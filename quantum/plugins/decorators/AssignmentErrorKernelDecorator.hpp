/*******************************************************************************
 * Copyright (c) 2018 UT-Battelle, LLC.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v1.0 which accompanies this
 * distribution. The Eclipse Public License is available at
 * http://www.eclipse.org/legal/epl-v10.html and the Eclipse Distribution
 *License is available at https://eclipse.org/org/documents/edl-v10.php
 *
 * Contributors:
 *   Tyler Kharazi - initial implementation
 *******************************************************************************/
#ifndef XACC_ASSIGNMENTERRORKERNELDECORATOR_HPP_
#define XACC_ASSIGNMENTERRORKERNELDECORATOR_HPP_

#include "AcceleratorDecorator.hpp"
#include <Eigen/Dense>
#include "xacc.hpp"
#include <unsupported/Eigen/KroneckerProduct>
#include <map>
#include <bitset>
#include <algorithm>

namespace xacc {
namespace quantum {
class AssignmentErrorKernelDecorator : public AcceleratorDecorator {
protected:
  //boolean designed to forgo excessive computation of kernels
  bool gen_kernel;
  Eigen::MatrixXd errorKernel;
  std::vector<std::size_t> layout;
  bool cumulant = false;
  bool multiplex;
  std::vector<std::vector<std::size_t>> cluster_map;
  bool spectators;
  int order = 0;
  std::vector<std::string> permutations;

  //Helper Functions
  //genPermutations ->  the set of all permutations of n_bits
  //genPartitions ->  the set of all partitions of layout set
  //genCombinations ->  the set of all combinations of layout with k <= order
  //sparsity ->  the number of elements with a 0 entry divided by total number of elements
  //trace -> distribution marginalized over the redundant qubits used by spectating


  std::vector<std::string> genPermutations(int num_bits) {
    int pow_bits = std::pow(2, num_bits);
    std::vector<std::string> bitstring(pow_bits);
    std::string str = "";
    std::string curr;
    int counter = 0;
    int j = num_bits;
    while (j--) {
      str.push_back('0');
    }
    for (int k = 0; k <= num_bits; k++) {
      str[num_bits - k] = '1';
      curr = str;
      do {
        bitstring[counter] = curr;
        counter++;
      } while (next_permutation(curr.begin(), curr.end()));
    }
    return bitstring;
  }

  std::vector<std::vector<std::vector<std::size_t>>> 
  genPartitions(const std::vector<std::size_t>& elements){
    std::vector<std::vector<std::vector<std::size_t>>> fList;
    std::vector<std::vector<std::size_t>> lists;
    std::vector<std::size_t> indexes(elements.size(), 0); // Allocate?
    lists.emplace_back(std::vector<std::size_t>());
    lists[0].insert(lists[0].end(), elements.begin(), elements.end());
    int counter = -1;
    for(;;){
        counter += 1;
        fList.emplace_back(lists);
        int i,index;
        bool obreak = false;
        for (i=indexes.size()-1;; --i) {
            if (i<=0){
                obreak = true;
                break;
              }
            index = indexes[i];
            lists[index].erase(lists[index].begin() + lists[index].size()-1);
            if (lists[index].size()>0)
                break;
            lists.erase(lists.begin() + index);
          }
        if(obreak) break;
        ++index;
        if (index >= lists.size())
            lists.emplace_back(std::vector<std::size_t>());
        for (;i<indexes.size();++i) {
            indexes[i]=index;
            lists[index].emplace_back(elements[i]);
            index=0;
          }
      }
    return fList;
  }

  std::vector<std::vector<std::size_t>> genCombinations(std::vector<std::size_t> list, int order){
    std::vector<std::vector<std::size_t>> combinations;
    int N = list.size();
    for(int i = 1; i <= order; i++){
      std::string bitmask(i, 1); // K leading 1's
      bitmask.resize(N, 0); // N-K trailing 0's
      // print integers and permute bitmask
      do {
          std::vector<std::size_t> comb_els;
          for (int i = 0; i < N; ++i) // [0..N-1] integers
          {
              if (bitmask[i]){
                comb_els.push_back(i);
              } 
          };
          std::vector<std::size_t> list_map;
          for(auto &x: comb_els){
            list_map.push_back(list[x]);
          }
          combinations.push_back(list_map);
      } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    }
    return combinations;
  }


  std::vector<std::shared_ptr<CompositeInstruction>> 
  genCircuits(std::vector<std::size_t> qubits){
    //generate all circuits needed for kernel on qubits
    if(spectators){
      std::cout<<"running with spectator qubits\n";
      return genCircuitsSpectators(qubits);
    }
    else{
      int num_bits = qubits.size();
      this->permutations = genPermutations(num_bits);
      std::vector<std::shared_ptr<CompositeInstruction>> circuits;
      auto provider = xacc::getIRProvider("quantum");
      int pow_bits = std::pow(2, qubits.size());
      for (int i = 0; i < pow_bits; i++) {
        auto circuit = provider->createComposite(permutations[i]);
        int j = num_bits-1;
        for (char c : permutations[i]) {
          if (c == '1') {
            auto x = provider->createInstruction("X", j);
            circuit->addInstruction(x);
          }
          j = (j-1)%num_bits;
        }
        for(int i = 0; i < num_bits; i++){
          circuit->addInstruction(provider->createInstruction("Measure", i));
        }
        if (!layout.empty()){
          circuit->mapBits(qubits);
        }
        circuits.push_back(circuit);
      }
      return circuits;
    }
  }

  std::vector<std::shared_ptr<CompositeInstruction>> 
  genCircuitsSpectators(std::vector<std::size_t> qubits){
    //generate circuits with spectator hadamards on qubits not in qubit, but in layout
    int num_bits = qubits.size();
    auto permutations = genPermutations(qubits.size());
    std::vector<int> qubit_idxs_in_layout;
    for(auto &x:qubits){
      auto it = std::find(layout.begin(), layout.end(), x);
      int idx = std::distance(layout.begin(), it);
      qubit_idxs_in_layout.push_back(idx);
    }
    std::vector<std::shared_ptr<CompositeInstruction>> circuits; 
    auto provider = xacc::getIRProvider("quantum");
    int pow_bits = std::pow(2, qubits.size());
    for (int i = 0; i < pow_bits; i++) {
        auto circuit = provider->createComposite(permutations[i]);
        int j = num_bits-1;
        for (char c : permutations[i]) {
          if (c == '1') {
            auto x = provider->createInstruction("X", qubit_idxs_in_layout[j]);
            circuit->addInstruction(x);
          }
          j = (j-1)%num_bits;
        }
        for(int i = 0; i < layout.size(); i++){
          if(std::find(qubits.begin(), qubits.end(), layout[i]) == qubits.end()){
            auto H = provider->createInstruction("H", i);
            circuit->addInstruction(H);
          } 
        }
        for(int i = 0; i < layout.size(); i++){
          circuit->addInstruction(provider->createInstruction("Measure", i));
        }
        if (!layout.empty()){
          circuit->mapBits(layout);
        }
        circuits.push_back(circuit);
      }
      return circuits;
    }

  std::string vecToString(std::vector<std::size_t> layout){
    std::string str = "";
    for(auto &x: layout){
      if(x == layout.back()){
        str += std::to_string(x);
      }
      else{
        str += std::to_string(x) + " ";
      }
    }
    return str;
  }

  template <typename T>
  void printVec(std::vector<T> vec){
    for(auto &x:vec){
      std::cout<<x<<" ";
    }
    std::cout<<"\n";
  }


  Eigen::MatrixXd trace(std::vector<std::vector<double>>& dists, std::vector<std::size_t>& idxs){
    //Return subspace of only qubits on idxs
    int new_size = std::pow(2, idxs.size());
    Eigen::MatrixXd new_mat(new_size, new_size);
    new_mat.setZero(new_size,new_size);
    int x = 0;
    std::reverse(idxs.begin(), idxs.end());
    for(auto &dist:dists){
      int y = 0;
      for(double &val:dist){
        std::string bits_y = std::bitset< 64 >( y ).to_string();
        std::string new_y = "";
        for(std::size_t &idx: idxs){
          new_y += bits_y[63-idx];
        }

        int y_idx = std::stoi(new_y, 0, 2);
        new_mat(x, y_idx) += val;
        y++;
      }
      x++;
    }
    return new_mat;
  }

  std::vector<std::shared_ptr<CompositeInstruction>>
  cumulantCircuits(){
    //each vector of compositeinstruction corresponds to the circuits needed to execute for
    //the kernel. The set of all these instructions represent the set of all circuits
    //that are needed to generate the cumulant
    std::map<std::string,std::vector<std::shared_ptr<CompositeInstruction>>> cum_circuits;
    std::vector<std::shared_ptr<CompositeInstruction>> all_circuits;
    if(cluster_map.size() > 0){
      for(auto &x: this -> cluster_map){
        printVec(x);
        for(auto &circ:genCircuits(x)){
          all_circuits.push_back(circ);
        }
      }
    }
    else{
      std::cout<<"generating cumulant circuits\n";
      auto combinations = genCombinations(layout, this->order);
      for(auto &x:combinations){
        for(auto &circ:genCircuits(x)){
          all_circuits.push_back(circ);
        }
      }
    }
    return all_circuits;
  }

  std::vector<std::shared_ptr<CompositeInstruction>> 
  kernelCircuits(std::shared_ptr<AcceleratorBuffer> buffer){
    int num_bits = buffer->size();
    int pow_bits = std::pow(2, num_bits);
    this->permutations = genPermutations(num_bits);
    std::vector<std::shared_ptr<CompositeInstruction>> circuits;
    auto provider = xacc::getIRProvider("quantum");
    for(int i = 0; i < pow_bits; i++){
      auto circuit = provider->createComposite(permutations[i]);
      int j = num_bits - 1;
      for(char c : permutations[i]){
        if(c == '1'){
          auto x = provider->createInstruction("X", j);
          circuit->addInstruction(x);
        }
        j = (j-1)%num_bits;
      }
      for(int i = 0; i < num_bits; i++){
        circuit->addInstruction(provider->createInstruction("Measure", i));
      }
      if (!layout.empty()){
        circuit->mapBits(layout);
      }
      circuits.push_back(circuit);
    }
    return circuits;
  }



  std::vector<std::shared_ptr<AcceleratorBuffer>> 
  genKernel(std::shared_ptr<AcceleratorBuffer> buffer){
    if(cumulant){
      std::cout<<"generating cumulant kernel from buffer result\n";
      return genCumulantKernel(buffer);
    }
    int num_bits = buffer->size();
    int pow_bits = std::pow(2,num_bits);
    this->permutations = genPermutations(num_bits);
    Eigen::MatrixXd kernel(pow_bits, pow_bits);
    kernel.setZero();
    auto buffers = buffer->getChildren();
    for(int i = 0; i < pow_bits; i++){
      int col = 0;
      for(auto &x: permutations){
        auto temp = buffers[i]->computeMeasurementProbability(x);
        kernel(i, col) = temp;
        col ++;
      }
    }
    this->errorKernel = kernel.inverse();
    std::cout<<"NAIVE KERNEL:\n";
    print_mat(errorKernel);

    buffers.erase(buffers.begin(), buffers.begin()+pow_bits);
    gen_kernel = false;
    return buffers;
  }

  Eigen::VectorXd normalize_vec(Eigen::VectorXd& vec){

  }

  Eigen::MatrixXd normalize_mat(Eigen::MatrixXd& mat){
    int i = 0;
    double run_sum = 0;
    for(int j = 0; j < mat.cols(); j++){
      run_sum += mat(i,j);
    }
    return 1/run_sum * mat;
  }
  void print_mat(Eigen::MatrixXd& mat){
    for(int i = 0; i < mat.rows(); i++){
      for(int j = 0; j < mat.cols(); j++){
        std::cout<<mat(i,j)<<" ";
      }
      std::cout<<"\n";
    }
  }

  std::vector<std::size_t> strToVec(std::string vec_str){
    std::stringstream iss( vec_str );
    int number;
    std::vector<std::size_t> myNumbers;
    while ( iss >> number ){
      myNumbers.push_back( number );
    }
      return myNumbers;
  }

  std::vector<std::shared_ptr<AcceleratorBuffer>>
   genCumulantKernel(std::shared_ptr<AcceleratorBuffer> buffer){
    auto buffers = buffer->getChildren();
    //go through buffers to get kernel expermiments
    std::map<std::string, Eigen::MatrixXd> kernels;
    //global index counter for buffers
    int glob_idx = 0;
    //putting all kernels into map object for later retrieval
    std::vector<std::vector<std::size_t>> combinations;
    if(cluster_map.size()>0){
      combinations = cluster_map;
    }
    else{
      combinations = genCombinations(this->layout, this->order);
    }

    int shot_size = 0;
    for(auto &x:buffers[0]->getMeasurementCounts()){
      shot_size += x.second;
    }

    for(auto &comb : combinations){
      int pow_size = std::pow(2, comb.size());
      //need to trace out dependence on spectators for kernel
      if(spectators){
        std::vector<std::vector<double>> kernel(std::pow(2, comb.size()), std::vector<double>(std::pow(2, layout.size())));
        for(auto &perm:genPermutations(comb.size())){
          for(auto &perm_y:genPermutations(this->layout.size())){
            int x = std::stoi(perm, 0, 2);
            int y = std::stoi(perm_y, 0, 2);
            kernel[x][y] = buffers[glob_idx]->computeMeasurementProbability(perm_y);
          }
          glob_idx++;
        }
        //generate indices of combination to trace
        std::vector<std::size_t> comb_idxs;
        for(auto &x: comb){
          auto it = std::find(layout.begin(), layout.end(), x);
          int idx = std::distance(layout.begin(), it);
          comb_idxs.push_back(idx);
        }
        auto mat = trace(kernel, comb_idxs);
        kernels[vecToString(comb)] = mat;
      }//if spectators
      else{
        Eigen::MatrixXd kernel(pow_size,pow_size);
        for(auto &perm : genPermutations(comb.size()) ){
          int x = stoi(perm, 0, 2);
          for(auto &perm_y : genPermutations(comb.size()) ){
            int y = stoi(perm_y, 0, 2);
            kernel(x,y) =  buffers[glob_idx]->computeMeasurementProbability(perm_y);
          }
          glob_idx++;
        }
        //kernels are being generated properly then
        kernels[vecToString(comb)] = kernel;

      }//else
    }//end iteration over combinations 
    //kernels map object generated with key: combination str, val: kernel mat

    //kernels saved now need to generate cumulant map object:
    std::map<std::string, Eigen::MatrixXd> cumulants;
    if(cluster_map.size()==0){
      for(int ord = 1; ord <= this->order; ord++){
        for(const auto &comb : genCombinations(layout, ord)){
          if(comb.size() == 1){
            cumulants[vecToString(comb)] = kernels[vecToString(comb)];
          }
          else{
            auto partitions = genPartitions(comb);
            auto kernel = kernels[vecToString(comb)];
            for(auto &partition:partitions){
              Eigen::MatrixXd cumulant = Eigen::MatrixXd::Identity(1,1);
              for(auto &set:partition){
                if(set.size() == comb.size()){
                  continue;
                  std::cout<<"getting next partition element\n";
                }
                else{
                  cumulant = Eigen::KroneckerProduct(cumulant, cumulants[vecToString(set)]).eval();
                  
                }
              }
              if(cumulant.size() == kernel.size()){
                kernel -= cumulant;
              }
            }
              cumulants[vecToString(comb)] = kernel;
          }
        }
      }
    }

    //cumulant map made
    //don't need circuits involving kernel generation anymore
    buffers.erase(buffers.begin(), buffers.begin()+glob_idx);

    int pow_size = std::pow(2, this->layout.size()); 
    Eigen::MatrixXd cumulant_kernel(pow_size, pow_size);
    cumulant_kernel.setZero(pow_size,pow_size);
    //do matrix inversion with cumulant as (A+D)^-1 -> A^-1 + A^-1DA^-1
    if(cluster_map.size() == 0){
      auto partitions = genPartitions(this->layout);
      Eigen::MatrixXd single_kernel = Eigen::MatrixXd::Identity(1,1);
      for(auto &x:layout){
          auto kernel = kernels[std::to_string(x)];
          single_kernel = Eigen::KroneckerProduct(single_kernel, kernel).eval();
      }
      single_kernel = normalize_mat(single_kernel);
      //working properly
      //generating cumulants of higher order
      for(auto &partition:partitions){
        Eigen::MatrixXd sub_kernel = Eigen::MatrixXd::Identity(1,1);
        if(partition.size() == layout.size()){
          continue;
        }
        else{
          for(auto &set:partition){
            if(set.size() > order){
              //getting a new partition
              break;
            }
            else{
              auto cumulant = cumulants[vecToString(set)];
              sub_kernel = Eigen::KroneckerProduct(sub_kernel, cumulant).eval();
            }
          }
          if(sub_kernel.size() != cumulant_kernel.size()){
            continue;
          }
          else{
            cumulant_kernel += sub_kernel;
          }
        }
      }//iterate over partitions
      auto inv_ker = single_kernel.inverse();
      this->errorKernel = inv_ker + inv_ker*cumulant_kernel*inv_ker;
      std::cout<<"ERROR KERNEL CUMULANT: \n";
      print_mat(errorKernel);
    }//end cumulant kernel gen
    else{
      Eigen::MatrixXd sub_kernel = Eigen::MatrixXd::Identity(1,1);
      for(auto &x:cluster_map){
        auto kernel = kernels[vecToString(x)];
        sub_kernel = Eigen::KroneckerProduct(sub_kernel, kernel.inverse()).eval();
      }
      cumulant_kernel = normalize_mat(sub_kernel);
      std::cout<<"CLUSTER KERNEL:\n";
      print_mat(cumulant_kernel);
      this->errorKernel = cumulant_kernel;
    }//end clustering kernel gen
    return buffers;
  }//genCumulantKernel

  

public:
  AssignmentErrorKernelDecorator() = default;
  
  const std::vector<std::string> configurationKeys() override {
    return {"gen-kernel"};
  }

  void initialize(const HeterogeneousMap &params = {}) override;

  void execute(std::shared_ptr<AcceleratorBuffer> buffer,
               const std::shared_ptr<CompositeInstruction> function) override;
  void execute(std::shared_ptr<AcceleratorBuffer> buffer,
               const std::vector<std::shared_ptr<CompositeInstruction>>
                   functions) override;

  const std::string name() const override { return "assignment-error-kernel"; }
  const std::string description() const override { return ""; }

  ~AssignmentErrorKernelDecorator() override {}
};
} // namespace quantum
} // namespace xacc
#endif
