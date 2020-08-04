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

  std::vector<std::vector<std::vector<std::size_t>>> genPartitions(const std::vector<std::size_t>& elements){
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


  std::vector<std::shared_ptr<CompositeInstruction>> genCircuits(std::vector<std::size_t> qubits){
    //generate all circuits needed for kernel on qubits
    if(spectators){
      std::cout<<"running with spectator qubits\n";
      return genCircuitsSpectators(qubits);
    }
    else{
      std::cout<<"\n";
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

  std::vector<std::shared_ptr<CompositeInstruction>> genCircuitsSpectators(std::vector<std::size_t> qubits){
    //generate circuits with spectator hadamards on qubits not in qubits, but in layout
    int num_bits = layout.size();
    this->permutations = genPermutations(qubits.size());
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
          else{
            auto x = provider->createInstruction("H", j);
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

  std::string vecToString(std::vector<std::size_t> layout){
    std::string str = "";
    for(auto &x: layout){
      str += std::to_string(x);
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
    int num_qubits = (this->layout).size();
    int new_size = std::pow(2, idxs.size());
    Eigen::MatrixXd new_mat(new_size, new_size);
    int x = 0;
    for(auto &dist:dists){
      int y = 0;
      for(double &val:dist){
        y++;
        std::string bits_x = std::bitset< 64 >( x ).to_string();
        std::string bits_y = std::bitset< 64 >( y ).to_string();
        std::string new_x = "";
        std::string new_y = "";
        for(std::size_t &idx: idxs){
          new_x += bits_x[idx];
          new_y += bits_y[idx];
        }
        new_mat(std::stoi(new_x, 0, 2), std::stoi(new_y, 0, 2)) += val;
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
      std::cout<<"cluster map invoked\n";
      for(auto &x: this -> cluster_map){
        for(auto &circ:genCircuits(x)){
          all_circuits.push_back(circ);
        }
      }
      std::cout<<"cluster circuits:\n";
      for(auto &x:all_circuits){
        std::cout<<x->toString()<<std::endl;
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
    std::cout<<"generated " <<circuits.size()<<" circuits\n";
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

    this->errorKernel = kernel;

    buffers.erase(buffers.begin(), buffers.begin()+pow_bits);
    std::cout<<"buffers for kernel generation removed \n ";
    std::cout<<"new size: "<<buffers.size()<<"\n";

    std::vector<double> vec(errorKernel.data(), errorKernel.data() + errorKernel.rows()*errorKernel.cols());
    buffer->addExtraInfo("error-kernel", vec);
    gen_kernel = false;
    std::cout<<"sucesful run of genKernel\n";
    return buffers;
  }
  std::vector<std::shared_ptr<AcceleratorBuffer>>
   genCumulantKernel(std::shared_ptr<AcceleratorBuffer> buffer){
    auto buffers = buffer->getChildren();

    //go through buffers to get kernel expermiments
    std::map<std::string, Eigen::MatrixXd> kernels;
    //global index counter for buffers
    int i = 0;
    //putting all kernels into map object for later retrieval
    std::vector<std::vector<std::size_t>> combinations;
    if(cluster_map.size()>0){
      combinations = cluster_map;
    }
    else{
      combinations = genCombinations(this->layout, this->order);
    }
    for(auto &comb : combinations){
      int pow_size = std::pow(2, comb.size());
      if(spectators){
        std::vector<std::vector<double>> kernel;
        for(auto &perm:genPermutations(comb.size())){
          for(auto &perm_y:genPermutations(this->layout.size())){
            int x = std::stoi(perm, 0, 2);
            int y = std::stoi(perm_y, 0, 2);
            kernel[x][y] = buffers[x]->computeMeasurementProbability(perm_y);
          }
          i++;
        }
        kernels[vecToString(comb)] = trace(kernel, comb);
      }
      else{
        Eigen::MatrixXd kernel(pow_size,pow_size);
        for(auto &perm : genPermutations(comb.size()) ){
          for(auto &perm_y : genPermutations(comb.size()) ){
            int x = stoi(perm, 0, 2);
            int y = stoi(perm_y, 0, 2);
            kernel(x,y) = buffers[i]->computeMeasurementProbability(perm_y);
          }
          i++;
        }
      kernels[vecToString(comb)] = kernel;
      }
    }
    //don't need circuits involving kernel generation anymore
    buffers.erase(buffers.begin(), buffers.begin()+i);
    std::cout<<"circuits not used for kernel generation erased\n";
    std::cout<<"num circuits left: "<<buffers.size()<<"\n";

    auto partitions = genPartitions(this->layout);
    int pow_size = std::pow(2, this->layout.size());
    Eigen::MatrixXd cumulant_kernel(pow_size, pow_size);
    for(auto &partition:partitions){
      Eigen::MatrixXd sub_kernel(1, 1);
      for(auto &set:partition){
        if(set.size() > this->order){
          std::cout<<"leaving this partition, set "<<vecToString(set)<<" too big \n";
          break;
        }
        else{
          auto kernel = kernels[vecToString(set)].inverse();
          sub_kernel = Eigen::KroneckerProduct(sub_kernel, kernel).eval();
          std::cout<<"sub kernel cols: " <<sub_kernel.cols() <<" rows: " << sub_kernel.rows()<<"\n";
        }
      }
      if(sub_kernel.size() == cumulant_kernel.size()){
        cumulant_kernel += sub_kernel;
      }
      else{
        std::cout<<"matrix sizes don't match\n";
      }
      
    }

    this->errorKernel = cumulant_kernel;

    return buffers;
  }
  

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
