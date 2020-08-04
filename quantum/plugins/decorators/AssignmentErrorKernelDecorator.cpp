/*******************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC.
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
#include "AssignmentErrorKernelDecorator.hpp"
#include "InstructionIterator.hpp"
#include "Utils.hpp"
#include "xacc.hpp"
#include <fstream>
#include <Eigen/Dense>


/*
const std::vector<std::vector<int>> mapData { { 1, 3}, {2, 4}, {5, 6} };
  xacc::HeterogeneousMap m{ {"myMap", mapData} };
  if (m.keyExists<std::vector<std::vector<int>>>("myMap")) {
    std::cout << "Has map\n";
    auto mapData = m.get<std::vector<std::vector<int>>>("myMap");
    std::cout << "Data length = " << mapData.size() << "\n";
  } 
*/

namespace xacc {
namespace quantum {
void AssignmentErrorKernelDecorator::initialize(
    const HeterogeneousMap &params) {

  if(params.keyExists<std::vector<std::vector<std::size_t>>>("cluster-map")){
    cluster_map = params.get<std::vector<std::vector<std::size_t>>>("cluster-map");
    }

  if(params.keyExists<bool>("spectators")){
    spectators = params.get<bool>("spectators");
  }
  if(params.keyExists<int>("order")){
    order = params.get<int>("order");
    std::cout<<"cumulant order specified: "<<order<<std::endl;
  }
  if(params.keyExists<bool>("cumulant")){
    std::cout<<"Generating Error kernel using cumulant method\n";
    cumulant = params.get<bool>("cumulant");
    if(order == 0){
      std::cout<<"order not specified, setting order = 2";
      order = 2;
    }
  }
  if(!cumulant && spectators){
    std::cout<<"must have cumulants == True in order to use spectators";
    spectators = false;
  }
  if (params.keyExists<bool>("gen-kernel")) {
    gen_kernel = params.get<bool>("gen-kernel");
  }
  if (params.keyExists<std::vector<std::size_t>>("layout")){
    layout = params.get<std::vector<std::size_t>>("layout");
    std::cout<<"layout recieved"<<std::endl;
  }
  if (params.keyExists<std::vector<int>>("layout")){
    auto tmp = params.get<std::vector<int>>("layout");
    std::cout<<"layout specified, num qubits: "<<tmp.size()<<std::endl;
    for (auto& a : tmp){
      layout.push_back(a);
    }
  }
} // initialize

void AssignmentErrorKernelDecorator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::shared_ptr<CompositeInstruction> function){
  int num_bits = buffer->size();
  std::cout<<"num qubits allocated: "<<buffer->size() <<std::endl;
  int pow_size = std::pow(2, num_bits);
  if (!layout.empty()) {
    std::cout<<"mapbits\n";
     function->mapBits(layout);
  }

  std::vector<std::shared_ptr<CompositeInstruction>> circuits;
  if(cumulant){
    std::cout<<"if cumulant\n";
    circuits = cumulantCircuits();
    circuits.push_back(function);
    std::cout<<"all circuits added to circuits vector\n";
  }
  else{
    circuits = kernelCircuits(buffer);
    std::cout<<"num kernel circuits: "<<circuits.size()<<"\n";
    circuits.push_back(function);
  }

  //vector of buffers the last one corresponding to the passed in user program
  decoratedAccelerator->execute(buffer, circuits);
  auto buffers = genKernel(buffer);
  int shots = 0;
  for(auto &x: buffers[0]->getMeasurementCounts()){
    shots += x.second;
  }
  Eigen::VectorXd init_dist(pow_size);
  int i = 0;
  for(auto &x:permutations){
    int idx = std::stoi(x, 0, 2);
    init_dist(idx) = buffers.back()->computeMeasurementProbability(x);
  }
  Eigen::VectorXd mit_dist;
  if(!cumulant){
    mit_dist = errorKernel.inverse()*init_dist;
  }
  else{
    mit_dist = errorKernel*init_d
  }
  

  std::map<std::string, double> orig_counts;
  for(auto &x:permutations){
    orig_counts[x] = (double)buffers.back()->getMeasurementCounts()[x];
    int count = floor(shots * mit_dist(i) + 0.5);
    buffers.back()->appendMeasurement(x, count);
    i++;
  }
  buffers.back()->addExtraInfo("unmitigated-counts", orig_counts);



  }

void AssignmentErrorKernelDecorator::execute(
    const std::shared_ptr<AcceleratorBuffer> buffer,
    const std::vector<std::shared_ptr<CompositeInstruction>> functions) {
  int num_bits = buffer->size();
  int pow_size = std::pow(2, num_bits);
  int num_circs = functions.size();
  std::cout<<"num circs passed: " <<num_circs<<"\n";
  if(!layout.empty()){
    for(auto &program:functions){
      program->mapBits(layout);
    }
  }
  std::vector<std::shared_ptr<CompositeInstruction>> circuits;
  if(cumulant){
    circuits = cumulantCircuits();
    for(auto & program:functions){
      circuits.push_back(program);
    }
  }
  else{
    circuits = kernelCircuits(buffer);
    for(auto & program:functions){
      circuits.push_back(program);
    }
    std::cout<<"all circuits size = "<<circuits.size()<<"\n";
  }


  // get the raw states
  decoratedAccelerator->execute(buffer, circuits);
  auto buffers = genKernel(buffer);
  int shots = 0;
  for(auto &x: buffers[0]->getMeasurementCounts()){
    shots += x.second;
  }
  for(int i = circuits.size()-functions.size(); i < buffers.size(); i++){
    if(pow_size < std::pow(2, buffers[i]->size())){
      std::cout<<"qubits allocated less than parent buffer, tracing error kernel over unused indices";
      auto bitmap = buffers[i]->getBitMap();
      //get values from bitmap
      //compare values from bitmap to values in layout
      //generate index list and pass that to trace, mitigate this circuit
      //with that specific kernel
    }
    Eigen::VectorXd(pow_size);
    int j = 0;
  }
  std::vector<Eigen::VectorXd> init_states;
  // compute number of shots;

  int i = 0;
  for (auto &b : buffers) {
    Eigen::VectorXd temp(pow_size);
    int j = 0;
    for (auto &x : permutations) {
      temp(j) =(double)b->computeMeasurementProbability(x);
      j++;
    }
    // std::cout << "seg fault right here " << std::endl;
    init_states.push_back(temp);
  }
//   std::cout << "unmitigated states: " << std::endl;
  std::vector<Eigen::VectorXd> EM_states;
  i = 0;
  for (auto &x : init_states) {
    std::cout <<"vector size: "<< x.size() << std::endl;
    std::cout<<"matrix size "<<errorKernel.cols()<<std::endl;
    EM_states.push_back(errorKernel * x);
    i++;
  }
//   std::cout << "mitigated states: " << std::endl;
//   for (auto &x : EM_states) {
//     std::cout << x << std::endl << std::endl;
//   }
  for (int i = 0; i < EM_states.size(); i++) {
    for (int j = 0; j < EM_states[i].size(); j++) {
      if (EM_states[i](j) < 0.0) {
        int count = floor(shots * EM_states[i](j) + 0.5);
        // std::cout << "found negative value, clipping and renorming "
        //           << std::endl;
        // std::cout << "removed " << abs(count) << " shots from total shots"
        //           << std::endl;
        shots += count;
        EM_states[i](j) = 0.0;
      }
    }
  }


  std::vector<std::map<std::string, double>> origCounts(buffer->nChildren());
  int total = 0;
  i = 0;
  for (auto &b : buffers) {
    int j = 0;
    for (auto &x : permutations) {
      origCounts[i][x] = (double)b->getMeasurementCounts()[x];
      int count = floor(shots * EM_states[i](j) + 0.5);
      j++;
      total += count;
      b->appendMeasurement(x, count);
      b->addExtraInfo("unmitigated-counts", origCounts[i]);
    }
    i++;
  }
  return;
} // execute (vectorized)

} // namespace quantum
} // namespace xacc
