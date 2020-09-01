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
    std::cout<<"Generating Error Kernel using cumulant method \n";
    cumulant = params.get<bool>("cumulant");
    if(!params.keyExists<int>("order")){
      std::cout<<"order not specified, setting order = 2 \n";
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
  }
  if (params.keyExists<std::vector<int>>("layout")){
    auto tmp = params.get<std::vector<int>>("layout");
    for (auto& a : tmp){
      layout.push_back(a);
    }
  }
} // initialize

void AssignmentErrorKernelDecorator::execute(
    std::shared_ptr<AcceleratorBuffer> buffer,
    const std::shared_ptr<CompositeInstruction> function){
  int num_bits = buffer->size();
  int pow_size = std::pow(2, num_bits);
  if (!layout.empty()) {
     function->mapBits(layout);
  }
  else{
    for(int i = 0; i < buffer->size(); i++){
      layout.push_back(i);
    }
  }

  std::vector<std::shared_ptr<CompositeInstruction>> circuits;
  if(cumulant){
    circuits = cumulantCircuits();
    circuits.push_back(function);
  }
  else{
    circuits = kernelCircuits(buffer);
    circuits.push_back(function);
  }

  decoratedAccelerator->execute(buffer, circuits);
  auto buffers = genKernel(buffer);
  //buffers are just the non-kernel buffers
  int shots = 0;
  for(auto &x: buffers[0]->getMeasurementCounts()){
    shots += x.second;
  }
  Eigen::VectorXd init_dist(pow_size);
  double sum = 0.0;
  std::cout<<"init dist:\n";
  for(auto &x:genPermutations(num_bits)){
    int idx = std::stoi(x, 0, 2);
    std::cout<<"bitstring: "<<x<<" ";
    init_dist(idx) = buffers.back()->computeMeasurementProbability(x);
    std::cout<<init_dist(idx)<<"\n";
    sum += init_dist(idx);
  }
  std::cout<<"\n";

  std::cout<<"sum of init_dist: \n";
  std::cout<<sum<<"\n";

  Eigen::VectorXd mit_dist;
  mit_dist = errorKernel*init_dist;
  std::cout<<"sum of mit_dist: \n";
  sum = 0.0;
  for(int i = 0; i < mit_dist.size(); i++){
    sum += mit_dist(i);
  }
  std::cout<<sum<<"\n";
  //finding shots that map below zero
  double new_shots = (double) shots;
  for(int i =0; i < mit_dist.size(); i++){
    if(mit_dist(i) < 0.0){
      new_shots += shots*floor(mit_dist(i) + 0.5);
      mit_dist(i) = 0.0;
    }
  }

  //updating distribution with new shot size
  mit_dist.normalize();

  std::map<std::string, double> orig_counts;
  this->permutations = genPermutations(num_bits);
  for(auto &x:permutations){
    orig_counts[x] = (double)buffers.back()->getMeasurementCounts()[x];
    int idx = std::stoi(x, 0, 2);
    int count = floor(new_shots * mit_dist(idx) + 0.5);
    std::cout<<"appending: "<<x<<" with: "<<count<<"\n";
    buffers.back()->appendMeasurement(x, count);
  }
  buffer->addExtraInfo("unmitigated-counts", orig_counts);
  std::cout<<"shot dropout is: "<<shots - new_shots<<"\n";
  gen_kernel = false;
  std::vector<double> vec(errorKernel.data(), errorKernel.data() + errorKernel.rows()*errorKernel.cols());
  buffer->addExtraInfo("error-kernel", vec);
  }//execute

void AssignmentErrorKernelDecorator::execute(
    const std::shared_ptr<AcceleratorBuffer> buffer,
    const std::vector<std::shared_ptr<CompositeInstruction>> functions) {
  int num_bits = buffer->size();
  int pow_size = std::pow(2, num_bits);
  int num_circs = functions.size();
  if(!layout.empty()){
    for(auto &program:functions){
      program->mapBits(layout);
    }
  }
  else{
    for(int i = 0; i < buffer->size(); i++){
      layout.push_back(i);
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
  //buffers of programs passed to the decorator (no kernel circuits)
  int shots = 0;
  for(auto &x: buffers[0]->getMeasurementCounts()){
    shots += x.second;
  }
  std::vector<Eigen::VectorXd> init_states;
  this->permutations = genPermutations(num_bits);
  for (auto &b : buffers) {
    Eigen::VectorXd temp(pow_size);
    for (auto &x : permutations) {
      int idx = std::stoi(x, 0, 2);
      temp(idx) =(double)b->computeMeasurementProbability(x);
    }
    init_states.push_back(temp);
  }
  std::vector<Eigen::VectorXd> mit_dist;
  for (auto &x : init_states) {
    mit_dist.push_back(errorKernel * x);
  }
  std::cout<<"distribution pre clipping:\n";
  for(int i = 0; i < mit_dist.size(); i++){
    std::cout<<"distribution: "<<i<<": ";
    for(int j = 0; j < mit_dist[i].size(); j++){
      std::cout<<mit_dist[i](j)<<" ";
    }
    std::cout<<"\n";
  }
  std::cout<<"\n";

  //finding shots that map below zero
  std::vector<std::map<std::string, double>> origCounts(buffer->nChildren());
  int b_idx = 0;
  for(auto &b: buffers){
    int new_shots = shots;
    for(auto &x: permutations){
      origCounts[b_idx][x] = (double)b->getMeasurementCounts()[x];
      int idx = std::stoi(x, 0, 2);
      if(mit_dist[b_idx](idx) < 0.0){
        new_shots += floor(shots*mit_dist[b_idx](idx) + 0.5);
        mit_dist[b_idx](idx) = 0;
      }
      int count = floor(new_shots*mit_dist[b_idx](idx) + 0.5);
      std::cout<<"appending: "<<x<<" with: "<< count<<"\n";
      b->appendMeasurement(x, mit_dist[b_idx](idx));
    }
    b->addExtraInfo("unmitigated-counts", origCounts[b_idx]);
    b_idx++;
  }

  std::vector<double> vec(errorKernel.data(), errorKernel.data() + errorKernel.rows()*errorKernel.cols());
  buffer->addExtraInfo("error-kernel", vec);
  gen_kernel = false;
  return;
} // execute (vectorized)

} // namespace quantum
} // namespace xacc
