#ifndef __POMPP_Markov_Chain_HPP__
#define __POMPP_Markov_Chain_HPP__

// Abstract class meant to be used for MCMC calculations

#include <R.h>
//#include "MCMCprogress.hpp"

class MarkovChain {
  double logPosterior;
  unsigned int iteration;
//  MCMCprogress* progressBar;

  virtual double applyTransitionKernel() = 0; // returns log posterior
  void doNothing();
public:
//  MarkovChain() : iteration(0), progressBar(new MCMCprogress()) {}
  MarkovChain() : iteration(0) {}
  virtual ~MarkovChain() {}

  // Update the chain
  void update() {
    GetRNGstate();
    logPosterior = applyTransitionKernel();
    PutRNGstate();
//    progressBar->increment(1);
    iteration++;
  }
  void update(unsigned int times) {
    for (unsigned int i = 0; i < times; i++) update();
  }

  // Getters
  double getLogPosterior() {return logPosterior;}
  unsigned int getIteration() {return iteration;}
};

#endif
