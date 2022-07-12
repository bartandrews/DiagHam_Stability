////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//      class for a basic Monte Carlo algorith for particles on a sphere      //
//                                                                            //
//                        last modification : 23/01/2008                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef BLOCKMONTECARLOONSPHEREALGORITHM_H
#define BLOCKMONTECARLOONSPHEREALGORITHM_H


#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "AbstractObservable.h"
#include "AbstractBlockObservable.h"
#include "AbstractMCBlockSamplingFunctionOnSphere.h"
#include "AbstractParticleCollectionOnSphere.h"

#include <iostream>
using std::ostream;

class OptionManager;

class BlockMonteCarloOnSphereAlgorithm
{
 protected:

  // number of particles
  int NbrParticles;

  // number of Blocks
  int NbrBlocks;

  // sampling function
  AbstractMCBlockSamplingFunctionOnSphere *TargetFunction;  
  
  // class holding the particle coordinates
  AbstractParticleCollectionOnSphere *System;

  // pointer to the option manager
  OptionManager* Options;

  // array with all observables in the system
  AbstractObservable **Observables;

  // array with frequencies associated to observables
  int *Frequencies;

  // size of array Observables
  int MaxNbrObservables;

  // number of observable
  int NbrObservables;

    // array with all observables in the system
  AbstractBlockObservable **BlockObservables;

  // array with frequencies associated to observables
  int *BlockFrequencies;

  // size of array Observables
  int MaxNbrBlockObservables;

  // number of observable
  int NbrBlockObservables;

  // number of attempted microsteps
  long NbrAttemptedMoves;

  // number of accepted microsteps
  long NbrAcceptedMoves;

  
  
 public:

  // default constructor
  BlockMonteCarloOnSphereAlgorithm();

  // set up for basic monte-carlo scheme
  // nbrParticles = number of particles in system
  // waveFunction = wavefunction to be simulated
  // samplingFunction = function to be used to generate samples
  // manager = pointer to option manager
  // maxNbrObservables = maximum number of observables to be assigned
  BlockMonteCarloOnSphereAlgorithm(int nbrParticles, AbstractMCBlockSamplingFunctionOnSphere *targetFunction,
				    OptionManager *manager, int maxNbrObservables = 10);

  
  // destructor
  ~BlockMonteCarloOnSphereAlgorithm();

  // add an observable
  // O = pointer to the observable to be added
  // frequency = integer value indicating on which microsteps observations are being made
  void AddObservable(AbstractObservable *O, int frequency = 1);

  // add a block observable for detailed views of the contribution of individual terms in the symmetric expansion
  // O = pointer to the observable to be added
  // frequency = integer value indicating on which microsteps observations are being made
  void AddBlockObservable(AbstractBlockObservable *O, int frequency = 1);

  // thermalize system with a number of microsteps
  // time = number of microsteps
  // startFromRandom = flag indicating if we want to restart from a random configuration
  void Thermalize(int time, bool startFromRandom = false);

  // renormalize wavefunction
  // time = number of points to average
  void NormalizePsi(int time=500);

  // run simulation
  // all options are included via the AddOptionGroup method
  // Output = stream where to direct output
  void Simulate(std::ostream &Output = std::cout);

  //
  // perform a number of Monte-Carlo microsteps
  // nbrSteps = number of steps
  //
  void PerformMicroSteps(int nbrSteps);
  
  // add an option group containing all options related to 
  // this Monte-Carlo module
  //
  // manager = pointer to the option manager
  static void AddOptionGroup(OptionManager* manager);

  // write measurements of all observables to a given stream
  // str = output stream
  //
  void WriteObservations(ostream &str=std::cout);

  
  
};

#endif // BLOCKMONTECARLOONSPHEREALGORITHM_H
