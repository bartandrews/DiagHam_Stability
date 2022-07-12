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


#ifndef ABSTRACTMCBLOCKSAMPLINGFUNCTIONONSPHERE_H
#define ABSTRACTMCBLOCKSAMPLINGFUNCTIONONSPHERE_H

#include "config.h"
#include "MathTools/Complex.h"
#include "AbstractMCBlockSamplingFunction.h"
#include "AbstractParticleCollectionOnSphere.h"


class AbstractMCBlockSamplingFunctionOnSphere : public AbstractMCBlockSamplingFunction
{
 protected:
  // collection of particles
  AbstractParticleCollectionOnSphere *System;
 public:
  // virtual destructor
  virtual ~AbstractMCBlockSamplingFunctionOnSphere();

  // method for ratio of probabilities with respect to the last configuration
  // allows for more rapid calculation due to cancellation of factors
  virtual double GetTransitionRatio()=0;

  // get the estimate of the full function value calculated over all Blocks for the given system of particles
  virtual Complex GetFunctionValue()=0;

  // get the estimate of the sampling function value calculated over the sampling block only, for the given system of particles
  virtual Complex GetSamplingBlockValue()=0;

  // call this method to scale the sampling function (needed to normalize the function)
  virtual void ScaleByFactor(double scale)=0;

  // accessor routine for NbrBlocks
  virtual int GetNbrBlocks()=0;

  // get number of flux quanta, or degree of underlying polynomial for simulation wavefunction
  virtual int GetNbrFluxQuanta() = 0;

  // get the Monte Carlo amplitude for the requested block with nbrPermute particles exchanged
  // nbrBlock = nbrPermute = number of particles to exchange between blocks
  // amplitude = reference of return argument
  virtual void GetBlockAmplitude(int nbrBlock, Complex &amplitude) = 0;

  // get the Monte Carlo amplitude for the requested block with nbrPermute particles exchanged
  // amplitudes = pointer to array of return arguments
  virtual void GetAllBlockAmplitudes(Complex *amplitudes) = 0;

  // query permutation of particles applied to given block
  // nbrBlock = block index
  // permutations = pointer to array to be filled with return argument
  virtual void GetBlockPermutations(int nbrBlock, int* permutations) = 0;

  // query permutation of particles applied to given block
  // nbrBlock = block index
  // permutations = pointer to array to be filled with return argument
  virtual void GetAllBlockPermutations(int** permutations) = 0;

  // query weights of individual blocks
  // weights = array reference to values of weighting factors
  virtual void GetBlockWeights(double *weights)=0;

    // register basic system of particles
  virtual void RegisterSystem(AbstractParticleCollection *system);

  // register basic system of particles
  virtual void RegisterSystem(AbstractParticleCollectionOnSphere *system) = 0;
  
 protected:
  // register basic system of particles
  virtual AbstractParticleCollection * GetSystem();  
  
};

#endif // ABSTRACTMCBLOCKSAMPLINGFUNCTION_H
