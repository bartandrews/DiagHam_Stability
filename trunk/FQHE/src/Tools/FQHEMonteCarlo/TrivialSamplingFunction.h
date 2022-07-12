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


#ifndef TRIVIALSAMPLINGFUNCTION_H
#define TRIVIALSAMPLINGFUNCTION_H

#include "AbstractMCSamplingFunction.h"
#include "Vector/RealVector.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"

class TrivialSamplingFunction : public AbstractMCSamplingFunction
{
 protected:
    // pointer to the ensemble of particles that shall be examined in MonteCarlo
  AbstractParticleCollection *System;
  
  // wavefunction that is being simulated
  Abstract1DComplexFunction *WaveFunction;

  // Memory of the last value of the wavefunction squared
  double LastAmplitude;

  // Memory of the last value for the wavefunction at an attempted move
  double TentativeNewValue;

  // private copy of vector describing particle positions
  RealVector Positions;
  
 public:
  // constructor
  // waveFunction = wavefunction to be sampled from
  TrivialSamplingFunction(Abstract1DComplexFunction *waveFunction);
  
  // virtual destructor
  virtual ~TrivialSamplingFunction();

  // register basic system of particles
  virtual void RegisterSystem(AbstractParticleCollection *system);

  // method for ratio of probabilities with respect to the last configuration
  // allows for more rapid calculation due to cancellation of factors
  virtual double GetTransitionRatio();

  // get the full function value for a system of particles
  virtual Complex GetFunctionValue();

  // signal that the last move was accepted
  virtual void AcceptedMove();

  // call this method to scale the sampling function (needed to normalize the function)
  // scale = total scaling factor
  virtual void ScaleByFactor(double scale);

 protected:
  // register basic system of particles
  virtual AbstractParticleCollection * GetSystem();  
  
};

#endif 
