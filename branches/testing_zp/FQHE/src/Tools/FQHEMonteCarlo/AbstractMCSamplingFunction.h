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


#ifndef ABSTRACTMCSAMPLINGFUNCTION_H
#define ABSTRACTMCSAMPLINGFUNCTION_H

#include "MathTools/Complex.h"
#include "AbstractParticleCollection.h"


class AbstractMCSamplingFunction
{
 protected:
  // pointer to the ensemble of particles that shall be examined in MonteCarlo
  AbstractParticleCollection *System;
  
 public:
  // virtual destructor
  virtual ~AbstractMCSamplingFunction();

  // register basic system of particles
  virtual void RegisterSystem(AbstractParticleCollection *system);

  // method for ratio of probabilities with respect to the last configuration
  // allows for more rapid calculation due to cancellation of factors
  virtual double GetTransitionRatio()=0;

  // get the full function value for a system of particles
  virtual Complex GetFunctionValue()=0;

  // set function value to one for present particle positions
  virtual void AdaptNorm();

  // signal that the last move was accepted
  virtual void AcceptedMove();
  
  // set function value for a typical average of particle positions in MC sampling
  //
  virtual void AdaptAverageMCNorm(int thermalize=500, int average=500);

  // call this method to scale the sampling function (needed to normalize the function)
  virtual void ScaleByFactor(double scale)=0;
  
};

#endif // ABSTRACTMCSAMPLINGFUNCTION_H
