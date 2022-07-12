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


#ifndef HALPERINSAMPLINGFUNCTION_H
#define HALPERINSAMPLINGFUNCTION_H

#include "config.h"
#include "AbstractMCSamplingFunctionOnSphere.h"

class HalperinSamplingFunction : public AbstractMCSamplingFunctionOnSphere
{
 protected:
  // number of particles: total and with spin up
  int NbrParticles;  

  int NbrUp;
  // exponent of Jastrow Factors
  int Exponent_K;
  int Exponent_L;
  int Exponent_M;

  // pointers to spinor coordinates (external)
  Complex *SpinorUCoordinates;
  Complex *SpinorVCoordinates;

  // for access to ParticleCollection
  Complex LastU;
  Complex LastV;

  // norm for total function value
  double ElementNorm;
  
 public:
  // constructor
  // nbrParticles = number of particles in system
  // exponentK = power raised for up-spins
  // exponentL = power raised for down-spins
  // exponentM = power raised for inter-spin factors
  // nbrUp = number of particles with up spin (negative number -> choose half the particles)
  HalperinSamplingFunction(int nbrParticles, int exponentK, int exponentL, int exponentM, int nbrUp=-1);
    
  // virtual destructor
  virtual ~HalperinSamplingFunction();

  // register basic system of particles
  virtual void RegisterSystem(AbstractParticleCollectionOnSphere *system);

  // method for ratio of probabilities with respect to the last configuration
  // allows for more rapid calculation due to cancellation of factors
  virtual double GetTransitionRatio();

  // get the full function value for a system of particles
  virtual Complex GetFunctionValue();

  // call this method to scale the sampling function (needed to normalize the function)
  // scale = total scaling factor
  virtual void ScaleByFactor(double scale);
  
};

#endif 
