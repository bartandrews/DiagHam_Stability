////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//             class of hamiltonian with particles on the 4D manifold         //
//      Cylinder x Cylinder with 4D delta interaction and truncated range     //
//                                                                            //
//                        last modification : 17/10/2016                      //
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


#ifndef PARTICLEONCYLINDERXCYLINDERDELTATRUNCATEDHAMILTONIAN_H
#define PARTICLEONCYLINDERXCYLINDERDELTATRUNCATEDHAMILTONIAN_H

#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnCylinderxCylinderDeltaHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnCylinderxCylinderDeltaTruncatedHamiltonian : public ParticleOnCylinderxCylinderDeltaHamiltonian
{

 protected:
  
 public:

  // default constructor
  //
  ParticleOnCylinderxCylinderDeltaTruncatedHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuanta1 = number of flux quanta for the first cylinder
  // nbrFluxQuanta2 = number of flux quanta for the second cylinder
  // ratio1 = ratio between the length in the x direction and the length in the y direction of the first cylinder
  // ratio2 = ratio between the length in the x direction and the length in the y direction of the second cylinder
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCylinderxCylinderDeltaTruncatedHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta1, int nbrFluxQuanta2, double ratio1, double ratio2,
						       AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnCylinderxCylinderDeltaTruncatedHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

};

#endif
