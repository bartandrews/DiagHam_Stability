////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//     class of hamiltonian with particles on the 4D manifold S2 x S2         //
//                     with 4D generic two-body interaction                   //
//                                                                            //
//                        last modification : 07/12/2016                      //
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


#ifndef PARTICLEONS2XS2GENERICTWOBODYHAMILTONIAN_H
#define PARTICLEONS2XS2GENERICTWOBODYHAMILTONIAN_H

#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnS2xS2DeltaHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnS2xS2GenericTwoBodyHamiltonian : public ParticleOnS2xS2DeltaHamiltonian
{

 protected:
  
  // number of pseudo-potentials
  int NbrPseudoPotentials;
  // pseudo-potential first sphere relative angular momenta 
  int* PseudoPotentialAngularMomentum1;
  // pseudo-potential second sphere relative angular momenta
  int* PseudoPotentialAngularMomentum2;
  //pseudo-potential coefficients
  double* PseudoPotentials;

 public:

  // default constructor
  //
  ParticleOnS2xS2GenericTwoBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuanta1 = number of flux quanta for the first sphere
  // nbrFluxQuanta2 = number of flux quanta for the second sphere
  // nbrPseudoPotentials = number of pseudo-potentials
  // pseudoPotentialAngularMomentum1= pseudo-potential first sphere relative angular momenta
  // pseudoPotentialAngularMomentum2 = pseudo-potential second sphere relative angular momenta
  // pseudoPotentials = pseudo-potential coefficients
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnS2xS2GenericTwoBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta1, int nbrFluxQuanta2, 
					   int nbrPseudoPotentials, int* pseudoPotentialAngularMomentum1, int* pseudoPotentialAngularMomentum2, 
					   double* pseudoPotentials, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnS2xS2GenericTwoBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

};

#endif
