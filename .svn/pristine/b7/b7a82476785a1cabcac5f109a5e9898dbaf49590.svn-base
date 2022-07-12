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
//            Cylinder x Cylinder with 4D two body generic interaction        //
//                                                                            //
//                        last modification : 08/12/2016                      //
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


#ifndef PARTICLEONCYLINDERXCYLINDERTWOBODYGENERICHAMILTONIAN_H
#define PARTICLEONCYLINDERXCYLINDERTWOBODYGENERICHAMILTONIAN_H

#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnS2xS2GenericTwoBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian : public ParticleOnS2xS2GenericTwoBodyHamiltonian
{

 protected:
  
  // ratio between the length in the x direction and the length in the y direction of the first cylinder
  double Ratio1; 
  // ratio between the length in the x direction and the length in the y direction of the second cylinder
  double Ratio2;

 public:

  // default constructor
  //
  ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuanta1 = number of flux quanta for the first cylinder
  // nbrFluxQuanta2 = number of flux quanta for the second cylinder
  // ratio1 = ratio between the length in the x direction and the length in the y direction of the first cylinder
  // ratio2 = ratio between the length in the x direction and the length in the y direction of the second cylinder
  // nbrPseudoPotentials = number of pseudo-potentials
  // pseudoPotentialAngularMomentum1= pseudo-potential first sphere relative angular momenta
  // pseudoPotentialAngularMomentum2 = pseudo-potential second sphere relative angular momenta
  // pseudoPotentials = pseudo-potential coefficients
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta1, int nbrFluxQuanta2, double ratio1, double ratio2, 
						       int nbrPseudoPotentials, int* pseudoPotentialAngularMomentum1, int* pseudoPotentialAngularMomentum2, 
						       double* pseudoPotentials, AbstractArchitecture* architecture, long memory = -1);
  
  // destructor
  //
  ~ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // evaluate the numerical coefficient  in front of the a+_(lz1,kz1) a+_(lz2,kz2) a_(lz3,kz3) a_(lz4,kz4) coupling term
  //
  // lz1 = first lz index
  // lz2 = second lz index
  // lz3 = third lz index
  // lz4 = fourth lz index
  // kz1 = first kz index
  // kz2 = second kz index
  // kz3 = third kz index
  // kz4 = fourth kz index
  // pseudopotentialIndex1 = pseudopotential index for the interaction on the first cylinder 
  // pseudopotentialIndex2 = pseudopotential index for the interaction on the second cylinder 
  // exponentialCoefficient1 = array that contains the precomputed exponential factors for the first cylinder  
  // exponentialCoefficient2 = array that contains the precomputed exponential factors for the second cylinder  
  // return value = numerical coefficient
  virtual double EvaluateInteractionCoefficient(int lz1, int lz2, int lz3, int lz4, 
						int kz1, int kz2, int kz3, int kz4, 
						int pseudopotentialIndex1, int pseudopotentialIndex2,
						double* exponentialCoefficient1, double* exponentialCoefficient2);

};

#endif
