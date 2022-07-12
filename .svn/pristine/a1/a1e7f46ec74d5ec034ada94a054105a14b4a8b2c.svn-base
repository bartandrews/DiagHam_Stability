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
//                           with 4D delta interaction                        //
//                                                                            //
//                        last modification : 27/09/2016                      //
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


#ifndef PARTICLEONS2XS2DELTAHAMILTONIAN_H
#define PARTICLEONS2XS2DELTAHAMILTONIAN_H

#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereFullHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnS2xS2DeltaHamiltonian : public AbstractQHEOnSphereFullHamiltonian
{

 protected:
  
  // number of flux quanta for the first sphere
  int NbrFluxQuanta1;
  // number of flux quanta for the second sphere
  int NbrFluxQuanta2;
  // number of orbitals for the first sphere
  int NbrLzValues1;
  // number of orbitals for the second sphere
  int NbrLzValues2;


 public:

  // default constructor
  //
  ParticleOnS2xS2DeltaHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrFluxQuanta1 = number of flux quanta for the first sphere
  // nbrFluxQuanta2 = number of flux quanta for the second sphere
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnS2xS2DeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta1, int nbrFluxQuanta2, 
				  AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnS2xS2DeltaHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // get a linearized index from the two angular momenta
  //
  // lz = z projection of the angular momentum for the first sphere
  // kz = z projection of the angular momentum for the second sphere
  // return value = linearized index 
  virtual int GetLinearizedIndex(int lz, int kz);

  // get the two angular momenta associated to a given linearized index
  //
  // index = linearized index 
  // lz = reference on the z projection of the angular momentum for the first sphere
  // kz = reference on the  z projection of the angular momentum for the second sphere
  virtual void GetLinearizedIndex(int index, int& lz, int& kz);

};

// get a linearized index from the two angular momenta
//
// lz = z projection of the angular momentum for the first sphere
// kz = z projection of the angular momentum for the second sphere
// return value = linearized index 

inline int ParticleOnS2xS2DeltaHamiltonian::GetLinearizedIndex(int lz, int kz)
{
  return ((lz * this->NbrLzValues2) + kz);
}

// get the two angular momenta associated to a given linearized index
//
// index = linearized index 
// lz = reference on the z projection of the angular momentum for the first sphere
// kz = reference on the  z projection of the angular momentum for the second sphere

inline void ParticleOnS2xS2DeltaHamiltonian::GetLinearizedIndex(int index, int& lz, int& kz)
{
  lz = index / this->NbrLzValues2;
  kz = index % this->NbrLzValues2;
}

#endif
