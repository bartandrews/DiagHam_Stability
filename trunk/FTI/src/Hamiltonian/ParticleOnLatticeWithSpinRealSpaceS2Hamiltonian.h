////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of S^2 hamiltonian for interacting spinful particles        //
//                       on lattice written in real space                     //
//                                                                            //
//                        last modification : 13/06/2015                      //
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


#ifndef PARTICLEONLATTICEWITHSPINREALSPACES2HAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINREALSPACES2HAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian : public ParticleOnLatticeWithSpinRealSpaceHamiltonian
{

 protected:

  // numerical factor in front of S^2
  double S2Factor;
  // true if the Hilbert space gas a fixed a total Sz value
  bool FixedSzFlag;
 
 public:

  // default constructor
  //
  ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSites = number of sites
  // s2Factor = numerical factor in front of S^2
  // fixedSzFlag = true if the Hilbert space gas a fixed a total Sz value
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSites, 
						  double s2Factor, bool fixedSzFlag,
						  AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  virtual ~ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian();

  
 protected:
 
  // evaluate the two body interaction factors from a generic density-density interaction
  //
  virtual void EvaluateInteractionFactors();

};


#endif
