////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of eta pairing J^2 hamiltonian for interacting spinful particles  //
//                             for the Hubbard model                          //
//                                                                            //
//                        last modification : 09/02/2016                      //
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


#ifndef PARTICLEONLATTICEWITHSPINREALSPACEETAPAIRINGJ2HAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINREALSPACEETAPAIRINGJ2HAMILTONIAN_H


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


class ParticleOnLatticeWithSpinRealSpaceEtaPairingJ2Hamiltonian : public ParticleOnLatticeWithSpinRealSpaceHamiltonian
{

 protected:

  // numerical factor in front of J^2
  double J2Factor;

  // number of sites along the x direction
  int NbrSitesX;
  // number of sites along the y direction
  int NbrSitesY;

 public:

  // default constructor
  //
  ParticleOnLatticeWithSpinRealSpaceEtaPairingJ2Hamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSitesX = number of sites along the x direction
  // nbrSitesY = number of sites along the y direction
  // j2Factor = numerical factor in front of J^2
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeWithSpinRealSpaceEtaPairingJ2Hamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSitesX, int nbrSitesY, double j2Factor,
							    AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  virtual ~ParticleOnLatticeWithSpinRealSpaceEtaPairingJ2Hamiltonian();

  
 protected:
 
  // evaluate the two body interaction factors from a generic density-density interaction
  //
  virtual void EvaluateInteractionFactors();

};


#endif
