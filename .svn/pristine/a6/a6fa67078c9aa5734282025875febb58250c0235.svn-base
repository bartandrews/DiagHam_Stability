////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated to particles on a thin annulus      //
//           with generic interaction defined by its pseudopotential          //
//                                                                            //
//                        last modification : 07/11/2021                      //
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


#ifndef PARTICLEONTHINANNULUSGENERICHAMILTONIAN_H
#define PARTICLEONTHINANNULUSGENERICHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Hamiltonian/AbstractQHEOnSphereHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class ParticleOnThinAnnulusGenericHamiltonian : public AbstractQHEOnSphereHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:

  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  double* PseudoPotential;

  // array with the coefficient in front of each one body term (ordered such that the first element corresponds to the one of a+_-s a_-s)
  double* OneBodyPotentials;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzMax = maximum angular momentum that a single particle can reach
  // architecture = architecture to use for precalculation
  // pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction, V_m=\int d^2r r^2m V(r) e^(-r^2/8) )
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnThinAnnulusGenericHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzMax, double* pseudoPotential,
				   AbstractArchitecture* architecture, long memory = -1,
				   bool onDiskCacheFlag = false, char* precalculationFileName = 0);

  // constructor with one body terms
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzMax = maximum angular momentum that a single particle can reach
  // architecture = architecture to use for precalculation
  // pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction, V_m=\int d^2r r^2m V(r) e^(-r^2/8) )
  // oneBodyPotentials = array with the coefficient in front of each one body term (ordered such that the first element corresponds to the one of a+_-s a_-s)
  // mMax = maximum momentum for single orbitals
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnThinAnnulusGenericHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzMax, double* pseudoPotential, double* oneBodyPotentials, int mMax,
				   AbstractArchitecture* architecture, long memory = -1,
				   bool onDiskCacheFlag = false, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnThinAnnulusGenericHamiltonian();

 
 protected:

   // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();


};

#endif
