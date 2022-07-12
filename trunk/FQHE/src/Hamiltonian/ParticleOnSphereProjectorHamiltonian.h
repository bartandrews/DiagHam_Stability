////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of hamiltonian defined as a projector                 //
//                       for particles on a sphere defined                    //
//                                                                            //
//                        last modification : 14/08/2009                      //
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


#ifndef PARTICLEONSPHEREPROJECTORHAMILTONIAN_H
#define PARTICLEONSPHEREPROJECTORHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereNBodyInteractionHamiltonian.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients;


class ParticleOnSphereProjectorHamiltonian : public AbstractQHEOnSphereNBodyInteractionHamiltonian
{

 protected:


  // number of projectors
  int NbrProjectors;
  // vector states describing each projector
  RealVector* ProjectorStates;
  // number of particles for the projector state
  int ProjectorNbrParticles;
  // spaces associated to the projector states
  ParticleOnSphere** ProjectorSpaces;
  

 public:

  // default constructor
  //
  ParticleOnSphereProjectorHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // projectorState = state describing the projector
  // projectorSpace = space associated to the projector state 
  // projectorNbrParticles = number of particles for the projector state
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereProjectorHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
				       RealVector& projectorState, ParticleOnSphere* projectorSpace, int projectorNbrParticles, double l2Factor, 
				       AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
				       char* precalculationFileName = 0);

  // constructor from multiple projectors
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // projectorState = states describing each projector
  // projectorSpace = spaces associated to the projector states
  // projectorNbrParticles = number of projectors
  // projectorNbrParticles = number of particles for the projector state
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereProjectorHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
				       RealVector* projectorStates, ParticleOnSphere** projectorSpaces, int NbrProjectors, int projectorNbrParticles, double l2Factor, 
				       AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
				       char* precalculationFileName);

  // destructor
  //
  ~ParticleOnSphereProjectorHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();


 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

};

#endif
