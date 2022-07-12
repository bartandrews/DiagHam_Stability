////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                         a generic n-body interaction                       //
//                                                                            //
//                        last modification : 23/06/2008                      //
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


#ifndef PARTICLEONSPHEREGENERICNBODYHAMILTONIAN_H
#define PARTICLEONSPHEREGENERICNBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnSphereNBodyHardCoreHamiltonian.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients;


class ParticleOnSphereGenericNBodyHamiltonian : public ParticleOnSphereNBodyHardCoreHamiltonian
{

 protected:

  // number of particle that interact simultaneously through the hard core interaction
  int NbrNbody;
  
  // weight of the different n-body interaction terms with respect to each other
  double* NBodyInteractionWeightFactors;

  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  double* PseudoPotential;

 public:

  // default constructor
  //
  ParticleOnSphereGenericNBodyHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // maxNbrBody = maximum number of particle that interact simultaneously through an interaction
  // maxRelativeMomenta = array that indicates the maximum relative angular momentum between n particles that may occur in the given n-body channel (negative if the channel is absent)
  // nBodyFactors = weight of the different n-body interaction terms with respect to each other (first entry being the n, second the relative angular momentum between n particles)
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereGenericNBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
					  int maxNbrBody, int* maxRelativeMomenta, double** nBodyFactors, double l2Factor = 0.0,
					  AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
					  char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereGenericNBodyHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();


 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();


};

#endif
