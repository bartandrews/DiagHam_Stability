////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                         n-body hard core interaction                       //
//                                                                            //
//                        last modification : 23/09/2004                      //
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


#ifndef PARTICLEONSPHERENBODYHARDCOREHAMILTONIAN_H
#define PARTICLEONSPHERENBODYHARDCOREHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereNBodyInteractionHamiltonian.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients;


class ParticleOnSphereNBodyHardCoreHamiltonian : public AbstractQHEOnSphereNBodyInteractionHamiltonian
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
  ParticleOnSphereNBodyHardCoreHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // nbrBody = number of particle that interact simultaneously through the hard core interaction
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int nbrBody, double l2Factor, 
					   AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
					   char* precalculationFileName = 0);

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // maxNbrBody = maximum number of particle that interact simultaneously through the hard core interaction
  // nBodyFactors = weight of the different n-body interaction terms with respect to each other
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
					   int maxNbrBody, double* nBodyFactors, double l2Factor,
					   AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
					   char* precalculationFileName = 0);

  // constructor from datas with a fully-defined two body interaction
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // maxNbrBody = maximum number of particle that interact simultaneously through the hard core interaction
  // nBodyFactors = weight of the different n-body interaction terms with respect to each other
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
					   int maxNbrBody, double* nBodyFactors, double l2Factor, double* pseudoPotential, 
					   AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
					   char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereNBodyHardCoreHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();


 protected:
 
  // compute all projector coefficient associated to a given relative angular momentum between k particles
  //
  // relativeMomentum = value of twice the relative angular momentum between the k particles
  // nbrIndices = number of indices per set
  // indices = array that contains all possible sets of indices (size of the array is nbrIndices * nbrIndexSets)
  // nbrIndexSets = number of sets
  double* ComputeProjectorCoefficients(int relativeMomentum, int nbrIndices, int* indices, int nbrIndexSets);

  // compute a given projector coefficient for the 4-body interaction 
  //
  // m1 = first index
  // m2 = second index
  // m3 = third inde
  // m4 = fourth index
  // jValue = total angular momentum
  // clebshArray = array that conatins Clebsch Gordan coefficients for each angular momentum value
  // return value = corresponding projector coefficient
  double ComputeProjectorCoefficients4Body(int m1, int m2, int m3, int m4, int jValue, 
					   ClebschGordanCoefficients* clebshArray);

  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();


};

#endif
