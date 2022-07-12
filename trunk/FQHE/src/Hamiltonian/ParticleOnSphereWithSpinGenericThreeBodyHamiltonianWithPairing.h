////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                     spin and generic 3-body interaction                    //
//                                                                            //
//                        last modification : 26/08/2008                      //
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


#ifndef PARTICLEONSPHEREWITHSPINGENERICTHREEBODYHAMILTONIANWITHPAIRING_H
#define PARTICLEONSPHEREWITHSPINGENERICTHREEBODYHAMILTONIANWITHPAIRING_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients;


class ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing : public AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonianWithPairing
{

 protected:

  // array with the three-body pseudo-potentials in the S=3/2 sector
  double* ThreeBodyPseudoPotentials32;
  // number of elements in the ThreeBodyPseudoPotential array in the S=3/2 sector
  int NbrThreeBodyPseudoPotentials32;
  // maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential32
  int MaxRelativeAngularMomentum32;

  // array with the three-body pseudo-potentials in the S=1/2 sector
  double* ThreeBodyPseudoPotentials12;
  // number of elements in the ThreeBodyPseudoPotential array in the S=3/2 sector
  int NbrThreeBodyPseudoPotentials12;
  // maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential32
  int MaxRelativeAngularMomentum12;

  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  // first index refered to the spin sector (sorted as up-up, down-down, up-down)
  double** PseudoPotentials;

  // magnitude of contact pairing terms
  double Alpha;

  
 public:

  // default constructor
  //
  ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // threeBodyPseudoPotential32 = array with the three-body pseudo-potentials in the S=3/2 sector
  // maxRelativeAngularMomentum32 =  maximum relative angular momentum that is used in ThreeBodyPseudoPotential in the S=3/2 sector
  // threeBodyPseudoPotential12 = array with the three-body pseudo-potentials in the S=1/2 sector
  // maxRelativeAngularMomentum12 =  maximum relative angular momentum that is used in ThreeBodyPseudoPotential in the S=1/2 sector
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, 
								 double* threeBodyPseudoPotential32, int maxRelativeAngularMomentum32,
								 double* threeBodyPseudoPotential12, int maxRelativeAngularMomentum12,
								 double alpha, double **pseudoPotential,
								 AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
								 char* precalculationFileName = 0);


  // constructor from datas with a fully-defined two body interaction
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // threeBodyPseudoPotential32 = array with the three-body pseudo-potentials in the S=3/2 sector
  // maxRelativeAngularMomentum32 =  maximum relative angular momentum that is used in ThreeBodyPseudoPotential in the S=3/2 sector
  // threeBodyPseudoPotential12 = array with the three-body pseudo-potentials in the S=1/2 sector
  // maxRelativeAngularMomentum12 =  maximum relative angular momentum that is used in ThreeBodyPseudoPotential in the S=1/2 sector
  // alpha = magnitude of the contact pairing terms
  // pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
  //                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
  // onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
  // onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, 
								 double* threeBodyPseudoPotential32, int maxRelativeAngularMomentum32,
								 double* threeBodyPseudoPotential12, int maxRelativeAngularMomentum12,
								 double alpha, double** pseudoPotential, double* onebodyPotentialUpUp, double* onebodyPotentialDownDown,
								 AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
								 char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();


 protected:
 

  // evaluate additional interaction factors for pairing
  //   
  void EvaluatePairingInteractionFactors();

    // compute all projector coefficient associated to a given relative angular momentum between 3 particles
  //
  // relativeMomentum = value of twice the relative angular momentum between the 3 particles
  // degeneracyIndex = optional degeneracy index for relative angular momentum greater than 5 for bosons (8 for fermions)
  // indices = array that contains all possible sets of indices (size of the array is 3 * nbrIndexSets)
  // nbrIndexSets = number of sets
  // maxJValue = twice the maximum total angular momentum two particles can have
  // spinIndex = indicate for which of three body operators coeeficients are computed (0 for up-up-up, 1 for up-up-down, 2 for up-down-up, 3 for down-up-up) 
  double* ComputeProjectorCoefficients(int relativeMomentum, int degeneracyIndex, int* indices, int nbrIndexSets, int maxJValue, int spinIndex);

  // compute all projector coefficient associated to a given relative angular momentum between 3 particles in the spin 1/2 sector along z_1^u + z_2^u - 2 z_3^d
  //
  // relativeMomentum = value of twice the relative angular momentum between the 3 particles
  // degeneracyIndex = optional degeneracy index for relative angular momentum greater than 5 for bosons (8 for fermions)
  // indices = array that contains all possible sets of indices (size of the array is 3 * nbrIndexSets)
  // nbrIndexSets = number of sets
  // maxJValue = twice the maximum total angular momentum two particles can have
  // spinIndex = indicate for which of three body operators coeeficients are computed (0 for up-up-up, 1 for up-up-down, 2 for up-down-up, 3 for down-up-up) 
  double* ComputeProjectorCoefficientsSpin12E112(int relativeMomentum, int degeneracyIndex, int* indices, int nbrIndexSets, int maxJValue, int spinIndex);

  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();




};

#endif


