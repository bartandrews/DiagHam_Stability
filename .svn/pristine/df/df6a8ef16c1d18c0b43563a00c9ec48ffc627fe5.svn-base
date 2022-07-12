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


#ifndef PARTICLEONSPHEREWITHSPINBASICTHREEBODYHAMILTONIAN_H
#define PARTICLEONSPHEREWITHSPINBASICTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients;


class ParticleOnSphereWithSpinBasicThreeBodyHamiltonian : public AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian
{

 protected:

  // array with the three-body pseudo-potentials, sorted with respect to the relative angular momentu (second index)
  // taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  // first index stands for the spin sector (0 up-up-up, 1 down-down-down, 2 up-up-down, 3 up-down-down)
  double** ThreeBodyPseudoPotentials;
  // number of elements in the ThreeBodyPseudoPotential array for each spin sector
  int* NbrThreeBodyPseudoPotentials;
  // maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential
  int* MaxRelativeAngularMomentum;

  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  // first index refered to the spin sector (sorted as up-up, down-down, up-down)
  double** PseudoPotentials;

 public:

  // default constructor
  //
  ParticleOnSphereWithSpinBasicThreeBodyHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // threeBodyPseudoPotential = array with the three-body pseudo-potentials sorted with respect to the relative angular momentum,
  //                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  //                            first index is the spin sector (0 up-up-up, 1 down-down-down, 2 up-up-down, 3 up-down-down)
  // maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential  for each spin sector
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereWithSpinBasicThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax,
                                                      double** threeBodyPseudoPotential, int* maxRelativeAngularMomentum,
                                                      AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false,
                                                      char* precalculationFileName = 0);

  // constructor from datas with a fully-defined two body interaction
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // threeBodyPseudoPotential = array with the three-body pseudo-potentials sorted with respect to the relative angular momentum,
  //                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  //                            first index is the spin sector (0 up-up-up, 1 down-down-down, 2 up-up-down, 3 up-down-down)
  // maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential  for each spin sector
  // pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
  //                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
  // onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
  // onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereWithSpinBasicThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax,
                                                      double** threeBodyPseudoPotential, int* maxRelativeAngularMomentum,
                                                      double** pseudoPotential, double* onebodyPotentialUpUp, double* onebodyPotentialDownDown,
                                                      AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false,
                                                      char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereWithSpinBasicThreeBodyHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();


 protected:
 
  // compute all projector coefficient associated to a given relative angular momentum between 3 particles
  //
  // relativeMomentum = value of twice the relative angular momentum between the 3 particles
  // degeneracyIndex = optional degeneracy index for relative angular momentum greater than 5 for bosons (8 for fermions)
  // indices = array that contains all possible sets of indices (size of the array is 3 * nbrIndexSets)
  // nbrIndexSets = number of sets
  // maxJValue = twice the maximum total angular momentum two particles can have
  // spinIndex = indicate for which of three body operators coeeficients are computed (0 for up-up-up, 1 for up-up-down, 2 for up-down-up, 3 for down-up-up)
  double* ComputeProjectorCoefficients(int relativeMomentum, int degeneracyIndex, int* indices, int nbrIndexSets, int maxJValue, int spinIndex);

  // evaluate all interaction factors
  //  
  void EvaluateInteractionFactors();


};

#endif
