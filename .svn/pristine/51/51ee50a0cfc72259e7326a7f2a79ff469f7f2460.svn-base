////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated to particles on a disk with         //
//                          generic 3-body interaction                        //
//                                                                            //
//                        last modification : 16/09/2010                      //
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


#ifndef PARTICLEONDISKGENERICTHREEBODYHAMILTONIAN_H
#define PARTICLEONDISKGENERICTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereNBodyInteractionHamiltonian.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients;


class ParticleOnDiskGenericThreeBodyHamiltonian : public AbstractQHEOnSphereNBodyInteractionHamiltonian
{

 protected:

  
  // array with the three-body pseudo-potentials sorted with respect to the relative angular momentum, taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  double* ThreeBodyPseudoPotential;
  // nuber of elements in the ThreeBodyPseudoPotential array
  int NbrThreeBodyPseudoPotential;
  // maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential
  int MaxRelativeAngularMomentum;

  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  double* PseudoPotential;

  // array with the onebody potentials 
  double* OneBodyPotentials;

  // flag whether to use normalized n-body eigenstates as the reference potentials
  bool NormalizeFlag;


 public:

  // default constructor
  //
  ParticleOnDiskGenericThreeBodyHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // threeBodyPseudoPotential = array with the three-body pseudo-potentials sorted with respect to the relative angular momentum, 
  //                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  // maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential
  // pseudoPotential = array of two-body pseudopotentials (optional)
  // onebodyPotential = array of one-body potentials (optional)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  // normalizePPs = use normalized n-body eigenstates as the reference potentials
  ParticleOnDiskGenericThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
					    double* threeBodyPseudoPotential, int maxRelativeAngularMomentum, double* pseudoPotential, double* onebodyPotentials,
					    AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
					    char* precalculationFileName = 0, bool normalizePPs = false);

  // destructor
  //
  ~ParticleOnDiskGenericThreeBodyHamiltonian();

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
  double* ComputeProjectorCoefficients(int relativeMomentum, int degeneracyIndex, int* indices, int nbrIndexSets);

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();


};

#endif
