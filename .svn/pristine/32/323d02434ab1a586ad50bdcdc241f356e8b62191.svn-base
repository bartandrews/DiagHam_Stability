////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                          generic 4-body interaction                        //
//                                                                            //
//                        last modification : 12/08/2009                      //
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


#ifndef PARTICLEONSPHEREGENERICFOURBODYHAMILTONIAN_H
#define PARTICLEONSPHEREGENERICFOURBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/AbstractQHEOnSphereNBodyInteractionHamiltonian.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients;


class ParticleOnSphereGenericFourBodyHamiltonian : public AbstractQHEOnSphereNBodyInteractionHamiltonian
{

 protected:

  
  // array with the four-body pseudo-potentials sorted with respect to the relative angular momentum, taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  double* FourBodyPseudoPotential;
  // nuber of elements in the FourBodyPseudoPotential array
  int NbrFourBodyPseudoPotential;
  // maxixmum relative angular momentum that is used in FourBodyPseudoPotential
  int MaxRelativeAngularMomentum;

  // array with the three-body pseudo-potentials sorted with respect to the relative angular momentum, taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  double* ThreeBodyPseudoPotential;
  // nuber of elements in the ThreeBodyPseudoPotential array
  int NbrThreeBodyPseudoPotential;
  // maxixmum relative angular momentum that is used in FourBodyPseudoPotential
  int ThreeBodyMaxRelativeAngularMomentum;

  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  double* PseudoPotential;

  // flag whether to use normalized n-body eigenstates as the reference potentials
  bool NormalizeFlag;

 public:

  // default constructor
  //
  ParticleOnSphereGenericFourBodyHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // fourBodyPseudoPotential = array with the four-body pseudo-potentials sorted with respect to the relative angular momentum, 
  //                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  // maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FourBodyPseudoPotential
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereGenericFourBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
					      double* fourBodyPseudoPotential, int maxRelativeAngularMomentum, double l2Factor, 
					      AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
					      char* precalculationFileName = 0);

  // constructor from default datas with three-body interaction
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // fourBodyPseudoPotential = array with the four-body pseudo-potentials sorted with respect to the relative angular momentum, 
  //                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  // maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FourBodyPseudoPotential
// threeBodyPseudoPotential = array with the three-body pseudo-potentials sorted with respect to the relative angular momentum, 
//                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
// threeBodymaxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereGenericFourBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
					      double* fourBodyPseudoPotential, int maxRelativeAngularMomentum, double* threeBodyPseudoPotential, int threeBodyMaxRelativeAngularMomentum, double l2Factor, 
					      AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
					      char* precalculationFileName = 0);

  // constructor from datas with a fully-defined two body interaction
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // fourBodyPseudoPotential = array with the four-body pseudo-potentials sorted with respect to the relative angular momentum, 
  //                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  // maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FourBodyPseudoPotential
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereGenericFourBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
					      double* fourBodyPseudoPotential, int maxRelativeAngularMomentum,
					      double l2Factor, double* pseudoPotential, 
					      AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
					      char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereGenericFourBodyHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();


 protected:
 
  // compute all projector coefficient associated to a given relative angular momentum between 4 particles
  //
  // relativeMomentum = value of twice the relative angular momentum between the 4 particles
  // indices = array that contains all possible sets of indices (size of the array is 4 * nbrIndexSets)
  // nbrIndexSets = number of sets
  double* Compute4BodyCoefficients(int relativeMomentum, int* indices, int nbrIndexSets);

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // evaluate all interaction factors for the four body interaction part
  //   
  void Evaluate4BodyInteractionFactors();

  // evaluate all interaction factors for 3body
  //   
  virtual void Evaluate3BodyInteractionFactors();

  // compute a given projector coefficient for the 4-body interaction 
  //
  // m1 = first index
  // m2 = second index
  // m3 = third inde
  // m4 = fourth index
  // jValue = total angular momentum
  // minJ = minimum angular momentum that can be reach by three particles
  // return value = corresponding projector coefficient
  double ComputeProjectorCoefficients4Body(int m1, int m2, int m3, int m4, int jValue, 
					   ClebschGordanCoefficients* clebshArray);

  // compute a given projector coefficient for the 4-body interaction 
  //
  // m1 = first index
  // m2 = second index
  // m3 = third inde
  // m4 = fourth index
  // jValue = total angular momentum
  // minJ = minimum angular momentum that can be reach by three particles
  // return value = corresponding projector coefficient
  double ComputeProjectorCoefficients4Body2(int m1, int m2, int m3, int m4, int jValue, 
					    ClebschGordanCoefficients* clebshArray);

  // compute all projector coefficient associated to a given relative angular momentum between 3 particles
  //
  // relativeMomentum = value of twice the relative angular momentum between the 3 particles
  // degeneracyIndex = optional degeneracy index for relative angular momentum greater than 5 for bosons (8 for fermions)
  // indices = array that contains all possible sets of indices (size of the array is 3 * nbrIndexSets)
  // nbrIndexSets = number of sets
  double* ComputeProjectorCoefficients3Body(int relativeMomentum, int degeneracyIndex, int* indices, int nbrIndexSets);
};

#endif
