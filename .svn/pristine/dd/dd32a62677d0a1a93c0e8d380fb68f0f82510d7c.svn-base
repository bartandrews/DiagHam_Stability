////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                          generic 5-body interaction                        //
//                                                                            //
//                        last modification : 29/08/2009                      //
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


#ifndef PARTICLEONSPHEREGENERICFIVEBODYHAMILTONIAN_H
#define PARTICLEONSPHEREGENERICFIVEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnSphereGenericFourBodyHamiltonian.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients;


class ParticleOnSphereGenericFiveBodyHamiltonian : public ParticleOnSphereGenericFourBodyHamiltonian
{

 protected:

  
  // array with the five-body pseudo-potentials sorted with respect to the relative angular momentum, taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  double* FiveBodyPseudoPotential;
  // nuber of elements in the FiveBodyPseudoPotential array
  int NbrFiveBodyPseudoPotential;
  // maxixmum relative angular momentum that is used in FiveBodyPseudoPotential
  int MaxRelativeAngularMomentum;

 public:

  // default constructor
  //
  ParticleOnSphereGenericFiveBodyHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // fiveBodyPseudoPotential = array with the five-body pseudo-potentials sorted with respect to the relative angular momentum, 
  //                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  // maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FiveBodyPseudoPotential
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereGenericFiveBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
					      double* fiveBodyPseudoPotential, int maxRelativeAngularMomentum, double l2Factor, 
					      AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
					      char* precalculationFileName = 0);

 // constructor from default datas with a full four-body interaction
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // fiveBodyPseudoPotential = array with the five-body pseudo-potentials sorted with respect to the relative angular momentum, 
  // maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FiveBodyPseudoPotential
  // fourBodyPseudoPotential = array with the four-body pseudo-potentials sorted with respect to the relative angular momentum, 
  // fourBodyMaxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FourBodyPseudoPotential
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereGenericFiveBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
					     double* fiveBodyPseudoPotential, int maxRelativeAngularMomentum, 
					     double* fourBodyPseudoPotential, int fourBodyMaxRelativeAngularMomentum, 
					     double l2Factor, AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
					     char* precalculationFileName);

  // constructor from datas with a fully-defined two body interaction
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // fiveBodyPseudoPotential = array with the five-body pseudo-potentials sorted with respect to the relative angular momentum, 
  //                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  // maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FiveBodyPseudoPotential
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereGenericFiveBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
					      double* fiveBodyPseudoPotential, int maxRelativeAngularMomentum,
					      double l2Factor, double* pseudoPotential, 
					      AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
					      char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereGenericFiveBodyHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();


 protected:
 
  // compute all projector coefficient associated to a given relative angular momentum between 5 particles
  //
  // relativeMomentum = value of twice the relative angular momentum between the 5 particles
  // indices = array that contains all possible sets of indices (size of the array is 5 * nbrIndexSets)
  // nbrIndexSets = number of sets
  double* Compute5BodyCoefficients(int relativeMomentum, int* indices, int nbrIndexSets);

  // compute all projector coefficient associated to a given relative angular momentum between 5 particles in a given direction
  //
  // relativeMomentum = value of twice the relative angular momentum between the 4 particles
  // indices = array that contains all possible sets of indices (size of the array is 4 * nbrIndexSets)
  // nbrIndexSets = number of sets
  // maxClosing = array that gives the maximum angular momentum when a particle approach the cluster of n+1 particles  (n being the index of MaxClosing)
  double* Compute5BodyCoefficientsWithDirection(int relativeMomentum, int* indices, int nbrIndexSets, int* maxClosing);

  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // evaluate all interaction factors for the five body interaction part
  //   
  void Evaluate5BodyInteractionFactors();

  // compute a given projector coefficient for the 5-body interaction 
  //
  // m1 = first index
  // m2 = second index
  // m3 = third inde
  // m4 = fourth index
  // m5 = fifth index
  // jValue = total angular momentum
  // minJ = minimum angular momentum that can be reach by three particles
  // return value = corresponding projector coefficient
  double ComputeProjectorCoefficients5Body(int m1, int m2, int m3, int m4, int m5, int jValue, 
					   ClebschGordanCoefficients* clebshArray);

  // compute a given projector coefficient for the 5-body interaction in a given direction
  //
  // m1 = first index
  // m2 = second index
  // m3 = third inde
  // m4 = fourth index
  // m5 = fifth index
  // jValue = total angular momentum
  // maxClosing = array that gives the maximum angular momentum when a particle approach the cluster of n+1 particles  (n being the index of MaxClosing)
  // return value = corresponding projector coefficient
  double ComputeProjectorCoefficients5BodyWithDirection(int m1, int m2, int m3, int m4, int m5, int jValue, 
							int* maxClosing, ClebschGordanCoefficients* clebshArray);

};

#endif
