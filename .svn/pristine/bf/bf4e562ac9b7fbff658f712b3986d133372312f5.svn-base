////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                           Class author : Cecile Repellin                   //
//                                                                            //
//    class of hamiltonian associated to particles on a CP2 sphere with       //
//                            delta 3-body interaction                        //
//                                                                            //
//                        last modification : 09/01/2013                      //
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


#ifndef PARTICLEONCP2THREEBODYDELTAHAMILTONIAN_H
#define PARTICLEONCP2THREEBODYDELTAHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/BosonOnCP2.h"
#include "Hamiltonian/AbstractQHEOnSphereNBodyInteractionHamiltonian.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients;


class ParticleOnCP2ThreeBodyDeltaHamiltonian : public AbstractQHEOnSphereNBodyInteractionHamiltonian
{

 protected:

   // array with the onebody potentials 
  double* OneBodyPotentials;

  // flag whether to use normalized n-body eigenstates as the reference potentials
  bool NormalizeFlag;
  
  // number of flux quanta
  int NbrFluxQuanta;
  
 public:

  // default constructor
  //
  ParticleOnCP2ThreeBodyDeltaHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // onebodyPotential = array of one-body potentials (optional)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  // normalizePPs = use normalized n-body eigenstates as the reference potentials
  ParticleOnCP2ThreeBodyDeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta, double* onebodyPotentials,
											 AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
											 char* precalculationFileName);

  // destructor
  //
  ~ParticleOnCP2ThreeBodyDeltaHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();


 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();
   
  // get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices
  //
  // nbrValues = number of different values an index can have
  // nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
  // sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
  // return value = total number of index groups

  void GetAllIndices (int nbrValues, int nbrFluxQuanta, int nbrSectorSum, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum, double**& interactionCoef);

 

};

#endif
