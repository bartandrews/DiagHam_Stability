////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of fermions on lattice with spin  and Gutzwiller          //
//                          projection in real space                          //
//                                                                            //
//                       class author: Nicolas Regnault                       //
//                                                                            //
//                        last modification : 17/06/2014                      //
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


#ifndef FERMIONONLATTICEWITHSPINANDGUTZWILLERPROJECTIONREALSPACE_H
#define FERMIONONLATTICEWITHSPINANDGUTZWILLERPROJECTIONREALSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"

#include <iostream>



class FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace : public FermionOnLatticeWithSpinRealSpace
{

  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;
  friend class FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation;

 protected:


 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (int nbrFermions, int nbrSite, unsigned long memory = 10000000);

  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the total spin value
  // nbrSite = number of sites in the x direction
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (int nbrFermions, int totalSpin, int nbrSite, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace(const FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace& operator = (const FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  virtual RealVector ForgeSU2FromU1(RealVector& upState, ParticleOnSphere* upStateSpace, RealVector& downState, ParticleOnSphere* downStateSpace);

  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  virtual ComplexVector ForgeSU2FromU1(ComplexVector& upState, ParticleOnSphere* upStateSpace, ComplexVector& downState, ParticleOnSphere* downStateSpace);

  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  virtual RealVector ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace);

  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  virtual ComplexVector ForgeSU2FromU1(ComplexVector& upState, FermionOnSphere& upStateSpace, ComplexVector& downState, FermionOnSphere& downStateSpace);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = kx sector in which the density matrix has to be evaluated 
  // kySector = kx sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in given momentum and Sz sectors.
  //
  // nbrParticleSector = number of particles that belong to the subsytem
  // szSector  = twice the total Sz value of the subsytem
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
  
  // evaluate the orbital cut entanglement matrix. The entanglement matrix is only evaluated for fixed number of particles
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // keptOrbitals = array of orbitals that have to be kept, should be sorted from the smallest index to the largest index 
  // nbrKeptOrbitals = array of orbitals that have to be kept
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrParticleSector, int nbrKeptOrbitals, int* keptOrbitals, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
  

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions);

  // evaluate Hilbert space dimension with a fixed number of fermions with spin up
  //
  // nbrFermions = number of fermions
  // nbrSpinUp = number of fermions with spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // nbrHoles = number of unoccupied sites
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, int nbrHoles, long pos);
  
  // generate all states corresponding to the constraints with a fixed number of fermions with spin up
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // nbrHoles = number of unoccupied sites
  //nbrSpinUp = number of fermions with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, int nbrHoles, int nbrSpinUp, long pos);
  
  // find state index (and checks state belongs to Hilbert space)
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);
  
  // find state index and checks that states belongs to Hilbert space
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int CarefulFindStateIndex(unsigned long stateDescription, int lzmax);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real space
  // nbrSpinUp = number of fermions with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
//   virtual long GenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
								  ComplexVector& groundState, HermitianMatrix* densityMatrix);


};


// find state index and checks that states belongs to Hilbert space
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index
inline int FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::CarefulFindStateIndex(unsigned long stateDescription, int lzmax)
{
 return this->FindStateIndex(stateDescription, lzmax); 
}

#endif


