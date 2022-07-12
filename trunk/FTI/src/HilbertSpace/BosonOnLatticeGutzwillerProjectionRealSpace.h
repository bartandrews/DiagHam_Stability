////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of hardcore boson on lattice in real space        //
//                                                                            //
//                        last modification : 10/09/2014                      //
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


#ifndef BOSONONLATTICEGUTZWILLERPROJECTIONREALSPACE_H
#define BOSONONLATTICEGUTZWILLERPROJECTIONREALSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"

#include <iostream>



class BosonOnLatticeGutzwillerProjectionRealSpace : public FermionOnLatticeRealSpace
{

  friend class BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation;

 protected:

  // number of bosons
  int NbrBosons;
  int IncNbrBosons;

  // indices of the orbitals that are kept when performing an orbital cut
  int* KeptOrbitals;

 public:

  // default constructor
  // 
  BosonOnLatticeGutzwillerProjectionRealSpace ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // memory = amount of memory granted for precalculations
  BosonOnLatticeGutzwillerProjectionRealSpace (int nbrBosons, int nbrSite, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  BosonOnLatticeGutzwillerProjectionRealSpace(const BosonOnLatticeGutzwillerProjectionRealSpace& bosons);

  // destructor
  //
  ~BosonOnLatticeGutzwillerProjectionRealSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnLatticeGutzwillerProjectionRealSpace & operator = (const BosonOnLatticeGutzwillerProjectionRealSpace& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // evaluate the orbital cut entanglement matrix. The entanglement matrix is only evaluated for fixed number of particles
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // keptOrbitals = array of orbitals that have to be kept, should be sorted from the smallest index to the largest index 
  // nbrKeptOrbitals = array of orbitals that have to be kept
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrParticleSector, int nbrKeptOrbitals, int* keptOrbitals, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  virtual void GetCompositeFermionWavefunction(ComplexVector & trialState, ComplexMatrix & jastrowEigenVecs,ComplexMatrix & cFEigenVecs);

 protected:

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in complementary Hilbert space
  // nbrIndex = number of indices to consider in complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  
								  ParticleOnSphere* destinationHilbertSpace,
								  ComplexVector& groundState,  HermitianMatrix* densityMatrix);
    
  // apply a_n  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A (int index, int n);

  // apply a^+_n  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator
  // return value =  multiplicative factor 
  virtual double Ad (int index, int n);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual long AdA (long index, int m, int n, double& coefficient);


  // apply a^+_m operator to the state produced using the A or Ad method (without destroying it)
  //
  // m = index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad (int m, double& coefficient);
  
  // apply a_n  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = index for first annihilation operator
  // n2 = index for first annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int n1, int n2);
  
  // apply a^+_m operator to the state produced using the A or Ad method (without destroying it)
  //
  // m = index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient);


  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();
};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnLatticeGutzwillerProjectionRealSpace::GetParticleStatistic()
{
  return ParticleOnSphere::BosonicStatistic;
}

#endif


