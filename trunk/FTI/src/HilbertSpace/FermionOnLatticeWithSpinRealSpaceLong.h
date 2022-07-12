////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//                  in real space supporting up to 64 sites                   //
//                                                                            //
//                       class author: Cecile Repellin                        //
//                                                                            //
//                        last modification : 15/07/2016                      //
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


#ifndef FERMIONONLATTICEWITHSPINREALSPACELONG_H
#define FERMIONONLATTICEWITHSPINREALSPACELONG_H

#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"

#include <iostream>



class FermionOnLatticeWithSpinRealSpaceLong : public FermionOnSphereWithSpinLong
{

  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;
  friend class FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong;
  friend class FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong;

 protected:

  // total number of sites
  int NbrSite;
  
  // flag to indicate that the Hilbert space should preserve Sz
  bool SzFlag;

  // indices of the orbitals that are kept when performing an orbital cut
  int* KeptOrbitals;

 public:

  // default constructor
  // 
  FermionOnLatticeWithSpinRealSpaceLong ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinRealSpaceLong (int nbrFermions, int nbrSite, unsigned long memory = 10000000);

  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the total spin value
  // nbrSite = number of sites in the x direction
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinRealSpaceLong (int nbrFermions, int totalSpin, int nbrSite, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinRealSpaceLong(const FermionOnLatticeWithSpinRealSpaceLong& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinRealSpaceLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinRealSpaceLong& operator = (const FermionOnLatticeWithSpinRealSpaceLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // find state index from a string
  //
  // stateDescription = string describing the state
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(char* stateDescription);

  // evaluate the orbital cut entanglement matrix. The entanglement matrix is only evaluated for fixed number of particles
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // keptOrbitals = array of orbitals that have to be kept, should be sorted from the smallest index to the largest index 
  // nbrKeptOrbitals = array of orbitals that have to be kept
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrParticleSector, int nbrKeptOrbitals, int* keptOrbitals, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
  
  // evaluate the orbital cut entanglement matrix. The entanglement matrix is only evaluated for fixed number of particles and Sz
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // szSector  = twice the total Sz value of the subsytem 
  // groundState = reference on the total system ground state
  // keptOrbitals = array of orbitals that have to be kept, should be sorted from the smallest index to the largest index 
  // nbrKeptOrbitals = array of orbitals that have to be kept
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrParticleSector, int szSector, int nbrKeptOrbitals, int* keptOrbitals, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
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
  
 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(ULONGLONG stateDescription, int lzmax);
  
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
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real space
  // nbrSpinUp = number of fermions with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos);

  // core part of the evaluation orbital cut entanglement matrix calculation
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialEntanglementMatrixCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
						      ComplexVector& groundState, ComplexMatrix* entanglementMatrix);

};


#endif


