////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of fermions on lattice in real space              //
//                                                                            //
//                        last modification : 09/09/2014                      //
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


#ifndef FERMIONONLATTICEREALSPACE_H
#define FERMIONONLATTICEREALSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnSphere.h"

#include <iostream>



class FermionOnLatticeRealSpace : public FermionOnSphere
{

 friend class FermionOnLatticeRealSpaceAnd2DTranslation;
 friend class FermionOnLatticeWithSpinRealSpaceAnd2DTranslation;
 friend class BosonOnLatticeRealSpaceAnd2DTranslation;

 protected:

  // total number of sites
  int NbrSite;
  
  // indices of the orbitals that are kept when performing an orbital cut
  int* KeptOrbitals;

 public:

  // default constructor
  // 
  FermionOnLatticeRealSpace ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // memory = amount of memory granted for precalculations
  FermionOnLatticeRealSpace (int nbrFermions, int nbrSite, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeRealSpace(const FermionOnLatticeRealSpace& fermions);

  // destructor
  //
  ~FermionOnLatticeRealSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeRealSpace& operator = (const FermionOnLatticeRealSpace& fermions);

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

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  //virtual int FindStateIndex(unsigned long stateDescription, int lzmax);
  
  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, long pos);

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
  virtual long EvaluatePartialEntanglementMatrixZeroParticuleCase (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace, ComplexVector& groundState, ComplexMatrix* entanglementMatrix);
  virtual long EvaluatePartialEntanglementMatrixMaxParticuleCase (ParticleOnSphere* destinationHilbertSpace, ComplexVector& groundState, ComplexMatrix* entanglementMatrix);

};


#endif


