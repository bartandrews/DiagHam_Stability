////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//        that has to handle (temporary) Hilbert space larger than 2^31       //
//                                                                            //
//                        last modification : 14/05/2009                      //
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


#ifndef FERMIONONSPHEREHALDANELARGEBASIS_H
#define FERMIONONSPHEREHALDANELARGEBASIS_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include <iostream>


class FermionOnSphereHaldaneLargeBasis :  public FermionOnSphereHaldaneBasis
{

 protected:
  
  // depth of the B-tree used to improve state index search speed
  int BTreeLookUpTableDepth;
  // B-tree used to improve state index search speed
  unsigned long* BTreeLookUpTable;
  // number of state indices covered by a leaf of the B-tree
  long BTreeLookUpTableNbrIndices;

 public:

  // default constructor
  //
  FermionOnSphereHaldaneLargeBasis();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = reference on twice the momentum total value, totalLz will be recomputed from referenceState and stored in totalLz
  // lzMax = maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  // referenceState = array that describes the reference state to start from
  FermionOnSphereHaldaneLargeBasis (int nbrFermions, int& totalLz, int lzMax, int* referenceState,
			       unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereHaldaneLargeBasis (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereHaldaneLargeBasis(const FermionOnSphereHaldaneLargeBasis& fermions);

  // destructor
  //
  virtual ~FermionOnSphereHaldaneLargeBasis ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereHaldaneLargeBasis& operator = (const FermionOnSphereHaldaneLargeBasis& fermions);

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // convert a given state from Haldane basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis);

  // convert a given state from the usual n-body basis to the Haldane basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertFromNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis);

  // create the Jack polynomial decomposition corresponding to the root partition
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  virtual RealVector& GenerateJackPolynomial(RealVector& jack, double alpha);

  // create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  virtual RealVector& GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
  // 
  // subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
  // nbrFermionSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  RealSymmetricMatrix EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // return value = corresponding index
  virtual long FindStateIndex(unsigned long stateDescription);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int lzMax, unsigned long referenceState, long pos, long& memory);

  // generate all states (i.e. all possible skew symmetric polynomials with fixed Lz)
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, long pos);

};


// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

inline int FermionOnSphereHaldaneLargeBasis::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  return ((int) this->FindStateIndex(stateDescription));
}

#endif


