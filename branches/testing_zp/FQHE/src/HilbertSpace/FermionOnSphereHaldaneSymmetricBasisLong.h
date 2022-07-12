////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//                           with Lz <-> -Lz symmetry                         //
//   that allow LzMax up to 127 (for systems with 128 bit integer support)    //
//               or 63 (on 32 bit systems without 128 bit integer support)    //
//                                                                            //
//                        last modification : 23/09/2007                      //
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


#ifndef FERMIONONSPHEREHALDANESYMMETRYBASISLONG_H
#define FERMIONONSPHEREHALDANESYMMETRYBASISLONG_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"

#include <iostream>


class FermionOnSphereHaldaneSymmetricBasisLong :  public FermionOnSphereSymmetricBasisLong
{

 protected:

  // topmost state 
  ULONGLONG ReferenceState;

  // three temporary arrays used during Hilbert space generation
  ULONGLONG* TmpGeneratedStates;
  int* TmpGeneratedStatesLzMax;
  unsigned long* KeepStateFlag;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // lzMax = maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  // referenceState = array that describes the reference state to start from
  FermionOnSphereHaldaneSymmetricBasisLong (int nbrFermions, int lzMax, int* referenceState,
					    unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereHaldaneSymmetricBasisLong (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereHaldaneSymmetricBasisLong(const FermionOnSphereHaldaneSymmetricBasisLong& fermions);

  // destructor
  //
  virtual ~FermionOnSphereHaldaneSymmetricBasisLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereHaldaneSymmetricBasisLong& operator = (const FermionOnSphereHaldaneSymmetricBasisLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // convert a gien state from Haldane basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphereLong& nbodyBasis);

  // convert a gien state from Haldane basis to the usual symmetric n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the  symmetric nbody-basis to use
  // return value = converted vector
  RealVector ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereSymmetricBasisLong& nbodyBasis);

  // convert a gien state from Lz-symmetric Haldane basis to the usual Haldane n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToHaldaneNbodyBasis(RealVector& state, FermionOnSphereHaldaneBasisLong& nbodyBasis);

  // convert a given state from the usual n-body Haldane basis to the Lz-symmetric Haldane basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToSymmetricHaldaneNbodyBasis(RealVector& state, FermionOnSphereLong& nbodyBasis);


 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(ULONGLONG stateDescription, int lzmax);

  // generate all states corresponding to the constraints
  // 
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual int GenerateStates(int lzMax, ULONGLONG referenceState, int pos, long& memory);

  // generate all states (i.e. all possible skew symmetric polynomials with fixed Lz)
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual int RawGenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int pos);

  // get canonical expression of a given state
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state
  virtual ULONGLONG GetCanonicalState (ULONGLONG initialState);

  // get symmetry of a given state 
  //
  // initialState = referennce state whose symmetry has to be computed
  virtual void GetStateSymmetry (ULONGLONG& initialState);

  // get canonical expression of a given state and its symmetry
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state (with symmetry bit)
  virtual ULONGLONG GetSignedCanonicalState (ULONGLONG initialState);

};


// get canonical expression of a given state
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state

inline ULONGLONG FermionOnSphereHaldaneSymmetricBasisLong::GetCanonicalState (ULONGLONG initialState)
{
  initialState <<= this->InvertShift;
#ifdef __128_BIT_LONGLONG__
  ULONGLONG TmpState = FermionOnSphereSymmetricBasisInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 120;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 112;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 104;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 96;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 88;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 80;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 72;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 56) & ((ULONGLONG) 0xfful)] << 64; 
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 64) & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 72) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 80) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 88) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 96) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 104) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 112) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[initialState >> 120]; 
#else
  ULONGLONG TmpState = FermionOnSphereSymmetricBasisInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[initialState >> 56]; 
#endif
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
  if (TmpState < initialState)
    return TmpState;
  else
    return initialState;
}

// get symmetry of a given state 
//
// initialState = reference on the state whose symmetry has to be computed

inline void FermionOnSphereHaldaneSymmetricBasisLong::GetStateSymmetry (ULONGLONG& initialState)
{
  initialState <<= this->InvertShift;
#ifdef __128_BIT_LONGLONG__
  ULONGLONG TmpState = FermionOnSphereSymmetricBasisInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 120;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 112;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 104;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 96;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 88;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 80;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 72;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 56) & ((ULONGLONG) 0xfful)] << 64; 
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 64) & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 72) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 80) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 88) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 96) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 104) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 112) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[initialState >> 120]; 
#else
  ULONGLONG TmpState = FermionOnSphereSymmetricBasisInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[initialState >> 56];  
#endif
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
  if (TmpState != initialState)    
    initialState |= FERMION_SPHERE_LONG_SYMMETRIC_BIT;
}

// get canonical expression of a given state and its symmetry
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state (with symmetry bit)

inline ULONGLONG FermionOnSphereHaldaneSymmetricBasisLong::GetSignedCanonicalState (ULONGLONG initialState)
{
  initialState <<= this->InvertShift;
#ifdef __128_BIT_LONGLONG__
  ULONGLONG TmpState = FermionOnSphereSymmetricBasisInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 120;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 112;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 104;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 96;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 88;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 80;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 72;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 56) & ((ULONGLONG) 0xfful)] << 64; 
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 64) & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 72) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 80) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 88) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 96) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 104) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 112) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[initialState >> 120]; 
#else
  ULONGLONG TmpState = FermionOnSphereSymmetricBasisInvertTableLong[initialState & ((ULONGLONG) 0xfful)] << 56;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 8) & ((ULONGLONG) 0xfful)] << 48;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 16) & ((ULONGLONG) 0xfful)] << 40;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 24) & ((ULONGLONG) 0xfful)] << 32;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 32) & ((ULONGLONG) 0xfful)] << 24;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 40) & ((ULONGLONG) 0xfful)] << 16;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[(initialState >> 48) & ((ULONGLONG) 0xfful)] << 8;
  TmpState |= FermionOnSphereSymmetricBasisInvertTableLong[initialState >> 56];  
#endif
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
  if (TmpState < initialState)
    {
      return (TmpState | FERMION_SPHERE_LONG_SYMMETRIC_BIT);
    }
  else
    if (TmpState != initialState)
      return (initialState | FERMION_SPHERE_LONG_SYMMETRIC_BIT);
    else
      return initialState;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

inline int FermionOnSphereHaldaneSymmetricBasisLong::FindStateIndex(ULONGLONG stateDescription, int lzmax)
{
  stateDescription &= FERMION_SPHERE_LONG_SYMMETRIC_MASK;
  long PosMax = (long) (stateDescription >> this->LookUpTableShift[lzmax]);
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  ULONGLONG CurrentState = (this->StateDescription[PosMid] & FERMION_SPHERE_LONG_SYMMETRIC_MASK);
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = (this->StateDescription[PosMid] & FERMION_SPHERE_LONG_SYMMETRIC_MASK);
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}



#endif


