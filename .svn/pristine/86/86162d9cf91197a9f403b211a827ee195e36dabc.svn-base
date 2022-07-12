////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of fermions on sphere using                     //
//                            the Lz <-> -Lz symmetry                         //
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


#ifndef FERMIONONSPHERESYMMETRICBASISLONG_H
#define FERMIONONSPHERESYMMETRICBASISLONG_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereLong.h"

#include <iostream>



class FermionOnSphereSymmetricBasisLong :  public FermionOnSphereLong
{

  friend class FermionOnSphereHaldaneSymmetricBasisLong;

 protected:

  // shift to apply to a state before inverting its expression
  int InvertShift;
  // shift to apply to a state after inverting its expression
  int InvertUnshift;


  // signature associated to temporary state used when applying ProdA operator
  ULONGLONG ProdASignature;

 public:

  // default constuctor
  //
  FermionOnSphereSymmetricBasisLong();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // lzMax = maximum Lz value reached by a fermion
  // memory = amount of memory granted for precalculations
  FermionOnSphereSymmetricBasisLong (int nbrFermions, int lzMax, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereSymmetricBasisLong (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereSymmetricBasisLong(const FermionOnSphereSymmetricBasisLong& fermions);

  // destructor
  //
  virtual ~FermionOnSphereSymmetricBasisLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereSymmetricBasisLong& operator = (const FermionOnSphereSymmetricBasisLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter);

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphere* targetSpace);

  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

  // get information about any additional symmetry of the Hilbert space
  //
  // return value = symmetry id
  virtual int GetHilbertSpaceAdditionalSymmetry();

  // convert a given state from Lz-symmetric basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphereLong& nbodyBasis);

  // convert a given state from the usual n-body basis to the Lz-symmetric basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereLong& nbodyBasis);

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
  //
  // index = index of the state on which the operator has to be applied
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient);

  // apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int n1, int n2);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int nbrIndices);

  // apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int nbrIndices, double& coefficient);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdA (int index, int m);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdA (int index, int m, int n, double& coefficient);

  // evaluate wave function in real space using a given basis and only for agiven range of components
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
					int firstComponent, int nbrComponent);                                
  
  // initialize evaluation of wave function in real space using a given basis and only for a given range of components and
  //
  // timeCoherence = true if time coherence has to be used
  virtual void InitializeWaveFunctionEvaluation (bool timeCoherence = false);
  
  // compute particule-hole symmetric state from a given state
  //
  // state = vector corresponding to the state to symmetrize
  // holeBasis = n-body basis on which the symmetrized state has to be expressed
  virtual RealVector ParticleHoleSymmetrize (RealVector& state, FermionOnSphereLong& holeBasis);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(ULONGLONG stateDescription, int lzmax);

  // get symmetry of a given state 
  //
  // initialState = referennce state whose symmetry has to be computed
  virtual void GetStateSymmetry (ULONGLONG& initialState);

  // get canonical expression of a given state
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state
  virtual ULONGLONG GetCanonicalState (ULONGLONG initialState);

  // get canonical expression of a given state and its symmetry
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state (with symmetry bit)
  virtual ULONGLONG GetSignedCanonicalState (ULONGLONG initialState);

};

// get symmetry of a given state 
//
// initialState = reference on the state whose symmetry has to be computed

inline void FermionOnSphereSymmetricBasisLong::GetStateSymmetry (ULONGLONG& initialState)
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

// get canonical expression of a given state
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state

inline ULONGLONG FermionOnSphereSymmetricBasisLong::GetCanonicalState (ULONGLONG initialState)
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

// get canonical expression of a given state and its symmetry
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state (with symmetry bit)

inline ULONGLONG FermionOnSphereSymmetricBasisLong::GetSignedCanonicalState (ULONGLONG initialState)
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

inline int FermionOnSphereSymmetricBasisLong::FindStateIndex(ULONGLONG stateDescription, int lzmax)
{
  stateDescription &= FERMION_SPHERE_LONG_SYMMETRIC_MASK;
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
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
    return PosMin;
}


#endif


 
