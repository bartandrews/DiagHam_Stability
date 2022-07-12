////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of fermions on sphere with spin with             //
//                                 Lz<->-Lz symmetry                          //
//                                                                            //
//                        last modification : 15/08/2007                      //
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


#ifndef FERMIONONSPHEREWITHSPINLZSYMMETRY_H
#define FERMIONONSPHEREWITHSPINLZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"

#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;



class FermionOnSphereWithSpinLzSymmetry :  public FermionOnSphereWithSpinLzSzSymmetry
{

 protected:


 public:

  // default constructor 
  //
  FermionOnSphereWithSpinLzSymmetry ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSz = twce the total spin value
  // minusParity = select the  Lz <-> -Lz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinLzSymmetry (int nbrFermions, int lzMax, int totalSz, bool minusParity, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinLzSymmetry (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinLzSymmetry(const FermionOnSphereWithSpinLzSymmetry& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinLzSymmetry ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinLzSymmetry& operator = (const FermionOnSphereWithSpinLzSymmetry& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // convert a given state from symmetric basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector  
  virtual RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis);

  // convert a given state from the usual n-body basis to the symmetric basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  virtual RealVector ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis);

  // apply a^+_m1_d a^+_m2_d a_n1_d a_n2_d operator to a given state (with m1+m2=n1+n2, only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m1_u a^+_m2_u a_n1_u a_n2_u operator to a given state (with m1+m2=n1+n2, only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_d_m1 a^+_u_m2 a_d_n1 a_u_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

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
  
  protected:

  // get canonical expression of a given state
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state
  virtual unsigned long GetCanonicalState (unsigned long initialState);

  // get symmetry of a given state 
  //
  // initialState = referennce state whose symmetry has to be computed
  virtual void GetStateSymmetry (unsigned long& initialState);

  // get canonical expression of a given state and its symmetry
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state (with symmetry bit)
  virtual unsigned long GetSignedCanonicalState (unsigned long initialState);

  // factorized code that is used to symmetrize result of the AdxAdy operations
  //
  // state = reference on the state that has been produced with the AdxAdy operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state
  virtual int SymmetrizeAdAdResult(unsigned long& state, double& coefficient);

};

// get canonical expression of a given state
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state

inline unsigned long FermionOnSphereWithSpinLzSymmetry::GetCanonicalState (unsigned long initialState)
{
  initialState <<= this->InvertShift;
//  cout << hex << initialState << dec << endl;
#ifdef __64_BITS__
  unsigned long TmpState = FermionOnSphereWithSpinLzInvertTable[initialState & 0xff] << 56;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 8) & 0xff] << 48;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 16) & 0xff] << 40;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 24) & 0xff] << 32;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 32) & 0xff] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 40) & 0xff] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 48) & 0xff] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[initialState >> 56]; 
#else
  unsigned long TmpState = FermionOnSphereWithSpinLzInvertTable[initialState & 0xff] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 8) & 0xff] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 16) & 0xff] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[initialState >> 24];
#endif	
//  cout << hex << TmpState << dec << endl;
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
//  cout << hex << TmpState << dec << endl;
  if (TmpState < initialState)
    return TmpState;
  else
    return initialState;
}

// get symmetry of a given state 
//
// initialState = reference on the state whose symmetry has to be computed

inline void FermionOnSphereWithSpinLzSymmetry::GetStateSymmetry (unsigned long& initialState)
{
  initialState <<= this->InvertShift;
#ifdef __64_BITS__
  unsigned long TmpState = FermionOnSphereWithSpinLzInvertTable[initialState & 0xff] << 56;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 8) & 0xff] << 48;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 16) & 0xff] << 40;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 24) & 0xff] << 32;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 32) & 0xff] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 40) & 0xff] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 48) & 0xff] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[initialState >> 56];  
#else
  unsigned long TmpState = FermionOnSphereWithSpinLzInvertTable[initialState & 0xff] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 8) & 0xff] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 16) & 0xff] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[initialState >> 24];
#endif
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
  if (TmpState != initialState)    
    initialState |= FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT;
}

// get canonical expression of a given state and its symmetry
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state (with symmetry bit)

inline unsigned long FermionOnSphereWithSpinLzSymmetry::GetSignedCanonicalState (unsigned long initialState)
{
  initialState <<= this->InvertShift;
#ifdef __64_BITS__
  unsigned long TmpState = FermionOnSphereWithSpinLzInvertTable[initialState & 0xff] << 56;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 8) & 0xff] << 48;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 16) & 0xff] << 40;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 24) & 0xff] << 32;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 32) & 0xff] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 40) & 0xff] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 48) & 0xff] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[initialState >> 56];  
#else
  unsigned long TmpState = FermionOnSphereWithSpinLzInvertTable[initialState & 0xff] << 24;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 8) & 0xff] << 16;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[(initialState >> 16) & 0xff] << 8;
  TmpState |= FermionOnSphereWithSpinLzInvertTable[initialState >> 24];
#endif
  initialState >>= this->InvertShift;
  TmpState >>= this->InvertUnshift;
  if (TmpState < initialState)
    {
      return (TmpState | FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT);
    }
  else
    if (TmpState != initialState)
      return (initialState | FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT);
    else
      return initialState;
}

// factorized code that is used to symmetrize result of the AdxAdy operations
//
// state = reference on the state that has been produced with the AdxAdy operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int FermionOnSphereWithSpinLzSymmetry::SymmetrizeAdAdResult(unsigned long& state, double& coefficient)
{
   unsigned long TmpState2 = state;
   state = this->GetSignedCanonicalState(state);
   if ((state & FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT) == 0x0ul)
     {
       this->GetStateSingletParity(TmpState2);
       if (((1.0 - 2.0 * ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul))) * this->LzParitySign) < 0.0)
	 return this->HilbertSpaceDimension;
      if ((this->ProdASignature & FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT) != 0x0ul)
	coefficient *= M_SQRT2;
      int NewLzMax = 1 + (this->LzMax << 1);
      while (((state & FERMION_SPHERE_SU2_SYMMETRIC_MASK) >> NewLzMax) == 0)
	--NewLzMax;
      return this->FindStateIndex(state, NewLzMax);
     }
  int NewLzMax = 1 + (this->LzMax << 1);
  while (((state & FERMION_SPHERE_SU2_SYMMETRIC_MASK) >> NewLzMax) == 0)
    --NewLzMax;
  if ((TmpState2 & FERMION_SPHERE_SU2_SYMMETRIC_MASK) != (state  & FERMION_SPHERE_SU2_SYMMETRIC_MASK))
    {
      this->GetStateSingletParity(TmpState2);
      coefficient *= (1.0 - 2.0 * ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
      coefficient *= this->LzParitySign;
    }
  if ((this->ProdASignature & FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT) == 0x0ul)
    coefficient *= M_SQRT1_2;
  return this->FindStateIndex(state, NewLzMax);
}

#endif


