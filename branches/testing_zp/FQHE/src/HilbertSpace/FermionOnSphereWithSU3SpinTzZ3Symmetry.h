////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//              class of fermions on sphere with SU(3) spin including         //
//                      the Z3 and Tz discrete symmetries                     //
//                                                                            //
//                        last modification : 20/02/2008                      //
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


#ifndef FERMIONONSPHEREWITHSU3SPINTZZ3SYMMETRY_H
#define FERMIONONSPHEREWITHSU3SPINTZZ3SYMMETRY_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinTzSymmetry.h"

#include <iostream>
#include <math.h>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


class FermionOnSphereWithSU3SpinTzZ3Symmetry :  public FermionOnSphereWithSU3SpinTzSymmetry
{


 public:

  // default constructor
  // 
  FermionOnSphereWithSU3SpinTzZ3Symmetry();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // minusTzParity = select the  Tz <-> -Tz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSU3SpinTzZ3Symmetry (int nbrFermions, int totalLz, int lzMax, bool minusTzParity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSU3SpinTzZ3Symmetry(const FermionOnSphereWithSU3SpinTzZ3Symmetry& fermions);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSU3SpinTzZ3Symmetry (char* fileName, unsigned long memory);

  // destructor
  //
  ~FermionOnSphereWithSU3SpinTzZ3Symmetry ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSU3SpinTzZ3Symmetry& operator = (const FermionOnSphereWithSU3SpinTzZ3Symmetry& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // convert a given state from symmetric basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector  
  virtual RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSU3Spin& nbodyBasis);

  // convert a given state from the usual n-body basis to the symmetric basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  virtual RealVector ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSU3Spin& nbodyBasis);

 protected:

  // get canonical expression of a given state
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state
  virtual unsigned long GetCanonicalState (unsigned long initialState);

  // get symmetry of a given state 
  //
  // initialState = referennce state whose symmetry has to be computed
  // return value = corresponding symmetry bit
  virtual unsigned long GetStateSymmetry (unsigned long& initialState);

  // compute the parity of the number of spin singlet aka the number of particle pairs (no more than two) per Lz value
  //
  // initialState = reference on the state whose parity has to be evaluated
  // return value = corresponding parity bit
  virtual unsigned long GetSignedCanonicalState (unsigned long& initialState);

  // compute the sign when applying Z3 rotation due to the fermionic statistics
  //
  // initialState = reference on the state whose sign has to be evaluated
  // rotationDirection = specify if the rotation has been done in the right or left direction
  // return value = corresponding sign 
  virtual double GetStateRotationSign(unsigned long& initialState, unsigned long rotationDirection);

  // compute the parity of the number of both doublet and triplet
  //
  // initialState = reference on the state whose parity has to be evaluated
  // return value = sign corresponding to the parity (+1 for even , -1 for odd)
  virtual double GetStateDoubletTripletParity(unsigned long& initialState);

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

inline unsigned long FermionOnSphereWithSU3SpinTzZ3Symmetry::GetCanonicalState (unsigned long initialState)
{
  unsigned long TmpState = (((initialState & FERMION_SPHERE_SU3_1_MASK) << 2) |
			    ((initialState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  if (TmpState == initialState)
    return initialState;
  unsigned long CanonicalState = initialState;
  TmpState = ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpState |= TmpState << 1;
  TmpState ^= initialState;
  if (TmpState < CanonicalState)
    CanonicalState = TmpState;
  TmpState = (((initialState & FERMION_SPHERE_SU3_1_MASK) << 2) |
	      ((initialState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  if (TmpState < CanonicalState)
    CanonicalState = TmpState;
  initialState = TmpState;
  TmpState = ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpState |= TmpState << 1;
  TmpState ^= initialState;
  if (TmpState < CanonicalState)
   CanonicalState  = TmpState;
  TmpState =  (((initialState & FERMION_SPHERE_SU3_1_MASK) << 2) |
	       ((initialState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  if (TmpState < CanonicalState)
    CanonicalState = TmpState;
  initialState = TmpState;
  TmpState = ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpState |= TmpState << 1;
  TmpState ^= initialState;
  if (TmpState < CanonicalState)
    return TmpState;
  else
    return CanonicalState;

}

// get symmetry of a given state 
//
// initialState = reference on the state whose symmetry has to be computed
// return value = corresponding symmetry bit

inline unsigned long FermionOnSphereWithSU3SpinTzZ3Symmetry::GetStateSymmetry (unsigned long& initialState)
{
  unsigned long TmpState =  (((initialState & FERMION_SPHERE_SU3_1_MASK) << 2) |
			     ((initialState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  if (TmpState == initialState)
    return (FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT | FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT);
  unsigned long TmpMask =  ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((initialState ^ TmpMask) == initialState)
    return FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT;
  TmpMask =  ((TmpState >> 1) ^ TmpState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((TmpState ^ TmpMask) == TmpState)
    return FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT;
  TmpState =  (((TmpState & FERMION_SPHERE_SU3_1_MASK) << 2) |
	       ((TmpState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  TmpMask =  ((TmpState >> 1) ^ TmpState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((TmpState ^ TmpMask) == TmpState)
    return FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT;
  return 0x0ul;
}

// get canonical expression of a given state and its symmetry
//
// initialState = reference on the state that has to be converted to its canonical expression
// return value = corresponding symmetry bit

inline unsigned long FermionOnSphereWithSU3SpinTzZ3Symmetry::GetSignedCanonicalState (unsigned long& initialState)
{ 
  unsigned long TmpSign = 0x0ul;
  unsigned long TmpRotationSign = 0x0ul;
  unsigned long CanonicalState = initialState;
  unsigned long TmpState = ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpState |= TmpState << 1;
  TmpState ^= initialState;
  if (TmpState < CanonicalState)
    {      
      TmpRotationSign = FERMION_SPHERE_SU3_TZ_FLIP_BIT;
      CanonicalState = TmpState;
    }
  else
    if (TmpState == initialState)
      TmpSign = FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT;
   TmpState = (((initialState & FERMION_SPHERE_SU3_1_MASK) << 2) |
	      ((initialState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  if (TmpState == initialState)
    return (FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT | FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT);
  if (TmpState <= CanonicalState)
    {
      TmpRotationSign = FERMION_SPHERE_SU3_Z3LEFTROTATION_BIT;
      CanonicalState = TmpState;
    }
  initialState = TmpState;
  TmpState = ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpState |= TmpState << 1;
  TmpState ^= initialState;
  if (TmpState < CanonicalState)
    {      
      TmpRotationSign = FERMION_SPHERE_SU3_TZ_FLIP_BIT |  FERMION_SPHERE_SU3_Z3LEFTROTATION_BIT;
      CanonicalState = TmpState;
    }
  else
    if (TmpState == initialState)
      TmpSign = FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT;
   TmpState = (((initialState & FERMION_SPHERE_SU3_1_MASK) << 2) |
	       ((initialState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  if (TmpState <= CanonicalState)
    {
      TmpRotationSign = FERMION_SPHERE_SU3_Z3RIGHTROTATION_BIT;
      CanonicalState = TmpState;
    }
  initialState = TmpState;
  TmpState = ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpState |= TmpState << 1;
  TmpState ^= initialState;
  if (TmpState < CanonicalState)
    {      
      TmpRotationSign = FERMION_SPHERE_SU3_TZ_FLIP_BIT |  FERMION_SPHERE_SU3_Z3RIGHTROTATION_BIT;
      CanonicalState = TmpState;
    }
  else
    if (TmpState == initialState)
      TmpSign = FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT;
  initialState = CanonicalState;
  return (TmpSign | TmpRotationSign);
}

// compute the sign when applying Z3 rotation due to the fermionic statistics
//
// initialState = reference on the state whose sign has to be evaluated
// rotationDirection = specify if the rotation has been done in the right or left direction
// return value = corresponding sign 

inline double FermionOnSphereWithSU3SpinTzZ3Symmetry::GetStateRotationSign(unsigned long& initialState, unsigned long rotationDirection)
{
  unsigned long TmpState = initialState;
  if ((rotationDirection &  FERMION_SPHERE_SU3_Z3LEFTROTATION_BIT) == 0x0ul)
    {
      initialState &= FERMION_SPHERE_SU3_3_MASK;
      TmpState &= initialState | (initialState  >> 1) | (initialState >> 2);
    }
  else
    {
      initialState &= FERMION_SPHERE_SU3_1_MASK;
      TmpState &= initialState | (initialState << 1) | (initialState << 2);   
    }
  initialState = TmpState;
  initialState ^= TmpState >> 1;
  initialState ^= TmpState >> 2;
  initialState ^= ((TmpState >> 2) | (TmpState >> 1) | TmpState);
  initialState &= FERMION_SPHERE_SU3_1_MASK;
#ifdef __64_BITS__
  initialState ^= (initialState >> 32);
#endif
  initialState ^= (initialState >> 16);
  initialState ^= (initialState >> 8);
  initialState ^= (initialState >> 4);
  initialState ^= (initialState >> 2);
  initialState ^= (initialState >> 1);
  return (1.0 - ((double) ((initialState & 1ul) << 1)));
}

// compute the parity of the number of both doublet and triplet
//
// initialState = reference on the state whose parity has to be evaluated
// return value = sign corresponding to the parity (+1 for even , -1 for odd)

inline double FermionOnSphereWithSU3SpinTzZ3Symmetry::GetStateDoubletTripletParity(unsigned long& initialState)
{
  unsigned long TmpState = initialState;
  TmpState &= (TmpState >> 1);
  TmpState |= (TmpState >> 1);
  unsigned long TmpState2 = initialState;
  TmpState2 &= (TmpState2 >> 2);
  TmpState |= TmpState2;
  TmpState &= FERMION_SPHERE_SU3_TZ_MASK;
#ifdef __64_BITS__
  TmpState ^= (TmpState >> 32);
#endif
  TmpState ^= (TmpState >> 16);
  TmpState ^= (TmpState >> 8);
  TmpState ^= (TmpState >> 4);
  TmpState ^= (TmpState >> 2);
  TmpState ^= (TmpState >> 1);
  return (1.0 - ((double) ((TmpState & 1ul) << 1)));
}

// factorized code that is used to symmetrize result of the AdxAdy operations
//
// state = reference on the state that has been produced with the AdxAdy operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int FermionOnSphereWithSU3SpinTzZ3Symmetry::SymmetrizeAdAdResult(unsigned long& state, double& coefficient)
{
  unsigned long TmpState2 = state;
  unsigned long TmpSign = this->GetSignedCanonicalState(state);
  if ((TmpSign & FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT) != 0x0ul)
    {
      if ((this->GetStateDoubletTripletParity(TmpState2) * this->TzParitySign) < 0.0)
	return this->HilbertSpaceDimension;
      if ((this->ProdASignature & FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT) == 0x0ul)
	coefficient *= MSQRT3;
      if ((this->ProdASignature & FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT) == 0x0ul)
	coefficient *= M_SQRT2;
      int NewLzMax = 2 + (this->LzMax * 3);
      while ((state >> NewLzMax) == 0x0ul)
	--NewLzMax;
      return this->FindStateIndex(state, NewLzMax);
    }
  int NewLzMax = 2 + (this->LzMax * 3);
  while ((state >> NewLzMax) == 0x0ul)
    --NewLzMax;
  if ((TmpSign  & FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT) != 0x0ul)
    {
      if ((this->GetStateDoubletTripletParity(TmpState2) * this->TzParitySign) < 0.0)
	{
	  return this->HilbertSpaceDimension;
	}
      if (TmpState2 != state)
	{
	  coefficient *= this->GetStateRotationSign(TmpState2, TmpSign);
	}
      if ((this->ProdASignature & FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT) != 0x0ul)
	coefficient *= MSQRT1_3;
      else
	if ((this->ProdASignature & FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT) == 0x0ul)
	  coefficient *= M_SQRT2;
      return this->FindStateIndex(state, NewLzMax);
   }
  if (TmpState2 != state)
    {
      if ((TmpSign &  FERMION_SPHERE_SU3_TZ_FLIP_BIT) != 0x0ul)
	{
	  coefficient *= this->GetStateSingletTzParity(state);
	  coefficient *= this->TzParitySign;
	}
      if ((TmpSign & FERMION_SPHERE_SU3_Z3ROTATION_BIT) != 0x0ul)
	{
	  coefficient *= this->GetStateRotationSign(TmpState2, TmpSign);
	}
    }
//  cout << endl;
  if ((this->ProdASignature & FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT) != 0x0ul)
    coefficient *= MSQRT1_3;
  if ((this->ProdASignature & FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT) != 0x0ul)
    coefficient *= M_SQRT1_2;
  return this->FindStateIndex(state, NewLzMax);
}
#endif


