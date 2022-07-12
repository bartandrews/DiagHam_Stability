////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//              class of fermions on sphere with SU(3) spin including         //
//                          Tz<->-Tz discrete symmetry                        //
//                                                                            //
//                        last modification : 04/02/2008                      //
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


#ifndef FERMIONONSPHEREWITHSU3SPINTZSYMMETRY_H
#define FERMIONONSPHEREWITHSU3SPINTZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"

#include <iostream>


class FermionOnSphere;


#define FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT    0x1ul
#define FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT     0x2ul
#define FERMION_SPHERE_SU3_TZZ3_SYMMETRIC_BIT     0x3ul
#define FERMION_SPHERE_SU3_TZSINGLETPARITY_BIT   0x10000ul
#define FERMION_SPHERE_SU3_TZSINGLETPARITY_SHIFT 16
#define FERMION_SPHERE_SU3_Z3ROTATIONSIGN_BIT   0x20000ul
#define FERMION_SPHERE_SU3_Z3ROTATIONSIGN_SHIFT 17
#define FERMION_SPHERE_SU3_Z3ROTATION_BIT       0xc0000ul
#define FERMION_SPHERE_SU3_Z3LEFTROTATION_BIT   0x40000ul
#define FERMION_SPHERE_SU3_Z3RIGHTROTATION_BIT  0x80000ul
#define FERMION_SPHERE_SU3_TZ_FLIP_BIT         0x100000ul
#ifdef __64_BITS__
#define FERMION_SPHERE_SU3_TZ_MASK 0x1249249249249249ul
#define FERMION_SPHERE_SU3_1_MASK  0x1249249249249249ul
#define FERMION_SPHERE_SU3_3_MASK  0x4924924924924924ul
#define FERMION_SPHERE_SU3_23_MASK 0x6db6db6db6db6db6ul
#define FERMION_SPHERE_SU3_TZ_UNFLIP_BIT 0xfffffffffefffffful
#else
#define FERMION_SPHERE_SU3_TZ_MASK 0x09249249ul
#define FERMION_SPHERE_SU3_1_MASK  0x09249249ul
#define FERMION_SPHERE_SU3_3_MASK  0x24924924ul
#define FERMION_SPHERE_SU3_23_MASK 0x36db6db6ul
#define FERMION_SPHERE_SU3_TZ_UNFLIP_BIT 0xfefffffful
#endif
#define MSQRT3   1.73205080756888
#define MSQRT1_3 0.577350269189626
#define MSQRT1_6 0.408248290463863


class FermionOnSphereWithSU3SpinTzSymmetry :  public FermionOnSphereWithSU3Spin
{

 protected:

  // additional sign due to the parity sector for the Tz<->-Tz symmetry
  double TzParitySign;
  // additional sign due to the parity sector for the Y<->-Y symmetry
  double YParitySign;
  // additional sign due to the parity sector for the Lz<->-Lz symmetry
  double LzParitySign;

  // symmetry of the temporary state
  int ProdASignature;

 public:

  // default constructor
  // 
  FermionOnSphereWithSU3SpinTzSymmetry();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalY = three time the total Y value
  // minusTzParity = select the  Tz <-> -Tz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSU3SpinTzSymmetry (int nbrFermions, int totalLz, int lzMax, int totalY, bool minusTzParity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSU3SpinTzSymmetry(const FermionOnSphereWithSU3SpinTzSymmetry& fermions);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSU3SpinTzSymmetry (char* fileName, unsigned long memory);

  // destructor
  //
  ~FermionOnSphereWithSU3SpinTzSymmetry ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSU3SpinTzSymmetry& operator = (const FermionOnSphereWithSU3SpinTzSymmetry& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  bool WriteHilbertSpace (char* fileName);

 // read Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description is stored
  // return value = true if no error occured
  virtual bool ReadHilbertSpace (char* fileName);

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

  // apply a^+_m_1 a_m_1 operator to a given state (only state 1 Tz=+1/2, Y=+1/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_1 a_m_1
  virtual double Ad1A1 (int index, int m);

  // apply a^+_m_2 a_m_2 operator to a given state (only state 2 Tz=-1/2, Y=+1/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_2 a_m_2
  virtual double Ad2A2 (int index, int m);

  // apply a^+_m_3 a_m_3 operator to a given state (only state 3 Tz=0, Y=-2/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_3 a_m_3
  virtual double Ad3A3 (int index, int m);

  // apply a_n1_1 a_n2_1 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A1A1 (int index, int n1, int n2);

  // apply a_n1_1 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A1A2 (int index, int n1, int n2);

  // apply a_n1_1 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A1A3 (int index, int n1, int n2);

  // apply a_n1_2 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A2A2 (int index, int n1, int n2);

  // apply a_n1_2 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A2A3 (int index, int n1, int n2);

  // apply a_n1_3 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A3A3 (int index, int n1, int n2);

  // apply a^+_m1_1 a^+_m2_1 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1Ad1 (int m1, int m2, double& coefficient);

  // apply a^+_m1_1 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1Ad2 (int m1, int m2, double& coefficient);

  // apply a^+_m1_1 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1Ad3 (int m1, int m2, double& coefficient);

  // apply a^+_m1_2 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2Ad2 (int m1, int m2, double& coefficient);

  // apply a^+_m1_2 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2Ad3 (int m1, int m2, double& coefficient);

  // apply a^+_m1_3 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad3Ad3 (int m1, int m2, double& coefficient);

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

  // get canonical expression of a given state and its symmetry
  //
  // initialState = reference on the state that has to be converted to its canonical expression
  // return value = corresponding symmetry bit
  virtual unsigned long GetSignedCanonicalState (unsigned long& initialState);

  // compute the parity of the number of spin singlet in the SU(2) sector
  //
  // initialState = reference on the state whose parity has to be evaluated
  // return value = sign corresponding to the parity (+1 for even , -1 for odd)
  virtual double GetStateSingletTzParity(unsigned long& initialState);

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

inline unsigned long FermionOnSphereWithSU3SpinTzSymmetry::GetCanonicalState (unsigned long initialState)
{
  unsigned long TmpMask =  ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((initialState ^ TmpMask) < initialState)
    return (initialState ^ TmpMask);
  else
    return initialState;
}

// get symmetry of a given state 
//
// initialState = reference on the state whose symmetry has to be computed
// return value = corresponding symmetry bit

inline unsigned long FermionOnSphereWithSU3SpinTzSymmetry::GetStateSymmetry (unsigned long& initialState)
{
  unsigned long TmpMask =  ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((initialState ^ TmpMask) != initialState)    
    return 0x0ul;
  else
    return FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT;
}

// get canonical expression of a given state and its symmetry
//
// initialState = reference on the state that has to be converted to its canonical expression
// return value = corresponding symmetry bit

inline unsigned long FermionOnSphereWithSU3SpinTzSymmetry::GetSignedCanonicalState (unsigned long& initialState)
{  
  unsigned long TmpMask =  ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((initialState ^ TmpMask) < initialState)
    {
      initialState ^= TmpMask;
      return 0x0ul;
    }
  else
    if ((initialState ^ TmpMask) == initialState)
      return FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT;
    else
      return 0x0ul;
}

// compute the parity of the number of spin singlet in the SU(2) sector
//
// initialState = reference on the state whose parity has to be evaluated
// return value = corresponding parity bit

inline double FermionOnSphereWithSU3SpinTzSymmetry::GetStateSingletTzParity(unsigned long& initialState)
{
  unsigned long TmpState = initialState;
  TmpState &= (TmpState >> 1);
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

inline int FermionOnSphereWithSU3SpinTzSymmetry::SymmetrizeAdAdResult(unsigned long& state, double& coefficient)
{
  unsigned long TmpState2 = state;
  if ((this->GetSignedCanonicalState(state) & FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT) != 0x0ul)
    {
      if ((this->GetStateSingletTzParity(TmpState2) * this->TzParitySign) < 0.0)
	return this->HilbertSpaceDimension;
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
  if (TmpState2 != state)
    {
      coefficient *= this->GetStateSingletTzParity(TmpState2);
      coefficient *= this->TzParitySign;
    }
  if ((this->ProdASignature & FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT) != 0x0ul)
    coefficient *= M_SQRT1_2;
  return this->FindStateIndex(state, NewLzMax);
}
#endif


