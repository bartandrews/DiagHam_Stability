////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of fermions on sphere with spin with             //
//                                 Sz<->-Sz symmetry                          //
//                                                                            //
//                        last modification : 13/08/2007                      //
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


#ifndef FERMIONONSPHEREWITHSPINSZSYMMETRY_H
#define FERMIONONSPHEREWITHSPINSZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"

#include <iostream>

using std::cout;
using std::endl;


class FermionOnSphereWithSpinSzSymmetry :  public FermionOnSphereWithSpinLzSzSymmetry
{


 public:

  // default constructor 
  //
  FermionOnSphereWithSpinSzSymmetry ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // minusParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinSzSymmetry (int nbrFermions, int totalLz, int lzMax, bool minusParity, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinSzSymmetry (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinSzSymmetry(const FermionOnSphereWithSpinSzSymmetry& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinSzSymmetry ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinSzSymmetry& operator = (const FermionOnSphereWithSpinSzSymmetry& fermions);

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
  
  // apply a_n_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next Adu call
  //
  // index = index of the state on which the operator has to be applied
  // n = first index for annihilation operator (spin up)
  // return value =  multiplicative factor 
  virtual double Au (int index, int n);
  
  // apply a_n_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next Add call
  //
  // index = index of the state on which the operator has to be applied
  // n = first index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double Ad (int index, int n);
  
  // apply a^+_m_u operator to the state produced using Au (or Ad) method (without destroying it)
  //
  // m = index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Adu (int m, double& coefficient);
  
  // apply a^+_m_u operator to the state produced using Au (or Ad) method (without destroying it)
  //
  // m = index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Add (int m, double& coefficient);
  
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

inline unsigned long FermionOnSphereWithSpinSzSymmetry::GetCanonicalState (unsigned long initialState)
{
  unsigned long TmpMask =  ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU2_SZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((initialState ^ TmpMask) < initialState)
    return (initialState ^ TmpMask);
  else
    return initialState;
}

// get symmetry of a given state 
//
// initialState = reference on the state whose symmetry has to be computed

inline void FermionOnSphereWithSpinSzSymmetry::GetStateSymmetry (unsigned long& initialState)
{
  unsigned long TmpMask =  ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU2_SZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((initialState ^ TmpMask) != initialState)    
    initialState |= FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT;
}

// get canonical expression of a given state and its symmetry
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state (with symmetry bit)

inline unsigned long FermionOnSphereWithSpinSzSymmetry::GetSignedCanonicalState (unsigned long initialState)
{
  unsigned long TmpMask =  ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU2_SZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((initialState ^ TmpMask) < initialState)
    {
      return ((initialState ^ TmpMask) | FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT);
    }
  else
    if ((initialState ^ TmpMask) != initialState)
      return (initialState | FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT);
    else
      return initialState;
}

// factorized code that is used to symmetrize result of the AdxAdy operations
//
// state = reference on the state that has been produced with the AdxAdy operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int FermionOnSphereWithSpinSzSymmetry::SymmetrizeAdAdResult(unsigned long& state, double& coefficient)
{
  unsigned long TmpState2 = state;
  state = this->GetSignedCanonicalState(state);
  if ((state & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul)
    {
      this->GetStateSingletParity(TmpState2);
      if (((1.0 - 2.0 * ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul))) * this->SzParitySign) < 0.0)
	return this->HilbertSpaceDimension;
      int NewLzMax = 1 + (this->LzMax << 1);
      while (((state & FERMION_SPHERE_SU2_SYMMETRIC_MASK) >> NewLzMax) == 0)
	--NewLzMax;
      if ((this->ProdASignature & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) != 0x0ul)
	coefficient *= M_SQRT2;
      return this->FindStateIndex(state, NewLzMax);
    }
  int NewLzMax = 1 + (this->LzMax << 1);
  while (((state & FERMION_SPHERE_SU2_SYMMETRIC_MASK) >> NewLzMax) == 0)
    --NewLzMax;
  if ((TmpState2 & FERMION_SPHERE_SU2_SYMMETRIC_MASK) != (state  & FERMION_SPHERE_SU2_SYMMETRIC_MASK))
    {
      this->GetStateSingletParity(TmpState2);
      coefficient *= (1.0 - 2.0 * ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
      coefficient *= this->SzParitySign;
    }
  if ((this->ProdASignature & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul)
    coefficient *= M_SQRT1_2;
  return this->FindStateIndex(state, NewLzMax);
}

#endif


