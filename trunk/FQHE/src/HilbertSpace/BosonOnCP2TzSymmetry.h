////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Cecile Repellin                 //
//                     with Tz <-> -Tz symmetry                               //
//                                                                            //
//                          class of bosons on CP2                            //
//                     including Tz <-> -Tz symmetry                          //
//                                                                            //
//                        last modification : 18/01/2013                      //
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


#ifndef BOSONONCP2TZSYMMETRY_H
#define BOSONONCP2TZSYMMETRY_H

#include "config.h"
#include "HilbertSpace/BosonOnCP2.h"

#include <iostream>


class BosonOnCP2TzSymmetry : public BosonOnCP2
{

 protected:
// additional sign due to the parity sector for the Tz<->-Tz symmetry
  double TzParitySign;  
// symmetry of the temporary state
  int ProdASignature;
  // Tz<->-Tz conjugate (bosonic representation) of a state
  unsigned long* TzStateBosonic;
 public:

  // default constructor
  // 
  BosonOnCP2TzSymmetry ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta = number of flux quanta (p)
  // totalJz = total value of jz
  // totalKz = total value of kz
  // minusTzParity = select the  Tz <-> -Tz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  BosonOnCP2TzSymmetry (int nbrBosons, int nbrFluxQuanta, int totalTz, int totalY, bool minusTzParity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnCP2TzSymmetry(const BosonOnCP2TzSymmetry& bosons);

  // destructor
  //
  ~BosonOnCP2TzSymmetry ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnCP2TzSymmetry& operator = (const BosonOnCP2TzSymmetry& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);
  
  // apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int n1, int n2);

  // apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient);
  
  // convert a given state from symmetric basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector  
  virtual RealVector ConvertToNbodyBasis(RealVector& state, BosonOnCP2& nbodyBasis);


 protected:
   
  // get canonical expression of a given state and its symmetry
  //
  // initialState = ID of the state (in the unsymmetrized basis table) that has to be converted to its canonical expression
  // tzSymmetry = reference to an integer that describes the symmetry of the state: 1 if the state is symmetric, 0 if not 
  // return value = corresponding canonical state
  virtual unsigned long GetCanonicalState (int initialState, int& tzSymmetry);
  
  // get canonical expression of a given state and its symmetry
  //
  // tzSymmetry = reference to an integer that describes the symmetry of the state: 1 if the state is symmetric, 0 if not 
  // canonicalFlag = reference on an integer that says if the state was canonical (0) or not (1)
  // return value = corresponding canonical state
  virtual unsigned long GetCanonicalStateFromTemporaryBosonicPartition (int& tzSymmetry, int& canonicalFlag);
  
  // get canonical expression of a given state and its symmetry
  //
  // tzSymmetry = reference to an integer that describes the symmetry of the state: 1 if the state is symmetric, 0 if not
  // return value = corresponding canonical state
  virtual unsigned long GetCanonicalStateFromFermionicPartition (unsigned long initialStateFermionic, int& tzSymmetry);
  
  // factorized code that is used to symmetrize result of the AdxAdy operations
  //
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state
  virtual int SymmetrizeAdAdResult(double& coefficient);

};

inline unsigned long BosonOnCP2TzSymmetry::GetCanonicalState (int initialState, int& tzSymmetry)
{
  
  this->FermionToBoson(this->FermionBasis->StateDescription[initialState], this->FermionBasis->StateLzMax[initialState], this->TemporaryState, this->TemporaryStateLzMax);
  unsigned long finalState;
  for (int i = 0; i < this->NbrLzValue; ++i)
  {
   this->TzStateBosonic[i] = 0x0ul; 
  }
  int TzStateBosonicLzMax = 0;
  for (int index = 0; index <= this->TemporaryStateLzMax; ++index)
  {
   if (this->TemporaryState[index] > 0)
    {
      int indexFinalState = this->GetLinearizedIndex(-this->quantumNumberTz[index], this->quantumNumberY[index], 1);
//       cout << index << "  " << indexFinalState << "  ;  " << this->quantumNumberTz[index] << "," << this->quantumNumberY[index] << endl;
      this->TzStateBosonic[indexFinalState] = this->TemporaryState[index];
      if ( indexFinalState > TzStateBosonicLzMax)
	TzStateBosonicLzMax = indexFinalState;
    }
  }
  finalState = this->BosonToFermion(this->TzStateBosonic, TzStateBosonicLzMax);
  tzSymmetry = 0;
  if (finalState == this->FermionBasis->StateDescription[initialState])
    tzSymmetry = 1;
  if (finalState > this->FermionBasis->StateDescription[initialState])
    return finalState;
  else
    return this->FermionBasis->StateDescription[initialState];
}

inline unsigned long BosonOnCP2TzSymmetry::GetCanonicalStateFromTemporaryBosonicPartition (int& tzSymmetry, int& canonicalFlag)
{
  unsigned long TzStateFermionic;
  unsigned long initialStateFermionic = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
  for (int i = 0; i < this->NbrLzValue; ++i)
  {
   this->TzStateBosonic[i] = 0x0ul; 
  }
  int TzStateBosonicLzMax = 0;
  for (int index = 0; index <= this->TemporaryStateLzMax; ++index)
  {
   if (this->TemporaryState[index] > 0)
    {
      int indexFinalState = this->GetLinearizedIndex(-this->quantumNumberTz[index], this->quantumNumberY[index], 1);
//       cout << index << "  " << indexFinalState << "  ;  " << this->quantumNumberTz[index] << "," << this->quantumNumberY[index] << endl;
      this->TzStateBosonic[indexFinalState] = this->TemporaryState[index];
      if ( indexFinalState > TzStateBosonicLzMax)
	TzStateBosonicLzMax = indexFinalState;
    }
  }
  TzStateFermionic = this->BosonToFermion(this->TzStateBosonic, TzStateBosonicLzMax);
  tzSymmetry = 0;
  canonicalFlag = 0;
  if (TzStateFermionic == initialStateFermionic)
    tzSymmetry = 1;
  if (TzStateFermionic > initialStateFermionic)
  {
    for (int i = 0; i <= this->LzMax; ++i)
    {
      this->TemporaryState[i] = this->TzStateBosonic[i];
    }
    
    canonicalFlag = 1;
    return TzStateFermionic;
  }
  else
    return initialStateFermionic;
}

inline unsigned long BosonOnCP2TzSymmetry::GetCanonicalStateFromFermionicPartition (unsigned long initialStateFermionic, int& tzSymmetry)
{
  int initialStateFermionicLzMax = this->LzMax + this->NbrBosons - 1;
  while ((initialStateFermionic >> initialStateFermionicLzMax) == 0x0ul)
	--initialStateFermionicLzMax;
  this->FermionToBoson(initialStateFermionic, initialStateFermionicLzMax, this->TemporaryState, this->TemporaryStateLzMax);
  unsigned long TzStateFermionic;
  for (int i = 0; i < this->NbrLzValue; ++i)
  {
   this->TzStateBosonic[i] = 0x0ul; 
  }
  int TzStateBosonicLzMax = 0;
  for (int index = 0; index <= this->TemporaryStateLzMax; ++index)
  {
   if (this->TemporaryState[index] > 0)
    {
//       cout << this->TemporaryStateLzMax << " " << this->quantumNumberTz[index] << " " << this->quantumNumberY[index] << " " << endl;
      int indexFinalState = this->GetLinearizedIndex(-this->quantumNumberTz[index], this->quantumNumberY[index], 1);
//       cout << index << "  " << indexFinalState << "  ;  " << this->quantumNumberTz[index] << "," << this->quantumNumberY[index] << endl;
      this->TzStateBosonic[indexFinalState] = this->TemporaryState[index];
      if ( indexFinalState > TzStateBosonicLzMax)
	TzStateBosonicLzMax = indexFinalState;
    }
  }
  TzStateFermionic = this->BosonToFermion(this->TzStateBosonic, TzStateBosonicLzMax);
  tzSymmetry = 0;
  if (TzStateFermionic == initialStateFermionic)
    tzSymmetry = 1;
  if (TzStateFermionic > initialStateFermionic)
  {
    for (int i = 0; i <= this->LzMax; ++i)
    {
      this->TemporaryState[i] = this->TzStateBosonic[i];
    }
    return TzStateFermionic;
  }
  else
    return initialStateFermionic;
}

// factorized code that is used to symmetrize result of the AdxAdy operations
//
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int BosonOnCP2TzSymmetry::SymmetrizeAdAdResult(double& coefficient)
{
  unsigned long TmpState;
  int tzSymmetry;
  int canonicalFlag;
  TmpState = this->GetCanonicalStateFromTemporaryBosonicPartition(tzSymmetry, canonicalFlag);
  int TmpStateLzMax = this->NbrLzValue + this->NbrBosons - 1;
  while (((TmpState >> TmpStateLzMax) & 1) == 0)
    --TmpStateLzMax;
  if (tzSymmetry != 0)
    {
      if (this->TzParitySign < 0.0)
	return this->HilbertSpaceDimension;
      if (this->ProdASignature == 0)
	coefficient *= M_SQRT2;
//       cout << this->FindStateIndex(TmpState, TmpStateLzMax) << endl;
      return this->FindStateIndex(TmpState, TmpStateLzMax);
    }
  if (canonicalFlag == 1)
    coefficient *= this->TzParitySign;
  if (this->ProdASignature != 0)
    coefficient *= M_SQRT1_2;
//   cout << this->FindStateIndex(TmpState, TmpStateLzMax) << endl;
  return this->FindStateIndex(TmpState, TmpStateLzMax);
}
#endif


