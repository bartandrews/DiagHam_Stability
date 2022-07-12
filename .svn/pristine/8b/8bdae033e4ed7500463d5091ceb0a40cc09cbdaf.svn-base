////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                             class of pair hopping p=2                      //
//            Hilbert space written as spin 1 chain with translations         //
//                    the combination of inversion and spin flip              //
//                             for more than 32 spins                         //
//                                                                            //
//                        last modification : 19/03/2019                      //
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


#ifndef PAIRHOPPINGP2SPIN1CHAINWITHTRANSLATIONSANDINVERSIONSZSYMMETRYLONG_H
#define PAIRHOPPINGP2SPIN1CHAINWITHTRANSLATIONSANDINVERSIONSZSYMMETRYLONG_H


#include "config.h"
#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslationsLong.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;
using std::hex;
using std::dec;



class PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong : public PairHoppingP2AsSpin1ChainWithTranslationsLong
{

 protected:

  // sign of the inversion symmetry sector
  double InversionSector;

  // bit shift to apply before performing the inversion symmetry
  int InversionShift;
  // bit shift to apply after performing the inversion symmetry
  int InversionUnshift;

  // mask needed for the Sz<->-Sz symmetry application
  ULONGLONG SzSymmetryMask;

 public:

  // default constructor
  //
  PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong ();

  // constructor for complete Hilbert space
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // inversionSector = inversion symmetry sector (can be either +1 or -1)
  // memory = amount of memory granted for precalculations
  PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong (int chainLength, int momentum, int inversionSector, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong (const PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong& chain);

  // destructor
  //
  ~PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong& operator = (const PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

 protected:

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // nbrStateInOrbit = original number of states in the orbit before the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeResult(ULONGLONG& state, int nbrStateInOrbit, double& coefficient, int& nbrTranslations);

  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
  // szSymmetrySign = reference on the additional sign coming from the inversion symmetry
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual ULONGLONG FindCanonicalForm(ULONGLONG stateDescription, int& nbrTranslations, double& szSymmetrySign);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(ULONGLONG stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(ULONGLONG stateDescription);

  // apply the inversion symmetry to a state description
  //
  // stateDescription = reference on the state on which the inversion symmetry has to be applied
  virtual void ApplyInversionSymmetry(ULONGLONG& stateDescription);

};


// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// nbrStateInOrbit = original number of states in the orbit before the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong::SymmetrizeResult(ULONGLONG& state, int nbrStateInOrbit, double& coefficient, 
											     int& nbrTranslations)
{
  double TmpSign;
  state = this->FindCanonicalForm(state, nbrTranslations, TmpSign);
  int TmpMaxMomentum = 2 * this->ChainLength - 1;
  while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= TmpSign;
      coefficient *= this->RescalingFactors[nbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslations = (this->MaxXMomentum - nbrTranslations) % this->MaxXMomentum;
     }
  return TmpIndex;
}

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// szSymmetrySign = reference on the additional sign coming from the inversion symmetry
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline ULONGLONG PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong::FindCanonicalForm(ULONGLONG stateDescription, int& nbrTranslations, double& szSymmetrySign)
{
  ULONGLONG CanonicalState = stateDescription;
  nbrTranslations = 0;
  szSymmetrySign = 1.0;
  ULONGLONG TmpStateDescription = stateDescription;  
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslations = n;	      
	}
    }
  TmpStateDescription = stateDescription;
  this->ApplyInversionSymmetry(TmpStateDescription);
  if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      nbrTranslations = 0;	      
      szSymmetrySign = this->InversionSector;      
    }
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);    
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslations = n;	      
	  szSymmetrySign = this->InversionSector;      
	}
    }
  return CanonicalState;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong::FindOrbitSize(ULONGLONG stateDescription)
{
  ULONGLONG TmpStateDescription = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  TmpStateDescription = stateDescription;
  this->ApplyInversionSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    return XSize;  
  int XSize2 = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  
  if (XSize2 != XSize)
    {
      return XSize;
    }
  return (2 * XSize);
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong::TestMomentumConstraint(ULONGLONG stateDescription)
{
  ULONGLONG TmpStateDescription = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);   
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  if (((this->Momentum * XSize) % this->MaxXMomentum) != 0)
    return false;
  TmpStateDescription = stateDescription;
  this->ApplyInversionSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    {
      if (this->InversionSector < 0.0)
	return false;
      else
	return true;
    }
  int XSize2 = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  

  if (XSize2 < XSize)
    {
      if (this->InversionSector < 0.0)
	{
	  if ((((this->Momentum * XSize2 * 2) + this->MaxXMomentum) % (2 * this->MaxXMomentum)) != 0)
	    return false;
	  else
	    return true;
	}
      if ((((this->Momentum * XSize2)) % this->MaxXMomentum) != 0)
	return false;
      else
	return true;  
    }
  return true;
}

// apply the inversion symmetry to a state description
//
// stateDescription = reference on the state on which the inversion symmetry has to be applied

inline void PairHoppingP2AsSpin1ChainWithTranslationsAndInversionSzSymmetryLong::ApplyInversionSymmetry(ULONGLONG& stateDescription)
{
  ULONGLONG InitialState = stateDescription;  
  InitialState <<= this->InversionShift;
#ifdef __128_BIT_LONGLONG__
  ULONGLONG TmpState = PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) InitialState) & 0xfful] << 120;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 8)) & 0xfful] << 112;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 16)) & 0xfful] << 104;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 24)) & 0xfful] << 96;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 32)) & 0xfful] << 88;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 40)) & 0xfful] << 80;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 48)) & 0xfful] << 72;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 56)) & 0xfful] << 64;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 64)) & 0xfful] << 56;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 72)) & 0xfful] << 48;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 80)) & 0xfful] << 40;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 88)) & 0xfful] << 32;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 96)) & 0xfful] << 24;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 104)) & 0xfful] << 16;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 112)) & 0xfful] << 8;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[(unsigned long) (InitialState >> 120)]; 
#else
  ULONGLONG TmpState = PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) InitialState) & 0xfful] << 56;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 8)) & 0xfful] << 48;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 16)) & 0xfful] << 40;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 24)) & 0xfful] << 32;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 32)) & 0xfful] << 24;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 40)) & 0xfful] << 16;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) (InitialState >> 48)) & 0xfful] << 8;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTableLong[((unsigned long) InitialState) >> 56];
#endif	
  TmpState >>= this->InversionUnshift;
  TmpState &= this->SzSymmetryMask;
  stateDescription = TmpState;
}

#endif


