////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of pair hopping p=1 (a.k.a PXP model)               //
//            Hilbert space written as spin 1 chain with translations         //
//                    the combination of inversion and spin flip              //
//                                                                            //
//                        last modification : 14/03/2019                      //
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


#ifndef PAIRHOPPINGP1SPIN1CHAINWITHTRANSLATIONSANDINVERSIONSZSYMMETRY_H
#define PAIRHOPPINGP1SPIN1CHAINWITHTRANSLATIONSANDINVERSIONSZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/PairHoppingP1AsSpin1ChainWithTranslations.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;
using std::hex;
using std::dec;


//  precalculation table used to apply the inversion+spin flip symmetry
//
static  unsigned long PairHoppingAsSpin1ChainWithTranslationsInversionTable[] = {
0xfful, 0x3ful, 0xbful, 0x3ful, 0xcful, 0x0ful, 0x8ful, 0x0ful, 0xeful, 0x2ful, 0xaful, 0x2ful, 0xcful, 0x0ful, 0x8ful, 0x0ful, 
0xf3ul, 0x33ul, 0xb3ul, 0x33ul, 0xc3ul, 0x03ul, 0x83ul, 0x03ul, 0xe3ul, 0x23ul, 0xa3ul, 0x23ul, 0xc3ul, 0x03ul, 0x83ul, 0x03ul, 
0xfbul, 0x3bul, 0xbbul, 0x3bul, 0xcbul, 0x0bul, 0x8bul, 0x0bul, 0xebul, 0x2bul, 0xabul, 0x2bul, 0xcbul, 0x0bul, 0x8bul, 0x0bul, 
0xf3ul, 0x33ul, 0xb3ul, 0x33ul, 0xc3ul, 0x03ul, 0x83ul, 0x03ul, 0xe3ul, 0x23ul, 0xa3ul, 0x23ul, 0xc3ul, 0x03ul, 0x83ul, 0x03ul, 
0xfcul, 0x3cul, 0xbcul, 0x3cul, 0xccul, 0x0cul, 0x8cul, 0x0cul, 0xecul, 0x2cul, 0xacul, 0x2cul, 0xccul, 0x0cul, 0x8cul, 0x0cul, 
0xf0ul, 0x30ul, 0xb0ul, 0x30ul, 0xc0ul, 0x00ul, 0x80ul, 0x00ul, 0xe0ul, 0x20ul, 0xa0ul, 0x20ul, 0xc0ul, 0x00ul, 0x80ul, 0x00ul, 
0xf8ul, 0x38ul, 0xb8ul, 0x38ul, 0xc8ul, 0x08ul, 0x88ul, 0x08ul, 0xe8ul, 0x28ul, 0xa8ul, 0x28ul, 0xc8ul, 0x08ul, 0x88ul, 0x08ul, 
0xf0ul, 0x30ul, 0xb0ul, 0x30ul, 0xc0ul, 0x00ul, 0x80ul, 0x00ul, 0xe0ul, 0x20ul, 0xa0ul, 0x20ul, 0xc0ul, 0x00ul, 0x80ul, 0x00ul, 
0xfeul, 0x3eul, 0xbeul, 0x3eul, 0xceul, 0x0eul, 0x8eul, 0x0eul, 0xeeul, 0x2eul, 0xaeul, 0x2eul, 0xceul, 0x0eul, 0x8eul, 0x0eul, 
0xf2ul, 0x32ul, 0xb2ul, 0x32ul, 0xc2ul, 0x02ul, 0x82ul, 0x02ul, 0xe2ul, 0x22ul, 0xa2ul, 0x22ul, 0xc2ul, 0x02ul, 0x82ul, 0x02ul, 
0xfaul, 0x3aul, 0xbaul, 0x3aul, 0xcaul, 0x0aul, 0x8aul, 0x0aul, 0xeaul, 0x2aul, 0xaaul, 0x2aul, 0xcaul, 0x0aul, 0x8aul, 0x0aul, 
0xf2ul, 0x32ul, 0xb2ul, 0x32ul, 0xc2ul, 0x02ul, 0x82ul, 0x02ul, 0xe2ul, 0x22ul, 0xa2ul, 0x22ul, 0xc2ul, 0x02ul, 0x82ul, 0x02ul, 
0xfcul, 0x3cul, 0xbcul, 0x3cul, 0xccul, 0x0cul, 0x8cul, 0x0cul, 0xecul, 0x2cul, 0xacul, 0x2cul, 0xccul, 0x0cul, 0x8cul, 0x0cul, 
0xf0ul, 0x30ul, 0xb0ul, 0x30ul, 0xc0ul, 0x00ul, 0x80ul, 0x00ul, 0xe0ul, 0x20ul, 0xa0ul, 0x20ul, 0xc0ul, 0x00ul, 0x80ul, 0x00ul, 
0xf8ul, 0x38ul, 0xb8ul, 0x38ul, 0xc8ul, 0x08ul, 0x88ul, 0x08ul, 0xe8ul, 0x28ul, 0xa8ul, 0x28ul, 0xc8ul, 0x08ul, 0x88ul, 0x08ul, 
0xf0ul, 0x30ul, 0xb0ul, 0x30ul, 0xc0ul, 0x00ul, 0x80ul, 0x00ul, 0xe0ul, 0x20ul, 0xa0ul, 0x20ul, 0xc0ul, 0x00ul, 0x80ul, 0x00ul};

//  precalculation table used to apply the Sz<->-Sz symmetry
//
static  unsigned long PairHoppingAsSpin1ChainWithTranslationsSzFlipTable[] = {0x3ul, 0x1ul, 0x2ul, 0x0ul};


class PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry : public PairHoppingP1AsSpin1ChainWithTranslations
{

 protected:

  // sign of the inversion symmetry sector
  double InversionSector;

  // bit shift to apply before performing the inversion symmetry
  int InversionShift;
  // bit shift to apply after performing the inversion symmetry
  int InversionUnshift;

  // mask needed for the Sz<->-Sz symmetry application
  unsigned long SzSymmetryMask;

 public:

  // default constructor
  //
  PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry ();

  // constructor for complete Hilbert space
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // inversionSector = inversion symmetry sector (can be either +1 or -1)
  // memory = amount of memory granted for precalculations
  PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry (int chainLength, int momentum, int inversionSector, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry (const PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry& chain);

  // destructor
  //
  ~PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry& operator = (const PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry& chain);

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
  virtual int SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, int& nbrTranslations);

 // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
  // szSymmetrySign = reference on the additional sign coming from the inversion symmetry
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations, double& szSymmetrySign);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(unsigned long stateDescription);

  // apply the inversion symmetry to a state description
  //
  // stateDescription = reference on the state on which the inversion symmetry has to be applied
  virtual void ApplyInversionSymmetry(unsigned long& stateDescription);

};


// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// nbrStateInOrbit = original number of states in the orbit before the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry::SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, 
											     int& nbrTranslations)
{
  double TmpSign;
  state = this->FindCanonicalForm(state, nbrTranslations, TmpSign);
  int TmpMaxMomentum = 2 * this->ChainLength;
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

inline unsigned long PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations, double& szSymmetrySign)
{
  unsigned long CanonicalState = stateDescription;
  nbrTranslations = 0;
  szSymmetrySign = 1.0;
  unsigned long TmpStateDescription = stateDescription;  
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

inline int PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry::FindOrbitSize(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
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

inline bool PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry::TestMomentumConstraint(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
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

inline void PairHoppingP1AsSpin1ChainWithTranslationsAndInversionSzSymmetry::ApplyInversionSymmetry(unsigned long& stateDescription)
{
#ifdef __64_BITS__
  unsigned long InitialState = stateDescription;  
#else
  unsigned long InitialState = stateDescription;  
#endif	
  InitialState <<= this->InversionShift;
#ifdef __64_BITS__
  unsigned long TmpState = PairHoppingAsSpin1ChainWithTranslationsInversionTable[InitialState & 0xfful] << 56;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTable[(InitialState >> 8) & 0xfful] << 48;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTable[(InitialState >> 16) & 0xfful] << 40;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTable[(InitialState >> 24) & 0xfful] << 32;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTable[(InitialState >> 32) & 0xfful] << 24;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTable[(InitialState >> 40) & 0xfful] << 16;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTable[(InitialState >> 48) & 0xfful] << 8;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTable[InitialState >> 56]; 
#else
  unsigned long TmpState =PairHoppingAsSpin1ChainWithTranslationsInversionTable[InitialState & 0xfful] << 24;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTable[(InitialState >> 8) & 0xfful] << 16;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTable[(InitialState >> 16) & 0xfful] << 8;
  TmpState |= PairHoppingAsSpin1ChainWithTranslationsInversionTable[InitialState >> 24];
#endif	
  TmpState >>= this->InversionUnshift;
  TmpState &= this->SzSymmetryMask;
  stateDescription = TmpState;
}

#endif


