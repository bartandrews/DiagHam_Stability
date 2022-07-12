////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of Potts-3 chain with the translation symmetry            //
//                            and inversion symmetry                          //
//                                                                            //
//                        last modification : 10/12/2019                      //
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


#ifndef POTTS3CHAINWITHTRANSLATIONSANDINVERSION_H
#define POTTS3CHAINWITHTRANSLATIONSANDINVERSION_H


#include "config.h"
#include "HilbertSpace/Potts3ChainWithTranslations.h"

#include <iostream>


using std::ostream;
using std::hex;
using std::dec;


class HermitianMatrix;
class RealMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


//  precalculation table used to apply the inversion symmetry
//
static  unsigned long Potts3ChainWithTranslationsInversionTable[] = {
0x00ul, 0x40ul, 0x80ul, 0xc0ul, 0x10ul, 0x50ul, 0x90ul, 0xd0ul, 0x20ul, 0x60ul, 0xa0ul, 0xe0ul, 0x30ul, 0x70ul, 0xb0ul, 0xf0ul, 
0x04ul, 0x44ul, 0x84ul, 0xc4ul, 0x14ul, 0x54ul, 0x94ul, 0xd4ul, 0x24ul, 0x64ul, 0xa4ul, 0xe4ul, 0x34ul, 0x74ul, 0xb4ul, 0xf4ul, 
0x08ul, 0x48ul, 0x88ul, 0xc8ul, 0x18ul, 0x58ul, 0x98ul, 0xd8ul, 0x28ul, 0x68ul, 0xa8ul, 0xe8ul, 0x38ul, 0x78ul, 0xb8ul, 0xf8ul, 
0x0cul, 0x4cul, 0x8cul, 0xccul, 0x1cul, 0x5cul, 0x9cul, 0xdcul, 0x2cul, 0x6cul, 0xacul, 0xecul, 0x3cul, 0x7cul, 0xbcul, 0xfcul, 
0x01ul, 0x41ul, 0x81ul, 0xc1ul, 0x11ul, 0x51ul, 0x91ul, 0xd1ul, 0x21ul, 0x61ul, 0xa1ul, 0xe1ul, 0x31ul, 0x71ul, 0xb1ul, 0xf1ul, 
0x05ul, 0x45ul, 0x85ul, 0xc5ul, 0x15ul, 0x55ul, 0x95ul, 0xd5ul, 0x25ul, 0x65ul, 0xa5ul, 0xe5ul, 0x35ul, 0x75ul, 0xb5ul, 0xf5ul, 
0x09ul, 0x49ul, 0x89ul, 0xc9ul, 0x19ul, 0x59ul, 0x99ul, 0xd9ul, 0x29ul, 0x69ul, 0xa9ul, 0xe9ul, 0x39ul, 0x79ul, 0xb9ul, 0xf9ul, 
0x0dul, 0x4dul, 0x8dul, 0xcdul, 0x1dul, 0x5dul, 0x9dul, 0xddul, 0x2dul, 0x6dul, 0xadul, 0xedul, 0x3dul, 0x7dul, 0xbdul, 0xfdul, 
0x02ul, 0x42ul, 0x82ul, 0xc2ul, 0x12ul, 0x52ul, 0x92ul, 0xd2ul, 0x22ul, 0x62ul, 0xa2ul, 0xe2ul, 0x32ul, 0x72ul, 0xb2ul, 0xf2ul, 
0x06ul, 0x46ul, 0x86ul, 0xc6ul, 0x16ul, 0x56ul, 0x96ul, 0xd6ul, 0x26ul, 0x66ul, 0xa6ul, 0xe6ul, 0x36ul, 0x76ul, 0xb6ul, 0xf6ul, 
0x0aul, 0x4aul, 0x8aul, 0xcaul, 0x1aul, 0x5aul, 0x9aul, 0xdaul, 0x2aul, 0x6aul, 0xaaul, 0xeaul, 0x3aul, 0x7aul, 0xbaul, 0xfaul, 
0x0eul, 0x4eul, 0x8eul, 0xceul, 0x1eul, 0x5eul, 0x9eul, 0xdeul, 0x2eul, 0x6eul, 0xaeul, 0xeeul, 0x3eul, 0x7eul, 0xbeul, 0xfeul, 
0x03ul, 0x43ul, 0x83ul, 0xc3ul, 0x13ul, 0x53ul, 0x93ul, 0xd3ul, 0x23ul, 0x63ul, 0xa3ul, 0xe3ul, 0x33ul, 0x73ul, 0xb3ul, 0xf3ul, 
0x07ul, 0x47ul, 0x87ul, 0xc7ul, 0x17ul, 0x57ul, 0x97ul, 0xd7ul, 0x27ul, 0x67ul, 0xa7ul, 0xe7ul, 0x37ul, 0x77ul, 0xb7ul, 0xf7ul, 
0x0bul, 0x4bul, 0x8bul, 0xcbul, 0x1bul, 0x5bul, 0x9bul, 0xdbul, 0x2bul, 0x6bul, 0xabul, 0xebul, 0x3bul, 0x7bul, 0xbbul, 0xfbul, 
0x0ful, 0x4ful, 0x8ful, 0xcful, 0x1ful, 0x5ful, 0x9ful, 0xdful, 0x2ful, 0x6ful, 0xaful, 0xeful, 0x3ful, 0x7ful, 0xbful, 0xfful};


class Potts3ChainWithTranslationsAndInversion : public Potts3ChainWithTranslations
{

 protected:

  // sign of the inversion symmetry sector
  double InversionSector;

  // bit shift to apply before performing the inversion symmetry
  int InversionShift;
  // bit shift to apply after performing the inversion symmetry
  int InversionUnshift;

 public:

  // default constructor
  //
  Potts3ChainWithTranslationsAndInversion  ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // momemtum = momentum sector
  // inversionSector = inversion symmetry sector (can be either +1 or -1)
  // memorySize = memory size in bytes allowed for look-up table
  Potts3ChainWithTranslationsAndInversion (int chainLength, int momentum, int inversionSector, int memorySize);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // sz = twice the value of total Sz component
  // momemtum = momentum sector
  // inversionSector = inversion symmetry sector (can be either +1 or -1)
  // memorySize = memory size in bytes allowed for look-up table
  Potts3ChainWithTranslationsAndInversion (int chainLength, int sz, int momentum, int inversionSector, int memorySize);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Potts3ChainWithTranslationsAndInversion (const Potts3ChainWithTranslationsAndInversion& chain);

  // destructor
  //
  ~Potts3ChainWithTranslationsAndInversion ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Potts3ChainWithTranslationsAndInversion& operator = (const Potts3ChainWithTranslationsAndInversion& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture = 0);

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

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(unsigned long stateDescription);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);

  // apply the inversion symmetry to a state description
  //
  // stateDescription = reference on the state on which the inversion symmetry has to be applied
  virtual void ApplyInversionSymmetry(unsigned long& stateDescription);

  // generate all states corresponding to a given momnetum
  //
  virtual void GenerateStates();

};


// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// nbrStateInOrbit = original number of states in the orbit before the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int Potts3ChainWithTranslationsAndInversion::SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, 
								     int& nbrTranslations)
{
  double TmpSign;
  state = this->FindCanonicalForm(state, nbrTranslations, TmpSign);
  // int TmpMaxMomentum = 2 * this->ChainLength;
  // while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
  //   --TmpMaxMomentum;
  // int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  int TmpIndex = this->FindStateIndex(state);
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

inline unsigned long Potts3ChainWithTranslationsAndInversion::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations, double& szSymmetrySign)
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

inline int Potts3ChainWithTranslationsAndInversion::FindOrbitSize(unsigned long stateDescription)
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

// apply the inversion symmetry to a state description
//
// stateDescription = reference on the state on which the inversion symmetry has to be applied

inline void Potts3ChainWithTranslationsAndInversion::ApplyInversionSymmetry(unsigned long& stateDescription)
{
#ifdef __64_BITS__
  unsigned long InitialState = stateDescription & 0xfffffffffffffffcul;  
#else
  unsigned long InitialState = stateDescription & 0xfffffffcul;  
#endif	
  InitialState <<= this->InversionShift;
#ifdef __64_BITS__
  unsigned long TmpState = Potts3ChainWithTranslationsInversionTable[InitialState & 0xfful] << 56;
  TmpState |=Potts3ChainWithTranslationsInversionTable[(InitialState >> 8) & 0xfful] << 48;
  TmpState |=Potts3ChainWithTranslationsInversionTable[(InitialState >> 16) & 0xfful] << 40;
  TmpState |=Potts3ChainWithTranslationsInversionTable[(InitialState >> 24) & 0xfful] << 32;
  TmpState |=Potts3ChainWithTranslationsInversionTable[(InitialState >> 32) & 0xfful] << 24;
  TmpState |=Potts3ChainWithTranslationsInversionTable[(InitialState >> 40) & 0xfful] << 16;
  TmpState |=Potts3ChainWithTranslationsInversionTable[(InitialState >> 48) & 0xfful] << 8;
  TmpState |=Potts3ChainWithTranslationsInversionTable[InitialState >> 56]; 
#else
  unsigned long TmpState =Potts3ChainWithTranslationsInversionTable[InitialState & 0xfful] << 24;
  TmpState |=Potts3ChainWithTranslationsInversionTable[(InitialState >> 8) & 0xfful] << 16;
  TmpState |=Potts3ChainWithTranslationsInversionTable[(InitialState >> 16) & 0xfful] << 8;
  TmpState |=Potts3ChainWithTranslationsInversionTable[InitialState >> 24];
#endif	
  TmpState >>= this->InversionUnshift;
  TmpState |= stateDescription & 0x3ul;
  stateDescription = TmpState;
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool Potts3ChainWithTranslationsAndInversion::TestMomentumConstraint(unsigned long stateDescription)
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

#endif


