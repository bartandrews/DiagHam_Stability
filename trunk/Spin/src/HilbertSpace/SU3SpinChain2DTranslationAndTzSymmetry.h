////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                        Class author Cecile Repellin                        //
//                                                                            //
//                                                                            //
//                class of SU(3) spin chain with Sz contraint                 //
//                       2d translations and Tz symmetry                      //
//                                                                            //
//                        last modification : 07/02/2018                      //
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


#ifndef SU3SPINCHAIN2DTRANSLATIONANDTZSYMMETRY_H
#define SU3SPINCHAIN2DTRANSLATIONANDTZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/SU3SpinChainAnd2DTranslation.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;

//  precalculation table used to apply the Tz<->-Tz symmetry
//
static  unsigned long Spin1ChainWithSzInversionTable[] = {
0xfful, 0xfdul, 0xfeul, 0xfcul, 0xf7ul, 0xf5ul, 0xf6ul, 0xf4ul, 0xfbul, 0xf9ul, 0xfaul, 0xf8ul, 0xf3ul, 0xf1ul, 0xf2ul, 0xf0ul, 
0xdful, 0xddul, 0xdeul, 0xdcul, 0xd7ul, 0xd5ul, 0xd6ul, 0xd4ul, 0xdbul, 0xd9ul, 0xdaul, 0xd8ul, 0xd3ul, 0xd1ul, 0xd2ul, 0xd0ul, 
0xeful, 0xedul, 0xeeul, 0xecul, 0xe7ul, 0xe5ul, 0xe6ul, 0xe4ul, 0xebul, 0xe9ul, 0xeaul, 0xe8ul, 0xe3ul, 0xe1ul, 0xe2ul, 0xe0ul, 
0xcful, 0xcdul, 0xceul, 0xccul, 0xc7ul, 0xc5ul, 0xc6ul, 0xc4ul, 0xcbul, 0xc9ul, 0xcaul, 0xc8ul, 0xc3ul, 0xc1ul, 0xc2ul, 0xc0ul, 
0x7ful, 0x7dul, 0x7eul, 0x7cul, 0x77ul, 0x75ul, 0x76ul, 0x74ul, 0x7bul, 0x79ul, 0x7aul, 0x78ul, 0x73ul, 0x71ul, 0x72ul, 0x70ul, 
0x5ful, 0x5dul, 0x5eul, 0x5cul, 0x57ul, 0x55ul, 0x56ul, 0x54ul, 0x5bul, 0x59ul, 0x5aul, 0x58ul, 0x53ul, 0x51ul, 0x52ul, 0x50ul, 
0x6ful, 0x6dul, 0x6eul, 0x6cul, 0x67ul, 0x65ul, 0x66ul, 0x64ul, 0x6bul, 0x69ul, 0x6aul, 0x68ul, 0x63ul, 0x61ul, 0x62ul, 0x60ul, 
0x4ful, 0x4dul, 0x4eul, 0x4cul, 0x47ul, 0x45ul, 0x46ul, 0x44ul, 0x4bul, 0x49ul, 0x4aul, 0x48ul, 0x43ul, 0x41ul, 0x42ul, 0x40ul, 
0xbful, 0xbdul, 0xbeul, 0xbcul, 0xb7ul, 0xb5ul, 0xb6ul, 0xb4ul, 0xbbul, 0xb9ul, 0xbaul, 0xb8ul, 0xb3ul, 0xb1ul, 0xb2ul, 0xb0ul, 
0x9ful, 0x9dul, 0x9eul, 0x9cul, 0x97ul, 0x95ul, 0x96ul, 0x94ul, 0x9bul, 0x99ul, 0x9aul, 0x98ul, 0x93ul, 0x91ul, 0x92ul, 0x90ul, 
0xaful, 0xadul, 0xaeul, 0xacul, 0xa7ul, 0xa5ul, 0xa6ul, 0xa4ul, 0xabul, 0xa9ul, 0xaaul, 0xa8ul, 0xa3ul, 0xa1ul, 0xa2ul, 0xa0ul, 
0x8ful, 0x8dul, 0x8eul, 0x8cul, 0x87ul, 0x85ul, 0x86ul, 0x84ul, 0x8bul, 0x89ul, 0x8aul, 0x88ul, 0x83ul, 0x81ul, 0x82ul, 0x80ul, 
0x3ful, 0x3dul, 0x3eul, 0x3cul, 0x37ul, 0x35ul, 0x36ul, 0x34ul, 0x3bul, 0x39ul, 0x3aul, 0x38ul, 0x33ul, 0x31ul, 0x32ul, 0x30ul, 
0x1ful, 0x1dul, 0x1eul, 0x1cul, 0x17ul, 0x15ul, 0x16ul, 0x14ul, 0x1bul, 0x19ul, 0x1aul, 0x18ul, 0x13ul, 0x11ul, 0x12ul, 0x10ul, 
0x2ful, 0x2dul, 0x2eul, 0x2cul, 0x27ul, 0x25ul, 0x26ul, 0x24ul, 0x2bul, 0x29ul, 0x2aul, 0x28ul, 0x23ul, 0x21ul, 0x22ul, 0x20ul, 
0x0ful, 0x0dul, 0x0eul, 0x0cul, 0x07ul, 0x05ul, 0x06ul, 0x04ul, 0x0bul, 0x09ul, 0x0aul, 0x08ul, 0x03ul, 0x01ul, 0x02ul, 0x00ul};


class SU3SpinChain2DTranslationAndTzSymmetry : public SU3SpinChainAnd2DTranslation
{

 protected:
  
   // sign of the Tz<->-Tz symmetry sector
  double TzSymmetrySector;

  // mask needed for the Tz<->-Tz symmetry application
  unsigned long TzSymmetryMask;
  
 public:

  // default constructor
  //
  SU3SpinChain2DTranslationAndTzSymmetry ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // nbrSite = total number or spins
  // sz = value of the total magnetization
  // xMomentum = momentum along the x direction
  // maxXMomentum = number of sites in the x direction
  // yMomentum = momentum along the y direction
  // maxYMomentum = number of sites in the y direction
  // memory = amount of memory granted for precalculations
  SU3SpinChain2DTranslationAndTzSymmetry (int nbrSite, int tz, int y, int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, int tzSymmetrySector, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  SU3SpinChain2DTranslationAndTzSymmetry (const SU3SpinChain2DTranslationAndTzSymmetry& chain);

  // destructor
  //
  ~SU3SpinChain2DTranslationAndTzSymmetry ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  SU3SpinChain2DTranslationAndTzSymmetry& operator = (const SU3SpinChain2DTranslationAndTzSymmetry& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();
    
 protected:

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  long GenerateStates();

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // nbrStateInOrbit = original number of states in the orbit before the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY, double& tzSymmetrySign);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  inline int FindOrbitSize(unsigned long stateDescription);

  // apply a single translation in the x direction for a state description
  //
  // stateDescription = reference on the state description
  virtual void ApplyTzSymmetry(unsigned long& stateDescription);
    
  // compute the rescaling factors
  //
  virtual void ComputeRescalingFactors();
  
};

// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// nbrStateInOrbit = original number of states in the orbit before the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int SU3SpinChain2DTranslationAndTzSymmetry::SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, 
							      int& nbrTranslationX, int& nbrTranslationY)
{
  double TmpSign;
  state = this->FindCanonicalForm(state, nbrTranslationX, nbrTranslationY, TmpSign);
  int TmpMaxMomentum = 2*this->ChainLength;
  while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= TmpSign;
      coefficient *= this->RescalingFactors[nbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslationX = (this->MaxXMomentum - nbrTranslationX) % this->MaxXMomentum;
      nbrTranslationY = (this->MaxYMomentum - nbrTranslationY) % this->MaxYMomentum;
     }
  return TmpIndex;
}

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long SU3SpinChain2DTranslationAndTzSymmetry::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY, double& tzSymmetrySign)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;  
  unsigned long TmpStateDescription;  
  nbrTranslationX = 0;
  nbrTranslationY = 0;
  tzSymmetrySign = 1.0;
  TmpStateDescription = stateDescription;
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslationX = n;	      
	  nbrTranslationY = 0;	      
	}
    }
  for (int m = 1; m < this->MaxYMomentum; ++m)
    {
      this->ApplySingleYTranslation(stateDescription);      
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslationX = 0;	      
	  nbrTranslationY = m;	      
	}
      TmpStateDescription = stateDescription;
      for (int n = 1; n < this->MaxXMomentum; ++n)
	{
	  this->ApplySingleXTranslation(TmpStateDescription);      
	  if (TmpStateDescription < CanonicalState)
	    {
	      CanonicalState = TmpStateDescription;
	      nbrTranslationX = n;	      
	      nbrTranslationY = m;	      
	    }
	}
    }
  stateDescription = stateDescriptionReference;
  TmpStateDescription = stateDescription;
  this->ApplyTzSymmetry(TmpStateDescription);
  if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      nbrTranslationX = 0;	      
      nbrTranslationY = 0;	      
      tzSymmetrySign = this->TzSymmetrySector;
    }
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslationX = n;	      
	  nbrTranslationY = 0;	      
	  tzSymmetrySign = this->TzSymmetrySector;
	}
    }
  this->ApplyTzSymmetry(stateDescription);
  for (int m = 1; m < this->MaxYMomentum; ++m)
    {
      this->ApplySingleYTranslation(stateDescription);      
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslationX = 0;	      
	  nbrTranslationY = m;	      
	  tzSymmetrySign = this->TzSymmetrySector;
	}
      TmpStateDescription = stateDescription;
      for (int n = 1; n < this->MaxXMomentum; ++n)
	{
	  this->ApplySingleXTranslation(TmpStateDescription);      
	  if (TmpStateDescription < CanonicalState)
	    {
	      CanonicalState = TmpStateDescription;
	      nbrTranslationX = n;	      
	      nbrTranslationY = m;	      
	      tzSymmetrySign = this->TzSymmetrySector;
	    }
	}
    }
  return CanonicalState;
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool SU3SpinChain2DTranslationAndTzSymmetry::TestMomentumConstraint(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  unsigned long TmpStateDescription2 = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);   
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  if (((this->XMomentum * XSize) % this->MaxXMomentum) != 0)
    return false;
  int YSize = this->MaxYMomentum;
  int TmpXSize = 0;
  TmpStateDescription2 = stateDescription;
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(TmpStateDescription2); 
      TmpStateDescription = TmpStateDescription2;
      TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	}
      else
	{
	  TmpXSize = 0;
	}
    } 
  if ((((this->YMomentum * YSize * this->MaxXMomentum)
	+ (this->XMomentum * TmpXSize * this->MaxYMomentum)) % (this->MaxXMomentum * this->MaxYMomentum)) != 0)
    return false;

  TmpStateDescription2 = stateDescription;
  this->ApplyTzSymmetry(TmpStateDescription2);
  if (stateDescription == TmpStateDescription2)
    {
      if (this->TzSymmetrySector < 0.0)
	return false;
      else
	return true;
    }

  int XSize2 = 1;
  TmpStateDescription = TmpStateDescription2;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  
  if (XSize2 < XSize)
    {
      if (this->TzSymmetrySector < 0.0)
	{
	  if ((((this->XMomentum * XSize2 * 2 * this->MaxYMomentum) + (this->MaxXMomentum * this->MaxYMomentum)) % (2 * this->MaxXMomentum * this->MaxYMomentum)) != 0)
	    return false;
	  else
	    return true;
	}
      if ((((this->XMomentum * XSize2 * this->MaxYMomentum)) % (this->MaxXMomentum * this->MaxYMomentum)) != 0)
	return false;
      else
	return true;  
   }
  int YSize2 = YSize;
  TmpXSize = 0;
  TmpStateDescription2 = stateDescription;
  this->ApplyTzSymmetry(TmpStateDescription2);
  for (int m = 1; m < YSize2; ++m)
    {
      this->ApplySingleYTranslation(TmpStateDescription2); 
      TmpStateDescription = TmpStateDescription2;
      TmpXSize = 0;
      while ((TmpXSize < XSize2) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize2)
	{
	  YSize2 = m;
	}
      else
	{
	  TmpXSize = 0;
	}
    } 

  if (YSize == YSize2)
    return true;

  if (this->TzSymmetrySector < 0.0)
    {
      if ((((this->YMomentum * YSize2 * 2 * this->MaxXMomentum)
	    + (this->XMomentum * TmpXSize * 2 * this->MaxYMomentum) + (this->MaxXMomentum * this->MaxYMomentum)) % (2 * this->MaxXMomentum * this->MaxYMomentum)) != 0)
	return false;
      else
	return true;
    }
  if ((((this->YMomentum * YSize2 * this->MaxXMomentum)
	+ (this->XMomentum * TmpXSize * this->MaxYMomentum)) % (this->MaxXMomentum * this->MaxYMomentum)) != 0)
    return false;
  else
    return true;  
  return true;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int SU3SpinChain2DTranslationAndTzSymmetry::FindOrbitSize(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  unsigned long TmpStateDescription2 = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  int YSize = this->MaxYMomentum;
  TmpStateDescription2 = stateDescription;
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(TmpStateDescription2); 
      TmpStateDescription = TmpStateDescription2;
      int TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	}
    }

  TmpStateDescription2 = stateDescription;
  this->ApplyTzSymmetry(TmpStateDescription2);
  if (stateDescription == TmpStateDescription2)
    return (XSize * YSize);  

  int XSize2 = 1;
  TmpStateDescription = TmpStateDescription2;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  
  if (XSize2 != XSize)
    {
      return (XSize * YSize);
    }
  TmpStateDescription2 = stateDescription;
  this->ApplyTzSymmetry(TmpStateDescription2);
  int YSize2 = YSize;
  int TmpXSize;
  for (int m = 1; m < YSize2; ++m)
    {
      this->ApplySingleYTranslation(TmpStateDescription2); 
      TmpStateDescription = TmpStateDescription2;
      int TmpXSize = 0;
      while ((TmpXSize < XSize2) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize2)
	{
	  return (XSize * YSize);
	}
    }
  return (2 * XSize * YSize);
}


// apply the Sz<->-Sz symmetry to a state description
//
// stateDescription = reference on the state on which the Sz<->-Sz symmetry has to be applied

inline void SU3SpinChain2DTranslationAndTzSymmetry::ApplyTzSymmetry(unsigned long& stateDescription)
{
#ifdef __64_BITS__
  unsigned long TmpState = Spin1ChainWithSzInversionTable[stateDescription & 0xfful];
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 8) & 0xfful] << 8;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 16) & 0xfful] << 16;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 24) & 0xfful] << 24;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 32) & 0xfful] << 32;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 40) & 0xfful] << 40;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 48) & 0xfful] << 48;
  TmpState |= Spin1ChainWithSzInversionTable[stateDescription >> 56] << 56; 
#else
  unsigned long TmpState = Spin1ChainWithSzInversionTable[stateDescription & 0xfful];
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 8) & 0xfful] << 8;
  TmpState |= Spin1ChainWithSzInversionTable[(stateDescription >> 16) & 0xfful] << 16;
  TmpState |= Spin1ChainWithSzInversionTable[stateDescription >> 24] << 24;
#endif	
  TmpState &= this->TzSymmetryMask;
  stateDescription = TmpState;
}
#endif