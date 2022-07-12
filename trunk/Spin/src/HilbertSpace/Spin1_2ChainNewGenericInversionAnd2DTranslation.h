////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of spin 1/2 chain with Sz contraint                //
//                and 2d translations any generic inversion symmetry          //
//                                                                            //
//                        last modification : 14/06/2017                      //
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


#ifndef SPIN1_2CHAINNEWGENERICINVERSIONAND2DTRANSLATION_H
#define SPIN1_2CHAINNEWGENERICINVERSIONAND2DTRANSLATION_H


#include "config.h"
#include "HilbertSpace/Spin1_2ChainNewSzSymmetryAnd2DTranslation.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;
using std::hex;
using std::dec;


class Spin1_2ChainNewGenericInversionAnd2DTranslation : public Spin1_2ChainNewSzSymmetryAnd2DTranslation
{

 protected:
  
  // sign of the inversion sector
  double InversionSymmetrySector;
  // table used to perform the inversion symmetry (mapping each site index) 
  int* InversionSymmetryTable;

 public:

  // default constructor
  //
  Spin1_2ChainNewGenericInversionAnd2DTranslation ();

  // constructor for the Hilbert space
  //
  // nbrSite = total number or spins
  // sz = value of the total magnetization
  // inversionSymmetrySectorSector = inversion symmetry sector (can be either +1 or -1)
  // inversionList = array that give the mapped indices under the inversion symmetry
  // xMomentum = momentum along the x direction
  // maxXMomentum = number of sites in the x direction
  // yMomentum = momentum along the y direction
  // maxYMomentum = number of sites in the y direction
  // memory = amount of memory granted for precalculations
  Spin1_2ChainNewGenericInversionAnd2DTranslation (int nbrSite, int sz, int inversionSymmetrySector, int* inversionList,
						   int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, unsigned long memory = 10000000);
  
  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  Spin1_2ChainNewGenericInversionAnd2DTranslation (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainNewGenericInversionAnd2DTranslation (const Spin1_2ChainNewGenericInversionAnd2DTranslation& chain);

  // destructor
  //
  ~Spin1_2ChainNewGenericInversionAnd2DTranslation ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainNewGenericInversionAnd2DTranslation& operator = (const Spin1_2ChainNewGenericInversionAnd2DTranslation& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);
    
 protected:

  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // inversionSign = reference on the additional sign coming from the inversion symmetry
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY, double& inversionSign);

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
  // stateDescription = reference on the state description  
  virtual void ApplyInversionSymmetry(unsigned long& stateDescription);

  // compute the rescaling factors
  //
  virtual void ComputeRescalingFactors();

};

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// inversionSign = reference on the additional sign coming from the inversion symmetry
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long Spin1_2ChainNewGenericInversionAnd2DTranslation::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, 
												  int& nbrTranslationY, double& inversionSign)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;  
  unsigned long TmpStateDescription;  
  nbrTranslationX = 0;
  nbrTranslationY = 0;
  inversionSign = 1.0;
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
  this->ApplyInversionSymmetry(TmpStateDescription);
  if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      nbrTranslationX = 0;	      
      nbrTranslationY = 0;	      
      inversionSign = this->InversionSymmetrySector;
    }
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslationX = n;	      
	  nbrTranslationY = 0;	      
	  inversionSign = this->InversionSymmetrySector;
	}
    }
  this->ApplyInversionSymmetry(stateDescription);
  for (int m = 1; m < this->MaxYMomentum; ++m)
    {
      this->ApplySingleYTranslation(stateDescription);      
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslationX = 0;	      
	  nbrTranslationY = m;	      
	  inversionSign = this->InversionSymmetrySector;
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
	      inversionSign = this->InversionSymmetrySector;
	    }
	}
    }
  return CanonicalState;
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool Spin1_2ChainNewGenericInversionAnd2DTranslation::TestMomentumConstraint(unsigned long stateDescription)
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
  this->ApplyInversionSymmetry(TmpStateDescription2);
  if (stateDescription == TmpStateDescription2)
    {
      if (this->InversionSymmetrySector < 0.0)
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
      if (this->InversionSymmetrySector < 0.0)
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
  this->ApplyInversionSymmetry(TmpStateDescription2);
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

  if (this->InversionSymmetrySector < 0.0)
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

inline int Spin1_2ChainNewGenericInversionAnd2DTranslation::FindOrbitSize(unsigned long stateDescription)
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
  this->ApplyInversionSymmetry(TmpStateDescription2);
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
  this->ApplyInversionSymmetry(TmpStateDescription2);
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

// apply the inversion symmetry to a state description
//
// stateDescription = reference on the state description  

inline void Spin1_2ChainNewGenericInversionAnd2DTranslation::ApplyInversionSymmetry(unsigned long& stateDescription)
{
  unsigned long TmpState = 0x0ul;
  for (int i = 0; i < this->NbrSite; ++i)
    {
      TmpState |=  (stateDescription & 0x1ul) << this->InversionSymmetryTable[i];
      stateDescription >>= 1;
    }
  stateDescription = TmpState;
}


#endif


