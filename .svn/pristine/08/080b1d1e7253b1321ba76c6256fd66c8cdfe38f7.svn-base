////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of spin 2 chain with translations                  //
//                  and the inversion and Sz<->-Sz symmetries                 //
//                                                                            //
//                        last modification : 12/11/2016                      //
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


#ifndef SPIN2CHAINWITHTRANSLATIONSANDSZINVERSIONSYMMETRIES_H
#define SPIN2CHAINWITHTRANSLATIONSANDSZINVERSIONSYMMETRIES_H


#include "config.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndInversionSymmetry.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;
using std::hex;
using std::dec;


class Spin2ChainWithTranslationsAndSzInversionSymmetries : public Spin2ChainWithTranslationsAndInversionSymmetry
{

 protected:

 public:

  // default constructor
  //
  Spin2ChainWithTranslationsAndSzInversionSymmetries ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // inversionSector = inversion symmetry sector (can be either +1 or -1)
  // szSymmetrySector = Sz<->-Sz symmetry sector (can be either +1 or -1)
  // memory = amount of memory granted for precalculations
  Spin2ChainWithTranslationsAndSzInversionSymmetries (int chainLength, int momentum, int inversionSector, int szSymmetrySector, unsigned long memory = 10000000);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momentum = total momentum of each state
  // sz = twice the value of total Sz component (should be equal to zero)
  // inversionSector = inversion symmetry sector (can be either +1 or -1)
  // szSymmetrySector = Sz<->-Sz symmetry sector (can be either +1 or -1)
  // memory = amount of memory granted for precalculations
  Spin2ChainWithTranslationsAndSzInversionSymmetries (int chainLength, int momentum, int inversionSector, int szSymmetrySector, int sz, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin2ChainWithTranslationsAndSzInversionSymmetries (const Spin2ChainWithTranslationsAndSzInversionSymmetries& chain);

  // destructor
  //
  ~Spin2ChainWithTranslationsAndSzInversionSymmetries ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin2ChainWithTranslationsAndSzInversionSymmetries& operator = (const Spin2ChainWithTranslationsAndSzInversionSymmetries& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

 protected:

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

  // compute the rescaling factors
  //
  virtual void ComputeRescalingFactors();

};

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// szSymmetrySign = reference on the additional sign coming from the inversion symmetry
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long Spin2ChainWithTranslationsAndSzInversionSymmetries::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations, double& szSymmetrySign)
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
  this->ApplySzSymmetry(TmpStateDescription);
  if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      nbrTranslations = 0;	      
      szSymmetrySign = this->SzSymmetrySector;      
    }
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);    
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslations = n;	      
	  szSymmetrySign = this->SzSymmetrySector;      
	}
    }
  this->ApplyInversionSymmetry(stateDescription);
  TmpStateDescription = stateDescription;
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
  TmpStateDescription = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription);
   if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      nbrTranslations = 0;	      
      szSymmetrySign = this->SzSymmetrySector * this->InversionSector;      
    }
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);    
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslations = n;	      
	  szSymmetrySign = this->SzSymmetrySector * this->InversionSector;      
	}
    }
 return CanonicalState;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int Spin2ChainWithTranslationsAndSzInversionSymmetries::FindOrbitSize(unsigned long stateDescription)
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
  this->ApplySzSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    {
      this->ApplyInversionSymmetry(TmpStateDescription);
      if (stateDescription == TmpStateDescription)
	{
	  return XSize;  
	}
      else
	{
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
    }
  int XSize2 = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  
  if (XSize2 != XSize)
    {
      this->ApplyInversionSymmetry(TmpStateDescription);
      if (stateDescription == TmpStateDescription)
	{
	  return XSize;  
	}
      else
	{
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
    }
  TmpStateDescription = stateDescription;
  this->ApplyInversionSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    {
      return (2 * XSize);  
    }
  XSize2 = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  
  if (XSize2 != XSize)
    {
      return (2 * XSize);
    }
  TmpStateDescription = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription);
  this->ApplyInversionSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    {
      return (2 * XSize);  
    }
  XSize2 = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  
  if (XSize2 != XSize)
    {
      return (2 * XSize);
    }
  return (4 * XSize);
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool Spin2ChainWithTranslationsAndSzInversionSymmetries::TestMomentumConstraint(unsigned long stateDescription)
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
  this->ApplySzSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    {
      if (this->SzSymmetrySector < 0.0)
	return false;
      else
	{
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
		    {
		      return true;
		    }
		}
	      if ((((this->Momentum * XSize2)) % this->MaxXMomentum) != 0)
		return false;
	      else
		{

		  return true;
		}
	    }
	  return true;
	}
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
      if (((this->SzSymmetrySector < 0.0) && ((((this->Momentum * XSize2 * 2) + this->MaxXMomentum) % (2 * this->MaxXMomentum)) != 0)) ||
	  ((this->SzSymmetrySector > 0.0) && ((((this->Momentum * XSize2)) % this->MaxXMomentum) != 0)))
	{
	  return false;
	}
    }

  // consider the sector where the inversion symmetry is applied
  TmpStateDescription = stateDescription;
  this->ApplyInversionSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    {
      if (this->InversionSector < 0.0)
	return false;
      else
	return true;
    }
  XSize2 = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  	  
  if (XSize2 < XSize)
    {
      if (((this->InversionSector < 0.0) && ((((this->Momentum * XSize2 * 2) + this->MaxXMomentum) % (2 * this->MaxXMomentum)) != 0)) ||
	  ((this->InversionSector > 0.0) && ((((this->Momentum * XSize2)) % this->MaxXMomentum) != 0)))
	{
	  return false;
	}
   }

  TmpStateDescription = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription);
  this->ApplyInversionSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    {
      if ((this->InversionSector * this->SzSymmetrySector) < 0.0)
	return false;
      else
	return true;
    }
  XSize2 = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  	  
  if (XSize2 < XSize)
    {
      if ((this->InversionSector * this->SzSymmetrySector) < 0.0)
	{
	  if ((((this->Momentum * XSize2 * 2) + this->MaxXMomentum) % (2 * this->MaxXMomentum)) != 0)
	    return false;
	  else
	    return true;
	}
      else
	{
	  if ((((this->Momentum * XSize2)) % this->MaxXMomentum) != 0)
	    return false;
	  else
	    return true;
	}
    }

  return true;
}


#endif


