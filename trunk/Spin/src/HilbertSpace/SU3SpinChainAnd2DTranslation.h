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
//                               and 2d translations                          //
//                                                                            //
//                        last modification : 06/02/2018                      //
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


#ifndef SU3SPINCHAINAND2DTRANSLATION_H
#define SU3SPINCHAINAND2DTRANSLATION_H


#include "config.h"
#include "HilbertSpace/SU3SpinChain.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;


class SU3SpinChainAnd2DTranslation : public SU3SpinChain
{

 protected:
  
  // number of sites in the x direction
  int MaxXMomentum;
  // number of sites in the y direction
  int MaxYMomentum;
  // momentum in the x direction
  int XMomentum;
  // momentum in the y direction
  int YMomentum;
  
  // bit shift that has to applied to perform a translation in the x direction 
  int StateXShift;
  // binary mask for the StateXShift first bits 
  unsigned long XMomentumMask;
  // bit shift to apply to move the first StateXShift bits at the end of a state description
  int ComplementaryStateXShift;

  // bit shift that has to applied to perform a translation in the y direction 
  int StateYShift;
 // binary mask for the StateYShift first bits 
  unsigned long YMomentumMask;
   // binary mask for the StateYShift first bits of each group
  unsigned long YMomentumFullMask;
  // binary mask for the ~YMomentumFullMask
  unsigned long ComplementaryYMomentumFullMask;
  // bit shift to apply to move the first StateYShift bits at the end of a state description
  int ComplementaryStateYShift;
  // number of bits that are related by a translation along the y direction 
  int YMomentumBlockSize;
  // binary mask corresponding to YMomentumBlockSize
  unsigned long YMomentumBlockMask;
  // number of independant blocks related by translations in the y direction 
  int NbrYMomentumBlocks;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

  unsigned long LookUpTableMask;
  int LookUpPosition;
  int LookUpTableSize;

 public:

  // default constructor
  //
  SU3SpinChainAnd2DTranslation ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // nbrSite = total number or spins
  // sz = value of the total magnetization
  // xMomentum = momentum along the x direction
  // maxXMomentum = number of sites in the x direction
  // yMomentum = momentum along the y direction
  // maxYMomentum = number of sites in the y direction
  // memory = amount of memory granted for precalculations
  SU3SpinChainAnd2DTranslation (int nbrSite, int tz, int y, int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  SU3SpinChainAnd2DTranslation (const SU3SpinChainAnd2DTranslation& chain);

  // destructor
  //
  ~SU3SpinChainAnd2DTranslation ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  SU3SpinChainAnd2DTranslation& operator = (const SU3SpinChainAnd2DTranslation& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // return index of resulting state from application of two-site exchange operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = index of resulting state
  virtual int Pij (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // return index of resulting state from application of four-site exchange operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = index of resulting state
  virtual int Pijkl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // convert a state defined in the real space basis into a state in the (Kx,Ky) basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertToKxKyBasis(ComplexVector& state, AbstractSpinChain* space);

  // convert a state defined in the (Kx,Ky) basis into a state in the real space basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space);
  
 protected:

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  long GenerateStates();

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = maximum bit set to one in stateDescription
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int maxMomentum);

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
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY);

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
  virtual void ApplySingleXTranslation(unsigned long& stateDescription);

  // apply a single translation in the y direction for a state description
  //
  // stateDescription = reference on the state description  
  virtual void ApplySingleYTranslation(unsigned long& stateDescription);
  
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

inline int SU3SpinChainAnd2DTranslation::SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, 
							      int& nbrTranslationX, int& nbrTranslationY)
{
  state = this->FindCanonicalForm(state, nbrTranslationX, nbrTranslationY);
  int TmpMaxMomentum = 2*this->ChainLength;
  while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
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

inline unsigned long SU3SpinChainAnd2DTranslation::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;  
  unsigned long TmpStateDescription;  
  nbrTranslationX = 0;
  nbrTranslationY = 0;
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
  return CanonicalState;
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool SU3SpinChainAnd2DTranslation::TestMomentumConstraint(unsigned long stateDescription)
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
  if (YSize == this->MaxYMomentum)
    {
      this->ApplySingleYTranslation(TmpStateDescription2); 
    }
  if ((((this->YMomentum * YSize * this->MaxXMomentum)
	+ (this->XMomentum * TmpXSize * this->MaxYMomentum)) % (this->MaxXMomentum * this->MaxYMomentum)) != 0)
    return false;
  return true;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int SU3SpinChainAnd2DTranslation::FindOrbitSize(unsigned long stateDescription)
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
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(stateDescription); 
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
  return (XSize * YSize);
}

// apply a single translation in the x direction for a state description
//
// stateDescription = reference on the state description

inline void SU3SpinChainAnd2DTranslation::ApplySingleXTranslation(unsigned long& stateDescription)
{
  stateDescription = (stateDescription >> this->StateXShift) | ((stateDescription & this->XMomentumMask) << this->ComplementaryStateXShift);
}

// apply a single translation in the y direction for a state description
//
// stateDescription = reference on the state description

inline void SU3SpinChainAnd2DTranslation::ApplySingleYTranslation(unsigned long& stateDescription)
{
  stateDescription = (((stateDescription & this->ComplementaryYMomentumFullMask) >> this->StateYShift) | ((stateDescription & this->YMomentumFullMask) << this->ComplementaryStateYShift));
}

#endif


