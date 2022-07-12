////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                    class author: Cecile Repellin                           //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with a Sz = 0                        //
//        a pseudospin 1/2 (not a conserved quantity), spin inversion         //
//                     symmetry and 2D translations                           //
//                                                                            //
//                        last modification : 11/12/2016                      //
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


#ifndef SPIN1_2CHAINWITHPSEUDOSPINSZSYMMETRYAND2DTRANSLATION_H
#define SPIN1_2CHAINWITHPSEUDOSPINSZSYMMETRYAND2DTRANSLATION_H


#include "config.h"
#include "HilbertSpace/Spin1_2ChainWithPseudospinAnd2DTranslation.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;


class Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation : public Spin1_2ChainWithPseudospinAnd2DTranslation
{
 protected:
   
    // sign of the inversion sector
  double SzSymmetrySector;
  //mask to perform the spin flip operations
  unsigned long SzSymmetryMask;
  //mask to perform the spin flip operations
  unsigned long SzSymmetryComplementaryMask;

 public:


  // default constructor
  //
  Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation ();

  // constructor for complete Hilbert space
  //
  // chainLength = number of spin 1/2
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation (int nbrSite, int sz, int szSymmetrySector, int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, unsigned long memory) ;

  
  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation (char* fileName, unsigned long memory);
  
  
  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation (const Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation& chain);

  // destructor
  //
  ~Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation& operator = (const Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  bool WriteHilbertSpace (char* fileName);
  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
//   virtual ostream& PrintState (ostream& Str, int state);
  
  // convert a state defined in the (Kx,Ky) basis into a state in the real space basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  ComplexVector ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space);
  
  // convert a state from a SU(2) basis to another one, transforming the one body basis in each momentum sector
  //
  // initialState = state to transform  
  // targetState = vector where the transformed state has to be stored
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // firstComponent = index of the first component to compute in initialState
  // nbrComponents = number of consecutive components to compute
  virtual void TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, RealMatrix oneBodyBasis, AbstractSpinChain* space, long firstComponent = 0l, long nbrComponents = 0l);

  // recursive part of the convertion from a SU(2) basis to another one, transforming the one body basis in each momentum sector
  //
  // targetState = vector where the transformed state has to be stored
  // coefficient = current coefficient to assign
  // position = current particle consider in the n-body state
  // momentumIndices = array that gives the momentum partition of the initial n-body state
  // initialSU2Indices = array that gives the spin dressing the initial n-body state
  // currentSU2Indices = array that gives the spin dressing the current transformed n-body state
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  virtual void TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
				      int position, int* spinIndices, int* initialPseudospinIndices, int* currentPseudospinIndices, RealMatrix oneBodyBasis, int nbrStateInOrbit, Complex** exponentialFactors);
  
  // factorized code that is used to symmetrize the result of any operator action (for external use)
  //
  // i = index of original state (to compute orbit size)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeResult(int i, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
 protected:
   
  // generate all states with constraints
  // 
  // return value = dimension of Hilbert space
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
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY, double& inversionSign);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  inline int FindOrbitSize(unsigned long stateDescription);

  // apply the inversion symmetry to a state description
  //
  // stateDescription = reference on the state description  
  virtual void ApplySzSymmetry(unsigned long& stateDescription);
  
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

inline int Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  double TmpSign;
  state = this->FindCanonicalForm(state, nbrTranslationX, nbrTranslationY, TmpSign);
//   int TmpMaxMomentum = this->ChainLength;
//   while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
//     --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= TmpSign;
      coefficient *= this->RescalingFactors[nbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslationX = (this->MaxXMomentum - nbrTranslationX) % this->MaxXMomentum;
      nbrTranslationY = (this->MaxYMomentum - nbrTranslationY) % this->MaxYMomentum;
     }
//   else
//     coefficient = 0.0;
  return TmpIndex;
}

// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// i = index of initial state before operations
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::SymmetrizeResult(int i, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  double TmpSign;
  unsigned long state = this->FindCanonicalForm(this->TransientState, nbrTranslationX, nbrTranslationY, TmpSign);
  int TmpIndex = this->FindStateIndex(state);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= TmpSign;
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[i]][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslationX = (this->MaxXMomentum - nbrTranslationX) % this->MaxXMomentum;
      nbrTranslationY = (this->MaxYMomentum - nbrTranslationY) % this->MaxYMomentum;
     }
//   else
//     coefficient = 0.0;
//   cout << "sym result " << nbrTranslationX << " " << nbrTranslationY << endl;
  return TmpIndex;
}



// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY, double& inversionSign)
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
  this->ApplySzSymmetry(TmpStateDescription);
  if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      nbrTranslationX = 0;	      
      nbrTranslationY = 0;	      
      inversionSign = this->SzSymmetrySector;
    }
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslationX = n;	      
	  nbrTranslationY = 0;	      
	  inversionSign = this->SzSymmetrySector;
	}
    }
  this->ApplySzSymmetry(stateDescription);
  for (int m = 1; m < this->MaxYMomentum; ++m)
    {
      this->ApplySingleYTranslation(stateDescription);      
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslationX = 0;	      
	  nbrTranslationY = m;	      
	  inversionSign = this->SzSymmetrySector;
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
	      inversionSign = this->SzSymmetrySector;
	    }
	}
    }  
  return CanonicalState;
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::TestMomentumConstraint(unsigned long stateDescription)
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
  this->ApplySzSymmetry(TmpStateDescription2);
  if (stateDescription == TmpStateDescription2)
    {
      if (this->SzSymmetrySector < 0.0)
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
      if (this->SzSymmetrySector < 0.0)
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
  this->ApplySzSymmetry(TmpStateDescription2);
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

  if (this->SzSymmetrySector < 0.0)
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

inline int Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::FindOrbitSize(unsigned long stateDescription)
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
  this->ApplySzSymmetry(TmpStateDescription2);
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
  this->ApplySzSymmetry(TmpStateDescription2);
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

inline void Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation::ApplySzSymmetry(unsigned long& stateDescription)
{
  stateDescription = ((~stateDescription) & (this->SzSymmetryMask)) | ((stateDescription) & (this->SzSymmetryComplementaryMask));
}


#endif


