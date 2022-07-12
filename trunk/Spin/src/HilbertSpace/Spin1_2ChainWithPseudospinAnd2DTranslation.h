////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                    class author: Cecile Repellin                           //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with a fixed Sz value                //
//               and a pseudospin 1/2 (not a conserved quantity)              //
//                        and 2D translations                                 //
//                                                                            //
//                        last modification : 25/11/2016                      //
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


#ifndef SPIN1_2CHAINWITHPSEUDOSPINAND2DTRANSLATION_H
#define SPIN1_2CHAINWITHPSEUDOSPINAND2DTRANSLATION_H


#include "config.h"
#include "HilbertSpace/Spin1_2ChainWithPseudospin.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;


class Spin1_2ChainWithPseudospinAnd2DTranslation : public Spin1_2ChainWithPseudospin
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
  
  // state description to be kept in cache to perform operations
  unsigned long TransientState;
  
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

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given maxMomentum sector
  int LookUpTableMemorySize;
  // shift used in each maxMomentum sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used maxMomentum value of the state an the second 
  int** LookUpTable;
  
  unsigned long LookUpTableMask;
  int LookUpPosition;
  int LookUpTableSize;

 public:


  // default constructor
  //
  Spin1_2ChainWithPseudospinAnd2DTranslation ();

  // constructor for complete Hilbert space
  //
  // chainLength = number of spin 1/2
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2ChainWithPseudospinAnd2DTranslation (int chainLength, int sz, int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, int memorySize);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainWithPseudospinAnd2DTranslation (const Spin1_2ChainWithPseudospinAnd2DTranslation& chain);

  // destructor
  //
  ~Spin1_2ChainWithPseudospinAnd2DTranslation ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainWithPseudospinAnd2DTranslation& operator = (const Spin1_2ChainWithPseudospinAnd2DTranslation& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  virtual int TotalSz (int index);
  
  // return eigenvalue of Sz_i Sz_j associated to a given state (acts only on spin part of many-body state)
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  virtual double SziSzj (int i, int j, int state);
  
  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // operator acting on pseudospin on site i (off-diagonal part)
  //
  // i = position of pseudospin operator
  // state = index of the state to be applied on JAi operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of the resulting state
  virtual int JOffDiagonali (int i, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // operator acting on pseudospin on site i and j(off-diagonal part)
  //
  // i = position of pseudospin operator
  // j = position of pseudospin operator
  // state = index of the state to be applied on JAi operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of the resulting state
  virtual int JoffiJoffj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // operator acting on pseudospin on site i (diagonal) and j(diagonal part)
  //
  // i = position of pseudospin operator
  // j = position of pseudospin operator
  // state = index of the state to be applied on JAi operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of the resulting state
  virtual int SmiSpjJiJj (int i, int j, int state, double* couplingI, double* couplingJ, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  
  // operator acting on pseudospin on site i (off-diagonal) and j(diagonal part)
  //
  // i = position of pseudospin operator
  // j = position of pseudospin operator
  // state = index of the state to be applied on JAi operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of the resulting state
  virtual int SmiSpjJoffiJj (int i, int j, int state, double* coupling, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
   // operator acting on pseudospin on site i (diagonal) and j(off-diagonal part)
  //
  // i = position of pseudospin operator
  // j = position of pseudospin operator
  // state = index of the state to be applied on JAi operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of the resulting state
  virtual int SmiSpjJiJoffj (int i, int j, int state, double* coupling, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // operator acting on pseudospin on site i (off-diagonal) and j(diagonal part)
  //
  // i = position of pseudospin operator
  // j = position of pseudospin operator
  // state = index of the state to be applied on JAi operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of the resulting state
  virtual int SmiSpjJoffiJoffj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // operator acting on pseudospin on site i (diagonal part)
  //
  // i = position of pseudospin operator
  // state = index of the state to be applied on JAi operator
  // coupling = array where the coupling coefficients are stored
  // return value = numerical coefficient
  virtual double JDiagonali (int i, int state, double* coupling);
  
  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSpjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SziSzjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
    
  // give a state description value to this->TransientState
  //
  // state = index of the state to convert
  virtual void InitializeTransientState (int state);
  
  // apply SziSzj to transient state and return corresponding coefficient
  //
  // i = position of S- operator
  // j = position of S+ operator
  // return value = corresponding eigenvalue
  virtual double SziSzj(int i, int j);
  
  // apply S+S- to transient state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // return value = corresponding eigenvalue
  virtual double SmiSpj(int i, int j);
  
  // apply off-diagonal part of pseudospin operator to the transient state
  //
  // i = position of the Joff operator
  // return value = corresponding eigenvalue
  virtual double JOffDiagonali(int i);
  
  // apply diagonal part of pseudospin operator to the transient state
  //
  // i = position of the JDiag operator
  // coupling = array of coupling coefficients
  // return value = corresponding eigenvalue
  virtual double JDiagonali(int i, double* coupling);
  
  // compute the parity (prod_i Sz_i) for a given state
  //
  // state = index of the state to be applied on Sz_i operator
  // return value = total Sz value
  virtual unsigned long GetParity (int state);
  
  // translate a state assuming the system have periodic boundary
  // conditions (increasing the site index)
  //
  // nbrTranslations = number of translations to apply
  // state = index of the state to translate 
  // return value = index of resulting state
  virtual int TranslateState (int nbrTranslations, int state);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);
  
  // convert a state defined on a lattice with a number of sites equals to a multiple of three
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual RealVector ProjectToEffectiveSubspaceThreeToOne(ComplexVector& state, AbstractSpinChain* space);
  
  // apply the mirror symmetry to a state
  //
  // stateIndex = index of the current state in the Hlibert space
  // coefficient = coefficient associated with the transformation
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  virtual inline int ApplyMirrorSymmetry(int stateIndex, double& coefficient,
                                             int& nbrTranslationX, int& nbrTranslationY);

  // factorized code that is used to symmetrize the result of any operator action (for external use)
  //
  // i = index of original state (to compute orbit size)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeResult(int i, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  

 protected:

  // find state index
  //
  // stateDescription = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription);

  // evaluate Hilbert space dimension
  //
  // sz = twice the Sz value
  // nbrSites = number of sites
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int sz, int nbrSites);

  // generate all states with constraints
  // 
  // return value = dimension of Hilbert space
  long GenerateStates();
  
  // generate all states
  // 
  // nbrSpinUp = number of spin up
  // currentPosition = current position to consider in the chain
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long RawGenerateStates(int nbrSpinUp, int currentPosition, long pos);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

//   // return the Bosonic Occupation of a given state in the basis
//   //
//   // index = index of the state in the basis
//   // finalState = reference on the array where the monomial representation has to be stored
//   virtual void GetBosonicOccupation (unsigned int index, int * finalState);
  
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
  
  // apply the mirror symmetry to a state description
  //
  // stateDescription = reference on the state description  
  virtual inline void ApplyMirrorSymmetry(unsigned long& stateDescription, int& TmpSign);
  
  // get a linearized position index from the 2d coordinates
  //
  // xPosition = position along the x direction
  // yPosition = position along the y direction
  // return value = linearized index
  virtual int GetLinearizedIndex(int xPosition, int yPosition);
  
  // fill up the table for the mirror symmetry
  //
  virtual void FillMirrorSymmetryTable();
  
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

inline int Spin1_2ChainWithPseudospinAnd2DTranslation::SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  state = this->FindCanonicalForm(state, nbrTranslationX, nbrTranslationY);
//   int TmpMaxMomentum = this->ChainLength;
//   while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
//     --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[nbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslationX = (this->MaxXMomentum - nbrTranslationX) % this->MaxXMomentum;
      nbrTranslationY = (this->MaxYMomentum - nbrTranslationY) % this->MaxYMomentum;
     }
//   else
//     coefficient = 0.0;
//   cout << "sym result " << nbrTranslationX << " " << nbrTranslationY << endl;
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

inline int Spin1_2ChainWithPseudospinAnd2DTranslation::SymmetrizeResult(int i, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long state = this->FindCanonicalForm(this->TransientState, nbrTranslationX, nbrTranslationY);
//   int TmpMaxMomentum = this->ChainLength;
//   while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
//     --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
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

inline unsigned long Spin1_2ChainWithPseudospinAnd2DTranslation::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY)
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

inline bool Spin1_2ChainWithPseudospinAnd2DTranslation::TestMomentumConstraint(unsigned long stateDescription)
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

inline int Spin1_2ChainWithPseudospinAnd2DTranslation::FindOrbitSize(unsigned long stateDescription)
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

inline void Spin1_2ChainWithPseudospinAnd2DTranslation::ApplySingleXTranslation(unsigned long& stateDescription)
{
  stateDescription = (stateDescription >> this->StateXShift) | ((stateDescription & this->XMomentumMask) << this->ComplementaryStateXShift);
}

// apply a single translation in the y direction for a state description
//
// stateDescription = reference on the state description

inline void Spin1_2ChainWithPseudospinAnd2DTranslation::ApplySingleYTranslation(unsigned long& stateDescription)
{
  stateDescription = (((stateDescription & this->ComplementaryYMomentumFullMask) >> this->StateYShift) | ((stateDescription & this->YMomentumFullMask) << this->ComplementaryStateYShift));
}

// apply the inversion symmetry to a state description
//
// stateDescription = reference on the state description


inline void Spin1_2ChainWithPseudospinAnd2DTranslation::ApplyMirrorSymmetry(unsigned long& stateDescription, int& TmpSign)
{
   unsigned long InitialState = stateDescription;
   stateDescription = 0x0ul;
   TmpSign = 1;
   for (int i =0; i < this->ChainLength; ++i)
   {
      stateDescription |= ((InitialState >> (2*i)) & 0x3ul) << (2 * this->MirrorTransformationTable[i]);
      if (((InitialState >> (2*i)) & 0x1ul) != 0x0ul)
	TmpSign *= -1;
   }
}



// apply the mirror symmetry to a state
//
// stateIndex = index of the current state in the Hlibert space
// coefficient = coefficient associated with the transformation
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
inline int Spin1_2ChainWithPseudospinAnd2DTranslation::ApplyMirrorSymmetry(int stateIndex, double& coefficient,
                                             int& nbrTranslationX, int& nbrTranslationY)
{
      coefficient = 1.0;
      unsigned long TmpState = this->StateDescription[stateIndex];
      int TmpNbrStateInOrbit = this->NbrStateInOrbit[stateIndex];
      int TmpSign;
      this->ApplyMirrorSymmetry(TmpState, TmpSign);
      int TmpIndex = this->SymmetrizeResult(TmpState, TmpNbrStateInOrbit, coefficient, nbrTranslationX, nbrTranslationY);
      coefficient *= TmpSign;
      return TmpIndex;
}



// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// return value = linearized index

inline int Spin1_2ChainWithPseudospinAnd2DTranslation::GetLinearizedIndex(int xPosition, int yPosition)
{
  if (xPosition < 0)
    xPosition += this->MaxXMomentum;
  if (xPosition >= this->MaxXMomentum)
    xPosition -= this->MaxXMomentum;
  if (yPosition < 0)
    yPosition += this->MaxYMomentum;
  if (yPosition >= this->MaxYMomentum)
    yPosition -= this->MaxYMomentum;
  return ((xPosition * this->MaxYMomentum) + yPosition);
}
#endif


