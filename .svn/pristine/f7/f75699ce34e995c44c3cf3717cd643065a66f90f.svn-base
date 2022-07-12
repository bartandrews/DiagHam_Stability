////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of doubled spin 0 +1/2 chain with translations              //
//                                                                            //
//                        last modification : 21/01/2016                      //
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


#ifndef SPIN0_1_2_CHAINWITHTRANSLATIONS_H
#define SPIN0_1_2_CHAINWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"


#include <iostream>


using std::ostream;
class DoubledSpin0_1_2_ChainWithTranslations;
class DoubledSpin0_1_2_ChainWithTranslationsStaggered;

class Spin0_1_2_ChainWithTranslations : public AbstractSpinChainWithTranslations
{
  friend class DoubledSpin0_1_2_ChainWithTranslations;
  friend class DoubledSpin0_1_2_ChainWithTranslationsStaggered;
  friend class DoubledSpin1_ChainWithTranslations;
   
 protected:

  // array containing falg indicating if a state beloging to an orbit with a given number of member is compatible with momentum constraint
  bool* CompatibilityWithMomentum;

  // number of sites in the x direction
  int MaxXMomentum;
  // shift to apply to move a group spin using an elementary translation
  int StateShift;
  // shift to apply to move the spin from one end to the other one
  int ComplementaryStateShift;
  // mask to apply to get a group of spins
  unsigned long StateMask;
  

  int DiffSz;
  bool FixedSpinProjectionFlag;
  
  // shift to apply to a state to obtain an index to the look-up table 
  int LookUpTableShift;
  // look-up table (LookUpTable[i] gives the index of the smallest state that greater than i <<  LookUpTableShift)
  long* LookUpTable;
  int ShiftLookUpTable;
  int ShiftLookUpTableNegativeSz;

  // look-up table (LookUpTable[i] gives the index of the smallest state that greater than i <<  LookUpTableShift)
  long* LookUpTablePositifSz;
  // look-up table (LookUpTable[i] gives the index of the smallest state that greater than i <<  LookUpTableShift)
  long* LookUpTableNegativeSz;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

  // array describing each n-nody state 
  unsigned long* ChainDescription;
  
  int ShiftNegativeDiffSz;
  int ShiftNegativeSz;

 public:

  // default constructor
  //
  Spin0_1_2_ChainWithTranslations ();

  // constructor for complete Hilbert space with no restriction on momentum
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // memorySize = memory size in bytes allowed for look-up table
  // memorySlice = maximum amount of memory that can be allocated to partially evalauted the states
  Spin0_1_2_ChainWithTranslations (int chainLength,  int diffSz, int memorySize, int memorySlice);
  Spin0_1_2_ChainWithTranslations (int chainLength,  int memorySize=100000, int memorySlice=100000);
  
  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin0_1_2_ChainWithTranslations (int chainLength, int momentum,  int translationStep, int sz, int memorySize, int memorySlice);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin0_1_2_ChainWithTranslations (const Spin0_1_2_ChainWithTranslations & chain);

  // destructor
  //
  ~Spin0_1_2_ChainWithTranslations ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin0_1_2_ChainWithTranslations & operator = (const Spin0_1_2_ChainWithTranslations & chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  virtual int FindStateIndex(unsigned long stateDescription);


  virtual int SmiSpj (int i, int j, int state, double& coefficient);
  virtual int GetLocalSpin(int site, int state);


 protected:

  // return value of twice spin projection on (Oz) for a given state
  //
  // stateDescription = state to which the spin projection has to be evaluated
  // return value = twice spin projection on (Oz)
  int GetTotalSz (unsigned long stateDescription);
  virtual inline int TotalSz(int); 
    
  // create precalculation tables
  //
  void CreatePrecalculationTable();

  // find the canonical form of a state
  //
  // state = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // return value = canonical form of the state
  inline void FindCanonicalForm ( unsigned long stateDescription, unsigned long & canonicalState, int& nbrTranslation);

  // find the canonical form of a state and find how many translations are needed to obtain the same state
  //
  // stateDescription = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
  // return value = canonical form of the state
  inline void FindCanonicalForm ( unsigned long stateDescription, unsigned long & canonicalState, int& nbrTranslation, int& nbrTranslationToIdentity);
  
  // find how many translations are needed to obtain the same state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of translation needed to obtain the same state
  int FindNumberTranslation(unsigned long stateDescription);


  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);
  
  // generate all states corresponding to the constraints
  // 
  // lengthBra = length of the chain to be decided for bra spins
  // lengthBra = length of the chain to be decided for ket spins
  // diffSz = difference of spin projection between bra and ket chain
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int length, int diffSz, long pos);
  long GenerateStates(int length, long pos);
  
  long ShiftedEvaluateHilbertSpaceDimension(int length, int diffSz);
  long ShiftedEvaluateHilbertSpaceDimension(int length);

  double TotalSzSz (int index);
  double SziSzj (int i, int j, int state);
  int SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation);
  int SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation);
  int SpiSpi (int i, int state, double& coefficient, int& nbrTranslation);
  int SmiSmi (int i, int state, double& coefficient, int& nbrTranslation);
  int SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation);
  int SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation);
  int SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation);
  int Spi (int i, int state, double& coefficient, int& nbrTranslation);
  int Smi (int i, int state, double& coefficient, int& nbrTranslation);

  virtual inline unsigned long FindCanonicalForm(unsigned long state, int& nbrTranslation);
  inline void ApplySingleXTranslation(unsigned long& stateDescription);
  
  virtual void GenerateLookUpTable();
};



inline int Spin0_1_2_ChainWithTranslations::TotalSz(int index)
{
  return this-> DiffSz;
}

inline unsigned long Spin0_1_2_ChainWithTranslations::FindCanonicalForm(unsigned long state, int& nbrTranslation)
{
  unsigned long CanonicalForm;
  this->FindCanonicalForm(state,CanonicalForm, nbrTranslation);
  return CanonicalForm;
}

// apply a single translation in the x direction for a state description
//
// stateDescription = reference on the state description

inline void Spin0_1_2_ChainWithTranslations::ApplySingleXTranslation(unsigned long& stateDescription)
{
  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->StateMask) << this->ComplementaryStateShift);
}



// find the canonical form of a state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// return value = canonical form of the state

inline void Spin0_1_2_ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription,unsigned long & canonicalState, int& nbrTranslation)
{
  nbrTranslation = 0;
  canonicalState = stateDescription;
  int index = 1;  
  while (index < this->ChainLength)
    {
      this->ApplySingleXTranslation(stateDescription);
      if (stateDescription < canonicalState)
	{
	  canonicalState = stateDescription;
	  nbrTranslation = index;
	}
      ++index;
    }
}

// find the canonical form of a state and find how many translations are needed to obtain the same state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
// return value = canonical form of the state

inline void Spin0_1_2_ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription, unsigned long & canonicalState, int& nbrTranslation, int& nbrTranslationToIdentity)
{
  nbrTranslation = 0;
  nbrTranslationToIdentity = 1;
  canonicalState = stateDescription;
  unsigned long ReferenceState = stateDescription;
  this->ApplySingleXTranslation(stateDescription);
  while ((ReferenceState != stateDescription)  && (nbrTranslationToIdentity < this->ChainLength))
    {
      if (stateDescription < canonicalState)
	{
	  canonicalState = stateDescription;
	  nbrTranslation = nbrTranslationToIdentity;
	}
      this->ApplySingleXTranslation(stateDescription);
      ++nbrTranslationToIdentity;
    }
}

// find how many translations are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

inline int Spin0_1_2_ChainWithTranslations::FindNumberTranslation(unsigned long stateDescription)
{
  unsigned long TmpState = stateDescription;
  this->ApplySingleXTranslation(TmpState);
  int index = 1;  
  while (TmpState != stateDescription)
    {     
      this->ApplySingleXTranslation(TmpState);
      ++index;
    }
  return index;
}



#endif

