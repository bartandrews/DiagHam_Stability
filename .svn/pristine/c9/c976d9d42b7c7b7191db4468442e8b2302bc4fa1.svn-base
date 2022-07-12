////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of doubled spin 1/2 chain with translations                //
//                                                                            //
//                        last modification : 01/05/2016                      //
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


#ifndef DOUBLEDSPIN1_2_CHAINWITHTRANSLATIONS_ALTERNATIVE_H
#define DOUBLEDSPIN1_2_CHAINWITHTRANSLATIONS_ALTERNATIVE_H


#include "config.h"
#include "HilbertSpace/AbstractDoubledSpinChainWithTranslations.h"


#include <iostream>


using std::ostream;

class DoubledSpin1_2_ChainWithTranslations_alternative : public AbstractDoubledSpinChainWithTranslations
{
  
 protected:
  int ShiftNegativeDiffSz;
  int BraShiftNegativeSz;
  int * PowerD;
  unsigned long * ChainDescription;
  Complex * TranslationPhase;
    
 public:

  // default constructor
  //
  DoubledSpin1_2_ChainWithTranslations_alternative ();

 
  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  DoubledSpin1_2_ChainWithTranslations_alternative (int chainLength, int momentum, int sz, int memorySize, int memorySlice);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  DoubledSpin1_2_ChainWithTranslations_alternative (const DoubledSpin1_2_ChainWithTranslations_alternative & chain);

  // destructor
  //
  ~DoubledSpin1_2_ChainWithTranslations_alternative ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  DoubledSpin1_2_ChainWithTranslations_alternative & operator = (const DoubledSpin1_2_ChainWithTranslations_alternative & chain);

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

  // evaluate a density matrix of a subsystem of the whole system described by a given real ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem (return a zero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int szSector, RealVector& groundState);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem (return a zero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int szSector, ComplexVector& groundState);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz and momentum sector.
  // 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // momentumSector = momentum sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem (return a zero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int szSector, int momentumSector, ComplexVector& groundState);

  
  virtual void AddConvertFromGeneralSpace(ComplexVector vSource,ComplexVector & vDestination);
  virtual void ConvertToGeneralSpace(ComplexVector vSource,ComplexVector & vDestination);
  
  virtual void AddConvertFromGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination);
  virtual void ConvertToGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination);
  
 protected:

  // return value of twice spin projection on (Oz) for a given state
  //
  // stateDescription = state to which the spin projection has to be evaluated
  // return value = twice spin projection on (Oz)
  int GetTotalSz (unsigned long stateDescriptionBra, unsigned long stateDescriptionKet);
  
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

  // generate all states corresponding to the constraints
  // 
  // lengthBra = length of the chain to be decided for bra spins
  // lengthBra = length of the chain to be decided for ket spins
  // diffSz = difference of spin projection between bra and ket chain
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int lengthBra, int lengthKet, int diffSz, long pos);
  
  long ShiftedEvaluateHilbertSpaceDimension(int lengthBra, int lengthKet, int diffSz);

  inline unsigned int GetCommonIndexFromBraAndKetIndices(unsigned int braIndex, unsigned int ketIndex );
  inline void GetBraAndKetIndicesFromCommonIndex(unsigned int & braIndex, unsigned int & ketIndex, unsigned long commonIndex);
  
  void GenerateLookUpTable(unsigned long memory);


  // evaluate all exponential factors
  //   
  virtual void EvaluateExponentialFactors();
  
};


// find the canonical form of a state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// return value = canonical form of the state

inline void DoubledSpin1_2_ChainWithTranslations_alternative::FindCanonicalForm(unsigned long stateDescription, unsigned long & canonicalState, int& nbrTranslation)
{
  nbrTranslation = 0;
  canonicalState = stateDescription;
  int index = 1;  
  while (index < this->ChainLength)
    {
      stateDescription = (stateDescription / 4) + (stateDescription % 4)*this->PowerD[this->ChainLength-1];
      
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

inline void DoubledSpin1_2_ChainWithTranslations_alternative::FindCanonicalForm(unsigned long stateDescription, unsigned long & canonicalState, int& nbrTranslation, int& nbrTranslationToIdentity)
{
  nbrTranslation = 0;
  nbrTranslationToIdentity = 1;
  canonicalState = stateDescription;

  unsigned long ReferenceState = stateDescription;

  stateDescription = (stateDescription / 4) + (stateDescription % 4)*this->PowerD[this->ChainLength-1];
  while ((ReferenceState != stateDescription) && (nbrTranslationToIdentity < this->ChainLength))
    {
      if (stateDescription < canonicalState)
	{
	  canonicalState = stateDescription;
	  nbrTranslation = nbrTranslationToIdentity;
	}
      
      stateDescription = (stateDescription / 4) + (stateDescription % 4)*this->PowerD[this->ChainLength-1];
      ++nbrTranslationToIdentity;
    }
}

// find how many translations are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

inline int DoubledSpin1_2_ChainWithTranslations_alternative::FindNumberTranslation(unsigned long stateDescription)
{
  int index = 1;
  unsigned long TmpState = (stateDescription / 4) + (stateDescription % 4)*this->PowerD[this->ChainLength-1];
  while (TmpState !=  stateDescription)
    {     
      TmpState = (TmpState / 4) + (TmpState % 4)*this->PowerD[this->ChainLength-1];
      ++index;
    }
  return index;
}



inline unsigned int DoubledSpin1_2_ChainWithTranslations_alternative::GetCommonIndexFromBraAndKetIndices(unsigned int braIndex, unsigned int ketIndex )
{
  return ketIndex * 2 + braIndex;
}


inline void DoubledSpin1_2_ChainWithTranslations_alternative::GetBraAndKetIndicesFromCommonIndex(unsigned int & braIndex, unsigned int & ketIndex, unsigned long commonIndex)
{
  braIndex = commonIndex % 2;
  ketIndex = commonIndex / 2;
}

#endif


