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


#ifndef VIRTUALSPACETRANSFERMATRIXWITHTRANSLATIONS_H
#define VIRTUALSPACETRANSFERMATRIXWITHTRANSLATIONS_H 


#include "config.h"
#include "HilbertSpace/AbstractDoubledSpinChainWithTranslations.h"


#include <iostream>


using std::ostream;

class VirtualSpaceTransferMatrixWithTranslations : public AbstractDoubledSpinChainWithTranslations
{
  
 protected:
  int BondDimension;
  int * PowerD;
  unsigned long * ChainDescription;
  Complex * TranslationPhase;
   
 public:

  // default constructor
  //
  VirtualSpaceTransferMatrixWithTranslations ();

  // constructor for complete Hilbert space with no restriction on momentum
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // memorySize = memory size in bytes allowed for look-up table
  // memorySlice = maximum amount of memory that can be allocated to partially evalauted the states
  VirtualSpaceTransferMatrixWithTranslations (int chainLength, int bondDimension, int memorySize, int memorySlice);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  VirtualSpaceTransferMatrixWithTranslations (int chainLength, int bondDimension,  int momentum,  int translationStep, int memorySize, int memorySlice);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  VirtualSpaceTransferMatrixWithTranslations (const VirtualSpaceTransferMatrixWithTranslations & chain);

  // destructor
  //
  ~VirtualSpaceTransferMatrixWithTranslations ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  VirtualSpaceTransferMatrixWithTranslations & operator = (const VirtualSpaceTransferMatrixWithTranslations & chain);

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

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (RealVector& groundState);
  virtual HermitianMatrix EvaluatePartialDensityMatrix (ComplexVector& groundState);
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int momentumSector, ComplexVector& groundState);

  virtual void AddConvertFromGeneralSpace(ComplexVector vSource,ComplexVector & vDestination);
  virtual void ConvertToGeneralSpace(ComplexVector vSource,ComplexVector & vDestination);
  
  virtual void AddConvertFromGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination);
  virtual void ConvertToGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination);
//  virtual void ApplyInversionSymmetry(ComplexVector & sourceVector,  ComplexVector & destinationVector);  
//  virtual void ApplyInversionSymmetry(ComplexVector & sourceVector,  ComplexVector & destinationVector, bool translationFlag); 

  virtual void NormalizeDensityMatrix(ComplexVector & sourceVector);
 protected:

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
  long GenerateStates(int length, long pos);
  
  long ShiftedEvaluateHilbertSpaceDimension(int length);

  inline unsigned int GetCommonIndexFromBraAndKetIndices(unsigned int braIndex, unsigned int ketIndex );
  inline void GetBraAndKetIndicesFromCommonIndex(unsigned int & braIndex, unsigned int & ketIndex, unsigned long commonIndex);
  
  void GenerateLookUpTable(unsigned long memory);


  // evaluate all exponential factors
  //   
  virtual void EvaluateExponentialFactors();

  virtual int FindNextInversionSymetricIndice(int);
  virtual  inline void ApplySingleXTranslation(unsigned long& stateDescription);

};


// find the canonical form of a state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// return value = canonical form of the state

inline void VirtualSpaceTransferMatrixWithTranslations::FindCanonicalForm(unsigned long stateDescription, unsigned long & canonicalState, int& nbrTranslation)
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

inline void VirtualSpaceTransferMatrixWithTranslations::FindCanonicalForm(unsigned long stateDescription, unsigned long & canonicalState, int& nbrTranslation, int& nbrTranslationToIdentity)
{
  nbrTranslation = 0;
  nbrTranslationToIdentity = 1;
  canonicalState = stateDescription;

  unsigned long ReferenceState = stateDescription;

  this->ApplySingleXTranslation(stateDescription);
  while ((ReferenceState != stateDescription) && (nbrTranslationToIdentity < this->ChainLength))
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

inline int VirtualSpaceTransferMatrixWithTranslations::FindNumberTranslation(unsigned long stateDescription)
{
  int index = 1;
  unsigned long TmpState = stateDescription;
  this->ApplySingleXTranslation(TmpState);
  while (TmpState !=  stateDescription)
    {     
      this->ApplySingleXTranslation(TmpState);
      ++index;
    }
  return index;
}


// apply a single translation in the x direction for a state description
//
// stateDescription = reference on the state description

inline void VirtualSpaceTransferMatrixWithTranslations::ApplySingleXTranslation(unsigned long& stateDescription)
{
  for(int i =0 ; i < this->ChainLength/this->MaxXMomentum; i++)
    stateDescription = (stateDescription/this->PowerD[1]) + (stateDescription%this->PowerD[1])*this->PowerD[this->ChainLength-1];
}

inline unsigned int VirtualSpaceTransferMatrixWithTranslations::GetCommonIndexFromBraAndKetIndices(unsigned int braIndex, unsigned int ketIndex )
{
  return ketIndex * this->BondDimension + braIndex;
}


inline void VirtualSpaceTransferMatrixWithTranslations::GetBraAndKetIndicesFromCommonIndex(unsigned int & braIndex, unsigned int & ketIndex, unsigned long commonIndex)
{
  braIndex = commonIndex%this->BondDimension;
  ketIndex = commonIndex/this->BondDimension;
}

#endif


