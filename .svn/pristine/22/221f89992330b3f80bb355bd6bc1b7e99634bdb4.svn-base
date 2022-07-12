////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of abstract doubled spin chain with translations            //
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


#ifndef ABSTRACTDOUBLEDSPINCHAINWITHTRANSLATIONS_H
#define ABSTRACTDOUBLEDSPINCHAINWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractDoubledSpinChain.h"

#include <iostream>


using std::ostream;

class AbstractDoubledSpinChainWithTranslations : public AbstractDoubledSpinChain
{
  
  friend class ComplexPEPSTransfertMatrixPBCWithTranslations;
  friend class ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations;
  
 protected:

  //  total difference sz component (if fixed)
  int Momentum;
  int MaxXMomentum;
  
  // array containing flag indicating if a state belonging to an orbit with a given number of member is compatible with momentum constraint
  bool* CompatibilityWithMomentum;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;

  // number of state in each orbit
  int* NbrStateInOrbit;

  // shift to apply to move the spin from one end to the other one
  int ComplementaryStateShift;

 public:

  // destructor
  //
  ~AbstractDoubledSpinChainWithTranslations ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  AbstractDoubledSpinChainWithTranslations& operator = (const AbstractDoubledSpinChainWithTranslations & chain);

  virtual void AddConvertFromGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination){};
  virtual void ConvertToGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination){};

 protected:

  // find the canonical form of a state
  //
  // state = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // return value = canonical form of the state
  virtual void FindCanonicalForm ( unsigned long stateDescriptionBra, unsigned long stateDescriptionKet, unsigned long & canonicalStateBra, unsigned long & canonicalStateKet, int& nbrTranslation){};

  // find the canonical form of a state and find how many translations are needed to obtain the same state
  //
  // stateDescription = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
  // return value = canonical form of the state
  virtual void FindCanonicalForm ( unsigned long stateDescriptionBra, unsigned long stateDescriptionKet, unsigned long & canonicalStateBra, unsigned long & canonicalStateKet, int& nbrTranslation, int& nbrTranslationToIdentity){};
  
  // find how many translations are needed to obtain the same state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of translation needed to obtain the same state
  virtual int FindNumberTranslation(unsigned long stateDescriptionBra, unsigned long stateDescriptionKet){};

  // create precalculation tables
  //
  void CreatePrecalculationTable();
  

  // generate all states corresponding to the constraints
  // 
  // lengthBra = length of the chain to be decided for bra spins
  // lengthBra = length of the chain to be decided for ket spins
  // diffSz = difference of spin projection between bra and ket chain
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int lengthBra, int lengthKet, int diffSz, long pos);
  
  inline int GetMomentum(){return this->Momentum;};


  //return the scaling factor when going from state i to state j
  inline double GetRescalingFactor(int i,int j) const {return this->RescalingFactors[this->NbrStateInOrbit[i]][this->NbrStateInOrbit[j]];};
  
  unsigned long FindCanonicalForm(unsigned long state, int& nbrTranslation) {cout<< "Using undefined function unsigned long FindCanonicalForm(unsigned long state, int& nbrTranslation) in DoubledSpin0_1_2_chainWithTranslations"<<endl; return 0; };

};



#endif


