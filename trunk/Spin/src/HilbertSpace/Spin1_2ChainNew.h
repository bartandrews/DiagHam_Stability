////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with a fixed Sz value                //
//                                                                            //
//                        last modification : 29/06/2015                      //
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


#ifndef SPIN1_2CHAINNEW_H
#define SPIN1_2CHAINNEW_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;


class Spin1_2ChainNew : public AbstractSpinChain
{

 friend class Spin1_2ChainNewAnd2DTranslation;
 friend class Spin1_2ChainNewSzSymmetryAnd2DTranslation;
 friend class Spin1_2ChainNewInversionAnd2DTranslation;
 friend class Spin1_2ChainNewSzSymmetryInversionAnd2DTranslation;
 friend class Spin1_2ChainWithPseudospin;
 friend class Spin1_2ChainWithPseudospinAnd2DTranslation;
 friend class Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation;
 friend class ComplexPEPSPBC;
 
 protected:

  // twice the Sz value
  int Sz;
  // true if Sz is set for the Hilbert space
  bool FixedQuantumNumberFlag;

  // array describing each n-nody state 
  unsigned long* StateDescription;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table
  unsigned long LookUpTableMemorySize;
  // shift used for configurations with a spin up at a given largest position
  int* LookUpTableShift;
  // look-up table
  int** LookUpTable;


 public:


  // default constructor
  //
  Spin1_2ChainNew ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin 1/2
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2ChainNew (int chainLength, int sz, int memorySize);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainNew (const Spin1_2ChainNew& chain);

  // destructor
  //
  ~Spin1_2ChainNew ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainNew& operator = (const Spin1_2ChainNew& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // get the value of the spin (i.e. S) at a given site
  // 
  // site = site index
  // return value = twice the spin
  virtual int GetLocalSpin(int site);

  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  virtual int TotalSz (int index);

  // return index of resulting state from application of S+_i operator on a given state
  //
  // i = position of S+ operator
  // state = index of the state to be applied on S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int Spi (int i, int state, double& coefficient);

  // return index of resulting state from application of S-_i operator on a given state
  //
  // i = position of S- operator
  // state = index of the state to be applied on S-_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int Smi (int i, int state, double& coefficient);

  // return index of resulting state from application of Sz_i operator on a given state
  //
  // i = position of Sz operator
  // state = index of the state to be applied on Sz_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int Szi (int i, int state, double& coefficient);

  // compute the parity (prod_i Sz_i) for a given state
  //
  // state = index of the state to be applied on Sz_i operator
  // return value = total Sz value
  virtual unsigned long GetParity (int state);

  // return index of resulting state from application of P_ij operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to be applied on P_ij operator
  // return value = index of resulting state
  virtual int Pij (int i, int j, int state);

  // return eigenvalue of Sz_i Sz_j associated to a given state
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
  virtual int SmiSpj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i S-_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of S- operator
  // state = index of the state to be applied on S+_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSmj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i S+_j operator on a given state
  //
  // i = position of first S+ operator
  // j = position of second S+ operator
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSpj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S-_i S-_j operator on a given state
  //
  // i = position of first S- operator
  // j = position of second S- operator
  // state = index of the state to be applied on S-_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSmj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i Sz_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of Sz operator
  // state = index of the state to be applied on S+_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSzj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S-_i Sz_j operator on a given state
  //
  // i = position of S- operator
  // j = position of Sz operator
  // state = index of the state to be applied on S-_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSzj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i S-_j Sz_k operator on a given state
  //
  // i = position of S+ operator
  // j = position of S- operator
  // k = position of Sz operator
  // state = index of the state to be applied on S+_i S-_j Sz_k operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSmjSzk (int i, int j, int k, int state, double& coefficient);
  
  // operate local isometry on three sites
  //
  // i = position of first site
  // j = position of second site 
  // k = position of third site
  // state = index of the state that the isometry has to be applied on
  // indices = reference to an array where the indices of the resulting states have to be stored
  // coefficients = reference to the array where the coefficients have to be stored
  // return value = number of non-zero coefficients
  virtual int ThreeSiteIsometry (int i, int j, int k, int state, int*& indices, double*& coefficients);
  

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

  // generate all states
  // 
  // nbrSpinUp = number of spin up
  // currentPosition = current position to consider in the chain
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrSpinUp, int currentPosition, long pos);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  // stateMask = an optional mask to apply to each state to focus on the relevant bits
#ifdef __64_BITS__
  virtual void GenerateLookUpTable(unsigned long memory, unsigned long stateMask = 0xfffffffffffffffful);
#else
  virtual void GenerateLookUpTable(unsigned long memory, unsigned long stateMask = 0xfffffffful);
#endif


  // convert the state on the site to its binary representation
  //
  // state = state to be stored
  // sitePosition = position on the chain of the state
  // return integer that code the state
  virtual unsigned long EncodeSiteState(int physicalState, int sitePosition);


  // return the Bosonic Occupation of a given state in the basis
  //
  // index = index of the state in the basis
  // finalState = reference on the array where the monomial representation has to be stored
  virtual void GetBosonicOccupation (unsigned int index, int * finalState);
  
};

// get the value of the spin (i.e. S) at a given site
// 
// site = site index
// return value = twice the spin

inline int Spin1_2ChainNew::GetLocalSpin(int site)
{
  return 1;
}

#endif


