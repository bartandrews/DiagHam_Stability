////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of Potts-3 chain                           //
//                                                                            //
//                        last modification : 04/06/2012                      //
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


#ifndef POTTS3CHAIN_H
#define POTTS3CHAIN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


using std::cout;
using std::endl;


class Potts3Chain : public AbstractSpinChain
{


  friend class Potts3ChainWithTranslations;
  friend class Potts3ChainWithTranslationsAndInversion;
  
 protected:

  // total quantum number for the chain (i.e. either 0, 1 or 2)
  int Sz;
  // flag to indicate if Sz is fixed 
  bool FixedQuantumNumberFlag;
  
  // look-up table
  int* LookUpTable;
  // look-up table size
  int LookUpTableSize;
  // shift to apply to a state to get the key in th look-up table
  int LookUpTableShift;

  // array containing the state descriptions
  unsigned long* StateDescription;

public:


  // default constructor
  //
  Potts3Chain ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // memorySize = memory size in bytes allowed for look-up table
  Potts3Chain (int chainLength, int memorySize);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Potts3Chain (int chainLength, int sz, int memorySize) ;

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Potts3Chain (const Potts3Chain& chain);

  // destructor
  //
  ~Potts3Chain ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Potts3Chain& operator = (const Potts3Chain& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  int TotalSz (int index);

  // return index of resulting state from application of S+_i operator on a given state
  //
  // i = position of S+ operator
  // state = index of the state to be applied on S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int Spi (int i, int state, double& coefficient);

  // return index of resulting state from application of S-_i operator on a given state
  //
  // i = position of S- operator
  // state = index of the state to be applied on S-_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int Smi (int i, int state, double& coefficient);

  // return index of resulting state from application of Sz_i operator on a given state
  //
  // i = position of Sz operator
  // state = index of the state to be applied on Sz_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int Szi (int i, int state, double& coefficient);

  // return eigenvalue of Sz_i Sz_j associated to a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  double SziSzj (int i, int j, int state);

  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int SmiSpj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i S+_j operator on a given state
  //
  // i = position of first S+ operator
  // j = position of second S+ operator
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int SpiSpj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i S+_j S+_k operator on a given state
  //
  // i = position of the first S+ operator
  // j = position of the second S+ operator
  // k = position of the third S+ operator
  // state = index of the state to be applied on S+_i S+_j S+_k operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int SpiSpjSpk (int i, int j, int k, int state, double& coefficient);

  // return index of resulting state from application of S-_i S-_j operator on a given state
  //
  // i = position of first S- operator
  // j = position of second S- operator
  // state = index of the state to be applied on S-_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int SmiSmj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S-_i S-_j S-_k operator on a given state
  //
  // i = position of the first S- operator
  // j = position of the second S- operator
  // k = position of the third S- operator
  // state = index of the state to be applied on S-_i S-_j S-_k operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int SmiSmjSmk (int i, int j, int k, int state, double& coefficient);

  // return index of resulting state from application of S+_i Sz_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of Sz operator
  // state = index of the state to be applied on S+_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int SpiSzj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S-_i Sz_j operator on a given state
  //
  // i = position of S- operator
  // j = position of Sz operator
  // state = index of the state to be applied on S-_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int SmiSzj (int i, int j, int state, double& coefficient);

  // return value the product of the tau_j operators (Q operator) on a given state
  //
  // index = index of the state 
  // return value = Q value (0 for 1, 1 for exp(i 2 \pi / 3), -1 1 for exp(i 2 \pi / 3)) 
  double QValue (int index);

  // translate a state assuming the system have periodic boundary conditions (increasing the site index)
  //
  // nbrTranslations = number of translations to apply
  // state = index of the state to translate 
  // return value = index of resulting state
  int TranslateState (int nbrTranslations, int state);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter);

  // find state index
  //
  // state = state description
  // return value = corresponding index
  int FindStateIndex(unsigned long state);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

 protected:
  
  // constructor from pre-constructed datas
  //
  // largeHilbertSpaceDimension = Hilbert space dimension
  // chainDescription = array describing states
  // chainLength = number of spin 1
  // sz = twice the value of total Sz component
  // fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
  // lookUpTable = look-up table
  // lookUpTableSize = look-Up table size
  // lookUpTableShift = shift to apply to a state to get the key in th look-up table
  Potts3Chain (long largeHilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
	       int sz, bool fixedQuantumNumberFlag, int* lookUpTable, int lookUpTableSize, 
	       int lookUpTableShift);
  
  // evaluate Hilbert space dimension
  //
  // currentSite = current site to occupy
  // currentSzValue = state current Sz value 
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int currentSite, int currentSzValue);

  // generate all states with no constraint on total Sz
  //
  // currentSite = current site to occupy
  // currentPosition = current position of the state that has to be considered
  // return value = number of generated states
  long GenerateStates(int currentSite, long currentPosition);

  // generate all states corresponding to a given total Sz
  //
  // currentSite = current site to occupy
  // currentSzValue = state current Sz value 
  // currentPosition = current position of the state that has to be considered
  // return value = number of generated states
  long GenerateStates(int currentSite, int currentSzValue, long currentPosition);

  
};

// return index of resulting state from application of Sz_i operator on a given state
//
// i = position of Sz operator
// state = index of the state to be applied on Sz_i operator
// coefficient = reference on double where numerical coefficient has to be stored (0 for 1.0, 1.0 for exp(i 2 \pi / 3), 2.0 for exp(i 4 \pi / 3)) 
// return value = index of resulting state

inline int Potts3Chain::Szi (int i, int state, double& coefficient)
{
  //  coefficient = (double) ((this->StateDescription[state] >> (i << 1)) & 0x3ul);
  coefficient = (((double) ((this->StateDescription[state] >> (i << 1)) & 0x1ul))
		 - ((double) ((this->StateDescription[state] >> ((i << 1) + 1)) & 0x1ul)));
  return state;
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

inline double Potts3Chain::SziSzj (int i, int j, int state)
{  
  return ((double) (((this->StateDescription[state] >> (i << 1)) & 0x3ul)
		    * ((this->StateDescription[state] >> (j << 1)) & 0x3ul)));
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

inline int Potts3Chain::SmiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = tmpState;
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << j);
  switch (tmpState2)
    {
    case 0x1ul:
      tmpState |= 0x2ul << j;
      break;
    case 0x0ul:
      tmpState |= 0x1ul << j;
      break;
    }	  
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x2ul:
      tmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x2ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return value the product of the tau_j operators (Q operator) on a given state
//
// index = index of the state 
// return value = Q value (0 for 1, 1 for exp(i 2 \pi / 3), -1 1 for exp(i 2 \pi / 3)) 

inline double Potts3Chain::QValue (int index)
{
  unsigned long Tmp = 0x0ul;
  unsigned long Tmp2 = this->StateDescription[index];
  for (int i = 0; i < this->ChainLength; ++i)
    {
      Tmp += Tmp2 & 0x3ul;
      Tmp2 >>= 2;
    }
  return ((double) (Tmp % 3));
}

// find state index
//
// state = state description
// return value = corresponding index

inline int Potts3Chain::FindStateIndex(unsigned long state)
{
  return SearchInArrayDownOrdering<unsigned long> (state, this->StateDescription, 
						   this->HilbertSpaceDimension);
/*   int TmpShift = this->LookUpTable[state >> this->LookUpTableShift]; */
/*   return (TmpShift + SearchInArrayDownOrdering<unsigned long> (state, this->StateDescription + TmpShift, */
/* 							       this->HilbertSpaceDimension - TmpShift)); */
}

#endif


