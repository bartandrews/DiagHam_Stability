////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of doubled spin 1/2 chain                      //
//                                                                            //
//                        last modification : 02/02/2016                      //
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


#ifndef _DOUBLEDSPIN1_2_CHAIN_H
#define _DOUBLEDSPIN1_2_CHAIN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class DoubledSpin1_2_Chain : public AbstractSpinChain
{

 protected:

  //  total difference sz component (if fixed)
  int DiffSz;

  // flag to indicate if the total sz component is fixed
  bool FixedSpinProjectionFlag;
  
  // shift to apply to move the spin from one end to the other one
  int ComplementaryStateShift;

  // shift to apply to a state to obtain an index to the look-up table 
  int LookUpTableShift;
  // look-up table (LookUpTable[i] gives the index of the smallest state that greater than i <<  LookUpTableShift)
  long* LookUpTable;

  // array describing each state
  unsigned long* ChainDescriptionBra;
  unsigned long* ChainDescriptionKet;


  // sorted array that contains each unique configuration for the type up particles
  unsigned long* UniqueStateDescriptionBra;
  // number of time each unique configuration for the type up particles appears in StateDescriptionUp
  int* UniqueStateDescriptionSubArraySizeBra;
  // number of unique configurations for the type up-plus particles
  long NbrUniqueStateDescriptionBra;
  // first time a type up appears in the Hilbert space
  int* FirstIndexUniqueStateDescriptionBra;


 public:

  // default constructor
  //
  DoubledSpin1_2_Chain ();

  // constructor for complete Hilbert space with no restriction on momentum
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // memorySize = memory size in bytes allowed for look-up table
  // memorySlice = maximum amount of memory that can be allocated to partially evalauted the states
  DoubledSpin1_2_Chain (int chainLength,  int diffSz, int memorySize, int memorySlice);


  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
   DoubledSpin1_2_Chain (const  DoubledSpin1_2_Chain & chain);

  // destructor
  //
  ~DoubledSpin1_2_Chain ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  DoubledSpin1_2_Chain & operator = (const  DoubledSpin1_2_Chain & chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  int TotalSz (int index);
  
  // find state index
  //
  // state = state description
  // return value = corresponding index
  int FindStateIndex(unsigned long stateBra, unsigned long stateKet);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // return the Bosonic Occupation of a given state in the basis
  //
  // index = index of the state in the basis
  // return value bosonic occupation 
  inline void GetBosonicOccupation (unsigned int index, int * finalStateBra,int * finalStateKet);

  // convert the state on the site to its binary representation
  //
  // state = state to be stored
  // sitePosition = position on the chain of the state
  // return integer that code the state
  inline unsigned long EncodeSiteStateBra(int physicalState, int sitePosition);

  // convert the state on the site to its binary representation
  //
  // state = state to be stored
  // sitePosition = position on the chain of the state
  // return integer that code the state
  inline unsigned long EncodeSiteStateKet(int physicalState, int sitePosition);

 private:

  // return value of twice spin projection on (Oz) for a given state
  //
  // stateDescription = state to which the spin projection has to be evaluated
  // return value = twice spin projection on (Oz)
  int GetTotalSz  (unsigned long stateDescriptionBra,unsigned long stateDescriptionKet);
  
  // generate all states corresponding to the constraints
  // 
  // lengthBra = length of the chain to be decided for bra spins
  // lengthBra = length of the chain to be decided for ket spins
  // diffSz = difference of spin projection between bra and ket chain
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int lengthBra, int lengthKet, int diffSz, long pos);
  
  // create look-up table used to speed up index search
  //
  void GenerateLookUpTable(unsigned long memory);
  
  long ShiftedEvaluateHilbertSpaceDimension(int lengthBra, int lengthKet, int diffSz);

  double TotalSzSz (int index);
  double SziSzj (int i, int j, int state);

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
  
  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state 
  virtual int SmiSpj (int i, int j, int state, double& coefficient);
  
  // return index of resulting state from application of S+_i operator on a given state
  //
  // i = position of S+ operator
  // state = index of the state to be applied on S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int Spi (int i, int state, double& coefficient);

};

// return the Bosonic Occupation of a given state in the basis
//
// index = index of the state in the basis
// finalState = reference on the array where the monomial representation has to be stored

inline void DoubledSpin1_2_Chain::GetBosonicOccupation (unsigned int index, int * finalStateBra,int * finalStateKet)
{
  for (unsigned long i = 0; i < this->ChainLength; i++)
    {
      finalStateBra[i] = (this->ChainDescriptionBra[index] >> i ) & 0x1ul;
      finalStateKet[i] = (this->ChainDescriptionKet[index] >> i ) & 0x1ul;
    }
}

// convert the state on the site to its binary representation
//
// state = state to be stored
// sitePosition = position on the chain of the state
// return integer that code the state

inline unsigned long DoubledSpin1_2_Chain::EncodeSiteStateBra(int physicalState, int sitePosition)
{
  return  physicalState << sitePosition;
}

// convert the state on the site to its binary representation
//
// state = state to be stored
// sitePosition = position on the chain of the state
// return integer that code the state

inline unsigned long DoubledSpin1_2_Chain::EncodeSiteStateKet(int physicalState, int sitePosition)
{
  return  physicalState << sitePosition;
}

#endif


