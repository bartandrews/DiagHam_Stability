////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of spin 1 chain with one 1/2 spin at the left end of the chain    //
//                                                                            //
//                        last modification : 23/02/2001                      //
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


#ifndef SPIN1AKLTCHAIN_H
#define SPIN1AKLTCHAIN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class Spin1AKLTChain : public AbstractSpinChain
{

 protected:

  int Sz;
  bool FixedQuantumNumberFlag;
  
  unsigned long* ChainDescription;

  int* LookUpTable;
  unsigned long LookUpTableMask;
  int LookUpPosition;
  int LookUpTableSize;

 public:

  // default constructor
  //
  Spin1AKLTChain ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // memorySize = memory size in bytes allowed for look-up table
  Spin1AKLTChain (int chainLength, int memorySize);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin1AKLTChain (int chainLength, int sz, int memorySize) ;

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1AKLTChain (const Spin1AKLTChain& chain);

  // destructor
  //
  ~Spin1AKLTChain ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1AKLTChain& operator = (const Spin1AKLTChain& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // return Hilbert space dimension
  //
  // return value = Hilbert space dimension
  int GetHilbertSpaceDimension();

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

  // return matrix representation of Sx
  //
  // i = operator position
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  Matrix& Sxi (int i, Matrix& M);

  // return matrix representation of i * Sy
  //
  // i = operator position
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  Matrix& Syi (int i, Matrix& M);

  // return matrix representation of Sz
  //
  // i = operator position
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  Matrix& Szi (int i, Matrix& M);

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

  // return index of resulting state from application of S-_i S-_j operator on a given state
  //
  // i = position of first S- operator
  // j = position of second S- operator
  // state = index of the state to be applied on S-_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int SmiSmj (int i, int j, int state, double& coefficient);

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

 private:

  // constructor from pre-constructed datas
  //
  // hilbertSpaceDimension = Hilbert space dimension
  // chainDescription = array describing states
  // chainLength = number of spin 1
  // sz = twice the value of total Sz component
  // fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
  // lookUpTable = look-up table
  // lookUpTableSize = look-Up table size
  // lookUpTablePosition = last position described by the look-Up table
  // lookUpTableMask = look-Up table mask
  Spin1AKLTChain (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
		  int sz, bool fixedQuantumNumberFlag, int* lookUpTable, int lookUpTableSize, 
		  int lookUpPosition, unsigned long lookUpTableMask);
  
  // generate all states
  //
  // statePosition = position for the new states
  // sitePosition = site on chain where spin has to be changed
  // currentStateDescription = description of current state
  // return value = number of generated states
  int GenerateStates(int statePosition, int sitePosition, unsigned long currentStateDescription);

  // generate all states corresponding to a given total Sz
  //
  // statePosition = position for the new states
  // sitePosition = site on chain where spin has to be changed
  // currentStateDescription = description of current state
  // currentSz = total Sz value of current state
  // return value = number of generated states
  int GenerateStates(int statePosition,int sitePosition, unsigned long currentStateDescription, int currentSz);

};

#endif


