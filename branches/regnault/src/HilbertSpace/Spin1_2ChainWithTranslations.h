////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin 1/2 chain with translastion invariance          //
//                                                                            //
//                        last modification : 29/01/2002                      //
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


#ifndef SPIN1_2CHAINWITHTRANSLATIONS_H
#define SPIN1_2CHAINWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"

#include <iostream>


using std::ostream;


class Spin1_2ChainWithTranslations : public AbstractSpinChainWithTranslations
{

 protected:

  int ChainLength;
  int ReducedChainLength;
  int HilbertSpaceDimension;

  int Sz;
  bool FixedSpinProjectionFlag;
  
  int Momentum;
  bool FixedMomentumFlag;
  
  int* LookUpTable;
  int LookUpTableShift;
  int LookUpTableSize;

  unsigned long* ChainDescription;
  int* OrbitSize;

 public:

  // default constructor
  //
  Spin1_2ChainWithTranslations ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin 1/2
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2ChainWithTranslations (int chainLength, int memorySize);

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin 1/2
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2ChainWithTranslations (int chainLength, int sz, int memorySize);

  // constructor for complete Hilbert space with a given total spin projection Sz an a given momentum 
  //
  // chainLength = number of spin 1/2
  // sz = twice the value of total Sz component
  // momentum = momentum
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2ChainWithTranslations (int chainLength, int sz, int momentum, int memorySize);
    
  // constructor from pre-constructed datas
  //
  // hilbertSpaceDimension = Hilbert space dimension
  // chainDescription = array describing states
  // orbitSize = number of state per orbit
  // chainLength = number of spin 1/2
  // sz = twice the value of total Sz component
  // fixedSpinPorjectionFlag = true if hilbert space is restricted to a given spin projection
  // momentum = momentum
  // fixedSpinPorjectionFlag = true if hilbert space is restricted to a given momentum
  // lookUpTableSize = look-Up table size
  Spin1_2ChainWithTranslations (int hilbertSpaceDimension, unsigned long* chainDescription, 
				int* orbitSize, int chainLength, 
				int sz, bool fixedSpinProjectionFlag, 
				int momentum, bool fixedMomentumFlag, 
				int lookUpTableSize);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainWithTranslations (const Spin1_2ChainWithTranslations& chain);

  // destructor
  //
  ~Spin1_2ChainWithTranslations ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainWithTranslations& operator = (const Spin1_2ChainWithTranslations& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // re-initialize chain with another total Sz component
  //
  // sz = twice the value of total Sz component
  // return value = reference on current chain  
  Spin1_2ChainWithTranslations& Reinitialize(int sz);

  // return Hilbert space dimension
  //
  // return value = Hilbert space dimension
  int GetHilbertSpaceDimension() {return this->HilbertSpaceDimension;};

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

  // return index of resulting state from application of P_ij operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to be applied on P_ij operator
  // return value = index of resulting state
  int Pij (int i, int j, int state);

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
  // translations = number of translations to apply to the resulting state to obtain the true resulting state
  // return value = index of resulting state (orbit index)
  int SmiSpj (int i, int j, int state, double& coefficient, int& translations);

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

  // find canonical state
  //
  // state = state description
  // nbrTranslation = reference on an integer used to store number of translations
  // return value = canonical state description
  unsigned long FindCanonicalState(unsigned long state, int& nbrTranslation);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

 private:

  // generate Spin 1/2 states
  //
  // sitePosition = site on chain where spin has to be changed
  // currentStateDescription = description of current state
  // return value = number of generated states
  void GenerateStates(int sitePosition, unsigned long currentStateDescription);

  // generate Spin 1/2 states for a given total spin projection Sz
  //
  // sitePosition = site on chain where spin has to be changed
  // currentStateDescription = description of current state
  // currentSz = total Sz value of current state
  // return value = number of generated states
  void GenerateStates(int sitePosition, unsigned long currentStateDescription, int currentSz);

  // store a state in its canonical form
  //
  // stateDescription = description of the stateto store
  void StoreCanonicalForm (unsigned long stateDescription);

  // build look-up table
  //
  // memorySize = memory size in bytes allowed for look-up table
  void BuildLookUpTable(int memorySize = 0);

  // Apply momentum criteria to keep allowed state 
  //
  void ApplyMomentumCriteria ();

};

#endif


