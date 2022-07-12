////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of spin 3/2 chain                          //
//                                                                            //
//                        last modification : 08/07/2018                      //
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


#ifndef SPIN3_2CHAIN_H
#define SPIN3_2CHAIN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Vector/RealVector.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class ComplexMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class Spin3_2Chain : public AbstractSpinChain
{

  friend class Spin3_2ChainWithTranslations;
  friend class Spin3_2ChainWithTranslationsAndSzSymmetry;
  friend class Spin3_2ChainWithTranslationsAndInversionSymmetry;
  friend class Spin3_2ChainAnd2DTranslation;
  friend class Spin3_2ChainSzSymmetryAnd2DTranslation;
  
 protected:

  // true if the Hilbert is defined with a specific total Sz value
  bool FixedQuantumNumberFlag;
  //  total sz component (if fixed)
  int Sz;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given maxMomentum sector
  int LookUpTableMemorySize;
  // shift used in each maxMomentum sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used maxMomentum value of the state an the second 
  int** LookUpTable;

  // array describing each state
  unsigned long* StateDescription;

 public:


  // default constructor
  //
  Spin3_2Chain ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // memorySize = memory size in bytes allowed for look-up table
  Spin3_2Chain (int chainLength, int memorySize);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin3_2Chain (int chainLength, int sz, int memorySize) ;

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin3_2Chain (const Spin3_2Chain& chain);

  // destructor
  //
  ~Spin3_2Chain ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin3_2Chain& operator = (const Spin3_2Chain& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // get the value of the spin (i.e. S) at a given site
  // 
  // site = site index
  // return value = twice the spin
  virtual int GetLocalSpin(int site);

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index);

  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  virtual int TotalSz (int index);

  // return matrix representation of Sx
  //
  // i = operator position
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  virtual Matrix& Sxi (int i, Matrix& M);

  // return matrix representation of i * Sy
  //
  // i = operator position
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  virtual Matrix& Syi (int i, Matrix& M);

  // return matrix representation of Sz
  //
  // i = operator position
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  virtual Matrix& Szi (int i, Matrix& M);

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

  // return index of resulting state from application of S-_i S+_j Sz_k operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // k = position of Sz operator
  // state = index of the state to be applied on S-_i S+_j Sz_k operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSpjSzk (int i, int j, int k, int state, double& coefficient);

  // return index of resulting state from application of S-_i1 S+_j1 S-_i2 S+_j2 operator on a given state
  //
  // i1 = position of leftmost S- operator
  // j1 = position of leftmost S+ operator
  // i2 = position of rightmost S- operator
  // j2 = position of rightmost S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient);

  // return index of resulting state from application of Sz_i1 Sz_j1 S-_i2 S+_j2 operator on a given state
  //
  // i1 = position of first Sz operator
  // j1 = position of second Sz operator
  // i2 = position of S- operator
  // j2 = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient);

  // translate a state assuming the system have periodic boundary conditions (increasing the site index)
  //
  // nbrTranslations = number of translations to apply
  // state = index of the state to translate 
  // return value = index of resulting state
  virtual int TranslateState (int nbrTranslations, int state);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture = 0);
	
  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);


 protected:


  // find state index
  //
  // stateDescription = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription);

  // find state index
  //
  // stateDescription = state description
  // maxBitPosition = maximum bit set to one in stateDescription
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int maxBitPosition);

  // evaluate Hilbert space dimension with no constraint on the total Sz
  //
  // nbrSites = number of sites
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrSites);

  // evaluate Hilbert space dimension
  //
  // sz = twice the Sz value
  // nbrSites = number of sites
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int sz, int nbrSites);

  // generate all states with no constraint on total Sz and no discrete symmtry constraint
  //
  // statePosition = position for the new states
  // sitePosition = site on chain where spin has to be changed
  // currentStateDescription = description of current state
  // return value = number of generated states
  virtual long RawGenerateStates(long statePosition, int sitePosition);

  // generate all states corresponding to a given total Sz and no discrete symmtry constraint
  //
  // statePosition = position for the new states
  // sitePosition = site on chain where spin has to be changed
  // currentSz = total Sz value of current state
  // return value = number of generated states
  virtual long RawGenerateStates(long statePosition, int sitePosition, int currentSz); 

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);
  
};

// get the value of the spin (i.e. S) at a given site
// 
// site = site index
// return value = twice the spin

inline int Spin3_2Chain::GetLocalSpin(int site)
{
  return 3;
}

// find state index
//
// stateDescription = state description
// return value = corresponding index

inline int Spin3_2Chain::FindStateIndex(unsigned long stateDescription)
{
  int TmpMaxBitPosition= 3 * this->ChainLength;
  while (((stateDescription >> TmpMaxBitPosition) == 0x0ul) && (TmpMaxBitPosition > 0))
    {
      --TmpMaxBitPosition;
    }
  return this->FindStateIndex(stateDescription, TmpMaxBitPosition);
}

#endif


