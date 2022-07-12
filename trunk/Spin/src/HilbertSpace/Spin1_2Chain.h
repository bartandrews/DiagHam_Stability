////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of spin 1/2 chain                          //
//                                                                            //
//                        last modification : 05/03/2001                      //
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


#ifndef SPIN1_2CHAIN_H
#define SPIN1_2CHAIN_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChain.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;


class Spin1_2Chain : public AbstractSpinChain
{
  
  friend class DoubledSpin1_2_ChainWithTranslations_alternative;
  friend class Spin1_2ChainWithTranslationsAndSzSymmetry;
  friend class Spin1_2ChainWithTranslationsAndInversionSymmetry;
  friend class Spin1_2ChainWithTranslations;

 protected:

  int Sz;
  bool FixedQuantumNumberFlag;
  
  int* LookUpTable;
  unsigned long LookUpTableMask;
  int LookUpPosition;
  int LookUpTableSize;

  // array describing each n-nody state 
  unsigned long* StateDescription;

 public:


  // default constructor
  //
  Spin1_2Chain ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin 1/2
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2Chain (int chainLength, int memorySize);

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin 1/2
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2Chain (int chainLength, int sz, int memorySize);

  // constructor from pre-constructed datas
  //
  // hilbertSpaceDimension = Hilbert space dimension
  // chainDescription = array describing states
  // chainLength = number of spin 1/2
  // sz = twice the value of total Sz component
  // fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
  // lookUpTable = look-up table
  // lookUpTableSize = look-Up table size
  // lookUpTablePosition = last position described by the look-Up table
  // lookUpTableMask = look-Up table mask
  Spin1_2Chain (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
		int sz, bool fixedQuantumNumberFlag, int* lookUpTable, int lookUpTableSize, 
		int lookUpPosition, unsigned long lookUpTableMask);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2Chain (const Spin1_2Chain& chain);

  // destructor
  //
  ~Spin1_2Chain ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2Chain& operator = (const Spin1_2Chain& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // re-initialize chain with another total Sz component
  //
  // sz = twice the value of total Sz component
  // return value = reference on current chain  
  virtual Spin1_2Chain& Reinitialize(int sz);

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
  Matrix& Syi (int i, Matrix& M);

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

  // compute the parity (prod_i Sz_i) for a given state
  //
  // state = index of the state to be applied on Sz_i operator
  // return value = total Sz value
  virtual unsigned long Parity (int state);

  // return index of resulting state from application of P_ij operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to be applied on P_ij operator
  // return value = index of resulting state
  virtual int Pij (int i, int j, int state);

  // return index of resulting state from application of a 3 sites permutation operator on a given state
  //
  // i = first position
  // j = second position (j > i)
  // k = third position  (k > j)
  // state = index of the state to be applied on P_ijk operator
  // return value = index of resulting state
  virtual int Pijk (int i, int j, int k, int state);

  // return index of resulting state from application of a 3 sites permutation inverse operator on a given state
  //
  // i = first position
  // j = second position (j > i)
  // k = third position  (k > j)
  // state = index of the state to be applied on P_ijk operator
  // return value = index of resulting state
  virtual int Pminusijk (int i, int j, int k, int state);

  // return eigenvalue of Sz_i Sz_j associated to a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  virtual double SziSzj (int i, int j, int state);
  
  // return the eigenvalue the product of consecutive Sz_i's 
  //
  // indexMin = index of the leftmost site 
  // indexMax = index of the righttmost site 
  // state = index of the state to consider
  // return value = corresponding eigenvalue (either -1 or +1)
  virtual int ProdSzj (int indexMin, int indexMax, int state);

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

  // translate a state assuming the system have periodic boundary
  // conditions (increasing the site index)
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

  // find state index
  //
  // state = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long state);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition.
  // 
  // nbrSpinUp = number of spin up that belong to the subsytem 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrSpinUpSector, RealVector& groundState);

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

  // get the normalization factor in front of each basis state (i.e. 1/sqrt(orbit size))
  //
  // return value = pointer to normalization factors
  virtual double* GetBasisNormalization();
 
 protected:

  // generate Spin 1/2 states
  //
  // statePosition = position for the new states
  // sitePosition = site on chain where spin has to be changed
  // currentStateDescription = description of current state
  // return value = number of generated states
  virtual int GenerateStates(int statePosition, int sitePosition, unsigned currentStateDescription);

  // generate Spin 1/2 states for a given total spin projection Sz
  //
  // statePosition = position for the new states
  // sitePosition = site on chain where spin has to be changed
  // currentStateDescription = description of current state
  // currentSz = total Sz value of current state
  // return value = number of generated states
  virtual int GenerateStates(int statePosition, int sitePosition, unsigned currentStateDescription, int currentSz);

  // evaluate Hilbert space dimension
  //
  // nbrSpins = number of spins
  // sz = twice the z projection of the total momentum
  // return value = Hilbert space dimension
  virtual int EvaluateHilbertSpaceDimension(int nbrSpins, int szMax);

};
 
// get the value of the spin (i.e. S) at a given site
// 
// site = site index
// return value = twice the spin

inline int Spin1_2Chain::GetLocalSpin(int site)
{
  return 1;
}

#endif


