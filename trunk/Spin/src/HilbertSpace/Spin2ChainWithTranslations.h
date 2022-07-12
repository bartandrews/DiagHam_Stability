////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of spin 2 chain with translations                  //
//                                                                            //
//                        last modification : 09/11/2016                      //
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


#ifndef SPIN2CHAINWITHTRANSLATIONS_H
#define SPIN2CHAINWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class Spin2ChainWithTranslations : public Spin1ChainWithTranslations
{

 protected:

  // number of coefficients per component for the spin 3 projector
  int* Spin3ProjectorNbrCoefficients;
  // number of coefficients per component for the spin 4 projector
  int* Spin4ProjectorNbrCoefficients;
  // coefficients per component for the spin 3 projector
  double** Spin3ProjectorCoefficients;
  // coefficients per component for the spin 4 projector
  double** Spin4ProjectorCoefficients;
  // configuration per component for the spin 3 projector
  unsigned long** Spin3ProjectorStates;
  // configuration per component for the spin 4 projector
  unsigned long** Spin4ProjectorStates;
  

 public:

  // default constructor
  //
  Spin2ChainWithTranslations ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // memory = amount of memory granted for precalculations
  Spin2ChainWithTranslations (int chainLength, int momentum, unsigned long memory = 10000000);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // sz = twice the value of total Sz component
  // memory = amount of memory granted for precalculations
  Spin2ChainWithTranslations (int chainLength, int momentum, int sz, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin2ChainWithTranslations (const Spin2ChainWithTranslations& chain);

  // destructor
  //
  ~Spin2ChainWithTranslations ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin2ChainWithTranslations& operator = (const Spin2ChainWithTranslations& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // get the normalization factor in front of each basis state (i.e. 1/sqrt(orbit size))
  //
  // return value = pointer to normalization factors
  virtual double* GetBasisNormalization();
  
  // get the value of the spin (i.e. S) at a given site
  // 
  // site = site index
  // return value = twice the spin
  virtual int GetLocalSpin(int site);

  // return value of the value of the sum of the square of spin projection on (Oz) 
  //
  // index = index of the state to test
  // return value = twice spin projection on (Oz)
  virtual double TotalSzSz (int index);

  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  virtual int TotalSz (int index);

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
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S+_i operator on a given state (only valid if there is no constraint on total Sz)
  //
  // i = operator position
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int Spi (int i, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i operator on a given state (only valid if there is no constraint on total Sz)
  //
  // i = operator position
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int Smi (int i, int state, double& coefficient, int& nbrTranslation);
    
  // return index of resulting state from application of S+_i S+_j operator on a given state
  //
  // i = position of first S+ operator
  // j = position of second S+ operator
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i S-_j operator on a given state
  //
  // i = position of first S- operator
  // j = position of second S- operator
  // state = index of the state to be applied on S-_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S+_i S+_i operator on a given state
  //
  // i = position of first S+ operator
  // state = index of the state to be applied on S+_i S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SpiSpi (int i, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i S-_i operator on a given state
  //
  // i = position of the S- operator
  // state = index of the state to be applied on S-_i S-_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SmiSmi (int i, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S+_i Sz_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of Sz operator
  // state = index of the state to be applied on S+_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i Sz_j operator on a given state
  //
  // i = position of S- operator
  // j = position of Sz operator
  // state = index of the state to be applied on S-_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i1 S+_j1 S-_i2 S+_j2 operator on a given state
  //
  // i1 = position of leftmost S- operator
  // j1 = position of leftmost S+ operator
  // i2 = position of rightmost S- operator
  // j2 = position of rightmost S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state (orbit index)
  virtual int SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of Sz_i1 Sz_j1 S-_i2 S+_j2 operator on a given state
  //
  // i1 = position of first Sz operator
  // j1 = position of second Sz operator
  // i2 = position of S- operator
  // j2 = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state (orbit index)
  virtual int SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient, int& nbrTranslation);

  // compute all the states connected to a single one by a two site spin 3 projector
  //
  // i = index of the first site
  // j = index of the first site
  // state = index of state on whcih the projector should be applied
  // indices = pointer to the array where the connected state indeices will be stored
  // coefficients = pointer to the array where coefficients related to each connected states will be stored
  // nbrTranslations  = pointer to the array where number of translations related to each connected states will be stored
  // return value = number of connected states
  virtual int Spin3Projector (int  i, int j, int state, int* indices, double* coefficients, int* nbrTranslations);

  // compute all the states connected to a single one by a two site spin 4 projector
  //
  // i = index of the first site
  // j = index of the first site
  // state = index of state on whcih the projector should be applied
  // indices = pointer to the array where the connected state indeices will be stored
  // coefficients = pointer to the array where coefficients related to each connected states will be stored
  // nbrTranslations  = pointer to the array where number of translations related to each connected states will be stored
  // return value = number of connected states
  virtual int Spin4Projector (int  i, int j, int state, int* indices, double* coefficients, int* nbrTranslations);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

 protected:

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // nbrStateInOrbit = original number of states in the orbit before the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, int& nbrTranslations);

  // return value of twice spin projection on (Oz) for a given state
  //
  // stateDescription = state to which the spin projection has to be evaluated
  // return value = twice spin projection on (Oz)
  virtual int GetTotalSz (unsigned long stateDescription);

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
  
  // build the projectors of two spins on a given total spin sector
  //
  virtual void BuildProjectors();

  // compute all the states connected to a single one by a two site spin projector
  //
  // i = index of the first site
  // j = index of the first site
  // state = index of state on whcih the projector should be applied
  // indices = pointer to the array where the connected state indeices will be stored
  // coefficients = pointer to the array where coefficients related to each connected states will be stored
  // nbrTranslations  = pointer to the array where number of translations related to each connected states will be stored
  // spinProjectorNbrCoefficients = number of coefficients per component for the spin projector
  // spinProjectorCoefficients = coefficients per component for the spin projector
  // spinProjectorStates = configuration per component for the spin projector
  // return value = number of connected states
  virtual int SpinProjector (int  i, int j, int state, int* indices, double* coefficients, int* nbrTranslations,
			     int* spinProjectorNbrCoefficients,  double** spinProjectorCoefficients, unsigned long** spinProjectorStates);

};

// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// nbrStateInOrbit = original number of states in the orbit before the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int Spin2ChainWithTranslations::SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, 
							int& nbrTranslations)
{
  state = this->FindCanonicalForm(state, nbrTranslations);
  int TmpMaxMomentum = 3 * this->ChainLength;
  while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[nbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslations = (this->MaxXMomentum - nbrTranslations) % this->MaxXMomentum;
     }
  return TmpIndex;
}

// get the value of the spin (i.e. S) at a given site
// 
// site = site index
// return value = twice the spin

inline int Spin2ChainWithTranslations::GetLocalSpin(int site)
{
  return 4;
}

// compute all the states connected to a single one by a two site spin 3 projector
//
// i = index of the first site
// j = index of the first site
// state = index of state on whcih the projector should be applied
// indices = pointer to the array where the connected state indeices will be stored
// coefficients = pointer to the array where coefficients related to each connected states will be stored
// nbrTranslations  = pointer to the array where number of translations related to each connected states will be stored
// return value = number of connected states

inline int Spin2ChainWithTranslations::Spin3Projector (int  i, int j, int state, int* indices, double* coefficients, int* nbrTranslations)
{
  return this->SpinProjector(i, j, state, indices, coefficients, nbrTranslations, this->Spin3ProjectorNbrCoefficients, 
			     this->Spin3ProjectorCoefficients, this->Spin3ProjectorStates);
}

// compute all the states connected to a single one by a two site spin 4 projector
//
// i = index of the first site
// j = index of the first site
// state = index of state on whcih the projector should be applied
// indices = pointer to the array where the connected state indeices will be stored
// coefficients = pointer to the array where coefficients related to each connected states will be stored
// nbrTranslations  = pointer to the array where number of translations related to each connected states will be stored
// return value = number of connected states

inline int Spin2ChainWithTranslations::Spin4Projector (int  i, int j, int state, int* indices, double* coefficients, int* nbrTranslations)
{
  return this->SpinProjector(i, j, state, indices, coefficients, nbrTranslations, this->Spin4ProjectorNbrCoefficients, 
			     this->Spin4ProjectorCoefficients, this->Spin4ProjectorStates);
}


#endif


