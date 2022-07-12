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
#include "Matrix/HermitianMatrix.h"

#include <iostream>


using std::ostream;


class Spin1_2ChainWithTranslations : public AbstractSpinChainWithTranslations
{
  friend class DoubledSpin1_2_ChainWithTranslations_alternative;

 protected: 

  // number of sites in the x direction
  int MaxXMomentum;

  // array containing falg indicating if a state beloging to an orbit with a given number of member is compatible with momentum constraint
  bool* CompatibilityWithMomentum;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

  // flag to indicate if the total sz component is fixed
  bool FixedSpinProjectionFlag;
  //  total sz component (if fixed)
  int Sz;
  
  // shift to apply to move a group spin using an elementary translation
  int StateShift;
  // shift to apply to move the spin from one end to the other one
  int ComplementaryStateShift;
  // mask to apply to get a group of spins
  unsigned long StateMask;

  // shift to apply to a state to obtain an index to the look-up table 
  int LookUpTableShift;
  // look-up table (LookUpTable[i] gives the index of the smallest state that greater than i <<  LookUpTableShift)
  long* LookUpTable;

  int ShiftLookUpTable;

  // array describing each state
  unsigned long* StateDescription;

 public:

  // default constructor
  //
  Spin1_2ChainWithTranslations ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // translationStep = indicates the step for an elementary translation
  // memorySize = memory size in bytes allowed for look-up table
  // memorySlice = maximum amount of memory that can be allocated to partially evaluate the states
  Spin1_2ChainWithTranslations (int chainLength, int momentum, int translationStep, int memorySize, int memorySlice);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // translationStep = indicates the step for an elementary translation
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2ChainWithTranslations (int chainLength, int momentum, int translationStep, int sz, int memorySize, int memorySlice);

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
  virtual AbstractHilbertSpace* Clone();

  // get the normalization factor in front of each basis state (i.e. 1/sqrt(orbit size))
  //
  // return value = pointer to normalization factors
  virtual double* GetBasisNormalization();
  
  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index);

  // return value of the value of the sum of the square of spin projection on (Oz) 
  //
  // index = index of the state to test
  // return value = twice spin projection on (Oz)
  virtual double TotalSzSz (int index);

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

  // return eigenvalue of Sz_i associated to a given state
  //
  // i = position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  virtual double Szi (int i, int state);

  // return eigenvalue of Sz_i Sz_j associated to a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  virtual double SziSzj (int i, int j, int state);

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
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  virtual int Pij (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of a 3 sites permutation operator on a given state
  //
  // i = first position
  // j = second position (j > i)
  // k = third position  (k > j)
  // state = index of the state to be applied on P_ijk operator
  // coefficient = reference on the numerical coefficient
  // nbrTranslations = reference on the number of translations to apply to the resulting state to obtain the canonical state
  // return value = index of resulting state
  virtual int Pijk (int i, int j, int k, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of a 3 sites permutation inverse operator on a given state
  //
  // i = first position
  // j = second position (j > i)
  // k = third position  (k > j)
  // state = index of the state to be applied on P_ijk operator
  // coefficient = reference on the numerical coefficient
  // nbrTranslations = reference on the number of translations to apply to the resulting state to obtain the canonical state
  // return value = index of resulting state
  virtual int Pminusijk (int i, int j, int kt, int state, double& coefficient, int& nbrTranslation);

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

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

  ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. Sz is not conserved
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

  ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, ComplexVector& groundState, AbstractArchitecture* architecture);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition.
  // 
  // nbrSpinUp = number of spin up that belong to the subsytem 
  // kSector = momentum of the subsystem
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrSpinUpSector, int kSector, RealVector& groundState);
  

 protected:

  // return value of twice spin projection on (Oz) for a given state
  //
  // stateDescription = state to which the spin projection has to be evaluated
  // return value = twice spin projection on (Oz)
  virtual int GetTotalSz (unsigned long stateDescription);

  // find the canonical form of a state
  //
  // state = state description
  // nbrTranslations = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // return value = canonical form of the state
  virtual unsigned long FindCanonicalForm(unsigned long state, int& nbrTranslations);

  // find the canonical form of a state and find how many translations are needed to obtain the same state
  //
  // stateDescription = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
  // return value = canonical form of the state
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation, int& nbrTranslationToIdentity);

  // find how many translations are needed to obtain the same state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of translation needed to obtain the same state
  virtual int FindNumberTranslation(unsigned long stateDescription);

  // constructor from pre-constructed datas
  //
  // hilbertSpaceDimension = Hilbert space dimension
  // chainDescription = array describing states
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // sz = twice the value of total Sz component
  // fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
  // lookUpTableShift = shift to apply to a state to obtain an index to the look-up table 
  // complementaryStateShift = shift to apply to move the spin from one end to the other one
  Spin1_2ChainWithTranslations (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
				int momentum, int sz, bool fixedQuantumNumberFlag, int lookUpTableShift, 
				int complementaryStateShift);

  // generate all states with constraint on total Sz
  //
  // position = current position in the Hilbert space basis
  // sitePosition = largest sit position that has to be filled 
  // currentNbrSpinUp = number of spin to put in each state
  // return value = new current position in the Hilbert space basis
  virtual long GenerateStates(long position, int sitePosition, int currentNbrSpinUp);

  // create precalculation tables
  //
  virtual void CreatePrecalculationTable();

  // create look-up table used to speed up index search
  //
  virtual void CreateLookUpTable();

  // evaluate Hilbert space dimension
  //
  // nbrSpins = number of spins
  // sz = twice the z projection of the total momentum
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrSpins, int szMax);

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
  
  //return the scaling factor when going from state i to state j
  virtual double GetRescalingFactor(int i,int j) const {return this->RescalingFactors[this->NbrStateInOrbit[i]][this->NbrStateInOrbit[j]];};
  
  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // nbrStateInOrbit = original number of states in the orbit before the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, int& nbrTranslations);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(unsigned long stateDescription);

  // apply a single translation in the x direction for a state description
  //
  // stateDescription = reference on the state description
  virtual void ApplySingleXTranslation(unsigned long& stateDescription);

};

// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// nbrStateInOrbit = original number of states in the orbit before the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int Spin1_2ChainWithTranslations::SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, 
							  int& nbrTranslations)
{
  state = this->FindCanonicalForm(state, nbrTranslations);
  int TmpIndex = this->FindStateIndex(state);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[nbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslations = (this->MaxXMomentum - nbrTranslations) % this->MaxXMomentum;
     }
  return TmpIndex;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int Spin1_2ChainWithTranslations::FindOrbitSize(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  return XSize;
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool Spin1_2ChainWithTranslations::TestMomentumConstraint(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);   
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  if (((this->Momentum * XSize) % this->MaxXMomentum) != 0)
    return false;
  return true;
}


// find the canonical form of a state
//
// stateDescription = state description
// nbrTranslations = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// return value = canonical form of the state

inline unsigned long Spin1_2ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations)
{
  nbrTranslations = 0;
  unsigned long CanonicalState = stateDescription;
  nbrTranslations = 0;
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(stateDescription);      
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslations = n;	      
	}
    }
  return CanonicalState;
}

// find the canonical form of a state and find how many translations are needed to obtain the same state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
// return value = canonical form of the state

inline unsigned long Spin1_2ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation, int& nbrTranslationToIdentity)
{
  nbrTranslation = 0;
  nbrTranslationToIdentity = 1;
  unsigned long CanonicalState = stateDescription;
  unsigned long ReferenceState = stateDescription;
  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->StateMask) << this->ComplementaryStateShift);
  while ((ReferenceState != stateDescription) && (nbrTranslationToIdentity < this->ChainLength))
    {
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = nbrTranslationToIdentity;
	}
      stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->StateMask) << this->ComplementaryStateShift);
      ++nbrTranslationToIdentity;
    }
  return CanonicalState;
}

// find how many translations are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

inline int Spin1_2ChainWithTranslations::FindNumberTranslation(unsigned long stateDescription)
{
  unsigned long TmpState = (stateDescription >> this->StateShift) | ((stateDescription & this->StateMask) << this->ComplementaryStateShift);
  int index = 1;  
  while (TmpState != stateDescription)
    {
      TmpState = (TmpState >> this->StateShift) | ((TmpState & this->StateMask) << this->ComplementaryStateShift);
      ++index;
    }
  return index;
}


// apply a single translation in the x direction for a state description
//
// stateDescription = reference on the state description

inline void Spin1_2ChainWithTranslations::ApplySingleXTranslation(unsigned long& stateDescription)
{
  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->StateMask) << this->ComplementaryStateShift);
}

// get the value of the spin (i.e. S) at a given site
// 
// site = site index
// return value = twice the spin

inline int Spin1_2ChainWithTranslations::GetLocalSpin(int site)
{
  return 1;
}

#endif


