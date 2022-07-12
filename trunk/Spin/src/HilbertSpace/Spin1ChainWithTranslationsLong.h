////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of spin 1 chain with translations for more than 32 spins      //
//                                                                            //
//                        last modification : 18/03/2019                      //
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


#ifndef SPIN1CHAINWITHTRANSLATIONSLONG_H
#define SPIN1CHAINWITHTRANSLATIONSLONG_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class Spin1ChainWithTranslationsLong : public AbstractSpinChainWithTranslations
{

 protected:

  // number of sites in the x direction
  int MaxXMomentum;

  // bit shift that has to applied to perform a translation in the x direction 
  int StateXShift;
  // binary mask for the StateXShift first bits 
  ULONGLONG XMomentumMask;
  // bit shift to apply to move the first StateXShift bits at the end of a state description
  int ComplementaryStateXShift;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

  // flag to indicate if the total sz component is fixed
  bool FixedSpinProjectionFlag;
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
  ULONGLONG* StateDescription;

 public:

  // default constructor
  //
  Spin1ChainWithTranslationsLong ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // memory = amount of memory granted for precalculations
  Spin1ChainWithTranslationsLong (int chainLength, int momentum, unsigned long memory = 10000000);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // sz = twice the value of total Sz component
  // memory = amount of memory granted for precalculations
  Spin1ChainWithTranslationsLong (int chainLength, int momentum, int sz, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1ChainWithTranslationsLong (const Spin1ChainWithTranslationsLong& chain);

  // destructor
  //
  ~Spin1ChainWithTranslationsLong ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1ChainWithTranslationsLong& operator = (const Spin1ChainWithTranslationsLong& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

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

  // get the momentum of each state in the current Hilbert space
  //
  // return value = momentum value
  virtual int GetMomentum();

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

  // find state index
  //
  // stateDescription = state description
  // maxBitPosition = maximum bit set to one in stateDescription
  // return value = corresponding index
  virtual int FindStateIndex(ULONGLONG stateDescription);

  // find state index
  //
  // stateDescription = state description
  // maxBitPosition = maximum bit set to one in stateDescription
  // return value = corresponding index
  virtual int FindStateIndex(ULONGLONG stateDescription, int maxBitPosition);

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

  // convert a state defined in the real space basis into a state in the (Kx,Ky) basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertToKxKyBasis(ComplexVector& state, AbstractSpinChain* space);

  // convert a state defined in the (Kx,Ky) basis into a state in the real space basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space);

 protected:

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // nbrStateInOrbit = original number of states in the orbit before the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeResult(ULONGLONG& state, int nbrStateInOrbit, double& coefficient, int& nbrTranslations);

  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual ULONGLONG FindCanonicalForm(ULONGLONG stateDescription, int& nbrTranslations);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(ULONGLONG stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(ULONGLONG stateDescription);

  // apply a single translation in the x direction for a state description
  //
  // stateDescription = reference on the state description
  virtual void ApplySingleXTranslation(ULONGLONG& stateDescription);

  // return value of twice spin projection on (Oz) for a given state
  //
  // stateDescription = state to which the spin projection has to be evaluated
  // return value = twice spin projection on (Oz)
  virtual int GetTotalSz (ULONGLONG stateDescription);

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

  // generate all states with constraints 
  //
  // return value = number of generated states
  virtual long GenerateStates();

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);
  
  // compute the rescaling factors
  //
  virtual void ComputeRescalingFactors();

};

// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// nbrStateInOrbit = original number of states in the orbit before the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int Spin1ChainWithTranslationsLong::SymmetrizeResult(ULONGLONG& state, int nbrStateInOrbit, double& coefficient, 
							    int& nbrTranslations)
{
  state = this->FindCanonicalForm(state, nbrTranslations);
  int TmpMaxMomentum = 2 * this->ChainLength;
  while (((state >> TmpMaxMomentum) == ((ULONGLONG) 0x0ul)) && (TmpMaxMomentum > 0))
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[nbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslations = (this->MaxXMomentum - nbrTranslations) % this->MaxXMomentum;
     }
  return TmpIndex;
}

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline ULONGLONG Spin1ChainWithTranslationsLong::FindCanonicalForm(ULONGLONG stateDescription, int& nbrTranslations)
{
  ULONGLONG CanonicalState = stateDescription;
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

// find the size of the orbit for a given state
//
// return value = orbit size

inline int Spin1ChainWithTranslationsLong::FindOrbitSize(ULONGLONG stateDescription)
{
  ULONGLONG TmpStateDescription = stateDescription;
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

inline bool Spin1ChainWithTranslationsLong::TestMomentumConstraint(ULONGLONG stateDescription)
{
  ULONGLONG TmpStateDescription = stateDescription;
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

// apply a single translation in the x direction for a state description
//
// stateDescription = reference on the state description

inline void Spin1ChainWithTranslationsLong::ApplySingleXTranslation(ULONGLONG& stateDescription)
{
  stateDescription = (stateDescription >> this->StateXShift) | ((stateDescription & this->XMomentumMask) << this->ComplementaryStateXShift);
}

// get the momentum of each state in the current Hilbert space
//
// return value = momentum value

inline int Spin1ChainWithTranslationsLong::GetMomentum()
{
  return this->Momentum;
}

// get the value of the spin (i.e. S) at a given site
// 
// site = site index
// return value = twice the spin

inline int Spin1ChainWithTranslationsLong::GetLocalSpin(int site)
{
  return 2;
}

#endif


