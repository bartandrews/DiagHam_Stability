////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of fermion on a torus taking into account magnetic translations    //
//                                                                            //
//                        last modification : 10/09/2003                      //
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


#ifndef FERMIONONTORUSWITHMAGNETICTRANSLATIONS_H
#define FERMIONONTORUSWITHMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnSphere.h"

using std::cout;
using std::endl;
using std::dec;
using std::hex;


class FermionOnTorusWithMagneticTranslations :  public ParticleOnTorusWithMagneticTranslations
{

 protected:

  // number of fermions
  int NbrFermions;
  // number of fermions plus 1
  int IncNbrFermions;

  // maximum momentum value reached by a fermion
  int MaxMomentum;
  // number of momentum values in a state (= MaxMomentum +1)
  int NbrMomentum;
  // GCD of MaxMomentum and NbrFermions (momemta are defined modulo MomentumModulo)
  int MomentumModulo;
  // momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
  int XMomentum;
  // momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)
  int YMomentum;

  // value that has to be substracted to the momentum for each translation of the canonical form research
  int MomentumIncrement;
  // shift that has to be done on a state for each translation of the canonical form research
  int StateShift;
  // complementary shift (with respect to MaxMomentum) to StateShift
  int ComplementaryStateShift;
  // mask corresponding to StateShift
  unsigned long MomentumMask;

  // shift to apply to a state before inverting its expression
  int InvertShift;
  // shift to apply to a state after inverting its expression
  int InvertUnshift;

  // array describing each state 
  unsigned long* StateDescription;
  // array giving maximum momentum value reached for a fermion in a given state
  int* StateMaxMomentum;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given maxMomentum sector
  int LookUpTableMemorySize;
  // shift used in each maxMomentum sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used maxMomentum value of the state an the second 
  int** LookUpTable;

  // a table containing ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  unsigned long* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;
  // a table containing parity of the sum of 1 bits for all integer ranging from 0 to 2^MaximumSignLookUp - 1 (1 if odd)
  int* NbrParticleLookUpTable;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

  // temporary variables when using AdAd / ProdAd operations
  unsigned long ProdATemporaryState;
  int ProdATemporaryStateMaxMomentum;
  int ProdATemporaryNbrStateInOrbit;

  // sign due to state reordering when applying translation operator 
  unsigned long* ReorderingSign;
  // array of unsigned long where each bit describes sign associated to each translation of the orbit representant (0 for +, 1 for -) with respect to N-body ordering convention
//  int* StateSignature;

  // array containing for each state the sign due to fermion reordering when translating state (1 bit to 0 if sign is negative)
//  unsigned long* TranslationSign;

 public:

  // default constructor
  // 
  FermionOnTorusWithMagneticTranslations ();

  // basic constructor
  // 
  // nbrFermions = number of fermions 
  // maxMomentum = momentum maximum value for a fermion
  // xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)
  FermionOnTorusWithMagneticTranslations (int nbrFermions, int maxMomentum, int xMomentum, int yMomentum);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnTorusWithMagneticTranslations(const FermionOnTorusWithMagneticTranslations& fermions);

  // destructor
  //
  ~FermionOnTorusWithMagneticTranslations ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnTorusWithMagneticTranslations& operator = (const FermionOnTorusWithMagneticTranslations& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

  // get momemtum value in the y direction of a given state
  //
  // index = state index
  // return value = state momentum in the y direction
  int GetYMomentumValue(int index);

  // get the number of orbitals
  //
  // return value = number of orbitals
  virtual int GetNbrOrbitals();

  // get the number of particles
  //
  // return value = number of particles
  virtual int GetNbrParticles();

  // get the momentum along the x axis
  // 
  // return avlue = momentum along the x axis
  virtual int GetKxMomentum();

  // get the momentum along the y axis
  // 
  // return avlue = momentum along the y axis
  virtual int GetKyMomentum();

  // get the maximum momentum along the x axis (i.e. the number of momentum sectors)
  // 
  // return avlue = maximum momentum along the x axis
  virtual int GetMaxXMomentum();
  
  // get the maximum momentum along the y axis (i.e. the number of momentum sectors)
  // 
  // return avlue = maximum momentum along the y axis
  virtual int GetMaxYMomentum();
  
  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
					 SubspaceSpaceConverter& converter);

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation);

  // apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int n1, int n2);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int nbrIndices);

  // apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int nbrIndices, double& coefficient, int& nbrTranslation);

  // apply a^+_m a_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for creation operator
  // return value =  resulting multiplicative factor 
  virtual double AdA (int index, int m);
  
  // apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A (int index, int n);
 
  // convert a state defined in the Ky basis into a state in the (Kx,Ky) basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertToKxKyBasis(ComplexVector& state, ParticleOnTorus* space);

  // convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertFromKxKyBasis(ComplexVector& state, ParticleOnTorus* space);
  
  // convert a state defined in the (Kx,Ky) basis into a state in the Ky basis, and change its kx quantum number
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // oldKx = vnew value of the relative quantum number kx
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertFromKxKyBasisAndModifyKx(ComplexVector& state, ParticleOnTorus* space, int oldKx);

  // get the C2 symmetric state of a given state 
  //
  // index = index of the state whose symmetric counterpart has to be computed
  // nbrTranslation = number of translations that has to be applied to C2 symmetric state to be in the canonical form
  virtual int GetC2SymmetricState (int index, int& nbrTranslation);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // lzSector = Lz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // request whether state with given index satisfies a general Pauli exclusion principle
  //
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  // return value = true if teh state satisfies the general Pauli exclusion principle
  bool HasPauliExclusions(int index, int pauliK, int pauliR);

 protected:

  // find canonical form of a state description
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
  // nbrTranslation = number of translation needed to obtain the canonical form
  // yMomentum = state momentum value in the y direction
  // return value = canonical form of a state description
  unsigned long FindCanonicalForm(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation);

  // safe version to find canonical form of a state description (with additional checks)
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
  // nbrTranslation = number of translation needed to obtain the canonical form
  // yMomentum = state momentum value in the y direction
  // return value = canonical form of a state description
  unsigned long SafeFindCanonicalForm(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation);

  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = canonical form of a state description and -1 in nbrTranslation if the state does not fit the x momentum constraint
  unsigned long FindCanonicalFormAndTestXMomentumConstraint(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation);

 // find how many translations on the x direction are needed to obtain the same state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of translation needed to obtain the same state
  int FindNumberXTranslation(unsigned long stateDescription);

  // test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = maximum momentum value that can be reached by a fermion in the stateDescription state
  // return value = true if the state satisfy the x momentum constraint
  bool TestXMomentumConstraint(unsigned long stateDescription, int maxMomentum);

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription, int maxMomentum);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrFermions, int maxMomentum);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(int memory);

  // generate look-up table associated to sign calculations
  // 
  void GenerateSignLookUpTable();

  // generate all states corresponding to the constraints
  // tmpDimension = max dimension of Hilbert space (to be reduced by symmetries)
  // return value = hilbert space dimension
  long GenerateStates(long tmpDimension);

  // generate all states corresponding to the constraints (without taking into the canonical form) 
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion in the state
  // currentMaxMomentum = momentum maximum value for fermions that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentYMomentum = current value of the momentum in the y direction
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, long pos, int currentUMomentum);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in complementary Hilbert space
  // nbrIndex = number of indices to consider in complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorusWithMagneticTranslations* complementaryHilbertSpace,  ParticleOnTorusWithMagneticTranslations* destinationHilbertSpace,
								  ComplexVector& groundState,  HermitianMatrix* densityMatrix);

  // core part of the C4 rotation
  //
  // inputState = reference on the state that has to be rotated
  // inputSpace = Hilbert space associated to the input state
  // outputState = reference on the rotated state
  // minIndex = minimum index that has to be computed
  // nbrIndices = number of indices that have to be computed
  // clockwise = the rotation is done clockwise
  // return value = reference on the rotated state
  virtual ComplexVector& CoreC4Rotation (ComplexVector& inputState, ParticleOnTorusWithMagneticTranslations* inputSpace, 
					 ComplexVector& outputState, int minIndex, int nbrIndices, bool clockwise);

  // convert a fermionic state to its monomial representation
  //
  // initialState = initial fermionic state in its fermionic representation
  // finalState = reference on the array where the monomial representation has to be stored
  virtual void ConvertToMonomial(unsigned long initialState, unsigned long*& finalState);

  // convert a fermionic state from its monomial representation
  //
  // initialState = array where the monomial representation is stored
  // return value = fermionic state in its fermionic representation
  virtual unsigned long ConvertFromMonomial(unsigned long* initialState);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnTorusWithMagneticTranslations::GetParticleStatistic()
{
  return ParticleOnTorusWithMagneticTranslations::FermionicStatistic;
}

// get the number of orbitals
//
// return value = number of orbitals

inline int FermionOnTorusWithMagneticTranslations::GetNbrOrbitals()
{
  return this->MaxMomentum;
}

// get the number of particles
//
// return value = number of particles

inline int FermionOnTorusWithMagneticTranslations::GetNbrParticles()
{
  return this->NbrFermions;
}

// get the momentum along the x axis
// 
// return avlue = momentum along the x axis

inline int FermionOnTorusWithMagneticTranslations::GetKxMomentum()
{
  return this->XMomentum;
}

// get the momentum along the y axis
// 
// return avlue = momentum along the y axis

inline int FermionOnTorusWithMagneticTranslations::GetKyMomentum()
{
  return this->YMomentum;
}

// get the maximum momentum along the x axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the x axis

inline int FermionOnTorusWithMagneticTranslations::GetMaxXMomentum()
{
  return this->MomentumModulo;
}
  
// get the maximum momentum along the y axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the y axis

inline int FermionOnTorusWithMagneticTranslations::GetMaxYMomentum()
{
  return this->MaxMomentum;
}
  
// find canonical form of a state description
//
// stateDescription = unsigned integer describing the state
// maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description

inline unsigned long FermionOnTorusWithMagneticTranslations::FindCanonicalForm(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;
  int index = 1;
  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
  while (stateDescriptionReference != stateDescription)
    {
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = index;
	}
      stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      ++index;
    }
  if (nbrTranslation != 0)
    {
      maxMomentum = this->MaxMomentum;
      stateDescription = 0x1ul << this->MaxMomentum;
      while ((CanonicalState & stateDescription) ==0)      
	{
	  --maxMomentum;
	  stateDescription >>= 1;
	}
      nbrTranslation = index - nbrTranslation;
    }
  return CanonicalState;
}

// find how many translations on the x direction are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

inline int FermionOnTorusWithMagneticTranslations::FindNumberXTranslation(unsigned long stateDescription)
{
  unsigned long TmpState = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
  int index = 1;  
  while (TmpState != stateDescription)
    {
      TmpState = (TmpState >> this->StateShift) | ((TmpState & this->MomentumMask) << this->ComplementaryStateShift);
      ++index;
    }
  return index;
}

// test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescription = unsigned integer describing the state
// maxMomentum = maximum momentum value that can be reached by a fermion in the stateDescription state
// return value = true if the state satisfy the x momentum constraint

inline bool FermionOnTorusWithMagneticTranslations::TestXMomentumConstraint(unsigned long stateDescription, int maxMomentum)
{
  if (this->NbrFermions & 1)
    {
      if (this->XMomentum == 0)
	return true;
      unsigned long TmpState = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      int index = 1;  
      while (TmpState != stateDescription)
	{
	  TmpState = (TmpState >> this->StateShift) | ((TmpState & this->MomentumMask) << this->ComplementaryStateShift);
	  ++index;
	}
      if (((this->XMomentum * index) % this->MomentumModulo) == 0)
	return true;
      else
	return false;
    }
  else
    {
      unsigned long TmpState2 = stateDescription & this->MomentumMask;
      unsigned long TmpState = (stateDescription >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
      int TmpSignature = 0;
      int index = 1;  
#ifndef  __64_BITS__
      TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 16) & 0xfffful]) & 1;
#else
      TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 16) & 0xfffful]
		       + this->NbrParticleLookUpTable[(TmpState2 >> 32) & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 48) & 0xfffful]) & 1;
#endif

      while (TmpState != stateDescription)
	{
	  TmpState2 = TmpState & this->MomentumMask;
	  TmpState = (TmpState >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
#ifndef  __64_BITS__
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 16) & 0xfffful]) & 1;
#else
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 16) & 0xfffful]
			   + this->NbrParticleLookUpTable[(TmpState2 >> 32) & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 48) & 0xfffful]) & 1;
#endif
	  ++index;
	}
      if ((((this->XMomentum * index) - ((this->MomentumModulo * TmpSignature) >> 1 )) % this->MomentumModulo) == 0)
	return true;
      else
	return false;
    }
}

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescription = unsigned integer describing the state
// maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description and -1 in nbrTranslation if the state does not fit the x momentum constraint

inline unsigned long FermionOnTorusWithMagneticTranslations::FindCanonicalFormAndTestXMomentumConstraint(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;
  int index = 1;  
  if (this->NbrFermions & 1)
    {
      stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      while (stateDescriptionReference != stateDescription)
	{
	  if (stateDescription < CanonicalState)
	    {
	      CanonicalState = stateDescription;
	      nbrTranslation = index;
	    }
	  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
	  ++index;
	}
      if (((this->XMomentum * index) % this->MomentumModulo) != 0)
	{
	  nbrTranslation = -1;
	  return CanonicalState;
	}
    }
  else
    {
      unsigned long TmpState = stateDescription & this->MomentumMask;
      stateDescription = (stateDescription >> this->StateShift) | (TmpState << this->ComplementaryStateShift);
      int TmpSignature = 0;
#ifndef  __64_BITS__
      TmpSignature += (this->NbrParticleLookUpTable[TmpState & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 16) & 0xfffful]) & 1;
#else
      TmpSignature += (this->NbrParticleLookUpTable[TmpState & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 16) & 0xfffful]
		       + this->NbrParticleLookUpTable[(TmpState >> 32) & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 48) & 0xfffful]) & 1;
#endif

      while (stateDescription != stateDescriptionReference)
	{
	  if (stateDescription < CanonicalState)
	    {
	      CanonicalState = stateDescription;
	      nbrTranslation = index;
	    }
	  TmpState = stateDescription & this->MomentumMask;
	  stateDescription = (stateDescription >> this->StateShift) | (TmpState << this->ComplementaryStateShift);
#ifndef  __64_BITS__
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 16) & 0xfffful]) & 1;
#else
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 16) & 0xfffful]
			   + this->NbrParticleLookUpTable[(TmpState >> 32) & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 48) & 0xfffful]) & 1;
#endif
	  ++index;
	}
      if ((((this->XMomentum * index) - ((this->MomentumModulo * TmpSignature) >> 1)) % this->MomentumModulo) != 0)
	{
	  nbrTranslation = -1;
	  return CanonicalState;
	}
    }
  if (nbrTranslation != 0)
    {
      maxMomentum = this->MaxMomentum;
      stateDescription = 0x1ul << this->MaxMomentum;
      while ((CanonicalState & stateDescription) ==0)      
	{
	  --maxMomentum;
	  stateDescription >>= 1;
	}
      nbrTranslation = index - nbrTranslation;
    }
  return CanonicalState;
}

// find state index
//
// stateDescription = unsigned longeger describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

inline int FermionOnTorusWithMagneticTranslations::FindStateIndex(unsigned long stateDescription, int maxMomentum)
{
  long PosMax = stateDescription >> this->LookUpTableShift[maxMomentum];
  long PosMin = this->LookUpTable[maxMomentum][PosMax];
  PosMax = this->LookUpTable[maxMomentum][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMin != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    return PosMax;
}

// get the C2 symmetric state of a given state 
//
// index = index of the state whose symmetric counterpart has to be computed
// nbrTranslation = number of translations that has to be applied to C2 symmetric state to be in the canonical form

inline int FermionOnTorusWithMagneticTranslations::GetC2SymmetricState (int index, int& nbrTranslation)
{
  //  cout << index << " : " << hex << this->StateDescription[index] << " ";
  unsigned long InitialState = this->StateDescription[index] << this->InvertShift;
  //  cout << hex << InitialState << " ";
#ifdef __64_BITS__
  unsigned long TmpState = FermionOnSphereInvertTable[InitialState & 0xff] << 56;
  TmpState |= FermionOnSphereInvertTable[(InitialState >> 8) & 0xff] << 48;
  TmpState |= FermionOnSphereInvertTable[(InitialState >> 16) & 0xff] << 40;
  TmpState |= FermionOnSphereInvertTable[(InitialState >> 24) & 0xff] << 32;
  TmpState |= FermionOnSphereInvertTable[(InitialState >> 32) & 0xff] << 24;
  TmpState |= FermionOnSphereInvertTable[(InitialState >> 40) & 0xff] << 16;
  TmpState |= FermionOnSphereInvertTable[(InitialState >> 48) & 0xff] << 8;
  TmpState |= FermionOnSphereInvertTable[InitialState >> 56]; 
#else
  unsigned long TmpState = FermionOnSphereInvertTable[InitialState & 0xff] << 24;
  TmpState |= FermionOnSphereInvertTable[(InitialState >> 8) & 0xff] << 16;
  TmpState |= FermionOnSphereInvertTable[(InitialState >> 16) & 0xff] << 8;
  TmpState |= FermionOnSphereInvertTable[InitialState >> 24];
#endif	
  //  cout << TmpState << " ";
  TmpState >>= this->InvertUnshift;
  TmpState |= TmpState >> this->MaxMomentum;
  TmpState &= ~(0x1ul << this->MaxMomentum);
  //  cout << TmpState  << dec;
  int TmpStateMaxMomentum = this->MaxMomentum;
  while ((TmpState >> TmpStateMaxMomentum) == 0x0ul)
    --TmpStateMaxMomentum;
  //  cout << " " << TmpStateMaxMomentum;
  TmpState = this->FindCanonicalForm(TmpState, TmpStateMaxMomentum, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, TmpStateMaxMomentum) == false)
    {
      return this->HilbertSpaceDimension;
    }
  //  cout << " " << hex << TmpState << " " << dec << TmpStateMaxMomentum << endl;;
  return this->FindStateIndex(TmpState, TmpStateMaxMomentum);
}

// convert a fermionic state to its monomial representation
//
// initialState = initial fermionic state in its fermionic representation
// finalState = reference on the array where the monomial representation has to be stored

inline void FermionOnTorusWithMagneticTranslations::ConvertToMonomial(unsigned long initialState, unsigned long*& finalState)
{
  int Index = 0;
  for (long j = this->MaxMomentum; j >= 0l; --j)
    if (((initialState >> j) & 1ul) != 0ul)
      finalState[Index++] = (unsigned long) j;
}


// convert a fermionic state from its monomial representation
//
// initialState = array where the monomial representation is stored
// return value = fermionic state in its fermionic representation

inline unsigned long FermionOnTorusWithMagneticTranslations::ConvertFromMonomial(unsigned long* initialState)
{
  unsigned long TmpState = 0x0ul;  
  for (int j = 0; j < this->NbrFermions; ++j)
    TmpState |= 0x1ul << initialState[j];
  return TmpState;
 }

#endif


