////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of bosons on torus with magnetic translations          //
//                         and for system size such that                      //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 30/01/2012                      //
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


#ifndef BOSONONTORUSWITHMAGNETICTRANSLATIONSSHORT_H
#define BOSONONTORUSWITHMAGNETICTRANSLATIONSSHORT_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Matrix/HermitianMatrix.h"

#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;


class BosonOnTorusWithMagneticTranslationsShort :  public ParticleOnTorusWithMagneticTranslations
{

 protected:

  // number of bosons
  int NbrBosons;
  // number of bosons plus 1
  int IncNbrBosons;
  // number of flux quanta
  int MaxMomentum;
  // number of flux quanta for the fermionic representation
  int FermionicMaxMomentum;
  // number of Ky values in a state
  int NbrKyValue;

  // momentum value in the x direction (modulo GCD of nbrBosons and maxMomentum)
  int KxMomentum;
  // momentum value in the y direction (modulo GCD of nbrBosons and maxMomentum)
  int KyMomentum;
  //  GCD of nbrBosons and maxMomentum
  int MomentumModulo;
  // translation step used for the magnetic translation along x 
  int XMomentumTranslationStep;

  // index of the momentum orbit
  int TotalKy;

  // array describing each state
  unsigned long* StateDescription;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given lzmax sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;

  // value that has to be substracted to the momentum for each translation of the canonical form research
  int MomentumIncrement;
  // shift that has to be done on a state for each translation of the canonical form research
  int StateShift;
  // complementary shift (with respect to MaxMomentum) to StateShift
  int ComplementaryStateShift;
  // mask that corresponds to last bit that can be set to one
  unsigned long LastMomentumMask;

  // temporary state used when applying operators
  unsigned long* TemporaryState;
  int TemporaryStateKyMax;
  // temporary state used when applying ProdA operator
  unsigned long* ProdATemporaryState;
  int ProdATemporaryStateKyMax;
  int ProdATemporaryStateNbrStateInOrbit;


   // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;
  
public:

  // default constructor
  // 
  BosonOnTorusWithMagneticTranslationsShort ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // kxMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // kyMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)
  BosonOnTorusWithMagneticTranslationsShort (int nbrBosons, int maxMomentum, int kxMomentum, int kyMomentum);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnTorusWithMagneticTranslationsShort(const BosonOnTorusWithMagneticTranslationsShort& bosons);

  // destructor
  //
  ~BosonOnTorusWithMagneticTranslationsShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnTorusWithMagneticTranslationsShort& operator = (const BosonOnTorusWithMagneticTranslationsShort& bosons);

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

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation);

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
  
   // apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next Ad call
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A (int index, int n);

  // apply a^+_m a_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for creation operator
  // return value =  resulting multiplicative factor 
  virtual double AdA (int index, int m);
  
  // convert a state to its occupation number representation
  //
  // index = index of the state
  // finalState = reference on the array where the occupation number representation has to be stored
  virtual void GetOccupationNumber(long index, unsigned long*& finalState);

  // return matrix representation of the annihilation operator a_i
  //
  // i = operator index
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  Matrix& A (int i, Matrix& M);

  // return matrix representation ofthw creation operator a^+_i
  //
  // i = operator index
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  Matrix& Ad (int i, Matrix& M);

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomial (ostream& Str, long state);
  
  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given (Kx,Ky) sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // kxSector = Kx sector in which the density matrix has to be evaluated 
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kxSector, int kySector, ComplexVector& groundState);

 protected:

  // find canonical form of a state description
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
  // nbrTranslation = number of translation needed to obtain the canonical form
  // yMomentum = state momentum value in the y direction
  // return value = canonical form of a state description
  //  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation);

  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = canonical form of a state description and -1 in nbrTranslation if the state does not fit the x momentum constraint
  //  virtual unsigned long FindCanonicalFormAndTestXMomentumConstraint(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation);

 // find how many translations on the x direction are needed to obtain the same state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of translation needed to obtain the same state
  virtual int FindNumberXTranslation(unsigned long stateDescription);

  // test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfy the x momentum constraint
  virtual bool TestXMomentumConstraint(unsigned long stateDescription);

  // find canonical form of a state description
  //
  // stateDescription = fermionic representation of the state
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = canonical form of a state description (fermionic representation)
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation);
 
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescription = fermionic representation of the state
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = canonical form of a state description (fermionic representation) and -1 in nbrTranslation if the state does not fit the x momentum constraint
  virtual unsigned long FindCanonicalFormAndTestXMomentumConstraint(unsigned long stateDescription, int& nbrTranslation);

  // apply a single translation to a bosonic state in its fermionic representation
  //
  // stateDescription = state to translate
  virtual void ApplySingleTranslation(unsigned long& stateDescription);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initialbosonic  state is stored
  // initialStateLzMax = reference on the initial bosonic state maximum Lz value
  // return value = corresponding fermionic state
  virtual unsigned long BosonToFermion(unsigned long*& initialState, int& initialStateLzMax);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax = initial fermionic state maximum Lz value
  // finalState = reference on the array where the bosonic state has to be stored
  // finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored
  virtual void FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax);

  // find state index
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKyMax = momentum maximum value for bosons that are still to be placed
  // currentMomentum = current value of the momentum
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKyMax, int currentMomentum);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(int memory);

  // generate all states with both the kx and ky constraint
  // 
  // return value = new dimension of the Hilbert space
  virtual long GenerateStates();

  // generate all states corresponding to the ky constraint  without taking care of the kx constraint
  // 
  // nbrBosons = number of bosons
  // currentKyMax = momentum maximum value for bosons that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentMomentum = current value of the momentum
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrBosons, int currentKyMax, long pos, int currentMomentum);

  // convert a bosonic state to its monomial representation
  //
  // initialState = initial  bosonic state
  // initialStateLzMax = initial bosonic state maximum Ky value
  // finalState = reference on the array where the monomial representation has to be stored
  virtual void ConvertToMonomial(unsigned long* initialState, int initialStateKyMax, unsigned long*& finalState);

  // convert a bosonic state to its monomial representation
  //
  // initialState = initial bosonic state in its fermionic representation
  // initialStateLzMax = initial bosonic state maximum Lz value
  // finalState = reference on the array where the monomial representation has to be stored
  virtual void ConvertToMonomial(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState);

  // convert a bosonic state from its monomial representation
  //
  // initialState = array where the monomial representation is stored
  // finalState = bosonic state
  // return value = maximum Ky value reached by a particle
  virtual int ConvertFromMonomial(unsigned long* initialState, unsigned long*& finalState);

  // convert a bosonic state from its monomial representation
  //
  // initialState = array where the monomial representation is stored
  // return value = bosonic state in its fermionic representation
  virtual unsigned long ConvertFromMonomial(unsigned long* initialState);

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

 protected:

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

};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnTorusWithMagneticTranslationsShort::GetParticleStatistic()
{
  return ParticleOnTorusWithMagneticTranslations::BosonicStatistic;
}

// get the number of orbitals
//
// return value = number of orbitals

inline int BosonOnTorusWithMagneticTranslationsShort::GetNbrOrbitals()
{
  return this->MaxMomentum;
}

// get the number of particles
//
// return value = number of particles

inline int BosonOnTorusWithMagneticTranslationsShort::GetNbrParticles()
{
  return this->NbrBosons;
}

// get the momentum along the x axis
// 
// return avlue = momentum along the x axis

inline int BosonOnTorusWithMagneticTranslationsShort::GetKxMomentum()
{
  return this->KxMomentum;
}

// get the momentum along the y axis
// 
// return avlue = momentum along the y axis

inline int BosonOnTorusWithMagneticTranslationsShort::GetKyMomentum()
{
  return this->KyMomentum;
}

// get the maximum momentum along the x axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the x axis

inline int BosonOnTorusWithMagneticTranslationsShort::GetMaxXMomentum()
{
  return this->MomentumModulo;
}
  
// get the maximum momentum along the y axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the y axis

inline int BosonOnTorusWithMagneticTranslationsShort::GetMaxYMomentum()
{
  return this->MaxMomentum;
}
  
// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initialbosonic  state is stored
// initialStateLzMax = reference on the initial bosonic state maximum Lz value
// return value = corresponding fermionic state

inline unsigned long BosonOnTorusWithMagneticTranslationsShort::BosonToFermion(unsigned long*& initialState, int& initialStateLzMax)
{
  unsigned long TmpState = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= initialStateLzMax; ++i)
    {
      TmpState |= ((1ul << initialState[i]) - 1ul) << Shift;
      Shift += initialState[i];
      ++Shift;
    }
  return TmpState;
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateLzMax = initial fermionic state maximum Lz value
// finalState = reference on the array where the bosonic state has to be stored
// finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored

inline void BosonOnTorusWithMagneticTranslationsShort::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax)
{
  finalStateLzMax = 0;
  while (initialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialState - 1ul) ^ (~initialState);
      TmpState &= ~(TmpState >> 1);
#ifdef  __64_BITS__
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0f0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ffff0000ul) != 0) << 4;      
      TmpPower |= ((TmpState & 0xffffffff00000000ul) != 0) << 5;      
#else
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ul) != 0) << 4;      
#endif
      finalState[finalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState >>= TmpPower;
      ++finalStateLzMax;
      initialStateLzMax -= TmpPower;
    }
  --finalStateLzMax;
}

// convert a bosonic state to its monomial representation
//
// initialState = initial  bosonic state
// initialStateKyMax = initial bosonic state maximum Ky value
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnTorusWithMagneticTranslationsShort::ConvertToMonomial(unsigned long* initialState, int initialStateKyMax, unsigned long*& finalState)
{
  int Index = 0;
  for (int i = initialStateKyMax; i >= 0; --i)
    for (unsigned long j = 0ul; j < initialState[i]; ++j)
      finalState[Index++] = i;
}

// convert a bosonic state to its monomial representation
//
// initialState = initial bosonic state in its fermionic representation
// initialStateKyMax = initial bosonic state maximum Ky value
// finalState = reference on the array where the monomial representation has to be stored

inline void BosonOnTorusWithMagneticTranslationsShort::ConvertToMonomial(unsigned long initialState, int initialStateKyMax, unsigned long*& finalState)
{
  int Index = 0;
  int TmpKy = initialStateKyMax - this->NbrBosons + 1;
  while (initialStateKyMax >= 0)
    {
      while ((initialStateKyMax >= 0) && (((initialState >> initialStateKyMax) & 0x1ul) != 0x0ul))
	{
	  finalState[Index++] = TmpKy;
	  --initialStateKyMax;
	}
      while ((initialStateKyMax >= 0) && (((initialState >> initialStateKyMax) & 0x1ul) == 0x0ul))
	{
	  --TmpKy;
	  --initialStateKyMax;
	}
    }
}

// convert a state to its occupation number representation
//
// index = index of the state
// finalState = reference on the array where the occupation number representation has to be stored

inline void BosonOnTorusWithMagneticTranslationsShort::GetOccupationNumber(long index, unsigned long*& finalState)
{
  int FinalStateLzMax;
  unsigned long InitialState = this->StateDescription[index];
  int InitialStateLzMax = this->FermionicMaxMomentum;
  while ((InitialState >> InitialStateLzMax) == 0x0ul)
    {
      --InitialStateLzMax;
    }
  this->FermionToBoson(InitialState, InitialStateLzMax, finalState, FinalStateLzMax);
  for (++FinalStateLzMax; FinalStateLzMax < this->MaxMomentum; ++FinalStateLzMax)
    finalState[FinalStateLzMax] = 0x0ul;
}

// find how many translations on the x direction are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

inline int BosonOnTorusWithMagneticTranslationsShort::FindNumberXTranslation(unsigned long stateDescription)
{
  unsigned long TmpState = stateDescription;
  this->ApplySingleTranslation(TmpState);
  int Index = 1;  
  while (TmpState != stateDescription)
    {
      this->ApplySingleTranslation(TmpState);
      ++Index;  
    }
  return Index;
}

// find canonical form of a state description
//
// stateDescription = fermionic representation of the state
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description (fermionic representation)

inline unsigned long BosonOnTorusWithMagneticTranslationsShort::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  for (int i = 0; i < this->MomentumModulo; ++i)
    {
      this->ApplySingleTranslation(stateDescription);
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = 1;
	}
    }
  return CanonicalState;
}

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescription = fermionic representation of the state
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description (fermionic representation) and -1 in nbrTranslation if the state does not fit the x momentum constraint

inline unsigned long BosonOnTorusWithMagneticTranslationsShort::FindCanonicalFormAndTestXMomentumConstraint(unsigned long stateDescription, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  unsigned long InitialStateDescription = stateDescription;
  int Index = 0;
  int OrbitSize = 0;
  while (Index < this->MomentumModulo)
    {
      this->ApplySingleTranslation(stateDescription);
      ++Index;
      ++OrbitSize;
      if (stateDescription != InitialStateDescription)
	{
	  if (stateDescription < CanonicalState)
	    {
	      CanonicalState = stateDescription;
	      nbrTranslation = Index;
	    }
	}
      else
	{
	  Index = this->MomentumModulo;
	}
    }

  if (((this->KxMomentum * OrbitSize) % this->MomentumModulo) != 0)
    {
      nbrTranslation = -1;
    }

  return CanonicalState;
}


// apply a single translation to a bosonic state in its fermionic representation
//
// stateDescription = state to translate

inline void BosonOnTorusWithMagneticTranslationsShort::ApplySingleTranslation(unsigned long& stateDescription)
{
  for (int i = 0; i < this->StateShift;)
    {
      while ((i < this->StateShift) && ((stateDescription & 0x1ul) == 0x0ul))
	{
	  stateDescription >>= 1;
	  ++i;
	}
      if (i < this->StateShift)
	{
	  while ((stateDescription & 0x1ul) == 0x1ul)
	    {
	      stateDescription >>= 1;
	      stateDescription |= this->LastMomentumMask;
	    }
	  stateDescription >>= 1;	  
	  ++i;
	}
    }
}

// test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfy the x momentum constraint

inline bool BosonOnTorusWithMagneticTranslationsShort::TestXMomentumConstraint(unsigned long stateDescription)
{
  if (((this->KxMomentum * this->FindNumberXTranslation(stateDescription)) % this->MomentumModulo) == 0)
    return true;
  else
    return false;
}

// convert a bosonic state from its monomial representation
//
// initialState = array where the monomial representation is stored
// finalState = bosonic state
// return value = maximum Ky value reached by a particle

inline int BosonOnTorusWithMagneticTranslationsShort::ConvertFromMonomial(unsigned long* initialState, unsigned long*& finalState)
{
  int TmpKyMax = initialState[0]; 
  for (int i = 0; i < TmpKyMax; ++i)
    finalState[i] = 0ul;
  for (int i = 0; i < this->NbrBosons; ++i)
    finalState[initialState[i]]++;
  return TmpKyMax;
}

// convert a bosonic state from its monomial representation
//
// initialState = array where the monomial representation is stored
// return value = bosonic state in its fermionic representation

inline unsigned long BosonOnTorusWithMagneticTranslationsShort::ConvertFromMonomial(unsigned long* initialState)
{
  unsigned long Tmp = 0x0ul;
  for (int i = 0; i < this->NbrBosons; ++i)
    Tmp |= 0x1ul << (initialState[i] + ((unsigned long) (this->NbrBosons - i)) - 1ul);
  return Tmp;
}

#endif


