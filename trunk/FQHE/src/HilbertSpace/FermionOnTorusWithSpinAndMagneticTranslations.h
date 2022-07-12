////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of fermion with spin on a torus taking               //
//                      into account magnetic translations                    //
//                                                                            //
//                        last modification : 26/11/2007                      //
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


#ifndef FERMIONONTORUSWITHSPINANDMAGNETICTRANSLATIONS_H
#define FERMIONONTORUSWITHSPINANDMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"

using std::cout;
using std::endl;
using std::dec;
using std::hex;


#define FERMION_TORUS_SU2_SZ_MASK 0x5555555555555555ul


class FermionOnTorusWithSpinAndMagneticTranslations :  public ParticleOnTorusWithSpinAndMagneticTranslations
{


 protected:

  // number of fermions
  int NbrFermions;

  // number of fermions with spin up
  int NbrFermionsUp;

  // number of fermions with spin down
  int NbrFermionsDown;
  
  // total number of fermions plus 1
  int IncNbrFermions;

  // total value of Spin
  int TotalSpin;

  // maximum momentum value reached by a fermion
  int MaxMomentum;
  // number of momentum values in a state (= MaxMomentum +1)
  int NbrMomentum;
  // number of different electron orbitals in a state (= 2*(MaxMomentum +1))
  int NbrFermionStates;
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

  // array describing each state 
  unsigned long* StateDescription;
  // array giving maximum momentum value reached for a fermion in a given state
  int* StateHighestBit;
  
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

  // temporary state used when applying ProdA operator
  unsigned long ProdATemporaryState;
  // Highest Bit value associated to temporary state used when applying ProdA operator
  int ProdAHighestBit;
  // Index of the initial state when applying ProdA operator
  int ProdAIndex;
  //  number of states in the orbit of the  temporary state used when applying ProdA operator
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
  FermionOnTorusWithSpinAndMagneticTranslations ();  

  // basic constructor
  // 
  // nbrFermions= number of fermions
  // totalSpin = twice the total spin value
  // maxMomentum = momentum maximum value for a fermion
  // xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)  
  FermionOnTorusWithSpinAndMagneticTranslations (int nbrFermions, int totalSpin, int maxMomentum, int xMomentum, int yMomentum);

  // basic constructor without constraint on Sz
  // 
  // nbrFermions= number of fermions
  // maxMomentum = momentum maximum value for a fermion
  // xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)  
  FermionOnTorusWithSpinAndMagneticTranslations (int nbrFermions, int maxMomentum, int xMomentum, int yMomentum);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnTorusWithSpinAndMagneticTranslations(const FermionOnTorusWithSpinAndMagneticTranslations& fermions);

  // destructor
  //
  ~FermionOnTorusWithSpinAndMagneticTranslations();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnTorusWithSpinAndMagneticTranslations& operator = (const FermionOnTorusWithSpinAndMagneticTranslations& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // get the number of orbitals
  //
  // return value = number of orbitals
  virtual int GetNbrOrbitals();

  // get the number of particles
  //
  // return value = number of particles
  virtual int GetNbrParticles();

  // get the total spin
  //
  //return value: total spin of the Hilbert space
  virtual int GetTotalSpin();

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
  
  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

  // get momemtum value in the y direction of a given state
  //
  // index = state index
  // return value = state momentum in the y direction
  int GetYMomentumValue(int index);

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

  // apply a^+_(d,m1) a^+_(d,m2) a_(d,n1) a_(d,n2) operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation);

  // apply a^+_(u,m1) a^+_(u,m2) a_(u,n1) a_(u,n2) operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation);
  // verbose version
  virtual int AduAduAuAuV (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation);

  // apply a^+_(u,m1) a^+_(d,m2) a_(d,n1) a_(u,n2) operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to be applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation);

  
  // apply a^+_m_d a_m_d operator to a given state (only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AddAd (int index, int m);

  // apply a^+_m_u a_m_u operator to a given state  (only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AduAu (int index, int m);

  // apply a^+_m_u a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AduAu (int index, int m, int n, double& coefficient, int& nbrTranslation);
  
  // apply a^+_m_d a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AddAd (int index, int m, int n, double& coefficient, int& nbrTranslation);
  
  // apply a^+_m_u a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation/annihilation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AduAd (int index, int m, int n, double& coefficient, int& nbrTranslation);

  // apply a^+_m_d a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation/annihilation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, int n, double& coefficient, int& nbrTranslation);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // return value =  multiplicative factor 
  virtual double AuAu (int index, int n1, int n2);
  virtual double AuAuV (int index, int n1, int n2);

  // apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AdAd (int index, int n1, int n2);

  // apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AuAd (int index, int n1, int n2);

  // apply a_n1_d a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AdAu (int index, int n1, int n2);

  
  // apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdu (int m1, int m2, double& coefficient, int& nbrTranslation);
  virtual int AduAduV (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAdd (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdd (int m1, int m2, double& coefficient, int& nbrTranslation);  
  
  // apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator (spin up)
  // return value =  multiplicative factor 
  virtual double Au (int index, int n);
  
  // apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n = index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double Ad (int index, int n);
 
  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // convert a state defined in the Ky basis into a state in the (Kx,Ky) basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertToKxKyBasis(ComplexVector& state, ParticleOnSphere* space);

  // convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space);

  // convert a given state from a generic basis to the current Sz subspace basis
  //
  // state = reference on the vector to convert
  // basis = reference on the basis associated to state
  // return value = converted vector
  virtual ComplexVector ConvertFromNbodyBasis(ComplexVector& state, FermionOnTorusWithSpinAndMagneticTranslations& basis);

 protected:

  // find canonical form of a state description
  //
  // stateDescription = unsigned integer describing the state
  // stateHighestBit = reference on the highest non-zero bit in the canonical representation
  // nbrTranslation = number of translation needed to obtain the canonical form
  // yMomentum = state momentum value in the y direction
  // return value = canonical form of a state description
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& stateHighestBit, int& nbrTranslation);

  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = canonical form of a state description and -1 in nbrTranslation if the state does not fit the x momentum constraint
  virtual unsigned long FindCanonicalFormAndTestXMomentumConstraint(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation);

 // find how many translations on the x direction are needed to obtain the same state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of translation needed to obtain the same state
  virtual int FindNumberXTranslation(unsigned long stateDescription);

  // test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = maximum momentum value that can be reached by a fermion in the stateDescription state
  // return value = true if the state satisfy the x momentum constraint
  virtual bool TestXMomentumConstraint(unsigned long stateDescription, int maxMomentum);

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int maxMomentum);

  // evaluate Hilbert space dimension for a given total momentum
  //
  // nbrFermions = number of fermions
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKy, int currentTotalKy);

  // evaluate Hilbert space dimension using recursive algorithm
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalMomentum = momentum total value
  // totalSpin = number of particles with spin up
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalMomentum, int totalSpin);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(int memory);

  // generate look-up table associated to sign calculations
  // 
  virtual void GenerateSignLookUpTable();

  // generate all states corresponding to the constraints
  // 
  // fullSzFlag = if true, does not apply the Sz contraint
  // fullKyFlag = if true, does not apply the Ky contraint
  // return value = hilbert space dimension
  virtual int GenerateStates(bool fullSzFlag = false, bool fullKyFlag = false);
 
  // evaluate Hilbert space dimension without the translation symmmetry along x, Sz or Ky quantu numbers
  //
  // nbrFermions = number of fermions
  // currentMomentum = current one-body momentum 
  // pos = position in StateDescription array where to store states
  // return value = Hilbert space dimension
  virtual long RawGenerateStates(int nbrFermions, int currentMomentum, long pos);

  // evaluate Hilbert space dimension without the translation symmmetry along x
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalMomentum = momentum total value
  // totalSpinUp = number of particles with spin up
  // pos = position in StateDescription array where to store states
  // return value = Hilbert space dimension
  virtual long RawGenerateStates(int nbrFermions, int lzMax, int totalMomentum, int totalSpinUp, long pos);

  // generate all states corresponding to the constraints without the translation symmmetry along x
  // 
  // nbrFermions = number of fermions
  // maxMomentum = momentum maximum value for a fermion in the state
  // currentMaxMomentum = momentum maximum value for fermions that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentMomentum = current value of the momentum
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, long pos, int currentMomentum);
  
  // generate all states corresponding to the constraints without the translation symmmetry along x
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // nbrSpinUp = number of fermions with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos);
  
  // compute the parity of the number of spin singlets 
  //
  // state = reference on the state whose parity has to be evaluated
  // return value = parity of the number of spin singlets (0 for even, 1 for odd)
  virtual unsigned long GetStateSingletParity(unsigned long state);

  // apply the spin flip symmetry on a given state
  //
  // stateDescription= state on which the spin flip symmetry has to be applied
  // return value = resulting state
  virtual unsigned long ApplySpinFlipSymmetry(unsigned long stateDescription);

  // apply a Gutzwiller projection (in the orbital space) to a given state
  //
  // state = reference on the state to project
  // space = pointer to the Hilbert space where state is defined
  // return value = Gutzwiller projected state
  virtual ComplexVector GutzwillerProjection(ComplexVector& state, ParticleOnSphere* space);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnTorusWithSpinAndMagneticTranslations::GetParticleStatistic()
{
  return ParticleOnTorusWithSpinAndMagneticTranslations::FermionicStatistic;
}

// compute the parity of the number of spin singlets 
//
// state = reference on the state whose parity has to be evaluated
// return value = parity of the number of spin singlets (0 for even, 1 for odd)

inline unsigned long FermionOnTorusWithSpinAndMagneticTranslations::GetStateSingletParity(unsigned long state)
{
  state &= (state >> 1);
  state &= FERMION_TORUS_SU2_SZ_MASK;
#ifdef __64_BITS__
  state ^= (state >> 32);
#endif
  state ^= (state >> 16);
  state ^= (state >> 8);
  state ^= (state >> 4);
  state ^= (state >> 2);
  return (state & 0x1ul);
}

// apply the spin flip symmetry on a given state
//
// stateDescription= state on which the spin flip symmetry has to be applied
// return value = resulting state

inline unsigned long FermionOnTorusWithSpinAndMagneticTranslations::ApplySpinFlipSymmetry(unsigned long stateDescription)
{
  unsigned long TmpState = ((stateDescription >> 1) ^ stateDescription) & FERMION_TORUS_SU2_SZ_MASK;
  TmpState |= TmpState << 1;
  TmpState ^= stateDescription; 
  return TmpState;
}

// get the number of orbitals
//
// return value = number of orbitals

inline int FermionOnTorusWithSpinAndMagneticTranslations::GetNbrOrbitals()
{
  return this->MaxMomentum;
}

// get the number of particles
//
// return value = number of particles

inline int FermionOnTorusWithSpinAndMagneticTranslations::GetNbrParticles()
{
  return this->NbrFermions;
}

// get the total spin
//
//return value: total spin of the Hilbert space

inline int FermionOnTorusWithSpinAndMagneticTranslations::GetTotalSpin()
{
  return this->TotalSpin;
}

// get the momentum along the x axis
// 
// return avlue = momentum along the x axis

inline int FermionOnTorusWithSpinAndMagneticTranslations::GetKxMomentum()
{
  return this->XMomentum;
}

// get the momentum along the y axis
// 
// return avlue = momentum along the y axis

inline int FermionOnTorusWithSpinAndMagneticTranslations::GetKyMomentum()
{
  return this->YMomentum;
}

// get the maximum momentum along the x axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the x axis

inline int FermionOnTorusWithSpinAndMagneticTranslations::GetMaxXMomentum()
{
  return this->MomentumModulo;
}

// get the maximum momentum along the y axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the y axis

inline int FermionOnTorusWithSpinAndMagneticTranslations::GetMaxYMomentum()
{
  return this->MaxMomentum;
}
  
#endif // FERMIONONTORUSWITHSPINANDMAGNETICTRANSLATIONS_H


