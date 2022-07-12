////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with SU(8) spin              //
//                                                                            //
//                        last modification : 11/05/2020                      //
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


#ifndef FERMIONONSPHEREWITHSU8SPIN_H
#define FERMIONONSPHEREWITHSU8SPIN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU8Spin.h"

#include <iostream>

#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

class FermionOnSphere;
class FermionOnSphereWithSpin;
class ParticleOnSphereWithSpin;


class FermionOnSphereWithSU8Spin :  public ParticleOnSphereWithSU8Spin
{

 protected:

  // number of fermions
  int NbrFermions;
  // number of fermions plus 1
  int IncNbrFermions;
  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a fermion
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;

  // highest bit in a given state description
  int HighestBit;

  // number of particles with sigma=1
  int NbrFermions1;
  // number of particles with sigma=2
  int NbrFermions2;
  // number of particles with sigma=3
  int NbrFermions3;
  // number of particles with sigma=4
  int NbrFermions4;
  // number of particles with sigma=5
  int NbrFermions5;
  // number of particles with sigma=6
  int NbrFermions6;
  // number of particles with sigma=7
  int NbrFermions7;
  // number of particles with sigma=8
  int NbrFermions8;

  // array describing each state
  unsigned long* StateDescription;
  // array giving maximum Lz value reached for a fermion in a given state
  int* StateHighestBit;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given lzmax sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;

  // a table containing ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  unsigned long* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;

  // temporary state used when applying ProdA operator
  unsigned long ProdATemporaryState;
  // Lz maximum value associated to temporary state used when applying ProdA operator
  int ProdALzMax;

 public:

  // default constructor
  //
  FermionOnSphereWithSU8Spin();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // nbrParticleSigma = array that provides the number of particles with a given internal degree of freedom
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSU8Spin (int nbrFermions, int totalLz, int lzMax, int* nbrParticleSigma, unsigned long memory = 10000000);

   // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSU8Spin(const FermionOnSphereWithSU8Spin& fermions);

  // destructor
  //
  ~FermionOnSphereWithSU8Spin ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSU8Spin& operator = (const FermionOnSphereWithSU8Spin& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

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

  // apply a^+_m_s a_m_s operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // sigma = internal degree of freedom label of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdsigmaAsigma (int index, int m, int sigma);

  // apply a^+_m_s a_m_s operator to a given state)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // sigma = internal degree of freedom label of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdsigmaAsigma (long index, int m, int sigma);

  // apply a^+_m1_s1 a_m2_s2 operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m1 = index of the creation operator
  // sigma1 = internal degree of freedom label of the creation operator
  // m2 = index of the annihilation operator
  // sigma2 = internal degree of freedom label of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdsigmaAsigma (int index, int m1, int sigma1, int m2, int sigma2, double& coefficient);

  // apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up, 1 for um, 2 for dp and 3 for dm 
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // sigma1 = SU(4) index for the first annihilation operator
  // sigma2 = SU(4) index for the second annihilation operator
  // return value =  multiplicative factor 
  virtual double AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2);

  // apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is 0 for up, 1 for um, 2 for dp and 3 for dm 
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // sigma1 = SU(4) index for the first creation operator
  // sigma2 = SU(4) index for the second creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // evaluate wave function in real space using a given basis and only for agiven range of components
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
				int firstComponent, int nbrComponent);                                
  
  // initialize evaluation of wave function in real space using a given basis and only for a given range of components and
  //
  // timeCoherence = true if time coherence has to be used
  void InitializeWaveFunctionEvaluation (bool timeCoherence = false);
  
  // create a U(1) state from an SU(4) state
  //
  // state = vector describing the SU(4) state
  // u1Space = reference on the Hilbert space associated to the U(1) state
  // return value = resulting U(1) state
  virtual RealVector ForgeU1FromSU8(RealVector& state, FermionOnSphere& u1Space);

  // create a SU(2) state from an SU(4) state (fusing same spin values,i.e symmetrizing over the isospin)
  //
  // state = vector describing the SU(4) state
  // su2Space = reference on the Hilbert space associated to the SU(2) state
  // return value = resulting SU(2) state
  virtual RealVector ForgeSU2FromSU8(RealVector& state, FermionOnSphereWithSpin& su2Space);

  // convert a state from one SU(4) basis to another, transforming the one body basis in each momentum sector
  //
  // initialState = state to transform  
  // targetState = vector where the transformed state has to be stored
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // firstComponent = index of the first component to compute in initialState
  // nbrComponents = number of consecutive components to compute
  virtual void TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, long firstComponent = 0l, long nbrComponents = 0l);

  // compute the transformation matrix from one SU(4) basis to another, transforming the one body basis in each momentum sector
  //
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // return value = transformation matrix
  virtual ComplexMatrix TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis);

  // compute the projection matrix from the SU(8) Hilbert space to an SU(2) Hilbert space
  // 
  // targetSpace = pointer to the SU(2) Hilbert space
  // spinUp = index of the component that has to be consider as a spin up
  // spinDown = index of the component that has to be consider as a spin down
  // return value = projection matrix
  virtual ComplexMatrix TransformationMatrixSU8ToSU2(ParticleOnSphereWithSpin* targetSpace, int spinUp = 0, int spinDown = 1);

  protected:

  // factorized code for any a^+_m_x a_n_y operator 
  //
  // index = index of the state on which the operator has to be applied
  // m = global index of the creation operator
  // n = global index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int GenericAdA(int index, int m, int n, double& coefficient);

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription, int lzmax);


  // evaluate Hilbert space dimension with shifted values for lzMax and totalLz
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
  // nbrParticles1 = number of particles with sigma=1
  // nbrParticles2 = number of particles with sigma=2
  // nbrParticles3 = number of particles with sigma=3
  // nbrParticles4 = number of particles with sigma=4
  // nbrParticles5 = number of particles with sigma=5
  // nbrParticles6 = number of particles with sigma=6
  // nbrParticles7 = number of particles with sigma=7
  // nbrParticles8 = number of particles with sigma=8
  // return value = Hilbert space dimension  
  long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int nbrParticles1, int nbrParticles2,
					    int nbrParticles3, int nbrParticles4, int nbrParticles5, int nbrParticles6,
					    int nbrParticles7, int nbrParticles8);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // nbrParticles1 = number of particles with sigma=1
  // nbrParticles2 = number of particles with sigma=2
  // nbrParticles3 = number of particles with sigma=3
  // nbrParticles4 = number of particles with sigma=4
  // nbrParticles5 = number of particles with sigma=5
  // nbrParticles6 = number of particles with sigma=6
  // nbrParticles7 = number of particles with sigma=7
  // nbrParticles8 = number of particles with sigma=8
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrFermions, int lzMax, int totalLz, int nbrParticles1, int nbrParticles2,
		      int nbrParticles3, int nbrParticles4, int nbrParticles5, int nbrParticles6,
		      int nbrParticles7, int nbrParticles8, long pos);

  // recursive part of the convertion from a state from one SU(4) basis to another, transforming the one body basis in each momentum sector
  //
  // targetState = vector where the transformed state has to be stored
  // coefficient = current coefficient to assign
  // position = current particle consider in the n-body state
  // momentumIndices = array that gives the momentum partition of the initial n-body state
  // initialSU8Indices = array that gives the spin dressing the initial n-body state
  // currentSU8Indices = array that gives the spin dressing the current transformed n-body state
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  void TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
				      int position, int* momentumIndices, int* initialSU8Indices, int* currentSU8Indices, ComplexMatrix* oneBodyBasis);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnSphereWithSU8Spin::GetParticleStatistic()
{
  return AbstractQHEParticle::FermionicStatistic;
}

// apply a^+_m1_s1 a_m2_s2 operator to a given state
//
// index = index of the state on which the operator has to be applied
// m1 = index of the creation operator
// sigma1 = internal degree of freedom label of the creation operator
// m2 = index of the annihilation operator
// sigma2 = internal degree of freedom label of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int FermionOnSphereWithSU8Spin::AdsigmaAsigma (int index, int m1, int sigma1, int m2, int sigma2, double& coefficient)
{
  return this->GenericAdA(index, (m1 << 3) + sigma1, (m2 << 3) + sigma2, coefficient);
}

// factorized code for any a^+_m_x a_n_y operator 
//
// index = index of the state on which the operator has to be applied
// m = global index of the creation operator
// n = global index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int FermionOnSphereWithSU8Spin::GenericAdA(int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0x0ul) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = -this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(0x1ul << n);
  if (State != 0x0ul)
    while ((State >> NewLargestBit) == 0x0ul)
      --NewLargestBit;

  if ((State & (0x1ul << m)) != 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= 0x1ul << m;
  return this->FindStateIndex(State, NewLargestBit);
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call.
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(8) index for the first annihilation operator
// sigma2 = SU(8) index for the second annihilation operator
// return value =  multiplicative factor 

inline double FermionOnSphereWithSU8Spin::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 3;
  n1 +=  sigma1;
  n2 <<= 3;
  n2 +=  sigma2;
  if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  if (this->ProdATemporaryState != 0x0ul)
    {
      while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
	--this->ProdALzMax;
    }
  else
    this->ProdALzMax = 0;
  return Coefficient;
}

// apply a^+_m_s a_m_s operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// sigma = internal degree of freedom label of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

inline double FermionOnSphereWithSU8Spin::AdsigmaAsigma (int index, int m, int sigma)
{
  m <<= 3;
  m += sigma;
  return ((double) ((this->StateDescription[index] >> m) & 0x1ul));
}

// apply a^+_m_s a_m_s operator to a given state)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// sigma = internal degree of freedom label of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

inline double FermionOnSphereWithSU8Spin::AdsigmaAsigma (long index, int m, int sigma)
{
  m <<= 3;
  m += sigma;
  return ((double) ((this->StateDescription[index] >> m) & 0x1ul));
}

// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(8) index for the first creation operator
// sigma2 = SU(8) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int FermionOnSphereWithSU8Spin::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 3;
  m1 +=  sigma1;
  m2 <<= 3;
  m2 +=  sigma2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

#endif


