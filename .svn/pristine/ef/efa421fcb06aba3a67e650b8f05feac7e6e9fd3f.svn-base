////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                    class of bosons on sphere with SU(3) spin               //
//                                                                            //
//                        last modification : 12/10/2011                      //
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


#ifndef BOSONONSPHEREWITHSU3SPIN_H
#define BOSONONSPHEREWITHSU3SPIN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU3Spin.h"

#include <iostream>


class BosonOnSphereShort;


class BosonOnSphereWithSU3Spin :  public ParticleOnSphereWithSU3Spin
{

 protected:

  // number of bosons
  int NbrBosons;
  // number of bosons plus 1
  int IncNbrBosons;
  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a boson
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;
  // twice the total Tz value
  int TotalTz;
  // three time the total Y value
  int TotalY;

  // lzmax value for the fermionic states associated to the type 1 particles
  int N1LzMax;
  // lzmax value for the fermionic states associated to the type 1 particles
  int N2LzMax;
  // lzmax value for the fermionic states associated to the type 1 particles
  int N3LzMax;
  // largest lzmax value for the fermionic states among N1LzMax, N2LzMax and N3LzMax
  int FermionicLzMax;

  // array that contains the state description, the first entry being StateDescription1, the second entry being StateDescription2, the third entry StateDescription3
  unsigned long* StateDescriptionSigma[3];

  // temporay array describing the type 1 particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescription1;
  // temporay array describing the type 2 particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescription2;
  // temporay array describing the type 3 particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescription3;

  // sorted array that contains each unique configuration for the type 1 particles
  unsigned long* UniqueStateDescription1;
  // number of time each unique configuration for the type 1 particles appears in StateDescription1
  int* UniqueStateDescriptionSubArraySize1;
  // number of unique configurations for the type 1 particles
  long NbrUniqueStateDescription1;
  // number of unique configurations for the type 2 particles per unique type 1 particle configuration 
  int* NbrUniqueStateDescription2;
  // unique configurations for the type 2 particles per unique type 1 particle configuration 
  unsigned long** UniqueStateDescription2;
  // number of time each unique configuration for the type 2 particles appears in StateDescription2 having the same type 1 particle configucation
  int** UniqueStateDescriptionSubArraySize2;
  // first time a unique combination of type 1 and type 2 particle configurations appears in the Hilbert space
  int** FirstIndexUniqueStateDescription2;

  // temporary state used when applying operators for type 1 particles
  unsigned long* TemporaryState1;
  // temporary state used when applying operators for type 2 particles
  unsigned long* TemporaryState2;
  // temporary state used when applying operators for type 3 particles
  unsigned long* TemporaryState3;
  // array that contains the temporary state, the first entry being TemporaryState1, the second entry being TemporaryState2 and the third entry TemporaryState3
  unsigned long* TemporaryStateSigma[3];

  // temporary state used when applying ProdA operator for type 1 particles
  unsigned long* ProdATemporaryState1;
  // temporary state used when applying ProdA operator for type 2 particles
  unsigned long* ProdATemporaryState2;
  // temporary state used when applying ProdA operator for type 3 particles
  unsigned long* ProdATemporaryState3;
  // array that contains the temporary state used when applying ProdA operator, the first entry being ProdATemporaryState1, the second entry being ProdATemporaryState2 and the third entry ProdATemporaryState3
  unsigned long* ProdATemporaryStateSigma[3];

 public:

  // default constructor
  // 
  BosonOnSphereWithSU3Spin ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // totalSpin = twice the total spin value
  // totalIsospin = twice the total isospin value
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU3Spin (int nbrBosons, int totalLz, int lzMax, int totalSpin, int totalIsospin, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereWithSU3Spin(const BosonOnSphereWithSU3Spin& bosons);

  // destructor
  //
  virtual ~BosonOnSphereWithSU3Spin ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereWithSU3Spin& operator = (const BosonOnSphereWithSU3Spin& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter);

  // apply a^+_m_1 a_m_1 operator to a given state (only state 1 Tz=+1/2, Y=+1/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_1 a_m_1
  virtual double Ad1A1 (int index, int m);

  // apply a^+_m_2 a_m_2 operator to a given state (only state 2 Tz=-1/2, Y=+1/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_2 a_m_2
  virtual double Ad2A2 (int index, int m);

  // apply a^+_m_3 a_m_3 operator to a given state (only state 3 Tz=0, Y=-2/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_3 a_m_3
  virtual double Ad3A3 (int index, int m);

  // apply a^+_m_1 a_n_1 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1A1 (int index, int m, int n, double& coefficient);

  // apply a^+_m_1 a_n_2 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1A2 (int index, int m, int n, double& coefficient);

  // apply a^+_m_1 a_n_3 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1A3 (int index, int m, int n, double& coefficient);

  // apply a^+_m_2 a_n_1 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2A1 (int index, int m, int n, double& coefficient);

  // apply a^+_m_2 a_n_2 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2A2 (int index, int m, int n, double& coefficient);

  // apply a^+_m_2 a_n_3 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2A3 (int index, int m, int n, double& coefficient);

  // apply a^+_m_3 a_n_1 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad3A1 (int index, int m, int n, double& coefficient);

  // apply a^+_m_3 a_n_2 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad3A2 (int index, int m, int n, double& coefficient);

  // apply a^+_m_3 a_n_3 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad3A3 (int index, int m, int n, double& coefficient);

  // apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0, 1 or 2 
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // sigma1 = SU(3) index for the first annihilation operator
  // sigma2 = SU(3) index for the second annihilation operator
  // return value =  multiplicative factor 
  virtual double AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2);

  // apply a_n1_1 a_n2_1 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A1A1 (int index, int n1, int n2);

  // apply a_n1_1 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A1A2 (int index, int n1, int n2);

  // apply a_n1_1 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A1A3 (int index, int n1, int n2);

  // apply a_n1_2 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A2A2 (int index, int n1, int n2);

  // apply a_n1_2 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A2A3 (int index, int n1, int n2);

  // apply a_n1_3 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A3A3 (int index, int n1, int n2);

  // apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is 0, 1 or 2
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // sigma1 = SU(3) index for the first creation operator
  // sigma2 = SU(3) index for the second creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient);

  // apply a^+_m1_1 a^+_m2_1 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1Ad1 (int m1, int m2, double& coefficient);

  // apply a^+_m1_1 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1Ad2 (int m1, int m2, double& coefficient);

  // apply a^+_m1_1 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1Ad3 (int m1, int m2, double& coefficient);

  // apply a^+_m1_2 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2Ad2 (int m1, int m2, double& coefficient);

  // apply a^+_m1_2 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2Ad3 (int m1, int m2, double& coefficient);

  // apply a^+_m1_3 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad3Ad3 (int m1, int m2, double& coefficient);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // convert a state from one SU(3) basis to another, transforming the one body basis in each momentum sector
  //
  // initialState = state to transform  
  // targetState = vector where the transformed state has to be stored
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // firstComponent = index of the first component to compute in initialState
  // nbrComponents = number of consecutive components to compute
  virtual void TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, long firstComponent = 0l, long nbrComponents = 0l);

  // compute the transformation matrix from one SU(3) basis to another, transforming the one body basis in each momentum sector
  //
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // return value = transformation matrix
  virtual ComplexMatrix TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis);

  // compute the projection matrix from the SU(3) Hilbert space to an U(1) Hilbert space
  // 
  // targetSpace = pointer to the U(1) Hilbert space
  // type = type of particles that has to be kept (0 for type 1, 1 for type 2, 2 for type 3
  // return value = projection matrix
  virtual ComplexMatrix TransformationMatrixSU3ToU1(BosonOnSphereShort* targetSpace, int type = 0);

  protected:

  // find state index
  //
  // stateDescription1 = unsigned integer describing the fermionic state for type 1 particles
  // stateDescription2 = unsigned integer describing the fermionic state for type 2 particles
  // stateDescription3 = unsigned integer describing the fermionic state for type 3 particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription1, unsigned long stateDescription2, unsigned long stateDescription3);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // totalTz = twice the total Tz value
  // totalY = three time the total Y value
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalTz, int totalY);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
  // nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
  // nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson in the state
  // totalLz = momentum total value
  // nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
  // nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
  // nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int lzMax1, int lzMax2, int lzMax3, int totalLz, int nbrN1, int nbrN2, int nbrN3, long pos);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initial bosonic state is stored
  // return value = corresponding fermionic state
  virtual unsigned long BosonToFermion(unsigned long*& initialState);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState1 = reference on the array where initial bosonic state for the type 1 particles is stored
  // initialState2 = reference on the array where initial bosonic state for the type 2 particles is stored
  // initialState3 = reference on the array where initial bosonic state for the type 3 particles is stored
  // finalState1 = reference on the corresponding fermionic state for the type 1 particles
  // finalState2 = reference on the corresponding fermionic state for the type 2 particles
  // finalState3 = reference on the corresponding fermionic state for the type 3 particles
  virtual void BosonToFermion(unsigned long*& initialState1, unsigned long*& initialState2, unsigned long*& initialState3, 
			      unsigned long& finalState1, unsigned long& finalState2, unsigned long& finalState3);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax= maximum lz value reached by the fermionic state
  // finalState = reference on the array where the bosonic state for the type 1 particles has to be stored
  virtual void FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState1 = initial fermionic state for the type 1 particles
  // initialState2 = initial fermionic state for the type 2 particles
  // initialState3 = initial fermionic state for the type 3 particles
  // finalState1 = reference on the array where the bosonic state for the type 1 particles has to be stored
  // finalState2 = reference on the array where the bosonic state for the type 2 particles has to be stored
  // finalState3 = reference on the array where the bosonic state for the type 3 particles has to be stored
  virtual void FermionToBoson(unsigned long initialState1, unsigned long initialState2, unsigned long initialState3,
			      unsigned long*& finalState1, unsigned long*& finalState2, unsigned long*& finalState3);

  // apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // temporaryStatei= reference on the temporary array for the type i particles
  // temporaryStatej= reference on the temporary array for the type j particles
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state
  virtual int AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej, double& coefficient);

  // find state index
  //
  // stateDescription1 = array describing the bosonic state for type 1 particles
  // stateDescription2 = array describing the bosonic state for type 2 particles
  // stateDescription3 = array describing the bosonic state for type 3 particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long*& stateDescription1, unsigned long*& stateDescription2, unsigned long*& stateDescription3);

  // recursive part of the convertion from a state from one SU(3) basis to another, transforming the one body basis in each momentum sector
  //
  // targetState = vector where the transformed state has to be stored
  // coefficient = current coefficient to assign
  // position = current particle consider in the n-body state
  // momentumIndices = array that gives the momentum partition of the initial n-body state
  // initialSU3Indices = array that gives the spin dressing the initial n-body state
  // currentSU3Indices = array that gives the spin dressing the current transformed n-body state
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
  // occupationCoefficientArray = array that provides 1/2 ln (N!)
  void TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
				      int position, int* momentumIndices, int* initialSU3Indices, int* currentSU3Indices, ComplexMatrix* oneBodyBasis, 
				      double occupationCoefficient, double* occupationCoefficientArray);

};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnSphereWithSU3Spin::GetParticleStatistic()
{
  return AbstractQHEParticle::BosonicStatistic;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initial bosonic state is stored
// return value = corresponding fermionic state

inline unsigned long BosonOnSphereWithSU3Spin::BosonToFermion(unsigned long*& initialState)
{
  unsigned long FinalState = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      FinalState |= ((1ul << initialState[i]) - 1ul) << Shift;
      Shift += initialState[i];
      ++Shift;
    }
  return FinalState;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState1 = reference on the array where initial bosonic state for the type 1 particles is stored
// initialState2 = reference on the array where initial bosonic state for the type 2 particles is stored
// initialState3 = reference on the array where initial bosonic state for the type 3 particles is stored
// finalState1 = reference on the corresponding fermionic state for the type 1 particles
// finalState2 = reference on the corresponding fermionic state for the type 2 particles
// finalState3 = reference on the corresponding fermionic state for the type 3 particles

inline void BosonOnSphereWithSU3Spin::BosonToFermion(unsigned long*& initialState1, unsigned long*& initialState2, unsigned long*& initialState3, 
						     unsigned long& finalState1, unsigned long& finalState2, unsigned long& finalState3)
{
  finalState1 = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalState1 |= ((1ul << initialState1[i]) - 1ul) << Shift;
      Shift += initialState1[i];
      ++Shift;
    }
  finalState2 = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalState2 |= ((1ul << initialState2[i]) - 1ul) << Shift;
      Shift += initialState2[i];
      ++Shift;
    }
  finalState3 = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalState3 |= ((1ul << initialState3[i]) - 1ul) << Shift;
      Shift += initialState3[i];
      ++Shift;
    }
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateLzMax= maximum lz value reached by the fermionic state
// finalState = reference on the array where the bosonic state for the type 1 particles has to be stored

inline void BosonOnSphereWithSU3Spin::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState)
{
  int FinalStateLzMax = 0;
  while ((initialStateLzMax >= 0) && ((initialState >> initialStateLzMax) == 0x0ul))
    --initialStateLzMax;
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
      finalState[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState >>= TmpPower;
      ++FinalStateLzMax;
      initialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalState[FinalStateLzMax] = 0x0ul;
}


// convert a fermionic state into its bosonic  counterpart
//
// initialState1 = initial fermionic state for the type 1 particles
// initialState2 = initial fermionic state for the type 2 particles
// initialState3 = initial fermionic state for the type 3 particles
// finalState1 = reference on the array where the bosonic state for the type 1 particles has to be stored
// finalState2 = reference on the array where the bosonic state for the type 2 particles has to be stored
// finalState3 = reference on the array where the bosonic state for the type 3 particles has to be stored

inline void BosonOnSphereWithSU3Spin::FermionToBoson(unsigned long initialState1, unsigned long initialState2, unsigned long initialState3,
						     unsigned long*& finalState1, unsigned long*& finalState2, unsigned long*& finalState3)
{
  int FinalStateLzMax = 0;
  int InitialStateLzMax = this->N1LzMax;
  while ((InitialStateLzMax >= 0) && ((initialState1 >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialState1 - 1ul) ^ (~initialState1);
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
      finalState1[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState1 >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalState1[FinalStateLzMax] = 0x0ul;

  FinalStateLzMax = 0;
  InitialStateLzMax = this->N2LzMax;
  while ((InitialStateLzMax >= 0) && ((initialState2 >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialState2 - 1ul) ^ (~initialState2);
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
      finalState2[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState2 >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalState2[FinalStateLzMax] = 0x0ul;

  FinalStateLzMax = 0;
  InitialStateLzMax = this->N3LzMax;
  while ((InitialStateLzMax >= 0) && ((initialState3 >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialState3 - 1ul) ^ (~initialState3);
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
      finalState3[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState3 >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalState3[FinalStateLzMax] = 0x0ul;

}

// apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// temporaryStatei= reference on the temporary array for the type i particles
// temporaryStatej= reference on the temporary array for the type j particles
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int BosonOnSphereWithSU3Spin::AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej,
					     double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      this->TemporaryState1[i] = this->ProdATemporaryState1[i];
      this->TemporaryState2[i] = this->ProdATemporaryState2[i];
      this->TemporaryState3[i] = this->ProdATemporaryState3[i];
    }
  ++temporaryStatej[m2];
  coefficient = temporaryStatej[m2];
  ++temporaryStatei[m1];
  coefficient *= temporaryStatei[m1];
  coefficient = sqrt(coefficient);
  return this->FindStateIndex(this->TemporaryState1, this->TemporaryState2, this->TemporaryState3);
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0, 1 or 2
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(3) index for the first annihilation operator
// sigma2 = SU(3) index for the second annihilation operator
// return value =  multiplicative factor 

inline double BosonOnSphereWithSU3Spin::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1LzMax, this->ProdATemporaryState1);
  this->FermionToBoson(this->StateDescription2[index], this->N2LzMax, this->ProdATemporaryState2);
  this->FermionToBoson(this->StateDescription3[index], this->N3LzMax, this->ProdATemporaryState3);
  if ((this->ProdATemporaryStateSigma[sigma1][n1] == 0) || (this->ProdATemporaryStateSigma[sigma2][n2] == 0) || 
      ((n1 == n2) && (sigma1 == sigma2) && (this->ProdATemporaryStateSigma[sigma1][n1] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryStateSigma[sigma2][n2];
  --this->ProdATemporaryStateSigma[sigma2][n2];
  Coefficient *= this->ProdATemporaryStateSigma[sigma1][n1];
  --this->ProdATemporaryStateSigma[sigma1][n1];
  return sqrt(Coefficient);
}

// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(3) index for the first creation operator
// sigma2 = SU(3) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int BosonOnSphereWithSU3Spin::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      this->TemporaryState1[i] = this->ProdATemporaryState1[i];
      this->TemporaryState2[i] = this->ProdATemporaryState2[i];
      this->TemporaryState3[i] = this->ProdATemporaryState3[i];
    }
  ++this->TemporaryStateSigma[sigma2][m2];
  coefficient = this->TemporaryStateSigma[sigma2][m2];
  ++this->TemporaryStateSigma[sigma1][m1];
  coefficient *= this->TemporaryStateSigma[sigma1][m1];
  coefficient = sqrt(coefficient);
  return this->FindStateIndex(this->TemporaryState1, this->TemporaryState2, 
			      this->TemporaryState3);
}

// find state index
//
// stateDescription1 = array describing the bosonic state for type 1 particles
// stateDescription2 = array describing the bosonic state for type 2 particles
// stateDescription3 = array describing the bosonic state for type 3 particles
// return value = corresponding index

inline int BosonOnSphereWithSU3Spin::FindStateIndex(unsigned long*& stateDescription1, unsigned long*& stateDescription2, unsigned long*& stateDescription3)
{
  unsigned long Tmp1;
  unsigned long Tmp2;
  unsigned long Tmp3;
  this->BosonToFermion(stateDescription1, stateDescription2, stateDescription3, Tmp1, Tmp2, Tmp3);
  return this->FindStateIndex(Tmp1, Tmp2, Tmp3);
}

#endif


