////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with SU(3) spin              //
//                                                                            //
//                        last modification : 20/01/2008                      //
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


#ifndef FERMIONONSPHEREWITHSU3SPIN_H
#define FERMIONONSPHEREWITHSU3SPIN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU3Spin.h"

#include <iostream>


class FermionOnSphere;


class FermionOnSphereWithSU3Spin :  public ParticleOnSphereWithSU3Spin
{

  friend class FermionOnSphereWithSU3SpinTzSymmetry;
  friend class FermionOnSphereWithSU3SpinZ3Symmetry;
  friend class FermionOnSphereWithSU3SpinTzZ3Symmetry;

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
  // twice the total Tz value
  int TotalTz;
  // three time the total Y value
  int TotalY;

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
  FermionOnSphereWithSU3Spin ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twice the total spin value
  // totalIsospin = twice the total isospin value
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSU3Spin (int nbrFermions, int totalLz, int lzMax, int totalSpin, int totalIsospin, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSU3Spin(const FermionOnSphereWithSU3Spin& fermions);

  // destructor
  //
  virtual ~FermionOnSphereWithSU3Spin ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSU3Spin& operator = (const FermionOnSphereWithSU3Spin& fermions);

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

  // evaluate wave function in real space using a given basis and only for agiven range of components
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
					int firstComponent, int nbrComponent);                                
  
  // initialize evaluation of wave function in real space using a given basis and only for a given range of components and
  //
  // timeCoherence = true if time coherence has to be used
  virtual void InitializeWaveFunctionEvaluation (bool timeCoherence = false);
  
  // create a U(1) state from an SU(3) state
  //
  // state = vector describing the SU(3) state
  // u1Space = reference on the Hilbert space associated to the U(1) state
  // return value = resulting U(1) state
  virtual RealVector ForgeU1FromSU3(RealVector& state, FermionOnSphere& u1Space);

  protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);


  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // totalTz = twice the total Tz value
  // totalY = three time the total Y value
  // return value = Hilbert space dimension
  virtual int EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalTz, int totalY);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
  // nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
  // nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
  // nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
  // nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3, long pos);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnSphereWithSU3Spin::GetParticleStatistic()
{
  return AbstractQHEParticle::FermionicStatistic;
}

#endif


