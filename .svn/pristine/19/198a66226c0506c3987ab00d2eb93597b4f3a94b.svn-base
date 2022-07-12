////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with SU(4) spin              //
//                                                                            //
//                        last modification : 11/10/2006                      //
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


#ifndef FERMIONONSPHEREWITHSU4SPIN_H
#define FERMIONONSPHEREWITHSU4SPIN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU4Spin.h"

#include <iostream>


class FermionOnSphere;
class FermionOnSphereWithSpin;


class FermionOnSphereWithSU4Spin :  public ParticleOnSphereWithSU4Spin
{

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
  // twice the total spin value
  int TotalSpin;
  // twice the total isospin value
  int TotalIsospin;
  // twice the total entanglement value (greater than NbrFermions if there is no constraint on total lEntanglement value)
  int TotalEntanglement;

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

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twice the total spin value
  // totalIsospin = twice the total isospin value
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSU4Spin (int nbrFermions, int totalLz, int lzMax, int totalSpin, int totalIsospin, unsigned long memory = 10000000);

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twice the total spin value
  // totalIsospin = twice the total isospin value
  // totalEntanglement = twice the total entanglement value
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSU4Spin (int nbrFermions, int totalLz, int lzMax, int totalSpin, int totalIsospin, 
			      int totzlEntanglement, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSU4Spin(const FermionOnSphereWithSU4Spin& fermions);

  // destructor
  //
  ~FermionOnSphereWithSU4Spin ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSU4Spin& operator = (const FermionOnSphereWithSU4Spin& fermions);

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

  // apply a^+_m_dp a_m_dp operator to a given state (only spin down isospin plus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_dm a_m_dm
  double AddpAdp (int index, int m);

  // apply a^+_m_up a_m_up operator to a given state  (only spin up isospin plus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_um a_m_um
  double AdupAup (int index, int m);

  // apply a^+_m_dm a_m_dm operator to a given state (only spin down isospin minus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_dm a_m_dm
  double AddmAdm (int index, int m);

  // apply a^+_m_um a_m_um operator to a given state  (only spin up isospin minus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_um a_m_um
  double AdumAum (int index, int m);

  // apply a_n1_up a_n2_up operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AupAup (int index, int n1, int n2);

  // apply a_n1_up a_n2_um operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AupAum (int index, int n1, int n2);

  // apply a_n1_up a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AupAdp (int index, int n1, int n2);

  // apply a_n1_up a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AupAdm (int index, int n1, int n2);

  // apply a_n1_um a_n2_um operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AumAum (int index, int n1, int n2);

  // apply a_n1_um a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AumAdp (int index, int n1, int n2);

  // apply a_n1_um a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AumAdm (int index, int n1, int n2);

  // apply a_n1_dp a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AdpAdp (int index, int n1, int n2);

  // apply a_n1_dp a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AdpAdm (int index, int n1, int n2);

  // apply a_n1_dm a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  double AdmAdm (int index, int n1, int n2);

  // apply a^+_m1_up a^+_m2_up operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdupAdup (int m1, int m2, double& coefficient);

  // apply a^+_m1_up a^+_m2_um operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdupAdum (int m1, int m2, double& coefficient);

  // apply a^+_m1_up a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdupAddp (int m1, int m2, double& coefficient);

  // apply a^+_m1_up a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdupAddm (int m1, int m2, double& coefficient);

  // apply a^+_m1_um a^+_m2_um operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdumAdum (int m1, int m2, double& coefficient);

  // apply a^+_m1_um a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdumAddp (int m1, int m2, double& coefficient);

  // apply a^+_m1_um a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdumAddm (int m1, int m2, double& coefficient);

  // apply a^+_m1_dp a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddpAddp (int m1, int m2, double& coefficient);

  // apply a^+_m1_dp a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddpAddm (int m1, int m2, double& coefficient);

  // apply a^+_m1_dm a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddmAddm (int m1, int m2, double& coefficient);

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
  virtual RealVector ForgeU1FromSU4(RealVector& state, FermionOnSphere& u1Space);

  // create a SU(2) state from an SU(4) state (fusing same spin values,i.e symmetrizing over the isospin)
  //
  // state = vector describing the SU(4) state
  // su2Space = reference on the Hilbert space associated to the SU(2) state
  // return value = resulting SU(2) state
  virtual RealVector ForgeSU2FromSU4(RealVector& state, FermionOnSphereWithSpin& su2Space);

  private:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription, int lzmax);


  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // totalSpin = twice the total spin value
  // totalIsospin = twice the total isospin value
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin);

  // evaluate Hilbert space dimension with shifted values for lzMax and totalLz
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
  // totalSpin = number of particles with spin up
  // totalIsospin = number of particles with isospin plus
  // return value = Hilbert space dimension  
  long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin);

  // evaluate Hilbert space dimension with shifted values for lzMax and totalLz
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
  // totalSpin = number of particles with spin up
  // totalIsospin = number of particles with isospin plus
  // entanglement = number of particles with entanglement plus
  // return value = Hilbert space dimension  
  long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin, int entanglement);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // totalSpin = number of particles with spin up
  // totalIsospin = number of particles with isospin plus
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // totalSpin = number of particles with spin up
  // totalIsospin = number of particles with isospin plus
  // totalEntanglement = number of particles with entanglement plus
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin, int totalEntanglement, long pos);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnSphereWithSU4Spin::GetParticleStatistic()
{
  return AbstractQHEParticle::FermionicStatistic;
}

#endif


