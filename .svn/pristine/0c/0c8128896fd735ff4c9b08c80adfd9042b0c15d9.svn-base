////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                    class of bosons on sphere with SU(4) spin               //
//                                                                            //
//                        last modification : 19/12/2011                      //
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


#ifndef BOSONONSPHEREWITHSU4SPINALLENTANGLEMENT_H
#define BOSONONSPHEREWITHSU4SPINALLENTANGLEMENT_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU4Spin.h"

#include <iostream>
using std::cout;
using std::endl;

class BosonOnSphereShort;


class BosonOnSphereWithSU4SpinAllEntanglement : public ParticleOnSphereWithSU4Spin
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
  // twice the total spin value
  int TotalSpin;
  // twice the total isospin value
  int TotalIsospin;

  // effective LzMax in fermionic coding of block I / plus states
  int NPlusLzMax;
  // effective LzMax in fermionic coding of block II / minus states
  int NMinusLzMax;

  
  // old stuffs
  // lzmax value for the fermionic states associated to the type up-plus particles
  int NUpPlusLzMax;
  // lzmax value for the fermionic states associated to the type up-minus particles
  int NUpMinusLzMax;
  // lzmax value for the fermionic states associated to the type down-plus particles
  int NDownPlusLzMax;
  // lzmax value for the fermionic states associated to the type down-minus particles
  int NDownMinusLzMax;
  // largest lzmax value for the fermionic states among N1LzMax, N2LzMax and N3LzMax
  int FermionicLzMax;

  // temporay array describing the type up-plus particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionPlus;
  // temporay array describing the type up-minus particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionMinus;

  // flag indicating whether full look-up table is present
  bool FullLookUp;
  // size of memory used for storage
  long LookUpMemorySize;
  // full look-up table for plus spins (used for small systems)
  unsigned long* LookUpTablePlus;
  // full look-up table for plus spins (used for small systems)
  unsigned long* LookUpTableMinus;

  // look-up for large systems:
  // sorted array that contains each unique configuration for the type plus particles
  unsigned long* UniqueStateDescriptionPlus;
  // number of time each unique configuration for the type plus particles appears in StateDescriptionUpPlus
  int* UniqueStateDescriptionSubArrayIndicesPlus;
  // number of unique configurations for the type plus particles
  long NbrUniqueStateDescriptionPlus;

  // temporary state used when applying operators for type up-plus particles (actual memory)
  unsigned long* TemporaryStatePlus;
  // temporary state used when applying operators for type up-minus particles
  unsigned long* TemporaryStateMinus;

  // temporary state used when applying operators for type up-plus particles (pointers to segments of above)
  unsigned long* TemporaryStateUpPlus;
  // temporary state used when applying operators for type up-minus particles
  unsigned long* TemporaryStateUpMinus;
  // temporary state used when applying operators for type down-plus particles
  unsigned long* TemporaryStateDownPlus;
  // temporary state used when applying operators for type down-minus particles
  unsigned long* TemporaryStateDownMinus;

  // temporary state used when applying ProdA operator for type up-plus particles (actual memory)
  unsigned long* ProdATemporaryStatePlus;
  // temporary state used when applying ProdA operator for type up-minus particles
  unsigned long* ProdATemporaryStateMinus;
  
  // temporary state used when applying ProdA operator for type up-plus particles (pointers to segments of above)
  unsigned long* ProdATemporaryStateUpPlus;
  // temporary state used when applying ProdA operator for type up-minus particles
  unsigned long* ProdATemporaryStateUpMinus;
  // temporary state used when applying ProdA operator for type down-plus particles
  unsigned long* ProdATemporaryStateDownPlus;
  // temporary state used when applying ProdA operator for type down-minus particles
  unsigned long* ProdATemporaryStateDownMinus;

 public:

  // default constructor
  // 
  BosonOnSphereWithSU4SpinAllEntanglement ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // totalSz = twice the total sz projection
  // totalIsospin = twice the total isospin value (number imbalance between Plus and Minus)
  // memory = amount of memory granted for precalculations
  
  BosonOnSphereWithSU4SpinAllEntanglement (int nbrBosons, int totalLz, int lzMax, int totalSz, int totalIsospin,
					 unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereWithSU4SpinAllEntanglement(const BosonOnSphereWithSU4SpinAllEntanglement& bosons);

  // destructor
  //
  virtual ~BosonOnSphereWithSU4SpinAllEntanglement ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereWithSU4SpinAllEntanglement& operator = (const BosonOnSphereWithSU4SpinAllEntanglement& bosons);

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

  // apply a^+_m_dp a_m_dp operator to a given state (only spin down isospin plus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_dm a_m_dm
  virtual double AddpAdp (int index, int m);

  // apply a^+_m_up a_m_up operator to a given state  (only spin up isospin plus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_um a_m_um
  virtual double AdupAup (int index, int m);

  // apply a^+_m_dm a_m_dm operator to a given state (only spin down isospin minus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_dm a_m_dm
  virtual double AddmAdm (int index, int m);

  // apply a^+_m_um a_m_um operator to a given state  (only spin up isospin minus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_um a_m_um
  virtual double AdumAum (int index, int m);

  // apply a^+_m_up a_n_up operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdupAup (int index, int m, int n, double& coefficient);

  // apply a^+_m_up a_n_um operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdupAum (int index, int m, int n, double& coefficient);

  // apply a^+_m_up a_n_dp operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdupAdp (int index, int m, int n, double& coefficient);

  // apply a^+_m_up a_n_dm operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdupAdm (int index, int m, int n, double& coefficient);

  // apply a^+_m_um a_n_up operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdumAup (int index, int m, int n, double& coefficient);

  // apply a^+_m_um a_n_um operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdumAum (int index, int m, int n, double& coefficient);

  // apply a^+_m_um a_n_dp operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdumAdp (int index, int m, int n, double& coefficient);

  // apply a^+_m_um a_n_dm operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdumAdm (int index, int m, int n, double& coefficient);

  // apply a^+_m_dp a_n_up operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddpAup (int index, int m, int n, double& coefficient);

  // apply a^+_m_dp a_n_um operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddpAum (int index, int m, int n, double& coefficient);

  // apply a^+_m_dp a_n_dp operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddpAdp (int index, int m, int n, double& coefficient);

  // apply a^+_m_dp a_n_dm operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddpAdm (int index, int m, int n, double& coefficient);

  // apply a^+_m_dm a_n_up operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddmAup (int index, int m, int n, double& coefficient);

  // apply a^+_m_dm a_n_um operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddmAum (int index, int m, int n, double& coefficient);

  // apply a^+_m_dm a_n_dp operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddmAdp (int index, int m, int n, double& coefficient);

  // apply a^+_m_dm a_n_dm operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddmAdm (int index, int m, int n, double& coefficient);

  // apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up, 1 for um, 2 for dp and 3 for dm 
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // sigma1 = SU(4) index for the first annihilation operator
  // sigma2 = SU(4) index for the second annihilation operator
  // return value =  multiplicative factor 
  virtual double AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2);

  // apply a_n1_up a_n2_up operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AupAup (int index, int n1, int n2);

  // apply a_n1_up a_n2_um operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AupAum (int index, int n1, int n2);

  // apply a_n1_up a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AupAdp (int index, int n1, int n2);

  // apply a_n1_up a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AupAdm (int index, int n1, int n2);

  // apply a_n1_um a_n2_um operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AumAum (int index, int n1, int n2);

  // apply a_n1_um a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AumAdp (int index, int n1, int n2);

  // apply a_n1_um a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AumAdm (int index, int n1, int n2);

  // apply a_n1_dp a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AdpAdp (int index, int n1, int n2);

  // apply a_n1_dp a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AdpAdm (int index, int n1, int n2);

  // apply a_n1_dm a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AdmAdm (int index, int n1, int n2);

  // apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is 0 for up, 1 for um, 2 for dp and 3 for dm 
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // sigma1 = SU(4) index for the first creation operator
  // sigma2 = SU(4) index for the second creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient);

  // apply a^+_m1_up a^+_m2_up operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdupAdup (int m1, int m2, double& coefficient);

  // apply a^+_m1_up a^+_m2_um operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdupAdum (int m1, int m2, double& coefficient);

  // apply a^+_m1_up a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdupAddp (int m1, int m2, double& coefficient);

  // apply a^+_m1_up a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdupAddm (int m1, int m2, double& coefficient);

  // apply a^+_m1_um a^+_m2_um operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdumAdum (int m1, int m2, double& coefficient);

  // apply a^+_m1_um a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdumAddp (int m1, int m2, double& coefficient);

  // apply a^+_m1_um a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdumAddm (int m1, int m2, double& coefficient);

  // apply a^+_m1_dp a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddpAddp (int m1, int m2, double& coefficient);

  // apply a^+_m1_dp a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddpAddm (int m1, int m2, double& coefficient);

  // apply a^+_m1_dm a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddmAddm (int m1, int m2, double& coefficient);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  /*
  // convert a state from one SU(4) basis to another, transforming the one body basis in each momentum sector
  //
  // initialState = state to transform  
  // targetState = vector where the transformed state has to be stored
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  virtual void TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis);

  // compute the transformation matrix from one SU(4) basis to another, transforming the one body basis in each momentum sector
  //
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // return value = transformation matrix
  virtual ComplexMatrix TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis);

  // compute the projection matrix from the SU(4) Hilbert space to an U(1) Hilbert space
  // 
  // targetSpace = pointer to the U(1) Hilbert space
  // type = type of particles that has to be kept (0 for type up-plus, 1 for type up-minus, 2 for type down-plus, 3 for type down-minus)
  // return value = projection matrix
  virtual ComplexMatrix TransformationMatrixSU4ToU1(BosonOnSphereShort* targetSpace, int type = 0);
  */

  //get bosonic state description for state with index index
  //PlusOccupation = reference on array where plus state occupation numbers will be stored
  //MinusOccupation = reference on array where minus state occupation numbers will be stored
  inline void GetBosonicStateDescription(int index, unsigned long *& PlusOccupation, unsigned long *& MinusOccupation );

  //get fermionic state description for state with index index
  //PlusDescription = unsigned long where plus state description will be stored
  //MinusDescription = unsigned long where minus state description will be stored
  inline void GetFermionicStateDescription(int index, unsigned long & PlusDescription, unsigned long & MinusDescription );

  //get total lz of spins for state index
  //TotalLzPlus = int where twice total lz of plus spins will be stored
  //TotalLzMinus = int where twice total lz of minus spins will be stored
  inline void GetTotalLz(int index, int & TotalLzPlus, int & TotalLzMinus);

  //get total sz of plus spins for state index 
  //TotalSzPlus = int where twice total sz of plus spins will be stored
  //TotalSzMinus = int where twice total sz of minus spins will be stored
  inline void GetTotalSz(int index, int & TotalSzPlus, int & TotalSzMinus);


  protected:

  // find state index
  //
  // stateDescriptionPlus = unsigned integer describing the fermionic state for type plus particles
  // stateDescriptionMinus = unsigned integer describing the fermionic state for type minus particles
  // return value = corresponding index
  virtual inline int FindStateIndex(unsigned long stateDescriptionPlus, unsigned long stateDescriptionMinus); 


  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // totalNPlus = total number of particles in block I
  // totalSz = total spin quantum number Sz
  // return value = Hilbert space dimension  
  long EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalNPlus, int totalLz, int totalSz);

  

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons remaining to be placed
  // lzMax = momentum maximum value for a boson in the state 
  // currentLzMax = momentum maximum value for bosons that are still to be placed
  //                (counting from 0 to lzMax for down and lzMax+1 to 2LzMax+1 for up)
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int currentLzMax, int nPlus, int totalLz, int totalSz, int level);


  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons remaining to be placed
  // lzMax = momentum maximum value for a boson in the state 
  // currentLzMax = momentum maximum value for bosons that are still to be placed
  //                (counting from 0 to lzMax for down and lzMax+1 to 2LzMax+1 for up)
  // totalLz = momentum total value
  // nPlus = remaining number of particles in plus / block I
  // nUp = remaining number of particles with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  
  long GenerateStates(int nbrBosons, int currentLzMax, int nPlus, int totalLz, int nUp, long pos);


  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initial bosonic state is stored
  // return value = corresponding fermionic state
  virtual unsigned long BosonToFermion(unsigned long*& initialState);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialStateUpPlus = reference on the array where initial bosonic state for the type up-plus particles is stored
  // initialStateUpMinus = reference on the array where initial bosonic state for the type up-minus particles is stored
  // initialStateDownPlus = reference on the array where initial bosonic state for the type down-plus particles is stored
  // initialStateDownMinus = reference on the array where initial bosonic state for the type down-minus particles is stored
  // finalStateUpPlus = reference on the corresponding fermionic state for the type up-plus particles
  // finalStateUpMinus = reference on the corresponding fermionic state for the type up-minus particles
  // finalStateDownPlus = reference on the corresponding fermionic state for the type down-plus particles
  // finalStateDownMinus = reference on the corresponding fermionic state for the type down-minus particles
  virtual void BosonToFermion(unsigned long*& initialStatePlus, unsigned long*& initialStateMinus, 
			      unsigned long& finalStatePlus, unsigned long& finalStateMinus);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax= maximum lz value reached by the fermionic state
  // finalState = reference on the array where the bosonic state for the type up-plus particles has to be stored
  virtual void FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialStateUpPlus = initial fermionic state for the type up-plus particles
  // initialStateUpMinus = initial fermionic state for the type up-minus particles
  // initialStateDownPlus = initial fermionic state for the type down-plus particles
  // initialStateDownMinus = initial fermionic state for the type down-minus particles
  // finalStateUpPlus = reference on the array where the bosonic state for the type up-plus particles has to be stored
  // finalStateUpMinus = reference on the array where the bosonic state for the type up-minus particles has to be stored
  // finalStateDownPlus = reference on the array where the bosonic state for the type down-plus particles has to be stored
  // finalStateDownMinus = reference on the array where the bosonic state for the type down-minus particles has to be stored
  virtual void FermionToBoson(unsigned long initialStatePlus, unsigned long initialStateMinus,
			      unsigned long*& finalStatePlus, unsigned long*& finalStateMinus);

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
  // stateDescriptionUpPlus = array describing the bosonic state for type up-plus particles
  // stateDescriptionUpMinus = array describing the bosonic state for type up-minus particles
  // stateDescriptionDownPlus = array describing the bosonic state for type down-plus particles
  // stateDescriptionDownMinus = array describing the bosonic state for type down-minus particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long*& stateDescriptionPlus, unsigned long*& stateDescriptionMinus);

  // recursive part of the convertion from a state from one SU(4) basis to another, transforming the one body basis in each momentum sector
  //
  // targetState = vector where the transformed state has to be stored
  // coefficient = current coefficient to assign
  // position = current particle consider in the n-body state
  // momentumIndices = array that gives the momentum partition of the initial n-body state
  // initialSU4Indices = array that gives the spin dressing the initial n-body state
  // currentSU4Indices = array that gives the spin dressing the current transformed n-body state
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
  // occupationCoefficientArray = array that provides 1/2 ln (N!)
  void TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
				      int position, int* momentumIndices, int* initialSU4Indices, int* currentSU4Indices, ComplexMatrix* oneBodyBasis, 
				      double occupationCoefficient, double* occupationCoefficientArray);

};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnSphereWithSU4SpinAllEntanglement::GetParticleStatistic()
{
  return AbstractQHEParticle::BosonicStatistic;
}


// old boson to fermion conversions:

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initial bosonic state is stored
// return value = corresponding fermionic state

inline unsigned long BosonOnSphereWithSU4SpinAllEntanglement::BosonToFermion(unsigned long*& initialState)
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
// initialStatePlus = reference on the array where initial bosonic state for the type up-plus particles is stored
// initialStateMinus = reference on the array where initial bosonic state for the type up-minus particles is stored
// finalStatePlus = reference on the corresponding fermionic state for the type up-plus particles
// finalStateMinus = reference on the corresponding fermionic state for the type up-minus particles

inline void BosonOnSphereWithSU4SpinAllEntanglement::BosonToFermion(unsigned long*& initialStatePlus, unsigned long*& initialStateMinus, 
								    unsigned long& finalStatePlus, unsigned long& finalStateMinus) 
{
  finalStatePlus = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i < 2*(this->NbrLzValue); ++i)
    {
      finalStatePlus |= ((1ul << initialStatePlus[i]) - 1ul) << Shift;
/*       cout << "initialStatePlus["<<i<<"]="<<initialStatePlus[i]<<endl; */
      Shift += initialStatePlus[i];
      ++Shift;
    }
  finalStateMinus = 0x0ul;
  Shift = 0;
  for (int i = 0; i < 2*(this->NbrLzValue); ++i)
    {
      finalStateMinus |= ((1ul << initialStateMinus[i]) - 1ul) << Shift;
/*       cout << "initialStateMinus["<<i<<"]="<<initialStateMinus[i]<<endl; */
      Shift += initialStateMinus[i];
      ++Shift;
    }
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateLzMax= maximum lz value reached by the fermionic state
// finalState = reference on the array where the bosonic state for the type up-plus particles has to be stored

inline void BosonOnSphereWithSU4SpinAllEntanglement::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState)
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
  for (; FinalStateLzMax < (this->NbrLzValue<<1); ++FinalStateLzMax)
    finalState[FinalStateLzMax] = 0x0ul;
}


// convert a fermionic state into its bosonic  counterpart
//
// initialStatePlus = initial fermionic state for the type up-plus particles
// initialStateMinus = initial fermionic state for the type up-minus particles
// finalStatePlus = reference on the array where the bosonic state for the type up-plus particles has to be stored
// finalStateMinus = reference on the array where the bosonic state for the type up-minus particles has to be stored

inline void BosonOnSphereWithSU4SpinAllEntanglement::FermionToBoson(unsigned long initialStatePlus, unsigned long initialStateMinus, 
								    unsigned long*& finalStatePlus, unsigned long*& finalStateMinus)
{
  int FinalStateLzMax = 0;
  int InitialStateLzMax = this->NPlusLzMax;
  while ((InitialStateLzMax >= 0) && ((initialStatePlus >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStatePlus - 1ul) ^ (~initialStatePlus);
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
      finalStatePlus[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStatePlus >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax < (this->NbrLzValue<<1); ++FinalStateLzMax)
    finalStatePlus[FinalStateLzMax] = 0x0ul;

  FinalStateLzMax = 0;
  InitialStateLzMax = this->NMinusLzMax;
  while ((InitialStateLzMax >= 0) && ((initialStateMinus >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateMinus - 1ul) ^ (~initialStateMinus);
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
      finalStateMinus[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateMinus >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax < (this->NbrLzValue<<1); ++FinalStateLzMax)
    finalStateMinus[FinalStateLzMax] = 0x0ul;
}

// apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// temporaryStatei= reference on the temporary array for the type i particles
// temporaryStatej= reference on the temporary array for the type j particles
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int BosonOnSphereWithSU4SpinAllEntanglement::AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej,
							    double& coefficient)
{
  for (int i = 0; i < (this->NbrLzValue<<1); ++i)
    {
      this->TemporaryStatePlus[i] = this->ProdATemporaryStatePlus[i];
      this->TemporaryStateMinus[i] = this->ProdATemporaryStateMinus[i];
    }
  ++temporaryStatei[m1];
  coefficient = temporaryStatei[m1];
  ++temporaryStatej[m2];
  coefficient *= temporaryStatej[m2];
  coefficient = sqrt(coefficient);
  //return this->FindStateIndex(this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, this->TemporaryStateDownPlus, this->TemporaryStateDownMinus);

  return this->FindStateIndex(this->TemporaryStatePlus, this->TemporaryStateMinus);
}

// find state index
//
// stateDescriptionUpPlus = array describing the bosonic state for type up-plus particles
// stateDescriptionUpMinus = array describing the bosonic state for type up-minus particles
// stateDescriptionDownPlus = array describing the bosonic state for type down-plus particles
// stateDescriptionDownMinus = array describing the bosonic state for type down-minus particles
// return value = corresponding index

inline int BosonOnSphereWithSU4SpinAllEntanglement::FindStateIndex(unsigned long*& stateDescriptionPlus, unsigned long*& stateDescriptionMinus)
{
  unsigned long Tmp1;
  unsigned long Tmp2;
/*   cout << "Current state:"<<endl; */
/*   for (int o=0; o<2*NbrLzValue; ++o) */
/*     cout << "("<<stateDescriptionPlus[o]<<","<< stateDescriptionMinus[o]<<") "; */
/*   cout << " done "<<endl; */
  this->BosonToFermion(stateDescriptionPlus, stateDescriptionMinus, Tmp1, Tmp2);
/*   cout << "Tmp1="<<std::hex<<Tmp1; */
/*   cout << "Tmp2="<<Tmp2<<std::dec<<endl; */
  return this->FindStateIndex(Tmp1, Tmp2);
}

// find state index
//
// stateDescriptionPlus = unsigned integer describing the fermionic state for type up-plus particles
// stateDescriptionMinus = unsigned integer describing the fermionic state for type up-minus particles
// return value = corresponding index

int BosonOnSphereWithSU4SpinAllEntanglement::FindStateIndex(unsigned long stateDescriptionPlus, unsigned long stateDescriptionMinus)
{
/*   std::cout << "searching "<< std::hex << stateDescriptionPlus << " " << stateDescriptionMinus << std::dec << std::endl; */
  if (FullLookUp)
    {
      int Result = LookUpTablePlus[stateDescriptionPlus];
      Result += LookUpTableMinus[stateDescriptionMinus];
      return Result;
    }
  else
    {
      int PosMin = 0;
      int PosMax = this->NbrUniqueStateDescriptionPlus - 1;
      int PosMid = (PosMin + PosMax) >> 1;
      //  cout << "entering " << PosMin << " " << PosMax << endl;
      unsigned long CurrentState = this->UniqueStateDescriptionPlus[PosMid];
      while ((PosMax > PosMin) && (CurrentState != stateDescriptionPlus))
	{
	  if (CurrentState > stateDescriptionPlus)
	    {
	      PosMin = PosMid + 1;
	    }
	  else
	    {
	      PosMax = PosMid - 1;
	    } 
	  PosMid = (PosMin + PosMax) >> 1;
	  CurrentState = this->UniqueStateDescriptionPlus[PosMid];
	}
      if (CurrentState != stateDescriptionPlus)
	PosMid = PosMax;
      //  cout << "pass 1 : " << PosMid << endl;
      PosMin = UniqueStateDescriptionSubArrayIndicesPlus[PosMid];
      PosMax = UniqueStateDescriptionSubArrayIndicesPlus[PosMid+1]-1;
      PosMid = (PosMin + PosMax) >> 1;
      //  cout << "entring pass 2 : " << PosMin << " " << PosMax << endl;
      CurrentState = StateDescriptionMinus[PosMid];
      while ((PosMax > PosMin) && (CurrentState != stateDescriptionMinus))
	{
	  if (CurrentState > stateDescriptionMinus)
	    {
	      PosMin = PosMid + 1;
	    }
	  else
	    {
	      PosMax = PosMid - 1;
	    } 
	  PosMid = (PosMin + PosMax) >> 1;
	  CurrentState = StateDescriptionMinus[PosMid];
	}
      if (CurrentState != stateDescriptionMinus)
	return PosMax;
      else
	return PosMid;
    }
}

//get Bosonic description of the state with index index
//PlusStateDescription = reference on array where Bosonic description of the plus state will be stored
//MinusStateDescription = reference on array where Bosonic description of the minus state will be stored
inline void BosonOnSphereWithSU4SpinAllEntanglement::GetBosonicStateDescription(int index, unsigned long *& PlusStateDescription, unsigned long *& MinusStateDescription)
{
  this->FermionToBoson(this->StateDescriptionPlus[index], this->StateDescriptionMinus[index], PlusStateDescription, MinusStateDescription);
}

//get Fermionic description of the state with index index
//PlusDescription = unsigned long where Fermionic description of the plus state will be stored
//MinusDescription = unsigned long where Fermionic description of the minus state will be stored
inline void BosonOnSphereWithSU4SpinAllEntanglement::GetFermionicStateDescription(int index, unsigned long & PlusDescription, unsigned long & MinusDescription)
{
  PlusDescription = this->StateDescriptionPlus[index];
  MinusDescription = this->StateDescriptionMinus[index];
}


//get total lz of plus spins for state index
//TotalLzPlus = int where twice total lz of plus spins will be stored
//TotalLzMinus = int where twice total lz of minus spins will be stored
inline void BosonOnSphereWithSU4SpinAllEntanglement::GetTotalLz(int index, int & TotalLzPlus, int & TotalLzMinus)
{
  for(int i=0; i<2*this->NbrLzValue; i++)
    {
      TemporaryStatePlus[i] = 0x0ul;
      TemporaryStateMinus[i] = 0x0ul;
    }
  this->FermionToBoson(this->StateDescriptionPlus[index], this->StateDescriptionMinus[index], TemporaryStatePlus, TemporaryStateMinus);
  TotalLzPlus=0;
  TotalLzMinus=0;
  for(int lz=0; lz<this->NbrLzValue; lz++)
    {
      TotalLzPlus += (TemporaryStateUpPlus[lz]+TemporaryStateDownPlus[lz])*(2*lz - this->LzMax);
      TotalLzMinus += (TemporaryStateUpMinus[lz]+TemporaryStateDownMinus[lz])*(2*lz - this->LzMax);
    }

}

//get total lz of minus spins for state index
//TotalSzPlus = int where twice total sz of plus spins will be stored
//TotalSzMinus = int where twice total sz of minus spins will be stored
inline void BosonOnSphereWithSU4SpinAllEntanglement::GetTotalSz(int index, int & TotalSzPlus, int & TotalSzMinus)
{
  for(int i=0; i<2*this->NbrLzValue; i++)
    {
      TemporaryStatePlus[i] = 0x0ul;
      TemporaryStateMinus[i] = 0x0ul;
    }
  this->FermionToBoson(this->StateDescriptionPlus[index], this->StateDescriptionMinus[index], TemporaryStatePlus, TemporaryStateMinus);
  TotalSzPlus = 0;
  TotalSzMinus = 0;
  for(int lz=0; lz<this->NbrLzValue; lz++)
    {
      TotalSzPlus += TemporaryStateUpPlus[lz] - TemporaryStateDownPlus[lz];
      TotalSzMinus += TemporaryStateUpMinus[lz] - TemporaryStateDownMinus[lz];
    }
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up, 1 for um, 2 for dp and 3 for dm 
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(4) index for the first annihilation operator
// sigma2 = SU(4) index for the second annihilation operator
// return value =  multiplicative factor 

inline double BosonOnSphereWithSU4SpinAllEntanglement::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  return 0.0;
}


// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is 0 for up, 1 for um, 2 for dp and 3 for dm 
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(4) index for the first creation operator
// sigma2 = SU(4) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int BosonOnSphereWithSU4SpinAllEntanglement::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  return -1;
}


#endif


