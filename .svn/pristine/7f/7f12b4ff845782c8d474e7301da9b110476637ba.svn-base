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


#ifndef BOSONONSPHEREWITHSU4SPIN_H
#define BOSONONSPHEREWITHSU4SPIN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU4Spin.h"

#include <iostream>


using std::cout;
using std::endl;


class BosonOnSphereShort;
class ParticleOnSphereWithSpin;


class BosonOnSphereWithSU4Spin :  public ParticleOnSphereWithSU4Spin
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
  // twice the total entanglement value
  int TotalEntanglement;

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

  // array that contains the state description, the first entry being StateDescriptionUpPlus, the second entry being StateDescriptionUpMinus, third entry StateDescriptionDownPlus and the fourth entry StateDescriptionDownMinus
  unsigned long* StateDescriptionSigma[4];

  // temporay array describing the type up-plus particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionUpPlus;
  // temporay array describing the type up-minus particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionUpMinus;
  // temporay array describing the type down-plus particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionDownPlus;
  // temporay array describing the type down-minus occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionDownMinus;

  // sorted array that contains each unique configuration for the type up-plus particles
  unsigned long* UniqueStateDescriptionUpPlus;
  // number of time each unique configuration for the type up-plus particles appears in StateDescriptionUpPlus
  int* UniqueStateDescriptionSubArraySizeUpPlus;
  // number of unique configurations for the type up-plus particles
  long NbrUniqueStateDescriptionUpPlus;
  // number of unique configurations for the type up-minus particles per unique type up-plus particle configuration 
  int* NbrUniqueStateDescriptionUpMinus;
  // unique configurations for the type up-minus particles per unique type up-plus particle configuration 
  unsigned long** UniqueStateDescriptionUpMinus;
  // number of time each unique configuration for the type up-minus particles appears in StateDescription2 having the same type up-plus particle configucation
  int** UniqueStateDescriptionSubArraySizeUpMinus;
  // first time a unique combination of type up-plus and type up-minus particle configurations appears in the Hilbert space
  int** FirstIndexUniqueStateDescriptionUpMinus;
  // number of unique configurations for the type down-plus particles per unique type up-plus / up-minus particle configuration 
  int** NbrUniqueStateDescriptionDownPlus;
  // unique configurations for the type down-plus particles per unique type up-plus / up-minus particle configuration 
  unsigned long*** UniqueStateDescriptionDownPlus;
  // number of time each unique configuration for the type down-plus particles appears in StateDescription2 having the same type up-plus / up-minus particle configuration
  int*** UniqueStateDescriptionSubArraySizeDownPlus;
  // first time a unique combination of type up-plus, type up-minus and type down minus particle configurations appears in the Hilbert space
  int*** FirstIndexUniqueStateDescriptionDownPlus;

  // temporary state used when applying operators for type up-plus particles
  unsigned long* TemporaryStateUpPlus;
  // temporary state used when applying operators for type up-minus particles
  unsigned long* TemporaryStateUpMinus;
  // temporary state used when applying operators for type down-plus particles
  unsigned long* TemporaryStateDownPlus;
  // temporary state used when applying operators for type down-minus particles
  unsigned long* TemporaryStateDownMinus;
  // array that contains the temporary state, the first entry being TemporaryStateUpPlus, the second entry being TemporaryStateUpMinus, third entry TemporaryStateDownPlus and the fourth entry TemporaryStateDownMinus
  unsigned long* TemporaryStateSigma[4];

  // temporary state used when applying ProdA operator for type up-plus particles
  unsigned long* ProdATemporaryStateUpPlus;
  // temporary state used when applying ProdA operator for type up-minus particles
  unsigned long* ProdATemporaryStateUpMinus;
  // temporary state used when applying ProdA operator for type down-plus particles
  unsigned long* ProdATemporaryStateDownPlus;
  // temporary state used when applying ProdA operator for type down-minus particles
  unsigned long* ProdATemporaryStateDownMinus;
  // array that contains the temporary state used when applying ProdA operator, the first entry being ProdATemporaryStateUpPlus, the second entry being ProdATemporaryStateUpMinus, third entry ProdATemporaryStateDownPlus and the fourth entry ProdATemporaryStateDownMinus
  unsigned long* ProdATemporaryStateSigma[4];

 public:

  // default constructor
  // 
  BosonOnSphereWithSU4Spin ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // totalSpin = twice the total spin value
  // totalIsospin = twice the total isospin value
  // totalEntanglement = twice the total entanglement value
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU4Spin (int nbrBosons, int totalLz, int lzMax, int totalSpin, int totalIsospin, int totalEntanglement, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereWithSU4Spin(const BosonOnSphereWithSU4Spin& bosons);

  // destructor
  //
  virtual ~BosonOnSphereWithSU4Spin ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereWithSU4Spin& operator = (const BosonOnSphereWithSU4Spin& bosons);

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

  // compute the projection matrix from the SU(4) Hilbert space to an SU(2) Hilbert space
  // 
  // targetSpace = pointer to the SU(2) Hilbert space
  // spinUp = index of the component that has to be consider as a spin up
  // spinDown = index of the component that has to be consider as a spin down
  // return value = projection matrix
  virtual ComplexMatrix TransformationMatrixSU4ToSU2(ParticleOnSphereWithSpin* targetSpace, int spinUp = 0, int spinDown = 1);

  // compute the projection matrix from the SU(4) Hilbert space to an U(1) Hilbert space
  // 
  // targetSpace = pointer to the U(1) Hilbert space
  // type = type of particles that has to be kept (0 for type up-plus, 1 for type up-minus, 2 for type down-plus, 3 for type down-minus)
  // return value = projection matrix
  virtual ComplexMatrix TransformationMatrixSU4ToU1(BosonOnSphereShort* targetSpace, int type = 0);

  protected:

  // find state index
  //
  // stateDescriptionUpPlus = unsigned integer describing the fermionic state for type up-plus particles
  // stateDescriptionUpMins = unsigned integer describing the fermionic state for type up-minus particles
  // stateDescriptionDownPlus = unsigned integer describing the fermionic state for type down-plus particles
  // stateDescriptionDownMinus = unsigned integer describing the fermionic state for type down-plus particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescriptionUpPlus, unsigned long stateDescriptionUpMinus, 
			     unsigned long stateDescriptionDownPlus, unsigned long stateDescriptionDownMinus);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // nbrNUpPlus = number of particles with quantum number up-plus
  // nbrNUpMinus = number of particles with quantum number up-minus
  // nbrNDownPlus = number of particles with quantum number down-plus
  // nbrNDownMinus = number of particles with quantum number down-minus
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, 
						    int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMax1 = momentum maximum value for a boson in the state with up-plus
  // lzMax2 = momentum maximum value for a boson in the state with up-minus
  // lzMax3 = momentum maximum value for a boson in the state with down-plus
  // lzMax4 = momentum maximum value for a boson in the state with down-minus
  // totalLz = momentum total value
  // nbrNUpPlus = number of particles with up-plus
  // nbrNUpMinus = number of particles with up-minus
  // nbrNDownPlus = number of particles with down-plus
  // nbrNDownMinus = number of particles with down-minus
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int lzMax1, int lzMax2, int lzMax3, int lzMax4, int totalLz, 
			      int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus, long pos);

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
  virtual void BosonToFermion(unsigned long*& initialStateUpPlus, unsigned long*& initialStateUpMinus, 
			      unsigned long*& initialStateDownPlus, unsigned long*& initialStateDownMius, 
			      unsigned long& finalStateUpPlus, unsigned long& finalStateUpMinus, 
			      unsigned long& finalStateDownPlus, unsigned long& finalStateDownMinus);

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
  virtual void FermionToBoson(unsigned long initialStateUpPlus, unsigned long initialStateUpMinus,
			      unsigned long initialStateDownPlus, unsigned long initialStateDownMius,
			      unsigned long*& finalStateUpPlus, unsigned long*& finalStateUpMinus, 
			      unsigned long*& finalStateDownPlus, unsigned long*& finalStateDownMinus);

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
  virtual int FindStateIndex(unsigned long*& stateDescriptionUpPlus, unsigned long*& stateDescriptionUpMinus, 
			     unsigned long*& stateDescriptionDownPlus, unsigned long*& stateDescriptionDownMinus);

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

inline int BosonOnSphereWithSU4Spin::GetParticleStatistic()
{
  return AbstractQHEParticle::BosonicStatistic;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initial bosonic state is stored
// return value = corresponding fermionic state

inline unsigned long BosonOnSphereWithSU4Spin::BosonToFermion(unsigned long*& initialState)
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
// initialStateUpPlus = reference on the array where initial bosonic state for the type up-plus particles is stored
// initialStateUpMinus = reference on the array where initial bosonic state for the type up-minus particles is stored
// initialStateDownPlus = reference on the array where initial bosonic state for the type down-plus particles is stored
// initialStateDownMinus = reference on the array where initial bosonic state for the type down-minus particles is stored
// finalStateUpPlus = reference on the corresponding fermionic state for the type up-plus particles
// finalStateUpMinus = reference on the corresponding fermionic state for the type up-minus particles
// finalStateDownPlus = reference on the corresponding fermionic state for the type down-plus particles
// finalStateDownMinus = reference on the corresponding fermionic state for the type down-minus particles

inline void BosonOnSphereWithSU4Spin::BosonToFermion(unsigned long*& initialStateUpPlus, unsigned long*& initialStateUpMinus, 
						     unsigned long*& initialStateDownPlus, unsigned long*& initialStateDownMinus, 
						     unsigned long& finalStateUpPlus, unsigned long& finalStateUpMinus, 
						     unsigned long& finalStateDownPlus, unsigned long& finalStateDownMinus)
{
  finalStateUpPlus = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalStateUpPlus |= ((1ul << initialStateUpPlus[i]) - 1ul) << Shift;
      Shift += initialStateUpPlus[i];
      ++Shift;
    }
  finalStateUpMinus = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalStateUpMinus |= ((1ul << initialStateUpMinus[i]) - 1ul) << Shift;
      Shift += initialStateUpMinus[i];
      ++Shift;
    }
  finalStateDownPlus = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalStateDownPlus |= ((1ul << initialStateDownPlus[i]) - 1ul) << Shift;
      Shift += initialStateDownPlus[i];
      ++Shift;
    }
  finalStateDownMinus = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalStateDownMinus |= ((1ul << initialStateDownMinus[i]) - 1ul) << Shift;
      Shift += initialStateDownMinus[i];
      ++Shift;
    }
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateLzMax= maximum lz value reached by the fermionic state
// finalState = reference on the array where the bosonic state for the type up-plus particles has to be stored

inline void BosonOnSphereWithSU4Spin::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState)
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
// initialStateUpPlus = initial fermionic state for the type up-plus particles
// initialStateUpMinus = initial fermionic state for the type up-minus particles
// initialStateDownPlus = initial fermionic state for the type down-plus particles
// initialStateDownMinus = initial fermionic state for the type down-minus particles
// finalStateUpPlus = reference on the array where the bosonic state for the type up-plus particles has to be stored
// finalStateUpMinus = reference on the array where the bosonic state for the type up-minus particles has to be stored
// finalStateDownPlus = reference on the array where the bosonic state for the type down-plus particles has to be stored
// finalStateDownMinus = reference on the array where the bosonic state for the type down-minus particles has to be stored

inline void BosonOnSphereWithSU4Spin::FermionToBoson(unsigned long initialStateUpPlus, unsigned long initialStateUpMinus, 
						     unsigned long initialStateDownPlus, unsigned long initialStateDownMinus,
						     unsigned long*& finalStateUpPlus, unsigned long*& finalStateUpMinus, 
						     unsigned long*& finalStateDownPlus, unsigned long*& finalStateDownMinus)
{
  int FinalStateLzMax = 0;
  int InitialStateLzMax = this->NUpPlusLzMax;
  while ((InitialStateLzMax >= 0) && ((initialStateUpPlus >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateUpPlus - 1ul) ^ (~initialStateUpPlus);
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
      finalStateUpPlus[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateUpPlus >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalStateUpPlus[FinalStateLzMax] = 0x0ul;

  FinalStateLzMax = 0;
  InitialStateLzMax = this->NUpMinusLzMax;
  while ((InitialStateLzMax >= 0) && ((initialStateUpMinus >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateUpMinus - 1ul) ^ (~initialStateUpMinus);
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
      finalStateUpMinus[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateUpMinus >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalStateUpMinus[FinalStateLzMax] = 0x0ul;

  FinalStateLzMax = 0;
  InitialStateLzMax = this->NDownPlusLzMax;
  while ((InitialStateLzMax >= 0) && ((initialStateDownPlus >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateDownPlus - 1ul) ^ (~initialStateDownPlus);
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
      finalStateDownPlus[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateDownPlus >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalStateDownPlus[FinalStateLzMax] = 0x0ul;

  FinalStateLzMax = 0;
  InitialStateLzMax = this->NDownMinusLzMax;
  while ((InitialStateLzMax >= 0) && ((initialStateDownMinus >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateDownMinus - 1ul) ^ (~initialStateDownMinus);
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
      finalStateDownMinus[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateDownMinus >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalStateDownMinus[FinalStateLzMax] = 0x0ul;
}

// apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// temporaryStatei= reference on the temporary array for the type i particles
// temporaryStatej= reference on the temporary array for the type j particles
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int BosonOnSphereWithSU4Spin::AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej,
					     double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      this->TemporaryStateUpPlus[i] = this->ProdATemporaryStateUpPlus[i];
      this->TemporaryStateUpMinus[i] = this->ProdATemporaryStateUpMinus[i];
      this->TemporaryStateDownPlus[i] = this->ProdATemporaryStateDownPlus[i];
      this->TemporaryStateDownMinus[i] = this->ProdATemporaryStateDownMinus[i];
    }
  ++temporaryStatei[m2];
  coefficient = temporaryStatei[m2];
  ++temporaryStatej[m1];
  coefficient *= temporaryStatej[m1];
  coefficient = sqrt(coefficient);
  return this->FindStateIndex(this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, 
			      this->TemporaryStateDownPlus, this->TemporaryStateDownMinus);
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up, 1 for um, 2 for dp and 3 for dm 
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(4) index for the first annihilation operator
// sigma2 = SU(4) index for the second annihilation operator
// return value =  multiplicative factor 

inline double BosonOnSphereWithSU4Spin::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->ProdATemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->ProdATemporaryStateDownMinus);
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
// sigma1 = SU(4) index for the first creation operator
// sigma2 = SU(4) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int BosonOnSphereWithSU4Spin::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      this->TemporaryStateUpPlus[i] = this->ProdATemporaryStateUpPlus[i];
      this->TemporaryStateUpMinus[i] = this->ProdATemporaryStateUpMinus[i];
      this->TemporaryStateDownPlus[i] = this->ProdATemporaryStateDownPlus[i];
      this->TemporaryStateDownMinus[i] = this->ProdATemporaryStateDownMinus[i];
    }
  ++this->TemporaryStateSigma[sigma2][m2];
  coefficient = this->TemporaryStateSigma[sigma2][m2];
  ++this->TemporaryStateSigma[sigma1][m1];
  coefficient *= this->TemporaryStateSigma[sigma1][m1];
  coefficient = sqrt(coefficient);
  return this->FindStateIndex(this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, 
			      this->TemporaryStateDownPlus, this->TemporaryStateDownMinus);
}

// find state index
//
// stateDescriptionUpPlus = array describing the bosonic state for type up-plus particles
// stateDescriptionUpMinus = array describing the bosonic state for type up-minus particles
// stateDescriptionDownPlus = array describing the bosonic state for type down-plus particles
// stateDescriptionDownMinus = array describing the bosonic state for type down-minus particles
// return value = corresponding index

inline int BosonOnSphereWithSU4Spin::FindStateIndex(unsigned long*& stateDescriptionUpPlus, unsigned long*& stateDescriptionUpMinus,
						    unsigned long*& stateDescriptionDownPlus, unsigned long*& stateDescriptionDownMinus)
{
  unsigned long Tmp1;
  unsigned long Tmp2;
  unsigned long Tmp3;
  unsigned long Tmp4;
  this->BosonToFermion(stateDescriptionUpPlus, stateDescriptionUpMinus, stateDescriptionDownPlus, stateDescriptionDownMinus,
		       Tmp1, Tmp2, Tmp3, Tmp4);
  return this->FindStateIndex(Tmp1, Tmp2, Tmp3, Tmp4);
}

#endif


