////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of boson with SU4 spin on a torus taking               //
//                    into account magnetic translations                      //
//                                                                            //
//                        last modification : 21/06/2012                      //
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


#ifndef BOSONONTORUSWITHSU4SPINANDMAGNETICTRANSLATIONS_H
#define BOSONONTORUSWITHSU4SPINANDMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSU4SpinAndMagneticTranslations.h"

#include <cmath>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class BosonOnTorusWithSU4SpinAndMagneticTranslations :  public ParticleOnTorusWithSU4SpinAndMagneticTranslations
{

 protected:

  // number of bosons
  int NbrBosons;
  // number of bosons plus 1
  int IncNbrBosons;
  // number of flux quanta
  int MaxMomentum;
  // maximum value for the Ky momentum
  int KyMax;
  // number of flux quanta for the fermionic representation
  int FermionicMaxMomentum;
  // number of Ky values in a state
  int NbrKyValue;
  // twice the total spin value
  int TotalSpin;
  // twice the total isospin value
  int TotalIsospin;
  // twice the total entanglement value
  int TotalEntanglement;

  // number of bosons of type up-plus
  int NbrBosonsUpPlus;
  // number of bosons of type up-minus
  int NbrBosonsUpMinus;
  // number of bosons of type down-plus
  int NbrBosonsDownPlus;
  // number of bosons of type down-minus
  int NbrBosonsDownMinus;

  // momentum value in the x direction (modulo GCD of nbrBosons and maxMomentum)
  int KxMomentum;
  // momentum value in the y direction (modulo GCD of nbrBosons and maxMomentum)
  int KyMomentum;
  //  GCD of nbrBosons and maxMomentum
  int MomentumModulo;
  // translation step used for the magnetic translation along x 
  int XMomentumTranslationStep;

  // value that has to be substracted to the momentum for each translation of the canonical form research
  int MomentumIncrement;
  // shift that has to be done on a state for each translation of the canonical form research
  int StateShift;
  // complementary shift (with respect to MaxMomentum) to StateShift
  int ComplementaryStateShift;
  // mask that corresponds to last bit that can be set to one for type up-plus particles 
  unsigned long LastMomentumMaskUpPlus;
  // mask that corresponds to last bit that can be set to one for type up-minus particles 
  unsigned long LastMomentumMaskUpMinus;
  // mask that corresponds to last bit that can be set to one for type down-plus particles 
  unsigned long LastMomentumMaskDownPlus;
  // mask that corresponds to last bit that can be set to one for type down-minus particles 
  unsigned long LastMomentumMaskDownMinus;
  // mask corresponding to StateShift
  unsigned long MomentumMask;

  // kymax value for the fermionic states associated to the type up-plus particles
  int NUpPlusKyMax;
  // kymax value for the fermionic states associated to the type up-minus particles
  int NUpMinusKyMax;
  // kymax value for the fermionic states associated to the type down-plus particles
  int NDownPlusKyMax;
  // kymax value for the fermionic states associated to the type down-minus particles
  int NDownMinusKyMax;
  // largest kymax value for the fermionic states among NUpPlusKymax, NUpMinusKymax and NDownPlusKymax
  int FermionicKyMax;

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

  // temporary state used when applying ProdA operator for type up-plus particles
  unsigned long* ProdATemporaryStateUpPlus;
  // temporary state used when applying ProdA operator for type up-minus particles
  unsigned long* ProdATemporaryStateUpMinus;
  // temporary state used when applying ProdA operator for type down-plus particles
  unsigned long* ProdATemporaryStateDownPlus;
  // temporary state used when applying ProdA operator for type down-minus particles
  unsigned long* ProdATemporaryStateDownMinus;
  // temporay pointer to the rescaling factor of the state currently processed by ProdA
  double* ProdANbrStateInOrbit;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

 public:

  
  // basic constructor
  // 
  // nbrBosons= number of bosons
  // totalSpin = twice the total spin value
  // totalIsospin = twice the total isospin value
  // totalEntanglement = twice the total entanglement value
  // maxMomentum = momentum maximum value for a boson
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)
  BosonOnTorusWithSU4SpinAndMagneticTranslations (int nbrBosons, int totalSpin, int totalIsospin, int totalEntanglement, int maxMomentum, int xMomentum, int yMomentum);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnTorusWithSU4SpinAndMagneticTranslations(const BosonOnTorusWithSU4SpinAndMagneticTranslations& bosons);

  // destructor
  //
  ~BosonOnTorusWithSU4SpinAndMagneticTranslations();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnTorusWithSU4SpinAndMagneticTranslations& operator = (const BosonOnTorusWithSU4SpinAndMagneticTranslations& bosons);

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

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // apply a^+_m_up a_m_up operator to a given state (only state 1 Tz=+1/2, Y=+1/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_up a_m_up
  virtual double AdupAup (int index, int m);

  // apply a^+_m_um a_m_um operator to a given state (only state 2 Tz=-1/2, Y=+1/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_um a_m_um
  virtual double AdumAum (int index, int m);

  // apply a^+_m_dp a_m_dp operator to a given state (only state 3 Tz=0, Y=-2/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_dp a_m_dp
  virtual double AddpAdp (int index, int m);

  // apply a^+_m_dm a_m_dm operator to a given state (only state 3 Tz=0, Y=-2/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_dp a_m_dp
  virtual double AddmAdm (int index, int m);

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

  // apply a^+_m1_up a^+_m2_up operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int AdupAdup (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_up a^+_m2_um operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int AdupAdum (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_up a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int AdupAddp (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_up a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int AdupAddm (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_um a^+_m2_um operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int AdumAdum (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_um a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int AdumAddp (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_um a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int AdumAddm (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_dp a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int AddpAddp (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_dp a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int AddpAddm (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_dm a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int AddmAddm (int m1, int m2, double& coefficient, int& nbrTranslation);

 protected:

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
			      unsigned long*& initialStateDownPlus, unsigned long*& initialStateDownMinus,
			      unsigned long& finalStateUpPlus, unsigned long& finalStateUpMinus, 
			      unsigned long& finalStateDownPlus, unsigned long& finalStateDownMius);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateKyMax= maximum lz value reached by the fermionic state
  // finalState = reference on the array where the bosonic state for the type up particles has to be stored
  virtual void FermionToBoson(unsigned long initialState, int initialStateKyMax, unsigned long*& finalState);

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
			      unsigned long initialStateDownPlus, unsigned long initialStateDownMinus,
			      unsigned long*& finalStateUpPlus, unsigned long*& finalStateUpMinus, 
			      unsigned long*& finalStateDownPlus, unsigned long*& finalStateDownMinus);

  // apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // temporaryStatei= reference on the temporary array for the type i particles
  // temporaryStatej= reference on the temporary array for the type j particles
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = index of the destination state
  virtual int AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej, double& coefficient, int& nbrTranslation);

  // find canonical form of a state description
  //
  // stateDescriptionUpPlus = reference on the unsigned integer describing the state related to type up-plus particles
  // stateDescriptionUpMinus = reference on the unsigned integer describing the state related to type up-minus particles
  // stateDescriptionDownPlus = reference on the unsigned integer describing the state related to type down-plus particles
  // stateDescriptionDownMinus = reference on the unsigned integer describing the state related to type down-minus particles
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = canonical form of a state description (fermionic representation)
  virtual void FindCanonicalForm(unsigned long& stateDescriptionUpPlus, unsigned long& stateDescriptionUpMinus, 
				 unsigned long& stateDescriptionDownPlus, unsigned long& stateDescriptionDownMinus, int& nbrTranslation);
  
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescriptionUpPlus = reference on the unsigned integer describing the state related to type up-plus particles
  // stateDescriptionUpMinus = reference on the unsigned integer describing the state related to type up-minus particles
  // stateDescriptionDownPlus = reference on the unsigned integer describing the state related to type down-plus particles
  // stateDescriptionDownMinus = reference on the unsigned integer describing the state related to type down-minus particles
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = true if the state does not fit the x momentum constraint
  virtual bool FindCanonicalFormAndTestXMomentumConstraint(unsigned long& stateDescriptionUpPlus, unsigned long& stateDescriptionUpMinus, 
							   unsigned long& stateDescriptionDownPlus, unsigned long& stateDescriptionDownMinus, int& nbrTranslation);
  
  // find how many translations on the x direction are needed to obtain the same state
  //
  // stateDescriptionUpPlus = unsigned integer describing the state for type up-plus particles
  // stateDescriptionUpMinus = unsigned integer describing the state for type up-minus particles
  // stateDescriptionDownPlus = unsigned integer describing the state for type down-plus particles
  // stateDescriptionDownMinus = unsigned integer describing the state for type down-minus particles
  // return value = number of translation needed to obtain the same state
  virtual int FindNumberXTranslation(unsigned long stateDescriptionUpPlus, unsigned long stateDescriptionUpMinus, 
				     unsigned long stateDescriptionDownPlus, unsigned long stateDescriptionDownMinus);
  
  // test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescriptionUpPlus = unsigned integer describing the state related to type up-plus particles
  // stateDescriptionUpMinus = unsigned integer describing the state related to type up-minus particles
  // stateDescriptionDownPlus = unsigned integer describing the state related to type down-plus particles
  // stateDescriptionDownMinus = unsigned integer describing the state related to type down-minus particles
  // return value = true if the state satisfy the x momentum constraint  
  virtual bool TestXMomentumConstraint(unsigned long stateDescriptionUpPlus, unsigned long stateDescriptionUpMinus, 
				       unsigned long stateDescriptionDownPlus, unsigned long stateDescriptionDownMinus);

  // apply a single translation to a bosonic state in its fermionic representation
  //
  // stateDescriptionUpPlus = reference on the unsigned integer describing the state related to type up-plus particles
  // stateDescriptionUpMinus = reference on the unsigned integer describing the state related to type up-minus particles
  // stateDescriptionDownPlus = reference on the unsigned integer describing the state related to type down-plus particles
  // stateDescriptionDownMinus = reference on the unsigned integer describing the state related to type down-minus particles
  virtual void ApplySingleTranslation(unsigned long& stateDescriptionUpPlus, unsigned long& stateDescriptionUpMinus, 
				      unsigned long& stateDescriptionDownPlus, unsigned long& stateDescriptionDownMinus);

  // find state index
  //
  // stateDescriptionUpPlus = unsigned integer describing the fermionic state for type up-plus particles
  // stateDescriptionUpMinus = unsigned integer describing the fermionic state for type up-minus particles
  // stateDescriptionDownPlus = unsigned integer describing the fermionic state for type down-plus particles
  // stateDescriptionDownMinus = unsigned integer describing the fermionic state for type down-minus particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescriptionUpPlus, unsigned long stateDescriptionUpMinus, 
			     unsigned long stateDescriptionDownPlus, unsigned long stateDescriptionDownMinus);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // return value = hilbert space dimension
  long GenerateStates();
 
  // evaluate Hilbert space dimension for a given total spin momentum
  //
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // nbrNUpPlus = number of particles with quantum number up-plus
  // nbrNUpMinus = number of particles with quantum number up-minus
  // nbrNDownPlus = number of particles with quantum number down-plus
  // nbrNDownMinus = number of particles with quantum number down-minus
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus);

  // generate all states corresponding to the constraints without the mangetic translations
  // 
  // nbrBosons = number of bosons
  // currentKyUpPlus = current momentum along y for a single type up-plus particle
  // currentKyUpMinus = current momentum along y for a single type up-minus particle
  // currentKyDownPlus = current momentum along y for a single type down-plus particle
  // currentKyDownMinus = current momentum along y for a single type down-minus particle
  // currentTotalKy = current total momentum along y
  // nbrNUpPlus = number of particles with quantum number up-plus
  // nbrNUpMinus = number of particles with quantum number up-minus
  // nbrNDownPlus = number of particles with quantum number down-plus
  // nbrNDownMinus = number of particles with quantum number down-minus
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long RawGenerateStates(int nbrBosons, int currentKyUpPlus, int currentKyUpMinus, int currentKyDownPlus, int currentKyDownMinus, int currentTotalKy, 
			 int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus, long pos);

};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnTorusWithSU4SpinAndMagneticTranslations::GetParticleStatistic()
{
  return ParticleOnTorusWithSU4SpinAndMagneticTranslations::BosonicStatistic;
}

// find how many translations on the x direction are needed to obtain the same state
//
// stateDescriptionUpPlus = unsigned integer describing the state for type up-plus particles
// stateDescriptionUpMinus = unsigned integer describing the state for type up-minus particles
// stateDescriptionDownPlus = unsigned integer describing the state for type down-plus particles
// stateDescriptionDownMinus = unsigned integer describing the state for type down-minus particles
// return value = number of translation needed to obtain the same state

inline int BosonOnTorusWithSU4SpinAndMagneticTranslations::FindNumberXTranslation(unsigned long stateDescriptionUpPlus, unsigned long stateDescriptionUpMinus, 
										  unsigned long stateDescriptionDownPlus, unsigned long stateDescriptionDownMinus)
{
  unsigned long TmpStateUpPlus = stateDescriptionUpPlus;
  unsigned long TmpStateUpMinus = stateDescriptionUpMinus;
  unsigned long TmpStateDownPlus = stateDescriptionDownPlus;
  unsigned long TmpStateDownMinus = stateDescriptionDownMinus;
  this->ApplySingleTranslation(TmpStateUpPlus, TmpStateUpMinus, TmpStateDownPlus, TmpStateDownMinus);
  int Index = 1;  
  while ((TmpStateUpPlus != stateDescriptionUpPlus) || (TmpStateUpMinus != stateDescriptionUpMinus) || 
	 (TmpStateDownPlus != stateDescriptionDownPlus) || (TmpStateDownMinus != stateDescriptionDownMinus))
    {
      this->ApplySingleTranslation(TmpStateUpPlus, TmpStateUpMinus, TmpStateDownPlus, TmpStateDownMinus);
      ++Index;  
    }
  return Index;
}

// find canonical form of a state description
//
// stateDescriptionUpPlus = reference on the unsigned integer describing the state related to type up-plus particles
// stateDescriptionUpMinus = reference on the unsigned integer describing the state related to type up-minus particles
// stateDescriptionDownPlus = reference on the unsigned integer describing the state related to type down-plus particles
// stateDescriptionDownMinus = reference on the unsigned integer describing the state related to type down-minus particles
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description (fermionic representation)

inline void BosonOnTorusWithSU4SpinAndMagneticTranslations::FindCanonicalForm(unsigned long& stateDescriptionUpPlus, unsigned long& stateDescriptionUpMinus, 
									      unsigned long& stateDescriptionDownPlus, unsigned long& stateDescriptionDownMinus,
									      int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalStateUpPlus = stateDescriptionUpPlus;
  unsigned long CanonicalStateUpMinus = stateDescriptionUpMinus;
  unsigned long CanonicalStateDownPlus = stateDescriptionDownPlus;
  unsigned long CanonicalStateDownMinus = stateDescriptionDownMinus;
  for (int i = 0; i < this->MomentumModulo; ++i)
    {
      this->ApplySingleTranslation(stateDescriptionUpPlus, stateDescriptionUpMinus, stateDescriptionDownPlus, stateDescriptionDownMinus);
      if ((stateDescriptionUpPlus < CanonicalStateUpPlus) || 
	  ((stateDescriptionUpPlus == CanonicalStateUpPlus) && ((stateDescriptionUpMinus < CanonicalStateUpMinus) || ((stateDescriptionUpMinus == CanonicalStateUpMinus) && 
														      ((stateDescriptionDownPlus < CanonicalStateDownPlus) || 
														       ((stateDescriptionDownPlus == CanonicalStateDownPlus) && (stateDescriptionDownMinus < CanonicalStateDownMinus)))))))
	{
	  CanonicalStateUpPlus = stateDescriptionUpPlus;
	  CanonicalStateUpMinus = stateDescriptionUpMinus;
	  CanonicalStateDownPlus = stateDescriptionDownPlus;
	  CanonicalStateDownMinus = stateDescriptionDownMinus;
	  nbrTranslation = i;
	}
    }
}

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescriptionUpPlus = reference on the unsigned integer describing the state related to type up-plus particles
// stateDescriptionUpMinus = reference on the unsigned integer describing the state related to type up-minus particles
// stateDescriptionDownPlus = reference on the unsigned integer describing the state related to type down-plus particles
// stateDescriptionDownMinus = reference on the unsigned integer describing the state related to type down-minus particles
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = true if the state does not fit the x momentum constraint

inline bool BosonOnTorusWithSU4SpinAndMagneticTranslations::FindCanonicalFormAndTestXMomentumConstraint(unsigned long& stateDescriptionUpPlus, unsigned long& stateDescriptionUpMinus, 
													unsigned long& stateDescriptionDownPlus, unsigned long& stateDescriptionDownMinus, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalStateUpPlus = stateDescriptionUpPlus;
  unsigned long InitialStateDescriptionUpPlus = stateDescriptionUpPlus;
  unsigned long CanonicalStateUpMinus = stateDescriptionUpMinus;
  unsigned long InitialStateDescriptionUpMinus = stateDescriptionUpMinus;
  unsigned long CanonicalStateDownPlus = stateDescriptionDownPlus;
  unsigned long InitialStateDescriptionDownPlus = stateDescriptionDownPlus;
  unsigned long CanonicalStateDownMinus = stateDescriptionDownMinus;
  unsigned long InitialStateDescriptionDownMinus = stateDescriptionDownMinus;
  int Index = 0;
  int OrbitSize = 0;
  while (Index < this->MomentumModulo)
    {
      this->ApplySingleTranslation(stateDescriptionUpPlus, stateDescriptionUpMinus, stateDescriptionDownPlus, stateDescriptionDownMinus);
      ++Index;
      ++OrbitSize;
      if ((stateDescriptionUpPlus != InitialStateDescriptionUpPlus) || (stateDescriptionUpMinus != InitialStateDescriptionUpMinus) || (stateDescriptionDownPlus != InitialStateDescriptionDownPlus) || (stateDescriptionDownMinus != InitialStateDescriptionDownMinus))
	{
	  if ((stateDescriptionUpPlus < CanonicalStateUpPlus) || 
	      ((stateDescriptionUpPlus == CanonicalStateUpPlus) && ((stateDescriptionUpMinus < CanonicalStateUpMinus) || ((stateDescriptionUpMinus == CanonicalStateUpMinus) && ((stateDescriptionDownPlus < CanonicalStateDownPlus) || ((stateDescriptionDownPlus == CanonicalStateDownPlus) && (stateDescriptionDownMinus < CanonicalStateDownMinus)))))))
	    {
	      CanonicalStateUpPlus = stateDescriptionUpPlus;
	      CanonicalStateUpMinus = stateDescriptionUpMinus;
	      CanonicalStateDownPlus = stateDescriptionDownPlus;
	      CanonicalStateDownMinus = stateDescriptionDownMinus;
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
      return false;
    }
  stateDescriptionUpPlus = CanonicalStateUpPlus;
  stateDescriptionUpMinus = CanonicalStateUpMinus;
  stateDescriptionDownPlus = CanonicalStateDownPlus;
  stateDescriptionDownMinus = CanonicalStateDownMinus;
  return true;
}


// apply a single translation to a bosonic state in its fermionic representation
//
// stateDescriptionUpPlus = reference on the unsigned integer describing the state related to type up-plus particles
// stateDescriptionUpMinus = reference on the unsigned integer describing the state related to type up-minus particles
// stateDescriptionDownPlus = reference on the unsigned integer describing the state related to type down-plus particles
// stateDescriptionDownMinus = reference on the unsigned integer describing the state related to type down-minus particles

inline void BosonOnTorusWithSU4SpinAndMagneticTranslations::ApplySingleTranslation(unsigned long& stateDescriptionUpPlus, unsigned long& stateDescriptionUpMinus,
										   unsigned long& stateDescriptionDownPlus, unsigned long& stateDescriptionDownMinus)
{
  for (int i = 0; i < this->StateShift;)
    {
      while ((i < this->StateShift) && ((stateDescriptionUpPlus & 0x1ul) == 0x0ul))
	{
	  stateDescriptionUpPlus >>= 1;
	  ++i;
	}
      if (i < this->StateShift)
	{
	  while ((stateDescriptionUpPlus & 0x1ul) == 0x1ul)
	    {
	      stateDescriptionUpPlus >>= 1;
	      stateDescriptionUpPlus |= this->LastMomentumMaskUpPlus;
	    }
	  stateDescriptionUpPlus >>= 1;	  
	  ++i;
	}
    }
  for (int i = 0; i < this->StateShift;)
    {
      while ((i < this->StateShift) && ((stateDescriptionUpMinus & 0x1ul) == 0x0ul))
	{
	  stateDescriptionUpMinus >>= 1;
	  ++i;
	}
      if (i < this->StateShift)
	{
	  while ((stateDescriptionUpMinus & 0x1ul) == 0x1ul)
	    {
	      stateDescriptionUpMinus >>= 1;
	      stateDescriptionUpMinus |= this->LastMomentumMaskUpMinus;
	    }
	  stateDescriptionUpMinus >>= 1;	  
	  ++i;
	}
    }
  for (int i = 0; i < this->StateShift;)
    {
      while ((i < this->StateShift) && ((stateDescriptionDownPlus & 0x1ul) == 0x0ul))
	{
	  stateDescriptionDownPlus>>= 1;
	  ++i;
	}
      if (i < this->StateShift)
	{
	  while ((stateDescriptionDownPlus & 0x1ul) == 0x1ul)
	    {
	      stateDescriptionDownPlus >>= 1;
	      stateDescriptionDownPlus |= this->LastMomentumMaskDownPlus;
	    }
	  stateDescriptionDownPlus >>= 1;	  
	  ++i;
	}
    }
  for (int i = 0; i < this->StateShift;)
    {
      while ((i < this->StateShift) && ((stateDescriptionDownMinus & 0x1ul) == 0x0ul))
	{
	  stateDescriptionDownMinus>>= 1;
	  ++i;
	}
      if (i < this->StateShift)
	{
	  while ((stateDescriptionDownMinus & 0x1ul) == 0x1ul)
	    {
	      stateDescriptionDownMinus >>= 1;
	      stateDescriptionDownMinus |= this->LastMomentumMaskDownMinus;
	    }
	  stateDescriptionDownMinus >>= 1;	  
	  ++i;
	}
    }
}

// test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescriptionUpPlus = unsigned integer describing the state related to type up-plus particles
// stateDescriptionUpMinus = unsigned integer describing the state related to type up-minus particles
// stateDescriptionDownPlus = unsigned integer describing the state related to type down-plus particles
// stateDescriptionDownMinus = unsigned integer describing the state related to type down-minus particles
// return value = true if the state satisfy the x momentum constraint

inline bool BosonOnTorusWithSU4SpinAndMagneticTranslations::TestXMomentumConstraint(unsigned long stateDescriptionUpPlus, unsigned long stateDescriptionUpMinus,
										    unsigned long stateDescriptionDownPlus, unsigned long stateDescriptionDownMinus)
{
  if (((this->KxMomentum * this->FindNumberXTranslation(stateDescriptionUpPlus, stateDescriptionUpMinus, stateDescriptionDownPlus, stateDescriptionDownMinus)) % this->MomentumModulo) == 0)
    return true;
  else
    return false;
}

// apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// temporaryStatei= reference on the temporary array for the type i particles
// temporaryStatej= reference on the temporary array for the type j particles
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = index of the destination state

inline int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej, double& coefficient, int& nbrTranslation)
{
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      this->TemporaryStateUpPlus[i] = this->ProdATemporaryStateUpPlus[i];
      this->TemporaryStateUpMinus[i] = this->ProdATemporaryStateUpMinus[i];
      this->TemporaryStateDownPlus[i] = this->ProdATemporaryStateDownPlus[i];
      this->TemporaryStateDownMinus[i] = this->ProdATemporaryStateDownMinus[i];
    }
  ++temporaryStatej[m2];
  coefficient = temporaryStatej[m2];
  ++temporaryStatei[m1];
  coefficient *= temporaryStatei[m1];
  coefficient = sqrt(coefficient);
  unsigned long TmpStateUpPlus;
  unsigned long TmpStateUpMinus;
  unsigned long TmpStateDownPlus;
  unsigned long TmpStateDownMinus;
  this->BosonToFermion(this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, this->TemporaryStateDownPlus, this->TemporaryStateDownMinus, 
		       TmpStateUpPlus, TmpStateUpMinus, TmpStateDownPlus, TmpStateDownMinus);
  if (this->FindCanonicalFormAndTestXMomentumConstraint(TmpStateUpPlus, TmpStateUpMinus, TmpStateDownPlus, TmpStateDownMinus, nbrTranslation) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpStateUpPlus, TmpStateUpMinus, TmpStateDownPlus, TmpStateDownMinus);
  coefficient *= this->ProdANbrStateInOrbit[this->NbrStateInOrbit[TmpIndex]];
  nbrTranslation *= this->StateShift;
  return TmpIndex;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initial bosonic state is stored
// return value = corresponding fermionic state

inline unsigned long BosonOnTorusWithSU4SpinAndMagneticTranslations::BosonToFermion(unsigned long*& initialState)
{
  unsigned long FinalState = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->KyMax; ++i)
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

inline void BosonOnTorusWithSU4SpinAndMagneticTranslations::BosonToFermion(unsigned long*& initialStateUpPlus, unsigned long*& initialStateUpMinus, 
									   unsigned long*& initialStateDownPlus, unsigned long*& initialStateDownMinus,
									   unsigned long& finalStateUpPlus, unsigned long& finalStateUpMinus, 
									   unsigned long& finalStateDownPlus, unsigned long& finalStateDownMinus)
{
  finalStateUpPlus = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->KyMax; ++i)
    {
      finalStateUpPlus |= ((1ul << initialStateUpPlus[i]) - 1ul) << Shift;
      Shift += initialStateUpPlus[i];
      ++Shift;
    }
  finalStateUpMinus = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->KyMax; ++i)
    {
      finalStateUpMinus |= ((1ul << initialStateUpMinus[i]) - 1ul) << Shift;
      Shift += initialStateUpMinus[i];
      ++Shift;
    }
  finalStateDownPlus = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->KyMax; ++i)
    {
      finalStateDownPlus |= ((1ul << initialStateDownPlus[i]) - 1ul) << Shift;
      Shift += initialStateDownPlus[i];
      ++Shift;
    }
  finalStateDownMinus = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->KyMax; ++i)
    {
      finalStateDownMinus |= ((1ul << initialStateDownMinus[i]) - 1ul) << Shift;
      Shift += initialStateDownMinus[i];
      ++Shift;
    }
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateKyMax= maximum lz value reached by the fermionic state
// finalState = reference on the array where the bosonic state for the type up particles has to be stored

inline void BosonOnTorusWithSU4SpinAndMagneticTranslations::FermionToBoson(unsigned long initialState, int initialStateKyMax, unsigned long*& finalState)
{
  int FinalStateKyMax = 0;
  while ((initialStateKyMax >= 0) && ((initialState >> initialStateKyMax) == 0x0ul))
    --initialStateKyMax;
  while (initialStateKyMax >= 0)
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
      finalState[FinalStateKyMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState >>= TmpPower;
      ++FinalStateKyMax;
      initialStateKyMax -= TmpPower;
    }
  for (; FinalStateKyMax <= this->KyMax; ++FinalStateKyMax)
    finalState[FinalStateKyMax] = 0x0ul;
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

inline void BosonOnTorusWithSU4SpinAndMagneticTranslations::FermionToBoson(unsigned long initialStateUpPlus, unsigned long initialStateUpMinus, 
									   unsigned long initialStateDownPlus, unsigned long initialStateDownMinus,
									   unsigned long*& finalStateUpPlus, unsigned long*& finalStateUpMinus, 
									   unsigned long*& finalStateDownPlus, unsigned long*& finalStateDownMinus)
{
  int FinalStateKyMax = 0;
  int InitialStateKyMax = this->NUpPlusKyMax;
  while ((InitialStateKyMax >= 0) && ((initialStateUpPlus >> InitialStateKyMax) == 0x0ul))
    --InitialStateKyMax;
  while (InitialStateKyMax >= 0)
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
      finalStateUpPlus[FinalStateKyMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateUpPlus >>= TmpPower;
      ++FinalStateKyMax;
      InitialStateKyMax -= TmpPower;
    }
  for (; FinalStateKyMax <= this->KyMax; ++FinalStateKyMax)
    finalStateUpPlus[FinalStateKyMax] = 0x0ul;

  FinalStateKyMax = 0;
  InitialStateKyMax = this->NUpMinusKyMax;
  while ((InitialStateKyMax >= 0) && ((initialStateUpMinus >> InitialStateKyMax) == 0x0ul))
    --InitialStateKyMax;
  while (InitialStateKyMax >= 0)
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
      finalStateUpMinus[FinalStateKyMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateUpMinus >>= TmpPower;
      ++FinalStateKyMax;
      InitialStateKyMax -= TmpPower;
    }
  for (; FinalStateKyMax <= this->KyMax; ++FinalStateKyMax)
    finalStateUpMinus[FinalStateKyMax] = 0x0ul;

  FinalStateKyMax = 0;
  InitialStateKyMax = this->NDownPlusKyMax;
  while ((InitialStateKyMax >= 0) && ((initialStateDownPlus >> InitialStateKyMax) == 0x0ul))
    --InitialStateKyMax;
  while (InitialStateKyMax >= 0)
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
      finalStateDownPlus[FinalStateKyMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateDownPlus >>= TmpPower;
      ++FinalStateKyMax;
      InitialStateKyMax -= TmpPower;
    }
  for (; FinalStateKyMax <= this->KyMax; ++FinalStateKyMax)
    finalStateDownPlus[FinalStateKyMax] = 0x0ul;

  FinalStateKyMax = 0;
  InitialStateKyMax = this->NDownMinusKyMax;
  while ((InitialStateKyMax >= 0) && ((initialStateDownMinus >> InitialStateKyMax) == 0x0ul))
    --InitialStateKyMax;
  while (InitialStateKyMax >= 0)
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
      finalStateDownMinus[FinalStateKyMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateDownMinus >>= TmpPower;
      ++FinalStateKyMax;
      InitialStateKyMax -= TmpPower;
    }
  for (; FinalStateKyMax <= this->KyMax; ++FinalStateKyMax)
    finalStateDownMinus[FinalStateKyMax] = 0x0ul;
}

#endif
