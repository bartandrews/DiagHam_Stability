////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of boson with SU3 spin on a torus taking               //
//                    into account magnetic translations                      //
//                                                                            //
//                        last modification : 18/06/2012                      //
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


#ifndef BOSONONTORUSWITHSU3SPINANDMAGNETICTRANSLATIONS_H
#define BOSONONTORUSWITHSU3SPINANDMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSU3SpinAndMagneticTranslations.h"

#include <cmath>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class BosonOnTorusWithSU3SpinAndMagneticTranslations :  public ParticleOnTorusWithSU3SpinAndMagneticTranslations
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

  // number of bosons of type 1
  int NbrBosons1;
  // number of bosons of type 2
  int NbrBosons2;
  // number of bosons of type 3
  int NbrBosons3;
  // twice the total Tz value
  int TotalTz;
  // three time the total Y value
  int TotalY;

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
  // mask that corresponds to last bit that can be set to one for type 1 particles 
  unsigned long LastMomentumMask1;
   // mask that corresponds to last bit that can be set to one for type 2 particles 
  unsigned long LastMomentumMask2;
   // mask that corresponds to last bit that can be set to one for type 3 particles 
  unsigned long LastMomentumMask3;
 // mask corresponding to StateShift
  unsigned long MomentumMask;

  // kymax value for the fermionic states associated to the type 1 particles
  int N1KyMax;
  // kymax value for the fermionic states associated to the type 2 particles
  int N2KyMax;
  // kymax value for the fermionic states associated to the type 3 particles
  int N3KyMax;
  // largest kymax value for the fermionic states among N1KyMax, N2KyMax and N3KyMax
  int FermionicKyMax;

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

  // temporary state used when applying ProdA operator for type 1 particles
  unsigned long* ProdATemporaryState1;
  // temporary state used when applying ProdA operator for type 2 particles
  unsigned long* ProdATemporaryState2;
  // temporary state used when applying ProdA operator for type 3 particles
  unsigned long* ProdATemporaryState3;
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
  // totalTz = twice the total Tz value
  // totalY = three time the total Y value
  // maxMomentum = momentum maximum value for a boson
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)
  BosonOnTorusWithSU3SpinAndMagneticTranslations (int nbrBosons, int totalTz, int totalY, int maxMomentum, int xMomentum, int yMomentum);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnTorusWithSU3SpinAndMagneticTranslations(const BosonOnTorusWithSU3SpinAndMagneticTranslations& bosons);

  // destructor
  //
  ~BosonOnTorusWithSU3SpinAndMagneticTranslations();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnTorusWithSU3SpinAndMagneticTranslations& operator = (const BosonOnTorusWithSU3SpinAndMagneticTranslations& bosons);

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
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad1Ad1 (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_1 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad1Ad2 (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_1 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad1Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_2 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad2Ad2 (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_2 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad2Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation);

  // apply a^+_m1_3 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad3Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation);

 protected:

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
  // initialStateKyMax= maximum lz value reached by the fermionic state
  // finalState = reference on the array where the bosonic state for the type up particles has to be stored
  virtual void FermionToBoson(unsigned long initialState, int initialStateKyMax, unsigned long*& finalState);

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
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = index of the destination state
  virtual int AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej, double& coefficient, int& nbrTranslation);

  // find canonical form of a state description
  //
  // stateDescription1 = reference on the unsigned integer describing the state related to type 1 particles
  // stateDescription2 = reference on the unsigned integer describing the state related to type 2 particles
  // stateDescription3 = reference on the unsigned integer describing the state related to type 3 particles
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = canonical form of a state description (fermionic representation)
  virtual void FindCanonicalForm(unsigned long& stateDescription1, unsigned long& stateDescription2, unsigned long& stateDescription3, int& nbrTranslation);
  
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescription1 = reference on the unsigned integer describing the state related to type 1 particles
  // stateDescription2 = reference on the unsigned integer describing the state related to type 2 particles
  // stateDescription3 = reference on the unsigned integer describing the state related to type 3 particles
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = true if the state does not fit the x momentum constraint
  virtual bool FindCanonicalFormAndTestXMomentumConstraint(unsigned long& stateDescription1, unsigned long& stateDescription2, 
							   unsigned long& stateDescription3, int& nbrTranslation);
  
  // find how many translations on the x direction are needed to obtain the same state
  //
  // stateDescription1 = unsigned integer describing the state for type 1 particles
  // stateDescription2 = unsigned integer describing the state for type 2 particles
  // stateDescription3 = unsigned integer describing the state for type 3 particles
  // return value = number of translation needed to obtain the same state
  virtual int FindNumberXTranslation(unsigned long stateDescription1, unsigned long stateDescription2, unsigned long stateDescription3);
  
  // test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescription1 = unsigned integer describing the state related to type 1 particles
  // stateDescription2 = unsigned integer describing the state related to type 2 particles
  // stateDescription3 = unsigned integer describing the state related to type 3 particles
  // return value = true if the state satisfy the x momentum constraint  
  virtual bool TestXMomentumConstraint(unsigned long stateDescription1, unsigned long stateDescription2, unsigned long stateDescription3);

  // apply a single translation to a bosonic state in its fermionic representation
  //
  // stateDescription1 = reference on the unsigned integer describing the state related to type 1 particles
  // stateDescription2 = reference on the unsigned integer describing the state related to type 2 particles
  // stateDescription3 = reference on the unsigned integer describing the state related to type 3 particles
  virtual void ApplySingleTranslation(unsigned long& stateDescription1, unsigned long& stateDescription2, unsigned long& stateDescription3);

  // find state index
  //
  // stateDescription1 = unsigned integer describing the fermionic state for type 1 particles
  // stateDescription2 = unsigned integer describing the fermionic state for type 2 particles
  // stateDescription3 = unsigned integer describing the fermionic state for type 3 particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription1, unsigned long stateDescription2, unsigned long stateDescription3);

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
  // nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
  // nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
  // nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrN1, int nbrN2, int nbrN3);

  // generate all states corresponding to the constraints without the mangetic translations
  // 
  // nbrBosons = number of bosons
  // currentKy1 = current momentum along y for a single type 1 particle
  // currentKy2 = current momentum along y for a single type 2 particle
  // currentKy3 = current momentum along y for a single type 3 particle
  // currentTotalKy = current total momentum along y
  // nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
  // nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
  // nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long RawGenerateStates(int nbrBosons, int currentKy1, int currentKy2, int currentKy3, int currentTotalKy, int nbrN1, int nbrN2, int nbrN3, long pos);

};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnTorusWithSU3SpinAndMagneticTranslations::GetParticleStatistic()
{
  return ParticleOnTorusWithSU3SpinAndMagneticTranslations::BosonicStatistic;
}

// find how many translations on the x direction are needed to obtain the same state
//
// stateDescription1 = unsigned integer describing the state for type 1 particles
// stateDescription2 = unsigned integer describing the state for type 2 particles
// stateDescription3 = unsigned integer describing the state for type 3 particles
// return value = number of translation needed to obtain the same state

inline int BosonOnTorusWithSU3SpinAndMagneticTranslations::FindNumberXTranslation(unsigned long stateDescription1, unsigned long stateDescription2, unsigned long stateDescription3)
{
  unsigned long TmpState1 = stateDescription1;
  unsigned long TmpState2 = stateDescription2;
  unsigned long TmpState3 = stateDescription3;
  this->ApplySingleTranslation(TmpState1, TmpState2, TmpState3);
  int Index = 1;  
  while ((TmpState1 != stateDescription1) || (TmpState2 != stateDescription2) || (TmpState3 != stateDescription3))
    {
      this->ApplySingleTranslation(TmpState1, TmpState2, TmpState3);
      ++Index;  
    }
  return Index;
}

// find canonical form of a state description
//
// stateDescription1 = reference on the unsigned integer describing the state related to type 1 particles
// stateDescription2 = reference on the unsigned integer describing the state related to type 2 particles
// stateDescription3 = reference on the unsigned integer describing the state related to type 3 particles
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description (fermionic representation)

inline void BosonOnTorusWithSU3SpinAndMagneticTranslations::FindCanonicalForm(unsigned long& stateDescription1, unsigned long& stateDescription2, unsigned long& stateDescription3,
									      int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState1 = stateDescription1;
  unsigned long CanonicalState2 = stateDescription2;
  unsigned long CanonicalState3 = stateDescription3;
  for (int i = 0; i < this->MomentumModulo; ++i)
    {
      this->ApplySingleTranslation(stateDescription1, stateDescription2, stateDescription3);
      if ((stateDescription1 < CanonicalState1) || 
	  ((stateDescription1 == CanonicalState1) && ((stateDescription2 < CanonicalState2) || ((stateDescription2 == CanonicalState2) && (stateDescription3 < CanonicalState3)))))
	{
	  CanonicalState1 = stateDescription1;
	  CanonicalState2 = stateDescription2;
	  CanonicalState3 = stateDescription3;
	  nbrTranslation = i;
	}
    }
}

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescription1 = reference on the unsigned integer describing the state related to type 1 particles
// stateDescription2 = reference on the unsigned integer describing the state related to type 2 particles
// stateDescription3 = reference on the unsigned integer describing the state related to type 3 particles
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = true if the state does not fit the x momentum constraint

inline bool BosonOnTorusWithSU3SpinAndMagneticTranslations::FindCanonicalFormAndTestXMomentumConstraint(unsigned long& stateDescription1, unsigned long& stateDescription2, 
													unsigned long& stateDescription3, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState1 = stateDescription1;
  unsigned long InitialStateDescription1 = stateDescription1;
  unsigned long CanonicalState2 = stateDescription2;
  unsigned long InitialStateDescription2 = stateDescription2;
  unsigned long CanonicalState3 = stateDescription3;
  unsigned long InitialStateDescription3 = stateDescription3;
  int Index = 0;
  int OrbitSize = 0;
  while (Index < this->MomentumModulo)
    {
      this->ApplySingleTranslation(stateDescription1, stateDescription2, stateDescription3);
      ++Index;
      ++OrbitSize;
      if ((stateDescription1 != InitialStateDescription1) || (stateDescription2 != InitialStateDescription2) || (stateDescription3 != InitialStateDescription3))
	{
	  if ((stateDescription1 < CanonicalState1) || 
	      ((stateDescription1 == CanonicalState1) && ((stateDescription2 < CanonicalState2) || ((stateDescription2 == CanonicalState2) && (stateDescription3 < CanonicalState3)))))
	    {
	      CanonicalState1 = stateDescription1;
	      CanonicalState2 = stateDescription2;
	      CanonicalState3 = stateDescription3;
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
  stateDescription1 = CanonicalState1;
  stateDescription2 = CanonicalState2;
  stateDescription3 = CanonicalState3;
  return true;
}


// apply a single translation to a bosonic state in its fermionic representation
//
// stateDescription1 = reference on the unsigned integer describing the state related to type 1 particles
// stateDescription2 = reference on the unsigned integer describing the state related to type 2 particles
// stateDescription3 = reference on the unsigned integer describing the state related to type 3 particles

inline void BosonOnTorusWithSU3SpinAndMagneticTranslations::ApplySingleTranslation(unsigned long& stateDescription1, unsigned long& stateDescription2,
										   unsigned long& stateDescription3)
{
  for (int i = 0; i < this->StateShift;)
    {
      while ((i < this->StateShift) && ((stateDescription1 & 0x1ul) == 0x0ul))
	{
	  stateDescription1 >>= 1;
	  ++i;
	}
      if (i < this->StateShift)
	{
	  while ((stateDescription1 & 0x1ul) == 0x1ul)
	    {
	      stateDescription1 >>= 1;
	      stateDescription1 |= this->LastMomentumMask1;
	    }
	  stateDescription1 >>= 1;	  
	  ++i;
	}
    }
  for (int i = 0; i < this->StateShift;)
    {
      while ((i < this->StateShift) && ((stateDescription2 & 0x1ul) == 0x0ul))
	{
	  stateDescription2 >>= 1;
	  ++i;
	}
      if (i < this->StateShift)
	{
	  while ((stateDescription2 & 0x1ul) == 0x1ul)
	    {
	      stateDescription2 >>= 1;
	      stateDescription2 |= this->LastMomentumMask2;
	    }
	  stateDescription2 >>= 1;	  
	  ++i;
	}
    }
  for (int i = 0; i < this->StateShift;)
    {
      while ((i < this->StateShift) && ((stateDescription3 & 0x1ul) == 0x0ul))
	{
	  stateDescription3>>= 1;
	  ++i;
	}
      if (i < this->StateShift)
	{
	  while ((stateDescription3 & 0x1ul) == 0x1ul)
	    {
	      stateDescription3 >>= 1;
	      stateDescription3 |= this->LastMomentumMask3;
	    }
	  stateDescription3 >>= 1;	  
	  ++i;
	}
    }
}

// test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescription1 = unsigned integer describing the state related to type 1 particles
// stateDescription2 = unsigned integer describing the state related to type 2 particles
// stateDescription3 = unsigned integer describing the state related to type 3 particles
// return value = true if the state satisfy the x momentum constraint

inline bool BosonOnTorusWithSU3SpinAndMagneticTranslations::TestXMomentumConstraint(unsigned long stateDescription1, unsigned long stateDescription2,
										    unsigned long stateDescription3)
{
  if (((this->KxMomentum * this->FindNumberXTranslation(stateDescription1, stateDescription2, stateDescription3)) % this->MomentumModulo) == 0)
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

inline int BosonOnTorusWithSU3SpinAndMagneticTranslations::AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej, double& coefficient, int& nbrTranslation)
{
  for (int i = 0; i < this->MaxMomentum; ++i)
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
  unsigned long TmpState1;
  unsigned long TmpState2;
  unsigned long TmpState3;
  this->BosonToFermion(this->TemporaryState1, this->TemporaryState2, this->TemporaryState3, TmpState1, TmpState2, TmpState3);
  if (this->FindCanonicalFormAndTestXMomentumConstraint(TmpState1, TmpState2, TmpState3, nbrTranslation) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState1, TmpState2, TmpState3);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->ProdANbrStateInOrbit[this->NbrStateInOrbit[TmpIndex]];
      nbrTranslation *= this->StateShift;
    }
  return TmpIndex;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initial bosonic state is stored
// return value = corresponding fermionic state

inline unsigned long BosonOnTorusWithSU3SpinAndMagneticTranslations::BosonToFermion(unsigned long*& initialState)
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
// initialState1 = reference on the array where initial bosonic state for the type 1 particles is stored
// initialState2 = reference on the array where initial bosonic state for the type 2 particles is stored
// initialState3 = reference on the array where initial bosonic state for the type 3 particles is stored
// finalState1 = reference on the corresponding fermionic state for the type 1 particles
// finalState2 = reference on the corresponding fermionic state for the type 2 particles
// finalState3 = reference on the corresponding fermionic state for the type 3 particles

inline void BosonOnTorusWithSU3SpinAndMagneticTranslations::BosonToFermion(unsigned long*& initialState1, unsigned long*& initialState2, unsigned long*& initialState3,
									   unsigned long& finalState1, unsigned long& finalState2, unsigned long& finalState3)
{
  finalState1 = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->KyMax; ++i)
    {
      finalState1 |= ((1ul << initialState1[i]) - 1ul) << Shift;
      Shift += initialState1[i];
      ++Shift;
    }
  finalState2 = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->KyMax; ++i)
    {
      finalState2 |= ((1ul << initialState2[i]) - 1ul) << Shift;
      Shift += initialState2[i];
      ++Shift;
    }
  finalState3 = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->KyMax; ++i)
    {
      finalState3 |= ((1ul << initialState3[i]) - 1ul) << Shift;
      Shift += initialState3[i];
      ++Shift;
    }
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateKyMax= maximum lz value reached by the fermionic state
// finalState = reference on the array where the bosonic state for the type up particles has to be stored

inline void BosonOnTorusWithSU3SpinAndMagneticTranslations::FermionToBoson(unsigned long initialState, int initialStateKyMax, unsigned long*& finalState)
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
// initialState1 = initial fermionic state for the type 1 particles
// initialState2 = initial fermionic state for the type 2 particles
// initialState3 = initial fermionic state for the type 3 particles
// finalState1 = reference on the array where the bosonic state for the type 1 particles has to be stored
// finalState2 = reference on the array where the bosonic state for the type 2 particles has to be stored
// finalState3 = reference on the array where the bosonic state for the type 3 particles has to be stored

inline void BosonOnTorusWithSU3SpinAndMagneticTranslations::FermionToBoson(unsigned long initialState1, unsigned long initialState2, unsigned long initialState3,
									   unsigned long*& finalState1, unsigned long*& finalState2, unsigned long*& finalState3)
{
  int FinalStateKyMax = 0;
  int InitialStateKyMax = this->N1KyMax;
  while ((InitialStateKyMax >= 0) && ((initialState1 >> InitialStateKyMax) == 0x0ul))
    --InitialStateKyMax;
  while (InitialStateKyMax >= 0)
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
      finalState1[FinalStateKyMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState1 >>= TmpPower;
      ++FinalStateKyMax;
      InitialStateKyMax -= TmpPower;
    }
  for (; FinalStateKyMax <= this->KyMax; ++FinalStateKyMax)
    finalState1[FinalStateKyMax] = 0x0ul;

  FinalStateKyMax = 0;
  InitialStateKyMax = this->N2KyMax;
  while ((InitialStateKyMax >= 0) && ((initialState2 >> InitialStateKyMax) == 0x0ul))
    --InitialStateKyMax;
  while (InitialStateKyMax >= 0)
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
      finalState2[FinalStateKyMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState2 >>= TmpPower;
      ++FinalStateKyMax;
      InitialStateKyMax -= TmpPower;
    }
  for (; FinalStateKyMax <= this->KyMax; ++FinalStateKyMax)
    finalState2[FinalStateKyMax] = 0x0ul;

  FinalStateKyMax = 0;
  InitialStateKyMax = this->N3KyMax;
  while ((InitialStateKyMax >= 0) && ((initialState3 >> InitialStateKyMax) == 0x0ul))
    --InitialStateKyMax;
  while (InitialStateKyMax >= 0)
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
      finalState3[FinalStateKyMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState3 >>= TmpPower;
      ++FinalStateKyMax;
      InitialStateKyMax -= TmpPower;
    }
  for (; FinalStateKyMax <= this->KyMax; ++FinalStateKyMax)
    finalState3[FinalStateKyMax] = 0x0ul;
}

#endif
