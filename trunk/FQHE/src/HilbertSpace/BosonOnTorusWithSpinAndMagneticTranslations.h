////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of boson with spin on a torus taking                //
//                    into account magnetic translations                      //
//                                                                            //
//                        last modification : 26/04/2012                      //
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


#ifndef BOSONONTORUSWITHSPINANDMAGNETICTRANSLATIONS_H
#define BOSONONTORUSWITHSPINANDMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"

#include <cmath>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class BosonOnTorusWithSpinAndMagneticTranslations :  public ParticleOnTorusWithSpinAndMagneticTranslations
{

  friend class BosonOnTorusWithSpinAllSzAndMagneticTranslations;

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

  // number of bosons with spin up
  int NbrBosonsUp;
  // number of bosons with spin down
  int NbrBosonsDown;
  // total value of Spin
  int TotalSpin;
  // flag to indicate that the Hilbert space should preserve Sz
  bool SzFlag;

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
  // mask that corresponds to last bit that can be set to one for spin up
  unsigned long LastMomentumMaskUp;
  // mask that corresponds to last bit that can be set to one for spin down
  unsigned long LastMomentumMaskDown;
  // mask corresponding to StateShift
  unsigned long MomentumMask;

  // lzmax value for the fermionic states associated to the type up particles
  int NUpKyMax;
  // lzmax value for the fermionic states associated to the type down particles
  int NDownKyMax;
  // largest lzmax value for the fermionic states among N1KyMax, N2KyMax and N3KyMax
  int FermionicKyMax;

  // temporay array describing the type up particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionUp;
  // temporay array describing the type down particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionDown;

  // sorted array that contains each unique configuration for the type up particles
  unsigned long* UniqueStateDescriptionUp;
  // number of time each unique configuration for the type up particles appears in StateDescriptionUp
  int* UniqueStateDescriptionSubArraySizeUp;
  // number of unique configurations for the type up-plus particles
  long NbrUniqueStateDescriptionUp;
  // first time a type up appears in the Hilbert space
  int* FirstIndexUniqueStateDescriptionUp;

  // temporary state used when applying operators for type up particles
  unsigned long* TemporaryStateUp;
  // temporary state used when applying operators for type down particles
  unsigned long* TemporaryStateDown;

  // temporary state used when applying ProdA operator for type up particles
  unsigned long* ProdATemporaryStateUp;
  // temporary state used when applying ProdA operator for type down particles
  unsigned long* ProdATemporaryStateDown;
  // temporay pointer to the rescaling factor of the state currently processed by ProdA
  double* ProdANbrStateInOrbit;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

  // target space for operations leaving the Hilbert-space
  BosonOnTorusWithSpinAndMagneticTranslations* TargetSpace;

 public:

  
  // default constructor
  // 
  BosonOnTorusWithSpinAndMagneticTranslations ();

  // basic constructor
  // 
  // nbrBosons= number of bosons
  // totalSpin = twice the total spin value
  // maxMomentum = momentum maximum value for a boson
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)
  BosonOnTorusWithSpinAndMagneticTranslations (int nbrBosons, int totalSpin, int maxMomentum, int xMomentum, int yMomentum);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnTorusWithSpinAndMagneticTranslations(const BosonOnTorusWithSpinAndMagneticTranslations& bosons);

  // destructor
  //
  ~BosonOnTorusWithSpinAndMagneticTranslations();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnTorusWithSpinAndMagneticTranslations& operator = (const BosonOnTorusWithSpinAndMagneticTranslations& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

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

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphereWithSpin* targetSpace);
  
  
  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

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

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // return value =  multiplicative factor 
  virtual double AuAu (int index, int n1, int n2);

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

  // apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdu (int m1, int m2, double& coefficient, int& nbrTranslation);

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

  // apply a^+_m_d a_m_u operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, double& coefficient, int& nbrTranslation);

  // apply a^+_m_u a_m_d operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AduAd (int index, int m, double& coefficient, int& nbrTranslation);

  // apply a^+_m_d a_n_u operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation/annihilation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, int n, double& coefficient, int& nbrTranslation);

  // apply a^+_m_u a_n_d operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation/annihilation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AduAd (int index, int m, int n, double& coefficient, int& nbrTranslation);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int* spinIndices, int nbrIndices);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient, int& nbrTranslation);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis  
  virtual ComplexVector ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space);

 protected:

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initial bosonic state is stored
  // return value = corresponding fermionic state
  virtual unsigned long BosonToFermion(unsigned long*& initialState);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialStateUp = reference on the array where initial bosonic state for the type up particles is stored
  // initialStateDown = reference on the array where initial bosonic state for the type down particles is stored
  // finalStateUp = reference on the corresponding fermionic state for the type up particles
  // finalStateDown = reference on the corresponding fermionic state for the type down particles
  virtual void BosonToFermion(unsigned long*& initialStateUp, unsigned long*& initialStateDown, 
			      unsigned long& finalStateUp, unsigned long& finalStateDown);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateKyMax= maximum lz value reached by the fermionic state
  // finalState = reference on the array where the bosonic state for the type up particles has to be stored
  virtual void FermionToBoson(unsigned long initialState, int initialStateKyMax, unsigned long*& finalState);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialStateUp = initial fermionic state for the type up particles
  // initialStateDown = initial fermionic state for the type down particles
  // finalStateUp = reference on the array where the bosonic state for the type up particles has to be stored
  // finalStateDown = reference on the array where the bosonic state for the type down particles has to be stored
  virtual void FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown, 
			      unsigned long*& finalStateUp, unsigned long*& finalStateDown);

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
  // stateDescriptionUp = reference on the unsigned integer describing the state for spin up
  // stateDescriptionDown = reference on the unsigned integer describing the state for spin down
  // nbrTranslation = number of translation needed to obtain the canonical form
  virtual void FindCanonicalForm(unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown, int& nbrTranslation);
  
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescriptionUp = reference on the unsigned integer describing the state for spin up
  // stateDescriptionDown = reference on the unsigned integer describing the state for spin down
  // nbrTranslation = number of translation needed to obtain the canonical form
  // return value = true if the state does not fit the x momentum constraint
  virtual bool FindCanonicalFormAndTestXMomentumConstraint(unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown, int& nbrTranslation);
  
  // find how many translations on the x direction are needed to obtain the same state
  //
  // stateDescriptionUp = unsigned integer describing the state for spin up
  // stateDescriptionDown = unsigned integer describing the state for spin down
  // return value = number of translation needed to obtain the same state
  virtual int FindNumberXTranslation(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown);
  
  // test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
  //
  // stateDescriptionUp = unsigned integer describing the state for spin up
  // stateDescriptionDown = unsigned integer describing the state for spin down
  // return value = true if the state satisfy the x momentum constraint
  virtual bool TestXMomentumConstraint(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown);
				       
  // apply a single translation to a bosonic state in its fermionic representation
  //
  // stateDescriptionUp = reference on the unsigned integer describing the state for spin up
  // stateDescriptionDown = reference on the unsigned integer describing the state for spin down
  virtual void ApplySingleTranslation(unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown);

  // find state index
  //
  // stateDescriptionUp = unsigned integer describing the fermionic state for type up particles
  // stateDescriptionDown = unsigned integer describing the fermionic state for type down particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // return value = hilbert space dimension
  virtual long GenerateStates();
 
  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy);

  // evaluate Hilbert space dimension for a given total spin momentum
  //
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // nbrSpinUp = number of particles with spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrSpinUp);

  // generate all states corresponding to the constraints without the mangetic translations
  // 
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // currentFermionicPositionUp = current fermionic position within the state description for the spin up
  // currentFermionicPositionDown = current fermionic position within the state description for the spin down
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrBosons, int currentKy, int currentTotalKy, int currentFermionicPositionUp, int currentFermionicPositionDown, long pos);

  // generate all states corresponding to the constraints without the mangetic translations
  // 
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // currentFermionicPositionUp = current fermionic position within the state description for the spin up
  // currentFermionicPositionDown = current fermionic position within the state description for the spin down
  // nbrSpinUp = number of particles with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrBosons, int currentKy, int currentTotalKy, 
				 int currentFermionicPositionUp, int currentFermionicPositionDown, int nbrSpinUp, long pos);

};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnTorusWithSpinAndMagneticTranslations::GetParticleStatistic()
{
  return ParticleOnTorusWithSpinAndMagneticTranslations::BosonicStatistic;
}

// find how many translations on the x direction are needed to obtain the same state
//
// stateDescriptionUp = unsigned integer describing the state for spin up
// stateDescriptionDown = unsigned integer describing the state for spin down
// return value = number of translation needed to obtain the same state

inline int BosonOnTorusWithSpinAndMagneticTranslations::FindNumberXTranslation(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
{
  unsigned long TmpStateUp = stateDescriptionUp;
  unsigned long TmpStateDown = stateDescriptionDown;
  this->ApplySingleTranslation(TmpStateUp, TmpStateDown);
  int Index = 1;  
  while ((TmpStateUp != stateDescriptionUp) || (TmpStateDown != stateDescriptionDown))
    {
      this->ApplySingleTranslation(TmpStateUp, TmpStateDown);
      ++Index;  
    }
  return Index;
}

// find canonical form of a state description
//
// stateDescriptionUp = reference on the unsigned integer describing the state for spin up
// stateDescriptionDown = reference on the unsigned integer describing the state for spin down
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description (fermionic representation)

inline void BosonOnTorusWithSpinAndMagneticTranslations::FindCanonicalForm(unsigned long& stateDescriptionUp, 
										unsigned long& stateDescriptionDown, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalStateUp = stateDescriptionUp;
  unsigned long CanonicalStateDown = stateDescriptionDown;
  for (int i = 0; i < this->MomentumModulo; ++i)
    {
      this->ApplySingleTranslation(stateDescriptionUp, stateDescriptionDown);
      if ((stateDescriptionUp < CanonicalStateUp) || 
	  ((stateDescriptionUp == CanonicalStateUp) && (stateDescriptionDown < CanonicalStateDown)))
	{
	  CanonicalStateUp = stateDescriptionUp;
	  CanonicalStateDown = stateDescriptionDown;
	  nbrTranslation = i;
	}
    }
}

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescriptionUp = reference on the unsigned integer describing the state for spin up
// stateDescriptionDown = reference on the unsigned integer describing the state for spin down
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = true if the state does not fit the x momentum constraint

inline bool BosonOnTorusWithSpinAndMagneticTranslations::FindCanonicalFormAndTestXMomentumConstraint(unsigned long& stateDescriptionUp,
												     unsigned long& stateDescriptionDown, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalStateUp = stateDescriptionUp;
  unsigned long InitialStateDescriptionUp = stateDescriptionUp;
  unsigned long CanonicalStateDown = stateDescriptionDown;
  unsigned long InitialStateDescriptionDown = stateDescriptionDown;
  int Index = 0;
  int OrbitSize = 0;
  while (Index < this->MomentumModulo)
    {
      this->ApplySingleTranslation(stateDescriptionUp, stateDescriptionDown);
      ++Index;
      ++OrbitSize;
      if ((stateDescriptionUp != InitialStateDescriptionUp) || (stateDescriptionDown != InitialStateDescriptionDown))
	{
	  if ((stateDescriptionUp < CanonicalStateUp) ||
	      ((stateDescriptionUp == CanonicalStateUp) && (stateDescriptionDown < CanonicalStateDown)))
	    {
	      CanonicalStateUp = stateDescriptionUp;
	      CanonicalStateDown = stateDescriptionDown;
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
  stateDescriptionUp = CanonicalStateUp;
  stateDescriptionDown = CanonicalStateDown;
  return true;
}


// apply a single translation to a bosonic state in its fermionic representation
//
// stateDescriptionUp = reference on the unsigned integer describing the state for spin up
// stateDescriptionDown = reference on the unsigned integer describing the state for spin down

inline void BosonOnTorusWithSpinAndMagneticTranslations::ApplySingleTranslation(unsigned long& stateDescriptionUp,
										unsigned long& stateDescriptionDown)
{
  for (int i = 0; i < this->StateShift;)
    {
      while ((i < this->StateShift) && ((stateDescriptionUp & 0x1ul) == 0x0ul))
	{
	  stateDescriptionUp >>= 1;
	  ++i;
	}
      if (i < this->StateShift)
	{
	  while ((stateDescriptionUp & 0x1ul) == 0x1ul)
	    {
	      stateDescriptionUp >>= 1;
	      stateDescriptionUp |= this->LastMomentumMaskUp;
	    }
	  stateDescriptionUp >>= 1;	  
	  ++i;
	}
    }
  for (int i = 0; i < this->StateShift;)
    {
      while ((i < this->StateShift) && ((stateDescriptionDown & 0x1ul) == 0x0ul))
	{
	  stateDescriptionDown >>= 1;
	  ++i;
	}
      if (i < this->StateShift)
	{
	  while ((stateDescriptionDown & 0x1ul) == 0x1ul)
	    {
	      stateDescriptionDown >>= 1;
	      stateDescriptionDown |= this->LastMomentumMaskDown;
	    }
	  stateDescriptionDown >>= 1;	  
	  ++i;
	}
    }
}

// test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescriptionUp = unsigned integer describing the state for spin up
// stateDescriptionDown = unsigned integer describing the state for spin down
// return value = true if the state satisfy the x momentum constraint

inline bool BosonOnTorusWithSpinAndMagneticTranslations::TestXMomentumConstraint(unsigned long stateDescriptionUp,
										      unsigned long stateDescriptionDown)
{
  if (((this->KxMomentum * this->FindNumberXTranslation(stateDescriptionUp, stateDescriptionDown)) % this->MomentumModulo) == 0)
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

inline int BosonOnTorusWithSpinAndMagneticTranslations::AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej, double& coefficient, int& nbrTranslation)
{
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      this->TemporaryStateUp[i] = this->ProdATemporaryStateUp[i];
      this->TemporaryStateDown[i] = this->ProdATemporaryStateDown[i];
    }
  ++temporaryStatej[m2];
  coefficient = temporaryStatej[m2];
  ++temporaryStatei[m1];
  coefficient *= temporaryStatei[m1];
  coefficient = sqrt(coefficient);
  unsigned long TmpStateUp;
  unsigned long TmpStateDown;
  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpStateUp, TmpStateDown);
  if (this->FindCanonicalFormAndTestXMomentumConstraint(TmpStateUp, TmpStateDown, nbrTranslation) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpStateUp, TmpStateDown);
  coefficient *= this->ProdANbrStateInOrbit[this->NbrStateInOrbit[TmpIndex]];
  nbrTranslation *= this->StateShift;
  return TmpIndex;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initial bosonic state is stored
// return value = corresponding fermionic state

inline unsigned long BosonOnTorusWithSpinAndMagneticTranslations::BosonToFermion(unsigned long*& initialState)
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
// initialStateUp = reference on the array where initial bosonic state for the type up particles is stored
// initialStateDown = reference on the array where initial bosonic state for the type down particles is stored
// finalStateUp = reference on the corresponding fermionic state for the type up particles
// finalStateDown = reference on the corresponding fermionic state for the type down particles

inline void BosonOnTorusWithSpinAndMagneticTranslations::BosonToFermion(unsigned long*& initialStateUp, unsigned long*& initialStateDown,
									unsigned long& finalStateUp, unsigned long& finalStateDown)
{
  finalStateUp = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->KyMax; ++i)
    {
      finalStateUp |= ((1ul << initialStateUp[i]) - 1ul) << Shift;
      Shift += initialStateUp[i];
      ++Shift;
    }
  finalStateDown = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->KyMax; ++i)
    {
      finalStateDown |= ((1ul << initialStateDown[i]) - 1ul) << Shift;
      Shift += initialStateDown[i];
      ++Shift;
    }
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateKyMax= maximum lz value reached by the fermionic state
// finalState = reference on the array where the bosonic state for the type up particles has to be stored

inline void BosonOnTorusWithSpinAndMagneticTranslations::FermionToBoson(unsigned long initialState, int initialStateKyMax, unsigned long*& finalState)
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
// initialStateUp = initial fermionic state for the type up particles
// initialStateDown = initial fermionic state for the type down particles
// finalStateUp = reference on the array where the bosonic state for the type up particles has to be stored
// finalStateDown = reference on the array where the bosonic state for the type down particles has to be stored

inline void BosonOnTorusWithSpinAndMagneticTranslations::FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown,
									unsigned long*& finalStateUp, unsigned long*& finalStateDown)
{
  int FinalStateKyMax = 0;
  int InitialStateKyMax = this->NUpKyMax;
  while ((InitialStateKyMax >= 0) && ((initialStateUp >> InitialStateKyMax) == 0x0ul))
    --InitialStateKyMax;
  while (InitialStateKyMax >= 0)
    {
      unsigned long TmpState = (~initialStateUp - 1ul) ^ (~initialStateUp);
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
      finalStateUp[FinalStateKyMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateUp >>= TmpPower;
      ++FinalStateKyMax;
      InitialStateKyMax -= TmpPower;
    }
  for (; FinalStateKyMax <= this->KyMax; ++FinalStateKyMax)
    finalStateUp[FinalStateKyMax] = 0x0ul;

  FinalStateKyMax = 0;
  InitialStateKyMax = this->NDownKyMax;
  while ((InitialStateKyMax >= 0) && ((initialStateDown >> InitialStateKyMax) == 0x0ul))
    --InitialStateKyMax;
  while (InitialStateKyMax >= 0)
    {
      unsigned long TmpState = (~initialStateDown - 1ul) ^ (~initialStateDown);
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
      finalStateDown[FinalStateKyMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateDown >>= TmpPower;
      ++FinalStateKyMax;
      InitialStateKyMax -= TmpPower;
    }
  for (; FinalStateKyMax <= this->KyMax; ++FinalStateKyMax)
    finalStateDown[FinalStateKyMax] = 0x0ul;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

inline int BosonOnTorusWithSpinAndMagneticTranslations::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// get the number of orbitals
//
// return value = number of orbitals

inline int BosonOnTorusWithSpinAndMagneticTranslations::GetNbrOrbitals()
{
  return this->MaxMomentum;
}

// get the number of particles
//
// return value = number of particles

inline int BosonOnTorusWithSpinAndMagneticTranslations::GetNbrParticles()
{
  return this->NbrBosons;
}

// get the total spin
//
//return value: total spin of the Hilbert space

inline int BosonOnTorusWithSpinAndMagneticTranslations::GetTotalSpin()
{
  return this->TotalSpin;
}

// get the momentum along the x axis
// 
// return avlue = momentum along the x axis

inline int BosonOnTorusWithSpinAndMagneticTranslations::GetKxMomentum()
{
  return this->KxMomentum;
}

// get the momentum along the y axis
// 
// return avlue = momentum along the y axis

inline int BosonOnTorusWithSpinAndMagneticTranslations::GetKyMomentum()
{
  return this->KyMomentum;
}

// get the maximum momentum along the x axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the x axis

inline int BosonOnTorusWithSpinAndMagneticTranslations::GetMaxXMomentum()
{
  return this->MomentumModulo;
}

// get the maximum momentum along the y axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the y axis

inline int BosonOnTorusWithSpinAndMagneticTranslations::GetMaxYMomentum()
{
  return this->MaxMomentum;
}
  
  
#endif
