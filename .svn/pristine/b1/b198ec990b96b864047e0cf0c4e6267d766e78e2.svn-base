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


#ifndef BOSONONTORUSWITHSPINALLSZANDMAGNETICTRANSLATIONS_H
#define BOSONONTORUSWITHSPINALLSZANDMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/BosonOnTorusWithSpinAndMagneticTranslations.h"

#include <cmath>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class BosonOnTorusWithSpinAllSzAndMagneticTranslations :  public BosonOnTorusWithSpinAndMagneticTranslations
{

 protected:

  
 public:

  
  // constructor when Sz is not conserved
  // 
  // nbrBosons= number of bosons
  // maxMomentum = momentum maximum value for a boson
  // kxMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // kyMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)
  BosonOnTorusWithSpinAllSzAndMagneticTranslations (int nbrBosons, int maxMomentum, int kxMomentum, int kyMomentum);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnTorusWithSpinAllSzAndMagneticTranslations(const BosonOnTorusWithSpinAllSzAndMagneticTranslations& bosons);

  // destructor
  //
  ~BosonOnTorusWithSpinAllSzAndMagneticTranslations();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnTorusWithSpinAllSzAndMagneticTranslations& operator = (const BosonOnTorusWithSpinAllSzAndMagneticTranslations& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphereWithSpin* targetSpace);
  
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

  // convert state of a SU(2) Hilbert space with fixed Sz to a SU(2) space with all sz sectors
  //
  // state = state that needs to be projected
  // su2space = SU(2) space with fixed sz of the input state
  // return value = input state expression in the SU(2) basis
  virtual ComplexVector SU2ToSU2AllSz(ComplexVector& state, ParticleOnSphereWithSpin* su2space);

  // convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // initialState = state to transform  
  // targetState = vector where the transformed state has to be stored
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // firstComponent = index of the first component to compute in initialState
  // nbrComponents = number of consecutive components to compute
  virtual void TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, long firstComponent = 0l, long nbrComponents = 0l);

 protected:

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
  // lastMomentumMaskUp = reference on themask that corresponds to last bit that can be set to one for spin up
  // lastMomentumMaskDown = mask that corresponds to last bit that can be set to one for spin down
  // return value = number of translation needed to obtain the same state
  virtual int FindNumberXTranslation(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown);
  
  // apply a single translation to a bosonic state in its fermionic representation
  //
  // stateDescriptionUp = reference on the unsigned integer describing the state for spin up
  // stateDescriptionDown = reference on the unsigned integer describing the state for spin down
  virtual void ApplySingleTranslation(unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown,
				      unsigned long& lastMomentumMaskUp, unsigned long& lastMomentumMaskDown);

  // generate all states corresponding to the constraints
  // 
  // return value = hilbert space dimension
  long GenerateStates();
 
  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy);

  // generate all states corresponding to the constraints without the mangetic translations
  // 
  // nbrBosons = number of bosons
  // currentKy = current momentum along y for a single particle
  // currentTotalKy = current total momentum along y
  // currentFermionicPositionUp = current fermionic position within the state description for the spin up
  // currentFermionicPositionDown = current fermionic position within the state description for the spin down
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long RawGenerateStates(int nbrBosons, int currentKy, int currentTotalKy, int currentFermionicPositionUp, int currentFermionicPositionDown, long pos);

  // compute the number of particles in a given state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of particles
  virtual int ComputeNbrParticles(unsigned long stateDescription);

  // recursive part of the convertion from a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // targetState = vector where the transformed state has to be stored
  // coefficient = current coefficient to assign
  // position = current particle consider in the n-body state
  // momentumIndices = array that gives the momentum partition of the initial n-body state
  // initialSU2Indices = array that gives the spin dressing the initial n-body state
  // currentSU2Indices = array that gives the spin dressing the current transformed n-body state
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
  // occupationCoefficientArray = array that provides 1/2 ln (N!)
  // fourrierCoefficients = array of Fourrier coefficients for the kx momentum
  virtual void TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
					      int position, int* momentumIndices, int* initialSU2Indices, int* currentSU2Indices, ComplexMatrix* oneBodyBasis, 
					      double occupationCoefficient, double* occupationCoefficientArray, Complex* fourrierCoefficients);

};

// find how many translations on the x direction are needed to obtain the same state
//
// stateDescriptionUp = unsigned integer describing the state for spin up
// stateDescriptionDown = unsigned integer describing the state for spin down
// return value = number of translation needed to obtain the same state

inline int BosonOnTorusWithSpinAllSzAndMagneticTranslations::FindNumberXTranslation(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
{
  unsigned long TmpStateUp = stateDescriptionUp;
  unsigned long TmpStateDown = stateDescriptionDown;
  int TmpNbrUps = this->ComputeNbrParticles(stateDescriptionUp);
  unsigned long TmpLastMomentumMaskUp = 0x1ul << (this->MaxMomentum + TmpNbrUps - 1);
  unsigned long TmpLastMomentumMaskDown = 0x1ul << (this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1);
  this->ApplySingleTranslation(TmpStateUp, TmpStateDown, TmpLastMomentumMaskUp, TmpLastMomentumMaskDown);
  int Index = 1;  
  while ((TmpStateUp != stateDescriptionUp) || (TmpStateDown != stateDescriptionDown))
    {
      this->ApplySingleTranslation(TmpStateUp, TmpStateDown, TmpLastMomentumMaskUp, TmpLastMomentumMaskDown);
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

inline void BosonOnTorusWithSpinAllSzAndMagneticTranslations::FindCanonicalForm(unsigned long& stateDescriptionUp, 
										unsigned long& stateDescriptionDown, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalStateUp = stateDescriptionUp;
  unsigned long CanonicalStateDown = stateDescriptionDown;
  int TmpNbrUps = this->ComputeNbrParticles(stateDescriptionUp);
  unsigned long TmpLastMomentumMaskUp = 0x1ul << (this->MaxMomentum + TmpNbrUps - 1);
  unsigned long TmpLastMomentumMaskDown = 0x1ul << (this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1);
  for (int i = 0; i < this->MomentumModulo; ++i)
    {
      this->ApplySingleTranslation(stateDescriptionUp, stateDescriptionDown, TmpLastMomentumMaskUp, TmpLastMomentumMaskDown);
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

inline bool BosonOnTorusWithSpinAllSzAndMagneticTranslations::FindCanonicalFormAndTestXMomentumConstraint(unsigned long& stateDescriptionUp,
													  unsigned long& stateDescriptionDown, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalStateUp = stateDescriptionUp;
  unsigned long InitialStateDescriptionUp = stateDescriptionUp;
  unsigned long CanonicalStateDown = stateDescriptionDown;
  unsigned long InitialStateDescriptionDown = stateDescriptionDown;
  int TmpNbrUps = this->ComputeNbrParticles(stateDescriptionUp);
  unsigned long TmpLastMomentumMaskUp = 0x1ul << (this->MaxMomentum + TmpNbrUps - 1);
  unsigned long TmpLastMomentumMaskDown = 0x1ul << (this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1);
  int Index = 0;
  int OrbitSize = 0;
  while (Index < this->MomentumModulo)
    {
      this->ApplySingleTranslation(stateDescriptionUp, stateDescriptionDown, TmpLastMomentumMaskUp, TmpLastMomentumMaskDown);
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
// lastMomentumMaskUp = reference on themask that corresponds to last bit that can be set to one for spin up
// lastMomentumMaskDown = mask that corresponds to last bit that can be set to one for spin down

inline void BosonOnTorusWithSpinAllSzAndMagneticTranslations::ApplySingleTranslation(unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown,
										     unsigned long& lastMomentumMaskUp, unsigned long& lastMomentumMaskDown)
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
	      stateDescriptionUp |= lastMomentumMaskUp;
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
	      stateDescriptionDown |= lastMomentumMaskDown;
	    }
	  stateDescriptionDown >>= 1;	  
	  ++i;
	}
    }
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateKyMax= maximum lz value reached by the fermionic state
// finalState = reference on the array where the bosonic state for the type up particles has to be stored

inline void BosonOnTorusWithSpinAllSzAndMagneticTranslations::FermionToBoson(unsigned long initialState, int initialStateKyMax, unsigned long*& finalState)
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

inline void BosonOnTorusWithSpinAllSzAndMagneticTranslations::FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown,
									     unsigned long*& finalStateUp, unsigned long*& finalStateDown)
{
  int FinalStateKyMax = 0;
  int TmpNbrUps = this->ComputeNbrParticles(initialStateUp);
  int InitialStateKyMax = this->MaxMomentum + TmpNbrUps - 1;
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
  InitialStateKyMax = this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1;
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

// compute the number of particles in a given state
//
// stateDescription = unsigned integer describing the state
// return value = number of particles

inline int BosonOnTorusWithSpinAllSzAndMagneticTranslations::ComputeNbrParticles(unsigned long stateDescription)
{
  int TmpNbrParticle = 0;
  while(stateDescription) 
    {
      ++TmpNbrParticle;
      stateDescription &= stateDescription - 0x1ul;
    }
  return TmpNbrParticle;
}


  
#endif
