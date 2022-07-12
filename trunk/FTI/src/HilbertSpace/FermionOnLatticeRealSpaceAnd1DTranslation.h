////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of fermions on lattice                        //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 09/09/2014                      //
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


#ifndef FERMIONONLATTICEREALSPACEAND1DTRANSLATION_H
#define FERMIONONLATTICEREALSPACEAND1DTRANSLATION_H

#include "config.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"

#include <iostream>
#include <bitset>

using std::bitset;


class  FermionOnLatticeRealSpaceAnd1DTranslation : public FermionOnTorusWithMagneticTranslations
{

 protected:

  // total number of sites
  int NbrSite;
  // number of sites per unit cell
  int NbrSitePerUnitCell;
  
  // number of momentum sectors in the x direction 
  int MaxXMomentum;
  // bit shift that has to applied to perform a translation in the x direction 
  int StateXShift;
  // binary mask for the StateXShift first bits 
  unsigned long XMomentumMask;
  // bit shift to apply to move the first StateXShift bits at the end of a state description
  int ComplementaryStateXShift;

  // parity of the number of fermions, 0x1ul if even, 0x0ul if odd
  unsigned long NbrFermionsParity;

 public:

  // default constructor
  // 
  FermionOnLatticeRealSpaceAnd1DTranslation ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // memory = amount of memory granted for precalculations
  FermionOnLatticeRealSpaceAnd1DTranslation (int nbrFermions, int nbrSite, int xMomentum, int maxXMomentum,
					     unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeRealSpaceAnd1DTranslation(const FermionOnLatticeRealSpaceAnd1DTranslation& fermions);

  // destructor
  //
  ~FermionOnLatticeRealSpaceAnd1DTranslation ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeRealSpaceAnd1DTranslation& operator = (const FermionOnLatticeRealSpaceAnd1DTranslation& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // apply a^+_m_u a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient, int& nbrTranslationX);
  
  // apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient, int& nbrTranslationX);
  
  // apply a^+_m operator to the state produced using AuAu method (without destroying it)
  //
  // m = first index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int Ad (int m, double& coefficient, int& nbrTranslationX);

  // find state index from a string
  //
  // stateDescription = string describing the state
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(char* stateDescription);

  // get the momentum along the x axis
  // 
  // return avlue = momentum along the x axis
  virtual int GetKxMomentum() const;

  // get the maximum momentum along the x axis (i.e. the number of momentum sectors)
  // 
  // return avlue = maximum momentum along the x axis
  virtual int GetMaxXMomentum() const;

  // get the maximum momentum along the x axis (i.e. the number of momentum sectors)
  // 
  // return avlue = maximum momentum along the x axis
  virtual int GetNbrSites() const;
  
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = subsystem momentum along the x direction
  // kySector = subsystem momentum along the x direction
  // groundState = reference on the total system ground state  
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  // virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position, sanitizing the input data first
  // 
  // xPosition = x coordinate of the unit cell
  // yPosition = y coordinate of the unit cell
  // orbitalIndex = index of the orbital within the unit cell
  // return value = linearized index
  virtual int GetLinearizedIndexSafe(int xPosition, int orbitalIndex);

  // get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position
  // 
  // xPosition = x coordinate of the unit cell
  // yPosition = y coordinate of the unit cell
  // orbitalIndex = index of the orbital within the unit cell
  // return value = linearized index
  virtual int GetLinearizedIndex(int xPosition, int orbitalIndex);

  // get the lattice position from the linearized index (e.g. used for the creation/annihilation operators)
  // 
  // index = linearized index
  // xPosition = reference on the x coordinate of the unit cell
  // yPosition =reference on the  y coordinate of the unit cell
  // orbitalIndex = reference on the index of the orbital within the unit cell
  virtual void GetLinearizedIndex(int index, int& xPosition, int& orbitalIndex);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int maxMomentum);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions);

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int currentSite, long pos);

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(unsigned long& state, double& coefficient, int& nbrTranslationX);

  
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  inline int FindOrbitSize(unsigned long stateDescription);

  // apply a single translation in the x direction for a state description
  //
  // stateDescription = reference on the state description
  virtual void ApplySingleXTranslation(unsigned long& stateDescription);

  // get the fermonic sign when performing a single translation in the x direction on a state description, and apply the single translation
  //
  // stateDescription = reference on state description
  // return value = 0 if the sign is +1, 1 if the sign is -1
  unsigned long GetSignAndApplySingleXTranslation(unsigned long& stateDescription);

   // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(unsigned long memory);
  
  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in complementary Hilbert space
  // nbrIndex = number of indices to consider in complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  //  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorusWithMagneticTranslations* complementaryHilbertSpace,  		  								  ParticleOnTorusWithMagneticTranslations* destinationHilbertSpace, ComplexVector& groundState,  HermitianMatrix* densityMatrix);
  
};


// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int FermionOnLatticeRealSpaceAnd1DTranslation::SymmetrizeAdAdResult(unsigned long& state, double& coefficient, int& nbrTranslationX)
{
  state = this->FindCanonicalForm(state, nbrTranslationX);
  int TmpMaxMomentum = this->NbrSite;
  while ((state >> TmpMaxMomentum) == 0x0ul)
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslationX = (this->MaxXMomentum - nbrTranslationX) % this->MaxXMomentum;
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >>  nbrTranslationX) & 0x1ul))); 
    }
  return TmpIndex;
}


// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long FermionOnLatticeRealSpaceAnd1DTranslation::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;  
  unsigned long TmpStateDescription;  
  nbrTranslationX = 0;
  TmpStateDescription = stateDescription;
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription); 
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslationX = n;
	}
    }

  return CanonicalState;
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool FermionOnLatticeRealSpaceAnd1DTranslation::TestMomentumConstraint(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  unsigned long TmpStateDescription2 = stateDescription;
  int XSize = 1;
  unsigned long TmpSign = this->GetSignAndApplySingleXTranslation(TmpStateDescription);   
  unsigned long TmpSign2 = 0x0ul;
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      TmpSign ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
    }
  if ((((this->XMomentum * XSize) + ((((int) TmpSign) * this->MaxXMomentum) >> 1)) % this->MaxXMomentum) != 0)
    return false;

  return true;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int FermionOnLatticeRealSpaceAnd1DTranslation::FindOrbitSize(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  return XSize;
}

// apply a single translation in the x direction for a state description
//
// stateDescription = reference on the state description

inline void FermionOnLatticeRealSpaceAnd1DTranslation::ApplySingleXTranslation(unsigned long& stateDescription)
{
  stateDescription = (stateDescription >> this->StateXShift) | ((stateDescription & this->XMomentumMask) << this->ComplementaryStateXShift);
}

// get the fermonic sign when performing a single translation in the x direction on a state description, and apply the single translation
//
// stateDescription = reference on state description
// return value = 0 if the sign is +1, 1 if the sign is -1

inline unsigned long FermionOnLatticeRealSpaceAnd1DTranslation::GetSignAndApplySingleXTranslation(unsigned long& stateDescription)
{
  unsigned long TmpSign =  stateDescription >> this->StateXShift;
  stateDescription = TmpSign | ((stateDescription & this->XMomentumMask) << this->ComplementaryStateXShift);
#ifdef __64_BITS__
  TmpSign ^= (TmpSign >> 32);
#endif
  TmpSign ^= (TmpSign >> 16);
  TmpSign ^= (TmpSign >> 8);
  TmpSign ^= (TmpSign >> 4);
  TmpSign ^= (TmpSign >> 2);
  TmpSign ^= (TmpSign >> 1);
  TmpSign &= this->NbrFermionsParity;
  return TmpSign;
}

// get the momentum along the x axis
// 
// return avlue = momentum along the x axis

inline int FermionOnLatticeRealSpaceAnd1DTranslation::GetKxMomentum() const
{
  return this->XMomentum;
}

// get the maximum momentum along the x axis (i.e. the number of momentum sectors)
// 
// return value = maximum momentum along the x axis

inline int FermionOnLatticeRealSpaceAnd1DTranslation::GetMaxXMomentum() const
{
  return this->MaxXMomentum;
}

// get the number of sites 
// 
// return value = number of sites

inline int FermionOnLatticeRealSpaceAnd1DTranslation::GetNbrSites() const
{
  return this->NbrSite;
}
  


// get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position, sanitizing the input data first
// 
// xPosition = x coordinate of the unit cell
// yPosition = y coordinate of the unit cell
// orbitalIndex = index of the orbital within the unit cell
// return value = linearized index

inline int FermionOnLatticeRealSpaceAnd1DTranslation::GetLinearizedIndexSafe(int xPosition,  int orbitalIndex)
{
  orbitalIndex %= this->NbrSitePerUnitCell;
  if (orbitalIndex < 0)
    orbitalIndex +=  this->NbrSitePerUnitCell;
  xPosition %= this->MaxXMomentum;
  if (xPosition < 0)
    xPosition +=  this->MaxXMomentum;
  return this->GetLinearizedIndex(xPosition, orbitalIndex); 
}

// get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position
// 
// xPosition = x coordinate of the unit cell
// yPosition = y coordinate of the unit cell
// orbitalIndex = index of the orbital within the unit cell
// return value = linearized index

inline int FermionOnLatticeRealSpaceAnd1DTranslation::GetLinearizedIndex(int xPosition,  int orbitalIndex)
{
  return (xPosition * this->NbrSitePerUnitCell) + orbitalIndex;
}

// get the lattice position from the linearized index (e.g. used for the creation/annihilation operators)
// 
// index = linearized index
// xPosition = reference on the x coordinate of the unit cell
// yPosition =reference on the  y coordinate of the unit cell
// orbitalIndex = reference on the index of the orbital within the unit cell

inline void FermionOnLatticeRealSpaceAnd1DTranslation::GetLinearizedIndex(int index, int& xPosition,  int& orbitalIndex)
{
  orbitalIndex = index % this->NbrSitePerUnitCell;
  xPosition = index / this->NbrSitePerUnitCell;
}

#endif


