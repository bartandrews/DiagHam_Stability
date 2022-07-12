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


#ifndef FERMIONONLATTICEREALSPACEAND2DTRANSLATION_H
#define FERMIONONLATTICEREALSPACEAND2DTRANSLATION_H

#include "config.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"

#include <iostream>
#include <bitset>


class  FermionOnLatticeRealSpaceAnd2DTranslation : public FermionOnTorusWithMagneticTranslations
{

  friend class FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation;
  
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

  // number of momentum sectors in the y direction 
  int MaxYMomentum;
  // bit shift that has to applied to perform a translation in the y direction 
  int StateYShift;
  // binary mask for the StateYShift first bits 
  unsigned long YMomentumMask;
  // binary mask for the StateYShift first bits of each group
  unsigned long YMomentumFullMask;
  // binary mask for the ~YMomentumFullMask
  unsigned long ComplementaryYMomentumFullMask;
  // bit shift to apply to move the first StateYShift bits at the end of a state description
  int ComplementaryStateYShift;
  // number of bits that are related by a translation along the y direction 
  int YMomentumBlockSize;
  // binary mask corresponding to YMomentumBlockSize
  unsigned long YMomentumBlockMask;
  // number of independant blockse related by translations in the y direction 
  int NbrYMomentumBlocks;

  // parity of the number of fermions, 0x1ul if even, 0x0ul if odd
  unsigned long NbrFermionsParity;

 public:

  // default constructor
  // 
  FermionOnLatticeRealSpaceAnd2DTranslation ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeRealSpaceAnd2DTranslation (int nbrFermions, int nbrSite, int xMomentum, int maxXMomentum,
					     int yMomentum, int maxYMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeRealSpaceAnd2DTranslation(const FermionOnLatticeRealSpaceAnd2DTranslation& fermions);

  // destructor
  //
  ~FermionOnLatticeRealSpaceAnd2DTranslation ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeRealSpaceAnd2DTranslation& operator = (const FermionOnLatticeRealSpaceAnd2DTranslation& fermions);

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
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // apply a^+_m operator to the state produced using AuAu method (without destroying it)
  //
  // m = first index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int Ad (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // find state index from a string
  //
  // stateDescription = string describing the state
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(char* stateDescription);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = subsystem momentum along the x direction
  // kySector = subsystem momentum along the x direction
  // groundState = reference on the total system ground state  
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // get an Hilbert space without the translation symmetries, but preserving any other properties
  //
  // return value = Hilbert space without translations
  virtual FermionOnLatticeRealSpace* GetHilbertSpaceWithoutTranslations();
  
  
  // get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position, sanitizing the input data first
  // 
  // xPosition = x coordinate of the unit cell
  // yPosition = y coordinate of the unit cell
  // orbitalIndex = index of the orbital within the unit cell
  // return value = linearized index
  virtual int GetLinearizedIndexSafe(int xPosition, int yPosition, int orbitalIndex);

  // get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position
  // 
  // xPosition = x coordinate of the unit cell
  // yPosition = y coordinate of the unit cell
  // orbitalIndex = index of the orbital within the unit cell
  // return value = linearized index
  virtual int GetLinearizedIndex(int xPosition, int yPosition, int orbitalIndex);

  // get the lattice position from the linearized index (e.g. used for the creation/annihilation operators)
  // 
  // index = linearized index
  // xPosition = reference on the x coordinate of the unit cell
  // yPosition =reference on the  y coordinate of the unit cell
  // orbitalIndex = reference on the index of the orbital within the unit cell
  virtual void GetLinearizedIndex(int index, int& xPosition, int& yPosition, int& orbitalIndex);

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
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(unsigned long& state, double& coefficient, 
				   int& nbrTranslationX, int& nbrTranslationY);

  
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY);

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

  // apply a single translation in the y direction for a state description
  //
  // stateDescription = reference on the state description  
  virtual void ApplySingleYTranslation(unsigned long& stateDescription);

  // get the fermonic sign when performing a single translation in the x direction on a state description, and apply the single translation
  //
  // stateDescription = reference on state description
  // return value = 0 if the sign is +1, 1 if the sign is -1
  unsigned long GetSignAndApplySingleXTranslation(unsigned long& stateDescription);

  // get the fermonic sign when performing a single translation in the y direction on a state description, and apply the single translation
  //
  // stateDescription = reference on state description
  // return value = 0 if the sign is +1, 1 if the sign is -1
  virtual unsigned long GetSignAndApplySingleYTranslation(unsigned long& stateDescription);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(unsigned long memory);
  
  // convert a state defined in the real space basis into a state in the (Kx,Ky) basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertToKxKyBasis(ComplexVector& state, ParticleOnSphere* space);

  // convert a state defined in the (Kx,Ky) basis into a state in the real space basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space);
  
    
  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in complementary Hilbert space
  // nbrIndex = number of indices to consider in complementary Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorusWithMagneticTranslations* complementaryHilbertSpace,  
								  ParticleOnTorusWithMagneticTranslations* destinationHilbertSpace,
								  ComplexVector& groundState,  HermitianMatrix* densityMatrix);

};


// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int FermionOnLatticeRealSpaceAnd2DTranslation::SymmetrizeAdAdResult(unsigned long& state, double& coefficient, 
									   int& nbrTranslationX, int& nbrTranslationY)
{
  state = this->FindCanonicalForm(state, nbrTranslationX, nbrTranslationY);
  int TmpMaxMomentum = this->NbrSite;
  while ((state >> TmpMaxMomentum) == 0x0ul)
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslationX = (this->MaxXMomentum - nbrTranslationX) % this->MaxXMomentum;
      nbrTranslationY = (this->MaxYMomentum - nbrTranslationY) % this->MaxYMomentum;
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> ((nbrTranslationY * this->MaxXMomentum) + nbrTranslationX)) & 0x1ul))); 
    }
  return TmpIndex;
}


// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long FermionOnLatticeRealSpaceAnd2DTranslation::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;  
  unsigned long TmpStateDescription;  
  nbrTranslationX = 0;
  nbrTranslationY = 0;
  TmpStateDescription = stateDescription;
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslationX = n;	      
	  nbrTranslationY = 0;	      
	}
    }
  for (int m = 1; m < this->MaxYMomentum; ++m)
    {
      this->ApplySingleYTranslation(stateDescription);      
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslationX = 0;	      
	  nbrTranslationY = m;	      
	}
      TmpStateDescription = stateDescription;
      for (int n = 1; n < this->MaxXMomentum; ++n)
	{
	  this->ApplySingleXTranslation(TmpStateDescription);      
	  if (TmpStateDescription < CanonicalState)
	    {
	      CanonicalState = TmpStateDescription;
	      nbrTranslationX = n;	      
	      nbrTranslationY = m;	      
	    }
	}
    }
  return CanonicalState;
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool FermionOnLatticeRealSpaceAnd2DTranslation::TestMomentumConstraint(unsigned long stateDescription)
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
  int YSize = this->MaxYMomentum;
  int TmpXSize = 0;
  TmpSign = 0x0ul;
  TmpStateDescription2 = stateDescription;
  for (int m = 1; m < YSize; ++m)
    {
      TmpSign ^= this->GetSignAndApplySingleYTranslation(TmpStateDescription2); 
      TmpSign2 = TmpSign;
      TmpStateDescription = TmpStateDescription2;
      TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  TmpSign2 ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	}
      else
	{
	  TmpXSize = 0;
	}
    } 
  if (YSize == this->MaxYMomentum)
    {
      TmpSign ^= this->GetSignAndApplySingleYTranslation(TmpStateDescription2); 
      TmpSign2 = TmpSign;
    }
  if ((((2 * this->YMomentum * YSize * this->MaxXMomentum)
	+ (2 * this->XMomentum * TmpXSize * this->MaxYMomentum)
	+ (((int) TmpSign2) * this->MaxXMomentum * this->MaxYMomentum)) % (2 * this->MaxXMomentum * this->MaxYMomentum)) != 0)
    return false;
  return true;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int FermionOnLatticeRealSpaceAnd2DTranslation::FindOrbitSize(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  unsigned long TmpStateDescription2 = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  int YSize = this->MaxYMomentum;
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(stateDescription); 
      TmpStateDescription = TmpStateDescription2;
      int TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	}
    }
  return (XSize * YSize);
}

// apply a single translation in the x direction for a state description
//
// stateDescription = reference on the state description

inline void FermionOnLatticeRealSpaceAnd2DTranslation::ApplySingleXTranslation(unsigned long& stateDescription)
{
  stateDescription = (stateDescription >> this->StateXShift) | ((stateDescription & this->XMomentumMask) << this->ComplementaryStateXShift);
}

// apply a single translation in the y direction for a state description
//
// stateDescription = reference on the state description

inline void FermionOnLatticeRealSpaceAnd2DTranslation::ApplySingleYTranslation(unsigned long& stateDescription)
{
//  cout << std::bitset<16>(stateDescription) << " " <<std::bitset<16>(this->ComplementaryYMomentumFullMask) << " " <<std::bitset<16>(this->StateYShift)<<" " <<std::bitset<16>(this->YMomentumFullMask)<< " " <<std::bitset<16>(this->ComplementaryStateYShift)<<endl;
  stateDescription = (((stateDescription & this->ComplementaryYMomentumFullMask) >> this->StateYShift) | ((stateDescription & this->YMomentumFullMask) << this->ComplementaryStateYShift));
}

// get the fermonic sign when performing a single translation in the x direction on a state description, and apply the single translation
//
// stateDescription = reference on state description
// return value = 0 if the sign is +1, 1 if the sign is -1

inline unsigned long FermionOnLatticeRealSpaceAnd2DTranslation::GetSignAndApplySingleXTranslation(unsigned long& stateDescription)
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

// get the fermonic sign when performing a single translation in the y direction on a state description, and apply the single translation
//
// stateDescription = reference on state description
// return value = 0 if the sign is +1, 1 if the sign is -1

inline unsigned long FermionOnLatticeRealSpaceAnd2DTranslation::GetSignAndApplySingleYTranslation(unsigned long& stateDescription)
{
  unsigned long TmpState = 0x0ul;
  unsigned long TmpSign =  0x0ul;
  unsigned long TmpSign2;
  unsigned long TmpSign3;
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      TmpSign2 = (stateDescription & this->YMomentumBlockMask) >> this->StateYShift;
      TmpSign3 = (stateDescription & this->YMomentumMask) << this->ComplementaryStateYShift;
      TmpState |= (TmpSign2 | TmpSign3) << (this->YMomentumBlockSize * i);
#ifdef __64_BITS__
      TmpSign2 ^= (TmpSign2 >> 32);
#endif
      TmpSign2 ^= (TmpSign2 >> 16);
      TmpSign2 ^= (TmpSign2 >> 8);
      TmpSign2 ^= (TmpSign2 >> 4);
      TmpSign2 ^= (TmpSign2 >> 2);
      TmpSign2 ^= (TmpSign2 >> 1);
#ifdef __64_BITS__
      TmpSign3 ^= (TmpSign3 >> 32);
#endif
      TmpSign3 ^= (TmpSign3 >> 16);
      TmpSign3 ^= (TmpSign3 >> 8);
      TmpSign3 ^= (TmpSign3 >> 4);
      TmpSign3 ^= (TmpSign3 >> 2);
      TmpSign3 ^= (TmpSign3 >> 1);
      TmpSign2 *= TmpSign3;
      TmpSign2 &= 0x1ul;
      TmpSign ^= TmpSign2;
      stateDescription >>= this->YMomentumBlockSize;
    }
  stateDescription = TmpState;
  return TmpSign;
}

// get an Hilbert space without the translation symmetries, but preserving any other properties
//
// return value = Hilbert space without translations

inline FermionOnLatticeRealSpace* FermionOnLatticeRealSpaceAnd2DTranslation::GetHilbertSpaceWithoutTranslations()
{
  return new FermionOnLatticeRealSpace(this->NbrFermions, this->NbrSite);
}


// get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position, sanitizing the input data first
// 
// xPosition = x coordinate of the unit cell
// yPosition = y coordinate of the unit cell
// orbitalIndex = index of the orbital within the unit cell
// return value = linearized index

inline int FermionOnLatticeRealSpaceAnd2DTranslation::GetLinearizedIndexSafe(int xPosition, int yPosition, int orbitalIndex)
{
  orbitalIndex %= this->NbrSitePerUnitCell;
  if (orbitalIndex < 0)
    orbitalIndex +=  this->NbrSitePerUnitCell;
  xPosition %= this->MaxXMomentum;
  if (xPosition < 0)
    xPosition +=  this->MaxXMomentum;
  yPosition %= this->MaxYMomentum;
  if (yPosition < 0)
    yPosition +=  this->MaxYMomentum;
  return this->GetLinearizedIndex(xPosition, yPosition, orbitalIndex); 
}

// get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position
// 
// xPosition = x coordinate of the unit cell
// yPosition = y coordinate of the unit cell
// orbitalIndex = index of the orbital within the unit cell
// return value = linearized index

inline int FermionOnLatticeRealSpaceAnd2DTranslation::GetLinearizedIndex(int xPosition, int yPosition, int orbitalIndex)
{
  return (((xPosition * this->MaxYMomentum) + yPosition) * this->NbrSitePerUnitCell) + orbitalIndex;
}

// get the lattice position from the linearized index (e.g. used for the creation/annihilation operators)
// 
// index = linearized index
// xPosition = reference on the x coordinate of the unit cell
// yPosition =reference on the  y coordinate of the unit cell
// orbitalIndex = reference on the index of the orbital within the unit cell

inline void FermionOnLatticeRealSpaceAnd2DTranslation::GetLinearizedIndex(int index, int& xPosition, int& yPosition, int& orbitalIndex)
{
  orbitalIndex = index % this->NbrSitePerUnitCell;
  index /= this->NbrSitePerUnitCell;
  xPosition = index / this->MaxYMomentum;
  yPosition = index % this->MaxYMomentum;
}

#endif


