////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of fermions on lattice                        //
//                       in real space with C4 symmetry                       //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 27/02/2018                      //
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


#ifndef FERMIONONSQUARELATTICEREALSPACEANDC4SYMMETRY_H
#define FERMIONONSQUARELATTICEREALSPACEANDC4SYMMETRY_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd1DTranslation.h"

#include <iostream>
#include <bitset>


class  FermionOnSquareLatticeRealSpaceAndC4Symmetry : public FermionOnLatticeRealSpaceAnd1DTranslation
{

 protected:

  // mask for the site corresponding to the C4 rotation center
  unsigned long CenterSiteMask;
  // complementary mask of CenterSiteMask
  unsigned long ComplementaryCenterSiteMask;
  // position of the site corresponding to the C4 rotation center
  int CenterSitePosition;

 public:

  // default constructor
  // 
  FermionOnSquareLatticeRealSpaceAndC4Symmetry ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSitesX = number of sites along the x (or y direction)
  // symmetrySector = C4 symmetry sector (either 0, 1, 2 or 3)
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeRealSpaceAndC4Symmetry (int nbrFermions, int nbrSitesX, int symmetrySector, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeRealSpaceAndC4Symmetry(const FermionOnSquareLatticeRealSpaceAndC4Symmetry& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeRealSpaceAndC4Symmetry ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeRealSpaceAndC4Symmetry& operator = (const FermionOnSquareLatticeRealSpaceAndC4Symmetry& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

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

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

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
  
};


// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int FermionOnSquareLatticeRealSpaceAndC4Symmetry::SymmetrizeAdAdResult(unsigned long& state, double& coefficient, int& nbrTranslationX)
{
  state = this->FindCanonicalForm(state, nbrTranslationX);
  int TmpMaxMomentum = this->NbrSite - 1;
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

inline unsigned long FermionOnSquareLatticeRealSpaceAndC4Symmetry::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX)
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

inline bool FermionOnSquareLatticeRealSpaceAndC4Symmetry::TestMomentumConstraint(unsigned long stateDescription)
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

inline int FermionOnSquareLatticeRealSpaceAndC4Symmetry::FindOrbitSize(unsigned long stateDescription)
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

inline void FermionOnSquareLatticeRealSpaceAndC4Symmetry::ApplySingleXTranslation(unsigned long& stateDescription)
{
  stateDescription = (((stateDescription & this->ComplementaryCenterSiteMask) >> this->StateXShift) 
		      | ((stateDescription & this->XMomentumMask) << this->ComplementaryStateXShift)
		      | (stateDescription & this->CenterSiteMask));
}

// get the fermonic sign when performing a single translation in the x direction on a state description, and apply the single translation
//
// stateDescription = reference on state description
// return value = 0 if the sign is +1, 1 if the sign is -1

inline unsigned long FermionOnSquareLatticeRealSpaceAndC4Symmetry::GetSignAndApplySingleXTranslation(unsigned long& stateDescription)
{
  unsigned long TmpSign =  (stateDescription & this->ComplementaryCenterSiteMask) >> this->StateXShift;
  unsigned long TmpSign2 = stateDescription & this->CenterSiteMask;
  stateDescription = TmpSign | ((stateDescription & this->XMomentumMask) << this->ComplementaryStateXShift) | TmpSign2;
#ifdef __64_BITS__
  TmpSign ^= (TmpSign >> 32);
#endif
  TmpSign ^= (TmpSign >> 16);
  TmpSign ^= (TmpSign >> 8);
  TmpSign ^= (TmpSign >> 4);
  TmpSign ^= (TmpSign >> 2);
  TmpSign ^= (TmpSign >> 1);
  TmpSign &= this->NbrFermionsParity ^ (TmpSign2 >> this->CenterSitePosition);
  return TmpSign;
}

// get the maximum momentum along the x axis (i.e. the number of momentum sectors)
// 
// return value = maximum momentum along the x axis

inline int FermionOnSquareLatticeRealSpaceAndC4Symmetry::GetMaxXMomentum() const
{
  return 4;
}

// get the number of sites 
// 
// return value = number of sites

inline int FermionOnSquareLatticeRealSpaceAndC4Symmetry::GetNbrSites() const
{
  return this->NbrSite;
}
  


// get the linearized index (e.g. used for the creation/annihilation operators) from the lattice position, sanitizing the input data first
// 
// xPosition = x coordinate of the unit cell
// yPosition = y coordinate of the unit cell
// orbitalIndex = index of the orbital within the unit cell
// return value = linearized index

inline int FermionOnSquareLatticeRealSpaceAndC4Symmetry::GetLinearizedIndexSafe(int xPosition,  int orbitalIndex)
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

inline int FermionOnSquareLatticeRealSpaceAndC4Symmetry::GetLinearizedIndex(int xPosition,  int orbitalIndex)
{
  return (xPosition * this->NbrSitePerUnitCell) + orbitalIndex;
}

// get the lattice position from the linearized index (e.g. used for the creation/annihilation operators)
// 
// index = linearized index
// xPosition = reference on the x coordinate of the unit cell
// yPosition =reference on the  y coordinate of the unit cell
// orbitalIndex = reference on the index of the orbital within the unit cell

inline void FermionOnSquareLatticeRealSpaceAndC4Symmetry::GetLinearizedIndex(int index, int& xPosition,  int& orbitalIndex)
{
  orbitalIndex = index % this->NbrSitePerUnitCell;
  xPosition = index / this->NbrSitePerUnitCell;
}

#endif


