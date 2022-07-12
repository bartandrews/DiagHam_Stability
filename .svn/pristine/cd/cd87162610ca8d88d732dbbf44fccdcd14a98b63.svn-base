////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//       in real space with translation invariance in two directions and      //
//                               Sz<->-Sz symmetry                            //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 03/11/2014                      //
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


#ifndef FERMIONONLATTICEWITHSPINSZSYMMETRYREALSPACEAND2DTRANSLATION_H
#define FERMIONONLATTICEWITHSPINSZSYMMETRYREALSPACEAND2DTRANSLATION_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"

#include <iostream>



class FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation : public FermionOnLatticeWithSpinRealSpaceAnd2DTranslation
{

  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;

 protected:

  // sign of the parity sector for the Sz<->-Sz symmetry
  double SzParitySign;

  // sign of the parity sector for the Sz<->-Sz symmetry, 0 for 1, 1 for -1
  unsigned long SzParity;

 public:

  // default constructor
  // 
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (int nbrFermions, int nbrSite, bool minusSzParity, int xMomentum, int maxXMomentum,
							       int yMomentum, int maxYMomentum, unsigned long memory = 10000000);
  
  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the value of Sz
  // nbrSite = number of sites
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (int nbrFermions, int totalSpin, int nbrSite, bool minusSzParity, int xMomentum, int maxXMomentum,
							       int yMomentum, int maxYMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation(const FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation& operator = (const FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();
    
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in given momentum and Sz sectors.
  //
  // nbrParticleSector = number of particles that belong to the subsytem
  // szParitySector = Sz parity sector (can be either -1 or +1)
  // kxSector = subsystem momentum along the x direction
  // kySector = subsystem momentum along the x direction
  // szSector  = twice the total Sz value of the subsytem
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int szSector, int szParitySector, int kxSector, int kySector, 
									 ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. 
  // The density matrix is only evaluated in given momentum and Sz sectors.
  //
  // nbrParticleSector = number of particles that belong to the subsytem
  // szSector  = twice the total Sz value of the subsytem
  // szParitySector = Sz parity sector (can be either -1 or +1)
  // kxSector = subsystem momentum along the x direction
  // kySector = subsystem momentum along the x direction
  // nbrGroundStates = number of projectors
  // groundStates = array of degenerate groundstates associated to each projector
  // weights = array of weights in front of each projector
  // architecture = pointer to the architecture to use parallelized algorithm
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int szSector, int szParitySector, int kxSector, int kySector, 
									 int nbrGroundStates, ComplexVector* groundStates, double* weights, AbstractArchitecture* architecture);

  // convert a given state from the n-body basis with a fized Sz parity to the full n-body basis
  //
  // state = reference on the vector to convert
  // targetNbodyBasis = reference on the nbody-basis where the final state will be expressed
  // return value = converted vector
  ComplexVector ConvertToNbodyBasis(ComplexVector& state, FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* targetNbodyBasis);
  
  // convert a given state from the full n-body basis to the current n-body basis with a fized Sz parity 
  //
  // state = reference on the vector to convert
  // inputNbodyBasis = reference on the nbody-basis where the inital state is expressed
  // return value = converted vector
  ComplexVector ConvertFromNbodyBasis(ComplexVector& state, FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* inputNbodyBasis);
  
 protected:

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

  // generate all states corresponding to the constraints (core part of the method)
  //
  // return value = Hilbert space dimension
  virtual long CoreGenerateStates();

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(unsigned long& state, double& coefficient, 
				   int& nbrTranslationX, int& nbrTranslationY);

  
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslationX = reference on the number of translations to apply in the x direction to the resulting state to obtain the return orbit describing state
  // nbrTranslationY = reference on the number of translations to apply in the y direction to the resulting state to obtain the return orbit describing state
  // additionalSign = reference on the additional sign coming from symmetries beyond translation
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY, int& nbrSzSymmetry);//double& additionalSign);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(unsigned long stateDescription);

  // find the reordering sign when applying a sequence of discrete symmetries
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslationX = number of translations to apply  in the x direction
  // nbrTranslationY = number of translations to apply in the y direction to
  // nbrSpinFlip = number of full spin flip to apply
  // return value = reordering sign (0 for +1, 1 for -1)
  virtual unsigned long FindReorderingSign(unsigned long stateDescription, int nbrTranslationX, int nbrTranslationY, int nbrSzSymmetry);

  // Apply the Sz operator to flip all the spins
  //
  // stateDescription = reference on state description
  // stateDescription = state that has to be converted to its canonical expression
  virtual void ApplySzSymmetry (unsigned long& stateDescription);
  
  // get the fermonic sign when performing a flip all the spins, and apply the flip sign 
  //
  // stateDescription = reference on state description
  // return value = 0 if the sign is +1, 1 if the sign is -1
  virtual unsigned long GetSignAndApplySzSymmetry (unsigned long& stateDescription);

  // core part of the evaluation density matrix particle partition calculation
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorusWithSpinAndMagneticTranslations* complementaryHilbertSpace,  
								  ParticleOnTorusWithSpinAndMagneticTranslations* destinationHilbertSpace,
								  ComplexVector& groundState, HermitianMatrix* densityMatrix);

  // core part of the evaluation density matrix particle partition calculation involving a sum of projectors
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // nbrGroundStates = number of projectors
  // groundStates = array of degenerate groundstates associated to each projector
  // weights = array of weights in front of each projector
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorusWithSpinAndMagneticTranslations* complementaryHilbertSpace,
								  ParticleOnTorusWithSpinAndMagneticTranslations* destinationHilbertSpace,
								  int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix* densityMatrix);

};


// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state  

inline int FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::SymmetrizeAdAdResult(unsigned long& state, double& coefficient, 
											     int& nbrTranslationX, int& nbrTranslationY)
{
  int NbrSzSymmetry;
  state = this->FindCanonicalForm(state, nbrTranslationX, nbrTranslationY, NbrSzSymmetry);
  int TmpMaxMomentum = 2 * this->NbrSite - 1;
  while ((state >> TmpMaxMomentum) == 0x0ul)
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> ((((nbrTranslationY * this->MaxXMomentum) + nbrTranslationX) * 2) + NbrSzSymmetry)) & 0x1u)));
      if (NbrSzSymmetry == 1) 
	coefficient *= this->SzParitySign;  
    }
  return TmpIndex;
}


// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations toapply  in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to apply in the y direction to the resulting state to obtain the return orbit describing state
// nbrSpinFlip = reference on the number of full spin flip to apply to the resulting state to obtain the return orbit describing state
// additionalSign = reference on the additional sign coming from symmetries beyond translation
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY,
												    int& nbrSzSymmetry)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;  
  unsigned long TmpStateDescription;  
  nbrTranslationX = 0;
  nbrTranslationY = 0;
  nbrSzSymmetry = 0;
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
  stateDescription = stateDescriptionReference;
  TmpStateDescription = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription);
  if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      nbrTranslationX = 0;	      
      nbrTranslationY = 0;	      
      nbrSzSymmetry = 1;
    }
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslationX = n;	      
	  nbrTranslationY = 0;	      
	  nbrSzSymmetry = 1;
	}
    }
  this->ApplySzSymmetry(stateDescription);
  for (int m = 1; m < this->MaxYMomentum; ++m)
    {
      this->ApplySingleYTranslation(stateDescription);      
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslationX = 0;	      
	  nbrTranslationY = m;	      
	  nbrSzSymmetry = 1;
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
	      nbrSzSymmetry = 1;
	    }
	}
    }
  return CanonicalState;
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::TestMomentumConstraint(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  unsigned long TmpStateDescription2 = stateDescription;
  unsigned long TmpStateDescription3 = stateDescription;
  int XSize = 1;
  unsigned long TmpSign = this->GetSignAndApplySingleXTranslation(TmpStateDescription);   
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      TmpSign ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
    }
  if (TmpSign == 0x0ul)
    {
      if (((this->XMomentum * XSize) % this->MaxXMomentum) != 0)
	return false;
    }
  else
    {
      if ((((this->XMomentum * XSize) + this->MaxXMomentum) % (2 * this->MaxXMomentum)) != 0)
	return false;
    }
  int YSize = this->MaxYMomentum;
  int TmpXSize = 0;
  TmpStateDescription2 = stateDescription;
  for (int m = 1; m < YSize; ++m)
    {
      TmpSign ^= this->GetSignAndApplySingleYTranslation(TmpStateDescription2); 
      TmpStateDescription = TmpStateDescription2;
      TmpXSize = 0;
      unsigned long TmpSign2 = TmpSign;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  TmpSign2 ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	  TmpSign = TmpSign2;
	}
      else
	{
	  TmpXSize = -1;
	}
    } 
  if (TmpXSize < 0)
    {
      TmpXSize = 0;
      TmpSign ^= this->GetSignAndApplySingleYTranslation(TmpStateDescription2); 
    }
   if (TmpSign == 0x0ul)
    {
      if ((((this->YMomentum * YSize * this->MaxXMomentum)
	    + (this->XMomentum * TmpXSize * this->MaxYMomentum)) % (this->MaxXMomentum * this->MaxYMomentum)) != 0)
	return false;
    }
   else
     {
       if (((((this->YMomentum * YSize * this->MaxXMomentum)
	      + (this->XMomentum * TmpXSize * this->MaxYMomentum)) + (this->MaxXMomentum * this->MaxYMomentum)) % (2 * this->MaxXMomentum * this->MaxYMomentum)) != 0)
	return false;
     }

  TmpStateDescription2 = stateDescription;
  TmpSign = this->GetSignAndApplySzSymmetry(TmpStateDescription2);
  if (stateDescription == TmpStateDescription2)
    {
      if ((this->SzParity ^ TmpSign) == 0x0ul)
	return true;
      else
	return false;
    }

  int XSize2 = 1;
  TmpStateDescription = TmpStateDescription2;
  unsigned long TmpSign2 = TmpSign;
  TmpSign2 ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      TmpSign2 ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
    }  
  if (XSize2 < XSize)
    {
      if ((this->SzParity ^ TmpSign2) != 0x0ul)
	{
	  if ((((this->XMomentum * XSize2 * 2 * this->MaxYMomentum) + (this->MaxXMomentum * this->MaxYMomentum)) % (2 * this->MaxXMomentum * this->MaxYMomentum)) != 0)
	    return false;
	  else
	    return true;
	}
      else
	{
	  if ((((this->XMomentum * XSize2 * this->MaxYMomentum)) % (this->MaxXMomentum * this->MaxYMomentum)) != 0)
	    return false;
	  else
	    return true;  
	}
   }
  int YSize2 = YSize;
  TmpXSize = 0;
  TmpStateDescription2 = stateDescription;
  TmpSign = this->GetSignAndApplySzSymmetry(TmpStateDescription2);
  for (int m = 1; m < YSize2; ++m)
    {      
      TmpSign ^= this->GetSignAndApplySingleYTranslation(TmpStateDescription2); 
      TmpStateDescription = TmpStateDescription2;
      TmpXSize = 0; 
      TmpSign2 = TmpSign;
      while ((TmpXSize < XSize2) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  TmpSign2 ^= this->GetSignAndApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize2)
	{
	  YSize2 = m;
	  TmpSign = TmpSign2;
	}
      else
	{
	  TmpXSize = 0;
	}
    } 

  if (YSize == YSize2)
    return true;

  if ((this->SzParity ^ TmpSign) != 0x0ul)
    {
      if ((((this->YMomentum * YSize2 * 2 * this->MaxXMomentum)
	    + (this->XMomentum * TmpXSize * 2 * this->MaxYMomentum) + (this->MaxXMomentum * this->MaxYMomentum)) % (2 * this->MaxXMomentum * this->MaxYMomentum)) != 0)
	return false;
      else
	return true;
    }
  if ((((this->YMomentum * YSize2 * this->MaxXMomentum)
	+ (this->XMomentum * TmpXSize * this->MaxYMomentum)) % (this->MaxXMomentum * this->MaxYMomentum)) != 0)
    return false;
  else
    return true;  
  return true;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::FindOrbitSize(unsigned long stateDescription)
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
  TmpStateDescription2 = stateDescription;
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(TmpStateDescription2); 
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

  TmpStateDescription2 = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription2);
  if (stateDescription == TmpStateDescription2)
    return (XSize * YSize);  

  int XSize2 = 1;
  TmpStateDescription = TmpStateDescription2;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  
  if (XSize2 != XSize)
    {
      return (XSize * YSize);
    }
  TmpStateDescription2 = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription2);
  int YSize2 = YSize;
  int TmpXSize;
  for (int m = 1; m < YSize2; ++m)
    {
      this->ApplySingleYTranslation(TmpStateDescription2); 
      TmpStateDescription = TmpStateDescription2;
      int TmpXSize = 0;
      while ((TmpXSize < XSize2) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(TmpStateDescription);      
	}
      if (TmpXSize < XSize2)
	{
	  return (XSize * YSize);
	}
    }
  return (2 * XSize * YSize);
}

// find the reordering sign when applying a sequence of discrete symmetries
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = number of translations to apply  in the x direction
// nbrTranslationY = number of translations to apply in the y direction to
// nbrSpinFlip = number of full spin flip to apply
// return value = reordering sign (0 for +1, 1 for -1)

inline unsigned long FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::FindReorderingSign(unsigned long stateDescription, int nbrTranslationX, int nbrTranslationY,
												     int nbrSzSymmetry)
{
  unsigned long TmpSign = 0x0ul;
  if (nbrSzSymmetry == 1)
    {
      TmpSign ^= this->GetSignAndApplySzSymmetry(stateDescription);
    }
  for (int m = 0; m < nbrTranslationY; ++m)
    {
      TmpSign ^=  this->GetSignAndApplySingleYTranslation(stateDescription); 
    }
  for (int n = 0; n < nbrTranslationX; ++n)
    {
      TmpSign ^=  this->GetSignAndApplySingleXTranslation(stateDescription); 
    }
  return TmpSign;
}

// Apply the Sz operator to flip all the spins
//
// stateDescription = reference on state description

inline void FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::ApplySzSymmetry (unsigned long& stateDescription)
{
  
  unsigned long TmpState = stateDescription;
  stateDescription = ((TmpState >> 1) ^ TmpState) & FERMION_LATTICE_REALSPACE_SU2_SZ_MASK;
  stateDescription |= stateDescription << 1;
  stateDescription ^= TmpState; 
}

// get the fermonic sign when performing a flip all the spins, and apply the flip sign 
//
// stateDescription = reference on state description
// return value = 0 if the sign is +1, 1 if the sign is -1

inline unsigned long FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation::GetSignAndApplySzSymmetry (unsigned long& stateDescription)
{
  unsigned long TmpState = stateDescription;
  stateDescription = ((TmpState >> 1) ^ TmpState) & FERMION_LATTICE_REALSPACE_SU2_SZ_MASK;
  stateDescription |= stateDescription << 1;
  stateDescription ^= TmpState; 
  // compute the parity of the pair number
  TmpState &= (TmpState >> 1);
  TmpState &= FERMION_LATTICE_REALSPACE_SU2_SZ_MASK;
#ifdef __64_BITS__
  TmpState ^= (TmpState >> 32);
#endif
  TmpState ^= (TmpState >> 16);
  TmpState ^= (TmpState >> 8);
  TmpState ^= (TmpState >> 4);
  TmpState ^= (TmpState >> 2);
  return (TmpState & 0x01ul);
}

#endif


