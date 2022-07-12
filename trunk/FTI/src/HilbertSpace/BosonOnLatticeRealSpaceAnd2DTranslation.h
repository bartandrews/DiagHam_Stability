////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of bosons on lattice                          //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Antoine Sterdyniak                    //
//                                                                            //
//                        last modification : 11/09/2014                      //
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


#ifndef BOSONONLATTICEREALSPACEAND2DTRANSLATION_H
#define BOSONONLATTICEREALSPACEAND2DTRANSLATION_H

#include "config.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"

#include <iostream>



class BosonOnLatticeRealSpaceAnd2DTranslation : public BosonOnTorusWithMagneticTranslationsShort
{

 protected:
  
  // total number of sites
  int NbrSite;
  // number of sites per unit cell
  int NbrSitePerUnitCell;
  
  
  int NbrMomentum;
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
  // bit shift to apply to move the first StateYShift bits at the end of a state description
  int ComplementaryStateYShift;
  // number of bits that are related by a translation along the y direction 
  int YMomentumBlockSize;
  // binary mask corresponding to YMomentumBlockSize
  unsigned long YMomentumBlockMask;
  // number of independant blockse related by translations in the y direction 
  int NbrYMomentumBlocks;
  
  // temporary state used when applying operators
  unsigned long* TemporaryStateOperators;
  int TemporaryStateOperatorsKyMax;

  
  public:

  // default constructor		
  // 
  BosonOnLatticeRealSpaceAnd2DTranslation ();

  // basic constructor
  // 
  // nbrBosons = number of fermions
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  BosonOnLatticeRealSpaceAnd2DTranslation (int nbrBosons, int nbrSite, int xMomentum, int maxXMomentum,
					   int yMomentum, int maxYMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnLatticeRealSpaceAnd2DTranslation (const BosonOnLatticeRealSpaceAnd2DTranslation& bosons);
  
  // destructor
  //
  ~BosonOnLatticeRealSpaceAnd2DTranslation ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnLatticeRealSpaceAnd2DTranslation& operator = (const BosonOnLatticeRealSpaceAnd2DTranslation& bosons);

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
  
  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream & PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // kxSector = subsystem momentum along the x direction
  // kySector = subsystem momentum along the x direction
  // groundState = reference on the total system ground state  
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

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
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrFermions, int currentSite, long pos);


  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(unsigned long& state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state 
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY);


  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalFormAndTestMomentumConstraint(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY);

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
  virtual void ApplySingleXTranslation(unsigned long * stateDescription);

  // apply a single translation in the y direction for a state description
  //
  // stateDescription = reference on the state description  
  virtual void ApplySingleYTranslation(unsigned long * stateDescription);

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

inline int BosonOnLatticeRealSpaceAnd2DTranslation::SymmetrizeAdAdResult(unsigned long& state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  state = this->FindCanonicalForm(state, nbrTranslationX, nbrTranslationY);
  
  int TmpMaxMomentum =  this->FermionicMaxMomentum;
  while ((state >> TmpMaxMomentum) == 0x0ul)
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->ProdATemporaryStateNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
    }
   return TmpIndex;
}


 
// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long BosonOnLatticeRealSpaceAnd2DTranslation::FindCanonicalFormAndTestMomentumConstraint(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;  
  this->FermionToBoson(stateDescription, this->FermionicMaxMomentum, this->TemporaryState, this->TemporaryStateKyMax);
  for(int i =  this->TemporaryStateKyMax + 1; i < this->NbrMomentum ; i++)
    {
      this->TemporaryState[i] = 0;
    }
  
  this->FermionToBoson(stateDescription, this->FermionicMaxMomentum, this->ProdATemporaryState, this->ProdATemporaryStateKyMax);
  for(int i =  this->ProdATemporaryStateKyMax + 1; i < this->NbrMomentum ; i++)
    {
      this->ProdATemporaryState[i] = 0;
    }
  
  unsigned long TmpStateDescription;  
  nbrTranslationX = 0;
  nbrTranslationY = 0;
  for (int m = 0; (m < this->MaxYMomentum) && (stateDescriptionReference != stateDescription) ; ++m)
    {
      for(int i =  this->TemporaryStateKyMax + 1; i < this->NbrMomentum ; i++)
	{
	  this->TemporaryState[i] =  this->ProdATemporaryState[i];
	}
      TmpStateDescription = stateDescription;
      
      for (int n = 0; (n < this->MaxXMomentum) && (TmpStateDescription != stateDescription) ; ++n)
	{
	  if (TmpStateDescription < CanonicalState)
	    {
	      CanonicalState = TmpStateDescription;
	      nbrTranslationX = n;	      
	      nbrTranslationY = m;	      
	    }
	  
	  this->ApplySingleXTranslation(this->TemporaryState);      
          int TmpMomentumMax =  this->NbrMomentum - 1;
	  while (this->TemporaryState[ TmpMomentumMax] == 0x0ul)
	    --TmpMomentumMax;
          TmpStateDescription = this->BosonToFermion( this->TemporaryState, TmpMomentumMax);
	}
      this->ApplySingleYTranslation(this->ProdATemporaryState);
      int ProdATmpMomentumMax =  this->NbrMomentum - 1;
      while (this->ProdATemporaryState[ProdATmpMomentumMax] == 0x0ul)
	--ProdATmpMomentumMax;
      stateDescription = this->BosonToFermion(this->ProdATemporaryState,ProdATmpMomentumMax);
    }
  return CanonicalState;  
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool BosonOnLatticeRealSpaceAnd2DTranslation::TestMomentumConstraint(unsigned long stateDescription)
{
  this->FermionToBoson(stateDescription, this->FermionicMaxMomentum, this->TemporaryState, this->TemporaryStateKyMax);
  for(int i =  this->TemporaryStateKyMax + 1; i < this->NbrMomentum ; i++)
    {
      this->TemporaryState[i] = 0;
    }
  
  unsigned long TmpStateDescription = stateDescription;
  unsigned long TmpStateDescription2 = stateDescription;
  int XSize = 1;
  
  this->ApplySingleXTranslation(this->TemporaryState);      
  int TmpMomentumMax =  this->NbrMomentum - 1;
  while (this->TemporaryState[ TmpMomentumMax] == 0x0ul)
    --TmpMomentumMax;
  TmpStateDescription = this->BosonToFermion( this->TemporaryState, TmpMomentumMax);
  
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(this->TemporaryState);      
      int TmpMomentumMax =  this->NbrMomentum - 1;
      while (this->TemporaryState[ TmpMomentumMax] == 0x0ul)
	--TmpMomentumMax;
      TmpStateDescription = this->BosonToFermion( this->TemporaryState, TmpMomentumMax);
    }
  if (((this->KxMomentum * XSize) % this->MaxXMomentum) != 0)
    return false;
  int YSize = this->MaxYMomentum;
  int TmpXSize = 0;
  this->FermionToBoson(stateDescription, this->FermionicMaxMomentum, this->ProdATemporaryState, this->ProdATemporaryStateKyMax);
  for(int i =  this->ProdATemporaryStateKyMax + 1 ; i < this->NbrMomentum ; i++)
    {
      this->ProdATemporaryState[i] = 0;
    }
  TmpStateDescription2 = stateDescription;
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(this->ProdATemporaryState);	 
      int ProdATmpMomentumMax =  this->NbrMomentum - 1;
      while (this->ProdATemporaryState[ProdATmpMomentumMax] == 0x0ul)
	--ProdATmpMomentumMax;
      TmpStateDescription2 = this->BosonToFermion( this->ProdATemporaryState,ProdATmpMomentumMax);
      TmpStateDescription = TmpStateDescription2;
      for(int i = 0 ; i < this->NbrMomentum ; i++)
	{
	  this->TemporaryState[i] = this->ProdATemporaryState[i] ;
	}
      
      TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(this->TemporaryState);      
          int TmpMomentumMax =  this->NbrMomentum - 1;
	  while (this->TemporaryState[ TmpMomentumMax] == 0x0ul)
	    --TmpMomentumMax;
          TmpStateDescription = this->BosonToFermion( this->TemporaryState, TmpMomentumMax);
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
  if ((((this->KyMomentum * YSize * this->MaxXMomentum)
	+ (this->KxMomentum * TmpXSize * this->MaxYMomentum)) % (this->MaxXMomentum * this->MaxYMomentum)) != 0)
    return false;
  return true;
}


// apply a single translation to a bosonic state in its fermionic representation
//
// stateDescription = state to translate

inline void BosonOnLatticeRealSpaceAnd2DTranslation::ApplySingleXTranslation(unsigned long * stateDescription)
{
  unsigned long Tmp;
  for (int i = 0; i <  this->StateXShift ; i++)
    {
      Tmp = stateDescription[i];
      for (int k = 1; k < this->MaxXMomentum ; k++)
	{
	  stateDescription[i+(k-1)*this->StateXShift] = this->TemporaryState[i+k*this->StateXShift];
	}
      stateDescription[i + (this->MaxXMomentum - 1)*this->StateXShift] = Tmp;
    }
}

// apply a single translation to a bosonic state in its fermionic representation
//
// stateDescription = state to translate

inline void BosonOnLatticeRealSpaceAnd2DTranslation::ApplySingleYTranslation(unsigned long * stateDescription)
{
  unsigned long Tmp;
  for (int i = 0; i < this->NbrYMomentumBlocks; i++)
    { 	
      for (int j = 0; j<this->StateYShift ;j++)
	{
	  Tmp = stateDescription[i*this->StateXShift + (this->MaxYMomentum - 1) * this->StateYShift + j];
	  for (int k = this->MaxYMomentum - 1; k > 0; k--)
	    {
	      stateDescription[i*this->StateXShift + k * this->StateYShift + j] = stateDescription[i*this->StateXShift + j + (k-1) * this->StateYShift];			
	    }
	  stateDescription[i*this->StateXShift + j] = Tmp;
	}
    }
}



// find the size of the orbit for a given state
//
// return value = orbit size

inline int BosonOnLatticeRealSpaceAnd2DTranslation::FindOrbitSize(unsigned long stateDescription)
{ 
  this->FermionToBoson(stateDescription, this->FermionicMaxMomentum, this->TemporaryState, this->TemporaryStateKyMax);
  for(int i =  this->TemporaryStateKyMax + 1; i < this->NbrMomentum ; i++)
    {
      this->TemporaryState[i] = 0;
    }
  
  unsigned long TmpStateDescription = stateDescription;
  this->FermionToBoson(stateDescription, this->FermionicMaxMomentum, this->ProdATemporaryState, this->ProdATemporaryStateKyMax);
  for(int i =  this->ProdATemporaryStateKyMax + 1 ; i < this->NbrMomentum ; i++)
    {
      this->ProdATemporaryState[i] = 0;
    }

  unsigned long TmpStateDescription2 = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(this->TemporaryState);      
  int TmpMomentumMax =  this->NbrMomentum-1;
  while (this->TemporaryState[TmpMomentumMax] == 0x0ul)
    --TmpMomentumMax;
  TmpStateDescription = this->BosonToFermion(this->TemporaryState, TmpMomentumMax);
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(this->TemporaryState);      
      int TmpMomentumMax =  this->NbrMomentum - 1;
      while (this->TemporaryState[TmpMomentumMax] == 0x0ul)
	--TmpMomentumMax;
      TmpStateDescription = this->BosonToFermion( this->TemporaryState, TmpMomentumMax);
    }
  int YSize = this->MaxYMomentum;
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(this->ProdATemporaryState);	 
      int ProdATmpMomentumMax =  this->NbrMomentum - 1;
      while (this->ProdATemporaryState[ProdATmpMomentumMax] == 0x0ul)
	--ProdATmpMomentumMax;
      TmpStateDescription2 = this->BosonToFermion( this->ProdATemporaryState,ProdATmpMomentumMax);
      this->FermionToBoson(stateDescription, this->FermionicMaxMomentum, this->TemporaryState, this->TemporaryStateKyMax);
      for(int i =  this->TemporaryStateKyMax + 1 ; i < this->NbrMomentum ; i++)
	{
	  this->TemporaryState[i] = 0;
 	}
      
      TmpStateDescription = stateDescription;
      int TmpXSize = 0;
      while ((TmpXSize < XSize) && (TmpStateDescription2 != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(this->TemporaryState);      
          int TmpMomentumMax =  this->NbrMomentum - 1;
	  while (this->TemporaryState[TmpMomentumMax] == 0x0ul)
	    --TmpMomentumMax;
          TmpStateDescription = this->BosonToFermion(this->TemporaryState, TmpMomentumMax);
	}
      if (TmpXSize < XSize)
	{
	  YSize = m;
	}
    }
  return (XSize * YSize);
}




// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long BosonOnLatticeRealSpaceAnd2DTranslation::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;  
  this->FermionToBoson(stateDescription, this->FermionicMaxMomentum, this->TemporaryState, this->TemporaryStateKyMax);
  for(int i =  this->TemporaryStateKyMax + 1 ; i < this->NbrMomentum ; i++)
    {
      this->TemporaryState[i] = 0;
    }
    
   for(int i =  0 ; i < this->NbrMomentum ; i++)
    {
      this->TemporaryStateOperators[i] = this->TemporaryState[i];
    }
  
  this->FermionToBoson(stateDescription, this->FermionicMaxMomentum, this->TemporaryStateOperators, this->TemporaryStateOperatorsKyMax);
  for(int i =  this->TemporaryStateOperatorsKyMax + 1 ; i < this->NbrMomentum ; i++)
    {
      this->TemporaryStateOperators[i] = 0;
    }
  
  
  unsigned long TmpStateDescription;  
  nbrTranslationX = 0;
  nbrTranslationY = 0;
   for (int n = 1; n < this->MaxXMomentum; ++n)
	{
	  this->ApplySingleXTranslation(this->TemporaryState);      
          int TmpMomentumMax =  this->NbrMomentum - 1;
	  while (this->TemporaryState[TmpMomentumMax] == 0x0ul)
	    --TmpMomentumMax;
          TmpStateDescription = this->BosonToFermion( this->TemporaryState, TmpMomentumMax);
	  
	  if (TmpStateDescription < CanonicalState)
	    {
	      CanonicalState = TmpStateDescription;
	      nbrTranslationX = n;	      
	      nbrTranslationY = 0;	      

	    }
	}

  for (int m = 1; (m < this->MaxYMomentum) && (stateDescription != 0x0ul) ; ++m)
    {
      TmpStateDescription = stateDescriptionReference;
      this->ApplySingleYTranslation(this->TemporaryStateOperators);	 
      int ProdATmpMomentumMax =  this->NbrMomentum - 1;
      while (this->TemporaryStateOperators[ProdATmpMomentumMax] == 0x0ul)
	--ProdATmpMomentumMax;
      stateDescriptionReference = this->BosonToFermion(this->TemporaryStateOperators,ProdATmpMomentumMax); 
      for(int i =  0; i < this->NbrMomentum ; i++)
	{
	  this->TemporaryState[i] = this->TemporaryStateOperators[i];
 	}
      
      for (int n = 1; n < this->MaxXMomentum; ++n)
	{
	  
	  this->ApplySingleXTranslation(this->TemporaryState);      
          int TmpMomentumMax =  this->NbrMomentum - 1;
	  while (this->TemporaryState[TmpMomentumMax] == 0x0ul)
	    --TmpMomentumMax;
          TmpStateDescription = this->BosonToFermion( this->TemporaryState, TmpMomentumMax);
	  
	  if (TmpStateDescription < CanonicalState)
	    {
	      CanonicalState = TmpStateDescription;
	      nbrTranslationX = n;	      
	      nbrTranslationY = m;	      
	    }
	  if  (TmpStateDescription == stateDescription)
	    n = this->MaxXMomentum;
	}
    
      
      if (stateDescription == stateDescriptionReference)
	{
	  m = this->MaxYMomentum;
	}
      else
	{
	  if (stateDescriptionReference < CanonicalState)
	    {
	      CanonicalState = stateDescriptionReference;
	      nbrTranslationX = 0;	      
	      nbrTranslationY = m;
	    }
	}
    } 
  return CanonicalState;
}


#endif


