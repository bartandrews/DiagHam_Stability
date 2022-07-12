////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of bosons hardcore on lattice                 //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Antoine Sterdyniak                    //
//                                                                            //
//                        last modification : 30/06/2016                      //
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


#ifndef BOSONONLATTICEGUTZWILLERPROJECTIONREALSPACEAND2DTRANSLATIONLONG_H
#define BOSONONLATTICEGUTZWILLERPROJECTIONREALSPACEAND2DTRANSLATIONLONG_H

#include "config.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslationsLong.h"

#include <iostream>



class BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong : public   FermionOnTorusWithMagneticTranslationsLong
{

 protected:
  
  // number of bosons
  int NbrBosons;	
  
  // number of bosons + 1
  int IncNbrBosons;
  
  // total number of sites
  int NbrSite;
  // number of sites per unit cell
  int NbrSitePerUnitCell;
  
  // number of momentum sectors in the x direction 
  int MaxXMomentum;
  // bit shift that has to applied to perform a translation in the x direction 
  int StateXShift;
  // binary mask for the StateXShift first bits 
  ULONGLONG XMomentumMask;
  // bit shift to apply to move the first StateXShift bits at the end of a state description
  int ComplementaryStateXShift;

  // number of momentum sectors in the y direction 
  int MaxYMomentum;
  // bit shift that has to applied to perform a translation in the y direction 
  int StateYShift;
  // binary mask for the StateYShift first bits 
  ULONGLONG YMomentumMask;
  // binary mask for the StateYShift first bits of each group
  ULONGLONG YMomentumFullMask;
  // binary mask for the ~YMomentumFullMask
  ULONGLONG ComplementaryYMomentumFullMask;
  // bit shift to apply to move the first StateYShift bits at the end of a state description
  int ComplementaryStateYShift;
  // number of bits that are related by a translation along the y direction 
  int YMomentumBlockSize;
  // binary mask corresponding to YMomentumBlockSize
  ULONGLONG YMomentumBlockMask;
  // number of independant blockse related by translations in the y direction 
  int NbrYMomentumBlocks;

 protected:
    
  // target space for operations leaving the Hilbert space
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong * TargetSpace;


 public:

  // default constructor		
  // 
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrSite = number of sites
  // xMomentum = momentum sector in the x direction
  // maxXMomentum = maximum momentum in the x direction
  // yMomentum = momentum sector in the y direction
  // maxYMomentum = maximum momentum in the y direction 
  // memory = amount of memory granted for precalculations
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong (int nbrBosons, int nbrSite, int xMomentum, int maxXMomentum,
							       int yMomentum, int maxYMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong (const BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong& bosons);

  // destructor
  //
  ~BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong& operator = (const BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong& bosons);

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

  // apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n = first index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A (int index, int n) ;

  // apply a^+_m operator to the state produced using AuAu method (without destroying it)
  //
  // m = first index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state 
  virtual int Ad (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n = first index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int n1, int n2) ;


  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);


 protected:


  // find canonical form of a state description
  //
  // stateDescription = unsigned integer describing the state
  // maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
  // nbrTranslation = number of translation needed to obtain the canonical form
  // yMomentum = state momentum value in the y direction
  // return value = canonical form of a state description
  virtual ULONGLONG FindCanonicalForm(ULONGLONG stateDescription, int& maxMomentum, int& nbrTranslation);
  

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(unsigned long memory);
  
  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(ULONGLONG & state, double& coefficient, 
				   int& nbrTranslationX, int& nbrTranslationY);

  
  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual ULONGLONG FindCanonicalFormAndTestMomentumConstraint(ULONGLONG stateDescription, int& nbrTranslationX, int& nbrTranslationY);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(ULONGLONG stateDescription);

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(ULONGLONG stateDescription);

  // apply a single translation in the x direction for a state description
  //
  // stateDescription = reference on the state description
  virtual void ApplySingleXTranslation(ULONGLONG& stateDescription);

  // apply a single translation in the y direction for a state description
  //
  // stateDescription = reference on the state description  
  virtual void ApplySingleYTranslation(ULONGLONG& stateDescription);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // currentSite = current site index
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long RawGenerateStates(int nbrBosons, int currentSite, long pos);
  
};


// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::SymmetrizeAdAdResult(ULONGLONG & state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  state = this->FindCanonicalFormAndTestMomentumConstraint(state, nbrTranslationX, nbrTranslationY);
  if (nbrTranslationX < 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpMaxMomentum = this->NbrSite;
  while ((state >> TmpMaxMomentum) == ((ULONGLONG) 0x0ul))
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
    }
  return TmpIndex;
}



// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline ULONGLONG BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::FindCanonicalForm(ULONGLONG stateDescription, int& nbrTranslationX, int& nbrTranslationY)
{
  ULONGLONG CanonicalState = stateDescription;
  ULONGLONG stateDescriptionReference = stateDescription;  
  ULONGLONG TmpStateDescription;  
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

inline bool BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::TestMomentumConstraint(ULONGLONG stateDescription)
{
  ULONGLONG  TmpStateDescription = stateDescription;
  ULONGLONG  TmpStateDescription2 = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
	this->ApplySingleXTranslation(TmpStateDescription);
    }
  if (((this->XMomentum * XSize) % this->MaxXMomentum) != 0)
    return false;
  int YSize = this->MaxYMomentum;
  int TmpXSize = 0;
  TmpStateDescription2 = stateDescription;
  for (int m = 1; m < YSize; ++m)
    {
      this->ApplySingleYTranslation(TmpStateDescription2);	
      TmpStateDescription = TmpStateDescription2;
      TmpXSize = 0;
      while ((TmpXSize < XSize) && (stateDescription != TmpStateDescription))
	{	  
	  ++TmpXSize;
	  this->ApplySingleXTranslation(TmpStateDescription);
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
  
  if ((((this->YMomentum * YSize * this->MaxXMomentum)
	+ (this->XMomentum * TmpXSize * this->MaxYMomentum)) % (this->MaxXMomentum * this->MaxYMomentum)) != 0)
    return false;
  return true;
}


// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::GetParticleStatistic()
{
  return ParticleOnSphere::BosonicStatistic;
}



// find the size of the orbit for a given state
//
// return value = orbit size

inline int  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::FindOrbitSize(ULONGLONG stateDescription)
{
  ULONGLONG TmpStateDescription = stateDescription;
  ULONGLONG TmpStateDescription2 = stateDescription;
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

inline void  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::ApplySingleXTranslation(ULONGLONG& stateDescription)
{
  stateDescription = (stateDescription >> this->StateXShift) | ((stateDescription & this->XMomentumMask) << this->ComplementaryStateXShift);
}

// apply a single translation in the y direction for a state description
//
// stateDescription = reference on the state description

inline void  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::ApplySingleYTranslation(ULONGLONG& stateDescription)
{
  stateDescription = (((stateDescription & this->ComplementaryYMomentumFullMask) >> this->StateYShift) | 
		      ((stateDescription & this->YMomentumFullMask) << this->ComplementaryStateYShift));
}


// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline ULONGLONG BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong::FindCanonicalFormAndTestMomentumConstraint(ULONGLONG stateDescription, int& nbrTranslationX, int& nbrTranslationY)
{
  ULONGLONG CanonicalState = this->FindCanonicalForm(stateDescription,nbrTranslationX,nbrTranslationY);
  if (this->TestMomentumConstraint(CanonicalState)  == false)
    nbrTranslationX = -1;

  /*  for (int i = 0; i < this->MaxMomentum; ++i)
    cout << (unsigned long)  ((stateDescription >> i) & ((ULONGLONG) 0x1ul)) << " ";

  cout <<endl;
  for (int i = 0; i < this->MaxMomentum; ++i)
    cout << (unsigned long)  ((CanonicalState >> i) & ((ULONGLONG) 0x1ul)) << " ";
  
    cout <<nbrTranslationX<<endl;*/
  return CanonicalState;
}


#endif


