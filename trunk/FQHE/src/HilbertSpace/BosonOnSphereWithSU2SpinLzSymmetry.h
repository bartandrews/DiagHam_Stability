////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//     class for bosons on sphere with SU(2) spin and Lz<->-Lz symmetry       //
//                                                                            //
//                        last modification : 24/09/2016                      //
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


#ifndef BOSONONSPHEREWITHSU2SPINLZSYMMETRY_H
#define BOSONONSPHEREWITHSU2SPINLZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphere.h"

#include <iostream>

using std::cout;
using std::endl;



class BosonOnSphereWithSU2SpinLzSymmetry :  public BosonOnSphereWithSU2SpinSzSymmetry
{


 protected:


  // sign of the parity sector for the Lz<->-Lz symmetry
  double LzParitySign;

  // maximum index to consider when swapping orbitals
  int LzSymmetryMaxSwapPosition;

 public:

  // default constructor
  // 
  BosonOnSphereWithSU2SpinLzSymmetry ();

  // basic constructor without any constraint on Lz
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
  BosonOnSphereWithSU2SpinLzSymmetry (int nbrBosons, int lzMax, bool minusLzParity);

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // totalSpin = twice the total spin value (not taken into account)
  // minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU2SpinLzSymmetry (int nbrBosons, int lzMax, int totalSpin, bool minusLzParity, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  // 
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU2SpinLzSymmetry (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereWithSU2SpinLzSymmetry(const BosonOnSphereWithSU2SpinLzSymmetry& bosons);

  // destructor
  //
  virtual ~BosonOnSphereWithSU2SpinLzSymmetry ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereWithSU2SpinLzSymmetry& operator = (const BosonOnSphereWithSU2SpinLzSymmetry& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // convert a given state from a generic basis to the current Sz subspace basis
  //
  // state = reference on the vector to convert
  // space = reference on the basis associated to state
  // return value = converted vector
  virtual RealVector ConvertFromNbodyBasis(RealVector& state, ParticleOnSphereWithSpin* space);
  
  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

 protected:

  // read Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description is stored
  // return value = true if no error occured
  virtual bool ReadHilbertSpace (char* fileName);

  // generate the Hilbert space with the discrete symmetry constraint
  //
  virtual void GenerateStatesWithDiscreteSymmetry();

  // factorized code that is used to symmetrize the result of any operator action
  //
  // stateDescriptionUp = reference on the state up part that has been produced with the operator action
  // stateDescriptionDown = reference on the state down part that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown, double& coefficient);

  // find the canonical form of a given state
  //
  // stateDescriptionUp = unsigned integer describing the state up spin component
  // stateDescriptionDown = unsigned integer describing the state down spin component
  // canonicalStateDescriptionUp = reference on the canonical state up spin component
  // canonicalStateDescriptionDown = reference on the canonical state down spin component
  // nbrLzSymmetry = reference on the number of Lz<->-Lz flip that has to be applied to obtained the canonical form
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual void FindCanonicalForm(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown,
				 unsigned long& canonicalStateDescriptionUp, unsigned long& canonicalStateDescriptionDown, int& nbrLzSymmetry);

  // find the size of the orbit for a given state
  //
  // stateDescriptionUp = unsigned integer describing the state up spin component
  // stateDescriptionDown = unsigned integer describing the state down spin component
  // return value = orbit size
  virtual int FindOrbitSize(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown);
  
  //  test if the state is allowed by the the symmetry constraint
  //
  // stateDescriptionUp = unsigned integer describing the state up spin component
  // stateDescriptionDown = unsigned integer describing the state down spin component
  // return value = true if the state satisfies the symmetry constraint
  virtual bool TestSymmetryConstraint(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown);

  // Apply the Lz<->-Lz operator to a single spin component
  //
  // stateDescription = reference on state description
  // stateLzMax = maximum lz value reached by the fermionic state
  virtual void ApplyLzSymmetry (unsigned long& stateDescription, int stateLzMax);

  // Apply the Lz<->-Lz operator to both spin components
  //
  // stateDescriptionUp = reference on up spin state description
  // stateDescriptionDown = reference on down spin state description
  virtual void ApplyLzSymmetry (unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown);
  
};


// factorized code that is used to symmetrize the result of any operator action
//
// stateDescriptionUp = reference on the state up part that has been produced with the operator action
// stateDescriptionDown = reference on the state down part that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state  

inline int BosonOnSphereWithSU2SpinLzSymmetry::SymmetrizeAdAdResult(unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown, double& coefficient)
{
  unsigned long TmpCanonicalStateDescriptionUp;
  unsigned long TmpCanonicalStateDescriptionDown;
  int NbrLzSymmetry;
  this->FindCanonicalForm(stateDescriptionUp, stateDescriptionDown, TmpCanonicalStateDescriptionUp, TmpCanonicalStateDescriptionDown, NbrLzSymmetry);
  int TmpIndex = this->FindStateIndex(TmpCanonicalStateDescriptionUp, TmpCanonicalStateDescriptionDown);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      if (NbrLzSymmetry == 1) 
	coefficient *= this->LzParitySign;  
    }
  return TmpIndex;
}

// find the canonical form of a given state
//
// stateDescriptionUp = unsigned integer describing the state up spin component
// stateDescriptionDown = unsigned integer describing the state down spin component
// canonicalStateDescriptionUp = reference on the canonical state up spin component
// canonicalStateDescriptionDown = reference on the canonical state down spin component
// nbrLzSymmetry = reference on the number of Lz<->-Lz flip that has to be applied to obtained the canonical form
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline void BosonOnSphereWithSU2SpinLzSymmetry::FindCanonicalForm(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown,
								  unsigned long& canonicalStateDescriptionUp, unsigned long& canonicalStateDescriptionDown, 
								  int& nbrLzSymmetry)
{
  nbrLzSymmetry = 0;
  canonicalStateDescriptionUp = stateDescriptionUp;
  this->ApplyLzSymmetry(canonicalStateDescriptionUp, this->NUpLzMax);
  if (stateDescriptionUp > canonicalStateDescriptionUp)
    {
      canonicalStateDescriptionUp = stateDescriptionUp;
      canonicalStateDescriptionDown = stateDescriptionDown;
      return;
    }
  canonicalStateDescriptionDown = stateDescriptionDown;
  this->ApplyLzSymmetry(canonicalStateDescriptionDown, this->NDownLzMax);
  if (stateDescriptionUp == canonicalStateDescriptionUp)
    {
      if (stateDescriptionDown > canonicalStateDescriptionDown)
	{
	  canonicalStateDescriptionUp = stateDescriptionUp;
	  canonicalStateDescriptionDown = stateDescriptionDown;
	  return;
	}
    } 
  nbrLzSymmetry = 1;
}

//  test if the state is allowed by the the symmetry constraint
//
// stateDescriptionUp = unsigned integer describing the state up spin component
// stateDescriptionDown = unsigned integer describing the state down spin component
// return value = true if the state satisfies the symmetry constraint

inline bool BosonOnSphereWithSU2SpinLzSymmetry::TestSymmetryConstraint(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
{
  unsigned long CanonicalStateDescriptionUp = stateDescriptionUp;
  this->ApplyLzSymmetry(CanonicalStateDescriptionUp, this->NUpLzMax);
  if (stateDescriptionUp > CanonicalStateDescriptionUp)
    {
      return true;
    }
  unsigned long CanonicalStateDescriptionDown = stateDescriptionDown;
  this->ApplyLzSymmetry(CanonicalStateDescriptionDown, this->NDownLzMax);
  if ((stateDescriptionUp == CanonicalStateDescriptionUp) && (stateDescriptionDown == CanonicalStateDescriptionDown))
    {
      if (this->LzParitySign > 0.0)
	{
	  return true;
	}
      else
	{
	  return false;
	}
    }
  return true;
}

// find the size of the orbit for a given state
//
// stateDescriptionUp = unsigned integer describing the state up spin component
// stateDescriptionDown = unsigned integer describing the state down spin component
// return value = orbit size

inline int BosonOnSphereWithSU2SpinLzSymmetry::FindOrbitSize(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
{
  unsigned long CanonicalStateDescriptionUp = stateDescriptionUp;
  this->ApplyLzSymmetry(CanonicalStateDescriptionUp, this->NUpLzMax);
  if (stateDescriptionUp > CanonicalStateDescriptionUp)
    {
      return 2;
    }
  unsigned long CanonicalStateDescriptionDown = stateDescriptionDown;
  this->ApplyLzSymmetry(CanonicalStateDescriptionDown, this->NDownLzMax);
  if ((stateDescriptionUp == CanonicalStateDescriptionUp) && (stateDescriptionDown == CanonicalStateDescriptionDown))
    {
      return 1;
    }
  return 2;
}

// Apply the Lz<->-Lz operator to a single spin component
//
// stateDescription = reference on state description
// stateLzMax = maximum lz value reached by the fermionic state

inline void BosonOnSphereWithSU2SpinLzSymmetry::ApplyLzSymmetry (unsigned long& stateDescription, int stateLzMax)
{
#ifdef __64_BITS__
  unsigned long TmpState = FermionOnSphereInvertTable[stateDescription & 0xff] << 56;
  TmpState |= FermionOnSphereInvertTable[(stateDescription >> 8) & 0xff] << 48;
  TmpState |= FermionOnSphereInvertTable[(stateDescription >> 16) & 0xff] << 40;
  TmpState |= FermionOnSphereInvertTable[(stateDescription >> 24) & 0xff] << 32;
  TmpState |= FermionOnSphereInvertTable[(stateDescription >> 32) & 0xff] << 24;
  TmpState |= FermionOnSphereInvertTable[(stateDescription >> 40) & 0xff] << 16;
  TmpState |= FermionOnSphereInvertTable[(stateDescription >> 48) & 0xff] << 8;
  TmpState |= FermionOnSphereInvertTable[stateDescription >> 56]; 
  TmpState >>= 63 - stateLzMax;
#else
  unsigned long TmpState = FermionOnSphereInvertTable[stateDescription & 0xff] << 24;
  TmpState |= FermionOnSphereInvertTable[(stateDescription >> 8) & 0xff] << 16;
  TmpState |= FermionOnSphereInvertTable[(stateDescription >> 16) & 0xff] << 8;
  TmpState |= FermionOnSphereInvertTable[stateDescription >> 24];
  TmpState >>= 31 - stateLzMax;
#endif	
  stateDescription = TmpState;
/*   this->FermionToBoson(stateDescription, stateLzMax, this->TemporaryStateUp); */
/*   for (int i = 0; i <= this->LzSymmetryMaxSwapPosition; ++i) */
/*     { */
/*       unsigned long Tmp = this->TemporaryStateUp[i]; */
/*       this->TemporaryStateUp[i] = this->TemporaryStateUp[this->LzMax - i]; */
/*       this->TemporaryStateUp[this->LzMax - i] = Tmp; */
/*     } */
/*   stateDescription = this->BosonToFermion(this->TemporaryStateUp); */
}

// Apply the Lz<->-Lz operator to both spin components
//
// stateDescriptionUp = reference on up spin state description
// stateDescriptionDown = reference on down spin state description

inline void BosonOnSphereWithSU2SpinLzSymmetry::ApplyLzSymmetry (unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown)
{
  this->FermionToBoson(stateDescriptionUp, this->NUpLzMax, this->TemporaryStateUp);
  this->FermionToBoson(stateDescriptionDown, this->NDownLzMax, this->TemporaryStateDown);
  for (int i = 0; i <= this->LzSymmetryMaxSwapPosition; ++i)
    {
      unsigned long Tmp = this->TemporaryStateUp[i];
      this->TemporaryStateUp[i] = this->TemporaryStateUp[this->LzMax - i];
      this->TemporaryStateUp[this->LzMax - i] = Tmp;
      Tmp = this->TemporaryStateDown[i];
      this->TemporaryStateDown[i] = this->TemporaryStateDown[this->LzMax - i];
      this->TemporaryStateDown[this->LzMax - i] = Tmp;
    }
  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, stateDescriptionUp, stateDescriptionDown);
}

#endif


