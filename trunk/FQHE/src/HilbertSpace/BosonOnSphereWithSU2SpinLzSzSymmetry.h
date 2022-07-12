////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//           class for bosons on sphere with SU(2) spin and both the          //
//                       Lz<->-Lz and Sz<->-Sz symmetries                     //
//                                                                            //
//                        last modification : 25/09/2016                      //
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


#ifndef BOSONONSPHEREWITHSU2SPINLZSZSYMMETRY_H
#define BOSONONSPHEREWITHSU2SPINLZSZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSymmetry.h"

#include <iostream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;



class BosonOnSphereWithSU2SpinLzSzSymmetry :  public BosonOnSphereWithSU2SpinLzSymmetry
{


 protected:


 public:

  // default constructor
  // 
  BosonOnSphereWithSU2SpinLzSzSymmetry ();

  // basic constructor without any constraint on Lz
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
  BosonOnSphereWithSU2SpinLzSzSymmetry (int nbrBosons, int lzMax, bool minusSzParity, bool minusLzParity);

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // totalSpin = twice the total spin value (not taken into account)
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU2SpinLzSzSymmetry (int nbrBosons, int lzMax, int totalSpin, bool minusSzParity, bool minusLzParity, unsigned long memory = 10000000);

  // constructor from a binary file that describes the Hilbert space
  // 
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU2SpinLzSzSymmetry (char* fileName, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereWithSU2SpinLzSzSymmetry(const BosonOnSphereWithSU2SpinLzSzSymmetry& bosons);

  // destructor
  //
  virtual ~BosonOnSphereWithSU2SpinLzSzSymmetry ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereWithSU2SpinLzSzSymmetry& operator = (const BosonOnSphereWithSU2SpinLzSzSymmetry& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

 protected:

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
  // nbrSzSymmetry = reference on the number of Sz<->-Sz flips that has to be applied to obtained the canonical form
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual void FindCanonicalForm(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown,
				 unsigned long& canonicalStateDescriptionUp, unsigned long& canonicalStateDescriptionDown, 
				 int& nbrLzSymmetry, int& nbrSzSymmetry);

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

};


// factorized code that is used to symmetrize the result of any operator action
//
// stateDescriptionUp = reference on the state up part that has been produced with the operator action
// stateDescriptionDown = reference on the state down part that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state  

inline int BosonOnSphereWithSU2SpinLzSzSymmetry::SymmetrizeAdAdResult(unsigned long& stateDescriptionUp, unsigned long& stateDescriptionDown, double& coefficient)
{
  unsigned long TmpCanonicalStateDescriptionUp;
  unsigned long TmpCanonicalStateDescriptionDown;
  int NbrLzSymmetry;
  int NbrSzSymmetry;
  this->FindCanonicalForm(stateDescriptionUp, stateDescriptionDown, TmpCanonicalStateDescriptionUp, TmpCanonicalStateDescriptionDown, NbrLzSymmetry, NbrSzSymmetry);
  int TmpIndex = this->FindStateIndex(TmpCanonicalStateDescriptionUp, TmpCanonicalStateDescriptionDown);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->ProdATemporaryNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      if (NbrLzSymmetry == 1) 
	coefficient *= this->LzParitySign;  
      if (NbrSzSymmetry == 1) 
	coefficient *= this->SzParitySign;  
    }
  return TmpIndex;
}

// find the canonical form of a given state
//
// stateDescriptionUp = unsigned integer describing the state up spin component
// stateDescriptionDown = unsigned integer describing the state down spin component
// canonicalStateDescriptionUp = reference on the canonical state up spin component
// canonicalStateDescriptionDown = reference on the canonical state down spin component
// nbrLzSymmetry = reference on the number of Lz<->-Lz flips that has to be applied to obtained the canonical form
// nbrSzSymmetry = reference on the number of Sz<->-Sz flips that has to be applied to obtained the canonical form
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline void BosonOnSphereWithSU2SpinLzSzSymmetry::FindCanonicalForm(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown,
								    unsigned long& canonicalStateDescriptionUp, unsigned long& canonicalStateDescriptionDown, 
								    int& nbrLzSymmetry, int& nbrSzSymmetry)
{
  unsigned long TmpStateDescriptionUp = stateDescriptionUp;
  unsigned long TmpStateDescriptionDown = stateDescriptionDown;  
  this->ApplyLzSymmetry(TmpStateDescriptionUp, this->NUpLzMax);
  this->ApplyLzSymmetry(TmpStateDescriptionDown, this->NDownLzMax);

  if ((stateDescriptionUp > TmpStateDescriptionUp) || 
      ((stateDescriptionUp == TmpStateDescriptionUp) && (stateDescriptionDown >= TmpStateDescriptionDown)))
    {
      nbrLzSymmetry = 0;
      nbrSzSymmetry = 0;
      canonicalStateDescriptionUp = stateDescriptionUp;
      canonicalStateDescriptionDown = stateDescriptionDown;
    }
  else
    {
      canonicalStateDescriptionUp = TmpStateDescriptionUp;
      canonicalStateDescriptionDown = TmpStateDescriptionDown;
      nbrSzSymmetry = 0;
      nbrLzSymmetry = 1;
    }

  if ((stateDescriptionDown > canonicalStateDescriptionUp) || 
      ((stateDescriptionDown == canonicalStateDescriptionUp) && (stateDescriptionUp > canonicalStateDescriptionDown)))
    {
      canonicalStateDescriptionUp = stateDescriptionDown;
      canonicalStateDescriptionDown = stateDescriptionUp;
      nbrSzSymmetry = 1;
      nbrLzSymmetry = 0;
    }

  if ((TmpStateDescriptionDown > canonicalStateDescriptionUp) || 
      ((TmpStateDescriptionDown == canonicalStateDescriptionUp) && (TmpStateDescriptionUp > canonicalStateDescriptionDown)))
    {
      canonicalStateDescriptionUp = TmpStateDescriptionDown;
      canonicalStateDescriptionDown = TmpStateDescriptionUp;
      nbrSzSymmetry = 1;
      nbrLzSymmetry = 1;
    }
}

//  test if the state is allowed by the the symmetry constraint
//
// stateDescriptionUp = unsigned integer describing the state up spin component
// stateDescriptionDown = unsigned integer describing the state down spin component
// return value = true if the state satisfies the symmetry constraint

inline bool BosonOnSphereWithSU2SpinLzSzSymmetry::TestSymmetryConstraint(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
{
  unsigned long TmpStateDescriptionUp = stateDescriptionUp;
  unsigned long TmpStateDescriptionDown = stateDescriptionDown;  
  this->ApplyLzSymmetry(TmpStateDescriptionUp, this->NUpLzMax);
  this->ApplyLzSymmetry(TmpStateDescriptionDown, this->NDownLzMax);
  if ((stateDescriptionUp == TmpStateDescriptionUp) && (stateDescriptionDown == TmpStateDescriptionDown))
    {
     if (this->LzParitySign < 0.0)
	{
	  return false;
	}
      if ((TmpStateDescriptionUp == TmpStateDescriptionDown) && (this->SzParitySign < 0.0))
	{
	  return false;
	}
      return true;
    }
  if (TmpStateDescriptionUp == TmpStateDescriptionDown)  
    {
      if (this->SzParitySign < 0.0)
	{
	  return false;
	}
      else
	{
	  return true;
	}
    }
  if ((stateDescriptionUp == TmpStateDescriptionDown) && (stateDescriptionDown == TmpStateDescriptionUp))
    {
      return (this->SzParitySign == this->LzParitySign);
    }
  return true;
}

// find the size of the orbit for a given state
//
// stateDescriptionUp = unsigned integer describing the state up spin component
// stateDescriptionDown = unsigned integer describing the state down spin component
// return value = orbit size

inline int BosonOnSphereWithSU2SpinLzSzSymmetry::FindOrbitSize(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
{
  unsigned long CanonicalStateDescriptionUp = stateDescriptionUp;
  unsigned long CanonicalStateDescriptionDown = stateDescriptionDown;
  this->ApplyLzSymmetry(CanonicalStateDescriptionUp, this->NUpLzMax);
  this->ApplyLzSymmetry(CanonicalStateDescriptionDown, this->NDownLzMax);
  if ((stateDescriptionUp == CanonicalStateDescriptionUp) && (stateDescriptionDown == CanonicalStateDescriptionDown))
    {
      if (stateDescriptionUp == stateDescriptionDown)
	{
	  return 1;
	}
      else
	{
	  return 2;
	}
    }
  else
    {
      if (stateDescriptionUp == stateDescriptionDown)
	{
	  return 2;
	}
      else
	{
	  if ((stateDescriptionUp == CanonicalStateDescriptionDown) && (stateDescriptionDown == CanonicalStateDescriptionUp))
	    {
	      return 2;
	    }
	  else
	    {
	      return 4;
	    }
	}
    }
}


#endif


