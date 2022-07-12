////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of spin 1 chain with the Sz<->Sz and                  //
//                              inversion symmetries                          //
//                                                                            //
//                        last modification : 11/05/2017                      //
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


#ifndef SPIN1CHAINWITHSZINVERSIONSYMMETRIES_H
#define SPIN1CHAINWITHSZINVERSIONSYMMETRIES_H


#include "config.h"
#include "HilbertSpace/Spin1ChainWithInversionSymmetry.h"
#include "Vector/RealVector.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class ComplexMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;



class Spin1ChainWithSzInversionSymmetries : public Spin1ChainWithInversionSymmetry
{

  friend class Spin1ChainWithTranslations;
  friend class Spin1ChainWithTranslationsAndSzSymmetry;
  friend class Spin1ChainWithTranslationsAndInversionSymmetry;

 protected:

 public:

  // default constructor
  //
  Spin1ChainWithSzInversionSymmetries ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // inversionSector = inversion symmetry sector (can be either +1 or -1)
  // szSymmetrySector = Sz<->-Sz symmetry sector (can be either +1 or -1)
  // memorySize = memory size in bytes allowed for look-up table
  Spin1ChainWithSzInversionSymmetries (int chainLength, int inversionSector, int szSymmetrySector, unsigned long memorySize);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // inversionSector = inversion symmetry sector (can be either +1 or -1)
  // szSymmetrySector = Sz<->-Sz symmetry sector (can be either +1 or -1)
  // sz = twice the value of total Sz component (should be equal to zero)
  // memorySize = memory size in bytes allowed for look-up table
  Spin1ChainWithSzInversionSymmetries (int chainLength, int inversionSector, int szSymmetrySector, int sz, unsigned long memory) ;

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1ChainWithSzInversionSymmetries (const Spin1ChainWithSzInversionSymmetries& chain);

  // destructor
  //
  ~Spin1ChainWithSzInversionSymmetries ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1ChainWithSzInversionSymmetries& operator = (const Spin1ChainWithSzInversionSymmetries& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

 protected:

  // find canonical form of a state description
  //
  // stateDescription = unsigned integer describing the state
  // szSymmetrySign = reference on the additional sign coming from the Sz<->-Sz symmetry
  // return value = canonical form of a state description 
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, double& szSymmetrySign);

  //  test if the state is compatible with the discrete symmetry sector
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the discrete symmetry sector
  virtual bool TestDiscreteSymmetryConstraint(unsigned long stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(unsigned long stateDescription);


};

// find canonical form of a state description 
//
// stateDescription = unsigned integer describing the state
// szSymmetrySign = reference on the additional sign coming from the Sz<->-Sz symmetry
// return value = canonical form of a state description

inline unsigned long Spin1ChainWithSzInversionSymmetries::FindCanonicalForm(unsigned long stateDescription, double& szSymmetrySign)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long TmpStateDescription = stateDescription;
  this->ApplyInversionSymmetry(TmpStateDescription);  
  szSymmetrySign = 1.0;
  if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      szSymmetrySign = this->InversionSector;      
    }
  this->ApplySzSymmetry(stateDescription);  
  if (stateDescription < CanonicalState)
    {
      CanonicalState = stateDescription;
      szSymmetrySign = this->SzSymmetrySector;      
    }
  this->ApplySzSymmetry(TmpStateDescription);  
  if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      szSymmetrySign = this->SzSymmetrySector * this->InversionSector;      
    } 
  return CanonicalState;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int Spin1ChainWithSzInversionSymmetries::FindOrbitSize(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  this->ApplyInversionSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    {
      this->ApplySzSymmetry(TmpStateDescription);
      if (stateDescription == TmpStateDescription)
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
      unsigned long TmpStateDescription2 = stateDescription;
      this->ApplySzSymmetry(TmpStateDescription2);
      if ((stateDescription == TmpStateDescription2) || (TmpStateDescription == TmpStateDescription2))
	{
	  return 2;
	}
      else
	{
	  return 4;
	}
    }
}

//  test if the state is compatible with the discrete symmetry sector
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the discrete symmetry sector

inline bool Spin1ChainWithSzInversionSymmetries::TestDiscreteSymmetryConstraint(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  this->ApplyInversionSymmetry(TmpStateDescription);
  if ((stateDescription == TmpStateDescription) && (this->InversionSector < 0.0))
    return false;
  unsigned long TmpStateDescription2 = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription2);
  if ((stateDescription == TmpStateDescription2) && (this->SzSymmetrySector < 0.0))
    return false;
  if ((TmpStateDescription2 == TmpStateDescription) && ((this->InversionSector * this->SzSymmetrySector) < 0.0))
    return false;
  else
    return true;
}

#endif


