////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of fermions on a square lattice with SU(4) spin           //
//       in momentum space with the Z2 symmetry swapping simultaneously       //
//                    Sz<->-Sz and Pz<->-Pz while preserving Ez               //
//                                                                            //
//                        last modification : 17/06/2020                      //
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


#ifndef FERMIONONSQUARELATTICEWITHSU4SPINMOMENTUMSPACESZPZPRESERVINGEZSYMMETRY_H
#define FERMIONONSQUARELATTICEWITHSU4SPINMOMENTUMSPACESZPZPRESERVINGEZSYMMETRY_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry.h"


#include <iostream>

static double FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[] = {1.0, 1.0, 1.0, -1.0,
									    1.0, -1.0, -1.0, -1.0,
									    1.0, -1.0, -1.0, -1.0,
									    1.0, 1.0, -1.0, 1.0};


class FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry : public FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry
{

 protected:

  // additional sign due to the parity sector for the Z2 symmetry swapping simultaneously Sz<->-Sz and Pz<->-Pz
  double SzPzParitySign;

public:

  // default constructor
  //
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // minusSzPzParity = select the Z2 symmetry swapping simultaneously Sz <-> -Sz and Pz<->-Pz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, bool minusSzPzParity, unsigned long memory = 10000000);

  // constructor when preserving only spin
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // totalSpin = twice the total spin value
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSpin, bool minusSzParity, unsigned long memory = 10000000);

  // constructor when preserving spin and isospin
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // totalSpin = twice the total spin value
  // totalIsospin = twice the total isospin value
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSpin, int totalIsospin, bool minusSzParity, unsigned long memory = 10000000);

  // constructor when preserving the three Cartan quantum numbers
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // totalSpin = twice the total spin value
  // totalIsospin = twice the total isospin value
  // totalEntanglement = twice the total entanglement value
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSpin, int totalIsospin,
							    int totalEntanglement, bool minusSzParity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry(const FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry& operator = (const FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:

  // factorized code that is used to compute any symmetry information of a given state
  //
  // state = reference on the state that has been produced with the operator action
  virtual void SymmetrizeAAInput(unsigned long& state);
  
  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // highestBit = highest bit set to one in state
  // return value = index of the destination state  
  virtual int SymmetrizeAdAdResult(unsigned long& state, double& coefficient, int highestBit);

  // generate all states corresponding to the constraints and discrete symmetries
  // 
  // return value = Hilbert space dimension
  virtual long GenerateStatesWithSymmetries();
  
  // get the fermionic sign when applying the Sz flip operator
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = fermionic sign
  virtual double GetSzPzSymmetrySign (unsigned long initialState);

  // Apply the Sz operator to flip all the spins
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = flipped configuration
  virtual unsigned long ApplySzPzSymmetry (unsigned long initialState);
  
};

// get the fermionic sign when applying the Sz flip operator
//
// initialState = state that has to be converted to its canonical expression
// return value = fermionic sign

inline double FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry::GetSzPzSymmetrySign (unsigned long initialState)
{
#ifdef __64_BITS__
  return (FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 60) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 56) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 52) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 48) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 44) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 40) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 36) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 32) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 28) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 24) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 20) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 16) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 12) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 8) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 4) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[initialState & 0xful]);
#else 
  return (FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 28) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 24) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 20) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 16) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 12) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 8) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[(initialState >> 4) & 0xful]
	  * FermionWithSU4SpinSzPzPreservingEzSymmetryReorderingSign[initialState & 0xful]);
#endif
}

// apply the Sz operator to flip all the spins
//
// initialState = state that has to be flipped
// return value = fermionic sign

inline unsigned long FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry::ApplySzPzSymmetry (unsigned long initialState)
{
  return (((initialState & 0x1111111111111111ul) << 3) | ((initialState & 0x2222222222222222ul) << 1)
	  | ((initialState & 0x4444444444444444ul) >> 1) | ((initialState & 0x8888888888888888ul) >> 3));
}


// factorized code that is used to compute any symmetry information of a given state
//
// state = reference on the state that has been produced with the operator action

inline void FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry::SymmetrizeAAInput(unsigned long& state)
{
  if (state != this->ApplySzPzSymmetry(state))
    {
      this->ProdASymmetryFactor = M_SQRT2;
    }
  else
    {
      this->ProdASymmetryFactor = 1.0;
    }
}
  

// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// highestBit = highest bit set to one in state
// return value = index of the destination state  

inline int FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzPzPreservingEzSymmetry::SymmetrizeAdAdResult(unsigned long& state, double& coefficient, int highestBit)
{
  unsigned long TmpState = this->ApplySzPzSymmetry(state);
  coefficient *= this->ProdASymmetryFactor;
  if (state == TmpState)
    {
      if (this->SzPzParitySign != this->GetSzPzSymmetrySign(state))
	{
	  return this->HilbertSpaceDimension;
	}
    }
  else
    {
      coefficient *= M_SQRT1_2;
      if (state > TmpState)
	{
	  TmpState = state;
	}
      else
	{
	  coefficient *= this->SzPzParitySign;
	  coefficient *= this->GetSzPzSymmetrySign(state);
	}
    }
  highestBit = this->HighestBit;
  while (((TmpState >> highestBit) & 0x1ul) == 0x0ul)
    {
      --highestBit;
    }
  return this->FindStateIndex(TmpState, highestBit);
}


#endif


