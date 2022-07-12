////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(4) spin         //
//                      in momentum space with Sz<->-Sz symmetry              //
//                                                                            //
//                        last modification : 16/06/2020                      //
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


#ifndef FERMIONONSQUARELATTICEWITHSU4SPINMOMENTUMSPACESZSYMMETRY_H
#define FERMIONONSQUARELATTICEWITHSU4SPINMOMENTUMSPACESZSYMMETRY_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"


#include <iostream>



class FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry : public FermionOnSquareLatticeWithSU4SpinMomentumSpace
{

 protected:

  // additional sign due to the parity sector for the Sz<->-Sz symmetry
  double SzParitySign;

public:

  // default constructor
  //
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, bool minusSzParity, unsigned long memory = 10000000);

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
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSpin, bool minusSzParity, unsigned long memory = 10000000);

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
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSpin, int totalIsospin, bool minusSzParity, unsigned long memory = 10000000);

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
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSpin, int totalIsospin,
							    int totalEntanglement, bool minusSzParity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry(const FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry& operator = (const FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry& fermions);

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
  virtual double GetSzSymmetrySign (unsigned long initialState);

  // Apply the Sz operator to flip all the spins
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = flipped configuration
  virtual unsigned long ApplySzSymmetry (unsigned long initialState);
  
};

// get the fermionic sign when applying the Sz flip operator
//
// initialState = state that has to be converted to its canonical expression
// return value = fermionic sign

inline double FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry::GetSzSymmetrySign (unsigned long initialState)
{
  initialState = (initialState >> 1) ^ initialState;
  initialState &= (initialState >> 2);
  initialState &= 0x1111111111111111ul;
#ifdef __64_BITS__
  initialState ^= (initialState >> 32);
#endif
  initialState ^= (initialState >> 16);
  initialState ^= (initialState >> 8);
  initialState ^= (initialState >> 4);
  return (((double) ((initialState & 0x1ul) << 1)) - 1.0);
}

// apply the Sz operator to flip all the spins
//
// initialState = state that has to be flipped
// return value = fermionic sign

inline unsigned long FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry::ApplySzSymmetry (unsigned long initialState)
{
  return (((initialState & 0x3333333333333333ul) << 2) | ((initialState & 0xccccccccccccccccul) >> 2));
}


// factorized code that is used to compute any symmetry information of a given state
//
// state = reference on the state that has been produced with the operator action

inline void FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry::SymmetrizeAAInput(unsigned long& state)
{
  if (state != this->ApplySzSymmetry(state))
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

inline int FermionOnSquareLatticeWithSU4SpinMomentumSpaceSzSymmetry::SymmetrizeAdAdResult(unsigned long& state, double& coefficient, int highestBit)
{
  unsigned long TmpState = this->ApplySzSymmetry(state);
  coefficient *= this->ProdASymmetryFactor;
  if (state == TmpState)
    {
      if (this->SzParitySign != this->GetSzSymmetrySign(state))
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
	  coefficient *= this->SzParitySign;
	  coefficient *= this->GetSzSymmetrySign(state);
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


