////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of fermions on lattice with spin  and Gutzwiller          //
//             projection in real space with the Sz<->-Sz symmetry            //
//                                                                            //
//                        last modification : 18/08/2014                      //
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


#ifndef FERMIONONLATTICEWITHSPINSZSYMMETRYANDGUTZWILLERPROJECTIONREALSPACE_H
#define FERMIONONLATTICEWITHSPINSZSYMMETRYANDGUTZWILLERPROJECTIONREALSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpace.h"

#include <iostream>



class FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace : public FermionOnLatticeWithSpinSzSymmetryRealSpace
{

  friend class FermionOnSquareLatticeWithSU4SpinMomentumSpace;

 protected:


 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace (int nbrFermions, int nbrSite, bool minusSzParity, unsigned long memory = 10000000);

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalSpin = twice the Sz projection
  // nbrSite = number of sites
  // minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace (int nbrFermions, int totalSpin, int nbrSite, bool minusSzParity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace(const FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace& fermions);

  // destructor
  //
  ~FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace& operator = (const FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();
  

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions);

  // evaluate Hilbert space dimension with a fixed number of fermions with spin up
  //
  // nbrFermions = number of fermions
  // nbrSpinUp = number of fermions with spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // nbrHoles = number of unoccupied sites
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, int nbrHoles, long pos);
  
  // generate all states corresponding to the constraints with a fixed number of fermions with spin up
  // 
  // nbrFermions = number of fermions
  // currentSite = current site index in real state
  // nbrHoles = number of unoccupied sites
  //nbrSpinUp = number of fermions with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentSite, int nbrHoles, int nbrSpinUp, long pos);
  
  // find state index (and checks state belongs to Hilbert space)
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);
  
  // find state index and checks that states belongs to Hilbert space
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int CarefulFindStateIndex(unsigned long stateDescription, int lzmax);



};


// find state index and checks that states belongs to Hilbert space
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index
inline int FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::CarefulFindStateIndex(unsigned long stateDescription, int lzmax)
{
 return this->FindStateIndex(stateDescription, lzmax); 
}

#endif


