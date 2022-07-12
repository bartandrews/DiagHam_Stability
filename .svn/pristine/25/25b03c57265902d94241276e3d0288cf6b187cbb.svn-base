////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                class of bosons on a square lattice with SU(3) spin         //
//                                in momentum space                           //
//                                                                            //
//                        last modification : 06/12/2011                      //
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


#ifndef BOSONONSQUARELATTICEWITHSU3SPINMOMENTUMSPACE_H
#define BOSONONSQUARELATTICEWITHSU3SPINMOMENTUMSPACE_H

#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSU3Spin.h"

#include <iostream>



class BosonOnSquareLatticeWithSU3SpinMomentumSpace : public BosonOnSphereWithSU3Spin
{

 protected:

  // number of sites in the x direction
  int NbrSiteX;
  // number of sites in the y direction
  int NbrSiteY;

  // momentum along the x direction
  int KxMomentum;
  // momentum along the y direction
  int KyMomentum;

 public:

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  BosonOnSquareLatticeWithSU3SpinMomentumSpace (int nbrBosons, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSquareLatticeWithSU3SpinMomentumSpace(const BosonOnSquareLatticeWithSU3SpinMomentumSpace& bosons);

  // destructor
  //
  ~BosonOnSquareLatticeWithSU3SpinMomentumSpace ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSquareLatticeWithSU3SpinMomentumSpace& operator = (const BosonOnSquareLatticeWithSU3SpinMomentumSpace& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentFermionicPosition1 = current fermionic position within the state description for the type 1 particles
  // currentFermionicPosition2 = current fermionic position within the state description for the type 2 particles
  // currentFermionicPosition3 = current fermionic position within the state description for the type 3 particles
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, 
			      int currentFermionicPosition1, int currentFermionicPosition2, int currentFermionicPosition3, long pos);


};


#endif


