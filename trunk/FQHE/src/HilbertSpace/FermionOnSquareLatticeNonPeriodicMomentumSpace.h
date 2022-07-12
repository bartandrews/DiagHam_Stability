////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of fermions on square lattice                 //
//                  in momentum space with open boundary conditions           //
//                                                                            //
//                        last modification : 28/02/2011                      //
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


#ifndef FERMIONONSQUARELATTICENONPERIODICMOMENTUMSPACE_H
#define FERMIONONSQUARELATTICENONPERIODICMOMENTUMSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

#include <iostream>



class FermionOnSquareLatticeNonPeriodicMomentumSpace : public FermionOnSquareLatticeMomentumSpace
{

  friend class FermionOnSquareLatticeMomentumSpace;

 protected:

  // number of x momenta allowe for a single particle
  int NbrAllowedKx;
  // number of y momenta allowe for a single particle
  int NbrAllowedKy;
  // minimal x momentum allowed for a single particle
  int MinKx;
  // minimal y momentum allowed for a single particle
  int MinKy;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrAllowedKx = number of kx momenta allowed in the space
  // nbrAllowedKy = number of kx momenta allowed in the space
  // minKx = minimum value of kx momenta allowed in the space
  // minKy = minimum value of ky momenta allowed in the space
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeNonPeriodicMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int nbrAllowedKx, int nbrAllowedKy, int minKx, int minKy,
                                                  int kxMomentum, int kyMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeNonPeriodicMomentumSpace(const FermionOnSquareLatticeNonPeriodicMomentumSpace& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeNonPeriodicMomentumSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeNonPeriodicMomentumSpace& operator = (const FermionOnSquareLatticeNonPeriodicMomentumSpace& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKxShifted = current momentum along x for a single particle, shifted to fit in [0,NbrAllowedKx)
  // currentKyShifted = current momentum along y for a single particle, shifted to fit in [0,NbrAllowedKy)
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKxShifted, int currentKyShifted, int currentTotalKx, int currentTotalKy);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKxShifted = current momentum along x for a single particle, shifted to fit in [0,NbrAllowedKx)
  // currentKyShifted = current momentum along y for a single particle, shifted to fit in [0,NbrAllowedKy)
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKxShifted, int currentKyShifted, int currentTotalKx, int currentTotalKy, long pos);

};


#endif


