////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on hyper cubic 4D lattice with spin         //
//                                  in momentum space                         //
//                                                                            //
//                        last modification : 12/07/2011                      //
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


#ifndef FERMIONONHYPERCUBICLATTICEWITHSPINMOMENTUMSPACE_H
#define FERMIONONHYPERCUBICLATTICEWITHSPINMOMENTUMSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

#include <iostream>



class FermionOnHyperCubicLatticeWithSpinMomentumSpace : public FermionOnSphereWithSpin
{

 protected:

  // number of sites in the x direction
  int NbrSiteX;
  // number of sites in the y direction
  int NbrSiteY;
  // number of sites in the z direction
  int NbrSiteZ;
  // number of sites in the t direction
  int NbrSiteT;
  // number of sites in the direction perpendicular to X
  int NbrSiteYZT;
  // number of sites in the direction perpendicular to XY
  int NbrSiteZT;

  // momentum along the x direction
  int KxMomentum;
  // momentum along the y direction
  int KyMomentum;
  // momentum along the z direction
  int KzMomentum;
  // momentum along the t direction
  int KtMomentum;

  // flag to indicate that the Hilbert space should preserve Sz
  bool SzFlag;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // nbrSiteT = number of sites in the t direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // kzMomentum = momentum along the z direction
  // ktMomentum = momentum along the t direction
  // memory = amount of memory granted for precalculations
  FermionOnHyperCubicLatticeWithSpinMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT, int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, unsigned long memory = 10000000);

  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // nbrSpinUp = number of particles with spin up
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // nbrSiteT = number of sites in the t direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // kzMomentum = momentum along the z direction
  // ktMomentum = momentum along the t direction
  // memory = amount of memory granted for precalculations
  FermionOnHyperCubicLatticeWithSpinMomentumSpace (int nbrFermions, int nbrSpinUp, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT, int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnHyperCubicLatticeWithSpinMomentumSpace(const FermionOnHyperCubicLatticeWithSpinMomentumSpace& fermions);

  // destructor
  //
  ~FermionOnHyperCubicLatticeWithSpinMomentumSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnHyperCubicLatticeWithSpinMomentumSpace& operator = (const FermionOnHyperCubicLatticeWithSpinMomentumSpace& fermions);

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
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentKt = current momentum along t for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // currentTotalKt = current total momentum along t
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt);

  // evaluate Hilbert space dimension with a fixed number of fermions with spin up
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentKt = current momentum along t for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // currentTotalKt = current total momentum along t
  // nbrSpinUp = number of fermions with spin up
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, int nbrSpinUp);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentKt = current momentum along t for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // currentTotalKt = current total momentum along t
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentKt = current momentum along t for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // currentTotalKt = current total momentum along t
  // nbrSpinUp = number of fermions with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, int nbrSpinUp, long pos);

};


#endif


