////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                class of fermions on cubic lattice with SU(4) spin          //
//               in momentum space that allows a number of sites up to        //
//                  31 (for systems with 128 bit integer support)             //
//               or 15 (on 32 bit systems without 128 bit integer support)    //                        
//                                                                            //
//                        last modification : 24/09/2011                      //
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


#ifndef FERMIONONCUBICLATTICEWITHSU4SPINMOMENTUMSPACELONG_H
#define FERMIONONCUBICLATTICEWITHSU4SPINMOMENTUMSPACELONG_H

#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSU4SpinLong.h"

#include <iostream>



class FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong : public FermionOnSphereWithSU4SpinLong
{

 protected:

  // number of sites in the x direction
  int NbrSiteX;
  // number of sites in the y direction
  int NbrSiteY;
  // number of sites in the z direction
  int NbrSiteZ;
  // number of sites in the direction perpendicular to X
  int NbrSiteYZ;

  // momentum along the x direction
  int KxMomentum;
  // momentum along the y direction
  int KyMomentum;
  // momentum along the z direction
  int KzMomentum;

  // flag to indicate that the Hilbert space should preserve Sz
  bool SzFlag;
  // flag to indicate that the Hilbert space should preserve Pz (Pz being the difference N_A - N_B)
  bool PzFlag;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // kzMomentum = momentum along the z direction
  // memory = amount of memory granted for precalculations
  FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int kxMomentum, int kyMomentum, int kzMomentum, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong(const FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong& fermions);

  // destructor
  //
  ~FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong& operator = (const FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong& fermions);

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
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentKz = current momentum along z for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentTotalKz = current total momentum along z
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz, long pos);


};


#endif


