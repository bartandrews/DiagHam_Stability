////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere including two                //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 19/05/2009                      //
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


#ifndef FERMIONONSPHERETWOLANDAULEVELS_H
#define FERMIONONSPHERETWOLANDAULEVELS_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

#include <iostream>



class FermionOnSphereTwoLandauLevels :  public FermionOnSphereWithSpin
{

 protected:

  // maximum Lz value reached by a fermion with a spin up
  int LzMaxUp;
  // maximum Lz value reached by a fermion with a spin down
  int LzMaxDown;
  // shift to apply on the spin up part
  int LzShiftUp;
  // shift to apply on the spin down part
  int LzShiftDown;
  // sum of LzShiftUp and LzShiftDown
  int LzTotalShift;

 public:

  // default constructor
  //
  FermionOnSphereTwoLandauLevels();

  // basic constructor with no contraint on the number of particles per spin component
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMaxUp = twice the maximum Lz value reached by a fermion with a spin up
  // lzMaxDown = twice the maximum Lz value reached by a fermion with a spin down
  // memory = amount of memory granted for precalculations
  FermionOnSphereTwoLandauLevels (int nbrFermions, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereTwoLandauLevels(const FermionOnSphereTwoLandauLevels& fermions);

  // destructor
  //
  ~FermionOnSphereTwoLandauLevels ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereTwoLandauLevels& operator = (const FermionOnSphereTwoLandauLevels& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  virtual RealVector ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermionsUp = number of fermions with spin up
  // nbrFermionsDown = number of fermions with spin down
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension      
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz);

  // evaluate Hilbert space dimension without constraint on the number of particles per level
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateFullHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // generate all states corresponding to the constraints
  // 
  // nbrFermionsUp = number of fermions with spin up
  // nbrFermionsDown = number of fermions with spin down
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMaxUp = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateFullStates(int nbrFermions, int lzMax, int totalLz, long pos);

};

#endif


