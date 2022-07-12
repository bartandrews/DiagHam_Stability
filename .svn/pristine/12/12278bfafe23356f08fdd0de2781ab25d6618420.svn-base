////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on sphere with spin and pairing             //
//      (i.e. Sz conservation but no conservation of the particle number)     //
//                  without any contraint on the momentum along z             //
//                                                                            //
//                        last modification : 19/09/2016                      //
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


#ifndef FERMIONONSPHEREWITHSPINANDPAIRINGALLLZ_H
#define FERMIONONSPHEREWITHSPINANDPAIRINGALLLZ_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpinAndPairing.h"

#include <iostream>


class FermionOnSphereWithSpinAndPairingAllLz :  public FermionOnSphereWithSpinAndPairing
{


 protected:


 public:

  // default constructor
  //
  FermionOnSphereWithSpinAndPairingAllLz();

  // basic constructor
  // 
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twice the total spin value
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinAndPairingAllLz (int lzMax, int totalSpin, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinAndPairingAllLz(const FermionOnSphereWithSpinAndPairingAllLz& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinAndPairingAllLz ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinAndPairingAllLz& operator = (const FermionOnSphereWithSpinAndPairingAllLz& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();



 protected:

  // evaluate Hilbert space dimension
  //
  // lzMax = momentum maximum value for a fermion
  // totalSpin = twice the total spin
  // return value = Hilbert space dimension      
  virtual long EvaluateHilbertSpaceDimension(int lzMax, int totalSpin);


  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion in the state
  // totalSpin = twice the total spin
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int lzMax, int totalSpin, long pos);
  
};

#endif


