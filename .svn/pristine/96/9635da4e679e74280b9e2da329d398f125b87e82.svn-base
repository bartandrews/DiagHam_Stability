////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2014 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons on the 4D manifold T2 x T2               //
//                                                                            //
//                        last modification : 15/02/2017                      //
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


#ifndef BOSONONT2XT2_H
#define BOSONONT2XT2_H


#include "config.h"
#include "HilbertSpace/BosonOnS2xS2.h"

#include <iostream>



class BosonOnT2xT2 : public BosonOnS2xS2
{
  
 protected:

 public:

  // default constructor
  // 
  BosonOnT2xT2 ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta1 = number of flux quanta for the first torus
  // nbrFluxQuanta2 = number of flux quanta for the second torus
  // totalKy1 = total momentum along the y of the first torus
  // totalKy2 = total momentum along the y of the second torus
  // memory = amount of memory granted for precalculations
  BosonOnT2xT2 (int nbrBosons, int nbrFluxQuanta1, int nbrFluxQuanta2, int totalKy1, int totalKy2, unsigned long memory = 10000000);

  // copy constructor (without duplicating data)
  //
  // space = reference on the hilbert space to copy
  BosonOnT2xT2(const BosonOnT2xT2& space);

  // destructor
  //
  ~BosonOnT2xT2 ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnT2xT2 & operator = (const BosonOnT2xT2 & bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

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
  // stateDescription = array that gives each state description
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentFermionicPosition = current fermionic position within the state description
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(unsigned long* stateDescription, int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentFermionicPosition, long pos);

};

#endif

