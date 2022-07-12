////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                  class of fermions on the 4D manifold S2 x S2              //
//                                                                            //
//                        last modification : 26/09/2016                      //
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


#ifndef FERMIONONS2XS2_H
#define FERMIONONS2XS2_H

#include "config.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

#include <iostream>



class FermionOnS2xS2 : public FermionOnSquareLatticeMomentumSpace
{

 protected:

  // number of flux quanta for the first sphere
  int NbrFluxQuanta1;
  // number of flux quanta for the second sphere
  int NbrFluxQuanta2;

  // total angular momentum along z for the first sphere
  int TotalLz1;
  // total angular momentum along z for the second sphere
  int TotalLz2;

 public:

  // default constructor
  // 
  FermionOnS2xS2 ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrFluxQuanta1 = number of flux quanta for the first sphere
  // nbrFluxQuanta2 = number of flux quanta for the second sphere
  // totalLz1 = total angular momentum for the first sphere
  // totalLz2 = total angular momentum for the second sphere
  // memory = amount of memory granted for precalculations
  FermionOnS2xS2 (int nbrFermions, int nbrFluxQuanta1, int nbrFluxQuanta2, int totalLz1, int totalLz2, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnS2xS2(const FermionOnS2xS2& fermions);

  // destructor
  //
  ~FermionOnS2xS2 ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnS2xS2& operator = (const FermionOnS2xS2& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos);

};


#endif


