////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of bosons on cylinder for system size such that            //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 11/12/2009                      //
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


#ifndef BOSONONCYLINDERSHORT_H
#define BOSONONCYLINDERSHORT_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereShort.h"


class BosonOnCylinderShort :  public BosonOnSphereShort
{


 public:

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalKy = momentum total value
  // kyMax = maximum momentum that a single particle can reach
  BosonOnCylinderShort (int nbrBosons, int totalKy, int kyMax);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnCylinderShort(const BosonOnCylinderShort& bosons);

  // destructor
  //
  ~BosonOnCylinderShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnCylinderShort& operator = (const BosonOnCylinderShort& bosons);

};

#endif


