////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of bosons on disk  using the Haldane basis            //
//                            for system size such that                       //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 08/07/2008                      //
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


#ifndef BOSONONDISKHALDANEBASISSHORT_H
#define BOSONONDISKHALDANEBASISSHORT_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"


class BosonOnDiskHaldaneBasisShort :  public BosonOnSphereHaldaneBasisShort
{


 public:

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum angular momentum that a single particle can reach
  // referenceState = array that describes the reference state to start from
  BosonOnDiskHaldaneBasisShort (int nbrBosons, int totalLz, int lzMax, int* referenceState);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnDiskHaldaneBasisShort(const BosonOnDiskHaldaneBasisShort& bosons);

  // destructor
  //
  ~BosonOnDiskHaldaneBasisShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnDiskHaldaneBasisShort& operator = (const BosonOnDiskHaldaneBasisShort& bosons);

};

#endif


