////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of bosons on disk with SU(2) spin                 //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 11/10/2012                      //
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


#ifndef BOSONONDISKWITHSU2SPIN_H
#define BOSONONDISKWITHSU2SPIN_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"


class BosonOnDiskWithSU2Spin :  public BosonOnSphereWithSU2Spin
{

 public:
  
  // basic constructor without any constraint on Sz
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum angular momentum that a single particle can reach (negative if it has to be deduced from nbrBosons and totalLz)
  BosonOnDiskWithSU2Spin (int nbrBosons, int totalLz, int lzMax = -1);

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // totalSpin = twice the total spin value
  // lzMax = twice the maximum Lz value reached by a boson
  BosonOnDiskWithSU2Spin (int nbrBosons, int totalLz, int totalSpin, int lzMax);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnDiskWithSU2Spin(const BosonOnDiskWithSU2Spin& bosons);

  // destructor
  //
  ~BosonOnDiskWithSU2Spin ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnDiskWithSU2Spin& operator = (const BosonOnDiskWithSU2Spin& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();	

};

#endif


