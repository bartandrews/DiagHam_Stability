////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of bosons on disk                         //
//                                                                            //
//                        last modification : 04/06/2002                      //
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


#ifndef BOSONONDISK_H
#define BOSONONDISK_H


#include "config.h"
#include "HilbertSpace/BosonOnSphere.h"


class BosonOnDisk :  public BosonOnSphere
{


 public:

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = momentum total value
  // lzMax = maximum angular momentum that a single particle can reach (negative if it has to be deduced from nbrBosons and totalLz)
  BosonOnDisk (int nbrBosons, int totalLz, int lzMax = -1);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnDisk(const BosonOnDisk& bosons);

  // destructor
  //
  ~BosonOnDisk ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnDisk& operator = (const BosonOnDisk& bosons);

  // forge an eigenstate from a description given by a file
  //
  // filename = name of the file that contains the state description
  // state = reference on the vector where the state has to be stored
  // return value = true if no error occured
  bool ForgeEigenstate(char* filename, RealVector& state);

};

#endif


