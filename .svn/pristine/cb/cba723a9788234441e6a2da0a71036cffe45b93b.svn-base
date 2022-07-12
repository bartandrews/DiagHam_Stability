////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of bosons on sphere  using the Haldane basis           //
//                            for system size such that                       //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                           with Lz <-> -Lz symmetry                         //
//                                                                            //
//                        last modification : 22/11/2010                      //
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


#ifndef BOSONONSPHEREHALDANESYMMETRICBASISSHORT_H
#define BOSONONSPHEREHALDANESYMMETRICBASISSHORT_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Vector/RationalVector.h"
#include "Vector/LongRationalVector.h"

#include <iostream>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class RationalPolynomial;
class LongRationalPolynomial;


class BosonOnSphereHaldaneSymmetricBasisShort :  public BosonOnSphereSymmetricBasisShort
{

 public:

  // default constructor
  //
  BosonOnSphereHaldaneSymmetricBasisShort ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // lzMax = maximum Lz value reached by a boson
  // referenceState = array that describes the reference state to start from
  BosonOnSphereHaldaneSymmetricBasisShort (int nbrBosons, int lzMax, int* referenceState);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  BosonOnSphereHaldaneSymmetricBasisShort (char* fileName);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereHaldaneSymmetricBasisShort(const BosonOnSphereHaldaneSymmetricBasisShort& bosons);

  // destructor
  //
  virtual ~BosonOnSphereHaldaneSymmetricBasisShort ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereHaldaneSymmetricBasisShort& operator = (const BosonOnSphereHaldaneSymmetricBasisShort& bosons);

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

};

#endif
