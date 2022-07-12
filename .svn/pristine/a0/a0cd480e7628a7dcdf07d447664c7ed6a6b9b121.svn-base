////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with translastion invariance         //  
//                             and fixed parity                               //
//                                                                            //
//                        last modification : 06/06/2014                      //
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


#ifndef SPIN1_2CHAINFIXEDPARITYWITHTRANSLATIONS_H
#define SPIN1_2CHAINFIXEDPARITYWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "Matrix/HermitianMatrix.h"

#include <iostream>


using std::ostream;


class Spin1_2ChainFixedParityWithTranslations : public Spin1_2ChainWithTranslations
{

 protected:

  // parity of the total (Sz + 1/2) (can be 0 or 1)
  int SzParity;

 public:

  // default constructor
  //
  Spin1_2ChainFixedParityWithTranslations ();

  // constructor for complete Hilbert space corresponding to a given parity
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // translationStep = indicates the step for an elementary translation
  // memorySize = memory size in bytes allowed for look-up table
  // memorySlice = maximum amount of memory that can be allocated to partially evalauted the states
  Spin1_2ChainFixedParityWithTranslations (int chainLength, int momentum, int translationStep, int parity, 
					   int memorySize = 10000000, int memorySlice = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainFixedParityWithTranslations (const Spin1_2ChainFixedParityWithTranslations& chain);

  // destructor
  //
  ~Spin1_2ChainFixedParityWithTranslations ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainFixedParityWithTranslations& operator = (const Spin1_2ChainFixedParityWithTranslations& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // compute the parity (prod_i Sz_i) for a given state
  //
  // state = index of the state to be applied on Sz_i operator
  // return value = 0 if prod_i Sz_i = 1, 1 if prod_i Sz_i = -1
  virtual unsigned long Parity (int state);

 protected:

  // constructor from pre-constructed datas
  //
  // hilbertSpaceDimension = Hilbert space dimension
  // chainDescription = array describing states
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // parity = parity of the total (Sz + 1/2) (can be 0 or 1)
  // lookUpTableShift = shift to apply to a state to obtain an index to the look-up table 
  // complementaryStateShift = shift to apply to move the spin from one end to the other one
  Spin1_2ChainFixedParityWithTranslations (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
					   int momentum, int parity, int lookUpTableShift, 
					   int complementaryStateShift);

};

// compute the parity (prod_i Sz_i) for a given state
//
// state = index of the state to be applied on Sz_i operator
// return value = 0 if prod_i Sz_i = 1, 1 if prod_i Sz_i = -1

inline unsigned long Spin1_2ChainFixedParityWithTranslations::Parity (int state)
{
  return (unsigned long) this->SzParity;
}


#endif


