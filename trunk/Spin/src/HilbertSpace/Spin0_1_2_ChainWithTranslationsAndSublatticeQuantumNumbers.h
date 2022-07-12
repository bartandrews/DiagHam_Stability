////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of doubled spin 0 +1/2 chain with translations              //
//                                                                            //
//                        last modification : 21/01/2016                      //
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


#ifndef SPIN0_1_2_CHAINWITHTRANSLATIONSANDSUBLATTICEQUANTUMNUMBERS_H
#define SPIN0_1_2_CHAINWITHTRANSLATIONSANDSUBLATTICEQUANTUMNUMBERS_H


#include "config.h"
#include "HilbertSpace/Spin0_1_2_ChainWithTranslations.h"

#include <iostream>


using std::ostream;

class Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers : public  Spin0_1_2_ChainWithTranslations
{
  friend class DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers;
  friend class DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers;
 protected:
  int SubLatticeDifference;

 public:

  // default constructor
  //
  Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers ();

  // constructor for complete Hilbert space with no restriction on momentum
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // memorySize = memory size in bytes allowed for look-up table
  // memorySlice = maximum amount of memory that can be allocated to partially evalauted the states
  Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers (int chainLength,  int diffSz, int subLatticeDifference, int memorySize, int memorySlice);
  Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers (int chainLength,  int subLatticeDifference, int memorySize=100000, int memorySlice=100000);
  
  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers (int chainLength, int momentum,  int translationStep, int sz,  int subLatticeDifference, int memorySize, int memorySlice);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers (const Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers & chain);

  // destructor
  //
  ~Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers & operator = (const Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers & chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:
  inline void ComputeDifferenceSubLatticeNumberZero(unsigned long state,  int & sublatticeNumber);

};

inline void Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers::ComputeDifferenceSubLatticeNumberZero(unsigned long state,  int & sublatticeNumber)
{
  unsigned long SourceState =  state;
  sublatticeNumber = 0;
  unsigned  int Tmp;
  for (int p = 0;p < this->ChainLength;p++)
    {
      Tmp = SourceState & 0x3ul;
      if( p%2 == 0 )
	{
	  if (Tmp == 2)
	      sublatticeNumber++;

	}
      else
	{
	  if (Tmp == 2)
	    sublatticeNumber--;
	}
      SourceState >>=2;
    }
}
#endif


