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
 

#ifndef DOUBLEDSPIN0_1_2_CHAINWITHTRANSLATIONSANDZZSYMMETRYANDSUBLATTICEQUANTUMNUMBERS_H
#define DOUBLEDSPIN0_1_2_CHAINWITHTRANSLATIONSANDZZSYMMETRYANDSUBLATTICEQUANTUMNUMBERS_H 


#include "config.h"
#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetry.h"


#include <iostream>


using std::ostream;

class DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers : public DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetry 
{
  friend class Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers;
 protected:
  int SubLatticeDifferenceKet;
  int SubLatticeDifferenceBra;
  int SubLatticeDifferenceProduct;

 public:

  // default constructor
  //
  DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers ();
  
  // constructor for complete Hilbert space with no restriction on momentum
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // memorySize = memory size in bytes allowed for look-up table
  // memorySlice = maximum amount of memory that can be allocated to partially evalauted the states
  DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers (int chainLength,  int diffSz, int subLatticeDifferenceKet, int subLatticeDifferenceBra, int subLatticeDifferenceProduct, int memorySize, int memorySlice);
   
  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers (int chainLength, int momentum,  int translationStep, int sz, int subLatticeDifferenceKet, int subLatticeDifferenceBra, int subLatticeDifferenceProduct, int memorySize, int memorySlice);
  
  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers (const DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers & chain);
  
  // destructor
  //
  ~DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers & operator = (const DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers & chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  HermitianMatrix  EvaluatePartialDensityMatrix (int szSector, int momentumSector, int sublatticeQuantumNumberSector,  int complementarySublatticeQuantumNumberSector, ComplexVector& groundState);
 protected:

};


#endif


