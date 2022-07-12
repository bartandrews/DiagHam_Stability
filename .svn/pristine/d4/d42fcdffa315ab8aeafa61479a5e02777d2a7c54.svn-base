////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin                    //
//                  including the Haldane squeezing technique                 //
//        allow LzMax up to 63 (for systems with 128 bit integer support)     //
//             or 31 (on 32 bit systems without 128 bit integer support)      //
//                                                                            //
//                        last modification : 18/03/2009                      //
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


#ifndef FERMIONONSPHEREWITHSPINHALDANEBASISLONG_H
#define FERMIONONSPHEREWITHSPINHALDANEBASISLONG_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"

#include <iostream>


class FermionOnSphereWithSpinHaldaneBasisLong :  public FermionOnSphereWithSpinLong
{

 protected:

  // array of root partitions describing the squeezed basis
  ULONGLONG* RootPartitions;
  // number of root partitions
  int NbrRootPartitions;

  // three temporary arrays used during Hilbert space generation
  ULONGLONG* TmpGeneratedStates;
  int* TmpGeneratedStatesLzMax;
  unsigned long* KeepStateFlag;

 public:

  // default constructor
  //
  FermionOnSphereWithSpinHaldaneBasisLong();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twce the total spin value
  // rootPartitions = array of root partitions describing the squeezed basis
  // nbrRootPartitions = number of root partitions
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinHaldaneBasisLong (int nbrFermions, int& totalLz, int lzMax, int& totalSpin, 
					   int** rootPartitions, int nbrRootPartitions, 
					   unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinHaldaneBasisLong(const FermionOnSphereWithSpinHaldaneBasisLong& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinHaldaneBasisLong ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinHaldaneBasisLong& operator = (const FermionOnSphereWithSpinHaldaneBasisLong& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // convert a gien state from Haldane basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  virtual RealVector ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpinLong& nbodyBasis);

  // convert a given state from the usual n-body basis to the Haldane basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  virtual RealVector ConvertFromNbodyBasis(RealVector& state, FermionOnSphereWithSpinLong& nbodyBasis);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(ULONGLONG stateDescription, int lzmax);


  // generate all squeezed states from a root partition
  // 
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateSqueezedStates(int lzMax, ULONGLONG referenceState, long pos, long& memory);

};


#endif


