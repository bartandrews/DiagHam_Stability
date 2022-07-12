////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin                    //
//    including the Haldane squeezing technique for large system size         //
//                                                                            //
//                        last modification : 05/09/2009                      //
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


#ifndef FERMIONONSPHEREWITHSPINHALDANELARGEBASIS_H
#define FERMIONONSPHEREWITHSPINHALDANELARGEBASIS_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasis.h"

#include <iostream>


class FermionOnSphereWithSpinHaldaneLargeBasis :  public FermionOnSphereWithSpinHaldaneBasis
{

 protected:

  // look-up table for large Hilbert space
  long** LargeLookUpTable;

 public:

  // default constructor
  //
  FermionOnSphereWithSpinHaldaneLargeBasis();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twce the total spin value
  // rootPartitions = array of root partitions describing the squeezed basis
  // nbrRootPartitions = number of root partitions
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinHaldaneLargeBasis (int nbrFermions, int& totalLz, int lzMax, int& totalSpin, 
					    int** rootPartitions, int nbrRootPartitions, 
					    unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinHaldaneLargeBasis(const FermionOnSphereWithSpinHaldaneLargeBasis& fermions);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinHaldaneLargeBasis (char* fileName, unsigned long memory = 10000000);

  // destructor
  //
  ~FermionOnSphereWithSpinHaldaneLargeBasis ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinHaldaneLargeBasis& operator = (const FermionOnSphereWithSpinHaldaneLargeBasis& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  bool WriteHilbertSpace (char* fileName);

  // create the Jack polynomial decomposition corresponding to the root partition
  //
  // jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
  // alpha = value of the Jack polynomial alpha coefficient
  // return value = decomposition of the corresponding Jack polynomial on the unnormalized basis
  virtual RealVector& GenerateJackPolynomial(RealVector& jack, double alpha);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // return value = corresponding index
  virtual long FindStateIndex(unsigned long stateDescription);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(unsigned long memory);

  // generate all squeezed states from a root partition
  // 
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateSqueezedStates(int lzMax, unsigned long referenceState, long pos, long& memory);

};


#endif


