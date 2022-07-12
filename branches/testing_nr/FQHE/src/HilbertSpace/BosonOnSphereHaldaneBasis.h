////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of bosons on sphere using the Haldane basis             //
//                                                                            //
//                        last modification : 01/10/2007                      //
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


#ifndef BOSONONSPHEREHALDANEBASIS_H
#define BOSONONSPHEREHALDANEBASIS_H


#include "config.h"
#include "HilbertSpace/BosonOnSphere.h"

#include <iostream>


class BosonOnSphereHaldaneBasis :  public BosonOnSphere
{

 protected:

  // topmost state 
  int* ReferenceState;

  // three temporary arrays used during Hilbert space generation
  int** TmpGeneratedStates;
  int* TmpGeneratedStatesLzMax;
  unsigned long* KeepStateFlag;

 public:

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = reference on twice the momentum total value, totalLz will be recomputed from referenceState and stored in totalLz
  // lzMax = maximum Lz value reached by a boson
  // referenceState = array that describes the reference state to start from
  BosonOnSphereHaldaneBasis (int nbrBosons, int& totalLz, int lzMax, int* referenceState);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  BosonOnSphereHaldaneBasis (char* fileName);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereHaldaneBasis(const BosonOnSphereHaldaneBasis& bosons);

  // destructor
  //
  virtual ~BosonOnSphereHaldaneBasis ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereHaldaneBasis& operator = (const BosonOnSphereHaldaneBasis& bosons);

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // convert a gien state from Haldane basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToNbodyBasis(RealVector& state, BosonOnSphere& nbodyBasis);

 protected:

  // find state index
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  int FindStateIndex(int* stateDescription, int lzmax);

  // generate all states corresponding to the constraints
  // 
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual int HaldaneGenerateStates(int lzMax, int* referenceState, int pos, long& memory);

};

#endif


