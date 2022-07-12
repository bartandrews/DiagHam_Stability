////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of abstract Hilbert space                      //
//                                                                            //
//                        last modification : 23/04/2001                      //
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


#ifndef ABSTRACTHILBERTSPACE_H
#define ABSTRACTHILBERTSPACE_H


#include "config.h"
#include "GeneralTools/List.h"
#include "GeneralTools/GarbageFlag.h"


#include <iostream>


using std::ostream;


class AbstractQuantumNumber;
class SubspaceSpaceConverter;


class AbstractHilbertSpace
{

 protected:

   // dimension of the hilbert space
  int HilbertSpaceDimension;

  // dimension of the hilbert space for large Hilbert space
  long LargeHilbertSpaceDimension;

  // garbage flag used to share datas from different copy of the same hilbert space
  GarbageFlag Flag;

 public:

  // virtual destructor
  //
  virtual ~AbstractHilbertSpace ();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone() = 0;

  // return Hilbert space dimension
  //
  // return value = Hilbert space dimension
  virtual int GetHilbertSpaceDimension();

  // return Hilbert space dimension for large Hilbert space
  //
  // return value = Hilbert space dimension
  virtual long GetLargeHilbertSpaceDimension();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers () = 0;

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index) = 0;

  // get information about any additional symmetry of the Hilbert space
  //
  // return value = symmetry id  
  int GetHilbertSpaceAdditionalSymmetry();

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter) = 0;

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state) = 0;

};

// return Hilbert space dimension
//
// return value = Hilbert space dimension

inline int AbstractHilbertSpace::GetHilbertSpaceDimension() 
{
  return this->HilbertSpaceDimension;
}

// return Hilbert space dimension for large Hilbert space
//
// return value = Hilbert space dimension

inline long AbstractHilbertSpace::GetLargeHilbertSpaceDimension()
{
  return this->LargeHilbertSpaceDimension;
}

#endif


