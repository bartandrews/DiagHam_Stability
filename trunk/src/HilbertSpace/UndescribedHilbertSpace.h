////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of Hilbert space with no state description             //
//                                                                            //
//                        last modification : 07/06/2001                      //
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


#ifndef UNDESCRIBEDHILBERTSPACE_H
#define UNDESCRIBEDHILBERTSPACE_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"


class UndescribedHilbertSpace : public AbstractHilbertSpace
{

 protected:

  AbstractQuantumNumber* QuantumNumber;

 public:

  // default constructor
  //
  UndescribedHilbertSpace ();

  // constructor from datas
  //
  // dimension = Hilbert space dimension
  UndescribedHilbertSpace (int dimension);

  // constructor from datas
  //
  // dimension = Hilbert space dimension
  // quantumNumber = reference on the quantum number associated to the Hilbert space
  UndescribedHilbertSpace (int dimension, AbstractQuantumNumber& quantumNumber);

  // copy constructor (without duplicating datas)
  //
  // space = reference on Hilbert space to copy
  UndescribedHilbertSpace (const UndescribedHilbertSpace& space);

  // destructor
  //
  ~UndescribedHilbertSpace ();

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // assignment (without duplicating datas)
  //
  // space = reference on Hilbert space to copy
  // return value = reference on current Hilbert space
  UndescribedHilbertSpace& operator = (const UndescribedHilbertSpace& space);

  // return Hilbert space dimension
  //
  // return value = Hilbert space dimension
  int GetHilbertSpaceDimension();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
					 SubspaceSpaceConverter& converter);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

};

#endif


