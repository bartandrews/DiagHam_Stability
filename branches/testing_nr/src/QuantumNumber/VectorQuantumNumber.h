////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of vector of quantum numbers                    //
//                                                                            //
//                        last modification : 11/05/2001                      //
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


#ifndef VECTORQUANTUMNUMBER_H
#define VECTORQUANTUMNUMBER_H


#include "config.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "GeneralTools/List.h"

#include <iostream>


using std::ostream;


class VectorQuantumNumber : public AbstractQuantumNumber
{

 protected:

  List<AbstractQuantumNumber*> QuantumNumbers;

 public:

  // default constructor
  //
  VectorQuantumNumber ();

  // constructor from a list of quantum numbers
  //
  // quantumNumbers = list of quantum numbers
  VectorQuantumNumber (List<AbstractQuantumNumber*> quantumNumbers);

  // copy constructor
  //
  // Q = quantum number to copy
  VectorQuantumNumber (const VectorQuantumNumber& Q);

  // destructor
  //
  ~VectorQuantumNumber ();

  // assignement
  //
  // Q = quantum number to copy
  // return value = reference on current quantum number
  VectorQuantumNumber& operator = (const VectorQuantumNumber& Q);

  // clone current quantum number
  //
  // return value = pointer on cloned quantum number  
  AbstractQuantumNumber* Clone ();

  // add a quantum number
  //
  // numberParticle = value to assign
  // return value = reference on current quantum number
  VectorQuantumNumber& operator += (AbstractQuantumNumber* quantumNumber);

  // Set i-th component of vector of quantum numbers
  //
  // index = index of the component to set 
  // quantumNumber = new quantum number
  // return value = reference on current quantum number
  VectorQuantumNumber& SetQuantumNumber(int index, AbstractQuantumNumber* quantumNumber);

  // Get i-th component of vector of quantum numbers
  //
  // index = index of the component to get 
  // return value = pointer to the quantum number
  AbstractQuantumNumber* operator [] (int index);

  // Get reference on the list of quantum numbers that are part of the vector
  //
  // return value = reference on the list of quantum numbers
  List<AbstractQuantumNumber*>& GetQuantumNumbers();

  // add a quantum nunber to the current one
  //
  // number = quantum number to add 
  // return value = reference to the current quantum number
  AbstractQuantumNumber& operator += (const AbstractQuantumNumber& number);

  // substract a quantum nunber to the current one
  //
  // number = quantum number to add 
  // return value = reference to the current quantum number
  AbstractQuantumNumber& operator -= (const AbstractQuantumNumber& number);

  // test if two quantum numbers are identical
  //
  // number = quantum number to compare to the current one
  bool IsEqual (const AbstractQuantumNumber& number);

  // test if two quantum numbers are different
  //
  // number = quantum number to compare to the current one
  bool IsDifferent (const AbstractQuantumNumber& number);

  // print quantum number
  //
  // str = reference on current output stream 
  // return value = reference on current output stream 
  ostream& PrintQuantumNumber (ostream& str);

};

// Get reference on the list of quantum numbers that are part of the vector
//
// return value = reference on the list of quantum numbers

inline List<AbstractQuantumNumber*>& VectorQuantumNumber::GetQuantumNumbers()
{
  return this->QuantumNumbers;
}

#endif


