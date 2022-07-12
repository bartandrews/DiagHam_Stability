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

  // add two quantum numbers
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  // return value = sum of the two quantum numbers
  friend VectorQuantumNumber operator + (VectorQuantumNumber& Q1, VectorQuantumNumber& Q2);

  // substract two quantum numbers
  //
  // Q1 = first quantum number
  // Q2 = quantum number to substract
  // return value = sum of the two quantum numbers
  friend VectorQuantumNumber operator - (VectorQuantumNumber& Q1, VectorQuantumNumber& Q2);

  // test if two quantum numbers are identical
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  // return value = true if quantum numbers are identical
  friend bool operator == (const VectorQuantumNumber& Q1, const VectorQuantumNumber& Q2);

  // test if two quantum numbers are different
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  // return value = true if quantum numbers are different
  friend bool operator != (const VectorQuantumNumber& Q1, const VectorQuantumNumber& Q2);

  // print quantum number
  //
  // Str = reference on current output stream 
  // Q = quantum number to print
  // return value = reference on current output stream 
  friend ostream& operator << (ostream& Str, const VectorQuantumNumber& Q);

};

#endif


