////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of Sz quantum number                        //
//                                                                            //
//                        last modification : 19/04/2001                      //
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


#ifndef SZQUANTUMNUMBER_H
#define SZQUANTUMNUMBER_H


#include "config.h"
#include "QuantumNumber/AbstractQuantumNumber.h"


#include <iostream>


using std::ostream;


class SzQuantumNumber : public AbstractQuantumNumber
{

 protected:

  int Sz;

 public:

  // default constructor
  //
  SzQuantumNumber ();

  // constructor from a Sz value
  //
  // sz = twice the value of  Sz component
  SzQuantumNumber (int sz);

  // copy constructor
  //
  // Q = quantum number to copy
  SzQuantumNumber (const SzQuantumNumber& Q);

  // destructor
  //
  ~SzQuantumNumber ();

  // assignement
  //
  // Q = quantum number to copy
  // return value = reference on current quantum number
  SzQuantumNumber& operator = (const SzQuantumNumber& Q);

  // set total Sz component
  //
  // sz = value to assign
  // return value = reference on current quantum number
  SzQuantumNumber& operator = (int sz);

  // clone current quantum number
  //
  // return value = pointer on cloned quantum number  
  AbstractQuantumNumber* Clone ();

  // Get total Sz component
  //
  // return value = total Sz component
  int GetSz ();

  // add two quantum numbers
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  // return value = sum of the two quantum numbers
  friend SzQuantumNumber operator + (const SzQuantumNumber& Q1, const SzQuantumNumber& Q2);

  // substract two quantum numbers
  //
  // Q1 = first quantum number
  // Q2 = quantum number to substract
  // return value = difference of the two quantum numbers
  friend SzQuantumNumber operator - (const SzQuantumNumber& Q1, const SzQuantumNumber& Q2);

  // test if two quantum numbers are identical
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  // return value = true if quantum numbers are identical
  friend bool operator == (const SzQuantumNumber& Q1, const SzQuantumNumber& Q2);

  // test if two quantum numbers are different
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  // return value = true if quantum numbers are different
  friend bool operator != (const SzQuantumNumber& Q1, const SzQuantumNumber& Q2);

  // print quantum number
  //
  // Str = reference on current output stream 
  // Q = quantum number to print
  // return value = reference on current output stream 
  friend ostream& operator << (ostream& Str, const SzQuantumNumber& Q);

};

#endif


