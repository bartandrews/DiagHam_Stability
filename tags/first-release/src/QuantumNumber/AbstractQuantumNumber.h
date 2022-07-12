////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of abstract quantum number                     //
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


#ifndef ABSTRACTQUANTUMNUMBER_H
#define ABSTRACTQUANTUMNUMBER_H


#include "config.h"

#include <iostream>


using std::ostream;


class AbstractQuantumNumber
{

 protected:

  int QuantumNumberType;

 public:

  enum QuantumNumberID
  {
    Void = 0x00000000,
    Sz = 0x00000001,
    NumberParticle = 0x00000002,
    Vector = 0x00000004,
    Momentum = 0x00000008,
    PeriodicMomentum = 0x00000010
  };

  // default constructor
  //
  AbstractQuantumNumber ();

  // copy constructor
  //
  // Q = quantum number to copy
  AbstractQuantumNumber (const AbstractQuantumNumber& Q);

  // virtual destructor
  //
  virtual ~AbstractQuantumNumber ();

  // assignement (without duplicating datas)
  //
  // Q = quantum number to copy
  // return value = reference on current quantum number
  virtual AbstractQuantumNumber& operator = (const AbstractQuantumNumber& Q);

  // clone current quantum number
  //
  // return value = pointer on cloned quantum number  
  virtual AbstractQuantumNumber* Clone ();

  // get quantum number type
  //
  // return value = Hilbert space dimension
  int GetQuantumNumberType();

  // add two quantum numbers
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  // return value = pointer to the sum of the two quantum numbers
  friend AbstractQuantumNumber* operator + (const AbstractQuantumNumber& Q1, 
					    const AbstractQuantumNumber& Q2);

  // substract two quantum numbers
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  // return value = pointer to the difference of the two quantum numbers
  friend AbstractQuantumNumber* operator - (const AbstractQuantumNumber& Q1, 
					    const AbstractQuantumNumber& Q2);

  // test if two quantum numbers are identical
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  friend bool operator == (const AbstractQuantumNumber& Q1, const AbstractQuantumNumber& Q2);

  // test if two quantum numbers are different
  //
  // Q1 = first quantum number
  // Q2 = second quantum number
  friend bool operator != (const AbstractQuantumNumber& Q1, const AbstractQuantumNumber& Q2);

  // print quantum number
  //
  // Str = reference on current output stream 
  // Q = quantum number to print
  // return value = reference on current output stream 
  friend ostream& operator << (ostream& Str, const AbstractQuantumNumber& Q1);

};

#endif


