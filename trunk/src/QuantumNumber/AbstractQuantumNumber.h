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

 public:

  int QuantumNumberType;

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
  virtual int GetQuantumNumberType();

  // add a quantum nunber to the current one
  //
  // number = quantum number to add 
  // return value = reference to the current quantum number
  virtual AbstractQuantumNumber& operator += (const AbstractQuantumNumber& number);

  // substract a quantum nunber to the current one
  //
  // number = quantum number to add 
  // return value = reference to the current quantum number
  virtual AbstractQuantumNumber& operator -= (const AbstractQuantumNumber& number);

  // add a quantum nunber to the current one and return a new quantum number (not modifying the current one)
  //
  // number = quantum number to add to the current one
  // return value = pointer on cloned quantum number corresponding to the sum 
  virtual AbstractQuantumNumber* Add (const AbstractQuantumNumber& number);

  // substract a quantum nunber to the current one and return a new quantum number (not modifying the current one)
  //
  // number = quantum number to substract to the current one
  // return value = pointer on cloned quantum number corresponding to the difference 
  virtual AbstractQuantumNumber* Sub (const AbstractQuantumNumber& number);

  // test if two quantum numbers are identical
  //
  // number = quantum number to compare to the current one
  virtual bool IsEqual (const AbstractQuantumNumber& number);

  // test if two quantum numbers are different
  //
  // number = quantum number to compare to the current one
  virtual bool IsDifferent (const AbstractQuantumNumber& number);

  // print quantum number
  //
  // str = reference on current output stream 
  // return value = reference on current output stream 
  virtual ostream& PrintQuantumNumber (ostream& str);

};

#endif


