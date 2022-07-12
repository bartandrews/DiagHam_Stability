////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of momenetum quantum number                     //
//                                                                            //
//                        last modification : 29/01/2002                      //
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


#ifndef MOMENTUMQUANTUMNUMBER_H
#define MOMENTUMQUANTUMNUMBER_H


#include "config.h"
#include "QuantumNumber/AbstractQuantumNumber.h"

#include <iostream>


using std::ostream;


class MomentumQuantumNumber : public AbstractQuantumNumber
{

 protected:

  int Momentum;

 public:

  // default constructor
  //
  MomentumQuantumNumber ();

  // constructor from a momentum value
  //
  // momentum = momentum value
  MomentumQuantumNumber (int momentum);

  // copy constructor
  //
  // Q = quantum number to copy
  MomentumQuantumNumber (const MomentumQuantumNumber& Q);

  // destructor
  //
  ~MomentumQuantumNumber ();

  // assignement
  //
  // Q = quantum number to copy
  // return value = reference on current quantum number
  MomentumQuantumNumber& operator = (const MomentumQuantumNumber& Q);

  // set momentum
  //
  // momentum = momentum value
  // return value = reference on current quantum number
  MomentumQuantumNumber& operator = (int momentum);

  // clone current quantum number
  //
  // return value = pointer on cloned quantum number  
  AbstractQuantumNumber* Clone ();

  // Get momentum
  //
  // return value = momentum
  int GetMomentum ();

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

#endif


