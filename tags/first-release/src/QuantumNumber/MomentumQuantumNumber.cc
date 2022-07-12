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


#include "QuantumNumber/MomentumQuantumNumber.h"
#include <iostream>


// default constructor
//

MomentumQuantumNumber::MomentumQuantumNumber () 
{
  this->Momentum = 0;
  this->QuantumNumberType = AbstractQuantumNumber::Momentum;
}

// constructor from a momentum value
//
// momentum = momentum value

MomentumQuantumNumber::MomentumQuantumNumber (int momentum) 
{
  this->Momentum = momentum;
  this->QuantumNumberType = AbstractQuantumNumber::Momentum;
}

// copy constructor
//
// Q = quantum number to copy

MomentumQuantumNumber::MomentumQuantumNumber (const MomentumQuantumNumber& Q) 
{
  this->Momentum = Q.Momentum;
  this->QuantumNumberType = Q.QuantumNumberType;
}

// destructor
//

MomentumQuantumNumber::~MomentumQuantumNumber () 
{
}

// assignement
//
// Q = quantum number to copy
// return value = reference on current quantum number

MomentumQuantumNumber& MomentumQuantumNumber::operator = (const MomentumQuantumNumber& Q) 
{
  this->Momentum = Q.Momentum;
  this->QuantumNumberType = Q.QuantumNumberType;
  return *this;
}

// clone current quantum number
//
// return value = pointer on cloned quantum number  

AbstractQuantumNumber* MomentumQuantumNumber::Clone ()
{
  return new MomentumQuantumNumber(this->Momentum);
}

// set momentum
//
// momentum = momentum value
// return value = reference on current quantum number

MomentumQuantumNumber& MomentumQuantumNumber::operator = (int momentum) 
{
  this->Momentum = momentum;
  return *this;
}

// Get momentum
//
// return value = momentum

int MomentumQuantumNumber::GetMomentum () 
{
  return this->Momentum;
}

// add two quantum numbers
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = sum of the two quantum numbers

MomentumQuantumNumber operator + (const MomentumQuantumNumber& Q1, const MomentumQuantumNumber& Q2) 
{
  return MomentumQuantumNumber(Q1.Momentum + Q2.Momentum);
}

// substract two quantum numbers
//
// Q1 = first quantum number
// Q2 = quantum number to substract
// return value = sum of the two quantum numbers

MomentumQuantumNumber operator - (const MomentumQuantumNumber& Q1, const MomentumQuantumNumber& Q2) 
{
  return MomentumQuantumNumber(Q1.Momentum - Q2.Momentum);
}

// test if two quantum numbers are identical
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = true if quantum numbers are identical

bool operator == (const MomentumQuantumNumber& Q1, const MomentumQuantumNumber& Q2) 
{
  if (Q1.Momentum == Q2.Momentum)
    return true;
  else
    return false;
}

// test if two quantum numbers are different
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = true if quantum numbers are different

bool operator != (const MomentumQuantumNumber& Q1, const MomentumQuantumNumber& Q2) 
{
  if (Q1.Momentum != Q2.Momentum)
    return true;
  else
    return false;
}

// print quantum number
//
// Str = reference on current output stream 
// Q = quantum number to print
// return value = reference on current output stream 

ostream& operator << (ostream& Str, const MomentumQuantumNumber& Q)
{
  Str << "P = " << Q.Momentum;
  return Str;
}

