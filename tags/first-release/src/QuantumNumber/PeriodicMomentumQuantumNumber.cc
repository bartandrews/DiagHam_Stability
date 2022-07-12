////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of periodic momentum quantum number                 //
//                                                                            //
//                        last modification : 08/10/2002                      //
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


#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include <iostream>


// default constructor
//

PeriodicMomentumQuantumNumber::PeriodicMomentumQuantumNumber () 
{
  this->Momentum = 0;
  this->Period = 1;
  this->QuantumNumberType = AbstractQuantumNumber::PeriodicMomentum;
}

// constructor from a momentum value
//
// momentum = momentum value

PeriodicMomentumQuantumNumber::PeriodicMomentumQuantumNumber (int momentum, int period) 
{
  this->Momentum = momentum;
  this->Period = period;
  this->QuantumNumberType = AbstractQuantumNumber::PeriodicMomentum;

}

// copy constructor
//
// Q = quantum number to copy

PeriodicMomentumQuantumNumber::PeriodicMomentumQuantumNumber (const PeriodicMomentumQuantumNumber& Q) 
{
  this->Momentum = Q.Momentum;
  this->Period = Q.Period;
  this->QuantumNumberType = Q.QuantumNumberType;
}

// destructor
//

PeriodicMomentumQuantumNumber::~PeriodicMomentumQuantumNumber () 
{
}

// assignement
//
// Q = quantum number to copy
// return value = reference on current quantum number

PeriodicMomentumQuantumNumber& PeriodicMomentumQuantumNumber::operator = (const PeriodicMomentumQuantumNumber& Q) 
{
  this->Momentum = Q.Momentum;
  this->Period = Q.Period;
  this->QuantumNumberType = Q.QuantumNumberType;
  return *this;
}

// clone current quantum number
//
// return value = pointer on cloned quantum number  

AbstractQuantumNumber* PeriodicMomentumQuantumNumber::Clone ()
{
  return new PeriodicMomentumQuantumNumber(*this);
}

// set momentum
//
// momentum = momentum value
// return value = reference on current quantum number

PeriodicMomentumQuantumNumber& PeriodicMomentumQuantumNumber::operator = (int momentum) 
{
  this->Momentum = momentum;
  return *this;
}

// Get momentum
//
// return value = momentum

int PeriodicMomentumQuantumNumber::GetMomentum () 
{
  return this->Momentum;
}

// add two quantum numbers
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = sum of the two quantum numbers

PeriodicMomentumQuantumNumber operator + (const PeriodicMomentumQuantumNumber& Q1, const PeriodicMomentumQuantumNumber& Q2) 
{
  return PeriodicMomentumQuantumNumber((Q1.Momentum + Q2.Momentum) % Q1.Period, Q1.Period); 
}

// substract two quantum numbers
//
// Q1 = first quantum number
// Q2 = quantum number to substract
// return value = sum of the two quantum numbers

PeriodicMomentumQuantumNumber operator - (const PeriodicMomentumQuantumNumber& Q1, const PeriodicMomentumQuantumNumber& Q2) 
{
  return PeriodicMomentumQuantumNumber((Q1.Momentum - Q2.Momentum) % Q1.Period, Q1.Period);
}

// test if two quantum numbers are identical
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = true if quantum numbers are identical

bool operator == (const PeriodicMomentumQuantumNumber& Q1, const PeriodicMomentumQuantumNumber& Q2) 
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

bool operator != (const PeriodicMomentumQuantumNumber& Q1, const PeriodicMomentumQuantumNumber& Q2) 
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

ostream& operator << (ostream& Str, const PeriodicMomentumQuantumNumber& Q)
{
  Str << "P = " << Q.Momentum;
  return Str;
}

