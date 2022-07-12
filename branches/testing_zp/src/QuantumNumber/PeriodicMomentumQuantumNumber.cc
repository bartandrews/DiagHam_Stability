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

// add a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& PeriodicMomentumQuantumNumber::operator += (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumberType == number.QuantumNumberType)
    {
      this->Momentum += ((PeriodicMomentumQuantumNumber&) number).Momentum;
      this->Momentum %= this->Period;
    }
  return *this;
}

// substract a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& PeriodicMomentumQuantumNumber::operator -= (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumberType == number.QuantumNumberType)
    {
      this->Momentum -= ((PeriodicMomentumQuantumNumber&) number).Momentum;
      this->Momentum %= this->Period;
    }
  return *this;
}

// test if two quantum numbers are identical
//
// number = quantum number to compare to the current one

bool PeriodicMomentumQuantumNumber::IsEqual (const AbstractQuantumNumber& number)
{
  if ((this->QuantumNumberType == number.QuantumNumberType) && 
      (this->Momentum == ((PeriodicMomentumQuantumNumber&) number).Momentum))
    return true;
  else
    return false;
}

// test if two quantum numbers are different
//
// number = quantum number to compare to the current one

bool PeriodicMomentumQuantumNumber::IsDifferent (const AbstractQuantumNumber& number)
{
  if ((this->QuantumNumberType != number.QuantumNumberType) || 
      (this->Momentum != ((PeriodicMomentumQuantumNumber&) number).Momentum))
    return true;
  else
    return false;
}

// print quantum number
//
// str = reference on current output stream 
// return value = reference on current output stream 

ostream& PeriodicMomentumQuantumNumber::PrintQuantumNumber (ostream& str)
{
  str << "P = " << this->Momentum;
  return str;
}

