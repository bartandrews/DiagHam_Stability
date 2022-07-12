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


#include "QuantumNumber/AbstractQuantumNumber.h"


// default constructor
//

AbstractQuantumNumber::AbstractQuantumNumber () 
{
  this->QuantumNumberType = 0;
}

// copy constructor
//
// Q = quantum number to copy

AbstractQuantumNumber::AbstractQuantumNumber (const AbstractQuantumNumber& Q) 
{
  this->QuantumNumberType = Q.QuantumNumberType;
}

// virtual destructor
//

AbstractQuantumNumber::~AbstractQuantumNumber () 
{
}

// assignement (without duplicating datas)
//
// Q = quantum number to copy
// return value = reference on current quantum number

AbstractQuantumNumber& AbstractQuantumNumber::operator = (const AbstractQuantumNumber& Q) 
{
  this->QuantumNumberType = Q.QuantumNumberType;
  return *this;
}

// clone current quantum number
//
// return value = pointer on cloned quantum number

AbstractQuantumNumber* AbstractQuantumNumber::Clone () 
{
  return new AbstractQuantumNumber ();
}

// get quantum number type
//
// return value = Hilbert space dimension

int AbstractQuantumNumber::GetQuantumNumberType() 
{
  return this->QuantumNumberType;
}

// add a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& AbstractQuantumNumber::operator += (const AbstractQuantumNumber& number)
{
  return *this;
}

// substract a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& AbstractQuantumNumber::operator -= (const AbstractQuantumNumber& number)
{
  return *this;
}

// add a quantum nunber to the current one and return a new quantum number (not modifying the current one)
//
// number = quantum number to add to the current one
  // return value = pointer on cloned quantum number corresponding to the sum 

AbstractQuantumNumber* AbstractQuantumNumber::Add (const AbstractQuantumNumber& number)
{
  return this->Clone();
}

// substract a quantum nunber to the current one and return a new quantum number (not modifying the current one)
//
// number = quantum number to substract to the current one
// return value = pointer on cloned quantum number corresponding to the difference 

AbstractQuantumNumber* AbstractQuantumNumber::Sub (const AbstractQuantumNumber& number)
{
  return this->Clone();
}

// test if two quantum numbers are identical
//
// number = quantum number to compare to the current one

bool AbstractQuantumNumber::IsEqual (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumberType == number.QuantumNumberType)
    return true;
  else
    return false;
}

// test if two quantum numbers are different
//
// number = quantum number to compare to the current one

bool AbstractQuantumNumber::IsDifferent (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumberType != number.QuantumNumberType)
    return true;
  else
    return false;
}

// print quantum number
//
// str = reference on current output stream 
// return value = reference on current output stream 

ostream& AbstractQuantumNumber::PrintQuantumNumber (ostream& str)
{
  return str;
}

