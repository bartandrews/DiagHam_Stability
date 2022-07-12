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


#include "QuantumNumber/SzQuantumNumber.h"
#include <iostream>


// default constructor
//

SzQuantumNumber::SzQuantumNumber () 
{
  this->Sz = 0;
  this->QuantumNumberType = AbstractQuantumNumber::Sz;
}

// constructor from a Sz value
//
// sz = twice the value of  Sz component

SzQuantumNumber::SzQuantumNumber (int sz) 
{
  this->Sz = sz;
  this->QuantumNumberType = AbstractQuantumNumber::Sz;
}

// copy constructor
//
// Q = quantum number to copy

SzQuantumNumber::SzQuantumNumber (const SzQuantumNumber& Q) 
{
  this->Sz = Q.Sz;
  this->QuantumNumberType = Q.QuantumNumberType;
}

// destructor
//

SzQuantumNumber::~SzQuantumNumber () 
{
}

// assignement
//
// Q = quantum number to copy
// return value = reference on current quantum number

SzQuantumNumber& SzQuantumNumber::operator = (const SzQuantumNumber& Q) 
{
  this->Sz = Q.Sz;
  this->QuantumNumberType = Q.QuantumNumberType;
  return *this;
}

// clone current quantum number
//
// return value = pointer on cloned quantum number  

AbstractQuantumNumber* SzQuantumNumber::Clone ()
{
  return new SzQuantumNumber(this->Sz);
}

// set total Sz component
//
// sz = value to assign
// return value = reference on current quantum number

SzQuantumNumber& SzQuantumNumber::operator = (int sz) 
{
  this->Sz = sz;
  return *this;
}

// Get total Sz component
//
// return value = total Sz component

int SzQuantumNumber::GetSz () 
{
  return this->Sz;
}

// add two quantum numbers
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = sum of the two quantum numbers

SzQuantumNumber operator + (const SzQuantumNumber& Q1, const SzQuantumNumber& Q2) 
{
  return SzQuantumNumber(Q1.Sz + Q2.Sz);
}

// substract two quantum numbers
//
// Q1 = first quantum number
// Q2 = quantum number to substract
// return value = sum of the two quantum numbers

SzQuantumNumber operator - (const SzQuantumNumber& Q1, const SzQuantumNumber& Q2) 
{
  return SzQuantumNumber(Q1.Sz - Q2.Sz);
}

// test if two quantum numbers are identical
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = true if quantum numbers are identical

bool operator == (const SzQuantumNumber& Q1, const SzQuantumNumber& Q2) 
{
  if (Q1.Sz == Q2.Sz)
    return true;
  else
    return false;
}

// test if two quantum numbers are different
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = true if quantum numbers are different

bool operator != (const SzQuantumNumber& Q1, const SzQuantumNumber& Q2) 
{
  if (Q1.Sz != Q2.Sz)
    return true;
  else
    return false;
}

// print quantum number
//
// Str = reference on current output stream 
// Q = quantum number to print
// return value = reference on current output stream 

ostream& operator << (ostream& Str, const SzQuantumNumber& Q)
{
  Str << "Sz = " << Q.Sz;
  return Str;
}

