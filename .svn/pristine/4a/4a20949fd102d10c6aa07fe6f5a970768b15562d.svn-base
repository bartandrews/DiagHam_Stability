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

// add a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& SzQuantumNumber::operator += (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumberType == number.QuantumNumberType)
    this->Sz += ((SzQuantumNumber&) number).Sz;
  return *this;
}

// substract a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& SzQuantumNumber::operator -= (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumberType == number.QuantumNumberType)
    this->Sz -= ((SzQuantumNumber&) number).Sz;
  return *this;
}

// test if two quantum numbers are identical
//
// number = quantum number to compare to the current one

bool SzQuantumNumber::IsEqual (const AbstractQuantumNumber& number)
{
  if ((this->QuantumNumberType == number.QuantumNumberType) && 
      (this->Sz == ((SzQuantumNumber&) number).Sz))
    return true;
  else
    return false;
}

// test if two quantum numbers are different
//
// number = quantum number to compare to the current one

bool SzQuantumNumber::IsDifferent (const AbstractQuantumNumber& number)
{
  if ((this->QuantumNumberType != number.QuantumNumberType) || 
      (this->Sz != ((SzQuantumNumber&) number).Sz))
    return true;
  else
    return false;
}

// print quantum number
//
// str = reference on current output stream 
// return value = reference on current output stream 

ostream& SzQuantumNumber::PrintQuantumNumber (ostream& str)
{
  str << "Sz = " << this->Sz;
  return str;
}

