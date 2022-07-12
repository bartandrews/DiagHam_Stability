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

// add a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& MomentumQuantumNumber::operator += (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumberType == number.QuantumNumberType)
    this->Momentum += ((MomentumQuantumNumber&) number).Momentum;
  return *this;
}

// substract a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& MomentumQuantumNumber::operator -= (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumberType == number.QuantumNumberType)
    this->Momentum -= ((MomentumQuantumNumber&) number).Momentum;
  return *this;
}

// test if two quantum numbers are identical
//
// number = quantum number to compare to the current one

bool MomentumQuantumNumber::IsEqual (const AbstractQuantumNumber& number)
{
  if ((this->QuantumNumberType == number.QuantumNumberType) && 
      (this->Momentum == ((MomentumQuantumNumber&) number).Momentum))
    return true;
  else
    return false;
}

// test if two quantum numbers are different
//
// number = quantum number to compare to the current one

bool MomentumQuantumNumber::IsDifferent (const AbstractQuantumNumber& number)
{
  if ((this->QuantumNumberType != number.QuantumNumberType) || 
      (this->Momentum != ((MomentumQuantumNumber&) number).Momentum))
    return true;
  else
    return false;
}

// print quantum number
//
// str = reference on current output stream 
// return value = reference on current output stream 

ostream& MomentumQuantumNumber::PrintQuantumNumber (ostream& str)
{
  str << "P = " << this->Momentum;
  return str;
}

