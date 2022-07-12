////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of number of particles quantum number               //
//                                                                            //
//                        last modification : 11/05/2001                      //
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


#include "QuantumNumber/NumberParticleQuantumNumber.h"
#include <iostream>


// default constructor
//

NumberParticleQuantumNumber::NumberParticleQuantumNumber () 
{
  this->NumberParticle = 0;
  this->QuantumNumberType = AbstractQuantumNumber::NumberParticle;
}

// constructor from a numberof particles
//
// numberParticle = number of particles

NumberParticleQuantumNumber::NumberParticleQuantumNumber (int numberParticle) 
{
  this->NumberParticle = numberParticle;
  this->QuantumNumberType = AbstractQuantumNumber::NumberParticle;
}

// copy constructor
//
// Q = quantum number to copy

NumberParticleQuantumNumber::NumberParticleQuantumNumber (const NumberParticleQuantumNumber& Q) 
{
  this->NumberParticle = Q.NumberParticle;
  this->QuantumNumberType = Q.QuantumNumberType;
}

// destructor
//

NumberParticleQuantumNumber::~NumberParticleQuantumNumber () 
{
}

// assignement
//
// Q = quantum number to copy
// return value = reference on current quantum number

NumberParticleQuantumNumber& NumberParticleQuantumNumber::operator = (const NumberParticleQuantumNumber& Q) 
{
  this->NumberParticle = Q.NumberParticle;
  this->QuantumNumberType = Q.QuantumNumberType;
  return *this;
}

// clone current quantum number
//
// return value = pointer on cloned quantum number  

AbstractQuantumNumber* NumberParticleQuantumNumber::Clone ()
{
  return new NumberParticleQuantumNumber (this->NumberParticle);
}

// set number of particles
//
// numberParticle = value to assign
// return value = reference on current quantum number

NumberParticleQuantumNumber& NumberParticleQuantumNumber::operator = (int numberParticle) 
{
  this->NumberParticle = numberParticle;
  return *this;
}

// Get number of particles
//
// return value = number of particles

int NumberParticleQuantumNumber::GetNumberParticle () 
{
  return this->NumberParticle;
}

// add a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& NumberParticleQuantumNumber::operator += (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumberType == number.QuantumNumberType)
    this->NumberParticle += ((NumberParticleQuantumNumber&) number).NumberParticle;
  return *this;
}

// substract a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& NumberParticleQuantumNumber::operator -= (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumberType == number.QuantumNumberType)
    this->NumberParticle -= ((NumberParticleQuantumNumber&) number).NumberParticle;
  return *this;
}

// test if two quantum numbers are identical
//
// number = quantum number to compare to the current one

bool NumberParticleQuantumNumber::IsEqual (const AbstractQuantumNumber& number)
{
  if ((this->QuantumNumberType == number.QuantumNumberType) && 
      (this->NumberParticle == ((NumberParticleQuantumNumber&) number).NumberParticle))
    return true;
  else
    return false;
}

// test if two quantum numbers are different
//
// number = quantum number to compare to the current one

bool NumberParticleQuantumNumber::IsDifferent (const AbstractQuantumNumber& number)
{
  if ((this->QuantumNumberType != number.QuantumNumberType) || 
      (this->NumberParticle != ((NumberParticleQuantumNumber&) number).NumberParticle))
    return true;
  else
    return false;
}

// print quantum number
//
// str = reference on current output stream 
// return value = reference on current output stream 

ostream& NumberParticleQuantumNumber::PrintQuantumNumber (ostream& str)
{
  str << "Number of particles = " << this->NumberParticle;
  return str;
}

