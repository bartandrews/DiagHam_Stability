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
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/NumberParticleQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "QuantumNumber/MomentumQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"


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

// add two quantum numbers
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = pointer to the sum of the two quantum numbers

AbstractQuantumNumber* operator + (const AbstractQuantumNumber& Q1, 
				   const AbstractQuantumNumber& Q2)
{
  if (Q1.QuantumNumberType != Q2.QuantumNumberType)
    return 0;
  switch (Q1.QuantumNumberType)
    {
    case AbstractQuantumNumber::Sz:
      {
	return new SzQuantumNumber (((SzQuantumNumber&) Q1) + ((SzQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::NumberParticle:
      {
	return new NumberParticleQuantumNumber (((NumberParticleQuantumNumber&) Q1) +
						((NumberParticleQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::Vector:
      {
	return new VectorQuantumNumber(((VectorQuantumNumber&) Q1) + ((VectorQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::Momentum:
      {
	return new  MomentumQuantumNumber (((MomentumQuantumNumber&) Q1) + ((MomentumQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::PeriodicMomentum:
      {
	return new  PeriodicMomentumQuantumNumber (((PeriodicMomentumQuantumNumber&) Q1) + ((PeriodicMomentumQuantumNumber&) Q2));  
      }
      break;
    }
  return 0;
}

// substract two quantum numbers
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = pointer to the difference of the two quantum numbers

AbstractQuantumNumber* operator - (const AbstractQuantumNumber& Q1, 
				   const AbstractQuantumNumber& Q2)
{
  if (Q1.QuantumNumberType != Q2.QuantumNumberType)
    return 0;
  switch (Q1.QuantumNumberType)
    {
    case AbstractQuantumNumber::Sz:
      {
	return new SzQuantumNumber (((SzQuantumNumber&) Q1) - ((SzQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::NumberParticle:
      {
	return new NumberParticleQuantumNumber (((NumberParticleQuantumNumber&) Q1) -
						((NumberParticleQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::Vector:
      {
	return new VectorQuantumNumber(((VectorQuantumNumber&) Q1) - ((VectorQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::Momentum:
      {
	return new  MomentumQuantumNumber (((MomentumQuantumNumber&) Q1) - ((MomentumQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::PeriodicMomentum:
      {
	return new  PeriodicMomentumQuantumNumber (((PeriodicMomentumQuantumNumber&) Q1) - ((PeriodicMomentumQuantumNumber&) Q2));  
      }
      break;
    }
  return 0;
}

// test if two quantum numbers are identical
//
// Q1 = first quantum number
// Q2 = second quantum number

bool operator == (const AbstractQuantumNumber& Q1, const AbstractQuantumNumber& Q2) 
{
  if (Q1.QuantumNumberType != Q2.QuantumNumberType)
    return false;
  switch (Q1.QuantumNumberType)
    {
    case AbstractQuantumNumber::Sz:
      {
	return (((SzQuantumNumber&) Q1) == ((SzQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::NumberParticle:
      {
	return (((NumberParticleQuantumNumber&) Q1) == ((NumberParticleQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::Vector:
      {
	return (((VectorQuantumNumber&) Q1) == ((VectorQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::Momentum:
      {
	return (((MomentumQuantumNumber&) Q1) == ((MomentumQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::PeriodicMomentum:
      {
	return (((PeriodicMomentumQuantumNumber&) Q1) == ((PeriodicMomentumQuantumNumber&) Q2));  
      }
      break;
    }
  return false;
}

// test if two quantum numbers are different
//
// Q1 = first quantum number
// Q2 = second quantum number

bool operator != (const AbstractQuantumNumber& Q1, const AbstractQuantumNumber& Q2) 
{
  if (Q1.QuantumNumberType != Q2.QuantumNumberType)
    return true;
  switch (Q1.QuantumNumberType)
    {
    case AbstractQuantumNumber::Sz:
      {
	return (((SzQuantumNumber&) Q1) != ((SzQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::NumberParticle:
      {
	return (((NumberParticleQuantumNumber&) Q1) != ((NumberParticleQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::Vector:
      {
	return (((VectorQuantumNumber&) Q1) != ((VectorQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::Momentum:
      {
	return (((MomentumQuantumNumber&) Q1) != ((MomentumQuantumNumber&) Q2));
      }
      break;
    case AbstractQuantumNumber::PeriodicMomentum:
      {
	return (((PeriodicMomentumQuantumNumber&) Q1) != ((PeriodicMomentumQuantumNumber&) Q2));  
      }
      break;
    }
  return true;
}

// print quantum number
//
// Str = reference on current output stream 
// Q = quantum number to print
// return value = reference on current output stream 

ostream& operator << (ostream& Str, const AbstractQuantumNumber& Q)
{
  if ((Q.QuantumNumberType & AbstractQuantumNumber::Vector) != 0)
    return (Str << ((VectorQuantumNumber&) Q));  
  switch (Q.QuantumNumberType)
    {
    case AbstractQuantumNumber::Sz:
      {
	return (Str << ((SzQuantumNumber&) Q));  
      }
      break;
    case AbstractQuantumNumber::NumberParticle:
      {
	return (Str << ((NumberParticleQuantumNumber&) Q));  
      }
      break;
    case AbstractQuantumNumber::Momentum:
      {
	return (Str << ((MomentumQuantumNumber&) Q));  
      }
      break;
    case AbstractQuantumNumber::PeriodicMomentum:
      {
	return (Str << ((PeriodicMomentumQuantumNumber&) Q));  
      }
      break;
    }
  return Str;
}

