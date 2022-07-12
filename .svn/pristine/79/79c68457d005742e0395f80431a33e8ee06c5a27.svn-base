////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of vector of quantum numbers                    //
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


#include "QuantumNumber/VectorQuantumNumber.h"
#include "GeneralTools/ListIterator.h"
#include <iostream>


// default constructor
//

VectorQuantumNumber::VectorQuantumNumber () 
{
  this->QuantumNumberType = AbstractQuantumNumber::Vector;
}

// constructor from a list of quantum numbers
//
// quantumNumbers = list of quantum numbers

VectorQuantumNumber::VectorQuantumNumber(List<AbstractQuantumNumber*> quantumNumbers)
{
  this->QuantumNumbers = quantumNumbers;
  this->QuantumNumberType = AbstractQuantumNumber::Vector;  
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber(quantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber;
  while ((TmpQuantumNumber = IterQuantumNumber()))
    this->QuantumNumberType |= (*TmpQuantumNumber)->GetQuantumNumberType();
}

// copy constructor
//
// Q = quantum number to copy

VectorQuantumNumber::VectorQuantumNumber (const VectorQuantumNumber& Q) 
{
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber(Q.QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber;
  while ((TmpQuantumNumber = IterQuantumNumber()))
    this->QuantumNumbers += (*TmpQuantumNumber)->Clone();
  this->QuantumNumberType = Q.QuantumNumberType;
}

// destructor
//

VectorQuantumNumber::~VectorQuantumNumber () 
{
}

// assignement
//
// Q = quantum number to copy
// return value = reference on current quantum number

VectorQuantumNumber& VectorQuantumNumber::operator = (const VectorQuantumNumber& Q) 
{
  AbstractQuantumNumber** TmpQuantumNumber;
  if (this->QuantumNumbers.GetNbrElement() != 0)
    {
      ListIterator<AbstractQuantumNumber*> IterQuantumNumber2(this->QuantumNumbers);
      while ((TmpQuantumNumber = IterQuantumNumber2()))
	delete  (*TmpQuantumNumber);
      this->QuantumNumbers.DeleteList();
    }
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber(Q.QuantumNumbers);
  while ((TmpQuantumNumber = IterQuantumNumber()))
    this->QuantumNumbers += (*TmpQuantumNumber)->Clone();
  this->QuantumNumberType = Q.QuantumNumberType;
  return *this;
}

// clone current quantum number
//
// return value = pointer on cloned quantum number  

AbstractQuantumNumber* VectorQuantumNumber::Clone ()
{
  List<AbstractQuantumNumber*> TmpList;
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber(this->QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber;
  while ((TmpQuantumNumber = IterQuantumNumber()))
    TmpList += (*TmpQuantumNumber)->Clone();
  return new VectorQuantumNumber (TmpList);
}

// add a quantum number
//
// numberParticle = value to assign
// return value = reference on current quantum number

VectorQuantumNumber& VectorQuantumNumber::operator += (AbstractQuantumNumber* quantumNumber)
{
  this->QuantumNumbers += quantumNumber;
  this->QuantumNumberType |= quantumNumber->GetQuantumNumberType();
  return *this;
}

// Set i-th component of vector of quantum numbers
//
// index = index of the component to set 
// quantumNumber = new quantum number
// return value = reference on current quantum number

VectorQuantumNumber& VectorQuantumNumber::SetQuantumNumber(int index, AbstractQuantumNumber* quantumNumber)
{
  delete this->QuantumNumbers[index];
  this->QuantumNumbers[index] = quantumNumber;
  return *this;
}

// Get i-th component of vector of quantum numbers
//
// index = index of the component to get 
// return value = pointer to the quantum number

AbstractQuantumNumber* VectorQuantumNumber::operator [] (int index)
{
  return this->QuantumNumbers[index];
}

// add two quantum numbers
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = sum of the two quantum numbers

VectorQuantumNumber operator + (VectorQuantumNumber& Q1, VectorQuantumNumber& Q2) 
{
  if (Q1.QuantumNumbers.GetNbrElement() != Q2.QuantumNumbers.GetNbrElement())
    return VectorQuantumNumber();
  List<AbstractQuantumNumber*> TmpQ;
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber1(Q1.QuantumNumbers);
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber2(Q2.QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber1;
  AbstractQuantumNumber** TmpQuantumNumber2;
  while (((TmpQuantumNumber1 = IterQuantumNumber1())) && ((TmpQuantumNumber2 = IterQuantumNumber2())))
    {
      TmpQ += ((**TmpQuantumNumber1) + (**TmpQuantumNumber2));
    }
  return VectorQuantumNumber(TmpQ);
}

// substract two quantum numbers
//
// Q1 = first quantum number
// Q2 = quantum number to substract
// return value = sum of the two quantum numbers

VectorQuantumNumber operator - (VectorQuantumNumber& Q1, VectorQuantumNumber& Q2) 
{
  if (Q1.QuantumNumbers.GetNbrElement() != Q2.QuantumNumbers.GetNbrElement())
    return VectorQuantumNumber();
  List<AbstractQuantumNumber*> TmpQ;
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber1(Q1.QuantumNumbers);
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber2(Q2.QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber1;
  AbstractQuantumNumber** TmpQuantumNumber2;
  while (((TmpQuantumNumber1 = IterQuantumNumber1())) && ((TmpQuantumNumber2 = IterQuantumNumber2())))
    {
      TmpQ += ((**TmpQuantumNumber1) - (**TmpQuantumNumber2));
    }
  return VectorQuantumNumber(TmpQ);
}

// test if two quantum numbers are identical
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = true if quantum numbers are identical

bool operator == (const VectorQuantumNumber& Q1, const VectorQuantumNumber& Q2) 
{
  bool Flag = true;
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber1(Q1.QuantumNumbers);
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber2(Q2.QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber1;
  AbstractQuantumNumber** TmpQuantumNumber2;
  while ((Flag == true) && ((TmpQuantumNumber1 = IterQuantumNumber1())) && ((TmpQuantumNumber2 = IterQuantumNumber2())))
  if ((**TmpQuantumNumber1) != (**TmpQuantumNumber2))
    Flag = false;
  return Flag;
}

// test if two quantum numbers are different
//
// Q1 = first quantum number
// Q2 = second quantum number
// return value = true if quantum numbers are different

bool operator != (const VectorQuantumNumber& Q1, const VectorQuantumNumber& Q2) 
{
  bool Flag = true;
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber1(Q1.QuantumNumbers);
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber2(Q2.QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber1;
  AbstractQuantumNumber** TmpQuantumNumber2;
  while ((Flag == true) && ((TmpQuantumNumber1 = IterQuantumNumber1())) && ((TmpQuantumNumber2 = IterQuantumNumber2())))
  if ((**TmpQuantumNumber1) != (**TmpQuantumNumber2))
    Flag = false;
  return !Flag;
}

// print quantum number
//
// Str = reference on current output stream 
// Q = quantum number to print
// return value = reference on current output stream 

ostream& operator << (ostream& Str, const VectorQuantumNumber& Q)
{
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber(Q.QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber;
  while ((TmpQuantumNumber = IterQuantumNumber()))
    Str << (**TmpQuantumNumber) << " ";
  return Str;
}

