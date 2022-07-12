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

// add a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& VectorQuantumNumber::operator += (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumbers.GetNbrElement() != ((VectorQuantumNumber&) number).QuantumNumbers.GetNbrElement())
    return *this;
  List<AbstractQuantumNumber*> TmpQ;
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber1(this->QuantumNumbers);
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber2(((VectorQuantumNumber&) number).QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber1;
  AbstractQuantumNumber** TmpQuantumNumber2;
  while (((TmpQuantumNumber1 = IterQuantumNumber1())) && ((TmpQuantumNumber2 = IterQuantumNumber2())))
    {
      (**TmpQuantumNumber1) += (**TmpQuantumNumber2);
    }
  return *this;
}

// substract a quantum nunber to the current one
//
// number = quantum number to add 
// return value = reference to the current quantum number

AbstractQuantumNumber& VectorQuantumNumber::operator -= (const AbstractQuantumNumber& number)
{
  if (this->QuantumNumbers.GetNbrElement() != ((VectorQuantumNumber&) number).QuantumNumbers.GetNbrElement())
    return *this;
  List<AbstractQuantumNumber*> TmpQ;
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber1(this->QuantumNumbers);
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber2(((VectorQuantumNumber&) number).QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber1;
  AbstractQuantumNumber** TmpQuantumNumber2;
  while (((TmpQuantumNumber1 = IterQuantumNumber1())) && ((TmpQuantumNumber2 = IterQuantumNumber2())))
    {
      (**TmpQuantumNumber1) -= (**TmpQuantumNumber2);
    }
  return *this;
}

// test if two quantum numbers are identical
//
// number = quantum number to compare to the current one

bool VectorQuantumNumber::IsEqual (const AbstractQuantumNumber& number)
{
  bool Flag = true;
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber1(this->QuantumNumbers);
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber2(((VectorQuantumNumber &) number).QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber1;
  AbstractQuantumNumber** TmpQuantumNumber2;
  while ((Flag == true) && ((TmpQuantumNumber1 = IterQuantumNumber1())) && ((TmpQuantumNumber2 = IterQuantumNumber2())))
    Flag = ((**TmpQuantumNumber1).IsEqual(**TmpQuantumNumber2));
  return Flag;
}

// test if two quantum numbers are different
//
// number = quantum number to compare to the current one

bool VectorQuantumNumber::IsDifferent (const AbstractQuantumNumber& number)
{
  bool Flag = true;
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber1(this->QuantumNumbers);
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber2(((VectorQuantumNumber &) number).QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber1;
  AbstractQuantumNumber** TmpQuantumNumber2;
  while ((Flag == true) && ((TmpQuantumNumber1 = IterQuantumNumber1())) && ((TmpQuantumNumber2 = IterQuantumNumber2())))
    Flag = (**TmpQuantumNumber1).IsEqual(**TmpQuantumNumber2);
  return !Flag;
}

// print quantum number
//
// str = reference on current output stream 
// return value = reference on current output stream 

ostream& VectorQuantumNumber::PrintQuantumNumber (ostream& str)
{
  ListIterator<AbstractQuantumNumber*> IterQuantumNumber(this->QuantumNumbers);
  AbstractQuantumNumber** TmpQuantumNumber;
  while ((TmpQuantumNumber = IterQuantumNumber()))
    (*TmpQuantumNumber)->PrintQuantumNumber(str) << " ";
  return str;
}

