////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of Hilbert space with no state description             //
//                                                                            //
//                        last modification : 07/06/2001                      //
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


#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"

#include <iostream.h>


// default constructor
//

UndescribedHilbertSpace::UndescribedHilbertSpace () 
{
  this->QuantumNumber = 0;
  this->HilbertSpaceDimension = 0;
}

// constructor from datas
//
// dimension = Hilbert space dimension

UndescribedHilbertSpace::UndescribedHilbertSpace (int dimension) 
{
  this->HilbertSpaceDimension = dimension;
  this->QuantumNumber = 0;
}

// constructor from datas
//
// dimension = Hilbert space dimension
// quantumNumber = reference on the quantum number associated to the Hilbert space

UndescribedHilbertSpace::UndescribedHilbertSpace (int dimension, AbstractQuantumNumber& quantumNumber) 
{
  this->HilbertSpaceDimension = dimension;
  this->QuantumNumber = quantumNumber.Clone();
}

// copy constructor (without duplicating datas)
//
// space = reference on Hilbert space to copy

UndescribedHilbertSpace::UndescribedHilbertSpace (const UndescribedHilbertSpace& space) 
{
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;
  if (space.QuantumNumber != 0)
    this->QuantumNumber = space.QuantumNumber->Clone();
  else
    this->QuantumNumber = 0;
}

// destructor
//

UndescribedHilbertSpace::~UndescribedHilbertSpace () 
{
  if (this->QuantumNumber != 0)
    delete this->QuantumNumber;
}

// assignment (without duplicating datas)
//
// space = reference on Hilbert space to copy
// return value = reference on current Hilbert space

UndescribedHilbertSpace& UndescribedHilbertSpace::operator = (const UndescribedHilbertSpace& space)
{
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;
  if (space.QuantumNumber != 0)
    this->QuantumNumber = space.QuantumNumber->Clone();
  else
    this->QuantumNumber = 0;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* UndescribedHilbertSpace::Clone()
{
  return new UndescribedHilbertSpace(*this);
}

// return Hilbert space dimension
//
// return value = Hilbert space dimension

int UndescribedHilbertSpace::GetHilbertSpaceDimension() 
{
  return this->HilbertSpaceDimension;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> UndescribedHilbertSpace::GetQuantumNumbers () 
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->QuantumNumber != 0)
    TmpList += this->QuantumNumber;
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* UndescribedHilbertSpace::GetQuantumNumber (int index) 
{
  return this->QuantumNumber;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* UndescribedHilbertSpace::ExtractSubspace (AbstractQuantumNumber& q, 
								SubspaceSpaceConverter& converter)
{
  if (this->QuantumNumber == 0)
    return 0;
  if ((q != (*(this->QuantumNumber))))
    return 0;
  else
    return this;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& UndescribedHilbertSpace::PrintState (ostream& Str, int state) 
{
  Str << state;
  return Str;
}
