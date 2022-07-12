////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            diagham  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of abstract Hilbert space                      //
//                                                                            //
//                        last modification : 23/04/2001                      //
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


#include "HilbertSpace/AbstractHilbertSpace.h"


// default constructor
//
AbstractHilbertSpace::AbstractHilbertSpace ()
{
  LargeHilbertSpaceDimension=0x0l;
}

// virtual destructor
//

AbstractHilbertSpace::~AbstractHilbertSpace () 
{
}


// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id  

int AbstractHilbertSpace::GetHilbertSpaceAdditionalSymmetry()
{
  return 0; // universal coding for no symmetry
}


// write the entire Hilbert-space basis to the given stream
//
// Str = stream to write on

ostream& AbstractHilbertSpace::ShowBasis (ostream& Str)
{
  for (int i=0; i<this->GetHilbertSpaceDimension(); ++i)
    this->PrintState(Str,i)<<std::endl;
  return Str;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> AbstractHilbertSpace::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* AbstractHilbertSpace::GetQuantumNumber (int index)
{
  return 0;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* AbstractHilbertSpace::ExtractSubspace (AbstractQuantumNumber& q, 
							     SubspaceSpaceConverter& converter)
{
  return 0;
}
