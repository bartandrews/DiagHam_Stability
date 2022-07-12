////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//                     class of hilbert space of two 3d particles             //
//                                                                            //
//                        last modification : 13/10/2004                      //
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


#include "config.h"
#include "HilbertSpace/ThreeDTwoParticles.h"

using std::cout;
using std::endl;

// default constructor
//

ThreeDTwoParticles::ThreeDTwoParticles ()
{
}

// constructor from two abstract particles
//
// firstParticle = pointer to the first abstract particle
// secondParticle = pointer to the second abstract particle

ThreeDTwoParticles::ThreeDTwoParticles (ThreeDOneParticle* firstParticle, ThreeDOneParticle* secondParticle)
{
  this->FirstParticle = firstParticle;
  this->SecondParticle = secondParticle;
  this->HilbertSpaceDimension = this->FirstParticle->GetHilbertSpaceDimension () * this->SecondParticle->GetHilbertSpaceDimension ();
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// copy constructor
//
// space = reference on Hilbert space to copy

ThreeDTwoParticles::ThreeDTwoParticles (const ThreeDTwoParticles& space)
{
  this->FirstParticle = space.FirstParticle;
  this->SecondParticle = space.SecondParticle;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

ThreeDTwoParticles::~ThreeDTwoParticles ()
{
  delete this->FirstParticle;
  delete this->SecondParticle;
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

ThreeDTwoParticles& ThreeDTwoParticles::operator = (const ThreeDTwoParticles& space)
{
  this->FirstParticle = space.FirstParticle;
  this->SecondParticle = space.SecondParticle;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* ThreeDTwoParticles::Clone ()
{
  return new ThreeDTwoParticles (*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> ThreeDTwoParticles::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* ThreeDTwoParticles::GetQuantumNumber (int index)
{
  return 0;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* ThreeDTwoParticles::ExtractSubspace (AbstractQuantumNumber& q, 
							      SubspaceSpaceConverter& converter)
{
  return 0;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& ThreeDTwoParticles::PrintState (ostream& Str, int state)
{ 
  int dimension2 = this->SecondParticle->GetHilbertSpaceDimension ();
  int state1 = state / dimension2; int state2 = state - state1 * dimension2;

  Str << "'" << this->FirstParticle->PrintState (Str, state1) << ", " << this->SecondParticle->PrintState (Str, state2) << ")";
  return Str;
}
