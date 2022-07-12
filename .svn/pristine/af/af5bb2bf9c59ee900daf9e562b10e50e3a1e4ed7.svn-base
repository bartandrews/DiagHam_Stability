////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//              class of hilbert space of one 1d periodic box particle        //
//                                                                            //
//                        last modification : 05/07/2004                      //

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
#include "HilbertSpace/PeriodicOneDOneParticle.h"

using std::cout;
using std::endl;

// default constructor
//

PeriodicOneDOneParticle::PeriodicOneDOneParticle ()
{
}

// basic constructor
//
// nbrState = wave function basis dimension
// low = lower impulsion

PeriodicOneDOneParticle::PeriodicOneDOneParticle (int nbrState, int low)
{
  this->NbrState = nbrState;
  this->LowerImpulsion = low;
  this->HilbertSpaceDimension = this->NbrState;
}

// copy constructor
//
// space = reference on Hilbert space to copy

PeriodicOneDOneParticle::PeriodicOneDOneParticle (const PeriodicOneDOneParticle& space)
{
  this->NbrState = space.NbrState;
  this->LowerImpulsion = space.LowerImpulsion;
}


// destructor
//

PeriodicOneDOneParticle::~PeriodicOneDOneParticle ()
{
}


// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

PeriodicOneDOneParticle& PeriodicOneDOneParticle::operator = (const PeriodicOneDOneParticle& space)
{
  this->NbrState = space.NbrState;
  this->LowerImpulsion = space.LowerImpulsion;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* PeriodicOneDOneParticle::Clone()
{
  return new PeriodicOneDOneParticle(*this);
}

