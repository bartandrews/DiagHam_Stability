////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//             class of hilbert space of one periodic 3d particle             //
//                                                                            //
//                        last modification : 11/08/2004                      //
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
#include "HilbertSpace/PeriodicThreeDOneParticle.h"

using std::cout;
using std::endl;

// default constructor
//

PeriodicThreeDOneParticle::PeriodicThreeDOneParticle ()
{
}

// constructor
//
// nbrStateX = wave function basis dimension in the x direction
// lowX = lower impulsion in X direction (in unit of 2 * Pi / Lx)
// nbrStateY = wave function basis dimension in the y direction
// lowY = lower impulsion in Y direction (in unit of 2 * Pi / Ly)
// nbrStateZ = wave function basis dimension in the z direction
// lowZ = lower impulsion in Z direction (in unit of 2 * Pi / Lz)

PeriodicThreeDOneParticle::PeriodicThreeDOneParticle (int nbrStateX, int lowX, int nbrStateY, int lowY, int nbrStateZ, int lowZ)
{
  this->StateX = new PeriodicOneDOneParticle (nbrStateX, lowX);
  this->StateY = new PeriodicOneDOneParticle (nbrStateY, lowY);
  this->StateZ = new PeriodicOneDOneParticle (nbrStateZ, lowZ);
  this->HilbertSpaceDimension = this->StateX->GetNbrState () * this->StateY->GetNbrState () * this->StateZ->GetNbrState ();
}

// copy constructor
//
// space = reference on Hilbert space to copy

PeriodicThreeDOneParticle::PeriodicThreeDOneParticle (const PeriodicThreeDOneParticle& space)
{
  this->StateX = (PeriodicOneDOneParticle*) space.StateX;
  this->StateY = (PeriodicOneDOneParticle*) space.StateY;
  this->StateZ = (PeriodicOneDOneParticle*) space.StateZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
}

// destructor
//

PeriodicThreeDOneParticle::~PeriodicThreeDOneParticle ()
{
  delete (PeriodicOneDOneParticle*) this->StateX;
  delete (PeriodicOneDOneParticle*) this->StateY;
  delete (PeriodicOneDOneParticle*) this->StateZ;
  this->StateX = NULL; this->StateY = NULL; this->StateZ = NULL; 
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

PeriodicThreeDOneParticle& PeriodicThreeDOneParticle::operator = (const PeriodicThreeDOneParticle& space)
{
  this->StateX = (PeriodicOneDOneParticle*) space.StateX;
  this->StateY = (PeriodicOneDOneParticle*) space.StateY;
  this->StateZ = (PeriodicOneDOneParticle*) space.StateZ;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* PeriodicThreeDOneParticle::Clone ()
{
  return new PeriodicThreeDOneParticle (*this);
}

