////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//class of hilbert space of one 3d periodic box particle with reflexion symmetry of X axis //
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
#include "HilbertSpace/PeriodicXReflexionYZPeriodicThreeDOneParticle.h"


// default constructor
//

PeriodicXReflexionYZPeriodicThreeDOneParticle::PeriodicXReflexionYZPeriodicThreeDOneParticle ()
{
}

// constructor
//
// nbrStateX = wave function basis dimension in the x direction without symmetry redundancy
// evenX = true if the wave functions in X is even, else odd
// nbrStateY = wave function basis dimension in the y direction
// lowY = lower impulsion in Y direction (in unit of 2 * Pi / Lz)
// nbrStateZ = wave function basis dimension in the z direction
// lowZ = lower impulsion in Z direction (in unit of 2 * Pi / Lz)

PeriodicXReflexionYZPeriodicThreeDOneParticle::PeriodicXReflexionYZPeriodicThreeDOneParticle (int nbrStateX, bool evenX, int nbrStateY, int lowY, int nbrStateZ, int lowZ)
{
  this->StateX = new PeriodicReflexionSymmetryOneDOneParticle (nbrStateX, evenX);
  this->StateY = new PeriodicOneDOneParticle (nbrStateY, lowY);
  this->StateZ = new PeriodicOneDOneParticle (nbrStateZ, lowZ);
  this->HilbertSpaceDimension = this->StateX->GetNbrState () * this->StateY->GetNbrState () * this->StateZ->GetNbrState ();
}

// copy constructor
//
// space = reference on Hilbert space to copy

PeriodicXReflexionYZPeriodicThreeDOneParticle::PeriodicXReflexionYZPeriodicThreeDOneParticle (const PeriodicXReflexionYZPeriodicThreeDOneParticle& space)
{
  this->StateX = (PeriodicReflexionSymmetryOneDOneParticle*) space.StateX;
  this->StateY = (PeriodicOneDOneParticle*) space.StateY;
  this->StateZ = (PeriodicOneDOneParticle*) space.StateZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
}

// destructor
//

PeriodicXReflexionYZPeriodicThreeDOneParticle::~PeriodicXReflexionYZPeriodicThreeDOneParticle ()
{
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

PeriodicXReflexionYZPeriodicThreeDOneParticle& PeriodicXReflexionYZPeriodicThreeDOneParticle::operator = (const PeriodicXReflexionYZPeriodicThreeDOneParticle& space)
{
  this->StateX = (PeriodicReflexionSymmetryOneDOneParticle*) space.StateX;
  this->StateY = (PeriodicOneDOneParticle*) space.StateY;
  this->StateZ = (PeriodicOneDOneParticle*) space.StateZ;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* PeriodicXReflexionYZPeriodicThreeDOneParticle::Clone ()
{
  return new PeriodicXReflexionYZPeriodicThreeDOneParticle (*this);
}

