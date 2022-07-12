////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
//class of hilbert space of one particle with z periodicity and planar rotation symmetry//
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
#include "HilbertSpace/PlanarRotationSymmetryZPeriodicOneParticle.h"
#include "HilbertSpace/PeriodicOneDOneParticle.h"


// constructor
//
// lz = quantum number of kinetic momentum in Z direction
// nbrStateR = number of states in plane 
// nbrStateZ = wave function basis dimension in the z direction
// lowerImpulsionZ = lower impulsion in Z direction (in unit of 2 * Pi / Lz)

PlanarRotationSymmetryZPeriodicOneParticle::PlanarRotationSymmetryZPeriodicOneParticle (int lz, int nbrStateR, int nbrStateZ, int lowerImpulsionZ)
{
  this->Lz = lz;
  this->StateR = new OneDOneParticle (nbrStateR);
  this->StateZ = new PeriodicOneDOneParticle (nbrStateZ, lowerImpulsionZ);
  this->HilbertSpaceDimension = this->StateR->GetNbrState () * this->StateZ->GetNbrState ();
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// copy constructor
//
// space = reference on Hilbert space to copy

PlanarRotationSymmetryZPeriodicOneParticle::PlanarRotationSymmetryZPeriodicOneParticle(const PlanarRotationSymmetryZPeriodicOneParticle& space)
{
  this->Lz = space.Lz;
  this->StateR = space.StateR;
  this->StateZ = space.StateZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

PlanarRotationSymmetryZPeriodicOneParticle::~PlanarRotationSymmetryZPeriodicOneParticle()
{
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

PlanarRotationSymmetryZPeriodicOneParticle& PlanarRotationSymmetryZPeriodicOneParticle::operator = (const PlanarRotationSymmetryZPeriodicOneParticle& space)
{
  this->Lz = space.Lz;
  this->StateR = space.StateR;
  this->StateZ = space.StateZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* PlanarRotationSymmetryZPeriodicOneParticle::Clone()
{
  return new PlanarRotationSymmetryZPeriodicOneParticle(*this);
}
