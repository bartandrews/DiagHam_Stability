////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2003-2004 Duc-Phuong Nguyen                 //
//                                                                            //
//                                                                            //
// class of hilbert space of one particle in 3d with planar rotation symmetry //
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
#include "HilbertSpace/PlanarRotationSymmetryThreeDOneParticle.h"


// default constructor
//

PlanarRotationSymmetryThreeDOneParticle::PlanarRotationSymmetryThreeDOneParticle ()
{
}

// constructor
//
// lz = quantum number of kinetic momentum in Z direction
// nbrStateR = wave function basis dimension of the radial coordinate
// nbrStateZ = wave function basis dimension in the z direction

PlanarRotationSymmetryThreeDOneParticle::PlanarRotationSymmetryThreeDOneParticle(int lz, int nbrStateR, int nbrStateZ)
{
  this->Lz = lz;
  this->StateR = new OneDOneParticle (nbrStateR);
  this->StateZ = new OneDOneParticle (nbrStateZ);
  this->HilbertSpaceDimension = this->StateR->GetNbrState () * this->StateZ->GetNbrState ();
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// copy constructor
//
// space = reference on Hilbert space to copy

PlanarRotationSymmetryThreeDOneParticle::PlanarRotationSymmetryThreeDOneParticle(const PlanarRotationSymmetryThreeDOneParticle& space)
{
  this->Lz = space.Lz;
  this->StateR = space.StateR;
  this->StateZ = space.StateZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

PlanarRotationSymmetryThreeDOneParticle::~PlanarRotationSymmetryThreeDOneParticle()
{
  delete this->StateR;
  delete this->StateZ;
}

// assignement
//
// space = reference on Hilbert space to assign
// return value = reference on current Hilbert space

PlanarRotationSymmetryThreeDOneParticle& PlanarRotationSymmetryThreeDOneParticle::operator = (const PlanarRotationSymmetryThreeDOneParticle& space)
{
  this->Lz = space.Lz;
  this->StateR = space.StateR;
  this->StateZ = space.StateZ;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;  
  return *this;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* PlanarRotationSymmetryThreeDOneParticle::Clone()
{
  return new PlanarRotationSymmetryThreeDOneParticle(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> PlanarRotationSymmetryThreeDOneParticle::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* PlanarRotationSymmetryThreeDOneParticle::GetQuantumNumber (int index)
{
  return 0;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* PlanarRotationSymmetryThreeDOneParticle::ExtractSubspace (AbstractQuantumNumber& q, 
							      SubspaceSpaceConverter& converter)
{
  return 0;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& PlanarRotationSymmetryThreeDOneParticle::PrintState (ostream& Str, int state)
{
  int NbrStateZ = this->StateZ->GetNbrState ();
  int i1 = state / NbrStateZ;
  Str << "(" << this->Lz << ", " << i1 << ", " << (state - (i1 * NbrStateZ)) << ")";
  return Str;
}
