////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of decribing subspace                       //
//                                                                            //
//                        last modification : 26/04/2001                      //
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


#include "DMRGAlgorithm/Subspace.h"


// default constructor
//

Subspace::Subspace () 
{
  this->SubspaceDimension = 0;
  this->QuantumNumber = 0;
  this->VectorsInSubspace = RealMatrix();
  this->Hamiltonian = RealDiagonalMatrix();
  this->Converter = SubspaceSpaceConverter();
}

// constructor for a not yet defined subspace coresponding to a given quantum number
// 
// Q = pointer to quantum number

Subspace::Subspace (AbstractQuantumNumber* Q) 
{
  this->SubspaceDimension = 0;
  this->QuantumNumber = Q;
  this->VectorsInSubspace = RealMatrix();
  this->Hamiltonian = RealDiagonalMatrix();
  this->Converter = SubspaceSpaceConverter();
}

// copy constructor
//
// subspace = subspace to copy

Subspace::Subspace (const Subspace& subspace) 
{
  this->SubspaceDimension = subspace.SubspaceDimension;
  this->QuantumNumber = subspace.QuantumNumber;
  this->Converter = subspace.Converter;
  this->Hamiltonian = subspace.Hamiltonian;
  this->VectorsInSubspace = subspace.VectorsInSubspace;
}

// destructor
//

Subspace::~Subspace() 
{
}

// assignment 
//
// subspace = subspace to assign
// return value = reference to current subspace

Subspace& Subspace::operator = (const Subspace& subspace) 
{
  this->SubspaceDimension = subspace.SubspaceDimension;
  this->QuantumNumber = subspace.QuantumNumber;
  this->Converter = subspace.Converter;
  this->Hamiltonian = subspace.Hamiltonian;
  this->VectorsInSubspace = subspace.VectorsInSubspace;
  return *this;
}

// get subspace dimension
//
// return value = subspace dimension

int Subspace::GetSubspaceDimension () 
{
  return this->Converter.GetSubspaceDimension();
}

// get reference to associated subspace-space converter
//
// rerturn value = reference to converter

SubspaceSpaceConverter& Subspace::GetConverter () 
{
  return this->Converter;
}

// get reference to diagonalized hamiltonian associated to the current subspace
//
// return value = reference to hamiltonian

RealDiagonalMatrix& Subspace::GetHamiltonian () 
{
  return this->Hamiltonian;
}
  
// get reference to Transformation Matrix from eigenvector base to canonical subspace base 
//
// return value = reference to coresponding matrix

RealMatrix& Subspace::GetMatrix () 
{
  return this->VectorsInSubspace;
}
