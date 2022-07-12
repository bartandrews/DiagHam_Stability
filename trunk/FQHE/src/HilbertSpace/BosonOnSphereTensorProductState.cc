////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of tensor product state for bosons on sphere           //  
//                                                                            //
//                        last modification : 16/03/2011                      //
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
#include "HilbertSpace/BosonOnSphereTensorProductState.h"

  
// default constructor
// 

BosonOnSphereTensorProductState::BosonOnSphereTensorProductState()
{
  this->NbrTensorProducts = 0;
  this->TensorProductStates = 0;
  this->TensorProductSpaces = 0;
  this->Coefficient = 0.0;
  this->Paddings = 0;
}

// constructor
// 
// nbrTensorProducts = number of states that makes the tensor product state
// tensorProductStates = states that makes the tensor product state
// tensorProductSpaces = Hilbert space associated to each state
// paddings = orbital shifts between two successive Hilbert spaces
// coefficient = an optional global multiplicative factor

BosonOnSphereTensorProductState::BosonOnSphereTensorProductState(int nbrTensorProducts, RealVector* tensorProductStates, BosonOnSphereShort** tensorProductSpaces, int* paddings, double coefficient)
{
  this->NbrTensorProducts = nbrTensorProducts;
  this->TensorProductStates = tensorProductStates;
  this->TensorProductSpaces = tensorProductSpaces;
  this->Paddings = paddings;
  this->Coefficient = coefficient;
}

// destructor
// 

BosonOnSphereTensorProductState::~BosonOnSphereTensorProductState()
{
  if (this->TensorProductSpaces == 0)
    {
      delete[] this->TensorProductSpaces;
      delete[] this->Paddings;  
    }
}

// check if the tensor product state has a component  on a given occupation number
//
// state = array describing the occupation number 
// return value = true if the tensor product state has a component on the given occupation number

bool BosonOnSphereTensorProductState::Belong(unsigned long* state)
{
//   for (int i = 0; i < this->NbrTensorProducts; ++i)
//     {
//       int TmpLz = 0;
//       for (int j = 0; < 
//     }
//   return true;
  return false;
}

// check if the tensor product state has a component on a given occupation number and return its value
//
// state = array describing the occupation number 
// return value = components (0 if the tensor product state has no component on the given occupation number)

double BosonOnSphereTensorProductState::GetComponent(unsigned long* state)
{
  return 0.0;
}

// convert the state associated to the fused basis into a vector defined in a full basis
//
// targetBasis = pointer to the full basis
// outputVector = reference on the vector where the state will be stored
// return value = true if no error occured

bool BosonOnSphereTensorProductState::ConvertToFullBasis(BosonOnSphereShort* targetBasis, RealVector& outputVector)
{
  targetBasis->FuseMultipleStates(outputVector, this->NbrTensorProducts, this->TensorProductStates, this->Paddings, (ParticleOnSphere**) this->TensorProductSpaces, false, this->Coefficient);
  return true;
}

