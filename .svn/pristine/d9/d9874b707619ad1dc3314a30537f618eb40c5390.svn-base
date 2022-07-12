////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of sum of tensor product states for bosons on sphere       //  
//                                                                            //
//                        last modification : 19/03/2011                      //
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
#include "HilbertSpace/BosonOnSphereSumTensorProductState.h"
#include "Vector/RealVector.h"


// default constructor
// 
BosonOnSphereSumTensorProductState::BosonOnSphereSumTensorProductState()
{
  this->NbrTensorProductsPerSum = 0;
  this->TensorProductStates = 0;
  this->TensorProductSpaces = 0;
  this->Paddings = 0;
  this->Coefficients = 0;
}

// constructor
// 
// nbrTensorProducts = number of states that makes the tensor product state
// tensorProductStates = states that makes the tensor product state
// tensorProductSpaces = Hilbert space associated to each state
// paddings = orbital shifts between two successive Hilbert spaces
// coefficient = an optional global multiplicative factor

BosonOnSphereSumTensorProductState::BosonOnSphereSumTensorProductState(int nbrSumTensorProducts, int* nbrTensorProducts, RealVector** tensorProductStates, BosonOnSphereShort*** tensorProductSpaces, int** paddings, double* coefficients)
{
  this->NbrTensorProducts = nbrSumTensorProducts;
  this->NbrTensorProductsPerSum = nbrTensorProducts;
  this->TensorProductStates = tensorProductStates;
  this->TensorProductSpaces = tensorProductSpaces;
  this->Paddings = paddings;
  this->Coefficients = coefficients;
}

// destructor
//
 
BosonOnSphereSumTensorProductState::~BosonOnSphereSumTensorProductState()
{
  if (this->TensorProductSpaces != 0)
    {
      for (int i = 0; i < this->NbrTensorProducts; ++i)
	{
	  delete[] this->TensorProductSpaces[i];
	  delete[] this->Paddings[i];
	}
      delete[] this->TensorProductStates;
      delete[] this->NbrTensorProductsPerSum;
      delete[] this->Paddings;
    }
}

// check if the tensor product state has a component  on a given occupation number
//
// state = array describing the occupation number 
// return value = true if the tensor product state has a component on the given occupation number

bool BosonOnSphereSumTensorProductState::Belong(unsigned long* state)
{
  return false;
}

// check if the tensor product state has a component on a given occupation number and return its value
//
// state = array describing the occupation number 
// return value = components (0 if the tensor product state has no component on the given occupation number)
double BosonOnSphereSumTensorProductState::GetComponent(unsigned long* state)
{
  return 0.0;
}

// convert the state associated to the fused basis into a vector defined in a full basis
//
// targetBasis = pointer to the full basis
// outputVector = reference on the vector where the state will be stored
// return value = true if no error occured
bool BosonOnSphereSumTensorProductState::ConvertToFullBasis(BosonOnSphereShort* targetBasis, RealVector& outputVector)
{
  RealVector TmpVector (outputVector.GetLargeVectorDimension(), true);
  for (int i = 0; i < this->NbrTensorProducts; ++i)
    {
      targetBasis->FuseMultipleStates(TmpVector, this->NbrTensorProductsPerSum[i], this->TensorProductStates[i], this->Paddings[i], (ParticleOnSphere**) (this->TensorProductSpaces[i]), false, this->Coefficients[i]);
      outputVector += TmpVector;
      TmpVector.ClearVector();
    }
  return true;
}
