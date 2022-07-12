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

#ifndef BOSONONSPHERESUMTENSORPRODUCTSTATE
#define BOSONONSPHERESUMTENSORPRODUCTSTATE

#include "config.h"
#include "HilbertSpace/BosonOnSphereTensorProductState.h"
#include "Vector/RealVector.h"


class BosonOnSphereSumTensorProductState : public BosonOnSphereTensorProductState
{

 protected:

  // number of states that makes each tensor product state
  int* NbrTensorProductsPerSum;
  // states that makes the tensor product state
  RealVector** TensorProductStates;
  // Hilbert space associated to each state
  BosonOnSphereShort*** TensorProductSpaces;
  // orbital shifts between two successive Hilbert spaces
  int** Paddings;
  // global multiplicative factor per tensor product state 
  double* Coefficients;

 public:

  // default constructor
  // 
  BosonOnSphereSumTensorProductState();

  // constructor
  // 
  // nbrTensorProducts = number of states that makes the tensor product state
  // tensorProductStates = states that makes the tensor product state
  // tensorProductSpaces = Hilbert space associated to each state
  // paddings = orbital shifts between two successive Hilbert spaces
  // coefficient = an optional global multiplicative factor
  BosonOnSphereSumTensorProductState(int NbrSumTensorProducts, int* nbrTensorProducts, RealVector** tensorProductStates, BosonOnSphereShort*** tensorProductSpaces, int** paddings, double* coefficients);

  // destructor
  // 
  virtual ~BosonOnSphereSumTensorProductState();

  // check if the tensor product state has a component  on a given occupation number
  //
  // state = array describing the occupation number 
  // return value = true if the tensor product state has a component on the given occupation number
  virtual bool Belong(unsigned long* state);

  // check if the tensor product state has a component on a given occupation number and return its value
  //
  // state = array describing the occupation number 
  // return value = components (0 if the tensor product state has no component on the given occupation number)
  virtual double GetComponent(unsigned long* state);

  // convert the state associated to the fused basis into a vector defined in a full basis
  //
  // targetBasis = pointer to the full basis
  // outputVector = reference on the vector where the state will be stored
  // return value = true if no error occured
  virtual bool ConvertToFullBasis(BosonOnSphereShort* targetBasis, RealVector& outputVector);

};

#endif
