////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of exact state wave function                    //
//                                                                            //
//                        last modification : 20/04/2005                      //
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
#include "Tools/FQHEWaveFunction/ExactWaveFunction.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "HilbertSpace/AbstractQHEParticle.h"

#ifdef USE_HILBERT_SPACE

// constructor
//
// components = vector that describes exact state components in ExactStateSpace basis
// exactStateSpace = Hilbert space associated to the exact state
// oneBodyBasis = one body real space basis to use 

ExactWaveFunction::ExactWaveFunction(RealVector& components, AbstractQHEParticle* exactStateSpace, AbstractFunctionBasis* oneBodyBasis)
{
  this->ExactState = components;
  this->ExactStateSpace = (AbstractQHEParticle*) (exactStateSpace->Clone());
  this->OneBodyBasis = oneBodyBasis;
  this->ExactStateIndex = -1;
}

// constructor for a one component exact state
//
// index = index that describes the one component exact state
// exactStateSpace = Hilbert space associated to the exact state
// oneBodyBasis = one body real space basis to use 

ExactWaveFunction::ExactWaveFunction(int index, AbstractQHEParticle* exactStateSpace, AbstractFunctionBasis* oneBodyBasis)
{
  this->ExactStateSpace = (AbstractQHEParticle*) (exactStateSpace->Clone());
  this->OneBodyBasis = oneBodyBasis;
  this->ExactStateIndex = index;
  this->ExactState = RealVector (this->ExactStateSpace->GetHilbertSpaceDimension(), true);
  if ((this->ExactStateIndex >= 0) && (this->ExactStateIndex < this->ExactStateSpace->GetHilbertSpaceDimension()))
    this->ExactState[this->ExactStateIndex] = 1.0;
}

// copy constructor
//
// function = reference on the wave function to copy

ExactWaveFunction::ExactWaveFunction(const ExactWaveFunction& function)
{
  this->ExactState = function.ExactState;
  this->OneBodyBasis = function.OneBodyBasis;
  this->ExactStateIndex = function.ExactStateIndex;
  this->ExactStateSpace = (AbstractQHEParticle*) (function.ExactStateSpace->Clone());
}

// destructor
//

ExactWaveFunction::~ExactWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* ExactWaveFunction::Clone ()
{
  return new ExactWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex ExactWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp;
  if (this->ExactStateIndex == -1)
    {
      Tmp = this->ExactStateSpace->EvaluateWaveFunction(this->ExactState, x, *(this->OneBodyBasis));
    }
  else
    {
      Tmp = this->ExactStateSpace->EvaluateWaveFunction(this->ExactState, x, *(this->OneBodyBasis), this->ExactStateIndex, 1);
    }
  return Tmp;
}

#endif
