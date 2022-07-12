////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of linear 1D real function                     //
//                                                                            //
//                        last modification : 11/08/2015                      //
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
#include "MathTools/NumericalAnalysis/Linear1DRealFunction.h"


// constructor
//
// linearCoefficient = linear coefficient
// offsetValue = value of the function at zerooffsetValue

Linear1DRealFunction::Linear1DRealFunction(double linearCoefficient, double offsetValue)
{
  this->LinearCoefficient = linearCoefficient;
  this->OffsetValue = offsetValue;
}

// copy constructor 
//
// function = function to copy

Linear1DRealFunction::Linear1DRealFunction(const Linear1DRealFunction& function)
{
  this->LinearCoefficient = function.LinearCoefficient;
  this->OffsetValue = function.OffsetValue;
}

// destructor
//

Linear1DRealFunction::~Linear1DRealFunction()
{ 
}

// clone function 
//
// return value = clone of the function 

Abstract1DRealFunction* Linear1DRealFunction::Clone ()
{
  return new Linear1DRealFunction(*this);
}

// get function derivative at a given point
//
// x = point where the function derivative has to be evaluated
// return value = function derivative

double Linear1DRealFunction::GetDerivative(const double& x)
{
  return 0.0;
}

// evaluate derivative of the function on the same interval that the function itself
//
// return value = function derivative

Abstract1DRealFunction* Linear1DRealFunction::GetDerivative()
{
  return 0;
}

// get function laplacian at a given point
//
// x = point where the function laplacian has to be evaluated
// return value = function laplacian

double Linear1DRealFunction::GetLaplacian(const double& x)
{
  return 0.0;
}

// evaluate laplacian of the function on the same interval that the function itself
//
// return value = function laplacian

Abstract1DRealFunction* Linear1DRealFunction::GetLaplacian()
{
  return 0;
}

// evaluate integral on the function of a given interval
// 
// interval = reference on the interval on which the integral has to be evaluated
// return value = integral value

double Linear1DRealFunction::GetIntegral(AbstractNumericalInterval& interval)
{
  return 0.0;
}

