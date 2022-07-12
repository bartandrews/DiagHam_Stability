////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of abstract function basis                     //
//                                                                            //
//                        last modification : 09/12/2002                      //
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
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"


// virtual destructor
//

AbstractFunctionBasis::~AbstractFunctionBasis ()
{
}

// get value of the i-th function at a given point (for functions which take values in R)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 

void AbstractFunctionBasis::GetFunctionValue(Vector& value, double& result, int index)
{
  if (value.GetVectorType() == Vector::RealDatas)
    return this->GetFunctionValue((RealVector&) value, result, index);
  if (value.GetVectorType() == Vector::ComplexDatas)
    return this->GetFunctionValue((ComplexVector&) value, result, index);
}

// get value of the i-th function at a given point (for functions which take values in R)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 

void AbstractFunctionBasis::GetFunctionValue(RealVector& value, double& result, int index)
{
  result = 0.0;
  return;
}

// get value of the i-th function at a given point (for functions which take values in R)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 

void AbstractFunctionBasis::GetFunctionValue(ComplexVector& value, double& result, int index)
{
  result = 0.0;
  return;
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 

void AbstractFunctionBasis::GetFunctionValue(Vector& value, Complex& result, int index)
{
  if (value.GetVectorType() == Vector::RealDatas)
    return this->GetFunctionValue((RealVector&) value, result, index);
  if (value.GetVectorType() == Vector::ComplexDatas)
    return this->GetFunctionValue((ComplexVector&) value, result, index);
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 

void AbstractFunctionBasis::GetFunctionValue(RealVector& value, Complex& result, int index)
{
  result.Re = 0.0;
  result.Im = 0.0;
  return;
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 

void AbstractFunctionBasis::GetFunctionValue(ComplexVector& value, Complex& result, int index)
{
  result.Re = 0.0;
  result.Im = 0.0;
  return;
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// x, y = coordinate where the function has to be evaluated
// index = function index 
Complex AbstractFunctionBasis::GetFunctionValue(double x, double y, int index)
{
  Complex result;
  result.Re = 0.0;
  result.Im = 0.0;
  return result;
}
