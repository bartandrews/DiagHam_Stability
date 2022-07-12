////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of abstract 1D function                       //
//                                                                            //
//                        last modification : 08/07/2004                      //
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
#include "MathTools/NumericalAnalysis/Abstract1DRealFunction.h"


// destructor
//

Abstract1DRealFunction::~Abstract1DRealFunction()
{
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

double Abstract1DRealFunction::operator ()(double x)
{
  return 0.0;
}

// get function derivative at a given point
//
// x = point where the function derivative has to be evaluated
// return value = function derivative

double Abstract1DRealFunction::GetDerivative(const double& x)
{
  return 0.0;
}

// evaluate derivative of the function on the same interval that the function itself
//
// return value = function derivative

Abstract1DRealFunction* Abstract1DRealFunction::GetDerivative()
{
  return 0;
}

// get function laplacian at a given point
//
// x = point where the function laplacian has to be evaluated
// return value = function laplacian

double Abstract1DRealFunction::GetLaplacian(const double& x)
{
  return 0.0;
}

// evaluate laplacian of the function on the same interval that the function itself
//
// return value = function laplacian

Abstract1DRealFunction* Abstract1DRealFunction::GetLaplacian()
{
  return 0;
}

// evaluate integral on the function of a given interval
// 
// interval = reference on the interval on which the integral has to be evaluated
// return value = integral value

double Abstract1DRealFunction::GetIntegral(AbstractNumericalInterval& interval)
{
  return 0.0;
}
