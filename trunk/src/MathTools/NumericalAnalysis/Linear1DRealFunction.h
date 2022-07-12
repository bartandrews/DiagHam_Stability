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


#ifndef LINEAR1DREALFUNCTION_H
#define LINEAR1DREALFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DRealFunction.h"
#include "GeneralTools/GarbageFlag.h"


class Linear1DRealFunction : public Abstract1DRealFunction
{

 protected:

   // linear coefficient
  double LinearCoefficient;

 // value of the function at zero
  double OffsetValue;

 public:

  // constructor
  //
  // linearCoefficient = linear coefficient
  // offsetValue = value of the function at zerooffsetValue
  Linear1DRealFunction(double linearCoefficient, double offsetValue = 0.0);

  // copy constructor 
  //
  // function = function to copy
  Linear1DRealFunction(const Linear1DRealFunction& function);

  // destructor
  //
  ~Linear1DRealFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DRealFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  double operator ()(double x);

  // get function derivative at a given point
  //
  // x = point where the function derivative has to be evaluated
  // return value = function derivative
  double GetDerivative(const double& x);

  // evaluate derivative of the function on the same interval that the function itself
  //
  // return value = function derivative
  Abstract1DRealFunction* GetDerivative();

  // get function laplacian at a given point
  //
  // x = point where the function laplacian has to be evaluated
  // return value = function laplacian
  double GetLaplacian(const double& x);

  // evaluate laplacian of the function on the same interval that the function itself
  //
  // return value = function laplacian
  Abstract1DRealFunction* GetLaplacian();

  // evaluate integral on the function of a given interval
  // 
  // interval = reference on the interval on which the integral has to be evaluated
  // return value = integral value
  double GetIntegral(AbstractNumericalInterval& interval);

 protected:

    
};

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

inline double Linear1DRealFunction::operator ()(double x)
{
  return ((this->LinearCoefficient * x) + this->OffsetValue);
}



#endif
