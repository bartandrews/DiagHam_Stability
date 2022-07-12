////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of constant 1D real function                    //
//                                                                            //
//                        last modification : 08/08/2015                      //
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


#ifndef CONSTANT1DREALFUNCTION_H
#define CONSTANT1DREALFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DRealFunction.h"
#include "GeneralTools/GarbageFlag.h"


class Constant1DRealFunction : public Abstract1DRealFunction
{

 protected:

  // value of the function
  double ConstantValue;

 public:

  // constructor
  //
  // constantValue  = value of the function
  Constant1DRealFunction(double constantValue);

  // copy constructor 
  //
  // function = function to copy
  Constant1DRealFunction(const Constant1DRealFunction& function);

  // destructor
  //
  ~Constant1DRealFunction();

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

inline double Constant1DRealFunction::operator ()(double x)
{
  return this->ConstantValue;
}


#endif
