////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of tabulated 1D real function                    //
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


#ifndef TABULATED1DREALFUNCTION_H
#define TABULATED1DREALFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DRealFunction.h"



class AbstractNumericalInterval;
class AbstractSubdividedInterval;


class Tabulated1DRealFunction : public Abstract1DRealFunction
{

 protected:

  // pointer to the interval description
  AbstractSubdividedInterval* Interval;

  // array containing function tabulated values 
  double* TabulatedValues;

 public:

  // constructor from raw datas
  //
  // interval = interval on which the function is defined (including the subdivision scheme)
  Tabulated1DRealFunction(AbstractSubdividedInterval* interval, double* tabulatedValues = 0);

  // constructor from a C function
  //
  // interval = interval on which the function is defined (including the subdivision scheme)
  // function = pointer to the C function that described the mathematical function that will be used to initialized the tabulated function
  Tabulated1DRealFunction(AbstractSubdividedInterval* interval, double (*function) (double));

  // copy constructor 
  //
  // function = function to copy
  Tabulated1DRealFunction(const Tabulated1DRealFunction& function);

  // destructor
  //
  ~Tabulated1DRealFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DRealFunction* Clone ();

  // get reference on a function tabulated value
  //
  // i = index of the tabulated value
  // return value = reference on the tabulated value
  double& operator [](unsigned long i);

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

  // use linear interpolation to obtain function value at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value obatined using linear interpolation
  double ValueFromLinearInterpolation (const double& x);
  
  // use spline interpolation to obtain function value at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value obatined using spline interpolation
  double ValueFromSplineInterpolation (const double& x);
  
    
};

// get reference on a function tabulated value
//
// i = index of the tabulated value
// return value = reference on the tabulated value

inline double& operator [](unsigned long i)
{
  return this->TabulatedValues[i];
}

#endif
