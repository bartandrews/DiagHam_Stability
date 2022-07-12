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


#ifndef ABSTRACT1DREALFUNCTION_H
#define ABSTRACT1DREALFUNCTION_H


#include "config.h"


class AbstractNumericalInterval;


class Abstract1DRealFunction
{

 public:

  // virtual destructor
  //
  virtual ~Abstract1DRealFunction();

  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DRealFunction* Clone () = 0;

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  virtual double operator ()(double x);

  // get function derivative at a given point
  //
  // x = point where the function derivative has to be evaluated
  // return value = function derivative
  virtual double GetDerivative(const double& x);

  // evaluate derivative of the function on the same interval that the function itself
  //
  // return value = function derivative
  virtual Abstract1DRealFunction* GetDerivative();

  // get function laplacian at a given point
  //
  // x = point where the function laplacian has to be evaluated
  // return value = function laplacian
  virtual double GetLaplacian(const double& x);

  // evaluate laplacian of the function on the same interval that the function itself
  //
  // return value = function laplacian
  virtual Abstract1DRealFunction* GetLaplacian();

  // evaluate integral on the function of a given interval
  // 
  // interval = reference on the interval on which the integral has to be evaluated
  // return value = integral value
  virtual double GetIntegral(AbstractNumericalInterval& interval);
    
};

#endif
