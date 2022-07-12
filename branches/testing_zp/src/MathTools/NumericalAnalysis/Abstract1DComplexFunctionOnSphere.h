////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of abstract 1D complex function                   //
//                                                                            //
//                        last modification : 01/09/2004                      //
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


#ifndef ABSTRACT1DCOMPLEXFUNCTIONONSPHERE_H
#define ABSTRACT1DCOMPLEXFUNCTIONONSPHERE_H


#include "config.h"
#include "Abstract1DComplexFunction.h"
#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"

class RealVector;


class Abstract1DComplexFunctionOnSphere : public Abstract1DComplexFunction
{

 public:

  // virtual destructor
  //
  virtual ~Abstract1DComplexFunctionOnSphere();

  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DComplexFunction* Clone () = 0;

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  virtual Complex operator ()(RealVector& x);

  // evaluate function at a given point
  //
  // uv = ensemble of spinor variables on sphere describing point
  //      where function has to be evaluated
  //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
  // return value = function value at (uv)
  virtual Complex CalculateFromSpinorVariables(ComplexVector& uv) = 0;

  // get function properties, and possible extensions of interface 
  // 
  virtual unsigned GetProperties();

};

#endif
