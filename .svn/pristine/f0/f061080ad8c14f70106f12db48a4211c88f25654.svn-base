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


#ifndef ABSTRACTFUNCTIONBASIS_H
#define ABSTRACTFUNCTIONBASIS_H


#include "config.h"
#include "MathTools/Complex.h"


class Vector;
class RealVector;
class ComplexVector;


class AbstractFunctionBasis
{

 protected:

  // dimension of the hilbert space spawned by the function basis
  int HilbertSpaceDimension;

 public:

  // virtual destructor
  //
  virtual ~AbstractFunctionBasis ();

  // return Hilbert space dimension
  //
  // return value = Hilbert space dimension
  virtual int GetHilbertSpaceDimension();

  // get value of the i-th function at a given point (for functions which take values in R)
  //
  // value = reference on the value where the function has to be evaluated
  // result = reference on the value where the result has to be stored
  // index = function index 
  virtual void GetFunctionValue(Vector& value, double& result, int index);

  // get value of the i-th function at a given point (for functions which take values in R)
  //
  // value = reference on the value where the function has to be evaluated
  // result = reference on the value where the result has to be stored
  // index = function index 
  virtual void GetFunctionValue(RealVector& value, double& result, int index);

  // get value of the i-th function at a given point (for functions which take values in R)
  //
  // value = reference on the value where the function has to be evaluated
  // result = reference on the value where the result has to be stored
  // index = function index 
  virtual void GetFunctionValue(ComplexVector& value, double& result, int index);

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // value = reference on the value where the function has to be evaluated
  // result = reference on the value where the result has to be stored
  // index = function index 
  virtual void GetFunctionValue(Vector& value, Complex& result, int index);

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // value = reference on the value where the function has to be evaluated
  // result = reference on the value where the result has to be stored
  // index = function index 
  virtual void GetFunctionValue(RealVector& value, Complex& result, int index);

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // value = reference on the value where the function has to be evaluated
  // result = reference on the value where the result has to be stored
  // index = function index 
  virtual void GetFunctionValue(ComplexVector& value, Complex& result, int index);

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // x, y = coordinate where the function has to be evaluated
  // index = function index 
  virtual Complex GetFunctionValue(double x, double y, int index);

};

// return Hilbert space dimension
//
// return value = Hilbert space dimension

inline int AbstractFunctionBasis::GetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}

#endif


