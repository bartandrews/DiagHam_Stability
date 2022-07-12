////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of function basis for particle on disk              //
//                         (without the gaussian factor)                      //
//                                                                            //
//                        last modification : 05/02/2004                      //
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


#ifndef PARTICLEONDISKFUNCTIONBASIS_H
#define PARTICLEONDISKFUNCTIONBASIS_H


#include "config.h"
#include "FunctionBasis/AbstractFunctionBasis.h"


class ParticleOnDiskFunctionBasis: public AbstractFunctionBasis
{

 protected:

  // maximum Lz value reached by a particle
  int LzMax;

  // array containing numerical prefactor of each function
  double* Prefactor;

 public:

  // constructor
  //
  // lzMax = twice the maximum Lz value reached by a particle
  ParticleOnDiskFunctionBasis(int lzMax);

  // destructor
  //
  ~ParticleOnDiskFunctionBasis ();

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // value = reference on the value where the function has to be evaluated
  // result = reference on the value where the result has to be stored
  // index = function index 
  void GetFunctionValue(RealVector& value, Complex& result, int index);

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // x, y = coordinate where the function has to be evaluated
  // index = function index 
  virtual Complex GetFunctionValue(double x, double y, int index);

};

#endif


