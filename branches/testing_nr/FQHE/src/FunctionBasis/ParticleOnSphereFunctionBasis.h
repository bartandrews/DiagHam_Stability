////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of function basis for particle on sphere             //
//                                                                            //
//                        last modification : 10/12/2002                      //
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


#ifndef PARTICLEONSPHEREFUNCTIONBASIS_H
#define PARTICLEONSPHEREFUNCTIONBASIS_H


#include "config.h"
#include "FunctionBasis/AbstractFunctionBasis.h"


class ParticleOnSphereFunctionBasis: public AbstractFunctionBasis
{

 protected:

  // twice the maximum Lz value reached by a particle
  int LzMax;
  
  // array containing numerical prefactor of each function
  double* Prefactor;

  // Chirality of the phase
  int Chirality;

 public:

  enum
    {
      RightHanded = +1,
      LeftHanded = -1
    };

  // constructor
  //
  // lzMax = twice the maximum Lz value reached by a particle
  // chirality = flag that allows to choose between either one of two conventions for
  // the phase of the orbitals
  ParticleOnSphereFunctionBasis(int lzMax, int chirality=ParticleOnSphereFunctionBasis::RightHanded);

  // destructor
  //
  ~ParticleOnSphereFunctionBasis ();

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // value = reference on the value where the function has to be evaluated
  // result = reference on the value where the result has to be stored
  // index = function index 
  void GetFunctionValue(RealVector& value, Complex& result, int index);

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // value = reference on the value where the function has to be evaluated
  // result = reference on the value where the result has to be stored
  // index = function index 
  // terrible performance, to be used only for testing
  void GetFunctionValue(Complex U, Complex V, Complex& result, int index);

  // get value of the i-th function at a given point (for functions which take values in C)
  //
  // (U, V) = spinor coordinates of point where basis has to be evaluated
  // result = reference on the vector with values of the basis
  void GetAllFunctionValues(Complex U, Complex V, Complex *result);

  // get value of the i-th function at a given point (theta, phi=0), which happens to be real
  //
  // theta = reference on the value where the function has to be evaluated
  // index = function index 

  double GetRealFunctionValue(double theta, int index);

  // calculate them all:
  void GetRealFunctionValues(double theta, double *result);


};

#endif


