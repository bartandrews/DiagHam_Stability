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


#include "config.h"
#include "FunctionBasis/ParticleOnDiskFunctionBasis.h"
#include "Vector/RealVector.h"

#include <math.h>


// constructor
//
// lzMax = twice the maximum Lz value reached by a particle

ParticleOnDiskFunctionBasis::ParticleOnDiskFunctionBasis(int lzMax)
{
  this->LzMax = lzMax;
  this->HilbertSpaceDimension = this->LzMax + 1;
  this->Prefactor = new double [this->HilbertSpaceDimension];
  this->Prefactor[0] = sqrt (1.0 / (2.0 * M_PI));
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      this->Prefactor[i] = this->Prefactor[i - 1] / sqrt (2.0 * ((double) i));
    }
}

// destructor
//

ParticleOnDiskFunctionBasis::~ParticleOnDiskFunctionBasis ()
{
  delete[] this->Prefactor;
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 

void ParticleOnDiskFunctionBasis::GetFunctionValue(RealVector& value, Complex& result, int index)
{
  result = this->Prefactor[index];
  Complex Tmp(value[0], value[1]);
  while (index > 0)
    {
      result *= Tmp;
      --index;
    }
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// x, y = coordinate where the function has to be evaluated
// index = function index 

Complex ParticleOnDiskFunctionBasis::GetFunctionValue(double x, double y, int index)
{
    Complex result = this->Prefactor[index];
    Complex Tmp(x, y);
    while (index > 0)
    {
        result *= Tmp;
        --index;
    }
    return result;
}

