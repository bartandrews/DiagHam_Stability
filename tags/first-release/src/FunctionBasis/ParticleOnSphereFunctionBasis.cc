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


#include "config.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "Vector/RealVector.h"

#include <math.h>


// constructor
//
// lzMax = twice the maximum Lz value reached by a particle

ParticleOnSphereFunctionBasis::ParticleOnSphereFunctionBasis(int lzMax)
{
  this->LzMax = lzMax;
  this->HilbertSpaceDimension = this->LzMax + 1;
  double TmpFactor = ((double) this->HilbertSpaceDimension) / (4.0 * M_PI);
  double TmpBinomial = 1.0;
  this->Prefactor = new double [this->HilbertSpaceDimension];
  this->Prefactor[0] = sqrt (TmpBinomial * TmpFactor);
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      TmpBinomial *= this->LzMax - ((double) i) + 1.0;
      TmpBinomial /= ((double) i);
      this->Prefactor[i] = sqrt (TmpBinomial * TmpFactor);
    }
}

// destructor
//

ParticleOnSphereFunctionBasis::~ParticleOnSphereFunctionBasis ()
{
  delete[] this->Prefactor;
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 

void ParticleOnSphereFunctionBasis::GetFunctionValue(RealVector& value, Complex& result, int index)
{
  double Arg = value[1] * (((double) index) - 0.5 * ((double) (this->LzMax)));
  result.Re = cos(Arg);
  result.Im = sin(Arg);
  result *= this->Prefactor[index] * pow(cos (0.5 * value[0]), (double) (index)) * pow(sin (0.5 * value[0]), (double) (this->LzMax - index));
}



