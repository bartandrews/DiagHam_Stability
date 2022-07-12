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

#include <iostream>
using std::cout;
using std::endl;

// constructor
//
// lzMax = twice the maximum Lz value reached by a particle
// chirality = flag that allows to choose between either one of two conventions for
// the phase of the orbitals
ParticleOnSphereFunctionBasis::ParticleOnSphereFunctionBasis(int lzMax, int chirality)
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
  this->Chirality = chirality;
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
  double Arg = value[1] * (((double) (this->Chirality*index)) - 0.5 * ((double) (this->Chirality*this->LzMax)));
  result.Re = cos(Arg);
  result.Im = sin(Arg);
  result *= this->Prefactor[index] * pow(cos (0.5 * value[0]), (double) (index)) * pow(sin (0.5 * value[0]), (double) (this->LzMax - index));
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 
// terrible performance, to be used only for testing
void ParticleOnSphereFunctionBasis::GetFunctionValue(Complex U, Complex V, Complex& result, int index)
{
  Complex tmpC=1.0;
  for (int i=0; i<index;++i) tmpC*=U;
  for (int i=0; i<this->LzMax-index;++i) tmpC*=V;
  result = this->Prefactor[index] * tmpC;
}


// get value of the i-th function at a given point (for functions which take values in C)
//
// (theta, phi) = coordinates of point where basis has to be evaluated
// result = reference on the vector with values of the basis
void ParticleOnSphereFunctionBasis::GetAllFunctionValues(Complex U, Complex V, Complex *result)
{
  result[0]=1.0;
  result[1]=U;
  for (int i=2; i<=LzMax; ++i)
    result[i]=U*result[i-1];
  Complex Tmp=1.0;
  result[LzMax]*=this->Prefactor[LzMax];
  for (int index=LzMax-1; index>=0; --index)
    {
      Tmp*=V;
      result[index]*=this->Prefactor[index]*Tmp;
    }
}


// get value of the i-th function at a given point (theta, phi=0), which happens to be real
//
// theta = reference on the value where the function has to be evaluated
// index = function index 

double ParticleOnSphereFunctionBasis::GetRealFunctionValue(double theta, int index)
{
  return this->Prefactor[index] * pow(cos (0.5 * theta), (double) (index)) * pow(sin (0.5 * theta), (double) (this->LzMax - index));
}


void ParticleOnSphereFunctionBasis::GetRealFunctionValues(double theta, double *result)
{
  double st=sin (0.5 * theta);
  double ct=cos (0.5 * theta);
  for (int index=0; index<=LzMax; ++index)
    result[index]=this->Prefactor[index] * pow(ct, (double) (index)) * pow(st, (double) (this->LzMax - index));
}
