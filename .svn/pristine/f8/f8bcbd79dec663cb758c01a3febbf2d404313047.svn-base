////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                          Class author Cecile Repellin                      //
//                                                                            //
//                 class of function basis for particle on CP2                //
//                                                                            //
//                        last modification : 07/02/2013                      //
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
#include "FunctionBasis/ParticleOnCP2FunctionBasis.h"
#include "MathTools/FactorialCoefficient.h"
#include "Vector/RealVector.h"

#include <math.h>

#include <iostream>
using std::cout;
using std::endl;

// constructor
//
// nbrFluxQuanta = number of flux quanta
// chirality = flag that allows to choose between either one of two conventions for
// the phase of the orbitals
ParticleOnCP2FunctionBasis::ParticleOnCP2FunctionBasis(int nbrFluxQuanta, int chirality)
{
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->HilbertSpaceDimension = (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 2)/2;
  this->LzMax = this->HilbertSpaceDimension - 1;
  this->quantumNumberR = new int [this->HilbertSpaceDimension];
  this->quantumNumberS = new int [this->HilbertSpaceDimension];
  this->GetQuantumNumbersFromLinearizedIndex(this->quantumNumberR, this->quantumNumberS);
  FactorialCoefficient Coef;
  Coef.SetToOne();
  Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
  double TmpFactor = sqrt(Coef.GetNumericalValue())/(M_SQRT_2* M_PI);
  this->Prefactor = new double [this->HilbertSpaceDimension];
  for (int index = 0; index <= this->LzMax; ++index)
  {
    Coef.SetToOne();
    Coef.FactorialDivide(this->quantumNumberR[index]);
    Coef.FactorialDivide(this->quantumNumberS[index]);
    Coef.FactorialDivide(this->NbrFluxQuanta - this->quantumNumberR[index] - this->quantumNumberS[index]);
    this->Prefactor[index] = TmpFactor * sqrt(Coef.GetNumericalValue());
    }
  this->Chirality = chirality;
}

// destructor
//

ParticleOnCP2FunctionBasis::~ParticleOnCP2FunctionBasis ()
{
  delete[] this->quantumNumberR;
  delete[] this->quantumNumberS;
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 

void ParticleOnCP2FunctionBasis::GetFunctionValue(RealVector& value, Complex& result, int index)
{
  int r = this->quantumNumberR[index];
  int s = this->quantumNumberS[index];
  int t = this->NbrFluxQuanta - this->quantumNumberR[index] - this->quantumNumberS[index];
  double Arg = value[2] * (double) (this->Chirality * s) + value[3] * (double) (this->Chirality * t);
  result.Re = cos(Arg);
  result.Im = sin(Arg);
  result *= this->Prefactor[index] * pow(cos (value[0]), (double) (r)) * pow(sin(value[0]) * cos(value[1]), (double) (s)) * pow( sin(value[0]) * sin(value[1]), (double) (t));
}


// get value of the i-th function at a given point (theta, phi=0), which happens to be real
//
// theta = reference on the value where the function has to be evaluated
// index = function index 

double ParticleOnCP2FunctionBasis::GetRealFunctionValue(double theta1, double theta2, int index)
{
  int r = this->quantumNumberR[index];
  int s = this->quantumNumberS[index];
  int t = this->NbrFluxQuanta - this->quantumNumberR[index] - this->quantumNumberS[index];
  return this->Prefactor[index] * pow(cos (theta1), (double) (r)) * pow(sin(theta1) * cos(theta2), (double) (s)) * pow( sin(theta1) * sin(theta2), (double) (t));
}


void ParticleOnCP2FunctionBasis::GetRealFunctionValues(double theta1, double theta2, double *result)
{
  int r;
  int s;
  int t;
  for (int index=0; index <= this->LzMax; ++index)
  {
    r = this->quantumNumberR[index];
    s = this->quantumNumberS[index];
    t = this->NbrFluxQuanta - this->quantumNumberR[index] - this->quantumNumberS[index];
    result[index]=this->Prefactor[index] * pow(cos (theta1), (double) (r)) * pow(sin(theta1) * cos(theta2), (double) (s)) * pow( sin(theta1) * sin(theta2), (double) (t));
  }
}
