////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of function basis for particle on sphere             //
//                            in a given Landau level                         //
//                                                                            //
//                        last modification : 29/11/2005                      //
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
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"

#include <math.h>


// constructor
//
// lzMax = twice the maximum Lz value reached by a particle in the lowest Landau level (i.e. number of flux quanta)
// landauLevel = index of the Landau level to consider (0 for the lowest Landau level)

ParticleOnSphereGenericLLFunctionBasis::ParticleOnSphereGenericLLFunctionBasis(int lzMax, int landauLevel)
{
  this->LzMax = lzMax;
  this->LandauLevel = landauLevel;
  this->HilbertSpaceDimension = this->LzMax + (2 * this->LandauLevel) + 1;
  this->EvaluateSumPrefactors();
  this->EvaluateNormalizationPrefactors();
}

// destructor
//

ParticleOnSphereGenericLLFunctionBasis::~ParticleOnSphereGenericLLFunctionBasis ()
{
  delete[] this->NormalizationPrefactors;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    delete[] this->SumPrefactors[i];
  delete[] this->SumPrefactors;
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// value = reference on the value where the function has to be evaluated
// result = reference on the value where the result has to be stored
// index = function index 

void ParticleOnSphereGenericLLFunctionBasis::GetFunctionValue(RealVector& value, Complex& result, int index)
{
  double Cos = cos (0.5 * value[0]);
  double Sin = sin (0.5 * value[0]);
  int Max = this->HilbertSpaceDimension - 1;
  result.Re = value[1] * (((double) index) - 0.5 * ((double) Max));
  result.Im = sin(result.Re);
  result.Re = cos(result.Re);
  double* TmpCoefficients = this->SumPrefactors[index];
  int Min = this->LandauLevel - index;
  if (Min < 0)
    Min = 0;
  Max -= index;
  double DMin = (double) (Max + this->LandauLevel - (2 * Min));  
  if (Max > this->LandauLevel)
    {
      Max = this->LandauLevel;
    }
  double Tmp = 0.0;
  double TmpCos = pow (Cos, (double) (index + (2 * Min) - this->LandauLevel));
  Cos *= Cos;
  Max -= Min;
  for (Min = 0; Min <= Max; ++Min)
    {
      Tmp += TmpCoefficients[Min] * TmpCos * pow(Sin, DMin);
      TmpCos *= Cos;
      DMin -= 2.0;
    }  
  result *= (this->NormalizationPrefactors[index] * Tmp);
}



// evaluate normalization factors of monopole harmonics
//

void ParticleOnSphereGenericLLFunctionBasis::EvaluateNormalizationPrefactors()
{
  int MaxMomentum = this->HilbertSpaceDimension - 1;
  FactorialCoefficient Coef;
  this->NormalizationPrefactors = new double[MaxMomentum + 1];
  double Factor = ((double) this->HilbertSpaceDimension) / (4.0  * M_PI);
  if (this->LandauLevel == 0)
    {
      for (int j = 0; j <= MaxMomentum; ++j)
	{
	  Coef.SetToOne();
	  Coef.PartialFactorialMultiply(this->HilbertSpaceDimension - j, MaxMomentum);
	  Coef.FactorialDivide(j);
	  this->NormalizationPrefactors[j] = sqrt(Factor * Coef.GetNumericalValue());
	  if ((j & 1) != 0)
	    {
	      this->NormalizationPrefactors[j] *= -1.0;
	    }
	}
    }
  else
    {
      double Sign = 1.0;
      if ((this->LzMax & 1) != 0)
	Sign = -1.0;
      for (int j = 0; j <= MaxMomentum; ++j)
	{	  
	  Coef.SetToOne();
	  Coef.FactorialMultiply(MaxMomentum - j);
	  Coef.FactorialDivide(this->LandauLevel);
	  Coef.FactorialMultiply(j);
	  Coef.FactorialDivide(MaxMomentum - this->LandauLevel);
	  this->NormalizationPrefactors[j] = sqrt(Factor * Coef.GetNumericalValue());
	  this->NormalizationPrefactors[j] *= Sign;
	  Sign *= -1.0;
	}
    }
}

// evaluate constant factors that appears in the sum of projected monopole harmonic
//

void ParticleOnSphereGenericLLFunctionBasis::EvaluateSumPrefactors()
{
  this->SumPrefactors = new double* [this->HilbertSpaceDimension];
  if (this->LandauLevel == 0)
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  this->SumPrefactors[i] = new double [1];
	  this->SumPrefactors[i][0] = 1.0;
	}
    }
  else
    {
      FactorialCoefficient Coef;
      int MaxMomentum = this->HilbertSpaceDimension - 1;      
      for (int j = 0; j <= MaxMomentum; ++j)
	{
	  double Factor = 1.0;
	  int Min = this->LandauLevel - j;
	  if (Min < 0)
	    Min = 0;
	  else
	    if ((Min & 1) != 0)
	      Factor *= -1.0;
	  int Max = MaxMomentum - j;
	  if (Max > this->LandauLevel)
	    Max = this->LandauLevel;
	  this->SumPrefactors[j] = new double [Max - Min + 1];
	  for (int k = Min; k <= Max; ++k)  
	    {
	      Coef.SetToOne();
	      Coef.PartialFactorialMultiply(k + 1, this->LandauLevel);
	      Coef.FactorialDivide(this->LandauLevel - k);
	      Coef.PartialFactorialMultiply(j + k + 1 - this->LandauLevel, this->LzMax + this->LandauLevel);	  
	      Coef.FactorialDivide(MaxMomentum - j - k);	  
	      this->SumPrefactors[j][k - Min] = Factor * Coef.GetNumericalValue();
	      Factor *= -1.0;	      
	    }
	}
    }
}


