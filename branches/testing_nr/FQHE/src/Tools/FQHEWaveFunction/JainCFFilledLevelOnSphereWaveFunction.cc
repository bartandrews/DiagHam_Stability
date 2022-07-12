////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 16/09/2004                      //
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
#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"
#include "Matrix/ComplexMatrix.h"
#include "MathTools/FactorialCoefficient.h"
#include "Vector/RealVector.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
//

JainCFFilledLevelOnSphereWaveFunction::JainCFFilledLevelOnSphereWaveFunction()
{
  this->NbrParticles = 0;
}

// constructor
//
// nbrParticles = number of particles
// nbrLandauLevel = number of Landau levels filled with composite fermions
// jastrowPower = power to which the Jastrow factor has to be raised

JainCFFilledLevelOnSphereWaveFunction::JainCFFilledLevelOnSphereWaveFunction(int nbrParticles, int nbrLandauLevels, int jastrowPower)
{
  this->NbrParticles = nbrParticles;
  this->NbrLandauLevels = nbrLandauLevels;
  this->TwiceS = (this->NbrParticles / this->NbrLandauLevels) - this->NbrLandauLevels;
  this->ActualJastrowPower = jastrowPower;
  if ((this->ActualJastrowPower & 1) == 1)
    {
      this->JastrowPower = (jastrowPower + 1) >> 1;
    }
  else
    {
      this->JastrowPower = jastrowPower >> 1;
    }
  this->JastrowPowerPowers = new double [this->NbrLandauLevels];
  this->JastrowPowerPowers[0] = 1.0;
  for (int i = 1; i < this->NbrLandauLevels; ++i)
    this->JastrowPowerPowers[i] = this->JastrowPowerPowers[i - 1] * ((double) this->JastrowPower);
  this->Flag.Initialize();
  this->EvaluateNormalizationPrefactors();  
  this->EvaluateSumPrefactors();
  this->JastrowFactorElements = new Complex*[this->NbrParticles];
  for (int i = 1; i < this->NbrParticles; ++i)
    this->JastrowFactorElements[i] = new Complex[i];
  this->DerivativeFactors = new Complex** [this->NbrParticles];
  this->DerivativeFactors2 = new Complex** [this->NbrParticles];
  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];
  this->SpinorUCoordinatePower = new Complex*[NbrParticles];
  this->SpinorVCoordinatePower = new Complex*[NbrParticles];
  int MaxPower = this->TwiceS + 2 * (this->NbrLandauLevels - 1);
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinatePower[i] = new Complex[MaxPower + 1];
      this->SpinorVCoordinatePower[i] = new Complex[MaxPower + 1];
      this->DerivativeFactors[i] = new Complex* [this->NbrLandauLevels];
      this->DerivativeFactors2[i] = new Complex* [this->NbrLandauLevels];
      for (int j = 0; j < this->NbrLandauLevels; ++j)
	{
	  this->DerivativeFactors[i][j] = new Complex [this->NbrLandauLevels];
	  this->DerivativeFactors2[i][j] = new Complex [this->NbrLandauLevels];
	}
    }
}

// copy constructor
//
// function = reference on the wave function to copy

JainCFFilledLevelOnSphereWaveFunction::JainCFFilledLevelOnSphereWaveFunction(const JainCFFilledLevelOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrLandauLevels = function.NbrLandauLevels;
  this->TwiceS = function.TwiceS;
  this->JastrowPower = function.JastrowPower;
  this->JastrowPowerPowers = function.JastrowPowerPowers;
  this->Flag = function.Flag;
  this->NormalizationPrefactors = function.NormalizationPrefactors;
  this->SumPrefactors = function.SumPrefactors;
  this->JastrowFactorElements = new Complex*[this->NbrParticles - 1];
  for (int i = 1; i < this->NbrParticles; ++i)
    this->JastrowFactorElements[i] = new Complex[i];
  this->DerivativeFactors = new Complex** [this->NbrParticles];
  this->DerivativeFactors2 = new Complex** [this->NbrParticles];
  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];
  this->SpinorUCoordinatePower = new Complex*[NbrParticles];
  this->SpinorVCoordinatePower = new Complex*[NbrParticles];
  int MaxPower = (this->NbrParticles / this->NbrLandauLevels) + this->NbrLandauLevels - 2;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinatePower[i] = new Complex[MaxPower + 1];
      this->SpinorVCoordinatePower[i] = new Complex[MaxPower + 1];
      this->DerivativeFactors[i] = new Complex* [this->NbrLandauLevels];
      this->DerivativeFactors2[i] = new Complex* [this->NbrLandauLevels];
      for (int j = 0; j < this->NbrLandauLevels; ++j)
	{
	  this->DerivativeFactors[i][j] = new Complex [this->NbrLandauLevels];
	  this->DerivativeFactors2[i][j] = new Complex [this->NbrLandauLevels];
	}
    }
}

// destructor
//

JainCFFilledLevelOnSphereWaveFunction::~JainCFFilledLevelOnSphereWaveFunction()
{
  if (this->Flag.Shared() == false)
    {
      delete[] this->NormalizationPrefactors;
      for (int i = 0; i < this->NbrLandauLevels; ++i)
	{
	  for (int j = 0; j < this->NbrLandauLevels; ++j)
	    delete[] this->SumPrefactors[i][j];
	  delete[] this->SumPrefactors[i];
	}
      delete[] this->SumPrefactors;
      delete[] this->JastrowPowerPowers;
    }
  for (int i = 0; i < this->NbrParticles; ++i)
    delete[] this->JastrowFactorElements[i];
  delete[] this->JastrowFactorElements;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      for (int j = 0; j < this->NbrLandauLevels; ++j)
	{
	  delete[] this->DerivativeFactors[i][j];
	  delete[] this->DerivativeFactors2[i][j];
	}
      delete[] this->DerivativeFactors[i];
      delete[] this->DerivativeFactors2[i];
      delete[] this->SpinorUCoordinatePower[i];
      delete[] this->SpinorVCoordinatePower[i];
    }
  delete[] this->SpinorUCoordinates;
  delete[] this->SpinorVCoordinates;
  delete[] this->SpinorUCoordinatePower;
  delete[] this->SpinorVCoordinatePower;
  delete[] this->DerivativeFactors;
  delete[] this->DerivativeFactors2;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* JainCFFilledLevelOnSphereWaveFunction::Clone ()
{
  return new JainCFFilledLevelOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex JainCFFilledLevelOnSphereWaveFunction::operator ()(RealVector& x)
{
  ComplexMatrix Slater (this->NbrParticles, this->NbrParticles);

  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = cos(0.5 * x[i << 1]);
      this->SpinorUCoordinates[i].Im = this->SpinorUCoordinates[i].Re;
      this->SpinorUCoordinates[i].Re *= cos(0.5 * x[1 + (i << 1)]);
      this->SpinorUCoordinates[i].Im *= sin(0.5 * x[1 + (i << 1)]);
      this->SpinorVCoordinates[i].Re = sin(0.5 * x[i << 1]);
      this->SpinorVCoordinates[i].Im = this->SpinorVCoordinates[i].Re;
      this->SpinorVCoordinates[i].Re *= cos(0.5 * x[1 + (i << 1)]);
      this->SpinorVCoordinates[i].Im *= -sin(0.5 * x[1 + (i << 1)]);
    }
  
  Complex JastrowFactor = this->EvaluateTables();

  for (int i = 0; i < this->NbrParticles; ++i)
    {
      int Index = 0;
      int MaxMomentum = this->TwiceS;
      for (int j = 0; j < this->NbrLandauLevels; ++j)
	{
	  for (int k = 0; k <= MaxMomentum; ++k)
	    {
	      Slater.SetMatrixElement(Index, i, EvaluateCFMonopoleHarmonic(i, k, j, MaxMomentum));
	      ++Index;
	    }
	  MaxMomentum += 2;
	}
    }  

  return JastrowFactor * Slater.Determinant();
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex JainCFFilledLevelOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  ComplexMatrix Slater (this->NbrParticles, this->NbrParticles);
  
  // Import from spinors
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = uv.Re(2*i);
      this->SpinorUCoordinates[i].Im = uv.Im(2*i);
      this->SpinorVCoordinates[i].Re = uv.Re(2*i+1);
      this->SpinorVCoordinates[i].Im = uv.Im(2*i+1);
    }

  Complex JastrowFactor = this->EvaluateTables();

  for (int i = 0; i < this->NbrParticles; ++i)
    {
      int Index = 0;
      int MaxMomentum = this->TwiceS;
      for (int j = 0; j < this->NbrLandauLevels; ++j)
	{
	  for (int k = 0; k <= MaxMomentum; ++k)
	    {
	      Slater.SetMatrixElement(Index, i, EvaluateCFMonopoleHarmonic(i, k, j, MaxMomentum));
	      ++Index;
	    }
	  MaxMomentum += 2;
	}
    }  

  return JastrowFactor * Slater.Determinant();
  
}

// evaluate precalculation tables used during wave function evaluation (called at each evaluation)
//
// requires SpinorUCoordinates and SpinorVCoordinates to be initialized prior to call
// derivativeFlag = indicate if precalculation tables invloved in derivative evaluation have to be calculated
// return value = value of the Jastrow factor

Complex JainCFFilledLevelOnSphereWaveFunction::EvaluateTables(bool derivativeFlag)
{  

  int MaxPower = this->TwiceS + 2 * (this->NbrLandauLevels - 1);
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinatePower[i][0] = 1.0;
      this->SpinorVCoordinatePower[i][0] = 1.0;
      Complex TmpU = this->SpinorUCoordinates[i];
      Complex TmpV = this->SpinorVCoordinates[i];
      for (int j = 1; j <= MaxPower; ++j)
	{
	  this->SpinorUCoordinatePower[i][j] = this->SpinorUCoordinatePower[i][j - 1] * TmpU;
	  this->SpinorVCoordinatePower[i][j] = this->SpinorVCoordinatePower[i][j - 1] * TmpV;
	}
    }

  Complex JastrowFactor(1.0);
  Complex Tmp;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      for (int j = 0; j < i; ++j)
	{
	  Tmp = ((this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]));
	  this->JastrowFactorElements[i][j] = 1.0 / Tmp;
	  JastrowFactor *= Tmp;
	}
    }
  Tmp = JastrowFactor;
  for (int i = 1; i < this->ActualJastrowPower; ++i)
    {
      JastrowFactor *= Tmp;
    }

  if (derivativeFlag == true)
    {
      Complex Tmp2;
      for (int i = 0; i < this->NbrParticles; ++i)
	{     
	  Complex** TmpDerivativeFactors = this->DerivativeFactors[i];
	  Complex** TmpDerivativeFactors2 = this->DerivativeFactors2[i];
	  for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
	    for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
	      TmpDerivativeFactors[k1][k2] = 0.0;
	  
	  int Index = 0;
	  for (int j = 1; j < this->NbrParticles; ++j)
	    {
	      if (Index == i)
		++Index;
	      if (Index > i)
		Tmp2 = this->SpinorVCoordinates[Index] * this->JastrowFactorElements[Index][i];
	      else
		Tmp2 = -this->SpinorVCoordinates[Index] * this->JastrowFactorElements[i][Index];
	      Tmp = 1.0;
	      for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
		{
		  for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
		    {
		      TmpDerivativeFactors2[k1][k2] = Tmp;
		    }
		  Tmp *= Tmp2;
		}
	      if (Index > i)
		Tmp2 = - this->SpinorUCoordinates[Index] * this->JastrowFactorElements[Index][i];
	      else
		Tmp2 = this->SpinorUCoordinates[Index] * this->JastrowFactorElements[i][Index];
	      Tmp = 1.0;
	      for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
		{
		  for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
		    {
		      TmpDerivativeFactors2[k2][k1] *= Tmp;
		    }
		  Tmp *= Tmp2;
		}
	      for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
		for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
		  TmpDerivativeFactors[k1][k2] += TmpDerivativeFactors2[k1][k2];
	      ++Index;
	    }
	}
    }
  return JastrowFactor;
}


// evaluate normalization factors of projected monopole harmonics
//

void JainCFFilledLevelOnSphereWaveFunction::EvaluateNormalizationPrefactors()
{
  this->NormalizationPrefactors = new double* [this->NbrLandauLevels];
  int MaxMomentum = this->TwiceS;
  FactorialCoefficient Coef;
  double Factor = ((double) (this->TwiceS + 1)) / (4.0  * M_PI);
  this->NormalizationPrefactors[0] = new double[MaxMomentum + 1];
  for (int j = 0; j <= MaxMomentum; ++j)
    {
      Coef.SetToOne();
      Coef.PartialFactorialDivide(this->TwiceS - j + 1, this->TwiceS);
      Coef.FactorialMultiply(j);
      this->NormalizationPrefactors[0][j] = sqrt(Factor * Coef.GetNumericalValue());
      if ((j & 1) != 0)
	{
	  this->NormalizationPrefactors[0][j] *= -1.0;
	}
    }
  MaxMomentum += 2;
  for (int i = 1; i < this->NbrLandauLevels; ++i)  
    {
      double Factor = ((double) (this->TwiceS + (2 * i) + 1)) / (4.0  * M_PI);
      this->NormalizationPrefactors[i] = new double[MaxMomentum + 1];
      double Sign = 1.0;
      if ((this->TwiceS & 1) != 0)
	Sign = -1.0;
      for (int j = 0; j <= MaxMomentum; ++j)
	{	  
	  Coef.SetToOne();
	  Coef.FactorialMultiply(this->TwiceS + 2 * i - j);
	  Coef.FactorialDivide(i);
	  Coef.FactorialMultiply(j);
	  Coef.FactorialDivide(this->TwiceS + i);
	  this->NormalizationPrefactors[i][j] = sqrt(Factor * Coef.GetNumericalValue());
	  this->NormalizationPrefactors[i][j] *= Sign;
	  Sign *= -1.0;
	}
      MaxMomentum += 2;
    }
}

// evaluate constant factors that appears in the sum of projected monopole harmonic (except LLL)
//

void JainCFFilledLevelOnSphereWaveFunction::EvaluateSumPrefactors()
{
  this->SumPrefactors = new double** [this->NbrLandauLevels];
  int MaxMomentum = this->TwiceS;
  FactorialCoefficient Coef;
  double Factor = 1.0;
  for (int i = 0; i < this->NbrLandauLevels; ++i)  
    {
      this->SumPrefactors[i] = new double* [i + 1];
      Factor = 1.0;
      for (int k = 0; k <= i; ++k)  
	{
	  this->SumPrefactors[i][k] = new double [MaxMomentum + 1];
	  for (int j = 0; j < (i - k); ++j)
	    this->SumPrefactors[i][k][j] = 0.0;
	  for (int j = i - k; j <= (MaxMomentum - k); ++j)
	    {
	      Coef.SetToOne();
	      Coef.PartialFactorialDivide(this->TwiceS + 2 * (1 + this->JastrowPower * (this->NbrParticles - 1)), 
					   this->TwiceS + 2 * (this->JastrowPower * (this->NbrParticles - 1)) + i + 1);
	      Coef.PartialFactorialMultiply(k + 1, i);
	      Coef.FactorialDivide(i - k);
	      Coef.PartialFactorialMultiply(j + k - i + 1, this->TwiceS + i);	  
	      Coef.FactorialDivide(this->TwiceS + (2 * i) -j - k);	  
	      this->SumPrefactors[i][k][j] = Factor * Coef.GetNumericalValue();
	    }
	  for (int j = MaxMomentum - k + 1; j <= MaxMomentum; ++j)
	    this->SumPrefactors[i][k][j] = 0.0;
	  Factor *= -1.0;
	}
      MaxMomentum += 2;
    }
}


// evaluate composite fermion monopole spherical harmonic 
//
// coordinate = index of the main coordinate (aka coordinate before project onto the lowest Landau level)
// momentum = monopole spherical harmonic Lz momentum (plus S shift)
// landauLevel = index of the pseudo Landau level
// maximumMomentum = maxixum momentum that can be reached in the current pseudo Landau level
// return value = value of the monopole spherical harmonic at the given point

Complex JainCFFilledLevelOnSphereWaveFunction::EvaluateCFMonopoleHarmonic (int coordinate, int momentum, int landauLevel, int maximumMomentum)
{
  Complex Tmp(0.0);
  int i = landauLevel - momentum;
  if (i < 0)
    i = 0;
  int UPower = momentum + i - landauLevel;
  int Max = maximumMomentum - momentum;
  int VPower = Max - i;
  if (Max > landauLevel)
    Max = landauLevel;
  for (; i <= Max; ++i)
    {
      Tmp += (this->SumPrefactors[landauLevel][i][momentum] * this->SpinorUCoordinatePower[coordinate][UPower] * 
	      this->SpinorVCoordinatePower[coordinate][VPower] * this->EvaluateCFMonopoleHarmonicDerivative(coordinate, i, landauLevel - i));
      ++UPower;
      --VPower;
    }
  return (this->NormalizationPrefactors[landauLevel][momentum] * Tmp);
}


// evaluate derivative part of the composite fermion monopole spherical harmonic 
//
// index = particle index
// alpha = number of (d\du) derivates
// beta = number of (d\dv) derivates
// return value = derivative contribution to the monopole spherical harmonic

Complex JainCFFilledLevelOnSphereWaveFunction::EvaluateCFMonopoleHarmonicDerivative(int index, int alpha, int beta)
{
  Complex Tmp;
  switch (alpha)
    {
    case 0:
      {
	switch (beta)
	  {
	  case 0:
	    Tmp = 1.0;
	    break;
	  case 1:
	    Tmp = this->JastrowPowerPowers[1] * this->DerivativeFactors[index][0][1];
	    break;
	  case 2:
	    Tmp = (this->JastrowPowerPowers[2] * this->DerivativeFactors[index][0][1] * this->DerivativeFactors[index][0][1] 
		   - this->JastrowPowerPowers[1] * this->DerivativeFactors[index][0][2]);
	    break;
	  case 3:
	    Tmp = (this->JastrowPowerPowers[3] * this->DerivativeFactors[index][0][1] * this->DerivativeFactors[index][0][1] * this->DerivativeFactors[index][0][1]
		   - 3 * this->JastrowPowerPowers[2] * this->DerivativeFactors[index][0][2] * this->DerivativeFactors[index][0][1]
		   + 2 * this->JastrowPowerPowers[1] * this->DerivativeFactors[index][0][3]);
	    break;
	  }
      }
      break;
    case 1:
      {
	switch (beta)
	  {
	  case 0:
	    Tmp = this->JastrowPowerPowers[1] * this->DerivativeFactors[index][1][0];
	    break;
	  case 1:
	    Tmp = (this->JastrowPowerPowers[2] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][0][1] 
		   - this->JastrowPowerPowers[1] * this->DerivativeFactors[index][1][1]);
	    break;
	  case 2:
	    Tmp = (this->JastrowPowerPowers[3] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][0][1] * this->DerivativeFactors[index][0][1]
		   - this->JastrowPowerPowers[2] * (2 * this->DerivativeFactors[index][1][1] * this->DerivativeFactors[index][0][1]
						    + this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][0][2])
		   - 2 * this->JastrowPowerPowers[1] * this->DerivativeFactors[index][1][2]);
	    break;
	  }
      }
      break;
    case 2:
      {
	switch (beta)
	  {
	  case 0:
	    Tmp = (this->JastrowPowerPowers[2] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][1][0] 
		   - this->JastrowPowerPowers[1] * this->DerivativeFactors[index][2][0]);
	    break;
	  case 1:
	    Tmp = (this->JastrowPowerPowers[3] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][0][1]
		   - this->JastrowPowerPowers[2] * (2 * this->DerivativeFactors[index][1][1] * this->DerivativeFactors[index][1][0]
						    + this->DerivativeFactors[index][0][1] * this->DerivativeFactors[index][2][0])
		   - 2 * this->JastrowPowerPowers[1] * this->DerivativeFactors[index][2][1]);
	    break;
	  }
      }
      break;
    case 3:
      {
	switch (beta)
	  {
	  case 0:
	    Tmp = (this->JastrowPowerPowers[3] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][1][0] * this->DerivativeFactors[index][1][0]
		   - 3 * this->JastrowPowerPowers[2] * this->DerivativeFactors[index][2][0] * this->DerivativeFactors[index][1][0]
		   + 2 * this->JastrowPowerPowers[1] * this->DerivativeFactors[index][3][0]);
	    break;
	  }
      }
      break;	    
    }
  return Tmp;
}
