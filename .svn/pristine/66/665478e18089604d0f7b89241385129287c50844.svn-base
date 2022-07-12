////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of unprojected Jain composite fermion wave function on sphere     //
//                                                                            //
//                        last modification : 13/04/2005                      //
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
#include "Tools/FQHEWaveFunction/UnprojectedJainCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/LandauSpectrumOnSphere.h"
#include "Matrix/ComplexMatrix.h"
#include "MathTools/FactorialCoefficient.h"
#include "Vector/RealVector.h"
#include "GeneralTools/ConfigurationParser.h"

#include <iostream>

using std::cout;
using std::endl;


// constructor
//
// filename = name of the file describing the occupation of the pseudo-Landau levels

UnprojectedJainCFOnSphereWaveFunction::UnprojectedJainCFOnSphereWaveFunction(char* filename)
{
  this->NbrParticles = 0;
  ConfigurationParser CFDefinition;
  if ((CFDefinition.Parse(filename) == true) && (this->ParseGeneralInformation(CFDefinition) == true))
    {
      if (this->NbrParticles != 0)
	{
	  this->JastrowPower = this->ActualJastrowPower;
	  this->Flag.Initialize();
	  this->EvaluateNormalizationPrefactors();  
	  this->EvaluateSumPrefactors();
	  this->JastrowFactorElements = new Complex*[this->NbrParticles];
	  for (int i = 1; i < this->NbrParticles; ++i)
	    this->JastrowFactorElements[i] = new Complex[i];
	  this->SpinorUCoordinates = new Complex[NbrParticles];
	  this->SpinorVCoordinates = new Complex[NbrParticles];
	  this->SpinorUCoordinatePower = new Complex*[NbrParticles];
	  this->SpinorVCoordinatePower = new Complex*[NbrParticles];
	  int MaxPower = this->TwiceS + 2 * (this->NbrLandauLevels - 1);
	  for (int i = 0; i < this->NbrParticles; ++i)
	    {
	      this->SpinorUCoordinatePower[i] = new Complex[MaxPower + 1];
	      this->SpinorVCoordinatePower[i] = new Complex[MaxPower + 1];
	    }
	}
    }
  else
    {
      this->NbrParticles = 0;
    }
}

// copy constructor
//
// function = reference on the wave function to copy

UnprojectedJainCFOnSphereWaveFunction::UnprojectedJainCFOnSphereWaveFunction(const UnprojectedJainCFOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  if (this->NbrParticles != 0)
    {
      this->NbrLandauLevels = function.NbrLandauLevels;
      this->TwiceS = function.TwiceS;
      this->JastrowPower = function.JastrowPower;
      this->Flag = function.Flag;
      this->ReverseFluxFlag = function.ReverseFluxFlag;
      this->ActualJastrowPower = function.ActualJastrowPower;
      this->LevelOccupation = function.LevelOccupation;
      this->LinearCombinationCoefficients =  function.LinearCombinationCoefficients;
      this->NbrLinearCombination = function.NbrLinearCombination;
      this->NormalizationPrefactors = function.NormalizationPrefactors;
      this->SumPrefactors = function.SumPrefactors;
      this->JastrowFactorElements = new Complex*[this->NbrParticles];
      for (int i = 1; i < this->NbrParticles; ++i)
	this->JastrowFactorElements[i] = new Complex[i];
      this->SpinorUCoordinates = new Complex[NbrParticles];
      this->SpinorVCoordinates = new Complex[NbrParticles];
      this->SpinorUCoordinatePower = new Complex*[NbrParticles];
      this->SpinorVCoordinatePower = new Complex*[NbrParticles];
      int MaxPower = (this->NbrParticles / this->NbrLandauLevels) + this->NbrLandauLevels - 2;
      for (int i = 0; i < this->NbrParticles; ++i)
	{
	  this->SpinorUCoordinatePower[i] = new Complex[MaxPower + 1];
	  this->SpinorVCoordinatePower[i] = new Complex[MaxPower + 1];
	}
    }
}

// destructor
//

UnprojectedJainCFOnSphereWaveFunction::~UnprojectedJainCFOnSphereWaveFunction()
{
  if (this->NbrParticles != 0)
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
	  delete[] this->LevelOccupation;
	  delete[] this->LinearCombinationCoefficients;
	}
      for (int i = 0; i < this->NbrParticles; ++i)
	delete[] this->JastrowFactorElements[i];
      delete[] this->JastrowFactorElements;
      for (int i = 0; i < this->NbrParticles; ++i)
	{
	  delete[] this->SpinorUCoordinatePower[i];
	  delete[] this->SpinorVCoordinatePower[i];
	}
      delete[] this->SpinorUCoordinates;
      delete[] this->SpinorVCoordinates;
      delete[] this->SpinorUCoordinatePower;
      delete[] this->SpinorVCoordinatePower;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* UnprojectedJainCFOnSphereWaveFunction::Clone ()
{
  return new UnprojectedJainCFOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex UnprojectedJainCFOnSphereWaveFunction::operator ()(RealVector& x)
{
  if (this->NbrParticles == 0)
    {
      return Complex(0.0);
    }
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
  
  Complex JastrowFactor = this->EvaluateTables(false);
  Complex TmpValue(0.0);
  for (int n = 0; n < this->NbrLinearCombination; ++n)
    {
      for (int i = 0; i < this->NbrParticles; ++i)
	{
	  int Index = 0;
	  int MaxMomentum = this->TwiceS;
	  for (int j = 0; j < this->NbrLandauLevels; ++j)
	    {
	      for (int k = 0; k <= MaxMomentum; ++k)
		{
		  if ((this->LevelOccupation[n])(j, k) == 1)
		    {
		      Slater.SetMatrixElement(Index, i, this->EvaluateMonopoleHarmonic(i, k, j, MaxMomentum));
		      ++Index;
		    }
		}
	      MaxMomentum += 2;
	    }
	} 
      TmpValue += this->LinearCombinationCoefficients[n] * Slater.Determinant();
    }

  return JastrowFactor * TmpValue;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex UnprojectedJainCFOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  if (this->NbrParticles == 0)
    {
      return Complex(0.0);
    }
  ComplexMatrix Slater (this->NbrParticles, this->NbrParticles);

  // Import from spinors
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = uv.Re(2*i);
      this->SpinorUCoordinates[i].Im = uv.Im(2*i);
      this->SpinorVCoordinates[i].Re = uv.Re(2*i+1);
      this->SpinorVCoordinates[i].Im = uv.Im(2*i+1);
    }
  
  Complex JastrowFactor = this->EvaluateTables(false);
  Complex TmpValue(0.0);
  for (int n = 0; n < this->NbrLinearCombination; ++n)
    {
      for (int i = 0; i < this->NbrParticles; ++i)
	{
	  int Index = 0;
	  int MaxMomentum = this->TwiceS;
	  for (int j = 0; j < this->NbrLandauLevels; ++j)
	    {
	      for (int k = 0; k <= MaxMomentum; ++k)
		{
		  if ((this->LevelOccupation[n])(j, k) == 1)
		    {
		      Slater.SetMatrixElement(Index, i, this->EvaluateMonopoleHarmonic(i, k, j, MaxMomentum));
		      ++Index;
		    }
		}
	      MaxMomentum += 2;
	    }
	} 
      TmpValue += this->LinearCombinationCoefficients[n] * Slater.Determinant();
    }

  return JastrowFactor * TmpValue;
}  

// evaluate monopole spherical harmonic 
//
// coordinate = index of the coordinate
// momentum = monopole spherical harmonic Lz momentum (plus S shift)
// landauLevel = index of the Landau level
// maximumMomentum = maxixum momentum that can be reached in the Landau level
// return value = value of the monopole spherical harmonic at the given point

Complex UnprojectedJainCFOnSphereWaveFunction::EvaluateMonopoleHarmonic (int coordinate, int momentum, int landauLevel, int maximumMomentum)
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
  if (this->ReverseFluxFlag == false)
    {
      for (; i <= Max; ++i)
	{
	  Tmp += (this->SumPrefactors[landauLevel][i][momentum] * this->SpinorUCoordinatePower[coordinate][UPower] * 
		  this->SpinorVCoordinatePower[coordinate][VPower] * Conj(this->SpinorUCoordinatePower[coordinate][i]) * 
		  Conj(this->SpinorVCoordinatePower[coordinate][landauLevel - i]));
	  ++UPower;
	  --VPower;
	}
    }
  else
    {
      for (; i <= Max; ++i)
	{
	  Tmp += (this->SumPrefactors[landauLevel][i][momentum] * Conj(this->SpinorUCoordinatePower[coordinate][UPower]) * 
		  Conj(this->SpinorVCoordinatePower[coordinate][VPower]) * this->SpinorUCoordinatePower[coordinate][i] * 
		  this->SpinorVCoordinatePower[coordinate][landauLevel - i]);
	  ++UPower;
	  --VPower;
	}
    }
  return (this->NormalizationPrefactors[landauLevel][momentum] * Tmp);
}


// evaluate constant factors that appears in the sum of projected monopole harmonic (except LLL)
//

void UnprojectedJainCFOnSphereWaveFunction::EvaluateSumPrefactors()
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
