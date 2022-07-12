////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                                                                            //
//                        last modification : 10/01/2005                      //
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
#include "Tools/FQHEWaveFunction/JainCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/LandauSpectrumOnSphere.h"
#include "Matrix/ComplexMatrix.h"
#include "MathTools/FactorialCoefficient.h"
#include "Vector/RealVector.h"
#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <cstdlib>
#include <cstring>

using std::cout;
using std::endl;


// default constructor
//

JainCFOnSphereWaveFunction::JainCFOnSphereWaveFunction()
{
  this->NbrParticles = 0;
}

// constructor
//
// filename = name of the file describing the occupation of the pseudo-Landau levels

JainCFOnSphereWaveFunction::JainCFOnSphereWaveFunction(char* filename)
{
  this->NbrParticles = 0;
  ConfigurationParser CFDefinition;
  if ((CFDefinition.Parse(filename) == true) && (this->ParseGeneralInformation(CFDefinition) == true))
    {
      if (this->NbrParticles != 0)
	{
	  if ((this->ActualJastrowPower & 1) == 1)
	    {
	      this->JastrowPower = (this->ActualJastrowPower + 1) >> 1;
	    }
	  else
	    {
	      this->JastrowPower = this->ActualJastrowPower >> 1;
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
    }
  else
    {
      this->NbrParticles = 0;
    }
}

// copy constructor
//
// function = reference on the wave function to copy

JainCFOnSphereWaveFunction::JainCFOnSphereWaveFunction(const JainCFOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  if (this->NbrParticles != 0)
    {
      this->NbrLandauLevels = function.NbrLandauLevels;
      this->TwiceS = function.TwiceS;
      this->JastrowPower = function.JastrowPower;
      this->JastrowPowerPowers = function.JastrowPowerPowers;
      this->Flag = function.Flag;
      this->LevelOccupation = function.LevelOccupation;
      this->LinearCombinationCoefficients =  function.LinearCombinationCoefficients;
      this->NbrLinearCombination = function.NbrLinearCombination;
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
}

// destructor
//

JainCFOnSphereWaveFunction::~JainCFOnSphereWaveFunction()
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
	  delete[] this->JastrowPowerPowers;
	  delete[] this->LevelOccupation;
	  delete[] this->LinearCombinationCoefficients;
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
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* JainCFOnSphereWaveFunction::Clone ()
{
  return new JainCFOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex JainCFOnSphereWaveFunction::operator ()(RealVector& x)
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
  
  Complex JastrowFactor = this->EvaluateTables();
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
		      Slater.SetMatrixElement(Index, i, this->EvaluateCFMonopoleHarmonic(i, k, j, MaxMomentum));
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
Complex JainCFOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
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
  
  Complex JastrowFactor = this->EvaluateTables();
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
		      Slater.SetMatrixElement(Index, i, this->EvaluateCFMonopoleHarmonic(i, k, j, MaxMomentum));
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

// parse general informations about the composite fermion state
// 
// state = reference on the configuration parser that contains the informations
// return value = false if an error occured

bool JainCFOnSphereWaveFunction::ParseGeneralInformation(ConfigurationParser& state)
{
  if ((state["Geometry"] == 0) || (strncmp(state["Geometry"], "sphere", 6) != 0))
    {
      return false;
    }
  if (state["LandauLevel"] != 0)
    {
      this->NbrLandauLevels = atoi(state["LandauLevel"]);
    }
  else
    {
      return false;
    }
  if (state["JastrowFactor"] != 0)
    {
      this->ActualJastrowPower = atoi(state["JastrowFactor"]);
      if (this->ActualJastrowPower < 0)
	{
	  return false;
	}
      if (state["Statistics"] == 0)
	{
	  return false;
	}
      if (((strncmp (state["Statistics"], "bosons", 6)  == 0) && ((this->ActualJastrowPower & 1) == 0)) ||
	  ((strncmp (state["Statistics"], "fermions", 8)  == 0) && ((this->ActualJastrowPower & 1) == 1)))
	{
	  return false;
	}
    }
  else
    {
      return false;
    }
  if (state["Flux"] != 0)
    {
      this->TwiceS = atoi(state["Flux"]);
      if (this->TwiceS < 0)
	{
	  this->TwiceS = - this->TwiceS;
	  this->ReverseFluxFlag = true;
	}
      else
	{
	  this->ReverseFluxFlag = false;
	}
    }
  else
    {
      return false;
    }
  if (state["Levels"] == 0)
    {
      return false;
    }
  if ((state["NbrLinear"] != 0) && (state["Linear"] != 0))
    {
      this->NbrLinearCombination = atoi (state["NbrLinear"]);
      if (this->NbrLinearCombination <= 0)
	{
	  return false;
	}
      this->LinearCombinationCoefficients = new double [this->NbrLinearCombination];
      int TmpNbrLinearCombination = 0;
      char* Start = state["Linear"];
      char* End = Start;
      while ((*Start) != '\0')
	{
	  while (((*End) != '\0') && ((*End) != '|'))
	    {
	      ++End;
	    }
	  while ((Start != End) && ((((*End) < '0') || ((*End) > '9')) && ((*End) != '.')))
	    {
	      --End;
	    }
	  ++End;
	  if (Start == End)
	    {
	      delete[] this->LinearCombinationCoefficients;
	      return false;		  
	    }
	  this->LinearCombinationCoefficients[TmpNbrLinearCombination] = strtod(Start, &End);
	  if ((this->LinearCombinationCoefficients[TmpNbrLinearCombination] == 0.0) && (Start == End))
	    {
	      delete[] this->LinearCombinationCoefficients;
	      return false;		  		  
	    }
	  ++TmpNbrLinearCombination;
	  while (((*End) != '\0') && ((((*End) < '0') || ((*End) > '9')) && ((*End) != '.') && ((*End) != '+') && ((*End) != '-')))
	    {
	      ++End;
	    }
	  Start = End;
	}
      if (TmpNbrLinearCombination != this->NbrLinearCombination)
	{
	  delete[] this->LinearCombinationCoefficients;
	  return false;		  		  	  
	}
      this->LevelOccupation = new LandauSpectrumOnSphere [this->NbrLinearCombination];
      char* TmpLevelDescription = new char [strlen(state["Levels"]) + 1];
      strcpy (TmpLevelDescription, state["Levels"]);
      Start = TmpLevelDescription;
      while (((*Start) != '\0') && ((*Start) != '(') && (((*Start) == ' ') || ((*Start) == '\t')))
	++Start;
      if ((*Start) != '(')
	{
	  delete[] this->LinearCombinationCoefficients;
	  delete[] TmpLevelDescription;
	  return false;		  		  	  
	}
      TmpNbrLinearCombination = 0;
      ++Start;
      while ((TmpNbrLinearCombination < this->NbrLinearCombination) &&  ((*Start) != '\0'))
	{
	  End = Start;
	  while (((*End) != '\0') && ((*End) != ')'))
	    ++End;
	  if ((*End) != ')')
	    {
	      delete[] this->LinearCombinationCoefficients;
	      delete[] TmpLevelDescription;
	      delete[] this->LevelOccupation;
	      return false;		  		  	  	      
	    }
	  (*End) = '\0';
	  this->LevelOccupation[TmpNbrLinearCombination] = LandauSpectrumOnSphere(this->NbrLandauLevels, this->TwiceS, Start);
	  if (this->LevelOccupation[TmpNbrLinearCombination].GetNbrParticles() <= 0)
	    {
	      delete[] this->LinearCombinationCoefficients;
	      delete[] TmpLevelDescription;
	      delete[] this->LevelOccupation;
	      return false;		  		  	  	      	      
	    }
	  ++TmpNbrLinearCombination;
	  ++End;
	  while (((*End) != '\0') && ((*End) != '(') && (((*End) == ' ') || ((*End) == '\t') || ((*End) == '\n')))
	    ++End;
	  if (((*End) != '(') && (TmpNbrLinearCombination != this->NbrLinearCombination))
	    {
	      delete[] this->LinearCombinationCoefficients;
	      delete[] TmpLevelDescription;
	      delete[] this->LevelOccupation;
	      return false;		  		  	  	      	      
	    }
	  else
	    {
	      Start = End + 1;	  
	    }
	}
      if (TmpNbrLinearCombination != this->NbrLinearCombination)
	{
	  delete[] TmpLevelDescription;
	  delete[] this->LinearCombinationCoefficients;
	  delete[] this->LevelOccupation;
	  return false;		  		  	  
	}
      this->NbrParticles = this->LevelOccupation[0].GetNbrParticles();
      for (int i = 1; i < this->NbrLinearCombination; ++i)
	if (this->NbrParticles != this->LevelOccupation[i].GetNbrParticles())
	  {
	    delete[] this->LinearCombinationCoefficients;
	    delete[] TmpLevelDescription;
	    delete[] this->LevelOccupation;
	    return false;		  		  	  	      	      
	  }
    }
  else
    {
      this->LevelOccupation = new LandauSpectrumOnSphere [1];
      this->LevelOccupation[0] = LandauSpectrumOnSphere(this->NbrLandauLevels, this->TwiceS, state["Levels"]);      
      this->NbrParticles = this->LevelOccupation[0].GetNbrParticles();
      if (this->NbrParticles == 0)
	{
	  delete[] this->LevelOccupation;
	  return false;
	}
      this->LinearCombinationCoefficients = new double [1];
      this->LinearCombinationCoefficients[0] = 1.0;
      this->NbrLinearCombination = 1;
    }
  return true;
}


