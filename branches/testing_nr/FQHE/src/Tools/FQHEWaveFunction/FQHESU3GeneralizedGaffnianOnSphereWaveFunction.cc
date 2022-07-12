////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of SU(3) generalized Halperine wave function on sphere       //
//                        times the generalized permament                     //
//                                                                            //
//                        last modification : 05/04/2008                      //
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
#include "Tools/FQHEWaveFunction/FQHESU3GeneralizedGaffnianOnSphereWaveFunction.h"
#include "GeneralTools/ConfigurationParser.h"
#include "Vector/RealVector.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>

using std::cout;
using std::endl;


// constructor
//
// nbrParticles = total number of particles 
// mIntra =  coefficient of the intra-component correlations
// mInter = coefficient of the inter-component correlations

FQHESU3GeneralizedGaffnianOnSphereWaveFunction::FQHESU3GeneralizedGaffnianOnSphereWaveFunction(int nbrParticles, int mIntra, int mInter)
{
  this->NbrParticles = nbrParticles;
  this->MIntra = mIntra;
  this->MInter = mInter;
  this->SUKWaveFunction = 0;
  this->KValue = 3;
  this->NbrParticlesPerColor = this->NbrParticles / this->KValue;
  this->FermionFlag = false;
  this->FullySymmetrize = true;  
  if ((this->MIntra & 1) == 1)
    {
      this->FermionFlag= true;
      --this->MIntra;
      --this->MInter;
    }
  this->EvaluatePermutations();
  this->Flag.Initialize();
  this->TemporaryUV = ComplexVector(this->NbrParticles * 2);
}

// copy constructor
//
// function = reference on the wave function to copy

FQHESU3GeneralizedGaffnianOnSphereWaveFunction::FQHESU3GeneralizedGaffnianOnSphereWaveFunction(const FQHESU3GeneralizedGaffnianOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->MIntra = function.MIntra;
  this->MInter = function.MInter;
  this->SUKWaveFunction = 0;
  this->KValue = 3;
  this->NbrParticlesPerColor = this->NbrParticles / this->KValue;
  this->FermionFlag = function.FermionFlag;
  this->TemporaryUV = ComplexVector(this->NbrParticles * 2);
  this->Flag = function.Flag;
  this->FullySymmetrize = true;  
}

// constructor from configuration file
//
// configuration = reference on the configuration file parser
// errorFlag = reference on the error flag that is set to false if an error occured while retriving the configuration
// nbrParticles = reference on the total number of particles (computed from the configuration file datas)
// lzMax = twice the maximum angular momentum for a particle (computed from the configuration file datas)

FQHESU3GeneralizedGaffnianOnSphereWaveFunction::FQHESU3GeneralizedGaffnianOnSphereWaveFunction(ConfigurationParser& configuration, bool& errorFlag, int& nbrParticles, int& lzMax)
{
  errorFlag = true;
  if ((configuration["WaveFunction"] == 0) || (strcmp ("SU3GeneralizedGaffnian", configuration["WaveFunction"]) != 0))
    {
      errorFlag = false;
      return;
    }
  if ((configuration.GetAsSingleInteger("NbrParticles", this->NbrParticles) == false))
    {
      errorFlag = false;
      return;
    }
  if ((configuration.GetAsSingleInteger("MIntra", this->MIntra) == false) ||
      (configuration.GetAsSingleInteger("MInter", this->MInter) == false) ||
      (this->MIntra < 0) || (this->MInter < 0))
    {
      errorFlag = false;
      return;
    }
  this->SUKWaveFunction = 0;
  this->KValue = 3;
  this->NbrParticlesPerColor = this->NbrParticles / this->KValue;
  this->TemporaryUV = ComplexVector(this->NbrParticles * 2);
  this->FermionFlag = false;
  this->FullySymmetrize = true;  
  if ((this->MIntra & 1) == 1)
    {
      this->FermionFlag= true;
      --this->MIntra;
      --this->MInter;
    }
  this->EvaluatePermutations();
  this->Flag.Initialize();
}

// destructor
//

FQHESU3GeneralizedGaffnianOnSphereWaveFunction::~FQHESU3GeneralizedGaffnianOnSphereWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* FQHESU3GeneralizedGaffnianOnSphereWaveFunction::Clone ()
{
  return new FQHESU3GeneralizedGaffnianOnSphereWaveFunction(*this);
}

// evaluate function at a given point(the first 2*N1 coordinates correspond to the position of the type 1 particles, 
//                                     the following 2*N2 coordinates correspond to the position of the type 2 particles,
//                                     last the 2*N3 coordinates correspond to the position of the type 3 particles)
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
//      this method is for only for internal class usage
// return value = function value at (uv)

Complex FQHESU3GeneralizedGaffnianOnSphereWaveFunction::LocalCalculateFromSpinorVariables(ComplexVector& uv)
{
  Complex Tmp;
  Complex TmpU;
  Complex TmpV;
  int NbrN1N2 = 2 * this->NbrParticlesPerColor;
  Complex TotalWaveFunction(1.0);
  if (this->MIntra > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrParticlesPerColor; ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->NbrParticlesPerColor; ++j)
	    WaveFunction *=  ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = this->NbrParticlesPerColor; i < (NbrN1N2); ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < (NbrN1N2); ++j)
	    WaveFunction *= ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = NbrN1N2; i < this->NbrParticles; ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->NbrParticles; ++j)
	    WaveFunction *= ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->MIntra; ++i)
	TotalWaveFunction *= WaveFunction;
    }

  if (this->MInter > 0)
    {
      Complex WaveFunction(1.0);
      Complex WaveFunction2(1.0);
      for (int i = 0; i < this->NbrParticlesPerColor; ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  int Lim = i + this->NbrParticlesPerColor;
	  int j = this->NbrParticlesPerColor;
	  for (; j < Lim; ++j)
	    {
	      Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	      WaveFunction *= Tmp;
	      WaveFunction2 *= Tmp;
	    }
	  Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	  WaveFunction *= Tmp;
	  ++j;
	  for (; j < (NbrN1N2); ++j)
	    {
	      Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	      WaveFunction *= Tmp;
	      WaveFunction2 *= Tmp;
	    }
	}
      for (int i = 0; i < this->NbrParticlesPerColor; ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  int Lim = i + NbrN1N2;
	  int j = NbrN1N2;
	  for (; j < Lim; ++j)
	    {
	      Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	      WaveFunction *= Tmp;
	      WaveFunction2 *= Tmp;
	    }
	  Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	  WaveFunction *= Tmp;
	  ++j;
	  for (; j < this->NbrParticles; ++j)
	    {
	      Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	      WaveFunction *= Tmp;
	      WaveFunction2 *= Tmp;
	    }
	}
      for (int i = this->NbrParticlesPerColor; i < NbrN1N2; ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  int Lim = i + this->NbrParticlesPerColor;
	  int j = NbrN1N2;
	  for (; j < Lim; ++j)
	    {
	      Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	      WaveFunction *= Tmp;
	      WaveFunction2 *= Tmp;
	    }
	  Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	  WaveFunction *= Tmp;
	  ++j;
	  for (; j < this->NbrParticles; ++j)
	    {
	      Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	      WaveFunction *= Tmp;
	      WaveFunction2 *= Tmp;
	    }
	}
      TotalWaveFunction *= WaveFunction2;
      for (int i = 0; i < (this->MInter - 1); ++i)
	TotalWaveFunction *= WaveFunction;

    }

  return (TotalWaveFunction);
}
