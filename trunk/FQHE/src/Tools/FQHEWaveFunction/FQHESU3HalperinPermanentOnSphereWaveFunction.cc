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
#include "Tools/FQHEWaveFunction/FQHESU3HalperinPermanentOnSphereWaveFunction.h"
#include "GeneralTools/ConfigurationParser.h"
#include "Vector/RealVector.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>
#include <cstring> 

using std::cout;
using std::endl;


// constructor
//
// nbrN1 = number of type 1 particles (i.e. Tz=+1, Y=+1)
// nbrN2 = number of type 2 particles (i.e. Tz=-1, Y=+1)
// nbrN3 = number of type 3 particles (i.e. Tz=0, Y=-2)
// m11 = coefficient of the intra-component correlations in the 1 - 1 sector
// m22 = coefficient of the intra-component correlations in the 2 - 2 sector
// m33 = coefficient of the intra-component correlations in the 3 - 3 sector
// m12 = coefficient of the inter-component correlations in the 1 - 2 sector
// m13 = coefficient of the inter-component correlations in the 1 - 3 sector
// m23 = coefficient of the inter-component correlations in the 2 - 3 sector
  // invertFlag = if true, use the invert of the matrix elements for the permanents

FQHESU3HalperinPermanentOnSphereWaveFunction::FQHESU3HalperinPermanentOnSphereWaveFunction(int nbrN1, int nbrN2, int nbrN3,
											   int m11, int m22, int m33, int m12, int m13, int m23,
											   bool invertFlag)
{
  this->NbrN1 = nbrN1;
  this->NbrN2 = nbrN2;
  this->NbrN3 = nbrN3;  
  this->TotalNbrParticles = (this->NbrN1 + this->NbrN2 + this->NbrN3);
  this->M11 = m11;
  this->M22 = m22;
  this->M33 = m33;
  this->M12 = m12;
  this->M13 = m13;
  this->M23 = m23;
  this->InvertFlag = invertFlag;
  this->Permanent12 = ComplexMatrix(this->NbrN1, this->NbrN1, true);
  this->Permanent13 = ComplexMatrix(this->NbrN1, this->NbrN1, true);
  this->Permanent23 = ComplexMatrix(this->NbrN1, this->NbrN1, true);

}

// copy constructor
//
// function = reference on the wave function to copy

FQHESU3HalperinPermanentOnSphereWaveFunction::FQHESU3HalperinPermanentOnSphereWaveFunction(const FQHESU3HalperinPermanentOnSphereWaveFunction& function)
{
  this->NbrN1 = function.NbrN1;
  this->NbrN2 = function.NbrN2;
  this->NbrN3 = function.NbrN3;  
  this->TotalNbrParticles = (this->NbrN1 + this->NbrN2 + this->NbrN3);
  this->M11 = function.M11;
  this->M22 = function.M22;
  this->M33 = function.M33;
  this->M12 = function.M12;
  this->M13 = function.M13;
  this->M23 = function.M23;
  this->InvertFlag = function.InvertFlag;
  this->Permanent12 = ComplexMatrix(this->NbrN1, this->NbrN1, true);
  this->Permanent13 = ComplexMatrix(this->NbrN1, this->NbrN1, true);
  this->Permanent23 = ComplexMatrix(this->NbrN1, this->NbrN1, true);
}

// constructor from configuration file
//
// configuration = reference on the configuration file parser
// errorFlag = reference on the error flag that is set to false if an error occured while retriving the configuration
// nbrParticles = reference on the total number of particles (computed from the configuration file datas)
// lzMax = twice the maximum angular momentum for a particle (computed from the configuration file datas)

FQHESU3HalperinPermanentOnSphereWaveFunction::FQHESU3HalperinPermanentOnSphereWaveFunction(ConfigurationParser& configuration, bool& errorFlag, int& nbrParticles, int& lzMax)
{
  errorFlag = true;
  if ((configuration["WaveFunction"] == 0) || (strcmp ("SU3HalperinPermament", configuration["WaveFunction"]) != 0))
    {
      errorFlag = false;
      return;
    }
  if ((configuration.GetAsSingleInteger("NbrN1", this->NbrN1) == false) ||
      (configuration.GetAsSingleInteger("NbrN2", this->NbrN2) == false) || 
      (configuration.GetAsSingleInteger("NbrN3", this->NbrN3) == false))
    {
      errorFlag = false;
      return;
    }
  this->TotalNbrParticles = (this->NbrN1 + this->NbrN2 + this->NbrN3);
  if ((configuration.GetAsSingleInteger("M11", this->M11) == false) ||
      (configuration.GetAsSingleInteger("M22", this->M22) == false) ||
      (configuration.GetAsSingleInteger("M33", this->M33) == false) ||
      (configuration.GetAsSingleInteger("M12", this->M12) == false) ||
      (configuration.GetAsSingleInteger("M13", this->M13) == false) ||
      (configuration.GetAsSingleInteger("M23", this->M23) == false) ||
      (this->M11 < 0) || (this->M22 < 0) || (this->M33 < 0) || (this->M12 < 0) || 
      (this->M13 < 0) || (this->M23 < 0))
    {
      errorFlag = false;
      return;
    }
  this->InvertFlag = false;
  if ((configuration.GetAsBoolean("Invert", this->InvertFlag)) == false)
    this->InvertFlag = false;
  this->Permanent12 = ComplexMatrix(this->NbrN1, this->NbrN1, true);
  this->Permanent13 = ComplexMatrix(this->NbrN1, this->NbrN1, true);
  this->Permanent23 = ComplexMatrix(this->NbrN1, this->NbrN1, true);
}

// destructor
//

FQHESU3HalperinPermanentOnSphereWaveFunction::~FQHESU3HalperinPermanentOnSphereWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* FQHESU3HalperinPermanentOnSphereWaveFunction::Clone ()
{
  return new FQHESU3HalperinPermanentOnSphereWaveFunction(*this);
}

// evaluate function at a given point(the first 2*N1 coordinates correspond to the position of the type 1 particles, 
//                                     the following 2*N2 coordinates correspond to the position of the type 2 particles,
//                                     last the 2*N3 coordinates correspond to the position of the type 3 particles)
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)

Complex FQHESU3HalperinPermanentOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  Complex Tmp;
  Complex TmpU;
  Complex TmpV;
  int NbrN1N2 = this->NbrN1 + this->NbrN2;
  Complex TotalWaveFunction(1.0);
  if (this->M11 > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrN1; ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->NbrN1; ++j)
	    WaveFunction *=  ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->M11; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  if (this->M22 > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = this->NbrN1; i < (NbrN1N2); ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < (NbrN1N2); ++j)
	    WaveFunction *= ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->M22; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  if (this->M33 > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = NbrN1N2; i < this->TotalNbrParticles; ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->TotalNbrParticles; ++j)
	    WaveFunction *= ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->M33; ++i)
	TotalWaveFunction *= WaveFunction;
    }

  if (this->M12 > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrN1; ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  for (int j = this->NbrN1; j < (NbrN1N2); ++j)
	    {
	      Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	      WaveFunction *= Tmp;
	      if (this->InvertFlag == false)
		Tmp = Inv(Tmp);
	      Permanent12.SetMatrixElement(i, j - this->NbrN1, Tmp);
	    }
	}
      for (int i = 0; i < this->M12; ++i)
	TotalWaveFunction *= WaveFunction;
    }

  if (this->M13 > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrN1; ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  for (int j = NbrN1N2; j < this->TotalNbrParticles; ++j)
	    {
	      Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	      WaveFunction *= Tmp;
	      if (this->InvertFlag == false)
		Tmp = Inv(Tmp);
	      int ShiftedJ = j - NbrN1N2;
	      Permanent13.SetMatrixElement(i, ShiftedJ, Tmp);
	    }
	}
      for (int i = 0; i < this->M13; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  if (this->M23 > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = this->NbrN1; i < (NbrN1N2); ++i)
	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
	  for (int j = NbrN1N2; j < this->TotalNbrParticles; ++j)
	    {
	      Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	      WaveFunction *= Tmp;
	      if (this->InvertFlag == false)
		Tmp = Inv(Tmp);
	      int ShiftedI = i - this->NbrN1;	      
	      int ShiftedJ = j - NbrN1N2;
	      Permanent23.SetMatrixElement(ShiftedI, ShiftedJ, Tmp);
	    }
	}
      for (int i = 0; i < this->M23; ++i)
	TotalWaveFunction *= WaveFunction;
    }

  return (TotalWaveFunction * Permanent12.Permanent() * Permanent13.Permanent() * Permanent23.Permanent());
}
