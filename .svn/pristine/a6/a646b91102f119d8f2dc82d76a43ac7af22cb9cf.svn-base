////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of SU(2) Halperin wave function on sphere             //
//                            times the permament                             //
//                                                                            //
//                        last modification : 09/04/2008                      //
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
#include "Tools/FQHEWaveFunction/FQHESU2HalperinPermanentOnSphereWaveFunction.h"
#include "Vector/RealVector.h"


// constructor
//
// nbrSpinUpParticles = number of particles with spin up
// nbrSpinDownParticles = number of particles with spin down
// m1Index = m1 index
// m2Index = m2 index
// nIndex = n index
// invertFlag = if true, use the invert of the matrix elements for the permanent

FQHESU2HalperinPermanentOnSphereWaveFunction::FQHESU2HalperinPermanentOnSphereWaveFunction(int nbrSpinUpParticles, int nbrSpinDownParticles, 
											   int m1Index, int m2Index, int nIndex, bool invertFlag)
{
  this->M1Index = m1Index;
  this->M2Index = m2Index;
  this->NIndex = nIndex;
  this->NbrSpinUpParticles = nbrSpinUpParticles;
  this->NbrSpinDownParticles = nbrSpinDownParticles;
  this->TotalNbrParticles = this->NbrSpinUpParticles + this->NbrSpinDownParticles;
  this->UCoordinates = new Complex[this->TotalNbrParticles];
  this->VCoordinates = new Complex[this->TotalNbrParticles];
  this->InvertFlag = invertFlag;
  this->Permanent12 = ComplexMatrix (this->NbrSpinUpParticles, this->NbrSpinUpParticles, true);
}

// copy constructor
//
// function = reference on the wave function to copy

FQHESU2HalperinPermanentOnSphereWaveFunction::FQHESU2HalperinPermanentOnSphereWaveFunction(const FQHESU2HalperinPermanentOnSphereWaveFunction& function)
{
  this->M1Index = function.M1Index;
  this->M2Index = function.M2Index;
  this->NIndex = function.NIndex;
  this->NbrSpinUpParticles = function.NbrSpinUpParticles;
  this->NbrSpinDownParticles = function.NbrSpinDownParticles;
  this->TotalNbrParticles = function.TotalNbrParticles;
  this->UCoordinates = new Complex[this->TotalNbrParticles];
  this->VCoordinates = new Complex[this->TotalNbrParticles];
  this->InvertFlag = function.InvertFlag;
  this->Permanent12 = ComplexMatrix (this->NbrSpinUpParticles, this->NbrSpinUpParticles, true);
}

// destructor
//

FQHESU2HalperinPermanentOnSphereWaveFunction::~FQHESU2HalperinPermanentOnSphereWaveFunction()
{
  delete[] this->UCoordinates;
  delete[] this->VCoordinates;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* FQHESU2HalperinPermanentOnSphereWaveFunction::Clone ()
{
  return new FQHESU2HalperinPermanentOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)

Complex FQHESU2HalperinPermanentOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  Complex TotalWaveFunction(1.0);
  Complex TmpU;
  Complex TmpV;
  Complex Tmp;
  if (this->M1Index > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrSpinUpParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->NbrSpinUpParticles; ++j)
	    WaveFunction *=  ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->M1Index; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  if (this->M2Index > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = this->NbrSpinUpParticles; i < this->TotalNbrParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->TotalNbrParticles; ++j)
	    WaveFunction *=  ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->M2Index; ++i)
	TotalWaveFunction *= WaveFunction;
    }
   if (this->NIndex > 0)
     {
       Complex WaveFunction(1.0);
       for (int i = 0; i < this->NbrSpinUpParticles; ++i)
 	{
 	  TmpU = uv[2 * i];
 	  TmpV = uv[2 * i + 1];
 	  for (int j = this->NbrSpinUpParticles; j < this->TotalNbrParticles; ++j)
	    {
	      Tmp = ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	      WaveFunction *= Tmp;
	      if (this->InvertFlag == false)
		Tmp = Inv(Tmp);
	      this->Permanent12.SetMatrixElement(i, j - this->NbrSpinUpParticles, Tmp);
	    }	  
 	}
       for (int i = 0; i < this->NIndex; ++i)
	 TotalWaveFunction *= WaveFunction;
       TotalWaveFunction *= this->Permanent12.Permanent();
     }
  return TotalWaveFunction;
}

