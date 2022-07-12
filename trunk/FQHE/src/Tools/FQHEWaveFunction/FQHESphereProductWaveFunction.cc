////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of wave function obtained from product of            //
//                            wave function on sphere                         //
//                                                                            //
//                        last modification : 28/08/2009                      //
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
#include "Tools/FQHEWaveFunction/FQHESphereProductWaveFunction.h"
#include "Vector/RealVector.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
//

FQHESphereProductWaveFunction::FQHESphereProductWaveFunction()
{
}

// constructor
//
// nbrParticles = number of particles
// wavefunctions = array to the wavefunctions that have to be multiplied
// nbrWaveFunctions = number of wavefunctions to multiply
// jastrowFactor = multiply(if positive) or divide (if negative) the wavefunction by an overall Jastrow factor
// separateCoordinates = if false use the same coordinates to evaluate wavefunctions, if true use coordinates 0 to nbrParticles/nbrWaveFunctions - 1 for the first wavefunctions and so on
 
FQHESphereProductWaveFunction::FQHESphereProductWaveFunction(int nbrParticles, Abstract1DComplexFunctionOnSphere** waveFunctions, int nbrWaveFunctions, int jastrowFactor, bool separateCoordinates)
{
  this->NbrParticles = nbrParticles;
  this->NbrWaveFunctions = nbrWaveFunctions;
  this->BaseWaveFunctions = new Abstract1DComplexFunctionOnSphere*[this->NbrWaveFunctions];
  for (int i = 0; i < this->NbrWaveFunctions; ++i)
    this->BaseWaveFunctions[i] =  (Abstract1DComplexFunctionOnSphere*) (waveFunctions[i]->Clone());
  this->JastrowFactor = jastrowFactor;
  this->SeparateCoordinateFlag = separateCoordinates;
  this->TemporaryUV = ComplexVector(this->NbrParticles * 2);
  if (this->SeparateCoordinateFlag == true)
    this->TemporaryUVSeparateCoordinates = ComplexVector((this->NbrParticles * 2) / this->NbrWaveFunctions);
  else
    this->TemporaryUVSeparateCoordinates = ComplexVector();
}

// copy constructor
//
// function = reference on the wave function to copy

FQHESphereProductWaveFunction::FQHESphereProductWaveFunction(const FQHESphereProductWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrWaveFunctions = function.NbrWaveFunctions;
  this->BaseWaveFunctions = new Abstract1DComplexFunctionOnSphere*[this->NbrWaveFunctions];
  for (int i = 0; i < this->NbrWaveFunctions; ++i)
    this->BaseWaveFunctions[i] =  (Abstract1DComplexFunctionOnSphere*) (function.BaseWaveFunctions[i]->Clone());
  this->JastrowFactor = function.JastrowFactor;
  this->SeparateCoordinateFlag = function.SeparateCoordinateFlag;
  this->TemporaryUV = ComplexVector(this->NbrParticles * 2);
  if (this->SeparateCoordinateFlag == true)
    this->TemporaryUVSeparateCoordinates = ComplexVector((this->NbrParticles * 2) / this->NbrWaveFunctions);
  else
    this->TemporaryUVSeparateCoordinates = ComplexVector();
}

// destructor
//

FQHESphereProductWaveFunction::~FQHESphereProductWaveFunction()
{
  delete[] this->BaseWaveFunctions;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* FQHESphereProductWaveFunction::Clone ()
{
  return new FQHESphereProductWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex FQHESphereProductWaveFunction::operator ()(RealVector& x)
{  
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->TemporaryUV[2 * i].Re = cos(0.5 * x[i << 1]);
      this->TemporaryUV[2 * i].Im = this->TemporaryUV[2 * i].Re;
      this->TemporaryUV[2 * i].Re *= cos(0.5 * x[1 + (i << 1)]);
      this->TemporaryUV[2 * i].Im *= sin(0.5 * x[1 + (i << 1)]);
      this->TemporaryUV[(2 * i) + 1].Re = sin(0.5 * x[i << 1]);
      this->TemporaryUV[(2 * i) + 1].Im = this->TemporaryUV[(2 * i) + 1].Re;
      this->TemporaryUV[(2 * i) + 1].Re *= cos(0.5 * x[1 + (i << 1)]);
      this->TemporaryUV[(2 * i) + 1].Im *= -sin(0.5 * x[1 + (i << 1)]);
    }
  return this->CalculateFromSpinorVariables(this->TemporaryUV);
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)

Complex FQHESphereProductWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{  
  Complex TotalValue = 1.0;
  if (this->SeparateCoordinateFlag == false)
    {
      for (int i = 0; i < this->NbrWaveFunctions; ++i)
	TotalValue *= this->BaseWaveFunctions[i]->CalculateFromSpinorVariables(uv);  
    }
  else
    {
      int Lim = (this->NbrParticles * 2) / this->NbrWaveFunctions;	
      for (int i = 0; i < this->NbrWaveFunctions; ++i)
	{
	  for (int j = 0; j < Lim; ++j)	    
	    this->TemporaryUVSeparateCoordinates[j] = uv[(Lim * i) + j];
	  TotalValue *= this->BaseWaveFunctions[i]->CalculateFromSpinorVariables(this->TemporaryUVSeparateCoordinates); 
	}
    }
  if (this->JastrowFactor != 0)
    {
      Complex TmpU;
      Complex TmpV;
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->NbrParticles; ++j)
	    WaveFunction *=  ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      if (this->JastrowFactor > 0)
	for (int i = 0; i < this->JastrowFactor; ++i)
	  TotalValue *= WaveFunction;
      else
	for (int i = 0; i > this->JastrowFactor; --i)
	  TotalValue /= WaveFunction;
    }
  return TotalValue;
}  


