////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of Laughlin wave function on sphere                 //
//                                                                            //
//                        last modification : 01/09/2004                      //
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
#include "Tools/FQHEWaveFunction/LaughlinOnSphereWaveFunction.h"
#include "Vector/RealVector.h"


// constructor
//
// nbrParticles = number of particles
// invFillingFactor = inverse value of the filling factor

LaughlinOnSphereWaveFunction::LaughlinOnSphereWaveFunction(int nbrParticles, int invFillingFactor)
{
  this->InvFillingFactor = invFillingFactor;
  this->NbrParticles = nbrParticles;
}

// copy constructor
//
// function = reference on the wave function to copy

LaughlinOnSphereWaveFunction::LaughlinOnSphereWaveFunction(const LaughlinOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->InvFillingFactor = function.InvFillingFactor;
}

// destructor
//

LaughlinOnSphereWaveFunction::~LaughlinOnSphereWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* LaughlinOnSphereWaveFunction::Clone ()
{
  return new LaughlinOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex LaughlinOnSphereWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp;
  Complex WaveFunction(1.0);
  double Theta;
  double Phi;
  double Factor = M_PI * 0.5;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      Theta = x[i << 1];
      Phi = x[1 + (i << 1)];
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp.Re = Factor * sin(0.5 * (x[j << 1] - Theta)) * cos(0.5 * (Phi - x[1 + (j << 1)]));
	  Tmp.Im = Factor * sin(0.5 * (Theta + x[j << 1])) * sin(0.5 * (Phi - x[1 + (j << 1)]));
	  WaveFunction *= Tmp;
	}
    }
  Tmp = WaveFunction;
  for (int i = 1; i < this->InvFillingFactor; ++i)
    {
      WaveFunction *= Tmp;
    }
  return WaveFunction;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex LaughlinOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  Complex Tmp;
  Complex WaveFunction(1.0);
  
  double Factor = M_PI * 0.5;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{	  
	  Tmp = Factor * ( uv[2*i] * uv[2*j+1] - uv[2*i+1] * uv[2*j] );
	  WaveFunction *= Tmp;
	}
    }
  Tmp = WaveFunction;
  for (int i = 1; i < this->InvFillingFactor; ++i)
    {
      WaveFunction *= Tmp;
    }
  return WaveFunction;
}
