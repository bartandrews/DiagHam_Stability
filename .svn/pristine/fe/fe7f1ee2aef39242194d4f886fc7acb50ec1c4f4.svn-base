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
  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];
}

// copy constructor
//
// function = reference on the wave function to copy

LaughlinOnSphereWaveFunction::LaughlinOnSphereWaveFunction(const LaughlinOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->InvFillingFactor = function.InvFillingFactor;
  this->SpinorUCoordinates = function.SpinorUCoordinates;
  this->SpinorVCoordinates = function.SpinorVCoordinates;

}

// destructor
//

LaughlinOnSphereWaveFunction::~LaughlinOnSphereWaveFunction()
{
  if (this->NbrParticles!=0)
    {
      delete [] SpinorUCoordinates;
      delete [] SpinorVCoordinates;
    }
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
//
Complex LaughlinOnSphereWaveFunction::operator ()(RealVector& x)
{
  double s,c;
  // Convert particle positions into spinor form:
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = cos(0.5 * x[i << 1]);
      this->SpinorUCoordinates[i].Im = this->SpinorUCoordinates[i].Re;
      this->SpinorUCoordinates[i].Re *= (c=cos(0.5 * x[1 + (i << 1)]));
      this->SpinorUCoordinates[i].Im *= -(s=sin(0.5 * x[1 + (i << 1)]));
      this->SpinorVCoordinates[i].Re = sin(0.5 * x[i << 1]);
      this->SpinorVCoordinates[i].Im = this->SpinorVCoordinates[i].Re;
      this->SpinorVCoordinates[i].Re *= c;
      this->SpinorVCoordinates[i].Im *= s;
    }

  Complex Result(1.0);
  
  for (int i=1;i<this->NbrParticles;i++)
    for (int j=0;j<i;j++)
    {
      Result *= (SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
    }
  Complex Base=Result;
  for (int i=1; i<this->InvFillingFactor; ++i)
    Result *= Base;
  return Result;
}


// evaluate function at a given point
// 
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex LaughlinOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  // Import from spinors
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = uv.Re(2*i);
      this->SpinorUCoordinates[i].Im = uv.Im(2*i);
      this->SpinorVCoordinates[i].Re = uv.Re(2*i+1);
      this->SpinorVCoordinates[i].Im = uv.Im(2*i+1);
    }
  Complex Result(1.0);
  
  for (int i=1;i<this->NbrParticles;i++)
    for (int j=0;j<i;j++)
    {
      Result *= (SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
    }
  Complex Base=Result;
  for (int i=1; i<this->InvFillingFactor; ++i)
    Result *= Base;
  return Result;
  
}

// Complex LaughlinOnSphereWaveFunction::operator ()(RealVector& x)
// {
//   Complex Tmp;
//   Complex WaveFunction(1.0);
//   double Theta;
//   double Phi;
//   double Factor = M_PI * 0.5;
//   for (int i = 0; i < this->NbrParticles; ++i)
//     {
//       Theta = x[i << 1];
//       Phi = x[1 + (i << 1)];
//       for (int j = i + 1; j < this->NbrParticles; ++j)
// 	{
// 	  Tmp.Re = Factor * sin(0.5 * (x[j << 1] - Theta)) * cos(0.5 * (Phi - x[1 + (j << 1)]));
// 	  Tmp.Im = Factor * sin(0.5 * (Theta + x[j << 1])) * sin(0.5 * (Phi - x[1 + (j << 1)]));
// 	  WaveFunction *= Tmp;
// 	}
//     }
//   Tmp = WaveFunction;
//   for (int i = 1; i < this->InvFillingFactor; ++i)
//     {
//       WaveFunction *= Tmp;
//     }
//   return WaveFunction;
// }

// // evaluate function at a given point
// //
// // uv = ensemble of spinor variables on sphere describing point
// //      where function has to be evaluated
// //      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// // return value = function value at (uv)
// Complex LaughlinOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
// {
//   Complex Tmp;
//   Complex WaveFunction(1.0);
  
//   double Factor = M_PI * 0.5;
//   for (int i = 0; i < this->NbrParticles; ++i)
//     {
//       for (int j = i + 1; j < this->NbrParticles; ++j)
// 	{	  
// 	  Tmp = Factor * ( uv[2*i] * uv[2*j+1] - uv[2*i+1] * uv[2*j] );
// 	  WaveFunction *= Tmp;
// 	}
//     }
//   Tmp = WaveFunction;
//   for (int i = 1; i < this->InvFillingFactor; ++i)
//     {
//       WaveFunction *= Tmp;
//     }
//   return WaveFunction;
// }
