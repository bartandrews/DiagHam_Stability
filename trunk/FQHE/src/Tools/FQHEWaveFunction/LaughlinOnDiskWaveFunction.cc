////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of Laughlin wave function on disk                  //
//                         (without the gaussian factor)                      //
//                                                                            //
//                        last modification : 10/10/2004                      //
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
#include "Tools/FQHEWaveFunction/LaughlinOnDiskWaveFunction.h"
#include "Vector/RealVector.h"
#include <cmath>


// constructor
//
// nbrParticles = number of particles
// invFillingFactor = inverse value of the filling factor
// scale = typical sytem size

LaughlinOnDiskWaveFunction::LaughlinOnDiskWaveFunction(int nbrParticles, int invFillingFactor, double scale, bool useExponentials)
{
  this->InvFillingFactor = invFillingFactor;
  this->NbrParticles = nbrParticles;
  if (scale==1.0)
    {
      this->InvScale=0.33/this->InvFillingFactor;
      std::cout << "Laughlin state: guessing inverse scaling factor: "<<this->InvScale<<std::endl;
    }
  else
    this->InvScale = 1.0 / scale;
  this->ExponentialFactors = useExponentials;
  this->LogScale = 0.25*this->NbrParticles*(this->NbrParticles-1)*this->InvFillingFactor;
  this->SumSqr = 2.0*this->LogScale;
  // std::cout << "InvScale=" << InvScale<<", Log-Scale="<<this->LogScale<<std::endl;
}

// copy constructor
//
// function = reference on the wave function to copy

LaughlinOnDiskWaveFunction::LaughlinOnDiskWaveFunction(const LaughlinOnDiskWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->InvFillingFactor = function.InvFillingFactor;
  this->InvScale = function.InvScale;
  this->ExponentialFactors = function.ExponentialFactors;
  this->LogScale = function.LogScale;
  this->SumSqr = function.SumSqr;
}

// destructor
//

LaughlinOnDiskWaveFunction::~LaughlinOnDiskWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* LaughlinOnDiskWaveFunction::Clone ()
{
  return new LaughlinOnDiskWaveFunction(*this);
}


// change the normalization of the funtion by a multiplicative factor
// factor = factor to be multiplied
void LaughlinOnDiskWaveFunction::Renormalize(double factor)
{
  // std::cout << "Renormalize ("<<factor<<")\n";
  this->InvScale *= std::exp(1.0/(this->NbrParticles*(this->NbrParticles-1))*(this->LogScale-0.5*this->SumSqr+std::log(factor)));
  this->LogScale = 0.5*this->SumSqr;
}


// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex LaughlinOnDiskWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp;
  Complex WaveFunction(1.0);
  double ZRe;
  double ZIm;
  this->SumSqr=0.0;
  //std::cout << "x="<<x<<"\n";
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      ZRe = x[i << 1];
      ZIm = x[1 + (i << 1)];
      SumSqr += ZRe*ZRe + ZIm*ZIm;
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp.Re = ZRe - x[j << 1];
	  Tmp.Im = ZIm - x[1 + (j << 1)];
	  Tmp *= this->InvScale;
	  WaveFunction *= Tmp;
	}
    }
  Tmp = WaveFunction;
  for (int i = 1; i < this->InvFillingFactor; ++i)
    {
      WaveFunction *= Tmp;
    }
  if (this->ExponentialFactors)
    {
      //std::cout << "SumSqr="<<SumSqr<<" exponent="<<-0.5*SumSqr+this->LogScale<<" WaveFunction="<<WaveFunction<<", Invscale="<<this->InvScale << std::endl;
      WaveFunction *= std::exp(-0.5*SumSqr+this->LogScale);
      //exit(1);

    }
  return WaveFunction;
}
