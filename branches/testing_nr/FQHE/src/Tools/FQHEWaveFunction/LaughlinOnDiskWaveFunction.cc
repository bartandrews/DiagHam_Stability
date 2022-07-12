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


// constructor
//
// nbrParticles = number of particles
// invFillingFactor = inverse value of the filling factor
// scale = typical sytem size

LaughlinOnDiskWaveFunction::LaughlinOnDiskWaveFunction(int nbrParticles, int invFillingFactor, double scale)
{
  this->InvFillingFactor = invFillingFactor;
  this->NbrParticles = nbrParticles;
  this->InvScale = 1.0 / scale;
}

// copy constructor
//
// function = reference on the wave function to copy

LaughlinOnDiskWaveFunction::LaughlinOnDiskWaveFunction(const LaughlinOnDiskWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->InvFillingFactor = function.InvFillingFactor;
  this->InvScale = function.InvScale;
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
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      ZRe = x[i << 1];
      ZIm = x[1 + (i << 1)];
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
  return WaveFunction;
}
