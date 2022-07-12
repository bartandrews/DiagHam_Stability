////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of Laughlin wave function on disk with one quasihole           //
//                         (without the gaussian factor)                      //
//                                                                            //
//                        last modification : 15/11/2008                      //
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
#include "Tools/FQHEWaveFunction/FQHEDiskLaughlinOneQuasiholeWaveFunction.h"
#include "Vector/RealVector.h"


// constructor
//
// nbrParticles = number of particles
// zHole = position of the quasihole
// invFillingFactor = inverse value of the filling factor

FQHEDiskLaughlinOneQuasiholeWaveFunction::FQHEDiskLaughlinOneQuasiholeWaveFunction(int nbrParticles, Complex zHole, int invFillingFactor)
{
  this->InvFillingFactor = invFillingFactor;
  this->NbrParticles = nbrParticles;
  this->ZHole = zHole;
}

// copy constructor
//
// function = reference on the wave function to copy

FQHEDiskLaughlinOneQuasiholeWaveFunction::FQHEDiskLaughlinOneQuasiholeWaveFunction(const FQHEDiskLaughlinOneQuasiholeWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->InvFillingFactor = function.InvFillingFactor;
  this->ZHole = function.ZHole;
}

// destructor
//

FQHEDiskLaughlinOneQuasiholeWaveFunction::~FQHEDiskLaughlinOneQuasiholeWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* FQHEDiskLaughlinOneQuasiholeWaveFunction::Clone ()
{
  return new FQHEDiskLaughlinOneQuasiholeWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex FQHEDiskLaughlinOneQuasiholeWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp;
  Complex WaveFunction(1.0);
  double ZRe;
  double ZIm;
  Complex HoleFactor (1.0);
  double Scale = 1.0 / 15;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      ZRe = x[i << 1];
      ZIm = x[1 + (i << 1)];
      Tmp.Re = ZRe - this->ZHole.Re;
      Tmp.Im = ZIm - this->ZHole.Im;
      HoleFactor *= Scale * Tmp;
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp.Re = ZRe - x[j << 1];
	  Tmp.Im = ZIm - x[1 + (j << 1)];
	  WaveFunction *= Scale * Tmp;
	}
    }
  Tmp = WaveFunction;
  for (int i = 1; i < this->InvFillingFactor; ++i)
    {
      WaveFunction *= Tmp;
    }
  WaveFunction *= HoleFactor;
  return WaveFunction;
}
