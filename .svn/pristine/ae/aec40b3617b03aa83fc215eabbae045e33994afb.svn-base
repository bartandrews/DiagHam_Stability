////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of Laughlin wave function on disk                  //
//              that evalute both the wave function and the wave              //
//                  function time the Coulomb interaction term                //
//                         (without the gaussian factor)                      //
//                                                                            //
//                        last modification : 23/10/2006                      //
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
#include "Tools/FQHEWaveFunction/LaughlinOnDiskWaveFunctionOneOverR.h"
#include "Vector/RealVector.h"


// constructor
//
// nbrParticles = number of particles
// invFillingFactor = inverse value of the filling factor

LaughlinOnDiskWaveFunctionOneOverR::LaughlinOnDiskWaveFunctionOneOverR(int nbrParticles, int invFillingFactor)
{
  this->InvFillingFactor = invFillingFactor;
  this->NbrParticles = nbrParticles;
  this->NbrFactors = ((this->NbrParticles - 1) * this->NbrParticles) >> 1;
  this->JastrowFactors = new Complex [this->NbrFactors];
  this->TemporaryFactors = new double [this->NbrFactors];
}

// copy constructor
//
// function = reference on the wave function to copy

LaughlinOnDiskWaveFunctionOneOverR::LaughlinOnDiskWaveFunctionOneOverR(const LaughlinOnDiskWaveFunctionOneOverR& function)
{
  this->NbrParticles = function.NbrParticles;
  this->InvFillingFactor = function.InvFillingFactor;
  this->NbrFactors = function.NbrFactors;
  this->JastrowFactors = new Complex [this->NbrFactors];
  for (int i = 0; i < this->NbrFactors; ++i)
    this->JastrowFactors[i] = function.JastrowFactors[i];
  this->TemporaryFactors = new double [this->NbrFactors];
}

// destructor
//

LaughlinOnDiskWaveFunctionOneOverR::~LaughlinOnDiskWaveFunctionOneOverR()
{
  delete[] this->JastrowFactors;
  delete[] this->TemporaryFactors;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* LaughlinOnDiskWaveFunctionOneOverR::Clone ()
{
  return new LaughlinOnDiskWaveFunctionOneOverR(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex LaughlinOnDiskWaveFunctionOneOverR::operator ()(RealVector& x)
{
  Complex Tmp;
  Complex WaveFunction(1.0);
  double ZRe;
  double ZIm;
  int Pos = 0;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      ZRe = x[i << 1];
      ZIm = x[1 + (i << 1)];
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp.Re = ZRe - x[j << 1];
	  Tmp.Im = ZIm - x[1 + (j << 1)];
	  this->JastrowFactors[Pos] = Tmp;
	  ++Pos;
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

// evaluate the norm to the square of the wave function at a given point time the coulomb term (assume the coordinates are those provides by the previous operator() method call)
//
// return value = corresponding numerical value

double LaughlinOnDiskWaveFunctionOneOverR::CoulombContribution()
{
  for (int i = 0; i < this->NbrFactors; ++i)
    this->TemporaryFactors[i] = 1.0;
  double TmpNorm;
  double Tmp;
  int k;  
  double TwiceInvFillingFactor = 2.0 * this->InvFillingFactor;
  for (int i = 0; i < this->NbrFactors; ++i)
    {
      TmpNorm = Norm(this->JastrowFactors[i]);
      Tmp = pow(TmpNorm, TwiceInvFillingFactor);
      for (k = 0; k < i; ++k)
	this->TemporaryFactors[k] *= Tmp;
      this->TemporaryFactors[k] *= pow(TmpNorm, TwiceInvFillingFactor - 1.0);
      ++k;
      for (; k < this->NbrFactors; ++k)
	this->TemporaryFactors[k] *= Tmp;	
    }
  Tmp =  this->TemporaryFactors[0];
  for (int i = 1; i < this->NbrFactors; ++i)
    Tmp += this->TemporaryFactors[i];
  return Tmp;
}

