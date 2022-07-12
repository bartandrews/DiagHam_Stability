////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of Halperin (m1, m2, n) wave function on disk             //
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
#include "Tools/FQHEWaveFunction/HalperinOnDiskWaveFunctionOneOverR.h"
#include "Vector/RealVector.h"


// constructor
//
// nbrSpinUpParticles = number of particles with spin up
// nbrSpinDownParticles = number of particles with spin down
// m1Index = m1 index
// m2Index = m2 index
// nIndex = n index

HalperinOnDiskWaveFunctionOneOverR::HalperinOnDiskWaveFunctionOneOverR(int nbrSpinUpParticles, int nbrSpinDownParticles, int m1Index, int m2Index, int nIndex)
{
  this->M1Index = m1Index;
  this->M2Index = m2Index;
  this->NIndex = nIndex;
  this->NbrSpinUpParticles = nbrSpinUpParticles;
  this->NbrSpinDownParticles = nbrSpinDownParticles;
  this->TotalNbrParticles = this->NbrSpinUpParticles + this->NbrSpinDownParticles;
  this->NbrFactors = ((this->TotalNbrParticles - 1) * this->TotalNbrParticles) >> 1;
  this->JastrowFactors = new Complex [this->NbrFactors];
  this->TemporaryFactors = new double [this->NbrFactors];
}

// copy constructor
//
// function = reference on the wave function to copy

HalperinOnDiskWaveFunctionOneOverR::HalperinOnDiskWaveFunctionOneOverR(const HalperinOnDiskWaveFunctionOneOverR& function)
{
  this->M1Index = function.M1Index;
  this->M2Index = function.M2Index;
  this->NIndex = function.NIndex;
  this->NbrSpinUpParticles = function.NbrSpinUpParticles;
  this->NbrSpinDownParticles = function.NbrSpinDownParticles;
  this->TotalNbrParticles = function.TotalNbrParticles;
  this->NbrFactors = function.NbrFactors;
  this->JastrowFactors = new Complex [this->NbrFactors];
  for (int i = 0; i < this->NbrFactors; ++i)
    this->JastrowFactors[i] = function.JastrowFactors[i];
  this->TemporaryFactors = new double [this->NbrFactors];
}

// destructor
//

HalperinOnDiskWaveFunctionOneOverR::~HalperinOnDiskWaveFunctionOneOverR()
{
  delete[] this->JastrowFactors;
  delete[] this->TemporaryFactors;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* HalperinOnDiskWaveFunctionOneOverR::Clone ()
{
  return new HalperinOnDiskWaveFunctionOneOverR(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex HalperinOnDiskWaveFunctionOneOverR::operator ()(RealVector& x)
{
  Complex Tmp;
  Complex WaveFunction(1.0);
  Complex WaveFunction2(1.0);
  double ZRe;
  double ZIm;
  int Pos = 0;
  for (int i = 0; i < this->NbrSpinUpParticles; ++i)
    {
      ZRe = x[i << 1];
      ZIm = x[1 + (i << 1)];
      for (int j = i + 1; j < this->NbrSpinUpParticles; ++j)
	{
	  Tmp.Re = ZRe - x[j << 1];
	  Tmp.Im = ZIm - x[1 + (j << 1)];
	  this->JastrowFactors[Pos] = Tmp;
	  ++Pos;
	  WaveFunction *= Tmp;
	}
    }
  for (int i = 0; i < this->M1Index; ++i)
    WaveFunction2 *= WaveFunction;

  for (int i = this->NbrSpinUpParticles; i < this->TotalNbrParticles; ++i)
    {
      ZRe = x[i << 1];
      ZIm = x[1 + (i << 1)];
      for (int j = i + 1; j < this->TotalNbrParticles; ++j)
	{
	  Tmp.Re = ZRe - x[j << 1];
	  Tmp.Im = ZIm - x[1 + (j << 1)];
	  this->JastrowFactors[Pos] = Tmp;
	  ++Pos;
	  WaveFunction *= Tmp;
	}
    }
  for (int i = 0; i < this->M2Index; ++i)
    WaveFunction2 *= WaveFunction;

  WaveFunction = 1.0;
  for (int i = 0; i < this->NbrSpinUpParticles; ++i)
    {
      ZRe = x[i << 1];
      ZIm = x[1 + (i << 1)];
      for (int j = this->NbrSpinUpParticles; j < this->TotalNbrParticles; ++j)
	{
	  Tmp.Re = ZRe - x[j << 1];
	  Tmp.Im = ZIm - x[1 + (j << 1)];
	  this->JastrowFactors[Pos] = Tmp;
	  ++Pos;
	  WaveFunction *= Tmp;
	}
    }
  for (int i = 1; i < this->NIndex; ++i)
    WaveFunction2 *= WaveFunction;

  return WaveFunction2;
}

// evaluate the norm to the square of the wave function at a given point time the coulomb term (assume the coordinates are those provides by the previous operator() method call)
//
// return value = corresponding numerical value

double HalperinOnDiskWaveFunctionOneOverR::CoulombContribution()
{
  for (int i = 0; i < this->NbrFactors; ++i)
    this->TemporaryFactors[i] = 1.0;
  double TmpNorm;
  double Tmp;
  int k;  
  double TwiceM1 = 2.0 * this->M1Index;
  double TwiceM2 = 2.0 * this->M2Index;
  double TwiceN = 2.0 * this->NIndex;
  int Lim1 = ((this->NbrSpinUpParticles - 1) * this->NbrSpinUpParticles) >> 1;
  int Lim2 = (((this->NbrSpinDownParticles - 1) * this->NbrSpinDownParticles) >> 1) + Lim2;

  int i = 0;
  for (; i < Lim1; ++i)
    {
      TmpNorm = Norm(this->JastrowFactors[i]);
      Tmp = pow(TmpNorm, TwiceM1);
      for (k = 0; k < i; ++k)
	this->TemporaryFactors[k] *= Tmp;
      this->TemporaryFactors[k] *= pow(TmpNorm, TwiceM1 - 1.0);
      ++k;
      for (; k < this->NbrFactors; ++k)
	this->TemporaryFactors[k] *= Tmp;	
    }
  for (; i < Lim2; ++i)
    {
      TmpNorm = Norm(this->JastrowFactors[i]);
      Tmp = pow(TmpNorm, TwiceM1);
      for (k = 0; k < i; ++k)
	this->TemporaryFactors[k] *= Tmp;
      this->TemporaryFactors[k] *= pow(TmpNorm, TwiceM2 - 1.0);
      ++k;
      for (; k < this->NbrFactors; ++k)
	this->TemporaryFactors[k] *= Tmp;	
    }
  for (; i < this->NbrFactors; ++i)
    {
      TmpNorm = Norm(this->JastrowFactors[i]);
      Tmp = pow(TmpNorm, TwiceM1);
      for (k = 0; k < i; ++k)
	this->TemporaryFactors[k] *= Tmp;
      this->TemporaryFactors[k] *= pow(TmpNorm, TwiceN - 1.0);
      ++k;
      for (; k < this->NbrFactors; ++k)
	this->TemporaryFactors[k] *= Tmp;	
    }
  Tmp =  this->TemporaryFactors[0];
  for (int i = 1; i < this->NbrFactors; ++i)
    Tmp += this->TemporaryFactors[i];
  return Tmp;
}

