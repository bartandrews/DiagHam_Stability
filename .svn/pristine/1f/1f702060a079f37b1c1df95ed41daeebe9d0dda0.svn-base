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
#include "Tools/FQHEWaveFunction/FQHEDiskLaughlinOneJainQuasielectronWaveFunction.h"
#include "Vector/RealVector.h"


// constructor
//
// nbrParticles = number of particles
// zElectron = quasielectron position
// scale = typical sytem size
// invFillingFactor = inverse value of the filling factor

FQHEDiskLaughlinOneJainQuasielectronWaveFunction::FQHEDiskLaughlinOneJainQuasielectronWaveFunction(int nbrParticles, Complex zElectron, double scale, int invFillingFactor)
{
  this->InvFillingFactor = invFillingFactor;
  this->NbrParticles = nbrParticles;
  this->InvScale = 2.0 / scale;
  this->ZElectron = zElectron;  
  this->ZElectron /= sqrt((double) invFillingFactor);
  this->ElectronWeight = exp (-0.5 * SqrNorm(this->ZElectron));
  this->TmpElectronWeight = new Complex [this->NbrParticles];
  this->TmpJastrow = new Complex*[this->NbrParticles];
  this->TmpSqrJastrow = new Complex*[this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->TmpJastrow[i] = new Complex [this->NbrParticles];
      this->TmpSqrJastrow[i] = new Complex [this->NbrParticles];
    }
}

// copy constructor
//
// function = reference on the wave function to copy

FQHEDiskLaughlinOneJainQuasielectronWaveFunction::FQHEDiskLaughlinOneJainQuasielectronWaveFunction(const FQHEDiskLaughlinOneJainQuasielectronWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->InvFillingFactor = function.InvFillingFactor;
  this->InvScale = function.InvScale;
  this->ZElectron = function.ZElectron;
  this->ElectronWeight = exp (-0.125 * SqrNorm(this->ZElectron));
  this->TmpElectronWeight = new Complex [this->NbrParticles];
  this->TmpJastrow = new Complex*[this->NbrParticles];
  this->TmpSqrJastrow = new Complex*[this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->TmpJastrow[i] = new Complex [this->NbrParticles];
      this->TmpSqrJastrow[i] = new Complex [this->NbrParticles];
    }
}

// destructor
//

FQHEDiskLaughlinOneJainQuasielectronWaveFunction::~FQHEDiskLaughlinOneJainQuasielectronWaveFunction()
{
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      delete[] this->TmpSqrJastrow[i];
      delete[] this->TmpJastrow[i];
    }
  delete[] this->TmpSqrJastrow;
  delete[] this->TmpJastrow;
  delete[] this->TmpElectronWeight;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* FQHEDiskLaughlinOneJainQuasielectronWaveFunction::Clone ()
{
  return new FQHEDiskLaughlinOneJainQuasielectronWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex FQHEDiskLaughlinOneJainQuasielectronWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp;
  Complex Tmp2;
  Complex Tmp3;
  Complex WaveFunction(1.0);
  double ZRe;
  double ZIm;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      ZRe = x[i << 1];
      ZIm = x[1 + (i << 1)];
      Tmp.Re = (this->ZElectron.Re * ZRe + this->ZElectron.Im * ZIm);
      Tmp.Im = (this->ZElectron.Re * ZIm - this->ZElectron.Im * ZRe);
      this->TmpElectronWeight[i] = exp (Tmp);
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp.Re = ZRe - x[j << 1];
	  Tmp.Im = ZIm - x[1 + (j << 1)];
	  Tmp *= this->InvScale;
	  this->TmpJastrow[i][j] = Tmp;
	  this->TmpJastrow[j][i] = -Tmp;
	  this->TmpSqrJastrow[i][j] = Tmp * Tmp;
	  this->TmpSqrJastrow[j][i] = this->TmpSqrJastrow[i][j];
	  WaveFunction *= Tmp;
	}
    }

  Complex WaveFunction2(0.0);  
  for (int i =0; i < this->NbrParticles; ++i)
    {
      Tmp = this->TmpElectronWeight[i];
      Tmp *= this->ElectronWeight;
      Tmp2 = 0.0;
      for (int j = 0; j < i; ++j)
	{
	  Tmp3 = 1.0;
	  for (int k = 0; k < j; ++k)
	    Tmp3 *= this->TmpJastrow[i][k];	    
	  for (int k = j + 1; k < i; ++k)	
	    {
	      Tmp *= this->TmpSqrJastrow[j][k];
	      Tmp3 *= this->TmpJastrow[i][k];	    
	    }
	  for (int k = i + 1; k < this->NbrParticles; ++k)
	    {	
	      Tmp *= this->TmpSqrJastrow[j][k];
	      Tmp3 *= this->TmpJastrow[i][k];
	    }
	  Tmp2 += Tmp3;
	}
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp3 = 1.0;
	  for (int k = 0; k < i; ++k)
	    Tmp3 *= this->TmpJastrow[i][k];	    
	  for (int k = i + 1; k < j; ++k)
	    Tmp3 *= this->TmpJastrow[i][k];	    
	  for (int k = j + 1; k < this->NbrParticles; ++k)	
	    {
	      Tmp *= this->TmpSqrJastrow[j][k];
	      Tmp3 *= this->TmpJastrow[i][k];
	    }
	  Tmp2 += Tmp3;
	}
      Tmp *= Tmp2;
      WaveFunction2 += Tmp;
    }
  for (int i = 2; i < this->InvFillingFactor; ++i)
    WaveFunction2 *= WaveFunction;
  return WaveFunction2;
}
