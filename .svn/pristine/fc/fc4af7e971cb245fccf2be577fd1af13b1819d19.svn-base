////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of Laughlin wave function with one quasihole on sphere        //
//                                                                            //
//                        last modification : 30/10/2008                      //
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
#include "Tools/FQHEWaveFunction/FQHESphereLaughlinOneQuasiholeWaveFunction.h"
#include "Vector/RealVector.h"
#include "GeneralTools/Endian.h"


#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;



// constructor
//
// nbrParticles = number of particles
// theta = position of the quasihole (spherical coordinates, theta angle)
// phi = position of the quasihole (spherical coordinates, phi angle)
// fermions = flag indicating whether to calculate bosonic or fermionic laughlin wave function

FQHESphereLaughlinOneQuasiholeWaveFunction::FQHESphereLaughlinOneQuasiholeWaveFunction(int nbrParticles, double theta, double phi, bool fermions)
{
  this->NbrParticles = nbrParticles;
  this->FermionFlag = fermions;

  this->UHole.Re = cos(0.5*phi);
  this->UHole.Im= -sin(0.5*phi);
  this->UHole *= cos(0.5*theta);

  this->VHole.Re = cos(0.5*phi);
  this->VHole.Im = sin(0.5*phi);
  this->VHole *= sin(0.5*theta);
}

// copy constructor
//
// function = reference on the wave function to copy

FQHESphereLaughlinOneQuasiholeWaveFunction::FQHESphereLaughlinOneQuasiholeWaveFunction(const FQHESphereLaughlinOneQuasiholeWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;

  this->UHole = function.UHole;
  this->VHole = function.VHole;

  this->FermionFlag=function.FermionFlag;
}

// destructor
//

FQHESphereLaughlinOneQuasiholeWaveFunction::~FQHESphereLaughlinOneQuasiholeWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* FQHESphereLaughlinOneQuasiholeWaveFunction::Clone ()
{
  return new FQHESphereLaughlinOneQuasiholeWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex FQHESphereLaughlinOneQuasiholeWaveFunction::operator ()(RealVector& x)
{
  Complex WaveFunction(1.0);
  return WaveFunction;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)

Complex FQHESphereLaughlinOneQuasiholeWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  Complex TmpU1;
  Complex TmpV1;
  Complex WaveFunction(1.0);
  Complex WaveFunction2(1.0);
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      TmpU1 = uv[i << 1];
      TmpV1 = uv[1 + (i << 1)];
      WaveFunction *= ((this->UHole * TmpV1) - (this->VHole * TmpU1));
      for (int j = i + 1; j < this->NbrParticles; ++j)
	WaveFunction2 *= (TmpU1 * uv[1 + (j << 1)]) - (uv[j << 1] * TmpV1);
    }

  if (this->FermionFlag == false)
    WaveFunction2 *= WaveFunction2;
  else
    WaveFunction2 *= WaveFunction2 * WaveFunction2;
  WaveFunction *= WaveFunction2;
  return WaveFunction;
}

