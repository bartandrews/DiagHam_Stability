////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of Pfaffian wave function with two quasiholes on disk          //
//                                                                            //
//                        last modification : 23/10/2008                      //
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
#include "Tools/FQHEWaveFunction/PfaffianOnDiskTwoQuasiholeWaveFunction.h"
#include "Vector/RealVector.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"


// constructor
//
// nbrParticles = number of particles
// zHole1 = position of the first quasihole
// zHole2 = position of the second quasihole (spherical coordinates, theta angle)
// fermions = flag indicating whether to calculate bosonic or fermionic pfaffian

PfaffianOnDiskTwoQuasiholeWaveFunction::PfaffianOnDiskTwoQuasiholeWaveFunction(int nbrParticles, Complex zHole1, Complex zHole2, bool fermions)
{
  this->NbrParticles = nbrParticles;
  this->ZHole1 = zHole1;
  this->ZHole2 = zHole2;
  this->TmpPfaffian = ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->FermionFlag = fermions;
}



// copy constructor
//
// function = reference on the wave function to copy

PfaffianOnDiskTwoQuasiholeWaveFunction::PfaffianOnDiskTwoQuasiholeWaveFunction(const PfaffianOnDiskTwoQuasiholeWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->ZHole1 = function.ZHole1;
  this->ZHole2 = function.ZHole2;
  this->TmpPfaffian = ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->FermionFlag=function.FermionFlag;
}

// destructor
//

PfaffianOnDiskTwoQuasiholeWaveFunction::~PfaffianOnDiskTwoQuasiholeWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PfaffianOnDiskTwoQuasiholeWaveFunction::Clone ()
{
  return new PfaffianOnDiskTwoQuasiholeWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PfaffianOnDiskTwoQuasiholeWaveFunction::operator ()(RealVector& x)
{
  Complex TmpZ;
  Complex TmpHole1;
  Complex TmpHole2;
  Complex TmpHole3;
  Complex TmpHole4;
  ComplexSkewSymmetricMatrix TmpPfaffian (this->NbrParticles);
  Complex WaveFunction(1.0);
  Complex Tmp;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      TmpZ.Re = x[i << 1];
      TmpZ.Im = x[1 + (i << 1)];
      TmpHole1 = TmpZ - this->ZHole1;
      TmpHole2 = TmpZ - this->ZHole2;
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp.Re = x[j << 1];
	  Tmp.Im = x[1 + (j << 1)];
	  TmpHole3 = Tmp - this->ZHole1;
	  TmpHole4 = Tmp - this->ZHole2;
	  Tmp -= TmpZ;
	  Tmp *= -1.0;
	  WaveFunction *= Tmp;	  
	  Tmp = (TmpHole1 *TmpHole4  + TmpHole2 * TmpHole3) / Tmp;
	  TmpPfaffian.SetMatrixElement (i , j, Tmp);
	}
    }
  if (this->FermionFlag == true)
    WaveFunction *= WaveFunction;
  WaveFunction *= TmpPfaffian.Pfaffian();
  return WaveFunction;
}

