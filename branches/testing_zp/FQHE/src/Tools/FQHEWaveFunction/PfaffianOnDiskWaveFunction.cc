////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of Pfaffian wave function on disk                  //
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
#include "Tools/FQHEWaveFunction/PfaffianOnDiskWaveFunction.h"
#include "Vector/RealVector.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"


// constructor
//
// nbrParticles = number of particles

PfaffianOnDiskWaveFunction::PfaffianOnDiskWaveFunction(int nbrParticles)
{
  this->NbrParticles = nbrParticles;
}

// copy constructor
//
// function = reference on the wave function to copy

PfaffianOnDiskWaveFunction::PfaffianOnDiskWaveFunction(const PfaffianOnDiskWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
}

// destructor
//

PfaffianOnDiskWaveFunction::~PfaffianOnDiskWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PfaffianOnDiskWaveFunction::Clone ()
{
  return new PfaffianOnDiskWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PfaffianOnDiskWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp;
  ComplexSkewSymmetricMatrix TmpPfaffian (this->NbrParticles);
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
	  WaveFunction *= Tmp;
	  Tmp = 1.0 / Tmp;
	  TmpPfaffian.SetMatrixElement (i , j, Tmp);
	}
   }
  WaveFunction *= TmpPfaffian.Pfaffian();
  return WaveFunction;
}
