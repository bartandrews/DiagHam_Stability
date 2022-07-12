////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of function basis for particle on sphere             //
//                                                                            //
//                        last modification : 10/12/2002                      //
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
#include "FunctionBasis/ParticleOnCylinderFunctionBasis.h"
#include "Vector/RealVector.h"

#include <math.h>

#include <iostream>
using std::cout;
using std::endl;

// constructor
//
// maxMomentum = maximum momentum reached by a particle
// landauLevel = Landau level index
// ratio = aspect ratio of the cylinder
// indexShiftFlag = true if apply a (maxMomentum / 2) shift to state index

ParticleOnCylinderFunctionBasis::ParticleOnCylinderFunctionBasis(int maxMomentum, int landauLevel, double ratio, bool indexShiftFlag)
{
  this->MaxMomentum = maxMomentum;
  this->LandauLevel = landauLevel;
  this->Ratio = ratio;
  this->HilbertSpaceDimension = this->MaxMomentum + 1;
  this->Perimeter = sqrt(2.0 * M_PI * (this->MaxMomentum + 1.0) * this->Ratio);
  this->Kappa = 2.0 * M_PI / this->Perimeter;
  this->Normalization = 1.0 / sqrt(this->Perimeter * sqrt(M_PI));
  if (this->LandauLevel == 1)
    {
      this->Normalization /= sqrt(2.0);
    }
  else 
    {
      if (this->LandauLevel > 1)
	{
	  cout << "LL >= 2 " << endl;
	  exit(1);
	}
    }
  if (indexShiftFlag)
      this->IndexShift = 0.5 * ((double) this->MaxMomentum);
  else
      this->IndexShift = 0;
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// x, y = coordinates where the function should be evaluated
// index = the function index 
// returns value of the wavefunction

Complex ParticleOnCylinderFunctionBasis::GetFunctionValue(double x, double y, double index)
{  
  Complex Phase;
  Phase.Re = cos(this->Kappa * index * y);
  Phase.Im = sin(this->Kappa * index * y);

  Complex Result = Phase * exp(-0.5 * pow(x - this->Kappa * index, 2.0));

  if (this->LandauLevel == 1)
    {
      Result *= (2.0 * (x - this->Kappa * index));    
    }
 
  return (Result * this->Normalization);
}

// get value of the i-th function at a given point (for functions which take values in C)
//
// x, y = coordinates where the function should be evaluated
// index = the function index 
// returns value of the wavefunction

Complex ParticleOnCylinderFunctionBasis::GetFunctionValue(double x, double y, int index)
{  
  double Tmp = this->Kappa * (((double) index) - this->IndexShift); 
  Complex Result = Phase(Tmp * y) * exp(-0.5 * pow(x - Tmp, 2.0));
  if (this->LandauLevel == 1)
    {
      Result *= (2.0 * (x - Tmp));    
    } 
  return (Result * this->Normalization);
}
