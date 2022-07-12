////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of abstract 1D complex function                   //
//                                                                            //
//                        last modification : 01/09/2004                      //
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
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunctionOnSphere.h"


// virtual destructor
//

Abstract1DComplexFunctionOnSphere::~Abstract1DComplexFunctionOnSphere()
{
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex Abstract1DComplexFunctionOnSphere::operator ()(RealVector& x)
{
  double Factor = M_PI * 0.5;
  int Max = x.GetVectorDimension();
  ComplexVector UVCoordinates (Max);    
  for (int i = 0; i < Max; ++i)
    {
      UVCoordinates.Re(i) = cos(0.5 * x[i << 1]);
      UVCoordinates.Im(i) = UVCoordinates.Re(i) * sin(0.5 * x[1 + (i << 1)]) * Factor;
      UVCoordinates.Re(i) *= cos(0.5 * x[1 + (i << 1)]) * Factor;
      ++i;
      UVCoordinates.Re(i) = sin(0.5 * x[i << 1]);
      UVCoordinates.Im(i) = - UVCoordinates.Re(i) * sin(0.5 * x[1 + (i << 1)]) * Factor;
      UVCoordinates.Re(i) *= cos(0.5 * x[1 + (i << 1)]) * Factor;
    }
  Complex Tmp = this->CalculateFromSpinorVariables(UVCoordinates);
  return Tmp;
}

