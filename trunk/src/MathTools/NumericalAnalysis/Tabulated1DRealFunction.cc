////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of tabulated 1D real function                    //
//                                                                            //
//                        last modification : 08/07/2004                      //
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
#include "MathTools/NumericalAnalysis/Tabulated1DRealFunction.h"

#include <iostream>


using std::cout;
using std::endl;


// constructor
//
// tabulatedCoordinates = array containing the coordinates where the tabluated function is defined
// tabulatedValues = array containing function tabulated values 
// nbrValues = number of tabluated values

Tabulated1DRealFunction::Tabulated1DRealFunction(double* tabulatedCoordinates, double* tabulatedValues, long nbrValues)
{
  this->TabulatedValues = tabulatedValues;
  this->TabulatedCoordinates = tabulatedCoordinates;
  this->NbrValues = nbrValues;
  this->Flag.Initialize();
}

// copy constructor 
//
// function = function to copy

Tabulated1DRealFunction::Tabulated1DRealFunction(const Tabulated1DRealFunction& function)
{
  this->TabulatedValues = function.TabulatedValues;
  this->TabulatedCoordinates = function.TabulatedCoordinates;
  this->NbrValues = function.NbrValues;
  this->Flag = function.Flag;
}

// destructor
//

Tabulated1DRealFunction::~Tabulated1DRealFunction()
{ 
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->TabulatedValues;
      delete[] this->TabulatedCoordinates;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DRealFunction* Tabulated1DRealFunction::Clone ()
{
  return new Tabulated1DRealFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

double Tabulated1DRealFunction::operator ()(double x)
{
  return this->ValueFromLinearInterpolation(x);
}

// get function derivative at a given point
//
// x = point where the function derivative has to be evaluated
// return value = function derivative

double Tabulated1DRealFunction::GetDerivative(const double& x)
{
  return 0.0;
}

// evaluate derivative of the function on the same interval that the function itself
//
// return value = function derivative

Abstract1DRealFunction* Tabulated1DRealFunction::GetDerivative()
{
  return 0;
}

// get function laplacian at a given point
//
// x = point where the function laplacian has to be evaluated
// return value = function laplacian

double Tabulated1DRealFunction::GetLaplacian(const double& x)
{
  return 0.0;
}

// evaluate laplacian of the function on the same interval that the function itself
//
// return value = function laplacian

Abstract1DRealFunction* Tabulated1DRealFunction::GetLaplacian()
{
  return 0;
}

// evaluate integral on the function of a given interval
// 
// interval = reference on the interval on which the integral has to be evaluated
// return value = integral value

double Tabulated1DRealFunction::GetIntegral(AbstractNumericalInterval& interval)
{
  return 0.0;
}

// evaluate the primitive of the function on the same interval that the function itself
//
// return value = function primitive

Abstract1DRealFunction* Tabulated1DRealFunction::GetPrimitive()
{
  double* TmpTabulatedValues = new double[this->NbrValues + 1];
  double* TmpTabulatedCoordinates = new double[this->NbrValues + 1];
  TmpTabulatedValues[0] = 0.0;
  TmpTabulatedCoordinates[0] = this->TabulatedCoordinates[0];
  for (int i = 1; i < this->NbrValues; ++i)
    {
      TmpTabulatedCoordinates[i] = 0.5 * (this->TabulatedCoordinates[i - 1] + this->TabulatedCoordinates[i]);
      TmpTabulatedValues[i] = (TmpTabulatedValues[i - 1] 
			       + (this->TabulatedValues[i - 1] * (this->TabulatedCoordinates[i] - this->TabulatedCoordinates[i - 1])));
    }
  TmpTabulatedCoordinates[this->NbrValues] =  this->TabulatedCoordinates[this->NbrValues - 1];
  TmpTabulatedValues[this->NbrValues] = 1.0;
  return new Tabulated1DRealFunction(TmpTabulatedCoordinates, TmpTabulatedValues, this->NbrValues + 1);
}

