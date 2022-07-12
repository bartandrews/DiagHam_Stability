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


#ifndef TABULATED1DREALFUNCTION_H
#define TABULATED1DREALFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Tabulated1DRealFunction.h"


// constructor from raw datas
//
// interval = interval on which the function is defined (including the subdivision scheme)

Tabulated1DRealFunction::Tabulated1DRealFunction(AbstractSubdividedInterval* interval, double* tabulatedValues)
{
  this->AbstractSubdividedInterval = AbstractSubdividedInterval->Clone();
  if (tabulatedValues == 0)
    {
      this->TabulatedValues = new double [this->AbstractSubdividedInterval->GetNbrSubdivision()];
    }
  else
    {
      this->TabulatedValues = tabulatedValues;
    }
}


// constructor from a C function
//
// interval = interval on which the function is defined (including the subdivision scheme)
// function = pointer to the C function that described the mathematical function that will be used to initialized the tabulated function

Tabulated1DRealFunction::Tabulated1DRealFunction(AbstractSubdividedInterval* interval, double (*function) (double))
{
}

// copy constructor 
//
// function = function to copy

Tabulated1DRealFunction::Tabulated1DRealFunction(const Tabulated1DRealFunction& function)
{
  this->AbstractSubdividedInterval = AbstractSubdividedInterval->Clone();
  if (tabulatedValues == 0)
    {
      this->TabulatedValues = new double [this->AbstractSubdividedInterval->GetNbrSubdivision()];
    }
  else
    {
      this->TabulatedValues = tabulatedValues;
    }
}


// destructor
//

Tabulated1DRealFunction::~Tabulated1DRealFunction()
{
  this->AbstractSubdividedInterval = AbstractSubdividedInterval->Clone();
  this->TabulatedValues = new double [this->AbstractSubdividedInterval->GetNbrSubdivision()];
  for (int i = 0; i < this->AbstractSubdividedInterval->GetNbrSubdivision(); ++i)
    {
      this->TabulatedValues[i] = function(x);
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
}

// get function derivative at a given point
//
// x = point where the function derivative has to be evaluated
// return value = function derivative

double Tabulated1DRealFunction::GetDerivative(const double& x)
{
}

// evaluate derivative of the function on the same interval that the function itself
//
// return value = function derivative

Abstract1DRealFunction* Tabulated1DRealFunction::GetDerivative()
{
  long Lim = this->AbstractSubdividedInterval->GetNbrSubdivision() - 1;
  double* TmpTabulatedValues = new double [this->AbstractSubdividedInterval->GetNbrSubdivision()];
  if ((this->AbstractSubdividedInterval->GetType() & AbstractNumericalInterval::RegularSubdivision) != 0)
    {
      double Coefficient = 0.5 / this->AbstractSubdividedInterval->GetSubdivisionStep(); 
      for (long i = 1; i < Lim; ++i)
	{
	  TmpTabulatedValues[i] = (this->TabulatedValues[i + 1] - this->TabulatedValues[i - 1]) * Coefficient;      
	}
      TmpTabulatedValues[0] = 2.0 * (this->TabulatedValues[1] - this->TabulatedValues[0]) * Coefficient;
      TmpTabulatedValues[Lim] = 2.0 * (this->TabulatedValues[Lim] - this->TabulatedValues[Lim - 1]) * Coefficient;
    }
  else
    {
      for (long i = 1; i < Lim; ++i)
	{
	  TmpTabulatedValues[i] = ((this->TabulatedValues[i + 1] - this->TabulatedValues[i - 1]) /
				   (this->AbstractSubdividedInterval.GetSubdivisionStep(i) + 
				    this->AbstractSubdividedInterval.GetSubdivisionStep(i - 1)));      
	}
      TmpTabulatedValues[0] = ((this->TabulatedValues[1] - this->TabulatedValues[0]) /
			       this->AbstractSubdividedInterval.GetSubdivisionStep(0));      
      TmpTabulatedValues[Lim] = ((this->TabulatedValues[Lim] - this->TabulatedValues[Lim - 1]) /
				 this->AbstractSubdividedInterval.GetSubdivisionStep(Lim - 1));      
    }
  return new Tabulated1DRealFunction(this->AbstractSubdividedInterval, TmpTabulatedValues);
}

// get function laplacian at a given point
//
// x = point where the function laplacian has to be evaluated
// return value = function laplacian

double Tabulated1DRealFunction::GetLaplacian(const double& x)
{
}

// evaluate laplacian of the function on the same interval that the function itself
//
// return value = function laplacian

Abstract1DRealFunction* Tabulated1DRealFunction::GetLaplacian()
{
  long Lim = this->AbstractSubdividedInterval->GetNbrSubdivision() - 1;
  double* TmpTabulatedValues = new double [this->AbstractSubdividedInterval->GetNbrSubdivision()];
  if ((this->AbstractSubdividedInterval->GetType() & AbstractNumericalInterval::RegularSubdivision) != 0)
    {
      double Coefficient = 1.0 / (this->AbstractSubdividedInterval->GetSubdivisionStep() * this->AbstractSubdividedInterval->GetSubdivisionStep()); 
      for (long i = 1; i < Lim; ++i)
	{
	  TmpTabulatedValues[i] = (this->TabulatedValues[i + 1] + this->TabulatedValues[i - 1] - (2.0 * this->TabulatedValues[i])) * Coefficient;      
	}
      TmpTabulatedValues[0] = ((this->TabulatedValues[0] + this->TabulatedValues[2]) - 2.0 * this->TabulatedValues[1]) * Coefficient;
      TmpTabulatedValues[Lim] = ((this->TabulatedValues[Lim] + this->TabulatedValues[Lim - 2]) - 2.0 * this->TabulatedValues[Lim - 1]) * Coefficient;
    }
  else
    {
      for (long i = 1; i < Lim; ++i)
	{
	  TmpTabulatedValues[i] = (((this->AbstractSubdividedInterval.GetSubdivisionStep(i) * (this->TabulatedValues[i + 1] - this->TabulatedValues[i]))
				    + (this->AbstractSubdividedInterval.GetSubdivisionStep(i + 1) * (this->TabulatedValues[i - 1] - this->TabulatedValues[i]))) / 
				   ((this->AbstractSubdividedInterval.GetSubdivisionStep(i) + this->AbstractSubdividedInterval.GetSubdivisionStep(i - 1)) *
				    this->AbstractSubdividedInterval.GetSubdivisionStep(i) * this->AbstractSubdividedInterval.GetSubdivisionStep(i - 1)));      
	}
      TmpTabulatedValues[0] = (((this->AbstractSubdividedInterval.GetSubdivisionStep(1) * (this->TabulatedValues[1] - this->TabulatedValues[0]))
				- (this->AbstractSubdividedInterval.GetSubdivisionStep(0) * (this->TabulatedValues[2] - this->TabulatedValues[0]))) / 
			       ((this->AbstractSubdividedInterval.GetSubdivisionStep(0) - this->AbstractSubdividedInterval.GetSubdivisionStep(1)) *
				this->AbstractSubdividedInterval.GetSubdivisionStep(0) * this->AbstractSubdividedInterval.GetSubdivisionStep(1)));      
      TmpTabulatedValues[Lim] = (((this->AbstractSubdividedInterval.GetSubdivisionStep(Lim - 2) * (this->TabulatedValues[Lim - 1] - this->TabulatedValues[Lim]))
				  - (this->AbstractSubdividedInterval.GetSubdivisionStep(Lim - 1) * (this->TabulatedValues[Lim - 2] - this->TabulatedValues[Lim]))) / 
				 ((this->AbstractSubdividedInterval.GetSubdivisionStep(Lim - 1) - this->AbstractSubdividedInterval.GetSubdivisionStep(Lim - 2)) *
				  this->AbstractSubdividedInterval.GetSubdivisionStep(Lim - 1) * this->AbstractSubdividedInterval.GetSubdivisionStep(Lim - 2)));      
    }
  return new Tabulated1DRealFunction(this->AbstractSubdividedInterval, TmpTabulatedValues);
}

// evaluate integral on the function of a given interval
// 
// interval = reference on the interval on which the integral has to be evaluated
// return value = integral value

double Tabulated1DRealFunction::GetIntegral(AbstractNumericalInterval& interval)
{
}

// use linear interpolation to obtain function value at a given point
//
// x = point where the function has to be evaluated
// return value = function value obatined using linear interpolation

inline double Tabulated1DRealFunction::ValueFromLinearInterpolation (const double& x)
{
  long Index = this->AbstractSubdividedInterval->GetFloorIndex(x);
  if (Index == (this->AbstractSubdividedInterval->GetNbrSubdivision() - 1))
    return this->TabulatedValues[Index];
  return ((this->TabulatedValues[Index + 1] - this->TabulatedValues[Index]) * (x - this->AbstractSubdividedInterval->GetPosistionfromIndex(Index)) /
	  this->AbstractSubdividedInterval.GetSubdivisionStep(Index)) + this->TabulatedValues[Index];
}
  
// use spline interpolation to obtain function value at a given point
//
// x = point where the function has to be evaluated
// return value = function value obatined using spline interpolation

inline double Tabulated1DRealFunction::ValueFromSplineInterpolation (const double& x)
{
  long Index = this->AbstractSubdividedInterval->GetFloorIndex(x);
  if (Index == 0)
    return this->TabulatedValues[Index];
  if (Index == (this->AbstractSubdividedInterval->GetNbrSubdivision() - 2))
    return this->TabulatedValues[Index];
  if (Index == (this->AbstractSubdividedInterval->GetNbrSubdivision() - 1))
    return this->TabulatedValues[Index];
  return ((this->TabulatedValues[Index + 1] - this->TabulatedValues[Index]) * (x - this->AbstractSubdividedInterval->GetPosistionfromIndex(Index)) /
	  this->AbstractSubdividedInterval.GetSubdivisionStep(Index)) + this->TabulatedValues[Index];
}
  

