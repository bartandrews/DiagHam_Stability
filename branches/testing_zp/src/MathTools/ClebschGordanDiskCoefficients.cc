////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of the equivalent of clebsch gordan coefficients          //
//          for addition of two angular monenta for Landau states in the LLL  //
//                              and on the disk geometry                      //
//                                                                            //
//                        last modification : 05/07/2008                      //
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
#include "MathTools/ClebschGordanDiskCoefficients.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"

#include <stdlib.h>
#include <math.h>
#include <iostream>


using std::cout;
using std::endl;


// default constructor
//

ClebschGordanDiskCoefficients::ClebschGordanDiskCoefficients()
{
  this->MaximumMomentum = -1;
  this->Coefficients = 0l;
}


// constructor 
//
// mmax = maximum angular momentum for a single particle

ClebschGordanDiskCoefficients::ClebschGordanDiskCoefficients(int mmax)
{
  this->MaximumMomentum = mmax;
  this->Flag.Initialize();
  this->EvaluateClebschGordanDiskCoefficients();
}

// copy constructor (without duplicating datas)
//
// coefficients = reference on Clebsch Gordan coefficients to copy

ClebschGordanDiskCoefficients::ClebschGordanDiskCoefficients (const ClebschGordanDiskCoefficients& coefficients)
{
  this->MaximumMomentum = coefficients.MaximumMomentum;
  this->Coefficients = coefficients.Coefficients;
  this->Flag = coefficients.Flag;
}

// destructor
//

ClebschGordanDiskCoefficients::~ClebschGordanDiskCoefficients ()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i <= this->MaximumMomentum; ++i)
	{
	  for (int j = 0; j <=  this->MaximumMomentum; ++j)
	    delete[] this->Coefficients[i][j];
	  delete[] this->Coefficients[i];
	}
      delete[] this->Coefficients;
    }
}

// assignment (without duplicating datas)
//
// coefficients = reference on Clebsch Gordan coefficients to assign
// return value = reference on current Clebsch Gordan coefficients

ClebschGordanDiskCoefficients& ClebschGordanDiskCoefficients::operator = (const ClebschGordanDiskCoefficients& coefficients)
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i <= this->MaximumMomentum; ++i)
	{
	  for (int j = 0; j <= this->MaximumMomentum; ++j)
	    delete[] this->Coefficients[i][j];
	  delete[] this->Coefficients[i];
	}
      delete[] this->Coefficients;
    }
  this->Flag = coefficients.Flag;
  this->MaximumMomentum = coefficients.MaximumMomentum;
  this->Coefficients = coefficients.Coefficients;
  return *this;
}

// get a particular coefficient (without testing if m1, m2 and j are valid)
//
// m1 = projection of first angular momentum 
// m2 = projection of second angular momentum 
// j = resulting angular momentum
// return value = corresponding Clebsch Gordan coefficient

double ClebschGordanDiskCoefficients::GetCoefficient (int m1, int m2, int j)
{
  return this->Coefficients[m1][m2][j];
}

// print a particular coefficient (without testing if m1, m2 and j are valid)
//
// str = reference on output stream
// m1 = projection of first angular momentum 
// m2 = projection of second angular momentum 
// j = resulting angular momentum
// return value = reference on output stream

ostream& ClebschGordanDiskCoefficients::PrintCoefficient (ostream& str, int m1, int m2, int j)
{
  str << " < "  << m1 << " " << m2 << " | " << j << " > = " << this->Coefficients[m1][m2][j];
  return str;
}

// evaluate all Clebsch Gordan coefficients using Schulten Gordon recursion algorithm
//

void ClebschGordanDiskCoefficients::EvaluateClebschGordanDiskCoefficients()
{
  int TotalM = (this->MaximumMomentum + 1);
  this->Coefficients = new double** [TotalM];
  for (int i = 0; i < TotalM; ++i)
    {
      this->Coefficients[i] = new double* [TotalM];
      for (int j = 0; j < TotalM; ++j)
	this->Coefficients[i][j] = new double [i + j + 1];
    }
  BinomialCoefficients BCoefficients (this->MaximumMomentum);
  FactorialCoefficient TmpCoef;

  for (int i = 0; i <= this->MaximumMomentum; ++i)
    for (int j = 0; j <= i; ++j)
      {
	int TmpMax = i + j;
	TmpCoef.SetToOne();
	TmpCoef.PartialFactorialMultiply(TmpMax - i + 1, TmpMax);
	TmpCoef.FactorialDivide(i);
	TmpCoef.Power2Divide(TmpMax);
	double TmpCoefficient = sqrt(TmpCoef.GetNumericalValue() / (4.0 * M_PI));       
	for (int k = 0; k <= TmpMax; ++k)
	  {
	    double TmpCoefficient2 = 0.0;
	    double Sign = 1.0;
	    int TmpMin2 = 0;
	    int TmpMax2 = i;
	    if (k < i)
	      TmpMax2 = k;
	    if (k > j)
	      TmpMin2 = k - j;
	    if (((k - TmpMin2) & 1) != 0)
	      Sign = -1.0;
	    for (int l = TmpMin2; l <= TmpMax2; ++l)
	      {
		TmpCoefficient2 += BCoefficients.GetNumericalCoefficient(i, l) * BCoefficients.GetNumericalCoefficient(j, k - l) * Sign;
		Sign *= -1.0;
	      }
	    TmpCoefficient2 *= TmpCoefficient;
	    this->Coefficients[i][j][k] = TmpCoefficient2;
	    if (j != i)
	      {
		Sign = 1.0;
		if ((k & 1) != 0)
		  Sign = -1.0;
		this->Coefficients[j][i][k] = Sign * TmpCoefficient2;
	      }
	    if (k != TmpMax)
	      TmpCoefficient /= sqrt(0.5 * ((double) (TmpMax - k)));	    
	  }	
      }

}

