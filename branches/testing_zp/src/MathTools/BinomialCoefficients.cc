////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        Class of binomial coefficients                      //
//                                                                            //
//                        last modification : 01/06/2005                      //
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
#include "MathTools/BinomialCoefficients.h"


// default constructor
//
// maximumIndex = new maximum major index (aka the (m) upper index in the [m n] notation)

BinomialCoefficients::BinomialCoefficients (int maximumIndex)
{
  if (maximumIndex < 2)
    maximumIndex = 2;
  this->MaximumIndex = 2;
  this->Coefficients = new long* [3];
  this->Coefficients[0] = new long [1];
  this->Coefficients[0][0] = 1l;
  this->Coefficients[1] = new long [1];
  this->Coefficients[1][0] = 1l;
  this->Coefficients[2] = new long [2];
  this->Coefficients[2][0] = 1l;
  this->Coefficients[2][1] = 2l;  
  this->NumericalCoefficients = new double* [3];
  this->NumericalCoefficients[0] = new double [1];
  this->NumericalCoefficients[0][0] = 1.0;
  this->NumericalCoefficients[1] = new double [1];
  this->NumericalCoefficients[1][0] = 1.0;
  this->NumericalCoefficients[2] = new double [2];
  this->NumericalCoefficients[2][0] = 1.0;
  this->NumericalCoefficients[2][1] = 2.0;
  this->Resize(maximumIndex);
  this->Flag.Initialize();
}

// copy constructor 
// 
// binomial = q-deformed binomial coefficients to copy

BinomialCoefficients::BinomialCoefficients (const BinomialCoefficients& binomial)
{
  this->MaximumIndex = binomial.MaximumIndex;
  this->Coefficients = new long* [this->MaximumIndex + 1];
  this->NumericalCoefficients = new double* [this->MaximumIndex + 1];
  for (int i = 0; i <= this->MaximumIndex; ++i)
    {
      this->Coefficients[i] = binomial.Coefficients[i];
      this->NumericalCoefficients[i] = binomial.NumericalCoefficients[i];
    }
  this->Flag = binomial.Flag;
}

// destructor
//

BinomialCoefficients::~BinomialCoefficients()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      for (int i = 0; i <= this->MaximumIndex; ++i)
	{
	  delete[] this->Coefficients[i];
          delete[] this->NumericalCoefficients[i];
	}
    }
  delete[] this->Coefficients; 
  delete[] this->NumericalCoefficients;
}

// modify the maximum major index that can be requested
//
// maximumIndex = new maximum major index
// return value = reference on current BinomialCoefficients

BinomialCoefficients& BinomialCoefficients::Resize(int maximumIndex)
{
  if (maximumIndex <= this->MaximumIndex)
    return *this;
  long** TmpCoefficients = new long* [maximumIndex + 1];
  double** TmpNumericalCoefficients = new double* [maximumIndex + 1];
  int i = 0;
  if (this->Flag.Shared() == false)
    for (; i <= this->MaximumIndex; ++i)
      {
	TmpCoefficients[i] = this->Coefficients[i];
	TmpNumericalCoefficients[i] = this->NumericalCoefficients[i];
      }
  else
    {
      TmpCoefficients[0] = new long [1];
      this->Coefficients[0][0] = 0l;
      TmpNumericalCoefficients[0] = new double [1];
      this->NumericalCoefficients[0][0] = 0.0;
      for (; i <= this->MaximumIndex; ++i)
	{	  
	  int Half = i >> 1;
	  TmpCoefficients[i] = new long [Half + 1];
          TmpCoefficients[i][0] = TmpCoefficients[0][0];
	  TmpNumericalCoefficients[i] = new double [Half + 1];
	  TmpNumericalCoefficients[i][0] = TmpNumericalCoefficients[0][0];
	  for (int j = 1; j <= Half; ++j)
	    {
	      TmpCoefficients[i][j] = this->Coefficients[i][j];     
              TmpNumericalCoefficients[i][j] = this->NumericalCoefficients[i][j];
	    }
	}
      this->Flag.Initialize();
    }
  for (; i <= maximumIndex; ++i)
    {
      int Half = i >> 1;
      TmpCoefficients[i] = new long [Half + 1];
      TmpCoefficients[i][0] = TmpCoefficients[0][0];
      TmpNumericalCoefficients[i] = new double [Half + 1];
      TmpNumericalCoefficients[i][0] = TmpNumericalCoefficients[0][0];
      if ((i & 1) != 0)
	for (int j = 1; j <= Half; ++j)
	  {
	    TmpCoefficients[i][j] = TmpCoefficients[i - 1][j - 1] + TmpCoefficients[i - 1][j];
	    TmpNumericalCoefficients[i][j] = TmpNumericalCoefficients[i - 1][j - 1] + TmpNumericalCoefficients[i - 1][j];
	  }
      else
	{
	  for (int j = 1; j < Half; ++j)
	    {
              TmpCoefficients[i][j] = TmpCoefficients[i - 1][j - 1] + TmpCoefficients[i - 1][j];
 	      TmpNumericalCoefficients[i][j] = TmpNumericalCoefficients[i - 1][j - 1] + TmpNumericalCoefficients[i - 1][j];
	    }
	  TmpCoefficients[i][Half] = TmpCoefficients[i - 1][Half - 1] + TmpCoefficients[i - 1][i - 1 - Half];
          TmpNumericalCoefficients[i][Half] = TmpNumericalCoefficients[i - 1][Half - 1] + TmpNumericalCoefficients[i - 1][i - 1 - Half];
	}
    }
  delete[] this->Coefficients;
  this->Coefficients = TmpCoefficients;
  delete[] this->NumericalCoefficients;
  this->NumericalCoefficients = TmpNumericalCoefficients;
  this->MaximumIndex = maximumIndex;
  return *this;
}
