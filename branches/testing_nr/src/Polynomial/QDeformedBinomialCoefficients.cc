////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   Class of q-deformed binomial coefficients                //
//                                                                            //
//                        last modification : 28/05/2005                      //
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
#include "Polynomial/QDeformedBinomialCoefficients.h"


// default constructor
//
// maximumIndex = new maximum major index (aka the (m) upper index in the [m n] notation)

QDeformedBinomialCoefficients::QDeformedBinomialCoefficients (int maximumIndex)
{
  if (maximumIndex < 2)
    maximumIndex = 2;
  this->MaximumIndex = 2;
  long Constant = 0l;
  this->NullPolynomial = new IntegerPolynomial(0, &Constant, false);
  this->Coefficients = new IntegerPolynomial** [3];
  this->Coefficients[0] = new IntegerPolynomial* [1];
  Constant = 1l;
  this->Coefficients[0][0] = new IntegerPolynomial(0, &Constant, false);
  this->Coefficients[1] = new IntegerPolynomial* [1];
  this->Coefficients[1][0] = this->Coefficients[0][0];
  this->Coefficients[2] = new IntegerPolynomial* [2];
  this->Coefficients[2][0] = this->Coefficients[0][0];
  this->Coefficients[2][1] = new IntegerPolynomial(*(this->Coefficients[0][0]), *(this->Coefficients[0][0]), 1);  
  this->Resize(maximumIndex);
  this->Flag.Initialize();
}

// copy constructor 
// 
// binomial = q-deformed binomial coefficients to copy

QDeformedBinomialCoefficients::QDeformedBinomialCoefficients (const QDeformedBinomialCoefficients& binomial)
{
  this->NullPolynomial = binomial.NullPolynomial;
  this->MaximumIndex = binomial.MaximumIndex;
  this->Coefficients = new IntegerPolynomial** [this->MaximumIndex + 1];
  for (int i = 0; i <= this->MaximumIndex; ++i)
    this->Coefficients[i] = binomial.Coefficients[i];
  this->Flag = binomial.Flag;
}

// destructor
//

QDeformedBinomialCoefficients::~QDeformedBinomialCoefficients()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete this->Coefficients[0][0];
      delete[] this->Coefficients[0];
      for (int i = 1; i <= this->MaximumIndex; ++i)
	{
	  int Half = i >> 1;
	  for (int j = 1; j <= Half; ++j)
	    delete this->Coefficients[i][j];
	  delete[] this->Coefficients[i];
	}
      delete this->NullPolynomial;
    }
  delete[] this->Coefficients;
}

// modify the maximum major index that can be requested
//
// maximumIndex = new maximum major index
// return value = reference on current QDeformedBinomialCoefficients

QDeformedBinomialCoefficients& QDeformedBinomialCoefficients::Resize(int maximumIndex)
{
  if (maximumIndex <= this->MaximumIndex)
    return *this;
  IntegerPolynomial*** TmpCoefficients = new IntegerPolynomial** [maximumIndex + 1];
  int i = 0;
  if (this->Flag.Shared() == false)
    for (; i <= this->MaximumIndex; ++i)
      TmpCoefficients[i] = this->Coefficients[i];
  else
    {
      long Constant = 1l;
      TmpCoefficients[0] = new IntegerPolynomial* [1];
      this->Coefficients[0][0] = new IntegerPolynomial(0, &Constant, false);
      for (; i <= this->MaximumIndex; ++i)
	{	  
	  int Half = i >> 1;
	  TmpCoefficients[i] = new IntegerPolynomial* [Half + 1];
	  TmpCoefficients[i][0] = TmpCoefficients[0][0];
	  for (int j = 1; j <= Half; ++j)
	    TmpCoefficients[i][j] = this->Coefficients[i][j];     
	}
      this->Flag.Initialize();
    }
  for (; i <= maximumIndex; ++i)
    {
      int Half = i >> 1;
      TmpCoefficients[i] = new IntegerPolynomial* [Half + 1];
      TmpCoefficients[i][0] = TmpCoefficients[0][0];
      if ((i & 1) != 0)
	for (int j = 1; j <= Half; ++j)
	  {
	    TmpCoefficients[i][j] = new IntegerPolynomial(*(TmpCoefficients[i - 1][j - 1]), *(TmpCoefficients[i - 1][j]), j);
	  }
      else
	{
	  for (int j = 1; j < Half; ++j)
	    {
	      TmpCoefficients[i][j] = new IntegerPolynomial(*(TmpCoefficients[i - 1][j - 1]), *(TmpCoefficients[i - 1][j]), j);
	    }
	  TmpCoefficients[i][Half] = new IntegerPolynomial(*(TmpCoefficients[i - 1][Half - 1]), *(TmpCoefficients[i - 1][i - 1 - Half]), Half);
	}
    }
  delete[] this->Coefficients;
  this->Coefficients = TmpCoefficients;
  this->MaximumIndex = maximumIndex;
  return *this;
}
