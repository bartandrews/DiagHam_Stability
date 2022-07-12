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


#ifndef BINOMIALCOEFFICIENTS_H
#define BINOMIALCOEFFICIENTS_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class BinomialCoefficients
{

 private:

  // maximum major index that can be requested
  int MaximumIndex;

  // array that contained all binomial coefficients (integer version)
  long** Coefficients;

  // array that contained all binomial coefficients (double version)
  double** NumericalCoefficients;

   // garbage flag used for the polynomials
  GarbageFlag Flag;

 public:

  // default constructor
  //
  // maximumIndex = new maximum major index (aka the (m) upper index in the [m n] notation)
  BinomialCoefficients (int maximumIndex = 2);

  // copy constructor 
  // 
  // binomial = binomial coefficients to copy
  BinomialCoefficients (const BinomialCoefficients& binomial);

  // destructor
  //
  ~BinomialCoefficients();

  // request a given binomial coefficient (referring to the [m n] notation)
  //
  // m = major index
  // n = minor index
  // return value = reference on the correponding binomial coefficient
  long operator () (int m, int n);

  // request a given binomial coefficient (referring to the [m n] notation) in double format
  //
  // m = major index
  // n = minor index
  // return value = reference on the correponding binomial coefficient
  double GetNumericalCoefficient (int m, int n);

  // modify the maximum major index that can be requested
  //
  // maximumIndex = new maximum major index
  // return value = reference on current BinomialCoefficients
  BinomialCoefficients& Resize(int maximumIndex);

 private:


};

// request a given binomial coefficient (referring to the [m n] notation)
//
// m = major index
// n = minor index
// return value = reference on the correponding binomial coefficient 

inline long BinomialCoefficients::operator ()(int m, int n)
{
  if ((m < 0) || (n < 0) || (m < n))
    return 0l;
  if (m > this->MaximumIndex)
    this->Resize(m);
  if ((n <= (m >> 1)))
    {
      return this->Coefficients[m][n];
    }
  else
    {
      return this->Coefficients[m][m - n];
    }
}

// request a given binomial coefficient (referring to the [m n] notation) in double format
//
// m = major index
// n = minor index
// return value = reference on the correponding binomial coefficient

inline double BinomialCoefficients::GetNumericalCoefficient (int m, int n)
{
  if ((m < 0) || (n < 0) || (m < n))
    return 0.0;
  if (m > this->MaximumIndex)
    this->Resize(m);
  if ((n <= (m >> 1)))
    {
      return this->NumericalCoefficients[m][n];
    }
  else
    {
      return this->NumericalCoefficients[m][m - n];
    }
}

#endif

