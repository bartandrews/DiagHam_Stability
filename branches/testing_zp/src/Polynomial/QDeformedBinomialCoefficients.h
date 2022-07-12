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


#ifndef QDEFORMEDBINOMIALCOEFFICIENTS_H
#define QDEFORMEDBINOMIALCOEFFICIENTS_H


#include "config.h"
#include "Polynomial/IntegerPolynomial.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class QDeformedBinomialCoefficients
{

 private:

  // polynomial used when an invalid q-deformed binomial coefficient is requested
  IntegerPolynomial* NullPolynomial;

  // maximum major index that can be requested
  int MaximumIndex;

  // array that contained all q-deformed binomial coefficients
  IntegerPolynomial*** Coefficients;

   // garbage flag used for the polynomials
  GarbageFlag Flag;

 public:

  // default constructor
  //
  // maximumIndex = new maximum major index (aka the (m) upper index in the [m n] notation)
  QDeformedBinomialCoefficients (int maximumIndex = 2);

  // copy constructor 
  // 
  // binomial = q-deformed binomial coefficients to copy
  QDeformedBinomialCoefficients (const QDeformedBinomialCoefficients& binomial);

  // destructor
  //
  ~QDeformedBinomialCoefficients();

  // request a given q-deformed binomial coefficient (referring to the [m n] notation)
  //
  // m = major index
  // n = minor index
  // return value = reference on the correponding q-deformed binomial coefficient (MUST NOT BE MODIFIED!)
  IntegerPolynomial& operator () (int m, int n);

  // modify the maximum major index that can be requested
  //
  // maximumIndex = new maximum major index
  // return value = reference on current QDeformedBinomialCoefficients
  QDeformedBinomialCoefficients& Resize(int maximumIndex);

 private:


};

// request a given q-deformed binomial coefficient (referring to the [m n] notation)
//
// m = major index
// n = minor index
// return value = reference on the correponding q-deformed binomial coefficient (MUST NOT BE MODIFIED!)

inline IntegerPolynomial& QDeformedBinomialCoefficients::operator ()(int m, int n)
{
  if ((m < 0) || (n < 0) || (m < n))
    return *(this->NullPolynomial);
  if (m > this->MaximumIndex)
    this->Resize(m);
  if ((n <= (m >> 1)))
    {
      return *(this->Coefficients[m][n]);
    }
  else
    {
      return *(this->Coefficients[m][m - n]);
    }
}

#endif

