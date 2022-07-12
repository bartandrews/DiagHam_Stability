////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                           RayTrace version  0.10                           //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   Class of polynomial with real coefficients               //
//                                                                            //
//                        last modification : 15/01/2001                      //
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


#ifndef SPECIALPOLYNOMIAL_H
#define SPECIALPOLYNOMIAL_H


#include "config.h"
#include "MathTools/Complex.h"
#include "Polynomial.h"
#include "Vector/RealVector.h"

#ifdef __GMP__
#include <gmp.h>
#endif


// return a Laguerre Polynomial of rank n (optionally associated Laguerre function)
// n = index
// alpha = modifier
Polynomial LaguerrePolynomial(int n, int alpha=0, bool verbose=false);



class LaguerrePolynomialRecursion
{
 public:
  // constructor
  // maximum order to be calculated
  LaguerrePolynomialRecursion(int nmax);

  // destructor
  ~LaguerrePolynomialRecursion(){}

  // return the value for a Laguerre Polynomial of rank n
  // n = index
  double EvaluateLaguerrePolynomial(int n, double x);

  // evaluate all Laguerre polynomials at the given coordinate
  void EvaluateLaguerrePolynomials(RealVector &, double x);

#ifdef __GMP__
  void CheckAgainstGMP(double xDouble, mp_bitcnt_t precision);
#endif

 protected:
  // temporary array for values to be calculated
  RealVector TmpValues;

  // maximum order
  int NMax;
};


#ifdef __GMP__

class LaguerrePolynomialRecursionAP
{
 public:
  // constructor
  // maximum order to be calculated
  LaguerrePolynomialRecursionAP(int nmax, mp_bitcnt_t precision);

  // destructor
  ~LaguerrePolynomialRecursionAP();

  // return the value for a Laguerre Polynomial of rank n
  // n = index
  void EvaluateLaguerrePolynomial(int N, mpf_t &result, const mpf_t &x);

  // evaluate all Laguerre polynomials at the given coordinate
  void EvaluateLaguerrePolynomials(int N, mpf_t *output, mpf_t &x);


 protected:
  // temporary array for values to be calculated
  mpf_t *TmpValues;
  mpf_t TmpValue;
  mpf_t TmpValue2;

  // maximum order
  int NMax;
};

#endif


#endif

