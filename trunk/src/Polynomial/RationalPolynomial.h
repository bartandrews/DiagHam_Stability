////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  Class of polynomial with rational coefficients            //
//                                                                            //
//                        last modification : 13/11/2010                      //
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


#ifndef RATIONALPOLYNOMIAL_H
#define RATIONALPOLYNOMIAL_H


#include "config.h"
#include "MathTools/Rational.h"
#include "MathTools/Complex.h"

#include <iostream>


using std::ostream;


class RationalPolynomial
{

private:

  // array where polynomial coefficients are stored
  Rational* Coefficient;
  

  // polynomial degree
  int Degree;

  // number of found roots
  int NbrRoot;
  // flag to indicate of roots have been evaluated
  bool RootFlag;
  // array where roots are stored
  Complex* Root;
  // array where rational roots are stored (only available if all roots are known)
  Rational* RationalRoots;

public:

  // default constructor
  //
  RationalPolynomial ();
 
  // constructor
  //
  // degree = polynomial degree
  RationalPolynomial (int degree);

  // constructor from raw datas
  //
  // degree = polynomial degree
  // coefficients = coefficients array ( first element is associated to the -power term)
  // flag = true if coefficients array has to be used directly and not duplicated
  RationalPolynomial (int degree, Rational* coefficients, bool flag = true);

  // constructor from raw datas, including all rational roots
  //
  // degree = polynomial degree
  // coefficients = coefficients array ( first element is associated to the -power term)
  // roots = array where the roots are stored
  // flag = true if coefficient and root arrays have to be used directly and not duplicated
  RationalPolynomial (int degree, Rational* coefficients, Rational* roots, bool flag = true);

  // copy constructor
  //
  // P = polynomial to copy
  RationalPolynomial (const RationalPolynomial& P);
  
  // constructor from P1+q^n P2
  //
  // P1 = first polynomial
  // P2 = second polynomial
  // degree = polynomial degree of the monomial with which P2 has to be multiplied 
  RationalPolynomial (const RationalPolynomial& P1, const RationalPolynomial& P2, int degree);

  // destructor
  //
  ~RationalPolynomial();

  // assignement
  //
  // P = polynomial to copy
  // return value = reference on current polynomial 
  RationalPolynomial& operator = (const RationalPolynomial& P);

  // test if the polynomial has been defined
  // 
  // return value = true if the polynomial hass been defined
  bool Defined();

  // Return polynomial value at a given point
  //
  // x = point where to evaluate polynomial
  // return value = polynomial value at x
  Rational PolynomialEvaluate (const Rational& x);
  
  // Return polynomial value at a given point
  //
  // x = point where to evaluate polynomial
  // return value = polynomial value at x
  double PolynomialEvaluate (double x);
  
  // Return polynomial value at a given point
  //
  // x = point where to evaluate polynomial
  // return value = polynomial value at x
  Complex PolynomialEvaluate (Complex x);
  
  // Evaluate polynomial derivative   
  //
  // x = position where to evaluate polynomial derivative
  // return value = polynomial derivative at x
  Rational DerivativeEvaluate (Rational x);

  // Evaluate polynomial derivative   
  //
  // x = position where to evaluate polynomial derivative
  // return value = polynomial derivative at x
  double DerivativeEvaluate (double x);

  // Evaluate polynomial derivative   
  //
  // x = position where to evaluate polynomial derivative
  // return value = polynomial derivative at x
  Complex DerivativeEvaluate (Complex x);

  // Evaluate polynomial n-th derivative 
  //
  // x = position where to evaluate polynomial n-th derivative
  // n = derivative order
  // return value = polynomial n-th derivative at x  
  double DerivativeEvaluate (double x, int n);

  // Evaluate polynomial n-th derivative 
  //
  // x = position where to evaluate polynomial n-th derivative
  // n = derivative order
  // return value = polynomial n-th derivative at x  
  Complex DerivativeEvaluate (Complex x, int n);

  // Return Derivative of the polynomial 
  //
  // return value = polynomial derivative
  RationalPolynomial DerivatePolynomial ();
  
  // Evaluate polynomial n-th derivative
  //
  // n = derivative order
  // return value = polynomial n-th derivative
  RationalPolynomial DerivatePolynomial (int n);

  // Find Roots of the polynomial
  void SolvePolynomial();
  
  // get the number of roots
  //
  // return value = number of roots
  int GetNbrRoots ();
  
  // get the degree of the polynomial
  //
  // return value = degree of the polynomial
  int GetPolynomialDegree ();
  
  // return reference on the coefficient corresponding to the nth-degree
  //
  // n = index of the corresponding coefficient
  // return value = reference on the coefficient
  Rational& operator [] (int n);

  // return coefficient corresponding to the nth-degree
  //
  // n = index of the corresponding coefficient
  // return value = coefficient
  Rational PolynomialCoefficient (int n);

  // Return the n-th raional root 
  //
  // n = index of the root 
  // return value = n-th root
  inline Rational PolynomialRationalRoot (int n);

  // Return the n-th root (with root(n-1) <= root(n))
  //
  // n = index of the root (0 being root with the lowest modulus)
  // return value = n-th root
  Complex PolynomialRoot (int n);

  // refine root value 
  void RefineRoot (double Epsilon, int MaxIter);

  // arithmetic operators
  friend RationalPolynomial operator - (const RationalPolynomial& P);
  friend RationalPolynomial operator + (const RationalPolynomial& P1, const RationalPolynomial& P2);
  friend RationalPolynomial operator - (const RationalPolynomial& P1, const RationalPolynomial& P2);
  friend RationalPolynomial operator * (const RationalPolynomial& P, const Rational& d);
  friend RationalPolynomial operator * (const Rational& d, const RationalPolynomial& P);
  friend RationalPolynomial operator * (const RationalPolynomial& P1, const RationalPolynomial& P2);
  friend RationalPolynomial operator / (const RationalPolynomial& P, const double& d);
  friend RationalPolynomial operator / (const RationalPolynomial& P1, const RationalPolynomial& P2); // Degree P1 > Degree P2
  friend RationalPolynomial operator % (const RationalPolynomial& P1, const RationalPolynomial& P2); // Degree P1 > Degree P2
  RationalPolynomial& operator += (const RationalPolynomial& P);
  RationalPolynomial& operator -= (const RationalPolynomial& P);
  RationalPolynomial& operator *= (const Rational& d);
  RationalPolynomial& operator *= (long d);
  RationalPolynomial& operator *= (const RationalPolynomial& P);
  RationalPolynomial& operator /= (const Rational& d);
  RationalPolynomial& operator /= (long d);

  // Divide polynomial by a monomial (z - z0) using Horner scheme
  //
  // z0 = monomial root
  // return value = result of polynomial division
  RationalPolynomial MonomialDivision (const Rational& z0);

  // Multiply polynomial by a monomial (z - z0)
  //
  // z0 = monomial root
  // return value = result of polynomial multiplication
  RationalPolynomial MonomialMultiplication (const Rational& z0);

  // shift all powers from a given value
  //
  // shift = shift to apply
  // return value = reference on the current polynomial
  RationalPolynomial& ShiftPowers(int shift);

  //Output Stream overload
  friend ostream& operator << (ostream& Str, const RationalPolynomial& P);
  
private:

  // Find Root of a linear polynomial
  void SolveLinear();

  // Sort roots of a polynomial
  //
  void SortRoots ();

};

// return number of roots
//
// return value = number of roots

inline int RationalPolynomial::GetNbrRoots ()
{
  return this->NbrRoot;
}
  
// get the degree of the polynomial
//
// return value = degree of the polynomial

inline int RationalPolynomial::GetPolynomialDegree () 
{
  return this->Degree;
}
  
// return coefficient corresponding to the nth-degree
//
// n = index of the corresponding coefficient
// return value = coefficient

inline Rational RationalPolynomial::PolynomialCoefficient (int n)
{
  return this->Coefficient[n];
}

// Return the n-th root (with root(n-1) <= root(n))
//
// n = index of the root (0 being root with the lowest modulus)
// return value = n-th root

inline Complex RationalPolynomial::PolynomialRoot (int n) 
{
  return this->Root[n];
}

// Return the n-th raional root 
//
// n = index of the root 
// return value = n-th root

inline Rational RationalPolynomial::PolynomialRationalRoot (int n) 
{
  return this->RationalRoots[n];
}

// return reference on the coefficient corresponding to the nth-degree
//
// n = index of the corresponding coefficient
// return value = reference on the coefficient

inline Rational& RationalPolynomial::operator [] (int n)
{
  return this->Coefficient[n];
}

// test if the polynomial has been defined
// 
// return value = true if the polynomial hass been defined

inline bool RationalPolynomial::Defined()
{
  return !(this->Coefficient == 0);
}

#endif

