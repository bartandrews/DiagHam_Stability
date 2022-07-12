////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                           RayTrace version  0.10                           //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   Class of polynomial with complex coefficients            //
//                                                                            //
//                        last modification : 23/03/2005                      //
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


#ifndef COMPLEXPOLYNOMIAL_H
#define COMPLEXPOLYNOMIAL_H


#include "config.h"
#include <iostream>
#include "MathTools/Complex.h"


using std::ostream;


class ComplexPolynomial
{

private:

  // array where polynomial coefficients are stored
  Complex* Coefficient;
  
  // array where roots are stored
  Complex* Root;

  // polynomial degree
  int Degree;

  // number of found roots
  int NbrRoot;

  bool RootFlag;
  bool NoRootFlag;

public:

  // default constructor
  //
  ComplexPolynomial ();
 
  // constructor
  //
  // degree = polynomial degree
  // coefficients = coefficients array ( first element is associated to the -power term)
  // flag = true if coefficients array has to be used directly and not duplicated
  ComplexPolynomial (int degree, Complex* coefficients, bool flag = true);

  // copy constructor
  //
  // P = polynomial to copy
  ComplexPolynomial (const ComplexPolynomial& P);
  
  // destructor
  //
  ~ComplexPolynomial();

  // assignement
  //
  // P = polynomial to copy
  // return value = reference on current polynomial 
  ComplexPolynomial& operator = (const ComplexPolynomial& P);

  // Return polynomial value at a given point
  //
  // x = point where to evaluate polynomial
  // return value = polynomial value at x
  Complex PolynomialEvaluate (double x);
  
  // Return polynomial value at a given point
  //
  // x = point where to evaluate polynomial
  // return value = polynomial value at x
  Complex PolynomialEvaluate (Complex x);
  
  // Evaluate polynomial derivative   
  //
  // x = position where to evaluate polynomial derivative
  // return value = polynomial derivative at x
  Complex DerivativeEvaluate (double x);

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
  Complex DerivativeEvaluate (double x, int n);

  // Evaluate polynomial n-th derivative 
  //
  // x = position where to evaluate polynomial n-th derivative
  // n = derivative order
  // return value = polynomial n-th derivative at x  
  Complex DerivativeEvaluate (Complex x, int n);

  // Return Derivative of the polynomial 
  //
  // return value = polynomial derivative
  ComplexPolynomial DerivatePolynomial ();
  
  // Evaluate polynomial n-th derivative
  //
  // n = derivative order
  // return value = polynomial n-th derivative
  ComplexPolynomial DerivatePolynomial (int n);

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
  
  // return coefficient corresponding to the nth-degree
  //
  // n = index of the corresponding coefficient
  // return value = coefficient
  Complex PolynomialCoefficient (int n);

  // Return the n-th root (with root(n-1) <= root(n))
  //
  // n = index of the root (0 being root with the lowest modulus)
  // return value = n-th root
  Complex PolynomialRoot (int n);

  // refine root value 
  void RefineRoot (double Epsilon, int MaxIter);

  // arithmetic operators
  friend ComplexPolynomial operator - (const ComplexPolynomial& P);
  friend ComplexPolynomial operator + (const ComplexPolynomial& P1, const ComplexPolynomial& P2);
  friend ComplexPolynomial operator - (const ComplexPolynomial& P1, const ComplexPolynomial& P2);
  friend ComplexPolynomial operator * (const ComplexPolynomial& P, const double& d);
  friend ComplexPolynomial operator * (const double& d, const ComplexPolynomial& P);
  friend ComplexPolynomial operator * (const ComplexPolynomial& P1, const ComplexPolynomial& P2);
  friend ComplexPolynomial operator / (const ComplexPolynomial& P, const double& d);
  friend ComplexPolynomial operator / (const ComplexPolynomial& P1, const ComplexPolynomial& P2); // Degree P1 > Degree P2
  friend ComplexPolynomial operator % (const ComplexPolynomial& P1, const ComplexPolynomial& P2); // Degree P1 > Degree P2
  ComplexPolynomial& operator += (const ComplexPolynomial& P);
  ComplexPolynomial& operator -= (const ComplexPolynomial& P);
  ComplexPolynomial& operator *= (const double& d);
  ComplexPolynomial& operator *= (const ComplexPolynomial& P);
  ComplexPolynomial& operator /= (const double& d);
  ComplexPolynomial& operator /= (const ComplexPolynomial& P);

  // Divide polynomial by a monomial (z - z0) using Horner scheme
  //
  // z0 = monomial root
  // return value = result of polynomial division
  ComplexPolynomial MonomialDivision (const Complex& z0);

  //Output Stream overload
  friend ostream& operator << (ostream& Str, const ComplexPolynomial& P);
  
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

inline int ComplexPolynomial::GetNbrRoots ()
{
  return this->NbrRoot;
}
  
// get the degree of the polynomial
//
// return value = degree of the polynomial

inline int ComplexPolynomial::GetPolynomialDegree () 
{
  return this->Degree;
}
  
// return coefficient corresponding to the nth-degree
//
// n = index of the corresponding coefficient
// return value = coefficient

inline Complex ComplexPolynomial::PolynomialCoefficient (int n)
{
  return this->Coefficient[n];
}

// Return the n-th root (with root(n-1) <= root(n))
//
// n = index of the root (0 being root with the lowest modulus)
// return value = n-th root

inline Complex ComplexPolynomial::PolynomialRoot (int n) 
{
  return this->Root[n];
}

#endif

