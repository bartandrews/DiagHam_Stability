////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             Class of polynomial with long rational coefficients            //
//                                                                            //
//                        last modification : 17/11/2010                      //
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


#ifndef LONGRATIONALPOLYNOMIAL_H
#define LONGRATIONALPOLYNOMIAL_H


#include "config.h"
#include "MathTools/LongRational.h"
#include "MathTools/Complex.h"

#include <iostream>


using std::ostream;


class LongRationalPolynomial
{

private:

  // array where polynomial coefficients are stored
  LongRational* Coefficient;
  

  // polynomial degree
  int Degree;
  // maximal polynomial degree that can be used without reallocating memory
  int MaxDegree;

  // number of found roots
  int NbrRoot;
  // flag to indicate of roots have been evaluated
  bool RootFlag;
  // array where roots are stored
  Complex* Root;
  // array where rational roots are stored (only available if all roots are known)
  LongRational* RationalRoots;

public:

  // default constructor
  //
  LongRationalPolynomial ();
 
  // constructor
  //
  // degree = polynomial degree
  LongRationalPolynomial (int degree);

  // constructor from raw datas
  //
  // degree = polynomial degree
  // coefficients = coefficients array ( first element is associated to the -power term)
  // flag = true if coefficients array has to be used directly and not duplicated
  LongRationalPolynomial (int degree, LongRational* coefficients, bool flag = true);

  // constructor from raw datas, including all rational roots
  //
  // degree = polynomial degree
  // coefficients = coefficients array ( first element is associated to the -power term)
  // roots = array where the roots are stored
  // flag = true if coefficient and root arrays have to be used directly and not duplicated
  LongRationalPolynomial (int degree, LongRational* coefficients, LongRational* roots, bool flag = true);

  // copy constructor
  //
  // P = polynomial to copy
  LongRationalPolynomial (const LongRationalPolynomial& P);
  
  // constructor from P1 * P2
  //
  // P1 = first polynomial
  // P2 = second polynomial
  LongRationalPolynomial (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2);

  // constructor from P1+q^n P2
  //
  // P1 = first polynomial
  // P2 = second polynomial
  // degree = polynomial degree of the monomial with which P2 has to be multiplied 
  LongRationalPolynomial (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2, int degree);

  // destructor
  //
  ~LongRationalPolynomial();

  // assignement
  //
  // P = polynomial to copy
  // return value = reference on current polynomial 
  LongRationalPolynomial& operator = (const LongRationalPolynomial& P);

  // test if the polynomial has been defined
  // 
  // return value = true if the polynomial hass been defined
  bool Defined();

  // Return polynomial value at a given point
  //
  // x = point where to evaluate polynomial
  // return value = polynomial value at x
  LongRational PolynomialEvaluate (const LongRational& x);
  
  // Return polynomial value at a given point
  //
  // x = point where to evaluate polynomial
  // result = variable where to store the result
  // return value = reference on the current polynomial
  LongRationalPolynomial& PolynomialEvaluate (const LongRational& x, LongRational& result);

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
  LongRational DerivativeEvaluate (LongRational x);

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
  LongRationalPolynomial DerivatePolynomial ();
  
  // Evaluate polynomial n-th derivative
  //
  // n = derivative order
  // return value = polynomial n-th derivative
  LongRationalPolynomial DerivatePolynomial (int n);

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
  LongRational& operator [] (int n);

  // return coefficient corresponding to the nth-degree
  //
  // n = index of the corresponding coefficient
  // return value = coefficient
  LongRational PolynomialCoefficient (int n);

  // Return the n-th root (with root(n-1) <= root(n))
  //
  // n = index of the root (0 being root with the lowest modulus)
  // return value = n-th root
  Complex PolynomialRoot (int n);

  // refine root value 
  void RefineRoot (double Epsilon, int MaxIter);

  // Return the n-th raional root 
  //
  // n = index of the root 
  // return value = n-th root
  inline LongRational PolynomialRationalRoot (int n);

  // arithmetic operators
  friend LongRationalPolynomial operator - (const LongRationalPolynomial& P);
  friend LongRationalPolynomial operator + (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2);
  friend LongRationalPolynomial operator - (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2);
  friend LongRationalPolynomial operator * (const LongRationalPolynomial& P, const LongRational& d);
  friend LongRationalPolynomial operator * (const LongRational& d, const LongRationalPolynomial& P);
  friend LongRationalPolynomial operator * (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2);
  friend LongRationalPolynomial operator / (const LongRationalPolynomial& P, const double& d);
  friend LongRationalPolynomial operator / (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2); // Degree P1 > Degree P2
  friend LongRationalPolynomial operator % (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2); // Degree P1 > Degree P2
  LongRationalPolynomial& operator += (const LongRationalPolynomial& P);
  LongRationalPolynomial& operator -= (const LongRationalPolynomial& P);
  LongRationalPolynomial& operator *= (const LongRational& d);
  LongRationalPolynomial& operator *= (long d);
  LongRationalPolynomial& operator *= (const LongRationalPolynomial& P);
  LongRationalPolynomial& operator /= (const LongRational& d);
  LongRationalPolynomial& operator /= (long d);
  //  LongRationalPolynomial&MultiplyAndAdd (const LongRationalPolynomial& P1 , const LongRationalPolynomial& P2);

  // Divide polynomial by a monomial (z - z0) using Horner scheme
  //
  // z0 = monomial root
  // return value = result of polynomial division
  LongRationalPolynomial MonomialDivision (const LongRational& z0);


  // Divide the current polynomial by a monomial (z - z0) using Horner scheme
  //
  // z0 = monomial root
  // return value = reference on the current polynomial
  LongRationalPolynomial& LocalMonomialDivision (const LongRational& z0);

  // Multiply polynomial by a monomial (z - z0)
  //
  // z0 = monomial root
  // return value = result of polynomial multiplication
  LongRationalPolynomial MonomialMultiplication (const LongRational& z0);

  // Multiply the current polynomial by a monomial (z - z0)
  //
  // z0 = monomial root
  // return value = reference on the current polynomial
  LongRationalPolynomial& LocalMonomialMultiplication (const LongRational& z0);

  // shift all powers from a given value
  //
  // shift = shift to apply
  // return value = reference on the current polynomial
  LongRationalPolynomial& ShiftPowers(int shift);

  //Output Stream overload
  friend ostream& operator << (ostream& Str, const LongRationalPolynomial& P);

  
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

inline int LongRationalPolynomial::GetNbrRoots ()
{
  return this->NbrRoot;
}
  
// get the degree of the polynomial
//
// return value = degree of the polynomial

inline int LongRationalPolynomial::GetPolynomialDegree () 
{
  return this->Degree;
}
  
// return coefficient corresponding to the nth-degree
//
// n = index of the corresponding coefficient
// return value = coefficient

inline LongRational LongRationalPolynomial::PolynomialCoefficient (int n)
{
  return this->Coefficient[n];
}

// Return the n-th root (with root(n-1) <= root(n))
//
// n = index of the root (0 being root with the lowest modulus)
// return value = n-th root

inline Complex LongRationalPolynomial::PolynomialRoot (int n) 
{
  return this->Root[n];
}

// Return the n-th raional root 
//
// n = index of the root 
// return value = n-th root

inline LongRational LongRationalPolynomial::PolynomialRationalRoot (int n) 
{
  return this->RationalRoots[n];
}

// return reference on the coefficient corresponding to the nth-degree
//
// n = index of the corresponding coefficient
// return value = reference on the coefficient

inline LongRational& LongRationalPolynomial::operator [] (int n)
{
  return this->Coefficient[n];
}

// test if the polynomial has been defined
// 
// return value = true if the polynomial hass been defined

inline bool LongRationalPolynomial::Defined()
{
  return !(this->Coefficient == 0);
}

// Return polynomial value at a given point
//
// x = point where to evaluate polynomial
// result = variable where to store the result
// return value = reference on the current polynomial

inline LongRationalPolynomial& LongRationalPolynomial::PolynomialEvaluate (const LongRational& x, LongRational& result)
{
  if (this->Coefficient != 0)
    {
      result =  this->Coefficient[this->Degree];
      for (int i = this->Degree - 1 ; i >= 0; i--)
	{
	  result *= x;
	  result += this->Coefficient[i];
	}
    }
  return *this;
}

#endif

