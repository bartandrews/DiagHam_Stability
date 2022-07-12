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


#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H


#include "config.h"
#include <iostream>
#include "Complex.h"


#define POLYERROR 0.0001
#define MC1_3 0.333333333333
#define MC1_9 0.111111111111
#define MC1_27 0.037037037037
#define cos2PI_3 -0.5
#define sin2PI_3 0.866025403785


using std::ostream;


class Polynomial
{

private:

  double* Coefficient;
  double* Root;
  int Degree;
  int NbrRoot;
  bool RootFlag;
  bool NoRootFlag;

public:
 
  // constructors
  Polynomial (int Deg, double* Coef);
  Polynomial (int Deg, double* Coef, bool flag); // flag = true if Coef has to be used directly
  Polynomial (const Polynomial& P);
  
  // destructor
  ~Polynomial();

  // assignement
  Polynomial& operator = (const Polynomial& P);

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
  Polynomial DerivatePolynomial ();
  
  // Evaluate polynomial n-th derivative
  //
  // n = derivative order
  // return value = polynomial n-th derivative
  Polynomial DerivatePolynomial (int n);

  // Find Roots of the polynomial
  void SolvePolynomial();
  
  // Return Number of root
  int NumberRoot () {return this->NbrRoot;};
  
  // Return Degree of the polynomial
  int PolynomialDegree () {return this->Degree;};
  
  // Return Coefficient corresponding to the nth-degree
  double PolynomialCoefficient (int n) {return this->Coefficient[n];};

  // Return the n-th root (with root(n-1) <= root(n))
  double PolynomialRoot (int n) {return this->Root[n];};

  //Refine root value 
  void RefineRoot (double Epsilon, int MaxIter);

  // arithmetic operators
  friend Polynomial operator - (const Polynomial& P);
  friend Polynomial operator + (const Polynomial& P1, const Polynomial& P2);
  friend Polynomial operator - (const Polynomial& P1, const Polynomial& P2);
  friend Polynomial operator * (const Polynomial& P, const double& d);
  friend Polynomial operator * (const double& d, const Polynomial& P);
  friend Polynomial operator * (const Polynomial& P1, const Polynomial& P2);
  friend Polynomial operator / (const Polynomial& P, const double& d);
  friend Polynomial operator / (const Polynomial& P1, const Polynomial& P2); // Degree P1 > Degree P2
  friend Polynomial operator % (const Polynomial& P1, const Polynomial& P2); // Degree P1 > Degree P2
  Polynomial& operator += (const Polynomial& P);
  Polynomial& operator -= (const Polynomial& P);
  Polynomial& operator *= (const double& d);
  Polynomial& operator *= (const Polynomial& P);
  Polynomial& operator /= (const double& d);
  Polynomial& operator /= (const Polynomial& P);

  // Divide polynomial by a monomial (x - x0) using Horner scheme
  //
  // x0 = monomial root
  // return value = result of polynomial division
  Polynomial MonomialDivision (double x0);

  //Output Stream overload
  friend ostream& operator << (ostream& Str, const Polynomial& P);
  
private:

  // Find Root of a linear polynomial
  void SolveLinear();

  // Find Roots of a quadratic polynomial
  void SolveQuadratic();

  // Find Roots of a cubic polynomial
  void SolveCubic();

  // Find Roots of a quartic polynomial
  void SolveQuartic();

  // Find Roots of a polynomial of degree > 4
  void SolveHigherOrder();

  //Find a root given local minimun and local maximum (with Xmin < Xmax)
  //
  // xmax = interval upper bound
  // xmin = interval lower bound
  // err = relative error on the root
  // absErr = absolute error on the root
  // return value = root belonging to the interval [xmin, xmax]
  double FindRootBottomUp (double xmin, double xmax, double err, double absErr);

  //Find a root given local minimun and local maximum (with Xmin > Xmax)
  //
  // xmax = interval upper bound
  // xmin = interval lower bound
  // err = relative error on the root
  // absErr = absolute error on the root
  // return value = root belonging to the interval [xmin, xmax]
  double FindRootTopDown (double xmax, double xmin, double err, double absErr);

  public:

  // use Laguerre Method to find roots
  //
  // err = relative error on the root
  // absErr = absolute error on the root
  // nbrIteration = maximum iteration number to find a root
  void LaguerreMethodSolver (double err, double absErr, int nbrIteration);

 private:

  // use Bauhuber Method to find a root (put NoRootFlag to true if no root was find)
  //
  // nbrIteration = maximum iteration number to find a root
  // err = relative error on the root
  // absErr = absolute error on the root
  // return value = root obtained
  Complex BauhuberMethod (double err, double absErr, int nbrIteration, Complex x0);

  // use Laguerre Method to find a root (put NoRootFlag to true if no root was find)
  //
  // nbrIteration = maximum iteration number to find a root
  // err = relative error on the root
  // absErr = absolute error on the root
  // return value = root obtained
  double LaguerreMethod (double err, double absErr, int nbrIteration, double x0 = 0.0);

  // use Laguerre Method to find a root (put NoRootFlag to true if no root was find)
  //
  // nbrIteration = maximum iteration number to find a root
  // err = relative error on the root
  // absErr = absolute error on the root
  // return value = root obtained
  Complex LaguerreMethod (double err, double absErr, int nbrIteration, Complex x0);

  // Sort roots of a polynomial
  //
  void SortRoots ();

};

#endif

