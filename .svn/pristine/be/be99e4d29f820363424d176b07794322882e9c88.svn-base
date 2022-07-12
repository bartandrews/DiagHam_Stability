////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of clebsch gordan coefficients                   //
//                                                                            //
//                        last modification : 19/06/2002                      //
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


#ifndef JACOBIPOLYNOMIALS_H
#define JACOBIPOLYNOMIALS_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"
#include <iostream>


using std::ostream;


class JacobiPolynomials
{
 private:

  // conjugation parameters alpha, beta
  double ParameterA;
  
  double ParameterB;

  // maximal degree
  int MaxDegreeN;

  // recursion coefficients
  double *CNPlusOne;
  double *CNConst;
  double *CNLinear;
  double *CNMinusOne;

  // coefficients of the first polynomial
  double A1Zero;
  double A1One;

  
  // array holding the function values of the degree n functions for the last argument
  double LastArgument;
  double *FunctionValues;

  // table of prefactors for explicit expansion
  double **ExplicitExpansionPrefactors;

  // (optional) last degree requested
  int LastN;
  
 public:

  // default constructor
  // creates object that knows about P0 and P1 with alpha=beta=0
  //
  JacobiPolynomials();

  // constructor for a general theta function \theta[^a_b](z|tau)
  // maxDegreeN = maximum degree of the function
  // param_a = parameter a
  // param_b = parameter b
  // tau = modulus
  // precision = required precision
  JacobiPolynomials(int maxDegreeN, double a, double b);

  // copy constructor
  JacobiPolynomials (const JacobiPolynomials& p);

  // destructor
  //
  ~JacobiPolynomials ();

  
  // get the value of the function for a given coordinate z
  //
  // n = degree (must be <=MaxDegreeN)
  // x = argument
  // return = function value of P_n(x)
  double GetValue(int n, double x);

  // get the value of the function from a direct series expansion rather than recursively
  // n = degree (must be <=MaxDegreeN)
  // x = argument
  // return = function value of P_n(x)
  double GetExplicitFunctionValue(int n, double x);

  // get the value of the function for a given coordinate z
  //
  // x = argument
  // return = function value of P_n(x)
  double* GetValues(double x);

  
  // pretty-print a function values
  // str = stream to print to
  // z = point where to evaluate
  ostream& PrintValues(ostream &str, double x);

 private:
  // evaluate recursively up to degree n
  void RunRecursion(int &n, double x);

  // initialize recursion coefficients
  void InitializeRecursion();

};

#endif
