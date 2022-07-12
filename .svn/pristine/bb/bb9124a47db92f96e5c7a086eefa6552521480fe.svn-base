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


#ifndef MONOPOLEHARMONICS_H
#define MONOPOLEHARMONICS_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"
#include "MathTools/JacobiPolynomials.h"
#include <iostream>


using std::ostream;


class MonopoleHarmonics
{

 private:

  // conjugation parameters alpha, beta
  int MomentumTwoM;

  // twice the monopole charge
  int ChargeTwoQ;

  // maximal degree
  int MaxOrderN;

  // recursion coefficients
  double *NormPrefactor;

  // Jacobi polynomial
  JacobiPolynomials *Jacobi;

  // offset of order of Jacobi polynomials w.r.t index of monopole harmonics
  int JacobiOffset;

  // table with the functino values
  double *FunctionValues;
  
  // array holding the function values of the degree n functions for the last argument
  double LastArgument;

  // (optional) last degree requested
  int LastN;
  
 public:

  // default constructor
  // creates object that knows about P0 and P1 with alpha=beta=0
  //
  MonopoleHarmonics();

  // constructor for a sequence of monopole harmonics Y_q,n,m with varying n (q,m fixed)
  // maxOrderN = maximum order of the function
  // twoQ = twice the monopole charge
  // twoM = twice the Lz component of angular momentum
  MonopoleHarmonics(int twoQ, int MaxOrderN, int twoM);

  // copy constructor
  MonopoleHarmonics (const MonopoleHarmonics& p);

  // destructor
  //
  ~MonopoleHarmonics ();

  
  // get the value of the function for a given coordinate z
  //
  // n = degree (must be <=MaxDegreeN)
  // x = argument
  // return = function value of P_n(x)
  double GetValue(int n, double x);

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
  // initialize recursion coefficients
  void InitializeNorms();

};

#endif
