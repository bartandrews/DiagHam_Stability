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


#ifndef JACOBITHETAFUNCTION_H
#define JACOBITHETAFUNCTION_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"
#include <iostream>


using std::ostream;


class JacobiThetaFunction
{

 protected:
  
  double ParameterA;
  
  double ParameterB;

  // modulus of the function
  Complex Tau;

  // overall phase from shifting B into standard interval [0,1]
  Complex ShiftPhase;

  // Offset of the sum over n 
  int SumOffset;

  
 public:

  // default constructor
  // creates \theta_1, with tau=I
  //
  JacobiThetaFunction();

  // constructor for a general theta function \theta[^a_b](z|tau)
  //
  // param_a = parameter a
  // param_b = parameter b
  // tau = modulus
  // precision = required precision
  JacobiThetaFunction(double a, double b, Complex tau=Complex(0.0,1.0));

  // constructor for one of the four Jacobi theta functions \theta_i(z|tau)
  //
  // type = index i of theta-function
  // tau = modulus
  // precision = required precision
  JacobiThetaFunction(int i, Complex tau=Complex(0.0,1.0));

  // copy constructor (without duplicating datas)
  //
  // coefficients = reference on Clebsch Gordan coefficients to copy
  JacobiThetaFunction (const JacobiThetaFunction& theta);

  // destructor
  //
  ~JacobiThetaFunction ();

  
  // get the value of the function for a given coordinate z
  //
  // z = complex coordinate
  // return = function value at z

  Complex GetValue(const Complex &z);

  // get the value of the function for a given coordinate z
  // values = complex vector where to store results
  // manyZ = complex vector of coordinates
  // return = function value at z
  void GetManyValues(ComplexVector &values, ComplexVector &manyZ);

  // pretty-print a function value
  // str = stream to print to
  // z = point where to evaluate
  ostream& PrintValue(ostream &str, const Complex &z);
  
  
};

#endif
