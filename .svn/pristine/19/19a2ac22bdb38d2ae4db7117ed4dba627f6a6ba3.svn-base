////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of factorial coefficient                      //
//                                                                            //
//                        last modification : 04/06/2002                      //
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


#ifndef FACTORIALCOEFFICIENT_H
#define FACTORIALCOEFFICIENT_H


#include "config.h"

#include <iostream>


class FactorialCoefficient
{

  // current position in numerator
  int NumeratorPosition;
  // current position in denominator
  int DenominatorPosition;

  // array representing the numerator
  long* Numerator;
  // array representing the denominator
  long* Denominator;

 public:

  // default constructor 
  //
  FactorialCoefficient();

  // constructor from an integer
  //
  // x = value to assign to the factorial coefficient
  FactorialCoefficient(long x);  

  // constructor from a rational number
  //
  // x = numerator to assign to the factorial coefficient
  // y = denominator to assign to the factorial coefficient
  FactorialCoefficient(long x, long y);  

  // destructor
  //
  ~FactorialCoefficient();

  // assignement
  //
  // factorial = factorial coefficient to assign
  // return value = reference on current factorial coefficient
  FactorialCoefficient& operator = (FactorialCoefficient& factorial);

  // set the coefficient to one
  //
  // return value = reference on current coefficient
  FactorialCoefficient& SetToOne();

  // multiply by an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  FactorialCoefficient& operator *= (long x);

  // divide by an integer
  //
  // y = integer to use
  // return value = reference on current coefficient
  FactorialCoefficient& operator /= (long y);

  // multiply the coefficient by a power of 2
  // 
  // power = power exponent (must be greater than 0)
  // return value = reference on current coefficient
  FactorialCoefficient& Power2Multiply (long power);

  // multiply the coefficient by a power of N
  // 
  // n = value whose  has to be risen to the power
  // power = power exponent (must be greater than 0)
  // return value = reference on current coefficient
  FactorialCoefficient& PowerNMultiply (long n, long power);

  // divide the coefficient by a power of 2
  // 
  // power = power exponent (must be greater than 0)
  // return value = reference on current coefficient
  FactorialCoefficient& Power2Divide (long power);

  // divide the coefficient by a power of N
  // 
  // n = value whose power has to be risen to the power
  // power = power exponent (must be greater than 0)
  // return value = reference on current coefficient
  FactorialCoefficient& PowerNDivide (long n, long power);

  // multiply the coefficient by the factorial of an integer
  // 
  // x = integer to use
  // return value = reference on current coefficient
  FactorialCoefficient& FactorialMultiply (long x);

  // multiply the coefficient by the partial factorial end! / (start - 1)!
  // 
  // start = first integer in the partial factorial product
  // end = last integer in the partial factorial product
  // return value = reference on current coefficient
  FactorialCoefficient& PartialFactorialMultiply (long start, long end);

  // multiply the coefficient by the factorial of an integer
  // 
  // x = integer to use
  // return value = reference on current coefficient
  FactorialCoefficient& FactorialDivide (long x);
  
  // divide the coefficient by the partial factorial end! / (start - 1)!
  // 
  // start = first integer in the partial factorial product
  // end = last integer in the partial factorial product
  // return value = reference on current coefficient
  FactorialCoefficient& PartialFactorialDivide (long start, long end);

  // return integer value associated to the coefficient numerator (0 if the coefficient can't be cast into an integer)
  //
  // return value = numerical value associated to the coefficient numerator  
  long GetIntegerNumeratorValue();

  // return integer value associated to the coefficient denominator (0 if the coefficient can't be cast into an integer)
  //
  // return value = numerical value associated to the coefficient denominator
  long GetIntegerDenominatorValue();

  // return numerical value associated to the coefficient
  //
  // return value = numerical value associated to the coefficient
  double GetNumericalValue();

  // return integer value associated to the coefficient (0 if the coefficient is not an integer, or can't be cast into an integer)
  //
  // return value = numerical value associated to the coefficient
  long GetIntegerValue();

 private:

  // find greatest common divider (recursive part of the method)
  //
  // m = first integer  
  // n = second integer (must be greater than m)
  // return value = GCD
  long FindGCD(long m, long n);

  // find greatest common divider (recurisive part of the method)
  //
  // m = first integer  
  // n = second integer (must be greater than m)
  // return value = GCD
  long RecursiveFindGCD(long m, long n);

};

#endif


