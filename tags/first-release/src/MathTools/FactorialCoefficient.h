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

#include <iostream.h>


class FactorialCoefficient
{

  // current position in numerator
  int NumeratorPosition;
  // current position in denominator
  int DenominatorPosition;

  // array representing the numerator
  int* Numerator;
  // array representing the denominator
  int* Denominator;

 public:

  // default constructor 
  //
  FactorialCoefficient();

  // constructor from an integer
  //
  // x = value to assign to the factorial coefficient
  FactorialCoefficient(int x);  

  // constructor from a rational number
  //
  // x = numerator to assign to the factorial coefficient
  // y = denominator to assign to the factorial coefficient
  FactorialCoefficient(int x, int y);  

  // destructor
  //
  ~FactorialCoefficient();

  // set the coefficient to one
  //
  // return value = reference on current coefficient
  FactorialCoefficient& SetToOne();

  // multiply by an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  FactorialCoefficient& operator *= (int x);

  // divide by an integer
  //
  // y = integer to use
  // return value = reference on current coefficient
  FactorialCoefficient& operator /= (int y);

  // multiply the coefficient by a power of 2
  // 
  // power = power exponent (must be greater than 0)
  // return value = reference on current coefficient
  FactorialCoefficient& Power2Multiply (int power);

  // divide the coefficient by a power of 2
  // 
  // power = power exponent (must be greater than 0)
  // return value = reference on current coefficient
  FactorialCoefficient& Power2Divide (int power);

  // multiply the coefficient by the factorial of an integer
  // 
  // x = integer to use
  // return value = reference on current coefficient
  FactorialCoefficient& FactorialMultiply (int x);

  // multiply the coefficient by the partial factorial end! / (start - 1)!
  // 
  // start = first integer in the partial factorial product
  // end = last integer in the partial factorial product
  // return value = reference on current coefficient
  FactorialCoefficient& PartialFactorialMultiply (int start, int end);

  // multiply the coefficient by the factorial of an integer
  // 
  // x = integer to use
  // return value = reference on current coefficient
  FactorialCoefficient& FactorialDivide (int x);
  
  // divide the coefficient by the partial factorial end! / (start - 1)!
  // 
  // start = first integer in the partial factorial product
  // end = last integer in the partial factorial product
  // return value = reference on current coefficient
  FactorialCoefficient& PartialFactorialDivide (int start, int end);

  // return numerical value associated to the coefficient
  //
  // return value = numerical value associated to the coefficient
  double GetNumericalValue();

  // return integer value associated to the coefficient (0 if the coefficient is not an integer, or can't be cast into an integer)
  //
  // return value = numerical value associated to the coefficient
  int GetIntegerValue();

 private:

  // find greatest common divider (recursive part of the method)
  //
  // m = first integer  
  // n = second integer (must be greater than m)
  // return value = GCD
  int FindGCD(int m, int n);

  // find greatest common divider (recurisive part of the method)
  //
  // m = first integer  
  // n = second integer (must be greater than m)
  // return value = GCD
  int RecursiveFindGCD(int m, int n);

};

#endif


