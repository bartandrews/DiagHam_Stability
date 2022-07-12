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


#include "config.h"
#include "MathTools/FactorialCoefficient.h"


#define MAX_INTEGER_SIZE 1000


// default constructor 
//

FactorialCoefficient::FactorialCoefficient()
{
  this->NumeratorPosition = 0;
  this->DenominatorPosition = 0;
  this->Numerator = new int [MAX_INTEGER_SIZE];
  this->Denominator = new int [MAX_INTEGER_SIZE];
  this->Numerator[0] = 1;
  this->Denominator[0] = 1;
}

// constructor from an integer
//
// x = value to assign to the factorial coefficient

FactorialCoefficient::FactorialCoefficient(int x)  
{
  this->NumeratorPosition = 0;
  this->DenominatorPosition = 0;
  this->Numerator = new int [MAX_INTEGER_SIZE];
  this->Denominator = new int [MAX_INTEGER_SIZE];
  this->Numerator[0] = x;
  this->Denominator[0] = 1;
}

// constructor from a rational number
//
// x = numerator to assign to the factorial coefficient
// y = denominator to assign to the factorial coefficient

FactorialCoefficient::FactorialCoefficient(int x, int y)  
{
  this->NumeratorPosition = 0;
  this->DenominatorPosition = 0;
  this->Numerator = new int [MAX_INTEGER_SIZE];
  this->Denominator = new int [MAX_INTEGER_SIZE];
  this->Numerator[0] = x;
  this->Denominator[0] = y;
}

// destructor
//

FactorialCoefficient::~FactorialCoefficient()
{
  delete[] this->Numerator;
  delete[] this->Denominator;
}

// set the coefficient to one
//
// return value = reference on current coefficient

FactorialCoefficient& FactorialCoefficient::SetToOne()
{
  this->NumeratorPosition = 0;
  this->DenominatorPosition = 0;
  this->Numerator[0] = 1;
  this->Denominator[0] = 1;
  return *this;
}

// multiply by an integer
//
// x = integer to use
// return value = reference on current coefficient

FactorialCoefficient& FactorialCoefficient::operator *= (int x)
{
  int GCD;
  int i;
  for (i = 0; i <= this->DenominatorPosition; ++i)
    {
      GCD = this->FindGCD(x, this->Denominator[i]);
      x /= GCD;
      this->Denominator[i] /= GCD;
    }
  i = 0;
  while (i <= this->NumeratorPosition)
    {
      if ((x * this->Numerator[i]) < this->Numerator[i])
	{
	  this->Numerator[i] *= x;
	  i = this->NumeratorPosition + 2;
	}
      ++i;
    }
  if (i == (this->NumeratorPosition + 1))
    {
      this->Numerator[i] = x;
      ++this->NumeratorPosition;
    }
  return *this;
}

// divide by an integer
//
// y = integer to use
// return value = reference on current coefficient

FactorialCoefficient& FactorialCoefficient::operator /= (int y)
{
  int GCD;
  int i;
  for (i = 0; i <= this->NumeratorPosition; ++i)
    {
      GCD = this->FindGCD(y, this->Numerator[i]);
      y /= GCD;
      this->Numerator[i] /= GCD;
    }
  i = 0;
  while (i <= this->DenominatorPosition)
    {
      if ((y * this->Denominator[i]) < this->Denominator[i])
	{
	  this->Denominator[i] *= y;
	  i = this->DenominatorPosition + 2;
	}
      ++i;
    }
  if (i == (this->DenominatorPosition + 1))
    {
      this->Denominator[i] = y;
      ++this->DenominatorPosition;
    }
  return *this;
}

// multiply the coefficient by a power of 2
// 
// power = power exponent (must be greater than 0)
// return value = reference on current coefficient

FactorialCoefficient& FactorialCoefficient::Power2Multiply (int power)
{
  int i;
  for (i = 0; (i <= this->DenominatorPosition) && (power > 0); ++i)
    {
      while (!(this->Denominator[i] & 0x1) && (power > 0))
	{
	  this->Denominator[i] >>= 1;
	  --power;
	}
    }
  for (i = 0; (i <= this->NumeratorPosition) && (power > 0); ++i)
    {
      while (((this->Numerator[i] << 1) > 0) && (power > 0))
	{
	  this->Numerator[i] <<= 1;
	  --power;
	}
    }
  while (power > 0)
    {
      i = power % 31;
      power -= i;
      ++this->NumeratorPosition;
      this->Numerator[this->NumeratorPosition] = 1 << i;
    }
  return *this;
}

// divide the coefficient by a power of 2
// 
// power = power exponent (must be greater than 0)
// return value = reference on current coefficient

FactorialCoefficient& FactorialCoefficient::Power2Divide (int power)
{
  int i;
  for (i = 0; (i <= this->NumeratorPosition) && (power > 0); ++i)
    {
      while (!(this->Numerator[i] & 0x1) && (power > 0))
	{
	  this->Numerator[i] >>= 1;
	  --power;
	}
    }
  for (i = 0; (i <= this->DenominatorPosition) && (power > 0); ++i)
    {
      while (((this->Denominator[i] << 1) > 0) && (power > 0))
	{
	  this->Denominator[i] <<= 1;
	  --power;
	}
    }
  while (power > 0)
    {
      i = power % 31;
      power -= i;
      ++this->DenominatorPosition;
      this->Denominator[this->DenominatorPosition] = 1 << i;
    }
  return *this;
}

// multiply the coefficient by the factorial of an integer
// 
// x = integer to use
// return value = reference on current coefficient

FactorialCoefficient& FactorialCoefficient::FactorialMultiply (int x)
{
  if (x <= 1)
    return *this;
  for (int i = 2; i <= x; ++i)
    {
      (*this) *= i;
    }
  return *this;
}

// multiply the coefficient by the partial factorial end! / (start - 1)!
// 
// start = first integer in the partial factorial product
// end = last integer in the partial factorial product
// return value = reference on current coefficient

FactorialCoefficient& FactorialCoefficient::PartialFactorialMultiply (int start, int end)
{
  if (end <= 1)
    return *this;
  if (start <= 1)
    start = 2;
  for (int i = start; i <= end; ++i)
    {
      (*this) *= i;
    }
  return *this;
}

// multiply the coefficient by the factorial of an integer
// 
// x = integer to use
// return value = reference on current coefficient

FactorialCoefficient& FactorialCoefficient::FactorialDivide (int x)
{
  if (x <= 1)
    return *this;
  for (int i = 2; i <= x; ++i)
    {
      (*this) /= i;
    }
  return *this;
}

// divide the coefficient by the partial factorial end! / (start - 1)!
// 
// start = first integer in the partial factorial product
// end = last integer in the partial factorial product
// return value = reference on current coefficient

FactorialCoefficient& FactorialCoefficient::PartialFactorialDivide (int start, int end)
{
  if (end <= 1)
    return *this;
  if (start <= 1)
    start = 2;
  for (int i = start; i <= end; ++i)
    {
      (*this) /= i;
    }
  return *this;
}

// return numerical value associated to the coefficient
//
// return value = numerical value associated to the coefficient

double FactorialCoefficient::GetNumericalValue()
{
  double x = 1.0;
  double y = 1.0;
  for (int i = 0; i <= this->NumeratorPosition; ++i)
    x *= (double) this->Numerator[i];
  for (int i = 0; i <= this->DenominatorPosition; ++i)
    y *= (double) this->Denominator[i];
  return (x / y);
}

// return integer value associated to the coefficient (0 if the coefficient is not an integer, or can't be cast into an integer)
//
// return value = numerical value associated to the coefficient

int  FactorialCoefficient::GetIntegerValue()
{
  for (int i = 0; i <= this->DenominatorPosition; ++i)
    if (this->Denominator[i] != 1)
      return 0;
  int x = 1;
  for (int i = 0; i <= this->NumeratorPosition; ++i)
    if ((x * this->Numerator[i]) >= x)
      x *= this->Numerator[i];
    else
      return 0;
  return x;
}

// find greatest common divider (recursive part of the method)
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD
 
int FactorialCoefficient::FindGCD(int m, int n)
{
  if (m < n)
    return this->RecursiveFindGCD (m, n);
  else
    return this->RecursiveFindGCD (n, m);
  return n;
}

// find greatest common divider (recurisive part of the method)
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD
 
int FactorialCoefficient::RecursiveFindGCD(int m, int n)
{
  if (m == 0)
    return n;
  else
    return this->FindGCD ((n % m), m);
}

