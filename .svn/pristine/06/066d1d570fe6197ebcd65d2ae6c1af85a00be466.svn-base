////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of rational numbers using long long               //
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


#ifndef LONGRATIONAL_H
#define LONGRATIONAL_H


#include "config.h"

#include <iostream>
#include <fstream>

#ifdef __GMP__
#include <gmp.h>
#endif

using std::ostream;
using std::ofstream;
using std::ifstream;


class LongRational
{

 private:

#ifdef __GMP__
  // GMP rational number
  mpq_t Value;
#else
  // numerator
  LONGLONG Numerator;
  // denominator
  LONGLONG Denominator;

#endif

 public:

  // default constructor 
  //
  LongRational();

  // constructor from an integer
  //
  // x = value to assign to the rational coefficient
  LongRational(long x);  

  // constructor from a rational number
  //
  // x = numerator to assign to the rational coefficient
  // y = denominator to assign to the rational coefficient
  LongRational(long x, long y);  

  // copy constructor from a rational number
  //
  // x =  rational coefficient to copy
  LongRational(const LongRational& x);  

  // destructor
  //
  ~LongRational();

  // assignement
  //
  // rational = rational coefficient to assign
  // return value = reference on current rational coefficient
  LongRational& operator = (const LongRational& rational);

  // assignement from integer number
  //
  // x = interger to assign
  // return value = reference on current rational coefficient
  LongRational& operator = (long x);

  // assignement from a rational number encoded as a string
  //
  // x = string 
  // return value = reference on current rational coefficient
  LongRational& operator = (char* x);

  // set the coefficient to one
  //
  // return value = reference on current coefficient
  LongRational& SetToOne();

  // compute the opposit number
  //
  // x = first rational
  // return value = ooposit number
  friend LongRational operator - (const LongRational& x);

  // multiply the current rational by -1
  //
  // return value = reference on the current rational
  LongRational& Neg();

  // add two rational numbers
  //
  // x = first rational
  // y = second rational
  // return value = sum
  friend LongRational operator + (const LongRational& x, const LongRational& y);

  // add a rational number and an integer
  //
  // x = rational number
  // y = integer
  // return value = sum
  friend LongRational operator + (const LongRational& x, long y);

  // add a rational number and an integer
  //
  // y = rational number
  // x = integer
  // return value = sum
  friend LongRational operator + (long y, const LongRational& x);

  // substract an rational number to another
  //
  // x = first rational
  // y = second rational
  // return value = substraction
  friend LongRational operator - (const LongRational& x, const LongRational& y);

  // substract an integer to a rational number
  //
  // x = rational number
  // y = integer
  // return value = substraction
  friend LongRational operator - (const LongRational& x, long y);

  // substract a rational number to an integer
  //
  // y = integer
  // x = rational number
  // return value = substraction
  friend LongRational operator - (long y, const LongRational& x);

  // multiply two rational numbers
  //
  // x = first rational
  // y = second rational
  // return value = product
  friend LongRational operator * (const LongRational& x, const LongRational& y);

  // multiply a rational number and an integer
  //
  // x = rational number
  // y = integer
  // return value = product
  friend LongRational operator * (const LongRational& x, long y);

  // multiply a rational number and an integer
  //
  // y = rational number
  // x = integer
  // return value = product
  friend LongRational operator * (long y, const LongRational& x);

  // divide two rational numbers
  //
  // x = first rational
  // y = second rational
  // return value = division
  friend LongRational operator / (const LongRational& x, const LongRational& y);

  // divide a rational number by an integer
  //
  // x = rational number
  // y = integer
  // return value = division
  friend LongRational operator / (const LongRational& x, long y);

  // divide an integer by a rational number
  //
  // y = rational number
  // x = integer
  // return value = division
  friend LongRational operator / (long y, const LongRational& x);

  // add a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  LongRational& operator += (const  LongRational& x);

  // add an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  LongRational& operator += (long x);

  // substract a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  LongRational& operator -= (const  LongRational& x);

  // substract an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  LongRational& operator -= (long x);

  // multiply by an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  LongRational& operator *= (long x);

  // multiply by a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  LongRational& operator *= (const LongRational& x);

  // divide by an integer
  //
  // y = integer to use
  // return value = reference on current coefficient
  LongRational& operator /= (long y);

  // divide by a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  LongRational& operator /= (const LongRational& x);

  // multiply by x! factorial
  //
  // x = value of x (x!)
  // return value = reference on current coefficient
  LongRational& FactorialMultiply (long x);

  // divide by x! factorial
  //
  // x = value of x (x!)
  // return value = reference on current coefficient
  LongRational& FactorialDivide (long x);

  // multiply by a binomial factor
  //
  // m = major index
  // n = minor index
  // return value = reference on current coefficient  
  LongRational& BinomialMultiply (long m, long n);

  // divide by a binomial factor
  //
  // m = major index
  // n = minor index
  // return value = reference on current coefficient
  LongRational& BinomialDivide (long m, long n);

  // multiply the current rational by 2^x
  // 
  // x = 2 power exponent
  // return value = referencce on the current rational
  LongRational& Power2Multiply (long x);

  // divide the current rational by 2^x
  // 
  // x = 2 power exponent
  // return value = referencce on the current rational
  LongRational& Power2Divide (long x);

  // swap two rationals in an effecient way
  //
  // x = first rational
  // y = second rational
  friend void Swap(LongRational& x, LongRational& y);

  // test is two rational numbers are identical
  // 
  // x = first rational number
  // y = first rational number
  // return value = true if the two numbers are identical
  friend bool operator == (const LongRational& x, const LongRational& y);

  // test is two rational numbers are different
  // 
  // x = first rational number
  // y = first rational number
  // return value = true if the two numbers are different
  friend bool operator != (const LongRational& x, const LongRational& y);

  // test is a rational numbers is zero
  // 
  // x = rational number
  // return value = true if the number is zero
  bool IsZero ();


#ifdef __GMP__

  // test is a rational number is equal to an integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the two numbers are identical
  friend bool operator == (const LongRational& x, long y);

  // test is a rational number is equal to an integer number
  // 
  // y = integer number
  // x = rational number
  // return value = true if the two numbers are identical
  friend bool operator == ( long y, const LongRational& x);

  // test is a rational number is different from a given integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the two numbers are different
  friend bool operator != (const LongRational& x, long y);

  // test is a rational number is different from a given integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the two numbers are different
  friend bool operator != (long y, const LongRational& x);

#else

  // test is a rational number is equal to an integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the two numbers are identical
  friend bool operator == (const LongRational& x, LONGLONG y);

  // test is a rational number is equal to an integer number
  // 
  // y = integer number
  // x = rational number
  // return value = true if the two numbers are identical
  friend bool operator == (LONGLONG y, const LongRational& x);

  // test is a rational number is different from a given integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the two numbers are different
  friend bool operator != (const LongRational& x, LONGLONG y);

  // test is a rational number is different from a given integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the two numbers are different
  friend bool operator != (LONGLONG y, const LongRational& x);

#endif

  // return integer value associated to the coefficient numerator (0 if the coefficient can't be cast into an integer)
  //
  // return value = numerical coefficient
  long Num();

  // return integer value associated to the coefficient denominator (0 if the coefficient can't be cast into an integer)
  //
  // return value = numerical value associated to the coefficient denominator
  long Den();

  // return numerical value associated to the coefficient
  //
  // return value = numerical value associated to the coefficient
  double GetNumericalValue();

  // find the closest rational to a given double with maximum value for the denominator
  //
  // x = double to compare with
  // maximumDenominator = maximum value for the denominator
  // return value = difference between the double number in the closest rational
  double GetClosestRational(double x, long maximumDenominator);

  // return string associated to the coefficient
  //
  // division = character to use instead of '/'
  // return value = string associated to the coefficient
  char* GetString(char division = '/');

  // Output stream overload
  //
  // str = reference on output stream
  // x = rational to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, LongRational& x);

  // write a rational in binary mode
  //
  // file = reference on the output file stream
  // return value = reference on the output file stream
  ofstream& Write(ofstream& file);

  // read a rational in binary mode
  //
  // file = reference on the output file stream
  // return value = reference on the output file stream
  ifstream& Read(ifstream& file);

 private:

  // Simplify the current rational number
  //
  void Simplify();

#ifndef __GMP__

  // find greatest common divider
  //
  // m = first integer  
  // n = second integer (must be greater than m)
  // return value = GCD
  long FindGCD(LONGLONG m, LONGLONG n);

#endif

};

#ifdef __GMP__

// return integer value associated to the coefficient numerator (0 if the coefficient can't be cast into an integer)
//
// return value = numerical coefficient

inline long LongRational::Num()
{
  return 0l;
}

// return integer value associated to the coefficient denominator (0 if the coefficient can't be cast into an integer)
//
// return value = numerical value associated to the coefficient denominator

inline long LongRational::Den()
{
  return 0l;
}

// return numerical value associated to the coefficient
//
// return value = numerical value associated to the coefficient

inline double LongRational::GetNumericalValue()
{
  return mpq_get_d(this->Value);
}

// swap two rationals in an effecient way
//
// x = first rational
// y = second rational

inline void Swap(LongRational& x, LongRational& y)
{
  mpq_swap(x.Value, y.Value);
}

// Simplify the current rational number
//

inline void LongRational::Simplify()
{
  mpq_canonicalize(this->Value);
}

// test is two rational numbers are identical
// 
// x = first rational number
// y = first rational number
// return value = true if the two numbers are identical

inline bool operator == (const LongRational& x, const LongRational& y)
{
  return mpq_equal(x.Value, y .Value);
}

// test is a rational number is equal to an integer number
// 
// x = rational number
// y = integer number
// return value = true if the two numbers are identical

inline bool operator == (const LongRational& x, long y)
{
  return (mpq_cmp_si (x.Value, y, 1ul) == 0);
}

// test is a rational number is equal to an integer number
// 
// y = integer number
// x = rational number
// return value = true if the two numbers are identical

inline bool operator == ( long y, const LongRational& x)
{
  return (mpq_cmp_si (x.Value, y, 1ul) == 0);
}


// test is two rational numbers are different
// 
// x = first rational number
// y = first rational number
// return value = true if the two numbers are different

inline bool operator != (const LongRational& x, const LongRational& y)
{
  return !(mpq_equal(x.Value, y .Value));
}

// test is a rational number is different from a given integer number
// 
// x = rational number
// y = integer number
// return value = true if the two numbers are different

inline bool operator != (const LongRational& x, long y)
{
  return (mpq_cmp_si (x.Value, y, 1ul) != 0);
}

// test is a rational number is different from a given integer number
// 
// x = rational number
// y = integer number
// return value = true if the two numbers are different

inline bool operator != (long y, const LongRational& x)
{
  return (mpq_cmp_si (x.Value, y, 1ul) != 0);
}

// test is a rational numbers is zero
// 
// x = rational number
// return value = true if the number is zero

inline bool LongRational::IsZero ()
{ 
  return (mpq_sgn(this->Value) == 0);
}

// compute the opposit number
//
// x = first rational
// return value = ooposit number

inline LongRational operator - (const LongRational& x)
{
  LongRational Tmp(x);
  mpq_neg(Tmp.Value, Tmp.Value);
  return Tmp;
}

// multiply the current rational by -1
//
// return value = reference on the current rational

inline LongRational& LongRational::Neg()
{
  mpq_neg(this->Value, this->Value);
  return *this;
}

// multiply two rational numbers
//
// x = first rational
// y = second rational
// return value = product

inline LongRational operator * (const LongRational& x, const LongRational& y)
{
  LongRational Tmp;
  mpq_mul(Tmp.Value, x.Value, y.Value);
  return Tmp;
}

// add a rational
//
// x = rational to use
// return value = reference on current coefficient

inline LongRational& LongRational::operator += (const  LongRational& x)
{
  mpq_add(this->Value, this->Value, x.Value);  
  return *this;
}

// multiply by a rational
//
// x = rational to use
// return value = reference on current coefficient

inline LongRational& LongRational::operator *= (const LongRational& x)
{
  mpq_mul(this->Value, this->Value, x.Value);  
  return *this;
}

// divide by a rational
//
// x = rational to use
// return value = reference on current coefficient

inline LongRational& LongRational::operator /= (const LongRational& x)
{
  mpq_div(this->Value, this->Value, x.Value);  
  return *this;
}

#else

// return integer value associated to the coefficient numerator (0 if the coefficient can't be cast into an integer)
//
// return value = numerical coefficient

inline long LongRational::Num()
{
  return this->Numerator;
}

// return integer value associated to the coefficient denominator (0 if the coefficient can't be cast into an integer)
//
// return value = numerical value associated to the coefficient denominator

inline long LongRational::Den()
{
  return this->Denominator;
}

// return numerical value associated to the coefficient
//
// return value = numerical value associated to the coefficient

inline double LongRational::GetNumericalValue()
{
  return (((double) this->Numerator) / ((double) this->Denominator));
}

// Simplify the current rational number
//

inline void LongRational::Simplify()
{
  long Tmp = this->FindGCD(this->Numerator, this->Denominator);
  this->Numerator /= Tmp;
  this->Denominator /= Tmp;
}

// find greatest common divider
//
// m = first integer  
// n = second integer (must be greater than m)
// return value = GCD

inline long LongRational::FindGCD(LONGLONG m, LONGLONG n)
{
  LONGLONG Tmp;
  while (n != 0l)
    {
      Tmp = n;
      n = m % n;
      m = Tmp;
    }
  return m;
}

// swap two rationals in an effecient way
//
// x = first rational
// y = second rational

inline void Swap(LongRational& x, LongRational& y)
{
  LongRational Tmp = x;
  x = y;
  y = Tmp;
}

// test is two rational numbers are identical
// 
// x = first rational number
// y = first rational number
// return value = true if the two numbers are identical

inline bool operator == (const LongRational& x, const LongRational& y)
{
  if (((x.Numerator == y.Numerator) && (x.Denominator == y.Denominator)) || ((x.Numerator == -y.Numerator) && (x.Denominator == -y.Denominator)))
    return true;
  else
    return false;    
}


// test is two rational numbers are different
// 
// x = first rational number
// y = first rational number
// return value = true if the two numbers are different

inline bool operator != (const LongRational& x, const LongRational& y)
{
  if (((x.Numerator != y.Numerator) || (x.Denominator != y.Denominator)) && ((x.Numerator != -y.Numerator) || (x.Denominator != -y.Denominator)))
    return true;
  else
    return false;    
}

// test is a rational number is equal to an integer number
// 
// x = rational number
// y = integer number
// return value = true if the two numbers are identical

inline bool operator == (const LongRational& x, LONGLONG y)
{
  if (((x.Denominator == 1l) && (x.Numerator == y)) || ((x.Denominator == -1l) && (x.Numerator == -y)))
    return true;    
  else
    return false;
}

// test is a rational number is equal to an integer number
// 
// x = rational number
// y = integer number
// return value = true if the two numbers are identical

inline bool operator == (LONGLONG y, const LongRational& x)
{
  if (((x.Denominator == 1l) && (x.Numerator == y)) || ((x.Denominator == -1l) && (x.Numerator == -y)))
    return true;    
  else
    return false;
}

// test is a rational number is different from a given integer number
// 
// x = rational number
// y = integer number
// return value = true if the two numbers are different

inline bool operator != (const LongRational& x, LONGLONG y)
{
  if (((x.Denominator == 1l) && (x.Numerator == y)) || ((x.Denominator == -1l) && (x.Numerator == -y)))
    return false;
  else
    return true;    
}

// test is a rational number is different from a given integer number
// 
// x = rational number
// y = integer number
// return value = true if the two numbers are different

inline bool operator != (LONGLONG y, const LongRational& x)
{
  if (((x.Denominator == 1l) && (x.Numerator == y)) || ((x.Denominator == -1l) && (x.Numerator == -y)))
    return false;
  else
    return true;    
}

// test is a rational numbers is zero
// 
// x = rational number
// return value = true if the number is zero

inline bool LongRational::IsZero ()
{ 
  return (this->Numerator == ((LONGLONG) 0l));
}

// compute the opposit number
//
// x = first rational
// return value = ooposit number

inline LongRational operator - (const LongRational& x)
{
  LongRational Tmp(x);
  Tmp.Numerator *= -1l;
  return Tmp;
}

// multiply the current rational by -1
//
// return value = reference on the current rational

inline LongRational& LongRational::Neg()
{
  this->Numerator *= -1l;
  return *this;
}

#endif

#endif


