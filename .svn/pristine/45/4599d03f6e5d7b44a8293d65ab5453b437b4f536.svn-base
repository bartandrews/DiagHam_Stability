////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of rational numbers                         //
//                                                                            //
//                        last modification : 12/11/2010                      //
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


#ifndef RATIONAL_H
#define RATIONAL_H


#include "config.h"

#include <iostream>
#include <fstream>


using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;


class Rational
{

 private:

  // numerator
  long Numerator;
  // denominator
  long Denominator;

 public:

  // default constructor 
  //
  Rational();

  // constructor from an integer
  //
  // x = value to assign to the rational coefficient
  Rational(long x);  

  // constructor from a rational number
  //
  // x = numerator to assign to the rational coefficient
  // y = denominator to assign to the rational coefficient
  Rational(long x, long y);  

  // destructor
  //
  ~Rational();

  // assignement
  //
  // rational = rational coefficient to assign
  // return value = reference on current rational coefficient
  Rational& operator = (const Rational& rational);

  // assignement from integer number
  //
  // x = interger to assign
  // return value = reference on current rational coefficient
  Rational& operator = (long x);

  // assignement from a rational number encoded as a string
  //
  // x = string 
  // return value = reference on current rational coefficient
  Rational& operator = (char* x);

  // set the coefficient to one
  //
  // return value = reference on current coefficient
  Rational& SetToOne();

  // compute the opposit number
  //
  // x = first rational
  // return value = ooposit number
  friend Rational operator - (const Rational& x);

  // add two rational numbers
  //
  // x = first rational
  // y = second rational
  // return value = sum
  friend Rational operator + (const Rational& x, const Rational& y);

  // add a rational number and an integer
  //
  // x = rational number
  // y = integer
  // return value = sum
  friend Rational operator + (const Rational& x, long y);

  // add a rational number and an integer
  //
  // y = rational number
  // x = integer
  // return value = sum
  friend Rational operator + (long y, const Rational& x);

  // substract an rational number to another
  //
  // x = first rational
  // y = second rational
  // return value = substraction
  friend Rational operator - (const Rational& x, const Rational& y);

  // substract an integer to a rational number
  //
  // x = rational number
  // y = integer
  // return value = substraction
  friend Rational operator - (const Rational& x, long y);

  // substract a rational number to an integer
  //
  // y = integer
  // x = rational number
  // return value = substraction
  friend Rational operator - (long y, const Rational& x);

  // multiply two rational numbers
  //
  // x = first rational
  // y = second rational
  // return value = product
  friend Rational operator * (const Rational& x, const Rational& y);

  // multiply a rational number and an integer
  //
  // x = rational number
  // y = integer
  // return value = product
  friend Rational operator * (const Rational& x, long y);

  // multiply a rational number and an integer
  //
  // y = rational number
  // x = integer
  // return value = product
  friend Rational operator * (long y, const Rational& x);

  // divide two rational numbers
  //
  // x = first rational
  // y = second rational
  // return value = division
  friend Rational operator / (const Rational& x, const Rational& y);

  // divide a rational number by an integer
  //
  // x = rational number
  // y = integer
  // return value = division
  friend Rational operator / (const Rational& x, long y);

  // divide an integer by a rational number
  //
  // y = rational number
  // x = integer
  // return value = division
  friend Rational operator / (long y, const Rational& x);

  // add a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  Rational& operator += (const  Rational& x);

  // add an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  Rational& operator += (long x);

  // substract a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  Rational& operator -= (const  Rational& x);

  // substract an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  Rational& operator -= (long x);

  // multiply by an integer
  //
  // x = integer to use
  // return value = reference on current coefficient
  Rational& operator *= (long x);

  // multiply by a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  Rational& operator *= (const Rational& x);

  // divide by an integer
  //
  // y = integer to use
  // return value = reference on current coefficient
  Rational& operator /= (long y);

  // divide by a rational
  //
  // x = rational to use
  // return value = reference on current coefficient
  Rational& operator /= (const Rational& x);

  // test is two rational numbers are identical
  // 
  // x = first rational number
  // y = first rational number
  // return value = true if the two numbers are identical
  friend bool operator == (const Rational& x, const Rational& y);

  // test is a rational number is equal to an integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the two numbers are identical
  friend bool operator == (const Rational& x, long y);

  // test is a rational number is equal to an integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the two numbers are identical
  friend bool operator == (long y, const Rational& x);

  // test is two rational numbers are different
  // 
  // x = first rational number
  // y = first rational number
  // return value = true if the two numbers are different
  friend bool operator != (const Rational& x, const Rational& y);

  // test is a rational number is different from a given integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the two numbers are different
  friend bool operator != (const Rational& x, long y);

  // test is a rational number is different from a given integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the two numbers are different
  friend bool operator != (long y, const Rational& x);

  // test is a rational number is less than integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the rational number is less than the integer number
  friend bool operator < (const Rational& x, long y);

  // test is a rational number is greater than integer number
  // 
  // x = rational number
  // y = integer number
  // return value = true if the rational number is greater than the integer number
  friend bool operator > (const Rational& x, long y);

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

  // Output stream overload
  //
  // str = reference on output stream
  // x = rational to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& str, Rational& x);

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

  // find greatest common divider
  //
  // m = first integer  
  // n = second integer (must be greater than m)
  // return value = GCD
  long FindGCD(long m, long n);


};

// return integer value associated to the coefficient numerator (0 if the coefficient can't be cast into an integer)
//
// return value = numerical coefficient

inline long Rational::Num()
{
  return this->Numerator;
}

// return integer value associated to the coefficient denominator (0 if the coefficient can't be cast into an integer)
//
// return value = numerical value associated to the coefficient denominator

inline long Rational::Den()
{
  return this->Denominator;
}

// return numerical value associated to the coefficient
//
// return value = numerical value associated to the coefficient

inline double Rational::GetNumericalValue()
{
  return (((double) this->Numerator) / ((double) this->Denominator));
}

// Simplify the current rational number
//

inline void Rational::Simplify()
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

inline long Rational::FindGCD(long m, long n)
{
  long Tmp;
  while (n != 0l)
    {
      Tmp = n;
      n = m % n;
      m = Tmp;
    }
  return m;
}

// test is two rational numbers are identical
// 
// x = first rational number
// y = first rational number
// return value = true if the two numbers are identical

inline bool operator == (const Rational& x, const Rational& y)
{
  if (((x.Numerator == y.Numerator) && (x.Denominator == y.Denominator)) || ((x.Numerator == -y.Numerator) && (x.Denominator == -y.Denominator)))
    return true;
  else
    return false;    
}

// test is a rational number is equal to an integer number
// 
// x = rational number
// y = integer number
// return value = true if the two numbers are identical

inline bool operator == (const Rational& x, long y)
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

inline bool operator == (long y, const Rational& x)
{
  if (((x.Denominator == 1l) && (x.Numerator == y)) || ((x.Denominator == -1l) && (x.Numerator == -y)))
    return true;    
  else
    return false;
}

// test is two rational numbers are different
// 
// x = first rational number
// y = first rational number
// return value = true if the two numbers are different

inline bool operator != (const Rational& x, const Rational& y)
{
  if (((x.Numerator != y.Numerator) || (x.Denominator != y.Denominator)) && ((x.Numerator != -y.Numerator) || (x.Denominator != -y.Denominator)))
    return true;
  else
    return false;    
}

// test is a rational number is different from a given integer number
// 
// x = rational number
// y = integer number
// return value = true if the two numbers are different

inline bool operator != (const Rational& x, long y)
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

inline bool operator != (long y, const Rational& x)
{
  if (((x.Denominator == 1l) && (x.Numerator == y)) || ((x.Denominator == -1l) && (x.Numerator == -y)))
    return false;
  else
    return true;    
}

// test is a rational number is less than integer number
// 
// x = rational number
// y = integer number
// return value = true if the rational number is less than the integer number

inline bool operator < (const Rational& x, long y)
{
  if (((x.Denominator > 0l) && (x.Numerator < (y * x.Denominator))) || ((x.Denominator < 0l) && (x.Numerator > (y * x.Denominator))))
    return true;
  else
    return false;    
}

// test is a rational number is greater than integer number
// 
// x = rational number
// y = integer number
// return value = true if the rational number is greater than the integer number

inline bool operator > (const Rational& x, long y)
{
  if (((x.Denominator > 0l) && (x.Numerator > (y * x.Denominator))) || ((x.Denominator < 0l) && (x.Numerator < (y * x.Denominator))))
    return true;
  else
    return false;    
}

// compute the opposit number
//
// x = first rational
// return value = ooposit number

inline Rational operator - (const Rational& x)
{
  Rational Tmp(x);
  Tmp.Numerator *= -1l;
  return Tmp;
}


#endif


