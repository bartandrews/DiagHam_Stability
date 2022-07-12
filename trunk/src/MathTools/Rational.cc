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


#include "config.h"
#include "MathTools/Rational.h"
#include "GeneralTools/Endian.h"

#include <iostream>
#include <cstring>
#include <cstdlib>


using std::cout;
using std::endl;


// default constructor 
//

Rational::Rational()
{
  this->Numerator = 0l;
  this->Denominator = 1l;
}

// constructor from an integer
//
// x = value to assign to the rational coefficient

Rational::Rational(long x)  
{
  this->Numerator = x;
  this->Denominator = 1l;
}

// constructor from a rational number
//
// x = numerator to assign to the rational coefficient
// y = denominator to assign to the rational coefficient

Rational::Rational(long x, long y)  
{
  this->Numerator = x;
  this->Denominator = y;
}

// destructor
//

Rational::~Rational()
{
}

// assignement
//
// rational = rational coefficient to assign
// return value = reference on current rational coefficient

Rational& Rational::operator = (const Rational& rational)
{
  this->Numerator = rational.Numerator;
  this->Denominator = rational.Denominator;
  return *this;
}

// assignement from integer number
//
// x = interger to assign
// return value = reference on current rational coefficient

Rational& Rational::operator = (long x)
{
  this->Numerator = x;
  this->Denominator = 1.0;
  return *this;
}

// assignement from a rational number encoded as a string
//
// x = string 
// return value = reference on current rational coefficient

Rational& Rational::operator = (char* x)
{
  char* TmpPos = strstr (x, "/");
  if (TmpPos == 0)
    {
      char* TmpError;
      this->Numerator =  strtol(x, &TmpError, 0);
      this->Denominator = 1l;
    }
  else
    {
      char* TmpError;
      (*TmpPos) = '\0';
      this->Numerator =  strtol(x, &TmpError, 0);
      this->Denominator = strtol(TmpPos + 1, &TmpError, 0);
    }
  this->Simplify();
  return *this;
}

// set the coefficient to one
//
// return value = reference on current coefficient

Rational& Rational::SetToOne()
{
  this->Numerator = 1l;
  this->Denominator = 1l;
  return *this;
}

// add two rational numbers
//
// x = first rational
// y = second rational
// return value = sum

Rational operator + (const Rational& x, const Rational& y)
{
  Rational Tmp;
  long Tmp2 = Tmp.FindGCD(x.Numerator, y.Denominator) ;
  long Tmp3 = Tmp.FindGCD(y.Numerator, x.Denominator) ;
  long Tmp4 = Tmp.FindGCD(Tmp2, Tmp3);
  Tmp.Numerator = ((x.Numerator / Tmp4) * (y.Denominator / Tmp4)) + ((y.Numerator / Tmp4) * (x.Denominator / Tmp4));  
  Tmp.Numerator *= Tmp4;
  Tmp.Numerator *= Tmp4;
  Tmp.Denominator = x.Denominator;
  Tmp.Simplify();  
  Tmp.Denominator *= y.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// add a rational number and an integer
//
// x = rational number
// y = integer
// return value = sum

Rational operator + (const Rational& x, long y)
{
  Rational Tmp;
  Tmp.Numerator = x.Numerator + (y * x.Denominator);
  Tmp.Denominator = x.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// add a rational number and an integer
//
// y = rational number
// x = integer
// return value = sum

Rational operator + (long y, const Rational& x)
{
  Rational Tmp;
  Tmp.Numerator = x.Numerator + (y * x.Denominator);
  Tmp.Denominator = x.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// substract an rational number to another
//
// x = first rational
// y = second rational
// return value = substraction

Rational operator - (const Rational& x, const Rational& y)
{
  Rational Tmp;
  Tmp.Numerator = (x.Numerator * y.Denominator) - (y.Numerator * x.Denominator);
  Tmp.Denominator = x.Denominator;
  Tmp.Simplify();
  Tmp.Denominator *= y.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// substract an integer to a rational number
//
// x = rational number
// y = integer
// return value = substraction

Rational operator - (const Rational& x, long y)
{
  Rational Tmp;
  Tmp.Numerator = x.Numerator - (y * x.Denominator);
  Tmp.Denominator = x.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// substract a rational number to an integer
//
// y = integer
// x = rational number
// return value = substraction

Rational operator - (long y, const Rational& x)
{
  Rational Tmp;
  Tmp.Numerator = (y * x.Denominator) - x.Numerator;
  Tmp.Denominator = x.Denominator;
  Tmp.Simplify();
  return Tmp;
}

// multiply two rational numbers
//
// x = first rational
// y = second rational
// return value = product

Rational operator * (const Rational& x, const Rational& y)
{
  Rational Tmp = x;
  long Tmp2 = y.Numerator;
  long Tmp3 = Tmp.FindGCD(Tmp2, Tmp.Denominator);
  long Tmp4 = y.Denominator;
  long Tmp5 = Tmp.FindGCD(Tmp4, Tmp.Numerator);
  Tmp.Denominator /= Tmp3;
  Tmp2 /= Tmp3;
  Tmp.Numerator /= Tmp5;
  Tmp4 /= Tmp5;
  Tmp.Numerator *= Tmp2;
  Tmp.Denominator *= Tmp4;
  return Tmp;
}

// multiply a rational number and an integer
//
// x = rational number
// y = integer
// return value = product

Rational operator * (const Rational& x, long y)
{
  Rational Tmp = x;
  long Tmp2 = Tmp.FindGCD(y, Tmp.Denominator);
  Tmp.Denominator /= Tmp2;
  y /= Tmp2;
  Tmp.Numerator *= y;
  return Tmp;
}

// multiply a rational number and an integer
//
// y = rational number
// x = integer
// return value = product

Rational operator * (long y, const Rational& x)
{
  Rational Tmp = x;
  long Tmp2 = Tmp.FindGCD(y, Tmp.Denominator);
  Tmp.Denominator /= Tmp2;
  y /= Tmp2;
  Tmp.Numerator *= y;
  return Tmp;
}

// divide two rational numbers
//
// x = first rational
// y = second rational
// return value = division

Rational operator / (const Rational& x, const Rational& y)
{
  Rational Tmp = x;
  long Tmp2 = y.Denominator;
  long Tmp3 = Tmp.FindGCD(Tmp2, Tmp.Denominator);
  long Tmp4 = y.Numerator;
  long Tmp5 = Tmp.FindGCD(Tmp4, Tmp.Numerator);
  Tmp.Denominator /= Tmp3;
  Tmp2 /= Tmp3;
  Tmp.Numerator /= Tmp5;
  Tmp4 /= Tmp5;
  Tmp.Numerator *= Tmp2;
  Tmp.Denominator *= Tmp4;
  return Tmp;
}

// divide a rational number ny an integer
//
// x = rational number
// y = integer
// return value = division

Rational operator / (const Rational& x, long y)
{
  Rational Tmp (x);
  long Tmp2 = Tmp.FindGCD(y, Tmp.Numerator);
  Tmp.Numerator /= Tmp2;
  y /= Tmp2;
  Tmp.Denominator *= y;
  return Tmp;
}

// divide an integer ny a rational number
//
// y = rational number
// x = integer
// return value = division

Rational operator / (long y, const Rational& x)
{
  Rational Tmp;
  Tmp.Denominator = x.Numerator;
  Tmp.Numerator = x.Denominator;
  long Tmp2 = Tmp.FindGCD(y, Tmp.Denominator);
  Tmp.Denominator /= Tmp2;
  y /= Tmp2;
  Tmp.Numerator *= y;
  return Tmp;
}

// add a rational
//
// x = rational to use
// return value = reference on current coefficient

Rational& Rational::operator += (const  Rational& x)
{
  this->Numerator *= x.Denominator;
  this->Numerator += x.Numerator * this->Denominator;
  this->Simplify() ;
  this->Denominator *= x.Denominator;
  this->Simplify() ;
  return *this;
}

// add an integer
//
// x = integer to use
// return value = reference on current coefficient

Rational& Rational::operator += (long x)
{
  this->Numerator += x * this->Denominator;
  this->Simplify();
  return *this;
}
 
// substract a rational
//
// x = rational to use
// return value = reference on current coefficient

Rational& Rational::operator -= (const  Rational& x)
{
  this->Numerator *= x.Denominator;
  this->Numerator -= x.Numerator * this->Denominator;
  this->Simplify();  
  this->Denominator *= x.Denominator;
  this->Simplify();  
  return *this;
}

// substract an integer
//
// x = integer to use
// return value = reference on current coefficient

Rational& Rational::operator -= (long x)
{
  this->Numerator -= x * this->Denominator;
  this->Simplify();
  return *this;
}

// multiply by an integer
//
// x = integer to use
// return value = reference on current coefficient

Rational& Rational::operator *= (long x)
{
  long Tmp = this->FindGCD(x, this->Denominator);
  x /= Tmp;
  this->Denominator /= Tmp;
  this->Numerator *= x;
  return *this;
}

// multiply by a rational
//
// x = rational to use
// return value = reference on current coefficient

Rational& Rational::operator *= (const Rational& x)
{
  long Tmp = this->FindGCD(this->Numerator, x.Denominator);
  long Tmp2 = this->FindGCD(x.Numerator, this->Denominator);
  this->Numerator /= Tmp;
  this->Denominator /= Tmp2;
  this->Numerator *=  x.Numerator / Tmp2;
  this->Denominator *= x.Denominator / Tmp;
  return *this;     
}

// divide by an integer
//
// y = integer to use
// return value = reference on current coefficient

Rational& Rational::operator /= (long y)
{
  long Tmp = this->FindGCD(y, this->Denominator);
  y /= Tmp;
  this->Denominator *= y;
  this->Numerator /= Tmp;
  return *this;
}

// divide by a rational
//
// x = rational to use
// return value = reference on current coefficient

Rational& Rational::operator /= (const Rational& x)
{
  long Tmp = this->FindGCD(this->Numerator, x.Numerator);
  long Tmp2 = this->FindGCD(x.Denominator, this->Denominator);
  this->Numerator /= Tmp;
  this->Denominator /= Tmp2;
  this->Numerator *=  x.Denominator / Tmp2;
  this->Denominator *= x.Numerator / Tmp;
  return *this;     
}

// Output stream overload
//
// str = reference on output stream
// x = rational to print
// return value = reference on output stream

ostream& operator << (ostream& str, Rational& x)
{
  if (x.Denominator > 0l)
    {
      if (x.Denominator == 1l)
	str << x.Numerator;
      else
	str << x.Numerator << "/" << x.Denominator;
    }
  else
    {
      if (x.Denominator == -1l)
	str << (-x.Numerator);
      else
	str << (-x.Numerator) << "/" << (-x.Denominator);
    }
  return str;
}


// write a rational in binary mode
//
// file = reference on the output file stream
// return value = reference on the output file stream

ofstream& Rational::Write(ofstream& file)
{
  WriteLittleEndian(file, this->Numerator);
  WriteLittleEndian(file, this->Denominator);
  return file;
}

// read a rational in binary mode
//
// file = reference on the output file stream
// return value = reference on the output file stream

ifstream& Rational::Read(ifstream& file)
{
  ReadLittleEndian(file, this->Numerator);
  ReadLittleEndian(file, this->Denominator);
  return file;
}
