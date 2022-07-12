////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                           RayTrace version  0.10                           //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            Class for Complex                               //
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


#include "config.h"
#ifdef USE_OUTPUT
#include "Output/MathematicaOutput.h"
#endif

#include <math.h>
#include <iostream>


#ifndef COMPLEX_H
#define COMPLEX_H


using std::ostream;


class Complex;

//define i
Complex I();


class Complex
{

public:

  double Re, Im;

  //constructors
  Complex();
  Complex (double x1, double x2);   
  Complex (double x1);   
  Complex (const Complex& z);

  //destructor
  ~Complex() {};

  //assignement
  Complex& operator = (const Complex& z);
  Complex& operator = (double x);
  
// basic operations on complex numbers
// return real part of z
  friend double Real(const Complex& z);
// return imaginary part of z
  friend double Imag(const Complex& z);
  //return norm of z
  friend double Norm(const Complex& z);
  // return square of the z modulus 
  friend double SqrNorm(const Complex& z);
  // return argument of z 
  friend double Arg(const Complex& z);
  //return complex conjugate
  friend Complex Conj(const Complex& z);
  // return complex corresponding to the polar definition
  friend Complex Polar(double r,double theta);
  //return  the invert of a given complex number
  friend Complex Inv(const Complex& z);

// basic arithmetic operations
  Complex& operator += (const Complex& z);
  Complex& operator -= (const Complex& z);
  Complex& operator *= (const Complex& z);
  Complex& operator /= (const Complex& z);
  Complex& operator += (const double x);
  Complex& operator -= (const double x);
  Complex& operator *= (const double x);
  Complex& operator /= (const double x);
  friend Complex operator + (const Complex& z1,const Complex& z2);
  friend Complex operator - (const Complex& z1,const Complex& z2);
  friend Complex operator * (const Complex& z1,const Complex& z2);
  friend Complex operator / (const Complex& z1,const Complex& z2);
  friend Complex operator + (const Complex& z,double x);
  friend Complex operator - (const Complex& z,double x);
  friend Complex operator * (const Complex& z,double x);
  friend Complex operator / (const Complex& z,double x);
  friend Complex operator + (double x,const Complex& z);
  friend Complex operator - (double x,const Complex& z);
  friend Complex operator * (double x,const Complex& z);
  friend Complex operator / (double x,const Complex& z);
  friend Complex operator + (const Complex& z);
  friend Complex operator - (const Complex& z);

  // multiply the current complex number by the conjugate of another complex number 
  // 
  // z = reference on the complex number to use for the multiplication
  // return value = reference on the current complex number
  Complex& ConjugateProduct(const Complex& z);

  Complex& AddMultiply(const Complex& z, double x);

// logical operations
  friend bool operator == (const Complex& z1,const Complex& z2);
  friend bool operator != (const Complex& z1,const Complex& z2);

// basic mathematical functions
// exponential
  friend Complex exp (const Complex& z);
// logarithme
  friend Complex ln (const Complex& z);
// decimal logarithme
  friend Complex log (const Complex& z);
//hyperbolic cosine
  friend Complex ch (const Complex& z);
//hyperbolic sine
  friend Complex sh (const Complex& z);
//hyperbolic tangent
  friend Complex th (const Complex& z);
//hyperbolic cotangent
  friend Complex coth (const Complex& z);
//cosine
  friend Complex cos (const Complex& z);
//sine
  friend Complex sin (const Complex& z);
//tangent
  friend Complex tan (const Complex& z);
//cotangent
  friend Complex cotan (const Complex& z);
// square root (returning those corresponding to argument arg(z)/2)
  friend Complex sqrt (const Complex& z);
// arccosine
  friend Complex arccos (const Complex& z);
// arcsine
  friend Complex arcsin (const Complex& z);
// arctangent
  friend Complex arctan (const Complex& z);
//hyperbolic arccosine
  friend Complex argch (const Complex& z);
//hyperbolic arcsine
  friend Complex argsh (const Complex& z);
//hyperbolic arctangent
  friend Complex argth (const Complex& z);
// z power y with y real
  friend Complex pow (const Complex& z,const double y);
// z power y
  friend Complex pow (const Complex& z,const Complex& y);
// y power z with y real
  friend Complex pow (const double y,const Complex& z);

//define i
  friend Complex I();

  // Output Stream overload
  //
  // Str = output stream
  // z = complex value to print
  // return value = reference on current output stream 
  friend ostream& operator << (ostream& Str, const Complex& z);

#ifdef USE_OUTPUT
  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // z = complex value to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, const Complex& z);
#endif
};


//constructor

inline Complex::Complex()
{
  this->Re = 0.0;
  this->Im = 0.0;
}
  
inline Complex::Complex(double x1, double x2)
{
  this->Re = x1;
  this->Im = x2;
}
  
inline Complex::Complex(double x1)
{
  this->Re = x1;
  this->Im = 0.0;
}
  
inline Complex::Complex (const Complex& z)
{
  this->Re = z.Re;
  this->Im = z.Im;
}

//assignement

inline Complex& Complex::operator  = (const Complex& z)
{
  this->Re = z.Re;
  this->Im = z.Im;
  return *this;
}

inline Complex& Complex::operator  = (double x)
{
  this->Re = x;
  this->Im = 0.0;
  return *this;
}

// basic operations on complex numbers

// return real part of z

inline double Real(const Complex& z)
{
  return z.Re;
}

// return imaginary part of z

inline double Imag(const Complex& z)
{
  return z.Im;
}

//return norm of z

inline double Norm(const Complex& z)
{
  return sqrt(z.Re * z.Re + z.Im * z.Im);
}

// return square of the z modulus 
inline double SqrNorm(const Complex& z)
{
  return (z.Re * z.Re + z.Im * z.Im);
}

// return argument of z 

inline double Arg(const Complex& z)
{
  if (Norm(z) == 0.l)
    return 0.l;
  else
    return atan2(z.Im, z.Re);
}

//return complex conjugate

inline Complex Conj(const Complex& z)
{
  return Complex (z.Re, -z.Im);
}

// return complex corresponding to the polar definition

inline Complex Polar(double r,double theta)
{
  return Complex (r * cos(theta),r * sin(theta));
}

//return  the invert of a given complex number

inline Complex Inv(const Complex& z)
{
  double Tmp = 1.0 / SqrNorm(z);
  return Complex(z.Re * Tmp, -z.Im * Tmp);
}


// basic arithmetic operations

inline Complex& Complex::operator +=(const Complex& z) 
{
  this->Re += z.Re;
  this->Im += z.Im;
  return *this;
}

inline Complex& Complex::operator *=(const Complex& z) 
{
  double x;
  x = this->Re * z.Re - this->Im * z.Im;
  this->Im = this->Re * z.Im + this->Im * z.Re;
  this->Re = x;
  return *this;
}

inline Complex& Complex::operator -=(const Complex& z)
{
  this->Re -= z.Re;
  this->Im -= z.Im;
  return *this;
}

inline Complex& Complex::operator /=(const Complex& z)
{
  double n = 1.0 / (z.Re * z.Re + z.Im * z.Im);
  double x = (this->Re * z.Re + this->Im * z.Im) * n;
  this->Im = (this->Im * z.Re - this->Re * z.Im) * n;
  this->Re = x;
  return *this;
}

inline Complex& Complex::operator +=(const double x)
{
  this->Re += x;
  return *this;
}

inline Complex& Complex::operator *=(const double x)
{
  this->Re *= x;
  this->Im *= x;
  return *this;
}

inline Complex& Complex::operator -=(const double x)
{
  this->Re -= x;
  return *this;
}

inline Complex& Complex::operator /=(const double x)
{
  this->Re = this->Re / x;
  this->Im = this->Im / x;
  return *this;
}

inline Complex operator + (const Complex& z1, const Complex& z2)
{
  return Complex  (z1.Re + z2.Re, z1.Im + z2.Im);
}

inline Complex operator - (const Complex& z1, const Complex& z2)
{
  return Complex  (z1.Re - z2.Re, z1.Im - z2.Im);
}

inline Complex operator * (const Complex& z1, const Complex& z2)
{
  return Complex  (z1.Re * z2.Re - z1.Im * z2.Im, z1.Re * z2.Im + z2.Re * z1.Im);
}

inline Complex operator / (const Complex& z1, const Complex& z2)
{
  double n = z2.Re * z2.Re + z2.Im * z2.Im;
  return Complex  ((z1.Re * z2.Re + z1.Im * z2.Im) / n, (z2.Re * z1.Im - z1.Re * z2.Im) / n);
}

inline Complex operator + (const Complex& z, double x)
{
   return Complex  (z.Re + x, z.Im);
}

inline Complex operator + (double x, const Complex& z)
{
  return Complex  (z.Re + x, z.Im);
}

inline Complex operator - (const Complex& z, double x)
{
  return Complex  (z.Re - x, z.Im);
}

inline Complex operator - (double x, const Complex& z)
{
  return Complex  (x - z.Re, -z.Im);
}

inline Complex operator * (const Complex& z, double x)
{
  return Complex  (z.Re * x, z.Im * x);
}

inline Complex operator * (double x, const Complex& z)
{
  return Complex  (z.Re * x, z.Im * x);
}

inline Complex operator / (const Complex& z, double x)
{
  return Complex  (z.Re / x, z.Im / x);
}

inline Complex operator / (double x, const Complex& z) 
{
  double n = 1.0 / ((z.Re * z.Re) + (z.Im * z.Im));
  return Complex  ((x * z.Re) * n , (-x * z.Im) * n);
}

inline Complex operator + (const Complex& z)
{
  return Complex  (z.Re, z.Im);
}

inline Complex operator - (const Complex& z)
{
   return Complex  (-z.Re, -z.Im);
}

// logical operations

inline bool operator == (const Complex& z1,  const Complex& z2)
{
  return ((z1.Re == z2.Re) && (z1.Im == z2.Im));
}

inline bool operator != (const Complex& z1,  const Complex& z2)
{
  return ((z1.Re != z2.Re) || (z1.Im != z2.Im));
}

// multiply the current complex number by the conjugate of another complex number 
// 
// z = reference on the complex number to use for the multiplication
// return value = reference on the current complex number

inline Complex& Complex::ConjugateProduct(const Complex& z)
{
  double x;
  x = this->Re * z.Re + this->Im * z.Im;
  this->Im = this->Im * z.Re - this->Re * z.Im;
  this->Re = x;  
  return *this;
}

inline Complex& Complex::AddMultiply(const Complex& z, double x)
{
  this->Re+=z.Re*x;
  this->Im+=z.Im*x;
  return *this;
}
  

// basic mathematical functions

// exponential

inline Complex exp (const Complex& z)
{
  return Complex (exp(z.Re) * cos(z.Im), exp(z.Re) * sin(z.Im));
}

// logarithme

inline Complex ln (const Complex& z)
{
  return Complex (log(z.Re * z.Re + z.Im * z.Im) * 0.5l, atan2(z.Im, z.Re));
}

// decimal logarithme

inline Complex log (const Complex& z)
{
  return Complex (log(z.Re * z.Re + z.Im * z.Im)/(2.0l * log(10.0l)), atan2(z.Im, z.Re) / log(10.0l));
}

//hyperbolic cosine

inline Complex ch (const Complex& z)
{
  return ((exp (z) + exp(-z)) * 0.5l);
}

//hyperbolic sine

inline Complex sh (const Complex& z)
{
  return ((exp (z) - exp(-z)) * 0.5l);
}

//hyperbolic tangent

inline Complex th (const Complex& z)
{
  return ((exp (z) - exp(-z)) / (exp (z) + exp(-z)));
}

//hyperbolic cotangent

inline Complex coth (const Complex& z)
{
  return ((exp (z) + exp(-z)) / (exp (z) - exp(-z)));
}

//cosine

inline Complex cos (const Complex& z)
{
  return ch(I() * z);
}

//sine

inline Complex sin (const Complex& z)
{
  return - sh(I() * z) * I();
}

//tangent

inline Complex tan (const Complex& z)
{
  return (sin(z) / cos(z));
}

//cotangent

inline Complex cotan (const Complex& z)
{
  return (cos(z) / sin(z));
}

// square root (returning those corresponding to argument arg(z)/2)

inline Complex sqrt (const Complex& z)
{
  double mod = sqrt(Norm(z));
  double theta = Arg(z) * 0.5l;
  return Polar(mod, theta);
}

// arccosine

inline Complex arccos (const Complex& z)
{
  return - I() * ln(z + I() * sqrt(1.0l- z*z));
}

// arcsine

inline Complex arcsin (const Complex& z)
{
  return - I() * ln(I() * z + sqrt(1.0l - z * z));
}

//arctangent

inline Complex arctan (const Complex& z)
{
  return I() * ln((z + I()) / (I() - z)) * 0.5l;
}

//hyperbolic arccosine

inline Complex argch (const Complex& z)
{
  return ln(z + sqrt(z * z -  1.0l));
}

//hyperbolic arcsine

inline Complex argsh (const Complex& z)
{
  return ln(z + sqrt(1.0l + z * z));
}

//hyperbolic arctangent

inline Complex argth (const Complex& z)
{
  return ln((z + 1.0l) / (1.0l - z)) * 0.5l;
}

// z power y with y real

inline Complex pow (const Complex& z, const double y) 
{
  if ((y == 0.0) && (z == Complex (0.0, 0.0)))
    return Complex (1.0, 0.0);
  if (y == 0.0)
    return Complex (0.0, 0.0);
  return exp(y * ln(z));
}

// y power z with y real

inline Complex pow (const double y, const Complex& z)
{
  if ((y == 0.0) && ( z == Complex (0.0, 0.0)))
    return Complex (1.0, 0.0);
  if (Complex (0.0, 0.0) == z)
    return Complex (0.0, 0.0);
  return exp(z * log(y));
}

// z power y

inline Complex pow (const Complex& z, const Complex& y)
{
  if ((y == Complex (0.0, 0.0)) && (z == Complex (0.0, 0.0)))
    return Complex (1.0, 0.0);
  if (y == Complex (0.0, 0.0))
    return Complex (0.0, 0.0);
  return exp(y * ln(z));
}

//define i

inline Complex I()
{
  return Complex(0.0, 1.0);
}

// Output Stream overload
//
// Str = output stream
// z = complex value to print
// return value = reference on current output stream 

inline ostream& operator << (ostream& Str, const Complex& z)
{
  Str << "(" << z.Re << "," << z.Im << ")";
  return Str;
}

#ifdef  USE_OUTPUT
// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// z = complex value to print
// return value = reference on output stream

inline MathematicaOutput& operator << (MathematicaOutput& Str, const Complex& z)
{
  if (z.Im < 0.0)
    Str << z.Re << z.Im << "I";
  else
    Str << z.Re << "+" << z.Im << "I";
  return Str;
}
#endif

#endif

