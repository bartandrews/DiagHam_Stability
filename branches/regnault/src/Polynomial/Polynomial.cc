////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                           RayTrace version  0.10                           //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   Class of polynomial with real coefficients               //
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


#include "Polynomial/Polynomial.h"
#include <math.h>


#define EPSILON 1E-6


using std::cout;
using std::endl;


// return absolute value of a given number
//
// x = number to proceed
// return value = absolute value of x
double Abs (double x);

// test if an integer is odd
//
// x = integer to test
// return value = true if x is odd
double TestOdd (int  x);

// constructors

Polynomial::Polynomial (int Deg, double* Coef)
{
  this->Degree = Deg;
  this->Coefficient = new double [this->Degree + 1];
  this->RootFlag = false;
  this->NoRootFlag = false;
  for (int i = 0; i <= this->Degree; i++)
    this->Coefficient[i] = Coef[i];
}

Polynomial::Polynomial (int Deg, double* Coef, bool flag)
{
  this->Degree = Deg;
  this->NoRootFlag = false;
  if (flag = true)
    this->Coefficient = Coef;
  else
    {
      this->Coefficient = new double [this->Degree + 1];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] = Coef[i];
    }
  this->RootFlag = false;
}

Polynomial::Polynomial (const Polynomial& P)
{
  this->Degree = P.Degree;
  this->Coefficient = new double [this->Degree + 1];
  this->RootFlag = P.RootFlag;
  for (int i = 0; i <= this->Degree; i++)
    this->Coefficient[i] = P.Coefficient[i];
  this->NoRootFlag = P.NoRootFlag;
  if (P.RootFlag == true)
    {
      this->NbrRoot = P.NbrRoot;
      this->Root = new double [this->NbrRoot];
      for (int i = 0; i < this->NbrRoot; i++)
	this->Root[i] = P.Root[i];
    }
}

// destructor
//

Polynomial::~Polynomial()
{
  if ((this->RootFlag == true) && (this->NbrRoot != 0))
    delete this->Root;
  delete[] this->Coefficient;
}

// assignement

Polynomial& Polynomial::operator = (const Polynomial& P)
{
  delete this->Coefficient;
  this->Degree = P.Degree;
  this->RootFlag = P.RootFlag;
  this->Coefficient = new double [this->Degree + 1];
  if ((this->RootFlag == true) && (this->NbrRoot != 0))
    delete this->Root;
  this->RootFlag = P.RootFlag;
  for (int i = 0; i <= this->Degree; i++)
    this->Coefficient[i] = P.Coefficient[i];
  if (P.RootFlag == true)
    {
      this->NbrRoot = P.NbrRoot;
      this->Root = new double [this->NbrRoot];
      for (int i = 0; i < this->NbrRoot; i++)
	this->Root[i] = P.Root[i];
    }
  return *this;
}

// Return polynomial value at a given point
//
// x = point where to evaluate polynomial
// return value = polynomial value at x

double Polynomial::PolynomialEvaluate (double x)
{
  double Res = this->Coefficient[this->Degree];
  for (int i = this->Degree - 1 ; i >= 0; i--)
    {
      Res *= x;
      Res += this->Coefficient[i];
    }
  return Res;
}
  
// Return polynomial value at a given point
//
// x = point where to evaluate polynomial
// return value = polynomial value at x

Complex Polynomial::PolynomialEvaluate (Complex x)
{
  Complex Res = this->Coefficient[this->Degree];
  for (int i = this->Degree - 1 ; i >= 0; i--)
    {
      Res *= x;
      Res += this->Coefficient[i];
    }
  return Res;
}
  
// Evaluate polynomial derivative 
//
// x = position where to evaluate polynomial derivative
// return value = polynomial derivative at x

double Polynomial::DerivativeEvaluate (double x)
{
  double Res = Degree * this->Coefficient[this->Degree];
  for (int i = this->Degree - 1 ; i > 0; i--)
    {
      Res *= x;
      Res += i * this->Coefficient[i];
    }
  return Res;
}
  
// Evaluate polynomial derivative 
//
// x = position where to evaluate polynomial derivative
// return value = polynomial derivative at x

Complex Polynomial::DerivativeEvaluate (Complex x)
{
  Complex Res = Degree * this->Coefficient[this->Degree];
  for (int i = this->Degree - 1 ; i > 0; i--)
    {
      Res *= x;
      Res += i * this->Coefficient[i];
    }
  return Res;
}
  
// Evaluate polynomial n-th derivative at a given value
//
// x = position where to evaluate polynomial n-th derivative
// n = derivative order
// return value = polynomial n-th derivative at x

double Polynomial::DerivativeEvaluate (double x, int n)
{
  if (n > this->Degree)
    return 0.0;
  double tmp = (double) 1.0;
  for (int j = this->Degree - n + 1; j <= this->Degree; j++)
    tmp *= (double) j;
  double Res = tmp * this->Coefficient[this->Degree];
  for (int i = this->Degree - 1; i >= n; i--)
    {
      tmp = (double) 1.0;
      for (int j = i - n + 1; j <= i; j++)
	tmp *= (double) j;      
      Res *= x;
      Res += tmp * this->Coefficient[i];
    }
  return Res;
}
  
// Evaluate polynomial n-th derivative at a given value
//
// x = position where to evaluate polynomial n-th derivative
// n = derivative order
// return value = polynomial n-th derivative at x

Complex Polynomial::DerivativeEvaluate (Complex x, int n)
{
  if (n > this->Degree)
    return Complex();
  double tmp = (double) 1.0;
  for (int j = this->Degree - n + 1; j <= this->Degree; j++)
    tmp *= (double) j;
  Complex Res (tmp * this->Coefficient[this->Degree]);
  for (int i = this->Degree - 1; i >= n; i--)
    {
      tmp = (double) 1.0;
      for (int j = i - n + 1; j <= i; j++)
	tmp *= (double) j;      
      Res *= x;
      Res += tmp * this->Coefficient[i];
    }
  return Res;
}
  
// Return Derivative of the polynomial 
//
// return value = polynomial derivative

Polynomial Polynomial::DerivatePolynomial ()
{
  double *Coef = new double [this->Degree];
  for (int i = 1; i <= this->Degree; i++)
    Coef[i-1] = (this->Coefficient[i] * i);
  Polynomial P = Polynomial(this->Degree - 1, Coef);
  delete Coef;
  return P;
}

// Evaluate polynomial n-th derivative
//
// n = derivative order
// return value = polynomial n-th derivative

Polynomial Polynomial::DerivatePolynomial (int n)
{
  if (n > this->Degree)
    {
      double* Coef = new double [1];
      Coef[0] = 0.0;
      return Polynomial(0, Coef, true);
    }
  int tmpDegree = this->Degree - n;
  double* Coef = new double [tmpDegree + 1];  
  double tmp = (double) 1.0;
  for (int j = this->Degree - n + 1; j <= this->Degree; j++)
    tmp *= (double) j;
  Coef[tmpDegree] = tmp * this->Coefficient[this->Degree];
  for (int i = this->Degree - 1; i >= n; i--)
    {
      tmp = (double) 1.0;
      for (int j = i - n + 1; j <= i; j++)
	tmp *= (double) j;      
      Coef[i - n] = tmp * this->Coefficient[i];
    }
  return  Polynomial(tmpDegree, Coef, true);
}

// arithmetic operators

Polynomial operator - (const Polynomial& P)
{
  Polynomial Pr(P);
  for (int i = 0; i <= Pr.Degree; i++)
    Pr.Coefficient[i] = -Pr.Coefficient[i];  
  return Pr;  
}

Polynomial operator + (const Polynomial& P1, const Polynomial& P2)
{
  if (P1.Degree > P2.Degree)    
    {
      double* Coef = new double [P1.Degree+1];
      for (int i = P1.Degree; i > P2.Degree; i--)
	Coef[i] = P1.Coefficient[i];
      for (int i = P2.Degree; i >= 0; i--)
	Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];	  
      return Polynomial(P1.Degree, Coef, true);
    }
  else
    if (P1.Degree < P2.Degree)    
      {
	double* Coef = new double [P2.Degree+1];
	for (int i = P2.Degree; i > P1.Degree; i--)
	  Coef[i] = P2.Coefficient[i];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];	  
	return Polynomial(P2.Degree, Coef, true);

      }
    else
      {
	double* Coef = new double [P1.Degree+1];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];
	int i = P1.Degree;
	while ((Coef[i] == 0) && (i!=0)) i--;
	if (i == P1.Degree)
	  return Polynomial(P2.Degree, Coef, true);
	else 
	  {
	    double* Coef2 = new double [i+1];
	    for (int j = 0; j <= i; j++)
	      Coef2[j] = Coef[j];
	    delete Coef;
	    return Polynomial(i, Coef2, true);
	  }
      }
}

Polynomial operator - (const Polynomial& P1, const Polynomial& P2)
{
  if (P1.Degree > P2.Degree)    
    {
      double* Coef = new double [P1.Degree+1];
      for (int i = P1.Degree; i > P2.Degree; i--)
	Coef[i] = P1.Coefficient[i];
      for (int i = P2.Degree; i >= 0; i--)
	Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];	  
      return Polynomial(P1.Degree, Coef, true);
    }
  else
    if (P1.Degree < P2.Degree)    
      {
	double* Coef = new double [P2.Degree+1];
	for (int i = P2.Degree; i > P1.Degree; i--)
	  Coef[i] = - P2.Coefficient[i];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];	  
	return Polynomial(P2.Degree, Coef, true);

      }
    else
      {
	double* Coef = new double [P1.Degree+1];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];
	int i = P1.Degree;
	while ((Coef[i] == 0) && (i!=0)) i--;
	if (i == P1.Degree)
	  return Polynomial(P2.Degree, Coef, true);
	else 
	  {
	    double* Coef2 = new double [i+1];
	    for (int j = 0; j <= i; j++)
	      Coef2[j] = Coef[j];
	    delete Coef;
	    return Polynomial(i, Coef2, true);
	  }
      }
}

Polynomial operator * (const Polynomial& P1, const Polynomial& P2)
{
  int Deg = P1.Degree+P2.Degree;
  double* Coef = new double [Deg + 1];
  for (int i = P1.Degree; i >= 0; i--)
    Coef[i+P2.Degree] = P1.Coefficient[i] * P2.Coefficient[P2.Degree];
  if (P2.Degree == 0)
    return Polynomial(Deg, Coef, true);
  for (int i = P2.Degree-1; i >= 0; i--)
    Coef[i] = 0;
  for (int i = P2.Degree - 1; i >= 0; i--)
    for (int j = P1.Degree; j >= 0; j--)
      Coef[i+j] += P1.Coefficient[j] * P2.Coefficient[i];
  return Polynomial(Deg, Coef, true);
}

Polynomial operator * (const Polynomial& P, const double& d)
{
  Polynomial Pr(P);
  for (int i = 0; i <= Pr.Degree; i++)
    Pr.Coefficient[i] *= d;
  return Pr;
}

Polynomial operator * (const double& d, const Polynomial& P)
{
  Polynomial Pr = P;
  for (int i = 0; i <= Pr.Degree; i++)
    Pr.Coefficient[i] *= d;
  return Pr;
}

Polynomial operator / (const Polynomial& P, const double& d)
{
  Polynomial Pr = P;
  for (int i = 0; i <= Pr.Degree; i++)
    Pr.Coefficient[i] /= d;
  return Pr;
}

Polynomial operator / (const Polynomial& P1, const Polynomial& P2)
{
  int Deg = P1.Degree - P2.Degree;
  double* Coef1 = new double [P1.Degree];
  for (int i = P1.Degree-1; i >= 0; i--)
    Coef1[i] = P1.Coefficient[i];
  double* Coef2 = new double [Deg+1];
  double c = 1 / P2.Coefficient[P2.Degree];
  Coef2[Deg] = P1.Coefficient[P1.Degree] * c;  
  for (int i = Deg; i > 0; i--)
    {
      for (int j = P2.Degree - 1; j >= 0; j--)
	Coef1[i + j] -= Coef2[i] * P2.Coefficient[j];
      Coef2[i-1] = c * Coef1[P2.Degree + i - 1]; 
    }
  delete[] Coef1;
  return Polynomial(Deg, Coef2, true);
}

Polynomial operator % (const Polynomial& P1, const Polynomial& P2)
{
  int Deg = P1.Degree - P2.Degree;
  double* Coef1 = new double [P1.Degree];
  for (int i = P1.Degree - 1; i >= 0; i--)
    Coef1[i] = P1.Coefficient[i];
  double* Coef2 = new double [Deg+1];
  double c = 1 / P2.Coefficient[P2.Degree];
  Coef2[Deg] = P1.Coefficient[P1.Degree] * c;  
  for (int i = Deg; i > 0; i--)
    {
      for (int j = P2.Degree - 1; j >= 0; j--)
	Coef1[i + j] -= Coef2[i] * P2.Coefficient[j];
      Coef2[i-1] = c * Coef1[P2.Degree + i - 1]; 
    }
  for (int j = P2.Degree - 1; j >= 0; j--)
    Coef1[j] -= Coef2[0] * P2.Coefficient[j];
  return Polynomial(P2.Degree - 1, Coef1, true);
}

Polynomial& Polynomial::operator += (const Polynomial& P)
{
  if ((this->RootFlag == true) && (this->NbrRoot != 0))
    {
      delete this->Root;
      this->RootFlag = false;
    }
  if (P.Degree < this->Degree)
    {
      for (int i = 0; i <= P.Degree; i++)
	this->Coefficient[i] += P.Coefficient[i];
      return *this;      
    }
  if (P.Degree > this->Degree) 
    {
      double* TmpCoef = new double [P.Degree];
      for (int i = 0; i <= this->Degree; i++)
	TmpCoef[i] += P.Coefficient[i];
      for (int i = this->Degree + 1; i <= P.Degree; i++)
	TmpCoef[i] = P.Coefficient[i];
      delete this->Coefficient;
      this->Coefficient = TmpCoef;
      this->Degree = P.Degree;
      return *this;      
    }
  else
    {
      double* TmpCoef = new double [this->Degree];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] += P.Coefficient[i];
      int i = this->Degree;
      while ((i > 0) && (TmpCoef[i] < EPSILON) && (TmpCoef[i] > -EPSILON))
	i--;
      if (i != this->Degree)
	{
	  double* TmpCoef = new double [i+1];
	  this->Degree = i;
	  for (int j = i; j >= 0; j--)
	    TmpCoef[j] = this->Coefficient [i]; 
	  delete this->Coefficient;
	  this->Coefficient = TmpCoef;
	}
      return *this;      
    }
}

Polynomial& Polynomial::operator -= (const Polynomial& P)
{
  if ((this->RootFlag == true) && (this->NbrRoot != 0))
    {
      delete this->Root;
      this->RootFlag = false;
    }
  if (P.Degree < this->Degree)
    {
      for (int i = 0; i <= P.Degree; i++)
	this->Coefficient[i] -= P.Coefficient[i];
      return *this;      
    }
  if (P.Degree > this->Degree) 
    {
      double* TmpCoef = new double [P.Degree];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] -= P.Coefficient[i];
      for (int i = this->Degree + 1; i <= P.Degree; i++)
	TmpCoef[i] = -P.Coefficient[i];
      delete this->Coefficient;
      this->Coefficient = TmpCoef;
      this->Degree = P.Degree;
      return *this;      
    }
  else
    {
      double* TmpCoef = new double [this->Degree];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] -= P.Coefficient[i];
      int i = this->Degree;
      while ((i > 0) && (TmpCoef[i] < EPSILON) && (TmpCoef[i] > -EPSILON))
	i--;
      if (i != this->Degree)
	{
	  double* TmpCoef = new double [i+1];
	  this->Degree = i;
	  for (int j = i; j >= 0; j--)
	    TmpCoef[j] = this->Coefficient [i]; 
	  delete this->Coefficient;
	  this->Coefficient = TmpCoef;
	}
      return *this;      
    }
}

Polynomial& Polynomial::operator *= (const double& d)
{
  for (int i = 0; i <= this->Degree; i++)
    this->Coefficient[i] *= d;
  return *this;
}

Polynomial& Polynomial::operator *= (const Polynomial& P)
{
  int Deg = P.Degree + this->Degree;
  double* Coef = new double [Deg + 1];
  if (P.Degree == 0)
    {
      delete this->Coefficient;
      this->Coefficient = Coef;
      return *this;
    }
  for (int i = this->Degree; i >= 0; i--)
    Coef[i + P.Degree] = this->Coefficient[i] * P.Coefficient[P.Degree];
  for (int i = P.Degree-1; i >= 0; i--)
    Coef[i] = 0;
  for (int i = P.Degree - 1; i >= 0; i--)
    for (int j = this->Degree; j >= 0; j--)
      Coef[i+j] += this->Coefficient[j] * P.Coefficient[i];
  delete this->Coefficient;
  this->Coefficient = Coef;
  if (this->RootFlag == true)
    {
      delete this->Root;
      this->RootFlag = false;
    }
  return *this;
}

Polynomial& Polynomial::operator /= (const double& d)
{
  for (int i = 0; i <= this->Degree; i++)
    this->Coefficient[i] /= d;
  return *this;
}

Polynomial& Polynomial::operator /= (const Polynomial& P)
{
  int Deg = this->Degree - P.Degree;
  double* Coef2 = new double [Deg+1];
  double c = 1 / P.Coefficient[P.Degree];
  Coef2[Deg] = this->Coefficient[this->Degree] * c;  
  for (int i = Deg; i > 0; i--)
    {
      for (int j = P.Degree - 1; j >= 0; j--)
	this->Coefficient[i + j] -= Coef2[i] * P.Coefficient[j];
      Coef2[i-1] = c * this->Coefficient[P.Degree + i - 1]; 
    }
  delete[] this->Coefficient;
  if (this->RootFlag == true)
    {
      delete this->Root;
      this->RootFlag = false;
    }
  this->Coefficient = Coef2;
  this->Degree = Deg;
  return *this;
}

// Divide polynomial by a monomial (x - x0) where x0 is a root of the given polynomial using Horner scheme
//
// x0 = root
// return value = result of polynomial division

Polynomial Polynomial::MonomialDivision (double x0)
{
  double* Coef1 = new double [this->Degree];
  Coef1[this->Degree - 1] = this->Coefficient[this->Degree];
  for (int i = this->Degree - 1; i > 0; i--)
    Coef1[i - 1] = Coef1[i] * x0 + this->Coefficient[i];
  return Polynomial(this->Degree - 1, Coef1, true);  
}


// Output Stream overload
//

ostream& operator << (ostream& Str, const Polynomial& P)
{
  if ((P.RootFlag == true) && (P.NbrRoot != 0))
    {
      Str << P.Root[P.Degree];
      for (int i = 0; i < P.Degree; i++)
	Str << "(x - " << P.Root[i] << ")";
      return Str;  
    }
  else
    {
      for (int i = P.Degree; i > 0; i--)
	Str << P.Coefficient[i] << " x^" << i << " + ";
      Str << P.Coefficient[0];
/*      Str << "[";
      for (int i = 0; i < P.Degree; i++)
	Str << P.Coefficient[i] << ",";
	  Str << P.Coefficient[P.Degree] << "]";*/
      return Str;  
    }
}

// Find Roots of the polynomial
//

void Polynomial::SolvePolynomial()
{
  switch(this->Degree)
    {
    case 1:
      this->SolveLinear();
      break;
    case 2:
      this->SolveQuadratic();
      break;
    case 3:
      this->SolveCubic();
      break;
      return;
    case 4:
      this->SolveQuartic();
      break;
    default:
      this->SolveHigherOrder();
      break;
    }
}
  
// Find Root of a linear polynomial
//

void Polynomial::SolveLinear()
{
  this->RootFlag = true;
  this->NbrRoot = 1;
  this->Root = new double [1];
  Root[0] = - this->Coefficient[0] / this->Coefficient[1];
}      

// Find Roots of a quadratic polynomial
//

void Polynomial::SolveQuadratic()
{
  this->RootFlag = true;
  double Delta = this->Coefficient[1] * this->Coefficient[1] - 4 * this->Coefficient[0] * this->Coefficient[2];
  if (Delta >= 0)
    {
      this->NbrRoot = 2;
      this->Root = new double [2];
      Delta = sqrt(Delta);
      double Inv = 1 / this->Coefficient[2];
      this->Root[0] = -0.5 * (this->Coefficient[1] + Delta) * Inv;
      this->Root[1] = this->Root[0] + Delta * Inv;
    }
  else
    this->NbrRoot = 0;
}      

// Find Roots of a cubic polynomial
//

void Polynomial::SolveCubic()
{
  this->RootFlag = true;
  long double a = this->Coefficient[2];
  long double b = this->Coefficient[1] ;
  long double a2 = a * a;
  long double q = MC1_9 * a2 - MC1_3 * b;
  long double r = 0.5 * (this->Coefficient[0] - a * b * MC1_3) + a * a2 * MC1_27;
  long double Delta = q * q * q - r * r;
  long double Inc = a * MC1_3;
  if (Delta >= 0)
    {
      this->NbrRoot = 3;
      this->Root = new double [3];
      long double sqq = sqrt (q);
      long double alpha = acos(- r / (q * sqq)) * MC1_3;      
      sqq *= 2;
      long double cosalpha = cos(alpha);
      long double sinalpha = sqq * sin(alpha);
      if (cosalpha >= 0.5)
	{
	  cosalpha *= sqq;
	  this->Root[2] = cosalpha - Inc;
	  this->Root[0] = cosalpha  * cos2PI_3 - sinalpha * sin2PI_3 - Inc;
	  this->Root[1] = cosalpha  * cos2PI_3 + sinalpha * sin2PI_3  - Inc;
	}
      else
	if (cosalpha >= -0.5)
	  {
	    cosalpha *= sqq;
	    this->Root[1] = cosalpha - Inc;
	    this->Root[0] = cosalpha  * cos2PI_3 - sinalpha * sin2PI_3 - Inc;
	    this->Root[2] = cosalpha  * cos2PI_3 + sinalpha * sin2PI_3  - Inc;
	  }
	else
	  {
	    cosalpha *= sqq;
	    this->Root[0] = cosalpha - Inc;
	    this->Root[1] = cosalpha  * cos2PI_3 - sinalpha * sin2PI_3 - Inc;
	    this->Root[2] = cosalpha  * cos2PI_3 + sinalpha * sin2PI_3  - Inc;
	  }
    }
  else
    {
      this->NbrRoot = 1;
      this->Root = new double [1];
      if (r > 0)
	{
	  long double expon = pow(r + sqrt (-Delta), MC1_3);
	  this->Root[0] = - (expon + q / expon) - Inc;
	}
      else
	{
	  long double expon = pow(- r + sqrt (-Delta), MC1_3);
	  this->Root[0] = (expon + q / expon) - Inc;
	}
    }	  
}      

// Find Roots of a quartic polynomial
//

void Polynomial::SolveQuartic()
{
  this->RootFlag = true;
  // solve resolvant
  long double a = - this->Coefficient[2];
  long double b = this->Coefficient[3] * this->Coefficient[1] - 4 * this->Coefficient[0];
  long double c = 4 * this->Coefficient[2] * this->Coefficient[0] - this->Coefficient[3] * this->Coefficient[3] * this->Coefficient[0]
    - this->Coefficient[1] * this->Coefficient[1];
  long double a2 = a * a;
  long double q = MC1_9 * a2 - MC1_3 * b;
  long double r = 0.5 * (c - a * b * MC1_3) + a * a2 * MC1_27;
  long double delta = q * q *q - r * r;
  long double inc = a * MC1_3;
  long double root3;
  if (delta >= 0)
    {
      long double sqq = sqrt (q);
      long double alpha = acos(- r / (q * sqq)) * MC1_3;      
      long double cosalpha = cos(alpha);
      sqq *= 2;
      root3 = sqq * cosalpha - inc;
    }
  else
    {
      if (r > 0)
	{
	  long double expon = pow(r + sqrt (-delta), MC1_3);
	  root3 = - (expon + q / expon) - inc;
	}
      else
	{
	  long double expon = pow(- r + sqrt (-delta), MC1_3);
	  root3 = (expon + q / expon) - inc;
	}
    }	  
  long double Fac = this->Coefficient[3] * this->Coefficient[3] * 0.25 + root3 - this->Coefficient[2];
  if (Fac <= -POLYERROR)
    {
      this->NbrRoot = 0;
      return;
    }
  long double d0, d1, e0, e1;
  if ((Fac > -POLYERROR) && (Fac < POLYERROR))
    {
      long double Fac2 = 0.25 * root3 * root3 - this->Coefficient[0];
      if (Fac2 < -POLYERROR)
	{
	  this->NbrRoot = 0;
	  return;
	}
      if (Fac2 < POLYERROR)
	Fac2 = 0;
      else
	Fac2 = sqrt(Fac2);
      d1 = 0.5 * this->Coefficient[3];
      d0 = 0.5 * root3 + Fac2;
      e1 = d1;
      e0 = 0.5 * root3 - Fac2;
    }
  else
    {
      Fac = sqrt(Fac);
      long double Fac2 = 0.5 * (0.5 * this->Coefficient[3] * root3 - this->Coefficient[1]) / Fac;
      d1 = 0.5 * this->Coefficient[3] - Fac;
      d0 = 0.5 * root3 - Fac2;
      e1 = 0.5 * this->Coefficient[3] + Fac;
      e0 = 0.5 * root3 + Fac2;
    }
  long double delta1 = d1 * d1 - 4 * d0;
  long double delta2 = e1 * e1 - 4 * e0;
  if ((delta1 < -POLYERROR) && (delta2 < -POLYERROR))
    {
      this->NbrRoot = 0;
      return;
    }
  if ((delta1 < -POLYERROR) && (delta2 >= -POLYERROR))
    {
      if (delta2 < POLYERROR)
	delta2 = 0;
      else
	delta2 = sqrt(delta2);	  
      this->NbrRoot = 2;
      this->Root = new double [2];
      this->Root[0] = - 0.5 * (delta2 + e1);
      this->Root[1] = this->Root[0] + delta2;
      return;
    }
  if ((delta2 < -POLYERROR) && (delta1 >= -POLYERROR))
    {
      if (delta1 < POLYERROR)
	delta1 = 0;
      else
	delta1 = sqrt(delta1);	  
      this->NbrRoot = 2;
      this->Root = new double [2];
      this->Root[0] = - 0.5 * (delta1 + d1);
      this->Root[1] = this->Root[0] + delta1;
      return;
    }
  if (delta1 < POLYERROR)
    delta1 = 0;
  else
    delta1 = sqrt(delta1);	  
  if (delta2 < POLYERROR)
    delta2 = 0;
  else
    delta2 = sqrt(delta2);	  
  this->NbrRoot = 4;
  this->Root = new double [4];
  double x1 = - 0.5 * (delta1 + d1);
  double x2 = - 0.5 * (delta2 + e1);
  if (x1 < x2)
    {
      this->Root[0] = x1;
      x1 += delta1;
      if (x1 < x2)
	{
	  this->Root[1] = x1;
	  this->Root[2] = x2;
	  this->Root[3] = x2 + delta2;   
	  return;
	}
      this->Root[1] = x2;
      x2 += delta2;
      if (x1 < x2)
	{
	  this->Root[2] = x1;
	  this->Root[3] = x2;   
	  return;
	}
      else
	{
	  this->Root[2] = x2;
	  this->Root[3] = x1;   
	  return;
	}
    }   
  else
    {
      this->Root[0] = x2;
      x2 += delta2;
      if (x2 < x1)
	{
	  this->Root[1] = x2;
	  this->Root[2] = x1;
	  this->Root[3] = x1 + delta1;   
	  return;
	}
      this->Root[1] = x1;
      x1 += delta1;
      if (x1 < x2)
	{
	  this->Root[2] = x1;
	  this->Root[3] = x2;   
	  return;
	}
      else
	{
	  this->Root[2] = x2;
	  this->Root[3] = x1;   
	  return;
	}
    }   
}

// Find Roots of a polynomial of degree > 4
//

void Polynomial::SolveHigherOrder()
{
  bool FlagEven = true;
 


  this->RootFlag = true;
  Polynomial DP = this->DerivatePolynomial();
  if (this->Degree > 4)
    DP.SolveHigherOrder();
  else
    DP.SolveQuartic();
  double* TmpRoot = new double [this->Degree];
  int NR = 0;
  if (DP.NbrRoot == 0)
    {
    }
  else
    {
      double f1, f2;
      if (FlagEven == true)
	
      if (DP.NbrRoot > 1)
	{
	  for (int i = 0; i < (DP.NbrRoot - 1); i++)
	    {
	      f1 = this->PolynomialEvaluate(DP.Root[i]);
	      f2 = this->PolynomialEvaluate(DP.Root[i+1]);
	      if ((f1 < f2) && (f1 <= 0) && (f2 >= 0))
		{
		  TmpRoot[NR] = this->FindRootBottomUp(Root[i], Root[i+1], (f1 - f2) * POLYERROR, POLYERROR);
		  NR++;
		}
	      else
		if ((f1 > f2) && (f1 >= 0) && (f2 <= 0))
		  {
		    TmpRoot[NR] = this->FindRootTopDown(Root[i], Root[i+1], (f1 - f2) * POLYERROR, POLYERROR);
		    NR++;
		  }
	      
	    }

	}
    }
}

// Refine root value 

void Polynomial::RefineRoot (double Epsilon, int MaxIter)
{
  if (this->NbrRoot == 0)
    return;
  double x, d, e;
  int k;
  for (int i = 0; i < this->NbrRoot; i++)
    {
      e = this->Root[i]; 
      k = 0;
      x = this->PolynomialEvaluate(e);
      while (((x> Epsilon) || (x < -Epsilon)) && (k++ < MaxIter))
	{
	  d = this->DerivativeEvaluate(e);
	  e = e - x / d;
	  x = this->PolynomialEvaluate(e);
	}
      if ((x> Epsilon) || (x < -Epsilon))
	if (i == (this->NbrRoot -1))
	  this->NbrRoot--;
	else
	  {
	    for (int j = i; j < (this->NbrRoot-1); j++)
	      this->Root[j] = this->Root[j+1];
	    this->NbrRoot--;
	    i--;
	  }
      else
	this->Root[i] = e;
    }
}

// Find a root given local minimun and local maximum (with P(Xmin) < P(Xmax)) using bisection method
//
// xmax = interval upper bound
// xmin = interval lower bound
// err = relative error on the root
// absErr = absolute error on the root
// return value = root belonging to the interval [xmin, xmax]

double Polynomial::FindRootBottomUp (double xmin, double xmax, double err, double absErr)
{
  double x = (xmin + xmax) * 0.5;
  if ((xmax - xmin) < err * Abs(xmax) + absErr)
    return x;
  double ErrBound = err;
  double P = this->PolynomialEvaluate(x);
  if ((P < absErr) && (P > -absErr))
    return x;
  if (P > 0)
    return this->FindRootTopDown(x, xmax, err, absErr);
  else 
    return this->FindRootTopDown(xmin, x, err, absErr);    
}
 
// Find a root given local minimun and local maximum (with P(Xmin) > P(Xmax)) using bisection method
//
// xmax = interval upper bound
// xmin = interval lower bound
// err = relative error on the root
// absErr = absolute error on the root
// return value = root belonging to the interval [xmin, xmax]

double Polynomial::FindRootTopDown (double xmax, double xmin, double err, double absErr)
{
  double x = (xmin + xmax) * 0.5;
  if ((xmax - xmin) < err * Abs(xmax) + absErr)
    return x;
  double ErrBound = err;
  double P = this->PolynomialEvaluate(x);
  if ((P < absErr) && (P > -absErr))
    return x;
  if (P > 0)
    return this->FindRootTopDown(x, xmax, err, absErr);
  else 
    return this->FindRootTopDown(xmin, x, err, absErr);    
}

// use Laguerre Method to find roots
//
// err = relative error on the root
// absErr = absolute error on the root
// nbrIteration = maximum iteration number to find a root

void Polynomial::LaguerreMethodSolver (double err, double absErr, int nbrIteration)
{
  int NbrComplexRoot = 0;  
  int NbrRealRoot = 0;
  double* TmpRoot = new double [this->Degree];
  Polynomial P = *this;
  while (NbrComplexRoot < this->Degree)
    {
      if (P.Degree == 1)
	{
	  TmpRoot[NbrRealRoot++] = - P.Coefficient[0] / P.Coefficient[1];
	  cout << "root = " << TmpRoot[NbrRealRoot - 1] << " " << this->PolynomialEvaluate(TmpRoot[NbrRealRoot - 1]) << endl;
	  NbrComplexRoot += 2;	  
	}
      else
	if (P.Degree == 2)
	  {
	    P.SolveQuadratic();
	    if (P.RootFlag == true)
	      {
	      TmpRoot[NbrRealRoot++] = P.Root[0];
	      cout << "root = " << TmpRoot[NbrRealRoot - 1] << " " << this->PolynomialEvaluate(TmpRoot[NbrRealRoot - 1]) << endl;
	      TmpRoot[NbrRealRoot++] = P.Root[1];
	      cout << "root = " << TmpRoot[NbrRealRoot - 1] << " " << this->PolynomialEvaluate(TmpRoot[NbrRealRoot - 1]) << endl;
	      NbrComplexRoot += 2;	  
	      }
	  }
	else
	  {
	    Complex Root = P.LaguerreMethod(err, absErr, nbrIteration, Complex());
	    if (P.NoRootFlag == true)
	      {
		this->NoRootFlag = true;
		return;
	      }
	    if (Abs (Root.Im) < (Norm (Root) * err + absErr))
	      {	      
		cout << "root = " << Root.Re << " " << this->PolynomialEvaluate(Root.Re) << endl;
		cout << Abs(2.0 * Root.Re * err * this->DerivativeEvaluate(Root.Re)) << endl;
		if (Abs(this->PolynomialEvaluate(Root.Re)) > absErr)
		  {
		    cout << "error root = " << Root << endl;
		    /*		  if ((this->DerivativeEvaluate(Root.Re) * this->PolynomialEvaluate(Root.Re)) > 0)
				  Root.Re *= (1.0 + 1e-12);//this->PolynomialEvaluate(Root.Re) / this->DerivativeEvaluate(Root.Re);
				  else
				  Root.Re *= (1.0 - 1e-12);
				  //		  Root = P.LaguerreMethod(err, absErr, nbrIteration, Complex(Root.Re * 100.00, 0.0));
				  Root = P.BauhuberMethod(err, absErr, nbrIteration, Root);
				  cout << "new root = " << Root << " " << this->PolynomialEvaluate(Root.Re) << endl;
				  cout << "end error root" << endl;*/
		  }
		TmpRoot[NbrRealRoot++] = Root.Re;	  
		P = P.MonomialDivision(Root.Re);	  
		NbrComplexRoot++;
	      }
	    else
	      {
		cout << "root = " << Root << " " << this->PolynomialEvaluate(Root) << endl;;
		cout << "root = " << Conj(Root) << " " << this->PolynomialEvaluate(Conj(Root)) << endl;;
		//	      TmpRoot[NbrRealRoot++] = Root.Re;	  
		double* Coef = new double[3];
		Coef[0] = Root.Re * Root.Re + Root.Im * Root.Im;
		Coef[1] =  -2.0 * Root.Re;
		Coef[2] = 1.0;
		Polynomial TmpP (2, Coef, true);
		P /= TmpP;
		NbrComplexRoot += 2;
	      }
	  }
    }
  this->NbrRoot = NbrRealRoot;
  this->RootFlag = true;
  if (this->NbrRoot > 0)
    {
      this->Root = new double [this->NbrRoot];
      for (int i = 0; i < this->NbrRoot; i++)
	this->Root[i] = TmpRoot[i];
      this->SortRoots();
      if ((TestOdd(this->Degree) == true) && (this->Coefficient[this->Degree] > 0))
	cout << "ok" << endl;
      cout << "Real roots : ";
      for (int i = 0; i < this->NbrRoot; i++)
	cout << this->Root[i] << " ";
      cout << endl;      
    }
  else
    {
      this->Root = 0;
    }
  delete[] TmpRoot;
  return;
}

// use Laguerre Method to find a root (put NoRootFlag to true if no root was find)
//
// nbrIteration = maximum iteration number to find a root
// err = relative error on the root
// absErr = absolute error on the root
// return value = root obtained

double Polynomial::LaguerreMethod (double err, double absErr, int nbrIteration, double x0)
{
  if (nbrIteration == 0)
    {
      this->NoRootFlag = true;
      return x0;
    }
  double nP = this->PolynomialEvaluate(x0) * (double) this->Degree;
  double DP = this->DerivativeEvaluate (x0);
  double tmpS = (this->Degree - 1) * ((this->Degree - 1) * DP * DP - nP * this->DerivativeEvaluate (x0, 2));
  Complex S;
  if (tmpS < 0)
    S = Complex (0.0, sqrt(-tmpS));
  else
    S = sqrt(tmpS);    
  Complex x = x0;
  Complex tmp = DP + S;
  Complex tmp2 = DP - S;
  if (Norm(tmp2) >= Norm(tmp))
    x -= nP / tmp2;
  else
    x -= nP / tmp;
  if (Norm(x - x0) < (err * Norm(x) + absErr))
    {
      return x.Re;
    }
  else
    return this->LaguerreMethod(err, absErr, nbrIteration - 1, x).Re;
}

// use Laguerre Method to find a root (put NoRootFlag to true if no root was find)
//
// nbrIteration = maximum iteration number to find a root
// err = relative error on the root
// absErr = absolute error on the root
// return value = root obtained

Complex Polynomial::LaguerreMethod (double err, double absErr, int nbrIteration, Complex x0)
{
  if (nbrIteration == 0)
    {
      this->NoRootFlag = true;
      return x0;
    }
  Complex nP = this->PolynomialEvaluate(x0) * (double) this->Degree;
  Complex DP = this->DerivativeEvaluate (x0);
  Complex S = sqrt( (this->Degree - 1) * ((this->Degree - 1) * DP * DP - nP * this->DerivativeEvaluate (x0, 2)));
  Complex x = x0;
  Complex tmp = DP + S;
  DP -= S;
  if (Norm(DP) >= Norm(tmp))
    x -= nP / DP;
  else
    x -= nP / tmp;
  if (Norm(x - x0) < (err * Norm(x) + absErr))
    {
      return x;
    }
  else
    return this->LaguerreMethod(err, absErr, nbrIteration - 1, x);
}

/*Polynomial* Polynomial::SturmSequence();
{
  
}*/

// return absolute value of a given number
//
// x = number to proceed
// return value = absolute value of x

double Abs (double x)
{
  if (x < 0)
    return -x;
  else
    return x;
}

// use Bauhuber Method to find a root (put NoRootFlag to true if no root was find)
//
// nbrIteration = maximum iteration number to find a root
// err = relative error on the root
// absErr = absolute error on the root
// return value = root obtained

Complex Polynomial::BauhuberMethod (double err, double absErr, int nbrIteration, Complex x0)
{
  Complex x1;
  Complex P = this->PolynomialEvaluate(x0);
  Complex DP;
  Complex D2P;
  double AbsP0 = Norm(P);
  double AbsP1;
  while ((nbrIteration > 0) && (AbsP0 > absErr))
    {
      DP = this->DerivativeEvaluate(x0);
      D2P = this->DerivativeEvaluate(x0, 2) / DP;
      P /= DP;
      x1 = x0 - (P / (1.0 - P * D2P));
      P = this->PolynomialEvaluate(x1);
      AbsP1 = Norm(P);
      while ((AbsP1 > AbsP0) && (nbrIteration > 0))
	{
	  // spiralization
	  x1 *= Complex (0.1, 0.9);
	  x1 += x0 * Complex (0.9, -0.9);
	  P = this->PolynomialEvaluate(x1);	  
	  AbsP1 = Norm(P);
	  nbrIteration--;
	}
      nbrIteration--;
    }
  if (nbrIteration == 0)
    this->NoRootFlag = true;
  return x0;
}

// Sort roots of a polynomial
//

void Polynomial::SortRoots ()
{
  if (this->NbrRoot <= 1)
    return;
  double tmp;
  for (int i = this->NbrRoot - 2; i >= 0; i--)
    for (int j = 0; j <= i; j++)
      if (this->Root[j] > this->Root[j+1])
	{
	  tmp = this->Root[j];
	  this->Root[j] = this->Root[j+1];
	  this->Root[j+1] = tmp;
	}
  return;
}

// test if an integer is odd
//
// x = integer to test
// return value = true if x is odd

double TestOdd (int  x)
{
  return (((x >> 1) << 1) != x);
}

