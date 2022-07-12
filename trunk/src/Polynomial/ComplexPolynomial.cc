////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                           RayTrace version  0.10                           //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   Class of polynomial with complex coefficients            //
//                                                                            //
//                        last modification : 23/03/2005                      //
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


#include "Polynomial/ComplexPolynomial.h"
#include <cmath>


#pragma optimize("", off)

#define MY_EPSILON 1E-6

#include <iostream>
using std::cout;
using std::endl;


// default constructor
//

ComplexPolynomial::ComplexPolynomial ()
{
  this->Degree = 0;
  this->NoRootFlag = false;
  this->Coefficient = 0;
  this->RootFlag = false; 
}
 
// constructor
//
// degree = polynomial degree
// coefficients = coefficients array ( first element is associated to the -power term)
// flag = true if coefficients array has to be used directly and not duplicated

ComplexPolynomial::ComplexPolynomial (int degree, Complex* coefficients, bool flag)
{
  this->Degree = degree;
  this->NoRootFlag = false;
  if (flag == true)
    this->Coefficient = coefficients;
  else
    {
      this->Coefficient = new Complex [this->Degree + 1];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] = coefficients[i];
    }
  this->RootFlag = false;
}

ComplexPolynomial::ComplexPolynomial (const ComplexPolynomial& P)
{
  if (P.Coefficient != 0)
    {
      this->Degree = P.Degree;
      this->Coefficient = new Complex [this->Degree + 1];
      this->RootFlag = P.RootFlag;
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] = P.Coefficient[i];
      this->NoRootFlag = P.NoRootFlag;
      if (P.RootFlag == true)
	{
	  this->NbrRoot = P.NbrRoot;
	  this->Root = new Complex [this->NbrRoot];
	  for (int i = 0; i < this->NbrRoot; i++)
	    this->Root[i] = P.Root[i];
	}
    }
  else
    {
      this->Degree = 0;
      this->NoRootFlag = false;
      this->Coefficient = 0;
      this->RootFlag = false; 
    }
}

// destructor
//

ComplexPolynomial::~ComplexPolynomial()
{
  if ((this->RootFlag == true) && (this->NbrRoot != 0))
    delete[] this->Root;
  if (this->Coefficient != 0)
    delete[] this->Coefficient;
}

// assignement

ComplexPolynomial& ComplexPolynomial::operator = (const ComplexPolynomial& P)
{
  if (this->Coefficient != 0)
    delete[] this->Coefficient;
  if ((this->RootFlag == true) && (this->NbrRoot != 0))
    delete[] this->Root;
  this->Degree = P.Degree;
  this->RootFlag = P.RootFlag;
  if (P.Coefficient != 0)
    {
      this->Coefficient = new Complex [this->Degree + 1];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] = P.Coefficient[i];
    }
  else
    {
      this->Coefficient = 0;
    }
  if (P.RootFlag == true)
    {
      this->NbrRoot = P.NbrRoot;
      this->Root = new Complex [this->NbrRoot];
      for (int i = 0; i < this->NbrRoot; i++)
	this->Root[i] = P.Root[i];
    }
  return *this;
}


// screen PolynomialEvaluate from compilation problems in optimization step by intel compilers
#ifdef __INTEL_COMPILER
#pragma optimize("", off)
#endif

// Return polynomial value at a given point
//
// x = point where to evaluate polynomial
// return value = polynomial value at x

Complex ComplexPolynomial::PolynomialEvaluate (double x)
{
  if (this->Coefficient != 0)
    {
      Complex Res (this->Coefficient[this->Degree]);
      for (int i = this->Degree - 1 ; i >= 0; i--)
	{
	  Res *= x;
	  Res += this->Coefficient[i];
	}
      return Res;
    }
  return Complex();
}
  
// Return polynomial value at a given point
//
// x = point where to evaluate polynomial
// return value = polynomial value at x

Complex ComplexPolynomial::PolynomialEvaluate (Complex x)
{
  if (this->Coefficient != 0)
    {
      Complex Res (this->Coefficient[this->Degree]);
      for (int i = this->Degree - 1 ; i >= 0; i--)
	{
	  Res *= x;
	  Res += this->Coefficient[i];
	}
      return Res;
    }
  return Complex();
}

#ifdef __INTEL_COMPILER
#pragma optimize("", on)
#endif

// Evaluate polynomial derivative 
//
// x = position where to evaluate polynomial derivative
// return value = polynomial derivative at x

Complex ComplexPolynomial::DerivativeEvaluate (double x)
{
  if (this->Coefficient != 0)
    {
      Complex Res = this->Degree * this->Coefficient[this->Degree];
      for (int i = this->Degree - 1 ; i > 0; i--)
	{
	  Res *= x;
	  Res += i * this->Coefficient[i];
	}
      return Res;
    }
  return Complex();
}
  
// Evaluate polynomial derivative 
//
// x = position where to evaluate polynomial derivative
// return value = polynomial derivative at x

Complex ComplexPolynomial::DerivativeEvaluate (Complex x)
{
  if (this->Coefficient != 0)
    {
      Complex Res = this->Degree * this->Coefficient[this->Degree];
      for (int i = this->Degree - 1 ; i > 0; i--)
	{
	  Res *= x;
	  Res += i * this->Coefficient[i];
	}
      return Res;
    }
  return Complex();
}
  
// Evaluate polynomial n-th derivative at a given value
//
// x = position where to evaluate polynomial n-th derivative
// n = derivative order
// return value = polynomial n-th derivative at x

Complex ComplexPolynomial::DerivativeEvaluate (double x, int n)
{
  if (this->Coefficient != 0)
    {
      if (n > this->Degree)
	return 0.0;
      Complex tmp = 1.0;
      for (int j = this->Degree - n + 1; j <= this->Degree; j++)
	tmp *= (double) j;
      Complex Res = tmp * this->Coefficient[this->Degree];
      for (int i = this->Degree - 1; i >= n; i--)
	{
	  tmp = 1.0;
	  for (int j = i - n + 1; j <= i; j++)
	    tmp *= (double) j;      
	  Res *= x;
	  Res += tmp * this->Coefficient[i];
	}
      return Res;
    }
  return Complex();
}
  
// Evaluate polynomial n-th derivative at a given value
//
// x = position where to evaluate polynomial n-th derivative
// n = derivative order
// return value = polynomial n-th derivative at x

Complex ComplexPolynomial::DerivativeEvaluate (Complex x, int n)
{
  if (this->Coefficient != 0)
    {
      if (n > this->Degree)
	return Complex();
      Complex tmp (1.0);
      for (int j = this->Degree - n + 1; j <= this->Degree; j++)
	tmp *= (double) j;
      Complex Res (tmp * this->Coefficient[this->Degree]);
      for (int i = this->Degree - 1; i >= n; i--)
	{
	  tmp = 1.0;
	  for (int j = i - n + 1; j <= i; j++)
	    tmp *= (double) j;      
	  Res *= x;
	  Res += tmp * this->Coefficient[i];
	}
      return Res;
    }
  return Complex();
}
  
// Return Derivative of the polynomial 
//
// return value = polynomial derivative

ComplexPolynomial ComplexPolynomial::DerivatePolynomial ()
{
  if (this->Coefficient != 0)
    {
      Complex* Coef = new Complex [this->Degree];
      for (int i = 1; i <= this->Degree; i++)
	Coef[i-1] = (this->Coefficient[i] * i);
      ComplexPolynomial P (this->Degree - 1, Coef);
      delete[] Coef;
      return P;
    }
  return ComplexPolynomial();
}

// Evaluate polynomial n-th derivative
//
// n = derivative order
// return value = polynomial n-th derivative

ComplexPolynomial ComplexPolynomial::DerivatePolynomial (int n)
{
  if (this->Coefficient != 0)
    {
      if (n > this->Degree)
	{
	  Complex* Coef = new Complex [1];
	  Coef[0] = 0.0;
	  return ComplexPolynomial(0, Coef, true);
	}
      int tmpDegree = this->Degree - n;
      Complex* Coef = new Complex [tmpDegree + 1];  
      Complex tmp (1.0);
      for (int j = this->Degree - n + 1; j <= this->Degree; j++)
	tmp *= (double) j;
      Coef[tmpDegree] = tmp * this->Coefficient[this->Degree];
      for (int i = this->Degree - 1; i >= n; i--)
	{
	  tmp = 1.0;
	  for (int j = i - n + 1; j <= i; j++)
	    tmp *= (double) j;      
	  Coef[i - n] = tmp * this->Coefficient[i];
	}
      return  ComplexPolynomial(tmpDegree, Coef, true);
    }
  return ComplexPolynomial();
}

// arithmetic operators

ComplexPolynomial operator - (const ComplexPolynomial& P)
{
  ComplexPolynomial Pr(P);
  if (Pr.Coefficient != 0)
    {
      for (int i = 0; i <= Pr.Degree; i++)
	Pr.Coefficient[i] = -Pr.Coefficient[i];  
    }
  return Pr;  
}

ComplexPolynomial operator + (const ComplexPolynomial& P1, const ComplexPolynomial& P2)
{
  if (P1.Degree > P2.Degree)    
    {
      Complex* Coef = new Complex [P1.Degree+1];
      for (int i = P1.Degree; i > P2.Degree; i--)
	Coef[i] = P1.Coefficient[i];
      for (int i = P2.Degree; i >= 0; i--)
	Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];	  
      return ComplexPolynomial(P1.Degree, Coef, true);
    }
  else
    if (P1.Degree < P2.Degree)    
      {
	Complex* Coef = new Complex [P2.Degree+1];
	for (int i = P2.Degree; i > P1.Degree; i--)
	  Coef[i] = P2.Coefficient[i];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];	  
	return ComplexPolynomial(P2.Degree, Coef, true);

      }
    else
      {
	Complex* Coef = new Complex [P1.Degree+1];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];
	int i = P1.Degree;
	while ((Coef[i] == 0) && (i!=0)) i--;
	if (i == P1.Degree)
	  return ComplexPolynomial(P2.Degree, Coef, true);
	else 
	  {
	    Complex* Coef2 = new Complex [i+1];
	    for (int j = 0; j <= i; j++)
	      Coef2[j] = Coef[j];
	    delete[] Coef;
	    return ComplexPolynomial(i, Coef2, true);
	  }
      }
}

ComplexPolynomial operator - (const ComplexPolynomial& P1, const ComplexPolynomial& P2)
{
  if (P1.Degree > P2.Degree)    
    {
      Complex* Coef = new Complex [P1.Degree+1];
      for (int i = P1.Degree; i > P2.Degree; i--)
	Coef[i] = P1.Coefficient[i];
      for (int i = P2.Degree; i >= 0; i--)
	Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];	  
      return ComplexPolynomial(P1.Degree, Coef, true);
    }
  else
    if (P1.Degree < P2.Degree)    
      {
	Complex* Coef = new Complex [P2.Degree+1];
	for (int i = P2.Degree; i > P1.Degree; i--)
	  Coef[i] = - P2.Coefficient[i];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];	  
	return ComplexPolynomial(P2.Degree, Coef, true);

      }
    else
      {
	Complex* Coef = new Complex [P1.Degree+1];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];
	int i = P1.Degree;
	while ((Coef[i] == 0) && (i!=0)) i--;
	if (i == P1.Degree)
	  return ComplexPolynomial(P2.Degree, Coef, true);
	else 
	  {
	    Complex* Coef2 = new Complex [i+1];
	    for (int j = 0; j <= i; j++)
	      Coef2[j] = Coef[j];
	    delete[] Coef;
	    return ComplexPolynomial(i, Coef2, true);
	  }
      }
}

ComplexPolynomial operator * (const ComplexPolynomial& P1, const ComplexPolynomial& P2)
{
  int Deg = P1.Degree+P2.Degree;
  Complex* Coef = new Complex [Deg + 1];
  for (int i = P1.Degree; i >= 0; i--)
    Coef[i+P2.Degree] = P1.Coefficient[i] * P2.Coefficient[P2.Degree];
  if (P2.Degree == 0)
    return ComplexPolynomial(Deg, Coef, true);
  for (int i = P2.Degree-1; i >= 0; i--)
    Coef[i] = 0;
  for (int i = P2.Degree - 1; i >= 0; i--)
    for (int j = P1.Degree; j >= 0; j--)
      Coef[i+j] += P1.Coefficient[j] * P2.Coefficient[i];
  return ComplexPolynomial(Deg, Coef, true);
}

ComplexPolynomial operator * (const ComplexPolynomial& P, const double& d)
{
  ComplexPolynomial Pr(P);
  if (P.Coefficient != 0)
    {
      for (int i = 0; i <= Pr.Degree; i++)
	Pr.Coefficient[i] *= d;
    }
  return Pr;
}

ComplexPolynomial operator * (const double& d, const ComplexPolynomial& P)
{
  ComplexPolynomial Pr (P);
  if (P.Coefficient != 0)
    {
      for (int i = 0; i <= Pr.Degree; i++)
	Pr.Coefficient[i] *= d;
    }
  return Pr;
}

ComplexPolynomial operator / (const ComplexPolynomial& P, const double& d)
{
  ComplexPolynomial Pr (P);
  if (Pr.Coefficient != 0)
    {
      for (int i = 0; i <= Pr.Degree; i++)
	Pr.Coefficient[i] /= d;
    }
  return Pr;
}

ComplexPolynomial operator / (const ComplexPolynomial& P1, const ComplexPolynomial& P2)
{
  int Deg = P1.Degree - P2.Degree;
  Complex* Coef1 = new Complex [P1.Degree];
  for (int i = P1.Degree-1; i >= 0; i--)
    Coef1[i] = P1.Coefficient[i];
  Complex* Coef2 = new Complex [Deg+1];
  Complex c = 1.0 / P2.Coefficient[P2.Degree];
  Coef2[Deg] = P1.Coefficient[P1.Degree] * c;  
  for (int i = Deg; i > 0; i--)
    {
      for (int j = P2.Degree - 1; j >= 0; j--)
	Coef1[i + j] -= Coef2[i] * P2.Coefficient[j];
      Coef2[i-1] = c * Coef1[P2.Degree + i - 1]; 
    }
  delete[] Coef1;
  return ComplexPolynomial(Deg, Coef2, true);
}

ComplexPolynomial operator % (const ComplexPolynomial& P1, const ComplexPolynomial& P2)
{
  int Deg = P1.Degree - P2.Degree;
  Complex* Coef1 = new Complex [P1.Degree];
  for (int i = P1.Degree - 1; i >= 0; i--)
    Coef1[i] = P1.Coefficient[i];
  Complex* Coef2 = new Complex [Deg+1];
  Complex c = 1.0 / P2.Coefficient[P2.Degree];
  Coef2[Deg] = P1.Coefficient[P1.Degree] * c;  
  for (int i = Deg; i > 0; i--)
    {
      for (int j = P2.Degree - 1; j >= 0; j--)
	Coef1[i + j] -= Coef2[i] * P2.Coefficient[j];
      Coef2[i-1] = c * Coef1[P2.Degree + i - 1]; 
    }
  for (int j = P2.Degree - 1; j >= 0; j--)
    Coef1[j] -= Coef2[0] * P2.Coefficient[j];
  return ComplexPolynomial(P2.Degree - 1, Coef1, true);
}

ComplexPolynomial& ComplexPolynomial::operator += (const ComplexPolynomial& P)
{
  if ((this->RootFlag == true) && (this->NbrRoot != 0))
    {
      delete[] this->Root;
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
      Complex* TmpCoef = new Complex [P.Degree];
      for (int i = 0; i <= this->Degree; i++)
	TmpCoef[i] += P.Coefficient[i];
      for (int i = this->Degree + 1; i <= P.Degree; i++)
	TmpCoef[i] = P.Coefficient[i];
      delete[] this->Coefficient;
      this->Coefficient = TmpCoef;
      this->Degree = P.Degree;
      return *this;      
    }
  else
    {
      Complex* TmpCoef = new Complex [this->Degree];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] += P.Coefficient[i];
      int i = this->Degree;
      while ((i > 0) && (Norm(TmpCoef[i]) < MY_EPSILON) && (Norm(TmpCoef[i]) > -MY_EPSILON))
	i--;
      if (i != this->Degree)
	{
	  Complex* TmpCoef = new Complex [i+1];
	  this->Degree = i;
	  for (int j = i; j >= 0; j--)
	    TmpCoef[j] = this->Coefficient [i]; 
	  delete[] this->Coefficient;
	  this->Coefficient = TmpCoef;
	}
      return *this;      
    }
}

ComplexPolynomial& ComplexPolynomial::operator -= (const ComplexPolynomial& P)
{
  if ((this->RootFlag == true) && (this->NbrRoot != 0))
    {
      delete[] this->Root;
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
      Complex* TmpCoef = new Complex [P.Degree];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] -= P.Coefficient[i];
      for (int i = this->Degree + 1; i <= P.Degree; i++)
	TmpCoef[i] = -P.Coefficient[i];
      delete[] this->Coefficient;
      this->Coefficient = TmpCoef;
      this->Degree = P.Degree;
      return *this;      
    }
  else
    {
      Complex* TmpCoef = new Complex [this->Degree];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] -= P.Coefficient[i];
      int i = this->Degree;
      while ((i > 0) && (Norm(TmpCoef[i]) < MY_EPSILON) && (Norm(TmpCoef[i]) > -MY_EPSILON))
	i--;
      if (i != this->Degree)
	{
	  Complex* TmpCoef = new Complex [i+1];
	  this->Degree = i;
	  for (int j = i; j >= 0; j--)
	    TmpCoef[j] = this->Coefficient [i]; 
	  delete[] this->Coefficient;
	  this->Coefficient = TmpCoef;
	}
      return *this;      
    }
}

ComplexPolynomial& ComplexPolynomial::operator *= (const double& d)
{
  if (this->Coefficient != 0)
    {
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] *= d;
    }
  return *this;
}

ComplexPolynomial& ComplexPolynomial::operator *= (const ComplexPolynomial& P)
{
  if (this->Coefficient != 0)
    {
      int Deg = P.Degree + this->Degree;
      Complex* Coef = new Complex [Deg + 1];
      if (P.Degree == 0)
	{
	  delete[] this->Coefficient;
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
      delete[] this->Coefficient;
      this->Coefficient = Coef;
      if (this->RootFlag == true)
	{
	  delete[] this->Root;
	  this->RootFlag = false;
	}
    }
  return *this;
}

ComplexPolynomial& ComplexPolynomial::operator /= (const double& d)
{
  if (this->Coefficient != 0)
    {
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] /= d;
    }
  return *this;
}

ComplexPolynomial& ComplexPolynomial::operator /= (const ComplexPolynomial& P)
{
  if (this->Coefficient != 0)
    {
      int Deg = this->Degree - P.Degree;
      Complex* Coef2 = new Complex [Deg+1];
      Complex c = 1.0 / P.Coefficient[P.Degree];
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
	  delete[] this->Root;
	  this->RootFlag = false;
	}
      this->Coefficient = Coef2;
      this->Degree = Deg;
    }
  return *this;
}

// Divide polynomial by a monomial (z - z0) where z0 is a root of the given polynomial using Horner scheme
//
// x0 = root
// return value = result of polynomial division

ComplexPolynomial ComplexPolynomial::MonomialDivision (const Complex& z0)
{
  if (this->Coefficient != 0)
    {
      Complex* Coef1 = new Complex [this->Degree];
      Coef1[this->Degree - 1] = this->Coefficient[this->Degree];
      for (int i = this->Degree - 1; i > 0; i--)
	Coef1[i - 1] = Coef1[i] * z0 + this->Coefficient[i];
      return ComplexPolynomial(this->Degree - 1, Coef1, true);  
    }
  return ComplexPolynomial();
}


// Output Stream overload
//

ostream& operator << (ostream& Str, const ComplexPolynomial& P)
{
  if (P.Coefficient != 0)
    {
      if ((P.RootFlag == true) && (P.NbrRoot != 0))
	{
	  Str << P.Root[P.Degree];
	  for (int i = 0; i < P.Degree; i++)
	    Str << "(z - " << P.Root[i] << ")";
	  return Str;  
	}
      else
	{
	  for (int i = P.Degree; i > 0; i--)
	    Str << P.Coefficient[i] << " z^" << i << " + ";
	  Str << P.Coefficient[0];
	  return Str;  
	}
    }
  return Str;  
}

// Find Roots of the polynomial
//

void ComplexPolynomial::SolvePolynomial()
{
  switch(this->Degree)
    {
    case 1:
      this->SolveLinear();
      break;
    }
}
  
// Find Root of a linear polynomial
//

void ComplexPolynomial::SolveLinear()
{
  this->RootFlag = true;
  this->NbrRoot = 1;
  this->Root = new Complex [1];
  Root[0] = - this->Coefficient[0] / this->Coefficient[1];
}      


// Sort roots of a polynomial
//

void ComplexPolynomial::SortRoots ()
{
  if (this->NbrRoot <= 1)
    return;
  Complex tmp;
  for (int i = this->NbrRoot - 2; i >= 0; i--)
    for (int j = 0; j <= i; j++)
      if (Norm(this->Root[j]) > Norm(this->Root[j+1]))
	{
	  tmp = this->Root[j];
	  this->Root[j] = this->Root[j+1];
	  this->Root[j+1] = tmp;
	}
  return;
}


#pragma optimize("", on)
