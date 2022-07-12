////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             Class of polynomial with long rational coefficients            //
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


#include "Polynomial/LongRationalPolynomial.h"
#include <math.h>


using std::cout;
using std::endl;


// default constructor
//

LongRationalPolynomial::LongRationalPolynomial ()
{
  this->Degree = 0;
  this->MaxDegree = 0;
  this->Coefficient = 0;
  this->RationalRoots = 0;
  this->RootFlag = false; 
}
 
// constructor
//
// degree = polynomial degree

LongRationalPolynomial::LongRationalPolynomial (int degree)
{
  this->Degree = degree;
  this->MaxDegree = degree;
  this->Coefficient = new LongRational [this->Degree + 1];
  for (int i = 0; i <= this->Degree; ++i)
    this->Coefficient[i] = 0l;
  this->RootFlag = false; 
  this->RationalRoots = 0;  
}

// constructor from raw datas
//
// degree = polynomial degree
// coefficients = coefficients array ( first element is associated to the -power term)
// flag = true if coefficients array has to be used directly and not duplicated

LongRationalPolynomial::LongRationalPolynomial (int degree, LongRational* coefficients, bool flag)
{
  this->Degree = degree;
  this->MaxDegree = degree;
  if (flag == true)
    this->Coefficient = coefficients;
  else
    {
      this->Coefficient = new LongRational [this->Degree + 1];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] = coefficients[i];
    }
  this->RootFlag = false;
  this->RationalRoots = 0;
}

// constructor from raw datas, including all rational roots
//
// degree = polynomial degree
// coefficients = coefficients array ( first element is associated to the -power term)
// roots = array where the roots are stored
// flag = true if coefficient and root arrays have to be used directly and not duplicated

LongRationalPolynomial::LongRationalPolynomial (int degree, LongRational* coefficients, LongRational* roots, bool flag)
{
  this->Degree = degree;
  this->MaxDegree = degree;
  if (flag == true)
    {
      this->Coefficient = coefficients;
      this->RationalRoots = roots;
    }
  else
    {
      this->Coefficient = new LongRational [this->Degree + 1];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] = coefficients[i];
      this->RationalRoots = new LongRational [this->Degree];
      for (int i = 0; i < this->Degree; i++)
	this->RationalRoots[i] = roots[i];
    }
  this->RootFlag = false;
}

// copy constructor
//
// P = polynomial to copy

LongRationalPolynomial::LongRationalPolynomial (const LongRationalPolynomial& P)
{
  if (P.Coefficient != 0)
    {
      this->Degree = P.Degree;
      this->MaxDegree = this->Degree;
      this->Coefficient = new LongRational [this->Degree + 1];
      this->RootFlag = P.RootFlag;
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] = P.Coefficient[i];
      if (P.RootFlag == true)
	{
	  this->NbrRoot = P.NbrRoot;
	  this->Root = new Complex [this->NbrRoot];
	  for (int i = 0; i < this->NbrRoot; i++)
	    this->Root[i] = P.Root[i];
	}
      if (P.RationalRoots == 0)
	{
	  this->RationalRoots = 0;
	}
      else
	{
	  this->RationalRoots = new LongRational [this->Degree];
	  for (int i = 0; i < this->Degree; ++i)
	    this->RationalRoots[i] = P.RationalRoots[i];
	}
    }
  else
    {
      this->Degree = 0;
      this->MaxDegree = 0;
      this->Coefficient = 0;
      this->RootFlag = false; 
      this->RationalRoots = 0;
    }
}

// constructor from P1*P2
//
// P1 = first polynomial
// P2 = second polynomial

LongRationalPolynomial::LongRationalPolynomial (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2)
{
  this->RootFlag = false; 
  this->Degree = P2.Degree + P1.Degree;
  this->RationalRoots = 0;
  if ((P2.Coefficient == 0) || (P1.Coefficient == 0))
    {
      this->Degree = 0;
      this->MaxDegree = 0;	  
      this->Coefficient = 0;
    }
  else
    {
      this->Coefficient = new LongRational [this->Degree + 1];
      for (int i = P1.Degree; i >= 0; i--)
	{
	  this->Coefficient[i + P2.Degree] = P1.Coefficient[i];
	  this->Coefficient[i + P2.Degree] *= P2.Coefficient[P2.Degree];
	}
      if (P2.Degree > 0)
	{
	  for (int i = P2.Degree-1; i >= 0; i--)
	    this->Coefficient[i] = 0l;
	  for (int i = P2.Degree - 1; i >= 0; i--)
	    for (int j = P1.Degree; j >= 0; j--)
	      this->Coefficient[i+j] += P1.Coefficient[j] * P2.Coefficient[i];
	}
      if ((P1.RationalRoots != 0) && (P2.RationalRoots != 0))
	{
	  this->RationalRoots = new LongRational[this->Degree];
	  for (int j = 0; j < P1.Degree; ++j)
	    this->RationalRoots[j] = P1.RationalRoots[j];
	  for (int j = 0; j < P2.Degree; ++j)
	    this->RationalRoots[P1.Degree + j] = P2.RationalRoots[j];	  
	}
    }
  this->MaxDegree = this->Degree;
}

// constructor from P1+q^n P2
//
// P1 = first polynomial
// P2 = second polynomial
// degree = polynomial degree of the monomial with which P2 has to be multiplied 

LongRationalPolynomial::LongRationalPolynomial (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2, int degree)
{
  this->RootFlag = false; 
  this->Degree = P2.Degree + degree;
  this->RationalRoots = 0;
  if (P2.Coefficient == 0)
    {
      if (P1.Coefficient == 0)
	{
	  this->Degree = 0;

	  this->Coefficient = 0;
	}
      else
	{
	  this->Degree = P1.Degree;
	  this->Coefficient = new LongRational [this->Degree + 1];
	  for (int i = 0; i <= this->Degree; ++i)
	    this->Coefficient[i] = P1.Coefficient[i];
	}
    }
  else
    if (P1.Degree > this->Degree)
      {
	this->Degree = P1.Degree;
	this->Coefficient = new LongRational [this->Degree + 1];
	for (int i = 0; i <= this->Degree; ++i)
	  this->Coefficient[i] = P1.Coefficient[i];
	for (int i = 0; i <= P2.Degree; ++i)
	  this->Coefficient[i + degree] += P2.Coefficient[i];      
      }
    else
      {
	this->Coefficient = new LongRational [this->Degree + 1];
	if (P1.Coefficient != 0)
	  for (int i = 0; i < degree; ++i)
	    this->Coefficient[i] = P1.Coefficient[i];
	else
	  for (int i = 0; i < degree; ++i)
	    this->Coefficient[i] = 0l;
	for (int i = 0; i <= P2.Degree; ++i)
	  this->Coefficient[i + degree] = P2.Coefficient[i];      
	if (P1.Coefficient != 0)
	  for (int i = degree; i <= P1.Degree; ++i)
	    this->Coefficient[i] += P1.Coefficient[i];      
      }
  this->MaxDegree = this->Degree;
}

// destructor
//

LongRationalPolynomial::~LongRationalPolynomial()
{
  if ((this->RootFlag == true) && (this->NbrRoot != 0))
    delete[] this->Root;
  if (this->Coefficient != 0)
    delete[] this->Coefficient;
  if (this->RationalRoots != 0)
    delete[] this->RationalRoots;
}

// assignement

LongRationalPolynomial& LongRationalPolynomial::operator = (const LongRationalPolynomial& P)
{
  if ((this->MaxDegree < P.Degree) || (this->MaxDegree == 0)) 
    {
      if (this->Coefficient != 0)
	delete[] this->Coefficient;
      if ((this->RootFlag == true) && (this->NbrRoot != 0))
	delete[] this->Root;
      if (this->RationalRoots != 0)
	delete[] this->RationalRoots;
      this->Degree = P.Degree;
      this->MaxDegree = this->Degree;
      if (P.Coefficient != 0)
	{
	  this->Coefficient = new LongRational [this->Degree + 1];
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
      if (P.RationalRoots == 0)
	{
	  this->RationalRoots = 0;
	}
      else
	{
	  this->RationalRoots = new LongRational [this->Degree];
	  for (int i = 0; i < this->Degree; ++i)
	    this->RationalRoots[i] = P.RationalRoots[i];
	}
    }
  else
    {
      this->Degree = P.Degree;
      if (P.Coefficient != 0)
	{
	  for (int i = 0; i <= this->Degree; i++)
	    {
	      this->Coefficient[i] = P.Coefficient[i];
	    }
	}
      else
	{
	  delete[] this->Coefficient;
	  this->Coefficient = 0;
	}
      if (P.RootFlag == true)
	{
	  this->NbrRoot = P.NbrRoot;
	  this->Root = new Complex [this->NbrRoot];
	  for (int i = 0; i < this->NbrRoot; i++)
	    this->Root[i] = P.Root[i];
	}
      if (P.RationalRoots == 0)
	{
	  if (this->RationalRoots != 0)
	    delete[] this->RationalRoots;
	  this->RationalRoots = 0;	  
	}
      else
	{
	  if (this->RationalRoots == 0)
	    {
	      this->RationalRoots = new LongRational[this->MaxDegree];
	    }
	  for (int i = 0; i < this->Degree; ++i)
	    this->RationalRoots[i] = P.RationalRoots[i];
	}
    }
  this->RootFlag = P.RootFlag;
  return *this;
}

// Return polynomial value at a given point
//
// x = point where to evaluate polynomial
// return value = polynomial value at x

LongRational LongRationalPolynomial::PolynomialEvaluate (const LongRational& x)
{
  if (this->Coefficient != 0)
    {
      LongRational Res = this->Coefficient[this->Degree];
      if (this->RationalRoots == 0)
	{
	  for (int i = this->Degree - 1 ; i >= 0; i--)
	    {
	      Res *= x;
	      Res += this->Coefficient[i];
	    }
	}
      else
	{
	  for (int i = this->Degree - 1 ; i >= 0; i--)
	    Res *= (x - this->RationalRoots[i]);
	}
      return Res;
    }
  return LongRational();
}
  
// Return polynomial value at a given point
//
// x = point where to evaluate polynomial
// return value = polynomial value at x

double LongRationalPolynomial::PolynomialEvaluate (double x)
{
  if (this->Coefficient != 0)
    {
      double Res = this->Coefficient[this->Degree].GetNumericalValue();
      for (int i = this->Degree - 1 ; i >= 0; i--)
	{
	  Res *= x;
	  Res += this->Coefficient[i].GetNumericalValue();
	}
      return Res;
    }
  return 0.0;
}
  
// Return polynomial value at a given point
//
// x = point where to evaluate polynomial
// return value = polynomial value at x

Complex LongRationalPolynomial::PolynomialEvaluate (Complex x)
{
  if (this->Coefficient != 0)
    {
      Complex Res (this->Coefficient[this->Degree].GetNumericalValue());
      for (int i = this->Degree - 1 ; i >= 0; i--)
	{
	  Res *= x;
	  Res += this->Coefficient[i].GetNumericalValue();
	}
      return Res;
    }
  return Complex();
}
  
// Evaluate polynomial derivative 
//
// x = position where to evaluate polynomial derivative
// return value = polynomial derivative at x

LongRational LongRationalPolynomial::DerivativeEvaluate (LongRational x)
{
  if (this->Coefficient != 0)
    {
      LongRational Res = this->Coefficient[this->Degree] * ((LongRational) this->Degree);
      for (long i = this->Degree - 1 ; i > 0; i--)
	{
	  Res *= x;
	  Res += i * this->Coefficient[i];
	}
      return Res;
    }
  return 0l;
}
  
// Evaluate polynomial derivative 
//
// x = position where to evaluate polynomial derivative
// return value = polynomial derivative at x

double LongRationalPolynomial::DerivativeEvaluate (double x)
{
  if (this->Coefficient != 0)
    {
      double Res = this->Coefficient[this->Degree].GetNumericalValue() * ((double) this->Degree);
      for (long i = this->Degree - 1 ; i > 0; i--)
	{
	  Res *= x;
	  Res += ((double) i) * this->Coefficient[i].GetNumericalValue();
	}
      return Res;
    }
  return 0.0;
}
  
// Evaluate polynomial derivative 
//
// x = position where to evaluate polynomial derivative
// return value = polynomial derivative at x

Complex LongRationalPolynomial::DerivativeEvaluate (Complex x)
{
  if (this->Coefficient != 0)
    {
      Complex Res = this->Coefficient[this->Degree].GetNumericalValue() * ((double) this->Degree);
      for (long i = this->Degree - 1 ; i > 0; i--)
	{
	  Res *= x;
	  Res += ((double) i) * this->Coefficient[i].GetNumericalValue();
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

double LongRationalPolynomial::DerivativeEvaluate (double x, int n)
{
  if (this->Coefficient != 0)
    {
      if (n > this->Degree)
	return 0.0;
      double tmp = 1.0;
      for (int j = this->Degree - n + 1; j <= this->Degree; j++)
	tmp *= (double) j;
      double Res = tmp * this->Coefficient[this->Degree].GetNumericalValue();
      for (int i = this->Degree - 1; i >= n; i--)
	{
	  tmp = 1.0;
	  for (int j = i - n + 1; j <= i; j++)
	    tmp *= (double) j;      
	  Res *= x;
	  Res += tmp * this->Coefficient[i].GetNumericalValue();
	}
      return Res;
    }
  return 0.0;
}
  
// Evaluate polynomial n-th derivative at a given value
//
// x = position where to evaluate polynomial n-th derivative
// n = derivative order
// return value = polynomial n-th derivative at x

Complex LongRationalPolynomial::DerivativeEvaluate (Complex x, int n)
{
  if (this->Coefficient != 0)
    {
      if (n > this->Degree)
	return Complex();
      double tmp = 1.0;
      for (int j = this->Degree - n + 1; j <= this->Degree; j++)
	tmp *= (double) j;
      Complex Res (tmp * this->Coefficient[this->Degree].GetNumericalValue());
      for (int i = this->Degree - 1; i >= n; i--)
	{
	  tmp = 1.0;
	  for (int j = i - n + 1; j <= i; j++)
	    tmp *= (double) j;      
	  Res *= x;
	  Res += tmp * this->Coefficient[i].GetNumericalValue();
	}
      return Res;
    }
  return Complex();
}
  
// Return Derivative of the polynomial 
//
// return value = polynomial derivative

LongRationalPolynomial LongRationalPolynomial::DerivatePolynomial ()
{
  if (this->Coefficient != 0)
    {
      LongRational* Coef = new LongRational [this->Degree];
      for (long i = 1; i <= this->Degree; i++)
	Coef[i-1] = (this->Coefficient[i] * i);
      LongRationalPolynomial P (this->Degree - 1, Coef);
      // delete Coef;
      return P;
    }
  return LongRationalPolynomial();
}

// Evaluate polynomial n-th derivative
//
// n = derivative order
// return value = polynomial n-th derivative

LongRationalPolynomial LongRationalPolynomial::DerivatePolynomial (int n)
{
  if (this->Coefficient != 0)
    {
      if (n > this->Degree)
	{
	  LongRational* Coef = new LongRational [1];
	  Coef[0] = 0l;
	  return LongRationalPolynomial(0, Coef, true);
	}
      int tmpDegree = this->Degree - n;
      LongRational* Coef = new LongRational [tmpDegree + 1];  
      LongRational tmp = 1l;
      for (long j = this->Degree - n + 1; j <= this->Degree; ++j)
	tmp *= j;
      Coef[tmpDegree] = tmp * this->Coefficient[this->Degree];
      for (long i = this->Degree - 1; i >= n; i--)
	{
	  tmp = 1l;
	  for (long j = i - n + 1; j <= i; j++)
	    tmp *= j;      
	  Coef[i - n] = tmp * this->Coefficient[i];
	}
      return  LongRationalPolynomial(tmpDegree, Coef, true);
    }
  return LongRationalPolynomial();
}

// arithmetic operators

LongRationalPolynomial operator - (const LongRationalPolynomial& P)
{
  LongRationalPolynomial Pr(P);
  if (Pr.Coefficient != 0)
    {
      for (int i = 0; i <= Pr.Degree; i++)
	Pr.Coefficient[i] *= -1l;  
    }
  return Pr;  
}

LongRationalPolynomial operator + (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2)
{
  if (P1.Degree > P2.Degree)    
    {
      LongRational* Coef = new LongRational [P1.Degree+1];
      for (int i = P1.Degree; i > P2.Degree; i--)
	Coef[i] = P1.Coefficient[i];
      for (int i = P2.Degree; i >= 0; i--)
	Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];	  
      return LongRationalPolynomial(P1.Degree, Coef, true);
    }
  else
    if (P1.Degree < P2.Degree)    
      {
	LongRational* Coef = new LongRational [P2.Degree+1];
	for (int i = P2.Degree; i > P1.Degree; i--)
	  Coef[i] = P2.Coefficient[i];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];	  
	return LongRationalPolynomial(P2.Degree, Coef, true);

      }
    else
      {
	LongRational* Coef = new LongRational [P1.Degree+1];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];
	int i = P1.Degree;
	while ((Coef[i]== 0l) && (i!=0)) i--;
	if (i == P1.Degree)
	  return LongRationalPolynomial(P2.Degree, Coef, true);
	else 
	  {
	    LongRational* Coef2 = new LongRational [i+1];
	    for (int j = 0; j <= i; j++)
	      Coef2[j] = Coef[j];
	    delete[] Coef;
	    return LongRationalPolynomial(i, Coef2, true);
	  }
      }
}

LongRationalPolynomial operator - (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2)
{
  if (P1.Degree > P2.Degree)    
    {
      LongRational* Coef = new LongRational [P1.Degree+1];
      for (int i = P1.Degree; i > P2.Degree; i--)
	Coef[i] = P1.Coefficient[i];
      for (int i = P2.Degree; i >= 0; i--)
	Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];	  
      return LongRationalPolynomial(P1.Degree, Coef, true);
    }
  else
    if (P1.Degree < P2.Degree)    
      {
	LongRational* Coef = new LongRational [P2.Degree+1];
	for (int i = P2.Degree; i > P1.Degree; i--)
	  Coef[i] = - P2.Coefficient[i];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];	  
	return LongRationalPolynomial(P2.Degree, Coef, true);

      }
    else
      {
	LongRational* Coef = new LongRational [P1.Degree+1];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];
	int i = P1.Degree;
	while ((Coef[i]== 0l) && (i != 0)) 
	  i--;
	if (i == P1.Degree)
	  return LongRationalPolynomial(P2.Degree, Coef, true);
	else 
	  {
	    LongRational* Coef2 = new LongRational [i+1];
	    for (int j = 0; j <= i; j++)
	      Coef2[j] = Coef[j];
	    delete[] Coef;
	    return LongRationalPolynomial(i, Coef2, true);
	  }
      }
}

LongRationalPolynomial operator * (const LongRationalPolynomial& P1, const LongRationalPolynomial& P2)
{
  int Deg = P1.Degree+P2.Degree;
  LongRational* Coef = new LongRational [Deg + 1];
  for (int i = P1.Degree; i >= 0; i--)
    Coef[i+P2.Degree] = P1.Coefficient[i] * P2.Coefficient[P2.Degree];
  if (P2.Degree == 0)
    return LongRationalPolynomial(Deg, Coef, true);
  for (int i = P2.Degree-1; i >= 0; i--)
    Coef[i] = 0l;
  for (int i = P2.Degree - 1; i >= 0; i--)
    for (int j = P1.Degree; j >= 0; j--)
      Coef[i+j] += P1.Coefficient[j] * P2.Coefficient[i];
  if (P1.RationalRoots != 0)
    {
      if (P2.RationalRoots != 0)
	{
	  LongRational* TmpRationalRoots = new LongRational[Deg];
	  for (int j = 0; j < P1.Degree; ++j)
	    TmpRationalRoots[j] = P1.RationalRoots[j];
	  for (int j = 0; j < P2.Degree; ++j)
	    TmpRationalRoots[P1.Degree + j] = P2.RationalRoots[j];
	  return LongRationalPolynomial(Deg, Coef, TmpRationalRoots, true);
	}
    }
  return LongRationalPolynomial(Deg, Coef, true);
}

LongRationalPolynomial operator * (const LongRationalPolynomial& P, const LongRational& d)
{
  LongRationalPolynomial Pr(P);
  if (P.Coefficient != 0)
    {
      for (int i = 0; i <= Pr.Degree; i++)
	Pr.Coefficient[i] *= d;
    }
  return Pr;
}

LongRationalPolynomial operator * (const LongRational& d, const LongRationalPolynomial& P)
{
  LongRationalPolynomial Pr (P);
  if (P.Coefficient != 0)
    {
      for (int i = 0; i <= Pr.Degree; i++)
	Pr.Coefficient[i] *= d;
    }
  return Pr;
}

LongRationalPolynomial& LongRationalPolynomial::operator += (const LongRationalPolynomial& P)
{
  if (P.Coefficient == 0)
    return *this;
  if (this->Coefficient == 0)
    {
      this->Coefficient = new LongRational [P.Degree + 1];
      this->Degree = P.Degree;
      this->MaxDegree = this->Degree;
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] = P.Coefficient[i];
      return *this;
   }
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
      if (P.Degree > this->MaxDegree)
	{
	  LongRational* TmpCoef = new LongRational [P.Degree + 1];
	  for (int i = 0; i <= this->Degree; i++)
	    {
	      TmpCoef[i] = this->Coefficient[i];
	      TmpCoef[i] += P.Coefficient[i];
	    }
	  for (int i = this->Degree + 1; i <= P.Degree; i++)
	    TmpCoef[i] = P.Coefficient[i];
	  delete[] this->Coefficient;
	  this->Coefficient = TmpCoef;
	  this->Degree = P.Degree;
	  this->MaxDegree = this->Degree;
	}
      else
	{
	  for (int i = 0; i <= this->Degree; i++)
	    this->Coefficient[i] += P.Coefficient[i];
	  for (int i = this->Degree + 1; i <= P.Degree; i++)
	    this->Coefficient[i] = P.Coefficient[i];
	  this->Degree = P.Degree;
	}
      return *this;      
    }
  else
    {
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] += P.Coefficient[i];
      int i = this->Degree;
      while ((i > 0) && (this->Coefficient[i]== 0l))
	i--;
      this->Degree = i;
      return *this;      
    }
}

LongRationalPolynomial& LongRationalPolynomial::operator -= (const LongRationalPolynomial& P)
{
  if (P.Coefficient == 0)
    return *this;
  if (this->Coefficient == 0)
    {
      this->Coefficient = new LongRational [P.Degree + 1];
      this->Degree = P.Degree;
      this->MaxDegree = this->Degree;
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] = -P.Coefficient[i];
      return *this;
   }
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
      if (P.Degree > this->MaxDegree)
	{
	  LongRational* TmpCoef = new LongRational [P.Degree + 1];
	  for (int i = 0; i <= this->Degree; i++)
	    this->Coefficient[i] -= P.Coefficient[i];
	  for (int i = this->Degree + 1; i <= P.Degree; i++)
	    TmpCoef[i] = -P.Coefficient[i];
	  delete this->Coefficient;
	  this->Coefficient = TmpCoef;
	  this->Degree = P.Degree;
	  this->MaxDegree = this->Degree;
	}
      else
	{
	  for (int i = 0; i <= this->Degree; i++)
	    this->Coefficient[i] += P.Coefficient[i];
	  for (int i = this->Degree + 1; i <= P.Degree; i++)
	    this->Coefficient[i] = -P.Coefficient[i];	      
	  this->Degree = P.Degree;
	}
      return *this;      
    }
  else
    {
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] -= P.Coefficient[i];
      int i = this->Degree;
      while ((i > 0) && (this->Coefficient[i] == 0l))
	i--;
      if (i != this->Degree)
	{
	  LongRational* TmpCoef = new LongRational [i+1];
	  this->Degree = i;
	  this->MaxDegree = this->Degree;
	  for (int j = i; j >= 0; j--)
	    TmpCoef[j] = this->Coefficient [i]; 
	  delete[] this->Coefficient;
	  this->Coefficient = TmpCoef;
	}
      return *this;      
    }
}

LongRationalPolynomial& LongRationalPolynomial::operator *= (const LongRational& d)
{
  if (this->Coefficient != 0)
    {
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] *= d;
    }
  return *this;
}

LongRationalPolynomial& LongRationalPolynomial::operator *= (long d)
{
  if (this->Coefficient != 0)
    {
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] *= d;
    }
  return *this;
}

LongRationalPolynomial& LongRationalPolynomial::operator *= (const LongRationalPolynomial& P)
{
  if (this->Coefficient != 0)
    {
      int Deg = P.Degree + this->Degree;
      if (Deg > this->MaxDegree)
	{
	  LongRational* Coef = new LongRational [Deg + 1];
	  LongRational Tmp = P.Coefficient[P.Degree];
	  for (int i = this->Degree; i >= 0; --i)
	    Coef[i + P.Degree] = this->Coefficient[i] * Tmp;
	  for (int i = 0; i < P.Degree; ++i)
	    Coef[i] = 0l;
	  if (P.Degree == 0)
	    {
	      this->MaxDegree = this->Degree;
	      delete[] this->Coefficient;
	      this->Coefficient = Coef;
	      return *this;
	    }
	  for (int i = 0; i < P.Degree; ++i)
	    {
	      Tmp = P.Coefficient[i];
	      for (int j = 0; j <= this->Degree; ++j)
		Coef[i + j] += this->Coefficient[j] * Tmp;
	    }
	  delete[] this->Coefficient;
	  this->Coefficient = Coef;
	  if (this->RationalRoots != 0)
	    {
	      if (P.RationalRoots != 0)
		{
		  LongRational* TmpRationalRoots = new LongRational[Deg];
		  for (int j = 0; j < this->Degree; ++j)
		    TmpRationalRoots[j] = this->RationalRoots[j];
		  for (int j = 0; j < P.Degree; ++j)
		    TmpRationalRoots[this->Degree + j] = P.RationalRoots[j];
		  delete[] this->RationalRoots;
		  this->RationalRoots = TmpRationalRoots;
		}
	      else
		{
		  delete[] this->RationalRoots;
		  this->RationalRoots = 0;
		}
	    }
	  else
	    {
	      if ((this->Degree == 0) && (P.RationalRoots != 0))
		{
		  this->RationalRoots = new LongRational[Deg];
		  for (int j = 0; j < Deg; ++j)
		    {
		      this->RationalRoots[j] = P.RationalRoots[j];
		    }
		}
	    }
	  this->Degree = Deg;
	  this->MaxDegree = this->Degree;
	  if (this->RootFlag == true)
	    {
	      delete[] this->Root;
	      this->RootFlag = false;
	    }
	}
      else
	{
	  if (P.Degree == 0)
	    {
	      for (int i = this->Degree; i >= 0; --i)
		this->Coefficient[i] *= P.Coefficient[0];
	      return *this;
	    }
	  LongRational* Coef = new LongRational [this->MaxDegree + 1];
	  LongRational Tmp = P.Coefficient[P.Degree];
	  for (int i = this->Degree; i >= 0; --i)
	    Coef[i + P.Degree] = this->Coefficient[i] * Tmp;
	  for (int i = 0; i < P.Degree; ++i)
	    Coef[i] = 0l;
	  for (int i = 0; i < P.Degree; ++i)
	    {
	      Tmp = P.Coefficient[i];
	      for (int j = 0; j <= this->Degree; ++j)
		Coef[i + j] += this->Coefficient[j] * Tmp;
	    }
	  delete[] this->Coefficient;
	  this->Coefficient = Coef;
	  if (this->RationalRoots != 0)
	    {
	      if (P.RationalRoots != 0)
		{
		  for (int j = 0; j < P.Degree; ++j)
		    this->RationalRoots[this->Degree + j] = P.RationalRoots[j];
		}
	      else
		{
		  delete[] this->RationalRoots;
		  this->RationalRoots = 0;
		}
	    }
	  else
	    {
	      if ((this->Degree == 0) && (P.RationalRoots != 0))
		{
		  this->RationalRoots = new LongRational[Deg];
		  for (int j = 0; j < Deg; ++j)
		    {
		      this->RationalRoots[j] = P.RationalRoots[j];
		    }
		}
	    }
	  this->Degree = Deg;
	  if (this->RootFlag == true)
	    {
	      delete[] this->Root;
	      this->RootFlag = false;
	    }
	}
    }
  return *this;
}

// LongRationalPolynomial& LongRationalPolynomial::MultiplyAndAdd (const LongRationalPolynomial& P1 , const LongRationalPolynomial& P2)
// {
//   if (this->Coefficient != 0)
//     {
//       int Deg = P1.Degree + this->Degree;
//       if (P2.Degree > Deg)
// 	Deg = P2.Degree;
//       if (Deg > this->MaxDegree)
// 	{
// 	  LongRational* Coef = new LongRational [Deg + 1];
// 	  LongRational Tmp = P1.Coefficient[P1.Degree];
// 	  for (int i = 0; i <= Deg; ++i)
// 	    Coef[i] = 0l;
// 	  for (int i = this->Degree; i >= 0; --i)
// 	    Coef[i + P1.Degree] = this->Coefficient[i] * Tmp;
// 	  if (P1.Degree > 0)
// 	    {
// 	      for (int i = 0; i < P1.Degree; ++i)
// 		{
// 		  Tmp = P1.Coefficient[i];
// 		  for (int j = 0; j <= this->Degree; ++j)
// 		    Coef[i + j] += this->Coefficient[j] * Tmp;
// 		}
// 	    }
// 	  for (int i = 0; i < P2.Degree; ++i)
// 	    Coef[i] += P2.Coefficient[i];
// 	  delete[] this->Coefficient;
// 	  this->Coefficient = Coef;
// 	  if (this->RationalRoots != 0)
// 	    {
// 	      delete[] this->RationalRoots;
// 	      this->RationalRoots = 0;
// 	    }
// 	  this->MaxDegree = Deg;
// 	}
//       else
// 	{
// 	  if (P1.Degree == 0)
// 	    {
// 	      for (int i = this->Degree; i >= 0; --i)
// 		this->Coefficient[i] *= P1.Coefficient[0];
// 	      for (int i = P2.Degree; i > this->Degree; --i)
// 		this->Coefficient[i] = P2.Coefficient[i];
// 	      for (int i = 0; i < P2.Degree; ++i)
// 		this->Coefficient[i] += P2.Coefficient[i];
// 	    }
// 	  else
// 	    {
// 	      LongRational* Coef = new LongRational [this->MaxDegree + 1];
// 	      LongRational Tmp = P1.Coefficient[P1.Degree];
// 	      for (int i = 0; i <= Deg; ++i)
// 		Coef[i] = 0l;
// 	      for (int i = this->Degree; i >= 0; --i)
// 		Coef[i + P1.Degree] = this->Coefficient[i] * Tmp;
// 	      for (int i = 0; i < P1.Degree; ++i)
// 		{
// 		  Tmp = P1.Coefficient[i];
// 		  for (int j = 0; j <= this->Degree; ++j)
// 		    Coef[i + j] += this->Coefficient[j] * Tmp;
// 		}
// 	      for (int i = 0; i < P2.Degree; ++i)
// 		Coef[i] += P2.Coefficient[i];
// 	      delete[] this->Coefficient;
// 	      this->Coefficient = Coef;
// 	      if (this->RationalRoots != 0)
// 		{
// 		  delete[] this->RationalRoots;
// 		  this->RationalRoots = 0;
// 		}	      
// 	    }
// 	}
//       this->Degree = Deg;
//       if (this->RootFlag == true)
// 	{
// 	  delete[] this->Root;
// 	  this->RootFlag = false;
// 	}
//     }
//   return *this;
// }

LongRationalPolynomial& LongRationalPolynomial::operator /= (const LongRational& d)
{
  if (this->Coefficient != 0)
    {
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] /= d;
    }
  return *this;
}

LongRationalPolynomial& LongRationalPolynomial::operator /= (long d)
{
  if (this->Coefficient != 0)
    {
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] /= d;
    }
  return *this;
}

// shift all powers from a given value
//
// shift = shift to apply
// return value = reference on the current polynomial

LongRationalPolynomial& LongRationalPolynomial::ShiftPowers(int shift)
{
  if ((shift == 0) || (this->Coefficient == 0))
    return *this;
  int TmpDegree = this->Degree + shift;
  if (TmpDegree > this->MaxDegree)
    {
      LongRational* TmpCoefficients = new LongRational [TmpDegree + 1];
      int i = 0;
      for (; i < shift; ++i)
	TmpCoefficients[i] = 0l;
      for (; i <= TmpDegree; ++i)
	TmpCoefficients[i] = this->Coefficient[i - shift];      
      delete[] this->Coefficient;
      this->Coefficient = TmpCoefficients;
      this->Degree = TmpDegree;
      this->MaxDegree = this->Degree;
    }
  else
    {
      int i = TmpDegree;
      for (; i >= shift; --i)
	this->Coefficient[i] = this->Coefficient[i - shift];      
      for (; i >=0; --i)
	this->Coefficient[i] = 0l;
      this->Degree = TmpDegree;
    }
  if (this->RationalRoots != 0)
    {
      delete[] this->RationalRoots;
      this->RationalRoots = 0;
    }
  return *this;
}

// Divide polynomial by a monomial (z - z0) using Horner scheme
//
// z0 = monomial root
// return value = result of polynomial division

LongRationalPolynomial LongRationalPolynomial::MonomialDivision (const LongRational& z0)
{
  if (this->Coefficient != 0)
    {
      LongRational* Coef1 = new LongRational [this->Degree];
      Coef1[this->Degree - 1] = this->Coefficient[this->Degree];
      for (int i = this->Degree - 1; i > 0; i--)
	{
	  Coef1[i - 1] = this->Coefficient[i];
	  Coef1[i - 1] += (Coef1[i] * z0);
	}
      if (this->RationalRoots == 0)
	{
	  return LongRationalPolynomial(this->Degree - 1, Coef1, true);  
	}
      else
	{
	  LongRational* TmpRoots = new LongRational [this->Degree - 1];
	  int i = 0;
	  for (; (i < this->Degree) && (this->RationalRoots[i] != z0); ++i)
	    {
	      TmpRoots[i] = this->RationalRoots[i];
	    }
	  ++i;
	  for (; i < this->Degree; ++i)
	    {
	      TmpRoots[i - 1] = this->RationalRoots[i];
	    }
	  return LongRationalPolynomial(this->Degree - 1, Coef1, TmpRoots, true);  
	}
   }
  return LongRationalPolynomial();
}

// Divide the current polynomial by a monomial (z - z0) using Horner scheme
//
// z0 = monomial root
// return value = reference on the current polynomial

LongRationalPolynomial& LongRationalPolynomial::LocalMonomialDivision (const LongRational& z0)
{
  if (this->Coefficient != 0)
    {
      LongRational Tmp = this->Coefficient[this->Degree - 1];
      this->Coefficient[this->Degree - 1] = this->Coefficient[this->Degree];
      for (int i = this->Degree - 1; i > 0; i--)
	{
	  Tmp += (this->Coefficient[i] * z0);
	  Swap(Tmp, this->Coefficient[i - 1]);
	}
      if (this->RationalRoots != 0)
	{
	  int i = 0;
	  while ((i < this->Degree) && (this->RationalRoots[i] != z0))
	    ++i;
	  ++i;
	  for (; i < this->Degree; ++i)
	    {
	      Swap(this->RationalRoots[i - 1], this->RationalRoots[i]);
	    }
	}
      --this->Degree;
   }
  return *this;
}

// Multiply polynomial by a monomial (z - z0)
//
// z0 = monomial root
// return value = result of polynomial multiplication

LongRationalPolynomial LongRationalPolynomial::MonomialMultiplication (const LongRational& z0)
{
  if (this->Coefficient != 0)
    {
      LongRational* Coef1 = new LongRational [this->Degree + 2];
      Coef1[this->Degree + 1] = this->Coefficient[this->Degree];
      for (int i = this->Degree - 1; i >= 0; i--)
	Coef1[i + 1] = this->Coefficient[i] - (this->Coefficient[i + 1] * z0);
      Coef1[0] = - (z0 * this->Coefficient[0]);
      if ((this->RationalRoots == 0) && (this->Degree > 0))
	{
	  return LongRationalPolynomial(this->Degree + 1, Coef1, true);  
	}
      else
	{
	  LongRational* TmpRoots = new LongRational [this->Degree + 1];
	  for (int i = 0; i < this->Degree; ++i)
	    {
	      TmpRoots[i] = this->RationalRoots[i];
	    }
	  TmpRoots[this->Degree] = z0;
	  return LongRationalPolynomial(this->Degree + 1, Coef1, TmpRoots, true);  
	}
    }
  return LongRationalPolynomial();
}

// Multiply the current polynomial by a monomial (z - z0)
//
// z0 = monomial root
// return value = reference on the current polynomial

LongRationalPolynomial& LongRationalPolynomial::LocalMonomialMultiplication (const LongRational& z0)
{
  if (this->Coefficient != 0)
    {
      if (this->MaxDegree == this->Degree)
	{
	  LongRational* Coef1 = new LongRational [this->Degree + 2];
	  Coef1[this->Degree + 1] = this->Coefficient[this->Degree];
	  for (int i = this->Degree - 1; i >= 0; i--)
	    Coef1[i + 1] = this->Coefficient[i] - (this->Coefficient[i + 1] * z0);
	  Coef1[0] = - (z0 * this->Coefficient[0]);
	  delete[] this->Coefficient;
	  this->Coefficient = Coef1;
	  if ((this->RationalRoots == 0) && (this->Degree > 0))
	    {
	      ++this->Degree;
	      ++this->MaxDegree;
	    }
	  else
	    {
	      LongRational* TmpRoots = new LongRational [this->Degree + 1];
	      for (int i = 0; i < this->Degree; ++i)
		{
		  TmpRoots[i] = this->RationalRoots[i];
		}
	      TmpRoots[this->Degree] = z0;
	      delete[] this->RationalRoots;
	      this->RationalRoots = TmpRoots;
	      ++this->Degree;
	      ++this->MaxDegree;
	    }
	  return *this;
	}
      else
	{
	  this->Coefficient[this->Degree + 1] = this->Coefficient[this->Degree];
	  LongRational Tmp = -z0;
	  for (int i = this->Degree - 1; i >= 0; i--)
	    {
	      this->Coefficient[i + 1] *= Tmp;
	      this->Coefficient[i + 1] += this->Coefficient[i];
	    }	  
	  this->Coefficient[0] *= Tmp;
	  if ((this->RationalRoots == 0) && (this->Degree > 0))
	    {
	      ++this->Degree;
	    }
	  else
	    {
	      if (this->Degree == 0)
		{
		  this->RationalRoots = new LongRational [this->MaxDegree];
		  this->RationalRoots[this->Degree] = z0;	  
		}
	      else
		{
		  this->RationalRoots[this->Degree] = z0;
		}
	      ++this->Degree;	      
	    }
	}
    }
  return *this;
}

// Output Stream overload
//

ostream& operator << (ostream& Str, const LongRationalPolynomial& P)
{
  if (P.Coefficient != 0)
    {
      if (P.RationalRoots != 0)
	{
	  if (P.Coefficient[P.Degree] != 1l)
	    Str << P.Coefficient[P.Degree];
	  for (int i = 0; i < P.Degree; i++)
	    {
	      Str << "(q - " << P.RationalRoots[i] << ")";
	    }
	  return Str;  
	}
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
	    if (P.Coefficient[i] == 1l)
	      Str <<"q^" << i << " + ";
	    else
	      Str << P.Coefficient[i] << " q^" << i << " + ";
	  Str << P.Coefficient[0];
	  return Str;  
	}
    }
  return Str;  
}

// Find Roots of the polynomial
//

void LongRationalPolynomial::SolvePolynomial()
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

void LongRationalPolynomial::SolveLinear()
{
  this->RootFlag = true;
  this->NbrRoot = 1;
  this->Root = new Complex [1];
  Root[0] = - this->Coefficient[0].GetNumericalValue() / this->Coefficient[1].GetNumericalValue();
}      


// Sort roots of a polynomial
//

void LongRationalPolynomial::SortRoots ()
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

