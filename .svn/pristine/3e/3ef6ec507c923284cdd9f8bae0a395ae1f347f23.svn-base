////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   Class of polynomial with integer coefficients            //
//                                                                            //
//                        last modification : 27/05/2005                      //
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


#include "Polynomial/IntegerPolynomial.h"
#include <math.h>


using std::cout;
using std::endl;


// default constructor
//

IntegerPolynomial::IntegerPolynomial ()
{
  this->Degree = 0;
  this->Coefficient = 0;
  this->RootFlag = false; 
}
 
// constructor
//
// degree = polynomial degree

IntegerPolynomial::IntegerPolynomial (int degree)
{
  this->Degree = degree;
  this->Coefficient = new long [this->Degree + 1];
  for (int i = 0; i <= this->Degree; ++i)
    this->Coefficient[i] = 0l;
  this->RootFlag = false; 
  
}

// constructor from raw datas
//
// degree = polynomial degree
// coefficients = coefficients array ( first element is associated to the -power term)
// flag = true if coefficients array has to be used directly and not duplicated

IntegerPolynomial::IntegerPolynomial (int degree, long* coefficients, bool flag)
{
  this->Degree = degree;
  if (flag == true)
    this->Coefficient = coefficients;
  else
    {
      this->Coefficient = new long [this->Degree + 1];
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] = coefficients[i];
    }
  this->RootFlag = false;
}

// copy constructor
//
// P = polynomial to copy

IntegerPolynomial::IntegerPolynomial (const IntegerPolynomial& P)
{
  if (P.Coefficient != 0)
    {
      this->Degree = P.Degree;
      this->Coefficient = new long [this->Degree + 1];
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
    }
  else
    {
      this->Degree = 0;
      this->Coefficient = 0;
      this->RootFlag = false; 
    }
}

// constructor from P1+q^n P2
//
// P1 = first polynomial
// P2 = second polynomial
// degree = polynomial degree of the monomial with which P2 has to be multiplied 

IntegerPolynomial::IntegerPolynomial (const IntegerPolynomial& P1, const IntegerPolynomial& P2, int degree)
{
  this->RootFlag = false; 
  this->Degree = P2.Degree + degree;
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
	  this->Coefficient = new long [this->Degree + 1];
	  for (int i = 0; i <= this->Degree; ++i)
	    this->Coefficient[i] = P1.Coefficient[i];
	}
    }
  else
    if (P1.Degree > this->Degree)
      {
	this->Degree = P1.Degree;
	this->Coefficient = new long [this->Degree + 1];
	for (int i = 0; i <= this->Degree; ++i)
	  this->Coefficient[i] = P1.Coefficient[i];
	for (int i = 0; i <= P2.Degree; ++i)
	  this->Coefficient[i + degree] += P2.Coefficient[i];      
      }
    else
      {
	this->Coefficient = new long [this->Degree + 1];
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
}

// destructor
//

IntegerPolynomial::~IntegerPolynomial()
{
  if ((this->RootFlag == true) && (this->NbrRoot != 0))
    delete[] this->Root;
  if (this->Coefficient != 0)
    delete[] this->Coefficient;
}

// assignement

IntegerPolynomial& IntegerPolynomial::operator = (const IntegerPolynomial& P)
{
  if (this->Coefficient != 0)
    delete[] this->Coefficient;
  if ((this->RootFlag == true) && (this->NbrRoot != 0))
    delete[] this->Root;
  this->Degree = P.Degree;
  this->RootFlag = P.RootFlag;
  if (P.Coefficient != 0)
    {
      this->Coefficient = new long [this->Degree + 1];
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

// Return polynomial value at a given point
//
// x = point where to evaluate polynomial
// return value = polynomial value at x

long IntegerPolynomial::PolynomialEvaluate (long x)
{
  if (this->Coefficient != 0)
    {
      long Res = this->Coefficient[this->Degree];
      for (int i = this->Degree - 1 ; i >= 0; i--)
	{
	  Res *= x;
	  Res += this->Coefficient[i];
	}
      return Res;
    }
  return 0l;
}
  
// Return polynomial value at a given point
//
// x = point where to evaluate polynomial
// return value = polynomial value at x

double IntegerPolynomial::PolynomialEvaluate (double x)
{
  if (this->Coefficient != 0)
    {
      double Res = (double) this->Coefficient[this->Degree];
      for (int i = this->Degree - 1 ; i >= 0; i--)
	{
	  Res *= x;
	  Res += (double) this->Coefficient[i];
	}
      return Res;
    }
  return 0.0;
}
  
// Return polynomial value at a given point
//
// x = point where to evaluate polynomial
// return value = polynomial value at x

Complex IntegerPolynomial::PolynomialEvaluate (Complex x)
{
  if (this->Coefficient != 0)
    {
      Complex Res ((double) this->Coefficient[this->Degree]);
      for (int i = this->Degree - 1 ; i >= 0; i--)
	{
	  Res *= x;
	  Res += (double) this->Coefficient[i];
	}
      return Res;
    }
  return Complex();
}
  
// Evaluate polynomial derivative 
//
// x = position where to evaluate polynomial derivative
// return value = polynomial derivative at x

long IntegerPolynomial::DerivativeEvaluate (long x)
{
  if (this->Coefficient != 0)
    {
      long Res = this->Coefficient[this->Degree] * ((long) this->Degree);
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

double IntegerPolynomial::DerivativeEvaluate (double x)
{
  if (this->Coefficient != 0)
    {
      double Res = (double) (this->Coefficient[this->Degree] * ((long) this->Degree));
      for (long i = this->Degree - 1 ; i > 0; i--)
	{
	  Res *= x;
	  Res += (double) (i * this->Coefficient[i]);
	}
      return Res;
    }
  return 0.0;
}
  
// Evaluate polynomial derivative 
//
// x = position where to evaluate polynomial derivative
// return value = polynomial derivative at x

Complex IntegerPolynomial::DerivativeEvaluate (Complex x)
{
  if (this->Coefficient != 0)
    {
      Complex Res = (double) (this->Coefficient[this->Degree] * ((long) this->Degree));
      for (long i = this->Degree - 1 ; i > 0; i--)
	{
	  Res *= x;
	  Res += (double) (i * this->Coefficient[i]);
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

double IntegerPolynomial::DerivativeEvaluate (double x, int n)
{
  if (this->Coefficient != 0)
    {
      if (n > this->Degree)
	return 0.0;
      double tmp = 1.0;
      for (int j = this->Degree - n + 1; j <= this->Degree; j++)
	tmp *= (double) j;
      double Res = tmp * (double) this->Coefficient[this->Degree];
      for (int i = this->Degree - 1; i >= n; i--)
	{
	  tmp = 1.0;
	  for (int j = i - n + 1; j <= i; j++)
	    tmp *= (double) j;      
	  Res *= x;
	  Res += tmp * (double) this->Coefficient[i];
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

Complex IntegerPolynomial::DerivativeEvaluate (Complex x, int n)
{
  if (this->Coefficient != 0)
    {
      if (n > this->Degree)
	return Complex();
      double tmp = 1.0;
      for (int j = this->Degree - n + 1; j <= this->Degree; j++)
	tmp *= (double) j;
      Complex Res (tmp * (double) this->Coefficient[this->Degree]);
      for (int i = this->Degree - 1; i >= n; i--)
	{
	  tmp = 1.0;
	  for (int j = i - n + 1; j <= i; j++)
	    tmp *= (double) j;      
	  Res *= x;
	  Res += tmp * (double) this->Coefficient[i];
	}
      return Res;
    }
  return Complex();
}
  
// Return Derivative of the polynomial 
//
// return value = polynomial derivative

IntegerPolynomial IntegerPolynomial::DerivatePolynomial ()
{
  if (this->Coefficient != 0)
    {
      long* Coef = new long [this->Degree];
      for (long i = 1; i <= this->Degree; i++)
	Coef[i-1] = (this->Coefficient[i] * i);
      IntegerPolynomial P (this->Degree - 1, Coef);
      delete Coef;
      return P;
    }
  return IntegerPolynomial();
}

// Evaluate polynomial n-th derivative
//
// n = derivative order
// return value = polynomial n-th derivative

IntegerPolynomial IntegerPolynomial::DerivatePolynomial (int n)
{
  if (this->Coefficient != 0)
    {
      if (n > this->Degree)
	{
	  long* Coef = new long [1];
	  Coef[0] = 0l;
	  return IntegerPolynomial(0, Coef, true);
	}
      int tmpDegree = this->Degree - n;
      long* Coef = new long [tmpDegree + 1];  
      long tmp = 1l;
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
      return  IntegerPolynomial(tmpDegree, Coef, true);
    }
  return IntegerPolynomial();
}

// arithmetic operators

IntegerPolynomial operator - (const IntegerPolynomial& P)
{
  IntegerPolynomial Pr(P);
  if (Pr.Coefficient != 0)
    {
      for (int i = 0; i <= Pr.Degree; i++)
	Pr.Coefficient[i] *= -1l;  
    }
  return Pr;  
}

IntegerPolynomial operator + (const IntegerPolynomial& P1, const IntegerPolynomial& P2)
{
  if (P1.Degree > P2.Degree)    
    {
      long* Coef = new long [P1.Degree+1];
      for (int i = P1.Degree; i > P2.Degree; i--)
	Coef[i] = P1.Coefficient[i];
      for (int i = P2.Degree; i >= 0; i--)
	Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];	  
      return IntegerPolynomial(P1.Degree, Coef, true);
    }
  else
    if (P1.Degree < P2.Degree)    
      {
	long* Coef = new long [P2.Degree+1];
	for (int i = P2.Degree; i > P1.Degree; i--)
	  Coef[i] = P2.Coefficient[i];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];	  
	return IntegerPolynomial(P2.Degree, Coef, true);

      }
    else
      {
	long* Coef = new long [P1.Degree+1];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] + P2.Coefficient[i];
	int i = P1.Degree;
	while ((Coef[i] == 0l) && (i!=0)) i--;
	if (i == P1.Degree)
	  return IntegerPolynomial(P2.Degree, Coef, true);
	else 
	  {
	    long* Coef2 = new long [i+1];
	    for (int j = 0; j <= i; j++)
	      Coef2[j] = Coef[j];
	    delete[] Coef;
	    return IntegerPolynomial(i, Coef2, true);
	  }
      }
}

IntegerPolynomial operator - (const IntegerPolynomial& P1, const IntegerPolynomial& P2)
{
  if (P1.Degree > P2.Degree)    
    {
      long* Coef = new long [P1.Degree+1];
      for (int i = P1.Degree; i > P2.Degree; i--)
	Coef[i] = P1.Coefficient[i];
      for (int i = P2.Degree; i >= 0; i--)
	Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];	  
      return IntegerPolynomial(P1.Degree, Coef, true);
    }
  else
    if (P1.Degree < P2.Degree)    
      {
	long* Coef = new long [P2.Degree+1];
	for (int i = P2.Degree; i > P1.Degree; i--)
	  Coef[i] = - P2.Coefficient[i];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];	  
	return IntegerPolynomial(P2.Degree, Coef, true);

      }
    else
      {
	long* Coef = new long [P1.Degree+1];
	for (int i = P1.Degree; i >= 0; i--)
	  Coef[i] = P1.Coefficient[i] - P2.Coefficient[i];
	int i = P1.Degree;
	while ((Coef[i] == 0l) && (i != 0)) 
	  i--;
	if (i == P1.Degree)
	  return IntegerPolynomial(P2.Degree, Coef, true);
	else 
	  {
	    long* Coef2 = new long [i+1];
	    for (int j = 0; j <= i; j++)
	      Coef2[j] = Coef[j];
	    delete[] Coef;
	    return IntegerPolynomial(i, Coef2, true);
	  }
      }
}

IntegerPolynomial operator * (const IntegerPolynomial& P1, const IntegerPolynomial& P2)
{
  int Deg = P1.Degree+P2.Degree;
  long* Coef = new long [Deg + 1];
  for (int i = P1.Degree; i >= 0; i--)
    Coef[i+P2.Degree] = P1.Coefficient[i] * P2.Coefficient[P2.Degree];
  if (P2.Degree == 0)
    return IntegerPolynomial(Deg, Coef, true);
  for (int i = P2.Degree-1; i >= 0; i--)
    Coef[i] = 0l;
  for (int i = P2.Degree - 1; i >= 0; i--)
    for (int j = P1.Degree; j >= 0; j--)
      Coef[i+j] += P1.Coefficient[j] * P2.Coefficient[i];
  return IntegerPolynomial(Deg, Coef, true);
}

IntegerPolynomial operator * (const IntegerPolynomial& P, const long& d)
{
  IntegerPolynomial Pr(P);
  if (P.Coefficient != 0)
    {
      for (int i = 0; i <= Pr.Degree; i++)
	Pr.Coefficient[i] *= d;
    }
  return Pr;
}

IntegerPolynomial operator * (const long& d, const IntegerPolynomial& P)
{
  IntegerPolynomial Pr (P);
  if (P.Coefficient != 0)
    {
      for (int i = 0; i <= Pr.Degree; i++)
	Pr.Coefficient[i] *= d;
    }
  return Pr;
}

IntegerPolynomial& IntegerPolynomial::operator += (const IntegerPolynomial& P)
{
  if (P.Coefficient == 0)
    return *this;
  if (this->Coefficient == 0)
    {
      this->Coefficient = new long [P.Degree + 1];
      this->Degree = P.Degree;
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
      long* TmpCoef = new long [P.Degree + 1];
      for (int i = 0; i <= this->Degree; i++)
	TmpCoef[i] = this->Coefficient[i] + P.Coefficient[i];
      for (int i = this->Degree + 1; i <= P.Degree; i++)
	TmpCoef[i] = P.Coefficient[i];
      delete this->Coefficient;
      this->Coefficient = TmpCoef;
      this->Degree = P.Degree;
      return *this;      
    }
  else
    {
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] += P.Coefficient[i];
      int i = this->Degree;
      while ((i > 0) && (this->Coefficient[i] == 0l))
	i--;
      if (i != this->Degree)
	{
	  long* TmpCoef = new long [i+1];
	  this->Degree = i;
	  for (int j = i; j >= 0; j--)
	    TmpCoef[j] = this->Coefficient [i]; 
	  delete[] this->Coefficient;
	  this->Coefficient = TmpCoef;
	}
      return *this;      
    }
}

IntegerPolynomial& IntegerPolynomial::operator -= (const IntegerPolynomial& P)
{
  if (P.Coefficient == 0)
    return *this;
  if (this->Coefficient == 0)
    {
      this->Coefficient = new long [P.Degree + 1];
      this->Degree = P.Degree;
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
      long* TmpCoef = new long [P.Degree + 1];
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
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] -= P.Coefficient[i];
      int i = this->Degree;
      while ((i > 0) && (this->Coefficient[i] == 0l))
	i--;
      if (i != this->Degree)
	{
	  long* TmpCoef = new long [i+1];
	  this->Degree = i;
	  for (int j = i; j >= 0; j--)
	    TmpCoef[j] = this->Coefficient [i]; 
	  delete[] this->Coefficient;
	  this->Coefficient = TmpCoef;
	}
      return *this;      
    }
}

IntegerPolynomial& IntegerPolynomial::operator *= (const long& d)
{
  if (this->Coefficient != 0)
    {
      for (int i = 0; i <= this->Degree; i++)
	this->Coefficient[i] *= d;
    }
  return *this;
}

IntegerPolynomial& IntegerPolynomial::operator *= (const IntegerPolynomial& P)
{
  if (this->Coefficient != 0)
    {
      int Deg = P.Degree + this->Degree;
      long* Coef = new long [Deg + 1];
      long Tmp = P.Coefficient[P.Degree];
      for (int i = this->Degree; i >= 0; --i)
	Coef[i + P.Degree] = this->Coefficient[i] * Tmp;
      for (int i = 0; i < P.Degree; ++i)
	Coef[i] = 0;
      if (P.Degree == 0)
	{
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
      this->Degree = Deg;
      if (this->RootFlag == true)
	{
	  delete[] this->Root;
	  this->RootFlag = false;
	}
    }
  return *this;
}

// shift all powers from a given value
//
// shift = shift to apply
// return value = reference on the current polynomial

IntegerPolynomial& IntegerPolynomial::ShiftPowers(int shift)
{
  if ((shift == 0) || (this->Coefficient == 0))
    return *this;
  int TmpDegree = this->Degree + shift;
  long* TmpCoefficients = new long [TmpDegree + 1];
  int i = 0;
  for (; i < shift; ++i)
    TmpCoefficients[i] = 0l;
  for (; i <= TmpDegree; ++i)
    TmpCoefficients[i] = this->Coefficient[i - shift];      
  delete[] this->Coefficient;
  this->Coefficient = TmpCoefficients;
  this->Degree = TmpDegree;
  return *this;
}

// Output Stream overload
//

ostream& operator << (ostream& Str, const IntegerPolynomial& P)
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

void IntegerPolynomial::SolvePolynomial()
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

void IntegerPolynomial::SolveLinear()
{
  this->RootFlag = true;
  this->NbrRoot = 1;
  this->Root = new Complex [1];
  Root[0] = - ((double) this->Coefficient[0]) / ((double) this->Coefficient[1]);
}      


// Sort roots of a polynomial
//

void IntegerPolynomial::SortRoots ()
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

