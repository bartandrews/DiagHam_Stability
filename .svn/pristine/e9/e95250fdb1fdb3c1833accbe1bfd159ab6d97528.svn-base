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

#include "SpecialPolynomial.h"
#include "MathTools/FactorialCoefficient.h"
#include <cassert>
using std::ostream;
using std::cout;
using std::endl;

// return a Laguerre Polynomial of rank n (optionally associated Laguerre function)
// n = index
// alpha = modifier
Polynomial LaguerrePolynomial(int n, int alpha, bool verbose)
{
  double *Coefficients;
  if (n<0)
    {
      Coefficients=new double[1];
      Coefficients[0]=0.0;
      n=0;
    }
  else
    {
      Coefficients=new double[n+1];
      FactorialCoefficient MyCoefficient; 
      for (int m=0; m<=n; ++m)
	{
	  MyCoefficient.SetToOne();
	  MyCoefficient.FactorialMultiply(n+alpha);
	  MyCoefficient.FactorialDivide(m);
	  MyCoefficient.FactorialDivide(alpha+m);
	  MyCoefficient.FactorialDivide(n-m);
	  Coefficients[m]=MyCoefficient.GetNumericalValue();
	  if (m&1)
	    Coefficients[m]*=-1.0;
	  if (verbose)
	    cout << "Coefficient m="<<m<<" for Laguerre_(n="<<n<<", alpha="<<alpha<<") is "<<Coefficients[m]<<endl;
	}
    }
  return Polynomial (n, Coefficients, true);
}



// constructor
// maximum order to be calculated
LaguerrePolynomialRecursion::LaguerrePolynomialRecursion(int nmax)
{
  TmpValues.Resize(nmax+1);
  NMax = nmax;
}
  
// return the value for a Laguerre Polynomial of rank n
// n = index
double LaguerrePolynomialRecursion::EvaluateLaguerrePolynomial(int N, double x)
{
  assert(N<=this->NMax);
  TmpValues[0]=1.0;
  TmpValues[1]=1.0-x;
  for (int n=2; n<=N; ++n)
    TmpValues[n] = ( (2.*n-1.-x)*TmpValues[n-1] - (n-1)*TmpValues[n-2]) / n;
  return TmpValues[N];							     
}

// evaluate all Laguerre polynomials up to the order of the length of rst at the given coordinate
// rst = return value for polynomials
// x = coordinate where to evaluate
void LaguerrePolynomialRecursion::EvaluateLaguerrePolynomials(RealVector& rst, double x)
{
  rst[0]=1.0;
  rst[1]=1.0-x;
  for (int n=2; n<rst.GetVectorDimension(); ++n)
    rst[n] = ( (2.*n-1.-x)*rst[n-1] - (n-1)*rst[n-2]) / n;
}


#ifdef __GMP__
void LaguerrePolynomialRecursion::CheckAgainstGMP(double xDouble, mp_bitcnt_t precision)
{
  LaguerrePolynomialRecursionAP RecursionAP(this->NMax, precision);
  mpf_t x;
  mpf_init2 (x, precision);

  mpf_t rst;
  mpf_init2 (rst, precision);

  mpf_set_d (x, xDouble);
  double rstDouble, rstAP;

  mpf_t TmpValue;
  mpf_t TmpValue2;
  mpf_init2 (TmpValue, precision);
  mpf_init2 (TmpValue2, precision);


  mpf_t *rsts = new mpf_t[this->NMax+1];
  for (int i=0;  i<=this->NMax; ++i)
     mpf_init2 (rsts[i], precision);

  RecursionAP.EvaluateLaguerrePolynomials(this->NMax, rsts, x);
    
  for (int n=0; n<=NMax; ++n)
    {
      RecursionAP.EvaluateLaguerrePolynomial(n, rst, x);
      rstAP=mpf_get_d(rst);
      rstDouble=this->EvaluateLaguerrePolynomial(n, xDouble);
      if (fabs(rstAP - rstDouble)>1e-10)
	cout << "For x="<<xDouble<<", the discrepancy between Double and AP at n="<<n<<" is "<< rstAP - rstDouble<<", AP="<<rstAP<<", Double="<< rstDouble<<endl;

      mpf_reldiff(TmpValue, rst, rsts[n]);
      mpf_set_d(TmpValue2,0.0);
      if (mpf_cmp (TmpValue, TmpValue2))
	cout << "For x="<<xDouble<<", have a relative discrepancy between AP1 and AP2 codes of "<<mpf_get_d(TmpValue)<<endl;
    }

    
}
#endif




#ifdef __GMP__

// constructor
// nmax = maximum order to be calculated
// precision = number of bits precision
LaguerrePolynomialRecursionAP::LaguerrePolynomialRecursionAP(int nmax, mp_bitcnt_t precision)
{
  mpf_set_default_prec(precision);
  TmpValues = new mpf_t [nmax+1];
  for (int n=0; n<nmax; ++n)
    mpf_init (TmpValues[n]);
  mpf_init2 (TmpValue, precision);
  mpf_init2 (TmpValue2, precision);
  this->NMax = nmax;
}


// destructor
LaguerrePolynomialRecursionAP::~LaguerrePolynomialRecursionAP()
{
  delete [] TmpValues;
}


// return the value for a Laguerre Polynomial of rank n
// n = index
void LaguerrePolynomialRecursionAP::EvaluateLaguerrePolynomial(int N, mpf_t &result, const mpf_t &x)
{
  assert(N<=this->NMax);
  mpf_set_d (TmpValues[0], 1.0);

  mpf_set_d (TmpValue,1.0);
  mpf_sub (TmpValues[1], TmpValue, x);
  for (int n=2; n<=N; ++n)
    {
      mpf_set_d (TmpValue,2.*n-1.);
      mpf_sub (TmpValue2, TmpValue, x);
      mpf_mul(TmpValues[n], TmpValue2, TmpValues[n-1]);
      mpf_set_d (TmpValue,n-1.);
      mpf_mul(TmpValue2, TmpValue, TmpValues[n-2]);
      mpf_sub (TmpValue, TmpValues[n], TmpValue2);
      mpf_set_d (TmpValue2,n);
      mpf_div (TmpValues[n], TmpValue, TmpValue2);
    }
  return mpf_set(result, TmpValues[N]); 
}

// evaluate all Laguerre polynomials up to the order of the length of rst at the given coordinate
// rst = return value for polynomials
// x = coordinate where to evaluate
void LaguerrePolynomialRecursionAP::EvaluateLaguerrePolynomials(int N, mpf_t* rst, mpf_t &x)
{
  mpf_set_d (rst[0], 1.0);

  mpf_set_d (TmpValue,1.0);
  mpf_sub (rst[1], TmpValue, x);
  for (int n=2; n<=N; ++n)
    {
      mpf_set_d (TmpValue,2.*n-1.);
      mpf_sub (TmpValue2, TmpValue, x);
      mpf_mul(rst[n], TmpValue2, rst[n-1]);
      mpf_set_d (TmpValue,n-1.);
      mpf_mul(TmpValue2, TmpValue, rst[n-2]);
      mpf_sub (TmpValue, rst[n], TmpValue2);
      mpf_set_d (TmpValue2,n);
      mpf_div (rst[n], TmpValue, TmpValue2);
    }
}

#endif
