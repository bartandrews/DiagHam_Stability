#include "LegendrePolynomials.h"
#include <iostream>
#include <cmath>
#include <cassert>

// precision below which abscissae are treated as equal on consecutive calls
double LegendrePolynomials::Precision = 1e-13;
   
/// Constructor for LegendrePolynomials
/// Class generating legendre polynomials on the standard interval [-1,1]
/// @param MaxOrder maximum order to be calculated
LegendrePolynomials::LegendrePolynomials(int MaxOrder) :
  NbrLegendrePolynomials (MaxOrder+1),
  Values((int)std::max((int) 3,MaxOrder+1), true) ,
  LastArgument(0.0),
  LastOrder(0)
{
}


/// evaluate a specific Legendre Polynomial at x
/// @param n requested order
/// @param x argument in [-1,1]
double LegendrePolynomials::GetPolynomial(int n, double x)
{
  if (IsStored(n,x))
    return this->Values[n];
  else
    {
      this->UpdatePolynomials(x,n);
      return this->Values[n];
    }
}


/// Update values for Legendre polynomials and return reference to internal storage
/// @param x argument in [-1,1]
/// @param[out] output array of output elements
/// @param requiredOrder order up to which it should be evaluated
RealVector const & LegendrePolynomials::GetPolynomials(double x)
{
  this->Values[0]=1.0;
  this->Values[1]=x;
  for (int n=2; n<NbrLegendrePolynomials; ++n)
      this->Values[n]=((2.0*n-1.0)*x*this->Values[n-1] - (n-1.0)*this->Values[n-2])/n;
  this->LastOrder=NbrLegendrePolynomials-1;
  this->LastArgument=x;
  return this->Values;
}



/// Calculate values for Legendre polynomials up to the required order and return to a supplied array
/// @param x argument in [-1,1]
/// @param[out] output array of output elements
/// @param requiredOrder order up to which it should be evaluated
void LegendrePolynomials::GetPolynomials(double x, double *output, int requiredOrder)
{
  if (requiredOrder<0) requiredOrder = this->NbrLegendrePolynomials-1;
  output[0]=1.0;
  output[1]=x;
  for (int n=2; n<=requiredOrder; ++n)
      output[n]=((2.0*n-1.0)*x*output[n-1] - (n-1.0)*output[n-2])/n;
}



/// Update values for Legendre polynomials up to the required order
/// @param x argument in [-1,1]
/// @param requiredOrder order up to which it should be evaluated
void LegendrePolynomials::UpdatePolynomials(double x, int requiredOrder)
{
  assert(requiredOrder < NbrLegendrePolynomials);
  this->Values[0]=1.0;
  this->Values[1]=x;
  for (int n=2; n<=requiredOrder; ++n)
    this->Values[n]=((2.0*n-1.0)*x*this->Values[n-1] - (n-1.0)*this->Values[n-2])/n;
  this->LastOrder=requiredOrder;
  this->LastArgument=x;
}

/// Update values for Legendre polynomials up to the given order, if recalculation required
/// @param x argument in [-1,1]
/// @param requiredOrder order up to which it should be evaluated
void LegendrePolynomials::ConditionalUpdatePolynomials(double x, int requiredOrder)
{
  assert(requiredOrder < NbrLegendrePolynomials);
  if (x<-1.0 || x> 1.0)
    std::cerr << "Argument out of bounds in LegendrePolynomials"<<std::endl;
  if (fabs(x-this->LastArgument) < Precision)
    {
      // only calculate required orders
      for (int n=LastOrder+1; n<=requiredOrder; ++n)
        this->Values[n]=((2.0*n-1.0)*x*this->Values[n-1] - (n-1.0)*this->Values[n-2])/n;
      if (requiredOrder > this->LastOrder) this->LastOrder=requiredOrder;
    }
  else
    {
      this->Values[0]=1.0;
      this->Values[1]=x;
      for (int n=2; n<=requiredOrder; ++n)
        this->Values[n]=((2.0*n-1.0)*x*this->Values[n-1] - (n-1.0)*this->Values[n-2])/n;
      this->LastOrder=requiredOrder;
      this->LastArgument=x;
    }
}


