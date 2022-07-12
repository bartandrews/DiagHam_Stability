#include "ArithmeticGeometricMean.h"

#include <cmath>
#include <cstdlib>

using std::sqrt;

// default constructor
//
ArithmeticGeometricMean::ArithmeticGeometricMean()
{
  this->Length=0;
}

// constructor 
//
// a0, b0, c0 = initial triple of numbers
// precision = required precision
ArithmeticGeometricMean::ArithmeticGeometricMean(double a0, double b0, double c0, double precision)
{
  this->Precision = precision;
  this->Flag.Initialize();

  double an=a0, bn=b0, cn=c0;

  this->Length=2;
  this->Recurse(an,bn,cn);
  while (fabs(cn) > this->Precision)
    {
      ++this->Length;
      this->Recurse(an,bn,cn);
    }

  this->A=new double[this->Length];
  this->B=new double[this->Length];
  this->C=new double[this->Length];  

  this->A[0]=a0;
  this->B[0]=b0;
  this->C[0]=c0;

  an=a0, bn=b0, cn=c0;

  this->Recurse(an,bn,cn);
  
  this->A[1]=an;
  this->B[1]=bn;
  this->C[1]=cn;
  
  this->Length=2;

  while (fabs(cn)>this->Precision)
    {
      this->Recurse(an,bn,cn);
      this->A[this->Length]=an;
      this->B[this->Length]=bn;
      this->C[this->Length]=cn;
      ++this->Length;      
    }  
}

// copy constructor (without duplicating datas)
//
// coefficients = reference on Clebsch Gordan coefficients to copy
ArithmeticGeometricMean::ArithmeticGeometricMean (const ArithmeticGeometricMean& agm)
{
  this->A=agm.A;
  this->B=agm.B;
  this->C=agm.C;
  this->Length=agm.Length;
  this->Flag=agm.Flag;
  this->Precision=agm.Precision;
}

// destructor
//
ArithmeticGeometricMean::~ArithmeticGeometricMean ()
{
  if ((this->Length!=0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] A;
      delete [] B;
      delete [] C;
    }
}

// do one step in the recursion
// a = on entry a_n on return: a_n+1
// b = on entry b_n ... 
// c = on entry c_n ...
void ArithmeticGeometricMean::Recurse(double &a, double &b, double &c)
{
  double an, bn;
  an=(a+b)/2.0;
  bn=sqrt(a*b);  
  c = (a-b)/2.0;
  a=an;
  b=bn;
}
