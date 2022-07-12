#include "JacobiPolynomials.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

using std::cout;
using std::endl;
using std::floor;


// default constructor
//
JacobiPolynomials::JacobiPolynomials()
{
  this->ParameterA = 0.0;
  this->ParameterB = 0.0;
  this->MaxDegreeN = 1;
  this->FunctionValues = new double[MaxDegreeN+1];
  this->FunctionValues[0]=1.0;
  this->CNPlusOne = NULL;
  this->CNConst = NULL;
  this->CNLinear = NULL;
  this->CNMinusOne = NULL;
  this->A1Zero = 0.5*(this->ParameterA-this->ParameterB);
  this->A1One = 0.5*(this->ParameterA+this->ParameterB)+1.0;
  this->ExplicitExpansionPrefactors=new double*[MaxDegreeN+1];
  for (int i=0; i<=MaxDegreeN; ++i)
    this->ExplicitExpansionPrefactors [i] = NULL; 
  this->LastN = -1;

}

// constructor for a general theta function \theta[^a_b](z|tau)
// maxDegreeN = maximum degree of the function
// param_a = parameter a
// param_b = parameter b
// tau = modulus
// precision = required precision
JacobiPolynomials::JacobiPolynomials(int maxDegreeN, double a, double b)
{
  if (maxDegreeN<0)
    {
      cout << "JacobiPolynomials not implemented for negative degree n"<<endl;
      exit(1);
    }
  this->ParameterA=a;
  this->ParameterB=b;
  this->MaxDegreeN=maxDegreeN;
  this->FunctionValues = new double[MaxDegreeN+1];
  
  this->CNPlusOne = new double[MaxDegreeN-1];
  this->CNConst = new double[MaxDegreeN-1];
  this->CNLinear = new double[MaxDegreeN-1];
  this->CNMinusOne = new double[MaxDegreeN-1];

  this->FunctionValues[0]=1.0;
  this->A1Zero = 0.5*(this->ParameterA-this->ParameterB);
  this->A1One = 0.5*(this->ParameterA+this->ParameterB)+1.0;
  this->ExplicitExpansionPrefactors=new double*[MaxDegreeN+1];
  for (int i=0; i<=MaxDegreeN; ++i)
    this->ExplicitExpansionPrefactors [i] = NULL; 

  this->InitializeRecursion();

  this->LastN = -1;  
}

// copy constructor (with duplicating of datas)
//
// 
JacobiPolynomials::JacobiPolynomials (const JacobiPolynomials& p)
{
  this->ParameterA=p.ParameterA;
  this->ParameterB=p.ParameterB;
  this->MaxDegreeN=p.MaxDegreeN;
  this->FunctionValues = new double[MaxDegreeN+1];
  for (int i=0; i<=MaxDegreeN; ++i)
    this->FunctionValues[i] = p.FunctionValues[i];
  if (MaxDegreeN>1)
    {
      this->CNPlusOne = new double[MaxDegreeN-1];
      this->CNConst = new double[MaxDegreeN-1];
      this->CNLinear = new double[MaxDegreeN-1];
      this->CNMinusOne = new double[MaxDegreeN-1];
      for (int i=0; i<MaxDegreeN-1; ++i)
	{
	  this->CNPlusOne[i] = p.CNPlusOne[i];
	  this->CNConst[i] = p.CNConst[i];
	  this->CNLinear[i] = p.CNLinear[i];
	  this->CNMinusOne[i] = p.CNMinusOne[i];
	}
    }
  this->ExplicitExpansionPrefactors=new double*[MaxDegreeN+1];
  for (int i=0; i<=MaxDegreeN; ++i)
    if (p.ExplicitExpansionPrefactors [i] != NULL)
      {
	this->ExplicitExpansionPrefactors[i]=new double[i+1];
	for (int j=0; j<=i; ++j)
	  this->ExplicitExpansionPrefactors[i][j]=p.ExplicitExpansionPrefactors[i][j];
      }
    else
      {
	this->ExplicitExpansionPrefactors[i]=NULL;
      }

  this->A1Zero = 0.5*(this->ParameterA-this->ParameterB);
  this->A1One = 0.5*(this->ParameterA+this->ParameterB)+1.0;
  this->LastArgument = p.LastArgument;
  this->LastN = p.LastN;

}

// destructor
//
JacobiPolynomials::~JacobiPolynomials ()
{
  if (MaxDegreeN>1)
    {
      delete [] this->CNPlusOne;
      delete [] this->CNConst;
      delete [] this->CNLinear;
      delete [] this->CNMinusOne;
    }
  for (int i=0; i<=MaxDegreeN; ++i)
    if (this->ExplicitExpansionPrefactors [i] != NULL)
      delete [] this->ExplicitExpansionPrefactors[i];
  delete [] this->ExplicitExpansionPrefactors;
  delete [] this->FunctionValues;
}


// get the value of the function for a given coordinate z
//
// n = degree (must be <=MaxDegreeN)
// x = argument
// return = function value of P_n(x)
double JacobiPolynomials::GetValue(int n, double x)
{
  if ((x!=LastArgument)||(n>LastN))
    {
      LastArgument=x;
      LastN=n;
      this->RunRecursion(n,x);
    }
  return FunctionValues[n];
}

// get the value of the function for a given coordinate z
//
// x = argument
// return = function value of P_n(x)
double* JacobiPolynomials::GetValues(double x)
{
  if ((x!=LastArgument)||(LastN<MaxDegreeN))
    {
      LastArgument=x;
      LastN=MaxDegreeN;
      this->RunRecursion(MaxDegreeN,x);
    }
  return this->FunctionValues;
}



  // pretty-print a function value
  // str = stream to print to
  // z = point where to evaluate
ostream& JacobiPolynomials::PrintValues(ostream &str, double x)
{
  this->GetValues(x);
  for (int i=0; i<=MaxDegreeN; ++i)
    str << "P_"<<i<<"^("<<this->ParameterA<<", "<<this->ParameterB<<") ("<<x<<")="<<this->FunctionValues[i]<<endl;
  return str;
}


// evaluate recursively up to degree n
void JacobiPolynomials::RunRecursion(int &limitN, double x)
{
  this->FunctionValues[0]=1.0;
  this->FunctionValues[1]=A1Zero+A1One*x;
  if (limitN>this->MaxDegreeN)
    {
      cout << "Attention: JacobiPolynomials object was initialized to calculate polynomials up to degree "<<this->MaxDegreeN<<endl;      
      limitN=this->MaxDegreeN;
    }
  for (int n=2; n<=limitN; ++n)
    {
      if (fabs(CNPlusOne[n-2])>1e-12)
	{
	  double A = (CNConst[n-2] + CNLinear[n-2]*x)*this->FunctionValues[n-1];
	  double B = CNMinusOne[n-2]*this->FunctionValues[n-2];
	  double C=A-B;
	  if (fabs(fabs(C)-fabs(A))<1e-8)
	    this->FunctionValues[n] = this->GetExplicitFunctionValue(n,x);
	  else
	    this->FunctionValues[n] = C/CNPlusOne[n-2];
	}
      else
	{
	  this->FunctionValues[n] = this->GetExplicitFunctionValue(n,x);
	}
    }
}

// initialize recursion coefficients
void JacobiPolynomials::InitializeRecursion()
{
  double SumAB = this->ParameterA+this->ParameterB;
  for (int n=1; n<this->MaxDegreeN; ++n)
    {
      CNPlusOne[n-1] = 2.0*(n+1.0)*(n+SumAB+1.0)*(2.0*n+SumAB);
      CNConst[n-1] = (2.0*n+SumAB+1.0)*(this->ParameterA*this->ParameterA-this->ParameterB*this->ParameterB);
      CNLinear[n-1] = (2.0*n+SumAB)*(2.0*n+SumAB+1.0)*(2.0*n+SumAB+2.0);
      CNMinusOne[n-1] = 2.0*(n+this->ParameterA)*(n+this->ParameterB)*(2.0*n+SumAB+2.0);
    }
}

// get the value of the function from a direct series expansion rather than recursively
// n = degree (must be <=MaxDegreeN)
// x = argument
// return = function value of P_n(x)
double JacobiPolynomials::GetExplicitFunctionValue(int n, double x)
{
  if (ExplicitExpansionPrefactors[n]==NULL)
    {
      if ((fabs(this->ParameterA-nearbyint(this->ParameterA))>1e-12)||(fabs(this->ParameterB-nearbyint(this->ParameterB))>1e-12))
	{
	  cout << "Cannot perform explicit evaluation except for integer values of alpha and beta in JacobiPolynomial"<<endl;
	  exit(1);
	}
      int A=(int)nearbyint(this->ParameterA);
      int B=(int)nearbyint(this->ParameterB);
      FactorialCoefficient Prefactor;
      this->ExplicitExpansionPrefactors[n]=new double[n+1];
      for (int m=0; m<-B; ++m)
	this->ExplicitExpansionPrefactors[n][m]=0.0;
      for (int m=n+A+1; m<=n; ++m)
	this->ExplicitExpansionPrefactors[n][m]=0.0;
      int minM=(B>=0?0:-B);
      int maxM=(A>=0?n:n+A);
      for (int m=minM; m<=maxM; ++m)
	{
	  Prefactor.SetToOne();

	  Prefactor.FactorialMultiply(n+A);
	  Prefactor.FactorialDivide(n+A-m);
	  Prefactor.FactorialDivide(m);

	  Prefactor.FactorialMultiply(n+B);
	  Prefactor.FactorialDivide(n-m);
	  Prefactor.FactorialDivide(B+m);

	  Prefactor.Power2Divide(n);

	  this->ExplicitExpansionPrefactors[n][m]=Prefactor.GetNumericalValue();
	}
    }
  double Result=0.0;
  for (int m=0; m<=n; ++m)
    {
      double XMinus1=x-1.0;
      double Term=1.0;
      for (int k=0; k<n-m; ++k)
	Term*=XMinus1;
      double XPlus1=x+1.0;
      for (int k=0; k<m; ++k)
	Term*=XPlus1;
      Result+=Term*this->ExplicitExpansionPrefactors[n][m];
    }
  return Result;
}
