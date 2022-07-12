#include "MonopoleHarmonics.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/JacobiPolynomials.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

using std::cout;
using std::endl;
using std::floor;


// default constructor
//
MonopoleHarmonics::MonopoleHarmonics()
{
  this->ChargeTwoQ=0;
  this->MomentumTwoM=0;
  this->MaxOrderN = 0;
  this->FunctionValues = new double[1];
  this->NormPrefactor = new double[1];
  this->InitializeNorms();
  this->FunctionValues[0]=1.0/sqrt(4.0*M_PI);
  this->JacobiOffset=0;
  this->Jacobi = new JacobiPolynomials(0, (-ChargeTwoQ+MomentumTwoM)/2,(ChargeTwoQ+MomentumTwoM)/2 );
  this->LastN = -1;

}

// constructor for a sequence of monopole harmonics Y_q,n,m with varying n (q,m fixed)
// maxOrderN = maximum order of the function
// twoQ = twice the monopole charge
// twoM = twice the Lz component of angular momentum
MonopoleHarmonics::MonopoleHarmonics(int twoQ, int maxOrderN, int twoM)
{
  this->ChargeTwoQ=twoQ;
  this->MomentumTwoM=twoM;
  
  this->MaxOrderN=maxOrderN;
  this->FunctionValues = new double[MaxOrderN+1];
  this->NormPrefactor = new double[MaxOrderN+1];
  
  this->InitializeNorms();
  int maxDegreeP = MaxOrderN + ((ChargeTwoQ<0?-ChargeTwoQ:ChargeTwoQ)-MomentumTwoM)/2;
  if (maxDegreeP<0)
    maxDegreeP = 0;
  this->JacobiOffset = ((ChargeTwoQ<0?-ChargeTwoQ:ChargeTwoQ)-MomentumTwoM)/2;
  this->Jacobi = new JacobiPolynomials(maxDegreeP, (-ChargeTwoQ+MomentumTwoM)/2,(ChargeTwoQ+MomentumTwoM)/2 );
}

// copy constructor (with duplicating of datas)
//
// 
MonopoleHarmonics::MonopoleHarmonics (const MonopoleHarmonics& p)
{
  
  this->ChargeTwoQ=p.ChargeTwoQ;
  this->MomentumTwoM=p.MomentumTwoM;
  
  this->MaxOrderN=p.MaxOrderN;
  this->FunctionValues = new double[MaxOrderN+1];
  this->NormPrefactor = new double[MaxOrderN+1];
  for (int i=0; i<=MaxOrderN; ++i)
    this->NormPrefactor[i]=p.NormPrefactor[i];
  int maxDegreeP = MaxOrderN + ((ChargeTwoQ<0?-ChargeTwoQ:ChargeTwoQ)-MomentumTwoM)/2;
  this->Jacobi = new JacobiPolynomials(maxDegreeP, (-ChargeTwoQ+MomentumTwoM)/2,(ChargeTwoQ+MomentumTwoM)/2 );
  this->JacobiOffset = p.JacobiOffset;
  this->LastArgument = p.LastArgument;
  this->LastN = p.LastN;

}

// destructor
//
MonopoleHarmonics::~MonopoleHarmonics ()
{
  delete [] this->FunctionValues;
  delete [] this->NormPrefactor;
  delete this->Jacobi;
}


// get the value of the function for a given coordinate z
//
// n = degree (must be <=MaxOrderN)
// x = argument
// return = function value of P_n(x)
double MonopoleHarmonics::GetValue(int n, double x)
{
  if ((n>0)&&((x==1.0)||(x==-1.0)))
    return 0.0;
  double XFactor = pow(1.0-x,(MomentumTwoM-ChargeTwoQ)/4.0)*pow(1.0+x,(MomentumTwoM+ChargeTwoQ)/4.0);
  double JacobiValue = Jacobi->GetValue(JacobiOffset+n,x);
  return NormPrefactor[n] * JacobiValue * XFactor;
}

// get function value
//
// x = argument
// return = function value of P_n(x)
double * MonopoleHarmonics::GetValues(double x)
{

  if ((x==1.0)||(x==-1.0))
    {
      FunctionValues[0] = NormPrefactor[0];
      for (int i=1; i<=MaxOrderN; ++i)
	FunctionValues[i]=0.0;
    }
  else
    {
      double XFactor = pow(1.0-x,(MomentumTwoM-ChargeTwoQ)/4.0)*pow(1.0+x,(MomentumTwoM+ChargeTwoQ)/4.0);
      double *JacobiValues = Jacobi->GetValues(x);
      int minN=(JacobiOffset>0?0:-JacobiOffset);
      for (int i=0; i<minN; ++i)
	FunctionValues[i] = 0.0;
      for (int i=minN; i<=MaxOrderN; ++i)
	{
	  FunctionValues[i] = NormPrefactor[i] * JacobiValues[JacobiOffset+i] * XFactor;
	  //cout << "NormPrefactor["<<i<<"]="<<NormPrefactor[i]<<" ,JacobiValues["<<JacobiOffset+i<<", "<< (-ChargeTwoQ+MomentumTwoM)/2
	  //     << ", " << (ChargeTwoQ+MomentumTwoM)/2 <<"; "<<x<<"]="<<JacobiValues[JacobiOffset+i]<<" Xfx="<<XFactor<<endl;
	}
    }
  return FunctionValues;
}


  // pretty-print a function value
  // str = stream to print to
  // z = point where to evaluate
ostream& MonopoleHarmonics::PrintValues(ostream &str, double x)
  {
  this->GetValues(x);
  for (int i=0; i<=MaxOrderN; ++i)
    str << "Y^"<<this->ChargeTwoQ<<"/2_("<<i<<", "<<MomentumTwoM<<"/2) ["<<x<<"]="<<this->FunctionValues[i]<<endl;
  return str;
}


// initialize recursion coefficients
void MonopoleHarmonics::InitializeNorms()
{
  int AbsTwoQ=(ChargeTwoQ<0?-ChargeTwoQ:ChargeTwoQ);
  FactorialCoefficient Factor;
  for (int n=0; n<=MaxOrderN; ++n)
    {
      Factor.SetToOne();
      Factor.FactorialMultiply((AbsTwoQ-MomentumTwoM)/2+n);
      Factor.FactorialDivide(n);
      Factor.FactorialMultiply((AbsTwoQ+MomentumTwoM)/2+n);
      Factor.FactorialDivide(AbsTwoQ+n);
      if (MomentumTwoM>0)
	{
	  Factor.Power2Divide(MomentumTwoM/2);
	  if (MomentumTwoM&1)
	    NormPrefactor[n]=1.0/sqrt(2.0);
	}
      else
	{
	  Factor.Power2Multiply(-MomentumTwoM/2);
	  if (MomentumTwoM&1)
	    NormPrefactor[n]=sqrt(2.0);
	}
      NormPrefactor[n]*=sqrt(Factor.GetNumericalValue())*sqrt((AbsTwoQ+2.0*n+1.0)/(4.0*M_PI));
    }
}


