#include "JacobiThetaDerivativeFunction.h"

#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::floor;


// maximum number of iterations before convergence has to be reached
#define ITERMAX 500


// default constructor
//
JacobiThetaDerivativeFunction::JacobiThetaDerivativeFunction() : JacobiThetaFunction()
{
}

// constructor for a general theta function \theta[^a_b](z|tau)
//
// param_a = parameter a
// param_b = parameter b
// tau = modulus
// precision = required precision
JacobiThetaDerivativeFunction::JacobiThetaDerivativeFunction(double a, double b, Complex tau) : JacobiThetaFunction(a,b,tau)
{
}

// constructor for one of the four Jacobi theta functions \theta_i(z|tau)
//
// type = index i of theta-function
// tau = modulus
// precision = required precision
JacobiThetaDerivativeFunction::JacobiThetaDerivativeFunction(int i, Complex tau ) : JacobiThetaFunction(i,tau)
{
}
      

// copy constructor (without duplicating datas)
//
// coefficients = reference on Clebsch Gordan coefficients to copy
JacobiThetaDerivativeFunction::JacobiThetaDerivativeFunction (const JacobiThetaDerivativeFunction& theta)
{
  this->ParameterA=theta.ParameterA;
  this->ParameterB=theta.ParameterB;
  this->Tau=theta.Tau;
  this->SumOffset=theta.SumOffset;
  this->ShiftPhase=theta.ShiftPhase;
}

// destructor
//
JacobiThetaDerivativeFunction::~JacobiThetaDerivativeFunction ()
{
}


// get the value of the function for a given coordinate z
//
// z = complex coordinate
// return = function value at z

Complex JacobiThetaDerivativeFunction::GetValue(const Complex &z)
{
  Complex IPiTau = Complex(0.0,M_PI)*this->Tau;
  Complex TwoPiI = Complex(0.0,2.0*M_PI);
  Complex TwoPiIZ = TwoPiI*(z+this->ParameterB);
  Complex Convergence=10.0;
  Complex Result=0.0;
  double NMinus, NPlus;
  cout <<this->Tau<<endl;
  if (SumOffset^1) // ParameterA < 0.5
    {
      Result = ParameterA * TwoPiI * exp(IPiTau*(ParameterA*ParameterA)+TwoPiIZ*ParameterA);
            cout<<"Adding n=0: "<<Result<<", argument of exp: "<<IPiTau*(ParameterA*ParameterA)+TwoPiIZ*ParameterA<<", 1: "<<IPiTau*(ParameterA*ParameterA)<<", 2:"<< TwoPiIZ*ParameterA <<", A:"<< ParameterA <<", (z+b)"<<z+ParameterB<<", 2pi(z+b)="<<TwoPiIZ<<endl;
    }

  int Iter=0;
  while ((SqrNorm(Convergence)>1e-29)&&(Iter++<ITERMAX))
    {
      NMinus = (-Iter+SumOffset)+this->ParameterA;
      NPlus = Iter+this->ParameterA;
      Convergence = NMinus *TwoPiI * exp(IPiTau*(NMinus*NMinus)+TwoPiIZ*NMinus); // term with negative n
            cout<<"Adding n="<<-Iter+SumOffset<<": "<<Convergence<<endl;
      Convergence+= NPlus * TwoPiI * exp(IPiTau*(NPlus*NPlus)+TwoPiIZ*NPlus); // term with positive n
            cout<<"Adding n="<<Iter<<": "<<NPlus * TwoPiI*exp(IPiTau*(NPlus*NPlus)+TwoPiIZ*NPlus)<<endl;
      Result+=Convergence;
    }
  return Result*this->ShiftPhase;
}

// get the value of the function for a given coordinate z
// values = complex vector where to store results
// manyZ = complex vector of coordinates
// return = function value at z
void JacobiThetaDerivativeFunction::GetManyValues(ComplexVector &values, ComplexVector &manyZ)
{
  int Dim = manyZ.GetVectorDimension();  
  Complex IPiTau = Complex(0.0,M_PI)*this->Tau;
  Complex TwoPiI = Complex(0.0,2.0*M_PI);
  Complex Term1M, Term1P;
  ComplexVector TwoPiIZ(Dim);
  Complex Convergence;
  Complex TotalConvergence=10.0;
  double NMinus, NPlus;
  
  if (SumOffset^1) // ParameterA < 0.5
    {
      values.Resize(Dim);
      for (int i=0; i<Dim; ++i)
	{
	  TwoPiIZ[i]= Complex(0.0,2.0*M_PI)*(manyZ[i]+this->ParameterB);
	  values[i]= ParameterA * TwoPiI * exp(IPiTau*(ParameterA*ParameterA)+TwoPiIZ[i]*ParameterA);
	}
    }
  else
    {
      values.ResizeAndClean(Dim);
      for (int i=0; i<Dim; ++i)
	TwoPiIZ[i]= Complex(0.0,2.0*M_PI)*(manyZ[i]+this->ParameterB);
    }

  int Iter=0;
  while ((SqrNorm(TotalConvergence)>1e-29)&&(Iter++<ITERMAX))
    {
      NMinus = (-Iter+SumOffset)+this->ParameterA;
      NPlus = Iter+this->ParameterA;
      Term1M=IPiTau*(NMinus*NMinus);
      Term1P=IPiTau*(NPlus*NPlus);
      for (int i=0; i<Dim; ++i)
	{
	  Convergence = NMinus *TwoPiI * exp(Term1M+TwoPiIZ[i]*NMinus); // term with negative n
	  Convergence+= NPlus *TwoPiI * exp(Term1P+TwoPiIZ[i]*NPlus);  // term with positive n
	  TotalConvergence+=Convergence;
	  values[i]+=Convergence;
	}
      
    }
  if (SqrNorm(this->ShiftPhase-1.0)>1e-24)
    for (int i=0; i<Dim; ++i)
      values[i]*=this->ShiftPhase;
}

  // pretty-print a function value
  // str = stream to print to
  // z = point where to evaluate
ostream& JacobiThetaDerivativeFunction::PrintValue(ostream &str, const Complex &z)
{
  str << "theta["<<this->ParameterA<<", "<<this->ParameterB<<"]("<<z.Re<<"+I*"<<z.Im<<"|"
      <<this->Tau.Re<<"+I*"<<this->Tau.Im<<")="<<this->GetValue(z);
  return str;
}
