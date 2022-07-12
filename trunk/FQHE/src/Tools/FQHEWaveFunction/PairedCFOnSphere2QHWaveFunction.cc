////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//           class implementing generalized Halperin states on the sphere          //
//                                                                            //
//                        last modification : 02/11/2007                      //
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


#include "config.h"
#include "PairedCFOnSphere2QHWaveFunction.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
//

PairedCFOnSphere2QHWaveFunction::PairedCFOnSphere2QHWaveFunction()
{
  this->NbrParticles = 0;
  this->NbrParameters = 0;
}

// constructor
//
// nbrParticles = number of particles
// nbrLandauLevel = number of Landau levels filled with composite fermions
// nbrEffectiveFlux = number of flux quanta of the magnetic monopole field experienced by CF's
// lValue = value of angular momentum for total state
// mValue = value of z component of total angular momentum
// MooreReadCoefficient = prefactor of singular 1/z term in pair-wave function
// CFCoefficients = prefactors of CF orbitals in shells 0, 1, 2, ... , nbrLandauLevel-1
// quasiParticleMCSteps = number of steps for internal integration over QP coordinates
// correctPrefactors = flag that enables the correction of prefactors to adopt the conventions of previous code
// jastrowPower = power to which the Jastrow factor has to be raised
//
PairedCFOnSphere2QHWaveFunction::PairedCFOnSphere2QHWaveFunction(int nbrParticles, int nbrLandauLevels, int nbrEffectiveFlux, int lValue, int mValue, double MooreReadCoefficient, double * givenCFCoefficients,
					 int quasiParticleMCSteps, bool correctPrefactors, int jastrowPower)
{
  this->NbrParticles = nbrParticles;
  if (nbrParticles&1)
    {
      cout << "QH States with odd number of electrons not defined, yet."<<endl;
      exit(1);
    }
  this->NbrLandauLevels = nbrLandauLevels;
  this->NbrParameters = this->NbrLandauLevels+1; // inherited field
  this->LValue = lValue;
  if ((2*lValue > nbrParticles) || ( (lValue&1)!= ((nbrParticles/2)&1)) )
    {
      cout << "Valid angular momenta of two quasiholes are L=[N/2, N/2-2,...]" << endl;
      exit(1);
    }
  this->MValue = mValue;
  if (abs(mValue)>lValue)
    {
      cout << "Valid z-component of the angular momentum has to satisfy -L<=M<=L" << endl;
      exit(1);
    }
  this->QuasiParticleMCSteps = quasiParticleMCSteps;
  this->AbsEffectiveFlux = abs(nbrEffectiveFlux);
  this->Orbitals = new JainCFOnSphereOrbitals(nbrParticles, nbrLandauLevels, nbrEffectiveFlux,jastrowPower);
  this->MooreReadCoefficient=MooreReadCoefficient;
  this->TrialParameters= new double [NbrParameters];
  for (int i=0; i<NbrLandauLevels; ++i) this->TrialParameters[i] = givenCFCoefficients[i];
  this->TrialParameters[this->NbrLandauLevels] = MooreReadCoefficient;
  this->ElementNorm=1.0;
  this->Slater = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->Gij = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->QuasiParticles = new ParticleOnSphereCollection(this->NbrParticles);
  
  QHBasis = new ParticleOnSphereFunctionBasis(this->NbrParticles/2,ParticleOnSphereFunctionBasis::LeftHanded);

  QHLzMax = this->NbrParticles/2;

  // calculate coupling coefficients:
  ClebschGordanCoefficients VectorCoupling(this->NbrParticles/2,this->NbrParticles/2);

  this->M1Values=new int[this->NbrParticles/2+1];
  this->M2Values=new int[this->NbrParticles/2+1];
  this->Couplings= new double[this->NbrParticles/2+1];
  this->NbrCouplings=0;  
  for (int m1=-this->NbrParticles/2; m1<=this->NbrParticles/2; m1+=2)
    {
      int m2 = 2*MValue - m1;
      if ((m2>=-this->NbrParticles/2) && (m2<=this->NbrParticles/2))
	{
	  if (VectorCoupling.GetCoefficient (m1, m2, 2*this->LValue) != 0.0)
	    {
	      M1Values[this->NbrCouplings] = (m1+this->NbrParticles/2)/2;
	      M2Values[this->NbrCouplings] = (m2+this->NbrParticles/2)/2;
	      Couplings[this->NbrCouplings]= VectorCoupling.GetCoefficient (m1, m2, 2*this->LValue);
	      cout << "Adding coupling : " <<Couplings[this->NbrCouplings]<<" * (m1="<<m1/2<<", m2="<<m2/2<<")"<<endl;
	      ++this->NbrCouplings;	      
	    }
	}
    }

  RealVector tmpVector(2);
  Complex tmpC;
  tmpVector[0]=QuasiParticles->Theta(0);
  tmpVector[1]=QuasiParticles->Phi(0);
  for (int i=0; i<NbrCouplings; ++i)
    QHBasis->GetFunctionValue(tmpVector, tmpC, M1Values[i]);
  tmpVector[0]=QuasiParticles->Theta(1);
  tmpVector[1]=QuasiParticles->Phi(1);
  for (int i=0; i<NbrCouplings; ++i)
    QHBasis->GetFunctionValue(tmpVector, tmpC, M2Values[i]);

  
  this->Flag.Initialize();
  this->Ji = new Complex[this->NbrParticles];
  this->Jqh1 = new Complex[this->NbrParticles];
  this->Jqh2 = new Complex[this->NbrParticles];
  this->QH1BasisValues = new Complex[this->NbrParticles/2+1];
  this->QH2BasisValues = new Complex[this->NbrParticles/2+1];
  this->gAlpha = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha[i]= new Complex[NbrParticles*NbrParticles];
  
  if (correctPrefactors)
    {
      int p=Orbitals->GetJastrowPower();
      FactorialCoefficient Coef;
      for (int n=0; n<NbrLandauLevels; ++n)
	{
	  Coef.SetToOne();
	  Coef.PartialFactorialDivide(nbrEffectiveFlux+p*(nbrParticles-1)+2,nbrEffectiveFlux+2*p*(nbrParticles-1)+1);
	  Coef.PartialFactorialMultiply(p*(nbrParticles-1)+n+2,2*p*(nbrParticles-1)+n+1);
	  this->TrialParameters[n]*=Coef.GetNumericalValue()*Coef.GetNumericalValue();
	  //cout << "Correction["<<n<<"]="<<Coef.GetNumericalValue()<<"from: " << nbrEffectiveFlux+p*(nbrParticles-1)+2<<","<<nbrEffectiveFlux+2*p*(nbrParticles-1)+1<<","<< p*(nbrParticles-1)+n+2<<","<<2*p*(nbrParticles-1)+n+1<<endl;
	}
    }  
}

// copy constructor
//
// function = reference on the wave function to copy

PairedCFOnSphere2QHWaveFunction::PairedCFOnSphere2QHWaveFunction(const PairedCFOnSphere2QHWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrLandauLevels = function.NbrLandauLevels;
  this->NbrParameters = function.NbrParameters;
  this->LValue = function.LValue;
  this->MValue = function.MValue;
  this->AbsEffectiveFlux = function.AbsEffectiveFlux;
  this->QuasiParticleMCSteps = function.QuasiParticleMCSteps;
  this->QHLzMax = function.QHLzMax;
  this->Flag = function.Flag;
  this->Orbitals = function.Orbitals;
  this->MooreReadCoefficient=function.MooreReadCoefficient;
  this->TrialParameters=function.TrialParameters;
  this->ElementNorm=function.ElementNorm;
  this->QHBasis = new ParticleOnSphereFunctionBasis(this->NbrParticles/2,ParticleOnSphereFunctionBasis::LeftHanded);
  this->M1Values = function.M1Values;
  this->M2Values = function.M2Values;
  this->Couplings = function.Couplings;
  this->NbrCouplings= function.NbrCouplings;

  RealVector tmpVector(2);
  Complex tmpC;
  tmpVector[0]=QuasiParticles->Theta(0);
  tmpVector[1]=QuasiParticles->Phi(0);
  for (int i=0; i<NbrCouplings; ++i)
    QHBasis->GetFunctionValue(tmpVector, tmpC, M1Values[i]);
  tmpVector[0]=QuasiParticles->Theta(1);
  tmpVector[1]=QuasiParticles->Phi(1);
  for (int i=0; i<NbrCouplings; ++i)
    QHBasis->GetFunctionValue(tmpVector, tmpC, M2Values[i]);
  
  this->Slater = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->Gij = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->QuasiParticles = new ParticleOnSphereCollection(this->NbrParticles);
  this->Ji = new Complex[this->NbrParticles];
  this->Jqh1 = new Complex[this->NbrParticles];
  this->Jqh2 = new Complex[this->NbrParticles];
  this->QH1BasisValues = new Complex[this->NbrParticles/2+1];
  this->QH2BasisValues = new Complex[this->NbrParticles/2+1];
  this->gAlpha = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha[i]= new Complex[NbrParticles*NbrParticles];
  
}

// destructor
//

PairedCFOnSphere2QHWaveFunction::~PairedCFOnSphere2QHWaveFunction()
{
  if ( (this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete Orbitals;
      delete [] TrialParameters;
      delete [] M1Values;
      delete [] M2Values;
      delete [] Couplings;
    }
  for (int i=0; i< this->NbrLandauLevels; ++i) delete [] gAlpha[i];
  delete [] gAlpha;
  delete [] Ji;
  delete [] Jqh1;
  delete [] Jqh2;
  delete [] this->QH1BasisValues;
  delete [] this->QH2BasisValues;
  delete [] Couplings;
  delete Slater;
  delete Gij;
  delete QHBasis;
  delete QuasiParticles;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PairedCFOnSphere2QHWaveFunction::Clone ()
{
  return new PairedCFOnSphere2QHWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PairedCFOnSphere2QHWaveFunction::operator ()(RealVector& x)
{
  this->OrbitalValues = (*Orbitals)(x);
  this->EvaluateTables();
  Complex tmp;
	      
  // initialize Slater determinant (or Pfaffian matrix)
  for (int i=0;i<this->NbrParticles;++i)
    {
      for(int j=0;j<i;++j)
	{
	  tmp=0.0;
	  for (int n=0; n<this->NbrLandauLevels; ++n)
	    tmp+=this->TrialParameters[n]*this->gAlpha[n][i*this->NbrParticles+j];	    
	  
	  Gij->SetMatrixElement(i,j, this->ElementNorm*this->Ji[i]*this->Ji[j]
				   *(MooreReadCoefficient/Orbitals->JastrowFactorElement(i,j) + tmp));
	}
    }  
  //cout << *Slater << endl;
  return MonteCarloEvaluate(); 
}


// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex PairedCFOnSphere2QHWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
 {
  this->OrbitalValues = Orbitals->CalculateFromSpinorVariables(uv);
  this->EvaluateTables();
  Complex tmp;
	      
  // initialize Slater determinant (or Pfaffian matrix)
  for (int i=0;i<this->NbrParticles;++i)
    {
      for(int j=0;j<i;++j)
	{
	  tmp=0.0;
	  for (int n=0; n<this->NbrLandauLevels; ++n)
	    tmp+=this->TrialParameters[n]*this->gAlpha[n][i*this->NbrParticles+j];	    
	  
	  Gij->SetMatrixElement(i,j, this->ElementNorm*this->Ji[i]*this->Ji[j]
				   *(MooreReadCoefficient/Orbitals->JastrowFactorElement(i,j) + tmp));
	}
    }  
  //cout << *Slater << endl;
  return this->MonteCarloEvaluate();
} 


// get a value of the wavefunction for the last set of coordinates, but with different variational coefficients
// coefficients: array of variational coefficients f_0, f_1, ...
//                    length has to be identical to the initial set of parameters!
// singular: new prefactor of the Moore-Read part
//
Complex PairedCFOnSphere2QHWaveFunction::GetForOtherParameters( double *coefficients)
{
  Complex tmp;	      
  // initialize Slater determinant (or Pfaffian matrix)
  for (int i=0;i<this->NbrParticles;++i)
    {
      for(int j=0;j<i;++j)
	{
	  tmp=0.0;
	  for (int n=0; n<this->NbrLandauLevels; ++n)
	    tmp+=coefficients[n]*this->gAlpha[n][i*this->NbrParticles+j];	    
	  Gij->SetMatrixElement(i,j, this->ElementNorm*this->Ji[i]*this->Ji[j]
				   *(coefficients[this->NbrLandauLevels]/Orbitals->JastrowFactorElement(i,j) + tmp));
	}
    }  
  return this->MonteCarloEvaluate();
}

// do many evaluations, storing the result in the vector results given in the call
// x: positions to evaluate the wavefuntion in
// format for passing parameters in the matrix coefficients: coefficients[nbrSet][LandauLevel],
// the entry [][NbrLandauLevels] corresponds to the MooreRead Term.
void PairedCFOnSphere2QHWaveFunction::GetForManyParameters(ComplexVector &results, RealVector& x, double **coefficients)
{
  this->OrbitalValues = (*Orbitals)(x);
  this->EvaluateTables();
  Complex tmp;
  int numParamSets=results.GetVectorDimension();
  double *tmpCoefficients;
  for (int s=0; s<numParamSets; ++s)
    {
      tmpCoefficients = coefficients[s];
      // initialize Slater determinant (or Pfaffian matrix)
      for (int i=0;i<this->NbrParticles;++i)
	{
	  for(int j=0;j<i;++j)
	    {
	      tmp=0.0;
	      for (int n=0; n<this->NbrLandauLevels; ++n)
		tmp+=tmpCoefficients[n]*this->gAlpha[n][i*this->NbrParticles+j];	    
	      Gij->SetMatrixElement(i,j, this->ElementNorm*this->Ji[i]*this->Ji[j]
				       *(tmpCoefficients[this->NbrLandauLevels]/Orbitals->JastrowFactorElement(i,j) + tmp));
	    }
	}
      tmp=this->MonteCarloEvaluate();
      results.Re(s)=Real(tmp);
      results.Im(s)=Imag(tmp);
    }
}


void PairedCFOnSphere2QHWaveFunction::SetTrialParameters(double * coefficients)
{
  for (int n=0; n<this->NbrParameters; ++n)
    this->TrialParameters[n]=coefficients[n];
  this->MooreReadCoefficient =  this->TrialParameters[this->NbrLandauLevels];
}

// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void PairedCFOnSphere2QHWaveFunction::AdaptNorm(RealVector& x)
{
  double det=Norm((*this)(x));
  while ((det<.1)||(det>50.0))
    {
      //cout <<"N'="<< this->ElementNorm << " det="<<det<<endl;
      if (det>1e300) 
	this->ElementNorm*= pow((double)1.0e-300,(double)2.0/this->NbrParticles);
      else if (det==0.0) 
	this->ElementNorm*= pow((double)1.0e300,(double)2.0/this->NbrParticles);
      else 
	this->ElementNorm*= pow(det,(double)-2.0/this->NbrParticles);
      det=Norm((*this)(x));
      //cout <<"N'="<< this->ElementNorm << endl;
    }
}


// normalize the wave-function over an average number of MC positions

void PairedCFOnSphere2QHWaveFunction::AdaptAverageMCNorm(int thermalize, int average)
{
  ParticleOnSphereCollection * Particles = new ParticleOnSphereCollection(this->NbrParticles);
  this->AdaptNorm(Particles->GetPositions());
  Complex TmpMetropolis, TrialValue = (*this)(Particles->GetPositions());  
  double PreviousSamplingAmplitude = SqrNorm(TrialValue);
  double CurrentSamplingAmplitude = PreviousSamplingAmplitude;
  int NextCoordinates=0;
  // do some MC moves: accept or reject move according to probability |Psi_new|^2  / |Psi_old|^2
  for (int i = 0; i < thermalize; ++i)
    {
      Particles->Move(NextCoordinates);
      TmpMetropolis = (*this)(Particles->GetPositions());
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((Particles->GetRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	}
      NextCoordinates = (int) (((double) NbrParticles) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrParticles) --NextCoordinates;      
    }
  this->AdaptNorm(Particles->GetPositions());
  double SumTrialValues=0.0;
  for (int i = 0; i < average; ++i)
    {
      Particles->Move(NextCoordinates);
      TmpMetropolis = (*this)(Particles->GetPositions());
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((Particles->GetRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	}
      NextCoordinates = (int) (((double) NbrParticles) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrParticles) --NextCoordinates;            
      SumTrialValues+=Norm(TmpMetropolis);
    }  
  this->ElementNorm*= pow(SumTrialValues/average,(double)-2.0/this->NbrParticles);
  cout << "averaged norm:"<<SumTrialValues/average<<endl;
  delete Particles;
}


// this is the main part of the calculation of the paired wavefunction:
void PairedCFOnSphere2QHWaveFunction::EvaluateTables()
{
  int i, j, offset, alpha;
  Complex tmp;
  // evaluate single particle Jastrow factors
  this->Interpolation=1.0;
  if (Orbitals->TestCriticality(Interpolation) == 0)
    {
      for (i=0;i<this->NbrParticles;i++)
	{
	  Ji[i]=1.0;
	  for(j=0;j<i;j++) Ji[i] *= Orbitals->JastrowFactorElement(i,j);
	  for(j=i+1;j<NbrParticles;j++) Ji[i] *= Orbitals->JastrowFactorElement(i,j);
	}
    }
  else // if some interpolation occurred, the true values of the Ji's have to be recalculated:
    {
      for (i=0;i<this->NbrParticles;i++)
	{
	  Ji[i]=1.0;
	  for(j=0;j<i;j++) Ji[i] *= ((Orbitals->SpinorU(i) * Orbitals->SpinorV(j)) - (Orbitals->SpinorU(j) * Orbitals->SpinorV(i)));
	  for(j=i+1;j<NbrParticles;j++) Ji[i] *= ((Orbitals->SpinorU(i) * Orbitals->SpinorV(j)) - (Orbitals->SpinorU(j) * Orbitals->SpinorV(i)));
	}
    }

  // evaluate sums over orbitals m for each LL:
  for (i=0;i<this->NbrParticles;i++)
    for(j=0;j<i;j++)
      {
	alpha=0;
	for (int n=0;n<this->NbrLandauLevels;n++)
	  {
	    tmp=0.0;
	    offset=2*n*(n+this->AbsEffectiveFlux+1)+this->AbsEffectiveFlux;	    
	    for (int m2=-AbsEffectiveFlux-2*n; m2<=AbsEffectiveFlux+2*n;m2+=2)
	      {
		//offset-alpha gives Phi[] with -m 
		tmp+=this->fsgn((m2+AbsEffectiveFlux)/2)*OrbitalValues[alpha][i]*OrbitalValues[offset-alpha][j];
		//cout << "matching up " << alpha << " with " << offset-alpha<<" sign: "<<this->fsgn((m2+AbsEffectiveFlux)/2)<<endl;
		alpha++;
	      }
	    this->gAlpha[n][i*this->NbrParticles+j] = tmp;
	  }
      }
}


// finish evaluation of the function after preparing all arrays with preceding function calls
// integrating over quasihole coordinates by Monte-Carlo
Complex PairedCFOnSphere2QHWaveFunction::MonteCarloEvaluate()
{
  Complex *Uqh, *Vqh, *NewBasis, Tmp;
  Complex Result=0.0;
  int toMove;
  int *NewMValues;
  RealVector tmpVector(2);
  Complex tmpC;
  Complex *AllBasis = new Complex[QHLzMax+1];
  QuasiParticles->GetSpinorCoordinates(Uqh, Vqh);
  for (int step=0; step<QuasiParticleMCSteps; ++step)
    {
      // integrating with sampling function 1 for the moment
      toMove = (int)(2.0*QuasiParticles->GetRandomNumber());
      if (toMove == 0)
	{
	  NewBasis=QH1BasisValues;
	  NewMValues=M1Values;
	}
      else
	{
	  toMove=1;
	  NewBasis=QH2BasisValues;
	  NewMValues=M2Values;
	}
      QuasiParticles->Move(toMove);
      tmpVector[0]=QuasiParticles->Theta(toMove);
      tmpVector[1]=QuasiParticles->Phi(toMove);
      // cout << "Uqh0="<<Uqh[0]<<" Vqh0="<<Vqh[0]<<endl<<"Uqh1="<<Uqh[1]<<"Vqh1="<<Vqh[1]<<" toMove="<<toMove<<endl;
      // QHBasis->GetAllFunctionValues(Uqh[toMove], Vqh[toMove], AllBasis);
      for (int i=0; i<NbrCouplings; ++i)
	{
	  // select one of the following lines
	  // QHBasis->GetFunctionValue(tmpVector, NewBasis[i], NewMValues[i]);
	  QHBasis->GetFunctionValue(Uqh[toMove], Vqh[toMove],NewBasis[i],NewMValues[i]);
	}

      
      // initialize Pfaffian part:
      for (int i=0;i<this->NbrParticles;++i)
	{
	  this->Jqh1[i]=(Orbitals->SpinorU(i) * Vqh[0]) - (Uqh[0] * Orbitals->SpinorV(i));
	  this->Jqh2[i]=(Orbitals->SpinorU(i) * Vqh[1]) - (Uqh[1] * Orbitals->SpinorV(i));
	}
      for (int i=0;i<this->NbrParticles;++i)
	for(int j=0;j<i;++j)
	  {
	    Gij->GetMatrixElement(i, j, Tmp);
	    Slater->SetMatrixElement(i,j, Tmp * (Jqh1[i]*Jqh2[j]-Jqh1[j]*Jqh2[i]));
	  }
      Tmp=0.0;
      // calculate effective quasihole wavefunction
      for (int i=0; i<NbrCouplings; ++i)
	Tmp+=Couplings[i]*QH1BasisValues[i]*QH2BasisValues[i];
      
      Result+=Tmp*Slater->Pfaffian();
    }
  delete [] AllBasis;
  return Result*Interpolation/QuasiParticleMCSteps;
}
