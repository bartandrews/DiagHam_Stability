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
#include "PairedCFOnSphereWaveFunction.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
//

PairedCFOnSphereWaveFunction::PairedCFOnSphereWaveFunction()
{
  this->NbrParticles = 0;
  this->NbrParameters = 0;
}

  // constructor
  //
  // nbrParticles = number of particles
  // nbrLandauLevel = number of Landau levels filled with composite fermions
  // nbrEffectiveFlux = number of flux quanta of the magnetic monopole field experienced by CF's
  // MooreReadCoefficient = prefactor of singular 1/z term in pair-wave function
  // CFCoefficients = prefactors of CF orbitals in shells 0, 1, 2, ... , nbrLandauLevel-1
  // correctPrefactors = flag that enables the correction of prefactors to adopt the conventions of previous code
  // jastrowPower = power to which the Jastrow factor has to be raised


PairedCFOnSphereWaveFunction::PairedCFOnSphereWaveFunction(int nbrParticles, int nbrLandauLevels, int nbrEffectiveFlux,
							   double MooreReadCoefficient, double * givenCFCoefficients,
							   bool correctPrefactors, int jastrowPower)
{
  this->NbrParticles = nbrParticles;
  this->NbrLandauLevels = nbrLandauLevels;
  this->NbrParameters = this->NbrLandauLevels+1; // inherited field
  this->AbsEffectiveFlux = abs(nbrEffectiveFlux);
  this->AddJastrowPower = (jastrowPower>=2?jastrowPower - 2:0);
  this->Orbitals = new JainCFOnSphereOrbitals(nbrParticles, nbrLandauLevels, nbrEffectiveFlux,jastrowPower);
  this->MooreReadCoefficient=MooreReadCoefficient;
  this->TrialParameters= new double [NbrParameters];
  for (int i=0; i<NbrLandauLevels; ++i) this->TrialParameters[i] = givenCFCoefficients[i];
  this->TrialParameters[this->NbrLandauLevels] = MooreReadCoefficient;
  this->ElementNorm=1.0;
  this->Slater = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->Flag.Initialize();
  this->Ji = new Complex[this->NbrParticles];

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

  cout << "trial parameters in Paired:";
  for (int i=0; i<NbrParameters; ++i) cout <<" "<<TrialParameters[i];
  cout << endl;

}

// copy constructor
//
// function = reference on the wave function to copy

PairedCFOnSphereWaveFunction::PairedCFOnSphereWaveFunction(const PairedCFOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrLandauLevels = function.NbrLandauLevels;
  this->NbrParameters = function.NbrParameters;
  this->AbsEffectiveFlux = function.AbsEffectiveFlux;
  this->Flag = function.Flag;
  this->AddJastrowPower = function.AddJastrowPower;
  this->Orbitals = function.Orbitals;
  this->MooreReadCoefficient=function.MooreReadCoefficient;
  this->TrialParameters=function.TrialParameters;
  this->ElementNorm=function.ElementNorm;
  
  this->Slater = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->Ji = new Complex[this->NbrParticles];
  this->gAlpha = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha[i]= new Complex[NbrParticles*NbrParticles];
  
}

// destructor
//

PairedCFOnSphereWaveFunction::~PairedCFOnSphereWaveFunction()
{
  if ( (this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete Orbitals;
      delete [] TrialParameters;
    }
  for (int i=0; i< this->NbrLandauLevels; ++i) delete [] gAlpha[i];
  delete [] gAlpha;
  delete [] Ji;
  delete Slater;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PairedCFOnSphereWaveFunction::Clone ()
{
  return new PairedCFOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PairedCFOnSphereWaveFunction::operator ()(RealVector& x)
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
	  
	  Slater->SetMatrixElement(i,j, this->ElementNorm*this->Ji[i]*this->Ji[j]
				   *(MooreReadCoefficient/Orbitals->JastrowFactorElement(i,j) + tmp));
	}
    }  
  //cout << *Slater << endl;
  return Slater->Pfaffian()*Interpolation*AdditionalJastrow;
}


// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex PairedCFOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
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
	  
	  Slater->SetMatrixElement(i,j, this->ElementNorm*this->Ji[i]*this->Ji[j]
				   *(MooreReadCoefficient/Orbitals->JastrowFactorElement(i,j) + tmp));
	}
    }  
  //cout << *Slater << endl;
  return Slater->Pfaffian()*Interpolation*AdditionalJastrow;
} 


// get a value of the wavefunction for the last set of coordinates, but with different variational coefficients
// coefficients: array of variational coefficients f_0, f_1, ...
//                    length has to be identical to the initial set of parameters!
// singular: new prefactor of the Moore-Read part
//
Complex PairedCFOnSphereWaveFunction::GetForOtherParameters( double *coefficients)
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
	  Slater->SetMatrixElement(i,j, this->ElementNorm*this->Ji[i]*this->Ji[j]
				   *(coefficients[this->NbrLandauLevels]/Orbitals->JastrowFactorElement(i,j) + tmp));
	}
    }  
  return Slater->Pfaffian()*Interpolation*AdditionalJastrow;
}

// do many evaluations, storing the result in the vector results given in the call
// x: positions to evaluate the wavefuntion in
// format for passing parameters in the matrix coefficients: coefficients[nbrSet][LandauLevel],
// the entry [][NbrLandauLevels] corresponds to the MooreRead Term.
void PairedCFOnSphereWaveFunction::GetForManyParameters(ComplexVector &results, RealVector& x, double **coefficients)
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
	      Slater->SetMatrixElement(i,j, this->ElementNorm*this->Ji[i]*this->Ji[j]
				       *(tmpCoefficients[this->NbrLandauLevels]/Orbitals->JastrowFactorElement(i,j) + tmp));
	    }
	}
      tmp=Slater->Pfaffian()*this->Interpolation*AdditionalJastrow;
      results.Re(s)=Real(tmp);
      results.Im(s)=Imag(tmp);
    }  
}


void PairedCFOnSphereWaveFunction::SetTrialParameters(double * coefficients)
{
  for (int n=0; n<this->NbrParameters; ++n)
    this->TrialParameters[n]=coefficients[n];
  this->MooreReadCoefficient =  this->TrialParameters[this->NbrLandauLevels];
}

// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void PairedCFOnSphereWaveFunction::AdaptNorm(RealVector& x)
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

void PairedCFOnSphereWaveFunction::AdaptAverageMCNorm(int thermalize, int average)
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
  delete Particles;
}


// this is the main part of the calculation of the paired wavefunction:
void PairedCFOnSphereWaveFunction::EvaluateTables()
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
  this->AdditionalJastrow=1.0;
  if (AddJastrowPower)
    {
      Complex Base=1.0;
      for (int i=0; i<this->NbrParticles;++i)
	Base *= Ji[i];
      for (int i=0; i<AddJastrowPower/2; ++i)
	this->AdditionalJastrow *= Base;
    }
}
