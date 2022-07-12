////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                      Copyright (C) 2009 Gunnar Möller                      //
//                                                                            //
//                                                                            //
//        class implementing a paired CF wave function on the sphere          //
//                                                                            //
//                        last modification : 04/25/2009                      //
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
#include "PairedCFOnSphereSpinSingletWaveFunction.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
//

PairedCFOnSphereSpinSingletWaveFunction::PairedCFOnSphereSpinSingletWaveFunction()
{
  this->NbrUp = 0;
  this->NbrParameters = 0;
}

// constructor
//
// nbrParticles= total number of particles 
// nbrLandauLevel = number of Landau levels filled with composite fermions
// nbrEffectiveFlux = number of flux quanta of the magnetic monopole field experienced by CF's
// cfCoefficients = prefactors of CF orbitals in shells 0, 1, 2, ... , nbrLandauLevel-1
// jastrowPower = power to which the Jastrow factor has to be raised
//
PairedCFOnSphereSpinSingletWaveFunction::PairedCFOnSphereSpinSingletWaveFunction(int nbrParticles,
					  int nbrLandauLevels, int nbrEffectiveFlux,
					  double * cfCoefficients, int jastrowPower)
{
  if (nbrParticles&1)
    {
      cout << "Error: Paired spin-singlet states require an even number of particles!"<<endl;
      exit(1);
    }
  this->NbrUp = nbrParticles/2;
  this->NbrLandauLevels = nbrLandauLevels;
  this->NbrParameters = this->NbrLandauLevels; // inherited field
  this->AbsEffectiveFlux = abs(nbrEffectiveFlux);
  this->JastrowPower=jastrowPower;
  this->Orbitals1 = new JainCFOnSphereOrbitals(NbrUp, nbrLandauLevels, nbrEffectiveFlux, jastrowPower);
  this->Orbitals2 = new JainCFOnSphereOrbitals(NbrUp, nbrLandauLevels, nbrEffectiveFlux, jastrowPower);
  this->TrialParameters= new double [NbrParameters];
  for (int i=0; i<NbrLandauLevels; ++i) this->TrialParameters[i] = cfCoefficients[i];
  this->ElementNorm=1.0;
  this->Flag.Initialize();
  
#ifdef USE_LAPACK_CFCB
  this->Matrix = new ComplexLapackDeterminant(this->NbrUp);
#else
  this->Matrix = new ComplexMatrix(this->NbrUp,this->NbrUp);
#endif
  
  this->J11 = new Complex[this->NbrUp];
  this->J12 = new Complex[this->NbrUp];
  this->J21 = new Complex[this->NbrUp];
  this->J22 = new Complex[this->NbrUp];
  this->gAlpha = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha[i]= new Complex[NbrUp*NbrUp];      
 
}

// copy constructor
//
// function = reference on the wave function to copy

PairedCFOnSphereSpinSingletWaveFunction::PairedCFOnSphereSpinSingletWaveFunction(const PairedCFOnSphereSpinSingletWaveFunction& function)
{
  this->NbrUp = function.NbrUp;
  this->NbrLandauLevels = function.NbrLandauLevels;
  this->NbrParameters = function.NbrParameters;
  this->AbsEffectiveFlux = function.AbsEffectiveFlux;
  this->JastrowPower=function.JastrowPower;
  this->Flag = function.Flag;  
  this->Orbitals1 = function.Orbitals1;
  this->Orbitals2 = function.Orbitals2;
  this->TrialParameters=function.TrialParameters;
  this->ElementNorm=function.ElementNorm;

#ifdef USE_LAPACK_CFCB
  this->Matrix = new ComplexLapackDeterminant(this->NbrUp);
#else
  this->Matrix = new ComplexMatrix(this->NbrUp,this->NbrUp);
#endif
  
  this->J11 = new Complex[this->NbrUp];
  this->J12 = new Complex[this->NbrUp];
  this->J21 = new Complex[this->NbrUp];
  this->J22 = new Complex[this->NbrUp];
  
  this->gAlpha = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha[i]= new Complex[NbrUp*NbrUp];
  
}

// destructor
//

PairedCFOnSphereSpinSingletWaveFunction::~PairedCFOnSphereSpinSingletWaveFunction()
{
  if ( (this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete Orbitals1;
      delete Orbitals2;
      delete [] TrialParameters;
    }
  for (int i=0; i< this->NbrLandauLevels; ++i)
    delete [] gAlpha[i];
  delete [] gAlpha;
  delete [] J11;
  delete [] J12;
  delete [] J21;
  delete [] J22;
  delete Matrix;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PairedCFOnSphereSpinSingletWaveFunction::Clone ()
{
  return new PairedCFOnSphereSpinSingletWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PairedCFOnSphereSpinSingletWaveFunction::operator ()(RealVector& x)
{
  RealVector part=x.Extract(0, 2*this->NbrUp-1);
  this->OrbitalValues1 = (*Orbitals1)(part);
  part = x.Extract(2*this->NbrUp, 4*this->NbrUp-1);
  this->OrbitalValues2 = (*Orbitals2)(part);
  this->EvaluateTables();
  Complex tmp;

  // initialize Slater determinant 
  for (int i=0;i<this->NbrUp;++i)
    {
      for(int j=0;j<this->NbrUp;++j)
	{
	  tmp=0.0;
	  for (int n=0; n<this->NbrLandauLevels; ++n)
	    tmp+=this->TrialParameters[n]*this->gAlpha[n][i*this->NbrUp+j];
	  tmp *=  - this->ElementNorm * J11[i]*J22[j];
#ifdef USE_LAPACK_CFCB
	  Matrix->SetMatrixElement(i,j,Real(tmp), Imag(tmp));
#else
	  (*Matrix)[i].Re(j) = Real(tmp);
	  (*Matrix)[i].Im(j) = Imag(tmp);
#endif
	}
    }
  //cout << *Matrix << endl;
  tmp = Matrix->Determinant()*Interpolation;

  return tmp*InterSpinJastrow;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex PairedCFOnSphereSpinSingletWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  ComplexVector part=uv.Extract(0, 2*this->NbrUp-1);
  this->OrbitalValues1 = Orbitals1->CalculateFromSpinorVariables(part);
  part = uv.Extract(2*this->NbrUp, 4*this->NbrUp-1);
  this->OrbitalValues2 = Orbitals2->CalculateFromSpinorVariables(part);
  
  this->EvaluateTables();
  Complex tmp;

  // initialize Slater determinant 
  for (int i=0;i<this->NbrUp;++i)
    {
      for(int j=0;j<this->NbrUp;++j)
	{
	  tmp=0.0;
	  for (int n=0; n<this->NbrLandauLevels; ++n)
	    tmp+=this->TrialParameters[n]*this->gAlpha[n][i*this->NbrUp+j];
	  tmp *=  - this->ElementNorm * J11[i]*J22[j];
#ifdef USE_LAPACK_CFCB
	  Matrix->SetMatrixElement(i,j,Real(tmp), Imag(tmp));
#else
	  (*Matrix)[i].Re(j) = Real(tmp);
	  (*Matrix)[i].Im(j) = Imag(tmp);
#endif
	}
    }
  //cout << *Matrix << endl;
  tmp = Matrix->Determinant()*Interpolation;
  
  return tmp*InterSpinJastrow;
}

// get a value of the wavefunction for the last set of coordinates, but with different variational coefficients
// coefficients: array of variational coefficients f_0, f_1, ...
//                    length has to be identical to the initial set of parameters!
// singular: new prefactor of the Moore-Read part
//
Complex PairedCFOnSphereSpinSingletWaveFunction::GetForOtherParameters( double *coefficients)
{
  Complex tmp;	      

  // initialize Slater determinant 
  for (int i=0;i<this->NbrUp;++i)
    {
      for(int j=0;j<this->NbrUp;++j)
	{
	  tmp=0.0;
	  for (int n=0; n<this->NbrLandauLevels; ++n)
	    tmp+=coefficients[n]*this->gAlpha[n][i*this->NbrUp+j];
	  tmp *=  - this->ElementNorm * J11[i]*J22[j];
#ifdef USE_LAPACK_CFCB
	  Matrix->SetMatrixElement(i,j,Real(tmp), Imag(tmp));
#else
	  (*Matrix)[i].Re(j) = Real(tmp);
	  (*Matrix)[i].Im(j) = Imag(tmp);
#endif
	}
    }
  //cout << *Matrix << endl;
  tmp= Matrix->Determinant()*Interpolation;
  
  return tmp*InterSpinJastrow;
}


// do many evaluations, storing the result in the vector results given in the call
// x: positions to evaluate the wavefuntion in
// format for passing parameters in the matrix coefficients: coefficients[nbrSet][LandauLevel],
// the entry [][NbrLandauLevels] corresponds to the MooreRead Term.
void PairedCFOnSphereSpinSingletWaveFunction::GetForManyParameters(ComplexVector &results, RealVector& x, double **coefficients)
{
  RealVector part=x.Extract(0, 2*this->NbrUp-1);
  this->OrbitalValues1 = (*Orbitals1)(part);
  part = x.Extract(2*this->NbrUp, 4*this->NbrUp-1);
  this->OrbitalValues2 = (*Orbitals2)(part);
  this->EvaluateTables();
  Complex tmp;
  int numParamSets=results.GetVectorDimension();
  double *tmpCoefficients;
  for (int s=0; s<numParamSets; ++s)
    {
      tmpCoefficients = coefficients[s];
      
      // initialize Slater determinant 
      for (int i=0;i<this->NbrUp;++i)
	{
	  for(int j=0;j<this->NbrUp;++j)
	    {
	      tmp=0.0;
	      for (int n=0; n<this->NbrLandauLevels; ++n)
		tmp+=tmpCoefficients[n]*this->gAlpha[n][i*this->NbrUp+j];
	      tmp *=  - this->ElementNorm * J11[i]*J22[j];
#ifdef USE_LAPACK_CFCB
	      Matrix->SetMatrixElement(i,j,Real(tmp), Imag(tmp));
#else
	      (*Matrix)[i].Re(j) = Real(tmp);
	      (*Matrix)[i].Im(j) = Imag(tmp);
#endif
	    }
	}
      //cout << *Matrix << endl;
      tmp = Matrix->Determinant()*Interpolation*InterSpinJastrow;

      results.Re(s)=Real(tmp);
      results.Im(s)=Imag(tmp);
    }
}


void PairedCFOnSphereSpinSingletWaveFunction::SetTrialParameters(double * coefficients)
{
  for (int n=0; n<this->NbrParameters; ++n)
    this->TrialParameters[n]=coefficients[n];
}


// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void PairedCFOnSphereSpinSingletWaveFunction::AdaptNorm(RealVector& x)
{
  double det=Norm((*this)(x));
  while ((det<.1)||(det>50.0))
    {
      //cout <<"N'="<< this->ElementNorm << " det="<<det<<endl;
      if (det>1e300) 
	this->ElementNorm*= pow((double)1.0e-300,(double)1.0/this->NbrUp);
      else if (det==0.0) 
	this->ElementNorm*= pow((double)1.0e300,(double)1.0/this->NbrUp);
      else 
	this->ElementNorm*= pow(det,(double)-1.0/this->NbrUp);
      det=Norm((*this)(x));
      //cout <<"N'="<< this->ElementNorm << endl;
    }
}


// normalize the wave-function over an average number of MC positions

void PairedCFOnSphereSpinSingletWaveFunction::AdaptAverageMCNorm(int thermalize, int average)
{
  ParticleOnSphereCollection * Particles = new ParticleOnSphereCollection(2*this->NbrUp);
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
      NextCoordinates = (int) (((double) NbrUp) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrUp) --NextCoordinates;      
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
      NextCoordinates = (int) (((double) NbrUp) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrUp) --NextCoordinates;            
      SumTrialValues+=Norm(TmpMetropolis);
    }  
  this->ElementNorm*= pow(SumTrialValues/average,(double)-2.0/this->NbrUp);
  delete Particles;
}


// this is the main part of the calculation of the paired wavefunction:
// assumes that Orbitals1 and Orbitals2 have already been initialized prior to call
void PairedCFOnSphereSpinSingletWaveFunction::EvaluateTables()
{
  int offset, alpha;  
  Complex tmp;
  // evaluate single particle Jastrow factors
  this->Interpolation=1.0;
  if (Orbitals1->TestCriticality(Interpolation) == 0)
    {
      for (int i=0;i<this->NbrUp;i++)
	{
	  J11[i]=1.0;
	  for(int j=0;j<i;j++) J11[i] *= Orbitals1->JastrowFactorElement(i,j);
	  for(int j=i+1;j<NbrUp;j++) J11[i] *= Orbitals1->JastrowFactorElement(i,j);
	}
    }
  else // if some interpolation occurred, the true values of the J11's have to be recalculated:
    {
      for (int i=0;i<this->NbrUp;i++)
	{
	  J11[i]=1.0;
	  for(int j=0;j<i;j++) J11[i] *= ((Orbitals1->SpinorU(i) * Orbitals1->SpinorV(j)) - (Orbitals1->SpinorU(j) * Orbitals1->SpinorV(i)));
	  for(int j=i+1;j<NbrUp;j++) J11[i] *= ((Orbitals1->SpinorU(i) * Orbitals1->SpinorV(j)) - (Orbitals1->SpinorU(j) * Orbitals1->SpinorV(i)));
	}
    }

  if (Orbitals2->TestCriticality(Interpolation) == 0)
    {
      for (int i=0;i<this->NbrUp;i++)
	{
	  J22[i]=1.0;
	  for(int j=0;j<i;j++) J22[i] *= Orbitals2->JastrowFactorElement(i,j);
	  for(int j=i+1;j<NbrUp;j++) J22[i] *= Orbitals2->JastrowFactorElement(i,j);
	}
    }
  else // if some interpolation occurred, the true values of the J22's have to be recalculated:
    {
      for (int i=0;i<this->NbrUp;i++)
	{
	  J22[i]=1.0;
	  for(int j=0;j<i;j++) J22[i] *= ((Orbitals2->SpinorU(i) * Orbitals2->SpinorV(j)) - (Orbitals2->SpinorU(j) * Orbitals2->SpinorV(i)));
	  for(int j=i+1;j<NbrUp;j++) J22[i] *= ((Orbitals2->SpinorU(i) * Orbitals2->SpinorV(j)) - (Orbitals2->SpinorU(j) * Orbitals2->SpinorV(i)));
	}
    }

  double Factor = M_PI * 0.5;
  Complex Base(1.0);
  
  for (int i=0; i<this->NbrUp; ++i)
    for (int j=0; j<this->NbrUp; ++j)
      Base *= Factor*((Orbitals1->SpinorU(i) * Orbitals2->SpinorV(j)) - (Orbitals2->SpinorU(j) * Orbitals1->SpinorV(i)));
  InterSpinJastrow=Base;
  for (int i=1; i<this->JastrowPower; ++i)
    InterSpinJastrow*=Base;
  
  // evaluate sums over orbitals m for each LL:
  for (int i=0;i<this->NbrUp;i++)
    for(int j=0;j<this->NbrUp;j++)
      {
	alpha=0;
	for (int n=0;n<this->NbrLandauLevels;n++)
	  {
	    tmp=0.0;
	    offset=2*n*(n+this->AbsEffectiveFlux+1)+this->AbsEffectiveFlux;	    
	    for (int m2=-AbsEffectiveFlux-2*n; m2<=AbsEffectiveFlux+2*n;m2+=2)
	      {
		//offset-alpha gives Phi[] with -m 
		tmp+=this->fsgn((m2+AbsEffectiveFlux)/2)*OrbitalValues1[alpha][i]*OrbitalValues2[offset-alpha][j];
		//tmp+=Conj(OrbitalValues1[alpha][i])*OrbitalValues2[offset-alpha][j];

		//cout << "matching up " << alpha << " with " << offset-alpha<<" sign: "<<this->fsgn((m2+AbsEffectiveFlux)/2)<<endl;
		alpha++;
	      }
	    this->gAlpha[n][i*this->NbrUp+j] = tmp;
	  }
      }
}
