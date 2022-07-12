////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//           class implementing a paired CF wave function on the sphere          //
//                                                                            //
//                        last modification : 18/05/2007                      //
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
#include "PairedCFOnSphereWithSpinWaveFunction.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
//

PairedCFOnSphereWithSpinWaveFunction::PairedCFOnSphereWithSpinWaveFunction()
{
  this->NbrParticlesPerLayer = 0;
  this->NbrParameters = 0;
}

  // constructor
  //
  // nbrParticles = number of particles
  // nbrLandauLevel = number of Landau levels filled with composite fermions
  // nbrEffectiveFlux = number of flux quanta of the magnetic monopole field experienced by CF's
  // haveBosons = indicates whether bosons are present
  // bosonCoefficient = prefactor of singular 1/z term in pair-wave function
  // CFCoefficients = prefactors of CF orbitals in shells 0, 1, 2, ... , nbrLandauLevel-1
  // correctPrefactors = flag that enables the correction of prefactors to adopt the conventions of previous code
  // jastrowPower = power to which the Jastrow factor has to be raised


PairedCFOnSphereWithSpinWaveFunction::PairedCFOnSphereWithSpinWaveFunction(int nbrParticles, int nbrLandauLevels,
									   int nbrEffectiveFlux, bool haveBosons, 
							   double bosonCoefficient, double * givenCFCoefficients,
							   bool correctPrefactors, int jastrowPower)
{
  if (nbrParticles&1)
    {
      cout << "Error: Paired bilayer states require an even number of particles!"<<endl;
      exit(1);
    }
  this->NbrParticlesPerLayer = nbrParticles/2;
  this->NbrLandauLevels = nbrLandauLevels;
  this->NbrParameters = this->NbrLandauLevels+1; // inherited field
  this->AbsEffectiveFlux = abs(nbrEffectiveFlux);
  this->HaveBosons = haveBosons;
  this->Orbitals1 = new JainCFOnSphereOrbitals(NbrParticlesPerLayer, nbrLandauLevels, nbrEffectiveFlux, jastrowPower);
  this->Orbitals2 = new JainCFOnSphereOrbitals(NbrParticlesPerLayer, nbrLandauLevels, nbrEffectiveFlux, jastrowPower);
  this->BosonCoefficient=bosonCoefficient;
  this->TrialParameters= new double [NbrParameters];
  for (int i=0; i<NbrLandauLevels; ++i) this->TrialParameters[i] = givenCFCoefficients[i];
  this->TrialParameters[this->NbrLandauLevels] = BosonCoefficient;
  this->ElementNorm=1.0;
  this->Flag.Initialize();
  
#ifdef USE_LAPACK_CFCB
  this->Matrix = new ComplexLapackDeterminant(this->NbrParticlesPerLayer);
#else
  this->Matrix = new ComplexMatrix(this->NbrParticlesPerLayer,this->NbrParticlesPerLayer);
#endif
  
  this->J11 = new Complex[this->NbrParticlesPerLayer];
  this->J12 = new Complex[this->NbrParticlesPerLayer];
  this->J21 = new Complex[this->NbrParticlesPerLayer];
  this->J22 = new Complex[this->NbrParticlesPerLayer];
  this->gAlpha = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha[i]= new Complex[NbrParticlesPerLayer*NbrParticlesPerLayer];      
  
  this->InterLayerDistances = new Complex*[this->NbrParticlesPerLayer];
  for (int i=0; i< this->NbrParticlesPerLayer; ++i)
    InterLayerDistances[i]= new Complex[NbrParticlesPerLayer];
  
  if (correctPrefactors) 
    {
      cout << "Option correctPrefactors in PairedCFOnSphereWithSpinWaveFunction needs to be checked!" << endl;
//       int p=Orbitals1->GetJastrowPower();
//       FactorialCoefficient Coef;
//       for (int n=0; n<NbrLandauLevels; ++n)
// 	{
// 	  Coef.SetToOne();
// 	  Coef.PartialFactorialDivide(nbrEffectiveFlux+p*(NbrParticlesPerLayer-1)+2,nbrEffectiveFlux+2*p*(NbrParticlesPerLayer-1)+1);
// 	  Coef.PartialFactorialMultiply(p*(NbrParticlesPerLayer-1)+n+2,2*p*(NbrParticlesPerLayer-1)+n+1);
// 	  this->TrialParameters[n]*=Coef.GetNumericalValue()*Coef.GetNumericalValue();
// 	  cout << "Correction["<<n<<"]="<<Coef.GetNumericalValue()<<"from: " << nbrEffectiveFlux+p*(NbrParticlesPerLayer-1)+2<<","<<nbrEffectiveFlux+2*p*(NbrParticlesPerLayer-1)+1<<","<< p*(NbrParticlesPerLayer-1)+n+2<<","<<2*p*(NbrParticlesPerLayer-1)+n+1<<endl;
//	}
    }  
}

// copy constructor
//
// function = reference on the wave function to copy

PairedCFOnSphereWithSpinWaveFunction::PairedCFOnSphereWithSpinWaveFunction(const PairedCFOnSphereWithSpinWaveFunction& function)
{
  this->NbrParticlesPerLayer = function.NbrParticlesPerLayer;
  this->NbrLandauLevels = function.NbrLandauLevels;
  this->NbrParameters = function.NbrParameters;
  this->AbsEffectiveFlux = function.AbsEffectiveFlux;
  this->HaveBosons = function.HaveBosons;
  this->Flag = function.Flag;
  this->Orbitals1 = function.Orbitals1;
  this->Orbitals2 = function.Orbitals2;
  this->BosonCoefficient=function.BosonCoefficient;
  this->TrialParameters=function.TrialParameters;
  this->ElementNorm=function.ElementNorm;

#ifdef USE_LAPACK_CFCB
  this->Matrix = new ComplexLapackDeterminant(this->NbrParticlesPerLayer);
#else
  this->Matrix = new ComplexMatrix(this->NbrParticlesPerLayer,this->NbrParticlesPerLayer);
#endif
  
  this->J11 = new Complex[this->NbrParticlesPerLayer];
  this->J12 = new Complex[this->NbrParticlesPerLayer];
  this->J21 = new Complex[this->NbrParticlesPerLayer];
  this->J22 = new Complex[this->NbrParticlesPerLayer];
  
  this->gAlpha = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha[i]= new Complex[NbrParticlesPerLayer*NbrParticlesPerLayer];
  
  this->InterLayerDistances = new Complex*[this->NbrParticlesPerLayer];
  for (int i=0; i< this->NbrParticlesPerLayer; ++i)
    InterLayerDistances[i]= new Complex[NbrParticlesPerLayer];
  
}

// destructor
//

PairedCFOnSphereWithSpinWaveFunction::~PairedCFOnSphereWithSpinWaveFunction()
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
  for (int i=0; i< this->NbrParticlesPerLayer; ++i)
    delete [] InterLayerDistances[i];
  delete [] InterLayerDistances;
  delete [] J11;
  delete [] J12;
  delete [] J21;
  delete [] J22;
  delete Matrix;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PairedCFOnSphereWithSpinWaveFunction::Clone ()
{
  return new PairedCFOnSphereWithSpinWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PairedCFOnSphereWithSpinWaveFunction::operator ()(RealVector& x)
{
  RealVector part=x.Extract(0, 2*this->NbrParticlesPerLayer-1);
  this->OrbitalValues1 = (*Orbitals1)(part);
  part = x.Extract(2*this->NbrParticlesPerLayer, 4*this->NbrParticlesPerLayer-1);
  this->OrbitalValues2 = (*Orbitals2)(part);
  this->EvaluateTables();
  Complex tmp;

  if (this->HaveBosons)
    {
      // initialize Slater determinant 
      for (int i=0;i<this->NbrParticlesPerLayer;++i)
	{
	  for(int j=0;j<this->NbrParticlesPerLayer;++j)
	    {
	      tmp=0.0;
	      for (int n=0; n<this->NbrLandauLevels; ++n)
		tmp+=this->TrialParameters[n]*this->gAlpha[n][i*this->NbrParticlesPerLayer+j];
	      tmp = this->ElementNorm*( this->BosonCoefficient*J12[i]*J21[j]/InterLayerDistances[i][j] - tmp*J11[i]*J22[j]);
#ifdef USE_LAPACK_CFCB
	      Matrix->SetMatrixElement(i,j, Real(tmp), Imag(tmp));
#else
	      (*Matrix)[i].Re(j) = Real(tmp);
	      (*Matrix)[i].Im(j) = Imag(tmp);
#endif
	    }
	}
      //cout << *Matrix << endl;
      tmp= Matrix->Determinant()*Interpolation;
    }
  else // no Bosons:
    {
      // initialize Slater determinant 
      for (int i=0;i<this->NbrParticlesPerLayer;++i)
	{
	  for(int j=0;j<this->NbrParticlesPerLayer;++j)
	    {
	      tmp=0.0;
	      for (int n=0; n<this->NbrLandauLevels; ++n)
		tmp+=this->TrialParameters[n]*this->gAlpha[n][i*this->NbrParticlesPerLayer+j];
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
    }
  return tmp;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex PairedCFOnSphereWithSpinWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  ComplexVector part=uv.Extract(0, 2*this->NbrParticlesPerLayer-1);
  this->OrbitalValues1 = Orbitals1->CalculateFromSpinorVariables(part);
  part = uv.Extract(2*this->NbrParticlesPerLayer, 4*this->NbrParticlesPerLayer-1);
  this->OrbitalValues2 = Orbitals2->CalculateFromSpinorVariables(part);
  
  this->EvaluateTables();
  Complex tmp;

  if (this->HaveBosons)
    {
      // initialize Slater determinant 
      for (int i=0;i<this->NbrParticlesPerLayer;++i)
	{
	  for(int j=0;j<this->NbrParticlesPerLayer;++j)
	    {
	      tmp=0.0;
	      for (int n=0; n<this->NbrLandauLevels; ++n)
		tmp+=this->TrialParameters[n]*this->gAlpha[n][i*this->NbrParticlesPerLayer+j];
	      tmp = this->ElementNorm*( this->BosonCoefficient*J12[i]*J21[j]/InterLayerDistances[i][j] - tmp*J11[i]*J22[j]);
#ifdef USE_LAPACK_CFCB
	      Matrix->SetMatrixElement(i,j, Real(tmp), Imag(tmp));
#else
	      (*Matrix)[i].Re(j) = Real(tmp);
	      (*Matrix)[i].Im(j) = Imag(tmp);
#endif
	    }
	}
      //cout << *Matrix << endl;
      tmp= Matrix->Determinant()*Interpolation;
    }
  else // no Bosons:
    {
      // initialize Slater determinant 
      for (int i=0;i<this->NbrParticlesPerLayer;++i)
	{
	  for(int j=0;j<this->NbrParticlesPerLayer;++j)
	    {
	      tmp=0.0;
	      for (int n=0; n<this->NbrLandauLevels; ++n)
		tmp+=this->TrialParameters[n]*this->gAlpha[n][i*this->NbrParticlesPerLayer+j];
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
    }
  return tmp;
}

// get a value of the wavefunction for the last set of coordinates, but with different variational coefficients
// coefficients: array of variational coefficients f_0, f_1, ...
//                    length has to be identical to the initial set of parameters!
// singular: new prefactor of the Moore-Read part
//
Complex PairedCFOnSphereWithSpinWaveFunction::GetForOtherParameters( double *coefficients)
{
  Complex tmp;	      
    if (this->HaveBosons)
    {
      // initialize Slater determinant 
      for (int i=0;i<this->NbrParticlesPerLayer;++i)
	{
	  for(int j=0;j<this->NbrParticlesPerLayer;++j)
	    {
	      tmp=0.0;
	      for (int n=0; n<this->NbrLandauLevels; ++n)
		tmp+=coefficients[n]*this->gAlpha[n][i*this->NbrParticlesPerLayer+j];
	      tmp = this->ElementNorm *(J12[i]*J21[j]/InterLayerDistances[i][j] - tmp*J11[i]*J22[j]);
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
    }
  else // no Bosons:
    {
      // initialize Slater determinant 
      for (int i=0;i<this->NbrParticlesPerLayer;++i)
	{
	  for(int j=0;j<this->NbrParticlesPerLayer;++j)
	    {
	      tmp=0.0;
	      for (int n=0; n<this->NbrLandauLevels; ++n)
		tmp+=coefficients[n]*this->gAlpha[n][i*this->NbrParticlesPerLayer+j];
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
    }
  return tmp;
}


// do many evaluations, storing the result in the vector results given in the call
// x: positions to evaluate the wavefuntion in
// format for passing parameters in the matrix coefficients: coefficients[nbrSet][LandauLevel],
// the entry [][NbrLandauLevels] corresponds to the MooreRead Term.
void PairedCFOnSphereWithSpinWaveFunction::GetForManyParameters(ComplexVector &results, RealVector& x, double **coefficients)
{
  RealVector part=x.Extract(0, 2*this->NbrParticlesPerLayer-1);
  this->OrbitalValues1 = (*Orbitals1)(part);
  part = x.Extract(2*this->NbrParticlesPerLayer, 4*this->NbrParticlesPerLayer-1);
  this->OrbitalValues2 = (*Orbitals2)(part);
  this->EvaluateTables();
  Complex tmp;
  int numParamSets=results.GetVectorDimension();
  double *tmpCoefficients;
  for (int s=0; s<numParamSets; ++s)
    {
      tmpCoefficients = coefficients[s];
      if (this->HaveBosons)
	{
	  // initialize Slater determinant 
	  for (int i=0;i<this->NbrParticlesPerLayer;++i)
	    {
	      for(int j=0;j<this->NbrParticlesPerLayer;++j)
		{
		  tmp=0.0;
		  for (int n=0; n<this->NbrLandauLevels; ++n)
		    tmp+= tmpCoefficients[n]*this->gAlpha[n][i*this->NbrParticlesPerLayer+j];
		  tmp = this->ElementNorm *(J12[i]*J21[j]/InterLayerDistances[i][j] - tmp*J11[i]*J22[j]);
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
	}
      else // no Bosons:
	{
	  // initialize Slater determinant 
	  for (int i=0;i<this->NbrParticlesPerLayer;++i)
	    {
	      for(int j=0;j<this->NbrParticlesPerLayer;++j)
		{
		  tmp=0.0;
		  for (int n=0; n<this->NbrLandauLevels; ++n)
		    tmp+=tmpCoefficients[n]*this->gAlpha[n][i*this->NbrParticlesPerLayer+j];
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
	}
      results.Re(s)=Real(tmp);
      results.Im(s)=Imag(tmp);
    }
}


void PairedCFOnSphereWithSpinWaveFunction::SetTrialParameters(double * coefficients)
{
  for (int n=0; n<this->NbrParameters; ++n)
    this->TrialParameters[n]=coefficients[n];
  this->BosonCoefficient =  this->TrialParameters[this->NbrLandauLevels];
}

// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void PairedCFOnSphereWithSpinWaveFunction::AdaptNorm(RealVector& x)
{
  double det=Norm((*this)(x));
  while ((det<.1)||(det>50.0))
    {
      //cout <<"N'="<< this->ElementNorm << " det="<<det<<endl;
      if (det>1e300) 
	this->ElementNorm*= pow((double)1.0e-300,(double)1.0/this->NbrParticlesPerLayer);
      else if (det==0.0) 
	this->ElementNorm*= pow((double)1.0e300,(double)1.0/this->NbrParticlesPerLayer);
      else 
	this->ElementNorm*= pow(det,(double)-1.0/this->NbrParticlesPerLayer);
      det=Norm((*this)(x));
      //cout <<"N'="<< this->ElementNorm << endl;
    }
}


// normalize the wave-function over an average number of MC positions

void PairedCFOnSphereWithSpinWaveFunction::AdaptAverageMCNorm(int thermalize, int average)
{
  ParticleOnSphereCollection * Particles = new ParticleOnSphereCollection(2*this->NbrParticlesPerLayer);
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
      NextCoordinates = (int) (((double) NbrParticlesPerLayer) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrParticlesPerLayer) --NextCoordinates;      
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
      NextCoordinates = (int) (((double) NbrParticlesPerLayer) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrParticlesPerLayer) --NextCoordinates;            
      SumTrialValues+=Norm(TmpMetropolis);
    }  
  this->ElementNorm*= pow(SumTrialValues/average,(double)-2.0/this->NbrParticlesPerLayer);
  delete Particles;
}


// this is the main part of the calculation of the paired wavefunction:
// assumes that Orbitals1 and Orbitals2 have already been initialized prior to call
void PairedCFOnSphereWithSpinWaveFunction::EvaluateTables()
{
  int offset, alpha;  
  Complex tmp;
  // evaluate single particle Jastrow factors
  this->Interpolation=1.0;
  if (Orbitals1->TestCriticality(Interpolation) == 0)
    {
      for (int i=0;i<this->NbrParticlesPerLayer;i++)
	{
	  J11[i]=1.0;
	  for(int j=0;j<i;j++) J11[i] *= Orbitals1->JastrowFactorElement(i,j);
	  for(int j=i+1;j<NbrParticlesPerLayer;j++) J11[i] *= Orbitals1->JastrowFactorElement(i,j);
	}
    }
  else // if some interpolation occurred, the true values of the J11's have to be recalculated:
    {
      for (int i=0;i<this->NbrParticlesPerLayer;i++)
	{
	  J11[i]=1.0;
	  for(int j=0;j<i;j++) J11[i] *= ((Orbitals1->SpinorU(i) * Orbitals1->SpinorV(j)) - (Orbitals1->SpinorU(j) * Orbitals1->SpinorV(i)));
	  for(int j=i+1;j<NbrParticlesPerLayer;j++) J11[i] *= ((Orbitals1->SpinorU(i) * Orbitals1->SpinorV(j)) - (Orbitals1->SpinorU(j) * Orbitals1->SpinorV(i)));
	}
    }

  if (Orbitals2->TestCriticality(Interpolation) == 0)
    {
      for (int i=0;i<this->NbrParticlesPerLayer;i++)
	{
	  J22[i]=1.0;
	  for(int j=0;j<i;j++) J22[i] *= Orbitals2->JastrowFactorElement(i,j);
	  for(int j=i+1;j<NbrParticlesPerLayer;j++) J22[i] *= Orbitals2->JastrowFactorElement(i,j);
	}
    }
  else // if some interpolation occurred, the true values of the J22's have to be recalculated:
    {
      for (int i=0;i<this->NbrParticlesPerLayer;i++)
	{
	  J22[i]=1.0;
	  for(int j=0;j<i;j++) J22[i] *= ((Orbitals2->SpinorU(i) * Orbitals2->SpinorV(j)) - (Orbitals2->SpinorU(j) * Orbitals2->SpinorV(i)));
	  for(int j=i+1;j<NbrParticlesPerLayer;j++) J22[i] *= ((Orbitals2->SpinorU(i) * Orbitals2->SpinorV(j)) - (Orbitals2->SpinorU(j) * Orbitals2->SpinorV(i)));
	}
    }

  if (this->HaveBosons) // need to calculate additional terms, then:
    {
      for (int i=0; i<this->NbrParticlesPerLayer; ++i)
	for (int j=0; j<this->NbrParticlesPerLayer; ++j)
	  InterLayerDistances[i][j] = ((Orbitals1->SpinorU(i) * Orbitals2->SpinorV(j)) - (Orbitals2->SpinorU(j) * Orbitals1->SpinorV(i)));

      for (int i=0;i<this->NbrParticlesPerLayer;++i)
	{
	  J12[i]=1.0;
	  J21[i]=1.0;
	  for(int j=0;j<this->NbrParticlesPerLayer;++j)
	    {
	      J12[i] *= InterLayerDistances[i][j];
	      J21[i] *= -InterLayerDistances[j][i];
	    }
	}
    }


  // evaluate sums over orbitals m for each LL:
  for (int i=0;i<this->NbrParticlesPerLayer;i++)
    for(int j=0;j<this->NbrParticlesPerLayer;j++)
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

		//cout << "matching up " << alpha << " with " << offset-alpha<<" sign: "<<this->fsgn((m2+AbsEffectiveFlux)/2)<<endl;
		alpha++;
	      }
	    this->gAlpha[n][i*this->NbrParticlesPerLayer+j] = tmp;
	  }
      }
}
