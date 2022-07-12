////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                          //
//                                                                            //
//                                                                            //
//   class implementing a CF bilayer wave function with parton tunnelling     //
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
#include "CFOnSphereWithSpinPartonTunnellingWaveFunction.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
//

CFOnSphereWithSpinPartonTunnellingWaveFunction::CFOnSphereWithSpinPartonTunnellingWaveFunction()
{
  this->NbrParticlesUp = 0;
  this->NbrParticlesDown = 0;
  this->NbrParticles = 0;
}

// constructor
//
// nbrParticles = total number of particles
// lFBonding = Fermi momentum in Up - layer
// lFAntiBonding = Fermi momentum in Down - layer
// nbrEffectiveFlux = number of flux quanta of the magnetic monopole field experienced by CF's
// jastrowPower = power to which the Jastrow factor has to be raised
// nbrParticlesUp = number of particles in upper layer (-1 if half filling each)
CFOnSphereWithSpinPartonTunnellingWaveFunction::CFOnSphereWithSpinPartonTunnellingWaveFunction(int nbrParticles, int lFBonding, int lFAntiBonding, int nbrEffectiveFlux, int jastrowPower, int nbrParticlesUp)
{
  if ((lFBonding<0)||(lFAntiBonding<0)||(nbrParticles<=0))
    {
      cout << "Need to fill at least one shell with particles"<<endl;
    }
  if (nbrParticlesUp<0)
    {
      if (nbrParticles&1)
	{
	  cout << "Please give explicit number of up-particles for any layer unbalanced configuration!"<<endl;
	  exit(1);
	}
      this->NbrParticlesUp = nbrParticles/2;
      this->NbrParticlesDown = nbrParticles/2;
    }
  else
    {
      this->NbrParticlesUp = nbrParticlesUp;
      this->NbrParticlesDown = nbrParticles-nbrParticlesUp;
    }
  this->LFBonding=lFBonding;
  this->LFAntiBonding=lFAntiBonding;
  this->MaxLF = (lFBonding>lFAntiBonding?lFBonding:lFAntiBonding);
  this->AbsEffectiveFlux = abs(nbrEffectiveFlux);
  this->NbrOrbsBonding = LFBonding*(AbsEffectiveFlux+1)+LFBonding*(LFBonding-1);
  this->NbrOrbsAntiBonding = LFAntiBonding*(AbsEffectiveFlux+1)+LFAntiBonding*(LFAntiBonding-1);
  this->NbrParticles = this->NbrOrbsBonding + this->NbrOrbsAntiBonding;  
  if (this->NbrParticles!=nbrParticles)
    {
      cout << "Total number of particles does not agree with parameters!"<<endl;
      exit(-1);
    }
  this->Orbitals1 = new JainCFOnSphereOrbitals(NbrParticlesUp, MaxLF, nbrEffectiveFlux, jastrowPower);
  this->Orbitals2 = new JainCFOnSphereOrbitals(NbrParticlesDown, MaxLF, nbrEffectiveFlux, jastrowPower);
  this->ElementNorm=1.0;
  this->Flag.Initialize();
  
#ifdef USE_LAPACK_CFCB
  this->Matrix = new ComplexLapackDeterminant(this->NbrParticles);
#else
  this->Matrix = new ComplexMatrix(this->NbrParticles, this->NbrParticles);
#endif
  
  this->J11 = new Complex[this->NbrParticlesUp];
//   this->J12 = new Complex[this->NbrParticlesUp];
//   this->J21 = new Complex[this->NbrParticlesDown];
  this->J22 = new Complex[this->NbrParticlesDown];
  
//   this->InterLayerDistances = new Complex*[this->NbrParticlesUp];
//   for (int i=0; i< this->NbrParticlesPerLayer; ++i)
//     InterLayerDistances[i]= new Complex[NbrParticlesDown];
  
}

// copy constructor
//
// function = reference on the wave function to copy

CFOnSphereWithSpinPartonTunnellingWaveFunction::CFOnSphereWithSpinPartonTunnellingWaveFunction(const CFOnSphereWithSpinPartonTunnellingWaveFunction& function)
{
  this->NbrParticlesUp = function.NbrParticlesUp;
  this->NbrParticlesDown = function.NbrParticlesDown;
  this->NbrOrbsBonding = function.NbrOrbsBonding;
  this->NbrOrbsAntiBonding = function.NbrOrbsAntiBonding;
  this->NbrParticles = function.NbrParticles;
  this->LFBonding = function.LFBonding;
  this->LFAntiBonding = function.LFAntiBonding;
  this->MaxLF = function.MaxLF;
  this->AbsEffectiveFlux = function.AbsEffectiveFlux;
  this->Flag = function.Flag;  
  this->Orbitals1 = function.Orbitals1;
  this->Orbitals2 = function.Orbitals2;
  this->ElementNorm=function.ElementNorm;

#ifdef USE_LAPACK_CFCB
  this->Matrix = new ComplexLapackDeterminant(this->NbrParticles);
#else
  this->Matrix = new ComplexMatrix(this->NbrParticles, this->NbrParticles);
#endif
  
  this->J11 = new Complex[this->NbrParticlesUp];
//   this->J12 = new Complex[this->NbrParticlesUp];
//   this->J21 = new Complex[this->NbrParticlesDown];
  this->J22 = new Complex[this->NbrParticlesDown];
  
//   this->InterLayerDistances = new Complex*[this->NbrParticlesUp];
//   for (int i=0; i< this->NbrParticlesPerLayer; ++i)
//     InterLayerDistances[i]= new Complex[NbrParticlesDown];
}

// destructor
//

CFOnSphereWithSpinPartonTunnellingWaveFunction::~CFOnSphereWithSpinPartonTunnellingWaveFunction()
{
  if ( (this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete Orbitals1;
      delete Orbitals2;
    }
//   for (int i=0; i< this->NbrParticlesUp; ++i)
//     delete [] InterLayerDistances[i];
//   delete [] InterLayerDistances;
  delete [] J11;
//   delete [] J12;
//   delete [] J21;
  delete [] J22;
  delete Matrix;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* CFOnSphereWithSpinPartonTunnellingWaveFunction::Clone ()
{
  return new CFOnSphereWithSpinPartonTunnellingWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex CFOnSphereWithSpinPartonTunnellingWaveFunction::operator ()(RealVector& x)
{
  RealVector part=x.Extract(0, 2*this->NbrParticlesUp-1);
  this->OrbitalValues1 = (*Orbitals1)(part);
  part = x.Extract(2*this->NbrParticlesUp, 2*(this->NbrParticles)-1);
  this->OrbitalValues2 = (*Orbitals2)(part);
  this->EvaluateTables();
  
  return Matrix->Determinant()*Interpolation;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex CFOnSphereWithSpinPartonTunnellingWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  ComplexVector part=uv.Extract(0, 2*this->NbrParticlesUp-1);
  this->OrbitalValues1 = Orbitals1->CalculateFromSpinorVariables(part);
  part = uv.Extract(2*this->NbrParticlesUp, 2*(this->NbrParticles)-1);
  this->OrbitalValues2 = Orbitals2->CalculateFromSpinorVariables(part);
  
  this->EvaluateTables();

  return Matrix->Determinant()*Interpolation;
}


// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void CFOnSphereWithSpinPartonTunnellingWaveFunction::AdaptNorm(RealVector& x)
{
  double det=Norm((*this)(x));
  while ((det<.1)||(det>50.0))
    {
      cout <<"N'="<< this->ElementNorm << " det="<<det<<endl;
      if (det>1e300) 
	this->ElementNorm*= pow((double)1.0e-300,(double)1.0/(this->NbrParticles));
      else if (det==0.0) 
	this->ElementNorm*= pow((double)1.0e300,(double)1.0/(this->NbrParticles));
      else 
	this->ElementNorm*= pow(det,(double)-1.0/(this->NbrParticles));
      det=Norm((*this)(x));
      cout <<"N'="<< this->ElementNorm << endl;
    }
}


// normalize the wave-function over an average number of MC positions

void CFOnSphereWithSpinPartonTunnellingWaveFunction::AdaptAverageMCNorm(int thermalize, int average)
{
  ParticleOnSphereCollection * Particles = new ParticleOnSphereCollection(2*this->NbrParticles);
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
  this->ElementNorm*= pow(SumTrialValues/average,(double)-2.0/(this->NbrParticles));
  delete Particles;
}


// this is the main part of the calculation of the paired wavefunction:
// assumes that Orbitals1 and Orbitals2 have already been initialized prior to call
void CFOnSphereWithSpinPartonTunnellingWaveFunction::EvaluateTables()
{
  int alpha;  
  Complex TmpC;
  // evaluate single particle Jastrow factors
  this->Interpolation=1.0;
  if (Orbitals1->TestCriticality(Interpolation) == 0)
    {
      for (int i=0;i<this->NbrParticlesUp;i++)
	{
	  J11[i]=this->ElementNorm;
	  for(int j=0;j<i;j++) J11[i] *= Orbitals1->JastrowFactorElement(i,j);
	  for(int j=i+1;j<NbrParticlesUp;j++) J11[i] *= Orbitals1->JastrowFactorElement(i,j);
	}
    }
  else // if some interpolation occurred, the true values of the J11's have to be recalculated:
    {
      for (int i=0;i<this->NbrParticlesUp;i++)
	{
	  J11[i]=this->ElementNorm;
	  for(int j=0;j<i;j++) J11[i] *= ((Orbitals1->SpinorU(i) * Orbitals1->SpinorV(j)) - (Orbitals1->SpinorU(j) * Orbitals1->SpinorV(i)));
	  for(int j=i+1;j<NbrParticlesUp;j++) J11[i] *= ((Orbitals1->SpinorU(i) * Orbitals1->SpinorV(j)) - (Orbitals1->SpinorU(j) * Orbitals1->SpinorV(i)));
	}
    }

  if (Orbitals2->TestCriticality(Interpolation) == 0)
    {
      for (int i=0;i<this->NbrParticlesDown;i++)
	{
	  J22[i]=this->ElementNorm;
	  for(int j=0;j<i;j++) J22[i] *= Orbitals2->JastrowFactorElement(i,j);
	  for(int j=i+1;j<NbrParticlesDown;j++) J22[i] *= Orbitals2->JastrowFactorElement(i,j);
	}
    }
  else // if some interpolation occurred, the true values of the J22's have to be recalculated:
    {
      for (int i=0;i<this->NbrParticlesDown;i++)
	{
	  J22[i]=this->ElementNorm;
	  for(int j=0;j<i;j++) J22[i] *= ((Orbitals2->SpinorU(i) * Orbitals2->SpinorV(j)) - (Orbitals2->SpinorU(j) * Orbitals2->SpinorV(i)));
	  for(int j=i+1;j<NbrParticlesDown;j++) J22[i] *= ((Orbitals2->SpinorU(i) * Orbitals2->SpinorV(j)) - (Orbitals2->SpinorU(j) * Orbitals2->SpinorV(i)));
	}
    }

  // initialize Slater determinant 

  // up particles: first Nup lines
  for (int i=0;i<this->NbrParticlesUp;++i)
    {
      alpha=0;
      for (int n=0;n<this->LFBonding;n++)
	for (int m2=-AbsEffectiveFlux-2*n; m2<=AbsEffectiveFlux+2*n;m2+=2)
	  {
	    TmpC=OrbitalValues1[alpha][i]*J11[i];
#ifdef USE_LAPACK_CFCB
	    Matrix->SetMatrixElement(i, alpha, Real(TmpC), Imag(TmpC));
#else
	    (*Matrix)[i].Re(alpha) = Real(TmpC);
	    (*Matrix)[i].Im(alpha) = Imag(TmpC);
#endif
	    ++alpha;
	  }
      alpha=0;
      for (int n=0;n<this->LFAntiBonding;n++)
	for (int m2=-AbsEffectiveFlux-2*n; m2<=AbsEffectiveFlux+2*n;m2+=2)
	  {
	    TmpC=OrbitalValues1[alpha][i]*J11[i];
#ifdef USE_LAPACK_CFCB
	    Matrix->SetMatrixElement(i, this->NbrOrbsBonding+alpha, Real(TmpC), Imag(TmpC));
#else
	    (*Matrix)[i].Re(this->NbrOrbsBonding+alpha) = Real(TmpC);
	    (*Matrix)[i].Im(this->NbrOrbsBonding+alpha) = Imag(TmpC);
#endif
	    ++alpha;
	  }
    }
  
  // down particles - lines Nup+1...N
  for (int i=0;i<this->NbrParticlesDown;++i)
    {
      alpha=0;
      for (int n=0;n<this->LFBonding;n++)
	for (int m2=-AbsEffectiveFlux-2*n; m2<=AbsEffectiveFlux+2*n;m2+=2)
	  {
	    TmpC=OrbitalValues2[alpha][i]*J22[i];
#ifdef USE_LAPACK_CFCB
	    Matrix->SetMatrixElement(this->NbrParticlesUp+i, alpha, Real(TmpC), Imag(TmpC));
#else
	    (*Matrix)[this->NbrParticlesUp+i].Re(alpha) = Real(TmpC);
	    (*Matrix)[this->NbrParticlesUp+i].Im(alpha) = Imag(TmpC);
#endif
	    ++alpha;
	  }
      alpha=0;
      for (int n=0;n<this->LFAntiBonding;n++)
	for (int m2=-AbsEffectiveFlux-2*n; m2<=AbsEffectiveFlux+2*n;m2+=2)
	  {
	    // minus sign here for asymmetric states!
	    TmpC=-OrbitalValues2[alpha][i]*J22[i];
#ifdef USE_LAPACK_CFCB
	    Matrix->SetMatrixElement(this->NbrParticlesUp+i, this->NbrOrbsBonding+alpha, Real(TmpC), Imag(TmpC));
#else
	    (*Matrix)[this->NbrParticlesUp+i].Re(this->NbrOrbsBonding+alpha) = Real(TmpC);
	    (*Matrix)[this->NbrParticlesUp+i].Im(this->NbrOrbsBonding+alpha) = Imag(TmpC);
#endif
	    ++alpha;
	  }
    }
  
//   for (int i=0; i<NbrParticles; ++i)
//     for (int j=0; j<NbrParticles; ++j)
//       {
// 	Matrix->GetMatrixElement(i,j,TmpC);
// 	cout << "M["<<i<<","<<j<<"]="<<TmpC<<endl;
//       }
}
