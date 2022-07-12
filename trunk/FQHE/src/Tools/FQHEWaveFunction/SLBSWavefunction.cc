////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                          //
//                                                                            //
//                                                                            //
//           class implementing generalized Halperin states on the sphere     //
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
#include "SLBSWavefunction.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <fstream>

using std::ios;
using std::ofstream;
using std::cout;
using std::endl;


// default constructor
//

SLBSWavefunction::SLBSWavefunction()
{
  this->NbrParticles = 0;
}

// constructor
//
// nbrParticles = number of particles
// negFluxFlag = indicating sign of Jain wavefunction component
//
SLBSWavefunction::SLBSWavefunction(int nbrParticles, bool negFluxFlag)
{
  if ((nbrParticles & 1)||(nbrParticles < 4))
    {
      cout << "This state requires an even number of fermions >= 4" << endl;
      exit(1);
    }
  this->NbrParticles = nbrParticles;
  // discard "bad" state with positive flux:
  negFluxFlag = true;
  this->NegativeFieldFlag = negFluxFlag;
  this->EffectiveFlux = ((NbrParticles -2)/2 - 1)*(1-2*NegativeFieldFlag);
  cout << "EffectiveFlux="<<EffectiveFlux<<endl;
  this->OrbitalFactory = new JainCFOnSphereOrbitals(NbrParticles, /*NbrLandauLevels*/ 2, EffectiveFlux, /*JastrowPower*/ 2);
  
  this->Flag.Initialize();
  this->Interpolation=1.0;
  this->Slater = new ComplexMatrix(this->NbrParticles,this->NbrParticles);
  this->Pfaffian = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->SlaterElementNorm = 1.0;

  this->Jij = new Complex*[this->NbrParticles];
  for (int i=0; i<NbrParticles; ++i)
    this->Jij[i] = new Complex[this->NbrParticles];

  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];
}


// constructor
//
// nbrParticles = number of particles
// levels = number of LL's in CF part
//
SLBSWavefunction::SLBSWavefunction(int nbrParticles, int levels)
{
  int Levels = abs(levels);
  if (((nbrParticles % Levels)!=0)||(nbrParticles < 2*Levels))
    {
      cout << "This state requires an multiple of "<<Levels<<" fermions >= "<<2*Levels << endl;
      exit(1);
    }
  this->NbrParticles = nbrParticles;
  if (levels<0)
    this->NegativeFieldFlag = true;
  else
    this->NegativeFieldFlag = false;
  this->EffectiveFlux = ((NbrParticles -Levels*(Levels-1))/Levels - 1)*(1-2*NegativeFieldFlag);
  cout << "EffectiveFlux="<<EffectiveFlux<<endl;
  this->OrbitalFactory = new JainCFOnSphereOrbitals(NbrParticles, /*NbrLandauLevels*/ Levels, EffectiveFlux, /*JastrowPower*/ 2);
  
  this->Flag.Initialize();
  this->Interpolation=1.0;
  this->Slater = new ComplexMatrix(this->NbrParticles,this->NbrParticles);
  this->Pfaffian = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->SlaterElementNorm = 1.0;

  this->Jij = new Complex*[this->NbrParticles];
  for (int i=0; i<NbrParticles; ++i)
    this->Jij[i] = new Complex[this->NbrParticles];

  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];
}

// copy constructor
//
// function = reference on the wave function to copy

SLBSWavefunction::SLBSWavefunction(const SLBSWavefunction& function)
{
  // to adapt
  this->NbrParticles = function.NbrParticles;
  this->EffectiveFlux = function.EffectiveFlux;
  this->OrbitalFactory = function.OrbitalFactory;
  this->SlaterElementNorm=function.SlaterElementNorm;
  this->Determinant = function.Determinant;
  this->NegativeFieldFlag = function.NegativeFieldFlag;
  this->EffectiveFlux = function.EffectiveFlux;
  this->Flag = function.Flag;
  this->Interpolation=function.Interpolation;
  this->Slater = new ComplexMatrix(this->NbrParticles,this->NbrParticles);
  this->Pfaffian = new ComplexSkewSymmetricMatrix(this->NbrParticles);
    
  this->Jij = new Complex*[this->NbrParticles];
  for (int i=0; i<NbrParticles; ++i)
    this->Jij[i] = new Complex[this->NbrParticles];

  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];
}

// destructor
//

SLBSWavefunction::~SLBSWavefunction()
{
  if (this->NbrParticles>0)
    {
      if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
	{
	  delete OrbitalFactory;
	}
      delete Slater;
      delete Pfaffian;
      delete [] SpinorUCoordinates;
      delete [] SpinorVCoordinates;
      for (int i=0; i<NbrParticles;++i)
	delete [] Jij[i];
      delete [] Jij;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* SLBSWavefunction::Clone ()
{
  return new SLBSWavefunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex SLBSWavefunction::operator ()(RealVector& x)
{
  double c, s;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      SpinorUCoordinates[i].Re = cos(0.5 * x[i << 1]);
      SpinorUCoordinates[i].Im = SpinorUCoordinates[i].Re;
      SpinorUCoordinates[i].Re *= (c=cos(0.5 * x[1 + (i << 1)]));
      SpinorUCoordinates[i].Im *= -(s=sin(0.5 * x[1 + (i << 1)]));
      SpinorVCoordinates[i].Re = sin(0.5 * x[i << 1]);
      SpinorVCoordinates[i].Im = SpinorVCoordinates[i].Re;
      SpinorVCoordinates[i].Re *= c;
      SpinorVCoordinates[i].Im *= s;
    }
  this->Orbitals = (*OrbitalFactory)(x);
  
  this->EvaluateTables();
  Complex Ji;
  Complex Tmp;

  // initialize Slater determinant and Pfaffian
  for (int i=0;i<this->NbrParticles;++i)
    {
      Ji=1.0;
      for(int j=0;j<i;++j)
	{
	  Ji*=Jij[i][j];
	  Pfaffian->SetMatrixElement(i,j,1.0/Jij[i][j]);
	}
      for(int j=i+1;j<NbrParticles;++j) Ji*=Jij[i][j];
      Ji*=SlaterElementNorm;
      for(int j=0;j<NbrParticles;++j)
	{
	  Orbitals.GetMatrixElement(i,j, Tmp);
	  // cout << i<<", "<<j<<": "<<Tmp<<endl;		  
	  Tmp*=Ji;	  
	  Slater->SetMatrixElement(i,j, Tmp);
	}
      // cout << endl;
    }
  this->Determinant=Slater->Determinant();
  this->PfaffianValue=Pfaffian->Pfaffian();

  return this->PfaffianValue*Interpolation*this->Determinant*this->Chi1;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex SLBSWavefunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = uv.Re(2*i);
      this->SpinorUCoordinates[i].Im = uv.Im(2*i);
      this->SpinorVCoordinates[i].Re = uv.Re(2*i+1);
      this->SpinorVCoordinates[i].Im = uv.Im(2*i+1);
    }
  this->Orbitals = OrbitalFactory->CalculateFromSpinorVariables(uv);
  this->EvaluateTables();

    Complex Ji;
  Complex Tmp;

  // initialize Slater determinant and Pfaffian
  for (int i=0;i<this->NbrParticles;++i)
    {
      Ji=1.0;
      for(int j=0;j<i;++j)
	{
	  Ji*=Jij[i][j];
	  Pfaffian->SetMatrixElement(i,j,1.0/Jij[i][j]);
	}
      for(int j=i+1;j<NbrParticles;++j) Ji*=Jij[i][j];
      Ji*=SlaterElementNorm;
      for(int j=0;j<NbrParticles;++j)
	{
	  Orbitals.GetMatrixElement(i,j, Tmp);
	  // cout << i<<", "<<j<<": "<<Tmp<<endl;		  
	  Tmp*=Ji;	  
	  Slater->SetMatrixElement(i,j, Tmp);
	}
      // cout << endl;
    }
  this->Determinant=Slater->Determinant();
  this->PfaffianValue=Pfaffian->Pfaffian();

  return this->PfaffianValue*Interpolation*this->Determinant*this->Chi1;
}



// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void SLBSWavefunction::AdaptNorm(RealVector& x)
{
  double TotalNorm=Norm((*this)(x));
  while ((TotalNorm<.1)||(TotalNorm>50.0))
    {
      cout <<"N="<< this->SlaterElementNorm << " Norm(Psi)="<<TotalNorm<<endl;
      if (TotalNorm>1e300) 
	this->SlaterElementNorm*= pow((double)1.0e-300,(double)1.0/this->NbrParticles);
      else if (TotalNorm==0.0) 
	this->SlaterElementNorm*= pow((double)1.0e300,(double)1.0/this->NbrParticles);
      else 
	this->SlaterElementNorm*= pow(TotalNorm,(double)-1.0/this->NbrParticles);
      TotalNorm=Norm((*this)(x));
      cout <<"N'="<< this->SlaterElementNorm << " Norm(Psi)="<<TotalNorm<<endl;
    }
}


// normalize the wave-function over an average number of MC positions

void SLBSWavefunction::AdaptAverageMCNorm(int thermalize, int average)
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
  this->SlaterElementNorm*= pow(SumTrialValues/average,(double)-1.0/this->NbrParticles);
  delete Particles;
}


// this is the main part of the calculation of the paired wavefunction:
// assumes Orbitals initialized
//
void SLBSWavefunction::EvaluateTables()
{
  int i, j;
  Complex tmp;
  // evaluate single particle Jastrow factors
  this->Interpolation=1.0;
  if (OrbitalFactory->TestCriticality(Interpolation) == 0)
    {
      for (i=0;i<this->NbrParticles;i++)
	{
	  for(j=0;j<i;j++)
	    {
	      Jij[i][j] = OrbitalFactory->JastrowFactorElement(i,j);
	      Jij[j][i] = - Jij[i][j];
	    }
	  Jij[i][j] = 0.0;
	  for(j=i+1;j<NbrParticles;j++)
	    {
	      Jij[i][j] = OrbitalFactory->JastrowFactorElement(i,j);
	      Jij[j][i] = - Jij[i][j];
	    }
	}
    }
  else // if some interpolation occurred, the true values of the Ji's have to be recalculated:
    {
      cout << "Critical Event" << endl;
      for (i=0;i<this->NbrParticles;i++)
       	{	  
       	  for(j=0;j<i;j++)
	    {
	      Jij[i][j] = ((OrbitalFactory->SpinorU(i) * OrbitalFactory->SpinorV(j)) - (OrbitalFactory->SpinorU(j) * OrbitalFactory->SpinorV(i)));
	      Jij[j][i] = - Jij[i][j];
	    }
       	  for(j=i+1;j<NbrParticles;j++)
	    {
	      Jij[i][j] = ((OrbitalFactory->SpinorU(i) * OrbitalFactory->SpinorV(j)) - (OrbitalFactory->SpinorU(j) * OrbitalFactory->SpinorV(i)));
	      Jij[j][i] = - Jij[i][j];
	    }
      	}
    }

  double Factor = M_PI * 0.5;
  this->Chi1=1.0;
  for (int i = 1; i < this->NbrParticles; ++i)
    for (int j = 0; j < i; ++j)
      Chi1*=(Jij[i][j]*Factor);

  if (!NegativeFieldFlag) // divide by Chi1 if positive flux, multiply if negative!
    {
      Chi1=1.0/Chi1;
    }
}


      
