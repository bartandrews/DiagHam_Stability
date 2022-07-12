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
#include "TwoThirdSingletState.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
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

TwoThirdSingletState::TwoThirdSingletState()
{
  this->NbrParticles = 0;
}

// constructor
//
// nbrParticles = number of particles
//
TwoThirdSingletState::TwoThirdSingletState(int nbrParticles)
{
  if (nbrParticles & 1)
    {
      cout << "This spin-singlet state requires an even number of fermions" << endl;
      exit(1);
    }
  this->NbrParticles = nbrParticles;
  this->EffectiveFlux = -NbrParticles/2 + 2;
  this->OrbitalFactory = new JainCFOnSphereOrbitals(NbrParticles, /*NbrLandauLevels*/ 2, EffectiveFlux, /*JastrowPower*/ 2);
  this->SlaterNorm=1.0;
  this->CauchyNorm=1.0;
  
  this->Flag.Initialize();
  this->Interpolation=1.0;
  this->Slater = new ComplexMatrix(this->NbrParticles,this->NbrParticles);
  this->CauchyPermanent = new ComplexMatrix(this->NbrParticles/2,this->NbrParticles/2);
  this->Jij = new Complex*[this->NbrParticles];
  for (int i=0; i<NbrParticles; ++i)
    this->Jij[i] = new Complex[this->NbrParticles];

  
}

// copy constructor
//
// function = reference on the wave function to copy

TwoThirdSingletState::TwoThirdSingletState(const TwoThirdSingletState& function)
{
  this->NbrParticles = function.NbrParticles;
  this->EffectiveFlux = function.EffectiveFlux;
  this->OrbitalFactory = function.OrbitalFactory;
  this->SlaterNorm=function.SlaterNorm;
  this->CauchyNorm=function.CauchyNorm;
  this->Flag = function.Flag;
  this->Interpolation=function.Interpolation;
  this->Slater = new ComplexMatrix(this->NbrParticles,this->NbrParticles);
  this->CauchyPermanent = new ComplexMatrix(this->NbrParticles/2,this->NbrParticles/2);
  this->Jij = new Complex*[this->NbrParticles];
  for (int i=0; i<NbrParticles; ++i)
    this->Jij[i] = new Complex[this->NbrParticles];
}

// destructor
//

TwoThirdSingletState::~TwoThirdSingletState()
{
  if (this->NbrParticles>0)
    {
      if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
	delete OrbitalFactory;
      delete Slater;
      delete CauchyPermanent;
      for (int i=0; i<NbrParticles;++i)
	delete [] Jij[i];
      delete [] Jij;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* TwoThirdSingletState::Clone ()
{
  return new TwoThirdSingletState(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex TwoThirdSingletState::operator ()(RealVector& x)
{
  this->Orbitals = (*OrbitalFactory)(x);
  this->EvaluateTables();
  Complex Ji;
  Complex Tmp;
  // initialize Slater determinant (or Pfaffian matrix)
  for (int i=0;i<this->NbrParticles;++i)
    {
      Ji=1.0;
      for(int j=0;j<i;++j) Ji*=Jij[i][j];
      for(int j=i+1;j<NbrParticles;++j) Ji*=Jij[i][j];
      Ji*=SlaterNorm;
      for(int j=0;j<NbrParticles;++j)
	{
	  Orbitals.GetMatrixElement(i,j,Tmp);
	  Tmp*=Ji;	  
	  Slater->SetMatrixElement(i,j, Tmp);
	}
    }
  this->DeterminantValue=Slater->Determinant();
  this->PermanentValue=CauchyPermanent->Permanent();
  return DeterminantValue*Interpolation*PermanentValue;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex TwoThirdSingletState::CalculateFromSpinorVariables(ComplexVector& uv)
{
  this->Orbitals = OrbitalFactory->CalculateFromSpinorVariables(uv);
  this->EvaluateTables();
  Complex Ji;
  Complex Tmp;
  // initialize Slater determinant (or Pfaffian matrix)
  for (int i=0;i<this->NbrParticles;++i)
    {
      Ji=1.0;
      for(int j=0;j<i;++j) Ji*=Jij[i][j];
      for(int j=i+1;j<NbrParticles;++j) Ji*=Jij[i][j];
      Ji*=SlaterNorm;
      for(int j=0;j<NbrParticles;++j)
	{
	  Orbitals.GetMatrixElement(i,j,Tmp);
	  Tmp*=Ji;	  
	  Slater->SetMatrixElement(i,j, Tmp);
	}
    }
  this->DeterminantValue=Slater->Determinant();
  this->PermanentValue=CauchyPermanent->Permanent();
  return DeterminantValue*Interpolation*PermanentValue;
}

Complex TwoThirdSingletState::GetTestValue(RealVector& x)
{
  this->Orbitals = (*OrbitalFactory)(x);
  this->EvaluateTables();
  Complex Ji;
  Complex Tmp;
  // initialize Slater determinant (or Pfaffian matrix)
  for (int i=0;i<this->NbrParticles;++i)
    {
      Ji=1.0;
      for(int j=0;j<i;++j) Ji*=Jij[i][j];
      for(int j=i+1;j<NbrParticles;++j) Ji*=Jij[i][j];
      Ji*=SlaterNorm;
      for(int j=0;j<NbrParticles;++j)
	{
	  Orbitals.GetMatrixElement(i,j,Tmp);
	  Tmp*=Ji;	  
	  Slater->SetMatrixElement(i,j, Tmp);
	}
    }
  cout << "Cauchy-Permanent:" << *CauchyPermanent << endl;
  this->DeterminantValue=Slater->Determinant();
  this->PermanentValue=CauchyPermanent->Permanent();
  cout << "DeterminantValue=" <<DeterminantValue<<endl;
  cout << "PermanentValue=" <<PermanentValue<<endl;
  return DeterminantValue*Interpolation*PermanentValue;
}



// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void TwoThirdSingletState::AdaptNorm(RealVector& x)
{
  double TotalNorm=Norm((*this)(x));
  double DeterminantNorm=Norm(DeterminantValue);  
  while ((DeterminantNorm<.1)||(DeterminantNorm>50.0))
    {
      cout <<"N="<< this->SlaterNorm << " DeterminantNorm="<<DeterminantNorm<<endl;
      if (DeterminantNorm>1e300) 
	this->SlaterNorm*= pow((double)1.0e-300,(double)1.0/this->NbrParticles);
      else if (DeterminantNorm==0.0) 
	this->SlaterNorm*= pow((double)1.0e300,(double)1.0/this->NbrParticles);
      else 
	this->SlaterNorm*= pow(DeterminantNorm,(double)-1.0/this->NbrParticles);
      TotalNorm=Norm((*this)(x));
      DeterminantNorm=Norm(DeterminantValue);
      cout <<"N'="<< this->SlaterNorm << " DeterminantNorm="<<DeterminantNorm<<endl;
    }
  double PermanentNorm=Norm(PermanentValue);
  while ((PermanentNorm<.1)||(PermanentNorm>50.0))
    {      
      cout <<"N="<< this->CauchyNorm << " PermanentNorm="<<PermanentNorm<<endl;
      if (PermanentNorm>1e300) 
	this->CauchyNorm*= pow((double)1.0e-300,(double)2.0/this->NbrParticles);
      else if (PermanentNorm==0.0) 
	this->CauchyNorm*= pow((double)1.0e300,(double)2.0/this->NbrParticles);
      else 
	this->CauchyNorm*= pow(PermanentNorm,(double)-2.0/this->NbrParticles);
      TotalNorm=Norm((*this)(x));
      PermanentNorm=Norm(PermanentValue);
      cout <<"N'="<< this->CauchyNorm << " PermanentNorm="<<PermanentNorm<<endl;
    }
  if ((TotalNorm<.1)||(TotalNorm>50.0))
    {
      cout << "Total Norm should also be well panned in, now... check code." << endl;
      exit(1);
    }
  
}


// normalize the wave-function over an average number of MC positions

void TwoThirdSingletState::AdaptAverageMCNorm(int thermalize, int average)
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
  this->SlaterNorm*= pow(SumTrialValues/average,(double)-1.0/this->NbrParticles);
  delete Particles;
}


// this is the main part of the calculation of the paired wavefunction:
// assumes Orbitals initialized
//
void TwoThirdSingletState::EvaluateTables()
{
  int i, j;
  Complex tmp;
  // evaluate single particle Jastrow factors
  this->Interpolation=1.0;
  if (OrbitalFactory->TestCriticality(Interpolation) == 0)
    {
      for (i=0;i<this->NbrParticles;i++)
	{
	  for(j=0;j<i;j++) Jij[i][j] = OrbitalFactory->JastrowFactorElement(i,j);
	  for(j=i+1;j<NbrParticles;j++) Jij[i][j] = OrbitalFactory->JastrowFactorElement(i,j);
	}
    }
  else // if some interpolation occurred, the true values of the Ji's have to be recalculated:
    {
      cout << "Critical Event" << endl;
      for (i=0;i<this->NbrParticles;i++)
       	{	  
       	  for(j=0;j<i;j++) Jij[i][j] = ((OrbitalFactory->SpinorU(i) * OrbitalFactory->SpinorV(j)) - (OrbitalFactory->SpinorU(j) * OrbitalFactory->SpinorV(i)));
       	  for(j=i+1;j<NbrParticles;j++) Jij[i][j] = ((OrbitalFactory->SpinorU(i) * OrbitalFactory->SpinorV(j)) - (OrbitalFactory->SpinorU(j) * OrbitalFactory->SpinorV(i)));
      	}
    }
  Complex Power,Tmp;
  int N1=NbrParticles/2;
  for (int i=0; i<N1;++i)
    for (int j=0; j<N1;++j)
      {
	Tmp = (Power = 1.0/Jij[i][j+N1]);
	//	Power *= Tmp;
	//	Power *= Tmp;
	Power *= this->CauchyNorm;
	CauchyPermanent->SetMatrixElement(i,j, Power);
      }
}


      
