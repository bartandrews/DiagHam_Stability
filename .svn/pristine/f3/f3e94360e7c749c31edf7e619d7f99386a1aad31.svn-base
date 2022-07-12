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
#include "TwoThirdUnpolarizedCF.h"
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

TwoThirdUnpolarizedCF::TwoThirdUnpolarizedCF()
{
  this->NbrParticles = 0;
}

// constructor
//
// nbrParticles = number of particles
//
TwoThirdUnpolarizedCF::TwoThirdUnpolarizedCF(int nbrParticles)
{
  if (nbrParticles & 1)
    {
      cout << "This spin-singlet state requires an even number of fermions" << endl;
      exit(1);
    }
  this->NbrParticles = nbrParticles;
  this->NbrParticlesPerSpin = nbrParticles/2;
  this->EffectiveFlux = -NbrParticlesPerSpin + 1;
  this->OrbitalFactoryUp = new JainCFOnSphereOrbitals(NbrParticlesPerSpin, /*NbrLandauLevels*/ 1, EffectiveFlux, /*JastrowPower*/ 2);
  this->OrbitalFactoryDown = new JainCFOnSphereOrbitals(NbrParticlesPerSpin, /*NbrLandauLevels*/ 1, EffectiveFlux, /*JastrowPower*/ 2);
  this->SlaterElementNorm=1.0;
  
  this->Flag.Initialize();
  this->InterpolationUp=1.0;
  this->InterpolationDown=1.0;
  this->SlaterUp = new ComplexMatrix(this->NbrParticlesPerSpin,this->NbrParticlesPerSpin);
  this->SlaterDown = new ComplexMatrix(this->NbrParticlesPerSpin,this->NbrParticlesPerSpin);

  this->Jij = new Complex*[this->NbrParticles];
  for (int i=0; i<NbrParticles; ++i)
    this->Jij[i] = new Complex[this->NbrParticles];

  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];

  
}

// copy constructor
//
// function = reference on the wave function to copy

TwoThirdUnpolarizedCF::TwoThirdUnpolarizedCF(const TwoThirdUnpolarizedCF& function)
{
  // to adapt
  this->NbrParticles = function.NbrParticles;
  this->EffectiveFlux = function.EffectiveFlux;
  this->OrbitalFactoryUp = function.OrbitalFactoryUp;
  this->OrbitalFactoryDown = function.OrbitalFactoryDown;
  this->SlaterElementNorm=function.SlaterElementNorm;
  this->DeterminantUpValue = function.DeterminantUpValue;
  this->DeterminantDownValue = function.DeterminantDownValue;
  this->InterSpinJastrow=function.InterSpinJastrow;
  this->Flag = function.Flag;
  this->InterpolationUp=function.InterpolationUp;
  this->InterpolationDown=function.InterpolationDown;
  this->SlaterUp = new ComplexMatrix(this->NbrParticlesPerSpin,this->NbrParticlesPerSpin);
  this->SlaterDown = new ComplexMatrix(this->NbrParticlesPerSpin,this->NbrParticlesPerSpin);
  this->Jij = new Complex*[this->NbrParticles];
  for (int i=0; i<NbrParticles; ++i)
    this->Jij[i] = new Complex[this->NbrParticles];

  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];
}

// destructor
//

TwoThirdUnpolarizedCF::~TwoThirdUnpolarizedCF()
{
  if (this->NbrParticles>0)
    {
      if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
	{
	  delete OrbitalFactoryUp;
	  delete OrbitalFactoryDown;
	}
      delete SlaterUp;
      delete SlaterDown;
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

Abstract1DComplexFunction* TwoThirdUnpolarizedCF::Clone ()
{
  return new TwoThirdUnpolarizedCF(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex TwoThirdUnpolarizedCF::operator ()(RealVector& x)
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
      //cout << "U["<<i<<"]="<<SpinorUCoordinates[i]<<", "<< "V["<<i<<"]="<<SpinorVCoordinates[i]<<endl;
    }
  
  RealVector part=x.Extract(0, 2*this->NbrParticlesPerSpin-1);
  this->OrbitalsUp = (*OrbitalFactoryUp)(part);
  part = x.Extract(2*this->NbrParticlesPerSpin, 4*this->NbrParticlesPerSpin-1);
  this->OrbitalsDown = (*OrbitalFactoryDown)(part);  
  
  this->EvaluateTables();
  Complex Ji;
  Complex Tmp;
  // initialize Slater determinants
  for (int i=0;i<this->NbrParticlesPerSpin;++i)
    {
      Ji=1.0;
      for(int j=0;j<i;++j) Ji*=Jij[i][j];
      for(int j=i+1;j<NbrParticlesPerSpin;++j) Ji*=Jij[i][j];
      Ji*=SlaterElementNorm;
      for(int j=0;j<NbrParticlesPerSpin;++j)
	{
	  OrbitalsUp.GetMatrixElement(i,j,Tmp);
	  Tmp*=Ji;	  
	  SlaterUp->SetMatrixElement(i,j, Tmp);
	}
    }
  this->DeterminantUpValue=SlaterUp->Determinant();

  // initialize Slater determinants
  for (int i=this->NbrParticlesPerSpin;i<2*this->NbrParticlesPerSpin;++i)
    {
      Ji=1.0;
      for(int j=this->NbrParticlesPerSpin;j<i;++j) Ji*=Jij[i][j];
      for(int j=i+1;j<NbrParticles;++j) Ji*=Jij[i][j];
      Ji*=SlaterElementNorm;
      for(int j=0;j<this->NbrParticlesPerSpin;++j)
	{
	  OrbitalsDown.GetMatrixElement(i-this->NbrParticlesPerSpin,j,Tmp);
	  Tmp*=Ji;	  
	  SlaterDown->SetMatrixElement(i-this->NbrParticlesPerSpin,j, Tmp);
	}
    }
  this->DeterminantDownValue=SlaterDown->Determinant();
  
  return DeterminantUpValue*InterpolationUp*DeterminantDownValue*InterpolationDown*InterSpinJastrow;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex TwoThirdUnpolarizedCF::CalculateFromSpinorVariables(ComplexVector& uv)
{

  ComplexVector part=uv.Extract(0, 2*this->NbrParticlesPerSpin-1);
  this->OrbitalsUp = OrbitalFactoryUp->CalculateFromSpinorVariables(part);
  part = uv.Extract(2*this->NbrParticlesPerSpin, 4*this->NbrParticlesPerSpin-1);
  this->OrbitalsDown = OrbitalFactoryDown->CalculateFromSpinorVariables(part);
  this->EvaluateTables();

  Complex Ji;
  Complex Tmp;
  // initialize Slater determinants
  for (int i=0;i<this->NbrParticlesPerSpin;++i)
    {
      Ji=1.0;
      for(int j=0;j<i;++j) Ji*=Jij[i][j];
      for(int j=i+1;j<NbrParticlesPerSpin;++j) Ji*=Jij[i][j];
      Ji*=SlaterElementNorm;
      for(int j=0;j<NbrParticlesPerSpin;++j)
	{
	  OrbitalsUp.GetMatrixElement(i,j,Tmp);
	  Tmp*=Ji;	  
	  SlaterUp->SetMatrixElement(i,j, Tmp);
	}
    }
  this->DeterminantUpValue=SlaterUp->Determinant();

  // initialize Slater determinants
  for (int i=this->NbrParticlesPerSpin;i<2*this->NbrParticlesPerSpin;++i)
    {
      Ji=1.0;
      for(int j=this->NbrParticlesPerSpin;j<i;++j) Ji*=Jij[i][j];
      for(int j=i+1;j<NbrParticles;++j) Ji*=Jij[i][j];
      Ji*=SlaterElementNorm;
      for(int j=0;j<this->NbrParticlesPerSpin;++j)
	{
	  OrbitalsDown.GetMatrixElement(i-this->NbrParticlesPerSpin,j,Tmp);
	  Tmp*=Ji;	  
	  SlaterDown->SetMatrixElement(i-this->NbrParticlesPerSpin,j, Tmp);
	}
    }
  this->DeterminantDownValue=SlaterDown->Determinant();
  
  return DeterminantUpValue*InterpolationUp*DeterminantDownValue*InterpolationDown*InterSpinJastrow;
}



// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void TwoThirdUnpolarizedCF::AdaptNorm(RealVector& x)
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

void TwoThirdUnpolarizedCF::AdaptAverageMCNorm(int thermalize, int average)
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
void TwoThirdUnpolarizedCF::EvaluateTables()
{
  int i, j;
  Complex tmp;
  // evaluate single particle Jastrow factors
  this->InterpolationUp=1.0;
  this->InterpolationDown=1.0;
  if (OrbitalFactoryUp->TestCriticality(InterpolationUp) == 0)
    {
      for (i=0;i<this->NbrParticlesPerSpin;i++)
	{
	  for(j=0;j<i;j++)
	    {
	      Jij[i][j] = OrbitalFactoryUp->JastrowFactorElement(i,j);
	      Jij[j][i] = - Jij[i][j];
	    }
	  Jij[i][j] = 0.0;
	  for(j=i+1;j<NbrParticlesPerSpin;j++)
	    {
	      Jij[i][j] = OrbitalFactoryUp->JastrowFactorElement(i,j);
	      Jij[j][i] = - Jij[i][j];
	    }
	}
    }
  else // if some interpolation occurred, the true values of the Ji's have to be recalculated:
    {
      cout << "Critical Event" << endl;
      for (i=0;i<this->NbrParticlesPerSpin;i++)
       	{	  
       	  for(j=0;j<i;j++)
	    {
	      Jij[i][j] = ((OrbitalFactoryUp->SpinorU(i) * OrbitalFactoryUp->SpinorV(j)) - (OrbitalFactoryUp->SpinorU(j) * OrbitalFactoryUp->SpinorV(i)));
	      Jij[j][i] = - Jij[i][j];
	    }
       	  for(j=i+1;j<NbrParticlesPerSpin;j++)
	    {
	      Jij[i][j] = ((OrbitalFactoryUp->SpinorU(i) * OrbitalFactoryUp->SpinorV(j)) - (OrbitalFactoryUp->SpinorU(j) * OrbitalFactoryUp->SpinorV(i)));
	      Jij[j][i] = - Jij[i][j];
	    }
      	}
    }

  if (OrbitalFactoryDown->TestCriticality(InterpolationDown) == 0)
    {
      for (i=0;i<this->NbrParticlesPerSpin;i++)
	{
	  for(j=0;j<i;j++)
	    {
	      Jij[i+NbrParticlesPerSpin][j+NbrParticlesPerSpin] = OrbitalFactoryDown->JastrowFactorElement(i,j);
	      Jij[j+NbrParticlesPerSpin][i+NbrParticlesPerSpin] = - Jij[i+NbrParticlesPerSpin][j+NbrParticlesPerSpin];
	    }
	  for(j=i+1;j<NbrParticlesPerSpin;j++)
	    {
	      Jij[i+NbrParticlesPerSpin][j+NbrParticlesPerSpin] = OrbitalFactoryDown->JastrowFactorElement(i,j);
	      Jij[j+NbrParticlesPerSpin][i+NbrParticlesPerSpin] = - Jij[i+NbrParticlesPerSpin][j+NbrParticlesPerSpin];
	    }
	}
    }
  else // if some interpolation occurred, the true values of the Ji's have to be recalculated:
    {
      cout << "Critical Event" << endl;
      for (i=0;i<this->NbrParticlesPerSpin;i++)
       	{	  
       	  for(j=0;j<i;j++)
	    {
	      Jij[i+NbrParticlesPerSpin][j+NbrParticlesPerSpin] = ((OrbitalFactoryDown->SpinorU(i) * OrbitalFactoryDown->SpinorV(j)) - (OrbitalFactoryDown->SpinorU(j) * OrbitalFactoryDown->SpinorV(i)));
	      Jij[j][i] = - Jij[i][j];
	    }
       	  for(j=i+1;j<NbrParticlesPerSpin;j++)
	    {
	      Jij[i+NbrParticlesPerSpin][j+NbrParticlesPerSpin] = ((OrbitalFactoryDown->SpinorU(i) * OrbitalFactoryDown->SpinorV(j)) - (OrbitalFactoryDown->SpinorU(j) * OrbitalFactoryDown->SpinorV(i)));
	      Jij[j][i] = - Jij[i][j];
	    }
      	}
    }

  for (i=0;i<this->NbrParticlesPerSpin;i++)
       	{
       	  for(j=0;j<NbrParticlesPerSpin;j++)
	    {
	      Jij[i][j+NbrParticlesPerSpin] = ((OrbitalFactoryUp->SpinorU(i) * OrbitalFactoryDown->SpinorV(j)) - (OrbitalFactoryDown->SpinorU(j) * OrbitalFactoryUp->SpinorV(i)));
	      Jij[j+NbrParticlesPerSpin][i] = - Jij[i][j+NbrParticlesPerSpin];
	    }
      	}


  double Factor = M_PI * 0.5;

  this->InterSpinJastrow=1.0;
  for (int i=0; i<this->NbrParticlesPerSpin;++i)
    for (int j=this->NbrParticlesPerSpin; j<2*this->NbrParticlesPerSpin;++j)
      {
	this->InterSpinJastrow*=(-Factor*Jij[i][j]);
      }
  this->InterSpinJastrow*=this->InterSpinJastrow; // take square

  // testing Jastrow-factors:
  Complex *SpinorU=new Complex[NbrParticles];
  Complex *SpinorV=new Complex[NbrParticles];
  for (int i=0; i<NbrParticlesPerSpin; ++i)
    {
      SpinorU[i]=OrbitalFactoryUp->SpinorU(i);
      if (Norm(SpinorU[i]-this->SpinorUCoordinates[i])>1e-10)
	{
	  cout << "problem with SpinorU["<<i<<"]"<<endl;
	}
      SpinorV[i]=OrbitalFactoryUp->SpinorV(i);
      if (Norm(SpinorV[i]-this->SpinorVCoordinates[i])>1e-10)
	{
	  cout << "problem with SpinorV["<<i<<"]"<<endl;
	}
      SpinorU[NbrParticlesPerSpin+i]=OrbitalFactoryDown->SpinorU(i);
      if (Norm(SpinorU[NbrParticlesPerSpin+i]-this->SpinorUCoordinates[NbrParticlesPerSpin+i])>1e-10)
	{
	  cout << "problem with SpinorU["<<NbrParticlesPerSpin+i<<"]"<<endl;
	}
      SpinorV[NbrParticlesPerSpin+i]=OrbitalFactoryDown->SpinorV(i);
      if (Norm(SpinorV[NbrParticlesPerSpin+i]-this->SpinorVCoordinates[NbrParticlesPerSpin+i])>1e-10)
	{
	  cout << "problem with SpinorV["<<NbrParticlesPerSpin+i<<"]"<<endl;
	}
    }
  for (int i=0; i<NbrParticles; ++i)
    for (int j=0; j<i; ++j)
      {
	Complex newJij=SpinorU[i]*SpinorV[j]-SpinorU[j]*SpinorV[i];
	if (Norm(newJij-Jij[i][j])>1e-10)
	  cout << "Jastrow factors "<<i<<","<<j<<" differ: "<< newJij << " vs "<<Jij[i][j]<<endl;
	
      }

  delete [] SpinorU;
  delete [] SpinorV;
    
}


      
