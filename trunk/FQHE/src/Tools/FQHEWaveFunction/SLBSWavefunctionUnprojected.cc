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
#include "SLBSWavefunctionUnprojected.h"
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

SLBSWavefunctionUnprojected::SLBSWavefunctionUnprojected()
{
  this->NbrParticles = 0;
}

// constructor
//
// nbrParticles = number of particles
// nbrLandauLevels = number of LL's of CF orbitals to fill (unprojected)
// negFluxFlag = indicating sign of Jain wavefunction component
//
SLBSWavefunctionUnprojected::SLBSWavefunctionUnprojected(int nbrParticles, int nbrLandauLevels, bool negFluxFlag)
{
  if ((nbrParticles & 1)||(nbrParticles < 4))
    {
      cout << "This state requires an even number of fermions >= 4" << endl;
      exit(1);
    }
  this->NbrParticles = nbrParticles;
  this->NbrLandauLevels=nbrLandauLevels;
  this->NegativeFieldFlag = negFluxFlag;
  this->ActualJastrowPower=2;
  
  this->TwiceS = (this->NbrParticles / this->NbrLandauLevels) - this->NbrLandauLevels;
  this->EffectiveFlux = this->TwiceS*(1-2*NegativeFieldFlag);
  cout << "EffectiveFlux="<<EffectiveFlux<<endl;

  this->Flag.Initialize();
  //this->Interpolation=1.0;
  this->Slater = new ComplexMatrix(this->NbrParticles,this->NbrParticles);
  this->Pfaffian = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->SlaterElementNorm = 1.0;
  
  
  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];

  // for unprojected orbitals:
  this->EvaluateNormalizationPrefactors();  
  this->EvaluateSumPrefactors();
  this->JastrowFactorElements = new Complex*[this->NbrParticles];
  for (int i = 1; i < this->NbrParticles; ++i)
    this->JastrowFactorElements[i] = new Complex[i];
  this->SpinorUCoordinatePower = new Complex*[NbrParticles];
  this->SpinorVCoordinatePower = new Complex*[NbrParticles];
  int MaxPower = this->TwiceS + 2 * (this->NbrLandauLevels - 1);
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinatePower[i] = new Complex[MaxPower + 1];
      this->SpinorVCoordinatePower[i] = new Complex[MaxPower + 1];
    }
}

// copy constructor
//
// function = reference on the wave function to copy

SLBSWavefunctionUnprojected::SLBSWavefunctionUnprojected(const SLBSWavefunctionUnprojected& function)
{
  // to adapt
  this->NbrParticles = function.NbrParticles;
  this->EffectiveFlux = function.EffectiveFlux;
  this->SlaterElementNorm=function.SlaterElementNorm;
  this->Determinant = function.Determinant;
  this->NegativeFieldFlag = function.NegativeFieldFlag;
  this->EffectiveFlux = function.EffectiveFlux;
  this->Flag = function.Flag;
  //this->Interpolation=function.Interpolation;
  this->Slater = new ComplexMatrix(this->NbrParticles,this->NbrParticles);
  this->Pfaffian = new ComplexSkewSymmetricMatrix(this->NbrParticles);
    

  this->SpinorUCoordinates = new Complex[NbrParticles];
  this->SpinorVCoordinates = new Complex[NbrParticles];
}

// destructor
//

SLBSWavefunctionUnprojected::~SLBSWavefunctionUnprojected()
{
  if (this->NbrParticles>0)
    {
      if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
	{
	  
	}
      delete Slater;
      delete Pfaffian;
      delete [] SpinorUCoordinates;
      delete [] SpinorVCoordinates;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* SLBSWavefunctionUnprojected::Clone ()
{
  return new SLBSWavefunctionUnprojected(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex SLBSWavefunctionUnprojected::operator ()(RealVector& x)
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
  this->EvaluateTables(false);
  
  Complex Tmp;

  // initialize Slater determinant and Pfaffian
  for (int i=0;i<this->NbrParticles;++i)
    for(int j=0;j<i;++j)
      Pfaffian->SetMatrixElement(i,j,JastrowFactorElements[i][j] /* testing*/ *SlaterElementNorm);
  
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      int Index = 0;
      int MaxMomentum = this->TwiceS;
      for (int j = 0; j < this->NbrLandauLevels; ++j)
	{
	  for (int k = 0; k <= MaxMomentum; ++k)
	    {	      
	      Slater->SetMatrixElement(Index, i, this->EvaluateMonopoleHarmonic(i, k, j, MaxMomentum)*SlaterElementNorm);
	      ++Index;
	    }
	  MaxMomentum += 2;
	}
    }  
  
  this->Determinant=Slater->Determinant();
  this->PfaffianValue=Pfaffian->Pfaffian();

  return this->PfaffianValue*this->Chi1*this->Determinant;
  
  // testing validity for 2/5 state:
  // return this->Determinant*this->Chi1*this->Chi1;

  // testing validity for MR state:
  // return this->PfaffianValue*this->Chi1*this->Chi1;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex SLBSWavefunctionUnprojected::CalculateFromSpinorVariables(ComplexVector& uv)
{
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = uv.Re(2*i);
      this->SpinorUCoordinates[i].Im = uv.Im(2*i);
      this->SpinorVCoordinates[i].Re = uv.Re(2*i+1);
      this->SpinorVCoordinates[i].Im = uv.Im(2*i+1);
    }

  this->EvaluateTables(false);
  
  Complex Tmp;

  // initialize Slater determinant and Pfaffian
  for (int i=0;i<this->NbrParticles;++i)
    for(int j=0;j<i;++j)
      Pfaffian->SetMatrixElement(i,j,1.0/Jij[i][j]);
  
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      int Index = 0;
      int MaxMomentum = this->TwiceS;
      for (int j = 0; j < this->NbrLandauLevels; ++j)
	{
	  for (int k = 0; k <= MaxMomentum; ++k)
	    {	      
	      Slater->SetMatrixElement(Index, i, this->EvaluateMonopoleHarmonic(i, k, j, MaxMomentum)*SlaterElementNorm);
	      ++Index;
	    }
	  MaxMomentum += 2;
	}
    }  
  
  this->Determinant=Slater->Determinant();
  this->PfaffianValue=Pfaffian->Pfaffian();

  return this->PfaffianValue*this->Chi1*this->Determinant;
  
  // testing validity for 2/5 state:
  // return this->Determinant*this->Chi1*this->Chi1;
}



// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void SLBSWavefunctionUnprojected::AdaptNorm(RealVector& x)
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

void SLBSWavefunctionUnprojected::AdaptAverageMCNorm(int thermalize, int average)
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


      

// evaluate monopole spherical harmonic 
//
// coordinate = index of the coordinate
// momentum = monopole spherical harmonic Lz momentum (plus S shift)
// landauLevel = index of the Landau level
// maximumMomentum = maxixum momentum that can be reached in the Landau level
// return value = value of the monopole spherical harmonic at the given point

Complex SLBSWavefunctionUnprojected::EvaluateMonopoleHarmonic (int coordinate, int momentum, int landauLevel, int maximumMomentum)
{
  Complex Tmp(0.0);
  int i = landauLevel - momentum;
  if (i < 0)
    i = 0;
  int UPower = momentum + i - landauLevel;
  int Max = maximumMomentum - momentum;
  int VPower = Max - i;
  if (Max > landauLevel)
    Max = landauLevel;
  if (this->NegativeFieldFlag == false)
    {
      for (; i <= Max; ++i)
	{
	  Tmp += (this->SumPrefactors[landauLevel][i][momentum] * this->SpinorUCoordinatePower[coordinate][UPower] * 
		  this->SpinorVCoordinatePower[coordinate][VPower] * Conj(this->SpinorUCoordinatePower[coordinate][i]) * 
		  Conj(this->SpinorVCoordinatePower[coordinate][landauLevel - i]));
	  ++UPower;
	  --VPower;
	}
    }
  else
    {
      for (; i <= Max; ++i)
	{
	  Tmp += (this->SumPrefactors[landauLevel][i][momentum] * Conj(this->SpinorUCoordinatePower[coordinate][UPower]) * 
		  Conj(this->SpinorVCoordinatePower[coordinate][VPower]) * this->SpinorUCoordinatePower[coordinate][i] * 
		  this->SpinorVCoordinatePower[coordinate][landauLevel - i]);
	  ++UPower;
	  --VPower;
	}
    }
  return (this->NormalizationPrefactors[landauLevel][momentum] * Tmp);
}


// evaluate constant factors that appears in the sum of projected monopole harmonic (except LLL)
//

void SLBSWavefunctionUnprojected::EvaluateSumPrefactors()
{
  this->SumPrefactors = new double** [this->NbrLandauLevels];
  int MaxMomentum = this->TwiceS;
  FactorialCoefficient Coef;
  double Factor = 1.0;
  for (int i = 0; i < this->NbrLandauLevels; ++i)  
    {
      this->SumPrefactors[i] = new double* [i + 1];
      Factor = 1.0;
      for (int k = 0; k <= i; ++k)  
	{
	  this->SumPrefactors[i][k] = new double [MaxMomentum + 1];
	  for (int j = 0; j < (i - k); ++j)
	    this->SumPrefactors[i][k][j] = 0.0;
	  for (int j = i - k; j <= (MaxMomentum - k); ++j)
	    {
	      Coef.SetToOne();
	      Coef.PartialFactorialMultiply(k + 1, i);
	      Coef.FactorialDivide(i - k);
	      Coef.PartialFactorialMultiply(j + k - i + 1, this->TwiceS + i);	  
	      Coef.FactorialDivide(this->TwiceS + (2 * i) -j - k);	  
	      this->SumPrefactors[i][k][j] = Factor * Coef.GetNumericalValue();
	    }
	  for (int j = MaxMomentum - k + 1; j <= MaxMomentum; ++j)
	    this->SumPrefactors[i][k][j] = 0.0;
	  Factor *= -1.0;
	}
      MaxMomentum += 2;
    }
}


// evaluate normalization factors of projected monopole harmonics
//

void SLBSWavefunctionUnprojected::EvaluateNormalizationPrefactors()
{
  this->NormalizationPrefactors = new double* [this->NbrLandauLevels];
  int MaxMomentum = this->TwiceS;
  FactorialCoefficient Coef;
  double Factor = ((double) (this->TwiceS + 1)) / (4.0  * M_PI);
  this->NormalizationPrefactors[0] = new double[MaxMomentum + 1];
  for (int j = 0; j <= MaxMomentum; ++j)
    {
      Coef.SetToOne();
      Coef.PartialFactorialDivide(this->TwiceS - j + 1, this->TwiceS);
      Coef.FactorialMultiply(j);
      this->NormalizationPrefactors[0][j] = sqrt(Factor * Coef.GetNumericalValue());
      if ((j & 1) != 0)
	{
	  this->NormalizationPrefactors[0][j] *= -1.0;
	}
    }
  MaxMomentum += 2;
  for (int i = 1; i < this->NbrLandauLevels; ++i)  
    {
      double Factor = ((double) (this->TwiceS + (2 * i) + 1)) / (4.0  * M_PI);
      this->NormalizationPrefactors[i] = new double[MaxMomentum + 1];
      double Sign = 1.0;
      if ((this->TwiceS & 1) != 0)
	Sign = -1.0;
      for (int j = 0; j <= MaxMomentum; ++j)
	{	  
	  Coef.SetToOne();
	  Coef.FactorialMultiply(this->TwiceS + 2 * i - j);
	  Coef.FactorialDivide(i);
	  Coef.FactorialMultiply(j);
	  Coef.FactorialDivide(this->TwiceS + i);
	  this->NormalizationPrefactors[i][j] = sqrt(Factor * Coef.GetNumericalValue());
	  this->NormalizationPrefactors[i][j] *= Sign;
	  Sign *= -1.0;
	}
      MaxMomentum += 2;
    }
}

// evaluate precalculation tables used during wave function evaluation (called at each evaluation)
//
// requires SpinorUCoordinates and SpinorVCoordinates to be initialized prior to call
// derivativeFlag = indicate if precalculation tables invloved in derivative evaluation have to be calculated
// return value = value of the Jastrow factor

Complex SLBSWavefunctionUnprojected::EvaluateTables(bool derivativeFlag)
{  

  int MaxPower = this->TwiceS + 2 * (this->NbrLandauLevels - 1);
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinatePower[i][0] = 1.0;
      this->SpinorVCoordinatePower[i][0] = 1.0;
      Complex TmpU = this->SpinorUCoordinates[i];
      Complex TmpV = this->SpinorVCoordinates[i];
      for (int j = 1; j <= MaxPower; ++j)
	{
	  this->SpinorUCoordinatePower[i][j] = this->SpinorUCoordinatePower[i][j - 1] * TmpU;
	  this->SpinorVCoordinatePower[i][j] = this->SpinorVCoordinatePower[i][j - 1] * TmpV;
	}
    }

  Complex JastrowFactor(1.0);
  Complex Tmp;
  for (int i = 1; i < this->NbrParticles; ++i)
    {
      for (int j = 0; j < i; ++j)
	{
	  Tmp = ((this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]));
	  this->JastrowFactorElements[i][j] = 1.0 / Tmp;
	  JastrowFactor *= Tmp;
	}
    }
  this->Chi1=JastrowFactor;
  for (int i = 1; i < this->ActualJastrowPower; ++i)
    {
      JastrowFactor *= this->Chi1;
    }

  if (derivativeFlag == true)
    {
      Complex Tmp2;
      for (int i = 0; i < this->NbrParticles; ++i)
	{     
	  Complex** TmpDerivativeFactors = this->DerivativeFactors[i];
	  Complex** TmpDerivativeFactors2 = this->DerivativeFactors2[i];
	  for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
	    for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
	      TmpDerivativeFactors[k1][k2] = 0.0;
	  
	  int Index = 0;
	  for (int j = 1; j < this->NbrParticles; ++j)
	    {
	      if (Index == i)
		++Index;
	      if (Index > i)
		Tmp2 = this->SpinorVCoordinates[Index] * this->JastrowFactorElements[Index][i];
	      else
		Tmp2 = -this->SpinorVCoordinates[Index] * this->JastrowFactorElements[i][Index];
	      Tmp = 1.0;
	      for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
		{
		  for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
		    {
		      TmpDerivativeFactors2[k1][k2] = Tmp;
		    }
		  Tmp *= Tmp2;
		}
	      if (Index > i)
		Tmp2 = - this->SpinorUCoordinates[Index] * this->JastrowFactorElements[Index][i];
	      else
		Tmp2 = this->SpinorUCoordinates[Index] * this->JastrowFactorElements[i][Index];
	      Tmp = 1.0;
	      for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
		{
		  for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
		    {
		      TmpDerivativeFactors2[k2][k1] *= Tmp;
		    }
		  Tmp *= Tmp2;
		}
	      for (int k1 = 0; k1 < this->NbrLandauLevels; ++k1)
		for (int k2 = 0; k2 < this->NbrLandauLevels; ++k2)
		  TmpDerivativeFactors[k1][k2] += TmpDerivativeFactors2[k1][k2];
	      ++Index;
	    }
	}
    }
  return JastrowFactor;
}
