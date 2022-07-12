////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//       generic two body interaction without assuming ky conservation        //
//                                                                            //
//                        last modification : 11/09/2016                      //
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


#include "Hamiltonian/ParticleOnTorusGenericHamiltonianAllMomenta.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Polynomial/SpecialPolynomial.h"

#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


#define M1_12 0.08333333333333333


// default constructor
//

ParticleOnTorusGenericHamiltonianAllMomenta::ParticleOnTorusGenericHamiltonianAllMomenta()
{
  this->LaguerrePolynomials = 0;
  this->Pseudopotentials = 0;
  this->OneBodyTermFlag = false;
  this->OneBodyInteractionFactors = 0;
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// ratio = ratio between the width in the x direction and the width in the y direction
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = pseudopotential coefficients
// oneBodyPotentials = array that contains the diagonal terms of the one-body potential
// xFlux = additional inserted flux along the x axis in flux quantum unit
// yFlux = additional inserted flux along the y axis in flux quantum unit
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusGenericHamiltonianAllMomenta::ParticleOnTorusGenericHamiltonianAllMomenta(ParticleOnTorus* particles, int nbrParticles, int maxMomentum, double ratio, 
											 int nbrPseudopotentials, double* pseudopotentials,
											 double* oneBodyPotentials,
											 double xFlux, double yFlux,
											 AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->HamiltonianShift = 0.0;
  if (oneBodyPotentials == 0)
    {
      this->OneBodyTermFlag = false;
      this->OneBodyInteractionFactors = 0;
    }
  else
    {
      this->OneBodyTermFlag = true;
      this->OneBodyInteractionFactors = new double[this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	{
	  this->OneBodyInteractionFactors[i] = oneBodyPotentials[i];
	}
      
    }
  this->FilterInteractionFlag = false;
  this->XFlux = xFlux;
  this->YFlux = yFlux;

  this->NbrPseudopotentials = nbrPseudopotentials;
  this->Pseudopotentials = new double[this->NbrPseudopotentials];
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    this->Pseudopotentials[i] = pseudopotentials[i];
  this->LaguerrePolynomials =new Polynomial[this->NbrPseudopotentials];
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

  this->EvaluateInteractionFactors();

  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  if (TmpMemory < 1024)
	    cout  << "fast = " <<  TmpMemory << "b ";
	  else
	    if (TmpMemory < (1 << 20))
	      cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	    else
	      if (TmpMemory < (1 << 30))
		cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	      else
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnTorusGenericHamiltonianAllMomenta::~ParticleOnTorusGenericHamiltonianAllMomenta() 
{
  if (this->Pseudopotentials != 0)
    delete[] this->Pseudopotentials;
  if (this->LaguerrePolynomials != 0)
    delete[] this->LaguerrePolynomials;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusGenericHamiltonianAllMomenta::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  if (this->FastMultiplicationFlag == true)
    {
      for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
    }
  this->Particles = (ParticleOnTorus*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// evaluate all interaction factors
//   

void ParticleOnTorusGenericHamiltonianAllMomenta::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
    {
      this->GetIndices();
      this->InteractionFactors = new double* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new double[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];

		  this->InteractionFactors[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->XFlux, this->XFlux, this->XFlux, this->XFlux)
							+ this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->XFlux, this->XFlux, this->XFlux, this->XFlux)
							- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->XFlux, this->XFlux, this->XFlux, this->XFlux)
							- this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->XFlux, this->XFlux, this->XFlux, this->XFlux));
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      this->GetIndices();
      this->InteractionFactors = new double* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new double[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];

		  this->InteractionFactors[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->XFlux, this->XFlux, this->XFlux, this->XFlux)
							+ this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->XFlux, this->XFlux, this->XFlux, this->XFlux)
							+ this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->XFlux, this->XFlux, this->XFlux, this->XFlux)
							+ this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->XFlux, this->XFlux, this->XFlux, this->XFlux));
		  if (m3 == m4)
		    this->InteractionFactors[i][Index] *= 0.5;
		  if (m1 == m2)
		    this->InteractionFactors[i][Index] *= 0.5;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// fluxM1 = additional inserted flux for m1
// fluxM2 = additional inserted flux for m2
// fluxM3 = additional inserted flux for m3
// fluxM4 = additional inserted flux for m4
// return value = numerical coefficient

double ParticleOnTorusGenericHamiltonianAllMomenta::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4,
										   double fluxM1, double fluxM2, double fluxM3, double fluxM4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->NbrLzValue);
  double TwoPIOnM = 2.0 * M_PI / ((double) this->NbrLzValue);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Precision;
  double TmpInteraction;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  TmpInteraction = 0.0;
	  for (int i = 0; i < this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(TwoPIOnM * Q2);
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
	  Precision = Coefficient;
	}
      else
	{
	  Precision = 1.0;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(TwoPIOnM * Q2);
	  Precision = 2.0 * exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->NbrLzValue;
    }
  N2 = (double) (m1 - m4 - this->NbrLzValue);
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  TmpInteraction = 0.0;
	  for (int i=0; i< this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(TwoPIOnM * Q2);
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
	  Precision = Coefficient;
	}
      else
	{
	  Precision = 1.0;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(TwoPIOnM * Q2);
	  Precision = 2.0 *  exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->NbrLzValue;
    }
  return (Sum / this->NbrLzValue);
}

