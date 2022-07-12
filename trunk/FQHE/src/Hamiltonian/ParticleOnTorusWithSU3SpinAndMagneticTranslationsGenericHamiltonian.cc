////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//      two body pseudopotential interaction and magnetic translations        //
//                                                                            //
//                        last modification : 07/03/2014                      //
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


#include "Hamiltonian/ParticleOnTorusWithSU3SpinAndMagneticTranslationsGenericHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Architecture/AbstractArchitecture.h"

#include "Polynomial/SpecialPolynomial.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


#define M1_12 0.08333333333333333


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// ratio = ratio between the width in the x direction and the width in the y direction
// nbrPseudoPotentials = array with the number of pseudo-potentials per interaction type
// pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as 11, 12, 13, 22, 23, 33)
// spinFlux1 = additional inserted flux for spin 1
// spinFlux2 = additional inserted flux for spin 2
// spinFlux3 = additional inserted flux for spin 3
// oneBodyPotentials = optional arry with the one body potentials (sorted as 11, 12, 13, 22, 23, 33)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusWithSU3SpinAndMagneticTranslationsGenericHamiltonian::ParticleOnTorusWithSU3SpinAndMagneticTranslationsGenericHamiltonian(ParticleOnTorusWithSU3SpinAndMagneticTranslations* particles, 
																	 int nbrParticles, int maxMomentum, int xMomentum, 
																	 double ratio, 
																	 int* nbrPseudopotentials, 
																	 double** pseudoPotentials, 
																	 double** oneBodyPotentials, 
																	 double spinFlux1, double spinFlux2, double spinFlux3,
																	 AbstractArchitecture* architecture, long memory, 
																	 char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->Ratio = ratio;  
  this->InvRatio = 1.0 / ratio;
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  

  this->SpinFlux1 = spinFlux1;
  this->SpinFlux2 = spinFlux2;
  this->SpinFlux3 = spinFlux3;

  this->NbrPseudopotentials  = new int [6];
  this->Pseudopotentials = new double*[6];
  this->MaxNbrPseudopotentials = 0;
  for (int i = 0; i < 6; ++i)
    {
      this->NbrPseudopotentials[i] = nbrPseudopotentials[i];
      if (this->NbrPseudopotentials[i] > 0)
	{
	  this->Pseudopotentials[i] =  new double[this->NbrPseudopotentials[i]];
	  for (int j = 0; j < this->NbrPseudopotentials[i]; ++j)
	    {
	      this->Pseudopotentials[i][j] = pseudoPotentials[i][j];
	    }
	}
      if (this->MaxNbrPseudopotentials < this->NbrPseudopotentials[i])
	this->MaxNbrPseudopotentials = this->NbrPseudopotentials[i];
    } 
  this->LaguerrePolynomials = new Polynomial[this->MaxNbrPseudopotentials];
  for (int i = 0; i < this->MaxNbrPseudopotentials; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

  this->OneBodyInteractionFactors11 = 0;
  this->OneBodyInteractionFactors12 = 0;
  this->OneBodyInteractionFactors13 = 0;
  this->OneBodyInteractionFactors22 = 0;
  this->OneBodyInteractionFactors23 = 0;
  this->OneBodyInteractionFactors33 = 0;
  if (oneBodyPotentials != 0)
    {
      if (oneBodyPotentials[0] != 0)
	{
	  this->OneBodyInteractionFactors11 = new double[this->NbrLzValue];
	  for (int i = 0; i < this->NbrLzValue; i++)
	    this->OneBodyInteractionFactors11[i] = oneBodyPotentials[0][i];
	}
      if (oneBodyPotentials[3] != 0)
	{
	  this->OneBodyInteractionFactors22 = new double[this->NbrLzValue];
	  for (int i = 0; i < this->NbrLzValue; i++)
	    this->OneBodyInteractionFactors22[i] = oneBodyPotentials[3][i];
	}
      if (oneBodyPotentials[5] != 0)
	{
	  this->OneBodyInteractionFactors33 = new double[this->NbrLzValue];
	  for (int i = 0; i < this->NbrLzValue; i++)
	    this->OneBodyInteractionFactors33[i] = oneBodyPotentials[5][i];
	}
      if(oneBodyPotentials[1] != 0)
	{
	  this->OneBodyInteractionFactors12 = new Complex[this->NbrLzValue];
	  for (int i = 0; i < this->NbrLzValue; i++)
	    {
	      this->OneBodyInteractionFactors12[i] = oneBodyPotentials[1][i];
	    }
	} 
      if(oneBodyPotentials[2] != 0)
	{
	  this->OneBodyInteractionFactors13 = new Complex[this->NbrLzValue];
	  for (int i = 0; i < this->NbrLzValue; i++)
	    {
	      this->OneBodyInteractionFactors13[i] = oneBodyPotentials[2][i];
	    }
	} 
      if(oneBodyPotentials[4] != 0)
	{
	  this->OneBodyInteractionFactors23 = new Complex[this->NbrLzValue];
	  for (int i = 0; i < this->NbrLzValue; i++)
	    {
	      this->OneBodyInteractionFactors23[i] = oneBodyPotentials[4][i];
	    }
	} 
    }
  this->EvaluateExponentialFactors();
  this->EvaluateInteractionFactors();

  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  if (TmpMemory < 1024l)
	    cout  << "fast = " <<  TmpMemory << "b ";
	  else
	    if (TmpMemory < (1l << 20))
	      cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	    else
	      if (TmpMemory < (1l << 30))
		cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	      else
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
	  cout << endl;
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

ParticleOnTorusWithSU3SpinAndMagneticTranslationsGenericHamiltonian::~ParticleOnTorusWithSU3SpinAndMagneticTranslationsGenericHamiltonian() 
{
  if (this->Pseudopotentials != 0)
    {
      for (int i = 0; i < 6; ++i)
	delete[] this->Pseudopotentials[i];
      delete[] this->Pseudopotentials;
    }
  if (this->LaguerrePolynomials != 0)
    delete[] this->LaguerrePolynomials;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusWithSU3SpinAndMagneticTranslationsGenericHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->FastMultiplicationFlag = false;
  this->Particles = (ParticleOnTorusWithSU3SpinAndMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusWithSU3SpinAndMagneticTranslationsGenericHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnTorusWithSU3SpinAndMagneticTranslationsGenericHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  long TotalNbrNonZeroInteractionFactors = 0;
  double MaxCoefficient = 0.0;
  this->GetIndices();
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      // 11-11
      this->InteractionFactors11 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors11[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										  this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										  this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors11[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors11[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // 22-22
      this->InteractionFactors22 = new Complex* [this->NbrIntraSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors22[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										  this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										  this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors22[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors22[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // 33-33
      this->InteractionFactors33 = new Complex* [this->NbrIntraSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors33[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										  this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										  this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors33[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors33[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // 12-12
      this->InteractionFactors12 = new Complex* [this->NbrInterSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors12[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[1], this->Pseudopotentials[1], 
										    this->SpinFlux1, this->SpinFlux2, this->SpinFlux2, this->SpinFlux1)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[1], this->Pseudopotentials[1], 
										    this->SpinFlux2, this->SpinFlux1, this->SpinFlux1, this->SpinFlux2));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[1], this->Pseudopotentials[1], 
										    this->SpinFlux1, this->SpinFlux2, this->SpinFlux2, this->SpinFlux1)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[1], this->Pseudopotentials[1], 
										    this->SpinFlux2, this->SpinFlux1, this->SpinFlux1, this->SpinFlux2));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors12[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors12[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // 13-13
      this->InteractionFactors13 = new Complex* [this->NbrInterSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors13[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[2], this->Pseudopotentials[2], 
										    this->SpinFlux1, this->SpinFlux3, this->SpinFlux3, this->SpinFlux1)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[2], this->Pseudopotentials[2], 
										    this->SpinFlux3, this->SpinFlux1, this->SpinFlux1, this->SpinFlux3));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[2], this->Pseudopotentials[2], 
										    this->SpinFlux1, this->SpinFlux2, this->SpinFlux3, this->SpinFlux1)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[2], this->Pseudopotentials[2], 
										    this->SpinFlux2, this->SpinFlux1, this->SpinFlux1, this->SpinFlux3));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors13[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors13[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // 23-23
      this->InteractionFactors23 = new Complex* [this->NbrInterSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors23[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[4], this->Pseudopotentials[4], 
										    this->SpinFlux2, this->SpinFlux3, this->SpinFlux3, this->SpinFlux2)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[4], this->Pseudopotentials[4], 
										    this->SpinFlux3, this->SpinFlux2, this->SpinFlux2, this->SpinFlux3));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[4], this->Pseudopotentials[4], 
										    this->SpinFlux2, this->SpinFlux3, this->SpinFlux3, this->SpinFlux2)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[4], this->Pseudopotentials[4], 
										    this->SpinFlux3, this->SpinFlux2, this->SpinFlux2, this->SpinFlux3));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors23[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors23[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      // 11-11 
      this->InteractionFactors11 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors11[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										  this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1));
		  if (m1 == m2)		    
		    TmpCoefficient *= 0.5;
		  if (m3 == m4)		    
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										  this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0], 
										    this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1));
		  if (m1 == m2)		    
		    TmpCoefficient *= 0.5;
		  if (m3 == m4)		    
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors11[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors11[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // 22-22
      this->InteractionFactors22 = new Complex* [this->NbrIntraSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors22[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										  this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2));
		  if (m1 == m2)		    
		    TmpCoefficient *= 0.5;
		  if (m3 == m4)		    
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										  this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3], 
										    this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2));
		  if (m1 == m2)		    
		    TmpCoefficient *= 0.5;
		  if (m3 == m4)		    
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors22[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors22[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // 33-33
      this->InteractionFactors33 = new Complex* [this->NbrIntraSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors33[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										  this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3));
		  if (m1 == m2)		    
		    TmpCoefficient *= 0.5;
		  if (m3 == m4)		    
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										  this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5], 
										    this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3));
		  if (m1 == m2)		    
		    TmpCoefficient *= 0.5;
		  if (m3 == m4)		    
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors33[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors33[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // 12-12
      this->InteractionFactors12 = new Complex* [this->NbrInterSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors12[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[1], this->Pseudopotentials[1], 
										    this->SpinFlux1, this->SpinFlux2, this->SpinFlux2, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[1], this->Pseudopotentials[1], 
										    this->SpinFlux2, this->SpinFlux1, this->SpinFlux1, this->SpinFlux2));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[1], this->Pseudopotentials[1], 
										    this->SpinFlux1, this->SpinFlux2, this->SpinFlux2, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[1], this->Pseudopotentials[1], 
										    this->SpinFlux2, this->SpinFlux1, this->SpinFlux1, this->SpinFlux2));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors12[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors12[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // 13-13
      this->InteractionFactors13 = new Complex* [this->NbrInterSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors13[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[2], this->Pseudopotentials[2], 
										    this->SpinFlux1, this->SpinFlux3, this->SpinFlux3, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[2], this->Pseudopotentials[2], 
										    this->SpinFlux3, this->SpinFlux1, this->SpinFlux1, this->SpinFlux3));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[2], this->Pseudopotentials[2], 
										    this->SpinFlux1, this->SpinFlux3, this->SpinFlux3, this->SpinFlux1)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[2], this->Pseudopotentials[2], 
										    this->SpinFlux3, this->SpinFlux1, this->SpinFlux1, this->SpinFlux3));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors13[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors13[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // 23-23
      this->InteractionFactors23 = new Complex* [this->NbrInterSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors23[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[4], this->Pseudopotentials[4], 
										    this->SpinFlux2, this->SpinFlux3, this->SpinFlux3, this->SpinFlux2)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[4], this->Pseudopotentials[4], 
										    this->SpinFlux3, this->SpinFlux2, this->SpinFlux2, this->SpinFlux3));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[4], this->Pseudopotentials[4], 
										    this->SpinFlux2, this->SpinFlux3, this->SpinFlux3, this->SpinFlux2)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[4], this->Pseudopotentials[4], 
										    this->SpinFlux3, this->SpinFlux2, this->SpinFlux2, this->SpinFlux3));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors23[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactors23[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "nbr non-zero interaction = " << TotalNbrNonZeroInteractionFactors << endl;
  cout << "====================================" << endl;
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = pseudopotential coefficients
// spinFluxM1 = additional inserted flux for m1
// spinFluxM2 = additional inserted flux for m2
// spinFluxM3 = additional inserted flux for m3
// spinFluxM4 = additional inserted flux for m4
// return value = numerical coefficient

double ParticleOnTorusWithSU3SpinAndMagneticTranslationsGenericHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
													   double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->NbrLzValue);
  double Factor =  - (((double) (m1-m3)) + spinFluxM1 - spinFluxM3) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = ((double) (m1 - m4)) + spinFluxM1 - spinFluxM4;
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
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0* PIOnM * Q2);
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
          if (fabs(Coefficient) != 0.0)
 	    Precision = Coefficient;
          else
            Precision = 1.0;
	}
       else
 	{
	  Precision = 1.0;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
	  Precision = 2.0 * exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->NbrLzValue;
    }
  N2 = (double) (m1 - m4 - this->NbrLzValue) + spinFluxM1 - spinFluxM4;
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  TmpInteraction = 0.0;
	  for (int i=0; i< nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
          if (fabs(Coefficient) != 0.0)
	    Precision = Coefficient;
          else
            Precision = 1.0;
	}
       else
 	{
	  Precision = 1.0;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
	  Precision = 2.0 *  exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->NbrLzValue;
    }
  //Normalize per flux (gives correct energy scale for 2-particle problem)
  return (Sum / ((double) this->NbrLzValue));
}

