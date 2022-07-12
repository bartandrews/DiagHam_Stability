////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//      two body pseudopotential interaction and magnetic translations        //
//                                                                            //
//                        last modification : 10/02/2015                      //
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


#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusWithSpinAndMagneticTranslationsGenericHamiltonian.h"
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
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


//default constructor

ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian::ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian (){
}

// constructor from pseudopotentials
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// ratio = ratio between the width in the x direction and the width in the y direction
// nbrPseudopotentialsUpUp = number of pseudopotentials for up-up interaction
// pseudopotentialsUpUp = pseudopotential coefficients for up-up interaction
// nbrPseudopotentialsDownDown = number of pseudopotentials for down-down interaction
// pseudopotentialsDownDown = pseudopotential coefficients for down-down interaction
// nbrPseudopotentialsUpDown = number of pseudopotentials for up-down interaction
// pseudopotentialsUpDown = pseudopotential coefficients for up-down interaction
// spinFluxUp = additional inserted flux for spin up
// spinFluxDown = additional inserted flux for spin down
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian::ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian (ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp, int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown, int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,  double spinFluxUp, double spinFluxDown,  AbstractArchitecture* architecture, long memory, char* precalculationFileName, double* oneBodyPotentielUpUp, double* oneBodyPotentielDownDown, double* oneBodyPotentielUpDown)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  if (this->Particles->GetTotalSpin() == 0)
    this->MomentumModulo = this->MaxMomentum;
  else
    this->MomentumModulo = FindGCD(this->Particles->GetTotalSpin(), this->MaxMomentum);
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = false;//true;
  this->Ratio = ratio;  
  this->InvRatio = 1.0 / ratio;
  this->SpinFluxUp = spinFluxUp;
  this->SpinFluxDown = spinFluxDown;
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  

  this->NbrPseudopotentialsUpUp = nbrPseudopotentialsUpUp;
  this->PseudopotentialsUpUp = new double[this->NbrPseudopotentialsUpUp];
  for (int i = 0; i < this->NbrPseudopotentialsUpUp; ++i)
    this->PseudopotentialsUpUp[i] = pseudopotentialsUpUp[i];
  this->NbrPseudopotentialsDownDown = nbrPseudopotentialsDownDown;
  this->PseudopotentialsDownDown = new double[this->NbrPseudopotentialsDownDown];
  for (int i = 0; i < this->NbrPseudopotentialsDownDown; ++i)
    this->PseudopotentialsDownDown[i] = pseudopotentialsDownDown[i];
  this->NbrPseudopotentialsUpDown = nbrPseudopotentialsUpDown;
  this->PseudopotentialsUpDown = new double[this->NbrPseudopotentialsUpDown];
  for (int i = 0; i < this->NbrPseudopotentialsUpDown; ++i)
    this->PseudopotentialsUpDown[i] = pseudopotentialsUpDown[i];

  this->MaxNbrPseudopotentials = this->NbrPseudopotentialsUpUp;
  if (this->NbrPseudopotentialsDownDown > this->MaxNbrPseudopotentials)
    this->MaxNbrPseudopotentials = this->NbrPseudopotentialsDownDown;
  if (this->NbrPseudopotentialsUpDown > this->MaxNbrPseudopotentials)
    this->MaxNbrPseudopotentials = this->NbrPseudopotentialsUpDown;
  this->LaguerrePolynomials = new Polynomial[this->MaxNbrPseudopotentials];
  for (int i = 0; i < this->MaxNbrPseudopotentials; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

  this->OneBodyInteractionFactorsupup = 0;
  if(oneBodyPotentielUpUp != 0)
    {
      this->OneBodyInteractionFactorsupup = new double[this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; i++)
	this->OneBodyInteractionFactorsupup[i] = oneBodyPotentielUpUp[i];
    }
  this->OneBodyInteractionFactorsdowndown = 0;
  if(oneBodyPotentielDownDown != 0)
    {
      this->OneBodyInteractionFactorsdowndown = new double[this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; i++)
	this->OneBodyInteractionFactorsdowndown[i] = oneBodyPotentielDownDown[i];
    } 
  this->OneBodyInteractionFactorsupdown = 0;
  if(oneBodyPotentielUpDown != 0)
    {
      this->OneBodyInteractionFactorsupdown = new Complex[this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; i++)
	{
	  this->OneBodyInteractionFactorsupdown[i] = oneBodyPotentielUpDown[i];
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

ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian::~ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian() 
{
}


  
// evaluate all interaction factors
//   

void ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  long TotalNbrNonZeroInteractionFactors = 0;
  double MaxCoefficient = 0.0;
  this->GetIndices();



  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      // upup-upup 
      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientUpUp(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										  this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     + this->EvaluateInteractionCoefficientUpUp(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     - this->EvaluateInteractionCoefficientUpUp(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     - this->EvaluateInteractionCoefficientUpUp(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
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

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientUpUp(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										  this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     + this->EvaluateInteractionCoefficientUpUp(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     - this->EvaluateInteractionCoefficientUpUp(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     - this->EvaluateInteractionCoefficientUpUp(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactorsupup[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsupup[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // downdown-downdown
      this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientDownDown(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										  this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     + this->EvaluateInteractionCoefficientDownDown(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     - this->EvaluateInteractionCoefficientDownDown(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     - this->EvaluateInteractionCoefficientDownDown(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
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

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientDownDown(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										  this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     + this->EvaluateInteractionCoefficientDownDown(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     - this->EvaluateInteractionCoefficientDownDown(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     - this->EvaluateInteractionCoefficientDownDown(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactorsdowndown[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsdowndown[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // updown-updown
      this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (
					      - this->EvaluateInteractionCoefficientUpDown(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
										    this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxUp)
					     - this->EvaluateInteractionCoefficientUpDown(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
										    this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown));
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

		  double TmpCoefficient   = (- this->EvaluateInteractionCoefficientUpDown(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
										    this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxUp)
					     - this->EvaluateInteractionCoefficientUpDown(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
										    this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactorsupdown[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsupdown[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      // upup-upup 
      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientUpUp(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										  this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     + this->EvaluateInteractionCoefficientUpUp(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     + this->EvaluateInteractionCoefficientUpUp(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     + this->EvaluateInteractionCoefficientUpUp(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
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

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientUpUp(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										  this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     + this->EvaluateInteractionCoefficientUpUp(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     + this->EvaluateInteractionCoefficientUpUp(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
					     + this->EvaluateInteractionCoefficientUpUp(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp, 
										    this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
		  if (m1 == m2)		    
		    TmpCoefficient *= 0.5;
		  if (m3 == m4)		    
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactorsupup[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsupup[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // downdown-downdown
      this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientDownDown(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										  this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     + this->EvaluateInteractionCoefficientDownDown(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     + this->EvaluateInteractionCoefficientDownDown(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     + this->EvaluateInteractionCoefficientDownDown(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
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

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientDownDown(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										  this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     + this->EvaluateInteractionCoefficientDownDown(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     + this->EvaluateInteractionCoefficientDownDown(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
					     + this->EvaluateInteractionCoefficientDownDown(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown, 
										    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
		  if (m1 == m2)		    
		    TmpCoefficient *= 0.5;
		  if (m3 == m4)		    
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactorsdowndown[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsdowndown[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // updown-updown
      this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientUpDown(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
										    this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxUp)
					     + this->EvaluateInteractionCoefficientUpDown(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
										    this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown));
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

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientUpDown(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
										    this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxUp)
					     + this->EvaluateInteractionCoefficientUpDown(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown, 
										    this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxDown));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactorsupdown[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsupdown[i][Index] = 0.0;
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

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 upup coupling term
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

double ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian::EvaluateInteractionCoefficientUpUp(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
													double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4)
{
  return this->EvaluateInteractionCoefficient (m1, m2, m3, m4, nbrPseudopotentials, pseudopotentials, spinFluxM1, spinFluxM2, spinFluxM3, spinFluxM4);
}



// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 downdown coupling term
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

double ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian::EvaluateInteractionCoefficientDownDown(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
													double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4)
{
  m1 = (this->NbrLzValue - m1) % this->NbrLzValue;
  m2 = (this->NbrLzValue - m2) % this->NbrLzValue;
  m3 = (this->NbrLzValue - m3) % this->NbrLzValue;
  m4 = (this->NbrLzValue - m4) % this->NbrLzValue;
  return this->EvaluateInteractionCoefficientUpUp (m1, m2, m3, m4, nbrPseudopotentials, pseudopotentials, spinFluxM1, spinFluxM2, spinFluxM3, spinFluxM4);
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 updown coupling term
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

double ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian::EvaluateInteractionCoefficientUpDown(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
													double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4)
{
//   m4 = (this->NbrLzValue - m4) % this->NbrLzValue;
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->NbrLzValue);
  double N2 = ((double) (m1 - m4)) + spinFluxM1 - spinFluxM4;
  double Factor =  - (((double) (m1-m2)) + spinFluxM1 - spinFluxM3) * PIOnM * 2.0;
  double Sum = 0.0;
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
//   cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << (Sum / ((double) this->NbrLzValue)) << endl;
  //Normalize per flux (gives correct energy scale for 2-particle problem)
  return (Sum / ((double) this->NbrLzValue));
  
}


// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->FastMultiplicationFlag = false;
  this->Particles = (ParticleOnTorusWithSpinAndMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}


// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}


// get all the indices that should appear in the annihilation/creation operators
//

void ParticleOnTorusWithSpinAndMagneticTranslationsTimeReversalSymmetricGenericHamiltonian::GetIndices()
{
  this->NbrInterSectorSums = this->NbrLzValue;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;      

      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = 0; m2 <= this->LzMax; ++m2)
	  ++this->NbrInterSectorIndicesPerSum[(m1 - m2 + this->NbrLzValue) % this->NbrLzValue];
      this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  if (this->NbrInterSectorIndicesPerSum[i]  > 0)
	    {
	      this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
	      this->NbrInterSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = 0; m2 <= this->LzMax; ++m2)
	  {
	    int TmpSum = (m1 - m2 + this->NbrLzValue) % this->NbrLzValue;
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = m1;
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = m2;
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	}


  this->NbrIntraSectorSums = this->NbrLzValue;
  this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    this->NbrIntraSectorIndicesPerSum[i] = 0;      
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    int TmpSum = (m1 + m2) % this->NbrLzValue;
	    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = m1;
	    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
	}
    }
  else
    {
      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  {
	    int TmpSum = (m1 + m2) % this->NbrLzValue;
	    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = m1;
	    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
	  }
    }
}
