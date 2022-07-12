////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of hamiltonian associated to particles on a cylinder with       //
//     SU(2) spin and a generic interaction defined by its pseudopotential    //
//                                                                            //
//                        last modification : 19/08/2016                      //
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


#include "Hamiltonian/ParticleOnCylinderWithSpinGenericHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"
#include "Polynomial/SpecialPolynomial.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

ParticleOnCylinderWithSpinGenericHamiltonian::ParticleOnCylinderWithSpinGenericHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// ratio = ratio between the width in the x direction and the width in the y direction
// architecture = architecture to use for precalculation
// nbrPseudopotentialsUpUp = number of pseudopotentials for up-up interaction
// pseudopotentialsUpUp = pseudopotential coefficients for up-up interaction
// nbrPseudopotentialsDownDown = number of pseudopotentials for down-down interaction
// pseudopotentialsDownDown = pseudopotential coefficients for down-down interaction
// nbrPseudopotentialsUpDown = number of pseudopotentials for up-down interaction
// pseudopotentialsUpDown = pseudopotential coefficients for up-down interaction
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialUpDown =  one-body tunnelling potential (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down

ParticleOnCylinderWithSpinGenericHamiltonian::ParticleOnCylinderWithSpinGenericHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, double ratio, 
											   int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
											   int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
											   int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
											   double spinFluxUp, double spinFluxDown,
											   AbstractArchitecture* architecture, long memory, char* precalculationFileName, 
											   double * oneBodyPotentialUpUp, double * oneBodyPotentialDownDown, double * oneBodyPotentialUpDown)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / this->Ratio;
  this->SpinFluxUp = 0.0;
  this->SpinFluxDown = 0.0;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = false;
  this->Architecture = architecture;
  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;

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
  this->LaguerrePolynomials =new Polynomial[this->MaxNbrPseudopotentials];
  for (int i = 0; i < this->MaxNbrPseudopotentials; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = false;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  if (oneBodyPotentialUpUp != 0)
    {
      this->OneBodyInteractionFactorsupup = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsupup[i] = oneBodyPotentialUpUp[i];
    }
  this->OneBodyInteractionFactorsdowndown = 0;
  if (oneBodyPotentialDownDown != 0)
    {
      this->OneBodyInteractionFactorsdowndown = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsdowndown[i] = oneBodyPotentialDownDown[i];
    }
  this->OneBodyInteractionFactorsupdown = 0;
  if (oneBodyPotentialUpDown != 0)
    {
      this->OneBodyInteractionFactorsupdown = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsupdown[i] = oneBodyPotentialUpDown[i];
    }
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
	    {
	      cout  << "fast = " << (TmpMemory >> 30) << ".";
	      TmpMemory -= ((TmpMemory >> 30) << 30);
	      TmpMemory *= 100l;
	      TmpMemory >>= 30;
	      if (TmpMemory < 10l)
		cout << "0";
	      cout  << TmpMemory << " Gb ";
	    }
	  if (this->DiskStorageFlag == false)
	    {
	      this->EnableFastMultiplication();
	    }
	  else
	    {
	      char* TmpFileName = this->Architecture->GetTemporaryFileName();
	      this->EnableFastMultiplicationWithDiskStorage(TmpFileName);	      
	      delete[] TmpFileName;
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnCylinderWithSpinGenericHamiltonian::~ParticleOnCylinderWithSpinGenericHamiltonian() 
{
}

// evaluate all interaction factors
//   

void ParticleOnCylinderWithSpinGenericHamiltonian::EvaluateInteractionFactors()
{
  // in absence of above code, initialize the following two fields to zero
  this->M1IntraValue = 0;
  this->M1InterValue = 0;
  
  // int Lim;
  // int Min;
  // int Pos = 0;
  // int m4;
  long TotalNbrInteractionFactors = 0;

  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  double TmpCoefficient = 0.0;

  this->NbrInterSectorSums = 2 * this->LzMax + 1;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      ++this->NbrInterSectorIndicesPerSum[m1 + m2];      
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
      this->NbrInterSectorIndicesPerSum[i] = 0;
    }
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      {
	this->InterSectorIndicesPerSum[(m1 + m2)][this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
	this->InterSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
	++this->NbrInterSectorIndicesPerSum[(m1 + m2)];
      }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NbrIntraSectorSums = 2 * this->LzMax - 1;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    this->IntraSectorIndicesPerSum[(m1 + m2) - 1][this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1] = m1;
	    this->IntraSectorIndicesPerSum[(m1 + m2) - 1][1 + (this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
	  }

      this->InteractionFactorsupup = new double* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new double* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndown[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
												 this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
							    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
												   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
							    - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
												   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
							    - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
												   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
		  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
												     this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
								+ this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
												       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
								- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
												       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
								- this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
												       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      this->InteractionFactorsupdown = new double* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  this->InteractionFactorsupdown[i][Index] = (-this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
												    this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxDown)
							      -this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
												    this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxUp));
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      this->NbrIntraSectorSums = 2 * this->LzMax+1;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[m1 + m2];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  {
	    this->IntraSectorIndicesPerSum[m1 + m2][this->NbrIntraSectorIndicesPerSum[m1 + m2] << 1] = m1;
	    this->IntraSectorIndicesPerSum[m1 + m2][1 + (this->NbrIntraSectorIndicesPerSum[m1 + m2] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[m1 + m2];
	  }

      this->InteractionFactorsupup = new double* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new double* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndown[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  if (m1 != m2)
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													 this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
								    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
			  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													     this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
									+ this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
									+ this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
									+ this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
			}
		      else
			{
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													 this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
			  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													     this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
									+ this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													 this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
			  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													     this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
									+ this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
			}
		      else
			{
			  this->InteractionFactorsupup[i][Index] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp);
			  this->InteractionFactorsdowndown[i][Index] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown);
			}
		    }
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      this->InteractionFactorsupdown = new double* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  this->InteractionFactorsupdown[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
												   this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxDown)
							      + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
												     this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxDown));
		  ++TotalNbrInteractionFactors;
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
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = pseudopotential coefficients
// spinFluxM1 = additional inserted flux for m1
// spinFluxM2 = additional inserted flux for m2
// spinFluxM3 = additional inserted flux for m3
// spinFluxM4 = additional inserted flux for m4
// return value = numerical coefficient

double ParticleOnCylinderWithSpinGenericHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
										    double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4)
{
  double Kappa2Factor = (2.0 * M_PI * this->Ratio / ((double) (this->LzMax + 1)));
  double Coefficient = 0.0;
  if ((nbrPseudopotentials > 0) && (pseudopotentials[0] != 0.0))
    {
      Coefficient += pseudopotentials[0] * exp (-0.25 * Kappa2Factor * (((double) ((m1 - m2) * (m1 - m2))) + ((double) ((m3 - m4) * (m3 - m4))))) / sqrt(2.0 * M_PI);
    }
  if (((nbrPseudopotentials > 1) && pseudopotentials[1] != 0.0))
    {
      Coefficient += pseudopotentials[1] * (((double) ((m1 - m2) * (m4 - m3)))
					    * exp (-0.25 * Kappa2Factor * (((double) ((m1 - m2) * (m1 - m2))) 
									  + ((double) ((m3 - m4) * (m3 - m4)))))) * Kappa2Factor / sqrt(this->Ratio * ((double) (this->LzMax + 1)));
    }
  return Coefficient;
}

