////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                         generic two body interaction                       //
//                                                                            //
//                        last modification : 17/04/2012                      //
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


#include "Hamiltonian/ParticleOnTorusWithSU3SpinGenericHamiltonian.h"
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
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusWithSU3SpinGenericHamiltonian::ParticleOnTorusWithSU3SpinGenericHamiltonian(ParticleOnSphereWithSU3Spin* particles, int nbrParticles, int maxMomentum, double ratio, 
											   int* nbrPseudopotentials, double** pseudoPotentials,
											   double spinFlux1, double spinFlux2, double spinFlux3,
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

  this->LaguerrePolynomials =new Polynomial[this->MaxNbrPseudopotentials];
  for (int i = 0; i < this->MaxNbrPseudopotentials; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  this->OneBodyInteractionFactors11 = 0;
  this->OneBodyInteractionFactors22 = 0;
  this->OneBodyInteractionFactors33 = 0;

  this->HamiltonianShift = 0.0;
  this->DiskStorageFlag = false;
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

ParticleOnTorusWithSU3SpinGenericHamiltonian::~ParticleOnTorusWithSU3SpinGenericHamiltonian() 
{
  for (int i = 0; i < 6; ++i)
    {
      if (this->NbrPseudopotentials[i] > 0)
	{
	  delete[] this->Pseudopotentials[i];
	}
    } 
  delete[] this->NbrPseudopotentials;
  delete[] this->Pseudopotentials;
  delete[] this->LaguerrePolynomials;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusWithSU3SpinGenericHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// evaluate all interaction factors
//   

void ParticleOnTorusWithSU3SpinGenericHamiltonian::EvaluateInteractionFactors()
{
  this->M1IntraValue = 0;
  this->M1InterValue = 0;
  
  long TotalNbrInteractionFactors = 0;

  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;

  this->NbrInterSectorSums = this->NbrLzValue;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      ++this->NbrInterSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];      
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
      this->NbrInterSectorIndicesPerSum[i] = 0;
    }
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      {
	int TmpIndex = (m1 + m2) % this->NbrLzValue;
	this->InterSectorIndicesPerSum[TmpIndex][this->NbrInterSectorIndicesPerSum[TmpIndex] << 1] = m1;
	this->InterSectorIndicesPerSum[TmpIndex][1 + (this->NbrInterSectorIndicesPerSum[TmpIndex] << 1)] = m2;
	++this->NbrInterSectorIndicesPerSum[TmpIndex];
      }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NbrIntraSectorSums = this->NbrLzValue;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    int TmpIndex = (m1 + m2) % this->NbrLzValue;
	    this->IntraSectorIndicesPerSum[TmpIndex][this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1] = m1;
	    this->IntraSectorIndicesPerSum[TmpIndex][1 + (this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[TmpIndex];
	  }

      this->InteractionFactors11 = new double* [this->NbrIntraSectorSums];
      this->InteractionFactors22 = new double* [this->NbrIntraSectorSums];
      this->InteractionFactors33 = new double* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors11[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactors22[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactors33[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  this->InteractionFactors11[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
											       this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
							    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
												   this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
							    - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
												   this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
							    - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
												   this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1));
		  this->InteractionFactors22[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
											       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
								+ this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
												       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
								- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
												       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
								- this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
												       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2));
		  this->InteractionFactors33[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
											       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
								+ this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
												       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
								- this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
												       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
								- this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
												       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3));
		  //		  this->InteractionFactors11[i][Index] = 0.0;
		  //		  this->InteractionFactors22[i][Index] = 0.0;		  
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      this->InteractionFactors12 = new double* [this->NbrInterSectorSums];
      this->InteractionFactors13 = new double* [this->NbrInterSectorSums];
      this->InteractionFactors23 = new double* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors12[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactors13[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactors23[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  this->InteractionFactors12[i][Index] = (-this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[1], this->Pseudopotentials[1],
												this->SpinFlux1, this->SpinFlux2, this->SpinFlux1, this->SpinFlux2)
							  -this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[1], this->Pseudopotentials[1],
												this->SpinFlux1, this->SpinFlux2, this->SpinFlux1, this->SpinFlux2));
		  this->InteractionFactors13[i][Index] = (-this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[2], this->Pseudopotentials[2],
												this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
							  -this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[2], this->Pseudopotentials[2],
												this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2));
		  this->InteractionFactors23[i][Index] = (-this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[4], this->Pseudopotentials[4],
												this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
							  -this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[4], this->Pseudopotentials[4],
												this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3));
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      this->NbrIntraSectorSums = this->NbrLzValue;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i] != 0)
	    this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  else
	    this->IntraSectorIndicesPerSum[i] = 0;
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  {
	    int TmpIndex = (m1 + m2) % this->NbrLzValue;
	    this->IntraSectorIndicesPerSum[TmpIndex][this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1] = m1;
	    this->IntraSectorIndicesPerSum[TmpIndex][1 + (this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[TmpIndex];
	  }

      this->InteractionFactors11 = new double* [this->NbrIntraSectorSums];
      this->InteractionFactors22 = new double* [this->NbrIntraSectorSums];
      this->InteractionFactors33 = new double* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors11[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactors22[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactors33[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
			  this->InteractionFactors11[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
												       this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
								    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
													   this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
													   this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
													   this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1));
			  this->InteractionFactors22[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
												       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
									+ this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
													       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
									+ this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
													       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
									+ this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
													       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2));
			  this->InteractionFactors33[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
												       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
									+ this->EvaluateInteractionCoefficient(m2, m1, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
													       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
									+ this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
													       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
									+ this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
													       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3));
			}
		      else
			{
			  this->InteractionFactors11[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
												       this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
													   this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1));
			  this->InteractionFactors22[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
												       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
									+ this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
													       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2));
			  this->InteractionFactors33[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
												       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
									+ this->EvaluateInteractionCoefficient(m2, m1, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
													       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactors11[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
												       this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1)
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
													   this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1));
			  this->InteractionFactors22[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
												       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2)
									+ this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
													       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2));
			  this->InteractionFactors33[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
												       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3)
									+ this->EvaluateInteractionCoefficient(m1, m2, m4, m3, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
													       this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3));
			}
		      else
			{
			  this->InteractionFactors11[i][Index] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[0], this->Pseudopotentials[0],
												      this->SpinFlux1, this->SpinFlux1, this->SpinFlux1, this->SpinFlux1);
			  this->InteractionFactors22[i][Index] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[3], this->Pseudopotentials[3],
												       this->SpinFlux2, this->SpinFlux2, this->SpinFlux2, this->SpinFlux2);
			  this->InteractionFactors33[i][Index] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[5], this->Pseudopotentials[5],
												      this->SpinFlux3, this->SpinFlux3, this->SpinFlux3, this->SpinFlux3);
			}
		    }
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      this->InteractionFactors12 = new double* [this->NbrInterSectorSums];
      this->InteractionFactors13 = new double* [this->NbrInterSectorSums];
      this->InteractionFactors23 = new double* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors12[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactors13[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactors23[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      double Factor = 2.0;
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

// 		  ParticleOnSphereWithSU3Spin* TmpHilbertTruc = (ParticleOnSphereWithSU3Spin*) this->Particles;
		    
// 		  cout << m1 << " " << m2 << " " << m3 << " " << m4 << " : " << endl;;
// 		  for (int truc = 0; truc < TmpHilbertTruc->GetHilbertSpaceDimension(); ++truc)
// 		    {
// 		      double TrucCoefficient = TmpHilbertTruc->A2A3(truc, m3, m4);
// 		      cout << truc << " : " << TrucCoefficient;
// 		      if (TrucCoefficient != 0.0)
// 			{
// 			  int TrucIndex = TmpHilbertTruc->Ad2Ad3(m1, m2, TrucCoefficient);
// 			  cout << " : " << TrucIndex << " " << TrucCoefficient;
// 			}
// 		      cout << endl;
// 		    }

		  this->InteractionFactors12[i][Index] = Factor * this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[1], this->Pseudopotentials[1],
												       this->SpinFlux1, this->SpinFlux2, this->SpinFlux1, this->SpinFlux2);
		  this->InteractionFactors13[i][Index] = Factor * this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[2], this->Pseudopotentials[2],
												       this->SpinFlux1, this->SpinFlux3, this->SpinFlux1, this->SpinFlux3);
		  this->InteractionFactors23[i][Index] = Factor * this->EvaluateInteractionCoefficient(m1, m2, m3, m4, this->NbrPseudopotentials[4], this->Pseudopotentials[4],
												       this->SpinFlux2, this->SpinFlux3, this->SpinFlux2, this->SpinFlux3);
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

double ParticleOnTorusWithSU3SpinGenericHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
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



