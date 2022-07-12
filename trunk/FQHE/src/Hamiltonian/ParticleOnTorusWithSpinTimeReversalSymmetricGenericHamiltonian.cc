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
//                   two body pseudopotential interaction                     //
//                                                                            //
//                        last modification : 11/02/2015                      //
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

#include "Hamiltonian/ParticleOnTorusWithSpinTimeReversalSymmetricGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnTorusWithSpinGenericHamiltonian.h"
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
// nbrPseudopotentialsUpUp = number of pseudopotentials for up-up interaction
// pseudopotentialsUpUp = pseudopotential coefficients for up-up interaction
// nbrPseudopotentialsDownDown = number of pseudopotentials for down-down interaction
// pseudopotentialsDownDown = pseudopotential coefficients for down-down interaction
// nbrPseudopotentialsUpDown = number of pseudopotentials for up-down interaction
// pseudopotentialsUpDown = pseudopotential coefficients for up-down interaction
// pseudopotentials = pseudopotential coefficients
// spinFluxUp = additional inserted flux for spin up
// spinFluxDown = additional inserted flux for spin down
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusWithSpinTimeReversalSymmetricGenericHamiltonian::ParticleOnTorusWithSpinTimeReversalSymmetricGenericHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int maxMomentum, double ratio, 
										     int nbrPseudopotentialsUpUp, double* pseudopotentialsUpUp,
										     int nbrPseudopotentialsDownDown, double* pseudopotentialsDownDown,
										     int nbrPseudopotentialsUpDown, double* pseudopotentialsUpDown,
										     double spinFluxUp, double spinFluxDown,
										     AbstractArchitecture* architecture, long memory, char* precalculationFileName, double * oneBodyPotentielUpUp, double * oneBodyPotentielDownDown, double * oneBodyPotentielUpDown )
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
  this->SpinFluxUp = spinFluxUp;
  this->SpinFluxDown = spinFluxDown;

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

  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->OneBodyInteractionFactorsupup = 0;
  if(oneBodyPotentielUpUp != 0)
    {
      cout << "One-body interaction factors up-up : " << endl;
      this->OneBodyInteractionFactorsupup = new double[this->NbrLzValue];
      for(int i = 0; i < this->NbrLzValue; i++)
	{
	  this->OneBodyInteractionFactorsupup[i] = oneBodyPotentielUpUp[i];
	  cout << this->OneBodyInteractionFactorsupup[i] << " ";
	}
      cout << endl;
    }
  this->OneBodyInteractionFactorsdowndown = 0;
  if(oneBodyPotentielDownDown != 0)
    {
      cout << "One-body interaction factors down-down : " << endl;
      this->OneBodyInteractionFactorsdowndown = new double[this->NbrLzValue];
      for(int i = 0; i < this->NbrLzValue; i++)
	{
	  this->OneBodyInteractionFactorsdowndown[i] = oneBodyPotentielDownDown[i];
	  cout << this->OneBodyInteractionFactorsdowndown[i] << " ";
	}
      cout << endl;
    } 
  this->OneBodyInteractionFactorsupdown = 0;
  if(oneBodyPotentielUpDown != 0)
    {
      cout << "One-body interaction factors up-down : " << endl;
      this->OneBodyInteractionFactorsupdown = new double[this->NbrLzValue];
      for(int i = 0; i < this->NbrLzValue; i++)
	{
	  this->OneBodyInteractionFactorsupdown[i] = oneBodyPotentielUpDown[i];
	  cout << this->OneBodyInteractionFactorsupdown[i] << " ";
	}
      cout <<endl;
    } 

 
  this->EvaluateInteractionFactors();
  this->S2Hamiltonian = 0;
  this->L2Hamiltonian = 0;
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
/*  
  
  for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
    for (int m2 = 0; m2 < this->NbrLzValue; ++m2)
      for (int m3 = 0; m3 < this->NbrLzValue; ++m3)
	for (int m4 = 0; m4 < this->NbrLzValue; ++m4)
	  cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << this->EvaluateInteractionCoefficientUpDown(m1, m2, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
												 this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp) << endl;

  cout << endl;
  for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
    for (int m2 = 0; m2 < this->NbrLzValue; ++m2)
      for (int m3 = 0; m3 < this->NbrLzValue; ++m3)
	for (int m4 = 0; m4 < this->NbrLzValue; ++m4)
	  cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << this->EvaluateInteractionCoefficientUpUp(m1, m2, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
												 this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp) << endl;*/
}

// destructor
//

ParticleOnTorusWithSpinTimeReversalSymmetricGenericHamiltonian::~ParticleOnTorusWithSpinTimeReversalSymmetricGenericHamiltonian() 
{
}


  
// evaluate all interaction factors
//   

void ParticleOnTorusWithSpinTimeReversalSymmetricGenericHamiltonian::EvaluateInteractionFactors()
{  
  /*
  this->QxValues = new double [this->NBodyValue];
  this->QyValues = new double [this->NBodyValue];
  this->Q2Values = new double [this->NBodyValue];
  this->CosineCoffients = new double [this->NBodyValue];*/
  
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
		  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficientUpUp(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
												 this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
							    + this->EvaluateInteractionCoefficientUpUp(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
												   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
							    - this->EvaluateInteractionCoefficientUpUp(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
												   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
							    - this->EvaluateInteractionCoefficientUpUp(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
												   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
		  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficientDownDown(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
												     this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
								+ this->EvaluateInteractionCoefficientDownDown(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
												       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
								- this->EvaluateInteractionCoefficientDownDown(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
												       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
								- this->EvaluateInteractionCoefficientDownDown(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
												       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
		  //		  this->InteractionFactorsupup[i][Index] = 0.0;
		  //		  this->InteractionFactorsdowndown[i][Index] = 0.0;		  
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
		  this->InteractionFactorsupdown[i][Index] = (-this->EvaluateInteractionCoefficientUpDown(m1, m2, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
												    this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxDown)
							      -this->EvaluateInteractionCoefficientUpDown(m2, m1, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
												    this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxUp));
          
//           this->InteractionFactorsupdown[i][Index] = -this->EvaluateInteractionCoefficientUpUp(m1, ((this->NbrLzValue-m2) % (this->NbrLzValue)), ((this->NbrLzValue-m4) % (this->NbrLzValue)), m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
// 												    this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxDown);
//           this->InteractionFactorsupdown[i][Index] -= this->EvaluateInteractionCoefficientUpUp(m2, ((this->NbrLzValue-m1) % (this->NbrLzValue)), ((this->NbrLzValue-m3) % (this->NbrLzValue)), m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
// 												    this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxDown);
          
// 		  cout << m1 << " " << m2 << " " << m4 << " " << m3 <<  " " << this->InteractionFactorsupdown[i][Index] << endl;
// 		  cout << m2 << " " << m1 << " " << m3 << " " << m4 << " " << this->InteractionFactorsupdown[i][Index] << endl;
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
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficientUpUp(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													 this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
								    + this->EvaluateInteractionCoefficientUpUp(m2, m1, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
								    + this->EvaluateInteractionCoefficientUpUp(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
								    + this->EvaluateInteractionCoefficientUpUp(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
			  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficientDownDown(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													     this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
									+ this->EvaluateInteractionCoefficientDownDown(m2, m1, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
									+ this->EvaluateInteractionCoefficientDownDown(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
									+ this->EvaluateInteractionCoefficientDownDown(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
			}
		      else
			{
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficientUpUp(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													 this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
								    + this->EvaluateInteractionCoefficientUpUp(m2, m1, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
			  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficientDownDown(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													     this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
									+ this->EvaluateInteractionCoefficientDownDown(m2, m1, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficientUpUp(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													 this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp)
								    + this->EvaluateInteractionCoefficientUpUp(m1, m2, m4, m3, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													   this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp));
			  this->InteractionFactorsdowndown[i][Index] = (this->EvaluateInteractionCoefficientDownDown(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													     this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown)
									+ this->EvaluateInteractionCoefficientDownDown(m1, m2, m4, m3, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													       this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown));
			}
		      else
			{
			  this->InteractionFactorsupup[i][Index] = this->EvaluateInteractionCoefficientUpUp(m1, m2, m3, m4, this->NbrPseudopotentialsUpUp, this->PseudopotentialsUpUp,
													this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp, this->SpinFluxUp);
			  this->InteractionFactorsdowndown[i][Index] = this->EvaluateInteractionCoefficientDownDown(m1, m2, m3, m4, this->NbrPseudopotentialsDownDown, this->PseudopotentialsDownDown,
													    this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown, this->SpinFluxDown);
			}
		    }
// 		  this->InteractionFactorsupup[i][Index] *= -1.0;
// 		  this->InteractionFactorsdowndown[i][Index] *= -1.0;		  
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
		  this->InteractionFactorsupdown[i][Index] = this->EvaluateInteractionCoefficientUpDown(m1, m2, m3, m4, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
													   this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxDown) 
		  + this->EvaluateInteractionCoefficientUpDown(m2, m1, m4, m3, this->NbrPseudopotentialsUpDown, this->PseudopotentialsUpDown,
													   this->SpinFluxUp, this->SpinFluxDown, this->SpinFluxUp, this->SpinFluxDown) ;
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }

  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
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

double ParticleOnTorusWithSpinTimeReversalSymmetricGenericHamiltonian::EvaluateInteractionCoefficientUpUp(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
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

double ParticleOnTorusWithSpinTimeReversalSymmetricGenericHamiltonian::EvaluateInteractionCoefficientDownDown(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
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

double ParticleOnTorusWithSpinTimeReversalSymmetricGenericHamiltonian::EvaluateInteractionCoefficientUpDown(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
										 double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->NbrLzValue);
  double Factor =  - (((double) (m1 + m2)) + spinFluxM1 + spinFluxM4) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = ((double) (m1 - m4)) + spinFluxM1 - spinFluxM4;
  double N1;
  double Q2x;
  double Q2y;
  double Q2;
  double Precision;
  double TmpInteraction;
  int count = 0;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2x = 0.0;
      Q2y = this->Ratio * N2 * N2;
      Q2 = Q2x + Q2y;
      if (N2 != 0.0)
	{
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	    {
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0* PIOnM * Q2);
// 	      TmpInteraction += pseudopotentials[i] * Q2 * exp(-4*2*M_PI*Q2);	
// 	      cout << i << " " << m1 << " " << m2 << " " << m3 << " " << m4 << " " << this->LaguerrePolynomials[i].PolynomialEvaluate(2.0* PIOnM * Q2) << endl;
	    }
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
// // 	  cout << "0" << " " << m1 << " " << m2 << " " << m3 << " " << m4 << " " << Coefficient << " " << Sum << endl;
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
	    {
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
// 	      TmpInteraction += pseudopotentials[i] * Q2* exp(-4*2*M_PI*Q2);
	    }
	  Coefficient = TmpInteraction;
// 	  cout << "0" << " " << m1 << " " << m2 << " " << m3 << " " << m4 << " " << Coefficient << " " << Sum << endl;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
// 	  cout << "0, " << N << " : " << m1 << " " << m2 << " " << m3 << " " << m4 << " " << Coefficient << " " << Sum << endl;	  
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	    {
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
// 	      TmpInteraction += pseudopotentials[i] * Q2 * exp(-4*2*M_PI*Q2);
	    }
	  Precision = 2.0 * exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
// 	  cout << "1" << " " << m1 << " " << m2 << " " << m3 << " " << m4 << " " << Coefficient << " " << Sum << endl;
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
	    {
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
// 	      TmpInteraction += pseudopotentials[i] * Q2 * exp(-4*2*M_PI*Q2);
	    }
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
// 	  cout << "2" << " " << m1 << " " << m2 << " " << m3 << " " << m4 << " " << Coefficient << " " << Sum << endl;
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
	    {
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
// 	      TmpInteraction += pseudopotentials[i] * Q2 * exp(-4*2*M_PI*Q2);
	    }
	  Coefficient = TmpInteraction;
// 	  cout << "2" << " " << m1 << " " << m2 << " " << m3 << " " << m4 << " " << Coefficient << " " << Sum << endl;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < nbrPseudopotentials; ++i)
	    if (pseudopotentials[i] != 0.0)
	    {
	      TmpInteraction += pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(2.0 * PIOnM * Q2);
// 	       TmpInteraction += pseudopotentials[i] * Q2 * exp(-4*2*M_PI*Q2);
	    }
	  Precision = 2.0 *  exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
// 	  cout << "3" << " " << m1 << " " << m2 << " " << m3 << " " << m4 << " " << Coefficient << " " << Sum << endl;
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->NbrLzValue;
    }
//   cout << "final" << " " << m1 << " " << m2 << " " << m3 << " " << m4 << " " << (Sum / ((double) this->NbrLzValue)) << endl;
  //Normalize per flux (gives correct energy scale for 2-particle problem)
  return (Sum / ((double) this->NbrLzValue));
}



// evaluate the numerical coefficient  in front of the Prod a^+_mi Prod a+_n coupling term
//
// mIndices = array containing the creation operator indices
// nIndices = array containing the annihilation operator indices
// return value = numerical coefficient  

// double ParticleOnTorusWithSpinTimeReversalSymmetricGenericHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrPseudopotentials, double* pseudopotentials,
// 										 double spinFluxM1, double spinFluxM2, double spinFluxM3, double spinFluxM4)
// {
//   int Tmp;
//   int* mIndices = new int[2];
//   int* nIndices = new int[2];
//   mIndices[0] = m1;
//   mIndices[1] = m2;
//   nIndices[0] = m3;
//   nIndices[1] = m4;
//   double Prefactor = powl(this->MaxMomentum, 1.0);
//   for (int i = 0; i < 2; ++i)
//     {
//       this->QxValues[i] = 0.0;
//       this->QyValues[i] = (double) (nIndices[i] - mIndices[i]);
//       this->CosineCoffients[i] = 2.0 * ((double) mIndices[i]);
//     }  
//   double CurrentPrecision;
//   double Coefficient = Prefactor * this->RecursiveEvaluateInteractionCoefficient(0, 0.0, 0.0, 0.0, 0.0, CurrentPrecision, nbrOperations);
//   return Coefficient;
// }
  

// double ParticleOnTorusWithSpinTimeReversalSymmetricGenericHamiltonian::RecursiveEvaluateInteractionCoefficient(int xPosition, double currentSumQx, double currentSumQy, double currentSumQ2, double currentSumPhase, double& currentPrecision, long& nbrOperations)
// {
//   if (xPosition < 1)
//     {
//       double TotalCoefficient  = 0.0;
//       double Coefficient  = 1.0;
//       int CurrentQy = this->QyValues[xPosition];
//       currentPrecision = 0.0;
//       double TmpPrecision;
//       double Tmp;
//       while ((fabs(Coefficient) + fabs(TotalCoefficient)) != fabs(TotalCoefficient))
// 	{	        
// 	  Tmp = (this->QyValues[xPosition] * this->QyValues[xPosition] * this->Ratio);
// 	  this->QxValues[xPosition] = 0.0;
// 	  this->Q2Values[xPosition] = Tmp;
// 	  Coefficient = this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
// 								      currentSumQy + this->QyValues[xPosition], 
// 								      currentSumQ2 + this->Q2Values[xPosition], 
// 								      currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
// 	  currentPrecision += TmpPrecision;
// 	  if (Coefficient == 0.0)
// 	    TmpPrecision = 1.0;
// 	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
// 	    {	  
// 	      ++this->QxValues[xPosition];
// 	      this->Q2Values[xPosition] = (this->QxValues[xPosition] * this->QxValues[xPosition] * this->InvRatio) + Tmp;
// 	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
// 									   currentSumQy + this->QyValues[xPosition], 
// 									   currentSumQ2 + this->Q2Values[xPosition], 
// 									   currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
// 	      currentPrecision += TmpPrecision;
// 	    }
// 	  this->QxValues[xPosition] = 0.0;
// 	  if (Coefficient == 0.0)
// 	    TmpPrecision = 1.0;
// 	  else
// 	    TmpPrecision = 2.0 * fabs(Coefficient);
// 	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
// 	    {	  
// 	      --this->QxValues[xPosition];
// 	      this->Q2Values[xPosition] = (this->QxValues[xPosition] * this->QxValues[xPosition] * this->InvRatio) + Tmp;
// 	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
// 									   currentSumQy + this->QyValues[xPosition], 
// 									   currentSumQ2 + this->Q2Values[xPosition], 
// 									   currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
// 	      currentPrecision += TmpPrecision;
// 	    }
// 	  this->QyValues[xPosition] +=  (double) this->MaxMomentum;
// 	  TotalCoefficient += Coefficient;
// 	}
//       this->QyValues[xPosition] = CurrentQy -  (double) this->MaxMomentum;
//       if (TotalCoefficient == 0.0)
// 	Coefficient = 1.0;
//       else
// 	Coefficient = 2.0 * TotalCoefficient;
//       while ((fabs(Coefficient) + fabs(TotalCoefficient)) != fabs(TotalCoefficient))
// 	{	        
// 	  this->QxValues[xPosition] = 0.0;
// 	  Tmp = (this->QyValues[xPosition] * this->QyValues[xPosition] * this->Ratio);
// 	  this->Q2Values[xPosition] = Tmp;
// 	  Coefficient = this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
// 								      currentSumQy + this->QyValues[xPosition], 
// 								      currentSumQ2 + this->Q2Values[xPosition], 
// 								      currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
// 	  currentPrecision += TmpPrecision;
// 	  if (Coefficient == 0.0)
// 	    TmpPrecision = 1.0;
// 	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
// 	    {	  
// 	      ++this->QxValues[xPosition];
// 	      this->Q2Values[xPosition] = (this->QxValues[xPosition] * this->QxValues[xPosition] * this->InvRatio) + Tmp;
// 	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
// 									   currentSumQy + this->QyValues[xPosition], 
// 									   currentSumQ2 + this->Q2Values[xPosition], 
// 									   currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
// 	      currentPrecision += TmpPrecision;
// 	    }
// 	  this->QxValues[xPosition] = 0.0;
// 	  if (Coefficient == 0.0)
// 	    TmpPrecision = 1.0;
// 	  else
// 	    TmpPrecision = 2.0 * fabs(Coefficient);
// 	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
// 	    {	  
// 	      --this->QxValues[xPosition];
// 	      this->Q2Values[xPosition] = (this->QxValues[xPosition] * this->QxValues[xPosition] * this->InvRatio) + Tmp;
// 	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
// 									   currentSumQy + this->QyValues[xPosition], 
// 									   currentSumQ2 + this->Q2Values[xPosition], 
// 									   currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
// 	      currentPrecision += TmpPrecision;
// 	    }
// 	  this->QyValues[xPosition] -=  (double) this->MaxMomentum;
// 	  TotalCoefficient += Coefficient;
// 	}      
//       this->QyValues[xPosition] = CurrentQy;
//       return TotalCoefficient;
//     }
//   else
//     {
//       double TmpExponentialFactor = M_PI / ((double) this->MaxMomentum);
//       this->QxValues[xPosition] = -currentSumQx;
//       this->QyValues[xPosition] = -currentSumQy;  
//       this->Q2Values[xPosition] = (this->QxValues[xPosition] * this->QxValues[xPosition] * this->InvRatio) + (this->QyValues[xPosition] * this->QyValues[xPosition] * this->Ratio);
//       currentSumPhase += this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]);
//       currentPrecision = exp(- 0.5 * TmpExponentialFactor * (this->Q2Values[xPosition] + currentSumQ2)) * this->VFactor(this->Q2Values);
//       ++nbrOperations;
//       return (cos(TmpExponentialFactor * currentSumPhase) * currentPrecision);
// // if we do not assume that it is invariant under {qx}<->{-qx}
// //      return (Complex(TmpExponentialFactor * Sum2) * currentPrecision);
//     }
// }
