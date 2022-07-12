////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//           class of bosons on the CP2 with delta interaction                //
//                                                                            // 
//                                                                            //
//                        last modification : 24/10/2012                      //
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
#include "Hamiltonian/ParticleOnCP2DeltaHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <cmath>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::sin;
using std::cos;



// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCP2DeltaHamiltonian::ParticleOnCP2DeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->NbrLzValue = (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 2)/2;
  this->LzMax = this->NbrLzValue - 1;
  this->Particles2 = (BosonOnCP2*) this->Particles;
  this->quantumNumberTz = new int [this->NbrLzValue];
  this->quantumNumberY = new int [this->NbrLzValue];
  this->quantumNumberR = new int [this->NbrLzValue];
  this->quantumNumberS = new int [this->NbrLzValue];
  this->OneBodyTermFlag = false;
  this->Particles2->GetQuantumNumbersFromLinearizedIndex(this->quantumNumberTz, this->quantumNumberY, this->quantumNumberR, this->quantumNumberS);
//   for (int i = 0; i<NbrLzValue; i++)
//   {
//    cout << i << " " << this->quantumNumberR[i] << " " << this->quantumNumberS[i] << endl; 
//   }
  this->HamiltonianShift = 0.0;
  
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  cout << "done" << endl;
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      cout << "fast = ";
      PrintMemorySize(cout, TmpMemory)<< endl;
      this->EnableFastMultiplication();
     }
}

// constructor with one body terms
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// oneBodyPotentials = array with the coefficient in front of each one body term (ordered such that the first element corresponds to the one of a+_-s a_-s)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCP2DeltaHamiltonian::ParticleOnCP2DeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta, double*  oneBodyPotentials, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->NbrLzValue = (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 2)/2;
  this->LzMax = this->NbrLzValue - 1;
  this->Particles2 = (BosonOnCP2*) this->Particles;
  this->quantumNumberTz = new int [this->NbrLzValue];
  this->quantumNumberY = new int [this->NbrLzValue];
  this->quantumNumberR = new int [this->NbrLzValue];
  this->quantumNumberS = new int [this->NbrLzValue];
  this->Particles2->GetQuantumNumbersFromLinearizedIndex(this->quantumNumberTz, this->quantumNumberY, this->quantumNumberR, this->quantumNumberS);
//   for (int i = 0; i<NbrLzValue; i++)
//   {
//    cout << i << " " << this->quantumNumberR[i] << " " << this->quantumNumberS[i] << endl; 
//   }
  this->HamiltonianShift = 0.0;
  
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyTermFlag = true;
  this->OneBodyPotentials = new double [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->OneBodyPotentials[i] = oneBodyPotentials[i];
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  cout << "done" << endl;
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      cout << "fast = ";
      PrintMemorySize(cout, TmpMemory)<< endl;
      this->EnableFastMultiplication();
     }
}

// destructor
//

ParticleOnCP2DeltaHamiltonian::~ParticleOnCP2DeltaHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnCP2DeltaHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  FactorialCoefficient Coef;
   
  this->NbrSectorSums = (2*this->NbrFluxQuanta + 1)*(2*this->NbrFluxQuanta + 2)/2;
  this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
  for (int i = 0; i < this->NbrSectorSums; ++i)
  this->NbrSectorIndicesPerSum[i] = 0;  
      
  for (int r1 = 0; r1 <= this->NbrFluxQuanta; ++r1)
    for (int r2 = 0; r2 <= this->NbrFluxQuanta; ++r2)
      for (int s1 = 0; s1 <= this->NbrFluxQuanta - r1; ++s1)
	for (int s2 = 0; s2 <= this->NbrFluxQuanta - r2; ++s2) 
	   {
	      int tz1 = r1 - s1;
	      int y1 = 3*(r1+s1) - 2*this->NbrFluxQuanta;
	      int tz2 = r2 - s2;
	      int y2 = 3*(r2 + s2) - 2*this->NbrFluxQuanta;
	      int Index1 = this->Particles2->GetLinearizedIndex(tz1, y1, 1);
	      int Index2 = this->Particles2->GetLinearizedIndex(tz2, y2, 1);
	      if (Index1 <= Index2)
		++this->NbrSectorIndicesPerSum[this->Particles2->GetLinearizedIndex(tz1 + tz2, y1 + y2 , 2)];    
	    }
		
  this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
  for (int i = 0; i < this->NbrSectorSums; ++i)
    {
      if (this->NbrSectorIndicesPerSum[i]  > 0)
	{
	   this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];      
	   this->NbrSectorIndicesPerSum[i] = 0;
	}
    }
  for (int r1 = 0; r1 <= this->NbrFluxQuanta; ++r1)
    for (int r2 = 0; r2 <= this->NbrFluxQuanta; ++r2)
      for (int s1 = 0; s1 <= this->NbrFluxQuanta - r1; ++s1)
	for (int s2 = 0; s2 <= this->NbrFluxQuanta - r2; ++s2) 
	  {
	      int tz1 = r1 - s1;
	      int y1 = 3*(r1+s1) - 2*this->NbrFluxQuanta;
	      int tz2 = r2 - s2;
	      int y2 = 3*(r2 + s2) - 2*this->NbrFluxQuanta;
	      int Index1 = this->Particles2->GetLinearizedIndex(tz1, y1, 1);
	      int Index2 = this->Particles2->GetLinearizedIndex(tz2, y2, 1);
	    if (Index1 <= Index2)
	      {
		int TmpSum = this->Particles2->GetLinearizedIndex(tz1 + tz2, y1 + y2 , 2);
		this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		++this->NbrSectorIndicesPerSum[TmpSum];    
	      }
	  }
      
      
      this->InteractionFactors = new double* [this->NbrSectorSums];
      //cout << this->NbrSectorSums << endl;
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new double[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int l1 = 0; l1 < this->NbrSectorIndicesPerSum[i]; ++l1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][l1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(l1 << 1) + 1];
	      int r1 = this->quantumNumberR[Index1];
	      int s1 = this->quantumNumberS[Index1];
	      int r2 = this->quantumNumberR[Index2];
	      int s2 = this->quantumNumberS[Index2];
	      for (int l2 = 0; l2 < this->NbrSectorIndicesPerSum[i]; ++l2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][l2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(l2 << 1) + 1];
		  int r3 = this->quantumNumberR[Index3];
		  int s3 = this->quantumNumberS[Index3];
		  int r4 = this->quantumNumberR[Index4];
		  int s4 = this->quantumNumberS[Index4];
		  
// 		  cout << j1 + j2 << "=" << j3 + j4 << endl;
// 		  cout << jz1 + jz2 << "=" << jz3 + jz4 << endl;
// 		  cout << kz1 + kz2 << "=" << kz3 + kz4 << endl;
		  
 		  this->InteractionFactors[i][Index] =  this->ComputeTwoBodyMatrixElement(r1,  s1, r2, s2, r3,  s3, r4, s4, Coef);
// 		  this->InteractionFactors[i][Index] =  0;
		  if (Index3 == Index4)
		    this->InteractionFactors[i][Index] *= 0.5;
		  if (Index1 == Index2)
		    this->InteractionFactors[i][Index] *= 0.5;
		  this->InteractionFactors[i][Index] *= 2.0;

		  TotalNbrInteractionFactors++;
		  ++Index;

		}
	    }
	}
  if (this->OneBodyTermFlag == true)
    {
      this->NbrOneBodyInteractionFactors = 0;
      for (int i = 0; i <= this->LzMax; ++i)
	if (this->OneBodyPotentials[i] != 0)
	  ++this->NbrOneBodyInteractionFactors;
      if (this->NbrOneBodyInteractionFactors != 0)
	{
	  this->OneBodyMValues = new int[this->NbrOneBodyInteractionFactors];
	  this->OneBodyNValues = new int[this->NbrOneBodyInteractionFactors];
	  this->OneBodyInteractionFactors = new double[this->NbrLzValue];
	  this->NbrOneBodyInteractionFactors = 0;
	  int TmpNbrOrbitals = 0;
	  for (int i = 0; i <= this->LzMax; ++i)
	  {
	    this->OneBodyInteractionFactors[TmpNbrOrbitals] = 0.5*this->OneBodyPotentials[i];
	    if (this->OneBodyPotentials[i] != 0)
	      {
		this->OneBodyMValues[this->NbrOneBodyInteractionFactors] = i;
		this->OneBodyNValues[this->NbrOneBodyInteractionFactors] = i;
		++this->NbrOneBodyInteractionFactors;
	      }	 
	    ++ TmpNbrOrbitals;
	  }
	}
      else
	{
	  delete[] this->OneBodyPotentials;
	  this->OneBodyTermFlag = false;
	}
    } 
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the two body delta interaction between two particles 
  //
  // r1 = creation r
  // s1 = creation s
  // r2 = annihilation r
  // s2 = annihilation s
  // return value = corresponding matrix element

double ParticleOnCP2DeltaHamiltonian::ComputeTwoBodyMatrixElement(int r1, int s1, int r2, int s2, int r3, int s3, int r4, int s4, FactorialCoefficient &Coef)
{
  int t1 = this->NbrFluxQuanta - r1 - s1;
  int t2 = this->NbrFluxQuanta - r2 - s2;
  int t3 = this->NbrFluxQuanta - r3 - s3;
  int t4 = this->NbrFluxQuanta - r4 - s4;
  
  Coef.SetToOne();
    
  Coef.FactorialMultiply(r1 + r2);
  Coef.FactorialDivide(r1);
  Coef.FactorialDivide(r2);
  Coef.FactorialMultiply(r1 + r2);
  Coef.FactorialDivide(r3);
  Coef.FactorialDivide(r4);
  Coef.FactorialMultiply(s1 + s2);
  Coef.FactorialDivide(s1);
  Coef.FactorialDivide(s2);
  Coef.FactorialMultiply(s1 + s2);
  Coef.FactorialDivide(s3);
  Coef.FactorialDivide(s4);
  Coef.FactorialMultiply(t1 + t2);
  Coef.FactorialDivide(t1);
  Coef.FactorialDivide(t2);
  Coef.FactorialMultiply(t1 + t2);
  Coef.FactorialDivide(t3);
  Coef.FactorialDivide(t4);
  Coef.FactorialDivide(2*this->NbrFluxQuanta + 2);
  Coef.FactorialDivide(2*this->NbrFluxQuanta + 2);
  Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
  Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
  Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
  Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
    
  double Tmp = sqrt (Coef.GetNumericalValue());
  
  return Tmp;
}
