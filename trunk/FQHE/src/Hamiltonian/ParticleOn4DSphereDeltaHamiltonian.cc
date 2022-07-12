////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Cecile Repellin                   //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//           class of ruby lattice model with interacting particles           //
//                       in the single band approximation                     // 
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
#include "Hamiltonian/ParticleOn4DSphereDeltaHamiltonian.h"
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

ParticleOn4DSphereDeltaHamiltonian::ParticleOn4DSphereDeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->NbrLzValue = - (this->NbrFluxQuanta*(this->NbrFluxQuanta+1)*(2*this->NbrFluxQuanta+1))/6 + this->NbrFluxQuanta*(this->NbrFluxQuanta*(this->NbrFluxQuanta+1))/2 + (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 1);
  this->LzMax = this->NbrLzValue - 1;
  this->quantumNumberJ = new int [this->NbrLzValue];
  this->quantumNumberJz = new int [this->NbrLzValue];
  this->quantumNumberKz = new int [this->NbrLzValue];

  this->GetQuantumNumbersFromLinearizedIndex(quantumNumberJ, quantumNumberJz, quantumNumberKz);

  
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

// destructor
//

ParticleOn4DSphereDeltaHamiltonian::~ParticleOn4DSphereDeltaHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOn4DSphereDeltaHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  FactorialCoefficient Coef;
   
      this->NbrSectorSums = - (2*this->NbrFluxQuanta*(2*this->NbrFluxQuanta+1)*(4*this->NbrFluxQuanta+1))/6
			     + 2*this->NbrFluxQuanta*(2*this->NbrFluxQuanta*(2*this->NbrFluxQuanta+1))/2 
			     + (2*this->NbrFluxQuanta + 1)*(2*this->NbrFluxQuanta + 1);
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;  
      
      for (int j1 = 0; j1 <= this->NbrFluxQuanta; ++j1)
	for (int j2 = 0; j2 <= this->NbrFluxQuanta; ++j2)
	  for (int jz1 = 0; jz1 <= j1; ++jz1)
	    for (int jz2 = 0; jz2 <= j2; ++jz2) 
	      for (int kz1 = 0; kz1 <= this->NbrFluxQuanta - j1; ++kz1)
		for (int kz2 = 0; kz2 <= this->NbrFluxQuanta - j2; ++kz2)
		{
		int Index1 = -((j1-1)*j1*(2*j1-1))/6 + this->NbrFluxQuanta*j1*(j1-1)/2 + (this->NbrFluxQuanta+1)*j1+ (this->NbrFluxQuanta + 1 - j1)*jz1 +kz1;
		int Index2 = -((j2-1)*j2*(2*j2-1))/6 + this->NbrFluxQuanta*j2*(j2-1)/2 + (this->NbrFluxQuanta+1)*j2+ (this->NbrFluxQuanta + 1 - j2)*jz2 +kz2;
		if (Index1 <= Index2)
		  ++this->NbrSectorIndicesPerSum[-(((j1 + j2)-1)*(j1 + j2)*(2*(j1 + j2)-1))/6 + 2*this->NbrFluxQuanta*(j1 + j2)*((j1 + j2)-1)/2 + (2*this->NbrFluxQuanta+1)*(j1 + j2)+ (2*this->NbrFluxQuanta + 1 - (j1 + j2))*(jz1 + jz2) + kz1 +kz2];    
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
      for (int j1 = 0; j1 <= this->NbrFluxQuanta; ++j1)
	for (int j2 = 0; j2 <= this->NbrFluxQuanta; ++j2)
	  for (int jz1 = 0; jz1 <= j1; ++jz1)
	    for (int jz2 = 0; jz2 <= j2; ++jz2) 
	      for (int kz1 = 0; kz1 <= this->NbrFluxQuanta - j1; ++kz1)
		for (int kz2 = 0; kz2 <= this->NbrFluxQuanta - j2; ++kz2)
		  {
		    int Index1 = -((j1-1)*j1*(2*j1-1))/6 + this->NbrFluxQuanta*j1*(j1-1)/2 + (this->NbrFluxQuanta+1)*j1+ (this->NbrFluxQuanta + 1 - j1)*jz1 +kz1;
		    int Index2 = -((j2-1)*j2*(2*j2-1))/6 + this->NbrFluxQuanta*j2*(j2-1)/2 + (this->NbrFluxQuanta+1)*j2+ (this->NbrFluxQuanta + 1 - j2)*jz2 +kz2;
		    if (Index1 <= Index2)
		      {
			int TmpSum = -(((j1 + j2)-1)*(j1 + j2)*(2*(j1 + j2)-1))/6 + 2*this->NbrFluxQuanta*(j1 + j2)*((j1 + j2)-1)/2 + (2*this->NbrFluxQuanta+1)*(j1 + j2)+ (2*this->NbrFluxQuanta + 1 - (j1 + j2))*(jz1 + jz2) + kz1 +kz2;
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
	      int j1 = quantumNumberJ[Index1];
	      int jz1 = quantumNumberJz[Index1];
	      int kz1 = quantumNumberKz[Index1];
	      int j2 = quantumNumberJ[Index2];
	      int jz2 = quantumNumberJz[Index2];
	      int kz2 = quantumNumberKz[Index2];
	      for (int l2 = 0; l2 < this->NbrSectorIndicesPerSum[i]; ++l2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][l2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(l2 << 1) + 1];
		  int j3 = quantumNumberJ[Index3];
		  int jz3 = quantumNumberJz[Index3];
		  int kz3 = quantumNumberKz[Index3];
		  int j4 = quantumNumberJ[Index4];
		  int jz4 = quantumNumberJz[Index4];
		  int kz4 = quantumNumberKz[Index4];
		  
// 		  cout << j1 + j2 << "=" << j3 + j4 << endl;
// 		  cout << jz1 + jz2 << "=" << jz3 + jz4 << endl;
// 		  cout << kz1 + kz2 << "=" << kz3 + kz4 << endl;
		  
 		  this->InteractionFactors[i][Index] =  this->ComputeTwoBodyMatrixElement(j1,  jz1, kz1, j2, jz2, kz2, j3,  jz3, kz3, j4, jz4, kz4, Coef);
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
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the two body delta interaction between two particles 
  //
  // j1 = creation j
  // jz1 = creation jz
  // kz1 = creation kz
  // j2 = annihilation j
  // jz2 = annihilation jz
  // kz2 = annihilation kz
  // return value = corresponding matrix element

double ParticleOn4DSphereDeltaHamiltonian::ComputeTwoBodyMatrixElement(int j1, int jz1, int kz1, int j2, int jz2, int kz2, int j3, int jz3, int kz3, int j4, int jz4, int kz4, FactorialCoefficient &Coef)
{
  Coef.SetToOne();
    
  Coef.FactorialMultiply(j1 + j2 - jz1 - jz2);
  Coef.FactorialDivide(j1 - jz1);
  Coef.FactorialDivide(jz1);
  Coef.FactorialMultiply(j1 + j2 - jz1 - jz2);
  Coef.FactorialDivide(this->NbrFluxQuanta - j1 - kz1);
  Coef.FactorialDivide(kz1);
  Coef.FactorialMultiply(jz1 + jz2);
  Coef.FactorialDivide(j2 - jz2);
  Coef.FactorialDivide(jz2);
  Coef.FactorialMultiply(jz1 + jz2);
  Coef.FactorialDivide(this->NbrFluxQuanta - j2 -kz2);
  Coef.FactorialDivide(kz2);
  Coef.FactorialMultiply(2*this->NbrFluxQuanta - j1 - j2 - kz1 - kz2);
  Coef.FactorialDivide(j3 - jz3);
  Coef.FactorialDivide(jz3);
  Coef.FactorialMultiply(2*this->NbrFluxQuanta - j1 - j2 - kz1 - kz2);
  Coef.FactorialDivide(this->NbrFluxQuanta - j3 -kz3);
  Coef.FactorialDivide(kz3);
  Coef.FactorialMultiply(kz1 + kz2);
  Coef.FactorialDivide(j4 - jz4);
  Coef.FactorialDivide(jz4);
  Coef.FactorialMultiply(kz1 + kz2);
  Coef.FactorialDivide(this->NbrFluxQuanta - j4 -kz4);
  Coef.FactorialDivide(kz4);
  Coef.FactorialDivide(2*this->NbrFluxQuanta + 3);
  Coef.FactorialMultiply(this->NbrFluxQuanta + 3);
  Coef.FactorialMultiply(this->NbrFluxQuanta + 3);
  Coef.FactorialDivide(2*this->NbrFluxQuanta + 3);
  Coef.FactorialMultiply(this->NbrFluxQuanta + 3);
  Coef.FactorialMultiply(this->NbrFluxQuanta + 3);
    
  double Tmp = sqrt (Coef.GetNumericalValue())/(this->NbrFluxQuanta);
  
  return Tmp;
}
