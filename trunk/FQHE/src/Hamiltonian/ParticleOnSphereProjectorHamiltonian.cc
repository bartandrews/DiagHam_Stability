////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of hamiltonian defined as a projector                 //
//                       for particles on a sphere defined                    //
//                                                                            //
//                        last modification : 14/08/2009                      //
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
#include "Hamiltonian/ParticleOnSphereProjectorHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"

  
#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// projectorState = state describing the projector
// projectorSpace = space associated to the projector state 
// projectorNbrParticles = number of particles for the projector state
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereProjectorHamiltonian::ParticleOnSphereProjectorHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
									   RealVector& projectorState, ParticleOnSphere* projectorSpace, int projectorNbrParticles, double l2Factor, 
									   AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
									   char* precalculationFileName )
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = false;
  this->MaxNBody = projectorNbrParticles;
  this->NbrProjectors = 1;
  this->ProjectorSpaces = new   ParticleOnSphere* [this->NbrProjectors];
  this->ProjectorStates = new RealVector [this->NbrProjectors];
  this->ProjectorSpaces[0] = projectorSpace;
  this->ProjectorStates[0] = projectorState;

  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double** [this->MaxNBody + 1];
  this->NbrSortedIndicesPerSum = new int* [this->MaxNBody + 1];
  this->SortedIndicesPerSum = new int** [this->MaxNBody + 1];
  this->MinSumIndices = new int [this->MaxNBody + 1];
  this->MaxSumIndices = new int [this->MaxNBody + 1];
  this->NBodySign = new double[this->MaxNBody + 1];

  this->NbrNIndices = new long[this->MaxNBody + 1];
  this->NIndices = new int*[this->MaxNBody + 1];
  this->NbrMIndices = new long*[this->MaxNBody + 1];
  this->MIndices = new int**[this->MaxNBody + 1];
  this->MNNBodyInteractionFactors = new double**[this->MaxNBody + 1];

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 1;
      this->MaxSumIndices[k] = 0;      
      this->NBodyFlags[k] = false;
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	{
	  this->NBodySign[k] = -1.0;
	}
      this->NbrNIndices[k] = 0;
      this->NIndices[k] = 0;
      this->NbrMIndices[k] = 0;
      this->MIndices[k] = 0;
      this->MNNBodyInteractionFactors[k] = 0;
    }

  this->NBodyFlags[this->MaxNBody] = true;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
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
      else
	{
	  this->FastMultiplicationFlag = false;
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);

  if (l2Factor != 0.0)
    {
      this->L2Operator = new ParticleOnSphereL2Hamiltonian(this->Particles, this->NbrParticles, this->LzMax, this->Particles->GetLzValue() , this->Architecture, l2Factor, memory); 
    }
  else
    {
      this->L2Operator = 0;
    }

}


// constructor from multiple projectors
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// projectorState = states describing each projector
// projectorSpace = spaces associated to the projector states
// projectorNbrParticles = number of projectors
// projectorNbrParticles = number of particles for the projector state
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereProjectorHamiltonian::ParticleOnSphereProjectorHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
									   RealVector* projectorStates, ParticleOnSphere** projectorSpaces, int nbrProjectors, int projectorNbrParticles, double l2Factor, 
									   AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
									   char* precalculationFileName )
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = false;
  this->MaxNBody = projectorNbrParticles;
  this->NbrProjectors = nbrProjectors;
  this->ProjectorSpaces = new   ParticleOnSphere* [this->NbrProjectors];
  this->ProjectorStates = new RealVector [this->NbrProjectors];
  for (int i = 0; i < this->NbrProjectors; ++i)
    {
      this->ProjectorSpaces[i] = projectorSpaces[i];
      this->ProjectorStates[i] = projectorStates[i];
    }

  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double** [this->MaxNBody + 1];
  this->NbrSortedIndicesPerSum = new int* [this->MaxNBody + 1];
  this->SortedIndicesPerSum = new int** [this->MaxNBody + 1];
  this->MinSumIndices = new int [this->MaxNBody + 1];
  this->MaxSumIndices = new int [this->MaxNBody + 1];
  this->NBodySign = new double[this->MaxNBody + 1];

  this->NbrNIndices = new long[this->MaxNBody + 1];
  this->NIndices = new int*[this->MaxNBody + 1];
  this->NbrMIndices = new long*[this->MaxNBody + 1];
  this->MIndices = new int**[this->MaxNBody + 1];
  this->MNNBodyInteractionFactors = new double**[this->MaxNBody + 1];

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 1;
      this->MaxSumIndices[k] = 0;      
      this->NBodyFlags[k] = false;
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	{
	  this->NBodySign[k] = -1.0;
	}
      this->NbrNIndices[k] = 0;
      this->NIndices[k] = 0;
      this->NbrMIndices[k] = 0;
      this->MIndices[k] = 0;
      this->MNNBodyInteractionFactors[k] = 0;
    }

  this->NBodyFlags[this->MaxNBody] = true;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
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
      else
	{
	  this->FastMultiplicationFlag = false;
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);

  if (l2Factor != 0.0)
    {
      this->L2Operator = new ParticleOnSphereL2Hamiltonian(this->Particles, this->NbrParticles, this->LzMax, this->Particles->GetLzValue() , this->Architecture, l2Factor); 
    }
  else
    {
      this->L2Operator = 0;
    }

}

// destructor
//

ParticleOnSphereProjectorHamiltonian::~ParticleOnSphereProjectorHamiltonian()
{
  for (int k = 1; k <= this->MaxNBody; ++k)
    if (this->NBodyFlags[k] == true)
      {
	for (int MinSum = this->MinSumIndices[k]; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	  {
	    delete[] this->SortedIndicesPerSum[k][MinSum];
	    if (this->MNNBodyInteractionFactors == 0)
	      delete[] this->NBodyInteractionFactors[k][MinSum];	      
	  }
	delete[] this->NbrSortedIndicesPerSum[k];
	delete[] this->SortedIndicesPerSum[k];
	if (this->MNNBodyInteractionFactors == 0)
	  delete[] this->NBodyInteractionFactors[k];
	else
	  {
	    for (int i = 0; i < this->NbrNIndices[k]; ++i)
	      {
		delete[] this->MNNBodyInteractionFactors[k][i];		
		delete[] this->MIndices[k][i];
	      }
	    delete[] this->NIndices[k];
	    delete[] this->NbrMIndices[k];
	    delete[] this->MIndices[k];
	    delete[] this->MNNBodyInteractionFactors[k];
	  }
      }


  delete[] this->ProjectorSpaces;
  delete[] this->ProjectorStates;

  delete[] this->NbrNIndices;
  delete[] this->NIndices;
  delete[] this->NbrMIndices;
  delete[] this->MIndices;
  delete[] this->MNNBodyInteractionFactors;

  delete[] this->NBodyFlags;
  delete[] this->NBodyInteractionFactors;
  delete[] this->SortedIndicesPerSum;
  delete[] this->NbrSortedIndicesPerSum;
  delete[] this->MinSumIndices;
  delete[] this->MaxSumIndices;
  delete[] this->NBodySign;
  if (this->L2Operator != 0)
    delete this->L2Operator;

}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* ParticleOnSphereProjectorHamiltonian::Clone ()
{
  return 0;
}


// evaluate all interaction factors
//   

void ParticleOnSphereProjectorHamiltonian::EvaluateInteractionFactors()
{
  unsigned long* TmpMonomial = new unsigned long[this->MaxNBody];
  int TmpMinSum = 0;
  int TmpMaxSum = 0;
  int* SumPerProjectorSpace = new int [this->NbrProjectors];
  double** IndexSymmetryFactor = new double* [this->NbrProjectors];
  for (int i = 0; i < this->NbrProjectors; ++i)
    IndexSymmetryFactor[i] = new double [this->ProjectorSpaces[i]->GetHilbertSpaceDimension()];
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int j = 0; j < this->NbrProjectors; ++j)
	{
	  FermionOnSphere* TmpSpace = (FermionOnSphere*) this->ProjectorSpaces[j];
	  TmpSpace->GetMonomial(0, TmpMonomial);
	  SumPerProjectorSpace[j] = 0;
	  for (int i = 0; i < this->MaxNBody; ++i)
	    SumPerProjectorSpace[j] += TmpMonomial[i];
	}
      TmpMaxSum = SumPerProjectorSpace[0];
      TmpMinSum = SumPerProjectorSpace[0];
      for (int j = 0; j < this->NbrProjectors; ++j)
	{
	  if (SumPerProjectorSpace[j] > TmpMaxSum)
	    TmpMaxSum = SumPerProjectorSpace[j];
	  if (SumPerProjectorSpace[j] < TmpMinSum)
	    TmpMinSum = SumPerProjectorSpace[j];
	}
      this->MaxSumIndices[this->MaxNBody] = TmpMaxSum;
      this->MinSumIndices[this->MaxNBody] = TmpMinSum;
      this->SortedIndicesPerSum[this->MaxNBody] = new int*[TmpMaxSum + 1];
      for (int i = 0; i <= TmpMaxSum; ++i)
	this->SortedIndicesPerSum[this->MaxNBody][i] = 0;
      for (int j = 0; j < this->NbrProjectors; ++j)     
	{
	  this->SortedIndicesPerSum[this->MaxNBody][SumPerProjectorSpace[j]] = new int [this->ProjectorSpaces[j]->GetHilbertSpaceDimension() * this->MaxNBody];
	  int* TmpSortedIndicesPerSum = this->SortedIndicesPerSum[this->MaxNBody][SumPerProjectorSpace[j]];
	  FermionOnSphere* TmpSpace = (FermionOnSphere*) this->ProjectorSpaces[j];
	  for (int i = 0; i < TmpSpace->GetHilbertSpaceDimension(); ++i)
	    {
	      TmpSpace->GetMonomial(i, TmpMonomial);	  
	      for (int k = 0; k < this->MaxNBody; ++k)
		{
		  (*TmpSortedIndicesPerSum) = TmpMonomial[k];
		  ++TmpSortedIndicesPerSum;
		}
	      IndexSymmetryFactor[j][i] = 1.0;
	    }
	}
    }
  else  
    {
      double* SymmetryFactor = new double [this->MaxNBody + 1];
      SymmetryFactor[0] = 1.0;
      SymmetryFactor[1] = 1.0;
      for (int i = 2; i <= this->MaxNBody; ++i)
	SymmetryFactor[i] = SymmetryFactor[i - 1] / sqrt ((double) i);
      this->MaxSumIndices[this->MaxNBody] = 0;
      BosonOnSphereShort* TmpSpace = 0;
      for (int i = 0; i < this->NbrProjectors; ++i)
	{
	  TmpSpace = (BosonOnSphereShort*) this->ProjectorSpaces[i];
	  TmpSpace->GetMonomial(0, TmpMonomial);
	  SumPerProjectorSpace[i] = 0;
	  for (int j = 0; j < this->MaxNBody; ++j)
	    SumPerProjectorSpace[i] += TmpMonomial[j];
	}
      TmpMaxSum = SumPerProjectorSpace[0];
      TmpMinSum = SumPerProjectorSpace[0];
      for (int j = 0; j < this->NbrProjectors; ++j)
	{
	  if (SumPerProjectorSpace[j] > TmpMaxSum)
	    TmpMaxSum = SumPerProjectorSpace[j];
	  if (SumPerProjectorSpace[j] < TmpMinSum)
	    TmpMinSum = SumPerProjectorSpace[j];
	}
      this->MaxSumIndices[this->MaxNBody] = TmpMaxSum;
      this->MinSumIndices[this->MaxNBody] = TmpMinSum;
      this->SortedIndicesPerSum[this->MaxNBody] = new int*[TmpMaxSum + 1];
      for (int i = 0; i <= TmpMaxSum; ++i)	
	this->SortedIndicesPerSum[this->MaxNBody][i] = 0;     
      for (int j = 0; j < this->NbrProjectors; ++j)
	{
	  TmpSpace = (BosonOnSphereShort*) this->ProjectorSpaces[j];
	  this->SortedIndicesPerSum[this->MaxNBody][SumPerProjectorSpace[j]] = new int [TmpSpace->GetHilbertSpaceDimension() * this->MaxNBody];
	  int* TmpSortedIndicesPerSum = this->SortedIndicesPerSum[this->MaxNBody][SumPerProjectorSpace[j]];
	  for (int i = 0; i < TmpSpace->GetHilbertSpaceDimension(); ++i)
	    {
	      TmpSpace->GetMonomial(i, TmpMonomial);	  
	      IndexSymmetryFactor[j][i] = 1.0;
	      unsigned long CurrentOccupation = 0ul;
	      int NbrOccupation = 0;
	      for (int k = 0; k < this->MaxNBody; ++k)
		{
		  (*TmpSortedIndicesPerSum) = TmpMonomial[k];
		  ++TmpSortedIndicesPerSum;
		  if (CurrentOccupation != TmpMonomial[k])
		    {
		      IndexSymmetryFactor[j][i] *= SymmetryFactor[NbrOccupation];
		      CurrentOccupation =  TmpMonomial[k];
		      NbrOccupation = 1;
		    }
		  else
		    ++NbrOccupation;
		}
	      IndexSymmetryFactor[j][i] *= SymmetryFactor[NbrOccupation];	  
	    }
	}
      delete[] SymmetryFactor;
    }

  this->NbrSortedIndicesPerSum[this->MaxNBody] = new int [TmpMaxSum + 1];
  for (int i = 0; i < TmpMaxSum; ++i)
    this->NbrSortedIndicesPerSum[this->MaxNBody][i] = 0;
  int TmpTotalNbrNIndices = 0;
  for (int i = 0; i < this->NbrProjectors; ++i)
    {
      this->NbrSortedIndicesPerSum[this->MaxNBody][SumPerProjectorSpace[i]] = this->ProjectorSpaces[i]->GetHilbertSpaceDimension();
      TmpTotalNbrNIndices += this->ProjectorSpaces[i]->GetHilbertSpaceDimension();
    }

  this->NbrNIndices[this->MaxNBody] = TmpTotalNbrNIndices;
  this->NIndices[this->MaxNBody] = new int[TmpTotalNbrNIndices * this->MaxNBody];
  this->NbrMIndices[this->MaxNBody] = new long[TmpTotalNbrNIndices];
  this->MIndices[this->MaxNBody] = new int*[TmpTotalNbrNIndices];
  this->MNNBodyInteractionFactors[this->MaxNBody] = new double* [TmpTotalNbrNIndices];

  int TmpNbrNIndices = 0;	 
  int* TmpNIndices = this->NIndices[this->MaxNBody];
  for (int k = 0; k < this->NbrProjectors; ++k)
    {
      int* TmpNIndices2 = this->SortedIndicesPerSum[this->MaxNBody][SumPerProjectorSpace[k]];
      for (int i = 0; i < this->ProjectorSpaces[k]->GetHilbertSpaceDimension(); ++i)
	{
	  this->NbrMIndices[this->MaxNBody][TmpNbrNIndices] = this->ProjectorSpaces[k]->GetHilbertSpaceDimension();		    
	  this->MIndices[this->MaxNBody][TmpNbrNIndices] = new int [this->ProjectorSpaces[k]->GetHilbertSpaceDimension() * this->MaxNBody];
	  int* TmpMIndices = this->MIndices[this->MaxNBody][TmpNbrNIndices];
	  int* TmpMIndices2 = this->SortedIndicesPerSum[this->MaxNBody][SumPerProjectorSpace[k]];
	  this->MNNBodyInteractionFactors[this->MaxNBody][TmpNbrNIndices] = new double [this->ProjectorSpaces[k]->GetHilbertSpaceDimension()];
	  double* TmpInteraction = this->MNNBodyInteractionFactors[this->MaxNBody][TmpNbrNIndices];
	  for (int j = 0; j < this->ProjectorSpaces[k]->GetHilbertSpaceDimension(); ++j)
	    {
	      for (int l = 0; l < this->MaxNBody; ++l)
		{
		  (*TmpMIndices) = (*TmpMIndices2);			
		  ++TmpMIndices;
		  ++TmpMIndices2;
		}			
	      TmpInteraction[j] = this->ProjectorStates[k][i] * this->ProjectorStates[k][j] * IndexSymmetryFactor[k][i] * IndexSymmetryFactor[k][j];
	    }
	  for (int j = 0; j < this->MaxNBody; ++j)
	    {
	      (*TmpNIndices) = (*TmpNIndices2);			
	      ++TmpNIndices;
	      ++TmpNIndices2;
	    }
	  ++TmpNbrNIndices;
	}
    }

}

