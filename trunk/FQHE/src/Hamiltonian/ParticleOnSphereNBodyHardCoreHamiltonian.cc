////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                         n-body hard core interaction                       //
//                                                                            //
//                        last modification : 23/09/2004                      //
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
#include "Hamiltonian/ParticleOnSphereNBodyHardCoreHamiltonian.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"

#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//

ParticleOnSphereNBodyHardCoreHamiltonian::ParticleOnSphereNBodyHardCoreHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// nbrBody = number of particle that interact simultaneously through the hard core interaction
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereNBodyHardCoreHamiltonian::ParticleOnSphereNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int nbrBody, double l2Factor,
										   AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag,
										   char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = false;
  this->NbrNbody = nbrBody;
  this->MaxNBody = this->NbrNbody;
  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double** [this->MaxNBody + 1];
  this->NBodyInteractionWeightFactors = new double [this->MaxNBody + 1];
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
      this->NBodyInteractionWeightFactors[k] = 0.0;
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	{
	  this->NBodySign[k] = -1.0;
	  if (k == 4)
	    this->NBodySign[k] = 1.0;
	}
      this->NbrNIndices[k] = 0;
      this->NIndices[k] = 0;
      this->NbrMIndices[k] = 0;
      this->MIndices[k] = 0;
      this->MNNBodyInteractionFactors[k] = 0;
    }

  this->NBodyFlags[this->NbrNbody] = true;
  this->NBodyInteractionWeightFactors[this->NbrNbody] = 1.0;
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

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// architecture = architecture to use for precalculation
// maxNbrBody = maximum number of particle that interact simultaneously through the hard core interaction
// nBodyFactors = weight of the different n-body interaction terms with respect to each other
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereNBodyHardCoreHamiltonian::ParticleOnSphereNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
										   int maxNbrBody, double* nBodyFactors, double l2Factor,
										   AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
										   char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = false;
  this->NbrNbody = maxNbrBody;
  this->MaxNBody = maxNbrBody;
  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double** [this->MaxNBody + 1];
  this->NBodyInteractionWeightFactors = new double [this->MaxNBody + 1];
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
      if (nBodyFactors[k] == 0.0)
	this->NBodyFlags[k] = false;
      else
	this->NBodyFlags[k] = true;
      this->NBodyInteractionWeightFactors[k] = nBodyFactors[k];
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	this->NBodySign[k] = -1.0;
      if (this->NBodyInteractionWeightFactors[k] < 0.0)
	{
	  this->NBodyInteractionWeightFactors[k] *= -1.0;
	  this->NBodySign[k] *= -1.0;
	}
    }
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
	    cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
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

// constructor from datas with a fully-defined two body interaction
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// maxNbrBody = maximum number of particle that interact simultaneously through the hard core interaction
// nBodyFactors = weight of the different n-body interaction terms with respect to each other
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereNBodyHardCoreHamiltonian::ParticleOnSphereNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
										   int maxNbrBody, double* nBodyFactors, double l2Factor, double* pseudoPotential, 
										   AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
										   char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  if (pseudoPotential != 0)
    {
      this->PseudoPotential = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->PseudoPotential[i] = pseudoPotential[this->LzMax - i];
      this->FullTwoBodyFlag = true;
    }
  else
    {
      this->FullTwoBodyFlag = false;
    }

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = true;
  this->NbrNbody = maxNbrBody;
  this->MaxNBody = maxNbrBody;
  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double** [this->MaxNBody + 1];
  this->NBodyInteractionWeightFactors = new double [this->MaxNBody + 1];
  this->NbrSortedIndicesPerSum = new int* [this->MaxNBody + 1];
  this->SortedIndicesPerSum = new int** [this->MaxNBody + 1];
  this->MinSumIndices = new int [this->MaxNBody + 1];
  this->MaxSumIndices = new int [this->MaxNBody + 1];
  this->NBodySign = new double[this->MaxNBody + 1];
  this->PseudoPotential = new double [this->NbrLzValue];

  this->NbrNIndices = new long[this->MaxNBody + 1];
  this->NIndices = new int*[this->MaxNBody + 1];
  this->NbrMIndices = new long*[this->MaxNBody + 1];
  this->MIndices = new int**[this->MaxNBody + 1];
  this->MNNBodyInteractionFactors = new double**[this->MaxNBody + 1];

  for (int i = 0; i < this->NbrLzValue; ++i)
    this->PseudoPotential[i] = pseudoPotential[this->LzMax - i];

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 1;
      this->MaxSumIndices[k] = 0;      
      if (nBodyFactors[k] == 0.0)
	this->NBodyFlags[k] = false;
      else
	this->NBodyFlags[k] = true;
      this->NBodyInteractionWeightFactors[k] = nBodyFactors[k];
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	this->NBodySign[k] = -1.0;
      if (this->NBodyInteractionWeightFactors[k] < 0.0)
	{
	  this->NBodyInteractionWeightFactors[k] *= -1.0;
	  this->NBodySign[k] *= -1.0;
	}
    }
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
	    cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
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

ParticleOnSphereNBodyHardCoreHamiltonian::~ParticleOnSphereNBodyHardCoreHamiltonian()
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

  if (this->MNNBodyInteractionFactors != 0)
    {
      delete[] this->NbrNIndices;
      delete[] this->NIndices;
      delete[] this->NbrMIndices;
      delete[] this->MIndices;
      delete[] this->MNNBodyInteractionFactors;
    }

  delete[] this->NBodyFlags;
  delete[] this->NBodyInteractionFactors;
  delete[] this->SortedIndicesPerSum;
  delete[] this->NbrSortedIndicesPerSum;
  delete[] this->NBodyInteractionWeightFactors;
  delete[] this->MinSumIndices;
  delete[] this->MaxSumIndices;
  delete[] this->NBodySign;
  if (this->L2Operator != 0)
    delete this->L2Operator;
  if (this->FullTwoBodyFlag == true)
    {
      delete[] this->InteractionFactors;
      delete[] this->M1Value;
      delete[] this->M2Value;
      delete[] this->M3Value;
      delete[] this->PseudoPotential;
    }
}

// evaluate all interaction factors
//   

void ParticleOnSphereNBodyHardCoreHamiltonian::EvaluateInteractionFactors()
{
  double* TmpNormalizationCoeffients = new double[this->NbrLzValue];
  double TmpFactor = ((double) this->NbrLzValue) / (4.0 * M_PI);
  double TmpBinomial = 1.0;
  TmpNormalizationCoeffients[0] = sqrt (TmpBinomial * TmpFactor);
  for (int i = 1; i < this->NbrLzValue; ++i)
    {
      TmpBinomial *= this->LzMax - ((double) i) + 1.0;
      TmpBinomial /= ((double) i);
      TmpNormalizationCoeffients[i] = sqrt (TmpBinomial * TmpFactor);
    }
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    double Coefficient;
	    GetAllSkewSymmetricIndices(this->NbrLzValue, k, this->NbrSortedIndicesPerSum[k], this->SortedIndicesPerSum[k]);
	    this->MaxSumIndices[k] = (this->LzMax * k) - (((k - 1) * (k - 2))/ 2);
	    this->MinSumIndices[k] = (k * (k - 1)) / 2;
	    double* TmpInteractionCoeffients = new double[this->MaxSumIndices[k] + 1];
	    TmpInteractionCoeffients[0] = 1.0;
	    TmpInteractionCoeffients[1] = 1.0;
	    for (int i = 2; i <= this->MaxSumIndices[k]; ++i)
	      {
		Coefficient = 1.0;
		for (int j = 1; j < i; ++j)
		  {
		    double Coefficient2 = TmpInteractionCoeffients[j];
		    TmpInteractionCoeffients[j] += Coefficient;
		    Coefficient = Coefficient2;
		  }
		TmpInteractionCoeffients[i] = 1.0;
	      }
	    Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[k]) + 1.0);
	    double Radius = 2.0 / ((double) this->LzMax);
	    for (int i = 2; i <= k; ++i)
	      {
		Coefficient *= (double) (i * i);	  
		Coefficient *= Radius;
	      }
	    for (int i = 0; i <= this->MaxSumIndices[k]; ++i)
	      TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);

	    long TmpNbrNIndices = 0;
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	      TmpNbrNIndices += this->NbrSortedIndicesPerSum[k][MinSum];
	    this->NbrNIndices[k] = TmpNbrNIndices;
	    this->NIndices[k] = new int[TmpNbrNIndices * k];
	    this->NbrMIndices[k] = new long[TmpNbrNIndices];
	    this->MIndices[k] = new int*[TmpNbrNIndices];
	    this->MNNBodyInteractionFactors[k] = new double* [TmpNbrNIndices];
	    TmpNbrNIndices = 0;	 
	    int* TmpNIndices = this->NIndices[k];
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	      {
		int Lim = this->NbrSortedIndicesPerSum[k][MinSum];
		int* TmpNIndices2 = this->SortedIndicesPerSum[k][MinSum];
		double* TmpProjectorCoefficients = this->ComputeProjectorCoefficients(k * (k - 1), k, TmpNIndices2, Lim);
		for (int i = 0; i < Lim; ++i)
		  {
		    this->NbrMIndices[k][TmpNbrNIndices] = Lim;		    
		    this->MIndices[k][TmpNbrNIndices] = new int [Lim * k];
		    this->MNNBodyInteractionFactors[k][TmpNbrNIndices] = new double [Lim];
		    int* TmpMIndices = this->MIndices[k][TmpNbrNIndices];
		    int* TmpMIndices2 = this->SortedIndicesPerSum[k][MinSum];
		    double* TmpInteraction = this->MNNBodyInteractionFactors[k][TmpNbrNIndices];
		    Coefficient = TmpProjectorCoefficients[i]* this->NBodyInteractionWeightFactors[k];
		    for (int j = 0; j < Lim; ++j)
		      {
			for (int l = 0; l < k; ++l)
			  {
			    (*TmpMIndices) = (*TmpMIndices2);			
			    ++TmpMIndices;
			    ++TmpMIndices2;
			  }			
			TmpInteraction[j] = Coefficient * TmpProjectorCoefficients[j];
		      }
		    for (int j = 0; j < k; ++j)
		      {
			(*TmpNIndices) = (*TmpNIndices2);			
			++TmpNIndices;
			++TmpNIndices2;
		      }
		    ++TmpNbrNIndices;
		  }
		delete[] TmpProjectorCoefficients;		
	      }


// 	    this->NBodyInteractionFactors[k] = new double* [this->MaxSumIndices[k] + 1];
// 	    int Lim;
// 	    for (int MinSum = this->MinSumIndices[k]; MinSum <= this->MaxSumIndices[k]; ++MinSum)
// 	      {
// 		Lim = this->NbrSortedIndicesPerSum[k][MinSum];
// 		this->NBodyInteractionFactors[k][MinSum] = new double [Lim];
// 		double* TmpNBodyInteractionFactors = this->NBodyInteractionFactors[k][MinSum];		  
// 		int* TmpMIndices = this->SortedIndicesPerSum[k][MinSum];
// 		double* TmpProjectorCoefficients = this->ComputeProjectorCoefficients(k * (k - 1), k, TmpMIndices, Lim);
// 		for (int i = 0; i < Lim; ++i)
// 		  {
// 		    Coefficient = TmpInteractionCoeffients[MinSum] * TmpProjectorCoefficients[i];
// 		    for (int l = 0; l < k; ++l)
// 		      {
// 			Coefficient *= TmpNormalizationCoeffients[TmpMIndices[l]];		    
// 		      }
// 		    TmpNBodyInteractionFactors[i] = TmpProjectorCoefficients[i] * sqrt(this->NBodyInteractionWeightFactors[k]);
// 		    TmpMIndices += k;
// 		  }
// 		delete[] TmpProjectorCoefficients;
// 	      }
	    delete[] TmpInteractionCoeffients;
	  }
    }
  else
    {
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    this->MinSumIndices[k] = 0;
	    this->MaxSumIndices[k] = this->LzMax * k;
	    double* TmpInteractionCoeffients = new double[this->MaxSumIndices[k] + 1];
	    double Coefficient;
	    TmpInteractionCoeffients[0] = 1.0;
	    TmpInteractionCoeffients[1] = 1.0;
	    for (int i = 2; i <= this->MaxSumIndices[k]; ++i)
	      {
		Coefficient = 1.0;
		for (int j = 1; j < i; ++j)
		  {
		    double Coefficient2 = TmpInteractionCoeffients[j];
		    TmpInteractionCoeffients[j] += Coefficient;
		    Coefficient = Coefficient2;
		  }
		TmpInteractionCoeffients[i] = 1.0;
	      }
	    Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[k]) + 1.0);
	    double Radius = 2.0 / ((double) this->LzMax);
	    for (int i = 2; i <= k; ++i)
	      {
		Coefficient *= (double) (i * i);	  
		Coefficient *= Radius;
	      }
	    for (int i = 0; i <= this->MaxSumIndices[k]; ++i)
	      TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);

	    double** SortedIndicesPerSumSymmetryFactor;
	    GetAllSymmetricIndices(this->NbrLzValue, k, this->NbrSortedIndicesPerSum[k], this->SortedIndicesPerSum[k],
				   SortedIndicesPerSumSymmetryFactor);


	    long TmpNbrNIndices = 0;
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	      TmpNbrNIndices += this->NbrSortedIndicesPerSum[k][MinSum];
	    this->NbrNIndices[k] = TmpNbrNIndices;
	    this->NIndices[k] = new int[TmpNbrNIndices * k];
	    this->NbrMIndices[k] = new long[TmpNbrNIndices];
	    this->MIndices[k] = new int*[TmpNbrNIndices];
	    this->MNNBodyInteractionFactors[k] = new double* [TmpNbrNIndices];
	    TmpNbrNIndices = 0;	 
	    int* TmpNIndices = this->NIndices[k];
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	      {
		int Lim = this->NbrSortedIndicesPerSum[k][MinSum];
		double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
		int* TmpNIndices2 = this->SortedIndicesPerSum[k][MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    this->NbrMIndices[k][TmpNbrNIndices] = Lim;		    
		    this->MIndices[k][TmpNbrNIndices] = new int [Lim * k];
		    this->MNNBodyInteractionFactors[k][TmpNbrNIndices] = new double [Lim];
		    int* TmpMIndices = this->MIndices[k][TmpNbrNIndices];
		    int* TmpMIndices2 = this->SortedIndicesPerSum[k][MinSum];
		    double* TmpInteraction = this->MNNBodyInteractionFactors[k][TmpNbrNIndices];
		    Coefficient = TmpSymmetryFactors[i] * this->NBodyInteractionWeightFactors[k] * TmpInteractionCoeffients[MinSum] * TmpInteractionCoeffients[MinSum];
		    for (int l = 0; l < k; ++l)
		      Coefficient *= TmpNormalizationCoeffients[TmpNIndices2[l]];		    		      
		    for (int j = 0; j < Lim; ++j)
		      {
			double Coefficient2 = TmpSymmetryFactors[j];
			for (int l = 0; l < k; ++l)
			  {
			    Coefficient2 *= TmpNormalizationCoeffients[(*TmpMIndices2)];		    
			    (*TmpMIndices) = (*TmpMIndices2);			
			    ++TmpMIndices;
			    ++TmpMIndices2;
			  }			
			TmpInteraction[j] = Coefficient * Coefficient2;
		      }
		    for (int j = 0; j < k; ++j)
		      {
			(*TmpNIndices) = (*TmpNIndices2);			
			++TmpNIndices;
			++TmpNIndices2;
		      }
		    ++TmpNbrNIndices;
		  }		
	      }

// 	    this->NBodyInteractionFactors[k] = new double* [this->MaxSumIndices[k] + 1];
//  	    int Lim;
// 	    for (int MinSum = 0; MinSum <= this->MaxSumIndices[k]; ++MinSum)
// 	      {
// 		Lim = this->NbrSortedIndicesPerSum[k][MinSum];
// 		double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
// 		this->NBodyInteractionFactors[k][MinSum] = new double [Lim];
// 		double* TmpNBodyInteractionFactors = this->NBodyInteractionFactors[k][MinSum];		  
// 		int* TmpMIndices = this->SortedIndicesPerSum[k][MinSum];
// 		for (int i = 0; i < Lim; ++i)
// 		  {
// 		    Coefficient = TmpSymmetryFactors[i] * TmpInteractionCoeffients[MinSum];
// 		    for (int l = 0; l < k; ++l)
// 		      Coefficient *= TmpNormalizationCoeffients[TmpMIndices[l]];		    
// 		    TmpNBodyInteractionFactors[i] = Coefficient * sqrt(this->NBodyInteractionWeightFactors[k]);
// 		    TmpMIndices += k;
// 		  }
// 	      }
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	      {
		delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	      }
	    delete[] SortedIndicesPerSumSymmetryFactor;
	    delete[] TmpInteractionCoeffients;
	  }
    }
  delete[] TmpNormalizationCoeffients;
  if (this->FullTwoBodyFlag == true)
    {
      int Lim;
      int Min;
      int Pos = 0;
      ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
      int J = 2 * this->LzMax - 2;
      int m4;
      double ClebschCoef;
      double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
      
      int Sign = 1;
      if (this->LzMax & 1)
	Sign = 0;
      double MaxCoefficient = 0.0;
      
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	{
	  for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
	    for (int m2 =  -this->LzMax; m2 < m1; m2 += 2)
	      {
		Lim = m1 + m2 + this->LzMax;
		if (Lim > this->LzMax)
		  Lim = this->LzMax;
		Min = m1 + m2 - this->LzMax;
		if (Min < -this->LzMax)
		  Min = -this->LzMax;
		for (int m3 = Min; m3 <= Lim; m3 += 2)
		  {
		    Clebsch.InitializeCoefficientIterator(m1, m2);
		    m4 = m1 + m2 - m3;
		    TmpCoefficient[Pos] = 0.0;
		    while (Clebsch.Iterate(J, ClebschCoef))
		      {
			if (((J >> 1) & 1) == Sign)
			  TmpCoefficient[Pos] += this->PseudoPotential[J >> 1] * ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		      }
		    if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		      MaxCoefficient = TmpCoefficient[Pos];
		    ++Pos;
		  }
	      }
	  this->NbrInteractionFactors = 0;
	  this->M1Value = new int [Pos];
	  this->M2Value = new int [Pos];
	  this->M3Value = new int [Pos];
	  this->InteractionFactors = new double [Pos];
	  cout << "nbr interaction = " << Pos << endl;
	  Pos = 0;
	  MaxCoefficient *= MACHINE_PRECISION;
	  double Factor = - 4.0;// / sqrt (0.5 * ((double) this->LzMax));
// 	  for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
// 	    for (int m2 = 0; m2 < m1; ++m2)
// 	      {
// 		Lim = m1 + m2;
// 		if (Lim > this->LzMax)
// 		  Lim = this->LzMax;
// 		Min = m1 + m2 - this->LzMax;
// 		if (Min < 0)
// 		  Min = 0;
// 		for (int m3 = Min; m3 <= Lim; ++m3)
// 		  {
// 		    if (/*(fabs(TmpCoefficient[Pos]) > MaxCoefficient) &&*/ ((2 * m3) > (m1 + m2)))
// 		      {
// 			this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
// 			this->M1Value[this->NbrInteractionFactors] = m1;
// 			this->M2Value[this->NbrInteractionFactors] = m2;
// 			this->M3Value[this->NbrInteractionFactors] = m3;
// 			++this->NbrInteractionFactors;
// 		      }
// 		    ++Pos;
// 		  }
// 	      }
	  this->NbrM12Indices = (this->NbrLzValue * (this->NbrLzValue - 1)) / 2;
	  this->M1Value = new int [this->NbrM12Indices];
	  this->M2Value = new int [this->NbrM12Indices];
	  this->NbrM3Values = new int [this->NbrM12Indices];
	  this->M3Values = new int* [this->NbrM12Indices];
	  int TotalIndex = 0;
	  Pos = 0;
	  for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	    for (int m2 = 0; m2 < m1; ++m2)
	      {
		Lim = m1 + m2;
		if (Lim > this->LzMax)
		  Lim = this->LzMax;
		Min = m1 + m2 - this->LzMax;
		if (Min < 0)
		  Min = 0;
		this->M1Value[TotalIndex] = m1;
		this->M2Value[TotalIndex] = m2;	    
		this->NbrM3Values[TotalIndex] = 0;
		for (int m3 = Min; m3 <= Lim; ++m3)
		  if ((2 * m3) > (m1 + m2))
		    ++this->NbrM3Values[TotalIndex];
		if (this->NbrM3Values[TotalIndex] > 0)
		  {
		    this->M3Values[TotalIndex] = new int [this->NbrM3Values[TotalIndex]];
		    int TmpIndex = 0;
		    for (int m3 = Min; m3 <= Lim; ++m3)
		      {
			if ((2 * m3) > (m1 + m2))
			  {
			    this->M3Values[TotalIndex][TmpIndex] = m3;
			    this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
			    ++this->NbrInteractionFactors;
			    ++TmpIndex;
			  }
			++Pos;
		      }
		  }
		++TotalIndex;
	      }
	}
      else
	{
	  for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
	    for (int m2 =  -this->LzMax; m2 <= m1; m2 += 2)
	      {
		Lim = m1 + m2 + this->LzMax;
		if (Lim > this->LzMax)
		  Lim = this->LzMax;
		Min = m1 + m2 - this->LzMax;
		if (Min < -this->LzMax)
		  Min = -this->LzMax;
		for (int m3 = Min; m3 <= Lim; m3 += 2)
		  {
		    Clebsch.InitializeCoefficientIterator(m1, m2);
		    m4 = m1 + m2 - m3;
		    TmpCoefficient[Pos] = 0.0;
		    while (Clebsch.Iterate(J, ClebschCoef))
		      {
			if (((J >> 1) & 1) != Sign)
			  TmpCoefficient[Pos] += this->PseudoPotential[J >> 1] * ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		      }
		    if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		      MaxCoefficient = TmpCoefficient[Pos];
		    ++Pos;
		  }
	      }
	  this->NbrInteractionFactors = 0;
	  this->M1Value = new int [Pos];
	  this->M2Value = new int [Pos];
	  this->M3Value = new int [Pos];
	  this->InteractionFactors = new double [Pos];
	  cout << "nbr interaction = " << Pos << endl;
	  Pos = 0;
	  MaxCoefficient *= MACHINE_PRECISION;
	  double Factor = 4.0 / sqrt (0.5 * ((double) this->LzMax));
	  for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	    {
	      for (int m2 = 0; m2 < m1; ++m2)
		{
		  Lim = m1 + m2;
		  if (Lim > this->LzMax)
		    Lim = this->LzMax;
		  Min = m1 + m2 - this->LzMax;
		  if (Min < 0)
		    Min = 0;
		  for (int m3 = Min; m3 <= Lim; ++m3)
		    {
		      if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
			{
			  if ((2 * m3) > (m1 + m2))
			    {
			      this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
			      this->M1Value[this->NbrInteractionFactors] = m1;
			      this->M2Value[this->NbrInteractionFactors] = m2;
			      this->M3Value[this->NbrInteractionFactors] = m3;
			      ++this->NbrInteractionFactors;
			    }
			  else
			    if ((2 * m3) == (m1 + m2))
			      {
				this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * Factor * TmpCoefficient[Pos];
				this->M1Value[this->NbrInteractionFactors] = m1;
				this->M2Value[this->NbrInteractionFactors] = m2;
				this->M3Value[this->NbrInteractionFactors] = m3;
				++this->NbrInteractionFactors;
			      }
			}
		      ++Pos;
		    }
		}	
	      Lim = 2 * m1;
	      if (Lim > this->LzMax)
		Lim = this->LzMax;
	      Min = 2 * m1 - this->LzMax;
	      if (Min < 0)
		Min = 0;
	      for (int m3 = Min; m3 <= Lim; ++m3)
		{
		  if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      if (m3 > m1)
			{
			  this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * Factor * TmpCoefficient[Pos];
			  this->M1Value[this->NbrInteractionFactors] = m1;
			  this->M2Value[this->NbrInteractionFactors] = m1;
			  this->M3Value[this->NbrInteractionFactors] = m3;
			  ++this->NbrInteractionFactors;
			}
		      else
			if (m3 == m1)
			  {
			    this->InteractionFactors[this->NbrInteractionFactors] = 0.25 * Factor * TmpCoefficient[Pos];
			    this->M1Value[this->NbrInteractionFactors] = m1;
			    this->M2Value[this->NbrInteractionFactors] = m1;
			    this->M3Value[this->NbrInteractionFactors] = m3;
			    ++this->NbrInteractionFactors;
			  }
		    }
		  ++Pos;
		}
	    }
	}
    }
}

// compute all projector coefficient associated to a given relative angular momentum between k particles
//
// relativeMomentum = value of twice the relative angular momentum between the k particles 
// nbrIndices = number of indices per set
// indices = array that contains all possible sets of indices (size of the array is nbrIndices * nbrIndexSets)
// nbrIndexSets = number of sets

double* ParticleOnSphereNBodyHardCoreHamiltonian::ComputeProjectorCoefficients(int relativeMomentum, int nbrIndices, int* indices, int nbrIndexSets)
{
  double* TmpCoefficients = new double [nbrIndexSets];
  int JValue = (nbrIndices * this->LzMax) - relativeMomentum;
  switch (nbrIndices)
    {
    case 2:
      {
	ClebschGordanCoefficients Clebsh (this->LzMax, this->LzMax);
	for (int i = 0; i < nbrIndexSets; ++i)
	  {
	    TmpCoefficients[i] = Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1) - this->LzMax), JValue);
	    indices += 2;
	  }
      }
      break;
    case 3:
      {
	int MaxJ = 2 * this->LzMax;
	if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	  MaxJ -= 2;
	ClebschGordanCoefficients Clebsh (this->LzMax, this->LzMax);
	ClebschGordanCoefficients* ClebshArray = new ClebschGordanCoefficients[MaxJ + 1];
	int MinJ = JValue - this->LzMax;
	if (MinJ < 0)
	  MinJ = 0;
	for (int j = MaxJ; j >= MinJ; j -= 4)
	  ClebshArray[j] = ClebschGordanCoefficients(j, this->LzMax);
	for (int i = 0; i < nbrIndexSets; ++i)
	  {
	    double Tmp = 0.0;
	    int Sum = ((indices[0] + indices[1]) << 1)  - (2 * this->LzMax);
	    int TmpMinJ = MinJ;
	    if (TmpMinJ < abs(Sum))
	      TmpMinJ = abs(Sum);
	    for (int j = MaxJ; j >= TmpMinJ; j -= 4)
	      {
		Tmp += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) * 
			ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue)); 
	      }
	    Sum = ((indices[1] + indices[2]) << 1)  - (2 * this->LzMax);
	    TmpMinJ = MinJ;
	    if (TmpMinJ < abs(Sum))
	      TmpMinJ = abs(Sum);
	    for (int j = MaxJ; j >= TmpMinJ; j -= 4)
	      {
		Tmp += (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) * 
			ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue)); 
	      }
	    Sum = ((indices[2] + indices[0]) << 1)  - (2 * this->LzMax);
	    TmpMinJ = MinJ;
	    if (TmpMinJ < abs(Sum))
	      TmpMinJ = abs(Sum);
	    for (int j = MaxJ; j >= TmpMinJ; j -= 4)
	      {
		Tmp += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) * 
			ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue)); 
	      }
	    TmpCoefficients[i] = Tmp;
	    indices += 3;
	  }
	delete[] ClebshArray;
      }
      break;
    case 4:
      {
	int TotalMaxJ = 3 * this->LzMax - 2;
	ClebschGordanCoefficients Clebsh (this->LzMax, this->LzMax);
	ClebschGordanCoefficients* ClebshArray = new ClebschGordanCoefficients[TotalMaxJ + 1];
	int MinJ = JValue - this->LzMax;
	if (MinJ < 0)
	  MinJ = 0;
	for (int j = TotalMaxJ; j >= 0; --j)
	  ClebshArray[j] = ClebschGordanCoefficients(j, this->LzMax);
	for (int i = 0; i < nbrIndexSets; ++i)
	  {
	    double Tmp = 0.0;
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[0], indices[1], indices[2], indices[3], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[0], indices[1], indices[3], indices[2], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[0], indices[2], indices[1], indices[3], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[0], indices[2], indices[3], indices[1], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[0], indices[3], indices[2], indices[1], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[0], indices[3], indices[1], indices[2], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[1], indices[2], indices[0], indices[3], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[1], indices[2], indices[3], indices[0], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[1], indices[3], indices[2], indices[0], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[1], indices[3], indices[0], indices[2], JValue, ClebshArray);	  
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[2], indices[3], indices[0], indices[1], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[2], indices[3], indices[1], indices[0], JValue, ClebshArray);
	    TmpCoefficients[i] = Tmp;
	    indices += 4;
	  }
	delete[] ClebshArray;
      }
      break;
    default:
      {
      }
    }
  return TmpCoefficients;
}

// compute a given projector coefficient for the 4-body interaction 
//
// m1 = first index
// m2 = second index
// m3 = third inde
// m4 = fourth index
// jValue = total angular momentum
// minJ = minimum angular momentum that can be reach by three particles
// return value = corresponding projector coefficient

double ParticleOnSphereNBodyHardCoreHamiltonian::ComputeProjectorCoefficients4Body(int m1, int m2, int m3, int m4, int jValue, 
										   ClebschGordanCoefficients* clebshArray)
{
  double Tmp = 0.0;
  double Tmp2 = 0.0;
  int TmpMinJ2;
  int Sum1 = ((m1 + m2) << 1)  - (2 * this->LzMax);
  int Sum2 = Sum1 + (m3 << 1)  - this->LzMax;
  int TmpMinJ = abs(Sum1);
  int MaxJ = 2 * this->LzMax - 2;
  for (int j1 = MaxJ; j1 >= TmpMinJ; j1 -= 4)
    {
      Tmp2 = clebshArray[this->LzMax].GetCoefficient(((m1 << 1) - this->LzMax), ((m2 << 1)- this->LzMax), j1);
      TmpMinJ2 = abs(jValue - this->LzMax);
      if (TmpMinJ2 < abs(Sum2))
	TmpMinJ2 = abs(Sum2);
      for (int j2 = j1 + this->LzMax; j2 >= TmpMinJ2; j2 -= 2)
	{
	  Tmp += Tmp2 * (clebshArray[j1].GetCoefficient(Sum1, ((m3 << 1) - this->LzMax), j2) * 
			 clebshArray[j2].GetCoefficient(Sum2, ((m4 << 1) - this->LzMax), jValue)); 
	}
    }
  return Tmp;
}
