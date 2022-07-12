////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated to particles on a disk with         //
//                         n-body hard core interaction                       //
//                                                                            //
//                        last modification : 04/09/2004                      //
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
#include "Hamiltonian/ParticleOnDiskNBodyHardCoreHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// architecture = architecture to use for precalculation
// nbrBody = number of particle that interact simultaneously through the hard core interaction
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnDiskNBodyHardCoreHamiltonian::ParticleOnDiskNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int nbrBody,
									       AbstractArchitecture* architecture, long memory, 
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
  this->DiskStorageFlag = false;
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
  this->L2Operator = 0;
}

// destructor
//

ParticleOnDiskNBodyHardCoreHamiltonian::~ParticleOnDiskNBodyHardCoreHamiltonian()
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

void ParticleOnDiskNBodyHardCoreHamiltonian::EvaluateInteractionFactors()
{
  double* TmpNormalizationCoeffients = new double[this->NbrLzValue];
  double TmpFactor = ((double) this->NbrLzValue) / (4.0 * M_PI);
  double TmpBinomial = 1.0;
  TmpNormalizationCoeffients[0] = sqrt (TmpBinomial * TmpFactor);
  for (int i = 1; i < this->NbrLzValue; ++i)
    {
      TmpBinomial *= this->NbrLzValue - ((double) i) + 1.0;
      TmpBinomial /= ((double) i);
      TmpNormalizationCoeffients[i] = sqrt (TmpBinomial * TmpFactor);
    }
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      // useless part (trivially equal to zero for fermions), to be used as an example for other interactions
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    double SumCoefficient = 0.0;
	    double Coefficient;
	    GetAllSkewSymmetricIndices(this->NbrLzValue, k, this->NbrSortedIndicesPerSum[k], this->SortedIndicesPerSum[k]);
	    this->MaxSumIndices[k] = (((this->NbrLzValue - 1) * this->NbrLzValue) - ((k - 1) * (k - 2)))/ 2;
	    this->MinSumIndices[k] = (k * (k - 1)) / 2;
	    this->NBodyInteractionFactors[k] = new double* [this->MaxSumIndices[k] + 1];
	    int Lim;
	    for (int MinSum = this->MinSumIndices[k]; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	      {
		Lim = this->NbrSortedIndicesPerSum[k][MinSum];
		this->NBodyInteractionFactors[k][MinSum] = new double [Lim];
		double* TmpNBodyInteractionFactors = this->NBodyInteractionFactors[k][MinSum];		  
		int* TmpMIndices = this->SortedIndicesPerSum[k][MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    Coefficient = SumCoefficient;
		    for (int l = 0; l < k; ++l)
		      Coefficient *= TmpNormalizationCoeffients[TmpMIndices[l]];		    
		    TmpNBodyInteractionFactors[i] = Coefficient;
		    TmpMIndices += k;
		  }
	      }
	  }
    }
  else
    {
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    this->MinSumIndices[k] = 0;
	    this->MaxSumIndices[k] = this->LzMax * k;
	    double** SortedIndicesPerSumSymmetryFactor;
	    GetAllSymmetricIndices(this->NbrLzValue, k, this->NbrSortedIndicesPerSum[k], this->SortedIndicesPerSum[k],
				   SortedIndicesPerSumSymmetryFactor);
	    this->NBodyInteractionFactors[k] = new double* [this->MaxSumIndices[k] + 1];
	    int Lim;
	    double Factor = 1.0;
	    for (int i = 1; i < k; ++i)
	      Factor *= 8.0 * M_PI;
	    Factor = 1.0 / Factor;
	    FactorialCoefficient Coef;


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
		Lim = this->NbrSortedIndicesPerSum[k][MinSum];
		double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
		double* TmpNBodyInteractionFactors = new double [Lim];
		this->NBodyInteractionFactors[k][MinSum] = TmpNBodyInteractionFactors;
		int* TmpNIndices2 = this->SortedIndicesPerSum[k][MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    Coef.SetToOne();
		    Coef.FactorialMultiply(MinSum);
		    for (int l = 0; l < k; ++l)
		      {
			Coef.FactorialDivide(TmpNIndices2[l]);
		      }
		    Coef.PowerNDivide(k, MinSum);
		    TmpNBodyInteractionFactors[i] = sqrt(Coef.GetNumericalValue() * Factor) * TmpSymmetryFactors[i];
		    TmpNIndices2 += k;
		  }
		TmpNIndices2 = this->SortedIndicesPerSum[k][MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    this->NbrMIndices[k][TmpNbrNIndices] = Lim;		    
		    this->MIndices[k][TmpNbrNIndices] = new int [Lim * k];
		    this->MNNBodyInteractionFactors[k][TmpNbrNIndices] = new double [Lim];
		    int* TmpMIndices = this->MIndices[k][TmpNbrNIndices];
		    int* TmpMIndices2 = this->SortedIndicesPerSum[k][MinSum];
		    double* TmpInteraction = this->MNNBodyInteractionFactors[k][TmpNbrNIndices];
		    double Coefficient = TmpNBodyInteractionFactors[i] * this->NBodyInteractionWeightFactors[k];
		    for (int j = 0; j < Lim; ++j)
		      {
			for (int l = 0; l < k; ++l)
			  {
			    (*TmpMIndices) = (*TmpMIndices2);			
			    ++TmpMIndices;
			    ++TmpMIndices2;
			  }			
			TmpInteraction[j] = Coefficient *  TmpNBodyInteractionFactors[j];
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
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	      {
		delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	      }
	    delete[] SortedIndicesPerSumSymmetryFactor;
	  }
    }
}

