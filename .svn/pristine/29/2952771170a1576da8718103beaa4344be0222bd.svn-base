////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                          generic 4-body interaction                        //
//                                                                            //
//                        last modification : 12/08/2009                      //
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
#include "Hamiltonian/ParticleOnSphereGenericFiveBodyHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"

  
#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//

ParticleOnSphereGenericFiveBodyHamiltonian::ParticleOnSphereGenericFiveBodyHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// fiveBodyPseudoPotential = array with the five-body pseudo-potentials sorted with respect to the relative angular momentum, 
//                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
// maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FiveBodyPseudoPotential
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereGenericFiveBodyHamiltonian::ParticleOnSphereGenericFiveBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
										       double* fiveBodyPseudoPotential, int maxRelativeAngularMomentum, double l2Factor, 
										       AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
										       char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = false;
  this->MaxNBody = 5;

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

  this->MaxRelativeAngularMomentum = maxRelativeAngularMomentum;
  this->NbrFiveBodyPseudoPotential = maxRelativeAngularMomentum + 1;
  this->FourBodyPseudoPotential = 0;
  this->FiveBodyPseudoPotential = new double[this->NbrFiveBodyPseudoPotential];
  for (int i = 0; i < this->NbrFiveBodyPseudoPotential; ++i)
    this->FiveBodyPseudoPotential[i] = fiveBodyPseudoPotential[i];

  this->NBodyFlags[5] = true;
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

// constructor from default datas with a full four-body interaction
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// fiveBodyPseudoPotential = array with the five-body pseudo-potentials sorted with respect to the relative angular momentum, 
// maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FiveBodyPseudoPotential
// fourBodyPseudoPotential = array with the four-body pseudo-potentials sorted with respect to the relative angular momentum, 
// fourBodyMaxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FourBodyPseudoPotential
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereGenericFiveBodyHamiltonian::ParticleOnSphereGenericFiveBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
										       double* fiveBodyPseudoPotential, int maxRelativeAngularMomentum, 
										       double* fourBodyPseudoPotential, int fourBodyMaxRelativeAngularMomentum, 
										       double l2Factor, AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
										       char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = false;
  this->MaxNBody = 5;

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

  this->MaxRelativeAngularMomentum = maxRelativeAngularMomentum;
  this->NbrFiveBodyPseudoPotential = maxRelativeAngularMomentum + 1;
  this->NbrFourBodyPseudoPotential = fourBodyMaxRelativeAngularMomentum + 1;
  this->FourBodyPseudoPotential = new double[this->NbrFourBodyPseudoPotential];
  for (int i = 0; i < this->NbrFourBodyPseudoPotential; ++i)
    {
      this->FourBodyPseudoPotential[i] = fourBodyPseudoPotential[i];
    }
  this->FiveBodyPseudoPotential = new double[this->NbrFiveBodyPseudoPotential];
  for (int i = 0; i < this->NbrFiveBodyPseudoPotential; ++i)
    {
      this->FiveBodyPseudoPotential[i] = fiveBodyPseudoPotential[i];
    }
  this->NBodyFlags[4] = true;
  this->NBodyFlags[5] = true;
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

// constructor from datas with a fully-defined two body interaction
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// fiveBodyPseudoPotential = array with the five-body pseudo-potentials sorted with respect to the relative angular momentum, 
//                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
// maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FiveBodyPseudoPotential
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereGenericFiveBodyHamiltonian::ParticleOnSphereGenericFiveBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
											 double* fiveBodyPseudoPotential, int maxRelativeAngularMomentum,
											 double l2Factor, double* pseudoPotential, 
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
  this->MaxNBody = 5;
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

  this->PseudoPotential = new double [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->PseudoPotential[i] = pseudoPotential[this->LzMax - i];

  this->MaxRelativeAngularMomentum = maxRelativeAngularMomentum;
  this->NbrFiveBodyPseudoPotential = maxRelativeAngularMomentum;
  this->FourBodyPseudoPotential = 0;
  this->FiveBodyPseudoPotential = new double[this->NbrFiveBodyPseudoPotential + 1];
  for (int i = 0; i <= this->NbrFiveBodyPseudoPotential; ++i)
    this->FiveBodyPseudoPotential[i] = fiveBodyPseudoPotential[i];


  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 1;
      this->MaxSumIndices[k] = 0;      
      this->NBodyFlags[k] = false;
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	this->NBodySign[k] = -1.0;
    }
  this->NBodyFlags[5] = true;
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

ParticleOnSphereGenericFiveBodyHamiltonian::~ParticleOnSphereGenericFiveBodyHamiltonian()
{
  if (this->FiveBodyPseudoPotential != 0)
    delete[] this->FiveBodyPseudoPotential;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* ParticleOnSphereGenericFiveBodyHamiltonian::Clone ()
{
  return 0;
}


// evaluate all interaction factors
//   

void ParticleOnSphereGenericFiveBodyHamiltonian::EvaluateInteractionFactors()
{
  this->Evaluate5BodyInteractionFactors();
  if (this->FourBodyPseudoPotential != 0)
    {
      this->Evaluate4BodyInteractionFactors();      
    }
}


// evaluate all interaction factors for the five body interaction part
//   

void ParticleOnSphereGenericFiveBodyHamiltonian::Evaluate5BodyInteractionFactors()
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
      double Coefficient;
      GetAllSkewSymmetricIndices(this->NbrLzValue, 5, this->NbrSortedIndicesPerSum[5], this->SortedIndicesPerSum[5]);
      this->MaxSumIndices[5] = (this->LzMax * 5) - 9;
      this->MinSumIndices[5] = 10;
      double* TmpInteractionCoeffients = new double[this->MaxSumIndices[5] + 1];
      TmpInteractionCoeffients[0] = 1.0;
      TmpInteractionCoeffients[1] = 1.0;
      for (int i = 2; i <= this->MaxSumIndices[5]; ++i)
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
      Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[5]) + 1.0);
      double Radius = 2.0 / ((double) this->LzMax);
      for (int i = 2; i <= 5; ++i)
	{
	  Coefficient *= (double) (i * i);	  
	  Coefficient *= Radius;
	}
      for (int i = 0; i <= this->MaxSumIndices[5]; ++i)
	TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);
      
      long TmpNbrNIndices = 0;
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[5]; ++MinSum)
	TmpNbrNIndices += this->NbrSortedIndicesPerSum[5][MinSum];
      this->NbrNIndices[5] = TmpNbrNIndices;
      this->NIndices[5] = new int[TmpNbrNIndices * 5];
      this->NbrMIndices[5] = new long[TmpNbrNIndices];
      this->MIndices[5] = new int*[TmpNbrNIndices];
      this->MNNBodyInteractionFactors[5] = new double* [TmpNbrNIndices];
      TmpNbrNIndices = 0;	 
      int* TmpNIndices = this->NIndices[5];
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[5]; ++MinSum)
	{
	  int Lim = this->NbrSortedIndicesPerSum[5][MinSum];
	  if (Lim > 0)
	    {
	      int* TmpNIndices2 = this->SortedIndicesPerSum[5][MinSum];
	      int TmpMaxRealtiveMonentum = this->NbrFiveBodyPseudoPotential - 1;
	      int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2] + TmpNIndices2[3] + TmpNIndices2[5];
	      while (((5 * this->LzMax) - TmpMaxRealtiveMonentum)  < TmpSum)
		--TmpMaxRealtiveMonentum;
	      double** TmpProjectorCoefficients = new double* [TmpMaxRealtiveMonentum + 1];
	      if (this->FiveBodyPseudoPotential[5] != 0.0)
		TmpProjectorCoefficients[5] = this->Compute5BodyCoefficients(12, TmpNIndices2, Lim);
	      for (int i = 5; i <= TmpMaxRealtiveMonentum; ++i)  
		if (this->FiveBodyPseudoPotential[i] != 0.0)
		  TmpProjectorCoefficients[i] = this->Compute5BodyCoefficients(2 * i, TmpNIndices2, Lim);
	      for (int i = 0; i < Lim; ++i)
		{
		  this->NbrMIndices[5][TmpNbrNIndices] = Lim;		    
		  this->MIndices[5][TmpNbrNIndices] = new int [Lim * 3];
		  this->MNNBodyInteractionFactors[5][TmpNbrNIndices] = new double [Lim];
		  int* TmpMIndices = this->MIndices[5][TmpNbrNIndices];
		  int* TmpMIndices2 = this->SortedIndicesPerSum[5][MinSum];
		  double* TmpInteraction = this->MNNBodyInteractionFactors[5][TmpNbrNIndices];
		  for (int j = 0; j < Lim; ++j)
		    {
		      for (int l = 0; l < 5; ++l)
			{
			  (*TmpMIndices) = (*TmpMIndices2);			
			  ++TmpMIndices;
			  ++TmpMIndices2;
			}			
		      double& TmpInteraction2 = TmpInteraction[j];
		      TmpInteraction2 = 0.0;
		      if (this->FiveBodyPseudoPotential[5] != 0.0)
			TmpInteraction2 += this->FiveBodyPseudoPotential[5] * TmpProjectorCoefficients[5][i] * TmpProjectorCoefficients[5][j];
		      for (int k = 5; k <= TmpMaxRealtiveMonentum; ++k)  
			if (this->FiveBodyPseudoPotential[k] != 0.0)
			  TmpInteraction2 += this->FiveBodyPseudoPotential[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
		    }
		  for (int j = 0; j < 5; ++j)
		    {
		      (*TmpNIndices) = (*TmpNIndices2);			
		      ++TmpNIndices;
		      ++TmpNIndices2;
		    }
		  ++TmpNbrNIndices;
		}
	      if (this->FiveBodyPseudoPotential[5] != 0.0)
		delete[] TmpProjectorCoefficients[5];
	      for (int i = 5; i <= TmpMaxRealtiveMonentum; ++i)  
		if (this->FiveBodyPseudoPotential[i] != 0.0)
		  delete[] TmpProjectorCoefficients[i];
	      delete[] TmpProjectorCoefficients;		
	    }
	}
      delete[] TmpInteractionCoeffients;
    }
  else
    {
      this->MinSumIndices[5] = 0;
      this->MaxSumIndices[5] = this->LzMax * 5;
      double* TmpInteractionCoeffients = new double[this->MaxSumIndices[5] + 1];
      double Coefficient;
      TmpInteractionCoeffients[0] = 1.0;
      TmpInteractionCoeffients[1] = 1.0;
      for (int i = 2; i <= this->MaxSumIndices[5]; ++i)
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
      Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[5]) + 1.0);
      double Radius = 2.0 / ((double) this->LzMax);
      for (int i = 2; i <= 5; ++i)
	{
	  Coefficient *= (double) (i * i);	  
	  Coefficient *= Radius;
	}
      for (int i = 0; i <= this->MaxSumIndices[5]; ++i)
	TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);
      
      double** SortedIndicesPerSumSymmetryFactor;
      GetAllSymmetricIndices(this->NbrLzValue, 5, this->NbrSortedIndicesPerSum[5], this->SortedIndicesPerSum[5],
			     SortedIndicesPerSumSymmetryFactor);
      
      
      long TmpNbrNIndices = 0;
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[5]; ++MinSum)
	TmpNbrNIndices += this->NbrSortedIndicesPerSum[5][MinSum];
      this->NbrNIndices[5] = TmpNbrNIndices;
      this->NIndices[5] = new int[TmpNbrNIndices * 5];
      this->NbrMIndices[5] = new long[TmpNbrNIndices];
      this->MIndices[5] = new int*[TmpNbrNIndices];
      this->MNNBodyInteractionFactors[5] = new double* [TmpNbrNIndices];
      TmpNbrNIndices = 0;	 
      int* TmpNIndices = this->NIndices[5];
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[5]; ++MinSum)
	{
	  int Lim = this->NbrSortedIndicesPerSum[5][MinSum];
	  double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
	  int* TmpNIndices2 = this->SortedIndicesPerSum[5][MinSum];
	  int TmpMaxRealtiveMonentum = this->NbrFiveBodyPseudoPotential - 1;
	  int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2] + TmpNIndices2[3] + TmpNIndices2[4];
	  while (((5 * this->LzMax) - TmpMaxRealtiveMonentum)  < TmpSum)
	    --TmpMaxRealtiveMonentum;
	  double** TmpProjectorCoefficients = new double* [ TmpMaxRealtiveMonentum+ 1];
	  double** TmpProjectorCoefficients2 = new double* [ TmpMaxRealtiveMonentum+ 1];
	  if (this->FiveBodyPseudoPotential[0] != 0.0)
	    TmpProjectorCoefficients[0] = this->Compute5BodyCoefficients(0, TmpNIndices2, Lim);
	  for (int i = 1; i <= TmpMaxRealtiveMonentum; ++i)  
	    if (this->FiveBodyPseudoPotential[i] != 0.0)
	      TmpProjectorCoefficients[i] = this->Compute5BodyCoefficients(2 * i, TmpNIndices2, Lim);	    
	  if ((TmpMaxRealtiveMonentum >= 4) && (this->FiveBodyPseudoPotential[4] != 0.0))
	    {
	      int MaxClosing[4];
	      MaxClosing[0] = 0;
	      MaxClosing[1] = 0;
	      MaxClosing[2] = 4;
	      TmpProjectorCoefficients2[4] = this->Compute5BodyCoefficientsWithDirection(8, TmpNIndices2, Lim, MaxClosing);	   
	    }
	  if ((TmpMaxRealtiveMonentum >= 5) && (this->FiveBodyPseudoPotential[5] != 0.0))
	    {
	      int MaxClosing[4];
	      MaxClosing[0] = 0;
	      MaxClosing[1] = 0;
	      MaxClosing[2] = 4;
	      TmpProjectorCoefficients2[5] = this->Compute5BodyCoefficientsWithDirection(10, TmpNIndices2, Lim, MaxClosing);	   
	    }

	  for (int i = 0; i < Lim; ++i)
	    {
	      this->NbrMIndices[5][TmpNbrNIndices] = Lim;		    
	      this->MIndices[5][TmpNbrNIndices] = new int [Lim * 5];
	      this->MNNBodyInteractionFactors[5][TmpNbrNIndices] = new double [Lim];
	      int* TmpMIndices = this->MIndices[5][TmpNbrNIndices];
	      int* TmpMIndices2 = this->SortedIndicesPerSum[5][MinSum];
	      double* TmpInteraction = this->MNNBodyInteractionFactors[5][TmpNbrNIndices];
	      for (int j = 0; j < Lim; ++j)
		{
		  double Coefficient2 = TmpSymmetryFactors[j];
		  for (int l = 0; l < 5; ++l)
		    {
		      Coefficient2 *= TmpNormalizationCoeffients[(*TmpMIndices2)];		    
		      (*TmpMIndices) = (*TmpMIndices2);			
		      ++TmpMIndices;
		      ++TmpMIndices2;
		    }			
		  double& TmpInteraction2 = TmpInteraction[j];
		  TmpInteraction2 = 0.0;
		  if (this->FiveBodyPseudoPotential[0] != 0.0)
		    TmpInteraction2 += this->FiveBodyPseudoPotential[0] * TmpProjectorCoefficients[0][i] * TmpProjectorCoefficients[0][j];
		  for (int k = 2; k <= TmpMaxRealtiveMonentum; ++k)  
		    if (this->FiveBodyPseudoPotential[k] != 0.0)
		      TmpInteraction2 += this->FiveBodyPseudoPotential[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
		  for (int k = 4; k <= TmpMaxRealtiveMonentum; ++k)  
		    if (this->FiveBodyPseudoPotential[k] != 0.0)
		      TmpInteraction2 += this->FiveBodyPseudoPotential[k] * TmpProjectorCoefficients2[k][i] * TmpProjectorCoefficients2[k][j];
		  TmpInteraction2 *= TmpSymmetryFactors[i] * TmpSymmetryFactors[j];
		}
	      for (int j = 0; j < 5; ++j)
		{
		  (*TmpNIndices) = (*TmpNIndices2);			
		  ++TmpNIndices;
		  ++TmpNIndices2;
		}
	      ++TmpNbrNIndices;
	    }		
	  if (this->FiveBodyPseudoPotential[0] != 0.0)
	    delete[] TmpProjectorCoefficients[0];
	  for (int i = 2; i <= TmpMaxRealtiveMonentum; ++i)  
	    if (this->FiveBodyPseudoPotential[i] != 0.0)
	      delete[] TmpProjectorCoefficients[i];

	  delete[] TmpProjectorCoefficients;		
	}
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[5]; ++MinSum)
	{
	  delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	}
      delete[] SortedIndicesPerSumSymmetryFactor;
      delete[] TmpInteractionCoeffients;
    }
  delete[] TmpNormalizationCoeffients;
}

// compute all projector coefficient associated to a given relative angular momentum between 5 particles
//
// relativeMomentum = value of twice the relative angular momentum between the 4 particles
// indices = array that contains all possible sets of indices (size of the array is 4 * nbrIndexSets)
// nbrIndexSets = number of sets

double* ParticleOnSphereGenericFiveBodyHamiltonian::Compute5BodyCoefficients(int relativeMomentum, int* indices, int nbrIndexSets)
{
  double* TmpCoefficients = new double [nbrIndexSets];
  int JValue = (5 * this->LzMax) - relativeMomentum;

  int MaxJ = 4 * this->LzMax;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    MaxJ -= 2;
  ClebschGordanCoefficients Clebsh (this->LzMax, this->LzMax);
  ClebschGordanCoefficients* ClebshArray = 0;
  ClebshArray = new ClebschGordanCoefficients[MaxJ + 1];
  for (int j = MaxJ; j >= 0; --j)
    ClebshArray[j] = ClebschGordanCoefficients(j, this->LzMax);

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    for (int i = 0; i < nbrIndexSets; ++i)
      {
	double Tmp = 0.0;
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[2], indices[3], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[2], indices[4], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[3], indices[2], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[3], indices[4], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[4], indices[2], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[4], indices[3], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[1], indices[3], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[1], indices[4], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[3], indices[1], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[3], indices[4], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[4], indices[1], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[4], indices[3], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[1], indices[2], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[1], indices[4], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[2], indices[1], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[2], indices[4], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[4], indices[1], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[4], indices[2], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[1], indices[2], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[1], indices[3], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[2], indices[1], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[2], indices[3], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[3], indices[1], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[3], indices[2], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[2], indices[3], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[2], indices[4], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[3], indices[2], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[3], indices[4], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[4], indices[2], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[4], indices[3], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[0], indices[3], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[0], indices[4], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[3], indices[0], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[3], indices[4], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[4], indices[0], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[4], indices[3], indices[0], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[0], indices[2], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[0], indices[4], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[2], indices[0], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[2], indices[4], indices[0], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[4], indices[0], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[4], indices[2], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[0], indices[2], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[0], indices[3], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[2], indices[0], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[2], indices[3], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[3], indices[0], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[3], indices[2], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[1], indices[3], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[1], indices[4], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[3], indices[1], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[3], indices[4], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[4], indices[1], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[4], indices[3], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[0], indices[3], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[0], indices[4], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[3], indices[0], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[3], indices[4], indices[0], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[4], indices[0], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[4], indices[3], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[0], indices[1], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[0], indices[4], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[1], indices[0], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[1], indices[4], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[4], indices[0], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[4], indices[1], indices[0], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[0], indices[1], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[0], indices[3], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[1], indices[0], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[1], indices[3], indices[0], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[3], indices[0], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[3], indices[1], indices[0], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[1], indices[2], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[1], indices[4], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[2], indices[1], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[2], indices[4], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[4], indices[1], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[4], indices[2], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[0], indices[2], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[0], indices[4], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[2], indices[0], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[2], indices[4], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[4], indices[0], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[4], indices[2], indices[0], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[0], indices[1], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[0], indices[4], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[1], indices[0], indices[4], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[1], indices[4], indices[0], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[4], indices[0], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[4], indices[1], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[0], indices[1], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[0], indices[2], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[1], indices[0], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[1], indices[2], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[2], indices[0], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[2], indices[1], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[1], indices[2], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[1], indices[3], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[2], indices[1], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[2], indices[3], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[3], indices[1], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[3], indices[2], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[0], indices[2], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[0], indices[3], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[2], indices[0], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[2], indices[3], indices[0], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[3], indices[0], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[3], indices[2], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[0], indices[1], indices[3], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[0], indices[3], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[1], indices[0], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[1], indices[3], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[3], indices[0], indices[1], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[3], indices[1], indices[0], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[0], indices[1], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[0], indices[2], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[1], indices[0], indices[2], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[1], indices[2], indices[0], JValue, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[2], indices[0], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[2], indices[1], indices[0], JValue, ClebshArray);
	TmpCoefficients[i] = Tmp;
	indices += 5;
      }
  else
    for (int i = 0; i < nbrIndexSets; ++i)
      {
	double Tmp = 0.0;
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[2], indices[3], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[2], indices[4], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[3], indices[2], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[3], indices[4], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[4], indices[2], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[4], indices[3], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[1], indices[3], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[1], indices[4], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[3], indices[1], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[3], indices[4], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[4], indices[1], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[2], indices[4], indices[3], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[1], indices[2], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[1], indices[4], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[2], indices[1], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[2], indices[4], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[4], indices[1], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[3], indices[4], indices[2], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[1], indices[2], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[1], indices[3], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[2], indices[1], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[2], indices[3], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[3], indices[1], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[4], indices[3], indices[2], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[2], indices[3], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[2], indices[4], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[3], indices[2], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[3], indices[4], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[4], indices[2], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[0], indices[4], indices[3], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[0], indices[3], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[0], indices[4], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[3], indices[0], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[3], indices[4], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[4], indices[0], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[2], indices[4], indices[3], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[0], indices[2], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[0], indices[4], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[2], indices[0], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[2], indices[4], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[4], indices[0], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[3], indices[4], indices[2], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[0], indices[2], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[0], indices[3], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[2], indices[0], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[2], indices[3], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[3], indices[0], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[1], indices[4], indices[3], indices[2], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[1], indices[3], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[1], indices[4], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[3], indices[1], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[3], indices[4], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[4], indices[1], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[0], indices[4], indices[3], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[0], indices[3], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[0], indices[4], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[3], indices[0], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[3], indices[4], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[4], indices[0], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[1], indices[4], indices[3], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[0], indices[1], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[0], indices[4], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[1], indices[0], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[1], indices[4], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[4], indices[0], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[3], indices[4], indices[1], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[0], indices[1], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[0], indices[3], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[1], indices[0], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[1], indices[3], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[3], indices[0], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[2], indices[4], indices[3], indices[1], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[1], indices[2], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[1], indices[4], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[2], indices[1], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[2], indices[4], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[4], indices[1], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[0], indices[4], indices[2], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[0], indices[2], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[0], indices[4], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[2], indices[0], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[2], indices[4], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[4], indices[0], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[1], indices[4], indices[2], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[0], indices[1], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[0], indices[4], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[1], indices[0], indices[4], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[1], indices[4], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[4], indices[0], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[2], indices[4], indices[1], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[0], indices[1], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[0], indices[2], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[1], indices[0], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[1], indices[2], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[2], indices[0], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[3], indices[4], indices[2], indices[1], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[1], indices[2], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[1], indices[3], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[2], indices[1], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[2], indices[3], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[3], indices[1], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[0], indices[3], indices[2], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[0], indices[2], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[0], indices[3], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[2], indices[0], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[2], indices[3], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[3], indices[0], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[1], indices[3], indices[2], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[0], indices[1], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[0], indices[3], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[1], indices[0], indices[3], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[1], indices[3], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[3], indices[0], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[2], indices[3], indices[1], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[0], indices[1], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[0], indices[2], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[1], indices[0], indices[2], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[1], indices[2], indices[0], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[2], indices[0], indices[1], JValue, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5Body(indices[4], indices[3], indices[2], indices[1], indices[0], JValue, ClebshArray);
	TmpCoefficients[i] = Tmp;
	indices += 5;
      }
  delete[] ClebshArray;
  return TmpCoefficients;
}

// compute all projector coefficient associated to a given relative angular momentum between 5 particles in a given direction
//
// relativeMomentum = value of twice the relative angular momentum between the 4 particles
// indices = array that contains all possible sets of indices (size of the array is 4 * nbrIndexSets)
// nbrIndexSets = number of sets
// maxClosing = array that gives the maximum angular momentum when a particle approach the cluster of n+1 particles  (n being the index of MaxClosing)

double* ParticleOnSphereGenericFiveBodyHamiltonian::Compute5BodyCoefficientsWithDirection(int relativeMomentum, int* indices, int nbrIndexSets, int* maxClosing)
{
  double* TmpCoefficients = new double [nbrIndexSets];
  int JValue = (5 * this->LzMax) - relativeMomentum;

  int MaxJ = 4 * this->LzMax;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    MaxJ -= 2;
  ClebschGordanCoefficients Clebsh (this->LzMax, this->LzMax);
  ClebschGordanCoefficients* ClebshArray = 0;
  ClebshArray = new ClebschGordanCoefficients[MaxJ + 1];
  for (int j = MaxJ; j >= 0; --j)
    ClebshArray[j] = ClebschGordanCoefficients(j, this->LzMax);

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    for (int i = 0; i < nbrIndexSets; ++i)
      {
	double Tmp = 0.0;
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[2], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[2], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[3], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[3], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[4], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[4], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[1], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[1], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[3], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[3], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[4], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[4], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[1], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[1], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[2], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[2], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[4], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[4], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[1], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[1], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[2], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[2], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[3], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[3], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[2], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[2], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[3], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[3], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[4], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[4], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[0], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[0], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[3], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[3], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[4], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[4], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[0], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[0], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[2], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[2], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[4], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[4], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[0], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[0], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[2], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[2], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[3], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[3], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[1], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[1], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[3], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[3], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[4], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[4], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[0], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[0], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[3], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[3], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[4], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[4], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[0], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[0], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[1], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[1], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[4], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[4], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[0], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[0], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[1], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[1], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[3], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[3], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[1], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[1], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[2], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[2], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[4], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[4], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[0], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[0], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[2], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[2], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[4], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[4], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[0], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[0], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[1], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[1], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[4], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[4], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[0], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[0], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[1], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[1], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[2], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[2], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[1], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[1], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[2], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[2], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[3], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[3], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[0], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[0], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[2], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[2], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[3], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[3], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[0], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[0], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[1], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[1], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[3], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[3], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[0], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[0], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[1], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[1], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp -= this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[2], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[2], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	TmpCoefficients[i] = Tmp;
	indices += 5;
      }
  else
    for (int i = 0; i < nbrIndexSets; ++i)
      {
	double Tmp = 0.0;
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[2], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[2], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[3], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[3], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[4], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[1], indices[4], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[1], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[1], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[3], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[3], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[4], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[2], indices[4], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[1], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[1], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[2], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[2], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[4], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[3], indices[4], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[1], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[1], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[2], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[2], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[3], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[0], indices[4], indices[3], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[2], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[2], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[3], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[3], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[4], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[0], indices[4], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[0], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[0], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[3], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[3], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[4], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[2], indices[4], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[0], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[0], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[2], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[2], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[4], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[3], indices[4], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[0], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[0], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[2], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[2], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[3], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[1], indices[4], indices[3], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[1], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[1], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[3], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[3], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[4], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[0], indices[4], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[0], indices[3], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[0], indices[4], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[3], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[3], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[4], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[1], indices[4], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[0], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[0], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[1], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[1], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[4], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[3], indices[4], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[0], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[0], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[1], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[1], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[3], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[2], indices[4], indices[3], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[1], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[1], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[2], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[2], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[4], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[0], indices[4], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[0], indices[2], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[0], indices[4], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[2], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[2], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[4], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[1], indices[4], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[0], indices[1], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[0], indices[4], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[1], indices[0], indices[4], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[1], indices[4], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[4], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[2], indices[4], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[0], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[0], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[1], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[1], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[2], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[3], indices[4], indices[2], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[1], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[1], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[2], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[2], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[3], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[0], indices[3], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[0], indices[2], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[0], indices[3], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[2], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[2], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[3], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[1], indices[3], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[0], indices[1], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[0], indices[3], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[1], indices[0], indices[3], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[1], indices[3], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[3], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[2], indices[3], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[0], indices[1], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[0], indices[2], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[1], indices[0], indices[2], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[1], indices[2], indices[0], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[2], indices[0], indices[1], JValue, maxClosing, ClebshArray);
	Tmp += this->ComputeProjectorCoefficients5BodyWithDirection(indices[4], indices[3], indices[2], indices[1], indices[0], JValue, maxClosing, ClebshArray);
	TmpCoefficients[i] = Tmp;
	indices += 5;
      }
  delete[] ClebshArray;
  return TmpCoefficients;
}

// compute a given projector coefficient for the 5-body interaction 
//
// m1 = first index
// m2 = second index
// m3 = third inde
// m4 = fourth index
// m5 = fifth index
// jValue = total angular momentum
// minJ = minimum angular momentum that can be reach by three particles
// return value = corresponding projector coefficient

double ParticleOnSphereGenericFiveBodyHamiltonian::ComputeProjectorCoefficients5Body(int m1, int m2, int m3, int m4, int m5, int jValue, 
										     ClebschGordanCoefficients* clebshArray)
{
  double Tmp = 0.0;
  double Tmp2 = 0.0;
  double Tmp3 = 0.0;
  int TmpMinJ2;
  int TmpMinJ3;
  int Sum1 = ((m1 + m2) << 1)  - (2 * this->LzMax);
  int Sum2 = Sum1 + (m3 << 1)  - this->LzMax;
  int Sum3 = Sum2 + (m4 << 1)  - this->LzMax;
  if (abs(Sum3 + ((m5 << 1) - this->LzMax)) > jValue)
    return 0.0;
  int TmpMinJ = abs(Sum1);
  int MaxJ = 2 * this->LzMax;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    MaxJ -= 2;
  for (int j1 = MaxJ; j1 >= TmpMinJ; j1 -= 4)
    {
      Tmp2 = clebshArray[this->LzMax].GetCoefficient(((m1 << 1) - this->LzMax), ((m2 << 1)- this->LzMax), j1);
      TmpMinJ2 = abs(j1 - this->LzMax);
      if (TmpMinJ2 < abs(Sum2))
	TmpMinJ2 = abs(Sum2);
      for (int j2 = j1 + this->LzMax; j2 >= TmpMinJ2; j2 -= 2)
	{
	  if (abs(Sum1 + ((m3 << 1) - this->LzMax)) <= j2)
	    {
	      Tmp3 = Tmp2 * clebshArray[j1].GetCoefficient(Sum1, ((m3 << 1) - this->LzMax), j2);
	      TmpMinJ3 = abs(jValue - this->LzMax);
	      if (TmpMinJ3 < abs(Sum3))
		TmpMinJ3 = abs(Sum3);	      
	      for (int j3 = j2 + this->LzMax; j3 >= TmpMinJ3; j3 -= 2)
		{
		  if ((abs(Sum2 + ((m4 << 1) - this->LzMax)) <= j3) && (abs(j3 - this->LzMax) <= jValue) && (abs(j2 - this->LzMax) <= j3))
		    Tmp += Tmp3 * (clebshArray[j2].GetCoefficient(Sum2, ((m4 << 1) - this->LzMax), j3) * 
				   clebshArray[j3].GetCoefficient(Sum3, ((m5 << 1) - this->LzMax), jValue)); 
		}
	    }
	}
    }
  return Tmp;
}

// compute a given projector coefficient for the 5-body interaction in a given direction
//
// m1 = first index
// m2 = second index
// m3 = third inde
// m4 = fourth index
// m5 = fifth index
// jValue = total angular momentum
// maxClosing = array that gives the maximum angular momentum when a particle approach the cluster of n+1 particles  (n being the index of MaxClosing)
// return value = corresponding projector coefficient

double ParticleOnSphereGenericFiveBodyHamiltonian::ComputeProjectorCoefficients5BodyWithDirection(int m1, int m2, int m3, int m4, int m5, int jValue, 
												  int* maxClosing, ClebschGordanCoefficients* clebshArray)
{
  double Tmp = 0.0;
  double Tmp2 = 0.0;
  double Tmp3 = 0.0;
  int TmpMinJ2;
  int TmpMinJ3;
  int Sum1 = ((m1 + m2) << 1)  - (2 * this->LzMax);
  int Sum2 = Sum1 + (m3 << 1)  - this->LzMax;
  int Sum3 = Sum2 + (m4 << 1)  - this->LzMax;
  if (abs(Sum3 + ((m5 << 1) - this->LzMax)) > jValue)
    return 0.0;
  int TmpMinJ = abs(Sum1);
  int MaxJ = 2 * this->LzMax;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    MaxJ -= 2;
  MaxJ -= maxClosing[0];
  for (int j1 = MaxJ; j1 >= TmpMinJ; j1 -= 4)
    {
      Tmp2 = clebshArray[this->LzMax].GetCoefficient(((m1 << 1) - this->LzMax), ((m2 << 1)- this->LzMax), j1);
      TmpMinJ2 = abs(j1 - this->LzMax);
      if (TmpMinJ2 < abs(Sum2))
	TmpMinJ2 = abs(Sum2);
      for (int j2 = j1 + this->LzMax - maxClosing[1]; j2 >= TmpMinJ2; j2 -= 2)
	{
	  if (abs(Sum1 + ((m3 << 1) - this->LzMax)) <= j2)
	    {
	      Tmp3 = Tmp2 * clebshArray[j1].GetCoefficient(Sum1, ((m3 << 1) - this->LzMax), j2);
	      TmpMinJ3 = abs(jValue - this->LzMax);
	      if (TmpMinJ3 < abs(Sum3))
		TmpMinJ3 = abs(Sum3);	      
	      for (int j3 = j2 + this->LzMax - maxClosing[2]; j3 >= TmpMinJ3; j3 -= 2)
		{
		  if ((abs(Sum2 + ((m4 << 1) - this->LzMax)) <= j3) && (abs(j3 - this->LzMax) <= jValue) && (abs(j2 - this->LzMax) <= j3))
		    Tmp += Tmp3 * (clebshArray[j2].GetCoefficient(Sum2, ((m4 << 1) - this->LzMax), j3) * 
				   clebshArray[j3].GetCoefficient(Sum3, ((m5 << 1) - this->LzMax), jValue)); 
		}
	    }
	}
    }
  return Tmp;
}

