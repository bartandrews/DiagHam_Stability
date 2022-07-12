////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                           Class author : Cecile Repellin                   //
//                                                                            //
//    class of hamiltonian associated to particles on a 4D sphere with        //
//                           delta 3-body interaction                         //
//                                                                            //
//                        last modification : 10/01/2013                      //
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
#include "Hamiltonian/ParticleOnCP2ThreeBodyDeltaHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"
#include "MathTools/FactorialCoefficient.h"
  
#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//

ParticleOnCP2ThreeBodyDeltaHamiltonian::ParticleOnCP2ThreeBodyDeltaHamiltonian()
{
}

// constructor from datas with three-body interaction and optional two-body and one-body
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
//  nbrFluxQuanta = maximum J value reached by a particle in the state
// onebodyPotential = array with the one-body potentials
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnCP2ThreeBodyDeltaHamiltonian::ParticleOnCP2ThreeBodyDeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta, double* onebodyPotentials,
											 AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
											 char* precalculationFileName)
{
  this->Particles = particles;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->NbrLzValue = (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 2) /2;
  this->LzMax = NbrLzValue - 1;
  this->NbrParticles = nbrParticles;


  this->FullTwoBodyFlag = false;
  this->OneBodyTermFlag = false;
  if (onebodyPotentials != 0)
    {
      this->OneBodyPotentials = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->OneBodyPotentials[i] = onebodyPotentials[this->LzMax - i];
      cout << "Setting up one-body potential" << endl;
      this->OneBodyTermFlag = true;
    }
  else
    {
      this->OneBodyTermFlag = false;
    }

  this->MaxNBody = 3;
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

//   this->MaxRelativeAngularMomentum = maxRelativeAngularMomentum;
//   this->NbrThreeBodyPseudoPotential = maxRelativeAngularMomentum;
//   this->ThreeBodyPseudoPotential = new double[this->NbrThreeBodyPseudoPotential + 1];
//   for (int i = 0; i <= this->NbrThreeBodyPseudoPotential; ++i)
//     this->ThreeBodyPseudoPotential[i] = threeBodyPseudoPotential[i];
//   this->NormalizeFlag=normalizePPs;

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 1;
      this->MaxSumIndices[k] = 0;      
      this->NBodyFlags[k] = false;
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	this->NBodySign[k] = -1.0;
    }
  this->NBodyFlags[3] = true;
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
  this->L2Operator = 0;
}

// destructor
//

ParticleOnCP2ThreeBodyDeltaHamiltonian::~ParticleOnCP2ThreeBodyDeltaHamiltonian()
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

  delete[] this->NbrNIndices;
  delete[] this->NIndices;
  delete[] this->NbrMIndices;
  delete[] this->MIndices;
  delete[] this->MNNBodyInteractionFactors;
//   delete[] this->ThreeBodyPseudoPotential;

  delete[] this->NBodyFlags;
  delete[] this->NBodyInteractionFactors;
  delete[] this->SortedIndicesPerSum;
  delete[] this->NbrSortedIndicesPerSum;
  delete[] this->MinSumIndices;
  delete[] this->MaxSumIndices;
  delete[] this->NBodySign;
//   if (this->L2Operator != 0)
//     delete this->L2Operator;
  if (this->FullTwoBodyFlag == true)
    {
      delete[] this->InteractionFactors;
      delete[] this->M1Value;
      delete[] this->M2Value;
      delete[] this->M3Value;
//       delete[] this->PseudoPotential;
    }
  if (this->OneBodyTermFlag == true)
    {
      delete[] this->OneBodyInteractionFactors;
      delete[] this->OneBodyPotentials;
      delete[] this->OneBodyMValues;
      delete[] this->OneBodyNValues;
    }
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* ParticleOnCP2ThreeBodyDeltaHamiltonian::Clone ()
{
  return 0;
}


// evaluate all interaction factors
//   

void ParticleOnCP2ThreeBodyDeltaHamiltonian::EvaluateInteractionFactors()
{
      // three-body interactions for bosons

   this->MinSumIndices[3] = 0;
   this->MaxSumIndices[3] = (3*this->NbrFluxQuanta +1)*(3*this->NbrFluxQuanta + 2)/2 - 1;

   double** InteractionFactor;
   GetAllIndices(this->NbrLzValue, this->NbrFluxQuanta, this->MaxSumIndices[3] + 1, this->NbrSortedIndicesPerSum[3], this->SortedIndicesPerSum[3], InteractionFactor); 
//         for (int MinSum = 0; MinSum <= this->MaxSumIndices[3]; ++MinSum)
// 	{
// 	      cout << "Sum = " << MinSum << endl;
// 	      int Lim = this->NbrSortedIndicesPerSum[3][MinSum];
// 	      cout << "nbr multiplets = " << Lim << endl;
// 	 	  for (int i = 0; i < Lim; ++i)
// 		      {
// 			cout << this->SortedIndicesPerSum[3][MinSum][3 * i] << "," <<  this->SortedIndicesPerSum[3][MinSum][3 * i + 1] << "," <<  this->SortedIndicesPerSum[3][MinSum][3 * i + 2] << endl;  
// 		      }
// 	}
	
    long TmpNbrNIndices = 0;
    for (int MinSum = 0; MinSum <= this->MaxSumIndices[3]; ++MinSum)
      TmpNbrNIndices += this->NbrSortedIndicesPerSum[3][MinSum];
    this->NbrNIndices[3] = TmpNbrNIndices;
    this->NIndices[3] = new int[TmpNbrNIndices * 3];
    this->NbrMIndices[3] = new long[TmpNbrNIndices];
    this->MIndices[3] = new int*[TmpNbrNIndices];
    this->MNNBodyInteractionFactors[3] = new double* [TmpNbrNIndices];
    TmpNbrNIndices = 0;	 
    int* TmpNIndices = this->NIndices[3];
    for (int MinSum = 0; MinSum <= this->MaxSumIndices[3]; ++MinSum)
      {
	int Lim = this->NbrSortedIndicesPerSum[3][MinSum];
 	int* TmpNIndices2 = this->SortedIndicesPerSum[3][MinSum];
	for (int i = 0; i < Lim; ++i)
	  {
	    this->NbrMIndices[3][TmpNbrNIndices] = Lim;		    
	    this->MIndices[3][TmpNbrNIndices] = new int [Lim * 3];
	    this->MNNBodyInteractionFactors[3][TmpNbrNIndices] = new double [Lim];
	    int* TmpMIndices = this->MIndices[3][TmpNbrNIndices];
	    int* TmpMIndices2 = this->SortedIndicesPerSum[3][MinSum];

	    for (int j = 0; j < Lim; ++j)
	      {
		for (int l = 0; l < 3; ++l)
		  {
		    (*TmpMIndices) = (*TmpMIndices2);			
		    ++TmpMIndices;
		    ++TmpMIndices2;
		  }			
		  this->MNNBodyInteractionFactors[3][TmpNbrNIndices][j] = InteractionFactor[MinSum][i]*InteractionFactor[MinSum][j];
	      }
		for (int j = 0; j < 3; ++j)
		  {
		    (*TmpNIndices) = (*TmpNIndices2);			
		    ++TmpNIndices;
		    ++TmpNIndices2;
		  }
	      ++TmpNbrNIndices;
	    }		
	}
     delete[] InteractionFactor;

  if (this->OneBodyTermFlag == true)
    {
      this->NbrOneBodyInteractionFactors = 0;
      for (int i = 0; i <= this->LzMax; ++i)
	if (this->OneBodyPotentials[i] != 0)
	  ++this->NbrOneBodyInteractionFactors;
      if (this->NbrOneBodyInteractionFactors != 0)
	{
	  double Sign = 1.0;
	  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	    Sign = -1.0;
	  this->OneBodyMValues = new int[this->NbrOneBodyInteractionFactors];
	  this->OneBodyNValues = new int[this->NbrOneBodyInteractionFactors];
	  this->OneBodyInteractionFactors = new double[this->NbrOneBodyInteractionFactors];
	  this->NbrOneBodyInteractionFactors = 0;
	  for (int i = 0; i <= this->LzMax; ++i)
	    if (this->OneBodyPotentials[i] != 0)
	      {
		this->OneBodyMValues[this->NbrOneBodyInteractionFactors] = i;
		this->OneBodyNValues[this->NbrOneBodyInteractionFactors] = i;
		cout << this->OneBodyPotentials[i] << endl;
		this->OneBodyInteractionFactors[this->NbrOneBodyInteractionFactors] = this->OneBodyPotentials[i] * Sign;
		++this->NbrOneBodyInteractionFactors;
	      }	  
	}
      else
	{
	  delete[] this->OneBodyPotentials;
	  this->OneBodyTermFlag = false;
	}
    }
}



// get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices and compute the corresponding interaction factor
//
// nbrValues = number of different values an index can have
// nbrFluxQuanta = number of flux quanta
// nbrSectorSum = number of triplets ((r1,s1), (r2,s2), (r3,s3))  with a different value of (r1 + r2 + r3, s1 + s2 + s3)
// nbrSortedIndicesPerSum = reference on an array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on an array where group of indices are stored (first array dimension corresponding to sum of the indices)
// interactionCoef = reference to an array of factorial coefficients  (interaction factors)


void ParticleOnCP2ThreeBodyDeltaHamiltonian::GetAllIndices (int nbrValues, int nbrFluxQuanta, int nbrSectorSum, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum, double**& interactionCoef)
{
  BosonOnCP2* TmpHilbertSpace = (BosonOnCP2*) this->Particles;
  int* QuantumNumberTz = new int [this->NbrLzValue];
  int* QuantumNumberY = new int [this->NbrLzValue];
  int* QuantumNumberR = new int [this->NbrLzValue];
  int* QuantumNumberS = new int [this->NbrLzValue];
  TmpHilbertSpace->GetQuantumNumbersFromLinearizedIndex(QuantumNumberTz, QuantumNumberY, QuantumNumberR, QuantumNumberS);
  nbrSortedIndicesPerSum = new int[nbrSectorSum];
  sortedIndicesPerSum = new int*[nbrSectorSum];
  interactionCoef = new double*[nbrSectorSum];
  FactorialCoefficient Coef;
  for (int i = 0; i < nbrSectorSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  
  for (int l = 0; l < nbrValues; ++l)
  {
    int r1 = QuantumNumberR[l];
    int s1 = QuantumNumberS[l];
    for (int m = 0; m < nbrValues; ++m)
    {
      int r2 = QuantumNumberR[m];
      int s2 = QuantumNumberS[m];
      for (int n = 0; n < nbrValues; ++n)
      {
	int r3 = QuantumNumberR[n];
	int s3= QuantumNumberS[n];
	int tz = r1 + r2 + r3 - (s1 + s2 + s3);
	int y = 3 * (r1 + r2 + r3 + s1 + s2 + s3) - 6 * this->NbrFluxQuanta;
	int IndexSum =  TmpHilbertSpace->GetLinearizedIndex(tz, y, 3); 
	++nbrSortedIndicesPerSum[IndexSum];
      }
    }
  }
  
  for (int i = 0; i < nbrSectorSum; ++i)
  {
    sortedIndicesPerSum[i] = new int[3*nbrSortedIndicesPerSum[i]];
    interactionCoef[i] = new double[nbrSortedIndicesPerSum[i]];
    nbrSortedIndicesPerSum[i] = 0;
   }

  for (int l = 0; l < nbrValues; ++l)
  {
    int r1 = QuantumNumberR[l];
    int s1 = QuantumNumberS[l];
    int t1 = this->NbrFluxQuanta - r1 -s1;
    for (int m = 0; m < nbrValues; ++m)
    {
      int r2 = QuantumNumberR[m];
      int s2 = QuantumNumberS[m];
      int t2 = this->NbrFluxQuanta - r2 -s2;
      for (int n = 0; n < nbrValues; ++n)
      {
	int r3 = QuantumNumberR[n];
	int s3= QuantumNumberS[n];
	int t3 = this->NbrFluxQuanta - r3 -s3;
	int tz = r1 + r2 + r3 - (s1 + s2 + s3);
	int y = 3 * (r1 + r2 + r3 + s1 + s2 + s3) - 6 * this->NbrFluxQuanta;
	int IndexSum = TmpHilbertSpace->GetLinearizedIndex(tz, y, 3); 
	sortedIndicesPerSum[IndexSum][3*nbrSortedIndicesPerSum[IndexSum]] = l;
	sortedIndicesPerSum[IndexSum][3*nbrSortedIndicesPerSum[IndexSum] + 1] = m;
	sortedIndicesPerSum[IndexSum][3*nbrSortedIndicesPerSum[IndexSum] + 2] = n;
	
	
	Coef.SetToOne();
	Coef.FactorialMultiply(r1 + r2 + r3);
	Coef.FactorialDivide(r1);
	Coef.FactorialDivide(r2);
	Coef.FactorialDivide(r3);
	Coef.FactorialMultiply(s1 + s2 + s3);
	Coef.FactorialDivide(s1);
	Coef.FactorialDivide(s2);
	Coef.FactorialDivide(s3);
	Coef.FactorialMultiply(t1 + t2 + t3);
	Coef.FactorialDivide(t1);
	Coef.FactorialDivide(t2);
	Coef.FactorialDivide(t3);
	Coef.FactorialDivide(3*this->NbrFluxQuanta + 2);
	Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
	Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
	Coef.FactorialMultiply(this->NbrFluxQuanta + 2);
	interactionCoef[IndexSum][nbrSortedIndicesPerSum[IndexSum]] = sqrt(Coef.GetNumericalValue());
	
	++nbrSortedIndicesPerSum[IndexSum];
      }
    }
  }
}
