////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                          generic 3-body interaction                        //
//                                                                            //
//                        last modification : 22/07/2008                      //
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
#include "Hamiltonian/ParticleOnSphereGenericThreeBodyHamiltonian.h"
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

ParticleOnSphereGenericThreeBodyHamiltonian::ParticleOnSphereGenericThreeBodyHamiltonian()
{
}

// constructor from datas with three-body interaction and optional two-body and one-body
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// threeBodyPseudoPotential = array with the three-body pseudo-potentials sorted with respect to the relative angular momentum, 
//                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
// maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction)
// onebodyPotential = array with the one-body potentials
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
// normalizePPs = use normalized n-body eigenstates as the reference potentials
ParticleOnSphereGenericThreeBodyHamiltonian::ParticleOnSphereGenericThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
											 double* threeBodyPseudoPotential, int maxRelativeAngularMomentum,
											 double l2Factor, double* pseudoPotential, double* onebodyPotentials,
											 AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
											 char* precalculationFileName, bool normalizePPs)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FullTwoBodyFlag = false;
  if (pseudoPotential != 0)
    {
      this->PseudoPotential = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->PseudoPotential[i] = pseudoPotential[this->LzMax - i];
      cout << "Setting up two-body potential" << endl;
      this->FullTwoBodyFlag = true;
    }
  else
    {
      this->FullTwoBodyFlag = false;
    }

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

  this->MaxRelativeAngularMomentum = maxRelativeAngularMomentum;
  this->NbrThreeBodyPseudoPotential = maxRelativeAngularMomentum;
  this->ThreeBodyPseudoPotential = new double[this->NbrThreeBodyPseudoPotential + 1];
  for (int i = 0; i <= this->NbrThreeBodyPseudoPotential; ++i)
    this->ThreeBodyPseudoPotential[i] = threeBodyPseudoPotential[i];
  this->NormalizeFlag=normalizePPs;

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

ParticleOnSphereGenericThreeBodyHamiltonian::~ParticleOnSphereGenericThreeBodyHamiltonian()
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
  delete[] this->ThreeBodyPseudoPotential;

  delete[] this->NBodyFlags;
  delete[] this->NBodyInteractionFactors;
  delete[] this->SortedIndicesPerSum;
  delete[] this->NbrSortedIndicesPerSum;
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

AbstractHamiltonian* ParticleOnSphereGenericThreeBodyHamiltonian::Clone ()
{
  return 0;
}


// evaluate all interaction factors
//   

void ParticleOnSphereGenericThreeBodyHamiltonian::EvaluateInteractionFactors()
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
      GetAllSkewSymmetricIndices(this->NbrLzValue, 3, this->NbrSortedIndicesPerSum[3], this->SortedIndicesPerSum[3]);
      this->MaxSumIndices[3] = (this->LzMax * 3) - 2;
      this->MinSumIndices[3] = 3;
      double* TmpInteractionCoeffients = new double[this->MaxSumIndices[3] + 1];
      TmpInteractionCoeffients[0] = 1.0;
      TmpInteractionCoeffients[1] = 1.0;
      for (int i = 2; i <= this->MaxSumIndices[3]; ++i)
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
      Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[3]) + 1.0);
      double Radius = 2.0 / ((double) this->LzMax);
      for (int i = 2; i <= 3; ++i)
	{
	  Coefficient *= (double) (i * i);	  
	  Coefficient *= Radius;
	}
      for (int i = 0; i <= this->MaxSumIndices[3]; ++i)
	TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);
      
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
	  if (Lim > 0)
	    {
// 	      cout << "--------------------------------------------" << endl;
//  	      cout << "index sum = " << MinSum << endl;
// 	      cout << "--------------------------------------------" << endl;
 	      int* TmpNIndices2 = this->SortedIndicesPerSum[3][MinSum];
	      int TmpMaxRelativeMomentum = 20;
	      if (this->MaxRelativeAngularMomentum <= TmpMaxRelativeMomentum)
		TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum;
	      int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2];
	      while (abs((TmpSum * 2) - (3 * this->LzMax)) > ((3 * this->LzMax) - (2 * TmpMaxRelativeMomentum)))
		--TmpMaxRelativeMomentum;
	      double** TmpProjectorCoefficients = new double* [TmpMaxRelativeMomentum + 1];
//	      cout << "V^3_3 : " << endl;
	      if (this->ThreeBodyPseudoPotential[3] != 0.0)
		TmpProjectorCoefficients[3] = this->ComputeProjectorCoefficients(6, 1, TmpNIndices2, Lim);
//	      cout << "--------------------------------------------" << endl;
	      for (int i = 5; i <= TmpMaxRelativeMomentum; ++i)  
		if (this->ThreeBodyPseudoPotential[i] != 0.0)
		  {
//		    cout << "V^3_" << i << " : " << endl;
		    TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim);
//		    cout << "--------------------------------------------" << endl;
 		  }
	      for (int i = 0; i < Lim; ++i)
		{
		  this->NbrMIndices[3][TmpNbrNIndices] = Lim;		    
		  this->MIndices[3][TmpNbrNIndices] = new int [Lim * 3];
		  this->MNNBodyInteractionFactors[3][TmpNbrNIndices] = new double [Lim];
		  int* TmpMIndices = this->MIndices[3][TmpNbrNIndices];
		  int* TmpMIndices2 = this->SortedIndicesPerSum[3][MinSum];
		  double* TmpInteraction = this->MNNBodyInteractionFactors[3][TmpNbrNIndices];
		  for (int j = 0; j < Lim; ++j)
		    {
		      for (int l = 0; l < 3; ++l)
			{
			  (*TmpMIndices) = (*TmpMIndices2);			
			  ++TmpMIndices;
			  ++TmpMIndices2;
			}			
		      double& TmpInteraction2 = TmpInteraction[j];
		      TmpInteraction2 = 0.0;
		      if (this->ThreeBodyPseudoPotential[3] != 0.0)
			TmpInteraction2 += this->ThreeBodyPseudoPotential[3] * TmpProjectorCoefficients[3][i] * TmpProjectorCoefficients[3][j];
		      for (int k = 5; k <= TmpMaxRelativeMomentum; ++k)  
			if (this->ThreeBodyPseudoPotential[k] != 0.0)
			  TmpInteraction2 += this->ThreeBodyPseudoPotential[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
		    }
		  for (int j = 0; j < 3; ++j)
		    {
		      (*TmpNIndices) = (*TmpNIndices2);			
		      ++TmpNIndices;
		      ++TmpNIndices2;
		    }
		  ++TmpNbrNIndices;
		}
	      if (this->ThreeBodyPseudoPotential[3] != 0.0)
		delete[] TmpProjectorCoefficients[3];
	      for (int i = 5; i <= TmpMaxRelativeMomentum; ++i)  
		if (this->ThreeBodyPseudoPotential[i] != 0.0)
		  delete[] TmpProjectorCoefficients[i];
	      delete[] TmpProjectorCoefficients;		
	    }
	}
      delete[] TmpInteractionCoeffients;
    }
  else // three-body interactions for bosons
    {
      this->MinSumIndices[3] = 0;
      this->MaxSumIndices[3] = this->LzMax * 3;
      double* TmpInteractionCoeffients = new double[this->MaxSumIndices[3] + 1];
      double Coefficient;
      TmpInteractionCoeffients[0] = 1.0;
      TmpInteractionCoeffients[1] = 1.0;
      for (int i = 2; i <= this->MaxSumIndices[3]; ++i)
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
      Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[3]) + 1.0);
      double Radius = 2.0 / ((double) this->LzMax);
      for (int i = 2; i <= 3; ++i)
	{
	  Coefficient *= (double) (i * i);	  
	  Coefficient *= Radius;
	}
      for (int i = 0; i <= this->MaxSumIndices[3]; ++i)
	TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);
      
      double** SortedIndicesPerSumSymmetryFactor;
      GetAllSymmetricIndices(this->NbrLzValue, 3, this->NbrSortedIndicesPerSum[3], this->SortedIndicesPerSum[3],
			     SortedIndicesPerSumSymmetryFactor);
      
      
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
	  double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
	  int* TmpNIndices2 = this->SortedIndicesPerSum[3][MinSum];
	  int TmpMaxRelativeMomentum = 10;
	  if (this->MaxRelativeAngularMomentum <= TmpMaxRelativeMomentum)
	    TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum;
	  int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2];
	  while (((3 * this->LzMax) - TmpMaxRelativeMomentum)  < TmpSum)
	    --TmpMaxRelativeMomentum;
	  double** TmpProjectorCoefficients = new double* [TmpMaxRelativeMomentum + 1];
	  if (this->ThreeBodyPseudoPotential[0] != 0.0)
	    TmpProjectorCoefficients[0] = this->ComputeProjectorCoefficients(0, 1, TmpNIndices2, Lim);
	  for (int i = 2; i <= TmpMaxRelativeMomentum; ++i)  
	    if (this->ThreeBodyPseudoPotential[i] != 0.0)
	      TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim);
	  for (int i = 0; i < Lim; ++i)
	    {
	      this->NbrMIndices[3][TmpNbrNIndices] = Lim;		    
	      this->MIndices[3][TmpNbrNIndices] = new int [Lim * 3];
	      this->MNNBodyInteractionFactors[3][TmpNbrNIndices] = new double [Lim];
	      int* TmpMIndices = this->MIndices[3][TmpNbrNIndices];
	      int* TmpMIndices2 = this->SortedIndicesPerSum[3][MinSum];
	      double* TmpInteraction = this->MNNBodyInteractionFactors[3][TmpNbrNIndices];
	      for (int j = 0; j < Lim; ++j)
		{
		  double Coefficient2 = TmpSymmetryFactors[j];
		  for (int l = 0; l < 3; ++l)
		    {
		      Coefficient2 *= TmpNormalizationCoeffients[(*TmpMIndices2)];		    
		      (*TmpMIndices) = (*TmpMIndices2);			
		      ++TmpMIndices;
		      ++TmpMIndices2;
		    }			
		  double& TmpInteraction2 = TmpInteraction[j];
		  TmpInteraction2 = 0.0;
		  if (this->ThreeBodyPseudoPotential[0] != 0.0)
		    TmpInteraction2 += this->ThreeBodyPseudoPotential[0] * TmpProjectorCoefficients[0][i] * TmpProjectorCoefficients[0][j];
		  for (int k = 2; k <= TmpMaxRelativeMomentum; ++k)  
		    if (this->ThreeBodyPseudoPotential[k] != 0.0)
		      TmpInteraction2 += this->ThreeBodyPseudoPotential[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
		  TmpInteraction2 *= TmpSymmetryFactors[i] * TmpSymmetryFactors[j];
		}
	      for (int j = 0; j < 3; ++j)
		{
		  (*TmpNIndices) = (*TmpNIndices2);			
		  ++TmpNIndices;
		  ++TmpNIndices2;
		}
	      ++TmpNbrNIndices;
	    }		
	  if (this->ThreeBodyPseudoPotential[0] != 0.0)
	    delete[] TmpProjectorCoefficients[0];
	  for (int i = 2; i <= TmpMaxRelativeMomentum; ++i)  
	    if (this->ThreeBodyPseudoPotential[i] != 0.0)
	      delete[] TmpProjectorCoefficients[i];
	  delete[] TmpProjectorCoefficients;		
	}
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[3]; ++MinSum)
	{
	  delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	}
      delete[] SortedIndicesPerSumSymmetryFactor;
      delete[] TmpInteractionCoeffients;
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
	  double Factor = - 4.0;
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
	  double Factor = 4.0;
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

// compute all projector coefficient associated to a given relative angular momentum between 3 particles
//
// relativeMomentum = value of twice the relative angular momentum between the 3 particles
// degeneracyIndex = optional degeneracy index for relative angular momentum greater than 5 for bosons (8 for fermions)
// indices = array that contains all possible sets of indices (size of the array is 3 * nbrIndexSets)
// nbrIndexSets = number of sets

double* ParticleOnSphereGenericThreeBodyHamiltonian::ComputeProjectorCoefficients(int relativeMomentum, int degeneracyIndex, int* indices, int nbrIndexSets)
{
  double* TmpCoefficients = new double [nbrIndexSets];
  int JValue = (3 * this->LzMax) - relativeMomentum;

  if (abs(((indices[0] + indices[1] + indices[2]) << 1)  - (3 * this->LzMax)) > JValue)
    {
      for (int i = 0; i < nbrIndexSets; ++i)
	TmpCoefficients[i] = 0.0;
      return TmpCoefficients;
    }

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
      double Tmp2 = 0.0;
      for (int j = MaxJ; j >= TmpMinJ; j -= 4)
	{
	  Tmp += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) * 
		  ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue)); 
	  Tmp2 += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) *
		   ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue));
	}
//      cout << indices[0] << " " << indices[1] << " " << indices[2] << " : " << Tmp2 << endl;
      Sum = ((indices[1] + indices[2]) << 1)  - (2 * this->LzMax);
      TmpMinJ = MinJ;
      if (TmpMinJ < abs(Sum))
	TmpMinJ = abs(Sum);
      Tmp2 = 0.0;
      for (int j = MaxJ; j >= TmpMinJ; j -= 4)
	{
	  Tmp += (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) * 
		  ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue)); 
	  Tmp2 += (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) *
                  ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue));
	}
//      cout << indices[1] << " " << indices[2] << " " << indices[0] << " : " << Tmp2 << endl;
      Sum = ((indices[2] + indices[0]) << 1)  - (2 * this->LzMax);
      TmpMinJ = MinJ;
      if (TmpMinJ < abs(Sum))
	TmpMinJ = abs(Sum);
      Tmp2 = 0.0;
      for (int j = MaxJ; j >= TmpMinJ; j -= 4)
	{
	  Tmp += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) * 
		  ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue)); 
	  Tmp2 += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) *
                  ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue));
	}
//      cout << indices[2] << " " << indices[0] << " " << indices[1] << " : " << Tmp2 << endl;
      TmpCoefficients[i] = Tmp;
      indices += 3;
    }
  delete[] ClebshArray;

  if (NormalizeFlag)
    {
      double sum=0.0;
      for (int i=0; i<nbrIndexSets; ++i)
	sum+=TmpCoefficients[i]*TmpCoefficients[i];
      if (sum>1e-10)
	for (int i=0; i<nbrIndexSets; ++i)
	  TmpCoefficients[i]/=sqrt(sum);
      else
	cout << "Attention: cannot normalize 3-body pseudopotential with relative angular momentum "<<relativeMomentum<<endl;
    }
  return TmpCoefficients;
}

