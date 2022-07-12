////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//               forcing vanishing properties with multiple groups            //
//                        of n-body hard core interaction                     //
//                                                                            //
//                        last modification : 04/06/2011                      //
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
#include "Hamiltonian/ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian.h"
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

ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian::ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// nbrGroups = number of groups that have a given vanishing property
// nbrNBodys = number of particle that interact simultaneously through the hard core interaction for the each group
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian::ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int nbrGroups, int* nbrNBodys, double l2Factor, AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag,  char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = false;

  this->NbrGroups = nbrGroups;
  this->NbrNBodys = new int [this->NbrGroups];
  
  this->NbrNbody = 0;
  for (int i = 0; i < this->NbrGroups; ++i)
    {
      this->NbrNBodys[i] = nbrNBodys[i];
      this->NbrNbody += this->NbrNBodys[i];
    }
    
  this->MaxNBody = this->NbrNbody;
  
//   this->GroupsNBodyFlags = new bool [this->MaxNBody + 1];
//   for (int k = 0; k <= this->MaxNBody; ++k)
//     {
//       this->GroupsNBodyFlags[k] = false;
//     }
//     
//     for (int i = 0; i < this->NbrGroups; ++i)
//     {
//       this->GroupsNBodyFlags[this->NbrNBodys[i]]= true;
//     }
    
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


// destructor
//

ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian::~ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian()
{
  delete[] this->NbrNBodys;
}

// evaluate all interaction factors
//   

void ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian::EvaluateInteractionFactors()
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
      cout << "fermionic case is not implemented" << endl;
    }
  else
    {
      double Coefficient;
      double InverseRadiusSquare = 2.0 / ((double) this->LzMax);
      double ** TmpInteractionCoeffients = new double * [this->NbrGroups];
      for (int i = 0; i < this->NbrGroups; i++)
      {
	cout <<"i = "<< i <<" Nbr Body: "<<this->NbrNBodys[i] <<endl;
	TmpInteractionCoeffients[i] = new double[this->LzMax * this->NbrNBodys[i] + 1];
	
	double TmpCoefficient;
	TmpInteractionCoeffients[i][0] = 1.0;
	TmpInteractionCoeffients[i][1] = 1.0;
	for (int k = 2; k <= this->LzMax * this->NbrNBodys[i]; ++k)
	    {
	      TmpCoefficient = 1.0;
	      for (int j = 1; j < k; ++j)
		 {
		   double Coefficient2 = TmpInteractionCoeffients[i][j];
		    TmpInteractionCoeffients[i][j] += TmpCoefficient;
		    TmpCoefficient = Coefficient2;
		  }
		TmpInteractionCoeffients[i][k] = 1.0;
	      }
	     
// 	     cout<< "TmpInteractionCoeffients[i] : ";
// 	     for (int k = 0; k <= this->LzMax * this->NbrNBodys[i]; ++k)
// 	    {
// 	     cout <<  TmpInteractionCoeffients[i][k] << " ";
// 	    }
// 	    cout <<endl;
	    
	    Coefficient = 4.0 * M_PI / (((double) this->NbrNBodys[i] * this->LzMax) + 1.0);

	    for (int k = 2; k <= this->NbrNBodys[i]; k++)
	      {	  
		
		Coefficient *= InverseRadiusSquare;
	      }
	      
	      
	    for (int p = 0; p <= this->LzMax * this->NbrNBodys[i]; p++)
	      TmpInteractionCoeffients[i][p] = (Coefficient / TmpInteractionCoeffients[i][p]);
      }
      
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    this->MinSumIndices[k] = 0;
	    this->MaxSumIndices[k] = this->LzMax * k;
	    
	      

	    double** SortedIndicesPerSumSymmetryFactor;
	    GetAllSymmetricIndices(this->NbrLzValue, k, this->NbrSortedIndicesPerSum[k], this->SortedIndicesPerSum[k],
				   SortedIndicesPerSumSymmetryFactor);
  
// 	    for (int l = 2; l <= this->NbrNbody; ++l)
// 	      {
// 		this->NBodyInteractionWeightFactors[k] *= (double) (l * l);
// 	      }
// 	    cout <<"this->NBodyInteractionWeightFactors[k] = "<< this->NBodyInteractionWeightFactors[k]<<endl;
	    
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
		    
		    
		    Coefficient = this->NBodyInteractionWeightFactors[k];//TmpSymmetryFactors[i] * this->NBodyInteractionWeightFactors[k];
		    for (int l = 0; l < k; ++l)
		      Coefficient *= TmpNormalizationCoeffients[TmpNIndices2[l]];
		    
		    int * TmpMIndicesTry = new int[k];
		    int * TmpNIndicesTry = new int[k];
		    
		    
		    for (int p = 0 ; p <k ; p++)
		    {
		      TmpNIndicesTry[p] = TmpNIndices2[p]; 
		      //cout << TmpNIndicesTry[p]<< " ";
		    }
		    //cout <<endl;
		    //cout <<"TmpSymmetryFactors[i] = "<<TmpSymmetryFactors[i]<<endl;
		      
		    //Coefficient *= TmpInteractionCoeffients[MinSum] * TmpInteractionCoeffients[MinSum];
		      
		    for (int j = 0; j < Lim; ++j)
		      {
			//cout <<"j ="<<j <<endl; 
			for (int p = 0 ; p <k ; p++)
		    {
		      TmpMIndicesTry[p] = TmpMIndices2[p]; 
		      //cout << TmpMIndicesTry[p]<< " ";
		    }
		    //cout <<endl;
		    //cout <<"TmpSymmetryFactors[j] = "<<TmpSymmetryFactors[j]<<endl;
		    
		    double UltimateCoef = 0.0;
		    
		    if (TmpMIndicesTry[0] != TmpMIndicesTry[1])
		    {
		      if (TmpMIndicesTry[1] != TmpMIndicesTry[2])
			{
			  if (TmpMIndicesTry[2] != TmpMIndicesTry[3])
			  {
			    
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[3], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[2], TmpMIndicesTry[1], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[2], TmpMIndicesTry[3], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[3], TmpMIndicesTry[1], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[3], TmpMIndicesTry[2], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[0], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[0], TmpMIndicesTry[3], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[2], TmpMIndicesTry[0], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[2], TmpMIndicesTry[3], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[3], TmpMIndicesTry[0], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[3], TmpMIndicesTry[2], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[1], TmpMIndicesTry[0], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[1], TmpMIndicesTry[3], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[0], TmpMIndicesTry[1], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[0], TmpMIndicesTry[3], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[3], TmpMIndicesTry[1], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[3], TmpMIndicesTry[0], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[1], TmpMIndicesTry[2], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[1], TmpMIndicesTry[0], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[2], TmpMIndicesTry[1], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[2], TmpMIndicesTry[0], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[0], TmpMIndicesTry[1], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[0], TmpMIndicesTry[2], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  }
			  else // TmpMIndicesTry[2] == TmpMIndicesTry[3]
			  {
			    UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[2], TmpMIndicesTry[1], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[2], TmpMIndicesTry[3], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[0], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[2], TmpMIndicesTry[0], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[2], TmpMIndicesTry[3], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[1], TmpMIndicesTry[0], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[1], TmpMIndicesTry[3], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[0], TmpMIndicesTry[1], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[0], TmpMIndicesTry[3], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[3], TmpMIndicesTry[1], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[3], TmpMIndicesTry[0], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  }
			}
			else
			{
			  if (TmpMIndicesTry[2] != TmpMIndicesTry[3]) //(TmpMIndices2[1] == TmpMIndices2[2]) 
			  {
			    
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[3], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[3], TmpMIndicesTry[1], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[0], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[0], TmpMIndicesTry[3], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[2], TmpMIndicesTry[0], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[2], TmpMIndicesTry[3], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[3], TmpMIndicesTry[0], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[3], TmpMIndicesTry[2], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[1], TmpMIndicesTry[2], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[1], TmpMIndicesTry[0], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[0], TmpMIndicesTry[1], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			    
			  }
			  else  //(TmpMIndices2[2] == TmpMIndices2[3]) && (TmpMIndices2[1] == TmpMIndices2[2])
			  {
			    UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[0], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[2], TmpMIndicesTry[0], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[1],TmpMIndicesTry[2], TmpMIndicesTry[3], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  }
			}
		    }
			  else //(TmpMIndices2[0] == TmpMIndices2[1])
			    {
			      if (TmpMIndicesTry[1] != TmpMIndicesTry[2])
			      {
				if (TmpMIndicesTry[2] != TmpMIndicesTry[3])
				{
				  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[3], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[2], TmpMIndicesTry[1], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[2], TmpMIndicesTry[3], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[3], TmpMIndicesTry[1], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[3], TmpMIndicesTry[2], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[1], TmpMIndicesTry[0], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[1], TmpMIndicesTry[3], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[3], TmpMIndicesTry[0], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[1], TmpMIndicesTry[2], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[1], TmpMIndicesTry[0], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[2], TmpMIndicesTry[0], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			      
				}
				else //(TmpMIndices2[0] == TmpMIndices2[1]) && (TmpMIndices2[2] == TmpMIndices2[3])
				{
				  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[2], TmpMIndicesTry[1], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[2], TmpMIndicesTry[3], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[1], TmpMIndicesTry[0], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[1], TmpMIndicesTry[3], TmpMIndicesTry[0],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[2],TmpMIndicesTry[3], TmpMIndicesTry[0], TmpMIndicesTry[1],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
				  
				}
			      }
			      else //(TmpMIndices2[0] == TmpMIndices2[1] == TmpMIndices2[2])
				{
				  if (TmpMIndicesTry[2] != TmpMIndicesTry[3])
				  {
				    UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[3], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[3], TmpMIndicesTry[1], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[3],TmpMIndicesTry[0], TmpMIndicesTry[1], TmpMIndicesTry[2],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
				    
				    
				  }
				else
				{
				   UltimateCoef += this->EvaluateSymmetrizedInteractionCoefficient(TmpMIndicesTry[0],TmpMIndicesTry[1], TmpMIndicesTry[2], TmpMIndicesTry[3],TmpNIndicesTry[0],TmpNIndicesTry[1], TmpNIndicesTry[2], TmpNIndicesTry[3], TmpInteractionCoeffients);
				}
				
				  
				}
			}
			  
			//cout <<"UltimateCoef = "<< UltimateCoef<<endl;
			double Coefficient2 = 1.0;//TmpSymmetryFactors[j];
			for (int l = 0; l < k; ++l)
			  {
			    Coefficient2 *= TmpNormalizationCoeffients[(*TmpMIndices2)];		    
			    (*TmpMIndices) = (*TmpMIndices2);			
			    ++TmpMIndices;
			    ++TmpMIndices2;
			  }
			  
			TmpInteraction[j] = Coefficient * Coefficient2 * UltimateCoef;
			//cout << TmpInteraction[j]<<endl;
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
	      for (int i = 0 ; i< this->NbrGroups; i++)
		delete[] TmpInteractionCoeffients[i];
	      
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



double ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3,int m4, int n1,int n2, int n3, int n4, double ** TmpInteractionCoeffients)
{
  double UltimateCoef=0.0;
  if ( m1+m2 == n1+n2) 
  {
	//cout <<"first if"<<endl;
	UltimateCoef += TmpInteractionCoeffients[0][m1+m2] * TmpInteractionCoeffients[0][m3+m4];
      }
      
     if  ( m1+m3 == n1+n3)
     {
	//cout <<"second if"<<endl;
	UltimateCoef += TmpInteractionCoeffients[0][m1+m3] * TmpInteractionCoeffients[0][m2+m4];
    }
		      
	if ( m1+m4 == n1+n4) 
	    {
	//cout <<"third if"<<endl;
	UltimateCoef += TmpInteractionCoeffients[0][m1+m4] * TmpInteractionCoeffients[0][m2+m3];
	    }
	    
  return UltimateCoef;	  
}