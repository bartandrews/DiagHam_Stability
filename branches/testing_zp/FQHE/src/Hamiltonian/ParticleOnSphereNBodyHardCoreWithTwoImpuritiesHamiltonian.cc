////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//  n-body hard core interaction and two localized impurities at the poles    //
//                                                                            //
//                        last modification : 04/05/2006                      //
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
#include "Hamiltonian/ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Vector/RealVector.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

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
// nbrBody = number of particle that interact simultaneously through the hard core interaction
// northPoleImpurityPotential = potential strength associted to the impurity at the north pole
// southPoleImpurityPotential = potential strength associted to the impurity at the south pole
// landauLevel = index of the Landau level where calculations have to be done (0 being the LLL)
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian::ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
														     int nbrBody, double northPoleImpurityPotential, 
														     double southPoleImpurityPotential, int landauLevel, double l2Factor,
														     AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag,
														     char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

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
  this->MNNBodyInteractionFactors = 0;

  this->OneBodyTermFlag = true;
  this->FullTwoBodyFlag = false;
  this->NorthPoleImpurityPotential = northPoleImpurityPotential;
  this->SouthPoleImpurityPotential = southPoleImpurityPotential;
  this->LandauLevel = landauLevel;

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 1;
      this->MaxSumIndices[k] = 0;      
      this->NBodyFlags[k] = false;
      this->NBodyInteractionWeightFactors[k] = 0.0;
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	this->NBodySign[k] = -1.0;
    }
  this->NBodyFlags[this->NbrNbody] = true;
  this->NBodyInteractionWeightFactors[this->NbrNbody] = 1.0;
  this->Architecture = architecture;
  this->EvaluateInteractionFactorsWithTwoImpurities();
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
      this->L2Operator = new ParticleOnSphereSquareTotalMomentumOperator(this->Particles, this->LzMax, l2Factor);
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
// maxNbrBody = maximum number of particle that interact simultaneously through the hard core interaction
// nBodyFactors = weight of the different n-body interaction terms with respect to each other
// northPoleImpurityPotential = potential strength associted to the impurity at the north pole
// southPoleImpurityPotential = potential strength associted to the impurity at the south pole
// landauLevel = index of the Landau level where calculations have to be done (0 being the LLL)
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian::ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
														     int maxNbrBody, double* nBodyFactors, 
														     double northPoleImpurityPotential, double southPoleImpurityPotential, 
														     int landauLevel, double l2Factor,
														     AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
														     char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

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

  this->OneBodyTermFlag = true;
  this->FullTwoBodyFlag = false;
  this->NorthPoleImpurityPotential = northPoleImpurityPotential;
  this->SouthPoleImpurityPotential = southPoleImpurityPotential;
  this->LandauLevel = landauLevel;

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
  this->EvaluateInteractionFactorsWithTwoImpurities();
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
      this->L2Operator = new ParticleOnSphereSquareTotalMomentumOperator(this->Particles, this->LzMax, l2Factor);
    }
  else
    {
      this->L2Operator = 0;
    }
}

// destructor
//

ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian::~ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian()
{
  if (this->OneBodyTermFlag == true)
    {
      delete[] this->OneBodyMValues;
      delete[] this->OneBodyNValues;
      delete[] this->OneBodyInteractionFactors;
    }
}

// evaluate all interaction factors (including those arising from impurities)
//   

void ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian::EvaluateInteractionFactorsWithTwoImpurities()
{
  this->EvaluateInteractionFactors();

  this->NbrOneBodyInteractionFactors = 2;
  this->OneBodyMValues = new int[this->NbrOneBodyInteractionFactors];
  this->OneBodyNValues = new int[this->NbrOneBodyInteractionFactors];
  this->OneBodyInteractionFactors = new double[this->NbrOneBodyInteractionFactors];
  ParticleOnSphereGenericLLFunctionBasis Basis (this->LzMax - (2 * this->LandauLevel), this->LandauLevel);
  RealVector Value(2, true);
  Complex TmpValue;
  Value[0] = M_PI;
  Basis.GetFunctionValue(Value, TmpValue, this->LandauLevel);
  this->OneBodyMValues[0] = this->LandauLevel;
  this->OneBodyNValues[0] = this->LandauLevel;
  this->OneBodyInteractionFactors[0] = -this->SouthPoleImpurityPotential * SqrNorm(TmpValue);
  Value[0] = 0.0;
  Basis.GetFunctionValue(Value, TmpValue, this->LzMax - this->LandauLevel);
  this->OneBodyMValues[1] = this->LzMax - this->LandauLevel;
  this->OneBodyNValues[1] = this->LzMax - this->LandauLevel;
  this->OneBodyInteractionFactors[1] = -this->NorthPoleImpurityPotential * SqrNorm(TmpValue);  
}



