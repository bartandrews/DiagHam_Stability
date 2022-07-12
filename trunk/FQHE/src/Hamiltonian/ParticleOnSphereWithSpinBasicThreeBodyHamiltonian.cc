////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                     spin and generic 3-body interaction                    //
//                                                                            //
//                        last modification : 26/08/2008                      //
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
#include "Hamiltonian/ParticleOnSphereWithSpinBasicThreeBodyHamiltonian.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"

 
#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//

ParticleOnSphereWithSpinBasicThreeBodyHamiltonian::ParticleOnSphereWithSpinBasicThreeBodyHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// threeBodyPseudoPotential = array with the three-body pseudo-potentials sorted with respect to the relative angular momentum,
//                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
//                            first index is the spin sector (0 up-up-up, 1 down-down-down, 2 up-up-down, 3 up-down-down)
// maxRelativeAngularMomentum =  maximum relative angular momentum that is used in ThreeBodyPseudoPotential  for each spin sector
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereWithSpinBasicThreeBodyHamiltonian::ParticleOnSphereWithSpinBasicThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax,
                                                                                                         double** threeBodyPseudoPotential, int* maxRelativeAngularMomentum,
                                                                                                         AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag,
                                                                                                         char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;

  this->OneBodyTermFlag = false;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->FullTwoBodyFlag = false;
  this->NbrIntraSectorSums = 0;
  this->NbrInterSectorSums = 0;
  this->M1IntraValue = 0;
  this->M1InterValue = 0;
  this->MaxNBody = 3;

  this->NBodyFlags = new bool [this->MaxNBody + 1];

  this->NbrSortedIndicesPerSum = new int** [this->MaxNBody + 1];
  this->SortedIndicesPerSum = new int*** [this->MaxNBody + 1];
  this->MinSumIndices = new int* [this->MaxNBody + 1];
  this->MaxSumIndices = new int* [this->MaxNBody + 1];

  this->NBodySign = new double*[this->MaxNBody + 1];
  this->SpinIndices = new int** [this->MaxNBody + 1];
  this->SpinIndicesShort = new int* [this->MaxNBody + 1];
  this->NbrNIndices = new long*[this->MaxNBody + 1];
  this->NIndices = new int**[this->MaxNBody + 1];
  this->NbrMIndices = new long**[this->MaxNBody + 1];
  this->MIndices = new int***[this->MaxNBody + 1];
  this->MNNBodyInteractionFactors = new double***[this->MaxNBody + 1];

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 0;
      this->MaxSumIndices[k] = 0;      
      this->NBodyFlags[k] = false;
      this->NBodySign[k] = 0;
      this->NbrNIndices[k] = 0;
      this->NIndices[k] = 0;
      this->NbrMIndices[k] = 0;
      this->MIndices[k] = 0;
      this->MNNBodyInteractionFactors[k] = 0;
      this->SpinIndices[k] = 0;
    }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NBodySign[3] = new double[4];
      this->NBodySign[3][0] = -1.0;
      this->NBodySign[3][1] = -1.0;
      this->NBodySign[3][2] = 1.0;
      this->NBodySign[3][3] = 1.0;
    }
  this->SpinIndices[3] = new int*[4];
  this->SpinIndices[3][0] = new int[3];
  this->SpinIndices[3][1] = new int[3];
  this->SpinIndices[3][2] = new int[3];
  this->SpinIndices[3][3] = new int[3];
  this->SpinIndices[3][0][0] = 0;
  this->SpinIndices[3][0][1] = 0;
  this->SpinIndices[3][0][2] = 0;
  this->SpinIndices[3][1][0] = 1;
  this->SpinIndices[3][1][1] = 1;
  this->SpinIndices[3][1][2] = 1;
  this->SpinIndices[3][2][0] = 0;
  this->SpinIndices[3][2][1] = 0;
  this->SpinIndices[3][2][2] = 1;
  this->SpinIndices[3][3][0] = 1;
  this->SpinIndices[3][3][1] = 1;
  this->SpinIndices[3][3][2] = 0;
  this->SpinIndicesShort[3] = new int[4];
  this->SpinIndicesShort[3][0] = 0x0;
  this->SpinIndicesShort[3][1] = 0x7;
  this->SpinIndicesShort[3][2] = 0x4;
  this->SpinIndicesShort[3][3] = 0x3;
 
  this->NbrSortedIndicesPerSum[3] = new int* [2];
  this->SortedIndicesPerSum[3] = new int** [2];
  this->MinSumIndices[3] = new int [2];
  this->MaxSumIndices[3] = new int [2];
  this->NBodyFlags[3] = true;
  this->NbrNIndices[3] = new long[2];
  this->NIndices[3] = new int*[2];
  this->NbrMIndices[3] = new long*[2];
  this->MIndices[3] = new int**[2];
  this->MNNBodyInteractionFactors[3] = new double**[4];
 

  this->MaxRelativeAngularMomentum = new int[4];
  for (int i = 0; i < 4; ++i)
    this->MaxRelativeAngularMomentum[i] = maxRelativeAngularMomentum[i] ;
  this->ThreeBodyPseudoPotentials = new double*[4];
  this->NbrThreeBodyPseudoPotentials = new int[4];
  for (int i = 0; i < 4; ++i)
    {
      this->NbrThreeBodyPseudoPotentials[i] = this->MaxRelativeAngularMomentum[i] + 1;
      this->ThreeBodyPseudoPotentials[i] = new double[this->NbrThreeBodyPseudoPotentials[i]];
      for (int k = 0; k < this->NbrThreeBodyPseudoPotentials[i]; ++k)
        this->ThreeBodyPseudoPotentials[i][k] = threeBodyPseudoPotential[i][k];
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

}

// constructor from datas with a fully-defined two body interaction
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// threeBodyPseudoPotential = array with the three-body pseudo-potentials sorted with respect to the relative angular momentum,
//                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
//                            first index is the spin sector (0 up-up-up, 1 down-down-down, 2 up-up-down, 3 up-down-down)
// maxRelativeAngularMomentum =  maximum relative angular momentum that is used in ThreeBodyPseudoPotential  for each spin sector
// pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereWithSpinBasicThreeBodyHamiltonian::ParticleOnSphereWithSpinBasicThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax,
                                                                                                         double** threeBodyPseudoPotential, int* maxRelativeAngularMomentum,
                                                                                                         double** pseudoPotential,
                                                                                                         double* onebodyPotentialUpUp, double* onebodyPotentialDownDown,
                                                                                                         AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag,
                                                                                                         char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;

  this->PseudoPotentials = new double* [3];
  for (int j = 0; j < 3; ++j)
    {
      this->PseudoPotentials[j] = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
        this->PseudoPotentials[j][i] = pseudoPotential[j][this->LzMax - i];
    }  
 
  if ((onebodyPotentialUpUp == 0) || (onebodyPotentialDownDown == 0))
    {
      this->OneBodyTermFlag = false;
      this->OneBodyInteractionFactorsupup = 0;
      this->OneBodyInteractionFactorsdowndown = 0;      
    }
  else
    {
      this->OneBodyTermFlag = true;
      this->OneBodyInteractionFactorsupup = new double [this->NbrLzValue];
      this->OneBodyInteractionFactorsdowndown = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
        {
          this->OneBodyInteractionFactorsupup[i] = onebodyPotentialUpUp[i];
          this->OneBodyInteractionFactorsdowndown[i] = onebodyPotentialDownDown[i];
        }
    }
  // don't have tunnelling implemented either way, here!
  this->OneBodyInteractionFactorsupdown = 0;
 
  this->FullTwoBodyFlag = true;

  this->MaxNBody = 3;

  this->NBodyFlags = new bool [this->MaxNBody + 1];

  this->NbrSortedIndicesPerSum = new int** [this->MaxNBody + 1];
  this->SortedIndicesPerSum = new int*** [this->MaxNBody + 1];
  this->MinSumIndices = new int* [this->MaxNBody + 1];
  this->MaxSumIndices = new int* [this->MaxNBody + 1];

  this->NBodySign = new double*[this->MaxNBody + 1];
  this->SpinIndices = new int** [this->MaxNBody + 1];
  this->SpinIndicesShort = new int* [this->MaxNBody + 1];
  this->NbrNIndices = new long*[this->MaxNBody + 1];
  this->NIndices = new int**[this->MaxNBody + 1];
  this->NbrMIndices = new long**[this->MaxNBody + 1];
  this->MIndices = new int***[this->MaxNBody + 1];
  this->MNNBodyInteractionFactors = new double***[this->MaxNBody + 1];

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 0;
      this->MaxSumIndices[k] = 0;      
      this->NBodyFlags[k] = false;
      this->NBodySign[k] = 0;
      this->NbrNIndices[k] = 0;
      this->NIndices[k] = 0;
      this->NbrMIndices[k] = 0;
      this->MIndices[k] = 0;
      this->MNNBodyInteractionFactors[k] = 0;
      this->SpinIndices[k] = 0;
    }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NBodySign[3] = new double[4];
      this->NBodySign[3][0] = -1.0;
      this->NBodySign[3][1] = -1.0;
      this->NBodySign[3][2] = 1.0;
      this->NBodySign[3][3] = 1.0;
    }
  this->SpinIndices[3] = new int*[4];
  this->SpinIndices[3][0] = new int[3];
  this->SpinIndices[3][1] = new int[3];
  this->SpinIndices[3][2] = new int[3];
  this->SpinIndices[3][3] = new int[3];
  this->SpinIndices[3][0][0] = 0;
  this->SpinIndices[3][0][1] = 0;
  this->SpinIndices[3][0][2] = 0;
  this->SpinIndices[3][1][0] = 1;
  this->SpinIndices[3][1][1] = 1;
  this->SpinIndices[3][1][2] = 1;
  this->SpinIndices[3][2][0] = 0;
  this->SpinIndices[3][2][1] = 0;
  this->SpinIndices[3][2][2] = 1;
  this->SpinIndices[3][3][0] = 1;
  this->SpinIndices[3][3][1] = 1;
  this->SpinIndices[3][3][2] = 0;
  this->SpinIndicesShort[3] = new int[4];
  this->SpinIndicesShort[3][0] = 0x0;
  this->SpinIndicesShort[3][1] = 0x7;
  this->SpinIndicesShort[3][2] = 0x4;
  this->SpinIndicesShort[3][3] = 0x3;

  this->NbrSortedIndicesPerSum[3] = new int* [2];
  this->SortedIndicesPerSum[3] = new int** [2];
  this->MinSumIndices[3] = new int [2];
  this->MaxSumIndices[3] = new int [2];
  this->NBodyFlags[3] = true;
  this->NbrNIndices[3] = new long[2];
  this->NIndices[3] = new int*[2];
  this->NbrMIndices[3] = new long*[2];
  this->MIndices[3] = new int**[2];
  this->MNNBodyInteractionFactors[3] = new double**[4];
 

  this->MaxRelativeAngularMomentum = new int[4];
  for (int i = 0; i < 4; ++i)
    this->MaxRelativeAngularMomentum[i] = maxRelativeAngularMomentum[i];
  this->ThreeBodyPseudoPotentials = new double*[4];
  this->NbrThreeBodyPseudoPotentials = new int[4];
  for (int i = 0; i < 4; ++i)
    {
      this->NbrThreeBodyPseudoPotentials[i] = this->MaxRelativeAngularMomentum[i] + 1;
      this->ThreeBodyPseudoPotentials[i] = new double[this->NbrThreeBodyPseudoPotentials[i]];
      for (int k = 0; k < this->NbrThreeBodyPseudoPotentials[i]; ++k)
        this->ThreeBodyPseudoPotentials[i][k] = threeBodyPseudoPotential[i][k];
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
}

// destructor
//

ParticleOnSphereWithSpinBasicThreeBodyHamiltonian::~ParticleOnSphereWithSpinBasicThreeBodyHamiltonian()
{
  for (int i = 0; i < 4; ++i)
    delete[] this->ThreeBodyPseudoPotentials[i];
  delete[] this->MaxRelativeAngularMomentum;
  delete[] this->NbrThreeBodyPseudoPotentials;
  delete[] this->ThreeBodyPseudoPotentials;
  if (this->FullTwoBodyFlag == true)
    {
      for (int i = 0; i < 3; ++i)
        delete[] this->PseudoPotentials[i];
      delete[] this->PseudoPotentials;
    }

  delete [] this->NBodyFlags;

  // rough job deleting only outer arrays...
  delete [] this->NbrSortedIndicesPerSum;
  delete [] this->SortedIndicesPerSum;
  delete [] this->MinSumIndices;
  delete [] this->MaxSumIndices;

  delete [] this->NBodySign;
  delete [] this->SpinIndices;
  delete [] this->SpinIndicesShort;
  delete [] this->NbrNIndices;
  delete [] this->NIndices; // multiple
  delete [] this->NbrMIndices; // multiple = new long**[this->MaxNBody + 1];
  delete [] this->MIndices; // multiple = new int***[this->MaxNBody + 1];
  delete [] this->MNNBodyInteractionFactors; // multiple = new double***[this->MaxNBody + 1];

}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* ParticleOnSphereWithSpinBasicThreeBodyHamiltonian::Clone ()
{
  return 0;
}

// evaluate all interaction factors
//  

void    ParticleOnSphereWithSpinBasicThreeBodyHamiltonian::EvaluateInteractionFactors()
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
      GetAllSkewSymmetricIndices(this->NbrLzValue, 3, this->NbrSortedIndicesPerSum[3][0], this->SortedIndicesPerSum[3][0]);
      this->MaxSumIndices[3][0] = ((this->NbrLzValue - 1) * 3) - 3;
      this->MinSumIndices[3][0] = 3;
      double* TmpInteractionCoeffients = new double[this->MaxSumIndices[3][0] + 1];
      TmpInteractionCoeffients[0] = 1.0;
      TmpInteractionCoeffients[1] = 1.0;
      for (int i = 2; i <= this->MaxSumIndices[3][0]; ++i)
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
      Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[3][0]) + 1.0);
      double Radius = 2.0 / ((double) this->LzMax);
      for (int i = 2; i <= 3; ++i)
        {
          Coefficient *= (double) (i * i);       
          Coefficient *= Radius;
        }
      for (int i = 0; i <= this->MaxSumIndices[3][0]; ++i)
        TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);
     
      long TmpNbrNIndices = 0;
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[3][0]; ++MinSum)
        TmpNbrNIndices += this->NbrSortedIndicesPerSum[3][0][MinSum];
      this->NbrNIndices[3][0] = TmpNbrNIndices;
      this->NIndices[3][0] = new int[TmpNbrNIndices * 3];
      this->NbrMIndices[3][0] = new long[TmpNbrNIndices];
      this->MIndices[3][0] = new int*[TmpNbrNIndices];
      this->MNNBodyInteractionFactors[3][0] = new double* [TmpNbrNIndices];
      this->MNNBodyInteractionFactors[3][1] = new double* [TmpNbrNIndices];
      TmpNbrNIndices = 0;        
      int* TmpNIndices = this->NIndices[3][0];
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[3][0]; ++MinSum)
        {
          int Lim = this->NbrSortedIndicesPerSum[3][0][MinSum];
          if (Lim > 0)
            {
              int* TmpNIndices2 = this->SortedIndicesPerSum[3][0][MinSum];
              int TmpMaxRelativeMomentum = 8;
              if (this->MaxRelativeAngularMomentum[0] <= TmpMaxRelativeMomentum)
                TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum[0];
              if (this->MaxRelativeAngularMomentum[1] > TmpMaxRelativeMomentum)
                {
                  if (this->MaxRelativeAngularMomentum[1] <= 8)
                    TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum[1];
                  else
                    TmpMaxRelativeMomentum = 8;            
                }
              int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2];
              while (((3 * this->LzMax) - TmpMaxRelativeMomentum)  < TmpSum)
                --TmpMaxRelativeMomentum;
              double** TmpProjectorCoefficients = new double* [TmpMaxRelativeMomentum + 1];
              if ((this->ThreeBodyPseudoPotentials[0][3] != 0.0) || (this->ThreeBodyPseudoPotentials[1][3] != 0.0))
                TmpProjectorCoefficients[3] = this->ComputeProjectorCoefficients(6, 1, TmpNIndices2, Lim, (2 * this->LzMax) - 2, 0);
              for (int i = 5; i <= TmpMaxRelativeMomentum; ++i)  
                if ((this->ThreeBodyPseudoPotentials[0][i] != 0.0) || (this->ThreeBodyPseudoPotentials[1][i] != 0.0))
                  TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim, (2 * this->LzMax) - 2, 0);
              for (int i = 0; i < Lim; ++i)
                {
                  this->NbrMIndices[3][0][TmpNbrNIndices] = Lim;                   
                  this->MIndices[3][0][TmpNbrNIndices] = new int [Lim * 3];
                  this->MNNBodyInteractionFactors[3][0][TmpNbrNIndices] = new double [Lim];
                  this->MNNBodyInteractionFactors[3][1][TmpNbrNIndices] = new double [Lim];
                  int* TmpMIndices = this->MIndices[3][0][TmpNbrNIndices];
                  int* TmpMIndices2 = this->SortedIndicesPerSum[3][0][MinSum];
                  double* TmpInteractionUpUpUp = this->MNNBodyInteractionFactors[3][0][TmpNbrNIndices];
                  double* TmpInteractionDownDownDown = this->MNNBodyInteractionFactors[3][1][TmpNbrNIndices];
                  for (int j = 0; j < Lim; ++j)
                    {
                      for (int l = 0; l < 3; ++l)
                        {
                          (*TmpMIndices) = (*TmpMIndices2);                    
                          ++TmpMIndices;
                          ++TmpMIndices2;
                        }                      
                      double& TmpInteraction2UpUp = TmpInteractionUpUpUp[j];
                      TmpInteraction2UpUp = 0.0;
                      double& TmpInteraction2DownDown = TmpInteractionDownDownDown[j];
                      TmpInteraction2DownDown = 0.0;
                      if (this->ThreeBodyPseudoPotentials[0][3] != 0.0)
                        TmpInteraction2UpUp += 4.0 * this->ThreeBodyPseudoPotentials[0][3] * TmpProjectorCoefficients[3][i] * TmpProjectorCoefficients[3][j];
                      if (this->ThreeBodyPseudoPotentials[1][3] != 0.0)
                        TmpInteraction2DownDown += 4.0 * this->ThreeBodyPseudoPotentials[1][3] * TmpProjectorCoefficients[3][i] * TmpProjectorCoefficients[3][j];
                      for (int k = 5; k <= TmpMaxRelativeMomentum; ++k)  
                        {
                          if (this->ThreeBodyPseudoPotentials[0][k] != 0.0)
                            TmpInteraction2UpUp += 4.0 * this->ThreeBodyPseudoPotentials[0][k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
                          if (this->ThreeBodyPseudoPotentials[1][k] != 0.0)
                            TmpInteraction2DownDown += 4.0 * this->ThreeBodyPseudoPotentials[1][k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
                        }
                    }
                  for (int j = 0; j < 3; ++j)
                    {
                      (*TmpNIndices) = (*TmpNIndices2);                
                      ++TmpNIndices;
                      ++TmpNIndices2;
                    }
                  ++TmpNbrNIndices;
                }
              if ((this->ThreeBodyPseudoPotentials[0][3] != 0.0) || (this->ThreeBodyPseudoPotentials[1][3]))
                delete[] TmpProjectorCoefficients[3];
              for (int i = 5; i <= TmpMaxRelativeMomentum; ++i)  
                if ((this->ThreeBodyPseudoPotentials[0][i] != 0.0) || (this->ThreeBodyPseudoPotentials[1][i]))
                  delete[] TmpProjectorCoefficients[i];
              delete[] TmpProjectorCoefficients;               
            }
        }
      delete[] TmpInteractionCoeffients;


      GetAllTwoSetSkewSymmetricIndices(this->NbrLzValue, 2, 1, this->NbrSortedIndicesPerSum[3][1], this->SortedIndicesPerSum[3][1]);
      this->MaxSumIndices[3][1] = this->LzMax * 3 - 1;
      this->MinSumIndices[3][1] = 0;
      TmpInteractionCoeffients = new double[this->MaxSumIndices[3][1] + 1];
      TmpInteractionCoeffients[0] = 1.0;
      TmpInteractionCoeffients[1] = 1.0;
      for (int i = 2; i <= this->MaxSumIndices[3][1]; ++i)
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
      Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[3][1]) + 1.0);
      for (int i = 2; i <= 3; ++i)
        {
          Coefficient *= (double) (i * i);       
          Coefficient *= Radius;
        }
      for (int i = 0; i <= this->MaxSumIndices[3][1]; ++i)
        TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);
     
      TmpNbrNIndices = 0;
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[3][1]; ++MinSum)
        TmpNbrNIndices += this->NbrSortedIndicesPerSum[3][1][MinSum];
      this->NbrNIndices[3][1] = TmpNbrNIndices;
      this->NIndices[3][1] = new int[TmpNbrNIndices * 3];
      this->NbrMIndices[3][1] = new long[TmpNbrNIndices];
      this->MIndices[3][1] = new int*[TmpNbrNIndices];
      this->MNNBodyInteractionFactors[3][2] = new double* [TmpNbrNIndices];
      this->MNNBodyInteractionFactors[3][3] = new double* [TmpNbrNIndices];
      TmpNbrNIndices = 0;        
      TmpNIndices = this->NIndices[3][1];
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[3][1]; ++MinSum)
        {
          int Lim = this->NbrSortedIndicesPerSum[3][1][MinSum];
          if (Lim > 0)
            {
              int* TmpNIndices2 = this->SortedIndicesPerSum[3][1][MinSum];
              int TmpMaxRelativeMomentum = 6;
              if (this->MaxRelativeAngularMomentum[2] <= TmpMaxRelativeMomentum)
                TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum[2];
              if (this->MaxRelativeAngularMomentum[3] > TmpMaxRelativeMomentum)
                {
                  if (this->MaxRelativeAngularMomentum[3] <= 6)
                    TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum[3];
                  else
                    TmpMaxRelativeMomentum = 6;            
                }
              int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2];
              while (((3 * this->LzMax) - TmpMaxRelativeMomentum)  < TmpSum)
                --TmpMaxRelativeMomentum;
              double** TmpProjectorCoefficients = new double* [TmpMaxRelativeMomentum + 1];
              for (int i = 1; i <= TmpMaxRelativeMomentum; ++i)  
                if ((this->ThreeBodyPseudoPotentials[2][i] != 0.0) || (this->ThreeBodyPseudoPotentials[3][i] != 0.0))
                  TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim, 2 * this->LzMax - 2, 1);
             
              int TmpNbrNIndices3 = TmpNbrNIndices;
              for (int i = 0; i < Lim; ++i)
                {
                  this->NbrMIndices[3][1][TmpNbrNIndices] = Lim;                   
                  this->MIndices[3][1][TmpNbrNIndices] = new int [Lim * 3];
                  this->MNNBodyInteractionFactors[3][2][TmpNbrNIndices] = new double [Lim];
                  this->MNNBodyInteractionFactors[3][3][TmpNbrNIndices] = new double [Lim];
                  int* TmpMIndices = this->MIndices[3][1][TmpNbrNIndices];
                  int* TmpMIndices2 = this->SortedIndicesPerSum[3][1][MinSum];
                  double* TmpInteractionUpUpDown = this->MNNBodyInteractionFactors[3][2][TmpNbrNIndices];
                  double* TmpInteractionDownDownUp = this->MNNBodyInteractionFactors[3][3][TmpNbrNIndices];
                  for (int j = 0; j < Lim; ++j)
                    {
                      for (int l = 0; l < 3; ++l)
                        {
                          (*TmpMIndices) = (*TmpMIndices2);                    
                          ++TmpMIndices;
                          ++TmpMIndices2;
                        }                      
                      double& TmpInteraction2UpUpDown = TmpInteractionUpUpDown[j];
                      TmpInteraction2UpUpDown = 0.0;
                      double& TmpInteraction2DownDownUp = TmpInteractionDownDownUp[j];
                      TmpInteraction2DownDownUp = 0.0;
                      for (int k = 0; k <= TmpMaxRelativeMomentum; ++k)  
                        {
                          if (this->ThreeBodyPseudoPotentials[2][k] != 0.0)
                            TmpInteraction2UpUpDown += this->ThreeBodyPseudoPotentials[2][k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
                          if (this->ThreeBodyPseudoPotentials[3][k] != 0.0)
                            TmpInteraction2DownDownUp += this->ThreeBodyPseudoPotentials[3][k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
                        }
                    }
                  for (int j = 0; j < 3; ++j)
                    {
                      (*TmpNIndices) = (*TmpNIndices2);                
                      ++TmpNIndices;
                      ++TmpNIndices2;
                    }
                  ++TmpNbrNIndices;
                }
              for (int i = 0; i <= TmpMaxRelativeMomentum; ++i)  
                if ((this->ThreeBodyPseudoPotentials[2][i] != 0.0) || (this->ThreeBodyPseudoPotentials[3][i] != 0.0))
                  delete[] TmpProjectorCoefficients[i];


              TmpNIndices2 = this->SortedIndicesPerSum[3][1][MinSum];
              TmpNbrNIndices = TmpNbrNIndices3;
              for (int i = 1; i <= TmpMaxRelativeMomentum; ++i)  
                if ((this->ThreeBodyPseudoPotentials[2][i] != 0.0) || (this->ThreeBodyPseudoPotentials[3][i] != 0.0))
                  TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim, 2 * this->LzMax, 2);
              for (int i = 0; i < Lim; ++i)
                {
                  double* TmpInteractionUpUpDown = this->MNNBodyInteractionFactors[3][2][TmpNbrNIndices];
                  double* TmpInteractionDownDownUp = this->MNNBodyInteractionFactors[3][3][TmpNbrNIndices];
                  for (int j = 0; j < Lim; ++j)
                    {
                      double& TmpInteraction2UpUpDown = TmpInteractionUpUpDown[j];
                      double& TmpInteraction2DownDownUp = TmpInteractionDownDownUp[j];
                      for (int k = 0; k <= TmpMaxRelativeMomentum; ++k)  
                        {
                          if (this->ThreeBodyPseudoPotentials[2][k] != 0.0)
                            TmpInteraction2UpUpDown += this->ThreeBodyPseudoPotentials[2][k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
                          if (this->ThreeBodyPseudoPotentials[3][k] != 0.0)
                            TmpInteraction2DownDownUp += this->ThreeBodyPseudoPotentials[3][k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
                        }
                    }
                  ++TmpNbrNIndices;
                }
              for (int i = 0; i <= TmpMaxRelativeMomentum; ++i)  
                if ((this->ThreeBodyPseudoPotentials[2][i] != 0.0) || (this->ThreeBodyPseudoPotentials[3][i] != 0.0))
                  delete[] TmpProjectorCoefficients[i];


              TmpNIndices2 = this->SortedIndicesPerSum[3][1][MinSum];
              TmpNbrNIndices = TmpNbrNIndices3;
              for (int i = 1; i <= TmpMaxRelativeMomentum; ++i)  
                if ((this->ThreeBodyPseudoPotentials[2][i] != 0.0) || (this->ThreeBodyPseudoPotentials[3][i] != 0.0))
                  TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim, 2 * this->LzMax, 3);
              for (int i = 0; i < Lim; ++i)
                {
                  double* TmpInteractionUpUpDown = this->MNNBodyInteractionFactors[3][2][TmpNbrNIndices];
                  double* TmpInteractionDownDownUp = this->MNNBodyInteractionFactors[3][3][TmpNbrNIndices];
                  for (int j = 0; j < Lim; ++j)
                    {
                      double& TmpInteraction2UpUpDown = TmpInteractionUpUpDown[j];
                      double& TmpInteraction2DownDownUp = TmpInteractionDownDownUp[j];
                      for (int k = 0; k <= TmpMaxRelativeMomentum; ++k)  
                        {
                          if (this->ThreeBodyPseudoPotentials[2][k] != 0.0)
                            TmpInteraction2UpUpDown += this->ThreeBodyPseudoPotentials[2][k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
                          if (this->ThreeBodyPseudoPotentials[3][k] != 0.0)
                            TmpInteraction2DownDownUp += this->ThreeBodyPseudoPotentials[3][k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
                        }
                    }
                  ++TmpNbrNIndices;
                }
              for (int i = 0; i <= TmpMaxRelativeMomentum; ++i)  
                if ((this->ThreeBodyPseudoPotentials[2][i] != 0.0) || (this->ThreeBodyPseudoPotentials[3][i] != 0.0))
                  delete[] TmpProjectorCoefficients[i];


              delete[] TmpProjectorCoefficients;               
            }
        }
      delete[] TmpInteractionCoeffients;
    }
  else
    { // xxx - case of bosons
      this->MinSumIndices[3][0] = 0;
      this->MaxSumIndices[3][0] = this->LzMax * 3;
      double* TmpInteractionCoeffients = new double[this->MaxSumIndices[3][0] + 1];
      double Coefficient;
      TmpInteractionCoeffients[0] = 1.0;
      TmpInteractionCoeffients[1] = 1.0;
      for (int i = 2; i <= this->MaxSumIndices[3][0]; ++i)
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
      Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[3][0]) + 1.0);
      double Radius = 2.0 / ((double) this->LzMax);
      for (int i = 2; i <= 3; ++i)
        {
          Coefficient *= (double) (i * i);       
          Coefficient *= Radius;
        }
      for (int i = 0; i <= this->MaxSumIndices[3][0]; ++i)
        TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);
     
      double** SortedIndicesPerSumSymmetryFactor;
      GetAllSymmetricIndices(this->NbrLzValue, 3, this->NbrSortedIndicesPerSum[3][0], this->SortedIndicesPerSum[3][0],
                             SortedIndicesPerSumSymmetryFactor);
     
     
      long TmpNbrNIndices = 0;
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[3][0]; ++MinSum)
        TmpNbrNIndices += this->NbrSortedIndicesPerSum[3][0][MinSum];
      this->NbrNIndices[3][0] = TmpNbrNIndices;
      this->NIndices[3][0] = new int[TmpNbrNIndices * 3];
      this->NbrMIndices[3][0] = new long[TmpNbrNIndices];
      this->MIndices[3][0] = new int*[TmpNbrNIndices];
      this->MNNBodyInteractionFactors[3][0] = new double* [TmpNbrNIndices];
      this->MNNBodyInteractionFactors[3][1] = new double* [TmpNbrNIndices];
      TmpNbrNIndices = 0;        
      int* TmpNIndices = this->NIndices[3][0];
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[3][0]; ++MinSum)
        {
          int Lim = this->NbrSortedIndicesPerSum[3][0][MinSum];
          double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
          int* TmpNIndices2 = this->SortedIndicesPerSum[3][0][MinSum];
          int TmpMaxRelativeMomentum = 10;
          if (this->MaxRelativeAngularMomentum[0] <= TmpMaxRelativeMomentum)
            TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum[0];
          int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2];
          while (((3 * this->LzMax) - TmpMaxRelativeMomentum)  < TmpSum)
            --TmpMaxRelativeMomentum;
          double** TmpProjectorCoefficients = new double* [TmpMaxRelativeMomentum + 1];
          if ((this->ThreeBodyPseudoPotentials[0][0] != 0.0)||(this->ThreeBodyPseudoPotentials[1][0] != 0.0))
            TmpProjectorCoefficients[0] = this->ComputeProjectorCoefficients(0, 1, TmpNIndices2, Lim, 2 * this->LzMax, 0);
          for (int i = 2; i <= TmpMaxRelativeMomentum; ++i)  
            if ((this->ThreeBodyPseudoPotentials[0][i] != 0.0)||(this->ThreeBodyPseudoPotentials[1][i] != 0.0))
              TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim, 2 * this->LzMax, 0);
          for (int i = 0; i < Lim; ++i)
            {
              this->NbrMIndices[3][0][TmpNbrNIndices] = Lim;               
              this->MIndices[3][0][TmpNbrNIndices] = new int [Lim * 3];
              this->MNNBodyInteractionFactors[3][0][TmpNbrNIndices] = new double [Lim];
              this->MNNBodyInteractionFactors[3][1][TmpNbrNIndices] = new double [Lim];
              int* TmpMIndices = this->MIndices[3][0][TmpNbrNIndices];
              int* TmpMIndices2 = this->SortedIndicesPerSum[3][0][MinSum];
              double* TmpInteractionUpUpUp = this->MNNBodyInteractionFactors[3][0][TmpNbrNIndices];
              double* TmpInteractionDownDownDown = this->MNNBodyInteractionFactors[3][1][TmpNbrNIndices];
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
                  double& TmpInteraction2UpUp = TmpInteractionUpUpUp[j];
                  TmpInteraction2UpUp = 0.0;
                  if (this->ThreeBodyPseudoPotentials[0][0] != 0.0)
                    TmpInteraction2UpUp += this->ThreeBodyPseudoPotentials[0][0] * TmpProjectorCoefficients[0][i] * TmpProjectorCoefficients[0][j];
                  for (int k = 2; k <= TmpMaxRelativeMomentum; ++k)
                    if (this->ThreeBodyPseudoPotentials[0][k] != 0.0)
                      TmpInteraction2UpUp += this->ThreeBodyPseudoPotentials[0][k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
                  TmpInteraction2UpUp *= TmpSymmetryFactors[i] * TmpSymmetryFactors[j];

                  double& TmpInteraction2DownDown = TmpInteractionDownDownDown[j];
                  TmpInteraction2DownDown = 0.0;
                  if (this->ThreeBodyPseudoPotentials[1][0] != 0.0)
                    TmpInteraction2DownDown += this->ThreeBodyPseudoPotentials[1][0] * TmpProjectorCoefficients[0][i] * TmpProjectorCoefficients[0][j];
                  for (int k = 2; k <= TmpMaxRelativeMomentum; ++k)
                    if (this->ThreeBodyPseudoPotentials[1][k] != 0.0)
                      TmpInteraction2DownDown += this->ThreeBodyPseudoPotentials[1][k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
                  TmpInteraction2DownDown *= TmpSymmetryFactors[i] * TmpSymmetryFactors[j];
                }
              for (int j = 0; j < 3; ++j)
                {
                  (*TmpNIndices) = (*TmpNIndices2);
                  ++TmpNIndices;
                  ++TmpNIndices2;
                }
              ++TmpNbrNIndices;
            }
          if ((this->ThreeBodyPseudoPotentials[0][0] != 0.0)||(this->ThreeBodyPseudoPotentials[1][0] != 0.0))
            delete[] TmpProjectorCoefficients[0];
          for (int i = 2; i <= TmpMaxRelativeMomentum; ++i)  
            if ((this->ThreeBodyPseudoPotentials[0][i] != 0.0)||(this->ThreeBodyPseudoPotentials[1][i] != 0.0))
              delete[] TmpProjectorCoefficients[i];
          delete[] TmpProjectorCoefficients;           
        }
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[3][0]; ++MinSum)
        {
          delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
        }
      delete[] SortedIndicesPerSumSymmetryFactor;
      delete[] TmpInteractionCoeffients;

      this->NbrNIndices[3][1] = 0;
      if ((this->NbrThreeBodyPseudoPotentials[2]>0)||(this->NbrThreeBodyPseudoPotentials[3]>0))
        cout << "No cross-spin terms implemented in Basic Three-Body Hamiltonian for Bosons!"<<endl;
    }
  delete[] TmpNormalizationCoeffients;
  if (this->FullTwoBodyFlag == true)
    {
      ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
      int J = 2 * this->LzMax - 2;
      double ClebschCoef;
      long TotalNbrInteractionFactors = 0;
     
      int Sign = 1;
      if (this->LzMax & 1)
        Sign = 0;
      double TmpCoefficient = 0.0;
     
      this->NbrInterSectorSums = 2 * this->LzMax + 1;
      this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
        this->NbrInterSectorIndicesPerSum[i] = 0;
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
        for (int m2 = 0; m2 <= this->LzMax; ++m2)
          ++this->NbrInterSectorIndicesPerSum[m1 + m2];      
      this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
        {
          this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
          this->NbrInterSectorIndicesPerSum[i] = 0;
        }
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
        for (int m2 = 0; m2 <= this->LzMax; ++m2)
          {
            this->InterSectorIndicesPerSum[(m1 + m2)][this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
            this->InterSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
            ++this->NbrInterSectorIndicesPerSum[(m1 + m2)];
          }
     
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
        {
          this->NbrIntraSectorSums = 2 * this->LzMax - 1;
          this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
          for (int i = 0; i < this->NbrIntraSectorSums; ++i)
            this->NbrIntraSectorIndicesPerSum[i] = 0;      
          for (int m1 = 0; m1 < this->LzMax; ++m1)
            for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
              ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
          this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
          for (int i = 0; i < this->NbrIntraSectorSums; ++i)
            {
              this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
              this->NbrIntraSectorIndicesPerSum[i] = 0;
            }
          for (int m1 = 0; m1 < this->LzMax; ++m1)
            for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
              {
                this->IntraSectorIndicesPerSum[(m1 + m2) - 1][this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1] = m1;
                this->IntraSectorIndicesPerSum[(m1 + m2) - 1][1 + (this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1)] = m2;
                ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
              }
         
          this->InteractionFactorsupup = new double* [this->NbrIntraSectorSums];
          this->InteractionFactorsdowndown = new double* [this->NbrIntraSectorSums];
          for (int i = 0; i < this->NbrIntraSectorSums; ++i)
            {
              this->InteractionFactorsupup[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
              this->InteractionFactorsdowndown[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
              int Index = 0;
              for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
                {
                  int m1 = (this->IntraSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
                  int m2 = (this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
                  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
                    {
                      int m3 = (this->IntraSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
                      int m4 = (this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
                      Clebsch.InitializeCoefficientIterator(m1, m2);
                      this->InteractionFactorsupup[i][Index] = 0.0;
                      this->InteractionFactorsdowndown[i][Index] = 0.0;
                      while (Clebsch.Iterate(J, ClebschCoef))
                        {
                          if (((J >> 1) & 1) == Sign)
                            {
                              TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
                              this->InteractionFactorsupup[i][Index] += this->PseudoPotentials[0][J >> 1] * TmpCoefficient;
                              this->InteractionFactorsdowndown[i][Index] += this->PseudoPotentials[1][J >> 1] * TmpCoefficient;
                            }
                        }
                      this->InteractionFactorsupup[i][Index] *= -4.0;
                      this->InteractionFactorsdowndown[i][Index] *= -4.0;
                      TotalNbrInteractionFactors += 2;
                      ++Index;
                    }
                }
            }
         
          this->InteractionFactorsupdown = new double* [this->NbrInterSectorSums];
          for (int i = 0; i < this->NbrInterSectorSums; ++i)
            {
              this->InteractionFactorsupdown[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
              int Index = 0;
              for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
                {
                  double Factor = 2.0;
                  int m1 = (this->InterSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
                  int m2 = (this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
                  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
                    {
                      int m3 = (this->InterSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
                      int m4 = (this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
                      Clebsch.InitializeCoefficientIterator(m1, m2);
                      this->InteractionFactorsupdown[i][Index] = 0.0;
                      while (Clebsch.Iterate(J, ClebschCoef))
                        {
                          TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
                          this->InteractionFactorsupdown[i][Index] += this->PseudoPotentials[2][J >> 1] * TmpCoefficient;
                        }
                      this->InteractionFactorsupdown[i][Index] *= -Factor;
                      ++TotalNbrInteractionFactors;
                      ++Index;
                    }
                }
            }
        }
      else
        { // two-body interaction for bosons needs to be tested
          this->NbrIntraSectorSums = 2 * this->LzMax + 1;
          this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
          for (int i = 0; i < this->NbrIntraSectorSums; ++i)
            this->NbrIntraSectorIndicesPerSum[i] = 0;
          for (int m1 = 0; m1 <= this->LzMax; ++m1)
            for (int m2 = m1; m2 <= this->LzMax; ++m2)
              ++this->NbrIntraSectorIndicesPerSum[m1 + m2];
          this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
          for (int i = 0; i < this->NbrIntraSectorSums; ++i)
            {
              this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
              this->NbrIntraSectorIndicesPerSum[i] = 0;
            }
          for (int m1 = 0; m1 <= this->LzMax; ++m1)
            for (int m2 = m1; m2 <= this->LzMax; ++m2)
              {
                this->IntraSectorIndicesPerSum[m1 + m2][this->NbrIntraSectorIndicesPerSum[m1 + m2] << 1] = m1;
                this->IntraSectorIndicesPerSum[m1 + m2][1 + (this->NbrIntraSectorIndicesPerSum[m1 + m2] << 1)] = m2;
                ++this->NbrIntraSectorIndicesPerSum[m1 + m2];
              }
         
          this->InteractionFactorsupup = new double* [this->NbrIntraSectorSums];
          this->InteractionFactorsdowndown = new double* [this->NbrIntraSectorSums];
          for (int i = 0; i < this->NbrIntraSectorSums; ++i)
            {
              this->InteractionFactorsupup[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
              this->InteractionFactorsdowndown[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
              int Index = 0;
              for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
                {
                  int m1 = (this->IntraSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
                  int m2 = (this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
                  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
                    {
                      int m3 = (this->IntraSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
                      int m4 = (this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
                      Clebsch.InitializeCoefficientIterator(m1, m2);
                      this->InteractionFactorsupup[i][Index] = 0.0;
                      this->InteractionFactorsdowndown[i][Index] = 0.0;
                      while (Clebsch.Iterate(J, ClebschCoef))
                        {
                          if (((J >> 1) & 1) != Sign)
                            {
                              TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
                              this->InteractionFactorsupup[i][Index] += this->PseudoPotentials[0][J >> 1] * TmpCoefficient;
                              this->InteractionFactorsdowndown[i][Index] += this->PseudoPotentials[1][J >> 1] * TmpCoefficient;
                            }
                        }
                      if (m1!=m2)
                        {
                          this->InteractionFactorsupup[i][Index] *= 2.0;
                          this->InteractionFactorsdowndown[i][Index] *= 2.0;
                        }
                      if (m3!=m4)
                        {
                          this->InteractionFactorsupup[i][Index] *= 2.0;
                          this->InteractionFactorsdowndown[i][Index] *= 2.0;
                        }
                      TotalNbrInteractionFactors += 2;
                      ++Index;
                    }
                }
            }
         
          this->InteractionFactorsupdown = new double* [this->NbrInterSectorSums];
          for (int i = 0; i < this->NbrInterSectorSums; ++i)
            {
              this->InteractionFactorsupdown[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
              int Index = 0;
              for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
                {
                  double Factor = 2.0;
                  int m1 = (this->InterSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
                  int m2 = (this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
                  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
                    {
                      int m3 = (this->InterSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
                      int m4 = (this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
                      Clebsch.InitializeCoefficientIterator(m1, m2);
                      this->InteractionFactorsupdown[i][Index] = 0.0;
                      while (Clebsch.Iterate(J, ClebschCoef))
                        {
                          TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
                          this->InteractionFactorsupdown[i][Index] += this->PseudoPotentials[2][J >> 1] * TmpCoefficient;
                        }
                      this->InteractionFactorsupdown[i][Index] *= Factor;
                      ++TotalNbrInteractionFactors;
                      ++Index;
                    }
                }
            }
        }
    }
}

// compute all projector coefficient associated to a given relative angular momentum between 3 particles
//
// relativeMomentum = value of twice the relative angular momentum between the 3 particles
// degeneracyIndex = optional degeneracy index for relative angular momentum greater than 5 for bosons (8 for fermions)
// indices = array that contains all possible sets of indices (size of the array is 3 * nbrIndexSets)
// nbrIndexSets = number of sets
// maxJValue = twice the maximum total angular momentum two particles can have
// spinIndex = indicate for which of three body operators coeeficients are computed (0 for up-up-up, 1 for up-up-down, 2 for up-down-up, 3 for down-up-up)

double* ParticleOnSphereWithSpinBasicThreeBodyHamiltonian::ComputeProjectorCoefficients(int relativeMomentum, int degeneracyIndex, int* indices, int nbrIndexSets, int maxJValue, int spinIndex)
{
  double* TmpCoefficients = new double [nbrIndexSets];
  int JValue = (3 * this->LzMax) - relativeMomentum;

  ClebschGordanCoefficients Clebsh (this->LzMax, this->LzMax);
  ClebschGordanCoefficients* ClebshArray = new ClebschGordanCoefficients[maxJValue + 1];
  int MinJ = JValue - this->LzMax;
  if (MinJ < 0)
    MinJ = 0;
  for (int j = maxJValue; j >= MinJ; j -= 4)
    ClebshArray[j] = ClebschGordanCoefficients(j, this->LzMax);
  for (int i = 0; i < nbrIndexSets; ++i)
    {
      double Tmp = 0.0;
      int Sum = ((indices[0] + indices[1]) << 1)  - (2 * this->LzMax);
      int TmpMinJ = MinJ;
      if (TmpMinJ < abs(Sum))
        TmpMinJ = abs(Sum);
      if (spinIndex == 0)
        {
          for (int j = maxJValue; j >= TmpMinJ; j -= 4)
            {
              Tmp += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) *
                      ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue));
            }
          Sum = ((indices[1] + indices[2]) << 1)  - (2 * this->LzMax);
          TmpMinJ = MinJ;
          if (TmpMinJ < abs(Sum))
            TmpMinJ = abs(Sum);
          for (int j = maxJValue; j >= TmpMinJ; j -= 4)
            {
              Tmp += (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) *
                      ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue));
            }
          Sum = ((indices[2] + indices[0]) << 1)  - (2 * this->LzMax);
          TmpMinJ = MinJ;
          if (TmpMinJ < abs(Sum))
            TmpMinJ = abs(Sum);
          for (int j = maxJValue; j >= TmpMinJ; j -= 4)
            {
              Tmp += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) *
                      ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue));
            }
        }
      else
        if (spinIndex == 1)
          {
            if (abs(Sum + ((indices[2] << 1) - this->LzMax)) <= JValue)
              for (int j = maxJValue; j >= TmpMinJ; j -= 4)
                {
                  Tmp += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) *
                          ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue));
                  Tmp -= (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) *
                          ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue));
                }
          }
        else
          if (spinIndex == 2)
            {
              Sum = ((indices[1] + indices[2]) << 1)  - (2 * this->LzMax);
              TmpMinJ = MinJ;
              if (TmpMinJ < abs(Sum))
                TmpMinJ = abs(Sum);
              if (abs(Sum + ((indices[0] << 1) - this->LzMax)) <= JValue)
                for (int j = maxJValue; j >= TmpMinJ; j -= 4)
                  {
                    Tmp += (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) *
                            ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue));
                  }
              Sum = ((indices[0] + indices[2]) << 1)  - (2 * this->LzMax);
              TmpMinJ = MinJ;
              if (TmpMinJ < abs(Sum))
                TmpMinJ = abs(Sum);
              if (abs(Sum + ((indices[1] << 1) - this->LzMax)) <= JValue)
                for (int j = maxJValue; j >= TmpMinJ; j -= 4)
                  {
                    Tmp -= (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) *
                            ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue));
                  }
            }
          else
            {
              Sum = ((indices[2] + indices[0]) << 1)  - (2 * this->LzMax);
              TmpMinJ = MinJ;
              if (TmpMinJ < abs(Sum))
                TmpMinJ = abs(Sum);
              if (abs(Sum + ((indices[1] << 1) - this->LzMax)) <= JValue)
                for (int j = maxJValue; j >= TmpMinJ; j -= 4)
                  {
                    Tmp += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) *
                            ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue));
                  }
              Sum = ((indices[2] + indices[1]) << 1)  - (2 * this->LzMax);
              TmpMinJ = MinJ;
              if (TmpMinJ < abs(Sum))
                TmpMinJ = abs(Sum);
              if (abs(Sum + ((indices[0] << 1) - this->LzMax)) <= JValue)
                for (int j = maxJValue; j >= TmpMinJ; j -= 4)
                  {
                    Tmp -= (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) *
                            ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue));
                  }
            }
      TmpCoefficients[i] = Tmp;
      indices += 3;
    }
  delete[] ClebshArray;

  return TmpCoefficients;
}
