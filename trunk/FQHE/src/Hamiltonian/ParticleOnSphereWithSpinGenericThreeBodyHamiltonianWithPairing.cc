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
#include "Hamiltonian/ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"

  
#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//

ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing::ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// threeBodyPseudoPotential32 = array with the three-body pseudo-potentials in the S=3/2 sector
// maxRelativeAngularMomentum32 =  maximum relative angular momentum that is used in ThreeBodyPseudoPotential in the S=3/2 sector
// threeBodyPseudoPotential12 = array with the three-body pseudo-potentials in the S=1/2 sector
// maxRelativeAngularMomentum12 =  maximum relative angular momentum that is used in ThreeBodyPseudoPotential in the S=1/2 sector
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing::ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, 
													 double* threeBodyPseudoPotential32, int maxRelativeAngularMomentum32,
															       double* threeBodyPseudoPotential12, int maxRelativeAngularMomentum12, double alpha, double **pseudoPotential,
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
  if (pseudoPotential!=NULL)
    this->FullTwoBodyFlag = true;
  else
    this->FullTwoBodyFlag = false;
  this->Alpha = alpha;
  this->PseudoPotentials = new double* [4];
  for (int j = 0; j < 4; ++j)
    {
      this->PseudoPotentials[j] = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	if (pseudoPotential!=NULL)
	  this->PseudoPotentials[j][i] = pseudoPotential[j][this->LzMax - i];
	else
	  this->PseudoPotentials[j][i] = 0.0;
    }
  
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
 
  this->MaxRelativeAngularMomentum32 = maxRelativeAngularMomentum32;
  this->NbrThreeBodyPseudoPotentials32 = this->MaxRelativeAngularMomentum32 + 1;
  if (this->NbrThreeBodyPseudoPotentials32 == 0)
    {
      this->MaxRelativeAngularMomentum32 = 0;
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	{
	  this->NbrThreeBodyPseudoPotentials32 = 4;
	  this->MaxRelativeAngularMomentum32 = 3;
	  this->ThreeBodyPseudoPotentials32 = new double[4];
	  this->ThreeBodyPseudoPotentials32[3] = 0.0;
	}
      else
	{
	  this->NbrThreeBodyPseudoPotentials32 = 1;
	  this->MaxRelativeAngularMomentum32 = 0;
	  this->ThreeBodyPseudoPotentials32 = new double[1];
	  this->ThreeBodyPseudoPotentials32[0] = 0.0;
	}
    }
  else
    {
      this->ThreeBodyPseudoPotentials32 = new double[this->NbrThreeBodyPseudoPotentials32];
      for (int k = 0; k < this->NbrThreeBodyPseudoPotentials32; ++k)
	this->ThreeBodyPseudoPotentials32[k] = threeBodyPseudoPotential32[k];
    }
  this->MaxRelativeAngularMomentum12 = maxRelativeAngularMomentum12;
  this->NbrThreeBodyPseudoPotentials12 = this->MaxRelativeAngularMomentum12 + 1;
  if (this->NbrThreeBodyPseudoPotentials12 > 0)
    {
      this->ThreeBodyPseudoPotentials12 = new double[this->NbrThreeBodyPseudoPotentials12];
      for (int k = 0; k < this->NbrThreeBodyPseudoPotentials12; ++k)
	this->ThreeBodyPseudoPotentials12[k] = threeBodyPseudoPotential12[k];
    }
  else
    {
      this->ThreeBodyPseudoPotentials12 = 0;
    }

  this->Architecture = architecture;

  this->EvaluateInteractionFactors();
  this->EvaluatePairingInteractionFactors();
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
// threeBodyPseudoPotential32 = array with the three-body pseudo-potentials in the S=3/2 sector
// maxRelativeAngularMomentum32 =  maximum relative angular momentum that is used in ThreeBodyPseudoPotential in the S=3/2 sector
// threeBodyPseudoPotential12 = array with the three-body pseudo-potentials in the S=1/2 sector
// maxRelativeAngularMomentum12 =  maximum relative angular momentum that is used in ThreeBodyPseudoPotential in the S=1/2 sector
// pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing::ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, 
													 double* threeBodyPseudoPotential32, int maxRelativeAngularMomentum32,
															       double* threeBodyPseudoPotential12, int maxRelativeAngularMomentum12, double alpha,
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
  this->Alpha = alpha;

  this->PseudoPotentials = new double* [4];
  for (int j = 0; j < 4; ++j)
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
 
 
  this->MaxRelativeAngularMomentum32 = maxRelativeAngularMomentum32;
  this->NbrThreeBodyPseudoPotentials32 = this->MaxRelativeAngularMomentum32 + 1;
  if (this->NbrThreeBodyPseudoPotentials32 == 0)
    {
      this->MaxRelativeAngularMomentum32 = 0;
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	{
	  this->NbrThreeBodyPseudoPotentials32 = 4;
	  this->MaxRelativeAngularMomentum32 = 3;
	  this->ThreeBodyPseudoPotentials32 = new double[4];
	  this->ThreeBodyPseudoPotentials32[3] = 0.0;
	}
      else
	{
	  this->NbrThreeBodyPseudoPotentials32 = 1;
	  this->MaxRelativeAngularMomentum32 = 0;
	  this->ThreeBodyPseudoPotentials32 = new double[1];
	  this->ThreeBodyPseudoPotentials32[0] = 0.0;
	}
    }
  else
    {
      this->ThreeBodyPseudoPotentials32 = new double[this->NbrThreeBodyPseudoPotentials32];
      for (int k = 0; k < this->NbrThreeBodyPseudoPotentials32; ++k)
	this->ThreeBodyPseudoPotentials32[k] = threeBodyPseudoPotential32[k];
    }
  this->MaxRelativeAngularMomentum12 = maxRelativeAngularMomentum12;
  this->NbrThreeBodyPseudoPotentials12 = this->MaxRelativeAngularMomentum12 + 1;
  if (this->NbrThreeBodyPseudoPotentials12 > 0)
    {
      this->ThreeBodyPseudoPotentials12 = new double[this->NbrThreeBodyPseudoPotentials12];
      for (int k = 0; k < this->NbrThreeBodyPseudoPotentials12; ++k)
	this->ThreeBodyPseudoPotentials12[k] = threeBodyPseudoPotential12[k];
    }
  else
    {
      this->ThreeBodyPseudoPotentials12 = 0;
    }

  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->EvaluatePairingInteractionFactors();
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

ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing::~ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing()
{
  delete[] this->ThreeBodyPseudoPotentials32;
  if (this->ThreeBodyPseudoPotentials12 != 0)
    delete[] this->ThreeBodyPseudoPotentials12;
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

AbstractHamiltonian* ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing::Clone ()
{
  return 0;
}


// evaluate all interaction factors to do with pairing
//   

void ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing::EvaluatePairingInteractionFactors()
{
  ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
  double ClebschCoef;
  long TotalNbrInteractionFactors = 0;
  int J;
  
  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  double TmpCoefficient = 0.0;


  //FERMIONS  *** this code is not really used for lattices, as we consider bosons
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      cout << "Fermions not defined in ThreeBodyHamiltonianWithPairing"<<endl;
      exit(1);
    }
  else // bosons
    {
      
      //************************************************************************************************************ /
      //********************************     P A I R I N G    T E R M S ******************************************** /
      //************************************************************************************************************ /
      // Pairing term indices have the same structure as intra terms
	      
      this->NbrPairingSectorSums = 2 * this->LzMax + 1;
      this->NbrPairingSectorIndicesPerSum = new int[this->NbrPairingSectorSums];
      for (int i = 0; i < this->NbrPairingSectorSums; ++i)
	this->NbrPairingSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  ++this->NbrPairingSectorIndicesPerSum[m1 + m2];
      this->PairingSectorIndicesPerSum = new int* [this->NbrPairingSectorSums];
      for (int i = 0; i < this->NbrPairingSectorSums; ++i)
	{
	  this->PairingSectorIndicesPerSum[i] = new int[2 * this->NbrPairingSectorIndicesPerSum[i]];      
	  this->NbrPairingSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  {
	    this->PairingSectorIndicesPerSum[m1 + m2][this->NbrPairingSectorIndicesPerSum[m1 + m2] << 1] = m1;
	    this->PairingSectorIndicesPerSum[m1 + m2][1 + (this->NbrPairingSectorIndicesPerSum[m1 + m2] << 1)] = m2;
	    ++this->NbrPairingSectorIndicesPerSum[m1 + m2];
	  }

      //Fill in the matrix elements

      this->InteractionFactorsPairing= new double* [this->NbrPairingSectorSums];
      for (int i = 0; i < this->NbrPairingSectorSums; ++i)
	{
	  this->InteractionFactorsPairing[i] = new double[this->NbrPairingSectorIndicesPerSum[i] * this->NbrPairingSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrPairingSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->PairingSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
	      int m2 = (this->PairingSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
	      for (int j2 = 0; j2 < this->NbrPairingSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->PairingSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		  int m4 = (this->PairingSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
		  Clebsch.InitializeCoefficientIterator(m1, m2);
		  this->InteractionFactorsPairing[i][Index] = 0.0;
		  while (Clebsch.Iterate(J, ClebschCoef))
		    {
		      if (((J >> 1) & 1) != Sign)
			{
			  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
			  this->InteractionFactorsPairing[i][Index] += this->PseudoPotentials[3][J >> 1] * TmpCoefficient;
 			  if (J==2*LzMax) // anomalous terms arising from effective lattice model in V0 channel
 			    this->InteractionFactorsPairing[i][Index] += -2.0*M_PI*this->Alpha * TmpCoefficient;
			}
		    }
		  if (m1 != m2)
		    this->InteractionFactorsPairing[i][Index] *= 2.0;
		  if (m3 != m4)
		    this->InteractionFactorsPairing[i][Index] *= 2.0;
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
        }
    }

  cout << "nbr pairing interactions = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}


// evaluate all interaction factors
//   

void ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing::EvaluateInteractionFactors()
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
	      if (this->MaxRelativeAngularMomentum32 <= TmpMaxRelativeMomentum)
		TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum32;
	      int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2];
	      while (((3 * this->LzMax) - TmpMaxRelativeMomentum)  < TmpSum)
		--TmpMaxRelativeMomentum;
	      double** TmpProjectorCoefficients = new double* [TmpMaxRelativeMomentum + 1];
	      if  ((this->NbrThreeBodyPseudoPotentials32 > 0) && (this->ThreeBodyPseudoPotentials32[3] != 0.0))
		TmpProjectorCoefficients[3] = this->ComputeProjectorCoefficients(6, 1, TmpNIndices2, Lim, (2 * this->LzMax) - 2, 0);
	      for (int i = 5; i <= TmpMaxRelativeMomentum; ++i)  
		if ((this->NbrThreeBodyPseudoPotentials32 > 0) && (this->ThreeBodyPseudoPotentials32[i] != 0.0))
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
		      if (this->ThreeBodyPseudoPotentials32[3] != 0.0)
			{
			  TmpInteraction2UpUp += 4.0 * this->ThreeBodyPseudoPotentials32[3] * TmpProjectorCoefficients[3][i] * TmpProjectorCoefficients[3][j];
			  TmpInteraction2DownDown += 4.0 * this->ThreeBodyPseudoPotentials32[3] * TmpProjectorCoefficients[3][i] * TmpProjectorCoefficients[3][j];
			}
		      for (int k = 5; k <= TmpMaxRelativeMomentum; ++k)  
			{
			  if (this->ThreeBodyPseudoPotentials32[k] != 0.0)
			    {
			      TmpInteraction2UpUp += 4.0 * this->ThreeBodyPseudoPotentials32[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
			      TmpInteraction2DownDown += 4.0 * this->ThreeBodyPseudoPotentials32[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
			    }
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
	      if (this->ThreeBodyPseudoPotentials32[3] != 0.0)
		delete[] TmpProjectorCoefficients[3];
	      for (int i = 5; i <= TmpMaxRelativeMomentum; ++i)  
		if (this->ThreeBodyPseudoPotentials32[i] != 0.0)
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
	      if (this->MaxRelativeAngularMomentum32 <= TmpMaxRelativeMomentum)
		TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum32;
	      int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2];
	      while (((3 * this->LzMax) - TmpMaxRelativeMomentum)  < TmpSum)
		--TmpMaxRelativeMomentum;
	      double** TmpProjectorCoefficients = new double* [TmpMaxRelativeMomentum + 1];
	      for (int i = 1; i <= TmpMaxRelativeMomentum; ++i)  
		if (this->ThreeBodyPseudoPotentials32[i] != 0.0)
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
			  if (this->ThreeBodyPseudoPotentials32[k] != 0.0)
			    {
			      TmpInteraction2UpUpDown += this->ThreeBodyPseudoPotentials32[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
			      TmpInteraction2DownDownUp += this->ThreeBodyPseudoPotentials32[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
			    }
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
		if ((this->ThreeBodyPseudoPotentials32[i] != 0.0) || (this->ThreeBodyPseudoPotentials32[i] != 0.0))
		  delete[] TmpProjectorCoefficients[i];


	      TmpNIndices2 = this->SortedIndicesPerSum[3][1][MinSum];
	      TmpNbrNIndices = TmpNbrNIndices3;
	      for (int i = 1; i <= TmpMaxRelativeMomentum; ++i)  
		if (this->ThreeBodyPseudoPotentials32[i] != 0.0)
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
			  if (this->ThreeBodyPseudoPotentials32[k] != 0.0)
			    {
			      TmpInteraction2UpUpDown += this->ThreeBodyPseudoPotentials32[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
			      TmpInteraction2DownDownUp += this->ThreeBodyPseudoPotentials32[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
			    }
			}
		    }
		  ++TmpNbrNIndices;
		}
	      for (int i = 0; i <= TmpMaxRelativeMomentum; ++i)  
		if (this->ThreeBodyPseudoPotentials32[i] != 0.0)
		  delete[] TmpProjectorCoefficients[i];


	      TmpNIndices2 = this->SortedIndicesPerSum[3][1][MinSum];
	      TmpNbrNIndices = TmpNbrNIndices3;
	      for (int i = 1; i <= TmpMaxRelativeMomentum; ++i)  
		if (this->ThreeBodyPseudoPotentials32[i] != 0.0)
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
			  if (this->ThreeBodyPseudoPotentials32[k] != 0.0)
			    {
			      TmpInteraction2UpUpDown += this->ThreeBodyPseudoPotentials32[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
			      TmpInteraction2DownDownUp += this->ThreeBodyPseudoPotentials32[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
			    }
			}
		    }
		  ++TmpNbrNIndices;
		}
	      for (int i = 0; i <= TmpMaxRelativeMomentum; ++i)  
		if (this->ThreeBodyPseudoPotentials32[i] != 0.0)
		  delete[] TmpProjectorCoefficients[i];


	      delete[] TmpProjectorCoefficients;		
	    }
	}
      delete[] TmpInteractionCoeffients;
    }
  else
    {
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
	  if (this->MaxRelativeAngularMomentum32 <= TmpMaxRelativeMomentum)
	    TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum32;
	  int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2];
	  while (((3 * this->LzMax) - TmpMaxRelativeMomentum)  < TmpSum)
	    --TmpMaxRelativeMomentum;
	  double** TmpProjectorCoefficients = new double* [TmpMaxRelativeMomentum + 1];
	  if (this->ThreeBodyPseudoPotentials32[0] != 0.0)
	    TmpProjectorCoefficients[0] = this->ComputeProjectorCoefficients(0, 1, TmpNIndices2, Lim, 2 * this->LzMax, 0);
	  for (int i = 2; i <= TmpMaxRelativeMomentum; ++i)  
	    if (this->ThreeBodyPseudoPotentials32[i] != 0.0)
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
		  double& TmpInteraction2DownDown = TmpInteractionDownDownDown[j];
		  TmpInteraction2UpUp = 0.0;
		  TmpInteraction2DownDown = 0.0;
		  if (this->ThreeBodyPseudoPotentials32[0] != 0.0)
		    {
		      TmpInteraction2UpUp += this->ThreeBodyPseudoPotentials32[0] * TmpProjectorCoefficients[0][i] * TmpProjectorCoefficients[0][j];
		    TmpInteraction2DownDown += this->ThreeBodyPseudoPotentials32[0] * TmpProjectorCoefficients[0][i] * TmpProjectorCoefficients[0][j];
		    }
		  for (int k = 2; k <= TmpMaxRelativeMomentum; ++k)
		    if (this->ThreeBodyPseudoPotentials32[k] != 0.0)
		      {
			TmpInteraction2UpUp += this->ThreeBodyPseudoPotentials32[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
			TmpInteraction2DownDown += this->ThreeBodyPseudoPotentials32[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
		      }
		  TmpInteraction2UpUp *= 2.0 * TmpSymmetryFactors[i] * TmpSymmetryFactors[j] / 3.0;
		  TmpInteraction2DownDown *= 2.0 * TmpSymmetryFactors[i] * TmpSymmetryFactors[j] / 3.0;
		}
	      for (int j = 0; j < 3; ++j)
		{
		  (*TmpNIndices) = (*TmpNIndices2);
		  ++TmpNIndices;
		  ++TmpNIndices2;
		}
	      ++TmpNbrNIndices;
	    }
	  if (this->ThreeBodyPseudoPotentials32[0] != 0.0)
	    delete[] TmpProjectorCoefficients[0];
	  for (int i = 2; i <= TmpMaxRelativeMomentum; ++i)  
	    if (this->ThreeBodyPseudoPotentials32[i] != 0.0)
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

      GetAllIndices(this->NbrLzValue, 3, this->NbrSortedIndicesPerSum[3][1], this->SortedIndicesPerSum[3][1]);

      this->MaxSumIndices[3][1] = this->LzMax * 3;
      this->MinSumIndices[3][1] = 0;
      TmpInteractionCoeffients = new double[this->MaxSumIndices[3][1] + 1];
      TmpInteractionCoeffients[0] = 1.0;
      TmpInteractionCoeffients[1] = 1.0;
      for (int i = 0; i <= this->MaxSumIndices[3][1]; ++i)
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
	      int TmpMaxRelativeMomentum = 10;
	      if (this->MaxRelativeAngularMomentum32 <= TmpMaxRelativeMomentum)
		TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum32;
	      int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2];
	      while (((3 * this->LzMax) - TmpMaxRelativeMomentum)  < TmpSum)
		--TmpMaxRelativeMomentum;
	      double** TmpProjectorCoefficients = new double* [TmpMaxRelativeMomentum + 1];
	      for (int i = 0; i <= TmpMaxRelativeMomentum; ++i)  
		if (this->ThreeBodyPseudoPotentials32[i] != 0.0)
		  TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim, 2 * this->LzMax, 1);
	      
	      
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
			  if (this->ThreeBodyPseudoPotentials32[k] != 0.0)
			    {
			      TmpInteraction2UpUpDown += this->ThreeBodyPseudoPotentials32[k] / 2.0 * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
			      TmpInteraction2DownDownUp += this->ThreeBodyPseudoPotentials32[k] / 2.0 * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
			    }
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
		if (this->ThreeBodyPseudoPotentials32[i] != 0.0)
		  delete[] TmpProjectorCoefficients[i];
	      delete[] TmpProjectorCoefficients;

	      // spin 1/2 channel
	      
	      TmpMaxRelativeMomentum = 10;
	      if (this->MaxRelativeAngularMomentum12 <= TmpMaxRelativeMomentum)
		TmpMaxRelativeMomentum = this->MaxRelativeAngularMomentum12;
	      if (TmpMaxRelativeMomentum >= 0)
		{
		  TmpNIndices2 = this->SortedIndicesPerSum[3][1][MinSum];
		  TmpNbrNIndices = TmpNbrNIndices3;
		  TmpProjectorCoefficients = new double* [TmpMaxRelativeMomentum + 1];
		  for (int i = 0; i <= TmpMaxRelativeMomentum; ++i)  
		    if (this->ThreeBodyPseudoPotentials12[i] != 0.0)
		      TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficientsSpin12E112(2 * i, 1, TmpNIndices2, Lim, 2 * this->LzMax, 1);
		  for (int i = 0; i < Lim; ++i)
		    {
		      double* TmpInteractionUpUpDown = this->MNNBodyInteractionFactors[3][2][TmpNbrNIndices];
		      double* TmpInteractionDownDownUp = this->MNNBodyInteractionFactors[3][3][TmpNbrNIndices];
		      for (int j = 0; j < Lim; ++j)
			{
			  double& TmpInteraction2UpUpDown = TmpInteractionUpUpDown[j];
			  double& TmpInteraction2DownDownUp = TmpInteractionDownDownUp[j];
			  for (int k = 1; (k == 1) && (k <= TmpMaxRelativeMomentum); ++k)  
			    {
			      if (this->ThreeBodyPseudoPotentials12[k] != 0.0)
				{
				  TmpInteraction2UpUpDown += this->ThreeBodyPseudoPotentials12[k] / 12.0 * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
				  TmpInteraction2DownDownUp += this->ThreeBodyPseudoPotentials12[k] / 12.0 * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
				}
			    }
			}
		      ++TmpNbrNIndices;
		    }
		  for (int i = 0; i <= TmpMaxRelativeMomentum; ++i)  
		    if (this->ThreeBodyPseudoPotentials12[i] != 0.0)
		      delete[] TmpProjectorCoefficients[i];
		  delete[] TmpProjectorCoefficients;		  		  
		}
	    }
	}
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
	{ 
	  this->NbrIntraSectorSums = 2 * this->LzMax+1;
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
		      if (m1 != m2)
			{
			  this->InteractionFactorsupup[i][Index] *= 2.0;
			  this->InteractionFactorsdowndown[i][Index] *= 2.0;
			}
		      if (m3 != m4)
			{
			  this->InteractionFactorsupup[i][Index] *= 2.0;
			  this->InteractionFactorsdowndown[i][Index] *= 2.0;
			}
		      this->InteractionFactorsupup[i][Index] *= 0.5;
		      this->InteractionFactorsdowndown[i][Index] *= 0.5;
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
		  double Factor = 1.0;
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

// compute all projector coefficient associated to a given relative angular momentum between 3 particles in the spin 3/2 sector
//
// relativeMomentum = value of twice the relative angular momentum between the 3 particles
// degeneracyIndex = optional degeneracy index for relative angular momentum greater than 5 for bosons (8 for fermions)
// indices = array that contains all possible sets of indices (size of the array is 3 * nbrIndexSets)
// nbrIndexSets = number of sets
// maxJValue = twice the maximum total angular momentum two particles can have
// spinIndex = indicate for which of three body operators coeeficients are computed (0 for up-up-up, 1 for up-up-down, 2 for up-down-up, 3 for down-up-up) 

double* ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing::ComputeProjectorCoefficients(int relativeMomentum, int degeneracyIndex, int* indices, int nbrIndexSets, int maxJValue, int spinIndex)
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
	  if (abs(Sum + ((indices[2] << 1) - this->LzMax)) <= JValue)
	    {
	      for (int j = maxJValue; j >= TmpMinJ; j -= 4)
		{
		  Tmp += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) *
			  ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue)); 
		}
	    }
	  Sum = ((indices[1] + indices[2]) << 1)  - (2 * this->LzMax);
	  TmpMinJ = MinJ;
	  if (TmpMinJ < abs(Sum))
	    TmpMinJ = abs(Sum);
	  if (abs(Sum + ((indices[0] << 1) - this->LzMax)) <= JValue)
	    {
	      for (int j = maxJValue; j >= TmpMinJ; j -= 4)
		{
		  Tmp += (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) * 
			  ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue)); 
		}
	    }
	  Sum = ((indices[2] + indices[0]) << 1)  - (2 * this->LzMax);
	  TmpMinJ = MinJ;
	  if (TmpMinJ < abs(Sum))
	    TmpMinJ = abs(Sum);
	  if (abs(Sum + ((indices[1] << 1) - this->LzMax)) <= JValue)
	    {
	      for (int j = maxJValue; j >= TmpMinJ; j -= 4)
		{
		  Tmp += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) * 
			  ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue)); 
		}
	    }
	}
      else
	{
	  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	    {
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
	    }
	  else
	    {
	      if (spinIndex == 1)
		{
		  if (abs(Sum + ((indices[2] << 1) - this->LzMax)) <= JValue)
		    for (int j = maxJValue; j >= TmpMinJ; j -= 4)
		      {
			Tmp += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) *
				ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue)); 
			Tmp += (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) *
				ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue)); 
		      }
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
			Tmp += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) * 
				ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue)); 
		      }
		  
		  Sum = ((indices[2] + indices[0]) << 1)  - (2 * this->LzMax);
		  TmpMinJ = MinJ;
		  if (TmpMinJ < abs(Sum))
		    TmpMinJ = abs(Sum);
		  if (abs(Sum + ((indices[1] << 1) - this->LzMax)) <= JValue)
		    for (int j = maxJValue; j >= TmpMinJ; j -= 4)
		      {
			Tmp += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[0] << 1) - this->LzMax), j) * 
				ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue)); 
		      }
		  Sum = ((indices[2] + indices[1]) << 1)  - (2 * this->LzMax);
		  TmpMinJ = MinJ;
		  if (TmpMinJ < abs(Sum))
		    TmpMinJ = abs(Sum);
		  if (abs(Sum + ((indices[0] << 1) - this->LzMax)) <= JValue)
		    for (int j = maxJValue; j >= TmpMinJ; j -= 4)
		      {
			Tmp += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[1] << 1) - this->LzMax), j) * 
				ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue)); 
		      }
		}
	      Tmp /= 6.0;
	    }
	}

      TmpCoefficients[i] = Tmp;
      indices += 3;
    }
  delete[] ClebshArray;

  return TmpCoefficients;
}

// compute all projector coefficient associated to a given relative angular momentum between 3 particles in the spin 1/2 sector along z_1^u + z_2^u - 2 z_3^d
//
// relativeMomentum = value of twice the relative angular momentum between the 3 particles
// degeneracyIndex = optional degeneracy index for relative angular momentum greater than 5 for bosons (8 for fermions)
// indices = array that contains all possible sets of indices (size of the array is 3 * nbrIndexSets)
// nbrIndexSets = number of sets
// maxJValue = twice the maximum total angular momentum two particles can have
// spinIndex = indicate for which of three body operators coeeficients are computed (0 for up-up-up, 1 for up-up-down, 2 for up-down-up, 3 for down-up-up) 

double* ParticleOnSphereWithSpinGenericThreeBodyHamiltonianWithPairing::ComputeProjectorCoefficientsSpin12E112(int relativeMomentum, int degeneracyIndex, int* indices, int nbrIndexSets, int maxJValue, int spinIndex)
{
  if (spinIndex == 0)
    {
      return 0;
    }
  double* TmpCoefficients = new double [nbrIndexSets];
  int JValue = (3 * this->LzMax) - relativeMomentum;

  ClebschGordanCoefficients Clebsh (this->LzMax, this->LzMax);
  ClebschGordanCoefficients* ClebshArray = new ClebschGordanCoefficients[maxJValue + 1];
  int MinJ = JValue - this->LzMax;
  if (MinJ < 0)
    MinJ = 0;
  for (int j = maxJValue; j >= MinJ; j -= 4)
    ClebshArray[j] = ClebschGordanCoefficients(j, this->LzMax);

  ClebschGordanCoefficients ClebshSpin (1, 1);
  ClebschGordanCoefficients* ClebshSpinArray = new ClebschGordanCoefficients[2];
  ClebshSpinArray[0] = ClebschGordanCoefficients(0, 1);
  ClebshSpinArray[1] = ClebschGordanCoefficients(2, 1);

  for (int i = 0; i < nbrIndexSets; ++i)
    {
      double Tmp = 0.0;
      int Sum = ((indices[0] + indices[1]) << 1)  - (2 * this->LzMax);
      int TmpMinJ = MinJ;
      if (TmpMinJ < abs(Sum))
	TmpMinJ = abs(Sum);
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	{
	  cout << "ComputeProjectorCoefficientsSpin12E112 is not implemented for fermions" << endl;
// 	  if (spinIndex == 1)
// 	    {
// 	      if (abs(Sum + ((indices[2] << 1) - this->LzMax)) <= JValue)
// 		for (int j = maxJValue; j >= TmpMinJ; j -= 4)
// 		  {
// 		    Tmp += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) *
// 			    ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue)); 
// 		    Tmp -= (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) *
// 			    ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue)); 
// 		  }
// 	    }
// 	  else
// 	    if (spinIndex == 2)
// 	      {
// 		Sum = ((indices[1] + indices[2]) << 1)  - (2 * this->LzMax);
// 		TmpMinJ = MinJ;
// 		if (TmpMinJ < abs(Sum))
// 		  TmpMinJ = abs(Sum);
// 		if (abs(Sum + ((indices[0] << 1) - this->LzMax)) <= JValue)
// 		  for (int j = maxJValue; j >= TmpMinJ; j -= 4)
// 		    {
// 		      Tmp += (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) * 
// 			      ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue)); 
// 		    }
// 		Sum = ((indices[0] + indices[2]) << 1)  - (2 * this->LzMax);
// 		TmpMinJ = MinJ;
// 		if (TmpMinJ < abs(Sum))
// 		  TmpMinJ = abs(Sum);
// 		if (abs(Sum + ((indices[1] << 1) - this->LzMax)) <= JValue)
// 		  for (int j = maxJValue; j >= TmpMinJ; j -= 4)
// 		    {
// 		      Tmp -= (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) * 
// 			      ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue)); 
// 		    }
// 	      }
// 	    else
// 	      {
// 		Sum = ((indices[2] + indices[0]) << 1)  - (2 * this->LzMax);
// 		TmpMinJ = MinJ;
// 		if (TmpMinJ < abs(Sum))
// 		  TmpMinJ = abs(Sum);
// 		if (abs(Sum + ((indices[1] << 1) - this->LzMax)) <= JValue)
// 		  for (int j = maxJValue; j >= TmpMinJ; j -= 4)
// 		    {
// 		      Tmp += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) * 
// 			      ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue)); 
// 		    }
// 		Sum = ((indices[2] + indices[1]) << 1)  - (2 * this->LzMax);
// 		TmpMinJ = MinJ;
// 		if (TmpMinJ < abs(Sum))
// 		      TmpMinJ = abs(Sum);
// 		if (abs(Sum + ((indices[0] << 1) - this->LzMax)) <= JValue)
// 		  for (int j = maxJValue; j >= TmpMinJ; j -= 4)
// 		    {
// 		      Tmp -= (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) * 
// 			      ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue)); 
// 			}
// 	      }
	}
      else
	{
	  double Tmp2 = 0.0;

	  // 1 - uud and udd
	  if (abs(Sum + ((indices[2] << 1) - this->LzMax)) <= JValue)
	    for (int j = maxJValue; j >= TmpMinJ; j -= 4)
	      {
		Tmp2 += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) *
			 ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue)); 
	      }
	  double TmpSpin = 0.0;
	  if (spinIndex == 1)
	    {
	      TmpSpin += (ClebshSpin.GetCoefficient(1, 1, 2) *
			  ClebshSpinArray[1].GetCoefficient(2, -1, 1)); 
	    }
	  else	      
	    {
	      TmpSpin += (ClebshSpin.GetCoefficient(1, -1, 0) *
			  ClebshSpinArray[0].GetCoefficient(0, -1, 1)); 
	      TmpSpin += (ClebshSpin.GetCoefficient(1, -1, 2) *
			  ClebshSpinArray[1].GetCoefficient(0, -1, 1)); 
	    }
	  Tmp += Tmp2 * TmpSpin;

	  // 2 - uud and dud
	  Tmp2 = 0.0;
	  if (abs(Sum + ((indices[2] << 1) - this->LzMax)) <= JValue)
	    for (int j = maxJValue; j >= TmpMinJ; j -= 4)
	      {
		Tmp2 += (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) *
			 ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue)); 
	      }
	  TmpSpin = 0.0;
	  if (spinIndex == 1)
	    {
	      TmpSpin += (ClebshSpin.GetCoefficient(1, 1, 2) *
			  ClebshSpinArray[1].GetCoefficient(2, -1, 1)); 
	    }
	  else	      
	    {
	      TmpSpin += (ClebshSpin.GetCoefficient(-1, 1, 0) *
			  ClebshSpinArray[0].GetCoefficient(0, -1, 1)); 
	      TmpSpin += (ClebshSpin.GetCoefficient(-1, 1, 2) *
			  ClebshSpinArray[1].GetCoefficient(0, -1, 1)); 
	    }
	  Tmp += Tmp2 * TmpSpin;

	  // 3 - udu and ddu
	  Tmp2 = 0.0;
	  Sum = ((indices[1] + indices[2]) << 1)  - (2 * this->LzMax);
	  TmpMinJ = MinJ;
	  if (TmpMinJ < abs(Sum))
	    TmpMinJ = abs(Sum);
	  if (abs(Sum + ((indices[0] << 1) - this->LzMax)) <= JValue)
	    for (int j = maxJValue; j >= TmpMinJ; j -= 4)
	      {
		Tmp2 += (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) * 
			 ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue)); 
	      }
	  TmpSpin = 0.0;
	  if (spinIndex == 1)
	    {
	      TmpSpin += (ClebshSpin.GetCoefficient(1, -1, 0) *
			  ClebshSpinArray[0].GetCoefficient(0, 1, 1)); 
	      TmpSpin += (ClebshSpin.GetCoefficient(1, -1, 2) *
			  ClebshSpinArray[1].GetCoefficient(0, 1, 1)); 
	    }
	  else	      
	    {
	      TmpSpin += (ClebshSpin.GetCoefficient(-1, -1, 2) *
			  ClebshSpinArray[1].GetCoefficient(-2, 1, 1)); 
	    }
	  Tmp += Tmp2 * TmpSpin;


	  // 4 - udu and udd
	  Tmp2 = 0.0;
	  Sum = ((indices[0] + indices[2]) << 1)  - (2 * this->LzMax);
	  TmpMinJ = MinJ;
	  if (TmpMinJ < abs(Sum))
	    TmpMinJ = abs(Sum);
	  if (abs(Sum + ((indices[1] << 1) - this->LzMax)) <= JValue)
	    for (int j = maxJValue; j >= TmpMinJ; j -= 4)
	      {
		Tmp2 += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) * 
			 ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue)); 
	      }
	  TmpSpin = 0.0;
	  if (spinIndex == 1)
	    {
	      TmpSpin += (ClebshSpin.GetCoefficient(1, -1, 0) *
			  ClebshSpinArray[0].GetCoefficient(0, 1, 1)); 
	      TmpSpin += (ClebshSpin.GetCoefficient(1, -1, 2) *
			  ClebshSpinArray[1].GetCoefficient(0, 1, 1)); 
	    }
	  else	      
	    {
	      TmpSpin += (ClebshSpin.GetCoefficient(1, -1, 0) *
			  ClebshSpinArray[0].GetCoefficient(0, -1, 1)); 
	      TmpSpin += (ClebshSpin.GetCoefficient(1, -1, 2) *
			  ClebshSpinArray[1].GetCoefficient(0, -1, 1)); 
	    }
	  Tmp += Tmp2 * TmpSpin;


	  // 5 - duu and dud	  
	  Tmp2 = 0.0;
	  Sum = ((indices[2] + indices[0]) << 1)  - (2 * this->LzMax);
	  TmpMinJ = MinJ;
	  if (TmpMinJ < abs(Sum))
	    TmpMinJ = abs(Sum);
	  if (abs(Sum + ((indices[1] << 1) - this->LzMax)) <= JValue)
	    for (int j = maxJValue; j >= TmpMinJ; j -= 4)
	      {
		Tmp2 += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[0] << 1) - this->LzMax), j) * 
			 ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue)); 
	      }
	  TmpSpin = 0.0;
	  if (spinIndex == 1)
	    {
	      TmpSpin += (ClebshSpin.GetCoefficient(-1, 1, 0) *
			  ClebshSpinArray[0].GetCoefficient(0, 1, 1)); 
	      TmpSpin += (ClebshSpin.GetCoefficient(-1, 1, 2) *
			  ClebshSpinArray[1].GetCoefficient(0, 1, 1)); 
	    }
	  else	      
	    {
	      TmpSpin += (ClebshSpin.GetCoefficient(-1, 1, 0) *
			  ClebshSpinArray[0].GetCoefficient(0, -1, 1)); 
	      TmpSpin += (ClebshSpin.GetCoefficient(-1, 1, 2) *
			  ClebshSpinArray[1].GetCoefficient(0, -1, 1)); 
	    }
	  Tmp += Tmp2 * TmpSpin;


	  // 6 - duu and ddu
	  Tmp2 = 0.0;	  
	  Sum = ((indices[2] + indices[1]) << 1)  - (2 * this->LzMax);
	  TmpMinJ = MinJ;
	  if (TmpMinJ < abs(Sum))
	    TmpMinJ = abs(Sum);
	  if (abs(Sum + ((indices[0] << 1) - this->LzMax)) <= JValue)
	    for (int j = maxJValue; j >= TmpMinJ; j -= 4)
	      {
		Tmp2 += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[1] << 1) - this->LzMax), j) * 
			 ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue)); 
	      }
	  TmpSpin = 0.0;
	  if (spinIndex == 1)
	    {
	      TmpSpin += (ClebshSpin.GetCoefficient(-1, 1, 0) *
			  ClebshSpinArray[0].GetCoefficient(0, 1, 1)); 
	      TmpSpin += (ClebshSpin.GetCoefficient(-1, 1, 2) *
			  ClebshSpinArray[1].GetCoefficient(0, 1, 1)); 
	    }
	  else	      
	    {
 	      TmpSpin += (ClebshSpin.GetCoefficient(-1, -1, 2) *
 			  ClebshSpinArray[1].GetCoefficient(-2, 1, 1)); 
	    }
	  Tmp += Tmp2 * TmpSpin;

	}
      TmpCoefficients[i] = Tmp;
      indices += 3;
    }
  delete[] ClebshArray;

  return TmpCoefficients;
}

