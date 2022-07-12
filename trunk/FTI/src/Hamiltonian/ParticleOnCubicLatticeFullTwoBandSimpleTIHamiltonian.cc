////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 3d topological insulator based on the simple TI        //
//              model and restricted to two bands (alternate version )        //
//                                                                            //
//                        last modification : 21/02/2013                      //
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
#include "Hamiltonian/ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// uPotential = repulsive on-site potential strength between different orbitals
// vPotential = repulsive on-site potential strength between opposite spins
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
													   int nbrSiteY, int nbrSiteZ, double uPotential, double vPotential, 
													   Abstract3DTightBindingModel* tightBindingModel, bool flatBandFlag, 
													   AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteYZ = this->NbrSiteY * this->NbrSiteZ;
  this->LzMax = nbrSiteX * nbrSiteY * nbrSiteZ - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->KzFactor = 2.0 * M_PI / ((double) this->NbrSiteZ);
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;

  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->AUpAUpPotential = 0.0;
  this->ADownADownPotential = 0.0;
  this->AUpADownPotential = 0.0;
  this->BUpBUpPotential = 0.0;
  this->BDownBDownPotential = 0.0;
  this->BUpBDownPotential = 0.0;
  this->AUpBUpPotential = 0.0;
  this->ADownBDownPotential = 0.0;
  this->AUpBDownPotential = 0.0;
  this->ADownBUpPotential = 0.0;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->AUpADownPotential = this->VPotential;
      this->BUpBDownPotential = this->VPotential;
      this->AUpBUpPotential = this->UPotential;
      this->ADownBDownPotential = this->UPotential;
      this->AUpBDownPotential = this->UPotential;
      this->ADownBUpPotential = this->UPotential;
    }
  else
    {
      this->AUpAUpPotential = this->UPotential;
      this->ADownADownPotential = this->UPotential;
      this->AUpADownPotential = this->VPotential;
      this->BUpBUpPotential = this->UPotential;
      this->BDownBDownPotential = this->UPotential;
      this->BUpBDownPotential = this->VPotential;
      this->AUpBUpPotential = this->UPotential;
      this->ADownBDownPotential = this->UPotential;
      this->AUpBDownPotential = this->UPotential;
      this->ADownBUpPotential = this->UPotential;
    }

  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
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
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::~ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->TightBindingModel->GetNbrStatePerBand()];
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;

  this->InteractionFactorsupupupup = 0;
  this->InteractionFactorsupupdowndown = 0;
  this->InteractionFactorsdowndownupup = 0;
  this->InteractionFactorsdowndowndowndown = 0;
  this->InteractionFactorsupdownupup = 0;
  this->InteractionFactorsupdowndowndown = 0;
  this->InteractionFactorsupupupdown = 0;
  this->InteractionFactorsdowndownupdown = 0;
  this->InteractionFactorsupdownupdown = 0;

  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactorsupup = new double [this->TightBindingModel->GetNbrStatePerBand()];
      this->OneBodyInteractionFactorsdowndown = new double [this->TightBindingModel->GetNbrStatePerBand()];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      for (int kz = 0; kz < this->NbrSiteZ; ++kz)
	{
	  int Index = ((Abstract3DTightBindingModel*) this->TightBindingModel)->GetLinearizedMomentumIndex(kx, ky, kz);
	  if (this->FlatBand == false)
	    {
	      this->OneBodyInteractionFactorsupup[Index] = this->TightBindingModel->GetEnergy(0, Index);
	      this->OneBodyInteractionFactorsdowndown[Index] = this->TightBindingModel->GetEnergy(1, Index);
	    }
	  OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
	}
  int NbrInternalIndices = 2;
  this->InteractionFactorsSigma = new Complex***** [NbrInternalIndices];
  for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
    {
      this->InteractionFactorsSigma[sigma3] = new Complex****  [NbrInternalIndices];
      for (int sigma4 = sigma3; sigma4 < NbrInternalIndices; ++sigma4)
	{
	  this->InteractionFactorsSigma[sigma3][sigma4] = new Complex***[NbrInternalIndices];
	  for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	    {
	      this->InteractionFactorsSigma[sigma3][sigma4][sigma1] = new Complex**[NbrInternalIndices];
	    }
	}
    }

 
  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
  this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    this->NbrIntraSectorIndicesPerSum[i] = 0;      

  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)      
	  for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
	    for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
	      ++this->NbrInterSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)) * this->NbrSiteZ + ((kz1 + kz2) % this->NbrSiteZ)];    
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      if (this->NbrInterSectorIndicesPerSum[i] > 0)
	{
	  this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
	  this->NbrInterSectorIndicesPerSum[i] = 0;
	}
    }
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)    
	  for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
	    for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
	      {
		int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)) * this->NbrSiteZ + ((kz1 + kz2) % this->NbrSiteZ);
		this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = ((kx1 * this->NbrSiteY) + ky1) * this->NbrSiteZ + kz1;
		this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = ((kx2 * this->NbrSiteY) + ky2) * this->NbrSiteZ + kz2;
		++this->NbrInterSectorIndicesPerSum[TmpSum];    
	      }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
		for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
		  {
		    int Index1 = ((kx1 * this->NbrSiteY) + ky1) * this->NbrSiteZ + kz1;
		    int Index2 = ((kx2 * this->NbrSiteY) + ky2) * this->NbrSiteZ + kz2;
		    if (Index1 < Index2)
		      ++this->NbrIntraSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)) * this->NbrSiteZ + ((kz1 + kz2) % this->NbrSiteZ)];    
		  }
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
		for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
		  {
		    int Index1 = ((kx1 * this->NbrSiteY) + ky1) * this->NbrSiteZ + kz1;
		    int Index2 = ((kx2 * this->NbrSiteY) + ky2) * this->NbrSiteZ + kz2;
		    if (Index1 < Index2)
		      {
			int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)) * this->NbrSiteZ + ((kz1 + kz2) % this->NbrSiteZ);
			this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
			this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
			++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		      }
		  }
      
      double Factor = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ));
      double FactorAUpADown = Factor * this->AUpADownPotential;
      double FactorBUpBDown = Factor * this->BUpBDownPotential;
      double FactorAUpBUp = Factor * this->AUpBUpPotential;
      double FactorADownBDown = Factor * this->ADownBDownPotential;
      double FactorAUpBDown = Factor * this->AUpBDownPotential;
      double FactorADownBUp = Factor * this->ADownBUpPotential;

      Complex Tmp;

      //  upup upup coefficient
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteYZ;
	      int ky1 = Index1 % this->NbrSiteYZ;
	      int kz1 = ky1 % this->NbrSiteZ;
	      ky1 /= this->NbrSiteZ;
	      int kx2 = Index2 / this->NbrSiteYZ;
	      int ky2 = Index2 % this->NbrSiteYZ;
	      int kz2 = ky2 % this->NbrSiteZ;
	      ky2 /= this->NbrSiteZ;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteYZ;
		  int ky3 = Index3 % this->NbrSiteYZ;
		  int kz3 = ky3 % this->NbrSiteZ;
		  ky3 /= this->NbrSiteZ;
		  int kx4 = Index4 / this->NbrSiteYZ;
		  int ky4 = Index4 % this->NbrSiteYZ;
		  int kz4 = ky4 % this->NbrSiteZ;
		  ky4 /= this->NbrSiteZ;
                  this->InteractionFactorsupupupup[i][Index] = 0.0;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupupup[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupupup[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupupup[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupupup[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupupup[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupupup[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  upup downdown coefficient
      this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteYZ;
	      int ky1 = Index1 % this->NbrSiteYZ;
	      int kz1 = ky1 % this->NbrSiteZ;
	      ky1 /= this->NbrSiteZ;
	      int kx2 = Index2 / this->NbrSiteYZ;
	      int ky2 = Index2 % this->NbrSiteYZ;
	      int kz2 = ky2 % this->NbrSiteZ;
	      ky2 /= this->NbrSiteZ;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteYZ;
		  int ky3 = Index3 % this->NbrSiteYZ;
		  int kz3 = ky3 % this->NbrSiteZ;
		  ky3 /= this->NbrSiteZ;
		  int kx4 = Index4 / this->NbrSiteYZ;
		  int ky4 = Index4 % this->NbrSiteYZ;
		  int kz4 = ky4 % this->NbrSiteZ;
		  ky4 /= this->NbrSiteZ;
                  this->InteractionFactorsupupdowndown[i][Index] = 0.0;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  downdown downdown coefficient
      this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndowndowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteYZ;
	      int ky1 = Index1 % this->NbrSiteYZ;
	      int kz1 = ky1 % this->NbrSiteZ;
	      ky1 /= this->NbrSiteZ;
	      int kx2 = Index2 / this->NbrSiteYZ;
	      int ky2 = Index2 % this->NbrSiteYZ;
	      int kz2 = ky2 % this->NbrSiteZ;
	      ky2 /= this->NbrSiteZ;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteYZ;
		  int ky3 = Index3 % this->NbrSiteYZ;
		  int kz3 = ky3 % this->NbrSiteZ;
		  ky3 /= this->NbrSiteZ;
		  int kx4 = Index4 / this->NbrSiteYZ;
		  int ky4 = Index4 % this->NbrSiteYZ;
		  int kz4 = ky4 % this->NbrSiteZ;
		  ky4 /= this->NbrSiteZ;
                  this->InteractionFactorsdowndowndowndown[i][Index] = 0.0;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  downdown upup coefficient
      this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndownupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteYZ;
	      int ky1 = Index1 % this->NbrSiteYZ;
	      int kz1 = ky1 % this->NbrSiteZ;
	      ky1 /= this->NbrSiteZ;
	      int kx2 = Index2 / this->NbrSiteYZ;
	      int ky2 = Index2 % this->NbrSiteYZ;
	      int kz2 = ky2 % this->NbrSiteZ;
	      ky2 /= this->NbrSiteZ;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteYZ;
		  int ky3 = Index3 % this->NbrSiteYZ;
		  int kz3 = ky3 % this->NbrSiteZ;
		  ky3 /= this->NbrSiteZ;
		  int kx4 = Index4 / this->NbrSiteYZ;
		  int ky4 = Index4 % this->NbrSiteYZ;
		  int kz4 = ky4 % this->NbrSiteZ;
		  ky4 /= this->NbrSiteZ;
                  this->InteractionFactorsdowndownupup[i][Index] = 0.0;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  updown upup coefficient
      this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupup[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteYZ;
	      int ky3 = Index3 % this->NbrSiteYZ;
	      int kz3 = ky3 % this->NbrSiteZ;
	      ky3 /= this->NbrSiteZ;
	      int kx4 = Index4 / this->NbrSiteYZ;
	      int ky4 = Index4 % this->NbrSiteYZ;
	      int kz4 = ky4 % this->NbrSiteZ;
	      ky4 /= this->NbrSiteZ;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteYZ;
		  int ky1 = Index1 % this->NbrSiteYZ;
		  int kz1 = ky1 % this->NbrSiteZ;
		  ky1 /= this->NbrSiteZ;
		  int kx2 = Index2 / this->NbrSiteYZ;
		  int ky2 = Index2 % this->NbrSiteYZ;
		  int kz2 = ky2 % this->NbrSiteZ;
		  ky2 /= this->NbrSiteZ;
                  this->InteractionFactorsupdownupup[i][Index] = 0.0;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  this->InteractionFactorsupdownupup[i][Index] = Conj(this->InteractionFactorsupdownupup[i][Index]);

		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      //  updown downdown coefficient
      this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupdowndowndown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteYZ;
	      int ky3 = Index3 % this->NbrSiteYZ;
	      int kz3 = ky3 % this->NbrSiteZ;
	      ky3 /= this->NbrSiteZ;
	      int kx4 = Index4 / this->NbrSiteYZ;
	      int ky4 = Index4 % this->NbrSiteYZ;
	      int kz4 = ky4 % this->NbrSiteZ;
	      ky4 /= this->NbrSiteZ;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteYZ;
		  int ky1 = Index1 % this->NbrSiteYZ;
		  int kz1 = ky1 % this->NbrSiteZ;
		  ky1 /= this->NbrSiteZ;
		  int kx2 = Index2 / this->NbrSiteYZ;
		  int ky2 = Index2 % this->NbrSiteYZ;
		  int kz2 = ky2 % this->NbrSiteZ;
		  ky2 /= this->NbrSiteZ;
                  this->InteractionFactorsupdowndowndown[i][Index] = 0.0;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
                  this->InteractionFactorsupdowndowndown[i][Index] = Conj(this->InteractionFactorsupdowndowndown[i][Index]);

		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  upup updown coefficient
      this->InteractionFactorsupupupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupupupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteYZ;
	      int ky3 = Index3 % this->NbrSiteYZ;
	      int kz3 = ky3 % this->NbrSiteZ;
	      ky3 /= this->NbrSiteZ;
	      int kx4 = Index4 / this->NbrSiteYZ;
	      int ky4 = Index4 % this->NbrSiteYZ;
	      int kz4 = ky4 % this->NbrSiteZ;
	      ky4 /= this->NbrSiteZ;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteYZ;
		  int ky1 = Index1 % this->NbrSiteYZ;
		  int kz1 = ky1 % this->NbrSiteZ;
		  ky1 /= this->NbrSiteZ;
		  int kx2 = Index2 / this->NbrSiteYZ;
		  int ky2 = Index2 % this->NbrSiteYZ;
		  int kz2 = ky2 % this->NbrSiteZ;
		  ky2 /= this->NbrSiteZ;
                  this->InteractionFactorsupupupdown[i][Index] = 0.0;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
                  this->InteractionFactorsupupupdown[i][Index] = Conj(this->InteractionFactorsupupupdown[i][Index]);

		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  downdown updown coefficient
      this->InteractionFactorsdowndownupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsdowndownupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteYZ;
	      int ky3 = Index3 % this->NbrSiteYZ;
	      int kz3 = ky3 % this->NbrSiteZ;
	      ky3 /= this->NbrSiteZ;
	      int kx4 = Index4 / this->NbrSiteYZ;
	      int ky4 = Index4 % this->NbrSiteYZ;
	      int kz4 = ky4 % this->NbrSiteZ;
	      ky4 /= this->NbrSiteZ;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteYZ;
		  int ky1 = Index1 % this->NbrSiteYZ;
		  int kz1 = ky1 % this->NbrSiteZ;
		  ky1 /= this->NbrSiteZ;
		  int kx2 = Index2 / this->NbrSiteYZ;
		  int ky2 = Index2 % this->NbrSiteYZ;
		  int kz2 = ky2 % this->NbrSiteZ;
		  ky2 /= this->NbrSiteZ;
                  this->InteractionFactorsdowndownupdown[i][Index] = 0.0;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  this->InteractionFactorsdowndownupdown[i][Index] = Conj(this->InteractionFactorsdowndownupdown[i][Index]);

		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}


      //  updown updown coefficient
      this->InteractionFactorsupdownupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteYZ;
	      int ky1 = Index1 % this->NbrSiteYZ;
	      int kz1 = ky1 % this->NbrSiteZ;
	      ky1 /= this->NbrSiteZ;
	      int kx2 = Index2 / this->NbrSiteYZ;
	      int ky2 = Index2 % this->NbrSiteYZ;
	      int kz2 = ky2 % this->NbrSiteZ;
	      ky2 /= this->NbrSiteZ;
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteYZ;
		  int ky3 = Index3 % this->NbrSiteYZ;
		  int kz3 = ky3 % this->NbrSiteZ;
		  ky3 /= this->NbrSiteZ;
		  int kx4 = Index4 / this->NbrSiteYZ;
		  int ky4 = Index4 % this->NbrSiteYZ;
		  int kz4 = ky4 % this->NbrSiteZ;
		  ky4 /= this->NbrSiteZ;
                  this->InteractionFactorsupdownupdown[i][Index] = 0.0;
//		  cout << Index1 << " " << Index2 << " " << Index3 << " " << Index4 << " : "; 
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorAUpADown * Tmp;
//		  cout << this->InteractionFactorsupdownupdown[i][Index] << " " << (-2.0 * FactorAUpADown * Tmp) << " | ";
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
//		  cout << this->InteractionFactorsupdownupdown[i][Index] << " " << (-2.0 * FactorBUpBDown * Tmp) << " | ";
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
//		  cout << this->InteractionFactorsupdownupdown[i][Index] << " " << (-2.0 * FactorAUpBUp * Tmp) << " | ";
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorADownBDown * Tmp;
//		  cout << this->InteractionFactorsupdownupdown[i][Index] << " " << (-2.0 * FactorADownBDown * Tmp) << " | ";
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
//		  cout << this->InteractionFactorsupdownupdown[i][Index] << " " << (-2.0 * FactorAUpBDown * Tmp) << " | ";
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorADownBUp * Tmp;
//		  cout << this->InteractionFactorsupdownupdown[i][Index] << " " << (-2.0 * FactorADownBUp * Tmp) << " " << endl;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      // bosonic interaction
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
		for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
		  {
		    int Index1 = ((kx1 * this->NbrSiteY) + ky1) * this->NbrSiteZ + kz1;
		    int Index2 = ((kx2 * this->NbrSiteY) + ky2) * this->NbrSiteZ + kz2;
		    if (Index1 <= Index2)
		      {
			int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
											 (ky1 + ky2) % this->NbrSiteY,
											 (kz1 + kz2) % this->NbrSiteZ);
			++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		      }
		  }
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
		for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
		  {
		    int Index1 = ((kx1 * this->NbrSiteY) + ky1) * this->NbrSiteZ + kz1;
		    int Index2 = ((kx2 * this->NbrSiteY) + ky2) * this->NbrSiteZ + kz2;
		    if (Index1 <= Index2)
		      {
			int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
											 (ky1 + ky2) % this->NbrSiteY,
											 (kz1 + kz2) % this->NbrSiteZ);
			this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
			this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
			++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		      }
		  }
      
      double Factor = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ));
      double FactorAUpADown = Factor * this->AUpADownPotential;
      double FactorBUpBDown = Factor * this->BUpBDownPotential;
      double FactorAUpAUp = Factor * this->AUpAUpPotential;
      double FactorADownADown = Factor * this->ADownADownPotential;
      double FactorBUpBUp = Factor * this->BUpBUpPotential;
      double FactorBDownBDown = Factor * this->BDownBDownPotential;
      double FactorAUpBUp = Factor * this->AUpBUpPotential;
      double FactorADownBDown = Factor * this->ADownBDownPotential;
      double FactorAUpBDown = Factor * this->AUpBDownPotential;
      double FactorADownBUp = Factor * this->ADownBUpPotential;

      Complex Tmp;
      Complex* TmpInteractionFactor;

      int* TmpIndices;
      int* TmpIndices2;
      for (int sigma1 = 0; sigma1 < 2; ++sigma1)
	{
	  for (int sigma3 = 0; sigma3 < 2; ++sigma3)
	    {
	      this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] = new Complex*[this->NbrIntraSectorSums];
	      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		{
		  this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j] = new Complex [this->NbrIntraSectorIndicesPerSum[j] * this->NbrIntraSectorIndicesPerSum[j]];
		  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		  TmpIndices = this->IntraSectorIndicesPerSum[j];
		  for (int i1 = 0; i1 < Lim; i1 += 2)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
		      int Index1 = TmpIndices[i1];
		      int Index2 = TmpIndices[i1 + 1];
		      int kx1, ky1, kz1;
		      this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1, kz1);
		      int kx2, ky2, kz2;
		      this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2, kz2);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  int Index3 = TmpIndices[i2];
			  int Index4 = TmpIndices[i2 + 1];
			  int kx3, ky3, kz3;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3, kz3);
			  int kx4, ky4, kz4;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4, kz4);

			  (*TmpInteractionFactor) = 0.0;
			  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			  (*TmpInteractionFactor) += 2.0 * FactorAUpAUp * Tmp;
			  
			  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			  (*TmpInteractionFactor) += 2.0 * FactorADownADown * Tmp;
			  
			  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			  (*TmpInteractionFactor) += 2.0 * FactorAUpADown * Tmp;
			  
			  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			  (*TmpInteractionFactor) += 2.0 * FactorBUpBUp * Tmp;
			  
			  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			  (*TmpInteractionFactor) += 2.0 * FactorBDownBDown * Tmp;
			  
			  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			  (*TmpInteractionFactor) += 2.0 * FactorBUpBDown * Tmp;
			  
			  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			  (*TmpInteractionFactor) += 2.0 * FactorAUpBUp * Tmp;
			  
			  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			  (*TmpInteractionFactor) += 2.0 * FactorADownBDown * Tmp;
			  
			  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			  (*TmpInteractionFactor) += 2.0 * FactorAUpBDown * Tmp;
			  
			  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma1, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			  (*TmpInteractionFactor) += 2.0 * FactorADownBUp * Tmp;
			  
			  if (Index1 == Index2)
			    (*TmpInteractionFactor) *= 0.5;
			  if (Index3 == Index4)
			    (*TmpInteractionFactor) *= 0.5;
			  
			  ++TmpInteractionFactor;
			}
		    }
		}
	    }

	  for (int sigma3 = 0; sigma3 < 2; ++sigma3)
	    {
	      for (int sigma4 = sigma3 + 1; sigma4 < 2; ++sigma4)
		{
		  this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] = new Complex*[this->NbrIntraSectorSums];
		  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		    {
		      this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j] = new Complex [this->NbrIntraSectorIndicesPerSum[j] * this->NbrInterSectorIndicesPerSum[j]];
		      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		      TmpIndices = this->IntraSectorIndicesPerSum[j];
		      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
		      TmpIndices2 = this->InterSectorIndicesPerSum[j];
		      for (int i1 = 0; i1 < Lim; i1 += 2)
			{
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim2) >> 2]);
			  int Index1 = TmpIndices[i1];
			  int Index2 = TmpIndices[i1 + 1];
			  int kx1, ky1, kz1;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1, kz1);
			  int kx2, ky2, kz2;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2, kz2);
			  for (int i2 = 0; i2 < Lim2; i2 += 2)
			    {
			      int Index3 = TmpIndices2[i2];
			      int Index4 = TmpIndices2[i2 + 1];
			      int kx3, ky3, kz3;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3, kz3);
			      int kx4, ky4, kz4;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4, kz4);

			      (*TmpInteractionFactor) = 0.0;
			      
			      Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma1, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma1, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorAUpAUp * Tmp;
			      
			      Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma1, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma1, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorADownADown * Tmp;
			      
			      Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma1, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma1, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorAUpADown * Tmp;
			      
			      Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma1, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma1, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorBUpBUp * Tmp;
			      
			      Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma1, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma1, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorBDownBDown * Tmp;
			      
			      Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma1, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma1, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorBUpBDown * Tmp;
			      
			      Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma1, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma1, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorAUpBUp * Tmp;
			      
			      Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma1, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma1, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorADownBDown * Tmp;
			      
			      Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma1, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma1, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorAUpBDown * Tmp;
			      
			      Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma1, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma1, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorADownBUp * Tmp;
			      
			      if (Index3 == Index4)
				(*TmpInteractionFactor) *= 0.5;

			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	    }			  
	}

      for (int sigma1 = 0; sigma1 < 2; ++sigma1)
	{
	  for (int sigma2 = sigma1 + 1; sigma2 < 2; ++sigma2)
	    {
	      for (int sigma3 = 0; sigma3 < 2; ++sigma3)
		{
		  this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
		  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		    {
		      this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j] = new Complex [this->NbrInterSectorIndicesPerSum[j] * this->NbrIntraSectorIndicesPerSum[j]];
		      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		      TmpIndices = this->IntraSectorIndicesPerSum[j];
		      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
		      TmpIndices2 = this->InterSectorIndicesPerSum[j];
		      for (int i1 = 0; i1 < Lim2; i1 += 2)
			{
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim) >> 2]);
			  int Index1 = TmpIndices2[i1];
			  int Index2 = TmpIndices2[i1 + 1];
			  int kx1, ky1, kz1;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1, kz1);
			  int kx2, ky2, kz2;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2, kz2);
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      int Index3 = TmpIndices[i2];
			      int Index4 = TmpIndices[i2 + 1];
			      int kx3, ky3, kz3;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3, kz3);
			      int kx4, ky4, kz4;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4, kz4);

			      (*TmpInteractionFactor) = 0.0;
			      
			      Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorAUpAUp * Tmp;
			      
			      Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorADownADown * Tmp;
			      
			      Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorAUpADown * Tmp;
			      
			      Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorBUpBUp * Tmp;
			      
			      Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorBDownBDown * Tmp;
			      
			      Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorBUpBDown * Tmp;
			      
			      Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorAUpBUp * Tmp;
			      
			      Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorADownBDown * Tmp;
			      
			      Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorAUpBDown * Tmp;
			      
			      Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma3, sigma3, sigma1, sigma2, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma3, sigma3, sigma2, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
			      (*TmpInteractionFactor) += 2.0 * FactorADownBUp * Tmp;
			      
			      if (Index1 == Index2)
				(*TmpInteractionFactor) *= 0.5;
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	      for (int sigma3 = 0; sigma3 < 2; ++sigma3)
		{
		  for (int sigma4 = sigma3 + 1; sigma4 < 2; ++sigma4)
		    {
		      this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
		      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
			{
			  this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j] = new Complex [this->NbrInterSectorIndicesPerSum[j] * this->NbrInterSectorIndicesPerSum[j]];
			  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
			  TmpIndices = this->IntraSectorIndicesPerSum[j];
			  int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
			  TmpIndices2 = this->InterSectorIndicesPerSum[j];
			  for (int i1 = 0; i1 < Lim2; i1 += 2)
			    {
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
			      int Index1 = TmpIndices2[i1];
			      int Index2 = TmpIndices2[i1 + 1];
			      int kx1, ky1, kz1;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1, kz1);
			      int kx2, ky2, kz2;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2, kz2);
			      for (int i2 = 0; i2 < Lim2; i2 += 2)
				{
				  int Index3 = TmpIndices2[i2];
				  int Index4 = TmpIndices2[i2 + 1];
				  int kx3, ky3, kz3;
				  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3, kz3);
				  int kx4, ky4, kz4;
				  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4, kz4);

				  (*TmpInteractionFactor) = 0.0;

				  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma2, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma2, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma2, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma2, sigma1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
				  (*TmpInteractionFactor) += 2.0 * FactorAUpAUp * Tmp;
				  
				  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma2, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma2, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma2, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma2, sigma1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementADownADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
				  (*TmpInteractionFactor) += 2.0 * FactorADownADown * Tmp;
				  
				  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma2, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma2, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma2, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma2, sigma1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
				  (*TmpInteractionFactor) += 2.0 * FactorAUpADown * Tmp;
				  
				  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma2, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma2, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma2, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma2, sigma1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
				  (*TmpInteractionFactor) += 2.0 * FactorBUpBUp * Tmp;
				  
				  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma2, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma2, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma2, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma2, sigma1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementBDownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
				  (*TmpInteractionFactor) += 2.0 * FactorBDownBDown * Tmp;
				  
				  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma2, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma2, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma2, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma2, sigma1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
				  (*TmpInteractionFactor) += 2.0 * FactorBUpBDown * Tmp;
				  
				  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma2, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma2, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma2, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma2, sigma1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
				  (*TmpInteractionFactor) += 2.0 * FactorAUpBUp * Tmp;
				  
				  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma2, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma2, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma2, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma2, sigma1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
				  (*TmpInteractionFactor) += 2.0 * FactorADownBDown * Tmp;
				  
				  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma2, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma2, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma2, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma2, sigma1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
				  (*TmpInteractionFactor) += 2.0 * FactorAUpBDown * Tmp;
				  
				  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma3, sigma4, sigma1, sigma2, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, sigma3, sigma4, sigma2, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, sigma4, sigma3, sigma1, sigma2, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, sigma4, sigma3, sigma2, sigma1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
				  (*TmpInteractionFactor) += 2.0 * FactorADownBUp * Tmp;

				  ++TmpInteractionFactor;
				}
			    }
			}
		    }
		}			  
	    }
	}
    }

  delete[] OneBodyBasis;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the two body interaction between two sites A with up spins
//
// kx1 = momentum along x for the creation operator on first A site with spin up
// ky1 = momentum along y for the creation operator on first A site with spin up
// kz1 = momentum along z for the creation operator on first A site with spin up
// kx2 = momentum along x for the creation operator on second A site with spin up
// ky2 = momentum along y for the creation operator on second A site with spin up
// kz2 = momentum along z for the creation operator on second A site with spin up
// kx3 = momentum along x for the annihilation operator on first A site with spin up
// ky3 = momentum along y for the annihilation operator on first A site with spin up
// kz3 = momentum along z for the annihilation operator on first A site with spin up
// kx4 = momentum along x for the annihilation operator on second A site with spin up
// ky4 = momentum along y for the annihilation operator on second A site with spin up
// kz4 = momentum along z for the annihilation operator on second A site with spin up
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementAUpAUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  Complex Tmp = 1.0 ;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A with down spins
//
// kx1 = momentum along x for the creation operator on first A site with spin down
// ky1 = momentum along y for the creation operator on first A site with spin down
// kz1 = momentum along z for the creation operator on first A site with spin down
// kx2 = momentum along x for the creation operator on second A site with spin down
// ky2 = momentum along y for the creation operator on second A site with spin down
// kz2 = momentum along z for the creation operator on second A site with spin down
// kx3 = momentum along x for the annihilation operator on first A site with spin down
// ky3 = momentum along y for the annihilation operator on first A site with spin down
// kz3 = momentum along z for the annihilation operator on first A site with spin down
// kx4 = momentum along x for the annihilation operator on second A site with spin down
// ky4 = momentum along y for the annihilation operator on second A site with spin down
// kz4 = momentum along z for the annihilation operator on second A site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementADownADown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
}

// compute the matrix element for the two body interaction between two sites B with up spins
//
// kx1 = momentum along x for the creation operator on first B site with spin up
// ky1 = momentum along y for the creation operator on first B site with spin up
// kz1 = momentum along z for the creation operator on first B site with spin up
// kx2 = momentum along x for the creation operator on second B site with spin up
// ky2 = momentum along y for the creation operator on second B site with spin up
// kz2 = momentum along z for the creation operator on second B site with spin up
// kx3 = momentum along x for the annihilation operator on first B site with spin up
// ky3 = momentum along y for the annihilation operator on first B site with spin up
// kz3 = momentum along z for the annihilation operator on first B site with spin up
// kx4 = momentum along x for the annihilation operator on second B site with spin up
// ky4 = momentum along y for the annihilation operator on second B site with spin up
// kz4 = momentum along z for the annihilation operator on second B site with spin up
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementBUpBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  Complex Tmp = 1.0 ;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B with down spins
//
// kx1 = momentum along x for the creation operator on first B site with spin down
// ky1 = momentum along y for the creation operator on first B site with spin down
// kz1 = momentum along z for the creation operator on first B site with spin down
// kx2 = momentum along x for the creation operator on second B site with spin down
// ky2 = momentum along y for the creation operator on second B site with spin down
// kz2 = momentum along z for the creation operator on second B site with spin down
// kx3 = momentum along x for the annihilation operator on first B site with spin down
// ky3 = momentum along y for the annihilation operator on first B site with spin down
// kz3 = momentum along z for the annihilation operator on first B site with spin down
// kx4 = momentum along x for the annihilation operator on second B site with spin down
// ky4 = momentum along y for the annihilation operator on second B site with spin down
// kz4 = momentum along z for the annihilation operator on second B site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementBDownBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return this->ComputeTwoBodyMatrixElementBUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
}

// compute the matrix element for the two body interaction between two sites A and B with up spins
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kz2 = momentum along z for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// kz4 = momentum along z for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  Complex Tmp = 1.0 ;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and B with down spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kz1 = momentum along z for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kz2 = momentum along z for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kz3 = momentum along z for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// kz4 = momentum along z for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
}

// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kz1 = momentum along z for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kz2 = momentum along z for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kz3 = momentum along z for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// kz4 = momentum along z for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementADownBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
}

// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kz2 = momentum along z for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// kz4 = momentum along z for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
}

// compute the matrix element for the two body interaction between two sites A with opposite spins 
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on A site with spin down
// ky2 = momentum along y for the creation operator on A site with spin down
// kz2 = momentum along z for the creation operator on A site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on A site with spin down
// ky4 = momentum along y for the annihilation operator on A site with spin down
// kz4 = momentum along z for the annihilation operator on A site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  Complex Tmp = 1.0;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B with opposite spins 
//
// kx1 = momentum along x for the creation operator on B site with spin up
// ky1 = momentum along y for the creation operator on B site with spin up
// kz1 = momentum along z for the creation operator on B site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kz2 = momentum along z for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on B site with spin up
// ky3 = momentum along y for the annihilation operator on B site with spin up
// kz3 = momentum along z for the annihilation operator on B site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// kz4 = momentum along z for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFullTwoBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  Complex Tmp = 1.0;
  return Tmp;
}

