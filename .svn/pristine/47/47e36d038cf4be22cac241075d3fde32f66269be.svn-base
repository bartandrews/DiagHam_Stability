////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 2d topological insulator based on the simple TI        //
//                       model and restricted to two bands                    //
//                                                                            //
//                        last modification : 27/09/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeTwoBandThirringHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"
#include "GeneralTools/StringTools.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

ParticleOnLatticeTwoBandThirringHamiltonian::ParticleOnLatticeTwoBandThirringHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = repulsive potential strength between orbitals A, B in same unit cell
// tightBindingModel = pointer to tight binding structure of model
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeTwoBandThirringHamiltonian::ParticleOnLatticeTwoBandThirringHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, 
												   AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
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
      PrintMemorySize(cout, TmpMemory);
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeTwoBandThirringHamiltonian::~ParticleOnLatticeTwoBandThirringHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeTwoBandThirringHamiltonian::EvaluateInteractionFactors()
{

  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactorsupup = new double [this->TightBindingModel->GetNbrStatePerBand()];
      this->OneBodyInteractionFactorsdowndown = new double [this->TightBindingModel->GetNbrStatePerBand()];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsupup[Index] = 0.5 * this->TightBindingModel->GetEnergy(0, Index);
	    this->OneBodyInteractionFactorsdowndown[Index] = 0.5 * this->TightBindingModel->GetEnergy(1, Index);
	  }
	OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
      }
 
  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    this->NbrIntraSectorIndicesPerSum[i] = 0;      

  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)      
	  ++this->NbrInterSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2)];
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
	  {
	    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2);
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	  }
  
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 < Index2)
		  ++this->NbrIntraSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2)];
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
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2);
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      
      double FactorAB = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if (this->FlatBand == false)
	FactorAB *= this->UPotential;
      
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
	      int kx1,ky1;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
	      int kx2,ky2;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);
	      
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3,ky3;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
		  int kx4,ky4;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);		  
		  
                  this->InteractionFactorsupupupup[i][Index] = -2.0 * FactorAB * Tmp;
		  
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
	      int kx1,ky1;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
	      int kx2,ky2;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3,ky3;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
		  int kx4,ky4;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);

                  this->InteractionFactorsupupdowndown[i][Index] = -2.0 * FactorAB * Tmp;
		  
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

	      int kx1,ky1;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
	      int kx2,ky2;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3,ky3;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
		  int kx4,ky4;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
		  
		  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);
		  
                  this->InteractionFactorsdowndowndowndown[i][Index] = -2.0 * FactorAB * Tmp;
		  
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
	      int kx1,ky1;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
	      int kx2,ky2;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3,ky3;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
		  int kx4,ky4;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

		  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);

                  this->InteractionFactorsdowndownupup[i][Index] = -2.0 * FactorAB * Tmp;
		  
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
	      int kx3,ky3;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
	      int kx4,ky4;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1,ky1;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
		  int kx2,ky2;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);

		  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);
		    
                  this->InteractionFactorsupdownupup[i][Index] = -2.0 * FactorAB * Tmp;
		  
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
	      int kx3,ky3;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
	      int kx4,ky4;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1,ky1;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
		  int kx2,ky2;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);

		  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);
		  
                  this->InteractionFactorsupdowndowndown[i][Index] = -2.0 * FactorAB * Tmp;
		  
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
	      	      int kx3,ky3;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
	      int kx4,ky4;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1,ky1;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
		  int kx2,ky2;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);

		  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);

                  this->InteractionFactorsupupupdown[i][Index] = -2.0 * FactorAB * Tmp;
                  
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
	      int kx3,ky3;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
	      int kx4,ky4;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1,ky1;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
		  int kx2,ky2;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);

		  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);
		  
                  this->InteractionFactorsdowndownupdown[i][Index] = -2.0 * FactorAB * Tmp;
                  
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
	      int kx1,ky1;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
	      int kx2,ky2;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3,ky3;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
		  int kx4,ky4;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);


		  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);
		  
                  this->InteractionFactorsupdownupdown[i][Index] = -2.0 * FactorAB * Tmp;
                  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      cout << "Bosonic statistics yet to be implemented in Thirring model"<<endl;

      exit(1);
      /*
      // bosonic interaction
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 <= Index2)
		  ++this->NbrIntraSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY))];    
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
	      {
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 <= Index2)
		  {
		    int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY));
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      
      double Factor = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorAUpADown = Factor * this->VPotential;
      double FactorBUpBDown = Factor * this->VPotential;
      if ((this->FlatBand == false) || (this->VPotential != 0.0))
	Factor *= this->UPotential;
      double FactorAUpAUp = Factor;
      double FactorBUpBUp = Factor;
      double FactorADownADown = Factor;
      double FactorBDownBDown = Factor;
      double FactorAUpBUp = Factor;
      double FactorADownBDown = Factor;
      double FactorAUpBDown = Factor;
      double FactorADownBUp = Factor;

      Complex Tmp;

      //  upup upup coefficient
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupupupup[i][Index] = 0.0;		  

                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupupupup[i][Index] += 2.0 * FactorADownBUp * Tmp;
		  
		  if (Index1 == Index2)
		    this->InteractionFactorsupupupup[i][Index] *= 0.5;
		  if (Index3 == Index4)
		    this->InteractionFactorsupupupup[i][Index] *= 0.5;

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
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupupdowndown[i][Index] = 0.0;

                   Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupupdowndown[i][Index] += 2.0 * FactorAUpAUp * Tmp;

		  if (Index1 == Index2)
		    this->InteractionFactorsupupdowndown[i][Index] *= 0.5;
		  if (Index3 == Index4)
		    this->InteractionFactorsupupdowndown[i][Index] *= 0.5;

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
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsdowndowndowndown[i][Index] = 0.0;

                   Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += 2.0 * FactorAUpAUp * Tmp;
		  
		  if (Index1 == Index2)
		    this->InteractionFactorsdowndowndowndown[i][Index] *= 0.5;
		  if (Index3 == Index4)
		    this->InteractionFactorsdowndowndowndown[i][Index] *= 0.5;

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
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsdowndownupup[i][Index] = 0.0;

                   Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndownupup[i][Index] += 2.0 * FactorAUpAUp * Tmp;

		  
		  if (Index1 == Index2)
		    this->InteractionFactorsdowndownupup[i][Index] *= 0.5;
		  if (Index3 == Index4)
		    this->InteractionFactorsdowndownupup[i][Index] *= 0.5;

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
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupdownupup[i][Index] = 0.0;

                   Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdownupup[i][Index] += 2.0 * FactorAUpAUp * Tmp;


		  if (Index3 == Index4)
		    this->InteractionFactorsupdownupup[i][Index] *= 0.5;

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
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupdowndowndown[i][Index] = 0.0;

		  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdowndowndown[i][Index] += 2.0 * FactorAUpAUp * Tmp;		  

		  if (Index3 == Index4)
		    this->InteractionFactorsupdowndowndown[i][Index] *= 0.5;

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
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupupupdown[i][Index] = 0.0;

		  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupupupdown[i][Index] += 2.0 * FactorAUpAUp * Tmp;


		  if (Index1 == Index2)
		    this->InteractionFactorsupupupdown[i][Index] *= 0.5;

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
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsdowndownupdown[i][Index] = 0.0;

		  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndownupdown[i][Index] += 2.0 * FactorAUpAUp * Tmp;

		  if (Index1 == Index2)
		    this->InteractionFactorsdowndownupdown[i][Index] *= 0.5;

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
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupdownupdown[i][Index] = 0.0;

		  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpAUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdownupdown[i][Index] += 2.0 * FactorAUpAUp * Tmp;

		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      */
    }
      

  delete[] OneBodyBasis;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}


// compute the matrix element for the two body interaction between and A sites and the adjacent B site 
//
// kx1 = annihilation momentum along x for the B site
// ky1 = annihilation momentum along y for the B site
// kx2 = creation momentum along x for the B site
// ky2 = creation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeTwoBandThirringHamiltonian::ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2)
{
  //double dx=((double)(kx1-kx2));
  //double dy=((double)(ky1-ky2));
  //Complex Tmp = Phase(-this->KxFactor * dx/sqrt(3.0));
  Complex Tmp = 1.0;
  return Tmp;
}
