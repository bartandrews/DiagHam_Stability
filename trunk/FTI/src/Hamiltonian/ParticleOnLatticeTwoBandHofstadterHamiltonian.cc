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
#include "Hamiltonian/ParticleOnLatticeTwoBandHofstadterHamiltonian.h"
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

ParticleOnLatticeTwoBandHofstadterHamiltonian::ParticleOnLatticeTwoBandHofstadterHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellsX = number of sites in the x direction
// nbrCellsY = number of sites in the y direction
// bandIndex = index of lower band n to consider (as well as band n+1)
// uPotential = strength of the repulsive two body neareast neighbor interaction
// vPotential = strength of the repulsive two body second nearest neighbor interactio
// tightBindingModel = pointer to the tight binding model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
ParticleOnLatticeTwoBandHofstadterHamiltonian::ParticleOnLatticeTwoBandHofstadterHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrCellsX, int nbrCellsY, int bandIndex1, double uPotential, double vPotential,  Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrCellsX;
  this->NbrSiteY = nbrCellsY;
  this->LzMax = NbrSiteX * NbrSiteY - 1;
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->LowerBandIndex = bandIndex1;
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
      cout  << "fast = " ;
      PrintMemorySize(cout, TmpMemory)<<endl;
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeTwoBandHofstadterHamiltonian::~ParticleOnLatticeTwoBandHofstadterHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeTwoBandHofstadterHamiltonian::EvaluateInteractionFactors()
{

  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;

  int NbrSubLattices = TightBindingModel->GetNbrBands();

  if (NbrSubLattices<2)
    {
      cout << "Hofstadter model has insufficient number of bands"<<endl;
      exit(1);
    }
  if (LowerBandIndex<0)
    LowerBandIndex=0;
  if (LowerBandIndex>NbrSubLattices-2)
    LowerBandIndex=NbrSubLattices-2;

  cout << "LowerBandIndex="<<LowerBandIndex<<endl;
    
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
	    // issue of overall Factor 1/2... to match old conventions, take sp energy twice the actual value
	    this->OneBodyInteractionFactorsupup[Index] = /* 0.5 * */ this->TightBindingModel->GetEnergy(LowerBandIndex, Index);  // mapping up <-> lower band
	    this->OneBodyInteractionFactorsdowndown[Index] = /* 0.5 * */ this->TightBindingModel->GetEnergy(LowerBandIndex+1, Index); // mapping down <-> upper band
	  }
	OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
      }

  if (this->FlatBand == false)
    for (int i=0; i<this->TightBindingModel->GetNbrStatePerBand(); ++i)
      {
	cout << "[" << this->OneBodyInteractionFactorsupup[i] << ", "<< this->OneBodyInteractionFactorsdowndown[i] << "]" << endl;
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

  // build array of InteractionFactors to simplify the creation of matrix elements:

  int NbrInternalIndices = 2;
  Complex ****** InteractionFactorsSigma = new Complex***** [NbrInternalIndices];
  for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
    {
      InteractionFactorsSigma[sigma3] = new Complex****  [NbrInternalIndices];
      for (int sigma4 = sigma3; sigma4 < NbrInternalIndices; ++sigma4)
	{
	  InteractionFactorsSigma[sigma3][sigma4] = new Complex***[NbrInternalIndices];
	  for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	    {
	      InteractionFactorsSigma[sigma3][sigma4][sigma1] = new Complex**[NbrInternalIndices];
	    }
	}
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
      
      double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if (this->FlatBand == false)
	FactorU *= this->UPotential;
      
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

		  cout << "Attention: Fermionic interactions need full implementation in ParticleOnLatticeHofstadterSingleBandHamiltonian!"<<endl;
		  exit(1);


                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 1, 1, 0) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);		  
		  
                  this->InteractionFactorsupupupup[i][Index] = -2.0 * FactorU * Tmp;
		  
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

                  this->InteractionFactorsupupdowndown[i][Index] = -2.0 * FactorU * Tmp;
		  
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
		  
                  this->InteractionFactorsdowndowndowndown[i][Index] = -2.0 * FactorU * Tmp;
		  
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

                  this->InteractionFactorsdowndownupup[i][Index] = -2.0 * FactorU * Tmp;
		  
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
		    
                  this->InteractionFactorsupdownupup[i][Index] = -2.0 * FactorU * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      //  updown downdown coefficient

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

                  this->InteractionFactorsupupupdown[i][Index] = -2.0 * FactorU * Tmp;
                  
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
		  
                  this->InteractionFactorsdowndownupdown[i][Index] = -2.0 * FactorU * Tmp;
                  
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
		  
                  this->InteractionFactorsupdownupdown[i][Index] = -2.0 * FactorU * Tmp;
                  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 <= Index2)
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
		if (Index1 <= Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2);
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }

      double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if (this->FlatBand == false)
	FactorU *= this->UPotential;
      //double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      
      Complex* TmpInteractionFactor;
      int* TmpIndices;
      int* TmpIndices2;
      for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	{
	  for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
	    {
	      InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] = new Complex*[this->NbrIntraSectorSums];
	      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		{
		  InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j] = new Complex [this->NbrIntraSectorIndicesPerSum[j] * this->NbrIntraSectorIndicesPerSum[j]];
		  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		  TmpIndices = this->IntraSectorIndicesPerSum[j];
		  for (int i1 = 0; i1 < Lim; i1 += 2)
		    {
		      TmpInteractionFactor = &(InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
		      int Index1 = TmpIndices[i1];
		      int Index2 = TmpIndices[i1 + 1];
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Complex Tmp = 0.0;
			  int Index3 = TmpIndices[i2];
			  int Index4 = TmpIndices[i2 + 1];
			  for (int s = 0; s < NbrSubLattices; ++s)
			    {
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, sigma3, sigma3, sigma1, sigma1, s, s, s, s);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, sigma3, sigma3, sigma1, sigma1, s, s, s, s);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, sigma3, sigma3, sigma1, sigma1, s, s, s, s);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, sigma3, sigma3, sigma1, sigma1, s, s, s, s);
			    }
			  Tmp *= FactorU;
			  if (Index1 == Index2)
			    Tmp *= 0.5;
			  if (Index3 == Index4)
			    Tmp *= 0.5;
			  (*TmpInteractionFactor) = 2.0*Tmp;
  
			  if (this->VPotential != 0.0)		    
			    {
			      cout<< "VPotential yet needs to be implemented!"<<endl;
			      exit(1);
			    }
			  
			  ++TmpInteractionFactor;
			}
		    }
		}
	    }
	  for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
	    {
	      for (int sigma4 = sigma3 + 1; sigma4 < NbrInternalIndices; ++sigma4)
		{
		  InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] = new Complex*[this->NbrIntraSectorSums];
		  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		    {
		      InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j] = new Complex [this->NbrIntraSectorIndicesPerSum[j] * this->NbrInterSectorIndicesPerSum[j]];
		      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		      TmpIndices = this->IntraSectorIndicesPerSum[j];
		      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
		      TmpIndices2 = this->InterSectorIndicesPerSum[j];
		      for (int i1 = 0; i1 < Lim; i1 += 2)
			{
			  TmpInteractionFactor = &(InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim2) >> 2]);
			  int Index1 = TmpIndices[i1];
			  int Index2 = TmpIndices[i1 + 1];
			  for (int i2 = 0; i2 < Lim2; i2 += 2)
			    {
			      Complex Tmp = 0.0;
			      int Index3 = TmpIndices2[i2];
			      int Index4 = TmpIndices2[i2 + 1];
			      for (int s = 0; s < NbrSubLattices; ++s)
				{
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, sigma3, sigma4, sigma1, sigma1, s, s, s, s);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, sigma4, sigma3, sigma1, sigma1, s, s, s, s);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, sigma3, sigma4, sigma1, sigma1, s, s, s, s);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, sigma4, sigma3, sigma1, sigma1, s, s, s, s);
				}
			      if (Index1 == Index2)
				Tmp *= 0.5;
			      Tmp *= FactorU;
			      (*TmpInteractionFactor) = 2.0*Tmp;
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	    }			  
	}
      for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	{
	  for (int sigma2 = sigma1 + 1; sigma2 < NbrInternalIndices; ++sigma2)
	    {
	      for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
		{
		  InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
		  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		    {
		      InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j] = new Complex [this->NbrInterSectorIndicesPerSum[j] * this->NbrIntraSectorIndicesPerSum[j]];
		      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		      TmpIndices = this->IntraSectorIndicesPerSum[j];
		      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
		      TmpIndices2 = this->InterSectorIndicesPerSum[j];
		      for (int i1 = 0; i1 < Lim2; i1 += 2)
			{
			  TmpInteractionFactor = &(InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim) >> 2]);
			  int Index1 = TmpIndices2[i1];
			  int Index2 = TmpIndices2[i1 + 1];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Complex Tmp = 0.0;
			      int Index3 = TmpIndices[i2];
			      int Index4 = TmpIndices[i2 + 1];
			      for (int s = 0; s < NbrSubLattices; ++s)
				{
 				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, sigma3, sigma3, sigma1, sigma2, s, s, s, s);
 				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, sigma3, sigma3, sigma1, sigma2, s, s, s, s);
 				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, sigma3, sigma3, sigma2, sigma1, s, s, s, s);
 				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, sigma3, sigma3, sigma2, sigma1, s, s, s, s);
				}
			      if (Index3 == Index4)
				Tmp *= 0.5;
			      Tmp *= FactorU;
			      (*TmpInteractionFactor) = 2.0*Tmp;
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	      for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
		{
		  for (int sigma4 = sigma3 + 1; sigma4 < NbrInternalIndices; ++sigma4)
		    {
		      InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
		      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
			{
			  InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j] = new Complex [this->NbrInterSectorIndicesPerSum[j] * this->NbrInterSectorIndicesPerSum[j]];
			  int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
			  TmpIndices2 = this->InterSectorIndicesPerSum[j];
			  for (int i1 = 0; i1 < Lim2; i1 += 2)
			    {
			      TmpInteractionFactor = &(InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
			      int Index1 = TmpIndices2[i1];
			      int Index2 = TmpIndices2[i1 + 1];
			      for (int i2 = 0; i2 < Lim2; i2 += 2)
				{
				  Complex Tmp = 0.0;
				  int Index3 = TmpIndices2[i2];
				  int Index4 = TmpIndices2[i2 + 1];
				  for (int s = 0; s < NbrSubLattices; ++s)
				    {
 				      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, sigma3, sigma4, sigma1, sigma2, s, s, s, s);
 				      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, sigma4, sigma3, sigma1, sigma2, s, s, s, s);
 				      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, sigma3, sigma4, sigma2, sigma1, s, s, s, s);
 				      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, sigma4, sigma3, sigma2, sigma1, s, s, s, s);
				    }
				  Tmp *= FactorU;
				  (*TmpInteractionFactor) = 2.0*Tmp;
				  ++TmpInteractionFactor;
				}
			    }
			}
		    }
		}			  
	    }
	}

      //  upup upup coefficient
      this->InteractionFactorsupupupup = InteractionFactorsSigma[0][0][0][0];
      //  upup downdown coefficient
      this->InteractionFactorsupupdowndown = InteractionFactorsSigma[0][0][1][1];
      //  downdown downdown coefficient
      this->InteractionFactorsdowndowndowndown = InteractionFactorsSigma[1][1][1][1];
      //  downdown upup coefficient
      this->InteractionFactorsdowndownupup = InteractionFactorsSigma[1][1][0][0];

      //  updown upup coefficient
      this->InteractionFactorsupdownupup = InteractionFactorsSigma[0][1][0][0];
      //  updown downdown coefficient
      this->InteractionFactorsupdowndowndown = InteractionFactorsSigma[0][1][1][1];

      //  upup updown coefficient
      this->InteractionFactorsupupupdown = InteractionFactorsSigma[0][0][0][1];
      //  downdown updown coefficient
      this->InteractionFactorsdowndownupdown = InteractionFactorsSigma[1][1][0][1];

      //  updown updown coefficient      
      this->InteractionFactorsupdownupdown = InteractionFactorsSigma[0][1][0][1];

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

Complex ParticleOnLatticeTwoBandHofstadterHamiltonian::ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2)
{
  //double dx=((double)(kx1-kx2));
  //double dy=((double)(ky1-ky2));
  //Complex Tmp = Phase(-this->KxFactor * dx/sqrt(3.0));
  Complex Tmp = 1.0;
  return Tmp;
}
