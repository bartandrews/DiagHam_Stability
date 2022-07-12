////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                          class author: Gunnar Möller                       //
//                                                                            //
//      class for a square Hofstadter model with interacting particles        //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 10/12/2015                      //
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
#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <cmath>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::sin;
using std::cos;

// constructor
//
ParticleOnLatticeHofstadterSingleBandHamiltonian::ParticleOnLatticeHofstadterSingleBandHamiltonian()
{
  this->BandIndex = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellX = number of sites in the x direction
// nbrCellY = number of sites in the y direction
// bandIndex = index of band to consider
// uPotential = strength of the repulsive two body neareast neighbor interaction
// vPotential = strength of the repulsive two body second nearest neighbor interaction
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeHofstadterSingleBandHamiltonian::ParticleOnLatticeHofstadterSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, int bandIndex, double uPotential, double vPotential,  Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrCellX;
  this->NbrSiteY = nbrCellY;
  this->LzMax = nbrCellX * nbrCellY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->OneBodyPotential = 0;
  this->BandIndex = bandIndex;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      cout << "fast = ";
      PrintMemorySize(cout, TmpMemory)<< endl;
      this->EnableFastMultiplication();
    }
}


// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellX = number of sites in the x direction
// nbrCellY = number of sites in the y direction
// bandIndex = index of band to consider
// uPotential = strength of the repulsive two body neareast neighbor interaction
// vPotential = strength of the repulsive two body second nearest neighbor interaction
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeHofstadterSingleBandHamiltonian::ParticleOnLatticeHofstadterSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, int bandIndex, double uPotential, double vPotential,  Abstract2DTightBindingModel* tightBindingModel, double** oneBodyPotential,  bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrCellX;
  this->NbrSiteY = nbrCellY;
  this->LzMax = nbrCellX * nbrCellY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->OneBodyPotential = oneBodyPotential;
  this->BandIndex = bandIndex;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      cout << "fast = ";
      PrintMemorySize(cout, TmpMemory)<< endl;
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeHofstadterSingleBandHamiltonian::~ParticleOnLatticeHofstadterSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeHofstadterSingleBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  int NbrSublattices = TightBindingModel->GetNbrBands();
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
  if ((this->FlatBand == false) || (this->OneBodyPotential != 0))
    {
      this->OneBodyInteractionFactors = new double [this->TightBindingModel->GetNbrStatePerBand()];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	if (this->FlatBand == false)
	  this->OneBodyInteractionFactors[Index] = this->TightBindingModel->GetEnergy(BandIndex, Index);
	if (this->OneBodyPotential != 0)
	  this->OneBodyInteractionFactors[Index] = this->OneBodyPotential[kx][ky] / ((double) (this->NbrSiteX * this->NbrSiteY));
	OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
      }

  if (this->FlatBand == false)
    for (int i=0; i<this->TightBindingModel->GetNbrStatePerBand(); ++i)
      {
	cout << "[" << this->OneBodyInteractionFactors[i] << "]" << endl;
      }

  double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
  if ((this->FlatBand == false) || (this->OneBodyPotential != 0))
    FactorU *= this->UPotential;
  double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

  if (FactorU==0.0 && FactorV==0.0 && (this->OneBodyPotential == 0))
    {
      std::cerr << "Error: HofstadterHamiltonian created with interaction zero - set non-zero --u-potential or --v-potential"<<std::endl;
      exit(1);
    }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {    
      this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;


      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 < Index2)
		  ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2)];
	      }
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  if (this->NbrSectorIndicesPerSum[i]  > 0)
	    {
	      this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];      
	      this->NbrSectorIndicesPerSum[i] = 0;
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
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }

      this->InteractionFactors = new Complex* [this->NbrSectorSums];

      Complex Tmp;

      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1,ky1;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
	      int kx2,ky2;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);
	      
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3,ky3;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
		  int kx4,ky4;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

		  this->InteractionFactors[i][Index] = 0.0;

		  if (this->UPotential != 0.0)		    
		    {
		      int xI, yI, dRx, dRy, sF;
		      Complex translationPhase;
		      int nbrNeighbors=4;
		      int dx[4]={1,-1,0,0};
		      int dy[4]={0,0,1,-1};
		      Tmp=0.0;
		      for (int s=0; s<NbrSublattices; ++s)
			{			  
			  TightBindingModel->DecodeSublatticeIndex(s, xI, yI);
			  
			  for (int n=0; n<nbrNeighbors; ++n)
			    {
			      sF = TightBindingModel->EncodeSublatticeIndex(xI + dx[n], yI + dy[n], dRx, dRy, translationPhase); // calculate final sublattice index.
			      
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx3, ky3);
			      Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx4, ky4);
			      Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx3, ky3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx4, ky4);
			    }
			}
                  	  
		      this->InteractionFactors[i][Index] += 2.0 * FactorU * Tmp;
		    }

		  if (this->VPotential != 0.0)		    
		    {
		      int xI, yI, dRx, dRy, sF;
		      Complex translationPhase;
		      int nbrNeighbors=4;
		      int dx[4]={1,1,-1,-1};
		      int dy[4]={1,-1,1,-1};
		      Tmp=0.0;
		      for (int s=0; s<NbrSublattices; ++s)
			{			  
			  TightBindingModel->DecodeSublatticeIndex(s, xI, yI);
			  
			  for (int n=0; n<nbrNeighbors; ++n)
			    {
			      sF = TightBindingModel->EncodeSublatticeIndex(xI + dx[n], yI + dy[n], dRx, dRy, translationPhase); // calculate final sublattice index.
			      
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, BandIndex, BandIndex, BandIndex, BandIndex, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx3, ky3);
			      Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, BandIndex, BandIndex, BandIndex, BandIndex, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx4, ky4);
			      Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, BandIndex, BandIndex, BandIndex, BandIndex, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx3, ky3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, BandIndex, BandIndex, BandIndex, BandIndex, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx4, ky4);
			    }
			}
                  	  
		      this->InteractionFactors[i][Index] += 2.0 * FactorV * Tmp;
		    }
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else // Bosonic Statistics
    {
      this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;


      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 <= Index2)
		  ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2)];
	      }
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  if (this->NbrSectorIndicesPerSum[i]  > 0)
	    {
	      this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];      
	      this->NbrSectorIndicesPerSum[i] = 0;
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
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }

      this->InteractionFactors = new Complex* [this->NbrSectorSums];

      Complex Tmp;

      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1,ky1;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
	      int kx2,ky2;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);
	      
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3,ky3;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
		  int kx4,ky4;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

		  Tmp=0.0;

		  for (int s=0; s<NbrSublattices; ++s)
		    {
		      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, BandIndex, BandIndex, BandIndex, BandIndex, s, s, s, s) * this->ComputeEmbeddingOnSite(s, kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, BandIndex, BandIndex, BandIndex, BandIndex, s, s, s, s) * this->ComputeEmbeddingOnSite(s, kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, BandIndex, BandIndex, BandIndex, BandIndex, s, s, s, s) * this->ComputeEmbeddingOnSite(s, kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, BandIndex, BandIndex, BandIndex, BandIndex, s, s, s, s) * this->ComputeEmbeddingOnSite(s, kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
		    }	  
                  	  
		  if (Index3 == Index4)
		    Tmp *= 0.5;
		  if (Index1 == Index2)
		    Tmp *= 0.5;

		  this->InteractionFactors[i][Index] = 2.0 * FactorU * Tmp;
		  if (this->VPotential != 0.0)		    
		    {
		      int xI, yI, dRx, dRy, sF;
		      Complex translationPhase;
		      int nbrNeighbors=4;
		      int dx[4]={1,-1,0,0};
		      int dy[4]={0,0,1,-1};
		      Tmp=0.0;
		      for (int s=0; s<NbrSublattices; ++s)
			{			  
			  TightBindingModel->DecodeSublatticeIndex(s, xI, yI);
			  
			  for (int n=0; n<nbrNeighbors; ++n)
			    {
			      sF = TightBindingModel->EncodeSublatticeIndex(xI + dx[n], yI + dy[n], dRx, dRy, translationPhase); // calculate final sublattice index.
			      
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, BandIndex, BandIndex, BandIndex, BandIndex, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx3, ky3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, BandIndex, BandIndex, BandIndex, BandIndex, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx4, ky4);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, BandIndex, BandIndex, BandIndex, BandIndex, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx3, ky3);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, BandIndex, BandIndex, BandIndex, BandIndex, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx4, ky4);
			    }
			}

		      if (Index3 == Index4)
			Tmp *= 0.5;
		      if (Index1 == Index2)
			Tmp *= 0.5;
                  	  
		      this->InteractionFactors[i][Index] += 2.0 * FactorV * Tmp;
		    }
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;

  delete [] OneBodyBasis;
}



// compute the phase factor for the embedding of four operators of an on-site two body interaction involving sites on a generic sublattic 
//
// subl = sublattice index
// kx1 = first creation momentum along x for the B site
// ky1 = first creation momentum along y for the B site
// kx2 = second creation momentum along x for the B site
// ky2 = second creation momentum along y for the B site
// kx3 = first annihilation momentum along x for the B site
// ky3 = first annihilation momentum along y for the B site
// kx4 = second annihilation momentum along x for the B site
// ky4 = second annihilation momentum along y for the B site
//
// return value = corresponding matrix element
Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeEmbeddingOnSite(int subl, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  double embeddingX, embeddingY;
  this->TightBindingModel->GetEmbedding(subl, embeddingX, embeddingY);
  double phase = (this->KxFactor * (kx3 + kx4 - kx1 - kx2) * embeddingX + this->KyFactor * (ky3 + ky4 - ky1 - ky2) * embeddingY);
  return Polar(phase);
}


// compute the matrix element for on-site two body interaction involving sites on generic sublattic 
//
// s1 = sublattice index for the first creation operator
// s2 = sublattice index for the second annihilation operator
// kx1 = first creation momentum along x on the first sublattice
// ky1 = first creation momentum along y on the first sublattice
// kx2 = second creation momentum along x on the second sublattice
// ky2 = second creation momentum along y on the second sublattice
// kx3 = first annihilation momentum along x on the second sublattice
// ky3 = first annihilation momentum along y on the second sublattice
// kx4 = second annihilation momentum along x on the first sublattice
// ky4 = second annihilation momentum along y on the first sublattice
//
// return value = corresponding matrix element
 Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeEmbeddingForTwoBodyOperator(int s1, int s2, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  double phase = this->TightBindingModel->GetEmbeddingPhase(s2, this->KxFactor * kx3, this->KyFactor * ky3);
  phase += this->TightBindingModel->GetEmbeddingPhase(s1, this->KxFactor * kx4, this->KyFactor * ky4);
  phase -= this->TightBindingModel->GetEmbeddingPhase(s1, this->KxFactor * kx1, this->KyFactor * ky1);
  phase -= this->TightBindingModel->GetEmbeddingPhase(s2, this->KxFactor * kx2, this->KyFactor * ky2);
  return Polar(phase);
}


// compute the matrix element for on-site two body interaction involving sites on generic sublattic 
//
// dRx = number of unit vector translations along x-direction from EncodeSublatticeIndex (translations back to unit cell)
// dRy = number of unit vector translations along y-direction from EncodeSublatticeIndex (translations back to unit cell)
// kx2 = second creation momentum along x for the translated site
// ky2 = second creation momentum along y for the translated site
// kx3 = first annihilation momentum along x for the translated site
// ky3 = first annihilation momentum along y for the translated site
//
// return value = corresponding matrix element
Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeBlochPhases(int dRx, int dRy, int kx2, int ky2, int kx3, int ky3)
{
  double phase = this->KxFactor * dRx * (kx2-kx3) + this->KyFactor * dRy * (ky2-ky3);
  return Polar(phase);
}

