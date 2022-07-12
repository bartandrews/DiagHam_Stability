////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Yang-Le Wu                            //
//                                                                            //
//            class of Kagome lattice model with interacting particles        //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 04/10/2012                      //
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
#include "Hamiltonian/ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/StringTools.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;

// default constructor
//
ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// tightBindingModel = pointer to the tight binding model
// uPotential = strength of the repulsive two body nearest neighbor interaction
// vPotential = strength of the repulsive two body second nearest neighbor interaction
// wPotential = strength of the repulsive two body third nearest neighbor interaction (or next nearest neighbor density-density potential for bosons)
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// bandIndex = index of the band that has to be partially filled
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX,
															       int nbrSiteY, Abstract2DTightBindingModel* tightBindingModel, 
															       double uPotential, double vPotential, double wPotential, 
															       int bandIndex, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
    this->Particles = particles;
    this->NbrParticles = nbrParticles;
    this->NbrSiteX = nbrSiteX;
    this->NbrSiteY = nbrSiteY;
    this->LzMax = nbrSiteX * nbrSiteY - 1;
    this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
    this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
    this->TightBindingModel = tightBindingModel;
    this->UPotential = uPotential;
    this->VPotential = vPotential;
    this->WPotential = wPotential;
    this->BandIndex = bandIndex;
    this->FlatBand = flatBandFlag;

    this->HamiltonianShift = 0.0;
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

ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::~ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::EvaluateInteractionFactors()
{
    long TotalNbrInteractionFactors = 0;
    double FactorU = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
    double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
    double FactorW = this->WPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

    ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
    if (this->FlatBand == false)
        this->OneBodyInteractionFactors = new double[this->TightBindingModel->GetNbrStatePerBand()];
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
            if (this->FlatBand == false)
                this->OneBodyInteractionFactors[Index] = 0.5 * this->TightBindingModel->GetEnergy(this->BandIndex, Index);
            OneBodyBasis[Index] = this->TightBindingModel->GetOneBodyMatrix(Index);
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
		  int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
		  int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
		  if (Index1 < Index2)
		    ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2)];
		}
        this->SectorIndicesPerSum = new int*[this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
	  {
            if (this->NbrSectorIndicesPerSum[i] > 0)
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
		  int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
		  int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
		  if (Index1 < Index2)
		    {
		      int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2);
		      this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		      this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		      ++this->NbrSectorIndicesPerSum[TmpSum];
		    }
		}
	
        this->InteractionFactors = new Complex* [this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
	  {
            this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
            int Index = 0;
            for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	      {
                int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
                int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
                int kx1, ky1, kx2, ky2;
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index1, kx1, ky1);
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index2, kx2, ky2);
                for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		  {
                    int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                    int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                    int kx3, ky3, kx4, ky4;
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index3, kx3, ky3);
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index4, kx4, ky4);
		    
                    // the InteractionFactors is supposed to be the coefficients to   A+_3 A_1 A+_4 A_2
                    // tricky part: OneBodyBasis[Index] stores the result of LapackDiagonalize
                    // and its [0][_] elements are the COMPLEX CONJUGATE of wave functions < _ |lower band>. (See the end of HermitianMatrix.cc)
		    
                    Complex sumU = 0.;
                    Complex sumV = 0.;
		    
                    sumU += Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
		      * this->ComputeTwoBodyMatrixElementNNAB(kx2, ky2, kx4, ky4);
                    sumU += Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
		      * this->ComputeTwoBodyMatrixElementNNBC(kx2, ky2, kx4, ky4);
                    sumU += Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
		      * this->ComputeTwoBodyMatrixElementNNCA(kx2, ky2, kx4, ky4);
                    sumU -= Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
		      * this->ComputeTwoBodyMatrixElementNNAB(kx2, ky2, kx3, ky3);
                    sumU -= Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
		      * this->ComputeTwoBodyMatrixElementNNBC(kx2, ky2, kx3, ky3);
                    sumU -= Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
		      * this->ComputeTwoBodyMatrixElementNNCA(kx2, ky2, kx3, ky3);
                    sumU -= Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
		      * this->ComputeTwoBodyMatrixElementNNAB(kx1, ky1, kx4, ky4);
                    sumU -= Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
		      * this->ComputeTwoBodyMatrixElementNNBC(kx1, ky1, kx4, ky4);
                    sumU -= Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
		      * this->ComputeTwoBodyMatrixElementNNCA(kx1, ky1, kx4, ky4);
                    sumU += Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
		      * this->ComputeTwoBodyMatrixElementNNAB(kx1, ky1, kx3, ky3);
                    sumU += Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
		      * this->ComputeTwoBodyMatrixElementNNBC(kx1, ky1, kx3, ky3);
                    sumU += Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
		      * this->ComputeTwoBodyMatrixElementNNCA(kx1, ky1, kx3, ky3);
		    
                    sumV += Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
		      * this->ComputeTwoBodyMatrixElementNNNAB(kx2, ky2, kx4, ky4);
                    sumV += Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
		      * this->ComputeTwoBodyMatrixElementNNNBC(kx2, ky2, kx4, ky4);
                    sumV += Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
		      * this->ComputeTwoBodyMatrixElementNNNCA(kx2, ky2, kx4, ky4);
                    sumV -= Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
		      * this->ComputeTwoBodyMatrixElementNNNAB(kx2, ky2, kx3, ky3);
                    sumV -= Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
		      * this->ComputeTwoBodyMatrixElementNNNBC(kx2, ky2, kx3, ky3);
                    sumV -= Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
		      * Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
		      * this->ComputeTwoBodyMatrixElementNNNCA(kx2, ky2, kx3, ky3);
                    sumV -= Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
		      * this->ComputeTwoBodyMatrixElementNNNAB(kx1, ky1, kx4, ky4);
                    sumV -= Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
		      * this->ComputeTwoBodyMatrixElementNNNBC(kx1, ky1, kx4, ky4);
                    sumV -= Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
		      * this->ComputeTwoBodyMatrixElementNNNCA(kx1, ky1, kx4, ky4);
                    sumV += Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
		      * this->ComputeTwoBodyMatrixElementNNNAB(kx1, ky1, kx3, ky3);
                    sumV += Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
		      * this->ComputeTwoBodyMatrixElementNNNBC(kx1, ky1, kx3, ky3);
                    sumV += Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
		      * Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
		      * this->ComputeTwoBodyMatrixElementNNNCA(kx1, ky1, kx3, ky3);
		    
                    this->InteractionFactors[i][Index] = -2.0 * (FactorU * sumU + FactorV * sumV);
		    
                    TotalNbrInteractionFactors++;
                    ++Index;
		  }
	      }
	  }
      }
    else // Boson
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
		  int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
		  int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
		  if (Index1 <= Index2)
		    ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2)];
		}
        this->SectorIndicesPerSum = new int*[this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
	  {
            if (this->NbrSectorIndicesPerSum[i] > 0)
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
		  int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
		  int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
		  if (Index1 <= Index2)
		    {
		      int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2);
		      this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		      this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		      ++this->NbrSectorIndicesPerSum[TmpSum];
		    }
		}
	
        this->InteractionFactors = new Complex* [this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
	  {
            this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
            int Index = 0;
            for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	      {
                int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
                int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
                int kx1, ky1, kx2, ky2;
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index1, kx1, ky1);
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index2, kx2, ky2);
                for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		  {
                    int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                    int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                    int kx3, ky3, kx4, ky4;
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index3, kx3, ky3);
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index4, kx4, ky4);
		    
                    // the InteractionFactors is supposed to be the coefficients to   A+_3 A_1 A+_4 A_2
                    // tricky part: OneBodyBasis[Index] stores the result of LapackDiagonalize
                    // and its [0][_] elements are the COMPLEX CONJUGATE of wave functions < _ |lower band>. (See the end of HermitianMatrix.cc)
		    
                    this->InteractionFactors[i][Index] = FactorU * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]);
                    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]);
                    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0]);
                    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0]);
		    
                    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]);
                    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]);
                    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1]);
                    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1]);
		    
                    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2]);
                    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2]);
                    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2]);
                    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2]);
		    
		    if (this->VPotential != 0.0)
		      {
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementNNAB(kx2, ky2, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementNNAB(kx1, ky1, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementNNAB(kx2, ky2, kx3, ky3);
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementNNAB(kx1, ky1, kx3, ky3);
			
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2] * Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0]) * this->ComputeTwoBodyMatrixElementNNCA(kx2, ky2, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2] * Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0]) * this->ComputeTwoBodyMatrixElementNNCA(kx1, ky1, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2] * Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0]) * this->ComputeTwoBodyMatrixElementNNCA(kx2, ky2, kx3, ky3);
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2] * Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0]) * this->ComputeTwoBodyMatrixElementNNCA(kx1, ky1, kx3, ky3);
			
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementNNBC(kx2, ky2, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementNNBC(kx1, ky1, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementNNBC(kx2, ky2, kx3, ky3);
			this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementNNBC(kx1, ky1, kx3, ky3);
		      }
		    
		    if (this->WPotential != 0.0)
		      {
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementNNNAB(kx2, ky2, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementNNNAB(kx1, ky1, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementNNNAB(kx2, ky2, kx3, ky3);
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementNNNAB(kx1, ky1, kx3, ky3);
			
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2] * Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0]) * this->ComputeTwoBodyMatrixElementNNNCA(kx2, ky2, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2] * Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0]) * this->ComputeTwoBodyMatrixElementNNNCA(kx1, ky1, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2] * Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0]) * this->ComputeTwoBodyMatrixElementNNNCA(kx2, ky2, kx3, ky3);
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2] * Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0]) * this->ComputeTwoBodyMatrixElementNNNCA(kx1, ky1, kx3, ky3);
			
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementNNNBC(kx2, ky2, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementNNNBC(kx1, ky1, kx4, ky4);
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementNNNBC(kx2, ky2, kx3, ky3);
			this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementNNNBC(kx1, ky1, kx3, ky3);
		      }
		    
		    if (Index3 == Index4)
		      this->InteractionFactors[i][Index] *= 0.5;
                    if (Index1 == Index2)
		      this->InteractionFactors[i][Index] *= 0.5;
                    this->InteractionFactors[i][Index] *= 2.0;
		    
                    TotalNbrInteractionFactors++;
                    ++Index;
		  }
	      }
	  }
      }
    cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
    cout << "====================================" << endl;
}

// compute the matrix element for the two body interaction between two sites A and B 
//
// kx2 = annihilation momentum along x for the B site
// ky2 = annihilation momentum along y for the B site
// kx4 = creation momentum along x for the B site
// ky4 = creation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNAB(int kx2, int ky2, int kx4, int ky4)
{
    double dx = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0));
    double dy = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1));
    Complex Tmp = 1 + Phase(-dx);
    return Tmp;
}

// compute the matrix element for the two body interaction between two sites B and C
//
// kx2 = annihilation momentum along x for the C site
// ky2 = annihilation momentum along y for the C site
// kx4 = creation momentum along x for the C site
// ky4 = creation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNBC(int kx2, int ky2, int kx4, int ky4)
{
    double dx = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0));
    double dy = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1));
    Complex Tmp = 1 + Phase(dx - dy);
    return Tmp;
}

// compute the matrix element for the two body interaction between two sites C and A
//
// kx2 = annihilation momentum along x for the A site
// ky2 = annihilation momentum along y for the A site
// kx4 = creation momentum along x for the A site
// ky4 = creation momentum along y for the A site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNCA(int kx2, int ky2, int kx4, int ky4)
{
    double dx = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0));
    double dy = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1));
    Complex Tmp = 1 + Phase(dy);
    return Tmp;
}

// compute the matrix element for the two body interaction between two A sites
//
// kx2 = annihilation momentum along x for the second site
// ky2 = annihilation momentum along y for the second site
// kx4 = creation momentum along x for the second site
// ky4 = creation momentum along y for the second site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNNAB(int kx2, int ky2, int kx4, int ky4)
{
    double dx = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0));
    double dy = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1));
    Complex Tmp = Phase(-dy) + Phase(dy - dx);
    return Tmp;
}

// compute the matrix element for the two body interaction between two B sites
//
// kx2 = annihilation momentum along x for the second site
// ky2 = annihilation momentum along y for the second site
// kx4 = creation momentum along x for the second site
// ky4 = creation momentum along y for the second site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNNBC(int kx2, int ky2, int kx4, int ky4)
{
    double dx = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0));
    double dy = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1));
    Complex Tmp = Phase(dx) + Phase(-dy);
    return Tmp;
}

// compute the matrix element for the two body interaction between two C sites
//
// kx2 = annihilation momentum along x for the second site
// ky2 = annihilation momentum along y for the second site
// kx4 = creation momentum along x for the second site
// ky4 = creation momentum along y for the second site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNNCA(int kx2, int ky2, int kx4, int ky4)
{
    double dx = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0));
    double dy = (double)(this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1) - this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1));
    Complex Tmp = Phase(dx) + Phase(dy - dx);
    return Tmp;
}
