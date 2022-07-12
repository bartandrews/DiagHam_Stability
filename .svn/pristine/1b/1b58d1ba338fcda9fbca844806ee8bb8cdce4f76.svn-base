 ////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//       class of checkerboard lattice model with interacting particles       //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 13/07/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
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

ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = strength of the repulsive three body neareast neighbor interaction
// vPotential = strength of the repulsive two body neareast neighbor interaction
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, 
																       int nbrSiteY, double uPotential, double vPotential, 
																       Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NBodyValue = 3;
  this->ComputePhaseArray();

  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;


  this->TwoBodyFlag = false;
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

ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::~ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian()
{
  delete[] this->XPhaseTable;
  delete[] this->XHalfPhaseTable;
  delete[] this->YPhaseTable;
  delete[] this->YHalfPhaseTable;
}

// evaluate all interaction factors
//   

void ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactors = new double [this->TightBindingModel->GetNbrStatePerBand()];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	if (this->FlatBand == false)
	  this->OneBodyInteractionFactors[Index] = this->TightBindingModel->GetEnergy(0, Index);
	OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
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
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 < Index2)
		  ++this->NbrSectorIndicesPerSum[(((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)];    
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
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 < Index2)
		  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }

      if (this->TwoBodyFlag == false)
	{
	  for (int i = 0; i < this->NbrSectorSums; ++i)	  
	    this->NbrSectorIndicesPerSum[i] = 0;
	  this->InteractionFactors = new Complex* [this->NbrSectorSums];
	}

      this->NbrNBodySectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrNBodySectorIndicesPerSum = new int[this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	this->NbrNBodySectorIndicesPerSum[i] = 0;      
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
	    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
		  {
		    int Index1 = (kx1 * this->NbrSiteY) + ky1;
		    int Index2 = (kx2 * this->NbrSiteY) + ky2;
		    int Index3 = (kx3 * this->NbrSiteY) + ky3;
		    if ((Index1 < Index2) && (Index2 < Index3))
		      ++this->NbrNBodySectorIndicesPerSum[(((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY)];    
		  }
      this->NBodySectorIndicesPerSum = new int* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  if (this->NbrNBodySectorIndicesPerSum[i]  > 0)
	    {
	      this->NBodySectorIndicesPerSum[i] = new int[this->NBodyValue * this->NbrNBodySectorIndicesPerSum[i]];      
	      this->NbrNBodySectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
	    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
		  {
		    int Index1 = (kx1 * this->NbrSiteY) + ky1;
		    int Index2 = (kx2 * this->NbrSiteY) + ky2;
		    int Index3 = (kx3 * this->NbrSiteY) + ky3;
		    if ((Index1 < Index2) && (Index2 < Index3))
		      {
			int TmpSum = (((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY);
			this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * 3] = Index1;
			this->NBodySectorIndicesPerSum[TmpSum][1 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = Index2;
			this->NBodySectorIndicesPerSum[TmpSum][2 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = Index3;
			++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		      }
		  }


      double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      double FactorV = this->VPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->NBodySectorIndicesPerSum[i][j1 * 3];
	      int Index2 = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 1];
	      int Index3 = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 2];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  int Index4 = this->NBodySectorIndicesPerSum[i][j2 * 3];
		  int Index5 = this->NBodySectorIndicesPerSum[i][(j2 * 3) + 1];
		  int Index6 = this->NBodySectorIndicesPerSum[i][(j2 * 3) + 2];
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
		  int kx5 = Index5 / this->NbrSiteY;
		  int ky5 = Index5 % this->NbrSiteY;
		  int kx6 = Index6 / this->NbrSiteY;
		  int ky6 = Index6 % this->NbrSiteY;


                  this->NBodyInteractionFactors[i][Index] = Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5);



                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6, kx5, ky5);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4, kx6, ky6);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5, kx4, ky4);
                  this->NBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4, kx5, ky5);


// 		  double Sign = 1.0;
// 		  if ((((kx1 + kx2 + kx3 - kx4 - kx5 - kx6) / this->NbrSiteX) % 2) == 1)
// 		    {
// 		      Sign *= -1.0;
// 		    }
// 		  if ((((ky1 + ky2 + ky3 - ky4 - ky5 - ky6) / this->NbrSiteY) % 2) == 1)
// 		    {
// 		      Sign *= -1.0;
// 		    }
		  
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] -= Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4);
//                   this->NBodyInteractionFactors[i][Index] += Sign * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5);

		  this->NBodyInteractionFactors[i][Index] *= 2.0 * FactorU;

		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  else
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
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 <= Index2)
		  ++this->NbrSectorIndicesPerSum[(((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)];    
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
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 <= Index2)
		  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }

      if (this->TwoBodyFlag == false)
	{
	  for (int i = 0; i < this->NbrSectorSums; ++i)	  
	    this->NbrSectorIndicesPerSum[i] = 0;
	  this->InteractionFactors = new Complex* [this->NbrSectorSums];
	}

      this->NbrNBodySectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrNBodySectorIndicesPerSum = new int[this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	this->NbrNBodySectorIndicesPerSum[i] = 0;      
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
	    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
		  {
		    int Index1 = (kx1 * this->NbrSiteY) + ky1;
		    int Index2 = (kx2 * this->NbrSiteY) + ky2;
		    int Index3 = (kx3 * this->NbrSiteY) + ky3;
		    if ((Index1 <= Index2) && (Index2 <= Index3))
		      ++this->NbrNBodySectorIndicesPerSum[(((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY)];    
		  }
      this->NBodySectorIndicesPerSum = new int* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  if (this->NbrNBodySectorIndicesPerSum[i]  > 0)
	    {
	      this->NBodySectorIndicesPerSum[i] = new int[this->NBodyValue * this->NbrNBodySectorIndicesPerSum[i]];      
	      this->NbrNBodySectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
	    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
		  {
		    int Index1 = (kx1 * this->NbrSiteY) + ky1;
		    int Index2 = (kx2 * this->NbrSiteY) + ky2;
		    int Index3 = (kx3 * this->NbrSiteY) + ky3;
		    if ((Index1 <= Index2) && (Index2 <= Index3))
		      {
			int TmpSum = (((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY);
			this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * 3] = Index1;
			this->NBodySectorIndicesPerSum[TmpSum][1 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = Index2;
			this->NBodySectorIndicesPerSum[TmpSum][2 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = Index3;
			++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		      }
		  }


      double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      double FactorV = this->VPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->NBodySectorIndicesPerSum[i][j1 * 3];
	      int Index2 = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 1];
	      int Index3 = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 2];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  int Index4 = this->NBodySectorIndicesPerSum[i][j2 * 3];
		  int Index5 = this->NBodySectorIndicesPerSum[i][(j2 * 3) + 1];
		  int Index6 = this->NBodySectorIndicesPerSum[i][(j2 * 3) + 2];
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
		  int kx5 = Index5 / this->NbrSiteY;
		  int ky5 = Index5 % this->NbrSiteY;
		  int kx6 = Index6 / this->NbrSiteY;
		  int ky6 = Index6 % this->NbrSiteY;

		  Complex Sum;
		  this->NBodyInteractionFactors[i][Index] = 0.0;

                  Sum = Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx1, ky1, kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx2, ky2, kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBB(kx3, ky3, kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4, kx5, ky5);

		  this->NBodyInteractionFactors[i][Index] += Sum;

                  Sum = Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx1, ky1, kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx2, ky2, kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAA(kx3, ky3, kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4, kx5, ky5);

		  this->NBodyInteractionFactors[i][Index] += Sum;		  

                  Sum = Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx1, ky1, kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx2, ky2, kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteAAB(kx3, ky3, kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4, kx5, ky5);

		  //		  this->NBodyInteractionFactors[i][Index] += Sum;		 

                  Sum = Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx1, ky1, kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx2, ky2, kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6, kx5, ky5);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4, kx6, ky6);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5, kx4, ky4);
                  Sum += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementOnSiteABB(kx3, ky3, kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4, kx5, ky5);


		  //		  this->NBodyInteractionFactors[i][Index] += Sum;		  

		  //		  this->NBodyInteractionFactors[i][Index] *= 2.0 * FactorU;

		  if ((Index1 == Index2) && (Index1 == Index3))
		    {
		      this->NBodyInteractionFactors[i][Index] /= 6.0;
		    }
		  else
		    {
		      if ((Index1 == Index2) || (Index1 == Index3))
			{
			  this->NBodyInteractionFactors[i][Index] *= 0.5;
			}
		    }

		  if ((Index4 == Index5) && (Index4 == Index6))
		    {
		      this->NBodyInteractionFactors[i][Index] /= 6.0;
		    }
		  else
		    {
		      if ((Index4 == Index5) || (Index4 == Index6))
			{
			  this->NBodyInteractionFactors[i][Index] *= 0.5;
			    }
		    }

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
// kx1 = momentum along x for the A site
// ky1 = momentum along y for the A site
// kx2 = momentum along x for the B site
// ky2 = momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp = 2.0 * (cos (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) - ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY))))) 
		       + cos (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) + ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY))))));
  return Tmp;
}

// compute the matrix element for the three body interaction between two sites B and one site A
//
// kx1 = creation momentum along x for the first B site
// ky1 = creation momentum along y for the first B site
// kx2 = creation momentum along x for the second B site
// ky2 = creation momentum along y for the second B site
// kx3 = annihilation momentum along x for the first B site
// ky3 = annihilation momentum along y for the first B site
// kx4 = annihilation momentum along x for the second B site
// ky4 = annihilation momentum along y for the second B site
// return value = corresponding matrix element

Complex ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 2.0 * (cos (M_PI * ((((double) ((kx1 - kx3) - (kx2 - kx4))) / ((double) this->NbrSiteX)) + ((((double) ((ky1 - ky3) + (ky2 - ky4))) / ((double) this->NbrSiteY))))) 
		       + cos (M_PI * ((((double) ((kx1 - kx3) + (kx2 - kx4))) / ((double) this->NbrSiteX)) - ((((double) ((ky1 - ky3) - (ky2 - ky4))) / ((double) this->NbrSiteY))))));
  //  Tmp = 0.0;
  return Tmp;
}

// compute the matrix element for the three body interaction between one site B and two sites A
//
// kx1 = creation momentum along x for the first A site
// ky1 = creation momentum along y for the first A site
// kx2 = creation momentum along x for the second A site
// ky2 = creation momentum along y for the second A site
// kx3 = annihilation momentum along x for the first A site
// ky3 = annihilation momentum along y for the first A site
// kx4 = annihilation momentum along x for the second A site
// ky4 = annihilation momentum along y for the second A site
// return value = corresponding matrix element

Complex ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementBAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  Complex Tmp = Phase (2.0 * M_PI * ((double) (kx3 - kx6)) / ((double) this->NbrSiteX)); 
  Tmp += Phase ((2.0 * M_PI * ((double) (kx2 - kx5 + kx3 - kx6)) / ((double) this->NbrSiteX)) +
		 (2.0 * M_PI * ((double) (ky3 - ky6)) / ((double) this->NbrSiteY))); 
  Tmp += Phase ((2.0 * M_PI * ((double) (ky2 - ky5 + ky3 - ky6)) / ((double) this->NbrSiteY)) +
		(2.0 * M_PI * ((double) (kx2 - kx5)) / ((double) this->NbrSiteX))); 
  Tmp += Phase (2.0 * M_PI * ((double) (ky2 - ky5)) / ((double) this->NbrSiteY)); 
  Tmp *= Phase ((M_PI * ((double) (kx1 - kx4)) / ((double) this->NbrSiteX)) +
		(M_PI * ((double) (ky1 - ky4)) / ((double) this->NbrSiteY))); 
  return Tmp;
}

// compute the matrix element for the on-site three body interaction related to sites A
//
// kx1 = first creation momentum along x for the first A site
// ky1 = first creation momentum along y for the first A site
// kx2 = second creation momentum along x for the second A site
// ky2 = second creation momentum along y for the second A site
// kx3 = third creation momentum along x for the second A site
// ky3 = third creation momentum along y for the second A site
// kx4 = first annihilation momentum along x for the first A site
// ky4 = first annihilation momentum along y for the first A site
// kx5 = second annihilation momentum along x for the second A site
// ky5 = second annihilation momentum along y for the second A site
// kx6 = third annihilation momentum along x for the second A site
// ky6 = third annihilation momentum along y for the second A site
// return value = corresponding matrix element

Complex ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the on-site three body interaction related to sites A
//
 // kx1 = first creation momentum along x for the first A site
// ky1 = first creation momentum along y for the first A site
// kx2 = second creation momentum along x for the second A site
// ky2 = second creation momentum along y for the second A site
// kx3 = third creation momentum along x for the second A site
// ky3 = third creation momentum along y for the second A site
// kx4 = first annihilation momentum along x for the first A site
// ky4 = first annihilation momentum along y for the first A site
// kx5 = second annihilation momentum along x for the second A site
// ky5 = second annihilation momentum along y for the second A site
// kx6 = third annihilation momentum along x for the second A site
// ky6 = third annihilation momentum along y for the second A site
// return value = corresponding matrix element

Complex ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  Complex Tmp = Phase (0.5 * ((((double) (kx6 + kx5 + kx4 - kx3 - kx2 - kx1)) * this->KxFactor)
			      + ((((double) (ky6 + ky5 + ky4 - ky3 - ky2 - ky1)) * this->KyFactor))));
  return Tmp;
}


// compute the matrix element for the on-site three body interaction with two particles on a A site and one on a B site
//
// kx1 = first creation momentum along x for the A site
// ky1 = first creation momentum along y for the A site
// kx2 = second creation momentum along x for the A site
// ky2 = second creation momentum along y for the A site
// kx3 = creation momentum along x for the B site
// ky3 = creation momentum along y for the B site
// kx4 = first annihilation momentum along x for the A site
// ky4 = first annihilation momentum along y for the A site
// kx5 = second annihilation momentum along x for the A site
// ky5 = second annihilation momentum along y for the sA site
// kx6 = annihilation momentum along x for the B site
// ky6 = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteAAB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  Complex Tmp = 1.0;
  Tmp += Phase (((double) (kx3 - kx6)) * this->KxFactor);
  Tmp += Phase (((double) (ky3 - ky6)) * this->KyFactor);
  Tmp += Phase ((((double) (kx3 - kx6)) * this->KxFactor) + (((double) (ky3 - ky6)) * this->KyFactor));
  Tmp *= Phase (0.5 * ((((double) (kx6 - kx3)) * this->KxFactor)
		       + ((((double) (ky6 - ky3)) * this->KyFactor))));
  return Tmp;
}

// compute the matrix element for the on-site three body interaction with two particles on a B site and one on a A site
//
// kx1 = creation momentum along x for the A site
// ky1 = creation momentum along y for the A site
// kx2 = first creation momentum along x for the B site
// ky2 = first creation momentum along y for the B site
// kx3 = second creation momentum along x for the B site
// ky3 = second creation momentum along y for the B site
// kx4 = annihilation momentum along x for the A site
// ky4 = annihilation momentum along y for the A site
// kx5 = first annihilation momentum along x for the B site
// ky5 = first annihilation momentum along y for the B site
// kx6 = second annihilation momentum along x for the B site
// ky6 = second annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  Complex Tmp = 1.0;
  Tmp += Phase (((double) (kx4 - kx1)) * this->KxFactor);
  Tmp += Phase (((double) (ky4 - ky1)) * this->KyFactor);
  Tmp += Phase ((((double) (kx4 - kx1)) * this->KxFactor) + (((double) (ky4 - ky1)) * this->KyFactor));
  Tmp *= Phase (0.5 * ((((double) (kx6 + kx5 - kx3 - kx2)) * this->KxFactor)
		       + ((((double) (ky6 + ky5 - ky3 - ky2)) * this->KyFactor))));
  return Tmp;
}




// compute all the phase precalculation arrays 
//

void ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ComputePhaseArray()
{
  this->XPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteX];
  this->XHalfPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteX];
  this->XPhaseTableShift = this->NBodyValue * this->NbrSiteX;
  for (int i = -this->XPhaseTableShift; i < this->XPhaseTableShift; ++i)
    {
      this->XPhaseTable[this->XPhaseTableShift + i] = Phase(this->KxFactor * ((double) i));
      this->XHalfPhaseTable[this->XPhaseTableShift + i] = Phase(0.5 * this->KxFactor * ((double) i));
    }
  this->YPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteY];
  this->YHalfPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteY];
  this->YPhaseTableShift = this->NBodyValue * this->NbrSiteY;
  for (int i = -this->YPhaseTableShift; i < this->YPhaseTableShift; ++i)
    {
      this->YPhaseTable[this->YPhaseTableShift + i] = Phase(this->KyFactor * ((double) i));
      this->YHalfPhaseTable[this->YPhaseTableShift + i] = Phase(0.5 * this->KyFactor * ((double) i));
    }
}
