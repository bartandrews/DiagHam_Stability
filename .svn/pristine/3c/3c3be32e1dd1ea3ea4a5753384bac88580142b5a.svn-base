 ////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                          class author: Yang-Le Wu                          //
//                                                                            //
//     class of 2-orbital square lattice model with interacting particles     //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 11/09/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandHamiltonian.h"
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



// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// tightBindingModel = pointer to the tight binding model
// uPotential = strength of the repulsive two body neareast neighbor interaction
// vPotential = strength of the repulsive two body next neareast neighbor interaction
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandHamiltonian::ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY,
        Abstract2DTightBindingModel* tightBindingModel, double uPotential, double vPotential, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);

  this->HamiltonianShift = 0.0;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;

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

ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandHamiltonian::~ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;

  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
  if (this->FlatBand == false)
      this->OneBodyInteractionFactors = new double[this->TightBindingModel->GetNbrStatePerBand()];
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
          int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
          if (this->FlatBand == false)
              this->OneBodyInteractionFactors[Index] = 0.5 * this->TightBindingModel->GetEnergy(0, Index);
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
      double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if (this->FlatBand == false)
	FactorU *= this->UPotential;
      double FactorV = 0.5 * this->VPotential / ((double) (this->NbrSiteX * this->NbrSiteY));

      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;

                  Complex sumU = 0.;
                  Complex sumV = 0.;

                  sumU += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1];
                  sumU -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1];
                  sumU -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1];
                  sumU += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1];


                  sumV += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1])
                         *(Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeTwoBodyMatrixElementNN(kx2, ky2, kx4, ky4);
                  sumV -= (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                         *(Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1])
                        * this->ComputeTwoBodyMatrixElementNN(kx2, ky2, kx3, ky3);
                  sumV -= (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1])
                         *(Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeTwoBodyMatrixElementNN(kx1, ky1, kx4, ky4);
                  sumV += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                         *(Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1])
                        * this->ComputeTwoBodyMatrixElementNN(kx1, ky1, kx3, ky3);

		  this->InteractionFactors[i][Index] = -2.0 * (FactorU * sumU + FactorV * sumV);
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the two body interaction between two NN sites
//
// kx2 = annihilation momentum along x for the second site
// ky2 = annihilation momentum along y for the second site
// kx4 = creation momentum along x for the second site
// ky4 = creation momentum along y for the second site
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2SquareLatticeTwoOrbitalSingleBandHamiltonian::ComputeTwoBodyMatrixElementNN(int kx2, int ky2, int kx4, int ky4)
{
  double dx = ((double)(kx2-kx4)) * this->KxFactor;
  double dy = ((double)(ky2-ky4)) * this->KyFactor;
  Complex Tmp = Phase(dx) + Phase(dy);
  return Tmp;
}
