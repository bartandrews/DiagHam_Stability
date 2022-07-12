///////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                           class author: Yang-Le Wu                         //
//                                                                            //
//     class of 2-orbital square lattice model with interacting particles     //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 12/09/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian.h"
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

ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = strength of the repulsive two body nearest neighbor interaction
// vPotential = strength of the repulsive two body second nearest neighbor interaction
// wPotential = strength of the repulsive three body nearest neighbor interaction
// sPotential = strength of the repulsive three body next-to-nearest neighbor interaction
// t1 = imag part of the inter-orbital hopping amplitude between nearest neighbors along the x direction
// t2 = the inter-orbital hopping amplitude between nearest neighbors along the y direction
// t3 = the intra-orbital hopping amplitude between nearest neighbors
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
        double uPotential,Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NBodyValue = 3;
  this->FlatBand = flatBandFlag;
this->TightBindingModel = tightBindingModel;
  this->HamiltonianShift = 0.0;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->UPotential = uPotential;
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
      cout << "fast = ";
      PrintMemorySize(cout, TmpMemory)<< endl;
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian::~ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian::EvaluateInteractionFactors()
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
	  this->OneBodyInteractionFactors[Index] = 0.5 * this->TightBindingModel->GetEnergy(0, Index);
	OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
      }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
  {
      /*this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;      
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      if (this->TwoBodyFlag == true)
      {
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
          double FactorU = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
          double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
          for (int i = 0; i < this->NbrSectorSums; ++i)
          {
              this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
              int Index = 0;
              for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1) // annihilation operators
              {
                  int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
                  int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
                  int kx1 = Index1 / this->NbrSiteY;
                  int ky1 = Index1 % this->NbrSiteY;
                  int kx2 = Index2 / this->NbrSiteY;
                  int ky2 = Index2 % this->NbrSiteY;
                  for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2) // creation operators
                  {
                      int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                      int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                      int kx3 = Index3 / this->NbrSiteY;
                      int ky3 = Index3 % this->NbrSiteY;
                      int kx4 = Index4 / this->NbrSiteY;
                      int ky4 = Index4 % this->NbrSiteY;
                      // the InteractionFactors is supposed to be the coefficients to   A+_3 A_1 A+_4 A_2
                      // tricky part: OneBodyBasis[Index] stores the result of LapackDiagonalize
                      // and its [0][_] elements are the COMPLEX CONJUGATE of wave functions < _ |lower band>. (See the end of HermitianMatrix.cc)

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
      cout << "nbr 2-body interaction = " << TotalNbrInteractionFactors << endl;
      TotalNbrInteractionFactors=0;


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

      double FactorW = this->WPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), this->NBodyValue - 1);
      double FactorS = this->SPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), this->NBodyValue - 1);
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1) // annihilation operators
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
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2) // creation operators
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

                  Complex sumW=0.;
                  Complex sumS=0.;

                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6);
                  sumW -= (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5);
                  sumW -= (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6);
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4);
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5);
                  sumW -= (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4);
                  sumW -= (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6);
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5);
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6);
                  sumW -= (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4);
                  sumW -= (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5);
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4);
                  sumW -= (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6);
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5);
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6);
                  sumW -= (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4);
                  sumW -= (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5);
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4);
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6);
                  sumW -= (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5);
                  sumW -= (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6);
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4);
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5);
                  sumW -= (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4);
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6);
                  sumW -= (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5);
                  sumW -= (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6);
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4);
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5);
                  sumW -= (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4);
                  sumW -= (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6);
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5);
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6);
                  sumW -= (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4);
                  sumW -= (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5);
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                         * this->ComputeThreeBodyMatrixElementNN(kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4);

                  this->NBodyInteractionFactors[i][Index] = 2.0 * (FactorW * sumW + FactorS * sumS);

		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
      cout << "nbr 3-body interaction = " << TotalNbrInteractionFactors << endl;
      cout << "====================================" << endl;
      */
  }
  else
  {
    int BandIndex = 0;
    this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;      
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      
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

      
      int KxIn[3];
      int KyIn[3];
      int IndexIn[3];

      int** Permutations = 0; 
      double* PermutationSign = 0; 
      int NbrPermutations = this->ComputePermutations(Permutations, PermutationSign);

      int TmpLargestSector = 0;
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	if (this->NbrNBodySectorIndicesPerSum[i] > TmpLargestSector)
	  TmpLargestSector = this->NbrNBodySectorIndicesPerSum[i];

      Complex** TmpIn  = new Complex*[this->TightBindingModel->GetNbrBands()];
      Complex** TmpOut  = new Complex*[this->TightBindingModel->GetNbrBands()];
      Complex* TmpIn2  = new Complex[this->TightBindingModel->GetNbrBands()];
      Complex* TmpOut2  = new Complex[this->TightBindingModel->GetNbrBands()];
      for (int i = 0; i < this->TightBindingModel->GetNbrBands(); ++i)
	{
	  TmpIn[i] = new Complex[TmpLargestSector];
	  TmpOut[i] = new Complex[TmpLargestSector];
	}

      double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      IndexIn[0] = this->NBodySectorIndicesPerSum[i][j1 * 3];
	      IndexIn[1] = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 1];
	      IndexIn[2] = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 2];
	      KxIn[0] = IndexIn[0] / this->NbrSiteY;
	      KyIn[0] = IndexIn[0] % this->NbrSiteY;
	      KxIn[1] = IndexIn[1] / this->NbrSiteY;
	      KyIn[1] = IndexIn[1] % this->NbrSiteY;
	      KxIn[2] = IndexIn[2] / this->NbrSiteY;
	      KyIn[2] = IndexIn[2] % this->NbrSiteY;
	      
	      for (int Site = 0; Site < this->TightBindingModel->GetNbrBands(); ++Site)
		{
		  TmpIn2[Site] = 0.0;
		  TmpOut2[Site] = 0.0;
		}
	      double SymmetryFactor = 1.0;
	      if ((IndexIn[0] == IndexIn[1]) && (IndexIn[1] == IndexIn[2]))
		{
		  SymmetryFactor = 1.0 / 6.0;
		}
	      else
		{
		  if ((IndexIn[0] == IndexIn[1]) || (IndexIn[1] == IndexIn[2]) || (IndexIn[0] == IndexIn[2]))
		    {
		      SymmetryFactor = 0.5;
		    }
		}
	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];
		  for (int Site = 0; Site < this->TightBindingModel->GetNbrBands(); ++Site)
		    {
		      TmpIn2[Site] +=  Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][BandIndex][Site] * OneBodyBasis[IndexIn[TmpPerm[1]]][BandIndex][Site] * OneBodyBasis[IndexIn[TmpPerm[2]]][BandIndex][Site]);
		      TmpOut2[Site] +=  OneBodyBasis[IndexIn[TmpPerm[0]]][BandIndex][Site] * OneBodyBasis[IndexIn[TmpPerm[1]]][BandIndex][Site] * OneBodyBasis[IndexIn[TmpPerm[2]]][BandIndex][Site];
		    }
		}
	      for (int Site = 0; Site < this->TightBindingModel->GetNbrBands(); ++Site)
		{
		  TmpIn[Site][j1] = TmpIn2[Site] * SymmetryFactor;			
		  TmpOut[Site][j1] = TmpOut2[Site] * SymmetryFactor;			
		}
	    }
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  Complex Tmp = 0.0;
		  for (int Site1 = 0; Site1 < this->TightBindingModel->GetNbrBands(); ++Site1)
		  {
		    for (int Site2 = 0; Site2 < this->TightBindingModel->GetNbrBands(); ++Site2)
		      Tmp += TmpIn[Site1][j1] * TmpOut[Site2][j2];
		  }
		  this->NBodyInteractionFactors[i][Index] = 2.0 * FactorU * Tmp;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
      delete[] TmpIn2;
      delete[] TmpOut2;
      for (int i = 0; i < this->TightBindingModel->GetNbrBands(); ++i)
	{
	  delete[] TmpIn[i];
	  delete[] TmpOut[i];
	}
      delete[] TmpIn;
      delete[] TmpOut;
    
    
      
     /* double SymmetryFactor1 = 1.0;
      double SymmetryFactor2 = 1.0;
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1) // annihilation operators
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
	      
	      if(Index3 == Index2)
	      {
		if(Index2 == Index1)
		  SymmetryFactor1 = 6.0;
		else
		  SymmetryFactor1 = 2.0;
	      }
	      else
	      {
		if(Index2 == Index1)
		  SymmetryFactor1 = 2.0;
		else
		  SymmetryFactor1 = 1.0;
	      }
	      
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2) // creation operators
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

                  Complex sumW=0;
		  Complex sumS=0;
		  
		  if(Index6 == Index5)
	      {
		if(Index5 == Index4)
		  SymmetryFactor2 = 6.0;
		else
		  SymmetryFactor2 = 2.0;
	      }
	      else
	      {
		if(Index5 == Index4)
		  SymmetryFactor2 = 2.0;
		else
		  SymmetryFactor2 = 1.0;
	      }
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1];
			 
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1];
			 
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1];
			  
			 
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] ;
			 
			
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] ;
			 
                  
		    sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] ;
                  
	      sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1];
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] ;
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1];
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] ;
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1];
                  sumW += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] ;
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] ;
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] ;
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1];
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1];
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] ;
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1];
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1];
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1];
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1];
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1];
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] ;
                  sumW += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1];
                  sumW += (Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                         + Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                         * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                         * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1];

		   sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]);
			   
			
                  this->NBodyInteractionFactors[i][Index] = 2.0 * FactorU*( 3*sumW +  sumS)/(SymmetryFactor1 * SymmetryFactor2);
		  
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}*/
      cout << "nbr 3-body interaction = " << TotalNbrInteractionFactors << endl;
      cout << "====================================" << endl;
  }

}
