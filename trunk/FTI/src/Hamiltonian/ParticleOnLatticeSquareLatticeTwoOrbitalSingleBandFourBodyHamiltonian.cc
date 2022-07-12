 ////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                           class author: Yang-Le Wu                         //
//                                                                            //
//     class of 2-orbital square lattice model with interacting particles     //
//         in the single band approximation and four body interaction         // 
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
#include "Hamiltonian/ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian.h"
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

ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian::ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian()
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
// wPotential = strength of the repulsive three body nearest neighbor interaction
// sPotential = strength of the repulsive three body next-to-nearest neighbor interaction
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian::ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY,
        Abstract2DTightBindingModel* tightBindingModel, double uPotential, double vPotential, double wPotential, double sPotential, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NBodyValue = 4;
  this->ComputePhaseArray();
  this->TightBindingModel = tightBindingModel;

  this->HamiltonianShift = 0.0;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->WPotential = wPotential;
  this->SPotential = sPotential;
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

ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian::~ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandFourBodyHamiltonian::EvaluateInteractionFactors()
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

      timeval TotalStartingTime;
      timeval TotalEndingTime;
      gettimeofday (&(TotalStartingTime), 0);

      this->NbrNBodySectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrNBodySectorIndicesPerSum = new int[this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	this->NbrNBodySectorIndicesPerSum[i] = 0;      
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	  {
	    int Index1 = (kx1 * this->NbrSiteY) + ky1;
	    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		{
		  int Index2 = (kx2 * this->NbrSiteY) + ky2;
		  if (Index1 < Index2)
		    {
		      for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
			for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
			  {
			    int Index3 = (kx3 * this->NbrSiteY) + ky3;
			    if (Index2 < Index3)
			      {
				for (int kx4 = 0; kx4 < this->NbrSiteX; ++kx4)
				  for (int ky4 = 0; ky4 < this->NbrSiteY; ++ky4) 
				    {
				      int Index4 = (kx4 * this->NbrSiteY) + ky4;
				      if (Index3 < Index4)
					++this->NbrNBodySectorIndicesPerSum[(((kx1 + kx2 + kx3 + kx4) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3 + ky4) % this->NbrSiteY)];    
				    }
			      }
			  }
		    }
		}
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
	for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	  {
	    int Index1 = (kx1 * this->NbrSiteY) + ky1;
	    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		{
		  int Index2 = (kx2 * this->NbrSiteY) + ky2;
		  if (Index1 < Index2)
		    {
		      for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
			for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
			  {
			    int Index3 = (kx3 * this->NbrSiteY) + ky3;
			    if (Index2 < Index3)
			      {
				for (int kx4 = 0; kx4 < this->NbrSiteX; ++kx4)
				  for (int ky4 = 0; ky4 < this->NbrSiteY; ++ky4) 
				    {
				      int Index4 = (kx4 * this->NbrSiteY) + ky4;
				      if (Index3 < Index4)
					{
					  int TmpSum = (((kx1 + kx2 + kx3 + kx4) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3 + ky4) % this->NbrSiteY);
					  this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * 4] = Index1;
					  this->NBodySectorIndicesPerSum[TmpSum][1 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 4)] = Index2;
					  this->NBodySectorIndicesPerSum[TmpSum][2 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 4)] = Index3;
					  this->NBodySectorIndicesPerSum[TmpSum][3 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 4)] = Index4;
					  ++this->NbrNBodySectorIndicesPerSum[TmpSum];    
					}
				    }
			      }
			  }
		    }
		}
	  }

      gettimeofday (&(TotalEndingTime), 0);
      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
      cout << "index generation done in  " << Dt << "s" << endl;
      gettimeofday (&(TotalStartingTime), 0);

      double FactorW = this->WPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), this->NBodyValue - 1);
      double FactorS = this->SPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), this->NBodyValue - 1);
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      int KxIn[4];
      int KyIn[4];
      int KxInPerm[4];
      int KyInPerm[4];
      int KxOut[4];
      int KyOut[4];
      int IndexIn[4];
      int IndexOut[4];
      
      int NbrPermutations = 1;
      for (int i = 1; i <= this->NBodyValue; ++i)
	NbrPermutations *= i;
      int** Permutations = new int*[NbrPermutations]; 
      double* PermutationSign = new double[NbrPermutations]; 
      Permutations[0] = new int [this->NBodyValue];
      for (int i = 0; i < this->NBodyValue; ++i)
	Permutations[0][i] = i;
      PermutationSign[0] = 1.0;
      double TmpSign = 1.0;
      for (int i = 1; i < NbrPermutations; ++i)
	{
	  Permutations[i] = new int [this->NBodyValue];
	  for (int j = 0; j < this->NBodyValue; ++j)
	    Permutations[i][j] = Permutations[i - 1][j];
	  int* TmpArrayPerm = Permutations[i];
	  int Pos1 = this->NBodyValue - 1;
	  while (TmpArrayPerm[Pos1 - 1] >= TmpArrayPerm[Pos1])
	    --Pos1;
	  --Pos1;
	  int Pos2 = this->NBodyValue - 1;      
	  while (TmpArrayPerm[Pos2] <= TmpArrayPerm[Pos1])
	    --Pos2;
	  int TmpIndex = TmpArrayPerm[Pos1];
	  TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	  TmpArrayPerm[Pos2] = TmpIndex;
	  TmpSign *= -1.0;
	  Pos2 = this->NBodyValue - 1;   
	  Pos1++;
	  while (Pos1 < Pos2)
	    {
	      TmpIndex = TmpArrayPerm[Pos1];
	      TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	      TmpArrayPerm[Pos2] = TmpIndex;
	      ++Pos1;
	      --Pos2;
	      TmpSign *= -1.0;
	    }
	  PermutationSign[i] = TmpSign;
	}
      
      int TmpLargestSector = 0;
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	if (this->NbrNBodySectorIndicesPerSum[i] > TmpLargestSector)
	  TmpLargestSector = this->NbrNBodySectorIndicesPerSum[i];

      Complex* TmpNNIn  = new Complex[TmpLargestSector];
      Complex* TmpNNOut = new Complex[TmpLargestSector];
      Complex* TmpNNNIn  = new Complex[TmpLargestSector];
      Complex* TmpNNNOut = new Complex[TmpLargestSector];

      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      IndexIn[0] = this->NBodySectorIndicesPerSum[i][j1 * 4];
	      IndexIn[1] = this->NBodySectorIndicesPerSum[i][(j1 * 4) + 1];
	      IndexIn[2] = this->NBodySectorIndicesPerSum[i][(j1 * 4) + 2];
	      IndexIn[3] = this->NBodySectorIndicesPerSum[i][(j1 * 4) + 3];
	      KxIn[0] = IndexIn[0] / this->NbrSiteY;
	      KyIn[0] = IndexIn[0] % this->NbrSiteY;
	      KxIn[1] = IndexIn[1] / this->NbrSiteY;
	      KyIn[1] = IndexIn[1] % this->NbrSiteY;
	      KxIn[2] = IndexIn[2] / this->NbrSiteY;
	      KyIn[2] = IndexIn[2] % this->NbrSiteY;
	      KxIn[3] = IndexIn[3] / this->NbrSiteY;
	      KyIn[3] = IndexIn[3] % this->NbrSiteY;	      

              Complex SumNNIn  = 0.0;
              Complex SumNNOut = 0.0;
              Complex SumNNNIn  = 0.0;
              Complex SumNNNOut = 0.0;

	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];

		  SumNNIn += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] 
                                                      * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][1])
                      * this->ComputeFourBodyMatrixElementNNIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]], KxIn[TmpPerm[3]], KyIn[TmpPerm[3]]);

		  SumNNNIn += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] 
                                                       * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][1])
                      * this->ComputeFourBodyMatrixElementNNNIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]], KxIn[TmpPerm[3]], KyIn[TmpPerm[3]]);

		  SumNNOut += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] 
                                                  * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][1]
                      * Conj(this->ComputeFourBodyMatrixElementNNIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]], KxIn[TmpPerm[3]], KyIn[TmpPerm[3]]));

		  SumNNNOut += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1]
                                                   * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][1]
                      * Conj(this->ComputeFourBodyMatrixElementNNNIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]], KxIn[TmpPerm[3]], KyIn[TmpPerm[3]]));

		}
	      TmpNNIn [j1] = SumNNIn;
	      TmpNNOut[j1] = SumNNOut;
	      TmpNNNIn [j1] = SumNNNIn;
	      TmpNNNOut[j1] = SumNNNOut;
	    }
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  this->NBodyInteractionFactors[i][Index] = 2.0 * (FactorW * TmpNNIn[j1] * TmpNNOut[j2] + FactorS * TmpNNNIn[j1] * TmpNNNOut[j2]);
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}

      delete[] TmpNNIn;
      delete[] TmpNNOut;
      delete[] TmpNNNIn;
      delete[] TmpNNNOut;

      delete[] PermutationSign;
      for (int i = 0; i < NbrPermutations; ++i)
	delete[] Permutations[i];
      delete[] Permutations;

      gettimeofday (&(TotalEndingTime), 0);
      Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
		     ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
      cout << "element generation done in  " << Dt << "s" << endl;

      cout << "nbr " << this->NBodyValue << "-body interaction = " << TotalNbrInteractionFactors << endl;
      cout << "====================================" << endl;
  }
}

