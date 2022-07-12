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
//         in the single band approximation and five body interaction         // 
//                                                                            //
//                        last modification : 07/08/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian.h"
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



// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = strength of the repulsive five body neareast neighbor interaction
// vPotential = strength of the repulsive two body neareast neighbor interaction
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian::ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, 
																     int nbrSiteY, double uPotential, double vPotential, 
																     Abstract2DTightBindingModel* tightBindingModel, 
																     bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->HamiltonianShift = 0.0;
  this->NBodyValue = 5;
  this->ComputePhaseArray();

  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->TightBindingModel = tightBindingModel;
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

ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian::~ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  if (this->FlatBand == false)
    this->OneBodyInteractionFactors = new double [this->NbrSiteX * this->NbrSiteY];
  this->ComputeOneBodyMatrices(OneBodyBasis);
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
					{
					  for (int kx5 = 0; kx5 < this->NbrSiteX; ++kx5)
					    for (int ky5 = 0; ky5 < this->NbrSiteY; ++ky5) 
					      {
						int Index5 = (kx5 * this->NbrSiteY) + ky5;
						if (Index4 < Index5)
						  ++this->NbrNBodySectorIndicesPerSum[(((kx1 + kx2 + kx3 + kx4 + kx5) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3 + ky4 + ky5) % this->NbrSiteY)]; 
					      }
					}
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
					  for (int kx5 = 0; kx5 < this->NbrSiteX; ++kx5)
					    for (int ky5 = 0; ky5 < this->NbrSiteY; ++ky5) 
					      {
						int Index5 = (kx5 * this->NbrSiteY) + ky5;
						if (Index4 < Index5)
						  {
						    int TmpSum = (((kx1 + kx2 + kx3 + kx4 + kx5) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3 + ky4 + ky5) % this->NbrSiteY);
						    this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * 5] = Index1;
						    this->NBodySectorIndicesPerSum[TmpSum][1 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 5)] = Index2;
						    this->NBodySectorIndicesPerSum[TmpSum][2 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 5)] = Index3;
						    this->NBodySectorIndicesPerSum[TmpSum][3 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 5)] = Index4;
						    this->NBodySectorIndicesPerSum[TmpSum][4 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 5)] = Index5;
						    ++this->NbrNBodySectorIndicesPerSum[TmpSum];    
						  }
					      }
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


      double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), this->NBodyValue - 1);
      double FactorV = this->VPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), this->NBodyValue - 1);
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      int KxIn[5];
      int KyIn[5];
      int KxInPerm[5];
      int KyInPerm[5];
      int KxOut[5];
      int KyOut[5];
      int IndexIn[5];
      int IndexOut[5];
      
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
      
      Complex TmpABBBB;
      Complex TmpBAAAA;



      int TmpLargestSector = 0;
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	if (this->NbrNBodySectorIndicesPerSum[i] > TmpLargestSector)
	  TmpLargestSector = this->NbrNBodySectorIndicesPerSum[i];

      Complex* TmpABBBBIn = new Complex[TmpLargestSector];
      Complex* TmpABBBBOut = new Complex[TmpLargestSector];
      Complex* TmpBAAAAIn = new Complex[TmpLargestSector];
      Complex* TmpBAAAAOut = new Complex[TmpLargestSector];

      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      IndexIn[0] = this->NBodySectorIndicesPerSum[i][j1 * 5];
	      IndexIn[1] = this->NBodySectorIndicesPerSum[i][(j1 * 5) + 1];
	      IndexIn[2] = this->NBodySectorIndicesPerSum[i][(j1 * 5) + 2];
	      IndexIn[3] = this->NBodySectorIndicesPerSum[i][(j1 * 5) + 3];
	      IndexIn[4] = this->NBodySectorIndicesPerSum[i][(j1 * 5) + 4];
	      KxIn[0] = IndexIn[0] / this->NbrSiteY;
	      KyIn[0] = IndexIn[0] % this->NbrSiteY;
	      KxIn[1] = IndexIn[1] / this->NbrSiteY;
	      KyIn[1] = IndexIn[1] % this->NbrSiteY;
	      KxIn[2] = IndexIn[2] / this->NbrSiteY;
	      KyIn[2] = IndexIn[2] % this->NbrSiteY;
	      KxIn[3] = IndexIn[3] / this->NbrSiteY;
	      KyIn[3] = IndexIn[3] % this->NbrSiteY;	      
	      KxIn[4] = IndexIn[4] / this->NbrSiteY;
	      KyIn[4] = IndexIn[4] % this->NbrSiteY;	      

	      Complex TmpABBBBIn2 = 0.0;
	      Complex TmpABBBBOut2 = 0.0;
	      Complex TmpBAAAAIn2 = 0.0;
	      Complex TmpBAAAAOut2 = 0.0;

	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];

		  TmpABBBBIn2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[4]]][0][1]) * this->ComputeFiveBodyMatrixElementABBBBIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]], KxIn[TmpPerm[3]], KyIn[TmpPerm[3]], KxIn[TmpPerm[4]], KyIn[TmpPerm[4]]);
		  TmpBAAAAIn2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[4]]][0][0]) * this->ComputeFiveBodyMatrixElementBAAAAIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]], KxIn[TmpPerm[3]], KyIn[TmpPerm[3]], KxIn[TmpPerm[4]], KyIn[TmpPerm[4]]);
		  
		  TmpABBBBOut2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[4]]][0][1] * this->ComputeFiveBodyMatrixElementABBBBOut(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]], KxIn[TmpPerm[3]], KyIn[TmpPerm[3]], KxIn[TmpPerm[4]], KyIn[TmpPerm[4]]);
		  TmpBAAAAOut2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[4]]][0][0] * this->ComputeFiveBodyMatrixElementBAAAAOut(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]], KxIn[TmpPerm[3]], KyIn[TmpPerm[3]], KxIn[TmpPerm[4]], KyIn[TmpPerm[4]]);
		}

	      TmpABBBBIn[j1] = TmpABBBBIn2;
	      TmpABBBBOut[j1] = TmpABBBBOut2;
	      TmpBAAAAIn[j1] = TmpBAAAAIn2;
	      TmpBAAAAOut[j1] = TmpBAAAAOut2;

	    }

	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  this->NBodyInteractionFactors[i][Index] = -2.0 * FactorU * ((TmpABBBBIn[j1] * TmpABBBBOut[j2]) + (TmpBAAAAIn[j1] * TmpBAAAAOut[j2]));
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}


      delete[] TmpABBBBIn;
      delete[] TmpABBBBOut;
      delete[] TmpBAAAAIn;
      delete[] TmpBAAAAOut;

// 	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
// 		{
// 		  IndexOut[0] = this->NBodySectorIndicesPerSum[i][j2 * 5];
// 		  IndexOut[1] = this->NBodySectorIndicesPerSum[i][(j2 * 5) + 1];
// 		  IndexOut[2] = this->NBodySectorIndicesPerSum[i][(j2 * 5) + 2];
// 		  IndexOut[3] = this->NBodySectorIndicesPerSum[i][(j2 * 5) + 3];
// 		  IndexOut[4] = this->NBodySectorIndicesPerSum[i][(j2 * 5) + 4];
// 		  KxOut[0] = IndexOut[0] / this->NbrSiteY;
// 		  KyOut[0] = IndexOut[0] % this->NbrSiteY;
// 		  KxOut[1] = IndexOut[1] / this->NbrSiteY;
// 		  KyOut[1] = IndexOut[1] % this->NbrSiteY;
// 		  KxOut[2] = IndexOut[2] / this->NbrSiteY;
// 		  KyOut[2] = IndexOut[2] % this->NbrSiteY;
// 		  KxOut[3] = IndexOut[3] / this->NbrSiteY;
// 		  KyOut[3] = IndexOut[3] % this->NbrSiteY;
// 		  KxOut[4] = IndexOut[4] / this->NbrSiteY;
// 		  KyOut[4] = IndexOut[4] % this->NbrSiteY;
// 		  Complex Tmp = 0.0;
// 		  TmpABBBBIn = 0.0;
// 		  TmpABBBBOut = 0.0;
// 		  TmpBAAAAIn = 0.0;
// 		  TmpBAAAAOut = 0.0;
//  		  for (int l1 = 0; l1 < NbrPermutations; ++l1)
//  		    {
//  		      int* TmpPerm = Permutations[l1];
//  		      TmpABBBBIn += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[4]]][0][1]) * this->ComputeFiveBodyMatrixElementABBBBIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]], KxIn[TmpPerm[3]], KyIn[TmpPerm[3]], KxIn[TmpPerm[4]], KyIn[TmpPerm[4]]);
// 		      TmpBAAAAIn += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[4]]][0][0]) * this->ComputeFiveBodyMatrixElementBAAAAIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]], KxIn[TmpPerm[3]], KyIn[TmpPerm[3]], KxIn[TmpPerm[4]], KyIn[TmpPerm[4]]);

// 		      TmpABBBBOut += PermutationSign[l1] * OneBodyBasis[IndexOut[TmpPerm[0]]][0][0] * OneBodyBasis[IndexOut[TmpPerm[1]]][0][1] * OneBodyBasis[IndexOut[TmpPerm[2]]][0][1] * OneBodyBasis[IndexOut[TmpPerm[3]]][0][1] * OneBodyBasis[IndexOut[TmpPerm[4]]][0][1] * this->ComputeFiveBodyMatrixElementABBBBOut(KxOut[TmpPerm[0]], KyOut[TmpPerm[0]], KxOut[TmpPerm[1]], KyOut[TmpPerm[1]], KxOut[TmpPerm[2]], KyOut[TmpPerm[2]], KxOut[TmpPerm[3]], KyOut[TmpPerm[3]], KxOut[TmpPerm[4]], KyOut[TmpPerm[4]]);
// 		      TmpBAAAAOut += PermutationSign[l1] * OneBodyBasis[IndexOut[TmpPerm[0]]][0][1] * OneBodyBasis[IndexOut[TmpPerm[1]]][0][0] * OneBodyBasis[IndexOut[TmpPerm[2]]][0][0] * OneBodyBasis[IndexOut[TmpPerm[3]]][0][0] * OneBodyBasis[IndexOut[TmpPerm[4]]][0][0] * this->ComputeFiveBodyMatrixElementBAAAAOut(KxOut[TmpPerm[0]], KyOut[TmpPerm[0]], KxOut[TmpPerm[1]], KyOut[TmpPerm[1]], KxOut[TmpPerm[2]], KyOut[TmpPerm[2]], KxOut[TmpPerm[3]], KyOut[TmpPerm[3]], KxOut[TmpPerm[4]], KyOut[TmpPerm[4]]);
		      
// 		      // 		      int* TmpPerm = Permutations[l1];
// // 		      KxInPerm[0] = KxIn[TmpPerm[0]];
// // 		      KyInPerm[0] = KyIn[TmpPerm[0]];
// // 		      KxInPerm[1] = KxIn[TmpPerm[1]];
// // 		      KyInPerm[1] = KyIn[TmpPerm[1]];
// // 		      KxInPerm[2] = KxIn[TmpPerm[2]];
// // 		      KyInPerm[2] = KyIn[TmpPerm[2]];
// // 		      KxInPerm[3] = KxIn[TmpPerm[3]];
// // 		      KyInPerm[3] = KyIn[TmpPerm[3]];
// // 		      KxInPerm[4] = KxIn[TmpPerm[4]];
// // 		      KyInPerm[4] = KyIn[TmpPerm[4]];
// // 		      TmpABBBB = PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[4]]][0][1]);
// // 		      TmpBAAAA = PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[3]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[4]]][0][0]);
		      
// // 		      for (int l2 = 0; l2 < NbrPermutations; ++l2)
// // 			{
// // 			  int* TmpPerm2 = Permutations[l2];			  

// // 			  Tmp += TmpABBBB * PermutationSign [l2] * OneBodyBasis[IndexOut[TmpPerm2[0]]][0][0] * OneBodyBasis[IndexOut[TmpPerm2[1]]][0][1] * OneBodyBasis[IndexOut[TmpPerm2[2]]][0][1] * OneBodyBasis[IndexOut[TmpPerm2[3]]][0][1] * OneBodyBasis[IndexOut[TmpPerm2[4]]][0][1] * this->ComputeFiveBodyMatrixElementABBBB(KxInPerm[0], KyInPerm[0], KxInPerm[1], KyInPerm[1], KxInPerm[2], KyInPerm[2], KxInPerm[3], KyInPerm[3], KxInPerm[4], KyInPerm[4], KxOut[TmpPerm2[0]], KyOut[TmpPerm2[0]], KxOut[TmpPerm2[1]], KyOut[TmpPerm2[1]], KxOut[TmpPerm2[2]], KyOut[TmpPerm2[2]], KxOut[TmpPerm2[3]], KyOut[TmpPerm2[3]], KxOut[TmpPerm2[4]], KyOut[TmpPerm2[4]]);
			  
// // 			  Tmp += TmpBAAAA * PermutationSign [l2] * OneBodyBasis[IndexOut[TmpPerm2[0]]][0][1] * OneBodyBasis[IndexOut[TmpPerm2[1]]][0][0] * OneBodyBasis[IndexOut[TmpPerm2[2]]][0][0] * OneBodyBasis[IndexOut[TmpPerm2[3]]][0][0] * OneBodyBasis[IndexOut[TmpPerm2[4]]][0][0] * this->ComputeFiveBodyMatrixElementBAAAA(KxInPerm[0], KyInPerm[0], KxInPerm[1], KyInPerm[1], KxInPerm[2], KyInPerm[2], KxInPerm[3], KyInPerm[3], KxInPerm[4], KyInPerm[4], KxOut[TmpPerm2[0]], KyOut[TmpPerm2[0]], KxOut[TmpPerm2[1]], KyOut[TmpPerm2[1]], KxOut[TmpPerm2[2]], KyOut[TmpPerm2[2]], KxOut[TmpPerm2[3]], KyOut[TmpPerm2[3]], KxOut[TmpPerm2[4]], KyOut[TmpPerm2[4]]);
// // 			}
// //		    }
// //		  
// //		  this->NBodyInteractionFactors[i][Index] = -2.0 * FactorU * Tmp;

// 		    }
// 		  this->NBodyInteractionFactors[i][Index] = -2.0 * FactorU * ((TmpABBBBIn * TmpABBBBOut) + (TmpBAAAAIn * TmpBAAAAOut));																																												
// 		  TotalNbrInteractionFactors++;
// 		  ++Index;
// 		}
// 	    }

      gettimeofday (&(TotalEndingTime), 0);
      Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
		     ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
      cout << "element generation done in  " << Dt << "s" << endl;

      delete[] PermutationSign;
      for (int i = 0; i < NbrPermutations; ++i)
	delete[] Permutations[i];
      delete[] Permutations;
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

