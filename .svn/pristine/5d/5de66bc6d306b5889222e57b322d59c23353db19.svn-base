////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//    class of pyrochlore slab lattice model with interacting particles       //
//      in the single band approximation and four body interaction            // 
//                                                                            //
//                        last modification : 18/06/2012                      //
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
#include "Hamiltonian/ParticleOnLatticePyrochloreSlabLatticeSingleBandFourBodyHamiltonian.h"
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
ParticleOnLatticePyrochloreSlabLatticeSingleBandFourBodyHamiltonian::ParticleOnLatticePyrochloreSlabLatticeSingleBandFourBodyHamiltonian()
{
  this->BandIndex = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellX = number of sites in the x direction
// nbrCellY = number of sites in the y direction
// uPotential = strength of the repulsive two body onsite interaction
// vPotential = strength of the repulsive two body neareast neighbor interaction
// bandIndex = index of the band to consider
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticePyrochloreSlabLatticeSingleBandFourBodyHamiltonian::ParticleOnLatticePyrochloreSlabLatticeSingleBandFourBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, 
															 double uPotential, double vPotential,  
															 Abstract2DTightBindingModel* tightBindingModel, int bandIndex, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrCellX;
  this->NbrSiteY = nbrCellY;
  this->NBodyValue = 4;
  this->LzMax = nbrCellX * nbrCellY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);

  this->HamiltonianShift = 0.0;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
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

ParticleOnLatticePyrochloreSlabLatticeSingleBandFourBodyHamiltonian::~ParticleOnLatticePyrochloreSlabLatticeSingleBandFourBodyHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticePyrochloreSlabLatticeSingleBandFourBodyHamiltonian::EvaluateInteractionFactors()
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

      double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), this->NBodyValue - 1);
      double FactorV = this->VPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), this->NBodyValue - 1);
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

      Complex* TmpAABBIn1 = new Complex[TmpLargestSector];
      Complex* TmpAABBIn2 = new Complex[TmpLargestSector];
      Complex* TmpAABBOut1 = new Complex[TmpLargestSector];
      Complex* TmpAABBOut2 = new Complex[TmpLargestSector];
      
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
	      
	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];
		}
	      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
		{
		  for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		    {
		      this->NBodyInteractionFactors[i][Index] = 2.0 * FactorU * ((TmpAABBIn1[j1] * TmpAABBOut1[j2]) + (TmpAABBIn2[j1] * TmpAABBOut2[j2]));
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
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
		  if (Index1 <= Index2)
		    {
		      for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
			for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
			  {
			    int Index3 = (kx3 * this->NbrSiteY) + ky3;
			    if (Index2 <= Index3)
			      {
				for (int kx4 = 0; kx4 < this->NbrSiteX; ++kx4)
				  for (int ky4 = 0; ky4 < this->NbrSiteY; ++ky4) 
				    {
				      int Index4 = (kx4 * this->NbrSiteY) + ky4;
				      if (Index3 <= Index4)
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
		  if (Index1 <= Index2)
		    {
		      for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
			for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
			  {
			    int Index3 = (kx3 * this->NbrSiteY) + ky3;
			    if (Index2 <= Index3)
			      {
				for (int kx4 = 0; kx4 < this->NbrSiteX; ++kx4)
				  for (int ky4 = 0; ky4 < this->NbrSiteY; ++ky4) 
				    {
				      int Index4 = (kx4 * this->NbrSiteY) + ky4;
				      if (Index3 <= Index4)
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

      int KxIn[4];
      int KyIn[4];
      int IndexIn[4];

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
	      
	      for (int Site = 0; Site < this->TightBindingModel->GetNbrBands(); ++Site)
		{
		  TmpIn2[Site] = 0.0;
		  TmpOut2[Site] = 0.0;
		}
	      double SymmetryFactor = 1.0;
	      int SymmetryIndex = 1;
	      if (IndexIn[0] == IndexIn[1])
		{
		  ++SymmetryIndex;
		}
	      if (IndexIn[1] == IndexIn[2])
		{
		  ++SymmetryIndex;
		}
	      if (IndexIn[2] == IndexIn[3])
		{
		  ++SymmetryIndex;
		}
	      for (int k = 1; k <= SymmetryIndex; ++k)
		SymmetryFactor *= (double) k;
	      SymmetryFactor = 1.0 / SymmetryFactor;
	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];
		  for (int Site = 0; Site < this->TightBindingModel->GetNbrBands(); ++Site)
		    {
		      TmpIn2[Site] +=  Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][BandIndex][Site] * OneBodyBasis[IndexIn[TmpPerm[1]]][BandIndex][Site] * OneBodyBasis[IndexIn[TmpPerm[2]]][BandIndex][Site] * OneBodyBasis[IndexIn[TmpPerm[3]]][BandIndex][Site]);
		      TmpOut2[Site] +=  OneBodyBasis[IndexIn[TmpPerm[0]]][BandIndex][Site] * OneBodyBasis[IndexIn[TmpPerm[1]]][BandIndex][Site] * OneBodyBasis[IndexIn[TmpPerm[2]]][BandIndex][Site] * OneBodyBasis[IndexIn[TmpPerm[3]]][BandIndex][Site];
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
		  for (int Site = 0; Site < this->TightBindingModel->GetNbrBands(); ++Site)
		    Tmp += TmpIn[Site][j1] * TmpOut[Site][j2];
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
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;

  delete [] OneBodyBasis;
}
