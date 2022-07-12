 ////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of ruby lattice model with interacting particles          //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 25/10/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian.h"
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

ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian()
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

ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX,   int nbrSiteY, double uPotential, double vPotential, 
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

  this->HamiltonianShift = 0.0;
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

ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::~ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::EvaluateInteractionFactors()
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

      int** Permutations = 0; 
      double* PermutationSign = 0; 
      int NbrPermutations = this->ComputePermutations(Permutations, PermutationSign);
      
      int TmpLargestSector = 0;
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	if (this->NbrNBodySectorIndicesPerSum[i] > TmpLargestSector)
	  TmpLargestSector = this->NbrNBodySectorIndicesPerSum[i];

      int KxIn[3];
      int KyIn[3];
      int IndexIn[3];

      double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      double FactorV = this->VPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];

      Complex* TmpA1A3A5In = new Complex[TmpLargestSector];
      Complex* TmpA1A3A5Out = new Complex[TmpLargestSector];
      Complex* TmpA2A4A6In = new Complex[TmpLargestSector];
      Complex* TmpA2A4A6Out = new Complex[TmpLargestSector];
      Complex* TmpA1A2A5In = new Complex[TmpLargestSector];
      Complex* TmpA1A2A5Out = new Complex[TmpLargestSector];
      Complex* TmpA1A4A5In = new Complex[TmpLargestSector];
      Complex* TmpA1A4A5Out = new Complex[TmpLargestSector];
      Complex* TmpA2A4A5In = new Complex[TmpLargestSector];
      Complex* TmpA2A4A5Out = new Complex[TmpLargestSector];
      Complex* TmpA1A2A4In = new Complex[TmpLargestSector];
      Complex* TmpA1A2A4Out = new Complex[TmpLargestSector];
      Complex* TmpA1A3A6In = new Complex[TmpLargestSector];
      Complex* TmpA1A3A6Out = new Complex[TmpLargestSector];
      Complex* TmpA1A3A4In = new Complex[TmpLargestSector];
      Complex* TmpA1A3A4Out = new Complex[TmpLargestSector];
      Complex* TmpA1A4A6In = new Complex[TmpLargestSector];
      Complex* TmpA1A4A6Out = new Complex[TmpLargestSector];
      Complex* TmpA3A4A6In = new Complex[TmpLargestSector];
      Complex* TmpA3A4A6Out = new Complex[TmpLargestSector];
      Complex* TmpA2A5A6In = new Complex[TmpLargestSector];
      Complex* TmpA2A5A6Out = new Complex[TmpLargestSector];
      Complex* TmpA2A3A6In = new Complex[TmpLargestSector];
      Complex* TmpA2A3A6Out = new Complex[TmpLargestSector];
      Complex* TmpA2A3A5In = new Complex[TmpLargestSector];
      Complex* TmpA2A3A5Out = new Complex[TmpLargestSector];
      Complex* TmpA6A3A5In = new Complex[TmpLargestSector];
      Complex* TmpA6A3A5Out = new Complex[TmpLargestSector];

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

    	      Complex TmpA1A3A5In2 = 0.0;
	      Complex TmpA1A3A5Out2 = 0.0;
    	      Complex TmpA2A4A6In2 = 0.0;
	      Complex TmpA2A4A6Out2 = 0.0;
    	      Complex TmpA1A2A5In2 = 0.0;
	      Complex TmpA1A2A5Out2 = 0.0;
    	      Complex TmpA1A4A5In2 = 0.0;
	      Complex TmpA1A4A5Out2 = 0.0;
    	      Complex TmpA2A4A5In2 = 0.0;
	      Complex TmpA2A4A5Out2 = 0.0;
    	      Complex TmpA1A2A4In2 = 0.0;
	      Complex TmpA1A2A4Out2 = 0.0;
    	      Complex TmpA1A3A6In2 = 0.0;
	      Complex TmpA1A3A6Out2 = 0.0;
    	      Complex TmpA1A3A4In2 = 0.0;
	      Complex TmpA1A3A4Out2 = 0.0;
    	      Complex TmpA1A4A6In2 = 0.0;
	      Complex TmpA1A4A6Out2 = 0.0;
    	      Complex TmpA3A4A6In2 = 0.0;
	      Complex TmpA3A4A6Out2 = 0.0;
    	      Complex TmpA2A5A6In2 = 0.0;
	      Complex TmpA2A5A6Out2 = 0.0;
    	      Complex TmpA2A3A6In2 = 0.0;
	      Complex TmpA2A3A6Out2 = 0.0;
    	      Complex TmpA2A3A5In2 = 0.0;
	      Complex TmpA2A3A5Out2 = 0.0;
    	      Complex TmpA6A3A5In2 = 0.0;
	      Complex TmpA6A3A5Out2 = 0.0;

	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];
 		  TmpA1A3A5In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4]) * this->ComputeThreeBodyMatrixElementA1A3A5In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A4A6In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5]) * this->ComputeThreeBodyMatrixElementA2A4A6In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A2A5In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4]) * this->ComputeThreeBodyMatrixElementA1A2A5In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A4A5In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4]) * this->ComputeThreeBodyMatrixElementA1A4A5In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A4A5In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4]) * this->ComputeThreeBodyMatrixElementA2A4A5In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A2A4In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][3]) * this->ComputeThreeBodyMatrixElementA1A2A4In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A3A6In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5]) * this->ComputeThreeBodyMatrixElementA1A3A6In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A3A4In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][3]) * this->ComputeThreeBodyMatrixElementA1A3A4In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A4A6In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5]) * this->ComputeThreeBodyMatrixElementA1A4A6In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA3A4A6In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5]) * this->ComputeThreeBodyMatrixElementA3A4A6In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A5A6In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][4] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5]) * this->ComputeThreeBodyMatrixElementA2A5A6In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A3A6In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5]) * this->ComputeThreeBodyMatrixElementA2A3A6In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A3A5In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4]) * this->ComputeThreeBodyMatrixElementA2A3A5In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA6A3A5In2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][5] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4]) * this->ComputeThreeBodyMatrixElementA6A3A5In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);

 		  TmpA1A3A5Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4] * this->ComputeThreeBodyMatrixElementA1A3A5Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A4A6Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5] * this->ComputeThreeBodyMatrixElementA2A4A6Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A2A5Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4] * this->ComputeThreeBodyMatrixElementA1A2A5Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A4A5Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4] * this->ComputeThreeBodyMatrixElementA1A4A5Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A4A5Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4] * this->ComputeThreeBodyMatrixElementA2A4A5Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A2A4Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][3] * this->ComputeThreeBodyMatrixElementA1A2A4Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A3A6Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5] * this->ComputeThreeBodyMatrixElementA1A3A6Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A3A4Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][3] * this->ComputeThreeBodyMatrixElementA1A3A4Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA1A4A6Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5] * this->ComputeThreeBodyMatrixElementA1A4A6Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA3A4A6Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5] * this->ComputeThreeBodyMatrixElementA3A4A6Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A5A6Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][4] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5] * this->ComputeThreeBodyMatrixElementA2A5A6Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A3A6Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5] * this->ComputeThreeBodyMatrixElementA2A3A6Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A3A5Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4] * this->ComputeThreeBodyMatrixElementA2A3A5Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA6A3A5Out2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][5] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4] * this->ComputeThreeBodyMatrixElementA6A3A5Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
		}
	      TmpA1A3A5In[j1] = TmpA1A3A5In2;
	      TmpA1A3A5Out[j1] = TmpA1A3A5Out2;
	      TmpA2A4A6In[j1] = TmpA2A4A6In2;
	      TmpA2A4A6Out[j1] = TmpA2A4A6Out2;
	      TmpA1A2A5In[j1] = TmpA1A2A5In2;
	      TmpA1A2A5Out[j1] = TmpA1A2A5Out2;
	      TmpA1A4A5In[j1] = TmpA1A4A5In2;
	      TmpA1A4A5Out[j1] = TmpA1A4A5Out2;
	      TmpA2A4A5In[j1] = TmpA2A4A5In2;
	      TmpA2A4A5Out[j1] = TmpA2A4A5Out2;
	      TmpA1A2A4In[j1] = TmpA1A2A4In2;
	      TmpA1A2A4Out[j1] = TmpA1A2A4Out2;
	      TmpA1A3A6In[j1] = TmpA1A3A6In2;
	      TmpA1A3A6Out[j1] = TmpA1A3A6Out2;
	      TmpA1A3A4In[j1] = TmpA1A3A4In2;
	      TmpA1A3A4Out[j1] = TmpA1A3A4Out2;
	      TmpA1A4A6In[j1] = TmpA1A4A6In2;
	      TmpA1A4A6Out[j1] = TmpA1A4A6Out2;
	      TmpA3A4A6In[j1] = TmpA3A4A6In2;
	      TmpA3A4A6Out[j1] = TmpA3A4A6Out2;
	      TmpA2A5A6In[j1] = TmpA2A5A6In2;
	      TmpA2A5A6Out[j1] = TmpA2A5A6Out2;
	      TmpA2A3A6In[j1] = TmpA2A3A6In2;
	      TmpA2A3A6Out[j1] = TmpA2A3A6Out2;
	      TmpA2A3A5In[j1] = TmpA2A3A5In2;
	      TmpA2A3A5Out[j1] = TmpA2A3A5Out2;
	      TmpA6A3A5In[j1] = TmpA6A3A5In2;
	      TmpA6A3A5Out[j1] = TmpA6A3A5Out2;
	    }
	  
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  this->NBodyInteractionFactors[i][Index] = 2.0 * FactorU * ((TmpA1A3A5In[j1] * TmpA1A3A5Out[j2])
									     + (TmpA2A4A6In[j1] * TmpA2A4A6Out[j2])
									     + (TmpA1A2A5In[j1] * TmpA1A2A5Out[j2])
									     + (TmpA1A4A5In[j1] * TmpA1A4A5Out[j2])
									     + (TmpA2A4A5In[j1] * TmpA2A4A5Out[j2])
									     + (TmpA1A2A4In[j1] * TmpA1A2A4Out[j2])
									     + (TmpA1A3A6In[j1] * TmpA1A3A6Out[j2])
									     + (TmpA1A3A4In[j1] * TmpA1A3A4Out[j2])
									     + (TmpA1A4A6In[j1] * TmpA1A4A6Out[j2])
									     + (TmpA3A4A6In[j1] * TmpA3A4A6Out[j2])
									     + (TmpA2A5A6In[j1] * TmpA2A5A6Out[j2])
									     + (TmpA2A3A6In[j1] * TmpA2A3A6Out[j2])
									     + (TmpA2A3A5In[j1] * TmpA2A3A5Out[j2])
									     + (TmpA6A3A5In[j1] * TmpA6A3A5Out[j2]));
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
      delete[] TmpA1A3A5In;
      delete[] TmpA1A3A5Out;
      delete[] TmpA2A4A6In;
      delete[] TmpA2A4A6Out;
      delete[] TmpA1A2A5In;
      delete[] TmpA1A2A5Out;
      delete[] TmpA1A4A5In;
      delete[] TmpA1A4A5Out;
      delete[] TmpA2A4A5In;
      delete[] TmpA2A4A5Out;
      delete[] TmpA1A2A4In;
      delete[] TmpA1A2A4Out;
      delete[] TmpA1A3A6In;
      delete[] TmpA1A3A6Out;
      delete[] TmpA1A3A4In;
      delete[] TmpA1A3A4Out;
      delete[] TmpA1A4A6In;
      delete[] TmpA1A4A6Out;
      delete[] TmpA3A4A6In;
      delete[] TmpA3A4A6Out;
      delete[] TmpA2A5A6In;
      delete[] TmpA2A5A6Out;
      delete[] TmpA2A3A6In;
      delete[] TmpA2A3A6Out;
      delete[] TmpA2A3A5In;
      delete[] TmpA2A3A5Out;
      delete[] TmpA6A3A5In;
      delete[] TmpA6A3A5Out;
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

      Complex* TmpA1A1A1In = new Complex[TmpLargestSector];
      Complex* TmpA1A1A1Out = new Complex[TmpLargestSector];
      Complex* TmpA2A2A2In = new Complex[TmpLargestSector];
      Complex* TmpA2A2A2Out = new Complex[TmpLargestSector];
      Complex* TmpA3A3A3In = new Complex[TmpLargestSector];
      Complex* TmpA3A3A3Out = new Complex[TmpLargestSector];
      Complex* TmpA4A4A4In = new Complex[TmpLargestSector];
      Complex* TmpA4A4A4Out = new Complex[TmpLargestSector];
      Complex* TmpA5A5A5In = new Complex[TmpLargestSector];
      Complex* TmpA5A5A5Out = new Complex[TmpLargestSector];
      Complex* TmpA6A6A6In = new Complex[TmpLargestSector];
      Complex* TmpA6A6A6Out = new Complex[TmpLargestSector];

      double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      double FactorV = this->VPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
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

	      Complex TmpA1A1A1In2 = 0.0;
	      Complex TmpA1A1A1Out2 = 0.0;
	      Complex TmpA2A2A2In2 = 0.0;
	      Complex TmpA2A2A2Out2 = 0.0;
	      Complex TmpA3A3A3In2 = 0.0;
	      Complex TmpA3A3A3Out2 = 0.0;
	      Complex TmpA4A4A4In2 = 0.0;
	      Complex TmpA4A4A4Out2 = 0.0;
	      Complex TmpA5A5A5In2 = 0.0;
	      Complex TmpA5A5A5Out2 = 0.0;
	      Complex TmpA6A6A6In2 = 0.0;
	      Complex TmpA6A6A6Out2 = 0.0;
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

 		  TmpA1A1A1In2 += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0]) * this->ComputeThreeBodyMatrixElementOnSiteA1A1A1In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A2A2In2 += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1]) * this->ComputeThreeBodyMatrixElementOnSiteA2A2A2In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA3A3A3In2 += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2]) * this->ComputeThreeBodyMatrixElementOnSiteA3A3A3In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA4A4A4In2 += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][3]) * this->ComputeThreeBodyMatrixElementOnSiteA4A4A4In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA5A5A5In2 += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][4] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][4] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4]) * this->ComputeThreeBodyMatrixElementOnSiteA5A5A5In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA6A6A6In2 += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][5] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][5] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5]) * this->ComputeThreeBodyMatrixElementOnSiteA6A6A6In(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
		  
 		  TmpA1A1A1Out2 += OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0] * this->ComputeThreeBodyMatrixElementOnSiteA1A1A1Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA2A2A2Out2 += OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1] * this->ComputeThreeBodyMatrixElementOnSiteA2A2A2Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA3A3A3Out2 += OneBodyBasis[IndexIn[TmpPerm[0]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2] * this->ComputeThreeBodyMatrixElementOnSiteA3A3A3Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA4A4A4Out2 += OneBodyBasis[IndexIn[TmpPerm[0]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][3] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][3] * this->ComputeThreeBodyMatrixElementOnSiteA4A4A4Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA5A5A5Out2 += OneBodyBasis[IndexIn[TmpPerm[0]]][0][4] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][4] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][4] * this->ComputeThreeBodyMatrixElementOnSiteA5A5A5Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpA6A6A6Out2 += OneBodyBasis[IndexIn[TmpPerm[0]]][0][5] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][5] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][5] * this->ComputeThreeBodyMatrixElementOnSiteA6A6A6Out(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);

		}

 	      TmpA1A1A1In[j1] =  TmpA1A1A1In2 * SymmetryFactor;
 	      TmpA1A1A1Out[j1] = TmpA1A1A1Out2 * SymmetryFactor;
 	      TmpA2A2A2In[j1] =  TmpA2A2A2In2 * SymmetryFactor;
 	      TmpA2A2A2Out[j1] = TmpA2A2A2Out2 * SymmetryFactor;
 	      TmpA3A3A3In[j1] =  TmpA3A3A3In2 * SymmetryFactor;
 	      TmpA3A3A3Out[j1] = TmpA3A3A3Out2 * SymmetryFactor;
 	      TmpA4A4A4In[j1] =  TmpA4A4A4In2 * SymmetryFactor;
 	      TmpA4A4A4Out[j1] = TmpA4A4A4Out2 * SymmetryFactor;
 	      TmpA5A5A5In[j1] =  TmpA5A5A5In2 * SymmetryFactor;
 	      TmpA5A5A5Out[j1] = TmpA5A5A5Out2 * SymmetryFactor;
 	      TmpA6A6A6In[j1] =  TmpA6A6A6In2 * SymmetryFactor;
 	      TmpA6A6A6Out[j1] = TmpA6A6A6Out2 * SymmetryFactor;

	    }

	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  this->NBodyInteractionFactors[i][Index] = 2.0 * FactorU * ((TmpA1A1A1In[j1] * TmpA1A1A1Out[j2]) 
									     + (TmpA2A2A2In[j1] * TmpA2A2A2Out[j2])
									     + (TmpA3A3A3In[j1] * TmpA3A3A3Out[j2])
									     + (TmpA4A4A4In[j1] * TmpA4A4A4Out[j2])
									     + (TmpA5A5A5In[j1] * TmpA5A5A5Out[j2])
									     + (TmpA6A6A6In[j1] * TmpA6A6A6Out[j2]));
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}


      delete[] TmpA1A1A1In;
      delete[] TmpA1A1A1Out;
      delete[] TmpA2A2A2In;
      delete[] TmpA2A2A2Out;
      delete[] TmpA3A3A3In;
      delete[] TmpA3A3A3Out;
      delete[] TmpA4A4A4In;
      delete[] TmpA4A4A4Out;
      delete[] TmpA5A5A5In;
      delete[] TmpA5A5A5Out;
      delete[] TmpA6A6A6In;
      delete[] TmpA6A6A6Out;
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

Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2)
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

Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
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

Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementBAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
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

Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
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

Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
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

Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteAAB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
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

Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  Complex Tmp = 1.0;
  Tmp += Phase (((double) (kx4 - kx1)) * this->KxFactor);
  Tmp += Phase (((double) (ky4 - ky1)) * this->KyFactor);
  Tmp += Phase ((((double) (kx4 - kx1)) * this->KxFactor) + (((double) (ky4 - ky1)) * this->KyFactor));
  Tmp *= Phase (0.5 * ((((double) (kx6 + kx5 - kx3 - kx2)) * this->KxFactor)
		       + ((((double) (ky6 + ky5 - ky3 - ky2)) * this->KyFactor))));
  return Tmp;
}


