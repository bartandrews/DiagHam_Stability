///////////////////////////////////////////////////////////////////////////////
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
#include "Hamiltonian/ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian.h"
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

ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// threeBodyPotential = strength of the repulsive three body neareast neighbor interaction
// uPotential = strength of the repulsive two body neareast neighbor interaction
// vPotential = strength of the repulsive two body neareast neighbor interaction
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, 
															   int nbrCellY, double threeBodyPotential, double uPotential, double vPotential ,  Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
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

  this->NBodyValue = 3;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->FlatBand = flatBandFlag;
  this->ThreeBodyPotential = threeBodyPotential;
  this->UPotential = uPotential;
  this->VPotential = vPotential;

  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;

  if (fabs(uPotential)<1e-15)
    this->TwoBodyFlag = false;
  else
    this->TwoBodyFlag = true;
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

ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::~ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::EvaluateInteractionFactors()
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
	    {
	      this->NbrSectorIndicesPerSum[i] = 0;
	      delete [] this->SectorIndicesPerSum[i];
	    }
	  this->InteractionFactors = new Complex* [this->NbrSectorSums];
	}
      else // this->TwoBodyFlag == true
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
	  this->InteractionFactors = new Complex* [this->NbrSectorSums];
      
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
		      
                      Complex sumU = 0.0;
                      Complex sumV = 0.0;

		      sumU = (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] *
			      Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
			* this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		      sumU -= (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] *
			       Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
			* this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);
		      sumU -= (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] *
			       Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1])
			* this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		      sumU += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] *
			       Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1])
			* this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		      
		      sumU += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] *
			       Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2])
			* this->ComputeTwoBodyMatrixElementAC(kx2, ky2, kx4, ky4);
		      sumU -= (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] *
			       Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2])
			* this->ComputeTwoBodyMatrixElementAC(kx1, ky1, kx4, ky4);
		      sumU -= (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] *
			       Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2])
			* this->ComputeTwoBodyMatrixElementAC(kx2, ky2, kx3, ky3);
		      sumU += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] *
			       Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2])
			* this->ComputeTwoBodyMatrixElementAC(kx1, ky1, kx3, ky3);
		      
		      sumU += (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1] *
			       Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2])
			* this->ComputeTwoBodyMatrixElementBC(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		      sumU -= (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1] *
			       Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2])
			* this->ComputeTwoBodyMatrixElementBC(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		      sumU -= (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] *
			       Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2])
			* this->ComputeTwoBodyMatrixElementBC(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		      sumU += (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] *
			       Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2])
			* this->ComputeTwoBodyMatrixElementBC(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
		      
                      this->InteractionFactors[i][Index] = -2.0 * (FactorU * sumU);// + FactorV * sumV);
                      TotalNbrInteractionFactors++;
                      ++Index;
		    }
		}
	    }
	
	  cout << "nbr 2-body interaction = " << TotalNbrInteractionFactors << endl;
	  TotalNbrInteractionFactors=0;
	  
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

      // double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      double ThreeBodyFactor = this->ThreeBodyPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];

      Complex* TmpABCUpIn = new Complex[TmpLargestSector];
      Complex* TmpABCUpOut = new Complex[TmpLargestSector];
      Complex* TmpABCDownIn = new Complex[TmpLargestSector];
      Complex* TmpABCDownOut = new Complex[TmpLargestSector];

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

    	      Complex TmpABCUpIn2 = 0.0;
	      Complex TmpABCUpOut2 = 0.0;
    	      Complex TmpABCDownIn2 = 0.0;
	      Complex TmpABCDownOut2 = 0.0;

	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];
 		  TmpABCUpIn2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2]) * this->ComputeThreeBodyMatrixElementABCUpIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
		  TmpABCDownIn2 += PermutationSign[l1] * Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2]) * this->ComputeThreeBodyMatrixElementABCDownIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);

		  TmpABCUpOut2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2] * this->ComputeThreeBodyMatrixElementABCUpOut(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
		  TmpABCDownOut2 += PermutationSign[l1] * OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2] * this->ComputeThreeBodyMatrixElementABCDownOut(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
		}

	      TmpABCUpIn[j1] = TmpABCUpIn2;
	      TmpABCUpOut[j1] = TmpABCUpOut2;
    	      TmpABCDownIn[j1] = TmpABCDownIn2;
	      TmpABCDownOut[j1] = TmpABCDownOut2;
	    }
	  
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  this->NBodyInteractionFactors[i][Index] = 2.0 * ThreeBodyFactor * ((TmpABCUpIn[j1] * TmpABCUpOut[j2])
										      + (TmpABCDownIn[j1] * TmpABCDownOut[j2]));
		  // cout << j1 << " "<<j2<<" "<< 2.0 * ThreeBodyFactor * ((TmpABCUpIn[j1] * TmpABCUpOut[j2]) + (TmpABCDownIn[j1] * TmpABCDownOut[j2])) << " " << TmpABCUpIn[j1] <<" "<< TmpABCUpOut[j2] <<" " << TmpABCDownIn[j1] <<" "<< TmpABCDownOut[j2] << endl;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
      delete[] TmpABCUpIn;
      delete[] TmpABCUpOut;
      delete[] TmpABCDownIn;
      delete[] TmpABCDownOut;
      
      for (int i=0; i<NbrPermutations; ++i)
	delete [] Permutations[i];
      delete [] Permutations;
      delete [] PermutationSign;
    }
  else // bosonic statistics
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
	    {
	      this->NbrSectorIndicesPerSum[i] = 0;
	      delete [] this->SectorIndicesPerSum[i];
	    }
	  this->InteractionFactors = new Complex* [this->NbrSectorSums];
	}
      else
	{
	  double FactorU = this->UPotential*0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
	  double FactorV = this->VPotential*0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
	  
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

		      Complex sumU=0.0, sumV=0.0;
		      sumU = (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		      sumU += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		      sumU += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		      sumU += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();

		      sumU += (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		      sumU += (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		      sumU += (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		      sumU += (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

		      sumU += (Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		      sumU += (Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		      sumU += (Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		      sumU += (Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

		      sumV = (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] *
			      Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
			* this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		      sumV += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] *
			       Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
			* this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);
		      sumV += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] *
			       Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1])
			* this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		      sumV += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] *
			       Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1])
			* this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		      
		      sumV += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] *
			       Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2])
			* this->ComputeTwoBodyMatrixElementAC(kx2, ky2, kx4, ky4);
		      sumV += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] *
			       Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2])
			* this->ComputeTwoBodyMatrixElementAC(kx1, ky1, kx4, ky4);
		      sumV += (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] *
			       Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2])
			* this->ComputeTwoBodyMatrixElementAC(kx2, ky2, kx3, ky3);
		      sumV += (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] *
			       Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2])
			* this->ComputeTwoBodyMatrixElementAC(kx1, ky1, kx3, ky3);
		      
		      sumV += (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1] *
			       Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2])
			* this->ComputeTwoBodyMatrixElementBC(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		      sumV += (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1] *
			       Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2])
			* this->ComputeTwoBodyMatrixElementBC(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		      sumV += (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] *
			       Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2])
			* this->ComputeTwoBodyMatrixElementBC(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		      sumV += (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] *
			       Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2])
			* this->ComputeTwoBodyMatrixElementBC(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);


		      if (Index3 == Index4)
			{
			  sumU *= 0.5;
			  sumV *= 0.5;
			}
		      if (Index1 == Index2)
			{
			  sumU *= 0.5;
			  sumV *= 0.5;
			}

		      this->InteractionFactors[i][Index] = 2.0 * FactorU * sumU + 2.0 * FactorV * sumV;

		      TotalNbrInteractionFactors++;
		      ++Index;

		    }
		}
	    }
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

      Complex* TmpAAAIn = new Complex[TmpLargestSector];
      Complex* TmpAAAOut = new Complex[TmpLargestSector];
      Complex* TmpBBBIn = new Complex[TmpLargestSector];
      Complex* TmpBBBOut = new Complex[TmpLargestSector];
      Complex* TmpCCCIn = new Complex[TmpLargestSector];
      Complex* TmpCCCOut = new Complex[TmpLargestSector];
      Complex* TmpABCUpIn = new Complex[TmpLargestSector];
      Complex* TmpABCUpOut = new Complex[TmpLargestSector];
      Complex* TmpABCDownIn = new Complex[TmpLargestSector];
      Complex* TmpABCDownOut = new Complex[TmpLargestSector];


      // double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      double ThreeBodyFactor = 1.0 * this->ThreeBodyPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      double VThreeBodyFactor = 0.0 * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);

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

	      Complex TmpAAAIn2 = 0.0;
	      Complex TmpAAAOut2 = 0.0;
	      Complex TmpBBBIn2 = 0.0;
	      Complex TmpBBBOut2 = 0.0;
	      Complex TmpCCCIn2 = 0.0;
	      Complex TmpCCCOut2 = 0.0;
    	      Complex TmpABCUpIn2 = 0.0;
	      Complex TmpABCUpOut2 = 0.0;
    	      Complex TmpABCDownIn2 = 0.0;
	      Complex TmpABCDownOut2 = 0.0;

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

 		  TmpAAAIn2 += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0]) * this->ComputeThreeBodyMatrixElementOnSiteAAAIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpBBBIn2 += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1]) * this->ComputeThreeBodyMatrixElementOnSiteBBBIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpCCCIn2 += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2]) * this->ComputeThreeBodyMatrixElementOnSiteCCCIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
		  
 		  TmpAAAOut2 += OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0] * this->ComputeThreeBodyMatrixElementOnSiteAAAOut(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpBBBOut2 += OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1] * this->ComputeThreeBodyMatrixElementOnSiteBBBOut(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
 		  TmpCCCOut2 += OneBodyBasis[IndexIn[TmpPerm[0]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][2] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2] * this->ComputeThreeBodyMatrixElementOnSiteCCCOut(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);

 		  TmpABCUpIn2 += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2]) * this->ComputeThreeBodyMatrixElementABCUpIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
		  TmpABCDownIn2 += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2]) * this->ComputeThreeBodyMatrixElementABCDownIn(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);

		  TmpABCUpOut2 += OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2] * this->ComputeThreeBodyMatrixElementABCUpOut(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);
		  TmpABCDownOut2 += OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][2] * this->ComputeThreeBodyMatrixElementABCDownOut(KxIn[TmpPerm[0]], KyIn[TmpPerm[0]], KxIn[TmpPerm[1]], KyIn[TmpPerm[1]], KxIn[TmpPerm[2]], KyIn[TmpPerm[2]]);

 		}

 	      TmpAAAIn[j1] =  TmpAAAIn2 * SymmetryFactor;
 	      TmpAAAOut[j1] = TmpAAAOut2 * SymmetryFactor;
 	      TmpBBBIn[j1] =  TmpBBBIn2 * SymmetryFactor;
 	      TmpBBBOut[j1] = TmpBBBOut2 * SymmetryFactor;
 	      TmpCCCIn[j1] =  TmpCCCIn2 * SymmetryFactor;
 	      TmpCCCOut[j1] = TmpCCCOut2 * SymmetryFactor;
	      TmpABCUpIn[j1] = TmpABCUpIn2 * SymmetryFactor;
	      TmpABCUpOut[j1] = TmpABCUpOut2 * SymmetryFactor;
    	      TmpABCDownIn[j1] = TmpABCDownIn2 * SymmetryFactor;
	      TmpABCDownOut[j1] = TmpABCDownOut2 * SymmetryFactor;

	    }

	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  this->NBodyInteractionFactors[i][Index] = 2.0 * ThreeBodyFactor * ((TmpAAAIn[j1] * TmpAAAOut[j2])
										     + (TmpBBBIn[j1] * TmpBBBOut[j2])
										     + (TmpCCCIn[j1] * TmpCCCOut[j2]));
		  this->NBodyInteractionFactors[i][Index] += 2.0 * VThreeBodyFactor * ((TmpABCUpIn[j1] * TmpABCUpOut[j2])
										       + (TmpABCDownIn[j1] * TmpABCDownOut[j2]));
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }

	}
      for (int i=0; i<NbrPermutations; ++i)
	delete [] Permutations[i];
      delete [] Permutations;
      delete [] PermutationSign;


      delete[] TmpAAAIn;
      delete[] TmpAAAOut;
      delete[] TmpBBBIn;
      delete[] TmpBBBOut;
      delete[] TmpCCCIn;
      delete[] TmpCCCOut;
      delete[] TmpABCUpIn;
      delete[] TmpABCUpOut;
      delete[] TmpABCDownIn;
      delete[] TmpABCDownOut;
    }

  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;

  delete [] OneBodyBasis;
}


// conventions adopted for matrix elements:
// modified with respect to Tang, Mei, Wen:
// unit cell is triangle standing on base. bottom left corner site A and bottom right corner B, tip is site C.


// compute the matrix element for the two body interaction between two sites A and B 
//
// k1a = creation momentum along x for the B site
// k1b = creation momentum along y for the B site
// k2a = annihilation momentum along x for the B site
// k2b = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAB(int k1a, int k1b, int k2a, int k2b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KxFactor * ((double) (k2a - k1a))));
  //Complex Tmp = Phase (0.5 * (this->KxFactor * ((double) (k2a - k1a))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and C 
//
// k1a = creation momentum along x for the C site
// k1b = creation momentum along y for the C site
// k2a = annihilation momentum along x for the C site
// k2b = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAC(int k1a, int k1b, int k2a, int k2b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KyFactor * ((double) (k2b - k1b))));
  //Complex Tmp = Phase (0.5 * (this->KyFactor * ((double) (k2b - k1b))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B and C 
//
// k1a = creation momentum along x for the B site
// k1b = creation momentum along y for the B site
// k2a = creation momentum along x for the C site
// k2b = creation momentum along y for the C site
// k3a = annihilation momentum along x for the B site
// k3b = annihilation momentum along y for the B site
// k4a = annihilation momentum along x for the C site
// k4b = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementBC(int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b)
{
  Complex Tmp = 2.0 * cos (0.5 * ((this->KxFactor * ((double) (k3a - k1a))) + (this->KyFactor * ((double) (k4b - k2b)))));
  //Complex Tmp = Phase(0.5 * ((this->KxFactor * ((double) (k3a - k1a))) + (this->KyFactor * ((double) (k4b - k2b)))));
  return Tmp;
}


// compute the matrix element for on-site two body interaction involving A sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementOnSiteAA()
{
  return 1.0;
}

// compute the matrix element for on-site two body interaction involving B sites
//
// kx1 = first creation momentum along x for the B site
// ky1 = first creation momentum along y for the B site
// kx2 = second creation momentum along x for the B site
// ky2 = second creation momentum along y for the B site
// kx3 = first annihilation momentum along x for the B site
// ky3 = first annihilation momentum along y for the B site
// kx4 = second annihilation momentum along x for the B site
// ky4 = second annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return Phase(0.5 * this->KxFactor * ((double) (kx4 + kx3 - kx2 -kx1)));
}

// compute the matrix element for on-site two body interaction involving C sites
//
// kx1 = first creation momentum along x for the C site
// ky1 = first creation momentum along y for the C site
// kx2 = second creation momentum along x for the C site
// ky2 = second creation momentum along y for the C site
// kx3 = first annihilation momentum along x for the C site
// ky3 = first annihilation momentum along y for the C site
// kx4 = second annihilation momentum along x for the C site
// ky4 = second annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementOnSiteCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return Phase(0.5 * this->KyFactor * ((double) (ky4 + ky3 - ky2 -ky1)));
}

/* END TWO-BODY TERMS */



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

Complex ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the on-site three body interaction related to sites B
//
// kx1 = first creation momentum along x for the first B site
// ky1 = first creation momentum along y for the first B site
// kx2 = second creation momentum along x for the second B site
// ky2 = second creation momentum along y for the second B site
// kx3 = third creation momentum along x for the second B site
// ky3 = third creation momentum along y for the second B site
// kx4 = first annihilation momentum along x for the first B site
// ky4 = first annihilation momentum along y for the first B site
// kx5 = second annihilation momentum along x for the second B site
// ky5 = second annihilation momentum along y for the second B site
// kx6 = third annihilation momentum along x for the second B site
// ky6 = third annihilation momentum along y for the second B site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  Complex Tmp = Phase (0.5 * ((((double) (kx6 + kx5 + kx4 - kx3 - kx2 - kx1)) * this->KxFactor)));
  return Tmp;
}

// compute the matrix element for the on-site three body interaction related to sites A
//
// kx1 = first creation momentum along x for the first C site
// ky1 = first creation momentum along y for the first C site
// kx2 = second creation momentum along x for the second C site
// ky2 = second creation momentum along y for the second C site
// kx3 = third creation momentum along x for the second C site
// ky3 = third creation momentum along y for the second C site
// kx4 = first annihilation momentum along x for the first C site
// ky4 = first annihilation momentum along y for the first C site
// kx5 = second annihilation momentum along x for the second C site
// ky5 = second annihilation momentum along y for the second C site
// kx6 = third annihilation momentum along x for the second C site
// ky6 = third annihilation momentum along y for the second C site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteCCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  Complex Tmp = Phase (0.5 * ((((double) (ky6 + ky5 + ky4 - ky3 - ky2 - ky1)) * this->KyFactor)));
  return Tmp;
}



