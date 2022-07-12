////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 2d topological insulator based on the simple TI        //
//                       model and restricted to two bands                    //
//                                                                            //
//                        last modification : 27/09/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeDiceLatticeTwoBandHamiltonian.h"
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

ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ParticleOnLatticeDiceLatticeTwoBandHamiltonian() 
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = repulsive on-site potential strength between different orbitals
// vPotential = repulsive on-site potential strength between opposite spins
// mass = mass term of the simple TI model
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)


//MY HEADER ETC 10/02/2018

ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ParticleOnLatticeDiceLatticeTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrCellX, int nbrCellY, 
													 double u3Potential, double u6Potential,  
													 Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
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
  this->U3Potential = u3Potential;
  this->U6Potential = u6Potential;
  //this->WPotential = wPotential;
  this->BandIndex = 0;
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

ParticleOnLatticeDiceLatticeTwoBandHamiltonian::~ParticleOnLatticeDiceLatticeTwoBandHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeDiceLatticeTwoBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactorsupup = new double [this->NbrSiteX * this->NbrSiteY];
      this->OneBodyInteractionFactorsdowndown = new double [this->NbrSiteX * this->NbrSiteY];
    }

  this->ComputeOneBodyMatrices(OneBodyBasis);
 
  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    this->NbrIntraSectorIndicesPerSum[i] = 0;      

  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)      
	  ++this->NbrInterSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY))];
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      if (this->NbrInterSectorIndicesPerSum[i] > 0)
	{
	  this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
	  this->NbrInterSectorIndicesPerSum[i] = 0;
	}
    }
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)    
	  {
	    int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY));
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = ((kx1 * this->NbrSiteY) + ky1);
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = ((kx2 * this->NbrSiteY) + ky2);
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	  }
  
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 < Index2)
		  ++this->NbrIntraSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY))];
	      }
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY));
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      
      double Factor = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorAUpADown = Factor * this->VPotential;
      double FactorBUpBDown = Factor * this->VPotential;
      double FactorAUpBDown = Factor * this->WPotential;
      double FactorADownBUp = Factor * this->WPotential;
      if (this->FlatBand == false)
	Factor *= this->UPotential;
      double FactorAUpBUp = Factor;
      double FactorADownBDown = Factor;

      Complex Tmp;

      //  upup upup coefficient
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;


		// In original file we had 4 numbers. here wel have 8. (see .h for details: first 4 energy index, other 4 - site index)
		
		//Maybe like this, changing site index only: (all zeros, all 1, all 2 etc)


                  this->InteractionFactorsupupupup[i][Index] = ComputeOnSiteContributions( OneBodyBasis, FactorU3, FactorU6,
							Index1,  Index2,  Index3, Index4, 
							0, 0, 0, 0); //!!! added OneBodyBasis
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  upup downdown coefficient
      this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
                  this->InteractionFactorsupupdowndown[i][Index] = 0.0;

		
//Fermions
		
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 0, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 0, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupupdowndown[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  downdown downdown coefficient
      this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndowndowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
                  this->InteractionFactorsdowndowndowndown[i][Index] = 0.0;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndowndowndown[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  downdown upup coefficient
      this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndownupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
                  this->InteractionFactorsdowndownupup[i][Index] = 0.0;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsdowndownupup[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  updown upup coefficient
      this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupup[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupdownupup[i][Index] = 0.0;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 0, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdownupup[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      //  updown downdown coefficient
      this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupdowndowndown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupdowndowndown[i][Index] = 0.0;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
                  this->InteractionFactorsupdowndowndown[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  upup updown coefficient
      this->InteractionFactorsupupupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupupupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupupupdown[i][Index] = 0.0;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 1, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 1, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupupupdown[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  downdown updown coefficient
      this->InteractionFactorsdowndownupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsdowndownupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsdowndownupdown[i][Index] = 0.0;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsdowndownupdown[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}


      //  updown updown coefficient
      this->InteractionFactorsupdownupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
                  this->InteractionFactorsupdownupdown[i][Index] = 0.0;
		  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 0, 0, 0, 0) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorAUpADown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                   Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 1, 1, 1, 1) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorBUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 2, 2, 2, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorAUpBUp * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 3, 3, 3, 3) * this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorADownBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                   Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                   Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 4, 4, 4, 4) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorAUpBDown * Tmp;
                  Tmp  = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 0, 1, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 0, 0, 1, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 0, 5, 5, 5, 5) * this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupdownupdown[i][Index] += -2.0 * FactorADownBUp * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      // bosonic interaction
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 <= Index2)
		  ++this->NbrIntraSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY))];    
	      }
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 <= Index2)
		  {
		    int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY));
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      
      double Factor = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorAUpADown = Factor * this->VPotential;
      double FactorBUpBDown = Factor * this->VPotential;
      if ((this->FlatBand == false) || (this->VPotential != 0.0))
	Factor *= this->UPotential;
      double FactorAUpAUp = Factor;
      double FactorBUpBUp = Factor;
      double FactorADownADown = Factor;
      double FactorBDownBDown = Factor;
      double FactorAUpBUp = Factor;
      double FactorADownBDown = Factor;
      double FactorAUpBDown = Factor;
      double FactorADownBUp = Factor;

      Complex Tmp;

      //  upup upup coefficient
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupupupup[i][Index] = ComputeOnSiteContributions( OneBodyBasis, FactorU3, FactorU6,
							Index1,  Index2,  Index3, Index4, 
							0, 0, 0, 0); //B!!! 
		}
	    }
	}

      //  upup downdown coefficient
      this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupupdowndown[i][Index] = ComputeOnSiteContributions( OneBodyBasis, FactorU3, FactorU6,
							Index1,  Index2,  Index3, Index4, 
							0, 0, 1, 1); //B!!!

                  
		}
	    }
	}

      //  downdown downdown coefficient
      this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndowndowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsdowndowndowndown[i][Index] = ComputeOnSiteContributions( OneBodyBasis, FactorU3, FactorU6,
							Index1,  Index2,  Index3, Index4, 
							1, 1, 1, 1);//B!!!

                 
		}
	    }
	}

      //  downdown upup coefficient
      this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndownupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsdowndownupup[i][Index] = ComputeOnSiteContributions( OneBodyBasis, FactorU3, FactorU6,
							Index1,  Index2,  Index3, Index4, 
							1, 1, 0, 0);
		}
	    }
	}

      //  updown upup coefficient
      this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupup[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupdownupup[i][Index] = ComputeOnSiteContributions( OneBodyBasis, FactorU3, FactorU6,
							Index1,  Index2,  Index3, Index4, 
							0, 1, 0, 0); //B!!!

                 
		}
	    }
	}
      //  updown downdown coefficient
      this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupdowndowndown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupdowndowndown[i][Index] = ComputeOnSiteContributions( OneBodyBasis, FactorU3, FactorU6,
							Index1,  Index2,  Index3, Index4, 
							0, 1, 1, 1);//B!!!

		
		}
	    }
	}

      //  upup updown coefficient
      this->InteractionFactorsupupupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupupupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupupupdown[i][Index] = ComputeOnSiteContributions( OneBodyBasis, FactorU3, FactorU6,
							Index1,  Index2,  Index3, Index4, 
							0, 0, 0, 1);//B!!!

		
		}
	    }
	}

      //  downdown updown coefficient
      this->InteractionFactorsdowndownupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsdowndownupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsdowndownupdown[i][Index] = ComputeOnSiteContributions( OneBodyBasis, FactorU3, FactorU6,
							Index1,  Index2,  Index3, Index4, 
							1, 1, 0, 1);//B!!!

		 
		}
	    }
	}


      //  updown updown coefficient
      this->InteractionFactorsupdownupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
                  this->InteractionFactorsupdownupdown[i][Index] = ComputeOnSiteContributions( OneBodyBasis, FactorU3, FactorU6,
							Index1,  Index2,  Index3, Index4, 
							0, 1, 0, 1);//B!!!

		  
		}
	    }
	}
    }

  delete[] OneBodyBasis;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the contributions to the interaction factor for on-site interactions (for bosons only)
// 
  // momentumIndex1 = compact momentum index of the first creation operator
  // momentumIndex2 = compact momentum index of the second creation operator
  // momentumIndex3 = compact momentum index of the first annihilation operator
  // momentumIndex4 = compact momentum index of the second annihiliation operator
  // energyIndex1 = energy index of the first creation operator
  // energyIndex2 = energy index of the second creation operator
  // energyIndex3 = energy index of the first annihilation operator
  // energyIndex4 = energy index of the second annihiliation operator
Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeOnSiteContributions(ComplexMatrix *OneBodyBasis, double FactorU3, double FactorU6,
							int Index1, int Index2, int Index3, int Index4, 
							int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4)





//Bosons
  




{
		  
  		  Complex Tmp, Result(0.0);

		//I've deleted minus before the factors of to be in correspondance with single band verison
  		  
		  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, energyIndex1, energyIndex2, energyIndex3, energyIndex4, 0, 0, 0, 0);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, energyIndex1, energyIndex2, energyIndex4, energyIndex3, 0, 0, 0, 0);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, energyIndex2, energyIndex1, energyIndex3, energyIndex4, 0, 0, 0, 0);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, energyIndex2, energyIndex1, energyIndex4, energyIndex3, 0, 0, 0, 0);
                  Result += 2.0 * FactorU6 * Tmp; //u6
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, energyIndex1, energyIndex2, energyIndex3, energyIndex4, 1, 1, 1, 1);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, energyIndex1, energyIndex2, energyIndex4, energyIndex3, 1, 1, 1, 1);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, energyIndex2, energyIndex1, energyIndex3, energyIndex4, 1, 1, 1, 1);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, energyIndex2, energyIndex1, energyIndex4, energyIndex3, 1, 1, 1, 1); 
                  Result += 2.0 * FactorU3 * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, energyIndex1, energyIndex2, energyIndex3, energyIndex4, 2, 2, 2, 2);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, energyIndex1, energyIndex2, energyIndex4, energyIndex3, 2, 2, 2, 2);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, energyIndex2, energyIndex1, energyIndex3, energyIndex4, 2, 2, 2, 2);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, energyIndex2, energyIndex1, energyIndex4, energyIndex3, 2, 2, 2, 2); 
                  Result += 2.0 * FactorU3 * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, energyIndex1, energyIndex2, energyIndex3, energyIndex4, 3, 3, 3, 3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, energyIndex1, energyIndex2, energyIndex4, energyIndex3, 3, 3, 3, 3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, energyIndex2, energyIndex1, energyIndex3, energyIndex4, 3, 3, 3, 3);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, energyIndex2, energyIndex1, energyIndex4, energyIndex3, 3, 3, 3, 3); 

                  Result += 2.0 * FactorU3 * Tmp;
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, energyIndex1, energyIndex2, energyIndex3, energyIndex4, 4, 4, 4, 4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, energyIndex1, energyIndex2, energyIndex4, energyIndex3, 4, 4, 4, 4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, energyIndex2, energyIndex1, energyIndex3, energyIndex4, 4, 4, 4, 4);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, energyIndex2, energyIndex1, energyIndex4, energyIndex3, 4, 4, 4, 4); 
                  Result += 2.0 * FactorU6 * Tmp; //u6
                  Tmp = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, energyIndex1, energyIndex2, energyIndex3, energyIndex4, 5, 5, 5, 5);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, energyIndex1, energyIndex2, energyIndex4, energyIndex3, 5, 5, 5, 5);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, energyIndex2, energyIndex1, energyIndex3, energyIndex4, 5, 5, 5, 5);
                  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, energyIndex2, energyIndex1, energyIndex4, energyIndex3, 5, 5, 5, 5);
                  Result += 2.0 * FactorU3 * Tmp;


		


   return Result;
}

    


// compute the matrix element for the two body interaction between two sites A with up spins
//
// kx1 = momentum along x for the creation operator on first A site with spin up
// ky1 = momentum along y for the creation operator on first A site with spin up
// kx2 = momentum along x for the creation operator on second A site with spin up
// ky2 = momentum along y for the creation operator on second A site with spin up
// kx3 = momentum along x for the annihilation operator on first A site with spin up
// ky3 = momentum along y for the annihilation operator on first A site with spin up
// kx4 = momentum along x for the annihilation operator on second A site with spin up
// ky4 = momentum along y for the annihilation operator on second A site with spin up
// return value = corresponding matrix element  

Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeTwoBodyMatrixElementAUpAUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0 ;
  return Tmp;
}
 
// compute the matrix element for the two body interaction between two sites A with down spins
//
// kx1 = momentum along x for the creation operator on first A site with spin down
// ky1 = momentum along y for the creation operator on first A site with spin down
// kx2 = momentum along x for the creation operator on second A site with spin down
// ky2 = momentum along y for the creation operator on second A site with spin down
// kx3 = momentum along x for the annihilation operator on first A site with spin down
// ky3 = momentum along y for the annihilation operator on first A site with spin down
// kx4 = momentum along x for the annihilation operator on second A site with spin down
// ky4 = momentum along y for the annihilation operator on second A site with spin down
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeTwoBodyMatrixElementADownADown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0 ;
  return Tmp;
}
 
// compute the matrix element for the two body interaction between two sites B with up spins
//
// kx1 = momentum along x for the creation operator on first B site with spin up
// ky1 = momentum along y for the creation operator on first B site with spin up
// kx2 = momentum along x for the creation operator on second B site with spin up
// ky2 = momentum along y for the creation operator on second B site with spin up
// kx3 = momentum along x for the annihilation operator on first B site with spin up
// ky3 = momentum along y for the annihilation operator on first B site with spin up
// kx4 = momentum along x for the annihilation operator on second B site with spin up
// ky4 = momentum along y for the annihilation operator on second B site with spin up
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeTwoBodyMatrixElementBUpBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0 ;
  return Tmp;
}
 
// compute the matrix element for the two body interaction between two sites B with down spins
//
// kx1 = momentum along x for the creation operator on first B site with spin down
// ky1 = momentum along y for the creation operator on first B site with spin down
// kx2 = momentum along x for the creation operator on second B site with spin down
// ky2 = momentum along y for the creation operator on second B site with spin down
// kx3 = momentum along x for the annihilation operator on first B site with spin down
// ky3 = momentum along y for the annihilation operator on first B site with spin down
// kx4 = momentum along x for the annihilation operator on second B site with spin down
// ky4 = momentum along y for the annihilation operator on second B site with spin down
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeTwoBodyMatrixElementBDownBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0 ;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and B with up spins
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kz2 = momentum along z for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// kz4 = momentum along z for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0 ;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and B with down spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kz1 = momentum along z for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kz2 = momentum along z for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kz3 = momentum along z for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// kz4 = momentum along z for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
}

// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kz1 = momentum along z for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kz2 = momentum along z for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kz3 = momentum along z for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// kz4 = momentum along z for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeTwoBodyMatrixElementADownBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
}

// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kz2 = momentum along z for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// kz4 = momentum along z for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
}

// compute the matrix element for the two body interaction between two sites A with opposite spins 
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on A site with spin down
// ky2 = momentum along y for the creation operator on A site with spin down
// kz2 = momentum along z for the creation operator on A site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on A site with spin down
// ky4 = momentum along y for the annihilation operator on A site with spin down
// kz4 = momentum along z for the annihilation operator on A site with spin down
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B with opposite spins 
//
// kx1 = momentum along x for the creation operator on B site with spin up
// ky1 = momentum along y for the creation operator on B site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on B site with spin up
// ky3 = momentum along y for the annihilation operator on B site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0;
  return Tmp;
}

// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices
//ff !!!
void ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	HermitianMatrix TmpOneBodyHamiltonian(4, true);
	int Index = ((kx * this->NbrSiteY) + ky);
	Complex d2 (sin (((double) ky) * this->KyFactor), 0.0);
	double d1 = sin (((double) kx) * this->KxFactor);
	double d3 = (this->Mass - cos (((double) kx) * this->KxFactor) - cos (((double) ky) * this->KyFactor));
	TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d3);
	TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -d3);
	TmpOneBodyHamiltonian.SetMatrixElement(2, 2, d3);
	TmpOneBodyHamiltonian.SetMatrixElement(3, 3, -d3);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 1, d1);
	TmpOneBodyHamiltonian.SetMatrixElement(2, 3, d1);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 3, d2);
	TmpOneBodyHamiltonian.SetMatrixElement(1, 2, d2);
	
	ComplexMatrix TmpMatrix(4, 4, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	TmpMatrix[2][2] = 1.0;
	TmpMatrix[3][3] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	oneBodyBasis[Index] = TmpMatrix;	
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsupup[Index] = TmpDiag(0, 0);
	    this->OneBodyInteractionFactorsdowndown[Index] = TmpDiag(1, 1);
	  }
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << " " << TmpDiag(3, 3) << endl;
      }
}
