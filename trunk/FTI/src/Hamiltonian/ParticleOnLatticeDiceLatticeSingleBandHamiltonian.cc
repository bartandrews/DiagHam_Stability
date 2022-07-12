////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//      class of checkerboard lattice model with interacting particles        //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 08/09/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeDiceLatticeSingleBandHamiltonian.h" //CHANGED HERE TO DICE FROM KAGOME
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
ParticleOnLatticeDiceLatticeSingleBandHamiltonian::ParticleOnLatticeDiceLatticeSingleBandHamiltonian() //CHANGED HERE TO DICE FROM KAGOME
{
  this->BandIndex = 0;
  this->U3Potential = 0.0;
  this->U6Potential = 0.0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellX = number of sites in the x direction
// nbrCellY = number of sites in the y direction
// u3Potential = strength of the repulsive two body neareast neighbor interaction (or on-site density-density potential for bosons)
// u6Potential = strength of the repulsive two body second nearest neighbor interaction (or nearest neighbor density-density potential for bosons)
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeDiceLatticeSingleBandHamiltonian::ParticleOnLatticeDiceLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, 
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

ParticleOnLatticeDiceLatticeSingleBandHamiltonian::~ParticleOnLatticeDiceLatticeSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeDiceLatticeSingleBandHamiltonian::EvaluateInteractionFactors()
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
      for (int k1a = 0; k1a < this->NbrSiteX; ++k1a)
	for (int k2a = 0; k2a < this->NbrSiteX; ++k2a)
	  for (int k1b = 0; k1b < this->NbrSiteY; ++k1b)
	    for (int k2b = 0; k2b < this->NbrSiteY; ++k2b) 
	      {
		int Index1 = (k1a * this->NbrSiteY) + k1b;
		int Index2 = (k2a * this->NbrSiteY) + k2b;
		if (Index1 < Index2)
		  ++this->NbrSectorIndicesPerSum[(((k1a + k2a) % this->NbrSiteX) *  this->NbrSiteY) + ((k1b + k2b) % this->NbrSiteY)];    
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
      for (int k1a = 0; k1a < this->NbrSiteX; ++k1a)
	for (int k2a = 0; k2a < this->NbrSiteX; ++k2a)
	  for (int k1b = 0; k1b < this->NbrSiteY; ++k1b)
	    for (int k2b = 0; k2b < this->NbrSiteY; ++k2b) 
	      {
		int Index1 = (k1a * this->NbrSiteY) + k1b;
		int Index2 = (k2a * this->NbrSiteY) + k2b;
		if (Index1 < Index2)
		  {
		    int TmpSum = (((k1a + k2a) % this->NbrSiteX) *  this->NbrSiteY) + ((k1b + k2b) % this->NbrSiteY);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }
      /*double FactorUAB = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUAC = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUBC = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY)); */   //I COMMENTED OUT THIS BIT BECAUSE WE DONT USE U & V POTENTIALS

      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int k1a = Index1 / this->NbrSiteY;
	      int k1b = Index1 % this->NbrSiteY;
	      int k2a = Index2 / this->NbrSiteY;
	      int k2b = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int k3a = Index3 / this->NbrSiteY;
		  int k3b = Index3 % this->NbrSiteY;
		  int k4a = Index4 / this->NbrSiteY;
		  int k4b = Index4 % this->NbrSiteY;


//MAYBE COMMENT OUT THIS SECTION AS WELL?	

/*	  
		  this->InteractionFactors[i][Index] = FactorUAB * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(k2a, k2b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUAB * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(k1a, k1b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUAB * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(k2a, k2b, k3a, k3b);
 		  this->InteractionFactors[i][Index] += FactorUAB * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(k1a, k1b, k3a, k3b);

 		  this->InteractionFactors[i][Index] += FactorUAC * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(k2a, k2b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUAC * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(k1a, k1b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUAC * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(k2a, k2b, k3a, k3b);
 		  this->InteractionFactors[i][Index] += FactorUAC * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(k1a, k1b, k3a, k3b);

 		  this->InteractionFactors[i][Index] += FactorUBC * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(k1a, k1b, k2a, k2b, k3a, k3b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUBC * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(k2a, k2b, k1a, k1b, k3a, k3b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUBC * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(k1a, k1b, k2a, k2b, k4a, k4b, k3a, k3b);
 		  this->InteractionFactors[i][Index] += FactorUBC * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(k2a, k2b, k1a, k1b, k4a, k4b, k3a, k3b); */

		  this->InteractionFactors[i][Index] *= -2.0;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  else  // code for bosons
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
      double FactorU3 = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if (this->FlatBand == false) //I WROTE U3Potential it was VPotential
	FactorU3 *= this->U3Potential;
      double FactorU6 = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if (this->FlatBand == false) //I WROTE U6Potential it was VPotential
	FactorU6 *= this->U6Potential;
      else
        if ( this->U3Potential != 0.0)
          FactorU6 *= this->U6Potential/this->U3Potential;
        // take U6Potential as the unit of energy if U3=0 and flat-band option is chosen.
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

// this one six-fold
 		  this->InteractionFactors[i][Index] = FactorU6 * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
 		  this->InteractionFactors[i][Index] += FactorU6 * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
 		  this->InteractionFactors[i][Index] += FactorU6 * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
 		  this->InteractionFactors[i][Index] += FactorU6 * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();

 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

// from here i add other parts which are the same but different numbers

                  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index1][BandIndex][3]) * OneBodyBasis[Index3][BandIndex][3] * Conj(OneBodyBasis[Index2][BandIndex][3]) * OneBodyBasis[Index4][BandIndex][3]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index2][BandIndex][3]) * OneBodyBasis[Index3][BandIndex][3] * Conj(OneBodyBasis[Index1][BandIndex][3]) * OneBodyBasis[Index4][BandIndex][3]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index1][BandIndex][3]) * OneBodyBasis[Index4][BandIndex][3] * Conj(OneBodyBasis[Index2][BandIndex][3]) * OneBodyBasis[Index3][BandIndex][3]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index2][BandIndex][3]) * OneBodyBasis[Index4][BandIndex][3] * Conj(OneBodyBasis[Index1][BandIndex][3]) * OneBodyBasis[Index3][BandIndex][3]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);


//this one six-fold
                  this->InteractionFactors[i][Index] += FactorU6 * (Conj(OneBodyBasis[Index1][BandIndex][4]) * OneBodyBasis[Index3][BandIndex][4] * Conj(OneBodyBasis[Index2][BandIndex][4]) * OneBodyBasis[Index4][BandIndex][4]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] += FactorU6 * (Conj(OneBodyBasis[Index2][BandIndex][4]) * OneBodyBasis[Index3][BandIndex][4] * Conj(OneBodyBasis[Index1][BandIndex][4]) * OneBodyBasis[Index4][BandIndex][4]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] += FactorU6 * (Conj(OneBodyBasis[Index1][BandIndex][4]) * OneBodyBasis[Index4][BandIndex][4] * Conj(OneBodyBasis[Index2][BandIndex][4]) * OneBodyBasis[Index3][BandIndex][4]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactors[i][Index] += FactorU6 * (Conj(OneBodyBasis[Index2][BandIndex][4]) * OneBodyBasis[Index4][BandIndex][4] * Conj(OneBodyBasis[Index1][BandIndex][4]) * OneBodyBasis[Index3][BandIndex][4]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);



                  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index1][BandIndex][5]) * OneBodyBasis[Index3][BandIndex][5] * Conj(OneBodyBasis[Index2][BandIndex][5]) * OneBodyBasis[Index4][BandIndex][5]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index2][BandIndex][5]) * OneBodyBasis[Index3][BandIndex][5] * Conj(OneBodyBasis[Index1][BandIndex][5]) * OneBodyBasis[Index4][BandIndex][5]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index1][BandIndex][5]) * OneBodyBasis[Index4][BandIndex][5] * Conj(OneBodyBasis[Index2][BandIndex][5]) * OneBodyBasis[Index3][BandIndex][5]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactors[i][Index] += FactorU3 * (Conj(OneBodyBasis[Index2][BandIndex][5]) * OneBodyBasis[Index4][BandIndex][5] * Conj(OneBodyBasis[Index1][BandIndex][5]) * OneBodyBasis[Index3][BandIndex][5]) * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);


/*
// we don't need this if statement?		  
		  if (this->VPotential != 0.0)
		    {
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);
		      
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(kx2, ky2, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(kx1, ky1, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(kx2, ky2, kx3, ky3);
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(kx1, ky1, kx3, ky3);
		      
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		      this->InteractionFactors[i][Index] += FactorV * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
		    }

		  if (this-> != 0.0)
		    {
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementABNNN(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementABNNN(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementABNNN(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementABNNN(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
		      
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementACNNN(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementACNNN(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementACNNN(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementACNNN(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
		      
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBCNNN(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBCNNN(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBCNNN(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		      this->InteractionFactors[i][Index] += FactorW * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBCNNN(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
		    }
*/
		  
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

Complex ParticleOnLatticeDiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementAB(int k1a, int k1b, int k2a, int k2b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->TightBindingModel->GetProjectedMomentum(k2a, k2b, 0) - this->TightBindingModel->GetProjectedMomentum(k1a, k1b, 0)));
//   Complex Tmp = 2.0 * cos (0.5 * (this->KxFactor * ((double) (k2a - k1a))));
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

Complex ParticleOnLatticeDiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementAC(int k1a, int k1b, int k2a, int k2b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->TightBindingModel->GetProjectedMomentum(k2a, k2b, 1) - this->TightBindingModel->GetProjectedMomentum(k1a, k1b, 1)));
//   Complex Tmp = 2.0 * cos (0.5 * (this->KyFactor * ((double) (k2b - k1b))));
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

Complex ParticleOnLatticeDiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementBC(int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->TightBindingModel->GetProjectedMomentum(k3a, k3b, 0) - this->TightBindingModel->GetProjectedMomentum(k1a, k1b, 0) + this->TightBindingModel->GetProjectedMomentum(k4a, k4b, 1) - this->TightBindingModel->GetProjectedMomentum(k2a, k2b, 1)));
//   Complex Tmp = 2.0 * cos (0.5 * ((this->KxFactor * ((double) (k3a - k1a))) + (this->KyFactor * ((double) (k4b - k2b)))));
  //Complex Tmp = Phase(0.5 * ((this->KxFactor * ((double) (k3a - k1a))) + (this->KyFactor * ((double) (k4b - k2b)))));
  return Tmp;
}


// compute the matrix element for the two body interaction between two sites A and B in the next nearest neighbor interaction
//
// k1a = creation momentum along x for the A site
// k1b = creation momentum along y for the A site
// k2a = creation momentum along x for the B site
// k2b = creation momentum along y for the B site
// k3a = annihilation momentum along x for the A site
// k3b = annihilation momentum along y for the A site
// k4a = annihilation momentum along x for the B site
// k4b = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementABNNN(int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->TightBindingModel->GetProjectedMomentum(k4a, k4b, 0) - this->TightBindingModel->GetProjectedMomentum(k2a, k2b, 0))
			   - (this->TightBindingModel->GetProjectedMomentum(k4a, k4b, 1) - this->TightBindingModel->GetProjectedMomentum(k2a, k2b, 1)));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and C in the next nearest neighbor interaction
//
// k1a = creation momentum along x for the A site
// k1b = creation momentum along y for the A site
// k2a = creation momentum along x for the C site
// k2b = creation momentum along y for the C site
// k3a = annihilation momentum along x for the A site
// k3b = annihilation momentum along y for the A site
// k4a = annihilation momentum along x for the C site
// k4b = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementACNNN(int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b)
{
  Complex Tmp = 2.0 * cos ((this->TightBindingModel->GetProjectedMomentum(k4a, k4b, 0) - this->TightBindingModel->GetProjectedMomentum(k2a, k2b, 0))
			   - 0.5 * (this->TightBindingModel->GetProjectedMomentum(k4a, k4b, 1) - this->TightBindingModel->GetProjectedMomentum(k2a, k2b, 1)));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B and C in the next nearest neighbor interaction
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

Complex ParticleOnLatticeDiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementBCNNN(int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b)
{
//   Complex Tmp = 2.0 * cos (0.5 * (this->TightBindingModel->GetProjectedMomentum(k4a, k4b, 0) - this->TightBindingModel->GetProjectedMomentum(k2a, k2b, 0)
// 				  + this->TightBindingModel->GetProjectedMomentum(k3a, k3b, 0) - this->TightBindingModel->GetProjectedMomentum(k1a, k1b, 0))
// 			   + 0.5 * (this->TightBindingModel->GetProjectedMomentum(k4a, k4b, 1) - this->TightBindingModel->GetProjectedMomentum(k2a, k2b, 1)
// 				    + this->TightBindingModel->GetProjectedMomentum(k3a, k3b, 1) - this->TightBindingModel->GetProjectedMomentum(k1a, k1b, 1)));
  Complex Tmp = 2.0 * cos ((this->TightBindingModel->GetProjectedMomentum(k4a, k4b, 0) - this->TightBindingModel->GetProjectedMomentum(k2a, k2b, 0))
			   + 0.5 * (this->TightBindingModel->GetProjectedMomentum(k4a, k4b, 1) - this->TightBindingModel->GetProjectedMomentum(k2a, k2b, 1))
			   + 0.5 * (this->TightBindingModel->GetProjectedMomentum(k3a, k3b, 0) - this->TightBindingModel->GetProjectedMomentum(k1a, k1b, 0)));
  return Tmp;
}


// compute the matrix element for on-site two body interaction involving A sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeDiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteAA()
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

Complex ParticleOnLatticeDiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return Phase(0.5 * (this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0) + this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 0) 
		      - this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0) - this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 0) ));
//   return Phase(0.5 * this->KxFactor * ((double) (kx4 + kx3 - kx2 -kx1)));
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

Complex ParticleOnLatticeDiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return Phase(0.5 * (this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1) + this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 1) - this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1) - this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 1) ));
//   return Phase(0.5 * this->KyFactor * ((double) (ky4 + ky3 - ky2 -ky1)));
}

