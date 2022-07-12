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
//                                                                            //
//                        last modification : 03/04/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;



// default constructor
//

ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian::ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = strength of the repulsive on-site interaction
// vPotential = strength of the repulsive two body neareast neighbor interaction
// flatBandFlag = use flat band model
// flatBandOneBodyGap = set the gap between the first band and the second band when using the flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian::ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
														 int nbrSiteY, Abstract2DTightBindingModel* tightBindingModel,
														 double uPotential, double vPotential, 
														 bool flatBandFlag, double flatBandOneBodyGap, 
														 AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;
  this->FlatBandOneBodyGap = flatBandOneBodyGap;
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;
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

ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian::~ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian::EvaluateInteractionFactors()
{
  double FactorU = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
  double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
  
  int NbrSites = this->NbrSiteX * this->NbrSiteY;
    
  long TotalNbrInteractionFactors = 0;
  
  int NbrInternalIndices = 2;
  this->InteractionFactorsSigma = new Complex***** [NbrInternalIndices];
  for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
    {
      this->InteractionFactorsSigma[sigma3] = new Complex****  [NbrInternalIndices];
      for (int sigma4 = sigma3; sigma4 < NbrInternalIndices; ++sigma4)
	{
	  this->InteractionFactorsSigma[sigma3][sigma4] = new Complex***[NbrInternalIndices];
	  for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	    {
	      this->InteractionFactorsSigma[sigma3][sigma4][sigma1] = new Complex**[NbrInternalIndices];
	    }
	}
    }

  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  
  this->InteractionFactorsupupupup = 0;
  this->InteractionFactorsupupupdown = 0;
  this->InteractionFactorsupupdowndown = 0;
  this->InteractionFactorsdowndownupup = 0;
  this->InteractionFactorsupdownupdown = 0;
  this->InteractionFactorsupdownupup = 0;
  this->InteractionFactorsupdowndowndown = 0;
  this->InteractionFactorsdowndownupdown = 0;
  this->InteractionFactorsdowndowndowndown = 0;
  
  
  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
  if ((this->FlatBand == false) || (this->FlatBandOneBodyGap != 0.0))
    {
      this->OneBodyInteractionFactorsupup = new double [NbrSites];
      this->OneBodyInteractionFactorsdowndown = new double [NbrSites];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	  if (this->FlatBand == false)
	    {
	      this->OneBodyInteractionFactorsupup[Index] = this->TightBindingModel->GetEnergy(0, Index);
	      this->OneBodyInteractionFactorsdowndown[Index] = this->TightBindingModel->GetEnergy(1, Index);
	    }
	  else
	    {
	      if (this->FlatBandOneBodyGap != 0.0)
		{
		  this->OneBodyInteractionFactorsupup[Index] = 0.0;
		  this->OneBodyInteractionFactorsdowndown[Index] = this->FlatBandOneBodyGap;		
		}
	    }
	  OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
	}
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)      
	  ++this->NbrInterSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndex(((kx1 + kx2) % this->NbrSiteX), ((ky1 + ky2) % this->NbrSiteY))];    
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
	    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex(((kx1 + kx2) % this->NbrSiteX), ((ky1 + ky2) % this->NbrSiteY));
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	  }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 < Index2)
		  ++this->NbrIntraSectorIndicesPerSum[(((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)];    
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
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 < Index2)
		  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int kx1 = this->IntraSectorIndicesPerSum[i][j1 << 1] / this->NbrSiteY;
	      int ky1 = this->IntraSectorIndicesPerSum[i][j1 << 1] % this->NbrSiteY;
	      int kx2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] / this->NbrSiteY;
	      int ky2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int kx3 = this->IntraSectorIndicesPerSum[i][j2 << 1] / this->NbrSiteY;
		  int ky3 = this->IntraSectorIndicesPerSum[i][j2 << 1] % this->NbrSiteY;
		  int kx4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] / this->NbrSiteY;
		  int ky4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] % this->NbrSiteY;
		  this->InteractionFactorsupup[i][Index] = 0.0;
		  this->InteractionFactorsdowndown[i][Index] = this->InteractionFactorsupup[i][Index];
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}
      this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int kx1 = this->InterSectorIndicesPerSum[i][j1 << 1] / this->NbrSiteY;
	      int ky1 = this->InterSectorIndicesPerSum[i][j1 << 1] % this->NbrSiteY;
	      int kx2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] / this->NbrSiteY;
	      int ky2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int kx3 = this->InterSectorIndicesPerSum[i][j2 << 1] / this->NbrSiteY;
		  int ky3 = this->InterSectorIndicesPerSum[i][j2 << 1] % this->NbrSiteY;
		  int kx4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] / this->NbrSiteY;
		  int ky4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] % this->NbrSiteY;
// 		  this->InteractionFactorsupdown[i][Index] = -2.0 * this->UPotential * (this->ComputeTwoBodyMatrixElementUpDown(kx2, ky2, kx4, ky4));// + this->ComputeTwoBodyMatrixElementUpDown(kx1, ky1, kx3, ky3));
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
    else //bosonic statistics
    {
      this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 <= Index2)
		  ++this->NbrIntraSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndex(((kx1 + kx2) % this->NbrSiteX), ((ky1 + ky2) % this->NbrSiteY))];    
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
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 <= Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex(((kx1 + kx2) % this->NbrSiteX), ((ky1 + ky2) % this->NbrSiteY));
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      
      
      Complex* TmpInteractionFactor;
      int* TmpIndices;
      int* TmpIndices2;
      for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	{
	  for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
	    {
	      this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] = new Complex*[this->NbrIntraSectorSums];
	      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		{
		  this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j] = new Complex [this->NbrIntraSectorIndicesPerSum[j] * this->NbrIntraSectorIndicesPerSum[j]];
		  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		  TmpIndices = this->IntraSectorIndicesPerSum[j];
		  for (int i1 = 0; i1 < Lim; i1 += 2)
		    {
		      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j][(i1 * Lim) >> 2]);
		      int Index1 = TmpIndices[i1];
		      int Index2 = TmpIndices[i1 + 1];
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  Complex Tmp = 0.0;
			  int Index3 = TmpIndices[i2];
			  int Index4 = TmpIndices[i2 + 1];
			  for (int i = 0; i < 2; ++i)
			    {
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, sigma3, sigma3, sigma1, sigma1, i, i, i, i);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, sigma3, sigma3, sigma1, sigma1, i, i, i, i);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, sigma3, sigma3, sigma1, sigma1, i, i, i, i);
			      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, sigma3, sigma3, sigma1, sigma1, i, i, i, i);
			    }
			  Tmp *= FactorU;
			  if (Index1 == Index2)
			    Tmp *= 0.5;
			  if (Index3 == Index4)
			    Tmp *= 0.5;
			  (*TmpInteractionFactor) = Tmp;
			  ++TmpInteractionFactor;
			}
		    }
		}
	    }
	  for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
	    {
	      for (int sigma4 = sigma3 + 1; sigma4 < NbrInternalIndices; ++sigma4)
		{
		  this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] = new Complex*[this->NbrIntraSectorSums];
		  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		    {
		      this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j] = new Complex [this->NbrIntraSectorIndicesPerSum[j] * this->NbrInterSectorIndicesPerSum[j]];
		      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		      TmpIndices = this->IntraSectorIndicesPerSum[j];
		      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
		      TmpIndices2 = this->InterSectorIndicesPerSum[j];
		      for (int i1 = 0; i1 < Lim; i1 += 2)
			{
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j][(i1 * Lim2) >> 2]);
			  int Index1 = TmpIndices[i1];
			  int Index2 = TmpIndices[i1 + 1];
			  for (int i2 = 0; i2 < Lim2; i2 += 2)
			    {
			      Complex Tmp = 0.0;
			      int Index3 = TmpIndices2[i2];
			      int Index4 = TmpIndices2[i2 + 1];
			      for (int i = 0; i < 2; ++i)
				{
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, sigma3, sigma4, sigma1, sigma1, i, i, i, i);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, sigma4, sigma3, sigma1, sigma1, i, i, i, i);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, sigma3, sigma4, sigma1, sigma1, i, i, i, i);
				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, sigma4, sigma3, sigma1, sigma1, i, i, i, i);
				}
			      if (Index1 == Index2)
				Tmp *= 0.5;
			      Tmp *= FactorU;
			      (*TmpInteractionFactor) = Tmp;
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	    }			  
	}
      for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	{
	  for (int sigma2 = sigma1 + 1; sigma2 < NbrInternalIndices; ++sigma2)
	    {
	      for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
		{
		  this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
		  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		    {
		      this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j] = new Complex [this->NbrInterSectorIndicesPerSum[j] * this->NbrIntraSectorIndicesPerSum[j]];
		      int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
		      TmpIndices = this->IntraSectorIndicesPerSum[j];
		      int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
		      TmpIndices2 = this->InterSectorIndicesPerSum[j];
		      for (int i1 = 0; i1 < Lim2; i1 += 2)
			{
			  TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j][(i1 * Lim) >> 2]);
			  int Index1 = TmpIndices2[i1];
			  int Index2 = TmpIndices2[i1 + 1];
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      Complex Tmp = 0.0;
			      int Index3 = TmpIndices[i2];
			      int Index4 = TmpIndices[i2 + 1];
			      for (int i = 0; i < 2; ++i)
				{
 				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, sigma3, sigma3, sigma1, sigma2, i, i, i, i);
 				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, sigma3, sigma3, sigma1, sigma2, i, i, i, i);
 				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, sigma3, sigma3, sigma2, sigma1, i, i, i, i);
 				  Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, sigma3, sigma3, sigma2, sigma1, i, i, i, i);
				}
			      if (Index3 == Index4)
				Tmp *= 0.5;
			      Tmp *= FactorU;
			      (*TmpInteractionFactor) = Tmp;
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	      for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
		{
		  for (int sigma4 = sigma3 + 1; sigma4 < NbrInternalIndices; ++sigma4)
		    {
		      this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
		      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
			{
			  this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j] = new Complex [this->NbrInterSectorIndicesPerSum[j] * this->NbrInterSectorIndicesPerSum[j]];
			  int Lim = 2 * this->NbrIntraSectorIndicesPerSum[j];
			  TmpIndices = this->IntraSectorIndicesPerSum[j];
			  int Lim2 = 2 * this->NbrInterSectorIndicesPerSum[j];
			  TmpIndices2 = this->InterSectorIndicesPerSum[j];
			  for (int i1 = 0; i1 < Lim2; i1 += 2)
			    {
			      TmpInteractionFactor = &(this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j][(i1 * Lim2) >> 2]);
			      int Index1 = TmpIndices2[i1];
			      int Index2 = TmpIndices2[i1 + 1];
			      for (int i2 = 0; i2 < Lim2; i2 += 2)
				{
				  Complex Tmp = 0.0;
				  int Index3 = TmpIndices2[i2];
				  int Index4 = TmpIndices2[i2 + 1];
				  for (int i = 0; i < 2; ++i)
				    {
 				      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, sigma3, sigma4, sigma1, sigma2, i, i, i, i);
 				      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, sigma4, sigma3, sigma1, sigma2, i, i, i, i);
 				      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, sigma3, sigma4, sigma2, sigma1, i, i, i, i);
 				      Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, sigma4, sigma3, sigma2, sigma1, i, i, i, i);
				    }
				  Tmp *= FactorU;
				  (*TmpInteractionFactor) = Tmp;
				  ++TmpInteractionFactor;
				}
			    }
			}
		    }
		}
	    }
	}
    }
  cout << this->InteractionFactorsSigma[0][0][0][1][0][0] << endl;
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

Complex ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian::ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp;
  double FactorX = 2.0 * M_PI / ((double) this->NbrSiteX);
  double FactorY = 2.0 * M_PI / ((double) this->NbrSiteY);
  
  Tmp = 0.5*Phase(0.5 * (FactorX * (double) (kx1 - kx2) + FactorY * (double) (ky1 - ky2)));
  Tmp += 0.5*Phase(0.5 * (FactorX * (double) (kx1 - kx2) - FactorY * (double) (ky1 - ky2)));
  Tmp += 0.5*Phase(0.5 * (-FactorX * (double) (kx1 - kx2) + FactorY * (double) (ky1 - ky2)));
  Tmp += 0.5*Phase(0.5 * (-FactorX * (double) (kx1 - kx2) - FactorY * (double) (ky1 - ky2)));
  
  return Tmp;
}
