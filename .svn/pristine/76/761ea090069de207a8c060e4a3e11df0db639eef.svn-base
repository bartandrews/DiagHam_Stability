////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 3d topological insulator based on the pryochlore       //
//                       model and restricted to four bands                   //
//                                                                            //
//                        last modification : 13/08/2012                      //
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
#include "Hamiltonian/ParticleOnCubicLatticeFourBandPyrochloreHamiltonian.h"
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

ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::ParticleOnCubicLatticeFourBandPyrochloreHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// uPotential = strength of the repulsive two body on site interactions
// vPotential = trength of the repulsive two body neareast neighbor interaction
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::ParticleOnCubicLatticeFourBandPyrochloreHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int nbrSiteX, 
													 int nbrSiteY, int nbrSiteZ, double uPotential, double vPotential, double wuPotential, double wvPotential, 
													 Abstract3DTightBindingModel* tightBindingModel, 
													 bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteYZ = this->NbrSiteY * this->NbrSiteZ;
  this->LzMax = nbrSiteX * nbrSiteY * nbrSiteZ - 1;
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;

  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->WUPotential = wuPotential;
  this->WVPotential = wvPotential;

  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsupum = 0;
  this->OneBodyInteractionFactorsupdp = 0;
  this->OneBodyInteractionFactorsupdm = 0;
  this->OneBodyInteractionFactorsumum = 0;
  this->OneBodyInteractionFactorsumdp = 0;
  this->OneBodyInteractionFactorsumdm = 0;
  this->OneBodyInteractionFactorsdpdp = 0;
  this->OneBodyInteractionFactorsdpdm = 0;
  this->OneBodyInteractionFactorsdmdm = 0;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();

  int Dim = this->Particles->GetHilbertSpaceDimension();
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

ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::~ParticleOnCubicLatticeFourBandPyrochloreHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  int NbrSites = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;

  int NbrInternalIndices = 4;
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
  this->InteractionFactorsupumupum = 0;
  this->InteractionFactorsupdpupdp = 0;
  this->InteractionFactorsupdmupdm = 0;
  this->InteractionFactorsumumumum = 0;
  this->InteractionFactorsumdpumdp = 0;
  this->InteractionFactorsumdmumdm = 0;
  this->InteractionFactorsdpdpdpdp = 0;
  this->InteractionFactorsdpdmdpdm = 0;
  this->InteractionFactorsdmdmdmdm = 0;

  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsumum = 0;
  this->OneBodyInteractionFactorsdpdp = 0;
  this->OneBodyInteractionFactorsdmdm = 0;
  this->OneBodyInteractionFactorsupum = 0;
  this->OneBodyInteractionFactorsupdp = 0;
  this->OneBodyInteractionFactorsupdm = 0;
  this->OneBodyInteractionFactorsumdp = 0;
  this->OneBodyInteractionFactorsumdm = 0;
  this->OneBodyInteractionFactorsdpdm = 0;

  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactorsupup = new double [NbrSites];
      this->OneBodyInteractionFactorsumum = new double [NbrSites];
      this->OneBodyInteractionFactorsdpdp = new double [NbrSites];
      this->OneBodyInteractionFactorsdmdm = new double [NbrSites];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      for (int kz = 0; kz < this->NbrSiteZ; ++kz)
	{
	  int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky, kz);
	  if (this->FlatBand == false)
	    {
	      this->OneBodyInteractionFactorsupup[Index] = 0.5 * this->TightBindingModel->GetEnergy(0, Index);
	      this->OneBodyInteractionFactorsumum[Index] = 0.5 * this->TightBindingModel->GetEnergy(1, Index);
	      this->OneBodyInteractionFactorsdpdp[Index] = 0.5 * this->TightBindingModel->GetEnergy(2, Index);
	      this->OneBodyInteractionFactorsdmdm[Index] = 0.5 * this->TightBindingModel->GetEnergy(3, Index);
	    }
	  OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
	}


  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
  this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    this->NbrIntraSectorIndicesPerSum[i] = 0;      

  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)      
	  for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
	    for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
	      {
		int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
										 (ky1 + ky2) % this->NbrSiteY,
										 (kz1 + kz2) % this->NbrSiteZ);
		++this->NbrInterSectorIndicesPerSum[TmpSum];    
	      }
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
	  for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
	    for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
	      {
		int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
										 (ky1 + ky2) % this->NbrSiteY,
										 (kz1 + kz2) % this->NbrSiteZ);
		this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1, kz1);
		this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2, kz2);
		++this->NbrInterSectorIndicesPerSum[TmpSum];    
	      }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
		for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
		  {
		    int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1, kz1);
		    int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2, kz2);
		    if (Index1 <= Index2)
		      {
			int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
											 (ky1 + ky2) % this->NbrSiteY,
											 (kz1 + kz2) % this->NbrSiteZ);
			++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		      }
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
	      for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
		for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
		  {
		    int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1, kz1);
		    int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2, kz2);
		    if (Index1 <= Index2)
		      {
			int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
											 (ky1 + ky2) % this->NbrSiteY,
											 (kz1 + kz2) % this->NbrSiteZ);
			this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
			this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
			++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		      }
		  }
      
      double Factor = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ))*2.0*this->VPotential;
      double FactorWU = 1.0 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ))*this->WUPotential;
      double FactorWV = 1.0 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ))*this->WVPotential;
      Complex* TmpInteractionFactor;
      int* TmpIndices;
      int* TmpIndices2;
      for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	{
	  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
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
			  int Index3 = TmpIndices[i2];
			  int Index4 = TmpIndices[i2 + 1];
			  (*TmpInteractionFactor) += this->ComputeOnSiteContributionOppositeSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma1, sigma3, sigma3, Factor, -1);
			  (*TmpInteractionFactor) += this->ComputeNearestNeighborInteractionContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma1, sigma3, sigma3, FactorWU, FactorWV, -1);
			  ++TmpInteractionFactor;
			}
		    }
		}
	    }
	  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
	    {
	      for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
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
			      int Index3 = TmpIndices2[i2];
			      int Index4 = TmpIndices2[i2 + 1];
			      (*TmpInteractionFactor) += this->ComputeOnSiteContributionOppositeSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma1, sigma3, sigma4, Factor, -1);
			      (*TmpInteractionFactor) += this->ComputeNearestNeighborInteractionContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma1, sigma3, sigma4, FactorWU, FactorWV, -1);
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	    }			  
	}
      for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	{
	  for (int sigma2 = sigma1 + 1; sigma2 < 4; ++sigma2)
	    {
	      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
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
			      int Index3 = TmpIndices[i2];
			      int Index4 = TmpIndices[i2 + 1];
			     (*TmpInteractionFactor) += this->ComputeOnSiteContributionOppositeSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma2, sigma3, sigma3, Factor, -1);
			     (*TmpInteractionFactor) += this->ComputeNearestNeighborInteractionContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma2, sigma3, sigma3, FactorWU, FactorWV, -1);
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
		{
		  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
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
				  int Index3 = TmpIndices2[i2];
				  int Index4 = TmpIndices2[i2 + 1];
				  (*TmpInteractionFactor) += this->ComputeOnSiteContributionOppositeSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma2, sigma3, sigma4, Factor, -1);
				  (*TmpInteractionFactor) += this->ComputeNearestNeighborInteractionContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma2, sigma3, sigma4, FactorWU, FactorWV, -1);
				  ++TmpInteractionFactor;
				}
			    }
			}
		    }
		}			  
	    }
	}
    }
  else
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
		for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
		  {
		    int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1, kz1);
		    int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2, kz2);
		    if (Index1 <= Index2)
		      {
			int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
											 (ky1 + ky2) % this->NbrSiteY,
											 (kz1 + kz2) % this->NbrSiteZ);
			++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		      }
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
	      for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
		for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
		  {
		    int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1, kz1);
		    int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2, kz2);
		    if (Index1 <= Index2)
		      {
			int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
											 (ky1 + ky2) % this->NbrSiteY,
											 (kz1 + kz2) % this->NbrSiteZ);
			this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
			this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
			++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		      }
		  }
      
      double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ)) * this->UPotential;
      double FactorV = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ)) * 2.0 * this->VPotential;
      double FactorWU = 1.0 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ)) * this->WUPotential;
      double FactorWV = 1.0 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ)) * this->WVPotential;
      Complex* TmpInteractionFactor;
      int* TmpIndices;
      int* TmpIndices2;
      for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	{
	  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
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
			  int Index3 = TmpIndices[i2];
			  int Index4 = TmpIndices[i2 + 1];
			  
			  (*TmpInteractionFactor) = this->ComputeOnSiteContributionSameSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma1, sigma3, sigma3, FactorU);
			  (*TmpInteractionFactor) += this->ComputeOnSiteContributionOppositeSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma1, sigma3, sigma3, FactorV, 1);
			  (*TmpInteractionFactor) += this->ComputeNearestNeighborInteractionContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma1, sigma3, sigma3, FactorWU, FactorWV, 1);
			  
			  ++TmpInteractionFactor;
			}
		    }
		}
	    }
	  for (int sigma3 = 0; sigma3 < 4; ++sigma3)
	    {
	      for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
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
			      int Index3 = TmpIndices2[i2];
			      int Index4 = TmpIndices2[i2 + 1];
			      (*TmpInteractionFactor) = this->ComputeOnSiteContributionSameSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma1, sigma3, sigma4, FactorU);
			      (*TmpInteractionFactor) += this->ComputeOnSiteContributionOppositeSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma1, sigma3, sigma4, FactorV, 1);
			      (*TmpInteractionFactor) += this->ComputeNearestNeighborInteractionContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma1, sigma3, sigma4, FactorWU, FactorWV, 1);
			      
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	    }			  
	}
      for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	{
	  for (int sigma2 = sigma1 + 1; sigma2 < 4; ++sigma2)
	    {
	      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
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
			      int Index3 = TmpIndices[i2];
			      int Index4 = TmpIndices[i2 + 1];
			      (*TmpInteractionFactor) = this->ComputeOnSiteContributionSameSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma2, sigma3, sigma3, FactorU);
			      (*TmpInteractionFactor) += this->ComputeOnSiteContributionOppositeSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma2, sigma3, sigma3, FactorV, 1);
			      (*TmpInteractionFactor) += this->ComputeNearestNeighborInteractionContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma2, sigma3, sigma3, FactorWU, FactorWV, 1);
			      ++TmpInteractionFactor;
			    }
			}
		    }
		}
	      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
		{
		  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
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
				  int Index3 = TmpIndices2[i2];
				  int Index4 = TmpIndices2[i2 + 1];
				  (*TmpInteractionFactor) = this->ComputeOnSiteContributionSameSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma2, sigma3, sigma4, FactorU);
				  (*TmpInteractionFactor) += this->ComputeOnSiteContributionOppositeSpin(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma2, sigma3, sigma4, FactorV, 1);
				  (*TmpInteractionFactor) += this->ComputeNearestNeighborInteractionContribution(OneBodyBasis, Index1, Index2, Index3, Index4, sigma1, sigma2, sigma3, sigma4, FactorWU, FactorWV, 1);
				  ++TmpInteractionFactor;
				}
			    }
			}
		    }
		}			  
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}


// compute the contribution to the onsite interaction matrix elements with the same spin
// 
// oneBodyBasis = array of transformation basis matrices
// momentumIndex1 = compact momentum index of the first creation operator
// momentumIndex2 = compact momentum index of the second creation operator
// momentumIndex3 = compact momentum index of the first annihilation operator
// momentumIndex4 = compact momentum index of the second annihiliation operator
// energyIndex1 = energy index of the first creation operator
// energyIndex2 = energy index of the second creation operator
// energyIndex3 = energy index of the first annihilation operator
// energyIndex4 = energy index of the second annihiliation operator
// factorU = repulsive on-site potential strength between identical spins
//return value = corresponding matrix element

Complex ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::ComputeOnSiteContributionSameSpin(ComplexMatrix* oneBodyBasis,
						       int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
						       int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,double factorU)
{
  Complex Tmp = 0;
  for (int i = 0; i < 8; ++i)
  {
    Tmp += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, i, i, i, i);
    Tmp += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, i, i, i, i);
    Tmp += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, i, i, i, i);
    Tmp += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, i, i, i, i);
   }
  Tmp *= factorU;
  if ((momentumIndex1 == momentumIndex2) && (energyIndex1 == energyIndex2))
    Tmp *= 0.5;
  if ((momentumIndex3 == momentumIndex4) && (energyIndex3 == energyIndex4))
    Tmp *= 0.5;
  return Tmp;
}

// compute the contribution to the onsite interaction matrix elements with the opposite spin
// 
// oneBodyBasis = array of transformation basis matrices
// momentumIndex1 = compact momentum index of the first creation operator
// momentumIndex2 = compact momentum index of the second creation operator
// momentumIndex3 = compact momentum index of the first annihilation operator
// momentumIndex4 = compact momentum index of the second annihiliation operator
// energyIndex1 = energy index of the first creation operator
// energyIndex2 = energy index of the second creation operator
// energyIndex3 = energy index of the first annihilation operator
// energyIndex4 = energy index of the second annihiliation operator
// factorV = repulsive on-site potential strength between identical spins
// epsilon = sign of the permutation +1 for bosons, -1 for fermions
//return value = corresponding matrix element

Complex ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::ComputeOnSiteContributionOppositeSpin(ComplexMatrix* oneBodyBasis,
						       int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
						       int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,double factorV, int epsilon)
{
  Complex Tmp = 0.0;
  for (int i = 0; i < 4; ++i)
  {
      Tmp += epsilon*this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, i, i + 4, i, i + 4);
      Tmp += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, i, i + 4, i, i + 4);
      Tmp += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, i, i + 4, i, i + 4);
      Tmp += epsilon*this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, i, i + 4, i, i + 4);
    }
  Tmp *= factorV;
  if (epsilon == 1)
  {
    if ((momentumIndex1 == momentumIndex2) && (energyIndex1 == energyIndex2))
      Tmp *= 0.5;
    if ((momentumIndex3 == momentumIndex4) && (energyIndex3 == energyIndex4))
      Tmp *= 0.5;
  }
  return Tmp;  
}

// compute the contribution to the nearest neighbor interaction 
// 
// oneBodyBasis = array of transformation basis matrices
// momentumIndex1 = compact momentum index of the first creation operator
// momentumIndex2 = compact momentum index of the second creation operator
// momentumIndex3 = compact momentum index of the first annihilation operator
// momentumIndex4 = compact momentum index of the second annihiliation operator
// energyIndex1 = energy index of the first creation operator
// energyIndex2 = energy index of the second creation operator
// energyIndex3 = energy index of the first annihilation operator
// energyIndex4 = energy index of the second annihiliation operator
// factorWU = repulsive nearest neighbor potential strength between identical spins
// factorWV = repulsive nearest neighbor potential strength between opposite spins
// epsilon = sign of the permutation +1 for bosons, -1 for fermions
//return value = corresponding matrix element
Complex ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::ComputeNearestNeighborInteractionContribution(ComplexMatrix* oneBodyBasis,
						       int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
						       int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,double factorWU, double factorWV, int epsilon)
{
 double TmpKx1;
 double TmpKy1; 
 double TmpKz1;
 double TmpKx2;
 double TmpKy2; 
 double TmpKz2;
 double TmpKx3;
 double TmpKy3; 
 double TmpKz3;
 double TmpKx4;
 double TmpKy4;
 double TmpKz4;
 int TmpKx; 
 int TmpKy;
 int TmpKz;
 double kxFactor = 2.0 * M_PI / this->NbrSiteX;
 double kyFactor = 2.0 * M_PI / this->NbrSiteY;
 double kzFactor = 2.0 * M_PI / this->NbrSiteZ;
 this->TightBindingModel->GetLinearizedMomentumIndex(momentumIndex1, TmpKx, TmpKy, TmpKz);
 TmpKx1 = TmpKx * kxFactor;
 TmpKy1 = TmpKy * kyFactor;
 TmpKz1 = TmpKz * kzFactor;
 this->TightBindingModel->GetLinearizedMomentumIndex(momentumIndex2, TmpKx, TmpKy, TmpKz);
 TmpKx2 = TmpKx * kxFactor;
 TmpKy2 = TmpKy * kyFactor;
 TmpKz2 = TmpKz * kzFactor;
 this->TightBindingModel->GetLinearizedMomentumIndex(momentumIndex3, TmpKx, TmpKy, TmpKz);
 TmpKx3 = TmpKx * kxFactor;
 TmpKy3 = TmpKy * kyFactor;
 TmpKz3 = TmpKz * kzFactor;
 this->TightBindingModel->GetLinearizedMomentumIndex(momentumIndex4, TmpKx, TmpKy, TmpKz);
 TmpKx4 = TmpKx * kxFactor;
 TmpKy4 = TmpKy * kyFactor;
 TmpKz4 = TmpKz * kzFactor;

 //Compute contribution to interaction factor of particles with identical spins
 Complex TmpWU = 0.0;
 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 0, 1, 0, 1)*(1.0 + Phase(TmpKx3 - TmpKx1 - (TmpKz3 - TmpKz1)))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 0, 1, 0, 1)*(1.0 + Phase(TmpKx4 - TmpKx1 - (TmpKz4 - TmpKz1)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 0, 1, 0, 1)*(1.0 + Phase(TmpKx3 - TmpKx2 - (TmpKz3 - TmpKz2)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 0, 1, 0, 1)*(1.0 + Phase(TmpKx4 - TmpKx2 - (TmpKz4 - TmpKz2)))*epsilon;
 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 4, 5, 4, 5)*(1.0 + Phase(TmpKx3 - TmpKx1 - (TmpKz3 - TmpKz1)))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 4, 5, 4, 5)*(1.0 + Phase(TmpKx4 - TmpKx1 - (TmpKz4 - TmpKz1)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 4, 5, 4, 5)*(1.0 + Phase(TmpKx3 - TmpKx2 - (TmpKz3 - TmpKz2)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 4, 5, 4, 5)*(1.0 + Phase(TmpKx4 - TmpKx2 - (TmpKz4 - TmpKz2)))*epsilon;
 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 0, 2, 0, 2)*(1.0 + Phase(TmpKx3 - TmpKx1))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 0, 2, 0, 2)*(1.0 + Phase(TmpKx4 - TmpKx1)); 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 0, 2, 0, 2)*(1.0 + Phase(TmpKx3 - TmpKx2));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 0, 2, 0, 2)*(1.0 + Phase(TmpKx4 - TmpKx2))*epsilon;
 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 4, 6, 4, 6)*(1.0 + Phase(TmpKx3 - TmpKx1))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 4, 6, 4, 6)*(1.0 + Phase(TmpKx4 - TmpKx1)); 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 4, 6, 4, 6)*(1.0 + Phase(TmpKx3 - TmpKx2));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 4, 6, 4, 6)*(1.0 + Phase(TmpKx4 - TmpKx2))*epsilon;

 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 0, 3, 0, 3)*(1.0 + Phase(TmpKx3 - TmpKx1 - (TmpKy3 - TmpKy1)))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 0, 3, 0, 3)*(1.0 + Phase(TmpKx4 - TmpKx1 - (TmpKy4 - TmpKy1)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 0, 3, 0, 3)*(1.0 + Phase(TmpKx3 - TmpKx2 - (TmpKy3 - TmpKy2)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 0, 3, 0, 3)*(1.0 + Phase(TmpKx4 - TmpKx2 - (TmpKy4 - TmpKy2)))*epsilon;
 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 4, 7, 4, 7)*(1.0 + Phase(TmpKx3 - TmpKx1 - (TmpKy3 - TmpKy1)))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 4, 7, 4, 7)*(1.0 + Phase(TmpKx4 - TmpKx1 - (TmpKy4 - TmpKy1)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 4, 7, 4, 7)*(1.0 + Phase(TmpKx3 - TmpKx2 - (TmpKy3 - TmpKy2)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 4, 7, 4, 7)*(1.0 + Phase(TmpKx4 - TmpKx2 - (TmpKy4 - TmpKy2)))*epsilon;

 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 1, 2, 1, 2)*(1.0 + Phase(TmpKz3 - TmpKz1))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 1, 2, 1, 2)*(1.0 + Phase(TmpKz4 - TmpKz1));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 1, 2, 1, 2)*(1.0 + Phase(TmpKz3 - TmpKz2));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 1, 2, 1, 2)*(1.0 + Phase(TmpKz4 - TmpKz2))*epsilon;
 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 5, 6, 5, 6)*(1.0 + Phase(TmpKz3 - TmpKz1))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 5, 6, 5, 6)*(1.0 + Phase(TmpKz4 - TmpKz1));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 5, 6, 5, 6)*(1.0 + Phase(TmpKz3 - TmpKz2));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 5, 6, 5, 6)*(1.0 + Phase(TmpKz4 - TmpKz2))*epsilon;
 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 1, 3, 1, 3)*(1.0 + Phase(TmpKz3 - TmpKz1 - (TmpKy3 - TmpKy1)))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 1, 3, 1, 3)*(1.0 + Phase(TmpKz4 - TmpKz1 - (TmpKy4 - TmpKy1)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 1, 3, 1, 3)*(1.0 + Phase(TmpKz3 - TmpKz2 - (TmpKy3 - TmpKy2)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 1, 3, 1, 3)*(1.0 + Phase(TmpKz4 - TmpKz2 - (TmpKy4 - TmpKy2)))*epsilon;
 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 5, 7, 5, 7)*(1.0 + Phase(TmpKz3 - TmpKz1 - (TmpKy3 - TmpKy1)))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 5, 7, 5, 7)*(1.0 + Phase(TmpKz4 - TmpKz1 - (TmpKy4 - TmpKy1)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 5, 7, 5, 7)*(1.0 + Phase(TmpKz3 - TmpKz2 - (TmpKy3 - TmpKy2)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 5, 7, 5, 7)*(1.0 + Phase(TmpKz4 - TmpKz2 - (TmpKy4 - TmpKy2)))*epsilon;
 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 2, 3, 2, 3)*(1.0 + Phase(-(TmpKy3 - TmpKy1)))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 2, 3, 2, 3)*(1.0 + Phase(-(TmpKy4 - TmpKy1)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 2, 3, 2, 3)*(1.0 + Phase(-(TmpKy3 - TmpKy2)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 2, 3, 2, 3)*(1.0 + Phase(-(TmpKy4 - TmpKy2)))*epsilon;
 
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 6, 7, 6, 7)*(1.0 + Phase(-(TmpKy3 - TmpKy1)))*epsilon;
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 6, 7, 6, 7)*(1.0 + Phase(-(TmpKy4 - TmpKy1)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 6, 7, 6, 7)*(1.0 + Phase(-(TmpKy3 - TmpKy2)));
 TmpWU += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 6, 7, 6, 7)*(1.0 + Phase(-(TmpKy4 - TmpKy2)))*epsilon;
 
//
 
 
 //Compute contribution to interaction factor of particles with opposite spins
 Complex TmpWV = 0.0;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 0, 5, 0, 5)*(1.0 + Phase(TmpKx3 - TmpKx1 - (TmpKz3 - TmpKz1)))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 0, 5, 0, 5)*(1.0 + Phase(TmpKx4 - TmpKx1 - (TmpKz4 - TmpKz1)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 0, 5, 0, 5)*(1.0 + Phase(TmpKx3 - TmpKx2 - (TmpKz3 - TmpKz2)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 0, 5, 0, 5)*(1.0 + Phase(TmpKx4 - TmpKx2 - (TmpKz4 - TmpKz2)))*epsilon;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 4, 1, 4, 1)*(1.0 + Phase(TmpKx3 - TmpKx1 - (TmpKz3 - TmpKz1)))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 4, 1, 4, 1)*(1.0 + Phase(TmpKx4 - TmpKx1 - (TmpKz4 - TmpKz1)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 4, 1, 4, 1)*(1.0 + Phase(TmpKx3 - TmpKx2 - (TmpKz3 - TmpKz2)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 4, 1, 4, 1)*(1.0 + Phase(TmpKx4 - TmpKx2 - (TmpKz4 - TmpKz2)))*epsilon;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 0, 6, 0, 6)*(1.0 + Phase(TmpKx3 - TmpKx1))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 0, 6, 0, 6)*(1.0 + Phase(TmpKx4 - TmpKx1)); 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 0, 6, 0, 6)*(1.0 + Phase(TmpKx3 - TmpKx2));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 0, 6, 0, 6)*(1.0 + Phase(TmpKx4 - TmpKx2))*epsilon;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 4, 2, 4, 2)*(1.0 + Phase(TmpKx3 - TmpKx1))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 4, 2, 4, 2)*(1.0 + Phase(TmpKx4 - TmpKx1)); 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 4, 2, 4, 2)*(1.0 + Phase(TmpKx3 - TmpKx2));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 4, 2, 4, 2)*(1.0 + Phase(TmpKx4 - TmpKx2))*epsilon;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 0, 7, 0, 7)*(1.0 + Phase(TmpKx3 - TmpKx1 - (TmpKy3 - TmpKy1)))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 0, 7, 0, 7)*(1.0 + Phase(TmpKx4 - TmpKx1 - (TmpKy4 - TmpKy1)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 0, 7, 0, 7)*(1.0 + Phase(TmpKx3 - TmpKx2 - (TmpKy3 - TmpKy2)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 0, 7, 0, 7)*(1.0 + Phase(TmpKx4 - TmpKx2 - (TmpKy4 - TmpKy2)))*epsilon;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 4, 3, 4, 3)*(1.0 + Phase(TmpKx3 - TmpKx1 - (TmpKy3 - TmpKy1)))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 4, 3, 4, 3)*(1.0 + Phase(TmpKx4 - TmpKx1 - (TmpKy4 - TmpKy1)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 4, 3, 4, 3)*(1.0 + Phase(TmpKx3 - TmpKx2 - (TmpKy3 - TmpKy2)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 4, 3, 4, 3)*(1.0 + Phase(TmpKx4 - TmpKx2 - (TmpKy4 - TmpKy2)))*epsilon;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 1, 6, 1, 6)*(1.0 + Phase(TmpKz3 - TmpKz1))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 1, 6, 1, 6)*(1.0 + Phase(TmpKz4 - TmpKz1));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 1, 6, 1, 6)*(1.0 + Phase(TmpKz3 - TmpKz2));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 1, 6, 1, 6)*(1.0 + Phase(TmpKz4 - TmpKz2))*epsilon;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 5, 2, 5, 2)*(1.0 + Phase(TmpKz3 - TmpKz1))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 5, 2, 5, 2)*(1.0 + Phase(TmpKz4 - TmpKz1));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 5, 2, 5, 2)*(1.0 + Phase(TmpKz3 - TmpKz2));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 5, 2, 5, 2)*(1.0 + Phase(TmpKz4 - TmpKz2))*epsilon;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 1, 7, 1, 7)*(1.0 + Phase(TmpKz3 - TmpKz1 - (TmpKy3 - TmpKy1)))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 1, 7, 1, 7)*(1.0 + Phase(TmpKz4 - TmpKz1 - (TmpKy4 - TmpKy1)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 1, 7, 1, 7)*(1.0 + Phase(TmpKz3 - TmpKz2 - (TmpKy3 - TmpKy2)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 1, 7, 1, 7)*(1.0 + Phase(TmpKz4 - TmpKz2 - (TmpKy4 - TmpKy2)))*epsilon;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 5, 3, 5, 3)*(1.0 + Phase(TmpKz3 - TmpKz1 - (TmpKy3 - TmpKy1)))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 5, 3, 5, 3)*(1.0 + Phase(TmpKz4 - TmpKz1 - (TmpKy4 - TmpKy1)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 5, 3, 5, 3)*(1.0 + Phase(TmpKz3 - TmpKz2 - (TmpKy3 - TmpKy2)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 5, 3, 5, 3)*(1.0 + Phase(TmpKz4 - TmpKz2 - (TmpKy4 - TmpKy2)))*epsilon;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 2, 7, 2, 7)*(1.0 + Phase(-(TmpKy3 - TmpKy1)))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 2, 7, 2, 7)*(1.0 + Phase(-(TmpKy4 - TmpKy1)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 2, 7, 2, 7)*(1.0 + Phase(-(TmpKy3 - TmpKy2)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 2, 7, 2, 7)*(1.0 + Phase(-(TmpKy4 - TmpKy2)))*epsilon;
 
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex1, momentumIndex2, energyIndex3, energyIndex4, energyIndex1, energyIndex2, 6, 3, 6, 3)*(1.0 + Phase(-(TmpKy3 - TmpKy1)))*epsilon;
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex1, momentumIndex2, energyIndex4, energyIndex3, energyIndex1, energyIndex2, 6, 3, 6, 3)*(1.0 + Phase(-(TmpKy4 - TmpKy1)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex3, momentumIndex4, momentumIndex2, momentumIndex1, energyIndex3, energyIndex4, energyIndex2, energyIndex1, 6, 3, 6, 3)*(1.0 + Phase(-(TmpKy3 - TmpKy2)));
 TmpWV += this->ComputeTransfomationBasisContribution(oneBodyBasis, momentumIndex4, momentumIndex3, momentumIndex2, momentumIndex1, energyIndex4, energyIndex3, energyIndex2, energyIndex1, 6, 3, 6, 3)*(1.0 + Phase(-(TmpKy4 - TmpKy2)))*epsilon;
 
 TmpWU *= factorWU;
 TmpWV *= factorWV;
 
 Complex Tmp;
 Tmp = TmpWU + TmpWV;
 if (epsilon == 1)
 {
  if ((momentumIndex1 == momentumIndex2) && (energyIndex1 == energyIndex2))
    Tmp *= 0.5;
  if ((momentumIndex3 == momentumIndex4) && (energyIndex3 == energyIndex4))
    Tmp *= 0.5;
 }
 return Tmp;
}