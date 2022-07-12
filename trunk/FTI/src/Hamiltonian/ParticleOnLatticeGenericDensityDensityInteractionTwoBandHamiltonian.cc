////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//          class of a generic density-density two body interaction           //
//                       projected onto two bands and                         //
//              assuming a Bloch form for the tight binding model             //
//                                                                            //
//                        last modification : 26/09/2014                      //
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
#include "Hamiltonian/ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian.h"
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

ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian::ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// bandIndex1 = index of the first band onto which the Hamiltonian is projected 
// bandIndex2 = index of the second band onto which the Hamiltonian is projected 
// nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsPotentials = intensity of each density-density term 
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// flatBandOneBodyGap = set the gap between the first band and the second band when using the flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian::ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
																	 int nbrSiteY, int bandIndex1, int bandIndex2,
																	 int* nbrInteractingOrbitals, int** interactingOrbitalsOrbitalIndices,
																	 int** interactingOrbitalsSpatialIndices, double** interactingOrbitalsPotentials,
																	 Abstract2DTightBindingModel* tightBindingModel, 
																	 bool flatBandFlag, double flatBandOneBodyGap,
																	 AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;
  this->FlatBandOneBodyGap = flatBandOneBodyGap;
  this->NbrInteractingOrbitals = new int[this->TightBindingModel->GetNbrBands()];
  this->InteractingOrbitalsOrbitalIndices = new int*[this->TightBindingModel->GetNbrBands()];
  this->InteractingOrbitalsSpatialIndices = new int*[this->TightBindingModel->GetNbrBands()];
  this->InteractingOrbitalsPotentials = new double*[this->TightBindingModel->GetNbrBands()];
  for (int i = 0; i < this->TightBindingModel->GetNbrBands(); ++i)
    {
      this->NbrInteractingOrbitals[i] = nbrInteractingOrbitals[i];
      if (this->NbrInteractingOrbitals[i] > 0)
	{
	  this->InteractingOrbitalsOrbitalIndices[i] = new int[this->NbrInteractingOrbitals[i]];
	  this->InteractingOrbitalsSpatialIndices[i] = new int[2 * this->NbrInteractingOrbitals[i]];
	  this->InteractingOrbitalsPotentials[i] = new double[this->NbrInteractingOrbitals[i]];
	  for (int j = 0; j < this->NbrInteractingOrbitals[i]; ++j)
	    {
	      this->InteractingOrbitalsOrbitalIndices[i][j] = interactingOrbitalsOrbitalIndices[i][j];
	      this->InteractingOrbitalsSpatialIndices[i][2 * j] = interactingOrbitalsSpatialIndices[i][2 * j];
	      this->InteractingOrbitalsSpatialIndices[i][(2 * j) + 1] = interactingOrbitalsSpatialIndices[i][(2 * j) + 1];
	      this->InteractingOrbitalsPotentials[i][j] = interactingOrbitalsPotentials[i][j] / ((double) (this->NbrSiteX * this->NbrSiteY));
	    }
	}
      else
	{
	  this->InteractingOrbitalsOrbitalIndices[i] = new int[1];
	  this->InteractingOrbitalsSpatialIndices[i] = new int[1];
	  this->InteractingOrbitalsPotentials[i] = new double[1];
	}
    }

  this->BandIndex1 = bandIndex1;
  this->BandIndex2 = bandIndex2;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = false;
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

ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian::~ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->TightBindingModel->GetNbrStatePerBand()];
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  this->InteractionFactorsupupupup = 0;
  this->InteractionFactorsupupdowndown = 0;
  this->InteractionFactorsdowndownupup = 0;
  this->InteractionFactorsdowndowndowndown = 0;
  this->InteractionFactorsupdownupup = 0;
  this->InteractionFactorsupdowndowndown = 0;
  this->InteractionFactorsupupupdown = 0;
  this->InteractionFactorsdowndownupdown = 0;
  this->InteractionFactorsupdownupdown = 0;

  int* SigmaToBand = new int[2];
  SigmaToBand[0] = this->BandIndex1;
  SigmaToBand[1] = this->BandIndex2;

  if ((this->FlatBand == false) || (this->FlatBandOneBodyGap != 0.0))
    {
      this->OneBodyInteractionFactorsupup = new double [this->TightBindingModel->GetNbrStatePerBand()];
      this->OneBodyInteractionFactorsdowndown = new double [this->TightBindingModel->GetNbrStatePerBand()];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsupup[Index] = this->TightBindingModel->GetEnergy(this->BandIndex1, Index);
	    this->OneBodyInteractionFactorsdowndown[Index] = this->TightBindingModel->GetEnergy(this->BandIndex2, Index);
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
	  ++this->NbrInterSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, (ky1 + ky2) % this->NbrSiteY)];    
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
	    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, (ky1 + ky2) % this->NbrSiteY);
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	  }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
    for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, (ky1 + ky2) % this->NbrSiteY);
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
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
										     (ky1 + ky2) % this->NbrSiteY);
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      

      Complex Tmp;
      Complex* TmpInteractionFactor;

      int* TmpIndices;
      int* TmpIndices2;
      
      double TmpKx1 = 0.0;
      double TmpKx2 = 0.0;
      double TmpKx3 = 0.0;
      double TmpKx4 = 0.0;
      double TmpKy1 = 0.0;
      double TmpKy2 = 0.0;
      double TmpKy3 = 0.0;
      double TmpKy4 = 0.0;
      
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
		      int kx1, ky1;
		      this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1);
		      int kx2, ky2;
		      this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2);
		      TmpKx1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 0);
		      TmpKx2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0);
		      TmpKy1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 1);
		      TmpKy2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  int Index3 = TmpIndices[i2];
			  int Index4 = TmpIndices[i2 + 1];
			  int kx3, ky3;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
			  int kx4, ky4;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
			  TmpKx3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 0);
			  TmpKx4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0);
			  TmpKy3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 1);
			  TmpKy4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1);
			  (*TmpInteractionFactor) = 0.0;
			  for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
			    {
			      for (int k = 0; k < this->NbrInteractingOrbitals[Orbital1]; ++k)
				{
				  int Orbital2 = this->InteractingOrbitalsOrbitalIndices[Orbital1][k];
				  int XOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][2 * k];
				  int YOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][(2 * k) + 1];
				  double TmpPotential = this->InteractingOrbitalsPotentials[Orbital1][k];
				  (*TmpInteractionFactor) -= (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, 
													  SigmaToBand[sigma3], SigmaToBand[sigma3], 
													  SigmaToBand[sigma1], SigmaToBand[sigma1], 
													  Orbital1, Orbital2, Orbital1, Orbital2)
							      * TmpPotential * Phase ((-(TmpKx2 - TmpKx4)* ((double) XOrbital2))
										      + (-(TmpKy2 - TmpKy4) * ((double) YOrbital2))));
				  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, 
													  SigmaToBand[sigma3], SigmaToBand[sigma3], 
													  SigmaToBand[sigma1], SigmaToBand[sigma1], 
													  Orbital1, Orbital2, Orbital1, Orbital2)
							      * TmpPotential * Phase ((-(TmpKx2 - TmpKx3) * ((double) XOrbital2))
										      + (-(TmpKy2 - TmpKy3) * ((double) YOrbital2))));
				  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, 
													  SigmaToBand[sigma3], SigmaToBand[sigma3], 
													  SigmaToBand[sigma1], SigmaToBand[sigma1], 
													  Orbital1, Orbital2, Orbital1, Orbital2)
							      * TmpPotential * Phase ((-(TmpKx1 - TmpKx4) * ((double) XOrbital2))
										      + (-(TmpKy1 - TmpKy4) * ((double) YOrbital2))));
				  (*TmpInteractionFactor) -= (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, 
													  SigmaToBand[sigma3], SigmaToBand[sigma3], 
													  SigmaToBand[sigma1], SigmaToBand[sigma1], 
													  Orbital1, Orbital2, Orbital1, Orbital2)
							      * TmpPotential * Phase ((-(TmpKx1 - TmpKx3) * ((double) XOrbital2))
										      + (-(TmpKy1 - TmpKy3) * ((double) YOrbital2))));
				}
			    }
			  
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
			  int kx1, ky1;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1);
			  int kx2, ky2;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2);
			  TmpKx1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 0);
			  TmpKx2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0);
			  TmpKy1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 1);
			  TmpKy2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1);
			  for (int i2 = 0; i2 < Lim2; i2 += 2)
			    {
			      int Index3 = TmpIndices2[i2];
			      int Index4 = TmpIndices2[i2 + 1];
			      int kx3, ky3;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
			      int kx4, ky4;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
			      TmpKx3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 0);
			      TmpKx4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0);
			      TmpKy3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 1);
			      TmpKy4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1);

			      (*TmpInteractionFactor) = 0.0;
			      for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
				{
				  for (int k = 0; k < this->NbrInteractingOrbitals[Orbital1]; ++k)
				    {
				      int Orbital2 = this->InteractingOrbitalsOrbitalIndices[Orbital1][k];
				      int XOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][2 * k];
				      int YOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][(2 * k) + 1];
				      double TmpPotential = this->InteractingOrbitalsPotentials[Orbital1][k];
				      (*TmpInteractionFactor) -= (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, 
													      SigmaToBand[sigma3], SigmaToBand[sigma4], 
													      SigmaToBand[sigma1], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx2 - TmpKx4) * ((double) XOrbital2))
											  + (-(TmpKy2 - TmpKy4) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, 
													      SigmaToBand[sigma4], SigmaToBand[sigma3], 
													      SigmaToBand[sigma1], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx2 - TmpKx3) * ((double) XOrbital2))
											  + (-(TmpKy2 - TmpKy3) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, 
													      SigmaToBand[sigma3], SigmaToBand[sigma4], 
													      SigmaToBand[sigma1], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx1 - TmpKx4) * ((double) XOrbital2))
											  + (-(TmpKy1 - TmpKy4) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) -= (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, 
													      SigmaToBand[sigma4], SigmaToBand[sigma3], 
													      SigmaToBand[sigma1], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx1 - TmpKx3) * ((double) XOrbital2))
											  + (-(TmpKy1 - TmpKy3) * ((double) YOrbital2))));
				    }
				}
			      
			     
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
			  int kx1, ky1;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1);
			  int kx2, ky2;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2);
			  TmpKx1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 0);
			  TmpKx2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0);
			  TmpKy1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 1);
			  TmpKy2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1);
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      int Index3 = TmpIndices[i2];
			      int Index4 = TmpIndices[i2 + 1];
			      int kx3, ky3;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
			      int kx4, ky4;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
			      TmpKx3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 0);
			      TmpKx4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0);
			      TmpKy3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 1);
			      TmpKy4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1);

			      (*TmpInteractionFactor) = 0.0;
			      for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
				{
				  for (int k = 0; k < this->NbrInteractingOrbitals[Orbital1]; ++k)
				    {
				      int Orbital2 = this->InteractingOrbitalsOrbitalIndices[Orbital1][k];
				      int XOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][2 * k];
				      int YOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][(2 * k) + 1];
				      double TmpPotential = this->InteractingOrbitalsPotentials[Orbital1][k];
				      (*TmpInteractionFactor) -= (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, 
													      SigmaToBand[sigma3], SigmaToBand[sigma3], 
													      SigmaToBand[sigma1], SigmaToBand[sigma2], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx2 - TmpKx4) * ((double) XOrbital2))
											  + (-(TmpKy2 - TmpKy4) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, 
													      SigmaToBand[sigma3], SigmaToBand[sigma3], 
													      SigmaToBand[sigma1], SigmaToBand[sigma2], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx2 - TmpKx3) * ((double)  XOrbital2))
											  + (-(TmpKy2 - TmpKy3) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, 
													      SigmaToBand[sigma3], SigmaToBand[sigma3], 
													      SigmaToBand[sigma2], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx1 - TmpKx4) * ((double) XOrbital2))
											  + (-(TmpKy1 - TmpKy4) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) -= (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, 
													      SigmaToBand[sigma3], SigmaToBand[sigma3], 
													      SigmaToBand[sigma2], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx1 - TmpKx3) * ((double) XOrbital2))
											  + (-(TmpKy1 - TmpKy3) * ((double) YOrbital2))));
				    }
				}			      
			      if (Index3 == Index4)
				(*TmpInteractionFactor) *= 0.5;
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
			      int kx1, ky1;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1);
			      int kx2, ky2;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2);
			      TmpKx1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 0);
			      TmpKx2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0);
			      TmpKy1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 1);
			      TmpKy2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1);
			      for (int i2 = 0; i2 < Lim2; i2 += 2)
				{
				  int Index3 = TmpIndices2[i2];
				  int Index4 = TmpIndices2[i2 + 1];
				  int kx3, ky3;
				  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
				  int kx4, ky4;
				  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
				  TmpKx3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 0);
				  TmpKx4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0);
				  TmpKy3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 1);
				  TmpKy4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1);

				  (*TmpInteractionFactor) = 0.0;				  
				  for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
				    {
				      for (int k = 0; k < this->NbrInteractingOrbitals[Orbital1]; ++k)
					{
					  int Orbital2 = this->InteractingOrbitalsOrbitalIndices[Orbital1][k];
					  int XOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][2 * k];
					  int YOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][(2 * k) + 1]; 
					  double TmpPotential = this->InteractingOrbitalsPotentials[Orbital1][k];
					  (*TmpInteractionFactor) -= (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, 
														  SigmaToBand[sigma3], SigmaToBand[sigma4], 
														  SigmaToBand[sigma1], SigmaToBand[sigma2], 
														  Orbital1, Orbital2, Orbital1, Orbital2)
								      * TmpPotential * Phase ((-(TmpKx2 - TmpKx4) * ((double)  XOrbital2))
											      + (-(TmpKy2 - TmpKy4) * ((double) YOrbital2))));
					  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, 
														  SigmaToBand[sigma4], SigmaToBand[sigma3], 
														  SigmaToBand[sigma1], SigmaToBand[sigma2], 
														  Orbital1, Orbital2, Orbital1, Orbital2)
								      * TmpPotential * Phase ((-(TmpKx2 - TmpKx3) * ((double) XOrbital2))
											      + (-(TmpKy2 - TmpKy3) * ((double) YOrbital2))));
					  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, 
														  SigmaToBand[sigma3], SigmaToBand[sigma4], 
														  SigmaToBand[sigma2], SigmaToBand[sigma1], 
														  Orbital1, Orbital2, Orbital1, Orbital2)
								      * TmpPotential * Phase ((-(TmpKx1 - TmpKx4) * ((double) XOrbital2))
											      + (-(TmpKy1 - TmpKy4) * ((double) YOrbital2))));
					  (*TmpInteractionFactor) -= (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, 
														  SigmaToBand[sigma4], SigmaToBand[sigma3], 
														  SigmaToBand[sigma2], SigmaToBand[sigma1], 
														  Orbital1, Orbital2, Orbital1, Orbital2)
								      * TmpPotential * Phase ((-(TmpKx1 - TmpKx3) * ((double) XOrbital2))
											      + (-(TmpKy1 - TmpKy3) * ((double) YOrbital2))));
					}
				    }
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
      // bosonic interaction
      
      double TmpKx1 = 0.0;
      double TmpKx2 = 0.0;
      double TmpKx3 = 0.0;
      double TmpKx4 = 0.0;
      double TmpKy1 = 0.0;
      double TmpKy2 = 0.0;
      double TmpKy3 = 0.0;
      double TmpKy4 = 0.0;
      
            
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 <= Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, (ky1 + ky2) % this->NbrSiteY);
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
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 <= Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
										     (ky1 + ky2) % this->NbrSiteY);
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      

      Complex Tmp;
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
		      int kx1, ky1;
		      this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1);
		      int kx2, ky2;
		      this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2);
		      TmpKx1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 0);
		      TmpKx2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0);
		      TmpKy1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 1);
		      TmpKy2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1);
		      for (int i2 = 0; i2 < Lim; i2 += 2)
			{
			  int Index3 = TmpIndices[i2];
			  int Index4 = TmpIndices[i2 + 1];
			  int kx3, ky3;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
			  int kx4, ky4;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
			  TmpKx3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 0);
			  TmpKx4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0);
			  TmpKy3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 1);
			  TmpKy4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1);

			  (*TmpInteractionFactor) = 0.0;
			  for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
			    {
			      for (int k = 0; k < this->NbrInteractingOrbitals[Orbital1]; ++k)
				{
				  int Orbital2 = this->InteractingOrbitalsOrbitalIndices[Orbital1][k];
				  int XOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][2 * k];
				  int YOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][(2 * k) + 1];
				  double TmpPotential = this->InteractingOrbitalsPotentials[Orbital1][k];
				  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, 
													  SigmaToBand[sigma3], SigmaToBand[sigma3], 
													  SigmaToBand[sigma1], SigmaToBand[sigma1], 
													  Orbital1, Orbital2, Orbital1, Orbital2)
							      * TmpPotential * Phase ((-(TmpKx2 - TmpKx4) * ((double) XOrbital2))
										      + (-(TmpKy2 - TmpKy4) * ((double) YOrbital2))));
				  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, 
													  SigmaToBand[sigma3], SigmaToBand[sigma3], 
													  SigmaToBand[sigma1], SigmaToBand[sigma1], 
													  Orbital1, Orbital2, Orbital1, Orbital2)
							      * TmpPotential * Phase ((-(TmpKx2 - TmpKx3) * ((double) XOrbital2))
										      + (-(TmpKy2 - TmpKy3) * ((double) YOrbital2))));
				  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, 
													  SigmaToBand[sigma3], SigmaToBand[sigma3], 
													  SigmaToBand[sigma1], SigmaToBand[sigma1], 
													  Orbital1, Orbital2, Orbital1, Orbital2)
							      * TmpPotential * Phase ((-(TmpKx1 - TmpKx4) * ((double) XOrbital2))
										      + (-(TmpKy1 - TmpKy4) * ((double) YOrbital2))));
				  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, 
													  SigmaToBand[sigma3], SigmaToBand[sigma3], 
													  SigmaToBand[sigma1], SigmaToBand[sigma1], 
													  Orbital1, Orbital2, Orbital1, Orbital2)
							      * TmpPotential * Phase ((-(TmpKx1 - TmpKx3) * ((double) XOrbital2))
										      + (-(TmpKy1 - TmpKy3) * ((double) YOrbital2))));
				}
			    }
			  if (Index1 == Index2)
			    (*TmpInteractionFactor) *= 0.5;
			  if (Index3 == Index4)
			    (*TmpInteractionFactor) *= 0.5;
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
			  int kx1, ky1;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1);
			  int kx2, ky2;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2);
			  TmpKx1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 0);
			  TmpKx2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0);
			  TmpKy1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 1);
			  TmpKy2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1);
			  for (int i2 = 0; i2 < Lim2; i2 += 2)
			    {
			      int Index3 = TmpIndices2[i2];
			      int Index4 = TmpIndices2[i2 + 1];
			      int kx3, ky3;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
			      int kx4, ky4;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
			      TmpKx3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 0);
			      TmpKx4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0);
			      TmpKy3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 1);
			      TmpKy4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1);


			      (*TmpInteractionFactor) = 0.0;
			      for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
				{
				  for (int k = 0; k < this->NbrInteractingOrbitals[Orbital1]; ++k)
				    {
				      int Orbital2 = this->InteractingOrbitalsOrbitalIndices[Orbital1][k];
				      int XOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][2 * k];
				      int YOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][(2 * k) + 1];
				      double TmpPotential = this->InteractingOrbitalsPotentials[Orbital1][k];
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, 
													      SigmaToBand[sigma3], SigmaToBand[sigma4], 
													      SigmaToBand[sigma1], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx2 - TmpKx4) * ((double) XOrbital2))
											  + (-(TmpKy2 - TmpKy4) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, 
													      SigmaToBand[sigma4], SigmaToBand[sigma3], 
													      SigmaToBand[sigma1], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx2 - TmpKx3) * ((double) XOrbital2))
											  + (-(TmpKy2 - TmpKy3) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, 
													      SigmaToBand[sigma3], SigmaToBand[sigma4], 
													      SigmaToBand[sigma1], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx1 - TmpKx4) * ((double) XOrbital2))
											  + (-(TmpKy1 - TmpKy4) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, 
													      SigmaToBand[sigma4], SigmaToBand[sigma3], 
													      SigmaToBand[sigma1], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx1 - TmpKx3) * ((double) XOrbital2))
											  + (-(TmpKy1 - TmpKy3) * ((double) YOrbital2))));
				    }
				}
			      
			      if (Index1 == Index2)
				(*TmpInteractionFactor) *= 0.5;
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
			  int kx1, ky1;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1);
			  int kx2, ky2;
			  this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2);
			  TmpKx1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 0);
			  TmpKx2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0);
			  TmpKy1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 1);
			  TmpKy2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1);
			  for (int i2 = 0; i2 < Lim; i2 += 2)
			    {
			      int Index3 = TmpIndices[i2];
			      int Index4 = TmpIndices[i2 + 1];
			      int kx3, ky3;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
			      int kx4, ky4;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
			      TmpKx3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 0);
			      TmpKx4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0);
			      TmpKy3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 1);
			      TmpKy4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1);

			      (*TmpInteractionFactor) = 0.0;
			      for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
				{
				  for (int k = 0; k < this->NbrInteractingOrbitals[Orbital1]; ++k)
				    {
				      int Orbital2 = this->InteractingOrbitalsOrbitalIndices[Orbital1][k];
				      int XOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][2 * k];
				      int YOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][(2 * k) + 1];
				      double TmpPotential = this->InteractingOrbitalsPotentials[Orbital1][k];
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, 
													      SigmaToBand[sigma3], SigmaToBand[sigma3], 
													      SigmaToBand[sigma1], SigmaToBand[sigma2], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx2 - TmpKx4) * ((double) XOrbital2))
											  + (-(TmpKy2 - TmpKy4) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, 
													      SigmaToBand[sigma3], SigmaToBand[sigma3], 
													      SigmaToBand[sigma1], SigmaToBand[sigma2], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx2 - TmpKx3) * ((double) XOrbital2))
											  + (-(TmpKy2 - TmpKy3) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, 
													      SigmaToBand[sigma3], SigmaToBand[sigma3], 
													      SigmaToBand[sigma2], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx1 - TmpKx4) * ((double) XOrbital2))
											  + (-(TmpKy1 - TmpKy4) * ((double) YOrbital2))));
				      (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, 
													      SigmaToBand[sigma3], SigmaToBand[sigma3], 
													      SigmaToBand[sigma2], SigmaToBand[sigma1], 
													      Orbital1, Orbital2, Orbital1, Orbital2)
								  * TmpPotential * Phase ((-(TmpKx1 - TmpKx3) * ((double) XOrbital2))
											  + (-(TmpKy1 - TmpKy3) * ((double) YOrbital2))));
				    }
				}			      
			      if (Index3 == Index4)
				(*TmpInteractionFactor) *= 0.5;
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
			      int kx1, ky1;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1);
			      int kx2, ky2;
			      this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2);
			      TmpKx1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 0);
			      TmpKx2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 0);
			      TmpKy1 = this->TightBindingModel->GetProjectedMomentum(kx1, ky1, 1);
			      TmpKy2 = this->TightBindingModel->GetProjectedMomentum(kx2, ky2, 1);
			      for (int i2 = 0; i2 < Lim2; i2 += 2)
				{
				  int Index3 = TmpIndices2[i2];
				  int Index4 = TmpIndices2[i2 + 1];
				  int kx3, ky3;
				  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
				  int kx4, ky4;
				  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
				  TmpKx3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 0);
				  TmpKx4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 0);
				  TmpKy3 = this->TightBindingModel->GetProjectedMomentum(kx3, ky3, 1);
				  TmpKy4 = this->TightBindingModel->GetProjectedMomentum(kx4, ky4, 1);

				  (*TmpInteractionFactor) = 0.0;				  
				  for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
				    {
				      for (int k = 0; k < this->NbrInteractingOrbitals[Orbital1]; ++k)
					{
					  int Orbital2 = this->InteractingOrbitalsOrbitalIndices[Orbital1][k];
					  int XOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][2 * k];
					  int YOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][(2 * k) + 1];
					  double TmpPotential = this->InteractingOrbitalsPotentials[Orbital1][k];
					  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index1, Index2, 
														  SigmaToBand[sigma3], SigmaToBand[sigma4], 
														  SigmaToBand[sigma1], SigmaToBand[sigma2], 
														  Orbital1, Orbital2, Orbital1, Orbital2)
								      * TmpPotential * Phase ((-(TmpKx2 - TmpKx4) * ((double) XOrbital2))
											      + (-(TmpKy2 - TmpKy4) * ((double) YOrbital2))));
					  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index1, Index2, 
														  SigmaToBand[sigma4], SigmaToBand[sigma3], 
														  SigmaToBand[sigma1], SigmaToBand[sigma2], 
														  Orbital1, Orbital2, Orbital1, Orbital2)
								      * TmpPotential * Phase ((-(TmpKx2 - TmpKx3) * ((double) XOrbital2))
											      + (-(TmpKy2 - TmpKy3) * ((double) YOrbital2))));
					  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index3, Index4, Index2, Index1, 
														  SigmaToBand[sigma3], SigmaToBand[sigma4], 
														  SigmaToBand[sigma2], SigmaToBand[sigma1], 
														  Orbital1, Orbital2, Orbital1, Orbital2)
								      * TmpPotential * Phase ((-(TmpKx1 - TmpKx4) * ((double) XOrbital2))
											      + (-(TmpKy1 - TmpKy4) * ((double) YOrbital2))));
					  (*TmpInteractionFactor) += (this->ComputeTransfomationBasisContribution(OneBodyBasis, Index4, Index3, Index2, Index1, 
														  SigmaToBand[sigma4], SigmaToBand[sigma3], 
														  SigmaToBand[sigma2], SigmaToBand[sigma1], 
														  Orbital1, Orbital2, Orbital1, Orbital2)
								      * TmpPotential * Phase ((-(TmpKx1 - TmpKx3) * ((double) XOrbital2)))
											      + (-(TmpKy1 - TmpKy3) * ((double) YOrbital2)));
					}
				    }
				  ++TmpInteractionFactor;
				}
			    }
			}
		    }
		}			  
	    }
	}
    }

  delete[] OneBodyBasis;
  delete[] SigmaToBand;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

