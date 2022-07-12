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
//                     projected onto a single band and                       //
//              assuming a Bloch form for the tight binding model             //
//                                                                            //
//                        last modification : 25/09/2014                      //
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
#include "Hamiltonian/ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian.h"
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

ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian::ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian()
{
  this->BandIndex = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellX = number of sites in the x direction
// nbrCellY = number of sites in the y direction
// bandIndex = index of the band in which the Hamiltonian is projected 
// nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsPotentials = intensity of each density-density term 
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian::ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, 
																 int bandIndex, int* nbrInteractingOrbitals, int** interactingOrbitalsOrbitalIndices,
																 int** interactingOrbitalsSpatialIndices, double** interactingOrbitalsPotentials,  
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

ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian::~ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian::EvaluateInteractionFactors()
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
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 < Index2)
		  ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2)];    
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
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1;
	      int ky1;
	      int kx2;
	      int ky2;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1);
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2);
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3;
		  int ky3;
		  int kx4;
		  int ky4;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
		  Complex& SumU = this->InteractionFactors[i][Index];
		  SumU = 0.0;
		  for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
		    {
		      for (int k = 0; k < this->NbrInteractingOrbitals[Orbital1]; ++k)
			{
			  int Orbital2 = this->InteractingOrbitalsOrbitalIndices[Orbital1][k];
			  int XOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][2 * k];
			  int YOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][(2 * k) + 1];
			  double TmpPotential = this->InteractingOrbitalsPotentials[Orbital1][k];
			  SumU -= (Conj(OneBodyBasis[Index1][this->BandIndex][Orbital1]) * OneBodyBasis[Index3][this->BandIndex][Orbital1]
				   * Conj(OneBodyBasis[Index2][this->BandIndex][Orbital2]) * OneBodyBasis[Index4][this->BandIndex][Orbital2]
				   * TmpPotential * Phase ((this->KxFactor * ((double) ((kx2 - kx4) * XOrbital2))) 
							   + (this->KyFactor * ((double) ((ky2 - ky4) * YOrbital2)))));
			  SumU += (Conj(OneBodyBasis[Index1][this->BandIndex][Orbital1]) * OneBodyBasis[Index4][this->BandIndex][Orbital1]
				   * Conj(OneBodyBasis[Index2][this->BandIndex][Orbital2]) * OneBodyBasis[Index3][this->BandIndex][Orbital2]
				   * TmpPotential * Phase ((this->KxFactor * ((double) ((kx2 - kx3) * XOrbital2))) 
							   + (this->KyFactor * ((double) ((ky2 - ky3) * YOrbital2)))));
			  SumU += (Conj(OneBodyBasis[Index2][this->BandIndex][Orbital1]) * OneBodyBasis[Index3][this->BandIndex][Orbital1]
				   * Conj(OneBodyBasis[Index1][this->BandIndex][Orbital2]) * OneBodyBasis[Index4][this->BandIndex][Orbital2]
				   * TmpPotential * Phase ((this->KxFactor * ((double) ((kx1 - kx4) * XOrbital2))) 
							   + (this->KyFactor * ((double) ((ky1 - ky4) * YOrbital2)))));
			  SumU -= (Conj(OneBodyBasis[Index2][this->BandIndex][Orbital1]) * OneBodyBasis[Index4][this->BandIndex][Orbital1]
				   * Conj(OneBodyBasis[Index1][this->BandIndex][Orbital2]) * OneBodyBasis[Index3][this->BandIndex][Orbital2]
				   * TmpPotential * Phase ((this->KxFactor * ((double) ((kx1 - kx3) * XOrbital2))) 
							   + (this->KyFactor * ((double) ((ky1 - ky3) * YOrbital2)))));
			}
		    }			
		  TotalNbrInteractionFactors++;
		  ++Index;
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
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 <= Index2)
		  ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2)];    
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
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 <= Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1;
	      int ky1;
	      int kx2;
	      int ky2;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1);
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2);
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3;
		  int ky3;
		  int kx4;
		  int ky4;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);
		  Complex& SumU = this->InteractionFactors[i][Index];
		  SumU = 0.0;
		  for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
		    {
		      for (int k = 0; k < this->NbrInteractingOrbitals[Orbital1]; ++k)
			{
			  int Orbital2 = this->InteractingOrbitalsOrbitalIndices[Orbital1][k];
			  int XOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][2 * k];
			  int YOrbital2 = this->InteractingOrbitalsSpatialIndices[Orbital1][(2 * k) + 1];
			  double TmpPotential = this->InteractingOrbitalsPotentials[Orbital1][k];
			  SumU += (Conj(OneBodyBasis[Index1][this->BandIndex][Orbital1]) * OneBodyBasis[Index3][this->BandIndex][Orbital1]
				   * Conj(OneBodyBasis[Index2][this->BandIndex][Orbital2]) * OneBodyBasis[Index4][this->BandIndex][Orbital2]
				   * TmpPotential * Phase ((this->KxFactor * ((double) ((kx2 - kx4) * XOrbital2))) 
							   + (this->KyFactor * ((double) ((ky2 - ky4) * YOrbital2)))));
			  SumU += (Conj(OneBodyBasis[Index1][this->BandIndex][Orbital1]) * OneBodyBasis[Index4][this->BandIndex][Orbital1]
				   * Conj(OneBodyBasis[Index2][this->BandIndex][Orbital2]) * OneBodyBasis[Index3][this->BandIndex][Orbital2]
				   * TmpPotential * Phase ((this->KxFactor * ((double) ((kx2 - kx3) * XOrbital2))) 
							   + (this->KyFactor * ((double) ((ky2 - ky3) * YOrbital2)))));
			  SumU += (Conj(OneBodyBasis[Index2][this->BandIndex][Orbital1]) * OneBodyBasis[Index3][this->BandIndex][Orbital1]
				   * Conj(OneBodyBasis[Index1][this->BandIndex][Orbital2]) * OneBodyBasis[Index4][this->BandIndex][Orbital2]
				   * TmpPotential * Phase ((this->KxFactor * ((double) ((kx1 - kx4) * XOrbital2))) 
							   + (this->KyFactor * ((double) ((ky1 - ky4) * YOrbital2)))));
			  SumU += (Conj(OneBodyBasis[Index2][this->BandIndex][Orbital1]) * OneBodyBasis[Index4][this->BandIndex][Orbital1]
				   * Conj(OneBodyBasis[Index1][this->BandIndex][Orbital2]) * OneBodyBasis[Index3][this->BandIndex][Orbital2]
				   * TmpPotential * Phase ((this->KxFactor * ((double) ((kx1 - kx3) * XOrbital2))) 
							   + (this->KyFactor * ((double) ((ky1 - ky3) * YOrbital2)))));
			}
		    }			
		  if (Index3 == Index4)
		    SumU *= 0.5;
		  if (Index1 == Index2)
		    SumU *= 0.5;

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

