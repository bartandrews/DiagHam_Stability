////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//      class of a generic interaction involving only interactions between    //
//        orbitals in the same unit cell, assuming a Bloch form for the       //
//                           the tight binding model                          //
//                                                                            //
//                        last modification : 30/06/2014                      //
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
#include "Hamiltonian/ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandHamiltonian.h"
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

ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandHamiltonian::ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandHamiltonian()
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
// uPotential = repulsive intra-orbital potential strength
// vPotential = repulsive inter-orbital potential strength
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandHamiltonian::ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, 
																	     int bandIndex, double uPotential,  double vPotential,  
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
  this->UFactors = RealMatrix (this->TightBindingModel->GetNbrStatePerBand(), this->TightBindingModel->GetNbrStatePerBand(), true);
  for (int i = 0; i < this->TightBindingModel->GetNbrStatePerBand(); ++i)
    {
      this->UFactors[i][i] = uPotential / ((double) (this->NbrSiteX * this->NbrSiteY));
    }
  if (vPotential != 0.0)
    {
      this->InterOrbitalInteractionFlag = true;
      for (int i = 0; i < this->TightBindingModel->GetNbrStatePerBand(); ++i)
	{
	  for (int j = i + 1; j < this->TightBindingModel->GetNbrStatePerBand(); ++j)
	    {
	      this->UFactors[i][j] = vPotential / ((double) (this->NbrSiteX * this->NbrSiteY));
	      this->UFactors[j][i] = vPotential / ((double) (this->NbrSiteX * this->NbrSiteY));
	    }
	}
    }
  else
    {
      this->InterOrbitalInteractionFlag = false;
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

ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandHamiltonian::~ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeGenericIntraUnitCellInteractionSingleBandHamiltonian::EvaluateInteractionFactors()
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
      cout << "warning, untested code!!!"  << endl;
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
		  Complex sumU = 0.0;
		  this->InteractionFactors[i][Index] = 0.0;
		  if (this->InterOrbitalInteractionFlag == true)
		    {
		      for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
			{
			  for (int Orbital2 = Orbital1 + 1; Orbital2 < this->TightBindingModel->GetNbrBands(); ++Orbital2)
			    {
			      sumU  = Conj(OneBodyBasis[Index1][this->BandIndex][Orbital1]) * OneBodyBasis[Index3][this->BandIndex][Orbital1]
				* Conj(OneBodyBasis[Index2][this->BandIndex][Orbital2]) * OneBodyBasis[Index4][this->BandIndex][Orbital2];
			      sumU -= Conj(OneBodyBasis[Index1][this->BandIndex][Orbital1]) * OneBodyBasis[Index4][this->BandIndex][Orbital1]
				* Conj(OneBodyBasis[Index2][this->BandIndex][Orbital2]) * OneBodyBasis[Index3][this->BandIndex][Orbital2];
			      sumU -= Conj(OneBodyBasis[Index2][this->BandIndex][Orbital1]) * OneBodyBasis[Index3][this->BandIndex][Orbital1]
				* Conj(OneBodyBasis[Index1][this->BandIndex][Orbital2]) * OneBodyBasis[Index4][this->BandIndex][Orbital2];
			      sumU += Conj(OneBodyBasis[Index2][this->BandIndex][Orbital1]) * OneBodyBasis[Index4][this->BandIndex][Orbital1]
				* Conj(OneBodyBasis[Index1][this->BandIndex][Orbital2]) * OneBodyBasis[Index3][this->BandIndex][Orbital2];
			      this->InteractionFactors[i][Index] += 2.0 * this->UFactors[Orbital1][Orbital2] *sumU;
			    }
			}
		    }
			
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
		  Complex sumU = 0.0;
		  this->InteractionFactors[i][Index] = 0.0;
		  for (int Orbital = 0; Orbital < this->TightBindingModel->GetNbrBands(); ++Orbital)
		    {
		      sumU  = Conj(OneBodyBasis[Index1][this->BandIndex][Orbital]) * OneBodyBasis[Index3][this->BandIndex][Orbital]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][Orbital]) * OneBodyBasis[Index4][this->BandIndex][Orbital];
		      sumU += Conj(OneBodyBasis[Index1][this->BandIndex][Orbital]) * OneBodyBasis[Index4][this->BandIndex][Orbital]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][Orbital]) * OneBodyBasis[Index3][this->BandIndex][Orbital];
		      sumU += Conj(OneBodyBasis[Index2][this->BandIndex][Orbital]) * OneBodyBasis[Index3][this->BandIndex][Orbital]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][Orbital]) * OneBodyBasis[Index4][this->BandIndex][Orbital];
		      sumU += Conj(OneBodyBasis[Index2][this->BandIndex][Orbital]) * OneBodyBasis[Index4][this->BandIndex][Orbital]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][Orbital]) * OneBodyBasis[Index3][this->BandIndex][Orbital];
		      this->InteractionFactors[i][Index] += this->UFactors[Orbital][Orbital] *sumU; 
		    }
		  if (this->InterOrbitalInteractionFlag == true)
		    {
		      for (int Orbital1 = 0; Orbital1 < this->TightBindingModel->GetNbrBands(); ++Orbital1)
			{
			  for (int Orbital2 = Orbital1 + 1; Orbital2 < this->TightBindingModel->GetNbrBands(); ++Orbital2)
			    {
			      sumU  = Conj(OneBodyBasis[Index1][this->BandIndex][Orbital1]) * OneBodyBasis[Index3][this->BandIndex][Orbital1]
				* Conj(OneBodyBasis[Index2][this->BandIndex][Orbital2]) * OneBodyBasis[Index4][this->BandIndex][Orbital2];
			      sumU += Conj(OneBodyBasis[Index1][this->BandIndex][Orbital1]) * OneBodyBasis[Index4][this->BandIndex][Orbital1]
				* Conj(OneBodyBasis[Index2][this->BandIndex][Orbital2]) * OneBodyBasis[Index3][this->BandIndex][Orbital2];
			      sumU += Conj(OneBodyBasis[Index2][this->BandIndex][Orbital1]) * OneBodyBasis[Index3][this->BandIndex][Orbital1]
				* Conj(OneBodyBasis[Index1][this->BandIndex][Orbital2]) * OneBodyBasis[Index4][this->BandIndex][Orbital2];
			      sumU += Conj(OneBodyBasis[Index2][this->BandIndex][Orbital1]) * OneBodyBasis[Index4][this->BandIndex][Orbital1]
				* Conj(OneBodyBasis[Index1][this->BandIndex][Orbital2]) * OneBodyBasis[Index3][this->BandIndex][Orbital2];
			      this->InteractionFactors[i][Index] += 2.0 * this->UFactors[Orbital1][Orbital2] *sumU;
			    }
			}
		    }
			
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

