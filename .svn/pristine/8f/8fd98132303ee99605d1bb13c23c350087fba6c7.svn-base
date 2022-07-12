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
#include "Hamiltonian/ParticleOnLatticeOFLNOrbitalTriangularLatticeSingleBandHamiltonian.h"
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
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellX = number of sites in the x direction
// nbrCellY = number of sites in the y direction
// nbrSpinValue = number of internal degrees of freedom
// nbrReciprocalVector = number of vectors of the reciprocal space in each direction to keep. It must be even!
// uPotential = strength of the repulsive two body neareast neighbor interaction
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeOFLNOrbitalTriangularLatticeSingleBandHamiltonian::ParticleOnLatticeOFLNOrbitalTriangularLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, int nbrSpinValue, int  nbrReciprocalVectors, double uPotential, Abstract2DTightBindingModel* tightBindingModel,  bool flatBandFlag, int bandIndex , AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrCellX;
  this->NbrSiteY = nbrCellY;
  this->LzMax = nbrCellX * nbrCellY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->HamiltonianShift = 0.0;
  this->NbrReciprocalVectors =  nbrReciprocalVectors + 1;
  this->NbrSpinValue =  nbrSpinValue;
  this->TightBindingModel = tightBindingModel;
  this->LocalTightBindingModel = ((TightBindingModelOFLNOrbitalTriangularLattice *) tightBindingModel);
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->BandIndex= bandIndex;
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

ParticleOnLatticeOFLNOrbitalTriangularLatticeSingleBandHamiltonian::~ParticleOnLatticeOFLNOrbitalTriangularLatticeSingleBandHamiltonian()
{
}




// evaluate all interaction factors
//   

void ParticleOnLatticeOFLNOrbitalTriangularLatticeSingleBandHamiltonian::EvaluateInteractionFactors()
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
	  this->OneBodyInteractionFactors[Index] = this->TightBindingModel->GetEnergy(0, Index);
	OneBodyBasis[Index] = this->TightBindingModel->GetOneBodyMatrix(Index);
      }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {

      cout <<"Fermionic case not implemented yet"<<endl;
      /*      this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
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
		  

		  this->InteractionFactors[i][Index] *= -2.0;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
      */
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
      double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if (this->FlatBand == false)
	FactorU *= this->UPotential;
      
      cout <<"get here without trouble"<<endl;

      Complex *** TmpCalculationTable = new Complex ** [this->NbrSiteX* this->NbrSiteY];
      
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	  {
	    int Index1 = (kx1 * this->NbrSiteY) + ky1;
	    TmpCalculationTable[Index1] = new Complex * [this->NbrSiteX*this->NbrSiteY];
	    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		{
		  int Index2 = (kx2 * this->NbrSiteY) + ky2;
		  TmpCalculationTable[Index1][Index2] = new Complex [this->NbrReciprocalVectors*this->NbrReciprocalVectors];
		  for (int Nx = 0; Nx < this->NbrReciprocalVectors; Nx++)
		    {
		      for (int Ny = 0; Ny < this->NbrReciprocalVectors; Ny++)
			{
			  int TmpIndex =  Nx * this->NbrReciprocalVectors + Ny;

			  TmpCalculationTable[Index1][Index2][TmpIndex] = 0;
			  
			  for (int Spin = 0; Spin < this->NbrSpinValue; Spin++)
			    {	 
			      for (int Nx1 = 0; Nx1 < this->NbrReciprocalVectors; Nx1++)
				{
				  for (int Ny1 = 0; Ny1 < this->NbrReciprocalVectors; Ny1++)
				    {
				      int IntermediateIndex1 =  this->LocalTightBindingModel->GetIntermediateLinearizedIndices(Nx1,Ny1,Spin);
				      int IntermediateIndex2 =  this->LocalTightBindingModel->GetIntermediateLinearizedIndices(Nx1 - Nx +  this->NbrReciprocalVectors, Ny1 - Ny +  this->NbrReciprocalVectors,Spin);
				      TmpCalculationTable[Index1][Index2][TmpIndex] += Conj(OneBodyBasis[Index1][BandIndex][IntermediateIndex1]) * OneBodyBasis[Index2][BandIndex][IntermediateIndex2];
				    }
				}
			    }
			}
		    }
		}
	  }
  
		
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] =  new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  this->InteractionFactors[i][Index] = 0;
	  
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
		  this->InteractionFactors[i][Index] = 0;
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
		  int ShiftMomentumX = (kx1 + kx2 - kx3 - kx4) /  this->NbrSiteX;
		  int ShiftMomentumY = (ky1 + ky2 - ky3 - ky4) /  this->NbrSiteY;
		  
		  for (int Nx = 0; Nx < this->NbrReciprocalVectors; Nx++)
		    {
		      for (int Ny = 0; Ny < this->NbrReciprocalVectors; Ny++)
			{
			  int TmpIndex =  Nx * this->NbrReciprocalVectors + Ny;
			  int TmpNx =  -ShiftMomentumX + this->NbrReciprocalVectors - Nx;
			  int TmpNy =  -ShiftMomentumY + this->NbrReciprocalVectors - Ny;
			  
			  if (TmpNx < 0)
			    TmpNx += this->NbrReciprocalVectors;

			  if (TmpNx >= this->NbrReciprocalVectors)
			    TmpNx -= this->NbrReciprocalVectors;
			  
			  if (TmpNy < 0)
			    TmpNy += this->NbrReciprocalVectors;
			  
			  if (TmpNy >= this->NbrReciprocalVectors)
			    TmpNy -= this->NbrReciprocalVectors;

			  int TmpIndex2 =  TmpNx *  this->NbrReciprocalVectors + TmpNy;

			  this->InteractionFactors[i][Index] +=  FactorU * TmpCalculationTable[Index1][Index4][TmpIndex] * TmpCalculationTable[Index2][Index3][TmpIndex2];
			  this->InteractionFactors[i][Index] +=  FactorU * TmpCalculationTable[Index2][Index4][TmpIndex] * TmpCalculationTable[Index1][Index3][TmpIndex2];

			  this->InteractionFactors[i][Index] +=  FactorU * TmpCalculationTable[Index1][Index3][TmpIndex] * TmpCalculationTable[Index2][Index4][TmpIndex2];
			  this->InteractionFactors[i][Index] +=  FactorU * TmpCalculationTable[Index2][Index3][TmpIndex] * TmpCalculationTable[Index1][Index4][TmpIndex2];
			}
		    }

		    
		    if (Index3 == Index4)
		      this->InteractionFactors[i][Index] *= 0.5;
		    if (Index1 == Index2)
		      this->InteractionFactors[i][Index] *= 0.5;

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



// evaluate all interaction factors
//   

void ParticleOnLatticeOFLNOrbitalTriangularLatticeSingleBandHamiltonian::OldEvaluateInteractionFactors()
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
	OneBodyBasis[Index] = this->TightBindingModel->GetOneBodyMatrix(Index);
      }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {

      cout <<"Fermionic case not implemented yet"<<endl;
      /*      this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
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
 		  this->InteractionFactors[i][Index] += FactorUBC * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(k2a, k2b, k1a, k1b, k4a, k4b, k3a, k3b);

		  this->InteractionFactors[i][Index] *= -2.0;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	    }
*/
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
		//if (Index1 <= Index2)
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
		//if (Index1 <= Index2)
		//  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		    // }
	      }
      double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if (this->FlatBand == false)
	FactorU *= this->UPotential;
      
      cout <<"get here without trouble"<<endl;
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  this->InteractionFactors[i][Index] = 0;
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
		  this->InteractionFactors[i][Index] = 0;
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
		  int ShiftMomentumX = (kx1 + kx2 - kx3 - kx4) /  this->NbrSiteX;
		  int ShiftMomentumY = (ky1 + ky2 - ky3 - ky4) /  this->NbrSiteY;

		  ///////////////////////////////////////////////////////////////////////////////::
		    
		    for (int Spin1 = 0; Spin1 < this->NbrSpinValue; Spin1++)
		      {
			for (int Spin2 = 0; Spin2 < this->NbrSpinValue; Spin2++)
			  {
			    for (int Nx1 = 0; Nx1 < this->NbrReciprocalVectors; Nx1++)
			      {
				for (int Ny1 = 0; Ny1 < this->NbrReciprocalVectors; Ny1++)
				  {
				    int IntermediateIndex1 =  this->LocalTightBindingModel->GetIntermediateLinearizedIndices(Nx1,Ny1,Spin1);					  
				    for (int Nx2 = 0; Nx2 < this->NbrReciprocalVectors; Nx2++)
				      {
					for (int Ny2 = 0; Ny2 < this->NbrReciprocalVectors; Ny2++)
					  {
					    int IntermediateIndex2 = this->LocalTightBindingModel->GetIntermediateLinearizedIndices(Nx2,Ny2,Spin2);
					    for (int Nx3 = 0; Nx3 < this->NbrReciprocalVectors; Nx3++)
					      {
						for (int Ny3 = 0; Ny3 < this->NbrReciprocalVectors; Ny3++)
						  {
						    int IntermediateIndex3 =  this->LocalTightBindingModel->GetIntermediateLinearizedIndices(Nx3,Ny3,Spin1);
						    
						    int Nx4 = Nx1 + Nx2 - Nx3 +  ShiftMomentumX;
						    int Ny4 = Ny1 + Ny2 - Ny3 +  ShiftMomentumY;
						    int IntermediateIndex4 =  this->LocalTightBindingModel->GetIntermediateLinearizedIndices(Nx4,Ny4,Spin2);
						    
						    this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][BandIndex][IntermediateIndex1]) * OneBodyBasis[Index3][BandIndex][IntermediateIndex3] * Conj(OneBodyBasis[Index2][BandIndex][IntermediateIndex2]) * OneBodyBasis[Index4][BandIndex][IntermediateIndex4]);
						    //this->InteractionFactors[i][Index] = 1;
						  }
					      }
					  }
				      }
				  }
			      }
			  }
		      }
		      /*	    
		    if (Index3 == Index4)
		      this->InteractionFactors[i][Index] *= 0.5;
		    if (Index1 == Index2)
		      this->InteractionFactors[i][Index] *= 0.5;
		    this->InteractionFactors[i][Index] *= 2.0;
*/
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

