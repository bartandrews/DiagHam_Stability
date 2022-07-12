////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                          class author: Gunnar Möller                       //
//                                                                            //
//      class of Hamiltonian of bosons (currently) in a generic optical       //
//                    flux lattice with extended Brillouin zone               //
//                                                                            //
//                        last modification : 01/07/2014                      //
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
#include "Hamiltonian/ParticleOnLatticeOFLGenericLatticeWithSymmetrySingleBandHamiltonian.h"
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
// tightBindingModel = Generic Tight Binding Model With Symmetry
// uPotential = spin-independent strength of the repulsive two body contact interactions
// flatBandFlag = use flat band model
// noDispersionFlag = remove the band dispersion and put a constant gap of 10 between the two lowest bands
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
ParticleOnLatticeOFLGenericLatticeWithSymmetrySingleBandHamiltonian::ParticleOnLatticeOFLGenericLatticeWithSymmetrySingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, TightBindingModelOFLGenericLatticeWithSymmetry* tightBindingModel, double uPotential, bool flatBandFlag, bool noDispersionFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->LocalTightBindingModel = tightBindingModel;
  this->TightBindingModel = ((Abstract2DTightBindingModel*) tightBindingModel);
  // these two entries not needed:
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  // 
  this->HamiltonianShift = 0.0;
  this->NbrReciprocalVectors1 = this->LocalTightBindingModel->GetNbrReciprocalVectors1();
  this->NbrReciprocalVectors2 = this->LocalTightBindingModel->GetNbrReciprocalVectors2();
  this->SymmetryMultiplier1 = this->LocalTightBindingModel->GetSymmetryMultiplier1();
  this->SymmetryMultiplier2 = this->LocalTightBindingModel->GetSymmetryMultiplier2();
  // Brillouin zone is enlarged in tight binding model, taking into account symmetry multipliers!
  this->NbrSiteX = this->LocalTightBindingModel->GetNbrSiteX();
  this->NbrSiteY = this->LocalTightBindingModel->GetNbrSiteY();
  this->LzMax = NbrSiteX * NbrSiteY - 1;
  this->FlatBand = flatBandFlag;
  this->NoDispersionFlag = noDispersionFlag;
  this->UPotentialMatrix.Resize(tightBindingModel->GetNbrSubLatticeFlavours(), tightBindingModel->GetNbrSubLatticeFlavours());
  this->UPotential=uPotential;
  UPotentialMatrix.SetAllEntries(uPotential);
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

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// tightBindingModel = Generic Tight Binding Model With Symmetry
// uPotentialMatrix = spin-dependent strength of the repulsive two body contact interactions
// flatBandFlag = use flat band model
// noDispersionFlag = remove the band dispersion and put a constant gap of 10 between the two lowest bands
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
ParticleOnLatticeOFLGenericLatticeWithSymmetrySingleBandHamiltonian::ParticleOnLatticeOFLGenericLatticeWithSymmetrySingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, TightBindingModelOFLGenericLatticeWithSymmetry* tightBindingModel, RealSymmetricMatrix &uPotentialMatrix, bool flatBandFlag, bool noDispersionFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->LocalTightBindingModel = tightBindingModel;
  this->TightBindingModel = ((Abstract2DTightBindingModel*) tightBindingModel);
  // these two entries not needed:
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  // 
  this->HamiltonianShift = 0.0;
  this->NbrReciprocalVectors1 = this->LocalTightBindingModel->GetNbrReciprocalVectors1();
  this->NbrReciprocalVectors2 = this->LocalTightBindingModel->GetNbrReciprocalVectors2();
  this->SymmetryMultiplier1 = this->LocalTightBindingModel->GetSymmetryMultiplier1();
  this->SymmetryMultiplier2 = this->LocalTightBindingModel->GetSymmetryMultiplier2();
  // Brillouin zone is enlarged in tight binding model, taking into account symmetry multipliers!
  this->NbrSiteX = this->LocalTightBindingModel->GetNbrSiteX();
  this->NbrSiteY = this->LocalTightBindingModel->GetNbrSiteY();
  this->LzMax = NbrSiteX * NbrSiteY - 1;
  this->FlatBand = flatBandFlag;
  this->NoDispersionFlag = noDispersionFlag;
  this->UPotentialMatrix = uPotentialMatrix;
  if (tightBindingModel->GetNbrSubLatticeFlavours()!=UPotentialMatrix.GetNbrRow())
    {
      cout << "Inconsistent size of interaction matrix uPotential"<<endl;
    }
  this->UPotential = UPotentialMatrix(0,0);
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

ParticleOnLatticeOFLGenericLatticeWithSymmetrySingleBandHamiltonian::~ParticleOnLatticeOFLGenericLatticeWithSymmetrySingleBandHamiltonian()
{
}




// evaluate all interaction factors
//   

void ParticleOnLatticeOFLGenericLatticeWithSymmetrySingleBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  
  this->InteractionFactors = 0;
  
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->LocalTightBindingModel->GetNbrStatePerBand()];

  if (this->FlatBand == false)
    this->OneBodyInteractionFactors = new double [this->LocalTightBindingModel->GetNbrStatePerBand()];

  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = this->LocalTightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	if (this->FlatBand == false)
	  {
	    if (this->NoDispersionFlag == true)
	      this->OneBodyInteractionFactors[Index] = 0;
	    else
	      this->OneBodyInteractionFactors[Index] = this->LocalTightBindingModel->GetEnergy(0, Index); 
	  }
	OneBodyBasis[Index] = this->LocalTightBindingModel->GetOneBodyMatrix(Index);
      }

  this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
  for (int i = 0; i < this->NbrSectorSums; ++i)
    this->NbrSectorIndicesPerSum[i] = 0;


  int NbrSpinFlavours = this->LocalTightBindingModel->GetNbrSubLatticeFlavours();
  int NbrFlavourSublattices = this->LocalTightBindingModel->GetNbrSubLatticesPerFlavour();

  Complex **** DensityOperator = new Complex *** [NbrSpinFlavours];
      
  for (int Spin = 0; Spin < NbrSpinFlavours; ++Spin)	
    {
      DensityOperator[Spin] = new Complex ** [this->NbrSiteX* this->NbrSiteY];
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	  {
	    int Index1 = (kx1 * this->NbrSiteY) + ky1;
	    DensityOperator[Spin][Index1] = new Complex * [this->NbrSiteX*this->NbrSiteY];
	      
	    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		{
		  int Index2 = (kx2 * this->NbrSiteY) + ky2;
		  DensityOperator[Spin][Index1][Index2] = new Complex [this->NbrReciprocalVectors1*this->NbrReciprocalVectors2];

		  for (int Nx = 0; Nx < this->NbrReciprocalVectors1; Nx++)
		    {
		      for (int Ny = 0; Ny < this->NbrReciprocalVectors2; Ny++)
			{
			  int TmpIndex =  Nx * this->NbrReciprocalVectors2 + Ny;
			  
			  DensityOperator[Spin][Index1][Index2][TmpIndex] = 0;
			  for (int Subl = 0; Subl < NbrFlavourSublattices; Subl++)
			    {	 
			      for (int Nx1 = 0; Nx1 < this->NbrReciprocalVectors1; Nx1++)
				{
				  for (int Ny1 = 0; Ny1 < this->NbrReciprocalVectors2; Ny1++)
				    {
				      int TotalSubl =  this->LocalTightBindingModel->TotalSublatticeIndex(Spin, Subl);
				      int ReciprocalVectorIndex1 =  this->LocalTightBindingModel->PeriodicLinearizedReciprocalSpaceIndex(Nx1,Ny1,TotalSubl);
				      int ReciprocalVectorIndex2 =  this->LocalTightBindingModel->PeriodicLinearizedReciprocalSpaceIndex(Nx1 - Nx +  this->NbrReciprocalVectors1, Ny1 - Ny +  this->NbrReciprocalVectors2,TotalSubl);
				      
				      DensityOperator[Spin][Index1][Index2][TmpIndex] += Conj(OneBodyBasis[Index1][0][ReciprocalVectorIndex1]) * OneBodyBasis[Index2][0][ReciprocalVectorIndex2];
				    }
				}
			    }
			}
		    }
		}
	  }
    }
  
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      cout <<"Fermionic case not tested, yet"<<endl;
      
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
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }

      if (this->FlatBand == true)
	{
	  if (this->UPotentialMatrix(0,0)!=0.0)
	    this->UPotentialMatrix *= 1.0/this->UPotentialMatrix(0,0);
	  else
	    this->UPotentialMatrix.SetAllEntries(1.0);
	}
      this->UPotentialMatrix *= this->SymmetryMultiplier1 * this->SymmetryMultiplier2 * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));      
      
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 <  this->NbrSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      
	      for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
		  
		  this->InteractionFactors[i][Index] = 0.0;

		  int ShiftMomentumX = (kx1 + kx2 - kx3 - kx4) /  this->NbrSiteX;
		  int ShiftMomentumY = (ky1 + ky2 - ky3 - ky4) /  this->NbrSiteY;

		  for (int Nx = 0; Nx < this->NbrReciprocalVectors1; Nx++)
		    {
		      for (int Ny = 0; Ny < this->NbrReciprocalVectors2; Ny++)
			{
			  int TmpIndex =  Nx * this->NbrReciprocalVectors2 + Ny;
			  int TmpNx =  -ShiftMomentumX + this->NbrReciprocalVectors1 - Nx;
			  int TmpNy =  -ShiftMomentumY + this->NbrReciprocalVectors2 - Ny;
			  if (TmpNx < 0)
			    TmpNx += this->NbrReciprocalVectors1;

			  if (TmpNx >= this->NbrReciprocalVectors1)
			    TmpNx -= this->NbrReciprocalVectors1;
			  
			  if (TmpNy < 0)
			    TmpNy += this->NbrReciprocalVectors2;
			  
			  if (TmpNy >= this->NbrReciprocalVectors2)
			    TmpNy -= this->NbrReciprocalVectors2;
			  
			  int TmpIndex2 =  TmpNx *  this->NbrReciprocalVectors2 + TmpNy;
			  for (int Spin1 = 0; Spin1 < NbrSpinFlavours; ++Spin1)
			    for (int Spin2 = 0; Spin2 < NbrSpinFlavours; ++Spin2)
			      if (Spin1!=Spin2) // equal spins don't interact due to Pauli exclusion
				{
				  this->InteractionFactors[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperator[Spin1][Index1][Index4][TmpIndex] * DensityOperator[Spin2][Index2][Index3][TmpIndex2];
				  this->InteractionFactors[i][Index] -=  UPotentialMatrix(Spin1,Spin2) * DensityOperator[Spin1][Index2][Index4][TmpIndex] * DensityOperator[Spin2][Index1][Index3][TmpIndex2];
				  this->InteractionFactors[i][Index] -=  UPotentialMatrix(Spin1,Spin2) * DensityOperator[Spin1][Index1][Index3][TmpIndex] * DensityOperator[Spin2][Index2][Index4][TmpIndex2];
				  this->InteractionFactors[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperator[Spin1][Index2][Index3][TmpIndex] * DensityOperator[Spin2][Index1][Index4][TmpIndex2];
				}
			}
		    }

		  this->InteractionFactors[i][Index] *= -2.0;
		  ++TotalNbrInteractionFactors;
		  ++Index;
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
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 <= Index2)
		  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }

      if (this->FlatBand == true)
	{
	  if (this->UPotentialMatrix(0,0)!=0.0)
	    this->UPotentialMatrix *= 1.0/this->UPotentialMatrix(0,0);
	  else
	    this->UPotentialMatrix.SetAllEntries(1.0);
	}
      this->UPotentialMatrix *= this->SymmetryMultiplier1 * this->SymmetryMultiplier2 * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));      
      
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j2 = 0; j2 <  this->NbrSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      
	      for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
		  
		  this->InteractionFactors[i][Index] = 0.0;

		  // int ShiftMomentumX;
		  // int ShiftMomentumY;
		  // this->LocalTightBindingModel->Periodize(kx1 + kx2 - kx3 - kx4, this->NbrSiteX, ShiftMomentumX);
		  // ShiftMomentumX /= -this->NbrSiteX;
		  
		  // this->LocalTightBindingModel->Periodize(ky1 + ky2 - ky3 - ky4, this->NbrSiteY, ShiftMomentumY);
		  // ShiftMomentumY /= -this->NbrSiteY;

		  int ShiftMomentumX = (kx1 + kx2 - kx3 - kx4) /  this->NbrSiteX;
		  int ShiftMomentumY = (ky1 + ky2 - ky3 - ky4) /  this->NbrSiteY;


		  for (int Nx = 0; Nx < this->NbrReciprocalVectors1; Nx++)
		    {
		      for (int Ny = 0; Ny < this->NbrReciprocalVectors2; Ny++)
			{
			  int TmpIndex =  Nx * this->NbrReciprocalVectors2 + Ny;
			  int TmpNx =  -ShiftMomentumX + this->NbrReciprocalVectors1 - Nx;
			  int TmpNy =  -ShiftMomentumY + this->NbrReciprocalVectors2 - Ny;
			  if (TmpNx < 0)
			    TmpNx += this->NbrReciprocalVectors1;

			  if (TmpNx >= this->NbrReciprocalVectors1)
			    TmpNx -= this->NbrReciprocalVectors1;
			  
			  if (TmpNy < 0)
			    TmpNy += this->NbrReciprocalVectors2;
			  
			  if (TmpNy >= this->NbrReciprocalVectors2)
			    TmpNy -= this->NbrReciprocalVectors2;
			  
			  int TmpIndex2 =  TmpNx *  this->NbrReciprocalVectors2 + TmpNy;
			  for (int Spin1 = 0; Spin1 < NbrSpinFlavours; ++Spin1)
			    for (int Spin2 = 0; Spin2 < NbrSpinFlavours; ++Spin2)
			      {
				this->InteractionFactors[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperator[Spin1][Index1][Index4][TmpIndex] * DensityOperator[Spin2][Index2][Index3][TmpIndex2];
				this->InteractionFactors[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperator[Spin1][Index2][Index4][TmpIndex] * DensityOperator[Spin2][Index1][Index3][TmpIndex2];
				this->InteractionFactors[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperator[Spin1][Index1][Index3][TmpIndex] * DensityOperator[Spin2][Index2][Index4][TmpIndex2];
				this->InteractionFactors[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperator[Spin1][Index2][Index3][TmpIndex] * DensityOperator[Spin2][Index1][Index4][TmpIndex2];
			      }
			 
			}
		    }

		  if (Index3 == Index4)
		    {
		      this->InteractionFactors[i][Index] *= 0.5;
		    }
		  if (Index1 == Index2)
		    {
		      this->InteractionFactors[i][Index] *= 0.5;
		    }
		  //this->InteractionFactors[i][Index] *= 2.0;
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
      

  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
  
  for (int Spin = 0; Spin < NbrSpinFlavours; ++Spin)
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	  {
	    int Index1 = (kx1 * this->NbrSiteY) + ky1;
	    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		{
		  int Index2 = (kx2 * this->NbrSiteY) + ky2;
		  delete [] DensityOperator[Spin][Index1][Index2];
		}
	    delete [] DensityOperator[Spin][Index1];
	  }
      delete [] DensityOperator[Spin];
    }
  delete [] DensityOperator;
    
  delete [] OneBodyBasis;
}

