////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                          class author: Gunnar Möller                       //
//                                                                            //
//      class for general Hofstadter models with interacting particles        //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 10/12/2015                      //
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
#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandGenericHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Tools/FTITightBinding/AbstractTightBindingInteraction.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/StringTools.h"
#include "MathTools/KahanSum.h"
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
ParticleOnLatticeHofstadterSingleBandGenericHamiltonian::ParticleOnLatticeHofstadterSingleBandGenericHamiltonian()
{
  this->BandIndex = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellX = number of sites in the x direction
// nbrCellY = number of sites in the y direction
// bandIndex = index of band to consider
// genericInteraction = pointer to object encoding interactions
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
ParticleOnLatticeHofstadterSingleBandGenericHamiltonian::ParticleOnLatticeHofstadterSingleBandGenericHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, int bandIndex, AbstractTightBindingInteraction* interaction,  Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
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
  this->BandIndex = bandIndex;
  this->Interaction = interaction;
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

ParticleOnLatticeHofstadterSingleBandGenericHamiltonian::~ParticleOnLatticeHofstadterSingleBandGenericHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeHofstadterSingleBandGenericHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  int NbrSublattices = TightBindingModel->GetNbrBands();
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
	  this->OneBodyInteractionFactors[Index] = this->TightBindingModel->GetEnergy(BandIndex, Index);
	OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
      }

  if (this->FlatBand == false)
    for (int i=0; i<this->TightBindingModel->GetNbrStatePerBand(); ++i)
      {
	cout << "[" << this->OneBodyInteractionFactors[i] << "]" << endl;
      }

  // double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
  // if (this->FlatBand == false)
  //   FactorU *= this->UPotential;
  // double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

  // if (FactorU==0.0 && FactorV==0.0)
  //   {
  //     std::cerr << "Error: HofstadterHamiltonian created with interaction zero - set non-zero --u-potential or --v-potential"<<std::endl;
  //     exit(1);
  //   }

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
		  ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2)];
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
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }

      this->InteractionFactors = new Complex* [this->NbrSectorSums];

      Complex Tmp, ATerm;
      // int numX, numY;
      // this->TightBindingModel->GetMUCDimensions(numX, numY);
      // Interaction tri(60,this->NbrSiteX*numX,this->NbrSiteY*numY,1,1,100);
      // cout << "Calculated Ewald summation" << endl;

      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1, ky1;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
	      int kx2, ky2;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);
	      
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3, ky3;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
		  int kx4, ky4;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

		  this->InteractionFactors[i][Index] = 0.0;

		  Tmp=0.0;
		  KahamSum<Complex> Tmp2;
            
		  for (int s=0; s<NbrSublattices; ++s)
		    {		      
		      for (int sF=0; sF<NbrSublattices; ++sF) // final sublattice
			{
			  for (int dRx = 0; dRx < this->NbrSiteX; ++dRx)
			    for (int dRy = 0; dRy < this->NbrSiteY; ++dRy)
			      {	
				double Vij = this->Interaction->GetAmplitude(s, dRx, dRy, sF);
                    
				ATerm = Vij * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx3, ky3);
				Tmp += ATerm;
				Tmp2 += ATerm;
				ATerm = Vij * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx4, ky4);
				Tmp -= ATerm;
				Tmp2 -= ATerm;
				ATerm = Vij * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx3, ky3);
				Tmp -= ATerm;
				Tmp2 -= ATerm;
				ATerm = Vij * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx4, ky4);
				Tmp += ATerm;
				Tmp2 += ATerm;
			      }
			}
		    }
		  cout << "Accuracy of sum: "<<Norm(Tmp - Tmp2.GetValue()) << endl;
		  this->InteractionFactors[i][Index] += 2.0 * Tmp2.GetValue();
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else // Bosonic Statistics
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
		  ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2)];
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
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }

      this->InteractionFactors = new Complex* [this->NbrSectorSums];

      Complex Tmp, ATerm;

      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1,ky1;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
	      int kx2,ky2;
	      this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);
	      
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3,ky3;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
		  int kx4,ky4;
		  this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

		  Tmp=0.0;
		  KahamSum<Complex> Tmp2;
            
		  for (int s=0; s<NbrSublattices; ++s)
		    {		      
		      for (int sF=0; sF<NbrSublattices; ++sF) // final sublattice
			{
			  for (int dRx = 0; dRx < this->NbrSiteX; ++dRx)
			    for (int dRy = 0; dRy < this->NbrSiteY; ++dRy)
			      {	
				double Vij = this->Interaction->GetAmplitude(s, dRx, dRy, sF);
                    
				ATerm = Vij * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx3, ky3);
				Tmp += ATerm;
				Tmp2 += ATerm;
				ATerm = Vij * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx4, ky4);
				Tmp += ATerm;
				Tmp2 += ATerm;
				ATerm = Vij * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx3, ky3);
				Tmp += ATerm;
				Tmp2 += ATerm;
				ATerm = Vij * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, sF, sF, s) * this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx4, ky4);
				Tmp += ATerm;
				Tmp2 += ATerm;
			      }
			}
                  	  
		      if (Index3 == Index4)
			{ 
			  Tmp *= 0.5; Tmp2 *= 0.5;
			}
		      if (Index1 == Index2)
			{
			  Tmp *= 0.5; Tmp2 *= 0.5;
			}

		      cout << "Accuracy of sum: "<<Norm(Tmp - Tmp2.GetValue()) << endl;

		      this->InteractionFactors[i][Index] = 2.0 * Tmp2.GetValue();
		  		  
		      ++TotalNbrInteractionFactors;
		      ++Index;
		    }
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;

  delete [] OneBodyBasis;
}



// compute the phase factor for the embedding of four operators of an on-site two body interaction involving sites on a generic sublattic 
//
// subl = sublattice index
// kx1 = first creation momentum along x for the B site
// ky1 = first creation momentum along y for the B site
// kx2 = second creation momentum along x for the B site
// ky2 = second creation momentum along y for the B site
// kx3 = first annihilation momentum along x for the B site
// ky3 = first annihilation momentum along y for the B site
// kx4 = second annihilation momentum along x for the B site
// ky4 = second annihilation momentum along y for the B site
//
// return value = corresponding matrix element
Complex ParticleOnLatticeHofstadterSingleBandGenericHamiltonian::ComputeEmbeddingOnSite(int subl, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  double embeddingX, embeddingY;
  this->TightBindingModel->GetEmbedding(subl, embeddingX, embeddingY);
  double phase = (this->KxFactor * (kx3 + kx4 - kx1 - kx2) * embeddingX + this->KyFactor * (ky3 + ky4 - ky1 - ky2) * embeddingY);
  return Polar(phase);
}


// compute the matrix element for on-site two body interaction involving sites on generic sublattic 
//
// s1 = sublattice index for the first creation operator
// s2 = sublattice index for the second annihilation operator
// kx1 = first creation momentum along x on the first sublattice
// ky1 = first creation momentum along y on the first sublattice
// kx2 = second creation momentum along x on the second sublattice
// ky2 = second creation momentum along y on the second sublattice
// kx3 = first annihilation momentum along x on the second sublattice
// ky3 = first annihilation momentum along y on the second sublattice
// kx4 = second annihilation momentum along x on the first sublattice
// ky4 = second annihilation momentum along y on the first sublattice
//
// return value = corresponding matrix element
Complex ParticleOnLatticeHofstadterSingleBandGenericHamiltonian::ComputeEmbeddingForTwoBodyOperator(int s1, int s2, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  double phase = this->TightBindingModel->GetEmbeddingPhase(s2, this->KxFactor * kx3, this->KyFactor * ky3);
  phase += this->TightBindingModel->GetEmbeddingPhase(s1, this->KxFactor * kx4, this->KyFactor * ky4);
  phase -= this->TightBindingModel->GetEmbeddingPhase(s1, this->KxFactor * kx1, this->KyFactor * ky1);
  phase -= this->TightBindingModel->GetEmbeddingPhase(s2, this->KxFactor * kx2, this->KyFactor * ky2);
  return Polar(phase);
}


// compute the matrix element for on-site two body interaction involving sites on generic sublattic 
//
// dRx = number of unit vector translations along x-direction from EncodeSublatticeIndex (translations back to unit cell)
// dRy = number of unit vector translations along y-direction from EncodeSublatticeIndex (translations back to unit cell)
// kx2 = second creation momentum along x for the translated site
// ky2 = second creation momentum along y for the translated site
// kx3 = first annihilation momentum along x for the translated site
// ky3 = first annihilation momentum along y for the translated site
//
// return value = corresponding matrix element
Complex ParticleOnLatticeHofstadterSingleBandGenericHamiltonian::ComputeBlochPhases(int dRx, int dRy, int kx2, int ky2, int kx3, int ky3)
{
  double phase = this->KxFactor * dRx * (kx2-kx3) + this->KyFactor * dRy * (ky2-ky3);
  return Polar(phase);
}

