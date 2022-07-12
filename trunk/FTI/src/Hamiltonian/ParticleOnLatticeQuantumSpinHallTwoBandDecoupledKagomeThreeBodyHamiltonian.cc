////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                class of quatum spin Hall restricted to two bands           //
//           using two decoupled kagome models and three body interaction     //
//                                                                            //
//                        last modification : 26/10/2013                      //
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
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
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


// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = strength of the repulsive two body neareast neighbor interaction
// vPotential = strength of the repulsive on site two body interaction between opposite spins
// tightBindingModel = pointer to the tight binding model of the first copy
// flatBandFlag = use flat band model
// timeReversalFlag = apply thge time reversal symmetry on the second copy of the tight binding model  (must be in Bloch form)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian::ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
																		       int nbrSiteY, double uPotential, double vPotential, 
																		       Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, bool timeReversalFlag,  
																		       AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->TightBindingModel = tightBindingModel;
  this->HamiltonianShift = 0.0;
  
  this->FlatBand = flatBandFlag;
  this->TimeReversal = timeReversalFlag;  

  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->TwoBodyUPotential = 0.0;
  this->TwoBodyVPotential = 0.0;

  this->Architecture = architecture;
  this->Memory = memory;

  this->FullTwoBodyFlag = false;
  this->S2Hamiltonian = 0;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;

  this->NbrIntraSectorSums = 0;
  this->NbrInterSectorSums = 0;
  this->MaxNBody = 3;


  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = true;
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

ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian::~ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian::EvaluateInteractionFactors()
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

  this->NBodyFlags = new bool [this->MaxNBody + 1];

  this->NbrSpinSectors = new int [this->MaxNBody + 1];

  this->NBodySign = new double*[this->MaxNBody + 1];
  this->SpinIndices = new int** [this->MaxNBody + 1];
  this->SpinIndicesShort = new int* [this->MaxNBody + 1];

  this->NbrNBodySpinMomentumSectorSum = new int* [this->MaxNBody + 1];
  this->NbrNBodySpinMomentumSectorIndicesPerSum = new int**[this->MaxNBody + 1];
  this->NBodySpinMomentumSectorIndicesPerSum = new int***[this->MaxNBody + 1];
  this->NBodyInteractionFactors = new Complex***[this->MaxNBody + 1];

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->NBodyFlags[k] = false;
      this->NbrSpinSectors[k] = 0;
      this->NBodySign[k] = 0;
      this->SpinIndices[k] = 0;
      this->SpinIndicesShort[k] = 0;
      this->NBodyInteractionFactors[k] = 0;
    }
  
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NBodySign[3] = new double[this->MaxNBody + 1];
      this->NBodySign[3][0] = -1.0;
      this->NBodySign[3][1] = -1.0;
      this->NBodySign[3][2] = 1.0;
      this->NBodySign[3][3] = 1.0;
    }

  this->NbrSpinSectors[this->MaxNBody] = 4;
  this->SpinIndices[this->MaxNBody] = new int*[this->NbrSpinSectors[this->MaxNBody]];
  this->SpinIndicesShort[this->MaxNBody] = new int[NbrSpinSectors[this->MaxNBody]];
  this->SpinIndicesShort[this->MaxNBody][0] = 0x0;
  this->SpinIndicesShort[this->MaxNBody][1] = 0x7 | (0x7<<3);
  this->SpinIndicesShort[this->MaxNBody][2] = 0x4 | (0x4<<3);
  this->SpinIndicesShort[this->MaxNBody][3] = 0x3 | (0x3<<3);
  for (int i = 0; i < this->NbrSpinSectors[this->MaxNBody]; ++i)
    {
      int Lim = 2 * this->MaxNBody;
      this->SpinIndices[this->MaxNBody][i] = new int[Lim];
      for (int k = 0; k < Lim; ++k)
	{	  
	  this->SpinIndices[this->MaxNBody][i][k] =  ((this->SpinIndicesShort[this->MaxNBody][i] >> k)) & 1;
	}
    }
  
  this->NBodyFlags[this->MaxNBody] = true;
  this->NbrNBodySpinMomentumSectorSum[this->MaxNBody] = new int[NbrSpinSectors[this->MaxNBody]];
  this->NbrNBodySpinMomentumSectorIndicesPerSum[this->MaxNBody] = new int* [NbrSpinSectors[this->MaxNBody]];
  this->NBodySpinMomentumSectorIndicesPerSum[this->MaxNBody] = new int** [NbrSpinSectors[this->MaxNBody]];
  this->NBodyInteractionFactors[this->MaxNBody] = new Complex**[NbrSpinSectors[this->MaxNBody]];


 
  if (this->FullTwoBodyFlag == true)
    {      
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
	      ++this->NbrInterSectorIndicesPerSum[(((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)];    
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
		int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = (kx1 * this->NbrSiteY) + ky1;
		this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = (kx2 * this->NbrSiteY) + ky2;
		++this->NbrInterSectorIndicesPerSum[TmpSum];    
	      }
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	{
	  cout << "ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian not implemented for fermions" << endl;
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
		    if (Index1 <= Index2)
		      {
			int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
			this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
			this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
			++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		      }
		  }
	  
	  double FactorU = 0.5 * this->TwoBodyUPotential / ((double) (this->NbrSiteX * this->NbrSiteY));
	  double FactorV = this->TwoBodyVPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      
	  this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
	  this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
	  int TmpNbrBands = 6;
	  int TmpHalfNbrBands = TmpNbrBands / 2;
	  
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	      this->InteractionFactorsdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
		      
		      // upup upup coefficient

		      Complex Tmp = 0.0;		      
		      for (int k = 0; k < TmpNbrBands; ++k)
			{
			  Tmp += FactorU * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, k, k, k, k); 
			  Tmp += FactorU * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, k, k, k, k);
			  Tmp += FactorU * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, k, k, k, k);
			  Tmp += FactorU * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, k, k, k, k);
			}		      
		      for (int k = 0; k < TmpHalfNbrBands; ++k)
			{
			  Tmp += FactorV * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, k, k + 3, k, k + 3);
			  Tmp += FactorV * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, k, k + 3, k, k + 3);
			  Tmp += FactorV * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, k, k + 3, k, k + 3);
			  Tmp += FactorV * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, k, k + 3, k, k + 3);
			}
		      if (Index1 == Index2)
			Tmp *= 0.5;
		      if (Index3 == Index4)
			Tmp *= 0.5;
		      this->InteractionFactorsupup[i][Index] = 2.0 * Tmp;

		      // downdown downdown coefficient
		      Tmp = 0.0;
		      for (int k = 0; k < TmpNbrBands; ++k)
			{
			  Tmp += FactorU * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, k, k, k, k); 
			  Tmp += FactorU * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, k, k, k, k);
			  Tmp += FactorU * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, k, k, k, k);
			  Tmp += FactorU * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, k, k, k, k);
			}		      
		      for (int k = 0; k < TmpHalfNbrBands; ++k)
			{
			  Tmp += FactorV * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, k, k + 3, k, k + 3);
			  Tmp += FactorV * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, k, k + 3, k, k + 3);
			  Tmp += FactorV * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, k, k + 3, k, k + 3);
			  Tmp += FactorV * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, k, k + 3, k, k + 3);
			}

		  
		      if (Index1 == Index2)
			Tmp *= 0.5;
		      if (Index3 == Index4)
			Tmp *= 0.5;
		      this->InteractionFactorsdowndown[i][Index] = 2.0 * Tmp;

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
		      
		      Complex Tmp = 0.0;

		      for (int k = 0; k < TmpNbrBands; ++k)
			{
			  Tmp += FactorU * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, k, k, k, k);
			}
		      for (int k = 0; k < TmpHalfNbrBands; ++k)
			{
			  Tmp += FactorV * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, k, k + 3, k, k + 3);
			}
		      this->InteractionFactorsupdown[i][Index] = 2.0 * Tmp;
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
		}
	    }
	}
    }
  else
    {
      this->NbrIntraSectorSums = 0;
      this->NbrInterSectorSums = 0;
    }
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      cout << "ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian not implemented for fermions" << endl;
    }
  else
    {
      for (int s = 0; s < this->NbrSpinSectors[this->MaxNBody]; ++s)
	{
	  this->NbrNBodySpinMomentumSectorSum[this->MaxNBody][s] = this->NbrSiteX * this->NbrSiteY;
	  this->NbrNBodySpinMomentumSectorIndicesPerSum[this->MaxNBody][s] = new int[this->NbrNBodySpinMomentumSectorSum[this->MaxNBody][s]];
	  this->NBodySpinMomentumSectorIndicesPerSum[this->MaxNBody][s] = new int*[this->NbrNBodySpinMomentumSectorSum[this->MaxNBody][s]];
	}
      int** Permutations = 0; 
      double* PermutationSign = 0; 
      int NbrPermutations = this->ComputePermutations(Permutations, PermutationSign, this->MaxNBody);
      double ThreeBodyUFactor = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      double ThreeBodyVFactor = 3.0 * this->VPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      int KxIn[3];
      int KyIn[3];
      int IndexIn[3];

      // up-up-up and down-down-down contribution
      for (int s = 0; s < 2; ++s)
	{
	  int* TmpSpinIndices = this->SpinIndices[this->MaxNBody][s];
	  cout << "Spin indices : " << TmpSpinIndices[0] << " "  << TmpSpinIndices[1] << " "  << TmpSpinIndices[2] << " " 
	       << TmpSpinIndices[3] << " "  << TmpSpinIndices[4] << " "  << TmpSpinIndices[5] << endl;
 	  NbrNBodySpinMomentumSectorSum[this->MaxNBody][s] = this->NbrSiteX * this->NbrSiteY;
	  int TmpNbrNBodySpinMomentumSectorSums = this->NbrNBodySpinMomentumSectorSum[this->MaxNBody][s];
	  int* TmpNbrNBodySpinMomentumSectorIndicesPerSum = this->NbrNBodySpinMomentumSectorIndicesPerSum[this->MaxNBody][s];

	  for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSums; ++i)
	    TmpNbrNBodySpinMomentumSectorIndicesPerSum[i] = 0;
	  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	      for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
		for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
		  for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		    for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
		      {
			int Index1 = (kx1 * this->NbrSiteY) + ky1;
			int Index2 = (kx2 * this->NbrSiteY) + ky2;
			int Index3 = (kx3 * this->NbrSiteY) + ky3;
			if ((Index1 <= Index2) && (Index2 <= Index3))
			  ++TmpNbrNBodySpinMomentumSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndex(kx1 + kx2 + kx3, ky1 + ky2 + ky3)];    

		      }
	  int ** TmpNBodySpinMomentumSectorIndicesPerSum = this->NBodySpinMomentumSectorIndicesPerSum[this->MaxNBody][s];
	  for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSums; ++i)
	    {
	      if (TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]  > 0)
		{
		  TmpNBodySpinMomentumSectorIndicesPerSum[i] = new int[this->MaxNBody * TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]];
		  TmpNbrNBodySpinMomentumSectorIndicesPerSum[i] = 0;
		}
	    }
	  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	      for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
		for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
		  for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		    for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
		      {
			int Index1 = (kx1 * this->NbrSiteY) + ky1;
			int Index2 = (kx2 * this->NbrSiteY) + ky2;
			int Index3 = (kx3 * this->NbrSiteY) + ky3;
			if ((Index1 <= Index2) && (Index2 <= Index3))
			  {
			    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex(kx1 + kx2 + kx3, ky1 + ky2 + ky3);
			    TmpNBodySpinMomentumSectorIndicesPerSum[TmpSum][TmpNbrNBodySpinMomentumSectorIndicesPerSum[TmpSum] * 3] = Index1;
			    TmpNBodySpinMomentumSectorIndicesPerSum[TmpSum][1 + (TmpNbrNBodySpinMomentumSectorIndicesPerSum[TmpSum] * 3)] = Index2;
			    TmpNBodySpinMomentumSectorIndicesPerSum[TmpSum][2 + (TmpNbrNBodySpinMomentumSectorIndicesPerSum[TmpSum] * 3)] = Index3;
			    ++TmpNbrNBodySpinMomentumSectorIndicesPerSum[TmpSum];
			  }
		      }
	  
	  int TmpLargestSector = 0;
	  for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSums; ++i)
	    if (TmpNbrNBodySpinMomentumSectorIndicesPerSum[i] > TmpLargestSector)
	      TmpLargestSector = TmpNbrNBodySpinMomentumSectorIndicesPerSum[i];

	  Complex** TmpAAAIn = new Complex*[TmpLargestSector];
	  Complex** TmpAAAOut = new Complex*[TmpLargestSector];
	  for (int k = 0; k < TmpLargestSector; ++k)
	    {
	      TmpAAAIn[k] = new Complex[6];
	      TmpAAAOut[k] = new Complex[6];
	    }

	  this->NBodyInteractionFactors[this->MaxNBody][s] = new Complex* [TmpNbrNBodySpinMomentumSectorSums];
	  Complex** TmpNBodyInteractionFactors = this->NBodyInteractionFactors[this->MaxNBody][s];

	  for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSums; ++i)
	    {
	      TmpNBodyInteractionFactors[i] = new Complex[TmpNbrNBodySpinMomentumSectorIndicesPerSum[i] * TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]; ++j1)
		{
		  IndexIn[0] = TmpNBodySpinMomentumSectorIndicesPerSum[i][j1 * this->MaxNBody];
		  IndexIn[1] = TmpNBodySpinMomentumSectorIndicesPerSum[i][(j1 * this->MaxNBody) + 1];
		  IndexIn[2] = TmpNBodySpinMomentumSectorIndicesPerSum[i][(j1 * this->MaxNBody) + 2];
		  KxIn[0] = IndexIn[0] / this->NbrSiteY;
		  KyIn[0] = IndexIn[0] % this->NbrSiteY;
		  KxIn[1] = IndexIn[1] / this->NbrSiteY;
		  KyIn[1] = IndexIn[1] % this->NbrSiteY;
		  KxIn[2] = IndexIn[2] / this->NbrSiteY;
		  KyIn[2] = IndexIn[2] % this->NbrSiteY;

		  Complex TmpAAAIn2[6];
		  Complex TmpAAAOut2[6];
		  for (int k = 0; k < 6; ++k)
		    {
		      TmpAAAIn2[k] = 0.0;
		      TmpAAAOut2[k] = 0.0;
		    }
		  double SymmetryFactor = 1.0;
		  if ((IndexIn[0] == IndexIn[1]) && (IndexIn[0] == IndexIn[2]))
		    {
		      SymmetryFactor = 1.0 / 6.0;
		    }
		  else
		    {
		      if ((IndexIn[0] == IndexIn[1]) || (IndexIn[0] == IndexIn[2]) || (IndexIn[1] == IndexIn[2]))
			{
			  SymmetryFactor = 0.5;
			}
		    }
		  for (int l1 = 0; l1 < NbrPermutations; ++l1)
		    {
		      int* TmpPerm = Permutations[l1];
		      for (int k = 0; k < 6; ++k)
			{
			  TmpAAAIn2[k] += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][s][k] * OneBodyBasis[IndexIn[TmpPerm[1]]][s][k] * OneBodyBasis[IndexIn[TmpPerm[2]]][s][k]);
			  TmpAAAOut2[k] += OneBodyBasis[IndexIn[TmpPerm[0]]][s][k] * OneBodyBasis[IndexIn[TmpPerm[1]]][s][k] * OneBodyBasis[IndexIn[TmpPerm[2]]][s][k];
			}		  
		    }
		  for (int k = 0; k < 6; ++k)
		    {
		      TmpAAAIn[j1][k] =  TmpAAAIn2[k] * SymmetryFactor;
		      TmpAAAOut[j1][k] = TmpAAAOut2[k] * SymmetryFactor;
		    }
		}

	      for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]; ++j1)
		{
		  for (int j2 = 0; j2 < TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]; ++j2)
		    {
		      TmpNBodyInteractionFactors[i][Index] = 0.0;
		      for (int k = 0; k < 6; ++k)
			{
			  TmpNBodyInteractionFactors[i][Index] += 2.0 * ThreeBodyUFactor * TmpAAAIn[j1][k] * TmpAAAOut[j2][k];
			}
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
		}	      
	    }
	  for (int k = 0; k < TmpLargestSector; ++k)
	    {
	      delete[] TmpAAAIn[k];
	      delete[] TmpAAAOut[k];
	    }
	  delete[] TmpAAAIn;
	  delete[] TmpAAAOut;
	} 
      for (int i = 0; i < NbrPermutations; ++i)
	delete[] Permutations[i];
      delete[] Permutations;
      delete[] PermutationSign;

      // up-up-down and up-down-down contribution

      this->NbrNBodySpinMomentumSectorSum[3][2] = this->NbrSiteX * this->NbrSiteY;
      int TmpNbrNBodySpinMomentumSectorSumsUpDownDown = this->NbrNBodySpinMomentumSectorSum[3][2];
      int * TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown = this->NbrNBodySpinMomentumSectorIndicesPerSum[3][2];
      
      this->NbrNBodySpinMomentumSectorSum[3][3] = this->NbrSiteX * this->NbrSiteY;
      int TmpNbrNBodySpinMomentumSectorSumsUpUpDown = this->NbrNBodySpinMomentumSectorSum[3][3];
      int * TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown = this->NbrNBodySpinMomentumSectorIndicesPerSum[3][3];

      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUpDownDown; ++i)
	{
	  TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[i] = 0;
	  TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[i] = 0;
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
	    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
		  {
		    int Index1 = (kx1 * this->NbrSiteY) + ky1;
		    int Index2 = (kx2 * this->NbrSiteY) + ky2;
		    int Index3 = (kx3 * this->NbrSiteY) + ky3;
		    if (Index2 <= Index3)
		      ++TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[this->TightBindingModel->GetLinearizedMomentumIndex(kx1 + kx2 + kx3, ky1 + ky2 + ky3)];
		    if (Index1 <= Index2)
		      ++TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[this->TightBindingModel->GetLinearizedMomentumIndex(kx1 + kx2 + kx3, ky1 + ky2 + ky3)];    
		  }
      
      int** TmpNBodySpinMomentumSectorIndicesPerSumUpDownDown = this->NBodySpinMomentumSectorIndicesPerSum[3][2];
      int** TmpNBodySpinMomentumSectorIndicesPerSumUpUpDown = this->NBodySpinMomentumSectorIndicesPerSum[3][3];

      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUpDownDown; ++i)
	{
	  if (TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[i]  > 0)
	    {
	      TmpNBodySpinMomentumSectorIndicesPerSumUpDownDown[i] = new int[this->MaxNBody * TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[i]];
	      TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[i] = 0;
	    }
	  if (TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[i]  > 0)
	    {
	      TmpNBodySpinMomentumSectorIndicesPerSumUpUpDown[i] = new int[this->MaxNBody * TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[i]];
	      TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
	    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
		  {
		    int Index1 = (kx1 * this->NbrSiteY) + ky1;
		    int Index2 = (kx2 * this->NbrSiteY) + ky2;
		    int Index3 = (kx3 * this->NbrSiteY) + ky3;
		    if (Index2 <= Index3)
		      {
			int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex(kx1 + kx2 + kx3, ky1 + ky2 + ky3);
			TmpNBodySpinMomentumSectorIndicesPerSumUpDownDown[TmpSum][TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[TmpSum] * 3] = Index1;
			TmpNBodySpinMomentumSectorIndicesPerSumUpDownDown[TmpSum][1 + (TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[TmpSum] * 3)] = Index2;
			TmpNBodySpinMomentumSectorIndicesPerSumUpDownDown[TmpSum][2 + (TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[TmpSum] * 3)] = Index3;
			++TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[TmpSum];
		      }
		    if (Index1 <= Index2)
		      {
			int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex(kx1 + kx2 + kx3, ky1 + ky2 + ky3);
			TmpNBodySpinMomentumSectorIndicesPerSumUpUpDown[TmpSum][TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[TmpSum] * 3] = Index1;
			TmpNBodySpinMomentumSectorIndicesPerSumUpUpDown[TmpSum][1 + (TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[TmpSum] * 3)] = Index2;
			TmpNBodySpinMomentumSectorIndicesPerSumUpUpDown[TmpSum][2 + (TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[TmpSum] * 3)] = Index3;
			++TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[TmpSum];
		      }
		  }

      int TmpLargestSector = 0;
      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUpUpDown; ++i)
	if (TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[i] > TmpLargestSector)
	  TmpLargestSector = TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[i];
      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUpDownDown; ++i)
	if (TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[i] > TmpLargestSector)
	  TmpLargestSector = TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[i];
      
      
      Complex** TmpAAAIn = new Complex*[TmpLargestSector];
      Complex** TmpAAAOut = new Complex*[TmpLargestSector];
      for (int k = 0; k < TmpLargestSector; ++k)
	{
	  TmpAAAIn[k] = new Complex[3];
	  TmpAAAOut[k] = new Complex[3];
	}

      // up-down-down contribution
      NbrPermutations = 2;
      Permutations = new int*[2];
      Permutations[0] = new int[3];
      Permutations[0][0] = 0;
      Permutations[0][1] = 1;
      Permutations[0][2] = 2;
      Permutations[1] = new int[3];
      Permutations[1][0] = 0;
      Permutations[1][1] = 2;
      Permutations[1][2] = 1;


      int* TmpSpinIndicesUpDownDown = this->SpinIndices[3][2];
      this->NBodyInteractionFactors[3][2] = new Complex* [TmpNbrNBodySpinMomentumSectorSumsUpDownDown];
      Complex** TmpNBodyInteractionFactorsUpDownDown = this->NBodyInteractionFactors[3][2];

      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUpDownDown; ++i)
	{
	  TmpNBodyInteractionFactorsUpDownDown[i] = new Complex[TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[i] * TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[i]; ++j1)
	    {
	      IndexIn[0] = TmpNBodySpinMomentumSectorIndicesPerSumUpDownDown[i][j1 * 3];
	      IndexIn[1] = TmpNBodySpinMomentumSectorIndicesPerSumUpDownDown[i][(j1 * 3) + 1];
	      IndexIn[2] = TmpNBodySpinMomentumSectorIndicesPerSumUpDownDown[i][(j1 * 3) + 2];
	      KxIn[0] = IndexIn[0] / this->NbrSiteY;
	      KyIn[0] = IndexIn[0] % this->NbrSiteY;
	      KxIn[1] = IndexIn[1] / this->NbrSiteY;
	      KyIn[1] = IndexIn[1] % this->NbrSiteY;
	      KxIn[2] = IndexIn[2] / this->NbrSiteY;
	      KyIn[2] = IndexIn[2] % this->NbrSiteY;

	      Complex TmpAAAIn2[3];
	      Complex TmpAAAOut2[3];
	      for (int k = 0; k < 3; ++k)
		{
		  TmpAAAIn2[k] = 0.0;
		  TmpAAAOut2[k] = 0.0;
		}
	      double SymmetryFactor = 1.0;
	      if (IndexIn[1] == IndexIn[2])
		{
		  SymmetryFactor = 0.5;
		}
	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];
		  for (int k = 0; k < 3; ++k)
		    {
		      TmpAAAIn2[k] += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][k] * 
					   OneBodyBasis[IndexIn[TmpPerm[1]]][1][k + 3] * 
					   OneBodyBasis[IndexIn[TmpPerm[2]]][1][k + 3]);
		      TmpAAAOut2[k] += (OneBodyBasis[IndexIn[TmpPerm[0]]][0][k] * 
					OneBodyBasis[IndexIn[TmpPerm[1]]][1][k + 3] * 
					OneBodyBasis[IndexIn[TmpPerm[2]]][1][k + 3]);
		    }
		}
	      for (int k = 0; k < 3; ++k)
		{
		  TmpAAAIn[j1][k] =  SymmetryFactor * TmpAAAIn2[k];
		  TmpAAAOut[j1][k] = SymmetryFactor * TmpAAAOut2[k];
		}
	    }

	  for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUpDownDown[i]; ++j2)
		{
		  TmpNBodyInteractionFactorsUpDownDown[i][Index] = 0.0;
		  for (int k = 0; k < 3; ++k)
		    {		  
		      TmpNBodyInteractionFactorsUpDownDown[i][Index] += 2.0 * ThreeBodyVFactor * (TmpAAAIn[j1][k] * TmpAAAOut[j2][k]);
		    }							    
		  cout << "upupdown :" << i << " " << j1 << " " << j2  << " " << TotalNbrInteractionFactors << " : " <<  TmpNBodyInteractionFactorsUpDownDown[i][Index] << endl;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }	      
	}

      // up-up-down contribution

      Permutations[0][0] = 0;
      Permutations[0][1] = 1;
      Permutations[0][2] = 2;

      Permutations[1][0] = 1;
      Permutations[1][1] = 0;
      Permutations[1][2] = 2;

      int* TmpSpinIndicesUpUpDown = this->SpinIndices[3][3];
      cout << "Spin indices : " << TmpSpinIndicesUpUpDown[0] << " "  << TmpSpinIndicesUpUpDown[1] << " "  << TmpSpinIndicesUpUpDown[2] << " " 
	   << TmpSpinIndicesUpUpDown[3] << " "  << TmpSpinIndicesUpUpDown[4] << " "  << TmpSpinIndicesUpUpDown[5] << endl;
      this->NBodyInteractionFactors[3][3] = new Complex* [TmpNbrNBodySpinMomentumSectorSumsUpUpDown];
      Complex** TmpNBodyInteractionFactorsUpUpDown = this->NBodyInteractionFactors[3][3]; 

      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUpUpDown; ++i)
	{
	  TmpNBodyInteractionFactorsUpUpDown[i] = new Complex[TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[i] * TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[i]; ++j1)
	    {
	      IndexIn[0] = TmpNBodySpinMomentumSectorIndicesPerSumUpUpDown[i][j1 * 3];
	      IndexIn[1] = TmpNBodySpinMomentumSectorIndicesPerSumUpUpDown[i][(j1 * 3) + 1];
	      IndexIn[2] = TmpNBodySpinMomentumSectorIndicesPerSumUpUpDown[i][(j1 * 3) + 2];
	      KxIn[0] = IndexIn[0] / this->NbrSiteY;
	      KyIn[0] = IndexIn[0] % this->NbrSiteY;
	      KxIn[1] = IndexIn[1] / this->NbrSiteY;
	      KyIn[1] = IndexIn[1] % this->NbrSiteY;
	      KxIn[2] = IndexIn[2] / this->NbrSiteY;
	      KyIn[2] = IndexIn[2] % this->NbrSiteY;

	      Complex TmpAAAIn2[3];
	      Complex TmpAAAOut2[3];
	      for (int k = 0; k < 3; ++k)
		{
		  TmpAAAIn2[k] = 0.0;
		  TmpAAAOut2[k] = 0.0;
		}
	      double SymmetryFactor = 1.0;
	      if (IndexIn[0] == IndexIn[1])
		{
		  SymmetryFactor = 0.5;
		}
	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];
		  for (int k = 0; k < 3; ++k)
		    {
		      TmpAAAIn2[k] += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][k] * 
					   OneBodyBasis[IndexIn[TmpPerm[1]]][0][k] * 
					   OneBodyBasis[IndexIn[TmpPerm[2]]][1][k + 3]);
		      TmpAAAOut2[k] += (OneBodyBasis[IndexIn[TmpPerm[0]]][0][k] * 
					OneBodyBasis[IndexIn[TmpPerm[1]]][0][k] * 
					OneBodyBasis[IndexIn[TmpPerm[2]]][1][k + 3]);
// 		      TmpAAAIn2[k] += Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][TmpSpinIndicesUpUpDown[3]][k] * 
// 					   OneBodyBasis[IndexIn[TmpPerm[1]]][TmpSpinIndicesUpUpDown[4]][k] * 
// 					   OneBodyBasis[IndexIn[TmpPerm[2]]][TmpSpinIndicesUpUpDown[5]][k + 3]);
// 		      TmpAAAOut2[k] += (OneBodyBasis[IndexIn[TmpPerm[0]]][TmpSpinIndicesUpUpDown[0]][k] * 
// 					OneBodyBasis[IndexIn[TmpPerm[1]]][TmpSpinIndicesUpUpDown[1]][k] * 
// 					OneBodyBasis[IndexIn[TmpPerm[2]]][TmpSpinIndicesUpUpDown[2]][k + 3]);
		    }
		}
	      for (int k = 0; k < 3; ++k)
		{
		  TmpAAAIn[j1][k] =  SymmetryFactor * TmpAAAIn2[k];
		  TmpAAAOut[j1][k] = SymmetryFactor * TmpAAAOut2[k];
		}
	    }
	  
	  for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUpUpDown[i]; ++j2)
		{
		  TmpNBodyInteractionFactorsUpUpDown[i][Index] = 0.0;
		  for (int k = 0; k < 3; ++k)
		    {		  
		      TmpNBodyInteractionFactorsUpUpDown[i][Index] += 2.0 * ThreeBodyVFactor * (TmpAAAIn[j1][k] * TmpAAAOut[j2][k]); 
		    }
		  cout << "upupdown :" << i << " " << j1 << " " << j2  << " " << TotalNbrInteractionFactors << " : " <<  TmpNBodyInteractionFactorsUpUpDown[i][Index] << endl;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }	      
	}
      
      for (int k = 0; k < TmpLargestSector; ++k)
	{
	  delete[] TmpAAAIn[k];
	  delete[] TmpAAAOut[k];
	}
      delete[] TmpAAAIn;
      delete[] TmpAAAOut;
      
      for (int i = 0; i < NbrPermutations; ++i)
	delete[] Permutations[i];
      delete[] Permutations;
      
    }
  
  delete[] OneBodyBasis;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices (the leftmost upper block for the spin up, the rightmsdt lower block for the spin down)

void ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	oneBodyBasis[Index] = ComplexMatrix(6, 6, true);
	for (int i = 0; i < 3; ++i)
	  for (int j = 0; j < 3; ++j)
	    oneBodyBasis[Index][2 * i][j] = this->TightBindingModel->GetOneBodyMatrix(Index)[i][j];
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsupup[Index] = 0.5 * this->TightBindingModel->GetEnergy(0, Index);
	  }
	cout << this->TightBindingModel->GetEnergy(0, Index) << " " << this->TightBindingModel->GetEnergy(1, Index) << " " << this->TightBindingModel->GetEnergy(2, Index) << endl;
		
	int IndexInv = this->TightBindingModel->GetLinearizedMomentumIndex(((this->NbrSiteX - kx) % this->NbrSiteX), ((this->NbrSiteY - ky) % this->NbrSiteY));
	for (int i = 0; i < 3; ++i)
	  for (int j = 0; j < 3; ++j)
	  {
	    if (this->TimeReversal == true)
	      oneBodyBasis[Index][2 * i + 1][3 + j] = Conj (this->TightBindingModel->GetOneBodyMatrix(IndexInv)[i][j]);
	    else
	      oneBodyBasis[Index][2 * i + 1][3 + j] = this->TightBindingModel->GetOneBodyMatrix(Index)[i][j];
	  }
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsdowndown[Index] = 0.5 * this->TightBindingModel->GetEnergy(0, IndexInv);
	  }
	cout << this->TightBindingModel->GetEnergy(0, IndexInv) << " " << this->TightBindingModel->GetEnergy(1, IndexInv) << " " << this->TightBindingModel->GetEnergy(2, IndexInv) << endl;
      }
}
