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
//                      using two decoupled checkerboard models               //
//                                                                            //
//                        last modification : 04/04/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian.h"
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
// wPotential = strength of the repulsive two body neareast neighbor interaction between opposite spins
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// t2p = hoping amplitude between second next neareast neighbor sites
// mus = sublattice staggered chemical potential 
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian::ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
																		 int nbrSiteY, double uPotential, double vPotential, double wPotential, double t1, double t2, double t2p, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->HamiltonianShift = 0.0;
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->SecondNextNNHoping = t2p;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->WPotential = wPotential;
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

ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian::~ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian::EvaluateInteractionFactors()
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

  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	HermitianMatrix TmpOneBodyHamiltonian(2, true);
	int Index = (kx * this->NbrSiteY) + ky;
	Complex B1 = 4.0 * this->NNHoping * Complex (cos (1.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (1.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY)) * cos(M_PI * 0.25), 
						     sin (1.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * sin (1.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY)) * sin(M_PI * 0.25));
	double d1 = 4.0 * this->SecondNextNNHoping * cos (2.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (2.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY));
	double d3 = this->MuS + (2.0 * this->NextNNHoping * (cos (2.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX))
							     - cos (2.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY))));
	TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
	TmpOneBodyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	ComplexMatrix TmpMatrix(2, 2, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
 	ComplexMatrix TmpMatrix2(4, 4, true);
	TmpMatrix2[0][0] = TmpMatrix[0][0];
	TmpMatrix2[0][1] = TmpMatrix[0][1];
	TmpMatrix2[2][0] = TmpMatrix[1][0];
	TmpMatrix2[2][1] = TmpMatrix[1][1];
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsupup[Index] = TmpDiag(0, 0);
	  }
// 	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " ";
	
	B1 = 4.0 * this->NNHoping * Complex (cos (1.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (1.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY)) * cos(M_PI * 0.25), 
					     sin (1.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX)) * sin (1.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY)) * sin(M_PI * 0.25));
	d1 = 4.0 * this->SecondNextNNHoping * cos (2.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (2.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY));
	d3 = this->MuS + (2.0 * this->NextNNHoping * (cos (2.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX))
						      - cos (2.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY))));
	TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 1, Conj(B1));
	TmpOneBodyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	TmpMatrix[0][1] = 0.0;
	TmpMatrix[1][0] = 0.0;
#ifdef __LAPACK__
	TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	TmpMatrix2[1][2] = TmpMatrix[0][0];
	TmpMatrix2[1][3] = TmpMatrix[0][1];
	TmpMatrix2[3][2] = TmpMatrix[1][0];
	TmpMatrix2[3][3] = TmpMatrix[1][1];
	OneBodyBasis[Index] = TmpMatrix2;	
// 	cout << OneBodyBasis[Index] << endl;
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsdowndown[Index] = TmpDiag(0, 0);
	  }
// 	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << endl;
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

      double Factor = -2.0 * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorAUpADown = Factor * this->VPotential;
      double FactorAUpBDown = Factor * this->WPotential;
      if (this->FlatBand == false)
	Factor *= this->UPotential;

      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];

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
 		  this->InteractionFactorsupup[i][Index] = Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsupup[i][Index] += Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx3, ky3);

 		  this->InteractionFactorsupup[i][Index] += Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsupup[i][Index] += Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx3, ky3);

 		  this->InteractionFactorsupup[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 2, 0, 2);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 2, 0, 2);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 2, 0, 2);
 		  this->InteractionFactorsupup[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 2, 0, 2);

 		  this->InteractionFactorsupup[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx4, ky4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx3, ky3);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx4, ky4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx3, ky3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx3, ky3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx4, ky4);

		  this->InteractionFactorsupup[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsupup[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx3, ky3);

 		  this->InteractionFactorsupup[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsupup[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx3, ky3);


// downdown downdown coefficient
 		  this->InteractionFactorsdowndown[i][Index] = Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsdowndown[i][Index] += Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx3, ky3);

 		  this->InteractionFactorsdowndown[i][Index] += Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsdowndown[i][Index] += Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx3, ky3);

 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 2, 0, 2);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 2, 0, 2);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 2, 0, 2);
 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 2, 0, 2);

 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx4, ky4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx3, ky3);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx4, ky4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx3, ky3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx3, ky3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx4, ky4);

 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx3, ky3);

 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementAUpBDown(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx3, ky3);


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

 		  this->InteractionFactorsupdown[i][Index] = Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);


 		  this->InteractionFactorsupdown[i][Index] += Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);

 		  this->InteractionFactorsupdown[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 2, 0, 2);

 		  this->InteractionFactorsupdown[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx4, ky4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx3, ky3);

 		  this->InteractionFactorsupdown[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);


		  // 		  this->InteractionFactorsupdown[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsupdown[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 1, 2, 1, 2) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx3, ky3);


		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }

  delete[] OneBodyBasis;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the two body interaction between two sites A and B belonging to the same layer
//
// kx1 = momentum along x for the A site
// ky1 = momentum along y for the A site
// kx2 = momentum along x for the B site
// ky2 = momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian::ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp = 2.0 * (cos (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) - ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY))))) 
		       + cos (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) + ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY))))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B with different layer indices 
//
// kx1 = momentum along x for the B site
// ky1 = momentum along y for the B site
// kx2 = momentum along x for the B site
// ky2 = momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian::ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp = Phase (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) + ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY)))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and B belonging to different layers
//
// kx1 = momentum along x for the A site
// ky1 = momentum along y for the A site
// kx2 = momentum along x for the B site
// ky2 = momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonian::ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp = 2.0 * (cos (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) - ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY))))) 
		       + cos (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) + ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY))))));
  return Tmp;
}

