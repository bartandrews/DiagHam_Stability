////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//          class of kagome lattice model with interacting particles          //
//                                                                            //
//                        last modification : 09/12/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeKagomeLatticeThreeBandHamiltonian.h"
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
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeKagomeLatticeThreeBandHamiltonian::ParticleOnLatticeKagomeLatticeThreeBandHamiltonian(ParticleOnSphereWithSU3Spin* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, double uPotential, 
												       double t1, double t2, double lambda1, double lambda2, double mus, double gammaX, double gammaY, 
												       bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->HamiltonianShift = 0.0;
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->NNSpinOrbit = lambda1;
  this->NextNNSpinOrbit = lambda2;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = 0.5 * uPotential;

  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors11 = 0;
  this->OneBodyInteractionFactors12 = 0;
  this->OneBodyInteractionFactors13 = 0;
  this->OneBodyInteractionFactors22 = 0;
  this->OneBodyInteractionFactors23 = 0;
  this->OneBodyInteractionFactors33 = 0;
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

ParticleOnLatticeKagomeLatticeThreeBandHamiltonian::~ParticleOnLatticeKagomeLatticeThreeBandHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeKagomeLatticeThreeBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  int NbrSites = this->NbrSiteX * this->NbrSiteY;
  HermitianMatrix*OneBodyHamiltonian  = new HermitianMatrix [NbrSites];
  this->ComputeOneBodyHamiltonian(OneBodyHamiltonian);

  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  this->OneBodyInteractionFactors11 = new double [NbrSites];
  this->OneBodyInteractionFactors22 = new double [NbrSites];
  this->OneBodyInteractionFactors33 = new double [NbrSites];
  this->OneBodyInteractionFactors12 = new Complex [NbrSites];
  this->OneBodyInteractionFactors13 = new Complex [NbrSites];
  this->OneBodyInteractionFactors23 = new Complex [NbrSites];

  for (int i = 0; i < NbrSites; ++i)
    {
      double Tmp1;
      Complex Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(0, 0, Tmp1);      
      this->OneBodyInteractionFactors11[i] = Tmp1;
      OneBodyHamiltonian[i].GetMatrixElement(1, 1, Tmp1);      
      this->OneBodyInteractionFactors22[i] = Tmp1;
      OneBodyHamiltonian[i].GetMatrixElement(2, 2, Tmp1);      
      this->OneBodyInteractionFactors33[i] = Tmp1;

      OneBodyHamiltonian[i].GetMatrixElement(0, 1, Tmp2);      
      this->OneBodyInteractionFactors12[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(0, 2, Tmp2);      
      this->OneBodyInteractionFactors13[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(1, 2, Tmp2);      
      this->OneBodyInteractionFactors23[i] = Tmp2;
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
	  ++this->NbrInterSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY))];    
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
	    int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY));
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = ((kx1 * this->NbrSiteY) + ky1);
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = ((kx2 * this->NbrSiteY) + ky2);
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	  }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 < Index2)
		  ++this->NbrIntraSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY))];    
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
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY));
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      
      double Factor = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorAUpADown = Factor;// * this->VPotential;
      double FactorBUpBDown = Factor;// * this->VPotential;
      Factor *= this->UPotential;
      double FactorAUpBUp = Factor;
      double FactorADownBDown = Factor;
      double FactorAUpBDown = Factor;
      double FactorADownBUp = Factor;

      Complex Tmp;

      //  11 11 coefficient
      this->InteractionFactors1111 = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactors2222 = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactors3333 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors1111[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactors2222[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactors3333[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
                  this->InteractionFactors1111[i][Index] = 0.0;
                  this->InteractionFactors2222[i][Index] = 0.0;
                  this->InteractionFactors3333[i][Index] = 0.0;
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}



      //  updown updown coefficient
      this->InteractionFactors1212 = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactors1313 = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactors2323 = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactors1212[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactors1313[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactors2323[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
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
//                   Tmp = this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
//                   this->InteractionFactors1212[i][Index] = -2.0 * FactorAUpBUp * Tmp;
//                   Tmp = this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
//                   this->InteractionFactors2323[i][Index] = -2.0 * FactorADownBUp * Tmp;
//                   Tmp = this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
//                   this->InteractionFactors1313[i][Index] = -2.0 * FactorAUpADown * Tmp;

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
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 <= Index2)
		  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
//      if (this->FlatBand == false)
      FactorU *= this->UPotential;

      this->InteractionFactors1111 = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactors2222 = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactors3333 = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactors1111[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactors2222[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactors3333[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
		  this->InteractionFactors1111[i][Index] = FactorU * this->ComputeTwoBodyMatrixElementOnSiteAA();
 		  this->InteractionFactors1111[i][Index] += FactorU * this->ComputeTwoBodyMatrixElementOnSiteAA();
 		  this->InteractionFactors1111[i][Index] += FactorU * this->ComputeTwoBodyMatrixElementOnSiteAA();
 		  this->InteractionFactors1111[i][Index] += FactorU * this->ComputeTwoBodyMatrixElementOnSiteAA();

 		  this->InteractionFactors2222[i][Index] = FactorU * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors2222[i][Index] += FactorU * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors2222[i][Index] += FactorU * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactors2222[i][Index] += FactorU * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

 		  this->InteractionFactors3333[i][Index] = FactorU * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors3333[i][Index] += FactorU * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors3333[i][Index] += FactorU * this->ComputeTwoBodyMatrixElementOnSiteCC(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactors3333[i][Index] += FactorU * this->ComputeTwoBodyMatrixElementOnSiteCC(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

		  if (Index3 == Index4)
		    {
		      this->InteractionFactors1111[i][Index] *= 0.5;
		      this->InteractionFactors2222[i][Index] *= 0.5;
		      this->InteractionFactors3333[i][Index] *= 0.5;
		    }
		  if (Index1 == Index2)
		    {
		      this->InteractionFactors1111[i][Index] *= 0.5;
		      this->InteractionFactors2222[i][Index] *= 0.5;
		      this->InteractionFactors3333[i][Index] *= 0.5;
		    }
		  this->InteractionFactors1111[i][Index] *= 2.0;
		  this->InteractionFactors2222[i][Index] *= 2.0;
		  this->InteractionFactors3333[i][Index] *= 2.0;
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      this->NbrInterSectorSums = 0;
    }

  delete[] OneBodyHamiltonian;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}


// compute the matrix element for the two body interaction between two sites A and B 
//
// kx1 = creation momentum along x for the B site
// ky1 = creation momentum along y for the B site
// k2a = annihilation momentum along x for the B site
// k2b = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeThreeBandHamiltonian::ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KxFactor * ((double) (kx2 - kx1))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and C 
//
// kx1 = creation momentum along x for the C site
// ky1 = creation momentum along y for the C site
// kx2 = annihilation momentum along x for the C site
// ky2 = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeThreeBandHamiltonian::ComputeTwoBodyMatrixElementAC(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KyFactor * ((double) (ky2 - ky1))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B and C 
//
// kx1 = creation momentum along x for the B site
// ky1 = creation momentum along y for the B site
// kx2 = creation momentum along x for the C site
// ky2 = creation momentum along y for the C site
// kx3 = annihilation momentum along x for the B site
// ky3 = annihilation momentum along y for the B site
// kx4 = annihilation momentum along x for the C site
// ky4 = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeThreeBandHamiltonian::ComputeTwoBodyMatrixElementBC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 2.0 * cos (0.5 * ((this->KxFactor * ((double) (kx3 - kx1))) + (this->KyFactor * ((double) (ky4 - ky2)))));
  return Tmp;
}


// compute the matrix element for on-site two body interaction involving A sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeThreeBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteAA()
{
  return 1.0;
}

// compute the matrix element for on-site two body interaction involving B sites
//
// kx1 = first creation momentum along x for the B site
// ky1 = first creation momentum along y for the B site
// kx2 = second creation momentum along x for the B site
// ky2 = second creation momentum along y for the B site
// kx3 = first annihilation momentum along x for the B site
// ky3 = first annihilation momentum along y for the B site
// kx4 = second annihilation momentum along x for the B site
// ky4 = second annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeThreeBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return Phase(0.5 * this->KxFactor * ((double) (kx4 + kx3 - kx2 -kx1)));
}

// compute the matrix element for on-site two body interaction involving C sites
//
// kx1 = first creation momentum along x for the C site
// ky1 = first creation momentum along y for the C site
// kx2 = second creation momentum along x for the C site
// ky2 = second creation momentum along y for the C site
// kx3 = first annihilation momentum along x for the C site
// ky3 = first annihilation momentum along y for the C site
// kx4 = second annihilation momentum along x for the C site
// ky4 = second annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeKagomeLatticeThreeBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return Phase(0.5 * this->KyFactor * ((double) (ky4 + ky3 - ky2 -ky1)));
}


// compute the one body hamiltonians related to the band stucture contribution
//
// oneBodyHamiltonians = array of one body hamiltonians

void ParticleOnLatticeKagomeLatticeThreeBandHamiltonian::ComputeOneBodyHamiltonian(HermitianMatrix* oneBodyHamiltonians)
{
  double KX, KY;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	KX = 0.5 * this->KxFactor * (((double) kx) + this->GammaX);
	KY = 0.5 * this->KyFactor * (((double) ky) + this->GammaY);
	int Index = (kx * this->NbrSiteY) + ky;
	Complex HAB (-2.0 * this->NNHopping, -2.0 * this->NNSpinOrbit);
	HAB *= cos (KX);
	Complex HAC(-2.0 * this->NNHopping, 2.0 * this->NNSpinOrbit);
	HAC *= cos (KY);
	Complex HBC(-2.0 * this->NNHopping, -2.0 * this->NNSpinOrbit);
	HBC *= cos(KX - KY);

	Complex HAB2 (-2.0 * this->NextNNHopping, 2.0 * this->NextNNSpinOrbit);
	HAB2 *= cos (KX - 2.0 * KY);
	Complex HAC2 (-2.0 * this->NextNNHopping, -2.0 * this->NextNNSpinOrbit);
	HAC2 *= cos (2.0 * KX - KY);
	Complex HBC2 (-2.0 * this->NextNNHopping, 2.0 * this->NextNNSpinOrbit);
	HBC2 *= cos (KX + KY);

	HAB += HAB2;
	HAC += HAC2;
	HBC += HBC2;
		
	HermitianMatrix TmpOneBodyHamiltonian(3, true);
	
	TmpOneBodyHamiltonian.SetMatrixElement(0, 1, HAB);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 2, HAC);
	TmpOneBodyHamiltonian.SetMatrixElement(1, 2, HBC);
	
// 	Atomic limit
// // 	TmpOneBodyHamiltonian.SetMatrixElement(0, 1, 0.0);
// // 	TmpOneBodyHamiltonian.SetMatrixElement(0, 2, 0.0);
// // 	TmpOneBodyHamiltonian.SetMatrixElement(1, 2, 0.0);
// // 	TmpOneBodyHamiltonian.SetMatrixElement(0,0,1.0);
// // 	TmpOneBodyHamiltonian.SetMatrixElement(1,1,-1.0);
// // 	TmpOneBodyHamiltonian.SetMatrixElement(2,2,-1.0);
	
	if (this->FlatBand == false)
	  {
	    oneBodyHamiltonians[Index] = TmpOneBodyHamiltonian;
	  }
	else
	  {
	    ComplexMatrix OneBodyBasis(3, 3);
	    OneBodyBasis.SetToIdentity();
	    oneBodyHamiltonians[Index].Copy(TmpOneBodyHamiltonian);
	    RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	    TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, OneBodyBasis);
#else
	    TmpOneBodyHamiltonian.Diagonalize(TmpDiag, OneBodyBasis);
#endif   
	    for (int i = 0; i < 3; ++i)
	      {
		cout << TmpDiag(i, i) << " ";
	      }
	    cout << endl;
	    TmpOneBodyHamiltonian.ClearMatrix();
	    double Tmp = 1.0;
	    TmpOneBodyHamiltonian.SetMatrixElement(1, 1, Tmp);
	    Tmp = 2.0;
	    TmpOneBodyHamiltonian.SetMatrixElement(2, 2, Tmp);
	    oneBodyHamiltonians[Index] = TmpOneBodyHamiltonian.InvConjugate(OneBodyBasis);
	  }
      }
}
