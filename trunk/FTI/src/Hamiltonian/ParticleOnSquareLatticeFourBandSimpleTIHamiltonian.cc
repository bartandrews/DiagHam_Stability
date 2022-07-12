////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 2d topological insulator based on the simple TI        //
//                       model and restricted to four bands                   //
//                                                                            //
//                        last modification : 27/09/2011                      //
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
#include "Hamiltonian/ParticleOnSquareLatticeFourBandSimpleTIHamiltonian.h"
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

ParticleOnSquareLatticeFourBandSimpleTIHamiltonian::ParticleOnSquareLatticeFourBandSimpleTIHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = strength of the repulsive two body on site interactions
// vPotential = trength of the repulsive two body neareast neighbor interaction
// kineticScale = global energy scale of the kinetic energy term (i.e t1 hopping term)
// mass = mass term of the simple TI model
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnSquareLatticeFourBandSimpleTIHamiltonian::ParticleOnSquareLatticeFourBandSimpleTIHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, double uPotential, double vPotential, 
												       double mass, double gammaX, double gammaY, 
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
  this->Mass = mass;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->WPotential = uPotential;
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

ParticleOnSquareLatticeFourBandSimpleTIHamiltonian::~ParticleOnSquareLatticeFourBandSimpleTIHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnSquareLatticeFourBandSimpleTIHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  int NbrSites = this->NbrSiteX * this->NbrSiteY;
  HermitianMatrix*OneBodyHamiltonian  = new HermitianMatrix [NbrSites];
  this->ComputeOneBodyHamiltonian(OneBodyHamiltonian);

  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  this->OneBodyInteractionFactorsupup = new double [NbrSites];
  this->OneBodyInteractionFactorsumum = new double [NbrSites];
  this->OneBodyInteractionFactorsdpdp = new double [NbrSites];
  this->OneBodyInteractionFactorsdmdm = new double [NbrSites];
  this->OneBodyInteractionFactorsupum = new Complex [NbrSites];
  this->OneBodyInteractionFactorsupdp = new Complex [NbrSites];
  this->OneBodyInteractionFactorsupdm = new Complex [NbrSites];
  this->OneBodyInteractionFactorsumdp = new Complex [NbrSites];
  this->OneBodyInteractionFactorsumdm = new Complex [NbrSites];
  this->OneBodyInteractionFactorsdpdm = new Complex [NbrSites];

  for (int i = 0; i < NbrSites; ++i)
    {
      double Tmp1;
      Complex Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(0, 0, Tmp1);      
      this->OneBodyInteractionFactorsupup[i] = Tmp1;
      OneBodyHamiltonian[i].GetMatrixElement(1, 1, Tmp1);      
      this->OneBodyInteractionFactorsumum[i] = Tmp1;
      OneBodyHamiltonian[i].GetMatrixElement(2, 2, Tmp1);      
      this->OneBodyInteractionFactorsdpdp[i] = Tmp1;
      OneBodyHamiltonian[i].GetMatrixElement(3, 3, Tmp1);      
      this->OneBodyInteractionFactorsdmdm[i] = Tmp1;

      OneBodyHamiltonian[i].GetMatrixElement(1, 0, Tmp2);      
      this->OneBodyInteractionFactorsupum[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(2, 0, Tmp2);      
      this->OneBodyInteractionFactorsupdp[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(3, 0, Tmp2);      
      this->OneBodyInteractionFactorsupdm[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(2, 1, Tmp2);      
      this->OneBodyInteractionFactorsumdp[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(3, 1, Tmp2);      
      this->OneBodyInteractionFactorsumdm[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(3, 2, Tmp2);      
      this->OneBodyInteractionFactorsdpdm[i] = Tmp2;
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
      double FactorAUpADown = Factor * this->VPotential;
      double FactorBUpBDown = Factor * this->VPotential;
      double FactorAUpBDown = Factor * this->WPotential;
      double FactorADownBUp = Factor * this->WPotential;
      if (this->FlatBand == false)
	Factor *= this->UPotential;
      double FactorAUpBUp = Factor;
      double FactorADownBDown = Factor;

      Complex Tmp;

      //  upup upup coefficient
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsumumumum = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdpdpdpdp = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdmdmdmdm = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsumumumum[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdpdpdpdp[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdmdmdmdm[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
                  this->InteractionFactorsupupupup[i][Index] = 0.0;
                  this->InteractionFactorsumumumum[i][Index] = 0.0;
                  this->InteractionFactorsdpdpdpdp[i][Index] = 0.0;
                  this->InteractionFactorsdmdmdmdm[i][Index] = 0.0;		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}



      //  updown updown coefficient
      this->InteractionFactorsupumupum = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsupdpupdp = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsupdmupdm = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsumdpumdp = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsumdmumdm = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsdpdmdpdm = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupumupum[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdpupdp[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdmupdm[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsumdpumdp[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsumdmumdm[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsdpdmdpdm[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
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
                  Tmp = this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupumupum[i][Index] = -2.0 * FactorAUpBUp * Tmp;
                  Tmp = this->ComputeTwoBodyMatrixElementADownBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsdpdmdpdm[i][Index] = -2.0 * FactorADownBDown * Tmp;
                  Tmp = this->ComputeTwoBodyMatrixElementAUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupdmupdm[i][Index] = -2.0 * FactorAUpBDown * Tmp;
                  Tmp = this->ComputeTwoBodyMatrixElementADownBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsumdpumdp[i][Index] = -2.0 * FactorADownBUp * Tmp;
                  Tmp = this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsupdpupdp[i][Index] = -2.0 * FactorAUpADown * Tmp;
                  Tmp = this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
                  this->InteractionFactorsumdmumdm[i][Index] = -2.0 * FactorBUpBDown * Tmp;

		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }

  delete[] OneBodyHamiltonian;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the two body interaction between two sites A and B with up spins
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnSquareLatticeFourBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0 ;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and B with down spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnSquareLatticeFourBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
}

// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnSquareLatticeFourBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementADownBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
}

// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnSquareLatticeFourBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
}

// compute the matrix element for the two body interaction between two sites A with opposite spins 
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on A site with spin down
// ky2 = momentum along y for the creation operator on A site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on A site with spin down
// ky4 = momentum along y for the annihilation operator on A site with spin down
// return value = corresponding matrix element

Complex ParticleOnSquareLatticeFourBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B with opposite spins 
//
// kx1 = momentum along x for the creation operator on B site with spin up
// ky1 = momentum along y for the creation operator on B site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on B site with spin up
// ky3 = momentum along y for the annihilation operator on B site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnSquareLatticeFourBandSimpleTIHamiltonian::ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0;
  return Tmp;
}

// compute the one body hamiltonians related to the band stucture contribution
//
// oneBodyHamiltonians = array of one body hamiltonians

void ParticleOnSquareLatticeFourBandSimpleTIHamiltonian::ComputeOneBodyHamiltonian(HermitianMatrix* oneBodyHamiltonians)
{
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	HermitianMatrix TmpOneBodyHamiltonian(4, true);
	int Index = (kx * this->NbrSiteY) + ky;
	Complex d2 (sin (((double) ky) * this->KyFactor), 0.0);
	double d1 = sin (((double) kx) * this->KxFactor);
	double d3 = (this->Mass - cos (((double) kx) * this->KxFactor) - cos (((double) ky) * this->KyFactor));
	TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d3);
	TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -d3);
	TmpOneBodyHamiltonian.SetMatrixElement(2, 2, d3);
	TmpOneBodyHamiltonian.SetMatrixElement(3, 3, -d3);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 1, d1);
	TmpOneBodyHamiltonian.SetMatrixElement(2, 3, d1);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 3, d2);
	TmpOneBodyHamiltonian.SetMatrixElement(1, 2, d2);
	
	if (this->FlatBand == false)
	  {
	    oneBodyHamiltonians[Index] = TmpOneBodyHamiltonian;
	  }
	else
	  {
	    ComplexMatrix OneBodyBasis(4, 4);
	    OneBodyBasis.SetToIdentity();
	    oneBodyHamiltonians[Index].Copy(TmpOneBodyHamiltonian);
	    RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	    TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, OneBodyBasis);
#else
	    TmpOneBodyHamiltonian.Diagonalize(TmpDiag, OneBodyBasis);
#endif   
	    for (int i = 0; i < 4; ++i)
	      {
		cout << TmpDiag(i, i) << " ";
	      }
	    cout << endl;
	    TmpOneBodyHamiltonian.ClearMatrix();
	    double Tmp = 1.0;
	    TmpOneBodyHamiltonian.SetMatrixElement(2, 2, Tmp);
	    TmpOneBodyHamiltonian.SetMatrixElement(3, 3, Tmp);
	    oneBodyHamiltonians[Index] = TmpOneBodyHamiltonian.InvConjugate(OneBodyBasis);
	  }
      }
}
