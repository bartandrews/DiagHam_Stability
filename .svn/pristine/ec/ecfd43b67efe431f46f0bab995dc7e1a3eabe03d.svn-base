////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                           class author: Yang-Le Wu                         //
//                                                                            //
//               class of Haldane model with interacting particles            //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 16/08/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"
#include "GeneralTools/StringTools.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


//debugging flag
//#define DEBUG_POLSHMTBDTBH


// default constructor
//

ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = strength of the repulsive two body nearest neighbor interaction
// vPotential = strength of the repulsive two body second nearest neighbor interaction
// wPotential = strength of the repulsive three body nearest neighbor interaction
// sPotential = strength of the repulsive three body next-to-nearest neighbor interaction
// t1 = hoping amplitude between nearest neighbor sites
// t2 = hoping amplitude between next nearest neighbor sites
// t3 = hoping amplitude between next to next nearest neighbor sites
// phi =  Haldane phase on nnn hopping
// mus = sublattice staggered chemical potential 
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
																	 double uPotential, double vPotential, double wPotential, double sPotential,
																		     double t1, double t2, double t3, double phi, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory, bool hermitianFlag)
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
  this->NextNextNNHopping = t3;
  this->HaldanePhase = phi;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->WPotential = wPotential;
  this->SPotential = sPotential;
  this->ThreeBodySpinAnisotropy = 1.0;
  
  this->Architecture = architecture;
  this->Memory = memory;
  this->NBodyValue = 3;

  this->ComputePhaseArray();

  this->FastMultiplicationFlag = false;
  this->FullTwoBodyFlag = false;
  this->S2Hamiltonian = 0;

  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->NbrIntraSectorSums = 0;
  this->NbrInterSectorSums = 0;
  this->MaxNBody = 3;

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
      this->NBodySign[3] = new double[4];
      this->NBodySign[3][0] = -1.0;
      this->NBodySign[3][1] = -1.0;
      this->NBodySign[3][2] = 1.0;
      this->NBodySign[3][3] = 1.0;
    }
  this->NbrSpinSectors[3] = 4;
  this->SpinIndices[3] = new int*[NbrSpinSectors[3]];
  this->SpinIndices[3][0] = new int[6];
  this->SpinIndices[3][1] = new int[6];
  this->SpinIndices[3][2] = new int[6];
  this->SpinIndices[3][3] = new int[6];
  // symmetry between annihilation indices (0,1,2) and creation indices is used in the current code, but abstract class allows for more general terms with non-conservation of spin, also
  this->SpinIndices[3][0][0] = 0;
  this->SpinIndices[3][0][1] = 0;
  this->SpinIndices[3][0][2] = 0;
  this->SpinIndices[3][0][3] = 0;
  this->SpinIndices[3][0][4] = 0;
  this->SpinIndices[3][0][5] = 0;
  this->SpinIndices[3][1][0] = 1;
  this->SpinIndices[3][1][1] = 1;
  this->SpinIndices[3][1][2] = 1;
  this->SpinIndices[3][1][3] = 1;
  this->SpinIndices[3][1][4] = 1;
  this->SpinIndices[3][1][5] = 1;
  
  this->SpinIndices[3][2][0] = 0;
  this->SpinIndices[3][2][1] = 1;
  this->SpinIndices[3][2][2] = 1;
  this->SpinIndices[3][2][3] = 0;
  this->SpinIndices[3][2][4] = 1;
  this->SpinIndices[3][2][5] = 1;

  this->SpinIndices[3][3][0] = 1;
  this->SpinIndices[3][3][1] = 0;
  this->SpinIndices[3][3][2] = 0;
  this->SpinIndices[3][3][3] = 1;
  this->SpinIndices[3][3][4] = 0;
  this->SpinIndices[3][3][5] = 0;
  this->SpinIndicesShort[3] = new int[NbrSpinSectors[3]];
  this->SpinIndicesShort[3][0] = 0x0;
  this->SpinIndicesShort[3][1] = 0x7 | (0x7<<3);
  this->SpinIndicesShort[3][2] = 0x4 | (0x4<<3);
  this->SpinIndicesShort[3][3] = 0x3 | (0x3<<3);
  
  this->NBodyFlags[3] = true;
  this->NbrNBodySpinMomentumSectorSum[3] = new int[NbrSpinSectors[3]];
  this->NbrNBodySpinMomentumSectorIndicesPerSum[3] = new int* [NbrSpinSectors[3]];
  this->NBodySpinMomentumSectorIndicesPerSum[3] = new int** [NbrSpinSectors[3]];
  this->NBodyInteractionFactors[3] = new Complex**[NbrSpinSectors[3]];
 

  this->Architecture = architecture;
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  

  this->Memory = memory;
  this->EvaluateInteractionFactors();

  if (hermitianFlag)
    this->HermitianSymmetrizeInteractionFactors();
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      PrintMemorySize(cout,TmpMemory);
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::~ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian()
{
  delete[] this->XPhaseTable;
  delete[] this->XHalfPhaseTable;
  delete[] this->YPhaseTable;
  delete[] this->YHalfPhaseTable;
}

// evaluate all interaction factors
//   

void ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix** OneBodyBasis = new ComplexMatrix *[2] ;
  OneBodyBasis[0] = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  OneBodyBasis[1] = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
   
  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactorsupup = new double [this->NbrSiteX * this->NbrSiteY];
      this->OneBodyInteractionFactorsdowndown = new double [this->NbrSiteX * this->NbrSiteY];
    }
  this->ComputeOneBodyMatrices(OneBodyBasis);

  // two-body inter-spin terms common for bosons and fermions:
  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
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
      cout << "Fermionic case not implemented, yet" << endl;
    }
  // Bosonic case
  else
    {
      double FactorU = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

      // Two-Body - intra spin terms
      this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
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
		      
		  Complex sumUUp=0.0, sumUDown;
		  sumUUp = (Conj(OneBodyBasis[0][Index1][0][0]) * OneBodyBasis[0][Index3][0][0] * Conj(OneBodyBasis[0][Index2][0][0]) * OneBodyBasis[0][Index4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		  sumUUp += (Conj(OneBodyBasis[0][Index2][0][0]) * OneBodyBasis[0][Index3][0][0] * Conj(OneBodyBasis[0][Index1][0][0]) * OneBodyBasis[0][Index4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		  sumUUp += (Conj(OneBodyBasis[0][Index1][0][0]) * OneBodyBasis[0][Index4][0][0] * Conj(OneBodyBasis[0][Index2][0][0]) * OneBodyBasis[0][Index3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		  sumUUp += (Conj(OneBodyBasis[0][Index2][0][0]) * OneBodyBasis[0][Index4][0][0] * Conj(OneBodyBasis[0][Index1][0][0]) * OneBodyBasis[0][Index3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();

		  sumUUp += (Conj(OneBodyBasis[0][Index1][0][1]) * OneBodyBasis[0][Index3][0][1] * Conj(OneBodyBasis[0][Index2][0][1]) * OneBodyBasis[0][Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  sumUUp += (Conj(OneBodyBasis[0][Index2][0][1]) * OneBodyBasis[0][Index3][0][1] * Conj(OneBodyBasis[0][Index1][0][1]) * OneBodyBasis[0][Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		  sumUUp += (Conj(OneBodyBasis[0][Index1][0][1]) * OneBodyBasis[0][Index4][0][1] * Conj(OneBodyBasis[0][Index2][0][1]) * OneBodyBasis[0][Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		  sumUUp += (Conj(OneBodyBasis[0][Index2][0][1]) * OneBodyBasis[0][Index4][0][1] * Conj(OneBodyBasis[0][Index1][0][1]) * OneBodyBasis[0][Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

		  sumUDown = (Conj(OneBodyBasis[0][Index1][0][0]) * OneBodyBasis[0][Index3][0][0] * Conj(OneBodyBasis[0][Index2][0][0]) * OneBodyBasis[0][Index4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		  sumUDown += (Conj(OneBodyBasis[0][Index2][0][0]) * OneBodyBasis[0][Index3][0][0] * Conj(OneBodyBasis[0][Index1][0][0]) * OneBodyBasis[0][Index4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		  sumUDown += (Conj(OneBodyBasis[0][Index1][0][0]) * OneBodyBasis[0][Index4][0][0] * Conj(OneBodyBasis[0][Index2][0][0]) * OneBodyBasis[0][Index3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		  sumUDown += (Conj(OneBodyBasis[0][Index2][0][0]) * OneBodyBasis[0][Index4][0][0] * Conj(OneBodyBasis[0][Index1][0][0]) * OneBodyBasis[0][Index3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();

		  sumUDown += (Conj(OneBodyBasis[0][Index1][0][1]) * OneBodyBasis[0][Index3][0][1] * Conj(OneBodyBasis[0][Index2][0][1]) * OneBodyBasis[0][Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  sumUDown += (Conj(OneBodyBasis[0][Index2][0][1]) * OneBodyBasis[0][Index3][0][1] * Conj(OneBodyBasis[0][Index1][0][1]) * OneBodyBasis[0][Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		  sumUDown += (Conj(OneBodyBasis[0][Index1][0][1]) * OneBodyBasis[0][Index4][0][1] * Conj(OneBodyBasis[0][Index2][0][1]) * OneBodyBasis[0][Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		  sumUDown += (Conj(OneBodyBasis[0][Index2][0][1]) * OneBodyBasis[0][Index4][0][1] * Conj(OneBodyBasis[0][Index1][0][1]) * OneBodyBasis[0][Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);



		  if (Index3 == Index4)
		    {
		      sumUUp *= 0.5;
		      sumUDown *= 0.5;
		    }
		  if (Index1 == Index2)
		    {
		      sumUUp *= 0.5;
		      sumUDown *= 0.5;
		    }

		  this->InteractionFactorsupup[i][Index] = 2.0 * FactorU * sumUUp ;
		  this->InteractionFactorsdowndown[i][Index] = 2.0 * FactorU * sumUDown ;

		  TotalNbrInteractionFactors++;
		  ++Index;

		}
	    }
	}
	  
      // Two-Body - inter-spin terms

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
		      
		  Complex sumU=0.0;
		  sumU = (Conj(OneBodyBasis[0][Index1][0][0]) * OneBodyBasis[0][Index3][0][0] * Conj(OneBodyBasis[1][Index2][0][0]) * OneBodyBasis[1][Index4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		  sumU += (Conj(OneBodyBasis[0][Index2][0][0]) * OneBodyBasis[0][Index3][0][0] * Conj(OneBodyBasis[1][Index1][0][0]) * OneBodyBasis[1][Index4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		  sumU += (Conj(OneBodyBasis[0][Index1][0][0]) * OneBodyBasis[0][Index4][0][0] * Conj(OneBodyBasis[1][Index2][0][0]) * OneBodyBasis[1][Index3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
		  sumU += (Conj(OneBodyBasis[0][Index2][0][0]) * OneBodyBasis[0][Index4][0][0] * Conj(OneBodyBasis[1][Index1][0][0]) * OneBodyBasis[1][Index3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();

		  sumU += (Conj(OneBodyBasis[0][Index1][0][1]) * OneBodyBasis[0][Index3][0][1] * Conj(OneBodyBasis[1][Index2][0][1]) * OneBodyBasis[1][Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
		  sumU += (Conj(OneBodyBasis[0][Index2][0][1]) * OneBodyBasis[0][Index3][0][1] * Conj(OneBodyBasis[1][Index1][0][1]) * OneBodyBasis[1][Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
		  sumU += (Conj(OneBodyBasis[0][Index1][0][1]) * OneBodyBasis[0][Index4][0][1] * Conj(OneBodyBasis[1][Index2][0][1]) * OneBodyBasis[1][Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
		  sumU += (Conj(OneBodyBasis[0][Index2][0][1]) * OneBodyBasis[0][Index4][0][1] * Conj(OneBodyBasis[1][Index1][0][1]) * OneBodyBasis[1][Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

		  this->InteractionFactorsupdown[i][Index] = 2.0 * FactorU * sumU;

		  TotalNbrInteractionFactors++;
		  ++Index;

		}
	    }
	}
      	



      // Three Body 
      for (int s=0; s<4; ++s)
	{
	  this->NbrNBodySpinMomentumSectorSum[3][s] = this->NbrSiteX * this->NbrSiteY;
	  this->NbrNBodySpinMomentumSectorIndicesPerSum[3][s] = new int[this->NbrNBodySpinMomentumSectorSum[3][s]];
	  this->NBodySpinMomentumSectorIndicesPerSum[3][s] = new int*[this->NbrNBodySpinMomentumSectorSum[3][s]];
	}
      // spin-polarized cases: channel [0] and [1]
      int KxIn[3];
      int KyIn[3];
      int IndexIn[3];
      
      int** Permutations = 0; 
      double* PermutationSign = 0; 
      int NbrPermutations = this->ComputePermutations(Permutations, PermutationSign, 3);

      for (int s=0; s<2; ++s)
	{
	  int * TmpSpinIndices = this->SpinIndices[3][s];
	  NbrNBodySpinMomentumSectorSum[3][s] = this->NbrSiteX * this->NbrSiteY;
	  int TmpNbrNBodySpinMomentumSectorSums = NbrNBodySpinMomentumSectorSum[3][s];
	  int *TmpNbrNBodySpinMomentumSectorIndicesPerSum = this->NbrNBodySpinMomentumSectorIndicesPerSum[3][s];
	  
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
			  ++TmpNbrNBodySpinMomentumSectorIndicesPerSum[(((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY)];    

		      }
	  int ** TmpNBodySpinMomentumSectorIndicesPerSum = this->NBodySpinMomentumSectorIndicesPerSum[3][s];
	  for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSums; ++i)
	    {
	      if (TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]  > 0)
		{
		  TmpNBodySpinMomentumSectorIndicesPerSum[i] = new int[this->NBodyValue * TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]];
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
			    int TmpSum = (((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY);
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

	  Complex* TmpAAAIn = new Complex[TmpLargestSector];
	  Complex* TmpAAAOut = new Complex[TmpLargestSector];
	  Complex* TmpBBBIn = new Complex[TmpLargestSector];
	  Complex* TmpBBBOut = new Complex[TmpLargestSector];

	  double ThreeBodyFactor = this->WPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);

	  this->NBodyInteractionFactors[3][s] = new Complex* [TmpNbrNBodySpinMomentumSectorSums];
	  Complex **TmpNBodyInteractionFactors = this->NBodyInteractionFactors[3][s];

	  for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSums; ++i)
	    {
	      TmpNBodyInteractionFactors[i] = new Complex[TmpNbrNBodySpinMomentumSectorIndicesPerSum[i] * TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]; ++j1)
		{
		  IndexIn[0] = TmpNBodySpinMomentumSectorIndicesPerSum[i][j1 * 3];
		  IndexIn[1] = TmpNBodySpinMomentumSectorIndicesPerSum[i][(j1 * 3) + 1];
		  IndexIn[2] = TmpNBodySpinMomentumSectorIndicesPerSum[i][(j1 * 3) + 2];
		  KxIn[0] = IndexIn[0] / this->NbrSiteY;
		  KyIn[0] = IndexIn[0] % this->NbrSiteY;
		  KxIn[1] = IndexIn[1] / this->NbrSiteY;
		  KyIn[1] = IndexIn[1] % this->NbrSiteY;
		  KxIn[2] = IndexIn[2] / this->NbrSiteY;
		  KyIn[2] = IndexIn[2] % this->NbrSiteY;

		  Complex TmpAAAIn2 = 0.0;
		  Complex TmpAAAOut2 = 0.0;
		  Complex TmpBBBIn2 = 0.0;
		  Complex TmpBBBOut2 = 0.0;


		  for (int l1 = 0; l1 < NbrPermutations; ++l1)
		    {
		      int* TmpPerm = Permutations[l1];
		      Complex SumTmpAAAIn2 = Conj(OneBodyBasis[TmpSpinIndices[3]][IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[TmpSpinIndices[4]][IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[TmpSpinIndices[5]][IndexIn[TmpPerm[2]]][0][0]) ;
		      Complex SumTmpBBBIn2 = Conj(OneBodyBasis[TmpSpinIndices[3]][IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[TmpSpinIndices[4]][IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[TmpSpinIndices[5]][IndexIn[TmpPerm[2]]][0][1]) ;
		      Complex SumTmpAAAOut2 = OneBodyBasis[TmpSpinIndices[0]][IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[TmpSpinIndices[1]][IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[TmpSpinIndices[2]][IndexIn[TmpPerm[2]]][0][0] ;
		      Complex SumTmpBBBOut2 = OneBodyBasis[TmpSpinIndices[0]][IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[TmpSpinIndices[1]][IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[TmpSpinIndices[2]][IndexIn[TmpPerm[2]]][0][1] ;

				  
		      if ((IndexIn[TmpPerm[0]] == IndexIn[TmpPerm[1]]) && (IndexIn[TmpPerm[1]] == IndexIn[TmpPerm[2]]))
			{
			  SumTmpAAAIn2 /= 6.0;
			  SumTmpBBBIn2 /= 6.0;
			  SumTmpAAAOut2 /= 6.0;
			  SumTmpBBBOut2 /= 6.0;
		      
			}
		      else
			{
			  if( (IndexIn[TmpPerm[0]] == IndexIn[TmpPerm[1]])  || (IndexIn[TmpPerm[1]] == IndexIn[TmpPerm[2]])  || (IndexIn[TmpPerm[0]] == IndexIn[TmpPerm[2]]) )
			    {
			      SumTmpAAAIn2 *= 0.5;
			      SumTmpBBBIn2 *= 0.5;
			      SumTmpAAAOut2 *= 0.5;
			      SumTmpBBBOut2 *= 0.5;
			    }
			}
		      TmpAAAIn2 += SumTmpAAAIn2 ;
		      TmpBBBIn2 += SumTmpBBBIn2 ;

		      TmpAAAOut2 += SumTmpAAAOut2 ;
		      TmpBBBOut2 += SumTmpBBBOut2 ;
		    }

		  TmpAAAIn[j1] =  TmpAAAIn2;
		  TmpAAAOut[j1] = TmpAAAOut2;
		  TmpBBBIn[j1] =  TmpBBBIn2;
		  TmpBBBOut[j1] = TmpBBBOut2;


		}

	      for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]; ++j1)
		{
		  for (int j2 = 0; j2 < TmpNbrNBodySpinMomentumSectorIndicesPerSum[i]; ++j2)
		    {
		      TmpNBodyInteractionFactors[i][Index] = 2.0 * ThreeBodyFactor * (TmpAAAIn[j1] * TmpAAAOut[j2]  + TmpBBBIn[j1] * TmpBBBOut[j2]);							    
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
		}	      
	    }
	  delete[] TmpAAAIn;
	  delete[] TmpAAAOut;
	  delete[] TmpBBBIn;
	  delete[] TmpBBBOut;
	  
	} // end of block: spin polarized case


      for (int i=0; i<NbrPermutations; ++i)
	delete [] Permutations[i];
      delete [] Permutations;
      delete [] PermutationSign;

      // terms for partial spin polarization
      // write explicitly for udd [2], and uud [3]
      
      NbrNBodySpinMomentumSectorSum[3][2] = this->NbrSiteX * this->NbrSiteY;
      int TmpNbrNBodySpinMomentumSectorSumsUDD = NbrNBodySpinMomentumSectorSum[3][2];
      int * TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD = this->NbrNBodySpinMomentumSectorIndicesPerSum[3][2];
      
      NbrNBodySpinMomentumSectorSum[3][3] = this->NbrSiteX * this->NbrSiteY;
      int TmpNbrNBodySpinMomentumSectorSumsUUD = NbrNBodySpinMomentumSectorSum[3][3];
      int * TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD = this->NbrNBodySpinMomentumSectorIndicesPerSum[3][3];

      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUDD; ++i)
	{
	  TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[i] = 0;
	  TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[i] = 0;
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
		      ++TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[(((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY)];
		    if (Index1 <= Index2)
		      ++TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[(((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY)];    
		  }
      
      int **TmpNBodySpinMomentumSectorIndicesPerSumUDD = this->NBodySpinMomentumSectorIndicesPerSum[3][2];
      int **TmpNBodySpinMomentumSectorIndicesPerSumUUD = this->NBodySpinMomentumSectorIndicesPerSum[3][3];

      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUDD; ++i)
	{
	  if (TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[i]  > 0)
	    {
	      TmpNBodySpinMomentumSectorIndicesPerSumUDD[i] = new int[this->NBodyValue * TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[i]];
	      TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[i] = 0;
	    }
	  if (TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[i]  > 0)
	    {
	      TmpNBodySpinMomentumSectorIndicesPerSumUUD[i] = new int[this->NBodyValue * TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[i]];
	      TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[i] = 0;
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
			int TmpSum = (((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY);
			TmpNBodySpinMomentumSectorIndicesPerSumUDD[TmpSum][TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[TmpSum] * 3] = Index1;
			TmpNBodySpinMomentumSectorIndicesPerSumUDD[TmpSum][1 + (TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[TmpSum] * 3)] = Index2;
			TmpNBodySpinMomentumSectorIndicesPerSumUDD[TmpSum][2 + (TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[TmpSum] * 3)] = Index3;
			++TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[TmpSum];
		      }
		    if (Index1 <= Index2)
		      {
			int TmpSum = (((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY);
			TmpNBodySpinMomentumSectorIndicesPerSumUUD[TmpSum][TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[TmpSum] * 3] = Index1;
			TmpNBodySpinMomentumSectorIndicesPerSumUUD[TmpSum][1 + (TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[TmpSum] * 3)] = Index2;
			TmpNBodySpinMomentumSectorIndicesPerSumUUD[TmpSum][2 + (TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[TmpSum] * 3)] = Index3;
			++TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[TmpSum];
		      }
		  }


      int TmpLargestSector = 0;
      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUUD; ++i)
	if (TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[i] > TmpLargestSector)
	  TmpLargestSector = TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[i];
      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUDD; ++i)
	if (TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[i] > TmpLargestSector)
	  TmpLargestSector = TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[i];


      Complex* TmpAAAIn = new Complex[TmpLargestSector];
      Complex* TmpAAAOut = new Complex[TmpLargestSector];
      Complex* TmpBBBIn = new Complex[TmpLargestSector];
      Complex* TmpBBBOut = new Complex[TmpLargestSector];

      double ThreeBodyFactor = 3.0*this->ThreeBodySpinAnisotropy*this->WPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);

      // UDD sector

      NbrPermutations = 2;
      Permutations = new int*[2];
      Permutations[0] = new int[3];
      Permutations[0][0]=0;
      Permutations[0][1]=1;
      Permutations[0][2]=2;
      Permutations[1] = new int[3];
      Permutations[1][0]=0;
      Permutations[1][1]=2;
      Permutations[1][2]=1;


      int * TmpSpinIndicesUDD = this->SpinIndices[3][2];
      this->NBodyInteractionFactors[3][2] = new Complex* [TmpNbrNBodySpinMomentumSectorSumsUDD];
      Complex ** TmpNBodyInteractionFactorsUDD = this->NBodyInteractionFactors[3][2];

      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUDD; ++i)
	{
	  TmpNBodyInteractionFactorsUDD[i] = new Complex[TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[i] * TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[i]; ++j1)
	    {
	      IndexIn[0] = TmpNBodySpinMomentumSectorIndicesPerSumUDD[i][j1 * 3];
	      IndexIn[1] = TmpNBodySpinMomentumSectorIndicesPerSumUDD[i][(j1 * 3) + 1];
	      IndexIn[2] = TmpNBodySpinMomentumSectorIndicesPerSumUDD[i][(j1 * 3) + 2];
	      KxIn[0] = IndexIn[0] / this->NbrSiteY;
   KyIn[0] = IndexIn[0] % this->NbrSiteY;
	      KxIn[1] = IndexIn[1] / this->NbrSiteY;
	      KyIn[1] = IndexIn[1] % this->NbrSiteY;
	      KxIn[2] = IndexIn[2] / this->NbrSiteY;
	      KyIn[2] = IndexIn[2] % this->NbrSiteY;

	      Complex TmpAAAIn2 = 0.0;
	      Complex TmpAAAOut2 = 0.0;
	      Complex TmpBBBIn2 = 0.0;
	      Complex TmpBBBOut2 = 0.0;

	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];
		  Complex SumTmpAAAIn2 = Conj(OneBodyBasis[TmpSpinIndicesUDD[3]][IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[TmpSpinIndicesUDD[4]][IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[TmpSpinIndicesUDD[5]][IndexIn[TmpPerm[2]]][0][0]) ;
		  Complex SumTmpBBBIn2 = Conj(OneBodyBasis[TmpSpinIndicesUDD[3]][IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[TmpSpinIndicesUDD[4]][IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[TmpSpinIndicesUDD[5]][IndexIn[TmpPerm[2]]][0][1]) ;
		  Complex SumTmpAAAOut2 = OneBodyBasis[TmpSpinIndicesUDD[0]][IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[TmpSpinIndicesUDD[1]][IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[TmpSpinIndicesUDD[2]][IndexIn[TmpPerm[2]]][0][0] ;
		  Complex SumTmpBBBOut2 = OneBodyBasis[TmpSpinIndicesUDD[0]][IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[TmpSpinIndicesUDD[1]][IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[TmpSpinIndicesUDD[2]][IndexIn[TmpPerm[2]]][0][1] ;

				  
		  if (IndexIn[TmpPerm[1]] == IndexIn[TmpPerm[2]])
		    {
		      SumTmpAAAIn2 *= 0.5;
		      SumTmpBBBIn2 *= 0.5;
		      SumTmpAAAOut2 *= 0.5;
		      SumTmpBBBOut2 *= 0.5;
		    }
		  TmpAAAIn2 += SumTmpAAAIn2 ;
		  TmpBBBIn2 += SumTmpBBBIn2 ;

		  TmpAAAOut2 += SumTmpAAAOut2 ;
		  TmpBBBOut2 += SumTmpBBBOut2 ;
		}

	      TmpAAAIn[j1] =  TmpAAAIn2;
	      TmpAAAOut[j1] = TmpAAAOut2;
	      TmpBBBIn[j1] =  TmpBBBIn2;
	      TmpBBBOut[j1] = TmpBBBOut2;


	    }

	  for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUDD[i]; ++j2)
		{
		  TmpNBodyInteractionFactorsUDD[i][Index] = 2.0 * ThreeBodyFactor * (TmpAAAIn[j1] * TmpAAAOut[j2]  + TmpBBBIn[j1] * TmpBBBOut[j2]);							    
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }	      
	}

      // UUD sector

      Permutations[0][0]=0;
      Permutations[0][1]=1;
      Permutations[0][2]=2;

      Permutations[1][0]=1;
      Permutations[1][1]=0;
      Permutations[1][2]=2;

      int * TmpSpinIndicesUUD = this->SpinIndices[3][3];
      this->NBodyInteractionFactors[3][3] = new Complex* [TmpNbrNBodySpinMomentumSectorSumsUUD];
      Complex ** TmpNBodyInteractionFactorsUUD = this->NBodyInteractionFactors[3][3]; 

      for (int i = 0; i < TmpNbrNBodySpinMomentumSectorSumsUUD; ++i)
	{
	  TmpNBodyInteractionFactorsUUD[i] = new Complex[TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[i] * TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[i]; ++j1)
	    {
	      IndexIn[0] = TmpNBodySpinMomentumSectorIndicesPerSumUUD[i][j1 * 3];
	      IndexIn[1] = TmpNBodySpinMomentumSectorIndicesPerSumUUD[i][(j1 * 3) + 1];
	      IndexIn[2] = TmpNBodySpinMomentumSectorIndicesPerSumUUD[i][(j1 * 3) + 2];
	      KxIn[0] = IndexIn[0] / this->NbrSiteY;
	      KyIn[0] = IndexIn[0] % this->NbrSiteY;
	      KxIn[1] = IndexIn[1] / this->NbrSiteY;
	      KyIn[1] = IndexIn[1] % this->NbrSiteY;
	      KxIn[2] = IndexIn[2] / this->NbrSiteY;
	      KyIn[2] = IndexIn[2] % this->NbrSiteY;

	      Complex TmpAAAIn2 = 0.0;
	      Complex TmpAAAOut2 = 0.0;
	      Complex TmpBBBIn2 = 0.0;
	      Complex TmpBBBOut2 = 0.0;


	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];
		  Complex SumTmpAAAIn2 = Conj(OneBodyBasis[TmpSpinIndicesUUD[3]][IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[TmpSpinIndicesUUD[4]][IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[TmpSpinIndicesUUD[5]][IndexIn[TmpPerm[2]]][0][0]) ;
		  Complex SumTmpBBBIn2 = Conj(OneBodyBasis[TmpSpinIndicesUUD[3]][IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[TmpSpinIndicesUUD[4]][IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[TmpSpinIndicesUUD[5]][IndexIn[TmpPerm[2]]][0][1]) ;
		  Complex SumTmpAAAOut2 = OneBodyBasis[TmpSpinIndicesUUD[0]][IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[TmpSpinIndicesUUD[1]][IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[TmpSpinIndicesUUD[2]][IndexIn[TmpPerm[2]]][0][0] ;
		  Complex SumTmpBBBOut2 = OneBodyBasis[TmpSpinIndicesUUD[0]][IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[TmpSpinIndicesUUD[1]][IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[TmpSpinIndicesUUD[2]][IndexIn[TmpPerm[2]]][0][1] ;

				  
		  if (IndexIn[TmpPerm[0]] == IndexIn[TmpPerm[1]])
		    {
		      SumTmpAAAIn2 *= 0.5;
		      SumTmpBBBIn2 *= 0.5;
		      SumTmpAAAOut2 *= 0.5;
		      SumTmpBBBOut2 *= 0.5;
		    }
		  TmpAAAIn2 += SumTmpAAAIn2 ;
		  TmpBBBIn2 += SumTmpBBBIn2 ;

		  TmpAAAOut2 += SumTmpAAAOut2 ;
		  TmpBBBOut2 += SumTmpBBBOut2 ;
		}

	      TmpAAAIn[j1] =  TmpAAAIn2;
	      TmpAAAOut[j1] = TmpAAAOut2;
	      TmpBBBIn[j1] =  TmpBBBIn2;
	      TmpBBBOut[j1] = TmpBBBOut2;

	    }
	  
	  for (int j1 = 0; j1 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < TmpNbrNBodySpinMomentumSectorIndicesPerSumUUD[i]; ++j2)
		{
		  TmpNBodyInteractionFactorsUUD[i][Index] = 2.0 * ThreeBodyFactor * (TmpAAAIn[j1] * TmpAAAOut[j2]  + TmpBBBIn[j1] * TmpBBBOut[j2]);							    
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }	      
	}
      
      delete[] TmpAAAIn;
      delete[] TmpAAAOut;
      delete[] TmpBBBIn;
      delete[] TmpBBBOut;
      
      for (int i=0; i<NbrPermutations; ++i)
	delete [] Permutations[i];
      delete [] Permutations;

    }


#ifdef DEBUG_POLSHMTBDTBH
  // testing output:

  for (int k=0; k<=MaxNBody; ++k)
    {
      if (NBodyFlags[k]==true)
	{
	  cout << "k="<<k<<endl;
	  for (int i=0; i<NbrSpinSectors[k];++i)
	    {
	      cout << "i="<<i<<endl;
	      for (int m=0; m<this->NbrNBodySpinMomentumSectorSum[k][i]; ++m)
		{
		  cout << "m="<<m<<endl;
		  for (int n=0; n<NbrNBodySpinMomentumSectorIndicesPerSum[k][i][m]; ++n)
		    cout << "NBodyInteractionFactors["<<k<<","<<i<<","<<m<<","<<n<<"]=("
			 <<NBodySpinMomentumSectorIndicesPerSum[k][i][m][3*n]<<","
			 <<NBodySpinMomentumSectorIndicesPerSum[k][i][m][3*n+1]<<","
			 <<NBodySpinMomentumSectorIndicesPerSum[k][i][m][3*n+2]<<"),"
			 <<NBodyInteractionFactors[k][i][m][n]<<endl;
		}
	    }
	}
    }
#endif
  delete [] OneBodyBasis[0];
  delete [] OneBodyBasis[1];
  delete [] OneBodyBasis;
}
  


// compute the matrix element for the two body interaction between two sites A and B 
//
// kx1 = annihilation momentum along x for the B site
// ky1 = annihilation momentum along y for the B site
// kx2 = creation momentum along x for the B site
// ky2 = creation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2)
{
  double dx = ((double)(kx1-kx2)) * this->KxFactor;
  double dy = ((double)(ky1-ky2)) * this->KyFactor;
  Complex Tmp = 1 + Phase(dx + dy) + Phase(dy);
  return Tmp;
}

// compute the matrix element for the two body interaction between two A sites (or two B sites) 
//
// kx1 = annihilation momentum along x for the second site
// ky1 = annihilation momentum along y for the second site
// kx2 = creation momentum along x for the second site
// ky2 = creation momentum along y for the second site
// return value = corresponding matrix element

Complex ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAA(int kx1, int ky1, int kx2, int ky2)
{
  double dx = ((double)(kx1-kx2)) * this->KxFactor;
  double dy = ((double)(ky1-ky2)) * this->KyFactor;
  Complex Tmp = Phase(dx) + Phase(dy) + Phase(dx + dy);
  return Tmp;
}

// compute the matrix element for the three body interaction between one site A and two sites B 
//
// kx2 = annihilation momentum along x for the first B site
// ky2 = annihilation momentum along y for the first B site
// kx3 = annihilation momentum along x for the second B site
// ky3 = annihilation momentum along y for the second B site
// kx5 = creation momentum along x for the first B site
// ky5 = creation momentum along y for the first B site
// kx6 = creation momentum along x for the second B site
// ky6 = creation momentum along y for the second B site
// return value = corresponding matrix element

Complex ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ComputeThreeBodyMatrixElementABB(int kx2, int ky2, int kx3, int ky3, int kx5, int ky5, int kx6, int ky6)
{
  //double dx2 = ((double)(kx2 - kx5)) * this->KxFactor;
  double dx3 = ((double)(kx3 - kx6)) * this->KxFactor;
  double dy2 = ((double)(ky2 - ky5)) * this->KyFactor;
  double dy3 = ((double)(ky3 - ky6)) * this->KyFactor;
  Complex Tmp = Phase(dy2 + dx3 + dy3) + Phase(dx3 + dy3) + Phase(dy3);
  return Tmp;
}

// compute the matrix element for the three body interaction between NNN sites 
//
// kx2 = annihilation momentum along x for the second site
// ky2 = annihilation momentum along y for the second site
// kx3 = annihilation momentum along x for the third site
// ky3 = annihilation momentum along y for the third site
// kx5 = creation momentum along x for the second site
// ky5 = creation momentum along y for the second site
// kx6 = creation momentum along x for the third site
// ky6 = creation momentum along y for the third site
// return value = corresponding matrix element

Complex ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ComputeThreeBodyMatrixElementAAA(int kx2, int ky2, int kx3, int ky3, int kx5, int ky5, int kx6, int ky6)
{
  double dx2 = ((double)(kx2 - kx5)) * this->KxFactor;
  double dx3 = ((double)(kx3 - kx6)) * this->KxFactor;
  double dy2 = ((double)(ky2 - ky5)) * this->KyFactor;
  double dy3 = ((double)(ky3 - ky6)) * this->KyFactor;
  Complex Tmp = Phase(dx2 + dy2) * (Phase(dx3) + Phase(dy3));
  return Tmp;
}

// compute the matrix element for on-site two body interaction involving A sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ComputeTwoBodyMatrixElementOnSiteAA()
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

Complex ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return 1.0;
}

// compute the matrix element for the on-site three body interaction related to sites A
//
// kx1 = first annihilation momentum along x for the first A site
// ky1 = first annihilation momentum along y for the first A site
// kx2 = second annihilation momentum along x for the second A site
// ky2 = second annihilation momentum along y for the second A site
// kx3 = third annihilation momentum along x for the second A site
// ky3 = third annihilation momentum along y for the second A site
// kx4 = first creation momentum along x for the first A site
// ky4 = first creation momentum along y for the first A site
// kx5 = second creation momentum along x for the second A site
// ky5 = second creation momentum along y for the second A site
// kx6 = third creation momentum along x for the second A site
// ky6 = third creation momentum along y for the second A site
// return value = corresponding matrix element

Complex ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the on-site three body interaction related to sites A
//
// kx1 = first annihilation momentum along x for the first A site
// ky1 = first annihilation momentum along y for the first A site
// kx2 = second annihilation momentum along x for the second A site
// ky2 = second annihilation momentum along y for the second A site
// kx3 = third annihilation momentum along x for the second A site
// ky3 = third annihilation momentum along y for the second A site
// kx4 = first creation momentum along x for the first A site
// ky4 = first creation momentum along y for the first A site
// kx5 = second creation momentum along x for the second A site
// ky5 = second creation momentum along y for the second A site
// kx6 = third creation momentum along x for the second A site
// ky6 = third creation momentum along y for the second A site
// return value = corresponding matrix element

Complex ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  //   Complex Tmp = Phase (0.5 * ((((double) (kx6 + kx5 + kx4 - kx3 - kx2 - kx1)) * this->KxFactor)
  // 			      + ((((double) (ky6 + ky5 + ky4 - ky3 - ky2 - ky1)) * this->KyFactor))));
  //   return Tmp;

  //   double dx = ((double)(-kx1-kx2-kx3+kx4+kx5+kx6)) * this->KxFactor;
  //    double dy = ((double)(-ky1-ky2-ky3+ky4+ky5+ky6)) * this->KyFactor;
  //    Complex Tmp = Phase((-dx-2.0*dy)/3.0);
  //return Tmp;
  
  return 1.0;
}


// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = vector of two matrices spin up, spin down : array of one body transformation matrices

void ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ComputeOneBodyMatrices(ComplexMatrix** oneBodyBasis)
{
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      double x=((double)kx + this->GammaX) * this->KxFactor;
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  double y=((double)ky + this->GammaY) * this->KyFactor;
	  int Index = (kx * this->NbrSiteY) + ky;

	  // My Convention
	  // Complex B1 = - this->NNHopping * Complex(1 + cos(x) + cos(y), - sin(x) - sin(y));
	  // Complex B2 = - this->NextNextNNHopping * Complex(cos(x+y)+2*cos(x-y),-sin(x+y) );
	  // double d0 = - 2.0 * this->NextNNHopping * cos(this->HaldanePhase) * (cos(x) + cos(y) + cos(x-y));
	  // double d3 = - 2.0 * this->NextNNHopping * sin(this->HaldanePhase) * (sin(x) - sin(y) - sin(x-y)) + this->MuS;

	  Complex B1 = - this->NNHopping * Complex(1 + cos(x+y) + cos(y), + sin(x+y) + sin(y));
	  Complex B2 = - this->NextNextNNHopping * Complex(2* cos(x) + cos(x+2*y),  sin(x+2*y));
	  double d0 = - 2.0 * this->NextNNHopping * cos(this->HaldanePhase) * (cos(x) + cos(y) + cos(x+y));
	  double d3 = - 2.0 * this->NextNNHopping * sin(this->HaldanePhase) * (sin(x) + sin(y) - sin(x+y)) + this->MuS;

	  HermitianMatrix TmpOneBodyHamiltonianUp(2, true);
	  TmpOneBodyHamiltonianUp.SetMatrixElement(0, 0, d0 + d3);
	  TmpOneBodyHamiltonianUp.SetMatrixElement(0, 1, B1+B2);
	  TmpOneBodyHamiltonianUp.SetMatrixElement(1, 1, d0 - d3);
	  ComplexMatrix TmpMatrix(2, 2, true);
	  TmpMatrix.SetToIdentity();
	  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	  TmpOneBodyHamiltonianUp.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	  TmpOneBodyHamiltonianUp.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	  oneBodyBasis[0][Index] = TmpMatrix;
	  if (this->FlatBand == false)
	    {
	      this->OneBodyInteractionFactorsupup[Index] = 0.5*TmpDiag(0, 0);
	      // The 1/2 factor is needed because this one body term is counted twice somewhere in ChernInsulatorSingleBandHamiltonian because it was written for spin half particles. This factor is also present for the Kagome Lattice
	    }
	  cout << "+ " << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << "  e1=[" << TmpMatrix[0][0] << ", " << TmpMatrix[0][1] << "]  e2=[" << TmpMatrix[1][0] << ", " << TmpMatrix[1][1] << "]" << endl;

	  // down spins:
	  //B1 = - this->NNHopping * Complex(1 + cos(x+y) + cos(y), - sin(x+y) - sin(y));
	  //B2 = - this->NextNextNNHopping * Complex(2* cos(x) + cos(x+2*y),  - sin(x+2*y));
	  d0 = - 2.0 * this->NextNNHopping * cos(- this->HaldanePhase) * (cos(x) + cos(y) + cos(x+y));
	  d3 = - 2.0 * this->NextNNHopping * sin(- this->HaldanePhase) * (sin(x) + sin(y) - sin(x+y)) + this->MuS;

	  HermitianMatrix TmpOneBodyHamiltonianDown(2, true);
	  TmpOneBodyHamiltonianDown.SetMatrixElement(0, 0, d0 + d3);
	  TmpOneBodyHamiltonianDown.SetMatrixElement(0, 1, B1+B2);
	  TmpOneBodyHamiltonianDown.SetMatrixElement(1, 1, d0 - d3);

	  TmpMatrix.SetToIdentity();
	  TmpDiag.ClearMatrix();
#ifdef __LAPACK__
	  TmpOneBodyHamiltonianDown.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	  TmpOneBodyHamiltonianDown.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	  oneBodyBasis[1][Index] = TmpMatrix;	
	  if (this->FlatBand == false)
	    {
	      this->OneBodyInteractionFactorsdowndown[Index] = 0.5*TmpDiag(0, 0);
	      // The 1/2 factor is needed because this one body term is counted twice somewhere in ChernInsulatorSingleBandHamiltonian because it was written for spin half particles. This factor is also present for the Kagome Lattice
	    }
	  cout << "- " << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << "  e1=[" << TmpMatrix[0][0] << ", " << TmpMatrix[0][1] << "]  e2=[" << TmpMatrix[1][0] << ", " << TmpMatrix[1][1] << "]" << endl;

	}
    }
}



// compute all the phase precalculation arrays 
//

void ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian::ComputePhaseArray()
{
  this->XPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteX];
  this->XHalfPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteX];
  this->XPhaseTableShift = this->NBodyValue * this->NbrSiteX;
  for (int i = -this->XPhaseTableShift; i < this->XPhaseTableShift; ++i)
    {
      this->XPhaseTable[this->XPhaseTableShift + i] = Phase(this->KxFactor * ((double) i));
      this->XHalfPhaseTable[this->XPhaseTableShift + i] = Phase(0.5 * this->KxFactor * ((double) i));
    }
  this->YPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteY];
  this->YHalfPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteY];
  this->YPhaseTableShift = this->NBodyValue * this->NbrSiteY;
  for (int i = -this->YPhaseTableShift; i < this->YPhaseTableShift; ++i)
    {
      this->YPhaseTable[this->YPhaseTableShift + i] = Phase(this->KyFactor * ((double) i));
      this->YHalfPhaseTable[this->YPhaseTableShift + i] = Phase(0.5 * this->KyFactor * ((double) i));
    }
}
 
