////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of quatum spin Hall restricted to two band         //
//                                                                            //
//                        last modification : 27/02/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian.h"
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

ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// bandParameter = band parameter
// szSymmetryBreaking = amplitude of the Sz symmetry breaking term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, double bandParameter, double szSymmetryBreaking, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->HamiltonianShift = 0.0;
  this->BandParameter = bandParameter;
  this->SzSymmetryBreaking = szSymmetryBreaking;
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

ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::~ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian()
{
  if (this->InteractionFactorsupupupup != 0)
    {
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsupupupup[i];
	}
      delete[] this->InteractionFactorsupupupup;
    }
  if (this->InteractionFactorsupupdowndown != 0)
    {
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsupupdowndown[i];
	}
      delete[] this->InteractionFactorsupupdowndown;
    }
  if (this->InteractionFactorsdowndownupup != 0)
    {
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsdowndownupup[i];
	}
      delete[] this->InteractionFactorsdowndownupup;
    }
  if (this->InteractionFactorsdowndowndowndown != 0)
    {
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsdowndowndowndown[i];
	}
      delete[] this->InteractionFactorsdowndowndowndown;
    }
  if (this->InteractionFactorsupdownupup != 0)
    {
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsupdownupup[i];
	}
      delete[] this->InteractionFactorsupdownupup;
    }
  if (this->InteractionFactorsupdowndowndown != 0)
    {
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsupdowndowndown[i];
	}
      delete[] this->InteractionFactorsupdowndowndown;
    }
  if (this->InteractionFactorsupupupdown != 0)
    {
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsupupupdown[i];
	}
      delete[] this->InteractionFactorsupupupdown;
    }
  if (this->InteractionFactorsdowndownupdown != 0)
    {
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
 	  delete[] this->InteractionFactorsdowndownupdown[i];
	}
     delete[] this->InteractionFactorsdowndownupdown;
    }
  if (this->InteractionFactorsupdownupdown != 0)
    {
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  delete[] this->InteractionFactorsupdownupdown[i];
	}
      delete[] this->InteractionFactorsupdownupdown;
    }
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
      {
	int Index = (kx1 * this->NbrSiteY) + ky1;
	double d1 = sin (2.0 * M_PI * ((double) kx1) / ((double) this->NbrSiteX));
	double d2 = sin (2.0 * M_PI * ((double) ky1) / ((double) this->NbrSiteY));
	double d3 = (this->BandParameter - cos (2.0 * M_PI * ((double) ky1) / ((double) this->NbrSiteY))
		     - cos (2.0 * M_PI * ((double) kx1) / ((double) this->NbrSiteX)));
	double TmpSzSymmetryBreaking = this->SzSymmetryBreaking * d1;
	HermitianMatrix TmpOneBobyHamiltonian(4, true);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d3);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, Complex(d1, -d2));
	TmpOneBobyHamiltonian.SetMatrixElement(1, 1, -d3);
	TmpOneBobyHamiltonian.SetMatrixElement(2, 2, d3);
	TmpOneBobyHamiltonian.SetMatrixElement(2, 3, Complex(-d1, -d2));
	TmpOneBobyHamiltonian.SetMatrixElement(0, 3, TmpSzSymmetryBreaking);
	TmpOneBobyHamiltonian.SetMatrixElement(1, 2, TmpSzSymmetryBreaking);

	TmpOneBobyHamiltonian.SetMatrixElement(3, 3, -d3);
	ComplexMatrix TmpMatrix(4, 4, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	TmpMatrix[2][2] = 1.0;
	TmpMatrix[3][3] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBobyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	OneBodyBasis[Index] = TmpMatrix;	
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1)  << " " << TmpDiag(2, 2)  << " " << TmpDiag(3, 3) << endl;
      }

  double** CosineTableX = new double*[this->NbrSiteX];
  for (int i = 0; i < this->NbrSiteX; ++i)
    {
      CosineTableX[i] = new double[this->NbrSiteX];
      for (int j = 0; j < this->NbrSiteX; ++j)
	{
	  CosineTableX[i][j] = cos (2.0 * M_PI * ((double) (i - j)) / ((double) this->NbrSiteX));
	}
    }
  double** CosineTableY = new double*[this->NbrSiteY];
  for (int i = 0; i < this->NbrSiteY; ++i)
    {
      CosineTableY[i] = new double[this->NbrSiteY];
      for (int j = 0; j < this->NbrSiteY; ++j)
	{
	  CosineTableY[i][j] = cos (2.0 * M_PI * ((double) (i - j)) / ((double) this->NbrSiteY));
	}
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
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupupupdown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndownupdown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdownupdown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];

      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndowndowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsupupdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndownupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdowndowndown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdownupup[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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
		  this->InteractionFactorsupupupup[i][Index] = -((this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 0, 0) * (CosineTableX[kx2][kx4] + CosineTableY[ky2][ky4]))
								 - (this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 0, 0) * (CosineTableX[kx1][kx4] + CosineTableY[ky1][ky4]))
								 - (this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 0, 0) * (CosineTableX[kx2][kx3] + CosineTableY[ky2][ky3]))
								 + (this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 0, 0) * (CosineTableX[kx1][kx3] + CosineTableY[ky1][ky3])));
		  this->InteractionFactorsdowndowndowndown[i][Index] = -((this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 1, 1) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 1, 1) * (CosineTableX[kx2][kx4] + CosineTableY[ky2][ky4]))
									- (this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 1, 1) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 1, 1) * (CosineTableX[kx1][kx4] + CosineTableY[ky1][ky4]))
									- (this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 1, 1) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 1, 1) * (CosineTableX[kx2][kx3] + CosineTableY[ky2][ky3]))
									+ (this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 1, 1) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 1, 1) * (CosineTableX[kx1][kx3] + CosineTableY[ky1][ky3])));
		  
		  this->InteractionFactorsdowndownupup[i][Index] = -((this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 1, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 1, 0) * (CosineTableX[kx2][kx4] + CosineTableY[ky2][ky4]))
								    - (this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 1, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 1, 0) * (CosineTableX[kx1][kx4] + CosineTableY[ky1][ky4]))
								    - (this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 1, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 1, 0) * (CosineTableX[kx2][kx3] + CosineTableY[ky2][ky3]))
								    + (this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 1, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 1, 0) * (CosineTableX[kx1][kx3] + CosineTableY[ky1][ky3])));
		  this->InteractionFactorsupupdowndown[i][Index] = -((this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 0, 1) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 0, 1) * (CosineTableX[kx2][kx4] + CosineTableY[ky2][ky4]))
								    - (this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 0, 1) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 0, 1) * (CosineTableX[kx1][kx4] + CosineTableY[ky1][ky4]))
								    - (this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 0, 1) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 0, 1) * (CosineTableX[kx2][kx3] + CosineTableY[ky2][ky3]))
								    + (this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 0, 1) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 0, 1) * (CosineTableX[kx1][kx3] + CosineTableY[ky1][ky3])));
		  TotalNbrInteractionFactors += 4;
		  ++Index;
		}
	    }
	  Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
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
		  this->InteractionFactorsupdownupup[i][Index] = -((this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 1, 0) * (CosineTableX[kx2][kx4] + CosineTableY[ky2][ky4]))
								   - (this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 1, 0) * (CosineTableX[kx1][kx4] + CosineTableY[ky1][ky4]))
								   - (this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 1, 0) * (CosineTableX[kx2][kx3] + CosineTableY[ky2][ky3]))
								   + (this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 1, 0) * (CosineTableX[kx1][kx3] + CosineTableY[ky1][ky3])));
		  this->InteractionFactorsupdowndowndown[i][Index] = -((this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 0, 1) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 1, 1) * (CosineTableX[kx2][kx4] + CosineTableY[ky2][ky4]))
								       - (this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 0, 1) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 1, 1) * (CosineTableX[kx1][kx4] + CosineTableY[ky1][ky4]))
								       - (this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 0, 1) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 1, 1) * (CosineTableX[kx2][kx3] + CosineTableY[ky2][ky3]))
								       + (this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 0, 1) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 1, 1) * (CosineTableX[kx1][kx3] + CosineTableY[ky1][ky3])));
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsupupupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndownupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
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
		  this->InteractionFactorsupupupdown[i][Index] = -((this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 0, 1) * (CosineTableX[kx2][kx4] + CosineTableY[ky2][ky4]))
								   - (this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 0, 1) * (CosineTableX[kx1][kx4] + CosineTableY[ky1][ky4]))
								   - (this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 0, 1) * (CosineTableX[kx2][kx3] + CosineTableY[ky2][ky3]))
								   + (this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 0, 1) * (CosineTableX[kx1][kx3] + CosineTableY[ky1][ky3])));
		  this->InteractionFactorsdowndownupdown[i][Index] = -((this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 1, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 1, 1) * (CosineTableX[kx2][kx4] + CosineTableY[ky2][ky4]))
								       - (this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 1, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 1, 1) * (CosineTableX[kx1][kx4] + CosineTableY[ky1][ky4]))
								       - (this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 1, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 1, 1) * (CosineTableX[kx2][kx3] + CosineTableY[ky2][ky3]))
								       + (this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 1, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 1, 1) * (CosineTableX[kx1][kx3] + CosineTableY[ky1][ky3])));
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	  Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
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
		  this->InteractionFactorsupdownupdown[i][Index] = -((this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 1, 1) * (CosineTableX[kx2][kx4] + CosineTableY[ky2][ky4]))
								     - (this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 1, 1) * (CosineTableX[kx1][kx4] + CosineTableY[ky1][ky4]))
								     - (this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 1, 1) * (CosineTableX[kx2][kx3] + CosineTableY[ky2][ky3]))
								     + (this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 1, 1) * (CosineTableX[kx1][kx3] + CosineTableY[ky1][ky3])));
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }

  for (int i = 0; i < this->NbrSiteX; ++i)
    {
      delete[] CosineTableX[i];
    }
  delete[] CosineTableX;
  for (int i = 0; i < this->NbrSiteY; ++i)
    {
      delete[] CosineTableY[i];
    }
  delete[] CosineTableY;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the part of the interaction coefficient coming from the two band truncated basis
//
// basisMatrices = array of basis matrices
// momentumIndex1 = momentum index for the first particle
// momentumIndex2 = momentum index for the second particle
// bandIndex1 = band index of the first particle
// bandIndex2 = band index of the second particle
// return value = coefficient

Complex ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::ComputeBasisContribution(ComplexMatrix* basisMatrices, int momentumIndex1, int momentumIndex2, int bandIndex1, int bandIndex2) 
{
  ComplexMatrix& TmpMatrix1 = basisMatrices[momentumIndex1];
  ComplexMatrix& TmpMatrix2 = basisMatrices[momentumIndex2];
  return Conj(TmpMatrix1[bandIndex1][0]) * TmpMatrix2[bandIndex2][0] 
    + Conj(TmpMatrix1[bandIndex1][1]) * TmpMatrix2[bandIndex2][1]
    + Conj(TmpMatrix1[bandIndex1][2]) * TmpMatrix2[bandIndex2][2]
    + Conj(TmpMatrix1[bandIndex1][3]) * TmpMatrix2[bandIndex2][3];
}
