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
#include "Hamiltonian/ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian.h"
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
ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian::ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, TightBindingModelOFLGenericLatticeWithSymmetry* tightBindingModel, double uPotential, bool flatBandFlag, bool noDispersionFlag, AbstractArchitecture* architecture, long memory)
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
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
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
ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian::ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, TightBindingModelOFLGenericLatticeWithSymmetry* tightBindingModel, RealSymmetricMatrix &uPotentialMatrix, bool flatBandFlag, bool noDispersionFlag, AbstractArchitecture* architecture, long memory)
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
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
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

ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian::~ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian()
{
}




// evaluate all interaction factors
//   

void ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->LocalTightBindingModel->GetNbrStatePerBand()];

  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactorsupup = new double [this->LocalTightBindingModel->GetNbrStatePerBand()];
      this->OneBodyInteractionFactorsdowndown = new double [this->LocalTightBindingModel->GetNbrStatePerBand()];
    }

  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = this->LocalTightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	if (this->FlatBand == false)
	  {
	    if (this->NoDispersionFlag == true)
	      {
		this->OneBodyInteractionFactorsupup[Index] = 0;
		this->OneBodyInteractionFactorsdowndown[Index] = 10;
	      }
	    else
	      {	  
		this->OneBodyInteractionFactorsupup[Index] = this->LocalTightBindingModel->GetEnergy(0, Index); 
		this->OneBodyInteractionFactorsdowndown[Index] = this->LocalTightBindingModel->GetEnergy(1, Index);
	      }
	  }
	OneBodyBasis[Index] = this->LocalTightBindingModel->GetOneBodyMatrix(Index);
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
  int NbrSpinFlavours = this->LocalTightBindingModel->GetNbrSubLatticeFlavours();
  int NbrFlavourSublattices = this->LocalTightBindingModel->GetNbrSubLatticesPerFlavour();

  Complex **** DensityOperatorUpUp = new Complex *** [NbrSpinFlavours];
  Complex **** DensityOperatorDownDown = new Complex *** [NbrSpinFlavours];
  Complex **** DensityOperatorUpDown = new Complex *** [NbrSpinFlavours];
  Complex **** DensityOperatorDownUp = new Complex *** [NbrSpinFlavours];
      
  for (int Spin = 0; Spin < NbrSpinFlavours; ++Spin)	
    {
      DensityOperatorUpUp[Spin] = new Complex ** [this->NbrSiteX* this->NbrSiteY];
      DensityOperatorDownDown[Spin] = new Complex ** [this->NbrSiteX* this->NbrSiteY];
      DensityOperatorUpDown[Spin] = new Complex ** [this->NbrSiteX* this->NbrSiteY];
      DensityOperatorDownUp[Spin] = new Complex ** [this->NbrSiteX* this->NbrSiteY];
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	  {
	    int Index1 = (kx1 * this->NbrSiteY) + ky1;
	    DensityOperatorUpUp[Spin][Index1] = new Complex * [this->NbrSiteX*this->NbrSiteY];
	    DensityOperatorDownDown[Spin][Index1] = new Complex * [this->NbrSiteX*this->NbrSiteY];
	    DensityOperatorUpDown[Spin][Index1] = new Complex * [this->NbrSiteX*this->NbrSiteY];
	    DensityOperatorDownUp[Spin][Index1] = new Complex * [this->NbrSiteX*this->NbrSiteY];
	      
	    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		{
		  int Index2 = (kx2 * this->NbrSiteY) + ky2;
		  DensityOperatorUpUp[Spin][Index1][Index2] = new Complex [this->NbrReciprocalVectors1*this->NbrReciprocalVectors2];
		  DensityOperatorDownDown[Spin][Index1][Index2] = new Complex [this->NbrReciprocalVectors1*this->NbrReciprocalVectors2];
		  DensityOperatorUpDown[Spin][Index1][Index2] = new Complex [this->NbrReciprocalVectors1*this->NbrReciprocalVectors2];
		  DensityOperatorDownUp[Spin][Index1][Index2] = new Complex [this->NbrReciprocalVectors1*this->NbrReciprocalVectors2];

		  for (int Nx = 0; Nx < this->NbrReciprocalVectors1; Nx++)
		    {
		      for (int Ny = 0; Ny < this->NbrReciprocalVectors2; Ny++)
			{
			  int TmpIndex =  Nx * this->NbrReciprocalVectors2 + Ny;
			  
			  DensityOperatorUpUp[Spin][Index1][Index2][TmpIndex] = 0;
			  DensityOperatorDownDown[Spin][Index1][Index2][TmpIndex] = 0;
			  DensityOperatorUpDown[Spin][Index1][Index2][TmpIndex] = 0;
			  DensityOperatorDownUp[Spin][Index1][Index2][TmpIndex] = 0;
			  for (int Subl = 0; Subl < NbrFlavourSublattices; Subl++)
			    {	 
			      for (int Nx1 = 0; Nx1 < this->NbrReciprocalVectors1; Nx1++)
				{
				  for (int Ny1 = 0; Ny1 < this->NbrReciprocalVectors2; Ny1++)
				    {
				      int TotalSubl =  this->LocalTightBindingModel->TotalSublatticeIndex(Spin, Subl);
				      int ReciprocalVectorIndex1 =  this->LocalTightBindingModel->PeriodicLinearizedReciprocalSpaceIndex(Nx1,Ny1,TotalSubl);
				      int ReciprocalVectorIndex2 =  this->LocalTightBindingModel->PeriodicLinearizedReciprocalSpaceIndex(Nx1 - Nx +  this->NbrReciprocalVectors1, Ny1 - Ny +  this->NbrReciprocalVectors2,TotalSubl);
				      
				      DensityOperatorUpUp[Spin][Index1][Index2][TmpIndex] += Conj(OneBodyBasis[Index1][0][ReciprocalVectorIndex1]) * OneBodyBasis[Index2][0][ReciprocalVectorIndex2];
				      DensityOperatorDownDown[Spin][Index1][Index2][TmpIndex] += Conj(OneBodyBasis[Index1][1][ReciprocalVectorIndex1]) * OneBodyBasis[Index2][1][ReciprocalVectorIndex2];
				      DensityOperatorUpDown[Spin][Index1][Index2][TmpIndex] += Conj(OneBodyBasis[Index1][0][ReciprocalVectorIndex1]) * OneBodyBasis[Index2][1][ReciprocalVectorIndex2];
				      DensityOperatorDownUp[Spin][Index1][Index2][TmpIndex] += Conj(OneBodyBasis[Index1][1][ReciprocalVectorIndex1]) * OneBodyBasis[Index2][0][ReciprocalVectorIndex2];
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
      cout <<"Fermionic case not implemented yet"<<endl;
      exit(-6);
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

      if (this->FlatBand == true)
	{
	  if (this->UPotentialMatrix(0,0)!=0.0)
	    this->UPotentialMatrix *= 1.0/this->UPotentialMatrix(0,0);
	  else
	    this->UPotentialMatrix.SetAllEntries(1.0);
	}
      this->UPotentialMatrix *= this->SymmetryMultiplier1 * this->SymmetryMultiplier2 * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
            
      //  upup upup coefficient
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  
	  int Index = 0;

	  for (int j2 = 0; j2 <  this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
		  
		  this->InteractionFactorsupupupup[i][Index] = 0.0;

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
				this->InteractionFactorsupupupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpUp[Spin1][Index1][Index4][TmpIndex] * DensityOperatorUpUp[Spin2][Index2][Index3][TmpIndex2];
				this->InteractionFactorsupupupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpUp[Spin1][Index2][Index4][TmpIndex] * DensityOperatorUpUp[Spin2][Index1][Index3][TmpIndex2];
				this->InteractionFactorsupupupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpUp[Spin1][Index1][Index3][TmpIndex] * DensityOperatorUpUp[Spin2][Index2][Index4][TmpIndex2];
				this->InteractionFactorsupupupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpUp[Spin1][Index2][Index3][TmpIndex] * DensityOperatorUpUp[Spin2][Index1][Index4][TmpIndex2];
			      }
			}
		    }

		  if (Index3 == Index4)
		    {
		      this->InteractionFactorsupupupup[i][Index] *= 0.5;
		    }
		  if (Index1 == Index2)
		    {
		      this->InteractionFactorsupupupup[i][Index] *= 0.5;
		    }
		  //this->InteractionFactorsupupupup[i][Index] *= 2.0;
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

       //  upup downdown coefficient
      
      this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  
	  for (int j2 = 0; j2 <  this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      
	      
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
		  
		  this->InteractionFactorsupupdowndown[i][Index] = 0.0;
		  
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
				this->InteractionFactorsupupdowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpDown[Spin1][Index1][Index4][TmpIndex] * DensityOperatorUpDown[Spin2][Index2][Index3][TmpIndex2];
				this->InteractionFactorsupupdowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpDown[Spin1][Index2][Index4][TmpIndex] * DensityOperatorUpDown[Spin2][Index1][Index3][TmpIndex2];
				
				this->InteractionFactorsupupdowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpDown[Spin1][Index1][Index3][TmpIndex] * DensityOperatorUpDown[Spin2][Index2][Index4][TmpIndex2];
				this->InteractionFactorsupupdowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpDown[Spin1][Index2][Index3][TmpIndex] * DensityOperatorUpDown[Spin2][Index1][Index4][TmpIndex2];
			      }
			}
		    }			  
		  if (Index1 == Index2)
		    this->InteractionFactorsupupdowndown[i][Index] *= 0.5;
		  if (Index3 == Index4)
		    this->InteractionFactorsupupdowndown[i][Index] *= 0.5;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      
      //  downdown downdown coefficient
      this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];

      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndowndowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;

	  for (int j2 = 0; j2 <  this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
	      
		  this->InteractionFactorsdowndowndowndown[i][Index] = 0.0;
		  
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
				this->InteractionFactorsdowndowndowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownDown[Spin1][Index1][Index4][TmpIndex] * DensityOperatorDownDown[Spin2][Index2][Index3][TmpIndex2];
				this->InteractionFactorsdowndowndowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownDown[Spin1][Index2][Index4][TmpIndex] * DensityOperatorDownDown[Spin2][Index1][Index3][TmpIndex2];
				
				this->InteractionFactorsdowndowndowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownDown[Spin1][Index1][Index3][TmpIndex] * DensityOperatorDownDown[Spin2][Index2][Index4][TmpIndex2];
				this->InteractionFactorsdowndowndowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownDown[Spin1][Index2][Index3][TmpIndex] * DensityOperatorDownDown[Spin2][Index1][Index4][TmpIndex2];
			      }
			}
		    }

		  if (Index1 == Index2)
		    this->InteractionFactorsdowndowndowndown[i][Index] *= 0.5;
		  if (Index3 == Index4)
		    this->InteractionFactorsdowndowndowndown[i][Index] *= 0.5;

		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}


      //  downdown upup coefficient
      this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndownupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  
	  for (int j2 = 0; j2 <  this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
		  
		  this->InteractionFactorsdowndownupup[i][Index] = 0.0;
		  
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
				this->InteractionFactorsdowndownupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownUp[Spin1][Index1][Index4][TmpIndex] * DensityOperatorDownUp[Spin2][Index2][Index3][TmpIndex2];
				this->InteractionFactorsdowndownupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownUp[Spin1][Index2][Index4][TmpIndex] * DensityOperatorDownUp[Spin2][Index1][Index3][TmpIndex2];
				
				this->InteractionFactorsdowndownupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownUp[Spin1][Index1][Index3][TmpIndex] * DensityOperatorDownUp[Spin2][Index2][Index4][TmpIndex2];
				this->InteractionFactorsdowndownupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownUp[Spin1][Index2][Index3][TmpIndex] * DensityOperatorDownUp[Spin2][Index1][Index4][TmpIndex2];
			      }
			}
		    }			  
		  if (Index1 == Index2)
		    this->InteractionFactorsdowndownupup[i][Index] *= 0.5;
		  if (Index3 == Index4)
		    this->InteractionFactorsdowndownupup[i][Index] *= 0.5;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}


      this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupup[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;

	  for (int j2 = 0; j2 <  this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
		  
		  this->InteractionFactorsupdownupup[i][Index] = 0.0;		
		  
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
				this->InteractionFactorsupdownupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpUp[Spin1][Index1][Index4][TmpIndex] * DensityOperatorDownUp[Spin2][Index2][Index3][TmpIndex2];
				this->InteractionFactorsupdownupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpUp[Spin1][Index1][Index3][TmpIndex] * DensityOperatorDownUp[Spin2][Index2][Index4][TmpIndex2];
				this->InteractionFactorsupdownupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpUp[Spin1][Index1][Index3][TmpIndex2] * DensityOperatorDownUp[Spin2][Index2][Index4][TmpIndex];
				this->InteractionFactorsupdownupup[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpUp[Spin1][Index1][Index4][TmpIndex2] * DensityOperatorDownUp[Spin2][Index2][Index3][TmpIndex];
			      }
			}
		    }			  
		  
		  if (Index3 == Index4)
		    this->InteractionFactorsupdownupup[i][Index] *= 0.5;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  updown downdown coefficient
      this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupdowndowndown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;

	  for (int j2 = 0; j2 <  this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
		  
		  this->InteractionFactorsupdowndowndown[i][Index] = 0.0;
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
				this->InteractionFactorsupdowndowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpDown[Spin1][Index1][Index4][TmpIndex] * DensityOperatorDownDown[Spin2][Index2][Index3][TmpIndex2];
				this->InteractionFactorsupdowndowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpDown[Spin1][Index1][Index3][TmpIndex] * DensityOperatorDownDown[Spin2][Index2][Index4][TmpIndex2];
				this->InteractionFactorsupdowndowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpDown[Spin1][Index1][Index3][TmpIndex2] * DensityOperatorDownDown[Spin2][Index2][Index4][TmpIndex];
				this->InteractionFactorsupdowndowndown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpDown[Spin1][Index1][Index4][TmpIndex2] * DensityOperatorDownDown[Spin2][Index2][Index3][TmpIndex];
			      }
			}
		    }			  
		  if (Index3 == Index4)
		    this->InteractionFactorsupdowndowndown[i][Index] *= 0.5;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
	  
      
      //  upup updown coefficient
      this->InteractionFactorsupupupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupupupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  
	  
	  for (int j2 = 0; j2 <  this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
		  
		  this->InteractionFactorsupupupdown[i][Index] = 0.0;
		  
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
				this->InteractionFactorsupupupdown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpDown[Spin1][Index1][Index4][TmpIndex] * DensityOperatorUpUp[Spin2][Index2][Index3][TmpIndex2];
				this->InteractionFactorsupupupdown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpDown[Spin1][Index2][Index4][TmpIndex] * DensityOperatorUpUp[Spin2][Index1][Index3][TmpIndex2];

				this->InteractionFactorsupupupdown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpUp[Spin1][Index1][Index3][TmpIndex] * DensityOperatorUpDown[Spin2][Index2][Index4][TmpIndex2];
				this->InteractionFactorsupupupdown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorUpUp[Spin1][Index2][Index3][TmpIndex] * DensityOperatorUpDown[Spin2][Index1][Index4][TmpIndex2];
			      }
			}
		    }			  
		  if (Index1 == Index2)
		     this->InteractionFactorsupupupdown[i][Index] *= 0.5;

		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //  downdown updown coefficient
      this->InteractionFactorsdowndownupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsdowndownupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  
	  for (int j2 = 0; j2 <  this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      
	      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
		  
		  this->InteractionFactorsdowndownupdown[i][Index] = 0;
		  
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
				this->InteractionFactorsdowndownupdown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownDown[Spin1][Index1][Index4][TmpIndex] * DensityOperatorDownUp[Spin2][Index2][Index3][TmpIndex2];
				this->InteractionFactorsdowndownupdown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownDown[Spin1][Index2][Index4][TmpIndex] * DensityOperatorDownUp[Spin2][Index1][Index3][TmpIndex2];
				this->InteractionFactorsdowndownupdown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownUp[Spin1][Index1][Index3][TmpIndex] * DensityOperatorDownDown[Spin2][Index2][Index4][TmpIndex2];
				this->InteractionFactorsdowndownupdown[i][Index] +=  UPotentialMatrix(Spin1,Spin2) * DensityOperatorDownUp[Spin1][Index2][Index3][TmpIndex] * DensityOperatorDownDown[Spin2][Index1][Index4][TmpIndex2];
			      }
			}
		    }			  
		  if (Index1 == Index2)
		     this->InteractionFactorsdowndownupdown[i][Index] *= 0.5;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      
      this->InteractionFactorsupdownupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	     
	      int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      int kx4 = Index4 / this->NbrSiteY;
	      int ky4 = Index4 % this->NbrSiteY;
	      
	      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
		  int kx1 = Index1 / this->NbrSiteY;
		  int ky1 = Index1 % this->NbrSiteY;
		  int kx2 = Index2 / this->NbrSiteY;
		  int ky2 = Index2 % this->NbrSiteY;
		  
		  this->InteractionFactorsupdownupdown[i][Index] = 0;
 
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
				this->InteractionFactorsupdownupdown[i][Index] += UPotentialMatrix(Spin1,Spin2) *   DensityOperatorUpDown[Spin1][Index1][Index4][TmpIndex] * DensityOperatorDownUp[Spin2][Index2][Index3][TmpIndex2];
				this->InteractionFactorsupdownupdown[i][Index] += UPotentialMatrix(Spin1,Spin2) *   DensityOperatorUpUp[Spin1][Index1][Index3][TmpIndex] * DensityOperatorDownDown[Spin2][Index2][Index4][TmpIndex2];
				this->InteractionFactorsupdownupdown[i][Index] += UPotentialMatrix(Spin1,Spin2) *   DensityOperatorUpUp[Spin1][Index1][Index3][TmpIndex2] * DensityOperatorDownDown[Spin2][Index2][Index4][TmpIndex];
				this->InteractionFactorsupdownupdown[i][Index] += UPotentialMatrix(Spin1,Spin2) *   DensityOperatorUpDown[Spin1][Index1][Index4][TmpIndex2] * DensityOperatorDownUp[Spin2][Index2][Index3][TmpIndex];
			      }
			}
		    }
		  
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
		  delete [] DensityOperatorUpUp[Spin][Index1][Index2];
		  delete [] DensityOperatorDownDown[Spin][Index1][Index2];
		  delete [] DensityOperatorUpDown[Spin][Index1][Index2];
		  delete [] DensityOperatorDownUp[Spin][Index1][Index2];
		}
	    delete [] DensityOperatorUpUp[Spin][Index1];
	    delete [] DensityOperatorDownDown[Spin][Index1];
	    delete [] DensityOperatorUpDown[Spin][Index1];
	    delete [] DensityOperatorDownUp[Spin][Index1];
	  }
      delete [] DensityOperatorUpUp[Spin];
      delete [] DensityOperatorDownDown[Spin];
      delete [] DensityOperatorUpDown[Spin];
      delete [] DensityOperatorDownUp[Spin];
    }
  delete [] DensityOperatorUpUp;
  delete [] DensityOperatorDownDown;
  delete [] DensityOperatorUpDown;
  delete [] DensityOperatorDownUp;
  delete [] OneBodyBasis;
}
