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
//                           with real matrix elements                        //
//                                                                            //
//                        last modification : 04/05/2020                      //
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
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandRealHamiltonian.h"
#include "Matrix/RealMatrix.h"
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

ParticleOnLatticeQuantumSpinHallTwoBandRealHamiltonian::ParticleOnLatticeQuantumSpinHallTwoBandRealHamiltonian()
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

ParticleOnLatticeQuantumSpinHallTwoBandRealHamiltonian::ParticleOnLatticeQuantumSpinHallTwoBandRealHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
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

ParticleOnLatticeQuantumSpinHallTwoBandRealHamiltonian::~ParticleOnLatticeQuantumSpinHallTwoBandRealHamiltonian()
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

void ParticleOnLatticeQuantumSpinHallTwoBandRealHamiltonian::EvaluateInteractionFactors()
{
  cout << "warning, using dummy method ParticleOnLatticeQuantumSpinHallTwoBandRealHamiltonian::EvaluateInteractionFactors" << endl;
}

// compute the part of the interaction coefficient coming from the two band truncated basis
//
// basisMatrices = array of basis matrices
// momentumIndex1 = momentum index for the first particle
// momentumIndex2 = momentum index for the second particle
// bandIndex1 = band index of the first particle
// bandIndex2 = band index of the second particle
// return value = coefficient

double ParticleOnLatticeQuantumSpinHallTwoBandRealHamiltonian::ComputeBasisContribution(RealMatrix* basisMatrices, int momentumIndex1, int momentumIndex2, int bandIndex1, int bandIndex2) 
{
  RealMatrix& TmpMatrix1 = basisMatrices[momentumIndex1];
  RealMatrix& TmpMatrix2 = basisMatrices[momentumIndex2];
  return ((TmpMatrix1[bandIndex1][0]) * TmpMatrix2[bandIndex2][0] 
	  + (TmpMatrix1[bandIndex1][1]) * TmpMatrix2[bandIndex2][1]
	  + (TmpMatrix1[bandIndex1][2]) * TmpMatrix2[bandIndex2][2]
	  + (TmpMatrix1[bandIndex1][3]) * TmpMatrix2[bandIndex2][3]);
}
