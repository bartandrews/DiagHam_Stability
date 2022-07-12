////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//        class of generic hamiltonian for interacting spinful particles      //
//                     on lattice without Sz conservation and                 //
//               written in real space and handling 2d translations           //
//                                                                            //
//                        last modification : 05/08/2015                      //
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
#include "Hamiltonian/ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationS2Hamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;



// default constructor
//

ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian::ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian()
{
  this->XMomentum = 0;
  this->MaxXMomentum = 1;
  this->YMomentum = 0;
  this->MaxYMomentum = 1;
  this->DiagonalElements = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSites = number of sites
// xMomentum = momentum sector in the x direction
// maxXMomentum = number of momentum sectors in the x direction
// yMomentum = momentum sector in the x direction
// maxYMomentum = number of momentum sectors in the x direction
// tightBinding = hamiltonian corresponding to the tight-binding model in real space, orbitals with even indices (resp. odd indices) are considered as spin up (resp. spin down)
// densityDensityupup = matrix that gives the amplitude of each density-density interaction term between particles with spin up
// densityDensitydowndown = matrix that gives the amplitude of each density-density interaction term between particles with spin down
// densityDensityupdown = matrix that gives the amplitude of each density-density interaction term between particles with spin up and down
// sxSx = matrix that gives the amplitude of each Sx_i Sx_j term
// sySy = matrix that gives the amplitude of each Sy_i Sy_j term
// szSz = matrix that gives the amplitude of each Sz_i Sz_j term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian::ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSites, 
																     int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, 
																     HermitianMatrix& tightBinding, 
																     RealSymmetricMatrix& densityDensityupup, 
																     RealSymmetricMatrix& densityDensitydowndown, 
																     RealSymmetricMatrix& densityDensityupdown, 
																     RealSymmetricMatrix& sxSx,
																     RealSymmetricMatrix& sySy, RealSymmetricMatrix& szSz, 
																     AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSites = nbrSites;
  this->XMomentum = xMomentum;
  this->MaxXMomentum = maxXMomentum;
  this->YMomentum = yMomentum;
  this->MaxYMomentum = maxYMomentum;
    
  this->LzMax = this->NbrSites - 1;
  this->HamiltonianShift = 0.0;
  this->DiagonalElements = 0;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = true;//true;
  
  this->EvaluateExponentialFactors();
  this->EvaluateOneBodyFactorsFromTightBingding(tightBinding);
  
  this->EvaluateInteractionFactorsFromDensityDensityAndHeisenberg(densityDensityupup, densityDensitydowndown, densityDensityupdown, sxSx, sySy, szSz);
  
//   cout << this->InteractionFactorsupup[0][0] << endl;
    
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      cout << TmpMemory << endl;
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


// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSites = number of sites
// xMomentum = momentum sector in the x direction
// maxXMomentum = number of momentum sectors in the x direction
// yMomentum = momentum sector in the x direction
// maxYMomentum = number of momentum sectors in the x direction
// tightBinding = hamiltonian corresponding to the tight-binding model in real space, orbitals with even indices (resp. odd indices) are considered as spin up (resp. spin down)
// densityDensityupup = matrix that gives the amplitude of each density-density interaction term between particles with spin up
// densityDensitydowndown = matrix that gives the amplitude of each density-density interaction term between particles with spin down
// densityDensityupdown = matrix that gives the amplitude of each density-density interaction term between particles with spin up and down
// sxSx = matrix that gives the amplitude of each Sx_i Sx_j term
// sySy = matrix that gives the amplitude of each Sy_i Sy_j term
// szSz = matrix that gives the amplitude of each Sz_i Sz_j term
// sxSy = matrix that gives the amplitude of each Sx_i Sy_j term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian::ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSites, 
																     int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, 
																     HermitianMatrix& tightBinding, 
																     RealSymmetricMatrix& densityDensityupup, 
																     RealSymmetricMatrix& densityDensitydowndown, 
																     RealSymmetricMatrix& densityDensityupdown, 
																     RealSymmetricMatrix& sxSx,
																     RealSymmetricMatrix& sySy, RealSymmetricMatrix& szSz, RealAntisymmetricMatrix& sxSy,
																     AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSites = nbrSites;
  this->XMomentum = xMomentum;
  this->MaxXMomentum = maxXMomentum;
  this->YMomentum = yMomentum;
  this->MaxYMomentum = maxYMomentum;
    
  this->LzMax = this->NbrSites - 1;
  this->HamiltonianShift = 0.0;
  this->DiagonalElements = 0;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = true;//true;
  
  this->EvaluateExponentialFactors();
  this->EvaluateOneBodyFactorsFromTightBingding(tightBinding);
  this->EvaluateInteractionFactorsFromDensityDensityAndHeisenberg(densityDensityupup, densityDensitydowndown, densityDensityupdown, sxSx, sySy, szSz, sxSy);
  
//   cout << this->InteractionFactorsupup[0][0] << endl;
    
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      cout << TmpMemory << endl;
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

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSites = number of sites
// xMomentum = momentum sector in the x direction
// maxXMomentum = number of momentum sectors in the x direction
// yMomentum = momentum sector in the x direction
// maxYMomentum = number of momentum sectors in the x direction
// tightBinding = hamiltonian corresponding to the tight-binding model in real space, orbitals with even indices (resp. odd indices) are considered as spin up (resp. spin down)
// densityDensityupup = matrix that gives the amplitude of each density-density interaction term between particles with spin up
// densityDensitydowndown = matrix that gives the amplitude of each density-density interaction term between particles with spin down
// densityDensityupdown = matrix that gives the amplitude of each density-density interaction term between particles with spin up and down
// sxSx = matrix that gives the amplitude of each Sx_i Sx_j term
// sySy = matrix that gives the amplitude of each Sy_i Sy_j term
// szSz = matrix that gives the amplitude of each Sz_i Sz_j term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian::ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSites, 
																     int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, 
																     HermitianMatrix& tightBinding, 
																     RealSymmetricMatrix& densityDensityupup, 
																     RealSymmetricMatrix& densityDensitydowndown, 
																     RealSymmetricMatrix& densityDensityupdown, 
																     RealSymmetricMatrix& sxSx,
																     RealSymmetricMatrix& sySy, RealSymmetricMatrix& szSz, RealSymmetricMatrix& sxSy, RealSymmetricMatrix& sySz, RealSymmetricMatrix& sxSz, 
																     AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSites = nbrSites;
  this->XMomentum = xMomentum;
  this->MaxXMomentum = maxXMomentum;
  this->YMomentum = yMomentum;
  this->MaxYMomentum = maxYMomentum;
    
  this->LzMax = this->NbrSites - 1;
  this->HamiltonianShift = 0.0;
  this->DiagonalElements = 0;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = true;//true;
  
  this->EvaluateExponentialFactors();
  this->EvaluateOneBodyFactorsFromTightBingding(tightBinding);
  this->EvaluateInteractionFactorsFromDensityDensityAndHeisenberg(densityDensityupup, densityDensitydowndown, densityDensityupdown, sxSx, sySy, szSz, sxSy, sySz, sxSz);
  
//   cout << this->InteractionFactorsupup[0][0] << endl;
    
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      cout << TmpMemory << endl;
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

ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian::~ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian()
{
  for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      delete[] this->ExponentialFactors[i];
    }
  delete[] this->ExponentialFactors;
}
  

// evaluate all exponential factors
//   

void ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex*[this->MaxXMomentum];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      this->ExponentialFactors[i] = new Complex[this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{ 
	  this->ExponentialFactors[i][j] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->MaxXMomentum))
							       + (this->YMomentum * ((double) j) / ((double) this->MaxYMomentum))));
	}
    }
}

// add an additional S^2 term to the Hamiltonian
//
// factor = factor in front of the S^2
// fixedSz = flag indicating whether Sz needs to be evaluated
// memory = amount of memory that can be used for S^2  precalculations

void ParticleOnLatticeWithSpinFullRealSpaceAnd2DTranslationHamiltonian::AddS2 (double factor, bool fixedSz, long memory)
{
  this->S2Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationS2Hamiltonian (this->Particles, this->NbrParticles, this->NbrSites, this->XMomentum, this->MaxXMomentum, this->YMomentum, this->MaxYMomentum, factor, fixedSz, this->Architecture, memory); 
}

