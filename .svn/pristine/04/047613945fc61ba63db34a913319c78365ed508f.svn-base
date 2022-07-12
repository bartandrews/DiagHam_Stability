////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//        class of generic hamiltonian for interacting spinless particles     //
//         on lattice written in real space and handling 2d translations      //
//                                                                            //
//                        last modification : 02/10/2014                      //
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
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian.h"
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

ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian::ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian()
{
  this->PhaseTranslationY = 0;
  this->PhaseTranslationX = 0;
  this->MagneticTranslationPhaseFactor = 0; 
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites
// xMomentum = momentum sector in the x direction
// maxXMomentum = number of momentum sectors in the x direction
// yMomentum = momentum sector in the x direction
// maxYMomentum = number of momentum sectors in the x direction
// tightBinding = hamiltonian corresponding to the tight-binding model in real space
// densityDensity = matrix that gives the amplitude of each density-density interaction term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian::ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSites, 
							int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum,
							 double phaseTranslationX, double phaseTranslationY,
							HermitianMatrix& tightBinding, RealSymmetricMatrix& densityDensity,
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
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = false;//true;
  this->PhaseTranslationX = phaseTranslationX;
  this->PhaseTranslationY = phaseTranslationY;  
  this->EvaluateExponentialFactors();
//  this->ComputeMagneticTranslationPhaseFactor();
  this->EvaluateOneBodyFactorsFromTightBingding(tightBinding);
  this->EvaluateInteractionFactorsFromDensityDensity(densityDensity);
    
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


  
ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian::~ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian()
{
}
