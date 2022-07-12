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
//                        last modification : 11/09/2014                      //
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
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian.h"
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

ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian::ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian()
{
  this->XMomentum = 0;
  this->MaxXMomentum = 1;
  this->DiagonalElements = 0;
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

ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian::ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSites, 
													     int xMomentum, int maxXMomentum, HermitianMatrix& tightBinding, RealSymmetricMatrix& densityDensity,
													     AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSites = nbrSites;
  this->XMomentum = xMomentum;
  this->MaxXMomentum = maxXMomentum; 
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
//  cout <<" MaxIndex = "<< MaxIndex<<" MinIndex  = "<< MinIndex <<endl;
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = true;
  
  this->EvaluateExponentialFactors();
  this->EvaluateOneBodyFactorsFromTightBingding(tightBinding);
  this->EvaluateInteractionFactorsFromDensityDensity(densityDensity);
    
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

ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian::~ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian()
{
  if (this->DiagonalElements != 0)
    delete[] this->DiagonalElements;
  delete[] this->ExponentialFactors;
}
  

// evaluate all exponential factors
//   

void ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex[this->MaxXMomentum];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      this->ExponentialFactors[i] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->MaxXMomentum))));
    }
}

