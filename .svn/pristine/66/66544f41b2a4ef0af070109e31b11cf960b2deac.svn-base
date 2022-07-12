////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                     class author: Cecile Repellin                          //
//                                                                            //
//        class of generic hamiltonian for interacting spinuful particles     //
//                       on lattice written in real space                     //
//                                                                            //
//                        last modification : 15/10/2014                      //
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
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian.h"
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

ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian::ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian()
{
  this->S2Factor = 0.0;
  this->FixedSzFlag = false;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSites = number of sites
// s2Factor = numerical factor in front of S^2
// fixedSzFlag = true if the Hilbert space gas a fixed a total Sz value
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian::ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSites, 
												 double s2Factor, bool fixedSzFlag,
												 AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSites = nbrSites;
    
  this->S2Factor = s2Factor;
  this->FixedSzFlag = fixedSzFlag;

  this->LzMax = this->NbrSites - 1;
  this->HamiltonianShift = 0.0;
  this->DiagonalElements = 0;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->OneBodyGenericNbrConnectedSitesupup = 0;
  this->OneBodyGenericConnectedSitesupup = 0;
  this->OneBodyGenericInteractionFactorsupup = 0;
  this->OneBodyGenericNbrConnectedSitesdowndown = 0;
  this->OneBodyGenericConnectedSitesdowndown = 0;
  this->OneBodyGenericInteractionFactorsdowndown = 0;
  
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  
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

ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian::~ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian()
{
}
  

// evaluate the two body interaction factors from a generic density-density interaction
//

void ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian::EvaluateInteractionFactors ()
{
  this->HamiltonianShift = 0.5 * this->S2Factor * ((double) this->Particles->GetNbrParticles()); 
  if (this->FixedSzFlag == true)
    {
      this->NbrIntraSectorSums = 0;    
      this->HamiltonianShift += 0.25 * this->S2Factor * ((double) (this->Particles->GetTotalSpin() * this->Particles->GetTotalSpin()));
      this->NbrInterSectorSums = (this->NbrSites * (this->NbrSites + 1)) / 2;
      this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
      this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
      this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
      this->NbrInterSectorSums = 0;  
      for (int i = 0; i < this->NbrSites; ++i)
	{
	  this->NbrInterSectorIndicesPerSum[this->NbrInterSectorSums] = 1;
	  this->InterSectorIndicesPerSum[this->NbrInterSectorSums] = new int [2];
	  this->InteractionFactorsupdown[this->NbrInterSectorSums] = new Complex [1];
	  ++this->NbrInterSectorSums;
	  for (int j =  i + 1; j < this->NbrSites; ++j)
	    {
	      this->NbrInterSectorIndicesPerSum[this->NbrInterSectorSums] = 2;
	      this->InterSectorIndicesPerSum[this->NbrInterSectorSums] = new int [4];
	      this->InteractionFactorsupdown[this->NbrInterSectorSums] = new Complex [4];
	      ++this->NbrInterSectorSums;
	    }
	}
      this->NbrInterSectorSums = 0;
      for (int i = 0; i < this->NbrSites; ++i)
	{
	  this->InterSectorIndicesPerSum[this->NbrInterSectorSums][0] = i;
	  this->InterSectorIndicesPerSum[this->NbrInterSectorSums][1] = i;
	  this->InteractionFactorsupdown[this->NbrInterSectorSums][0] =  this->S2Factor;
	  ++this->NbrInterSectorSums;
	  for (int j = i + 1; j < this->NbrSites; ++j)
	    {
	      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][0] = i;
	      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][1] = j;
	      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][2] = j;
	      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][3] = i;
	      this->InteractionFactorsupdown[this->NbrInterSectorSums][0] = 0.0;
	      this->InteractionFactorsupdown[this->NbrInterSectorSums][1] = this->S2Factor;
	      this->InteractionFactorsupdown[this->NbrInterSectorSums][2] = this->S2Factor;
	      this->InteractionFactorsupdown[this->NbrInterSectorSums][3] = 0.0;
	      ++this->NbrInterSectorSums;
	    }
	} 
    }
  else
    {
      cout << "error, ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian requires a fixes Sz Hilbert space" << endl;
    }
}


