////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//        class of generic hamiltonian for interacting spinless particles     //
//              on lattice written in real space and p-wave pairing           //
//                                                                            //
//                        last modification : 11/06/2016                      //
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
#include "Hamiltonian/ParticleOnLatticeRealSpacePairingHamiltonian.h"
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

ParticleOnLatticeRealSpacePairingHamiltonian::ParticleOnLatticeRealSpacePairingHamiltonian()
{
}

// constructor
//
// nbrSiteX = number of sites
// tightBinding = hamiltonian corresponding to the tight-binding model in real space
// densityDensity = matrix that gives the amplitude of each density-density interaction term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeRealSpacePairingHamiltonian::ParticleOnLatticeRealSpacePairingHamiltonian(ParticleOnSphere* particles, int nbrSites, 
											   HermitianMatrix& tightBinding, RealSymmetricMatrix& densityDensity,
											   AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = 0;
  this->NbrSites = nbrSites;    
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
  this->HermitianSymmetryFlag = false;
  
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


// destructor
//

ParticleOnLatticeRealSpacePairingHamiltonian::~ParticleOnLatticeRealSpacePairingHamiltonian()
{
}
  


// evaluate the one body interaction factors from a tight-binding matrix
//
// tightBinding = hamiltonian corresponding to the tight-binding model in real space

void ParticleOnLatticeRealSpacePairingHamiltonian::EvaluateOneBodyFactorsFromTightBingding (HermitianMatrix& tightBinding)
{
  if (tightBinding.GetNbrRow() != (2 * this->NbrSites))
    {
      cout << "error, the dimension of the tight binding matrix does not match the number of sites" << endl;
      return;
    }
  this->OneBodyGenericNbrConnectedSites = new int [this->NbrSites];
  this->OneBodyGenericConnectedSites = new int* [this->NbrSites];
  this->OneBodyGenericInteractionFactors = new Complex* [this->NbrSites];
  this->OneBodyGenericPairingNbrConnectedSites = new int [this->NbrSites];
  this->OneBodyGenericPairingConnectedSites = new int* [this->NbrSites];
  this->OneBodyGenericPairingInteractionFactors = new Complex* [this->NbrSites];
  Complex Tmp;
  double Sign = 1.0;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;

  for (int i = 0; i < this->NbrSites; ++i)
    {
      this->OneBodyGenericNbrConnectedSites[i] = 0;
      for (int j = 0; j <  this->NbrSites; ++j)
	{
	  tightBinding.GetMatrixElement(2 * i, 2 * j, Tmp);
	  if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
	    ++this->OneBodyGenericNbrConnectedSites[i];
	}
      this->OneBodyGenericPairingNbrConnectedSites[i] = 0;
      for (int j = 0; j <  this->NbrSites; ++j)
	{
	  tightBinding.GetMatrixElement(2 * i, (2 * j) + 1, Tmp);
	  if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
	    ++this->OneBodyGenericPairingNbrConnectedSites[i];
	}
      if (this->OneBodyGenericNbrConnectedSites[i] > 0)
	{
	  this->OneBodyGenericConnectedSites[i] = new int [this->OneBodyGenericNbrConnectedSites[i]];
	  this->OneBodyGenericInteractionFactors[i] = new Complex [this->OneBodyGenericNbrConnectedSites[i]];
	  this->OneBodyGenericNbrConnectedSites[i] = 0;
 	  for (int j = 0; j <  this->NbrSites; ++j)
	    {
	      tightBinding.GetMatrixElement(2 * i, 2 * j, Tmp);
	      if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
		{
		  this->OneBodyGenericConnectedSites[i][this->OneBodyGenericNbrConnectedSites[i]] = j;
		  this->OneBodyGenericInteractionFactors[i][this->OneBodyGenericNbrConnectedSites[i]] = Sign*Tmp;

		  ++this->OneBodyGenericNbrConnectedSites[i];
		}
	    }
	}
      if (this->OneBodyGenericPairingNbrConnectedSites[i] > 0)
	{
	  this->OneBodyGenericPairingConnectedSites[i] = new int [this->OneBodyGenericPairingNbrConnectedSites[i]];
	  this->OneBodyGenericPairingInteractionFactors[i] = new Complex [this->OneBodyGenericPairingNbrConnectedSites[i]];
	  this->OneBodyGenericPairingNbrConnectedSites[i] = 0;
 	  for (int j = 0; j <  this->NbrSites; ++j)
	    {
	      tightBinding.GetMatrixElement(2 * i, (2 * j) + 1, Tmp);
	      if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
		{
		  this->OneBodyGenericPairingConnectedSites[i][this->OneBodyGenericPairingNbrConnectedSites[i]] = j;
		  this->OneBodyGenericPairingInteractionFactors[i][this->OneBodyGenericPairingNbrConnectedSites[i]] = Sign * Tmp;

		  ++this->OneBodyGenericPairingNbrConnectedSites[i];
		}
	    }
	}
    }
}

