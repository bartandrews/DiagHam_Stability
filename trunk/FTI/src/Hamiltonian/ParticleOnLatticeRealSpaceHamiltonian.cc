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
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"
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

ParticleOnLatticeRealSpaceHamiltonian::ParticleOnLatticeRealSpaceHamiltonian()
{
  this->DiagonalElements = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites
// tightBinding = hamiltonian corresponding to the tight-binding model in real space
// densityDensity = matrix that gives the amplitude of each density-density interaction term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeRealSpaceHamiltonian::ParticleOnLatticeRealSpaceHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSites, 
									     HermitianMatrix& tightBinding, RealSymmetricMatrix& densityDensity,
									     AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
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

ParticleOnLatticeRealSpaceHamiltonian::~ParticleOnLatticeRealSpaceHamiltonian()
{
  if (this->DiagonalElements != 0)
    delete[] this->DiagonalElements;
}
  


// evaluate the one body interaction factors from a tight-binding matrix
//
// tightBinding = hamiltonian corresponding to the tight-binding model in real space

void ParticleOnLatticeRealSpaceHamiltonian::EvaluateOneBodyFactorsFromTightBingding (HermitianMatrix& tightBinding)
{
  if (tightBinding.GetNbrRow() != this->NbrSites)
    {
      cout << "error, the dimension of the tight binding matrix does not match the number of sites" << endl;
      return;
    }
  this->OneBodyGenericNbrConnectedSites = new int [this->NbrSites];
  this->OneBodyGenericConnectedSites = new int* [this->NbrSites];
  this->OneBodyGenericInteractionFactors = new Complex* [this->NbrSites];
  Complex Tmp;
  double Sign = 1.0;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;

  for (int i = 0; i < this->NbrSites; ++i)
    {
      this->OneBodyGenericNbrConnectedSites[i] = 0;
      for (int j = 0; j <  this->NbrSites; ++j)
	{
	  tightBinding.GetMatrixElement(i, j, Tmp);
	  if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
	    ++this->OneBodyGenericNbrConnectedSites[i];
	}
      if (this->OneBodyGenericNbrConnectedSites[i] > 0)
	{
	  this->OneBodyGenericConnectedSites[i] = new int [this->OneBodyGenericNbrConnectedSites[i]];
	  this->OneBodyGenericInteractionFactors[i] = new Complex [this->OneBodyGenericNbrConnectedSites[i]];
	  this->OneBodyGenericNbrConnectedSites[i] = 0;
 	  for (int j = 0; j <  this->NbrSites; ++j)
	    {
	      tightBinding.GetMatrixElement(i, j, Tmp);
	      if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
		{
		  this->OneBodyGenericConnectedSites[i][this->OneBodyGenericNbrConnectedSites[i]] = j;
		  this->OneBodyGenericInteractionFactors[i][this->OneBodyGenericNbrConnectedSites[i]] = Sign*Tmp;
		  ++this->OneBodyGenericNbrConnectedSites[i];
		}
	    }
	}
    }
}

// evaluate the two body interaction factors from a generic density-density interaction
//
// densityDensity = matrix that gives the amplitude of each density-density interaction term

void ParticleOnLatticeRealSpaceHamiltonian::EvaluateInteractionFactorsFromDensityDensity (RealSymmetricMatrix& densityDensity)
{
  this->NbrSectorSums = 0;  
  for (int i = 0; i < densityDensity.GetNbrRow(); ++i)
    {
      for (int j = i; j < densityDensity.GetNbrRow(); ++j)
	{
	  double Tmp;
	  densityDensity.GetMatrixElement(i, j, Tmp);
	  if (Tmp != 0.0)
	    {
	      ++this->NbrSectorSums;
	    }
	}
    }

  if (this->NbrSectorSums == 0)
    {
      cout << " I am returning as all interactions are zero! "<<endl;
      return;
    }
  this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
  this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
  this->InteractionFactors = new Complex* [this->NbrSectorSums];
  for (int i = 0; i < this->NbrSectorSums; ++i)
    {
      this->NbrSectorIndicesPerSum[i] = 1;
      this->SectorIndicesPerSum[i] = new int [2];
      this->InteractionFactors[i] = new Complex [1];
    }
  this->NbrSectorSums = 0;

  double Sign = 1.0;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;
  for (int i = 0; i < densityDensity.GetNbrRow(); ++i)
    {
      for (int j = i; j < densityDensity.GetNbrRow(); ++j)
	{
	  double Tmp;
	  densityDensity.GetMatrixElement(i, j, Tmp);
	  if (Tmp != 0.0)
	    {
	      this->SectorIndicesPerSum[this->NbrSectorSums][0] = i;
	      this->SectorIndicesPerSum[this->NbrSectorSums][1] = j;
	      this->InteractionFactors[this->NbrSectorSums][0] =  Sign * Tmp;
	      ++this->NbrSectorSums;
	    }
	}
    } 
}
  
