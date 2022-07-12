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
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceHamiltonian.h"
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

ParticleOnLatticeWithSpinRealSpaceHamiltonian::ParticleOnLatticeWithSpinRealSpaceHamiltonian()
{
  this->DiagonalElements = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSites = number of sites
  // tightBindingupup = hamiltonian corresponding to the tight-binding model in real space for particles with spin up
  // tightBindingdowndown = hamiltonian corresponding to the tight-binding model in real space for particles with spin down
  // densityDensityupup = matrix that gives the amplitude of each density-density interaction term between particles with spin up
  // densityDensitydowndown = matrix that gives the amplitude of each density-density interaction term between particles with spin down
  // densityDensityupdown = matrix that gives the amplitude of each density-density interaction term between particles with spin up and down
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinRealSpaceHamiltonian::ParticleOnLatticeWithSpinRealSpaceHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSites, 
													     HermitianMatrix& tightBindingupup,HermitianMatrix& tightBindingdowndown, RealSymmetricMatrix& densityDensityupup, RealSymmetricMatrix& densityDensitydowndown, RealSymmetricMatrix& densityDensityupdown,
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
  this->HermitianSymmetryFlag = true;
  
  this->EvaluateOneBodyFactorsFromTightBingding(tightBindingupup, tightBindingdowndown);
  
  this->EvaluateInteractionFactorsFromDensityDensity(densityDensityupup, densityDensitydowndown, densityDensityupdown);
    
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

ParticleOnLatticeWithSpinRealSpaceHamiltonian::~ParticleOnLatticeWithSpinRealSpaceHamiltonian()
{
  if (this->DiagonalElements != 0)
    delete[] this->DiagonalElements;
}
  


// evaluate the one body interaction factors from a tight-binding matrix
//
// tightBindingupup = hamiltonian corresponding to the tight-binding model in real space for particles with spin up
// tightBindingdowndown = hamiltonian corresponding to the tight-binding model in real space for particles with spin down

void ParticleOnLatticeWithSpinRealSpaceHamiltonian::EvaluateOneBodyFactorsFromTightBingding (HermitianMatrix& tightBindingupup, HermitianMatrix& tightBindingdowndown)
{
  if ((tightBindingupup.GetNbrRow() != this->NbrSites) || (tightBindingdowndown.GetNbrRow() != this->NbrSites))
    {
      cout << "error, the dimension of the tight binding matrix does not match the number of sites" << endl;
      return;
    }
  this->OneBodyGenericNbrConnectedSitesupup = new int [this->NbrSites];
  this->OneBodyGenericConnectedSitesupup = new int* [this->NbrSites];
  this->OneBodyGenericInteractionFactorsupup = new Complex* [this->NbrSites];
  
  this->OneBodyGenericNbrConnectedSitesdowndown = new int [this->NbrSites];
  this->OneBodyGenericConnectedSitesdowndown = new int* [this->NbrSites];
  this->OneBodyGenericInteractionFactorsdowndown = new Complex* [this->NbrSites];
  
  Complex Tmp;
  double Sign = 1.0;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;
  for (int i = 0; i < this->NbrSites; ++i)
    {
      this->OneBodyGenericNbrConnectedSitesupup[i] = 0;
      this->OneBodyGenericNbrConnectedSitesdowndown[i] = 0;
      for (int j = 0; j <  this->NbrSites; ++j)
	{
	  tightBindingupup.GetMatrixElement(i, j, Tmp);
	  if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
	    ++this->OneBodyGenericNbrConnectedSitesupup[i];
	  tightBindingdowndown.GetMatrixElement(i, j, Tmp);
	  if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
	    ++this->OneBodyGenericNbrConnectedSitesdowndown[i];
	}
      if (this->OneBodyGenericNbrConnectedSitesupup[i] > 0)
	{
	  this->OneBodyGenericConnectedSitesupup[i] = new int [this->OneBodyGenericNbrConnectedSitesupup[i]];
	  this->OneBodyGenericInteractionFactorsupup[i] = new Complex [this->OneBodyGenericNbrConnectedSitesupup[i]];
	  this->OneBodyGenericNbrConnectedSitesupup[i] = 0;
 	  for (int j = 0; j <  this->NbrSites; ++j)
	    {
	      tightBindingupup.GetMatrixElement(i, j, Tmp);
	      if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
		{
		  this->OneBodyGenericConnectedSitesupup[i][this->OneBodyGenericNbrConnectedSitesupup[i]] = j;
		  this->OneBodyGenericInteractionFactorsupup[i][this->OneBodyGenericNbrConnectedSitesupup[i]] = Sign*Tmp;
		  ++this->OneBodyGenericNbrConnectedSitesupup[i];
		}
	    }
	}
      if (this->OneBodyGenericNbrConnectedSitesdowndown[i] > 0)
	{
	  this->OneBodyGenericConnectedSitesdowndown[i] = new int [this->OneBodyGenericNbrConnectedSitesdowndown[i]];
	  this->OneBodyGenericInteractionFactorsdowndown[i] = new Complex [this->OneBodyGenericNbrConnectedSitesdowndown[i]];
	  this->OneBodyGenericNbrConnectedSitesdowndown[i] = 0;
 	  for (int j = 0; j <  this->NbrSites; ++j)
	    {
	      tightBindingdowndown.GetMatrixElement(i, j, Tmp);
	      if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
		{
		  this->OneBodyGenericConnectedSitesdowndown[i][this->OneBodyGenericNbrConnectedSitesdowndown[i]] = j;
		  this->OneBodyGenericInteractionFactorsdowndown[i][this->OneBodyGenericNbrConnectedSitesdowndown[i]] = Sign*Tmp;
		  ++this->OneBodyGenericNbrConnectedSitesdowndown[i];
		}
	    }
	}
    }
}

// evaluate the two body interaction factors from a generic density-density interaction
//
// densityDensityIntra = matrix that gives the amplitude of each density-density interaction term for particles with the same spin
// densityDensityInter = matrix that gives the amplitude of each density-density interaction term for particles with opposite spins

void ParticleOnLatticeWithSpinRealSpaceHamiltonian::EvaluateInteractionFactorsFromDensityDensity (RealSymmetricMatrix& densityDensityupup, RealSymmetricMatrix& densityDensitydowndown, RealSymmetricMatrix& densityDensityupdown)
{
  this->NbrIntraSectorSums = 0;  
  for (int i = 0; i < densityDensityupup.GetNbrRow(); ++i)
    {
      for (int j = i; j < densityDensityupup.GetNbrRow(); ++j)
	{
	  double Tmp;
	  double Tmp1;
	  densityDensityupup.GetMatrixElement(i, j, Tmp);
	  densityDensitydowndown.GetMatrixElement(i, j, Tmp1);
	  if ((Tmp != 0.0) || (Tmp1 != 0.0))
	    {
	      ++this->NbrIntraSectorSums;
	    }
	}
    }
    
    
  this->NbrInterSectorSums = 0;
  for (int i = 0; i < densityDensityupdown.GetNbrRow(); ++i)
    {
      for (int j = 0; j < densityDensityupdown.GetNbrRow(); ++j)
	{
	  double Tmp;
	  densityDensityupdown.GetMatrixElement(i, j, Tmp);
	  if (Tmp != 0.0)
	    {
	      ++this->NbrInterSectorSums;
	    }
	}
    }
  if ((this->NbrInterSectorSums == 0) && (this->NbrIntraSectorSums == 0))
    {
      return;
    }
    
  double Sign = 1.0;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;
  
  if (this->NbrInterSectorSums != 0)
    {
      this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
      this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
      this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->NbrInterSectorIndicesPerSum[i] = 1;
	  this->InterSectorIndicesPerSum[i] = new int [2];
	  this->InteractionFactorsupdown[i] = new Complex [1];
	}
      this->NbrInterSectorSums = 0;
      for (int i = 0; i < densityDensityupdown.GetNbrRow(); ++i)
	{
	  for (int j = 0; j < densityDensityupdown.GetNbrRow(); ++j)
	    {
	      double Tmp;
	      densityDensityupdown.GetMatrixElement(i, j, Tmp);
	      if (Tmp != 0.0)
		{
		  this->InterSectorIndicesPerSum[this->NbrInterSectorSums][0] = i;
		  this->InterSectorIndicesPerSum[this->NbrInterSectorSums][1] = j;
		  this->InteractionFactorsupdown[this->NbrInterSectorSums][0] =  Sign * Tmp;
		  ++this->NbrInterSectorSums;
		}
	    }
	}
    }
  
  if (this->NbrIntraSectorSums != 0)
    {
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
      
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->NbrIntraSectorIndicesPerSum[i] = 1;
	  this->IntraSectorIndicesPerSum[i] = new int [2];
	  this->InteractionFactorsupup[i] = new Complex [1];
	  this->InteractionFactorsdowndown[i] = new Complex [1];
	}
      this->NbrIntraSectorSums = 0;
      for (int i = 0; i < densityDensityupup.GetNbrRow(); ++i)
	{
	  for (int j = i; j < densityDensityupup.GetNbrRow(); ++j)
	    {
	      double Tmp;
	      double Tmp1;
	      densityDensityupup.GetMatrixElement(i, j, Tmp);
	      densityDensitydowndown.GetMatrixElement(i, j, Tmp1);
	      if ((Tmp != 0.0) || (Tmp1 != 0.0))
		{
		  this->IntraSectorIndicesPerSum[this->NbrIntraSectorSums][0] = i;
		  this->IntraSectorIndicesPerSum[this->NbrIntraSectorSums][1] = j;
		  this->InteractionFactorsupup[this->NbrIntraSectorSums][0] =  Sign * Tmp;  
		  this->InteractionFactorsdowndown[this->NbrIntraSectorSums][0] =  Sign * Tmp1;  
		  
		  ++this->NbrIntraSectorSums;
		}
	    }
	}
    }
}
  
// add an additional S^2 term to the Hamiltonian
//
// factor = factor in front of the S^2
// fixedSz = flag indicating whether Sz needs to be evaluated
// memory = amount of memory that can be used for S^2  precalculations

void ParticleOnLatticeWithSpinRealSpaceHamiltonian::AddS2 (double factor, bool fixedSz, long memory)
{
  this->S2Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian (this->Particles, this->NbrParticles, this->NbrSites, factor, fixedSz, this->Architecture, memory); 
}
