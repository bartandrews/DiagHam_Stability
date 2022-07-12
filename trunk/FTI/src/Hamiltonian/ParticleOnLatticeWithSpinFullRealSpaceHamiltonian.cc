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
//           on lattice written in real space without Sz conservation         //
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
#include "Hamiltonian/ParticleOnLatticeWithSpinFullRealSpaceHamiltonian.h"
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

ParticleOnLatticeWithSpinFullRealSpaceHamiltonian::ParticleOnLatticeWithSpinFullRealSpaceHamiltonian()
{
  this->DiagonalElements = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSites = number of sites
// tightBinding = hamiltonian corresponding to the tight-binding model in real space, orbitals with even indices (resp. odd indices) are considered as spin up (resp. spin down)
// densityDensityupup = matrix that gives the amplitude of each density-density interaction term between particles with spin up
// densityDensitydowndown = matrix that gives the amplitude of each density-density interaction term between particles with spin down
// densityDensityupdown = matrix that gives the amplitude of each density-density interaction term between particles with spin up and down
// sxSx = matrix that gives the amplitude of each Sx_i Sx_j term
// sySy = matrix that gives the amplitude of each Sy_i Sy_j term
// szSz = matrix that gives the amplitude of each Sz_i Sz_j term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinFullRealSpaceHamiltonian::ParticleOnLatticeWithSpinFullRealSpaceHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSites, 
												     HermitianMatrix& tightBinding, 
												     RealSymmetricMatrix& densityDensityupup, RealSymmetricMatrix& densityDensitydowndown, 
												     RealSymmetricMatrix& densityDensityupdown, RealSymmetricMatrix& sxSx,
												     RealSymmetricMatrix& sySy, RealSymmetricMatrix& szSz,
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
  
  this->EvaluateOneBodyFactorsFromTightBingding(tightBinding);
  
  this->EvaluateInteractionFactorsFromDensityDensityAndHeisenberg(densityDensityupup, densityDensitydowndown, densityDensityupdown, sxSx, sySy, szSz);
    
  cout << (this->NbrInterSectorSums) << " " << (this->NbrIntraSectorSums) << endl;
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


// constructor with SxSy terms
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSites = number of sites
// tightBinding = hamiltonian corresponding to the tight-binding model in real space, orbitals with even indices (resp. odd indices) are considered as spin up (resp. spin down)
// densityDensityupup = matrix that gives the amplitude of each density-density interaction term between particles with spin up
// densityDensitydowndown = matrix that gives the amplitude of each density-density interaction term between particles with spin down
// densityDensityupdown = matrix that gives the amplitude of each density-density interaction term between particles with spin up and down
// sxSx = matrix that gives the amplitude of each Sx_i Sx_j term
// sySy = matrix that gives the amplitude of each Sy_i Sy_j term
// szSz = matrix that gives the amplitude of each Sz_i Sz_j term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinFullRealSpaceHamiltonian::ParticleOnLatticeWithSpinFullRealSpaceHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSites, 
												     HermitianMatrix& tightBinding, 
												     RealSymmetricMatrix& densityDensityupup, RealSymmetricMatrix& densityDensitydowndown, 
												     RealSymmetricMatrix& densityDensityupdown, RealSymmetricMatrix& sxSx,
												     RealSymmetricMatrix& sySy, RealSymmetricMatrix& szSz, RealAntisymmetricMatrix& sxSy,
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
  
  this->EvaluateOneBodyFactorsFromTightBingding(tightBinding);
  
  this->EvaluateInteractionFactorsFromDensityDensityAndHeisenberg(densityDensityupup, densityDensitydowndown, densityDensityupdown, sxSx, sySy, szSz, sxSy);
    
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

// constructor with symmetric SxSy, SySz and SxSz terms
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSites = number of sites
// tightBinding = hamiltonian corresponding to the tight-binding model in real space, orbitals with even indices (resp. odd indices) are considered as spin up (resp. spin down)
// densityDensityupup = matrix that gives the amplitude of each density-density interaction term between particles with spin up
// densityDensitydowndown = matrix that gives the amplitude of each density-density interaction term between particles with spin down
// densityDensityupdown = matrix that gives the amplitude of each density-density interaction term between particles with spin up and down
// sxSx = matrix that gives the amplitude of each Sx_i Sx_j term
// sySy = matrix that gives the amplitude of each Sy_i Sy_j term
// szSz = matrix that gives the amplitude of each Sz_i Sz_j term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
ParticleOnLatticeWithSpinFullRealSpaceHamiltonian::ParticleOnLatticeWithSpinFullRealSpaceHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSites, 
						    HermitianMatrix& tightBinding, RealSymmetricMatrix& densityDensityupup, RealSymmetricMatrix& densityDensitydowndown, 
						    RealSymmetricMatrix& densityDensityupdown, RealSymmetricMatrix& sxSx,
						    RealSymmetricMatrix& sySy, RealSymmetricMatrix& szSz, RealSymmetricMatrix& sxSy, RealSymmetricMatrix& sySz, RealSymmetricMatrix& sxSz,
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
  
  this->EvaluateOneBodyFactorsFromTightBingding(tightBinding);
  
  this->EvaluateInteractionFactorsFromDensityDensityAndHeisenberg(densityDensityupup, densityDensitydowndown, densityDensityupdown, sxSx, sySy, szSz, sxSy, sySz, sxSz);
    
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

ParticleOnLatticeWithSpinFullRealSpaceHamiltonian::~ParticleOnLatticeWithSpinFullRealSpaceHamiltonian()
{
  if (this->DiagonalElements != 0)
    delete[] this->DiagonalElements;
}
  


// evaluate the one body interaction factors from a tight-binding matrix
//
// tightBinding = hamiltonian corresponding to the tight-binding model in real space, orbitals with even indices (resp. odd indices) are considered as spin up (resp. spin down)

void ParticleOnLatticeWithSpinFullRealSpaceHamiltonian::EvaluateOneBodyFactorsFromTightBingding (HermitianMatrix& tightBinding)
{
  if (tightBinding.GetNbrRow() != (2 *this->NbrSites))
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
  
  this->OneBodyGenericNbrConnectedSitesupdown = new int [this->NbrSites];
  this->OneBodyGenericConnectedSitesupdown = new int* [this->NbrSites];
  this->OneBodyGenericInteractionFactorsupdown = new Complex* [this->NbrSites];
  
  this->OneBodyGenericNbrConnectedSitesdownup = new int [this->NbrSites];
  this->OneBodyGenericConnectedSitesdownup = new int* [this->NbrSites];
  this->OneBodyGenericInteractionFactorsdownup = new Complex* [this->NbrSites];
  
  Complex Tmp;
  double Sign = 1.0;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;
  for (int i = 0; i < this->NbrSites; ++i)
    {
      this->OneBodyGenericNbrConnectedSitesupup[i] = 0;
      this->OneBodyGenericNbrConnectedSitesdowndown[i] = 0;
      this->OneBodyGenericNbrConnectedSitesupdown[i] = 0;
      this->OneBodyGenericNbrConnectedSitesdownup[i] = 0;
      for (int j = 0; j <  this->NbrSites; ++j)
	{
	  tightBinding.GetMatrixElement(2 * i, 2 * j, Tmp);
	  if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
	    ++this->OneBodyGenericNbrConnectedSitesupup[i];
	  tightBinding.GetMatrixElement((2 * i) + 1, (2 * j) + 1, Tmp);
	  if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
	    ++this->OneBodyGenericNbrConnectedSitesdowndown[i];
	  tightBinding.GetMatrixElement(2 * i, (2 * j) + 1, Tmp);
	  if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
	    ++this->OneBodyGenericNbrConnectedSitesupdown[i];
	  tightBinding.GetMatrixElement((2 * i) + 1, 2 * j, Tmp);
	  if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
	    ++this->OneBodyGenericNbrConnectedSitesdownup[i];
	}
      if (this->OneBodyGenericNbrConnectedSitesupup[i] > 0)
	{
	  this->OneBodyGenericConnectedSitesupup[i] = new int [this->OneBodyGenericNbrConnectedSitesupup[i]];
	  this->OneBodyGenericInteractionFactorsupup[i] = new Complex [this->OneBodyGenericNbrConnectedSitesupup[i]];
	  this->OneBodyGenericNbrConnectedSitesupup[i] = 0;
 	  for (int j = 0; j <  this->NbrSites; ++j)
	    {
	      tightBinding.GetMatrixElement(2 * i, 2 * j, Tmp);
	      if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
		{
		  this->OneBodyGenericConnectedSitesupup[i][this->OneBodyGenericNbrConnectedSitesupup[i]] = j;
		  this->OneBodyGenericInteractionFactorsupup[i][this->OneBodyGenericNbrConnectedSitesupup[i]] = Sign * Tmp;
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
	      tightBinding.GetMatrixElement((2 * i) + 1, (2 * j) + 1, Tmp);
	      if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
		{
		  this->OneBodyGenericConnectedSitesdowndown[i][this->OneBodyGenericNbrConnectedSitesdowndown[i]] = j;
		  this->OneBodyGenericInteractionFactorsdowndown[i][this->OneBodyGenericNbrConnectedSitesdowndown[i]] = Sign * Tmp;
		  ++this->OneBodyGenericNbrConnectedSitesdowndown[i];
		}
	    }
	}
      if (this->OneBodyGenericNbrConnectedSitesupdown[i] > 0)
	{
	  this->OneBodyGenericConnectedSitesupdown[i] = new int [this->OneBodyGenericNbrConnectedSitesupdown[i]];
	  this->OneBodyGenericInteractionFactorsupdown[i] = new Complex [this->OneBodyGenericNbrConnectedSitesupdown[i]];
	  this->OneBodyGenericNbrConnectedSitesupdown[i] = 0;
 	  for (int j = 0; j <  this->NbrSites; ++j)
	    {
	      tightBinding.GetMatrixElement(2 * i, (2 * j) + 1, Tmp);
	      if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
		{
		  this->OneBodyGenericConnectedSitesupdown[i][this->OneBodyGenericNbrConnectedSitesupdown[i]] = j;
		  this->OneBodyGenericInteractionFactorsupdown[i][this->OneBodyGenericNbrConnectedSitesupdown[i]] = Sign * Tmp;
		  ++this->OneBodyGenericNbrConnectedSitesupdown[i];
		}
	    }
	}
      if (this->OneBodyGenericNbrConnectedSitesdownup[i] > 0)
	{
	  this->OneBodyGenericConnectedSitesdownup[i] = new int [this->OneBodyGenericNbrConnectedSitesdownup[i]];
	  this->OneBodyGenericInteractionFactorsdownup[i] = new Complex [this->OneBodyGenericNbrConnectedSitesdownup[i]];
	  this->OneBodyGenericNbrConnectedSitesdownup[i] = 0;
 	  for (int j = 0; j <  this->NbrSites; ++j)
	    {
	      tightBinding.GetMatrixElement((2 * i) + 1, 2 * j, Tmp);
	      if ((Tmp.Re != 0.0) || (Tmp.Im != 0.0))
		{
		  this->OneBodyGenericConnectedSitesdownup[i][this->OneBodyGenericNbrConnectedSitesdownup[i]] = j;
		  this->OneBodyGenericInteractionFactorsdownup[i][this->OneBodyGenericNbrConnectedSitesdownup[i]] = Sign * Tmp;
		  ++this->OneBodyGenericNbrConnectedSitesdownup[i];
		}
	    }
	}
    }
}

// evaluate the two body interaction factors from a generic density-density interaction and a generic anisotropic Heisenberg interaction
//
// densityDensityIntra = matrix that gives the amplitude of each density-density interaction term for particles with the same spin
// densityDensityInter = matrix that gives the amplitude of each density-density interaction term for particles with opposite spins
// sxSx = matrix that gives the amplitude of each Sx_i Sx_j term
// sySy = matrix that gives the amplitude of each Sy_i Sy_j term
// szSz = matrix that gives the amplitude of each Sz_i Sz_j term

void ParticleOnLatticeWithSpinFullRealSpaceHamiltonian::EvaluateInteractionFactorsFromDensityDensityAndHeisenberg (RealSymmetricMatrix& densityDensityupup, 
														   RealSymmetricMatrix& densityDensitydowndown,
														   RealSymmetricMatrix& densityDensityupdown, RealSymmetricMatrix& sxSx,
														   RealSymmetricMatrix& sySy, RealSymmetricMatrix& szSz)
{
  long TotalNbrInteractionFactors = 0;
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;

  this->InteractionFactorsupupupup = 0;
  this->InteractionFactorsdowndowndowndown = 0;
  this->InteractionFactorsupupdowndown = 0;
  this->InteractionFactorsdowndownupup = 0;
  this->InteractionFactorsupdownupup = 0;
  this->InteractionFactorsupdowndowndown = 0;
  this->InteractionFactorsupupupdown = 0;
  this->InteractionFactorsdowndownupdown = 0;
  this->InteractionFactorsupdownupdown = 0;
  
  this->NbrIntraSectorSums = 0;  
  for (int i = 0; i < densityDensityupup.GetNbrRow(); ++i)
    {
      for (int j = i; j < densityDensityupup.GetNbrRow(); ++j)
	{
	  double Tmpupup;
	  double Tmpdowndown;
	  double TmpSxSx;
	  double TmpSySy;
	  double TmpSzSz;
	  densityDensityupup.GetMatrixElement(i, j, Tmpupup);
	  densityDensitydowndown.GetMatrixElement(i, j, Tmpdowndown);
	  sxSx.GetMatrixElement(i, j, TmpSxSx);
	  sySy.GetMatrixElement(i, j, TmpSySy);
	  szSz.GetMatrixElement(i, j, TmpSzSz);
	  if ((Tmpupup != 0.0) || (Tmpdowndown != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0))
	    {
	      ++this->NbrIntraSectorSums;
	    }
	}
    }
    
    
  this->NbrInterSectorSums = 0;
  for (int i = 0; i < densityDensityupdown.GetNbrRow(); ++i)
    {
      for (int j = i; j < densityDensityupdown.GetNbrRow(); ++j)
	{
	  double Tmpupdown;
	  double TmpSxSx;
	  double TmpSySy;
	  double TmpSzSz;
	  densityDensityupdown.GetMatrixElement(i, j, Tmpupdown);
	  sxSx.GetMatrixElement(i, j, TmpSxSx);
	  sySy.GetMatrixElement(i, j, TmpSySy);
	  szSz.GetMatrixElement(i, j, TmpSzSz);
	  if ((Tmpupdown != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0))
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
      this->InteractionFactorsupdownupdown = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsupupupdown = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsdowndownupdown = new Complex* [this->NbrInterSectorSums];
      this->NbrInterSectorSums = 0;
      for (int i = 0; i < densityDensityupdown.GetNbrRow(); ++i)
	{
	  for (int j = i; j < densityDensityupdown.GetNbrRow(); ++j)
	    {
	      double Tmpupdown;
	      double TmpSxSx;
	      double TmpSySy;
	      double TmpSzSz;
	      densityDensityupdown.GetMatrixElement(i, j, Tmpupdown);
	      sxSx.GetMatrixElement(i, j, TmpSxSx);
	      sySy.GetMatrixElement(i, j, TmpSySy);
	      szSz.GetMatrixElement(i, j, TmpSzSz);
	      if ((Tmpupdown != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0))
 		{
// 		  if ()
// 		    {
// 		      this->NbrInterSectorIndicesPerSum[this->NbrInterSectorSums] = 1;
// 		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums] = new int [2];
// 		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums] = new Complex [1];
// 		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums] = new Complex [1];
// 		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums] = new Complex [1];
// 		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][0] = i;
// 		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][1] = j;
// 		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][0] =  Sign * (Tmpupdown - (0.5 * TmpSzSz));
// 		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums][0] = 0.0;
// 		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums][0] = 0.0;
// 		    }
// 		  else
// 		    {
		      this->NbrInterSectorIndicesPerSum[this->NbrInterSectorSums] = 2;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums] = new int [4];
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums] = new Complex [4];
		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums] = new Complex [2];
		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums] = new Complex [2];
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][0] = i;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][1] = j;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][2] = j;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][3] = i;
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][0] =  Sign * (Tmpupdown - (0.25 * TmpSzSz));
		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums][0] = 0.0;
		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums][0] = 0.0;
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][1] =  -0.25 * Sign * (TmpSxSx + TmpSySy);
		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums][1] = 0.0;
		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums][1] = 0.0;
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][2] =  -0.25 * Sign * (TmpSxSx + TmpSySy);
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][3] =  Sign * (Tmpupdown -0.25 * TmpSzSz);
// 		    }
		  ++this->NbrInterSectorSums;
		}
	    }
	}
    }

  
  if (this->NbrIntraSectorSums != 0)
    {
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];
      this->NbrIntraSectorSums = 0;
      for (int i = 0; i < densityDensityupup.GetNbrRow(); ++i)
	{
	  for (int j = i; j < densityDensityupup.GetNbrRow(); ++j)
	    {
	      double Tmp;
	      double Tmp1;
	      double TmpSxSx;
	      double TmpSySy;
	      double TmpSzSz;
	      densityDensityupup.GetMatrixElement(i, j, Tmp);
	      densityDensitydowndown.GetMatrixElement(i, j, Tmp1);
	      sxSx.GetMatrixElement(i, j, TmpSxSx);
	      sySy.GetMatrixElement(i, j, TmpSySy);
	      szSz.GetMatrixElement(i, j, TmpSzSz);
	      if ((Tmp != 0.0) || (Tmp1 != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0))
		{
		  this->NbrIntraSectorIndicesPerSum[this->NbrIntraSectorSums] = 1;
		  this->IntraSectorIndicesPerSum[this->NbrIntraSectorSums] = new int [2];
		  this->InteractionFactorsupupupup[this->NbrIntraSectorSums] = new Complex [1];
		  this->InteractionFactorsdowndowndowndown[this->NbrIntraSectorSums] = new Complex [1];
		  this->InteractionFactorsupupdowndown[this->NbrIntraSectorSums] = new Complex [1];
		  this->InteractionFactorsdowndownupup[this->NbrIntraSectorSums] = new Complex [1];
		  this->IntraSectorIndicesPerSum[this->NbrIntraSectorSums][0] = i;
		  this->IntraSectorIndicesPerSum[this->NbrIntraSectorSums][1] = j;
		  this->InteractionFactorsupupupup[this->NbrIntraSectorSums][0] =  Sign * (Tmp + 0.25 * TmpSzSz);  
		  this->InteractionFactorsdowndowndowndown[this->NbrIntraSectorSums][0] =  Sign * (Tmp1 + 0.25 * TmpSzSz);  
		  this->InteractionFactorsupupdowndown[this->NbrIntraSectorSums][0] = 0.25 * Sign * (TmpSxSx - TmpSySy);
		  this->InteractionFactorsdowndownupup[this->NbrIntraSectorSums][0] = 0.25 * Sign * (TmpSxSx - TmpSySy);
// 		  if ((TmpSxSx == 0.0) && (TmpSySy == 0.0))
// 		    {
// 		      this->InteractionFactorsupdownupup[this->NbrIntraSectorSums] = new Complex [1];
// 		      this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums] = new Complex [1];
// 		      this->InteractionFactorsupdownupup[this->NbrIntraSectorSums][0] = 0.0;
// 		      this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums][0] = 0.0;
// 		    }
// 		  else
// 		    {
		      this->InteractionFactorsupdownupup[this->NbrIntraSectorSums] = new Complex [2];
		      this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums] = new Complex [2];
		      this->InteractionFactorsupdownupup[this->NbrIntraSectorSums][0] = 0.0;
		      this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums][0] = 0.0;
		      this->InteractionFactorsupdownupup[this->NbrIntraSectorSums][1] = 0.0;
		      this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums][1] = 0.0;
		      //		    }
		  ++this->NbrIntraSectorSums;
		}
	    }
	}
    }
  else
    {
      this->InteractionFactorsupupupup = 0;
      this->InteractionFactorsdowndowndowndown = 0;
      this->InteractionFactorsupupdowndown = 0;
      this->InteractionFactorsdowndownupup = 0;
      this->InteractionFactorsupdownupup = 0;
      this->InteractionFactorsupdowndowndown = 0;
    }
}

// evaluate the two body interaction factors from a generic density-density interaction and a generic anisotropic Heisenberg interaction
//
// densityDensityIntra = matrix that gives the amplitude of each density-density interaction term for particles with the same spin
// densityDensityInter = matrix that gives the amplitude of each density-density interaction term for particles with opposite spins
// sxSx = matrix that gives the amplitude of each Sx_i Sx_j term
// sySy = matrix that gives the amplitude of each Sy_i Sy_j term
// szSz = matrix that gives the amplitude of each Sz_i Sz_j term
// sxSy = matrix that gives the amplitude of each Sx_i Sy_j term

void ParticleOnLatticeWithSpinFullRealSpaceHamiltonian::EvaluateInteractionFactorsFromDensityDensityAndHeisenberg (RealSymmetricMatrix& densityDensityupup, 
														   RealSymmetricMatrix& densityDensitydowndown,
														   RealSymmetricMatrix& densityDensityupdown, RealSymmetricMatrix& sxSx,
														   RealSymmetricMatrix& sySy, RealSymmetricMatrix& szSz, RealAntisymmetricMatrix& sxSy)
{
  
  this->InteractionFactorsupupupup = 0;
  this->InteractionFactorsdowndowndowndown = 0;
  this->InteractionFactorsupupdowndown = 0;
  this->InteractionFactorsdowndownupup = 0;
  this->InteractionFactorsupdownupup = 0;
  this->InteractionFactorsupdowndowndown = 0;
  this->InteractionFactorsupupupdown = 0;
  this->InteractionFactorsdowndownupdown = 0;
  this->InteractionFactorsupdownupdown = 0;
  
  long TotalNbrInteractionFactors = 0;
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;

  this->NbrIntraSectorSums = 0;  
  for (int i = 0; i < densityDensityupup.GetNbrRow(); ++i)
    {
      for (int j = i; j < densityDensityupup.GetNbrRow(); ++j)
	{
	  double Tmpupup;
	  double Tmpdowndown;
	  double TmpSxSx;
	  double TmpSySy;
	  double TmpSzSz;
	  double TmpSxSy;
	  densityDensityupup.GetMatrixElement(i, j, Tmpupup);
	  densityDensitydowndown.GetMatrixElement(i, j, Tmpdowndown);
	  sxSx.GetMatrixElement(i, j, TmpSxSx);
	  sySy.GetMatrixElement(i, j, TmpSySy);
	  szSz.GetMatrixElement(i, j, TmpSzSz);
	  sxSy.GetMatrixElement(i, j, TmpSxSy);
	  if ((Tmpupup != 0.0) || (Tmpdowndown != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0) || (TmpSxSy != 0.0))
	    {
	      ++this->NbrIntraSectorSums;
	    }
	}
    }
    
    
  this->NbrInterSectorSums = 0;
  for (int i = 0; i < densityDensityupdown.GetNbrRow(); ++i)
    {
      for (int j = i; j < densityDensityupdown.GetNbrRow(); ++j)
	{
	  double Tmpupdown;
	  double TmpSxSx;
	  double TmpSySy;
	  double TmpSzSz;
	  double TmpSxSy;
	  densityDensityupdown.GetMatrixElement(i, j, Tmpupdown);
	  sxSx.GetMatrixElement(i, j, TmpSxSx);
	  sySy.GetMatrixElement(i, j, TmpSySy);
	  szSz.GetMatrixElement(i, j, TmpSzSz);
	  sxSy.GetMatrixElement(i, j, TmpSxSy);
	  if ((Tmpupdown != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0) || (TmpSxSy != 0.0))
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
      this->InteractionFactorsupdownupdown = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsupupupdown = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsdowndownupdown = new Complex* [this->NbrInterSectorSums];
      this->NbrInterSectorSums = 0;
      for (int i = 0; i < densityDensityupdown.GetNbrRow(); ++i)
	{
	  for (int j = i; j < densityDensityupdown.GetNbrRow(); ++j)
	    {
	      double Tmpupdown;
	      double TmpSxSx;
	      double TmpSySy;
	      double TmpSzSz;
	      double TmpSxSy;
	      densityDensityupdown.GetMatrixElement(i, j, Tmpupdown);
	      sxSx.GetMatrixElement(i, j, TmpSxSx);
	      sySy.GetMatrixElement(i, j, TmpSySy);
	      szSz.GetMatrixElement(i, j, TmpSzSz);
	      sxSy.GetMatrixElement(i, j, TmpSxSy);
	      if ((Tmpupdown != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0) || (TmpSxSy != 0.0))
 		{
// 		  if ()
// 		    {
// 		      this->NbrInterSectorIndicesPerSum[this->NbrInterSectorSums] = 1;
// 		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums] = new int [2];
// 		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums] = new Complex [1];
// 		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums] = new Complex [1];
// 		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums] = new Complex [1];
// 		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][0] = i;
// 		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][1] = j;
// 		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][0] =  Sign * (Tmpupdown - (0.5 * TmpSzSz));
// 		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums][0] = 0.0;
// 		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums][0] = 0.0;
// 		    }
// 		  else
// 		    {
		      this->NbrInterSectorIndicesPerSum[this->NbrInterSectorSums] = 2;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums] = new int [4];
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums] = new Complex [4];
		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums] = new Complex [2];
		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums] = new Complex [2];
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][0] = i;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][1] = j;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][2] = j;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][3] = i;
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][0] =  Sign * (Tmpupdown - (0.25 * TmpSzSz));
		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums][0] = 0.0;
		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums][0] = 0.0;
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][1] =  Complex (-0.25 * Sign * (TmpSxSx + TmpSySy), 0.5 * Sign * TmpSxSy);
		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums][1] = 0.0;
		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums][1] = 0.0;
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][2] =  Complex (-0.25 * Sign * (TmpSxSx + TmpSySy), -0.5 * Sign * TmpSxSy);
// 		      cout << this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][2] << endl;
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][3] =  Sign * (Tmpupdown -0.25 * TmpSzSz);
// 		    }
		  ++this->NbrInterSectorSums;
		}
	    }
	}
    }

  
  if (this->NbrIntraSectorSums != 0)
    {
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];
      this->NbrIntraSectorSums = 0;
      for (int i = 0; i < densityDensityupup.GetNbrRow(); ++i)
	{
	  for (int j = i; j < densityDensityupup.GetNbrRow(); ++j)
	    {
	      double Tmp;
	      double Tmp1;
	      double TmpSxSx;
	      double TmpSySy;
	      double TmpSzSz;
	      double TmpSxSy;
	      densityDensityupup.GetMatrixElement(i, j, Tmp);
	      densityDensitydowndown.GetMatrixElement(i, j, Tmp1);
	      sxSx.GetMatrixElement(i, j, TmpSxSx);
	      sySy.GetMatrixElement(i, j, TmpSySy);
	      szSz.GetMatrixElement(i, j, TmpSzSz);
	      sxSy.GetMatrixElement(i, j, TmpSxSy);
	      if ((Tmp != 0.0) || (Tmp1 != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0) || (TmpSxSy != 0.0))
		{
		  this->NbrIntraSectorIndicesPerSum[this->NbrIntraSectorSums] = 1;
		  this->IntraSectorIndicesPerSum[this->NbrIntraSectorSums] = new int [2];
		  this->InteractionFactorsupupupup[this->NbrIntraSectorSums] = new Complex [1];
		  this->InteractionFactorsdowndowndowndown[this->NbrIntraSectorSums] = new Complex [1];
		  this->InteractionFactorsupupdowndown[this->NbrIntraSectorSums] = new Complex [1];
		  this->InteractionFactorsdowndownupup[this->NbrIntraSectorSums] = new Complex [1];
		  this->IntraSectorIndicesPerSum[this->NbrIntraSectorSums][0] = i;
		  this->IntraSectorIndicesPerSum[this->NbrIntraSectorSums][1] = j;
		  this->InteractionFactorsupupupup[this->NbrIntraSectorSums][0] =  Sign * (Tmp + 0.25 * TmpSzSz);  
		  this->InteractionFactorsdowndowndowndown[this->NbrIntraSectorSums][0] =  Sign * (Tmp1 + 0.25 * TmpSzSz);  
		  this->InteractionFactorsupupdowndown[this->NbrIntraSectorSums][0] = 0.25 * Sign * (TmpSxSx - TmpSySy);
		  this->InteractionFactorsdowndownupup[this->NbrIntraSectorSums][0] = 0.25 * Sign * (TmpSxSx - TmpSySy);
// 		  if ((TmpSxSx == 0.0) && (TmpSySy == 0.0))
// 		    {
// 		      this->InteractionFactorsupdownupup[this->NbrIntraSectorSums] = new Complex [1];
// 		      this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums] = new Complex [1];
// 		      this->InteractionFactorsupdownupup[this->NbrIntraSectorSums][0] = 0.0;
// 		      this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums][0] = 0.0;
// 		    }
// 		  else
// 		    {
		      this->InteractionFactorsupdownupup[this->NbrIntraSectorSums] = new Complex [2];
		      this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums] = new Complex [2];
		      this->InteractionFactorsupdownupup[this->NbrIntraSectorSums][0] = 0.0;
		      this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums][0] = 0.0;
		      this->InteractionFactorsupdownupup[this->NbrIntraSectorSums][1] = 0.0;
		      this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums][1] = 0.0;
		      //		    }
		  ++this->NbrIntraSectorSums;
		}
	    }
	}
    }
  else
    {
      this->InteractionFactorsupupupup = 0;
      this->InteractionFactorsdowndowndowndown = 0;
      this->InteractionFactorsupupdowndown = 0;
      this->InteractionFactorsdowndownupup = 0;
      this->InteractionFactorsupdownupup = 0;
      this->InteractionFactorsupdowndowndown = 0;
    }
}


// evaluate the two body interaction factors from a generic density-density interaction and a generic anisotropic Heisenberg interaction
//
// densityDensityIntra = matrix that gives the amplitude of each density-density interaction term for particles with the same spin
// densityDensityInter = matrix that gives the amplitude of each density-density interaction term for particles with opposite spins
// sxSx = matrix that gives the amplitude of each Sx_i Sx_j term
// sySy = matrix that gives the amplitude of each Sy_i Sy_j term
// szSz = matrix that gives the amplitude of each Sz_i Sz_j term
// sxSy = matrix that gives the amplitude of each Sx_i Sy_j term
// sySz = matrix that gives the amplitude of each Sy_i Sz_j term
// sxSz = matrix that gives the amplitude of each Sx_i Sz_j term

void ParticleOnLatticeWithSpinFullRealSpaceHamiltonian::EvaluateInteractionFactorsFromDensityDensityAndHeisenberg (RealSymmetricMatrix& densityDensityupup, 
														   RealSymmetricMatrix& densityDensitydowndown,
														   RealSymmetricMatrix& densityDensityupdown, RealSymmetricMatrix& sxSx,
														   RealSymmetricMatrix& sySy, RealSymmetricMatrix& szSz, RealSymmetricMatrix& sxSy, RealSymmetricMatrix& sySz, RealSymmetricMatrix& sxSz)
{
  long TotalNbrInteractionFactors = 0;
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  
  
  this->InteractionFactorsupupupup = 0;
  this->InteractionFactorsdowndowndowndown = 0;
  this->InteractionFactorsupupdowndown = 0;
  this->InteractionFactorsdowndownupup = 0;
  this->InteractionFactorsupdownupup = 0;
  this->InteractionFactorsupdowndowndown = 0;
  this->InteractionFactorsupupupdown = 0;
  this->InteractionFactorsdowndownupdown = 0;
  this->InteractionFactorsupdownupdown = 0;

  this->NbrIntraSectorSums = 0;  
  for (int i = 0; i < densityDensityupup.GetNbrRow(); ++i)
    {
      for (int j = i; j < densityDensityupup.GetNbrRow(); ++j)
	{
	  double Tmpupup;
	  double Tmpdowndown;
	  double TmpSxSx;
	  double TmpSySy;
	  double TmpSzSz;
	  double TmpSxSy;
	  double TmpSySz;
	  double TmpSxSz;
	  densityDensityupup.GetMatrixElement(i, j, Tmpupup);
	  densityDensitydowndown.GetMatrixElement(i, j, Tmpdowndown);
	  sxSx.GetMatrixElement(i, j, TmpSxSx);
	  sySy.GetMatrixElement(i, j, TmpSySy);
	  szSz.GetMatrixElement(i, j, TmpSzSz);
	  sxSy.GetMatrixElement(i, j, TmpSxSy);
	  sySz.GetMatrixElement(i, j, TmpSySz);
	  sxSz.GetMatrixElement(i, j, TmpSxSz);
	  if ((Tmpupup != 0.0) || (Tmpdowndown != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0) || (TmpSxSy != 0.0) || (TmpSySz != 0.0) || (TmpSxSz != 0.0))
	    {
	      ++this->NbrIntraSectorSums;
	    }
	}
    }
    
    
  this->NbrInterSectorSums = 0;
  for (int i = 0; i < densityDensityupdown.GetNbrRow(); ++i)
    {
      for (int j = i; j < densityDensityupdown.GetNbrRow(); ++j)
	{
	  double Tmpupdown;
	  double TmpSxSx;
	  double TmpSySy;
	  double TmpSzSz;
	  double TmpSxSy;
	  double TmpSySz;
	  double TmpSxSz;
	  densityDensityupdown.GetMatrixElement(i, j, Tmpupdown);
	  sxSx.GetMatrixElement(i, j, TmpSxSx);
	  sySy.GetMatrixElement(i, j, TmpSySy);
	  szSz.GetMatrixElement(i, j, TmpSzSz);
	  sxSy.GetMatrixElement(i, j, TmpSxSy);
	  sySz.GetMatrixElement(i, j, TmpSySz);
	  sxSz.GetMatrixElement(i, j, TmpSxSz);
	  if ((Tmpupdown != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0) || (TmpSxSy != 0.0) || (TmpSySz != 0.0) || (TmpSxSz != 0.0))
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
      this->InteractionFactorsupdownupdown = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsupupupdown = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsdowndownupdown = new Complex* [this->NbrInterSectorSums];
      this->NbrInterSectorSums = 0;
      for (int i = 0; i < densityDensityupdown.GetNbrRow(); ++i)
	{
	  for (int j = i; j < densityDensityupdown.GetNbrRow(); ++j)
	    {
	      double Tmpupdown;
	      double TmpSxSx;
	      double TmpSySy;
	      double TmpSzSz;
	      double TmpSxSy;
	      double TmpSySz;
	      double TmpSxSz;
	      densityDensityupdown.GetMatrixElement(i, j, Tmpupdown);
	      sxSx.GetMatrixElement(i, j, TmpSxSx);
	      sySy.GetMatrixElement(i, j, TmpSySy);
	      szSz.GetMatrixElement(i, j, TmpSzSz);
	      sxSy.GetMatrixElement(i, j, TmpSxSy);
	      sySz.GetMatrixElement(i, j, TmpSySz);
	      sxSz.GetMatrixElement(i, j, TmpSxSz);
	      if ((Tmpupdown != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0) || (TmpSxSy != 0.0) || (TmpSySz != 0.0) || (TmpSxSz != 0.0))
 		{
		      this->NbrInterSectorIndicesPerSum[this->NbrInterSectorSums] = 2;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums] = new int [4];
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums] = new Complex [4];
		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums] = new Complex [2];
		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums] = new Complex [2];
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][0] = i;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][1] = j;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][2] = j;
		      this->InterSectorIndicesPerSum[this->NbrInterSectorSums][3] = i;
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][0] =  Sign * (Tmpupdown - (0.25 * TmpSzSz));
		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums][0] = 0.25 * Sign * Complex (TmpSxSz, -TmpSySz);
		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums][0] = 0.25 * Sign * Complex (-TmpSxSz, -TmpSySz);
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][1] =  -0.25 * Sign * (TmpSxSx + TmpSySy);
		      this->InteractionFactorsupupupdown[this->NbrInterSectorSums][1] = 0.25 * Sign * Complex (TmpSxSz, -TmpSySz);
		      this->InteractionFactorsdowndownupdown[this->NbrInterSectorSums][1] = 0.25 * Sign * Complex (-TmpSxSz, -TmpSySz);
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][2] =  -0.25 * Sign * (TmpSxSx + TmpSySy);
// 		      cout << this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][2] << endl;
		      this->InteractionFactorsupdownupdown[this->NbrInterSectorSums][3] =  Sign * (Tmpupdown -0.25 * TmpSzSz);
// 		    }
		  ++this->NbrInterSectorSums;
		}
	    }
	}
    }

  
  if (this->NbrIntraSectorSums != 0)
    {
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];
      this->NbrIntraSectorSums = 0;
      for (int i = 0; i < densityDensityupup.GetNbrRow(); ++i)
	{
	  for (int j = i; j < densityDensityupup.GetNbrRow(); ++j)
	    {
	      double Tmp;
	      double Tmp1;
	      double TmpSxSx;
	      double TmpSySy;
	      double TmpSzSz;
	      double TmpSxSy;
	      double TmpSySz;
	      double TmpSxSz;
	      densityDensityupup.GetMatrixElement(i, j, Tmp);
	      densityDensitydowndown.GetMatrixElement(i, j, Tmp1);
	      sxSx.GetMatrixElement(i, j, TmpSxSx);
	      sySy.GetMatrixElement(i, j, TmpSySy);
	      szSz.GetMatrixElement(i, j, TmpSzSz);
	      sxSy.GetMatrixElement(i, j, TmpSxSy);
	      sySz.GetMatrixElement(i, j, TmpSySz);
	      sxSz.GetMatrixElement(i, j, TmpSxSz);
	      if ((Tmp != 0.0) || (Tmp1 != 0.0) || (TmpSxSx != 0.0) || (TmpSySy != 0.0) || (TmpSzSz != 0.0) || (TmpSxSy != 0.0))
		{
		  this->NbrIntraSectorIndicesPerSum[this->NbrIntraSectorSums] = 1;
		  this->IntraSectorIndicesPerSum[this->NbrIntraSectorSums] = new int [2];
		  this->InteractionFactorsupupupup[this->NbrIntraSectorSums] = new Complex [1];
		  this->InteractionFactorsdowndowndowndown[this->NbrIntraSectorSums] = new Complex [1];
		  this->InteractionFactorsupupdowndown[this->NbrIntraSectorSums] = new Complex [1];
		  this->InteractionFactorsdowndownupup[this->NbrIntraSectorSums] = new Complex [1];
		  this->IntraSectorIndicesPerSum[this->NbrIntraSectorSums][0] = i;
		  this->IntraSectorIndicesPerSum[this->NbrIntraSectorSums][1] = j;
		  this->InteractionFactorsupupupup[this->NbrIntraSectorSums][0] =  Sign * (Tmp + 0.25 * TmpSzSz);  
		  this->InteractionFactorsdowndowndowndown[this->NbrIntraSectorSums][0] =  Sign * (Tmp1 + 0.25 * TmpSzSz);  
		  this->InteractionFactorsupupdowndown[this->NbrIntraSectorSums][0] = 0.25 * Sign * Complex (TmpSxSx - TmpSySy, - 2.0 * TmpSxSy);
		  this->InteractionFactorsdowndownupup[this->NbrIntraSectorSums][0] = 0.25 * Sign * Complex (TmpSxSx - TmpSySy, 2.0 * TmpSxSy);
		  
		  this->InteractionFactorsupdownupup[this->NbrIntraSectorSums] = new Complex [2];
		  this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums] = new Complex [2];
		  this->InteractionFactorsupdownupup[this->NbrIntraSectorSums][0] = Sign * 0.25 * Complex (TmpSxSz, TmpSySz);
		  this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums][0] = Sign * 0.25 * Complex (-TmpSxSz, TmpSySz);
		  this->InteractionFactorsupdownupup[this->NbrIntraSectorSums][1] = Sign * 0.25 * Complex (TmpSxSz, TmpSySz);
		  this->InteractionFactorsupdowndowndown[this->NbrIntraSectorSums][1] = Sign * 0.25 * Complex (-TmpSxSz, TmpSySz);
		   
		  ++this->NbrIntraSectorSums;
		}
	    }
	}
    }
  else
    {
      this->InteractionFactorsupupupup = 0;
      this->InteractionFactorsdowndowndowndown = 0;
      this->InteractionFactorsupupdowndown = 0;
      this->InteractionFactorsdowndownupup = 0;
      this->InteractionFactorsupdownupup = 0;
      this->InteractionFactorsupdowndowndown = 0;
    }
}


// add an additional S^2 term to the Hamiltonian
//
// factor = factor in front of the S^2
// fixedSz = flag indicating whether Sz needs to be evaluated
// memory = amount of memory that can be used for S^2  precalculations

void ParticleOnLatticeWithSpinFullRealSpaceHamiltonian::AddS2 (double factor, bool fixedSz, long memory)
{
  this->S2Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceS2Hamiltonian (this->Particles, this->NbrParticles, this->NbrSites, factor, fixedSz, this->Architecture, memory); 
}
