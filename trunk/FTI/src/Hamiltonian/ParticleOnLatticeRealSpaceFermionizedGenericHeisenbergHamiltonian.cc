////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                class of generic heisenberg hamiltonian written             //
//                             in fermionic language                          //
//                                                                            //
//                        last modification : 09/08/2015                      //
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
#include "Hamiltonian/ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian.h"
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

ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian::ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian()
{
  this->DiagonalElements = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// parity = parity of the particle number (0 for even or 1 for odd)
// jxx = coupling in front of the SxSx term
// jyy = coupling in front of the SySy term
// jzz = coupling in front of the SzSz term
// hField = array that gives the on-site magnetic field (zero pointer if none)
// periodicBoundaryConditionFlag = true if periodic bounday conditions should be used
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian::ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian(ParticleOnSphere* particles,
																     int parity, double jxx, double jyy, double jzz, 
																     double* hField, bool periodicBoundaryConditionFlag,
																     AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrSites = this->Particles->GetNbrOrbitals();    
  this->Jxx = jxx;
  this->Jyy = jyy;
  this->Jzz = jzz;
  this->PeriodicBoundaryConditionFlag = periodicBoundaryConditionFlag;
  this->LzMax = this->NbrSites - 1;
  this->Parity = parity;
  if (this->PeriodicBoundaryConditionFlag== true)
    {
      this->HamiltonianShift = 0.25 * this->Jzz * this->NbrSites;
    }
  else
    {
      this->HamiltonianShift = 0.25 * this->Jzz * (this->NbrSites - 1);
    }
  this->HField = new double [this->NbrSites];
  if (hField == 0)
    {
      for (int i = 0; i < this->NbrSites; ++i)
	this->HField[i] = 0.0;      
    }
  else
    {
      for (int i = 0; i < this->NbrSites; ++i)
	{
	  this->HField[i] = hField[i];      
	  this->HamiltonianShift += 0.5 * hField[i];
	}
    }

  this->DiagonalElements = 0;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyTermFlag = false;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = true;
  
  this->EvaluateOneBodyFactors();
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

ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian::~ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian()
{
  if (this->DiagonalElements != 0)
    delete[] this->DiagonalElements;
  if (this->OneBodyGenericNbrConnectedSites != 0)
    {
      for (int i =0 ; i < this->NbrSites; ++i)
	{
	  if (this->OneBodyGenericNbrConnectedSites[i] > 0)
	    {
	      delete[] this->OneBodyGenericConnectedSites[i];
	      delete[] this->OneBodyGenericInteractionFactors[i];
	    }
	}
      delete[] this->OneBodyGenericNbrConnectedSites;
      delete[] this->OneBodyGenericConnectedSites;
      delete[] this->OneBodyGenericInteractionFactors;
    }

  if (this->OneBodyPairingGenericNbrConnectedSites != 0)
    {
      for (int i =0 ; i < this->NbrSites; ++i)
	{
	  if (this->OneBodyPairingGenericNbrConnectedSites[i] > 0)
	    {
	      delete[] this->OneBodyPairingGenericConnectedSites[i];
	      delete[] this->OneBodyPairingGenericInteractionFactors[i];
	    }
	}
      delete[] this->OneBodyPairingGenericNbrConnectedSites;
      delete[] this->OneBodyPairingGenericConnectedSites;
      delete[] this->OneBodyPairingGenericInteractionFactors;
    }


}
  


// evaluate the one body interaction factors from a tight-binding matrix
//
// tightBinding = hamiltonian corresponding to the tight-binding model in real space

void ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian::EvaluateOneBodyFactors()
{
   this->OneBodyGenericNbrConnectedSites = new int [this->NbrSites];
   this->OneBodyGenericConnectedSites = new int* [this->NbrSites];
   this->OneBodyGenericInteractionFactors = new double* [this->NbrSites];
   this->OneBodyPairingGenericNbrConnectedSites = new int [this->NbrSites];
   this->OneBodyPairingGenericConnectedSites = new int* [this->NbrSites];
   this->OneBodyPairingGenericInteractionFactors = new double* [this->NbrSites];
   double Tmp;
   double Sign = 1.0;
   if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
     Sign = -1.0;

  if (this->PeriodicBoundaryConditionFlag == true)
    {
      this->OneBodyGenericNbrConnectedSites[0] = 3;
      this->OneBodyGenericConnectedSites[0] = new int [this->OneBodyGenericNbrConnectedSites[0]];
      this->OneBodyGenericInteractionFactors[0] = new double [this->OneBodyGenericNbrConnectedSites[0]];
      this->OneBodyGenericConnectedSites[0][0] = 1;
      this->OneBodyGenericInteractionFactors[0][0] = Sign * 0.25 * (this->Jxx + this->Jyy);
      this->OneBodyGenericConnectedSites[0][1] = this->NbrSites - 1;
      if (this->Parity == 0)
	this->OneBodyGenericInteractionFactors[0][1] = -Sign * 0.25 * (this->Jxx + this->Jyy);
      else
	this->OneBodyGenericInteractionFactors[0][1] = Sign * 0.25 * (this->Jxx + this->Jyy);
      this->OneBodyGenericConnectedSites[0][2] = 0;
      this->OneBodyGenericInteractionFactors[0][2] = this->Jzz + this->HField[0];

      this->OneBodyPairingGenericNbrConnectedSites[0] = 1;
      this->OneBodyPairingGenericConnectedSites[0] = new int [this->OneBodyPairingGenericNbrConnectedSites[0]];
      this->OneBodyPairingGenericInteractionFactors[0] = new double [this->OneBodyPairingGenericNbrConnectedSites[0]];
      this->OneBodyPairingGenericConnectedSites[0][0] = 1;
      this->OneBodyPairingGenericInteractionFactors[0][0] = Sign * 0.25 * (this->Jxx - this->Jyy);
//       this->OneBodyPairingGenericConnectedSites[0][1] = this->NbrSites - 1;
//       if (this->Parity == 0)
// 	this->OneBodyPairingGenericInteractionFactors[0][1] = -Sign * 0.25 * (this->Jxx - this->Jyy);
//       else
// 	this->OneBodyPairingGenericInteractionFactors[0][1] = Sign * 0.25 * (this->Jxx - this->Jyy);
    }
  else
    {
      this->OneBodyGenericNbrConnectedSites[0] = 2;
      this->OneBodyGenericConnectedSites[0] = new int [this->OneBodyGenericNbrConnectedSites[0]];
      this->OneBodyGenericInteractionFactors[0] = new double [this->OneBodyGenericNbrConnectedSites[0]];
      this->OneBodyGenericConnectedSites[0][0] = 1;
      this->OneBodyGenericInteractionFactors[0][0] = Sign * 0.25 * (this->Jxx + this->Jyy);
      this->OneBodyGenericConnectedSites[0][1] = 0;
      this->OneBodyGenericInteractionFactors[0][1] = (0.5 * this->Jzz) + this->HField[0];

      this->OneBodyPairingGenericNbrConnectedSites[0] = 1;
      this->OneBodyPairingGenericConnectedSites[0] = new int [this->OneBodyPairingGenericNbrConnectedSites[0]];
      this->OneBodyPairingGenericInteractionFactors[0] = new double [this->OneBodyPairingGenericNbrConnectedSites[0]];
      this->OneBodyPairingGenericConnectedSites[0][0] = 1;
      this->OneBodyPairingGenericInteractionFactors[0][0] = Sign * 0.25 * (this->Jxx - this->Jyy);
    }


  for (int i = 1; i < (this->NbrSites - 1); ++i)
     {
       this->OneBodyGenericNbrConnectedSites[i] = 3;
       this->OneBodyGenericConnectedSites[i] = new int [this->OneBodyGenericNbrConnectedSites[i]];
       this->OneBodyGenericInteractionFactors[i] = new double [this->OneBodyGenericNbrConnectedSites[i]];
       this->OneBodyGenericConnectedSites[i][0] = i - 1;
       this->OneBodyGenericInteractionFactors[i][0] = Sign * 0.25 * (this->Jxx + this->Jyy);
       this->OneBodyGenericConnectedSites[i][1] = i + 1;
       this->OneBodyGenericInteractionFactors[i][1] = Sign * 0.25 * (this->Jxx + this->Jyy);
       this->OneBodyGenericConnectedSites[i][2] = i;
       this->OneBodyGenericInteractionFactors[i][2] = this->Jzz + this->HField[i];

       this->OneBodyPairingGenericNbrConnectedSites[i] = 1;
       this->OneBodyPairingGenericConnectedSites[i] = new int [this->OneBodyPairingGenericNbrConnectedSites[i]];
       this->OneBodyPairingGenericInteractionFactors[i] = new double [this->OneBodyPairingGenericNbrConnectedSites[i]];
//        this->OneBodyPairingGenericConnectedSites[i][0] = i - 1;
//        this->OneBodyPairingGenericInteractionFactors[i][0] = Sign * 0.25 * (this->Jxx - this->Jyy);
       this->OneBodyPairingGenericConnectedSites[i][0] = i + 1;
       this->OneBodyPairingGenericInteractionFactors[i][0] = Sign * 0.25 * (this->Jxx - this->Jyy);
    }


  if (this->PeriodicBoundaryConditionFlag == true)
    {
      this->OneBodyGenericNbrConnectedSites[this->NbrSites - 1] = 3;
      this->OneBodyGenericConnectedSites[this->NbrSites - 1] = new int [this->OneBodyGenericNbrConnectedSites[this->NbrSites - 1]];
      this->OneBodyGenericInteractionFactors[this->NbrSites - 1] = new double [this->OneBodyGenericNbrConnectedSites[this->NbrSites - 1]];
      this->OneBodyGenericConnectedSites[this->NbrSites - 1][0] = this->NbrSites - 2;
      this->OneBodyGenericInteractionFactors[this->NbrSites - 1][0] = Sign * 0.25 * (this->Jxx + this->Jyy);
      this->OneBodyGenericConnectedSites[this->NbrSites - 1][1] = 0;
      if (this->Parity == 0)
	this->OneBodyGenericInteractionFactors[this->NbrSites - 1][1] = -Sign * 0.25 * (this->Jxx + this->Jyy);
      else
	this->OneBodyGenericInteractionFactors[this->NbrSites - 1][1] = Sign * 0.25 * (this->Jxx + this->Jyy);
      this->OneBodyGenericConnectedSites[this->NbrSites - 1][2] = this->NbrSites - 1;
      this->OneBodyGenericInteractionFactors[this->NbrSites - 1][2] = this->Jzz + this->HField[this->NbrSites - 1];
 
      this->OneBodyPairingGenericNbrConnectedSites[this->NbrSites - 1] = 1;
      this->OneBodyPairingGenericConnectedSites[this->NbrSites - 1] = new int [this->OneBodyPairingGenericNbrConnectedSites[this->NbrSites - 1]];
      this->OneBodyPairingGenericInteractionFactors[this->NbrSites - 1] = new double [this->OneBodyPairingGenericNbrConnectedSites[this->NbrSites - 1]];
      this->OneBodyPairingGenericConnectedSites[this->NbrSites - 1][0] = 0;
      if (this->Parity == 0)
	this->OneBodyPairingGenericInteractionFactors[this->NbrSites - 1][0] = -Sign * 0.25 * (this->Jxx - this->Jyy);
      else
	this->OneBodyPairingGenericInteractionFactors[this->NbrSites - 1][0] = Sign * 0.25 * (this->Jxx - this->Jyy);
//       this->OneBodyPairingGenericNbrConnectedSites[this->NbrSites - 1] = 2;
//       this->OneBodyPairingGenericConnectedSites[this->NbrSites - 1] = new int [this->OneBodyPairingGenericNbrConnectedSites[this->NbrSites - 1]];
//       this->OneBodyPairingGenericInteractionFactors[this->NbrSites - 1] = new double [this->OneBodyPairingGenericNbrConnectedSites[this->NbrSites - 1]];
//       this->OneBodyPairingGenericConnectedSites[this->NbrSites - 1][0] = this->NbrSites - 2;
//       this->OneBodyPairingGenericInteractionFactors[this->NbrSites - 1][0] = Sign * 0.25 * (this->Jxx - this->Jyy);
//       this->OneBodyPairingGenericConnectedSites[this->NbrSites - 1][1] = 0;
//       if (this->Parity == 0)
// 	this->OneBodyPairingGenericInteractionFactors[this->NbrSites - 1][1] = -Sign * 0.25 * (this->Jxx - this->Jyy);
//       else
// 	this->OneBodyPairingGenericInteractionFactors[this->NbrSites - 1][1] = Sign * 0.25 * (this->Jxx - this->Jyy);
    }
  else
    {
      this->OneBodyGenericNbrConnectedSites[this->NbrSites - 1] = 2;
      this->OneBodyGenericConnectedSites[this->NbrSites - 1] = new int [this->OneBodyGenericNbrConnectedSites[this->NbrSites - 1]];
      this->OneBodyGenericInteractionFactors[this->NbrSites - 1] = new double [this->OneBodyGenericNbrConnectedSites[this->NbrSites - 1]];
      this->OneBodyGenericConnectedSites[this->NbrSites - 1][0] = this->NbrSites - 2;
      this->OneBodyGenericInteractionFactors[this->NbrSites - 1][0] = Sign * 0.25 * (this->Jxx + this->Jyy);
      this->OneBodyGenericConnectedSites[this->NbrSites - 1][1] = this->NbrSites - 1;
      this->OneBodyGenericInteractionFactors[this->NbrSites - 1][1] = (0.5 * this->Jzz) + this->HField[this->NbrSites - 1];

      this->OneBodyPairingGenericNbrConnectedSites[this->NbrSites - 1] = 0;
      this->OneBodyPairingGenericConnectedSites[this->NbrSites - 1] = 0;
      this->OneBodyPairingGenericInteractionFactors[this->NbrSites - 1] = 0;
//       this->OneBodyPairingGenericNbrConnectedSites[this->NbrSites - 1] = 1;
//       this->OneBodyPairingGenericConnectedSites[this->NbrSites - 1] = new int [this->OneBodyPairingGenericNbrConnectedSites[this->NbrSites - 1]];
//       this->OneBodyPairingGenericInteractionFactors[this->NbrSites - 1] = new double [this->OneBodyPairingGenericNbrConnectedSites[this->NbrSites - 1]];
//       this->OneBodyPairingGenericConnectedSites[this->NbrSites - 1][0] = this->NbrSites - 2;
//       this->OneBodyPairingGenericInteractionFactors[this->NbrSites - 1][0] = Sign * 0.25 * (this->Jxx - this->Jyy);
    }
}

// evaluate the two body interaction factors from a generic density-density interaction
//
// densityDensity = matrix that gives the amplitude of each density-density interaction term

void ParticleOnLatticeRealSpaceFermionizedGenericHeisenbergHamiltonian::EvaluateInteractionFactors ()
{
  if (this->PeriodicBoundaryConditionFlag == true)
    this->NbrSectorSums = this->NbrSites;
  else
    this->NbrSectorSums = this->NbrSites - 1;  
  
  this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
  this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
  this->InteractionFactors = new double* [this->NbrSectorSums];
  for (int i = 0; i < this->NbrSectorSums; ++i)
    {
      this->NbrSectorIndicesPerSum[i] = 1;
      this->SectorIndicesPerSum[i] = new int [2];
      this->InteractionFactors[i] = new double [1];
    }

  this->NbrSectorSums = 0;
  double Sign = 1.0;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;
  for (int i = 1; i < this->NbrSites; ++i)
    {
      this->SectorIndicesPerSum[this->NbrSectorSums][0] = i - 1;
      this->SectorIndicesPerSum[this->NbrSectorSums][1] = i;
      this->InteractionFactors[this->NbrSectorSums][0] =  Sign * this->Jzz;
      ++this->NbrSectorSums;
    }
  if (this->PeriodicBoundaryConditionFlag == true)
    {
       this->SectorIndicesPerSum[this->NbrSectorSums][0] = this->NbrSites - 1;
       this->SectorIndicesPerSum[this->NbrSectorSums][1] = 0;
       this->InteractionFactors[this->NbrSectorSums][0] =  Sign * this->Jzz;
       ++this->NbrSectorSums;
   } 
}
  
