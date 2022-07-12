////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//                   class of Hubbard model hamiltonian associated            //
//                to particles on a lattice                                   //
//                                                                            //
//                        last modification : 19/06/2014                      //
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
#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian.h"
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

ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian()
{
  this->HermitianSymmetryFlag = false;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSite = number of sites
// kineticFactor = multiplicative factor in front of the kinetic term
// uPotential = Hubbard potential strength
// bandParameter = band parameter
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSite, int** geometryDescription, int nbrBonds, 
													   double kineticFactorIsotropic, double kineticFactorAnisotropic, double uPotential, double j1Factor, double j2Factor, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSite = nbrSite;
  this->NbrBonds = 0;
  this->SitesA = 0;
  this->SitesB = 0;
  this->Bonds = 0;
  if (geometryDescription != 0)
    {
      this->NbrBonds = nbrBonds;
      this->SitesA = geometryDescription[0];
      this->SitesB = geometryDescription[1];
      this->Bonds = geometryDescription[2];
    }
  
  this->LzMax = this->NbrSite - 1;
  this->UPotential = uPotential;
  this->HamiltonianShift = 0.0;
  this->KineticFactorIsotropic = kineticFactorIsotropic;
  this->KineticFactorAnisotropic = kineticFactorAnisotropic;
  this->J1Factor = j1Factor;
  this->J2Factor = j2Factor;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->OneBodyGenericInteractionFactorsupup = 0;
  this->OneBodyGenericInteractionFactorsdowndown = 0;
  this->OneBodyGenericInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->PlotMapNearestNeighborBonds();
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

// constructor from the ecplicit the bond description
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSite = number of sites
// nbrBonds = number of bonds
// sitesA = array of A sites for each bond
// sitesB = array of B sites for each bond
// bondTypes = array that describe each type of bond (0 for x, 1 fo y, 2 for z)
// kineticFactor = multiplicative factor in front of the kinetic term
// uPotential = Hubbard potential strength
// j1Factor = strength of the isotropic nearest neighbor spin interaction
// j2Factor = strength of the anisotropic nearest neighbor spin interaction
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSite, 
													   int nbrBonds,  int* sitesA, int* sitesB, int* bondTypes, 
													   double kineticFactorIsotropic, double kineticFactorAnisotropic, 
													   double uPotential, double j1Factor, double j2Factor, 
													   AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSite = nbrSite;
  this->NbrBonds = nbrBonds;
  this->SitesA = new int [this->NbrBonds];
  this->SitesB = new int [this->NbrBonds];
  this->Bonds = new int [this->NbrBonds];
  for (int i = 0; i < this->NbrBonds; ++i)
    {
      this->SitesA[i] = sitesA[i];
      this->SitesB[i] = sitesB[i];
      this->Bonds[i] = bondTypes[i];
    }
  this->LzMax = this->NbrSite - 1;
  this->UPotential = uPotential;
  this->HamiltonianShift = 0.0;
  this->KineticFactorIsotropic = kineticFactorIsotropic;
  this->KineticFactorAnisotropic = kineticFactorAnisotropic;
  this->J1Factor = j1Factor;
  this->J2Factor = j2Factor;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->OneBodyGenericInteractionFactorsupup = 0;
  this->OneBodyGenericInteractionFactorsdowndown = 0;
  this->OneBodyGenericInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->PlotMapNearestNeighborBonds();
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

ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::~ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian()
{
  if (this->MapNearestNeighborBonds != 0)
    {
      for (int i = 0; i < this->NbrSite; ++i)
	delete[] this->MapNearestNeighborBonds[i];
      delete[] this->MapNearestNeighborBonds; 
    }
//   if (this->InteractionFactorsupup != 0)
//     {
//       for (int i = 0; i < this->NbrSite; ++i)
// 	{
// 	  delete[] this->InteractionFactorsupup[i];
// 	  delete[] this->InteractionFactorsdowndown[i];
// 	  delete[] this->InteractionFactorsupdown[i];
// 	  delete[] this->InteractionFactorsupupdowndown[i];
// 	}
//       delete[] this->InteractionFactorsupup;
//       delete[] this->InteractionFactorsupdown;
//       delete[] this->InteractionFactorsdowndown;
//       delete[] this->InteractionFactorsupupdowndown;
//     }
  if (this->OneBodyGenericInteractionFactorsupup != 0)
    {
      for (int i = 0; i < this->NbrSite; ++i)
	{
	  delete[] this->OneBodyGenericInteractionFactorsupup[i];
	  delete[] this->OneBodyGenericInteractionFactorsdowndown[i];
	  delete[] this->OneBodyGenericInteractionFactorsupdown[i];
	}
	  
      delete[] this->OneBodyGenericInteractionFactorsupup;
      delete[] this->OneBodyGenericInteractionFactorsdowndown;
      delete[] this->OneBodyGenericInteractionFactorsupdown;
    }
  
    
//   if (this->FastMultiplicationFlag == true)
//     {
//       long MinIndex;
//       long MaxIndex;
//       this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
//       int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
//       int ReducedDim = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
//       if ((ReducedDim * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
// 	++ReducedDim;
//       for (int i = 0; i < ReducedDim; ++i)
// 	{
// 	  delete[] this->InteractionPerComponentIndex[i];
// 	  delete[] this->InteractionPerComponentCoefficient[i];
// 	}
//       delete[] this->InteractionPerComponentIndex;
//       delete[] this->InteractionPerComponentCoefficient;
//       delete[] this->NbrInteractionPerComponent;
//       this->FastMultiplicationFlag = false;
//     }
}
  


// evaluate all nearest neighbor interaction factors
//   

void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  
  this->NbrInterSectorSums = 0;
  for (int j = 0; j < this->NbrSite; ++j)
    {
      for (int k = 0; k < 3; ++k)
	if (this->MapNearestNeighborBonds[j][k] < this->NbrSite)
	  ++this->NbrInterSectorSums;
      ++this->NbrInterSectorSums;      
    }
   
  
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  
  
  this->NbrIntraSectorSums = this->NbrInterSectorSums;
  this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    this->NbrIntraSectorIndicesPerSum[i] = 0;      
  
  int index = 0;
  for (int j = 0; j < this->NbrSite; ++j)
    {
      for (int k = 0; k < 3; ++k)  
	if (this->MapNearestNeighborBonds[j][k] < this->NbrSite)
	  {
	    this->NbrInterSectorIndicesPerSum[index] += 2;
	    ++index;
	  }
      ++this->NbrInterSectorIndicesPerSum[index];
      ++index;
    }
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      if (this->NbrInterSectorIndicesPerSum[i] > 0)
	{
	  this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
	  this->NbrInterSectorIndicesPerSum[i] = 0;
	}
    }
  int TmpSum = 0;
  for (int j = 0; j < this->NbrSite; ++j)
    {
      for (int k = 0; k < 3; ++k)
	if (this->MapNearestNeighborBonds[j][k] < this->NbrSite)
	  {
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = j;
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = this->MapNearestNeighborBonds[j][k];
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	    
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = this->MapNearestNeighborBonds[j][k];
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = j;
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	    ++TmpSum;
	  }
      this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = j;
      this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = j;
      ++this->NbrInterSectorIndicesPerSum[TmpSum];    
      ++TmpSum;      
    }
 
  TmpSum = 0;
  for (int j = 0; j < this->NbrSite; ++j)
    for (int k = 0; k < 3; ++k)
      {
	int Index1 = j;
	int Index2 = this->MapNearestNeighborBonds[j][k];
	if (Index2 < this->NbrSite)
	{
	  if (Index1 < Index2)
	    ++this->NbrIntraSectorIndicesPerSum[TmpSum];
	  ++TmpSum;
	}
      }
  this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    {
      if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
    }
  TmpSum = 0;
  for (int j = 0;j < this->NbrSite; ++j)
    for (int k = 0; k < 3; ++k)
      {
	int Index1 = j;
	int Index2 = this->MapNearestNeighborBonds[j][k];
	if (Index2 < this->NbrSite)
	  {
	    if (Index1 < Index2)
	      {
		this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		++this->NbrIntraSectorIndicesPerSum[TmpSum];    
	      }
	    ++TmpSum;
	  }
      }
  this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
  this->InteractionFactorsupupupdown = new Complex* [this->NbrInterSectorSums];
  this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
  this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
  this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];
  this->InteractionFactorsdowndownupdown = new Complex* [this->NbrInterSectorSums];
  this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
  this->InteractionFactorsupdownupdown = new Complex* [this->NbrInterSectorSums];
  this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];
  
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    {
      this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
      this->InteractionFactorsdowndowndowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
      this->InteractionFactorsupupdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
      this->InteractionFactorsdowndownupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
      this->InteractionFactorsupdownupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
      this->InteractionFactorsupdowndowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
      int Index = 0;
      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	{
	  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
	      if (this->FindBondType(Index1, Index2) == 0)
		{
		  this->InteractionFactorsupupupup[i][Index] = -(this->J1Factor - this->J2Factor);
		  this->InteractionFactorsdowndowndowndown[i][Index] = -(this->J1Factor - this->J2Factor);
		  this->InteractionFactorsdowndownupup[i][Index] = -2.0 * this->J2Factor;
		  this->InteractionFactorsupupdowndown[i][Index] = -2.0 * this->J2Factor;
		}
	      
	      if (this->FindBondType(Index1, Index2) == 1)
		{
		  this->InteractionFactorsupupupup[i][Index] = -(this->J1Factor - this->J2Factor);
		  this->InteractionFactorsdowndowndowndown[i][Index] = -(this->J1Factor - this->J2Factor);
		  this->InteractionFactorsdowndownupup[i][Index] = 2.0 * this->J2Factor;
		  this->InteractionFactorsupupdowndown[i][Index] = 2.0 * this->J2Factor;
		}
	      
	      if (this->FindBondType(Index1, Index2) == 2)
		{
		  this->InteractionFactorsupupupup[i][Index] = -(this->J1Factor + this->J2Factor);
		  this->InteractionFactorsdowndowndowndown[i][Index] = -(this->J1Factor + this->J2Factor);
		  this->InteractionFactorsdowndownupup[i][Index] = 0.0;
		  this->InteractionFactorsupupdowndown[i][Index] = 0.0;
		}	      	      
// 	      cout << Index1 << " " << Index2 << " " << Index3 << " " << Index4 << " t=" << this->FindBondType(Index1, Index2)
// 	      << endl;
// 		   << " : uuuu=" << this->InteractionFactorsupupupup[i][Index]
// 		   << " dddd=" << this->InteractionFactorsdowndowndowndown[i][Index]
// 		   << " uudd=" << this->InteractionFactorsupupdowndown[i][Index]
// 		   << " dduu=" << this->InteractionFactorsdowndownupup[i][Index] << endl;
	      TotalNbrInteractionFactors += 4;
	      ++Index;
	    }
	}
      Index = 0;
      for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	{
	  for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
	    {
	      this->InteractionFactorsupdownupup[i][Index] = 0.0;
	      this->InteractionFactorsupdowndowndown[i][Index] = 0.0;
	      TotalNbrInteractionFactors += 2;	      
	      ++Index;
	    }
	}
    }
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      this->InteractionFactorsupdownupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
      this->InteractionFactorsupupupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
      this->InteractionFactorsdowndownupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
      int Index = 0;
      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	{
	  int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	  int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	  if (Index1 != Index2)	      
	    {
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  
// 		  cout << Index1 << " " << Index2 << " " << Index3 << " " << Index4 << endl;
		  if (Index3 != Index4)
		    {
		      
		      if (this->FindBondType(Index1, Index2) == 0)
			{
			  if (Index1 == Index3)
			    this->InteractionFactorsupdownupdown[i][Index] = 0.5 * (this->J1Factor - this->J2Factor);
			  else
			    this->InteractionFactorsupdownupdown[i][Index] = 0.5 * 2.0 * this->J1Factor;
			}
		      if (this->FindBondType(Index1, Index2) == 1)
			{
			  if (Index1 == Index3)
			    this->InteractionFactorsupdownupdown[i][Index] = 0.5 * (this->J1Factor - this->J2Factor);
			  else
			    this->InteractionFactorsupdownupdown[i][Index] = 0.5 * 2.0 * this->J1Factor;
			}
		      if (this->FindBondType(Index1, Index2) == 2)
			{		    
			  if (Index1 == Index3)
			    this->InteractionFactorsupdownupdown[i][Index] = 0.5 * (this->J1Factor + this->J2Factor);  
			  else
			    this->InteractionFactorsupdownupdown[i][Index] = 0.5 * 2.0 * (this->J1Factor - this->J2Factor);
			}
// 		      cout << Index1 << " " << Index2 << " " << Index3 << " " << Index4 << " : " << this->InteractionFactorsupdownupdown[i][Index] << endl;
		      TotalNbrInteractionFactors += 1;
		      ++Index;
		    }
		  else
		    {
		      this->InteractionFactorsupdownupdown[i][Index] = 0.0;    
		      TotalNbrInteractionFactors += 1;
		      ++Index;
		    }
		}
	    }
	  else
	    {
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  if (Index3 == Index4)
		    {
		      this->InteractionFactorsupdownupdown[i][Index] = -this->UPotential;
		    }
		  else
		    {
		      this->InteractionFactorsupdownupdown[i][Index] = 0.0;
		    }
		  TotalNbrInteractionFactors += 1;
		  ++Index;
		}
	    }
	}
      Index = 0;
      for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	{
	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
	    {
	      this->InteractionFactorsupupupdown[i][Index] = 0.0;
	      this->InteractionFactorsdowndownupdown[i][Index] = 0.0;
	      TotalNbrInteractionFactors += 2;	      
	      ++Index;
	    }
	}
    }
  
  
  this->OneBodyGenericInteractionFactorsupup = new double*[this->NbrSite];
  this->OneBodyGenericInteractionFactorsdowndown = new double*[this->NbrSite];
  this->OneBodyGenericInteractionFactorsupdown = new Complex*[this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
    {
      this->OneBodyGenericInteractionFactorsupup[i] = new double [3];
      this->OneBodyGenericInteractionFactorsdowndown[i] = new double [3];
      this->OneBodyGenericInteractionFactorsupdown[i] = new Complex [3];
      for (int j = 0; j < 3; ++j)
	{   
	  this->OneBodyGenericInteractionFactorsupup[i][j] = 0.0;
	  this->OneBodyGenericInteractionFactorsdowndown[i][j] = 0.0;
	  this->OneBodyGenericInteractionFactorsupdown[i][j] = 0.0;
	}
    }
  
  for (int i = 0; i < this->NbrSite; ++i)
    {
      int jx = this->MapNearestNeighborBonds[i][0];
      int jy = this->MapNearestNeighborBonds[i][1];
      int jz = this->MapNearestNeighborBonds[i][2];
      
      if (jx < this->NbrSite)
	{
	  this->OneBodyGenericInteractionFactorsupup[i][0] += -this->KineticFactorIsotropic;
	  this->OneBodyGenericInteractionFactorsdowndown[i][0] += -this->KineticFactorIsotropic;
	  this->OneBodyGenericInteractionFactorsupdown[i][0] += Complex(- this->KineticFactorAnisotropic, 0.0);
	}
      
      if (jy < this->NbrSite)
	{
	  this->OneBodyGenericInteractionFactorsupup[i][1] += -this->KineticFactorIsotropic;
	  this->OneBodyGenericInteractionFactorsdowndown[i][1] += -this->KineticFactorIsotropic;
	  this->OneBodyGenericInteractionFactorsupdown[i][1] += Complex(0.0,  this->KineticFactorAnisotropic);
	}
      
      if (jz < this->NbrSite)
	{
 	  this->OneBodyGenericInteractionFactorsupup[i][2] += -this->KineticFactorIsotropic - this->KineticFactorAnisotropic;
 	  this->OneBodyGenericInteractionFactorsdowndown[i][2] += -this->KineticFactorIsotropic + this->KineticFactorAnisotropic;	  
	}
      //    cout << i << " " << this->InteractionFactorsupup[i][0] << " " << this->InteractionFactorsupup[i][1] << " " << this->InteractionFactorsupup[i][2] << endl;
    }
  
//   for (int i = 0; i < this->NbrIntraSectorSums; ++i)
//   {
//     int Index = 0;
//       for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
// 	{
// 	  int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
// 	  int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
// 	  for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
// 	    {
// 	      int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
// 	      int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
// 	      
// 	      cout << Index1 << " " << Index2 << " " << this->FindBondType(Index1, Index2) << " " << this->InteractionFactorsupupupup[i][Index] << " " << this->InteractionFactorsdowndowndowndown[i][Index] << " " << this->InteractionFactorsupupdowndown[i][Index] << " " << this->InteractionFactorsdowndownupup[i][Index] << endl;
// 	      ++Index;
// 	    }
// 	}
//   }
  
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

//fill up the positions of the three types of nearest neighbor bonds
//

void ParticleOnLatticeWithSpinKitaevHeisenbergHamiltonian::PlotMapNearestNeighborBonds()
{
  int Index;
  this->MapNearestNeighborBonds = new int* [this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
    this->MapNearestNeighborBonds[i] = new int [3];
  
  if (this->Bonds == 0)
    {
      if ((this->NbrSite % 4) == 2)
	{
	  for (int i = 0; i <= (this->NbrSite - 2)/4 ; ++i)
	    {
	      Index = 4*i;
	      if (i == 0)
		this->MapNearestNeighborBonds[Index][0] = this->NbrSite;
	      else
		this->MapNearestNeighborBonds[Index][0] = Index - 2;
	      if (i < (this->NbrSite - 2)/4)
		this->MapNearestNeighborBonds[Index][1] = Index + 2;
	      else
		this->MapNearestNeighborBonds[Index][1] = this->NbrSite;
	      this->MapNearestNeighborBonds[Index][2] = Index + 1;
	      
	      Index = 4*i + 1;
	      if (i < (this->NbrSite - 2)/4)
		this->MapNearestNeighborBonds[Index][0] = Index + 2;
	      else
		this->MapNearestNeighborBonds[Index][0] = this->NbrSite;
	      if (i == 0)
		this->MapNearestNeighborBonds[Index][1] = this->NbrSite;
	      else
		this->MapNearestNeighborBonds[Index][1] = Index - 2;
	      this->MapNearestNeighborBonds[Index][2] = Index - 1;
	      
	      if (i < (this->NbrSite - 2)/4)
		{
		  Index = 4*i + 2;
		  this->MapNearestNeighborBonds[Index][0] = Index + 2;
		  this->MapNearestNeighborBonds[Index][1] = Index - 2;
		  this->MapNearestNeighborBonds[Index][2] = this->NbrSite;
      
		  Index = 4*i + 3;
		  this->MapNearestNeighborBonds[Index][0] = Index - 2;
		  this->MapNearestNeighborBonds[Index][1] = Index + 2;
		  this->MapNearestNeighborBonds[Index][2] = this->NbrSite;      
		}
	    }
	}
      else
	{
	  if (((this->NbrSite % 4) == 0) && (this->NbrSite > 4))
	    {
	      for (int i = 0; i < this->NbrSite /4 ; ++i)
		{
		  Index = 4*i;
		  if (i == 0)
		    this->MapNearestNeighborBonds[Index][0] = this->NbrSite - 2;
		  else
		    this->MapNearestNeighborBonds[Index][0] = Index - 2;
		  this->MapNearestNeighborBonds[Index][1] = Index + 2;
		  this->MapNearestNeighborBonds[Index][2] = Index + 1;
		  
		  Index = 4*i + 1;
		  this->MapNearestNeighborBonds[Index][0] = Index + 2;
		  if (i == 0)
		    this->MapNearestNeighborBonds[Index][1] = this->NbrSite - 1;
		  else
		    this->MapNearestNeighborBonds[Index][1] = Index - 2;
		  this->MapNearestNeighborBonds[Index][2] = Index - 1;
		  
		  
		  Index = 4*i + 2;
		  if (i == (this->NbrSite /4 - 1))
		    this->MapNearestNeighborBonds[Index][0] = 0;
		  else
		    this->MapNearestNeighborBonds[Index][0] = Index + 2;
		  this->MapNearestNeighborBonds[Index][1] = Index - 2;
		  this->MapNearestNeighborBonds[Index][2] = this->NbrSite;
		  
		  Index = 4*i + 3;
		  this->MapNearestNeighborBonds[Index][0] = Index - 2;
		  if (i == (this->NbrSite /4 - 1))
		    this->MapNearestNeighborBonds[Index][1] = 1;
		  else
		    this->MapNearestNeighborBonds[Index][1] = Index + 2;
		  this->MapNearestNeighborBonds[Index][2] = this->NbrSite;      
		  
		}
	    }
	  else
	    {
	      if (this->NbrSite == 4)
		{
		  this->MapNearestNeighborBonds[0][0] = this->NbrSite;
		  this->MapNearestNeighborBonds[0][1] = 1;
		  this->MapNearestNeighborBonds[0][2] = this->NbrSite;
		  this->MapNearestNeighborBonds[1][0] = 2;
		  this->MapNearestNeighborBonds[1][1] = 0;
		  this->MapNearestNeighborBonds[1][2] = this->NbrSite;
		  this->MapNearestNeighborBonds[2][0] = 1;
		  this->MapNearestNeighborBonds[2][1] = this->NbrSite;
		  this->MapNearestNeighborBonds[2][2] = 3;
		  this->MapNearestNeighborBonds[3][0] = this->NbrSite;
		  this->MapNearestNeighborBonds[3][1] = this->NbrSite;
		  this->MapNearestNeighborBonds[3][2] = 2;
		}
	      else
		{
		  for (int i = 0; i < this->NbrSite; ++i) 
		    {
		      if ((i % 2) != 0)
			{
			  if (i < this->NbrSite - 1)
			    this->MapNearestNeighborBonds[i][0] = i + 1;
			  else
			    this->MapNearestNeighborBonds[i][0] = this->NbrSite;
			  
			  this->MapNearestNeighborBonds[i][1] = i - 1;
			}
		      else
			{
			  if (i > 0)
			    this->MapNearestNeighborBonds[i][0] = i - 1;
			  else
			    this->MapNearestNeighborBonds[i][0] = this->NbrSite; 
			  if (i < this->NbrSite - 1)
			    this->MapNearestNeighborBonds[i][1] = i + 1;
			  else
			    this->MapNearestNeighborBonds[i][1] = this->NbrSite;
			}
		      this->MapNearestNeighborBonds[i][2] = this->NbrSite;
		      
		      //   cout << i << " " << this->MapNearestNeighborBonds[i][0] << " " << this->MapNearestNeighborBonds[i][1] << " " << this->MapNearestNeighborBonds[i][2] << endl;
		    }  
		}
	    }
	}
    }
  else
    {
      for (int i = 0; i < this->NbrSite; ++i)
	for (int k = 0; k < 3; ++k)
	  this->MapNearestNeighborBonds[i][k] = this->NbrSite;
      for (int i = 0; i < this->NbrSite; ++i)
	{
	  for (int k = 0; k < 3; ++k)
	    {
	      for (int j = 0; ((j < this->NbrBonds) && (this->MapNearestNeighborBonds[i][k] == this->NbrSite)); ++j)
		{
		  if ((this->SitesA[j] == i) && (this->Bonds[j] == k))
		    {
		      this->MapNearestNeighborBonds[i][k] = this->SitesB[j];
		      this->MapNearestNeighborBonds[this->SitesB[j]][k] = i;
		      // 	    cout << this->SitesA[j] << " " << this->SitesB[j] << " " << this->Bonds[j] << " " << k << endl;
		    }
		}
	    }
	}
    }
}



