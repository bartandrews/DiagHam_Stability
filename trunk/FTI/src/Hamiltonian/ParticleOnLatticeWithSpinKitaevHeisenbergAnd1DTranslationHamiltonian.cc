////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//                   class of Hubbard hamiltonian associated                  //
//  to particles on a lattice with spin and translation in the x direction    //
//                                                                            //
//                        last modification : 13/07/2014                      //
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
#include "Hamiltonian/ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian.h"
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

ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian()
{
  this->XMomentum = 0;
  this->MaxMomentum = 1;
  this->XIncrement = 1;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites
// momentum = momentum sector
// periodicity = periodicity with respect to site numbering 
// kineticFactor = multiplicative factor in front of the kinetic term
// uPotential = Hubbard potential strength
// bandParameter = band parameter
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSite, int momentum, int periodicity, char* geometryFile, double kineticFactorIsotropic, double kineticFactorAnisotropic, double uPotential, double j1Factor, double j2Factor, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSite = nbrSite;
  this->XMomentum = momentum;
  this->XIncrement = periodicity;
  this->MaxMomentum = this->XMomentum / this->XIncrement;
  this->NbrBonds = 0;
  this->SitesA = 0;
  this->SitesB = 0;
  this->Bonds = 0;
  if (geometryFile != 0)
    {
      MultiColumnASCIIFile LatticeFile;
      if (LatticeFile.Parse(geometryFile) == false)
	{
	  LatticeFile.DumpErrors(cout);
	}
      this->NbrBonds = LatticeFile.GetNbrLines();
      this->SitesA = LatticeFile.GetAsIntegerArray(0);
      this->SitesB = LatticeFile.GetAsIntegerArray(1);
      this->Bonds = LatticeFile.GetAsIntegerArray(2);
    }
    
  this->LzMax = this->NbrSite - 1;
  this->UPotential = uPotential; //2.0 * uPotential / ((double) this->NbrParticles);
  this->HamiltonianShift = 0.0;//4.0 * uPotential;
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
  
  this->EvaluateExponentialFactors();
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

// constructor from the explicit the bond description
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSites = number of sites
// nbrBonds = number of bonds
// sitesA = array of A sites for each bond
// sitesB = array of B sites for each bond
// bondTypes = array that describe each type of bond (0 for x, 1 fo y, 2 for z)
// momentum = momentum sector
// periodicity = periodicity with respect to site numbering 
// kineticFactorIntra = multiplicative factor in front of the intraspin kinetic term
// uPotential = Hubbard potential strength
// j1Factor = strength of the isotropic nearest neighbor spin interaction
// j2Factor = strength of the anisotropic nearest neighbor spin interaction
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSite, 
																	   int nbrBonds,  int* sitesA, int* sitesB, int* bondTypes,
																	   int momentum, int periodicity, 
																	   double kineticFactorIsotropic, double kineticFactorAnisotropic, double uPotential, double j1Factor, double j2Factor, 
																	   AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSite = nbrSite;
  this->XMomentum = momentum;
  this->XIncrement = periodicity;
  this->MaxMomentum = this->XMomentum / this->XIncrement;
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
  
  this->EvaluateExponentialFactors();
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

ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::~ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian()
{
}
  

// evaluate all exponential factors
//   

void ParticleOnLatticeWithSpinKitaevHeisenbergAnd1DTranslationHamiltonian::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex[this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
    {
      this->ExponentialFactors[i] = Phase(2.0 * M_PI * this->XMomentum * ((double) i) / ((double) this->NbrSite));
    }
}



