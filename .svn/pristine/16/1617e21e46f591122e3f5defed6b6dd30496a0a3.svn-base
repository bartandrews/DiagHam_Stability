////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of tight binding model for the Kitaev-Heisenberg honeycomb lattice  //
//                                                                            //
//                        last modification : 01/08/2015                      //
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
#include "Tools/FTITightBinding/TightBindingModelKitaevHeisenbergHoneycombLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <iostream>

using std::cout;
using std::endl;
using std::ostream;



// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t = isotropic hoping amplitude between neareast neighbor sites
// tPrime = anisotropic hoping amplitude between neareast neighbor sites
// tPrimePhase = additional phase for the anisotropic hoping amplitude between neareast neighbor sitesc(in pi units)
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelKitaevHeisenbergHoneycombLattice::TightBindingModelKitaevHeisenbergHoneycombLattice(int nbrSiteX, int nbrSiteY, double t, double tPrime, double tPrimePhase,
												     AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->Nx1 = this->NbrSiteX;
  this->Ny1 = 0;
  this->Nx2 = 0;
  this->Ny2 = this->NbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->KineticFactorIsotropic = t;
  this->KineticFactorAnisotropic = tPrime;
  this->KineticFactorAnisotropicPhase = tPrimePhase;
  this->HFieldX = 0.0;
  this->HFieldY = 0.0;
  this->HFieldZ = 0.0;
  this->GammaX = 0.0;
  this->GammaY = 0.0;
  this->NbrBands = 4;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  
  this->ComputeAllProjectedMomenta();
  
  if (storeOneBodyMatrices == true)
    {
      this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
    }
  else
    {
      this->OneBodyBasis = 0;
    }
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    }
  this->FindConnectedOrbitals();
  this->ComputeBandStructure();
}


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t = isotropic hoping amplitude between neareast neighbor sites
// tPrime = anisotropic hoping amplitude between neareast neighbor sites
// tPrimePhase = additional phase for the anisotropic hoping amplitude between neareast neighbor sitesc(in pi units)
// hX = amplitude of the magnetic field along the x direction
// hY = amplitude of the magnetic field along the y direction
// hZ = amplitude of the magnetic field along the z direction
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelKitaevHeisenbergHoneycombLattice::TightBindingModelKitaevHeisenbergHoneycombLattice(int nbrSiteX, int nbrSiteY, double t, double tPrime, double tPrimePhase, double hX, double hY, double hZ,
												     AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->Nx1 = this->NbrSiteX;
  this->Ny1 = 0;
  this->Nx2 = 0;
  this->Ny2 = this->NbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->KineticFactorIsotropic = t;
  this->KineticFactorAnisotropic = tPrime;
  this->KineticFactorAnisotropicPhase = tPrimePhase;
  this->HFieldX = hX;
  this->HFieldY = hY;
  this->HFieldZ = hZ;
  this->GammaX = 0.0;
  this->GammaY = 0.0;
  this->NbrBands = 4;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  
  this->ComputeAllProjectedMomenta();
  
  if (storeOneBodyMatrices == true)
    {
      this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
    }
  else
    {
      this->OneBodyBasis = 0;
    }
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    }
  this->FindConnectedOrbitals();
  this->ComputeBandStructure();
}


// destructor
//

TightBindingModelKitaevHeisenbergHoneycombLattice::~TightBindingModelKitaevHeisenbergHoneycombLattice()
{
}

// find the orbitals connected to those located at the origin unit cell
// 
  
void TightBindingModelKitaevHeisenbergHoneycombLattice::FindConnectedOrbitals()
{
  if (this->NbrConnectedOrbitals == 0)
    {
      this->NbrConnectedOrbitals = new int [this->NbrBands];
      this->ConnectedOrbitalIndices = new int* [this->NbrBands];
      this->ConnectedOrbitalSpatialIndices = new int* [this->NbrBands];
      this->ConnectedOrbitalHoppingAmplitudes = new Complex* [this->NbrBands];
      for (int i = 0; i < 4; ++i)
	this->NbrConnectedOrbitals[i] = 5;
      if ((this->HFieldX != 0.0) || (this->HFieldY != 0.0))
	for (int i = 0; i < 4; ++i)
	  this->NbrConnectedOrbitals[i] += 1;
      if (this->HFieldZ != 0.0)
	for (int i = 0; i < 4; ++i)
	  this->NbrConnectedOrbitals[i] += 1;
	
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  this->ConnectedOrbitalIndices[i] = new int[this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalSpatialIndices[i] = new int[2 * this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalHoppingAmplitudes[i] = new Complex[this->NbrConnectedOrbitals[i]];
	}
      Complex TmpAnisotropicFactor = this->KineticFactorAnisotropic * Phase(this->KineticFactorAnisotropicPhase * M_PI);
       int TmpIndex = 0;
       // site A, up spin
       this->ConnectedOrbitalIndices[0][TmpIndex] = 2;
       this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
       this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->KineticFactorIsotropic;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[0][TmpIndex] = 2;
       this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = -1;
       this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->KineticFactorIsotropic;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[0][TmpIndex] = 2;
       this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = -1;
       this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 1;
       this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->KineticFactorIsotropic + TmpAnisotropicFactor;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[0][TmpIndex] = 3;
       this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
       this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = I() * TmpAnisotropicFactor;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[0][TmpIndex] = 3;
       this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = -1;
       this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = TmpAnisotropicFactor;
       ++TmpIndex;
       if ((this->HFieldX != 0.0) || (this->HFieldY != 0.0))
       {
	this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
	this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
	this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
	this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = 0.5 * Complex(this->HFieldX, this->HFieldY);
	++TmpIndex;
       }
       if (this->HFieldZ != 0.0)
       {
	this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
	this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
	this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
	this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = 0.5 * this->HFieldZ;
	++TmpIndex;
       }

       // site A, down spin
       TmpIndex = 0;
       this->ConnectedOrbitalIndices[1][TmpIndex] = 3;
       this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
       this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->KineticFactorIsotropic;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[1][TmpIndex] = 3;
       this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = -1;
       this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->KineticFactorIsotropic;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[1][TmpIndex] = 3;
       this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = -1;
       this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 1;
       this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->KineticFactorIsotropic - TmpAnisotropicFactor;;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[1][TmpIndex] = 2;
       this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
       this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = -I() * TmpAnisotropicFactor;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[1][TmpIndex] = 2;
       this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = -1;
       this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = TmpAnisotropicFactor;
       ++TmpIndex;
       if ((this->HFieldX != 0.0) || (this->HFieldY != 0.0))
       {
	this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
	this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
	this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
	this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = 0.5 * Complex(this->HFieldX, -this->HFieldY);
	++TmpIndex;
       }
       if (this->HFieldZ != 0.0)
       {
	this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
	this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
	this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
	this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = -0.5 * this->HFieldZ;
	++TmpIndex;
       }

       // site B, up spin
       TmpIndex = 0;
       this->ConnectedOrbitalIndices[2][TmpIndex] = 0;
       this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 0;
       this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = this->KineticFactorIsotropic;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[2][TmpIndex] = 0;
       this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 1;
       this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = this->KineticFactorIsotropic;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[2][TmpIndex] = 0;
       this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 1;
       this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = -1;
       this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = this->KineticFactorIsotropic + Conj(TmpAnisotropicFactor);
       ++TmpIndex;
       this->ConnectedOrbitalIndices[2][TmpIndex] = 1;
       this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 0;
       this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = Conj(-I() * TmpAnisotropicFactor);
       ++TmpIndex;
       this->ConnectedOrbitalIndices[2][TmpIndex] = 1;
       this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 1;
       this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = Conj(TmpAnisotropicFactor);
       ++TmpIndex;
       if ((this->HFieldX != 0.0) || (this->HFieldY != 0.0))
       {
	this->ConnectedOrbitalIndices[2][TmpIndex] = 3;
	this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 0;
	this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
	this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = 0.5 * Complex(this->HFieldX, this->HFieldY);
	++TmpIndex;
       }
       if (this->HFieldZ != 0.0)
       {
	this->ConnectedOrbitalIndices[2][TmpIndex] = 2;
	this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 0;
	this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
	this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = 0.5 * this->HFieldZ;
	++TmpIndex;
       }

       // site B, down spin
       TmpIndex = 0;
       this->ConnectedOrbitalIndices[3][TmpIndex] = 1;
       this->ConnectedOrbitalSpatialIndices[3][TmpIndex * 2] = 0;
       this->ConnectedOrbitalSpatialIndices[3][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[3][TmpIndex] = this->KineticFactorIsotropic;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[3][TmpIndex] = 1;
       this->ConnectedOrbitalSpatialIndices[3][TmpIndex * 2] = 1;
       this->ConnectedOrbitalSpatialIndices[3][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[3][TmpIndex] = this->KineticFactorIsotropic;
       ++TmpIndex;
       this->ConnectedOrbitalIndices[3][TmpIndex] = 1;
       this->ConnectedOrbitalSpatialIndices[3][TmpIndex * 2] = 1;
       this->ConnectedOrbitalSpatialIndices[3][(TmpIndex * 2) + 1] = -1;
       this->ConnectedOrbitalHoppingAmplitudes[3][TmpIndex] = this->KineticFactorIsotropic - Conj(TmpAnisotropicFactor);
       ++TmpIndex;
       this->ConnectedOrbitalIndices[3][TmpIndex] = 0;
       this->ConnectedOrbitalSpatialIndices[3][TmpIndex * 2] = 0;
       this->ConnectedOrbitalSpatialIndices[3][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[3][TmpIndex] = Conj(I() * TmpAnisotropicFactor);
       ++TmpIndex;
       this->ConnectedOrbitalIndices[3][TmpIndex] = 0;
       this->ConnectedOrbitalSpatialIndices[3][TmpIndex * 2] = 1;
       this->ConnectedOrbitalSpatialIndices[3][(TmpIndex * 2) + 1] = 0;
       this->ConnectedOrbitalHoppingAmplitudes[3][TmpIndex] = Conj(TmpAnisotropicFactor);
       ++TmpIndex;
       if ((this->HFieldX != 0.0) || (this->HFieldY != 0.0))
       {
	this->ConnectedOrbitalIndices[3][TmpIndex] = 2;
	this->ConnectedOrbitalSpatialIndices[3][TmpIndex * 2] = 0;
	this->ConnectedOrbitalSpatialIndices[3][(TmpIndex * 2) + 1] = 0;
	this->ConnectedOrbitalHoppingAmplitudes[3][TmpIndex] = 0.5 * Complex(this->HFieldX, -this->HFieldY);
	++TmpIndex;
       }
       if (this->HFieldZ != 0.0)
       {
	this->ConnectedOrbitalIndices[3][TmpIndex] = 3;
	this->ConnectedOrbitalSpatialIndices[3][TmpIndex * 2] = 0;
	this->ConnectedOrbitalSpatialIndices[3][(TmpIndex * 2) + 1] = 0;
	this->ConnectedOrbitalHoppingAmplitudes[3][TmpIndex] = -0.5 * this->HFieldZ;
	++TmpIndex;
       }

    }
}

