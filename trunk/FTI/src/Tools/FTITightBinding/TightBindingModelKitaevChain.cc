////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//              class of tight binding model for the Kitaev chain             //
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
#include "Tools/FTITightBinding/TightBindingModelKitaevChain.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "MathTools/BinomialCoefficients.h"
#include "Matrix/HermitianMatrix.h"

#include <iostream>


using std::cout;
using std::endl;


// default constructor
//
// nbrSiteX = number of sites
// hopping = hopping between nearest neighboring sites
// delta = superconducting order parameter 
// mus = on site chemical potential
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelKitaevChain::TightBindingModelKitaevChain(int nbrSiteX, double hopping, double delta, double mus, 
							   AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
    this->NbrSiteX = nbrSiteX;
    this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
    this->NNHopping = hopping;
    this->Delta = delta;
    this->MuS = mus;
    this->GammaX = 0.0;
    this->NbrBands = 2;
    this->NbrStatePerBand = nbrSiteX;

    this->EmbeddingX = RealVector(this->NbrBands, true);

    this->Architecture = architecture;

    if (storeOneBodyMatrices == true)
      this->OneBodyBasis = new ComplexMatrix[this->NbrStatePerBand];
    else
      this->OneBodyBasis = 0;
    this->EnergyBandStructure = new double*[this->NbrBands];
    for (int i = 0; i < this->NbrBands; ++i)
      this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    this->ComputeBandStructure();
}

// destructor
//

TightBindingModelKitaevChain::~TightBindingModelKitaevChain()
{
}

// find the orbitals connected to those located at the origin unit cell
// 
  
void TightBindingModelKitaevChain::FindConnectedOrbitals()
{
  if (this->NbrConnectedOrbitals == 0)
    {
      this->NbrConnectedOrbitals = new int [this->NbrBands];
      this->ConnectedOrbitalIndices = new int* [this->NbrBands];
      this->ConnectedOrbitalSpatialIndices = new int* [this->NbrBands];
      this->ConnectedOrbitalHoppingAmplitudes = new Complex* [this->NbrBands];
      if (this->MuS != 0.0)
	{
	  this->NbrConnectedOrbitals[0] = 5; 
	  this->NbrConnectedOrbitals[1] = 5;      
	} 
      else
	{
	  this->NbrConnectedOrbitals[0] = 4; 
	  this->NbrConnectedOrbitals[1] = 4;
	}
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  this->ConnectedOrbitalIndices[i] = new int[this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalSpatialIndices[i] = new int[this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalHoppingAmplitudes[i] = new Complex[this->NbrConnectedOrbitals[i]];
	}
      
      int TmpIndex = 0;
      if (this->MuS != 0.0)
	{
	  this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
	  this->ConnectedOrbitalSpatialIndices[0][TmpIndex] = 0;
	  this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = 0.5 * this->MuS;
	  this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
	  this->ConnectedOrbitalSpatialIndices[1][TmpIndex] = 0;
	  this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = 0.5 * this->MuS;
	  ++TmpIndex;
	}
      
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NNHopping;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = -this->NNHopping;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NNHopping;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = -this->NNHopping;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = 0.5 * this->Delta;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = -0.5 * this->Delta;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = 0.5 * this->Delta;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = -0.5 * this->Delta;
      ++TmpIndex;
    }
}


