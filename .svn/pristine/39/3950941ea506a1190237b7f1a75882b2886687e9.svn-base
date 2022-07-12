////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of tight binding model for the simple C4 quadrupole insulator    //
//                                                                            //
//                        last modification : 10/03/2018                      //
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
#include "Tools/FTITightBinding/TightBindingModelSimpleC4Quadrupole.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <iostream>

using std::cout;
using std::endl;
using std::ostream;



// constructor
//
// nbrSiteX = number of unit cells in the x direction
// nbrSiteY = number of unit cells in the y direction
// t1 = hoping amplitude between neareast neighbor sites within the unit cell
// phi1 = phase (in pi units) for the the hoping between neareast neighbor sites within the unit cell
// t2 = hoping amplitude between neareast neighbor sites between unit cells
// phi2 = phase (in pi units) for the the hoping between neareast neighbor sites between unit cells
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelSimpleC4Quadrupole::TightBindingModelSimpleC4Quadrupole(int nbrSiteX, int nbrSiteY, double t1, double phi1, double t2, double phi2,
									 double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->Nx1 = this->NbrSiteX;
  this->Ny1 = 0;
  this->Nx2 = 0;
  this->Ny2 = this->NbrSiteY;
  this->NNHopping1 = t1;
  this->NNHoppingPhase1 = phi1 * M_PI;
  this->NNHopping2 = t2;
  this->NNHoppingPhase2 = phi2 * M_PI;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 4;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  this->NbrConnectedOrbitals = 0 ;

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

TightBindingModelSimpleC4Quadrupole::~TightBindingModelSimpleC4Quadrupole()
{
}

// find the orbitals connected to those located at the origin unit cell
// 
  
void TightBindingModelSimpleC4Quadrupole::FindConnectedOrbitals()
{
  if (this->NbrConnectedOrbitals == 0)
    {
      Complex TmpPhaseFactor1 = this->NNHopping1 * Phase(this->NNHoppingPhase1);
      Complex TmpPhaseFactor2 = this->NNHopping2 * Phase(this->NNHoppingPhase2);
      Complex TmpConjPhaseFactor1 = this->NNHopping1 * Phase(-this->NNHoppingPhase1);
      Complex TmpConjPhaseFactor2 = this->NNHopping2 * Phase(-this->NNHoppingPhase2);

      this->NbrConnectedOrbitals = new int [this->NbrBands];
      this->ConnectedOrbitalIndices = new int* [this->NbrBands];
      this->ConnectedOrbitalSpatialIndices = new int* [this->NbrBands];
      this->ConnectedOrbitalHoppingAmplitudes = new Complex* [this->NbrBands];
      this->NbrConnectedOrbitals[0] = 4; 
      this->NbrConnectedOrbitals[1] = 4;      
      this->NbrConnectedOrbitals[2] = 4; 
      this->NbrConnectedOrbitals[3] = 4;      
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  this->ConnectedOrbitalIndices[i] = new int[this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalSpatialIndices[i] = new int[2 * this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalHoppingAmplitudes[i] = new Complex[this->NbrConnectedOrbitals[i]];
	}
      
      int TmpIndex = 0;
     
      TmpIndex = 0;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = TmpConjPhaseFactor1;
      ++TmpIndex;
       this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = this->NbrSiteX - 1;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = TmpConjPhaseFactor2;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 3;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = TmpPhaseFactor1;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 3;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = this->NbrSiteY - 1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = TmpPhaseFactor2;
      ++TmpIndex;

      TmpIndex = 0;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = TmpPhaseFactor1;
      ++TmpIndex;
       this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 1;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = TmpPhaseFactor2;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 2;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = TmpConjPhaseFactor1;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 2;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = this->NbrSiteY - 1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = TmpConjPhaseFactor2;
      ++TmpIndex;

      TmpIndex = 0;
      this->ConnectedOrbitalIndices[2][TmpIndex] = 3;
      this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = TmpConjPhaseFactor1;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[2][TmpIndex] = 3;
      this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 1;
      this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = TmpConjPhaseFactor2;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[2][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = TmpPhaseFactor1;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[2][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) +1] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = TmpPhaseFactor2;
      ++TmpIndex;

      TmpIndex = 0;
      this->ConnectedOrbitalIndices[3][TmpIndex] = 2;
      this->ConnectedOrbitalSpatialIndices[3][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[3][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[3][TmpIndex] = TmpPhaseFactor1;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[3][TmpIndex] = 2;
      this->ConnectedOrbitalSpatialIndices[3][TmpIndex * 2] =  this->NbrSiteX - 1;
      this->ConnectedOrbitalSpatialIndices[3][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[3][TmpIndex] = TmpPhaseFactor2;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[3][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[3][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[3][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[3][TmpIndex] = TmpConjPhaseFactor1;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[3][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[3][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[3][(TmpIndex * 2) +1] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[3][TmpIndex] = TmpConjPhaseFactor2;
      ++TmpIndex;
    }
}

