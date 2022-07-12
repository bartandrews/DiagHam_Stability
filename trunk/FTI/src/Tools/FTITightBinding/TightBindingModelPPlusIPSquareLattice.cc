////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the Chern insulator of         //
//                 PHYSICAL REVIEW B 84, 075128 (2011) equivalent             //
//                           to a p+ip superconductor                         //
//                                                                            //
//                        last modification : 19/02/2020                      //
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
#include "Tools/FTITightBinding/TightBindingModelPPlusIPSquareLattice.h"
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
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelPPlusIPSquareLattice::TightBindingModelPPlusIPSquareLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double mus, 
									     double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->Nx1 = this->NbrSiteX;
  this->Ny1 = 0;
  this->Nx2 = 0;
  this->Ny2 = this->NbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->SecondNextNNHoping = 0.0;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 2;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  this->BlochFormFlag = true;
  
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
  if (this->BlochFormFlag == true)
    {
      this->FindConnectedOrbitals();
      this->ComputeAllProjectedMomenta();
    }
  else
    {
      this->ComputeAllProjectedMomenta();
    }

  this->ComputeBandStructure();
}



// constructor for a tilted lattice
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nx1 = first coordinate of the first spanning vector of the tilted lattice
// nx2 = second coordinate of the first spanning vector of the tilted lattice
// ny1 = first coordinate of the second spanning vector of the tilted lattice
// ny2 = second coordinate of the second spanning vector of the tilted lattice
// offset = second coordinate in momentum space of the second spanning vector of the reciprocal lattice (0 if lattice is untilted or if Ny = 1)
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelPPlusIPSquareLattice::TightBindingModelPPlusIPSquareLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double mus, 
									   double gammaX, double gammaY, AbstractArchitecture* architecture, int offsetReal, bool storeOneBodyMatrices)
{
  
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->Nx1 = nx1;
  this->Ny1 = ny1;
  this->Nx2 = nx2;
  this->Ny2 = ny2;
  this->Offset = offset;
  this->OffsetReal = offsetReal;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->SecondNextNNHoping = 0.0;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 2;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  this->BlochFormFlag = true;
  
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
  if (this->BlochFormFlag == true)
  {
    this->FindConnectedOrbitals();
  }
  this->ComputeBandStructure();
}

// destructor
//

TightBindingModelPPlusIPSquareLattice::~TightBindingModelPPlusIPSquareLattice()
{
}

// find the orbitals connected to those located at the origin unit cell
// 
  
void TightBindingModelPPlusIPSquareLattice::FindConnectedOrbitals()
{
  double* EmbeddingX = new double [2];
  double* EmbeddingY = new double [2];
  EmbeddingX [0] = 0.0;
  EmbeddingX [1] = 0.0;
  EmbeddingY [0] = 0.0;
  EmbeddingY [1] = 0.5;
  int p;
  int q;
  if (this->NbrConnectedOrbitals == 0)
    {
      this->NbrConnectedOrbitals = new int [this->NbrBands];
      this->ConnectedOrbitalIndices = new int* [this->NbrBands];
      this->ConnectedOrbitalSpatialIndices = new int* [this->NbrBands];
      this->ConnectedOrbitalHoppingAmplitudes = new Complex* [this->NbrBands];
      this->NbrConnectedOrbitals[0] = 9; 
      this->NbrConnectedOrbitals[1] = 9;

      for (int i = 0; i < this->NbrBands; ++i)
	{
	  this->ConnectedOrbitalIndices[i] = new int[this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalSpatialIndices[i] = new int[2 * this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalHoppingAmplitudes[i] = new Complex[this->NbrConnectedOrbitals[i]];
	}
      
      int TmpIndex = 0;
      this->GetRealSpaceIndex(0, 0, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->MuS;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = -this->MuS;
      ++TmpIndex;
      
      this->GetRealSpaceIndex(0, 0, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NNHoping * Phase(-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase(-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
      ++TmpIndex;
      this->GetRealSpaceIndex(1, 0, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = -this->NNHoping * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2])) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1])));
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
     this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2])) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1])));
      ++TmpIndex;
      this->GetRealSpaceIndex(-1, 0, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = -this->NNHoping * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2])) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1])));
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
       this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2])) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1])));
      ++TmpIndex;
      this->GetRealSpaceIndex(0, -1, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NNHoping * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
      this->GetRealSpaceIndex(0, 1, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
      ++TmpIndex;

      
      this->GetRealSpaceIndex(1, 0, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NextNNHoping * Phase (M_PI * 0.5) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
      this->GetRealSpaceIndex(1, 1, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NextNNHoping * Phase (-M_PI * 0.5) * Phase (-this->KxFactor * this->GammaX *  ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
      ++TmpIndex;
      this->GetRealSpaceIndex(-1, -1, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NextNNHoping * Phase (M_PI * 0.5) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
      this->GetRealSpaceIndex(-1, 0, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NextNNHoping * Phase (-M_PI * 0.5) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
      ++TmpIndex;
      this->GetRealSpaceIndex(-1, 0, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NextNNHoping * Phase (-M_PI * 0.5) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
      this->GetRealSpaceIndex(-1, 1, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NextNNHoping * Phase (M_PI * 0.5) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
      ++TmpIndex;
      this->GetRealSpaceIndex(1, -1, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NextNNHoping * Phase (-M_PI * 0.5) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
      this->GetRealSpaceIndex(1, 0, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NextNNHoping * Phase (M_PI * 0.5) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
      ++TmpIndex;
      
    }
  delete[] EmbeddingX;
  delete[] EmbeddingY;
}

