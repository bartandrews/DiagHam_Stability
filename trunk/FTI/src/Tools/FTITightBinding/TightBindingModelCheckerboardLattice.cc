////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of tight binding model for the checkerboard lattice        //
//                                                                            //
//                        last modification : 08/05/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelCheckerboardLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <iostream>

using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

TightBindingModelCheckerboardLattice::TightBindingModelCheckerboardLattice()
{
  this->NbrSiteX = 1;
  this->NbrSiteY = 1;
  this->Nx1 = this->NbrSiteX;
  this->Ny1 = 0;
  this->Nx2 = 0;
  this->Ny2 = this->NbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHoping = 0.0;
  this->NextNNHoping = 0.0;
  this->SecondNextNNHoping = 0.0;
  this->MuS = 0.0;
  this->GammaX = 0.0;
  this->GammaY = 0.0;
  this->NbrBands = 2;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = 0;
  this->BlochFormFlag = true;
}

// constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// t2p = hoping amplitude between second next neareast neighbor sites
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
// blochFormFlag = use the Bloch form instead of the the traditional form

TightBindingModelCheckerboardLattice::TightBindingModelCheckerboardLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double t2p, double mus, 
									   double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices, bool blochFormFlag)
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
  this->SecondNextNNHoping = t2p;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 2;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  this->BlochFormFlag = blochFormFlag;
  
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
    this->FindConnectedOrbitals();
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
// t2p = hoping amplitude between second next neareast neighbor sites
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
// blochFormFlag = use the Bloch form instead of the the traditional form

TightBindingModelCheckerboardLattice::TightBindingModelCheckerboardLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double t2p, double mus, 
									   double gammaX, double gammaY, AbstractArchitecture* architecture, int offsetReal, bool storeOneBodyMatrices, bool blochFormFlag)
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
  this->SecondNextNNHoping = t2p;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 2;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  this->BlochFormFlag = blochFormFlag;
  
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

TightBindingModelCheckerboardLattice::~TightBindingModelCheckerboardLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelCheckerboardLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (this->BlochFormFlag == false)
    {
      if (nbrStates == 0l)
	nbrStates = this->NbrStatePerBand;
      long MaxStateIndex = minStateIndex + nbrStates;
      double KX;
      double KY;
      for (int kx = 0; kx < this->NbrSiteX; ++kx)
	{
	  for (int ky = 0; ky < this->NbrSiteY; ++ky)
	    {
	      int Index = this->GetLinearizedMomentumIndex(kx, ky);
	      if ((Index >= minStateIndex) && (Index < MaxStateIndex))
		{
		  Complex B1 = 4.0 * this->NNHoping * Complex (cos (0.5 * this->ProjectedMomenta[Index][0]) * cos (0.5 * this->ProjectedMomenta[Index][1]) * cos(M_PI * 0.25), 
							       sin (0.5 * this->ProjectedMomenta[Index][0]) * sin (0.5 * this->ProjectedMomenta[Index][1]) * sin(M_PI * 0.25));
		  double d1 = 4.0 * this->SecondNextNNHoping * cos (this->ProjectedMomenta[Index][0]) * cos (this->ProjectedMomenta[Index][1]);
		  double d3 =  this->MuS + (2.0 * this->NextNNHoping * (cos (this->ProjectedMomenta[Index][0])
									- cos (this->ProjectedMomenta[Index][1])));
		  
		  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
		  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
		  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
		  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
		  if (this->OneBodyBasis != 0)
		    {
		      ComplexMatrix TmpMatrix(this->NbrBands, this->NbrBands, true);
		      TmpMatrix.SetToIdentity();
		      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		      TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
		      TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
		      this->OneBodyBasis[Index] = TmpMatrix;

		      for (int i = 0; i < this->NbrBands; ++i)
			this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		    }
		  else
		    {
		      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		      TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
		      TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif
		      for (int i = 0; i < this->NbrBands; ++i)
		      {
			this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
// 			cout << i << " " << Index << " " << TmpDiag(i,i) << endl;
		      }
		    }
		}
	    }
	}
    }
  else
    {
      this->Abstract2DTightBindingModel::CoreComputeBandStructure(minStateIndex, nbrStates);
    }
}

// find the orbitals connected to those located at the origin unit cell
// 
  
void TightBindingModelCheckerboardLattice::FindConnectedOrbitals()
{
  double* EmbeddingX = new double [2];
  double* EmbeddingY = new double [2];
  EmbeddingX [0] = 0.0;
  EmbeddingX [1] = 0.5;
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
      if (this->SecondNextNNHoping != 0.0)
	{
	  this->NbrConnectedOrbitals[0] = 13; 
	  this->NbrConnectedOrbitals[1] = 13;      
	} 
      else
	{
	  this->NbrConnectedOrbitals[0] = 9; 
	  this->NbrConnectedOrbitals[1] = 9;
	}
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  this->ConnectedOrbitalIndices[i] = new int[this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalSpatialIndices[i] = new int[2 * this->NbrConnectedOrbitals[i]];
	  this->ConnectedOrbitalHoppingAmplitudes[i] = new Complex[this->NbrConnectedOrbitals[i]];
	}
      
      int TmpIndex = 0;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->MuS;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = -this->MuS;
      ++TmpIndex;
      
      this->GetRealSpaceIndex(0, 0, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NNHoping * Phase (M_PI * 0.25) * Phase(-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase (-M_PI * 0.25) * Phase(-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
      ++TmpIndex;
      this->GetRealSpaceIndex(-1, -1, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NNHoping * Phase (M_PI * 0.25) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
      this->GetRealSpaceIndex(1, 1, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase (-M_PI * 0.25) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
      ++TmpIndex;
      this->GetRealSpaceIndex(-1, 0, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NNHoping * Phase (-M_PI * 0.25) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
      this->GetRealSpaceIndex(0, 1, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase (M_PI * 0.25) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
            
      ++TmpIndex;
      this->GetRealSpaceIndex(0, -1, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NNHoping * Phase (-M_PI * 0.25) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
      this->GetRealSpaceIndex(1, 0, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NNHoping * Phase (M_PI * 0.25) * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
      ++TmpIndex;
      
      this->GetRealSpaceIndex(1, 0, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
      this->GetRealSpaceIndex(1, 0, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = -this->NextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->GetRealSpaceIndex(-1, 0, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
      this->GetRealSpaceIndex(-1, 0, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = -this->NextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->GetRealSpaceIndex(0, 1, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = -this->NextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
      this->GetRealSpaceIndex(0, 1, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->GetRealSpaceIndex(0, -1, p , q);
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = -this->NextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
      this->GetRealSpaceIndex(0, -1, p , q);
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = p;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1] = q;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]);
      ++TmpIndex;
      
      if (this->SecondNextNNHoping != 0.0)
	{
	  for (int i = 0; i < this->NbrBands; ++i)
	    {
	      this->GetRealSpaceIndex(1, 1, p , q);
	      this->ConnectedOrbitalIndices[i][TmpIndex] = i;
	      this->ConnectedOrbitalSpatialIndices[i][TmpIndex * 2] = p;
	      this->ConnectedOrbitalSpatialIndices[i][(TmpIndex * 2) + 1] = q;
	      this->ConnectedOrbitalHoppingAmplitudes[i][TmpIndex] = this->SecondNextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
	    }
	  ++TmpIndex;
	  for (int i = 0; i < this->NbrBands; ++i)
	    {
	      this->GetRealSpaceIndex(-1, 1, p , q);
	      this->ConnectedOrbitalIndices[i][TmpIndex] = i;
	      this->ConnectedOrbitalSpatialIndices[i][TmpIndex * 2] = p;
	      this->ConnectedOrbitalSpatialIndices[i][(TmpIndex * 2) +1] = q;
	      this->ConnectedOrbitalHoppingAmplitudes[i][TmpIndex] = this->SecondNextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
	    }
	  ++TmpIndex;
	  for (int i = 0; i < this->NbrBands; ++i)
	    {
	      this->GetRealSpaceIndex(1, -1, p , q);
	      this->ConnectedOrbitalIndices[i][TmpIndex] = i;
	      this->ConnectedOrbitalSpatialIndices[i][TmpIndex * 2] = p;
	      this->ConnectedOrbitalSpatialIndices[i][(TmpIndex * 2) +1] = q;
	      this->ConnectedOrbitalHoppingAmplitudes[i][TmpIndex] = this->SecondNextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
	    }
	  ++TmpIndex;
	  for (int i = 0; i < this->NbrBands; ++i)
	    {
	      this->GetRealSpaceIndex(-1, -1, p , q);
	      this->ConnectedOrbitalIndices[i][TmpIndex] = i;
	      this->ConnectedOrbitalSpatialIndices[i][TmpIndex * 2] = p;
	      this->ConnectedOrbitalSpatialIndices[i][(TmpIndex * 2) +1] = q;
	      this->ConnectedOrbitalHoppingAmplitudes[i][TmpIndex] = this->SecondNextNNHoping * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
	    }
	  ++TmpIndex;
	}
    }
   delete[] EmbeddingX;
   delete[] EmbeddingY;
}

