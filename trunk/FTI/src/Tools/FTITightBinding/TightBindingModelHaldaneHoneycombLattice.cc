////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of tight binding model for the Haldane honeycomb lattice       //
//                                                                            //
//                        last modification : 30/05/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelHaldaneHoneycombLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// phi =  Haldane phase on next neareast neighbor hopping
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelHaldaneHoneycombLattice::TightBindingModelHaldaneHoneycombLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double phi, double mus, 
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
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->HaldanePhase = phi;
  this->MuS = mus;
  this->TwistAngle = M_PI / 3;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 2;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->ComputeAllProjectedMomenta();
  this->EmbeddingX = RealVector(this->NbrBands, true);
  this->EmbeddingY = RealVector(this->NbrBands, true);

  this->EmbeddingX[0] = 0; 
  this->EmbeddingX[1] = 1.0 / 2.0;
  this->EmbeddingY[0] = 0; 
  this->EmbeddingY[1] = 1.0 / 2.0;

  this->EmbeddingX[1] = 0.0;
  this->EmbeddingY[1] = 0.0;
  
/*  this->EmbeddingX[0] = 1.0 / 6;
  this->EmbeddingX[1] = -1.0 / 6;
  this->EmbeddingY[0] = 1.0 / 3;
  this->EmbeddingY[1] = - 1.0 / 3;
*/
  this->Inversion = ComplexMatrix(this->NbrBands, this->NbrBands, true);
  for (int i = 0; i < this->NbrBands; ++i)
      this->Inversion[i][1 - i] = 1.0;

  this->Architecture = architecture;

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
  this->ComputeBandStructure();
}


// constructor for the tilted lattice
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
//nx1 = first coordinate of the first spanning vector for a tilted lattice
//ny1 = second coordinate of the first spanning vector for a tilted lattice
//nx2 = first coordinate of the second spanning vector for a tilted lattice
//ny2 = second coordinate of the second spanning vector for a tilted lattice
//offset = second coordinate in momentum space of the second spanning vector of the reciprocal lattice for a tilted lattice
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// phi =  Haldane phase on next neareast neighbor hopping
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelHaldaneHoneycombLattice::TightBindingModelHaldaneHoneycombLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double phi, double mus, 
										   double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->Nx1 = nx1;
  this->Ny1 = ny1;
  this->Nx2 = nx2;
  this->Ny2 = ny2;
  this->Offset = offset;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->HaldanePhase = phi;
  this->MuS = mus;
  this->TwistAngle = M_PI / 3;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 2;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;

  this->EmbeddingX = RealVector(this->NbrBands, true);
  this->EmbeddingX[0] = 1.0 / 6;
  this->EmbeddingX[1] = -1.0 / 6;
  this->EmbeddingY = RealVector(this->NbrBands, true);
  this->EmbeddingY[0] = 1.0 / 3;
  this->EmbeddingY[1] = - 1.0 / 3;
  this->Inversion = ComplexMatrix(this->NbrBands, this->NbrBands, true);
  for (int i = 0; i < this->NbrBands; ++i)
      this->Inversion[i][1 - i] = 1.0;

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
  this->ComputeBandStructure();
}


// destructor
//

TightBindingModelHaldaneHoneycombLattice::~TightBindingModelHaldaneHoneycombLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelHaldaneHoneycombLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
	      double x = this->GetProjectedMomentum(kx, ky, 0);
              double y = this->GetProjectedMomentum(kx, ky, 1);

	      Complex B1 = this->NNHopping * Complex(1 + cos(x+y) + cos(y), + sin(x+y) + sin(y));
	      double d0 = 2.0 * this->NextNNHopping * cos(this->HaldanePhase) * (cos(x) + cos(y) + cos(x+y));
	      double d3 = 2.0 * this->NextNNHopping * sin(this->HaldanePhase) * (sin(x) + sin(y) - sin(x+y)) + this->MuS;

	      HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d0 + d3);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 1, d0 - d3);
	      
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
		    this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		}
	    }
	}
    }
}

// find the orbitals connected to those located at the origin unit cell
// 
  
void TightBindingModelHaldaneHoneycombLattice::FindConnectedOrbitals()
{
  if (this->NbrConnectedOrbitals == 0)
  {
    this->NbrConnectedOrbitals = new int [this->NbrBands];
    this->ConnectedOrbitalIndices = new int* [this->NbrBands];
    this->ConnectedOrbitalSpatialIndices = new int* [this->NbrBands];
    this->ConnectedOrbitalHoppingAmplitudes = new Complex* [this->NbrBands];
  
    this->NbrConnectedOrbitals[0] = 3; 
    this->NbrConnectedOrbitals[1] = 3;      

    if ( this->NextNNHopping != 0.0 ) 
      {       
	this->NbrConnectedOrbitals[0] += 6; 
	this->NbrConnectedOrbitals[1] += 6;      

      }
    if (this->MuS != 0.0)
      {
	++this->NbrConnectedOrbitals[0];
      }
    for (int i = 0; i < this->NbrBands; ++i)
      {
	this->ConnectedOrbitalIndices[i] = new int[this->NbrConnectedOrbitals[i]];
	this->ConnectedOrbitalSpatialIndices[i] = new int[2 * this->NbrConnectedOrbitals[i]];
	this->ConnectedOrbitalHoppingAmplitudes[i] = new Complex[this->NbrConnectedOrbitals[i]];
      }

  
    Complex Lambda1 = this->NextNNHopping * Phase(this->HaldanePhase); 
    
    int TmpIndex = 0;
    
    // links starting from A
    this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
    this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->NNHopping * Phase(-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
    ++TmpIndex;

    this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
    this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = -1;
    this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] =  this->NNHopping * Phase(-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0]));
    ++TmpIndex;
    this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
    this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = -1;
    this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] =  this->NNHopping * Phase(-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2]) + EmbeddingX[1] - EmbeddingX[0]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]) + EmbeddingY[1] - EmbeddingY[0])); 
    ++TmpIndex;


    if (this->NextNNHopping != 0.0)
    {
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 1;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Lambda1 * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Conj(Lambda1) *  Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = -1;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Lambda1 *  Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = -1;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Conj(Lambda1) *  Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Lambda1 *  Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 1;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Conj(Lambda1) *  Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1]);
      ++TmpIndex;
    }
    
    if (this->MuS != 0.0)
      {
	this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
	this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
	this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = 0;
	this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->MuS;
      }

    TmpIndex = 0;

    // links starting from B
    this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
    this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NNHopping * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
    ++TmpIndex;
    this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
    this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 1;
    this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NNHopping * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));;
    ++TmpIndex;
    this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
    this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 1;
    this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = this->NNHopping * Phase (-this->KxFactor * this->GammaX * ((double) (this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2]) + EmbeddingX[0] - EmbeddingX[1]) - this->KyFactor * this->GammaY * ((double) (this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]) + EmbeddingY[0] - EmbeddingY[1]));
    ++TmpIndex;
    
    if (this->NextNNHopping != 0.0)
    {
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 1;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Conj(Lambda1) * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]);
      ++TmpIndex; 
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Lambda1 * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = -1;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Conj(Lambda1) * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = -1;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Lambda1 * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Conj(Lambda1) * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 1;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Lambda1 * Phase (-this->KxFactor * this->GammaX * this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] - this->KyFactor * this->GammaY * this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) +1]);
      ++TmpIndex;
    }
  }
}
