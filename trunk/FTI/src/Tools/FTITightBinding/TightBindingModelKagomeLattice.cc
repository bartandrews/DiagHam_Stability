////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of tight binding model for the kagome lattice          //
//                                                                            //
//                        last modification : 07/05/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelKagomeLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A1 sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelKagomeLattice::TightBindingModelKagomeLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double lambda1, double lambda2, double mus, 
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
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->NNSpinOrbit = lambda1;
  this->NextNNSpinOrbit = lambda2;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 3;
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
// ny1 = second coordinate of the first spanning vector of the tilted lattice
// nx2 = first coordinate of the second spanning vector of the tilted lattice
// ny2 = second coordinate of the second spanning vector of the tilted lattice
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A1 sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelKagomeLattice::TightBindingModelKagomeLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double lambda1, double lambda2, double mus, 
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
  this->NNSpinOrbit = lambda1;
  this->NextNNSpinOrbit = lambda2;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 3;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  this->BlochFormFlag = false;
  
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

TightBindingModelKagomeLattice::~TightBindingModelKagomeLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelKagomeLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
	      KX = this->ProjectedMomenta[Index] [0];
	      KY = this->ProjectedMomenta[Index] [1];

	      HermitianMatrix TmpOneBodyHamiltonian = this->ComputeBlochHamiltonian(KX, KY);

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
  else
  {
      this->Abstract2DTightBindingModel::CoreComputeBandStructure(minStateIndex, nbrStates);
  }
}

// compute the Bloch hamiltonian at a point of the Brillouin zone
//
// kx = momentum along the x axis
// ky = momentum along the x axis
// return value = Bloch hamiltonian

HermitianMatrix TightBindingModelKagomeLattice::ComputeBlochHamiltonian(double kx, double ky)
{
  Complex HAB (-2.0 * this->NNHopping, -2.0 * this->NNSpinOrbit);
  HAB *= cos (kx * 0.5);
  Complex HAC(-2.0 * this->NNHopping, 2.0 * this->NNSpinOrbit);
  HAC *= cos (ky * 0.5);
  Complex HBC(-2.0 * this->NNHopping, -2.0 * this->NNSpinOrbit);
  HBC *= cos((kx - ky) * 0.5);
  
  Complex HAB2 (-2.0 * this->NextNNHopping, 2.0 * this->NextNNSpinOrbit);
  HAB2 *= cos ((kx - 2.0 * ky) * 0.5);
  Complex HAC2 (-2.0 * this->NextNNHopping, -2.0 * this->NextNNSpinOrbit);
  HAC2 *= cos ((2.0 * kx - ky) * 0.5);
  Complex HBC2 (-2.0 * this->NextNNHopping, 2.0  *  this->NextNNSpinOrbit);
  HBC2 *= cos ((kx + ky) * 0.5);
  
  HAB += HAB2;
  HAC += HAC2;
  HBC += HBC2;
  
  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, this->MuS);
  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, HAB);
  TmpOneBodyHamiltonian.SetMatrixElement(0, 2, HAC);
  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, HBC);
  return TmpOneBodyHamiltonian;
}

// get the high symmetry points 
//
// pointNames = name of each high symmetry point
// pointCoordinates = coordinates in the first Brillouin zone of the high symmetry points
// return value = number of high symmetry points

int TightBindingModelKagomeLattice::GetHighSymmetryPoints(char**& pointNames, double**& pointCoordinates)
{
  int NbrHighSymmetryPoints = 3;
  pointNames = new char*[NbrHighSymmetryPoints];
  pointCoordinates = new double*[NbrHighSymmetryPoints];

  pointNames[0] = new char[16];
  sprintf (pointNames[0], "Gamma");
  pointCoordinates[0] = new double[2];
  pointCoordinates[0][0] = 0.0; 
  pointCoordinates[0][1] = 0.0; 

  pointNames[1] = new char[16];
  sprintf (pointNames[1], "K");
  pointCoordinates[1] = new double[2];
  pointCoordinates[1][0] = 4.0 * M_PI / 3.0; 
  pointCoordinates[1][1] = 2.0 * M_PI / 3.0; 

  pointNames[2] = new char[16];
  sprintf (pointNames[2], "M");
  pointCoordinates[2] = new double[2];
  pointCoordinates[2][0] = M_PI; 
  pointCoordinates[2][1] = 0.0; 

  return NbrHighSymmetryPoints;
}

// compute the distance between two points in the first Brillouin zone, changing the coordinates the second one by a reciprocal lattice vector if needed
//
// kx1 = momentum of the first point along the x axis
// ky1 = momentum of the first point along the y axis
// kx2 = reference on the momentum of the second point along the x axis
// ky2 = reference on the momentum of the second point along the y axis
// return value = distance between the two points

double TightBindingModelKagomeLattice::GetDistanceReciprocalSpace(double kx1, double ky1, double& kx2, double& ky2)
{
  double AngleFactor = 2.0 * cos(2.0 * M_PI / 3.0);
  double DiffKx = kx1 - kx2;
  double DiffKy = ky1 - ky2;
  double MinDistance = sqrt ((DiffKx * DiffKx) + (DiffKy * DiffKy) + (AngleFactor * DiffKx * DiffKy));
  double MinKx2 = kx2;
  double MinKy2 = ky2;
  for (int i = -1; i <= 1; ++i)
    {
      double TmpKx2  = kx2 + (2.0 * ((double) i) * M_PI);
      for (int j = -1; j <= 1; ++j)
	{
	  double TmpKy2  = ky2 + (2.0 * ((double) j) * M_PI);	  
	  double DiffKx = kx1 - TmpKx2;
	  double DiffKy = ky1 - TmpKy2;
	  double TmpDistance =  sqrt ((DiffKx * DiffKx) + (DiffKy * DiffKy) + (AngleFactor * DiffKx * DiffKy));
	  if (TmpDistance < MinDistance)
	    {
	      MinDistance = TmpDistance;
	      MinKx2 = TmpKx2;
	      MinKy2 = TmpKy2;
	    }
	}
    }
  kx2 = MinKx2;
  ky2 = MinKy2;
  return MinDistance;
}


// find the orbitals connected to those located at the origin unit cell
// 
  
void TightBindingModelKagomeLattice::FindConnectedOrbitals()
{
  if (this->NbrConnectedOrbitals == 0)
  {
    this->NbrConnectedOrbitals = new int [this->NbrBands];
    this->ConnectedOrbitalIndices = new int* [this->NbrBands];
    this->ConnectedOrbitalSpatialIndices = new int* [this->NbrBands];
    this->ConnectedOrbitalHoppingAmplitudes = new Complex* [this->NbrBands];
  
    this->NbrConnectedOrbitals[0] = 4; 
    this->NbrConnectedOrbitals[1] = 4;      
    this->NbrConnectedOrbitals[2] = 4;      
    if ((this->NextNNHopping != 0.0) || (this->NextNNSpinOrbit != 0.0))
      {       
	this->NbrConnectedOrbitals[0] += 4; 
	this->NbrConnectedOrbitals[1] += 4;      
	this->NbrConnectedOrbitals[2] += 4;      
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

    Complex Lambda1 (-this->NNHopping, -this->NNSpinOrbit);
    Complex Lambda2 (-this->NextNNHopping, -this->NextNNSpinOrbit);

    int TmpIndex = 0;

  // links starting from A
    this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
    this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Lambda1;
    ++TmpIndex;
    this->ConnectedOrbitalIndices[0][TmpIndex] = 2;
    this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Conj(Lambda1);
    ++TmpIndex;
    this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
    this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = -1;
    this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Lambda1;
    ++TmpIndex;
    this->ConnectedOrbitalIndices[0][TmpIndex] = 2;
    this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = -1;
    this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Conj(Lambda1);
    ++TmpIndex;
    
    if ((this->NextNNHopping != 0.0) || (this->NextNNSpinOrbit != 0.0))
    {
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = -1;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Conj(Lambda2);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Conj(Lambda2);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 2;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 1;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Lambda2;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[0][TmpIndex] = 2;
      this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = -1;
      this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) + 1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = Lambda2;
    }

    TmpIndex = 0;;

    // links starting from B
    this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
    this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Conj(Lambda1);
    ++TmpIndex;
    this->ConnectedOrbitalIndices[1][TmpIndex] = 2;
    this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Lambda1;
    ++TmpIndex;
    this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
    this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 1;
    this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Conj(Lambda1);
    ++TmpIndex;
    this->ConnectedOrbitalIndices[1][TmpIndex] = 2;
    this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 1;
    this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = -1;
    this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Lambda1;
    ++TmpIndex;
    
    if ((this->NextNNHopping != 0.0) || (this->NextNNSpinOrbit != 0.0))
    {
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 1;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Lambda2;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Lambda2;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 2;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 1;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Conj(Lambda2);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[1][TmpIndex] = 2;
      this->ConnectedOrbitalSpatialIndices[1][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[1][(TmpIndex * 2) + 1] = -1;
      this->ConnectedOrbitalHoppingAmplitudes[1][TmpIndex] = Conj(Lambda2);
      ++TmpIndex;
    }

    TmpIndex = 0;

    // links starting from C
    this->ConnectedOrbitalIndices[2][TmpIndex] = 0;
    this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = Lambda1;
    ++TmpIndex;
    this->ConnectedOrbitalIndices[2][TmpIndex] = 1;
    this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
    this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = Conj(Lambda1);
    ++TmpIndex;
    this->ConnectedOrbitalIndices[2][TmpIndex] = 0;
    this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 0;
    this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
    this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = Lambda1;
    ++TmpIndex;
    this->ConnectedOrbitalIndices[2][TmpIndex] = 1;
    this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = -1;
    this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
    this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = Conj(Lambda1);
    ++TmpIndex;
    
    if ((this->NextNNHopping != 0.0) || (this->NextNNSpinOrbit != 0.0))
    {
      this->ConnectedOrbitalIndices[2][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 1;
      this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = Conj(Lambda2);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[2][TmpIndex] = 0;
      this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = -1;
      this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = Conj(Lambda2);
      ++TmpIndex;
      this->ConnectedOrbitalIndices[2][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = 0;
      this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 1;
      this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = Lambda2;
      ++TmpIndex;
      this->ConnectedOrbitalIndices[2][TmpIndex] = 1;
      this->ConnectedOrbitalSpatialIndices[2][TmpIndex * 2] = -1;
      this->ConnectedOrbitalSpatialIndices[2][(TmpIndex * 2) + 1] = 0;
      this->ConnectedOrbitalHoppingAmplitudes[2][TmpIndex] = Lambda2;
      ++TmpIndex;
    }

    
    if (this->MuS != 0.0)
      {
	this->ConnectedOrbitalIndices[0][TmpIndex] = 0;
	this->ConnectedOrbitalSpatialIndices[0][TmpIndex * 2] = 0;
	this->ConnectedOrbitalSpatialIndices[0][(TmpIndex * 2) +1] = 0;
	this->ConnectedOrbitalHoppingAmplitudes[0][TmpIndex] = this->MuS;
      }
   }
}