////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the checkerboard lattice       //
//      with full open boundary conditions and direct implentation of the C4  //
//  symmetry implentation (WANRING: not compatible with may-body hamitonians  //
//  use TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry)         //
//                                                                            //
//                        last modification : 02/03/2018                      //
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
#include "Tools/FTITightBinding/TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>

using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry::TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry()
{
  this->NbrSiteY = 0;
  this->HalfNbrSiteY = 0;
  this->NNHopping = 0.0;
  this->NextNNHopping = 0.0;
  this->NbrBands = 0;
  this->NbrStatePerBand = 0;
  this->Architecture = 0;
  
  this->EnergyBandStructure = 0;
  this->OneBodyBasis = 0;
  this->NbrConfiningPotentials = 0;
  this->ConfiningPotentialCoordinates = 0;
  this->ConfiningPotentialAmplitudes = 0;

  this->ConnectedOrbitalIndices = 0;
  this->ConnectedOrbitalSpatialIndices = 0;
  this->ConnectedOrbitalHoppingAmplitudes = 0;
  this->NbrConnectedOrbitals = 0;
}



// basic constructor
//
// nbrSiteX = number of sites in the x (or y) direction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry::TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry(int nbrSiteX, double t1, double t2, 
															   AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteY = nbrSiteX;
  this->HalfNbrSiteY = this->NbrSiteY / 2;
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->NbrBands = (this->NbrSiteY * this->NbrSiteY) / 4;
  this->NbrStatePerBand = 4;
  this->Architecture = architecture;
  
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[4];
    }
  if (storeOneBodyMatrices == true)
    {
      this->OneBodyBasis = new ComplexMatrix[4];
    }
  this->NbrConfiningPotentials = 0;
  this->ConfiningPotentialCoordinates = 0;
  this->ConfiningPotentialAmplitudes = 0;

  this->ConnectedOrbitalIndices = 0;
  this->ConnectedOrbitalSpatialIndices = 0;
  this->ConnectedOrbitalHoppingAmplitudes = 0;
  this->NbrConnectedOrbitals = 0;
//  this->FindConnectedOrbitals();
  this->ComputeBandStructure();
}

// constructor with an additional confining potential
//
// nbrSiteX = number of sites in the x (or y) direction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// confiningPotentialXCoordinates = x coordiantes of the confining potential
// confiningPotentialYCoordinates = y coordiantes of the confining potential
// confiningPotentialAmplitudes = amplitudes of the confining potential on each sites
// nbrConfiningPotentials = number of sites where there the confining potential has a non-zero amplitude
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry::TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry(int nbrSiteX, double t1, double t2,
															   int* confiningPotentialXCoordinates, int* confiningPotentialYCoordinates, 
															   double* confiningPotentialAmplitudes, int nbrConfiningPotentials,
															   AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteY = nbrSiteX;
  this->HalfNbrSiteY = this->NbrSiteY / 2;
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->NbrBands = (this->NbrSiteY * this->NbrSiteY) / 4;
  this->NbrStatePerBand = 4;
  this->Architecture = architecture;
  
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[4];
    }
  if (storeOneBodyMatrices == true)
    {
      this->OneBodyBasis = new ComplexMatrix[4];
    }

  this->NbrConfiningPotentials = nbrConfiningPotentials;
  this->ConfiningPotentialCoordinates = new int [this->NbrConfiningPotentials];
  this->ConfiningPotentialAmplitudes = new double [this->NbrConfiningPotentials];
  for (int i = 0; i < this->NbrConfiningPotentials; ++i)
    {
      this->ConfiningPotentialCoordinates[i] = this->GetRealSpaceTightBindingLinearizedIndex(confiningPotentialXCoordinates[i], confiningPotentialYCoordinates[i]);
      this->ConfiningPotentialAmplitudes[i] = confiningPotentialAmplitudes[i];
    }
  SortArrayDownOrdering<int>(this->ConfiningPotentialCoordinates, this->ConfiningPotentialAmplitudes, this->NbrConfiningPotentials);

  this->ConnectedOrbitalIndices = 0;
  this->ConnectedOrbitalSpatialIndices = 0;
  this->ConnectedOrbitalHoppingAmplitudes = 0;
  this->NbrConnectedOrbitals = 0;
//  this->FindConnectedOrbitals();
  this->ComputeBandStructure();
}

// destructor
//

TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry::~TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry()
{
  if (this->NbrConfiningPotentials != 0)
    {
      delete[] this->ConfiningPotentialCoordinates;
      delete[] this->ConfiningPotentialAmplitudes;
    }
}

// find the orbitals connected to those located at the origin unit cell in a given discrete symmetry sector
// 
  
void TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry::FindConnectedOrbitals(int sector)
{
  Complex TmpPhaseFactor = this->NNHopping * Phase(M_PI * 0.25);
  Complex TmpConjPhaseFactor = this->NNHopping * Phase(-M_PI * 0.25);

  Complex* TmpC4PhaseFactors = new Complex[4];
  for (int i = 0; i < 4; ++i)
    {
      TmpC4PhaseFactors[i] = Phase(-M_PI * 0.5 * ((double) (sector * i)));
    }

  if (this->NbrConnectedOrbitals != 0)
    {
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  delete[] this->ConnectedOrbitalIndices[i];
	  delete[] this->ConnectedOrbitalSpatialIndices[i];
	  delete[] this->ConnectedOrbitalHoppingAmplitudes[i];
	}
      delete[] this->ConnectedOrbitalIndices;
      delete[] this->ConnectedOrbitalSpatialIndices;
      delete[] this->ConnectedOrbitalHoppingAmplitudes;
      delete[] this->NbrConnectedOrbitals;
    }


  this->NbrConnectedOrbitals = new int [this->NbrBands];
  this->ConnectedOrbitalIndices = new int* [this->NbrBands];
  this->ConnectedOrbitalSpatialIndices = new int* [this->NbrBands];
  this->ConnectedOrbitalHoppingAmplitudes = new Complex* [this->NbrBands];
  for (int x = 0; x < this->HalfNbrSiteY; ++x)
    {
      for (int y = 0; y < this->HalfNbrSiteY; ++y)
	{
	  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndex(x, y);
	  this->NbrConnectedOrbitals[TmpLinearizedCoordinate] = 0;
	}
    }
  int TmpC4Rotation;
  int TmpC4RotationAnnihilation;

  for (int x = 0; x < this->HalfNbrSiteY; ++x)
    {
      for (int y = 0; y < this->HalfNbrSiteY; ++y)
	{
	  int TmpLinearizedCoordinateAnnihilation = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetry(x, y, TmpC4RotationAnnihilation);
	  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
	  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y + 1, TmpC4Rotation);
	  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
	    }
	  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y - 1, TmpC4Rotation);
	  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
	    }
	  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y + 1, TmpC4Rotation);
	  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
	    }
	  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y - 1, TmpC4Rotation);
	  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
	    }
	  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y + 1, TmpC4Rotation);
	  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
	    }
	  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y - 1, TmpC4Rotation);
	  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
	    }
	  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y, TmpC4Rotation);
	  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
	    }
	  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y, TmpC4Rotation);
	  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
	    }
	  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation] = new int[this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]];	  
	  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation] = new int[this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]];
	  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation] = new Complex[this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]];	  
	}
    }
  for (int x = 0; x < this->HalfNbrSiteY; ++x)
    {
      for (int y = 0; y < this->HalfNbrSiteY; ++y)
	{
	  int TmpLinearizedCoordinateAnnihilation = this->GetRealSpaceTightBindingLinearizedIndex(x, y);
	  for (int k = 0; k < this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]; ++k)
	    {
	      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][k] = 0.0;
	    }
	  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation] = 0;
	}
    }


  for (int x = 0; x < this->HalfNbrSiteY; ++x)
    {
      for (int y = 0; y < this->HalfNbrSiteY; ++y)
	{
	  if (((x + y) & 1) == 0)
	    {
	      int TmpLinearizedCoordinateAnnihilation = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetry(x, y, TmpC4RotationAnnihilation);
	      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinateAnnihilation;
	      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
	      if (this->NbrConfiningPotentials > 0)
		{
		  int TmpConfiningPosition = SearchInArrayDownOrdering<int>(TmpLinearizedCoordinateAnnihilation, this->ConfiningPotentialCoordinates, this->NbrConfiningPotentials);
		  if (TmpConfiningPosition != this->NbrConfiningPotentials)
		    {
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += this->ConfiningPotentialAmplitudes[TmpConfiningPosition];
		    }
		}
	      ++this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation];
	      int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y + 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += this->NextNNHopping * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y - 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += this->NextNNHopping * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y + 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += -this->NextNNHopping * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y - 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += -this->NextNNHopping * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y + 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y - 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	    }
	  else
	    {
	      int TmpLinearizedCoordinateAnnihilation = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetry(x, y, TmpC4RotationAnnihilation);
	      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinateAnnihilation;
	      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
	      if (this->NbrConfiningPotentials > 0)
		{
		  int TmpConfiningPosition = SearchInArrayDownOrdering<int>(TmpLinearizedCoordinateAnnihilation, this->ConfiningPotentialCoordinates, this->NbrConfiningPotentials);
		  if (TmpConfiningPosition != this->NbrConfiningPotentials)
		    {
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += this->ConfiningPotentialAmplitudes[TmpConfiningPosition];
		    }
		}
	      ++this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation];
	      int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y + 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += -this->NextNNHopping * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y - 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += -this->NextNNHopping * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y + 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += this->NextNNHopping * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y - 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += this->NextNNHopping * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y + 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y - 1, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	      TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y, TmpC4Rotation);
	      if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor * TmpC4PhaseFactors[TmpC4Rotation];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}
	    }
	}
    }
  
  delete[] TmpC4PhaseFactors;
}

