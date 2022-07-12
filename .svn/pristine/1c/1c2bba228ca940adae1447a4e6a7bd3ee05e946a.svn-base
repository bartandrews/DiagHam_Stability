////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the simple C4 quadrupole       //
//                  insulator with full open boundary conditions              //
//                  and direct implentation of the C4 symmetry                //
//                                                                            //
//                        last modification : 20/03/2018                      //
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
#include "Tools/FTITightBinding/TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry.h"
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

TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry::TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry()
{
  this->NbrSiteY = 0;
  this->HalfNbrSiteY = 0;
  this->NNHopping1 = 0.0;
  this->NNHoppingPhase1 = 0.0;
  this->NNHopping2 = 0.0;
  this->NNHoppingPhase2 = 0.0;
  this->NbrBands = 0;
  this->NbrStatePerBand = 0;
  this->Architecture = 0;
  
  this->EnergyBandStructure = 0;
  this->OneBodyBasis = 0;
  this->NbrConfiningPotentials = 0;
  this->ConfiningPotentialCoordinates = 0;
  this->ConfiningPotentialAmplitudes = 0;
}



// basic constructor
//
// nbrSiteX = number of sites in the x (or y) direction
// t1 = hoping amplitude between neareast neighbor sites within the unit cell
// phi1 = phase (in pi units) for the the hoping between neareast neighbor sites within the unit cell
// t2 = hoping amplitude between neareast neighbor sites between unit cells
// phi2 = phase (in pi units) for the the hoping between neareast neighbor sites between unit cells
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry::TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry(int nbrSiteX, double t1, double phi1, double t2, double phi2,
															 AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteY = nbrSiteX;
  this->HalfNbrSiteY = this->NbrSiteY / 2;
  this->NNHopping1 = t1;
  this->NNHoppingPhase1 = phi1 * M_PI;
  this->NNHopping2 = t2;
  this->NNHoppingPhase2 = phi2 * M_PI;
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

  this->ComputeBandStructure();
}

// constructor with an additional confining potential
//
// nbrSiteX = number of sites in the x (or y) direction
// t1 = hoping amplitude between neareast neighbor sites within the unit cell
// phi1 = phase (in pi units) for the the hoping between neareast neighbor sites within the unit cell
// t2 = hoping amplitude between neareast neighbor sites between unit cells
// phi2 = phase (in pi units) for the the hoping between neareast neighbor sites between unit cells
// confiningPotentialXCoordinates = x coordiantes of the confining potential
// confiningPotentialYCoordinates = y coordiantes of the confining potential
// confiningPotentialAmplitudes = amplitudes of the confining potential on each sites
// nbrConfiningPotentials = number of sites where there the confining potential has a non-zero amplitude
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry::TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry(int nbrSiteX, double t1, double phi1, double t2, double phi2,
															 int* confiningPotentialXCoordinates, int* confiningPotentialYCoordinates, 
															 double* confiningPotentialAmplitudes, int nbrConfiningPotentials,
															 AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteY = nbrSiteX;
  this->HalfNbrSiteY = this->NbrSiteY / 2;
  this->NNHopping1 = t1;
  this->NNHoppingPhase1 = phi1 * M_PI;
  this->NNHopping2 = t2;
  this->NNHoppingPhase2 = phi2 * M_PI;
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
  this->ComputeBandStructure();
}

// destructor
//

TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry::~TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry()
{
  if (this->NbrConfiningPotentials != 0)
    {
      delete[] this->ConfiningPotentialCoordinates;
      delete[] this->ConfiningPotentialAmplitudes;
    }
}

// find the orbitals connected to those located at the origin unit cell in a given discrete symmetry sector
// 

void TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry::FindConnectedOrbitals(int sector)
{
  Complex TmpPhaseFactor1 = this->NNHopping1 * Phase(this->NNHoppingPhase1);
  Complex TmpPhaseFactor2 = this->NNHopping2 * Phase(this->NNHoppingPhase2);
  Complex TmpConjPhaseFactor1 = this->NNHopping1 * Phase(-this->NNHoppingPhase1);
  Complex TmpConjPhaseFactor2 = this->NNHopping2 * Phase(-this->NNHoppingPhase2);

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
	  if (this->NbrConfiningPotentials > 0)
	    {
	      int TmpConfiningPosition = SearchInArrayDownOrdering<int>(TmpLinearizedCoordinateAnnihilation, this->ConfiningPotentialCoordinates, this->NbrConfiningPotentials);
	      if (TmpConfiningPosition != this->NbrConfiningPotentials)
		{
		   this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}	      
	    }
	  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y, TmpC4Rotation);
	  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
	    }
	  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y, TmpC4Rotation);
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
	  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation] = new int[this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]];	  
	  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation] = new int[this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]];
	  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation] = new Complex[this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]];	  
	}
    }
  
  for (int x = 0; x < this->HalfNbrSiteY; ++x)
    {
      for (int y = 0; y < this->HalfNbrSiteY; ++y)
	{
	  int TmpLinearizedCoordinateAnnihilation = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetry(x, y, TmpC4RotationAnnihilation);
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
	  int TmpLinearizedCoordinateAnnihilation = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetry(x, y, TmpC4RotationAnnihilation);
	  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation] = 0;
	  if (this->NbrConfiningPotentials > 0)
	    {
	      int TmpConfiningPosition = SearchInArrayDownOrdering<int>(TmpLinearizedCoordinateAnnihilation, this->ConfiningPotentialCoordinates, this->NbrConfiningPotentials);
	      if (TmpConfiningPosition != this->NbrConfiningPotentials)
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinateAnnihilation;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += this->ConfiningPotentialAmplitudes[TmpConfiningPosition];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		}	      
	    }
	  if ((x & 1) == 0)
	    {
	      if ((y & 1) == 0)
		{		  
		  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor1 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor2 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y + 1, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor1 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y - 1, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor2 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		}
	      else
		{
		  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor1 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor2 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y - 1, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor1 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y + 1, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor2 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		}	      
	    }
	  else
	    {
	      if ((y & 1) == 0)
		{		  
		  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor1 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor2 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y + 1, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor1 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y - 1, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor2 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		}
	      else
		{
		  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x - 1, y, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor1 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x + 1, y, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpConjPhaseFactor2 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y - 1, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor1 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		  TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndexAndDiscreteSymmetrySafe(x, y + 1, TmpC4Rotation);
		  if ((TmpLinearizedCoordinate >= 0) && (TmpLinearizedCoordinate <= TmpLinearizedCoordinateAnnihilation))
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = TmpLinearizedCoordinate;
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinateAnnihilation][this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]] += TmpPhaseFactor2 * TmpC4PhaseFactors[TmpC4Rotation];
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinateAnnihilation]++;
		    }
		}	      
	    }
	}
    }
  delete[] TmpC4PhaseFactors;
}

