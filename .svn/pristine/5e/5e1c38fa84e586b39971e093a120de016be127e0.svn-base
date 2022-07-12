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
#include "Tools/FTITightBinding/TightBindingModelSimpleC4QuadrupoleFullOBC.h"
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

TightBindingModelSimpleC4QuadrupoleFullOBC::TightBindingModelSimpleC4QuadrupoleFullOBC()
{
  this->NbrSiteX = 0;
  this->NbrSiteY = 0;
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
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = hoping amplitude between neareast neighbor sites within the unit cell
// phi1 = phase (in pi units) for the the hoping between neareast neighbor sites within the unit cell
// t2 = hoping amplitude between neareast neighbor sites between unit cells
// phi2 = phase (in pi units) for the the hoping between neareast neighbor sites between unit cells
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelSimpleC4QuadrupoleFullOBC::TightBindingModelSimpleC4QuadrupoleFullOBC(int nbrSiteX, int nbrSiteY, double t1, double phi1, double t2, double phi2,
										       AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NNHopping1 = t1;
  this->NNHoppingPhase1 = phi1 * M_PI;
  this->NNHopping2 = t2;
  this->NNHoppingPhase2 = phi2 * M_PI;
  this->NbrBands = this->NbrSiteX * this->NbrSiteY;
  this->NbrStatePerBand = 1;
  this->Architecture = architecture;
  
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[1];
    }
  if (storeOneBodyMatrices == true)
    {
      this->OneBodyBasis = new ComplexMatrix[1];
    }
  this->NbrConfiningPotentials = 0;
  this->ConfiningPotentialCoordinates = 0;
  this->ConfiningPotentialAmplitudes = 0;

  this->FindConnectedOrbitals();
  this->ComputeBandStructure();
}

// constructor with an additional confining potential
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction 
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

TightBindingModelSimpleC4QuadrupoleFullOBC::TightBindingModelSimpleC4QuadrupoleFullOBC(int nbrSiteX, int nbrSiteY, double t1, double phi1, double t2, double phi2,
										       int* confiningPotentialXCoordinates, int* confiningPotentialYCoordinates, 
										       double* confiningPotentialAmplitudes, int nbrConfiningPotentials,
										       AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NNHopping1 = t1;
  this->NNHoppingPhase1 = phi1 * M_PI;
  this->NNHopping2 = t2;
  this->NNHoppingPhase2 = phi2 * M_PI;
  this->NbrBands = this->NbrSiteX * this->NbrSiteY;
  this->NbrStatePerBand = 1;
  this->Architecture = architecture;
  
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[1];
    }
  if (storeOneBodyMatrices == true)
    {
      this->OneBodyBasis = new ComplexMatrix[1];
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
  this->FindConnectedOrbitals();
  this->ComputeBandStructure();
}

// destructor
//

TightBindingModelSimpleC4QuadrupoleFullOBC::~TightBindingModelSimpleC4QuadrupoleFullOBC()
{
  if (this->NbrConfiningPotentials != 0)
    {
      delete[] this->ConfiningPotentialCoordinates;
      delete[] this->ConfiningPotentialAmplitudes;
    }
}

// find the orbitals connected to those located at the origin unit cell
// 
  
void TightBindingModelSimpleC4QuadrupoleFullOBC::FindConnectedOrbitals()
{
  Complex TmpPhaseFactor1 = this->NNHopping1 * Phase(this->NNHoppingPhase1);
  Complex TmpPhaseFactor2 = this->NNHopping2 * Phase(this->NNHoppingPhase2);
  Complex TmpConjPhaseFactor1 = this->NNHopping1 * Phase(-this->NNHoppingPhase1);
  Complex TmpConjPhaseFactor2 = this->NNHopping2 * Phase(-this->NNHoppingPhase2);

  this->NbrConnectedOrbitals = new int [this->NbrBands];
  this->ConnectedOrbitalIndices = new int* [this->NbrBands];
  this->ConnectedOrbitalSpatialIndices = new int* [this->NbrBands];
  this->ConnectedOrbitalHoppingAmplitudes = new Complex* [this->NbrBands];
  for (int x = 0; x < this->NbrSiteX; ++x)
    {
      for (int y = 0; y < this->NbrSiteY; ++y)
	{
	  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndex(x, y);
	  this->NbrConnectedOrbitals[TmpLinearizedCoordinate] = 0;
	}
    }
  for (int x = 0; x < this->NbrSiteX; ++x)
    {
      for (int y = 0; y < this->NbrSiteY; ++y)
	{
	  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndex(x, y);
	  if (this->NbrConfiningPotentials > 0)
	    {
	      int TmpConfiningPosition = SearchInArrayDownOrdering<int>(TmpLinearizedCoordinate, this->ConfiningPotentialCoordinates, this->NbrConfiningPotentials);
	      if (TmpConfiningPosition != this->NbrConfiningPotentials)
		{
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		}
	    }
	  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y) >= 0)
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
	    }
	  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x - 1, y) >= 0)
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
	    }
	  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x, y + 1) >= 0)
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
	    }
	  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x, y - 1) >= 0)
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
	    }
	  this->ConnectedOrbitalIndices[TmpLinearizedCoordinate] = new int[this->NbrConnectedOrbitals[TmpLinearizedCoordinate]];	  
	  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate] = new int[this->NbrConnectedOrbitals[TmpLinearizedCoordinate]];
	  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate] = new Complex[this->NbrConnectedOrbitals[TmpLinearizedCoordinate]];	  
	}
    }


  for (int x = 0; x < this->NbrSiteX; ++x)
    {
      for (int y = 0; y < this->NbrSiteY; ++y)
	{
	  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndex(x, y);
	  this->NbrConnectedOrbitals[TmpLinearizedCoordinate] = 0;
	  if (this->NbrConfiningPotentials > 0)
	    {
	      int TmpConfiningPosition = SearchInArrayDownOrdering<int>(TmpLinearizedCoordinate, this->ConfiningPotentialCoordinates, this->NbrConfiningPotentials);
	      if (TmpConfiningPosition != this->NbrConfiningPotentials)
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpLinearizedCoordinate;
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->ConfiningPotentialAmplitudes[TmpConfiningPosition];
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		}
	    }
	  if ((x & 1) == 0)
	    {
	      if ((y & 1) == 0)
		{		  
		  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y) >= 0)
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x + 1, y);
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpConjPhaseFactor1;
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		    }
		  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x, y + 1) >= 0)
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x, y + 1);
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpPhaseFactor1;
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		    }
		}
	      else
		{
		  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y) >= 0)
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x + 1, y);
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpPhaseFactor1;
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		    }
		  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x, y + 1) >= 0)
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x, y + 1);
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpConjPhaseFactor2;
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		    }
		}	      
	    }
	  else
	    {
	      if ((y & 1) == 0)
		{		  
		  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y) >= 0)
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x + 1, y);
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpPhaseFactor2;
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		    }
		  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x, y + 1) >= 0)
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x, y + 1);
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpConjPhaseFactor1;
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		    }
		}
	      else
		{
		  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y) >= 0)
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x + 1, y);
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpConjPhaseFactor2;
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		    }
		  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x, y + 1) >= 0)
		    {
		      this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x, y + 1);
		      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpPhaseFactor2;
		      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		    }
		}	      
	    }
	}
    }
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool TightBindingModelSimpleC4QuadrupoleFullOBC::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  File.precision(14);
  this->WriteASCIIHeader(File, '#');
  File << "# Energy";
  File << endl;
  for (int i = 0; i < this->NbrBands; ++i)
    File << this->GetEnergy(i, 0) << endl;
  File.close();
  return true;
}

