////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the simple C4 quadrupole       //
//             insulator with full open boundary conditions compatible        //
//                         with the C4 symmetry implentation                  //
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
#include "Tools/FTITightBinding/TightBindingModelSimpleC4QuadrupoleFullOBCC4Symmetry.h"
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

TightBindingModelSimpleC4QuadrupoleFullOBCC4Symmetry::TightBindingModelSimpleC4QuadrupoleFullOBCC4Symmetry()
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
// nbrSiteX = number of sites in the x (or y) direction
// t1 = hoping amplitude between neareast neighbor sites within the unit cell
// phi1 = phase (in pi units) for the the hoping between neareast neighbor sites within the unit cell
// t2 = hoping amplitude between neareast neighbor sites between unit cells
// phi2 = phase (in pi units) for the the hoping between neareast neighbor sites between unit cells
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelSimpleC4QuadrupoleFullOBCC4Symmetry::TightBindingModelSimpleC4QuadrupoleFullOBCC4Symmetry(int nbrSiteX, double t1, double phi1, double t2, double phi2,
													   AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteX;
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

  this->CreateLinearizedIndexMap();
  this->NbrConfiningPotentials = 0;
  this->ConfiningPotentialCoordinates = 0;
  this->ConfiningPotentialAmplitudes = 0;

  this->FindConnectedOrbitals();
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

TightBindingModelSimpleC4QuadrupoleFullOBCC4Symmetry::TightBindingModelSimpleC4QuadrupoleFullOBCC4Symmetry(int nbrSiteX, double t1, double phi1, double t2, double phi2,
													   int* confiningPotentialXCoordinates, int* confiningPotentialYCoordinates, 
													   double* confiningPotentialAmplitudes, int nbrConfiningPotentials,
													   AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteX;
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

  this->CreateLinearizedIndexMap();

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

TightBindingModelSimpleC4QuadrupoleFullOBCC4Symmetry::~TightBindingModelSimpleC4QuadrupoleFullOBCC4Symmetry()
{
  delete[] this->LinearizedIndexMap;
  delete[] this->InvertLinearizedIndexMap;
}

// create the linearized index map to go from the canonical linearized index to the C4 compatible one
//

void TightBindingModelSimpleC4QuadrupoleFullOBCC4Symmetry::CreateLinearizedIndexMap()
{
  this->LinearizedIndexMap = new int[this->NbrBands];
  this->InvertLinearizedIndexMap = new int[this->NbrBands];
  if ((this->NbrSiteX & 1) == 0)
    {
      int TmpHalf = this->NbrSiteX / 2;
      int TotalIndex = 0;
      for (int x = 0; x < TmpHalf; ++x)
	{
	  for (int y = 0; y < TmpHalf; ++y)
	    {
	      this->LinearizedIndexMap[y  + (x * this->NbrSiteY)] = TotalIndex;	      
	      ++TotalIndex;
	    }
	}
      for (int x = 0; x < TmpHalf; ++x)
	{
	  for (int y = 0; y < TmpHalf; ++y)
	    {
	      this->LinearizedIndexMap[x  + ((this->NbrSiteX - 1 - y) * this->NbrSiteY)] = TotalIndex;	      
	      ++TotalIndex;
	    }
	}
      for (int x = 0; x < TmpHalf; ++x)
	{
	  for (int y = 0; y < TmpHalf; ++y)
	    {
	      this->LinearizedIndexMap[(this->NbrSiteX - 1 - y)  + ((this->NbrSiteX - 1 - x) * this->NbrSiteY)] = TotalIndex;	      
	      ++TotalIndex;
	    }
	}
      for (int x = 0; x < TmpHalf; ++x)
	{
	  for (int y = 0; y < TmpHalf; ++y)
	    {
	      this->LinearizedIndexMap[(this->NbrSiteX - 1 - x)  + (y * this->NbrSiteY)] = TotalIndex;	      
	      ++TotalIndex;	      
	    }
	}
    }
  else
    {
      int TmpHalf = this->NbrSiteX / 2;
      int TotalIndex = 0;
      for (int x = 0; x <= TmpHalf; ++x)
	{
	  for (int y = 0; y < TmpHalf; ++y)
	    {
	      this->LinearizedIndexMap[y  + (x * this->NbrSiteY)] = TotalIndex;	      
	      ++TotalIndex;
	    }
	}
      for (int x = 0; x <= TmpHalf; ++x)
	{
	  for (int y = 0; y < TmpHalf; ++y)
	    {
	      this->LinearizedIndexMap[x  + ((this->NbrSiteX - 1 - y) * this->NbrSiteY)] = TotalIndex;	      
	      ++TotalIndex;
	    }
	}
      for (int x = 0; x <= TmpHalf; ++x)
	{
	  for (int y = 0; y < TmpHalf; ++y)
	    {
	      this->LinearizedIndexMap[(this->NbrSiteX - 1 - y)  + ((this->NbrSiteX - 1 - x) * this->NbrSiteY)] = TotalIndex;	      
	      ++TotalIndex;
	    }
	}
      for (int x = 0; x <= TmpHalf; ++x)
	{
	  for (int y = 0; y < TmpHalf; ++y)
	    {
	      this->LinearizedIndexMap[(this->NbrSiteX - 1 - x)  + (y * this->NbrSiteY)] = TotalIndex;	      
	      ++TotalIndex;	      
	    }
	}
      this->LinearizedIndexMap[TmpHalf  + (TmpHalf * this->NbrSiteY)] = TotalIndex;	      
      ++TotalIndex;
    }
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->InvertLinearizedIndexMap[i] = SearchInUnsortedArray<int>(i, this->LinearizedIndexMap, this->NbrBands);
    }  
}
