////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of generic tight binding model for the square lattice        //
//                with full open boundary conditions and direct               //
//                implentation of the C4 symmetry implentation                //
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
#include "Tools/FTITightBinding/TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry.h"
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

TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry()
{
  this->NbrSiteY = 0;
  this->HalfNbrSiteY = 0;
  this->NbrBands = 0;
  this->NbrStatePerBand = 0;
  this->Architecture = 0;
  
  this->EnergyBandStructure = 0;
  this->OneBodyBasis = 0;

  this->ConnectedOrbitalIndices = 0;
  this->ConnectedOrbitalSpatialIndices = 0;
  this->ConnectedOrbitalHoppingAmplitudes = 0;
  this->NbrConnectedOrbitals = 0;
}

// destructor
//

TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::~TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry()
{
}

// find the orbitals connected to those located at the origin unit cell
// 
  
void TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::FindConnectedOrbitals()
{
  this->FindConnectedOrbitals(0);
}

// find the orbitals connected to those located at the origin unit cell in a given discrete symmetry sector
// 
  
void TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::FindConnectedOrbitals(int sector)
{
  cout << "warning, using dummy method TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::FindConnectedOrbitals" << endl;
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = 4;
  int maxStateIndex = minStateIndex + nbrStates;
  for (int j = minStateIndex; j < maxStateIndex; ++j)
    {
      this->FindConnectedOrbitals(j);
      
      HermitianMatrix TmpOneBodyHamiltonian = this->BuildTightBindingHamiltonianRealSpace(this->NbrConnectedOrbitals, this->ConnectedOrbitalIndices,
											  this->ConnectedOrbitalSpatialIndices, this->ConnectedOrbitalHoppingAmplitudes);
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
	  this->OneBodyBasis[j] = TmpMatrix;
	  for (int i = 0; i < this->NbrBands; ++i)
	    {
	      this->EnergyBandStructure[i][j] = TmpDiag(i, i);
	    }
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
	    this->EnergyBandStructure[i][j] = TmpDiag(i, i);
	}
    }
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  File.precision(14);
  this->WriteASCIIHeader(File, '#');
  File << "# C4 Energy";
  File << endl;
  for (int j= 0; j < 4; ++j)
    {
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  File << j << " " << this->GetEnergy(i, j) << endl;
	}
    }
  File.close();
  return true;
}

