////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of dummy tight binding model for the 2D lattice        //
//                  where the band structure is explicitly provided           //
//                                                                            //
//                        last modification : 21/05/2020                      //
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
#include "Tools/FTITightBinding/TightBindingModel2DExplicitBandStructure.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrBands = number of bands (identical to the number of orbitals per unit cell)
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// kxMomenta = array providing the momenta along x for the band structure
// kyMomenta = array providing the momenta along y for the band structure
// energies = array for the band structure energies (first index is the same than kxMomenta/kyMomenta, second index is the band index)
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModel2DExplicitBandStructure::TightBindingModel2DExplicitBandStructure(int nbrSiteX, int nbrSiteY, int nbrBands,
										   double gammaX, double gammaY,
										   int* kxMomenta, int* kyMomenta, double** energies,
										   AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = nbrBands;
  this->NbrStatePerBand =  this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;

  this->ComputeBandStructure();
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

  for (int i = 0; i < this->NbrStatePerBand; ++i)
    {
      int Index = this->GetLinearizedMomentumIndex(kxMomenta[i], kyMomenta[i]);
      for (int j = 0; j < this->NbrBands; ++j)
	{
	  this->EnergyBandStructure[j][Index] = energies[i][j];
	  if (storeOneBodyMatrices == true)
	    {
	      this->OneBodyBasis[Index] = ComplexMatrix(this->NbrBands, this->NbrBands, true);
	      this->OneBodyBasis[Index].SetToIdentity();
	    }
	}
   }
}

// destructor
//

TightBindingModel2DExplicitBandStructure::~TightBindingModel2DExplicitBandStructure()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModel2DExplicitBandStructure::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
}


// compute the Bloch hamiltonian at a point of the Brillouin zone
//
// kx = momentum along the x axis
// ky = momentum along the x axis
// return value = Bloch hamiltonian

HermitianMatrix TightBindingModel2DExplicitBandStructure::ComputeBlochHamiltonian(double kx, double ky)
{
  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
  int Index = this->GetLinearizedMomentumIndex(kx, ky);
  for (int i = 0; i < this->NbrBands; ++i)
    {
      TmpOneBodyHamiltonian.SetMatrixElement(i, i, this->EnergyBandStructure[i][Index]);
    }  
  return TmpOneBodyHamiltonian;
}
