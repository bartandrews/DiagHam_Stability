////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of tight binding model for the 3D atomic limit         //
//                                                                            //
//                        last modification : 25/08/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModel3DAtomicLimitLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// nbrSitesUnitCell = number of sites per unit cell
// chemicalPotentials = on site chemical potentials
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// gammaZ = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModel3DAtomicLimitLattice::TightBindingModel3DAtomicLimitLattice(int nbrSiteX, int nbrSiteY, int nbrSiteZ,
									     int nbrSitesUnitCell, double* chemicalPotentials,
									     double gammaX, double gammaY, double gammaZ, 
									     AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteYZ = this->NbrSiteY * this->NbrSiteZ;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->KzFactor = 2.0 * M_PI / ((double) this->NbrSiteZ);
  this->NbrSitesUnitCell = nbrSitesUnitCell;
  this->ChemicalPotentials = new double[this->NbrSitesUnitCell];
  for (int i = 0; i < this->NbrSitesUnitCell; ++i)
    this->ChemicalPotentials[i] = chemicalPotentials[i];
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->GammaZ = gammaZ;
  this->NbrBands = this->NbrSitesUnitCell;
  this->NbrStatePerBand =  this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
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

// constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// nbrSitesUnitCell = number of sites per unit cell
// NbrLowChemicalPotentials = number of sites with chemical potential -1, the other sites having chemical potential +1
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// gammaZ = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModel3DAtomicLimitLattice::TightBindingModel3DAtomicLimitLattice(int nbrSiteX, int nbrSiteY, int nbrSiteZ,
									     int nbrSitesUnitCell, int nbrLowChemicalPotentials,
									     double gammaX, double gammaY, double gammaZ, 
									     AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->KzFactor = 2.0 * M_PI / ((double) this->NbrSiteZ);
  this->NbrSitesUnitCell = nbrSitesUnitCell;
  this->ChemicalPotentials = new double[this->NbrSitesUnitCell];
  for (int i = 0; i < nbrLowChemicalPotentials; ++i)
    this->ChemicalPotentials[i] = -1.0;
  for (int i = nbrLowChemicalPotentials; i < this->NbrSitesUnitCell; ++i)
    this->ChemicalPotentials[i] = 1.0;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->GammaZ = gammaZ;
  this->NbrBands = this->NbrSitesUnitCell;
  this->NbrStatePerBand =  this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;

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

TightBindingModel3DAtomicLimitLattice::~TightBindingModel3DAtomicLimitLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModel3DAtomicLimitLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  double KX;
  double KY;
  double KZ;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  for (int kz = 0; kz < this->NbrSiteZ; ++kz)
	    {
	      int Index = this->GetLinearizedMomentumIndex(kx, ky, kz);
	      if ((Index >= minStateIndex) && (Index < MaxStateIndex))
		{
		  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
		  for (int i = 0; i < this->NbrBands; ++i)
		    TmpOneBodyHamiltonian.SetMatrixElement(i, i, this->ChemicalPotentials[i]);

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
}
