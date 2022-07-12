////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of tight binding model for the 3D  simple TI lattice         //
//                                                                            //
//                        last modification : 27/09/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModel3DSimpleTILattice.h"
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
// nbrSiteZ = number of sites in the z direction
// mass = mass term of the simple TI model
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// gammaZ = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModel3DSimpleTILattice::TightBindingModel3DSimpleTILattice(int nbrSiteX, int nbrSiteY, int nbrSiteZ,
								       double mass,
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
  this->Mass = mass;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->GammaZ = gammaZ;
  this->NbrBands = 4;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
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

// destructor
//

TightBindingModel3DSimpleTILattice::~TightBindingModel3DSimpleTILattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModel3DSimpleTILattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
		  KX = this->KxFactor * (((double) kx) + this->GammaX);
		  KY = this->KyFactor * (((double) ky) + this->GammaY);
		  KZ = this->KzFactor * (((double) kz) + this->GammaZ);
		  
		  Complex d2 (sin (KY), -sin (KZ));
		  double d1 = sin (KX);
		  double d3 = (this->Mass - cos (KX) - cos (KY) - cos (KZ));
		  
		  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
		  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d3);
		  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -d3);
		  TmpOneBodyHamiltonian.SetMatrixElement(2, 2, d3);
		  TmpOneBodyHamiltonian.SetMatrixElement(3, 3, -d3);
		  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, d1);
		  TmpOneBodyHamiltonian.SetMatrixElement(2, 3, d1);
		  TmpOneBodyHamiltonian.SetMatrixElement(0, 3, d2);
		  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, d2);


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
