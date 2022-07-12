////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of tight binding model for the C=2 dice lattice         //
//                                                                            //
//                        last modification : 08/05/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelChern2TriangularLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = nearest neighbor hopping amplitude
// mus = sublattice chemical potential on A1 sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelChern2TriangularLattice::TightBindingModelChern2TriangularLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double phi, double mus, 
										   double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->Phi = 2.0 * M_PI * phi;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 3;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
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

TightBindingModelChern2TriangularLattice::~TightBindingModelChern2TriangularLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelChern2TriangularLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
	      KX = this->KxFactor * (((double) kx) + this->GammaX);
	      KY = this->KyFactor * (((double) ky) + this->GammaY);
	      
	      Complex HAB = this->NNHopping * (1.0 + Phase(KX) + Phase(KY));
	      Complex HAC = -this->NNHopping * (Phase(-2.0 * this->Phi) + Phase(KY) + Phase(2.0 * this->Phi + KY - KX));
	      Complex HBC = this->NNHopping * (Phase(2.0 * this->Phi) + Phase(KY - KX - 2.0 * this->Phi) + Phase(-KX));
	      
	      Complex HAA = 2.0 * this->NextNNHopping * (cos(KY + this->Phi) + cos(KX - this->Phi) + cos(KX - KY + this->Phi)) - this->MuS;
	      Complex HBB = 2.0 * this->NextNNHopping * (cos(KY - this->Phi) + cos(KX + this->Phi) + cos(KX - KY - this->Phi));
	      Complex HCC = -2.0 * this->NextNNHopping * (cos(KY) + cos(KX) + cos(KX - KY));
	      
	      HermitianMatrix TmpOneBodyHamiltonian(3, true);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 1, HAB);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 2, HAC);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 2, HBC);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, HAA);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 1, HBB);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 2, HCC);

	      
	      if (this->OneBodyBasis != 0)
		{
		  ComplexMatrix TmpMatrix(3, 3, true);
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

