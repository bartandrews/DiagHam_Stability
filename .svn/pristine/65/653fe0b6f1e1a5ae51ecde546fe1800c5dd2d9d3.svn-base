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
#include "Tools/FTITightBinding/TightBindingModelChern2DiceLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t = nearest neighbor hopping amplitude
// epsilon = on site energy for site 3
// lambda = Rashba spin orbit coupling strength
// bfield1 = magnetic field strength on sites 1 and 2
// bfield3 = magnetic field strength on site 3
// mus = sublattice chemical potential on A1 sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelChern2DiceLattice::TightBindingModelChern2DiceLattice(int nbrSiteX, int nbrSiteY, double t, double epsilon, double lambda, double bfield1, double bfield3, double mus, 
								       double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->THopping = t;
  this->Epsilon = epsilon;
  this->Lambda = lambda;
  this->BField1 = bfield1;
  this->BField3 = bfield3;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 6;
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

TightBindingModelChern2DiceLattice::~TightBindingModelChern2DiceLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelChern2DiceLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
	      
	      double angle = 2.0 * M_PI / 3.0;
	      Complex GammaK = this->THopping * (1.0 + Phase(-KX) + Phase(-KY));
	      Complex GammaP = this->Lambda * (1.0 + Phase(-KX + angle) + Phase(-KY - angle));
	      Complex GammaM = this->Lambda * (1.0 + Phase(-KX - angle) + Phase(-KY + angle));
	      
	      HermitianMatrix TmpOneBodyHamiltonian(6, true);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 4, Conj(GammaK));
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 5, Conj(GammaK));
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 5, I()*Conj(GammaP));
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 4, I()*Conj(GammaM));
	      
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 4, GammaK);
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 5, GammaK);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 5, I()*GammaM);
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 4, I()*GammaP);
	      
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, this->BField1);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -this->BField1);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 2, this->BField1);
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 3, -this->BField1);
	      TmpOneBodyHamiltonian.SetMatrixElement(4, 4, this->Epsilon + this->BField3);
	      TmpOneBodyHamiltonian.SetMatrixElement(5, 5, this->Epsilon - this->BField3);
	      
	      TmpOneBodyHamiltonian *= -1.0;
	      
	      if (this->OneBodyBasis != 0)
		{
		  ComplexMatrix TmpMatrix(6, 6, true);
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

