////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                         Class author : Cecile Repellin                     //
//                                                                            //
//         class of tight binding model for the Checkerboard lattice          //
//                     Time Reversal Invariant Model                          //
//                   last modification : 17/04/2013                           //
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
#include "Tools/FTITightBinding/TightBindingModelTimeReversalCheckerboardLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// t2p = hoping amplitude between second next neareast neighbor sites
// mus = sublattice chemical potential on A sites
// mixingTerm = mixing term coupling the two copies of the checkerboard lattice
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelTimeReversalCheckerboardLattice::TightBindingModelTimeReversalCheckerboardLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double t2p, double mus, double mixingTermNorm, double mixingTermNormArg, double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->SecondNextNNHoping = t2p;
  this->MuS = mus;
  this->MixingTermNorm = mixingTermNorm;
  this->MixingTermArg = mixingTermNormArg;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 4;
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

TightBindingModelTimeReversalCheckerboardLattice::~TightBindingModelTimeReversalCheckerboardLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelTimeReversalCheckerboardLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
	      Complex B1 = 4.0 * this->NNHoping * Complex (cos (1.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (1.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY)) * cos(M_PI * 0.25), 
							   sin (1.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * sin (1.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY)) * sin(M_PI * 0.25));
	      double d1 = 4.0 * this->SecondNextNNHoping * cos (2.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (2.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY));
	      double d3 =  this->MuS + (2.0 * this->NextNNHoping * (cos (2.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX))
								    - cos (2.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY))));
	      
	      
	      Complex B1Conj = 4.0 * this->NNHoping * Complex (cos (1.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (1.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY)) * cos(M_PI * 0.25), 
							   -sin (1.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX)) * sin (1.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY)) * sin(M_PI * 0.25));
	      double d1Conj = 4.0 * this->SecondNextNNHoping * cos (2.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (2.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY));
	      double d3Conj =  this->MuS + (2.0 * this->NextNNHoping * (cos (2.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX))
								    - cos (2.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY))));

	      
	      
	      HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	      
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 2, d1Conj + d3Conj);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 3, B1Conj);
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 3, d1Conj - d3Conj);
	      
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 3, this->MixingTermNorm * Phase(this->MixingTermArg));
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 2, -this->MixingTermNorm * Phase(this->MixingTermArg));
	      

	      
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



