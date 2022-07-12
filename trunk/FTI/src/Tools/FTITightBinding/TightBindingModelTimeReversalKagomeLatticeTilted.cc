////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                         Class author : Cecile Repellin                     //
//                                                                            //
//         class of tight binding model for the Kagome lattice                //
//      Time Reversal Invariant Model with tilted boudary conditions          //
//                   last modification : 03/06/2013                           //
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
#include "Tools/FTITightBinding/TightBindingModelTimeReversalKagomeLatticeTilted.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

using std::cout;
using std::endl;

// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mixingTerm12 = mixing term coupling the two copies of the kagome lattice (sites 1 and 2)
// mixingTerm13 = mixing term coupling the two copies of the kagome lattice (sites 1 and 3)
// mixingTerm23 = mixing term coupling the two copies of the kagome lattice (sites 2 and 3)
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelTimeReversalKagomeLatticeTilted::TightBindingModelTimeReversalKagomeLatticeTilted(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double lambda1, double lambda2, double mixingTerm12, double mixingTerm13, double mixingTerm23, double gammaX, double gammaY, AbstractArchitecture* architecture, bool timeReversalFlag, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->Nx1 = nx1;
  this->Ny1 = ny1;
  this->Nx2 = nx2;
  this->Ny2 = ny2;
  this->Offset = offset;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->NNSpinOrbit = lambda1;
  this->NextNNSpinOrbit = lambda2;
  this->MixingTerm12 = mixingTerm12;
  this->MixingTerm13 = mixingTerm13;
  this->MixingTerm23 = mixingTerm23;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 6;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  this->TimeReversal = timeReversalFlag;
  
  this->ComputeAllProjectedMomenta();
  

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

TightBindingModelTimeReversalKagomeLatticeTilted::~TightBindingModelTimeReversalKagomeLatticeTilted()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelTimeReversalKagomeLatticeTilted::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
	  
	  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	  
	  KX = this->ProjectedMomenta[Index][0];
	  KY = this->ProjectedMomenta[Index][1];
	  Complex HAB (-this->NNHopping, -this->NNSpinOrbit);
	  HAB *= 1 + Phase(KX);
	  Complex HAC(-this->NNHopping, this->NNSpinOrbit);
	  HAC *= 1 + Phase(KY);
	  Complex HBC(- this->NNHopping, - this->NNSpinOrbit);
	  HBC *= 1 + Phase(KY - KX);
	  
	  
	  double InvKX = -this->ProjectedMomenta[Index][0];
	  double InvKY = -this->ProjectedMomenta[Index][1];
	  Complex InvHAB = Complex(- this->NNHopping, - this->NNSpinOrbit);
	  InvHAB *= 1 + Phase(InvKX);
	  Complex InvHAC = Complex(- this->NNHopping,  this->NNSpinOrbit);
	  InvHAC *= 1 + Phase(InvKY);
	  Complex InvHBC = Complex(- this->NNHopping, - this->NNSpinOrbit);
	  InvHBC *= 1 + Phase(InvKY - InvKX);
	  
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, HAB);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 2, HAC);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, HBC);
	  if (this->TimeReversal == true)
	    {
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 4, Conj(InvHAB));
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 5, Conj(InvHAC));
	      TmpOneBodyHamiltonian.SetMatrixElement(4, 5, Conj(InvHBC));
	    }
	  else
	    {
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 4, HAB);
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 5, HAC);
	      TmpOneBodyHamiltonian.SetMatrixElement(4, 5, HBC);
	    }
	  
	  
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 4, this->MixingTerm12);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 5, this->MixingTerm13);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 3, -this->MixingTerm12);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 5, this->MixingTerm23);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 3, -this->MixingTerm13);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 4, -this->MixingTerm23);
	  
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



