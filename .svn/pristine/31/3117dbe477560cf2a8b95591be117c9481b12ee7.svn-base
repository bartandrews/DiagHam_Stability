////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of tight binding model for the pyrochlore slab lattice        //
//                                                                            //
//                        last modification : 17/05/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelNOrbitalSquareLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrLayers= number of Kagome lattice layer
// t1 = real part of the hopping amplitude between neareast neighbor sites in a given Kagome lattice layer
// t2 = real part of the hopping amplitude between next neareast neighbor sites in a given Kagome lattice layer
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites in a given Kagome lattice layer
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites in a given Kagome lattice layer
// tPerp = hopping term between a triangular lattice layer and a kagome lattice layer
// mus = sublattice chemical potential on A1 sites
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelNOrbitalSquareLattice::TightBindingModelNOrbitalSquareLattice(int nbrSiteX, int nbrSiteY, int nbrLayers,
									       double t1, double t2, double phi, double mus, 
									       double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NbrLayers = nbrLayers;
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->Phi = phi;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = this->NbrLayers;
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

TightBindingModelNOrbitalSquareLattice::~TightBindingModelNOrbitalSquareLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelNOrbitalSquareLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
	      
	       HermitianMatrix TmpOneBodyHamiltonian(this->NbrLayers, true);
	      
	       for (int i = 0; i < this->NbrLayers - 2; ++i)
		{
		  int l=i+1;
		  TmpOneBodyHamiltonian.SetMatrixElement(i,i,2*this->NextNNHopping*cos(KX+KY -(2*l-1)*this->Phi));
		  TmpOneBodyHamiltonian.SetMatrixElement(i,i+1,this->NNHopping*(Phase(-KX) + Phase((KY-2*l*this->Phi))));
		  TmpOneBodyHamiltonian.SetMatrixElement(i,i+2,this->NextNNHopping*(Phase( -(KX-KY+(2*l+1)*this->Phi))));
		}
		int i = (this->NbrLayers-2);
	      int l=i+1;
		  TmpOneBodyHamiltonian.SetMatrixElement(i,i,2*this->NextNNHopping*cos(KX+KY -(2*l-1)*this->Phi));
		  TmpOneBodyHamiltonian.SetMatrixElement(i,i+1,this->NNHopping*(Phase(-KX) + Phase((KY-2*l*this->Phi)) ));
		  TmpOneBodyHamiltonian.SetMatrixElement(i,0,this->NextNNHopping*Phase(-(KX-KY+(2*l+1)*this->Phi)));
		i++;
		 l = this->NbrLayers;
	      TmpOneBodyHamiltonian.SetMatrixElement(i,i,2*this->NextNNHopping*cos(KX+KY -(2*l-1)*this->Phi));
		  TmpOneBodyHamiltonian.SetMatrixElement(i,0,this->NNHopping*(Phase(-KX) + Phase((KY-2*l*this->Phi))));
		  TmpOneBodyHamiltonian.SetMatrixElement(i,1,this->NextNNHopping*(Phase( -(KX-KY+(2*l+1)*this->Phi))));
	      TmpOneBodyHamiltonian.AddToMatrixElement(0,0,-this->MuS);
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
