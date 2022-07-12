////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of tight binding model for the ruby lattice           //
//                                                                            //
//                        last modification : 01/05/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelRubyLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// tr = real part of the hopping amplitude between neareast neighbor sites with same parity
// ti = imaginary part of the hopping amplitude between neareast neighbor sites with same parity
// t1r = real part of the hopping amplitude next neareast neighbor sites with different parity
// t1i = real part of the hopping amplitude next neareast neighbor sites with different parity
// t4 = hopping amplitude along square diagonal
// mus = sublattice chemical potential on A1 sites
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelRubyLattice::TightBindingModelRubyLattice(int nbrSiteX, int nbrSiteY, double tr, double ti, double t1r, double t1i, double t4, double mus, 
							   double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->Nx1 = this->NbrSiteX;
  this->Ny1 = 0;
  this->Nx2 = 0;
  this->Ny2 = this->NbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->TrHopping = tr;
  this->TiHopping = ti;
  this->T1rHopping = t1r;
  this->T1iHopping = t1i;
  this->T4Hopping = t4;
  this->MuS = mus;
  this->TwistAngle = 2 * M_PI / 3;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 6;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  this->ComputeAllProjectedMomenta();

  double a = 1 + sqrt(3);

  this->EmbeddingX = RealVector(this->NbrBands, true);
  this->EmbeddingX[0] = -sqrt(3) / (3 * a);
  this->EmbeddingX[1] = sqrt(3) / (3 * a);
  this->EmbeddingX[2] = - 1 / a - sqrt(3) / (3 * a);
  this->EmbeddingX[3] = sqrt(3) / (3 * a);
  this->EmbeddingX[4] = -sqrt(3) / (3 * a);
  this->EmbeddingX[5] = sqrt(3) / (3 * a) + 1 / a;

  this->EmbeddingY = RealVector(this->NbrBands, true);
  this->EmbeddingY[0] = -sqrt(3) / (6 * a) + 1 / (2 * a);
  this->EmbeddingY[1] = sqrt(3) / (6 * a) + 1 / (2 * a);
  this->EmbeddingY[2] = -1 / (2 * a) - sqrt(3) / (6 * a);
  this->EmbeddingY[3] = -1 / (2 * a) + sqrt(3) / (6 * a);
  this->EmbeddingY[4] = -1 / (2 * a) - sqrt(3) / (6 * a);
  this->EmbeddingY[5] = sqrt(3) / (6 * a) + 1 / (2 * a);

  this->Inversion = ComplexMatrix(this->NbrBands, this->NbrBands, true);
  for (int i = 0; i < this->NbrBands; ++i)
      this->Inversion[i][i] = 1.0;

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


// constructor for the tilted lattice
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// tr = real part of the hopping amplitude between neareast neighbor sites with same parity
// ti = imaginary part of the hopping amplitude between neareast neighbor sites with same parity
// t1r = real part of the hopping amplitude next neareast neighbor sites with different parity
// t1i = real part of the hopping amplitude next neareast neighbor sites with different parity
// t4 = hopping amplitude along square diagonal
// mus = sublattice chemical potential on A1 sites
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelRubyLattice::TightBindingModelRubyLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double tr, double ti, double t1r, double t1i, double t4, double mus, 
							   double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
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
  this->TrHopping = tr;
  this->TiHopping = ti;
  this->T1rHopping = t1r;
  this->T1iHopping = t1i;
  this->T4Hopping = t4;
  this->MuS = mus;
  this->TwistAngle = 2 * M_PI / 3;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 6;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  
  this->ComputeAllProjectedMomenta();

  double a = 1 + sqrt(3);

  this->EmbeddingX = RealVector(this->NbrBands, true);
  this->EmbeddingX[0] = -sqrt(3) / (3 * a);
  this->EmbeddingX[1] = sqrt(3) / (3 * a);
  this->EmbeddingX[2] = - 1 / a - sqrt(3) / (3 * a);
  this->EmbeddingX[3] = sqrt(3) / (3 * a);
  this->EmbeddingX[4] = -sqrt(3) / (3 * a);
  this->EmbeddingX[5] = sqrt(3) / (3 * a) + 1 / a;

  this->EmbeddingY = RealVector(this->NbrBands, true);
  this->EmbeddingY[0] = -sqrt(3) / (6 * a) + 1 / (2 * a);
  this->EmbeddingY[1] = sqrt(3) / (6 * a) + 1 / (2 * a);
  this->EmbeddingY[2] = -1 / (2 * a) - sqrt(3) / (6 * a);
  this->EmbeddingY[3] = -1 / (2 * a) + sqrt(3) / (6 * a);
  this->EmbeddingY[4] = -1 / (2 * a) - sqrt(3) / (6 * a);
  this->EmbeddingY[5] = sqrt(3) / (6 * a) + 1 / (2 * a);

  this->Inversion = ComplexMatrix(this->NbrBands, this->NbrBands, true);
  for (int i = 0; i < this->NbrBands; ++i)
      this->Inversion[i][i] = 1.0;

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

TightBindingModelRubyLattice::~TightBindingModelRubyLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelRubyLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  Complex CT (this->TrHopping, this->TiHopping);
  Complex CT1 (this->T1rHopping, this->T1iHopping);
  double KX;
  double KY;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int Index = this->GetLinearizedMomentumIndex(kx, ky);
	  if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	    {
	      KX = this->GetProjectedMomentum(kx, ky, 0);
	      KY = this->GetProjectedMomentum(kx, ky, 1);
	      
	      Complex PhaseX = Phase(KX);
	      Complex PhaseXY = Phase(KX + KY);
	      
	      HermitianMatrix TmpOneBodyHamiltonian(6, true);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, -this->MuS);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 2, Conj(CT));
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 4, CT);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 4, Conj(CT));
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 3, Conj(CT));
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 5, CT);
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 5, Conj(CT));
	      
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 1, CT1);
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 4, CT1);
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 5, Conj(CT1) * Conj(PhaseX));
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 2, CT1 * PhaseXY);
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 3, CT1 * Conj(PhaseX));	
	      TmpOneBodyHamiltonian.SetMatrixElement(4, 5, CT1 * Conj(PhaseXY));
	      
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 3, this->T4Hopping * (1.0 + Conj(PhaseX)));
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 4, this->T4Hopping * (1.0 + PhaseXY));
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 5, this->T4Hopping * (Conj(PhaseX) + Conj(PhaseXY)));
	      
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

