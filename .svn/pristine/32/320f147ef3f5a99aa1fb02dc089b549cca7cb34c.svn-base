////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//             class of tight binding model for the Kagome lattice            //
//   linearly combining a copy of the model and its time reversal conjugate   //
//                                                                            //
//                        last modification : 11/04/2014                      //
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
#include "Tools/FTITightBinding/TightBindingModelMixedKagomeLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A1 sites
// coefficient = coefficient of the linear combination (0 begin the same model than TightBindingAlternativeMixedKagomeLattice and 1 being its time reversal conjugate)
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelMixedKagomeLattice::TightBindingModelMixedKagomeLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double lambda1, double lambda2, double mus, double coefficient,
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
    this->NNHopping = t1;
    this->NextNNHopping = t2;
    this->NNSpinOrbit = lambda1;
    this->NextNNSpinOrbit = lambda2;
    this->MuS = mus;
    this->TwistAngle = M_PI / 3;
    this->GammaX = gammaX;
    this->GammaY = gammaY;
    this->Coefficient = coefficient;

    this->NbrBands = 3;
    this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
    
    this->ComputeAllProjectedMomenta();

    this->EmbeddingX = RealVector(this->NbrBands, true);
    this->EmbeddingX[1] = 0.5;
    this->EmbeddingY = RealVector(this->NbrBands, true);
    this->EmbeddingY[2] = 0.5;
    this->Inversion = ComplexMatrix(this->NbrBands, this->NbrBands, true);
    for (int i = 0; i < this->NbrBands; ++i)
        this->Inversion[i][i] = 1.0;

    this->Architecture = architecture;

    if (storeOneBodyMatrices == true)
        this->OneBodyBasis = new ComplexMatrix[this->NbrStatePerBand];
    else
        this->OneBodyBasis = 0;
    this->EnergyBandStructure = new double*[this->NbrBands];
    for (int i = 0; i < this->NbrBands; ++i)
        this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    this->ComputeBandStructure();
}


// constructor for tilted lattice
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A1 sites
// coefficient = coefficient of the linear combination (0 begin the same model than TightBindingAlternativeMixedKagomeLattice and 1 being its time reversal conjugate)
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelMixedKagomeLattice::TightBindingModelMixedKagomeLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double lambda1, double lambda2, double mus, double coefficient,
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
    this->NNHopping = t1;
    this->NextNNHopping = t2;
    this->NNSpinOrbit = lambda1;
    this->NextNNSpinOrbit = lambda2;
    this->MuS = mus;
    this->TwistAngle = M_PI / 3;
    this->GammaX = gammaX;
    this->GammaY = gammaY;
    this->Coefficient = coefficient;

    this->NbrBands = 3;
    this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
     
    this->ComputeAllProjectedMomenta();

    this->EmbeddingX = RealVector(this->NbrBands, true);
    this->EmbeddingX[1] = 0.5;
    this->EmbeddingY = RealVector(this->NbrBands, true);
    this->EmbeddingY[2] = 0.5;
    this->Inversion = ComplexMatrix(this->NbrBands, this->NbrBands, true);
    for (int i = 0; i < this->NbrBands; ++i)
        this->Inversion[i][i] = 1.0;

    this->Architecture = architecture;

    if (storeOneBodyMatrices == true)
        this->OneBodyBasis = new ComplexMatrix[this->NbrStatePerBand];
    else
        this->OneBodyBasis = 0;
    this->EnergyBandStructure = new double*[this->NbrBands];
    for (int i = 0; i < this->NbrBands; ++i)
        this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    this->ComputeBandStructure();
}

// destructor
//

TightBindingModelMixedKagomeLattice::~TightBindingModelMixedKagomeLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelMixedKagomeLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
    if (nbrStates == 0l)
        nbrStates = this->NbrStatePerBand;
    long MaxStateIndex = minStateIndex + nbrStates;
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
      {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
	  {
            int Index = this->GetLinearizedMomentumIndex(kx, ky);
	    double x = this->ProjectedMomenta[Index][0];
	    double y = this->ProjectedMomenta[Index][1];
            int InvIndex = this->GetLinearizedMomentumIndex((this->NbrSiteX - kx) % this->NbrSiteX, 
							    (this->NbrSiteY - ky) % this->NbrSiteY);
	    double InvX = this->ProjectedMomenta[InvIndex][0];
	    double InvY = this->ProjectedMomenta[InvIndex][1];
            if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	      {
                Complex nnBA = Complex(-this->NNHopping, -this->NNSpinOrbit) * (1.0 + Phase(x));
                Complex nnCA = Complex(-this->NNHopping, +this->NNSpinOrbit) * (1.0 + Phase(y));
                Complex nnCB = Complex(-this->NNHopping, -this->NNSpinOrbit) * (1.0 + Phase(y-x));
                Complex nnnBA = Complex(-this->NextNNHopping, +this->NextNNSpinOrbit) * (Phase(y) + Phase(x-y));
                Complex nnnCA = Complex(-this->NextNNHopping, -this->NextNNSpinOrbit) * (Phase(x) + Phase(y-x));
                Complex nnnCB = Complex(-this->NextNNHopping, +this->NextNNSpinOrbit) * (Phase(-x) + Phase(y));

                Complex InvnnBA = Conj(Complex(-this->NNHopping, -this->NNSpinOrbit) * (1.0 + Phase(InvX)));
                Complex InvnnCA = Conj(Complex(-this->NNHopping, +this->NNSpinOrbit) * (1.0 + Phase(InvY)));
                Complex InvnnCB = Conj(Complex(-this->NNHopping, -this->NNSpinOrbit) * (1.0 + Phase(InvY - InvX)));
                Complex InvnnnBA = Conj(Complex(-this->NextNNHopping, +this->NextNNSpinOrbit) * (Phase(InvY) + Phase(InvX - InvY)));
                Complex InvnnnCA = Conj(Complex(-this->NextNNHopping, -this->NextNNSpinOrbit) * (Phase(InvX) + Phase(InvY - InvX)));
                Complex InvnnnCB = Conj(Complex(-this->NextNNHopping, +this->NextNNSpinOrbit) * (Phase(-InvX) + Phase(InvY)));

                HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
                TmpOneBodyHamiltonian.SetMatrixElement(0, 0, this->MuS);
                TmpOneBodyHamiltonian.SetMatrixElement(1, 0, this->Coefficient * (nnBA + nnnBA) + (1.0 - this->Coefficient) * (InvnnBA + InvnnnBA));
                TmpOneBodyHamiltonian.SetMatrixElement(2, 0, this->Coefficient * (nnCA + nnnCA) + (1.0 - this->Coefficient) * (InvnnCA + InvnnnCA));
                TmpOneBodyHamiltonian.SetMatrixElement(2, 1, this->Coefficient * (nnCB + nnnCB) + (1.0 - this->Coefficient) * (InvnnCB + InvnnnCB));

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

