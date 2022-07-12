////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                          class author: Yang-Le Wu                          //
//                                                                            //
//      class of tight binding model for chiral model on cubic lattice        //
//                                                                            //
//                      last modification : 09/09/2012                        //
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
#include "Tools/FTITightBinding/TightBindingModelChiralCubicLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// mus = lambda_7 coefficient
// bandIndex = index of the fractionally filled band
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// gammaZ = boundary condition twisting angle along y
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelChiralCubicLattice::TightBindingModelChiralCubicLattice(int nbrSiteX, int nbrSiteY, int nbrSiteZ, double mus, int bandIndex,
        double gammaX, double gammaY, double gammaZ, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
    this->NbrSiteX = nbrSiteX;
    this->NbrSiteY = nbrSiteY;
    this->NbrSiteZ = nbrSiteZ;
    this->NbrSiteYZ = nbrSiteY * nbrSiteZ;
    this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
    this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
    this->KzFactor = 2.0 * M_PI / ((double) this->NbrSiteZ);
    this->MuS = mus;
    this->BandIndex = bandIndex;
    this->GammaX = gammaX;
    this->GammaY = gammaY;
    this->GammaZ = gammaZ;
    this->NbrBands = 3;
    this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;

    if (storeOneBodyMatrices == true)
    {
        this->OneBodyBasis = new ComplexMatrix[this->NbrStatePerBand];
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

TightBindingModelChiralCubicLattice::~TightBindingModelChiralCubicLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelChiralCubicLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
    if (nbrStates == 0l)
        nbrStates = this->NbrStatePerBand;
    long MaxStateIndex = minStateIndex + nbrStates;
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            for (int kz = 0; kz < this->NbrSiteZ; ++kz)
            {
                int Index = this->GetLinearizedMomentumIndex(kx, ky, kz);
                if ((Index >= minStateIndex) && (Index < MaxStateIndex))
                {
                    double KX = this->KxFactor * (((double) kx) + this->GammaX);
                    double KY = this->KyFactor * (((double) ky) + this->GammaY);
                    double KZ = this->KzFactor * (((double) kz) + this->GammaZ);

                    double d1 = sin(KX);
                    double d2 = sin(KY);
                    double d3 = sin(KZ);
                    double d4 = this->MuS - (cos(KX) + cos(KY) + cos(KZ));

                    HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);

                    // d1 * lambda_4 + d2 * lambda_5
                    TmpOneBodyHamiltonian.SetMatrixElement(0, 2, Complex(d1, -d2));
                    // d3 * lambda_6 + d4 * lambda_7
                    TmpOneBodyHamiltonian.SetMatrixElement(1, 2, Complex(d3, -d4));

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
