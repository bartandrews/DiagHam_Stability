////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of tight binding model for the 3D pyrochlore lattice         //
//                                                                            //
//                        last modification : 10/08/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelPyrochloreLattice.h"
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
// lambdaNN = spin orbit coupling to neareast neighbor sites
// lambdaNNN = spin orbit coupling to next neareast neighbor sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// gammaZ = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelPyrochloreLattice::TightBindingModelPyrochloreLattice(int nbrSiteX, int nbrSiteY, int nbrSiteZ,
								       double lambdaNN, double lambdaNNN,
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
  this->NNSpinOrbit = lambdaNN;
  this->NextNNSpinOrbit = lambdaNNN;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->GammaZ = gammaZ;
  this->NbrBands = 8;
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

TightBindingModelPyrochloreLattice::~TightBindingModelPyrochloreLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelPyrochloreLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
    if (nbrStates == 0l)
        nbrStates = this->NbrStatePerBand;
    long MaxStateIndex = minStateIndex + nbrStates;

    // b[i][j] = e_j * b_i / a
    double r = 1.0 / sqrt(2.0);
    double b[4][3] = {
        {r, 0, r},
        {0, r, r},
        {0, 0, 0},
        {r, r, 0}
    };

    // bb1[i][m][n] = e_i * ((b_n - b_m) x (2 * b_n + 2 * b_m - \sum b))
    double bb1[3][4][4];
    for (int m = 0; m < 4; ++m)
    {
        for (int n = 0; n < 4; ++n)
        {
            double d1[3];
            double d2[3];
            for (int i = 0; i < 3; ++i)
            {
                d1[i] = b[n][i] - b[m][i];
                d2[i] = 2 * (b[n][i] + b[m][i]) - 2 * r;
            }
            bb1[0][m][n] = d1[1] * d2[2] - d1[2] * d2[1];
            bb1[1][m][n] = d1[2] * d2[0] - d1[0] * d2[2];
            bb1[2][m][n] = d1[0] * d2[1] - d1[1] * d2[0];
        }
    }

    // bb2[i][m][n][l] = e_i * ((b_m - b_l) x (b_n - b_l))
    double bb2[3][4][4][4];
    for (int m = 0; m < 4; ++m)
    {
        for (int n = 0; n < 4; ++n)
        {
            for (int l = 0; l < 4; ++l)
            {
                double d1[3];
                double d2[3];
                for (int i = 0; i < 3; ++i)
                {
                    d1[i] = b[m][i] - b[l][i];
                    d2[i] = b[n][i] - b[l][i];
                }
                bb2[0][m][n][l] = d1[1] * d2[2] - d1[2] * d2[1];
                bb2[1][m][n][l] = d1[2] * d2[0] - d1[0] * d2[2];
                bb2[2][m][n][l] = d1[0] * d2[1] - d1[1] * d2[0];
            }
        }
    }

    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            for (int kz = 0; kz < this->NbrSiteZ; ++kz)
            {
                int Index = this->GetLinearizedMomentumIndex(kx, ky, kz);
                if ((Index >= minStateIndex) && (Index < MaxStateIndex))
                {
                    double x = this->KxFactor * (((double) kx) + this->GammaX);
                    double y = this->KyFactor * (((double) ky) + this->GammaY);
                    double z = this->KzFactor * (((double) kz) + this->GammaZ);
                    double kb[4] = {x, z, 0, y}; // kb[orbital i] = 2 * k * b_i

                    HermitianMatrix H(this->NbrBands, true); // spin structure (u u u u d d d d)

                    // fill the lower triangle H[n, m]
                    for (int m = 0; m < 4; ++m)
                    {
                        for (int n = 0; n < 4; ++n) // need both n < m and n > m to get all the spin flip terms
                        {
                            if (m == n)
                                continue;
                            Complex t = (1 + Phase(kb[n] - kb[m])); // use +t rather than -t as in Guo&Franz
                            Complex s[3];
                            for (int i = 0; i < 3; ++i)
                            {
                                s[i] = I() * this->NNSpinOrbit / sqrt(2.0) * bb1[i][m][n] * (1 + Phase(kb[n] - kb[m]));
                                for (int l = 0; l < 4; ++l)
                                    if ((l != m) && (l != n))
                                        s[i] += I() * this->NextNNSpinOrbit * bb2[i][m][n][l] * (Phase(kb[n] - kb[l]) + Phase(kb[l] - kb[m]));
                            }
                            H.SetMatrixElement(n, m, t + s[2]);
                            H.SetMatrixElement(n + 4, m + 4, t - s[2]);
                            H.SetMatrixElement(n + 4, m, s[0] + I() * s[1]);
                        }
                    }

                    if (this->OneBodyBasis != 0)
                    {
                        ComplexMatrix TmpMatrix(this->NbrBands, this->NbrBands, true);
                        TmpMatrix.SetToIdentity();
                        RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
                        H.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
                        H.Diagonalize(TmpDiag, TmpMatrix);
#endif
                        this->OneBodyBasis[Index] = TmpMatrix;
                        for (int i = 0; i < this->NbrBands; ++i)
                            this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
                    }
                    else
                    {
                        RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
                        H.LapackDiagonalize(TmpDiag);
#else
                        H.Diagonalize(TmpDiag);
#endif
                        for (int i = 0; i < this->NbrBands; ++i)
                            this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
                    }
                }
            }
        }
    }
}
