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
//                        last modification : 30/01/2013                      //
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
#include "Tools/FTITightBinding/TightBindingModelRandom3D.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <iostream>
#include <ctime>


using std::cout;
using std::endl;


// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// nbrBands = number of bands
// frozen = use the same Hamiltonian at every k
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelRandom3D::TightBindingModelRandom3D(int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrBands, bool frozen,
        AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
    this->NbrSiteX = nbrSiteX;
    this->NbrSiteY = nbrSiteY;
    this->NbrSiteZ = nbrSiteZ;
    this->NbrSiteYZ = this->NbrSiteY * this->NbrSiteZ;
    this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
    this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
    this->KzFactor = 2.0 * M_PI / ((double) this->NbrSiteZ);
    this->GammaX = 0.0;
    this->GammaY = 0.0;
    this->GammaZ = 0.0;
    this->NbrBands = nbrBands;
    this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
    this->Architecture = architecture;

    srand((unsigned)time(NULL));
    this->Frozen = frozen;

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

TightBindingModelRandom3D::~TightBindingModelRandom3D()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelRandom3D::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
    if (nbrStates == 0l)
        nbrStates = this->NbrStatePerBand;
    long MaxStateIndex = minStateIndex + nbrStates;
    HermitianMatrix H(this->NbrBands, true);
    if (this->Frozen)
    {
        for (int i = 0; i < this->NbrBands; ++i)
            for (int j = i; j < this->NbrBands; ++j)
                H.SetMatrixElement(i, j, Complex(((double)rand()/(double)RAND_MAX), ((double)rand()/(double)RAND_MAX)));
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
                    if (!this->Frozen)
                    {
                        for (int i = 0; i < this->NbrBands; ++i)
                            for (int j = i; j < this->NbrBands; ++j)
                                H.SetMatrixElement(i, j, Complex(((double)rand()/(double)RAND_MAX), ((double)rand()/(double)RAND_MAX)));
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
