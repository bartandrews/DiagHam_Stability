////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of tight binding model for 3D frozen h(k)              //
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
#include "Tools/FTITightBinding/TightBindingModelFrozen3D.h"
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
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelFrozen3D::TightBindingModelFrozen3D(int nbrSiteX, int nbrSiteY, int nbrSiteZ,
        ComplexMatrix oneBodyBasis, double* oneBodyEnergy,
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
    this->NbrBands = oneBodyBasis.GetNbrRow();
    this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
    this->Architecture = architecture;

    this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
    this->EnergyBandStructure = new double*[this->NbrBands];

    for (int i = 0; i < this->NbrBands; ++i)
    {
        this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
        for (int j = 0; j < this->NbrStatePerBand; ++j)
            this->EnergyBandStructure[i][j] = oneBodyEnergy[i];
    }

    for (int j = 0; j < this->NbrStatePerBand; ++j)
        this->OneBodyBasis[j] = oneBodyBasis;
}

// destructor
//

TightBindingModelFrozen3D::~TightBindingModelFrozen3D()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelFrozen3D::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
}
