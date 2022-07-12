////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                          class author: Yang-Le Wu                          //
//                                                                            //
//           class of tight binding model for flux model on FCC lattice       //
//                                                                            //
//                      last modification : 14/11/2012                        //
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


#ifndef TIGHTBINDINGMODEL3DFLUXPYROCHLORELATTICE_H
#define TIGHTBINDINGMODEL3DFLUXPYROCHLORELATTICE_H

#include "config.h"
#include "Tools/FTITightBinding/Abstract3DTightBindingModel.h"

class TightBindingModelFCCLattice : public Abstract3DTightBindingModel
{

protected:

    // kz flux
    double Flux;

    // index of the fractionally filled band
    int BandIndex;

public:
    // default constructor
    //
    // nbrSiteX = number of sites in the x direction
    // nbrSiteY = number of sites in the y direction
    // nbrSiteZ = number of sites in the z direction
    // flux = flux parameter for kz
    // bandIndex = index of the fractionally filled band
    // gammaX = boundary condition twisting angle along x
    // gammaY = boundary condition twisting angle along y
    // gammaZ = boundary condition twisting angle along y
    // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
    TightBindingModelFCCLattice(int nbrSiteX, int nbrSiteY, int nbrSiteZ, double flux, int bandIndex,
            double gammaX, double gammaY, double gammaZ, AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);

    // destructor
    //
    ~TightBindingModelFCCLattice();

protected :

    // core part that compute the band structure
    //
    // minStateIndex = minimum index of the state to compute
    // nbrStates = number of states to compute
    virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif

