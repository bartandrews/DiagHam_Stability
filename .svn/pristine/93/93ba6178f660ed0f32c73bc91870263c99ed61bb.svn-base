////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                        class author: Yang-Le Wu                            //
//                                                                            //
//  class of tight binding model for the two-orbital model on square lattice  //
//                                                                            //
//                        last modification : 02/10/2012                      //
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


#ifndef TIGHTBINDINGMODELTWOORBITALSQUARELATTICE_H
#define TIGHTBINDINGMODELTWOORBITALSQUARELATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelTwoOrbitalSquareLattice : public Abstract2DTightBindingModel
{
protected:

    // imag part of the inter-orbital hopping amplitude between nearest neighbors along the x direction
    double NNHoppingInterX;
    // the inter-orbital hopping amplitude between nearest neighbors along the y direction
    double NNHoppingInterY;
    // the intra-orbital hopping amplitude between nearest neighbors
    double NNHoppingIntra;
    // folding factor for the momenta along sigma_x and sigma_y
    double FoldingFactor;  
    // four times the sublattice staggered chemical potential 
    double MuS;

public:

    // default constructor
    //
    // nbrSiteX = number of sites in the x direction
    // nbrSiteY = number of sites in the y direction
    // t1 = imag part of the inter-orbital hopping amplitude between nearest neighbors along the x direction
    // t2 = the inter-orbital hopping amplitude between nearest neighbors along the y direction
    // t3 = the intra-orbital hopping amplitude between nearest neighbors
    // foldingFactor = folding factor for the momenta along sigma_x and sigma_y
    // mus = sublattice chemical potential on A sites
    // gammaX = boundary condition twisting angle along x
    // gammaY = boundary condition twisting angle along y
    // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
    TightBindingModelTwoOrbitalSquareLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double t3, int foldingFactor, double mus, double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);

    // destructor
    //
    ~TightBindingModelTwoOrbitalSquareLattice();

protected :

    // compute the band structure
    //
    // minStateIndex = minimum index of the state to compute
    // nbrStates = number of states to compute
    virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif

