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


#ifndef TIGHTBINDINGMODELPYROCHLORELATTICE_H
#define TIGHTBINDINGMODELPYROCHLORELATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract3DTightBindingModel.h"


class TightBindingModelPyrochloreLattice : public Abstract3DTightBindingModel
{

 protected:

  // spin orbit coupling to neareast neighbor sites
  double NNSpinOrbit;
  // spin orbit coupling to next neareast neighbor sites
  double NextNNSpinOrbit;

 public:

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
  TightBindingModelPyrochloreLattice(int nbrSiteX, int nbrSiteY, int nbrSiteZ,
				     double lambdaNN, double lambdaNNN,
				     double gammaX, double gammaY, double gammaZ, 
				     AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);

  // destructor
  //
  ~TightBindingModelPyrochloreLattice();

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif
