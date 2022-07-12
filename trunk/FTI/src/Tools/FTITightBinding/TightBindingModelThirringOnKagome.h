#////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the Checkerboard lattice       //
//                                                                            //
//                        last modification : 08/05/2012                      //
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


#ifndef TIGHTBINDINGMODELTHIRRINGONKAGOME_H
#define TIGHTBINDINGMODELTHIRRINGONKAGOME_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelThirringOnKagome : public Abstract2DTightBindingModel
{

 protected:

  // hopping amplitude between neareast neighbor sites
  double NNHopping;
  // hopping amplitude between next neareast neighbor sites on A sublattice
  double NextNNHoppingA;
  // hopping amplitude between next neareast neighbor sites on B sublattice
  double NextNNHoppingB;
  // phase of hopping A
  double NextNNHoppingPhase;
  
 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // tNN = hopping amplitude between neareast neighbor sites
  // tNNNa = hopping amplitude between next neareast neighbor sites on A sublattice
  // tNNNb = hopping amplitude between next neareast neighbor sites on B sublattice
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelThirringOnKagome(int nbrSiteX, int nbrSiteY, double tNN, double tNNNa, double tNNNb, double theta, 
				       double gammaX, double gammaY, 
				       AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);

  // destructor
  //
  ~TightBindingModelThirringOnKagome();

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif
