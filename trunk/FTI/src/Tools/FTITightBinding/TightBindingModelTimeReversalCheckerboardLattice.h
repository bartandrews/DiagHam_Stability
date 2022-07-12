////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                         Class author : Cecile Repellin                     //
//                                                                            //
//         class of tight binding model for the Checkerboard lattice          //
//                     Time Reversal Invariant Model                          //
//                   last modification : 17/04/2013                           //
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


#ifndef TIGHTBINDINGMODELTIMEREVERSALCHECKERBOARDLATTICE_H
#define TIGHTBINDINGMODELTIMEREVERSALCHECKERBOARDLATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelTimeReversalCheckerboardLattice : public Abstract2DTightBindingModel
{

 protected:

  // hoping amplitude between neareast neighbor sites
  double NNHoping;
  // hoping amplitude between next neareast neighbor sites
  double NextNNHoping;
  // hoping amplitude between second next neareast neighbor sites
  double SecondNextNNHoping;
  // mixingTermNorm = norm of the mixing term coupling the two copies of the checkerboard lattice
  double MixingTermNorm ;
  // mixingTermArg = argument of the mixing term coupling the two copies of the checkerboard lattice
  double MixingTermArg ;

  // four times the sublattice staggered chemical potential 
  double MuS;
  // boundary condition twisting angle along x
  double GammaX;
  // boundary condition twisting angle along y
  double GammaY;
  // nearest neighbor density-density potential strength
  

 public:

// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// t2p = hoping amplitude between second next neareast neighbor sites
// mus = sublattice chemical potential on A sites
// mixingTerm = mixing term coupling the two copies of the checkerboard lattice
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelTimeReversalCheckerboardLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double t2p, double mus, double mixingTermNorm, double mixingTermArg, double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices);

  // destructor
  //
  ~TightBindingModelTimeReversalCheckerboardLattice();

 protected :

  // core part that computes the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif
