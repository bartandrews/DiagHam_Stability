////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of tight binding model for the Haldane honeycomb lattice       //
//                                                                            //
//                        last modification : 30/05/2012                      //
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


#ifndef TIGHTBINDINGMODELHALDANEHONEYCOMBLATTICE_H
#define TIGHTBINDINGMODELHALDANEHONEYCOMBLATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelHaldaneHoneycombLattice : public Abstract2DTightBindingModel
{

 protected:

  // hoping amplitude between neareast neighbor sites
  double NNHopping;
  // hoping amplitude between next neareast neighbor sites
  double NextNNHopping;
  // phase on  next neareast neighbor hopping
  double HaldanePhase;
  
  // four times the sublattice staggered chemical potential 
  double MuS;

 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // phi =  Haldane phase on next neareast neighbor hopping
  // mus = sublattice chemical potential on A sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelHaldaneHoneycombLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double phi, double mus, 
					   double gammaX, double gammaY, 
					   AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  
  
  // constructor for the tilted lattice
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  //nx1 = first coordinate of the first spanning vector for a tilted lattice
  //ny1 = second coordinate of the first spanning vector for a tilted lattice
  //nx2 = first coordinate of the second spanning vector for a tilted lattice
  //ny2 = second coordinate of the second spanning vector for a tilted lattice
  //offset = second coordinate in momentum space of the second spanning vector of the reciprocal lattice for a tilted lattice
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // phi =  Haldane phase on next neareast neighbor hopping
  // mus = sublattice chemical potential on A sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelHaldaneHoneycombLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double phi, double mus, 
					   double gammaX, double gammaY, 
					   AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  // destructor
  //
  ~TightBindingModelHaldaneHoneycombLattice();

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

  // find the orbitals connected to those located at the origin unit cell
  // 
  virtual void FindConnectedOrbitals();

};


#endif
