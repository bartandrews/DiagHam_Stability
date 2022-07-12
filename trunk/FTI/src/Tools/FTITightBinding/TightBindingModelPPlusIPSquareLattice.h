////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the Chern insulator of         //
//                 PHYSICAL REVIEW B 84, 075128 (2011) equivalent             //
//                           to a p+ip superconductor                         //
//                                                                            //
//                        last modification : 19/02/2020                      //
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


#ifndef TIGHTBINDINGMODELPPLUSIPSQUARELATTICE_H
#define TIGHTBINDINGMODELPPLUSIPSQUARELATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/TightBindingModelCheckerboardLattice.h"


class TightBindingModelPPlusIPSquareLattice : public TightBindingModelCheckerboardLattice
{

 protected:

  public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // mus = sublattice chemical potential on A sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelPPlusIPSquareLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double mus, 
					double gammaX, double gammaY, 
					AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  
  // constructor for a tilted lattice
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nx1 = first coordinate of the first spanning vector of the tilted lattice
  // nx2 = second coordinate of the first spanning vector of the tilted lattice
  // ny1 = first coordinate of the second spanning vector of the tilted lattice
  // ny2 = second coordinate of the second spanning vector of the tilted lattice
  // offset = second coordinate in momentum space of the second spanning vector of the reciprocal lattice (0 if lattice is untilted or if Ny = 1)
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // mus = sublattice chemical potential on A sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelPPlusIPSquareLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double mus, 
					double gammaX, double gammaY, 
					AbstractArchitecture* architecture, int offsetReal = 0, bool storeOneBodyMatrices = true);

  // destructor
  //
  ~TightBindingModelPPlusIPSquareLattice();

 protected :

  // find the orbitals connected to those located at the origin unit cell
  // 
  virtual void FindConnectedOrbitals();

};


#endif
