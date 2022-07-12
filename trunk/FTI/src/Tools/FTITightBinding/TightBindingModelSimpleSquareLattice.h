////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of tight binding model for the simple square lattice       //
//                                                                            //
//                        last modification : 01/01/2016                      //
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


#ifndef TIGHTBINDINGMODELSINGLESQUARELATTICE_H
#define TIGHTBINDINGMODELSINGLESQUARELATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelSimpleSquareLattice : public Abstract2DTightBindingModel
{

 protected:

  // hoping amplitude between neareast neighbor sites
  double NNHoping;
  // hoping amplitude between next neareast neighbor sites
  double NextNNHoping;

 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelSimpleSquareLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double gammaX, double gammaY, 
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
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelSimpleSquareLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, 
				       double gammaX, double gammaY, 
				       AbstractArchitecture* architecture, int offsetReal = 0, bool storeOneBodyMatrices = true);

  // destructor
  //
  ~TightBindingModelSimpleSquareLattice();

  // evaluate the two point correlation function 
  //
  // x = linearized position index of the first point
  // y = linearized position index of the second point
  // occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
  // nbrOccupiedMomenta = number of occupied momenta
  // bandIndex = index of the band to consider
  // return value = value of the two point correlation function 
  virtual Complex EvaluateTwoPointCorrelationFunction(int x, int y, int* occupiedMomenta, int nbrOccupiedMomenta, int bandIndex);

  // evaluate the two point correlation function in a given region
  //
  // maxX = x coordinate of the region upper right corner 
  // maxY = y coordinate of the region upper right corner 
  // occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
  // nbrOccupiedMomenta = number of occupied momenta
  // bandIndex = index of the band to consider
  // return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)
  virtual HermitianMatrix EvaluateFullTwoPointCorrelationFunction(int maxX, int maxY, int* occupiedMomenta, int nbrOccupiedMomenta, int bandIndex);

 protected :

  // find the orbitals connected to those located at the origin unit cell
  // 
  virtual void FindConnectedOrbitals();

};


#endif
