////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of tight binding model for a simple chain with confinementlattice  //
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


#ifndef TIGHTBINDINGMODELSSH_H
#define TIGHTBINDINGMODELSSH_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract1DTightBindingModel.h"


class TightBindingModelSSH : public Abstract1DTightBindingModel
{

 protected:

  // hoping amplitude between neareast neighbor sites
  double Delta;

  // square confinement amplitude 
  bool CylinderFlag;

 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // t1 = hoping amplitude between neareast neighbor sites
  // epsilon = square confinement amplitude 
  // gammaX = boundary condition twisting angle along x
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelSSH (int nbrSite, double delta, bool cylinder, AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  
  // destructor
  //
  ~TightBindingModelSSH();

  // evaluate the two point correlation function in a given region
  //
  // maxX = x coordinate of the region upper right corner 
  // maxY = y coordinate of the region upper right corner 
  // occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
  // nbrOccupiedMomenta = number of occupied momenta
  // bandIndex = index of the band to consider
  // return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)
  virtual HermitianMatrix GetRealSpaceTightBindingHamiltonian();

};


#endif
