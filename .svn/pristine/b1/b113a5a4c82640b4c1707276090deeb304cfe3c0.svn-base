////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the checkerboard lattice       //
//      with full open boundary conditions and direct implentation of the C4  //
//  symmetry implentation (WANRING: not compatible with may-body hamitonians  //
//  use TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry)         //
//                                                                            //
//                        last modification : 02/03/2018                      //
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


#ifndef TIGHTBINDINGMODELCHECKERBOARDLATTICEFULLOBCWITHFULLC4SYMMETRY_H
#define TIGHTBINDINGMODELCHECKERBOARDLATTICEFULLOBCWITHFULLC4SYMMETRY_H


#include "config.h"
#include "Tools/FTITightBinding/TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry.h"


class TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry : public TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry
{

 protected:

  // hoping amplitude between neareast neighbor sites
  double NNHopping;
  // hoping amplitude between next neareast neighbor sites
  double NextNNHopping;
  
  // linearized coordiantes of the confining potential
  int* ConfiningPotentialCoordinates;
  // amplitudes of the confining potential on each sites
  double* ConfiningPotentialAmplitudes;
  // number of sites where there the confining potential has a non-zero amplitude
  int NbrConfiningPotentials;

 public:

  // default constructor
  //
  TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry();
  
  // basic constructor
  //
  // nbrSiteX = number of sites in the x (or y) direction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry(int nbrSiteX, double t1, double t2, 
							       AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  

  // constructor with an additional confining potential
  //
  // nbrSiteX = number of sites in the x (or y) direction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // confiningPotentialXCoordinates = x coordiantes of the confining potential
  // confiningPotentialYCoordinates = y coordiantes of the confining potential
  // confiningPotentialAmplitudes = amplitudes of the confining potential on each sites
  // nbrConfiningPotentials = number of sites where there the confining potential has a non-zero amplitude
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry(int nbrSiteX, double t1, double t2,
							int* confiningPotentialXCoordinates, int* confiningPotentialYCoordinates, 
							double* confiningPotentialAmplitudes, int nbrConfiningPotentials,
							AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  

  // destructor
  //
  ~TightBindingModelCheckerboardLatticeFullOBCAndFullC4Symmetry();

 protected :

  // find the orbitals connected to those located at the origin unit cell in a given discrete symmetry sector
  // 
  virtual void FindConnectedOrbitals(int sector);

};


#endif
