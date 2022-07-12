////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the simple C4 quadrupole       //
//                  insulator with full open boundary conditions              //
//                  and direct implentation of the C4 symmetry                //
//                                                                            //
//                        last modification : 20/03/2018                      //
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


#ifndef TIGHTBINDINGMODELSIMPLEC4QUADRUPOLEFULLOBCANDFULLC4SYMMETRY_H
#define TIGHTBINDINGMODELSIMPLEC4QUADRUPOLEFULLOBCANDFULLC4SYMMETRY_H


#include "config.h"
#include "Tools/FTITightBinding/TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry.h"


class TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry : public TightBindingModelSquareLatticeFullOBCAndFullC4Symmetry
{

 protected:

  // hoping amplitude between neareast neighbor sites within the unit cell
  double NNHopping1;
  // phase (for the the hoping between neareast neighbor sites within the unit cell
  double NNHoppingPhase1;

  // hoping amplitude between neareast neighbor sites between unit cells
  double NNHopping2;
  // phase for the the hoping between neareast neighbor sites between unit cells
  double NNHoppingPhase2;
  
  // linearized coordiantes of the confining potential
  int* ConfiningPotentialCoordinates;
  // amplitudes of the confining potential on each sites
  double* ConfiningPotentialAmplitudes;
  // number of sites where there the confining potential has a non-zero amplitude
  int NbrConfiningPotentials;


 public:

  // default constructor
  //
  TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry();

  // basic constructor
  //
  // nbrSiteX = number of sites in the x (or y) direction
  // t1 = hoping amplitude between neareast neighbor sites within the unit cell
  // phi1 = phase (in pi units) for the the hoping between neareast neighbor sites within the unit cell
  // t2 = hoping amplitude between neareast neighbor sites between unit cells
  // phi2 = phase (in pi units) for the the hoping between neareast neighbor sites between unit cells
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry(int nbrSiteX, double t1, double phi1, double t2, double phi2, 
							      AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  

  // constructor with an additional confining potential
  //
  // nbrSiteX = number of sites in the x (or y) direction
  // t1 = hoping amplitude between neareast neighbor sites within the unit cell
  // phi1 = phase (in pi units) for the the hoping between neareast neighbor sites within the unit cell
  // t2 = hoping amplitude between neareast neighbor sites between unit cells
  // phi2 = phase (in pi units) for the the hoping between neareast neighbor sites between unit cells
  // confiningPotentialXCoordinates = x coordiantes of the confining potential
  // confiningPotentialYCoordinates = y coordiantes of the confining potential
  // confiningPotentialAmplitudes = amplitudes of the confining potential on each sites
  // nbrConfiningPotentials = number of sites where there the confining potential has a non-zero amplitude
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry(int nbrSiteX, double t1, double phi1, double t2, double phi2,
							      int* confiningPotentialXCoordinates, int* confiningPotentialYCoordinates, 
							      double* confiningPotentialAmplitudes, int nbrConfiningPotentials,
							      AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  

  // destructor
  //
  ~TightBindingModelSimpleC4QuadrupoleFullOBCAndFullC4Symmetry();


 protected :

  // find the orbitals connected to those located at the origin unit cell in a given discrete symmetry sector
  // 
  virtual void FindConnectedOrbitals(int sector);

};


#endif
