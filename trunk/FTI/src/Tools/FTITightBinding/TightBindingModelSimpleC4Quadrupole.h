////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of tight binding model for the simple C4 quadrupole insulator    //
//                                                                            //
//                        last modification : 10/03/2018                      //
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


#ifndef TIGHTBINDINGMODELSIMPLEC4QUADRUPOLE_H
#define TIGHTBINDINGMODELSIMPLEC4QUADRUPOLE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelSimpleC4Quadrupole : public Abstract2DTightBindingModel
{

 protected:

  // hoping amplitude between neareast neighbor sites within the unit cell
  double NNHopping1;
  // phase for the the hoping between neareast neighbor sites within the unit cell
  double NNHoppingPhase1;

  // hoping amplitude between neareast neighbor sites between unit cells
  double NNHopping2;
  // phase for the the hoping between neareast neighbor sites between unit cells
  double NNHoppingPhase2;
  
 public:

  // constructor
  //
  // nbrSiteX = number of unit cells in the x direction
  // nbrSiteY = number of unit cells in the y direction
  // t1 = hoping amplitude between neareast neighbor sites within the unit cell
  // phi1 = phase (in pi units) for the the hoping between neareast neighbor sites within the unit cell
  // t2 = hoping amplitude between neareast neighbor sites between unit cells
  // phi2 = phase (in pi units) for the the hoping between neareast neighbor sites between unit cells
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelSimpleC4Quadrupole(int nbrSiteX, int nbrSiteY, double t1, double phi1, double t2, double phi2,
				      double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  

  // destructor
  //
  ~TightBindingModelSimpleC4Quadrupole();

 protected :

  // find the orbitals connected to those located at the origin unit cell
  // 
  virtual void FindConnectedOrbitals();

};


#endif
