////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of tight binding model for the Zhang-Qi C=2 model        //
//                         described in arxiv:1403.0164                       //
//                                                                            //
//                        last modification : 30/06/2014                      //
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


#ifndef TIGHTBINDINGMODELZHANGQILATTICE_H
#define TIGHTBINDINGMODELZHANGQILATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelZhangQiLattice : public Abstract2DTightBindingModel
{

 protected:

  // mixing angle between the two layers (0 for decoupled)
  double Theta;
  // azimuthal angle for the vectors connecting next nearest neighboring sites
  double AzimuthalAngle;
  // sublattice chemical potential on A sites
  double MuS;

 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // theta = mixing angle between the two layers in pi unit (0 for decoupled)
  // azimuthalAngle = azimuthal angle (in pi unit) for the vectors connecting next nearest neighboring sites
  // mus = sublattice chemical potential on A sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelZhangQiLattice(int nbrSiteX, int nbrSiteY, double theta, double azimuthalAngle, double muS,
				  double gammaX, double gammaY, 
				  AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);

  // destructor
  //
  ~TightBindingModelZhangQiLattice();

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif
