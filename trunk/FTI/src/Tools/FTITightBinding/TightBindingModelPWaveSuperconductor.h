////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of tight binding model for p-wave superconductor on a lattice     //
//                                                                            //
//                        last modification : 27/04/2015                      //
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


#ifndef TIGHTBINDINGMODELPWAVESUPERCONDUCTOR_H
#define TIGHTBINDINGMODELPWAVESUPERCONDUCTOR_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelPWaveSuperconductor : public Abstract2DTightBindingModel
{

 protected:

  // hoping amplitude between neareast neighbor sites
  double NNHoping;
  // amplitude of the superconducting order parameter
  double Delta;
  // chemical potential 
  double Mu;

 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // t = hoping amplitude between neareast neighbor sites
  // delta = amplitude pf the superconducting order parameter
  // mu = chemical potential
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelPWaveSuperconductor(int nbrSiteX, int nbrSiteY, double t, double delta, double mu, 
				       double gammaX, double gammaY, 
				       AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  
  
  // destructor
  //
  ~TightBindingModelPWaveSuperconductor();

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif
