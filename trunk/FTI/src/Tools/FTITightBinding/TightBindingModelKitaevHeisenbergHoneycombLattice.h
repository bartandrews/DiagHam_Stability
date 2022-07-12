////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of tight binding model for the Kitaev-Heisenberg honeycomb lattice  //
//                                                                            //
//                        last modification : 01/08/2015                      //
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


#ifndef TIGHTBINDINGMODELKITAEVHEISENBERGHONEYCOMBLATTICE_H
#define TIGHTBINDINGMODELKITAEVHEISENBERGHONEYCOMBLATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelKitaevHeisenbergHoneycombLattice : public Abstract2DTightBindingModel
{

 protected:

  // amplitude of the isotropic nearest neighbor hopping term
  double KineticFactorIsotropic;
  // amplitude of the anisotropic nearest neighbor hopping term
  double KineticFactorAnisotropic;
  // phase of the anisotropic hoping amplitude between neareast neighbor sitesc(in pi units)
  double KineticFactorAnisotropicPhase;
  
  // amplitude of the magnetic field along the x direction
  double HFieldX;
  // amplitude of the magnetic field along the y direction
  double HFieldY;
  // amplitude of the magnetic field along the z direction
  double HFieldZ;

 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // t = isotropic hoping amplitude between neareast neighbor sites
  // tPrime = anisotropic hoping amplitude between next neareast neighbor sites
  // tPrimePhase = additional phase for the anisotropic hoping amplitude between neareast neighbor sitesc(in pi units)
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelKitaevHeisenbergHoneycombLattice(int nbrSiteX, int nbrSiteY, double t, double tPrime, double tPrimePhase, 
						    AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  
  // constructor with magnetic field
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // t = isotropic hoping amplitude between neareast neighbor sites
  // tPrime = anisotropic hoping amplitude between next neareast neighbor sites
  // tPrimePhase = additional phase for the anisotropic hoping amplitude between neareast neighbor sitesc(in pi units)
  // hX = amplitude of the magnetic field along the x direction
  // hY = amplitude of the magnetic field along the y direction
  // hZ = amplitude of the magnetic field along the z direction
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelKitaevHeisenbergHoneycombLattice(int nbrSiteX, int nbrSiteY, double t, double tPrime, double tPrimePhase, double hX, double hY, double hZ,
						    AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  
  // destructor
  //
  ~TightBindingModelKitaevHeisenbergHoneycombLattice();

 protected :

  // find the orbitals connected to those located at the origin unit cell
  // 
  virtual void FindConnectedOrbitals();

};


#endif
