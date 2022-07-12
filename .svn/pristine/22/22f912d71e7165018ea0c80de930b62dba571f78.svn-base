////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//              class of tight binding model for the Kitaev chain             //
//                                                                            //
//                        last modification : 11/06/2016                      //
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


#ifndef TIGHTBINDINGMODELKITAEVCHAIN_H
#define TIGHTBINDINGMODELKITAEVCHAIN_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract1DTightBindingModel.h"


class TightBindingModelKitaevChain : public Abstract1DTightBindingModel
{

 protected:
  
  // hopping between nearest neighboring sites
  double NNHopping;
  // superconducting order parameter 
  double Delta;
  // on site chemical potential
  double MuS;
  
  
 public:
  
  // default constructor
  //
  // nbrSiteX = number of sites
  // hopping = hopping between nearest neighboring sites
  // delta = superconducting order parameter 
  // mus = on site chemical potential
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelKitaevChain(int nbrSiteX, double hopping, double delta, double mus, 
			       AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  
  // destructor
  //
  ~TightBindingModelKitaevChain();
  
 protected :
  
  // find the orbitals connected to those located at the origin unit cell
  // 
  virtual void FindConnectedOrbitals();

};



#endif

