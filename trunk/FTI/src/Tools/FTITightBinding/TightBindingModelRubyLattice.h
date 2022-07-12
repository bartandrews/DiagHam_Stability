////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of tight binding model for the ruby lattice           //
//                                                                            //
//                        last modification : 01/05/2012                      //
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


#ifndef TIGHTBINDINGMODELRUBYLATTICE_H
#define TIGHTBINDINGMODELRUBYLATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelRubyLattice : public Abstract2DTightBindingModel
{

 protected:

  // real part of the hopping amplitude between neareast neighbor sites with same parity
  double TrHopping;
  // imaginary part of the hopping amplitude between neareast neighbor sites with same parity
  double TiHopping;
  // real part of the hopping amplitude next neareast neighbor sites with different parity
  double T1rHopping;
  // real part of the hopping amplitude next neareast neighbor sites with different parity
  double T1iHopping;
  // t4 = hopping amplitude along square diagonal
  double T4Hopping;

  // four times the sublattice staggered chemical potential 
  double MuS;

 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // tr = real part of the hopping amplitude between neareast neighbor sites with same parity
  // ti = imaginary part of the hopping amplitude between neareast neighbor sites with same parity
  // t1r = real part of the hopping amplitude next neareast neighbor sites with different parity
  // t1i = real part of the hopping amplitude next neareast neighbor sites with different parity
  // t4 = hopping amplitude along square diagonal
  // mus = sublattice chemical potential on A1 sites
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelRubyLattice(int nbrSiteX, int nbrSiteY, double tr, double ti, double t1r, double t1i, double t4, double mus, 
			       double gammaX, double gammaY, 
			       AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  
  // constructor for tilted lattice
  //
  //nx1 = first coordinate of the first spanning vector for a tilted lattice
  //ny1 = second coordinate of the first spanning vector for a tilted lattice
  //nx2 = first coordinate of the second spanning vector for a tilted lattice
  //ny2 = second coordinate of the second spanning vector for a tilted lattice
  //offset = second coordinate in momentum space of the second spanning vector of the reciprocal lattice for a tilted lattice
  TightBindingModelRubyLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double tr, double ti, double t1r, double t1i, double t4, double mus, double gammaX, double gammaY, AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);

  // destructor
  //
  ~TightBindingModelRubyLattice();

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif
