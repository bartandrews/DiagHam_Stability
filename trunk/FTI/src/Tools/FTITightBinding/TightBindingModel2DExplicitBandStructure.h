////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of dummy tight binding model for the 2D lattice        //
//                  where the band structure is explicitly provided           //
//                                                                            //
//                        last modification : 21/05/2020                      //
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


#ifndef TIGHTBINDINGMODEL2DEXPLICITBANDSTRUCTURE_H
#define TIGHTBINDINGMODEL2DEXPLICITBANDSTRUCTURE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModel2DExplicitBandStructure : public Abstract2DTightBindingModel
{

 protected:


 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrBands = number of bands (identical to the number of orbitals per unit cell)
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // kxMomenta = array providing the momenta along x for the band structure
  // kyMomenta = array providing the momenta along y for the band structure
  // energies = array for the band structure energies (first index is the same than kxMomenta/kyMomenta, second index is the band index)
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModel2DExplicitBandStructure(int nbrSiteX, int nbrSiteY, int nbrBands,
					   double gammaX, double gammaY,
					   int* kxMomenta, int* kyMomenta, double** energies,
					   AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);

  // destructor
  //
  ~TightBindingModel2DExplicitBandStructure();

  // compute the Bloch hamiltonian at a point of the Brillouin zone
  //
  // kx = momentum along the x axis
  // ky = momentum along the x axis
  // return value = Bloch hamiltonian
  virtual HermitianMatrix ComputeBlochHamiltonian(double kx, double ky);

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif
