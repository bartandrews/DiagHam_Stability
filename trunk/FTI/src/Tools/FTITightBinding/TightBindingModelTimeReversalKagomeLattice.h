////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                         Class author : Cecile Repellin                     //
//                                                                            //
//         class of tight binding model for the Kagome lattice                //
//                     Time Reversal Invariant Model                          //
//                   last modification : 08/04/2013                           //
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


#ifndef TIGHTBINDINGMODELTIMEREVERSALKAGOMELATTICE_H
#define TIGHTBINDINGMODELTIMEREVERSALKAGOMELATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelTimeReversalKagomeLattice : public Abstract2DTightBindingModel
{

 protected:

  // hopping amplitude between neareast neighbor sites
  double NNHopping;
  // hopping amplitude between next neareast neighbor sites
  double NextNNHopping;
  // spin orbit coupling to neareast neighbor sites
  double NNSpinOrbit;
  // spin orbit coupling to next neareast neighbor sites
  double NextNNSpinOrbit;
  
  // mixingTerm12 = mixing term coupling the two copies of the kagome lattice (sites 1 and 2)
  Complex MixingTerm12;
  // mixingTerm13 = mixing term coupling the two copies of the kagome lattice (sites 1 and 3)
  Complex MixingTerm13;
  // mixingTerm23 = mixing term coupling the two copies of the kagome lattice (sites 2 and 3)
  Complex MixingTerm23;
  
  // four times the sublattice staggered chemical potential 
  double MuS;
  // boundary condition twisting angle along x
  double GammaX;
  // boundary condition twisting angle along y
  double GammaY;
  // nearest neighbor density-density potential strength
  
  // use model with time reversal symmetry
  bool TimeReversal;

 public:

// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mixingTerm12 = mixing term coupling the two copies of the kagome lattice (sites 1 and 2)
// mixingTerm13 = mixing term coupling the two copies of the kagome lattice (sites 1 and 3)
// mixingTerm23 = mixing term coupling the two copies of the kagome lattice (sites 2 and 3)
// mixingTerm11 = mixing term coupling the two copies of the kagome lattice (sites 1up and 1down, 2up and 2down, 3up and 3down)
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelTimeReversalKagomeLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double lambda1, double lambda2, double mixingTerm12, double mixingTerm13, double mixingTerm23, double gammaX, double gammaY, AbstractArchitecture* architecture, bool timeReversalFlag, bool storeOneBodyMatrices);

  // destructor
  //
  ~TightBindingModelTimeReversalKagomeLattice();

  // compute the Bloch hamiltonian at a point of the Brillouin zone
  //
  // kx = momentum along the x axis
  // ky = momentum along the x axis
  // return value = Bloch hamiltonian
  virtual HermitianMatrix ComputeBlochHamiltonian(double kx, double ky);

  // get the high symmetry points 
  //
  // pointNames = name of each high symmetry point
  // pointCoordinates = coordinates in the first Brillouin zone of the high symmetry points
  // return value = number of high symmetry points
  virtual int GetHighSymmetryPoints(char**& pointNames, double**& pointCoordinates);

  // compute the distance between two points in the first Brillouin zone, changing the coordinates the second one by a reciprocal lattice vector if needed
  //
  // kx1 = momentum of the first point along the x axis
  // ky1 = momentum of the first point along the y axis
  // kx2 = reference on the momentum of the second point along the x axis
  // ky2 = reference on the momentum of the second point along the y axis
  // return value = distance between the two points
  virtual double GetDistanceReciprocalSpace(double kx1, double ky1, double& kx2, double& ky2);

 protected :

  // core part that computes the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);

};


#endif
