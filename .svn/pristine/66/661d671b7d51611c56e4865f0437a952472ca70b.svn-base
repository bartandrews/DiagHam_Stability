////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                         Class author : Cecile Repellin                     //
//                                                                            //
//         class of tight binding model for the Kagome lattice                //
//      Time Reversal Invariant Model with tilted boudary conditions          //
//                   last modification : 03/06/2013                           //
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


#ifndef TIGHTBINDINGMODELTIMEREVERSALKAGOMELATTICETILTED_H
#define TIGHTBINDINGMODELTIMEREVERSALKAGOMELATTICETILTED_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelTimeReversalKagomeLatticeTilted : public Abstract2DTightBindingModel
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
  double MixingTerm12;
  // mixingTerm13 = mixing term coupling the two copies of the kagome lattice (sites 1 and 3)
  double MixingTerm13;
  // mixingTerm23 = mixing term coupling the two copies of the kagome lattice (sites 2 and 3)
  double MixingTerm23;

  // four times the sublattice staggered chemical potential 
  double MuS;
  
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
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelTimeReversalKagomeLatticeTilted(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset,double t1, double t2, double lambda1, double lambda2, double mixingTerm12, double mixingTerm13, double mixingTerm23, double gammaX, double gammaY, AbstractArchitecture* architecture, bool timeReversalFlag, bool storeOneBodyMatrices);

  // destructor
  //
  ~TightBindingModelTimeReversalKagomeLatticeTilted();

 protected :

  // core part that computes the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);
  
  //Computes value of projected momentum along the lattice directions
  //
  //kx = first coordinate of the given point in the Brillouin zone
  //ky = second coordinate of the given point in the Brillouin zone
  //latticeComponent = index of the lattice vector along which the projection is to be performed
  //return value = projected momentum
  virtual double GetProjectedMomentum(int kx, int ky, int latticeComponent);

};

inline double TightBindingModelTimeReversalKagomeLatticeTilted::GetProjectedMomentum(int kx, int ky, int latticeComponent)
  {
    return this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky)][latticeComponent];
  }
#endif
