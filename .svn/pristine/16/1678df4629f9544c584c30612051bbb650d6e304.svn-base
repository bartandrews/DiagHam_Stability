////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//                class of quatum spin Hall restricted to two bands           //
//                      using two decoupled checkerboard models               //
//                      and tilted periodic boundary conditions               //
//                                                                            //
//                        last modification : 02/12/2014                      //
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


#ifndef PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDDECOUPLEDCHECKERBOARDHAMILTONIANTILTED_H
#define PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDDECOUPLEDCHECKERBOARDHAMILTONIANTILTED_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonianTilted.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonianTilted : public ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonianTilted
{

 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // vPotential = strength of the repulsive on site two body interaction between opposite spins
  // wPotential = strength of the repulsive two body neareast neighbor interaction between opposite spins
  // tightBindingModel = pointer to the tight binding model of the first copy (must be in Bloch form)
  // flatBandFlag = use flat band model
  // timeReversalFlag = apply thge time reversal symmetry on the second copy of the tight binding model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonianTilted (ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
									  double uPotential, double vPotential, double wPotential,
									  Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, bool timeReversalFlag, 
									  AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeQuantumSpinHallTwoBandDecoupledCheckerboardHamiltonianTilted();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();
  
  // compute the one body transformation matrices and the optional one body band stucture contribution
  //
  // oneBodyBasis = array of one body transformation matrices (the leftmost upper block for the spin up, the rightmost lower block for the spin down)
  virtual void ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis);
  

  // compute the matrix element for the two body interaction between two sites A and B  belonging to the same layer
  //
  // kx1 = momentum along x for the B site
  // ky1 = momentum along y for the B site
  // kx2 = momentum along x for the B site
  // ky2 = momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the two body interaction between two sites B with different layer indices 
  //
  // kx1 = momentum along x for the A site
  // ky1 = momentum along y for the A site
  // kx2 = momentum along x for the B site
  // ky2 = momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2);
    
  // compute the matrix element for the two body interaction between two sites A and B belonging to different layers
  //
  // kx1 = momentum along x for the A site
  // ky1 = momentum along y for the A site
  // kx2 = momentum along x for the B site
  // ky2 = momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kx2, int ky2);

};



#endif
