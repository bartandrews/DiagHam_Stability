////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                class of quatum spin Hall restricted to two bands           //
//                           using the kagome model                           //
//                                                                            //
//                        last modification : 17/03/2012                      //
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


#ifndef PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDKAGOMEHAMILTONIAN_H
#define PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDKAGOMEHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnSquareLatticeTwoBandSimpleTIHamiltonian.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian : public ParticleOnSquareLatticeTwoBandSimpleTIHamiltonian
{

 protected:
  
  // nearest neighbor density-density potential strength
  double UPotential;
  // strength of the repulsive on site two body interaction between opposite spins
  double VPotential;
  // strength of the repulsive two body neareast neighbor interaction between opposite spins
  double WPotential;
  // use flat band model
  bool FlatBand;

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
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, double wPotential,Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two spin up on site A
  //
  // kx1 = momentum along x for the creation operator on first A site with spin up
  // ky1 = momentum along y for the creation operator on first A site with spin up
  // kx2 = momentum along x for the creation operator on second A site with spin up
  // ky2 = momentum along y for the creation operator on second A site with spin up
  // kx3 = momentum along x for the annihilation operator on first A site with spin up
  // ky3 = momentum along y for the annihilation operator on first A site with spin up
  // kx4 = momentum along x for the annihilation operator on second A site with spin up
  // ky4 = momentum along y for the annihilation operator on second A site with spin up
  // return value = corresponding matrix element  
  virtual Complex ComputeTwoBodyMatrixElementAUpAUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two spin down on site A
  //
  // kx1 = momentum along x for the creation operator on first A site with spin down
  // ky1 = momentum along y for the creation operator on first A site with spin down
  // kx2 = momentum along x for the creation operator on second A site with spin down
  // ky2 = momentum along y for the creation operator on second A site with spin down
  // kx3 = momentum along x for the annihilation operator on first A site with spin down
  // ky3 = momentum along y for the annihilation operator on first A site with spin down
  // kx4 = momentum along x for the annihilation operator on second A site with spin down
  // ky4 = momentum along y for the annihilation operator on second A site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementADownADown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two opposite spins on site A
  //
  // kx1 = momentum along x for the creation operator on first A site with spin up
  // ky1 = momentum along y for the creation operator on first A site with spin up
  // kx2 = momentum along x for the creation operator on second A site with spin down
  // ky2 = momentum along y for the creation operator on second A site with spin down
  // kx3 = momentum along x for the annihilation operator on first A site with spin up
  // ky3 = momentum along y for the annihilation operator on first A site with spin up
  // kx4 = momentum along x for the annihilation operator on second A site with spin down
  // ky4 = momentum along y for the annihilation operator on second A site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two spin up on site B
  //
  // kx1 = momentum along x for the creation operator on first B site with spin up
  // ky1 = momentum along y for the creation operator on first B site with spin up
  // kx2 = momentum along x for the creation operator on second B site with spin up
  // ky2 = momentum along y for the creation operator on second B site with spin up
  // kx3 = momentum along x for the annihilation operator on first B site with spin up
  // ky3 = momentum along y for the annihilation operator on first B site with spin up
  // kx4 = momentum along x for the annihilation operator on second B site with spin up
  // ky4 = momentum along y for the annihilation operator on second B site with spin up
  // return value = corresponding matrix element  
  virtual Complex ComputeTwoBodyMatrixElementBUpBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two spin down on site B
  //
  // kx1 = momentum along x for the creation operator on first B site with spin down
  // ky1 = momentum along y for the creation operator on first B site with spin down
  // kx2 = momentum along x for the creation operator on second B site with spin down
  // ky2 = momentum along y for the creation operator on second B site with spin down
  // kx3 = momentum along x for the annihilation operator on first B site with spin down
  // ky3 = momentum along y for the annihilation operator on first B site with spin down
  // kx4 = momentum along x for the annihilation operator on second B site with spin down
  // ky4 = momentum along y for the annihilation operator on second B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementBDownBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two opposite spins on site B
  //
  // kx1 = momentum along x for the creation operator on first B site with spin up
  // ky1 = momentum along y for the creation operator on first B site with spin up
  // kx2 = momentum along x for the creation operator on second B site with spin down
  // ky2 = momentum along y for the creation operator on second B site with spin down
  // kx3 = momentum along x for the annihilation operator on first B site with spin up
  // ky3 = momentum along y for the annihilation operator on first B site with spin up
  // kx4 = momentum along x for the annihilation operator on second B site with spin down
  // ky4 = momentum along y for the annihilation operator on second B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two spin up on site C
  //
  // kx1 = momentum along x for the creation operator on first C site with spin up
  // ky1 = momentum along y for the creation operator on first C site with spin up
  // kx2 = momentum along x for the creation operator on second C site with spin up
  // ky2 = momentum along y for the creation operator on second C site with spin up
  // kx3 = momentum along x for the annihilation operator on first C site with spin up
  // ky3 = momentum along y for the annihilation operator on first C site with spin up
  // kx4 = momentum along x for the annihilation operator on second C site with spin up
  // ky4 = momentum along y for the annihilation operator on second C site with spin up
  // return value = corresponding matrix element  
  virtual Complex ComputeTwoBodyMatrixElementCUpCUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two spin down on site C
  //
  // kx1 = momentum along x for the creation operator on first C site with spin down
  // ky1 = momentum along y for the creation operator on first C site with spin down
  // kx2 = momentum along x for the creation operator on second C site with spin down
  // ky2 = momentum along y for the creation operator on second C site with spin down
  // kx3 = momentum along x for the annihilation operator on first C site with spin down
  // ky3 = momentum along y for the annihilation operator on first C site with spin down
  // kx4 = momentum along x for the annihilation operator on second C site with spin down
  // ky4 = momentum along y for the annihilation operator on second C site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementCDownCDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two opposite spins on site C
  //
  // kx1 = momentum along x for the creation operator on first C site with spin up
  // ky1 = momentum along y for the creation operator on first C site with spin up
  // kx2 = momentum along x for the creation operator on second C site with spin down
  // ky2 = momentum along y for the creation operator on second C site with spin down
  // kx3 = momentum along x for the annihilation operator on first C site with spin up
  // ky3 = momentum along y for the annihilation operator on first C site with spin up
  // kx4 = momentum along x for the annihilation operator on second C site with spin down
  // ky4 = momentum along y for the annihilation operator on second C site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementCUpCDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  
};


#endif
