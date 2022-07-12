////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//       class of checkerboard lattice model with interacting particles       //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 13/07/2011                      //
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


#ifndef PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICESINGLEBANDTHREEBODYHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICESINGLEBANDTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian
{

 protected:
  
  // nearest neighbor density-density-density potential strength
  double UPotential;
  // nearest neighbor density-density potential strength
  double VPotential;

  // use flat band model
  bool FlatBand;
  
  // precalculation tables for cosine and sine factors
  Complex* XPhaseTable;
  Complex* YPhaseTable;
  Complex* XHalfPhaseTable;
  Complex* YHalfPhaseTable;
  int XPhaseTableShift;
  int YPhaseTableShift;

 public:

  // default constructor
  //
  ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive three body neareast neighbor interaction
  // vPotential = strength of the repulsive two body neareast neighbor interaction
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, 
								     Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute all the phase precalculation arrays 
  //
  virtual void ComputePhaseArray();

  // compute the matrix element for the two body interaction between two sites A and B 
  //
  // kx1 = momentum along x for the A site
  // ky1 = momentum along y for the A site
  // kx2 = momentum along x for the B site
  // ky2 = momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the three body interaction between two sites A and one site B 
  //
  // kx1 = creation momentum along x for the first A site
  // ky1 = creation momentum along y for the first A site
  // kx2 = creation momentum along x for the second A site
  // ky2 = creation momentum along y for the second A site
  // kx3 = annihilation momentum along x for the first A site
  // ky3 = annihilation momentum along y for the first A site
  // kx4 = annihilation momentum along x for the second A site
  // ky4 = annihilation momentum along y for the second A site
  // return value = corresponding matrix element
  Complex ComputeThreeBodyMatrixElementBAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the three body interaction between one site A and two sites B 
  //
  // kx1 = creation momentum along x for the first B site
  // ky1 = creation momentum along y for the first B site
  // kx2 = creation momentum along x for the second B site
  // ky2 = creation momentum along y for the second B site
  // kx3 = annihilation momentum along x for the first B site
  // ky3 = annihilation momentum along y for the first B site
  // kx4 = annihilation momentum along x for the second B site
  // ky4 = annihilation momentum along y for the second B site
  // return value = corresponding matrix element
  Complex ComputeThreeBodyMatrixElementABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the on-site three body interaction related to sites A
  //
  // kx1 = first creation momentum along x for the first A site
  // ky1 = first creation momentum along y for the first A site
  // kx2 = second creation momentum along x for the second A site
  // ky2 = second creation momentum along y for the second A site
  // kx3 = third creation momentum along x for the second A site
  // ky3 = third creation momentum along y for the second A site
  // kx4 = first annihilation momentum along x for the first A site
  // ky4 = first annihilation momentum along y for the first A site
  // kx5 = second annihilation momentum along x for the second A site
  // ky5 = second annihilation momentum along y for the second A site
  // kx6 = third annihilation momentum along x for the second A site
  // ky6 = third annihilation momentum along y for the second A site
  // return value = corresponding matrix element
  Complex ComputeThreeBodyMatrixElementOnSiteAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the on-site three body interaction related to sites A
  //
  // kx1 = first creation momentum along x for the first A site
  // ky1 = first creation momentum along y for the first A site
  // kx2 = second creation momentum along x for the second A site
  // ky2 = second creation momentum along y for the second A site
  // kx3 = third creation momentum along x for the second A site
  // ky3 = third creation momentum along y for the second A site
  // kx4 = first annihilation momentum along x for the first A site
  // ky4 = first annihilation momentum along y for the first A site
  // kx5 = second annihilation momentum along x for the second A site
  // ky5 = second annihilation momentum along y for the second A site
  // kx6 = third annihilation momentum along x for the second A site
  // ky6 = third annihilation momentum along y for the second A site
  // return value = corresponding matrix element
  Complex ComputeThreeBodyMatrixElementOnSiteBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the on-site three body interaction with two particles on a A site and one on a B site
  //
  // kx1 = first creation momentum along x for the A site
  // ky1 = first creation momentum along y for the A site
  // kx2 = second creation momentum along x for the A site
  // ky2 = second creation momentum along y for the A site
  // kx3 = creation momentum along x for the B site
  // ky3 = creation momentum along y for the B site
  // kx4 = first annihilation momentum along x for the A site
  // ky4 = first annihilation momentum along y for the A site
  // kx5 = second annihilation momentum along x for the A site
  // ky5 = second annihilation momentum along y for the sA site
  // kx6 = annihilation momentum along x for the B site
  // ky6 = annihilation momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeThreeBodyMatrixElementOnSiteAAB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the on-site three body interaction with two particles on a B site and one on a A site
  //
  // kx1 = creation momentum along x for the A site
  // ky1 = creation momentum along y for the A site
  // kx2 = first creation momentum along x for the B site
  // ky2 = first creation momentum along y for the B site
  // kx3 = second creation momentum along x for the B site
  // ky3 = second creation momentum along y for the B site
  // kx4 = annihilation momentum along x for the A site
  // ky4 = annihilation momentum along y for the A site
  // kx5 = first annihilation momentum along x for the B site
  // ky5 = first annihilation momentum along y for the B site
  // kx6 = second annihilation momentum along x for the B site
  // ky6 = second annihilation momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeThreeBodyMatrixElementOnSiteABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the two body interaction between two A sites (or two B sites) 
  //
  // kx1 = momentum along x for the first A site
  // ky1 = momentum along y for the first A site
  // kx2 = momentum along x for the second A site
  // ky2 = momentum along y for the second A site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAA(int kx1, int ky1, int kx2, int ky2);


};



#endif
