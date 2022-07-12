////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of ruby lattice model with interacting particles          //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 25/10/2011                      //
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



#ifndef PARTICLEONLATTICECHERN2TRIANGULARLATTICESINGLEBANDTHREEBODYHAMILTONIAN_H
#define PARTICLEONLATTICECHERN2TRIANGULARLATTICESINGLEBANDTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian : public ParticleOnLatticeKagomeLatticeSingleBandThreeBodyHamiltonian
{

 protected:
  
 public:

  // default constructor
  //
  ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // threeBodyPotential = strength of the repulsive three body neareast neighbor interaction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // vPotential = strength of the repulsive two body next neareast neighbor interaction
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, double threeBodyPotential, double uPotential, double vPotential, Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);
  
  // destructor
  //
  ~ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian();
  

 protected:
 
  // compute the matrix element for the two body interaction between two sites A and B 
  //
  // kx1 = momentum along x for the A site
  // ky1 = momentum along y for the A site
  // kx2 = momentum along x for the B site
  // ky2 = momentum along y for the B site
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2);

    // compute the matrix element for the two body interaction between two sites A and C
  //
  // kx1 = creation momentum along x for the C site
  // ky1 = creation momentum along y for the C site
  // kx2 = annihilation momentum along x for the C site
  // ky2 = annihilation momentum along y for the C site
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAC(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the two body interaction between two sites B and C 
  //
  // k1a = creation momentum along e_a for the B site
  // k1b = creation momentum along e_b for the B site
  // k2a = creation momentum along e_a for the A site
  // k2b = creation momentum along e_b for the A site
  // k3a = annihilation momentum along e_a for the B site
  // k3b = annihilation momentum along e_b for the B site
  // k4a = annihilation momentum along e_a for the A site
  // k4b = annihilation momentum along e_b for the A site
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementBC(int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b);

    // compute the matrix element for on-site two body interaction involving A sites
  //
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementOnSiteAA();

  // compute the matrix element for on-site two body interaction involving B sites
  //
  // kx1 = first creation momentum along x for the B site
  // ky1 = first creation momentum along y for the B site
  // kx2 = second creation momentum along x for the B site
  // ky2 = second creation momentum along y for the B site
  // kx3 = first annihilation momentum along x for the B site
  // ky3 = first annihilation momentum along y for the B site
  // kx4 = second annihilation momentum along x for the B site
  // ky4 = second annihilation momentum along y for the B site
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for on-site two body interaction involving C sites
  //
  // kx1 = first creation momentum along x for the C site
  // ky1 = first creation momentum along y for the C site
  // kx2 = second creation momentum along x for the C site
  // ky2 = second creation momentum along y for the C site
  // kx3 = first annihilation momentum along x for the C site
  // ky3 = first annihilation momentum along y for the C site
  // kx4 = second annihilation momentum along x for the C site
  // ky4 = second annihilation momentum along y for the C site
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementOnSiteCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

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
  virtual Complex ComputeThreeBodyMatrixElementOnSiteAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

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
  virtual Complex ComputeThreeBodyMatrixElementOnSiteBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  virtual Complex ComputeThreeBodyMatrixElementOnSiteCCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);
    

  // compute the matrix element for the creation part of the three body on site interaction for the A1 sites 
  //
  // kx1 = momentum along x of the first creation operator
  // ky1 = momentum along y of the first creation operator
  // kx2 = momentum along x of the second creation operator
  // ky2 = momentum along y of the second creation operator
  // kx3 = momentum along x of the third creation operator
  // ky3 = momentum along y of the third creation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteAAAIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body on site interaction for the A1 sites 
  //
  // kx4 = momentum along x of the first annihilation operator
  // ky4 = momentum along y of the first annihilation operator
  // kx5 = momentum along x of the second annihilation operator
  // ky5 = momentum along y of the secondannihilation operator
  // kx6 = momentum along x of the third annihilation operator
  // ky6 = momentum along y of the third annihilation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteAAAOut(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body on site interaction for the A2 sites 
  //
  // kx1 = momentum along x of the first creation operator
  // ky1 = momentum along y of the first creation operator
  // kx2 = momentum along x of the second creation operator
  // ky2 = momentum along y of the second creation operator
  // kx3 = momentum along x of the third creation operator
  // ky3 = momentum along y of the third creation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteBBBIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body on site interaction for the A2 sites 
  //
  // kx4 = momentum along x of the first annihilation operator
  // ky4 = momentum along y of the first annihilation operator
  // kx5 = momentum along x of the second annihilation operator
  // ky5 = momentum along y of the secondannihilation operator
  // kx6 = momentum along x of the third annihilation operator
  // ky6 = momentum along y of the third annihilation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteBBBOut(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body on site interaction for the A3 sites 
  //
  // kx1 = momentum along x of the first creation operator
  // ky1 = momentum along y of the first creation operator
  // kx2 = momentum along x of the second creation operator
  // ky2 = momentum along y of the second creation operator
  // kx3 = momentum along x of the third creation operator
  // ky3 = momentum along y of the third creation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteCCCIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body on site interaction for the A3 sites 
  //
  // kx4 = momentum along x of the first annihilation operator
  // ky4 = momentum along y of the first annihilation operator
  // kx5 = momentum along x of the second annihilation operator
  // ky5 = momentum along y of the secondannihilation operator
  // kx6 = momentum along x of the third annihilation operator
  // ky6 = momentum along y of the third annihilation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteCCCOut(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);


  // compute the matrix element for the creation part of the three body interaction between sites A, B and C, for triangle standing up 
  //
  // kx1 = momentum along x of the creation operator of the A site
  // ky1 = momentum along y of the creation operator of the A site
  // kx2 = momentum along x of the creation operator of the B site
  // ky2 = momentum along y of the creation operator of the B site
  // kx3 = momentum along x of the creation operator of the C site
  // ky3 = momentum along y of the creation operator of the C site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementABCUpIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);
  // same for triangle pointing down
  virtual Complex ComputeThreeBodyMatrixElementABCDownIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A, B and C , for triangle standing up 
  //
  // kx4 = momentum along x of the annihilation operator of the A site
  // ky4 = momentum along y of the annihilation operator of the A site
  // kx5 = momentum along x of the annihilation operator of the B site
  // ky5 = momentum along y of the annihilation operator of the B site
  // kx6 = momentum along x of the annihilation operator of the C site
  // ky6 = momentum along y of the annihilation operator of the C site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementABCUpOut(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);
  // same for triangle pointing down
  virtual Complex ComputeThreeBodyMatrixElementABCDownOut(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

};


// compute the matrix element for the creation part of the three body interaction between sites A, B and C, for triangle standing up 
//
// kx1 = momentum along x of the creation operator of the A site
// ky1 = momentum along y of the creation operator of the A site
// kx2 = momentum along x of the creation operator of the B site
// ky2 = momentum along y of the creation operator of the B site
// kx3 = momentum along x of the creation operator of the C site
// ky3 = momentum along y of the creation operator of the C site
// return value = corresponding matrix element
inline Complex ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementABCUpIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 0.0;
}



// compute the matrix element for the creation part of the three body interaction between sites A, B and C, for triangle pointing down
//
// kx1 = momentum along x of the creation operator of the A site
// ky1 = momentum along y of the creation operator of the A site
// kx2 = momentum along x of the creation operator of the B site
// ky2 = momentum along y of the creation operator of the B site
// kx3 = momentum along x of the creation operator of the C site
// ky3 = momentum along y of the creation operator of the C site
// return value = corresponding matrix element
inline Complex ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementABCDownIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 0.0;
}


// compute the matrix element for the annihilation part of the three body interaction between sites A, B and C, for triangle standing up 
//
// kx1 = momentum along x of the creation operator of the A site
// ky1 = momentum along y of the creation operator of the A site
// kx2 = momentum along x of the creation operator of the B site
// ky2 = momentum along y of the creation operator of the B site
// kx3 = momentum along x of the creation operator of the C site
// ky3 = momentum along y of the creation operator of the C site
// return value = corresponding matrix element
inline Complex ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementABCUpOut(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 0.0;
}



// compute the matrix element for the annihilation part of the three body interaction between sites A, B and C, for triangle pointing down
//
// kx1 = momentum along x of the creation operator of the A site
// ky1 = momentum along y of the creation operator of the A site
// kx2 = momentum along x of the creation operator of the B site
// ky2 = momentum along y of the creation operator of the B site
// kx3 = momentum along x of the creation operator of the C site
// ky3 = momentum along y of the creation operator of the C site
// return value = corresponding matrix element
inline Complex ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementABCDownOut(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 0.0;
}



// compute the matrix element for the creation part of the three body on site interaction for the A1 sites 
//
// kx1 = momentum along x of the first creation operator
// ky1 = momentum along y of the first creation operator
// kx2 = momentum along x of the second creation operator
// ky2 = momentum along y of the second creation operator
// kx3 = momentum along x of the third creation operator
// ky3 = momentum along y of the third creation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteAAAIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body on site interaction for the A1 sites 
//
// kx4 = momentum along x of the first annihilation operator
// ky4 = momentum along y of the first annihilation operator
// kx5 = momentum along x of the second annihilation operator
// ky5 = momentum along y of the secondannihilation operator
// kx6 = momentum along x of the third annihilation operator
// ky6 = momentum along y of the third annihilation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteAAAOut(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the creation part of the three body on site interaction for the A2 sites 
//
// kx1 = momentum along x of the first creation operator
// ky1 = momentum along y of the first creation operator
// kx2 = momentum along x of the second creation operator
// ky2 = momentum along y of the second creation operator
// kx3 = momentum along x of the third creation operator
// ky3 = momentum along y of the third creation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteBBBIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body on site interaction for the A2 sites 
//
// kx4 = momentum along x of the first annihilation operator
// ky4 = momentum along y of the first annihilation operator
// kx5 = momentum along x of the second annihilation operator
// ky5 = momentum along y of the secondannihilation operator
// kx6 = momentum along x of the third annihilation operator
// ky6 = momentum along y of the third annihilation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteBBBOut(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  //return Phase(-0.5*this->KxFactor * ((double) (kx4+kx5+kx6)));
  return 1.0;
}

// compute the matrix element for the creation part of the three body on site interaction for the A3 sites 
//
// kx1 = momentum along x of the first creation operator
// ky1 = momentum along y of the first creation operator
// kx2 = momentum along x of the second creation operator
// ky2 = momentum along y of the second creation operator
// kx3 = momentum along x of the third creation operator
// ky3 = momentum along y of the third creation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteCCCIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  //return Phase(0.5*this->KyFactor * ((double) (ky1+ky2+ky3)));
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body on site interaction for the A3 sites 
//
// kx4 = momentum along x of the first annihilation operator
// ky4 = momentum along y of the first annihilation operator
// kx5 = momentum along x of the second annihilation operator
// ky5 = momentum along y of the secondannihilation operator
// kx6 = momentum along x of the third annihilation operator
// ky6 = momentum along y of the third annihilation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeChern2TriangularLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteCCCOut(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  //return Phase(-0.5*this->KyFactor * ((double) (ky4+ky5+ky6)));
  return 1.0;
}


#endif
