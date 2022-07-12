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


#ifndef PARTICLEONLATTICEWITHSPINRUBYLATTICESINGLEBANDTHREEBODYHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINRUBYLATTICESINGLEBANDTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian
{

 protected:
  
  // nearest neighbor density-density-density potential strength
  double UPotential;
  // nearest neighbor density-density potential strength
  double VPotential;


  // use flat band model
  bool FlatBand;
  

 public:

  // default constructor
  //
  ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian();

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
  ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, 
							     Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two sites A and B 
  //
  // kx1 = momentum along x for the A site
  // ky1 = momentum along y for the A site
  // kx2 = momentum along x for the B site
  // ky2 = momentum along y for the B site
  // return value = corresponding matrix element
  virtual  Complex ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2);

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
  virtual  Complex ComputeThreeBodyMatrixElementBAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

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
  virtual  Complex ComputeThreeBodyMatrixElementABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

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
  virtual  Complex ComputeThreeBodyMatrixElementOnSiteAAB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

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
  virtual Complex ComputeThreeBodyMatrixElementOnSiteABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the two body interaction between two A sites (or two B sites) 
  //
  // kx1 = momentum along x for the first A site
  // ky1 = momentum along y for the first A site
  // kx2 = momentum along x for the second A site
  // ky2 = momentum along y for the second A site
  // return value = corresponding matrix element
  //  virtual Complex ComputeTwoBodyMatrixElementAA(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the creation part of the three body on site interaction for the A1 sites 
  //
  // kx1 = momentum along x of the first creation operator
  // ky1 = momentum along y of the first creation operator
  // kx2 = momentum along x of the second creation operator
  // ky2 = momentum along y of the second creation operator
  // kx3 = momentum along x of the third creation operator
  // ky3 = momentum along y of the third creation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteA1A1A1In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body on site interaction for the A1 sites 
  //
  // kx4 = momentum along x of the first annihilation operator
  // ky4 = momentum along y of the first annihilation operator
  // kx5 = momentum along x of the second annihilation operator
  // ky5 = momentum along y of the secondannihilation operator
  // kx6 = momentum along x of the third annihilation operator
  // ky6 = momentum along y of the third annihilation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteA1A1A1Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body on site interaction for the A2 sites 
  //
  // kx1 = momentum along x of the first creation operator
  // ky1 = momentum along y of the first creation operator
  // kx2 = momentum along x of the second creation operator
  // ky2 = momentum along y of the second creation operator
  // kx3 = momentum along x of the third creation operator
  // ky3 = momentum along y of the third creation operator
  // return value = corresponding matrix element
  virtual  Complex ComputeThreeBodyMatrixElementOnSiteA2A2A2In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body on site interaction for the A2 sites 
  //
  // kx4 = momentum along x of the first annihilation operator
  // ky4 = momentum along y of the first annihilation operator
  // kx5 = momentum along x of the second annihilation operator
  // ky5 = momentum along y of the secondannihilation operator
  // kx6 = momentum along x of the third annihilation operator
  // ky6 = momentum along y of the third annihilation operator
  // return value = corresponding matrix element
  virtual  Complex ComputeThreeBodyMatrixElementOnSiteA2A2A2Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body on site interaction for the A3 sites 
  //
  // kx1 = momentum along x of the first creation operator
  // ky1 = momentum along y of the first creation operator
  // kx2 = momentum along x of the second creation operator
  // ky2 = momentum along y of the second creation operator
  // kx3 = momentum along x of the third creation operator
  // ky3 = momentum along y of the third creation operator
  // return value = corresponding matrix element
  virtual  Complex ComputeThreeBodyMatrixElementOnSiteA3A3A3In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body on site interaction for the A3 sites 
  //
  // kx4 = momentum along x of the first annihilation operator
  // ky4 = momentum along y of the first annihilation operator
  // kx5 = momentum along x of the second annihilation operator
  // ky5 = momentum along y of the secondannihilation operator
  // kx6 = momentum along x of the third annihilation operator
  // ky6 = momentum along y of the third annihilation operator
  // return value = corresponding matrix element
  virtual  Complex ComputeThreeBodyMatrixElementOnSiteA3A3A3Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body on site interaction for the A4 sites 
  //
  // kx1 = momentum along x of the first creation operator
  // ky1 = momentum along y of the first creation operator
  // kx2 = momentum along x of the second creation operator
  // ky2 = momentum along y of the second creation operator
  // kx3 = momentum along x of the third creation operator
  // ky3 = momentum along y of the third creation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteA4A4A4In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body on site interaction for the A4 sites 
  //
  // kx4 = momentum along x of the first annihilation operator
  // ky4 = momentum along y of the first annihilation operator
  // kx5 = momentum along x of the second annihilation operator
  // ky5 = momentum along y of the secondannihilation operator
  // kx6 = momentum along x of the third annihilation operator
  // ky6 = momentum along y of the third annihilation operator
  // return value = corresponding matrix element
  virtual  Complex ComputeThreeBodyMatrixElementOnSiteA4A4A4Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body on site interaction for the A5 sites 
  //
  // kx1 = momentum along x of the first creation operator
  // ky1 = momentum along y of the first creation operator
  // kx2 = momentum along x of the second creation operator
  // ky2 = momentum along y of the second creation operator
  // kx3 = momentum along x of the third creation operator
  // ky3 = momentum along y of the third creation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteA5A5A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body on site interaction for the A5 sites 
  //
  // kx4 = momentum along x of the first annihilation operator
  // ky4 = momentum along y of the first annihilation operator
  // kx5 = momentum along x of the second annihilation operator
  // ky5 = momentum along y of the secondannihilation operator
  // kx6 = momentum along x of the third annihilation operator
  // ky6 = momentum along y of the third annihilation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteA5A5A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body on site interaction for the A6 sites 
  //
  // kx1 = momentum along x of the first creation operator
  // ky1 = momentum along y of the first creation operator
  // kx2 = momentum along x of the second creation operator
  // ky2 = momentum along y of the second creation operator
  // kx3 = momentum along x of the third creation operator
  // ky3 = momentum along y of the third creation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteA6A6A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body on site interaction for the A6 sites 
  //
  // kx4 = momentum along x of the first annihilation operator
  // ky4 = momentum along y of the first annihilation operator
  // kx5 = momentum along x of the second annihilation operator
  // ky5 = momentum along y of the secondannihilation operator
  // kx6 = momentum along x of the third annihilation operator
  // ky6 = momentum along y of the third annihilation operator
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementOnSiteA6A6A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A1, A3 and A5
  //
  // kx1 = momentum along x of the creation operator of the A1 site
  // ky1 = momentum along y of the creation operator of the A1 site
  // kx2 = momentum along x of the creation operator of the A3 site
  // ky2 = momentum along y of the creation operator of the A3 site
  // kx3 = momentum along x of the creation operator of the A5 site
  // ky3 = momentum along y of the creation operator of the A5 site
  // return value = corresponding matrix element
   virtual Complex ComputeThreeBodyMatrixElementA1A3A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A1, A3 and A5 
  //
  // kx4 = momentum along x of the annihilation operator of the A1 site
  // ky4 = momentum along y of the annihilation operator of the A1 site
  // kx5 = momentum along x of the annihilation operator of the A3 site
  // ky5 = momentum along y of the annihilation operator of the A3 site
  // kx6 = momentum along x of the annihilation operator of the A5 site
  // ky6 = momentum along y of the annihilation operator of the A5 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A3A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A2, A4 and A6
  //
  // kx1 = momentum along x of the creation operator of the A2 site
  // ky1 = momentum along y of the creation operator of the A2 site
  // kx2 = momentum along x of the creation operator of the A4 site
  // ky2 = momentum along y of the creation operator of the A4 site
  // kx3 = momentum along x of the creation operator of the A6 site
  // ky3 = momentum along y of the creation operator of the A6 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA2A4A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A2, A4 and A6 
  //
  // kx4 = momentum along x of the annihilation operator of the A2 site
  // ky4 = momentum along y of the annihilation operator of the A2 site
  // kx5 = momentum along x of the annihilation operator of the A4 site
  // ky5 = momentum along y of the annihilation operator of the A4 site
  // kx6 = momentum along x of the annihilation operator of the A6 site
  // ky6 = momentum along y of the annihilation operator of the A6 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA2A4A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A1, A2 and A5
  //
  // kx1 = momentum along x of the creation operator of the A1 site
  // ky1 = momentum along y of the creation operator of the A1 site
  // kx2 = momentum along x of the creation operator of the A2 site
  // ky2 = momentum along y of the creation operator of the A2 site
  // kx3 = momentum along x of the creation operator of the A5 site
  // ky3 = momentum along y of the creation operator of the A5 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A2A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A1, A2 and A5 
  //
  // kx4 = momentum along x of the annihilation operator of the A1 site
  // ky4 = momentum along y of the annihilation operator of the A1 site
  // kx5 = momentum along x of the annihilation operator of the A2 site
  // ky5 = momentum along y of the annihilation operator of the A2 site
  // kx6 = momentum along x of the annihilation operator of the A5 site
  // ky6 = momentum along y of the annihilation operator of the A5 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A2A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A1, A4 and A5
  //
  // kx1 = momentum along x of the creation operator of the A1 site
  // ky1 = momentum along y of the creation operator of the A1 site
  // kx2 = momentum along x of the creation operator of the A4 site
  // ky2 = momentum along y of the creation operator of the A4 site
  // kx3 = momentum along x of the creation operator of the A5 site
  // ky3 = momentum along y of the creation operator of the A5 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A4A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A1, A4 and A5 
  //
  // kx4 = momentum along x of the annihilation operator of the A1 site
  // ky4 = momentum along y of the annihilation operator of the A1 site
  // kx5 = momentum along x of the annihilation operator of the A4 site
  // ky5 = momentum along y of the annihilation operator of the A4 site
  // kx6 = momentum along x of the annihilation operator of the A5 site
  // ky6 = momentum along y of the annihilation operator of the A5 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A4A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A2, A4 and A5
  //
  // kx1 = momentum along x of the creation operator of the A2 site
  // ky1 = momentum along y of the creation operator of the A2 site
  // kx2 = momentum along x of the creation operator of the A4 site
  // ky2 = momentum along y of the creation operator of the A4 site
  // kx3 = momentum along x of the creation operator of the A5 site
  // ky3 = momentum along y of the creation operator of the A5 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA2A4A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A2, A4 and A5 
  //
  // kx4 = momentum along x of the annihilation operator of the A2 site
  // ky4 = momentum along y of the annihilation operator of the A2 site
  // kx5 = momentum along x of the annihilation operator of the A4 site
  // ky5 = momentum along y of the annihilation operator of the A4 site
  // kx6 = momentum along x of the annihilation operator of the A5 site
  // ky6 = momentum along y of the annihilation operator of the A5 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA2A4A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A1, A2 and A4
  //
  // kx1 = momentum along x of the creation operator of the A1 site
  // ky1 = momentum along y of the creation operator of the A1 site
  // kx2 = momentum along x of the creation operator of the A2 site
  // ky2 = momentum along y of the creation operator of the A2 site
  // kx3 = momentum along x of the creation operator of the A4 site
  // ky3 = momentum along y of the creation operator of the A4 site
  // return value = corresponding matrix element
   virtual Complex ComputeThreeBodyMatrixElementA1A2A4In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A1, A2 and A4 
  //
  // kx4 = momentum along x of the annihilation operator of the A1 site
  // ky4 = momentum along y of the annihilation operator of the A1 site
  // kx5 = momentum along x of the annihilation operator of the A2 site
  // ky5 = momentum along y of the annihilation operator of the A2 site
  // kx6 = momentum along x of the annihilation operator of the A4 site
  // ky6 = momentum along y of the annihilation operator of the A4 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A2A4Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A1, A2 and A6
  //
  // kx1 = momentum along x of the creation operator of the A1 site
  // ky1 = momentum along y of the creation operator of the A1 site
  // kx2 = momentum along x of the creation operator of the A2 site
  // ky2 = momentum along y of the creation operator of the A2 site
  // kx3 = momentum along x of the creation operator of the A6 site
  // ky3 = momentum along y of the creation operator of the A6 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A3A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A1, A2 and A6 
  //
  // kx4 = momentum along x of the annihilation operator of the A1 site
  // ky4 = momentum along y of the annihilation operator of the A1 site
  // kx5 = momentum along x of the annihilation operator of the A2 site
  // ky5 = momentum along y of the annihilation operator of the A2 site
  // kx6 = momentum along x of the annihilation operator of the A6 site
  // ky6 = momentum along y of the annihilation operator of the A6 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A3A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A1, A3 and A4
  //
  // kx1 = momentum along x of the creation operator of the A1 site
  // ky1 = momentum along y of the creation operator of the A1 site
  // kx2 = momentum along x of the creation operator of the A3 site
  // ky2 = momentum along y of the creation operator of the A3 site
  // kx3 = momentum along x of the creation operator of the A4 site
  // ky3 = momentum along y of the creation operator of the A4 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A3A4In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A1, A3 and A4 
  //
  // kx4 = momentum along x of the annihilation operator of the A1 site
  // ky4 = momentum along y of the annihilation operator of the A1 site
  // kx5 = momentum along x of the annihilation operator of the A3 site
  // ky5 = momentum along y of the annihilation operator of the A3 site
  // kx6 = momentum along x of the annihilation operator of the A4 site
  // ky6 = momentum along y of the annihilation operator of the A4 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A3A4Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A1, A4 and A6
  //
  // kx1 = momentum along x of the creation operator of the A1 site
  // ky1 = momentum along y of the creation operator of the A1 site
  // kx2 = momentum along x of the creation operator of the A4 site
  // ky2 = momentum along y of the creation operator of the A4 site
  // kx3 = momentum along x of the creation operator of the A6 site
  // ky3 = momentum along y of the creation operator of the A6 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A4A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A1, A4 and A6 
  //
  // kx4 = momentum along x of the annihilation operator of the A1 site
  // ky4 = momentum along y of the annihilation operator of the A1 site
  // kx5 = momentum along x of the annihilation operator of the A4 site
  // ky5 = momentum along y of the annihilation operator of the A4 site
  // kx6 = momentum along x of the annihilation operator of the A6 site
  // ky6 = momentum along y of the annihilation operator of the A6 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA1A4A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A3, A4 and A6
  //
  // kx1 = momentum along x of the creation operator of the A3 site
  // ky1 = momentum along y of the creation operator of the A3 site
  // kx2 = momentum along x of the creation operator of the A4 site
  // ky2 = momentum along y of the creation operator of the A4 site
  // kx3 = momentum along x of the creation operator of the A6 site
  // ky3 = momentum along y of the creation operator of the A6 site
  // return value = corresponding matrix element
   virtual Complex ComputeThreeBodyMatrixElementA3A4A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A3, A4 and A6 
  //
  // kx4 = momentum along x of the annihilation operator of the A3 site
  // ky4 = momentum along y of the annihilation operator of the A3 site
  // kx5 = momentum along x of the annihilation operator of the A4 site
  // ky5 = momentum along y of the annihilation operator of the A4 site
  // kx6 = momentum along x of the annihilation operator of the A6 site
  // ky6 = momentum along y of the annihilation operator of the A6 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA3A4A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A2, A5 and A6
  //
  // kx1 = momentum along x of the creation operator of the A2 site
  // ky1 = momentum along y of the creation operator of the A2 site
  // kx2 = momentum along x of the creation operator of the A5 site
  // ky2 = momentum along y of the creation operator of the A5 site
  // kx3 = momentum along x of the creation operator of the A6 site
  // ky3 = momentum along y of the creation operator of the A6 site
  // return value = corresponding matrix element
   virtual Complex ComputeThreeBodyMatrixElementA2A5A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A2, A5 and A6 
  //
  // kx4 = momentum along x of the annihilation operator of the A2 site
  // ky4 = momentum along y of the annihilation operator of the A2 site
  // kx5 = momentum along x of the annihilation operator of the A5 site
  // ky5 = momentum along y of the annihilation operator of the A5 site
  // kx6 = momentum along x of the annihilation operator of the A6 site
  // ky6 = momentum along y of the annihilation operator of the A6 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA2A5A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A2, A3 and A6
  //
  // kx1 = momentum along x of the creation operator of the A2 site
  // ky1 = momentum along y of the creation operator of the A2 site
  // kx2 = momentum along x of the creation operator of the A3 site
  // ky2 = momentum along y of the creation operator of the A3 site
  // kx3 = momentum along x of the creation operator of the A6 site
  // ky3 = momentum along y of the creation operator of the A6 site
  // return value = corresponding matrix element
   virtual Complex ComputeThreeBodyMatrixElementA2A3A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A2, A3 and A6 
  //
  // kx4 = momentum along x of the annihilation operator of the A2 site
  // ky4 = momentum along y of the annihilation operator of the A2 site
  // kx5 = momentum along x of the annihilation operator of the A3 site
  // ky5 = momentum along y of the annihilation operator of the A3 site
  // kx6 = momentum along x of the annihilation operator of the A6 site
  // ky6 = momentum along y of the annihilation operator of the A6 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA2A3A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A2, A3 and A5
  //
  // kx1 = momentum along x of the creation operator of the A2 site
  // ky1 = momentum along y of the creation operator of the A2 site
  // kx2 = momentum along x of the creation operator of the A3 site
  // ky2 = momentum along y of the creation operator of the A3 site
  // kx3 = momentum along x of the creation operator of the A5 site
  // ky3 = momentum along y of the creation operator of the A5 site
  // return value = corresponding matrix element
   virtual Complex ComputeThreeBodyMatrixElementA2A3A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A2, A3 and A5 
  //
  // kx4 = momentum along x of the annihilation operator of the A2 site
  // ky4 = momentum along y of the annihilation operator of the A2 site
  // kx5 = momentum along x of the annihilation operator of the A3 site
  // ky5 = momentum along y of the annihilation operator of the A3 site
  // kx6 = momentum along x of the annihilation operator of the A5 site
  // ky6 = momentum along y of the annihilation operator of the A5 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA2A3A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the creation part of the three body interaction between sites A6, A3 and A5
  //
  // kx1 = momentum along x of the creation operator of the A6 site
  // ky1 = momentum along y of the creation operator of the A6 site
  // kx2 = momentum along x of the creation operator of the A3 site
  // ky2 = momentum along y of the creation operator of the A3 site
  // kx3 = momentum along x of the creation operator of the A5 site
  // ky3 = momentum along y of the creation operator of the A5 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA6A3A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3);

  // compute the matrix element for the annihilation part of the three body interaction between sites A6, A3 and A5 
  //
  // kx4 = momentum along x of the annihilation operator of the A6 site
  // ky4 = momentum along y of the annihilation operator of the A6 site
  // kx5 = momentum along x of the annihilation operator of the A3 site
  // ky5 = momentum along y of the annihilation operator of the A3 site
  // kx6 = momentum along x of the annihilation operator of the A5 site
  // ky6 = momentum along y of the annihilation operator of the A5 site
  // return value = corresponding matrix element
  virtual Complex ComputeThreeBodyMatrixElementA6A3A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

};

// compute the matrix element for the creation part of the three body interaction between sites A1, A3 and A5
//
// kx1 = momentum along x of the creation operator of the A1 site
// ky1 = momentum along y of the creation operator of the A1 site
// kx2 = momentum along x of the creation operator of the A3 site
// ky2 = momentum along y of the creation operator of the A3 site
// kx3 = momentum along x of the creation operator of the A5 site
// ky3 = momentum along y of the creation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A3A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body interaction between sites A1, A3 and A5 
//
// kx6 = momentum along x of the annihilation operator of the A1 site
// ky4 = momentum along y of the annihilation operator of the A1 site
// kx5 = momentum along x of the annihilation operator of the A3 site
// ky5 = momentum along y of the annihilation operator of the A3 site
// kx6 = momentum along x of the annihilation operator of the A5 site
// ky6 = momentum along y of the annihilation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A3A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the creation part of the three body interaction between sites A2, A4 and A6
//
// kx1 = momentum along x of the creation operator of the A2 site
// ky1 = momentum along y of the creation operator of the A2 site
// kx2 = momentum along x of the creation operator of the A4 site
// ky2 = momentum along y of the creation operator of the A4 site
// kx3 = momentum along x of the creation operator of the A6 site
// ky3 = momentum along y of the creation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA2A4A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body interaction between sites A2, A4 and A6 
//
// kx4 = momentum along x of the annihilation operator of the A2 site
// ky4 = momentum along y of the annihilation operator of the A2 site
// kx5 = momentum along x of the annihilation operator of the A4 site
// ky5 = momentum along y of the annihilation operator of the A4 site
// kx6 = momentum along x of the annihilation operator of the A6 site
// ky6 = momentum along y of the annihilation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA2A4A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the creation part of the three body interaction between sites A1, A2 and A5
//
// kx1 = momentum along x of the creation operator of the A1 site
// ky1 = momentum along y of the creation operator of the A1 site
// kx2 = momentum along x of the creation operator of the A2 site
// ky2 = momentum along y of the creation operator of the A2 site
// kx3 = momentum along x of the creation operator of the A5 site
// ky3 = momentum along y of the creation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A2A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body interaction between sites A1, A2 and A5 
//
// kx4 = momentum along x of the annihilation operator of the A1 site
// ky4 = momentum along y of the annihilation operator of the A1 site
// kx5 = momentum along x of the annihilation operator of the A2 site
// ky5 = momentum along y of the annihilation operator of the A2 site
// kx6 = momentum along x of the annihilation operator of the A5 site
// ky6 = momentum along y of the annihilation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A2A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the creation part of the three body interaction between sites A1, A4 and A5
//
// kx1 = momentum along x of the creation operator of the A1 site
// ky1 = momentum along y of the creation operator of the A1 site
// kx2 = momentum along x of the creation operator of the A4 site
// ky2 = momentum along y of the creation operator of the A4 site
// kx3 = momentum along x of the creation operator of the A5 site
// ky3 = momentum along y of the creation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A4A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body interaction between sites A1, A4 and A5 
//
// kx4 = momentum along x of the annihilation operator of the A1 site
// ky4 = momentum along y of the annihilation operator of the A1 site
// kx5 = momentum along x of the annihilation operator of the A4 site
// ky5 = momentum along y of the annihilation operator of the A4 site
// kx6 = momentum along x of the annihilation operator of the A5 site
// ky6 = momentum along y of the annihilation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A4A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the creation part of the three body interaction between sites A2, A4 and A5
//
// kx1 = momentum along x of the creation operator of the A2 site
// ky1 = momentum along y of the creation operator of the A2 site
// kx2 = momentum along x of the creation operator of the A4 site
// ky2 = momentum along y of the creation operator of the A4 site
// kx3 = momentum along x of the creation operator of the A5 site
// ky3 = momentum along y of the creation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA2A4A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body interaction between sites A2, A4 and A5 
//
// kx4 = momentum along x of the annihilation operator of the A2 site
// ky4 = momentum along y of the annihilation operator of the A2 site
// kx5 = momentum along x of the annihilation operator of the A4 site
// ky5 = momentum along y of the annihilation operator of the A4 site
// kx6 = momentum along x of the annihilation operator of the A5 site
// ky6 = momentum along y of the annihilation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA2A4A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the creation part of the three body interaction between sites A1, A2 and A4
//
// kx1 = momentum along x of the creation operator of the A1 site
// ky1 = momentum along y of the creation operator of the A1 site
// kx2 = momentum along x of the creation operator of the A2 site
// ky2 = momentum along y of the creation operator of the A2 site
// kx3 = momentum along x of the creation operator of the A4 site
// ky3 = momentum along y of the creation operator of the A4 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A2A4In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body interaction between sites A1, A2 and A4 
//
// kx4 = momentum along x of the annihilation operator of the A1 site
// ky4 = momentum along y of the annihilation operator of the A1 site
// kx5 = momentum along x of the annihilation operator of the A2 site
// ky5 = momentum along y of the annihilation operator of the A2 site
// kx6 = momentum along x of the annihilation operator of the A4 site
// ky6 = momentum along y of the annihilation operator of the A4 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A2A4Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the creation part of the three body interaction between sites A1, A2 and A6
//
// kx1 = momentum along x of the creation operator of the A1 site
// ky1 = momentum along y of the creation operator of the A1 site
// kx2 = momentum along x of the creation operator of the A2 site
// ky2 = momentum along y of the creation operator of the A2 site
// kx3 = momentum along x of the creation operator of the A6 site
// ky3 = momentum along y of the creation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A3A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return Phase(-this->KxFactor * ((double) kx3));
}

// compute the matrix element for the annihilation part of the three body interaction between sites A1, A2 and A6 
//
// kx4 = momentum along x of the annihilation operator of the A1 site
// ky4 = momentum along y of the annihilation operator of the A1 site
// kx5 = momentum along x of the annihilation operator of the A2 site
// ky5 = momentum along y of the annihilation operator of the A2 site
// kx6 = momentum along x of the annihilation operator of the A6 site
// ky6 = momentum along y of the annihilation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A3A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return Phase(this->KxFactor * ((double) kx6));
}

// compute the matrix element for the creation part of the three body interaction between sites A1, A3 and A4
//
// kx1 = momentum along x of the creation operator of the A1 site
// ky1 = momentum along y of the creation operator of the A1 site
// kx2 = momentum along x of the creation operator of the A3 site
// ky2 = momentum along y of the creation operator of the A3 site
// kx3 = momentum along x of the creation operator of the A4 site
// ky3 = momentum along y of the creation operator of the A4 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A3A4In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return Phase(-this->KxFactor * ((double) kx3));
}

// compute the matrix element for the annihilation part of the three body interaction between sites A1, A3 and A4 
//
// kx4 = momentum along x of the annihilation operator of the A1 site
// ky4 = momentum along y of the annihilation operator of the A1 site
// kx5 = momentum along x of the annihilation operator of the A3 site
// ky5 = momentum along y of the annihilation operator of the A3 site
// kx6 = momentum along x of the annihilation operator of the A4 site
// ky6 = momentum along y of the annihilation operator of the A4 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A3A4Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return Phase(this->KxFactor * ((double) kx6));
}

// compute the matrix element for the creation part of the three body interaction between sites A1, A4 and A6
//
// kx1 = momentum along x of the creation operator of the A1 site
// ky1 = momentum along y of the creation operator of the A1 site
// kx2 = momentum along x of the creation operator of the A4 site
// ky2 = momentum along y of the creation operator of the A4 site
// kx3 = momentum along x of the creation operator of the A6 site
// ky3 = momentum along y of the creation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A4A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return Phase(-this->KxFactor * ((double) (kx2 + kx3)));
}

// compute the matrix element for the annihilation part of the three body interaction between sites A1, A4 and A6 
//
// kx4 = momentum along x of the annihilation operator of the A1 site
// ky4 = momentum along y of the annihilation operator of the A1 site
// kx5 = momentum along x of the annihilation operator of the A4 site
// ky5 = momentum along y of the annihilation operator of the A4 site
// kx6 = momentum along x of the annihilation operator of the A6 site
// ky6 = momentum along y of the annihilation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA1A4A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return Phase(this->KxFactor * ((double) (kx5 + kx6)));
}

// compute the matrix element for the creation part of the three body interaction between sites A3, A4 and A6
//
// kx1 = momentum along x of the creation operator of the A3 site
// ky1 = momentum along y of the creation operator of the A3 site
// kx2 = momentum along x of the creation operator of the A4 site
// ky2 = momentum along y of the creation operator of the A4 site
// kx3 = momentum along x of the creation operator of the A6 site
// ky3 = momentum along y of the creation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA3A4A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return Phase(-this->KxFactor * ((double) (kx2 + kx3)));
}

// compute the matrix element for the annihilation part of the three body interaction between sites A3, A4 and A6 
//
// kx4 = momentum along x of the annihilation operator of the A3 site
// ky4 = momentum along y of the annihilation operator of the A3 site
// kx5 = momentum along x of the annihilation operator of the A4 site
// ky5 = momentum along y of the annihilation operator of the A4 site
// kx6 = momentum along x of the annihilation operator of the A6 site
// ky6 = momentum along y of the annihilation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA3A4A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return Phase(this->KxFactor * ((double) (kx5 + kx6)));
}

// compute the matrix element for the creation part of the three body interaction between sites A2, A5 and A6
//
// kx1 = momentum along x of the creation operator of the A2 site
// ky1 = momentum along y of the creation operator of the A2 site
// kx2 = momentum along x of the creation operator of the A5 site
// ky2 = momentum along y of the creation operator of the A5 site
// kx3 = momentum along x of the creation operator of the A6 site
// ky3 = momentum along y of the creation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA2A5A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return Phase(this->KxFactor * ((double) kx2) + this->KyFactor * ((double) ky2));
}

// compute the matrix element for the annihilation part of the three body interaction between sites A2, A5 and A6 
//
// kx4 = momentum along x of the annihilation operator of the A2 site
// ky4 = momentum along y of the annihilation operator of the A2 site
// kx5 = momentum along x of the annihilation operator of the A5 site
// ky5 = momentum along y of the annihilation operator of the A5 site
// kx6 = momentum along x of the annihilation operator of the A6 site
// ky6 = momentum along y of the annihilation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA2A5A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return Phase(-this->KxFactor * ((double) kx5) - this->KyFactor * ((double) ky5));
}

// compute the matrix element for the creation part of the three body interaction between sites A2, A3 and A6
//
// kx1 = momentum along x of the creation operator of the A2 site
// ky1 = momentum along y of the creation operator of the A2 site
// kx2 = momentum along x of the creation operator of the A3 site
// ky2 = momentum along y of the creation operator of the A3 site
// kx3 = momentum along x of the creation operator of the A6 site
// ky3 = momentum along y of the creation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA2A3A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return Phase(this->KxFactor * ((double) kx2) + this->KyFactor * ((double) ky2));
}

// compute the matrix element for the annihilation part of the three body interaction between sites A2, A3 and A6 
//
// kx4 = momentum along x of the annihilation operator of the A2 site
// ky4 = momentum along y of the annihilation operator of the A2 site
// kx5 = momentum along x of the annihilation operator of the A3 site
// ky5 = momentum along y of the annihilation operator of the A3 site
// kx6 = momentum along x of the annihilation operator of the A6 site
// ky6 = momentum along y of the annihilation operator of the A6 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA2A3A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return Phase(-this->KxFactor * ((double) kx5) - this->KyFactor * ((double) ky5));
}

// compute the matrix element for the creation part of the three body interaction between sites A2, A3 and A5
//
// kx1 = momentum along x of the creation operator of the A2 site
// ky1 = momentum along y of the creation operator of the A2 site
// kx2 = momentum along x of the creation operator of the A3 site
// ky2 = momentum along y of the creation operator of the A3 site
// kx3 = momentum along x of the creation operator of the A5 site
// ky3 = momentum along y of the creation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA2A3A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return Phase(this->KxFactor * ((double) (kx2 + kx3)) + this->KyFactor * ((double) (ky2 + ky3)));
}

// compute the matrix element for the annihilation part of the three body interaction between sites A2, A3 and A5 
//
// kx4 = momentum along x of the annihilation operator of the A2 site
// ky4 = momentum along y of the annihilation operator of the A2 site
// kx5 = momentum along x of the annihilation operator of the A3 site
// ky5 = momentum along y of the annihilation operator of the A3 site
// kx6 = momentum along x of the annihilation operator of the A5 site
// ky6 = momentum along y of the annihilation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA2A3A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return Phase(-this->KxFactor * ((double) (kx5 + kx6)) - this->KyFactor * ((double) (ky5 + ky6)));
}

// compute the matrix element for the creation part of the three body interaction between sites A6, A3 and A5
//
// kx1 = momentum along x of the creation operator of the A6 site
// ky1 = momentum along y of the creation operator of the A6 site
// kx2 = momentum along x of the creation operator of the A3 site
// ky2 = momentum along y of the creation operator of the A3 site
// kx3 = momentum along x of the creation operator of the A5 site
// ky3 = momentum along y of the creation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA6A3A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return Phase(this->KxFactor * ((double) (kx2 + kx3)) + this->KyFactor * ((double) (ky2 + ky3)));
}

// compute the matrix element for the annihilation part of the three body interaction between sites A6, A3 and A5 
//
// kx4 = momentum along x of the annihilation operator of the A6 site
// ky4 = momentum along y of the annihilation operator of the A6 site
// kx5 = momentum along x of the annihilation operator of the A3 site
// ky5 = momentum along y of the annihilation operator of the A3 site
// kx6 = momentum along x of the annihilation operator of the A5 site
// ky6 = momentum along y of the annihilation operator of the A5 site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementA6A3A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return Phase(-this->KxFactor * ((double) (kx5 + kx6)) - this->KyFactor * ((double) (ky5 + ky6)));
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

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA1A1A1In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
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

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA1A1A1Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
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

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA2A2A2In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
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

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA2A2A2Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
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

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA3A3A3In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
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

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA3A3A3Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the creation part of the three body on site interaction for the A4 sites 
//
// kx1 = momentum along x of the first creation operator
// ky1 = momentum along y of the first creation operator
// kx2 = momentum along x of the second creation operator
// ky2 = momentum along y of the second creation operator
// kx3 = momentum along x of the third creation operator
// ky3 = momentum along y of the third creation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA4A4A4In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body on site interaction for the A4 sites 
//
// kx4 = momentum along x of the first annihilation operator
// ky4 = momentum along y of the first annihilation operator
// kx5 = momentum along x of the second annihilation operator
// ky5 = momentum along y of the secondannihilation operator
// kx6 = momentum along x of the third annihilation operator
// ky6 = momentum along y of the third annihilation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA4A4A4Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the creation part of the three body on site interaction for the A5 sites 
//
// kx1 = momentum along x of the first creation operator
// ky1 = momentum along y of the first creation operator
// kx2 = momentum along x of the second creation operator
// ky2 = momentum along y of the second creation operator
// kx3 = momentum along x of the third creation operator
// ky3 = momentum along y of the third creation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA5A5A5In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body on site interaction for the A5 sites 
//
// kx4 = momentum along x of the first annihilation operator
// ky4 = momentum along y of the first annihilation operator
// kx5 = momentum along x of the second annihilation operator
// ky5 = momentum along y of the secondannihilation operator
// kx6 = momentum along x of the third annihilation operator
// ky6 = momentum along y of the third annihilation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA5A5A5Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the creation part of the three body on site interaction for the A6 sites 
//
// kx1 = momentum along x of the first creation operator
// ky1 = momentum along y of the first creation operator
// kx2 = momentum along x of the second creation operator
// ky2 = momentum along y of the second creation operator
// kx3 = momentum along x of the third creation operator
// ky3 = momentum along y of the third creation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA6A6A6In(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3)
{
  return 1.0;
}

// compute the matrix element for the annihilation part of the three body on site interaction for the A6 sites 
//
// kx4 = momentum along x of the first annihilation operator
// ky4 = momentum along y of the first annihilation operator
// kx5 = momentum along x of the second annihilation operator
// ky5 = momentum along y of the secondannihilation operator
// kx6 = momentum along x of the third annihilation operator
// ky6 = momentum along y of the third annihilation operator
// return value = corresponding matrix element

inline Complex ParticleOnLatticeRubyLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteA6A6A6Out(int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}



#endif
