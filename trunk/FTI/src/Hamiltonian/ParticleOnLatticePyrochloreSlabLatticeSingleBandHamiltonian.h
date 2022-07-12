////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//    class of pyrochlore slab lattice model with interacting particles       //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 17/05/2012                      //
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


#ifndef PARTICLEONLATTICEPYROCHLORESLABLATTICESINGLEBANDHAMILTONIAN_H
#define PARTICLEONLATTICEPYROCHLORESLABLATTICESINGLEBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticePyrochloreSlabLatticeSingleBandHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{

 protected:
  
  
  // strength of the repulsive two body onsite interaction
  double UPotential;
  
  // strength of the repulsive two body neareast neighbor interaction
  double VPotential;

  // index of the band to be filled
  int BandIndex;

  // use flat band model
  bool FlatBand;
  
 public:

  // constructor
  //
  ParticleOnLatticePyrochloreSlabLatticeSingleBandHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body onsite interaction
  // vPotential = strength of the repulsive two body neareast neighbor interaction
  // tightBindingModel = pointer to the tight binding model
  // bandIndex = index of the band to consider
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticePyrochloreSlabLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, 
							      int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, 
							      Abstract2DTightBindingModel* tightBindingModel, int bandIndex,
							      bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticePyrochloreSlabLatticeSingleBandHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two sites A and B 
  //
  // kx1 = creation momentum along x for the B site
  // ky1 = creation momentum along y for the B site
  // kx2 = annihilation momentum along x for the B site
  // ky2 = annihilation momentum along y for the B site
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

};



#endif
