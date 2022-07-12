////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//          class of kagome lattice model with interacting particles          //
//                                                                            //
//                        last modification : 09/12/2011                      //
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


#ifndef PARTICLEONLATTICEWITHSPINKAGOMELATTICETHREEBANDHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINKAGOMELATTICETHREEBANDHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallThreeBandHamiltonian.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeKagomeLatticeThreeBandHamiltonian : public ParticleOnLatticeQuantumSpinHallThreeBandHamiltonian
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
  
  // four times the sublattice staggered chemical potential 
  double MuS;
  // nearest neighbor density-density potential strength
  double UPotential;

  // boundary condition twisting angle along x
  double GammaX;
  // boundary condition twisting angle along y
  double GammaY;

  // use flat band model
  bool FlatBand;

  // numerical factor for momentum along x
  double KxFactor;
  // numerical factor for momentum along y
  double KyFactor;
  

 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // t1 = real part of the hopping amplitude between neareast neighbor sites
  // t2 = real part of the hopping amplitude between next neareast neighbor sites
  // lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
  // lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
  // mus = sublattice chemical potential on A sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeKagomeLatticeThreeBandHamiltonian(ParticleOnSphereWithSU3Spin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double t1, double t2, double lambda1, double lambda2, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeKagomeLatticeThreeBandHamiltonian();
  
 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two sites A and B 
  //
  // kx1 = creation momentum along x for the B site
  // ky1 = creation momentum along y for the B site
  // k2a = annihilation momentum along x for the B site
  // k2b = annihilation momentum along y for the B site
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
  // kx1 = creation momentum along x for the B site
  // ky1 = creation momentum along y for the B site
  // kx2 = creation momentum along x for the C site
  // ky2 = creation momentum along y for the C site
  // kx3 = annihilation momentum along x for the B site
  // ky3 = annihilation momentum along y for the B site
  // kx4 = annihilation momentum along x for the C site
  // ky4 = annihilation momentum along y for the C site
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementBC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

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

  // compute the one body hamiltonians related to the band stucture contribution
  //
  // oneBodyHamiltonians = array of one body hamiltonians
  virtual void ComputeOneBodyHamiltonian(HermitianMatrix* oneBodyHamiltonians);


};

#endif
