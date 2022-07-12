////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 3d topological insulator based on the Fu-Kane-Mele     //
//                      model and restricted to four bands                    //
//                                                                            //
//                        last modification : 26/08/2011                      //
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


#ifndef PARTICLEONCUBICLATTICEFOURBANDFUKANEMELECHECKERBOARDHAMILTONIAN_H
#define PARTICLEONCUBICLATTICEFOURBANDFUKANEMELECHECKERBOARDHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnCubicLatticeFourBandSimpleTIHamiltonian.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian : public ParticleOnCubicLatticeFourBandSimpleTIHamiltonian
{

 protected:
 
  // global energy scale of the kinetic energy term (i.e t1 hopping term)
  double KineticScale;
  // distortion of nearest neighbor hoping amplitude in the (111) direction
  double NNHopingDistortion111;
  // amplitude of the spin orbit coupling
  double SpinOrbitCoupling;
  // hoping amplitude between neareast neighbor sites
  double NNHoping;
  // hoping amplitude between next neareast neighbor sites
  double NextNNHoping;
  // hoping amplitude between second next neareast neighbor sites
  double SecondNextNNHoping;

 public:

  // default constructor
  //
  ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // vPotential = strength of the repulsive two body on site interaction
  // kineticScale = global energy scale of the kinetic energy term (i.e t1 hopping term)
  // nnHopingDistortion111 = distortion of nearest neighbor hoping amplitude in the (111) direction
  // spinOrbitCoupling = amplitude of the spin orbit coupling
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // gammaZ = boundary condition twisting angle along z
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, int nbrSiteZ, double uPotential, double vPotential, double kineticScale, double nnHopingDistortion111, double spinOrbitCoupling, double gammaX, double gammaY, double gammaZ, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian();
  

 protected:
 
  // compute the matrix element for the two body on site interaction for site A and up spins
  //
  // kx1 = momentum along x for the first creation operator on A site with spin up
  // ky1 = momentum along y for the first creation operator on A site with spin up
  // kz1 = momentum along z for the first creation operator on A site with spin up
  // kx2 = momentum along x for the creation operator on A site with spin up
  // ky2 = momentum along y for the creation operator on A site with spin up
  // kz2 = momentum along z for the creation operator on A site with spin up
  // kx3 = momentum along x for the first annihilation operator on A site with spin up
  // ky3 = momentum along y for the first annihilation operator on A site with spin up
  // kz3 = momentum along z for the first annihilation operator on A site with spin up
  // kx4 = momentum along x for the creation annihilation operator on B site with spin up
  // ky4 = momentum along y for the creation annihilation operator on B site with spin up
  // kz4 = momentum along z for the creation annihilation operator on B site with spin up
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAUpAUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4);


  // compute the matrix element for the two body on site interaction for site A and down spins
  //
  // kx1 = momentum along x for the first creation operator on A site with spin down
  // ky1 = momentum along y for the first creation operator on A site with spin down
  // kz1 = momentum along z for the first creation operator on A site with spin down
  // kx2 = momentum along x for the creation operator on A site with spin down
  // ky2 = momentum along y for the creation operator on A site with spin down
  // kz2 = momentum along z for the creation operator on A site with spin down
  // kx3 = momentum along x for the first annihilation operator on A site with spin down
  // ky3 = momentum along y for the first annihilation operator on A site with spin down
  // kz3 = momentum along z for the first annihilation operator on A site with spin down
  // kx4 = momentum along x for the creation annihilation operator on B site with spin down
  // ky4 = momentum along y for the creation annihilation operator on B site with spin down
  // kz4 = momentum along z for the creation annihilation operator on B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementADownADown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4);

  // compute the matrix element for the two body on site interaction for site B and up spins
  //
  // kx1 = momentum along x for the first creation operator on B site with spin up
  // ky1 = momentum along y for the first creation operator on B site with spin up
  // kz1 = momentum along z for the first creation operator on B site with spin up
  // kx2 = momentum along x for the creation operator on B site with spin up
  // ky2 = momentum along y for the creation operator on B site with spin up
  // kz2 = momentum along z for the creation operator on B site with spin up
  // kx3 = momentum along x for the first annihilation operator on B site with spin up
  // ky3 = momentum along y for the first annihilation operator on B site with spin up
  // kz3 = momentum along z for the first annihilation operator on B site with spin up
  // kx4 = momentum along x for the creation annihilation operator on B site with spin up
  // ky4 = momentum along y for the creation annihilation operator on B site with spin up
  // kz4 = momentum along z for the creation annihilation operator on B site with spin up
  // return value = corresponding matrix element  
  virtual Complex ComputeTwoBodyMatrixElementBUpBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4);

  // compute the matrix element for the two body on site interaction for site B and down spins
  //
  // kx1 = momentum along x for the first creation operator on B site with spin down
  // ky1 = momentum along y for the first creation operator on B site with spin down
  // kz1 = momentum along z for the first creation operator on B site with spin down
  // kx2 = momentum along x for the creation operator on B site with spin down
  // ky2 = momentum along y for the creation operator on B site with spin down
  // kz2 = momentum along z for the creation operator on B site with spin down
  // kx3 = momentum along x for the first annihilation operator on B site with spin down
  // ky3 = momentum along y for the first annihilation operator on B site with spin down
  // kz3 = momentum along z for the first annihilation operator on B site with spin down
  // kx4 = momentum along x for the creation annihilation operator on B site with spin down
  // ky4 = momentum along y for the creation annihilation operator on B site with spin down
  // kz4 = momentum along z for the creation annihilation operator on B site with spin down
  // return value = corresponding matrix element  
  virtual Complex ComputeTwoBodyMatrixElementBDownBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4);

  // compute the matrix element for the two body interaction between two sites A and B  belonging to the same layer
  //
  // kx1 = momentum along x for the creation operator on A site with spin up
  // ky1 = momentum along y for the creation operator on A site with spin up
  // kz1 = momentum along z for the creation operator on A site with spin up
  // kx2 = momentum along x for the creation operator on B site with spin up
  // ky2 = momentum along y for the creation operator on B site with spin up
  // kz2 = momentum along z for the creation operator on B site with spin up
  // kx3 = momentum along x for the annihilation operator on A site with spin up
  // ky3 = momentum along y for the annihilation operator on A site with spin up
  // kz3 = momentum along z for the annihilation operator on A site with spin up
  // kx4 = momentum along x for the annihilation operator on B site with spin up
  // ky4 = momentum along y for the annihilation operator on B site with spin up
  // kz4 = momentum along z for the annihilation operator on B site with spin up
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4);

  // compute the matrix element for the two body interaction between two sites A and B with down spins
  //
  // kx1 = momentum along x for the creation operator on A site with spin down
  // ky1 = momentum along y for the creation operator on A site with spin down
  // kz1 = momentum along z for the creation operator on A site with spin down
  // kx2 = momentum along x for the creation operator on B site with spin down
  // ky2 = momentum along y for the creation operator on B site with spin down
  // kz2 = momentum along z for the creation operator on B site with spin down
  // kx3 = momentum along x for the annihilation operator on A site with spin down
  // ky3 = momentum along y for the annihilation operator on A site with spin down
  // kz3 = momentum along z for the annihilation operator on A site with spin down
  // kx4 = momentum along x for the annihilation operator on B site with spin down
  // ky4 = momentum along y for the annihilation operator on B site with spin down
  // kz4 = momentum along z for the annihilation operator on B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4);
  
  // compute the matrix element for the two body interaction between two sites A and B with opposite spins
  //
  // kx1 = momentum along x for the creation operator on A site with spin down
  // ky1 = momentum along y for the creation operator on A site with spin down
  // kz1 = momentum along z for the creation operator on A site with spin down
  // kx2 = momentum along x for the creation operator on B site with spin up
  // ky2 = momentum along y for the creation operator on B site with spin up
  // kz2 = momentum along z for the creation operator on B site with spin up
  // kx3 = momentum along x for the annihilation operator on A site with spin down
  // ky3 = momentum along y for the annihilation operator on A site with spin down
  // kz3 = momentum along z for the annihilation operator on A site with spin down
  // kx4 = momentum along x for the annihilation operator on B site with spin up
  // ky4 = momentum along y for the annihilation operator on B site with spin up
  // kz4 = momentum along z for the annihilation operator on B site with spin up
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementADownBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4);

  // compute the matrix element for the two body interaction between two sites A and B with opposite spins
  //
  // kx1 = momentum along x for the creation operator on A site with spin up
  // ky1 = momentum along y for the creation operator on A site with spin up
  // kz1 = momentum along z for the creation operator on A site with spin up
  // kx2 = momentum along x for the creation operator on B site with spin down
  // ky2 = momentum along y for the creation operator on B site with spin down
  // kz2 = momentum along z for the creation operator on B site with spin down
  // kx3 = momentum along x for the annihilation operator on A site with spin up
  // ky3 = momentum along y for the annihilation operator on A site with spin up
  // kz3 = momentum along z for the annihilation operator on A site with spin up
  // kx4 = momentum along x for the annihilation operator on B site with spin down
  // ky4 = momentum along y for the annihilation operator on B site with spin down
  // kz4 = momentum along z for the annihilation operator on B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4);

  // compute the matrix element for the two body interaction between two sites A with opposite spins 
  //
  // kx1 = momentum along x for the creation operator on A site with spin up
  // ky1 = momentum along y for the creation operator on A site with spin up
  // kz1 = momentum along z for the creation operator on A site with spin up
  // kx2 = momentum along x for the creation operator on A site with spin down
  // ky2 = momentum along y for the creation operator on A site with spin down
  // kz2 = momentum along z for the creation operator on A site with spin down
  // kx3 = momentum along x for the annihilation operator on A site with spin up
  // ky3 = momentum along y for the annihilation operator on A site with spin up
  // kz3 = momentum along z for the annihilation operator on A site with spin up
  // kx4 = momentum along x for the annihilation operator on A site with spin down
  // ky4 = momentum along y for the annihilation operator on A site with spin down
  // kz4 = momentum along z for the annihilation operator on A site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4);

  // compute the matrix element for the two body interaction between two sites B with opposite spins 
  //
  // kx1 = momentum along x for the creation operator on B site with spin up
  // ky1 = momentum along y for the creation operator on B site with spin up
  // kz1 = momentum along z for the creation operator on B site with spin up
  // kx2 = momentum along x for the creation operator on B site with spin down
  // ky2 = momentum along y for the creation operator on B site with spin down
  // kz2 = momentum along z for the creation operator on B site with spin down
  // kx3 = momentum along x for the annihilation operator on B site with spin up
  // ky3 = momentum along y for the annihilation operator on B site with spin up
  // kz3 = momentum along z for the annihilation operator on B site with spin up
  // kx4 = momentum along x for the annihilation operator on B site with spin down
  // ky4 = momentum along y for the annihilation operator on B site with spin down
  // kz4 = momentum along z for the annihilation operator on B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4);
    

  // compute the one body hamiltonians related to the band stucture contribution
  //
  // oneBodyHamiltonians = array of one body hamiltonians
  virtual void ComputeOneBodyHamiltonian(HermitianMatrix* oneBodyHamiltonians);


};

#endif
