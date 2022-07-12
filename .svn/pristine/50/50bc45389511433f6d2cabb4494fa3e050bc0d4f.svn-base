////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 2d topological insulator based on the simple TI        //
//                       model and restricted to two bands                    //
//                                                                            //
//                        last modification : 27/09/2011                      //
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


#ifndef PARTICLEONSQUARELATTICETWOBANDSIMPLETICHECKERBOARDHAMILTONIAN_H
#define PARTICLEONSQUARELATTICETWOBANDSIMPLETICHECKERBOARDHAMILTONIAN_H


//added these two
//#ifndef PARTICLEONLATTICEDICELATTICETWOBANDHAMILTONIAN_H
//#define PARTICLEONLATTICEDICELATTICETWOBANDHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian.h" //this one is different from what we use. change?
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeDiceLatticeTwoBandHamiltonian : public ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian
{

 protected:
 
  // mass term of the simple TI model
  double Mass;

  // numerical factor for momentum along x
  double KxFactor;
  // numerical factor for momentum along y
  double KyFactor;
  
  // boundary condition twisting angle along x
  double GammaX;
  // boundary condition twisting angle along y
  double GammaY;
  // nearest neighbor density-density potential strength
  double WPotential;
  // strength of the repulsive two body on site interaction
  double VPotential;
  // strength of the repulsive two body different site opposite spin interaction
  //double WPotential;
  //Band Index
  int BandIndex;
  // use flat band model
  bool FlatBand;


  //Added these:

  double U3Potential;
  double U6Potential;
  int OneBodyInteractionFactors;
  double FactorU3;
  double FactorU6;

  

 public:

  // default constructor
  //
  ParticleOnLatticeDiceLatticeTwoBandHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = repulsive on-site potential strength between different orbitals
  // vPotential = repulsive on-site potential strength between opposite spins
  // mass = mass term of the simple TI model
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

  //Here I changed uPotential to u3Potential and vPotential to u6Potential
  ParticleOnLatticeDiceLatticeTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double u3Potential, double u6Potential, Abstract2DTightBindingModel* tightBindingModel,/*double mass,
						 double gammaX, double gammaY,*/ bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeDiceLatticeTwoBandHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two sites A with up spins
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

  // compute the matrix element for the two body interaction between two sites A with down spins
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

  // compute the matrix element for the two body interaction between two sites B with up spins
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

  // compute the matrix element for the two body interaction between two sites B with down spins
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

  // compute the matrix element for the two body interaction between two sites A and B  belonging to the same layer
  //
  // kx1 = momentum along x for the creation operator on A site with spin up
  // ky1 = momentum along y for the creation operator on A site with spin up
  // kx2 = momentum along x for the creation operator on B site with spin up
  // ky2 = momentum along y for the creation operator on B site with spin up
  // kx3 = momentum along x for the annihilation operator on A site with spin up
  // ky3 = momentum along y for the annihilation operator on A site with spin up
  // kx4 = momentum along x for the annihilation operator on B site with spin up
  // ky4 = momentum along y for the annihilation operator on B site with spin up
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two sites A and B with down spins
  //
  // kx1 = momentum along x for the creation operator on A site with spin down
  // ky1 = momentum along y for the creation operator on A site with spin down
  // kx2 = momentum along x for the creation operator on B site with spin down
  // ky2 = momentum along y for the creation operator on B site with spin down
  // kx3 = momentum along x for the annihilation operator on A site with spin down
  // ky3 = momentum along y for the annihilation operator on A site with spin down
  // kx4 = momentum along x for the annihilation operator on B site with spin down
  // ky4 = momentum along y for the annihilation operator on B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);
  
  // compute the matrix element for the two body interaction between two sites A and B with opposite spins
  //
  // kx1 = momentum along x for the creation operator on A site with spin down
  // ky1 = momentum along y for the creation operator on A site with spin down
  // kx2 = momentum along x for the creation operator on B site with spin up
  // ky2 = momentum along y for the creation operator on B site with spin up
  // kx3 = momentum along x for the annihilation operator on A site with spin down
  // ky3 = momentum along y for the annihilation operator on A site with spin down
  // kx4 = momentum along x for the annihilation operator on B site with spin up
  // ky4 = momentum along y for the annihilation operator on B site with spin up
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementADownBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two sites A and B with opposite spins
  //
  // kx1 = momentum along x for the creation operator on A site with spin up
  // ky1 = momentum along y for the creation operator on A site with spin up
  // kx2 = momentum along x for the creation operator on B site with spin down
  // ky2 = momentum along y for the creation operator on B site with spin down
  // kx3 = momentum along x for the annihilation operator on A site with spin up
  // ky3 = momentum along y for the annihilation operator on A site with spin up
  // kx4 = momentum along x for the annihilation operator on B site with spin down
  // ky4 = momentum along y for the annihilation operator on B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two sites A with opposite spins 
  //
  // kx1 = momentum along x for the creation operator on A site with spin up
  // ky1 = momentum along y for the creation operator on A site with spin up
  // kx2 = momentum along x for the creation operator on A site with spin down
  // ky2 = momentum along y for the creation operator on A site with spin down
  // kx3 = momentum along x for the annihilation operator on A site with spin up
  // ky3 = momentum along y for the annihilation operator on A site with spin up
  // kx4 = momentum along x for the annihilation operator on A site with spin down
  // ky4 = momentum along y for the annihilation operator on A site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two sites B with opposite spins 
  //
  // kx1 = momentum along x for the creation operator on B site with spin up
  // ky1 = momentum along y for the creation operator on B site with spin up
  // kx2 = momentum along x for the creation operator on B site with spin down
  // ky2 = momentum along y for the creation operator on B site with spin down
  // kx3 = momentum along x for the annihilation operator on B site with spin up
  // ky3 = momentum along y for the annihilation operator on B site with spin up
  // kx4 = momentum along x for the annihilation operator on B site with spin down
  // ky4 = momentum along y for the annihilation operator on B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the contributions to the interaction factor for on-site interactions
  // 
  // momentumIndex1 = compact momentum index of the first creation operator
  // momentumIndex2 = compact momentum index of the second creation operator
  // momentumIndex3 = compact momentum index of the first annihilation operator
  // momentumIndex4 = compact momentum index of the second annihiliation operator
  // energyIndex1 = energy index of the first creation operator
  // energyIndex2 = energy index of the second creation operator
  // energyIndex3 = energy index of the first annihilation operator
  // energyIndex4 = energy index of the second annihiliation operator
  Complex ComputeOnSiteContributions(ComplexMatrix* OneBodyBasis, double factorU3, double factorU6,
							int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
							int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4
							);
    

  // compute the transformation basis contribution to the interaction matrix element
  // 
  // oneBodyBasis = array of transformation basis matrices
  // momentumIndex1 = compact momentum index of the first creation operator
  // momentumIndex2 = compact momentum index of the second creation operator
  // momentumIndex3 = compact momentum index of the first annihilation operator
  // momentumIndex4 = compact momentum index of the second annihiliation operator
  // energyIndex1 = energy index of the first creation operator
  // energyIndex2 = energy index of the second creation operator
  // energyIndex3 = energy index of the first annihilation operator
  // energyIndex4 = energy index of the second annihiliation operator
  // siteIndex1 = site index of the first creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  // siteIndex2 = site index of the second creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  // siteIndex3 = site index of the first annihilation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  // siteIndex4 = site index of the second annihiliation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  virtual Complex ComputeTransfomationBasisContribution(ComplexMatrix* OneBodyBasis,
							int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
							int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
							int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4);

  // compute the one body transformation matrices and the optional one body band stucture contribution
  //
  // oneBodyBasis = array of one body transformation matrices
  virtual void ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis);

};

// compute the transformation basis contribution to the interaction matrix element
// 
// oneBodyBasis = array of transformation basis matrices
// momentumIndex1 = compact momentum index of the first creation operator
// momentumIndex2 = compact momentum index of the second creation operator
// momentumIndex3 = compact momentum index of the first annihilation operator
// momentumIndex4 = compact momentum index of the second annihiliation operator
// energyIndex1 = energy index of the first creation operator
// energyIndex2 = energy index of the second creation operator
// energyIndex3 = energy index of the first annihilation operator
// energyIndex4 = energy index of the second annihiliation operator
// siteIndex1 = site index of the first creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
// siteIndex2 = site index of the second creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
// siteIndex3 = site index of the first annihilation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
// siteIndex4 = site index of the second annihiliation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)

inline Complex ParticleOnLatticeDiceLatticeTwoBandHamiltonian::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
													int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
													int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
													int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1]) * Conj(oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3] * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4]);
}


#endif
