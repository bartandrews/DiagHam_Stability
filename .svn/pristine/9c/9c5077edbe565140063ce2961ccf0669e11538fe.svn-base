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
//                       using two decoupled kagome models                    //
//                                                                            //
//                        last modification : 04/06/2013                      //
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


#ifndef PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDDECOUPLEDKAGOMEHAMILTONIANTILTED_H
#define PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDDECOUPLEDKAGOMEHAMILTONIANTILTED_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonianTilted : public ParticleOnLatticeWithSpinChernInsulatorHamiltonian
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
  // use two copies of Kagome with time reversal symmetry
  bool TimeReversal;
  // use a linear combination of the tight binding model and its time reversal conjugate for the down copy
  bool LinearCombinationFlag;
  // coefficietn of the linear combination of the tight binding model and its time reversal conjugate for the down copy
  double LinearCombinationCoefficient;  
  // pointer to the tight binding model of the first copy (must be in Bloch form)
  Abstract2DTightBindingModel* TightBindingModelUp;
  // pointer to the tight binding model of the second copy (must be in Bloch form)
  Abstract2DTightBindingModel* TightBindingModelDown;

 public:
   
  // default constructor
  //
  ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonianTilted();

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
  ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonianTilted(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
									  double uPotential, double vPotential, double wPotential,
									  Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, bool timeReversalFlag, 
									  AbstractArchitecture* architecture, long memory = -1);

  // constructor with the one-body Hamiltonian for the down copy is explicitly provided
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // vPotential = strength of the repulsive on site two body interaction between opposite spins
  // wPotential = strength of the repulsive two body neareast neighbor interaction between opposite spins
  // tightBindingModelUp = pointer to the tight binding model of the first copy (must be in Bloch form)
  // tightBindingModelDown = pointer to the tight binding model of the second copy (must be in Bloch form)
  // flatBandFlag = use flat band model
  // timeReversalFlag = apply thge time reversal symmetry on the second copy of the tight binding model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)  
  ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonianTilted(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
									  double uPotential, double vPotential, double wPotential,
									  Abstract2DTightBindingModel* tightBindingModelUp, Abstract2DTightBindingModel* tightBindingModelDown, bool flatBandFlag, bool timeReversalFlag, 
									  AbstractArchitecture* architecture, long memory = -1);
  
  // destructor
  //
  ~ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonianTilted();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two sites A and B belonging to the same up layer
  //
  // kx1 = creation momentum along x for the B site
  // ky1 = creation momentum along y for the B site
  // kx2 = annihilation momentum along x for the B site
  // ky2 = annihilation momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kx2, int ky2);


  // compute the matrix element for the two body interaction between two sites A and C belonging to the same up layer
  //
  // kx1 = creation momentum along x for the C site
  // ky1 = creation momentum along y for the C site
  // kx2 = annihilation momentum along x for the C site
  // ky2 = annihilation momentum along y for the C site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAUpCUp(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the two body interaction between two sites B and C belonging to the same up layer
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
  Complex ComputeTwoBodyMatrixElementBUpCUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two sites A and B belonging to the same down layer
  //
  // kx1 = creation momentum along x for the B site
  // ky1 = creation momentum along y for the B site
  // kx2 = annihilation momentum along x for the B site
  // ky2 = annihilation momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the two body interaction between two sites A and C belonging to the same down layer
  //
  // kx1 = creation momentum along x for the C site
  // ky1 = creation momentum along y for the C site
  // kx2 = annihilation momentum along x for the C site
  // ky2 = annihilation momentum along y for the C site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementADownCDown(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the two body interaction between two sites B and C belonging to the same down layer
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
  Complex ComputeTwoBodyMatrixElementBDownCDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two sites A and B belonging to different layers
  //
  // kx1 = creation momentum along x for the B site
  // ky1 = creation momentum along y for the B site
  // kx2 = annihilation momentum along x for the B site
  // ky2 = annihilation momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the two body interaction between two sites A and C belonging to different layers
  //
  // kx1 = creation momentum along x for the C site
  // ky1 = creation momentum along y for the C site
  // kx2 = annihilation momentum along x for the C site
  // ky2 = annihilation momentum along y for the C site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAUpCDown(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the two body interaction between two sites B and C belonging to different layers
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
  Complex ComputeTwoBodyMatrixElementBUpCDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two sites A with different layer indices 
  //
  // kx1 = first creation momentum along x for the A site
  // ky1 = first creation momentum along y for the A site
  // kx2 = second creation momentum along x for the A site
  // ky2 = second creation momentum along y for the A site
  // kx3 = first annihilation momentum along x for the A site
  // ky3 = first annihilation momentum along y for the A site
  // kx4 = second annihilation momentum along x for the A site
  // ky4 = second annihilation momentum along y for the A site
  Complex ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);
    
  // compute the matrix element for the two body interaction between two sites B with different layer indices 
  //
  // kx1 = first creation momentum along x for the B site
  // ky1 = first creation momentum along y for the B site
  // kx2 = second creation momentum along x for the B site
  // ky2 = second creation momentum along y for the B site
  // kx3 = first annihilation momentum along x for the B site
  // ky3 = first annihilation momentum along y for the B site
  // kx4 = second annihilation momentum along x for the B site
  // ky4 = second annihilation momentum along y for the B site
  Complex ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);
    
  // compute the matrix element for the two body interaction between two sites C with different layer indices 
  //
  // kx1 = first creation momentum along x for the C site
  // ky1 = first creation momentum along y for the C site
  // kx2 = second creation momentum along x for the C site
  // ky2 = second creation momentum along y for the C site
  // kx3 = first annihilation momentum along x for the C site
  // ky3 = first annihilation momentum along y for the C site
  // kx4 = second annihilation momentum along x for the C site
  // ky4 = second annihilation momentum along y for the C site
  Complex ComputeTwoBodyMatrixElementCUpCDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);
    
  // compute the matrix element for on-site two body interaction involving A sites and spin up-up
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementOnSiteAUpAUp();

  // compute the matrix element for on-site two body interaction involving A sites and spin up-down
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementOnSiteAUpADown();
  
  // compute the matrix element for on-site two body interaction involving A sites and spin down-down
  //
  // return value = corresponding matrix element  
  Complex ComputeTwoBodyMatrixElementOnSiteADownADown();

  // compute the matrix element for on-site two body interaction involving B sites and spin up-up
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
  Complex ComputeTwoBodyMatrixElementOnSiteBUpBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for on-site two body interaction involving B sites and spin up-down
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
  Complex ComputeTwoBodyMatrixElementOnSiteBUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for on-site two body interaction involving B sites and spin down-down
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
  Complex ComputeTwoBodyMatrixElementOnSiteBDownBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for on-site two body interaction involving C sites and spin up-up
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
  Complex ComputeTwoBodyMatrixElementOnSiteCUpCUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for on-site two body interaction involving C sites and spin up-down
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
  Complex ComputeTwoBodyMatrixElementOnSiteCUpCDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for on-site two body interaction involving C sites and spin down-down
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
  Complex ComputeTwoBodyMatrixElementOnSiteCDownCDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the one body transformation matrices and the optional one body band stucture contribution
  //
  // oneBodyBasis = array of one body transformation matrices (the leftmost upper block for the spin up, the rightmsdt lower block for the spin down)
  void ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis);

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
  Complex ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
						int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
						int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
						int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4);

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

inline Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonianTilted::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
														     int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
														     int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
														     int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1]) * oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3] * Conj(oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4]);
}


#endif
