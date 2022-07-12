////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//      class of checkerboard lattice model with interacting particles        //
//                                                                            //
//                        last modification : 03/04/2011                      //
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


#ifndef PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICEHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICEHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class AbstractArchitecture;


class ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian : public ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian
{

 protected:
  
  //strength of the NN interaction
  double VPotential;
  
  //  gap between the first band and the second band when using the flat band model   
  double FlatBandOneBodyGap;

 public:

  // default constructor
  //
  ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // t2p = hoping amplitude between second next neareast neighbor sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // flatBandFlag = use flat band model
  // flatBandOneBodyGap = set the gap between the first band and the second band when using the flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
							  Abstract2DTightBindingModel* tightBindingModel, double uPotential, double vPotential, 
							  bool flatBandFlag, double flatBandOneBodyGap, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian();
  
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
  Complex ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2);
  
  
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
// siteIndex1 = site index of the first creation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex2 = site index of the second creation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex3 = site index of the first annihilation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex4 = site index of the second annihiliation operator (0,...,3 = up, 4,...,7 = down)
  Complex ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis, int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,  int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4);

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
// siteIndex1 = site index of the first creation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex2 = site index of the second creation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex3 = site index of the first annihilation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex4 = site index of the second annihiliation operator (0,...,3 = up, 4,...,7 = down)

inline Complex ParticleOnLatticeWithSpinCheckerboardLatticeHamiltonian::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
													    int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
													    int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
													    int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1]) * Conj(oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3] * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4]);
}

#endif
