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


#ifndef PARTICLEONLATTICETWOBANDHOFSTADTERHAMILTONIAN_H
#define PARTICLEONLATTICETWOBANDHOFSTADTERHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeTwoBandHofstadterHamiltonian : public ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian
{

 protected:
 
  // nearest neighbor density-density potential strength
  double UPotential;
  // second nearest neighbor density-density potential strength
  double VPotential;

  // use flat band model
  bool FlatBand;

  // index of lower band
  int LowerBandIndex;

 public:

  // default constructor
  //
  ParticleOnLatticeTwoBandHofstadterHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrCellsX = number of sites in the x direction
  // nbrCellsY = number of sites in the y direction
  // bandIndex = index of lower band n to consider (as well as band n+1)
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // vPotential = strength of the repulsive two body second nearest neighbor interactio
  // tightBindingModel = pointer to the tight binding model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeTwoBandHofstadterHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrCellsX, int nbrCellsY, int bandIndex1, double uPotential, double vPotential,  Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeTwoBandHofstadterHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();
    

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
  virtual Complex ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
							int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
							int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
							int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4);

  // compute the matrix element for the two body interaction between and A sites and the adjacent B site 
  //
  // kx1 = annihilation momentum along x for the B site
  // ky1 = annihilation momentum along y for the B site
  // kx2 = creation momentum along x for the B site
  // ky2 = creation momentum along y for the B site
  // return value = corresponding matrix element
  
  Complex ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2);


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

inline Complex ParticleOnLatticeTwoBandHofstadterHamiltonian::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
													int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
													int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
													int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1]) * Conj(oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3] * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4]);
}


#endif
