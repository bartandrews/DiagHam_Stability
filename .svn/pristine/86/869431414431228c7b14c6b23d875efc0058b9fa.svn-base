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
//         in the single band approximation and four body interaction         // 
//                                                                            //
//                        last modification : 03/08/2011                      //
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


#ifndef PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICESINGLEBANDFOURBODYHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICESINGLEBANDFOURBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian : public ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian
{

 protected:
  
  
 public:

  // default constructor
  //
  ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive four body neareast neighbor interaction
  // vPotential = strength of the repulsive two body neareast neighbor interaction
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, 
								    Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the four body interaction between one site A and three sites B 
  //
  // kx1 = momentum along x of the creation operator of the A site
  // ky1 = momentum along y of the creation operator of the A site
  // kx2 = momentum along x of the creation operator of the first B site
  // ky2 = momentum along y of the creation operator of the first B site
  // kx3 = momentum along x of the creation operator of the second B site
  // ky3 = momentum along y of the creation operator of the second B site
  // kx4 = momentum along x of the creation operator of the third B site
  // ky4 = momentum along y of the creation operator of the third B site
  // kx5 = momentum along x of the annihilation operator of the A site
  // ky5 = momentum along y of the annihilation operator of the A site
  // kx6 = momentum along x of the annihilation operator of the first B site
  // ky6 = momentum along y of the annihilation operator of the first B site
  // kx7 = momentum along x of the annihilation operator of the second B site
  // ky7 = momentum along y of the annihilation operator of the second B site
  // kx8 = momentum along x of the annihilation operator of the third B site
  // ky8 = momentum along y of the annihilation operator of the third B site
  // return value = corresponding matrix element
  Complex ComputeFourBodyMatrixElementABBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8);

  // compute the matrix element for the four body interaction between three sites A and one site B 
  //
  // kx1 = momentum along x of the creation operator of the B site
  // ky1 = momentum along y of the creation operator of the B site
  // kx2 = momentum along x of the creation operator of the first A site
  // ky2 = momentum along y of the creation operator of the first A site
  // kx3 = momentum along x of the creation operator of the second A site
  // ky3 = momentum along y of the creation operator of the second A site
  // kx4 = momentum along x of the creation operator of the third A site
  // ky4 = momentum along y of the creation operator of the third A site
  // kx5 = momentum along x of the annihilation operator of the B site
  // ky5 = momentum along y of the annihilation operator of the B site
  // kx6 = momentum along x of the annihilation operator of the first A site
  // ky6 = momentum along y of the annihilation operator of the first A site
  // kx7 = momentum along x of the annihilation operator of the second A site
  // ky7 = momentum along y of the annihilation operator of the second A site
  // kx8 = momentum along x of the annihilation operator of the third A site
  // ky8 = momentum along y of the annihilation operator of the third A site
  // return value = corresponding matrix element
  Complex ComputeFourBodyMatrixElementBAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8);

  // compute the matrix element for the four body interaction between two sites A and two sites B 
  //
  // kx1 = momentum along x of the creation operator of the first A site
  // ky1 = momentum along y of the creation operator of the first A site
  // kx2 = momentum along x of the creation operator of the second A site
  // ky2 = momentum along y of the creation operator of the second A site
  // kx3 = momentum along x of the creation operator of the first B site
  // ky3 = momentum along y of the creation operator of the first B site
  // kx4 = momentum along x of the creation operator of the second B site
  // ky4 = momentum along y of the creation operator of the secnd B site
  // kx5 = momentum along x of the annihilation operator of the first A site
  // ky5 = momentum along y of the annihilation operator of the first A site
  // kx6 = momentum along x of the annihilation operator of the second B site
  // ky6 = momentum along y of the annihilation operator of the second B site
  // kx7 = momentum along x of the annihilation operator of the first A site
  // ky7 = momentum along y of the annihilation operator of the first A site
  // kx8 = momentum along x of the annihilation operator of the second A site
  // ky8 = momentum along y of the annihilation operator of the second A site
  // return value = corresponding matrix element
  Complex ComputeFourBodyMatrixElementAABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8);

  // compute the matrix element for the creation part of the four body interaction between two sites A and two sites B 
  //
  // kx1 = momentum along x of the creation operator of the first A site
  // ky1 = momentum along y of the creation operator of the first A site
  // kx2 = momentum along x of the creation operator of the second A site
  // ky2 = momentum along y of the creation operator of the second A site
  // kx3 = momentum along x of the creation operator of the first B site
  // ky3 = momentum along y of the creation operator of the first B site
  // kx4 = momentum along x of the creation operator of the second B site
  // ky4 = momentum along y of the creation operator of the secnd B site
  // return value = corresponding matrix element
  Complex ComputeFourBodyMatrixElementAABBIn1(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the creation part of the four body interaction between two sites A and two sites B 
  //
  // kx1 = momentum along x of the creation operator of the first A site
  // ky1 = momentum along y of the creation operator of the first A site
  // kx2 = momentum along x of the creation operator of the second A site
  // ky2 = momentum along y of the creation operator of the second A site
  // kx3 = momentum along x of the creation operator of the first B site
  // ky3 = momentum along y of the creation operator of the first B site
  // kx4 = momentum along x of the creation operator of the second B site
  // ky4 = momentum along y of the creation operator of the secnd B site
  // return value = corresponding matrix element  
  Complex ComputeFourBodyMatrixElementAABBIn2(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the annihilation part of the four body interaction between two sites A and two sites B 
  //
  // kx5 = momentum along x of the annihilation operator of the first A site
  // ky5 = momentum along y of the annihilation operator of the first A site
  // kx6 = momentum along x of the annihilation operator of the second B site
  // ky6 = momentum along y of the annihilation operator of the second B site
  // kx7 = momentum along x of the annihilation operator of the first A site
  // ky7 = momentum along y of the annihilation operator of the first A site
  // kx8 = momentum along x of the annihilation operator of the second A site
  // ky8 = momentum along y of the annihilation operator of the second A site
  // return value = corresponding matrix element
  Complex ComputeFourBodyMatrixElementAABBOut1(int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8);

  // compute the matrix element for the annihilation part of the four body interaction between two sites A and two sites B 
  //
  // kx5 = momentum along x of the annihilation operator of the first A site
  // ky5 = momentum along y of the annihilation operator of the first A site
  // kx6 = momentum along x of the annihilation operator of the second B site
  // ky6 = momentum along y of the annihilation operator of the second B site
  // kx7 = momentum along x of the annihilation operator of the first A site
  // ky7 = momentum along y of the annihilation operator of the first A site
  // kx8 = momentum along x of the annihilation operator of the second A site
  // ky8 = momentum along y of the annihilation operator of the second A site
  // return value = corresponding matrix element
  Complex ComputeFourBodyMatrixElementAABBOut2(int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8);

};

// compute the matrix element for the four body interaction between one site A and three sites B 
//
// kx1 = momentum along x of the creation operator of the A site
// ky1 = momentum along y of the creation operator of the A site
// kx2 = momentum along x of the creation operator of the first B site
// ky2 = momentum along y of the creation operator of the first B site
// kx3 = momentum along x of the creation operator of the second B site
// ky3 = momentum along y of the creation operator of the second B site
// kx4 = momentum along x of the creation operator of the third B site
// ky4 = momentum along y of the creation operator of the third B site
// kx5 = momentum along x of the annihilation operator of the A site
// ky5 = momentum along y of the annihilation operator of the A site
// kx6 = momentum along x of the annihilation operator of the first B site
// ky6 = momentum along y of the annihilation operator of the first B site
// kx7 = momentum along x of the annihilation operator of the second B site
// ky7 = momentum along y of the annihilation operator of the second B site
// kx8 = momentum along x of the annihilation operator of the third B site
// ky8 = momentum along y of the annihilation operator of the third B site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian::ComputeFourBodyMatrixElementABBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8)
{
  return 0.0;
  Complex Tmp = this->XPhaseTable[(kx8 - kx4) + this->XPhaseTableShift] * this->YPhaseTable[(ky6 - ky2) + this->YPhaseTableShift];
  Tmp += this->XPhaseTable[(kx7 + kx8 - kx3 - kx4) + this->XPhaseTableShift] * this->YPhaseTable[(ky8 - ky4) + this->YPhaseTableShift];
  Tmp += this->XPhaseTable[(kx6 + kx7 - kx2 - kx3) + this->XPhaseTableShift] * this->YPhaseTable[(ky7 + ky8 - ky3 - ky4) + this->YPhaseTableShift];
  Tmp += this->XPhaseTable[(kx6 - kx2) + this->XPhaseTableShift] * this->YPhaseTable[(ky6 + ky7 - ky2 - ky3) + this->YPhaseTableShift];
  Tmp *= this->XHalfPhaseTable[(kx2 + kx3 + kx4 - kx6 - kx7 - kx8) + this->XPhaseTableShift];
  Tmp *= this->YHalfPhaseTable[(ky2 + ky3 + ky4 - ky6 - ky7 - ky8) + this->YPhaseTableShift];
  return Tmp;
}

// compute the matrix element for the four body interaction between three sites A and one site B 
//
// kx1 = momentum along x of the creation operator of the B site
// ky1 = momentum along y of the creation operator of the B site
// kx2 = momentum along x of the creation operator of the first A site
// ky2 = momentum along y of the creation operator of the first A site
// kx3 = momentum along x of the creation operator of the second A site
// ky3 = momentum along y of the creation operator of the second A site
// kx4 = momentum along x of the creation operator of the third A site
// ky4 = momentum along y of the creation operator of the third A site
// kx5 = momentum along x of the annihilation operator of the B site
// ky5 = momentum along y of the annihilation operator of the B site
// kx6 = momentum along x of the annihilation operator of the first A site
// ky6 = momentum along y of the annihilation operator of the first A site
// kx7 = momentum along x of the annihilation operator of the second A site
// ky7 = momentum along y of the annihilation operator of the second A site
// kx8 = momentum along x of the annihilation operator of the third A site
// ky8 = momentum along y of the annihilation operator of the third A site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian::ComputeFourBodyMatrixElementBAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8)
{
  return 0.0;
  Complex Tmp = this->XPhaseTable[(kx3 + kx4 - kx7 - kx8) + this->XPhaseTableShift] * this->YPhaseTable[(ky4 - ky8) + this->YPhaseTableShift];
  Tmp += this->XPhaseTable[(kx2 + kx3 - kx6 - kx7) + this->XPhaseTableShift] * this->YPhaseTable[(ky3 + ky4 - ky7 -ky8) + this->YPhaseTableShift];
  Tmp += this->XPhaseTable[(kx2 - kx6) + this->XPhaseTableShift] * this->YPhaseTable[(ky2 + ky3 - ky6 - ky7) + this->YPhaseTableShift];
  Tmp += this->XPhaseTable[(kx4 - kx8) + this->XPhaseTableShift] * this->YPhaseTable[(ky2 - ky6) + this->YPhaseTableShift];
  Tmp *= this->XHalfPhaseTable[(kx1 - kx5) + this->XPhaseTableShift];
  Tmp *= this->YHalfPhaseTable[(ky1 - ky5) + this->YPhaseTableShift];
  return Tmp;
}

// compute the matrix element for the four body interaction between two sites A and two sites B 
//
// kx1 = momentum along x of the creation operator of the first A site
// ky1 = momentum along y of the creation operator of the first A site
// kx2 = momentum along x of the creation operator of the second A site
// ky2 = momentum along y of the creation operator of the second A site
// kx3 = momentum along x of the creation operator of the first B site
// ky3 = momentum along y of the creation operator of the first B site
// kx4 = momentum along x of the creation operator of the second B site
// ky4 = momentum along y of the creation operator of the secnd B site
// kx5 = momentum along x of the annihilation operator of the first A site
// ky5 = momentum along y of the annihilation operator of the first A site
// kx6 = momentum along x of the annihilation operator of the second B site
// ky6 = momentum along y of the annihilation operator of the second B site
// kx7 = momentum along x of the annihilation operator of the first A site
// ky7 = momentum along y of the annihilation operator of the first A site
// kx8 = momentum along x of the annihilation operator of the second A site
// ky8 = momentum along y of the annihilation operator of the second A site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian::ComputeFourBodyMatrixElementAABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8)
{
  Complex Tmp = this->XPhaseTable[(kx2 - kx6) + this->XPhaseTableShift] * this->YPhaseTable[(ky8 - ky4) + this->YPhaseTableShift];
  Tmp += this->XPhaseTable[(kx8 - kx4) + this->XPhaseTableShift] * this->YPhaseTable[(ky2 - ky6) + this->YPhaseTableShift];

  // interaction for ABAB losange configuration 
/*   Tmp += this->YPhaseTable[(ky2 + ky4 - ky6 - ky8) + this->YPhaseTableShift]; */
/*   Tmp += this->YPhaseTable[(ky2 + ky8 - ky6 - ky4) + this->YPhaseTableShift]; */
/*   Tmp += this->XPhaseTable[(kx7 + kx8 - kx3 - kx4) + this->XPhaseTableShift] * this->YPhaseTable[(ky2 + ky3 - ky6 - ky7) + this->YPhaseTableShift]; */
/*   Tmp += this->XPhaseTable[(kx7 + kx8 - kx3 - kx4) + this->XPhaseTableShift] * this->YPhaseTable[(ky2 + ky7 - ky6 - ky3) + this->YPhaseTableShift]; */

/*   Tmp += this->XPhaseTable[(kx2 + kx4 - kx6 - kx8) + this->XPhaseTableShift]; */
/*   Tmp += this->XPhaseTable[(kx2 + kx8 - kx6 - kx4) + this->XPhaseTableShift]; */
/*   Tmp += this->XPhaseTable[(kx2 + kx3 - kx6 - kx7) + this->XPhaseTableShift] * this->YPhaseTable[(ky7 + ky8 - ky3 - ky4) + this->YPhaseTableShift]; */
/*   Tmp += this->XPhaseTable[(kx2 + kx7 - kx6 - kx3) + this->XPhaseTableShift] * this->YPhaseTable[(ky7 + ky8 - ky3 - ky4) + this->YPhaseTableShift]; */

  Tmp *= this->XHalfPhaseTable[(kx3 + kx4 - kx7 - kx8) + this->XPhaseTableShift];
  Tmp *= this->YHalfPhaseTable[(ky3 + ky4 - ky7 - ky8) + this->YPhaseTableShift];
  return Tmp;
}

// compute the matrix element for the creation part of the four body interaction between two sites A and two sites B 
//
// kx1 = momentum along x of the creation operator of the first A site
// ky1 = momentum along y of the creation operator of the first A site
// kx2 = momentum along x of the creation operator of the second A site
// ky2 = momentum along y of the creation operator of the second A site
// kx3 = momentum along x of the creation operator of the first B site
// ky3 = momentum along y of the creation operator of the first B site
// kx4 = momentum along x of the creation operator of the second B site
// ky4 = momentum along y of the creation operator of the secnd B site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian::ComputeFourBodyMatrixElementAABBIn1(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = this->XPhaseTable[kx2 + this->XPhaseTableShift] * this->YPhaseTable[this->YPhaseTableShift - ky4];
  Tmp *= this->XHalfPhaseTable[(kx3 + kx4) + this->XPhaseTableShift];
  Tmp *= this->YHalfPhaseTable[(ky3 + ky4) + this->YPhaseTableShift];
  return Tmp;
}

// compute the matrix element for the creation part of the four body interaction between two sites A and two sites B 
//
// kx1 = momentum along x of the creation operator of the first A site
// ky1 = momentum along y of the creation operator of the first A site
// kx2 = momentum along x of the creation operator of the second A site
// ky2 = momentum along y of the creation operator of the second A site
// kx3 = momentum along x of the creation operator of the first B site
// ky3 = momentum along y of the creation operator of the first B site
// kx4 = momentum along x of the creation operator of the second B site
// ky4 = momentum along y of the creation operator of the secnd B site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian::ComputeFourBodyMatrixElementAABBIn2(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = this->XPhaseTable[this->XPhaseTableShift - kx4] * this->YPhaseTable[ky2 + this->YPhaseTableShift];
  Tmp *= this->XHalfPhaseTable[(kx3 + kx4) + this->XPhaseTableShift];
  Tmp *= this->YHalfPhaseTable[(ky3 + ky4) + this->YPhaseTableShift];
  return Tmp;
}

// compute the matrix element for the annihilation part of the four body interaction between two sites A and two sites B 
//
// kx5 = momentum along x of the annihilation operator of the first A site
// ky5 = momentum along y of the annihilation operator of the first A site
// kx6 = momentum along x of the annihilation operator of the second B site
// ky6 = momentum along y of the annihilation operator of the second B site
// kx7 = momentum along x of the annihilation operator of the first A site
// ky7 = momentum along y of the annihilation operator of the first A site
// kx8 = momentum along x of the annihilation operator of the second A site
// ky8 = momentum along y of the annihilation operator of the second A site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian::ComputeFourBodyMatrixElementAABBOut1(int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8)
{
  Complex Tmp = this->XPhaseTable[this->XPhaseTableShift - kx6] * this->YPhaseTable[ky8 + this->YPhaseTableShift];
  Tmp *= this->XHalfPhaseTable[this->XPhaseTableShift - kx7 - kx8];
  Tmp *= this->YHalfPhaseTable[this->YPhaseTableShift - ky7 - ky8];
  return Tmp;
}

// compute the matrix element for the annihilation part of the four body interaction between two sites A and two sites B 
//
// kx5 = momentum along x of the annihilation operator of the first A site
// ky5 = momentum along y of the annihilation operator of the first A site
// kx6 = momentum along x of the annihilation operator of the second B site
// ky6 = momentum along y of the annihilation operator of the second B site
// kx7 = momentum along x of the annihilation operator of the first A site
// ky7 = momentum along y of the annihilation operator of the first A site
// kx8 = momentum along x of the annihilation operator of the second A site
// ky8 = momentum along y of the annihilation operator of the second A site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian::ComputeFourBodyMatrixElementAABBOut2(int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8)
{
  Complex Tmp = this->XPhaseTable[kx8 + this->XPhaseTableShift] * this->YPhaseTable[this->YPhaseTableShift - ky6];
  Tmp *= this->XHalfPhaseTable[this->XPhaseTableShift - kx7 - kx8];
  Tmp *= this->YHalfPhaseTable[this->YPhaseTableShift - ky7 - ky8];
  return Tmp;
}


#endif
