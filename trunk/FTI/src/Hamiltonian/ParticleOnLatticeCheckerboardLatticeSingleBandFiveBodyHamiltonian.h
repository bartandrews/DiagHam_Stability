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
//         in the single band approximation and five body interaction         // 
//                                                                            //
//                        last modification : 07/08/2011                      //
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


#ifndef PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICESINGLEBANDFIVEBODYHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICESINGLEBANDFIVEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian : public ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian
{

 protected:
  
  
 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive five body neareast neighbor interaction
  // vPotential = strength of the repulsive two body neareast neighbor interaction
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, 
								    Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the five body interaction between one four A and three sites B 
  //
  // kx1 = momentum along x of the creation operator of the A site
  // ky1 = momentum along y of the creation operator of the A site
  // kx2 = momentum along x of the creation operator of the first B site
  // ky2 = momentum along y of the creation operator of the first B site
  // kx3 = momentum along x of the creation operator of the second B site
  // ky3 = momentum along y of the creation operator of the second B site
  // kx4 = momentum along x of the creation operator of the third B site
  // ky4 = momentum along y of the creation operator of the third B site
  // kx5 = momentum along x of the creation operator of the fourth B site
  // ky5 = momentum along y of the creation operator of the fourth B site
  // kx6 = momentum along x of the annihilation operator of the A site
  // ky6 = momentum along y of the annihilation operator of the A site
  // kx7 = momentum along x of the annihilation operator of the first B site
  // ky7 = momentum along y of the annihilation operator of the first B site
  // kx8 = momentum along x of the annihilation operator of the second B site
  // ky8 = momentum along y of the annihilation operator of the second B site
  // kx9 = momentum along x of the annihilation operator of the third B site
  // ky9 = momentum along y of the annihilation operator of the third B site
  // kx10 = momentum along x of the annihilation operator of the fourth B site
  // ky10 = momentum along y of the annihilation operator of the fourth B site
  // return value = corresponding matrix element
  Complex ComputeFiveBodyMatrixElementABBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8, int kx9, int ky9, int kx10, int ky10);

  // compute the matrix element for the creation part of the five body interaction between one site A and four sites B 
  //
  // kx1 = momentum along x of the creation operator of the A site
  // ky1 = momentum along y of the creation operator of the A site
  // kx2 = momentum along x of the creation operator of the first B site
  // ky2 = momentum along y of the creation operator of the first B site
  // kx3 = momentum along x of the creation operator of the second B site
  // ky3 = momentum along y of the creation operator of the second B site
  // kx4 = momentum along x of the creation operator of the third B site
  // ky4 = momentum along y of the creation operator of the third B site
  // kx5 = momentum along x of the creation operator of the fourth B site
  // ky5 = momentum along y of the creation operator of the fourth B site
  // kx6 = momentum along x of the creation operator of the A site
  // return value = corresponding matrix element
  inline Complex ComputeFiveBodyMatrixElementABBBBIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5);

  // compute the matrix element for the annihilation part of the five body interaction between one site A and four sites B 
  //
  // kx6 = momentum along x of the annihilation operator of the A site
  // ky6 = momentum along y of the annihilation operator of the A site
  // kx7 = momentum along x of the annihilation operator of the first B site
  // ky7 = momentum along y of the annihilation operator of the first B site
  // kx8 = momentum along x of the annihilation operator of the second B site
  // ky8 = momentum along y of the annihilation operator of the second B site
  // kx9 = momentum along x of the annihilation operator of the third B site
  // ky9 = momentum along y of the annihilation operator of the third B site
  // kx10 = momentum along x of the annihilation operator of the fourth B site
  // ky10 = momentum along y of the annihilation operator of the fourth B site
  // return value = corresponding matrix element
  Complex ComputeFiveBodyMatrixElementABBBBOut(int kx6, int ky6, int kx7, int ky7, int kx8, int ky8, int kx9, int ky9, int kx10, int ky10);

  // compute the matrix element for the five body interaction between four sites A and one site B 
  //
  // kx1 = momentum along x of the creation operator of the B site
  // ky1 = momentum along y of the creation operator of the B site
  // kx2 = momentum along x of the creation operator of the first A site
  // ky2 = momentum along y of the creation operator of the first A site
  // kx3 = momentum along x of the creation operator of the second A site
  // ky3 = momentum along y of the creation operator of the second A site
  // kx4 = momentum along x of the creation operator of the third A site
  // ky4 = momentum along y of the creation operator of the third A site
  // kx5 = momentum along x of the creation operator of the fourth A site
  // ky5 = momentum along y of the creation operator of the fourth A site
  // kx6 = momentum along x of the annihilation operator of the B site
  // ky6 = momentum along y of the annihilation operator of the B site
  // kx7 = momentum along x of the annihilation operator of the first A site
  // ky7 = momentum along y of the annihilation operator of the first A site
  // kx8 = momentum along x of the annihilation operator of the second A site
  // ky8 = momentum along y of the annihilation operator of the second A site
  // kx9 = momentum along x of the annihilation operator of the third A site
  // ky9 = momentum along y of the annihilation operator of the third A site
  // kx10 = momentum along x of the annihilation operator of the fourth A site
  // ky10 = momentum along y of the annihilation operator of the fourth A site
  // return value = corresponding matrix element
  Complex ComputeFiveBodyMatrixElementBAAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8, int kx9, int ky9, int kx10, int ky10);

  // compute the matrix element for creation part of the five body interaction between four sites A and one site B 
  //
  // kx1 = momentum along x of the creation operator of the B site
  // ky1 = momentum along y of the creation operator of the B site
  // kx2 = momentum along x of the creation operator of the first A site
  // ky2 = momentum along y of the creation operator of the first A site
  // kx3 = momentum along x of the creation operator of the second A site
  // ky3 = momentum along y of the creation operator of the second A site
  // kx4 = momentum along x of the creation operator of the third A site
  // ky4 = momentum along y of the creation operator of the third A site
  // kx5 = momentum along x of the creation operator of the fourth A site
  // ky5 = momentum along y of the creation operator of the fourth A site
  // return value = corresponding matrix element
  Complex ComputeFiveBodyMatrixElementBAAAAIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5);

  // compute the matrix element for the annihilation part of the five body interaction between four sites A and one site B 
  //
  // kx6 = momentum along x of the annihilation operator of the B site
  // ky6 = momentum along y of the annihilation operator of the B site
  // kx7 = momentum along x of the annihilation operator of the first A site
  // ky7 = momentum along y of the annihilation operator of the first A site
  // kx8 = momentum along x of the annihilation operator of the second A site
  // ky8 = momentum along y of the annihilation operator of the second A site
  // kx9 = momentum along x of the annihilation operator of the third A site
  // ky9 = momentum along y of the annihilation operator of the third A site
  // kx10 = momentum along x of the annihilation operator of the fourth A site
  // ky10 = momentum along y of the annihilation operator of the fourth A site
  // return value = corresponding matrix element
  Complex ComputeFiveBodyMatrixElementBAAAAOut(int kx6, int ky6, int kx7, int ky7, int kx8, int ky8, int kx9, int ky9, int kx10, int ky10);

};

// compute the matrix element for the five body interaction between one site A and four sites B 
//
// kx1 = momentum along x of the creation operator of the A site
// ky1 = momentum along y of the creation operator of the A site
// kx2 = momentum along x of the creation operator of the first B site
// ky2 = momentum along y of the creation operator of the first B site
// kx3 = momentum along x of the creation operator of the second B site
// ky3 = momentum along y of the creation operator of the second B site
// kx4 = momentum along x of the creation operator of the third B site
// ky4 = momentum along y of the creation operator of the third B site
// kx5 = momentum along x of the creation operator of the fourth B site
// ky5 = momentum along y of the creation operator of the fourth B site
// kx6 = momentum along x of the annihilation operator of the A site
// ky6 = momentum along y of the annihilation operator of the A site
// kx7 = momentum along x of the annihilation operator of the first B site
// ky7 = momentum along y of the annihilation operator of the first B site
// kx8 = momentum along x of the annihilation operator of the second B site
// ky8 = momentum along y of the annihilation operator of the second B site
// kx9 = momentum along x of the annihilation operator of the third B site
// ky9 = momentum along y of the annihilation operator of the third B site
// kx10 = momentum along x of the annihilation operator of the fourth B site
// ky10 = momentum along y of the annihilation operator of the fourth B site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian::ComputeFiveBodyMatrixElementABBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8, int kx9, int ky9, int kx10, int ky10)
{
  Complex Tmp = this->XPhaseTable[(kx8 + kx9 - kx3 - kx4) + this->XPhaseTableShift] * this->YPhaseTable[(ky9 + ky10 - ky4 - ky5) + this->YPhaseTableShift];;
  Tmp *= this->XHalfPhaseTable[(kx2 + kx3 + kx4 + kx5 - kx7 - kx8 - kx9 - kx10) + this->XPhaseTableShift];
  Tmp *= this->YHalfPhaseTable[(ky2 + ky3 + ky4 + ky5 - ky7 - ky8 - ky9 - ky10) + this->YPhaseTableShift];
  return Tmp;
}

// compute the matrix element for the creation part of the five body interaction between one site A and four sites B 
//
// kx1 = momentum along x of the creation operator of the A site
// ky1 = momentum along y of the creation operator of the A site
// kx2 = momentum along x of the creation operator of the first B site
// ky2 = momentum along y of the creation operator of the first B site
// kx3 = momentum along x of the creation operator of the second B site
// ky3 = momentum along y of the creation operator of the second B site
// kx4 = momentum along x of the creation operator of the third B site
// ky4 = momentum along y of the creation operator of the third B site
// kx5 = momentum along x of the creation operator of the fourth B site
// ky5 = momentum along y of the creation operator of the fourth B site
// kx6 = momentum along x of the creation operator of the A site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian::ComputeFiveBodyMatrixElementABBBBIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5)
{
  Complex Tmp = this->XPhaseTable[(- kx3 - kx4) + this->XPhaseTableShift] * this->YPhaseTable[(- ky4 - ky5) + this->YPhaseTableShift];;
  Tmp *= this->XHalfPhaseTable[(kx2 + kx3 + kx4 + kx5) + this->XPhaseTableShift];
  Tmp *= this->YHalfPhaseTable[(ky2 + ky3 + ky4 + ky5) + this->YPhaseTableShift];
  return Tmp;
}

// compute the matrix element for the annihilation part of the five body interaction between one site A and four sites B 
//
// kx6 = momentum along x of the annihilation operator of the A site
// ky6 = momentum along y of the annihilation operator of the A site
// kx7 = momentum along x of the annihilation operator of the first B site
// ky7 = momentum along y of the annihilation operator of the first B site
// kx8 = momentum along x of the annihilation operator of the second B site
// ky8 = momentum along y of the annihilation operator of the second B site
// kx9 = momentum along x of the annihilation operator of the third B site
// ky9 = momentum along y of the annihilation operator of the third B site
// kx10 = momentum along x of the annihilation operator of the fourth B site
// ky10 = momentum along y of the annihilation operator of the fourth B site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian::ComputeFiveBodyMatrixElementABBBBOut(int kx6, int ky6, int kx7, int ky7, int kx8, int ky8, int kx9, int ky9, int kx10, int ky10)
{
  Complex Tmp = this->XPhaseTable[(kx8 + kx9) + this->XPhaseTableShift] * this->YPhaseTable[(ky9 + ky10) + this->YPhaseTableShift];;
  Tmp *= this->XHalfPhaseTable[(- kx7 - kx8 - kx9 - kx10) + this->XPhaseTableShift];
  Tmp *= this->YHalfPhaseTable[(- ky7 - ky8 - ky9 - ky10) + this->YPhaseTableShift];
  return Tmp;
}

// compute the matrix element for the five body interaction between four sites A and one site B 
//
// kx1 = momentum along x of the creation operator of the B site
// ky1 = momentum along y of the creation operator of the B site
// kx2 = momentum along x of the creation operator of the first A site
// ky2 = momentum along y of the creation operator of the first A site
// kx3 = momentum along x of the creation operator of the second A site
// ky3 = momentum along y of the creation operator of the second A site
// kx4 = momentum along x of the creation operator of the third A site
// ky4 = momentum along y of the creation operator of the third A site
// kx5 = momentum along x of the creation operator of the fourth A site
// ky5 = momentum along y of the creation operator of the fourth A site
// kx6 = momentum along x of the annihilation operator of the B site
// ky6 = momentum along y of the annihilation operator of the B site
// kx7 = momentum along x of the annihilation operator of the first A site
// ky7 = momentum along y of the annihilation operator of the first A site
// kx8 = momentum along x of the annihilation operator of the second A site
// ky8 = momentum along y of the annihilation operator of the second A site
// kx9 = momentum along x of the annihilation operator of the third A site
// ky9 = momentum along y of the annihilation operator of the third A site
// kx10 = momentum along x of the annihilation operator of the fourth A site
// ky10 = momentum along y of the annihilation operator of the fourth A site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian::ComputeFiveBodyMatrixElementBAAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8, int kx9, int ky9, int kx10, int ky10)
{
  Complex Tmp = this->XPhaseTable[(kx3 + kx4 - kx8 - kx9) + this->XPhaseTableShift] * this->YPhaseTable[(ky4 + ky5 - ky9 -ky10) + this->YPhaseTableShift];;
  Tmp *= this->XHalfPhaseTable[(kx1 - kx6) + this->XPhaseTableShift];
  Tmp *= this->YHalfPhaseTable[(ky1 - ky6) + this->YPhaseTableShift];
  return Tmp;
}

// compute the matrix element for creation part of the five body interaction between four sites A and one site B 
//
// kx1 = momentum along x of the creation operator of the B site
// ky1 = momentum along y of the creation operator of the B site
// kx2 = momentum along x of the creation operator of the first A site
// ky2 = momentum along y of the creation operator of the first A site
// kx3 = momentum along x of the creation operator of the second A site
// ky3 = momentum along y of the creation operator of the second A site
// kx4 = momentum along x of the creation operator of the third A site
// ky4 = momentum along y of the creation operator of the third A site
// kx5 = momentum along x of the creation operator of the fourth A site
// ky5 = momentum along y of the creation operator of the fourth A site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian::ComputeFiveBodyMatrixElementBAAAAIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5)
{
  Complex Tmp = this->XPhaseTable[kx3 + kx4 + this->XPhaseTableShift] * this->YPhaseTable[ky4 + ky5 + this->YPhaseTableShift];;
  Tmp *= this->XHalfPhaseTable[kx1 + this->XPhaseTableShift];
  Tmp *= this->YHalfPhaseTable[ky1 + this->YPhaseTableShift];
  return Tmp;
}

// compute the matrix element for the annihilation part of the five body interaction between four sites A and one site B 
//
// kx6 = momentum along x of the annihilation operator of the B site
// ky6 = momentum along y of the annihilation operator of the B site
// kx7 = momentum along x of the annihilation operator of the first A site
// ky7 = momentum along y of the annihilation operator of the first A site
// kx8 = momentum along x of the annihilation operator of the second A site
// ky8 = momentum along y of the annihilation operator of the second A site
// kx9 = momentum along x of the annihilation operator of the third A site
// ky9 = momentum along y of the annihilation operator of the third A site
// kx10 = momentum along x of the annihilation operator of the fourth A site
// ky10 = momentum along y of the annihilation operator of the fourth A site
// return value = corresponding matrix element

inline Complex ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian::ComputeFiveBodyMatrixElementBAAAAOut(int kx6, int ky6, int kx7, int ky7, int kx8, int ky8, int kx9, int ky9, int kx10, int ky10)
{
  Complex Tmp = this->XPhaseTable[(- kx8 - kx9) + this->XPhaseTableShift] * this->YPhaseTable[(- ky9 -ky10) + this->YPhaseTableShift];;
  Tmp *= this->XHalfPhaseTable[(- kx6) + this->XPhaseTableShift];
  Tmp *= this->YHalfPhaseTable[(- ky6) + this->YPhaseTableShift];
  return Tmp;
}


#endif
