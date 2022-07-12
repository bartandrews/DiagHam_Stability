////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                           class author: Yang-Le Wu                         //
//                                                                            //
//               class of Haldane model with interacting particles            //
//         in the single band approximation and four body interaction         // 
//                                                                            //
//                        last modification : 16/08/2011                      //
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


#ifndef PARTICLEONLATTICEHALDANEMODELSINGLEBANDFOURBODYHAMILTONIAN_H
#define PARTICLEONLATTICEHALDANEMODELSINGLEBANDFOURBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeHaldaneModelSingleBandFourBodyHamiltonian : public ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian
{

 protected:
  
  
 public:

  // default constructor
  //
  ParticleOnLatticeHaldaneModelSingleBandFourBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // vPotential = strength of the repulsive two body second neareast neighbor interaction
  // wPotential = strength of the repulsive three body neareast neighbor interaction
  // sPotential = strength of the repulsive three body next-to-nearest neighbor interaction
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeHaldaneModelSingleBandFourBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
							     double uPotential, double vPotential, double wPotential, double sPotential,
							     Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeHaldaneModelSingleBandFourBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the annihilation part of the four body interaction between one site A and three sites B 
  //
  // kx1 = momentum along x of the annihilation operator of the A1 site
  // ky1 = momentum along y of the annihilation operator of the A1 site
  // kx2 = momentum along x of the annihilation operator of the B1 site
  // ky2 = momentum along y of the annihilation operator of the B1 site
  // kx3 = momentum along x of the annihilation operator of the B2 site
  // ky3 = momentum along y of the annihilation operator of the B2 site
  // kx4 = momentum along x of the annihilation operator of the B3 site
  // ky4 = momentum along y of the annihilation operator of the B3 site
  // return value = corresponding matrix element
  Complex ComputeFourBodyMatrixElementABBBIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the annihilation part of the four body interaction between one site A and three sites B 
  //
  // kx5 = momentum along x of the creation operator of the A1 site
  // ky5 = momentum along y of the creation operator of the A1 site
  // kx6 = momentum along x of the creation operator of the B1 site
  // ky6 = momentum along y of the creation operator of the B1 site
  // kx7 = momentum along x of the creation operator of the B2 site
  // ky7 = momentum along y of the creation operator of the B2 site
  // kx8 = momentum along x of the creation operator of the B3 site
  // ky8 = momentum along y of the creation operator of the B3 site
  // return value = corresponding matrix element
  Complex ComputeFourBodyMatrixElementABBBOut(int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8);

};

// compute the matrix element for the annihilation part of the four body interaction between one site A and three sites B 
//
// kx1 = momentum along x of the annihilation operator of the A1 site
// ky1 = momentum along y of the annihilation operator of the A1 site
// kx2 = momentum along x of the annihilation operator of the B1 site
// ky2 = momentum along y of the annihilation operator of the B1 site
// kx3 = momentum along x of the annihilation operator of the B2 site
// ky3 = momentum along y of the annihilation operator of the B2 site
// kx4 = momentum along x of the annihilation operator of the B3 site
// ky4 = momentum along y of the annihilation operator of the B3 site
// return value = corresponding matrix element
inline Complex ParticleOnLatticeHaldaneModelSingleBandFourBodyHamiltonian::ComputeFourBodyMatrixElementABBBIn(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
    Complex Tmp = this->XPhaseTable[this->XPhaseTableShift + kx2] * this->YPhaseTable[this->YPhaseTableShift + (ky2 + ky3)];
    return Tmp;
}

// compute the matrix element for the annihilation part of the four body interaction between one site A and three sites B 
//
// kx5 = momentum along x of the creation operator of the A1 site
// ky5 = momentum along y of the creation operator of the A1 site
// kx6 = momentum along x of the creation operator of the B1 site
// ky6 = momentum along y of the creation operator of the B1 site
// kx7 = momentum along x of the creation operator of the B2 site
// ky7 = momentum along y of the creation operator of the B2 site
// kx8 = momentum along x of the creation operator of the B3 site
// ky8 = momentum along y of the creation operator of the B3 site
// return value = corresponding matrix element
inline Complex ParticleOnLatticeHaldaneModelSingleBandFourBodyHamiltonian::ComputeFourBodyMatrixElementABBBOut(int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8)
{
    Complex Tmp = this->XPhaseTable[this->XPhaseTableShift - kx6] * this->YPhaseTable[this->YPhaseTableShift - (ky6 + ky7)];
    return Tmp;
}

#endif
