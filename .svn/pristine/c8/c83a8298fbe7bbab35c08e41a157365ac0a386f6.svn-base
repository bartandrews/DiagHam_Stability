///////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Yang-Le Wu                            //
//                                                                            //
//               class of Haldane model with interacting particles            //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 06/07/2011                      //
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


#ifndef PARTICLEONLATTICECHERN3TWOORBITALTRIANGULARLATTICESINGLEBANDHAMILTONIAN_H
#define PARTICLEONLATTICECHERN3TWOORBITALTRIANGULARLATTICESINGLEBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{

 protected:
  
  // nearest neighbor density-density potential strength
  double UPotential;
  // boundary condition twisting angle along x
  double GammaX;
  // boundary condition twisting angle along y
  double GammaY;

  // use flat band model
  bool FlatBand;
  
 public:

  // default constructor
  //
  ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // vPotential = strength of the repulsive two body second neareast neighbor interaction
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential,  
						     Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two sites A and B 
  //
  // kx1 = annihilation momentum along x for the B site
  // ky1 = annihilation momentum along y for the B site
  // kx2 = creation momentum along x for the B site
  // ky2 = creation momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the two body interaction between two A sites (or two B sites) 
  //
  // kx1 = annihilation momentum along x for the second site
  // ky1 = annihilation momentum along y for the second site
  // kx2 = creation momentum along x for the second site
  // ky2 = creation momentum along y for the second site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAA(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for on-site two body interaction involving A sites
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementOnSiteAA();

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
  Complex ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

};



#endif
