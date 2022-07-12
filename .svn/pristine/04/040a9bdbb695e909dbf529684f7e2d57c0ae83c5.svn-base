////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//              class of ruby lattice model with interacting particles        //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 17/10/2011                      //
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


#ifndef PARTICLEONLATTICEWITHSPINRUBYLATTICESINGLEBANDHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINRUBYLATTICESINGLEBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeRubyLatticeSingleBandHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{

 protected:
  
  // nearest neighbor density-density potential strength
  double UPotential;

  // next nearest neighbor density-density potential strength
  double VPotential;

  // use flat band model
  bool FlatBand;
  
 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // vPotential = strength of the repulsive two body next neareast neighbor interaction
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeRubyLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, 
						    Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeRubyLatticeSingleBandHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two sites A1 and A2 
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementA1A2();

  // compute the matrix element for the two body interaction between two sites A1 and A3 
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementA1A3();

  // compute the matrix element for the two body interaction between two sites A1 and A5
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementA1A5();

  // compute the matrix element for the two body interaction between two sites A1 and A6 
  //
  // k1x = creation momentum along x for the A6 site
  // k1y = creation momentum along y for the A6 site
  // k2x = annihilation momentum along x for the A6 site
  // k2y = annihilation momentum along y for the A6 site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementA1A6(int k1x, int k1y, int k2x, int k2y);

  // compute the matrix element for the two body interaction between two sites A2 and A3 
  //
  // k1x = creation momentum along x for the A3 site
  // k1y = creation momentum along y for the A3 site
  // k2x = annihilation momentum along x for the A3 site
  // k2y = annihilation momentum along y for the A3 site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementA2A3(int k1x, int k1y, int k2x, int k2y);

  // compute the matrix element for the two body interaction between two sites A2 and A4
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementA2A4();

  // compute the matrix element for the two body interaction between two sites A2 and A6
  //
  // return value = corresponding matrix element  
  Complex ComputeTwoBodyMatrixElementA2A6();

  // compute the matrix element for the two body interaction between two sites A3 and A4 
  //
  // k1x = creation momentum along x for the A4 site
  // k1y = creation momentum along y for the A4 site
  // k2x = annihilation momentum along x for the A4 site
  // k2y = annihilation momentum along y for the A4 site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementA3A4(int k1x, int k1y, int k2x, int k2y);

  // compute the matrix element for the two body interaction between two sites A3 and A5
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementA3A5();

  // compute the matrix element for the two body interaction between two sites A4 and A5
  //
  // return value = corresponding matrix element  
  Complex ComputeTwoBodyMatrixElementA4A5();

  // compute the matrix element for the two body interaction between two sites A4 and A6
  //
  // return value = corresponding matrix element  
  Complex ComputeTwoBodyMatrixElementA4A6();

  // compute the matrix element for the two body interaction between two sites A5 and A6 
  //
  // k1x = creation momentum along x for the A6 site
  // k1y = creation momentum along y for the A6 site
  // k2x = annihilation momentum along x for the A6 site
  // k2y = annihilation momentum along y for the A6 site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementA5A6(int k1x, int k1y, int k2x, int k2y);

  // compute the matrix element for on-site two body interaction involving A1 sites
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementOnSiteA1A1();

  // compute the matrix element for on-site two body interaction involving A2 sites
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementOnSiteA2A2();

  // compute the matrix element for on-site two body interaction involving A3 sites
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementOnSiteA3A3();

  // compute the matrix element for on-site two body interaction involving A4 sites
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementOnSiteA4A4();

  // compute the matrix element for on-site two body interaction involving A5 sites
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementOnSiteA5A5();

  // compute the matrix element for on-site two body interaction involving A6 sites
  //
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementOnSiteA6A6();


};



#endif
