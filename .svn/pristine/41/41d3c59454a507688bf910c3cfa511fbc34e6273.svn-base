////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                          class author: Yang-Le Wu                          //
//                                                                            //
//     class of 2-orbital square lattice model with interacting particles     //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 11/09/2011                      //
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


#ifndef PARTICLEONLATTICESQUARELATTICETWOORBITALSINGLEBANDHAMILTONIAN_H
#define PARTICLEONLATTICESQUARELATTICETWOORBITALSINGLEBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{

 protected:
  
  // nearest neighbor density-density potential strength
  double UPotential;
  // nearest neighbor density-density potential strength AB for bosons
  double UABPotential;
  // second nearest neighbor density-density potential strength
  double VPotential;

  // use flat band model
  bool FlatBand;
  
 public:

  // default constructor
  //
  ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // tightBindingModel = pointer to the tight binding model
  // uPotential = strength of the repulsive two body nearest neighbor interaction
  // vPotential = strength of the repulsive two body next nearest neighbor interaction
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, Abstract2DTightBindingModel* tightBindingModel,
          double uPotential, double uabPotential, double vPotential, bool flatBandFlag, AbstractArchitecture* architecture, long memory=-1);

  // destructor
  //
  ~ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two NN sites
  //
  // kx2 = annihilation momentum along x for the second site
  // ky2 = annihilation momentum along y for the second site
  // kx4 = creation momentum along x for the second site
  // ky4 = creation momentum along y for the second site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementNN(int kx2, int ky2, int kx4, int ky4);

};



#endif
