////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//       class of a two body interaction projected onto two bands with        //
//      spin-like degree of freedom (requiring only U(1) conservation)        //
//      from an ASCII file providing the two body matrix real elements        //
//                                                                            //
//                        last modification : 07/05/2020                      //
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


#ifndef PARTICLEONLATTICEFROMFILEINTERACTIONTWOBANDWITHSPINREALHAMILTONIAN_H
#define PARTICLEONLATTICEFROMFILEINTERACTIONTWOBANDWITHSPINREALHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeFromFileInteractionTwoBandWithSpinRealHamiltonian : public ParticleOnLatticeFromFileInteractionTwoBandRealHamiltonian
{

 protected:
 
 public:

  // default constructor
  //
  ParticleOnLatticeFromFileInteractionTwoBandWithSpinRealHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // matrixElementsInteractionFile = name of the ASCII file containing the matrix element for the generic two body interaction term
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // interactionRescalingFactor = global rescaling factor for the two-body interaction term
  // additionalSpinFlag = include an additional spin 1/2 degree of freedom, building an SU(2) invariant interaction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeFromFileInteractionTwoBandWithSpinRealHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY,
								     char* matrixElementsInteractionFile,
								     Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag,
								     double interactionRescalingFactor, 
								     bool additionalSpinFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeFromFileInteractionTwoBandWithSpinRealHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();
  

};


#endif
