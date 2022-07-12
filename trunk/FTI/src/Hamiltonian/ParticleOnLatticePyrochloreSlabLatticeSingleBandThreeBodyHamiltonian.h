////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//    class of pyrochlore slab lattice model with interacting particles       //
//      in the single band approximation and three body interaction           // 
//                                                                            //
//                        last modification : 26/05/2012                      //
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


#ifndef PARTICLEONLATTICEPYROCHLORESLABLATTICESINGLEBANDTHREEBODYHAMILTONIAN_H
#define PARTICLEONLATTICEPYROCHLORESLABLATTICESINGLEBANDTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticePyrochloreSlabLatticeSingleBandThreeBodyHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian
{

 protected:
  
  
  // strength of the repulsive two body onsite interaction
  double UPotential;
  // strength of the repulsive two body neareast neighbor interaction
  double VPotential;

  // index of the band to be filled
  int BandIndex;

  // use flat band model
  bool FlatBand;
  
 public:

  // constructor
  //
  ParticleOnLatticePyrochloreSlabLatticeSingleBandThreeBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body onsite interaction
  // vPotential = strength of the repulsive two body neareast neighbor interaction
  // tightBindingModel = pointer to the tight binding model
  // bandIndex = index of the band to consider
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticePyrochloreSlabLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, 
								       int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, 
								       Abstract2DTightBindingModel* tightBindingModel, int bandIndex,
								       bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticePyrochloreSlabLatticeSingleBandThreeBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

};



#endif
