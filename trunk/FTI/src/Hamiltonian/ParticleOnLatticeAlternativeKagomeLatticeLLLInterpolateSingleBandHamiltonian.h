////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Yang-Le Wu                            //
//                                                                            //
//            class of Kagome lattice model with interacting particles        //
//                       in the single band approximation                     //
//                           with interpolation to LLL                        //
//                                                                            //
//                        last modification : 14/10/2012                      //
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


#ifndef PARTICLEONLATTICEKAGOMELATTICELLLINTERPOLATESINGLEBANDHAMILTONIAN_H
#define PARTICLEONLATTICEKAGOMELATTICELLLINTERPOLATESINGLEBANDHAMILTONIAN_H

#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeAlternativeKagomeLatticeLLLInterpolateSingleBandHamiltonian : public ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian
{
public:

    // constructor
    //
    // particles = Hilbert space associated to the system
    // nbrParticles = number of particles
    // nbrSiteX = number of sites in the x direction
    // nbrSiteY = number of sites in the y direction
    // tightBindingModel = pointer to the tight binding model
    // uPotential = strength of the repulsive two body nearest neighbor interaction
    // vPotential = strength of the repulsive two body second nearest neighbor interaction
    // gammaX = boundary condition twisting angle along x
    // gammaY = boundary condition twisting angle along y
    // bandIndex = index of the band that has to be partially filled
    // flatBandFlag = use flat band model
    // interpolationToFQH = the interpolation parameter lambda between FQH and FCI, as in H = lambda * H_FQH + (1 - lambda) * H_FCI
    // twistAngle = angle between the two fundamental cycles of the torus in Radians
    // aspectRatio = aspect ratio of torus, Lx / Ly
    // nbrPseudopotentials = array of the number of pseudo-potentials
    // pseudoPotentials = array of the pseudo-potentials
    // architecture = architecture to use for precalculation
    // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
    ParticleOnLatticeAlternativeKagomeLatticeLLLInterpolateSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, Abstract2DTightBindingModel* tightBindingModel,
            double uPotential, double vPotential, int bandIndex, bool flatBandFlag,
            double interpolationToFQH, double twistAngle, double aspectRatio, int nbrPseudopotentials, double* pseudoPotentials,
            AbstractArchitecture* architecture, long memory = -1);

    // destructor
    //
    ~ParticleOnLatticeAlternativeKagomeLatticeLLLInterpolateSingleBandHamiltonian();
};

#endif

