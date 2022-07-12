////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Yang-Le Wu                            //
//                                                                            //
//     class of 2-orbital square lattice model with interacting particles     //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 12/09/2011                      //
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


#ifndef PARTICLEONLATTICECHERN3TWOORBITALSINGLEBANDTHREEBODYHAMILTONIAN_H
#define PARTICLEONLATTICECHERN3TWOORBITALSINGLEBANDTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian
{

 protected:
  
  // imag part of the inter-orbital hopping amplitude between nearest neighbors along the x direction
  double NNHoppingInterX;
  // the inter-orbital hopping amplitude between nearest neighbors along the y direction
  double NNHoppingInterY;
  // the intra-orbital hopping amplitude between nearest neighbors
  double NNHoppingIntra;
  // four times the sublattice staggered chemical potential 
  double MuS;
  // nearest neighbor density-density potential strength
  double UPotential;
  // second nearest neighbor density-density potential strength
  double VPotential;
  // nearest neighbor density-density-density potential strength
  double WPotential;
  // next-to-nearest neighbor density-density-density potential strength
  double SPotential;
  // boundary condition twisting angle along x
  double GammaX;
  // boundary condition twisting angle along y
  double GammaY;

  // use flat band model
  bool FlatBand;
  
  // precalculation tables for cosine and sine factors
  Complex* XPhaseTable;
  Complex* YPhaseTable;
  Complex* XHalfPhaseTable;
  Complex* YHalfPhaseTable;
  int XPhaseTableShift;
  int YPhaseTableShift;

 public:

  // default constructor
  //
  ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body nearest neighbor interaction
  // vPotential = strength of the repulsive two body next nearest neighbor interaction
  // wPotential = strength of the repulsive three body neareast neighbor interaction
  // sPotential = strength of the repulsive three body next-to-nearest neighbor interaction
  // t1 = imag part of the inter-orbital hopping amplitude between nearest neighbors along the x direction
  // t2 = the inter-orbital hopping amplitude between nearest neighbors along the y direction
  // t3 = the intra-orbital hopping amplitude between nearest neighbors
  // mus = sublattice chemical potential on A sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
          double uPotential, Abstract2DTightBindingModel* tightBindingModel,bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeChern3TwoOrbitalTriangularLatticeSingleBandThreeBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();


};

#endif
