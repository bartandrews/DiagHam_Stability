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


#ifndef PARTICLEONLATTICESQUARELATTICETWOORBITALSINGLEBANDTHREEBODYHAMILTONIAN_H
#define PARTICLEONLATTICESQUARELATTICETWOORBITALSINGLEBANDTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandThreeBodyHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian
{

 protected:

  // nearest neighbor density-density potential strength
  double UPotential;
  // second nearest neighbor density-density potential strength
  double VPotential;
  // nearest neighbor density-density-density potential strength
  double WPotential;
  // next-to-nearest neighbor density-density-density potential strength
  double SPotential;

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
  ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandThreeBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // tightBindingModel = pointer to the tight binding model
  // uPotential = strength of the repulsive two body nearest neighbor interaction
  // vPotential = strength of the repulsive two body next nearest neighbor interaction
  // wPotential = strength of the repulsive three body neareast neighbor interaction
  // sPotential = strength of the repulsive three body next-to-nearest neighbor interaction
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, Abstract2DTightBindingModel* tightBindingModel, 
          double uPotential, double vPotential, double wPotential, double sPotential, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeSquareLatticeTwoOrbitalSingleBandThreeBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute all the phase precalculation arrays 
  //
  virtual void ComputePhaseArray();

  // compute the matrix element for the two body interaction between two NN sites
  //
  // kx2 = annihilation momentum along x for the second site
  // ky2 = annihilation momentum along y for the second site
  // kx4 = creation momentum along x for the second site
  // ky4 = creation momentum along y for the second site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementNN(int kx2, int ky2, int kx4, int ky4);

  // compute the matrix element for the three body interaction between three orbitals on two nearest neighbor sites
  //
  // kx2 = annihilation momentum along x for the A orbital on the second site
  // ky2 = annihilation momentum along y for the A orbital on the second site
  // kx3 = annihilation momentum along x for the B orbital on the second site
  // ky3 = annihilation momentum along y for the B orbital on the second site
  // kx5 = creation momentum along x for the A orbital on the second site
  // ky5 = creation momentum along y for the A orbital on the second site
  // kx6 = creation momentum along x for the B orbital on the second site
  // ky6 = creation momentum along y for the B orbital on the second site
  // return value = corresponding matrix element
  Complex ComputeThreeBodyMatrixElementNN(int kx2, int ky2, int kx3, int ky3, int kx5, int ky5, int kx6, int ky6);

};

#endif
