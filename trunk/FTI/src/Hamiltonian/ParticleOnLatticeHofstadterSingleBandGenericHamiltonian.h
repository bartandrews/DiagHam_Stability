////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                  class author: Gunnar Möller & Ben Huddart                 //
//                                                                            //
//      class of Hofstadter model with particles interacting with a generic   //
//           two-body interaction in the single band approximation            // 
//                                                                            //
//                        last modification : 16/01/2017                      //
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


#ifndef PARTICLEONLATTICEHOFSTADTERSINGLEBANDGENERICHAMILTONIAN_H
#define PARTICLEONLATTICEHOFSTADTERSINGLEBANDGENERICHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"

#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"


#include <iostream>


using std::ostream;
using std::cout;
using std::endl;

class AbstractTightBindingInteraction;

class ParticleOnLatticeHofstadterSingleBandGenericHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{

 protected:

  // Object encoding shape of interactions
  AbstractTightBindingInteraction* Interaction;

  // index of the band to be filled
  int BandIndex;

  // use flat band model
  bool FlatBand;
  
 public:

  // constructor
  //
  ParticleOnLatticeHofstadterSingleBandGenericHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrCellsX = number of sites in the x direction
  // nbrCellsY = number of sites in the y direction
  // bandIndex = index of band to consider
  // genericInteraction = pointer to object encoding interactions
  // tightBindingModel = pointer to the tight binding model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeHofstadterSingleBandGenericHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellsX, int nbrCellsY, int bandIndex, AbstractTightBindingInteraction* interaction,  Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeHofstadterSingleBandGenericHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for on-site two body interaction involving sites on generic sublattic 
  //
  // subl = sublattice index
  // kx1 = first creation momentum along x for the B site
  // ky1 = first creation momentum along y for the B site
  // kx2 = second creation momentum along x for the B site
  // ky2 = second creation momentum along y for the B site
  // kx3 = first annihilation momentum along x for the B site
  // ky3 = first annihilation momentum along y for the B site
  // kx4 = second annihilation momentum along x for the B site
  // ky4 = second annihilation momentum along y for the B site
  //
  // return value = corresponding matrix element
  Complex ComputeEmbeddingOnSite(int subl, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for on-site two body interaction involving sites on generic sublattic 
  //
  // kx1 = first creation momentum along x for the B site
  // ky1 = first creation momentum along y for the B site
  // kx2 = second creation momentum along x for the B site
  // ky2 = second creation momentum along y for the B site
  // kx3 = first annihilation momentum along x for the B site
  // ky3 = first annihilation momentum along y for the B site
  // kx4 = second annihilation momentum along x for the B site
  // ky4 = second annihilation momentum along y for the B site
  // s1 = sublattice index for the first creation operator
  // s2 = sublattice index for the second annihilation operator
  //
  // return value = corresponding matrix element
  Complex ComputeEmbeddingForTwoBodyOperator(int s1, int s2, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);
    
  // compute the matrix element for on-site two body interaction involving sites on generic sublattic
  //
  // dRx = number of unit vector translations along x-direction from EncodeSublatticeIndex (translations back to unit cell)
  // dRy = number of unit vector translations along y-direction from EncodeSublatticeIndex (translations back to unit cell)
  // kx2 = second creation momentum along x for the translated site
  // ky2 = second creation momentum along y for the translated site
  // kx3 = first annihilation momentum along x for the translated site
  // ky3 = first annihilation momentum along y for the translated site
  //
  // return value = corresponding matrix element
  Complex ComputeBlochPhases(int dRx, int dRy, int kx2, int ky2, int kx3, int ky3);

  // compute the transformation basis contribution to the interaction matrix element
  // 
  // oneBodyBasis = array of transformation basis matrices
  // momentumIndex1 = compact momentum index of the first creation operator
  // momentumIndex2 = compact momentum index of the second creation operator
  // momentumIndex3 = compact momentum index of the first annihilation operator
  // momentumIndex4 = compact momentum index of the second annihiliation operator
  // energyIndex1 = energy index of the first creation operator
  // energyIndex2 = energy index of the second creation operator
  // energyIndex3 = energy index of the first annihilation operator
  // energyIndex4 = energy index of the second annihiliation operator
  // siteIndex1 = site index of the first creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  // siteIndex2 = site index of the second creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  // siteIndex3 = site index of the first annihilation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  // siteIndex4 = site index of the second annihiliation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  
  virtual Complex ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
							int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
							int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
							int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4);

};


// compute the transformation basis contribution to the interaction matrix element
// 
// oneBodyBasis = array of transformation basis matrices
// momentumIndex1 = compact momentum index of the first creation operator
// momentumIndex2 = compact momentum index of the second creation operator
// momentumIndex3 = compact momentum index of the first annihilation operator
// momentumIndex4 = compact momentum index of the second annihiliation operator
// energyIndex1 = energy index of the first creation operator
// energyIndex2 = energy index of the second creation operator
// energyIndex3 = energy index of the first annihilation operator
// energyIndex4 = energy index of the second annihiliation operator
// siteIndex1 = site index of the first creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
// siteIndex2 = site index of the second creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
// siteIndex3 = site index of the first annihilation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
// siteIndex4 = site index of the second annihiliation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)

inline Complex ParticleOnLatticeHofstadterSingleBandGenericHamiltonian::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
													      int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
													      int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
													      int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1]) * Conj(oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3] * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4]);
}


#endif
