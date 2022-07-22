////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//      class of checkerboard lattice model with interacting particles        //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 08/09/2011                      //
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



#ifndef PARTICLEONLATTICEHOFSTADTERSINGLEBANDHAMILTONIAN_H
#define PARTICLEONLATTICEHOFSTADTERSINGLEBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"

#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"


#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeHofstadterSingleBandHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{

 protected:
  
  
  // nearest neighbor density-density potential strength
  double UPotential;
  // second nearest neighbor density-density potential strength
  double VPotential;
  // 2nd nearest neighbor potential for bosonic interactions
  double WPotential;

  // bosonic onsite three body potential
  double ThreeBodyPotential;

  // index of the band to be filled
  int BandIndex;

  TightBindingModelHofstadterSquare* TightBindingModel;  // comment

  // use flat band model
  bool FlatBand;
  bool RSquaredInteraction;
  bool FourFoldSmoothInteraction;
  
 public:

  // constructor
  //
  ParticleOnLatticeHofstadterSingleBandHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrCellsX = number of sites in the x direction
  // nbrCellsY = number of sites in the y direction
  // bandIndex = index of band to consider
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // vPotential = strength of the repulsive two body second nearest neighbor interactio
  // tightBindingModel = pointer to the tight binding model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeHofstadterSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellsX, int nbrCellsY, int bandIndex, double uPotential, double vPotential, double wPotential, TightBindingModelHofstadterSquare* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);  // TightBindingModelHofstadterSquare

  // destructor
  //
  ~ParticleOnLatticeHofstadterSingleBandHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two sites B and C 
  //
  // subA = sublattice index of the first site
  // subB = sublattice index of the second site
  // kx1 = first creation momentum along x for the B site
  // ky1 = first creation momentum along y for the B site
  // kx2 = second creation momentum along x for the B site
  // ky2 = second creation momentum along y for the B site
  // kx3 = first annihilation momentum along x for the B site
  // ky3 = first annihilation momentum along y for the B site
  // kx4 = second annihilation momentum along x for the B site
  // ky4 = second annihilation momentum along y for the B site
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementGenericAB(int subA, int subB, int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b);

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
  virtual Complex ComputeTwoBodyMatrixElementOnSite(int subl, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

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
  // dRx = number of unit vector translations along x-direction
  // dRy = number of unit vector translations along y-direction
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

  //Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeMatrixElement(ComplexMatrix* oneBodyBasis, int k1x, int k1y, int k2x, int k2y,
//                                                                               int k3x, int k3y, int k4x, int k4y, int siteIndex);


  int getXMomentumFromIndex(int index);
  int getYMomentumFromIndex(int index);
  int getBandFromIndex(int index);
private:
  int tS(int index, int xtrans, int ytrans);

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

inline Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
												       int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
												       int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
												       int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1]) * Conj(oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3] * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4]);
}

inline int ParticleOnLatticeHofstadterSingleBandHamiltonian::getXMomentumFromIndex(int index)
{
  
  return (int) index / this->NbrSiteY;
}

inline int ParticleOnLatticeHofstadterSingleBandHamiltonian::getYMomentumFromIndex(int index)
{
  
  return (int) index % this->NbrSiteY;
}



#endif
