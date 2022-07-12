////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 3d topological insulator based on the pryochlore       //
//                       model and restricted to four bands                   //
//                                                                            //
//                        last modification : 13/08/2012                      //
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


#ifndef PARTICLEONCUBICLATTICEFOURBANDPYROCHLOREHAMILTONIAN_H
#define PARTICLEONCUBICLATTICEFOURBANDPYROCHLOREHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian.h"
#include "Tools/FTITightBinding/Abstract3DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnCubicLatticeFourBandPyrochloreHamiltonian : public ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian
{

 protected:
 
  // number of sites in the z direction
  int NbrSiteZ;
  // number of sites in the direction perpendicular to X
  int NbrSiteYZ;

  // repulsive on-site potential strength between identical spins
  double UPotential;
  // repulsive on-site potential strength between opposite spins
  double VPotential;
  // repulsive nearest neighbor potential strength between identical spins
  double WUPotential;
  // repulsive nearest neighbor potential strength between opposite spins
  double WVPotential;

  // pointer to the tight binding model
  Abstract3DTightBindingModel* TightBindingModel;

  // use flat band model
  bool FlatBand;


 public:

  // default constructor
  //
  ParticleOnCubicLatticeFourBandPyrochloreHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // uPotential = repulsive on-site potential strength between identical spins
  // vPotential = repulsive on-site potential strength between opposite spins
  // wuPotential = repulsive nearest neighbor strength between identical spins
  // wvPotential = repulsive nearest neighbor strength between opposite spins
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCubicLatticeFourBandPyrochloreHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, int nbrSiteZ, double uPotential, double vPotential, double wuPotential, double wvPotential, Abstract3DTightBindingModel* tightBindingModel, 
						      bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnCubicLatticeFourBandPyrochloreHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();
  
  // compute the transformation basis contribution to the interaction matrix element
  // 
  // oneBodyBasis = array of transformation basis matrices
  // momentumIndex3 = compact momentum index of the first creation operator
  // momentumIndex4 = compact momentum index of the second creation operator
  // momentumIndex1 = compact momentum index of the first annihilation operator
  // momentumIndex2 = compact momentum index of the second annihiliation operator
  // energyIndex3 = energy index of the first creation operator
  // energyIndex4 = energy index of the second creation operator
  // energyIndex1 = energy index of the first annihilation operator
  // energyIndex2 = energy index of the second annihiliation operator
  // siteIndex3 = site index of the first creation operator (0,...,3 = up, 4,...,7 = down)
  // siteIndex4 = site index of the second creation operator (0,...,3 = up, 4,...,7 = down)
  // siteIndex1 = site index of the first annihilation operator (0,...,3 = up, 4,...,7 = down)
  // siteIndex2 = site index of the second annihiliation operator (0,...,3 = up, 4,...,7 = down)
  inline Complex ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
                                                       int momentumIndex3, int momentumIndex4, int momentumIndex1, int momentumIndex2,
                                                       int energyIndex3, int energyIndex4, int energyIndex1, int energyIndex2,
                                                       int siteIndex3, int siteIndex4, int siteIndex1, int siteIndex2);
  
  // compute the contribution to the onsite interaction matrix elements with the same spin
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
  // factorU = repulsive on-site potential strength between identical spins
  //return value = corresponding matrix element
  Complex ComputeOnSiteContributionSameSpin(ComplexMatrix* oneBodyBasis,
						       int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
						       int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,double factorU);
  
  // compute the contribution to the onsite interaction matrix elements with opposite spin
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
  // factorV = repulsive on-site potential strength between opposite spins
  // epsilon = sign of the permutation +1 for bosons, -1 for fermions
  //return value = corresponding matrix element
  Complex ComputeOnSiteContributionOppositeSpin(ComplexMatrix* oneBodyBasis,
						       int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
						       int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,double factorV, int epsilon);
  
  
  // compute the contribution to the nearest neighbor interaction 
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
  // factorWU = repulsive nearest neighbor potential strength between identical spins
  // factorWV = repulsive nearest neighbor potential strength between opposite spins
  // epsilon = sign of the permutation +1 for bosons, -1 for fermions
  //return value = corresponding matrix element
  Complex ComputeNearestNeighborInteractionContribution(ComplexMatrix* oneBodyBasis,
						       int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
						       int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4, double factorWU, double factorWV, int epsilon);

};


// compute the transformation basis contribution to the interaction matrix element
// 
// oneBodyBasis = array of transformation basis matrices
// momentumIndex3 = compact momentum index of the first creation operator
// momentumIndex4 = compact momentum index of the second creation operator
// momentumIndex1 = compact momentum index of the first annihilation operator
// momentumIndex2 = compact momentum index of the second annihiliation operator
// energyIndex3 = energy index of the first creation operator
// energyIndex4 = energy index of the second creation operator
// energyIndex1 = energy index of the first annihilation operator
// energyIndex2 = energy index of the second annihiliation operator
// siteIndex3 = site index of the first creation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex4 = site index of the second creation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex1 = site index of the first annihilation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex2 = site index of the second annihiliation operator (0,...,3 = up, 4,...,7 = down)

inline Complex ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
													 int momentumIndex3, int momentumIndex4, int momentumIndex1, int momentumIndex2,
													 int energyIndex3, int energyIndex4, int energyIndex1, int energyIndex2,
													 int siteIndex3, int siteIndex4, int siteIndex1, int siteIndex2)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1]) * Conj(oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3] * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4]);
}


#endif
