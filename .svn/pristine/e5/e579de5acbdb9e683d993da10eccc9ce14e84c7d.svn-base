////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//         class of 3d topological insulator based on the simple TI model     //
//                          and full four band support                        //
//                                                                            //
//                        last modification : 28/09/2012                      //
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


#ifndef PARTICLEONCUBICLATTICEFULLFOURBANDSIMPLETIHAMILTONIAN_H
#define PARTICLEONCUBICLATTICEFULLFOURBANDSIMPLETIHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian.h"
#include "Tools/FTITightBinding/Abstract3DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnCubicLatticeFullFourBandSimpleTIHamiltonian : public ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian
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

  // pointer to the tight binding model
  Abstract3DTightBindingModel* TightBindingModel;

  // use flat band model
  bool FlatBand;


 public:

  // default constructor
  //
  ParticleOnCubicLatticeFullFourBandSimpleTIHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // uPotential = repulsive on-site potential strength between identical spins
  // vPotential = repulsive on-site potential strength between opposite spins
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCubicLatticeFullFourBandSimpleTIHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, int nbrSiteZ, double uPotential, double vPotential,  Abstract3DTightBindingModel* tightBindingModel, 
							bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnCubicLatticeFullFourBandSimpleTIHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();
  
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
  // siteIndex1 = site index of the first creation operator (0,...,3 = up, 4,...,7 = down)
  // siteIndex2 = site index of the second creation operator (0,...,3 = up, 4,...,7 = down)
  // siteIndex3 = site index of the first annihilation operator (0,...,3 = up, 4,...,7 = down)
  // siteIndex4 = site index of the second annihiliation operator (0,...,3 = up, 4,...,7 = down)
  inline Complex ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
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
// siteIndex1 = site index of the first creation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex2 = site index of the second creation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex3 = site index of the first annihilation operator (0,...,3 = up, 4,...,7 = down)
// siteIndex4 = site index of the second annihiliation operator (0,...,3 = up, 4,...,7 = down)

inline Complex ParticleOnCubicLatticeFullFourBandSimpleTIHamiltonian::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
													    int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
													    int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
													    int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1]) * Conj(oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3] * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4]);
}


#endif
