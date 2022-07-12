////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//          class of a generic density-density two body interaction           //
//                       projected onto two bands and                         //
//              assuming a Bloch form for the tight binding model             //
//                                                                            //
//                        last modification : 26/09/2014                      //
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


#ifndef PARTICLEONLATTICEDENSITYDENSITYINTERACTIONTWOBANDHAMILTONIAN_H
#define PARTICLEONLATTICEDENSITYDENSITYINTERACTIONTWOBANDHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian : public ParticleOnLatticeQuantumSpinHallFullTwoBandHamiltonian
{

 protected:
 
  // number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
  int* NbrInteractingOrbitals;
  // orbital indices of the orbitals interacting with each orbital within the unit cell 
  // at the origin through a density-density term  
  int** InteractingOrbitalsOrbitalIndices;  
  // spatial indices (sorted as 2 consecutive integers) of the orbitals interacting 
  // with each orbital within the unit cell at the origin through a density-density term
  int** InteractingOrbitalsSpatialIndices;
  // intensity of each density-density term 
  double** InteractingOrbitalsPotentials; 
  
  // index of the first band to be filled
  int BandIndex1;
  // index of the second band to be filled
  int BandIndex2;
  
  // numerical factor for momentum along x
  double KxFactor;
  // numerical factor for momentum along y
  double KyFactor;
  
  // use flat band model
  bool FlatBand;
  //  gap between the first band and the second band when using the flat band model   
  double FlatBandOneBodyGap;

 public:

  // default constructor
  //
  ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // bandIndex1 = index of the first band onto which the Hamiltonian is projected 
  // bandIndex2 = index of the second band onto which the Hamiltonian is projected 
  // nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
  // interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
  // interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
  // interactingOrbitalsPotentials = intensity of each density-density term 
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // flatBandOneBodyGap = set the gap between the first band and the second band when using the flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY,
								      int bandIndex1, int bandIndex2, 
								      int* nbrInteractingOrbitals, int** interactingOrbitalsOrbitalIndices,
								      int** interactingOrbitalsSpatialIndices, double** interactingOrbitalsPotentials,
								      Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, double flatBandOneBodyGap, 
								      AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian();
  

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
  // siteIndex1 = site index of the first creation operator
  // siteIndex2 = site index of the second creation operator
  // siteIndex3 = site index of the first annihilation operator
  // siteIndex4 = site index of the second annihiliation operator 
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
// siteIndex1 = site index of the first creation operator
// siteIndex2 = site index of the second creation operator
// siteIndex3 = site index of the first annihilation operator
// siteIndex4 = site index of the second annihiliation operator 

inline Complex ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
												       int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
												       int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
												       int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1] * oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3] * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4]);
}


#endif
