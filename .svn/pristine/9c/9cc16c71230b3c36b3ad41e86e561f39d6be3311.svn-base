////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                class of quatum spin Hall restricted to two bands           //
//           using two decoupled kagome models and three body interaction     //
//                                                                            //
//                        last modification : 26/10/2013                      //
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


#ifndef PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDDECOUPLEDKAGOMETHREEBODYHAMILTONIAN_H
#define PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDDECOUPLEDKAGOMETHREEBODYHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian : public ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian
{

 protected:
  
  // strength of the repulsive on site three body interaction between identical spins
  double UPotential;
  // strength of the repulsive on site three body interaction between opposite spins
  double VPotential;

  // strength of the repulsive on site two body interaction between identical spins
  double TwoBodyUPotential;
  // strength of the repulsive on site two body interaction between opposite spins
  double TwoBodyVPotential;

  // use flat band model
  bool FlatBand;
  // use two copies of Kagome with time reversal symmetry
  bool TimeReversal;
  // pointer to the tight binding model of the first copy (must be in Bloch form)
  Abstract2DTightBindingModel* TightBindingModel;

 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // vPotential = strength of the repulsive on site two body interaction between opposite spins
  // tightBindingModel = pointer to the tight binding model of the first copy
  // flatBandFlag = use flat band model
  // timeReversalFlag = apply thge time reversal symmetry on the second copy of the tight binding model  (must be in Bloch form)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
									     double uPotential, double vPotential, 
									     Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, bool timeReversalFlag, 
									     AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the one body transformation matrices and the optional one body band stucture contribution
  //
  // oneBodyBasis = array of one body transformation matrices (the leftmost upper block for the spin up, the rightmsdt lower block for the spin down)
  void ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis);

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
  Complex ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
						int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
						int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
						int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4);
  
  // compute the transformation basis contribution to the interaction matrix element
  // 
  // oneBodyBasis = array of transformation basis matrices
  // momentumIndex1 = compact momentum index of the first creation operator
  // momentumIndex2 = compact momentum index of the second creation operator
  // momentumIndex3 = compact momentum index of the third creation operator
  // momentumIndex4 = compact momentum index of the first annihilation operator
  // momentumIndex5 = compact momentum index of the second annihiliation operator
  // momentumIndex6 = compact momentum index of the third annihiliation operator
  // energyIndex1 = energy index of the first creation operator
  // energyIndex2 = energy index of the second creation operator
  // energyIndex3 = energy index of the third creation operator
  // energyIndex4 = energy index of the first annihilation operator
  // energyIndex5 = energy index of the second annihiliation operator
  // energyIndex6 = energy index of the third annihiliation operator
  // siteIndex1 = site index of the first creation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)
  // siteIndex2 = site index of the second creation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)
  // siteIndex3 = site index of the third creation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)
  // siteIndex4 = site index of the first annihilation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)
  // siteIndex5 = site index of the second annihiliation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)
  // siteIndex6 = site index of the third annihiliation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)
  Complex ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
						int momentumIndex1, int momentumIndex2, int momentumIndex3, 
						int momentumIndex4, int momentumIndex5, int momentumIndex6, 
						int energyIndex1, int energyIndex2, int energyIndex3, 
						int energyIndex4, int energyIndex5, int energyIndex6,
						int siteIndex1, int siteIndex2, int siteIndex3, 
						int siteIndex4, int siteIndex5, int siteIndex6);


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

inline Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
																 int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
																 int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
																 int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1]) * oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3] * Conj(oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4]);
}

// compute the transformation basis contribution to the interaction matrix element
// 
// oneBodyBasis = array of transformation basis matrices
// momentumIndex1 = compact momentum index of the first creation operator
// momentumIndex2 = compact momentum index of the second creation operator
// momentumIndex3 = compact momentum index of the third creation operator
// momentumIndex4 = compact momentum index of the first annihilation operator
// momentumIndex5 = compact momentum index of the second annihiliation operator
// momentumIndex6 = compact momentum index of the third annihiliation operator
// energyIndex1 = energy index of the first creation operator
// energyIndex2 = energy index of the second creation operator
// energyIndex3 = energy index of the third creation operator
// energyIndex4 = energy index of the first annihilation operator
// energyIndex5 = energy index of the second annihiliation operator
// energyIndex6 = energy index of the third annihiliation operator
// siteIndex1 = site index of the first creation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)
// siteIndex2 = site index of the second creation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)
// siteIndex3 = site index of the third creation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)
// siteIndex4 = site index of the first annihilation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)
// siteIndex5 = site index of the second annihiliation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)
// siteIndex6 = site index of the third annihiliation operator (0=Aup, 1=Bup, 2=Cup, 3=Adown, 4=Bdown, 5=Cdown)

inline Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeThreeBodyHamiltonian::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
																 int momentumIndex1, int momentumIndex2, int momentumIndex3, 
																 int momentumIndex4, int momentumIndex5, int momentumIndex6, 
																 int energyIndex1, int energyIndex2, int energyIndex3, 
																 int energyIndex4, int energyIndex5, int energyIndex6,
																 int siteIndex1, int siteIndex2, int siteIndex3, 
																 int siteIndex4, int siteIndex5, int siteIndex6)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1]) * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4] * 
	  Conj(oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex5][energyIndex5][siteIndex5] * 
	  Conj(oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3]) * oneBodyBasis[momentumIndex6][energyIndex6][siteIndex6]);
}


#endif
