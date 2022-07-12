////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Yang-Le Wu                            //
//                                                                            //
//               class of Haldane model with interacting particles            //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 16/08/2011                      //
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


#ifndef PARTICLEONLATTICEWITHSPINHALDANEMODELTWOBANDDECOUPLEDTHREEBODYHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINHALDANEMODELTWOBANDDECOUPLEDTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian : public ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian
{

 protected:
  
  // hoping amplitude between neareast neighbor sites
  double NNHopping;
  // hoping amplitude between next neareast neighbor sites
  double NextNNHopping;
  // hoping amplitude between next to next neareast neighbor sites
  double NextNextNNHopping;
  // phase on nnn hopping
  double HaldanePhase;
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

  // anisotropy for spin interactions:
  double ThreeBodySpinAnisotropy;
  
  // precalculation tables for cosine and sine factors
  int NBodyValue;	
  Complex* XPhaseTable;
  Complex* YPhaseTable;
  Complex* XHalfPhaseTable;
  Complex* YHalfPhaseTable;
  int XPhaseTableShift;
  int YPhaseTableShift;

 public:

  // default constructor
  //
  ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // vPotential = strength of the repulsive two body second neareast neighbor interaction
  // wPotential = strength of the repulsive three body neareast neighbor interaction
  // sPotential = strength of the repulsive three body next-to-nearest neighbor interaction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // t3 = hoping amplitude between next to next neareast neighbor sites
  // phi =  Haldane phase on nnn hopping
  // mus = sublattice staggered chemical potential 
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
          double uPotential, double vPotential, double wPotential, double sPotential,
									    double t1, double t2, double t3, double phi, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1, bool hermitianFlag=false);

  // destructor
  //
  virtual ~ParticleOnLatticeWithSpinHaldaneModelTwoBandDecoupledThreeBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the one body transformation matrices and the optional one body band stucture contribution
  //
  // oneBodyBasis = array of one body transformation matrices (spin up, spin down)
  virtual void ComputeOneBodyMatrices(ComplexMatrix** oneBodyBasis);

  // compute all the phase precalculation arrays 
  //
  virtual void ComputePhaseArray();

  // compute the matrix element for the two body interaction between two sites A and B 
  //
  // kx1 = annihilation momentum along x for the B site
  // ky1 = annihilation momentum along y for the B site
  // kx2 = creation momentum along x for the B site
  // ky2 = creation momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the two body interaction between two A sites (or two B sites) 
  //
  // kx1 = annihilation momentum along x for the second site
  // ky1 = annihilation momentum along y for the second site
  // kx2 = creation momentum along x for the second site
  // ky2 = creation momentum along y for the second site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAA(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the three body interaction between one site A and two sites B 
  //
  // kx2 = annihilation momentum along x for the first B site
  // ky2 = annihilation momentum along y for the first B site
  // kx3 = annihilation momentum along x for the second B site
  // ky3 = annihilation momentum along y for the second B site
  // kx5 = creation momentum along x for the first B site
  // ky5 = creation momentum along y for the first B site
  // kx6 = creation momentum along x for the second B site
  // ky6 = creation momentum along y for the second B site
  // return value = corresponding matrix element
  Complex ComputeThreeBodyMatrixElementABB(int kx2, int ky2, int kx3, int ky3, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the three body interaction between NNN sites 
  //
  // kx2 = annihilation momentum along x for the second site
  // ky2 = annihilation momentum along y for the second site
  // kx3 = annihilation momentum along x for the third site
  // ky3 = annihilation momentum along y for the third site
  // kx5 = creation momentum along x for the second site
  // ky5 = creation momentum along y for the second site
  // kx6 = creation momentum along x for the third site
  // ky6 = creation momentum along y for the third site
  // return value = corresponding matrix element
  Complex ComputeThreeBodyMatrixElementAAA(int kx2, int ky2, int kx3, int ky3, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for on-site two body interaction involving A sites
  //
  // return value = corresponding matrix element

  Complex ComputeTwoBodyMatrixElementOnSiteAA();


  // compute the matrix element for on-site two body interaction involving B sites
  //
  // kx1 = first creation momentum along x for the B site
  // ky1 = first creation momentum along y for the B site
  // kx2 = second creation momentum along x for the B site
  // ky2 = second creation momentum along y for the B site
  // kx3 = first annihilation momentum along x for the B site
  // ky3 = first annihilation momentum along y for the B site
  // kx4 = second annihilation momentum along x for the B site
  // ky4 = second annihilation momentum along y for the B site
  // return value = corresponding matrix element

  Complex ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the on-site three body interaction related to sites A
  //
  // kx1 = first creation momentum along x for the first A site
  // ky1 = first creation momentum along y for the first A site
  // kx2 = second creation momentum along x for the second A site
  // ky2 = second creation momentum along y for the second A site
  // kx3 = third creation momentum along x for the second A site
  // ky3 = third creation momentum along y for the second A site
  // kx4 = first annihilation momentum along x for the first A site
  // ky4 = first annihilation momentum along y for the first A site
  // kx5 = second annihilation momentum along x for the second A site
  // ky5 = second annihilation momentum along y for the second A site
  // kx6 = third annihilation momentum along x for the second A site
  // ky6 = third annihilation momentum along y for the second A site
  // return value = corresponding matrix element
  Complex ComputeThreeBodyMatrixElementOnSiteAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);

  // compute the matrix element for the on-site three body interaction related to sites A
  //
  // kx1 = first creation momentum along x for the first A site
  // ky1 = first creation momentum along y for the first A site
  // kx2 = second creation momentum along x for the second A site
  // ky2 = second creation momentum along y for the second A site
  // kx3 = third creation momentum along x for the second A site
  // ky3 = third creation momentum along y for the second A site
  // kx4 = first annihilation momentum along x for the first A site
  // ky4 = first annihilation momentum along y for the first A site
  // kx5 = second annihilation momentum along x for the second A site
  // ky5 = second annihilation momentum along y for the second A site
  // kx6 = third annihilation momentum along x for the second A site
  // ky6 = third annihilation momentum along y for the second A site
  // return value = corresponding matrix element
  Complex ComputeThreeBodyMatrixElementOnSiteBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6);


};



#endif
