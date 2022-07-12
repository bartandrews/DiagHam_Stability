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
//                                                                            //
//                        last modification : 04/10/2012                      //
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


#ifndef PARTICLEONLATTICEKAGOMELATTICESINGLEBANDHAMILTONIAN_H
#define PARTICLEONLATTICEKAGOMELATTICESINGLEBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{

protected:

    // index of the band that has to be partially filled
    int BandIndex;

    // nearest neighbor density-density potential strength (or on-site density-density potential for bosons)
    double UPotential;
    // second nearest neighbor density-density potential strength (or nearest neighbor density-density potential for bosons)
    double VPotential;
    // third nearest neighbor density-density potential strength (or next nearest neighbor density-density potential for bosons)
    double WPotential;

    // use flat band model
    bool FlatBand;

public:

    // default constructor
    //
    ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian();

    // constructor
    //
    // particles = Hilbert space associated to the system
    // nbrParticles = number of particles
    // nbrSiteX = number of sites in the x direction
    // nbrSiteY = number of sites in the y direction
    // tightBindingModel = pointer to the tight binding model
    // uPotential = strength of the repulsive two body neareast neighbor interaction
    // vPotential = strength of the repulsive two body second nearest neighbor interactio
    // wPotential = strength of the repulsive two body third nearest neighbor interaction (or next nearest neighbor density-density potential for bosons)
    // gammaX = boundary condition twisting angle along x
    // gammaY = boundary condition twisting angle along y
    // bandIndex = index of the band that has to be partially filled
    // flatBandFlag = use flat band model
    // architecture = architecture to use for precalculation
    // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
    ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
								   Abstract2DTightBindingModel* tightBindingModel,
								   double uPotential, double vPotential, double wPotential, int bandIndex, bool flatBandFlag, 
								   AbstractArchitecture* architecture, long memory = -1);

    // destructor
    //
    ~ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian();

protected:

    // evaluate all interaction factors
    //
    virtual void EvaluateInteractionFactors();

    // compute the matrix element for the two body interaction between two sites A and B
    //
    // kx2 = annihilation momentum along x for the B site
    // ky2 = annihilation momentum along y for the B site
    // kx4 = creation momentum along x for the B site
    // ky4 = creation momentum along y for the B site
    // return value = corresponding matrix element
    Complex ComputeTwoBodyMatrixElementNNAB(int kx2, int ky2, int kx4, int ky4);

    // compute the matrix element for the two body interaction between two sites B and C
    //
    // kx2 = annihilation momentum along x for the C site
    // ky2 = annihilation momentum along y for the C site
    // kx4 = creation momentum along x for the C site
    // ky4 = creation momentum along y for the C site
    // return value = corresponding matrix element

    Complex ComputeTwoBodyMatrixElementNNBC(int kx2, int ky2, int kx4, int ky4);

    // compute the matrix element for the two body interaction between two sites C and A
    //
    // kx2 = annihilation momentum along x for the A site
    // ky2 = annihilation momentum along y for the A site
    // kx4 = creation momentum along x for the A site
    // ky4 = creation momentum along y for the A site
    // return value = corresponding matrix element
    Complex ComputeTwoBodyMatrixElementNNCA(int kx2, int ky2, int kx4, int ky4);

    // compute the matrix element for the two body interaction between two A sites
    //
    // kx2 = annihilation momentum along x for the second site
    // ky2 = annihilation momentum along y for the second site
    // kx4 = creation momentum along x for the second site
    // ky4 = creation momentum along y for the second site
    // return value = corresponding matrix element
    Complex ComputeTwoBodyMatrixElementNNNAB(int kx2, int ky2, int kx4, int ky4);

    // compute the matrix element for the two body interaction between two B sites
    //
    // kx2 = annihilation momentum along x for the second site
    // ky2 = annihilation momentum along y for the second site
    // kx4 = creation momentum along x for the second site
    // ky4 = creation momentum along y for the second site
    // return value = corresponding matrix element
    Complex ComputeTwoBodyMatrixElementNNNBC(int kx2, int ky2, int kx4, int ky4);

    // compute the matrix element for the two body interaction between two C sites
    //
    // kx2 = annihilation momentum along x for the second site
    // ky2 = annihilation momentum along y for the second site
    // kx4 = creation momentum along x for the second site
    // ky4 = creation momentum along y for the second site
    // return value = corresponding matrix element
    Complex ComputeTwoBodyMatrixElementNNNCA(int kx2, int ky2, int kx4, int ky4);

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

    // compute the matrix element for on-site two body interaction involving C sites
    //
    // kx1 = first creation momentum along x for the C site
    // ky1 = first creation momentum along y for the C site
    // kx2 = second creation momentum along x for the C site
    // ky2 = second creation momentum along y for the C site
    // kx3 = first annihilation momentum along x for the C site
    // ky3 = first annihilation momentum along y for the C site
    // kx4 = second annihilation momentum along x for the C site
    // ky4 = second annihilation momentum along y for the C site
    // return value = corresponding matrix element
    Complex ComputeTwoBodyMatrixElementOnSiteCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

};

#endif
