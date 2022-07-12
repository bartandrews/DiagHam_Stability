////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                          class author: Yang-Le Wu                          //
//                                                                            //
//       class of three-orbital chiral model with interacting particles       //
//                       in the single band approximation                     //
//                                                                            //
//                        last modification : 09/09/2012                      //
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


#ifndef PARTICLEONCUBICLATTICESINGLEBANDHAMILTONIAN_H
#define PARTICLEONCUBICLATTICESINGLEBANDHAMILTONIAN_H

#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
#include "Tools/FTITightBinding/Abstract3DTightBindingModel.h"

class ParticleOnCubicLatticeSingleBandHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{

protected:

    // number of sites in the z direction
    int NbrSiteZ;
    // number of sites in the direction perpendicular to X
    int NbrSiteYZ;

    // numerical factor for momentum along z
    double KzFactor;

    // index of the band to be filled
    int BandIndex;

    // on-site density-density repulsion between like orbitals
    double UPotential;
    // on-site density-density repulsion between different orbitals
    double UABPotential;
    // nearest-neighbor density-density repulsion (between all orbitals)
    double VPotential;

    // pointer to the tight binding model
    Abstract3DTightBindingModel* TightBindingModel;

    // use flat band model
    bool FlatBand;

    // add self energy from non-normal ordering
    bool SelfEnergy;

    // add on-site intra-orbital interaction even for fermions
    bool OnSiteIntra;

public:

    // default constructor
    //
    ParticleOnCubicLatticeSingleBandHamiltonian();

    // constructor
    //
    // particles = Hilbert space associated to the system
    // nbrParticles = number of particles
    // nbrSiteX = number of sites in the x direction
    // nbrSiteY = number of sites in the y direction
    // nbrSiteZ = number of sites in the z direction
    // uPotential = repulsive on-site potential strength between like orbitals
    // uABPotential = repulsive on-site potential strength between different orbitals
    // vPotential = repulsive potential strength between NN sites
    // tightBindingModel = pointer to the tight binding model
    // bandIndex = index of the band to consider
    // flatBandFlag = use flat band model
    // selfEnergyFlag = add self energy from non-normal ordering
    // intraFlag = add on-site intra-orbital interaction even for fermions
    // architecture = architecture to use for precalculation
    // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
    ParticleOnCubicLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles,
            int nbrSiteX, int nbrSiteY, int nbrSiteZ, double uPotential, double uABPotential, double vPotential,
            Abstract3DTightBindingModel* tightBindingModel, int bandIndex, bool flatBandFlag, bool selfEnergyFlag, bool intraFlag, AbstractArchitecture* architecture, long memory = -1);

    // destructor
    //
    ~ParticleOnCubicLatticeSingleBandHamiltonian();

protected:

    // evaluate all interaction factors
    //
    virtual void EvaluateInteractionFactors();

    // compute the band-basis matrix element of the density-density interaction between like orbitals on the same site
    //
    // oneBodyBasis = array of transformation basis matrices
    // Index1 = compact momentum index of the first annihilation operator
    // Index2 = compact momentum index of the second annihilation operator
    // Index3 = compact momentum index of the first creation operator
    // Index4 = compact momentum index of the second creation operator
    virtual Complex ComputeTwoBodyMatrixElementOnSiteU(ComplexMatrix* OneBodyBasis, int Index1, int Index2, int Index3, int Index4);

    // compute the band-basis matrix element of the density-density interaction between different orbitals on the same site
    //
    // oneBodyBasis = array of transformation basis matrices
    // Index1 = compact momentum index of the first annihilation operator
    // Index2 = compact momentum index of the second annihilation operator
    // Index3 = compact momentum index of the first creation operator
    // Index4 = compact momentum index of the second creation operator
    virtual Complex ComputeTwoBodyMatrixElementOnSiteUAB(ComplexMatrix* OneBodyBasis, int Index1, int Index2, int Index3, int Index4);

    // compute the band-basis matrix element of the density-density interaction between orbitals on adjacent sites
    //
    // oneBodyBasis = array of transformation basis matrices
    // Index1 = compact momentum index of the first annihilation operator
    // Index2 = compact momentum index of the second annihilation operator
    // Index3 = compact momentum index of the first creation operator
    // Index4 = compact momentum index of the second creation operator
    virtual Complex ComputeTwoBodyMatrixElementNN(ComplexMatrix* OneBodyBasis, int Index1, int Index2, int Index3, int Index4);
};

// compute the band-basis matrix element of the density-density interaction between like orbitals on the same site
//
// oneBodyBasis = array of transformation basis matrices
// Index1 = compact momentum index of the first annihilation operator
// Index2 = compact momentum index of the second annihilation operator
// Index3 = compact momentum index of the first creation operator
// Index4 = compact momentum index of the second creation operator
inline Complex ParticleOnCubicLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteU(ComplexMatrix* OneBodyBasis, int Index1, int Index2, int Index3, int Index4)
{
    Complex Tmp = 0.0;
    for (int i = 0; i < 3; ++i)
        Tmp += Conj(OneBodyBasis[Index1][this->BandIndex][i]) * OneBodyBasis[Index3][this->BandIndex][i]
             * Conj(OneBodyBasis[Index2][this->BandIndex][i]) * OneBodyBasis[Index4][this->BandIndex][i];
    return Tmp;
}

// compute the band-basis matrix element of the density-density interaction between different orbitals on the same site
//
// oneBodyBasis = array of transformation basis matrices
// Index1 = compact momentum index of the first annihilation operator
// Index2 = compact momentum index of the second annihilation operator
// Index3 = compact momentum index of the first creation operator
// Index4 = compact momentum index of the second creation operator
inline Complex ParticleOnCubicLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteUAB(ComplexMatrix* OneBodyBasis, int Index1, int Index2, int Index3, int Index4)
{
    Complex Tmp = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = i + 1; j < 3; ++j)
            Tmp += Conj(OneBodyBasis[Index1][this->BandIndex][i]) * OneBodyBasis[Index3][this->BandIndex][i]
                 * Conj(OneBodyBasis[Index2][this->BandIndex][j]) * OneBodyBasis[Index4][this->BandIndex][j];
    return Tmp;
}

// compute the band-basis matrix element of the density-density interaction between orbitals on adjacent sites
//
// oneBodyBasis = array of transformation basis matrices
// Index1 = compact momentum index of the first annihilation operator
// Index2 = compact momentum index of the second annihilation creation operator
// Index3 = compact momentum index of the first creation operator
// Index4 = compact momentum index of the second creation operator
inline Complex ParticleOnCubicLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNN(ComplexMatrix* OneBodyBasis, int Index1, int Index2, int Index3, int Index4)
{
    int kx2, ky2, kz2;
    this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2, kz2);
    int kx4, ky4, kz4;
    this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4, kz4);

    Complex m = (Phase(this->KxFactor * (kx2 - kx4)) + Phase(this->KyFactor * (ky2 - ky4)) + Phase(this->KzFactor * (kz2 - kz4)));
    return (  Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
            + Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
            + Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2])
         * (  Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
            + Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
            + Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]) * m;
}


#endif
