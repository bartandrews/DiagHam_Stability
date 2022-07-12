////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Yang-Le Wu                            //
//                                                                            //
//            class of trianglular lattice model with two orbitals            //
//                       in the single band approximation                     // 
//                           with interpolation to LLL                        //
//                                                                            //
//                        last modification : 14/10/2012                      //
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


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeChern3TwoOrbitalTriangularLatticeLLLInterpolateSingleBandHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Polynomial/SpecialPolynomial.h"
#include "GeneralTools/StringTools.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// tightBindingModel = pointer to the tight binding model
// uPotential = strength of the repulsive two body nearest neighbor interaction
// flatBandFlag = use flat band model
// interpolationToFQH = the interpolation parameter lambda between FQH and FCI, as in H = lambda * H_FQH + (1 - lambda) * H_FCI
// twistAngle = angle between the two fundamental cycles of the torus in Radians
// aspectRatio = aspect ratio of torus, Lx / Ly
// nbrPseudopotentials = array of the number of pseudo-potentials
// pseudoPotentials = array of the pseudo-potentials
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeChern3TwoOrbitalTriangularLatticeLLLInterpolateSingleBandHamiltonian::ParticleOnLatticeChern3TwoOrbitalTriangularLatticeLLLInterpolateSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY,
        Abstract2DTightBindingModel* tightBindingModel, double uPotential, bool flatBandFlag,
        double interpolationToFQH, double twistAngle, double aspectRatio, int nbrPseudopotentials, double* pseudoPotentials,
        AbstractArchitecture* architecture, long memory)
{
    this->Particles = particles;
    this->NbrParticles = nbrParticles;
    this->NbrSiteX = nbrSiteX;
    this->NbrSiteY = nbrSiteY;
    this->LzMax = nbrSiteX * nbrSiteY - 1;
    this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
    this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
    this->TightBindingModel = tightBindingModel;
    this->UPotential = uPotential;
    this->BandIndex = 0;
    this->FlatBand = flatBandFlag;

    this->InterpolationToFQH = interpolationToFQH;
    this->TwistAngle = twistAngle;
    this->AspectRatio = aspectRatio;
    this->LLLGammaX = 0.0;
    this->LLLGammaY = 0.0;
    this->NbrPseudopotentials = nbrPseudopotentials;
    if (this->NbrPseudopotentials > 0)
    {
        this->Pseudopotentials = new double[this->NbrPseudopotentials];
        this->LaguerrePolynomials = new Polynomial[this->NbrPseudopotentials];
        for (int i = 0; i < this->NbrPseudopotentials; ++i)
        {
            this->Pseudopotentials[i] = pseudoPotentials[i];
            this->LaguerrePolynomials[i] = LaguerrePolynomial(i);
        }
    }
    else
    {
        this->Pseudopotentials = NULL;
        this->LaguerrePolynomials = NULL;
    }

    this->HamiltonianShift = 0.0;
    this->Architecture = architecture;
    this->Memory = memory;
    this->OneBodyInteractionFactors = 0;
    this->FastMultiplicationFlag = false;
    long MinIndex;
    long MaxIndex;
    this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
    this->PrecalculationShift = (int) MinIndex;  

    this->NbrColor = this->TightBindingModel->GetChernNumber(this->BandIndex);
    double gammaX = this->TightBindingModel->GetLLLGammaX(this->BandIndex);
    double gammaY = this->TightBindingModel->GetLLLGammaY(this->BandIndex);
    ComplexMatrix gauge(this->NbrSiteX, this->NbrSiteY, true);
    if (this->TightBindingModel->BuildGeneralizedLandauGauge(this->BandIndex, gauge))
    {
        cout << "failed to build the generalized Landau gauge!" << endl;
        exit(1);
    }

    if (this->InterpolationToFQH != 1.0)
    {
        this->EvaluateInteractionFactors();
        this->TransformToLLLGauge(gauge, gammaX, gammaY);
    }
    if (this->InterpolationToFQH != 0.0)
        this->EvaluateFQHInteractionFactors();

    if (memory > 0)
    {
        long TmpMemory = this->FastMultiplicationMemory(memory);
        cout << "fast = ";
        PrintMemorySize(cout, TmpMemory)<< endl;
        this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeChern3TwoOrbitalTriangularLatticeLLLInterpolateSingleBandHamiltonian::~ParticleOnLatticeChern3TwoOrbitalTriangularLatticeLLLInterpolateSingleBandHamiltonian()
{
}

