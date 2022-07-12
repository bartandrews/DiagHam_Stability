////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                          class author: Yang-Le Wu                          //
//                                                                            //
//    class of face-centered cubic lattice model with interacting particles   //
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


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeFCCLatticeSingleBandHamiltonian.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


// default constructor
//

ParticleOnLatticeFCCLatticeSingleBandHamiltonian::ParticleOnLatticeFCCLatticeSingleBandHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// uPotential = repulsive on-site potential strength between like orbitals
// vPotential = repulsive potential strength between NN sites
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeFCCLatticeSingleBandHamiltonian::ParticleOnLatticeFCCLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles,
        int nbrSiteX, int nbrSiteY, int nbrSiteZ, double uPotential, double vPotential,
        Abstract3DTightBindingModel* tightBindingModel, int bandIndex, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
    this->Particles = particles;
    this->NbrParticles = nbrParticles;
    this->NbrSiteX = nbrSiteX;
    this->NbrSiteY = nbrSiteY;
    this->NbrSiteZ = nbrSiteZ;
    this->NbrSiteYZ = this->NbrSiteY * this->NbrSiteZ;
    this->LzMax = nbrSiteX * nbrSiteY * nbrSiteZ - 1;
    this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
    this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
    this->KzFactor = 2.0 * M_PI / ((double) this->NbrSiteZ);

    this->HamiltonianShift = 0.0;
    this->TightBindingModel = tightBindingModel;
    this->BandIndex = bandIndex;
    this->FlatBand = flatBandFlag;

    this->UPotential = uPotential;
    this->VPotential = vPotential;

    this->Architecture = architecture;
    this->Memory = memory;
    this->OneBodyInteractionFactors = 0;
    this->FastMultiplicationFlag = false;
    this->HermitianSymmetryFlag = true;

    long MinIndex;
    long MaxIndex;
    this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
    this->PrecalculationShift = (int) MinIndex;
    this->EvaluateInteractionFactors();

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

ParticleOnLatticeFCCLatticeSingleBandHamiltonian::~ParticleOnLatticeFCCLatticeSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//

void ParticleOnLatticeFCCLatticeSingleBandHamiltonian::EvaluateInteractionFactors()
{
    long TotalNbrInteractionFactors = 0;
    int NbrStatePerBand = this->TightBindingModel->GetNbrStatePerBand();
    double FactorU = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ));
    double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ));

    ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
    if (this->FlatBand == false)
    {
        this->OneBodyInteractionFactors = new double [this->TightBindingModel->GetNbrStatePerBand()];
        for (int k = 0; k < NbrStatePerBand; ++k)
            this->OneBodyInteractionFactors[k] = 0.0;
    }

    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
            for (int kz = 0; kz < this->NbrSiteZ; ++kz)
            {
                int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky, kz);
                if (this->FlatBand == false) // FIXME: Why "0.5 *" ?
                    this->OneBodyInteractionFactors[Index] = 0.5 * this->TightBindingModel->GetEnergy(this->BandIndex, Index);
                OneBodyBasis[Index] = this->TightBindingModel->GetOneBodyMatrix(Index);
            }

    if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
        this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
        this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
            this->NbrSectorIndicesPerSum[i] = 0;
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                        for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
                            for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)
                            {
                                int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1, kz1);
                                int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2, kz2);
                                if (Index1 < Index2)
                                    ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2, kz1 + kz2)];
                            }
        this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
        {
            if (this->NbrSectorIndicesPerSum[i]  > 0)
            {
                this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];
                this->NbrSectorIndicesPerSum[i] = 0;
            }
        }
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                        for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
                            for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)
                            {
                                int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1, kz1);
                                int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2, kz2);
                                if (Index1 < Index2)
                                {
                                    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2, kz1 + kz2);
                                    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
                                    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
                                    ++this->NbrSectorIndicesPerSum[TmpSum];
                                }
                            }

        this->InteractionFactors = new Complex* [this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
        {
            this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
            int Index = 0;
            for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
            {
                int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
                int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
                int kx1, ky1, kz1, kx2, ky2, kz2;
                this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1, kz1);
                this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2, kz2);
                for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
                {
                    int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                    int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                    int kx3, ky3, kz3, kx4, ky4, kz4;
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3, kz3);
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4, kz4);

                    // the InteractionFactors are coefficients to A+_3 A_1 A+_4 A_2
                    // tricky part: OneBodyBasis[Index] stores the result of LapackDiagonalize
                    // and its [0][_] elements are the COMPLEX CONJUGATE of wave functions < _ |lower band>. (See the end of HermitianMatrix.cc)

                    // NN
                    Complex SumV = 0.0;
                    if (this->VPotential != 0.0)
                    {
                        SumV += this->ComputeTwoBodyMatrixElementNN(OneBodyBasis, Index1, Index2, Index3, Index4);
                        SumV -= this->ComputeTwoBodyMatrixElementNN(OneBodyBasis, Index2, Index1, Index3, Index4);
                        SumV -= this->ComputeTwoBodyMatrixElementNN(OneBodyBasis, Index1, Index2, Index4, Index3);
                        SumV += this->ComputeTwoBodyMatrixElementNN(OneBodyBasis, Index2, Index1, Index4, Index3);
                        SumV *= FactorV;
                    }

                    this->InteractionFactors[i][Index] = SumV;
                    this->InteractionFactors[i][Index] *= -1.0;

                    TotalNbrInteractionFactors++;
                    ++Index;
                }
            }
        }
    }
    else
    {
        this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
        this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
            this->NbrSectorIndicesPerSum[i] = 0;
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                        for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
                            for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)
                            {
                                int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1, kz1);
                                int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2, kz2);
                                if (Index1 <= Index2)
                                    ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2, kz1 + kz2)];
                            }
        this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
        {
            if (this->NbrSectorIndicesPerSum[i]  > 0)
            {
                this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];
                this->NbrSectorIndicesPerSum[i] = 0;
            }
        }
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                        for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
                            for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)
                            {
                                int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1, kz1);
                                int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2, kz2);
                                if (Index1 <= Index2)
                                {
                                    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2, kz1 + kz2);
                                    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
                                    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
                                    ++this->NbrSectorIndicesPerSum[TmpSum];
                                }
                            }

        this->InteractionFactors = new Complex* [this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
        {
            this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
            int Index = 0;
            for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
            {
                int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
                int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
                int kx1, ky1, kz1, kx2, ky2, kz2;
                this->TightBindingModel->GetLinearizedMomentumIndex(Index1, kx1, ky1, kz1);
                this->TightBindingModel->GetLinearizedMomentumIndex(Index2, kx2, ky2, kz2);
                for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
                {
                    int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                    int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                    int kx3, ky3, kz3, kx4, ky4, kz4;
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3, kz3);
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4, kz4);

                    // the InteractionFactors are coefficients to A+_3 A_1 A+_4 A_2
                    // tricky part: OneBodyBasis[Index] stores the result of LapackDiagonalize
                    // and its [0][_] elements are the COMPLEX CONJUGATE of wave functions < _ |lower band>. (See the end of HermitianMatrix.cc)

                    // on-site U
                    Complex SumU = 0.0;
                    SumU += this->ComputeTwoBodyMatrixElementOnSiteU(OneBodyBasis, Index1, Index2, Index3, Index4);
                    SumU += this->ComputeTwoBodyMatrixElementOnSiteU(OneBodyBasis, Index2, Index1, Index3, Index4);
                    SumU += this->ComputeTwoBodyMatrixElementOnSiteU(OneBodyBasis, Index1, Index2, Index4, Index3);
                    SumU += this->ComputeTwoBodyMatrixElementOnSiteU(OneBodyBasis, Index2, Index1, Index4, Index3);
                    SumU *= FactorU;

                    // NN
                    Complex SumV = 0.0;
                    if (this->VPotential != 0.0)
                    {
                        SumV += this->ComputeTwoBodyMatrixElementNN(OneBodyBasis, Index1, Index2, Index3, Index4);
                        SumV += this->ComputeTwoBodyMatrixElementNN(OneBodyBasis, Index2, Index1, Index3, Index4);
                        SumV += this->ComputeTwoBodyMatrixElementNN(OneBodyBasis, Index1, Index2, Index4, Index3);
                        SumV += this->ComputeTwoBodyMatrixElementNN(OneBodyBasis, Index2, Index1, Index4, Index3);
                        SumV *= FactorV;
                    }


                    this->InteractionFactors[i][Index] = SumU + SumV;
                    if (Index3 == Index4)
                        this->InteractionFactors[i][Index] *= 0.5;
                    if (Index1 == Index2)
                        this->InteractionFactors[i][Index] *= 0.5;

                    TotalNbrInteractionFactors++;
                    ++Index;
                }
            }
        }
    }
    cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
    cout << "====================================" << endl;

    delete [] OneBodyBasis;
}
