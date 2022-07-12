////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of hamiltonian with particles on                   //
//               Chern insulator in the single band approximation             //
//                                                                            //
//                        last modification : 23/02/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
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



// default constructor
//

ParticleOnLatticeChernInsulatorSingleBandHamiltonian::ParticleOnLatticeChernInsulatorSingleBandHamiltonian()
{
  this->InterpolationToFQH = 0.0;
  this->NbrColor = 1;
  this->TwistAngle = M_PI / 2;
  this->AspectRatio = 1e5; // this default value should never be used
  this->LLLGammaX = 0.0;
  this->LLLGammaY = 0.0;
  this->NbrPseudopotentials = 0;
  this->Pseudopotentials = NULL;
  this->LaguerrePolynomials = NULL;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// tightBindingModel = pointer to the tight binding model
// nbrColor = dimension of the internal (color) space. this should be equal to the Chern number of the occupid band
// twistAngle = angle between the two fundamental cycles of the torus in Radians
// aspectRatio = aspect ratio of torus, Lx / Ly
// lLLGammaX = color-entangled LLL boundary condition twisting angle along x
// lLLGammaY = color-entangled LLL boundary condition twisting angle along y
// nbrPseudopotentials = array of the number of pseudo-potentials
// pseudoPotentials = array of the pseudo-potentials
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeChernInsulatorSingleBandHamiltonian::ParticleOnLatticeChernInsulatorSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, Abstract2DTightBindingModel* tightBindingModel,
        int nbrColor, double twistAngle, double aspectRatio, double lLLGammaX, double lLLGammaY, int nbrPseudopotentials, double* pseudoPotentials,
        AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;

  this->NbrColor = nbrColor; // specified externally (ignore the FCI content in this class)
  this->TwistAngle = twistAngle;
  this->LLLGammaX = lLLGammaX;
  this->LLLGammaY = lLLGammaY;
  if (aspectRatio < 0)
      this->AspectRatio = ((double) nbrSiteX) / nbrSiteY;
  else
      this->AspectRatio = aspectRatio;
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

  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  

  // instance of this class is supposed to be used for colorful LLL calculations. NOT FCI-FQH interpolations
  this->InterpolationToFQH = 1.0;
  this->EvaluateFQHInteractionFactors();

  this->HermitianSymmetryFlag = true;
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

ParticleOnLatticeChernInsulatorSingleBandHamiltonian::~ParticleOnLatticeChernInsulatorSingleBandHamiltonian()
{
  delete[] this->Pseudopotentials;
  delete[] this->LaguerrePolynomials;
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeChernInsulatorSingleBandHamiltonian::EvaluateInteractionFactors()
{
    long TotalNbrInteractionFactors = 0;

    ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
    for (int kx = 0; kx < this->NbrSiteX; ++kx)
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
            OneBodyBasis[Index] = this->TightBindingModel->GetOneBodyMatrix(Index);
        }

    if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
        this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
        this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
            this->NbrSectorIndicesPerSum[i] = 0;
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                    {
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
                        if (Index1 < Index2)
                            ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2)];
                    }
        this->SectorIndicesPerSum = new int*[this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
        {
            if (this->NbrSectorIndicesPerSum[i] > 0)
            {
                this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];
                this->NbrSectorIndicesPerSum[i] = 0;
            }
        }
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                    {
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
                        if (Index1 < Index2)
                        {
                            int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2);
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
                int kx1, ky1, kx2, ky2;
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index1, kx1, ky1);
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index2, kx2, ky2);
                for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
                {
                    int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                    int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                    int kx3, ky3, kx4, ky4;
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index3, kx3, ky3);
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index4, kx4, ky4);

                    // the InteractionFactors are supposed to be the coefficients to  A+_3 A_1 A+_4 A_2
                    this->InteractionFactors[i][Index] = 0.0;
                    TotalNbrInteractionFactors++;
                    ++Index;
                }
            }
        }
    }
    else // Boson
    {
        this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
        this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
            this->NbrSectorIndicesPerSum[i] = 0;
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                    {
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
                        if (Index1 <= Index2)
                            ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2)];
                    }
        this->SectorIndicesPerSum = new int*[this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
        {
            if (this->NbrSectorIndicesPerSum[i] > 0)
            {
                this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];
                this->NbrSectorIndicesPerSum[i] = 0;
            }
        }
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                    {
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
                        if (Index1 <= Index2)
                        {
                            int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2);
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
                int kx1, ky1, kx2, ky2;
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index1, kx1, ky1);
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index2, kx2, ky2);
                for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
                {
                    int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                    int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                    int kx3, ky3, kx4, ky4;
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index3, kx3, ky3);
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index4, kx4, ky4);

                    // the InteractionFactors are supposed to be the coefficients to  A+_3 A_1 A+_4 A_2
                    this->InteractionFactors[i][Index] = 0.0;

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
}

// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices

void ParticleOnLatticeChernInsulatorSingleBandHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
}

// perform gauge transform to the LLL gauge
//

void ParticleOnLatticeChernInsulatorSingleBandHamiltonian::TransformToLLLGauge(ComplexMatrix& gauge, double gammaX, double gammaY)
{
    double NbrFlux = ((double)(this->NbrSiteX * this->NbrSiteY)) / this->NbrColor; // not necessarily an integer!
    for (int i = 0; i < this->NbrSectorSums; ++i)
    {
        int Index = 0;
        for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
        {
            int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
            int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
            int kx1, ky1, kx2, ky2;
            this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index1, kx1, ky1);
            this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index2, kx2, ky2);
            for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
            {
                int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                int kx3, ky3, kx4, ky4;
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index3, kx3, ky3);
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index4, kx4, ky4);

                // InteractionFactor = coefficient to  A+_3 A_1 A+_4 A_2
                this->InteractionFactors[i][Index] /= gauge[ky3][kx3] * Phase(- 2 * M_PI * ky3 * gammaX / NbrFlux);
                this->InteractionFactors[i][Index] /= gauge[ky4][kx4] * Phase(- 2 * M_PI * ky4 * gammaX / NbrFlux);
                this->InteractionFactors[i][Index] *= gauge[ky1][kx1] * Phase(- 2 * M_PI * ky1 * gammaX / NbrFlux);
                this->InteractionFactors[i][Index] *= gauge[ky2][kx2] * Phase(- 2 * M_PI * ky2 * gammaX / NbrFlux);
                ++Index;
            }
        }
    }
}

// evaluate all FQH interaction factors
//

void ParticleOnLatticeChernInsulatorSingleBandHamiltonian::EvaluateFQHInteractionFactors()
{
    long TotalNbrInteractionFactors = 0;

    int nbrSectorSums;
    int* nbrSectorIndicesPerSum;
    int** sectorIndicesPerSum;
    Complex** interactionFactors;

    if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
        nbrSectorSums = this->NbrSiteX * this->NbrSiteY;
        nbrSectorIndicesPerSum = new int[nbrSectorSums];
        for (int i = 0; i < nbrSectorSums; ++i)
            nbrSectorIndicesPerSum[i] = 0;
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                    {
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
                        if (Index1 < Index2)
                            ++nbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2)];
                    }
        sectorIndicesPerSum = new int*[nbrSectorSums];
        for (int i = 0; i < nbrSectorSums; ++i)
        {
            if (nbrSectorIndicesPerSum[i] > 0)
            {
                sectorIndicesPerSum[i] = new int[2 * nbrSectorIndicesPerSum[i]];
                nbrSectorIndicesPerSum[i] = 0;
            }
        }
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                    {
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
                        if (Index1 < Index2)
                        {
                            int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2);
                            sectorIndicesPerSum[TmpSum][nbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
                            sectorIndicesPerSum[TmpSum][1 + (nbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
                            ++nbrSectorIndicesPerSum[TmpSum];
                        }
                    }

        interactionFactors = new Complex* [nbrSectorSums];
        for (int i = 0; i < nbrSectorSums; ++i)
        {
            interactionFactors[i] = new Complex[nbrSectorIndicesPerSum[i] * nbrSectorIndicesPerSum[i]];
            int Index = 0;
            for (int j1 = 0; j1 < nbrSectorIndicesPerSum[i]; ++j1)
            {
                int Index1 = sectorIndicesPerSum[i][j1 << 1];
                int Index2 = sectorIndicesPerSum[i][(j1 << 1) + 1];
                int kx1, ky1, kx2, ky2;
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index1, kx1, ky1);
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index2, kx2, ky2);
                for (int j2 = 0; j2 < nbrSectorIndicesPerSum[i]; ++j2)
                {
                    int Index3 = sectorIndicesPerSum[i][j2 << 1];
                    int Index4 = sectorIndicesPerSum[i][(j2 << 1) + 1];
                    int kx3, ky3, kx4, ky4;
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index3, kx3, ky3);
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index4, kx4, ky4);

                    // the InteractionFactors are supposed to be the coefficients to  A+_3 A_1 A+_4 A_2
                    interactionFactors[i][Index] = 0.0;
                    interactionFactors[i][Index] += this->EvaluateFQHInteractionCoefficient(kx3, ky3, kx4, ky4, kx2, ky2, kx1, ky1);
                    interactionFactors[i][Index] -= this->EvaluateFQHInteractionCoefficient(kx4, ky4, kx3, ky3, kx2, ky2, kx1, ky1);
                    interactionFactors[i][Index] -= this->EvaluateFQHInteractionCoefficient(kx3, ky3, kx4, ky4, kx1, ky1, kx2, ky2);
                    interactionFactors[i][Index] += this->EvaluateFQHInteractionCoefficient(kx4, ky4, kx3, ky3, kx1, ky1, kx2, ky2);
                    interactionFactors[i][Index] *= -1.0;

                    TotalNbrInteractionFactors++;
                    ++Index;
                }
            }
        }
    }
    else // Boson
    {
        nbrSectorSums = this->NbrSiteX * this->NbrSiteY;
        nbrSectorIndicesPerSum = new int[nbrSectorSums];
        for (int i = 0; i < nbrSectorSums; ++i)
            nbrSectorIndicesPerSum[i] = 0;
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                    {
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
                        if (Index1 <= Index2)
                            ++nbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2)];
                    }
        sectorIndicesPerSum = new int*[nbrSectorSums];
        for (int i = 0; i < nbrSectorSums; ++i)
        {
            if (nbrSectorIndicesPerSum[i] > 0)
            {
                sectorIndicesPerSum[i] = new int[2 * nbrSectorIndicesPerSum[i]];
                nbrSectorIndicesPerSum[i] = 0;
            }
        }
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                    {
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx2, ky2);
                        if (Index1 <= Index2)
                        {
                            int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1 + kx2, ky1 + ky2);
                            sectorIndicesPerSum[TmpSum][nbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
                            sectorIndicesPerSum[TmpSum][1 + (nbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
                            ++nbrSectorIndicesPerSum[TmpSum];
                        }
                    }

        interactionFactors = new Complex* [nbrSectorSums];
        for (int i = 0; i < nbrSectorSums; ++i)
        {
            interactionFactors[i] = new Complex[nbrSectorIndicesPerSum[i] * nbrSectorIndicesPerSum[i]];
            int Index = 0;
            for (int j1 = 0; j1 < nbrSectorIndicesPerSum[i]; ++j1)
            {
                int Index1 = sectorIndicesPerSum[i][j1 << 1];
                int Index2 = sectorIndicesPerSum[i][(j1 << 1) + 1];
                int kx1, ky1, kx2, ky2;
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index1, kx1, ky1);
                this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index2, kx2, ky2);
                for (int j2 = 0; j2 < nbrSectorIndicesPerSum[i]; ++j2)
                {
                    int Index3 = sectorIndicesPerSum[i][j2 << 1];
                    int Index4 = sectorIndicesPerSum[i][(j2 << 1) + 1];
                    int kx3, ky3, kx4, ky4;
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index3, kx3, ky3);
                    this->TightBindingModel->GetLinearizedMomentumIndexSafe(Index4, kx4, ky4);

                    // the InteractionFactors are supposed to be the coefficients to  A+_3 A_1 A+_4 A_2
                    interactionFactors[i][Index] = 0.0;
                    interactionFactors[i][Index] += this->EvaluateFQHInteractionCoefficient(kx3, ky3, kx4, ky4, kx2, ky2, kx1, ky1);
                    interactionFactors[i][Index] += this->EvaluateFQHInteractionCoefficient(kx4, ky4, kx3, ky3, kx2, ky2, kx1, ky1);
                    interactionFactors[i][Index] += this->EvaluateFQHInteractionCoefficient(kx3, ky3, kx4, ky4, kx1, ky1, kx2, ky2);
                    interactionFactors[i][Index] += this->EvaluateFQHInteractionCoefficient(kx4, ky4, kx3, ky3, kx1, ky1, kx2, ky2);

                    if (Index3 == Index4)
                        interactionFactors[i][Index] *= 0.5;
                    if (Index1 == Index2)
                        interactionFactors[i][Index] *= 0.5;

                    TotalNbrInteractionFactors++;
                    ++Index;
                }
            }
        }
    }

    if (this->InterpolationToFQH == 1.0)
    {
        this->NbrSectorSums = nbrSectorSums;
        this->NbrSectorIndicesPerSum = new int[nbrSectorSums];
        this->SectorIndicesPerSum = new int*[nbrSectorSums];
        this->InteractionFactors = new Complex*[nbrSectorSums];
        for (int i = 0; i < nbrSectorSums; ++i)
        {
            this->NbrSectorIndicesPerSum[i] = nbrSectorIndicesPerSum[i];
            this->InteractionFactors[i] = NULL;
            if (nbrSectorIndicesPerSum[i] > 0)
            {
                this->SectorIndicesPerSum[i] = new int[2 * nbrSectorIndicesPerSum[i]];
                for (int j = 0; j < 2 * nbrSectorIndicesPerSum[i]; ++j)
                    this->SectorIndicesPerSum[i][j] = sectorIndicesPerSum[i][j];
                this->InteractionFactors[i] = new Complex[nbrSectorIndicesPerSum[i] * nbrSectorIndicesPerSum[i]];
                for (int Index = 0; Index < nbrSectorIndicesPerSum[i] * nbrSectorIndicesPerSum[i]; ++Index)
                    this->InteractionFactors[i][Index] = interactionFactors[i][Index];
            }
        }
    }
    else
    {
        bool Good = true;
        Good = Good && (nbrSectorSums == this->NbrSectorSums);
        for (int i = 0; Good && (i < nbrSectorSums); ++i)
        {
            Good = Good && (nbrSectorIndicesPerSum[i] == this->NbrSectorIndicesPerSum[i]);
            for (int j = 0; Good && (j < 2 * nbrSectorIndicesPerSum[i]); ++j)
                Good = Good && (sectorIndicesPerSum[i][j] == this->SectorIndicesPerSum[i][j]);
        }
        if (!Good)
        {
            cout << "Mismatch of InteractionFactors structure between FQH and FCI!" << endl;
            exit(1);
        }

        for (int i = 0; i < nbrSectorSums; ++i)
        {
            if (nbrSectorIndicesPerSum[i] > 0)
            {
                for (int Index = 0; Index < nbrSectorIndicesPerSum[i] * nbrSectorIndicesPerSum[i]; ++Index)
                {
                    this->InteractionFactors[i][Index] *= (1.0 - this->InterpolationToFQH);
                    this->InteractionFactors[i][Index] += this->InterpolationToFQH * interactionFactors[i][Index];
                }
            }
        }
    }

    for (int i = 0; i < nbrSectorSums; ++i)
    {
        if (nbrSectorIndicesPerSum[i] > 0)
        {
            delete[] sectorIndicesPerSum[i];
            delete[] interactionFactors[i];
        }
    }
    if (this->NbrSectorSums != 0)
    {
        delete[] interactionFactors;
        delete[] nbrSectorIndicesPerSum;
        delete[] sectorIndicesPerSum;
    }
    cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
    cout << "====================================" << endl;
}

// evaluate matrix element of the FQH pseudopotential Hamiltonian, namely the numerical coefficient in front of the a+_k1 a+_k2 a_k3 a_k4
// k1x = first kx index
// k1y = first ky index
// k2x = second kx index
// k2y = second ky index
// k3x = third kx index
// k3y = third ky index
// k4x = fourth kx index
// k4y = fourth ky index
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = pseudopotential coefficients
// return value = numerical coefficient
Complex ParticleOnLatticeChernInsulatorSingleBandHamiltonian::EvaluateFQHInteractionCoefficient(int k1x, int k1y, int k2x, int k2y, int k3x, int k3y, int k4x, int k4y)
{
    Complex Sum = 0.0;
    Complex Coefficient;
    Complex Term;
    double NbrFlux = ((double)(this->NbrSiteX * this->NbrSiteY)) / this->NbrColor; // not necessarily an integer!
    double Sine = sin(this->TwistAngle);
    double Cosine = cos(this->TwistAngle);
    double p1x = k1x + this->LLLGammaX;
    double p1y = k1y + this->LLLGammaY;
    double p2x = k2x + this->LLLGammaX;
    double p2y = k2y + this->LLLGammaY;
    double p3x = k3x + this->LLLGammaX;
    double p3y = k3y + this->LLLGammaY;
    double p4x = k4x + this->LLLGammaX;
    double p4y = k4y + this->LLLGammaY;

    // sum over all q = k3 - k2 mod (Nx, Ny)
    // need to properly poising the split when calculating the sum by splitting it into four quadrants
    // namely, we want to put the split to the point with smallest q²
    // this becomes significant for huge / tiny aspect ratio, and for large twisting angle
    int Range = 5;
    double MinQ2 = 1e10;
    int qxi = (k3x - k2x + this->NbrSiteX) % this->NbrSiteX;
    int qyi = (k3y - k2y + this->NbrSiteY) % this->NbrSiteY;
    int qx0 = qxi;
    int qy0 = qyi;
    for (int i = - Range + 1; i < Range; ++i)
    {
        int qx = qxi + i * this->NbrSiteX;
        for (int j = - Range + 1; j < Range; ++j)
        {
            int qy = qyi + j * this->NbrSiteY;
            double Q2 = M_PI * (qx * qx / this->AspectRatio - 2 * qx * qy * Cosine + qy * qy * this->AspectRatio) / fabs(NbrFlux * Sine);
            if (Q2 < MinQ2)
            {
                MinQ2 = Q2;
                qx0 = qx;
                qy0 = qy;
            }
        }
    }

    Coefficient = 1.0;
    int mx = 0;
    //for (int ix = 0; (ix % 2 == 0) || (SqrNorm(Sum + Coefficient) != SqrNorm(Sum)); ++ix) // not accurate enough for extremal aspect ratios
    for (int ix = 0; ((ix - 1) % 6 != 0) || ((Norm(Sum) + Norm(Coefficient)) != Norm(Sum)); ++ix)
    {
        Coefficient = 0.0;
        int qx = qx0 + mx * this->NbrSiteX;

        int qycenter = qy0 + int(Cosine * mx * this->NbrSiteX / (this->AspectRatio * this->NbrSiteY)) * this->NbrSiteY;

        Term = 1.0;
        int my = 0;
        //for (int iy = 0; (iy % 2 == 0) || (SqrNorm(Coefficient + Term) != SqrNorm(Coefficient)); ++iy)
        for (int iy = 0; ((iy - 1) % 6 != 0) || (Norm(Coefficient) + Norm(Term)) != Norm(Coefficient); ++iy)
        {
            int qy = qycenter + my * this->NbrSiteY;
            Term = this->GetVofQ(M_PI * (qx * qx / this->AspectRatio - 2 * qx * qy * Cosine + qy * qy * this->AspectRatio) / fabs(NbrFlux * Sine));
            Term *= Phase(2 * (M_PI / NbrFlux) * (qx * p1y - p4x * (p4y - p1y + qy) - qx * p2y - p3x * (p3y - p2y - qy) - qx * qy));
            Coefficient += Term;
            my += (iy + 1) * ((iy % 2 == 0)?1:(-1));
        }

        Sum += Coefficient;
        mx += (ix + 1) * ((ix % 2 == 0)?1:(-1));
    }

    // TODO: need to do the gauge-fixing here for adiabatic continuity.
    return Sum / (2 * fabs(NbrFlux));
}

// get fourier transform of interaction times the exponential supression
// Q2_half = one half of q² value
double ParticleOnLatticeChernInsulatorSingleBandHamiltonian::GetVofQ(double Q2_half)
{
    double Result = 0.0;
    double Q2 = 2.0 * Q2_half;
    for (int i = 0; i < this->NbrPseudopotentials; ++i)
        if (this->Pseudopotentials[i] != 0.0)
            Result += 2.0 * this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(Q2);
    return Result * exp(-Q2_half);
}


