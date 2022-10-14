////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//      class of checkerboard lattice model with interacting particles        //
//                       in the single band approximation                     //
//                                                                            //
//                        last modification : 08/09/2011                      //
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
// Modified by David Bauer -- Last updated 9/18/18
// modifications include:
// - fermions with NN and NNN interactions
// - bosons with NN and NNN interaction
// - e^{-r^4} interaction

#include "config.h"
#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"

#include <iostream>
#include <cmath>
#include <sys/time.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ostream;
using std::sin;
using std::cos;

// constructor
//
ParticleOnLatticeHofstadterSingleBandHamiltonian::ParticleOnLatticeHofstadterSingleBandHamiltonian()
{
    this->BandIndex = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellX = number of sites in the x direction
// nbrCellY = number of sites in the y direction
// bandIndex = index of band to consider
// uPotential = strength of the repulsive two body neareast neighbor interaction
// vPotential = strength of the repulsive two body second nearest neighbor interaction
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
//
ParticleOnLatticeHofstadterSingleBandHamiltonian::ParticleOnLatticeHofstadterSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, int bandIndex, double uPotential, double vPotential, double wPotential, TightBindingModelHofstadterSquare* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory)  // TightBindingModelHofstadterSquare
{
    this->Particles = particles;
    this->NbrParticles = nbrParticles;
    this->NbrSiteX = nbrCellX;
    this->NbrSiteY = nbrCellY;
    this->LzMax = nbrCellX * nbrCellY - 1;
    this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
    this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);

    this->HamiltonianShift = 0.0;
    this->TightBindingModel = tightBindingModel;
    this->FlatBand = flatBandFlag;
    this->UPotential = uPotential;
    this->VPotential = vPotential;
    this->WPotential = wPotential;
    this->RSquaredInteraction = false;
    this->FourFoldSmoothInteraction = false;
    this->BandIndex = bandIndex;
    this->Architecture = architecture;
    this->Memory = memory;
    this->OneBodyInteractionFactors = 0;
    this->FastMultiplicationFlag = false;
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
ParticleOnLatticeHofstadterSingleBandHamiltonian::~ParticleOnLatticeHofstadterSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//
void ParticleOnLatticeHofstadterSingleBandHamiltonian::EvaluateInteractionFactors()
{
    long TotalNbrInteractionFactors = 0;

    int NbrSublattices = TightBindingModel->GetNbrBands();

    ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];

    if (this->FlatBand == false)
    {
        this->OneBodyInteractionFactors = new double [this->TightBindingModel->GetNbrStatePerBand()];
	  }

    for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
        for (int ky = 0; ky < this->NbrSiteY; ++ky)
        {
            int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
            if (this->FlatBand == false)
                this->OneBodyInteractionFactors[Index] = 0.5 * this->TightBindingModel->GetEnergy(BandIndex, Index);  // 0.5* (reduces error in energies, especially for bosons NN)
            OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
        }
    }

    if (this->FlatBand == false)
    {
        for (int i=0; i<this->TightBindingModel->GetNbrStatePerBand(); ++i)
        {
            cout << "[" << this->OneBodyInteractionFactors[i] << "]" << endl;
        }
    }

	double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
	if (this->FlatBand == false)
		FactorU *= this->UPotential;
	double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
	double FactorW = this->WPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

	if (FactorU==0.0 && FactorV==0.0)
    {
      std::cerr << "Error: HofstadterHamiltonian created with interaction zero - set non-zero --u-potential or --v-potential"<<std::endl;
      exit(1);
    }

    /* ***** FERMIONS ***** */
    if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
        this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
        this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
        for (int i = 0; i < this->NbrSectorSums; ++i)
            this->NbrSectorIndicesPerSum[i] = 0;

        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
        {
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
            {
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                {
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                    {
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
                        if (Index1 < Index2)
                            ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2)];
                    }
                }
            }
        }

        cout << NbrSectorIndicesPerSum[0] << endl;
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
        {
            for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
            {
                for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                {
                    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
                    {
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
                        if (Index1 < Index2)
                        {
                            int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2);
                            this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
                            this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
                            ++this->NbrSectorIndicesPerSum[TmpSum];
                        }
                    }
                }
            }
        }

        this->InteractionFactors = new Complex* [this->NbrSectorSums];

        Complex Tmp;

        if(RSquaredInteraction==false)
        {
            for (int i = 0; i < this->NbrSectorSums; ++i)
            {
                this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
                int Index = 0;
                for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
                {
                    int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
                    int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
                    int kx1,ky1;
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
                    int kx2,ky2;
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);

                    for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
                    {
                        int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                        int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                        int kx3,ky3;
                        this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
                        int kx4,ky4;
                        this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

                        if (this->UPotential != 0.0)
                        {
													Tmp = 0;
													// NN interaction for fermions
													for (int s=0; s<NbrSublattices; ++s)
													{
														int xpos = s % this->TightBindingModel->GetXSitesInUC();
														int ypos = s / this->TightBindingModel->GetXSitesInUC();
														float lxfactor = (xpos+1) / (this->TightBindingModel->GetXSitesInUC());
														float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

														Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx2-kx4)*lxfactor);
														Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky2-ky4)*lyfactor);

														Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx2-kx3)*lxfactor);
														Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky2-ky3)*lyfactor);

														Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx1-kx4)*lxfactor);
														Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky1-ky4)*lyfactor);

														Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx1-kx3)*lxfactor);
														Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky1-ky3)*lyfactor);

													}
													if(Index3==Index4)
														Tmp *= 0.5;
													if(Index2==Index1)
														Tmp *=0.5;

													this->InteractionFactors[i][Index] = -2.0 * FactorU * Tmp;  // = -2.0
												}

//												if (this->UPotential != 0.0)
//												{
//													this->InteractionFactors[i][Index] = 0.0;
//													int xI, yI, dRx, dRy, sF;
//													Complex translationPhase;
//													int nbrNeighbors=4;
//													int dx[4]={1,-1,0,0};
//													int dy[4]={0,0,1,-1};
//													Tmp=0.0;
//													for (int s=0; s<NbrSublattices; ++s)
//													{
//														TightBindingModel->DecodeSublatticeIndex(s, xI, yI);
//
//														for (int n=0; n<nbrNeighbors; ++n)
//														{
//															sF = TightBindingModel->EncodeSublatticeIndex(xI + dx[n], yI + dy[n], dRx, dRy, translationPhase); // calculate final sublattice index.
//
//															Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, sF, sF, s)
//															* this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx3, ky3);
//															Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, sF, sF, s)
//															* this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx2, ky2, kx4, ky4);
//															Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, sF, sF, s)
//															* this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx3, ky3);
//															Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, sF, sF, s)
//															* this->ComputeEmbeddingForTwoBodyOperator(s, sF, kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3) * ComputeBlochPhases(dRx, dRy, kx1, ky1, kx4, ky4);
//														}
//													}
//
//													this->InteractionFactors[i][Index] += 2.0 * FactorU * Tmp;
//												}

                        if (this->VPotential != 0.0)
                        {
                            // NNN Interaction for fermions
                            Tmp =0.0;
                            int NbrSublattices = TightBindingModel->GetNbrBands();

														double OneOverSqrtTwo = 0.707107;  // 0.707107
                            for (int s=0; s<NbrSublattices; ++s)
                            {
                                int xpos = s % this->TightBindingModel->GetXSitesInUC();
                                int ypos = s / this->TightBindingModel->GetXSitesInUC();
                                float lxfactor = (xpos+1) / (this->TightBindingModel->GetXSitesInUC());
                                float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

																// cout<<Index1<<" "<<Index2<<" "<<Index3<<" "<<Index4<<" 0 "<<s<<" "<<tS(s,1,1)<<" "<<s<<" "<<tS(s,1,1)<<endl;
																// cout<<lxfactor<<" "<<lyfactor<<endl;

																Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky2-ky4)*lyfactor);
																Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky2-ky4)*lyfactor);

																Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky2-ky3)*lyfactor);
																Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky2-ky3)*lyfactor);

																Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky4)*lyfactor);
																Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky4)*lyfactor);

																Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky3)*lyfactor);
																Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky3)*lyfactor);

                            }
                            if (Index3 == Index4)
                                Tmp *= 0.5;
                            if (Index1 == Index2)
                                Tmp *= 0.5;
                            //std::cout << "NN interaction factor: " << Tmp << "\n";
                            this->InteractionFactors[i][Index] -= 2.0 * FactorV * Tmp;
                        }

                        ++TotalNbrInteractionFactors;
                        ++Index;
                    }
                }
            }
        }
        else if (RSquaredInteraction==true)
        {
        	cout << "RSquaredInteraction==true" << endl;
            double RSquaredFactorNN = 1;  // 0.367879
            double RSquaredFactorNNNDiag = 10;  // 0.0183156
            double RSquaredFactorNNNStraight = 0;  // 0.000000112535
            double RSquaredFactorNNNNDiag = 0;  // 0.0000000000138879

			// multiply by interaction factor for consistency
            if (this->FlatBand == false)
            	RSquaredFactorNN *= 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
            RSquaredFactorNNNDiag *= 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
            RSquaredFactorNNNStraight *= 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
            RSquaredFactorNNNNDiag *= 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

            for (int i = 0; i < this->NbrSectorSums; ++i)
            {
                this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
                int Index = 0;
                for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
                {
                    int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
                    int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
                    //  std::cout << Index1 << "\n";
                    int kx1 = this->getXMomentumFromIndex(Index1);
                    int ky1 = this->getYMomentumFromIndex(Index1);
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
                    int kx2= this->getXMomentumFromIndex(Index2);
                    int ky2 = this->getYMomentumFromIndex(Index2);
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);

                    for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
                    {
                        int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                        int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                        int kx3=this->getXMomentumFromIndex(Index3);
                        int ky3=this->getYMomentumFromIndex(Index3);
                        this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
                        int kx4 = this->getXMomentumFromIndex(Index4);
                        int ky4 = this->getYMomentumFromIndex(Index4);
                        this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

                        int NbrSublattices = TightBindingModel->GetNbrBands();
                        Tmp = 0;

                        // exp(-r^4) NN interaction for FERMIONS
                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            int xpos = s % this->TightBindingModel->GetXSitesInUC();
                            int ypos = s / this->TightBindingModel->GetXSitesInUC();
                            float lxfactor = (xpos+1) / (this->TightBindingModel->GetXSitesInUC());
                            float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx2-kx4)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky2-ky4)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx2-kx3)*lxfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky2-ky3)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx1-kx4)*lxfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky1-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx1-kx3)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky1-ky3)*lyfactor);

                        }
                        if(Index3==Index4)
                            Tmp *= 0.5;
                        if(Index2==Index1)
                            Tmp *=0.5;

                        this->InteractionFactors[i][Index] = -2.0 * RSquaredFactorNN * Tmp;

                        // exp(-r^4) NNN diagonal interaction for FERMIONS
                        Tmp =0.0;
                        double OneOverSqrtTwo = 0.707107;
                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            int xpos = s % this->TightBindingModel->GetXSitesInUC();
                            int ypos = s / this->TightBindingModel->GetXSitesInUC();
                            float lxfactor = (xpos+1) / (this->TightBindingModel->GetXSitesInUC());
                            float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky2-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx4)*lxfactor - OneOverSqrtTwo*this->KyFactor*(ky2-ky4)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky2-ky3)*lyfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx3)*lxfactor - OneOverSqrtTwo*this->KyFactor*(ky2-ky3)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky4)*lyfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx4)*lxfactor - OneOverSqrtTwo*this->KyFactor*(ky1-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky3)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx3)*lxfactor - OneOverSqrtTwo*this->KyFactor*(ky1-ky3)*lyfactor);
                        }

                        if (Index3 == Index4)
                            Tmp *= 0.5;
                        if (Index1 == Index2)
                            Tmp *= 0.5;
                        // cout << "NN interaction factor: " << Tmp << endl;
                        this->InteractionFactors[i][Index] += 2.0 * RSquaredFactorNNNDiag * Tmp;


                        Tmp =0.0;
                        // exp(-r^4) NNN straight-line, FERMIONS
                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            int xpos = s % this->TightBindingModel->GetXSitesInUC();
                            int ypos = s / this->TightBindingModel->GetXSitesInUC();
                            float lxfactor = (xpos +1) / (this->TightBindingModel->GetXSitesInUC());
                            float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(2.0*this->KxFactor*(kx2-kx4)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(2.0*this->KyFactor*(ky2-ky4)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(2.0*this->KxFactor*(kx2-kx3)*lxfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(2.0*this->KyFactor*(ky2-ky3)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(2.0*this->KxFactor*(kx1-kx4)*lxfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(2.0*this->KyFactor*(ky1-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(2.0*this->KxFactor*(kx1-kx3)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(2.0*this->KyFactor*(ky1-ky3)*lyfactor);

                        }

                        if (Index3 == Index4)
                            Tmp *= 0.5;
                        if (Index1 == Index2)
                            Tmp *= 0.5;
                        //std::cout << "NN interaction factor: " << Tmp << "\n";
                        this->InteractionFactors[i][Index] -= 2.0 * RSquaredFactorNNNStraight * Tmp;

                        Tmp =0.0;
                        double OneOverSqrtFive = 0.447214;

                        // exp(-r^4) NNNN interaction, FERMIONS
                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            int xpos = s % this->TightBindingModel->GetXSitesInUC();
                            int ypos = s / this->TightBindingModel->GetXSitesInUC();
                            float lxfactor = (xpos+1) / (this->TightBindingModel->GetXSitesInUC());
                            float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,2), s, tS(s,1,2))*Phase(OneOverSqrtFive*this->KxFactor*(kx2-kx4)*lxfactor + 2*OneOverSqrtFive*this->KyFactor*(ky2-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,1), s, tS(s,2,1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx2-kx4)*lxfactor + OneOverSqrtFive*this->KyFactor*(ky2-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,-1), s, tS(s,2,-1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx2-kx4)*lxfactor - OneOverSqrtFive*this->KyFactor*(ky2-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-2), s, tS(s,1,-2))*Phase(OneOverSqrtFive*this->KxFactor*(kx2-kx4)*lxfactor - 2*OneOverSqrtFive*this->KyFactor*(ky2-ky4)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,2), s, tS(s,1,2))*Phase(OneOverSqrtFive*this->KxFactor*(kx2-kx3)*lxfactor + 2*OneOverSqrtFive*this->KyFactor*(ky2-ky3)*lyfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,1), s, tS(s,2,1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx2-kx3)*lxfactor + OneOverSqrtFive*this->KyFactor*(ky2-ky3)*lyfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,-1), s, tS(s,2,-1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx2-kx3)*lxfactor - OneOverSqrtFive*this->KyFactor*(ky2-ky3)*lyfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-2), s, tS(s,1,-2))*Phase(OneOverSqrtFive*this->KxFactor*(kx2-kx3)*lxfactor - 2*OneOverSqrtFive*this->KyFactor*(ky2-ky3)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,2), s, tS(s,1,2))*Phase(OneOverSqrtFive*this->KxFactor*(kx1-kx4)*lxfactor + 2*OneOverSqrtFive*this->KyFactor*(ky1-ky4)*lyfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,1), s, tS(s,2,1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx1-kx4)*lxfactor + OneOverSqrtFive*this->KyFactor*(ky1-ky4)*lyfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,-1), s, tS(s,2,-1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx1-kx4)*lxfactor - OneOverSqrtFive*this->KyFactor*(ky1-ky4)*lyfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-2), s, tS(s,1,-2))*Phase(OneOverSqrtFive*this->KxFactor*(kx1-kx4)*lxfactor - 2*OneOverSqrtFive*this->KyFactor*(ky1-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,2), s, tS(s,1,2))*Phase(OneOverSqrtFive*this->KxFactor*(kx1-kx3)*lxfactor + 2*OneOverSqrtFive*this->KyFactor*(ky1-ky3)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,1), s, tS(s,2,1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx1-kx3)*lxfactor + 2*OneOverSqrtFive*this->KyFactor*(ky1-ky3)*lyfactor);  // possible mistake
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,-1), s, tS(s,2,-1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx1-kx3)*lxfactor - 2*OneOverSqrtFive*this->KyFactor*(ky1-ky3)*lyfactor);  // possible mistake
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-2), s, tS(s,1,-2))*Phase(OneOverSqrtFive*this->KxFactor*(kx1-kx3)*lxfactor - 2*OneOverSqrtFive*this->KyFactor*(ky1-ky3)*lyfactor);

                        }

                        if (Index3 == Index4)
                            Tmp *= 0.5;
                        if (Index1 == Index2)
                            Tmp *= 0.5;
                        //std::cout << "NN interaction factor: " << Tmp << "\n";
                        this->InteractionFactors[i][Index] -= 2.0 * RSquaredFactorNNNNDiag * Tmp;


                        ++TotalNbrInteractionFactors;
                        ++Index;
                    }
                }
            }

        }
        else if (FourFoldSmoothInteraction==true)
        {
        	cout << "FourFoldSmoothInteraction==true" << endl;
            double OneByE = 0.367879;
            for (int i = 0; i < this->NbrSectorSums; ++i)
            {
                this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
                int Index = 0;
                for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
                {
                    int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
                    int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
                    //  std::cout << Index1 << "\n";
                    int kx1 = this->getXMomentumFromIndex(Index1);
                    int ky1 = this->getYMomentumFromIndex(Index1);
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
                    int kx2= this->getXMomentumFromIndex(Index2);
                    int ky2 = this->getYMomentumFromIndex(Index2);
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);

                    for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
                    {
                        int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                        int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                        int kx3=this->getXMomentumFromIndex(Index3);
                        int ky3=this->getYMomentumFromIndex(Index3);
                        this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
                        int kx4 = this->getXMomentumFromIndex(Index4);
                        int ky4 = this->getYMomentumFromIndex(Index4);
                        this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

                        int NbrSublattices = TightBindingModel->GetNbrBands();
                        Tmp = 0;

                        // smooth C4 NN interaction for FERMIONS
                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            int xpos = s % this->TightBindingModel->GetXSitesInUC();
                            int ypos = s / this->TightBindingModel->GetXSitesInUC();
                            float lxfactor = (xpos +1) / (this->TightBindingModel->GetXSitesInUC());
                            float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx2-kx4)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky2-ky4)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx2-kx3)*lxfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky2-ky3)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx1-kx4)*lxfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky1-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx1-kx3)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky1-ky3)*lyfactor);

                        }
                        if(Index3==Index4)
                            Tmp *= 0.5;
                        if(Index2==Index1)
                            Tmp *=0.5;

                        this->InteractionFactors[i][Index] = -2.0 * Tmp;

                        // smooth C4 NNN diagonal interaction for FERMIONS
                        Tmp =0.0;
                        double OneOverSqrtTwo = 0.707107;
                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            int xpos = s % this->TightBindingModel->GetXSitesInUC();
                            int ypos = s / this->TightBindingModel->GetXSitesInUC();
                            float lxfactor = (xpos +1) / (this->TightBindingModel->GetXSitesInUC());
                            float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());


                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky2-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx4)*lxfactor - OneOverSqrtTwo*this->KyFactor*(ky2-ky4)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky2-ky3)*lyfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx3)*lxfactor - OneOverSqrtTwo*this->KyFactor*(ky2-ky3)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky4)*lyfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky3)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky3)*lyfactor);

                        }

                        if (Index3 == Index4)
                            Tmp *= 0.5;
                        if (Index1 == Index2)
                            Tmp *= 0.5;
                        //std::cout << "NN interaction factor: " << Tmp << "\n";
                        this->InteractionFactors[i][Index] -= 2.0 * 2.0 * OneByE * Tmp;


                        Tmp =0.0;
                        // smooth C4 NNN straight-line, FERMIONS
                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            int xpos = s % this->TightBindingModel->GetXSitesInUC();
                            int ypos = s / this->TightBindingModel->GetXSitesInUC();
                            float lxfactor = (xpos +1) / (this->TightBindingModel->GetXSitesInUC());
                            float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(2.0*this->KxFactor*(kx2-kx4)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(2.0*this->KyFactor*(ky2-ky4)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(2.0*this->KxFactor*(kx2-kx3)*lxfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(2.0*this->KyFactor*(ky2-ky3)*lyfactor);

                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(2.0*this->KxFactor*(kx1-kx4)*lxfactor);
                            Tmp -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(2.0*this->KyFactor*(ky1-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(2.0*this->KxFactor*(kx1-kx3)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(2.0*this->KyFactor*(ky1-ky3)*lyfactor);

                        }

                        if (Index3 == Index4)
                            Tmp *= 0.5;
                        if (Index1 == Index2)
                            Tmp *= 0.5;
                        //std::cout << "NN interaction factor: " << Tmp << "\n";
                        this->InteractionFactors[i][Index] -= 2.0 * OneByE * Tmp;

                        Tmp =0.0;
                        double OneOverSqrtFive = 0.447214;

                        ++TotalNbrInteractionFactors;
                        ++Index;
                    }
                }
            }

        }

        /* ***** END FERMIONS ***** */
    }
       
    /* ***** BOSONS ***** */
    else
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
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
                        if (Index1 <= Index2)
                            ++this->NbrSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2)];
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
                    {
                        int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
                        int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
                        if (Index1 <= Index2)
                        {
                            int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx1+kx2, ky1+ky2);
                            this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
                            this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
                            ++this->NbrSectorIndicesPerSum[TmpSum];
                        }
                    }

        double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
        if (this->FlatBand == false)
            FactorU *= this->UPotential;

        double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
        double FactorW = this->WPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

        this->InteractionFactors = new Complex* [this->NbrSectorSums];

        Complex Tmp;
        int NbrSublattices = TightBindingModel->GetNbrBands();
        if(RSquaredInteraction==false)
        {
        	cout << "RSquaredInteraction==false" << endl;
            for (int i = 0; i < this->NbrSectorSums; ++i)
            {
                //std::cout << "Momentum sector: " << i << "\n";
                //std::cout << "Number of states in this sector " << this->NbrSectorIndicesPerSum[i] << "\n";
                this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
                int Index = 0;
                for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
                {
                    int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
                    int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
                    int kx1,ky1;
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
                    int kx2,ky2;
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);

                    for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
                    {
                        int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                        int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                        int kx3,ky3;
                        this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
                        int kx4,ky4;
                        this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

                        Tmp=0.0;
                        //std::cout << "Calculating onsite interaction factors for bosons.\n";
                        //Onsite interaction factors for bosons
                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, s, s, s);// * this->ComputeTwoBodyMatrixElementOnSite(kx2, ky2, kx3, ky3);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, s, s, s);// * this->ComputeTwoBodyMatrixElementOnSite(kx2, ky2, kx4, ky4);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, s, s, s);// * this->ComputeTwoBodyMatrixElementOnSite(kx1, ky1, kx3, ky3);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, s, s, s);// * this->ComputeTwoBodyMatrixElementOnSite(kx1, ky1, kx4, ky4);
                        }

                        if (Index3 == Index4)
                            Tmp *= 0.5;
                        if (Index1 == Index2)
                            Tmp *= 0.5;
                        //std::cout << "Onsite interaction factor: " << Tmp << "\n";
                        this->InteractionFactors[i][Index] = 2.0 * FactorU * Tmp;
                        //interaction_output << Index << " " << Index1 << " " << Index2 << " " << Index3 << " " << Index4 << " " << this->InteractionFactors[i][Index] << endl;

                        if (this->VPotential != 0.0)
                        {
                            // NN Interaction for bosons
                            //std::cout << "Calculating NN interaction factors for bosons.\n";
                            Tmp =0.0;
                            int NbrSublattices = TightBindingModel->GetNbrBands();


                            for (int s=0; s<NbrSublattices; ++s)
                            {
                                int xpos = s % this->TightBindingModel->GetXSitesInUC();
                                int ypos = s / this->TightBindingModel->GetXSitesInUC();
                                float lxfactor = (xpos +1) / (this->TightBindingModel->GetXSitesInUC());
                                float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx2-kx4)*lxfactor);
                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky2-ky4)*lyfactor);

                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx2-kx3)*lxfactor);
                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky2-ky3)*lyfactor);

                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx1-kx4)*lxfactor);
                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky1-ky4)*lyfactor);

                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx1-kx3)*lxfactor);
                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky1-ky3)*lyfactor);

                            }

                            if (Index3 == Index4)
                                Tmp *= 0.5;
                            if (Index1 == Index2)
                                Tmp *= 0.5;
                            //std::cout << "NN interaction factor: " << Tmp << "\n";
                            this->InteractionFactors[i][Index] += 2.0 * FactorV * Tmp;
                        }

                        if(this->WPotential != 0.0)
                        {
                            //NNN interaction for bosons
                            //std::cout << "Calculating NNN interaction factors for bosons.\n";
                            Tmp =0.0;
                            int NbrSublattices = TightBindingModel->GetNbrBands();

                            for (int s=0; s<NbrSublattices; ++s)
                            {
                                int xpos = s % this->TightBindingModel->GetXSitesInUC();
                                int ypos = s / this->TightBindingModel->GetXSitesInUC();
                                float lxfactor = (xpos +1) / (this->TightBindingModel->GetXSitesInUC());
                                float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(1.0*this->KxFactor*(kx2-kx4)*lxfactor);
                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(1.0*this->KyFactor*(ky2-ky4)*lyfactor);

                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(1.0*this->KxFactor*(kx2-kx3)*lxfactor);
                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(1.0*this->KyFactor*(ky2-ky3)*lyfactor);

                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(1.0*this->KxFactor*(kx1-kx4)*lxfactor);
                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(1.0*this->KyFactor*(ky1-ky4)*lyfactor);

                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(1.0*this->KxFactor*(kx1-kx3)*lxfactor);
                                Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(1.0*this->KyFactor*(ky1-ky3)*lyfactor);

                            }

                            if (Index3 == Index4)
                                Tmp *= 0.5;
                            if (Index1 == Index2)
                                Tmp *= 0.5;
                            //std::cout << "NNN interaction factor: " << Tmp << "\n";
                            this->InteractionFactors[i][Index] += 2.0 * FactorW * Tmp;
                        }

                        ++TotalNbrInteractionFactors;
                        ++Index;
                    }

                }
                std::cout << "i: " << i << "\n";
                std::cout << "Index: " << Index << "\n";
            }
        }
        else
        {
        	cout << "RSquaredInteraction==true" << endl;

        	double RSquaredFactorOnsite = 1;  // 1
            double RSquaredFactorNN = 10;  // 0.367879
            double RSquaredFactorNNNDiag = 0;  // 0.0183156
            double RSquaredFactorNNNStraight = 0;  // 0.000000112535
            double RSquaredFactorNNNNDiag = 0;  // 0.0000000000138879

            // multiply by interaction factor for consistency
            if (this->FlatBand == false)
            	RSquaredFactorOnsite *= 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
            RSquaredFactorNN *= 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
            RSquaredFactorNNNDiag *= 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
            RSquaredFactorNNNStraight *= 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
            RSquaredFactorNNNNDiag *= 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

            double OneOverSqrtTwo = 0.707107;
            double OneOverSqrtFive = 0.447214;
            for (int i = 0; i < this->NbrSectorSums; ++i)
            {
                //std::cout << "Momentum sector: " << i << "\n";
                //std::cout << "Number of states in this sector " << this->NbrSectorIndicesPerSum[i] << "\n";
                this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
                int Index = 0;
                for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
                {
                    int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
                    int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
                    int kx1,ky1;
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index1,kx1, ky1);
                    int kx2,ky2;
                    this->TightBindingModel->GetLinearizedMomentumIndex(Index2,kx2, ky2);

                    for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
                    {
                        int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
                        int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
                        int kx3,ky3;
                        this->TightBindingModel->GetLinearizedMomentumIndex(Index3, kx3, ky3);
                        int kx4,ky4;
                        this->TightBindingModel->GetLinearizedMomentumIndex(Index4, kx4, ky4);

                        Tmp=0.0;
                        //std::cout << "Calculating onsite interaction factors for bosons.\n";
                        // exp(-r^4) Onsite interaction factors for bosons
                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, s, s, s);// * this->ComputeTwoBodyMatrixElementOnSite(kx2, ky2, kx3, ky3);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, s, s, s);// * this->ComputeTwoBodyMatrixElementOnSite(kx2, ky2, kx4, ky4);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, s, s, s);// * this->ComputeTwoBodyMatrixElementOnSite(kx1, ky1, kx3, ky3);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, s, s, s);// * this->ComputeTwoBodyMatrixElementOnSite(kx1, ky1, kx4, ky4);
                        }

                        if (Index3 == Index4)
                            Tmp *= 0.5;
                        if (Index1 == Index2)
                            Tmp *= 0.5;
                        //std::cout << "Onsite interaction factor: " << Tmp << "\n";
                        this->InteractionFactors[i][Index] = 2.0 * RSquaredFactorOnsite * Tmp;
                        //interaction_output << Index << " " << Index1 << " " << Index2 << " " << Index3 << " " << Index4 << " " << this->InteractionFactors[i][Index] << endl;

                        // exp(-r^4) NN Interaction for BOSONS
                        Tmp = 0.0;

                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            int xpos = s % this->TightBindingModel->GetXSitesInUC();
                            int ypos = s / this->TightBindingModel->GetXSitesInUC();
                            float lxfactor = (xpos +1) / (this->TightBindingModel->GetXSitesInUC());
                            float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx2-kx4)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky2-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx2-kx3)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky2-ky3)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx1-kx4)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky1-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,0), s, tS(s,1,0))*Phase(1.0*this->KxFactor*(kx1-kx3)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,1), s, tS(s,0,1))*Phase(1.0*this->KyFactor*(ky1-ky3)*lyfactor);

                        }

                        if (Index3 == Index4)
                            Tmp *= 0.5;
                        if (Index1 == Index2)
                            Tmp *= 0.5;
                        //std::cout << "NN interaction factor: " << Tmp << "\n";
                        this->InteractionFactors[i][Index] += 2.0 * RSquaredFactorNN * Tmp;
                        // END NN interaction for BOSONS

                        //NNN straight-line interaction for BOSONS
                        //std::cout << "Calculating NNN interaction factors for bosons.\n";
                        Tmp =0.0;

                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            int xpos = s % this->TightBindingModel->GetXSitesInUC();
                            int ypos = s / this->TightBindingModel->GetXSitesInUC();
                            float lxfactor = (xpos +1) / (this->TightBindingModel->GetXSitesInUC());
                            float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(1.0*this->KxFactor*(kx2-kx4)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(1.0*this->KyFactor*(ky2-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(1.0*this->KxFactor*(kx2-kx3)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(1.0*this->KyFactor*(ky2-ky3)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(1.0*this->KxFactor*(kx1-kx4)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(1.0*this->KyFactor*(ky1-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,0), s, tS(s,2,0))*Phase(1.0*this->KxFactor*(kx1-kx3)*lxfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,0,2), s, tS(s,0,2))*Phase(1.0*this->KyFactor*(ky1-ky3)*lyfactor);

                        }

                        if (Index3 == Index4)
                            Tmp *= 0.5;
                        if (Index1 == Index2)
                            Tmp *= 0.5;
                        //std::cout << "NNN interaction factor: " << Tmp << "\n";
                        this->InteractionFactors[i][Index] += 2.0 * RSquaredFactorNNNStraight * Tmp;

                        // END NNN straight-line for BOSONS

                        // NNN diagonal interaction for BOSONS
                        Tmp =0.0;
                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            int xpos = s % this->TightBindingModel->GetXSitesInUC();
                            int ypos = s / this->TightBindingModel->GetXSitesInUC();
                            float lxfactor = (xpos +1) / (this->TightBindingModel->GetXSitesInUC());
                            float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());


                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky2-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx4)*lxfactor - OneOverSqrtTwo*this->KyFactor*(ky2-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky2-ky3)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx2-kx3)*lxfactor - OneOverSqrtTwo*this->KyFactor*(ky2-ky3)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx4)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,1), s, tS(s,1,1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky3)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-1), s, tS(s,1,-1))*Phase(OneOverSqrtTwo*this->KxFactor*(kx1-kx3)*lxfactor + OneOverSqrtTwo*this->KyFactor*(ky1-ky3)*lyfactor);

                        }

                        if (Index3 == Index4)
                            Tmp *= 0.5;
                        if (Index1 == Index2)
                            Tmp *= 0.5;
                        //std::cout << "NN interaction factor: " << Tmp << "\n";
                        this->InteractionFactors[i][Index] += 2.0 * RSquaredFactorNNNDiag * Tmp;

                        // exp(-r^4) NNNN interaction for BOSONS
                        for (int s=0; s<NbrSublattices; ++s)
                        {
                            int xpos = s % this->TightBindingModel->GetXSitesInUC();
                            int ypos = s / this->TightBindingModel->GetXSitesInUC();
                            float lxfactor = (xpos +1) / (this->TightBindingModel->GetXSitesInUC());
                            float lyfactor = (ypos+1) / (this->TightBindingModel->GetYSitesInUC());

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,2), s, tS(s,1,2))*Phase(OneOverSqrtFive*this->KxFactor*(kx2-kx4)*lxfactor + 2*OneOverSqrtFive*this->KyFactor*(ky2-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,1), s, tS(s,2,1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx2-kx4)*lxfactor + OneOverSqrtFive*this->KyFactor*(ky2-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,-1), s, tS(s,2,-1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx2-kx4)*lxfactor - OneOverSqrtFive*this->KyFactor*(ky2-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-2), s, tS(s,1,-2))*Phase(OneOverSqrtFive*this->KxFactor*(kx2-kx4)*lxfactor - 2*OneOverSqrtFive*this->KyFactor*(ky2-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,2), s, tS(s,1,2))*Phase(OneOverSqrtFive*this->KxFactor*(kx2-kx3)*lxfactor + 2*OneOverSqrtFive*this->KyFactor*(ky2-ky3)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,1), s, tS(s,2,1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx2-kx3)*lxfactor + OneOverSqrtFive*this->KyFactor*(ky2-ky3)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,-1), s, tS(s,2,-1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx2-kx3)*lxfactor - OneOverSqrtFive*this->KyFactor*(ky2-ky3)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-2), s, tS(s,1,-2))*Phase(OneOverSqrtFive*this->KxFactor*(kx2-kx3)*lxfactor - 2*OneOverSqrtFive*this->KyFactor*(ky2-ky3)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,2), s, tS(s,1,2))*Phase(OneOverSqrtFive*this->KxFactor*(kx1-kx4)*lxfactor + 2*OneOverSqrtFive*this->KyFactor*(ky1-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,1), s, tS(s,2,1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx1-kx4)*lxfactor + OneOverSqrtFive*this->KyFactor*(ky1-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,2,-1), s, tS(s,2,-1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx1-kx4)*lxfactor - OneOverSqrtFive*this->KyFactor*(ky1-ky4)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, s, tS(s,1,-2), s, tS(s,1,-2))*Phase(OneOverSqrtFive*this->KxFactor*(kx1-kx4)*lxfactor - 2*OneOverSqrtFive*this->KyFactor*(ky1-ky4)*lyfactor);

                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,2), s, tS(s,1,2))*Phase(OneOverSqrtFive*this->KxFactor*(kx1-kx3)*lxfactor + 2*OneOverSqrtFive*this->KyFactor*(ky1-ky3)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,1), s, tS(s,2,1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx1-kx3)*lxfactor + 2*OneOverSqrtFive*this->KyFactor*(ky1-ky3)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,2,-1), s, tS(s,2,-1))*Phase(2*OneOverSqrtFive*this->KxFactor*(kx1-kx3)*lxfactor - 2*OneOverSqrtFive*this->KyFactor*(ky1-ky3)*lyfactor);
                            Tmp += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, s, tS(s,1,-2), s, tS(s,1,-2))*Phase(OneOverSqrtFive*this->KxFactor*(kx1-kx3)*lxfactor - 2*OneOverSqrtFive*this->KyFactor*(ky1-ky3)*lyfactor);

                        }

                        if (Index3 == Index4)
                            Tmp *= 0.5;
                        if (Index1 == Index2)
                            Tmp *= 0.5;
                        //std::cout << "NN interaction factor: " << Tmp << "\n";
                        this->InteractionFactors[i][Index] += 2.0 * RSquaredFactorNNNNDiag * Tmp;

                        ++TotalNbrInteractionFactors;
                        ++Index;
                    }

                }
                std::cout << "i: " << i << "\n";
                std::cout << "Index: " << Index << "\n";
            }
        }

        //End two body case
    }
    /* ***** END BOSONS ***** */

std::cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
std::cout << "====================================" << endl;

delete [] OneBodyBasis;
//interaction_output.close();

}


/*
// conventions adopted for matrix elements:
// modified with respect to Tang, Mei, Wen:
// unit cell is triangle standing on base. bottom left corner site A and bottom right corner B, tip is site C.


// compute the matrix element for the two body interaction between two sites A and B
//
// k1a = creation momentum along x for the B site
// k1b = creation momentum along y for the B site
// k2a = annihilation momentum along x for the B site
// k2b = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeTwoBodyMatrixElementAB(int k1a, int k1b, int k2a, int k2b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KxFactor * ((double) (k2a - k1a))));
  //Complex Tmp = Phase (0.5 * (this->KxFactor * ((double) (k2a - k1a))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and C
//
// k1a = creation momentum along x for the C site
// k1b = creation momentum along y for the C site
// k2a = annihilation momentum along x for the C site
// k2b = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeTwoBodyMatrixElementAC(int k1a, int k1b, int k2a, int k2b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KyFactor * ((double) (k2b - k1b))));
  //Complex Tmp = Phase (0.5 * (this->KyFactor * ((double) (k2b - k1b))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B and C
//
// k1a = creation momentum along x for the B site
// k1b = creation momentum along y for the B site
// k2a = creation momentum along x for the C site
// k2b = creation momentum along y for the C site
// k3a = annihilation momentum along x for the B site
// k3b = annihilation momentum along y for the B site
// k4a = annihilation momentum along x for the C site
// k4b = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeTwoBodyMatrixElementBC(int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b)
{
  Complex Tmp = 2.0 * cos (0.5 * ((this->KxFactor * ((double) (k3a - k1a))) + (this->KyFactor * ((double) (k4b - k2b)))));
  //Complex Tmp = Phase(0.5 * ((this->KxFactor * ((double) (k3a - k1a))) + (this->KyFactor * ((double) (k4b - k2b)))));
  return Tmp;
}

*/

// compute the matrix element for the two body interaction between two sites B and C
//
// subA = sublattice index of the first site
// subB = sublattice index of the second site
// kx1 = first creation momentum along x for the B site
// ky1 = first creation momentum along y for the B site
// kx2 = second creation momentum along x for the B site
// ky2 = second creation momentum along y for the B site
// kx3 = first annihilation momentum along x for the B site
// ky3 = first annihilation momentum along y for the B site
// kx4 = second annihilation momentum along x for the B site
// ky4 = second annihilation momentum along y for the B site
// return value = corresponding matrix element
Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeTwoBodyMatrixElementGenericAB(int subA, int subB, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
    return 1.0;
}




// compute the matrix element for on-site two body interaction involving sites on generic sublattic
//
// subl = sublattice index
// kx1 = first creation momentum along x for the B site
// ky1 = first creation momentum along y for the B site
// kx2 = second creation momentum along x for the B site
// ky2 = second creation momentum along y for the B site
// kx3 = first annihilation momentum along x for the B site
// ky3 = first annihilation momentum along y for the B site
// kx4 = second annihilation momentum along x for the B site
// ky4 = second annihilation momentum along y for the B site
//
// return value = corresponding matrix element
Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSite(int subl, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
    return 1.0;
}

/*Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeMatrixElement(ComplexMatrix* oneBodyBasis, int k1x, int k1y, int k2x, int k2y,
                                                                               int k3x, int k3y, int k4x, int k4y, int siteIndex)
{
    return 0.0;
}
*/
// compute the matrix element for on-site two body interaction involving A sites
//
// return value = corresponding matrix element

// Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteAA()
// {
//   return 1.0;
// }

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

// Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
// {
//   return Phase(0.5 * this->KxFactor * ((double) (kx4 + kx3 - kx2 -kx1)));
// }

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

// Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
// {
//   return Phase(0.5 * this->KyFactor * ((double) (ky4 + ky3 - ky2 -ky1)));
// }

int ParticleOnLatticeHofstadterSingleBandHamiltonian::tS(int index, int xtrans, int ytrans)
{
    int UCSitesX = this->TightBindingModel->GetXSitesInUC();
    int UCSitesY = this->TightBindingModel->GetYSitesInUC();

    int xpos = index % UCSitesX;
    int ypos = index / UCSitesX;

    int xpos2,ypos2;
    if(xtrans > 0)
        xpos2 = (xpos + xtrans) % UCSitesX;
    else
        xpos2 = (xpos + xtrans + UCSitesX) % UCSitesX;
        // xpos2 = UCSitesX + (xpos + xtrans) % UCSitesX;

    if(ytrans > 0)
        ypos2 = (ypos + ytrans) % UCSitesY;
    else
        ypos2 = (ypos + ytrans + UCSitesY) % UCSitesY;
		    // ypos2 = UCSitesY + (ypos + ytrans) % UCSitesY;

    int rval = xpos2 + ypos2*UCSitesX;
    // int rval = ypos2 + xpos2*UCSitesY;

    return rval;

}

// compute the matrix element for on-site two body interaction involving sites on generic sublattic
//
// s1 = sublattice index for the first creation operator
// s2 = sublattice index for the second annihilation operator
// kx1 = first creation momentum along x on the first sublattice
// ky1 = first creation momentum along y on the first sublattice
// kx2 = second creation momentum along x on the second sublattice
// ky2 = second creation momentum along y on the second sublattice
// kx3 = first annihilation momentum along x on the second sublattice
// ky3 = first annihilation momentum along y on the second sublattice
// kx4 = second annihilation momentum along x on the first sublattice
// ky4 = second annihilation momentum along y on the first sublattice
//
// return value = corresponding matrix element
 Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeEmbeddingForTwoBodyOperator(int s1, int s2, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
	// cout << s1 << " " << s2 << " " << this->KxFactor << " " << this->KyFactor << endl;
	double phase = this->TightBindingModel->GetEmbeddingPhase(s2, this->KxFactor * kx3, this->KyFactor * ky3);
  phase += this->TightBindingModel->GetEmbeddingPhase(s1, this->KxFactor * kx4, this->KyFactor * ky4);
  phase -= this->TightBindingModel->GetEmbeddingPhase(s1, this->KxFactor * kx1, this->KyFactor * ky1);
  phase -= this->TightBindingModel->GetEmbeddingPhase(s2, this->KxFactor * kx2, this->KyFactor * ky2);

  return Polar(phase);
}

// compute the matrix element for on-site two body interaction involving sites on generic sublattice
//
// dRx = number of unit vector translations along x-direction from EncodeSublatticeIndex (translations back to unit cell)
// dRy = number of unit vector translations along y-direction from EncodeSublatticeIndex (translations back to unit cell)
// kx2 = second creation momentum along x for the translated site
// ky2 = second creation momentum along y for the translated site
// kx3 = first annihilation momentum along x for the translated site
// ky3 = first annihilation momentum along y for the translated site
//
// return value = corresponding matrix element
Complex ParticleOnLatticeHofstadterSingleBandHamiltonian::ComputeBlochPhases(int dRx, int dRy, int kx2, int ky2, int kx3, int ky3)
{
  double phase = this->KxFactor * dRx * (kx2-kx3) + this->KyFactor * dRy * (ky2-ky3);
  return Polar(phase);
}
