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
// Modified by DB, 8/28/2014 -- Implemented fermions with NN interaction

#include "config.h"
#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian.h"
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
ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian::ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian()
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

ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian::ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, int bandIndex, double uPotential, double vPotential, double wPotential, double threeBodyPotential, TightBindingModelHofstadterSquare* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
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
  this->ThreeBodyPotential = threeBodyPotential;
  this->BandIndex = bandIndex;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  this->NBodyValue = 3;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
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

ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian::~ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian::EvaluateInteractionFactors()
{
    std::cout << "Entering EvaluateInterationFactors()." << endl;
 // std::ofstream interaction_output;
   //interaction_output.open("SingleBandInteraction.txt");
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactors = new double [this->TightBindingModel->GetNbrStatePerBand()];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	if (this->FlatBand == false)
	  this->OneBodyInteractionFactors[Index] = 0.5 * this->TightBindingModel->GetEnergy(BandIndex, Index);
	OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
      }

  if (this->FlatBand == false)
    for (int i=0; i<this->TightBindingModel->GetNbrStatePerBand(); ++i)
      {
	cout << "[" << this->OneBodyInteractionFactors[i] << "]" << endl;
      }
      //Fermions
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
      {    
          // Do nothing.
      // End Fermionic case
      }
      else // Bosonic Statistics
      {
        // Two body
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
         // std::cout << "i: " << i << "\n";
        // std::cout << "Index: " << Index << "\n";
        }
        //Three Body
        
        this->NbrNBodySectorSums = this->NbrSiteX * this->NbrSiteY;
        this->NbrNBodySectorIndicesPerSum = new int[this->NbrNBodySectorSums];
        
        for (int i = 0; i < this->NbrNBodySectorSums; ++i)
          this->NbrNBodySectorIndicesPerSum[i] = 0;      
       
        for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
          for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
            for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
              for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
                for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
                  for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
                  {
                    int Index1 = (kx1 * this->NbrSiteY) + ky1;
                    int Index2 = (kx2 * this->NbrSiteY) + ky2;
                    int Index3 = (kx3 * this->NbrSiteY) + ky3;
                    if ((Index1 <= Index2) && (Index2 <= Index3))
                     ++this->NbrNBodySectorIndicesPerSum[(((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY)];    
                  }
      
      this->NBodySectorIndicesPerSum = new int* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
      {
        if (this->NbrNBodySectorIndicesPerSum[i]  > 0)
        {
          this->NBodySectorIndicesPerSum[i] = new int[this->NBodyValue * this->NbrNBodySectorIndicesPerSum[i]];      
          this->NbrNBodySectorIndicesPerSum[i] = 0;
        }
      }
      
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
        for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
          for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
            for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
              for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
                for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
                {
                  int Index1 = (kx1 * this->NbrSiteY) + ky1;
                  int Index2 = (kx2 * this->NbrSiteY) + ky2;
                  int Index3 = (kx3 * this->NbrSiteY) + ky3;
                  if ((Index1 <= Index2) && (Index2 <= Index3))
                  {
                    int TmpSum = (((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY);
                    this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * 3] = Index1;
                    this->NBodySectorIndicesPerSum[TmpSum][1 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = Index2;
                    this->NBodySectorIndicesPerSum[TmpSum][2 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = Index3;
                    ++this->NbrNBodySectorIndicesPerSum[TmpSum];    
                  }
                }

      
      int KxIn[3];
      int KyIn[3];
      int IndexIn[3];

      int** Permutations = 0; 
      double* PermutationSign = 0; 
      int NbrPermutations = this->ComputePermutations(Permutations, PermutationSign);

      int TmpLargestSector = 0;
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
        if (this->NbrNBodySectorIndicesPerSum[i] > TmpLargestSector)
          TmpLargestSector = this->NbrNBodySectorIndicesPerSum[i];

      //const int NbrSublattices2 = 1;
      //int NbrSublattices = TightBindingModel->GetNbrBands();
      Complex** TmpIn = new Complex*[TmpLargestSector];
      Complex** TmpOut = new Complex*[TmpLargestSector];
      for(int i = 0; i < TmpLargestSector;++i)
      {
        TmpIn[i] = new Complex[NbrSublattices];
        TmpOut[i] = new Complex[NbrSublattices];
      }
      //Complex* TmpBBBIn = new Complex[TmpLargestSector];
      //Complex* TmpBBBOut = new Complex[TmpLargestSector];

      // double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      double ThreeBodyFactor = this->ThreeBodyPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
      {
        this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
        int Index = 0;
        for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
        {
          IndexIn[0] = this->NBodySectorIndicesPerSum[i][j1 * 3];
          IndexIn[1] = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 1];
          IndexIn[2] = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 2];
          KxIn[0] = IndexIn[0] / this->NbrSiteY;
          KyIn[0] = IndexIn[0] % this->NbrSiteY;
          KxIn[1] = IndexIn[1] / this->NbrSiteY;
          KyIn[1] = IndexIn[1] % this->NbrSiteY;
          KxIn[2] = IndexIn[2] / this->NbrSiteY;
          KyIn[2] = IndexIn[2] % this->NbrSiteY;

        Complex TmpIn2 = 0.0;
        Complex TmpOut2 = 0.0;
        //Complex TmpBBBIn2 = 0.0;
        //Complex TmpBBBOut2 = 0.0;
            for(int s = 0; s < NbrSublattices; ++s)
            {
               TmpIn[j1][s] =  0.0;
              TmpOut[j1][s] = 0.0;
            }

          for (int l1 = 0; l1 < NbrPermutations; ++l1)
          {
            int* TmpPerm = Permutations[l1];
   

            for(int s = 0; s < NbrSublattices; ++s)
            {
              Complex SumTmpIn2 = Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][s] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][s] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][s]) ;
              //Complex SumTmpBBBIn2 = Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1]) ;
      
              Complex SumTmpOut2 = OneBodyBasis[IndexIn[TmpPerm[0]]][0][s] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][s] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][s] ;
              //Complex SumTmpBBBOut2 = OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1] ;

          
              if ((IndexIn[TmpPerm[0]] == IndexIn[TmpPerm[1]]) && (IndexIn[TmpPerm[1]] == IndexIn[TmpPerm[2]]))
              {
                SumTmpIn2 *= 1.0/6.0;
                //SumTmpIn2 *= 1.0/6.0;
                SumTmpOut2 *= 1.0/6.0;
                //SumTmpBBBOut2 *= 1.0/6.0;
          
              }
              else
              {
                if( (IndexIn[TmpPerm[0]] == IndexIn[TmpPerm[1]])  || (IndexIn[TmpPerm[1]] == IndexIn[TmpPerm[2]])  || (IndexIn[TmpPerm[0]] == IndexIn[TmpPerm[2]]) )
                {
                  SumTmpIn2 *= 0.5;
                  //SumTmpBBBIn2 *= 0.5;
                  SumTmpOut2 *= 0.5;
                  //SumTmpBBBOut2 *= 0.5;
                }
              }
      
              TmpIn[j1][s] += SumTmpIn2 ;
              //TmpBBBIn2 += SumTmpBBBIn2 ;
      
              TmpOut[j1][s] += SumTmpOut2 ;
              //TmpBBBOut2 += SumTmpBBBOut2 ;


            }

          //TmpIn[j1][s] =  TmpIn2;
          //TmpOut[j1][s] = TmpOut2;
          //TmpBBBIn[j1] =  TmpBBBIn2;
          //TmpBBBOut[j1] = TmpBBBOut2;

      } }
   // } } /*
      
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
      {
        for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
        {
          Complex TmpNBodyInteractionFactor =0.0;
          for(int s=0;s<NbrSublattices;++s)
          {
            TmpNBodyInteractionFactor += TmpIn[j1][s] * TmpOut[j2][s];
          }
          this->NBodyInteractionFactors[i][Index] = 2.0 * ThreeBodyFactor * (TmpNBodyInteractionFactor);                 
          TotalNbrInteractionFactors++;
          ++Index;
        }
      }

  
}
     // for (int i=0; i<NbrPermutations; ++i)
     //   delete [] Permutations[i];
     // delete [] Permutations;
      //delete [] PermutationSign;


      //delete[] TmpIn;
      //delete[] TmpOut;
      //delete[] TmpBBBIn;
      //delete[] TmpBBBOut;
      //*/
       
      //End Bosonic statistics
      }
  
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
  
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
Complex ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementGenericAB(int subA, int subB, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
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
Complex ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementOnSite(int subl, int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
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

int ParticleOnLatticeHofstadterSingleBandThreeBodyHamiltonian::tS(int index, int xtrans, int ytrans)
{
  int UCSitesX = this->TightBindingModel->GetXSitesInUC();
  int UCSitesY = this->TightBindingModel->GetYSitesInUC();

  int xpos = index % UCSitesX;
  int ypos = index / UCSitesX;

  int xpos2,ypos2;
  if(xtrans >0)
     xpos2 = (xpos + xtrans) % UCSitesX;
  else
     xpos2 = (xpos + xtrans + UCSitesX) % UCSitesX;

  if(ytrans >0)
    ypos2 = (ypos + ytrans) % UCSitesY;
  else
    ypos2 = (ypos + ytrans + UCSitesY) % UCSitesY;

  int rval = xpos2 + ypos2*UCSitesX;
  return rval;

}

