////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//     class of hamiltonian with particles on the 4D manifold S2 x S2         //
//                     with 4D generic two-body interaction                   //
//                                                                            //
//                        last modification : 07/12/2016                      //
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
#include "Hamiltonian/ParticleOnS2xS2GenericTwoBodyHamiltonian.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


// default constructor
//

ParticleOnS2xS2GenericTwoBodyHamiltonian::ParticleOnS2xS2GenericTwoBodyHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrFluxQuanta1 = number of flux quanta for the first sphere
// nbrFluxQuanta2 = number of flux quanta for the second sphere
// nbrPseudoPotentials = number of pseudo-potentials
// pseudoPotentialAngularMomentum1= pseudo-potential first sphere relative angular momenta
// pseudoPotentialAngularMomentum2 = pseudo-potential second sphere relative angular momenta
// pseudoPotentials = pseudo-potential coefficients
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnS2xS2GenericTwoBodyHamiltonian::ParticleOnS2xS2GenericTwoBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta1, int nbrFluxQuanta2, 
										   int nbrPseudoPotentials, int* pseudoPotentialAngularMomentum1, int* pseudoPotentialAngularMomentum2, 
										   double* pseudoPotentials, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrFluxQuanta1 = nbrFluxQuanta1;
  this->NbrFluxQuanta2 = nbrFluxQuanta2;
  this->LzMax = (this->NbrFluxQuanta1 + 1) * (this->NbrFluxQuanta2 + 1) - 1;
  this->NbrLzValues1 = this->NbrFluxQuanta1 + 1;
  this->NbrLzValues2 = this->NbrFluxQuanta2 + 1;
  this->NbrPseudoPotentials = nbrPseudoPotentials;
  this->PseudoPotentialAngularMomentum1 = new int [this->NbrPseudoPotentials];
  this->PseudoPotentialAngularMomentum2 = new int [this->NbrPseudoPotentials];
  this->PseudoPotentials = new double [this->NbrPseudoPotentials];
  for (int i = 0; i < this->NbrPseudoPotentials; ++i)
    {
      this->PseudoPotentialAngularMomentum1[i] = pseudoPotentialAngularMomentum1[i];
      this->PseudoPotentialAngularMomentum2[i] = pseudoPotentialAngularMomentum2[i];
      this->PseudoPotentials[i] = pseudoPotentials[i];
    }

  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyTermFlag = false;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  

  this->EvaluateInteractionFactors();

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

ParticleOnS2xS2GenericTwoBodyHamiltonian::~ParticleOnS2xS2GenericTwoBodyHamiltonian()
{
  delete[] this->PseudoPotentialAngularMomentum1;
  delete[] this->PseudoPotentialAngularMomentum2;
  delete[] this->PseudoPotentials;
}
  

// evaluate all interaction factors
//   

void ParticleOnS2xS2GenericTwoBodyHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ClebschGordanCoefficients Clebsch1 (this->NbrFluxQuanta1, this->NbrFluxQuanta1);
  ClebschGordanCoefficients Clebsch2 (this->NbrFluxQuanta2, this->NbrFluxQuanta2);
  int AngularMaxMomentum1 = 2 * this->NbrFluxQuanta1;
  int AngularMaxMomentum2 = 2 * this->NbrFluxQuanta2;
  
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      int TmpNbrLzSumPerSphere = 2 * this->NbrFluxQuanta2 + 1;
      this->NbrSectorSums = TmpNbrLzSumPerSphere * (2 * this->NbrFluxQuanta1 + 1);
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;
      for (int lz1 = 0; lz1 <= this->NbrFluxQuanta1; ++lz1)
	for (int lz2 = 0; lz2 <= this->NbrFluxQuanta1; ++lz2)
	  for (int kz1 = 0; kz1 <= this->NbrFluxQuanta2; ++kz1)
	    for (int kz2 = 0; kz2 <= this->NbrFluxQuanta2; ++kz2)
	      {
		int Index1 = this->GetLinearizedIndex(lz1, kz1);
		int Index2 = this->GetLinearizedIndex(lz2, kz2);
		if (Index1 < Index2)
		  {
		    ++this->NbrSectorIndicesPerSum[(lz1 + lz2) * (TmpNbrLzSumPerSphere) + (kz1 + kz2)];
		  }
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
       for (int lz1 = 0; lz1 <= this->NbrFluxQuanta1; ++lz1)
	for (int lz2 = 0; lz2 <= this->NbrFluxQuanta1; ++lz2)
	  for (int kz1 = 0; kz1 <= this->NbrFluxQuanta2; ++kz1)
	    for (int kz2 = 0; kz2 <= this->NbrFluxQuanta2; ++kz2)
	      {
		int Index1 = this->GetLinearizedIndex(lz1, kz1);
		int Index2 = this->GetLinearizedIndex(lz2, kz2);
		if (Index1 < Index2)
		  {
 		    int TmpSum = (lz1 + lz2) * (TmpNbrLzSumPerSphere) + (kz1 + kz2);
 		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
 		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
 		    ++this->NbrSectorIndicesPerSum[TmpSum];
 		  }
 	      }
      
       int lz1;
       int lz2;
       int lz3;
       int lz4;
       int kz1;
       int kz2;
       int kz3;
       int kz4;
       this->InteractionFactors = new double* [this->NbrSectorSums];
       for (int i = 0; i < this->NbrSectorSums; ++i)
        {
	  if (this->NbrSectorIndicesPerSum[i] > 0)
	    {
	      this->InteractionFactors[i] = new double[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
		  this->GetLinearizedIndex(Index1, lz1, kz1);
		  this->GetLinearizedIndex(Index2, lz2, kz2);
		  for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		    {
		      int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		      int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		      this->GetLinearizedIndex(Index3, lz3, kz3);
		      this->GetLinearizedIndex(Index4, lz4, kz4);
		      this->InteractionFactors[i][Index] = 0.0;
		      for (int TmpPseudoPotentialIndex = 0; TmpPseudoPotentialIndex < this->NbrPseudoPotentials; ++TmpPseudoPotentialIndex)
			{
			  if ((AngularMaxMomentum1 >= (2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex])) && (AngularMaxMomentum2 >= (2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex]))
			      && (abs(2 * (lz1 + lz2 - this->NbrFluxQuanta1)) <= (AngularMaxMomentum1 - 2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex]))
			      && (abs(2 * (kz1 + kz2 - this->NbrFluxQuanta2)) <= (AngularMaxMomentum2 - 2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex])))
			    {
			      double TmpFactorAnnihilation = (Clebsch1.GetCoefficient((2 * lz1) - this->NbrFluxQuanta1, (2 * lz2) - this->NbrFluxQuanta1, AngularMaxMomentum1 - 2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex]) 
							      * Clebsch2.GetCoefficient((2 * kz1) - this->NbrFluxQuanta2, (2 * kz2) - this->NbrFluxQuanta2, AngularMaxMomentum2 - 2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex]));
			      double TmpFactorAnnihilationPerm = (Clebsch1.GetCoefficient((2 * lz2) - this->NbrFluxQuanta1, (2 * lz1) - this->NbrFluxQuanta1, AngularMaxMomentum1 - 2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex]) 
								  * Clebsch2.GetCoefficient((2 * kz2) - this->NbrFluxQuanta2, (2 * kz1) - this->NbrFluxQuanta2, AngularMaxMomentum2 - 2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex]));
			      double TmpFactorCreation = (Clebsch1.GetCoefficient((2 * lz3) - this->NbrFluxQuanta1, (2 * lz4) - this->NbrFluxQuanta1, AngularMaxMomentum1 - 2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex]) 
							  * Clebsch2.GetCoefficient((2 * kz3) - this->NbrFluxQuanta2, (2 * kz4) - this->NbrFluxQuanta2, AngularMaxMomentum2 - 2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex]));
			      double TmpFactorCreationPerm = (Clebsch1.GetCoefficient((2 * lz4) - this->NbrFluxQuanta1, (2 * lz3) - this->NbrFluxQuanta1, AngularMaxMomentum1 - 2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex]) 
							      * Clebsch2.GetCoefficient((2 * kz4) - this->NbrFluxQuanta2, (2 * kz3) - this->NbrFluxQuanta2, AngularMaxMomentum2 - 2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex]));
			      
			      this->InteractionFactors[i][Index] -= 0.5 * this->PseudoPotentials[TmpPseudoPotentialIndex] * ((TmpFactorAnnihilation * TmpFactorCreation) - (TmpFactorAnnihilationPerm * TmpFactorCreation)
															     - (TmpFactorAnnihilation * TmpFactorCreationPerm) + (TmpFactorAnnihilationPerm * TmpFactorCreationPerm));
			    }
			}
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
		}
	    }
	}
    } 
  else 
    {
      int TmpNbrLzSumPerSphere = 2 * this->NbrFluxQuanta2 + 1;
      this->NbrSectorSums = TmpNbrLzSumPerSphere * (2 * this->NbrFluxQuanta1 + 1);
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;
      for (int lz1 = 0; lz1 <= this->NbrFluxQuanta1; ++lz1)
	for (int lz2 = 0; lz2 <= this->NbrFluxQuanta1; ++lz2)
	  for (int kz1 = 0; kz1 <= this->NbrFluxQuanta2; ++kz1)
	    for (int kz2 = 0; kz2 <= this->NbrFluxQuanta2; ++kz2)
	      {
		int Index1 = this->GetLinearizedIndex(lz1, kz1);
		int Index2 = this->GetLinearizedIndex(lz2, kz2);
		if (Index1 <= Index2)
		  {
		    ++this->NbrSectorIndicesPerSum[(lz1 + lz2) * (TmpNbrLzSumPerSphere) + (kz1 + kz2)];
		  }
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
       for (int lz1 = 0; lz1 <= this->NbrFluxQuanta1; ++lz1)
	for (int lz2 = 0; lz2 <= this->NbrFluxQuanta1; ++lz2)
	  for (int kz1 = 0; kz1 <= this->NbrFluxQuanta2; ++kz1)
	    for (int kz2 = 0; kz2 <= this->NbrFluxQuanta2; ++kz2)
	      {
		int Index1 = this->GetLinearizedIndex(lz1, kz1);
		int Index2 = this->GetLinearizedIndex(lz2, kz2);
		if (Index1 <= Index2)
		  {
 		    int TmpSum = (lz1 + lz2) * (TmpNbrLzSumPerSphere) + (kz1 + kz2);
 		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
 		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
 		    ++this->NbrSectorIndicesPerSum[TmpSum];
 		  }
 	      }
      
       int lz1;
       int lz2;
       int lz3;
       int lz4;
       int kz1;
       int kz2;
       int kz3;
       int kz4;
       this->InteractionFactors = new double* [this->NbrSectorSums];
       for (int i = 0; i < this->NbrSectorSums; ++i)
        {
	  if (this->NbrSectorIndicesPerSum[i] > 0)
	    {
	      this->InteractionFactors[i] = new double[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
		{
		  int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
		  int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
		  this->GetLinearizedIndex(Index1, lz1, kz1);
		  this->GetLinearizedIndex(Index2, lz2, kz2);
		  for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		    {
		      int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		      int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		      this->GetLinearizedIndex(Index3, lz3, kz3);
		      this->GetLinearizedIndex(Index4, lz4, kz4);
		      this->InteractionFactors[i][Index] = 0.0;
		      for (int TmpPseudoPotentialIndex = 0; TmpPseudoPotentialIndex < this->NbrPseudoPotentials; ++TmpPseudoPotentialIndex)
			{
			  if ((AngularMaxMomentum1 >= (2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex])) && (AngularMaxMomentum2 >= (2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex]))
			      && (abs(2 * (lz1 + lz2 - this->NbrFluxQuanta1)) <= (AngularMaxMomentum1 - 2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex]))
			      && (abs(2 * (kz1 + kz2 - this->NbrFluxQuanta2)) <= (AngularMaxMomentum2 - 2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex])))
			    {
			      double TmpFactorAnnihilation = (Clebsch1.GetCoefficient((2 * lz1) - this->NbrFluxQuanta1, (2 * lz2) - this->NbrFluxQuanta1, AngularMaxMomentum1 - 2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex]) 
							      * Clebsch2.GetCoefficient((2 * kz1) - this->NbrFluxQuanta2, (2 * kz2) - this->NbrFluxQuanta2, AngularMaxMomentum2 - 2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex]));
			      double TmpFactorAnnihilationPerm = (Clebsch1.GetCoefficient((2 * lz2) - this->NbrFluxQuanta1, (2 * lz1) - this->NbrFluxQuanta1, AngularMaxMomentum1 - 2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex]) 
								  * Clebsch2.GetCoefficient((2 * kz2) - this->NbrFluxQuanta2, (2 * kz1) - this->NbrFluxQuanta2, AngularMaxMomentum2 - 2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex]));
			      double TmpFactorCreation = (Clebsch1.GetCoefficient((2 * lz3) - this->NbrFluxQuanta1, (2 * lz4) - this->NbrFluxQuanta1, AngularMaxMomentum1 - 2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex]) 
							  * Clebsch2.GetCoefficient((2 * kz3) - this->NbrFluxQuanta2, (2 * kz4) - this->NbrFluxQuanta2, AngularMaxMomentum2 - 2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex]));
			      double TmpFactorCreationPerm = (Clebsch1.GetCoefficient((2 * lz4) - this->NbrFluxQuanta1, (2 * lz3) - this->NbrFluxQuanta1, AngularMaxMomentum1 - 2 * this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex]) 
							      * Clebsch2.GetCoefficient((2 * kz4) - this->NbrFluxQuanta2, (2 * kz3) - this->NbrFluxQuanta2, AngularMaxMomentum2 - 2 * this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex]));
			      
			      this->InteractionFactors[i][Index] += 0.5 * this->PseudoPotentials[TmpPseudoPotentialIndex] * ((TmpFactorAnnihilation * TmpFactorCreation) + (TmpFactorAnnihilationPerm * TmpFactorCreation)
															     + (TmpFactorAnnihilation * TmpFactorCreationPerm) + (TmpFactorAnnihilationPerm * TmpFactorCreationPerm));
			    }
			}
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
    }
  
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}


