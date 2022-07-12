////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//     class of bosons on the CP2 with generic two body interaction           //
//                                                                            // 
//                                                                            //
//                        last modification : 14/02/2013                      //
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
#include "Hamiltonian/ParticleOnCP2GenericTwoBodyHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"
#include "MathTools/SU3ClebschGordanCoefficients.h"
#include "MathTools/SU3IrreducibleRepresentations.h"

#include <iostream>
#include <cmath>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::sin;
using std::cos;



// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// pseudoPotential = pointer to an array containing the SU(3) pseudopotentials describing the interaction
// nbrPseudoPotentials = number of pseudo-potentials
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCP2GenericTwoBodyHamiltonian::ParticleOnCP2GenericTwoBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta, 
									       double* pseudoPotential, int nbrPseudoPotentials,
									       AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->NbrLzValue = (this->NbrFluxQuanta + 1)*(this->NbrFluxQuanta + 2)/2;
  this->LzMax = this->NbrLzValue - 1;
  this->PseudoPotentialIndexMax = 0;
  this->PseudoPotential = new double[this->NbrFluxQuanta + 1];
  for (int i = 0; i < nbrPseudoPotentials; ++i)
    {
      this->PseudoPotential[i] = pseudoPotential[i];
    }    
  for (int i = nbrPseudoPotentials; i <= this->NbrFluxQuanta; ++i)
    {
      this->PseudoPotential[i] = 0.0;
    }
  this->PseudoPotentialIndexMax = nbrPseudoPotentials - 1;
  this->quantumNumberTz = new int [this->NbrLzValue];
  this->quantumNumberY = new int [this->NbrLzValue];
  this->quantumNumberR = new int [this->NbrLzValue];
  this->quantumNumberS = new int [this->NbrLzValue];
  this->OneBodyTermFlag = false;
  
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
  {
    this->ParticlesFermions = (FermionOnCP2*) this->Particles;
    this->ParticlesFermions->GetQuantumNumbersFromLinearizedIndex(this->quantumNumberTz, this->quantumNumberY, this->quantumNumberR, this->quantumNumberS);
  }
  else
    {
    this->ParticlesBosons = (BosonOnCP2*) this->Particles;
    this->ParticlesBosons->GetQuantumNumbersFromLinearizedIndex(this->quantumNumberTz, this->quantumNumberY, this->quantumNumberR, this->quantumNumberS);
    }
  
//   for (int i = 0; i<NbrLzValue; i++)
//   {
//    cout << i << " " << this->quantumNumberR[i] << " " << this->quantumNumberS[i] << endl; 
//   }
  this->HamiltonianShift = 0.0;
  
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  cout << "done" << endl;
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

ParticleOnCP2GenericTwoBodyHamiltonian::~ParticleOnCP2GenericTwoBodyHamiltonian()
{
  delete[] this->PseudoPotential;
  delete[] this->quantumNumberR;
  delete[] this->quantumNumberS;
  delete[] this->quantumNumberTz;
  delete[] this->quantumNumberY;
}

// evaluate all interaction factors
//   

void ParticleOnCP2GenericTwoBodyHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  
  SU3ClebschGordanCoefficients Clebsch(this->NbrFluxQuanta, 0, this->NbrFluxQuanta, 0);
  SU3IrreducibleRepresentations Representation(this->NbrFluxQuanta, 0);	  
  
  this->NbrSectorSums = (2*this->NbrFluxQuanta + 1)*(2*this->NbrFluxQuanta + 2)/2;
  this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
  for (int i = 0; i < this->NbrSectorSums; ++i)
  this->NbrSectorIndicesPerSum[i] = 0;  
      
  
  
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
  {
    for (int r1 = 0; r1 <= this->NbrFluxQuanta; ++r1)
    for (int r2 = 0; r2 <= this->NbrFluxQuanta; ++r2)
      for (int s1 = 0; s1 <= this->NbrFluxQuanta - r1; ++s1)
	for (int s2 = 0; s2 <= this->NbrFluxQuanta - r2; ++s2) 
	   {
	      int tz1 = r1 - s1;
	      int y1 = 3*(r1+s1) - 2*this->NbrFluxQuanta;
	      int tz2 = r2 - s2;
	      int y2 = 3*(r2 + s2) - 2*this->NbrFluxQuanta;
	      int Index1 = this->ParticlesFermions->GetLinearizedIndex(tz1, y1, 1);
	      int Index2 = this->ParticlesFermions->GetLinearizedIndex(tz2, y2, 1);
	      if (Index1 < Index2)
		++this->NbrSectorIndicesPerSum[this->ParticlesFermions->GetLinearizedIndex(tz1 + tz2, y1 + y2 , 2)];    
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
  for (int r1 = 0; r1 <= this->NbrFluxQuanta; ++r1)
    for (int r2 = 0; r2 <= this->NbrFluxQuanta; ++r2)
      for (int s1 = 0; s1 <= this->NbrFluxQuanta - r1; ++s1)
	for (int s2 = 0; s2 <= this->NbrFluxQuanta - r2; ++s2) 
	  {
	      int tz1 = r1 - s1;
	      int y1 = 3*(r1+s1) - 2*this->NbrFluxQuanta;
	      int tz2 = r2 - s2;
	      int y2 = 3*(r2 + s2) - 2*this->NbrFluxQuanta;
	      int Index1 = this->ParticlesFermions->GetLinearizedIndex(tz1, y1, 1);
	      int Index2 = this->ParticlesFermions->GetLinearizedIndex(tz2, y2, 1);
	      if (Index1 < Index2)
		{
		   int TmpSum = this->ParticlesFermions->GetLinearizedIndex(tz1 + tz2, y1 + y2 , 2);
		   this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		++this->NbrSectorIndicesPerSum[TmpSum];    
	      }
	  }
      
      
      this->InteractionFactors = new double* [this->NbrSectorSums];
      
      //cout << this->NbrSectorSums << endl;
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new double[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int l1 = 0; l1 < this->NbrSectorIndicesPerSum[i]; ++l1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][l1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(l1 << 1) + 1];
	      int tz1 = this->quantumNumberTz[Index1];
	      int y1 = this->quantumNumberY[Index1];
	      int tz2 = this->quantumNumberTz[Index2];
	      int y2 = this->quantumNumberY[Index2];
	      int qIndice1 = Representation.GetQIndices(tz1, y1, 0);
	      int qIndice2 = Representation.GetQIndices(tz2, y2, 0);
	      for (int l2 = 0; l2 < this->NbrSectorIndicesPerSum[i]; ++l2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][l2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(l2 << 1) + 1];
		  int tz3 = this->quantumNumberTz[Index3];
		  int y3 = this->quantumNumberY[Index3];
		  int tz4 = this->quantumNumberTz[Index4];
		  int y4 = this->quantumNumberY[Index4];
		  int qIndice3 = Representation.GetQIndices(tz3, y3, 0);
		  int qIndice4 = Representation.GetQIndices(tz4, y4, 0);
		  
// 		  cout << tz1 << " " << y1 << " " << qIndice1 << " ; " << tz2 << " " << y2 << " " << qIndice2 << " ; " << tz3 << " " << y3 << " " << qIndice3 << " ; " << tz4 << " " << y4 << " " << qIndice4 << endl;
		 
		  double TmpInteractionFactor;
		  double TmpFactor = 0;
		  
		  for (int j = 0; j < Clebsch.GetNbrPQRepresentations(); ++j)
		  {
		    TmpInteractionFactor = 0;
		    if ((this->PseudoPotential[j] != 0) && (Clebsch.GetDegeneracy(j, qIndice1, qIndice2) != 0))
		      {
// 			cout << j << " " << Clebsch.GetDegeneracy(j, qIndice1, qIndice2) << " " << Clebsch.GetRepresentationDimension(j) << endl;
		      
		      for (int qIndex = 0; qIndex < Clebsch.GetRepresentationDimension(j); ++qIndex)
			{
// 			  cout << qIndex << " " << Clebsch.GetClebschGordanCoefficient(j, qIndice1, qIndice2, qIndex) << endl;
			 TmpInteractionFactor -= Clebsch.GetClebschGordanCoefficient(j, qIndice1, qIndice2, qIndex)*Clebsch.GetClebschGordanCoefficient(j, qIndice3, qIndice4, qIndex);
			 TmpInteractionFactor += Clebsch.GetClebschGordanCoefficient(j, qIndice2, qIndice1, qIndex)*Clebsch.GetClebschGordanCoefficient(j, qIndice3, qIndice4, qIndex);
			 TmpInteractionFactor += Clebsch.GetClebschGordanCoefficient(j, qIndice1, qIndice2, qIndex)*Clebsch.GetClebschGordanCoefficient(j, qIndice4, qIndice3, qIndex);
			 TmpInteractionFactor -= Clebsch.GetClebschGordanCoefficient(j, qIndice2, qIndice1, qIndex)*Clebsch.GetClebschGordanCoefficient(j, qIndice4, qIndice3, qIndex);
			}
		      TmpInteractionFactor *= PseudoPotential[j];
// 		      cout << TmpInteractionFactor << endl;
		      }
		    
		    TmpFactor += TmpInteractionFactor;
// 		    cout << j << " " << i << " " << Index << " = " << TmpInteractionFactor << endl;
		  }
	      this->InteractionFactors[i][Index] =  0.5*TmpFactor;
	      
	      TotalNbrInteractionFactors++;
	      ++Index;
	    }
	   }
	}
  }
  else
  {
    for (int r1 = 0; r1 <= this->NbrFluxQuanta; ++r1)
      for (int r2 = 0; r2 <= this->NbrFluxQuanta; ++r2)
	for (int s1 = 0; s1 <= this->NbrFluxQuanta - r1; ++s1)
	  for (int s2 = 0; s2 <= this->NbrFluxQuanta - r2; ++s2) 
	    {
	      int tz1 = r1 - s1;
	      int y1 = 3*(r1+s1) - 2*this->NbrFluxQuanta;
	      int tz2 = r2 - s2;
	      int y2 = 3*(r2 + s2) - 2*this->NbrFluxQuanta;
	      int Index1 = this->ParticlesBosons->GetLinearizedIndex(tz1, y1, 1);
	      int Index2 = this->ParticlesBosons->GetLinearizedIndex(tz2, y2, 1);
	      if (Index1 <= Index2)
		++this->NbrSectorIndicesPerSum[this->ParticlesBosons->GetLinearizedIndex(tz1 + tz2, y1 + y2 , 2)];    
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
    for (int r1 = 0; r1 <= this->NbrFluxQuanta; ++r1)
      for (int r2 = 0; r2 <= this->NbrFluxQuanta; ++r2)
	for (int s1 = 0; s1 <= this->NbrFluxQuanta - r1; ++s1)
	  for (int s2 = 0; s2 <= this->NbrFluxQuanta - r2; ++s2) 
	    {
	      int tz1 = r1 - s1;
	      int y1 = 3*(r1+s1) - 2*this->NbrFluxQuanta;
	      int tz2 = r2 - s2;
	      int y2 = 3*(r2 + s2) - 2*this->NbrFluxQuanta;
	      int Index1 = this->ParticlesBosons->GetLinearizedIndex(tz1, y1, 1);
	      int Index2 = this->ParticlesBosons->GetLinearizedIndex(tz2, y2, 1);
	      if (Index1 <= Index2)
		{
		   int TmpSum = this->ParticlesBosons->GetLinearizedIndex(tz1 + tz2, y1 + y2 , 2);
		   this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		++this->NbrSectorIndicesPerSum[TmpSum];    
		}
	    }
      
    this->InteractionFactors = new double* [this->NbrSectorSums];
      
      //cout << this->NbrSectorSums << endl;
    for (int i = 0; i < this->NbrSectorSums; ++i)
      {
	 this->InteractionFactors[i] = new double[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	 int Index = 0;
	 for (int l1 = 0; l1 < this->NbrSectorIndicesPerSum[i]; ++l1)
	   {
	     int Index1 = this->SectorIndicesPerSum[i][l1 << 1];
	     int Index2 = this->SectorIndicesPerSum[i][(l1 << 1) + 1];
	     int tz1 = this->quantumNumberTz[Index1];
	     int y1 = this->quantumNumberY[Index1];
	     int tz2 = this->quantumNumberTz[Index2];
	     int y2 = this->quantumNumberY[Index2];
	     int qIndice1 = Representation.GetQIndices(tz1, y1, 0);
	     int qIndice2 = Representation.GetQIndices(tz2, y2, 0);
	     for (int l2 = 0; l2 < this->NbrSectorIndicesPerSum[i]; ++l2)
	      {
		 int Index3 = this->SectorIndicesPerSum[i][l2 << 1];
		 int Index4 = this->SectorIndicesPerSum[i][(l2 << 1) + 1];
		 int tz3 = this->quantumNumberTz[Index3];
		 int y3 = this->quantumNumberY[Index3];
		 int tz4 = this->quantumNumberTz[Index4];
		 int y4 = this->quantumNumberY[Index4];
		 int qIndice3 = Representation.GetQIndices(tz3, y3, 0);
		 int qIndice4 = Representation.GetQIndices(tz4, y4, 0);
		  
// 		  cout << tz1 << " " << y1 << " " << qIndice1 << " ; " << tz2 << " " << y2 << " " << qIndice2 << " ; " << tz3 << " " << y3 << " " << qIndice3 << " ; " << tz4 << " " << y4 << " " << qIndice4 << endl;
		 
		  double TmpInteractionFactor;
		  double TmpFactor = 0;
		  
		  for (int j = 0; j < Clebsch.GetNbrPQRepresentations(); ++j)
		  {
		    TmpInteractionFactor = 0;
		    if ((this->PseudoPotential[j] != 0) && (Clebsch.GetDegeneracy(j, qIndice1, qIndice2) != 0))
		      {
// 			cout << j << " " << Clebsch.GetDegeneracy(j, qIndice1, qIndice2) << " " << Clebsch.GetRepresentationDimension(j) << endl;
		      
		      for (int qIndex = 0; qIndex < Clebsch.GetRepresentationDimension(j); ++qIndex)
			{
// 			  cout << qIndex << " " << Clebsch.GetClebschGordanCoefficient(j, qIndice1, qIndice2, qIndex) << endl;
			 TmpInteractionFactor += Clebsch.GetClebschGordanCoefficient(j, qIndice1, qIndice2, qIndex)*Clebsch.GetClebschGordanCoefficient(j, qIndice3, qIndice4, qIndex);
			 TmpInteractionFactor += Clebsch.GetClebschGordanCoefficient(j, qIndice2, qIndice1, qIndex)*Clebsch.GetClebschGordanCoefficient(j, qIndice3, qIndice4, qIndex);
			 TmpInteractionFactor += Clebsch.GetClebschGordanCoefficient(j, qIndice1, qIndice2, qIndex)*Clebsch.GetClebschGordanCoefficient(j, qIndice4, qIndice3, qIndex);
			 TmpInteractionFactor += Clebsch.GetClebschGordanCoefficient(j, qIndice2, qIndice1, qIndex)*Clebsch.GetClebschGordanCoefficient(j, qIndice4, qIndice3, qIndex);
			}
		      TmpInteractionFactor *= PseudoPotential[j];
// 		      cout << TmpInteractionFactor << endl;
		     }
		    
		    TmpFactor += TmpInteractionFactor;
// 		    cout << j << " " << i << " " << Index << " = " << TmpInteractionFactor << endl;
		  }
	      this->InteractionFactors[i][Index] =  TmpFactor;
	      if (Index3 == Index4)
		this->InteractionFactors[i][Index] *= 0.5;
	      if (Index1 == Index2)
		this->InteractionFactors[i][Index] *= 0.5;
	      this->InteractionFactors[i][Index] *= 0.5;

	      TotalNbrInteractionFactors++;
	      ++Index;
	   }
	}
      }
    }
  if (this->OneBodyTermFlag == true)
    {
      this->NbrOneBodyInteractionFactors = 0;
      for (int i = 0; i <= this->LzMax; ++i)
	if (this->OneBodyPotentials[i] != 0)
	  ++this->NbrOneBodyInteractionFactors;
      if (this->NbrOneBodyInteractionFactors != 0)
	{
	  this->OneBodyMValues = new int[this->NbrOneBodyInteractionFactors];
	  this->OneBodyNValues = new int[this->NbrOneBodyInteractionFactors];
	  this->OneBodyInteractionFactors = new double[this->NbrLzValue];
	  this->NbrOneBodyInteractionFactors = 0;
	  int TmpNbrOrbitals = 0;
	  for (int i = 0; i <= this->LzMax; ++i)
	  {
	    this->OneBodyInteractionFactors[TmpNbrOrbitals] = 0.5*this->OneBodyPotentials[i];
	    if (this->OneBodyPotentials[i] != 0)
	      {
		this->OneBodyMValues[this->NbrOneBodyInteractionFactors] = i;
		this->OneBodyNValues[this->NbrOneBodyInteractionFactors] = i;
		++this->NbrOneBodyInteractionFactors;
	      }	 
	    ++ TmpNbrOrbitals;
	  }
	}
      else
	{
	  delete[] this->OneBodyPotentials;
	  this->OneBodyTermFlag = false;
	}
    } 
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}
