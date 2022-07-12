////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//             class of hamiltonian with particles on the 4D manifold         //
//            Cylinder x Cylinder with 4D two body generic interaction        //
//                                                                            //
//                        last modification : 08/12/2016                      //
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
#include "Hamiltonian/ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


// default constructor
//

ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian::ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrFluxQuanta1 = number of flux quanta for the first sphere
// nbrFluxQuanta2 = number of flux quanta for the second sphere
// ratio1 = ratio between the length in the x direction and the length in the y direction of the first cylinder
// ratio2 = ratio between the length in the x direction and the length in the y direction of the second cylinder
// nbrPseudoPotentials = number of pseudo-potentials
// pseudoPotentialAngularMomentum1= pseudo-potential first sphere relative angular momenta
// pseudoPotentialAngularMomentum2 = pseudo-potential second sphere relative angular momenta
// pseudoPotentials = pseudo-potential coefficients
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian::ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, 
													   int nbrFluxQuanta1, int nbrFluxQuanta2, double ratio1, double ratio2, 
													   int nbrPseudoPotentials, int* pseudoPotentialAngularMomentum1, int* pseudoPotentialAngularMomentum2, 
													   double* pseudoPotentials, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrFluxQuanta1 = nbrFluxQuanta1;
  this->NbrFluxQuanta2 = nbrFluxQuanta2;
  this->LzMax = (this->NbrFluxQuanta1 + 1) * (this->NbrFluxQuanta2 + 1) - 1;
  this->NbrLzValues1 = this->NbrFluxQuanta1 + 1;
  this->NbrLzValues2 = this->NbrFluxQuanta2 + 1;
  this->Ratio1 = ratio1;
  this->Ratio2 = ratio2;
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

ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian::~ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian()
{
}
  

// evaluate all interaction factors
//   

void ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  double* TmpCoefficientCylinder1 = new double[2 * this->NbrFluxQuanta1 +1];
  for (int m = 0; m <= (2 * this->NbrFluxQuanta1); ++m)
    {
      TmpCoefficientCylinder1[m] = (exp(- 0.25 * 2.0 * M_PI * this->Ratio1 / ((double) (this->NbrFluxQuanta1 + 1)) * 
					(((double) ((m -  this->NbrFluxQuanta1)) * (m -  this->NbrFluxQuanta1))))
				    * pow(2.0 * M_PI, -0.25)); 
    }
  double* TmpCoefficientCylinder2 = new double[2 * this->NbrFluxQuanta2 +1];
  for (int m = 0; m <= (2 * this->NbrFluxQuanta2); ++m)
    {
      TmpCoefficientCylinder2[m] = (exp(- 0.25 * 2.0 * M_PI * this->Ratio2 / ((double) (this->NbrFluxQuanta2 + 1)) * 
					(((double) ((m -  this->NbrFluxQuanta2)) * (m -  this->NbrFluxQuanta2))))
				    * pow(2.0 * M_PI, -0.25)); 
    }
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
		  double TmpFactorAnnihilation = TmpCoefficientCylinder1[lz1 - lz2 + this->NbrFluxQuanta1] * TmpCoefficientCylinder2[kz1 - kz2 + this->NbrFluxQuanta2];
		  double TmpFactorAnnihilationPerm = TmpCoefficientCylinder1[lz2 - lz1 + this->NbrFluxQuanta1] * TmpCoefficientCylinder2[kz2 - kz1 + this->NbrFluxQuanta2];
		  for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		    {
		      int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		      int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		      this->GetLinearizedIndex(Index3, lz3, kz3);
		      this->GetLinearizedIndex(Index4, lz4, kz4);
		      this->InteractionFactors[i][Index] = 0.0;
		      this->InteractionFactors[i][Index] = 0.0;
		      for (int TmpPseudoPotentialIndex = 0; TmpPseudoPotentialIndex < this->NbrPseudoPotentials; ++TmpPseudoPotentialIndex)
			{
			  this->InteractionFactors[i][Index] += 0.5 * (this->EvaluateInteractionCoefficient(lz1, lz2, lz3, lz4, kz1, kz2, kz3, kz4, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], TmpCoefficientCylinder1, TmpCoefficientCylinder2)
								       - this->EvaluateInteractionCoefficient(lz2, lz1, lz3, lz4, kz2, kz1, kz3, kz4, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], TmpCoefficientCylinder1, TmpCoefficientCylinder2)
								       - this->EvaluateInteractionCoefficient(lz1, lz2, lz4, lz3, kz1, kz2, kz4, kz3, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], TmpCoefficientCylinder1, TmpCoefficientCylinder2)
								       + this->EvaluateInteractionCoefficient(lz2, lz1, lz4, lz3, kz2, kz1, kz4, kz3, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], TmpCoefficientCylinder1, TmpCoefficientCylinder2));
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
		  double TmpFactorAnnihilation = TmpCoefficientCylinder1[lz1 - lz2 + this->NbrFluxQuanta1] * TmpCoefficientCylinder2[kz1 - kz2 + this->NbrFluxQuanta2];
		  double TmpFactorAnnihilationPerm = TmpCoefficientCylinder1[lz2 - lz1 + this->NbrFluxQuanta1] * TmpCoefficientCylinder2[kz2 - kz1 + this->NbrFluxQuanta2];
		  for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		    {
		      int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		      int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		      this->GetLinearizedIndex(Index3, lz3, kz3);
		      this->GetLinearizedIndex(Index4, lz4, kz4);
		      this->InteractionFactors[i][Index] = 0.0;
		      for (int TmpPseudoPotentialIndex = 0; TmpPseudoPotentialIndex < this->NbrPseudoPotentials; ++TmpPseudoPotentialIndex)
			{			  
			  this->InteractionFactors[i][Index] += 0.5 * (this->EvaluateInteractionCoefficient(lz1, lz2, lz3, lz4, kz1, kz2, kz3, kz4, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], TmpCoefficientCylinder1, TmpCoefficientCylinder2)
								       + this->EvaluateInteractionCoefficient(lz2, lz2, lz3, lz4, kz2, kz1, kz3, kz4, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], TmpCoefficientCylinder1, TmpCoefficientCylinder2)
								       + this->EvaluateInteractionCoefficient(lz1, lz2, lz4, lz3, kz1, kz2, kz4, kz3, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], TmpCoefficientCylinder1, TmpCoefficientCylinder2)
								       + this->EvaluateInteractionCoefficient(lz2, lz1, lz4, lz3, kz2, kz1, kz4, kz3, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], TmpCoefficientCylinder1, TmpCoefficientCylinder2));
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
  delete[] TmpCoefficientCylinder1;
  delete[] TmpCoefficientCylinder2;
  
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}


// evaluate the numerical coefficient  in front of the a+_(lz1,kz1) a+_(lz2,kz2) a_(lz3,kz3) a_(lz4,kz4) coupling term
//
// lz1 = first lz index
// lz2 = second lz index
// lz3 = third lz index
// lz4 = fourth lz index
// kz1 = first kz index
// kz2 = second kz index
// kz3 = third kz index
// kz4 = fourth kz index
// pseudopotentialIndex1 = pseudopotential index for the interaction on the first cylinder 
// pseudopotentialIndex2 = pseudopotential index for the interaction on the second cylinder 
// exponentialCoefficient1 = array that contains the precomputed exponential factors for the first cylinder  
// exponentialCoefficient2 = array that contains the precomputed exponential factors for the second cylinder  
// return value = numerical coefficient

double ParticleOnCylinderxCylinderGenericTwoBodyHamiltonian::EvaluateInteractionCoefficient(int lz1, int lz2, int lz3, int lz4, 
											    int kz1, int kz2, int kz3, int kz4, 
											    int pseudopotentialIndex1, int pseudopotentialIndex2,
											    double* exponentialCoefficient1, double* exponentialCoefficient2)
{
  double Kappa2Factor1 = (2.0 * M_PI * this->Ratio1 / ((double) (this->NbrFluxQuanta1 + 1)));
  double Kappa2Factor2 = (2.0 * M_PI * this->Ratio2 / ((double) (this->NbrFluxQuanta2 + 1)));
  double TmpFactor = (exponentialCoefficient1[lz1 - lz2 + this->NbrFluxQuanta1] * exponentialCoefficient2[kz1 - kz2 + this->NbrFluxQuanta2] * 
		      exponentialCoefficient1[lz4 - lz3 + this->NbrFluxQuanta1] * exponentialCoefficient2[kz4 - kz3 + this->NbrFluxQuanta2]);
  if ((pseudopotentialIndex1 == 0) && (pseudopotentialIndex2 == 0))
    {
      return TmpFactor / sqrt(2.0 * M_PI);
    }
  if ((pseudopotentialIndex1 == 1) && (pseudopotentialIndex2 == 0))
    {
      return (((double) ((lz1 - lz2) * (lz4 - lz3))) * TmpFactor * Kappa2Factor1 / sqrt(this->Ratio1 * ((double) (this->NbrFluxQuanta1 + 1))));
    }
  if ((pseudopotentialIndex1 == 0) && (pseudopotentialIndex2 == 1))
    {
      return (((double) ((kz1 - kz2) * (kz4 - kz3))) * TmpFactor * Kappa2Factor2 / sqrt(this->Ratio2 * ((double) (this->NbrFluxQuanta2 + 1))));
    }
  return 0.0;
}

