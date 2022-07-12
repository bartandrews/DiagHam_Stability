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
#include "Hamiltonian/ParticleOnT2xT2GenericTwoBodyHamiltonian.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Polynomial/SpecialPolynomial.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


// default constructor
//

ParticleOnT2xT2GenericTwoBodyHamiltonian::ParticleOnT2xT2GenericTwoBodyHamiltonian()
{
  this->LaguerrePolynomials = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrFluxQuanta1 = number of flux quanta for the first torus
// nbrFluxQuanta2 = number of flux quanta for the second torus
// ratio1 = ratio between the width in the x direction and the width in the y direction for the first torus
// ratio2 = ratio between the width in the x direction and the width in the y direction for the second torus
// nbrPseudoPotentials = number of pseudo-potentials
// pseudoPotentialMomentum1= first torus pseudo-potential relative angular momenta
// pseudoPotentialMomentum2 = second torus pseudo-potential relative angular momenta
// pseudoPotentials = pseudo-potential coefficients
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnT2xT2GenericTwoBodyHamiltonian::ParticleOnT2xT2GenericTwoBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrFluxQuanta1, int nbrFluxQuanta2, 
										   double ratio1, double ratio2, 
										   int nbrPseudoPotentials, int* pseudoPotentialMomentum1, int* pseudoPotentialMomentum2, 
										   double* pseudoPotentials, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrFluxQuanta1 = nbrFluxQuanta1;
  this->NbrFluxQuanta2 = nbrFluxQuanta2;
  this->Ratio1 = ratio1;
  this->Ratio2 = ratio2;
  this->InvRatio1 = 1.0 / this->Ratio1;
  this->InvRatio2 = 1.0 / this->Ratio2;
  this->LzMax = this->NbrFluxQuanta1 * this->NbrFluxQuanta2 - 1;
  this->NbrPseudoPotentials = nbrPseudoPotentials;
  this->PseudoPotentialMomentum1 = new int [this->NbrPseudoPotentials];
  this->PseudoPotentialMomentum2 = new int [this->NbrPseudoPotentials];
  this->PseudoPotentials = new double [this->NbrPseudoPotentials];
  int TmpLargestRelativeMomentum = -1;
  for (int i = 0; i < this->NbrPseudoPotentials; ++i)
    {
      this->PseudoPotentialMomentum1[i] = pseudoPotentialMomentum1[i];
      if (this->PseudoPotentialMomentum1[i] > TmpLargestRelativeMomentum)
	{
	  TmpLargestRelativeMomentum = this->PseudoPotentialMomentum1[i];
	}
      this->PseudoPotentialMomentum2[i] = pseudoPotentialMomentum2[i];
      if (this->PseudoPotentialMomentum2[i] > TmpLargestRelativeMomentum)
	{
	  TmpLargestRelativeMomentum = this->PseudoPotentialMomentum2[i];
	}
      this->PseudoPotentials[i] = pseudoPotentials[i];
    }
  this->LaguerrePolynomials = new Polynomial[TmpLargestRelativeMomentum + 1];
  for (int i = 0; i <= TmpLargestRelativeMomentum; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

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

ParticleOnT2xT2GenericTwoBodyHamiltonian::~ParticleOnT2xT2GenericTwoBodyHamiltonian()
{
  delete[] this->PseudoPotentialMomentum1;
  delete[] this->PseudoPotentialMomentum2;
  delete[] this->PseudoPotentials;
  if (this->LaguerrePolynomials != 0)
    delete[] this->LaguerrePolynomials;
}
  

// evaluate all interaction factors
//   

void ParticleOnT2xT2GenericTwoBodyHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NbrSectorSums = this->NbrFluxQuanta1 * this->NbrFluxQuanta2;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;
      for (int lz1 = 0; lz1 < this->NbrFluxQuanta1; ++lz1)
	for (int lz2 = 0; lz2 < this->NbrFluxQuanta1; ++lz2)
	  for (int kz1 = 0; kz1 < this->NbrFluxQuanta2; ++kz1)
	    for (int kz2 = 0; kz2 < this->NbrFluxQuanta2; ++kz2)
	      {
		int Index1 = this->GetLinearizedIndex(lz1, kz1);
		int Index2 = this->GetLinearizedIndex(lz2, kz2);
		if (Index1 < Index2)
		  {
		    ++this->NbrSectorIndicesPerSum[((lz1 + lz2) % this->NbrFluxQuanta1) * this->NbrFluxQuanta2 + ((kz1 + kz2) % this->NbrFluxQuanta2)];
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
       for (int lz1 = 0; lz1 < this->NbrFluxQuanta1; ++lz1)
	for (int lz2 = 0; lz2 < this->NbrFluxQuanta1; ++lz2)
	  for (int kz1 = 0; kz1 < this->NbrFluxQuanta2; ++kz1)
	    for (int kz2 = 0; kz2 < this->NbrFluxQuanta2; ++kz2)
	      {
		int Index1 = this->GetLinearizedIndex(lz1, kz1);
		int Index2 = this->GetLinearizedIndex(lz2, kz2);
		if (Index1 < Index2)
		  {
 		    int TmpSum = ((lz1 + lz2) % this->NbrFluxQuanta1) * this->NbrFluxQuanta2 + ((kz1 + kz2) % this->NbrFluxQuanta2);
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
			  this->InteractionFactors[i][Index] += 0.5 * (this->EvaluateInteractionCoefficient(lz1, lz2, lz3, lz4, kz1, kz2, kz3, kz4, this->PseudoPotentialMomentum1[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       - this->EvaluateInteractionCoefficient(lz2, lz1, lz3, lz4, kz2, kz1, kz3, kz4, this->PseudoPotentialMomentum1[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       - this->EvaluateInteractionCoefficient(lz1, lz2, lz4, lz3, kz1, kz2, kz4, kz3, this->PseudoPotentialMomentum1[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       + this->EvaluateInteractionCoefficient(lz2, lz1, lz4, lz3, kz2, kz1, kz4, kz3, this->PseudoPotentialMomentum1[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex]));

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
      this->NbrSectorSums = this->NbrFluxQuanta1 * this->NbrFluxQuanta2;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;
      for (int lz1 = 0; lz1 < this->NbrFluxQuanta1; ++lz1)
	for (int lz2 = 0; lz2 < this->NbrFluxQuanta1; ++lz2)
	  for (int kz1 = 0; kz1 < this->NbrFluxQuanta2; ++kz1)
	    for (int kz2 = 0; kz2 < this->NbrFluxQuanta2; ++kz2)
	      {
		int Index1 = this->GetLinearizedIndex(lz1, kz1);
		int Index2 = this->GetLinearizedIndex(lz2, kz2);
		if (Index1 <= Index2)
		  {
		    ++this->NbrSectorIndicesPerSum[((lz1 + lz2) % this->NbrFluxQuanta1) * this->NbrFluxQuanta2 + ((kz1 + kz2) % this->NbrFluxQuanta2)];
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
       for (int lz1 = 0; lz1 < this->NbrFluxQuanta1; ++lz1)
	for (int lz2 = 0; lz2 < this->NbrFluxQuanta1; ++lz2)
	  for (int kz1 = 0; kz1 < this->NbrFluxQuanta2; ++kz1)
	    for (int kz2 = 0; kz2 < this->NbrFluxQuanta2; ++kz2)
	      {
		int Index1 = this->GetLinearizedIndex(lz1, kz1);
		int Index2 = this->GetLinearizedIndex(lz2, kz2);
		if (Index1 <= Index2)
		  {
 		    int TmpSum = ((lz1 + lz2) % this->NbrFluxQuanta1) * this->NbrFluxQuanta2 + ((kz1 + kz2) % this->NbrFluxQuanta2);
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
			  this->InteractionFactors[i][Index] += 0.5 * (this->EvaluateInteractionCoefficient(lz1, lz2, lz3, lz4, kz1, kz2, kz3, kz4, this->PseudoPotentialMomentum1[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       + this->EvaluateInteractionCoefficient(lz2, lz2, lz3, lz4, kz2, kz1, kz3, kz4, this->PseudoPotentialMomentum1[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       + this->EvaluateInteractionCoefficient(lz1, lz2, lz4, lz3, kz1, kz2, kz4, kz3, this->PseudoPotentialMomentum1[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       + this->EvaluateInteractionCoefficient(lz2, lz1, lz4, lz3, kz2, kz1, kz4, kz3, this->PseudoPotentialMomentum1[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex]));
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
// pseudoPotential = pseudopotential amplitude
// return value = numerical coefficient

double ParticleOnT2xT2GenericTwoBodyHamiltonian::EvaluateInteractionCoefficient(int lz1, int lz2, int lz3, int lz4, 
										int kz1, int kz2, int kz3, int kz4, 
										int pseudopotentialIndex1, int pseudopotentialIndex2, double pseudoPotential)
{
  return (pseudoPotential * this->EvaluateInteractionCoefficient(lz1, lz2, lz3, lz4, this->NbrFluxQuanta1, this->Ratio1, pseudopotentialIndex1)
	  * this->EvaluateInteractionCoefficient(kz1, kz2, kz3, kz4, this->NbrFluxQuanta2, this->Ratio2, pseudopotentialIndex2));
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// nbrFluxQuanta = number of flux quanta
// ratio = ratio between the width in the x direction and the width in the y direction
// pseudopotentialMomentum = pseudo-potential relative angular momenta
// return value = numerical coefficient

double ParticleOnT2xT2GenericTwoBodyHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrFluxQuanta, 
										double ratio, int pseudopotentialMomentum)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) nbrFluxQuanta);
  double TwoPIOnM = 2.0 * M_PI / ((double) nbrFluxQuanta);
  double InvRatio = 1.0 / ratio;
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Precision;
  double TmpInteraction;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 =ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  TmpInteraction = this->LaguerrePolynomials[pseudopotentialMomentum].PolynomialEvaluate(TwoPIOnM * Q2);
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
	  Precision = Coefficient;
	}
      else
	{
	  Precision = 1.0;
	  TmpInteraction = this->LaguerrePolynomials[pseudopotentialMomentum].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = InvRatio * N1 * N1 +ratio * N2 * N2;
	  TmpInteraction = this->LaguerrePolynomials[pseudopotentialMomentum].PolynomialEvaluate(TwoPIOnM * Q2);
	  Precision = 2.0 * exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += nbrFluxQuanta;
    }
  N2 = (double) (m1 - m4 - nbrFluxQuanta);
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 =ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  TmpInteraction = this->LaguerrePolynomials[pseudopotentialMomentum].PolynomialEvaluate(TwoPIOnM * Q2);
	  Coefficient = exp(- PIOnM * Q2) * TmpInteraction;
	  Precision = Coefficient;
	}
      else
	{
	  Precision = 1.0;
	  TmpInteraction = this->LaguerrePolynomials[pseudopotentialMomentum].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = InvRatio * N1 * N1 + ratio * N2 * N2;
	  TmpInteraction = this->LaguerrePolynomials[pseudopotentialMomentum].PolynomialEvaluate(TwoPIOnM * Q2);
	  Precision = 2.0 *  exp(- PIOnM * Q2) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= nbrFluxQuanta;
    }
  return (Sum / nbrFluxQuanta);
}

