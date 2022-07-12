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
//                   T2 x S2 with 4D two body generic interaction             //
//                                                                            //
//                        last modification : 09/03/2017                      //
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
#include "Hamiltonian/ParticleOnT2xS2GenericTwoBodyHamiltonian.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Polynomial/SpecialPolynomial.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


// default constructor
//

ParticleOnT2xS2GenericTwoBodyHamiltonian::ParticleOnT2xS2GenericTwoBodyHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrFluxQuantumTorus = number of flux quanta piercing the torus
// ratio = ratio between the length in the x direction and the length in the y direction of the torus
// nbrFluxQuantumSphere = number of flux quanta piercing the sphere
// nbrPseudoPotentials = number of pseudo-potentials
// pseudoPotentialMomentumTorus = torus pseudo-potential relative momenta
// pseudoPotentialMomentumSphere = sphere pseudo-potential relative angular momenta
// pseudoPotentials = pseudo-potential coefficients
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnT2xS2GenericTwoBodyHamiltonian::ParticleOnT2xS2GenericTwoBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, 
													   int nbrFluxQuantumTorus, double ratio, int nbrFluxQuantumSphere, 
													   int nbrPseudoPotentials, int* pseudoPotentialMomentumTorus, int* pseudoPotentialMomentumSphere, 
													   double* pseudoPotentials, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrFluxQuanta1 = nbrFluxQuantumTorus;
  this->NbrFluxQuanta2 = nbrFluxQuantumSphere;
  this->LzMax = (this->NbrFluxQuanta1 * (this->NbrFluxQuanta2 + 1)) - 1;
  this->NbrLzValues1 = this->NbrFluxQuanta1;
  this->NbrLzValues2 = this->NbrFluxQuanta2 + 1;
  this->Ratio = ratio;
  this->NbrPseudoPotentials = nbrPseudoPotentials;
  this->PseudoPotentialAngularMomentum1 = new int [this->NbrPseudoPotentials];
  this->PseudoPotentialAngularMomentum2 = new int [this->NbrPseudoPotentials];
  this->PseudoPotentials = new double [this->NbrPseudoPotentials];
  int TmpLargestRelativeMomentum = -1;
  for (int i = 0; i < this->NbrPseudoPotentials; ++i)
    {
      this->PseudoPotentialAngularMomentum1[i] = pseudoPotentialMomentumTorus[i];
      if (this->PseudoPotentialAngularMomentum1[i] > TmpLargestRelativeMomentum)
	{
	  TmpLargestRelativeMomentum = this->PseudoPotentialAngularMomentum1[i];
	}
      this->PseudoPotentialAngularMomentum2[i] = pseudoPotentialMomentumSphere[i];
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

ParticleOnT2xS2GenericTwoBodyHamiltonian::~ParticleOnT2xS2GenericTwoBodyHamiltonian()
{
}
  

// evaluate all interaction factors
//   

void ParticleOnT2xS2GenericTwoBodyHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ClebschGordanCoefficients Clebsch (this->NbrFluxQuanta2, this->NbrFluxQuanta2);
  this->GetIndices();
  this->InteractionFactors = new double* [this->NbrSectorSums];
  int lz1;
  int lz2;
  int lz3;
  int lz4;
  int kz1;
  int kz2;
  int kz3;
  int kz4;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
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
		      this->InteractionFactors[i][Index] = 0.0;
		      for (int TmpPseudoPotentialIndex = 0; TmpPseudoPotentialIndex < this->NbrPseudoPotentials; ++TmpPseudoPotentialIndex)
			{
			  this->InteractionFactors[i][Index] += 0.5 * (this->EvaluateInteractionCoefficient(lz1, lz2, lz3, lz4, kz1, kz2, kz3, kz4, Clebsch, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       - this->EvaluateInteractionCoefficient(lz2, lz1, lz3, lz4, kz2, kz1, kz3, kz4, Clebsch, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       - this->EvaluateInteractionCoefficient(lz1, lz2, lz4, lz3, kz1, kz2, kz4, kz3, Clebsch, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       + this->EvaluateInteractionCoefficient(lz2, lz1, lz4, lz3, kz2, kz1, kz4, kz3, Clebsch, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex]));
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
		      this->InteractionFactors[i][Index] += 0.5 * (this->EvaluateInteractionCoefficient(lz1, lz2, lz3, lz4, kz1, kz2, kz3, kz4, Clebsch, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								   + this->EvaluateInteractionCoefficient(lz2, lz2, lz3, lz4, kz2, kz1, kz3, kz4, Clebsch, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													  this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								   + this->EvaluateInteractionCoefficient(lz1, lz2, lz4, lz3, kz1, kz2, kz4, kz3, Clebsch, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													  this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								   + this->EvaluateInteractionCoefficient(lz2, lz1, lz4, lz3, kz2, kz1, kz4, kz3, Clebsch, this->PseudoPotentialAngularMomentum1[TmpPseudoPotentialIndex], 
													  this->PseudoPotentialAngularMomentum2[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex]));
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
// ky1 = first ky index
// ky2 = second ky index
// ky3 = third ky index
// ky4 = fourth ky index
// lz1 = first lz index
// lz2 = second lz index
// lz3 = third lz index
// lz4 = fourth lz index
// clebsch = reference on the Clebsch-Gordan coefficients for the sphere geometry
// pseudopotentialIndex1 = pseudopotential index for the interaction on the torus
// pseudopotentialIndex2 = pseudopotential index for the interaction on the sphere 
// pseudoPotential = pseudopotential amplitude
// return value = numerical coefficient

double ParticleOnT2xS2GenericTwoBodyHamiltonian::EvaluateInteractionCoefficient(int ky1, int ky2, int ky3, int ky4, 
										int lz1, int lz2, int lz3, int lz4, 
										ClebschGordanCoefficients& clebsch,
										int pseudopotentialIndex1, int pseudopotentialIndex2, double pseudoPotential)
{
  return (pseudoPotential * this->EvaluateTorusInteractionCoefficient(ky1, ky2, ky3, ky4, this->NbrFluxQuanta1, this->Ratio, pseudopotentialIndex1)
	  * this->EvaluateSphereInteractionCoefficient(lz1, lz2, lz3, lz4, clebsch, this->NbrFluxQuanta2, pseudopotentialIndex2));
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term for the torus
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// nbrFluxQuanta = number of flux quanta
// ratio = ratio between the width in the x direction and the width in the y direction
// pseudopotentialMomentum = pseudo-potential relative momentum
// return value = numerical coefficient

double ParticleOnT2xS2GenericTwoBodyHamiltonian::EvaluateTorusInteractionCoefficient(int m1, int m2, int m3, int m4, int nbrFluxQuanta, 
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

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term for the sphere part
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// clebsch = reference on the Clebsch-Gordan coefficients for the sphere geometry
// nbrFluxQuanta = number of flux quanta
// pseudopotentialMomentum = pseudo-potential relative momentum
// return value = numerical coefficient

double ParticleOnT2xS2GenericTwoBodyHamiltonian::EvaluateSphereInteractionCoefficient(int m1, int m2, int m3, int m4, ClebschGordanCoefficients& clebsch,
										      int nbrFluxQuanta,int pseudopotentialMomentum)
{
  int AngularMaxMomentum = 2 * nbrFluxQuanta;
  if ((AngularMaxMomentum >= (2 * pseudopotentialMomentum)) && 
      (abs(2 * (m1 + m2 - nbrFluxQuanta)) <= (AngularMaxMomentum - 2 * pseudopotentialMomentum)))
    {
      return (clebsch.GetCoefficient((2 * m1) - nbrFluxQuanta, (2 * m2) - nbrFluxQuanta, AngularMaxMomentum - 2 * pseudopotentialMomentum) 
	      * clebsch.GetCoefficient((2 * m3) - nbrFluxQuanta, (2 * m4) - nbrFluxQuanta, AngularMaxMomentum - 2 * pseudopotentialMomentum));
    }
  else
    {
      return 0;
    }
}

// get all the indices that should appear in the annihilation/creation operators
//

void ParticleOnT2xS2GenericTwoBodyHamiltonian::GetIndices()
{
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      int TmpNbrLzSumPerSphere = 2 * this->NbrFluxQuanta2 + 1;
      this->NbrSectorSums = TmpNbrLzSumPerSphere * this->NbrFluxQuanta1;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;      
      for (int ky1 = 0; ky1 < this->NbrFluxQuanta1; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrFluxQuanta1; ++ky2)
	  for (int lz1 = 0; lz1 <= this->NbrFluxQuanta2; ++lz1)
	    for (int lz2 = 0; lz2 <= this->NbrFluxQuanta2; ++lz2)
	      {
		int Index1 = this->GetLinearizedIndex(ky1, lz1);
		int Index2 = this->GetLinearizedIndex(ky2, lz2);
		if (Index1 < Index2)
		  {
		    ++this->NbrSectorIndicesPerSum[(((ky1 + ky2) % this->NbrFluxQuanta1) * TmpNbrLzSumPerSphere) + (lz1 + lz2)];
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
      for (int ky1 = 0; ky1 < this->NbrFluxQuanta1; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrFluxQuanta1; ++ky2)
	  for (int lz1 = 0; lz1 <= this->NbrFluxQuanta2; ++lz1)
	    for (int lz2 = 0; lz2 <= this->NbrFluxQuanta2; ++lz2)
	      {
		int Index1 = this->GetLinearizedIndex(ky1, lz1);
		int Index2 = this->GetLinearizedIndex(ky2, lz2);
		if (Index1 < Index2)
		  {
 		    int TmpSum = (((ky1 + ky2) % this->NbrFluxQuanta1) * TmpNbrLzSumPerSphere) + (lz1 + lz2);
 		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
 		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
 		    ++this->NbrSectorIndicesPerSum[TmpSum];
		  }
	      }
    }
  else
    {
      int TmpNbrLzSumPerSphere = 2 * this->NbrFluxQuanta2 + 1;
      this->NbrSectorSums = TmpNbrLzSumPerSphere * this->NbrFluxQuanta1;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;      
      for (int ky1 = 0; ky1 < this->NbrFluxQuanta1; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrFluxQuanta1; ++ky2)
	  for (int lz1 = 0; lz1 <= this->NbrFluxQuanta2; ++lz1)
	    for (int lz2 = 0; lz2 <= this->NbrFluxQuanta2; ++lz2)
	      {
		int Index1 = this->GetLinearizedIndex(ky1, lz1);
		int Index2 = this->GetLinearizedIndex(ky2, lz2);
		if (Index1 <= Index2)
		  {
		    ++this->NbrSectorIndicesPerSum[(((ky1 + ky2) % this->NbrFluxQuanta1) * TmpNbrLzSumPerSphere) + (lz1 + lz2)];
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
      for (int ky1 = 0; ky1 < this->NbrFluxQuanta1; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrFluxQuanta1; ++ky2)
	  for (int lz1 = 0; lz1 <= this->NbrFluxQuanta2; ++lz1)
	    for (int lz2 = 0; lz2 <= this->NbrFluxQuanta2; ++lz2)
	      {
		int Index1 = this->GetLinearizedIndex(ky1, lz1);
		int Index2 = this->GetLinearizedIndex(ky2, lz2);
		if (Index1 <= Index2)
		  {
 		    int TmpSum = ((ky1 + ky2) % this->NbrFluxQuanta1) * (TmpNbrLzSumPerSphere) + (lz1 + lz2);
 		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
 		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
 		    ++this->NbrSectorIndicesPerSum[TmpSum];
		  }
	      }
    }
}
