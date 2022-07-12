////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of hamiltonian associated to particles on a             //
//        4D space  torus x cylinder with a generic two-body interaction      //
//                          and magnetic translations                         //
//                                                                            //
//                        last modification : 07/03/2017                      //
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


#include "Hamiltonian/ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "GeneralTools/StringTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Polynomial/SpecialPolynomial.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian::ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrFluxQuantumTorus = number of flux quanta piercing the torus
// kxMomentum = momentum in the x direction for the torus
// maxMomentum = maximum Lz value reached by a particle in the state
// xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// nbrFluxQuantumCylinder = number of flux quanta piercing the cylinder
// cylinderRatio = ratio between the length in the x direction and the length in the y direction of the cylinder
// nbrPseudoPotentials = number of pseudo-potentials
// pseudoPotentialMomentumTorus = torus pseudo-potential relative momenta
// pseudoPotentialMomentumCylinder = cylinder pseudo-potential relative angular momenta
// pseudoPotentials = pseudo-potential coefficients
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian::ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, 
																	       int nbrParticles, int nbrFluxQuantumTorus, int xMomentum,
																	       double ratio, int nbrFluxQuantumCylinder,  double cylinderRatio, int nbrPseudoPotentials, 
																	       int* pseudoPotentialMomentumTorus, int* pseudoPotentialMomentumCylinder, 
																	       double* pseudoPotentials, AbstractArchitecture* architecture, long memory, 
																	       char* precalculationFileName)
{
  this->Particles = particles;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrFluxQuantumTorus = nbrFluxQuantumTorus;
  this->NbrFluxQuantumSphere = nbrFluxQuantumCylinder;
  this->NbrLzValues =  this->NbrFluxQuantumSphere + 1;
  this->LzMax = (this->NbrLzValues) * this->NbrFluxQuantumTorus - 1;  
  this->XMomentum = xMomentum;
  this->NbrParticles = nbrParticles;
  this->MaxMomentum = (this->NbrFluxQuantumSphere + 1) * this->NbrFluxQuantumTorus;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->NbrFluxQuantumTorus);
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->OneBodyInteractionFactors = 0;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->CylinderRatio = cylinderRatio;

  this->NbrPseudoPotentials = nbrPseudoPotentials;
  this->PseudoPotentialMomentumTorus = new int [this->NbrPseudoPotentials];
  this->PseudoPotentialMomentumSphere = new int [this->NbrPseudoPotentials];
  this->PseudoPotentials = new double [this->NbrPseudoPotentials];
  int TmpLargestRelativeMomentum = -1;
  for (int i = 0; i < this->NbrPseudoPotentials; ++i)
    {
      this->PseudoPotentialMomentumTorus[i] = pseudoPotentialMomentumTorus[i];
      if (this->PseudoPotentialMomentumTorus[i] > TmpLargestRelativeMomentum)
	{
	  TmpLargestRelativeMomentum = this->PseudoPotentialMomentumTorus[i];
	}
      this->PseudoPotentialMomentumSphere[i] = pseudoPotentialMomentumCylinder[i];
      if (this->PseudoPotentialMomentumSphere[i] > TmpLargestRelativeMomentum)
	{
	  TmpLargestRelativeMomentum = this->PseudoPotentialMomentumSphere[i];
	}
      this->PseudoPotentials[i] = pseudoPotentials[i];
    }
  this->LaguerrePolynomials = new Polynomial[TmpLargestRelativeMomentum + 1];
  for (int i = 0; i <= TmpLargestRelativeMomentum; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  this->EvaluateExponentialFactors();
  this->HamiltonianShift = 0.0;
  this->EvaluateInteractionFactors();
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout << "fast memory = ";
	  PrintMemorySize(cout,TmpMemory)<<endl;
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
      else
	{
	  if (this->Architecture->HasAutoLoadBalancing())
	    {
	      this->FastMultiplicationMemory(0l);
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian::~ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian() 
{
}

// evaluate all interaction factors
//   

void ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian::EvaluateInteractionFactors()
{
  double* TmpCoefficientCylinder = new double[2 * this->NbrFluxQuantumSphere +1];
  for (int m = 0; m <= (2 * this->NbrFluxQuantumSphere); ++m)
    {
      TmpCoefficientCylinder[m] = (exp(- 0.25 * 2.0 * M_PI * this->CylinderRatio / ((double) (this->NbrFluxQuantumSphere + 1)) * 
					(((double) ((m -  this->NbrFluxQuantumSphere)) * (m -  this->NbrFluxQuantumSphere))))
				    * pow(2.0 * M_PI, -0.25)); 
    }
  long TotalNbrInteractionFactors = 0;
  long TotalNbrNonZeroInteractionFactors = 0;
  double MaxCoefficient = 0.0;
  this->GetIndices();
  this->InteractionFactors = new Complex* [this->NbrSectorSums];
  int lz1;
  int lz2;
  int lz3;
  int lz4;
  int kz1;
  int kz2;
  int kz3;
  int kz4;
  if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
    {
       for (int i = 0; i < this->NbrSectorSums; ++i)
        {
	  if (this->NbrSectorIndicesPerSum[i] > 0)
	    {
	      this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
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
			  this->InteractionFactors[i][Index] += 0.5 * (this->EvaluateInteractionCoefficient(lz1, lz2, lz3, lz4, kz1, kz2, kz3, kz4, TmpCoefficientCylinder, this->PseudoPotentialMomentumTorus[TmpPseudoPotentialIndex], 
													    this->PseudoPotentialMomentumSphere[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       - this->EvaluateInteractionCoefficient(lz2, lz1, lz3, lz4, kz2, kz1, kz3, kz4, TmpCoefficientCylinder, this->PseudoPotentialMomentumTorus[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialMomentumSphere[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       - this->EvaluateInteractionCoefficient(lz1, lz2, lz4, lz3, kz1, kz2, kz4, kz3, TmpCoefficientCylinder, this->PseudoPotentialMomentumTorus[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialMomentumSphere[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								       + this->EvaluateInteractionCoefficient(lz2, lz1, lz4, lz3, kz2, kz1, kz4, kz3, TmpCoefficientCylinder, this->PseudoPotentialMomentumTorus[TmpPseudoPotentialIndex], 
													      this->PseudoPotentialMomentumSphere[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex]));

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
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
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
		      this->InteractionFactors[i][Index] += 0.5 * (this->EvaluateInteractionCoefficient(lz1, lz2, lz3, lz4, kz1, kz2, kz3, kz4, TmpCoefficientCylinder, this->PseudoPotentialMomentumTorus[TmpPseudoPotentialIndex], 
													this->PseudoPotentialMomentumSphere[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								   + this->EvaluateInteractionCoefficient(lz2, lz2, lz3, lz4, kz2, kz1, kz3, kz4, TmpCoefficientCylinder, this->PseudoPotentialMomentumTorus[TmpPseudoPotentialIndex], 
													  this->PseudoPotentialMomentumSphere[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								   + this->EvaluateInteractionCoefficient(lz1, lz2, lz4, lz3, kz1, kz2, kz4, kz3, TmpCoefficientCylinder, this->PseudoPotentialMomentumTorus[TmpPseudoPotentialIndex], 
													  this->PseudoPotentialMomentumSphere[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex])
								   + this->EvaluateInteractionCoefficient(lz2, lz1, lz4, lz3, kz2, kz1, kz4, kz3, TmpCoefficientCylinder, this->PseudoPotentialMomentumTorus[TmpPseudoPotentialIndex], 
													  this->PseudoPotentialMomentumSphere[TmpPseudoPotentialIndex], this->PseudoPotentials[TmpPseudoPotentialIndex]));
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
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "nbr non-zero interaction = " << TotalNbrNonZeroInteractionFactors << endl;
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
// exponentialCoefficient = array that contains the precomputed exponential factors for the cylinder  
// pseudopotentialIndex1 = pseudopotential index for the interaction on the torus
// pseudopotentialIndex2 = pseudopotential index for the interaction on the sphere 
// pseudoPotential = pseudopotential amplitude
// return value = numerical coefficient

double ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian::EvaluateInteractionCoefficient(int ky1, int ky2, int ky3, int ky4, 
													      int lz1, int lz2, int lz3, int lz4, 
													      double* exponentialCoefficient,
													      int pseudopotentialIndex1, int pseudopotentialIndex2, double pseudoPotential)
{
  return (pseudoPotential * this->EvaluateTorusInteractionCoefficient(ky1, ky2, ky3, ky4, this->NbrFluxQuantumTorus, this->Ratio, pseudopotentialIndex1)
	  * this->EvaluateCylinderInteractionCoefficient(lz1, lz2, lz3, lz4, pseudopotentialIndex2, exponentialCoefficient));
}

// evaluate the numerical coefficient  in front of the a+_lz1 a+_lz2 a_lz3 a_lz4 coupling term for the cylinder part
//
// lz1 = first lz index
// lz2 = second lz index
// lz3 = third lz index
// lz4 = fourth lz index
// pseudopotentialIndex = pseudopotential index for the interaction on the cylinder 
// exponentialCoefficient = array that contains the precomputed exponential factors for the cylinder  
// return value = numerical coefficient

double ParticleOnT2xCylinderWithMagneticTranslationsGenericTwoBodyHamiltonian::EvaluateCylinderInteractionCoefficient(int lz1, int lz2, int lz3, int lz4, int pseudopotentialIndex, double* exponentialCoefficient)
{
  double Kappa2Factor = (2.0 * M_PI * this->CylinderRatio / ((double) (this->NbrFluxQuantumSphere + 1)));
  double TmpFactor = (exponentialCoefficient[lz1 - lz2 + this->NbrFluxQuantumSphere] * exponentialCoefficient[lz4 - lz3 + this->NbrFluxQuantumSphere]);
  if (pseudopotentialIndex == 0)
    {
      return TmpFactor / sqrt(2.0 * M_PI);
    }
  if (pseudopotentialIndex == 1)
    {
      return (((double) ((lz1 - lz2) * (lz4 - lz3))) * TmpFactor * Kappa2Factor / sqrt(this->CylinderRatio * ((double) (this->NbrFluxQuantumSphere + 1))));
    }
  return 0.0;
}

