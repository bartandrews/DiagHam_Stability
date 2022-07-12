////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                             coulombian interaction                         //
//                                                                            //
//                        last modification : 18/07/2002                      //
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

#include "GeneralTools/StringTools.h"
#include "Hamiltonian/ParticleOnTorusPerturbedCoulombHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Polynomial/SpecialPolynomial.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>

using std::cout;
using std::endl;
using std::ostream;

#define M1_12 0.08333333333333333


// constructor from default data
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// ratio = ratio between the width in the x direction and the width in the y direction
// landauLevel = landauLevel to be simulated
// coulombPrefactor = prefactor multiplying Coulomb terms
// yukawaParameter = modify Coulomb interaction to a Yukawa form with exponential decay exp(-\lambda r) added
// nbrPseudopotentials = number of perturbation in terms of additional pseudopotentials
// pseudopotentials = additional pseudopotential coefficients
// perturbationStrength = prefactor multiplying perturbation terms
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
ParticleOnTorusPerturbedCoulombHamiltonian::ParticleOnTorusPerturbedCoulombHamiltonian(ParticleOnTorus* particles, int nbrParticles, int maxMomentum, double ratio, int landauLevel, double coulombPrefactor, double yukawaParameter, int nbrPseudopotentials, double* pseudopotentials, double perturbationStrength, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->LandauLevel = landauLevel;
  if (this->LandauLevel>=0)
    {
      this->FormFactor=LaguerrePolynomial(this->LandauLevel);
    }
  else
    {
      this->FormFactor=0.5*(LaguerrePolynomial(abs(this->LandauLevel))+LaguerrePolynomial(abs(this->LandauLevel)-1));
    }
  this->CoulombPrefactor=coulombPrefactor;
  this->YukawaParameterSqr=yukawaParameter*yukawaParameter;
  cout << "FormFactor=" << this->FormFactor << endl;
  this->NbrPseudopotentials = nbrPseudopotentials;
  this->Pseudopotentials = new double[this->NbrPseudopotentials];
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    this->Pseudopotentials[i] = perturbationStrength*pseudopotentials[i];
  int max=this->NbrPseudopotentials-1;
  while (this->Pseudopotentials[max] == 0.0 && max>2)
    --max;
  
#ifdef __GMP__
  this->ArraySize = max+1;
  // this->LaguerrePolynomialsStandard = new LaguerrePolynomialRecursion(max+1);
  // this->LaguerrePolynomialsStandard->CheckAgainstGMP(1000.0, 512);
  this->LaguerrePolynomials = new LaguerrePolynomialRecursionAP(this->ArraySize, 128);
  this->TmpLaguerreArray = new mpf_t[this->ArraySize];
  for (int i = 0; i < this->ArraySize; i++)
    {
      mpf_init(TmpLaguerreArray[i]);
      mpf_set_d(TmpLaguerreArray[i], 0);
    }
#else
  this->LaguerrePolynomials = new LaguerrePolynomialRecursion(max+1);
  this->TmpLaguerre.Resize(max+1);
#endif
  
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	    {
	      long TmpMemory = this->FastMultiplicationMemory(memory);
	      PrintMemorySize(cout,TmpMemory)<<endl;
	      if (memory > 0)
	        {
	          this->EnableFastMultiplication();
	        }
	    }
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//
ParticleOnTorusPerturbedCoulombHamiltonian::~ParticleOnTorusPerturbedCoulombHamiltonian() 
{
  delete [] this->Pseudopotentials;
#ifdef __GMP__
  delete this->LaguerrePolynomials;
  delete [] this->TmpLaguerreArray;
#endif
}

// evaluate the numerical coefficient in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient
//
double ParticleOnTorusPerturbedCoulombHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->NbrLzValue);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Precision;
  
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
        {
          Coefficient = this->GetVofQ(PIOnM*Q2, Precision); // Coulomb terms
          cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << N1 << " " << N2 << " " << Coefficient << endl;
        }
      else
        {
          Coefficient = this->GetVofQ(PIOnM*Q2, Precision);  // yields non-zero terms only for non-singular interactions
          cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << N1 << " " << N2 << " " << Coefficient << endl;
        }
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
        {
          Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
          Coefficient += 2.0 * this->GetVofQ(PIOnM*Q2, Precision) * cos (N1 * Factor);
          cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << N1 << " " << N2 << " " << Precision * cos (N1 * Factor) << endl;
          N1 += 1.0;
        }
      Sum += Coefficient;
      N2 += this->NbrLzValue;
    }
  N2 = (double) (m1 - m4 - this->NbrLzValue);
  Coefficient = 1.0;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
        {
          Coefficient =  this->GetVofQ(PIOnM*Q2, Precision);
          cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << N1 << " " << N2 << " " << Coefficient << endl;
        }
      else
        {
	      Coefficient = this->GetVofQ(PIOnM*Q2, Precision); // yields non-zero terms only for non-singular interactions
          cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << N1 << " " << N2 << " " << Coefficient << endl;
	    }
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	    {
	      Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	      Coefficient += 2.0 * this->GetVofQ(PIOnM*Q2, Precision) * cos (N1 * Factor);
          cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << N1 << " " << N2 << " " << Precision * cos (N1 * Factor) << endl;
	      N1 += 1.0;
	    }
      Sum += Coefficient;
      N2 -= this->NbrLzValue;
    }
  cout << "Interaction coefficient for m = (" << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ") = " << (Sum / (2.0 * this->NbrLzValue)) << endl;
  return (Sum / (2.0 * this->NbrLzValue));
}


// get Fourier transform of interaction
//
// Q2_half = one half of q^2 value
// return value = Fourier tranform
//
double ParticleOnTorusPerturbedCoulombHamiltonian::GetVofQ(double Q2_half, double &Precision)
{
  // 1) Coulomb term
  
  double CoulombTerm = 1.0;
  double Q2 = 2.0*Q2_half;
  if ((Q2_half != 0.0) && (this->CoulombPrefactor != 0.0))
    {
      CoulombTerm = this->FormFactor(Q2_half);
      CoulombTerm *= this->CoulombPrefactor*CoulombTerm;
      CoulombTerm /= sqrt(Q2+this->YukawaParameterSqr);
    }
  else CoulombTerm=0.0;
  
  // 2) Additional pseudopotential terms
  
#ifdef __GMP__
  mpf_t TmpInteractionMPF;
  mpf_init_set_d(TmpInteractionMPF, 0.0);
  mpf_t Q2MPF;
  mpf_init_set_d(Q2MPF, Q2);
  this->LaguerrePolynomials->EvaluateLaguerrePolynomials(this->ArraySize, TmpLaguerreArray, Q2MPF);
  mpf_t TmpValue, TmpValue2;  // declare and initialize values for next if statement
  mpf_init(TmpValue);
  mpf_init(TmpValue2);
#else
  double TmpInteraction = 0.0;
  this->LaguerrePolynomials->EvaluateLaguerrePolynomials(TmpLaguerre, Q2);
#endif

  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    {
#ifdef __GMP__
      if (this->Pseudopotentials[i] != 0.0)
        {
          mpf_set_d(TmpValue, 2.0*this->Pseudopotentials[i]);
          mpf_mul(TmpValue2, TmpValue, *(TmpLaguerreArray+i));
          mpf_add(TmpInteractionMPF, TmpInteractionMPF, TmpValue2);
        }
#else
      if (this->Pseudopotentials[i] != 0.0)
        TmpInteraction += 2.0*this->Pseudopotentials[i] * TmpLaguerre[i]; // extra factor of 2 relative to Coulomb terms
#endif
    }
#ifdef __GMP__
    double TmpInteraction = mpf_get_d(TmpInteractionMPF);
#endif
  if (TmpInteraction+CoulombTerm==0.0)
    {
      Precision=1.0;
      return 0.0;
    }
  else
    {
      Precision= exp(-Q2_half) * (TmpInteraction+CoulombTerm);
      return Precision;
    }
}
