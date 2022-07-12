////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                  coulomb interaction and magnetic translations             //
//                                                                            //
//                        last modification : 02/10/2003                      //
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


#include "Hamiltonian/ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian.h"
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


#define M1_12 0.08333333333333333

static double MySqrArg;
#define GETSQR(a) ((MySqrArg=(a)) == 1.0 ? 1.0 : MySqrArg*MySqrArg)


// default constructor
//

ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian::ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// xMomentum = momentum in the x direction (modulo GCD of nbrParticles and maxMomentum)
// ratio = ratio between the width in the x direction and the width in the y direction
// haveCoulomb = flag indicating whether a coulomb term is present
// landauLevel = landauLevel to be simulated (GaAs (>=0) or graphene (<0))
// nbrPseudopotentials = number of pseudopotentials indicated
// pseudopotentials = pseudopotential coefficients
// noWignerEnergy = do not consider the energy contribution from the Wigner crystal 
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian::ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, 
														     int nbrParticles, int maxMomentum, 
														     int xMomentum, double ratio, bool haveCoulomb, int landauLevel, int nbrPseudopotentials, double* pseudopotentials, bool noWignerEnergy,
														     AbstractArchitecture* architecture, long memory, 
														     char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->OneBodyInteractionFactors = 0;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->LandauLevel = landauLevel;
  this->NbrPseudopotentials = nbrPseudopotentials;  
  if (this->NbrPseudopotentials>0)
    {
      this->Pseudopotentials = pseudopotentials;
      this->LaguerreM = new Polynomial[NbrPseudopotentials];
      for (int i = 0; i < NbrPseudopotentials; ++i)
	this->LaguerreM[i] = LaguerrePolynomial(i);
    }
  else
    {
      this->Pseudopotentials = NULL;
      this->LaguerreM = NULL;
    }
  this->HaveCoulomb = haveCoulomb;
  if (HaveCoulomb)
    {
      if (this->LandauLevel>=0)
	{
	  // simple coulomb interactions
	  this->FormFactor=LaguerrePolynomial(this->LandauLevel);
	}
      else
	{
	  // coulomb interactions in graphene
	  this->FormFactor=0.5*(LaguerrePolynomial(abs(this->LandauLevel))+LaguerrePolynomial(abs(this->LandauLevel)-1));
	}
    }
  cout << "FormFactor=" << this->FormFactor << endl;
  if ((particles->GetHilbertSpaceDimension() > 0) && (noWignerEnergy == false))
    this->WignerEnergy = ((double) this->NbrParticles) * this->EvaluateWignerCrystalEnergy() / 2.0;
  else 
    this->WignerEnergy = 0.0;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  cout << "Wigner Energy = " << WignerEnergy << endl;  
  this->EvaluateExponentialFactors();
  this->HamiltonianShift = ((double) this->NbrParticles)*WignerEnergy;
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

ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian::~ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian() 
{
  if (this->NbrPseudopotentials>0)
    delete [] this->LaguerreM;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift + this->WignerEnergy;
}
  
// evaluate all interaction factors
//   

void ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  long TotalNbrNonZeroInteractionFactors = 0;
  double MaxCoefficient = 0.0;
  this->GetIndices();
  this->InteractionFactors = new Complex* [this->NbrSectorSums];
  if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
    {
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					   + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					   - this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					   - this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors[i][Index] = TmpCoefficient;
		      TotalNbrNonZeroInteractionFactors++;
		    }
		  else
		    {
		      this->InteractionFactors[i][Index] = 0.0;
		    }
		  TotalNbrInteractionFactors++;
		  ++Index;
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
	      int m1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					   + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					   + this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					   + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (m3 == m4)
		    TmpCoefficient *= 0.5;
		  if (m1 == m2)
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					   + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					   + this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					   + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (m3 == m4)
		    TmpCoefficient *= 0.5;
		  if (m1 == m2)
		    TmpCoefficient *= 0.5;
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactors[i][Index] = TmpCoefficient;
		      TotalNbrNonZeroInteractionFactors++;
		    }
		  else
		    {
		      this->InteractionFactors[i][Index] = 0.0;
		    }
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

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->MaxMomentum);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Precision;
//  cout << "new coef====================================" << m1 << " "  << m2 << " "  << m3 << " "  << m4 << endl;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = this->GetVofQ(PIOnM*Q2);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnM*Q2); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM*Q2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->MaxMomentum;
    }
  N2 = (double) (m1 - m4 - this->MaxMomentum);
  Coefficient = 1.0;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient =  this->GetVofQ(PIOnM*Q2);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnM*Q2); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM*Q2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
    }
  return (Sum / (2.0 * this->MaxMomentum));
}


// get fourier transform of interaction
// Q2_half = one half of q² value

double ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian::GetVofQ(double Q2_half)
{
  double Result;
  double Q2=2.0*Q2_half;
  if ((this->HaveCoulomb)&&(Q2_half!=0.0))
    {
      //cout << "branch 1 : Ln="<<this->FormFactor.GetValue(Q2_half)<<" Ln2="<<GETSQR(this->FormFactor(Q2_half))<<", exp="<<exp(-Q2_half)<<" 1/Q="<<1.0/sqrt(Q2)<<" ";
      //this->FormFactor.PrintValue(cout, Q2_half)<<" ";
      Result=GETSQR(this->FormFactor(Q2_half)) / sqrt(Q2);
    }
  else
    Result=0.0;
  for (int i=0; i<NbrPseudopotentials; ++i)
    if (this->Pseudopotentials[i]!=0.0)
      Result += 2.0*this->Pseudopotentials[i]*this->LaguerreM[i].PolynomialEvaluate(Q2);
  //cout <<"V("<<2*Q2_half<<")="<<Result<<" LL="<<this->LandauLevel<<endl;
  return Result * exp(-Q2_half);
}

// evaluate Wigner crystal energy per particle
//
// return value = Wigner crystal energy per particle

double ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian::EvaluateWignerCrystalEnergy ()
{
  double TmpRatio = M_PI * this->Ratio;
  double TmpInvRatio = M_PI * this->InvRatio;
  double Energy = this->MisraFunction(-0.5, TmpRatio);
  double Precision = Energy;
  int L1 = 2;
  while ((Energy + Precision) > Energy)
    {
      Precision = this->MisraFunction(-0.5, TmpRatio * L1 * L1);
      Energy += Precision;
      ++L1;
    }
  Energy *= 2.0;
  int L2 = 1;
  double PartialEnergy = Energy;
  while ((PartialEnergy + Energy) > Energy)
    {
      PartialEnergy = 2.0 * this->MisraFunction(-0.5, TmpInvRatio * L2 * L2);
      Precision = PartialEnergy;
      L1 = 1;
      while (((PartialEnergy + Precision) > PartialEnergy))// && ((fabs(PartialEnergy - Precision) + Energy) > Energy))
	{
	  Precision = 4.0 * this->MisraFunction(-0.5, TmpRatio * L1 * L1 + TmpInvRatio * L2 * L2);
	  PartialEnergy += Precision;
	  ++L1;	  
	}
      Energy += PartialEnergy;
      ++L2;
    }
  return 2.0 * (Energy - 2.0) / sqrt (2.0 * M_PI * this->MaxMomentum);
}

// evaluate Misra function (integral of t^n exp (-xt) between 1 and +inf)
//
// n = index of the Misra function
// x = point where the function has to be evaluated (> 0)
// return value = value of the n-Misra function at x

double ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian::MisraFunction (double n, double x)
{
  int Count=0;
  int NbrSubdivision = 100000;
  double PreviousSum = this->PartialMisraFunction(n, x, 0.0, 1.0, NbrSubdivision);
  double NewSum = PreviousSum;
  PreviousSum *= 2.0;
  while (((fabs(PreviousSum - NewSum) / PreviousSum) > MACHINE_PRECISION) && (Count<5))
    {
      if ((fabs(PreviousSum - NewSum) / PreviousSum) < 1e-11)
	++Count;
      PreviousSum = NewSum;
      NbrSubdivision += 10000;
      NewSum = this->PartialMisraFunction(n, x, 0.0, 1.0, NbrSubdivision);
      //cout << " PreviousSum = " << PreviousSum << "   NewSum = " << NewSum << "  diff="<<PreviousSum-NewSum<<endl;
    }
  return 2.0 * (sqrt(M_PI * 0.25 / x) - NewSum);
}

// evaluate part of the integral needed in the Misra function (integral of t^n exp (-xt) between min and max)
//
// n = index of the Misra function
// x = point where the function has to be evaluated (> 0)
// min = lower bound of the integral
// max = upper bound of the integral
// nbrSubdivision = number of subdivision used for the integral
// return value = value of the integral

double ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian::PartialMisraFunction (double n, double x, double min, double max, int nbrSubdivision)
{
  double Sum = 0.0;
  x *= -1.0;
  --nbrSubdivision;
  max  = (max - min) / ((double) nbrSubdivision);
  Sum += (0.5 + M1_12 * 2.0 * x * min * max )* exp(min * min * x);
  min += max;
  --nbrSubdivision;
  while (nbrSubdivision > 0)
    {
      Sum += exp(min * min * x);
      min += max;
      --nbrSubdivision;
    }
  Sum += (0.5 - M1_12 * 2.0 * x * min * max) * exp(min * min * x);
  Sum *= max;
  return Sum;
}

