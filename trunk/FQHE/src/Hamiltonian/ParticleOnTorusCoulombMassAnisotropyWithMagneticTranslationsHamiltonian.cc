////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//      coulomb interaction, magnetic translations, and mass anisotropy       //
//                                                                            //
//                        last modification : 14/05/2014                      //
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


#include "Hamiltonian/ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian.h"
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


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// xMomentum = momentum in the x direction (modulo GCD of nbrParticles and maxMomentum)
// ratio = ratio between the width in the x direction and the width in the y direction
// anisotropy = anisotropy parameter alpha (q_g^2 = alpha q_x^2 + q_y^2 / alpha)
// haveCoulomb = flag indicating whether a coulomb term is present
// landauLevel = landauLevel to be simulated (GaAs (>=0) or graphene (<0))
// nbrPseudopotentials = number of pseudopotentials indicated
// pseudopotentials = pseudopotential coefficients
// noWignerEnergy = do not consider the energy contribution from the Wigner crystal 
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian::ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, 
																		 int nbrParticles, int maxMomentum,

																		 int xMomentum, double ratio, double anisotropy,
																		 bool haveCoulomb, int landauLevel, int nbrPseudopotentials, double* pseudopotentials, bool noWignerEnergy,
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
  this->Anisotropy = anisotropy;
  this->InvAnisotropy = 1.0/ this->Anisotropy;
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
    this->WignerEnergy = this->EvaluateWignerCrystalEnergy() / 2.0;
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
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian::~ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian() 
{
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->MaxMomentum);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Q2Anistropy;
  double Precision;
//  cout << "new coef====================================" << m1 << " "  << m2 << " "  << m3 << " "  << m4 << endl;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      Q2Anistropy = this->Ratio * N2 * N2 * this->InvAnisotropy;
      if (N2 != 0.0)
	{
	  Coefficient = this->GetVofQ(PIOnM * Q2Anistropy, PIOnM * Q2Anistropy);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnM * Q2Anistropy, PIOnM * Q2Anistropy); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Q2Anistropy = this->InvRatio * N1 * N1  * this->Anisotropy + this->Ratio * N2 * N2  * this->InvAnisotropy;
	  Precision = 2.0 * this->GetVofQ(PIOnM * Q2Anistropy, PIOnM * Q2Anistropy);
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
      Q2Anistropy = this->Ratio * N2 * N2 * this->InvAnisotropy;
      if (N2 != 0.0)
	{
	  Coefficient =  this->GetVofQ(PIOnM * Q2Anistropy, PIOnM * Q2Anistropy);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnM * Q2Anistropy, PIOnM * Q2Anistropy); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Q2Anistropy = this->InvRatio * N1 * N1  * this->Anisotropy + this->Ratio * N2 * N2  * this->InvAnisotropy;
	  Precision = 2.0 * this->GetVofQ(PIOnM * Q2Anistropy, PIOnM * Q2Anistropy);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
    }
  return (Sum / (2.0 * this->MaxMomentum));
}


// get fourier transform of interaction
//
// Q2_half = one half of q² value
// Q2_halfgaussian = one half of q_g² value for the one body part (i.e. including the anisotropy)

double ParticleOnTorusCoulombMassAnisotropyWithMagneticTranslationsHamiltonian::GetVofQ(double Q2_half, double Q2_halfgaussian)
{
  double Result;
  double Q2 = 2.0*Q2_half;
  if ((this->HaveCoulomb)&&(Q2_half!=0.0))
    {
      //cout << "branch 1 : Ln="<<this->FormFactor.GetValue(Q2_half)<<" Ln2="<<GETSQR(this->FormFactor(Q2_half))<<", exp="<<exp(-Q2_half)<<" 1/Q="<<1.0/sqrt(Q2)<<" ";
      //this->FormFactor.PrintValue(cout, Q2_half)<<" ";
      Result=GETSQR(this->FormFactor(Q2_halfgaussian)) / sqrt(Q2);
    }
  else
    Result=0.0;
  for (int i=0; i<NbrPseudopotentials; ++i)
    if (this->Pseudopotentials[i]!=0.0)
      Result += 2.0*this->Pseudopotentials[i]*this->LaguerreM[i].PolynomialEvaluate(Q2);
  //cout <<"V("<<2*Q2_half<<")="<<Result<<" LL="<<this->LandauLevel<<endl;
  return Result * exp(-Q2_halfgaussian);
}

