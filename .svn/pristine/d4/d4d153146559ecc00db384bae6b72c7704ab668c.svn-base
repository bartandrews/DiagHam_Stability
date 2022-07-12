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


#include "Hamiltonian/ParticleOnTorusNematicParameterWithMagneticTranslationsHamiltonian.h"
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

ParticleOnTorusNematicParameterWithMagneticTranslationsHamiltonian::ParticleOnTorusNematicParameterWithMagneticTranslationsHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// xMomentum = momentum in the x direction (modulo GCD of nbrParticles and maxMomentum)
// ratio = ratio between the width in the x direction and the width in the y direction
// length = characteristic length (in magnetic lentgh unit) used in the order parameter, i.e. (cos q_xl - cos q_yl)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusNematicParameterWithMagneticTranslationsHamiltonian::ParticleOnTorusNematicParameterWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, 
														     int nbrParticles, int maxMomentum, 
														     int xMomentum, double ratio, double length, 
														     AbstractArchitecture* architecture, long memory, 
														     char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->NbrParticles = nbrParticles;
  this->Length = length;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->OneBodyInteractionFactors = 0;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->LandauLevel = 0;
  this->NbrPseudopotentials = 0;  
  this->Pseudopotentials = NULL;
  this->LaguerreM = NULL;
  this->HaveCoulomb = false;
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

ParticleOnTorusNematicParameterWithMagneticTranslationsHamiltonian::~ParticleOnTorusNematicParameterWithMagneticTranslationsHamiltonian() 
{
}


// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnTorusNematicParameterWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->MaxMomentum);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double FactorQx =  sqrt (2.0 * PIOnM * this->InvRatio) * this->Length;
  double FactorQy =  sqrt (2.0 * PIOnM * this->Ratio) * this->Length;
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
	  Coefficient = this->GetVofQ(PIOnM * Q2);
	  Precision = Coefficient;
	  Coefficient *= (1.0 - cos(FactorQy * N2));
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnM * Q2) * (1.0 - cos(FactorQy * N2)); 
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM * Q2);
	  Coefficient += Precision * cos (N1 * Factor) * (cos (FactorQx * N1) - cos(FactorQy * N2));
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
	  Coefficient =  this->GetVofQ(PIOnM * Q2); 
	  Precision = Coefficient;
	  Coefficient *= (1.0 - cos(FactorQy * N2));
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnM * Q2) * (1.0 - cos(FactorQy * N2));
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM * Q2);
	  Coefficient += Precision * cos (N1 * Factor) * (cos (FactorQx * N1) - cos(FactorQy * N2));
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

double ParticleOnTorusNematicParameterWithMagneticTranslationsHamiltonian::GetVofQ(double Q2_half)
{
  double Result = 2.0;
  double Q2 = 2.0 * Q2_half;
  return (Result * exp(-Q2_half));
}

