////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//       coulomb interaction, a double gate and magnetic translations         //
//                                                                            //
//                        last modification : 30/10/2019                      //
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


#include "Hamiltonian/ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian.h"
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

static double MySqrArg;
#define GETSQR(a) ((MySqrArg=(a)) == 1.0 ? 1.0 : MySqrArg*MySqrArg)


// default constructor
//

ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian::ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// xMomentum = momentum in the x direction (modulo GCD of nbrParticles and maxMomentum)
// ratio = ratio between the width in the x direction and the width in the y direction
// gateDistance = distance of any of the two gates to the electron gas (in lb units)
// haveCoulomb = flag indicating whether a coulomb term is present
// landauLevel = landauLevel to be simulated (GaAs (>=0) or graphene (<0))
// nbrPseudopotentials = number of pseudopotentials indicated
// pseudopotentials = pseudopotential coefficients
// noWignerEnergy = do not consider the energy contribution from the Wigner crystal 
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian::ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, 
														     int nbrParticles, int maxMomentum, 
														     int xMomentum, double ratio, double gateDistance, bool haveCoulomb, int landauLevel, int nbrPseudopotentials, double* pseudopotentials, bool noWignerEnergy,
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
  this->GateDistance = gateDistance;
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

ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian::~ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian() 
{
}

// get fourier transform of interaction
// Q2_half = one half of q² value

double ParticleOnTorusDoubleGatedCoulombWithMagneticTranslationsHamiltonian::GetVofQ(double Q2_half)
{
  double Result;
  double Q2 = 2.0 * Q2_half;
  if ((this->HaveCoulomb) && (Q2_half != 0.0))
    {
      double TmpSqrtQ2 = sqrt(Q2);
      Result = GETSQR(this->FormFactor(Q2_half)) / (TmpSqrtQ2 * tanh(TmpSqrtQ2 * this->GateDistance));
    }
  else
    Result = 0.0;
  for (int i = 0; i < NbrPseudopotentials; ++i)
    if (this->Pseudopotentials[i]!=0.0)
      Result += 2.0*this->Pseudopotentials[i]*this->LaguerreM[i].PolynomialEvaluate(Q2);
  return (Result * exp(-Q2_half));
}

