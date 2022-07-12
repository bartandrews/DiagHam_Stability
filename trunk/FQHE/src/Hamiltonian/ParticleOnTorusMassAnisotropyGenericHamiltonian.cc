////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//        generic two body interaction with an anisotropic mass term          //
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


#include "Hamiltonian/ParticleOnTorusMassAnisotropyGenericHamiltonian.h"
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


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = maximum Lz value reached by a particle in the state
// ratio = ratio between the width in the x direction and the width in the y direction
// anisotropy = anisotropy parameter alpha (q_g^2 = alpha q_x^2 + q_y^2 / alpha)
// nbrPseudopotentials = number of pseudopotentials
// pseudopotentials = pseudopotential coefficients
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusMassAnisotropyGenericHamiltonian::ParticleOnTorusMassAnisotropyGenericHamiltonian(ParticleOnTorus* particles, int nbrParticles, int maxMomentum, double ratio, 
												   double anisotropy, int nbrPseudopotentials, double* pseudopotentials,
												   AbstractArchitecture* architecture, long memory, char* precalculationFileName)
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
  this->HamiltonianShift = 0.0;
  this->OneBodyTermFlag = false;
  this->OneBodyInteractionFactors = 0;
  this->FilterInteractionFlag = false;
  this->Anisotropy = anisotropy;
  this->InvAnisotropy = 1.0/ this->Anisotropy;
  this->NbrPseudopotentials = nbrPseudopotentials;
  this->Pseudopotentials = new double[this->NbrPseudopotentials];
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    this->Pseudopotentials[i] = pseudopotentials[i];
  this->LaguerrePolynomials =new Polynomial[this->NbrPseudopotentials];
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    this->LaguerrePolynomials[i] = LaguerrePolynomial(i);

  this->EvaluateInteractionFactors();

  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  if (TmpMemory < 1024)
	    cout  << "fast = " <<  TmpMemory << "b ";
	  else
	    if (TmpMemory < (1 << 20))
	      cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	    else
	      if (TmpMemory < (1 << 30))
		cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	      else
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
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

ParticleOnTorusMassAnisotropyGenericHamiltonian::~ParticleOnTorusMassAnisotropyGenericHamiltonian() 
{
}


// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnTorusMassAnisotropyGenericHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->NbrLzValue);
  double TwoPIOnM = 2.0 * M_PI / ((double) this->NbrLzValue);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Q2Anistropy;
  double Precision;
  double TmpInteraction;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      Q2Anistropy = this->Ratio * N2 * N2 * this->InvAnisotropy;
      if (N2 != 0.0)
	{
	  TmpInteraction = 0.0;
	  for (int i = 0; i < this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(TwoPIOnM * Q2Anistropy);
	  Coefficient = exp(- PIOnM * Q2Anistropy) * TmpInteraction;
	  Precision = Coefficient;
	}
      else
	{
	  Precision = 1.0;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Q2Anistropy = this->InvRatio * N1 * N1  * this->Anisotropy + this->Ratio * N2 * N2  * this->InvAnisotropy;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(TwoPIOnM * Q2Anistropy);
	  Precision = 2.0 * exp(- PIOnM * Q2Anistropy) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->NbrLzValue;
    }
  N2 = (double) (m1 - m4 - this->NbrLzValue);
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      Q2Anistropy = this->Ratio * N2 * N2  * this->InvAnisotropy;
      if (N2 != 0.0)
	{
	  TmpInteraction = 0.0;
	  for (int i=0; i< this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(TwoPIOnM * Q2Anistropy);
	  Coefficient = exp(- PIOnM * Q2Anistropy) * TmpInteraction;
	  Precision = Coefficient;
	}
      else
	{
	  Precision = 1.0;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(0.0);
	  Coefficient = TmpInteraction;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Q2Anistropy = this->InvRatio * N1 * N1  * this->Anisotropy + this->Ratio * N2 * N2  * this->InvAnisotropy;
	  TmpInteraction = 0.0;
	  for (int i = 0; i < this->NbrPseudopotentials; ++i)
	    if (this->Pseudopotentials[i] != 0.0)
	      TmpInteraction += this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(TwoPIOnM * Q2Anistropy);
	  Precision = 2.0 *  exp(- PIOnM * Q2Anistropy) * TmpInteraction;
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->NbrLzValue;
    }
  return (Sum / this->NbrLzValue);
}

