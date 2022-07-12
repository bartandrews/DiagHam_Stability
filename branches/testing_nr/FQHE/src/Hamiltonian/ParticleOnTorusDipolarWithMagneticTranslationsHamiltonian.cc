////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                  dipolar interaction and magnetic translations             //
//                                                                            //
//                        last modification : 30/11/2004                      //
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


#include "Hamiltonian/ParticleOnTorusDipolarWithMagneticTranslationsHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"

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
// xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// ratio = ratio between the width in the x direction and the width in the y direction
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusDipolarWithMagneticTranslationsHamiltonian::ParticleOnTorusDipolarWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, 
														     int nbrParticles, 
														     int maxMomentum,  int xMomentum, double ratio, 
														     AbstractArchitecture* architecture, int memory, 
														     char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->EnergyShift = 0.0;
  this->CosinusTable = new double [this->MaxMomentum];
  this->SinusTable = new double [this->MaxMomentum];
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      this->CosinusTable[i] = cos(2.0 * M_PI * this->XMomentum * ((double) i) / ((double) this->MaxMomentum));
      this->SinusTable[i] = sin(2.0 * M_PI * this->XMomentum * ((double) i) / ((double) this->MaxMomentum));
    }
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

ParticleOnTorusDipolarWithMagneticTranslationsHamiltonian::~ParticleOnTorusDipolarWithMagneticTranslationsHamiltonian() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
  delete[] this->M4Value;
  delete[] this->CosinusTable;
  delete[] this->SinusTable;
  if (this->FastMultiplicationFlag == true)
    {
      int ReducedDim = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
      if ((ReducedDim * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
	++ReducedDim;
      for (int i = 0; i < ReducedDim; ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	  delete[] this->InteractionPerComponentNbrTranslation[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
      delete[] this->InteractionPerComponentNbrTranslation;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusDipolarWithMagneticTranslationsHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  if (this->FastMultiplicationFlag == true)
    {
      for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	  delete[] this->InteractionPerComponentNbrTranslation[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
      delete[] this->InteractionPerComponentNbrTranslation;
    }
  this->Particles = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusDipolarWithMagneticTranslationsHamiltonian::ShiftHamiltonian (double shift)
{
}
  
// evaluate all interaction factors
//   

void ParticleOnTorusDipolarWithMagneticTranslationsHamiltonian::EvaluateInteractionFactors()
{
  int Pos = 0;
  int m4;
  double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
  double MaxCoefficient = 0.0;

  if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
    {
     for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if (m4 < 0)
		m4 += this->MaxMomentum;
	      else
		if (m4 >= this->MaxMomentum)
		  m4 -= this->MaxMomentum;
	      if (m3 > m4)
		{
		  TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					 - this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					 - this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		    MaxCoefficient = fabs(TmpCoefficient[Pos]);
		  ++Pos;
		}
	    }
      this->NbrInteractionFactors = 0;
      this->M1Value = new int [Pos];
      this->M2Value = new int [Pos];
      this->M3Value = new int [Pos];
      this->M4Value = new int [Pos];
      this->InteractionFactors = new double [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if (m4 < 0)
		m4 += this->MaxMomentum;
	      else
		if (m4 >= this->MaxMomentum)
		  m4 -= this->MaxMomentum;
	      if (m3 > m4)
		{
		  if  (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
		      this->M1Value[this->NbrInteractionFactors] = m1;
		      this->M2Value[this->NbrInteractionFactors] = m2;
		      this->M3Value[this->NbrInteractionFactors] = m3;
		      this->M4Value[this->NbrInteractionFactors] = m4;
		      ++this->NbrInteractionFactors;
		    }
		  ++Pos;
		}
	    }
    }
  else
    {
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if (m4 < 0)
		m4 += this->MaxMomentum;
	      else
		if (m4 >= this->MaxMomentum)
		  m4 -= this->MaxMomentum;
	      if (m3 > m4)
		{
		  if (m1 != m2)
		    {
		      TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					     + this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		    }
		  else
		    TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					   + this->EvaluateInteractionCoefficient(m1, m2, m4, m3));
		  if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		    MaxCoefficient = fabs(TmpCoefficient[Pos]);
		  ++Pos;
		}
	      else
		if (m3 == m4)
		  {
		    if (m1 != m2)
		      TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		    else
		      TmpCoefficient[Pos] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4);
		    if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		      MaxCoefficient = fabs(TmpCoefficient[Pos]);
		    ++Pos;
		  }
	    }
      this->NbrInteractionFactors = 0;
      this->M1Value = new int [Pos];
      this->M2Value = new int [Pos];
      this->M3Value = new int [Pos];
      this->M4Value = new int [Pos];
      this->InteractionFactors = new double [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      for (int m1 = 0; m1 < this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 < this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if (m4 < 0)
		m4 += this->MaxMomentum;
	      else
		if (m4 >= this->MaxMomentum)
		  m4 -= this->MaxMomentum;
	      if (m3 >= m4)
		{
		  if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
		      this->M1Value[this->NbrInteractionFactors] = m1;
		      this->M2Value[this->NbrInteractionFactors] = m2;
		      this->M3Value[this->NbrInteractionFactors] = m3;
		      this->M4Value[this->NbrInteractionFactors] = m4;
		      ++this->NbrInteractionFactors;
		    }
		  ++Pos;
		}
	    }
    }
  cout << "nbr interaction = " << this->NbrInteractionFactors << endl;
  cout << "====================================" << endl;
  delete[] TmpCoefficient;
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnTorusDipolarWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
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
      Coefficient = exp(- PIOnM * Q2) * sqrt(Q2);
      Precision = Coefficient;
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * exp(- PIOnM * Q2) * sqrt(Q2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 += this->MaxMomentum;
    }
  N2 = (double) (m1 - m4 - this->MaxMomentum);
  Coefficient = Sum;	    
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      Coefficient = exp(- PIOnM * Q2) * sqrt(Q2);
      Precision = Coefficient;
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 *  exp(- PIOnM * Q2) * sqrt(Q2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
    }
  return (Sum / (2.0 * sqrt(2.0 * M_PI * this->MaxMomentum)));
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, ParticleOnTorusDipolarWithMagneticTranslationsHamiltonian& H) 
{
  RealVector TmpV2 (H.Particles->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Particles->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Particles->GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j][i] << "    ";
	}
      Str << endl;
    }
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnTorusDipolarWithMagneticTranslationsHamiltonian& H) 
{
  RealVector TmpV2 (H.Particles->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Particles->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Particles->GetHilbertSpaceDimension() - 1); i++)
    {
      Str << "{";
      for (int j = 0; j < (H.Particles->GetHilbertSpaceDimension() - 1); j++)
	{
	  Str << TmpV[j][i] << ",";
	}
      Str << TmpV[H.Particles->GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Particles->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][H.Particles->GetHilbertSpaceDimension() - 1] << ",";
    }
  Str << TmpV[H.Particles->GetHilbertSpaceDimension() - 1][H.Particles->GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}

