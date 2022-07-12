////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
// antiparallel magnetic field, coulomb interaction and magnetic translations //
//                                                                            //
//                        last modification : 30/05/2018                      //
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


#include "Hamiltonian/ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Architecture/AbstractArchitecture.h"

#include "Polynomial/SpecialPolynomial.h"

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
// xMomentum = momentum in the x direction (modulo GCD of nbrParticles and maxMomentum)
// ratio = ratio between the width in the x direction and the width in the y direction
// layerSeparation = layer separation in units of magnetic length
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian::ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian (ParticleOnTorusWithSpinAndMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio,
 double layerSeparation, AbstractArchitecture* architecture, int memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->Ratio = ratio;  
  this->InvRatio = 1.0 / ratio;
  this->LayerSeparation=layerSeparation;
//   double WignerEnergy = this->EvaluateWignerCrystalEnergy() / 2.0;
  double WignerEnergy = 0.0;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  

  cout << "Wigner Energy = " << WignerEnergy << endl;  
  
  this->EvaluateExponentialFactors();
  this->EvaluateInteractionFactors();

  this->HamiltonianShift = ((double) this->NbrParticles) * WignerEnergy;

  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  int TmpMemory = this->FastMultiplicationMemory(memory);
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
	  cout << endl;
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

ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian::~ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->FastMultiplicationFlag = false;
  this->Particles = (ParticleOnTorusWithSpinAndMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift += shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  long TotalNbrNonZeroInteractionFactors = 0;
  double MaxCoefficient = 0.0;
  this->GetIndices();

//   unsigned L16Mask = (1u<<16)-1;
//   int Pos = 0;
//   int M12Index = 0;
//   int m4;
//   double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
//   double MaxCoefficient = 0.0;  

  if (this->Particles->GetParticleStatistic() == ParticleOnTorusWithSpinAndMagneticTranslations::FermionicStatistic)
    {
      // upup-upup 
      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, 0.0)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, 0.0));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, 0.0)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, 0.0));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactorsupup[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsupup[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      // downdown-downdown
      this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientDownDown(m1, m2, m3, m4, 0.0)
					     + this->EvaluateInteractionCoefficientDownDown(m2, m1, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficientDownDown(m1, m2, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficientDownDown(m2, m1, m3, m4, 0.0));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (this->EvaluateInteractionCoefficientDownDown(m1, m2, m3, m4, 0.0)
					     + this->EvaluateInteractionCoefficientDownDown(m2, m1, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficientDownDown(m1, m2, m4, m3, 0.0)
					     - this->EvaluateInteractionCoefficientDownDown(m2, m1, m3, m4, 0.0));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    {
		      this->InteractionFactorsdowndown[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsdowndown[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      
      
      
      
      
      
      
      // updown-updown
//       this->NbrInterSectorSums = 0;
      this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
      MaxCoefficient = 0.0;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (- this->EvaluateInteractionCoefficientUpDown(m1, m2, m4, m3, this->LayerSeparation)
					     - this->EvaluateInteractionCoefficientUpDown(m2, m1, m3, m4, this->LayerSeparation));
		  if (fabs(TmpCoefficient) > MaxCoefficient)
		    MaxCoefficient = fabs(TmpCoefficient);
		}
	    }
	}
      MaxCoefficient *= MACHINE_PRECISION;
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  double TmpCoefficient   = (- this->EvaluateInteractionCoefficientUpDown(m1, m2, m4, m3, this->LayerSeparation)
					     - this->EvaluateInteractionCoefficientUpDown(m2, m1, m3, m4, this->LayerSeparation));
          
		  if ((fabs(TmpCoefficient) > MaxCoefficient) && (fabs(TmpCoefficient) > MACHINE_PRECISION))
		    {
		      this->InteractionFactorsupdown[i][Index] = TmpCoefficient;
		      ++TotalNbrNonZeroInteractionFactors;
		    }
		  else
		    {
		      this->InteractionFactorsupdown[i][Index] = 0.0;
		    }
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}


    }
  else
    {
      cout << "Bosonic statistics not defined yet for ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian"<<endl;
      exit(1);
    }
  cout << "====================================" << endl;
//  delete[] TmpCoefficient;
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// layerSeparation = separation of layers
// return value = numerical coefficient

double ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian::EvaluateInteractionCoefficientDownDown(int m1, int m2, int m3, int m4, double layerSeparation)
{
    m1 = (this->NbrLzValue - m1) % this->NbrLzValue;
    m2 = (this->NbrLzValue - m2) % this->NbrLzValue;
    m3 = (this->NbrLzValue - m3) % this->NbrLzValue;
    m4 = (this->NbrLzValue - m4) % this->NbrLzValue;
    return this->EvaluateInteractionCoefficient(m1, m2, m3, m4, layerSeparation);
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// layerSeparation = separation of layers
// return value = numerical coefficient

double ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian::EvaluateInteractionCoefficientUpDown(int m1, int m2, int m3, int m4, double layerSeparation)
{
//   if ((m1==((this->NbrLzValue-m2) % this->NbrLzValue))&&(m1==((this->NbrLzValue-m3) % this->NbrLzValue))&&(m1==m4)) 
//       return 0.5;
    
//   if ((m1==m4) && (m2 == m3) && (m1==(this->NbrLzValue - m2) % this->NbrLzValue)) 
//       return 0.5;
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->MaxMomentum);
  double Factor =  - ((double) (m1-m2)) * PIOnM * 2.0;
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
	  Coefficient = this->GetVofQ(PIOnM*Q2, layerSeparation);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = 0.0;// this->GetVofQ(PIOnM*Q2, layerSeparation); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM*Q2, layerSeparation);
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
	  Coefficient =  this->GetVofQ(PIOnM*Q2, layerSeparation);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = 0.0; // this->GetVofQ(PIOnM*Q2, layerSeparation); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM*Q2, layerSeparation);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
    }
  return (Sum / (2.0 * this->MaxMomentum));
}



// get all the indices that should appear in the annihilation/creation operators
//

void ParticleOnTorusCoulombWithSpinAndMagneticTranslationsTimeReversalSymmetricHamiltonian::GetIndices()
{
  this->NbrInterSectorSums = this->NbrLzValue;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;      

      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = 0; m2 <= this->LzMax; ++m2)
	  ++this->NbrInterSectorIndicesPerSum[(m1 - m2 + this->NbrLzValue) % this->NbrLzValue];
      this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  if (this->NbrInterSectorIndicesPerSum[i]  > 0)
	    {
	      this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
	      this->NbrInterSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = 0; m2 <= this->LzMax; ++m2)
	  {
	    int TmpSum = (m1 - m2 + this->NbrLzValue) % this->NbrLzValue;
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = m1;
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = m2;
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	}


  this->NbrIntraSectorSums = this->NbrLzValue;
  this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    this->NbrIntraSectorIndicesPerSum[i] = 0;      
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    int TmpSum = (m1 + m2) % this->NbrLzValue;
	    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = m1;
	    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
	}
    }
  else
    {
      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  {
	    int TmpSum = (m1 + m2) % this->NbrLzValue;
	    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = m1;
	    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
	  }
    }
}
