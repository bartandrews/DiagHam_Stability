////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated to particles on a disk with         //
//             generic interaction defined by its pseudopotential             //
//                                                                            //
//                        last modification : 03/07/2008                      //
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


#include "Hamiltonian/ParticleOnDiskGenericHamiltonian.h"

#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"

#include "MathTools/Complex.h"
#include "MathTools/ClebschGordanDiskCoefficients.h"
#include "MathTools/FactorialCoefficient.h"

#include "Output/MathematicaOutput.h"

#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzMax = maximum angular momentum that a single particle can reach
// pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction, V_m=\int d^2r r^2m V(r) e^(-r^2/8) )
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnDiskGenericHamiltonian::ParticleOnDiskGenericHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzMax, double* pseudoPotential,
								   AbstractArchitecture* architecture, long memory, 
								   bool onDiskCacheFlag, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = false;
  this->Architecture = architecture;
  this->PseudoPotential = new double [this->LzMax + this->NbrLzValue];
  for (int i = 0; i <= (2 * this->LzMax); ++i)
    this->PseudoPotential[i] = pseudoPotential[i];
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
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
	    {
	      cout  << "fast = " << (TmpMemory >> 30) << ".";
	      TmpMemory -= ((TmpMemory >> 30) << 30);
	      TmpMemory *= 100l;
	      TmpMemory >>= 30;
	      if (TmpMemory < 10l)
		cout << "0";
	      cout  << TmpMemory << " Gb ";
	    }
	  if (this->DiskStorageFlag == false)
	    {
	      this->EnableFastMultiplication();
	    }
	  else
	    {
	      char* TmpFileName = this->Architecture->GetTemporaryFileName();
	      this->EnableFastMultiplicationWithDiskStorage(TmpFileName);	      
	      delete[] TmpFileName;
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
  this->L2Operator = 0;
}

// constructor with one body terms
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzMax = maximum angular momentum that a single particle can reach
// architecture = architecture to use for precalculation
// pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction, V_m=\int d^2r r^2m V(r) e^(-r^2/8) )
// oneBodyPotentials = array with the coefficient in front of each one body term (ordered such that the first element corresponds to the one of a+_-s a_-s)
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnDiskGenericHamiltonian::ParticleOnDiskGenericHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzMax, double* pseudoPotential, double* oneBodyPotentials,
								   AbstractArchitecture* architecture, long memory,
								   bool onDiskCacheFlag, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Architecture = architecture;
  this->PseudoPotential = new double [this->LzMax + this->NbrLzValue];
  for (int i = 0; i <= (2 * this->LzMax); ++i)
    this->PseudoPotential[i] = pseudoPotential[i];
  this->OneBodyTermFlag = true;
  this->OneBodyPotentials = new double [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->OneBodyPotentials[i] = oneBodyPotentials[i];
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
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
	    {
	      cout  << "fast = " << (TmpMemory >> 30) << ".";
	      TmpMemory -= ((TmpMemory >> 30) << 30);
	      TmpMemory *= 100l;
	      TmpMemory >>= 30;
	      if (TmpMemory < 10l)
		cout << "0";
	      cout  << TmpMemory << " Gb ";
	    }
	  if (this->DiskStorageFlag == false)
	    {
	      this->EnableFastMultiplication();
	    }
	  else
	    {
	      char* TmpFileName = this->Architecture->GetTemporaryFileName();
	      this->EnableFastMultiplicationWithDiskStorage(TmpFileName);	      
	      delete[] TmpFileName;
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
  this->L2Operator = 0;
}

// destructor
//

ParticleOnDiskGenericHamiltonian::~ParticleOnDiskGenericHamiltonian() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
  delete[] this->PseudoPotential;
  if (this->OneBodyTermFlag == true)
    delete[] this->OneBodyPotentials;
  if (this->FastMultiplicationFlag == true)
    {
       if (this->DiskStorageFlag == false)
	{
	  long MinIndex;
	  long MaxIndex;
	  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
	  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
	  int ReducedDim = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if ((ReducedDim * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	    ++ReducedDim;
	  for (int i = 0; i < ReducedDim; ++i)
	    {
	      delete[] this->InteractionPerComponentIndex[i];
	      delete[] this->InteractionPerComponentCoefficient[i];
	    }
	  delete[] this->InteractionPerComponentIndex;
	  delete[] this->InteractionPerComponentCoefficient;
	}
       else
	 {
	  remove (this->DiskStorageFileName);
	  delete[] this->DiskStorageFileName;
	 }
       delete[] this->NbrInteractionPerComponent;
    }
}

// evaluate all interaction factors
//   

void ParticleOnDiskGenericHamiltonian::EvaluateInteractionFactors()
{
  int Lim;
  int Min;
  int Pos = 0;
  int m4;
  double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
  double MaxCoefficient = 0.0;

  ClebschGordanDiskCoefficients CGCoefficients(this->LzMax);

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  {
	    Lim = m1 + m2;
	    if (Lim > this->LzMax)
	      Lim = this->LzMax;
	    Min = m1 + m2 - this->LzMax;
	    if (Min < 0)
	      Min = 0;
	    int MaxSum = m1 + m2;
	    for (int m3 = Min; m3 <= Lim; ++m3)
	      {
		m4 = m1 + m2 - m3;
		double TmpCoef = 0.0;
		for (int j = 1; j <= MaxSum; j += 2)
		  {
		    TmpCoef -= 4.0 * this->PseudoPotential[j] * (CGCoefficients.GetCoefficient(m1, m2, j) * CGCoefficients.GetCoefficient(m3, m4, j));
		  }
		TmpCoefficient[Pos] = TmpCoef;
		if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		  MaxCoefficient = fabs(TmpCoefficient[Pos]);
		++Pos;
	      }
	  }
      this->NbrInteractionFactors = 0;
      this->M1Value = new int [Pos];
      this->M2Value = new int [Pos];
      this->M3Value = new int [Pos];

      this->InteractionFactors = new double [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;

      this->NbrM12Indices = (this->NbrLzValue * (this->NbrLzValue - 1)) / 2;
      this->M1Value = new int [this->NbrM12Indices];
      this->M2Value = new int [this->NbrM12Indices];
      this->NbrM3Values = new int [this->NbrM12Indices];
      this->M3Values = new int* [this->NbrM12Indices];
      int TotalIndex = 0;
      Pos = 0;


      for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  {
 	    Lim = m1 + m2;
 	    if (Lim > this->LzMax)
 	      Lim = this->LzMax;
 	    Min = m1 + m2 - this->LzMax;
 	    if (Min < 0)
 	      Min = 0;
	    this->M1Value[TotalIndex] = m1;
	    this->M2Value[TotalIndex] = m2;	    
	    this->NbrM3Values[TotalIndex] = 0;
 	    for (int m3 = Min; m3 <= Lim; ++m3)
	      if ((2 * m3) > (m1 + m2))
		++this->NbrM3Values[TotalIndex];
	    if (this->NbrM3Values[TotalIndex] > 0)
	      {
		this->M3Values[TotalIndex] = new int [this->NbrM3Values[TotalIndex]];
		int TmpIndex = 0;
		for (int m3 = Min; m3 <= Lim; ++m3)
		  {
		    if ((2 * m3) > (m1 + m2))
		      {
			this->M3Values[TotalIndex][TmpIndex] = m3;
			this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
			++this->NbrInteractionFactors;
			++TmpIndex;
		      }
		    ++Pos;
		  }
	      }
	    ++TotalIndex;
	  }
    }
  else
    {
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  {
	    Lim = m1 + m2;
	    if (Lim > this->LzMax)
	      Lim = this->LzMax;
	    Min = m1 + m2 - this->LzMax;
	    if (Min < 0)
	      Min = 0;
	    int MaxSum = m1 + m2;
	    for (int m3 = Min; m3 <= Lim; ++m3)
	      {
		m4 = m1 + m2 - m3;
		if (m3 > m4)
		  {
		    if (m1 != m2)
		      {
			double TmpCoef = 0.0;
			for (int j = 0; j <= MaxSum; j += 2)
			  {
			    TmpCoef += 4.0 * this->PseudoPotential[j] * (CGCoefficients.GetCoefficient(m1, m2, j) * CGCoefficients.GetCoefficient(m3, m4, j));
			  }
			TmpCoefficient[Pos] = TmpCoef;
		      }
		    else
		      {
			double TmpCoef = 0.0;
			for (int j = 0; j <= MaxSum; j += 2)
			  {
			    TmpCoef += 2.0 * this->PseudoPotential[j] * (CGCoefficients.GetCoefficient(m1, m2, j) * CGCoefficients.GetCoefficient(m3, m4, j));
			  }
			TmpCoefficient[Pos] = TmpCoef;
		      }
		    if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		      MaxCoefficient = fabs(TmpCoefficient[Pos]);
		    ++Pos;
		  }
		else
		  if (m3 == m4)
		    {
		      if (m1 != m2)
			{
			  double TmpCoef = 0.0;
			  for (int j = 0; j <= MaxSum; j += 2)
			    {
			      TmpCoef += 2.0 * this->PseudoPotential[j] * (CGCoefficients.GetCoefficient(m1, m2, j) * CGCoefficients.GetCoefficient(m3, m4, j));
			    }
			  TmpCoefficient[Pos] = TmpCoef;
			}
		      else
			{
			  double TmpCoef = 0.0;
			  for (int j = 0; j <= MaxSum; j += 2)
			    {
			      TmpCoef += this->PseudoPotential[j] * (CGCoefficients.GetCoefficient(m1, m2, j) * CGCoefficients.GetCoefficient(m3, m4, j));
			    }
			  TmpCoefficient[Pos] = TmpCoef;
			}
		      if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
			MaxCoefficient = fabs(TmpCoefficient[Pos]);
		      ++Pos;
		    }
	      }
	  }
      this->NbrInteractionFactors = 0;
      this->M1Value = new int [Pos];
      this->M2Value = new int [Pos];
      this->M3Value = new int [Pos];
      this->InteractionFactors = new double [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
 
      this->NbrM12Indices = (this->NbrLzValue * (this->NbrLzValue + 1)) / 2;
      this->M1Value = new int [this->NbrM12Indices];
      this->M2Value = new int [this->NbrM12Indices];
      this->NbrM3Values = new int [this->NbrM12Indices];
      this->M3Values = new int* [this->NbrM12Indices];
      int TotalIndex = 0;
      Pos = 0;
      for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	{
	  for (int m2 = 0; m2 < m1; ++m2)
	    {
	      Lim = m1 + m2;
	      if (Lim > this->LzMax)
		Lim = this->LzMax;
	      Min = m1 + m2 - this->LzMax;
	      if (Min < 0)
		Min = 0;
	      this->M1Value[TotalIndex] = m1;
	      this->M2Value[TotalIndex] = m2;	    
	      this->NbrM3Values[TotalIndex] = 0;
	      for (int m3 = Min; m3 <= Lim; ++m3)
		if ((2 * m3) >= (m1 + m2))
		  ++this->NbrM3Values[TotalIndex];
	      if (this->NbrM3Values[TotalIndex] > 0)
		{
		  this->M3Values[TotalIndex] = new int [this->NbrM3Values[TotalIndex]];
		  int TmpIndex = 0;
		  for (int m3 = Min; m3 <= Lim; ++m3)
		    {
		      if ((2 * m3) > (m1 + m2))
			{
			  this->M3Values[TotalIndex][TmpIndex] = m3;
			  this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
			  ++this->NbrInteractionFactors;
			  ++TmpIndex;
			}
		      else
			if ((2 * m3) == (m1 + m2))
			  {
			    this->M3Values[TotalIndex][TmpIndex] = m3;
			    this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * TmpCoefficient[Pos];
			    ++this->NbrInteractionFactors;
			    ++TmpIndex;
			  }			
		      ++Pos;
		    }
		}
	      ++TotalIndex;
	    }
	  Lim = 2 * m1;
	  if (Lim > this->LzMax)
	    Lim = this->LzMax;
	  Min = 2 * m1 - this->LzMax;
	  if (Min < 0)
	    Min = 0;
	  this->M1Value[TotalIndex] = m1;
	  this->M2Value[TotalIndex] = m1;	    
	  this->NbrM3Values[TotalIndex] = 0;
	  for (int m3 = Min; m3 <= Lim; ++m3)
	    if (m3 >= m1)
	      ++this->NbrM3Values[TotalIndex];
	  if (this->NbrM3Values[TotalIndex] > 0)
	    {
	      this->M3Values[TotalIndex] = new int [this->NbrM3Values[TotalIndex]];
	      int TmpIndex = 0;
	      for (int m3 = Min; m3 <= Lim; ++m3)
		{
		  if (m3 > m1)
		    {
		      this->M3Values[TotalIndex][TmpIndex] = m3;
		      this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * TmpCoefficient[Pos];
		      ++this->NbrInteractionFactors;
		      ++TmpIndex;
		    }
		  else
		    if (m3 == m1)
		      {
			this->M3Values[TotalIndex][TmpIndex] = m3;
			this->InteractionFactors[this->NbrInteractionFactors] = 0.25 * TmpCoefficient[Pos];
			++this->NbrInteractionFactors;
			++TmpIndex;
		      }
		  ++Pos;
		}
	    }
	  ++TotalIndex;
	}      
//      for (int m1 = 0; m1 <= this->LzMax; ++m1)
// 	for (int m2 = 0; m2 <= m1; ++m2)
// 	  for (int m3 = 0; m3 <= this->LzMax; ++m3)
// 	    {
// 	      m4 = m1 + m2 - m3;
// 	      if ((m4 >= 0) && (m3 > m4))
// 		{
// 		  if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
// 		    {
// 		      this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
// 		      this->M1Value[this->NbrInteractionFactors] = m1;
// 		      this->M2Value[this->NbrInteractionFactors] = m2;
// 		      this->M3Value[this->NbrInteractionFactors] = m3;
// 		      this->M4Value[this->NbrInteractionFactors] = m4;
// 		      ++this->NbrInteractionFactors;
// 		    }
// 		  ++Pos;
// 		}
// 	    }
    }
  if (this->OneBodyTermFlag == true)
    {
      this->NbrOneBodyInteractionFactors = 0;
      for (int i = 0; i <= this->LzMax; ++i)
	if (this->OneBodyPotentials[i] != 0)
	  ++this->NbrOneBodyInteractionFactors;
      if (this->NbrOneBodyInteractionFactors != 0)
	{
	  double Sign = 1.0;
	  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	    Sign = -1.0;
	  this->OneBodyMValues = new int[this->NbrOneBodyInteractionFactors];
	  this->OneBodyNValues = new int[this->NbrOneBodyInteractionFactors];
	  this->OneBodyInteractionFactors = new double[this->NbrOneBodyInteractionFactors];
	  this->NbrOneBodyInteractionFactors = 0;
	  for (int i = 0; i <= this->LzMax; ++i)
	    if (this->OneBodyPotentials[i] != 0)
	      {
		this->OneBodyMValues[this->NbrOneBodyInteractionFactors] = i;
		this->OneBodyNValues[this->NbrOneBodyInteractionFactors] = i;
		cout << this->OneBodyPotentials[i] << endl;
		this->OneBodyInteractionFactors[this->NbrOneBodyInteractionFactors] = this->OneBodyPotentials[i] * Sign;
		++this->NbrOneBodyInteractionFactors;
	      }	  
	}
      else
	{
	  delete[] this->OneBodyPotentials;
	  this->OneBodyTermFlag = false;
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

double ParticleOnDiskGenericHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  if ((m1 == m2) || (m3 == m4))
    return 0.0;
  FactorialCoefficient Coef;
  Coef.SetToOne();
  if (m2 > 1)
    {
      Coef.PartialFactorialMultiply(m1 + 1, m1 + m2 - 1);
      Coef.FactorialDivide(m2);
    }
  else
    {
      if (m2 == 0)
	Coef /= m1;	
    }
  if (m4 > 1)
    {
      Coef.PartialFactorialMultiply(m3 + 1, m3 + m4 - 1);
      Coef.FactorialDivide(m4);
    }
  else
    {
      if (m4 == 0)
	Coef /= m3;	
    }
  Coef.Power2Divide(2 * (m1 + m2));
  return (sqrt(Coef.GetNumericalValue()) * ((double) ((m2 - m1) * (m3 - m4))) / M_PI);
}

