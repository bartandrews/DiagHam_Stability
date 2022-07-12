////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                         laplacian delta interaction                        //
//                                                                            //
//                        last modification : 14/10/2003                      //
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


#include "Hamiltonian/ParticleOnSphereLaplacianDeltaHamiltonian.h"
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
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereLaplacianDeltaHamiltonian::ParticleOnSphereLaplacianDeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax,
										     AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = false;
  this->OneBodyTermFlag = false;
  this->L2Operator = 0;
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

ParticleOnSphereLaplacianDeltaHamiltonian::~ParticleOnSphereLaplacianDeltaHamiltonian() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
  if (this->FastMultiplicationFlag == true)
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
      delete[] this->NbrInteractionPerComponent;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereLaplacianDeltaHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereLaplacianDeltaHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereLaplacianDeltaHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnSphereLaplacianDeltaHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereLaplacianDeltaHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Particles->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereLaplacianDeltaHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ParticleOnSphereLaplacianDeltaHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> ParticleOnSphereLaplacianDeltaHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all interaction factors
//   

void ParticleOnSphereLaplacianDeltaHamiltonian::EvaluateInteractionFactors()
{
  int Lim;
  int Min;
  int Pos = 0;
  double TmpV = ((((double) this->LzMax) + 1.0) * (((double) this->LzMax) + 1.0)) / 
		 ((0.5 * ((double) this->LzMax)) * ((2.0 * (double) this->LzMax) - 1.0));
  ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
  int J = 2 * this->LzMax - 2;
  int m4;
  double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];

  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  double MaxCoefficient = 0.0;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
	for (int m2 =  -this->LzMax; m2 < m1; m2 += 2)
	  {
	    Lim = m1 + m2 + this->LzMax;
	    if (Lim > this->LzMax)
	      Lim = this->LzMax;
	    Min = m1 + m2 - this->LzMax;
	    if (Min < -this->LzMax)
	      Min = -this->LzMax;
	    for (int m3 = Min; m3 <= Lim; m3 += 2)
	      {
		Clebsch.InitializeCoefficientIterator(m1, m2);
		m4 = m1 + m2 - m3;
		TmpCoefficient[Pos] = 0.0;
		if (((J >> 1) & 1) == Sign)
		  TmpCoefficient[Pos] += TmpV * Clebsch.GetCoefficient(m1, m2, J) * Clebsch.GetCoefficient(m3, m4, J);
		if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		  MaxCoefficient = TmpCoefficient[Pos];
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
      double Factor = - 4.0;
      for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  {
	    Lim = m1 + m2;
	    if (Lim > this->LzMax)
	      Lim = this->LzMax;
	    Min = m1 + m2 - this->LzMax;
	    if (Min < 0)
	      Min = 0;
	    for (int m3 = Min; m3 <= Lim; ++m3)
	      {
		if (((2 * m3) > (m1 + m2)))
		  {
		    this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
		    this->M1Value[this->NbrInteractionFactors] = m1;
		    this->M2Value[this->NbrInteractionFactors] = m2;
		    this->M3Value[this->NbrInteractionFactors] = m3;
		    ++this->NbrInteractionFactors;
		  }
		++Pos;
	      }
	  }
    }
  else
    {
       for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
	for (int m2 =  -this->LzMax; m2 <= m1; m2 += 2)
	  {
	    Lim = m1 + m2 + this->LzMax;
	    if (Lim > this->LzMax)
	      Lim = this->LzMax;
	    Min = m1 + m2 - this->LzMax;
	    if (Min < -this->LzMax)
	      Min = -this->LzMax;
	    for (int m3 = Min; m3 <= Lim; m3 += 2)
	      {
		Clebsch.InitializeCoefficientIterator(m1, m2);
		m4 = m1 + m2 - m3;
		TmpCoefficient[Pos] = 0.0;
		if (((J >> 1) & 1) != Sign)
		  TmpCoefficient[Pos] += TmpV * Clebsch.GetCoefficient(m1, m2, J) * Clebsch.GetCoefficient(m3, m4, J);
		if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		  MaxCoefficient = TmpCoefficient[Pos];
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
      double Factor = 4.0;
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
	      for (int m3 = Min; m3 <= Lim; ++m3)
		{
		  if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      if ((2 * m3) > (m1 + m2))
			{
			  this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
			  this->M1Value[this->NbrInteractionFactors] = m1;
			  this->M2Value[this->NbrInteractionFactors] = m2;
			  this->M3Value[this->NbrInteractionFactors] = m3;
			  ++this->NbrInteractionFactors;
			}
		      else
			if ((2 * m3) == (m1 + m2))
			  {
			    this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * Factor * TmpCoefficient[Pos];
			    this->M1Value[this->NbrInteractionFactors] = m1;
			    this->M2Value[this->NbrInteractionFactors] = m2;
			    this->M3Value[this->NbrInteractionFactors] = m3;
			  ++this->NbrInteractionFactors;
			  }
		    }
		  ++Pos;
		}
	    }	
	  Lim = 2 * m1;
	  if (Lim > this->LzMax)
	    Lim = this->LzMax;
	  Min = 2 * m1 - this->LzMax;
	  if (Min < 0)
	    Min = 0;
	  for (int m3 = Min; m3 <= Lim; ++m3)
	    {
	      if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		{
		  if (m3 > m1)
		    {
		      this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * Factor * TmpCoefficient[Pos];
		      this->M1Value[this->NbrInteractionFactors] = m1;
		      this->M2Value[this->NbrInteractionFactors] = m1;
		      this->M3Value[this->NbrInteractionFactors] = m3;
		      ++this->NbrInteractionFactors;
		    }
		  else
		    if (m3 == m1)
		      {
			this->InteractionFactors[this->NbrInteractionFactors] = 0.25 * Factor * TmpCoefficient[Pos];
			this->M1Value[this->NbrInteractionFactors] = m1;
			this->M2Value[this->NbrInteractionFactors] = m1;
			this->M3Value[this->NbrInteractionFactors] = m3;
			++this->NbrInteractionFactors;
		      }
		}
	      ++Pos;
	    }
	}
    }
  cout << "nbr interaction = " << this->NbrInteractionFactors << endl;
  cout << "====================================" << endl;
  delete[] TmpCoefficient;
}

