////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                       coulombian and delta interaction                     //
//                                                                            //
//                        last modification : 24/03/2003                      //
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


#include "Hamiltonian/ParticleOnSphereCoulombDeltaHamiltonian.h"
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
// ratio = ratio between coulombian interaction and delta interaction
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnSphereCoulombDeltaHamiltonian::ParticleOnSphereCoulombDeltaHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax,
										 double ratio, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Ratio = ratio;
  this->HamiltonianShift = -0.5 * ((double) (this->NbrParticles * this->NbrParticles)) / (0.5 * ((double) this->LzMax));
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = false;
  this->OneBodyTermFlag = false;
  this->L2Operator = 0;
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

// destructor
//

ParticleOnSphereCoulombDeltaHamiltonian::~ParticleOnSphereCoulombDeltaHamiltonian() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
}

// evaluate all interaction factors
//   

void ParticleOnSphereCoulombDeltaHamiltonian::EvaluateInteractionFactors()
{
  int Lim;
  int Min;
  int Pos = 0;
//  int NbrNonZero = 0;
  double* TmpV = new double [this->NbrLzValue];
  FactorialCoefficient Coef;
  for (int j = 0; j <= this->LzMax; ++j)
    {
      Coef.SetToOne();
      Coef.PartialFactorialMultiply(this->LzMax - j + 1, 2 * this->LzMax - 2 * j);
      Coef.FactorialDivide(this->LzMax - j);
      Coef.PartialFactorialMultiply(this->LzMax + j + 2, 2 * this->LzMax + 2 * j + 2);
      Coef.FactorialDivide(this->LzMax + j + 1);
      Coef.PartialFactorialDivide(this->LzMax + 2, 2 * this->LzMax + 2);
      Coef.FactorialMultiply(this->LzMax + 1);
      Coef.PartialFactorialDivide(this->LzMax + 2, 2 * this->LzMax + 2);
      Coef.FactorialMultiply(this->LzMax + 1);
      TmpV[j] = Coef.GetNumericalValue();
      cout << "TmpV[" << j << "] = " << TmpV[j]<< endl;
    }
  double TmpVDelta = (((double) this->LzMax) + 1.0);
  TmpVDelta = (TmpVDelta * TmpVDelta) / (4.0 * M_PI * (2.0 * ((double) this->LzMax) + 1.0));
  for (int j = 0; j <= this->LzMax; ++j)
    {
      TmpV[j] *= this->Ratio * TmpVDelta / TmpV[this->LzMax];
    }
  TmpV[this->LzMax] = TmpVDelta;
  ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
  int J;
  int m4;
  double ClebschCoef;
  double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];

  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  double MaxCoefficient = 0.0;
//  double Factor = 4.0 / sqrt (0.5 * ((double) this->LzMax)) * this->Ratio;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
	for (int m2 =  -this->LzMax; m2 < m1; m2 += 2)//= this->LzMax; m2 += 2)// m1; m2 += 2)
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
		while (Clebsch.Iterate(J, ClebschCoef))
		  {
		    if (((J >> 1) & 1) == Sign)
		      TmpCoefficient[Pos] += TmpV[J >> 1] * ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		  }
/*		TmpCoefficient[Pos] *= Factor;
		J = 2 * this->LzMax;
		if (((J >> 1) & 1) == Sign)
		  TmpCoefficient[Pos] += TmpVDelta * Clebsch.GetCoefficient(m1, m2, J) * Clebsch.GetCoefficient(m3, m4, J);
		J -= 2;
		if (((J >> 1) & 1) == Sign)
		  TmpCoefficient[Pos] += TmpVDelta * this->Ratio * Clebsch.GetCoefficient(m1, m2, J) * Clebsch.GetCoefficient(m3, m4, J);*/
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
      double Factor = - 4.0 / sqrt (0.5 * ((double) this->LzMax));
      for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)//this->NbrLzValue; ++m2)// m1; ++m2)
	  {
	    Lim = m1 + m2;
	    if (Lim > this->LzMax)
	      Lim = this->LzMax;
	    Min = m1 + m2 - this->LzMax;
	    if (Min < 0)
	      Min = 0;
	    for (int m3 = Min; m3 <= Lim; ++m3)
	      {
		if ((fabs(TmpCoefficient[Pos]) > MaxCoefficient) && ((2 * m3) > (m1 + m2)))
	      {
		this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
		this->M1Value[this->NbrInteractionFactors] = m1;
		this->M2Value[this->NbrInteractionFactors] = m2;
		this->M3Value[this->NbrInteractionFactors] = m3;
		/*		cout << this->M1Value[this->NbrInteractionFactors] << " " << this->M2Value[this->NbrInteractionFactors] 
				<< " " << this->M3Value[this->NbrInteractionFactors] 
				<< " " << this->InteractionFactors[this->NbrInteractionFactors] << endl;*/
		++this->NbrInteractionFactors;
	      }
		++Pos;
	      }
	  }
    }
  else
    {
       for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
	for (int m2 =  -this->LzMax; m2 <= m1; m2 += 2) // this->LzMax; m2 += 2)// 
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
		while (Clebsch.Iterate(J, ClebschCoef))
		  {
		    if (((J >> 1) & 1) != Sign)
		      TmpCoefficient[Pos] += TmpV[J >> 1] * ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		  }
/*		TmpCoefficient[Pos] *= Factor;
		J = 2 * this->LzMax;
		if (((J >> 1) & 1) != Sign)
		  TmpCoefficient[Pos] += this->Ratio * TmpVDelta * Clebsch.GetCoefficient(m1, m2, J) * Clebsch.GetCoefficient(m3, m4, J);
		J -= 4;
		if (((J >> 1) & 1) != Sign)
		  TmpCoefficient[Pos] += TmpVDelta * Clebsch.GetCoefficient(m1, m2, J) * Clebsch.GetCoefficient(m3, m4, J);*/
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
      for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	{
	  for (int m2 = 0; m2 < m1; ++m2)//this->NbrLzValue; ++m2)//
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
			  this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
			  this->M1Value[this->NbrInteractionFactors] = m1;
			  this->M2Value[this->NbrInteractionFactors] = m2;
			  this->M3Value[this->NbrInteractionFactors] = m3;
			  ++this->NbrInteractionFactors;
			}
		      else
			if ((2 * m3) == (m1 + m2))
			  {
			    this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * TmpCoefficient[Pos];
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
		      this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * TmpCoefficient[Pos];
		      this->M1Value[this->NbrInteractionFactors] = m1;
		      this->M2Value[this->NbrInteractionFactors] = m1;
		      this->M3Value[this->NbrInteractionFactors] = m3;
		      ++this->NbrInteractionFactors;
		    }
		  else
		    if (m3 == m1)
		      {
			this->InteractionFactors[this->NbrInteractionFactors] = 0.25 * TmpCoefficient[Pos];
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
  delete[] TmpV;
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, ParticleOnSphereCoulombDeltaHamiltonian& H) 
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

MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnSphereCoulombDeltaHamiltonian& H) 
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

