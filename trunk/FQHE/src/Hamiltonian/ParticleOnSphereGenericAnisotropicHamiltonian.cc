////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//             generic interaction defined by its pseudopotential             //
//                                                                            //
//                        last modification : 03/06/2004                      //
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


#include "Hamiltonian/ParticleOnSphereGenericAnisotropicHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"

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
// pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction)
// anisotropicPseudoPotentialAlpha2 = array with the anisotropic pseudo-potentials alpha=2 (ordered such that the first element corresponds to the delta interaction)
// anisotropicPseudoPotentialAlpha4 = array with the anisotropic pseudo-potentials alpha=4 (ordered such that the first element corresponds to the delta interaction)
// chirality = chirality (+/-1) of the anisotropic pseudopotential
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
// hermitianFlag = flag to indicate if hermitian symmetry of Hamiltonian shall be used
ParticleOnSphereGenericAnisotropicHamiltonian::ParticleOnSphereGenericAnisotropicHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, double* pseudoPotential, double* anisotropicPseudoPotentialAlpha2, double* anisotropicPseudoPotentialAlpha4, int chirality, double l2Factor,
								       AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag,
								       char* precalculationFileName, bool hermitianFlag)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = false;
  this->Architecture = architecture;

  this->PseudoPotential = new double [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->PseudoPotential[i] = pseudoPotential[this->LzMax - i];

  this->AnisotropicPseudoPotentialAlpha2 = new double [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->AnisotropicPseudoPotentialAlpha2[i] = anisotropicPseudoPotentialAlpha2[i];

  this->AnisotropicPseudoPotentialAlpha4 = new double [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->AnisotropicPseudoPotentialAlpha4[i] = anisotropicPseudoPotentialAlpha4[i];
 
  this->Chirality = chirality;
  if (this->Chirality != 0)
     cout << "Using chiral anisotropic pseudopotential " << this->Chirality << endl;

  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
  if (hermitianFlag)
    this->HermitianSymmetrizeInteractionFactors();
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

  if (l2Factor != 0.0)
    { 
      this->L2Operator = new ParticleOnSphereL2Hamiltonian(this->Particles, this->NbrParticles, this->LzMax,
							   this->Particles->GetLzValue() , this->Architecture, l2Factor,
							   /*memory*/ -1, /* fixedLz */ true, /* onDiskCacheFlag */ false,
							   /*precalculationFileName */ 0, hermitianFlag);
    }
  else
    {
      this->L2Operator = 0;
    }
}

// constructor with one body terms
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// architecture = architecture to use for precalculation
// pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction)
// oneBodyPotentials = array with the coefficient in front of each one body term (ordered such that the first element corresponds to the one of a+_-s a_-s)
// anisotropicPseudoPotentialAlpha2 = array with the anisotropic pseudo-potentials alpha=2 (ordered such that the first element corresponds to the delta interaction)
// anisotropicPseudoPotentialAlpha4 = array with the anisotropic pseudo-potentials alpha=4 (ordered such that the first element corresponds to the delta interaction)
// chirality = chirality (+/-1) of the anisotropic pseudopotential
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
// hermitianFlag = flag to indicate if hermitian symmetry of Hamiltonian shall be used
ParticleOnSphereGenericAnisotropicHamiltonian::ParticleOnSphereGenericAnisotropicHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
								       double* pseudoPotential, double* oneBodyPotentials, double* anisotropicPseudoPotentialAlpha2, double* anisotropicPseudoPotentialAlpha4, int chirality, double l2Factor,
								       AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag,
								       char* precalculationFileName, bool hermitianFlag)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Architecture = architecture;

  this->PseudoPotential = new double [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->PseudoPotential[i] = pseudoPotential[this->LzMax - i];

  this->AnisotropicPseudoPotentialAlpha2 = new double [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->AnisotropicPseudoPotentialAlpha2[i] = anisotropicPseudoPotentialAlpha2[i];

  this->AnisotropicPseudoPotentialAlpha4 = new double [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->AnisotropicPseudoPotentialAlpha4[i] = anisotropicPseudoPotentialAlpha4[i];

  this->Chirality = chirality;
  if (this->Chirality != 0)
     cout << "Using chiral anisotropic pseudopotential " << this->Chirality << endl;

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
  if (hermitianFlag)
    this->HermitianSymmetrizeInteractionFactors();
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

  if (l2Factor != 0.0)
    {
      this->L2Operator = new ParticleOnSphereL2Hamiltonian(this->Particles, this->NbrParticles, this->LzMax,
							   this->Particles->GetLzValue() , this->Architecture, l2Factor,
							   /*memory*/ -1, /* fixedLz */ true, /* onDiskCacheFlag */ false,
							   /*precalculationFileName */ 0, hermitianFlag);
    }
  else
    {
      this->L2Operator = 0;
    }
}

// destructor
//

ParticleOnSphereGenericAnisotropicHamiltonian::~ParticleOnSphereGenericAnisotropicHamiltonian() 
{
  delete[] this->InteractionFactors;  
  delete[] this->M1Value;
  delete[] this->M2Value;
  if (M3Value!=0)
    {
      delete[] this->M3Value;
    }
  else
    delete[] this->M3Value;
  delete[] this->PseudoPotential;
  delete[] this->AnisotropicPseudoPotentialAlpha2;
  delete[] this->AnisotropicPseudoPotentialAlpha4;
  if (this->OneBodyTermFlag == true)
    delete[] this->OneBodyPotentials;
  if (this->L2Operator != 0)
    delete this->L2Operator;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereGenericAnisotropicHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereGenericAnisotropicHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereGenericAnisotropicHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnSphereGenericAnisotropicHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereGenericAnisotropicHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
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

Complex ParticleOnSphereGenericAnisotropicHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ParticleOnSphereGenericAnisotropicHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> ParticleOnSphereGenericAnisotropicHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all interaction factors
//   

void ParticleOnSphereGenericAnisotropicHamiltonian::EvaluateInteractionFactors()
{
  int Lim;
  int Min;
  int Pos = 0;
  ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
  int J;
  int m4;
  double ClebschCoef;
  double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
  for (int i = 0; i < (this->NbrLzValue * this->NbrLzValue * this->NbrLzValue); i++) TmpCoefficient[i] = 0.0;

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
		m4 = m1 + m2 - m3;

		TmpCoefficient[Pos] = 0.0;

		// add the standard pseudopotential

		J = 2 * this->LzMax - 2;
		Clebsch.InitializeCoefficientIterator(m1, m2);
		while (Clebsch.Iterate(J, ClebschCoef))
		  {
		    if (((J >> 1) & 1) == Sign)
		    {
		      TmpCoefficient[Pos] += this->PseudoPotential[J >> 1] * ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
// 		      cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << (this->LzMax - (J >> 1)) << " " << TmpCoefficient[Pos] << endl;
		    }
		  }


		// add the anisotropic pseudopotential alpha = 2
		double TmpAnisoPos = 0.0;		
		double TmpAnisoNeg = 0.0;		
		int Alpha = 2;

		int Mom = 1;
		while (Mom <= this->LzMax)
		  {
		    if (this->AnisotropicPseudoPotentialAlpha2[Mom] != 0.0)
		      {
			int L12 = this->LzMax - Mom;
			int L34 = this->LzMax - Mom - Alpha;
			int tmpmin;
			if (L12 < L34)	
				tmpmin = L12;
			else
				tmpmin = L34;

			for (int ML = -tmpmin; ML <= tmpmin; ML++)
				if ( (2*ML==(m1+m2)) && (2*ML==(m3+m4)) )
					TmpAnisoPos +=  this->AnisotropicPseudoPotentialAlpha2[Mom]  * Clebsch.GetCoefficient(m1, m2, 2*L12) * Clebsch.GetCoefficient(m3, m4, 2*L34);

			L12 = this->LzMax - Mom - Alpha;
			L34 = this->LzMax - Mom;
			if (L12 < L34)	
				tmpmin = L12;
			else
				tmpmin = L34;

			for (int ML = -tmpmin; ML <= tmpmin; ML++)
				if ( (2*ML==(m1+m2)) && (2*ML==(m3+m4)) )
					TmpAnisoNeg += this->AnisotropicPseudoPotentialAlpha2[Mom] * Clebsch.GetCoefficient(m1, m2, 2*L12) * Clebsch.GetCoefficient(m3, m4, 2*L34);
		     }
		    Mom = Mom + 2;
		  }

		if (this->Chirality == 1)
		    TmpCoefficient[Pos] += TmpAnisoPos;
		else if (this->Chirality == -1)
		    TmpCoefficient[Pos] += TmpAnisoNeg;
                else                      
	    	    TmpCoefficient[Pos] += 0.5 * (TmpAnisoPos + TmpAnisoNeg);


		// add the anisotropic pseudopotential alpha = 4
		TmpAnisoPos = 0.0;		
		TmpAnisoNeg = 0.0;		
		Alpha = 4;

		Mom = 1;
		while (Mom <= this->LzMax)
		  {
		    if (this->AnisotropicPseudoPotentialAlpha4[Mom] != 0.0)
		      {
			int L12 = this->LzMax - Mom;
			int L34 = this->LzMax - Mom - Alpha;
			int tmpmin;
			if (L12 < L34)	
				tmpmin = L12;
			else
				tmpmin = L34;

			for (int ML = -tmpmin; ML <= tmpmin; ML++)
				if ( (2*ML==(m1+m2)) && (2*ML==(m3+m4)) )
					TmpAnisoPos +=  this->AnisotropicPseudoPotentialAlpha4[Mom]  * Clebsch.GetCoefficient(m1, m2, 2*L12) * Clebsch.GetCoefficient(m3, m4, 2*L34);

			L12 = this->LzMax - Mom - Alpha;
			L34 = this->LzMax - Mom;
			if (L12 < L34)	
				tmpmin = L12;
			else
				tmpmin = L34;

			for (int ML = -tmpmin; ML <= tmpmin; ML++)
				if ( (2*ML==(m1+m2)) && (2*ML==(m3+m4)) )
					TmpAnisoNeg += this->AnisotropicPseudoPotentialAlpha4[Mom] * Clebsch.GetCoefficient(m1, m2, 2*L12) * Clebsch.GetCoefficient(m3, m4, 2*L34);
		     }
		    Mom = Mom + 2;
		  }

		if (this->Chirality == 1)
		    TmpCoefficient[Pos] += TmpAnisoPos;
		else if (this->Chirality == -1)
		    TmpCoefficient[Pos] += TmpAnisoNeg;
                else                      
	    	    TmpCoefficient[Pos] += 0.5 * (TmpAnisoPos + TmpAnisoNeg);

		//cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << TmpCoefficient[Pos] << endl;

		if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		  MaxCoefficient = TmpCoefficient[Pos];
		++Pos;
	      }
	  }
      this->NbrInteractionFactors = 0;      
      this->InteractionFactors = new double [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      double Factor = - 4.0;
//       this->NbrM3Values=0;
//       this->M1Value = new int [Pos];
//       this->M2Value = new int [Pos];
//       this->M3Value = new int [Pos];      
//       for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
// 	for (int m2 = 0; m2 < m1; ++m2)
// 	  {
// 	    Lim = m1 + m2;
// 	    if (Lim > this->LzMax)
// 	      Lim = this->LzMax;
// 	    Min = m1 + m2 - this->LzMax;
// 	    if (Min < 0)
// 	      Min = 0;
// 	    for (int m3 = Min; m3 <= Lim; ++m3)
// 	      {
// 		if (((2 * m3) > (m1 + m2)))
// 		  {
// 		    this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
// 		    this->M1Value[this->NbrInteractionFactors] = m1;
// 		    this->M2Value[this->NbrInteractionFactors] = m2;
// 		    this->M3Value[this->NbrInteractionFactors] = m3;
// 		    ++this->NbrInteractionFactors;
// 		  }
// 		++Pos;
// 	      }
// 	  }
//       this->NbrInteractionFactors = 0;

      this->NbrM12Indices = (this->NbrLzValue * (this->NbrLzValue - 1)) / 2;
      this->M1Value = new int [this->NbrM12Indices];
      this->M2Value = new int [this->NbrM12Indices];
      this->NbrM3Values = new int [this->NbrM12Indices];
      this->M3Values = new int* [this->NbrM12Indices];
      this->M3Value = 0;
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
			this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
			++this->NbrInteractionFactors;
			++TmpIndex;
		      }
		    ++Pos;
		  }
	      }
	    ++TotalIndex;
	  }
    }
  else //BOSONS
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

		TmpCoefficient[Pos] = 0.0;
		m4 = m1 + m2 - m3;

		// add the standard pseudopotential

		J = 2 * this->LzMax;
		Clebsch.InitializeCoefficientIterator(m1, m2);
		while (Clebsch.Iterate(J, ClebschCoef))
		  {
		    if (((J >> 1) & 1) != Sign)
		    {
		      TmpCoefficient[Pos] += this->PseudoPotential[J >> 1] * ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
// 		      cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << (this->LzMax - (J >> 1)) << " " << TmpCoefficient[Pos] << endl;
		    }
		  }


		// add the anisotropic pseudopotential alpha = 2
		double TmpAnisoPos = 0.0;		
		double TmpAnisoNeg = 0.0;		
		int Alpha = 2;

		int Mom = 0;
		while (Mom <= this->LzMax)
		  {
		    if (this->AnisotropicPseudoPotentialAlpha2[Mom] != 0.0)
		      {
			int L12 = this->LzMax - Mom;
			int L34 = this->LzMax - Mom - Alpha;
			int tmpmin;
			if (L12 < L34)	
				tmpmin = L12;
			else
				tmpmin = L34;

			for (int ML = -tmpmin; ML <= tmpmin; ML++)
				if ( (2*ML==(m1+m2)) && (2*ML==(m3+m4)) )
					TmpAnisoPos +=  this->AnisotropicPseudoPotentialAlpha2[Mom]  * Clebsch.GetCoefficient(m1, m2, 2*L12) * Clebsch.GetCoefficient(m3, m4, 2*L34);

			L12 = this->LzMax - Mom - Alpha;
			L34 = this->LzMax - Mom;
			if (L12 < L34)	
				tmpmin = L12;
			else
				tmpmin = L34;

			for (int ML = -tmpmin; ML <= tmpmin; ML++)
				if ( (2*ML==(m1+m2)) && (2*ML==(m3+m4)) )
					TmpAnisoNeg += this->AnisotropicPseudoPotentialAlpha2[Mom] * Clebsch.GetCoefficient(m1, m2, 2*L12) * Clebsch.GetCoefficient(m3, m4, 2*L34);
		     }
		    Mom = Mom + 2;
		  }

		if (this->Chirality == 1)
		    TmpCoefficient[Pos] += TmpAnisoPos;
		else if (this->Chirality == -1)
		    TmpCoefficient[Pos] += TmpAnisoNeg;
                else                      
	    	    TmpCoefficient[Pos] += 0.5 * (TmpAnisoPos + TmpAnisoNeg);


		// add the anisotropic pseudopotential alpha = 4
		TmpAnisoPos = 0.0;		
		TmpAnisoNeg = 0.0;		
		Alpha = 4;

		Mom = 0;
		while (Mom <= this->LzMax)
		  {
		    if (this->AnisotropicPseudoPotentialAlpha4[Mom] != 0.0)
		      {
			int L12 = this->LzMax - Mom;
			int L34 = this->LzMax - Mom - Alpha;
			int tmpmin;
			if (L12 < L34)	
				tmpmin = L12;
			else
				tmpmin = L34;

			for (int ML = -tmpmin; ML <= tmpmin; ML++)
				if ( (2*ML==(m1+m2)) && (2*ML==(m3+m4)) )
					TmpAnisoPos +=  this->AnisotropicPseudoPotentialAlpha4[Mom]  * Clebsch.GetCoefficient(m1, m2, 2*L12) * Clebsch.GetCoefficient(m3, m4, 2*L34);

			L12 = this->LzMax - Mom - Alpha;
			L34 = this->LzMax - Mom;
			if (L12 < L34)	
				tmpmin = L12;
			else
				tmpmin = L34;

			for (int ML = -tmpmin; ML <= tmpmin; ML++)
				if ( (2*ML==(m1+m2)) && (2*ML==(m3+m4)) )
					TmpAnisoNeg += this->AnisotropicPseudoPotentialAlpha4[Mom] * Clebsch.GetCoefficient(m1, m2, 2*L12) * Clebsch.GetCoefficient(m3, m4, 2*L34);
		     }
		    Mom = Mom + 2;
		  }

		if (this->Chirality == 1)
		    TmpCoefficient[Pos] += TmpAnisoPos;
		else if (this->Chirality == -1)
		    TmpCoefficient[Pos] += TmpAnisoNeg;
                else                      
	    	    TmpCoefficient[Pos] += 0.5 * (TmpAnisoPos + TmpAnisoNeg);

		//cout << m1 << " " << m2 << " " << m3 << " " << m4 << " " << TmpCoefficient[Pos] << endl;

		if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		  MaxCoefficient = TmpCoefficient[Pos];
		++Pos;



	      }
	  }
      this->NbrInteractionFactors = 0;
      this->InteractionFactors = new double [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      double Factor = 4.0;

//       this->NbrM12Indices = 0;
//       this->M1Value = new int [Pos];
//       this->M2Value = new int [Pos];
//       this->M3Value = new int [Pos];
//       for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
// 	{
// 	  for (int m2 = 0; m2 < m1; ++m2)
// 	    {
// 	      Lim = m1 + m2;
// 	      if (Lim > this->LzMax)
// 		Lim = this->LzMax;
// 	      Min = m1 + m2 - this->LzMax;
// 	      if (Min < 0)
// 		Min = 0;
// 	      for (int m3 = Min; m3 <= Lim; ++m3)
// 		{
// 		  if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
// 		    {
// 		      if ((2 * m3) > (m1 + m2))
// 			{
// 			  this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
// 			  this->M1Value[this->NbrInteractionFactors] = m1;
// 			  this->M2Value[this->NbrInteractionFactors] = m2;
// 			  this->M3Value[this->NbrInteractionFactors] = m3;
// 			  ++this->NbrInteractionFactors;
// 			}
// 		      else
// 			if ((2 * m3) == (m1 + m2))
// 			  {
// 			    this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * Factor * TmpCoefficient[Pos];
// 			    this->M1Value[this->NbrInteractionFactors] = m1;
// 			    this->M2Value[this->NbrInteractionFactors] = m2;
// 			    this->M3Value[this->NbrInteractionFactors] = m3;
// 			  ++this->NbrInteractionFactors;
// 			  }
// 		    }
// 		  ++Pos;
// 		}
// 	    }	
// 	  Lim = 2 * m1;
// 	  if (Lim > this->LzMax)
// 	    Lim = this->LzMax;
// 	  Min = 2 * m1 - this->LzMax;
// 	  if (Min < 0)
// 	    Min = 0;
// 	  for (int m3 = Min; m3 <= Lim; ++m3)
// 	    {
// 	      if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
// 		{
// 		  if (m3 > m1)
// 		    {
// 		      this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * Factor * TmpCoefficient[Pos];
// 		      this->M1Value[this->NbrInteractionFactors] = m1;
// 		      this->M2Value[this->NbrInteractionFactors] = m1;
// 		      this->M3Value[this->NbrInteractionFactors] = m3;
// 		      ++this->NbrInteractionFactors;
// 		    }
// 		  else
// 		    if (m3 == m1)
// 		      {
// 			this->InteractionFactors[this->NbrInteractionFactors] = 0.25 * Factor * TmpCoefficient[Pos];
// 			this->M1Value[this->NbrInteractionFactors] = m1;
// 			this->M2Value[this->NbrInteractionFactors] = m1;
// 			this->M3Value[this->NbrInteractionFactors] = m3;
// 			++this->NbrInteractionFactors;
// 		      }
// 		}
// 	      ++Pos;
// 	    }
// 	}

      this->NbrM12Indices = (this->NbrLzValue * (this->NbrLzValue + 1)) / 2;
      this->M1Value = new int [this->NbrM12Indices];
      this->M2Value = new int [this->NbrM12Indices];
      this->NbrM3Values = new int [this->NbrM12Indices];
      this->M3Value = 0;
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
			  this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
			  ++this->NbrInteractionFactors;
			  ++TmpIndex;
			}
		      else
			if ((2 * m3) == (m1 + m2))
			  {
			    this->M3Values[TotalIndex][TmpIndex] = m3;
			    this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * Factor * TmpCoefficient[Pos];
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
		      this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * Factor * TmpCoefficient[Pos];
		      ++this->NbrInteractionFactors;
		      ++TmpIndex;
		    }
		  else
		    if (m3 == m1)
		      {
			this->M3Values[TotalIndex][TmpIndex] = m3;
			this->InteractionFactors[this->NbrInteractionFactors] = 0.25 * Factor * TmpCoefficient[Pos];
			++this->NbrInteractionFactors;
			++TmpIndex;
		      }
		  ++Pos;
		}
	    }
	  ++TotalIndex;
	}     
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
