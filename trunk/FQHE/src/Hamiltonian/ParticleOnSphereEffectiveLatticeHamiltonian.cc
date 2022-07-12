////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//     SU(2) spin and a generic interaction defined by its pseudopotential    //
//                                                                            //
//                        last modification : 07/06/2007                      //
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


#include "Hamiltonian/ParticleOnSphereEffectiveLatticeHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include "GeneralTools/StringTools.h"

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
// alpha = parameter encoding deviation of flux-density from half filling on the lattice
// pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialUpDown =  one-body tunnelling potential (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereEffectiveLatticeHamiltonian::ParticleOnSphereEffectiveLatticeHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, double alpha, double** pseudoPotential, double* onebodyPotentialUpUp, double* onebodyPotentialDownDown, double* onebodyPotentialUpDown, AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = false;
  this->Architecture = architecture;
  this->Alpha = alpha;
  this->PseudoPotentials = new double* [4];
  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;
  if (pseudoPotential!=NULL)
    HaveOtherMixed=true;
  else
    HaveOtherMixed=false;
  for (int j = 0; j < 4; ++j)
    {
      this->PseudoPotentials[j] = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	if (pseudoPotential!=NULL)
	  this->PseudoPotentials[j][i] = pseudoPotential[j][this->LzMax - i];
	else
	  this->PseudoPotentials[j][i] = 0.0;
    }
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  if (onebodyPotentialUpUp != 0)
    {
      this->OneBodyInteractionFactorsupup = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsupup[i] = onebodyPotentialUpUp[i];
    }
  this->OneBodyInteractionFactorsdowndown = 0;
  if (onebodyPotentialDownDown != 0)
    {
      this->OneBodyInteractionFactorsdowndown = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsdowndown[i] = onebodyPotentialDownDown[i];
    }
  this->OneBodyInteractionFactorsupdown = 0;
  if (onebodyPotentialUpDown != 0)
    {
      this->OneBodyInteractionFactorsupdown = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsupdown[i] = onebodyPotentialUpDown[i];
    }
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout  << "fast = ";
	  PrintMemorySize(cout, TmpMemory)<<endl;
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

}

// destructor
//

ParticleOnSphereEffectiveLatticeHamiltonian::~ParticleOnSphereEffectiveLatticeHamiltonian() 
{
  for (int j = 0; j < 4; ++j)
    delete[] this->PseudoPotentials[j];
  delete[] this->PseudoPotentials;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereEffectiveLatticeHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereEffectiveLatticeHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereEffectiveLatticeHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnSphereEffectiveLatticeHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereEffectiveLatticeHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
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

Complex ParticleOnSphereEffectiveLatticeHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> ParticleOnSphereEffectiveLatticeHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> ParticleOnSphereEffectiveLatticeHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// evaluate all interaction factors
//   

void ParticleOnSphereEffectiveLatticeHamiltonian::EvaluateInteractionFactors()
{
  // int Lim;
  // int Min;
  // int Pos = 0;
  ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
  // int m4;
  double ClebschCoef;
  long TotalNbrInteractionFactors = 0;
  int J;
  
  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  double TmpCoefficient = 0.0;

  //INTER
  this->NbrInterSectorSums = 2 * this->LzMax + 1;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      ++this->NbrInterSectorIndicesPerSum[m1 + m2];      
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
      this->NbrInterSectorIndicesPerSum[i] = 0;
    }
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      {
	this->InterSectorIndicesPerSum[(m1 + m2)][this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
	this->InterSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
	++this->NbrInterSectorIndicesPerSum[(m1 + m2)];
      }

  //FERMIONS  *** this code is not really used for lattices, as we first considered bosons
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      //INTRA
      this->NbrIntraSectorSums = 2 * this->LzMax - 1;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    this->IntraSectorIndicesPerSum[(m1 + m2) - 1][this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1] = m1;
	    this->IntraSectorIndicesPerSum[(m1 + m2) - 1][1 + (this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
	  }

      this->InteractionFactorsupup = new double* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new double* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndown[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->IntraSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
	      int m2 = (this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->IntraSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		  int m4 = (this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
		  Clebsch.InitializeCoefficientIterator(m1, m2);
		  this->InteractionFactorsupup[i][Index] = 0.0;
		  this->InteractionFactorsdowndown[i][Index] = 0.0;
		  while (Clebsch.Iterate(J, ClebschCoef))
		    {
		      if (((J >> 1) & 1) == Sign)
			{
			  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
			  this->InteractionFactorsupup[i][Index] += this->PseudoPotentials[0][J >> 1] * TmpCoefficient;
			  this->InteractionFactorsdowndown[i][Index] += this->PseudoPotentials[1][J >> 1] * TmpCoefficient;
			}
		    }
		  this->InteractionFactorsupup[i][Index] *= -4.0;
		  this->InteractionFactorsdowndown[i][Index] *= -4.0;
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      //FILL IN INTER

      this->InteractionFactorsupdown = new double* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      double Factor = 2.0;
	      int m1 = (this->InterSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
	      int m2 = (this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->InterSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		  int m4 = (this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
		  Clebsch.InitializeCoefficientIterator(m1, m2);
		  this->InteractionFactorsupdown[i][Index] = 0.0;
		  while (Clebsch.Iterate(J, ClebschCoef))
		    {
		      TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		      this->InteractionFactorsupdown[i][Index] += this->PseudoPotentials[2][J >> 1] * TmpCoefficient;
		    }
		  this->InteractionFactorsupdown[i][Index] *= -Factor;
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

      //************************************************************************************************************
      //********************************     M I X E D     T E R M S ***********************************************
      //************************************************************************************************************
      // Mixed term indices that have the same structure as intra terms

      this->NbrMixedIntraSectorSums = 2 * this->LzMax - 1;
      this->NbrMixedIntraSectorIndicesPerSum = new int[this->NbrMixedIntraSectorSums];
      for (int i = 0; i < this->NbrMixedIntraSectorSums; ++i)
	this->NbrMixedIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrMixedIntraSectorIndicesPerSum[(m1 + m2) - 1];
      this->MixedIntraSectorIndicesPerSum = new int* [this->NbrMixedIntraSectorSums];
      for (int i = 0; i < this->NbrMixedIntraSectorSums; ++i)
	{
	  this->MixedIntraSectorIndicesPerSum[i] = new int[2 * this->NbrMixedIntraSectorIndicesPerSum[i]];      
	  this->NbrMixedIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    this->MixedIntraSectorIndicesPerSum[(m1 + m2) - 1][this->NbrMixedIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1] = m1;
	    this->MixedIntraSectorIndicesPerSum[(m1 + m2) - 1][1 + (this->NbrMixedIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1)] = m2;
	    ++this->NbrMixedIntraSectorIndicesPerSum[(m1 + m2) - 1];
	  }

      //Fill in the matrix elements

      this->InteractionFactorsmixedintra= new double* [this->NbrMixedIntraSectorSums];
      for (int i = 0; i < this->NbrMixedIntraSectorSums; ++i)
	{
	  this->InteractionFactorsmixedintra[i] = new double[this->NbrMixedIntraSectorIndicesPerSum[i] * this->NbrMixedIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrMixedIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->MixedIntraSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
	      int m2 = (this->MixedIntraSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
	      for (int j2 = 0; j2 < this->NbrMixedIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->MixedIntraSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		  int m4 = (this->MixedIntraSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
		  Clebsch.InitializeCoefficientIterator(m1, m2);
		  this->InteractionFactorsmixedintra[i][Index] = 0.0;
		  while (Clebsch.Iterate(J, ClebschCoef))
		    {
		      if (((J >> 1) & 1) == Sign)
			{
			  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
			  this->InteractionFactorsmixedintra[i][Index] += this->PseudoPotentials[3][J >> 1] * TmpCoefficient;
			}
		    }
		  this->InteractionFactorsmixedintra[i][Index] *= -4.0;
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
        }


      //Mixed term indices that have the same structure as inter terms


      this->NbrMixedInterSectorSums = 2 * this->LzMax + 1;
      this->NbrMixedInterSectorIndicesPerSum = new int[this->NbrMixedInterSectorSums];
      for (int i = 0; i < this->NbrMixedInterSectorSums; ++i)
	this->NbrMixedInterSectorIndicesPerSum[i] = 0;
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = 0; m2 <= this->LzMax; ++m2)
	  ++this->NbrMixedInterSectorIndicesPerSum[m1 + m2];      
      this->MixedInterSectorIndicesPerSum = new int* [this->NbrMixedInterSectorSums];
      for (int i = 0; i < this->NbrMixedInterSectorSums; ++i)
	{
	  this->MixedInterSectorIndicesPerSum[i] = new int[2 * this->NbrMixedInterSectorIndicesPerSum[i]];      
	  this->NbrMixedInterSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = 0; m2 <= this->LzMax; ++m2)
	  {
	    this->MixedInterSectorIndicesPerSum[(m1 + m2)][this->NbrMixedInterSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
	    this->MixedInterSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrMixedInterSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
	    ++this->NbrMixedInterSectorIndicesPerSum[(m1 + m2)];
	  }


      //Fill in the matrix elements


      this->InteractionFactorsmixedinter = new double* [this->NbrMixedInterSectorSums];
      for (int i = 0; i < this->NbrMixedInterSectorSums; ++i)
	{
	  this->InteractionFactorsmixedinter[i] = new double[this->NbrMixedInterSectorIndicesPerSum[i] * this->NbrMixedInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrMixedInterSectorIndicesPerSum[i]; ++j1)
	    {
	      double Factor = 2.0;
	      int m1 = (this->MixedInterSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
	      int m2 = (this->MixedInterSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
	      for (int j2 = 0; j2 < this->NbrMixedInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->MixedInterSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		  int m4 = (this->MixedInterSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
		  Clebsch.InitializeCoefficientIterator(m1, m2);
		  this->InteractionFactorsmixedinter[i][Index] = 0.0;
		  while (Clebsch.Iterate(J, ClebschCoef))
		    {
		      TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		      this->InteractionFactorsmixedinter[i][Index] += this->PseudoPotentials[3][J >> 1] * TmpCoefficient;
		    }
		  //The sign of the Factor is opposite to UpDown case
		  this->InteractionFactorsmixedinter[i][Index] *= Factor;
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}




      //*********************************************************************************************************
      //********************************     M I X E D     T E R M S ***********************************************
      //*********************************************************************************************************

    }
  else // bosons
    {
      //INTRA
      this->NbrIntraSectorSums = 2 * this->LzMax+1;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2)];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  {
	    this->IntraSectorIndicesPerSum[m1 + m2][this->NbrIntraSectorIndicesPerSum[m1 + m2] << 1] = m1;
	    this->IntraSectorIndicesPerSum[m1 + m2][1 + (this->NbrIntraSectorIndicesPerSum[m1 + m2] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[m1 + m2];
	  }

      this->InteractionFactorsupup = new double* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new double* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndown[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->IntraSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
	      int m2 = (this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->IntraSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		  int m4 = (this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
		  Clebsch.InitializeCoefficientIterator(m1, m2);
		  this->InteractionFactorsupup[i][Index] = 0.0;
		  this->InteractionFactorsdowndown[i][Index] = 0.0;
		  while (Clebsch.Iterate(J, ClebschCoef))
		    {
		      if (((J >> 1) & 1) != Sign)
			{
			  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
			  this->InteractionFactorsupup[i][Index] += this->PseudoPotentials[0][J >> 1] * TmpCoefficient;
			  this->InteractionFactorsdowndown[i][Index] += this->PseudoPotentials[1][J >> 1] * TmpCoefficient;
			  if (J==2*LzMax) // Delta interactions on lattice
			    {
			      this->InteractionFactorsupup[i][Index] += 1.0 * TmpCoefficient;
			      this->InteractionFactorsdowndown[i][Index] += 1.0 * TmpCoefficient;
			    }
			}
		    }
		  if (m1 != m2)
		    {
		      this->InteractionFactorsupup[i][Index] *= 2.0;
		      this->InteractionFactorsdowndown[i][Index] *= 2.0;
		    }
		  if (m3 != m4)
		    {
		      this->InteractionFactorsupup[i][Index] *= 2.0;
		      this->InteractionFactorsdowndown[i][Index] *= 2.0;
		    }
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      //FILL IN INTER

      this->InteractionFactorsupdown = new double* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      double Factor = 2.0;
	      int m1 = (this->InterSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
	      int m2 = (this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->InterSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		  int m4 = (this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
		  Clebsch.InitializeCoefficientIterator(m1, m2);
		  this->InteractionFactorsupdown[i][Index] = 0.0;
		  while (Clebsch.Iterate(J, ClebschCoef))
		    {
		      TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		      this->InteractionFactorsupdown[i][Index] += this->PseudoPotentials[2][J >> 1] * TmpCoefficient;
 		      if (J==2*LzMax) // Delta interactions on lattice
 			this->InteractionFactorsupdown[i][Index] += 1.0*TmpCoefficient;
		    }
		  this->InteractionFactorsupdown[i][Index] *= Factor;
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
      
      //************************************************************************************************************ /
      //********************************     M I X E D     T E R M S *********************************************** /
      //************************************************************************************************************ /
      // Mixed term indices that have the same structure as intra terms
	      
      this->NbrMixedIntraSectorSums = 2 * this->LzMax + 1;
      this->NbrMixedIntraSectorIndicesPerSum = new int[this->NbrMixedIntraSectorSums];
      for (int i = 0; i < this->NbrMixedIntraSectorSums; ++i)
	this->NbrMixedIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  ++this->NbrMixedIntraSectorIndicesPerSum[m1 + m2];
      this->MixedIntraSectorIndicesPerSum = new int* [this->NbrMixedIntraSectorSums];
      for (int i = 0; i < this->NbrMixedIntraSectorSums; ++i)
	{
	  this->MixedIntraSectorIndicesPerSum[i] = new int[2 * this->NbrMixedIntraSectorIndicesPerSum[i]];      
	  this->NbrMixedIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  {
	    this->MixedIntraSectorIndicesPerSum[m1 + m2][this->NbrMixedIntraSectorIndicesPerSum[m1 + m2] << 1] = m1;
	    this->MixedIntraSectorIndicesPerSum[m1 + m2][1 + (this->NbrMixedIntraSectorIndicesPerSum[m1 + m2] << 1)] = m2;
	    ++this->NbrMixedIntraSectorIndicesPerSum[m1 + m2];
	  }

      //Fill in the matrix elements

      this->InteractionFactorsmixedintra= new double* [this->NbrMixedIntraSectorSums];
      for (int i = 0; i < this->NbrMixedIntraSectorSums; ++i)
	{
	  this->InteractionFactorsmixedintra[i] = new double[this->NbrMixedIntraSectorIndicesPerSum[i] * this->NbrMixedIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrMixedIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = (this->MixedIntraSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
	      int m2 = (this->MixedIntraSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
	      for (int j2 = 0; j2 < this->NbrMixedIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = (this->MixedIntraSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		  int m4 = (this->MixedIntraSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
		  Clebsch.InitializeCoefficientIterator(m1, m2);
		  this->InteractionFactorsmixedintra[i][Index] = 0.0;
		  while (Clebsch.Iterate(J, ClebschCoef))
		    {
		      if (((J >> 1) & 1) != Sign)
			{
			  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
			  this->InteractionFactorsmixedintra[i][Index] += this->PseudoPotentials[3][J >> 1] * TmpCoefficient;
 			  if (J==2*LzMax) // anomalous terms arising from effective lattice model in V0 channel
 			    this->InteractionFactorsmixedintra[i][Index] += -2.0*M_PI*this->Alpha * TmpCoefficient;
			}
		    }
		  if (m1 != m2)
		    this->InteractionFactorsmixedintra[i][Index] *= 2.0;
		  if (m3 != m4)
		    this->InteractionFactorsmixedintra[i][Index] *= 2.0;
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
        }
      
      if (this->HaveOtherMixed)
	{
	  //Mixed term indices that have the same structure as inter terms
	  this->NbrMixedInterSectorSums = 2 * this->LzMax + 1;
	  this->NbrMixedInterSectorIndicesPerSum = new int[this->NbrMixedInterSectorSums];
	  for (int i = 0; i < this->NbrMixedInterSectorSums; ++i)
	    this->NbrMixedInterSectorIndicesPerSum[i] = 0;
	  for (int m1 = 0; m1 <= this->LzMax; ++m1)
	    for (int m2 = 0; m2 <= this->LzMax; ++m2)
	      ++this->NbrMixedInterSectorIndicesPerSum[m1 + m2];      
	  this->MixedInterSectorIndicesPerSum = new int* [this->NbrMixedInterSectorSums];
	  for (int i = 0; i < this->NbrMixedInterSectorSums; ++i)
	    {
	      this->MixedInterSectorIndicesPerSum[i] = new int[2 * this->NbrMixedInterSectorIndicesPerSum[i]];      
	      this->NbrMixedInterSectorIndicesPerSum[i] = 0;
	    }
	  for (int m1 = 0; m1 <= this->LzMax; ++m1)
	    for (int m2 = 0; m2 <= this->LzMax; ++m2)
	      {
		this->MixedInterSectorIndicesPerSum[(m1 + m2)][this->NbrMixedInterSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
		this->MixedInterSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrMixedInterSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
		++this->NbrMixedInterSectorIndicesPerSum[(m1 + m2)];
	      }


	  //Fill in the matrix elements


	  this->InteractionFactorsmixedinter = new double* [this->NbrMixedInterSectorSums];
	  for (int i = 0; i < this->NbrMixedInterSectorSums; ++i)
	    {
	      this->InteractionFactorsmixedinter[i] = new double[this->NbrMixedInterSectorIndicesPerSum[i] * this->NbrMixedInterSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrMixedInterSectorIndicesPerSum[i]; ++j1)
		{
		  double Factor = 2.0;
		  int m1 = (this->MixedInterSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
		  int m2 = (this->MixedInterSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
		  for (int j2 = 0; j2 < this->NbrMixedInterSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = (this->MixedInterSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
		      int m4 = (this->MixedInterSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
		      Clebsch.InitializeCoefficientIterator(m1, m2);
		      this->InteractionFactorsmixedinter[i][Index] = 0.0;
		      while (Clebsch.Iterate(J, ClebschCoef))
			{
			  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
			  this->InteractionFactorsmixedinter[i][Index] += this->PseudoPotentials[3][J >> 1] * TmpCoefficient;
			}
		      //For bosons, the sign of the Factor is probably same as UpDown case
		      this->InteractionFactorsmixedinter[i][Index] *= Factor;
		      ++TotalNbrInteractionFactors;
		      ++Index;
		    }
		}
	    }
	}
      else
	{
	  this->NbrMixedInterSectorSums = 0;
	}
    }

  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

