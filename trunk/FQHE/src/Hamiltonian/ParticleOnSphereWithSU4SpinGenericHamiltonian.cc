////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//     SU(4) spin and a generic interaction defined by its pseudopotential    //
//                                                                            //
//                        last modification : 28/11/2006                      //
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


#include "Hamiltonian/ParticleOnSphereWithSU4SpinGenericHamiltonian.h"
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
// pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin/isospin sector (sorted as up-up, up-um, up-dp, up-dm, um-um, um-dp, um-dm, dp-dp, dp-dm and dm-dm)
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereWithSU4SpinGenericHamiltonian::ParticleOnSphereWithSU4SpinGenericHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int lzmax, 
											     double** pseudoPotential, double zeeman,
											     AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
											     char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = false;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Pseudopotentials = new double* [10];
  for (int j = 0; j < 10; ++j)
    {
      this->Pseudopotentials[j] = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->Pseudopotentials[j][i] = pseudoPotential[j][this->LzMax - i];
    }

  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();

  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsumum = 0;
  this->OneBodyInteractionFactorsdpdp = 0;
  this->OneBodyInteractionFactorsdmdm = 0;
 
  if (zeeman != 0.0)
    {
      this->OneBodyTermFlag = true;
      this->OneBodyInteractionFactorsupup = new double [this->NbrLzValue];
      this->OneBodyInteractionFactorsumum = new double [this->NbrLzValue];
      this->OneBodyInteractionFactorsdpdp = new double [this->NbrLzValue];
      this->OneBodyInteractionFactorsdmdm = new double [this->NbrLzValue];
      for (int j = 0; j <= this->LzMax; ++j)
	{ 
	  this->OneBodyInteractionFactorsupup[j] = - zeeman;
	  this->OneBodyInteractionFactorsumum[j] = - zeeman;
	  this->OneBodyInteractionFactorsdpdp[j] = zeeman;
	  this->OneBodyInteractionFactorsdmdm[j] = zeeman;	
	}
    }

  this->EvaluateInteractionFactors();
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;

  this->HamiltonianShift = 0.0;
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

}

// destructor
//

ParticleOnSphereWithSU4SpinGenericHamiltonian::~ParticleOnSphereWithSU4SpinGenericHamiltonian() 
{
  for (int j = 0; j < 10; ++j)
    delete[] this->Pseudopotentials[j];
  delete[] this->Pseudopotentials;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSU4SpinGenericHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// evaluate all interaction factors
//   

void ParticleOnSphereWithSU4SpinGenericHamiltonian::EvaluateInteractionFactors()
{
  this->M1IntraValue = 0;
  this->M1InterValue = 0;
  
  long TotalNbrInteractionFactors = 0;

  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);

  this->NbrInterSectorSums = (2 * this->NbrLzValue) + 1;
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
	int TmpIndex = m1 + m2;
	this->InterSectorIndicesPerSum[TmpIndex][this->NbrInterSectorIndicesPerSum[TmpIndex] << 1] = m1;
	this->InterSectorIndicesPerSum[TmpIndex][1 + (this->NbrInterSectorIndicesPerSum[TmpIndex] << 1)] = m2;
	++this->NbrInterSectorIndicesPerSum[TmpIndex];
      }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NbrIntraSectorSums = 2 * this->NbrLzValue + 1;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[m1 + m2];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    int TmpIndex = m1 + m2;
	    this->IntraSectorIndicesPerSum[TmpIndex][this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1] = m1;
	    this->IntraSectorIndicesPerSum[TmpIndex][1 + (this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[TmpIndex];
	  }

      this->InteractionFactorsupup = new double* [this->NbrIntraSectorSums];
      this->InteractionFactorsumum = new double* [this->NbrIntraSectorSums];
      this->InteractionFactorsdpdp = new double* [this->NbrIntraSectorSums];
      this->InteractionFactorsdmdm = new double* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsumum[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdpdp[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdmdm[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];

		  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[0])
							    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, Clebsch, this->Pseudopotentials[0])
							    - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[0])
							    - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[0]));
		  this->InteractionFactorsumum[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[4])
							    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, Clebsch, this->Pseudopotentials[4])
							    - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[4])
							    - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[4]));
		  this->InteractionFactorsdpdp[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[7])
							    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, Clebsch, this->Pseudopotentials[7])
							    - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[7])
							    - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[7]));
		  this->InteractionFactorsdmdm[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[9])
							    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, Clebsch, this->Pseudopotentials[9])
							    - this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[9])
							    - this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[9]));
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      this->InteractionFactorsupum = new double* [this->NbrInterSectorSums];
      this->InteractionFactorsupdp = new double* [this->NbrInterSectorSums];
      this->InteractionFactorsupdm = new double* [this->NbrInterSectorSums];
      this->InteractionFactorsumdp = new double* [this->NbrInterSectorSums];
      this->InteractionFactorsumdm = new double* [this->NbrInterSectorSums];
      this->InteractionFactorsdpdm = new double* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupum[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdp[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdm[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsumdp[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsumdm[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsdpdm[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  this->InteractionFactorsupum[i][Index] = (-this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[1])
							    -this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[1]));
		  this->InteractionFactorsupdp[i][Index] = (-this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[2])
							    -this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[2]));
		  this->InteractionFactorsupdm[i][Index] = (-this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[3])
							    -this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[3]));
		  this->InteractionFactorsumdp[i][Index] = (-this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[5])
							    -this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[5]));
		  this->InteractionFactorsumdm[i][Index] = (-this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[6])
							    -this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[6]));
		  this->InteractionFactorsdpdm[i][Index] = (-this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[8])
							    -this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[8]));
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      this->NbrIntraSectorSums = 2 * this->NbrLzValue + 1;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  ++this->NbrIntraSectorIndicesPerSum[m1 + m2];
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i] != 0)
	    this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  else
	    this->IntraSectorIndicesPerSum[i] = 0;
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  {
	    int TmpIndex = (m1 + m2);
	    this->IntraSectorIndicesPerSum[TmpIndex][this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1] = m1;
	    this->IntraSectorIndicesPerSum[TmpIndex][1 + (this->NbrIntraSectorIndicesPerSum[TmpIndex] << 1)] = m2;
	    ++this->NbrIntraSectorIndicesPerSum[TmpIndex];
	  }

      this->InteractionFactorsupup = new double* [this->NbrIntraSectorSums];
      this->InteractionFactorsumum = new double* [this->NbrIntraSectorSums];
      this->InteractionFactorsdpdp = new double* [this->NbrIntraSectorSums];
      this->InteractionFactorsdmdm = new double* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsumum[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdpdp[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdmdm[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];


		  if (m1 != m2)
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[0])
								    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, Clebsch, this->Pseudopotentials[0])
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[0])
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[0]));
			  this->InteractionFactorsumum[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[4])
								    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, Clebsch, this->Pseudopotentials[4])
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[4])
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[4]));
			  this->InteractionFactorsdpdp[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[7])
								    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, Clebsch, this->Pseudopotentials[7])
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[7])
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[7]));
			  this->InteractionFactorsdmdm[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[9])
								    + this->EvaluateInteractionCoefficient(m2, m1, m4, m3, Clebsch, this->Pseudopotentials[9])
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[9])
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[9]));
			}
		      else
			{
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[0])
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[0]));
			  this->InteractionFactorsumum[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[4])
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[4]));
			  this->InteractionFactorsdpdp[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[7])
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[7]));
			  this->InteractionFactorsdmdm[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[9])
								    + this->EvaluateInteractionCoefficient(m2, m1, m3, m4, Clebsch, this->Pseudopotentials[9]));
			}
		    }
		  else
		    {
		      if (m3 != m4)
			{
			  this->InteractionFactorsupup[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[0])
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[0]));
			  this->InteractionFactorsumum[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[4])
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[4]));
			  this->InteractionFactorsdpdp[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[7])
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[7]));
			  this->InteractionFactorsdmdm[i][Index] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[9])
								    + this->EvaluateInteractionCoefficient(m1, m2, m4, m3, Clebsch, this->Pseudopotentials[9]));
			}
		      else
			{
			  this->InteractionFactorsupup[i][Index] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[0]);
			  this->InteractionFactorsumum[i][Index] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[4]);
			  this->InteractionFactorsdpdp[i][Index] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[7]);
			  this->InteractionFactorsdmdm[i][Index] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[9]);
			}
		    }
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      this->InteractionFactorsupum = new double* [this->NbrInterSectorSums];
      this->InteractionFactorsupdp = new double* [this->NbrInterSectorSums];
      this->InteractionFactorsupdm = new double* [this->NbrInterSectorSums];
      this->InteractionFactorsumdp = new double* [this->NbrInterSectorSums];
      this->InteractionFactorsumdm = new double* [this->NbrInterSectorSums];
      this->InteractionFactorsdpdm = new double* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupum[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdp[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdm[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsumdp[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsumdm[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsdpdm[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      double Factor = 2.0;
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];

		  this->InteractionFactorsupum[i][Index] = Factor * this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[1]);
		  this->InteractionFactorsupdp[i][Index] = Factor * this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[2]);
		  this->InteractionFactorsupdm[i][Index] = Factor * this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[3]);
		  this->InteractionFactorsumdp[i][Index] = Factor * this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[5]);
		  this->InteractionFactorsumdm[i][Index] = Factor * this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[6]);
		  this->InteractionFactorsdpdm[i][Index] = Factor * this->EvaluateInteractionCoefficient(m1, m2, m3, m4, Clebsch, this->Pseudopotentials[8]);
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }

  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
//   int Lim;
//   int Min;
//   int Pos = 0;
//   ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
//   int J = 2 * this->LzMax - 2;
//   int m4;
//   double ClebschCoef;
//   long TotalNbrInteractionFactors = 0;

//   int Sign = 1;
//   if (this->LzMax & 1)
//     Sign = 0;
//   double TmpCoefficient = 0.0;

//   this->NbrInterSectorSums = 2 * this->LzMax + 1;
//   this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
//   for (int i = 0; i < this->NbrInterSectorSums; ++i)
//     this->NbrInterSectorIndicesPerSum[i] = 0;
//   for (int m1 = 0; m1 <= this->LzMax; ++m1)
//     for (int m2 = 0; m2 <= this->LzMax; ++m2)
//       ++this->NbrInterSectorIndicesPerSum[m1 + m2];      
//   this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
//   for (int i = 0; i < this->NbrInterSectorSums; ++i)
//     {
//       this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
//       this->NbrInterSectorIndicesPerSum[i] = 0;
//     }
//   for (int m1 = 0; m1 <= this->LzMax; ++m1)
//     for (int m2 = 0; m2 <= this->LzMax; ++m2)
//       {
// 	this->InterSectorIndicesPerSum[(m1 + m2)][this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
// 	this->InterSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
// 	++this->NbrInterSectorIndicesPerSum[(m1 + m2)];
//       }

//   if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
//     {
//       this->NbrIntraSectorSums = 2 * this->LzMax - 1;
//       this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
//       for (int i = 0; i < this->NbrIntraSectorSums; ++i)
// 	this->NbrIntraSectorIndicesPerSum[i] = 0;      
//       for (int m1 = 0; m1 < this->LzMax; ++m1)
// 	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
// 	  ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
//       this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
//       for (int i = 0; i < this->NbrIntraSectorSums; ++i)
// 	{
// 	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
// 	  this->NbrIntraSectorIndicesPerSum[i] = 0;
// 	}
//       for (int m1 = 0; m1 < this->LzMax; ++m1)
// 	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
// 	  {
// 	    this->IntraSectorIndicesPerSum[(m1 + m2) - 1][this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1] = m1;
// 	    this->IntraSectorIndicesPerSum[(m1 + m2) - 1][1 + (this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1)] = m2;
// 	    ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
// 	  }

//       this->InteractionFactorsupup = new double* [this->NbrIntraSectorSums];
//       this->InteractionFactorsumum = new double* [this->NbrIntraSectorSums];
//       this->InteractionFactorsdpdp = new double* [this->NbrIntraSectorSums];
//       this->InteractionFactorsdmdm = new double* [this->NbrIntraSectorSums];
//       for (int i = 0; i < this->NbrIntraSectorSums; ++i)
// 	{
// 	  this->InteractionFactorsupup[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
// 	  this->InteractionFactorsumum[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
// 	  this->InteractionFactorsdpdp[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
// 	  this->InteractionFactorsdmdm[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
// 	  int Index = 0;
// 	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
// 	    {
// 	      int m1 = (this->IntraSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
// 	      int m2 = (this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
// 	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
// 		{
// 		  int m3 = (this->IntraSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
// 		  int m4 = (this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
// 		  Clebsch.InitializeCoefficientIterator(m1, m2);
// 		  this->InteractionFactorsupup[i][Index] = 0.0;
// 		  this->InteractionFactorsumum[i][Index] = 0.0;
// 		  this->InteractionFactorsdpdp[i][Index] = 0.0;
// 		  this->InteractionFactorsdmdm[i][Index] = 0.0;
// 		  while (Clebsch.Iterate(J, ClebschCoef))
// 		    {
// 		      if (((J >> 1) & 1) == Sign)
// 			{
// 			  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
// 			  this->InteractionFactorsupup[i][Index] += this->Pseudopotentials[0][J >> 1] * TmpCoefficient;
// 			  this->InteractionFactorsumum[i][Index] += this->Pseudopotentials[4][J >> 1] * TmpCoefficient;
// 			  this->InteractionFactorsdpdp[i][Index] += this->Pseudopotentials[7][J >> 1] * TmpCoefficient;
// 			  this->InteractionFactorsdmdm[i][Index] += this->Pseudopotentials[9][J >> 1] * TmpCoefficient;
// 			}
// 		    }
// 		  this->InteractionFactorsupup[i][Index] *= -4.0;
// 		  this->InteractionFactorsumum[i][Index] *= -4.0;
// 		  this->InteractionFactorsdpdp[i][Index] *= -4.0;
// 		  this->InteractionFactorsdmdm[i][Index] *= -4.0;
// 		  TotalNbrInteractionFactors += 4;
// 		  ++Index;
// 		}
// 	    }
// 	}

//       this->InteractionFactorsupum = new double* [this->NbrInterSectorSums];
//       this->InteractionFactorsupdp = new double* [this->NbrInterSectorSums];
//       this->InteractionFactorsupdm = new double* [this->NbrInterSectorSums];
//       this->InteractionFactorsumdp = new double* [this->NbrInterSectorSums];
//       this->InteractionFactorsumdm = new double* [this->NbrInterSectorSums];
//       this->InteractionFactorsdpdm = new double* [this->NbrInterSectorSums];
//       for (int i = 0; i < this->NbrInterSectorSums; ++i)
// 	{
// 	  this->InteractionFactorsupum[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
// 	  this->InteractionFactorsupdp[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
// 	  this->InteractionFactorsupdm[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
// 	  this->InteractionFactorsumdp[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
// 	  this->InteractionFactorsumdm[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
// 	  this->InteractionFactorsdpdm[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
// 	  int Index = 0;
// 	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
// 	    {
// 	      double Factor = 2.0;
// 	      int m1 = (this->InterSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
// 	      int m2 = (this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
// 	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
// 		{
// 		  int m3 = (this->InterSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
// 		  int m4 = (this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
// 		  Clebsch.InitializeCoefficientIterator(m1, m2);
// 		  this->InteractionFactorsupum[i][Index] = 0.0;
// 		  this->InteractionFactorsupdp[i][Index] = 0.0;
// 		  this->InteractionFactorsupdm[i][Index] = 0.0;
// 		  this->InteractionFactorsumdp[i][Index] = 0.0;
// 		  this->InteractionFactorsumdm[i][Index] = 0.0;
// 		  this->InteractionFactorsdpdm[i][Index] = 0.0;
// 		  while (Clebsch.Iterate(J, ClebschCoef))
// 		    {
// 		      TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
// 		      this->InteractionFactorsupum[i][Index] += this->Pseudopotentials[1][J >> 1] * TmpCoefficient;
// 		      this->InteractionFactorsupdp[i][Index] += this->Pseudopotentials[2][J >> 1] * TmpCoefficient;
// 		      this->InteractionFactorsupdm[i][Index] += this->Pseudopotentials[3][J >> 1] * TmpCoefficient;
// 		      this->InteractionFactorsumdp[i][Index] += this->Pseudopotentials[5][J >> 1] * TmpCoefficient;
// 		      this->InteractionFactorsumdm[i][Index] += this->Pseudopotentials[6][J >> 1] * TmpCoefficient;
// 		      this->InteractionFactorsdpdm[i][Index] += this->Pseudopotentials[8][J >> 1] * TmpCoefficient;
// 		    }
// 		  this->InteractionFactorsupum[i][Index] *= -Factor;
// 		  this->InteractionFactorsupdp[i][Index] *= -Factor;
// 		  this->InteractionFactorsupdm[i][Index] *= -Factor;
// 		  this->InteractionFactorsumdp[i][Index] *= -Factor;
// 		  this->InteractionFactorsumdm[i][Index] *= -Factor;
// 		  this->InteractionFactorsdpdm[i][Index] *= -Factor;
// 		  TotalNbrInteractionFactors += 6;
// 		  ++Index;
// 		}
// 	    }
// 	}
//     }
//   else
//     {
//     }
//   cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
//   cout << "====================================" << endl;
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// clebsch = reference to the Clebsch-Gordan coefficients
// pseudopotentials = pseudopotential coefficients
// return value = numerical coefficient

double ParticleOnSphereWithSU4SpinGenericHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4, ClebschGordanCoefficients& clebsch,
										    double* pseudopotentials)
{
  double Coefficient = 0.0;
  double ClebschCoef;
  int J;
  clebsch.InitializeCoefficientIterator((m1 << 1) - this->LzMax, (m2 << 1) - this->LzMax);
  while (clebsch.Iterate(J, ClebschCoef))
    {
      Coefficient += pseudopotentials[J >> 1] * ClebschCoef * clebsch.GetCoefficient((m3 << 1) - this->LzMax, (m4 << 1) - this->LzMax, J);
    }
  return Coefficient;
}


