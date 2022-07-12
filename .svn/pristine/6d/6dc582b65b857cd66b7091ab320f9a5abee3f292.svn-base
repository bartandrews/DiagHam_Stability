////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated to particles on a disk with         //
//     SU(2) spin and a generic interaction defined by its pseudopotential    //
//                                                                            //
//                        last modification : 26/11/2012                      //
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


#include "Hamiltonian/ParticleOnDiskWithSpinGenericHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanDiskCoefficients.h"

#include "Hamiltonian/ParticleOnDiskGenericHamiltonian.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

ParticleOnDiskWithSpinGenericHamiltonian::ParticleOnDiskWithSpinGenericHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// totalLz = system total value
// architecture = architecture to use for precalculation
// pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialUpDown =  one-body tunnelling potential (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnDiskWithSpinGenericHamiltonian::ParticleOnDiskWithSpinGenericHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, int totalLz,
										   double** pseudoPotential, double* onebodyPotentialUpUp, 
										   double* onebodyPotentialDownDown, double* onebodyPotentialUpDown,
										   AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
										   char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->TwoParticleLzMax = 2 * this->LzMax;
  if (totalLz < this->TwoParticleLzMax)
    this->TwoParticleLzMax = totalLz;
  this->OneBodyTermFlag = false;
  this->Architecture = architecture;
  this->PseudoPotentials = new double* [3];
  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;
  this->MaxPseudoPotentials = new int[3];
  for (int j = 0; j < 3; ++j)
    {
      this->PseudoPotentials[j] = new double [this->LzMax + this->NbrLzValue];
      this->MaxPseudoPotentials[j] = 0;
      for (int i = 0; i <= (2 * this->LzMax); ++i)
	{
	  this->PseudoPotentials[j][i] = pseudoPotential[j][i];
	  if (this->PseudoPotentials[j][i] != 0.0)
	    this->MaxPseudoPotentials[j] = i;
	}
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

ParticleOnDiskWithSpinGenericHamiltonian::~ParticleOnDiskWithSpinGenericHamiltonian() 
{
  for (int j = 0; j < 3; ++j)
    delete[] this->PseudoPotentials[j];
  delete[] this->PseudoPotentials;
  delete[] this->MaxPseudoPotentials;
}

// evaluate all interaction factors
//   

void ParticleOnDiskWithSpinGenericHamiltonian::EvaluateInteractionFactors()
{
  this->M1IntraValue = 0;
  this->M1InterValue = 0;
  
  ClebschGordanDiskCoefficients CGCoefficients(this->LzMax);

  long TotalNbrInteractionFactors = 0;

  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  double TmpCoefficient = 0.0;

  this->NbrInterSectorSums = 2 * this->LzMax + 1;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  for (int m1 = 0; m1 <= this->LzMax; ++m1)
    for (int m2 = 0; m2 <= this->LzMax; ++m2)
      if ((m1 + m2) <= this->TwoParticleLzMax)
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
	if ((m1 + m2) <= this->TwoParticleLzMax)
	  {
	    this->InterSectorIndicesPerSum[(m1 + m2)][this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
	    this->InterSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrInterSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
	    ++this->NbrInterSectorIndicesPerSum[(m1 + m2)];
	  }
      }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NbrIntraSectorSums = 2 * this->LzMax - 1;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  if ((m1 + m2) <= this->TwoParticleLzMax)
	    {
	      ++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
	    }
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 < this->LzMax; ++m1)
	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	  {
	    if ((m1 + m2) <= this->TwoParticleLzMax)
	      {
		this->IntraSectorIndicesPerSum[(m1 + m2) - 1][this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1] = m1;
		this->IntraSectorIndicesPerSum[(m1 + m2) - 1][1 + (this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1)] = m2;
		++this->NbrIntraSectorIndicesPerSum[(m1 + m2) - 1];
	      }
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
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int MaxSum = m1 + m2;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  this->InteractionFactorsupup[i][Index] = 0.0;
		  this->InteractionFactorsdowndown[i][Index] = 0.0;
		  for (int j = 1; j <= MaxSum; j += 2)
		    {
		      TmpCoefficient =  (CGCoefficients.GetCoefficient(m1, m2, j) * CGCoefficients.GetCoefficient(m3, m4, j));
		      this->InteractionFactorsupup[i][Index] += this->PseudoPotentials[0][j] * TmpCoefficient;
		      this->InteractionFactorsdowndown[i][Index] += this->PseudoPotentials[1][j] * TmpCoefficient;
		    }
		  this->InteractionFactorsupup[i][Index] *= -4.0;
		  this->InteractionFactorsdowndown[i][Index] *= -4.0;
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      this->InteractionFactorsupdown = new double* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      double Factor = 2.0;
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int MaxSum = m1 + m2;
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  this->InteractionFactorsupdown[i][Index] = 0.0;
		  for (int j = 0; j <= MaxSum; ++j)
		    {
		      TmpCoefficient =  (CGCoefficients.GetCoefficient(m1, m2, j) * CGCoefficients.GetCoefficient(m3, m4, j));
		      this->InteractionFactorsupdown[i][Index] += this->PseudoPotentials[2][j] * TmpCoefficient;
		    }
		  this->InteractionFactorsupdown[i][Index] *= -Factor;
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      this->NbrIntraSectorSums = 2 * this->LzMax+1;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  if ((m1 + m2) <= this->TwoParticleLzMax)
	    {
	      ++this->NbrIntraSectorIndicesPerSum[m1 + m2];
	    }
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	  this->NbrIntraSectorIndicesPerSum[i] = 0;
	}
      for (int m1 = 0; m1 <= this->LzMax; ++m1)
	for (int m2 = m1; m2 <= this->LzMax; ++m2)
	  {
	    if ((m1 + m2) <= this->TwoParticleLzMax)
	      {
		this->IntraSectorIndicesPerSum[m1 + m2][this->NbrIntraSectorIndicesPerSum[m1 + m2] << 1] = m1;
		this->IntraSectorIndicesPerSum[m1 + m2][1 + (this->NbrIntraSectorIndicesPerSum[m1 + m2] << 1)] = m2;
		++this->NbrIntraSectorIndicesPerSum[m1 + m2];
	      }
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
	      int m1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int MaxSum = m1 + m2;
	      if ((this->MaxPseudoPotentials[0] < MaxSum) || (this->MaxPseudoPotentials[1] < MaxSum))
		{
		  if (this->MaxPseudoPotentials[0] >= this->MaxPseudoPotentials[1])
		    MaxSum = this->MaxPseudoPotentials[0];
		  else
		    MaxSum = this->MaxPseudoPotentials[1];
		}
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  this->InteractionFactorsupup[i][Index] = 0.0;
		  this->InteractionFactorsdowndown[i][Index] = 0.0;
		  for (int j = 0; j <= MaxSum; j += 2)
		    {
		      TmpCoefficient =  (CGCoefficients.GetCoefficient(m1, m2, j) * CGCoefficients.GetCoefficient(m3, m4, j));
		      this->InteractionFactorsupup[i][Index] += this->PseudoPotentials[0][j] * TmpCoefficient;
		      this->InteractionFactorsdowndown[i][Index] += this->PseudoPotentials[1][j] * TmpCoefficient;
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

      this->InteractionFactorsupdown = new double* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      double Factor = 2.0;
	      int m1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int m2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int MaxSum = m1 + m2;
	      if (this->MaxPseudoPotentials[2] < MaxSum)
		MaxSum = this->MaxPseudoPotentials[2];
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int m3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int m4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  this->InteractionFactorsupdown[i][Index] = 0.0;
		  for (int j = 0; j <= MaxSum; ++j)
		    {
		      TmpCoefficient =  (CGCoefficients.GetCoefficient(m1, m2, j) * CGCoefficients.GetCoefficient(m3, m4, j));
		      this->InteractionFactorsupdown[i][Index] += this->PseudoPotentials[2][j] * TmpCoefficient;
		    }
		  this->InteractionFactorsupdown[i][Index] *= Factor;
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }

 cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
 cout << "====================================" << endl;
}

