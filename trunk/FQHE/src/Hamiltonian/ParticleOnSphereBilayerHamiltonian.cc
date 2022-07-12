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


#include "Hamiltonian/ParticleOnSphereBilayerHamiltonian.h"
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
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialUpDown =  one-body tunnelling potential (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
// chargingEnergy = charging energy for bilayers with d > 0
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereBilayerHamiltonian::ParticleOnSphereBilayerHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, 
										       double** pseudoPotential, double* onebodyPotentialUpUp, double* onebodyPotentialDownDown, double* onebodyPotentialUpDown, double chargingEnergy, AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = false;
  this->Architecture = architecture;
  this->PseudoPotentials = new double* [3];
  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;
  for (int j = 0; j < 3; ++j)
    {
      this->PseudoPotentials[j] = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->PseudoPotentials[j][i] = pseudoPotential[j][this->LzMax - i];
    }
  this->ChargingEnergy = chargingEnergy;
  cout << "Shifting the bilayer energies by " << this->ChargingEnergy << endl;
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

ParticleOnSphereBilayerHamiltonian::~ParticleOnSphereBilayerHamiltonian() 
{
}


void ParticleOnSphereBilayerHamiltonian::EvaluateInteractionFactors()
{
// this part of the code has been tested and is working (27/07/2007) but seems slower than the other method. This part is kept for testing purpose only

//   this->NbrIntraSectorSums = 0;
//   this->NbrInterSectorSums = 0;
//   int Lim;
//   int Min;
//   int Pos = 0;
//   ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
//   int J = 2 * this->LzMax - 2;
//   int m4;
//   double ClebschCoef;
//   double* TmpCoefficientupup = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
//   double* TmpCoefficientdowndown = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
//   double* TmpCoefficientupdown = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];

//   int Sign = 1;
//   if (this->LzMax & 1)
//     Sign = 0;
//   double MaxCoefficient = 0.0;

//   for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
//     for (int m2 =  -this->LzMax; m2 < m1; m2 += 2)
//       {
// 	Lim = m1 + m2 + this->LzMax;
// 	if (Lim > this->LzMax)
// 	  Lim = this->LzMax;
// 	Min = m1 + m2 - this->LzMax;
// 	if (Min < -this->LzMax)
// 	  Min = -this->LzMax;
// 	for (int m3 = Min; m3 <= Lim; m3 += 2)
// 	  {
// 	    Clebsch.InitializeCoefficientIterator(m1, m2);
// 	    m4 = m1 + m2 - m3;
// 	    TmpCoefficientupup[Pos] = 0.0;
// 	    TmpCoefficientdowndown[Pos] = 0.0;
// 	    while (Clebsch.Iterate(J, ClebschCoef))
// 	      {
// 		if (((J >> 1) & 1) == Sign)
// 		  {
// 		    TmpCoefficientupup[Pos] += this->PseudoPotentials[0][J >> 1] * ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
// 		    TmpCoefficientdowndown[Pos] += this->PseudoPotentials[1][J >> 1] * ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
// 		  }
// 	      }
// 	    ++Pos;
// 	  }
//       }

//   this->M12InteractionFactorsupup = new double [Pos];
//   this->M12InteractionFactorsdowndown = new double [Pos];

//   Pos = 0;
//   for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
//     for (int m2 =  -this->LzMax; m2 <= this->LzMax; m2 += 2)
//       {
// 	Lim = m1 + m2 + this->LzMax;
// 	if (Lim > this->LzMax)
// 	  Lim = this->LzMax;
// 	Min = m1 + m2 - this->LzMax;
// 	if (Min < -this->LzMax)
// 	  Min = -this->LzMax;
// 	for (int m3 = Min; m3 <= Lim; m3 += 2)
// 	  {
// 	    Clebsch.InitializeCoefficientIterator(m1, m2);
// 	    m4 = m1 + m2 - m3;
// 	    TmpCoefficientupdown[Pos] = 0.0;
// 	    while (Clebsch.Iterate(J, ClebschCoef))
// 	      {
// 		TmpCoefficientupdown[Pos] += this->PseudoPotentials[2][J >> 1] * ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
// 	      }
// 	    ++Pos;
// 	  }
//       }
//   this->M12InteractionFactorsupdown = new double [Pos];
  
//   this->NbrM12IntraIndices = (this->NbrLzValue * (this->NbrLzValue - 1)) / 2;
//   this->M1IntraValue = new int [this->NbrM12IntraIndices];
//   this->M2IntraValue = new int [this->NbrM12IntraIndices];
//   this->NbrM3IntraValues = new int [this->NbrM12IntraIndices];
//   this->M3IntraValues = new int* [this->NbrM12IntraIndices];
//   int TotalIndex = 0;
//   Pos = 0;
//   int TmpNbrInteractionFactors = 0;
//   double Factor = - 4.0;
//   for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
//     for (int m2 = 0; m2 < m1; ++m2)
//       {
// 	Lim = m1 + m2;
// 	if (Lim > this->LzMax)
// 	  Lim = this->LzMax;
// 	Min = m1 + m2 - this->LzMax;
// 	if (Min < 0)
// 	  Min = 0;
// 	this->M1IntraValue[TotalIndex] = m1;
// 	this->M2IntraValue[TotalIndex] = m2;	    
// 	this->NbrM3IntraValues[TotalIndex] = 0;
// 	for (int m3 = Min; m3 <= Lim; ++m3)
// 	  if ((2 * m3) > (m1 + m2))
// 	    ++this->NbrM3IntraValues[TotalIndex];
// 	if (this->NbrM3IntraValues[TotalIndex] > 0)
// 	  {
// 	    this->M3IntraValues[TotalIndex] = new int [this->NbrM3IntraValues[TotalIndex]];
// 	    int TmpIndex = 0;
// 	    for (int m3 = Min; m3 <= Lim; ++m3)
// 	      {
// 		if ((2 * m3) > (m1 + m2))
// 		  {
// 		    this->M3IntraValues[TotalIndex][TmpIndex] = m3;
// 		    this->M12InteractionFactorsupup[TmpNbrInteractionFactors] = Factor * TmpCoefficientupup[Pos];
// 		    this->M12InteractionFactorsdowndown[TmpNbrInteractionFactors] = Factor * TmpCoefficientdowndown[Pos];
// 		    ++TmpIndex;
// 		    ++TmpNbrInteractionFactors;
// 		  }
// 		++Pos;
// 	      }
// 	  }
// 	++TotalIndex;
//       }

//   this->NbrM12InterIndices = (this->NbrLzValue * this->NbrLzValue);
//   this->M1InterValue = new int [this->NbrM12InterIndices];
//   this->M2InterValue = new int [this->NbrM12InterIndices];
//   this->NbrM3InterValues = new int [this->NbrM12InterIndices];
//   this->M3InterValues = new int* [this->NbrM12InterIndices];
//   Factor = -2.0;
//   Pos = 0;
//   TotalIndex = 0;
//   TmpNbrInteractionFactors = 0;
//   for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
//     for (int m2 = 0; m2 < this->NbrLzValue; ++m2)
//       {
// 	Lim = m1 + m2;
// 	if (Lim > this->LzMax)
// 	  Lim = this->LzMax;
// 	Min = m1 + m2 - this->LzMax;
// 	if (Min < 0)
// 	  Min = 0;
// 	this->M1InterValue[TotalIndex] = m1;
// 	this->M2InterValue[TotalIndex] = m2;	    
// 	this->NbrM3InterValues[TotalIndex] = 0;
// 	for (int m3 = Min; m3 <= Lim; ++m3)
// 	  ++this->NbrM3InterValues[TotalIndex];
// 	if (this->NbrM3InterValues[TotalIndex] > 0)
// 	  {
// 	    this->M3InterValues[TotalIndex] = new int [this->NbrM3InterValues[TotalIndex]];
// 	    int TmpIndex = 0;
// 	    for (int m3 = Min; m3 <= Lim; ++m3)
// 	      {
// 		this->M3InterValues[TotalIndex][TmpIndex] = m3;
// 		this->M12InteractionFactorsupdown[TmpNbrInteractionFactors] = Factor * TmpCoefficientupdown[Pos];
// 		++TmpNbrInteractionFactors;
// 		++TmpIndex;
// 		++Pos;
// 	      }
// 	  }
// 	++TotalIndex;
//       }


  // int Lim;
  // int Min;
  // int Pos = 0;
  ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
  int J = 2 * this->LzMax - 2;
  // int m4;
  double ClebschCoef;
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

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
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
    }
  else
    {
    }


 cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
 cout << "====================================" << endl;
}

