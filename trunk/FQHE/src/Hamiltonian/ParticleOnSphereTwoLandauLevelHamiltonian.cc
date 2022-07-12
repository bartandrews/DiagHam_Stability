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


#include "Hamiltonian/ParticleOnSphereTwoLandauLevelHamiltonian.h"
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
// lzmax = maximum Lz value reached by a particle in the state in the lower Landau level
// landauLevelIndexDifference = difference of indices between the lower and higher Landau levels
// pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
// cyclotronEnergy = cyclotron energy in e^2/(epsilon l_b) unit
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereTwoLandauLevelHamiltonian::ParticleOnSphereTwoLandauLevelHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, int landauLevelIndexDifference, 
										     double** pseudoPotential, double cyclotronEnergy,
										     AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = false;
  this->Architecture = architecture;
  this->TotalCyclotronEnergy = ((double) landauLevelIndexDifference) * cyclotronEnergy;
  this->LandauLevelIndexDifference = landauLevelIndexDifference;
  this->LzMaxDown = this->LzMax;
  this->LzMaxUp = this->LzMax + (2 * this->LandauLevelIndexDifference);
  this->PseudoPotentials = new double* [4];
  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;
  for (int j = 0; j < 2; ++j)
    {
      this->PseudoPotentials[j] = new double [this->NbrLzValue];
      cout << "Pseudo type "<<j<<endl;
      for (int i = 0; i < this->NbrLzValue; ++i)
	{
          this->PseudoPotentials[j][i] = pseudoPotential[j][this->LzMax - i];
          cout << this->PseudoPotentials[j][i] << " ";
        }
      cout<<endl;
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
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;

  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  cout << "done" << endl;
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout << "done" << endl;
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

ParticleOnSphereTwoLandauLevelHamiltonian::~ParticleOnSphereTwoLandauLevelHamiltonian() 
{
  for (int j = 0; j < 4; ++j)
    delete[] this->PseudoPotentials[j];
  delete[] this->PseudoPotentials;
}

// evaluate all interaction factors
//   

void ParticleOnSphereTwoLandauLevelHamiltonian::EvaluateInteractionFactors()
{
  ClebschGordanCoefficients ClebschDownDown (this->LzMaxDown, this->LzMaxDown);
  ClebschGordanCoefficients ClebschUpUp (this->LzMaxUp, this->LzMaxUp);
  ClebschGordanCoefficients ClebschUpDown (this->LzMaxUp, this->LzMaxDown);
  ClebschGordanCoefficients ClebshDownUp (this->LzMaxDown, this->LzMaxUp);

  int J = 2 * this->LzMax - 2;
  // int m4;
  double ClebschCoef;
  long TotalNbrInteractionFactors = 0;

  int Sign = 1;
  if (this->LzMax & 1)
    Sign = 0;
  double TmpCoefficient = 0.0;

  int MaxSum = 0;
  this->NbrUpUpSectorSums = 2 * this->LzMaxUp + 1;
  this->NbrDownDownSectorSums = 2 * this->LzMaxDown + 1;
  this->NbrDownDownSectorSums += 2 * this->LandauLevelIndexDifference;
  this->NbrUpDownSectorSums = this->LzMaxUp + this->LzMaxDown + this->LandauLevelIndexDifference + 1;
  this->NbrUpUpSectorIndicesPerSum = new int[this->NbrUpUpSectorSums];
  this->NbrUpDownSectorIndicesPerSum = new int [this->NbrUpDownSectorSums];
  this->NbrDownDownSectorIndicesPerSum = new int [this->NbrDownDownSectorSums];
  for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
    this->NbrUpUpSectorIndicesPerSum[i] = 0;
  for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
    this->NbrUpDownSectorIndicesPerSum[i] = 0;
  for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
    this->NbrDownDownSectorIndicesPerSum[i] = 0;

  for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
    for (int m2 = 0; m2 <= this->LzMaxUp; ++m2)
      this->NbrUpUpSectorIndicesPerSum[m1 + m2]++;
  for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
    for (int m2 = 0; m2 <= this->LzMaxDown; ++m2)
      this->NbrUpDownSectorIndicesPerSum[m1 + m2 + this->LandauLevelIndexDifference]++;
  for (int m1 = 0; m1 <= this->LzMaxDown; ++m1)
    for (int m2 = 0; m2 <= this->LzMaxDown; ++m2)
      this->NbrDownDownSectorIndicesPerSum[m1 + m2 + (2 * this->LandauLevelIndexDifference)]++;

  this->UpUpSectorIndicesPerSum = new int* [this->NbrUpUpSectorSums];
  for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
    if (this->NbrUpUpSectorIndicesPerSum[i] > 0)
    {
      this->UpUpSectorIndicesPerSum[i] = new int[2 * this->NbrUpUpSectorIndicesPerSum[i]];      
      this->NbrUpUpSectorIndicesPerSum[i] = 0;
    }
  this->UpDownSectorIndicesPerSum = new int* [this->NbrUpDownSectorSums];
  for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
    if (this->NbrUpDownSectorIndicesPerSum[i] > 0)
    {
      this->UpDownSectorIndicesPerSum[i] = new int[2 * this->NbrUpDownSectorIndicesPerSum[i]];      
      this->NbrUpDownSectorIndicesPerSum[i] = 0;
    }
  this->DownDownSectorIndicesPerSum = new int* [this->NbrUpDownSectorSums];
  for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
    if (this->NbrDownDownSectorIndicesPerSum[i] > 0)
    {
      this->DownDownSectorIndicesPerSum[i] = new int[2 * this->NbrDownDownSectorIndicesPerSum[i]];      
      this->NbrDownDownSectorIndicesPerSum[i] = 0;
    }


  for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
    for (int m2 = 0; m2 <= this->LzMaxUp; ++m2)
      {
 	this->UpUpSectorIndicesPerSum[(m1 + m2)][this->NbrUpUpSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
 	this->UpUpSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrUpUpSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
 	++this->NbrUpUpSectorIndicesPerSum[(m1 + m2)];
       }
  for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
    for (int m2 = 0; m2 <= this->LzMaxDown; ++m2)
      {
 	this->UpDownSectorIndicesPerSum[m1 + m2 + this->LandauLevelIndexDifference][this->NbrUpDownSectorIndicesPerSum[(m1 + m2 + this->LandauLevelIndexDifference)] << 1] = m1;
 	this->UpDownSectorIndicesPerSum[m1 + m2 + this->LandauLevelIndexDifference][1 + (this->NbrUpDownSectorIndicesPerSum[(m1 + m2 + this->LandauLevelIndexDifference)] << 1)] = m2 + this->LandauLevelIndexDifference;
 	++this->NbrUpDownSectorIndicesPerSum[(m1 + m2 + this->LandauLevelIndexDifference)];
       }
  for (int m1 = 0; m1 <= this->LzMaxDown; ++m1)
    for (int m2 = 0; m2 <= this->LzMaxDown; ++m2)
      {
 	this->DownDownSectorIndicesPerSum[m1 + m2 + (2 * this->LandauLevelIndexDifference)][this->NbrDownDownSectorIndicesPerSum[m1 + m2 + (2 * this->LandauLevelIndexDifference)] << 1] = m1 + this->LandauLevelIndexDifference;
 	this->DownDownSectorIndicesPerSum[m1 + m2 + (2 * this->LandauLevelIndexDifference)][1 + (this->NbrDownDownSectorIndicesPerSum[m1 + m2 + (2 * this->LandauLevelIndexDifference)] << 1)] = m2 + this->LandauLevelIndexDifference;
 	++this->NbrDownDownSectorIndicesPerSum[m1 + m2 + (2 * this->LandauLevelIndexDifference)];
       }

  this->InteractionFactorsUpUpUpUp = new double* [this->NbrUpUpSectorSums];
  this->InteractionFactorsUpDownUpUp = new double* [this->NbrUpUpSectorSums];
  this->InteractionFactorsDownDownUpUp = new double* [this->NbrUpUpSectorSums];
  for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
    {
      this->InteractionFactorsUpUpUpUp = 0;
      this->InteractionFactorsUpDownUpUp = 0;
      this->InteractionFactorsDownDownUpUp = 0;
    }
  this->InteractionFactorsDownDownDownDown = new double* [this->NbrDownDownSectorSums];
  this->InteractionFactorsUpDownDownDown = new double* [this->NbrDownDownSectorSums];
  this->InteractionFactorsUpUpDownDown = new double* [this->NbrDownDownSectorSums];
  for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
    {
      this->InteractionFactorsUpUpUpDown = 0;
      this->InteractionFactorsUpDownUpDown = 0;
      this->InteractionFactorsDownDownUpDown = 0;
    }
  this->InteractionFactorsUpUpUpDown = new double* [this->NbrUpDownSectorSums];
  this->InteractionFactorsUpDownUpDown = new double* [this->NbrUpDownSectorSums];
  this->InteractionFactorsDownDownUpDown = new double* [this->NbrUpDownSectorSums];
  for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
    {
      this->InteractionFactorsUpUpUpDown = 0;
      this->InteractionFactorsUpDownUpDown = 0;
      this->InteractionFactorsDownDownUpDown = 0;
    }


  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      //upup-upup term
       for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	 {
	  if (this->NbrUpUpSectorIndicesPerSum[i] > 0)
	    {
	      this->InteractionFactorsUpUpUpUp[i] = new double[this->NbrUpUpSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = (this->UpUpSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp;
		  int m2 = (this->UpUpSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxUp;
		  for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
		      int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxUp;
		      ClebschUpUp.InitializeCoefficientIterator(m1, m2);
		      this->InteractionFactorsUpUpUpUp[i][Index] = 0.0;
		      while (ClebschUpUp.Iterate(J, ClebschCoef))
			{
			  if (((J >> 1) & 1) == Sign)
			    {
			      TmpCoefficient = ClebschCoef * ClebschUpUp.GetCoefficient(m3, m4, J);
			      this->InteractionFactorsUpUpUpUp[i][Index] += this->PseudoPotentials[0][J >> 1] * TmpCoefficient;
			    }
			}
		      this->InteractionFactorsUpUpUpUp[i][Index] *= -4.0;
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
 		}
 	    }
        }
       // downdown-downdown term
       for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
 	{
	  if (this->NbrDownDownSectorIndicesPerSum[i] > 0)
	    {
	      this->InteractionFactorsDownDownDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrDownDownSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - this->LandauLevelIndexDifference;
		  int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - this->LandauLevelIndexDifference;
		  for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxDown - this->LandauLevelIndexDifference;
		      int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - this->LandauLevelIndexDifference;
		      ClebschDownDown.InitializeCoefficientIterator(m1, m2);
		      this->InteractionFactorsDownDownDownDown[i][Index] = 0.0;
		      while (ClebschDownDown.Iterate(J, ClebschCoef))
			{
			  if (((J >> 1) & 1) == Sign)
			    {
			      TmpCoefficient = ClebschCoef * ClebschDownDown.GetCoefficient(m3, m4, J);
			      this->InteractionFactorsDownDownDownDown[i][Index] += this->PseudoPotentials[1][J >> 1] * TmpCoefficient;
			    }
			}
		      this->InteractionFactorsDownDownDownDown[i][Index] *= -4.0;
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
 		}
 	    }
        }
     }

//FILL IN INTER

//       this->InteractionFactorsupdown = new double* [this->NbrInterSectorSums];
//       for (int i = 0; i < this->NbrInterSectorSums; ++i)
// 	{
// 	  this->InteractionFactorsupdown[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
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
// 		  this->InteractionFactorsupdown[i][Index] = 0.0;
// 		  while (Clebsch.Iterate(J, ClebschCoef))
// 		    {
// 		      TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
// 		      this->InteractionFactorsupdown[i][Index] += this->PseudoPotentials[2][J >> 1] * TmpCoefficient;
// 		    }
// 		  this->InteractionFactorsupdown[i][Index] *= -Factor;
// 		  ++TotalNbrInteractionFactors;
// 		  ++Index;
// 		}
// 	    }
// 	}

// //*********************************************************************************************************
// //********************************     M I X E D     T E R M S ***********************************************
// //*********************************************************************************************************
// // Mixed term indices that have the same structure as intra terms

//      this->NbrMixedIntraSectorSums = 2 * this->LzMax - 1;
//       this->NbrMixedIntraSectorIndicesPerSum = new int[this->NbrMixedIntraSectorSums];
//       for (int i = 0; i < this->NbrMixedIntraSectorSums; ++i)
// 	this->NbrMixedIntraSectorIndicesPerSum[i] = 0;      
//       for (int m1 = 0; m1 < this->LzMax; ++m1)
// 	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
// 	  ++this->NbrMixedIntraSectorIndicesPerSum[(m1 + m2) - 1];
//       this->MixedIntraSectorIndicesPerSum = new int* [this->NbrMixedIntraSectorSums];
//       for (int i = 0; i < this->NbrMixedIntraSectorSums; ++i)
// 	{
// 	  this->MixedIntraSectorIndicesPerSum[i] = new int[2 * this->NbrMixedIntraSectorIndicesPerSum[i]];      
// 	  this->NbrMixedIntraSectorIndicesPerSum[i] = 0;
// 	}
//       for (int m1 = 0; m1 < this->LzMax; ++m1)
// 	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
// 	  {
// 	    this->MixedIntraSectorIndicesPerSum[(m1 + m2) - 1][this->NbrMixedIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1] = m1;
// 	    this->MixedIntraSectorIndicesPerSum[(m1 + m2) - 1][1 + (this->NbrMixedIntraSectorIndicesPerSum[(m1 + m2) - 1] << 1)] = m2;
// 	    ++this->NbrMixedIntraSectorIndicesPerSum[(m1 + m2) - 1];
// 	  }

// //Fill in the matrix elements

//       this->InteractionFactorsmixedintra= new double* [this->NbrMixedIntraSectorSums];
//       for (int i = 0; i < this->NbrMixedIntraSectorSums; ++i)
// 	{
// 	  this->InteractionFactorsmixedintra[i] = new double[this->NbrMixedIntraSectorIndicesPerSum[i] * this->NbrMixedIntraSectorIndicesPerSum[i]];
// 	  int Index = 0;
// 	  for (int j1 = 0; j1 < this->NbrMixedIntraSectorIndicesPerSum[i]; ++j1)
// 	    {
// 	      int m1 = (this->MixedIntraSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
// 	      int m2 = (this->MixedIntraSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
// 	      for (int j2 = 0; j2 < this->NbrMixedIntraSectorIndicesPerSum[i]; ++j2)
// 		{
// 		  int m3 = (this->MixedIntraSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
// 		  int m4 = (this->MixedIntraSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
// 		  Clebsch.InitializeCoefficientIterator(m1, m2);
// 		  this->InteractionFactorsmixedintra[i][Index] = 0.0;
// 		  while (Clebsch.Iterate(J, ClebschCoef))
// 		    {
// 		      if (((J >> 1) & 1) == Sign)
// 			{
// 			  TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
// 			  this->InteractionFactorsmixedintra[i][Index] += this->PseudoPotentials[3][J >> 1] * TmpCoefficient;
// 			}
// 		    }
// 		  this->InteractionFactorsmixedintra[i][Index] *= -4.0;
// 		  TotalNbrInteractionFactors += 2;
// 		  ++Index;
// 		}
// 	    }
//         }


// //Mixed term indices that have the same structure as inter terms


//     this->NbrMixedInterSectorSums = 2 * this->LzMax + 1;
//     this->NbrMixedInterSectorIndicesPerSum = new int[this->NbrMixedInterSectorSums];
//     for (int i = 0; i < this->NbrMixedInterSectorSums; ++i)
//      this->NbrMixedInterSectorIndicesPerSum[i] = 0;
//     for (int m1 = 0; m1 <= this->LzMax; ++m1)
//      for (int m2 = 0; m2 <= this->LzMax; ++m2)
//        ++this->NbrMixedInterSectorIndicesPerSum[m1 + m2];      
//     this->MixedInterSectorIndicesPerSum = new int* [this->NbrMixedInterSectorSums];
//     for (int i = 0; i < this->NbrMixedInterSectorSums; ++i)
//      {
//        this->MixedInterSectorIndicesPerSum[i] = new int[2 * this->NbrMixedInterSectorIndicesPerSum[i]];      
//        this->NbrMixedInterSectorIndicesPerSum[i] = 0;
//      }
//     for (int m1 = 0; m1 <= this->LzMax; ++m1)
//       for (int m2 = 0; m2 <= this->LzMax; ++m2)
//        {
// 	 this->MixedInterSectorIndicesPerSum[(m1 + m2)][this->NbrMixedInterSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
// 	 this->MixedInterSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrMixedInterSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
// 	 ++this->NbrMixedInterSectorIndicesPerSum[(m1 + m2)];
//        }


// //Fill in the matrix elements


//       this->InteractionFactorsmixedinter = new double* [this->NbrMixedInterSectorSums];
//       for (int i = 0; i < this->NbrMixedInterSectorSums; ++i)
// 	{
// 	  this->InteractionFactorsmixedinter[i] = new double[this->NbrMixedInterSectorIndicesPerSum[i] * this->NbrMixedInterSectorIndicesPerSum[i]];
// 	  int Index = 0;
// 	  for (int j1 = 0; j1 < this->NbrMixedInterSectorIndicesPerSum[i]; ++j1)
// 	    {
// 	      double Factor = 2.0;
// 	      int m1 = (this->MixedInterSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMax;
// 	      int m2 = (this->MixedInterSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMax;
// 	      for (int j2 = 0; j2 < this->NbrMixedInterSectorIndicesPerSum[i]; ++j2)
// 		{
// 		  int m3 = (this->MixedInterSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMax;
// 		  int m4 = (this->MixedInterSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMax;
// 		  Clebsch.InitializeCoefficientIterator(m1, m2);
// 		  this->InteractionFactorsmixedinter[i][Index] = 0.0;
// 		  while (Clebsch.Iterate(J, ClebschCoef))
// 		    {
// 		      TmpCoefficient = ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
// 		      this->InteractionFactorsmixedinter[i][Index] += this->PseudoPotentials[3][J >> 1] * TmpCoefficient;
// 		    }
// 		  //The sign of the Factor is opposite to UpDown case
// 		  this->InteractionFactorsmixedinter[i][Index] *= Factor;
// 		  ++TotalNbrInteractionFactors;
// 		  ++Index;
// 		}
// 	    }
// 	}




// //*********************************************************************************************************
// //********************************     M I X E D     T E R M S ***********************************************
// //*********************************************************************************************************

//     }
//   else
//     {
//     }

  this->NbrOneBodyInteractionFactorsUpUp = this->LzMaxUp + 1;
  this->NbrOneBodyInteractionFactorsUpDown = 0;
  this->NbrOneBodyInteractionFactorsDownUp = 0;
  this->NbrOneBodyInteractionFactorsDownDown = 0;
  this->OneBodyInteractionFactorsUpUp = new double [this->NbrOneBodyInteractionFactorsUpUp];
  this->OneBodyInteractionFactorsUpDown = 0;
  this->OneBodyInteractionFactorsDownUp = 0;
  this->OneBodyInteractionFactorsDownDown = 0;
  this->OneBodyMValuesUpUp = new int[this->NbrOneBodyInteractionFactorsUpUp];
  this->OneBodyMValuesUpDown = 0;
  this->OneBodyMValuesDownUp = 0;
  this->OneBodyMValuesDownDown = 0;
  for (int i = 0; i <= this->LzMaxUp; ++i)
    {
      this->OneBodyMValuesUpUp[i] = i;
      this->OneBodyInteractionFactorsUpUp[i] = this->TotalCyclotronEnergy;
    }

 cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
 cout << "====================================" << endl;
}

