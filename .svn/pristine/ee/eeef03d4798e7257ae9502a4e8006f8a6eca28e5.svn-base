////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          hardcore n-body interaction                       //
//                                                                            //
//                        last modification : 29/07/2014                      //
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


#include "config.h"
#include "Hamiltonian/ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <algorithm>
#include <sys/time.h>
#include <set>


// default constructor
//

ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian()
{
  this->NBodyPrefactor = 1.0;
  this->NbrPseudopotentials = 0;
  this->LaguerrePolynomials = 0;
  this->Pseudopotentials = 0;
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// nbrNBody = value of the n (i.e. the n-body interaction)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, int nbrNBody,  AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->NBodyValue = nbrNBody;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->TwoBodyFlag = false;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->NBodyPrefactor = 1.0;
  this->NbrPseudopotentials = 0;
  this->LaguerrePolynomials = 0;
  this->Pseudopotentials = 0;
  this->OneBodyInteractionFactors = 0;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  this->EvaluateExponentialFactors();
  
  
  this->NbrEntryPrecalculatedInteractionCoefficients1 = this->MaxMomentum;
  for (int i = 1; i < this->NBodyValue; ++i)
    this->NbrEntryPrecalculatedInteractionCoefficients1 *= this->MaxMomentum;
  this->PrecalculatedInteractionCoefficients = new double* [this->NbrEntryPrecalculatedInteractionCoefficients1];
  this->NbrEntryPrecalculatedInteractionCoefficients2 = this->NBodyValue*(2*this->NBodyValue - 1);
  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
    {      
      this->PrecalculatedInteractionCoefficients[m1] = new double [this->NbrEntryPrecalculatedInteractionCoefficients2];
    }
  char* InteractionCoefficientFileName = new char [512];
  sprintf (InteractionCoefficientFileName, "%dbodydelta_interactioncoefficient_2s_%d_ratio_%.10f.dat", this->NBodyValue, this->NbrLzValue, this->Ratio);
  if (IsFile(InteractionCoefficientFileName))
    {
      ifstream File;
      File.open(InteractionCoefficientFileName, ios::binary | ios::in);
      if (!File.is_open())
	{
	  cout << "cannot open " << InteractionCoefficientFileName << endl;
	}
      else
	{
	  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
	    ReadBlockLittleEndian(File, this->PrecalculatedInteractionCoefficients[m1], this->NbrEntryPrecalculatedInteractionCoefficients2);
	  File.close();
	}
    }
  else
    {
      ofstream File;
      File.open(InteractionCoefficientFileName, ios::binary | ios::out);
      if (!File.is_open())
	{
	  cout << "cannot create " << InteractionCoefficientFileName << endl;
	}
      else
	{  
	    
	  int** TmpIndices = new int* [this->NbrEntryPrecalculatedInteractionCoefficients1];
	  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
	    {
	      TmpIndices[m1] = new int [this->NBodyValue];
	      int Tmp = m1;
	      for (int i = 0; i < this->NBodyValue; ++i)
		{
		  TmpIndices[m1][i] = Tmp % this->MaxMomentum;
		  Tmp /= this->MaxMomentum;
		}
	    }
	  	  
	  cout << "this->NbrEntryPrecalculatedInteractionCoefficients1 = "  << this->NbrEntryPrecalculatedInteractionCoefficients1 << endl;
	  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
	    {
	      int TmpIndex = 0;
	      for (int j = 0; j < 2*this->NBodyValue - 1; ++j)
	      {
		int momentumTransfer = j - this->NBodyValue + 1;
		double* Coefficient = this->EvaluateInteractionCoefficientCreation(TmpIndices[m1], momentumTransfer);
		for (int g = 0; g < this->NBodyValue; ++g)
		{
		  this->PrecalculatedInteractionCoefficients[m1][TmpIndex] = Coefficient[g];
		  if (TmpIndex != (this->NBodyValue*j + g))
		    cout << TmpIndex << " " << (this->NBodyValue*j + g) << endl;
		  ++TmpIndex;
		}
		delete[] Coefficient;
	      }
	      WriteBlockLittleEndian(File, this->PrecalculatedInteractionCoefficients[m1], this->NbrEntryPrecalculatedInteractionCoefficients2);
	    }

	  
	  delete[] InteractionCoefficientFileName;
	  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
	    delete[] TmpIndices[m1];
	  delete[] TmpIndices;
	  File.close();
	}
    }
  
  
  
  this->HamiltonianShift = 0.0;
  this->EvaluateInteractionFactors();
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout << "fast memory = ";
	  PrintMemorySize(cout,TmpMemory)<<endl;
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

ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::~ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian()
{
  for (int i = 0; i < this->NbrEntryPrecalculatedInteractionCoefficients1; ++i)
    delete[] this->PrecalculatedInteractionCoefficients[i];  
  delete[] this->PrecalculatedInteractionCoefficients;
  if (this->NbrPseudopotentials != 0)
    {
      delete[] this->LaguerrePolynomials;
      delete[] this->Pseudopotentials;
    }
}
  
// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}


void ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0l;
  this->GetIndices();
  
  int** LinearizedNBodySectorIndicesPerSum = new int* [this->NbrNBodySectorSums];
  for (int i = 0; i < this->NbrNBodySectorSums; ++i)
    {
      LinearizedNBodySectorIndicesPerSum[i] = new int [NbrNBodySectorIndicesPerSum[i]];
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  int Tmp = 0;
	  for (int k = this->NBodyValue - 1; k >= 0; --k)
	    {
	      Tmp *= this->MaxMomentum;
	      Tmp += this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
	    }
	  LinearizedNBodySectorIndicesPerSum[i][j1] = Tmp;
	}
      SortArrayUpOrdering<int>(LinearizedNBodySectorIndicesPerSum[i], this->NbrNBodySectorIndicesPerSum[i]);
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  int Tmp = LinearizedNBodySectorIndicesPerSum[i][j1];
	  for (int k = 0; k < this->NBodyValue; ++k)
	    { 
	      this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k] = Tmp % this->MaxMomentum;
	      Tmp /= this->MaxMomentum;
	    }	  
	}
    }

  if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
    {
      int* TmpNIndices =  new int [this->NBodyValue];
      int* TmpMIndices =  new int [this->NBodyValue];

      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      int TmpCanonicalMIndices;
      int TmpCanonicalNIndices;
      int TmpCanonicalSum;
      double TmpCanonicalPermutationCoefficient;

      long TmpNbrMatrixElements = 0l;

      int sign = 1;
      if ((this->NBodyValue == 2) || (this->NBodyValue == 5))
	sign = -1;

      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 <= j1; ++j2)
		{
		  for (int k = 0; k < this->NBodyValue; ++k)
		    {
		      TmpNIndices[k]  = this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
		      TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + k];
		    }
 		  if (this->FindCanonicalIndices(TmpMIndices, TmpNIndices, LinearizedNBodySectorIndicesPerSum[i][j2], LinearizedNBodySectorIndicesPerSum[i][j1], i, TmpCanonicalMIndices, TmpCanonicalNIndices, TmpCanonicalSum, TmpCanonicalPermutationCoefficient) == true)
 		    {
		      ++TmpNbrMatrixElements;
		    }
		}
	    }
	}
      int* TmpMatrixElementIIndices = new int [TmpNbrMatrixElements];
      int* TmpMatrixElementJ2Indices = new int [TmpNbrMatrixElements];
      int* TmpMatrixElementJ1Indices = new int [TmpNbrMatrixElements];
      Complex* TmpMatrixElement = new Complex[TmpNbrMatrixElements];
      TmpNbrMatrixElements = 0l;

      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 <= j1; ++j2)
		{
		  for (int k = 0; k < this->NBodyValue; ++k)
		    {
		      TmpNIndices[k] = this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
		      TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + k];
		    }
 		  if (this->FindCanonicalIndices(TmpMIndices, TmpNIndices, LinearizedNBodySectorIndicesPerSum[i][j2], LinearizedNBodySectorIndicesPerSum[i][j1], i, TmpCanonicalMIndices, TmpCanonicalNIndices, TmpCanonicalSum, TmpCanonicalPermutationCoefficient) == true)
 		    {
		      TmpMatrixElementIIndices[TmpNbrMatrixElements] = i;
		      TmpMatrixElementJ1Indices[TmpNbrMatrixElements] = j1;
		      TmpMatrixElementJ2Indices[TmpNbrMatrixElements] = j2;
		      ++TmpNbrMatrixElements;
		    }
		}
	    }
	}

      cout << "generating " << TmpNbrMatrixElements << " unique matrix elements" << endl;
      timeval TotalStartingTime;
      gettimeofday (&TotalStartingTime, 0);

      for (long i = 0l; i < TmpNbrMatrixElements; ++i)
	{
	  int TmpJ1 = TmpMatrixElementJ1Indices[i];
	  int TmpJ2 = TmpMatrixElementJ2Indices[i];
	  int MomentumSector = TmpMatrixElementIIndices[i];
	  for (int k = 0; k < this->NBodyValue; ++k)
	    {
	      TmpNIndices[k]  = this->NBodySectorIndicesPerSum[MomentumSector][(TmpJ1 * this->NBodyValue) + k];
	      TmpMIndices[k]  = this->NBodySectorIndicesPerSum[MomentumSector][(TmpJ2 * this->NBodyValue) + k];
	    }
	  TmpMatrixElement[i] = sign * this->EvaluateInteractionNIndexAntiSymmetrizedCoefficient(TmpNIndices, TmpMIndices);	
	}
      timeval TotalEndingTime;
      gettimeofday (&TotalEndingTime, 0);
      double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
		    (((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
      cout << "generation done in " << Dt << endl;
      
      TmpNbrMatrixElements = 0l;
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      Index = j1 * this->NbrNBodySectorIndicesPerSum[i];
	      for (int j2 = 0; j2 <= j1; ++j2)
		{
		  for (int k = 0; k < this->NBodyValue; ++k)
		    {
		      TmpNIndices[k]  = this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
		      TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + k];
		    }
 		  if (this->FindCanonicalIndices(TmpMIndices, TmpNIndices, LinearizedNBodySectorIndicesPerSum[i][j2], LinearizedNBodySectorIndicesPerSum[i][j1], i, TmpCanonicalMIndices, TmpCanonicalNIndices, TmpCanonicalSum, TmpCanonicalPermutationCoefficient) == true)
 		    {
		      this->NBodyInteractionFactors[i][Index] = this->NBodyPrefactor * TmpMatrixElement[TmpNbrMatrixElements];
		      ++TmpNbrMatrixElements;
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
		  else
		    {
		      int TmpJ1 = SearchInArray<int>(TmpCanonicalNIndices, LinearizedNBodySectorIndicesPerSum[TmpCanonicalSum], this->NbrNBodySectorIndicesPerSum[TmpCanonicalSum]);
		      int TmpJ2 = SearchInArray<int>(TmpCanonicalMIndices, LinearizedNBodySectorIndicesPerSum[TmpCanonicalSum], this->NbrNBodySectorIndicesPerSum[TmpCanonicalSum]);
		      this->NBodyInteractionFactors[i][Index] = TmpCanonicalPermutationCoefficient * this->NBodyInteractionFactors[TmpCanonicalSum][TmpJ1 * this->NbrNBodySectorIndicesPerSum[TmpCanonicalSum] + TmpJ2];
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
		}
	    }
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = j1 + 1; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  this->NBodyInteractionFactors[i][j1 * this->NbrNBodySectorIndicesPerSum[i] + j2] = this->NBodyInteractionFactors[i][j2 * this->NbrNBodySectorIndicesPerSum[i] + j1];
		  TotalNbrInteractionFactors++;		  
		}
	    }
	}
    }
  else
    {      
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;  
	  	  
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      int* TmpNIndices = &(this->NBodySectorIndicesPerSum[i][j1 * this->NBodyValue]);
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  int* TmpMIndices = &(this->NBodySectorIndicesPerSum[i][j2 * this->NBodyValue]);
		  this->NBodyInteractionFactors[i][Index] = this->NBodyPrefactor * this->EvaluateInteractionNIndexSymmetrizedCoefficient(TmpNIndices, TmpMIndices);
		  TotalNbrInteractionFactors++;
		  ++Index;
		  
		  
		}
	    }
	}
    }

  if (this->TwoBodyFlag == true)
    {
      double MaxCoefficient = 0.0;
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
	{
	  for (int i = 0; i < this->NbrSectorSums; ++i)
	    {
	      this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->SectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient   = (this->EvaluateTwoBodyInteractionCoefficient(m1, m2, m3, m4)
						 + this->EvaluateTwoBodyInteractionCoefficient(m2, m1, m4, m3)
						 - this->EvaluateTwoBodyInteractionCoefficient(m1, m2, m4, m3)
						 - this->EvaluateTwoBodyInteractionCoefficient(m2, m1, m3, m4));
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			MaxCoefficient = fabs(TmpCoefficient);
		}
		}
	    }
	  MaxCoefficient *= MACHINE_PRECISION;
	  for (int i = 0; i < this->NbrSectorSums; ++i)
	    {
	      this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->SectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient = (this->EvaluateTwoBodyInteractionCoefficient(m1, m2, m3, m4)
					       + this->EvaluateTwoBodyInteractionCoefficient(m2, m1, m4, m3)
					       - this->EvaluateTwoBodyInteractionCoefficient(m1, m2, m4, m3)
					       - this->EvaluateTwoBodyInteractionCoefficient(m2, m1, m3, m4));
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			{
			  this->InteractionFactors[i][Index] = TmpCoefficient;
			}
		      else
			{
			  this->InteractionFactors[i][Index] = 0.0;
			}
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
		}
	    }
	}
      else
	{
	  for (int i = 0; i < this->NbrSectorSums; ++i)
	    {
	      this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->SectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient = (this->EvaluateTwoBodyInteractionCoefficient(m1, m2, m3, m4)
					       + this->EvaluateTwoBodyInteractionCoefficient(m2, m1, m4, m3)
					       + this->EvaluateTwoBodyInteractionCoefficient(m1, m2, m4, m3)
					       + this->EvaluateTwoBodyInteractionCoefficient(m2, m1, m3, m4));
		      if (m3 == m4)
			TmpCoefficient *= 0.5;
		      if (m1 == m2)
			TmpCoefficient *= 0.5;
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			MaxCoefficient = fabs(TmpCoefficient);
		    }
		}
	    }
	  MaxCoefficient *= MACHINE_PRECISION;
	  for (int i = 0; i < this->NbrSectorSums; ++i)
	    {
	      this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = this->SectorIndicesPerSum[i][j1 << 1];
		  int m2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
		  for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = this->SectorIndicesPerSum[i][j2 << 1];
		      int m4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		      
		      double TmpCoefficient = (this->EvaluateTwoBodyInteractionCoefficient(m1, m2, m3, m4)
					       + this->EvaluateTwoBodyInteractionCoefficient(m2, m1, m4, m3)
					       + this->EvaluateTwoBodyInteractionCoefficient(m1, m2, m4, m3)
					       + this->EvaluateTwoBodyInteractionCoefficient(m2, m1, m3, m4));
		      if (m3 == m4)
			TmpCoefficient *= 0.5;
		      if (m1 == m2)
			TmpCoefficient *= 0.5;
		      if (fabs(TmpCoefficient) > MaxCoefficient)
			{
			  this->InteractionFactors[i][Index] = TmpCoefficient;
			}
		      else
			{
			  this->InteractionFactors[i][Index] = 0.0;
			}
		      TotalNbrInteractionFactors++;
		      ++Index;
		    }
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}




// evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term (factor corresponding to the creation or the annihilation operators only) for each integer modulo the NBodyValue
//
// mIndices = array that contains the creation indices
// tmpSum = sum of all creation indices
// momentumTransfer = momentum transfer operated by the \prod_i a+_mi \prod_j a_nj operator, in units of the number of flux quanta
// return value = array of numerical coefficients  

double* ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficientCreation(int* mIndices, int momentumTransfer)
{
  double* Coefficient = new double [this->NBodyValue];
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  int tmpSum = 0;
  for (int i = 0; i < this->NBodyValue; ++i)
    tmpSum += mIndices[i];
  
  int* momFactor = new int [this->NBodyValue - 1];
  int* TmpIndices = new int [(this->NBodyValue - 1)];
  int* countIter = new int[(this->NBodyValue - 1)];
  for (int i = 0; i < this->NBodyValue - 1; ++i)
    momFactor[i] = this->NBodyValue * mIndices[i] -  tmpSum; 
  
  for (int i = 0; i < this->NBodyValue; ++i)
  {
    for (int j = 0; j < this->NBodyValue - 1; ++j)
      {
	int coef1 = (momFactor[j] / this->NbrLzValue) - (momFactor[j] / this->NbrLzValue) % this->NBodyValue;   
	TmpIndices[j] = (double) (((i + momentumTransfer) % this->NBodyValue) - coef1);	
	countIter[j] = 0;
      }
    Coefficient[i] = this->EvaluateGaussianSum(0, TmpIndices, 0, countIter, momFactor);
  }
  delete[] TmpIndices;
  delete[] momFactor;
  delete[] countIter;
  
  return Coefficient;
}



// evaluate the N nested infinite sums of EvaluateInteractionCoefficientCreation
//
// nBodyValue = current index being incremented 
//TmpIndices = array containing the current value of all indices
// Sum = current value of the sum
// countIter = array of integers giving, for each TmpIndices[i], the number of times it has been incremented
// momFactor = array of indices that contains the information of the creation (or annihilation) indices
//return value = value of the coefficient
double ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::EvaluateGaussianSum(int nBodyValue, int* TmpIndices, double Sum, int* countIter, int* momFactor)
{
  if (nBodyValue == this->NBodyValue - 1)
    return 0.0;
  
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  double normalizationCoefficient = pow(DoubleNbrLzValue,((double) (this->NBodyValue + 1)) / 4.0);
  normalizationCoefficient = 1.0;
  double PIOnM = M_PI / DoubleNbrLzValue ;
  double Factor = 2.0*M_PI*this->Ratio / (DoubleNbrLzValue * ((double)(this->NBodyValue))* ((double)(this->NBodyValue)));
  int MinIter = 3;
  int TmpIndex = TmpIndices[0];
  double ExpFactor;

  
  for (int i = nBodyValue; i < this->NBodyValue - 1; ++i)
    {
      TmpIndices[i] = TmpIndex;
      countIter[i] = 0;
    }
  ExpFactor = 0.0;
  for (int j = 0; j < this->NBodyValue - 1; ++j)
    for (int k = j; k < this->NBodyValue - 1; ++k)
      ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
  double Coefficient = normalizationCoefficient * exp(-Factor*ExpFactor);
  
  while ((Coefficient + Sum != Sum) || (countIter[nBodyValue] < MinIter))
  {
    countIter[nBodyValue] += 1;
    Sum += this->EvaluateGaussianSum(nBodyValue + 1, TmpIndices, 0.0, countIter, momFactor);
    if (nBodyValue == this->NBodyValue - 2)
      Sum += Coefficient;
    TmpIndices[nBodyValue] += this->NBodyValue;
    for (int i = nBodyValue + 1; i < this->NBodyValue - 1; ++i)
      TmpIndices[i] = TmpIndex;
    ExpFactor = 0;
    for (int j = 0; j < this->NBodyValue - 1; ++j)
      for (int k = j; k < this->NBodyValue - 1; ++k)
	ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
    Coefficient = normalizationCoefficient * exp(-Factor*ExpFactor);
    
  }
  
  TmpIndices[nBodyValue] = TmpIndex - this->NBodyValue;
  countIter[nBodyValue] = 0;
  for (int i = nBodyValue + 1; i < this->NBodyValue - 1; ++i)
  {
    TmpIndices[i] = TmpIndex;
    countIter[i] = 0;
  }
    
  ExpFactor = 0.0;
  for (int j = 0; j < this->NBodyValue - 1; ++j)
    for (int k = j; k < this->NBodyValue - 1; ++k)
      ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
  Coefficient = normalizationCoefficient * exp(-Factor*ExpFactor);
  
  
  while ((Coefficient + Sum != Sum) || (countIter[nBodyValue] < MinIter))
  {
    countIter[nBodyValue] += 1;
    Sum += this->EvaluateGaussianSum(nBodyValue + 1, TmpIndices, 0.0, countIter, momFactor);
    if (nBodyValue == this->NBodyValue - 2)
      Sum += Coefficient;
    TmpIndices[nBodyValue] -= this->NBodyValue;
    ExpFactor = 0;
    for (int i = nBodyValue + 1; i < this->NBodyValue - 1; ++i)
      TmpIndices[i] = TmpIndex;
    for (int j = 0; j < this->NBodyValue - 1; ++j)
      for (int k = j; k < this->NBodyValue - 1; ++k)
	ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
    Coefficient = normalizationCoefficient * exp(-Factor*ExpFactor);    
  }
  
  return Sum;
  
}
  

// Evaluate gaussian sum for the three body interaction (for test purposes only)
//
// momFactor = array of indices that contains the information of the creation (or annihilation) indices
// TmpIndices = array of indices that gives the initial indices that will be incremented in the sum
// return value = value of the sum

double ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::EvaluateGaussianSum(int* momFactor, int* TmpIndices)
{
  double Sum = 0;
  double Precision = 1.0;
  int MinIter = 5;
  int countIter1 = 0;
  int countIter2 = 0;
  
  
  int TmpIndex1 = TmpIndices[0];
  int TmpIndex2 = TmpIndices[1];
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  double PIOnM = M_PI / DoubleNbrLzValue ;
  double Factor = 2.0*M_PI*this->Ratio / (DoubleNbrLzValue * ((double)(this->NBodyValue))* ((double)(this->NBodyValue)));
  double ExpFactor1;
  double ExpFactor2;
  
  ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
  ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
      
  double Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
  
  
  while ((Sum + Precision*Coefficient != Sum) || (countIter1 < MinIter))
  {
    countIter1 += 1;
    countIter2 = 0;
    while ((Sum + Precision*Coefficient != Sum) || (countIter2 < MinIter))
    {
      countIter2 += 1;
      Sum += Coefficient;
      TmpIndices[0] += this->NBodyValue;
      ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
      ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
      Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
      
    }
    
    TmpIndices[0] = TmpIndex1 - this->NBodyValue;
    countIter2 = 0;
    ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
    ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
    Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
    
    while ((Sum + Precision*Coefficient != Sum) || (countIter2 < MinIter))
    {
      countIter2 += 1;
      Sum += Coefficient;
      TmpIndices[0] -= this->NBodyValue;
      ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
      ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
      Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
      
    }
    
    TmpIndices[0] = TmpIndex1;
    TmpIndices[1] += this->NBodyValue;
    ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
    ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
    Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
  }
  
  countIter1 = 0;
  TmpIndices[0] = TmpIndex1;
  TmpIndices[1] = TmpIndex2 - this->NBodyValue;
  ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
  ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
  Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
  
  while ((Sum + Precision*Coefficient != Sum) || (countIter1 < MinIter))
  {
    countIter1 += 1;
    countIter2 = 0;
    while ((Sum + Precision*Coefficient != Sum) || (countIter2 < MinIter))
    {
      countIter2 += 1;
      Sum += Coefficient;
      TmpIndices[0] += this->NBodyValue;
      ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
      ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
      Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
      
    }
    
    TmpIndices[0] = TmpIndex1 - this->NBodyValue;
    countIter2 = 0;
    ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
    ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
    Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
    
    while ((Sum + Precision*Coefficient != Sum) || (countIter2 < MinIter))
    {
      countIter2 += 1;
      Sum += Coefficient;
      TmpIndices[0] -= this->NBodyValue;
      ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
      ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
      Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
      
    }
    
    TmpIndices[0] = TmpIndex1;
    TmpIndices[1] -= this->NBodyValue;
    ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
    ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
    Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
  }
  
  return Sum;
}


	
// evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term
//
// creationCoefficients = array that contains the creation coefficients
// annihilationCoefficients = array that contains the annihilation coefficients
// nbrPermutations1 = number of permutations of the creation indexes
// nbrPermutations2 = number of permutations of the annihilation indexes
// return value = numerical coefficient  

double ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(double* creationCoefficients, double* annihilationCoefficients)
{
  double Coefficient = 0.0;
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  FactorialCoefficient FactorialNBody;
  FactorialNBody.SetToOne();
  FactorialNBody.FactorialMultiply(this->NBodyValue);

 
  for (int j = 0; j < this->NBodyValue; ++j)
    Coefficient += creationCoefficients[j]*annihilationCoefficients[j];

  for (int i = 0; i < this->NBodyValue; ++i)
    Coefficient /= (M_PI * DoubleNbrLzValue);
  
  return (Coefficient/ (FactorialNBody.GetNumericalValue()));
}

  
// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::EvaluateTwoBodyInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  double Coefficient = 1.0;
  double PIOnM = M_PI / ((double) this->MaxMomentum);
  double Factor =  - ((double) (m1-m3)) * PIOnM * 2.0;
  double Sum = 0.0;
  double N2 = (double) (m1 - m4);
  double N1;
  double Q2;
  double Precision;
//  cout << "new coef====================================" << m1 << " "  << m2 << " "  << m3 << " "  << m4 << endl;
  while ((fabs(Sum) + fabs(Coefficient)) != fabs(Sum))
    {
      N1 = 1.0;
      Q2 = this->Ratio * N2 * N2;
      if (N2 != 0.0)
	{
	  Coefficient = this->GetVofQ(PIOnM*Q2);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnM*Q2); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM*Q2);
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
	  Coefficient =  this->GetVofQ(PIOnM*Q2);
	  Precision = Coefficient;
	}
      else
	{
	  Coefficient = this->GetVofQ(PIOnM*Q2); // yields non-zero terms only for non-singular interactions
	  Precision = 1.0;
	}
      while ((fabs(Coefficient) + Precision) != fabs(Coefficient))
	{
	  Q2 = this->InvRatio * N1 * N1 + this->Ratio * N2 * N2;
	  Precision = 2.0 * this->GetVofQ(PIOnM*Q2);
	  Coefficient += Precision * cos (N1 * Factor);
	  N1 += 1.0;
	}
      Sum += Coefficient;
      N2 -= this->MaxMomentum;
    }
  return (Sum / (2.0 * this->MaxMomentum));
}

// get Fourier transform of the interaction
//
// Q2_half = one half of q² value
// retrun value =Fourier transform of ithe nteraction

double ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::GetVofQ(double Q2_half)
{
  double Result = 0.0;
  double Q2 = 2.0 * Q2_half;
//   if ((this->HaveCoulomb)&&(Q2_half!=0.0))
//     {
//       Result=GETSQR(this->FormFactor(Q2_half)) / sqrt(Q2);
//     }
//   else
//     Result=0.0;
  for (int i = 0; i < this->NbrPseudopotentials; ++i)
    if (this->Pseudopotentials[i]!=0.0)
      Result += 2.0 * this->Pseudopotentials[i] * this->LaguerrePolynomials[i].PolynomialEvaluate(Q2);
  return Result * exp(-Q2_half);
}




// find the canonical index corresponding to a giving index
//
// index = linearized index 
// sign = reference on an optional sign resulting when finiding the canonical index
// indices = temporary array to store the full index
// nbrPermutations  = temporary array to store the numebr of permutations
// return value = canonical index

int ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::GetCanonicalIndex (int index, int& sign, int* indices, int* nbrPermutations)
{
  sign = 1;
  int g;
  int MomentumTransfer;
  this->GetIndicesFromLinearizedIndex (index, g, MomentumTransfer, indices);
  int TmpIndex = index;  
  int CountModulo;
  int TmpSign;
  nbrPermutations[0] = 0;
  SortArrayDownOrderingPermutationBubbleSort<int>(indices, this->NBodyValue, nbrPermutations[0]);
  int TmpSumGMomentumTransfer = (g + MomentumTransfer + this->NBodyValue) % this->NBodyValue;
  int SumGMomentumTransfer = TmpSumGMomentumTransfer;
  MomentumTransfer = - this->NBodyValue + 1;
  g = (TmpSumGMomentumTransfer + this->NBodyValue - 1) % this->NBodyValue;
  int CanonicalIndex = this->EvaluateLinearizedIndex(indices, g, MomentumTransfer);
  int IndexForCanonical = 0;
  for (int j = 0; j <  this->NbrLzValue - 1; ++j)
    {
      nbrPermutations[j + 1] = nbrPermutations[j];
      CountModulo = 0;
      for (int i = 0; i < this->NBodyValue; ++i)
	{
	  indices[i]++;
	  if (indices[i] == this->NbrLzValue)
	    {
	      CountModulo++;
	      indices[i] = 0;
	    }   
	}      
      
      if (CountModulo != 0)
	{
	  SortArrayDownOrderingPermutationBubbleSort<int>(indices, this->NBodyValue, nbrPermutations[j + 1]);	  
	}
      
      TmpSumGMomentumTransfer = (TmpSumGMomentumTransfer - CountModulo) % this->NBodyValue;
      MomentumTransfer = - this->NBodyValue + 1;
      g = (TmpSumGMomentumTransfer + this->NBodyValue - 1) % this->NBodyValue;
      TmpIndex = this->EvaluateLinearizedIndex(indices, g, MomentumTransfer);
      if (TmpIndex < CanonicalIndex)
	{
	  CanonicalIndex = TmpIndex;
	  IndexForCanonical = j + 1;
	}
    }
  
     
  if ((nbrPermutations[IndexForCanonical] % 2) == 0)
    sign = 1;
  else
    sign = -1;

  return CanonicalIndex;   
}


// test if creation/annihilation are in canonical form
//
// mIndices = array that contains the creation indices (will be modified)
// nIndices = array that contains the annihilation indices (will be modified)
// linearizedMIndices = linearized creation index
// linearizedNIndices = linearized annihilation index
// totalMomentum = momentum sector of the creation/annihilation indices
// canonicalMIndices = reference on the linearized creation index of the canonical form
// canonicalMIndices = reference on the linearized annihilation index of the canonical form
// canonicalTotalMomentum = reference on the momentum sector of the creation/annihilation indices of the canonical form
// canonicalPermutationCoefficient = additional permutation coefficient to get the canonical form
// return value = true if the indices are in canonical form

bool ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::FindCanonicalIndices(int* mIndices, int* nIndices, int linearizedMIndices, int linearizedNIndices, 
											   int totalMomentum, int& canonicalMIndices, int& canonicalNIndices, 
											   int& canonicalTotalMomentum, double& canonicalPermutationCoefficient)
{
//  return true;
  int* TmpMIndices = new int [this->NBodyValue];
  int* TmpNIndices = new int [this->NBodyValue];
  canonicalPermutationCoefficient = 1.0;
  canonicalMIndices = linearizedMIndices;
  canonicalNIndices = linearizedNIndices;
  canonicalTotalMomentum = totalMomentum;
  canonicalPermutationCoefficient = 1.0;
  bool IsCanonical = true;
  for (int k = 0; k < this->NBodyValue; ++k)
    {
      if (mIndices[k] > 0)
	{
	  TmpMIndices[k] = this->MaxMomentum - mIndices[k];
	}
      else
	{
	  TmpMIndices[k] = 0;
	}
      if (nIndices[k] > 0)
	{
	  TmpNIndices[k] = this->MaxMomentum - nIndices[k];
	}
      else
	{
	  TmpNIndices[k] = 0;
	}
    }
  int TmpNbrPermutation = 0;
  SortArrayDownOrderingPermutationBubbleSort<int>(TmpMIndices, this->NBodyValue, TmpNbrPermutation);
  SortArrayDownOrderingPermutationBubbleSort<int>(TmpNIndices, this->NBodyValue, TmpNbrPermutation);
  int TmpCanonicalTotalMomentum = 0;
  int TmpCanonicalMIndices = 0;
  int TmpCanonicalNIndices = 0;
  for (int k = this->NBodyValue - 1; k >= 0; --k)
    {
      TmpCanonicalMIndices *= this->MaxMomentum;
      TmpCanonicalMIndices += TmpMIndices[k];
      TmpCanonicalNIndices *= this->MaxMomentum;
      TmpCanonicalNIndices += TmpNIndices[k];
      TmpCanonicalTotalMomentum += TmpNIndices[k];
    }
  TmpCanonicalTotalMomentum %= this->MaxMomentum;
  if (TmpCanonicalMIndices > TmpCanonicalNIndices)
    {
      int Tmp = TmpCanonicalMIndices;
      TmpCanonicalMIndices = TmpCanonicalNIndices;
      TmpCanonicalNIndices = Tmp;
    }
  if ((TmpCanonicalTotalMomentum < canonicalTotalMomentum) || 
      ((TmpCanonicalTotalMomentum == canonicalTotalMomentum) && ((TmpCanonicalNIndices < canonicalNIndices) ||
								 ((TmpCanonicalNIndices == canonicalNIndices) && (TmpCanonicalMIndices < canonicalMIndices)))))
    {
      IsCanonical = false;
      canonicalMIndices = TmpCanonicalMIndices;
      canonicalNIndices = TmpCanonicalNIndices;
      canonicalTotalMomentum = TmpCanonicalTotalMomentum;
      if ((TmpNbrPermutation & 1) == 0)
	{
	  canonicalPermutationCoefficient = 1.0;
	}
      else
	{
	  canonicalPermutationCoefficient = -1.0;
	}
    }
  for (int i = 1; i <  this->MaxMomentum; ++i)
    {
      bool ResuffleMFlag = false;
      bool ResuffleNFlag = false;
      for (int k = 0; k < this->NBodyValue; ++k)
	{
	  TmpMIndices[k] = mIndices[k] - i;
	  if (TmpMIndices[k] < 0)
	    {
	      TmpMIndices[k] += this->MaxMomentum;
	      ResuffleMFlag = true;
	    }
	  TmpNIndices[k] = nIndices[k] - i;
	  if (TmpNIndices[k] < 0)
	    {
	      TmpNIndices[k] += this->MaxMomentum;
	      ResuffleNFlag = true;
	    }
	}
      TmpNbrPermutation = 0;
      if (ResuffleMFlag == true)
	SortArrayDownOrderingPermutationBubbleSort<int>(TmpMIndices, this->NBodyValue, TmpNbrPermutation);
      if (ResuffleNFlag == true)
	SortArrayDownOrderingPermutationBubbleSort<int>(TmpNIndices, this->NBodyValue, TmpNbrPermutation);
      TmpCanonicalTotalMomentum = 0;
      TmpCanonicalMIndices = 0;
      TmpCanonicalNIndices = 0;
      for (int k = this->NBodyValue - 1; k >= 0; --k)
	{
	  TmpCanonicalMIndices *= this->MaxMomentum;
	  TmpCanonicalMIndices += TmpMIndices[k];
	  TmpCanonicalNIndices *= this->MaxMomentum;
	  TmpCanonicalNIndices += TmpNIndices[k];
	  TmpCanonicalTotalMomentum += TmpNIndices[k];
	}
      TmpCanonicalTotalMomentum %= this->MaxMomentum;
      if (TmpCanonicalMIndices > TmpCanonicalNIndices)
	{
	  int Tmp = TmpCanonicalMIndices;
	  TmpCanonicalMIndices = TmpCanonicalNIndices;
	  TmpCanonicalNIndices = Tmp;
	}
      if ((TmpCanonicalTotalMomentum < canonicalTotalMomentum) || 
	  ((TmpCanonicalTotalMomentum == canonicalTotalMomentum) && ((TmpCanonicalNIndices < canonicalNIndices) ||
								     ((TmpCanonicalNIndices == canonicalNIndices) && (TmpCanonicalMIndices < canonicalMIndices)))))
	{
	  IsCanonical = false;
	  canonicalMIndices = TmpCanonicalMIndices;
	  canonicalNIndices = TmpCanonicalNIndices;
	  canonicalTotalMomentum = TmpCanonicalTotalMomentum;
	  if ((TmpNbrPermutation & 1) == 0)
	    {
	      canonicalPermutationCoefficient = 1.0;
	    }
	  else
	    {
	      canonicalPermutationCoefficient = -1.0;
	    }
 	}
      for (int k = 0; k < this->NBodyValue; ++k)
	{
	  if (TmpMIndices[k] > 0)
	    {
	      TmpMIndices[k] = this->MaxMomentum - TmpMIndices[k];
	    }
	  else
	    {
	      TmpMIndices[k] = 0;
	    }
	  if (TmpNIndices[k] > 0)
	    {
	      TmpNIndices[k] = this->MaxMomentum - TmpNIndices[k];
	    }
	  else
	    {
	      TmpNIndices[k] = 0;
	    }
	}
      // warning do not set TmpNbrPermutation to zero
      SortArrayDownOrderingPermutationBubbleSort<int>(TmpMIndices, this->NBodyValue, TmpNbrPermutation);
      SortArrayDownOrderingPermutationBubbleSort<int>(TmpNIndices, this->NBodyValue, TmpNbrPermutation);
      TmpCanonicalTotalMomentum = 0;
      TmpCanonicalMIndices = 0;
      TmpCanonicalNIndices = 0;
      for (int k = this->NBodyValue - 1; k >= 0; --k)
	{
	  TmpCanonicalMIndices *= this->MaxMomentum;
	  TmpCanonicalMIndices += TmpMIndices[k];
	  TmpCanonicalNIndices *= this->MaxMomentum;
	  TmpCanonicalNIndices += TmpNIndices[k];
	  TmpCanonicalTotalMomentum += TmpNIndices[k];
	}
      TmpCanonicalTotalMomentum %= this->MaxMomentum;
      if (TmpCanonicalMIndices > TmpCanonicalNIndices)
	{
	  int Tmp = TmpCanonicalMIndices;
	  TmpCanonicalMIndices = TmpCanonicalNIndices;
	  TmpCanonicalNIndices = Tmp;
	}
      if ((TmpCanonicalTotalMomentum < canonicalTotalMomentum) || 
	  ((TmpCanonicalTotalMomentum == canonicalTotalMomentum) && ((TmpCanonicalNIndices < canonicalNIndices) ||
								     ((TmpCanonicalNIndices == canonicalNIndices) && (TmpCanonicalMIndices < canonicalMIndices)))))
	{
	  IsCanonical = false;
	  canonicalMIndices = TmpCanonicalMIndices;
	  canonicalNIndices = TmpCanonicalNIndices;
	  canonicalTotalMomentum = TmpCanonicalTotalMomentum;
	  if ((TmpNbrPermutation & 1) == 0)
	    {
	      canonicalPermutationCoefficient = 1.0;
	    }
	  else
	    {
	      canonicalPermutationCoefficient = -1.0;
	    }
	}
    }
  delete[] TmpMIndices;
  delete[] TmpNIndices;
  return IsCanonical;
}

