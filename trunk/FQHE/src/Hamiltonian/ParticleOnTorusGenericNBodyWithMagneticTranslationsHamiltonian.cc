////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          generic n-body interaction                        //
//                                                                            //
//                        last modification : 24/01/2015                      //
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
#include "Hamiltonian/ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/FQHETorusComputeMatrixElementOperation.h"

#include <iostream>
#include <algorithm>
#include <sys/time.h>


// default constructor
//

ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// maxMomentum = number of flux quanta
// xMomentum = relative angular momentum along x 
// ratio = torus aspect ratio (Lx/Ly)
// nbrNBody = type of interaction i.e. the number of density operators that are involved in the interaction
// interactionName = name of the interaction, will be use to generate the interaction matrix element output file name
// nbrMonomials = number of monomials in the Fourier transformed interaction
// monomialCoefficients = coefficients in front of each monomial in the Fourier transformed interaction
// monomialDescription = description of each monomial in the Fourier transformed interaction
// regenerateElementFlag = regenerate th interaction matrix elements instead of reading them from the harddrive
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio,
															       int nbrNBody, char* interactionName, 
															       int nbrMonomials, double* monomialCoefficients, int** monomialDescription, 
															       bool regenerateElementFlag, AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
															       char* precalculationFileName)
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
  this->OneBodyInteractionFactors = 0;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  this->EvaluateExponentialFactors();
  this->HamiltonianShift = 0.0;
  char* MatrixElementFileName = new char [strlen(interactionName) + 256];
  sprintf (MatrixElementFileName, "%dbody_%s_2s_%d_ratio_%.14f.dat", this->NBodyValue, interactionName, this->MaxMomentum, this->Ratio);
  this->NbrMonomials = nbrMonomials;
  this->MonomialCoefficients = monomialCoefficients;
  this->MonomialDescription = monomialDescription;
  if ((regenerateElementFlag == true) || (!(IsFile(MatrixElementFileName))))
    {
      this->EvaluateInteractionFactors();
      this->WriteInteractionFactors(MatrixElementFileName);
    }
  else
    {
      this->ReadInteractionFactors(MatrixElementFileName);
    }
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

ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::~ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian()
{
}
  
// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// evaluate all interaction factors
//   

void ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0l;
  this->GetIndices();
  double* QxValues = new double [this->NBodyValue];
  double* QyValues = new double [this->NBodyValue];
  double* Q2Values = new double [this->NBodyValue];
  double* CosineCoefficients = new double [this->NBodyValue];
  double Factor = pow((double) this->MaxMomentum, -0.5 * ((double) (this->NBodyValue)));
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
      int NbrPermutations = 1;
      for (int i = 1; i <= this->NBodyValue; ++i)
	NbrPermutations *= i;
      int** Permutations = new int*[NbrPermutations]; 
      double* PermutationSign = new double[NbrPermutations]; 
      Permutations[0] = new int [this->NBodyValue];
      for (int i = 0; i < this->NBodyValue; ++i)
	Permutations[0][i] = i;
      PermutationSign[0] = 1.0;
      double TmpSign = 1.0;
      for (int i = 1; i < NbrPermutations; ++i)
	{
	  Permutations[i] = new int [this->NBodyValue];
	  for (int j = 0; j < this->NBodyValue; ++j)
	    Permutations[i][j] = Permutations[i - 1][j];
	  int* TmpArrayPerm = Permutations[i];
	  int Pos1 = this->NBodyValue - 1;
	  while (TmpArrayPerm[Pos1 - 1] >= TmpArrayPerm[Pos1])
	    --Pos1;
	  --Pos1;
	  int Pos2 = this->NBodyValue - 1;      
	  while (TmpArrayPerm[Pos2] <= TmpArrayPerm[Pos1])
	    --Pos2;
	  int TmpIndex = TmpArrayPerm[Pos1];
	  TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	  TmpArrayPerm[Pos2] = TmpIndex;
	  TmpSign *= -1.0;
	  Pos2 = this->NBodyValue - 1;   
	  Pos1++;
	  while (Pos1 < Pos2)
	    {
	      TmpIndex = TmpArrayPerm[Pos1];
	      TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	      TmpArrayPerm[Pos2] = TmpIndex;
	      ++Pos1;
	      --Pos2;
	      TmpSign *= -1.0;
	    }
	  PermutationSign[i] = TmpSign;
	}


      int* TmpNIndices =  new int [this->NBodyValue];
      int* TmpMIndices =  new int [this->NBodyValue];

      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      int TmpCanonicalMIndices;
      int TmpCanonicalNIndices;
      int TmpCanonicalSum;
      double TmpCanonicalPermutationCoefficient;

      long TmpNbrMatrixElements = 0l;
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
		      TmpNIndices[k]  = this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
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
       FQHETorusComputeMatrixElementOperation Operation(this, TmpNbrMatrixElements, TmpMatrixElementIIndices, TmpMatrixElementJ1Indices, TmpMatrixElementJ2Indices, TmpMatrixElement);
      Operation.ApplyOperation(this->Architecture);
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
		      this->NBodyInteractionFactors[i][Index] = TmpMatrixElement[TmpNbrMatrixElements] * Factor;
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
      int* TmpNIndices =  new int [this->NBodyValue];
      int* TmpMIndices =  new int [this->NBodyValue];
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
		      TmpNIndices[k] = this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
		    }
		  double TmpInteraction = 0.0;
		  for (int k = 0; k < this->NBodyValue; ++k)
		    {
		      TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + k];
		    }
		  TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, 
									 QxValues, QyValues, Q2Values, CosineCoefficients);
		  while (std::prev_permutation(TmpMIndices, TmpMIndices  + this->NBodyValue))
		    {
		      TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, 
									     QxValues, QyValues, Q2Values, CosineCoefficients);		    
		    }
		  while (std::prev_permutation(TmpNIndices, TmpNIndices  + this->NBodyValue))
		    {
		      for (int k = 0; k < this->NBodyValue; ++k)
			{
			  TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + k];
			}
		      TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, 
									     QxValues, QyValues, Q2Values, CosineCoefficients);
		      while (std::prev_permutation(TmpMIndices, TmpMIndices  + this->NBodyValue))
			{
			  TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, 
										 QxValues, QyValues, Q2Values, CosineCoefficients);
			}
		    }
		  this->NBodyInteractionFactors[i][Index] = TmpInteraction * Factor;
		  for (int k = 0; k < this->NBodyValue; ++k)
		    cout << TmpMIndices[k] << " ";
		  cout << "| ";
		  for (int k = 0; k < this->NBodyValue; ++k)
		    cout << TmpNIndices[k] << " ";
		  TotalNbrInteractionFactors++;
		  ++Index;
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
  delete[] QxValues;
  delete[] QyValues;
  delete[] Q2Values;
  delete[] CosineCoefficients;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}
	
  
// evaluate the numerical coefficient  in front of the Prod a^+_mi Prod a+_n coupling term
//
// mIndices = array containing the creation operator indices
// nIndices = array containing the annihilation operator indices
// return value = numerical coefficient  

double ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int* mIndices, int* nIndices, double* qxValues, double* qyValues, double* q2Values, double* cosineCoefficients)
{
  int Tmp;
  double Prefactor = powl(this->MaxMomentum, -this->NBodyValue + 1.0);
  for (int i = 0; i < this->NBodyValue; ++i)
    {
      qxValues[i] = 0.0;
      qyValues[i] = (double) (nIndices[i] - mIndices[i]);
      cosineCoefficients[i] = 2.0 * ((double) mIndices[i]);
    }  
  double CurrentPrecision;
  double Coefficient = Prefactor * this->RecursiveEvaluateInteractionCoefficient(0.0, 0.0, 0.0, 0.0, CurrentPrecision, 
										 qxValues, qyValues, q2Values, cosineCoefficients);
  return Coefficient;
}
  
double ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::RecursiveEvaluateInteractionCoefficient(double currentSumQx, double currentSumQy, double currentSumQ2, double currentSumPhase, double& currentPrecision, double* qxValues, double* qyValues, double* q2Values, double* cosineCoefficients)
{
  double TotalCoefficient  = 0.0;
  double Coefficient  = 1.0;
  double CurrentQy = qyValues[0];
  currentPrecision = 0.0;
  double TmpPrecision;
  double Tmp;
  while ((fabs(Coefficient) + fabs(TotalCoefficient)) != fabs(TotalCoefficient))
    {	        
      Tmp = (qyValues[0] * qyValues[0] * this->Ratio);
      qxValues[0] = 0.0;
      q2Values[0] = Tmp;
      Coefficient = this->RecursiveEvaluateInteractionCoefficient(1, currentSumQx + qxValues[0], 
								  currentSumQy + qyValues[0], 
								  currentSumQ2 + q2Values[0], 
								  currentSumPhase +  qxValues[0] * (cosineCoefficients[0] + qyValues[0]), TmpPrecision, qxValues, qyValues, q2Values, cosineCoefficients);
      currentPrecision += TmpPrecision;
      if (Coefficient == 0.0)
	TmpPrecision = 1.0;
      while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
	{	  
	  ++qxValues[0];
	  q2Values[0] = (qxValues[0] * qxValues[0] * this->InvRatio) + Tmp;
	  Coefficient += 2.0 * this->RecursiveEvaluateInteractionCoefficient(1, currentSumQx + qxValues[0], 
									     currentSumQy + qyValues[0], 
									     currentSumQ2 + q2Values[0], 
									     currentSumPhase +  qxValues[0] * (cosineCoefficients[0] + qyValues[0]), TmpPrecision, qxValues, qyValues, q2Values, cosineCoefficients);
	  currentPrecision += 2.0 * TmpPrecision;
	}
      qyValues[0] +=  (double) this->MaxMomentum;
      TotalCoefficient += Coefficient;
    }
  qyValues[0] = CurrentQy -  (double) this->MaxMomentum;
  if (TotalCoefficient == 0.0)
    Coefficient = 1.0;
  else
    Coefficient = 2.0 * TotalCoefficient;
  while ((fabs(Coefficient) + fabs(TotalCoefficient)) != fabs(TotalCoefficient))
    {	        
      qxValues[0] = 0.0;
      Tmp = (qyValues[0] * qyValues[0] * this->Ratio);
      q2Values[0] = Tmp;
      Coefficient = this->RecursiveEvaluateInteractionCoefficient(1, currentSumQx + qxValues[0], 
								  currentSumQy + qyValues[0], 
								  currentSumQ2 + q2Values[0], 
								  currentSumPhase +  qxValues[0] * (cosineCoefficients[0] + qyValues[0]), TmpPrecision, qxValues, qyValues, q2Values, cosineCoefficients);
      currentPrecision += TmpPrecision;
      if (Coefficient == 0.0)
	TmpPrecision = 1.0;
      while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
	{	  
	  ++qxValues[0];
	  q2Values[0] = (qxValues[0] * qxValues[0] * this->InvRatio) + Tmp;
	  Coefficient += 2.0 * this->RecursiveEvaluateInteractionCoefficient(1, currentSumQx + qxValues[0], 
									     currentSumQy + qyValues[0], 
									     currentSumQ2 + q2Values[0], 
									     currentSumPhase +  qxValues[0] * (cosineCoefficients[0] + qyValues[0]), TmpPrecision, qxValues, qyValues, q2Values, cosineCoefficients);
	  currentPrecision += 2.0 * TmpPrecision;
	}
      qyValues[0] -=  (double) this->MaxMomentum;
      TotalCoefficient += Coefficient;
    }      
  qyValues[0] = CurrentQy;
  return TotalCoefficient;
}


double ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::RecursiveEvaluateInteractionCoefficient(int xPosition, double currentSumQx, double currentSumQy, double currentSumQ2, double currentSumPhase, double& currentPrecision, double* qxValues, double* qyValues, double* q2Values, double* cosineCoefficients)
{
  if (xPosition < (this->NBodyValue - 1))
    {
      double TotalCoefficient  = 0.0;
      double Coefficient  = 1.0;
      double CurrentQy = qyValues[xPosition];
      currentPrecision = 0.0;
      double TmpPrecision;
      double Tmp;
      while ((fabs(Coefficient) + fabs(TotalCoefficient)) != fabs(TotalCoefficient))
	{	        
	  Tmp = (qyValues[xPosition] * qyValues[xPosition] * this->Ratio);
	  qxValues[xPosition] = 0.0;
	  q2Values[xPosition] = Tmp;
	  Coefficient = this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + qxValues[xPosition], 
								      currentSumQy + qyValues[xPosition], 
								      currentSumQ2 + q2Values[xPosition], 
								      currentSumPhase +  qxValues[xPosition] * (cosineCoefficients[xPosition] + qyValues[xPosition]), TmpPrecision, qxValues, qyValues, q2Values, cosineCoefficients);
	  currentPrecision += TmpPrecision;
	  if (Coefficient == 0.0)
	    TmpPrecision = 1.0;
	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
	    {	  
	      ++qxValues[xPosition];
	      q2Values[xPosition] = (qxValues[xPosition] * qxValues[xPosition] * this->InvRatio) + Tmp;
	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + qxValues[xPosition], 
									   currentSumQy + qyValues[xPosition], 
									   currentSumQ2 + q2Values[xPosition], 
									   currentSumPhase +  qxValues[xPosition] * (cosineCoefficients[xPosition] + qyValues[xPosition]), TmpPrecision, qxValues, qyValues, q2Values, cosineCoefficients);
	      currentPrecision += TmpPrecision;
	    }
	  qxValues[xPosition] = 0.0;
	  if (Coefficient == 0.0)
	    TmpPrecision = 1.0;
	  else
	    TmpPrecision = 2.0 * fabs(Coefficient);
	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
	    {	  
	      --qxValues[xPosition];
	      q2Values[xPosition] = (qxValues[xPosition] * qxValues[xPosition] * this->InvRatio) + Tmp;
	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + qxValues[xPosition], 
									   currentSumQy + qyValues[xPosition], 
									   currentSumQ2 + q2Values[xPosition], 
									   currentSumPhase +  qxValues[xPosition] * (cosineCoefficients[xPosition] + qyValues[xPosition]), TmpPrecision, qxValues, qyValues, q2Values, cosineCoefficients);
	      currentPrecision += TmpPrecision;
	    }
	  qyValues[xPosition] +=  (double) this->MaxMomentum;
	  TotalCoefficient += Coefficient;
	}
      qyValues[xPosition] = CurrentQy -  (double) this->MaxMomentum;
      if (TotalCoefficient == 0.0)
	Coefficient = 1.0;
      else
	Coefficient = 2.0 * TotalCoefficient;
      while ((fabs(Coefficient) + fabs(TotalCoefficient)) != fabs(TotalCoefficient))
	{	        
	  qxValues[xPosition] = 0.0;
	  Tmp = (qyValues[xPosition] * qyValues[xPosition] * this->Ratio);
	  q2Values[xPosition] = Tmp;
	  Coefficient = this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + qxValues[xPosition], 
								      currentSumQy + qyValues[xPosition], 
								      currentSumQ2 + q2Values[xPosition], 
								      currentSumPhase +  qxValues[xPosition] * (cosineCoefficients[xPosition] + qyValues[xPosition]), TmpPrecision, qxValues, qyValues, q2Values, cosineCoefficients);
	  currentPrecision += TmpPrecision;
	  if (Coefficient == 0.0)
	    TmpPrecision = 1.0;
	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
	    {	  
	      ++qxValues[xPosition];
	      q2Values[xPosition] = (qxValues[xPosition] * qxValues[xPosition] * this->InvRatio) + Tmp;
	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + qxValues[xPosition], 
									   currentSumQy + qyValues[xPosition], 
									   currentSumQ2 + q2Values[xPosition], 
									   currentSumPhase +  qxValues[xPosition] * (cosineCoefficients[xPosition] + qyValues[xPosition]), TmpPrecision, qxValues, qyValues, q2Values, cosineCoefficients);
	      currentPrecision += TmpPrecision;
	    }
	  qxValues[xPosition] = 0.0;
	  if (Coefficient == 0.0)
	    TmpPrecision = 1.0;
	  else
	    TmpPrecision = 2.0 * fabs(Coefficient);
	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
	    {	  
	      --qxValues[xPosition];
	      q2Values[xPosition] = (qxValues[xPosition] * qxValues[xPosition] * this->InvRatio) + Tmp;
	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + qxValues[xPosition], 
									   currentSumQy + qyValues[xPosition], 
									   currentSumQ2 + q2Values[xPosition], 
									   currentSumPhase +  qxValues[xPosition] * (cosineCoefficients[xPosition] + qyValues[xPosition]), TmpPrecision, qxValues, qyValues, q2Values, cosineCoefficients);
	      currentPrecision += TmpPrecision;
	    }
	  qyValues[xPosition] -=  (double) this->MaxMomentum;
	  TotalCoefficient += Coefficient;
	}      
      qyValues[xPosition] = CurrentQy;
      return TotalCoefficient;
    }
  else
    {
      double TmpExponentialFactor = M_PI / ((double) this->MaxMomentum);
      qxValues[xPosition] = -currentSumQx;
      qyValues[xPosition] = -currentSumQy;  
      q2Values[xPosition] = (qxValues[xPosition] * qxValues[xPosition] * this->InvRatio) + (qyValues[xPosition] * qyValues[xPosition] * this->Ratio);
      currentSumPhase += qxValues[xPosition] * (cosineCoefficients[xPosition] + qyValues[xPosition]);
      currentPrecision = exp(- 0.5 * TmpExponentialFactor * (q2Values[xPosition] + currentSumQ2)) * this->VFactor(q2Values);
      return (cos(TmpExponentialFactor * currentSumPhase) * currentPrecision);
    }
}



double ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::VFactor(double* q2Values)
{
  double Tmp = 0.0;
  double Tmp2;
  for (int i = 0; i < this->NbrMonomials; ++i)
    {
      Tmp2 = this->MonomialCoefficients[i];
      for (int j = 0; j < this->NBodyValue; ++j)
	{
  	  for (int k = 0; k < this->MonomialDescription[i][j]; ++k)
  	    Tmp2 *= q2Values[j];
	}
      Tmp += Tmp2;
    }
  return Tmp;
}


// read the interaction matrix elements from disk
//
// fileName = name of the file where the interaction matrix elements are stored
// return value = true if no error occured

bool ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::ReadInteractionFactors(char* fileName)
{
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

  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "cannot open " << fileName << endl;
      return false;
    }
  this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
  for (int i = 0; i < this->NbrNBodySectorSums; ++i)
    {
      this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  int Index = j1 * this->NbrNBodySectorIndicesPerSum[i];
	  for (int j2 = 0; j2 <= j1; ++j2)
	    {
	      ReadLittleEndian(File, this->NBodyInteractionFactors[i][Index].Re);
	      ReadLittleEndian(File, this->NBodyInteractionFactors[i][Index].Im);
	      ++Index;
	    }
	}
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  for (int j2 = j1 + 1; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
	    {
	      this->NBodyInteractionFactors[i][j1 * this->NbrNBodySectorIndicesPerSum[i] + j2] = this->NBodyInteractionFactors[i][j2 * this->NbrNBodySectorIndicesPerSum[i] + j1];
	    }
	}
    }
  File.close();
  return true;
}

// write the interaction matrix elements from disk
//
// fileName = name of the file where the interaction matrix elements are stored
// return value = true if no error occured

bool ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::WriteInteractionFactors(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "cannot create " << fileName << endl;
      return false;
    }
  for (int i = 0; i < this->NbrNBodySectorSums; ++i)
    {
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  int Index = j1 * this->NbrNBodySectorIndicesPerSum[i];
	  for (int j2 = 0; j2 <= j1; ++j2)
	    {
	      WriteLittleEndian(File, this->NBodyInteractionFactors[i][Index].Re);
	      WriteLittleEndian(File, this->NBodyInteractionFactors[i][Index].Im);
	      ++Index;
	    }
	}
    }
  File.close();
  return true;
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

bool ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::FindCanonicalIndices(int* mIndices, int* nIndices, int linearizedMIndices, int linearizedNIndices, 
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

