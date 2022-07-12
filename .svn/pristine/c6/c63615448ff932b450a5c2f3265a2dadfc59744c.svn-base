////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of abstract quantum Hall hamiltonian associated            //
//     to particles on a sphere with spin and n-body interaction terms        //
//                                                                            //
//                        last modification : 26/08/2008                      //
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
#include "Hamiltonian/AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::ios;


// destructor
//

AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian::~AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian()
{
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
											int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
      if (this->FullTwoBodyFlag == true)
	for (int i = firstComponent; i < LastComponent; ++i)
	  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination) ;
      for (int k = 3; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	 {
	   for (int i = firstComponent; i < LastComponent; ++i)
	     this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination, k) ;
	 }	
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      Coefficient = vSource[k];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      vDestination[k++] += this->HamiltonianShift * Coefficient;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
	      int* TmpIndexArray;
	      double* TmpCoefficientArray; 
	      int j;
	      int TmpNbrInteraction;
	      firstComponent -= this->PrecalculationShift;
	      LastComponent -= this->PrecalculationShift;
	      int Pos = firstComponent / this->FastMultiplicationStep; 
	      int PosMod = firstComponent % this->FastMultiplicationStep;
	      if (PosMod != 0)
		{
		  ++Pos;
		  PosMod = this->FastMultiplicationStep - PosMod;
		}
	      int l =  PosMod + firstComponent + this->PrecalculationShift;
	      for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
		{
		  TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
		  TmpIndexArray = this->InteractionPerComponentIndex[Pos];
		  TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
		  Coefficient = vSource[l];
		  for (j = 0; j < TmpNbrInteraction; ++j)
		    vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		  vDestination[l] += this->HamiltonianShift * Coefficient;
		  l += this->FastMultiplicationStep;
		  ++Pos;
		}
	      firstComponent += this->PrecalculationShift;
	      LastComponent += this->PrecalculationShift;
	      for (l = 0; l < this->FastMultiplicationStep; ++l)
		if (PosMod != l)
		  {	
		    
		    if (this->FullTwoBodyFlag == true)
		      for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
		    for (int k = 2; k <= this->MaxNBody; ++k)
		      if (this->NBodyFlags[k] == true)
			{
			  for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			    this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination, k) ;
			}
		    this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSource, vDestination);
		  }
	      delete TmpParticles;
	    }
	  else
	    this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);	  
	}
    }
  if (this->L2Hamiltonian != 0)
    this->L2Hamiltonian->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  if (this->S2Hamiltonian != 0)
    this->S2Hamiltonian->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
												int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      double* Coefficient2 = new double [nbrVectors];
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
      if (this->FullTwoBodyFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2) ;
	}
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	  {
	    for (int i = firstComponent; i < LastComponent; ++i)
	      this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, k, Coefficient2);
	  }	    
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
      delete TmpParticles;
      delete[] Coefficient2;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  double* Coefficient2 = new double [nbrVectors];
	  int j;
	  int Pos;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  Coefficient2[l] = vSources[l][k];
		  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		    }
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
	      int* TmpIndexArray;
	      double* TmpCoefficientArray; 
	      double* Coefficient2 = new double [nbrVectors];
	      int j;
	      int Pos2;
	      int TmpNbrInteraction;
	      firstComponent -= this->PrecalculationShift;
	      LastComponent -= this->PrecalculationShift;
	      int Pos = firstComponent / this->FastMultiplicationStep; 
	      int PosMod = firstComponent % this->FastMultiplicationStep;
	      if (PosMod != 0)
		{
		  ++Pos;
		  PosMod = this->FastMultiplicationStep - PosMod;
		}
	      int l =  PosMod + firstComponent + this->PrecalculationShift;
	      for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
		{
		  TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
		  TmpIndexArray = this->InteractionPerComponentIndex[Pos];
		  TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      Coefficient2[k] = vSources[k][l];
		    }
		  for (j = 0; j < TmpNbrInteraction; ++j)
		    {
		      Pos2 = TmpIndexArray[j];
		      Coefficient = TmpCoefficientArray[j];
		      for (int k = 0; k < nbrVectors; ++k)
			{
			  vDestinations[k][Pos2] += Coefficient  * Coefficient2[k];
			}
		    }
		  l += this->FastMultiplicationStep;
		  ++Pos;
		}
	      firstComponent += this->PrecalculationShift;
	      LastComponent += this->PrecalculationShift;
	      for (l = 0; l < this->FastMultiplicationStep; ++l)
		if (PosMod != l)
		  {	
		    if (this->FullTwoBodyFlag == true)
		      {
			for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2) ;
 		      }
		    for (int k = 2; k <= this->MaxNBody; ++k)
		      if (this->NBodyFlags[k] == true)
			{
			  for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			    this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, k, Coefficient2) ;
			}
		    this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSources, vDestinations, nbrVectors);
		  }
	      delete[] Coefficient2;
	      delete TmpParticles;
	    }
	  else
	    this->LowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	}
    }
  if (this->L2Hamiltonian != 0)
    this->L2Hamiltonian->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
  if (this->S2Hamiltonian != 0)
    this->S2Hamiltonian->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
  return vDestinations;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long Memory = 0;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;
  
  if ((this->OneBodyInteractionFactorsdowndown != 0) || (this->OneBodyInteractionFactorsupup != 0))
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  ++Memory;
	  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
	}
    }
  
  if (this->FullTwoBodyFlag == true)
    EvaluateMNTwoBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);

  EvaluateMNOneBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
      
  for (int k = 3; k <= this->MaxNBody; ++k)
    if (this->NBodyFlags[k] == true)
      EvaluateMNNBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, k, Memory);
  delete TmpParticles;
  
  return Memory;
}


// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that have to be precalcualted
//
void AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{
  int LastComponent = nbrComponent + firstComponent;

  int* TmpIndexArray;
  double* TmpCoefficientArray;
  long ColumnIndex;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();

  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  long Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
      ColumnIndex = 0l;
      if (this->FullTwoBodyFlag == true)
	this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, ColumnIndex);
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	  {
	    this->EvaluateMNNBodyFastMultiplicationComponent(TmpParticles, i, k, TmpIndexArray, TmpCoefficientArray, ColumnIndex);
	  }
      if ((this->OneBodyInteractionFactorsdowndown != 0) || (this->OneBodyInteractionFactorsupup != 0))
	{
	  double TmpDiagonal = 0.0;
	  if (this->OneBodyInteractionFactorsupup != 0)
	    for (int j = 0; j <= this->LzMax; ++j) 
	      TmpDiagonal += this->OneBodyInteractionFactorsupup[j] * TmpParticles->AduAu(i + this->PrecalculationShift, j);
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	    for (int j = 0; j <= this->LzMax; ++j) 
	      TmpDiagonal += this->OneBodyInteractionFactorsdowndown[j] * TmpParticles->AddAd(i + this->PrecalculationShift, j);
	  TmpIndexArray[ColumnIndex] = i + this->PrecalculationShift;
	  TmpCoefficientArray[ColumnIndex] = TmpDiagonal;
	  ++ColumnIndex;	  
	}
      ++Pos;
    }
  delete TmpParticles;
}


// enable fast multiplication algorithm using on disk cache 
//
// fileName = prefix of the name of the file where temporary matrix elements will be stored

void AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
{
  if (this->FastMultiplicationStep == 1)
    {
      this->DiskStorageFlag = false;
      this->DiskStorageFileName = 0;
      this->EnableFastMultiplication();
      return;
    }

  this->DiskStorageFlag = true;
  this->DiskStorageFileName = new char [strlen(fileName) + 8];
  sprintf (this->DiskStorageFileName, "%s.ham", fileName);

  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  this->DiskStorageStart = (int) MinIndex;
  int DiskStorageEnd = 1 + (int) MaxIndex;

  int Index;
  double Coefficient;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  long Pos;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  this->InteractionPerComponentIndex = 0;
  this->InteractionPerComponentCoefficient = 0;
  double* TmpInteraction;
  int* TmpMIndices;
  int* TmpNIndices;
  int m1;
  int m2;
  int m3;
  this->MaxNbrInteractionPerComponent = 0;

  long TotalPos = 0l;
  ofstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::out);
 
  File.write((char*) &(EffectiveHilbertSpaceDimension), sizeof(int));
  File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
  File.write((char*) this->NbrInteractionPerComponent, sizeof(int) * EffectiveHilbertSpaceDimension);

  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      FileJump += (long) this->NbrInteractionPerComponent[i];
      if (this->MaxNbrInteractionPerComponent < this->NbrInteractionPerComponent[i])
	this->MaxNbrInteractionPerComponent = this->NbrInteractionPerComponent[i];
    }
  FileJump *= sizeof(int);

  TmpIndexArray = new int [this->MaxNbrInteractionPerComponent];
  TmpCoefficientArray = new double [this->MaxNbrInteractionPerComponent];      

  for (int i = this->DiskStorageStart; i < DiskStorageEnd; ++i)
    {
      if (this->NbrInteractionPerComponent[TotalPos] > 0)
	{
	  Pos = 0;
	  if (this->OneBodyTermFlag == true)
	    {
	      for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
		{
		  int m1 = this->OneBodyMValues[j];
		  int m2 = this->OneBodyNValues[j];
		  Index = this->Particles->AdA(i, m1, m2, Coefficient);
		  if (Index < this->Particles->GetHilbertSpaceDimension())
		    {
		      TmpIndexArray[Pos] = Index;
		      TmpCoefficientArray[Pos] = Coefficient * this->OneBodyInteractionFactors[j];
		      ++Pos;
		    }
		}
	    }
	  if (this->NBodyFlags[1] == true)
	    {
	      double Sign = this->NBodySign[1][0];
	      int TmpMinSumIndices = this->MinSumIndices[1][0];
	      int TmpMaxSumIndices = this->MaxSumIndices[1][0];	      
	      for (int j = TmpMinSumIndices; j <= TmpMaxSumIndices; ++j)
		{
		  int Lim = NbrSortedIndicesPerSum[1][0][j];
		  TmpNIndices = this->SortedIndicesPerSum[1][0][j];
		  TmpInteraction = this->NBodyInteractionFactors[1][0][j];
		  for (int i1 = 0; i1 < Lim; ++i1)
		    {
		      TmpMIndices = this->SortedIndicesPerSum[1][0][j];
		      for (int i2 = 0; i2 < Lim; ++i2)
			{
			  Index = this->Particles->AdA(i, TmpMIndices[i1], TmpNIndices[i2], Coefficient);
			  if (Index < this->Particles->GetHilbertSpaceDimension())
			    {
			      TmpIndexArray[Pos] = Index;
			      TmpCoefficientArray[Pos] = Sign * Coefficient * TmpInteraction[i1] *  TmpInteraction[i2];
			      ++Pos;
			    }
			  ++TmpMIndices;
			}
		    }
		  ++TmpNIndices;
		}
	    }
	  if (this->FullTwoBodyFlag == true)
	    {
	      if (this->NbrM12Indices == 0)
		{
		  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
		    {
		      m1 = this->M1Value[j];
		      m2 = this->M2Value[j];
		      m3 = this->M3Value[j];
		      Index = this->Particles->AdAdAA(i, m1, m2, m3, m1 + m2 - m3, Coefficient);
		      if (Index < this->Particles->GetHilbertSpaceDimension())
			{
			  TmpIndexArray[Pos] = Index;
			  TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
			  ++Pos;
			}
		    }
		}
	      else
		{
		  double Coefficient2;
		  int SumIndices;
		  int TmpNbrM3Values;
		  int* TmpM3Values;
		  int ReducedNbrInteractionFactors = 0;
		  for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		    {
		      Coefficient = this->Particles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		      if (Coefficient != 0.0)
			{
			  SumIndices = this->M1Value[m1] + this->M2Value[m1];
			  TmpM3Values = this->M3Values[m1];
			  TmpNbrM3Values = this->NbrM3Values[m1];
			  for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			    {
			      Index = this->Particles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			      if (Index < this->Particles->GetHilbertSpaceDimension())
				{
				  TmpIndexArray[Pos] = Index;
				  TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
				  ++Pos;
				}
			      ++ReducedNbrInteractionFactors;
			    }    
			}
		      else
			ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		    }
		}
	    }
	  for (int k = 2; k <= this->MaxNBody; ++k)
	    if (this->NBodyFlags[k] == true)
	      this->EvaluateMNNBodyFastMultiplicationComponent(TmpParticles, i, k, TmpIndexArray, TmpCoefficientArray, Pos);
	  File.write((char*) TmpIndexArray, sizeof(int) * this->NbrInteractionPerComponent[TotalPos]);
	  FileJump -= sizeof(int) * this->NbrInteractionPerComponent[TotalPos];
	  File.seekp(FileJump, ios::cur);
	  File.write((char*) TmpCoefficientArray, sizeof(double) * this->NbrInteractionPerComponent[TotalPos]);
	  FileJump += sizeof(double) * this->NbrInteractionPerComponent[TotalPos];
	  File.seekp(-FileJump, ios::cur);	  
	}
      ++TotalPos;
    }
  delete[] TmpIndexArray;
  delete[] TmpCoefficientArray;
  File.close();

  this->FastMultiplicationFlag = true;
  this->BufferSize = this->Memory / ((this->MaxNbrInteractionPerComponent * (sizeof(int) + sizeof(double))) + sizeof(int*) + sizeof(double*));

  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// get all indices needed to characterize a completly skew symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long  AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian::GetAllSkewSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
											  int**& sortedIndicesPerSum)
{
  long** BinomialCoefficients = GetBinomialCoefficients(nbrValues);
  long NbrElements = BinomialCoefficients[nbrValues][nbrIndices];
  int** Indices = new int* [NbrElements];
  int* Sum = new int [NbrElements];
  int Min = nbrIndices - 1;
  int Max;
  int Step;
  int Pos = 0;
  for (int i = nbrValues - 1; i >= Min; --i)
    {
      Step = BinomialCoefficients[i][nbrIndices - 1];
      for (int j = 0; j < Step; ++j)
	{
	  Indices[Pos] = new int [nbrIndices];
	  Indices[Pos][0] = i;
	  Sum[Pos] = i;
	  ++Pos;
	}
    }
  for (int i = 1; i < nbrIndices; ++i)
    {
      int Pos = 0;
      Min = nbrIndices - i - 1;
      while (Pos < NbrElements)
	{
	  Max = Indices[Pos][i - 1] - 1;
	  for (; Max >= Min; --Max)
	    {
	      Step = BinomialCoefficients[Max][Min];
	      for (int j = 0; j < Step; ++j)
		{
		  Indices[Pos][i] = Max;
		  Sum[Pos] += Max;
		  ++Pos;
		}
	    }
	}
    }

  int MaxSum = ((nbrValues - 1) * nbrIndices) - (((nbrIndices - 1) * (nbrIndices)) / 2);
  int MinSum = (nbrIndices * (nbrIndices - 1)) / 2;
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = 0; i < NbrElements; ++i)
    ++nbrSortedIndicesPerSum[Sum[i]];
  long* TmpPos = new long [MaxSum + 1];
  for (int i = MinSum; i <= MaxSum; ++i)
    {
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * nbrIndices];
      nbrSortedIndicesPerSum[i] = 0;
      TmpPos[i] = 0l;      
    }
  for (int i = 0; i < NbrElements; ++i)
    {   
      Pos = Sum[i];
      Max = nbrSortedIndicesPerSum[Pos];
      for (int j = 0; j < nbrIndices; ++j)
	{
	  sortedIndicesPerSum[Pos][TmpPos[Pos]] = Indices[i][j];
	  ++TmpPos[Pos];
	}
      ++nbrSortedIndicesPerSum[Pos];
      delete[] Indices[i];
    }
  delete[] TmpPos;
  delete[] Sum;
  delete[]Indices;
  for (int i = 0; i <= nbrValues; ++i)
    delete[] BinomialCoefficients[i];
  delete[] BinomialCoefficients;
  return NbrElements;
}

// get all indices sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian::GetAllIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
									    int**& sortedIndicesPerSum)
{
  long NbrElements = nbrValues;
  for (int i = 1; i < nbrIndices; ++i)
    NbrElements *= (long) nbrValues;
  int** Indices = new int* [NbrElements];
  int* Sum = new int [NbrElements];
  int Max;
  int Step;
  int Pos = 0;
  int TmpNbrIndices;
  for (long i = 0l; i < NbrElements; ++i)
    {
      Indices[Pos] = new int [nbrIndices];
      Sum[Pos] = 0;
      long TmpIndex = i;
      for (int j = 0; j < nbrIndices; ++j)
	{	  
	  Indices[Pos][j] = (int) (TmpIndex % nbrValues);
	  TmpIndex /= (long) nbrValues;
	  Sum[Pos] += Indices[Pos][j];
	}
      ++Pos;
    }
  int MaxSum = (nbrValues - 1) * nbrIndices;
  long* TmpPos = new long [MaxSum + 1];
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = 0; i < NbrElements; ++i)
    {
      ++nbrSortedIndicesPerSum[Sum[i]];
    }
  for (int i = 0; i <= MaxSum; ++i)
    {
      TmpPos[i] = 0l;
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * nbrIndices];
      nbrSortedIndicesPerSum[i] = 0;
    }
  for (long i = 0l; i < NbrElements; ++i)
    {   
      Pos = Sum[i];
      Max = nbrSortedIndicesPerSum[Pos];
      int* TmpIndices = Indices[i];
      for (int j = 0; j < nbrIndices; ++j)
	{
	  sortedIndicesPerSum[Pos][TmpPos[Pos]] = TmpIndices[j];
	  ++TmpPos[Pos];
	}
      delete[] TmpIndices;
      ++nbrSortedIndicesPerSum[Pos];
    }

  delete[] TmpPos;
  delete[] Sum;
  delete[] Indices;
  return NbrElements;
}

// get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
//                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian::GetAllSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
										     int**& sortedIndicesPerSum,
										     double**& sortedIndicesPerSumSymmetryFactor)
{
  long** DimensionSymmetricGroup;
  if (nbrValues >= nbrIndices)
    DimensionSymmetricGroup = GetIrreducibleRepresentationDimensionSymmetricGroup(nbrValues);
  else
    DimensionSymmetricGroup = GetIrreducibleRepresentationDimensionSymmetricGroup(nbrIndices);
  long NbrElements = DimensionSymmetricGroup[nbrValues][nbrIndices];

  int** Indices = new int* [NbrElements];
  int* Sum = new int [NbrElements];
  int Max;
  int Step;
  int Pos = 0;
  int TmpNbrIndices;
  for (int i = nbrValues - 1; i >= 0; --i)
    {
      Step = DimensionSymmetricGroup[i + 1][nbrIndices - 1];
      for (int j = 0; j < Step; ++j)
	{
	  Indices[Pos] = new int [nbrIndices];
	  Indices[Pos][0] = i;
	  Sum[Pos] = i;
	  ++Pos;
	}
    }
  for (int i = 1; i < nbrIndices; ++i)
    {
      int Pos = 0;
      TmpNbrIndices = nbrIndices - i - 1;
      while (Pos < NbrElements)
	{
	  Max = Indices[Pos][i - 1];
	  Step = DimensionSymmetricGroup[Max + 1][TmpNbrIndices];
	  for (int j = 0; j < Step; ++j)
	    {
	      Indices[Pos][i] = Max;
	      Sum[Pos] += Max;
	      ++Pos;
	    }
	  --Max;
	  for (; Max >= 0; --Max)
	    {
	      Step = DimensionSymmetricGroup[Max + 1][TmpNbrIndices];
	      for (int j = 0; j < Step; ++j)
		{
		  Indices[Pos][i] = Max;
		  Sum[Pos] += Max;
		  ++Pos;
		}
	    }
	}
    }

  int MaxSum = (nbrValues - 1) * nbrIndices;
  long* TmpPos = new long [MaxSum + 1];
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  sortedIndicesPerSumSymmetryFactor = new double* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = 0; i < NbrElements; ++i)
    {
      ++nbrSortedIndicesPerSum[Sum[i]];
    }
  for (int i = 0; i <= MaxSum; ++i)
    {
      TmpPos[i] = 0l;
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * nbrIndices];
      sortedIndicesPerSumSymmetryFactor[i] = new double [nbrSortedIndicesPerSum[i]];
      nbrSortedIndicesPerSum[i] = 0;
    }
  for (int i = 0; i < NbrElements; ++i)
    {   
      Pos = Sum[i];
      Max = nbrSortedIndicesPerSum[Pos];
      int* TmpIndices = Indices[i];
      for (int j = 0; j < nbrIndices; ++j)
	{
	  sortedIndicesPerSum[Pos][TmpPos[Pos]] = TmpIndices[j];
	  ++TmpPos[Pos];
	}
      double& SymmetryFactor = sortedIndicesPerSumSymmetryFactor[Pos][Max];
      SymmetryFactor = 1.0;
      for (int j = 1; j < nbrIndices; ++j)
	{
	  int TmpSymmetryFactor = 1;
	  while ((j < nbrIndices) && (TmpIndices[j - 1] == TmpIndices[j]))
	    {
	      ++TmpSymmetryFactor;
	      ++j;
	    }
	  if (TmpSymmetryFactor != 1)
	    for (int k = 2; k <= TmpSymmetryFactor; ++k)
	      SymmetryFactor *= (double) k;
	}
      delete[] TmpIndices;
      SymmetryFactor = 1.0 / SymmetryFactor;
      ++nbrSortedIndicesPerSum[Pos];
    }

  delete[] TmpPos;
  delete[] Sum;
  delete[]Indices;
  for (int i = 0; i <= nbrValues; ++i)
    delete[] DimensionSymmetricGroup[i];
  delete[] DimensionSymmetricGroup;
  return NbrElements;
}

// get all indices needed to characterize a  tensor made of two completly skew symmetric  sets of indices, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndicesUp = number of indices for the first set of indices (i.e. spin up)
// nbrIndicesDown = number of indices for the first set of indices (i.e. spin down), warning nbrIndicesDown should lower of equal to nbrIndicesUp
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long  AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian::GetAllTwoSetSkewSymmetricIndices (int nbrValues, int nbrIndicesUp, int nbrIndicesDown, int*& nbrSortedIndicesPerSum, 
												int**& sortedIndicesPerSum)
{
  int* TmpNbrSortedIndicesPerSumUp;
  int** TmpSortedIndicesPerSumUp;
  this->GetAllSkewSymmetricIndices(nbrValues, nbrIndicesUp, TmpNbrSortedIndicesPerSumUp, TmpSortedIndicesPerSumUp);
  int MaxSumUp = ((nbrValues - 1) * nbrIndicesUp) - (((nbrIndicesUp - 1) * nbrIndicesUp) / 2);
  int MinSumUp = (nbrIndicesUp * (nbrIndicesUp - 1)) / 2;
  long NbrElements = 0l;
  if (nbrIndicesDown == 1)
    {
      int MaxSum = MaxSumUp + nbrValues - 1;
      int MinSum = MinSumUp;
      nbrSortedIndicesPerSum = new int [MaxSum + 1];
      sortedIndicesPerSum = new int* [MaxSum + 1];
      for (int i = 0; i <= MaxSum; ++i)
	nbrSortedIndicesPerSum[i] = 0;
      nbrSortedIndicesPerSum[MinSum] = TmpNbrSortedIndicesPerSumUp[MinSumUp];
      sortedIndicesPerSum[MinSum] = new int [nbrSortedIndicesPerSum[MinSum] * (nbrIndicesUp + 1)];     
      int Lim = nbrSortedIndicesPerSum[MinSum];
      int* TmpSortedIndicesPerSum = sortedIndicesPerSum[MinSum];
      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[MinSum];
      int Pos = 0;
      int Pos2 = 0;
      for (int i = 0; i < Lim; ++i)
	{
	  for (int j = 0; j < nbrIndicesUp; ++j)
	    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
	  TmpSortedIndicesPerSum[Pos2++] = 0;
	}
      nbrSortedIndicesPerSum[MaxSum] = TmpNbrSortedIndicesPerSumUp[MaxSumUp];
      sortedIndicesPerSum[MaxSum] = new int [nbrSortedIndicesPerSum[MaxSum] * (nbrIndicesUp + 1)];     
      Lim = nbrSortedIndicesPerSum[MaxSum];
      TmpSortedIndicesPerSum = sortedIndicesPerSum[MaxSum];
      TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[MaxSumUp];
      Pos = 0;
      Pos2 = 0;
      for (int i = 0; i < Lim; ++i)
	{
	  for (int j = 0; j < nbrIndicesUp; ++j)
	    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
	  TmpSortedIndicesPerSum[Pos2++] = nbrValues - 1;
	}
      for (int i = MinSum + 1; i < MaxSum; ++i)
	{
	  int TmpMinSumUp = i - nbrValues + 1;
	  if (TmpMinSumUp < MinSumUp)
	    TmpMinSumUp = MinSumUp;
	  int TmpMaxSumUp = i;
	  if (TmpMaxSumUp > MaxSumUp)
	    TmpMaxSumUp = MaxSumUp;
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    nbrSortedIndicesPerSum[i] += TmpNbrSortedIndicesPerSumUp[j];
	  NbrElements += (long) nbrSortedIndicesPerSum[i];
	  sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * (nbrIndicesUp + 1)];
	  int Pos2 = 0;
	  int* TmpSortedIndicesPerSum = sortedIndicesPerSum[i];
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    {
	      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[j];
	      int Lim = TmpNbrSortedIndicesPerSumUp[j];
	      int Pos = 0;
	      for (int k = 0; k < Lim; ++k)
		{
		  for (int l = 0; l < nbrIndicesUp; ++l)
		    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
		  TmpSortedIndicesPerSum[Pos2++] = i - j;
		}
	    }
	  
	}
    }
  else
    {
      int* TmpNbrSortedIndicesPerSumDown;
      int** TmpSortedIndicesPerSumDown;
      this->GetAllSkewSymmetricIndices(nbrValues, nbrIndicesDown, TmpNbrSortedIndicesPerSumDown, TmpSortedIndicesPerSumDown);
      int MaxSumDown = ((nbrValues - 1) * nbrIndicesDown) - (((nbrIndicesDown - 1) * nbrIndicesDown) / 2);
      int MinSumDown = (nbrIndicesDown * (nbrIndicesDown - 1)) / 2;
      int MaxSum = MaxSumUp + MaxSumDown;
      int MinSum = MinSumUp + MinSumDown;
      sortedIndicesPerSum = new int* [MaxSum + 1];
      for (int i = 0; i <= MaxSum; ++i)
	nbrSortedIndicesPerSum[i] = 0;
      for (int i = MinSum; i <= MaxSum; ++i)
	{
	  int TmpMinSumUp = i - MaxSumDown;
	  if (TmpMinSumUp < MinSumUp)
	    TmpMinSumUp = MinSumUp;
	  int TmpMaxSumUp = i - MinSumDown;
	  if (TmpMaxSumUp > MaxSumUp)
	    TmpMaxSumUp = MaxSumUp;
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    nbrSortedIndicesPerSum[i] += TmpNbrSortedIndicesPerSumUp[j] * TmpNbrSortedIndicesPerSumDown[i - j];
	  NbrElements += (long) nbrSortedIndicesPerSum[i];
	  sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * (nbrIndicesUp + nbrIndicesDown)];
	  int Pos2 = 0;
	  int* TmpSortedIndicesPerSum = sortedIndicesPerSum[i];
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    {
	      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[j];
	      int Lim = TmpNbrSortedIndicesPerSumUp[j];
	      int Pos = 0;
	      for (int k = 0; k < Lim; ++k)
		{
		  int Pos3 = 0;
		  int* TmpSortedIndicesPerSumDown2 = TmpSortedIndicesPerSumDown[i - j];
		  int Lim2 = TmpNbrSortedIndicesPerSumDown[i - j];
		  for (int m = 0; m < Lim2; ++m)
		    {
		      for (int l = 0; l < nbrIndicesUp; ++l)
			TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos + l];
		      for (int l = 0; l < nbrIndicesDown; ++l)
			TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumDown2[Pos3++];		      
		    }
		  Pos += nbrIndicesUp;
		}
	    }
	}
      for (long i = MinSumDown; i <= MaxSumDown; ++i)
	delete[] TmpSortedIndicesPerSumDown[i];
      delete[] TmpNbrSortedIndicesPerSumDown;
      delete[] TmpSortedIndicesPerSumDown;
    }
  for (long i = MinSumUp; i <= MaxSumUp; ++i)
    delete[] TmpSortedIndicesPerSumUp[i];
  delete[] TmpNbrSortedIndicesPerSumUp;
  delete[] TmpSortedIndicesPerSumUp;
  return NbrElements;
}

// get all indices needed to characterize a  tensor made of two completly symmetric  sets of indices, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndicesUp = number of indices for the first set of indices (i.e. spin up)
// nbrIndicesDown = number of indices for the first set of indices (i.e. spin down), warning nbrIndicesDown should lower of equal to nbrIndicesUp
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
//                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian::GetAllTwoSetSymmetricIndices (int nbrValues, int nbrIndicesUp, int nbrIndicesDown, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum,
											   double**& sortedIndicesPerSumSymmetryFactor)
{
  int* TmpNbrSortedIndicesPerSumUp;
  int** TmpSortedIndicesPerSumUp;
  double** TmpSortedIndicesPerSumSymmetryFactorUp;
  this->GetAllSymmetricIndices(nbrValues, nbrIndicesUp, TmpNbrSortedIndicesPerSumUp, TmpSortedIndicesPerSumUp,
			       TmpSortedIndicesPerSumSymmetryFactorUp);
  int MaxSumUp = (nbrValues - 1) * nbrIndicesUp;
  int MinSumUp = 0;
  long NbrElements = 0l;
  if (nbrIndicesDown == 1)
    {
      int MaxSum = MaxSumUp + nbrValues - 1;
      int MinSum = MinSumUp;
      nbrSortedIndicesPerSum = new int [MaxSum + 1];
      sortedIndicesPerSum = new int* [MaxSum + 1];
      sortedIndicesPerSumSymmetryFactor = new double* [MaxSum + 1];
      for (int i = 0; i <= MaxSum; ++i)
	nbrSortedIndicesPerSum[i] = 0;
      nbrSortedIndicesPerSum[MinSum] = TmpNbrSortedIndicesPerSumUp[MinSumUp];
      sortedIndicesPerSum[MinSum] = new int [nbrSortedIndicesPerSum[MinSum] * (nbrIndicesUp + 1)];     
      sortedIndicesPerSumSymmetryFactor[MinSum] = new double [nbrSortedIndicesPerSum[MinSum]];
      int Lim = nbrSortedIndicesPerSum[MinSum];
      int* TmpSortedIndicesPerSum = sortedIndicesPerSum[MinSum];
      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[MinSum];
      double* TmpSortedIndicesPerSumSymmetryFactor = sortedIndicesPerSumSymmetryFactor[MinSum];
      double* TmpSortedIndicesPerSumSymmetryFactorUp2 = TmpSortedIndicesPerSumSymmetryFactorUp[MinSum];
	int Pos = 0;
      int Pos2 = 0;
      for (int i = 0; i < Lim; ++i)
	{
	  for (int j = 0; j < nbrIndicesUp; ++j)
	    {
	      TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
	    }
	  TmpSortedIndicesPerSum[Pos2++] = 0;
	  TmpSortedIndicesPerSumSymmetryFactor[i] = TmpSortedIndicesPerSumSymmetryFactorUp2[i];
	}
      nbrSortedIndicesPerSum[MaxSum] = TmpNbrSortedIndicesPerSumUp[MaxSumUp];
      sortedIndicesPerSum[MaxSum] = new int [nbrSortedIndicesPerSum[MaxSum] * (nbrIndicesUp + 1)];     
      sortedIndicesPerSumSymmetryFactor[MaxSum] = new double [nbrSortedIndicesPerSum[MaxSum]];
      Lim = nbrSortedIndicesPerSum[MaxSum];
      TmpSortedIndicesPerSum = sortedIndicesPerSum[MaxSum];
      TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[MaxSumUp];
      TmpSortedIndicesPerSumSymmetryFactor = sortedIndicesPerSumSymmetryFactor[MaxSum];
      TmpSortedIndicesPerSumSymmetryFactorUp2 = TmpSortedIndicesPerSumSymmetryFactorUp[MaxSumUp];
      Pos = 0;
      Pos2 = 0;
      for (int i = 0; i < Lim; ++i)
	{
	  for (int j = 0; j < nbrIndicesUp; ++j)
	    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
	  TmpSortedIndicesPerSum[Pos2++] = nbrValues - 1;
	  TmpSortedIndicesPerSumSymmetryFactor[i] = TmpSortedIndicesPerSumSymmetryFactorUp2[i];
	}
      for (int i = MinSum + 1; i < MaxSum; ++i)
	{
	  int TmpMinSumUp = i - nbrValues + 1;
	  if (TmpMinSumUp < MinSumUp)
	    TmpMinSumUp = MinSumUp;
	  int TmpMaxSumUp = i;
	  if (TmpMaxSumUp > MaxSumUp)
	    TmpMaxSumUp = MaxSumUp;
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    nbrSortedIndicesPerSum[i] += TmpNbrSortedIndicesPerSumUp[j];
	  NbrElements += (long) nbrSortedIndicesPerSum[i];
	  sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * (nbrIndicesUp + 1)];
	  sortedIndicesPerSumSymmetryFactor[i] = new double [nbrSortedIndicesPerSum[i]];
	  int Pos2 = 0;
	  int* TmpSortedIndicesPerSum = sortedIndicesPerSum[i];
	  double* TmpSortedIndicesPerSumSymmetryFactor = sortedIndicesPerSumSymmetryFactor[i];
	  int Pos3 = 0;
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    {
	      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[j];
	      double* TmpSortedIndicesPerSumSymmetryFactorUp2 = TmpSortedIndicesPerSumSymmetryFactorUp[j];
	      int Lim = TmpNbrSortedIndicesPerSumUp[j];
	      int Pos = 0;
	      for (int k = 0; k < Lim; ++k)
		{
		  for (int l = 0; l < nbrIndicesUp; ++l)
		    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
		  TmpSortedIndicesPerSum[Pos2++] = i - j;
		  TmpSortedIndicesPerSumSymmetryFactor[Pos3++] = TmpSortedIndicesPerSumSymmetryFactorUp2[k];
		}
	    }
	}
    }
  else
    {
      int* TmpNbrSortedIndicesPerSumDown;
      int** TmpSortedIndicesPerSumDown;
      double** TmpSortedIndicesPerSumSymmetryFactorDown;
      this->GetAllSymmetricIndices(nbrValues, nbrIndicesDown, TmpNbrSortedIndicesPerSumDown, TmpSortedIndicesPerSumDown,
				   TmpSortedIndicesPerSumSymmetryFactorDown);
      int MaxSumDown = (nbrValues - 1) * nbrIndicesDown;
      int MinSumDown = 0;
      int MaxSum = MaxSumUp + MaxSumDown;
      int MinSum = MinSumUp + MinSumDown;
      sortedIndicesPerSum = new int* [MaxSum + 1];
      sortedIndicesPerSumSymmetryFactor = new double* [MaxSum + 1];
      for (int i = 0; i <= MaxSum; ++i)
	nbrSortedIndicesPerSum[i] = 0;
      for (int i = MinSum; i <= MaxSum; ++i)
	{
	  int TmpMinSumUp = i - MaxSumDown;
	  if (TmpMinSumUp < MinSumUp)
	    TmpMinSumUp = MinSumUp;
	  int TmpMaxSumUp = i - MinSumDown;
	  if (TmpMaxSumUp > MaxSumUp)
	    TmpMaxSumUp = MaxSumUp;
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    nbrSortedIndicesPerSum[i] += TmpNbrSortedIndicesPerSumUp[j] * TmpNbrSortedIndicesPerSumDown[i - j];
	  NbrElements += (long) nbrSortedIndicesPerSum[i];
	  sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * (nbrIndicesUp + nbrIndicesDown)];
	  sortedIndicesPerSumSymmetryFactor[i] = new double [nbrSortedIndicesPerSum[i]];
	  int Pos2 = 0;
	  int* TmpSortedIndicesPerSum = sortedIndicesPerSum[i];
	  double* TmpSortedIndicesPerSumSymmetryFactor = sortedIndicesPerSumSymmetryFactor[i];
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    {
	      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[j];
	      double* TmpSortedIndicesPerSumSymmetryFactorUp2 = TmpSortedIndicesPerSumSymmetryFactorUp[j];
	      int Lim = TmpNbrSortedIndicesPerSumUp[j];
	      int Pos = 0;
	      int Pos4 = 0;
	      for (int k = 0; k < Lim; ++k)
		{
		  int Pos3 = 0;
		  int* TmpSortedIndicesPerSumDown2 = TmpSortedIndicesPerSumDown[i - j];
		  double* TmpSortedIndicesPerSumSymmetryFactorDown2 = TmpSortedIndicesPerSumSymmetryFactorDown[i - j];
		  int Lim2 = TmpNbrSortedIndicesPerSumDown[i - j];
		  for (int m = 0; m < Lim2; ++m)
		    {
		      for (int l = 0; l < nbrIndicesUp; ++l)
			TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos + l];
		      for (int l = 0; l < nbrIndicesDown; ++l)
			TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumDown2[Pos3++];		      
		      TmpSortedIndicesPerSumSymmetryFactor[Pos4] = (TmpSortedIndicesPerSumSymmetryFactorUp2[k] * 
								    TmpSortedIndicesPerSumSymmetryFactorDown2[m]);
		      ++Pos4;
		    }
		  Pos += nbrIndicesUp;
		}
	    }
	}
      for (long i = MinSumDown; i <= MaxSumDown; ++i)
	delete[] TmpSortedIndicesPerSumDown[i];
      delete[] TmpNbrSortedIndicesPerSumDown;
      delete[] TmpSortedIndicesPerSumDown;
    }
  for (long i = MinSumUp; i <= MaxSumUp; ++i)
    delete[] TmpSortedIndicesPerSumUp[i];
  delete[] TmpNbrSortedIndicesPerSumUp;
  delete[] TmpSortedIndicesPerSumUp;
  return NbrElements;
}
