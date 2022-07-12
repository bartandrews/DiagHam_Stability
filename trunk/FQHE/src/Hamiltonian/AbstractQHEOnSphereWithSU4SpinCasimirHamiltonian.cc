////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract fractional quantum Hall hamiltonian          //
//              associated to particles with SU(4) spin on a sphere           //
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


#include "config.h"
#include "Hamiltonian/AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSU4SpinL2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSU4SpinS2Hamiltonian.h"
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


AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian()
{
  this->L2Hamiltonian=0;
  this->S2Hamiltonian=0;
}

// destructor
//

AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::~AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian()
{
      if ((this->M1IntraValue != 0) || (this->M1InterValue != 0))
      {
	delete[] this->M1IntraValue;
	delete[] this->M2IntraValue;
	for (int i = 0; i < this->NbrM12IntraIndices; ++i)
	  delete[] this->M3IntraValues[i];
	delete[] this->M3IntraValues;
	delete[] this->NbrM3IntraValues;
	
	delete[] this->M1InterValue;
	delete[] this->M2InterValue;
	for (int i = 0; i < this->NbrM12InterIndices; ++i)
	  delete[] this->M3InterValues[i];
	delete[] this->M3InterValues;
	delete[] this->NbrM3InterValues;
	if (M12InteractionFactorsupup)
	  delete[] M12InteractionFactorsupup;
	if (M12InteractionFactorsumum)
	  delete[] M12InteractionFactorsumum;
	delete[] M12InteractionFactorsdpdp;
	if (M12InteractionFactorsdpdp)
	  delete[] M12InteractionFactorsdmdm;
	delete[] M12InteractionFactorsupum;
	if (M12InteractionFactorsupum)
	  delete[] M12InteractionFactorsupdp;
	if (M12InteractionFactorsupdm)
	delete[] M12InteractionFactorsupdm;
	if (M12InteractionFactorsumdp)
	  delete[] M12InteractionFactorsumdp;
	if (M12InteractionFactorsumdm)
	  delete[] M12InteractionFactorsumdm;
	if (M12InteractionFactorsdpdm)
	  delete[] M12InteractionFactorsdpdm;
	if (M12InteractionFactorsumdpEnt)
	  delete[] M12InteractionFactorsumdpEnt;
	if (M12InteractionFactorsupdmEnt)
	  delete[] M12InteractionFactorsupdmEnt;
      }
  if (this->OneBodyInteractionFactorsupup != 0)
    {
      delete[] this->OneBodyInteractionFactorsupup;
      delete[] this->OneBodyInteractionFactorsumum;
      delete[] this->OneBodyInteractionFactorsdpdp;
      delete[] this->OneBodyInteractionFactorsdmdm;
    }
  if (this->S2Hamiltonian != 0)
    delete this->S2Hamiltonian;
  if (this->L2Hamiltonian != 0)
    delete this->L2Hamiltonian;
}

// ask if Hamiltonian implements hermitian symmetry operations
//
bool AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::IsHermitian()
{
  return false;
}

// ask if Hamiltonian implements conjugate methods
//
bool AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::IsConjugate()
{
  return false;
}

// add an additional S^2 term to the Hamiltonian
//
// totalLz = twice the projected momentum total value
// totalSz = twice the projected spin total value
// factor = factor in front of the S^2
// memory = amount of memory that can be used for S^2  precalculations 

void AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::AddS2 (int totalLz, int totalSz, double factor, long memory)
{
  this->S2Hamiltonian = new ParticleOnSphereWithSU4SpinS2Hamiltonian((ParticleOnSphereWithSU4Spin*)this->Particles, this->NbrParticles, this->LzMax, totalLz, totalSz, this->Architecture, factor, memory);
}

// add an additional L^2 term to the Hamiltonian
//
// totalLz = twice the projected momentum total value
// totalSz = twice the projected spin total value
// factor = factor in front of the L^2
// memory = amount of memory that can be used for L^2  precalculations 

void AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::AddL2 (int totalLz, int totalSz, double factor, long memory)
{
  this->L2Hamiltonian = new ParticleOnSphereWithSU4SpinL2Hamiltonian((ParticleOnSphereWithSU4Spin*)this->Particles, this->NbrParticles, this->LzMax, totalLz, this->Architecture, factor, memory);
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnSphereWithSU4Spin* TmpParticles = (ParticleOnSphereWithSU4Spin*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
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
	      ParticleOnSphereWithSU4Spin* TmpParticles = (ParticleOnSphereWithSU4Spin*) this->Particles->Clone();
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
		    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
		      this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
		    this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSource, vDestination);
		  }
	      delete TmpParticles;
	    }
	  else
	    {
	      int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
	      double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
	      int TmpNbrIteration = nbrComponent / this->BufferSize;
	      int* TmpIndexArray;
	      double* TmpCoefficientArray;
	      int TmpNbrInteraction;
	      int k = firstComponent;
	      int EffectiveHilbertSpaceDimension;
	      firstComponent -= this->PrecalculationShift;

	      ifstream File;
	      File.open(this->DiskStorageFileName, ios::binary | ios::in);
	      File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
	      long FileJump = 0;
	      for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
		FileJump += (long) this->NbrInteractionPerComponent[i];
	      FileJump *= sizeof(int);
	      long FileOffset = 0;
	      for (int i = this->DiskStorageStart; i < firstComponent; ++i)
		FileOffset += this->NbrInteractionPerComponent[i];
	      File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
	      FileJump += (sizeof(double) - sizeof(int)) * FileOffset;

	      for (int i = 0; i < TmpNbrIteration; ++i)
		{
		  int TmpPos = firstComponent;
		  long ReadBlockSize = 0;
		  for (int j = 0; j < this->BufferSize; ++j)
		    {
		      ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
		      ++TmpPos;
		    }		  
		  File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
		  FileJump -= sizeof(int) * ReadBlockSize;
		  File.seekg (FileJump, ios::cur);
		  File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
		  FileJump += sizeof(double) * ReadBlockSize;
		  File.seekg (-FileJump, ios::cur);
		  
		  TmpIndexArray = BufferIndexArray;
		  TmpCoefficientArray = BufferCoefficientArray;
		  for (int l = 0; l < this->BufferSize; ++l)
		    {
		      TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
		      Coefficient = vSource[k];
		      if (TmpNbrInteraction > 0)
			{
			  for (int j = 0; j < TmpNbrInteraction; ++j)
			    vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
			  TmpIndexArray += TmpNbrInteraction;
			  TmpCoefficientArray += TmpNbrInteraction;
			}
		      vDestination[k] += this->HamiltonianShift * Coefficient;
		      ++k;
		      ++firstComponent;
		    }
		}

	      if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
		{
		  int TmpPos = firstComponent;
		  int Lim =  nbrComponent % this->BufferSize;
		  long ReadBlockSize = 0;
		  for (int j = 0; j < Lim; ++j)
		    {
		      ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
		      ++TmpPos;
		    }		  
		  File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
		  FileJump -= sizeof(int) * ReadBlockSize;
		  File.seekg (FileJump, ios::cur);
		  File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
		  FileJump += sizeof(double) * ReadBlockSize;
		  File.seekg (-FileJump, ios::cur);

		  TmpIndexArray = BufferIndexArray;
		  TmpCoefficientArray = BufferCoefficientArray;
		  for (int i = 0; i < Lim; ++i)
		    {
		      TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
		      Coefficient = vSource[k];
		      if (TmpNbrInteraction > 0)
			{
			  for (int j = 0; j < TmpNbrInteraction; ++j)
			    vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
			  TmpIndexArray += TmpNbrInteraction;
			  TmpCoefficientArray += TmpNbrInteraction;
			}
		      vDestination[k] += this->HamiltonianShift * Coefficient;
		      ++k;
		      ++firstComponent;
		    }
		}

	      File.close();
	      delete[] BufferIndexArray;
	      delete[] BufferCoefficientArray;
	    }
	  
	}
    }

  if (this->L2Hamiltonian != 0)
    this->L2Hamiltonian->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  if (this->S2Hamiltonian != 0)
    this->S2Hamiltonian->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);


  return vDestination;
}

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      double* Coefficient2 = new double [nbrVectors];
      ParticleOnSphereWithSU4Spin* TmpParticles = (ParticleOnSphereWithSU4Spin*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
      delete[] Coefficient2;
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* Coefficient2 = new double [nbrVectors];
	  double* TmpCoefficientArray; 
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
		    vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
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

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
												      int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  ParticleOnSphereWithSU4Spin* TmpParticles = (ParticleOnSphereWithSU4Spin*) this->Particles->Clone();
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int Pos2;
  double Coefficient;
  int PosMod = firstComponent % this->FastMultiplicationStep;
  double* Coefficient2 = new double [nbrVectors];
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
	  vDestinations[k][l] += this->HamiltonianShift * Coefficient2[k];
	}
      for (j = 0; j < TmpNbrInteraction; ++j)
	{
	  Pos2 = TmpIndexArray[j];
	  Coefficient = TmpCoefficientArray[j];
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient2[k];
	}
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (l = 0; l < this->FastMultiplicationStep; ++l)
    if (PosMod != l)
      {
	for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
	  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
	this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSources, vDestinations, nbrVectors);
      }
  delete[] Coefficient2;
  delete TmpParticles;
  return vDestinations;
}


// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long Memory = 0;
  ParticleOnSphereWithSU4Spin* TmpParticles = (ParticleOnSphereWithSU4Spin*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;
  this->EvaluateMNTwoBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory);
  delete TmpParticles;
  return Memory;
}


// enable fast multiplication algorithm
//

void AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{
  int LastComponent = nbrComponent + firstComponent;

  int* TmpIndexArray;
  double* TmpCoefficientArray;
  long ColumnIndex;
  ParticleOnSphereWithSU4Spin* TmpParticles = (ParticleOnSphereWithSU4Spin*) this->Particles->Clone();

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
      this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, ColumnIndex);
      this->EvaluateMNOneBodyFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, ColumnIndex);
      ++Pos;
    }
  
  delete TmpParticles;
}
  

// enable fast multiplication algorithm using on disk cache 
//
// fileName = prefix of the name of the file where temporary matrix elements will be stored

void AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
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

  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  this->DiskStorageStart = (int) MinIndex;
  int DiskStorageEnd = 1 + (int) MaxIndex;

  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  this->InteractionPerComponentIndex = 0;
  this->InteractionPerComponentCoefficient = 0;
  this->MaxNbrInteractionPerComponent = 0;

  int TotalPos = 0;
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
// 	  for (int k = 2; k <= this->MaxNBody; ++k)
// 	    if (this->NBodyFlags[k] == true)
// 	      {
// 		double Sign = this->NBodySign[k];
// 		int TmpMinSumIndices = this->MinSumIndices[k];
// 		int TmpMaxSumIndices = this->MaxSumIndices[k];	      
// 		for (int j = TmpMinSumIndices; j <= TmpMaxSumIndices; ++j)
// 		  {
// 		    TmpInteraction = this->NBodyInteractionFactors[k][j];
// 		    int Lim = this->NbrSortedIndicesPerSum[k][j];
// 		    NIndices = this->SortedIndicesPerSum[k][j];
// 		    for (int i1 = 0; i1 < Lim; ++i1)
// 		      {
// 			Coefficient2 = Sign * this->Particles->ProdA(i, NIndices, k);
// 			if (Coefficient2 != 0.0)
// 			  {
// 			    MIndices = this->SortedIndicesPerSum[k][j];
// 			    for (int i2 = 0; i2 < Lim; ++i2)
// 			      {
// 				Index = this->Particles->ProdAd(MIndices, k, Coefficient);
// 				if (Index < this->Particles->GetHilbertSpaceDimension())
// 				  {
// 				    TmpIndexArray[Pos] = Index;
// 				    TmpCoefficientArray[Pos] = Coefficient2 * Coefficient * TmpInteraction[i1] *  TmpInteraction[i2];
// 				    ++Pos;
// 				  }
// 				MIndices += k;
// 			      }
// 			  }
// 			NIndices += k;
// 		      }
// 		  }
// 	      }
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

