////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                          to particles on a sphere with                     //
//                                                                            //
//                        last modification : 24/03/2003                      //
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
#include "Hamiltonian/AbstractSUNSpinOnLatticeHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "HilbertSpace/GenericSUNSpinCollection.h"

#include "GeneralTools/RealUniqueArray.h"
#include "GeneralTools/ComplexUniqueArray.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/SUNSpinPrecalculationOperation.h"

#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#include <limits.h>
#include <fstream>
#include <cstring>

using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ostream;


// default constructor
//
AbstractSUNSpinOnLatticeHamiltonian::AbstractSUNSpinOnLatticeHamiltonian()
{
  this->NbrPermutationTerms = 0;
  this->NbrCyclicPermutations = 0;
}

// destructor
//

AbstractSUNSpinOnLatticeHamiltonian::~AbstractSUNSpinOnLatticeHamiltonian()
{
  if (this->FastMultiplicationFlag == true)
    {
      long MinIndex;
      long MaxIndex;
      this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
      int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      int ReducedDim = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      if ((ReducedDim * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	++ReducedDim;
      for (int i = 0; i < ReducedDim; ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficientIndex[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficientIndex;
      delete[] this->NbrRealInteractionPerComponent;
      if (this->HaveComplexInteractions)
	delete[] this->NbrComplexInteractionPerComponent;
    }
  if (this->NbrPermutationTerms != 0)
    {
      delete[] this->PermutationPrefactors;
      delete[] this->PermutationI;
      delete[] this->PermutationJ;
    }
  if (this->NbrCyclicPermutations != 0)
    {
      delete[] this->CyclicPermutationPrefactors;
      delete[] this->CyclicPermutationLength;
      for (int i=0; i<this->NbrCyclicPermutations; ++i)
	delete[] this->CyclicPermutationIndices[i];
      delete[] this->CyclicPermutationIndices;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractSUNSpinOnLatticeHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  if (this->NbrPermutationTerms != 0)
    {
      delete[] this->PermutationPrefactors;
      delete[] this->PermutationI;
      delete[] this->PermutationJ;
    }
  if (this->NbrCyclicPermutations != 0)
    {
      delete[] this->CyclicPermutationPrefactors;
      delete[] this->CyclicPermutationLength;
      for (int i=0; i<this->NbrCyclicPermutations; ++i)
	delete[] this->CyclicPermutationIndices[i];
      delete[] this->CyclicPermutationIndices;
    }
  this->Spins = (GenericSUNSpinCollection*) hilbertSpace;
  this->EvaluateInteractionTerms();
  if (this->FastMultiplicationFlag == true)
    {
      long MinIndex;
      long MaxIndex;
      this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
      int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      int ReducedDim = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      if ((ReducedDim * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	++ReducedDim;
      for (int i = 0; i < ReducedDim; ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficientIndex[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficientIndex;
      delete[] this->NbrRealInteractionPerComponent;
      if (this->HaveComplexInteractions)
	delete[] this->NbrComplexInteractionPerComponent;
      this->EnableFastMultiplication();
    }
  
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractSUNSpinOnLatticeHamiltonian::GetHilbertSpace ()
{
  return this->Spins;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractSUNSpinOnLatticeHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Spins->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractSUNSpinOnLatticeHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractSUNSpinOnLatticeHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Spins->GetHilbertSpaceDimension();
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

Complex AbstractSUNSpinOnLatticeHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Spins->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
							       int firstComponent, int nbrComponent) 
{
  vDestination.ClearVector();
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and store result in another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors where result has to be stored
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractSUNSpinOnLatticeHamiltonian::LowLevelMultipleMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
								     int firstComponent, int nbrComponent)
{
  for (int i = 0; i < nbrVectors; ++i)
    vDestinations[i].ClearVector();
  return LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Spins->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Spins->GetHilbertSpaceDimension();
  int ReducedNbrPermutationTerms = NbrPermutationTerms-1;
  if (this->FastMultiplicationFlag == false)
    {
      if (this->HaveComplexInteractions)
	{
	  cout << "Attention, this is a complex Hamiltonian applied to RealVectors."<<endl;
	}
      int Index;
      int s1;
      int s2;
      double TmpInteraction;
      GenericSUNSpinCollection* TmpSpins = (GenericSUNSpinCollection*) this->Spins->Clone();
      for (int p=0; p<ReducedNbrPermutationTerms; ++p)
	{
	  s1 = this->PermutationI[p];
	  s2 = this->PermutationJ[p];
	  TmpInteraction = this->PermutationPrefactors[p];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpSpins->SpinPermutation(i, s1, s2);
	      if (Index < Dim)
		vDestination[Index] += TmpInteraction * vSource[i];
	    }
	}
      s1 = this->PermutationI[ReducedNbrPermutationTerms];
      s2 = this->PermutationJ[ReducedNbrPermutationTerms];
      TmpInteraction = this->PermutationPrefactors[ReducedNbrPermutationTerms];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  Index = TmpSpins->SpinPermutation(i, s1, s2);
	  if (Index < Dim)
	    vDestination[Index] += TmpInteraction * vSource[i];
	  vDestination[i] += this->HamiltonianShift * vSource[i];
	}
      delete TmpSpins;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  unsigned short* TmpCoefficientIndexArray;
	  double TmpSource;
	  unsigned short TmpNbrRealInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
	      TmpSource = vSource[k];
	      int Pos=0;
	      for (; Pos < TmpNbrRealInteraction; ++Pos)
		{
		  vDestination[TmpIndexArray[Pos]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpSource;
		}
	      vDestination[k] += this->HamiltonianShift * TmpSource;
	      ++k;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Spins->GetHilbertSpaceDimension();
  double TmpSource;
  GenericSUNSpinCollection* TmpSpins = (GenericSUNSpinCollection*) this->Spins->Clone();
  int* TmpIndexArray;
  short unsigned int* TmpCoefficientIndexArray;
  int TmpNbrRealInteraction;
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
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpSource = vSource[l];
      for (int j = 0; j < TmpNbrRealInteraction; ++j)
	vDestination[TmpIndexArray[j]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[j]]*TmpSource;
      vDestination[l] += this->HamiltonianShift * TmpSource;
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  int Index;
  int s1;
  int s2;
  double TmpInteraction;
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  int ReducedNbrPermutationTerms=NbrPermutationTerms-1;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	for (int p=0; p<ReducedNbrPermutationTerms; ++p)
	  {
	    s1 = this->PermutationI[p];
	    s2 = this->PermutationJ[p];
	    TmpInteraction = this->PermutationPrefactors[p];
	    for (int i = firstComponent+k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Index = TmpSpins->SpinPermutation(i, s1, s2);
		if (Index < Dim)
		  vDestination[Index] += TmpInteraction * vSource[i];
	      }
	  }
	s1 = this->PermutationI[ReducedNbrPermutationTerms];
	s2 = this->PermutationJ[ReducedNbrPermutationTerms];
	TmpInteraction = this->PermutationPrefactors[ReducedNbrPermutationTerms];
	for (int i = firstComponent+k; i < LastComponent; i += this->FastMultiplicationStep)
	  {
	    Index = TmpSpins->SpinPermutation(i, s1, s2);
	    if (Index < Dim)
	      vDestination[Index] += TmpInteraction * vSource[i];
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	  }
      }
  delete TmpSpins;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  cout << "Disk storage in AbstractSUNSpinOnLatticeHamiltonian not implemented yet."<<endl;
  /*
  double Coefficient;
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
  */
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

RealVector* AbstractSUNSpinOnLatticeHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Spins->GetHilbertSpaceDimension();
  int ReducedNbrPermutationTerms = NbrPermutationTerms-1;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int s1;
      int s2;
      double TmpInteraction;
      GenericSUNSpinCollection* TmpSpins = (GenericSUNSpinCollection*) this->Spins->Clone();
      // check : faster to invert loop order?
      for (int p=0; p<ReducedNbrPermutationTerms; ++p)
	{
	  s1 = this->PermutationI[p];
	  s2 = this->PermutationJ[p];
	  TmpInteraction = this->PermutationPrefactors[p];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpSpins->SpinPermutation(i, s1, s2);
	      if (Index < Dim)
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][Index] += TmpInteraction * vSources[l][i];
	    }
	}
      s1 = this->PermutationI[ReducedNbrPermutationTerms];
      s2 = this->PermutationJ[ReducedNbrPermutationTerms];
      TmpInteraction = this->PermutationPrefactors[ReducedNbrPermutationTerms];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  Index = TmpSpins->SpinPermutation(i, s1, s2);
	  if (Index < Dim)
	    for (int l = 0; l < nbrVectors; ++l)
	      vDestinations[l][Index] += TmpInteraction * vSources[l][i];
	  for (int l = 0; l < nbrVectors; ++l)
	    vDestinations[l][i] += this->HamiltonianShift * vSources[l][i];
	}
      delete TmpSpins;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;	  
	  unsigned short* TmpCoefficientIndexArray;
	  unsigned short TmpNbrRealInteraction;
	  int k = firstComponent;
	  double *TmpSources = new double[nbrVectors];
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
	      for (int l = 0; l < nbrVectors; ++l)
		TmpSources[l] = vSources[l][k];
	      int Pos=0;
	      for (; Pos < TmpNbrRealInteraction; ++Pos)
		{
		  
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][TmpIndexArray[Pos]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpSources[l];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][k] += this->HamiltonianShift * TmpSources[l];
	      ++k;
	    }
	  delete [] TmpSources;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	}
    }
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

RealVector* AbstractSUNSpinOnLatticeHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Spins->GetHilbertSpaceDimension();
  double* TmpSources = new double[nbrVectors];
  GenericSUNSpinCollection* TmpSpins = (GenericSUNSpinCollection*) this->Spins->Clone();
  int* TmpIndexArray;
  short unsigned int* TmpCoefficientIndexArray;
  int TmpNbrRealInteraction;
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
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      for (int k=0; k<nbrVectors; ++k)
	TmpSources[k] = vSources[k][l];
      for (int j = 0; j < TmpNbrRealInteraction; ++j)
	for (int k=0; k<nbrVectors; ++k)
	  vDestinations[k][TmpIndexArray[j]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[j]]*TmpSources[k];
      for (int k=0; k<nbrVectors; ++k)
	vDestinations[k][l] += this->HamiltonianShift * TmpSources[k];
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  int Index;
  int s1;
  int s2;
  double TmpInteraction;
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  int ReducedNbrPermutationTerms=NbrPermutationTerms-1;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	for (int p=0; p<ReducedNbrPermutationTerms; ++p)
	  {
	    s1 = this->PermutationI[p];
	    s2 = this->PermutationJ[p];
	    TmpInteraction = this->PermutationPrefactors[p];
	    for (int i = firstComponent+k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Index = TmpSpins->SpinPermutation(i, s1, s2);
		if (Index < Dim)
		  for (int v=0; v<nbrVectors; ++v)
		    vDestinations[v][Index] += TmpInteraction * vSources[v][i];
	      }
	  }
	s1 = this->PermutationI[ReducedNbrPermutationTerms];
	s2 = this->PermutationJ[ReducedNbrPermutationTerms];
	TmpInteraction = this->PermutationPrefactors[ReducedNbrPermutationTerms];
	for (int i = firstComponent+k; i < LastComponent; i += this->FastMultiplicationStep)
	  {
	    Index = TmpSpins->SpinPermutation(i, s1, s2);
	    if (Index < Dim)
	      for (int v=0; v<nbrVectors; ++v)
		vDestinations[v][Index] += TmpInteraction * vSources[v][i];
	    for (int v=0; v<nbrVectors; ++v)
	      vDestinations[v][i] += this->HamiltonianShift * vSources[v][i];
	  }
      }
  delete[] TmpSources;
  delete TmpSpins;
  return vDestinations;
}

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using disk storage option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractSUNSpinOnLatticeHamiltonian::LowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										   int firstComponent, int nbrComponent)
{
  cout << "AbstractSUNSpinOnLatticeHamiltonian::LowLevelMultipleAddMultiplyDiskStorage not defined, yet"<<endl;
  /* old code
  double Coefficient;
  int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* Coefficient2 = new double [nbrVectors];
  int TmpNbrIteration = nbrComponent / this->BufferSize;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int TmpNbrInteraction;
  int k = firstComponent;
  int EffectiveHilbertSpaceDimension;
  int Pos;
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
      for (int m = 0; m < this->BufferSize; ++m)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  for (int l = 0; l < nbrVectors; ++l)
	    {
	      Coefficient2[l] = vSources[l][k];
	      vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
	    }
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		    }
		}
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
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
      for (int m = 0; m < Lim; ++m)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  for (int l = 0; l < nbrVectors; ++l)
	    {
	      Coefficient2[l] = vSources[l][k];
	      vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
	    }
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		    }
		}
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  ++k;
	  ++firstComponent;
	}
    }
  
  File.close();
  delete[] BufferIndexArray;
  delete[] BufferCoefficientArray;
  delete[] Coefficient2;
  */
  return vDestinations;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Spins->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							 int firstComponent, int nbrComponent)
{
  vDestination.ClearVector();
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Spins->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
ComplexVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Spins->GetHilbertSpaceDimension();
  int ReducedNbrPermutationTerms = NbrPermutationTerms-1;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int s1;
      int s2;
      double TmpInteraction;
      double TmpInteractionRe;
      double TmpInteractionIm;
      GenericSUNSpinCollection* TmpSpins = (GenericSUNSpinCollection*) this->Spins->Clone();
      if (HaveComplexInteractions)
	{
	  int ReducedNbrComplexPermutationTerms = NbrComplexPermutationTerms-1;
	  for (int p=0; p<ReducedNbrComplexPermutationTerms; ++p)
	    {
	      s1 = this->PermutationI[p];
	      s2 = this->PermutationJ[p];
	      TmpInteractionRe = this->ComplexPermutationPrefactors[p].Re;
	      TmpInteractionIm = this->ComplexPermutationPrefactors[p].Im;
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpSpins->SpinPermutation(i, s1, s2);
		  if (Index < Dim)
		    {
		      vDestination.Re(Index) += TmpInteractionRe * vSource[i].Re - TmpInteractionIm * vSource[i].Im;
		      vDestination.Im(Index) += TmpInteractionIm * vSource[i].Re + TmpInteractionRe * vSource[i].Im;
		    }
		}
	    }
	  s1 = this->PermutationI[ReducedNbrComplexPermutationTerms];
	  s2 = this->PermutationJ[ReducedNbrComplexPermutationTerms];
	  TmpInteractionRe = this->ComplexPermutationPrefactors[ReducedNbrComplexPermutationTerms].Re;
	  TmpInteractionIm = this->ComplexPermutationPrefactors[ReducedNbrComplexPermutationTerms].Im;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpSpins->SpinPermutation(i, s1, s2);
	      if (Index < Dim)
		{
		  vDestination.Re(Index) += TmpInteractionRe * vSource[i].Re - TmpInteractionIm * vSource[i].Im;
		  vDestination.Im(Index) += TmpInteractionIm * vSource[i].Re + TmpInteractionRe * vSource[i].Im;
		}
	      vDestination.Re(i) += this->HamiltonianShift * vSource[i].Re;
	      vDestination.Im(i) += this->HamiltonianShift * vSource[i].Im;
	    }
	}
      
      for (int p=0; p<ReducedNbrPermutationTerms; ++p)
	{
	  s1 = this->PermutationI[p];
	  s2 = this->PermutationJ[p];
	  TmpInteraction = this->PermutationPrefactors[p];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpSpins->SpinPermutation(i, s1, s2);
	      if (Index < Dim)
		{
		  vDestination.Re(Index) += TmpInteraction * vSource[i].Re;
		  vDestination.Im(Index) += TmpInteraction * vSource[i].Im;
		}
	    }
	}
      if (NbrPermutationTerms>0)
	{
	  s1 = this->PermutationI[ReducedNbrPermutationTerms];
	  s2 = this->PermutationJ[ReducedNbrPermutationTerms];
	  TmpInteraction = this->PermutationPrefactors[ReducedNbrPermutationTerms];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpSpins->SpinPermutation(i, s1, s2);
	      if (Index < Dim)
		{
		  vDestination.Re(Index) += TmpInteraction * vSource[i].Re;
		  vDestination.Im(Index) += TmpInteraction * vSource[i].Im;
		}
	      vDestination.Re(i) += this->HamiltonianShift * vSource[i].Re;
	      vDestination.Im(i) += this->HamiltonianShift * vSource[i].Im;
	    }
	}
      else
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      vDestination.Re(i) += this->HamiltonianShift * vSource[i].Re;
	      vDestination.Im(i) += this->HamiltonianShift * vSource[i].Im;
	    }
	}
      // complex cyclic permutations
      int *s,l;
      for (int p=0; p<NbrCyclicPermutations; ++p)
	{
	  s = this->CyclicPermutationIndices[p];
	  l = this->CyclicPermutationLength[p];
	  TmpInteractionRe = this->CyclicPermutationPrefactors[p].Re;
	  TmpInteractionIm = this->CyclicPermutationPrefactors[p].Im;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpSpins->CyclicSpinPermutation(i, l, s);
	      if (Index < Dim)
		{
		  vDestination.Re(Index) += TmpInteractionRe * vSource[i].Re - TmpInteractionIm * vSource[i].Im;
		  vDestination.Re(Index) += TmpInteractionRe * vSource[i].Im + TmpInteractionIm * vSource[i].Re;
		}
	    }
	}
      delete TmpSpins;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  unsigned short* TmpCoefficientIndexArray;
	  double TmpRe;
	  double TmpIm;
	  Complex *TmpCPtr;
	  unsigned short TmpNbrRealInteraction;
	  unsigned short TmpNbrComplexInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
	      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
	      TmpRe = vSource[k].Re;
	      TmpIm = vSource[k].Im;
	      int Pos=0;
	      for (; Pos < TmpNbrRealInteraction; ++Pos)
		{
		  vDestination.Re(TmpIndexArray[Pos]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpRe;
		  vDestination.Im(TmpIndexArray[Pos]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpIm;
		}
	      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
		{
		  TmpCPtr= &(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]]);
		  vDestination.Re(TmpIndexArray[Pos]) +=  TmpCPtr->Re*TmpRe-TmpCPtr->Im*TmpIm;
		  vDestination.Im(TmpIndexArray[Pos]) +=  TmpCPtr->Re*TmpIm+TmpCPtr->Im*TmpRe;		  
		}
	      vDestination.Re(k) += this->HamiltonianShift * TmpRe;
	      vDestination.Im(k) += this->HamiltonianShift * TmpIm;	      
	      ++k;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Spins->GetHilbertSpaceDimension();
  double TmpRe;
  double TmpIm;
  Complex *TmpCPtr;
  GenericSUNSpinCollection* TmpSpins = (GenericSUNSpinCollection*) this->Spins->Clone();
  int* TmpIndexArray;
  short unsigned int* TmpCoefficientIndexArray;
  int TmpNbrRealInteraction;
  int TmpNbrComplexInteraction;
  double TmpInteractionRe;
  double TmpInteractionIm;
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
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpRe = vSource[l].Re;
      TmpIm = vSource[l].Im;
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	{
	  vDestination.Re(TmpIndexArray[Pos2]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpRe;
	  vDestination.Im(TmpIndexArray[Pos2]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpIm;
	}
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos2)
	{
	  TmpCPtr= &(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]);
	  vDestination.Re(TmpIndexArray[Pos2]) +=  TmpCPtr->Re*TmpRe-TmpCPtr->Im*TmpIm;
	  vDestination.Im(TmpIndexArray[Pos2]) +=  TmpCPtr->Re*TmpIm+TmpCPtr->Im*TmpRe;		  
	}
      vDestination.Re(l) += this->HamiltonianShift * TmpRe;
      vDestination.Im(l) += this->HamiltonianShift * TmpIm;	      
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  int Index;
  int s1;
  int s2;
  double TmpInteraction;
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  int ReducedNbrPermutationTerms=NbrPermutationTerms-1;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	if (HaveComplexInteractions)
	{
	  int ReducedNbrComplexPermutationTerms = NbrComplexPermutationTerms-1;
	  for (int p=0; p<ReducedNbrComplexPermutationTerms; ++p)
	    {
	      s1 = this->PermutationI[p];
	      s2 = this->PermutationJ[p];
	      TmpInteractionRe = this->ComplexPermutationPrefactors[p].Re;
	      TmpInteractionIm = this->ComplexPermutationPrefactors[p].Im;
	      for (int i = firstComponent+k; i < LastComponent; i += this->FastMultiplicationStep)
		{
		  Index = TmpSpins->SpinPermutation(i, s1, s2);
		  if (Index < Dim)
		    {
		      vDestination.Re(Index) += TmpInteractionRe * vSource[i].Re - TmpInteractionIm * vSource[i].Im;
		      vDestination.Im(Index) += TmpInteractionIm * vSource[i].Re + TmpInteractionRe * vSource[i].Im;
		    }
		}
	    }
	  s1 = this->PermutationI[ReducedNbrComplexPermutationTerms];
	  s2 = this->PermutationJ[ReducedNbrComplexPermutationTerms];
	  TmpInteractionRe = this->ComplexPermutationPrefactors[ReducedNbrComplexPermutationTerms].Re;
	  TmpInteractionIm = this->ComplexPermutationPrefactors[ReducedNbrComplexPermutationTerms].Im;
	  for (int i = firstComponent+k; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      Index = TmpSpins->SpinPermutation(i, s1, s2);
	      if (Index < Dim)
		{
		  vDestination.Re(Index) += TmpInteractionRe * vSource[i].Re - TmpInteractionIm * vSource[i].Im;
		  vDestination.Im(Index) += TmpInteractionIm * vSource[i].Re + TmpInteractionRe * vSource[i].Im;
		}
	      vDestination.Re(i) += this->HamiltonianShift * vSource[i].Re;
	      vDestination.Im(i) += this->HamiltonianShift * vSource[i].Im;
	    }
	}
      
      for (int p=0; p<ReducedNbrPermutationTerms; ++p)
	{
	  s1 = this->PermutationI[p];
	  s2 = this->PermutationJ[p];
	  TmpInteraction = this->PermutationPrefactors[p];
	  for (int i = firstComponent+k; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      Index = TmpSpins->SpinPermutation(i, s1, s2);
	      if (Index < Dim)
		{
		  vDestination.Re(Index) += TmpInteraction * vSource[i].Re;
		  vDestination.Im(Index) += TmpInteraction * vSource[i].Im;
		}
	    }
	}
      if (NbrPermutationTerms>0)
	{
	  s1 = this->PermutationI[ReducedNbrPermutationTerms];
	  s2 = this->PermutationJ[ReducedNbrPermutationTerms];
	  TmpInteraction = this->PermutationPrefactors[ReducedNbrPermutationTerms];
	  for (int i = firstComponent+k; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      Index = TmpSpins->SpinPermutation(i, s1, s2);
	      if (Index < Dim)
		{
		  vDestination.Re(Index) += TmpInteraction * vSource[i].Re;
		  vDestination.Im(Index) += TmpInteraction * vSource[i].Im;
		}
	      vDestination.Re(i) += this->HamiltonianShift * vSource[i].Re;
	      vDestination.Im(i) += this->HamiltonianShift * vSource[i].Im;
	    }
	}
      else
	{
	  for (int i = firstComponent+k; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      vDestination.Re(i) += this->HamiltonianShift * vSource[i].Re;
	      vDestination.Im(i) += this->HamiltonianShift * vSource[i].Im;
	    }
	}
      // complex cyclic permutations
      int *s,l;
      for (int p=0; p<NbrCyclicPermutations; ++p)
	{
	  s = this->CyclicPermutationIndices[p];
	  l = this->CyclicPermutationLength[p];
	  TmpInteractionRe = this->CyclicPermutationPrefactors[p].Re;
	  TmpInteractionIm = this->CyclicPermutationPrefactors[p].Im;
	  for (int i = firstComponent+k; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      Index = TmpSpins->CyclicSpinPermutation(i, l, s);
	      if (Index < Dim)
		{
		  vDestination.Re(Index) += TmpInteractionRe * vSource[i].Re - TmpInteractionIm * vSource[i].Im;
		  vDestination.Re(Index) += TmpInteractionRe * vSource[i].Im + TmpInteractionIm * vSource[i].Re;
		}
	    }
	}
      }
  delete TmpSpins;
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
ComplexVector& AbstractSUNSpinOnLatticeHamiltonian::LowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  cout << "Attention, disk storage not implemented in AbstractSUNSpinOnLatticeHamiltonian"<<endl;
  return vDestination;
}


 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AbstractSUNSpinOnLatticeHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AbstractSUNSpinOnLatticeHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}


// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long AbstractSUNSpinOnLatticeHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  this->LoadedPrecalculation=false;
  if (allowedMemory>0)
    this->AllowedMemory = allowedMemory;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  
  // storage per basis state:
  // 1 x array of unsigned short NbrRealInteractionPerComponent
  // 1 x array of unsigned short NbrComplexInteractionPerComponent (if with complex)
  // 1 x array of int* InteractionPerComponentIndex
  // 1 x array of unsigned short* InteractionPerComponentCoefficientIndex
  int StoragePerBasisState;
  if ((this->HaveComplexInteractions) || (NbrCyclicPermutations>0))
    StoragePerBasisState = (2*sizeof (unsigned short) + sizeof (int*) + sizeof(unsigned short*));
  else
    StoragePerBasisState = (sizeof (unsigned short) + sizeof (int*) + sizeof(unsigned short*));
  // storage per non-zero matrix element
  // 1 x target index ~ int
  // 1 x value index ~ unsigned short
  int StoragePerMatrixElement = sizeof (int) + sizeof(unsigned short);

  this->NbrRealInteractionPerComponent = new unsigned short [EffectiveHilbertSpaceDimension];
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    this->NbrRealInteractionPerComponent[i] = 0x0;
  this->NbrComplexInteractionPerComponent = NULL;
  if ((this->HaveComplexInteractions)||(NbrCyclicPermutations>0))
    {
      this->NbrComplexInteractionPerComponent = new unsigned short [EffectiveHilbertSpaceDimension];   
      for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
	this->NbrComplexInteractionPerComponent[i] = 0x0;
    }
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start memory" << endl;
  
  
  SUNSpinPrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);

  long Memory = 0;
  if ((this->HaveComplexInteractions)||(NbrCyclicPermutations>0))
    for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
      {
	Memory += this->NbrRealInteractionPerComponent[i];
	Memory += this->NbrComplexInteractionPerComponent[i];
      }
  else
    for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
      Memory += this->NbrRealInteractionPerComponent[i];
  
  cout << "nbr interaction = " << Memory << endl;

  // memory requirement, ignoring the actual storage size of the values of matrix
  // elements, which is assumed small (maybe need to add an estimate, at least)
  long TmpMemory = AllowedMemory - StoragePerBasisState * EffectiveHilbertSpaceDimension;
  if (TmpMemory / ((int) (StoragePerMatrixElement)) < Memory)
    cout << "of which can be stored: "<<(TmpMemory / ((int) (StoragePerMatrixElement)))<<endl;
  // else cout << "all of which can be stored"<<endl;  
  
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (StoragePerMatrixElement))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      //this->FastMultiplicationSubStep = 0;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || (TmpMemory / ((int) (StoragePerMatrixElement)) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (this->Spins->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  // memory requirement, ignoring the actual storage size of the values of matrix
	  // elements, which is assumed small (maybe need to add an estimate, at least, again!)
	  TmpMemory = AllowedMemory - StoragePerBasisState * ReducedSpaceDimension;
	  Memory = 0;
	  if ((this->HaveComplexInteractions)||(NbrCyclicPermutations>0))
	    for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	      {
		Memory += this->NbrRealInteractionPerComponent[i];
		Memory += this->NbrComplexInteractionPerComponent[i];
	      }
	  else
	    for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	      Memory += this->NbrRealInteractionPerComponent[i];
	}
      Memory = (StoragePerBasisState * ReducedSpaceDimension) + (Memory * StoragePerMatrixElement);

      if (this->DiskStorageFlag == false)
	{
	  int TotalReducedSpaceDimension = ReducedSpaceDimension;
	  unsigned short* TmpNbrRealInteractionPerComponent = new unsigned short [TotalReducedSpaceDimension];
	  unsigned short* TmpNbrComplexInteractionPerComponent = NULL;
	  if ((this->HaveComplexInteractions)||(NbrCyclicPermutations>0))
	    TmpNbrComplexInteractionPerComponent = new unsigned short [TotalReducedSpaceDimension];	  
	  int Pos = 0;
	  if ((this->HaveComplexInteractions)||(NbrCyclicPermutations>0))
	    for (int i = 0; i < ReducedSpaceDimension; ++i)
	      {
		TmpNbrRealInteractionPerComponent[i] = this->NbrRealInteractionPerComponent[Pos];
		TmpNbrComplexInteractionPerComponent[i] = this->NbrComplexInteractionPerComponent[Pos];
		Pos += this->FastMultiplicationStep;
	      }
	  else
	    for (int i = 0; i < ReducedSpaceDimension; ++i)
	      {
		TmpNbrRealInteractionPerComponent[i] = this->NbrRealInteractionPerComponent[Pos];
		Pos += this->FastMultiplicationStep;
	      }
	  delete[] this->NbrRealInteractionPerComponent;
	  if ((this->HaveComplexInteractions)||(NbrCyclicPermutations>0))
	    delete[] this->NbrComplexInteractionPerComponent;
	  this->NbrRealInteractionPerComponent = TmpNbrRealInteractionPerComponent;
	  this->NbrComplexInteractionPerComponent = TmpNbrComplexInteractionPerComponent;
	}
    }
  else
    {
      Memory = (StoragePerBasisState * EffectiveHilbertSpaceDimension) + (Memory * StoragePerMatrixElement);
      this->FastMultiplicationStep = 1;
      // this->FastMultiplicationSubStep = 0;
    }

  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  return Memory;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalculated
// nbrComponent  = number of components that have to be precalculated
// return value = number of non-zero matrix element

long AbstractSUNSpinOnLatticeHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  long Memory = 0;
  int Dim = this->Spins->GetHilbertSpaceDimension();
  int Index;
  int s1;
  int s2;
  double TmpInteraction;
  GenericSUNSpinCollection* TmpSpins = (GenericSUNSpinCollection*) this->Spins->Clone();
  for (int p=0; p<NbrPermutationTerms; ++p)
    {
      s1 = this->PermutationI[p];
      s2 = this->PermutationJ[p];
      TmpInteraction = this->PermutationPrefactors[p];
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  Index = TmpSpins->SpinPermutation(i, s1, s2);
	  if (Index < Dim)
	    {
	      ++Memory;		
	      ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];		  
	    }
	}
    }
  if (this->HaveComplexInteractions)
    {
      Complex TmpComplexInteraction;
      int FinalNbrPermutationTerms = NbrPermutationTerms + NbrComplexPermutationTerms;
      for (int p=NbrPermutationTerms, f=0; p<FinalNbrPermutationTerms; ++p, ++f)
	{
	  s1 = this->PermutationI[p];
	  s2 = this->PermutationJ[p];
	  TmpComplexInteraction = this->ComplexPermutationPrefactors[f];	  
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpSpins->SpinPermutation(i, s1, s2);
	      if (Index < Dim)
		{
		  ++Memory;		
		  ++this->NbrComplexInteractionPerComponent[i - this->PrecalculationShift];		  
		}
	    }
	}
    }
  if (NbrCyclicPermutations>0)
    {
      int *s,l;
      Complex TmpComplexInteraction;
      for (int p=0; p<NbrCyclicPermutations; ++p)
	{
	  s = this->CyclicPermutationIndices[p];
	  l = this->CyclicPermutationLength[p];
	  TmpComplexInteraction = this->CyclicPermutationPrefactors[p];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpSpins->CyclicSpinPermutation(i, l, s);
	      if (Index < Dim)
		{
		  ++Memory;		
		  ++this->NbrComplexInteractionPerComponent[i - this->PrecalculationShift];		  
		}
	    }
	}
    }
  delete TmpSpins;
  return Memory;
}

// enable fast multiplication algorithm
//

void AbstractSUNSpinOnLatticeHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  int Dim = this->Spins->GetHilbertSpaceDimension();
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  int Index;
  int s1;
  int s2;
  GenericSUNSpinCollection* TmpSpins = (GenericSUNSpinCollection*) this->Spins->Clone();
  int tmpElementPos;
  int* TmpIndexArray;
  unsigned short* TmpCoefficientIndexArray;
  int PosR, PosC;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficientIndex = new unsigned short* [ReducedSpaceDimension];
  
  int TotalPos = 0;
  int TmpInteractionPerComponent;
  double TmpInteraction;
  Complex TmpComplexInteraction;
  
  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
    {
      TmpInteractionPerComponent = this->NbrRealInteractionPerComponent[TotalPos];
      if ((HaveComplexInteractions)||(NbrCyclicPermutations>0))
	TmpInteractionPerComponent+=this->NbrComplexInteractionPerComponent[TotalPos];
      this->InteractionPerComponentIndex[TotalPos] = new int [TmpInteractionPerComponent];
      this->InteractionPerComponentCoefficientIndex[TotalPos] = new unsigned short [TmpInteractionPerComponent];      
      TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[TotalPos];
      PosR = 0;  // counter for position of real matrix elements
      PosC = this->NbrRealInteractionPerComponent[TotalPos];  // counter for position of complex matrix elements
      
      // deal with kinetic energy terms first!
      for (int p=0; p<NbrPermutationTerms; ++p)
	{
	  s1 = this->PermutationI[p];
	  s2 = this->PermutationJ[p];
	  TmpInteraction = this->PermutationPrefactors[p];
	  Index = TmpSpins->SpinPermutation(i, s1, s2);
	  if (Index < Dim)
	    {
	      TmpIndexArray[PosR] = Index;
	      tmpElementPos = RealInteractionCoefficients.InsertElement(TmpInteraction);
	      if (tmpElementPos > USHRT_MAX )
		{
		  cout << "Error: too many different real matrix elements for fast storage"<<endl;
		  exit(1);
		}
	      TmpCoefficientIndexArray[PosR] = (unsigned short) tmpElementPos;
	      ++PosR;
	    }
	}
      if (HaveComplexInteractions)
	{
	  int FinalNbrPermutationTerms = NbrPermutationTerms + NbrComplexPermutationTerms;
	  for (int p=NbrPermutationTerms, f=0; p<FinalNbrPermutationTerms; ++p, ++f)
	    {
	      s1 = this->PermutationI[p];
	      s2 = this->PermutationJ[p];
	      TmpComplexInteraction = this->ComplexPermutationPrefactors[f];	  
	      Index = TmpSpins->SpinPermutation(i, s1, s2);
	      if (Index < Dim)
		{
		  TmpIndexArray[PosC] = Index;
		  tmpElementPos = ComplexInteractionCoefficients.InsertElement(TmpComplexInteraction);
		  if (tmpElementPos > USHRT_MAX )
		    {
		      cout << "Error: too many different complex matrix elements for fast storage"<<endl;
		      exit(1);
		    }
		  TmpCoefficientIndexArray[PosC] = (unsigned short) tmpElementPos;
		  ++PosC;
		}
	    }
	}
      if (NbrCyclicPermutations>0)
	{
	  int *s,l;
	  Complex TmpComplexInteraction;
	  for (int p=0; p<NbrCyclicPermutations; ++p)
	    {
	      s = this->CyclicPermutationIndices[p];
	      l = this->CyclicPermutationLength[p];
	      TmpComplexInteraction = this->CyclicPermutationPrefactors[p];
	      Index = TmpSpins->CyclicSpinPermutation(i, l, s);
	      if (Index < Dim)
		{
		  TmpIndexArray[PosC] = Index;
		  tmpElementPos = ComplexInteractionCoefficients.InsertElement(TmpComplexInteraction);
		  if (tmpElementPos > USHRT_MAX )
		    {
		      cout << "Error: too many different complex matrix elements for fast storage"<<endl;
		      exit(1);
		    }
		  TmpCoefficientIndexArray[PosC] = (unsigned short) tmpElementPos;
		  ++PosC;
		}
	    }
	}
      ++TotalPos;
    }
  delete TmpSpins;

  cout << "Nbr distinct matrix elements: "<<RealInteractionCoefficients.GetNbrElements()<<" real, "
       << ComplexInteractionCoefficients.GetNbrElements()<<" complex"<<endl;
   
  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted

void AbstractSUNSpinOnLatticeHamiltonian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
  int Dim = this->Spins->GetHilbertSpaceDimension();
  int Index;
  int s1;
  int s2;
  GenericSUNSpinCollection* TmpSpins = (GenericSUNSpinCollection*) this->Spins->Clone();
  int tmpElementPos;
  int* TmpIndexArray;
  unsigned short* TmpCoefficientIndexArray;
  int PosR, PosC;
  
  int Min = firstComponent / this->FastMultiplicationStep;
  int Max = lastComponent / this->FastMultiplicationStep;
  
  int TmpInteractionPerComponent;
  double TmpInteraction;
  Complex TmpComplexInteraction;

  for (int i = Min; i < Max; ++i)
    {
      TmpInteractionPerComponent = this->NbrRealInteractionPerComponent[i];
      if ((HaveComplexInteractions)||(NbrCyclicPermutations>0))
	TmpInteractionPerComponent+=this->NbrComplexInteractionPerComponent[i];
      this->InteractionPerComponentIndex[i] = new int [TmpInteractionPerComponent];
      this->InteractionPerComponentCoefficientIndex[i] = new unsigned short [TmpInteractionPerComponent];      
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
      PosR = 0;  // counter for position of real matrix elements
      PosC = this->NbrRealInteractionPerComponent[i];  // counter for position of complex matrix elements
      
      // deal with kinetic energy terms first!
      for (int p=0; p<NbrPermutationTerms; ++p)
	{
	  s1 = this->PermutationI[p];
	  s2 = this->PermutationJ[p];
	  TmpInteraction = this->PermutationPrefactors[p];
	  Index = TmpSpins->SpinPermutation(i, s1, s2);
	  if (Index < Dim)
	    {
	      TmpIndexArray[PosR] = Index;
	      tmpElementPos = RealInteractionCoefficients.InsertElement(TmpInteraction);
	      if (tmpElementPos > USHRT_MAX )
		{
		  cout << "Error: too many different real matrix elements for fast storage"<<endl;
		  exit(1);
		}
	      TmpCoefficientIndexArray[PosR] = (unsigned short) tmpElementPos;
	      ++PosR;
	    }
	}
      if (HaveComplexInteractions)
	{
	  int FinalNbrPermutationTerms = NbrPermutationTerms + NbrComplexPermutationTerms;
	  for (int p=NbrPermutationTerms, f=0; p<FinalNbrPermutationTerms; ++p, ++f)
	    {
	      s1 = this->PermutationI[p];
	      s2 = this->PermutationJ[p];
	      TmpComplexInteraction = this->ComplexPermutationPrefactors[f];	  
	      Index = TmpSpins->SpinPermutation(i, s1, s2);
	      if (Index < Dim)
		{
		  TmpIndexArray[PosC] = Index;
		  tmpElementPos = ComplexInteractionCoefficients.InsertElement(TmpComplexInteraction);
		  if (tmpElementPos > USHRT_MAX )
		    {
		      cout << "Error: too many different complex matrix elements for fast storage"<<endl;
		      exit(1);
		    }
		  TmpCoefficientIndexArray[PosC] = (unsigned short) tmpElementPos;
		  ++PosC;
		}
	    }
	}
      if (NbrCyclicPermutations>0)
	{
	  int *s,l;
	  Complex TmpComplexInteraction;
	  for (int p=0; p<NbrCyclicPermutations; ++p)
	    {
	      s = this->CyclicPermutationIndices[p];
	      l = this->CyclicPermutationLength[p];
	      TmpComplexInteraction = this->CyclicPermutationPrefactors[p];
	      Index = TmpSpins->CyclicSpinPermutation(i, l, s);
	      if (Index < Dim)
		{
		  TmpIndexArray[PosC] = Index;
		  tmpElementPos = ComplexInteractionCoefficients.InsertElement(TmpComplexInteraction);
		  if (tmpElementPos > USHRT_MAX )
		    {
		      cout << "Error: too many different complex matrix elements for fast storage"<<endl;
		      exit(1);
		    }
		  TmpCoefficientIndexArray[PosC] = (unsigned short) tmpElementPos;
		  ++PosC;
		}
	    }
	}
    }
  delete TmpSpins;
}

// enable fast multiplication algorithm using on disk cache 
//
// fileName = prefix of the name of the file where temporary matrix elements will be stored

void AbstractSUNSpinOnLatticeHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
{
  cout << "Disk storage in AbstractSUNSpinOnLatticeHamiltonian not implemented yet."<<endl;
  /*
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

  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  double Coefficient;
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
  double Coefficient2;
  int SumIndices;
  int TmpNbrM3Values;
  int* TmpM3Values;
  int ReducedNbrInteractionFactors;

  for (int i = this->DiskStorageStart; i < DiskStorageEnd; ++i)
    {
      if (this->NbrInteractionPerComponent[TotalPos] > 0)
	{
	  Pos = 0;
	  if (this->NbrM12Indices == 0)
	    {
	      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
		{
		  m1 = this->M1Value[j];
		  m2 = this->M2Value[j];
		  m3 = this->M3Value[j];
		  m4 = m1 + m2 - m3;
		  Index = this->Spins->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < this->Spins->GetHilbertSpaceDimension())
		    {
		      TmpIndexArray[Pos] = Index;
		      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
		      ++Pos;
		    }
		}
	    }
	  else
	    {
	      ReducedNbrInteractionFactors = 0;
	      for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		{
		  Coefficient = this->Spins->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1Value[m1] + this->M2Value[m1];
		      TmpM3Values = this->M3Values[m1];
		      TmpNbrM3Values = this->NbrM3Values[m1];
		      for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = this->Spins->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < this->Spins->GetHilbertSpaceDimension())
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
	  if (this->OneBodyTermFlag == true)
	    {
	      for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
		{
		  m1 = this->OneBodyMValues[j];
		  m2 = this->OneBodyNValues[j];
		  Index = this->Spins->AdA(i, m1, m2, Coefficient);
		  if (Index < this->Spins->GetHilbertSpaceDimension())
		    {
		      TmpIndexArray[Pos] = Index;
		      TmpCoefficientArray[Pos] = Coefficient * this->OneBodyInteractionFactors[j];
		      ++Pos;
		    }
		}
	    }
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
  */
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs

bool AbstractSUNSpinOnLatticeHamiltonian::SavePrecalculation (char* fileName)
{
  if (this->FastMultiplicationFlag)
    {
      ofstream File;
      File.open(fileName, ios::binary | ios::out);
      int Tmp = this->Spins->GetHilbertSpaceDimension();
      File.write((char*) &(Tmp), sizeof(int));
      File.write((char*) &(HaveComplexInteractions), sizeof(bool));
      File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
      Tmp /= this->FastMultiplicationStep;
      if ((Tmp * this->FastMultiplicationStep) != this->Spins->GetHilbertSpaceDimension())
	++Tmp;
      File.write((char*) this->NbrRealInteractionPerComponent, sizeof(unsigned short) * Tmp);
      if (this->HaveComplexInteractions)
	File.write((char*) this->NbrComplexInteractionPerComponent, sizeof(unsigned short) * Tmp);
      int TmpNbr;
      for (int i = 0; i < Tmp; ++i)
	{
	  TmpNbr = this->NbrRealInteractionPerComponent[i];
	  if (this->HaveComplexInteractions)
	    TmpNbr += this->NbrComplexInteractionPerComponent[i];
	  File.write((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * TmpNbr);
	}
      for (int i = 0; i < Tmp; ++i)
	{
	  TmpNbr = this->NbrRealInteractionPerComponent[i];
	  if (this->HaveComplexInteractions)
	    TmpNbr += this->NbrComplexInteractionPerComponent[i];
	  File.write((char*) (this->InteractionPerComponentCoefficientIndex[i]), sizeof(unsigned short) * TmpNbr);
	}
      RealInteractionCoefficients.WriteArray(File);
      if (this->HaveComplexInteractions)
	ComplexInteractionCoefficients.WriteArray(File);
      File.close();
      return true;
    }
  else
    {
      return false;
    }
}

// load precalculations from a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be read
// return value = true if no error occurs

bool AbstractSUNSpinOnLatticeHamiltonian::LoadPrecalculation (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  int Tmp, TmpNbr;
  File.read((char*) &(Tmp), sizeof(int));
  if (Tmp != this->Spins->GetHilbertSpaceDimension())
    {
      File.close();
      return false;
    }
  File.read((char*) &(this->HaveComplexInteractions), sizeof(bool));
  File.read((char*) &(this->FastMultiplicationStep), sizeof(int));
  Tmp /= this->FastMultiplicationStep;
  if ((Tmp * this->FastMultiplicationStep) != this->Spins->GetHilbertSpaceDimension())
    ++Tmp;
  this->NbrRealInteractionPerComponent = new unsigned short [Tmp];
  File.read((char*) this->NbrRealInteractionPerComponent, sizeof(unsigned short) * Tmp);
  if (this->HaveComplexInteractions)
    {
      this->NbrComplexInteractionPerComponent = new unsigned short [Tmp];
      File.read((char*) this->NbrComplexInteractionPerComponent, sizeof(unsigned short) * Tmp);
    }
  this->InteractionPerComponentIndex = new int* [Tmp];
  this->InteractionPerComponentCoefficientIndex = new unsigned short* [Tmp];
  for (int i = 0; i < Tmp; ++i)
    {
      TmpNbr = this->NbrRealInteractionPerComponent[i];
      if (this->HaveComplexInteractions)
	TmpNbr += this->NbrComplexInteractionPerComponent[i];
      this->InteractionPerComponentIndex[i] = new int [TmpNbr];
      File.read((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * TmpNbr);	  
    }
  for (int i = 0; i < Tmp; ++i)
    {
      TmpNbr = this->NbrRealInteractionPerComponent[i];
      if (this->HaveComplexInteractions)
	TmpNbr += this->NbrComplexInteractionPerComponent[i];
      this->InteractionPerComponentCoefficientIndex[i]=new unsigned short[TmpNbr];
      File.read((char*) (this->InteractionPerComponentCoefficientIndex[i]), sizeof(unsigned short) * TmpNbr);
    }
  RealInteractionCoefficients.ReadArray(File);
  if (this->HaveComplexInteractions)
    ComplexInteractionCoefficients.ReadArray(File);
  File.close();
  this->FastMultiplicationFlag = true;
  return true;
}

