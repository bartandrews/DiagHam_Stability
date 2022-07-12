////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of bosons on sphere truncated to a given P level            //
//  such that LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)    //
//                                                                            //
//                        last modification : 17/10/2012                      //
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
#include "HilbertSpace/BosonOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncated.h"

#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <map>


using std::cout;
using std::endl;
using std::dec;
using std::hex;
using std::map;
using std::pair;

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson
// pLevel = truncation level
// maximumOccupation = maximum occupation for a single orbital
// referenceState = array that describes the root configuration

BosonOnSpherePTruncated::BosonOnSpherePTruncated (int nbrBosons, int& totalLz, int lzMax, int pLevel,
						  int maximumOccupation, int* referenceState)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->TotalLz = 0;
  this->PLevel = pLevel;
  this->MaximumOccupation = maximumOccupation;
  this->ReferenceStateMonomialBasis = new int [this->NbrBosons];
  this->ReferenceState = 0x0ul;
  int ShiftedLzMax = this->LzMax + this->NbrBosons - 1;
  int* FermionicReferenceState =  new int [ShiftedLzMax + 1];
  for (int i = 0; i <= ShiftedLzMax; ++i)
    FermionicReferenceState[i] = 0;   
  int TmpIndex = 0; 
  int TmpFermionIndex = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      for (int j = 0; j < referenceState[i]; ++j)
	{
	  this->ReferenceStateMonomialBasis[TmpIndex] = i;
	  FermionicReferenceState[TmpFermionIndex] = 1;
	  ++TmpFermionIndex;
	  this->TotalLz += i;
	  ++TmpIndex;
	}
      ++TmpFermionIndex;
    }
  this->TotalLz = ((this->TotalLz << 1) - (this->LzMax * this->NbrBosons));
  totalLz = this->TotalLz;
  if (nbrBosons > 0)
    this->FermionBasis = new FermionOnSpherePTruncated(nbrBosons, totalLz, ShiftedLzMax, this->PLevel, FermionicReferenceState);
  else
    this->FermionBasis = new FermionOnSpherePTruncated(nbrBosons, totalLz, lzMax, this->PLevel, FermionicReferenceState);
  this->HilbertSpaceDimension = this->FermionBasis->HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = this->FermionBasis->LargeHilbertSpaceDimension;
  if (this->MaximumOccupation < this->NbrBosons)
    {
      long NbrPreservedStates = 0l;
      bool* PreservedStates =  new bool[this->LargeHilbertSpaceDimension];
      unsigned int TmpMaximumOccupation = (unsigned int) this->MaximumOccupation;
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if (this->CheckMaximumOccupation(this->FermionBasis->StateDescription[i], TmpMaximumOccupation))
	    {
	      PreservedStates[i] = true;
	      ++NbrPreservedStates;
	    }
	  else
	    {
	      PreservedStates[i] = false;	      
	    }
	}
      FermionOnSpherePTruncated* TmpBasis = new FermionOnSpherePTruncated(*((FermionOnSpherePTruncated*) this->FermionBasis), NbrPreservedStates, PreservedStates);
      delete this->FermionBasis;
      delete[] PreservedStates;
      this->FermionBasis = TmpBasis;
      this->HilbertSpaceDimension = this->FermionBasis->HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = this->FermionBasis->LargeHilbertSpaceDimension;
    }

  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->Flag.Initialize();
  int TmpLzMax = this->LzMax;
  this->TargetSpace = this;

  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;

}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSpherePTruncated::BosonOnSpherePTruncated(const BosonOnSpherePTruncated& bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->PLevel = bosons.PLevel;
  this->MaximumOccupation = bosons.MaximumOccupation;
  this->ReferenceStateMonomialBasis = new int [this->NbrBosons];
  for (int i = 0; i < this->NbrBosons; ++i)
    this->ReferenceStateMonomialBasis[i] = bosons.ReferenceStateMonomialBasis[i];
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = (FermionOnSpherePTruncated*) bosons.FermionBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
}

// destructor
//

BosonOnSpherePTruncated::~BosonOnSpherePTruncated ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSpherePTruncated& BosonOnSpherePTruncated::operator = (const BosonOnSpherePTruncated& bosons)
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->PLevel = bosons.PLevel;
  this->MaximumOccupation = bosons.MaximumOccupation;
  this->ReferenceStateMonomialBasis = new int [this->NbrBosons];
  for (int i = 0; i < this->NbrBosons; ++i)
    this->ReferenceStateMonomialBasis[i] = bosons.ReferenceStateMonomialBasis[i];
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = (FermionOnSpherePTruncated*) bosons.FermionBasis->Clone();
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSpherePTruncated::Clone()
{
  return new BosonOnSpherePTruncated(*this);
}

// create a state from its MPS description, assuming the resulting state is real
//
// bMatrices = array that gives the B matrices 
// state = reference to vector that will contain the state description
// mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// mPSColumnIndex = column index of the MPS element that has to be evaluated
// memory = amount of memory that can be use to precompute matrix multiplications  
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void BosonOnSpherePTruncated::CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, RealVector& state, int mPSRowIndex, int mPSColumnIndex, 
							     long memory, long initialIndex, long nbrComponents)
{
  SparseRealMatrix TmpMatrix;
  long MaxIndex = initialIndex + nbrComponents;
  if ((nbrComponents == 0l) || (MaxIndex > this->LargeHilbertSpaceDimension))
    {
      MaxIndex = this->LargeHilbertSpaceDimension;
    }
  long TmpMemory = (((long) bMatrices[0].GetNbrRow()) * ((long)bMatrices[0].GetNbrColumn()));
  if (TmpMemory > (1l << 28))
    TmpMemory = 1l << 28;
  double* TmpMatrixElements = new double [TmpMemory];
  int* TmpColumnIndices = new int [TmpMemory];
  double* TmpElements = new double [bMatrices[0].GetNbrRow()];

  if (memory <= 1l)
    {
      for (long i = initialIndex; i < MaxIndex; ++i)
	{
	  if (((i - initialIndex) % 10000) == 0)
	    cout << "Completed " << (i - initialIndex) << " out of " << (MaxIndex - initialIndex) << endl; 
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
	  TmpMatrix.Copy(bMatrices[this->ProdATemporaryState[0]]);
	  for (int j = 1; j <= this->ProdATemporaryStateLzMax; ++j)
	    {
	      TmpMatrix.Multiply(bMatrices[this->ProdATemporaryState[j]], TmpMatrixElements, TmpColumnIndices, TmpElements);
	    } 
	  for (int j = this->ProdATemporaryStateLzMax + 1; j <= this->LzMax; ++j)
	    TmpMatrix.Multiply(bMatrices[0], TmpMatrixElements, TmpColumnIndices, TmpElements);
	  if (mPSRowIndex < 0)
	    state[i] = TmpMatrix.Tr();
	  else
	    TmpMatrix.GetMatrixElement(mPSRowIndex, mPSColumnIndex, state[i]);
	}
    }
  else
    {
//       int PrecalculationBlockSize = (int) memory;
//       unsigned long PrecalculationBlockMask = (0x1ul  << PrecalculationBlockSize) - 0x1ul;
//       int PrecalculationBlockLength  = this->LzMax - ((this->LzMax + 1) % PrecalculationBlockSize);
//       int RemaingOrbtals = (this->LzMax + 1) % PrecalculationBlockSize;
//       SparseComplexMatrix* TmpBlockbMatrices = new SparseComplexMatrix[1 << PrecalculationBlockSize];
//       for (unsigned long TmpState = 0x0ul; TmpState <= PrecalculationBlockMask; ++TmpState)
// 	{
// 	  unsigned long TmpStateDescription = TmpState;
// 	  TmpBlockbMatrices[TmpState].Copy(bMatrices[TmpStateDescription & 0x1ul]);
// 	  TmpStateDescription >>= 1;
// 	  for (int j = 1; j < PrecalculationBlockSize; ++j)
// 	    {
// 	      TmpBlockbMatrices[TmpState].Multiply(bMatrices[TmpStateDescription & 0x1ul], TmpMatrixElements, TmpColumnIndices, TmpElements);
// 	      TmpStateDescription >>= 1;
// 	    }      
// 	}
//       for (long i = initialIndex; i < MaxIndex; ++i)
// 	{
// 	  if (((i - initialIndex) % 10000) == 0)
// 	    cout << "Completed " << (i - initialIndex) << " out of " << (MaxIndex - initialIndex) << endl; 
// 	  unsigned long TmpStateDescription = this->StateDescription[i];
// 	  TmpMatrix.Copy(TmpBlockbMatrices[TmpStateDescription & PrecalculationBlockMask]);
// 	  TmpStateDescription >>= PrecalculationBlockSize;
// 	  int j = PrecalculationBlockSize;
// 	  for (; j < PrecalculationBlockLength; j += PrecalculationBlockSize)
// 	    {
// 	      TmpMatrix.Multiply(TmpBlockbMatrices[TmpStateDescription & PrecalculationBlockMask], TmpMatrixElements, TmpColumnIndices, TmpElements);
// 	      TmpStateDescription >>= PrecalculationBlockSize;
// 	    } 
// 	  for (; j <= this->LzMax; ++j)
// 	    {
// 	      TmpMatrix.Multiply(bMatrices[TmpStateDescription & 0x1ul], TmpMatrixElements, TmpColumnIndices, TmpElements);
// 	      TmpStateDescription >>= 1;
// 	    } 
// 	  if (mPSRowIndex < 0)
// 	    state[i] = TmpMatrix.Tr();
// 	  else
// 	    TmpMatrix.GetMatrixElement(mPSRowIndex, mPSColumnIndex, state[i]);
// 	}
    }
  delete[] TmpMatrixElements;
  delete[] TmpColumnIndices;
  delete[] TmpElements;
}

