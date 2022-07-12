////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of fermions on sphere truncated to a given P level           //
//      with  LzMax up to 127 (for systems with 128 bit integer support)      //
//           or 63 (on 32 bit systems without 128 bit integer support)        //
//                                                                            //
//                        last modification : 02/10/2012                      //
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
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/FactorialCoefficient.h" 
#include "MathTools/BinomialCoefficients.h"

#include <math.h>


using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constuctor
//

FermionOnSpherePTruncatedLong::FermionOnSpherePTruncatedLong()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// pLevel = truncation level
// referenceState = array that describes the root configuration
// memory = amount of memory granted for precalculations

FermionOnSpherePTruncatedLong::FermionOnSpherePTruncatedLong (int nbrFermions, int& totalLz, int lzMax, int pLevel,
							      int* referenceState, unsigned long memory)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->PLevel = pLevel;
  this->ReferenceStateMonomialBasis = new int [this->NbrFermions];
  int TmpIndex = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      this->ReferenceState |= ((unsigned long) (referenceState[i] & 1)) << i;
      if (referenceState[i] == 1)
	{
	  this->ReferenceStateMonomialBasis[TmpIndex] = i;
	  this->TotalLz += i;
	  ++TmpIndex;
	}
    }
  this->TotalLz = ((this->TotalLz << 1) - (this->LzMax * this->NbrFermions));
  totalLz = this->TotalLz;
#ifdef __128_BIT_LONGLONG__
  this->InvertShift = 64 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 32 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;
  if (this->NbrFermions > 0)
    this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz, 0);
  else
    this->LargeHilbertSpaceDimension = 1l;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new ULONGLONG [this->LargeHilbertSpaceDimension];
  this->StateLzMax = new int [this->LargeHilbertSpaceDimension];
  if (this->NbrFermions > 0)
    {
      long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, this->LzMax, (this->TotalLz + this->NbrFermions * this->LzMax) >> 1, 0, 0);
      if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating Hilbert space, generation = " << TmpHilbertSpaceDimension << " states, dimension = " << this->LargeHilbertSpaceDimension << endl;
	}
    }
  else
    {
      this->StateDescription[0] = 0x0ul; 
      this->StateLzMax[0] = 0;
    }
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(ULONGLONG) + sizeof(int));
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
  UsedMemory +=  (1 << this->MaximumSignLookUp) * sizeof(double);
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSpherePTruncatedLong::FermionOnSpherePTruncatedLong(const FermionOnSpherePTruncatedLong& fermions)
{
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->PLevel = fermions.PLevel;
  this->ReferenceStateMonomialBasis = new int [this->NbrFermions];
  for (int i = 0; i < this->NbrFermions; ++i)
    this->ReferenceStateMonomialBasis[i] = fermions.ReferenceStateMonomialBasis[i];
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
}

// destructor
//

FermionOnSpherePTruncatedLong::~FermionOnSpherePTruncatedLong ()
{
  delete[] this->ReferenceStateMonomialBasis;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSpherePTruncatedLong& FermionOnSpherePTruncatedLong::operator = (const FermionOnSpherePTruncatedLong& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      delete[] this->ReferenceStateMonomialBasis;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->PLevel = fermions.PLevel;
  this->ReferenceStateMonomialBasis = new int [this->NbrFermions];
  for (int i = 0; i < this->NbrFermions; ++i)
    this->ReferenceStateMonomialBasis[i] = fermions.ReferenceStateMonomialBasis[i];
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSpherePTruncatedLong::Clone()
{
  return new FermionOnSpherePTruncatedLong(*this);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// currentLzMax = momentum maximum value for fermions that are still to be placed
// totalLz = momentum total value
// level = current level for truncation
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSpherePTruncatedLong::GenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int level, long pos)
{
  if ((nbrFermions == 0) || (totalLz < 0) || (currentLzMax < (nbrFermions - 1)) || (level < 0))
    return pos;
  int LzTotalMax = ((2 * currentLzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (totalLz > LzTotalMax)
    return pos;
  if ((nbrFermions == 1) && (currentLzMax >= totalLz))
    {
      if (((level + this->ReferenceStateMonomialBasis[0] - totalLz) > this->PLevel) 
	  || ((level + this->ReferenceStateMonomialBasis[0] - totalLz) < 0))
	return pos;
      this->StateDescription[pos] = ((ULONGLONG) 0x1) << totalLz;
      this->StateLzMax[pos] = lzMax;
      return pos + 1l;
    }

  int ReducedCurrentLzMax = currentLzMax - 1;
  if ((level + this->ReferenceStateMonomialBasis[nbrFermions - 1] - currentLzMax) <= this->PLevel)
    {
      long TmpPos = this->GenerateStates(nbrFermions - 1, lzMax, ReducedCurrentLzMax, totalLz - currentLzMax, 
					 level + this->ReferenceStateMonomialBasis[nbrFermions - 1] - currentLzMax, pos);
      ULONGLONG Mask = ((ULONGLONG) 1) << currentLzMax;
      for (long i = pos; i < TmpPos; ++i)
	this->StateDescription[i] |= Mask;
      if (lzMax == currentLzMax)
	return this->GenerateStates(nbrFermions, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, level, TmpPos);
      else
	return this->GenerateStates(nbrFermions, lzMax, ReducedCurrentLzMax, totalLz, level, TmpPos);
    }
  else
    {
      if (lzMax == currentLzMax)
	return this->GenerateStates(nbrFermions, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, level, pos);
      else
	return this->GenerateStates(nbrFermions, lzMax, ReducedCurrentLzMax, totalLz, level, pos);
    }
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// level = current level for truncation
// return value = Hilbert space dimension

long FermionOnSpherePTruncatedLong::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int level)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + nbrFermions * lzMax) >> 1, level);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// level = current level for truncation
// return value = Hilbert space dimension

long FermionOnSpherePTruncatedLong::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int level)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)) || (level < 0))
    return 0l;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return 0l;
  if ((nbrFermions == 1) && (lzMax >= totalLz) && ((level + this->ReferenceStateMonomialBasis[0] - lzMax) <= this->PLevel)
      && ((level + this->ReferenceStateMonomialBasis[0] - lzMax) >= 0))
    return 1l;
  long TmpPos = 0l;
  if ((level + this->ReferenceStateMonomialBasis[nbrFermions - 1] - lzMax) <= this->PLevel)
    {
      TmpPos = this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax,
							  level + this->ReferenceStateMonomialBasis[nbrFermions - 1] - lzMax);
    }  
  return TmpPos + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, level);
}

