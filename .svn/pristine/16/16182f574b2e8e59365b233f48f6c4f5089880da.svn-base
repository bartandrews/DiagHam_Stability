////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system size such that            //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                     and forbidding orbital multiple occupancy              //
//                                                                            //
//                        last modification : 12/10/2016                      //
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
#include "HilbertSpace/BosonOnSphereShortHardcore.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h" 

#include <math.h>
#include <stdlib.h>
#include <fstream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constuctor
//

BosonOnSphereShortHardcore::BosonOnSphereShortHardcore()
{
  this->LookUpTableShift = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// memory = amount of memory granted for precalculations

BosonOnSphereShortHardcore::BosonOnSphereShortHardcore (int nbrFermions, int totalLz, int lzMax, unsigned long memory)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;
  if (this->NbrFermions > 0)
    this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  else
    this->LargeHilbertSpaceDimension = 1l;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  if (this->NbrFermions > 0)
    {
      this->GenerateStates(this->NbrFermions, this->LzMax, this->LzMax, (this->TotalLz + this->NbrFermions * this->LzMax) >> 1, 0);
      if ((this->StateDescription[0l] >> this->LzMax) == 0x0ul)
	{
	  int TmpLzMax = this->LzMax;
	  for  (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      while ((this->StateDescription[i] >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      this->StateLzMax[i] = TmpLzMax;
	    }
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
  UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
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

// constructor using an external array for state description
// 
// nbrFermions = number of fermions
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a fermion
// stateDescription = array that gives each state description (data are not duplicated)
// hilbertSpaceDimension = Hilbert space dimension
// memory = amount of memory granted for precalculations

BosonOnSphereShortHardcore::BosonOnSphereShortHardcore (int nbrFermions, int totalLz, int lzMax, unsigned long* stateDescription, long hilbertSpaceDimension, unsigned long memory)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;
  this->LargeHilbertSpaceDimension = hilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateLzMax = new int [this->LargeHilbertSpaceDimension];
  this->StateDescription = stateDescription;
  if (this->NbrFermions > 0)
    {
      int TmpLzMax = this->LzMax;
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  while ((this->StateDescription[i] >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  this->StateLzMax[i] = TmpLzMax;
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
  UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
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

BosonOnSphereShortHardcore::BosonOnSphereShortHardcore(const BosonOnSphereShortHardcore& fermions)
{
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->ShiftedTotalLz = fermions.ShiftedTotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
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
}

// copy constructor, preversing only some specific states 
//
// fermions = reference on the hilbert space to copy to copy
// nbrPreservedStates = number of preserved states
// preservedStates = array of flags that indicates if the corresponding state has to be preserved 
//                   (dimension of the array should the one of the original Hilbert space)

BosonOnSphereShortHardcore::BosonOnSphereShortHardcore(const BosonOnSphereShortHardcore& fermions, long nbrPreservedStates, bool* preservedStates)
{
  this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->ShiftedTotalLz = fermions.ShiftedTotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->LargeHilbertSpaceDimension = nbrPreservedStates;
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateLzMax = new int [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < fermions.LargeHilbertSpaceDimension; ++i)
    {
      if (preservedStates[i] == true)
	{
	  this->StateDescription[this->LargeHilbertSpaceDimension] =  fermions.StateDescription[i];
	  this->StateLzMax[this->LargeHilbertSpaceDimension] = fermions.StateLzMax[i];
	  ++this->LargeHilbertSpaceDimension;
	}
    }
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable();
#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
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
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
}

// destructor
//

BosonOnSphereShortHardcore::~BosonOnSphereShortHardcore ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereShortHardcore& BosonOnSphereShortHardcore::operator = (const BosonOnSphereShortHardcore& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
  if (this->TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->InitializeWaveFunctionEvaluation();
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereShortHardcore::Clone()
{
  return new BosonOnSphereShortHardcore(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnSphereShortHardcore::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->TargetSpace = (BosonOnSphereShortHardcore*) targetSpace;
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShortHardcore::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  if (((TmpState & (0x1ul << n1)) == 0x0ul) || ((TmpState & (0x1ul << n2)) == 0x0ul) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  TmpState &= ~(0x1ul << n2);
  TmpState &= ~(0x1ul << n1);
  if ((TmpState & (0x1ul << m2)) != 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  TmpState |= (0x1ul << m2);
  if ((TmpState & (0x1ul << m1)) != 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  TmpState |= (0x1ul << m1);
  int NewLzMax = this->LzMax;
  while ((TmpState >> NewLzMax) == 0x0ul)
    --NewLzMax;
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShortHardcore::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  --nbrIndices;
  for (int i = 0; i < nbrIndices; ++i)
    {
      if ((TmpState & (0x1ul << n[i])) == 0x0ul)
	{
	  coefficient = 0.0;
	  return this->TargetSpace->HilbertSpaceDimension;
	}
      for (int j = i + 1; j <= nbrIndices; ++j)
	if ((n[i] == n[j]) || (m[i] == m[j]))
	  {
	    coefficient = 0.0;
	    return this->TargetSpace->HilbertSpaceDimension; 	    
	  }
    }
  if ((TmpState & (0x1ul << n[nbrIndices])) == 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }

  coefficient = 1.0;
  for (int i = nbrIndices; i >= 0; --i)
    {
      TmpState &= ~(0x1ul << n[i]);
    }
  for (int i = nbrIndices; i >= 0; --i)
    {
      if ((TmpState & (0x1ul << m[i]))!= 0x0ul)
	{
	  coefficient = 0.0;
	  return this->TargetSpace->HilbertSpaceDimension;
	}
      TmpState |= (0x1ul << m[i]);
    }
  int NewLzMax = this->LzMax;
  while ((TmpState >> NewLzMax) == 0)
    --NewLzMax;
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereShortHardcore::ProdA (int index, int* n, int nbrIndices)
{
  this->ProdALzMax = this->LzMax;
  this->ProdATemporaryState = this->StateDescription[index];
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if ((this->ProdATemporaryState & (0x1l << n[i])) == 0)
	{
	  return 0.0;
	}
      this->ProdATemporaryState &= ~(0x1l << n[i]);
    }
  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdALzMax = 0;
      return 1.0;      
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0x0ul)
    --this->ProdALzMax;
  return 1.0;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnSphereShortHardcore::AA (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];

  if (((ProdATemporaryState & (0x1ul << n1)) == 0x0ul) 
      || ((ProdATemporaryState & (0x1ul << n2)) == 0x0ul) || (n1 == n2))
    return 0.0;

  this->ProdALzMax = this->LzMax;

  this->ProdATemporaryState &= ~(0x1ul << n2);
  this->ProdATemporaryState &= ~(0x1ul << n1);

  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdALzMax = 0;
      return 1.0;      
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0x0ul)
    --this->ProdALzMax;
  return 1.0;
}

// apply a_n1 a_n2 operator to a given state without keeping it in cache
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShortHardcore::AA (int index, int n1, int n2, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];

  if (((TmpState & (0x1ul << n1)) == 0x0ul) || ((TmpState & (0x1ul << n2)) == 0x0ul) || (n1 == n2))
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  TmpState &= ~(0x1ul << n2);
  TmpState &= ~(0x1ul << n1);
  if (TmpState == 0x0ul)
    {
      return this->TargetSpace->FindStateIndex(TmpState, 0);      
    }
  int TmpLzMax = this->LzMax;
  while ((TmpState >> TmpLzMax) == 0x0ul)
    --TmpLzMax;
  return this->TargetSpace->FindStateIndex(TmpState, TmpLzMax);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShortHardcore::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if ((TmpState & (0x1l << m[i])) != 0x0ul)
	{
	  coefficient = 0.0;
	  return this->TargetSpace->HilbertSpaceDimension;
	}
      TmpState |= (0x1l << m[i]);
    }
  coefficient = 1.0;
  int NewLzMax = this->LzMax;
  while ((TmpState >> NewLzMax) == 0x0ul)
    --NewLzMax;
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShortHardcore::AdAd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if ((TmpState & (0x1ul << m2))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  TmpState |= (0x1ul << m2);
  if ((TmpState & (0x1ul << m1))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  TmpState |= (0x1ul << m1);
  coefficient = 1.0;
  int NewLzMax = this->LzMax;
  while ((TmpState >> NewLzMax) == 0x0ul)
    --NewLzMax;
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1 a^+_m2 operator to the state 
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShortHardcore::AdAd (int index, int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  if ((TmpState & (0x1ul << m2))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  TmpState |= (0x1ul << m2);
  if ((TmpState & (0x1ul << m1))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  TmpState |= (0x1ul << m1);
  coefficient = 1.0;
  int NewLzMax = this->LzMax;
  while ((TmpState >> NewLzMax) == 0x0ul)
    --NewLzMax;
  return this->FindStateIndex(TmpState, NewLzMax);
}



// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShortHardcore::AdA (int index, int m, int n, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  if ((TmpState & (0x1ul << n)) == 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  TmpState &= ~(0x1ul << n);
  if ((TmpState & (0x1ul << m))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  TmpState |= (0x1ul << m);
  int NewLzMax = this->ProdALzMax;
  while ((TmpState >> NewLzMax) == 0x0ul)
    --NewLzMax;
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

long BosonOnSphereShortHardcore::AdA (long index, int m, int n, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  if ((TmpState & (0x1ul << n)) == 0x0ul)
    {
      coefficient = 0.0;
      return this->LargeHilbertSpaceDimension;
    }
  int NewLzMax = this->LzMax;
  TmpState &= ~(0x1ul << n);
  if (TmpState != 0x0ul)
    {
      while ((TmpState >> NewLzMax) == 0)
	--NewLzMax;
    }
  else
    {
      NewLzMax = 0;
    }
  if ((TmpState & (0x1ul << m))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->LargeHilbertSpaceDimension;
    }
  if (m > NewLzMax)
    {
      NewLzMax = m;
    }
  TmpState |= (0x1ul << m);
  coefficient = 1.0;
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}
 
// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// m = Lz value of particle to be added
// coefficient = reference on the double where the multiplicative factor has to be stored

unsigned long BosonOnSphereShortHardcore::Ad (unsigned long state, int m, double& coefficient)
{
  if ((state & (0x1ul << m)) != 0x0ul)
    {
      coefficient=0.0;
      return 0x0l;
    }
  int NewLzMax = getHighestBit(state)-1;
  coefficient = 1.0;
  state |= (0x1ul << m);
  return state;
}

// apply a_n  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator
// return value =  multiplicative factor 

double BosonOnSphereShortHardcore::A (int index, int n)
{
  this->ProdATemporaryState = this->StateDescription[index];
  if ((this->ProdATemporaryState & (0x1ul << n)) == 0x0ul)
    return 0.0;

  this->ProdALzMax = this->LzMax;

  this->ProdATemporaryState &= ~(0x1ul << n);
  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdALzMax = 0;
      return 1.0;      
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0x0ul)
    --this->ProdALzMax;
  return 1.0;
}

// apply a^+_n1  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
//
// index = index of the state on which the operator has to be applied
// m = index for annihilation operator
// return value =  multiplicative factor 

double BosonOnSphereShortHardcore::Ad (int index, int m)
{
  this->ProdATemporaryState = this->StateDescription[index];
  if ((this->ProdATemporaryState & (0x1ul << m)) != 0x0ul)
    return 0.0;
  this->ProdATemporaryState |= (0x1ul << m);
  this->ProdALzMax = this->LzMax;
  if (m > this->ProdALzMax)
    this->ProdALzMax = m;
  return 1.0;
}

// apply a_n operator to the state produced using the A or Ad method (without destroying it)
//
// n = first index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShortHardcore::A (int n, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if ((TmpState & (0x1ul << n))== 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  TmpState &= ~(0x1ul << n);
  if (TmpState == 0x0ul)
    {
      NewLzMax = 0;
    }
  else
    {
      while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;
    }
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m operator to the state produced using the A or Ad method (without destroying it)
//
// m = index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShortHardcore::Ad (int m, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
    
  if ((TmpState & (0x1ul << m))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m > NewLzMax)
    NewLzMax = m;
  TmpState |= (0x1ul << m);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a_n  operator to a given state. 
//
// index = index of the state on which the operator has to be applied
// n = index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value =  index of the resulting state 

int BosonOnSphereShortHardcore::A (int index, int n, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  if ((TmpState & (0x1ul << n))== 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  TmpState &= ~(0x1ul << n);  
  int NewLzMax;
  if (TmpState == 0x0ul)
    {
      NewLzMax = 0;
    }
  else
    {
      NewLzMax = this->LzMax;
      while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;
    }
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}


// generate look-up table for sign calculation
// 

void BosonOnSphereShortHardcore::GenerateSignLookUpTable()
{
  // look-up tables for evaluating sign when applying creation/annihilation operators
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      Count = 0;
      TmpNbr = j;
      while (TmpNbr != 0)
	{
	  if (TmpNbr & 0x1)
	    ++Count;
	  TmpNbr >>= 1;
	}
      if (Count & 1)
	this->SignLookUpTable[j] = -1.0;
      else
	this->SignLookUpTable[j] = 1.0;
    }
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#endif
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int BosonOnSphereShortHardcore::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->LargeHilbertSpaceDimension - 1l]))
    {
      return this->HilbertSpaceDimension;
    }
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereShortHardcore::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->NbrLzValue);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrLzValue)
    this->MaximumLookUpShift = this->NbrLzValue;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrLzValue];
  this->LookUpTableShift = new int [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentLzMax = this->StateLzMax[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentLzMax];
  if (CurrentLzMax < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentLzMax] = 0;
  else
    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentLzMax];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
 	  CurrentLzMax = this->StateLzMax[i];
	  TmpLookUpTable = this->LookUpTable[CurrentLzMax];
	  if (CurrentLzMax < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentLzMax] = 0;
	  else
	    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentLzMax];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      --CurrentLookUpTableValue;
    }

  // look-up tables for evaluating sign when applying creation/annihilation operators
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      this->SignLookUpTable[j] = 1.0;
    }
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#endif
}

