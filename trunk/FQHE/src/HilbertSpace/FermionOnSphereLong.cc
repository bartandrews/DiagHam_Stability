////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of fermions on sphere that allow LzMax up to           //
//                  127 (for systems with 128 bit integer support)            //
//               or 63 (on 32 bit systems without 128 bit integer support)    //
//                                                                            //
//                        last modification : 13/09/2007                      //
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
#include "HilbertSpace/FermionOnSphereLong.h"
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

FermionOnSphereLong::FermionOnSphereLong()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// memory = amount of memory granted for precalculations

FermionOnSphereLong::FermionOnSphereLong (int nbrFermions, int totalLz, int lzMax, unsigned long memory)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
#ifdef __128_BIT_LONGLONG__
  this->InvertShift = 64 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 32 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  this->Flag.Initialize();
  this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrFermions, this->LzMax, this->LzMax, (this->TotalLz + this->NbrFermions * this->LzMax) >> 1, 0);
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
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

// constructor using an external array for state description
// 
// nbrFermions = number of fermions
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a fermion
// stateDescription = array that gives each state description (data are not duplicated)
// hilbertSpaceDimension = Hilbert space dimension
// memory = amount of memory granted for precalculations

FermionOnSphereLong::FermionOnSphereLong (int nbrFermions, int totalLz, int lzMax, ULONGLONG* stateDescription, long hilbertSpaceDimension, unsigned long memory)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
#ifdef __128_BIT_LONGLONG__
  this->InvertShift = 64 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 32 - ((this->LzMax + 1 ) >> 1);
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
	  while ((this->StateDescription[i] >> TmpLzMax) == ((ULONGLONG) 0x0ul))
	    --TmpLzMax;
	  this->StateLzMax[i] = TmpLzMax;
	}
    }
  else
    {
      this->StateDescription[0] = ((ULONGLONG) 0x0ul); 
      this->StateLzMax[0] = 0;
    }
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(ULONGLONG) + sizeof(int));
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

FermionOnSphereLong::FermionOnSphereLong(const FermionOnSphereLong& fermions)
{
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnSphereLong::~FermionOnSphereLong ()
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
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereLong& FermionOnSphereLong::operator = (const FermionOnSphereLong& fermions)
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
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereLong::Clone()
{
  return new FermionOnSphereLong(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnSphereLong::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnSphereLong::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphereLong::ExtractSubspace (AbstractQuantumNumber& q, 
							SubspaceSpaceConverter& converter)
{
  return 0;
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnSphereLong::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->TargetSpace = (FermionOnSphereLong*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int FermionOnSphereLong::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
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

int FermionOnSphereLong::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  ULONGLONG State = this->StateDescription[index];
  if ((n1 > StateLzMax) || (n2 > StateLzMax) || ((State & (((ULONGLONG) (0x1)) << n1)) == ((ULONGLONG) 0)) 
      || ((State & (((ULONGLONG) (0x1)) << n2)) == ((ULONGLONG) 0)) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = StateLzMax;
  ULONGLONG TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 80))  & this->SignLookUpTableMask[n2 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  TmpState &= ~(((ULONGLONG) (0x1)) << n2);
  if (NewLzMax == n2)
    while ((TmpState >> NewLzMax) == 0)
      --NewLzMax;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 80))  & this->SignLookUpTableMask[n1 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  TmpState &= ~(((ULONGLONG) (0x1)) << n1);
  if (NewLzMax == n1)
    while ((TmpState >> NewLzMax) == 0)
      --NewLzMax;
  if ((TmpState & (((ULONGLONG) (0x1)) << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewLzMax)
    {
      NewLzMax = m2;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80))  & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) (0x1)) << m2);
  if ((TmpState & (((ULONGLONG) (0x1)) << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewLzMax)
    {
      NewLzMax = m1;
    }
  else
    {      
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80))  & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) (0x1)) << m1);
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

int FermionOnSphereLong::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  ULONGLONG State = this->StateDescription[index];
  --nbrIndices;
  for (int i = 0; i < nbrIndices; ++i)
    {
      if ((n[i] > StateLzMax) || ((State & (((ULONGLONG) (0x1)) << n[i])) == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      for (int j = i + 1; j <= nbrIndices; ++j)
	if ((n[i] == n[j]) || (m[i] == m[j]))
	  {
	    coefficient = 0.0;
	    return this->HilbertSpaceDimension; 	    
	  }
    }
  if (n[nbrIndices] > StateLzMax)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }

  int NewLzMax = StateLzMax;
  ULONGLONG TmpState = State;

  int Index;
  coefficient = 1.0;
  for (int i = nbrIndices; i >= 0; --i)
    {
      Index = n[i];
      coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      TmpState &= ~(((ULONGLONG) (0x1)) << Index);
      if (NewLzMax == Index)
	while ((TmpState >> NewLzMax) == 0)
	  --NewLzMax;
    }
  for (int i = nbrIndices; i >= 0; --i)
    {
      Index = m[i];
      if ((TmpState & (((ULONGLONG) (0x1)) << Index))!= 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (((ULONGLONG) (0x1)) << Index);
    }
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereLong::ProdA (int index, int* n, int nbrIndices)
{
  this->ProdALzMax = this->StateLzMax[index];
  this->ProdATemporaryState = this->StateDescription[index];
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = n[i];
      if ((this->ProdATemporaryState & (((ULONGLONG) (0x1)) << Index)) == 0)
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#ifdef __128_BIT_LONGLONG__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 64)) & this->SignLookUpTableMask[Index + 64]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 80))  & this->SignLookUpTableMask[Index + 80]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 96)) & this->SignLookUpTableMask[Index + 96]];
      Coefficient *= this->SignLookUpTable[( this->ProdATemporaryState >> (Index + 112)) & this->SignLookUpTableMask[Index + 112]];
#endif
      this->ProdATemporaryState &= ~(((ULONGLONG) (0x1)) << Index);
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;

  return Coefficient;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereLong::AA (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];

  if (((this->ProdATemporaryState & (((ULONGLONG) (0x1)) << n1)) == ((ULONGLONG) 0)) 
      || ((this->ProdATemporaryState & (((ULONGLONG) (0x1)) << n2)) == ((ULONGLONG) 0)) || (n1 == n2))
    return 0.0;

  this->ProdALzMax = this->StateLzMax[index];

  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 80))  & this->SignLookUpTableMask[n2 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  Coefficient *= this->SignLookUpTable[( this->ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) (0x1)) << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 80))  & this->SignLookUpTableMask[n1 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) (0x1)) << n1);

  if (this->ProdATemporaryState == ((ULONGLONG) 0x0ul))
    {
      this->ProdALzMax = 0;
      return Coefficient;      
    }

  while ((this->ProdATemporaryState >> this->ProdALzMax) == ((ULONGLONG) 0))
    --this->ProdALzMax;
  return Coefficient;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereLong::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  ULONGLONG TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = m[i];
      if ((TmpState & (((ULONGLONG) (0x1)) << Index)) != ((ULONGLONG) 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#ifdef __128_BIT_LONGLONG__
	  coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 64)) & this->SignLookUpTableMask[Index + 64]];
	  coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 80))  & this->SignLookUpTableMask[Index + 80]];
	  coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 96)) & this->SignLookUpTableMask[Index + 96]];
	  coefficient *= this->SignLookUpTable[( this->ProdATemporaryState >> (Index + 112)) & this->SignLookUpTableMask[Index + 112]];
#endif
	}
      TmpState |= (((ULONGLONG) (0x1)) << Index);
    }
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereLong::AdAd (int m1, int m2, double& coefficient)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  if ((TmpState & (((ULONGLONG) (0x1)) << m2))!= ((ULONGLONG) 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80))  & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) (0x1)) << m2);
  if ((TmpState & (((ULONGLONG) (0x1)) << m1))!= ((ULONGLONG) 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {      
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80))  & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) (0x1)) << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereLong::AdA (int index, int m)
{
  if ((this->StateDescription[index] & (((ULONGLONG) (0x1)) << m)) != ((ULONGLONG) 0x0))
    return 1.0;
  else
    return 0.0;
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereLong::AdA (long index, int m)
{
  if ((this->StateDescription[index] & (((ULONGLONG) (0x1)) << m)) != ((ULONGLONG) 0x0))
    return 1.0;
  else
    return 0.0;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
  
int FermionOnSphereLong::AdA (int index, int m, int n, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  ULONGLONG State = this->StateDescription[index];
  if ((n > StateLzMax) || ((State & (((ULONGLONG) (0x1)) << n)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = StateLzMax;
  ULONGLONG TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 64)) & this->SignLookUpTableMask[n + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 80))  & this->SignLookUpTableMask[n + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 96)) & this->SignLookUpTableMask[n + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 112)) & this->SignLookUpTableMask[n + 112]];
#endif
  TmpState &= ~(((ULONGLONG) (0x1)) << n);
  if (NewLzMax == n)
    while ((TmpState >> NewLzMax) == 0)
      --NewLzMax;
  if ((TmpState & (((ULONGLONG) (0x1)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLzMax)
    {
      NewLzMax = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 64)) & this->SignLookUpTableMask[m + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 80))  & this->SignLookUpTableMask[m + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 96)) & this->SignLookUpTableMask[m + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 112)) & this->SignLookUpTableMask[m + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) (0x1)) << m);
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereLong::FindStateIndex(ULONGLONG stateDescription, int lzmax)
{
  ULONGLONG Tmp = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][Tmp];
  long PosMax = this->LookUpTable[lzmax][Tmp + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  ULONGLONG CurrentState = this->StateDescription[PosMid];
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
    return PosMin;
}

// find state index from an array of occupied orbitals
//
// stateDescription = array describing the state (stored as k1,k2,k3,...)
// return value = corresponding index, -1 if an error occured

int FermionOnSphereLong::FindStateIndex(int* stateDescription)
{
  ULONGLONG TmpState = ((ULONGLONG) 0x0ul);
  for (int i = 0; i < this->NbrFermions; ++i)
    TmpState |= ((ULONGLONG) 0x1ul) << (stateDescription[i]);
  int TmpLzMax = this->LzMax;
  while ((TmpState >> TmpLzMax) == ((ULONGLONG) 0x0ul))
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
}

// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the  component
int FermionOnSphereLong::GetLzValue(int j)
{
  return this->TotalLz;
}


// get the list of occupied orbitals in a given state
//
// state = ID of the state
// orbitals = list of orbitals to be filled

void FermionOnSphereLong::GetOccupied(int state, int* orbitals)
{
  ULONGLONG TmpState = this->StateDescription[state];
  int i = 0;
  for (int l = 0; l < this->NbrLzValue; ++l)
    {
      if (((TmpState >> l) & ((ULONGLONG) 0x1l)) == ((ULONGLONG) 0x1))
	orbitals[i++] = l;
    }
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereLong::PrintState (ostream& Str, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  for (int i = 0; i < this->NbrLzValue; ++i)
    if (((TmpState >> i) & ((ULONGLONG) 0x1)) == ((ULONGLONG) 0x1))
      Str << "1 ";
    else
      Str << "0 ";
//  Str << " key = " << this->Keys[state] << " lzmax position = " << this->LzMaxPosition[Max * (this->NbrFermions + 1) + TmpState[Max]]
//  Str << " position = " << this->FindStateIndex(TmpState, this->StateLzMax[state]);
//  if (state !=  this->FindStateIndex(TmpState, this->StateLzMax[state]))
//    Str << " error! ";
  return Str;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereLong::PrintStateMonomial (ostream& Str, long state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  Str << "[";
  int i = this->LzMax;
  while (((TmpState >> i) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
    --i;
  Str << i;
  --i;
  for (; i >=0; --i)
    if (((TmpState >> i) & ((ULONGLONG) 0x1ul)) != ((ULONGLONG) 0x0ul))
      Str << "," << i;
  Str << "]";
  return Str;
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// currentLzMax = momentum maximum value for fermions that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int FermionOnSphereLong::GenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int pos)
{
  if ((nbrFermions == 0) || (totalLz < 0) || (currentLzMax < (nbrFermions - 1)))
    return pos;
  int LzTotalMax = ((2 * currentLzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (totalLz > LzTotalMax)
    return pos;
  if ((nbrFermions == 1) && (currentLzMax >= totalLz))
    {
      this->StateDescription[pos] = ((ULONGLONG) 0x1) << totalLz;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }
  if (LzTotalMax == totalLz)
    {
      ULONGLONG Mask = 0;
      for (int i = currentLzMax - nbrFermions + 1; i <= currentLzMax; ++i)
	Mask |= (((ULONGLONG) 1) << i);
      this->StateDescription[pos] = Mask;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }

  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, lzMax, ReducedCurrentLzMax, totalLz - currentLzMax, pos);
  ULONGLONG Mask = ((ULONGLONG) 1) << currentLzMax;
  for (int i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (lzMax == currentLzMax)
    return this->GenerateStates(nbrFermions, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, TmpPos);
  else
    return this->GenerateStates(nbrFermions, lzMax, ReducedCurrentLzMax, totalLz, TmpPos);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereLong::GenerateLookUpTable(unsigned long memory)
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
  ULONGLONG CurrentLookUpTableValue = this->LookUpTableMemorySize;
  ULONGLONG TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  /*	  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
	    cout << TmpLookUpTable[j] << " ";
	    cout << endl << "-------------------------------------------" << endl;*/
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
//	      CurrentLookUpTableValue = TmpLookUpTableValue;
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;
  /*  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
    cout << TmpLookUpTable[j] << " ";
    cout << endl << "-------------------------------------------" << endl;*/
  
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
#ifdef __128_BIT_LONGLONG__
  this->SignLookUpTableMask = new ULONGLONG [256];
  for (int i = 0; i < 112; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0xffff;
  for (int i = 112; i < 128; ++i)
    this->SignLookUpTableMask[i] = ((ULONGLONG) 0xffff) >> (i - 112);
  for (int i = 128; i < 256; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0;  
#else
  this->SignLookUpTableMask = new ULONGLONG [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((ULONGLONG) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0;
#endif
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereLong::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  cout << "test " << nbrFermions << " " << lzMax << " " << ((totalLz + (nbrFermions * lzMax)) >> 1) << " " << totalLz << endl;
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + (nbrFermions * lzMax)) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long FermionOnSphereLong::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)))
    return 0l;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return 0l;
  if ((nbrFermions == 1) && (lzMax >= totalLz))
    return 1l;
  if (LzTotalMax == totalLz)
    return 1l;
  return  (this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz));
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphereLong::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
					       int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
#ifdef __LAPACK__
  ComplexLapackDeterminant Slater(this->NbrFermions);
#else
  ComplexMatrix Slater(this->NbrFermions, this->NbrFermions);
#endif
  ComplexMatrix Functions(this->LzMax + 1, this->NbrFermions);
  RealVector TmpCoordinates(2);
  int* Indices = new int [this->NbrFermions];
  int Pos;
  int Lz;
  for (int j = 0; j < this->NbrFermions; ++j)
    {
      TmpCoordinates[0] = position[j << 1];
      TmpCoordinates[1] = position[1 + (j << 1)];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	  Functions[j].Re(i) = Tmp.Re;
	  Functions[j].Im(i) = Tmp.Im;
	}
    }
  double Factor = 1.0;
  for (int i = 2; i <= this->NbrFermions; ++i)
    Factor *= (double) i;
  Factor = 1.0 / sqrt(Factor);
  ULONGLONG TmpStateDescription;
  int LastComponent = firstComponent + nbrComponent;
  for (int k = firstComponent; k < LastComponent; ++k)
    {
      Pos = 0;
      Lz = 0;
      TmpStateDescription = this->StateDescription[k];
      while (Pos < this->NbrFermions)
	{
	  if ((TmpStateDescription & ((ULONGLONG) 1)) == ((ULONGLONG) 1))
	    {
	      Indices[Pos] = Lz;
	      ++Pos;
	    }
	  ++Lz;
	  TmpStateDescription >>= 1;
	}
      for (int i = 0; i < this->NbrFermions; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < this->NbrFermions; ++j)
	    {
#ifdef __LAPACK__
	      Slater.SetMatrixElement(i,j,TmpColum2.Re(Indices[j]), TmpColum2.Im(Indices[j]));
#else
	      Slater[i].Re(j) = TmpColum2.Re(Indices[j]);
	      Slater[i].Im(j) = TmpColum2.Im(Indices[j]);
#endif
	    }
	}
      //cout << Slater << endl;

      // can calculate with lapack for a regular ComplexMatrix by Complex SlaterDet = Slater.LapackDeterminant();

      Complex SlaterDet = Slater.Determinant();
      Value += SlaterDet * (state[k] * Factor);
    }
  delete[] Indices;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereLong::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}

  
// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  FermionOnSphereLong::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState)
{  
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrFermionSector == 0))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrFermionSector == this->NbrFermions))
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * groundState[j]);
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  ULONGLONG TmpMask = (((((ULONGLONG) 0x1ul) << (this->LzMax + 2)) - ((ULONGLONG) 0x1ul)) >> subsytemSize) << subsytemSize;
  ULONGLONG TmpSubsystemMask = (((ULONGLONG) 0x1ul) << subsytemSize) - ((ULONGLONG) 0x1ul);
  int TmpIndex;
  int ShiftedTotalLz = (this->TotalLz + this->NbrFermions * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrFermionSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int TmpStateMaxLz = ShiftedLzComplementarySector - (((NbrFermionsComplementarySector - 2 + (subsytemSize << 1)) * (NbrFermionsComplementarySector - 1)) >> 1);
  int MinIndex = 0;
  int MaxIndex = this->HilbertSpaceDimension - 1;
  if ((NbrFermionsComplementarySector > 0) && ((NbrFermionsComplementarySector + subsytemSize - 2) > this->StateLzMax[MaxIndex]))
    MaxIndex = this->LookUpTable[NbrFermionsComplementarySector + subsytemSize - 2][0];
  if ((TmpStateMaxLz < this->StateLzMax[0]) && ((TmpStateMaxLz + 1) >  this->StateLzMax[MaxIndex]) && (TmpStateMaxLz >= subsytemSize))
    MinIndex = this->LookUpTable[TmpStateMaxLz + 1][0];
  
  ULONGLONG TmpComplementarySubsystem;
  int TmpNbrFermions;
  int TmpTotalLz;
  int TmpNbrOne[] = {  
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
  int TmpSumOccupation[] = {
    0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6,
    4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 10, 10,
    5, 5, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 10, 10, 11, 11,
    9, 9, 10, 10, 11, 11, 12, 12, 12, 12, 13, 13, 14, 14, 15, 15,
    6, 6, 7, 7, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 12, 12,
    10, 10, 11, 11, 12, 12, 13, 13, 13, 13, 14, 14, 15, 15, 16, 16,
    11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17,
    15, 15, 16, 16, 17, 17, 18, 18, 18, 18, 19, 19, 20, 20, 21, 21,
    7, 7, 8, 8, 9, 9, 10, 10, 10, 10, 11, 11, 12, 12, 13, 13,
    11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17,
    12, 12, 13, 13, 14, 14, 15, 15, 15, 15, 16, 16, 17, 17, 18, 18,
    16, 16, 17, 17, 18, 18, 19, 19, 19, 19, 20, 20, 21, 21, 22, 22,
    13, 13, 14, 14, 15, 15, 16, 16, 16, 16, 17, 17, 18, 18, 19, 19,
    17, 17, 18, 18, 19, 19, 20, 20, 20, 20, 21, 21, 22, 22, 23, 23,
    18, 18, 19, 19, 20, 20, 21, 21, 21, 21, 22, 22, 23, 23, 24, 24,
    22, 22, 23, 23, 24, 24, 25, 25, 25, 25, 26, 26, 27, 27, 28, 28};
  int TmpPartialNbrOne;
  if (nbrFermionSector <= 1)
    {
      ULONGLONG Key = ((ULONGLONG) 0x0ul);
      if (nbrFermionSector == 1)
	Key = ((ULONGLONG) 0x1ul) << ShiftedLzSector;
      double TmpValue = 0.0;
      while (MinIndex <= MaxIndex)
	{
	  if ((this->StateDescription[MinIndex] & TmpSubsystemMask) == Key)
	    TmpValue += groundState[MinIndex] * groundState[MinIndex];	    
	  ++MinIndex;
	}
      RealSymmetricMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  else
    {
      FermionOnSphereLong TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
      RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      long TmpNbrNonZeroElements = 0;
      while (MinIndex <= MaxIndex)
	{
	  TmpIndex = MinIndex;
	  TmpComplementarySubsystem = this->StateDescription[TmpIndex] & TmpMask;
	  ++TmpIndex;
	  while ((TmpIndex <= MaxIndex) && ((this->StateDescription[TmpIndex] & TmpMask) == TmpComplementarySubsystem))
	    ++TmpIndex;
	  TmpPartialNbrOne = TmpNbrOne[TmpComplementarySubsystem & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions = TmpPartialNbrOne;
	  TmpTotalLz = TmpSumOccupation[TmpComplementarySubsystem & ((ULONGLONG) 0xfful)];
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 8) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 8) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne << 3;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 16) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 16) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne << 4;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 24) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 24) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne * 24;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 32) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 32) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne << 5;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 40) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 40) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne * 40;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 48) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 48) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne * 48;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 56) & ((ULONGLONG) 0xfful)];      
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 56) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne * 56;
#ifdef __128_BIT_LONGLONG__
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 64) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 64) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne << 6;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 72) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 72) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne * 72;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 80) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 80) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne * 80;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 88) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 88) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne * 88;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 96) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 96) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne * 96;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 104) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 104) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne * 104;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 112) & ((ULONGLONG) 0xfful)];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 112) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne * 112;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 120) & ((ULONGLONG) 0xfful)];      
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 120) & ((ULONGLONG) 0xfful)];
	  TmpTotalLz += TmpPartialNbrOne * 120;
#endif
	  if ((TmpNbrFermions == NbrFermionsComplementarySector) && (ShiftedLzComplementarySector == TmpTotalLz))
	    {
	      int Pos = 0;
	      for (int i = MinIndex; i < TmpIndex; ++i)
		{
		  ULONGLONG TmpState = this->StateDescription[i] & TmpSubsystemMask;
		  int TmpLzMax = subsytemSize - 1;
		  while ((TmpState & (0x1ul << TmpLzMax)) == 0x0ul)
		    --TmpLzMax;
		  TmpStatePosition[Pos] = TmpDestinationHilbertSpace.FindStateIndex(TmpState, TmpLzMax);
		  ++Pos;
		}
	      int Pos2;
	      Pos = 0;
	      int Pos3;
	      double TmpValue;
	      for (int i = MinIndex; i < TmpIndex; ++i)
		{
		  Pos2 = 0;
		  Pos3 = TmpStatePosition[Pos];
		  TmpValue = groundState[i];
		  for (int j = MinIndex; j < TmpIndex; ++j)
		    {
		      if (Pos3 <=  TmpStatePosition[Pos2])
			{
			  TmpDensityMatrix.AddToMatrixElement(Pos3, TmpStatePosition[Pos2], TmpValue * groundState[j]);
			  ++TmpNbrNonZeroElements;
			}
		      ++Pos2;
		    }
		  ++Pos;
		}
	    }
	  MinIndex = TmpIndex;
	}
      delete[] TmpStatePosition;
      if (TmpNbrNonZeroElements > 0)	
	return TmpDensityMatrix;
      else
	{
	  RealSymmetricMatrix TmpDensityMatrixZero;
	  return TmpDensityMatrixZero;
	}
    }
}

// compute particule-hole symmetric state from a given state
//
// state = vector corresponding to the state to symmetrize
// holeBasis = n-body basis on which the symmetrized state has to be expressed

RealVector FermionOnSphereLong::ParticleHoleSymmetrize (RealVector& state, FermionOnSphereLong& holeBasis)
{
  RealVector TmpVector(holeBasis.HilbertSpaceDimension, true);
  ULONGLONG TmpMask = (0x1ul << (this->LzMax + 1)) - 1;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      ULONGLONG TmpState = (~this->StateDescription[i]) & TmpMask;
      int TmpLzMax = this->LzMax;
      while ((TmpState & (0x1ul << TmpLzMax)) == 0x0l)
	--TmpLzMax;
      TmpVector[holeBasis.FindStateIndex(TmpState, TmpLzMax)] = state[i];
    }
  return TmpVector;
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& FermionOnSphereLong::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  double* SqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  if (reference >= 0l)
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0 / state[reference];
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	  InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
	}
      ULONGLONG TmpState = this->StateDescription[reference];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & ((ULONGLONG) 1ul)) != ((ULONGLONG) 0ul))
	  TmpMonomialReference[Index++] = j;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & ((ULONGLONG) 1ul)) != ((ULONGLONG) 0ul))
	      TmpMonomial[Index++] = j;
	  int Index1 = 0;
	  int Index2 = 0;
	  double Coefficient = Factor;
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	    {
	      while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      state[reference] = 1.0;
      delete[] InvSqrtCoefficients;
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
   }
  else
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0;
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	  InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
	}
      ULONGLONG TmpState = this->StateDescription[0l];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & 1ul) != 0ul)
	  TmpMonomialReference[Index++] = j;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & ((ULONGLONG) 1ul)) != ((ULONGLONG) 0ul))
	      TmpMonomial[Index++] = j;
	  int Index1 = 0;
	  int Index2 = 0;
	  double Coefficient = Factor;
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	    {
	      while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      delete[] InvSqrtCoefficients;
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
   }
  delete[] SqrtCoefficients;
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& FermionOnSphereLong::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  BinomialCoefficients Binomials(this->LzMax);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  if (reference >= 0l)
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0 / state[reference];
      state[reference] = 1.0;
      double* SqrtCoefficients = new double [this->LzMax + 1];
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	  SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
	}
      ULONGLONG TmpState = this->StateDescription[reference];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & 1ul) != 0ul)
	  TmpMonomialReference[Index++] = j;
      for (int i = 1; i < this->HilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & ((ULONGLONG) 1ul)) != ((ULONGLONG) 0ul))
	      TmpMonomial[Index++] = j;
	  int Index1 = 0;
	  int Index2 = 0;
	  double Coefficient = Factor;
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	    {
	      while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
      delete[] InvSqrtCoefficients;
      state /= state.Norm();
    }
  else
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0;
      double* SqrtCoefficients = new double [this->LzMax + 1];
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	  SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
	}
      ULONGLONG TmpState = this->StateDescription[0l];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & ((ULONGLONG) 1ul)) != ((ULONGLONG) 0ul))
	  TmpMonomialReference[Index++] = j;
      for (int i = 1; i < this->HilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & ((ULONGLONG) 1ul)) != ((ULONGLONG) 0ul))
	      TmpMonomial[Index++] = j;
	  int Index1 = 0;
	  int Index2 = 0;
	  double Coefficient = Factor;
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	    {
	      while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
      delete[] InvSqrtCoefficients;
    }
  delete[] SqrtCoefficients;
  return state;
}

