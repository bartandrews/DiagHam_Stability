////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with SU(3) spin              //
//                                                                            //
//                        last modification : 20/01/2008                      //
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
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <math.h>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
// 

FermionOnSphereWithSU3Spin::FermionOnSphereWithSU3Spin ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// totalTz = twice the total Tz value
// totalY = three time the total Y value
// memory = amount of memory granted for precalculations

FermionOnSphereWithSU3Spin::FermionOnSphereWithSU3Spin (int nbrFermions, int totalLz, int lzMax, int totalTz, int totalY,
							unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalY = totalY;
  this->TotalTz = totalTz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  int N1 = (2 * nbrFermions) + totalY + (3 * totalTz);
  int N2 = (2 * nbrFermions) + totalY - (3 * totalTz);
  int N3 = nbrFermions - totalY;
  if ((N1 < 0) || (N2 < 0) || (N3 < 0) || ((N1 % 6) != 0) || ((N2 % 6) != 0) || ((N3 % 3) != 0))
    this->HilbertSpaceDimension = 0;
  else
    {
      N1 /= 6;
      N2 /= 6;
      N3 /= 3;
      this->HilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, N1, N2, N3);
    }
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						       N1, N2, N3, 0l);
  if (TmpHilbertSpaceDimension != this->HilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->HilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in FermionOnSphereWithSU3Spin!" << endl;
      exit(1);
    }
  this->HilbertSpaceDimension = TmpHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
  this->GenerateLookUpTable(memory);
//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)	
//     this->PrintState(cout, i) << endl;
// #ifdef __DEBUG__
//   int UsedMemory = 0;
//   UsedMemory += this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
//   cout << "memory requested for Hilbert space = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;
//   UsedMemory = this->NbrLzValue * sizeof(int);
//   UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
//   cout << "memory requested for lookup table = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;

// #endif
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereWithSU3Spin::FermionOnSphereWithSU3Spin(const FermionOnSphereWithSU3Spin& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalTz = fermions.TotalTz;
  this->TotalY = fermions.TotalY;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
}

// destructor
//

FermionOnSphereWithSU3Spin::~FermionOnSphereWithSU3Spin ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < 2*this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSU3Spin& FermionOnSphereWithSU3Spin::operator = (const FermionOnSphereWithSU3Spin& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->TotalTz = fermions.TotalTz;
  this->TotalY = fermions.TotalY;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSU3Spin::Clone()
{
  return new FermionOnSphereWithSU3Spin(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnSphereWithSU3Spin::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnSphereWithSU3Spin::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphereWithSU3Spin::ExtractSubspace (AbstractQuantumNumber& q, 
							SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m_1 a_m_1 operator to a given state (only state 1 Tz=+1/2, Y=+1/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_1 a_m_1

double  FermionOnSphereWithSU3Spin::Ad1A1 (int index, int m)
{
  if ((this->StateDescription[index] & (0x1l << (m * 3))) != 0)
    return 1.0;
  else
    return 0.0;
}

// apply a^+_m_2 a_m_2 operator to a given state (only state 2 Tz=-1/2, Y=+1/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_2 a_m_2

double FermionOnSphereWithSU3Spin::Ad2A2 (int index, int m)
{
  if ((this->StateDescription[index] & (0x2l << (m * 3))) != 0)
    return 1.0;
  else
    return 0.0;
}
 
// apply a^+_m_3 a_m_3 operator to a given state (only state 3 Tz=0, Y=-2/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_3 a_m_3

double FermionOnSphereWithSU3Spin::Ad3A3 (int index, int m)
{
  if ((this->StateDescription[index] & (0x4l << (m * 3))) != 0)
    return 1.0;
  else
    return 0.0;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereWithSU3Spin::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  unsigned long CurrentState = stateDescription >> this->LookUpTableShift[lzmax];
  int PosMin = this->LookUpTable[lzmax][CurrentState];
  int PosMax = this->LookUpTable[lzmax][CurrentState+ 1];
  int PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->StateDescription[PosMid];
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



// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereWithSU3Spin::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << " | ";
  for (int i = this->NbrLzValue-1; i >=0 ; --i)
    {
      Tmp = ((TmpState >> (i * 3)) & ((unsigned long) 0x7ul));
      switch (Tmp)
	{
	case 0x7ul :
	  Str << "123 ";
	  break;
	case 0x6ul :
	  Str << " 23 ";
	  break;
	case 0x5ul :
	  Str << "1 3 ";
	  break;
	case 0x4ul :
	  Str << "  3 ";
	  break;
	case 0x3ul :
	  Str << "12  ";
	  break;
	case 0x2ul :
	  Str << " 2  ";
	  break;
	case 0x1ul :
	  Str << "1   ";
	  break;
	case 0x0ul :
	  Str << " 0  ";
	  break;
	}
      Str << "| ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSU3Spin::GenerateStates(int nbrFermions, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3, long pos)
{
//   if (pos == 0)
//     cout << nbrFermions << " " << lzMax << " " << totalLz << " " <<  nbrN1 << " " <<  nbrN2 << " " <<  nbrN3 << endl;
  if ((nbrFermions < 0) || (totalLz < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0))
    return pos;
//   if (pos == 0)
//     cout << "check" << endl;
  if ((nbrFermions == 0) && (totalLz == 0))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
  if ((lzMax < 0) || ((nbrN1 - 1)> lzMax) || ((nbrN2 - 1)> lzMax) || ((nbrN3 - 1)> lzMax) ||
      ((nbrFermions * lzMax - (((nbrN1 * nbrN1) + (nbrN2 * nbrN2) + (nbrN3 * nbrN3) - nbrFermions) >> 1)) < totalLz))
    return pos;

  if (nbrFermions == 1) 
    {
      if (lzMax >= totalLz)
	{
	  this->StateDescription[pos] = 0x1ul << ((totalLz * 3) + nbrN2 + (nbrN3 << 1));
	  return (pos + 1l);
	}
      else
	return pos;
    }

  long TmpPos;
  unsigned long Mask;
  if (nbrFermions >= 3)
    {
      TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (lzMax * 3), nbrN1 - 1, nbrN2 - 1, nbrN3 - 1, pos);
      Mask = 0x7ul << (lzMax * 3);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), nbrN1, nbrN2 - 1, nbrN3 - 1, pos);
  Mask = 0x6ul << (lzMax * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), nbrN1 - 1, nbrN2, nbrN3 - 1, pos);
  Mask = 0x5ul << (lzMax * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrN1, nbrN2, nbrN3 - 1, pos);
  Mask = 0x4ul << (lzMax * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), nbrN1 - 1, nbrN2 - 1, nbrN3, pos);
  Mask = 0x3ul << (lzMax * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrN1, nbrN2 - 1, nbrN3, pos);
  Mask = 0x2ul << (lzMax * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrN1 - 1, nbrN2, nbrN3, pos);
  Mask = 0x1ul << (lzMax * 3);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  return this->GenerateStates(nbrFermions, lzMax - 1, totalLz, nbrN1, nbrN2, nbrN3, pos);
};


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereWithSU3Spin::GenerateLookUpTable(unsigned long memory)
{
  // get every highest bit poisition
  unsigned long TmpPosition = this->StateDescription[0];
  int CurrentHighestBit = (this->LzMax + 1) * 3 - 1;
  while ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
    --CurrentHighestBit;  
  int MaxHighestBit = CurrentHighestBit;
  this->StateHighestBit[0] = CurrentHighestBit;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      TmpPosition = this->StateDescription[i];
      while ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
	--CurrentHighestBit;  
      this->StateHighestBit[i] = CurrentHighestBit;
   }

  // evaluate look-up table size
  memory /= (sizeof(int*) * MaxHighestBit);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > MaxHighestBit)
    this->MaximumLookUpShift = MaxHighestBit;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [MaxHighestBit + 1];
  this->LookUpTableShift = new int [MaxHighestBit + 1];
  for (int i = 0; i <= MaxHighestBit; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];

  CurrentHighestBit = this->StateHighestBit[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
  if (CurrentHighestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentHighestBit] = 0;
  else
    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentHighestBit];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentHighestBit != this->StateHighestBit[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
 	  CurrentHighestBit = this->StateHighestBit[i];
	  TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
	  if (CurrentHighestBit < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentHighestBit] = 0;
	  else
	    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentHighestBit];
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
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;

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

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// return value = Hilbert space dimension

long FermionOnSphereWithSU3Spin::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3)
{
  if ((nbrFermions < 0) || (totalLz < 0) || (lzMax < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0) || (lzMax < 0) || 
      ((nbrN1 - 1)> lzMax) || ((nbrN2 - 1)> lzMax) || ((nbrN3 - 1)> lzMax) ||
      ((nbrFermions * lzMax - (((nbrN1 * nbrN1) + (nbrN2 * nbrN2) + (nbrN3 * nbrN3) - nbrFermions) >> 1)) < totalLz))
    return 0l;
  if ((nbrFermions == 0) && (totalLz == 0))
    return 1l;
  if (nbrFermions == 1) 
    if (lzMax >= totalLz)
      return 1l;
    else
      return 0l;
  unsigned long Tmp = 0l;
  if (nbrFermions >= 3)
    {
      Tmp += (this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), nbrN1 - 1, nbrN2 - 1, nbrN3)
	      + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), nbrN1 - 1, nbrN2, nbrN3 - 1)
	      + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), nbrN1, nbrN2 - 1, nbrN3 - 1));
      if (nbrFermions == 3)
	{
	  if ((totalLz == (3 * lzMax)) && (nbrN1 == 1) && (nbrN2 == 1) && (nbrN3 == 1))
	    ++Tmp;      
	}
      else
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (lzMax * 3), nbrN1 - 1, nbrN2 - 1, nbrN3 -1);
    }
  else
    {
      if ((totalLz == (2 * lzMax)) && (nbrN1 <= 1) && (nbrN2 <= 1) && (nbrN3 <= 1))
	++Tmp;
    }
  return  (Tmp + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrN1 - 1, nbrN2, nbrN3)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrN1, nbrN2 - 1, nbrN3)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrN1, nbrN2, nbrN3 - 1)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, nbrN1, nbrN2, nbrN3));
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalTz = twice the total Tz value
// totalY = three time the total Y value
// return value = Hilbert space dimension

int FermionOnSphereWithSU3Spin::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalTz, int totalY)
{
  int N1 = (2 * nbrFermions) + totalY + (3 * totalTz);
  int N2 = (2 * nbrFermions) + totalY - (3 * totalTz);
  int N3 = nbrFermions - totalY;
  if ((N1 >= 0) && (N2 >= 0) && (N3 >= 0) && ((N1 % 6) == 0) && ((N2 % 6) == 0) && ((N3 % 3) == 0))
    {
      N1 /= 6;
      N2 /= 6;
      N3 /= 3;
      return this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + (nbrFermions * lzMax)) >> 1, N1, N2, N3);
    }
  else
    return 0;
}


// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphereWithSU3Spin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
							  int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
  int N1 = (2 * this->NbrFermions) + this->TotalY + (3 * this->TotalTz);
  int N2 = (2 * this->NbrFermions) + this->TotalY - (3 * this->TotalTz);
  int N3 = this->NbrFermions - this->TotalY;
  N1 /= 6;
  N2 /= 6;
  N3 /= 3;
  ComplexMatrix Slater1(N1, N1);
  ComplexMatrix Slater2(N2, N2);
  ComplexMatrix Slater3(N3, N3);
  ComplexMatrix Functions(this->LzMax + 1, this->NbrFermions);
  RealVector TmpCoordinates(2);
  int* Indices1 = new int [N1];
  int* Indices2 = new int [N2];
  int* Indices3 = new int [N3];
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
  int Pos1;
  int Pos2;
  int Pos3;
  int Lz;
  unsigned long TmpStateDescription;
  int LastComponent = firstComponent + nbrComponent;
  for (int k = firstComponent; k < LastComponent; ++k)
    {
      Pos1 = 0;
      Pos2 = 0;
      Pos3 = 0;
      Lz = 0;
      TmpStateDescription = this->StateDescription[k];
      while (Lz <= this->LzMax)
	{
	  if ((TmpStateDescription & 0x1l) != 0x0l)
	    {
	      Indices1[Pos1] = Lz;
	      ++Pos1;
	    }
	  if ((TmpStateDescription & 0x2l) != 0x0l)
	    {
	      Indices2[Pos2] = Lz;
	      ++Pos2;
	    }
	  if ((TmpStateDescription & 0x4l) != 0x0l)
	    {
	      Indices3[Pos3] = Lz;
	      ++Pos3;
	    }
	  ++Lz;
	  TmpStateDescription >>= 3;
	}
      for (int i = 0; i < N1; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < N1; ++j)
	    {
#ifdef __USE_LAPACK_HERE__
	      Slater1.SetMatrixElement(i,j,TmpColum2.Re(Indices1[j]), TmpColum2.Im(Indices1[j]));
#else
	      Slater1[i].Re(j) = TmpColum2.Re(Indices1[j]);
	      Slater1[i].Im(j) = TmpColum2.Im(Indices1[j]);
#endif
	    }
	}
      for (int i = 0; i < N2; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < N2; ++j)
	    {
#ifdef __USE_LAPACK_HERE__
	      Slater2.SetMatrixElement(i,j,TmpColum2.Re(Indices2[j]), TmpColum2.Im(Indices2[j]));
#else
	      Slater2[i].Re(j) = TmpColum2.Re(Indices2[j]);
	      Slater2[i].Im(j) = TmpColum2.Im(Indices2[j]);
#endif
	    }
	}
      for (int i = 0; i < N3; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < N3; ++j)
	    {
#ifdef __USE_LAPACK_HERE__
	      Slater3.SetMatrixElement(i,j,TmpColum2.Re(Indices3[j]), TmpColum2.Im(Indices3[j]));
#else
	      Slater3[i].Re(j) = TmpColum2.Re(Indices3[j]);
	      Slater3[i].Im(j) = TmpColum2.Im(Indices3[j]);
#endif
	    }
	}
      Complex SlaterDet1 = Slater1.Determinant();
      Complex SlaterDet2 = Slater2.Determinant();
      Complex SlaterDet3 = Slater3.Determinant();
//      Value += SlaterDet1 * SlaterDet2 * SlaterDet3 * (state[k] * Factor) * this->GetStateSign(TmpStateDescription, IndicesDown);
    }
  delete[] Indices1;
  delete[] Indices2;
  delete[] Indices3;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereWithSU3Spin::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}

// create a U(1) state from an SU(3) state
//
// state = vector describing the SU(3) state
// u1Space = reference on the Hilbert space associated to the U(1) state
// return value = resulting U(1) state

RealVector FermionOnSphereWithSU3Spin::ForgeU1FromSU3(RealVector& state, FermionOnSphere& u1Space)
{
  RealVector FinalState(u1Space.GetHilbertSpaceDimension(), true);
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      unsigned long TmpState = this->StateDescription[j];
      unsigned long TmpState2 = TmpState; 
      int TmpPos = this->LzMax * 3;
      while (TmpPos >=0)
	{
	  unsigned long  TmpNbrParticles = TmpState2 & 0x1ul;
	  TmpState2 >>= 1;
	  TmpNbrParticles += TmpState2 & 0x1ul;
	  TmpState2 >>= 1;
	  TmpNbrParticles += TmpState2 & 0x1ul;
	  TmpState2 >>= 1;
	  if (TmpNbrParticles > 0x1ul)
	    TmpPos = 1;
	  TmpPos -= 3;
	}
      if (TmpPos != -2)
	{ 
	  TmpPos = 0;
	  TmpState2 = 0x0ul; 
	  while (TmpPos <= this->LzMax)
	    {
	      TmpState2 |= ((TmpState & 0x1ul) | ((TmpState & 0x2ul) >> 1) | ((TmpState & 0x4ul) >> 2)) << TmpPos;
	      TmpState >>= 3;
	      ++TmpPos;
	    }
	  while ((TmpState2 >> TmpPos) == 0x0ul)
	    --TmpPos;
	  FinalState[u1Space.FindStateIndex(TmpState2, TmpPos)] += state[j];
	}
    }
  FinalState /= FinalState.Norm();
  return FinalState;  
}

// apply a_n1_1 a_n2_1 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double FermionOnSphereWithSU3Spin::A1A1 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  n2 *= 3;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_1 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double FermionOnSphereWithSU3Spin::A1A2 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  n2 *= 3;
  ++n2;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_1 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double FermionOnSphereWithSU3Spin::A1A3 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  n2 *= 3;
  n2 += 2;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_2 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double FermionOnSphereWithSU3Spin::A2A2 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  ++n1;
  n2 *= 3;
  ++n2;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_2 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double FermionOnSphereWithSU3Spin::A2A3 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  ++n1;
  n2 *= 3;
  n2 += 2;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_3 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double FermionOnSphereWithSU3Spin::A3A3 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  n1 += 2;
  n2 *= 3;
  n2 += 2;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a^+_m1_1 a^+_m2_1 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int FermionOnSphereWithSU3Spin::Ad1Ad1 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  m2 *= 3;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_1 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int FermionOnSphereWithSU3Spin::Ad1Ad2 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  m2 *= 3;
  ++m2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_1 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int FermionOnSphereWithSU3Spin::Ad1Ad3 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  m2 *= 3;
  m2 += 2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_2 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int FermionOnSphereWithSU3Spin::Ad2Ad2 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  ++m1;
  m2 *= 3;
  ++m2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_2 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int FermionOnSphereWithSU3Spin::Ad2Ad3 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  ++m1;
  m2 *= 3;
  m2 += 2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_3 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int FermionOnSphereWithSU3Spin::Ad3Ad3 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  m1 += 2;
  m2 *= 3;
  m2 += 2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

