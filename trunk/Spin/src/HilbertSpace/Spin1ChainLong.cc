////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of spin 1 chain with translations for more than 32 spins      //
//                                                                            //
//                        last modification : 18/03/2019                      //
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


#include "HilbertSpace/Spin1ChainLong.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;



// default constructor
//

Spin1ChainLong::Spin1ChainLong () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->MaximumLookUpShift = 0;
  this->LookUpTableMemorySize = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->ChainLength = 0;
  this->Sz = 0;
  this->FixedSpinProjectionFlag = false;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// memorySize = memory size in bytes allowed for look-up table

Spin1ChainLong::Spin1ChainLong (int chainLength, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;

  this->LargeHilbertSpaceDimension = 3l;
  for (int i = 1; i < this->ChainLength; i++)
    this->LargeHilbertSpaceDimension *= 3l;

  this->StateDescription = new ULONGLONG [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory +=  this->LargeHilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->ChainLength * sizeof(int);
      UsedMemory += this->ChainLength * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}

// constructor for complete Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Spin1ChainLong::Spin1ChainLong (int chainLength, int sz, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz;
  this->FixedSpinProjectionFlag = true;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(0, this->ChainLength);
  this->StateDescription = new ULONGLONG [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1, 0);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory +=  this->LargeHilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = 2 * this->ChainLength * sizeof(int);
      UsedMemory += 2 * this->ChainLength * (this->LookUpTableMemorySize + 1) * sizeof(int);
      cout << "memory requested for lookup table = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1ChainLong::Spin1ChainLong (const Spin1ChainLong& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->LookUpTable = chain.LookUpTable;
      this->MaximumLookUpShift = chain.MaximumLookUpShift;
      this->LookUpTableMemorySize = chain.LookUpTableMemorySize;
      this->LookUpTableShift = chain.LookUpTableShift;

      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
    }
  else
    {
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0l;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;

      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;
    }
}

// destructor
//

Spin1ChainLong::~Spin1ChainLong () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true) && (this->LargeHilbertSpaceDimension != 0l))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTableShift;
      int TmpMaxBitPosition = 2 * this->ChainLength;
      for (int i = 0; i < TmpMaxBitPosition; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1ChainLong& Spin1ChainLong::operator = (const Spin1ChainLong& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTableShift;
      int TmpMaxBitPosition = 2 * this->ChainLength;
      for (int i = 0; i < TmpMaxBitPosition; ++i)
	this->LookUpTable[i];
      delete[] this->LookUpTable;
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->StateDescription = chain.StateDescription;
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->Sz = chain.Sz;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;

      this->LookUpTable = chain.LookUpTable;
      this->MaximumLookUpShift = chain.MaximumLookUpShift;
      this->LookUpTableMemorySize = chain.LookUpTableMemorySize;
      this->LookUpTableShift = chain.LookUpTableShift;
    }
  else
    {
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0l;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
 
      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1ChainLong::Clone()
{
  return new Spin1ChainLong (*this);
}

// evaluate Hilbert space dimension
//
// sz = twice the Sz value
// nbrSites = number of sites
// return value = Hilbert space dimension

long Spin1ChainLong::EvaluateHilbertSpaceDimension(int sz, int nbrSites)
{
  if (nbrSites == 0)
    {
      if (sz == this->Sz)
	{
	  return 1l;	  
	}
      else
	{
	  return 0l;	  
	}
    }
  if (((sz + 2 * nbrSites) < this->Sz) || ((sz - 2 * nbrSites) > this->Sz))
    {
      return 0l;
    }
  long TmpDimension = this->EvaluateHilbertSpaceDimension(sz - 2, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(sz, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(sz + 2, nbrSites - 1);
  return TmpDimension;
}

// generate all states with no constraint on total Sz and no discrete symmetry constraint
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// return value = number of generated states

long Spin1ChainLong::RawGenerateStates(long statePosition, int sitePosition) 
{
  if (sitePosition < 0)
    {
      this->StateDescription[statePosition] = ((ULONGLONG) 0x0ul);
      return (statePosition + 1l);
    }

  ULONGLONG TmpMask = ((ULONGLONG) 0x3ul) << (sitePosition * 2);
  long TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = ((ULONGLONG) 0x2ul) << (sitePosition * 2);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  return this->RawGenerateStates(statePosition, sitePosition - 1);  
}

// generate all states corresponding to a given total Sz and no discrete symmetry constraint
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentSz = total Sz value of current state
// return value = number of generated states

long Spin1ChainLong::RawGenerateStates(long statePosition, int sitePosition, int currentSz) 
{
  if (sitePosition < 0)
    {
      if (currentSz == this->Sz)
	{
	  this->StateDescription[statePosition] = ((ULONGLONG) 0x0ul);
	  return (statePosition + 1l);
	}
      else
	{
	  return statePosition;
	}
    }
  if (((currentSz + 2 * sitePosition + 2) < this->Sz) || ((currentSz - 2 * sitePosition - 2) > this->Sz))
    {
      return statePosition;
    }
  ULONGLONG TmpMask = ((ULONGLONG) 0x3ul) << (sitePosition * 2);
  long TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1, currentSz + 2);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = ((ULONGLONG) 0x2ul) << (sitePosition * 2);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1, currentSz);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  return this->RawGenerateStates(statePosition, sitePosition - 1, currentSz - 2);  
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void Spin1ChainLong::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->ChainLength);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  int TmpMaxBitPosition = 2 * this->ChainLength;
  if (this->MaximumLookUpShift > TmpMaxBitPosition)
    this->MaximumLookUpShift = TmpMaxBitPosition;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [TmpMaxBitPosition];
  this->LookUpTableShift = new int [TmpMaxBitPosition];
  for (int i = 0; i < TmpMaxBitPosition; ++i)
    {
      this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
    }
  int CurrentMaxMomentum = TmpMaxBitPosition;
  while (((this->StateDescription[0] >> CurrentMaxMomentum) == ((ULONGLONG) 0x0ul)) && (CurrentMaxMomentum > 0))
    --CurrentMaxMomentum;
  int* TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
  if (CurrentMaxMomentum < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentMaxMomentum] = 0;
  else
    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
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
      int TmpMaxMomentum = CurrentMaxMomentum;
      while (((this->StateDescription[i] >> TmpMaxMomentum) == ((ULONGLONG) 0x0ul)) && (TmpMaxMomentum > 0))
	--TmpMaxMomentum;
      if (CurrentMaxMomentum != TmpMaxMomentum)
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  --CurrentMaxMomentum;
	  while (CurrentMaxMomentum > TmpMaxMomentum)
	    {
	      this->LookUpTableShift[CurrentMaxMomentum] = -1;
	      --CurrentMaxMomentum;
	    }
 	  CurrentMaxMomentum = TmpMaxMomentum;
	  TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
	  if (CurrentMaxMomentum < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentMaxMomentum] = 0;
	  else
	    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
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
}

// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long Spin1ChainLong::GenerateStates()
{
  return this->LargeHilbertSpaceDimension;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Spin1ChainLong::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedSpinProjectionFlag == true)
    {
      TmpList += new SzQuantumNumber (this->Sz);
    }
  else
    {
      int TmpSz = - 2 * this->ChainLength;
      for (int i = 0; i <= (2 * this->ChainLength); i++)
	{
	  TmpList += new SzQuantumNumber (TmpSz);
	  TmpSz += 2;
	}
    }
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* Spin1ChainLong::GetQuantumNumber (int index)
{ 
  return new SzQuantumNumber (this->TotalSz(index));
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin1ChainLong::TotalSz (int index)
{
  if (this->FixedSpinProjectionFlag == true)
    return this->Sz;
  ULONGLONG State = this->StateDescription[index];
  int TmpSz = 0;
  ULONGLONG TmpState;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpState = State & ((ULONGLONG) 0x3ul);
      switch (TmpState)
	{
	case ((ULONGLONG) 0x3ul):
	  TmpSz += 2;
	  break;
	case ((ULONGLONG) 0x0ul):
	  TmpSz -= 2;
	  break;
	}
      State >>= 2;
    }
  return TmpSz;
}

// return matrix representation of Sx
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin1ChainLong::Sxi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  i <<= 1;
  ULONGLONG tmpState;
  ULONGLONG State;
  double Factor = M_SQRT2 * 0.5;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      tmpState = this->StateDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= ((ULONGLONG) 0x3ul);
      switch (tmpState)
	{
	case ((ULONGLONG) 0x3ul):
	  M(this->FindStateIndex((State & ~(((ULONGLONG) 0x1ul) << i))), j) = Factor;
	  break;
	case ((ULONGLONG) 0x2ul):
	  M(this->FindStateIndex((State | (((ULONGLONG) 0x1ul) << i))), j) = Factor;
	  M(this->FindStateIndex((State & ~(((ULONGLONG) 0x2ul) << i))), j) = Factor;
	  break;
	case ((ULONGLONG) 0x0ul):
	  M(this->FindStateIndex((State | (((ULONGLONG) 0x2ul) << i))), j) = Factor;
	  break;
	}
    }
  return M;
}

// return matrix representation of i * Sy
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin1ChainLong::Syi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  i <<= 1;
  ULONGLONG tmpState;
  ULONGLONG State;
  double Factor = M_SQRT2 * 0.5;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      tmpState = this->StateDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= ((ULONGLONG) 0x3ul);
      switch (tmpState)
	{
	case ((ULONGLONG) 0x3ul):
	  M(this->FindStateIndex((State & ~(((ULONGLONG) 0x1ul) << i))), j) = -Factor;
	  break;
	case ((ULONGLONG) 0x2ul):
	  M(this->FindStateIndex((State | (((ULONGLONG) 0x1ul) << i))), j) = Factor;
	  M(this->FindStateIndex((State & ~(((ULONGLONG) 0x2ul) << i))), j) = -Factor;
	  break;
	case ((ULONGLONG) 0x0ul):
	  M(this->FindStateIndex((State | (((ULONGLONG) 0x2ul) << i))), j) = Factor;
	  break;
	}
    }
  return M;
}

// return matrix representation of Sz
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin1ChainLong::Szi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  i <<= 1;
  ULONGLONG tmpState;
  ULONGLONG State;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      tmpState = this->StateDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= ((ULONGLONG) 0x3ul);
      switch (tmpState)
	{
	case ((ULONGLONG) 0x3ul):
	  M(j, j) = 1.0;
	  break;
	case ((ULONGLONG) 0x0ul):
	  M(j, j) = -1.0;
	  break;
	}
    }
  return M;
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainLong::Spi (int i, int state, double& coefficient)
{
  ULONGLONG State = (this->StateDescription[state] >> (i << 1));
  ULONGLONG tmpState = State;
  tmpState &= ((ULONGLONG) 0x3ul);
  switch (tmpState)
    {
    case ((ULONGLONG) 0x2ul):
      {
	coefficient = M_SQRT2;
	return this->SymmetrizeResult(State | (((ULONGLONG) 0x1ul) << (i << 1)), coefficient);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	coefficient = M_SQRT2;
	return this->SymmetrizeResult(State | (((ULONGLONG) 0x2ul) << (i << 1)), coefficient);
      }
      break;
    case ((ULONGLONG) 0x3ul):
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
      break;
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainLong::Smi (int i, int state, double& coefficient)
{
  ULONGLONG State = (this->StateDescription[state] >> (i << 1));
  ULONGLONG tmpState = State;
  tmpState &= ((ULONGLONG) 0x3ul);
  switch (tmpState)
    {
    case ((ULONGLONG) 0x2ul):
      {
	coefficient = M_SQRT2;
	return this->SymmetrizeResult(State & ~(((ULONGLONG) 0x2ul) << (i << 1)), coefficient);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
      break;
    case ((ULONGLONG) 0x3ul):
      {
	coefficient = M_SQRT2;
	return this->SymmetrizeResult(State & ~(((ULONGLONG) 0x1ul) << (i << 1)), coefficient);
      }
      break;
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i operator on a given state
//
// i = position of Sz operator
// state = index of the state to be applied on Sz_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainLong::Szi (int i, int state, double& coefficient)
{
  ULONGLONG tmpState = (this->StateDescription[state] >> (i << 1));
  tmpState &= ((ULONGLONG) 0x3ul);
  switch (tmpState)
    {
    case ((ULONGLONG) 0x2ul):
      coefficient = 0.0;
      break;
    case ((ULONGLONG) 0x0ul):
      coefficient = -1.0;
      break;
    case ((ULONGLONG) 0x3ul):
      coefficient = 1.0;
      break;
    }
  return state;
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double Spin1ChainLong::SziSzj (int i, int j, int state)
{  
  ULONGLONG tmpState = this->StateDescription[state];
  ULONGLONG tmpState2 = (tmpState >> (j << 1)) & ((ULONGLONG) 0x3ul);
  tmpState >>= (i << 1);
  tmpState &= ((ULONGLONG) 0x3ul);
  if ((tmpState == ((ULONGLONG) 0x2ul)) || (tmpState2 == ((ULONGLONG) 0x2ul)))
    return 0.0;
  if (tmpState == tmpState2)
    return 1.0;
  else
    return -1.0;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainLong::SmiSpj (int i, int j, int state, double& coefficient)
{  
  ULONGLONG tmpState = this->StateDescription[state];
  ULONGLONG State = tmpState;
  ULONGLONG tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= ((ULONGLONG) 0x3ul);
  switch (tmpState)
    {
    case ((ULONGLONG) 0x3ul):
      coefficient = M_SQRT2;
      State &= ~(((ULONGLONG) 0x1ul) << i);
      break;
    case ((ULONGLONG) 0x2ul):
      coefficient = M_SQRT2;
      State&= ~(((ULONGLONG) 0x2ul) << i);
      break;
    case ((ULONGLONG) 0x0ul):
      return this->HilbertSpaceDimension;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= ((ULONGLONG) 0x3ul);
  switch (tmpState2)
    {
    case ((ULONGLONG) 0x3ul):
      return this->HilbertSpaceDimension;
    case ((ULONGLONG) 0x2ul):
      coefficient *= M_SQRT2;
      return this->SymmetrizeResult(State | (((ULONGLONG) 0x1ul) << j), coefficient);
    case ((ULONGLONG) 0x0ul):
      coefficient *= M_SQRT2;
      return this->SymmetrizeResult(State | (((ULONGLONG) 0x2ul) << j), coefficient);
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainLong::SpiSpj (int i, int j, int state, double& coefficient)
{
  ULONGLONG tmpState = this->StateDescription[state];
  ULONGLONG State = tmpState;
  ULONGLONG tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= ((ULONGLONG) 0x3ul);
  switch (tmpState)
    {
    case ((ULONGLONG) 0x3ul):
      return this->HilbertSpaceDimension;
      break;
    case ((ULONGLONG) 0x2ul):
      coefficient = M_SQRT2;
      State |= (((ULONGLONG) 0x1ul) << i);
      break;
    case ((ULONGLONG) 0x0ul):
      coefficient = M_SQRT2;
      State |= (((ULONGLONG) 0x2ul) << i);
      break;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= ((ULONGLONG) 0x3ul);
  switch (tmpState2)
    {
    case ((ULONGLONG) 0x3ul):
      return this->HilbertSpaceDimension;
    case ((ULONGLONG) 0x2ul):
      coefficient *= M_SQRT2;
      return this->SymmetrizeResult(State | (((ULONGLONG) 0x1ul) << j), coefficient);
    case ((ULONGLONG) 0x0ul):
      coefficient *= M_SQRT2;
      return this->SymmetrizeResult(State | (((ULONGLONG) 0x2ul) << j), coefficient);
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainLong::SmiSmj (int i, int j, int state, double& coefficient)
{
  ULONGLONG tmpState = this->StateDescription[state];
  ULONGLONG State = tmpState;
  ULONGLONG tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= ((ULONGLONG) 0x3ul);
  switch (tmpState)
    {
    case ((ULONGLONG) 0x3ul):
      coefficient = M_SQRT2;
      State &= ~(((ULONGLONG) 0x1ul) << i);
      break;
    case ((ULONGLONG) 0x2ul):
      coefficient = M_SQRT2;
      State&= ~(((ULONGLONG) 0x2ul) << i);
      break;
    case ((ULONGLONG) 0x0ul):
      return this->HilbertSpaceDimension;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= ((ULONGLONG) 0x3ul);
  switch (tmpState2)
    {
    case ((ULONGLONG) 0x3ul):
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State & ~(((ULONGLONG) 0x1ul) << j), coefficient);
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State & ~(((ULONGLONG) 0x2ul) << j), coefficient);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	coefficient = 0;
	return this->HilbertSpaceDimension;
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainLong::SpiSzj (int i, int j, int state, double& coefficient)
{
  ULONGLONG tmpState = this->StateDescription[state];
  ULONGLONG State = tmpState;
  ULONGLONG tmpState2 = (tmpState >> (j << 1)) & ((ULONGLONG) 0x3ul);
  if (tmpState2 == ((ULONGLONG) 0x2ul))
    {
      coefficient = 0;
      return this->HilbertSpaceDimension;
    }
  if (tmpState2 == ((ULONGLONG) 0x3ul))
    {
      coefficient = 1.0;
    }
  else
    {
      coefficient = -1.0;
    }
  i <<= 1;
  tmpState >>= i;
  tmpState &= ((ULONGLONG) 0x3ul);
  switch (tmpState)
    {
    case ((ULONGLONG) 0x3ul):
      {
	coefficient = 0;
	return this->HilbertSpaceDimension;
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State | (((ULONGLONG) 0x1ul) << i), coefficient);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State | (((ULONGLONG) 0x2ul) << i), coefficient);
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainLong::SmiSzj (int i, int j, int state, double& coefficient)
{
  ULONGLONG tmpState = this->StateDescription[state];
  ULONGLONG State = tmpState;
  ULONGLONG tmpState2 = (tmpState >> (j << 1)) & ((ULONGLONG) 0x3ul);
  if (tmpState2 == ((ULONGLONG) 0x2ul))
    {
      coefficient = 0;
      return this->HilbertSpaceDimension;
    }
  if (tmpState2 == ((ULONGLONG) 0x3ul))
    {
      coefficient = 1.0;
    }
  else
    {
      coefficient = -1.0;
    }
  i <<= 1;
  tmpState >>= i;
  tmpState &= ((ULONGLONG) 0x3ul);
  switch (tmpState)
    {
    case ((ULONGLONG) 0x3ul):
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State & ~(((ULONGLONG) 0x1ul) << i), coefficient);
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State & ~(((ULONGLONG) 0x1ul) << i), coefficient);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	coefficient *= 0;
	return this->HilbertSpaceDimension;
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i1 S+_j1 S-_i2 S+_j2 operator on a given state
//
// i1 = position of leftmost S- operator
// j1 = position of leftmost S+ operator
// i2 = position of rightmost S- operator
// j2 = position of rightmost S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainLong::SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient)
{
  ULONGLONG TmpState = this->StateDescription[state];
  j2 <<= 1;
  switch ((TmpState >> j2) & ((ULONGLONG) 0x3ul))
    {
    case ((ULONGLONG) 0x3ul):
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient = M_SQRT2;
	TmpState |= (((ULONGLONG) 0x1ul) << j2);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	coefficient = M_SQRT2;
 	TmpState |= (((ULONGLONG) 0x2ul) << j2);
      }
      break;
    }
  i2 <<= 1;
  switch ((TmpState >> i2) & ((ULONGLONG) 0x3ul))
    {
    case ((ULONGLONG) 0x3ul):
      {
	coefficient *= M_SQRT2;
 	TmpState &= ~(((ULONGLONG) 0x1ul) << i2);
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient *= M_SQRT2;
	TmpState &= ~(((ULONGLONG) 0x2ul) << i2);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	return this->HilbertSpaceDimension;
      }
      break;
    }
  j1 <<= 1;
  switch ((TmpState >> j1) & ((ULONGLONG) 0x3ul))
    {
    case ((ULONGLONG) 0x3ul):
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient *= M_SQRT2;
	TmpState |= (((ULONGLONG) 0x1ul) << j1);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	coefficient *= M_SQRT2;
 	TmpState |= (((ULONGLONG) 0x2ul) << j1);
      }
      break;
    }
  i1 <<= 1;
  switch ((TmpState >> i1) & ((ULONGLONG) 0x3ul))
    {
    case ((ULONGLONG) 0x3ul):
      {
	coefficient *= M_SQRT2;
	TmpState &= ~(((ULONGLONG) 0x1ul) << i1);
	return this->SymmetrizeResult(TmpState, coefficient);
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient *= M_SQRT2;
	TmpState &= ~(((ULONGLONG) 0x2ul) << i1);
	return this->SymmetrizeResult(TmpState, coefficient);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	return this->HilbertSpaceDimension;
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i1 Sz_j1 S-_i2 S+_j2 operator on a given state
//
// i1 = position of first Sz operator
// j1 = position of second Sz operator
// i2 = position of S- operator
// j2 = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainLong::SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient)
{
  ULONGLONG TmpState = this->StateDescription[state];
  j2 <<= 1;
  switch ((TmpState >> j2) & ((ULONGLONG) 0x3ul))
    {
    case ((ULONGLONG) 0x3ul):
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient = M_SQRT2;
	TmpState |= (((ULONGLONG) 0x1ul) << j2);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	coefficient = M_SQRT2;
 	TmpState |= (((ULONGLONG) 0x2ul) << j2);
      }
      break;
    }
  i2 <<= 1;
  switch ((TmpState >> i2) & ((ULONGLONG) 0x3ul))
    {
    case ((ULONGLONG) 0x3ul):
      {
	TmpState &= ~(((ULONGLONG) 0x1ul) << i2);
	ULONGLONG TmpState2 = (TmpState >> (i1 << 1)) & ((ULONGLONG) 0x3ul);
	ULONGLONG TmpState3 = (TmpState >> (j1 << 1)) & ((ULONGLONG) 0x3ul);
	if ((TmpState2 == ((ULONGLONG) 0x2ul)) || (TmpState3 == ((ULONGLONG) 0x2ul)))
	  {
	    coefficient = 0.0;
	    return this->HilbertSpaceDimension;
	  }
	if (TmpState2 == TmpState3)
	  coefficient *= M_SQRT2;
	else
	  coefficient *= -M_SQRT2;
	return this->SymmetrizeResult(TmpState, coefficient);
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	TmpState &= ~(((ULONGLONG) 0x2ul) << i2);
	ULONGLONG TmpState2 = (TmpState >> (i1 << 1)) & ((ULONGLONG) 0x3ul);
	ULONGLONG TmpState3 = (TmpState >> (j1 << 1)) & ((ULONGLONG) 0x3ul);
	if ((TmpState2 == ((ULONGLONG) 0x2ul)) || (TmpState3 == ((ULONGLONG) 0x2ul)))
	  {
	    coefficient = 0.0;
	    return this->HilbertSpaceDimension;
	  }
	if (TmpState2 == TmpState3)
	  coefficient *= M_SQRT2;
	else
	  coefficient *= -M_SQRT2;
	return this->SymmetrizeResult(TmpState, coefficient);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	return this->HilbertSpaceDimension;
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// translate a state assuming the system have periodic boundary conditions (increasing the site index)
//
// nbrTranslations = number of translations to apply
// state = index of the state to translate 
// return value = index of resulting state

int Spin1ChainLong::TranslateState (int nbrTranslations, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  TmpState = (((TmpState & (((ULONGLONG) 0x3ul) << ((this->ChainLength - nbrTranslations) - 1ul)) << 1) << nbrTranslations)
	      | (TmpState >> ((this->ChainLength - nbrTranslations) << 1)));
  return this->FindStateIndex(TmpState);
}


// find state index
//
// stateDescription = state description
// maxBitPosition = maximum bit set to one in stateDescription
// return value = corresponding index

int Spin1ChainLong::FindStateIndex(ULONGLONG stateDescription, int maxBitPosition)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    return this->HilbertSpaceDimension;
  if (this->LookUpTableShift[maxBitPosition] < 0)
    return this->HilbertSpaceDimension;
  long PosMax = stateDescription >> this->LookUpTableShift[maxBitPosition];
  long PosMin = this->LookUpTable[maxBitPosition][PosMax];
  PosMax = this->LookUpTable[maxBitPosition][PosMax + 1];
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
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin1ChainLong::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  ULONGLONG tmp;
  ULONGLONG StateDescription = this->StateDescription[state];  
  //Str << this->FindStateIndex(StateDescription) << " : ";
  for (int j = 0; j < this->ChainLength; j++)
    {
      tmp = StateDescription & ((ULONGLONG) 0x3ul);
      if (tmp == ((ULONGLONG) 0x0ul))
	Str << "-1 ";
      else
	if (tmp == ((ULONGLONG) 0x2ul))
	  Str << "0 ";
	else
	  Str << "1 ";
      StateDescription >>= 2;
    }
  return Str;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix Spin1ChainLong::EvaluatePartialDensityMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
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
  if (nbrSites == this->ChainLength)
    {
      if (szSector == this->Sz)
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
  Spin1ChainLong TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1ChainLong TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = nbrSites * 2;
  ULONGLONG Mask = (((ULONGLONG) 0x1ul) << Shift) - ((ULONGLONG) 0x1ul);
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  long TmpNbrNonZeroElements = 0l;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      int Pos = 0;
      ULONGLONG TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  ULONGLONG TmpState2 = TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask);
	  int TmpPos = this->FindStateIndex(TmpState2);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]];
	      for (int k = 0; k < Pos; ++k)
		{
		  if (TmpStatePosition2[k] >= Pos2)
		    {
		      TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
		    }
		}
	    }
	}
    }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  return TmpDensityMatrix;
}
	
// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix Spin1ChainLong::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
	{
	  RealMatrix TmpEntanglementMatrix(1, 1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
      
    }
  if (nbrSites == this->ChainLength)
    {
      if (szSector == this->Sz)
	{
	  RealMatrix TmpEntanglementMatrix(1, 1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}      
    }
  Spin1ChainLong TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1ChainLong TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);
  RealMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int Shift = 2 * nbrSites;
//  ULONGLONG Mask = (((ULONGLONG) 0x1ul) << Shift) - ((ULONGLONG) 0x1ul);
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
//      ULONGLONG TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift) & 0xfffffffful;
      ULONGLONG TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
//	  ULONGLONG TmpState2 = TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask);
	  ULONGLONG TmpState2 = TmpState | TmpDestinationHilbertSpace.StateDescription[j];
	  int TmpPos = this->FindStateIndex(TmpState2);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	       TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos]);
	    }
	}
    }

   return TmpEntanglementMatrix;
}
	
// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix Spin1ChainLong::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, 1);
          Complex Tmp(1.0, 0.0);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
      
    }
  if (nbrSites == this->ChainLength)
    {
      if (szSector == this->Sz)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, 1);
          Complex Tmp(1.0, 0.0);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}      
    }
  Spin1ChainLong TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1ChainLong TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);

  ComplexMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = 2 * nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      ULONGLONG TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  ULONGLONG TmpState2 = TmpState | TmpDestinationHilbertSpace.StateDescription[j];
	  int TmpPos = this->FindStateIndex(TmpState2);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	       TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos]);
	    }
	}
    }

   return TmpEntanglementMatrix;
}
