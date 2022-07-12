////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                        Class author: Cecile Repellin                       //
//                                                                            //
//                                                                            //
//                     class of SU(3) (3 color) spin chain                    //
//                                                                            //
//                        last modification : 13/11/2017                      //
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


#include "HilbertSpace/SU3SpinChain.h"
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

SU3SpinChain::SU3SpinChain () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->MaximumLookUpShift = 0;
  this->LookUpTableMemorySize = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->ChainLength = 0;
  this->Tz = 0;
  this->Y = 0;
  this->FixedSpinProjectionFlag = false;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// constructor for complete Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// tz = twice the value of total tz component
// y = value of the total y component
// memorySize = memory size in bytes allowed for look-up table

SU3SpinChain::SU3SpinChain (int chainLength, int tz, int y, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Tz = tz;
  this->Y = y;
  this->FixedSpinProjectionFlag = true;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(0, 0, this->ChainLength);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//   {
//     cout << (i) << " : " ;
//     this->PrintState(cout, i);
//     cout << endl;
//   }
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory +=  this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
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

SU3SpinChain::SU3SpinChain (const SU3SpinChain& chain) 
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
      this->Tz = chain.Tz;
      this->Y = chain.Y;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
    }
  else
    {
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0l;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Tz = 0;
      this->Y = 0;
      this->FixedSpinProjectionFlag = false;

      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;
    }
}

// destructor
//

SU3SpinChain::~SU3SpinChain () 
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

SU3SpinChain& SU3SpinChain::operator = (const SU3SpinChain& chain)
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

      this->Tz = chain.Tz;
      this->Y = chain.Y;
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
      this->Tz = 0;
      this->Y = 0;
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

AbstractHilbertSpace* SU3SpinChain::Clone()
{
  return new SU3SpinChain (*this);
}

// evaluate Hilbert space dimension
//
// sz = twice the Sz value
// nbrSites = number of sites
// return value = Hilbert space dimension

long SU3SpinChain::EvaluateHilbertSpaceDimension(int currentTotalTz, int currentTotalY, int nbrSites)
{
  if (((currentTotalTz + 2 * nbrSites) < this->Tz) || ((currentTotalTz - 2 * nbrSites) > this->Tz) || ((currentTotalY + 3 * nbrSites) < this->Y) || ((currentTotalY - 3 * nbrSites) > this->Y ))
    {
      return 0l;
    }
    
  if (nbrSites == 0)
    {
      if ((currentTotalTz == this->Tz) && (currentTotalY == this->Y))
	{
	  return 1l;  
	}
      else
	{
	  return 0l;
	}
    }
  long TmpDimension = this->EvaluateHilbertSpaceDimension(currentTotalTz, currentTotalY - 2, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(currentTotalTz - 1, currentTotalY + 1, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(currentTotalTz + 1, currentTotalY + 1, nbrSites - 1);
  
  return TmpDimension;
}
// generate all states corresponding to a given total Sz and no discrete symmetry constraint
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentSz = total Sz value of current state
// return value = number of generated states

long SU3SpinChain::RawGenerateStates(long statePosition, int sitePosition, int currentTotalTz, int currentTotalY) 
{  
  if (sitePosition < 0)
    {
//       cout << currentTotalTz << " " << currentTotalY << " " << sitePosition << endl;
      if ((currentTotalTz == this->Tz) && (currentTotalY == this->Y))
	{
	  this->StateDescription[statePosition] = 0x0ul;
	  return (statePosition + 1l);
	}
      else
	{
	  return statePosition;
	}
    }
  if (((currentTotalTz + 2 * (sitePosition + 1)) < this->Tz) || ((currentTotalTz - 2 * (sitePosition + 1)) > this->Tz) || ((currentTotalY + 3 * (sitePosition + 1)) < this->Y) || ((currentTotalY - 3 * (sitePosition + 1)) > this->Y ))
    {
      return statePosition;
    }
  
  unsigned long TmpMask = 0x3ul << (sitePosition * 2);
  long TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1, currentTotalTz + 1, currentTotalY + 1);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
    
  TmpMask = 0x2ul << (sitePosition * 2);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1, currentTotalTz, currentTotalY - 2);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  
  return this->RawGenerateStates(statePosition, sitePosition - 1, currentTotalTz - 1, currentTotalY + 1);  
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void SU3SpinChain::GenerateLookUpTable(unsigned long memory)
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
//   cout << (this->StateDescription[0]) << endl;
  while (((this->StateDescription[0] >> CurrentMaxMomentum) == 0x0ul) && (CurrentMaxMomentum > 0))
    --CurrentMaxMomentum;
  int* TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
  if (CurrentMaxMomentum < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentMaxMomentum] = 0;
  else
    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
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
      int TmpMaxMomentum = CurrentMaxMomentum;
      while (((this->StateDescription[i] >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
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

long SU3SpinChain::GenerateStates()
{
  int TmpHilbertSpaceDimension = this->RawGenerateStates(0l, this->ChainLength - 1, 0, 0);
  return TmpHilbertSpaceDimension;
}


// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int SU3SpinChain::TotalTz (int index)
{
//   if (this->FixedSpinProjectionFlag == true)
//     return this->Tz;
  unsigned long State = this->StateDescription[index];
  int TmpSz = 0;
  unsigned long TmpState;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpState = State & 0x3ul;
      switch (TmpState)
	{
	case 0x3ul:
	  TmpSz += 1;
	  break;
	case 0x0ul:
	  TmpSz -= 1;
	  break;
	}
      State >>= 2;
    }
  return TmpSz;
}


// return index of resulting state from application of two-site exchange operator on a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = index of resulting state

int SU3SpinChain::Pij (int i, int j, int state)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long TmpMask1 = ((tmpState >> (i << 1)) & 0x3ul) ;
  unsigned long TmpMask2 = ((tmpState >> (j << 1)) & 0x3ul);
  
  if (TmpMask1 == TmpMask2)
    return state;
  
  unsigned long TmpMask3 = 0x3ul << (i << 1);
  unsigned long TmpMask4 = 0x3ul << (j << 1);
  
  
  tmpState &= (~(TmpMask3));
  tmpState &= (~(TmpMask4));
  
  tmpState |= (TmpMask1 << (j << 1));
  tmpState |= (TmpMask2 << (i << 1));
  
  double TmpCoefficient;
//   cout << (this->StateDescription[state]) << " " << tmpState << " " ;
  return this->SymmetrizeResult(tmpState, TmpCoefficient);
}


// return index of resulting state from application of four-site exchange operator on a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = index of resulting state

int SU3SpinChain::Pijkl (int i, int j, int k, int l, int state)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long TmpMask1 = ((tmpState >> (i << 1)) & 0x3ul);
  unsigned long TmpMask2 = ((tmpState >> (j << 1)) & 0x3ul);
  unsigned long TmpMask3 = ((tmpState >> (k << 1)) & 0x3ul) ;
  unsigned long TmpMask4 = ((tmpState >> (l << 1)) & 0x3ul);
  
  if ((TmpMask1 == TmpMask2) && (TmpMask2 == TmpMask3) && (TmpMask3 == TmpMask4))
    return state;
  
  
  unsigned long TmpMask = (0x3ul << (i << 1)) | (0x3ul << (j << 1)) | (0x3ul << (k << 1)) | (0x3ul << (l << 1)); 
  
  tmpState &= (~TmpMask);
  
  tmpState |= (TmpMask1 << (j << 1));
  tmpState |= (TmpMask2 << (k << 1));
  tmpState |= (TmpMask3 << (l << 1));
  tmpState |= (TmpMask4 << (i << 1));
  
  double TmpCoefficient;
//   cout << (this->StateDescription[state]) << " " << tmpState << " " ;
  return this->SymmetrizeResult(tmpState, TmpCoefficient);
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double SU3SpinChain::Szi (int i, int state)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = (tmpState >> (i << 1)) & 0x3ul;
  if (tmpState2 == 0x2ul)
    return 0.0;
  if (tmpState2 == 0x1ul)
    return 1.0;
  else
    return -1.0;
}


// translate a state assuming the system have periodic boundary conditions (increasing the site index)
//
// nbrTranslations = number of translations to apply
// state = index of the state to translate 
// return value = index of resulting state

int SU3SpinChain::TranslateState (int nbrTranslations, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  TmpState = (((TmpState & (0x3ul << ((this->ChainLength - nbrTranslations) - 1ul)) << 1) << nbrTranslations)
	      | (TmpState >> ((this->ChainLength - nbrTranslations) << 1)));
  return this->FindStateIndex(TmpState);
}


// find state index
//
// stateDescription = state description
// maxBitPosition = maximum bit set to one in stateDescription
// return value = corresponding index

int SU3SpinChain::FindStateIndex(unsigned long stateDescription, int maxBitPosition)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    return this->HilbertSpaceDimension;
  if (this->LookUpTableShift[maxBitPosition] < 0)
    return this->HilbertSpaceDimension;
  long PosMax = stateDescription >> this->LookUpTableShift[maxBitPosition];
  long PosMin = this->LookUpTable[maxBitPosition][PosMax];
//   cout << PosMax << " " << PosMin << " " ;
  PosMax = this->LookUpTable[maxBitPosition][PosMax + 1];
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

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& SU3SpinChain::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long tmp;
  unsigned long StateDescription = this->StateDescription[state];  
  //Str << this->FindStateIndex(StateDescription) << " : ";
  for (int j = 0; j < this->ChainLength; j++)
    {
      tmp = StateDescription & 0x3ul;
      if (tmp == 0x0ul)
	Str << "-1 ";
      else
	if (tmp == 0x2ul)
	  Str << "0 ";
	else
	  Str << "1 ";
      StateDescription >>= 2;
    }
  return Str;
}
