////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                    class author: Cecile Repellin                           //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with a fixed Sz value                //
//               and a pseudospin 1/2 (not a conserved quantity)              //
//                           and 2D translations                              //
//                                                                            //
//                        last modification : 25/11/2016                      //
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


#include "HilbertSpace/Spin1_2ChainWithPseudospinAnd2DTranslation.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"

#include <math.h>
#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;



// default constructor
//

Spin1_2ChainWithPseudospinAnd2DTranslation::Spin1_2ChainWithPseudospinAnd2DTranslation ()
{
}

// constructor for complete Hilbert space with a given total spin projection Sz
//
// chainLength = number of spin 1/2
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Spin1_2ChainWithPseudospinAnd2DTranslation::Spin1_2ChainWithPseudospinAnd2DTranslation (int chainLength, int sz, int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
#ifdef __64_BITS__
  if (this->ChainLength  > 64)
#else
  if (this->ChainLength  > 32)
#endif
    {
      this->ChainLength = 1;
    }
    
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
    
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  
  this->FixedQuantumNumberFlag = true;
  this->Sz = sz;
  
  this->StateXShift = 2 * this->ChainLength / this->MaxXMomentum;
  this->ComplementaryStateXShift = 2 * this->ChainLength- this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->NbrYMomentumBlocks = (2 * this->ChainLength) / this->StateXShift;
  this->StateYShift = ((2 * this->ChainLength) / (this->MaxYMomentum * this->MaxXMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
  this->YMomentumFullMask = 0x0ul;
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 
  
  this->FillMirrorSymmetryTable();
  
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->ChainLength);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->HilbertSpaceDimension > 0)
  {
    this->GenerateLookUpTable(memorySize);
  
//   double TmpCoef;
//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//   {
//     this->PrintState(cout,i) << endl ;
//     cout << (this->JOffDiagonali(1, i, TmpCoef)) << " ";
// //     if ((this->SmiSpj(0, 1, i, TmpCoef)) < this->HilbertSpaceDimension)
// //       cout << TmpCoef;
//     cout << endl;
//   }
//   cout << endl;
    this->Flag.Initialize();
//     for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//       this->PrintState(cout, i) << endl;
    
  }
  
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainWithPseudospinAnd2DTranslation::Spin1_2ChainWithPseudospinAnd2DTranslation (const Spin1_2ChainWithPseudospinAnd2DTranslation& chain)
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
      this->Sz = chain.Sz;
      this->StateDescription = chain.StateDescription;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->LookUpTableMemorySize = chain.LookUpTableMemorySize;
      this->MaximumLookUpShift = chain.MaximumLookUpShift;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
      this->Flag = chain.Flag;
      
      
      this->MaxXMomentum = chain.MaxXMomentum;
      this->XMomentum = chain.XMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;
      this->MaxYMomentum = chain.MaxYMomentum;
      this->YMomentum = chain.YMomentum;
      this->NbrYMomentumBlocks = chain.NbrYMomentumBlocks;
      this->StateYShift = chain.StateYShift;
      this->YMomentumBlockSize = chain.YMomentumBlockSize;
      this->ComplementaryStateYShift = chain.ComplementaryStateYShift;
      this->YMomentumMask = chain.YMomentumMask;
      this->YMomentumBlockMask = chain.YMomentumBlockMask;  
      this->YMomentumFullMask = chain.YMomentumFullMask;
      this->ComplementaryYMomentumFullMask = chain.ComplementaryYMomentumFullMask; 
      
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
    }
  else
    {
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->LookUpTableMemorySize = 0;
      this->MaximumLookUpShift = 0;
      
      this->MaxXMomentum = 0;
      this->XMomentum = 0;
      this->StateXShift = 0;
      this->ComplementaryStateXShift = 0;
      this->XMomentumMask = 0;
      this->MaxYMomentum = 0;
      this->YMomentum = 0;
      this->NbrYMomentumBlocks = 0;
      this->StateYShift = 0;
      this->YMomentumBlockSize = 0;
      this->ComplementaryStateYShift = 0;
      this->YMomentumMask = 0;
      this->YMomentumBlockMask = 0;  
      this->YMomentumFullMask = 0;
      this->ComplementaryYMomentumFullMask = 0; 
      
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
    }
}

// destructor
//

Spin1_2ChainWithPseudospinAnd2DTranslation::~Spin1_2ChainWithPseudospinAnd2DTranslation () 
{
  delete[] this->LookUpTableShift;
  
  int CurrentMaximumUpPosition = 2 * this->ChainLength - 1;
  while (CurrentMaximumUpPosition >=0)
  {
    delete[] this->LookUpTable[CurrentMaximumUpPosition];
    --CurrentMaximumUpPosition;
  }
  delete[] this->LookUpTable;
  
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainWithPseudospinAnd2DTranslation& Spin1_2ChainWithPseudospinAnd2DTranslation::operator = (const Spin1_2ChainWithPseudospinAnd2DTranslation& chain)
{
  if ((this->Flag.Used() == true) && (this->ChainLength != 0) && (this->Flag.Shared() == false))
    {
      delete[] this->StateDescription;
    }
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
      this->Sz = chain.Sz;
      this->StateDescription = chain.StateDescription;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->LookUpTableMemorySize = chain.LookUpTableMemorySize;
      this->MaximumLookUpShift = chain.MaximumLookUpShift;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
      this->Flag = chain.Flag;
    }
  else
    {
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->LookUpTableMemorySize = 0;
      this->MaximumLookUpShift = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainWithPseudospinAnd2DTranslation::Clone()
{
  return new Spin1_2ChainWithPseudospin (*this);
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin1_2ChainWithPseudospinAnd2DTranslation::TotalSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->Sz;
  unsigned long State = this->StateDescription[index];
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpSz += ((State & 0x1ul) << 1);
      State >>= 2;
    }
  TmpSz -= this->ChainLength;
  return TmpSz;
}

// compute the parity (prod_i Sz_i) for a given state
//
// state = index of the state to be applied on Sz_i operator
// return value = 0 if prod_i Sz_i = 1, 1 if prod_i Sz_i = -1

unsigned long Spin1_2ChainWithPseudospinAnd2DTranslation::GetParity (int state)
{
  unsigned long TmpState = this->StateDescription[state];
#ifdef __64_BITS__
  TmpState ^= TmpState >> 32;
#endif
  TmpState ^= TmpState >> 16;
  TmpState ^= TmpState >> 8;
  TmpState ^= TmpState >> 4;
  TmpState ^= TmpState >> 2;
  TmpState ^= TmpState >> 1;
  return (TmpState & 0x1ul);
  
}

// translate a state assuming the system have periodic boundary conditions (increasing the site index)
//
// nbrTranslations = number of translations to apply
// state = index of the state to translate 
// return value = index of resulting state

int Spin1_2ChainWithPseudospinAnd2DTranslation::TranslateState (int nbrTranslations, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  TmpState = (((TmpState & ((0x1ul << (this->ChainLength - nbrTranslations)) - 1ul)) << nbrTranslations)
	      | (TmpState >> (this->ChainLength - nbrTranslations)));
  return this->FindStateIndex(TmpState);
}


// find state index
//
// stateDescription = state description
// return value = corresponding index

int Spin1_2ChainWithPseudospinAnd2DTranslation::FindStateIndex(unsigned long stateDescription)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    return this->HilbertSpaceDimension;
  
  
  int CurrentMaximumUpPosition = 2 * this->ChainLength - 1;
  while ((((stateDescription >> CurrentMaximumUpPosition) & 0x1ul) == 0x0ul) && (CurrentMaximumUpPosition > 0))
    --CurrentMaximumUpPosition;
  long PosMax = stateDescription >> this->LookUpTableShift[CurrentMaximumUpPosition];
  long PosMin = this->LookUpTable[CurrentMaximumUpPosition][PosMax];
  PosMax = this->LookUpTable[CurrentMaximumUpPosition][PosMax + 1];
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
    {
      if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
	return this->HilbertSpaceDimension;
      else
	return PosMin;
    }
}

// evaluate Hilbert space dimension
//
// sz = twice the Sz value
// nbrSites = number of sites
// return value = Hilbert space dimension

long Spin1_2ChainWithPseudospinAnd2DTranslation::EvaluateHilbertSpaceDimension(int sz, int nbrSites)
{
  BinomialCoefficients TmpCoefficients (nbrSites);
  long PseudospinContribution = 1l << ((long) this->ChainLength);
  return (TmpCoefficients(nbrSites, (nbrSites + sz) >> 1) * PseudospinContribution);
}



// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long Spin1_2ChainWithPseudospinAnd2DTranslation::GenerateStates()
{
  cout << "Intermediary Hilbert space dimension = " << (this->LargeHilbertSpaceDimension) << endl;
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates((this->Sz + this->ChainLength) >> 1, this->ChainLength - 1, 0l);
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
#ifdef  __64_BITS__
  unsigned long Discard = 0xfffffffffffffffful;
#else
  unsigned long Discard = 0xfffffffful;
#endif
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY) == this->StateDescription[i]))
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      ++TmpLargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescription[i] = Discard;
	    }
	}
      else
	{
	  this->StateDescription[i] = Discard;
	}
    }
  if (TmpLargeHilbertSpaceDimension > 0)
    {
      unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
      this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
      TmpLargeHilbertSpaceDimension = 0l;
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if (this->StateDescription[i] != Discard)
	    {
	      TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	      this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = this->FindOrbitSize(this->StateDescription[i]);
	      ++TmpLargeHilbertSpaceDimension;
	    }
	}
      delete[] this->StateDescription;
      this->StateDescription = TmpStateDescription;
    }
  else
    {
      delete[] this->StateDescription;
    }
  return TmpLargeHilbertSpaceDimension;
}



// generate all states
// 
// nbrSpinUp = number of spin up
// currentPosition = current position to consider in the chain
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long Spin1_2ChainWithPseudospinAnd2DTranslation::RawGenerateStates(int nbrSpinUp, int currentPosition, long pos)
{
  if ((nbrSpinUp == 1) && (currentPosition == 0))
  {
//     this->StateDescription[pos + 1] = 0x1ul;
//     this->StateDescription[pos] = 0x3ul;
    this->StateDescription[pos + 1] = 0x2ul;
    this->StateDescription[pos] = 0x3ul;
    pos += 2;
    return pos;
  }
    
  if (nbrSpinUp == 0)
    {
      if (currentPosition < 0)
	return pos;
      else
      {
	if (currentPosition == 0)
	{
// 	  this->StateDescription[pos + 1] = 0x0ul;
// 	  this->StateDescription[pos] = (0x1ul << 1l);
	  this->StateDescription[pos + 1] = 0x0ul;
	  this->StateDescription[pos] = 0x1ul;
	  pos += 2;
	  return pos;
	}
	else
	{
	  long TmpPos = this->RawGenerateStates(nbrSpinUp, currentPosition - 1, pos);
// 	  unsigned long Mask = 0x1ul << (2 * currentPosition + 1);
	  unsigned long Mask = 0x1ul << (2 * currentPosition);
	  for (long i = pos; i < TmpPos; i++)
	  {
	    this->StateDescription[i + (TmpPos - pos)] = this->StateDescription[i];
	    this->StateDescription[i] |= Mask;
	  }
	  return (2 * TmpPos - pos);
	}
      }
    }
  if (currentPosition < (nbrSpinUp - 1))
    return pos;
  int ReducedCurrentPosition = currentPosition - 1;
  long TmpPos = this->RawGenerateStates(nbrSpinUp - 1, ReducedCurrentPosition, pos);
//   unsigned long Mask = 0x1ul << (2 * currentPosition);
  unsigned long Mask = 0x1ul << (2 * currentPosition);
  for (long i = pos; i < TmpPos; i++)
  {
//     this->StateDescription[i + (TmpPos - pos)] = this->StateDescription[i];
//     this->StateDescription[i] |= Mask;
//     this->StateDescription[i] |= (Mask << 0x1ul);
//     this->StateDescription[i + (TmpPos - pos)] |= Mask;
    this->StateDescription[i + (TmpPos - pos)] = this->StateDescription[i];
    this->StateDescription[i] |= Mask;
    this->StateDescription[i] |= (Mask << 0x1ul);
    this->StateDescription[i + (TmpPos - pos)] |= (Mask << 0x1ul);
  }
  pos = 2 * TmpPos - pos;
  TmpPos =  this->RawGenerateStates(nbrSpinUp, ReducedCurrentPosition, pos);
  for (long i = pos; i < TmpPos; i++)
  {
    this->StateDescription[i + (TmpPos - pos)] = this->StateDescription[i];
//     this->StateDescription[i] |= (Mask << 0x1ul);
    this->StateDescription[i] |= (Mask);
  }
  return (2 * TmpPos - pos);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin1_2ChainWithPseudospinAnd2DTranslation::PrintState (ostream& Str, int state)
{  
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long Mask = 0x1ul;
  for (int k = 0; k < this->ChainLength; k++)    
    {
      if ((this->StateDescription[state] & Mask) == 0x0ul)
	Str << "0 ";
      else
	Str << "1 ";
      Mask <<= 1;
      if ((this->StateDescription[state] & Mask) == 0x0ul)
	Str << "- ";
      else
	Str << "+ ";
      Mask <<= 1;
    }
//   Str << " " << hex << this->StateDescription[state] << dec;
//   Str << " " << this->FindStateIndex(this->StateDescription[state]);
  return Str;
}


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void Spin1_2ChainWithPseudospinAnd2DTranslation::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * 2 * this->ChainLength);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > 2 * this->ChainLength)
    this->MaximumLookUpShift = 2 * this->ChainLength;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [2 * this->ChainLength];
  this->LookUpTableShift = new int [2 * this->ChainLength];
  for (int i = 0; i < 2 * this->ChainLength; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentMaxMomentum = 2 * this->ChainLength - 1;
  int TmpMaxMomentum = 2 * this->ChainLength - 1;
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
  this->ComputeRescalingFactors();
}


// compute the rescaling factors
//

void Spin1_2ChainWithPseudospinAnd2DTranslation::ComputeRescalingFactors()
{
  this->RescalingFactors = new double* [this->ChainLength + 1];
  for (int i = 1; i <= this->ChainLength; ++i)
    {
      this->RescalingFactors[i] = new double [this->ChainLength + 1];
      for (int j = 1; j <= this->ChainLength; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}


// return the Bosonic Occupation of a given state in the basis
//
// index = index of the state in the basis
// finalState = reference on the array where the monomial representation has to be stored

// void Spin1_2ChainWithPseudospinAnd2DTranslation::GetBosonicOccupation (unsigned int index, int * finalState)
// {
//   for (int i = 0; i < this->ChainLength; i++)
//     {
//       finalState[i] = (this->StateDescription[index] >> ((unsigned long) (2 * i)) )& 0x1ul;
//     }
// }


// return eigenvalue of Sz_i Sz_j associated to a given state (acts only on spin part of many-body state)
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double Spin1_2ChainWithPseudospinAnd2DTranslation::SziSzj (int i, int j, int state)
{  
  unsigned long Mask = ((0x1ul << (2*i + 1)) | (0x1ul << (2*j + 1)));
  unsigned long tmpState = this->StateDescription[state] & Mask;
  if ((tmpState == 0x0ul) || (tmpState == Mask))
    return 0.25;
  else
    return -0.25;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainWithPseudospinAnd2DTranslation::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  int tmpOrbitSize = this->NbrStateInOrbit[state];
  tmpState >>= (2*i + 1);
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= (2*j + 1); 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  coefficient = 1.0;
	  State = ((State | (0x1ul << (2*j + 1))) & ~(0x1ul << (2*i + 1)));
	  return this->SymmetrizeResult(State, tmpOrbitSize, coefficient, nbrTranslationX, nbrTranslationY);
// 	  return this->FindStateIndex((State | (0x1ul << (2*j + 1))) & ~(0x1ul << (2*i + 1)));
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  coefficient = 0.5;
  return state;
//   if (tmpState == 0x0ul)
//     {
//       coefficient = 1.0;
//       return state;
//     }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainWithPseudospinAnd2DTranslation::SziSzjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  int tmpOrbitSize = this->NbrStateInOrbit[state];
  tmpState >>= (2*k + 1);
  tmpState &= 0x1ul;
  if (k != l)
    {
      tmpState2 >>= (2*l + 1); 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  coefficient = 1.0;
	  State = ((State | (0x1ul << (2*l + 1))) & ~(0x1ul << (2*k + 1)));
	  unsigned long Mask = ((0x1ul << (2*i + 1)) | (0x1ul << (2*j + 1)));
	  unsigned long tmpState = State & Mask;
	  if ((tmpState == 0x0ul) || (tmpState == Mask))
	    coefficient = 0.25;
	  else
	    coefficient = -0.25;
	  return this->SymmetrizeResult(State, tmpOrbitSize, coefficient, nbrTranslationX, nbrTranslationY);
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  coefficient = 0.5 * this->SziSzj(i, j, state);
  return state;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainWithPseudospinAnd2DTranslation::SmiSpjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{  
//   cout << i << " " << j << " " << k << " " << l << endl;
  if (((i == k) && (i !=j)) || ((j == l) && (k != l)))
    return this->HilbertSpaceDimension;
  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  unsigned long tmpState3 = tmpState;
  unsigned long tmpState4 = tmpState;
  int tmpOrbitSize = this->NbrStateInOrbit[state];
  tmpState >>= (2*i + 1);
  tmpState &= 0x1ul;
  
  tmpState3 >>= (2*k + 1);
  tmpState3 &= 0x1ul;
  
  tmpState4 >>= (2*l + 1); 
  tmpState4 &= 0x1ul;
  tmpState4 <<= 1;
  tmpState4 |= tmpState3;
  if (i != j)
    {
      if (k != l)
      {
	tmpState2 >>= (2*j + 1); 
	tmpState2 &= 0x1ul;
	tmpState2 <<= 1;
	tmpState2 |= tmpState;
	if ((tmpState2 == 0x1ul) && (tmpState4 == 0x1ul))
	{
	  coefficient = 1.0;
	  State = ((State | (0x1ul << (2*j + 1))) & ~(0x1ul << (2*i + 1)));
	  State = ((State | (0x1ul << (2*l + 1))) & ~(0x1ul << (2*k + 1)));
	  return this->SymmetrizeResult(State, tmpOrbitSize, coefficient, nbrTranslationX, nbrTranslationY);
	}
	else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}  
      }
      
      if ((tmpState2 == 0x1ul) && (tmpState3 == 0x0ul))
      {
	coefficient = 1.0;
	State = ((State | (0x1ul << (2*j + 1))) & ~(0x1ul << (2*i + 1)));
	return this->SymmetrizeResult(State, tmpOrbitSize, coefficient, nbrTranslationX, nbrTranslationY);
      }      
      else
	return this->HilbertSpaceDimension;
    }
  if (tmpState == 0x0ul)
    {
      return this->SmiSpj(k, l, state, coefficient, nbrTranslationX, nbrTranslationY);
    }
  return this->HilbertSpaceDimension;
}

// operator acting on pseudospin on site i (off-diagonal part)
//
// i = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospinAnd2DTranslation::JOffDiagonali (int i, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long State = this->StateDescription[state];
  coefficient = 1.0;
  unsigned long Tmp = (State >> (2*i)) & 0x1ul;
  
  if (Tmp == 0x0ul)    
  {
    State = State | (0x1ul << (2*i));
    return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
  }
  else
  {
    State = State & ~(0x1ul << (2*i));
    this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
  }
}


// operator acting on pseudospin on site i (off-diagonal part)
//
// i = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospinAnd2DTranslation::JoffiJoffj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  if (i == j)
  {
    coefficient = 1.0;
    nbrTranslationX = 0;
    nbrTranslationY = 0;
    return state;
  }
  unsigned long State = this->StateDescription[state];
  coefficient = 1.0;
  unsigned long Tmp = (State >> (2*i)) & 0x1ul;
  unsigned long  Tmp1 = (State >> (2*j)) & 0x1ul;
  unsigned long tmpState;
  if (Tmp == 0x0ul)    
  {
    tmpState = State | (0x1ul << (2*i));
//     cout << "*" << tmpState << "*";
    if (Tmp1 == 0x0ul)
      State =  tmpState | (0x1ul << (2*j));
    else
      State = tmpState & ~(0x1ul << (2*j));
  }
  else
  {
    tmpState = State & ~(0x1ul << (2*i));
//     cout << "*" << tmpState << "*";
    if (Tmp1 == 0x0ul)
      State =  tmpState | (0x1ul << (2*j));
    else
      State = tmpState & ~(0x1ul << (2*j));
  }
//   cout << "<" << State << "> ";
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
}


// operator acting on pseudospin on site i (off-diagonal part)
//
// i = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coupling = array where the coupling coefficients are stored
// return value = numerical coefficient

double Spin1_2ChainWithPseudospinAnd2DTranslation::JDiagonali (int i, int state, double* coupling)
{
  int Tmp = (int) ((this->StateDescription[state] >> (2*i)) & 0x1ul);
  return coupling[Tmp];
}


// operator acting on pseudospin on site i (off-diagonal) and j(diagonal part)
//
// i = position of pseudospin operator
// j = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospinAnd2DTranslation::SmiSpjJoffiJj (int i, int j, int state, double* coupling, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  int tmpOrbitSize = this->NbrStateInOrbit[state];
  tmpState >>= (2*i + 1);
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= (2*j + 1); 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  coefficient = 1.0;
	  State = ((State | (0x1ul << (2*j + 1))) & ~(0x1ul << (2*i + 1)));
	  int Tmp1 = (int) ((State >> (2*j)) & 0x1ul);
	  coefficient *= coupling[Tmp1];
	  
	  unsigned long Tmp = (State >> (2*i)) & 0x1ul;  
	  if (Tmp == 0x0ul)    
	    State = State | (0x1ul << (2*i));
	  else
	    State = State & ~(0x1ul << (2*i));	  
	  return this->SymmetrizeResult(State, tmpOrbitSize, coefficient, nbrTranslationX, nbrTranslationY);
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0x0ul)
    {
      coefficient = -0.25;
      int Tmp = (int) ((State >> (2*j)) & 0x1ul);
      coefficient *= coupling[Tmp];
      if (Tmp == 0)    
	State = State | (0x1ul << (2*i));
      else
	State = State & ~(0x1ul << (2*i));	  
      return this->SymmetrizeResult(State, tmpOrbitSize, coefficient, nbrTranslationX, nbrTranslationY);
    }
  return this->HilbertSpaceDimension;
}



// operator acting on pseudospin on site i (off-diagonal) and j(diagonal part)
//
// i = position of pseudospin operator
// j = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospinAnd2DTranslation::SmiSpjJiJoffj (int i, int j, int state, double* coupling, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  int tmpOrbitSize = this->NbrStateInOrbit[state];
  tmpState >>= (2*i + 1);
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= (2*j + 1); 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  coefficient = 1.0;
	  State = ((State | (0x1ul << (2*j + 1))) & ~(0x1ul << (2*i + 1)));
	  int Tmp1 = (int) ((State >> (2*i)) & 0x1ul);
	  coefficient *= coupling[Tmp1];
	  
	  unsigned long Tmp = (State >> (2*j)) & 0x1ul;  
	  if (Tmp == 0x0ul)    
	    State = State | (0x1ul << (2*j));
	  else
	    State = State & ~(0x1ul << (2*j));	  
	  return this->SymmetrizeResult(State, tmpOrbitSize, coefficient, nbrTranslationX, nbrTranslationY);
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0x0ul)
    {
      coefficient = -0.25;
      int Tmp = (int) ((State >> (2*j)) & 0x1ul);
      coefficient *= coupling[Tmp];
      if (Tmp == 0)    
	State = State | (0x1ul << (2*i));
      else
	State = State & ~(0x1ul << (2*i));	  
      return this->SymmetrizeResult(State, tmpOrbitSize, coefficient, nbrTranslationX, nbrTranslationY);
    }
  return this->HilbertSpaceDimension;
}



// operator acting on pseudospin on site i (off-diagonal) and j(diagonal part)
//
// i = position of pseudospin operator
// j = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospinAnd2DTranslation::SmiSpjJoffiJoffj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  unsigned long Tmp = (State >> (2*i)) & 0x1ul;  
  unsigned long Tmp1 = (State >> (2*j)) & 0x1ul;  
  int tmpOrbitSize = this->NbrStateInOrbit[state];
  tmpState >>= (2*i + 1);
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= (2*j + 1); 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  coefficient = 1.0;
	  State = ((State | (0x1ul << (2*j + 1))) & ~(0x1ul << (2*i + 1)));
	  
	  if (Tmp == 0x0ul)
	    State = State | (0x1ul << (2*i));
	  else
	    State = State & ~(0x1ul << (2*i));
	  
	  if (Tmp1 == 0x0ul)
	    State = State | (0x1ul << (2*j));
	  else
	    State = State & ~(0x1ul << (2*j));
	  return this->SymmetrizeResult(State, tmpOrbitSize, coefficient, nbrTranslationX, nbrTranslationY);
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  
  return this->HilbertSpaceDimension;
}

// operator acting on pseudospin on site i (off-diagonal) and j(diagonal part)
//
// i = position of pseudospin operator
// j = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospinAnd2DTranslation::SmiSpjJiJj (int i, int j, int state, double* couplingI, double* couplingJ, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  int tmpOrbitSize = this->NbrStateInOrbit[state];
  int Tmp = (int) ((State >> (2*i)) & 0x1ul);
  int Tmp1 = (int) ((State >> (2*j)) & 0x1ul);
  
  tmpState >>= (2*i + 1);
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= (2*j + 1); 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  coefficient = couplingI[Tmp] * couplingJ[Tmp1];
	  State = ((State | (0x1ul << (2*j + 1))) & ~(0x1ul << (2*i + 1))); 
	  return this->SymmetrizeResult(State, tmpOrbitSize, coefficient, nbrTranslationX, nbrTranslationY);
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0x0ul)
    {
      coefficient = -0.25 * couplingI[Tmp] * couplingJ[Tmp1];
      return this->SymmetrizeResult(State, tmpOrbitSize, coefficient, nbrTranslationX, nbrTranslationY);
    }
  return this->HilbertSpaceDimension;
}

// give a state description value to this->TransientState
//
// state = index of the state to convert

void Spin1_2ChainWithPseudospinAnd2DTranslation::InitializeTransientState (int state)
{
  this->TransientState = this->StateDescription[state];
  
//   cout << "0 " << (this->TransientState) ;
}


// apply SziSzj to transient state and return corresponding coefficient
//
// i = position of S- operator
// j = position of S+ operator
// return value = corresponding eigenvalue

double Spin1_2ChainWithPseudospinAnd2DTranslation::SziSzj(int i, int j)
{
  unsigned long Mask = ((0x1ul << (2*i + 1)) | (0x1ul << (2*j + 1)));
  unsigned long tmpState = this->TransientState & Mask;
//   cout << "szsz " << (this->TransientState);
  if ((tmpState == 0x0ul) || (tmpState == Mask))
    return 0.25;
  else
    return -0.25;
}
  
// apply S+S- to transient state
//
// i = position of S- operator
// j = position of S+ operator
// return value = corresponding eigenvalue

double Spin1_2ChainWithPseudospinAnd2DTranslation::SmiSpj(int i, int j)
{
  unsigned long tmpState = this->TransientState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= (2*i + 1);
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= (2*j + 1); 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  this->TransientState = ((this->TransientState | (0x1ul << (2*j + 1))) & ~(0x1ul << (2*i + 1)));
//   cout << "smsp " << (this->TransientState) << endl;
	  
	  return 1.0;
	}
      else
	{
	  return 0.0;
	}
    }
  return 0.5;
//   if (tmpState == 0x0ul)
//     {
//       coefficient = 1.0;
//       return state;
//     }
}
  
  
// apply off-diagonal part of pseudospin operator to the transient state
//
// i = position of the Joff operator
// return value = corresponding eigenvalue

double Spin1_2ChainWithPseudospinAnd2DTranslation::JOffDiagonali(int i)
{
  unsigned long Tmp = (this->TransientState >> (2*i)) & 0x1ul;
  
  if (Tmp == 0x0ul)    
    this->TransientState = this->TransientState | (0x1ul << (2*i));
  
  else
    this->TransientState = this->TransientState & ~(0x1ul << (2*i));

  
//   cout << "joff " << (this->TransientState);
  return 1.0;
}
  
  
// apply diagonal part of pseudospin operator to the transient state
//
// i = position of the JDiag operator
// coupling = array of coupling coefficients
// return value = corresponding eigenvalue

double Spin1_2ChainWithPseudospinAnd2DTranslation::JDiagonali(int i, double* coupling)
{
  int Tmp = (int) ((this->TransientState >> (2*i)) & 0x1ul);
  return coupling[Tmp];
}

// convert a state defined on a lattice with a number of sites equals to a multiple of three
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

RealVector Spin1_2ChainWithPseudospinAnd2DTranslation::ProjectToEffectiveSubspaceThreeToOne(ComplexVector& state, AbstractSpinChain* space)
{
  Spin1_2ChainNew* TmpSpace = (Spin1_2ChainNew*) space;
  RealVector TmpVector (this->HilbertSpaceDimension, true);
  RealVector TmpVector1 (this->HilbertSpaceDimension, true);
  
  cout << "Dimension full space = " << (TmpSpace->HilbertSpaceDimension) << endl;
  cout << "Dimension projected space = " << (this->HilbertSpaceDimension) << endl;
  unsigned long TmpState;
  unsigned long TmpStateEffective;
  
  int tmpStateI;
  int tmpStateJ;
  int tmpStateK;
  
  int tmpStateLocalSz;
  int tmpStateLocalPseudospin;
  int TmpIndexEffective;
  int localSz;
  int localSpinStructure;
  double Coefficient;
  
  int* TmpFlag = new int[this->HilbertSpaceDimension];
  double** Isometry = new double* [2];
  for (int i = 0; i < 2; ++i)
    Isometry[i] = new double [3];
  Isometry[1][0] = 1.0 / sqrt(6.0);
  Isometry[1][1] = 1.0 / sqrt(6.0);
  Isometry[1][2] = -2.0 / sqrt(6.0);
  Isometry[0][0] = 1.0 / sqrt(2.0);
  Isometry[0][1] = -1.0 / sqrt(2.0);
  Isometry[0][2] = 0.0;
  
  for (int j = 0; j < TmpSpace->HilbertSpaceDimension; ++j)
  {
    Coefficient = state[j].Re;
//     cout << j << " " << Coefficient << " " << state[j].Im << endl;
    TmpState = TmpSpace->StateDescription[j];
    int k = TmpSpace->ChainLength - 1;
    TmpVector1.ClearVector();
    for (int i = 0; i < this->HilbertSpaceDimension; ++i)
      TmpFlag[i] = 0;
    while ((fabs(Coefficient) > 1.0e-15) && (k > 0))
    {
      tmpStateK = (int) ((TmpState >> k) & 0x1ul);
      tmpStateJ = (int) ((TmpState >> (k - 1)) & 0x1ul);
      tmpStateI = (int) ((TmpState >> (k - 2)) & 0x1ul);
      
      localSz = tmpStateI + tmpStateJ + tmpStateK - 1;
      
      if ((localSz == 0) || (localSz == 1))
      {
	TmpIndexEffective = (k + 1)/3 - 1;
	
	for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	  {
	    TmpStateEffective = this->StateDescription[i];
	    tmpStateLocalSz = (int) ((TmpStateEffective >> (2 * TmpIndexEffective + 1)) & 0x1ul);
	    if (tmpStateLocalSz == localSz)
	    {
	      tmpStateLocalPseudospin = (int) ((TmpStateEffective >> (2 * TmpIndexEffective)) & 0x1ul);
	      if (tmpStateI == tmpStateJ)
		localSpinStructure = 0;
	      else
	      {
		if (tmpStateI == tmpStateK)
		  localSpinStructure = 1;
		else
		  localSpinStructure = 2;
	      }
// 	      cout << i << " " << (Coefficient * Isometry[tmpStateLocalPseudospin][localSpinStructure]) << endl;
	      if (TmpFlag[i] == 0)
	      {
		TmpVector1[i] = Isometry[tmpStateLocalPseudospin][localSpinStructure];
		TmpFlag[i] = 1;
	      }
	      else
		TmpVector1[i] *= Isometry[tmpStateLocalPseudospin][localSpinStructure];
	    }
	    else
	    {
	      TmpVector1[i] = 0.0;
	      TmpFlag[i] = 1;
	    }
	  }
      }
      else
      {
	k = 0;
	Coefficient = 0.0;
      }
      k -= 3;
    }
   TmpVector.AddLinearCombination(Coefficient, TmpVector1);
  }

//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//     cout << TmpVector[i] << endl;
  for (int i = 0; i < 2; ++i)
    delete[] Isometry[i];
  delete[] Isometry;
  delete[] TmpFlag;
  return TmpVector;
}

 
// initialize the mirror symmetry table
//
void Spin1_2ChainWithPseudospinAnd2DTranslation::FillMirrorSymmetryTable()
{
  this->MirrorTransformationTable = new int[this->ChainLength];
  if (this->MaxXMomentum == this->MaxYMomentum)
  {
    for (int i = 0; i < this->MaxXMomentum; ++i)
      for (int j = 0; j < this->MaxYMomentum; ++j)
	this->MirrorTransformationTable[this->GetLinearizedIndex(i, j)] = this->GetLinearizedIndex(j, i);
  }
    
  if ((this->ChainLength == 12) && (this->MaxXMomentum == 6) && (this->MaxYMomentum == 2))
  {
    this->MirrorTransformationTable[this->GetLinearizedIndex(0, 0)] = this->GetLinearizedIndex(0, 0);
    this->MirrorTransformationTable[this->GetLinearizedIndex(1, 0)] = this->GetLinearizedIndex(2, 1);
    this->MirrorTransformationTable[this->GetLinearizedIndex(2, 0)] = this->GetLinearizedIndex(4, 0);
    this->MirrorTransformationTable[this->GetLinearizedIndex(3, 0)] = this->GetLinearizedIndex(0, 1);
    this->MirrorTransformationTable[this->GetLinearizedIndex(4, 0)] = this->GetLinearizedIndex(2, 0);
    this->MirrorTransformationTable[this->GetLinearizedIndex(5, 0)] = this->GetLinearizedIndex(4, 1);
    this->MirrorTransformationTable[this->GetLinearizedIndex(0, 1)] = this->GetLinearizedIndex(3, 0);
    this->MirrorTransformationTable[this->GetLinearizedIndex(1, 1)] = this->GetLinearizedIndex(5, 1);
    this->MirrorTransformationTable[this->GetLinearizedIndex(2, 1)] = this->GetLinearizedIndex(1, 0);
    this->MirrorTransformationTable[this->GetLinearizedIndex(3, 1)] = this->GetLinearizedIndex(3, 1);
    this->MirrorTransformationTable[this->GetLinearizedIndex(4, 1)] = this->GetLinearizedIndex(5, 0);
    this->MirrorTransformationTable[this->GetLinearizedIndex(5, 1)] = this->GetLinearizedIndex(1, 1);
  }
}