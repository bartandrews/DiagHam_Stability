////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                        Class author Cecile Repellin                        //
//                                                                            //
//                                                                            //
//                class of spin 1/2 chain with Sz contraint                   //
//                               and 2d translations                          //
//                                                                            //
//                        last modification : 19/06/2016                      //
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


#include "HilbertSpace/Spin1_2ChainNewAnd2DTranslation.h"
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
using std::dec;
using std::hex;
using std::endl;



// default constructor
//

Spin1_2ChainNewAnd2DTranslation::Spin1_2ChainNewAnd2DTranslation ()
{
  this->Flag.Initialize();
  this->ChainLength = 0;
  this->HilbertSpaceDimension = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->Sz = 0;
  this->StateDescription = 0;
  this->LookUpTable = 0;
  this->LookUpTableMask = 0;
  this->LookUpPosition = 0;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// xMomentum = momentum along the x direction
// maxXMomentum = number of sites in the x direction
// yMomentum = momentum along the y direction
// maxYMomentum = number of sites in the y direction
// memory = amount of memory granted for precalculations

Spin1_2ChainNewAnd2DTranslation::Spin1_2ChainNewAnd2DTranslation (int nbrSite, int sz, int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, unsigned long memory) 
{
  this->Flag.Initialize();
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
  this->NbrSite = nbrSite;
  this->ChainLength = this->NbrSite;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;

  this->Sz = sz;
  this->FixedQuantumNumberFlag = true;
  
  this->StateXShift = this->NbrSite / this->MaxXMomentum;
  this->ComplementaryStateXShift = this->NbrSite - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->NbrYMomentumBlocks = this->NbrSite / this->StateXShift;
  this->StateYShift = (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
  this->YMomentumFullMask = 0x0ul;
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
    }
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->NbrSite);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
 
  if (this->LargeHilbertSpaceDimension > 0)
    {
      this->GenerateLookUpTable(memory);
      
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
      UsedMemory = this->NbrSite * sizeof(int);
      UsedMemory += this->NbrSite * this->LookUpTableMemorySize * sizeof(int);
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
    
//      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//   {
//     unsigned long state = this->StateDescription[i];
//     int TmpMaxMomentum = this->NbrSite;
//     while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
//       --TmpMaxMomentum;
//     this->PrintState(cout, i) << " " << i << " " << this->FindStateIndex(state, TmpMaxMomentum) << endl;
//   }
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainNewAnd2DTranslation::Spin1_2ChainNewAnd2DTranslation (const Spin1_2ChainNewAnd2DTranslation& chain)
{
  this->Flag = chain.Flag;
  if (chain.NbrSite != 0)
    {
      this->NbrSite = chain.NbrSite;
      this->ChainLength = chain.NbrSite;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

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

      this->Sz = chain.Sz;

      this->StateDescription = chain.StateDescription;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
    }
  else
    {
      this->NbrSite = 0;
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->LookUpTableSize = 0;
    }
}

// destructor
//

Spin1_2ChainNewAnd2DTranslation::~Spin1_2ChainNewAnd2DTranslation () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainNewAnd2DTranslation& Spin1_2ChainNewAnd2DTranslation::operator = (const Spin1_2ChainNewAnd2DTranslation& chain)
{
  if ((this->Flag.Used() == true) && (this->NbrSite != 0) && (this->Flag.Shared() == false))
    {
      delete[] this->StateDescription;
    }
  this->Flag = chain.Flag;
  if (chain.NbrSite != 0)
    {
      this->NbrSite = chain.NbrSite;
      this->ChainLength = chain.NbrSite;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

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

      this->Sz = chain.Sz;

      this->StateDescription = chain.StateDescription;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
    }
  else
    {
      this->NbrSite = 0;
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainNewAnd2DTranslation::Clone()
{
  return new Spin1_2ChainNewAnd2DTranslation (*this);
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int Spin1_2ChainNewAnd2DTranslation::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{  
  unsigned long State = this->StateDescription[state];
  unsigned long Mask = (0x1ul << i) | (0x1ul << j);
  if (i != j)
    {
      if (((State & (0x1ul << j)) != 0x0ul) || ((State & (0x1ul << i)) == 0x0ul))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      State ^= 0x1ul << i;
      State |= 0x1ul << j;
      coefficient = 1.0;
      return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
    }
  else
    {
      nbrTranslationX = 0;
      nbrTranslationY = 0;
      coefficient = 0.5;
      return state;
    }
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int Spin1_2ChainNewAnd2DTranslation::SziSzjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{  
  unsigned long State = this->StateDescription[state];
  if (k != l)
    {
      if (((State & (0x1ul << k)) != 0x0ul) || ((State & (0x1ul << l)) == 0x0ul))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      State ^= 0x1ul << l;
      State |= 0x1ul << k;
      
      unsigned long Mask = ((0x1ul << i) | (0x1ul << j));
      unsigned long tmpState = State & Mask;
      if ((tmpState == 0x0ul) || (tmpState == Mask))
	coefficient = 0.25;
      else
	coefficient = -0.25;
      return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
    }
  else
    {
      nbrTranslationX = 0;
      nbrTranslationY = 0;
      coefficient = 0.5 * this->SziSzj(i, j, state);
      return state;
    }
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int Spin1_2ChainNewAnd2DTranslation::SmiSpjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{  
  unsigned long State = this->StateDescription[state];
  unsigned long Mask = (0x1ul << i) | (0x1ul << j);
  if ((i != j) && (k != l))
    {
      if (((State & (0x1ul << j)) != 0x0ul) || ((State & (0x1ul << i)) == 0x0ul))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      State ^= 0x1ul << i;
      State |= 0x1ul << j;
      if (((State & (0x1ul << l)) != 0x0ul) || ((State & (0x1ul << k)) == 0x0ul))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      State ^= 0x1ul << k;
      State |= 0x1ul << l;
      coefficient = 1.0;
      return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
    }
  if (k == l)
  {
    int TmpIndex;
    if (i != j)
    {
      if (((State & (0x1ul << j)) != 0x0ul) || ((State & (0x1ul << i)) == 0x0ul))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      State ^= 0x1ul << i;
      State |= 0x1ul << j;
      coefficient = 1.0;
      TmpIndex = this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
      coefficient *= 0.5;
      return TmpIndex;
    }
  else
    {
      coefficient = 0.25;
      return state;
    }
  }
  if (i == j)
  {
    int TmpIndex;
    if (k != l)
    {
      if (((State & (0x1ul << l)) != 0x0ul) || ((State & (0x1ul << k)) == 0x0ul))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      State ^= 0x1ul << k;
      State |= 0x1ul << l;
      coefficient = 1.0;
      TmpIndex = this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
      coefficient *= 0.5;
      return TmpIndex;
    }
  else
    {
      coefficient = 0.25;
      return state;
    }
  }
}


// find state index
//
// stateDescription = unsigned longeger describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int Spin1_2ChainNewAnd2DTranslation::FindStateIndex(unsigned long stateDescription, int maxMomentum)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    return this->HilbertSpaceDimension;
  if (this->LookUpTableShift[maxMomentum] < 0)
    return this->HilbertSpaceDimension;
  long PosMax = stateDescription >> this->LookUpTableShift[maxMomentum];
  long PosMin = this->LookUpTable[maxMomentum][PosMax];
  PosMax = this->LookUpTable[maxMomentum][PosMax + 1];
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

// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long Spin1_2ChainNewAnd2DTranslation::GenerateStates()
{
  
  cout << (this->LargeHilbertSpaceDimension) << endl;
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

// generate all states without the momentum constraint
//
// return value = Hilbert space dimension

long Spin1_2ChainNewAnd2DTranslation::RawGenerateStates(int nbrSpinUp, int currentPosition, long pos)
{  
  if (nbrSpinUp == 0)
    {
      this->StateDescription[pos] = 0x0ul;
      ++pos;
      return pos;
    }
  if (currentPosition < (nbrSpinUp - 1))
    return pos;
  int ReducedCurrentPosition = currentPosition - 1;
  long TmpPos = this->RawGenerateStates(nbrSpinUp - 1, ReducedCurrentPosition, pos);
  unsigned long Mask = 0x1ul << currentPosition;
  for (long i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  return this->RawGenerateStates(nbrSpinUp, ReducedCurrentPosition, TmpPos);
}
 
// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void Spin1_2ChainNewAnd2DTranslation::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->NbrSite);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrSite)
    this->MaximumLookUpShift = this->NbrSite;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrSite];
  this->LookUpTableShift = new int [this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentMaxMomentum = this->ChainLength;
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

void Spin1_2ChainNewAnd2DTranslation::ComputeRescalingFactors()
{
  this->RescalingFactors = new double* [this->NbrSite + 1];
  for (int i = 1; i <= this->NbrSite; ++i)
    {
      this->RescalingFactors[i] = new double [this->NbrSite + 1];
      for (int j = 1; j <= this->NbrSite; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}

// convert a state defined in the real space basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector Spin1_2ChainNewAnd2DTranslation::ConvertToKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  Spin1_2ChainNew* TmpSpace = (Spin1_2ChainNew*) space;
  ComplexVector TmpVector (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      int Pos = TmpSpace->FindStateIndex(TmpState);
      if (Pos < TmpSpace->HilbertSpaceDimension)
	{
	  TmpVector[i] =  state[Pos] * sqrt((double) this->NbrStateInOrbit[i]);
	}
    }
  return TmpVector;
}

// convert a state defined in the (Kx,Ky) basis into a state in the real space basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector Spin1_2ChainNewAnd2DTranslation::ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  Spin1_2ChainNew* TmpSpace = (Spin1_2ChainNew*) space;
  ComplexVector TmpVector (TmpSpace->HilbertSpaceDimension, true);
  Complex** FourrierCoefficients = new Complex* [this->MaxXMomentum];
  
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      FourrierCoefficients[i] = new Complex [this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{
	  FourrierCoefficients[i][j] = Phase (-2.0 * M_PI * ((double) (i * this->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * this->YMomentum) / ((double) this->MaxYMomentum)));
	}
    }
  
  for (int i = 0; i < TmpSpace->HilbertSpaceDimension; ++i)
    {
      int NbrTranslationX;
      int NbrTranslationY;
      unsigned long TmpState = TmpSpace->StateDescription[i];
      TmpState = this->FindCanonicalForm(TmpState, NbrTranslationX, NbrTranslationY);
      NbrTranslationX = (this->MaxXMomentum - NbrTranslationX) % this->MaxXMomentum;
      NbrTranslationY = (this->MaxYMomentum - NbrTranslationY) % this->MaxYMomentum;
      int TmpMaxMomentum = this->NbrSite;
      while (((TmpState >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
	--TmpMaxMomentum;
      
      int Pos = this->FindStateIndex(TmpState, TmpMaxMomentum);
      if (Pos < this->HilbertSpaceDimension)
	{
	  TmpVector[i] =  (state[Pos] * FourrierCoefficients[NbrTranslationX][NbrTranslationY] / sqrt((double) this->NbrStateInOrbit[Pos]));
	}
    }
  delete[] FourrierCoefficients;
  return TmpVector;
}
  

