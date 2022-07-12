////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of Potts-3 chain with the translation symmetry            //
//                                                                            //
//                        last modification : 01/01/2014                      //
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


#include "HilbertSpace/Potts3ChainWithTranslations.h"
#include "HilbertSpace/Potts3Chain.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include <iostream>


using std::cout;
using std::endl;


#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif


// default constructor
//

Potts3ChainWithTranslations::Potts3ChainWithTranslations () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableSize = 0;
  this->LookUpTableMask = 0;
  this->LookUpPosition = 0;
  this->HilbertSpaceDimension = 0;
  this->ChainDescription = 0;
  this->ChainLength = 0;
  this->Momentum = 0;
  this->MaxXMomentum = 0;
  this->ComplementaryStateShift = 0;
  this->Sz = 0;
  this->FixedQuantumNumberFlag = false;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = momentum sector
// memorySize = memory size in bytes allowed for look-up table

Potts3ChainWithTranslations::Potts3ChainWithTranslations (int chainLength, int momentum, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->ComplementaryStateShift = (this->ChainLength - 1) << 1;
  this->FixedQuantumNumberFlag = false;
  this->Momentum = momentum;
  this->MaxXMomentum = this->ChainLength;
  this->LargeHilbertSpaceDimension = 3l;
  for (int i = 1; i < chainLength; i++)
    this->LargeHilbertSpaceDimension *= 3l;
  

  this->LookUpPosition = 0;
  this->LookUpTableSize = 4;
  memorySize >>= 2;
  this->LookUpTableMask = 0xfffffffc;
  while ((this->LookUpPosition <= this->ChainLength) && (memorySize >=  4))
    {
      this->LookUpTableMask <<= 2;
      this->LookUpTableSize <<= 2;
      memorySize >>= 2;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];

  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = this->RawGenerateStates (this->ChainLength - 1, 0);
  this->GenerateStates();
}

// constructor for complete Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// sz = twice the value of total Sz component
// momemtum = momentum sector
// memorySize = memory size in bytes allowed for look-up table

Potts3ChainWithTranslations::Potts3ChainWithTranslations (int chainLength, int sz, int momentum, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->ComplementaryStateShift = (this->ChainLength - 1) << 1;
  this->Sz = sz % 3;
  this->FixedQuantumNumberFlag = true;
  this->Momentum = momentum;
  this->MaxXMomentum = this->ChainLength;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->ChainLength - 1, this->Sz);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->LookUpPosition = 0;
  this->LookUpTableSize = 4;
  memorySize >>= 1;
  this->LookUpTableMask = 0xfffffffc;
  while ((this->LookUpPosition < this->ChainLength) && (memorySize >=  4))
    {
      this->LookUpTableMask <<= 2;
      this->LookUpTableSize <<= 2;
      memorySize >>= 2;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];

  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = this->RawGenerateStates (this->ChainLength - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->GenerateStates ();
}


// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Potts3ChainWithTranslations::Potts3ChainWithTranslations (const Potts3ChainWithTranslations& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->ChainDescription = chain.ChainDescription;
      this->Sz = chain.Sz;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->Momentum = chain.Momentum;
      this->MaxXMomentum = chain.MaxXMomentum;
      this->RescalingFactors = chain.RescalingFactors;
   }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->NbrStateInOrbit = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
      this->LargeHilbertSpaceDimension = 0l;
      this->ComplementaryStateShift = 0;
      this->Momentum = 0;
    }
}

// destructor
//

Potts3ChainWithTranslations::~Potts3ChainWithTranslations () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true) && (this->LargeHilbertSpaceDimension != 0l))
    {
      delete[] this->ChainDescription;
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Potts3ChainWithTranslations& Potts3ChainWithTranslations::operator = (const Potts3ChainWithTranslations& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->ChainDescription;
      delete[] this->LookUpTable;
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainDescription = chain.ChainDescription;
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->LookUpPosition = chain.LookUpPosition;
      this->Sz = chain.Sz;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->Momentum = chain.Momentum;
      this->MaxXMomentum = chain.MaxXMomentum;
      this->RescalingFactors = chain.RescalingFactors;
   }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->NbrStateInOrbit = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
      this->LargeHilbertSpaceDimension = 0l;
      this->ComplementaryStateShift = 0;
      this->Momentum = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Potts3ChainWithTranslations::Clone()
{
  return new Potts3ChainWithTranslations (*this);
}

// evaluate Hilbert space dimension
//
// currentSite = current site to occupy
// currentSzValue = state current Sz value 
// return value = Hilbert space dimension

long Potts3ChainWithTranslations::EvaluateHilbertSpaceDimension(int currentSite, int currentSzValue)
{
  if (currentSite < 0)
    {
      if ((currentSzValue %3) == this->Sz)
	return 1l;
      else
	return 0l;
    }
  long TmpDimension = this->EvaluateHilbertSpaceDimension(currentSite - 1, currentSzValue + 2);
  TmpDimension += this->EvaluateHilbertSpaceDimension(currentSite - 1, currentSzValue + 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(currentSite - 1, currentSzValue);
  return TmpDimension;
}


// generate all states with no constraint on total Sz
//
// currentSite = current site to occupy
// currentPosition = current position of the state that has to be considered
// return value = number of generated states

long Potts3ChainWithTranslations::RawGenerateStates(int currentSite, long currentPosition) 
{
  if (currentSite < 0)
    {
      this->ChainDescription[currentPosition] = 0x0l;
      return (currentPosition + 1l);
    }
  long TmpPosition = this->RawGenerateStates(currentSite - 1, currentPosition);
  unsigned long TmpMask = 0x2ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    this->ChainDescription[currentPosition] |= TmpMask;
  TmpPosition = this->RawGenerateStates(currentSite - 1, currentPosition);
  TmpMask = 0x1ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    this->ChainDescription[currentPosition] |= TmpMask;
  return this->RawGenerateStates(currentSite - 1, currentPosition);
}

// generate all states corresponding to a given total Sz
//
// currentSite = current site to occupy
// currentSzValue = state current Sz value 
// currentPosition = current position of the state that has to be considered
// return value = number of generated states

long Potts3ChainWithTranslations::RawGenerateStates(int currentSite, int currentSzValue, long currentPosition) 
{
  if (currentSite < 0)
    {
      if ((currentSzValue %3) == this->Sz)
	{
	  this->ChainDescription[currentPosition] = 0x0l;
	  return (currentPosition + 1l);
	}
      else
	return currentPosition;
    }
  long TmpPosition = this->RawGenerateStates(currentSite - 1, currentSzValue + 2, currentPosition);
  unsigned long TmpMask = 0x2ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    {
      this->ChainDescription[currentPosition] |= TmpMask;
    }
  TmpPosition = this->RawGenerateStates(currentSite - 1, currentSzValue + 1, currentPosition);
  TmpMask = 0x1ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    {
      this->ChainDescription[currentPosition] |= TmpMask;
    }
  return this->RawGenerateStates(currentSite - 1, currentSzValue, currentPosition);
}

// generate all states corresponding to a given momnetum
//

void Potts3ChainWithTranslations::GenerateStates()
{
  unsigned long TmpLargeHilbertSpaceDimension = 0l;
  int TmpOrbitSize;
  int TmpNbrTranslation;
  unsigned long TmpCanonicalForm;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      TmpCanonicalForm = this->FindCanonicalForm(this->ChainDescription[i], TmpNbrTranslation, TmpOrbitSize);
      if ((TmpCanonicalForm == this->ChainDescription[i]) && 
	  (((TmpOrbitSize * this->Momentum) % this->ChainLength) == 0))
	++TmpLargeHilbertSpaceDimension;
    }
  unsigned long* TmpChainDescription = new unsigned long [TmpLargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      TmpCanonicalForm = this->FindCanonicalForm(this->ChainDescription[i], TmpNbrTranslation, TmpOrbitSize);
      if ((TmpCanonicalForm == this->ChainDescription[i]) && 
	  (((TmpOrbitSize * this->Momentum) % this->ChainLength) == 0))
	{
	  TmpChainDescription[TmpLargeHilbertSpaceDimension] = TmpCanonicalForm;
	  this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = TmpOrbitSize;
	  ++TmpLargeHilbertSpaceDimension;
	}
    }
  delete[] this->ChainDescription;
  this->ChainDescription = TmpChainDescription;
  this->LargeHilbertSpaceDimension = TmpLargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
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

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Potts3ChainWithTranslations::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedQuantumNumberFlag == true)
    {
      TmpList += new SzQuantumNumber (this->Sz);
    }
  else
    {
      for (int i = 0; i < 3; i++)
	{
	  TmpList += new SzQuantumNumber (i);
	}
    }
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* Potts3ChainWithTranslations::GetQuantumNumber (int index)
{ 
  return new SzQuantumNumber (this->TotalSz(index));
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Potts3ChainWithTranslations::TotalSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->Sz;
  unsigned long State = this->ChainDescription[index];
  unsigned long TmpSz = 0l;
  for (int i = 0; i < this->ChainLength; ++i)
    {
      TmpSz += (State & 0x3ul);
      State >>= 2;
    }
  TmpSz %= 3;
  return ((double) TmpSz);
}

// return value of the value of the sum of the square of spin projection on (Oz) 
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

double Potts3ChainWithTranslations::TotalSzSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return (this->Sz * this->Sz);
  unsigned long State = this->ChainDescription[index];
  unsigned long TmpSz = 0l;
  for (int i = 0; i < this->ChainLength; ++i)
    {
      TmpSz += (State & 0x3ul);
      State >>= 2;
    }
  TmpSz %= 3;
  return ((double) (TmpSz * TmpSz));
}

// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Potts3ChainWithTranslations::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->ChainDescription[state];
  unsigned long TmpState2 = TmpState;
  j <<= 1;
  TmpState2 >>= j;
  TmpState2 &= 0x3ul;
  TmpState &= ~(0x3ul << j);
  switch (TmpState2)
    {
    case 0x1ul:
      TmpState |= 0x2ul << j;
      break;
    case 0x0ul:
      TmpState |= 0x1ul << j;
      break;
    }	  
  i <<= 1;
  TmpState2 = TmpState;
  TmpState2 >>= i;
  TmpState2 &= 0x3ul;
  TmpState &= ~(0x3ul << i);
  switch (TmpState2)
    {
    case 0x1ul:
      TmpState |= 0x2ul << i;      
      break;
    case 0x0ul:
      TmpState |= 0x1ul << i;      
      break;
    }
  coefficient = 1.0;
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
  //  return this->FindStateIndexAndTransaltion(state, TmpState, nbrTranslation, i, coefficient);
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Potts3ChainWithTranslations::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->ChainDescription[state];
  unsigned long TmpState2 = TmpState;
  j <<= 1;
  TmpState2 >>= j;
  TmpState2 &= 0x3ul;
  TmpState &= ~(0x3ul << j);
  switch (TmpState2)
    {
    case 0x2ul:
      TmpState |= 0x1ul << j;
      break;
    case 0x0ul:
      TmpState |= 0x2ul << j;
      break;
    }	  
  i <<= 1;
  TmpState2 = TmpState;
  TmpState2 >>= i;
  TmpState2 &= 0x3ul;
  TmpState &= ~(0x3ul << i);
  switch (TmpState2)
    {
    case 0x2ul:
      TmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      TmpState |= 0x2ul << i;      
      break;
    }
  coefficient = 1.0;
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
  //  return this->FindStateIndexAndTransaltion(state, TmpState, nbrTranslation, i, coefficient);
}

// return index of resulting state from application of S+_i S+_i operator on a given state
//
// i = position of first S+ operator
// state = index of the state to be applied on S+_i S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Potts3ChainWithTranslations::SpiSpi (int i, int state, double& coefficient, int& nbrTranslation)
{
  return this->SpiSpj (i, i, state, coefficient, nbrTranslation);
}

// return index of resulting state from application of S-_i S-_i operator on a given state
//
// i = position of the S- operator
// state = index of the state to be applied on S-_i S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Potts3ChainWithTranslations::SmiSmi (int i, int state, double& coefficient, int& nbrTranslation)
{
  return this->SmiSmj (i, i, state, coefficient, nbrTranslation);
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Potts3ChainWithTranslations::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->ChainDescription[state];
  unsigned long TmpState2 = TmpState;
  j <<= 1;
  TmpState2 >>= j;
  TmpState2 &= 0x3ul;
  if (TmpState2 == 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = (double) TmpState2;
  i <<= 1;
  TmpState2 = TmpState;
  TmpState2 >>= i;
  TmpState2 &= 0x3ul;
  TmpState &= ~(0x3ul << i);
  switch (TmpState2)
    {
    case 0x1ul:
      TmpState |= 0x2ul << i;      
      break;
    case 0x0ul:
      TmpState |= 0x1ul << i;      
      break;
    }	  
  coefficient = 1.0;
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
  //  return this->FindStateIndexAndTransaltion(state, TmpState, nbrTranslation, i, coefficient);
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Potts3ChainWithTranslations::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->ChainDescription[state];
  unsigned long TmpState2 = TmpState;
  j <<= 1;
  TmpState2 >>= j;
  TmpState2 &= 0x3ul;
  if (TmpState2 == 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = (double) TmpState2;
  i <<= 1;
  TmpState2 = TmpState;
  TmpState2 >>= i;
  TmpState2 &= 0x3ul;
  TmpState &= ~(0x3ul << i);
  switch (TmpState2)
    {
    case 0x2ul:
      TmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      TmpState |= 0x2ul << i;      
      break;
    }	  
  coefficient = 1.0;
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
  //  return this->FindStateIndexAndTransaltion(state, TmpState, nbrTranslation, i, coefficient);
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Potts3ChainWithTranslations::Spi (int i, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->ChainDescription[state];
  unsigned long TmpState2 = TmpState;
  i <<= 1;
  TmpState2 = TmpState;
  TmpState2 >>= i;
  TmpState2 &= 0x3ul;
  TmpState &= ~(0x3ul << i);
  switch (TmpState2)
    {
    case 0x1ul:
      TmpState |= 0x2ul << i;      
      break;
    case 0x0ul:
      TmpState |= 0x1ul << i;      
      break;
    }	  
  coefficient = 1.0;
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
  //  return this->FindStateIndexAndTransaltion(state, TmpState, nbrTranslation, i, coefficient);
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Potts3ChainWithTranslations::Smi (int i, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->ChainDescription[state];
  unsigned long TmpState2 = TmpState;
  i <<= 1;
  TmpState2 = TmpState;
  TmpState2 >>= i;
  TmpState2 &= 0x3ul;
  TmpState &= ~(0x3ul << i);
  switch (TmpState2)
    {
    case 0x2ul:
      TmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      TmpState |= 0x2ul << i;      
      break;
    }	  
  coefficient = 1.0;
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
  //  return this->FindStateIndexAndTransaltion(state, TmpState, nbrTranslation, i, coefficient);
}

// translate a state assuming the system have periodic boundary conditions (increasing the site index)
//
// nbrTranslations = number of translations to apply
// state = index of the state to translate 
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Potts3ChainWithTranslations::TranslateState (int nbrTranslations, int state)
{
  unsigned long TmpState = this->ChainDescription[state];
  TmpState = (((TmpState & (0x3ul << ((this->ChainLength - nbrTranslations) - 1ul)) << 1) << nbrTranslations)
	      | (TmpState >> ((this->ChainLength - nbrTranslations) << 1)));
  return this->FindStateIndex(TmpState);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Potts3ChainWithTranslations::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
//   if (q.GetQuantumNumberType() != AbstractQuantumNumber::Sz)
//     return new Potts3ChainWithTranslations();
//   int TmpSz = ((SzQuantumNumber&) q).GetSz();
//   long LargeHilbertSubspaceDimension = 0l;
//   int* TmpConvArray = new int [this->LargeHilbertSpaceDimension];
//   for (int i = 0; i < this->LargeHilbertSpaceDimension; i++)
//     {
//       if (this->TotalSz(i) == TmpSz)
// 	{
// 	  TmpConvArray[LargeHilbertSubspaceDimension] = i;
// 	  LargeHilbertSubspaceDimension++;	  
// 	}
//     }
//   int* ConvArray = new int [LargeHilbertSubspaceDimension];
//   unsigned long* SubspaceDescription = new unsigned long [LargeHilbertSubspaceDimension];
//   int* SubspaceLookUpTable = new int [this->LookUpTableSize];
//   unsigned long TestMask = this->ChainDescription[TmpConvArray[0]] & this->LookUpTableMask;
//   SubspaceLookUpTable[TestMask] = 0;
//   SubspaceDescription[0] = this->ChainDescription[TmpConvArray[0]];
//   ConvArray[0] = TmpConvArray[0];
//   for (long i = 1; i < LargeHilbertSubspaceDimension; i++)
//     {
//       if ((this->ChainDescription[TmpConvArray[i]] & this->LookUpTableMask) != TestMask)
// 	{
// 	  TestMask = this->ChainDescription[TmpConvArray[i]] & this->LookUpTableMask;
// 	  SubspaceLookUpTable[TestMask] = i;
// 	}
//       SubspaceDescription[i] = this->ChainDescription[TmpConvArray[i]];
//       ConvArray[i] = TmpConvArray[i];
//     }
//   converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, (int) LargeHilbertSubspaceDimension, ConvArray);
//   return (AbstractSpinChain*) new Potts3ChainWithTranslations (LargeHilbertSubspaceDimension, SubspaceDescription, this->ChainLength,
// 							       TmpSz, true, SubspaceLookUpTable, this->LookUpTableSize, 
// 							       this->LookUpPosition, this->LookUpTableMask);
  return 0;
}

// find state index
//
// state = state description
// return value = corresponding index

int Potts3ChainWithTranslations::FindStateIndex(unsigned long state)
{
  int index = 0;//this->LookUpTable[state & this->LookUpTableMask];
  unsigned long* TmpState = &(this->ChainDescription[index]);
  while ((index < this->HilbertSpaceDimension) && (state != *(TmpState++)))
    ++index;
  return index;   
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Potts3ChainWithTranslations::PrintState (ostream& Str, int state)
{
  unsigned long StateDescription = this->ChainDescription[state];  
  for (int j = 0; j < this->ChainLength; ++j)
    {
      Str << (StateDescription & 0x3ul) << " ";
      StateDescription >>= 2;
    }
  return Str;
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix Potts3ChainWithTranslations::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
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
  Potts3Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  int QB = this->Sz - szSector;
  if (QB < 0)
    {
      QB += 3;
    }
  Potts3Chain TmpHilbertSpace(this->ChainLength - nbrSites, QB, 1000000);

  ComplexMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = 2 * nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  int TmpNbrTranslation;
  int TmpNbrTranslationToIdentity;
  Complex* TmpPhases = new Complex [2 * this->ChainLength];
  double Coef = 2.0 * M_PI * ((double) this->Momentum) / ((double) this->ChainLength);
  for (int i = 0; i < (2 * this->ChainLength); ++i)
    {
      TmpPhases[i] = Phase(Coef * ((double) i));
    }

  unsigned long Mask1 = (0x1ul << Shift) - 0x1ul;
  unsigned long Mask2 = (0x1ul << (2 * this->ChainLength)) - 0x1ul;
  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = (TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask1)) & Mask2;
	  double Coefficient = 1.0;
	  int TmpPos = this->SymmetrizeResult(TmpState2, 1, Coefficient, TmpNbrTranslation);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos] * TmpPhases[TmpNbrTranslation] / sqrt((double) this->NbrStateInOrbit[TmpPos]));
	    }
	}
    }
  delete[] TmpPhases;
  return TmpEntanglementMatrix;
}

