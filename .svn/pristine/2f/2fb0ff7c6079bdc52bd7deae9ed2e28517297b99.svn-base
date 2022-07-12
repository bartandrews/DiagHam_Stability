////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with a fixed Sz value                //
//                                                                            //
//                        last modification : 29/06/2015                      //
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


#include "HilbertSpace/Spin1_2ChainNew.h"
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

Spin1_2ChainNew::Spin1_2ChainNew ()
{
}

// constructor for complete Hilbert space with a given total spin projection Sz
//
// chainLength = number of spin 1/2
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Spin1_2ChainNew::Spin1_2ChainNew (int chainLength, int sz, int memorySize) 
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
  this->FixedQuantumNumberFlag = true;
  this->Sz = sz;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->ChainLength);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = this->GenerateStates ((this->Sz + this->ChainLength) >> 1, this->ChainLength - 1, 0l);
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->GenerateLookUpTable(memorySize);
  this->Flag.Initialize();
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainNew::Spin1_2ChainNew (const Spin1_2ChainNew& chain)
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
    }
}

// destructor
//

Spin1_2ChainNew::~Spin1_2ChainNew () 
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->LookUpTableShift != 0)
	{
	  delete[] this->LookUpTableShift;
	  for (int i = 0; i < this->ChainLength; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;
	}
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainNew& Spin1_2ChainNew::operator = (const Spin1_2ChainNew& chain)
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

AbstractHilbertSpace* Spin1_2ChainNew::Clone()
{
  return new Spin1_2ChainNew (*this);
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin1_2ChainNew::TotalSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->Sz;
  unsigned long State = this->StateDescription[index];
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpSz += ((State & 0x1ul) << 1);
      State >>= 1;
    }
  TmpSz -= this->ChainLength;
  return TmpSz;
}
// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainNew::Spi (int i, int state, double& coefficient)
{
  unsigned long State = this->StateDescription[state];
  unsigned long tmpState = (State >> i) & 0x1ul;
  if (tmpState == 0x1ul)
    {
      coefficient = 1.0;
      return this->FindStateIndex(State & ~(0x1ul << i));
    }
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of P_ij operator on a given state
//
// i = first position
// j = second position
// state = index of the state to be applied on P_ij operator
// return value = index of resulting state

int Spin1_2ChainNew::Pij (int i, int j, int state)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpMask = (0x1ul << i) | (0x1ul << j);
  unsigned long tmpState2 = tmpState & tmpMask;
  unsigned long tmpState3 = ~tmpState & tmpMask;
  if ((tmpState2 == 0x0ul) || (tmpState3 == 0x0ul))
    return this->HilbertSpaceDimension;
  else
    return this->FindStateIndex((tmpState & ~tmpMask) | tmpState3);
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double Spin1_2ChainNew::SziSzj (int i, int j, int state)
{  
  unsigned long Mask = ((0x1ul << i) | (0x1ul << j));
  unsigned long tmpState = this->StateDescription[state] & Mask;
  if ((tmpState == 0x0ul) || (tmpState == Mask))
    return 0.25;
  else
    return -0.25;
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainNew::Smi (int i, int state, double& coefficient)
{
  unsigned long State = this->StateDescription[state];
  unsigned long tmpState = (State >> i) & 0x1ul;
  if (tmpState != 0x0ul)
    {
      coefficient = 1.0;
      return this->FindStateIndex(State ^ (0x1ul << i));
    }
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i operator on a given state
//
// i = position of Sz operator
// state = index of the state to be applied on Sz_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainNew::Szi (int i, int state, double& coefficient)
{
  unsigned long tmpState = (this->StateDescription[state] >> i) & 0x1ul;
  if (tmpState == 0x0ul)
    coefficient = -0.5;
  else
    coefficient = 0.5;
  return state;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainNew::SmiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex((State | (0x1ul << j)) & ~(0x1ul << i));
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
      return state;
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i S-_j operator on a given state
//
// i = position of S+ operator
// j = position of S- operator
// state = index of the state to be applied on S+_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainNew::SpiSmj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= j;
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= i; 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex((State | (0x1ul << i)) & ~(0x1ul << j));
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0x1ul)
    {
      coefficient = -0.25;
      return state;
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

int Spin1_2ChainNew::SpiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x3ul)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex(State | ((0x1ul << j) | (0x1ul << i)));
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainNew::SmiSmj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x0ul)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex(State & ~((0x1ul << j) | (0x1ul << i)));
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainNew::SpiSzj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x1ul;
  if (tmpState == 0x0ul)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x1ul;
      if (tmpState2 == 0x0ul)
	coefficient = -0.5;
      else
	coefficient = 0.5;
      return this->FindStateIndex(State | (0x1ul << i));
      
    }
  else
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
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

int Spin1_2ChainNew::SmiSzj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x1ul;
  if (tmpState == 0x1ul)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x1ul;
      if (tmpState2 == 0x0ul)
	coefficient = -0.5;
      else
	coefficient = 0.5;
      return this->FindStateIndex(State & ~(0x1ul << i));
      
    }
  else
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i S-_j Sz_k operator on a given state
//
// i = position of S+ operator
// j = position of S- operator
// k = position of Sz operator
// state = index of the state to be applied on S+_i S-_j Sz_k operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainNew::SpiSmjSzk (int i, int j, int k, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  if (((tmpState >> k) & 0x1ul) == 0x0ul)
    coefficient = -0.5;
  else
    coefficient = 0.5;
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= j;
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= i; 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  return this->FindStateIndex((State | (0x1ul << i)) & ~(0x1ul << j));
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0x1ul)
    {
      coefficient *= -0.25;
      return state;
    }
  return this->HilbertSpaceDimension;
}

// operate local isometry on three sites
//
// i = position of first site
// j = position of second site 
// k = position of third site
// state = index of the state that the isometry has to be applied on
// indices = reference to an array where the indices of the resulting states have to be stored
// coefficients = reference to the array where the coefficients have to be stored
// return value = number of non-zero coefficients
  
int Spin1_2ChainNew::ThreeSiteIsometry (int i, int j, int k, int state, int*& indices, double*& coefficients)
{
  unsigned long tmpStateI = (this->StateDescription[state] >> i) & 0x1ul;
  unsigned long tmpStateJ = (this->StateDescription[state] >> j) & 0x1ul;
  unsigned long tmpStateK = (this->StateDescription[state] >> k) & 0x1ul;
  
  int localSz = (int) (tmpStateI + tmpStateJ + tmpStateK);
  if ((localSz == 0) || (localSz == 3))
    return 0;
  
  if (tmpStateJ == tmpStateK)
  {
    if (localSz == 2)
      indices[0] = this->SpiSmj (i, j, state, coefficients[0]);
    else
      indices[0] = this->SpiSmj (j, i, state, coefficients[0]);
    coefficients[0] = -2.0 / sqrt(6.0);    
    return 1;
  }
  
  if (tmpStateI == tmpStateJ)
  {
    indices[0] = state;
    coefficients[0] = 1.0 / sqrt(2.0);
    if (localSz == 2)
      indices[1] = this->SpiSmj (k, j, state, coefficients[1]);
    else
      indices[1] = this->SpiSmj (j, k, state, coefficients[1]);
    coefficients[1] = 1.0 / sqrt(6.0);
  }
  else
  {
    indices[0] = state;
    coefficients[0] = 1.0 / sqrt(6.0);
    if (localSz == 2)
      indices[1] = this->SpiSmj (j, k, state, coefficients[1]);
    else
      indices[1] = this->SpiSmj (k, j, state, coefficients[1]);
    coefficients[1] = -1.0 / sqrt(2.0);
  }
  return 2;
  
}



// compute the parity (prod_i Sz_i) for a given state
//
// state = index of the state to be applied on Sz_i operator
// return value = 0 if prod_i Sz_i = 1, 1 if prod_i Sz_i = -1

unsigned long Spin1_2ChainNew::GetParity (int state)
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

int Spin1_2ChainNew::TranslateState (int nbrTranslations, int state)
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

int Spin1_2ChainNew::FindStateIndex(unsigned long stateDescription)
{
  int CurrentMaximumUpPosition = this->ChainLength - 1;
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

long Spin1_2ChainNew::EvaluateHilbertSpaceDimension(int sz, int nbrSites)
{
  BinomialCoefficients TmpCoefficients (nbrSites);
  return TmpCoefficients(nbrSites, (nbrSites + sz) >> 1);
}

// generate all states
// 
// nbrSpinUp = number of spin up
// currentPosition = current position to consider in the chain
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long Spin1_2ChainNew::GenerateStates(int nbrSpinUp, int currentPosition, long pos)
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
  long TmpPos = this->GenerateStates(nbrSpinUp - 1, ReducedCurrentPosition, pos);
  unsigned long Mask = 0x1ul << currentPosition;
  for (long i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  return this->GenerateStates(nbrSpinUp, ReducedCurrentPosition, TmpPos);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin1_2ChainNew::PrintState (ostream& Str, int state)
{  
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long Mask = 0x1ul;
  for (int k = 0; k < this->ChainLength; k++)    
    {
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
// stateMask = an optional mask to apply to each state to focus on the relevant bits

void Spin1_2ChainNew::GenerateLookUpTable(unsigned long memory, unsigned long stateMask)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->ChainLength);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->ChainLength)
    this->MaximumLookUpShift = this->ChainLength;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->ChainLength];
  this->LookUpTableShift = new int [this->ChainLength];
  for (int i = 0; i < this->ChainLength; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentPosition = this->ChainLength - 1;
  while ((((this->StateDescription[0] >> CurrentPosition) & 0x1ul) == 0x0ul) && (CurrentPosition > 0))
    --CurrentPosition;
  int* TmpLookUpTable = this->LookUpTable[CurrentPosition];
  if (CurrentPosition < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentPosition] = 0;
  else
    this->LookUpTableShift[CurrentPosition] = CurrentPosition + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentPosition];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = (this->StateDescription[0] & stateMask) >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      int TmpCurrentPosition = CurrentPosition;
      while ((((this->StateDescription[i] >> TmpCurrentPosition) & 0x1ul) == 0x0ul) && (TmpCurrentPosition > 0))
	--TmpCurrentPosition;
      if (CurrentPosition != TmpCurrentPosition)
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
 	  CurrentPosition = TmpCurrentPosition;
	  TmpLookUpTable = this->LookUpTable[CurrentPosition];
	  if (CurrentPosition < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentPosition] = 0;
	  else
	    this->LookUpTableShift[CurrentPosition] = CurrentPosition + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentPosition];
	  TmpLookUpTableValue = (this->StateDescription[i] & stateMask) >> CurrentShift;
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
	  TmpLookUpTableValue = (this->StateDescription[i] & stateMask) >> CurrentShift;
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



// return the Bosonic Occupation of a given state in the basis
//
// index = index of the state in the basis
// finalState = reference on the array where the monomial representation has to be stored

void Spin1_2ChainNew::GetBosonicOccupation (unsigned int index, int * finalState)
{
  for (int i = 0; i < this->ChainLength; i++)
    {
      finalState[i] = (this->StateDescription[index] >> ((unsigned long) i) )& 0x1ul;
    }
}

// convert the state on the site to its binary representation
//
// state = state to be stored
// sitePosition = position on the chain of the state
// return integer that code the state

unsigned long Spin1_2ChainNew::EncodeSiteState(int physicalState, int sitePosition)
{
  return  physicalState << sitePosition;
}
