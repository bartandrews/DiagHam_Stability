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
//                and a fixed parity under the mirror symmetry                //
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


#include "HilbertSpace/Spin1_2ChainMirrorSymmetry.h"
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

Spin1_2ChainMirrorSymmetry::Spin1_2ChainMirrorSymmetry ()
{
  this->MirrorParity = 0;
  this->MirrorParitySign = 1.0;
}

// constructor for complete Hilbert space with a given total spin projection Sz
//
// chainLength = number of spin 1/2
// sz = twice the value of total Sz component
// parity = parity of the total (Sz + 1/2)
// memorySize = memory size in bytes allowed for look-up table

Spin1_2ChainMirrorSymmetry::Spin1_2ChainMirrorSymmetry (int chainLength, int sz, int parity, int memorySize) 
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
  this->MirrorParity = parity;
  if (this->MirrorParity == 0)
    this->MirrorParitySign = 1.0;
  else
    this->MirrorParitySign = -1.0;
#ifdef __64_BITS__
  if ((this->ChainLength & 1) == 0)
    {
      this->MirrorShift = 32 - (this->ChainLength >> 1);
      this->MirrorUnshift = this->MirrorShift;
    }
  else
    {
      this->MirrorShift = 32 - ((this->ChainLength - 1)  >> 1);
      this->MirrorUnshift = this->MirrorShift - 1;
    }
#else
  if ((this->ChainLength & 1) == 0)
    {
      this->MirrorShift = 16 - (this->ChainLength >> 1);
      this->MirrorUnshift = this->MirrorShift;
    }
  else
    {
      this->MirrorShift = 16 - ((this->ChainLength - 1)  >> 1);
      this->MirrorUnshift = this->MirrorShift - 1;
    }
#endif
  unsigned long TmpLargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->ChainLength);
  this->StateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = this->GenerateStates ((this->Sz + this->ChainLength) >> 1, this->ChainLength - 1, 0l);
  this->LargeHilbertSpaceDimension = 0ul;
  if (this->MirrorParity == 0)
    {
      for (long i = 0l; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->ApplyMirrorSymmetry(this->StateDescription[i]);
	  if (TmpState >= this->StateDescription[i])
	    {
	      ++this->LargeHilbertSpaceDimension;
	    }
	}
      unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->LargeHilbertSpaceDimension = 0ul;
      for (long i = 0l; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->ApplyMirrorSymmetry(this->StateDescription[i]);
	  if (TmpState >= this->StateDescription[i])
	    {
	      if (TmpState != this->StateDescription[i])
		{
		  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->StateDescription[i] | SPIN1_2CHAIN_MIRRORSYMMETRIC_BIT;
		}
	      else
		{
		  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->StateDescription[i];
		}
	      ++this->LargeHilbertSpaceDimension;	  
	    }
	}
      delete[] this->StateDescription;
      this->StateDescription = TmpStateDescription;
    }
  else
    {
      for (long i = 0l; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->ApplyMirrorSymmetry(this->StateDescription[i]);
	  if (TmpState > this->StateDescription[i])
	    {
	      ++this->LargeHilbertSpaceDimension;
	    }
	}
      if (this->LargeHilbertSpaceDimension > 0)
	{
	  unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
	  this->LargeHilbertSpaceDimension = 0ul;
	  for (long i = 0l; i < TmpLargeHilbertSpaceDimension; ++i)
	    {
	      unsigned long TmpState = this->ApplyMirrorSymmetry(this->StateDescription[i]);
	      if (TmpState > this->StateDescription[i])
		{
		  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->StateDescription[i] | SPIN1_2CHAIN_MIRRORSYMMETRIC_BIT;
		  ++this->LargeHilbertSpaceDimension;	  
		}
	    }
	  delete[] this->StateDescription;
	  this->StateDescription = TmpStateDescription;
	}
      else
	{
	  this->StateDescription = 0;
	}
    }
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->MaxStateDescription = this->StateDescription[0l] & SPIN1_2CHAIN_FULLSYMMETRY_MASK;
      this->MinStateDescription = this->StateDescription[this->LargeHilbertSpaceDimension - 1l] & SPIN1_2CHAIN_FULLSYMMETRY_MASK;
      this->GenerateLookUpTable(memorySize, SPIN1_2CHAIN_FULLSYMMETRY_MASK);
    }
  this->Flag.Initialize();
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainMirrorSymmetry::Spin1_2ChainMirrorSymmetry (const Spin1_2ChainMirrorSymmetry& chain)
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
      this->MirrorParity = chain.MirrorParity;
      this->MirrorParitySign = chain.MirrorParitySign;
      this->MirrorShift = chain.MirrorShift;
      this->MirrorUnshift = chain.MirrorUnshift;
      this->MaxStateDescription = chain.MaxStateDescription;
      this->MinStateDescription = chain.MinStateDescription;
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
      this->MirrorParity = 0;
      this->MirrorParitySign = 1.0;
      this->MirrorShift = 0;
      this->MirrorUnshift = 0;
      this->MaxStateDescription = 0x0ul;
      this->MinStateDescription = 0x0ul;
    }
}

// destructor
//

Spin1_2ChainMirrorSymmetry::~Spin1_2ChainMirrorSymmetry () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainMirrorSymmetry& Spin1_2ChainMirrorSymmetry::operator = (const Spin1_2ChainMirrorSymmetry& chain)
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
      this->MirrorParity = chain.MirrorParity;
      this->MirrorParitySign = chain.MirrorParitySign;
      this->MirrorShift = chain.MirrorShift;
      this->MirrorUnshift = chain.MirrorUnshift;
      this->MaxStateDescription = chain.MaxStateDescription;
      this->MinStateDescription = chain.MinStateDescription;
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
      this->MirrorParity = 0;
      this->MirrorParitySign = 1.0;
      this->MirrorShift = 0;
      this->MirrorUnshift = 0;
      this->MaxStateDescription = 0x0ul;
      this->MinStateDescription = 0x0ul;
   }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainMirrorSymmetry::Clone()
{
  return new Spin1_2ChainMirrorSymmetry (*this);
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainMirrorSymmetry::Spi (int i, int state, double& coefficient)
{
  unsigned long State = this->StateDescription[state];
  unsigned long tmpState = (State >> i) & 0x1ul;
  if (tmpState == 0x1ul)
    {
      this->SOperatorSignature = State & SPIN1_2CHAIN_FULLSYMMETRY_BIT;
      State &= SPIN1_2CHAIN_FULLSYMMETRY_MASK;
      coefficient = 1.0;
      return this->SymmetrizeResult(State & ~(0x1ul << i), coefficient);
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

int Spin1_2ChainMirrorSymmetry::Pij (int i, int j, int state)
{  
  cout << "warning, Spin1_2ChainMirrorSymmetry::Pij  should be used" << endl;
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpMask = (0x1ul << i) | (0x1ul << j);
  unsigned long tmpState2 = tmpState & tmpMask;
  unsigned long tmpState3 = ~tmpState & tmpMask;
  if ((tmpState2 == 0x0ul) || (tmpState3 == 0x0ul))
    {
      return this->HilbertSpaceDimension;
    }
  else
    {
      this->SOperatorSignature = tmpState & SPIN1_2CHAIN_FULLSYMMETRY_BIT;
      tmpState &= SPIN1_2CHAIN_FULLSYMMETRY_MASK;
      double Coefficient = 1.0;
      return this->SymmetrizeResult((tmpState & ~tmpMask) | tmpState3, Coefficient);
    }
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainMirrorSymmetry::Smi (int i, int state, double& coefficient)
{
  unsigned long State = this->StateDescription[state];
  unsigned long tmpState = (State >> i) & 0x1ul;
  if (tmpState != 0x0ul)
    {
      this->SOperatorSignature = State & SPIN1_2CHAIN_FULLSYMMETRY_BIT;
      State &= SPIN1_2CHAIN_FULLSYMMETRY_MASK;
      coefficient = 1.0;
      return this->SymmetrizeResult(State ^ (0x1ul << i), coefficient);
    }
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainMirrorSymmetry::SmiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  this->SOperatorSignature = tmpState & SPIN1_2CHAIN_FULLSYMMETRY_BIT;
  tmpState &= SPIN1_2CHAIN_FULLSYMMETRY_MASK;
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
	  return this->SymmetrizeResult((State | (0x1ul << j)) & ~(0x1ul << i), coefficient);
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

int Spin1_2ChainMirrorSymmetry::SpiSmj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  this->SOperatorSignature = tmpState & SPIN1_2CHAIN_FULLSYMMETRY_BIT;
  tmpState &= SPIN1_2CHAIN_FULLSYMMETRY_MASK;
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
	  return this->SymmetrizeResult((State | (0x1ul << i)) & ~(0x1ul << j), coefficient);
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

int Spin1_2ChainMirrorSymmetry::SpiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  this->SOperatorSignature = tmpState & SPIN1_2CHAIN_FULLSYMMETRY_BIT;
  tmpState &= SPIN1_2CHAIN_FULLSYMMETRY_MASK;
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
	  return this->SymmetrizeResult(State | ((0x1ul << j) | (0x1ul << i)), coefficient);
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

int Spin1_2ChainMirrorSymmetry::SmiSmj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  this->SOperatorSignature = tmpState & SPIN1_2CHAIN_FULLSYMMETRY_BIT;
  tmpState &= SPIN1_2CHAIN_FULLSYMMETRY_MASK;
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
	  return this->SymmetrizeResult(State & ~((0x1ul << j) | (0x1ul << i)), coefficient);
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

int Spin1_2ChainMirrorSymmetry::SpiSzj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  this->SOperatorSignature = tmpState & SPIN1_2CHAIN_FULLSYMMETRY_BIT;
  tmpState &= SPIN1_2CHAIN_FULLSYMMETRY_MASK;
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
      return this->SymmetrizeResult(State | (0x1ul << i), coefficient);
      
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

int Spin1_2ChainMirrorSymmetry::SmiSzj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  this->SOperatorSignature = tmpState & SPIN1_2CHAIN_FULLSYMMETRY_BIT;
  tmpState &= SPIN1_2CHAIN_FULLSYMMETRY_MASK;
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
      return this->SymmetrizeResult(State & ~(0x1ul << i), coefficient);
      
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

int Spin1_2ChainMirrorSymmetry::SpiSmjSzk (int i, int j, int k, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  this->SOperatorSignature = tmpState & SPIN1_2CHAIN_FULLSYMMETRY_BIT;
  tmpState &= SPIN1_2CHAIN_FULLSYMMETRY_MASK;
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
	  return this->SymmetrizeResult((State | (0x1ul << i)) & ~(0x1ul << j), coefficient);
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

// compute the parity (prod_i Sz_i) for a given state
//
// state = index of the state to be applied on Sz_i operator
// return value = 0 if prod_i Sz_i = 1, 1 if prod_i Sz_i = -1

unsigned long Spin1_2ChainMirrorSymmetry::GetParity (int state)
{
  unsigned long TmpState = this->StateDescription[state] & SPIN1_2CHAIN_FULLSYMMETRY_MASK;
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

int Spin1_2ChainMirrorSymmetry::TranslateState (int nbrTranslations, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  TmpState = (((TmpState & ((0x1ul << (this->ChainLength - nbrTranslations)) - 1ul)) << nbrTranslations)
	      | (TmpState >> (this->ChainLength - nbrTranslations)));
  double Coefficient = 1.0;
  return this->SymmetrizeResult(TmpState, Coefficient);
}


// find state index
//
// stateDescription = state description
// return value = corresponding index

int Spin1_2ChainMirrorSymmetry::FindStateIndex(unsigned long stateDescription)
{
  if ((stateDescription > this->MaxStateDescription) || (stateDescription < this->MinStateDescription))
    {
      return this->HilbertSpaceDimension;
    }
  int CurrentMaximumUpPosition = this->ChainLength - 1;
  while ((((stateDescription >> CurrentMaximumUpPosition) & 0x1ul) == 0x0ul) && (CurrentMaximumUpPosition > 0))
    --CurrentMaximumUpPosition;
  long PosMax = stateDescription >> this->LookUpTableShift[CurrentMaximumUpPosition];
  long PosMin = this->LookUpTable[CurrentMaximumUpPosition][PosMax];
  PosMax = this->LookUpTable[CurrentMaximumUpPosition][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid] & SPIN1_2CHAIN_FULLSYMMETRY_MASK;
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
      CurrentState = this->StateDescription[PosMid] & SPIN1_2CHAIN_FULLSYMMETRY_MASK;
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    {
      if (((this->StateDescription[PosMin] & SPIN1_2CHAIN_FULLSYMMETRY_MASK) != stateDescription) 
	  && ((this->StateDescription[PosMax] & SPIN1_2CHAIN_FULLSYMMETRY_MASK) != stateDescription))
	return this->HilbertSpaceDimension;
      else
	return PosMin;
    }
}

// get the normalization factor in front of each basis state (i.e. 1/sqrt(orbit size))
//
// return value = pointer to normalization factors

double* Spin1_2ChainMirrorSymmetry::GetBasisNormalization()
{
  double* TmpNorm = new double[this->HilbertSpaceDimension];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if ((this->StateDescription[i] & SPIN1_2CHAIN_MIRRORSYMMETRIC_BIT) != 0x0ul)
	{
	  TmpNorm[i] = 1.0 / M_SQRT2;
	}
      else
	{
	  TmpNorm[i] = 1.0;
	}
    }
  return TmpNorm;
}
 
