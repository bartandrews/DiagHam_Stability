////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of spin 1 chain with the Sz<->-Sz symmetry             //
//                                                                            //
//                        last modification : 11/05/2017                      //
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


#include "HilbertSpace/Spin1ChainWithSzSymmetry.h"
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

Spin1ChainWithSzSymmetry::Spin1ChainWithSzSymmetry () 
{
  this->SzSymmetrySector = 1.0;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// memorySize = memory size in bytes allowed for look-up table

Spin1ChainWithSzSymmetry::Spin1ChainWithSzSymmetry (int chainLength, int szSymmetrySector, unsigned long memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;

  this->LargeHilbertSpaceDimension = 3l;
  for (int i = 1; i < this->ChainLength; i++)
    this->LargeHilbertSpaceDimension *= 3l;

  if (szSymmetrySector == 1)
    {
      this->SzSymmetrySector = 1.0;
    }
  else
    {
      this->SzSymmetrySector = -1.0;
    }
  this->SzSymmetryMask = (0x1ul << (2 * this-> ChainLength)) - 0x1ul;

  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
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

Spin1ChainWithSzSymmetry::Spin1ChainWithSzSymmetry (int chainLength, int szSymmetrySector, int sz, unsigned long memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz;
  this->FixedSpinProjectionFlag = true;

  if (szSymmetrySector == 1)
    {
      this->SzSymmetrySector = 1.0;
    }
  else
    {
      this->SzSymmetrySector = -1.0;
    }
  this->SzSymmetryMask = (0x1ul << (2 * this-> ChainLength)) - 0x1ul;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(0, this->ChainLength);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1, 0);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
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

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1ChainWithSzSymmetry::Spin1ChainWithSzSymmetry (const Spin1ChainWithSzSymmetry& chain) 
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

      this->SzSymmetryMask = chain.SzSymmetryMask;
      this->SzSymmetrySector = chain.SzSymmetrySector;

      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;

      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
    }
  else
    {
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0l;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;

      this->SzSymmetryMask = 0x0ul;
      this->SzSymmetrySector = 0.0;

      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;

      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
    }
}

// destructor
//

Spin1ChainWithSzSymmetry::~Spin1ChainWithSzSymmetry () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1ChainWithSzSymmetry& Spin1ChainWithSzSymmetry::operator = (const Spin1ChainWithSzSymmetry& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTable;
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
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

      this->SzSymmetryMask = chain.SzSymmetryMask;
      this->SzSymmetrySector = chain.SzSymmetrySector;

      this->LookUpTable = chain.LookUpTable;
      this->MaximumLookUpShift = chain.MaximumLookUpShift;
      this->LookUpTableMemorySize = chain.LookUpTableMemorySize;
      this->LookUpTableShift = chain.LookUpTableShift;

      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
    }
  else
    {
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0l;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
 
      this->SzSymmetryMask = 0x0ul;
      this->SzSymmetrySector = 0.0;

      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;

      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1ChainWithSzSymmetry::Clone()
{
  return new Spin1ChainWithSzSymmetry (*this);
}

// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long Spin1ChainWithSzSymmetry::GenerateStates()
{
  long TmpLargeHilbertSpaceDimension = 0l;
#ifdef  __64_BITS__
  unsigned long Discard = 0xfffffffffffffffful;
#else
  unsigned long Discard = 0xfffffffful;
#endif
  double TmpSign;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], TmpSign) == this->StateDescription[i]))
	{
	  if (this->TestDiscreteSymmetryConstraint(this->StateDescription[i]) == true)
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
  this->ComputeRescalingFactors();
  return TmpLargeHilbertSpaceDimension;
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainWithSzSymmetry::Spi (int i, int state, double& coefficient)
{
  unsigned long State = (this->StateDescription[state] >> (i << 1));
  unsigned long tmpState = State;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x2ul:
      {
	coefficient = M_SQRT2;
	return this->SymmetrizeResult(State | (0x1ul << (i << 1)), this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT2;
	return this->SymmetrizeResult(State | (0x2ul << (i << 1)), this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x3ul:
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

int Spin1ChainWithSzSymmetry::Smi (int i, int state, double& coefficient)
{
  unsigned long State = (this->StateDescription[state] >> (i << 1));
  unsigned long tmpState = State;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x2ul:
      {
	coefficient = M_SQRT2;
	return this->SymmetrizeResult(State & ~(0x2ul << (i << 1)), this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x0ul:
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x3ul:
      {
	coefficient = M_SQRT2;
	return this->SymmetrizeResult(State & ~(0x1ul << (i << 1)), this->NbrStateInOrbit[state], coefficient);
      }
      break;
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1ChainWithSzSymmetry::SmiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      coefficient = M_SQRT2;
      State &= ~(0x1ul << i);
      break;
    case 0x2ul:
      coefficient = M_SQRT2;
      State&= ~(0x2ul << i);
      break;
    case 0x0ul:
      return this->HilbertSpaceDimension;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  switch (tmpState2)
    {
    case 0x3ul:
      return this->HilbertSpaceDimension;
    case 0x2ul:
      coefficient *= M_SQRT2;
      return this->SymmetrizeResult(State | (0x1ul << j), this->NbrStateInOrbit[state], coefficient);
    case 0x0ul:
      coefficient *= M_SQRT2;
      return this->SymmetrizeResult(State | (0x2ul << j), this->NbrStateInOrbit[state], coefficient);
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

int Spin1ChainWithSzSymmetry::SpiSpj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      return this->HilbertSpaceDimension;
      break;
    case 0x2ul:
      coefficient = M_SQRT2;
      State |= (0x1ul << i);
      break;
    case 0x0ul:
      coefficient = M_SQRT2;
      State |= (0x2ul << i);
      break;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  switch (tmpState2)
    {
    case 0x3ul:
      return this->HilbertSpaceDimension;
    case 0x2ul:
      coefficient *= M_SQRT2;
      return this->SymmetrizeResult(State | (0x1ul << j), this->NbrStateInOrbit[state], coefficient);
    case 0x0ul:
      coefficient *= M_SQRT2;
      return this->SymmetrizeResult(State | (0x2ul << j), this->NbrStateInOrbit[state], coefficient);
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

int Spin1ChainWithSzSymmetry::SmiSmj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      coefficient = M_SQRT2;
      State &= ~(0x1ul << i);
      break;
    case 0x2ul:
      coefficient = M_SQRT2;
      State&= ~(0x2ul << i);
      break;
    case 0x0ul:
      return this->HilbertSpaceDimension;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  switch (tmpState2)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State & ~(0x1ul << j), this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State & ~(0x2ul << j), this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x0ul:
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

int Spin1ChainWithSzSymmetry::SpiSzj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = (tmpState >> (j << 1)) & 0x3ul;
  if (tmpState2 == 0x2ul)
    {
      coefficient = 0;
      return this->HilbertSpaceDimension;
    }
  if (tmpState2 == 0x3ul)
    {
      coefficient = 1.0;
    }
  else
    {
      coefficient = -1.0;
    }
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      {
	coefficient = 0;
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State | (0x1ul << i), this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State | (0x2ul << i), this->NbrStateInOrbit[state], coefficient);
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

int Spin1ChainWithSzSymmetry::SmiSzj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = (tmpState >> (j << 1)) & 0x3ul;
  if (tmpState2 == 0x2ul)
    {
      coefficient = 0;
      return this->HilbertSpaceDimension;
    }
  if (tmpState2 == 0x3ul)
    {
      coefficient = 1.0;
    }
  else
    {
      coefficient = -1.0;
    }
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State & ~(0x1ul << i), this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT2;
	return this->SymmetrizeResult(State & ~(0x1ul << i), this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x0ul:
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
// return value = index of resulting state (orbit index)

int Spin1ChainWithSzSymmetry::SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[state];
  j2 <<= 1;
  switch ((TmpState >> j2) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT2;
	TmpState |= (0x1ul << j2);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT2;
 	TmpState |= (0x2ul << j2);
      }
      break;
    }
  i2 <<= 1;
  switch ((TmpState >> i2) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT2;
 	TmpState &= ~(0x1ul << i2);
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT2;
	TmpState &= ~(0x2ul << i2);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    }
  j1 <<= 1;
  switch ((TmpState >> j1) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT2;
	TmpState |= (0x1ul << j1);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT2;
 	TmpState |= (0x2ul << j1);
      }
      break;
    }
  i1 <<= 1;
  switch ((TmpState >> i1) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT2;
	TmpState &= ~(0x1ul << i1);
	return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT2;
	TmpState &= ~(0x2ul << i1);
	return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x0ul:
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
// return value = index of resulting state (orbit index)

int Spin1ChainWithSzSymmetry::SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[state];
  j2 <<= 1;
  switch ((TmpState >> j2) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT2;
	TmpState |= (0x1ul << j2);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT2;
 	TmpState |= (0x2ul << j2);
      }
      break;
    }
  i2 <<= 1;
  switch ((TmpState >> i2) & 0x3ul)
    {
    case 0x3ul:
      {
	TmpState &= ~(0x1ul << i2);
	unsigned long TmpState2 = (TmpState >> (i1 << 1)) & 0x3ul;
	unsigned long TmpState3 = (TmpState >> (j1 << 1)) & 0x3ul;
	if ((TmpState2 == 0x2ul) || (TmpState3 == 0x2ul))
	  {
	    coefficient = 0.0;
	    return this->HilbertSpaceDimension;
	  }
	if (TmpState2 == TmpState3)
	  coefficient *= M_SQRT2;
	else
	  coefficient *= -M_SQRT2;
	return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x2ul:
      {
	TmpState &= ~(0x2ul << i2);
	unsigned long TmpState2 = (TmpState >> (i1 << 1)) & 0x3ul;
	unsigned long TmpState3 = (TmpState >> (j1 << 1)) & 0x3ul;
	if ((TmpState2 == 0x2ul) || (TmpState3 == 0x2ul))
	  {
	    coefficient = 0.0;
	    return this->HilbertSpaceDimension;
	  }
	if (TmpState2 == TmpState3)
	  coefficient *= M_SQRT2;
	else
	  coefficient *= -M_SQRT2;
	return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// compute the rescaling factors
//

void Spin1ChainWithSzSymmetry::ComputeRescalingFactors()
{
  int Tmp = 4;
  this->RescalingFactors = new double* [Tmp + 1];
  for (int i = 1; i <= Tmp; ++i)
    {
      this->RescalingFactors[i] = new double [Tmp + 1];
      for (int j = 1; j <= Tmp; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}


// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix Spin1ChainWithSzSymmetry::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
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
  Spin1Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);
  RealMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int Shift = 2 * nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | TmpDestinationHilbertSpace.StateDescription[j];
	  double Coefficient = 1.0;
	  int TmpPos = this->SymmetrizeResult(TmpState2, 1, Coefficient);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	       TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos] * Coefficient);
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

ComplexMatrix Spin1ChainWithSzSymmetry::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
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
  Spin1Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);

  ComplexMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = 2 * nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | TmpDestinationHilbertSpace.StateDescription[j];
	  double Coefficient = 1.0;
	  int TmpPos = this->SymmetrizeResult(TmpState2, 1, Coefficient);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	       TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos] * Coefficient);
	    }
	}
    }

   return TmpEntanglementMatrix;
}
