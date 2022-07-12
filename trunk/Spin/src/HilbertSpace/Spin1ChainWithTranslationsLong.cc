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


#include "HilbertSpace/Spin1ChainWithTranslationsLong.h"
#include "HilbertSpace/Spin1ChainLong.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <math.h>

using std::cout;
using std::endl;
using std::dec;
using std::hex;



// default constructor
//

Spin1ChainWithTranslationsLong::Spin1ChainWithTranslationsLong () 
{
  this->Flag.Initialize();
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->ChainLength = 0;
  this->Momentum = 0;
  this->Sz = 0;
  this->FixedSpinProjectionFlag = false;
  this->RescalingFactors = 0;
  this->NbrStateInOrbit = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  this->MaxXMomentum = 0;
  this->StateXShift = 0;
  this->XMomentumMask = ((ULONGLONG) 0x0ul);
  this->ComplementaryStateXShift = 0;
}

// constructor for Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// memory = amount of memory granted for precalculations

Spin1ChainWithTranslationsLong::Spin1ChainWithTranslationsLong (int chainLength, int momentum, unsigned long memory) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->Momentum = momentum;

  this->MaxXMomentum = this->ChainLength;
  this->StateXShift = 2 * (this->ChainLength / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this-> ChainLength) - this->StateXShift;
  this->XMomentumMask = (((ULONGLONG) 0x1ul) << this->StateXShift) - ((ULONGLONG) 0x1ul);

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->ChainLength);
  this->StateDescription = new ULONGLONG [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memory);
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

// constructor for Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memory = amount of memory granted for precalculations

Spin1ChainWithTranslationsLong::Spin1ChainWithTranslationsLong (int chainLength, int momentum, int sz, unsigned long memory) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz;
  this->FixedSpinProjectionFlag = true;
  this->Momentum = momentum;

  this->MaxXMomentum = this->ChainLength;
  this->StateXShift = 2 * (this->ChainLength / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this-> ChainLength) - this->StateXShift;
  this->XMomentumMask = (((ULONGLONG) 0x1ul) << this->StateXShift) - ((ULONGLONG) 0x1ul);

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(0, this->ChainLength);
  this->StateDescription = new ULONGLONG [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1, 0);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memory);
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

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1ChainWithTranslationsLong::Spin1ChainWithTranslationsLong (const Spin1ChainWithTranslationsLong& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->MaxXMomentum = chain.MaxXMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;

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
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;

      this->MaxXMomentum = 0;
      this->StateXShift = 0;
      this->ComplementaryStateXShift = 0;
      this->XMomentumMask = ((ULONGLONG) 0x0ul);

      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;
    }
}

// destructor
//

Spin1ChainWithTranslationsLong::~Spin1ChainWithTranslationsLong () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  delete[] this->LookUpTableShift;
	  int TmpMaxBitPosition = 2 * this->ChainLength;
	  for (int i = 0; i < TmpMaxBitPosition; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;
	  delete[] this->StateDescription;
	}
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1ChainWithTranslationsLong& Spin1ChainWithTranslationsLong::operator = (const Spin1ChainWithTranslationsLong& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  delete[] this->StateDescription;
	  delete[] this->LookUpTable;
	  for (int i = 1; i <= this->ChainLength; ++i)
	    {
	      delete[] this->RescalingFactors[i];
	    } 
	  delete[] this->RescalingFactors;
	  delete[] this->NbrStateInOrbit;
	}
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->MaxXMomentum = chain.MaxXMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;

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
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;

      this->MaxXMomentum = 0;
      this->StateXShift = 0;
      this->ComplementaryStateXShift = 0;
      this->XMomentumMask = ((ULONGLONG) 0x0ul);

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

AbstractHilbertSpace* Spin1ChainWithTranslationsLong::Clone()
{
  return new Spin1ChainWithTranslationsLong (*this);
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin1ChainWithTranslationsLong::TotalSz (int index)
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

// return value of the value of the sum of the square of spin projection on (Oz) 
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

double Spin1ChainWithTranslationsLong::TotalSzSz (int index)
{  
  ULONGLONG State = this->StateDescription[index];
  double TmpSzSz = 0.0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      if ((State & ((ULONGLONG) 0x3ul)) != ((ULONGLONG) 0x2ul))
	TmpSzSz += 1.0;
      State >>= 2;
    }
  return TmpSzSz;
}

// return value of twice spin projection on (Oz) for a given state
//
// stateDescription = state to which the spin projection has to be evaluated
// return value = twice spin projection on (Oz)

inline int Spin1ChainWithTranslationsLong::GetTotalSz (ULONGLONG stateDescription)
{
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      switch (stateDescription & ((ULONGLONG) 0x3ul))
	{
	case ((ULONGLONG) 0x3ul):
	  TmpSz += 2;
	  break;
	case ((ULONGLONG) 0x0ul):
	  TmpSz -= 2;
	  break;
	}
      stateDescription >>= 2;
    }
  return TmpSz;
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double Spin1ChainWithTranslationsLong::SziSzj (int i, int j, int state)
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
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslationsLong::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
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
	return this->HilbertSpaceDimension;
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient *= M_SQRT2;
	State |= (((ULONGLONG) 0x1ul) << j);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
     }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	coefficient *= M_SQRT2;
	State |= (((ULONGLONG) 0x2ul) << j);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i operator on a given state (only valid if there is no constraint on total Sz)
//
// i = operator position
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslationsLong::Spi (int i, int state, double& coefficient, int& nbrTranslation)
{
  ULONGLONG tmpState = this->StateDescription[state];
  ULONGLONG State = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= ((ULONGLONG) 0x3ul);
  switch (tmpState)
    {
    case ((ULONGLONG) 0x3ul):
      return this->HilbertSpaceDimension;
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient = M_SQRT2;
	State |= (((ULONGLONG) 0x1ul) << i);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	coefficient = M_SQRT2;
	State |= (((ULONGLONG) 0x2ul) << i);
      }
      break;
    }	  
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}

// return index of resulting state from application of S-_i operator on a given state (only valid if there is no constraint on total Sz)
//
// i = operator position
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslationsLong::Smi (int i, int state, double& coefficient, int& nbrTranslation)
{
  ULONGLONG tmpState = this->StateDescription[state];
  ULONGLONG State = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= ((ULONGLONG) 0x3ul);
  switch (tmpState)
    {
    case 0x03ul:
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
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}
    
// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslationsLong::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
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
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient *= M_SQRT2;
	State |= (((ULONGLONG) 0x1ul) << j);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
    case ((ULONGLONG) 0x0ul):
      {
	coefficient *= M_SQRT2;
	State |= (((ULONGLONG) 0x2ul) << j);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslationsLong::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation)
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
	State &= ~(((ULONGLONG) 0x1ul) << j);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient *= M_SQRT2;
	State &= ~(((ULONGLONG) 0x2ul) << j);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
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

// return index of resulting state from application of S+_i S+_i operator on a given state
//
// i = position of first S+ operator
// state = index of the state to be applied on S+_i S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslationsLong::SpiSpi (int i, int state, double& coefficient, int& nbrTranslation)
{
  ULONGLONG tmpState = this->StateDescription[state];
  if (((tmpState >> (i << 1)) & ((ULONGLONG) 0x3ul)) !=  ((ULONGLONG) 0x0ul))
    {
      return this->HilbertSpaceDimension;
    }
  coefficient = 2.0;
  tmpState |= (((ULONGLONG) 0x3ul) << (i << 1));
  return this->SymmetrizeResult(tmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}

// return index of resulting state from application of S-_i S-_i operator on a given state
//
// i = position of the S- operator
// state = index of the state to be applied on S-_i S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslationsLong::SmiSmi (int i, int state, double& coefficient, int& nbrTranslation)
{
  ULONGLONG tmpState = this->StateDescription[state];
  if (((tmpState >> (i << 1)) & ((ULONGLONG) 0x3ul)) !=  ((ULONGLONG) 0x3ul))
    {
      return this->HilbertSpaceDimension;
    }
  coefficient = 2.0;
  tmpState &= ~(((ULONGLONG) 0x3ul) << (i << 1));
  return this->SymmetrizeResult(tmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslationsLong::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  ULONGLONG tmpState = this->StateDescription[state];
  ULONGLONG State = tmpState;
  ULONGLONG tmpState2 = (tmpState >> (j << 1)) & ((ULONGLONG) 0x3ul);
  if (tmpState2 == ((ULONGLONG) 0x2ul))
    {
      coefficient = 0.0;
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
	State |= (((ULONGLONG) 0x1ul) << i);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    case ((ULONGLONG) 0x0ul):
      {
	coefficient *= M_SQRT2;
	State |= (((ULONGLONG) 0x2ul) << i);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
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
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslationsLong::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
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
	State &= ~(((ULONGLONG) 0x1ul) << i);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient *= M_SQRT2;
	State &= ~(((ULONGLONG) 0x2ul) << i);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
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

// return index of resulting state from application of S-_i1 S+_j1 S-_i2 S+_j2 operator on a given state
//
// i1 = position of leftmost S- operator
// j1 = position of leftmost S+ operator
// i2 = position of rightmost S- operator
// j2 = position of rightmost S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state (orbit index)

int Spin1ChainWithTranslationsLong::SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient, int& nbrTranslation)
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
	return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    case ((ULONGLONG) 0x2ul):
      {
	coefficient *= M_SQRT2;
	TmpState &= ~(((ULONGLONG) 0x2ul) << i1);
	return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
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
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state (orbit index)

int Spin1ChainWithTranslationsLong::SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient, int& nbrTranslation)
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
	return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
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
	return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
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


// find state index
//
// stateDescription = state description
// return value = corresponding index

int Spin1ChainWithTranslationsLong::FindStateIndex(ULONGLONG stateDescription)
{
  int TmpMaxBitPosition = (this->ChainLength << 1) - 1;
  while ((TmpMaxBitPosition >= 0) && (((stateDescription >> TmpMaxBitPosition) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul)))
    {
      --TmpMaxBitPosition;
    }
  return this->FindStateIndex(stateDescription, TmpMaxBitPosition);
}

// find state index
//
// stateDescription = state description
// maxBitPosition = maximum bit set to one in stateDescription
// return value = corresponding index

int Spin1ChainWithTranslationsLong::FindStateIndex(ULONGLONG stateDescription, int maxBitPosition)
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

ostream& Spin1ChainWithTranslationsLong::PrintState (ostream& Str, int state)
{
  ULONGLONG tmp;
  ULONGLONG TmpStateDescription = this->StateDescription[state];  
  for (int j = 0; j < this->ChainLength; j++)
    {
      tmp = TmpStateDescription & ((ULONGLONG) 0x3ul);
      if (tmp == 0)
	Str << "-1 ";
      else
	if (tmp == ((ULONGLONG) 0x2ul))
	  Str << "0 ";
	else
	  Str << "1 ";
      TmpStateDescription >>= 2;
    }
  Str << " (" << this->NbrStateInOrbit[state] << ")";
  return Str;
}


// evaluate Hilbert space dimension with no constraint on the total Sz
//
// nbrSites = number of sites
// return value = Hilbert space dimension

long Spin1ChainWithTranslationsLong::EvaluateHilbertSpaceDimension(int nbrSites)
{
  long Tmp = 3l;
  for (int i = 1; i < this->ChainLength; ++i)
    {
      Tmp *= 3l;
    }
  return Tmp;
}


// evaluate Hilbert space dimension
//
// sz = twice the Sz value
// nbrSites = number of sites
// return value = Hilbert space dimension

long Spin1ChainWithTranslationsLong::EvaluateHilbertSpaceDimension(int sz, int nbrSites)
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

// generate all states with no constraint on total Sz and no discrete symmtry constraint
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// return value = number of generated states

long Spin1ChainWithTranslationsLong::RawGenerateStates(long statePosition, int sitePosition) 
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

// generate all states corresponding to a given total Sz and no discrete symmtry constraint
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentSz = total Sz value of current state
// return value = number of generated states

long Spin1ChainWithTranslationsLong::RawGenerateStates(long statePosition, int sitePosition, int currentSz) 
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


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long Spin1ChainWithTranslationsLong::GenerateStates()
{
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  ULONGLONG Discard = ~((ULONGLONG) 0x0ul);
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX) == this->StateDescription[i]))
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
      ULONGLONG* TmpStateDescription = new ULONGLONG [TmpLargeHilbertSpaceDimension];  
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

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void Spin1ChainWithTranslationsLong::GenerateLookUpTable(unsigned long memory)
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
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentMaxMomentum = TmpMaxBitPosition - 1;
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
  this->ComputeRescalingFactors();
}

// compute the rescaling factors
//

void Spin1ChainWithTranslationsLong::ComputeRescalingFactors()
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

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix Spin1ChainWithTranslationsLong::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
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
  int TmpNbrTranslation;
  int TmpNbrTranslationToIdentity;
  Complex* TmpPhases = new Complex [2 * this->ChainLength];
  double Coef = 2.0 * M_PI * ((double) this->Momentum) / ((double) this->ChainLength);
  for (int i = 0; i < (2 * this->ChainLength); ++i)
    {
      TmpPhases[i] = Phase(Coef * ((double) i));
    }

  ULONGLONG Mask1 = (((ULONGLONG) 0x1ul) << Shift) - ((ULONGLONG) 0x1ul);
  ULONGLONG Mask2 = (((ULONGLONG) 0x1ul) << (2 * this->ChainLength)) - ((ULONGLONG) 0x1ul);
  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      ULONGLONG TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  ULONGLONG TmpState2 = (TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask1)) & Mask2;
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

// convert a state defined in the real space basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector Spin1ChainWithTranslationsLong::ConvertToKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  Spin1ChainLong* TmpSpace = (Spin1ChainLong*) space;
  ComplexVector TmpVector (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      ULONGLONG TmpState = this->StateDescription[i];
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

ComplexVector Spin1ChainWithTranslationsLong::ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  Spin1ChainLong* TmpSpace = (Spin1ChainLong*) space;
  ComplexVector TmpVector (TmpSpace->HilbertSpaceDimension, true);
  Complex* FourrierCoefficients = new Complex [this->MaxXMomentum];  
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      FourrierCoefficients[i] = Phase (-2.0 * M_PI * ((double) (i * this->Momentum) / ((double) this->MaxXMomentum)));
    }
   
  int NbrTranslations;
  double TmpCoefficient;
  for (int i = 0; i < TmpSpace->HilbertSpaceDimension; ++i)
    {
      ULONGLONG TmpState = TmpSpace->StateDescription[i];
      TmpCoefficient = 1.0;
      int Pos = this->SymmetrizeResult(TmpState, 1, TmpCoefficient, NbrTranslations);
      if (Pos < this->HilbertSpaceDimension)
	{
	  TmpVector[i] =  (state[Pos] * TmpCoefficient * FourrierCoefficients[NbrTranslations]);
	}
    }
  delete[] FourrierCoefficients;
  return TmpVector;
}
  
