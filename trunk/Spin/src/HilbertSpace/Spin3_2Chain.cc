////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of spin 3/2 chain                          //
//                                                                            //
//                        last modification : 08/07/2018                      //
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


#include "HilbertSpace/Spin3_2Chain.h"
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


#ifndef M_SQRT3
#define M_SQRT3 1.73205080756887729
#endif

// default constructor
//

Spin3_2Chain::Spin3_2Chain () 
{
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// memory = memory size in bytes allowed for look-up table

Spin3_2Chain::Spin3_2Chain (int chainLength, int memory) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedQuantumNumberFlag = false;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->ChainLength);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1);

  if (this->LargeHilbertSpaceDimension > 0l)
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
// memory = memory size in bytes allowed for look-up table

Spin3_2Chain::Spin3_2Chain (int chainLength, int sz, int memory) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz;
  this->FixedQuantumNumberFlag = true;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->ChainLength);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1, 0);
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
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

Spin3_2Chain::Spin3_2Chain (const Spin3_2Chain& chain) 
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
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
    }
  else
    {
      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;

      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0l;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

Spin3_2Chain::~Spin3_2Chain () 
{
  if ((this->ChainLength != 0) && (this->LargeHilbertSpaceDimension > 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin3_2Chain& Spin3_2Chain::operator = (const Spin3_2Chain& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTable;
    }  
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
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
    }
  else
    {
      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;

      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0l;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin3_2Chain::Clone()
{
  return new Spin3_2Chain (*this);
}

// evaluate Hilbert space dimension with no constraint on the total Sz
//
// nbrSites = number of sites
// return value = Hilbert space dimension

long Spin3_2Chain::EvaluateHilbertSpaceDimension(int nbrSites)
{
  long Tmp = 4l;
  for (int i = 1; i < this->ChainLength; ++i)
    {
      Tmp *= 4l;
    }
  return Tmp;
}


// evaluate Hilbert space dimension
//
// sz = twice the Sz value
// nbrSites = number of sites
// return value = Hilbert space dimension

long Spin3_2Chain::EvaluateHilbertSpaceDimension(int sz, int nbrSites)
{
  if (nbrSites == 0)
    {
      if (sz == 0)
	{
	  return 1l;	  
	}
      else
	{
	  return 0l;	  
	}
    }
  long TmpDimension = this->EvaluateHilbertSpaceDimension(sz - 3, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(sz - 1, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(sz + 1, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(sz + 3, nbrSites - 1);
  return TmpDimension;
}

// generate all states with no constraint on total Sz and no discrete symmtry constraint
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// return value = number of generated states

long Spin3_2Chain::RawGenerateStates(long statePosition, int sitePosition) 
{
  if (sitePosition < 0)
    {
      this->StateDescription[statePosition] = 0x0ul;
      return (statePosition + 1l);
    }

  unsigned long TmpMask = 0x3ul << (sitePosition * 2);
  long TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = 0x2ul << (sitePosition * 2);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = 0x1ul << (sitePosition * 2);
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

long Spin3_2Chain::RawGenerateStates(long statePosition, int sitePosition, int currentSz) 
{
  if (sitePosition < 0)
    {
      if (currentSz == this->Sz)
	{
	  this->StateDescription[statePosition] = 0x0ul;
	  return (statePosition + 1l);
	}
      else
	{
	  return statePosition;
	}
    }
  unsigned long TmpMask = 0x3ul << (sitePosition * 2);
  long TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1, currentSz + 3);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = 0x2ul << (sitePosition * 2);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1, currentSz + 1);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = 0x1ul << (sitePosition * 2);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1, currentSz - 1);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  return this->RawGenerateStates(statePosition, sitePosition - 1, currentSz - 3);  
}


// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Spin3_2Chain::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedQuantumNumberFlag == true)
    {
      TmpList += new SzQuantumNumber (this->Sz);
    }
  else
    {
      int TmpSz = - 3 * this->ChainLength;
      for (int i = 0; i <= (6 * this->ChainLength); i++)
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

AbstractQuantumNumber* Spin3_2Chain::GetQuantumNumber (int index)
{ 
  return new SzQuantumNumber (this->TotalSz(index));
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin3_2Chain::TotalSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->Sz;
  unsigned long State = this->StateDescription[index];
  int TmpSz = 0;
  unsigned long TmpState;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpState = State & 0x3ul;
      switch (TmpState)
	{
	case 0x3ul:
	  TmpSz += 3;
	  break;
	case 0x2ul:
	  TmpSz++;
	  break;
	case 0x1ul:
	  TmpSz--;
	  break;
	case 0x0ul:
	  TmpSz -= 3;
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

Matrix& Spin3_2Chain::Sxi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  i *= 2;
  cout << "error, Spin3_2Chain::Sxi is not implemented" << endl;
  return M;
}

// return matrix representation of i * Sy
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin3_2Chain::Syi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  i *= 2;
  cout << "error, Spin3_2Chain::Syi is not implemented" << endl;
  return M;
}

// return matrix representation of Sz
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin3_2Chain::Szi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  i *= 2;
  cout << "error, Spin3_2Chain::Szi is not implemented" << endl;
  return M;
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin3_2Chain::Spi (int i, int state, double& coefficient)
{
  i <<= 1;
  unsigned long State = (this->StateDescription[state] >> i);
  unsigned long tmpState = State;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	State &= ~(0x3ul << i);
	State |= (0x2ul << i);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	State |= (0x1ul << i);
      }
      break;
    }
  return this->FindStateIndex(State);
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin3_2Chain::Smi (int i, int state, double& coefficient)
{
  i <<= 1;
  unsigned long State = (this->StateDescription[state] >> i);
  unsigned long tmpState = State;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      {
	coefficient = M_SQRT3;
	State &= ~(0x2ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient = 2.0;
	State &= ~(0x2ul << i);
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient = M_SQRT3;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x0ul:
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
      break;
    }
  return this->FindStateIndex(State);
}

// return index of resulting state from application of Sz_i operator on a given state
//
// i = position of Sz operator
// state = index of the state to be applied on Sz_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin3_2Chain::Szi (int i, int state, double& coefficient)
{
  unsigned long tmpState = (this->StateDescription[state] >> (i * 2));
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      coefficient = 1.5;
      break;
    case 0x2ul:
      coefficient = 0.5;
      break;
    case 0x1ul:
      coefficient = -0.5;
      break;
     case 0x0ul:
      coefficient = -1.5;
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

double Spin3_2Chain::SziSzj (int i, int j, int state)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = (tmpState >> (j * 2)) & 0x3ul;
  double Tmp = 0.0;
  switch (tmpState2)
    {
    case 0x3ul:
      Tmp = 1.5;
      break;
    case 0x2ul:
      Tmp = 0.5;
      break;
    case 0x1ul:
      Tmp = -0.5;
      break;
    case 0x0ul:
      Tmp = -1.5;
      break;
    }
  tmpState2 = (tmpState >> (i * 2)) & 0x3ul;
  switch (tmpState2)
    {
    case 0x3ul:
      Tmp *= 1.5;
      break;
    case 0x2ul:
      Tmp *= 0.5;
      break;
    case 0x1ul:
      Tmp *= -0.5;
      break;
    case 0x0ul:
      Tmp *= -1.5;
      break;
    }
  return Tmp;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin3_2Chain::SmiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i *= 2;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      {
	coefficient = M_SQRT3;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient = 2.0;
	State &= ~(0x2ul << i);
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient = M_SQRT3;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  j *= 2;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  switch (tmpState2)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	State |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	State &= ~(0x1ul << j);
	State |= (0x2ul << j);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	State |= (0x1ul << j);
      }
      break;
    }	  
  return this->FindStateIndex(State);
}

// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin3_2Chain::SpiSpj (int i, int j, int state, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long State = TmpState;
  i *= 2;
  TmpState >>= i;
  TmpState &= 0x3ul;
  switch (TmpState)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT3;
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	State &= ~(0x1ul << i);
	State |= (0x2ul << i);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT3;
	State |= (0x1ul << i);
      }
      break;
    }	  
  TmpState = State;
  j *= 2;
  TmpState >>= j;
  TmpState &= 0x3ul;
  switch (TmpState)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	State |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	State &= ~(0x1ul << j);
	State |= (0x2ul << j);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	State |= (0x1ul << j);
      }
      break;
    }	  
  return this->FindStateIndex(State);
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin3_2Chain::SmiSmj (int i, int j, int state, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long State = TmpState;
  i *= 2;
  TmpState >>= i;
  TmpState &= 0x3ul;
  switch (TmpState)
    {
    case 0x3ul:
      {
	coefficient = M_SQRT3;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient = 2.0;
	State &= ~(0x2ul << i);
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient = M_SQRT3;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  TmpState = State;
  j *= 2;
  TmpState >>= j;
  TmpState &= 0x3ul;
  switch (TmpState)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	State &= ~(0x1ul << j);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	State &= ~(0x2ul << j);
	State |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	State &= ~(0x1ul << j);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  return this->FindStateIndex(State);
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin3_2Chain::SpiSzj (int i, int j, int state, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long State = TmpState;
  TmpState >>= (j * 2);
  TmpState &= 0x3ul;
  switch (TmpState)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }
  TmpState = State;
  i *= 2;
  TmpState >>= i;
  TmpState &= 0x3ul;
  switch (TmpState)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	State &= ~(0x1ul << i);
	State |= (0x2ul << i);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	State |= (0x1ul << i);
      }
      break;
    }	  
  return this->FindStateIndex(State);
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin3_2Chain::SmiSzj (int i, int j, int state, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long State = TmpState;
  TmpState >>= (j * 2);
  TmpState &= 0x3ul;
  switch (TmpState)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }
  TmpState = State;
  i *= 2;
  TmpState >>= i;
  TmpState &= 0x3ul;
  switch (TmpState)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	State &= ~(0x1ul << j);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	State &= ~(0x2ul << j);
	State |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	State &= ~(0x1ul << j);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  return this->FindStateIndex(State);
}

// return index of resulting state from application of S+_i S-_j Sz_k operator on a given state
//
// i = position of S+ operator
// j = position of S- operator
// k = position of Sz operator
// state = index of the state to be applied on S+_i S-_j Sz_k operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin3_2Chain::SpiSmjSzk (int i, int j, int k, int state, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[state];
  k *= 2;
  switch ((TmpState >> k) & 0x3ul)
    {
    case 0x3ul:
      coefficient = 1.5;
      break;
    case 0x2ul:
      coefficient = 0.5;
      break;
    case 0x1ul:
      coefficient = -0.5;
      break;
    case 0x0ul:
      coefficient = -1.5;
      break;
    }  
  j *= 2;
  switch ((TmpState >> j) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << j);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << j);
	TmpState |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << j);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  i *= 2;
  switch ((TmpState >> i) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x1ul << i);
	TmpState |= (0x2ul << i);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << i);
      }
      break;
    }	  
  return this->FindStateIndex(TmpState);
}


// return index of resulting state from application of S-_i S+_j Sz_k operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// k = position of Sz operator
// state = index of the state to be applied on S-_i S+_j Sz_k operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin3_2Chain::SmiSpjSzk (int i, int j, int k, int state, double& coefficient)
{
  double Tmp;
  switch ((this->StateDescription[state] >> (k << 1)) & 0x3ul)
    {
    case 0x3ul:
      Tmp = 1.5;
      break;
    case 0x2ul:
      Tmp = 0.5;
      break;
    case 0x1ul:
      Tmp = -0.5;
      break;
    case 0x0ul:
      Tmp = -1.5;
      break;
    }  
  int TmpPos = this->SpiSmj(i, j, state, coefficient);
  coefficient *= Tmp;
  return TmpPos;
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

int Spin3_2Chain::SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient)
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
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << j2);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	TmpState &= ~(0x1ul << j2);
	TmpState |= (0x2ul << j2);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << j2);
      }
      break;
    }
  i2 <<= 1;
  switch ((TmpState >> i2) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i2);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << i2);
	TmpState |= (0x1ul << i2);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i2);
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
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << j1);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x1ul << j1);
	TmpState |= (0x2ul << j1);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << j1);
      }
      break;
    }
  i1 <<= 1;
  switch ((TmpState >> i1) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i1);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << i1);
	TmpState |= (0x1ul << i1);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i1);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    }	  
  return this->FindStateIndex(TmpState);
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

int Spin3_2Chain::SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient)
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
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << j2);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	TmpState &= ~(0x1ul << j2);
	TmpState |= (0x2ul << j2);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << j2);
      }
      break;
    }
  i2 <<= 1;
  switch ((TmpState >> i2) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i2);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << i2);
	TmpState |= (0x1ul << i2);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i2);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    }	  
  switch ((TmpState >> (j1 * 2)) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }
  switch ((TmpState >> (i1 * 2)) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }
  return this->FindStateIndex(TmpState);
}

// translate a state assuming the system have periodic boundary conditions (increasing the site index)
//
// nbrTranslations = number of translations to apply
// state = index of the state to translate 
// return value = index of resulting state

int Spin3_2Chain::TranslateState (int nbrTranslations, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  TmpState = (((TmpState & (0x3ul << ((this->ChainLength - nbrTranslations) - 1ul)) * 2) << nbrTranslations)
	      | (TmpState >> ((this->ChainLength - nbrTranslations) * 2)));
  return this->FindStateIndex(TmpState);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Spin3_2Chain::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
  return 0;
}

// find state index
//
// stateDescription = state description
// maxBitPosition = maximum bit set to one in stateDescription
// return value = corresponding index

int Spin3_2Chain::FindStateIndex(unsigned long stateDescription, int maxBitPosition)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    return this->HilbertSpaceDimension;
  if (this->LookUpTableShift[maxBitPosition] < 0)
    return this->HilbertSpaceDimension;
  long PosMax = stateDescription >> this->LookUpTableShift[maxBitPosition];
  long PosMin = this->LookUpTable[maxBitPosition][PosMax];
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

ostream& Spin3_2Chain::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long tmp;
  unsigned long TmpStateDescription = this->StateDescription[state];  
  for (int j = 0; j < this->ChainLength; j++)
    {
      switch (TmpStateDescription & 0x3ul)
	{
	case 0x0ul:
	  Str << "-3 ";
	  break;
	case 0x1ul:
	  Str << "-1 ";
	  break;
	case 0x2ul:
	  Str << "+1 ";
	  break;
	case 0x3ul:
	  Str << "+3 ";
	  break;
	}
      TmpStateDescription >>= 2;
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

RealSymmetricMatrix Spin3_2Chain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
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
  Spin3_2Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin3_2Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = nbrSites * 2;
  unsigned long Mask = (0x1ul << Shift) - 0x1ul;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  long TmpNbrNonZeroElements = 0l;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift) & 0xfffffffful;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask);
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

RealMatrix Spin3_2Chain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
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
  Spin3_2Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin3_2Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);
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

ComplexMatrix Spin3_2Chain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
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
  Spin3_2Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin3_2Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);

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
	  int TmpPos = this->FindStateIndex(TmpState2);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	       TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos]);
	    }
	}
    }

   return TmpEntanglementMatrix;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void Spin3_2Chain::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * (2 * this->ChainLength));
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
  int CurrentMaxMomentum = TmpMaxBitPosition;
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

