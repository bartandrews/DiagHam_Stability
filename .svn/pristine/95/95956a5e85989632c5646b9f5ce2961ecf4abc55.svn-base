////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of spin 1/2 chain                          //
//                                                                            //
//                        last modification : 05/03/2001                      //
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


#include "HilbertSpace/Spin1_2Chain.h"
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


#define M_SQRT3 1.73205080756888
#define M_SQRT6 2.44948974278318


// default constructor
//

Spin1_2Chain::Spin1_2Chain ()
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
// chainLength = number of spin 1/2
// memorySize = memory size in bytes allowed for look-up table

Spin1_2Chain::Spin1_2Chain (int chainLength, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  if (this->ChainLength  > 32)
    {
      this->ChainLength = 1;
    }
  this->FixedQuantumNumberFlag = false;

  this->HilbertSpaceDimension = 1 << this->ChainLength;

  this->LookUpPosition = 0;
  this->LookUpTableSize = 1;
  memorySize >>= 2;
  this->LookUpTableMask = 0xffffffff;
  while ((this->LookUpPosition < this->ChainLength) && (memorySize >= 2))
    {
      this->LookUpTableMask <<= 1;
      this->LookUpTableSize *= 2;
      memorySize >>= 1;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateDescription[0] = 0xffffffff; 
  this->LargeHilbertSpaceDimension = this->GenerateStates (0, 0, this->StateDescription[0]);
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
}

// constructor for complete Hilbert space with a given total spin projection Sz
//
// chainLength = number of spin 1/2
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Spin1_2Chain::Spin1_2Chain (int chainLength, int sz, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  if (this->ChainLength  > 32)
    {
      this->ChainLength = 1;
    }
  this->FixedQuantumNumberFlag = true;
  this->Sz = sz;

  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->ChainLength, this->Sz);

  this->LookUpPosition = 0;
  this->LookUpTableSize = 1;
  memorySize >>= 2;
  this->LookUpTableMask = 0xffffffff;
  while ((this->LookUpPosition < this->ChainLength) && (memorySize >= 2))
    {
      this->LookUpTableMask <<= 1;
      this->LookUpTableSize *= 2;
      memorySize >>= 1;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateDescription[0] = 0xfffffffful; 
  this->LargeHilbertSpaceDimension = this->GenerateStates (0, 0, 0xfffffffful, this->ChainLength);
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
}

// constructor from pre-constructed datas
//
// hilbertSpaceDimension = Hilbert space dimension
// chainDescription = array describing states
// chainLength = number of spin 1/2
// sz = twice the value of total Sz component
// fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
// lookUpTable = look-up table
// lookUpTableSize = look-Up table size
// lookUpTablePosition = last position described by the look-Up table
// lookUpTableMask = look-Up table mask

Spin1_2Chain::Spin1_2Chain (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
			    int sz, bool fixedQuantumNumberFlag, int* lookUpTable, int lookUpTableSize, 
			    int lookUpPosition, unsigned long lookUpTableMask)
{
  this->Flag.Initialize();
  this->LookUpTable = lookUpTable;
  this->LookUpTableMask = lookUpTableMask;
  this->LookUpPosition = lookUpPosition;
  this->LookUpTableSize = lookUpTableSize;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->LargeHilbertSpaceDimension =  this->HilbertSpaceDimension;
  this->StateDescription = chainDescription;
  this->Sz = sz;
  this->FixedQuantumNumberFlag = fixedQuantumNumberFlag;
  this->ChainLength = chainLength;
}
  
// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2Chain::Spin1_2Chain (const Spin1_2Chain& chain)
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
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
    }
  else
    {
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

Spin1_2Chain::~Spin1_2Chain () 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true) && (this->ChainLength != 0) && (this->LargeHilbertSpaceDimension > 0))
    {
//      delete[] this->StateDescription;
//      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2Chain& Spin1_2Chain::operator = (const Spin1_2Chain& chain)
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
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
    }
  else
    {
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

AbstractHilbertSpace* Spin1_2Chain::Clone()
{
  return new Spin1_2Chain (*this);
}

// re-initialize chain with another total Sz component
//
// sz = twice the value of total Sz component
// return value = reference on current chain

Spin1_2Chain& Spin1_2Chain::Reinitialize(int sz)
{ 
  this->Sz = sz;
  this->HilbertSpaceDimension = this->GenerateStates (0, 0, 0xffffffff, this->ChainLength);
  return *this;
}

// generate Spin 1/2 states
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentStateDescription = description of current state
// return value = number of generated states

int Spin1_2Chain::GenerateStates(int statePosition, int sitePosition, unsigned currentStateDescription) 
{
  int NbrGeneratedState = -statePosition;
  int NextSitePosition = sitePosition + 1;
  unsigned long mask;
  sitePosition &= 0x0000001f;
  if (NextSitePosition != this->ChainLength)
    {
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription);

      mask = ~(0x00000001 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & mask);      
    }
  else
    {
      if (NextSitePosition == this->LookUpPosition)
	{
	  this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
	}
      this->StateDescription[statePosition] = currentStateDescription;
      statePosition++;
      mask = ~(0x00000001 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	{
	  this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
	}
      this->StateDescription[statePosition] = (currentStateDescription & mask);
      statePosition++;
    }
  NbrGeneratedState += statePosition;
  return NbrGeneratedState;
}

// generate Spin 1/2 states for a given total spin projection Sz
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentStateDescription = description of current state
// currentSz = total Sz value of current state
// return value = number of generated states

int Spin1_2Chain::GenerateStates(int statePosition, int sitePosition, unsigned currentStateDescription, int currentSz) 
{
  int MaxSz = (this->ChainLength - sitePosition) * 2;
  int DiffSz = currentSz - this->Sz;
  if (DiffSz > MaxSz)
    return 0;
  int NbrGeneratedState = -statePosition;
  int NextSitePosition = sitePosition + 1;
  unsigned long mask;
  sitePosition &= 0x0000001f;
  if (NextSitePosition != this->ChainLength)
    {
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription, currentSz);

      mask = ~(0x1ul << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 2);      
    }
  else
    {
      if (currentSz == this->Sz)
	{
	  if (NextSitePosition == this->LookUpPosition)
	    {
	      this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
	    }
	  this->StateDescription[statePosition] = currentStateDescription;
	  statePosition++;
	}
      else
	{
	  currentSz -= 2;
	  if (currentSz == this->Sz)
	    {
	      mask = ~(0x1ul << sitePosition);
	      if (NextSitePosition == this->LookUpPosition)
		{
		  this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
		}
	      this->StateDescription[statePosition] = (currentStateDescription & mask);
	      statePosition++;
	    }
	}
    }
  NbrGeneratedState += statePosition;
  return NbrGeneratedState;
}

// long Spin1_2Chain::GenerateStates(long statePosition, int sitePosition, int currentNbrSpinUp) 
// {
//   if (currentNbrSpinUp > sitePosition)
//     return 0l;
//   if (currentNbrSpinUp == 0)
//     {
//       this->StateDescription[statePosition] = 0x0ul;
//       return statePosition + 1;
//     }
//   unsigned long Mask = 0x1ul << sitePosition;
//   long TmpPosition = this->statePosition(statePosition, sitePosition  - 1, currentNbrSpinUp - 1);
//   for (; statePosition < TmpPosition; ++statePosition)
//     this->StateDescription[statePosition] |= Mask;
//   return this->GenerateStates(statePosition, sitePosition - 1, CurrentNbrSpinUp);
// }

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Spin1_2Chain::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedQuantumNumberFlag == true)
    {
      TmpList += new SzQuantumNumber (this->Sz);
    }
  else
    {
      int TmpSz = - this->ChainLength;
      for (int i = 0; i <= this->ChainLength; i++)
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

AbstractQuantumNumber* Spin1_2Chain::GetQuantumNumber (int index)
{ 
  return new SzQuantumNumber (this->TotalSz(index));
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin1_2Chain::TotalSz (int index)
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

// return matrix representation of Sx
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin1_2Chain::Sxi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  unsigned long tmpState;
  unsigned long State;
  unsigned long Mask = 0x1ul << i;
  unsigned long NotMask = ~Mask;
//  double Factor = M_SQRT2 * 0.5;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      tmpState = this->StateDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= 0x00000001;
      if (tmpState == 0x0ul)
	M(this->FindStateIndex(State | Mask), j) = 0.5;
      else
	M(this->FindStateIndex(State & NotMask), j) = 0.5;
    }
  return M;
}

// return matrix representation of i * Sy
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin1_2Chain::Syi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  unsigned long tmpState;
  unsigned long State;
  unsigned long Mask = 0x1ul << i;
  unsigned long NotMask = ~Mask;
//  double Factor = 0.5;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      tmpState = this->StateDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= 0x00000001;
      if (tmpState == 0x0ul)
	M(this->FindStateIndex(State | Mask), j) = 0.5;
      else
	M(this->FindStateIndex(State & NotMask), j) = -0.5;
    }
  return M;
}

// return matrix representation of Sz
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin1_2Chain::Szi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  unsigned long tmpState;
  unsigned long State;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      tmpState = this->StateDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= 0x1ul;
      if (tmpState == 0x0ul)
	M(j, j) = -0.5;
      else
	M(j, j) = 0.5;
    }
  return M;
}

// return index of resulting state from application of P_ij operator on a given state
//
// i = first position
// j = second position
// state = index of the state to be applied on P_ij operator
// return value = index of resulting state

int Spin1_2Chain::Pij (int i, int j, int state)
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

// return index of resulting state from application of a 3 sites permutation operator on a given state
//
// i = first position
// j = second position (j > i)
// k = third position  (k > j)
// state = index of the state to be applied on P_ijk operator
// return value = index of resulting state

int Spin1_2Chain::Pijk (int i, int j, int k, int state)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = tmpState & (0x1ul << i);
  unsigned long tmpState3 = tmpState & (0x1ul << j);
  unsigned long tmpState4 = tmpState & (0x1ul << k);
  tmpState &= ~((0x1ul << i) | (0x1ul << j) | (0x1ul << k));
  tmpState |= tmpState2 << (j - i);
  tmpState |= tmpState3 << (k - j);
  tmpState |= tmpState4 >> (k - i);
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of a 3 sites permutation inverse operator on a given state
//
// i = first position
// j = second position (j > i)
// k = third position  (k > j)
// state = index of the state to be applied on P_ijk operator
// return value = index of resulting state

int Spin1_2Chain::Pminusijk (int i, int j, int k, int state)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = tmpState & (0x1ul << i);
  unsigned long tmpState3 = tmpState & (0x1ul << j);
  unsigned long tmpState4 = tmpState & (0x1ul << k);
  tmpState &= ~((0x1ul << i) | (0x1ul << j) | (0x1ul << k));
  tmpState |= tmpState3 >> (j - i);
  tmpState |= tmpState4 >> (k - j);
  tmpState |= tmpState2 << (k - i);
  return this->FindStateIndex(tmpState);
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double Spin1_2Chain::SziSzj (int i, int j, int state)
{  
  unsigned long Mask = ((0x1ul << i) | (0x1ul << j));
  unsigned long tmpState = this->StateDescription[state] & Mask;
  if ((tmpState == 0x0ul) || (tmpState == Mask))
    return 0.25;
  else
    return -0.25;
}

// return the eigenvalue the product of consecutive Sz_i's 
//
// indexMin = index of the leftmost site 
// indexMax = index of the righttmost site 
// state = index of the state to consider
// return value = corresponding eigenvalue (either -1 or +1)

int Spin1_2Chain::ProdSzj (int indexMin, int indexMax, int state)
{
  int TmpSz = 0;
  unsigned long TmpState = this->StateDescription[state];
  for (int i = indexMin; i <= indexMax; i++)
    {
      TmpSz ^= (int) ((TmpState >> i) & 0x1ul);
    }
  return ((2 * TmpSz) - 1);
}


// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2Chain::Spi (int i, int state, double& coefficient)
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

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2Chain::Smi (int i, int state, double& coefficient)
{
  unsigned long State = this->StateDescription[state];
  unsigned long tmpState = (State >> i) & 0x1ul;
  if (tmpState == 0x0)
    {
      coefficient = 1.0;
      return this->FindStateIndex(State | (0x1ul << i));
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

int Spin1_2Chain::Szi (int i, int state, double& coefficient)
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

int Spin1_2Chain::SmiSpj (int i, int j, int state, double& coefficient)
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

int Spin1_2Chain::SpiSmj (int i, int j, int state, double& coefficient)
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

int Spin1_2Chain::SpiSpj (int i, int j, int state, double& coefficient)
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
	  return this->FindStateIndex(State & ~((0x1ul << j) | (0x1ul << i)));
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
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

int Spin1_2Chain::SmiSmj (int i, int j, int state, double& coefficient)
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
	  return this->FindStateIndex(State | ((0x1ul << j) | (0x1ul << i)));
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
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

int Spin1_2Chain::SpiSzj (int i, int j, int state, double& coefficient)
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
	coefficient = 0.5;
      else
	coefficient = -0.5;
      return this->FindStateIndex(State & ~(0x1ul << i));
      
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

int Spin1_2Chain::SmiSzj (int i, int j, int state, double& coefficient)
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
	coefficient = 0.5;
      else
	coefficient = -0.5;
      return this->FindStateIndex(State | (0x1ul << i));
      
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

int Spin1_2Chain::SpiSmjSzk (int i, int j, int k, int state, double& coefficient)
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

// compute the parity (prod_i Sz_i) for a given state
//
// state = index of the state to be applied on Sz_i operator
// return value = 0 if prod_i Sz_i = 1, 1 if prod_i Sz_i = -1

unsigned long Spin1_2Chain::Parity (int state)
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

int Spin1_2Chain::TranslateState (int nbrTranslations, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  TmpState = (((TmpState & ((0x1ul << (this->ChainLength - nbrTranslations)) - 1ul)) << nbrTranslations)
	      | (TmpState >> (this->ChainLength - nbrTranslations)));
  return this->FindStateIndex(TmpState);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Spin1_2Chain::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != AbstractQuantumNumber::Sz)
    return new Spin1_2Chain();
  int TmpSz = ((SzQuantumNumber&) q).GetSz();
  int HilbertSubspaceDimension = 0;
  int* TmpConvArray = new int [this->HilbertSpaceDimension];
  for (int i = 0; i < this->HilbertSpaceDimension; i++)
    {
      if (this->TotalSz(i) == TmpSz)
	{
	  TmpConvArray[HilbertSubspaceDimension] = i;
	  HilbertSubspaceDimension++;	  
	}
    }
  int* ConvArray = new int [HilbertSubspaceDimension];
  unsigned long* SubspaceDescription = new unsigned long [HilbertSubspaceDimension];
  int* SubspaceLookUpTable = new int [this->LookUpTableSize];
  unsigned long TestMask = this->StateDescription[TmpConvArray[0]] & this->LookUpTableMask;
  SubspaceLookUpTable[TestMask] = 0;
  SubspaceDescription[0] = this->StateDescription[TmpConvArray[0]];
  ConvArray[0] = TmpConvArray[0];
  for (int i = 1; i < HilbertSubspaceDimension; i++)
    {
      if ((this->StateDescription[TmpConvArray[i]] & this->LookUpTableMask) != TestMask)
	{
	  TestMask = this->StateDescription[TmpConvArray[i]] & this->LookUpTableMask;
	  SubspaceLookUpTable[TestMask] = i;
	}
      SubspaceDescription[i] = this->StateDescription[TmpConvArray[i]];
      ConvArray[i] = TmpConvArray[i];
    }
  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, HilbertSubspaceDimension, ConvArray);
  return (AbstractSpinChain*) new Spin1_2Chain (HilbertSubspaceDimension, SubspaceDescription, this->ChainLength,
						TmpSz, true, SubspaceLookUpTable, this->LookUpTableSize, 
						this->LookUpPosition, this->LookUpTableMask);
}

// find state index
//
// state = state description
// return value = corresponding index

int Spin1_2Chain::FindStateIndex(unsigned long state)
{
  int index = this->LookUpTable[state & this->LookUpTableMask];
  unsigned long* tmpState = &(this->StateDescription[index]);
  while ((index < this->HilbertSpaceDimension) && (state != *(tmpState++)))
    index++;
  return index;    
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin1_2Chain::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long Mask = 0x1ul;
  for (int k = 0; k < this->ChainLength; k++)    
    {
      if ((this->StateDescription[state] & Mask) == 0x0ul)
	Str << "+ ";
      else
	Str << "- ";
      Mask <<= 1;
    }
  return Str;
}

// evaluate Hilbert space dimension
//
// nbrSpins = number of spins
// sz = twice the z projection of the total momentum
// return value = Hilbert space dimension

int Spin1_2Chain::EvaluateHilbertSpaceDimension(int nbrSpins, int szMax)
{
   FactorialCoefficient Coef;
   Coef.SetToOne();
   Coef.FactorialMultiply(nbrSpins);
   Coef.FactorialDivide((nbrSpins + szMax) / 2);
   Coef.FactorialDivide((nbrSpins - szMax) / 2);
   return Coef.GetIntegerValue();
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition.
// 
// nbrSpinUp = number of spin up that belong to the subsytem 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix Spin1_2Chain::EvaluatePartialDensityMatrixParticlePartition (int nbrSpinUpSector, RealVector& groundState)
{
  if (nbrSpinUpSector == 0)
    {
      RealSymmetricMatrix TmpDensityMatrix(1);
      TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
      return TmpDensityMatrix;
    }

  int TmpTotalNbrSpinUp = (this->ChainLength + this->Sz) >> 1;

  if (nbrSpinUpSector == TmpTotalNbrSpinUp)
    {
      RealSymmetricMatrix TmpDensityMatrix(1);
      TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
      return TmpDensityMatrix;
    }

  int ComplementaryNbrSpinUpSector = TmpTotalNbrSpinUp - nbrSpinUpSector;
  BinomialCoefficients TmpBinomial (TmpTotalNbrSpinUp);
  double TmpInvBinomial = 1.0 / (TmpBinomial(TmpTotalNbrSpinUp, nbrSpinUpSector));

  Spin1_2Chain TmpDestinationHilbertSpace(this->ChainLength, nbrSpinUpSector, 1 << 18);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;
  Spin1_2Chain TmpHilbertSpace(this->ChainLength, ComplementaryNbrSpinUpSector, 1 << 18);
  TmpInvBinomial = sqrt(TmpInvBinomial);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex];
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace.StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
	      int TmpPos = this->FindStateIndex(TmpState | TmpState2);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  TmpStateCoefficient[Pos] = TmpInvBinomial;
		  ++Pos;
		}
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  delete[] TmpStateCoefficient;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }  
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix Spin1_2Chain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
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
  Spin1_2Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1_2Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = nbrSites;
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

RealMatrix Spin1_2Chain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
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
  Spin1_2Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1_2Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);

  RealMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int Shift = nbrSites;
  unsigned long Mask = (0x1ul << Shift) - 0x1ul;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift) & 0xfffffffful;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask);
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

ComplexMatrix Spin1_2Chain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
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
  Spin1_2Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1_2Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);

  ComplexMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int Shift = nbrSites;
  unsigned long Mask = (0x1ul << Shift) - 0x1ul;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift) & 0xfffffffful;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask);
	  int TmpPos = this->FindStateIndex(TmpState2);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	       TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos]);
	    }
	}
    }

   return TmpEntanglementMatrix;
}


// return the Bosonic Occupation of a given state in the basis
//
// index = index of the state in the basis
// finalState = reference on the array where the monomial representation has to be stored

void Spin1_2Chain::GetBosonicOccupation (unsigned int index, int * finalState)
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

unsigned long Spin1_2Chain::EncodeSiteState(int physicalState, int sitePosition)
{
  return  physicalState << sitePosition;
}

// get the normalization factor in front of each basis state (i.e. 1/sqrt(orbit size))
//
// return value = pointer to normalization factors

double* Spin1_2Chain::GetBasisNormalization()
{
  double* TmpNorm = new double[this->HilbertSpaceDimension];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      TmpNorm[i] = 1.0;
    }
  return TmpNorm;
}
 
