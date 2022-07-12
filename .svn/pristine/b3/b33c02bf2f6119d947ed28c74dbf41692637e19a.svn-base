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
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include <iostream.h>


#define M_SQRT3 1.73205080756888
#define M_SQRT6 2.44948974278318
#define MAXHILBERTSPACEDIMENSION 4194304


// default constructor
//

Spin1_2Chain::Spin1_2Chain ()
{
  this->Flag.Initialize();
  this->ChainLength = 0;
  this->HilbertSpaceDimension = 0;
  this->Sz = 0;
  this->ChainDescription = 0;
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
  this->ChainDescription = new unsigned long [this->HilbertSpaceDimension];
  this->ChainDescription[0] = 0xffffffff; 
  this->HilbertSpaceDimension = this->GenerateStates (0, 0, this->ChainDescription[0]);
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

  this->HilbertSpaceDimension = MAXHILBERTSPACEDIMENSION;

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
  this->ChainDescription = new unsigned long [this->HilbertSpaceDimension];
  this->ChainDescription[0] = 0xffffffff; 
  this->HilbertSpaceDimension = this->GenerateStates (0, 0, 0xffffffff, this->ChainLength);
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
  this->ChainDescription = chainDescription;
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
      this->Sz = chain.Sz;
      this->ChainDescription = chain.ChainDescription;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      
    }
  else
    {
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->ChainDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
    }
}

// destructor
//

Spin1_2Chain::~Spin1_2Chain () 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true) && (this->ChainLength != 0))
    {
      delete[] this->ChainDescription;
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
      delete[] this->ChainDescription;
    }
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->Sz = chain.Sz;
      this->ChainDescription = chain.ChainDescription;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
    }
  else
    {
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->ChainDescription = 0;
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
      this->ChainDescription[statePosition] = currentStateDescription;
      statePosition++;
      mask = ~(0x00000001 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	{
	  this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
	}
      this->ChainDescription[statePosition] = (currentStateDescription & mask);
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

      mask = ~(0x00000001 << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 2);      
    }
  else
    {
      if (currentSz == this-> Sz)
	{
	  if (NextSitePosition == this->LookUpPosition)
	    {
	      this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
	    }
	  this->ChainDescription[statePosition] = currentStateDescription;
	  statePosition++;
	}
      else
	{
	  currentSz -= 2;
	  if (currentSz == this-> Sz)
	    {
	      mask = ~(0x00000001 << sitePosition);
	      if (NextSitePosition == this->LookUpPosition)
		{
		  this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
		}
	      this->ChainDescription[statePosition] = (currentStateDescription & mask);
	      statePosition++;
	    }
	}
    }
  NbrGeneratedState += statePosition;
  return NbrGeneratedState;
}

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
  unsigned long State = this->ChainDescription[index];
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpSz += ((State & 0x00000001) << 1);
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
  unsigned long Mask = 0x00000001 << i;
  unsigned long NotMask = ~Mask;
//  double Factor = M_SQRT2 * 0.5;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      tmpState = this->ChainDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= 0x00000001;
      if (tmpState == 0)
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
  unsigned long Mask = 0x00000001 << i;
  unsigned long NotMask = ~Mask;
//  double Factor = 0.5;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      tmpState = this->ChainDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= 0x00000001;
      if (tmpState == 0)
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
      tmpState = this->ChainDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= 0x00000001;
      if (tmpState == 0)
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
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpMask = (0x00000001 << i) | (0x00000001 << j);
  unsigned long tmpState2 = tmpState & tmpMask;
  unsigned long tmpState3 = ~tmpState & tmpMask;
  if ((tmpState2 == 0x00000000) || (tmpState3 == 0x00000000))
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

double Spin1_2Chain::SziSzj (int i, int j, int state)
{  
  unsigned long Mask = ((0x00000001 << i) | (0x00000001 << j));
  unsigned long tmpState = this->ChainDescription[state] & Mask;
  if ((tmpState == 0x00000000) || (tmpState == Mask))
    return 0.25;
  else
    return -0.25;
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2Chain::Spi (int i, int state, double& coefficient)
{
  unsigned long State = this->ChainDescription[state];
  unsigned long tmpState = (State >> i) & 0x00000001;
  if (tmpState == 0x00000001)
    {
      coefficient = 1.0;
      return this->FindStateIndex(State & ~(0x00000001 << i));
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
  unsigned long State = this->ChainDescription[state];
  unsigned long tmpState = (State >> i) & 0x00000001;
  if (tmpState == 0x00000000)
    {
      coefficient = 1.0;
      return this->FindStateIndex(State | (0x00000001 << i));
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
  unsigned long tmpState = (this->ChainDescription[state] >> i) & 0x00000001;
  if (tmpState == 0x00000000)
    coefficient = 0.5;
  else
    coefficient = -0.5;
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
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x00000001;
  if (i != j)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x00000001;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x00000001)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex((State | (0x00000001 << j)) & ~(0x00000001 << i));
	}
      else
	{
	  coefficient = 0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0)
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
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x00000001;
  if (i != j)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x00000001;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x00000003)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex(State & ~((0x00000001 << j) | (0x00000001 << i)));
	}
      else
	{
	  coefficient = 0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0)
    {
      coefficient = 0;
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
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x00000001;
  if (i != j)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x00000001;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x00000000)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex(State | ((0x00000001 << j) | (0x00000001 << i)));
	}
      else
	{
	  coefficient = 0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0)
    {
      coefficient = 0;
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
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x00000001;
  if (tmpState == 0x00000001)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x00000001;
      if (tmpState2 == 0x00000000)
	coefficient = 0.5;
      else
	coefficient = -0.5;
      return this->FindStateIndex(State & ~(0x00000001 << i));
      
    }
  else
    {
      coefficient = 0;
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
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x00000001;
  if (tmpState == 0x00000000)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x00000001;
      if (tmpState2 == 0x00000000)
	coefficient = 0.5;
      else
	coefficient = -0.5;
      return this->FindStateIndex(State | (0x00000001 << i));
      
    }
  else
    {
      coefficient = 0;
      return this->HilbertSpaceDimension;
    }
  return this->HilbertSpaceDimension;
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
  unsigned long TestMask = this->ChainDescription[TmpConvArray[0]] & this->LookUpTableMask;
  SubspaceLookUpTable[TestMask] = 0;
  SubspaceDescription[0] = this->ChainDescription[TmpConvArray[0]];
  ConvArray[0] = TmpConvArray[0];
  for (int i = 1; i < HilbertSubspaceDimension; i++)
    {
      if ((this->ChainDescription[TmpConvArray[i]] & this->LookUpTableMask) != TestMask)
	{
	  TestMask = this->ChainDescription[TmpConvArray[i]] & this->LookUpTableMask;
	  SubspaceLookUpTable[TestMask] = i;
	}
      SubspaceDescription[i] = this->ChainDescription[TmpConvArray[i]];
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
//  if (this->FixedQuantumNumberFlag == false)
//    return ~state;
  int index = this->LookUpTable[state & this->LookUpTableMask];
  unsigned long* tmpState = &(this->ChainDescription[index]);
  while ((index < this->HilbertSpaceDimension) && (state != *(tmpState++)))
    index++;
//  if (index != (~state))
//    cout << "error" << ~state << " " << index << endl;
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
  unsigned long Mask = 0x00000001;
  for (int k = 0; k < this->ChainLength; k++)    
    {
      if ((this->ChainDescription[state] & Mask) == 0x00000000)
	Str << "+ ";
      else
	Str << "- ";
      Mask <<= 1;
    }
  return Str;
}
