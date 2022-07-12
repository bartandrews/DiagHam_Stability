////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//    class of spin 1 chain with one 1/2 spin at the left end of the chain    //
//                                                                            //
//                        last modification : 14/03/2001                      //
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


#include "HilbertSpace/Spin1AKLTChain.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include <iostream.h>

#define MAX_HILBERTSPACE_DIMENSION 32768
#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif

// default constructor
//

Spin1AKLTChain::Spin1AKLTChain () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableMask = 0;
  this->LookUpPosition = 0;
  this->LookUpTableSize = 0;
  this->HilbertSpaceDimension = 0;
  this->ChainDescription = 0;
  this->ChainLength = 0;
  this->Sz = 0;
  this->FixedQuantumNumberFlag = false;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// memorySize = memory size in bytes allowed for look-up table

Spin1AKLTChain::Spin1AKLTChain (int chainLength, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedQuantumNumberFlag = false;
  this->HilbertSpaceDimension = 6;
  for (int i = 1; i < chainLength; i++)
    this->HilbertSpaceDimension *= 3;

  this->LookUpPosition = 1;
  this->LookUpTableSize = 2;
  memorySize >>= 1;
  this->LookUpTableMask = 0xfffffffe;
  while ((this->LookUpPosition < this->ChainLength) && (memorySize >=  4))
    {
      this->LookUpTableMask <<= 2;
      this->LookUpTableSize <<= 2;
      memorySize >>= 2;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];

  this->ChainDescription = new unsigned long [this->HilbertSpaceDimension];
  this->GenerateStates (0, 0, 0xffffffff);
}

// constructor for complete Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Spin1AKLTChain::Spin1AKLTChain (int chainLength, int sz, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->HilbertSpaceDimension = MAX_HILBERTSPACE_DIMENSION;
  this->Sz = sz;
  this->FixedQuantumNumberFlag = true;

  this->LookUpPosition = 1;
  this->LookUpTableSize = 2;
  memorySize >>= 1;
  this->LookUpTableMask = 0xfffffffe;
  while ((this->LookUpPosition < this->ChainLength) && (memorySize >=  4))
    {
      this->LookUpTableMask <<= 2;
      this->LookUpTableSize <<= 2;
      memorySize >>= 2;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];
  cout << this->LookUpPosition << endl;

  this->ChainDescription = new unsigned long [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates (0, 0, 0xffffffff, this->ChainLength * 2 + 1);
}

// constructor from pre-constructed datas
//
// hilbertSpaceDimension = Hilbert space dimension
// chainDescription = array describing states
// chainLength = number of spin 1
// sz = twice the value of total Sz component
// lookUpTable = look-up table
// lookUpTableSize = look-Up table size
// lookUpTablePosition = last position described by the look-Up table
// lookUpTableMask = look-Up table mask

Spin1AKLTChain::Spin1AKLTChain (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength,
				int sz, bool fixedQuantumNumber, int* lookUpTable, int lookUpTableSize, 
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
  this->FixedQuantumNumberFlag = fixedQuantumNumber;
  this->ChainLength = chainLength;
}
  

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1AKLTChain::Spin1AKLTChain (const Spin1AKLTChain& chain) 
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
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->LookUpTableSize = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
      this->ChainLength = 0;
    }
}

// destructor
//

Spin1AKLTChain::~Spin1AKLTChain () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->ChainDescription;
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1AKLTChain& Spin1AKLTChain::operator = (const Spin1AKLTChain& chain)
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
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->Sz = chain.Sz;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->LookUpTableSize = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1AKLTChain::Clone()
{
  return new Spin1AKLTChain (*this);
}

// generate all states with no constraint on total Sz
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentStateDescription = description of current state
// return value = number of generated states

int Spin1AKLTChain::GenerateStates(int statePosition,int sitePosition, unsigned long currentStateDescription) 
{
  int NextSitePosition = sitePosition + 1;     
  int NbrGeneratedState = -statePosition;
  if (sitePosition == 0)
    {
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription);

      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & 0xfffffffe) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & 0xfffffffe);
    }
  else
    {
      sitePosition--;
      sitePosition &= 0x0000000f;
      sitePosition <<= 1;
      sitePosition++;
      unsigned long mask;
      if (NextSitePosition != (this->ChainLength + 1))
	{
	  if (NextSitePosition == this->LookUpPosition)
	    this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
	  statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription);

	  mask = ~(0x00000001 << sitePosition);
	  if (NextSitePosition == this->LookUpPosition)
	    this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
	  statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & mask);
	  
	  mask = ~(0x00000003 << sitePosition);
	  if (NextSitePosition == this->LookUpPosition)
	    this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
	  statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & mask);
	}
      else
	{
	  if (NextSitePosition == this->LookUpPosition)
	    this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
	  this->ChainDescription[statePosition] = currentStateDescription;
	  statePosition++;
	  
	  mask = ~(0x00000001 << sitePosition);
	  if (NextSitePosition == this->LookUpPosition)
	    this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
	  this->ChainDescription[statePosition] = (currentStateDescription & mask);
	  statePosition++;
	  
	  mask = ~(0x00000003 << sitePosition);
	  if (NextSitePosition == this->LookUpPosition)
	    this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
	  this->ChainDescription[statePosition] = (currentStateDescription & mask);
	  statePosition++;
	}
    }
  NbrGeneratedState += statePosition;
  return NbrGeneratedState;
}

// generate all states corresponding to a given total Sz
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentStateDescription = description of current state
// currentSz = total Sz value of current state
// return value = number of generated states

int Spin1AKLTChain::GenerateStates(int statePosition,int sitePosition, unsigned long currentStateDescription, int currentSz) 
{
  int MaxSz = (this->ChainLength - sitePosition + 1) * 4;
  if (sitePosition == 0)
    MaxSz += 2;
  if ((currentSz - this->Sz) > MaxSz)
    return 0;
  int NextSitePosition = sitePosition + 1;     
  int NbrGeneratedState = -statePosition;
  if (sitePosition == 0)
    {
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription, currentSz);

      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & 0xfffffffe) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & 0xfffffffe, currentSz - 2);
    }
  else
    {
      sitePosition--;
      sitePosition &= 0x0000000f;
      sitePosition <<= 1;
      sitePosition++;
      unsigned long mask;
      if (NextSitePosition != (this->ChainLength + 1))
	{
	  if (NextSitePosition == this->LookUpPosition)
	    this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
	  statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription, currentSz);

	  mask = ~(0x00000001 << sitePosition);
	  if (NextSitePosition == this->LookUpPosition)
	    this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
	  statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 2);
	  
	  mask = ~(0x00000003 << sitePosition);
	  if (NextSitePosition == this->LookUpPosition)
	    this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
	  statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 4);
	}
      else
	{
	  if (currentSz == this->Sz)
	    {
	      if (NextSitePosition == this->LookUpPosition)
		this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
	      this->ChainDescription[statePosition] = currentStateDescription;
	      statePosition++;
	    }
	  else
	    {
	      currentSz -= 2;
	      if (currentSz == this->Sz)
		{
		  mask = ~(0x00000001 << sitePosition);
		  if (NextSitePosition == this->LookUpPosition)
		    this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
		  this->ChainDescription[statePosition] = (currentStateDescription & mask);
		  statePosition++;
		}
	      else
		{
		  currentSz -= 2;
		  if (currentSz == this->Sz)
		    {
		      mask = ~(0x00000003 << sitePosition);
		      if (NextSitePosition == this->LookUpPosition)
			this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
		      this->ChainDescription[statePosition] = (currentStateDescription & mask);
		      statePosition++;
		    }
		}
	    }
	}
    }
  NbrGeneratedState += statePosition;
  return NbrGeneratedState;
}

// return Hilbert space dimension
//
// return value = Hilbert space dimension

int Spin1AKLTChain::GetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Spin1AKLTChain::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedQuantumNumberFlag == true)
    {
      TmpList += new SzQuantumNumber (this->Sz);
    }
  else
    {
      int TmpSz = -(2 * this->ChainLength + 1);
      for (int i = 0; i <= (2 * this->ChainLength + 1); i++)
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

AbstractQuantumNumber* Spin1AKLTChain::GetQuantumNumber (int index)
{ 
  return new SzQuantumNumber (this->TotalSz(index));
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin1AKLTChain::TotalSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->Sz;
  unsigned long State = this->ChainDescription[index];
  int TmpSz = 0;
  unsigned long TmpState = State;
  State >>= 1;
  if ((TmpState & 0x00000001) != 0x00000000)
    TmpSz += 1;
  else
    TmpSz -= 1;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpState = State & 0x00000003;
      switch (TmpState)
	{
	case 0x00000003:
	  TmpSz += 2;
	  break;
	case 0x00000000:
	  TmpSz -= 2;
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

Matrix& Spin1AKLTChain::Sxi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  int pos = i;
  pos <<= 1;
  pos--;
  
  unsigned long tmpState;
  unsigned long State;
  if (i == 0)
    {
      for (int j = 0; j < this->HilbertSpaceDimension; j++)
	{
	  tmpState = this->ChainDescription[j];
	  if ((tmpState & 0x00000001) == 0x00000000)
	    {
	      M(this->FindStateIndex (tmpState | 0x00000001), j) = 1.0;
	    }
	  else
	    {
	      M(this->FindStateIndex (tmpState & 0xfffffffe), j) = 1.0;
	    }
	}
    }
  else
    {
      double Factor = M_SQRT2 * 0.5;
      for (int j = 0; j < this->HilbertSpaceDimension; j++)
	{
	  tmpState = this->ChainDescription[j];
	  State = tmpState;
	  tmpState >>= pos;
	  tmpState &= 0x00000003;
	  switch (tmpState)
	    {
	    case 0x00000003:
	      M(this->FindStateIndex((State & ~(0x00000001 << pos))), j) = Factor;
	      break;
	    case 0x00000002:
	      M(this->FindStateIndex((State | (0x00000001 << pos))), j) = Factor;
	      M(this->FindStateIndex((State & ~(0x00000002 << pos))), j) = Factor;
	      break;
	    case 0x00000000:
	      M(this->FindStateIndex((State | (0x00000002 << pos))), j) = Factor;
	      break;
	    }
	}
    }
  return M;
}

// return matrix representation of i * Sy
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin1AKLTChain::Syi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  int pos = i;
  pos <<= 1;
  pos--;
  
  unsigned long tmpState;
  unsigned long State;
  if (i == 0)
    {
      for (int j = 0; j < this->HilbertSpaceDimension; j++)
	{
	  tmpState = this->ChainDescription[j];
	  if ((tmpState & 0x00000001) == 0x00000000)
	    {
	      M(this->FindStateIndex (tmpState | 0x00000001), j) = 1.0;
	    }
	  else
	    {
	      M(this->FindStateIndex (tmpState & 0xfffffffe), j) = -1.0;
	    }
	}
    }
  else
    {
      double Factor = M_SQRT2 * 0.5;
      for (int j = 0; j < this->HilbertSpaceDimension; j++)
	{
	  tmpState = this->ChainDescription[j];
	  State = tmpState;
	  tmpState >>= pos;
	  tmpState &= 0x00000003;
	  switch (tmpState)
	    {
	    case 0x00000003:
	      M(this->FindStateIndex((State & ~(0x00000001 << pos))), j) = -Factor;
	      break;
	    case 0x00000002:
	      M(this->FindStateIndex((State | (0x00000001 << pos))), j) = Factor;
	      M(this->FindStateIndex((State & ~(0x00000002 << pos))), j) = -Factor;
	      break;
	    case 0x00000000:
	      M(this->FindStateIndex((State | (0x00000002 << pos))), j) = Factor;
	      break;
	    }
	}
    }
  return M;
}

// return matrix representation of Sz
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin1AKLTChain::Szi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  int pos = i;
  pos <<= 1;
  pos--;
  unsigned long tmpState;
  unsigned long State;
  if (i == 0)
    {
      for (int j = 0; j < this->HilbertSpaceDimension; j++)
	{
	  tmpState = this->ChainDescription[j];
	  if ((tmpState & 0x00000001) == 0x00000001)
	    {
	      M(j, j) = 0.5;
	    }
	  else
	    M(j, j) = -0.5;
	}
    }
  else
    {
      for (int j = 0; j < this->HilbertSpaceDimension; j++)
	{
	  tmpState = this->ChainDescription[j];
	  State = tmpState;
	  tmpState >>= pos;
	  tmpState &= 0x00000003;
	  switch (tmpState)
	    {
	    case 0x00000003:
	      M(j, j) = 1.0;
	      break;
	    case 0x00000000:
	      M(j, j) = -1.0;
	      break;
	    }
	}
    }
  return M;
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1AKLTChain::Spi (int i, int state, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1AKLTChain::Smi (int i, int state, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i operator on a given state
//
// i = position of Sz operator
// state = index of the state to be applied on Sz_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1AKLTChain::Szi (int i, int state, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double Spin1AKLTChain::SziSzj (int i, int j, int state)
{  
  if ((i == 0) && (j == 0))
    return 0.25;
  unsigned long tmpState = this->ChainDescription[state];
  if (i == 0)
    {
      if ((tmpState & 0x00000001) == 0x00000000)
	{
	  tmpState >>= ((j << 1) - 1);
	  tmpState &= 0x00000003;
	  if (tmpState == 0x00000002)
	    return 0.0;
	  if (tmpState == 0x00000000)
	    return 0.5;
	  else
	    return -0.5;	    
	}
      else
	{
	  tmpState >>= ((j << 1) - 1);
	  tmpState &= 0x00000003;
	  if (tmpState == 0x00000002)
	    return 0.0;
	  if (tmpState == 0x00000000)
	    return -0.5;
	  else
	    return 0.5;	    
	}
    }
  if (j == 0)
    {
      if ((tmpState & 0x00000001) == 0x00000000)
	{
	  tmpState >>= ((i << 1) - 1);
	  tmpState &= 0x00000003;
	  if (tmpState == 0x00000002)
	    return 0.0;
	  if (tmpState == 0x00000000)
	    return 0.5;
	  else
	    return -0.5;	    
	}
      else
	{
	  tmpState >>= ((i << 1) - 1);
	  tmpState &= 0x00000003;
	  if (tmpState == 0x00000002)
	    return 0.0;
	  if (tmpState == 0x00000000)
	    return -0.5;
	  else
	    return 0.5;	    
	}
    }
  unsigned long tmpState2 = (tmpState >> ((j << 1) - 1)) & 0x00000003;
  tmpState >>= ((i << 1) - 1);
  tmpState &= 0x00000003;
  if ((tmpState == 0x00000002) || (tmpState2 == 0x00000002))
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
// return value = index of resulting state

int Spin1AKLTChain::SmiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->ChainDescription[state];
  if ((i == 0) && (j == 0))
    return this->HilbertSpaceDimension;
  if (i == 0)
    {
      if ((tmpState & 0x00000001) == 0x00000000)
	return this->HilbertSpaceDimension;
      unsigned long State = tmpState;
      j <<= 1;
      j--;
      tmpState >>= j;
      tmpState &= 0x00000003;
      switch (tmpState)
	{
	case 0x00000003:
	  return this->HilbertSpaceDimension;
	case 0x00000002:
	  coefficient = M_SQRT2;
	  return this->FindStateIndex((State | (0x00000001 << j)) & 0xfffffffe);
	case 0x00000000:
	  coefficient = M_SQRT2;
	  return this->FindStateIndex((State | (0x00000002 << j)) & 0xfffffffe);
	}	  
    }
  if (j == 0)
    {
      if ((tmpState & 0x00000001) == 0x00000001)
	return this->HilbertSpaceDimension;
      unsigned long State = tmpState;
      i <<= 1;
      i--;
      tmpState >>= i;
      tmpState &= 0x00000003;
      switch (tmpState)
	{
	case 0x00000003:
	  coefficient = M_SQRT2;
	  return this->FindStateIndex((State & ~(0x00000001 << i)) | 0x00000001);
	case 0x00000002:
	  coefficient = M_SQRT2;
	  return this->FindStateIndex((State & ~(0x00000002 << i)) | 0x00000001);
	case 0x00000000:
	  return this->HilbertSpaceDimension;
	}	  
    }
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  i--;
  tmpState >>= i;
  tmpState &= 0x00000003;
  switch (tmpState)
    {
    case 0x00000003:
      coefficient = M_SQRT2;
      State &= ~(0x00000001 << i);
      break;
    case 0x00000002:
      coefficient = M_SQRT2;
      State&= ~(0x00000002 << i);
      break;
    case 0x00000000:
      return this->HilbertSpaceDimension;
    }	  
  j <<= 1;
  j--;
  tmpState2 >>= j;
  tmpState2 &= 0x00000003;
  switch (tmpState2)
    {
    case 0x00000003:
      return this->HilbertSpaceDimension;
    case 0x00000002:
      coefficient *= M_SQRT2;
      return this->FindStateIndex(State | (0x00000001 << j));
    case 0x00000000:
      coefficient *= M_SQRT2;
      return this->FindStateIndex(State | (0x00000002 << j));
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

int Spin1AKLTChain::SpiSpj (int i, int j, int state, double& coefficient)
{
  coefficient = 0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1AKLTChain::SmiSmj (int i, int j, int state, double& coefficient)
{
  coefficient = 0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1AKLTChain::SpiSzj (int i, int j, int state, double& coefficient)
{
  coefficient = 0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1AKLTChain::SmiSzj (int i, int j, int state, double& coefficient)
{
  coefficient = 0;
  return this->HilbertSpaceDimension;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Spin1AKLTChain::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != AbstractQuantumNumber::Sz)
    return new Spin1AKLTChain();
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
  return (AbstractSpinChain*) new Spin1AKLTChain (HilbertSubspaceDimension, SubspaceDescription, this->ChainLength,
						  TmpSz, true, SubspaceLookUpTable, this->LookUpTableSize, 
						  this->LookUpPosition, this->LookUpTableMask);
}

// find state index
//
// state = state description
// return value = corresponding index

int Spin1AKLTChain::FindStateIndex(unsigned long state)
{
  int index = this->LookUpTable[state & this->LookUpTableMask];
  unsigned long* tmpState = &(this->ChainDescription[index]);
  while ((index < this->HilbertSpaceDimension) && (state != *(tmpState++)))
    index++;
  return index;   
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin1AKLTChain::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long tmp;
  unsigned long StateDescription = this->ChainDescription[state];
  Str << this->FindStateIndex(StateDescription) << " : ";
  if (StateDescription & 0x00000001)
    Str << "+1/2 ";
  else
    Str << "-1/2 ";
  StateDescription >>= 1;
  for (int j = 0; j < this->ChainLength; j++)
    {
      tmp = StateDescription & 0x00000003;
      if (tmp == 0x00000000)
	Str << "-1 ";
      else
	if (tmp == 0x00000002)
	  Str << "0 ";
	else
	  Str << "1 ";
      StateDescription >>= 2;
    }
  return Str;
}
