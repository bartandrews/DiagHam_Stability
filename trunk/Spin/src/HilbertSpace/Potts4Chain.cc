////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of Potts-4 chain                           //
//                                                                            //
//                        last modification : 05/10/2014                      //
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


#include "HilbertSpace/Potts4Chain.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/HermitianMatrix.h"
#include <iostream>


using std::cout;
using std::endl;


#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif


// default constructor
//

Potts4Chain::Potts4Chain () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableSize = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->ChainLength = 0;
  this->Sz = 0;
  this->FixedQuantumNumberFlag = false;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// memorySize = memory size in bytes allowed for look-up table

Potts4Chain::Potts4Chain (int chainLength, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedQuantumNumberFlag = false;
  this->LargeHilbertSpaceDimension = 4l;
  for (int i = 1; i < chainLength; i++)
    this->LargeHilbertSpaceDimension *= 4l;
  
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = this->GenerateStates (this->ChainLength - 1, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

  this->LookUpTableSize = this->ChainLength;
  if (this->LookUpTableSize > 16)
    {
      this->LookUpTableSize = 16;
    }
  this->LookUpTableShift = (2 * this->ChainLength) - this->LookUpTableSize;
  this->LookUpTableSize = 1 << this->LookUpTableSize;
  this->LookUpTable = new int [this->LookUpTableSize + 1];
  for (int i = 0; i <= this->LookUpTableSize; ++i)
    {
      this->LookUpTable[i] = this->HilbertSpaceDimension;
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long Tmp = this->StateDescription[i] >> this->LookUpTableShift;
      if (this->LookUpTable[Tmp] == this->HilbertSpaceDimension)
	this->LookUpTable[Tmp] = i;
    }
}

// constructor for complete Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Potts4Chain::Potts4Chain (int chainLength, int sz, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz & 2;
  this->FixedQuantumNumberFlag = true;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->ChainLength - 1, this->Sz);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = this->GenerateStates (this->ChainLength - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

  this->LookUpTableSize = this->ChainLength;
  if (this->LookUpTableSize > 16)
    {
      this->LookUpTableSize = 16;
    }
  this->LookUpTableShift = (2 * this->ChainLength) - this->LookUpTableSize;
  this->LookUpTableSize = 1 << this->LookUpTableSize;
  this->LookUpTable = new int [this->LookUpTableSize + 1];
  for (int i = 0; i <= this->LookUpTableSize; ++i)
    {
      this->LookUpTable[i] = this->HilbertSpaceDimension;
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long Tmp = this->StateDescription[i] >> this->LookUpTableShift;
      if (this->LookUpTable[Tmp] == this->HilbertSpaceDimension)
	this->LookUpTable[Tmp] = i;
    }
}

// constructor from pre-constructed datas
//
// largehilbertSpaceDimension = Hilbert space dimension
// chainDescription = array describing states
// chainLength = number of spin 1
// sz = twice the value of total Sz component
// fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
// lookUpTable = look-up table
// lookUpTableSize = look-up table size
// lookUpTableShift = shift to apply to a state to get the key in th look-up table

Potts4Chain::Potts4Chain (long largeHilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
			  int sz, bool fixedQuantumNumberFlag, int* lookUpTable, int lookUpTableSize, 
			  int lookUpTableShift)
{
  this->Flag.Initialize();
  this->LookUpTable = lookUpTable;
  this->LookUpTableShift = lookUpTableShift;
  this->LookUpTableSize = lookUpTableSize;
  this->StateDescription = chainDescription;
  this->Sz = sz;
  this->FixedQuantumNumberFlag = fixedQuantumNumberFlag;
  this->ChainLength = chainLength;
  this->LargeHilbertSpaceDimension = largeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
}
  
// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Potts4Chain::Potts4Chain (const Potts4Chain& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableShift = 0;
      this->HilbertSpaceDimension = 0;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
      this->LargeHilbertSpaceDimension = 0l;
    }
}

// destructor
//

Potts4Chain::~Potts4Chain () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Potts4Chain& Potts4Chain::operator = (const Potts4Chain& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTable;
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->StateDescription = chain.StateDescription;
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->Sz = chain.Sz;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
   }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableShift = 0;
      this->HilbertSpaceDimension = 0;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
      this->LargeHilbertSpaceDimension = 0l;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Potts4Chain::Clone()
{
  return new Potts4Chain (*this);
}

// evaluate Hilbert space dimension
//
// currentSite = current site to occupy
// currentSzValue = state current Sz value 
// return value = Hilbert space dimension

long Potts4Chain::EvaluateHilbertSpaceDimension(int currentSite, int currentSzValue)
{
  if (currentSite < 0)
    {
      if ((currentSzValue & 2) == this->Sz)
	return 1l;
      else
	return 0l;
    }
  long TmpDimension = this->EvaluateHilbertSpaceDimension(currentSite - 1, currentSzValue + 3);
  TmpDimension += this->EvaluateHilbertSpaceDimension(currentSite - 1, currentSzValue + 2);
  TmpDimension += this->EvaluateHilbertSpaceDimension(currentSite - 1, currentSzValue + 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(currentSite - 1, currentSzValue);
  return TmpDimension;
}


// generate all states with no constraint on total Sz
//
// currentSite = current site to occupy
// currentPosition = current position of the state that has to be considered
// return value = number of generated states

long Potts4Chain::GenerateStates(int currentSite, long currentPosition) 
{
  if (currentSite < 0)
    {
      this->StateDescription[currentPosition] = 0x0l;
      return (currentPosition + 1l);
    }
  long TmpPosition = this->GenerateStates(currentSite - 1, currentPosition);
  unsigned long TmpMask = 0x3ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    this->StateDescription[currentPosition] |= TmpMask;
  TmpPosition = this->GenerateStates(currentSite - 1, currentPosition);
  TmpMask = 0x2ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    this->StateDescription[currentPosition] |= TmpMask;
  TmpPosition = this->GenerateStates(currentSite - 1, currentPosition);
  TmpMask = 0x1ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    this->StateDescription[currentPosition] |= TmpMask;
  return this->GenerateStates(currentSite - 1, currentPosition);
}

// generate all states corresponding to a given total Sz
//
// currentSite = current site to occupy
// currentSzValue = state current Sz value 
// currentPosition = current position of the state that has to be considered
// return value = number of generated states

long Potts4Chain::GenerateStates(int currentSite, int currentSzValue, long currentPosition) 
{
  if (currentSite < 0)
    {
      if ((currentSzValue & 2) == this->Sz)
	{
	  this->StateDescription[currentPosition] = 0x0l;
	  return (currentPosition + 1l);
	}
      else
	return currentPosition;
    }
  long TmpPosition = this->GenerateStates(currentSite - 1, currentSzValue + 3, currentPosition);
  unsigned long TmpMask = 0x3ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    {
      this->StateDescription[currentPosition] |= TmpMask;
    }
  TmpPosition = this->GenerateStates(currentSite - 1, currentSzValue + 2, currentPosition);
  TmpMask = 0x2ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    {
      this->StateDescription[currentPosition] |= TmpMask;
    }
  TmpPosition = this->GenerateStates(currentSite - 1, currentSzValue + 1, currentPosition);
  TmpMask = 0x1ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    {
      this->StateDescription[currentPosition] |= TmpMask;
    }
  return this->GenerateStates(currentSite - 1, currentSzValue, currentPosition);
}


// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Potts4Chain::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedQuantumNumberFlag == true)
    {
      TmpList += new SzQuantumNumber (this->Sz);
    }
  else
    {
      for (int i = 0; i < 4; i++)
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

AbstractQuantumNumber* Potts4Chain::GetQuantumNumber (int index)
{ 
  return new SzQuantumNumber (this->TotalSz(index));
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Potts4Chain::TotalSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->Sz;
  unsigned long State = this->StateDescription[index];
  unsigned long TmpSz = 0l;
  for (int i = 0; i < this->ChainLength; ++i)
    {
      TmpSz += (State & 0x3ul);
      State >>= 2;
    }
  return (((int) TmpSz) & 2);
}

// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts4Chain::SpiSpj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = tmpState;
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << j);
  switch (tmpState2)
    {
    case 0x2ul:
      tmpState |= 0x3ul << j;
      break;
    case 0x1ul:
      tmpState |= 0x2ul << j;
      break;
    case 0x0ul:
      tmpState |= 0x1ul << j;
      break;
    }	  
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x2ul:
      tmpState |= 0x3ul << i;      
      break;
    case 0x1ul:
      tmpState |= 0x2ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x1ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts4Chain::SmiSmj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = tmpState;
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << j);
  switch (tmpState2)
    {
    case 0x3ul:
      tmpState |= 0x2ul << j;
      break;
    case 0x2ul:
      tmpState |= 0x1ul << j;
      break;
    case 0x0ul:
      tmpState |= 0x3ul << j;
      break;
    }	  
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x3ul:
      tmpState |= 0x2ul << i;      
      break;
    case 0x2ul:
      tmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x3ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts4Chain::SpiSzj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = tmpState;
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  if (tmpState2 == 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = (double) tmpState2;
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x2ul:
      tmpState |= 0x3ul << i;      
      break;
    case 0x1ul:
      tmpState |= 0x2ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x1ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts4Chain::SmiSzj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = tmpState;
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  if (tmpState2 == 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = (double) tmpState2;
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x3ul:
      tmpState |= 0x2ul << i;      
      break;
    case 0x2ul:
      tmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x3ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts4Chain::Spi (int i, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x2ul:
      tmpState |= 0x3ul << i;      
      break;
    case 0x1ul:
      tmpState |= 0x2ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x1ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts4Chain::Smi (int i, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x3ul:
      tmpState |= 0x2ul << i;      
      break;
    case 0x2ul:
      tmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x3ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// translate a state assuming the system have periodic boundary conditions (increasing the site index)
//
// nbrTranslations = number of translations to apply
// state = index of the state to translate 
// return value = index of resulting state

int Potts4Chain::TranslateState (int nbrTranslations, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  TmpState = (((TmpState & (0x3ul << ((this->ChainLength - nbrTranslations) - 1ul)) << 1) << nbrTranslations)
	      | (TmpState >> ((this->ChainLength - nbrTranslations) << 1)));
  return this->FindStateIndex(TmpState);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Potts4Chain::PrintState (ostream& Str, int state)
{
  unsigned long StateDescription = this->StateDescription[state];  
  for (int j = 0; j < this->ChainLength; ++j)
    {
      Str << (StateDescription & 0x3ul) << " ";
      StateDescription >>= 2;
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

HermitianMatrix Potts4Chain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{ 
  if (nbrSites == 0)
    {
      if (szSector == 0)
	{
	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
      
    }
  if (nbrSites == this->ChainLength)
    {
      if (szSector == this->Sz)
	{
	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}      
    }
  Potts4Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  int ComplementaryQSector = this->Sz - szSector;
  if (ComplementaryQSector < 0)
    ComplementaryQSector += 4;
   if (ComplementaryQSector >= 4)
    ComplementaryQSector -= 4;
  Potts4Chain TmpHilbertSpace(this->ChainLength - nbrSites, ComplementaryQSector, 1000000);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = nbrSites * 2;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  long TmpNbrNonZeroElements = 0l;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | TmpDestinationHilbertSpace.StateDescription[j];
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
	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]);
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

ComplexMatrix Potts4Chain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, 1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
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
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}      
    }
  Potts4Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  int ComplementaryQSector = this->Sz - szSector;
  if (ComplementaryQSector < 0)
    ComplementaryQSector += 3;
  if (ComplementaryQSector >= 3)
    ComplementaryQSector -= 3;
  Potts4Chain TmpHilbertSpace(this->ChainLength - nbrSites, ComplementaryQSector, 1000000);
  ComplexMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpHilbertSpace.HilbertSpaceDimension, true);
  int Shift = nbrSites * 2;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  long TmpNbrNonZeroElements = 0l;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | TmpDestinationHilbertSpace.StateDescription[j];
	  int TmpPos = this->FindStateIndex(TmpState2);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpEntanglementMatrix.AddToMatrixElement(j, MinIndex, groundState[TmpPos]);
	    }
	}
    }
  return TmpEntanglementMatrix;
}
