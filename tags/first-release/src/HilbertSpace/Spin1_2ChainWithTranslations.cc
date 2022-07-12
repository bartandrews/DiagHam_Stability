////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin 1/2 chain with translastion invariance          //
//                                                                            //
//                        last modification : 29/01/2002                      //
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


#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/MomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include <iostream.h>


#define M_SQRT3 1.73205080756888
#define M_SQRT6 2.44948974278318
#define MAXHILBERTSPACEDIMENSION 4194304


// default constructor
//

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations ()
{
  this->Flag.Initialize();
  this->ChainLength = 0;
  this->ReducedChainLength = 0;
  this->HilbertSpaceDimension = 0;
  this->Sz = 0;
  this->ChainDescription = 0;  
  this->OrbitSize = 0;
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->FixedSpinProjectionFlag = false;
  this->FixedMomentumFlag = false;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1/2
// memorySize = memory size in bytes allowed for look-up table

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (int chainLength, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  if (this->ChainLength  > 32)
    {
      this->ChainLength = 1;
    }
  this->ReducedChainLength = this->ChainLength - 1;

  this->FixedSpinProjectionFlag = false;
  this->FixedMomentumFlag = false;

  this->HilbertSpaceDimension = 1 << this->ChainLength;

  this->ChainDescription = new unsigned long [this->HilbertSpaceDimension];
  this->OrbitSize = new int [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = 0;
  this->GenerateStates (0, (0xffffffff >> (32 - this->ChainLength)));
  this->BuildLookUpTable(memorySize);
}

// constructor for complete Hilbert space with a given total spin projection Sz
//
// chainLength = number of spin 1/2
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (int chainLength, int sz, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  if (this->ChainLength  > 32)
    {
      this->ChainLength = 1;
    }
  this->ReducedChainLength = this->ChainLength - 1;
  this->FixedSpinProjectionFlag = true;
  this->Sz = sz;

  this->HilbertSpaceDimension = MAXHILBERTSPACEDIMENSION;

  this->ChainDescription = new unsigned long [this->HilbertSpaceDimension];
  this->OrbitSize = new int [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = 0;
  this->GenerateStates (0, (0xffffffff >> (32 - this->ChainLength)), this->ChainLength);
  this->BuildLookUpTable(memorySize);
}

// constructor for complete Hilbert space with a given total spin projection Sz an a given momentum 
//
// chainLength = number of spin 1/2
// sz = twice the value of total Sz component
// momentum = momentum
// memorySize = memory size in bytes allowed for look-up table

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (int chainLength, int sz, int momentum, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  if (this->ChainLength  > 32)
    {
      this->ChainLength = 1;
    }
  this->ReducedChainLength = this->ChainLength - 1;

  this->FixedSpinProjectionFlag = true;
  this->Sz = sz;
  this->FixedMomentumFlag = true;
  this->Momentum = momentum;

  this->HilbertSpaceDimension = MAXHILBERTSPACEDIMENSION;

  this->ChainDescription = new unsigned long [this->HilbertSpaceDimension];
  this->OrbitSize = new int [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = 0;
  this->GenerateStates (0, (0xffffffff >> (32 - this->ChainLength)), this->ChainLength);
  this->BuildLookUpTable(memorySize);
}

// constructor from pre-constructed datas
//
// hilbertSpaceDimension = Hilbert space dimension
// chainDescription = array describing states
// orbitSize = number of state per orbit
// chainLength = number of spin 1/2
// sz = twice the value of total Sz component
// fixedSpinPorjectionFlag = true if hilbert space is restricted to a given spin projection
// momentum = momentum
// fixedSpinPorjectionFlag = true if hilbert space is restricted to a given momentum
// lookUpTableSize = look-Up table size

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (int hilbertSpaceDimension, unsigned long* chainDescription, 
							    int* orbitSize, int chainLength, 
							    int sz, bool fixedSpinProjectionFlag, 
							    int momentum, bool fixedMomentumFlag,
							    int lookUpTableSize)
{
  this->Flag.Initialize();
  this->LookUpTableSize = lookUpTableSize;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->ChainDescription = chainDescription;
  this->OrbitSize = orbitSize;
  this->Sz = sz;
  this->Momentum = momentum;
  this->FixedSpinProjectionFlag = fixedSpinProjectionFlag;
  this->FixedMomentumFlag = fixedMomentumFlag;
  this->ChainLength = chainLength;
  this->ReducedChainLength = this->ChainLength - 1;
  this->BuildLookUpTable(this->LookUpTableSize);
}
  
// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (const Spin1_2ChainWithTranslations& chain)
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->ReducedChainLength = chain.ReducedChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->Sz = chain.Sz;
      this->ChainDescription = chain.ChainDescription;
      this->OrbitSize = chain.OrbitSize;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      
    }
  else
    {
      this->ChainLength = 0;
      this->ReducedChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->ChainDescription = 0;
      this->OrbitSize = 0;
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
    }
}

// destructor
//

Spin1_2ChainWithTranslations::~Spin1_2ChainWithTranslations () 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true) && (this->ChainLength != 0))
    {
      delete[] this->ChainDescription;
      delete[] this->OrbitSize;
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainWithTranslations& Spin1_2ChainWithTranslations::operator = (const Spin1_2ChainWithTranslations& chain)
{
  if ((this->Flag.Used() == true) && (this->ChainLength != 0) && (this->Flag.Shared() == false))
    {
      delete[] this->ChainDescription;
      delete[] this->OrbitSize;
      delete[] this->LookUpTable;
    }
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->ReducedChainLength = chain.ReducedChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->Sz = chain.Sz;
      this->ChainDescription = chain.ChainDescription;
      this->OrbitSize = chain.OrbitSize;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
    }
  else
    {
      this->ChainLength = 0;
      this->ReducedChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->ChainDescription = 0;
      this->OrbitSize = 0;
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainWithTranslations::Clone()
{
  return new Spin1_2ChainWithTranslations (*this);
}

// re-initialize chain with another total Sz component
//
// sz = twice the value of total Sz component
// return value = reference on current chain

Spin1_2ChainWithTranslations& Spin1_2ChainWithTranslations::Reinitialize(int sz)
{ 
  this->Sz = sz;
  this->HilbertSpaceDimension = 0;
  this->GenerateStates (0, (0xffffffff >> (32 - this->ChainLength)), this->ChainLength);
  this->BuildLookUpTable();
  return *this;
}

// generate Spin 1/2 states
//
// sitePosition = site on chain where spin has to be changed
// currentStateDescription = description of current state
// return value = number of generated states

void Spin1_2ChainWithTranslations::GenerateStates(int sitePosition, unsigned long currentStateDescription) 
{
  int NextSitePosition = sitePosition + 1;
  unsigned long mask;
  sitePosition &= 0x0000001f;
  if (NextSitePosition != this->ChainLength)
    {
      this->GenerateStates(NextSitePosition, currentStateDescription);
      mask = ~(0x00000001 << sitePosition);
      this->GenerateStates(NextSitePosition, currentStateDescription & mask);      
    }
  else
    {
      this->StoreCanonicalForm(currentStateDescription);
      mask = ~(0x00000001 << sitePosition);
      this->StoreCanonicalForm(currentStateDescription & mask);
    }
}

// generate Spin 1/2 states for a given total spin projection Sz
//
// sitePosition = site on chain where spin has to be changed
// currentStateDescription = description of current state
// currentSz = total Sz value of current state
// return value = number of generated states

void Spin1_2ChainWithTranslations::GenerateStates(int sitePosition, unsigned long currentStateDescription, int currentSz) 
{
  int MaxSz = (this->ChainLength - sitePosition) * 2;
  int DiffSz = currentSz - this->Sz;
  if (DiffSz > MaxSz)
    return;
  int NextSitePosition = sitePosition + 1;
  unsigned long mask;
  sitePosition &= 0x0000001f;
  if (NextSitePosition != this->ChainLength)
    {
      this->GenerateStates(NextSitePosition, currentStateDescription, currentSz);
      mask = ~(0x00000001 << sitePosition);
      this->GenerateStates(NextSitePosition, currentStateDescription & mask, currentSz - 2);      
    }
  else
    {
      if (currentSz == this-> Sz)
	{
	  this->StoreCanonicalForm(currentStateDescription);
	}
      else
	{
	  currentSz -= 2;
	  if (currentSz == this-> Sz)
	    {
	      mask = ~(0x00000001 << sitePosition);
	      this->StoreCanonicalForm(currentStateDescription & mask);
	    }
	}
    }
}

// store a state in its canonical form
//
// stateDescription = description of the stateto store

void Spin1_2ChainWithTranslations::StoreCanonicalForm (unsigned long stateDescription)
{
  unsigned long CanonicalState = stateDescription;
  int index = 1;
  while (index < this->ChainLength)
    {
      stateDescription = (stateDescription >> 1) | ((stateDescription & 0x1) << this->ReducedChainLength);
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	}
      ++index;
    }
  index = 0;
  while ((index < this->HilbertSpaceDimension) && (this->ChainDescription[index] < CanonicalState))
    {
      ++index;
    }
  if (index == this->HilbertSpaceDimension)
    {
      this->ChainDescription[index] = CanonicalState;
      this->OrbitSize[index] = 1;
      this->HilbertSpaceDimension++;
    }
  else
    {
      if (this->ChainDescription[index] != CanonicalState)
	{
	  for (int i = this->HilbertSpaceDimension; i > index; --i)
	    {
	      this->ChainDescription[i] = this->ChainDescription[i - 1];
	      this->OrbitSize[i] = this->OrbitSize[i - 1];
	    }
	  this->ChainDescription[index] = CanonicalState;
	  this->OrbitSize[index] = 1;
	  ++this->HilbertSpaceDimension;	  
	}
      else
	{
	  ++this->OrbitSize[index];
	}
    }
}

// build look-up table
//
// memorySize = memory size in bytes allowed for look-up table

void Spin1_2ChainWithTranslations::BuildLookUpTable(int memorySize)
{
  if (memorySize != 0)
    {
      this->LookUpTableSize = 2;
      memorySize >>= 2;
      this->LookUpTableShift = this->ChainLength - 1;
      while ((this->LookUpTableShift > 0) && (memorySize >= 2))
	{
	  memorySize >>= 1;
	  --this->LookUpTableShift;
	  this->LookUpTableSize <<= 1;
	}
      this->LookUpTable = new int [this->LookUpTableSize];
    }

  unsigned long TmpLookUp = this->ChainDescription[0] >> this->LookUpTableShift;
  unsigned long TmpOldLookUp = TmpLookUp;
  this->LookUpTable[TmpLookUp] = 0;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      TmpLookUp = this->ChainDescription[i] >> this->LookUpTableShift;
      if (TmpLookUp != TmpOldLookUp)
	{
	  this->LookUpTable[TmpLookUp] = i;
	  TmpOldLookUp = TmpLookUp;
	}
    }
}

// Apply momentum criteria to keep allowed state 
//

void Spin1_2ChainWithTranslations::ApplyMomentumCriteria ()
{
  int TmpHilbertSpaceDimension = 0;
  unsigned long* TmpChainDescription = new unsigned long [this->HilbertSpaceDimension];
  int* TmpOrbitSize = new int [this->HilbertSpaceDimension];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if ((this->Momentum % (this->ChainLength / this->OrbitSize[i])) == 0)
	{
	  TmpChainDescription[TmpHilbertSpaceDimension] = this->ChainDescription[i];
	  TmpOrbitSize[TmpHilbertSpaceDimension] = this->OrbitSize[i];
	  TmpHilbertSpaceDimension++;
	}
    }
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true) && (this->ChainLength != 0))
    {
      delete[] this->ChainDescription;
      delete[] this->OrbitSize;
    }
  this->ChainDescription = TmpChainDescription;
  this->OrbitSize = TmpOrbitSize;
  this->HilbertSpaceDimension = TmpHilbertSpaceDimension;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Spin1_2ChainWithTranslations::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedSpinProjectionFlag == true)
    if (this->FixedMomentumFlag == true)
      {
	VectorQuantumNumber* TmpQ = new VectorQuantumNumber();
	(*TmpQ) += new SzQuantumNumber (this->Sz);
	(*TmpQ) += new MomentumQuantumNumber (this->Momentum);
	TmpList += TmpQ;
	return TmpList;
      }
    else
      {
	for (int j = 0; j < this->ChainLength; j++)
	  {
	    VectorQuantumNumber* TmpQ = new VectorQuantumNumber();
	    (*TmpQ) += new SzQuantumNumber (this->Sz);
	    (*TmpQ) += new MomentumQuantumNumber (j);
	    TmpList += TmpQ;
	  }
	return TmpList;
      }
  if (this->FixedSpinProjectionFlag == true)
    {
      int TmpSz = - this->ChainLength;
      for (int i = 0; i <= this->ChainLength; i++)
	{
	  VectorQuantumNumber* TmpQ = new VectorQuantumNumber();
	  (*TmpQ) += new SzQuantumNumber (TmpSz);
	  (*TmpQ) += new MomentumQuantumNumber (this->Momentum);
	  TmpList += TmpQ;
	  TmpSz += 2;
	}
      return TmpList;
    }
  for (int j = 0; j < this->ChainLength; j++)
    {
      int TmpSz = - this->ChainLength;
      for (int i = 0; i <= this->ChainLength; i++)
	{
	  VectorQuantumNumber* TmpQ = new VectorQuantumNumber();
	  (*TmpQ) += new SzQuantumNumber (TmpSz);
	  (*TmpQ) += new MomentumQuantumNumber (j);
	  TmpList += TmpQ;
	  TmpSz += 2;
	}
    }
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* Spin1_2ChainWithTranslations::GetQuantumNumber (int index)
{ 
  if (this->FixedMomentumFlag == false)
    return new SzQuantumNumber (this->TotalSz(index));
  else
    {
      VectorQuantumNumber* TmpQ = new VectorQuantumNumber();
      (*TmpQ) += new SzQuantumNumber (this->TotalSz(index));
      (*TmpQ) += new MomentumQuantumNumber (this->Momentum);
      return TmpQ;
    }
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin1_2ChainWithTranslations::TotalSz (int index)
{
  if (this->FixedSpinProjectionFlag == true)
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

Matrix& Spin1_2ChainWithTranslations::Sxi (int i, Matrix& M)
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

Matrix& Spin1_2ChainWithTranslations::Syi (int i, Matrix& M)
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

Matrix& Spin1_2ChainWithTranslations::Szi (int i, Matrix& M)
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

int Spin1_2ChainWithTranslations::Pij (int i, int j, int state)
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

double Spin1_2ChainWithTranslations::SziSzj (int i, int j, int state)
{  
  unsigned long Mask = ((0x00000001 << i) | (0x00000001 << j));
  unsigned long tmpState = this->ChainDescription[state] & Mask;
  if ((tmpState == 0x00000000) || (tmpState == Mask))
    return 0.25;
  else
    return -0.25;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// translations = number of translations to apply to the resulting state to obtain the true resulting state
// return value = index of resulting state (orbit index)

int Spin1_2ChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& translation)
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
	  return this->FindStateIndex(this->FindCanonicalState((State | (0x00000001 << j)) & ~(0x00000001 << i), translation));
	}
      else
	{
	  coefficient = 0.0;
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

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Spin1_2ChainWithTranslations::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
  switch (q.GetQuantumNumberType())
    {
    case AbstractQuantumNumber::Sz:
      {
	int TmpSz = ((SzQuantumNumber&) q).GetSz();
	if (this->FixedSpinProjectionFlag == true)
	  if (TmpSz == this->Sz)
	    return this->Clone();
	  else
	    return 0;
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
	int* SubspaceOrbitSize = new int [HilbertSubspaceDimension];
	SubspaceDescription[0] = this->ChainDescription[TmpConvArray[0]];
	ConvArray[0] = TmpConvArray[0];
	SubspaceOrbitSize[0] = this->OrbitSize[TmpConvArray[0]];
	for (int i = 1; i < HilbertSubspaceDimension; i++)
	  {
	    SubspaceDescription[i] = this->ChainDescription[TmpConvArray[i]];
	    SubspaceOrbitSize[i] = this->OrbitSize[TmpConvArray[i]];
	    ConvArray[i] = TmpConvArray[i];
	  }
	delete[] TmpConvArray;
	converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, HilbertSubspaceDimension, ConvArray);
	return (AbstractSpinChainWithTranslations*) new Spin1_2ChainWithTranslations (HilbertSubspaceDimension, SubspaceDescription, 
										      SubspaceOrbitSize,
										      this->ChainLength,
										      TmpSz, true, 0, false, this->LookUpTableSize);
      }
    break;
    case AbstractQuantumNumber::Momentum:
      {
      }
      break;
    case (AbstractQuantumNumber::Vector | AbstractQuantumNumber::Momentum | AbstractQuantumNumber::Sz):
      {
	int TmpSz = ((SzQuantumNumber*) ((VectorQuantumNumber&) q)[0])->GetSz();
	int TmpP = ((MomentumQuantumNumber*) ((VectorQuantumNumber&) q)[1])->GetMomentum();
	if (this->FixedSpinProjectionFlag == true)
	  if (TmpSz == this->Sz)
	    return this->Clone();
	  else
	    return 0;
	int HilbertSubspaceDimension = 0;
	int* TmpConvArray = new int [this->HilbertSpaceDimension];
	for (int i = 0; i < this->HilbertSpaceDimension; i++)
	  {
	    if ((this->TotalSz(i) == TmpSz) && ((TmpP % (this->ChainLength / this->OrbitSize[i])) == 0))
	      {
		TmpConvArray[HilbertSubspaceDimension] = i;
		HilbertSubspaceDimension++;	  
	      }
	  }
	if (HilbertSubspaceDimension == 0)
	  {
	    return (AbstractSpinChainWithTranslations*) new Spin1_2ChainWithTranslations ();
	  }
	int* ConvArray = new int [HilbertSubspaceDimension];
	unsigned long* SubspaceDescription = new unsigned long [HilbertSubspaceDimension];
	int* SubspaceOrbitSize = new int [HilbertSubspaceDimension];
	SubspaceDescription[0] = this->ChainDescription[TmpConvArray[0]];
	ConvArray[0] = TmpConvArray[0];
	SubspaceOrbitSize[0] = 1;
	for (int i = 1; i < HilbertSubspaceDimension; i++)
	  {
	    SubspaceDescription[i] = this->ChainDescription[TmpConvArray[i]];
	    SubspaceOrbitSize[i] = 1;
	    ConvArray[i] = TmpConvArray[i];
	  }
	delete[] TmpConvArray;
	converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, HilbertSubspaceDimension, ConvArray);
	return (AbstractSpinChainWithTranslations*) new Spin1_2ChainWithTranslations (HilbertSubspaceDimension, 
										      SubspaceDescription, SubspaceOrbitSize,
										      this->ChainLength,
										      TmpSz, true, 0, false, this->LookUpTableSize);
      }
    default:
      return 0;
    }
  return 0;
}

// find state index
//
// state = state description
// return value = corresponding index

int Spin1_2ChainWithTranslations::FindStateIndex(unsigned long state)
{
  int index = this->LookUpTable[state >> this->LookUpTableShift];
  unsigned long* tmpState = &(this->ChainDescription[index]);
  while ((index < this->HilbertSpaceDimension) && (state != *(tmpState++)))
    index++;
  return index;
    
}

// find canonical state
//
// state = state description
// nbrTranslation = reference on an integer used to store number of translations
// return value = canonical state description

unsigned long Spin1_2ChainWithTranslations::FindCanonicalState(unsigned long state, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = state;
  unsigned long TmpState = state;
  int index = 1;
  while (index < this->ChainLength)
    {
      TmpState = (TmpState >> 1) | ((TmpState & 0x1) << this->ReducedChainLength);
      if (TmpState < CanonicalState)
	{
	  CanonicalState = TmpState;
	  nbrTranslation = index;
	}
      ++index;
    }
  return CanonicalState;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin1_2ChainWithTranslations::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long Mask = 0x00000001;
  for (int k = 0; k < this->ChainLength; k++)    
    {
      if ((this->ChainDescription[state] & Mask) == 0x00000000)
	Str << "- ";
      else
	Str << "+ ";
      Mask <<= 1;
    }
  Str << "  orbit size = "  << this->OrbitSize[state];
  return Str;
}
