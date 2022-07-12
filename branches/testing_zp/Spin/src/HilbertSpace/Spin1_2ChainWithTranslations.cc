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
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;


#define M_SQRT3 1.73205080756888
#define M_SQRT6 2.44948974278318


// default constructor
//

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->ChainDescription = 0;
  this->ChainLength = 0;
  this->Momentum = 0;
  this->ComplementaryStateShift = 0;
  this->Sz = 0;
  this->FixedSpinProjectionFlag = false;
  this->CompatibilityWithMomentum = 0;
  this->RescalingFactors = 0;
  this->NbrStateInOrbit = 0;
}

// constructor for Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (int chainLength, int momentum, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->ComplementaryStateShift = (this->ChainLength - 1) << 1;
  this->Momentum = momentum;

  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->CreatePrecalculationTable();
  this->HilbertSpaceDimension = this->GenerateStates (memorySlice >> 3);
  this->CreateLookUpTable();
}

// constructor for Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (int chainLength, int momentum, int sz, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz;
  this->FixedSpinProjectionFlag = true;
  this->ComplementaryStateShift = (this->ChainLength - 1) << 1;
  this->Momentum = momentum;

  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->CreatePrecalculationTable();
  this->HilbertSpaceDimension = this->GenerateStates (this->Sz, memorySlice >> 3);
  this->CreateLookUpTable();
}

// constructor from pre-constructed datas
//
// hilbertSpaceDimension = Hilbert space dimension
// chainDescription = array describing states
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
// lookUpTableShift = shift to apply to a state to obtain an index to the look-up table 
// complementaryStateShift = shift to apply to move the spin from one end to the other one

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
							int momentum, int sz, bool fixedQuantumNumberFlag, int lookUpTableShift, 
							int complementaryStateShift)
{
  this->Flag.Initialize();
  this->ComplementaryStateShift = complementaryStateShift;
  this->LookUpTableShift = lookUpTableShift;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->ChainDescription = chainDescription;
  this->Momentum = momentum;
  this->Sz = sz;
  this->FixedSpinProjectionFlag = fixedQuantumNumberFlag;
  this->ChainLength = chainLength;
  this->CreatePrecalculationTable();
  this->CreateLookUpTable();
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
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescription = chain.ChainDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
    }
}

// destructor
//

Spin1_2ChainWithTranslations::~Spin1_2ChainWithTranslations () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->ChainDescription;
      delete[] this->LookUpTable;
      delete[] this->CompatibilityWithMomentum;
      for (int i = 1; i <= this->ChainLength; ++i)
	{
	  delete[] this->RescalingFactors[i];
	} 
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainWithTranslations& Spin1_2ChainWithTranslations::operator = (const Spin1_2ChainWithTranslations& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->ChainDescription;
      delete[] this->LookUpTable;
      delete[] this->CompatibilityWithMomentum;
      for (int i = 1; i <= this->ChainLength; ++i)
	{
	  delete[] this->RescalingFactors[i];
	} 
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescription = chain.ChainDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
   }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
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

// return Hilbert space dimension
//
// return value = Hilbert space dimension

int Spin1_2ChainWithTranslations::GetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Spin1_2ChainWithTranslations::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedSpinProjectionFlag == true)
    {
      List<AbstractQuantumNumber*> TmpList2;
      TmpList2 += new PeriodicMomentumQuantumNumber (this->Momentum, this->ChainLength);
      TmpList2 += new SzQuantumNumber (this->Sz);
      TmpList += new VectorQuantumNumber (TmpList2);
    }
  else
    {
      int TmpSz = - 2 * this->ChainLength;
      for (int i = 0; i <= (2 * this->ChainLength); i++)
	{
	  List<AbstractQuantumNumber*> TmpList2;
	  TmpList2 += new PeriodicMomentumQuantumNumber (this->Momentum, this->ChainLength);
	  TmpList2 += new SzQuantumNumber (TmpSz);
	  TmpList += new VectorQuantumNumber (TmpList2);
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
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->Momentum, this->ChainLength);
  TmpList += new SzQuantumNumber (this->TotalSz(index));
  return new VectorQuantumNumber (TmpList);
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
      TmpSz += (State & 0x00000001) << 1;
      State >>= 1;
    }
  return ((2 * TmpSz) - this->ChainLength);
}

// return value of the value of the sum of the square of spin projection on (Oz) 
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

inline double Spin1_2ChainWithTranslations::TotalSzSz (int index)
{  
  return (((double) this->ChainLength) * 0.25);
}

// return value of twice spin projection on (Oz) for a given state
//
// stateDescription = state to which the spin projection has to be evaluated
// return value = twice spin projection on (Oz)

inline int Spin1_2ChainWithTranslations::GetTotalSz (unsigned long stateDescription)
{
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpSz += (stateDescription & 0x00000001);
      stateDescription >>= 1;
    }
  return ((2 * TmpSz) - this->ChainLength);
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

// return index of resulting state from application of P_ij operator on a given state
//
// i = first position
// j = second position
// state = index of the state to be applied on P_ij operator
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::Pij (int i, int j, int state, double& coefficient, int& nbrTranslation)
{  
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpMask = (0x00000001 << i) | (0x00000001 << j);
  unsigned long tmpState2 = tmpState & tmpMask;
  unsigned long tmpState3 = ~tmpState & tmpMask;
  if ((tmpState2 == 0x00000000) || (tmpState3 == 0x00000000))
    return this->HilbertSpaceDimension;
  else
    {
      tmpState = this->FindCanonicalForm((tmpState & ~tmpMask) | tmpState3, nbrTranslation, i);
      if (this->CompatibilityWithMomentum[i] == false)
	return this->HilbertSpaceDimension;
      j = this->FindStateIndex(tmpState);
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[state]][i];
      return j;
    }
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
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
	  State = this->FindCanonicalForm((State | (0x00000001 << j)) & ~(0x00000001 << i), nbrTranslation, i);
	  if (this->CompatibilityWithMomentum[i] == false)
	    return this->HilbertSpaceDimension;
	  j = this->FindStateIndex(State);
	  coefficient = this->RescalingFactors[this->NbrStateInOrbit[state]][i];
	  return j;
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

// return index of resulting state from application of S+_i operator on a given state (only valid if there is no constraint on total Sz)
//
// i = operator position
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::Spi (int i, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long State = this->ChainDescription[state];
  unsigned long TmpMask = (((unsigned long) 0x1) << i);
  if ((State & TmpMask) == 0)
    {
      State = this->FindCanonicalForm((State | TmpMask), nbrTranslation, i);
      if (this->CompatibilityWithMomentum[i] == false)
	return this->HilbertSpaceDimension;
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[state]][i];
      return this->FindStateIndex(State);
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i operator on a given state (only valid if there is no constraint on total Sz)
//
// i = operator position
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::Smi (int i, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long State = this->ChainDescription[state];
  unsigned long TmpMask = (((unsigned long) 0x1) << i);
  if ((State & TmpMask) != 0)
    {
      State = this->FindCanonicalForm((State & (~TmpMask)), nbrTranslation, i);
      if (this->CompatibilityWithMomentum[i] == false)
	return this->HilbertSpaceDimension;
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[state]][i];
      return this->FindStateIndex(State);
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  if (i == j)
    return this->HilbertSpaceDimension;
  unsigned long State = this->ChainDescription[state];
  unsigned long TmpMask = (((unsigned long) 0x1) << j) | (((unsigned long) 0x1) << i);
  if ((State & TmpMask) == 0)
    {
      State = this->FindCanonicalForm((State | TmpMask), nbrTranslation, i);
      if (this->CompatibilityWithMomentum[i] == false)
	return this->HilbertSpaceDimension;
      j = this->FindStateIndex(State);
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[state]][i];
      return j;
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

int Spin1_2ChainWithTranslations::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  if (i == j)
    return this->HilbertSpaceDimension;
  unsigned long State = this->ChainDescription[state];
  unsigned long TmpMask = (((unsigned long) 0x1) << j) | (((unsigned long) 0x1) << i);
  if ((State & TmpMask) == TmpMask)
    {
      State = this->FindCanonicalForm((State & (~TmpMask)), nbrTranslation, i);
      if (this->CompatibilityWithMomentum[i] == false)
	return this->HilbertSpaceDimension;
      j = this->FindStateIndex(State);
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[state]][i];
      return j;
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

int Spin1_2ChainWithTranslations::SpiSpi (int i, int state, double& coefficient, int& nbrTranslation)
{
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_i operator on a given state
//
// i = position of the S- operator
// state = index of the state to be applied on S-_i S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::SmiSmi (int i, int state, double& coefficient, int& nbrTranslation)
{
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long State = this->ChainDescription[state];
  if (((State >> j) & (unsigned long) 0x1) == 0)
    coefficient = 0.5;
  else
    coefficient = -0.5;
  if ((State & (((unsigned long) 0x1) << i)) != 0)
    return this->HilbertSpaceDimension;
  State = this->FindCanonicalForm(State | (((unsigned long) 0x1) << i), nbrTranslation, i);
  if (this->CompatibilityWithMomentum[i] == false)
    return this->HilbertSpaceDimension;
  j = this->FindStateIndex(State);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[state]][i];
  return j;
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long State = this->ChainDescription[state];
  if (((State >> j) & (unsigned long) 0x1) == 0)
    coefficient = 0.5;
  else
    coefficient = -0.5;
  if ((State & (((unsigned long) 0x1) << i)) == 0)
    return this->HilbertSpaceDimension;
  State = this->FindCanonicalForm(State & (~(((unsigned long) 0x1) << i)), nbrTranslation, i);
  if (this->CompatibilityWithMomentum[i] == false)
    return this->HilbertSpaceDimension;
  j = this->FindStateIndex(State);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[state]][i];
  return j;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Spin1_2ChainWithTranslations::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
  if (((VectorQuantumNumber&) q).GetQuantumNumbers().GetNbrElement() != 2)
    return 0;
  if (((((VectorQuantumNumber&) q)[0]))->GetQuantumNumberType() != AbstractQuantumNumber::PeriodicMomentum)
    return 0;
  if (((((VectorQuantumNumber&) q)[1]))->GetQuantumNumberType() != AbstractQuantumNumber::Sz)
    return 0;
  if (this->Momentum != ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[0]))->GetMomentum())
    return 0;
  if (this->FixedSpinProjectionFlag == true)
    if (this->Sz != ((SzQuantumNumber*) (((VectorQuantumNumber&) q)[1]))->GetSz())
      return 0;
    else
      return this;
  int TmpSz = ((SzQuantumNumber*) (((VectorQuantumNumber&) q)[1]))->GetSz();
  if ((TmpSz < (-2 * this->ChainLength)) || (TmpSz > (2 * this->ChainLength)))
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
  SubspaceDescription[0] = this->ChainDescription[TmpConvArray[0]];
  ConvArray[0] = TmpConvArray[0];
//  unsigned long TestMask = this->LookUpTableMask;
  for (int i = 1; i < HilbertSubspaceDimension; i++)
    {
      SubspaceDescription[i] = this->ChainDescription[TmpConvArray[i]];
      ConvArray[i] = TmpConvArray[i];
    }
  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, HilbertSubspaceDimension, ConvArray);
  return new Spin1_2ChainWithTranslations (HilbertSubspaceDimension, SubspaceDescription, this->ChainLength,
					 this->Momentum, TmpSz, true, 
					 this->LookUpTableShift, this->ComplementaryStateShift);
}

// find the canonical form of a state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// return value = canonical form of the state

inline unsigned long Spin1_2ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  int index = 1;  
  while (index < this->ChainLength)
    {
      stateDescription = (stateDescription >> 2) | ((stateDescription & 0x3) << this->ComplementaryStateShift);
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = index;
	}
      ++index;
    }
  return CanonicalState;
}

// find the canonical form of a state and find how many translations are needed to obtain the same state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
// return value = canonical form of the state

inline unsigned long Spin1_2ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation, int& nbrTranslationToIdentity)
{
  nbrTranslation = 0;
  nbrTranslationToIdentity = 1;
  unsigned long CanonicalState = stateDescription;
  unsigned long ReferenceState = stateDescription;
  stateDescription = (stateDescription >> 2) | ((stateDescription & 0x3) << this->ComplementaryStateShift);
  while ((ReferenceState != stateDescription) && (nbrTranslationToIdentity < this->ChainLength))
    {
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = nbrTranslationToIdentity;
	}
      stateDescription = (stateDescription >> 2) | ((stateDescription & 0x3) << this->ComplementaryStateShift);
      ++nbrTranslationToIdentity;
    }
  return CanonicalState;
}

// find how many translations are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

inline int Spin1_2ChainWithTranslations::FindNumberTranslation(unsigned long stateDescription)
{
  unsigned long TmpState = (stateDescription >> 2) | ((stateDescription & 0x3) << this->ComplementaryStateShift);
  int index = 1;  
  while (TmpState != stateDescription)
    {
      TmpState = (TmpState >> 2) | ((TmpState & 0x3) << this->ComplementaryStateShift);
      ++index;
    }
  return index;
}

// find state index
//
// state = state description
// return value = corresponding index

inline int Spin1_2ChainWithTranslations::FindStateIndex(unsigned long state)
{
  unsigned long MidPos = state >> this->LookUpTableShift;
  unsigned long LowPos = this->LookUpTable[MidPos];
  unsigned long HighPos = this->LookUpTable[MidPos + 1];
  while (this->ChainDescription[HighPos] != state)
    {
      MidPos = (HighPos + LowPos) >> 1;
      if (this->ChainDescription[MidPos] >= state)
	HighPos = MidPos;
      else
	LowPos = MidPos;
    }
  return HighPos;   
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
  return Str;
}

// generate all states with no constraint on total Sz
//
// memorySlice = maximum amount of memory (in unsigned long unit) that can be allocated to partially evalauted the states
// return value = number of generated states

int Spin1_2ChainWithTranslations::GenerateStates(long memorySlice) 
{
  long MaximumNbrState = 1 << this->ChainLength;
  List<unsigned long*> TmpGeneratedStateList;
  List<int*> TmpNbrStateInOrbitList;
  List<long> TmpNbrGeneratedStateList;
  int CurrentNbrStateInOrbit;
  unsigned long TmpState = 0;
  unsigned long TmpState2;
  int NbrTranslation = 0;
  int TmpHilbertSpaceDimension = 0;
  while (MaximumNbrState > 0)
    {
      unsigned long* TmpGeneratedStates = new unsigned long [memorySlice];
      int* TmpNbrStateInOrbit = new int [memorySlice];
      long Pos = 0;
      while ((Pos < memorySlice) && (MaximumNbrState > 0))
	{
	  TmpState2 = this->FindCanonicalForm(TmpState, NbrTranslation);
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpState2);
	  if ((NbrTranslation == 0) && (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true))
	    {
	      TmpGeneratedStates[Pos] = TmpState2;
	      TmpNbrStateInOrbit[Pos] = CurrentNbrStateInOrbit;
	      ++Pos;
	    }
	  ++TmpState;
	  --MaximumNbrState;
	}
      if (Pos > 0)
	{
	  TmpGeneratedStateList += TmpGeneratedStates;
	  TmpNbrStateInOrbitList += TmpNbrStateInOrbit;
	  TmpNbrGeneratedStateList += Pos;
	  TmpHilbertSpaceDimension += Pos;
	}
    }
  this->ChainDescription = SmartMergeArrayListIntoArray(TmpGeneratedStateList, TmpNbrGeneratedStateList);
  this->NbrStateInOrbit = SmartMergeArrayListIntoArray(TmpNbrStateInOrbitList, TmpNbrGeneratedStateList);

  return TmpHilbertSpaceDimension;
}

// generate all states with constraint on total Sz
//
// sz = twice the sz value
// memorySlice = maximum amount of memory (in unsigned long unit) that can be allocated to partially evalauted the states
// return value = number of generated states

int Spin1_2ChainWithTranslations::GenerateStates(int sz, long memorySlice) 
{
  long MaximumNbrState = 1 << this->ChainLength;
  List<unsigned long*> TmpGeneratedStateList;
  List<int*> TmpNbrStateInOrbitList;
  List<long> TmpNbrGeneratedStateList;
  int CurrentNbrStateInOrbit;
  unsigned long TmpState = 0;
  unsigned long TmpState2;
  int NbrTranslation = 0;
  int TmpHilbertSpaceDimension = 0;
  while (MaximumNbrState > 0)
    {
      unsigned long* TmpGeneratedStates = new unsigned long [memorySlice];
      int* TmpNbrStateInOrbit = new int [memorySlice];
      long Pos = 0;
      while ((Pos < memorySlice) && (MaximumNbrState > 0))
	{
	  TmpState2 = this->FindCanonicalForm(TmpState, NbrTranslation);
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpState2);
	  if ((NbrTranslation == 0) && (this->GetTotalSz(TmpState2) == sz) && (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true))
	    {
	      TmpGeneratedStates[Pos] = TmpState2;
	      TmpNbrStateInOrbit[Pos] = CurrentNbrStateInOrbit;
	      ++Pos;
	    }
	  ++TmpState;
	  --MaximumNbrState;
	}
      if (Pos > 0)
	{
	  TmpGeneratedStateList += TmpGeneratedStates;
	  TmpNbrStateInOrbitList += TmpNbrStateInOrbit;
	  TmpNbrGeneratedStateList += Pos;
	  TmpHilbertSpaceDimension += Pos;
	}
    }
  this->ChainDescription = SmartMergeArrayListIntoArray(TmpGeneratedStateList, TmpNbrGeneratedStateList);
  this->NbrStateInOrbit = SmartMergeArrayListIntoArray(TmpNbrStateInOrbitList, TmpNbrGeneratedStateList);
  return TmpHilbertSpaceDimension;
}


// create precalculation tables
//

void Spin1_2ChainWithTranslations::CreatePrecalculationTable()
{
  this->CompatibilityWithMomentum = new bool [this->ChainLength + 1];
  for (int i = 0; i <= this->ChainLength; ++i)
    if (((i * this->Momentum) % this->ChainLength) == 0)
      this->CompatibilityWithMomentum[i] = true;
    else
      this->CompatibilityWithMomentum[i] = false;

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


// create look-up table used to speed up index search
//

void Spin1_2ChainWithTranslations::CreateLookUpTable()
{
  int TmpHilbertSpaceDimension = this->HilbertSpaceDimension;
  // create the look-up table
  unsigned long Max = ((unsigned long) 1) << ((this->ChainLength << 1) - this->LookUpTableShift + 1);
  this->LookUpTable = new long [Max + 1];
  long LowPos;
  long MidPos;
  long HighPos;
  unsigned long Max2 = (this->ChainDescription[TmpHilbertSpaceDimension - 1]) >> this->LookUpTableShift;
  for (unsigned long i = 0; i <= Max2; ++i)
    {
      LowPos = 0;
      HighPos = TmpHilbertSpaceDimension - 1;
      while ((HighPos - LowPos) > 1)
	{
	  MidPos = (HighPos + LowPos) >> 1;
	  if (this->ChainDescription[MidPos] >= (i << this->LookUpTableShift))
	    HighPos = MidPos;
	  else
	    LowPos = MidPos;
	}      
      this->LookUpTable[i] = LowPos;
    }
  --TmpHilbertSpaceDimension;
  for (unsigned long i = Max2 + 1; i <= Max; ++i)    
    this->LookUpTable[i] = TmpHilbertSpaceDimension;
  ++TmpHilbertSpaceDimension;
}

