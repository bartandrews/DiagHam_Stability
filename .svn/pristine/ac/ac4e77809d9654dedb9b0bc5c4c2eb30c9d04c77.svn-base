////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                    class author: Cecile Repellin                           //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with a fixed Sz value                //
//               and a pseudospin 1/2 (not a conserved quantity)              //
//                                                                            //
//                        last modification : 11/06/2016                      //
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


#include "HilbertSpace/Spin1_2ChainWithPseudospin.h"
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

Spin1_2ChainWithPseudospin::Spin1_2ChainWithPseudospin ()
{
}

// constructor for complete Hilbert space with a given total spin projection Sz
//
// chainLength = number of spin 1/2
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Spin1_2ChainWithPseudospin::Spin1_2ChainWithPseudospin (int chainLength, int sz, int memorySize) 
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

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->ChainLength);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = this->GenerateStates ((this->Sz + this->ChainLength) >> 1, this->ChainLength - 1, 0l);
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->GenerateLookUpTable(memorySize);
  
//   double TmpCoef;
//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//   {
//     this->PrintState(cout,i) << endl ;
//     cout << (this->JOffDiagonali(1, i, TmpCoef)) << " ";
// //     if ((this->SmiSpj(0, 1, i, TmpCoef)) < this->HilbertSpaceDimension)
// //       cout << TmpCoef;
//     cout << endl;
//   }
//   cout << endl;
  this->Flag.Initialize();
  
//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//     this->PrintState(cout, i) << endl;
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainWithPseudospin::Spin1_2ChainWithPseudospin (const Spin1_2ChainWithPseudospin& chain)
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
    }
}

// destructor
//

Spin1_2ChainWithPseudospin::~Spin1_2ChainWithPseudospin () 
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->LookUpTableShift != 0)
	{
	  delete[] this->LookUpTableShift;
	  for (int i = 0; i < this->ChainLength; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;
	}
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainWithPseudospin& Spin1_2ChainWithPseudospin::operator = (const Spin1_2ChainWithPseudospin& chain)
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
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainWithPseudospin::Clone()
{
  return new Spin1_2ChainWithPseudospin (*this);
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin1_2ChainWithPseudospin::TotalSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->Sz;
  unsigned long State = this->StateDescription[index];
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpSz += ((State & 0x1ul) << 1);
      State >>= 2;
    }
  TmpSz -= this->ChainLength;
  return TmpSz;
}

// compute the parity (prod_i Sz_i) for a given state
//
// state = index of the state to be applied on Sz_i operator
// return value = 0 if prod_i Sz_i = 1, 1 if prod_i Sz_i = -1

unsigned long Spin1_2ChainWithPseudospin::GetParity (int state)
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

int Spin1_2ChainWithPseudospin::TranslateState (int nbrTranslations, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  TmpState = (((TmpState & ((0x1ul << (this->ChainLength - nbrTranslations)) - 1ul)) << nbrTranslations)
	      | (TmpState >> (this->ChainLength - nbrTranslations)));
  return this->FindStateIndex(TmpState);
}


// find state index
//
// stateDescription = state description
// return value = corresponding index

int Spin1_2ChainWithPseudospin::FindStateIndex(unsigned long stateDescription)
{
  int CurrentMaximumUpPosition = 2 * this->ChainLength - 1;
  while ((((stateDescription >> CurrentMaximumUpPosition) & 0x1ul) == 0x0ul) && (CurrentMaximumUpPosition > 0))
    --CurrentMaximumUpPosition;
  long PosMax = stateDescription >> this->LookUpTableShift[CurrentMaximumUpPosition];
  long PosMin = this->LookUpTable[CurrentMaximumUpPosition][PosMax];
  PosMax = this->LookUpTable[CurrentMaximumUpPosition][PosMax + 1];
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
    {
      if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
	return this->HilbertSpaceDimension;
      else
	return PosMin;
    }
}

// evaluate Hilbert space dimension
//
// sz = twice the Sz value
// nbrSites = number of sites
// return value = Hilbert space dimension

long Spin1_2ChainWithPseudospin::EvaluateHilbertSpaceDimension(int sz, int nbrSites)
{
  BinomialCoefficients TmpCoefficients (nbrSites);
  long PseudospinContribution = 1l << ((long) this->ChainLength);
  return (TmpCoefficients(nbrSites, (nbrSites + sz) >> 1) * PseudospinContribution);
}

// generate all states
// 
// nbrSpinUp = number of spin up
// currentPosition = current position to consider in the chain
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long Spin1_2ChainWithPseudospin::GenerateStates(int nbrSpinUp, int currentPosition, long pos)
{
  if ((nbrSpinUp == 1) && (currentPosition == 0))
  {
//     this->StateDescription[pos + 1] = 0x1ul;
//     this->StateDescription[pos] = 0x3ul;
    this->StateDescription[pos + 1] = 0x2ul;
    this->StateDescription[pos] = 0x3ul;
    pos += 2;
    return pos;
  }
    
  if (nbrSpinUp == 0)
    {
      if (currentPosition < 0)
	return pos;
      else
      {
	if (currentPosition == 0)
	{
// 	  this->StateDescription[pos + 1] = 0x0ul;
// 	  this->StateDescription[pos] = (0x1ul << 1l);
	  this->StateDescription[pos + 1] = 0x0ul;
	  this->StateDescription[pos] = 0x1ul;
	  pos += 2;
	  return pos;
	}
	else
	{
	  long TmpPos = this->GenerateStates(nbrSpinUp, currentPosition - 1, pos);
// 	  unsigned long Mask = 0x1ul << (2 * currentPosition + 1);
	  unsigned long Mask = 0x1ul << (2 * currentPosition);
	  for (long i = pos; i < TmpPos; i++)
	  {
	    this->StateDescription[i + (TmpPos - pos)] = this->StateDescription[i];
	    this->StateDescription[i] |= Mask;
	  }
	  return (2 * TmpPos - pos);
	}
      }
    }
  if (currentPosition < (nbrSpinUp - 1))
    return pos;
  int ReducedCurrentPosition = currentPosition - 1;
  long TmpPos = this->GenerateStates(nbrSpinUp - 1, ReducedCurrentPosition, pos);
//   unsigned long Mask = 0x1ul << (2 * currentPosition);
  unsigned long Mask = 0x1ul << (2 * currentPosition);
  for (long i = pos; i < TmpPos; i++)
  {
//     this->StateDescription[i + (TmpPos - pos)] = this->StateDescription[i];
//     this->StateDescription[i] |= Mask;
//     this->StateDescription[i] |= (Mask << 0x1ul);
//     this->StateDescription[i + (TmpPos - pos)] |= Mask;
    this->StateDescription[i + (TmpPos - pos)] = this->StateDescription[i];
    this->StateDescription[i] |= Mask;
    this->StateDescription[i] |= (Mask << 0x1ul);
    this->StateDescription[i + (TmpPos - pos)] |= (Mask << 0x1ul);
  }
  pos = 2 * TmpPos - pos;
  TmpPos =  this->GenerateStates(nbrSpinUp, ReducedCurrentPosition, pos);
  for (long i = pos; i < TmpPos; i++)
  {
    this->StateDescription[i + (TmpPos - pos)] = this->StateDescription[i];
//     this->StateDescription[i] |= (Mask << 0x1ul);
    this->StateDescription[i] |= (Mask);
  }
  return (2 * TmpPos - pos);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin1_2ChainWithPseudospin::PrintState (ostream& Str, int state)
{  
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long Mask = 0x1ul;
  for (int k = 0; k < this->ChainLength; k++)    
    {
      if ((this->StateDescription[state] & Mask) == 0x0ul)
	Str << "0 ";
      else
	Str << "1 ";
      Mask <<= 1;
      if ((this->StateDescription[state] & Mask) == 0x0ul)
	Str << "- ";
      else
	Str << "+ ";
      Mask <<= 1;
    }
//   Str << " " << hex << this->StateDescription[state] << dec;
//   Str << " " << this->FindStateIndex(this->StateDescription[state]);
  return Str;
}


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table
// stateMask = an optional mask to apply to each state to focus on the relevant bits

void Spin1_2ChainWithPseudospin::GenerateLookUpTable(unsigned long memory, unsigned long stateMask)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * (2 * this->ChainLength));
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > 2 * this->ChainLength)
    this->MaximumLookUpShift = 2 * this->ChainLength;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [2 * this->ChainLength];
  this->LookUpTableShift = new int [2 * this->ChainLength];
  for (int i = 0; i < 2 * this->ChainLength; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentPosition = 2 * this->ChainLength - 1;
  while ((((this->StateDescription[0] >> CurrentPosition) & 0x1ul) == 0x0ul) && (CurrentPosition > 0))
    --CurrentPosition;
  int* TmpLookUpTable = this->LookUpTable[CurrentPosition];
  if (CurrentPosition < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentPosition] = 0;
  else
    this->LookUpTableShift[CurrentPosition] = CurrentPosition + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentPosition];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = (this->StateDescription[0] & stateMask) >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      int TmpCurrentPosition = CurrentPosition;
      while ((((this->StateDescription[i] >> TmpCurrentPosition) & 0x1ul) == 0x0ul) && (TmpCurrentPosition > 0))
	--TmpCurrentPosition;
      if (CurrentPosition != TmpCurrentPosition)
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
 	  CurrentPosition = TmpCurrentPosition;
	  TmpLookUpTable = this->LookUpTable[CurrentPosition];
	  if (CurrentPosition < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentPosition] = 0;
	  else
	    this->LookUpTableShift[CurrentPosition] = CurrentPosition + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentPosition];
	  TmpLookUpTableValue = (this->StateDescription[i] & stateMask) >> CurrentShift;
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
	  TmpLookUpTableValue = (this->StateDescription[i] & stateMask) >> CurrentShift;
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



// return the Bosonic Occupation of a given state in the basis
//
// index = index of the state in the basis
// finalState = reference on the array where the monomial representation has to be stored

void Spin1_2ChainWithPseudospin::GetBosonicOccupation (unsigned int index, int * finalState)
{
  for (int i = 0; i < this->ChainLength; i++)
    {
      finalState[i] = (this->StateDescription[index] >> ((unsigned long) (2 * i)) )& 0x1ul;
    }
}


// return eigenvalue of Sz_i Sz_j associated to a given state (acts only on spin part of many-body state)
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double Spin1_2ChainWithPseudospin::SziSzj (int i, int j, int state)
{  
  unsigned long Mask = ((0x1ul << (2*i + 1)) | (0x1ul << (2*j + 1)));
  unsigned long tmpState = this->StateDescription[state] & Mask;
  if ((tmpState == 0x0ul) || (tmpState == Mask))
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
// return value = index of resulting state

int Spin1_2ChainWithPseudospin::SmiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= (2*i + 1);
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= (2*j + 1); 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex((State | (0x1ul << (2*j + 1))) & ~(0x1ul << (2*i + 1)));
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

int Spin1_2ChainWithPseudospin::SpiSmj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= (2*j + 1);
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= (2*i + 1); 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  coefficient = 1.0;
	  return this->FindStateIndex((State | (0x1ul << (2*i + 1))) & ~(0x1ul << (2*j + 1)));
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

// operator acting on pseudospin on site i (off-diagonal part)
//
// i = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospin::JOffDiagonali (int i, int state, double& coefficient)
{
  unsigned long State = this->StateDescription[state];
  coefficient = 1.0;
  unsigned long Tmp = (State >> (2*i)) & 0x1ul;
  if (Tmp == 0x0ul)
    return this->FindStateIndex(State | (0x1ul << (2*i)));
  else
    return this->FindStateIndex(State & ~(0x1ul << (2*i)));
}

// operator acting on pseudospin on site i (off-diagonal part)
//
// i = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coupling = array where the coupling coefficients are stored
// return value = numerical coefficient

double Spin1_2ChainWithPseudospin::JDiagonali (int i, int state, double* coupling)
{
  int Tmp = (int) ((this->StateDescription[state] >> (2*i)) & 0x1ul);
  return coupling[Tmp];
}

// convert a state defined on a lattice with a number of sites equals to a multiple of three
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

RealVector Spin1_2ChainWithPseudospin::ProjectToEffectiveSubspaceThreeToOne(ComplexVector& state, AbstractSpinChain* space)
{
  Spin1_2ChainNew* TmpSpace = (Spin1_2ChainNew*) space;
  RealVector TmpVector (this->HilbertSpaceDimension, true);
  RealVector TmpVector1 (this->HilbertSpaceDimension, true);
  
  cout << "Dimension full space = " << (TmpSpace->HilbertSpaceDimension) << endl;
  cout << "Dimension projected space = " << (this->HilbertSpaceDimension) << endl;
  unsigned long TmpState;
  unsigned long TmpStateEffective;
  
  int tmpStateI;
  int tmpStateJ;
  int tmpStateK;
  
  int tmpStateLocalSz;
  int tmpStateLocalPseudospin;
  int TmpIndexEffective;
  int localSz;
  int localSpinStructure;
  double Coefficient;
  
  int* TmpFlag = new int[this->HilbertSpaceDimension];
  double** Isometry = new double* [2];
  for (int i = 0; i < 2; ++i)
    Isometry[i] = new double [3];
  Isometry[1][0] = 1.0 / sqrt(6.0);
  Isometry[1][1] = 1.0 / sqrt(6.0);
  Isometry[1][2] = -2.0 / sqrt(6.0);
  Isometry[0][0] = 1.0 / sqrt(2.0);
  Isometry[0][1] = -1.0 / sqrt(2.0);
  Isometry[0][2] = 0.0;
  
  for (int j = 0; j < TmpSpace->HilbertSpaceDimension; ++j)
  {
    Coefficient = state[j].Re;
//     cout << j << " " << Coefficient << " " << state[j].Im << endl;
    TmpState = TmpSpace->StateDescription[j];
    int k = TmpSpace->ChainLength - 1;
    TmpVector1.ClearVector();
    for (int i = 0; i < this->HilbertSpaceDimension; ++i)
      TmpFlag[i] = 0;
    while ((fabs(Coefficient) > 1.0e-15) && (k > 0))
    {
      tmpStateK = (int) ((TmpState >> k) & 0x1ul);
      tmpStateJ = (int) ((TmpState >> (k - 1)) & 0x1ul);
      tmpStateI = (int) ((TmpState >> (k - 2)) & 0x1ul);
      
      localSz = tmpStateI + tmpStateJ + tmpStateK - 1;
      
      if ((localSz == 0) || (localSz == 1))
      {
	TmpIndexEffective = (k + 1)/3 - 1;
	
	for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	  {
	    TmpStateEffective = this->StateDescription[i];
	    tmpStateLocalSz = (int) ((TmpStateEffective >> (2 * TmpIndexEffective + 1)) & 0x1ul);
	    if (tmpStateLocalSz == localSz)
	    {
	      tmpStateLocalPseudospin = (int) ((TmpStateEffective >> (2 * TmpIndexEffective)) & 0x1ul);
	      if (tmpStateI == tmpStateJ)
		localSpinStructure = 0;
	      else
	      {
		if (tmpStateI == tmpStateK)
		  localSpinStructure = 1;
		else
		  localSpinStructure = 2;
	      }
// 	      cout << i << " " << (Coefficient * Isometry[tmpStateLocalPseudospin][localSpinStructure]) << endl;
	      if (TmpFlag[i] == 0)
	      {
		TmpVector1[i] = Isometry[tmpStateLocalPseudospin][localSpinStructure];
		TmpFlag[i] = 1;
	      }
	      else
		TmpVector1[i] *= Isometry[tmpStateLocalPseudospin][localSpinStructure];
	    }
	    else
	    {
	      TmpVector1[i] = 0.0;
	      TmpFlag[i] = 1;
	    }
	  }
      }
      else
      {
	k = 0;
	Coefficient = 0.0;
      }
      k -= 3;
    }
   TmpVector.AddLinearCombination(Coefficient, TmpVector1);
  }

//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//     cout << TmpVector[i] << endl;
  for (int i = 0; i < 2; ++i)
    delete[] Isometry[i];
  delete[] Isometry;
  delete[] TmpFlag;
  return TmpVector;
}


// convert a state from a SU(2) basis to another one, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void Spin1_2ChainWithPseudospin::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, RealMatrix oneBodyBasis, AbstractSpinChain* space, 
							 long firstComponent, long nbrComponents)
{
  Spin1_2ChainNew* TmpSpace = (Spin1_2ChainNew*) space;
  int* TmpSpinIndices = new int [this->ChainLength];
  int* TmpPseudospinIndices = new int [this->ChainLength];
  int* TmpPseudospinIndices2 = new int [this->ChainLength];
  targetState.ClearVector();
  for (long i = 0; i < (TmpSpace->HilbertSpaceDimension); ++i)
    {
      unsigned long TmpState = TmpSpace->StateDescription[i];
      unsigned long Tmp;
      unsigned long Tmp1;
      int TmpIndex = 0;
      bool flag = true;
      int j;
      while ((TmpIndex < this->ChainLength) && (flag = true))
	{
	  j = this->ChainLength - TmpIndex - 1;
	  Tmp = (TmpState >> (3 * TmpIndex)) & 0x7ul;
	  Tmp1 = (Tmp & 0x1ul) + ((Tmp >> 1) & 0x1ul) + ((Tmp >> 2) & 0x1ul);
// 	  cout << TmpState << " " << (3*j) << " " << (TmpState >> (3 * j)) << " " << Tmp << " " << Tmp1 << endl;
	  if ((Tmp == 0x7ul) || (Tmp == 0x0ul))
	  {
	    flag = false;
	    TmpIndex = this->ChainLength;
	  }
	  else
	  {
	    if (Tmp1 == 0x1ul)
	      TmpSpinIndices[TmpIndex] = 0;
	    if (Tmp1 == 0x2ul)
	      TmpSpinIndices[TmpIndex] = 1;
	    
	    if ((Tmp & 0x1ul) == ((Tmp >>1) & 0x1ul))
	      TmpPseudospinIndices[TmpIndex] = 0;
	    if ((Tmp & 0x1ul) == ((Tmp >>2) & 0x1ul))
	      TmpPseudospinIndices[TmpIndex] = 1;
	    if (((Tmp >>1) & 0x1ul) == ((Tmp >>2) & 0x1ul))
	      TmpPseudospinIndices[TmpIndex] = 2;
	    
	    ++TmpIndex;
	  }
	  
	  
// 	  cout << TmpState << " " << j << " " << TmpIndex << " " << flag << endl;
	}

      if (flag == true)
      {
// 	for (int k = 0; k < this->ChainLength; ++k)
// 	  cout <<  TmpPseudospinIndices[k] << endl;
	this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpSpinIndices, TmpPseudospinIndices, TmpPseudospinIndices2, oneBodyBasis);
      }
//     cout << endl;
    }
  delete[] TmpSpinIndices;
  delete[] TmpPseudospinIndices;
  delete[] TmpPseudospinIndices2;
}

// recursive part of the convertion from a state from a SU(2) basis to another one, transforming the one body basis in each momentum sector
//
// targetState = vector where the transformed state has to be stored
// coefficient = current coefficient to assign
// position = current particle consider in the n-body state
// momentumIndices = array that gives the momentum partition of the initial n-body state
// initialSpinIndices = array that gives the spin dressing the initial n-body state
// currentSpinIndices = array that gives the spin dressing the current transformed n-body state
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector

void Spin1_2ChainWithPseudospin::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
								int position, int* spinIndices, int* initialPseudospinIndices, int* currentPseudospinIndices, RealMatrix oneBodyBasis) 
{
//   cout << position << " : " << endl;
//   for (int i = 0; i < position; ++i)
//     cout << currentSpinIndices[i] << " ";
//   cout << endl;
  if (position == this->ChainLength)
    {
      unsigned long TmpState = 0x0ul;
      unsigned long Mask = 0x0ul;
      for (int i = 0; i < this->ChainLength; ++i)
	{
	  Mask =  (((currentPseudospinIndices[i]) << (2*i)) | ((spinIndices[i]) << (2*i + 1))); // Mask = 00...0100...0 : one fermion state in the second quantized basis
// 	  cout << i << " mask = " << ((currentPseudospinIndices[i]) << (2*i)) << " " << ((spinIndices[i]) << (2*i + 1)) << " " << Mask << endl;
	  if ((TmpState & Mask) != 0x0ul)
	    return;
	  TmpState |= Mask; //set bit corresponding to the current fermion state to 1 in TmpState
	}
      int Index = this->FindStateIndex(TmpState);
      if (Index < this->HilbertSpaceDimension)
	{
	  targetState[Index] += coefficient;
	}
      return;      
    }
  else
    {
      currentPseudospinIndices[position] = 0;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[0][initialPseudospinIndices[position]]), position + 1, spinIndices, initialPseudospinIndices, currentPseudospinIndices, oneBodyBasis);
      currentPseudospinIndices[position] = 1;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[1][initialPseudospinIndices[position]]), position + 1, spinIndices, initialPseudospinIndices, currentPseudospinIndices, oneBodyBasis);
    }
}


// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1_2ChainWithPseudospin::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Caution: using dummy operator Spin1_2ChainWithPseudospin::SmiSpj" << endl;
  return 0;
}


// operator acting on pseudospin on site i (off-diagonal part)
//
// i = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospin::JOffDiagonali (int i, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Caution: using dummy operator Spin1_2ChainWithPseudospin::JOffDiagonali" << endl;
  return 0;
}

// operator acting on pseudospin on site i (off-diagonal part)
//
// i = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospin::JoffiJoffj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Caution: using dummy operator Spin1_2ChainWithPseudospin::JoffiJoffj" << endl;
  return 0;
}


// operator acting on pseudospin on site i (off-diagonal) and j(diagonal part)
//
// i = position of pseudospin operator
// j = position of pseudospin operator
// state = index of the state to be applied on JAi operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospin::JoffiJj (int i, int j, int state, double* coupling, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Caution: using dummy operator Spin1_2ChainWithPseudospin::JoffiJj" << endl;
  return 0;
}

  
// operator acting on pseudospin on site i (diagonal) and j (diagonal part) and spin (SpSm)
//
// i = position of first spin*pseudospin operator
// j = position of spin*pseudospin operator
// state = index of the state to which operator has to be applied
// couplingI = array of coefficients characterizing the diagonal pseudospin coupling on site i
// couplingJ = array of coefficients characterizing the diagonal pseudospin coupling on site J
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospin::SmiSpjJiJj (int i, int j, int state, double* couplingI, double* couplingJ, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Caution: using dummy operator Spin1_2ChainWithPseudospin::SmiSpjJiJj" << endl;
  return 0;
}
  
  
// operator acting on pseudospin on site i (off-diagonal) and j (diagonal part) and spin (SpSm)
//
// i = position of first spin*pseudospin operator
// j = position of spin*pseudospin operator
// state = index of the state to which operator has to be applied
// coupling = array of coefficients characterizing the diagonal pseudospin coupling on site j
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting state

int Spin1_2ChainWithPseudospin::SmiSpjJoffiJj (int i, int j, int state, double* coupling, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Caution: using dummy operator Spin1_2ChainWithPseudospin::SmiSpjJoffiJj" << endl;
  return 0;
}
  
// operator acting on pseudospin on site i (diagonal) and j (off-diagonal part) and spin (SpSm)
//
// i = position of first spin*pseudospin operator
// j = position of spin*pseudospin operator
// state = index of the state to which operator has to be applied
// coupling = array of coefficients characterizing the diagonal pseudospin coupling on site i
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting stateDescription

int Spin1_2ChainWithPseudospin::SmiSpjJiJoffj (int i, int j, int state, double* coupling, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Caution: using dummy operator Spin1_2ChainWithPseudospin::SmiSpjJiJoffj" << endl;
  return 0;
}
  
// operator acting on pseudospin on site i (off-diagonal) and j (off-diagonal part) and spin (SpSm)
//
// i = position of first spin*pseudospin operator
// j = position of spin*pseudospin operator
// state = index of the state to which operator has to be applied
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of the resulting stateDescription

int Spin1_2ChainWithPseudospin::SmiSpjJoffiJoffj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Caution: using dummy operator Spin1_2ChainWithPseudospin::SmiSpjJoffiJoffj" << endl;
  return 0;
}