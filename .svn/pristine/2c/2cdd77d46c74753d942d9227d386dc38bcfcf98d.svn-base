////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of trapped bosons                         //
//                                                                            //
//                        last modification : 04/06/2002                      //
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


#include "config.h"
#include "HilbertSpace/TrappedBosons.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"

#include <math.h>
#include <iostream>


using std::cout;
using std::endl;


// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value

TrappedBosons::TrappedBosons (int nbrBosons, int totalLz)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->TotalLz, this->TotalLz);
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrBosons, this->TotalLz, this->TotalLz, this->TotalLz, 0);
  this->KeyMultiplicationTable = new int [this->TotalLz + 1];
//  this->ErasthothenesSlieve(10000000, this->TotalLz + 1, this->IncNbrBosons);
  this->GenerateLookUpTable(0);
#ifdef __DEBUG__
  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += (this->StateLzMax[i] + 1) * sizeof(int) + sizeof(int*);
  UsedMemory += (this->TotalLz + 1) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += ((this->TotalLz + 1) * this->IncNbrBosons) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

TrappedBosons::TrappedBosons(const TrappedBosons& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->Flag = bosons.Flag;
}

// destructor
//

TrappedBosons::~TrappedBosons ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

TrappedBosons& TrappedBosons::operator = (const TrappedBosons& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
    }
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->Flag = bosons.Flag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* TrappedBosons::Clone()
{
  return new TrappedBosons(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> TrappedBosons::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* TrappedBosons::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* TrappedBosons::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int TrappedBosons::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int LzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  if ((n1 > LzMax) || (n2 > LzMax) || (State[n1] == 0) || (State[n2] == 0) || ((n1 == n2) && (State[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = LzMax;
  if (LzMax < m1)
    NewLzMax = m1;
  if (LzMax < m2)
    NewLzMax = m2;
  int* TmpState = new int [NewLzMax + 1];
  int i = 0;
  for (; i <= LzMax; ++i)
    TmpState[i] = State[i];
  for (; i <= NewLzMax; ++i)
    TmpState[i] = 0;
  coefficient = TmpState[n2];
  --TmpState[n2];
  coefficient *= TmpState[n1];
  --TmpState[n1];
  ++TmpState[m2];
  coefficient *= TmpState[m2];
  ++TmpState[m1];
  coefficient *= TmpState[m1];
  coefficient = sqrt(coefficient);
  while (TmpState[NewLzMax] == 0)
    --NewLzMax;
  int DestIndex = this->FindStateIndex(TmpState, NewLzMax);
  delete[] TmpState;
  return DestIndex;
}

// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int TrappedBosons::FindStateIndex(int* stateDescription, int lzmax)
{
  int TmpKey = this->GenerateKey(stateDescription, lzmax);
  int i;
  int* TmpStateDescription;
  int Start = this->LzMaxPosition[lzmax * this->IncNbrBosons + stateDescription[lzmax]];
  int* TmpKey2 = &(this->Keys[Start]);
  for (; Start < this->HilbertSpaceDimension; ++Start)    
    {
      if ((*TmpKey2) == TmpKey)
	{
	  i = 0;
	  TmpStateDescription = this->StateDescription[Start];
	  while (i <= lzmax)
	    {
	      if (stateDescription[i] != TmpStateDescription[i])
		i = lzmax + 2;
	      ++i;
	    }
	  if (i == (lzmax + 1))
	    return Start;
	}
      ++TmpKey2;
    }
  return this->HilbertSpaceDimension;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& TrappedBosons::PrintState (ostream& Str, int state)
{
  int* TmpState = this->StateDescription[state];
  int Max = this->StateLzMax[state];
  int i = 0;
  for (; i <= Max; ++i)
    Str << TmpState[i] << " ";
  for (; i <= this->TotalLz; ++i)
    Str << "0 ";
  Str << " key = " << this->Keys[state] << " lzmax position = " << this->LzMaxPosition[Max * (this->NbrBosons + 1) + TmpState[Max]]
      << " position = " << FindStateIndex(TmpState, Max);
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson in the state
// currentLzMax = momentum maximum value for bosons that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int TrappedBosons::GenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int pos)
{
  if ((nbrBosons == 0) || ((nbrBosons * currentLzMax) < totalLz) || (pos == this->HilbertSpaceDimension))
    {
//      cout << "reject "<< nbrBosons << " " << lzMax << " " << currentLzMax << " " << totalLz << " " << pos << " " << endl;
      return pos;
    }
  if ((nbrBosons * currentLzMax) == totalLz)
    {
//      cout << "check ! " << nbrBosons << " " << lzMax << " " << currentLzMax << " " << totalLz << " " << pos << " " << endl;
      this->StateDescription[pos] = new int [lzMax + 1];
      int* TmpState = this->StateDescription[pos];
      for (int i = 0; i <= lzMax; ++i)
	TmpState[i] = 0;
      TmpState[currentLzMax] = nbrBosons;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }
  if ((currentLzMax == 0) || (totalLz == 0))
    {
//      cout << "check " << nbrBosons << " " << lzMax << " " << currentLzMax << " " << totalLz << " " << pos << " " << endl;
      this->StateDescription[pos] = new int [lzMax + 1];
      int* TmpState = this->StateDescription[pos];
      for (int i = 1; i <= lzMax; ++i)
	TmpState[i] = 0;
      TmpState[0] = nbrBosons;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }
//  cout << "normal " << nbrBosons << " " << lzMax << " " << currentLzMax << " " << totalLz << " " << pos << " " << endl;

  int TmpTotalLz = totalLz / currentLzMax;
  int TmpNbrBosons = nbrBosons - TmpTotalLz;
  TmpTotalLz = totalLz - TmpTotalLz * currentLzMax;
  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = pos;
  while (TmpNbrBosons < nbrBosons)
    {
      TmpPos = this->GenerateStates(TmpNbrBosons, lzMax, ReducedCurrentLzMax, TmpTotalLz, pos);
      for (int i = pos; i < TmpPos; i++)
	this->StateDescription[i][currentLzMax] = nbrBosons - TmpNbrBosons;
      ++TmpNbrBosons;
      pos = TmpPos;
      TmpTotalLz += currentLzMax;
    }
  if (lzMax == currentLzMax)
    return this->GenerateStates(nbrBosons, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, pos);
  else
    return this->GenerateStates(nbrBosons, lzMax, ReducedCurrentLzMax, totalLz, pos);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void TrappedBosons::GenerateLookUpTable(int memeory)
{
  this->Keys = new int [this->HilbertSpaceDimension];
  for (int i = 0; i <= this->TotalLz; ++i)
    this->KeyMultiplicationTable[i] = i * i * this->IncNbrBosons;

  this->LzMaxPosition = new int [(this->TotalLz + 1) * this->IncNbrBosons];
  int CurrentLzMax = this->TotalLz;
  int CurrentNbrLzMax = 1;
  this->LzMaxPosition[CurrentLzMax * (this->NbrBosons +1) + CurrentNbrLzMax] = 0; 
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->Keys[i] = this->GenerateKey(this->StateDescription[i], this->StateLzMax[i]);
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  CurrentLzMax = this->StateLzMax[i];
	  CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	  this->LzMaxPosition[CurrentLzMax * this->IncNbrBosons + CurrentNbrLzMax] = i; 
	}
      else
	if (this->StateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	    this->LzMaxPosition[CurrentLzMax * this->IncNbrBosons + CurrentNbrLzMax] = i; 
	  }
    }
}

// generate look-up table associated to current Hilbert space
// 
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = key associated to the state

int TrappedBosons::GenerateKey(int* stateDescription, int lzmax)
{
  int Key = 0;
  for (int i = 0; i <= lzmax; ++i)
    {
      Key += this->KeyMultiplicationTable[i] * stateDescription[i];
    }
  return Key;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

int TrappedBosons::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons == 0) || ((nbrBosons * lzMax) < totalLz))
    return 0;
  if (((nbrBosons * lzMax) == totalLz) || (lzMax == 0) || (totalLz == 0))
    return 1;
  int TmpDim = 0;
  while (totalLz >= 0)
    {
      TmpDim += this->EvaluateHilbertSpaceDimension(nbrBosons, lzMax - 1, totalLz);
      --nbrBosons;
      totalLz -= lzMax;
    }
  return TmpDim;
}

void TrappedBosons::ErasthothenesSlieve(int maxNumber, int nbrWantedPrime, int factor)
{
   maxNumber &= ~((int) 0x1);
  int Size = (int) sqrt(maxNumber);
  int NbrPrime = 1;
  int* PrimeNumbers = new int [Size];
  int NbrKeyFactor = 1;
  this->KeyMultiplicationTable[0] = 0;
  PrimeNumbers[0] = 3;
  int Lim;
  int Lim2 = factor;
  int j;  
  for (int i = 5; i < maxNumber; i +=2)
    {
      Lim = (int) sqrt((double) i);
      for (j = 0; j < NbrPrime; ++j)
	if ((j > Lim) || ((i % PrimeNumbers[j]) == 0))
	  {
	    j = NbrPrime + 1;
	  }
      if (j == NbrPrime)
	{
	  if(i < Size)
	    {
	      PrimeNumbers[NbrPrime] = i;
	      ++NbrPrime;
	    }
	  if (i > Lim2)
	    {
	      this->KeyMultiplicationTable[NbrKeyFactor] = i;
	      ++NbrKeyFactor;
	      Lim2 += i * factor;
	      if (NbrKeyFactor == nbrWantedPrime)
		i = maxNumber;
	    }
	}
    }
  for (int i = 0; i <= nbrWantedPrime; ++i)
    cout << this->KeyMultiplicationTable[i] << " ";
  cout << endl;
  delete[] PrimeNumbers;
}
