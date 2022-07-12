////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of bosons in Aharonov Bohm cages                  //
//                                                                            //
//                        last modification : 09/07/2002                      //
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
#include "HilbertSpace/ManyBodyHilbertSpace/AharonovBohmCages.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"

#include <math.h>


using std::cout;
using std::endl;


// basic constructor
// 
// nbrBosons = number of bosons
// nbrCages = number of cages

AharonovBohmCages::AharonovBohmCages (int nbrBosons, int nbrCages)
{
  this->NbrBosons = nbrBosons;
  this->NbrCages = nbrCages;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrCages);
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateMaxPosition = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrBosons, this->NbrCages, this->NbrCages, 0);
  this->KeyMultiplicationTable = new int [this->NbrCages + 1];
//  this->ErasthothenesSlieve(10000000, this->NbrCages + 1, this->IncNbrBosons);
  this->GenerateLookUpTable(0);
#ifdef __DEBUG__
  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += (this->StateMaxPosition[i] + 1) * sizeof(int) + sizeof(int*);
  UsedMemory += (this->NbrCages + 1) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += ((this->NbrCages + 1) * this->IncNbrBosons) * sizeof(int);
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
// space = reference on the hilbert space to copy to copy

AharonovBohmCages::AharonovBohmCages(const AharonovBohmCages& space)
{
  this->NbrBosons = space.NbrBosons;
  this->IncNbrBosons = space.IncNbrBosons;
  this->NbrCages = space.NbrCages;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;
  this->StateDescription = space.StateDescription;
  this->StateMaxPosition = space.StateMaxPosition;
  this->Flag = space.Flag;
}

// destructor
//

AharonovBohmCages::~AharonovBohmCages ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateMaxPosition;
    }
}

// assignement (without duplicating datas)
//
// space = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

AharonovBohmCages& AharonovBohmCages::operator = (const AharonovBohmCages& space)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateMaxPosition;
    }
  this->NbrBosons = space.NbrBosons;
  this->IncNbrBosons = space.IncNbrBosons;
  this->NbrCages = space.NbrCages;
  this->HilbertSpaceDimension = space.HilbertSpaceDimension;
  this->StateDescription = space.StateDescription;
  this->StateMaxPosition = space.StateMaxPosition;
  this->Flag = space.Flag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* AharonovBohmCages::Clone()
{
  return new AharonovBohmCages(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> AharonovBohmCages::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->NbrCages);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* AharonovBohmCages::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->NbrCages);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* AharonovBohmCages::ExtractSubspace (AbstractQuantumNumber& q, 
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

int AharonovBohmCages::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int LzMax = this->StateMaxPosition[index];
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

int AharonovBohmCages::FindStateIndex(int* stateDescription, int lzmax)
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

ostream& AharonovBohmCages::PrintState (ostream& Str, int state)
{
  int* TmpState = this->StateDescription[state];
  int Max = this->StateMaxPosition[state];
  int i = 0;
  for (; i <= Max; ++i)
    Str << TmpState[i] << " ";
  for (; i <= this->NbrCages; ++i)
    Str << "0 ";
  Str << " key = " << this->Keys[state] << " lzmax position = " << this->LzMaxPosition[Max * (this->NbrBosons + 1) + TmpState[Max]]
      << " position = " << FindStateIndex(TmpState, Max);
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson in the state
// currentCage = position of the current cage that has to be filled
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int AharonovBohmCages::GenerateStates(int nbrBosons, int currentMaxCage, int currentCage, int pos)
{
  if ((nbrBosons == 0) || (currentCage < 0) || (pos == this->HilbertSpaceDimension))
    {
//      cout << "reject "<< nbrBosons << " " << lzMax << " " << currentLzMax << " " << totalLz << " " << pos << " " << endl;
      return pos;
    }
  int TmpPos = pos;
  if (currentMaxCage == currentCage)
    TmpPos += GenerateStates(nbrBosons, currentCage - 1, currentCage - 1, TmpPos);
  else
    TmpPos += GenerateStates(nbrBosons, currentMaxCage, currentCage - 1, TmpPos);    
  return this->GenerateStates(nbrBosons, currentMaxCage, currentCage - 1, pos);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void AharonovBohmCages::GenerateLookUpTable(int memeory)
{
  this->Keys = new int [this->HilbertSpaceDimension];
  for (int i = 0; i <= this->NbrCages; ++i)
    this->KeyMultiplicationTable[i] = i * i * this->IncNbrBosons;

  this->LzMaxPosition = new int [(this->NbrCages + 1) * this->IncNbrBosons];
  int CurrentLzMax = this->NbrCages;
  int CurrentNbrLzMax = 1;
  this->LzMaxPosition[CurrentLzMax * (this->NbrBosons +1) + CurrentNbrLzMax] = 0; 
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->Keys[i] = this->GenerateKey(this->StateDescription[i], this->StateMaxPosition[i]);
      if (CurrentLzMax != this->StateMaxPosition[i])
	{
	  CurrentLzMax = this->StateMaxPosition[i];
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

int AharonovBohmCages::GenerateKey(int* stateDescription, int lzmax)
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
// nbrCages = number of cages
// return value = Hilbert space dimension

int AharonovBohmCages::EvaluateHilbertSpaceDimension(int nbrBosons, int nbrCages)
{
  int TmpDim = 3 * nbrCages;
  for (int i = 1; i < nbrBosons; ++i)
    TmpDim *= 3 * nbrCages;
  return TmpDim;
}

