////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of bosons on torus                         //
//                                                                            //
//                        last modification : 03/09/2002                      //
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
#include "HilbertSpace/BosonOnTorus.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"

#include <math.h>


// basic constructor
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson

BosonOnTorus::BosonOnTorus (int nbrBosons, int maxMomentum)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;

  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->MaxMomentum);
  cout << "dim = " << this->HilbertSpaceDimension << endl;
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
  cout << this->GenerateStates(this->NbrBosons, this->MaxMomentum - 1, this->MaxMomentum - 1, 0) << endl;
  this->KeyMultiplicationTable = new int [this->MaxMomentum + 1];
  this->GenerateLookUpTable(0);
#ifdef __DEBUG__
  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += (this->StateMaxMomentum[i] + 1) * sizeof(int) + sizeof(int*) + sizeof(int);
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

BosonOnTorus::BosonOnTorus(const BosonOnTorus& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrLzValue = bosons.NbrLzValue;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;
  this->Flag = bosons.Flag;
}

// destructor
//

BosonOnTorus::~BosonOnTorus ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorus& BosonOnTorus::operator = (const BosonOnTorus& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
    }
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrLzValue = bosons.NbrLzValue;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;
  this->Flag = bosons.Flag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnTorus::Clone()
{
  return new BosonOnTorus(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnTorus::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (0);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnTorus::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (0);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnTorus::ExtractSubspace (AbstractQuantumNumber& q, 
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

int BosonOnTorus::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int CurrentLzMax = this->StateMaxMomentum[index];
  int* State = this->StateDescription[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || (State[n1] == 0) || (State[n2] == 0) || ((n1 == n2) && (State[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = CurrentLzMax;
  if (NewLzMax < m1)
    NewLzMax = m1;
  if (NewLzMax < m2)
    NewLzMax = m2;
  int i = 0;
  int* TemporaryState = new int [NewLzMax + 1];
  for (; i <= CurrentLzMax; ++i)
    TemporaryState[i] = State[i];
  for (; i <= NewLzMax; ++i)
    TemporaryState[i] = 0;
  coefficient = TemporaryState[n2];
  --TemporaryState[n2];
  coefficient *= TemporaryState[n1];
  --TemporaryState[n1];
  ++TemporaryState[m2];
  coefficient *= TemporaryState[m2];
  ++TemporaryState[m1];
  coefficient *= TemporaryState[m1];
  coefficient = sqrt(coefficient);
  while (TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  int DestIndex = this->FindStateIndex(TemporaryState, NewLzMax);
  delete[] TemporaryState;
  return DestIndex;
}

// return matrix representation of the annihilation operator a_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& BosonOnTorus::A (int i, Matrix& M)
{
  return M;
}

// return matrix representation of the creation operator a^+_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& BosonOnTorus::Ad (int i, Matrix& M)
{
  return M;
}

// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnTorus::FindStateIndex(int* stateDescription, int lzmax)
{
  int TmpKey = this->GenerateKey(stateDescription, lzmax);
  int i;
  int* TmpStateDescription;
  int Start = this->MomentumMaxPosition[lzmax * this->IncNbrBosons + stateDescription[lzmax]];
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

ostream& BosonOnTorus::PrintState (ostream& Str, int state)
{
  int* TmpState = this->StateDescription[state];
  int Max = this->StateMaxMomentum[state];
  int i = 0;
  for (; i <= Max; ++i)
    Str << TmpState[i] << " ";
  for (; i < this->MaxMomentum; ++i)
    Str << "0 ";
  Str << " key = " << this->Keys[state] << " lzmax position = " << this->MomentumMaxPosition[Max * (this->NbrBosons + 1) + TmpState[Max]]
      << " position = " << FindStateIndex(TmpState, Max);
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int BosonOnTorus::GenerateStates(int nbrBosons, int maxMomentum, int currentMaxMomentum, int pos)
{
//  cout << nbrBosons << " " << lzMax << " " << currentLzMax << " " << totalLz << " " << pos << endl;
  if ((nbrBosons == 0) || (currentMaxMomentum < 0))
    {
      return pos;
    }
  if (currentMaxMomentum == 0)
    {
      this->StateDescription[pos] = new int [maxMomentum + 1];
      int* TmpState = this->StateDescription[pos];
      for (int i = 1; i <= maxMomentum; ++i)
	TmpState[i] = 0;
      TmpState[0] = nbrBosons;
      this->StateMaxMomentum[pos] = maxMomentum;
      return pos + 1;
    }

  this->StateDescription[pos] = new int [maxMomentum + 1];
  int* TmpState = this->StateDescription[pos];
  for (int i = 0; i <= maxMomentum; ++i)
    TmpState[i] = 0;
  TmpState[currentMaxMomentum] = nbrBosons;
  this->StateMaxMomentum[pos] = maxMomentum;
  ++pos;

  int TmpNbrBosons = 1;
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int TmpPos = pos;
  while (TmpNbrBosons < nbrBosons)
    {
      TmpPos = this->GenerateStates(TmpNbrBosons, maxMomentum, ReducedCurrentMaxMomentum, pos);
      for (int i = pos; i < TmpPos; i++)
	this->StateDescription[i][currentMaxMomentum] = nbrBosons - TmpNbrBosons;
      ++TmpNbrBosons;
      pos = TmpPos;
    }
  if (maxMomentum == currentMaxMomentum)
    return this->GenerateStates(nbrBosons, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, pos);
  else
    return this->GenerateStates(nbrBosons, maxMomentum, ReducedCurrentMaxMomentum, pos);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnTorus::GenerateLookUpTable(int memory)
{
  this->Keys = new int [this->HilbertSpaceDimension];
  for (int i = 0; i <= this->MaxMomentum; ++i)
    this->KeyMultiplicationTable[i] = i * i * this->IncNbrBosons;

  this->MomentumMaxPosition = new int [(this->MaxMomentum + 2) * this->IncNbrBosons];
  int CurrentMaxMomentum = this->MaxMomentum;
  int CurrentNbrMaxMomentum = 1;
  this->MomentumMaxPosition[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum] = 0; 
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->Keys[i] = this->GenerateKey(this->StateDescription[i], this->StateMaxMomentum[i]);
      if (CurrentMaxMomentum != this->StateMaxMomentum[i])
	{
	  CurrentMaxMomentum = this->StateMaxMomentum[i];
	  CurrentNbrMaxMomentum = this->StateDescription[i][CurrentMaxMomentum];
	  this->MomentumMaxPosition[CurrentMaxMomentum * this->IncNbrBosons + CurrentNbrMaxMomentum] = i; 
	}
      else
	if (this->StateDescription[i][CurrentMaxMomentum] != CurrentNbrMaxMomentum)
	  {
	    CurrentNbrMaxMomentum = this->StateDescription[i][CurrentMaxMomentum];
	    this->MomentumMaxPosition[CurrentMaxMomentum * this->IncNbrBosons + CurrentNbrMaxMomentum] = i; 
	  }
    }
}

// generate look-up table associated to current Hilbert space
// 
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = key associated to the state

int BosonOnTorus::GenerateKey(int* stateDescription, int lzmax)
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
// maxMomentum = momentum maximum value for a boson
// return value = Hilbert space dimension

int BosonOnTorus::EvaluateHilbertSpaceDimension(int nbrBosons, int maxMomentum)
{
  FactorialCoefficient Dimension; 
  Dimension.PartialFactorialMultiply(maxMomentum, maxMomentum + nbrBosons - 1); 
  Dimension.FactorialDivide(nbrBosons);
  return Dimension.GetIntegerValue();
}
