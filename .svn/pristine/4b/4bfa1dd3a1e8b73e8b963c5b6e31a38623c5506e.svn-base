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
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include <math.h>


using std::cout;
using std::endl;


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
  this->MomentumConstraintFlag = false;

  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->MaxMomentum);
  cout << "dim = " << this->HilbertSpaceDimension << endl;
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
  cout << this->GenerateStates(this->NbrBosons, this->MaxMomentum - 1, this->MaxMomentum - 1, 0) << endl;
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

// constructor with a constraint of the total momentum of states
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// momentumConstraint = index of the momentum orbit

BosonOnTorus::BosonOnTorus (int nbrBosons, int maxMomentum, int momentumConstraint)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->MomentumConstraint = momentumConstraint;
  this->MomentumConstraintFlag = true;

  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->MaxMomentum);
  cout << "dim = " << this->HilbertSpaceDimension << endl;
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->MaxMomentum - 1, this->MaxMomentum - 1, 0, 0);
  this->GenerateLookUpTable(1000000);
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

// constructor from full datas (with no constraint on the total momentum)
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// hilbertSpaceDimension = Hilbert space dimension
// stateDescription = array describing each state
// stateMaxMomentum = array giving maximum Lz value reached for a fermion in a given state

BosonOnTorus::BosonOnTorus (int nbrBosons, int maxMomentum, int hilbertSpaceDimension, 
			    int** stateDescription, int* stateMaxMomentum)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->MomentumConstraintFlag = false;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = stateDescription;
  this->StateMaxMomentum = stateMaxMomentum;
  this->GenerateLookUpTable(1000000);
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

// constructor from full datas
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// momentumConstraint = index of the momentum orbit
// hilbertSpaceDimension = Hilbert space dimension
// stateDescription = array describing each state
// stateMaxMomentum = array giving maximum Lz value reached for a boson in a given state

BosonOnTorus::BosonOnTorus (int nbrBosons, int maxMomentum, int momentumConstraint, int hilbertSpaceDimension, 
			    int** stateDescription, int* stateMaxMomentum)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->MomentumConstraint = momentumConstraint;
  this->MomentumConstraintFlag = true;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = stateDescription;
  this->StateMaxMomentum = stateMaxMomentum;
  this->GenerateLookUpTable(1000000);
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
      delete[] this->Keys;
      delete[] this->KeyMultiplicationTable;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
      int Size = (this->MaxMomentum + 2) * this->IncNbrBosons;
      for (int i = 0; i < Size; ++i)
	{
	  if (this->KeyInvertSectorSize[i] > 0)
	    {
	      for (int j= 0; j < this->KeyInvertSectorSize[i]; ++j)
		delete[] this->KeyInvertIndices[i][j];
	      delete[] this->KeyInvertTable[i];
	      delete[] this->KeyInvertTableNbrIndices[i];
	      delete[] this->KeyInvertIndices[i];
	    }
	}
      
      delete[] this->KeyInvertSectorSize;
      delete[] this->KeyInvertTable;
      delete[] this->KeyInvertTableNbrIndices;
      delete[] this->KeyInvertIndices;
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

// get momemtum value of a given state
//
// index = state index
// return value = state momentum

int BosonOnTorus::GetMomentumValue(int index)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  int Momentum = 0;
  int* TmpStateDescription = this->StateDescription[index];
  for (int i = 0; i <= StateMaxMomentum; ++i)
    {
      Momentum += (TmpStateDescription[i] * i);
    }
  return (Momentum % this->MaxMomentum);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnTorus::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->MomentumConstraintFlag == false)
    {
      for (int i = 0; i < this->MaxMomentum; ++i)
	TmpList += new PeriodicMomentumQuantumNumber (i, this->MaxMomentum);
    }
  else
    TmpList += new PeriodicMomentumQuantumNumber (this->MomentumConstraint, this->MaxMomentum);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnTorus::GetQuantumNumber (int index)
{
  if (this->MomentumConstraintFlag == false)
    {
      return  new PeriodicMomentumQuantumNumber (this->GetMomentumValue(index), this->MaxMomentum);
    }
  else
    return new PeriodicMomentumQuantumNumber (this->MomentumConstraint, this->MaxMomentum);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnTorus::ExtractSubspace (AbstractQuantumNumber& q, 
						     SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != AbstractQuantumNumber::PeriodicMomentum)
    return 0;
  if (this->MomentumConstraintFlag == true)
    if (this->MomentumConstraint == ((PeriodicMomentumQuantumNumber&) q).GetMomentum())
      return this;
    else 
      return 0;
  int Momentum = ((PeriodicMomentumQuantumNumber&) q).GetMomentum();
  int* TmpConvArray = new int [this->HilbertSpaceDimension];
  int SubspaceHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (this->GetMomentumValue(i) == Momentum)
	{
	  TmpConvArray[SubspaceHilbertSpaceDimension] = i;
	  ++SubspaceHilbertSpaceDimension;
	}
    }
  int** SubspaceStateDescription = new int* [SubspaceHilbertSpaceDimension];
  int* SubspaceStateMaxMomentum = new int [SubspaceHilbertSpaceDimension];
  int* ConvArray = new int [SubspaceHilbertSpaceDimension];
  for (int i = 0; i < SubspaceHilbertSpaceDimension; ++i)
    {
      ConvArray[i] = TmpConvArray[i];
      SubspaceStateMaxMomentum[i] = this->StateMaxMomentum[TmpConvArray[i]];
      SubspaceStateDescription[i] = new int [SubspaceStateMaxMomentum[i] + 1];
      for (int j = 0; j <= SubspaceStateMaxMomentum[i]; ++j)
	SubspaceStateDescription[i][j] = this->StateDescription[TmpConvArray[i]][j];
    }
  delete[] TmpConvArray;
  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, SubspaceHilbertSpaceDimension, ConvArray);
  return new BosonOnTorus(this->NbrBosons, this->MaxMomentum, Momentum, SubspaceHilbertSpaceDimension, 
			  SubspaceStateDescription, SubspaceStateMaxMomentum);  
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
/*  int TmpKey = this->GenerateKey(stateDescription, lzmax);
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
  return this->HilbertSpaceDimension;*/
  int TmpKey = this->GenerateKey(stateDescription, lzmax);
  int Sector = lzmax * this->IncNbrBosons + stateDescription[lzmax];
  int TmpPos = 0;
  int TmpPos2 = this->KeyInvertSectorSize[Sector] - 1;
  int TmpPos3;
  int* TmpKeyInvertTable = this->KeyInvertTable[Sector];
  while (TmpPos2 != TmpPos)
    {
      TmpPos3 = (TmpPos2 + TmpPos) >> 1;
      if (TmpKey < TmpKeyInvertTable[TmpPos3])
	{
	  TmpPos2 = TmpPos3 - 1;
	}
      else
	if (TmpKey > TmpKeyInvertTable[TmpPos3])
	  {
	    TmpPos = TmpPos3 + 1;
	  }
	else
	  {
	    TmpPos2 = TmpPos3;
	    TmpPos = TmpPos3;		    
	  }
    }
  int i;
  int* TmpStateDescription;
  int Start;
  TmpPos3 = this->KeyInvertTableNbrIndices[Sector][TmpPos];
  int* TmpKeyInvertIndices = this->KeyInvertIndices[Sector][TmpPos];
  for (int j = 0; j < TmpPos3; ++j)
    {
      Start = TmpKeyInvertIndices[j];
      i = 0;
      TmpStateDescription = this->StateDescription[Start];
      while (i <= lzmax)
	{
	  if (stateDescription[i] != TmpStateDescription[i])
	    i = lzmax + 1;
	  ++i;
	}
      if (i == (lzmax + 1))
	return Start;
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
  Str << " key = " << this->Keys[state] << " max momentum = " << this->MomentumMaxPosition[Max * (this->NbrBosons + 1) + TmpState[Max]]
      << " position = " << FindStateIndex(TmpState, Max) << " momentum = " << this->GetMomentumValue(state);;
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

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson in the state
// currentMaxMomentum = momentum maximum value for bosons that are still to be placed
// pos = position in StateDescription array where to store states
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

int BosonOnTorus::GenerateStates(int nbrBosons, int maxMomentum, int currentMaxMomentum, int pos, int currentMomentum)
{
  if (nbrBosons == 0)
    {
      if ((currentMomentum % this->MaxMomentum) == this->MomentumConstraint)
	{
	  this->StateDescription[pos] = new int [maxMomentum + 1];
	  int* TmpState = this->StateDescription[pos];
	  for (int i = 0; i <= maxMomentum; ++i)
	    TmpState[i] = 0;
	  this->StateMaxMomentum[pos] = maxMomentum;
	  return pos + 1;
	}
      else
	{
	  return pos;
	}
    }
  if (currentMaxMomentum == 0)
    {
      if ((currentMomentum % this->MaxMomentum) == this->MomentumConstraint)
	{
	  this->StateDescription[pos] = new int [maxMomentum + 1];
	  int* TmpState = this->StateDescription[pos];
	  for (int i = 1; i <= maxMomentum; ++i)
	    TmpState[i] = 0;
	  TmpState[0] = nbrBosons;
	  this->StateMaxMomentum[pos] = maxMomentum;
	  return pos + 1;
	}
      else
	return pos;
    }

  int TmpNbrBosons = 0;
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int TmpPos = pos;
  while (TmpNbrBosons < nbrBosons)
    {
      TmpPos = this->GenerateStates(TmpNbrBosons, maxMomentum, ReducedCurrentMaxMomentum, pos, currentMomentum + (nbrBosons - TmpNbrBosons) * currentMaxMomentum);
      for (int i = pos; i < TmpPos; i++)
	this->StateDescription[i][currentMaxMomentum] = nbrBosons - TmpNbrBosons;
      ++TmpNbrBosons;
      pos = TmpPos;
    }
  if (maxMomentum == currentMaxMomentum)
    return this->GenerateStates(nbrBosons, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, pos, currentMomentum);
  else
    return this->GenerateStates(nbrBosons, maxMomentum, ReducedCurrentMaxMomentum, pos, currentMomentum);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnTorus::GenerateLookUpTable(int memory)
{
  this->Keys = new int [this->HilbertSpaceDimension];
  this->KeyMultiplicationTable = new int [this->MaxMomentum + 1];
  for (int i = 0; i <= this->MaxMomentum; ++i)
    this->KeyMultiplicationTable[i] = i * i * this->IncNbrBosons;

  int Size = (this->MaxMomentum + 2) * this->IncNbrBosons;
  this->MomentumMaxPosition = new int [Size];
  this->KeyInvertSectorSize = new int [Size];
  this->KeyInvertTable = new int* [Size];
  this->KeyInvertTableNbrIndices = new int* [Size];
  this->KeyInvertIndices = new int** [Size];
  for (int i = 0; i < Size; ++i)
    this->KeyInvertSectorSize[i] =0;
  int CurrentMaxMomentum = this->MaxMomentum;
  int CurrentNbrMaxMomentum = 1;
  this->MomentumMaxPosition[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum] = 0; 
  this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum] = 0;   
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->Keys[i] = this->GenerateKey(this->StateDescription[i], this->StateMaxMomentum[i]);
      if (CurrentMaxMomentum != this->StateMaxMomentum[i])
	{
	  CurrentMaxMomentum = this->StateMaxMomentum[i];
	  CurrentNbrMaxMomentum = this->StateDescription[i][CurrentMaxMomentum];
	  this->MomentumMaxPosition[CurrentMaxMomentum * this->IncNbrBosons + CurrentNbrMaxMomentum] = i; 
	  this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum] = 1; 
	}
      else
	if (this->StateDescription[i][CurrentMaxMomentum] != CurrentNbrMaxMomentum)
	  {
	    CurrentNbrMaxMomentum = this->StateDescription[i][CurrentMaxMomentum];
	    this->MomentumMaxPosition[CurrentMaxMomentum * this->IncNbrBosons + CurrentNbrMaxMomentum] = i; 
	    this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum] = 1; 
	  }
	else
	  {
	    ++this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum]; 
	  }
    }

  for (int i = 0; i < Size; ++i)
    {
      if (this->KeyInvertSectorSize[i] > 0)
	{
	  this->KeyInvertTable[i] = new int [this->KeyInvertSectorSize[i]];
	  this->KeyInvertTableNbrIndices[i] = new int [this->KeyInvertSectorSize[i]];
	}
      else
	{
	  this->KeyInvertTable[i] = 0; 
	  this->KeyInvertTableNbrIndices[i] = 0;
	}
    }

  CurrentMaxMomentum = this->StateMaxMomentum[0];
  CurrentNbrMaxMomentum = this->StateDescription[0][CurrentMaxMomentum];
  int CurrentKeyInvertSectorSize = this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
  this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum] = 1;
  int* TmpKeyInvertTable = this->KeyInvertTable[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
  int* TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
  TmpKeyInvertTable[0] = this->Keys[0];
  TmpKeyInvertTableNbrIndices[0] = 1;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentMaxMomentum != this->StateMaxMomentum[i])
	{
	  cout << "sector " << CurrentNbrMaxMomentum << "/" << CurrentMaxMomentum << ": " << CurrentKeyInvertSectorSize  
	       << " " << this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum]<< endl;
	  CurrentMaxMomentum = this->StateMaxMomentum[i];
	  CurrentNbrMaxMomentum = this->StateDescription[i][CurrentMaxMomentum];
	  CurrentKeyInvertSectorSize = this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
	  this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum] = 1;
	  TmpKeyInvertTable = this->KeyInvertTable[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
	  TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
	  TmpKeyInvertTable[0] = this->Keys[i];
	  TmpKeyInvertTableNbrIndices[0] = 1;
	}
      else
	if (this->StateDescription[i][CurrentMaxMomentum] != CurrentNbrMaxMomentum)
	  {
	    cout << "sector " << CurrentNbrMaxMomentum << "/" << CurrentMaxMomentum << ": " <<  CurrentKeyInvertSectorSize
		 << " " << this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum] << endl;
	    CurrentNbrMaxMomentum = this->StateDescription[i][CurrentMaxMomentum];
	    CurrentKeyInvertSectorSize = this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
	    this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum] = 1;
	    TmpKeyInvertTable = this->KeyInvertTable[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
	    TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
	    TmpKeyInvertTable[0] = this->Keys[i];
	    TmpKeyInvertTableNbrIndices[0] = 1;
	  }
	else
	  {
	    int Lim = this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
	    int j = 0;
	    bool Flag = false;
	    int TmpKey = this->Keys[i];
	    for (; ((j < Lim) && (Flag == false)); ++j)
	      if (TmpKeyInvertTable[j] == TmpKey)
		{
		  Flag = true;
		  --j;
		}
	    if (Flag == true)
	      {
		++TmpKeyInvertTableNbrIndices[j];
	      }
	    else
	      {
		TmpKeyInvertTable[this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum]] = TmpKey;
		TmpKeyInvertTableNbrIndices[this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum]] = 1;
		++this->KeyInvertSectorSize[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
	      }
	  }
    }
 
  // sort key invert table and resize arrays
  int TmpPos;
  int TmpSize;
  for (int i = 0; i < Size; ++i)
    {
      if (this->KeyInvertSectorSize[i] > 0)
	{
	  int Lim = this->KeyInvertSectorSize[i];
	  TmpKeyInvertTable = this->KeyInvertTable[i];
	  TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[i];
	  int* TmpKeyInvertTable2 = new int [Lim];
	  int* TmpKeyInvertTableNbrIndices2 = new int [Lim];
	  int Tmp;
/*	  cout << "sector size = " << Lim << endl;
	  cout << "before sort" << endl;
	  for (int j = 0; j < Lim; ++j)
	    {
	      cout << TmpKeyInvertTable[j] << " " <<  TmpKeyInvertTableNbrIndices[j] << endl;
	    }*/
	  for (int j = 0; j < Lim; ++j)
	    {
	      Tmp = TmpKeyInvertTable[j];
	      TmpPos = j;
	      for (int k = j + 1; k < Lim; ++k)
		if (Tmp > TmpKeyInvertTable[k])
		  {
		    Tmp = TmpKeyInvertTable[k];
		    TmpPos = k;
		  }
	      TmpKeyInvertTable2[j] = Tmp;
	      TmpKeyInvertTableNbrIndices2[j] = TmpKeyInvertTableNbrIndices[TmpPos];
	      TmpKeyInvertTableNbrIndices[TmpPos] = TmpKeyInvertTableNbrIndices[j];
	      TmpKeyInvertTable[TmpPos] = TmpKeyInvertTable[j];
	    }
	  delete[] TmpKeyInvertTable;
	  delete[] TmpKeyInvertTableNbrIndices;
	  this->KeyInvertTable[i] = TmpKeyInvertTable2;
	  this->KeyInvertTableNbrIndices[i] = TmpKeyInvertTableNbrIndices2;
/*	  cout << "after sort" << endl;
	  for (int j = 0; j < Lim; ++j)
	    {
	      cout << TmpKeyInvertTable2[j] << " " <<  TmpKeyInvertTableNbrIndices2[j] << endl;
	    }*/
	  this->KeyInvertIndices[i] = new int* [Lim];
	  for (int j = 0; j < Lim; ++j)
	    {
	      this->KeyInvertIndices[i][j] = new int [this->KeyInvertTableNbrIndices[i][j]];
	      this->KeyInvertTableNbrIndices[i][j] = 0;
	    }
	}
    }

  // find all hilbert space indices that have the same key in each sector
  CurrentMaxMomentum = this->MaxMomentum;
  CurrentNbrMaxMomentum = 1;
  TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
  int** TmpKeyInvertIndices = this->KeyInvertIndices[CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum];
  int TmpPos2;
  int TmpPos3;
  int TmpPos4 = CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum;
  TmpSize = this->KeyInvertSectorSize[TmpPos4];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      int TmpKey = this->Keys[i];
      if (CurrentMaxMomentum != this->StateMaxMomentum[i])
	{
	  CurrentMaxMomentum = this->StateMaxMomentum[i];
	  CurrentNbrMaxMomentum = this->StateDescription[i][CurrentMaxMomentum];
	  TmpPos4 = CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum;
	  TmpKeyInvertIndices = this->KeyInvertIndices[TmpPos4];
	  TmpKeyInvertTable = this->KeyInvertTable[TmpPos4];
	  TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[TmpPos4];
	  TmpSize = this->KeyInvertSectorSize[TmpPos4];
	  TmpPos = 0;
	  TmpPos2 = TmpSize - 1;
	  while (TmpPos2 != TmpPos)
	    {
	      TmpPos3 = (TmpPos2 + TmpPos) >> 1;
	      if (TmpKey < TmpKeyInvertTable[TmpPos3])
		{
		  TmpPos2 = TmpPos3 - 1;
		}
	      else
		if (TmpKey > TmpKeyInvertTable[TmpPos3])
		  {
		    TmpPos = TmpPos3 + 1;
		  }
		else
		  {
		    TmpPos2 = TmpPos3;
		    TmpPos = TmpPos3;		    
		  }
	    }
	  TmpKeyInvertIndices[TmpPos][0] = i;
	  TmpKeyInvertTableNbrIndices[TmpPos] = 1;
	}
      else
	if (this->StateDescription[i][CurrentMaxMomentum] != CurrentNbrMaxMomentum)
	  {
	    CurrentNbrMaxMomentum = this->StateDescription[i][CurrentMaxMomentum];
	    TmpPos4 = CurrentMaxMomentum * (this->IncNbrBosons) + CurrentNbrMaxMomentum;
	    TmpKeyInvertIndices = this->KeyInvertIndices[TmpPos4];
	    TmpKeyInvertTable = this->KeyInvertTable[TmpPos4];
	    TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[TmpPos4];
	    TmpSize = this->KeyInvertSectorSize[TmpPos4];
	    TmpPos = 0;
	    TmpPos2 = TmpSize - 1;
	    while (TmpPos2 != TmpPos)
	      {
		TmpPos3 = (TmpPos2 + TmpPos) >> 1;
		if (TmpKey < TmpKeyInvertTable[TmpPos3])
		  {
		    TmpPos2 = TmpPos3 - 1;
		  }
	      else
		if (TmpKey > TmpKeyInvertTable[TmpPos3])
		  {
		    TmpPos = TmpPos3 + 1;
		  }
		else
		  {
		    TmpPos2 = TmpPos3;
		    TmpPos = TmpPos3;		    
		  }
	      }
	    TmpKeyInvertIndices[TmpPos][0] = i;
	    TmpKeyInvertTableNbrIndices[TmpPos] = 1;
	  }
	else
	  {
	    TmpPos = 0;
	    TmpPos2 = TmpSize - 1;
	    while (TmpPos2 != TmpPos)
	      {
		TmpPos3 = (TmpPos2 + TmpPos) >> 1;
		if (TmpKey < TmpKeyInvertTable[TmpPos3])
		  {
		    TmpPos2 = TmpPos3 - 1;
		  }
	      else
		if (TmpKey > TmpKeyInvertTable[TmpPos3])
		  {
		    TmpPos = TmpPos3 + 1;
		  }
		else
		  {
		    TmpPos2 = TmpPos3;
		    TmpPos = TmpPos3;		    
		  }
	      }
	    TmpKeyInvertIndices[TmpPos][TmpKeyInvertTableNbrIndices[TmpPos]] = i;
	    ++TmpKeyInvertTableNbrIndices[TmpPos];
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

