////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of bosons on sphere                        //
//                                                                            //
//                        last modification : 05/07/2002                      //
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
#include "HilbertSpace/BosonOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"

#include <math.h>


// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson

BosonOnSphere::BosonOnSphere (int nbrBosons, int totalLz, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, this->TotalLz);
//  cout << "dim = " << this->HilbertSpaceDimension << endl;
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrBosons, this->LzMax, this->LzMax, this->ShiftedTotalLz, 0);
  this->KeyMultiplicationTable = new int [this->LzMax + 1];
  this->GenerateLookUpTable(0);

  this->RatioSearchTotalDim = 0.0;
  this->RatioSearchCoarseDeep = 0.0;
  this->SqrRatioSearchTotalDim = 0.0;
  this->SqrRatioSearchCoarseDeep = 0.0;
  this->NbrSearch = 0;

#ifdef __DEBUG__
  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += (this->StateLzMax[i] + 1) * sizeof(int) + sizeof(int*);
  UsedMemory += (this->TotalLz + 1) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += ((this->TotalLz + 1) * this->IncNbrBosons) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
/*  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;*/
#endif
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSphere::BosonOnSphere(const BosonOnSphere& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->Flag = bosons.Flag;

  this->RatioSearchTotalDim = 0.0;
  this->RatioSearchCoarseDeep = 0.0;
  this->SqrRatioSearchTotalDim = 0.0;
  this->SqrRatioSearchCoarseDeep = 0.0;
  this->NbrSearch = 0;
}

// destructor
//

BosonOnSphere::~BosonOnSphere ()
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

BosonOnSphere& BosonOnSphere::operator = (const BosonOnSphere& bosons)
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
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->Flag = bosons.Flag;

  this->RatioSearchTotalDim = 0.0;
  this->RatioSearchCoarseDeep = 0.0;
  this->SqrRatioSearchTotalDim = 0.0;
  this->SqrRatioSearchCoarseDeep = 0.0;
  this->NbrSearch = 0;

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphere::Clone()
{
  return new BosonOnSphere(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphere::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphere::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphere::ExtractSubspace (AbstractQuantumNumber& q, 
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

int BosonOnSphere::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int CurrentLzMax = this->StateLzMax[index];
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

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphere::AdA (int index, int m)
{
  if (this->StateLzMax[index] < m)  
    return 0.0;
  return (double) (this->StateDescription[index][m]);  
}

// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnSphere::FindStateIndex(int* stateDescription, int lzmax)
{
/*  int TmpKey = this->GenerateKey(stateDescription, lzmax);
  int i;
  int* TmpStateDescription;
  int Start = this->LzMaxPosition[lzmax * this->IncNbrBosons + stateDescription[lzmax]];
  int* TmpKey2 = &(this->Keys[Start]);
  int NbrCoarseSearch = 0;
  int NbrDeepSearch = 0;
  for (; Start < this->HilbertSpaceDimension; ++Start)    
    {
      ++NbrCoarseSearch;
      if ((*TmpKey2) == TmpKey)
	{
	  ++NbrDeepSearch;
	  i = 0;
	  TmpStateDescription = this->StateDescription[Start];
	  while (i <= lzmax)
	    {
	      if (stateDescription[i] != TmpStateDescription[i])
		i = lzmax + 2;
	      ++i;
	    }
	  if (i == (lzmax + 1))
	    {
	      ++this->NbrSearch;
	      this->RatioSearchTotalDim += ((double) NbrCoarseSearch) / ((double) HilbertSpaceDimension);
	      this->RatioSearchCoarseDeep += ((double) NbrDeepSearch) / ((double) NbrCoarseSearch);
	      this->SqrRatioSearchTotalDim += (((double) NbrCoarseSearch) / ((double) HilbertSpaceDimension)) * (((double) NbrCoarseSearch) / ((double) HilbertSpaceDimension));
	      this->SqrRatioSearchCoarseDeep += (((double) NbrDeepSearch) / ((double) NbrCoarseSearch)) * (((double) NbrDeepSearch) / ((double) NbrCoarseSearch));
	      
	      if (this->NbrSearch >= 100000)
		{
		  cout << "search/total dim = " << (this->RatioSearchTotalDim / this->NbrSearch) << " +/- "  
		       << sqrt (this->SqrRatioSearchTotalDim / this->NbrSearch - (this->RatioSearchTotalDim / this->NbrSearch) * (this->RatioSearchTotalDim / this->NbrSearch))<< endl;
		  cout << "deep/coarse =" << (this->RatioSearchCoarseDeep / this->NbrSearch) << " +/- "  
		       << sqrt (this->SqrRatioSearchCoarseDeep / this->NbrSearch - (this->RatioSearchCoarseDeep / this->NbrSearch) * (this->RatioSearchCoarseDeep / this->NbrSearch))<< endl;
		  this->RatioSearchTotalDim = 0.0;
		  this->RatioSearchCoarseDeep = 0.0;
		  this->SqrRatioSearchTotalDim = 0.0;
		  this->SqrRatioSearchCoarseDeep = 0.0;
		  this->NbrSearch = 0;
		}
	      return Start;
	    }
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

ostream& BosonOnSphere::PrintState (ostream& Str, int state)
{
  int* TmpState = this->StateDescription[state];
  int Max = this->StateLzMax[state];
  int i = 0;
  for (; i <= Max; ++i)
    Str << TmpState[i] << " ";
  for (; i <= this->LzMax; ++i)
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

int BosonOnSphere::GenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int pos)
{
//  cout << nbrBosons << " " << lzMax << " " << currentLzMax << " " << totalLz << " " << pos << endl;
  if ((nbrBosons == 0) || ((nbrBosons * currentLzMax) < totalLz) || (pos == this->HilbertSpaceDimension))
    {
      return pos;
    }
  if ((nbrBosons * currentLzMax) == totalLz)
    {
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
      this->StateDescription[pos] = new int [lzMax + 1];
      int* TmpState = this->StateDescription[pos];
      for (int i = 1; i <= lzMax; ++i)
	TmpState[i] = 0;
      TmpState[0] = nbrBosons;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }

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

void BosonOnSphere::GenerateLookUpTable(int memory)
{
  this->Keys = new int [this->HilbertSpaceDimension];
  for (int i = 0; i <= this->LzMax; ++i)
    this->KeyMultiplicationTable[i] = i * i * i * i;//this->IncNbrBosons;
  int Size = (this->LzMax + 2) * this->IncNbrBosons;
  this->LzMaxPosition = new int [Size];
  this->KeyInvertSectorSize = new int [Size];
  this->KeyInvertTable = new int* [Size];
  this->KeyInvertTableNbrIndices = new int* [Size];
  this->KeyInvertIndices = new int** [Size];
  for (int i = 0; i < Size; ++i)
    this->KeyInvertSectorSize[i] =0;
  int CurrentLzMax = this->LzMax;
  int CurrentNbrLzMax = 1;
  this->LzMaxPosition[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 0; 
  this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 0;   
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->Keys[i] = this->GenerateKey(this->StateDescription[i], this->StateLzMax[i]);
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  CurrentLzMax = this->StateLzMax[i];
	  CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	  this->LzMaxPosition[CurrentLzMax * this->IncNbrBosons + CurrentNbrLzMax] = i; 
	  this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	}
      else
	if (this->StateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	    this->LzMaxPosition[CurrentLzMax * this->IncNbrBosons + CurrentNbrLzMax] = i;
	    this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	  }
	else
	  {
	    ++this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax]; 
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

  CurrentLzMax = this->StateLzMax[0];
  CurrentNbrLzMax = this->StateDescription[0][CurrentLzMax];
  int CurrentKeyInvertSectorSize = this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
  int* TmpKeyInvertTable = this->KeyInvertTable[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  int* TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  TmpKeyInvertTable[0] = this->Keys[0];
  TmpKeyInvertTableNbrIndices[0] = 1;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  cout << "sector " << CurrentNbrLzMax << "/" << CurrentLzMax << ": " << CurrentKeyInvertSectorSize  
	       << " " << this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax]<< endl;
	  CurrentLzMax = this->StateLzMax[i];
	  CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	  CurrentKeyInvertSectorSize = this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	  TmpKeyInvertTable = this->KeyInvertTable[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  TmpKeyInvertTable[0] = this->Keys[i];
	  TmpKeyInvertTableNbrIndices[0] = 1;
	}
      else
	if (this->StateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    cout << "sector " << CurrentNbrLzMax << "/" << CurrentLzMax << ": " <<  CurrentKeyInvertSectorSize
		 << " " << this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] << endl;
	    CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	    CurrentKeyInvertSectorSize = this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	    this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	    TmpKeyInvertTable = this->KeyInvertTable[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	    TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	    TmpKeyInvertTable[0] = this->Keys[i];
	    TmpKeyInvertTableNbrIndices[0] = 1;
	  }
	else
	  {
	    int Lim = this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
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
		TmpKeyInvertTable[this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax]] = TmpKey;
		TmpKeyInvertTableNbrIndices[this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax]] = 1;
		++this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
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
  CurrentLzMax = this->LzMax;
  CurrentNbrLzMax = 1;
  TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  int** TmpKeyInvertIndices = this->KeyInvertIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  int TmpPos2;
  int TmpPos3;
  int TmpPos4 = CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax;
  TmpSize = this->KeyInvertSectorSize[TmpPos4];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      int TmpKey = this->Keys[i];
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  CurrentLzMax = this->StateLzMax[i];
	  CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	  TmpPos4 = CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax;
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
	if (this->StateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	    TmpPos4 = CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax;
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

int BosonOnSphere::GenerateKey(int* stateDescription, int lzmax)
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

int BosonOnSphere::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, (totalLz + lzMax * nbrBosons) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrBosons = number of bosons
// lzMax = two times momentum maximum value for a boson plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

int BosonOnSphere::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons == 0) || ((nbrBosons * lzMax) < totalLz))
    return 0;
  if (((nbrBosons * lzMax) == totalLz) || (lzMax == 0) || (totalLz == 0))
    return 1;
  int TmpDim = 0;
  while ((totalLz >= 0) && (nbrBosons > 0))
    {
      TmpDim += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax - 1, totalLz);
      --nbrBosons;
      totalLz -= lzMax;
    }
  return TmpDim;
}
