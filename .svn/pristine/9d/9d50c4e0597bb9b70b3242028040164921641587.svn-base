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
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/StringTools.h"

#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//

BosonOnSphere::BosonOnSphere ()
{
}

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
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  int TmpLzMax = this->LzMax;
  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  this->GenerateStates(this->NbrBosons, TmpLzMax, TmpLzMax, this->ShiftedTotalLz, 0);
  this->KeyMultiplicationTable = new int [this->LzMax + 1];
  this->GenerateLookUpTable(0);
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

//   for (int i=0; i<HilbertSpaceDimension; ++i)
//     {
//       PrintState(cout,i)<<endl;      
//     }

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
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->Flag = bosons.Flag;
  this->Keys = bosons.Keys;
  this->KeyMultiplicationTable = bosons.KeyMultiplicationTable;
  this->LzMaxPosition = bosons.LzMaxPosition;
  this->KeyInvertSectorSize = bosons.KeyInvertSectorSize;
  this->KeyInvertTable = bosons.KeyInvertTable;
  this->KeyInvertTableNbrIndices = bosons.KeyInvertTableNbrIndices;
  this->KeyInvertIndices = bosons.KeyInvertIndices;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

BosonOnSphere::~BosonOnSphere ()
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Keys;
      delete[] this->KeyMultiplicationTable;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      int Size = (this->LzMax + 2) * this->IncNbrBosons;
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
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
      delete[] this->LzMaxPosition;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphere& BosonOnSphere::operator = (const BosonOnSphere& bosons)
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Keys;
      delete[] this->KeyMultiplicationTable;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      int Size = (this->LzMax + 2) * this->IncNbrBosons;
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
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
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
  this->Keys = bosons.Keys;
  this->KeyMultiplicationTable = bosons.KeyMultiplicationTable;
  this->LzMaxPosition = bosons.LzMaxPosition;
  this->KeyInvertSectorSize = bosons.KeyInvertSectorSize;
  this->KeyInvertTable = bosons.KeyInvertTable;
  this->KeyInvertTableNbrIndices = bosons.KeyInvertTableNbrIndices;
  this->KeyInvertIndices = bosons.KeyInvertIndices;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

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
  for (; i <= CurrentLzMax; ++i)
    this->TemporaryState[i] = State[i];
  for (; i <= NewLzMax; ++i)
    this->TemporaryState[i] = 0;
  coefficient = this->TemporaryState[n2];
  --this->TemporaryState[n2];
  coefficient *= this->TemporaryState[n1];
  --this->TemporaryState[n1];
  ++this->TemporaryState[m2];
  coefficient *= this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  int DestIndex = this->FindStateIndex(this->TemporaryState, NewLzMax);
  return DestIndex;
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphere::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  --nbrIndices;
  for (int i = 0; i <= nbrIndices; ++i)
    {
      if ((n[i] > CurrentLzMax) || (State[n[i]] == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
 
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->TemporaryState[i] = State[i];
  for (; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = 0;
  coefficient = 1.0;
  for (i = nbrIndices; i >= 0; --i)
    {
      if (this->TemporaryState[n[i]] == 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension; 	    
	}
      coefficient *= (double) this->TemporaryState[n[i]];
      --this->TemporaryState[n[i]];
    }
  for (i = nbrIndices; i >= 0; --i)
    {
      ++this->TemporaryState[m[i]];
      coefficient *= (double) this->TemporaryState[m[i]];
    }
  coefficient = sqrt(coefficient);
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  int DestIndex = this->FindStateIndex(this->TemporaryState, NewLzMax);
  return DestIndex;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnSphere::AA (int index, int n1, int n2)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || (State[n1] == 0) || (State[n2] == 0) || ((n1 == n2) && (State[n1] == 1)))
    {
      return 0.0;
    }
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  double Coefficient = this->ProdATemporaryState[n2];
  --this->ProdATemporaryState[n2];
  Coefficient *= this->ProdATemporaryState[n1];
  --this->ProdATemporaryState[n1];
  for (i = CurrentLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt(Coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphere::ProdA (int index, int* n, int nbrIndices)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  int TmpCoefficient = 1;
  for (i = nbrIndices - 1; i >= 0; --i)
    {
      if (n[i] > CurrentLzMax)
        {
          return 0.0;
        }
      int& Tmp = this->ProdATemporaryState[n[i]];
      if (Tmp == 0)
        {
          return 0.0;
        }
      TmpCoefficient *= Tmp;
      --Tmp;
    }
  for (i = CurrentLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt((double) TmpCoefficient);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphere::AdAd (int m1, int m2, double& coefficient)
{
  int i = 0;
  for (; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m2];
  coefficient = this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphere::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  int i = 0;
  for (; i < this->NbrLzValue; ++i)
    {
      this->TemporaryState[i] = this->ProdATemporaryState[i];
    }
  int TmpCoefficient = 1;
  for (i = 0; i < nbrIndices; ++i)
    TmpCoefficient *= ++this->TemporaryState[m[i]];
  coefficient = sqrt((double) TmpCoefficient);
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
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

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphere::AdA (long index, int m)
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
  int* TmpKeyInvertIndices = this->KeyInvertIndices[Sector][TmpPos];

  TmpPos2 =this->KeyInvertTableNbrIndices[Sector][TmpPos] - 1;
  TmpPos = 0;
  while (TmpPos2 != TmpPos)
    {
      TmpPos3 = (TmpPos2 + TmpPos) >> 1;
      Start = TmpKeyInvertIndices[TmpPos3];
      TmpStateDescription = this->StateDescription[Start];
      i = lzmax;
      while ((i >= 0) && (stateDescription[i] ==  TmpStateDescription[i]))
        --i;
      if (i == -1)
        {
          return Start;
        }
      if (stateDescription[i] < TmpStateDescription[i])
        {
          TmpPos = TmpPos3 + 1;
        }
      else
        {
           TmpPos2 = TmpPos3 - 1;
        }
    }
  return TmpKeyInvertIndices[TmpPos];
}

// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int BosonOnSphere::FindStateIndex(char* stateDescription)
{
  char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != (this->LzMax + 1))
    return -1;
  int TmpNbrParticles = 0;
  int TmpTotalLz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      int Tmp = atoi(TmpDescription[i]);
      this->TemporaryState[i] = Tmp;
      TmpTotalLz += (i * Tmp);
      TmpNbrParticles += Tmp;
      delete[] TmpDescription[i];
    }
  delete[] TmpDescription;
  if ((TmpNbrParticles != this->NbrBosons) || (TmpTotalLz != ((this->TotalLz + this->NbrBosons * this->LzMax) >> 1)))
    return -1;
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the Lz component
int BosonOnSphere::GetLzValue(int j)
{
  return this->TotalLz;
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
  Str << " key = " << this->Keys[state] << " lzmax  = " << this->StateLzMax[state]<< " position = " << FindStateIndex(TmpState, Max);
//   Str << " key = " << this->Keys[state] << " lzmax position = " << this->LzMaxPosition[Max * (this->NbrBosons + 1) + TmpState[Max]]
//       << " position = " << FindStateIndex(TmpState, Max);
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
  this->CoreGenerateLookUpTable(this->HilbertSpaceDimension, this->LzMax, this->StateDescription, this->StateLzMax, this->Keys, this->LzMaxPosition, this->KeyInvertSectorSize, 
				this->KeyInvertTable, this->KeyInvertTableNbrIndices, this->KeyInvertIndices);
}

// generate look-up table associated to current Hilbert space (core part of the look-up table generation)
// 
// dimension = Hilbert space dimension
// lzMax = maximum Lz value that can be reached by a particle
// stateDescription = array that contains state description
// stateLzMax = array giving maximum Lz value reached for a boson in a given state
// keys = keys associated to each state
// lzMaxPosition = indicate position of the first state with a given number of boson having a given maximum Lz value
// keyInvertSectorSize = array that indicates how many different states are store for each sector
// keyInvertTable = array that contains sorted possible key for each sector
// keyInvertTableNbrIndices = array that contains number of indices that have the same key per sector 
// keyInvertIndices = array that contains state index per sector and per key
// indexShift = optional shift to apply before storing any index

void BosonOnSphere::CoreGenerateLookUpTable(int dimension, int lzMax, int** stateDescription, int* stateLzMax, int* keys, int* lzMaxPosition, int* keyInvertSectorSize, 
					    int** keyInvertTable, int** keyInvertTableNbrIndices, int*** keyInvertIndices, int indexShift)
{
  int Size = (lzMax + 2) * this->IncNbrBosons;
  for (int i = 0; i < Size; ++i)
    keyInvertSectorSize[i] =0;
  int CurrentLzMax = stateLzMax[0];
  int CurrentNbrLzMax = stateDescription[0][CurrentLzMax];
  lzMaxPosition[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 0; 
  keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 0;   
  for (int i = 0; i < dimension; ++i)
    {
      keys[i] = this->GenerateKey(stateDescription[i], stateLzMax[i]);
      if (CurrentLzMax != stateLzMax[i])
	{
	  CurrentLzMax = stateLzMax[i];
	  CurrentNbrLzMax = stateDescription[i][CurrentLzMax];
	  lzMaxPosition[CurrentLzMax * this->IncNbrBosons + CurrentNbrLzMax] = indexShift + i; 
	  keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	}
      else
	if (stateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    CurrentNbrLzMax = stateDescription[i][CurrentLzMax];
	    lzMaxPosition[CurrentLzMax * this->IncNbrBosons + CurrentNbrLzMax] = indexShift + i;
	    keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	  }
	else
	  {
	    ++keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax]; 
	  }
    }

  for (int i = 0; i < Size; ++i)
    {
      if (keyInvertSectorSize[i] > 0)
	{
	  keyInvertTable[i] = new int [keyInvertSectorSize[i]];
	  keyInvertTableNbrIndices[i] = new int [keyInvertSectorSize[i]];
	}
      else
	{
	  keyInvertTable[i] = 0; 
	  keyInvertTableNbrIndices[i] = 0;
	}
    }

  CurrentLzMax = stateLzMax[0];
  CurrentNbrLzMax = stateDescription[0][CurrentLzMax];
  int CurrentKeyInvertSectorSize = keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
  int* TmpKeyInvertTable = keyInvertTable[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  int* TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  TmpKeyInvertTable[0] = keys[0];
  TmpKeyInvertTableNbrIndices[0] = 1;
  for (int i = 1; i < dimension; ++i)
    {
      if (CurrentLzMax != stateLzMax[i])
	{
	  CurrentLzMax = stateLzMax[i];
	  CurrentNbrLzMax = stateDescription[i][CurrentLzMax];
	  CurrentKeyInvertSectorSize = keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	  TmpKeyInvertTable = keyInvertTable[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  TmpKeyInvertTable[0] = keys[i];
	  TmpKeyInvertTableNbrIndices[0] = 1;
	}
      else
	if (stateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    CurrentNbrLzMax = stateDescription[i][CurrentLzMax];
	    CurrentKeyInvertSectorSize = keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	    keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	    TmpKeyInvertTable = keyInvertTable[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	    TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	    TmpKeyInvertTable[0] = keys[i];
	    TmpKeyInvertTableNbrIndices[0] = 1;
	  }
	else
	  {
	    int Lim = keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	    int j = 0;
	    bool Flag = false;
	    int TmpKey = keys[i];
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
		TmpKeyInvertTable[keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax]] = TmpKey;
		TmpKeyInvertTableNbrIndices[keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax]] = 1;
		++keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	      }
	  }
    }
 
  // sort key invert table and resize arrays
  int TmpPos;
  int TmpSize;
  for (int i = 0; i < Size; ++i)
    {
      if (keyInvertSectorSize[i] > 0)
	{
	  int Lim = keyInvertSectorSize[i];
	  TmpKeyInvertTable = keyInvertTable[i];
	  TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[i];
	  int* TmpKeyInvertTable2 = new int [Lim];
	  int* TmpKeyInvertTableNbrIndices2 = new int [Lim];
	  int Tmp;
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
	  keyInvertTable[i] = TmpKeyInvertTable2;
	  keyInvertTableNbrIndices[i] = TmpKeyInvertTableNbrIndices2;
	  keyInvertIndices[i] = new int* [Lim];
	  for (int j = 0; j < Lim; ++j)
	    {
	      keyInvertIndices[i][j] = new int [keyInvertTableNbrIndices[i][j]];
	      keyInvertTableNbrIndices[i][j] = 0;
	    }
	}
    }

  // find all hilbert space indices that have the same key in each sector
  CurrentLzMax = lzMax;
  CurrentNbrLzMax = 1;
  TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  int** TmpKeyInvertIndices = keyInvertIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  int TmpPos2;
  int TmpPos3;
  int TmpPos4 = CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax;
  TmpSize = keyInvertSectorSize[TmpPos4];
  TmpKeyInvertIndices = keyInvertIndices[TmpPos4];
  TmpKeyInvertTable = keyInvertTable[TmpPos4];
  TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[TmpPos4];
  for (int i = 0; i < dimension; ++i)
    {
      int TmpKey = keys[i];
      if (CurrentLzMax != stateLzMax[i])
	{
	  CurrentLzMax = stateLzMax[i];
	  CurrentNbrLzMax = stateDescription[i][CurrentLzMax];
	  TmpPos4 = CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax;
	  TmpKeyInvertIndices = keyInvertIndices[TmpPos4];
	  TmpKeyInvertTable = keyInvertTable[TmpPos4];
	  TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[TmpPos4];
	  TmpSize = keyInvertSectorSize[TmpPos4];
	  TmpPos = 0;
	  TmpPos2 = TmpSize - 1;
	  while ((TmpPos2 - TmpPos) > 1)
	    {
	      TmpPos3 = (TmpPos2 + TmpPos) >> 1;
	      if (TmpKey < TmpKeyInvertTable[TmpPos3])
		TmpPos2 = TmpPos3;
	      else
		TmpPos = TmpPos3;
	    }
	  if (TmpKey == TmpKeyInvertTable[TmpPos])
	    {
	      TmpKeyInvertIndices[TmpPos][0] = indexShift + i;
	      TmpKeyInvertTableNbrIndices[TmpPos] = 1;
	    }
	  else
	    {
	      TmpKeyInvertIndices[TmpPos2][0] = indexShift + i;
	      TmpKeyInvertTableNbrIndices[TmpPos2] = 1;
	    }
	}
      else
	if (stateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    CurrentNbrLzMax = stateDescription[i][CurrentLzMax];
	    TmpPos4 = CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax;
	    TmpKeyInvertIndices = keyInvertIndices[TmpPos4];
	    TmpKeyInvertTable = keyInvertTable[TmpPos4];
	    TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[TmpPos4];
	    TmpSize = keyInvertSectorSize[TmpPos4];
	    TmpPos = 0;
	    TmpPos2 = TmpSize - 1;
	    while ((TmpPos2 - TmpPos) > 1)
	      {
		TmpPos3 = (TmpPos2 + TmpPos) >> 1;
		if (TmpKey < TmpKeyInvertTable[TmpPos3])
		  TmpPos2 = TmpPos3;
		else
		  TmpPos = TmpPos3;
	      }
	    if (TmpKey == TmpKeyInvertTable[TmpPos])
	      {
		TmpKeyInvertIndices[TmpPos][0] = indexShift + i;
		TmpKeyInvertTableNbrIndices[TmpPos] = 1;
	      }
	    else
	      {
		TmpKeyInvertIndices[TmpPos2][0] = indexShift + i;
		TmpKeyInvertTableNbrIndices[TmpPos2] = 1;
	      }
	  }
	else
	  {
	    TmpPos = 0;
	    TmpPos2 = TmpSize - 1;
	    while ((TmpPos2 - TmpPos) > 1)
	      {
		TmpPos3 = (TmpPos2 + TmpPos) >> 1;
		if (TmpKey < TmpKeyInvertTable[TmpPos3])
		  {
		    TmpPos2 = TmpPos3;
		  }
		else
		  TmpPos = TmpPos3;
	      }
	    if (TmpKey == TmpKeyInvertTable[TmpPos])
	      {
		TmpKeyInvertIndices[TmpPos][TmpKeyInvertTableNbrIndices[TmpPos]] = indexShift + i;
		++TmpKeyInvertTableNbrIndices[TmpPos];
	      }
	    else
	      {
		TmpKeyInvertIndices[TmpPos2][TmpKeyInvertTableNbrIndices[TmpPos2]] = indexShift + i;
		++TmpKeyInvertTableNbrIndices[TmpPos2];
	      }
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


// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex BosonOnSphere::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
{
  return this->EvaluateWaveFunction(state, position, basis, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// return value = wave function evaluated at the given location

Complex BosonOnSphere::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
							      int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex BosonOnSphere::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
					     int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
  ComplexMatrix Perm(this->NbrBosons, this->NbrBosons);
  ComplexMatrix Functions(this->LzMax + 1, this->NbrBosons);
  RealVector TmpCoordinates(2);
  int* Indices = new int [this->NbrBosons];
  int Pos;
  int Lz;
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      TmpCoordinates[0] = position[j << 1];
      TmpCoordinates[1] = position[1 + (j << 1)];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	  Functions[j].Re(i) = Tmp.Re;
	  Functions[j].Im(i) = Tmp.Im;
	}
    }
  double* Factors = new double [this->NbrBosons + 1];
  Factors[0] = 1.0;
  Factors[1] = 1.0;
  for (int i = 2; i <= this->NbrBosons; ++i)
    Factors[i] = Factors[i - 1] / sqrt((double) i);
  double TmpFactor;
  int* ChangeBitSign;
  int* ChangeBit;
  int TmpStateDescription;
  Perm.EvaluateFastPermanentPrecalculationArray(ChangeBit, ChangeBitSign);
  int LastComponent = firstComponent + nbrComponent;
  for (int k = firstComponent; k < LastComponent; ++k)
    {
      Pos = 0;
      Lz = 0;
      TmpFactor = state[k] * Factors[this->NbrBosons];
      while (Pos < this->NbrBosons)
	{
	  TmpStateDescription = this->StateDescription[k][Lz];
	  if (TmpStateDescription != 0)
	    {
	      TmpFactor *= Factors[TmpStateDescription];
	      for (int j = 0; j < TmpStateDescription; ++j)
		{
		  Indices[Pos] = Lz;
		  ++Pos;
		}
	    }
	  ++Lz;
	}
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < this->NbrBosons; ++j)
	    {
	      Perm[i].Re(j) = TmpColum2.Re(Indices[j]);
	      Perm[i].Im(j) = TmpColum2.Im(Indices[j]);
	    }
	}
      Value += Perm.FastPermanent(ChangeBit, ChangeBitSign) * TmpFactor;
    }
  delete[] ChangeBitSign;
  delete[] ChangeBit;
  delete[] Factors;
  delete[] Indices;
  return Value;
}

// evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex BosonOnSphere::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
							      int nextCoordinates, int firstComponent, int nbrComponent)
{
  double* Factors = new double [this->NbrBosons + 1];
  Factors[0] = 1.0;
  Factors[1] = 1.0;
  for (int i = 2; i <= this->NbrBosons; ++i)
    Factors[i] = Factors[i - 1] / sqrt((double) i);
  double TmpFactor;
  Complex Value;
  Complex Tmp;
  Complex TmpPerm;
  int* Indices = new int [this->NbrBosons];
  int Pos;
  int Lz;
  int TmpStateDescription;
  int LastComponent = firstComponent + nbrComponent;
  if ((*(this->KeptCoordinates)) == -1)
    {
      ComplexMatrix Perm(this->NbrBosons, this->NbrBosons);
      ComplexMatrix Functions(this->LzMax + 1, this->NbrBosons);
      RealVector TmpCoordinates(2);
      for (int j = 0; j < this->NbrBosons; ++j)
	{
	  TmpCoordinates[0] = position[j << 1];
	  TmpCoordinates[1] = position[1 + (j << 1)];
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	      Functions[j].Re(i) = Tmp.Re;
	      Functions[j].Im(i) = Tmp.Im;
	    }
	}
      int* ChangeBitSign;
      int* ChangeBit;
      int TmpStateDescription;
      Perm.EvaluateFastPermanentPrecalculationArray(ChangeBit, ChangeBitSign, true);
      int LastComponent = firstComponent + nbrComponent;
      for (int k = firstComponent; k < LastComponent; ++k)
	{
	  Pos = 0;
	  Lz = 0;
	  TmpFactor = state[k] * Factors[this->NbrBosons];
	  while (Pos < this->NbrBosons)
	    {
	      TmpStateDescription = this->StateDescription[k][Lz];
	      if (TmpStateDescription != 0)
		{
		  TmpFactor *= Factors[TmpStateDescription];
		  for (int j = 0; j < TmpStateDescription; ++j)
		    {
		      Indices[Pos] = Lz;
		      ++Pos;
		    }
		}
	      ++Lz;
	    }
	  for (int i = 0; i < this->NbrBosons; ++i)
	    {
	      ComplexVector& TmpColum2 = Functions[i];	  
	      for (int j = 0; j < this->NbrBosons; ++j)
		{
		  Perm[i].Re(j) = TmpColum2.Re(Indices[j]);
		  Perm[i].Im(j) = TmpColum2.Im(Indices[j]);
		}
	    }
	  if (this->Minors[k] == 0)
	    {
	      this->Minors[k] = new Complex [this->NbrBosons];
	    }
	  Perm.FastPermanentMinorDevelopment(ChangeBit, ChangeBitSign, nextCoordinates, this->Minors[k]);
	  TmpPerm = 0.0;
	  for (int i = 0; i < this->NbrBosons; ++i)
	    TmpPerm += this->Minors[k][i] * Complex (Perm[nextCoordinates].Re(i), 
						     Perm[nextCoordinates].Im(i));
	  Value += TmpPerm * TmpFactor;
	}
      delete[] ChangeBitSign;
      delete[] ChangeBit;
      (*(this->KeptCoordinates)) = nextCoordinates;
    }
  else
    {
      Complex* Functions = new Complex[this->LzMax + 1];
      RealVector TmpCoordinates(2);
      TmpCoordinates[0] = position[(*(this->KeptCoordinates)) << 1];
      TmpCoordinates[1] = position[1 + ((*(this->KeptCoordinates)) << 1)];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  basis.GetFunctionValue(TmpCoordinates, Functions[i], i);
	}
      for (int k = firstComponent; k < LastComponent; ++k)
	{
	  Pos = 0;
	  Lz = 0;
	  TmpFactor = Factors[this->NbrBosons] * state[k];
	  while (Pos < this->NbrBosons)
	    {
	      TmpStateDescription = this->StateDescription[k][Lz];
	      if (TmpStateDescription != 0)
		{
		  TmpFactor *= Factors[TmpStateDescription];
		  for (int j = 0; j < TmpStateDescription; ++j)
		    {
		      Indices[Pos] = Lz;
		      ++Pos;
		    }
		}
	      ++Lz;
	    }
	  Complex* TmpMinors = this->Minors[k];
	  TmpPerm = 0.0;
	  for (int i = 0; i < this->NbrBosons; ++i)
	    TmpPerm += TmpMinors[i] * Functions[Indices[i]];
	  Value += TmpPerm * TmpFactor;
	}
      delete[] Functions;
      (*(this->KeptCoordinates)) = -1;
    }
  delete[] Factors;
  delete[] Indices;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void BosonOnSphere::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
  if ((timeCoherence == true) && (this->Minors == 0))
    {
      this->Minors = new Complex* [this->HilbertSpaceDimension];
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	this->Minors[i] = 0;
    }
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  BosonOnSphere::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState)
{  
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrBosonSector == 0))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrBosonSector == this->NbrBosons))
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * groundState[j]);
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0;
	  for (int MinIndex = 0; MinIndex < this->HilbertSpaceDimension; ++MinIndex)
	    if (this->StateDescription[MinIndex][0] == nbrBosonSector)
	       TmpValue += groundState[MinIndex] * groundState[MinIndex];
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}      
    }
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0;
	  int MinIndex = 0;
	  while ((MinIndex < this->HilbertSpaceDimension) && (this->StateLzMax[MinIndex] >= subsytemSize))
	    {
	      int TmpPos = 0;
	      int* TmpState = this->StateDescription[MinIndex];
	      while ((TmpPos < subsytemSize) && (TmpState[TmpPos] == 0))
		++TmpPos;
	      if (TmpPos == subsytemSize)
		{
		  TmpValue += groundState[MinIndex] * groundState[MinIndex];
		}
	      ++MinIndex;
	    }
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  if (ShiftedLzComplementarySector < 0)
    {
      RealSymmetricMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  int MinIndex = 0;
  int MaxIndex = this->HilbertSpaceDimension - 1;
  if (nbrBosonSector == 1)
    {
      double TmpValue = 0.0;
      int TmpLzMax = this->StateLzMax[MinIndex];
      while ((MinIndex <= MaxIndex) && (subsytemSize <= TmpLzMax))
	{
	  int* TmpStateDescription = this->StateDescription[MinIndex];
	  if (TmpStateDescription[ShiftedLzSector] == 1)
	    {	      
	      int TmpPos = 0;
	      int TmpNbrBosons = 0;
	      while (TmpPos < subsytemSize)
		TmpNbrBosons += TmpStateDescription[TmpPos++];
	      if (TmpNbrBosons == 1)
		TmpValue += groundState[MinIndex] * groundState[MinIndex];	    
	    }
	  ++MinIndex;
	  if (MinIndex <= MaxIndex)
	    TmpLzMax = this->StateLzMax[MinIndex];
	}
      RealSymmetricMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  if (NbrBosonsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
      BosonOnSphere TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      double TmpValue;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpValue = groundState[MinIndex + i];
	  for (int j = i; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    TmpDensityMatrix.SetMatrixElement(i, j, TmpValue * groundState[MinIndex + j]);
	}
      return TmpDensityMatrix;
    }


  int TmpNbrBosons;
  int TmpTotalLz;
  int TmpIndex;
  BosonOnSphere TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  long TmpNbrNonZeroElements = 0;
  int TmpComplementarySubsystemLzMax = this->StateLzMax[MinIndex];
  while ((MinIndex <= MaxIndex) && (TmpComplementarySubsystemLzMax >= subsytemSize))
    {
      int* TmpComplementarySubsystem = this->StateDescription[MinIndex];
      TmpIndex = MinIndex + 1;
      int TmpPos = TmpComplementarySubsystemLzMax;	  
      int TmpLzMax = TmpPos;
      while ((TmpIndex <= MaxIndex) && (TmpPos == TmpLzMax))
	{
	  if (TmpLzMax == this->StateLzMax[TmpIndex])
	    {
	      TmpPos = subsytemSize;
	      while ((TmpPos <= TmpLzMax) && (this->StateDescription[TmpIndex][TmpPos] == TmpComplementarySubsystem[TmpPos]))
		++TmpPos;
	      if (TmpPos > TmpLzMax)
		{
		  ++TmpIndex;
		  --TmpPos;
		}
	      else
		{
		  TmpPos = -1;
		}
	    }
	  else
	    {
	      TmpPos = -1;
	    }
	}
      TmpNbrBosons = 0;
      TmpTotalLz = 0;
      TmpPos = subsytemSize;	  
      while (TmpPos <= TmpComplementarySubsystemLzMax)
	{
	  TmpNbrBosons += TmpComplementarySubsystem[TmpPos];
	  TmpTotalLz += TmpComplementarySubsystem[TmpPos] * TmpPos;
	  ++TmpPos;
	}
      if ((TmpNbrBosons == NbrBosonsComplementarySector) && (ShiftedLzComplementarySector == TmpTotalLz))
	{
	  int Pos = 0;
	  for (int i = MinIndex; i < TmpIndex; ++i)
	    {
	      int* TmpState = this->StateDescription[i];
	      int TmpLzMax = subsytemSize - 1;
	      while (TmpState[TmpLzMax] == 0) 
		--TmpLzMax;
	      TmpStatePosition[Pos] = TmpDestinationHilbertSpace.FindStateIndex(TmpState, TmpLzMax);
	      ++Pos;
	    }
	  int Pos2;
	  Pos = 0;
	  int Pos3;
	  double TmpValue;
	  for (int i = MinIndex; i < TmpIndex; ++i)
	    {
	      Pos2 = 0;
	      Pos3 = TmpStatePosition[Pos];
	      TmpValue = groundState[i];
	      for (int j = MinIndex; j < TmpIndex; ++j)
		{
		  if (Pos3 <=  TmpStatePosition[Pos2])
		    {
		      TmpDensityMatrix.AddToMatrixElement(Pos3, TmpStatePosition[Pos2], TmpValue * groundState[j]);
		      ++TmpNbrNonZeroElements;
		    }
		  ++Pos2;
		}
	      ++Pos;
	    }
	}
      MinIndex = TmpIndex;
      if (MinIndex <= MaxIndex)
	TmpComplementarySubsystemLzMax = StateLzMax[MinIndex];
    }
  delete[] TmpStatePosition;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}
