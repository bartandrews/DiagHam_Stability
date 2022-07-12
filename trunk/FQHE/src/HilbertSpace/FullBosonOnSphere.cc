////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere without total Lz restriction         //
//                                                                            //
//                        last modification : 16/03/2006                      //
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
#include "HilbertSpace/FullBosonOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"

#include <math.h>


using std::cout;
using std::endl;


// basic constructor
// 
// nbrBosons = number of bosons
// lzMax = maximum Lz value reached by a boson

FullBosonOnSphere::FullBosonOnSphere (int nbrBosons, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->LzMax = lzMax;
  this->MaxTotalLz = this->NbrBosons * this->LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax);
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->StateTotalLz = new int [this->HilbertSpaceDimension];
  this->TotalLzPosition = new int [this->MaxTotalLz + 2];
  this->FullGenerateStates(this->NbrBosons, this->LzMax);
  this->KeyMultiplicationTable = new int [this->LzMax + 1];
  this->GenerateLookUpTable(0);
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

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

FullBosonOnSphere::FullBosonOnSphere(const FullBosonOnSphere& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue  = bosons.NbrLzValue;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->StateTotalLz = bosons.StateTotalLz;
  this->TotalLzPosition = bosons.TotalLzPosition;
  this->MaxTotalLz = bosons.MaxTotalLz;
  this->Flag = bosons.Flag;
  this->Keys = bosons.Keys;
  this->KeyMultiplicationTable = bosons.KeyMultiplicationTable;
  this->FullLzMaxPosition = bosons.FullLzMaxPosition;
  this->FullKeyInvertSectorSize = bosons.FullKeyInvertSectorSize;
  this->FullKeyInvertTable = bosons.FullKeyInvertTable;
  this->FullKeyInvertTableNbrIndices = bosons.FullKeyInvertTableNbrIndices;
  this->FullKeyInvertIndices = bosons.FullKeyInvertIndices;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FullBosonOnSphere::~FullBosonOnSphere ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Keys;
      delete[] this->KeyMultiplicationTable;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->StateTotalLz;
      delete[] this->TotalLzPosition;
      for (int k = 0; k <= this->MaxTotalLz; ++k)
	{
	  int TmpLzMax = this->LzMax;
	  if (k < TmpLzMax)
	    {
	      TmpLzMax = k;	  
	    }
	  int Size = (TmpLzMax + 2) * this->IncNbrBosons;
	  for (int i = 0; i < Size; ++i)
	    {
	      if (this->FullKeyInvertSectorSize[k][i] > 0)
		{
		  for (int j= 0; j < this->FullKeyInvertSectorSize[k][i]; ++j)
		    delete[] this->FullKeyInvertIndices[k][i][j];
		  delete[] this->FullKeyInvertTable[k][i];
		  delete[] this->FullKeyInvertTableNbrIndices[k][i];
		  delete[] this->FullKeyInvertIndices[k][i];
		}
	    }
	  delete[] this->FullKeyInvertSectorSize[k];
	  delete[] this->FullKeyInvertTable[k];
	  delete[] this->FullKeyInvertTableNbrIndices[k];
	  delete[] this->FullKeyInvertIndices[k];
	  delete[] this->FullLzMaxPosition[k];
	}      
      delete[] this->FullKeyInvertSectorSize;
      delete[] this->FullKeyInvertTable;
      delete[] this->FullKeyInvertTableNbrIndices;
      delete[] this->FullKeyInvertIndices;
      delete[] this->FullLzMaxPosition;

      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
  this->HilbertSpaceDimension = 0;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FullBosonOnSphere& FullBosonOnSphere::operator = (const FullBosonOnSphere& bosons)
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
      delete[] this->StateTotalLz;
      delete[] this->TotalLzPosition;
      for (int k = 0; k <= this->MaxTotalLz; ++k)
	{
	  int TmpLzMax = this->LzMax;
	  if (k < TmpLzMax)
	    {
	      TmpLzMax = k;	  
	    }
	  int Size = (TmpLzMax + 2) * this->IncNbrBosons;
	  for (int i = 0; i < Size; ++i)
	    {
	      if (this->FullKeyInvertSectorSize[k][i] > 0)
		{
		  for (int j= 0; j < this->FullKeyInvertSectorSize[k][i]; ++j)
		    delete[] this->FullKeyInvertIndices[k][i][j];
		  delete[] this->FullKeyInvertTable[k][i];
		  delete[] this->FullKeyInvertTableNbrIndices[k][i];
		  delete[] this->FullKeyInvertIndices[k][i];
		}
	    }
	  delete[] this->FullKeyInvertSectorSize[k];
	  delete[] this->FullKeyInvertTable[k];
	  delete[] this->FullKeyInvertTableNbrIndices[k];
	  delete[] this->FullKeyInvertIndices[k];
	  delete[] this->FullLzMaxPosition[k];
	}      
      delete[] this->FullKeyInvertSectorSize;
      delete[] this->FullKeyInvertTable;
      delete[] this->FullKeyInvertTableNbrIndices;
      delete[] this->FullKeyInvertIndices;
      delete[] this->FullLzMaxPosition;

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
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->StateTotalLz = bosons.StateTotalLz;
  this->TotalLzPosition = bosons.TotalLzPosition;
  this->MaxTotalLz = bosons.MaxTotalLz;
  this->Flag = bosons.Flag;
  this->Keys = bosons.Keys;
  this->KeyMultiplicationTable = bosons.KeyMultiplicationTable;
  this->FullLzMaxPosition = bosons.FullLzMaxPosition;
  this->FullKeyInvertSectorSize = bosons.FullKeyInvertSectorSize;
  this->FullKeyInvertTable = bosons.FullKeyInvertTable;
  this->FullKeyInvertTableNbrIndices = bosons.FullKeyInvertTableNbrIndices;
  this->FullKeyInvertIndices = bosons.FullKeyInvertIndices;
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

AbstractHilbertSpace* FullBosonOnSphere::Clone()
{
  return new FullBosonOnSphere(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FullBosonOnSphere::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  for (int i = -this->MaxTotalLz; i <= this->MaxTotalLz; i += 2)
    TmpList += new SzQuantumNumber (i);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FullBosonOnSphere::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->StateTotalLz[index]);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FullBosonOnSphere::ExtractSubspace (AbstractQuantumNumber& q, 
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

int FullBosonOnSphere::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
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
  int DestIndex = this->FindStateIndex(this->TemporaryState, NewLzMax, this->StateTotalLz[index] + m1 + m2 - n1 - n2);
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

int FullBosonOnSphere::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  --nbrIndices;
  int CurrentTotalLz = this->StateTotalLz[index];
  for (int i = 0; i <= nbrIndices; ++i)
    {
      if ((n[i] > CurrentLzMax) || (State[n[i]] == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      CurrentTotalLz += m[i];
      CurrentTotalLz -= n[i];
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
  int DestIndex = this->FindStateIndex(this->TemporaryState, NewLzMax, CurrentTotalLz);
  return DestIndex;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FullBosonOnSphere::AA (int index, int n1, int n2)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  this->ProdATemporaryStateTotalLz = this->StateTotalLz[index] - n1 - n2;
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
  return Coefficient;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FullBosonOnSphere::ProdA (int index, int* n, int nbrIndices)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  int TmpCoefficient = 1;
  this->ProdATemporaryStateTotalLz = this->StateTotalLz[index];
  for (i = nbrIndices - 1; i >= 0; --i)
    {
      if (n[i] > CurrentLzMax)
        {
          return 0.0;
        }
      this->ProdATemporaryStateTotalLz -= n[i];
      if (this->ProdATemporaryStateTotalLz < 0)
	return 0.0;
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

int FullBosonOnSphere::AdAd (int m1, int m2, double& coefficient)
{
  int i = 0;  
  for (; i < this->NbrLzValue; ++i)
    {
      this->TemporaryState[i] = this->ProdATemporaryState[i];
    }
  ++this->TemporaryState[m2];
  coefficient *= this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax, this->ProdATemporaryStateTotalLz + m1 + m2);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FullBosonOnSphere::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  int i = 0;
  int TmpTotalLz =  this->ProdATemporaryStateTotalLz;
  for (; i < this->NbrLzValue; ++i)
    {
      this->TemporaryState[i] = this->ProdATemporaryState[i];
    }
  int TmpCoefficient = 1;
  for (i = 0; i < nbrIndices; ++i)
    {
      TmpCoefficient *= ++this->TemporaryState[m[i]];
      TmpTotalLz += m[i];
    }
  coefficient = sqrt((double) TmpCoefficient);
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax, TmpTotalLz);
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FullBosonOnSphere::AdA (int index, int m, int n, double& coefficient)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  if ((n > CurrentLzMax) || (State[n] == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = CurrentLzMax;
  if (NewLzMax < m)
    NewLzMax = m;
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->TemporaryState[i] = State[i];
  for (; i <= NewLzMax; ++i)
    this->TemporaryState[i] = 0;
  coefficient = this->TemporaryState[n];
  --this->TemporaryState[n];
  ++this->TemporaryState[m];
  coefficient *= this->TemporaryState[m];
  coefficient = sqrt(coefficient);
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax, this->StateTotalLz[index] + m - n);
}

// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// totalLz = state total Lz value
// return value = corresponding index

int FullBosonOnSphere::FindStateIndex(int* stateDescription, int lzmax, int totalLz)
{
  int TmpKey = this->GenerateKey(stateDescription, lzmax);
  int Sector = lzmax * this->IncNbrBosons + stateDescription[lzmax];
  int TmpPos = 0;
  int TmpPos2 = this->FullKeyInvertSectorSize[totalLz][Sector] - 1;
  int TmpPos3;
  int* TmpKeyInvertTable = this->FullKeyInvertTable[totalLz][Sector];
  while ((TmpPos2 - TmpPos) > 1)
    {
      TmpPos3 = (TmpPos2 + TmpPos) >> 1;
      if (TmpKey < TmpKeyInvertTable[TmpPos3])
	TmpPos2 = TmpPos3;
      else
	TmpPos = TmpPos3;
    }
  if (TmpKey == TmpKeyInvertTable[TmpPos2])
    TmpPos = TmpPos2;
  int i;
  int* TmpStateDescription;
  int Start;
  int* TmpKeyInvertIndices = this->FullKeyInvertIndices[totalLz][Sector][TmpPos];

  TmpPos2 = this->FullKeyInvertTableNbrIndices[totalLz][Sector][TmpPos] - 1;
  TmpPos = 0;
  while ((TmpPos2 - TmpPos) > 1)
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
	TmpPos = TmpPos3;
      else
	TmpPos2 = TmpPos3;
    }
  Start = TmpKeyInvertIndices[TmpPos];
  TmpStateDescription = this->StateDescription[Start];
  i = lzmax;
  while ((i >= 0) && (stateDescription[i] ==  TmpStateDescription[i]))
    --i;
  if (i == -1)
    {
      return TmpKeyInvertIndices[TmpPos];
    }
  return TmpKeyInvertIndices[TmpPos2];
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FullBosonOnSphere::PrintState (ostream& Str, int state)
{
  int* TmpState = this->StateDescription[state];
  int Max = this->StateLzMax[state];
  int i = 0;
  for (; i <= Max; ++i)
    Str << TmpState[i] << " ";
  for (; i <= this->LzMax; ++i)
    Str << "0 ";
  Str << " key = " << this->Keys[state] << " lzmax  = " << this->StateLzMax[state] << " totallz  = " << this->StateTotalLz[state] << " position = " << FindStateIndex(TmpState, Max, this->StateTotalLz[state]);
//    Str << " key = " << this->Keys[state] << " lzmax position = " << this->LzMaxPosition[Max * (this->NbrBosons + 1) + TmpState[Max]]
//        << " position = " << FindStateIndex(TmpState, Max);
  return Str;
}

// generate all states
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson in the state
// return value = Hilbert space size

int FullBosonOnSphere::FullGenerateStates(int nbrBosons, int lzMax)
{
  int Pos = 0;
  for (int i = 0; i <= this->MaxTotalLz; ++i)
    {
      int TmpLzMax = lzMax;
      if (i < TmpLzMax)
	{
	  TmpLzMax = i;	  
	}
      int TmpPos = this->GenerateStates(nbrBosons, TmpLzMax, TmpLzMax, i, Pos);
      for (int j = Pos; j < TmpPos; ++j)
	{
	  this->StateTotalLz[j] = i;
	}
      this->TotalLzPosition[i] = Pos;
      Pos = TmpPos;
    }
  this->TotalLzPosition[this->MaxTotalLz + 1] = Pos;
  return Pos;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FullBosonOnSphere::GenerateLookUpTable(int memory)
{
  this->Keys = new int [this->HilbertSpaceDimension];
  for (int i = 0; i <= this->LzMax; ++i)
    this->KeyMultiplicationTable[i] = i * i * i * i;//this->IncNbrBosons;
  this->FullLzMaxPosition = new int* [this->MaxTotalLz + 1];
  this->FullKeyInvertSectorSize = new int* [this->MaxTotalLz + 1];
  this->FullKeyInvertTable = new int** [this->MaxTotalLz + 1];
  this->FullKeyInvertTableNbrIndices = new int** [this->MaxTotalLz + 1];
  this->FullKeyInvertIndices = new int*** [this->MaxTotalLz + 1];

  for (int i = 0; i <= this->MaxTotalLz; ++i)
    {
      int TmpLzMax = this->LzMax;
      if (i < TmpLzMax)
	{
	  TmpLzMax = i;	  
	}
      int Size = (TmpLzMax + 2) * this->IncNbrBosons;
      this->FullLzMaxPosition[i] = new int [Size];
      this->FullKeyInvertSectorSize[i] = new int [Size];
      this->FullKeyInvertTable[i] = new int* [Size];
      this->FullKeyInvertTableNbrIndices[i] = new int* [Size];
      this->FullKeyInvertIndices[i] = new int** [Size];
      this->CoreGenerateLookUpTable(this->TotalLzPosition[i + 1] - this->TotalLzPosition[i], TmpLzMax,
				    this->StateDescription + this->TotalLzPosition[i], this->StateLzMax + this->TotalLzPosition[i], this->Keys + this->TotalLzPosition[i], 
				    this->FullLzMaxPosition[i], this->FullKeyInvertSectorSize[i], this->FullKeyInvertTable[i], 
				    this->FullKeyInvertTableNbrIndices[i], this->FullKeyInvertIndices[i], this->TotalLzPosition[i]);
    }
}


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

int FullBosonOnSphere::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax)
{
  BinomialCoefficients TmpCoefficients (lzMax + nbrBosons);
  return TmpCoefficients(lzMax + nbrBosons, nbrBosons);
}



