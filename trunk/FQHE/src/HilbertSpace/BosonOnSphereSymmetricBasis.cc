////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of bosons on sphere using                     //
//                            the Lz <-> -Lz symmetry                         //
//                                                                            //
//                        last modification : 07/10/2007                      //
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
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"

#include <math.h>


using std::cout;
using std::endl;


// default constructor
//

BosonOnSphereSymmetricBasis::BosonOnSphereSymmetricBasis ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// lzMax = maximum Lz value reached by a boson

BosonOnSphereSymmetricBasis::BosonOnSphereSymmetricBasis (int nbrBosons, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
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

  int TmpHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    if (this->IsCanonicalState(this->StateDescription[i], this->StateLzMax[i]))
      ++TmpHilbertSpaceDimension;
    else
      this->StateLzMax[i] = -1;
  int** TmpStateDescription = new int* [TmpHilbertSpaceDimension];
  int* TmpStateLzMax = new int [TmpHilbertSpaceDimension];
  TmpHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    if (this->StateLzMax[i] >= 0)
      {
	TmpStateDescription[TmpHilbertSpaceDimension] = this->StateDescription[i];
	TmpStateLzMax[TmpHilbertSpaceDimension] = this->StateLzMax[i];
	++TmpHilbertSpaceDimension;
      }
    else
      delete[] this->StateDescription[i];
  delete[] this->StateDescription;
  delete[] this->StateLzMax;
  this->StateDescription = TmpStateDescription;
  this->StateLzMax = TmpStateLzMax;
  this->HilbertSpaceDimension = TmpHilbertSpaceDimension;  
  this->KeyMultiplicationTable = new int [this->LzMax + 1];
  this->GenerateLookUpTable(0);
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    this->GetStateSymmetry(this->StateDescription[i], this->StateLzMax[i]);
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

BosonOnSphereSymmetricBasis::BosonOnSphereSymmetricBasis(const BosonOnSphereSymmetricBasis& bosons)
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

BosonOnSphereSymmetricBasis::~BosonOnSphereSymmetricBasis ()
{  
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereSymmetricBasis& BosonOnSphereSymmetricBasis::operator = (const BosonOnSphereSymmetricBasis& bosons)
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

AbstractHilbertSpace* BosonOnSphereSymmetricBasis::Clone()
{
  return new BosonOnSphereSymmetricBasis(*this);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereSymmetricBasis::ExtractSubspace (AbstractQuantumNumber& q, 
								    SubspaceSpaceConverter& converter)
{
  return 0;
}

// convert a given state from Lz-symmetric basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector BosonOnSphereSymmetricBasis::ConvertToNbodyBasis(RealVector& state, BosonOnSphere& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      NewLzMax = nbodyBasis.StateLzMax[i];
      int* TmpState = nbodyBasis.StateDescription[i];
      for (int j = 0; j <= NewLzMax; ++j)
	this->TemporaryState[j] = TmpState[j];
      this->GetSignedCanonicalState(this->TemporaryState, NewLzMax);
      if ((NewLzMax & BOSON_SPHERE_SYMMETRIC_BIT) != 0x0ul)	
	TmpVector[i] = state[this->FindStateIndex(this->TemporaryState, NewLzMax)] * M_SQRT1_2;
      else
	TmpVector[i] = state[this->FindStateIndex(this->TemporaryState, NewLzMax)];
    }
  return TmpVector;  
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

int BosonOnSphereSymmetricBasis::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int CurrentLzMax = this->StateLzMax[index];
  int TmpSignature = CurrentLzMax & BOSON_SPHERE_SYMMETRIC_BIT;
  CurrentLzMax &= BOSON_SPHERE_SYMMETRIC_MASK;
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
  this->GetSignedCanonicalState(this->TemporaryState, NewLzMax);
  if ((NewLzMax & BOSON_SPHERE_SYMMETRIC_BIT) != TmpSignature)
    {
      if (TmpSignature != 0)
	coefficient *= M_SQRT2;
      else
	coefficient *= M_SQRT1_2;
    }
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereSymmetricBasis::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  int CurrentLzMax = this->StateLzMax[index];
  int TmpSignature = CurrentLzMax & BOSON_SPHERE_SYMMETRIC_BIT;
  CurrentLzMax &= BOSON_SPHERE_SYMMETRIC_MASK;
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
  this->GetSignedCanonicalState(this->TemporaryState, NewLzMax);
  if ((NewLzMax & BOSON_SPHERE_SYMMETRIC_BIT) != TmpSignature)
    {
      if (TmpSignature != 0)
	coefficient *= M_SQRT2;
      else
	coefficient *= M_SQRT1_2;
    }
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnSphereSymmetricBasis::AA (int index, int n1, int n2)
{
  int CurrentLzMax = this->StateLzMax[index];
  this->ProdASignature = CurrentLzMax & BOSON_SPHERE_SYMMETRIC_BIT;
  CurrentLzMax &= BOSON_SPHERE_SYMMETRIC_MASK;
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

double BosonOnSphereSymmetricBasis::ProdA (int index, int* n, int nbrIndices)
{
  int CurrentLzMax = this->StateLzMax[index];
  this->ProdASignature = CurrentLzMax & BOSON_SPHERE_SYMMETRIC_BIT;
  CurrentLzMax &= BOSON_SPHERE_SYMMETRIC_MASK;
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

int BosonOnSphereSymmetricBasis::AdAd (int m1, int m2, double& coefficient)
{
  int i = 0;
  for (; i < this->NbrLzValue; ++i)
    {
      this->TemporaryState[i] = this->ProdATemporaryState[i];
    }
  ++this->TemporaryState[m2];
  coefficient = this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  this->GetSignedCanonicalState(this->TemporaryState, NewLzMax);
  if ((NewLzMax & BOSON_SPHERE_SYMMETRIC_BIT) != this->ProdASignature)
    {
      if (this->ProdASignature != 0)
	coefficient *= M_SQRT2;
      else
	coefficient *= M_SQRT1_2;
    }
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereSymmetricBasis::ProdAd (int* m, int nbrIndices, double& coefficient)
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
  this->GetSignedCanonicalState(this->TemporaryState, NewLzMax);
  if ((NewLzMax & BOSON_SPHERE_SYMMETRIC_BIT) != this->ProdASignature)
    {
      if (this->ProdASignature != 0)
	coefficient *= M_SQRT2;
      else
	coefficient *= M_SQRT1_2;
    }
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereSymmetricBasis::AdA (int index, int m)
{
  if ((this->StateLzMax[index] & BOSON_SPHERE_SYMMETRIC_MASK) < m)  
    return 0.0;
  return (double) (this->StateDescription[index][m]);  
}

// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnSphereSymmetricBasis::FindStateIndex(int* stateDescription, int lzmax)
{
  lzmax &= BOSON_SPHERE_SYMMETRIC_MASK;
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

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereSymmetricBasis::PrintState (ostream& Str, int state)
{
  int* TmpState = this->StateDescription[state];
  int Max = (this->StateLzMax[state] & BOSON_SPHERE_SYMMETRIC_MASK);
  int i = 0;
  for (; i <= Max; ++i)
    Str << TmpState[i] << " ";
  for (; i <= this->LzMax; ++i)
    Str << "0 ";
  Str << " lzmax  = " << this->StateLzMax[state];
//  Str << " key = " << this->Keys[state] << " lzmax  = " << (this->StateLzMax[state] & BOSON_SPHERE_SYMMETRIC_MASK) << " position = " << FindStateIndex(TmpState, Max);
  return Str;
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex BosonOnSphereSymmetricBasis::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
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

Complex BosonOnSphereSymmetricBasis::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
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

void BosonOnSphereSymmetricBasis::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
  if ((timeCoherence == true) && (this->Minors == 0))
    {
      this->Minors = new Complex* [this->HilbertSpaceDimension];
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	this->Minors[i] = 0;
    }
}


