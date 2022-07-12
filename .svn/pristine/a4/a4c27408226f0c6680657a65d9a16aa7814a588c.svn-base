////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of bosons on sphere wih SU(2) spin                 //
//                                                                            //
//                        last modification : 10/10/2008                      //
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
#include "HilbertSpace/BosonOnSphereWithSpin.h"
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

BosonOnSphereWithSpin::BosonOnSphereWithSpin ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson
// totalSpin = twce the total spin value

BosonOnSphereWithSpin::BosonOnSphereWithSpin (int nbrBosons, int totalLz, int lzMax, int totalSpin)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->TotalSpin = totalSpin;
  this->NbrBosonsUp = (this->NbrBosons + this->TotalSpin)/2;
  this->NbrBosonsDown = (this->NbrBosons - this->TotalSpin)/2;
  this->HilbertSpaceDimension = (int) this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, this->TotalLz, this->TotalSpin);
  cout << "dim = " << this->HilbertSpaceDimension << endl;
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
  int TmpDim = this->GenerateStates(this->NbrBosonsUp, this->NbrBosonsDown, TmpLzMax, TmpLzMax, this->ShiftedTotalLz, 0);
  cout << "dim 2 = " << TmpDim << endl;
  for (int i = 0; i < TmpDim; ++i)
    this->PrintState(cout, i) << endl;
  this->KeyMultiplicationTable = new int [this->LzMax + 1];
  this->GenerateLookUpTable(0);
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;

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

BosonOnSphereWithSpin::BosonOnSphereWithSpin(const BosonOnSphereWithSpin& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->TotalSpin = bosons.TotalSpin;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
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
}

// destructor
//

BosonOnSphereWithSpin::~BosonOnSphereWithSpin ()
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

BosonOnSphereWithSpin& BosonOnSphereWithSpin::operator = (const BosonOnSphereWithSpin& bosons)
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
  this->TotalSpin = bosons.TotalSpin;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
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

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereWithSpin::Clone()
{
  return new BosonOnSphereWithSpin(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereWithSpin::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereWithSpin::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereWithSpin::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m1_d a^+_m2_d a_n1_d a_n2_d operator to a given state (with m1+m2=n1+n2, only spin down)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  cout << "WARNING :  AddAduAdAu method is nopt implemented" << endl;
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_u a^+_m2_u a_n1_u a_n2_u operator to a given state (with m1+m2=n1+n2, only spin up)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  cout << "WARNING :  AddAduAdAu method is nopt implemented" << endl;
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// apply a^+_d_m1 a^+_u_m2 a_d_n1 a_u_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  cout << "WARNING :  AddAduAdAu method is nopt implemented" << endl;
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}


// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 

double BosonOnSphereWithSpin::AuAu (int index, int n1, int n2)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || ((State[n1] >> 16) == 0) || ((State[n2] >> 16) == 0) || ((n1 == n2) && ((State[n1] >> 16) == 1)))
    {
      return 0.0;
    }
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  int& TmpDown = this->ProdATemporaryState[n2];
  int TmpUp = TmpDown >> 16;
  TmpDown &= 0xffff;
  double Coefficient = TmpUp;
  --TmpUp;
  TmpDown |= TmpUp;
  int& TmpDown2 = this->ProdATemporaryState[n1];
  TmpUp = TmpDown2 >> 16;
  TmpDown2 &= 0xffff;
  Coefficient *= TmpUp;
  --TmpUp;
  TmpDown2 |= TmpUp;
  for (i = CurrentLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt(Coefficient);
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double BosonOnSphereWithSpin::AdAd (int index, int n1, int n2)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || ((State[n1] & 0xffff) == 0) || ((State[n2] & 0xffff) == 0) || ((n1 == n2) && ((State[n1] & 0xffff) == 1)))
    {
      return 0.0;
    }
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  int& TmpDown = this->ProdATemporaryState[n2];
  int TmpUp = TmpDown >> 16;
  TmpDown &= 0xffff;
  double Coefficient = TmpDown;
  --TmpDown;
  TmpDown |= TmpUp;
  int& TmpDown2 = this->ProdATemporaryState[n1];
  TmpUp = TmpDown2 >> 16;
  TmpDown2 &= 0xffff;
  Coefficient *= TmpDown2;
  --TmpDown2;
  TmpDown2 |= TmpUp;
  for (i = CurrentLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt(Coefficient);
}

// apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double BosonOnSphereWithSpin::AuAd (int index, int n1, int n2)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || ((State[n1] >> 16) == 0) || ((State[n2] & 0xffff) == 0))
    {
      return 0.0;
    }
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  int& TmpDown = this->ProdATemporaryState[n2];
  int TmpUp = TmpDown >> 16;
  TmpDown &= 0xffff;
  double Coefficient = TmpUp;
  --TmpUp;
  TmpDown |= TmpUp;
  int& TmpDown2 = this->ProdATemporaryState[n1];
  TmpUp = TmpDown2 >> 16;
  TmpDown2 &= 0xffff;
  Coefficient *= TmpDown2;
  --TmpDown2;
  TmpDown2 |= TmpUp;
  for (i = CurrentLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt(Coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereWithSpin::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
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
      int& TmpDown = this->ProdATemporaryState[n[i]];
      int TmpUp = TmpDown >> 16;
      TmpDown &= 0xffff;
      if (spinIndices[i] == 0)
	{
	  if (TmpDown == 0)
	    return 0.0;
	  TmpCoefficient *= TmpDown;
	  --TmpDown;
	}
      else 
	{
	  if (TmpUp == 0)
	    return 0.0;
	  TmpCoefficient *= TmpUp;
	  --TmpUp;
	}
      TmpDown |= TmpUp << 16;
    }
  for (i = CurrentLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt((double) TmpCoefficient);
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AduAdu (int m1, int m2, double& coefficient)
{
  int i = 0;
  for (; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  int& TmpDown = this->TemporaryState[m2];
  int TmpUp = TmpDown >> 16;
  TmpDown &= 0xffff;
  ++TmpUp;
  coefficient = TmpUp;
  TmpDown |= TmpUp;
  int& TmpDown2 = this->TemporaryState[m1];
  TmpUp = TmpDown2 >> 16;
  TmpDown2 &= 0xffff;
  ++TmpUp;
  coefficient *= TmpUp;
  TmpDown2 |= TmpUp;
  coefficient = sqrt(coefficient);
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AddAdd (int m1, int m2, double& coefficient)
{
  int i = 0;
  for (; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  int& TmpDown = this->TemporaryState[m2];
  int TmpUp = TmpDown >> 16;
  TmpDown &= 0xffff;
  ++TmpDown;
  coefficient = TmpDown;
  TmpDown |= TmpUp;
  int& TmpDown2 = this->TemporaryState[m1];
  TmpUp = TmpDown2 >> 16;
  TmpDown2 &= 0xffff;
  ++TmpDown2;
  coefficient *= TmpDown2;
  TmpDown2 |= TmpUp;
  coefficient = sqrt(coefficient);
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// apply a^+_m1_d a^+_m2_u operator to the state produced using AuAd method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AduAdd (int m1, int m2, double& coefficient)
{
  int i = 0;
  for (; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  int& TmpDown = this->TemporaryState[m2];
  int TmpUp = TmpDown >> 16;
  TmpDown &= 0xffff;
  ++TmpUp;
  coefficient = TmpUp;
  TmpDown |= TmpUp;
  int& TmpDown2 = this->TemporaryState[m1];
  TmpUp = TmpDown2 >> 16;
  TmpDown2 &= 0xffff;
  ++TmpDown2;
  coefficient *= TmpDown2;
  TmpDown2 |= TmpUp;
  coefficient = sqrt(coefficient);
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  int i = 0;
  for (; i < this->NbrLzValue; ++i)
    {
      this->TemporaryState[i] = this->ProdATemporaryState[i];
    }
  int TmpCoefficient = 1;
  for (i = 0; i < nbrIndices; ++i)
    {
      int& TmpDown = this->TemporaryState[m[i]];
      int TmpUp = TmpDown >> 16;
      TmpDown &= 0xffff;
      if (spinIndices[i] == 0)
	{
	  TmpCoefficient *= TmpDown;
	  ++TmpDown;
	}
      else 
	{
	  TmpCoefficient *= TmpUp;
	  ++TmpUp;
	}
      TmpDown |= TmpUp << 16;
    }
  coefficient = sqrt((double) TmpCoefficient);
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereWithSpin::AduAu (int index, int m)
{
  if (this->StateLzMax[index] < m)  
    return 0.0;
  return (double) ((this->StateDescription[index][m]) >> 16);  
}

// apply a^+_m_d a_m_ operator to a given state  (only spin down)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereWithSpin::AddAd (int index, int m)
{
  if (this->StateLzMax[index] < m)  
    return 0.0;
  return (double) ((this->StateDescription[index][m]) & 0xffff);  
}

// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnSphereWithSpin::FindStateIndex(int* stateDescription, int lzmax)
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

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSpin::PrintState (ostream& Str, int state)
{
  int* TmpState = this->StateDescription[state];
  int Max = this->StateLzMax[state];
  int i = 0;
  Str << "|";
  for (; i <= Max; ++i)
    {
      if (TmpState[i] == 0)
	Str << "  0 ";
      else
	{
	  if ((TmpState[i] & 0xffff0000) != 0) 
	    Str << (TmpState[i] >> 16) << "u";
	  else 
	    Str << "  ";
	  if ((TmpState[i] & 0xffff) != 0) 
	    Str << (TmpState[i] & 0xffff) << "d";
	  else
	    Str << "  ";
	}
      Str << "|";
    }
  for (; i <= this->LzMax; ++i)
    Str << "  0 |";
  //  Str << " key = " << this->Keys[state] << " lzmax  = " << this->StateLzMax[state]<< " position = " << FindStateIndex(TmpState, Max);
//    Str << " key = " << this->Keys[state] << " lzmax position = " << this->LzMaxPosition[Max * (this->NbrBosons + 1) + TmpState[Max]]
//        << " position = " << FindStateIndex(TmpState, Max);
  Str << " *** " << this->StateLzMax[state] << endl;
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosonsUp = number of bosons with spin up
// nbrBosonsDown = number of bosons with spin down
// lzMax = momentum maximum value for a boson in the state
// currentLzMax = momentum maximum value for bosons that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int BosonOnSphereWithSpin::GenerateStates(int nbrBosonsUp, int nbrBosonsDown, int lzMax, int currentLzMax, int totalLz, int pos)
{
  if (((nbrBosonsUp + nbrBosonsDown) == 0) || (totalLz < 0) || (((nbrBosonsUp + nbrBosonsDown) * currentLzMax) < totalLz))
    {
      return pos;
    }
  if (((nbrBosonsUp + nbrBosonsDown) * currentLzMax) == totalLz)
    {
      this->StateDescription[pos] = new int [lzMax + 1];
      int* TmpState = this->StateDescription[pos];
      for (int i = 0; i <= lzMax; ++i)
	TmpState[i] = 0;
      TmpState[currentLzMax] = (nbrBosonsUp << 16) | nbrBosonsDown;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }
  if ((currentLzMax == 0) || (totalLz == 0))
    {
      this->StateDescription[pos] = new int [lzMax + 1];
      int* TmpState = this->StateDescription[pos];
      for (int i = 1; i <= lzMax; ++i)
	TmpState[i] = 0;
      TmpState[0] = (nbrBosonsUp << 16) | nbrBosonsDown;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }
  if ((nbrBosonsDown + nbrBosonsUp) == 1)
    if (nbrBosonsUp == 1)
      {
	this->StateDescription[pos] = new int [lzMax + 1];
	int* TmpState = this->StateDescription[pos];
	for (int i = 0; i <= lzMax; ++i)
	  TmpState[i] = 0;
	TmpState[totalLz] = 1 << 16;
	this->StateLzMax[pos] = lzMax;
	return pos + 1;
      }
    else
      {
	this->StateDescription[pos] = new int [lzMax + 1];
	int* TmpState = this->StateDescription[pos];
	for (int i = 0; i <= lzMax; ++i)
	  TmpState[i] = 0;
	TmpState[totalLz] = 1;
	this->StateLzMax[pos] = lzMax;
	return pos + 1;
      }
  int TotalNbrBosons = nbrBosonsDown + nbrBosonsUp;
  int TmpTotalLz = totalLz - (TotalNbrBosons * currentLzMax);
  TmpTotalLz += currentLzMax;
  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = pos;
  for (int TmpTotalNbrBosons = 1; TmpTotalNbrBosons < TotalNbrBosons; ++TmpTotalNbrBosons)
    {
      int TmpNbrBosonsUp = nbrBosonsUp;
      if (TmpNbrBosonsUp > TmpTotalNbrBosons)
	TmpNbrBosonsUp = TmpTotalNbrBosons;
      for (int TmpNbrBosonsDown = 0; (TmpNbrBosonsDown <= nbrBosonsDown) && (TmpNbrBosonsDown <= TmpTotalNbrBosons) && (TmpNbrBosonsUp >= 0); ++TmpNbrBosonsDown)
	{
	  TmpPos = this->GenerateStates(nbrBosonsUp - TmpNbrBosonsUp, TmpNbrBosonsDown, lzMax, ReducedCurrentLzMax, TmpTotalLz, pos);
	  for (int i = pos; i < TmpPos; i++)
	    this->StateDescription[i][currentLzMax] = (TmpNbrBosonsUp << 16) | (nbrBosonsDown - TmpNbrBosonsDown);
	  pos = TmpPos;	  
	  --TmpNbrBosonsUp;
	}
      TmpTotalLz += currentLzMax;
    }


  if (lzMax == currentLzMax)
    return this->GenerateStates(nbrBosonsUp, nbrBosonsDown, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, pos);
  else
    return this->GenerateStates(nbrBosonsUp, nbrBosonsDown, lzMax, ReducedCurrentLzMax, totalLz, pos);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereWithSpin::GenerateLookUpTable(int memory)
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

void BosonOnSphereWithSpin::CoreGenerateLookUpTable(int dimension, int lzMax, int** stateDescription, int* stateLzMax, int* keys, int* lzMaxPosition, int* keyInvertSectorSize, 
					    int** keyInvertTable, int** keyInvertTableNbrIndices, int*** keyInvertIndices, int indexShift)
{
  int Size = (lzMax + 2) * this->IncNbrBosons;
  for (int i = 0; i < Size; ++i)
    keyInvertSectorSize[i] =0;
  int CurrentLzMax = stateLzMax[0];
  int CurrentNbrLzMax = (stateDescription[0][CurrentLzMax] & 0xffff) + (stateDescription[0][CurrentLzMax] >> 16);
  lzMaxPosition[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 0; 
  keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 0;   
  for (int i = 0; i < dimension; ++i)
    {
      keys[i] = this->GenerateKey(stateDescription[i], stateLzMax[i]);
      if (CurrentLzMax != stateLzMax[i])
	{
	  CurrentLzMax = stateLzMax[i];
	  CurrentNbrLzMax = (stateDescription[i][CurrentLzMax] >> 16) + (stateDescription[i][CurrentLzMax] & 0xffff);
	  lzMaxPosition[CurrentLzMax * this->IncNbrBosons + CurrentNbrLzMax] = indexShift + i; 
	  keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	}
      else
	if (stateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    CurrentNbrLzMax = (stateDescription[i][CurrentLzMax] >> 16) + (stateDescription[i][CurrentLzMax] & 0xffff);
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
  CurrentNbrLzMax = (stateDescription[0][CurrentLzMax] >> 16) + (stateDescription[0][CurrentLzMax] & 0xffff);
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
	  CurrentNbrLzMax = (stateDescription[i][CurrentLzMax] >> 16) + (stateDescription[i][CurrentLzMax] & 0xffff);
	  CurrentKeyInvertSectorSize = keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  keyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	  TmpKeyInvertTable = keyInvertTable[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  TmpKeyInvertTable[0] = keys[i];
	  TmpKeyInvertTableNbrIndices[0] = 1;
	}
      else
	{
	  int TmpNbrLzMax = (stateDescription[i][CurrentLzMax] >> 16) + (stateDescription[i][CurrentLzMax] & 0xffff);
	  if (TmpNbrLzMax != CurrentNbrLzMax)
	    {
	      CurrentNbrLzMax = TmpNbrLzMax;
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
	  CurrentNbrLzMax = (stateDescription[i][CurrentLzMax] >> 16) + (stateDescription[i][CurrentLzMax] & 0xffff);
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
	    CurrentNbrLzMax = (stateDescription[i][CurrentLzMax] >> 16) + (stateDescription[i][CurrentLzMax] & 0xffff);
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

int BosonOnSphereWithSpin::GenerateKey(int* stateDescription, int lzmax)
{
  int Key = 0;
  for (int i = 0; i <= lzmax; ++i)
    {
      Key += this->KeyMultiplicationTable[i] * ((stateDescription[i] & 0xffff) + (stateDescription[i] >> 16));
    }
  return Key;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// totalSpin = twice the total spin value
// return value = Hilbert space dimension

long BosonOnSphereWithSpin::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalSpin)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, 
						    (totalLz + lzMax * nbrBosons) >> 1, (totalSpin + nbrBosons) >> 1);
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// return value = Hilbert space dimension      

long BosonOnSphereWithSpin::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalSpin)
{
  if ((nbrBosons < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalSpin > nbrBosons) || ((lzMax * nbrBosons) < totalLz))
    return 0l;
    
  if (nbrBosons == 1) 
    if (lzMax >= totalLz)
      return 1l;
    else
      return 0l;
  if (totalLz == 0)
    return 1l;

  unsigned long Tmp = 0l;  
  for (int i = totalSpin; i >= 0; --i)
    for (int j = (nbrBosons - totalSpin); j >= 0; --j)
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - (i + j), lzMax - 1, totalLz - (lzMax * (i + j)), totalSpin - i);
  return Tmp;
  //  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - 1, lzMax, totalLz - lzMax, totalSpin - 1);
  //  return  (Tmp + this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax - 1, totalLz, totalSpin));
}

