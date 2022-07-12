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
  if ((nbrBosons^totalSpin)&1)
    {
      cout << "NbrBosons and TotalSpin need to have the same parity"<<endl;
      exit(1);
    }
  this->NbrBosons = nbrBosons;  
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->TotalSpin = totalSpin;
  this->NbrBosonsUp = (this->NbrBosons + this->TotalSpin) >> 1;
  this->NbrBosonsDown = (this->NbrBosons - this->TotalSpin) >> 1;
  if (NbrBosonsUp>NbrBosonsDown)
    this->IncMaxNbrBosons = this->NbrBosonsUp + 1;
  else
    this->IncMaxNbrBosons = this->NbrBosonsDown + 1;
  this->HilbertSpaceDimension = (int) this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, this->TotalLz, this->TotalSpin);
  cout << "dim = " << this->HilbertSpaceDimension << endl;
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateLzSzMax = new int [this->HilbertSpaceDimension];
  int TmpLzMax = this->LzMax;
  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  this->GenerateStates(this->NbrBosonsUp, this->NbrBosonsDown, TmpLzMax, TmpLzMax, this->ShiftedTotalLz, 0);  
  this->KeyMultiplicationTable = new int [2*(this->LzMax + 1)];
  this->GenerateLookUpTable(0);
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->MinorsUp = 0;
  this->MinorsDown = 0;

  //  for (int i=0; i<HilbertSpaceDimension; ++i)
  //    PrintState(cout,i)<<endl;

#ifdef __DEBUG__
  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += (this->StateLzSzMax[i] + 1) * sizeof(int) + sizeof(int*);
  UsedMemory += (this->TotalLz + 1) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += ((this->TotalLz + 1) * this->IncMaxNbrBosons) * sizeof(int);
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
  this->IncMaxNbrBosons = bosons.IncMaxNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->TotalSpin = bosons.TotalSpin;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzSzMax = bosons.StateLzSzMax;
  this->Flag = bosons.Flag;
  this->Keys = bosons.Keys;
  this->KeyMultiplicationTable = bosons.KeyMultiplicationTable;
  this->LzSzMaxPosition = bosons.LzSzMaxPosition;
  this->KeyInvertSectorSize = bosons.KeyInvertSectorSize;
  this->KeyInvertTable = bosons.KeyInvertTable;
  this->KeyInvertTableNbrIndices = bosons.KeyInvertTableNbrIndices;
  this->KeyInvertIndices = bosons.KeyInvertIndices;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->MinorsUp = bosons.MinorsUp;
  this->MinorsDown = bosons.MinorsDown;
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
      delete[] this->StateLzSzMax;
      int Size = (2*(this->LzMax+1) + 1) * this->IncMaxNbrBosons;
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

      if (this->MinorsUp != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      if (this->MinorsUp[i] != 0)
		delete[] this->MinorsUp[i];	      
	      delete[] this->MinorsUp;
	    }
	}
      if (this->MinorsDown != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      if (this->MinorsDown[i] != 0)
		delete[] this->MinorsDown[i];	      
	      delete[] this->MinorsDown;
	    }
	}
      delete this->KeptCoordinates;
      delete[] this->LzSzMaxPosition;
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
      delete[] this->StateLzSzMax;
      int Size = (this->LzMax + 2) * this->IncMaxNbrBosons;
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
      if (this->MinorsUp != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      if (this->MinorsUp[i] != 0)
		delete[] this->MinorsUp[i];	      
	      delete[] this->MinorsUp;
	    }
	}
      if (this->MinorsDown != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      if (this->MinorsDown[i] != 0)
		delete[] this->MinorsDown[i];	      
	      delete[] this->MinorsDown;
	    }
	}
      delete this->KeptCoordinates;
    }
  this->NbrBosons = bosons.NbrBosons;
  this->IncMaxNbrBosons = bosons.IncMaxNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzSzMax = bosons.StateLzSzMax;
  this->Flag = bosons.Flag;
  this->Keys = bosons.Keys;
  this->KeyMultiplicationTable = bosons.KeyMultiplicationTable;
  this->LzSzMaxPosition = bosons.LzSzMaxPosition;
  this->KeyInvertSectorSize = bosons.KeyInvertSectorSize;
  this->KeyInvertTable = bosons.KeyInvertTable;
  this->KeyInvertTableNbrIndices = bosons.KeyInvertTableNbrIndices;
  this->KeyInvertIndices = bosons.KeyInvertIndices;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->MinorsUp = bosons.MinorsUp;
  this->MinorsDown = bosons.MinorsDown;
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
  return this->HilbertSpaceDimension;
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
  this->TemporaryState[m2] += 0x10000;
  coefficient = (this->TemporaryState[m2] >> 16);
  this->TemporaryState[m1] += 0x10000;
  coefficient *= (this->TemporaryState[m1] >> 16);
  coefficient = sqrt(coefficient);
  int NewLzSzMax = this->LzMax;
  while (this->TemporaryState[NewLzSzMax] == 0)
    --NewLzSzMax;
  if (this->TemporaryState[NewLzSzMax]&0xffff0000)
    {
      NewLzSzMax<<=1;
      NewLzSzMax+=1;
    }
  else NewLzSzMax<<=1;
  return this->FindStateIndex(this->TemporaryState, NewLzSzMax);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AdAd method (without destroying it)
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
  ++this->TemporaryState[m2];
  coefficient = (this->TemporaryState[m2] & 0xffff);
  ++this->TemporaryState[m1];
  coefficient *= (this->TemporaryState[m1] & 0xffff);
  coefficient = sqrt(coefficient);
  int NewLzSzMax = this->LzMax;
  while (this->TemporaryState[NewLzSzMax] == 0)
    --NewLzSzMax;  
  if (this->TemporaryState[NewLzSzMax]&0xffff0000)
    {
      NewLzSzMax<<=1;
      NewLzSzMax+=1;
    }
  else NewLzSzMax<<=1;
  return this->FindStateIndex(this->TemporaryState, NewLzSzMax);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAd method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AduAdd (int m1, int m2, double& coefficient)
{
  int i = 0;
  for (; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m2];
  coefficient = (this->TemporaryState[m2] & 0xffff);
  this->TemporaryState[m1] += 0x10000;
  coefficient *= (this->TemporaryState[m1] >> 16);
  coefficient = sqrt(coefficient);
  int NewLzSzMax = this->LzMax;
  while (this->TemporaryState[NewLzSzMax] == 0)
    --NewLzSzMax;  
  if (this->TemporaryState[NewLzSzMax]&0xffff0000)
    {
      NewLzSzMax<<=1;
      NewLzSzMax+=1;
    }
  else NewLzSzMax<<=1;    
  return this->FindStateIndex(this->TemporaryState, NewLzSzMax);
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 

double BosonOnSphereWithSpin::AuAu (int index, int n1, int n2)
{
  int CurrentLzMax = this->StateLzSzMax[index]>>1;
  int* State = this->StateDescription[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || ((State[n1] >> 16) == 0) || ((State[n2] >> 16) == 0) || ((n1 == n2) && ((State[n1] >> 16) == 1)))
    {
      return 0.0;
    }
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  double Coefficient = (this->ProdATemporaryState[n2] >> 16);
  this->ProdATemporaryState[n2] -= 0x10000;
  Coefficient *= (this->ProdATemporaryState[n1] >> 16);
  this->ProdATemporaryState[n1] -= 0x10000;
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
  int CurrentLzMax = this->StateLzSzMax[index]>>1;
  int* State = this->StateDescription[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || ((State[n1] & 0xffff) == 0) || ((State[n2] & 0xffff) == 0) || ((n1 == n2) && ((State[n1] & 0xffff) == 1)))
    {
      return 0.0;
    }
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  double Coefficient = (this->ProdATemporaryState[n2] & 0xffff);
  --this->ProdATemporaryState[n2];
  Coefficient *= (this->ProdATemporaryState[n1] & 0xffff);
  --this->ProdATemporaryState[n1];
  for (i = CurrentLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt(Coefficient);
}

// apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double BosonOnSphereWithSpin::AuAd (int index, int n1, int n2)
{
  int CurrentLzMax = this->StateLzSzMax[index]>>1;
  int* State = this->StateDescription[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || ((State[n1] >> 16) == 0) || ((State[n2] & 0xffff) == 0))
    {
      return 0.0;
    }
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  double Coefficient = (this->ProdATemporaryState[n2] & 0xffff);
  --this->ProdATemporaryState[n2];
  Coefficient *= (this->ProdATemporaryState[n1] >> 16);
  this->ProdATemporaryState[n1] -= 0x10000;
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
  int CurrentLzMax = this->StateLzSzMax[index]>>1;
  int* State = this->StateDescription[index];
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  int TmpCoefficient = 1;
  for (i = nbrIndices - 1; i >= 0; --i)
    {
      if (n[i] > CurrentLzMax)
	return 0.0;
      int& Tmp = this->ProdATemporaryState[n[i]];
      if (spinIndices[i] == 0)
	{
	  if ((Tmp & 0xffff) == 0)
	    return 0.0;
	  TmpCoefficient *= (Tmp & 0xffff);
	  --Tmp;
	}
      else
	{
	  if ((Tmp >> 16) == 0)
	    return 0.0;
	  TmpCoefficient *= (Tmp >> 16);
	  Tmp -= 0x10000;
	}
    }
  for (i = CurrentLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt((double) TmpCoefficient);
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
    if (spinIndices[i] == 0)
      {
	++this->TemporaryState[m[i]];
	TmpCoefficient *= (this->TemporaryState[m[i]] &0xffff);
      }
    else
      {
	this->TemporaryState[m[i]] += 0x10000;
	TmpCoefficient *= (this->TemporaryState[m[i]] >> 16);
      }
  coefficient = sqrt((double) TmpCoefficient);
  int NewLzSzMax = this->LzMax;
  while (this->TemporaryState[NewLzSzMax] == 0)
    --NewLzSzMax;
  if (this->TemporaryState[NewLzSzMax]&0xffff0000)
    {
      NewLzSzMax<<=1;
      NewLzSzMax+=1;
    }
  else NewLzSzMax<<=1;
  return this->FindStateIndex(this->TemporaryState, NewLzSzMax);
}


// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereWithSpin::AduAu (int index, int m)
{
  if (this->StateLzSzMax[index]>>1 < m)  
    return 0.0;
  return (double) ((this->StateDescription[index][m] >> 16));  
}

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereWithSpin::AddAd (int index, int m)
{
  if (this->StateLzSzMax[index]>>1 < m)  
    return 0.0;
  return (double) ((this->StateDescription[index][m] & 0xffff));  
}

// find state index
//
// stateDescription = array describing the state
// lzszmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnSphereWithSpin::FindStateIndex(int* stateDescription, int lzszmax)
{
//   cout << "Searching: ";
//   PrintState(cout, stateDescription)<<" ";
  int TmpKey = this->GenerateKey(stateDescription, lzszmax);
//   cout << "key="<<TmpKey<<endl;
  int Offset;
  if (lzszmax&1)
    Offset=stateDescription[lzszmax>>1]>>16;
  else
    Offset=stateDescription[lzszmax>>1]&0xffff;
  int Sector = lzszmax * this->IncMaxNbrBosons + Offset;
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
      i = lzszmax;
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

int BosonOnSphereWithSpin::FindStateIndex(char* stateDescription)
{
  char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != (this->LzMax + 1))
    return -1;
  int TmpNbrParticles = 0;
  int TmpTotalLz = 0;
  //int TmpTotalSz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      cout << "Need to fix BosonOnSphereWithSpin::FindStateIndex(char* stateDescription)"<<endl;
      int Tmp = atoi(TmpDescription[i]);
      this->TemporaryState[i] = Tmp;
      TmpTotalLz += (i * Tmp);
      TmpNbrParticles += Tmp;
      delete[] TmpDescription[i];
    }
  delete[] TmpDescription;
  if ((TmpNbrParticles != this->NbrBosons) || (TmpTotalLz != ((this->TotalLz + this->NbrBosons * this->LzMax) >> 1)))
    return -1;
  int NewLzSzMax = this->LzMax;
  while (this->TemporaryState[NewLzSzMax] == 0)
    --NewLzSzMax;
  if (this->TemporaryState[NewLzSzMax]&0xffff0000)
    {
      NewLzSzMax<<=1;
      NewLzSzMax+=1;
    }
  else NewLzSzMax<<=1;
  return this->FindStateIndex(this->TemporaryState, NewLzSzMax);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSpin::PrintState (ostream& Str, int state)
{
  int* TmpState = this->StateDescription[state];
  int Max = (this->StateLzSzMax[state]>>1);
  int i = 0;
  Str << state << " = " << "|";
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
//   Str << " lzszmax  = " << this->StateLzSzMax[state];
//   Str << " key = " << GenerateKey(TmpState, this->StateLzSzMax[state]);
//   Str << " position = " << FindStateIndex(TmpState, this->StateLzSzMax[state]);
  return Str;
}

// print a given State
//
// Str = reference on current output stream 
// stateDesc = array describing the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSpin::PrintState (ostream& Str, int* stateDesc)
{
  int Max = this->LzMax;
  while((stateDesc[Max]==0)&&(Max>0)) --Max;
  int i = 0;
  Str << " |";
  for (; i <= Max; ++i)
    {
      if (stateDesc[i] == 0)
	Str << "  0 ";
      else
	{
	  if ((stateDesc[i] & 0xffff0000) != 0) 
	    Str << (stateDesc[i] >> 16) << "u";
	  else 
	    Str << "  ";
	  if ((stateDesc[i] & 0xffff) != 0) 
	    Str << (stateDesc[i] & 0xffff) << "d";
	  else
	    Str << "  ";
	}
      Str << "|";
    }
  for (; i <= this->LzMax; ++i)
    Str << "  0 |";
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
  //cout << "GenerateStates(Up:"<< nbrBosonsUp<<", Down:"<< nbrBosonsDown<<", lzMax:"<<lzMax<<", currLz:"<<currentLzMax<< ", totalLz:" <<totalLz<<", pos:"<< pos<<")"<<endl;
  if ((nbrBosonsUp < 0) || (nbrBosonsDown < 0) || (totalLz < 0) || (currentLzMax < 0) || (((nbrBosonsUp + nbrBosonsDown) * currentLzMax) < totalLz))
    return pos;
  if ((nbrBosonsUp + nbrBosonsDown) == 0)
    {
      if (totalLz == 0)
	{	  
	  this->StateDescription[pos] = new int [lzMax + 1];
	  int* TmpState = this->StateDescription[pos];
	  for (int i = 0; i <= lzMax; ++i)
	    TmpState[i] = 0;
	  this->StateLzSzMax[pos] = 0;
	  return pos + 1;
	}
      else return pos;
    }
  if (((nbrBosonsUp + nbrBosonsDown) * currentLzMax) == totalLz)
    {
      this->StateDescription[pos] = new int [lzMax + 1];
      int* TmpState = this->StateDescription[pos];
      for (int i = 0; i <= lzMax; ++i)
	TmpState[i] = 0;
      TmpState[currentLzMax] = (nbrBosonsUp << 16) | nbrBosonsDown;
      if (nbrBosonsUp>0)
	this->StateLzSzMax[pos] = (currentLzMax<<1)+1;
      else
	this->StateLzSzMax[pos] = currentLzMax<<1;
      return pos + 1;
    }
  if ((nbrBosonsDown + nbrBosonsUp) == 1)
    if (lzMax >= totalLz)
      {
	this->StateDescription[pos] = new int [lzMax + 1];
	int* TmpState = this->StateDescription[pos];
	for (int i = 0; i <= lzMax; ++i)
	  TmpState[i] = 0;
	if (nbrBosonsUp == 1)
	  {
	    TmpState[totalLz] = 1 << 16; // place spin up
	    this->StateLzSzMax[pos] = (totalLz<<1)+1;
	  }
	else
	  {
	    TmpState[totalLz] = 1; // place spin down
	    this->StateLzSzMax[pos] = totalLz<<1;
	  }
	return pos + 1;	
      }
    else return pos;
  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = pos;
  int TmpLzSzMax;
  for (int UpToPlace=nbrBosonsUp; UpToPlace>=0; --UpToPlace)
    for (int DownToPlace=nbrBosonsDown; DownToPlace>=0; --DownToPlace)    
      {
	TmpPos = this->GenerateStates(nbrBosonsUp-UpToPlace, nbrBosonsDown-DownToPlace, lzMax, ReducedCurrentLzMax, totalLz-(UpToPlace+DownToPlace)*currentLzMax, pos);
	TmpLzSzMax = (currentLzMax<<1);
	if (UpToPlace>0) ++TmpLzSzMax;
	for (int i = pos; i < TmpPos; i++)
	  {
	    this->StateDescription[i][currentLzMax] = (UpToPlace << 16) | DownToPlace;
	    if ((UpToPlace>0)||(DownToPlace>0))
	      this->StateLzSzMax[i] = TmpLzSzMax;
	  }	
	pos = TmpPos;
      }
  return pos;
}


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereWithSpin::GenerateLookUpTable(int memory)
{
  this->Keys = new int [this->HilbertSpaceDimension];
  this->KeyMultiplicationTable[0]=this->IncMaxNbrBosons;
  for (int i = 1; i < 2*(this->LzMax+1); ++i)
    this->KeyMultiplicationTable[i] = this->IncMaxNbrBosons*this->KeyMultiplicationTable[i-1];
  int Size = (2*(this->LzMax+1) + 1) * this->IncMaxNbrBosons;
  this->LzSzMaxPosition = new int [Size];
  this->KeyInvertSectorSize = new int [Size];
  this->KeyInvertTable = new int* [Size];
  this->KeyInvertTableNbrIndices = new int* [Size];
  this->KeyInvertIndices = new int** [Size];
  this->CoreGenerateLookUpTable(this->HilbertSpaceDimension, this->LzMax, this->StateDescription, this->StateLzSzMax, this->Keys, this->LzSzMaxPosition, this->KeyInvertSectorSize, 
				this->KeyInvertTable, this->KeyInvertTableNbrIndices, this->KeyInvertIndices);
}

// generate look-up table associated to current Hilbert space (core part of the look-up table generation)
// 
// dimension = Hilbert space dimension
// lzMax = maximum Lz value that can be reached by a particle
// stateDescription = array that contains state description
// stateLzMax = array giving maximum Lz value reached for a boson in a given state
// keys = keys associated to each state
// lzSzMaxPosition = indicate position of the first state with a given number of boson having a given maximum Lz value
// keyInvertSectorSize = array that indicates how many different states are store for each sector
// keyInvertTable = array that contains sorted possible key for each sector
// keyInvertTableNbrIndices = array that contains number of indices that have the same key per sector 
// keyInvertIndices = array that contains state index per sector and per key
// indexShift = optional shift to apply before storing any index

void BosonOnSphereWithSpin::CoreGenerateLookUpTable(int dimension, int lzMax, int** stateDescription, int* stateLzSzMax, int* keys, int* lzSzMaxPosition, int* keyInvertSectorSize, 
					    int** keyInvertTable, int** keyInvertTableNbrIndices, int*** keyInvertIndices, int indexShift)
{
  int Size = (2*(this->LzMax+1) + 1) * this->IncMaxNbrBosons;
  for (int i = 0; i < Size; ++i)
    keyInvertSectorSize[i] = 0;
  int CurrentLzSzMax = stateLzSzMax[0];
  int CurrentNbrLzSzMax, NewNbrLzSzMax;
  if (CurrentLzSzMax&1)
    CurrentNbrLzSzMax = (stateDescription[0][CurrentLzSzMax>>1])>>16;
  else
    CurrentNbrLzSzMax = (stateDescription[0][CurrentLzSzMax>>1])&0xffff;
  lzSzMaxPosition[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax] = 0; 
  keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax] = 0;   
  for (int i = 0; i < dimension; ++i)
    {
      keys[i] = this->GenerateKey(stateDescription[i], stateLzSzMax[i]);
      if (CurrentLzSzMax != stateLzSzMax[i])
	{
	  CurrentLzSzMax = stateLzSzMax[i];
	  if (CurrentLzSzMax&1)
	    CurrentNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])>>16;
	  else
	    CurrentNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])&0xffff;
	  lzSzMaxPosition[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax] = indexShift + i; 
	  keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax] = 1;
	}
      else
	{
	  if (stateLzSzMax[i]&1)
	    NewNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])>>16;
	  else
	    NewNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])&0xffff;	  
	  if (NewNbrLzSzMax != CurrentNbrLzSzMax)
	    {
	      CurrentNbrLzSzMax = NewNbrLzSzMax;
	      lzSzMaxPosition[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax] = indexShift + i;
	      keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax] = 1;
	    }
	  else
	    {
	      ++keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax]; 
	    }
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

  CurrentLzSzMax = stateLzSzMax[0];
  if (CurrentLzSzMax&1)
    CurrentNbrLzSzMax = (stateDescription[0][CurrentLzSzMax>>1])>>16;
  else
    CurrentNbrLzSzMax = (stateDescription[0][CurrentLzSzMax>>1])&0xffff;
  int CurrentKeyInvertSectorSize = keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
  keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax] = 1;
  int* TmpKeyInvertTable = keyInvertTable[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
  int* TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
  //cout << "keys="<<keys<<", CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax="<<CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax<<" TmpKeyInvertTable="<<TmpKeyInvertTable<<endl;
  TmpKeyInvertTable[0] = keys[0];
  TmpKeyInvertTableNbrIndices[0] = 1;
  for (int i = 1; i < dimension; ++i)
    {
      if (CurrentLzSzMax != stateLzSzMax[i])
	{
	  CurrentLzSzMax = stateLzSzMax[i];
	  if (CurrentLzSzMax&1)
	    CurrentNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])>>16;
	  else
	    CurrentNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])&0xffff;
	  CurrentKeyInvertSectorSize = keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
	  keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax] = 1;
	  TmpKeyInvertTable = keyInvertTable[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
	  TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
	  TmpKeyInvertTable[0] = keys[i];
	  TmpKeyInvertTableNbrIndices[0] = 1;
	}
      else
	{
	  if (stateLzSzMax[i]&1)
	    NewNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])>>16;
	  else
	    NewNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])&0xffff;	  
	  if (NewNbrLzSzMax != CurrentNbrLzSzMax)
	    {
	      CurrentNbrLzSzMax = NewNbrLzSzMax;
	      CurrentKeyInvertSectorSize = keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
	      keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax] = 1;
	      TmpKeyInvertTable = keyInvertTable[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
	      TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
	      TmpKeyInvertTable[0] = keys[i];
	      TmpKeyInvertTableNbrIndices[0] = 1;
	    }
	  else
	    {
	      int Lim = keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
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
		  TmpKeyInvertTable[keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax]] = TmpKey;
		  TmpKeyInvertTableNbrIndices[keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax]] = 1;
		  ++keyInvertSectorSize[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
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
  CurrentLzSzMax = (lzMax<<1)+1;
  CurrentNbrLzSzMax = 1;
  TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
  int** TmpKeyInvertIndices = keyInvertIndices[CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax];
  int TmpPos2;
  int TmpPos3;
  int TmpPos4 = CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax;
  TmpSize = keyInvertSectorSize[TmpPos4];
  TmpKeyInvertIndices = keyInvertIndices[TmpPos4];
  TmpKeyInvertTable = keyInvertTable[TmpPos4];
  TmpKeyInvertTableNbrIndices = keyInvertTableNbrIndices[TmpPos4];
  for (int i = 0; i < dimension; ++i)
    {
      int TmpKey = keys[i];
      if (CurrentLzSzMax != stateLzSzMax[i])
	{
	  CurrentLzSzMax = stateLzSzMax[i];
	  if (CurrentLzSzMax&1)
	    CurrentNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])>>16;
	  else
	    CurrentNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])&0xffff;
	  TmpPos4 = CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax;
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
	  if (stateLzSzMax[i]&1)
	    NewNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])>>16;
	  else
	    NewNbrLzSzMax = (stateDescription[i][CurrentLzSzMax>>1])&0xffff;	  
	  if (NewNbrLzSzMax != CurrentNbrLzSzMax)
	    {
	      CurrentNbrLzSzMax = NewNbrLzSzMax;
	      TmpPos4 = CurrentLzSzMax * this->IncMaxNbrBosons + CurrentNbrLzSzMax;
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
}

// generate look-up table associated to current Hilbert space
// 
// stateDescription = array describing the state
// lzszmax = maximum Lz/Sz value lzMax<<1 (spin down) / lzMax<<1 + 1 (spin up) reached by a boson in the state
// return value = key associated to the state

int BosonOnSphereWithSpin::GenerateKey(int* stateDescription, int lzszmax)
{
  int Key = 0;
  lzszmax>>=1;
  //cout << "GenerateKey: lzszmax="<<lzszmax<<endl;
  for (int i = 0; i <= lzszmax; ++i)
    {
      Key += this->KeyMultiplicationTable[i<<1] * (stateDescription[i] & 0xffff) +
	this->KeyMultiplicationTable[(i<<1)+1] * (stateDescription[i] >>16);
      //cout << "i="<<i<<"Key="<<Key<<" KeyMultiplicationTable[i<<1]="<<KeyMultiplicationTable[i<<1]<<" KeyMultiplicationTable[(i<<1)+1]="<<KeyMultiplicationTable[(i<<1)+1]<<endl;
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
}


// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex BosonOnSphereWithSpin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
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

Complex BosonOnSphereWithSpin::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
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

Complex BosonOnSphereWithSpin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
					     int firstComponent, int nbrComponent)
{  
  Complex Value(0.0, 0.0), TmpValue(0.0, 0.0);
  Complex Tmp;
  ComplexMatrix PermUp(this->NbrBosonsUp, this->NbrBosonsUp);  
  ComplexMatrix FunctionsUp(this->LzMax + 1, this->NbrBosonsUp);
  ComplexMatrix PermDown(this->NbrBosonsDown, this->NbrBosonsDown);  
  ComplexMatrix FunctionsDown(this->LzMax + 1, this->NbrBosonsDown);
  RealVector TmpCoordinates(2);
  int* Indices = new int [this->NbrBosons];
  int Pos;
  int Lz;
  for (int j = 0; j < this->NbrBosonsUp; ++j)
    {
      TmpCoordinates[0] = position[j << 1];
      TmpCoordinates[1] = position[1 + (j << 1)];
      for (int i = 0; i <= this->LzMax; ++i)
	{	  
	  basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	  FunctionsUp[j].Re(i) = Tmp.Re;
	  FunctionsUp[j].Im(i) = Tmp.Im;
	}
    }
  for (int j = this->NbrBosonsUp; j < this->NbrBosons; ++j)
    {
      TmpCoordinates[0] = position[j << 1];
      TmpCoordinates[1] = position[1 + (j << 1)];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	  FunctionsDown[j-this->NbrBosonsUp].Re(i) = Tmp.Re;
	  FunctionsDown[j-this->NbrBosonsUp].Im(i) = Tmp.Im;
	}
    }
  double* Factors = new double [this->IncMaxNbrBosons];  
  Factors[0] = 1.0;
  Factors[1] = 1.0;
  for (int i = 2; i < this->IncMaxNbrBosons; ++i)
    Factors[i] = Factors[i - 1] / sqrt((double) i);
  double TmpFactor, TmpComponent;
  int TmpStateDescription;
  int LastComponent = firstComponent + nbrComponent;  
  int* ChangeBitSignUp;
  int* ChangeBitUp;
  PermUp.EvaluateFastPermanentPrecalculationArray(ChangeBitUp, ChangeBitSignUp);
  int* ChangeBitSignDown;
  int* ChangeBitDown;
  PermDown.EvaluateFastPermanentPrecalculationArray(ChangeBitDown, ChangeBitSignDown);  
  for (int k = firstComponent; k < LastComponent; ++k)
    {
      TmpComponent = state[k];
      if (TmpComponent!=0.0)
	{
	  TmpFactor = Factors[this->NbrBosonsUp];
	  Pos = 0;
	  Lz = 0;
	  while (Pos < this->NbrBosonsUp)
	    {	      
	      TmpStateDescription = (this->StateDescription[k][Lz])>>16;
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
	  for (int i = 0; i < this->NbrBosonsUp; ++i)
	    {
	      ComplexVector& TmpColum2 = FunctionsUp[i];	  
	      for (int j = 0; j < this->NbrBosonsUp; ++j)
		{
		  PermUp[i].Re(j) = TmpColum2.Re(Indices[j]);
		  PermUp[i].Im(j) = TmpColum2.Im(Indices[j]);
		}
	    }
	  TmpValue = PermUp.FastPermanent(ChangeBitUp, ChangeBitSignUp) * TmpFactor;
	  Pos = 0;
	  Lz = 0;
	  TmpFactor = Factors[this->NbrBosonsDown];
	  while (Pos < this->NbrBosonsDown)
	    {
	      TmpStateDescription = (this->StateDescription[k][Lz])&0xffff;
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
	  for (int i = 0; i < this->NbrBosonsDown; ++i)
	    {
	      ComplexVector& TmpColum2 = FunctionsDown[i];	  
	      for (int j = 0; j < this->NbrBosonsDown; ++j)
		{
		  PermDown[i].Re(j) = TmpColum2.Re(Indices[j]);
		  PermDown[i].Im(j) = TmpColum2.Im(Indices[j]);
		}
	    }
	  TmpValue *= PermDown.FastPermanent(ChangeBitDown, ChangeBitSignDown) * TmpFactor;
	  Value+=TmpComponent*TmpValue;
	}
    }
  delete[] ChangeBitSignUp;
  delete[] ChangeBitUp;
  delete[] ChangeBitSignDown;
  delete[] ChangeBitDown;
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

Complex BosonOnSphereWithSpin::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
							      int nextCoordinates, int firstComponent, int nbrComponent)
{
  cout << "Attention, TimeCoherence not fully implemented, yet!"<<endl;
  double* Factors = new double [this->IncMaxNbrBosons];
  Factors[0] = 1.0;
  Factors[1] = 1.0;
  for (int i = 2; i <= this->IncMaxNbrBosons; ++i)
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
      ComplexMatrix PermUp(this->NbrBosonsUp, this->NbrBosonsUp);  
      ComplexMatrix FunctionsUp(this->LzMax + 1, this->NbrBosonsUp);
      ComplexMatrix PermDown(this->NbrBosonsDown, this->NbrBosonsDown);  
      ComplexMatrix FunctionsDown(this->LzMax + 1, this->NbrBosonsDown);
      RealVector TmpCoordinates(2);
      int* Indices = new int [this->NbrBosons];
      int Pos;
      int Lz;
      for (int j = 0; j < this->NbrBosonsUp; ++j)
	{
	  TmpCoordinates[0] = position[j << 1];
	  TmpCoordinates[1] = position[1 + (j << 1)];
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	      FunctionsUp[j].Re(i) = Tmp.Re;
	      FunctionsUp[j].Im(i) = Tmp.Im;
	    }
	}
      for (int j = this->NbrBosonsUp; j < this->NbrBosons; ++j)
	{
	  TmpCoordinates[0] = position[j << 1];
	  TmpCoordinates[1] = position[1 + (j << 1)];
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	      FunctionsDown[j-this->NbrBosonsUp].Re(i) = Tmp.Re;
	      FunctionsDown[j-this->NbrBosonsUp].Im(i) = Tmp.Im;
	    }
	}
      int* ChangeBitSignUp;
      int* ChangeBitUp;
      PermUp.EvaluateFastPermanentPrecalculationArray(ChangeBitUp, ChangeBitSignUp);
      int* ChangeBitSignDown;
      int* ChangeBitDown;
      PermDown.EvaluateFastPermanentPrecalculationArray(ChangeBitDown, ChangeBitSignDown);  
      int TmpStateDescription;
      int LastComponent = firstComponent + nbrComponent;
      for (int k = firstComponent; k < LastComponent; ++k)
	{
	  TmpFactor = state[k] * Factors[this->NbrBosonsUp]* Factors[this->NbrBosonsDown];

	  if (TmpFactor!=0.0)
	    {
	      Pos = 0;
	      Lz = 0;
	      while (Pos < this->NbrBosonsUp)
		{	      
		  TmpStateDescription = (this->StateDescription[k][Lz])>>16;
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
	      for (int i = 0; i < this->NbrBosonsUp; ++i)
		{
		  ComplexVector& TmpColum2 = FunctionsUp[i];	  
		  for (int j = 0; j < this->NbrBosonsUp; ++j)
		    {
		      PermUp[i].Re(j) = TmpColum2.Re(Indices[j]);
		      PermUp[i].Im(j) = TmpColum2.Im(Indices[j]);
		    }
		}
	      
	      Pos = 0;
	      Lz = 0;
	      while (Pos < this->NbrBosonsDown)
		{
		  TmpStateDescription = (this->StateDescription[k][Lz])&0xffff;
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
	      for (int i = 0; i < this->NbrBosonsDown; ++i)
		{
		  ComplexVector& TmpColum2 = FunctionsDown[i];	  
		  for (int j = 0; j < this->NbrBosonsDown; ++j)
		    {
		      PermDown[i].Re(j) = TmpColum2.Re(Indices[j]);
		      PermDown[i].Im(j) = TmpColum2.Im(Indices[j]);
		    }
		}
	    }
	  // attention, need to decide how time coherence works, and edit from here...
	  if (this->MinorsUp[k] == 0)
	    {
	      this->MinorsUp[k] = new Complex [this->NbrBosonsUp];
	    }
	  PermUp.FastPermanentMinorDevelopment(ChangeBitUp, ChangeBitSignUp, nextCoordinates, this->MinorsUp[k]);
	  TmpPerm = 0.0;
	  for (int i = 0; i < this->NbrBosons; ++i)
	    TmpPerm += this->MinorsUp[k][i] * Complex (PermUp[nextCoordinates].Re(i), 
						     PermUp[nextCoordinates].Im(i));
	  Value += TmpPerm * TmpFactor;
	}
      delete[] ChangeBitSignUp;
      delete[] ChangeBitUp;
      delete[] ChangeBitSignDown;
      delete[] ChangeBitDown;
      (*(this->KeptCoordinates)) = nextCoordinates;
    }
  else
    {
      int StartCoordinates, EndCoordinates, EffectiveNbrBosons;
      if (nextCoordinates<NbrBosonsUp)
	{
	  StartCoordinates=0;
	  EndCoordinates=NbrBosonsUp;
	  EffectiveNbrBosons=NbrBosonsUp;
	}
      else
	{
	  StartCoordinates=NbrBosonsUp;
	  EndCoordinates=NbrBosons;
	  EffectiveNbrBosons=NbrBosonsDown;
	}
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
	  Complex* TmpMinors = this->MinorsUp[k];
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

void BosonOnSphereWithSpin::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
  if ((timeCoherence == true) && (this->MinorsUp == 0))
    {
      this->MinorsUp = new Complex* [this->HilbertSpaceDimension];
      this->MinorsDown = new Complex* [this->HilbertSpaceDimension];
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  this->MinorsUp[i] = 0;
	  this->MinorsDown[i] = 0;
	}
    }
}
