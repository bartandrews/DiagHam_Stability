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
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/BosonOnSphereFullShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Matrix/SparseRealMatrix.h"


#include <math.h>
#include <algorithm>
#include <set>
#include <limits>


using std::cout;
using std::endl;
using std::hex;
using std::dec;

// print some test outputs
//#define TEST_BOSONONTORUS_SHORT


// default constructor
// 

BosonOnTorusShort::BosonOnTorusShort ()
{
  this->MaxOccupation = -1;
  this->LookUpTable = 0;
  this->LargeLookUpTable = 0;
}

// basic constructor
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson

BosonOnTorusShort::BosonOnTorusShort (int nbrBosons, int maxMomentum)
{
  this->TargetSpace = this;
  this->NbrBosons = nbrBosons;
  this->MaxOccupation = this->NbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->KyMax = maxMomentum;
  this->NbrKyValue = this->KyMax + 1;
  this->TotalKyFlag = false;
  this->TemporaryState = new unsigned long [this->KyMax + 1];
  this->ProdATemporaryState = new unsigned long [this->KyMax + 1];

  this->MomentumModulo = FindGCD(this->NbrBosons, this->KyMax);
  this->KyMomentumModulo = this->KyMax;
  this->StateShift = this->KyMax / this->MomentumModulo;
  this->LastMomentumMask = 0x1ul << (this->KyMax + this->NbrBosons - 1);

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->KyMax - 1, this->KyMax - 1, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateKyMax = new int [this->LargeHilbertSpaceDimension];
  long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->KyMax - 1, this->KyMax - 1, 0);  
  if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << "error while generating the Hilbert space " << TmpLargeHilbertSpaceDimension << ", should be " 
	   <<  this->LargeHilbertSpaceDimension << endl;
    }
  
  this->GenerateLookUpTable(0);
#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
  UsedMemory += this->NbrKyValue * sizeof(int);
  UsedMemory += this->NbrKyValue * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
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

BosonOnTorusShort::BosonOnTorusShort (int nbrBosons, int maxMomentum, int momentumConstraint)
{
  this->TargetSpace = this;
  this->NbrBosons = nbrBosons;
  this->MaxOccupation = this->NbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->KyMax = maxMomentum;
  this->NbrKyValue = this->KyMax + 1;
  this->KyMomentumModulo = this->KyMax;
  this->TotalKy = momentumConstraint % this->KyMomentumModulo;
  this->TotalKyFlag = true;

  this->MomentumModulo = FindGCD(this->NbrBosons, this->KyMax);
  this->StateShift = this->KyMax / this->MomentumModulo;
  this->LastMomentumMask = 0x1ul << (this->KyMax + this->NbrBosons - 1);

  this->TemporaryState = new unsigned long [this->KyMax + 1];
  this->ProdATemporaryState = new unsigned long [this->KyMax + 1];

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->KyMax);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateKyMax = new int [this->LargeHilbertSpaceDimension];
  long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->KyMax - 1, this->KyMax - 1, 0, 0);
  if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << "error while generating the Hilbert space " << TmpLargeHilbertSpaceDimension << ", should be " 
	   <<  this->LargeHilbertSpaceDimension << endl;
    }
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
  UsedMemory += this->NbrKyValue * sizeof(int);
  UsedMemory += this->NbrKyValue * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
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

// constructor with a constraint of the total momentum of states and the maximum occupation per orbital
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// momentumConstraint = index of the momentum orbit
// maxOccupation = maximum bosonic occupation per orbital

BosonOnTorusShort::BosonOnTorusShort (int nbrBosons, int maxMomentum, int momentumConstraint, int maxOccupation)
{
  this->TargetSpace = this;
  this->NbrBosons = nbrBosons;
  this->MaxOccupation = maxOccupation;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->KyMax = maxMomentum;
  this->NbrKyValue = this->KyMax + 1;
  this->KyMomentumModulo = this->KyMax;
  this->TotalKy = momentumConstraint % this->KyMomentumModulo;
  this->TotalKyFlag = true;

  this->MomentumModulo = FindGCD(this->NbrBosons, this->KyMax);
  this->StateShift = this->KyMax / this->MomentumModulo;
  this->LastMomentumMask = 0x1ul << (this->KyMax + this->NbrBosons - 1);

  this->TemporaryState = new unsigned long [this->KyMax + 1];
  this->ProdATemporaryState = new unsigned long [this->KyMax + 1];

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimensionTruncatedOccupation(this->NbrBosons, this->KyMax - 1, this->KyMax - 1, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateKyMax = new int [this->LargeHilbertSpaceDimension];
  long TmpLargeHilbertSpaceDimension = this->GenerateStatesTruncatedOccupation(this->NbrBosons, this->KyMax - 1, this->KyMax - 1, 0, 0);
  if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << "error while generating the Hilbert space " << TmpLargeHilbertSpaceDimension << ", should be " 
	   <<  this->LargeHilbertSpaceDimension << endl;
    }
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
  UsedMemory += this->NbrKyValue * sizeof(int);
  UsedMemory += this->NbrKyValue * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
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

BosonOnTorusShort::BosonOnTorusShort(const BosonOnTorusShort& bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->KyMax = bosons.KyMax;
  this->NbrKyValue = bosons.NbrKyValue;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->MomentumModulo = bosons.MomentumModulo;
  this->StateKyMax = bosons.StateKyMax;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMask = bosons.LastMomentumMask;
  this->Flag = bosons.Flag;
  this->TemporaryState = new unsigned long [this->KyMax + 1];
  this->ProdATemporaryState = new unsigned long [this->KyMax + 1];
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->LargeLookUpTable = bosons.LargeLookUpTable;
}

// destructor
//

BosonOnTorusShort::~BosonOnTorusShort ()
{
  if ((this->LargeHilbertSpaceDimension != 0l) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateKyMax;
      if (this->LookUpTable != 0)
	{
	  for (int i=0;i<this->KyMax + this->NbrBosons; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;
	}
      else
	{
	  for (int i=0;i<this->KyMax + this->NbrBosons; ++i)
	    delete[] this->LargeLookUpTable[i];
	  delete[] this->LargeLookUpTable;
	}
      delete[] this->LookUpTableShift;
    }
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorusShort& BosonOnTorusShort::operator = (const BosonOnTorusShort& bosons)
{
  if ((this->LargeHilbertSpaceDimension != 0l) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateKyMax;
      for (int i=0;i<this->KyMax + this->NbrBosons; ++i)
	delete [] this->LookUpTable[i];
      delete[] this->LookUpTable;
      delete[] this->LookUpTableShift;
    }
  if (this->TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->KyMax = bosons.KyMax;
  this->NbrKyValue = bosons.NbrKyValue;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateKyMax = bosons.StateKyMax;
  this->MomentumModulo = bosons.MomentumModulo;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMask = bosons.LastMomentumMask;
  this->TemporaryState = new unsigned long [this->KyMax + 1];
  this->ProdATemporaryState = new unsigned long [this->KyMax + 1];
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->LargeLookUpTable = bosons.LargeLookUpTable;
  this->Flag = bosons.Flag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnTorusShort::Clone()
{
  return new BosonOnTorusShort(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnTorusShort::SetTargetSpace(ParticleOnTorus* targetSpace)
{
  this->TargetSpace = (BosonOnTorusShort*) targetSpace;
}


// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int BosonOnTorusShort::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}


// return number of particles of the target space
//
// return value = Hilbert space dimension

int BosonOnTorusShort::GetTargetNbrParticles()
{
  return this->TargetSpace->GetNbrParticles();
}

// get momemtum value of a given state
//
// index = state index
// return value = state momentum

int BosonOnTorusShort::GetMomentumValue(int index)
{
  this->FermionToBoson(this->StateDescription[index], this->StateKyMax[index]  + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
  int Momentum = 0;
  for (int i = 0; i <= this->TemporaryStateKyMax; ++i)
    {
      Momentum += (this->TemporaryState[i] * i);
    }
  return (Momentum % this->KyMomentumModulo);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnTorusShort::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->TotalKyFlag == false)
    {
      for (int i = 0; i < this->KyMax; ++i)
	TmpList += new PeriodicMomentumQuantumNumber (i, this->KyMax);
    }
  else
    TmpList += new PeriodicMomentumQuantumNumber (this->TotalKy, this->KyMax);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnTorusShort::GetQuantumNumber (int index)
{
  if (this->TotalKyFlag == false)
    {
      return  new PeriodicMomentumQuantumNumber (this->GetMomentumValue(index), this->KyMax);
    }
  else
    return new PeriodicMomentumQuantumNumber (this->TotalKy, this->KyMax);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnTorusShort::ExtractSubspace (AbstractQuantumNumber& q, 
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

int BosonOnTorusShort::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int CurrentLzMax = this->StateKyMax[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }    
  this->FermionToBoson(this->StateDescription[index], CurrentLzMax + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || (this->TemporaryState[n1] == 0) || (this->TemporaryState[n2] == 0) || ((n1 == n2) && (this->TemporaryState[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  this->TemporaryStateKyMax = CurrentLzMax;
  if (this->TemporaryStateKyMax < m1)
    this->TemporaryStateKyMax = m1;
  if (this->TemporaryStateKyMax < m2)
    this->TemporaryStateKyMax = m2;
  for (int i = CurrentLzMax + 1; i <= this->TemporaryStateKyMax; ++i)
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
  while (this->TemporaryState[this->TemporaryStateKyMax] == 0x0ul)
    --this->TemporaryStateKyMax;
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
}

// return matrix representation of the annihilation operator a_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& BosonOnTorusShort::A (int i, Matrix& M)
{
  return M;
}

// return matrix representation of the creation operator a^+_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& BosonOnTorusShort::Ad (int i, Matrix& M)
{
  return M;
}

// apply an annihilation operator a_i and return the index in the target space
//
// index = state index
// n = index of annihilation operator
// coefficient = will be multiplied by the prefactor of the bosonic ladder operator
// return value = index in the target space
int BosonOnTorusShort::A (int index, int n, double &coefficient)
{
  int CurrentLzMax = this->StateKyMax[index];
  if (n > CurrentLzMax)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }    
  this->FermionToBoson(this->StateDescription[index], CurrentLzMax + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
  if (this->TemporaryState[n] == 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  this->TemporaryStateKyMax = CurrentLzMax;
  for (int i = CurrentLzMax + 1; i <= this->TemporaryStateKyMax; ++i)
    this->TemporaryState[i] = 0;
  coefficient *= sqrt(this->TemporaryState[n]);
  --this->TemporaryState[n];
  while (this->TemporaryState[this->TemporaryStateKyMax] == 0x0ul)
    --this->TemporaryStateKyMax;
  return this->TargetSpace->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
}

// apply a creation operator a_i and return the index in the target space
//
// index = state index
// m = index of creation operator
// coefficient = will be multiplied by the prefactor of the bosonic ladder operator
// return value = index in the target space
int BosonOnTorusShort::Ad (int index, int m, double &coefficient)
{
  int CurrentLzMax = this->StateKyMax[index];
  this->FermionToBoson(this->StateDescription[index], CurrentLzMax + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
  this->TemporaryStateKyMax = CurrentLzMax;
  if (this->TemporaryStateKyMax < m)
    this->TemporaryStateKyMax = m;
  for (int i = CurrentLzMax + 1; i <= this->TemporaryStateKyMax; ++i)
    this->TemporaryState[i] = 0;
  ++this->TemporaryState[m];
  coefficient *= sqrt(this->TemporaryState[m]);
  return this->TargetSpace->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
}


// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnTorusShort::ProdA (int index, int* n, int nbrIndices)
{
  this->FermionToBoson(this->StateDescription[index], this->StateKyMax[index]  + this->NbrBosons - 1, this->ProdATemporaryState, this->ProdATemporaryStateKyMax);
  int TmpCoefficient = 1;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (n[i] > this->ProdATemporaryStateKyMax)
	return 0.0;
      unsigned long& Tmp = this->ProdATemporaryState[n[i]];
      if (Tmp == 0)
	return 0.0;
      TmpCoefficient *= Tmp;
      --Tmp;
    }
  for (int i = this->ProdATemporaryStateKyMax + 1; i < this->NbrKyValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt((double) TmpCoefficient);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusShort::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  this->TemporaryStateKyMax=this->ProdATemporaryStateKyMax;
  for (int i = 0; i < this->NbrKyValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  
  int TmpCoefficient = 1;
  for (int i = 0; i < nbrIndices; ++i)
    {
      if(m[i]>this->TemporaryStateKyMax)
	this->TemporaryStateKyMax=m[i];
      TmpCoefficient *= ++this->TemporaryState[m[i]];
    }
  coefficient = sqrt((double) TmpCoefficient);
  while (this->TemporaryState[this->TemporaryStateKyMax] == 0)
    --this->TemporaryStateKyMax;
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnTorusShort::AdA (int index, int m)
{
  this->FermionToBoson(this->StateDescription[index], this->StateKyMax[index]  + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
  if (this->TemporaryStateKyMax < m)  
    return 0.0;
  return (double) (this->TemporaryState[m]);  
}


// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusShort::AdA (int index, int m, int n, double& coefficient)
{
  int CurrentLzMax = this->StateKyMax[index];
  if (n > CurrentLzMax)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }    
  this->FermionToBoson(this->StateDescription[index], CurrentLzMax + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
  if ((n > CurrentLzMax) || (this->TemporaryState[n] == 0))
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  this->TemporaryStateKyMax = CurrentLzMax;
  if (this->TemporaryStateKyMax < m)
    this->TemporaryStateKyMax = m;
  for (int i = CurrentLzMax + 1; i <= this->TemporaryStateKyMax; ++i)
    this->TemporaryState[i] = 0;
  coefficient = this->TemporaryState[n];
  --this->TemporaryState[n];
  ++this->TemporaryState[m];
  coefficient *= this->TemporaryState[m];
  coefficient = sqrt(coefficient);
  while (this->TemporaryState[this->TemporaryStateKyMax] == 0x0ul)
    --this->TemporaryStateKyMax;
  return this->TargetSpace->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusShort::AdAd (int m1, int m2, double& coefficient)
{
  for (int i = 0; i < this->NbrKyValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m2];
  coefficient = this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  this->TemporaryStateKyMax = this->KyMax;
  while (this->TemporaryState[this->TemporaryStateKyMax] == 0)
    --this->TemporaryStateKyMax;

#ifdef DEBUG_BosonOnTorusShort
  cout << "AdAd("<<m1<<", "<<m2<<") result: ";
  int i = 0;
  for (; i <= this->TemporaryStateKyMax; ++i)
    cout << this->TemporaryState[i] << " ";
  for (; i < this->NbrKyValue; ++i)
    cout << "0 ";
  cout << "(KyMax="<< this->TemporaryStateKyMax<<endl;
#endif
  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnTorusShort::AA (int index, int n1, int n2)
{
#ifdef DEBUG_BosonOnTorusShort
  cout << "AA("<<n1<<", "<<n2<<") on ";
  this->PrintState(cout,index);
#endif
  int CurrentLzMax = this->StateKyMax[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription[index], CurrentLzMax + this->NbrBosons - 1, this->ProdATemporaryState, this->ProdATemporaryStateKyMax);
#ifdef DEBUG_BosonOnTorusShort
  cout << " (KyMax=" << this->ProdATemporaryStateKyMax<<" vs "<<CurrentLzMax<<")"<<endl;
#endif
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || 
      (this->ProdATemporaryState[n1] == 0) || (this->ProdATemporaryState[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryState[n1] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryState[n2];
  --this->ProdATemporaryState[n2];
  Coefficient *= this->ProdATemporaryState[n1];
  --this->ProdATemporaryState[n1];
  for (int i = this->ProdATemporaryStateKyMax + 1; i < this->NbrKyValue; ++i)
    this->ProdATemporaryState[i] = 0ul;

#ifdef DEBUG_BosonOnTorusShort  
  int i = 0;
  for (; i <= this->ProdATemporaryStateKyMax; ++i)
    cout << this->ProdATemporaryState[i] << " ";
  for (; i < this->NbrKyValue; ++i)
    cout << "0 ";
  cout << "(KyMax="<< this->ProdATemporaryStateKyMax<<")"<<endl;
#endif
  
  return sqrt(Coefficient);
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusShort::PrintStateMonomial (ostream& Str, long state)
{
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  cout << this->StateDescription[state] << endl;
  this->ConvertToMonomial(this->StateDescription[state], this->StateKyMax[state] + this->NbrBosons - 1, TmpMonomial);
  Str << "[";
  if (TmpMonomial[0] != 0)
    Str << TmpMonomial[0];
  for (int i = 1; (i < this->NbrBosons) && (TmpMonomial[i] > 0); ++i)
    Str << "," << TmpMonomial[i];
  Str << "]";
  delete[] TmpMonomial;
  return Str;
}


// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnTorusShort::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    return PosMin;
}

// find state index when the Hilbert space is larger than 2^31
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

long BosonOnTorusShort::FindStateLargeIndex(unsigned long stateDescription, int lzmax)
{
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LargeLookUpTable[lzmax][PosMax];
  PosMax = this->LargeLookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    return PosMin;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusShort::PrintState (ostream& Str, int state)
{
  int Max = this->StateKyMax[state];
  this->FermionToBoson(this->StateDescription[state], Max + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
  int i = 0;
  for (; i <= this->TemporaryStateKyMax; ++i)
    Str << this->TemporaryState[i] << " ";
  for (; i < this->KyMax; ++i)
    Str << "0 ";
 return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a fermion in the state
// currentKyMax = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnTorusShort::GenerateStates(int nbrBosons, int maxMomentum, int currentKyMax, long pos)
{
  if ((nbrBosons == 0) || (currentKyMax < 0))
    {
      return pos;
    }
  if (currentKyMax == 0)
    {
      this->StateDescription[pos] = 0ul;      
      this->StateDescription[pos] = ((1ul << (nbrBosons + 1)) - 1ul);
      this->StateKyMax[pos] = maxMomentum;
      return pos + 1l;
    }

  this->StateDescription[pos] = 0ul;
  this->StateKyMax[pos] = maxMomentum;
  ++pos;

  int TmpNbrBosons = 1;
  int ReducedCurrentKyMax = currentKyMax - 1;
  long TmpPos = pos;
  while (TmpNbrBosons < nbrBosons)
    {
      unsigned long Mask = ((0x1ul << (nbrBosons - TmpNbrBosons)) - 0x1ul) << (currentKyMax + TmpNbrBosons - 1);
      TmpPos = this->GenerateStates(TmpNbrBosons, maxMomentum, ReducedCurrentKyMax, pos);
      for (; pos <TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
      ++TmpNbrBosons;
    }
  if (maxMomentum == currentKyMax)
    return this->GenerateStates(nbrBosons, ReducedCurrentKyMax, ReducedCurrentKyMax, pos);
  else
    return this->GenerateStates(nbrBosons, maxMomentum, ReducedCurrentKyMax, pos);
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson in the state
// currentKyMax = momentum maximum value for bosons that are still to be placed
// pos = position in StateDescription array where to store states
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

long BosonOnTorusShort::GenerateStates(int nbrBosons, int maxMomentum, int currentKyMax, long pos, int currentMomentum)
{
  if (nbrBosons == 0)
    {
      if ((currentMomentum % this->KyMomentumModulo) == this->TotalKy)
	{
	  this->StateDescription[pos] = 0x0ul;
	  this->StateKyMax[pos] = maxMomentum;
	  return pos + 1l;
	}
      else
	{
	  return pos;
	}
    }
  if (currentKyMax == 0)
    {
      if ((currentMomentum % this->KyMomentumModulo) == this->TotalKy)
	{
	  this->StateDescription[pos] = (0x1ul << nbrBosons) - 0x1ul;
	  this->StateKyMax[pos] = maxMomentum;
	  return pos + 1l;
	}
      else
	return pos;
    }

  int TmpNbrBosons = 0;
  int ReducedCurrentKyMax = currentKyMax - 1;
  long TmpPos = pos;
  while (TmpNbrBosons < nbrBosons)
    {
      TmpPos = this->GenerateStates(TmpNbrBosons, maxMomentum, ReducedCurrentKyMax, pos, currentMomentum + (nbrBosons - TmpNbrBosons) * currentKyMax);
      unsigned long Mask = ((0x1ul << (nbrBosons - TmpNbrBosons)) - 0x1ul) << (currentKyMax + TmpNbrBosons);

      for (; pos <TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
      ++TmpNbrBosons;
    }
  if (maxMomentum == currentKyMax)
    return this->GenerateStates(nbrBosons, ReducedCurrentKyMax, ReducedCurrentKyMax, pos, currentMomentum);
  else
    return this->GenerateStates(nbrBosons, maxMomentum, ReducedCurrentKyMax, pos, currentMomentum);
}

// generate all states corresponding to the constraints, including a truncation on the occupation numbers
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson in the state
// currentKyMax = momentum maximum value for bosons that are still to be placed
// pos = position in StateDescription array where to store states
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

long BosonOnTorusShort::GenerateStatesTruncatedOccupation(int nbrBosons, int maxMomentum, int currentKyMax, long pos, int currentMomentum)
{
  if (nbrBosons == 0)
    {
      if ((currentMomentum % this->KyMomentumModulo) == this->TotalKy)
	{
	  this->StateDescription[pos] = 0x0ul;
	  this->StateKyMax[pos] = maxMomentum;
	  return pos + 1l;
	}
      else
	{
	  return pos;
	}
    }
  if (currentKyMax == 0)
    {
      if (((currentMomentum % this->KyMomentumModulo) == this->TotalKy) && (nbrBosons <= this->MaxOccupation))
	{
	  this->StateDescription[pos] = (0x1ul << nbrBosons) - 0x1ul;
	  this->StateKyMax[pos] = maxMomentum;
	  return pos + 1l;
	}
      else
	return pos;
    }

  int TmpNbrBosons = nbrBosons -  this->MaxOccupation;
  if (TmpNbrBosons < 0)
    TmpNbrBosons = 0;
  int ReducedCurrentKyMax = currentKyMax - 1;
  long TmpPos = pos;
  while (TmpNbrBosons < nbrBosons)
    {
      TmpPos = this->GenerateStatesTruncatedOccupation(TmpNbrBosons, maxMomentum, ReducedCurrentKyMax, pos, currentMomentum + (nbrBosons - TmpNbrBosons) * currentKyMax);
      unsigned long Mask = ((0x1ul << (nbrBosons - TmpNbrBosons)) - 0x1ul) << (currentKyMax + TmpNbrBosons);

      for (; pos <TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
      ++TmpNbrBosons;
    }
  if (maxMomentum == currentKyMax)
    return this->GenerateStatesTruncatedOccupation(nbrBosons, ReducedCurrentKyMax, ReducedCurrentKyMax, pos, currentMomentum);
  else
    return this->GenerateStatesTruncatedOccupation(nbrBosons, maxMomentum, ReducedCurrentKyMax, pos, currentMomentum);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnTorusShort::GenerateLookUpTable(int memory)
{
  // evaluate look-up table size
  int TmpNbrKyValue =  this->KyMax + this->NbrBosons;
  memory /= (sizeof(int*) * TmpNbrKyValue);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > TmpNbrKyValue)
    this->MaximumLookUpShift = TmpNbrKyValue;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = 0;
  this->LargeLookUpTable = 0;
  this->LookUpTableShift = new int [TmpNbrKyValue];
  int CurrentLzMax = this->KyMax + this->NbrBosons - 1;
  while ((this->StateDescription[0] >> CurrentLzMax) == 0x0ul)
    --CurrentLzMax;
  
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    {
      this->LargeLookUpTable = new long* [TmpNbrKyValue];
      for (int i = 0; i < TmpNbrKyValue; ++i)
	this->LargeLookUpTable[i] = new long [this->LookUpTableMemorySize + 1];
      long* TmpLookUpTable = this->LargeLookUpTable[CurrentLzMax];
      if (CurrentLzMax < this->MaximumLookUpShift)
	this->LookUpTableShift[CurrentLzMax] = 0;
      else
	this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
      int CurrentShift = this->LookUpTableShift[CurrentLzMax];
      unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
      unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
      while (CurrentLookUpTableValue > TmpLookUpTableValue)
	{
	  TmpLookUpTable[CurrentLookUpTableValue] = 0;
	  --CurrentLookUpTableValue;
	}
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->StateDescription[i];
	  int TmpLzMax = CurrentLzMax;
	  while ((TmpState >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  if (CurrentLzMax != TmpLzMax)
	    {
	      while (CurrentLookUpTableValue > 0)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[0] = i;
	      CurrentLzMax = TmpLzMax;
	      TmpLookUpTable = this->LargeLookUpTable[CurrentLzMax];
	      if (CurrentLzMax < this->MaximumLookUpShift)
		this->LookUpTableShift[CurrentLzMax] = 0;
	      else
		this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
	      CurrentShift = this->LookUpTableShift[CurrentLzMax];
	      TmpLookUpTableValue = TmpState >> CurrentShift;
	      CurrentLookUpTableValue = this->LookUpTableMemorySize;
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	  else
	    {
	      TmpLookUpTableValue = TmpState >> CurrentShift;
	      if (TmpLookUpTableValue != CurrentLookUpTableValue)
		{
		  while (CurrentLookUpTableValue > TmpLookUpTableValue)
		    {
		      TmpLookUpTable[CurrentLookUpTableValue] = i;
		      --CurrentLookUpTableValue;
		    }
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		}
	    }
	}
      while (CurrentLookUpTableValue > 0)
	{
	  TmpLookUpTable[CurrentLookUpTableValue] = this->LargeHilbertSpaceDimension - 1;
	  --CurrentLookUpTableValue;
	}
      TmpLookUpTable[0] = this->LargeHilbertSpaceDimension - 1;
   }
  else
    {
      this->LookUpTable = new int* [TmpNbrKyValue];
      for (int i = 0; i < TmpNbrKyValue; ++i)
	this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
      
      int* TmpLookUpTable = this->LookUpTable[CurrentLzMax];
      if (CurrentLzMax < this->MaximumLookUpShift)
	this->LookUpTableShift[CurrentLzMax] = 0;
      else
	this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
      int CurrentShift = this->LookUpTableShift[CurrentLzMax];
      unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
      unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
      while (CurrentLookUpTableValue > TmpLookUpTableValue)
	{
	  TmpLookUpTable[CurrentLookUpTableValue] = 0;
	  --CurrentLookUpTableValue;
	}
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->StateDescription[i];
	  int TmpLzMax = CurrentLzMax;
	  while ((TmpState >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  if (CurrentLzMax != TmpLzMax)
	    {
	      while (CurrentLookUpTableValue > 0)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[0] = i;
	      CurrentLzMax = TmpLzMax;
	      TmpLookUpTable = this->LookUpTable[CurrentLzMax];
	      if (CurrentLzMax < this->MaximumLookUpShift)
		this->LookUpTableShift[CurrentLzMax] = 0;
	      else
		this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
	      CurrentShift = this->LookUpTableShift[CurrentLzMax];
	      TmpLookUpTableValue = TmpState >> CurrentShift;
	      CurrentLookUpTableValue = this->LookUpTableMemorySize;
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	  else
	    {
	      TmpLookUpTableValue = TmpState >> CurrentShift;
	      if (TmpLookUpTableValue != CurrentLookUpTableValue)
		{
		  while (CurrentLookUpTableValue > TmpLookUpTableValue)
		    {
		      TmpLookUpTable[CurrentLookUpTableValue] = i;
		      --CurrentLookUpTableValue;
		    }
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		}
	    }
	}
      while (CurrentLookUpTableValue > 0)
	{
	  TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
	  --CurrentLookUpTableValue;
	}
      TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;
    }
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// return value = Hilbert space dimension

long BosonOnTorusShort::EvaluateHilbertSpaceDimension(int nbrBosons, int maxMomentum)
{
  return this->EvaluateHilbertSpaceDimension(nbrBosons, maxMomentum - 1, maxMomentum - 1, 0);
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson in the state
// currentKyMax = momentum maximum value for bosons that are still to be placed
// currentMomentum = current value of the momentum
// return value = Hilbert space dimension

long BosonOnTorusShort::EvaluateHilbertSpaceDimension(int nbrBosons, int maxMomentum, int currentKyMax, int currentMomentum)
{
  if (nbrBosons == 0)
    {
      if ((currentMomentum % this->KyMomentumModulo) == this->TotalKy)
	{
	  return 1l;
	}
      else
	{
	  return 0l;
	}
    }
  if (currentKyMax == 0)
    {
      if ((currentMomentum % this->KyMomentumModulo) == this->TotalKy)
	{
	  return 1l;
	}
      else
	return 0l;
    }

  int TmpNbrBosons = 0;
  int ReducedCurrentKyMax = currentKyMax - 1;
  long TmpNbrStates = 0l;
  while (TmpNbrBosons < nbrBosons)
    {
     TmpNbrStates += this->EvaluateHilbertSpaceDimension(TmpNbrBosons, maxMomentum, ReducedCurrentKyMax, currentMomentum + (nbrBosons - TmpNbrBosons) * currentKyMax);
     ++TmpNbrBosons;
    }
  if (maxMomentum == currentKyMax)
    return TmpNbrStates + this->EvaluateHilbertSpaceDimension(nbrBosons, ReducedCurrentKyMax, ReducedCurrentKyMax, currentMomentum);
  else
    return TmpNbrStates + this->EvaluateHilbertSpaceDimension(nbrBosons, maxMomentum, ReducedCurrentKyMax, currentMomentum);
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson in the state
// currentKyMax = momentum maximum value for bosons that are still to be placed
// currentMomentum = current value of the momentum
// return value = Hilbert space dimension

long BosonOnTorusShort::EvaluateHilbertSpaceDimensionTruncatedOccupation(int nbrBosons, int maxMomentum, int currentKyMax, int currentMomentum)
{
  if (nbrBosons == 0)
    {
      if ((currentMomentum % this->KyMomentumModulo) == this->TotalKy)
	{
	  return 1l;
	}
      else
	{
	  return 0l;
	}
    }
  if (currentKyMax == 0)
    {
      if (((currentMomentum % this->KyMomentumModulo) == this->TotalKy) && (nbrBosons <= this->MaxOccupation))
	{
	  return 1l;
	}
      else
	return 0l;
    }

  int TmpNbrBosons = nbrBosons -  this->MaxOccupation;
  if (TmpNbrBosons < 0)
    TmpNbrBosons = 0;
  int ReducedCurrentKyMax = currentKyMax - 1;
  long TmpNbrStates = 0l;
  while (TmpNbrBosons < nbrBosons)
    {
     TmpNbrStates += this->EvaluateHilbertSpaceDimensionTruncatedOccupation(TmpNbrBosons, maxMomentum, ReducedCurrentKyMax, currentMomentum + (nbrBosons - TmpNbrBosons) * currentKyMax);
     ++TmpNbrBosons;
    }
  if (maxMomentum == currentKyMax)
    return TmpNbrStates + this->EvaluateHilbertSpaceDimensionTruncatedOccupation(nbrBosons, ReducedCurrentKyMax, ReducedCurrentKyMax, currentMomentum);
  else
    return TmpNbrStates + this->EvaluateHilbertSpaceDimensionTruncatedOccupation(nbrBosons, maxMomentum, ReducedCurrentKyMax, currentMomentum);
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Ky sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// kySector = Ky sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  BosonOnTorusShort::EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, RealVector& groundState)
{  
  if (nbrBosonSector == 0)
    {
      if (kySector == 0)
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

  if (nbrBosonSector == this->NbrBosons)
    {
      if (kySector == this->TotalKy)
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

  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;
  int ComplementaryKySector = this->TotalKy - kySector;
  if (ComplementaryKySector < 0)
    ComplementaryKySector += this->KyMax;
  if (ComplementaryKySector >= this->KyMax)
    ComplementaryKySector -= this->KyMax;

  if (nbrBosonSector == 1)
    {
      double TmpValue = 0.0;
      unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
      unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
      cout << "ComplementaryKySector = " << ComplementaryKySector << endl;
      BosonOnTorusShort TmpHilbertSpace(this->NbrBosons - 1, this->KyMax, ComplementaryKySector);
      FactorialCoefficient Factorial;
//       for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
// 	TmpHilbertSpace.PrintState(cout, MinIndex) << " | " << TmpHilbertSpace.StateKyMax[MinIndex] << " | " << hex << TmpHilbertSpace.StateDescription[MinIndex] << dec << endl;
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.StateKyMax[MinIndex] + TmpHilbertSpace.NbrBosons - 1, TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateKyMax);
	  TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.TemporaryStateKyMax + TmpHilbertSpace.NbrBosons - 1, TmpMonomial1);
	  int TmpIndex2 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial1[TmpIndex2] >= kySector))
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  TmpMonomial3[TmpIndex4] = kySector;
	  ++TmpIndex4;
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FindStateIndex(TmpState, TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
	      Factorial.SetToOne();
	      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateKyMax; ++k)
		if (TmpHilbertSpace.TemporaryState[k] > 1)
		  Factorial.FactorialDivide(TmpHilbertSpace.TemporaryState[k]);
	      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
		if (this->TemporaryState[k] > 1)
		  Factorial.FactorialMultiply(this->TemporaryState[k]);
	      Factorial.BinomialDivide(this->NbrBosons, nbrBosonSector);
	      TmpValue += groundState[TmpPos] * groundState[TmpPos] * (Factorial.GetNumericalValue());	
	    }
	}
      
      RealSymmetricMatrix TmpDensityMatrix(1);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
      return TmpDensityMatrix;
    }


  BosonOnTorusShort TmpDestinationHilbertSpace(nbrBosonSector, this->KyMax, kySector);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;
  unsigned long* TmpMonomial2 = new unsigned long [nbrBosonSector];
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
  BosonOnTorusShort TmpHilbertSpace(ComplementaryNbrBosonSector, this->KyMax, ComplementaryKySector);
  FactorialCoefficient Factorial;

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.StateKyMax[MinIndex] + TmpHilbertSpace.NbrBosons - 1, TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateKyMax);
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.TemporaryStateKyMax + TmpHilbertSpace.NbrBosons - 1 , TmpMonomial1);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  TmpDestinationHilbertSpace.FermionToBoson(TmpDestinationHilbertSpace.StateDescription[j], TmpDestinationHilbertSpace.StateKyMax[j] + TmpDestinationHilbertSpace.NbrBosons - 1, TmpDestinationHilbertSpace.TemporaryState, TmpDestinationHilbertSpace.TemporaryStateKyMax);
	  TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.StateDescription[j], TmpDestinationHilbertSpace.TemporaryStateKyMax + TmpDestinationHilbertSpace.NbrBosons - 1, TmpMonomial2);

	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < nbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonSector)
		{
		  while ((TmpIndex3 < nbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
		    {
		      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < nbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FindStateIndex(TmpState, TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
 	      Factorial.SetToOne();
 	      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateKyMax; ++k)
		{
		  if (TmpHilbertSpace.TemporaryState[k] > 1)
		    Factorial.FactorialDivide(TmpHilbertSpace.TemporaryState[k]);
		}
  	      for (int k = 0; k <= TmpDestinationHilbertSpace.TemporaryStateKyMax; ++k)
		{
		  if (TmpDestinationHilbertSpace.TemporaryState[k] > 1)
		    Factorial.FactorialDivide(TmpDestinationHilbertSpace.TemporaryState[k]);
		}
 	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
 	      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
		if (this->TemporaryState[k] > 1)
 		  Factorial.FactorialMultiply(this->TemporaryState[k]);
	      Factorial.BinomialDivide(this->NbrBosons, nbrBosonSector);
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      TmpStateCoefficient[Pos] = sqrt(Factorial.GetNumericalValue());
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  delete[] TmpStateCoefficient;
  delete[] TmpMonomial1;
  delete[] TmpMonomial2;
  delete[] TmpMonomial3;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Ky sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// kySector = Ky sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix  BosonOnTorusShort::EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, ComplexVector& groundState)
{  
  if (nbrBosonSector == 0)
    {
      if (kySector == 0)
	{
	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  if (nbrBosonSector == this->NbrBosons)
    {
      if (kySector == this->TotalKy)
	{
	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;
  int ComplementaryKySector = this->TotalKy - kySector;
  if (ComplementaryKySector < 0)
    ComplementaryKySector += this->KyMax;
  if (ComplementaryKySector >= this->KyMax)
    ComplementaryKySector -= this->KyMax;

  if (nbrBosonSector == 1)
    {
      Complex TmpValue = 0.0;
      unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
      unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
      cout << "ComplementaryKySector = " << ComplementaryKySector << endl;
      BosonOnTorusShort TmpHilbertSpace(this->NbrBosons - 1, this->KyMax, ComplementaryKySector);
      FactorialCoefficient Factorial;
//       for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
// 	TmpHilbertSpace.PrintState(cout, MinIndex) << " | " << TmpHilbertSpace.StateKyMax[MinIndex] << " | " << hex << TmpHilbertSpace.StateDescription[MinIndex] << dec << endl;
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.StateKyMax[MinIndex] + TmpHilbertSpace.NbrBosons - 1, TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateKyMax);
	  TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.TemporaryStateKyMax + TmpHilbertSpace.NbrBosons - 1, TmpMonomial1);
	  int TmpIndex2 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial1[TmpIndex2] >= kySector))
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  TmpMonomial3[TmpIndex4] = kySector;
	  ++TmpIndex4;
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FindStateIndex(TmpState, TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
	      Factorial.SetToOne();
	      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateKyMax; ++k)
		if (TmpHilbertSpace.TemporaryState[k] > 1)
		  Factorial.FactorialDivide(TmpHilbertSpace.TemporaryState[k]);
	      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
		if (this->TemporaryState[k] > 1)
		  Factorial.FactorialMultiply(this->TemporaryState[k]);
	      Factorial.BinomialDivide(this->NbrBosons, nbrBosonSector);
	      TmpValue += Conj(groundState[TmpPos]) * groundState[TmpPos] * (Factorial.GetNumericalValue());	
	    }
	}
      
      HermitianMatrix TmpDensityMatrix(1);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
      return TmpDensityMatrix;
    }


  BosonOnTorusShort TmpDestinationHilbertSpace(nbrBosonSector, this->KyMax, kySector);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;
  unsigned long* TmpMonomial2 = new unsigned long [nbrBosonSector];
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
  BosonOnTorusShort TmpHilbertSpace(ComplementaryNbrBosonSector, this->KyMax, ComplementaryKySector);
  FactorialCoefficient Factorial;

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.StateKyMax[MinIndex] + TmpHilbertSpace.NbrBosons - 1, TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateKyMax);
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.TemporaryStateKyMax + TmpHilbertSpace.NbrBosons - 1 , TmpMonomial1);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  TmpDestinationHilbertSpace.FermionToBoson(TmpDestinationHilbertSpace.StateDescription[j], TmpDestinationHilbertSpace.StateKyMax[j] + TmpDestinationHilbertSpace.NbrBosons - 1, TmpDestinationHilbertSpace.TemporaryState, TmpDestinationHilbertSpace.TemporaryStateKyMax);
	  TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.StateDescription[j], TmpDestinationHilbertSpace.TemporaryStateKyMax + TmpDestinationHilbertSpace.NbrBosons - 1, TmpMonomial2);

	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < nbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonSector)
		{
		  while ((TmpIndex3 < nbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
		    {
		      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < nbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FindStateIndex(TmpState, TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
 	      Factorial.SetToOne();
 	      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateKyMax; ++k)
		{
		  if (TmpHilbertSpace.TemporaryState[k] > 1)
		    Factorial.FactorialDivide(TmpHilbertSpace.TemporaryState[k]);
		}
  	      for (int k = 0; k <= TmpDestinationHilbertSpace.TemporaryStateKyMax; ++k)
		{
		  if (TmpDestinationHilbertSpace.TemporaryState[k] > 1)
		    Factorial.FactorialDivide(TmpDestinationHilbertSpace.TemporaryState[k]);
		}
 	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
 	      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
		if (this->TemporaryState[k] > 1)
 		  Factorial.FactorialMultiply(this->TemporaryState[k]);
	      Factorial.BinomialDivide(this->NbrBosons, nbrBosonSector);
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      TmpStateCoefficient[Pos] = sqrt(Factorial.GetNumericalValue());
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      Complex TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * Conj(groundState[TmpStatePosition[k]]) * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  delete[] TmpStateCoefficient;
  delete[] TmpMonomial1;
  delete[] TmpMonomial2;
  delete[] TmpMonomial3;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnTorusShort::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int kySector, RealVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrBosonSector == 0))
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
  if (subsytemSize > this->KyMax)
    {
      if (((kySector % this->KyMomentumModulo) == this->TotalKy) && (nbrBosonSector == this->NbrBosons))
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

  if (nbrBosonSector == 0)
    {
      if (kySector == 0)
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  double Coefficient = 0.0;
	  unsigned long Mask  = (0x1ul << (subsytemSize + nbrBosonSector)) - 1ul;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
                Coefficient += groundState[i] * groundState[i];
            }
	  TmpDensityMatrix.SetMatrixElement(0, 0, Coefficient);
//	  cout << "check " << TmpDensityMatrix(0, 0) << endl;
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  if (nbrBosonSector == this->NbrBosons)
    {
      if ((kySector % this->KyMomentumModulo) == this->TotalKy)
	{
	  BosonOnDiskShort TmpDestinationHilbertSpace(nbrBosonSector, kySector, subsytemSize - 1);
	  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
	  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	    {
	      int TmpPos = this->FindStateIndex(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i]);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpDensityMatrix.AddToMatrixElement(i, i, groundState[TmpPos] * groundState[TmpPos]);
		  for (int j = i + 1; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
		    {
		      int TmpPos2 = this->FindStateIndex(TmpDestinationHilbertSpace.FermionBasis->StateDescription[j], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[j]);
		      if (TmpPos2 != this->HilbertSpaceDimension)
			{
			  TmpDensityMatrix.AddToMatrixElement(i, j, groundState[TmpPos] * groundState[TmpPos2]);
			}
		    }
		}
 	    }
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  long TmpNbrNonZeroElements = 0;
  BosonOnDiskShort TmpDestinationHilbertSpace(nbrBosonSector, kySector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrBosons - nbrBosonSector) * subsytemSize);
  while (ComplementaryTotalKy < 0)
    ComplementaryTotalKy += this->KyMax;
  ComplementaryTotalKy = ComplementaryTotalKy % this->KyMomentumModulo;
  int ComplementaryMaxTotalKy = (this->KyMax - subsytemSize - 1) * (this->NbrBosons - nbrBosonSector);    

  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      BosonOnDiskShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  int Pos = 0;
	  unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax + this->NbrBosons - 1;
	      while (((TmpState >> TmpKyMax) & 0x1ul) == 0x0ul)
		--TmpKyMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpKyMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  ++Pos;
		}
	    }
	  if (Pos != 0)
	    {
	      ++TmpNbrNonZeroElements;
	      for (int j = 0; j < Pos; ++j)
		{
		  int Pos2 = TmpStatePosition2[j];
		  double TmpValue = groundState[TmpStatePosition[j]];
		  for (int k = 0; k < Pos; ++k)
		    if (TmpStatePosition2[k] >= Pos2)
		      TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
		}
	    }
	}
      ComplementaryTotalKy += this->KyMomentumModulo;
    }
  
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnTorusShort::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int kySector, ComplexVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrBosonSector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  if (subsytemSize > this->KyMax)
    {
      if (((kySector % this->KyMomentumModulo) == this->TotalKy) && (nbrBosonSector == this->NbrBosons))
	{
	  HermitianMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * Conj(groundState[j]));
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  if (nbrBosonSector == 0)
    {
      if (kySector == 0)
	{
	  HermitianMatrix TmpDensityMatrix(1);
	  Complex Coefficient = 0.0;
	  unsigned long Mask  = (0x1ul << (subsytemSize + nbrBosonSector)) - 1ul;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
                Coefficient += groundState[i] * Conj(groundState[i]);
            }
	  TmpDensityMatrix.SetMatrixElement(0, 0, Coefficient);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  if (nbrBosonSector == this->NbrBosons)
    {
      if ((kySector % this->KyMomentumModulo) == this->TotalKy)
	{
	  BosonOnDiskShort TmpDestinationHilbertSpace(nbrBosonSector, kySector, subsytemSize - 1);
	  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
	  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	    {
	      int TmpPos = this->FindStateIndex(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i]);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpDensityMatrix.AddToMatrixElement(i, i, groundState[TmpPos] * Conj(groundState[TmpPos]));
		  for (int j = i + 1; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
		    {
		      int TmpPos2 = this->FindStateIndex(TmpDestinationHilbertSpace.FermionBasis->StateDescription[j], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[j]);
		      if (TmpPos2 != this->HilbertSpaceDimension)
			{
			  TmpDensityMatrix.AddToMatrixElement(i, j, groundState[TmpPos] * Conj(groundState[TmpPos2]));
			}
		    }
		}
 	    }
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  long TmpNbrNonZeroElements = 0;
  BosonOnDiskShort TmpDestinationHilbertSpace(nbrBosonSector, kySector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrBosons - nbrBosonSector) * subsytemSize);
  while (ComplementaryTotalKy < 0)
    ComplementaryTotalKy += this->KyMax;
  ComplementaryTotalKy = ComplementaryTotalKy % this->KyMomentumModulo;
  int ComplementaryMaxTotalKy = (this->KyMax - subsytemSize - 1) * (this->NbrBosons - nbrBosonSector);    

  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      BosonOnDiskShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  int Pos = 0;
	  unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax + this->NbrBosons - 1;
	      while (((TmpState >> TmpKyMax) & 0x1ul) == 0x0ul)
		--TmpKyMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpKyMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  ++Pos;
		}
	    }
	  if (Pos != 0)
	    {
	      ++TmpNbrNonZeroElements;
	      for (int j = 0; j < Pos; ++j)
		{
		  int Pos2 = TmpStatePosition2[j];
		  Complex TmpValue = groundState[TmpStatePosition[j]];
		  for (int k = 0; k < Pos; ++k)
		    if (TmpStatePosition2[k] >= Pos2)
		      TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * Conj(groundState[TmpStatePosition[k]]));
		}
	    }
	}
      ComplementaryTotalKy += this->KyMomentumModulo;
    }
  
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// apply a magnetic translation along x to a given state
//
// index = state index 
// return value = translated state index

int BosonOnTorusShort::ApplyXMagneticTranslation(int index)
{
  unsigned long TmpState = this->StateDescription[index];
  for (int i = 0; i < this->StateShift;)
    {
      while ((i < this->StateShift) && ((TmpState & 0x1ul) == 0x0ul))
	{
	  TmpState >>= 1;
	  ++i;
	}
      if (i < this->StateShift)
	{
	  while ((TmpState & 0x1ul) == 0x1ul)
	    {
	      TmpState >>= 1;
	      TmpState |= this->LastMomentumMask;
	    }
	  TmpState >>= 1;	  
	  ++i;
	}
    }
  int TmpKyMax = this->NbrBosons + this->KyMax - 1;
  while ((TmpState >> TmpKyMax) == 0x0ul)
    --TmpKyMax;
  return this->FindStateIndex(TmpState, TmpKyMax);
}




// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem

RealMatrix BosonOnTorusShort::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrBosonSector, int kySector, RealVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrBosonSector == 0))
	{
	  RealMatrix TmpEntanglementMatrix(1,1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  
  if (subsytemSize == this->KyMax)
    {
      if ((kySector == this->TotalKy) && (nbrBosonSector == this->NbrBosons))
	{
	  RealMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[i]);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  long TmpNbrNonZeroElements = 0;
  if (nbrBosonSector == 0)
    {
      if (kySector == 0)
	{
	  double Coefficient = 0.0;
	  unsigned long Mask  = (0x1ul << (subsytemSize + nbrBosonSector)) - 1ul;
	  int TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
                ++TmpSize;
            }
	  if (TmpSize == 0)
	    {
	      RealMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
	    }
	  RealMatrix TmpEntanglementMatrix(1, TmpSize, true);
	  TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
		{
 		  TmpEntanglementMatrix.AddToMatrixElement(0, TmpSize, groundState[i]);
		  ++TmpSize;
		}
	      Coefficient += groundState[i] * groundState[i];
            }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  int MinIndex = 0;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrBosons - nbrBosonSector) * subsytemSize);
  while (ComplementaryTotalKy < 0)
    ComplementaryTotalKy += this->KyMomentumModulo;
  int ComplementaryMaxTotalKy = (this->KyMax - subsytemSize - 1) * (this->NbrBosons - nbrBosonSector);    
  int NbrComplementarySpaces = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMomentumModulo;
    }
  if (NbrComplementarySpaces == 0)
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;	  
    }
  BosonOnDiskShort** ComplementaryHilbertSpaces = new BosonOnDiskShort*[NbrComplementarySpaces];
  NbrComplementarySpaces = 0;
  ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrBosons - nbrBosonSector) * subsytemSize);
  while (ComplementaryTotalKy < 0)
    ComplementaryTotalKy += this->KyMomentumModulo;
  int TotalComplementaryHilbertSpaceDimension = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ComplementaryHilbertSpaces[NbrComplementarySpaces] = new BosonOnDiskShort(this->NbrBosons - nbrBosonSector, 
										ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[NbrComplementarySpaces]->GetHilbertSpaceDimension();
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMomentumModulo;
    }
 
  BosonOnDiskShort TmpDestinationHilbertSpace(nbrBosonSector, kySector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,
				   TotalComplementaryHilbertSpaceDimension, true);
  TotalComplementaryHilbertSpaceDimension = 0;
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      for (int MinIndex = 0; MinIndex < ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpComplementaryState = ComplementaryHilbertSpaces[SpaceIndex]->FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax + this->NbrBosons - 1;
	      while (((TmpState >> TmpKyMax) & 0x1ul) == 0x0ul)
		--TmpKyMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpKyMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  ++TmpNbrNonZeroElements;
		  TmpEntanglementMatrix.AddToMatrixElement(j, TotalComplementaryHilbertSpaceDimension + MinIndex, 
							   groundState[TmpPos]);
		}
	    }
	}
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension;
    }
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      delete ComplementaryHilbertSpaces[SpaceIndex];
    }
  delete[] ComplementaryHilbertSpaces;
  if (TmpNbrNonZeroElements > 0)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
}
  
// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem

ComplexMatrix BosonOnTorusShort::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrBosonSector, int kySector, ComplexVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrBosonSector == 0))
	{
	  ComplexMatrix TmpEntanglementMatrix(1,1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  
  if (subsytemSize == this->KyMax)
    {
      if ((kySector == this->TotalKy) && (nbrBosonSector == this->NbrBosons))
	{
	  ComplexMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[i]);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  long TmpNbrNonZeroElements = 0;
  if (nbrBosonSector == 0)
    {
      if (kySector == 0)
	{
	  Complex Coefficient = 0.0;
	  unsigned long Mask  = (0x1ul << (subsytemSize + nbrBosonSector)) - 1ul;
	  int TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
                ++TmpSize;
            }
	  if (TmpSize == 0)
	    {
	      ComplexMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
	    }
	  ComplexMatrix TmpEntanglementMatrix(1, TmpSize, true);
	  TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
		{
 		  TmpEntanglementMatrix.AddToMatrixElement(0, TmpSize, groundState[i]);
		  ++TmpSize;
		}
	      Coefficient += groundState[i] * groundState[i];
            }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  int MinIndex = 0;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrBosons - nbrBosonSector) * subsytemSize);
  while (ComplementaryTotalKy < 0)
    ComplementaryTotalKy += this->KyMomentumModulo;
  int ComplementaryMaxTotalKy = (this->KyMax - subsytemSize - 1) * (this->NbrBosons - nbrBosonSector);    
  int NbrComplementarySpaces = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMomentumModulo;
    }
  if (NbrComplementarySpaces == 0)
    {
      ComplexMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;	  
    }
  BosonOnDiskShort** ComplementaryHilbertSpaces = new BosonOnDiskShort*[NbrComplementarySpaces];
  NbrComplementarySpaces = 0;
  ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrBosons - nbrBosonSector) * subsytemSize);
  while (ComplementaryTotalKy < 0)
    ComplementaryTotalKy += this->KyMomentumModulo;
  int TotalComplementaryHilbertSpaceDimension = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ComplementaryHilbertSpaces[NbrComplementarySpaces] = new BosonOnDiskShort(this->NbrBosons - nbrBosonSector, 
										ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[NbrComplementarySpaces]->GetHilbertSpaceDimension();
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMomentumModulo;
    }
 
  BosonOnDiskShort TmpDestinationHilbertSpace(nbrBosonSector, kySector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  ComplexMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,
				   TotalComplementaryHilbertSpaceDimension, true);
  TotalComplementaryHilbertSpaceDimension = 0;
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      for (int MinIndex = 0; MinIndex < ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpComplementaryState = ComplementaryHilbertSpaces[SpaceIndex]->FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax + this->NbrBosons - 1;
	      while (((TmpState >> TmpKyMax) & 0x1ul) == 0x0ul)
		--TmpKyMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpKyMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  ++TmpNbrNonZeroElements;
		  TmpEntanglementMatrix.AddToMatrixElement(j, TotalComplementaryHilbertSpaceDimension + MinIndex, 
							   groundState[TmpPos]);
		}
	    }
	}
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension;
    }
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      delete ComplementaryHilbertSpaces[SpaceIndex];
    }
  delete[] ComplementaryHilbertSpaces;
  if (TmpNbrNonZeroElements > 0)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      ComplexMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
}
  
// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. 
// The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles for the part A 
// but without the Ky constraint for the part B
// 
// subsytemSize = number of states that belong to the subsytem
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem

RealMatrix BosonOnTorusShort::EvaluatePartialEntanglementMatrixFullKyPartB (int subsytemSize, int nbrBosonSector, int kySector, RealVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrBosonSector == 0))
	{
	  RealMatrix TmpEntanglementMatrix(1,1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  
  if (subsytemSize == this->KyMax)
    {
      if ((kySector == this->TotalKy) && (nbrBosonSector == this->NbrBosons))
	{
	  RealMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[i]);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  long TmpNbrNonZeroElements = 0;
  if (nbrBosonSector == 0)
    {
      if (kySector == 0)
	{
	  double Coefficient = 0.0;
	  unsigned long Mask  = (0x1ul << (subsytemSize + nbrBosonSector)) - 1ul;
	  int TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
                ++TmpSize;
            }
	  if (TmpSize == 0)
	    {
	      RealMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
	    }
	  BosonOnSphereFullShort TmpFullComplementaryHilbertSpace(this->NbrBosons - nbrBosonSector, this->KyMax - subsytemSize - 1);
	  RealMatrix TmpEntanglementMatrix(1, TmpFullComplementaryHilbertSpace.HilbertSpaceDimension, true);
	  TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
		{
		  unsigned long TmpState = this->StateDescription[i] >> (subsytemSize + nbrBosonSector);
		  int TmpKyMax = this->KyMax - subsytemSize - 1;
		  while ((TmpState >> TmpKyMax) == 0x0ul)
		    --TmpKyMax;
 		  TmpEntanglementMatrix.AddToMatrixElement(0, TmpFullComplementaryHilbertSpace.FindStateIndex(TmpState, TmpKyMax), groundState[i]);
		  ++TmpSize;
		}
	      Coefficient += groundState[i] * groundState[i];
            }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  int MinIndex = 0;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrBosons - nbrBosonSector) * subsytemSize);
  while (ComplementaryTotalKy < 0)
    ComplementaryTotalKy += this->KyMomentumModulo;
  int ComplementaryMaxTotalKy = (this->KyMax - subsytemSize - 1) * (this->NbrBosons - nbrBosonSector);    
  int NbrComplementarySpaces = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMomentumModulo;
    }
  if (NbrComplementarySpaces == 0)
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;	  
    }
  BosonOnDiskShort** ComplementaryHilbertSpaces = new BosonOnDiskShort*[NbrComplementarySpaces];
  NbrComplementarySpaces = 0;
  ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrBosons - nbrBosonSector) * subsytemSize);
  while (ComplementaryTotalKy < 0)
    ComplementaryTotalKy += this->KyMomentumModulo;
  int TotalComplementaryHilbertSpaceDimension = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ComplementaryHilbertSpaces[NbrComplementarySpaces] = new BosonOnDiskShort(this->NbrBosons - nbrBosonSector, 
										ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[NbrComplementarySpaces]->GetHilbertSpaceDimension();
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMomentumModulo;
    }
 
  BosonOnDiskShort TmpDestinationHilbertSpace(nbrBosonSector, kySector, subsytemSize - 1);
  BosonOnSphereFullShort TmpFullComplementaryHilbertSpace(this->NbrBosons - nbrBosonSector, this->KyMax - subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,
				   TmpFullComplementaryHilbertSpace.HilbertSpaceDimension, true);
  TotalComplementaryHilbertSpaceDimension = 0;
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      for (int MinIndex = 0; MinIndex < ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpComplementaryState = ComplementaryHilbertSpaces[SpaceIndex]->FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	  int ComplementaryIndex = TmpFullComplementaryHilbertSpace.FindStateIndex(ComplementaryHilbertSpaces[SpaceIndex]->FermionBasis->StateDescription[MinIndex], 
										   ComplementaryHilbertSpaces[SpaceIndex]->FermionBasis->StateLzMax[MinIndex]);
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax + this->NbrBosons - 1;
	      while (((TmpState >> TmpKyMax) & 0x1ul) == 0x0ul)
		--TmpKyMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpKyMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  ++TmpNbrNonZeroElements;
		  TmpEntanglementMatrix.AddToMatrixElement(j, ComplementaryIndex, 
							   groundState[TmpPos]);
		}
	    }
	}
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension;
    }
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      delete ComplementaryHilbertSpaces[SpaceIndex];
    }
  delete[] ComplementaryHilbertSpaces;
  if (TmpNbrNonZeroElements > 0)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
}
  
// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. 
// The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles for the part A 
// but without the Ky constraint for the part B
// 
// subsytemSize = number of states that belong to the subsytem
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem

ComplexMatrix BosonOnTorusShort::EvaluatePartialEntanglementMatrixFullKyPartB (int subsytemSize, int nbrBosonSector, int kySector, ComplexVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrBosonSector == 0))
	{
	  ComplexMatrix TmpEntanglementMatrix(1,1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  
  if (subsytemSize == this->KyMax)
    {
      if ((kySector == this->TotalKy) && (nbrBosonSector == this->NbrBosons))
	{
	  ComplexMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[i]);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  long TmpNbrNonZeroElements = 0;
  if (nbrBosonSector == 0)
    {
      if (kySector == 0)
	{
	  Complex Coefficient = 0.0;
	  unsigned long Mask  = (0x1ul << (subsytemSize + nbrBosonSector)) - 1ul;
	  int TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
                ++TmpSize;
            }
	  if (TmpSize == 0)
	    {
	      ComplexMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
	    }
	  BosonOnSphereFullShort TmpFullComplementaryHilbertSpace(this->NbrBosons - nbrBosonSector, this->KyMax - subsytemSize - 1);
	  ComplexMatrix TmpEntanglementMatrix(1, TmpFullComplementaryHilbertSpace.HilbertSpaceDimension, true);
	  TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
		{
		  unsigned long TmpState = this->StateDescription[i] >> (subsytemSize + nbrBosonSector);
		  int TmpKyMax = this->KyMax - subsytemSize - 1;
		  while ((TmpState >> TmpKyMax) == 0x0ul)
		    --TmpKyMax;
 		  TmpEntanglementMatrix.AddToMatrixElement(0, TmpFullComplementaryHilbertSpace.FindStateIndex(TmpState, TmpKyMax), groundState[i]);
		  ++TmpSize;
		}
	      Coefficient += groundState[i] * groundState[i];
            }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  int MinIndex = 0;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrBosons - nbrBosonSector) * subsytemSize);
  while (ComplementaryTotalKy < 0)
    ComplementaryTotalKy += this->KyMomentumModulo;
  int ComplementaryMaxTotalKy = (this->KyMax - subsytemSize - 1) * (this->NbrBosons - nbrBosonSector);    
  int NbrComplementarySpaces = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMomentumModulo;
    }
  if (NbrComplementarySpaces == 0)
    {
      ComplexMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;	  
    }
  BosonOnDiskShort** ComplementaryHilbertSpaces = new BosonOnDiskShort*[NbrComplementarySpaces];
  NbrComplementarySpaces = 0;
  ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrBosons - nbrBosonSector) * subsytemSize);
  while (ComplementaryTotalKy < 0)
    ComplementaryTotalKy += this->KyMomentumModulo;
  int TotalComplementaryHilbertSpaceDimension = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ComplementaryHilbertSpaces[NbrComplementarySpaces] = new BosonOnDiskShort(this->NbrBosons - nbrBosonSector, 
										ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[NbrComplementarySpaces]->GetHilbertSpaceDimension();
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMomentumModulo;
    }
 
  BosonOnDiskShort TmpDestinationHilbertSpace(nbrBosonSector, kySector, subsytemSize - 1);
  BosonOnSphereFullShort TmpFullComplementaryHilbertSpace(this->NbrBosons - nbrBosonSector, this->KyMax - subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  ComplexMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,
				   TmpFullComplementaryHilbertSpace.HilbertSpaceDimension, true);
  TotalComplementaryHilbertSpaceDimension = 0;
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      for (int MinIndex = 0; MinIndex < ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpComplementaryState = ComplementaryHilbertSpaces[SpaceIndex]->FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	  int ComplementaryIndex = TmpFullComplementaryHilbertSpace.FindStateIndex(ComplementaryHilbertSpaces[SpaceIndex]->FermionBasis->StateDescription[MinIndex], 
										   ComplementaryHilbertSpaces[SpaceIndex]->FermionBasis->StateLzMax[MinIndex]);
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax + this->NbrBosons - 1;
	      while (((TmpState >> TmpKyMax) & 0x1ul) == 0x0ul)
		--TmpKyMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpKyMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  ++TmpNbrNonZeroElements;
		  TmpEntanglementMatrix.AddToMatrixElement(j, ComplementaryIndex, 
							   groundState[TmpPos]);
		}
	    }
	}
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension;
    }
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      delete ComplementaryHilbertSpaces[SpaceIndex];
    }
  delete[] ComplementaryHilbertSpaces;
  if (TmpNbrNonZeroElements > 0)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      ComplexMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
}
  
// remove part of each Fock state, discarding component if the Fock state does not a given pattern
//
// inputVector = state to truncate
// reducedSpace = Hilbert space where the truncated state will lie
// pattern = array describing the pattern 
// patternSize = pattern size
// patternShift = indicate where the pattern has to be applied
// return value = trucated state

RealVector BosonOnTorusShort::TruncateStateWithPatternConstraint(RealVector& inputVector, ParticleOnTorus* reducedSpace, int* pattern, int patternSize, int patternShift)
{
  BosonOnTorusShort* TmpSpace = (BosonOnTorusShort*) reducedSpace; 
  unsigned long PatternMask =  (0x1ul << patternSize) - 0x1ul;
  PatternMask <<= patternShift;
  unsigned long Pattern = 0x0ul;
  for (int i = 0; i < patternSize; ++i)
    {
      if (pattern[i] == 1)
	Pattern |= 0x1ul << i;
    }
  Pattern <<= patternShift;
  unsigned long Mask = (0x1ul << patternShift) - 0x1ul;
  RealVector TmpVector (TmpSpace->LargeHilbertSpaceDimension, true);
  unsigned long Tmp;
//   for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//     {      
//       Tmp = this->StateDescription[i]; 
//       if ((Tmp & PatternMask) == Pattern)
// 	{	  	  
// 	  Tmp = ((Tmp >> patternSize) & (~Mask)) | (Tmp & Mask);
// 	  int TmpKyMax = TmpSpace->LzMax;
// 	  while (((Tmp >> TmpLzMax) & 0x1ul) == 0x0ul)
// 	    --TmpLzMax;
// 	  int TmpIndex = TmpSpace->FindStateIndex(Tmp, TmpLzMax);
// 	  if (TmpIndex < TmpSpace->HilbertSpaceDimension)
// 	    TmpVector[TmpIndex] = inputVector[i];
// 	}
//     }
  return TmpVector;
 
}


// symmetrized a product of two uncoupled states 
//
// outputVector = reference on the vector which will contain the symmetrozed state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

void BosonOnTorusShort::SymmetrizeU1U1StateCore (RealVector& symmetrizedVector, RealVector& leftVector, RealVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents)
{
  if (unnormalizedBasisFlag == true)
    {
      cout << "Unnormalized basis not implemented" << endl;
    }
  else
    {
      long LastComponent = (long) (firstComponent + nbrComponents);
      BosonOnTorusShort* TmpLeftSpace = (BosonOnTorusShort*) leftSpace;
      BosonOnTorusShort* TmpRightSpace = (BosonOnTorusShort*) rightSpace;
     
      FactorialCoefficient Factorial1;
      FactorialCoefficient Factorial2;
      for (long i = (long) firstComponent; i < LastComponent; ++i)
	{
	  this->FermionToBoson(TmpLeftSpace->StateDescription[i], TmpLeftSpace->StateKyMax[i] + TmpLeftSpace->NbrBosons - 1, 
			       TmpLeftSpace->TemporaryState, TmpLeftSpace->TemporaryStateKyMax);
	  
	  for (int k = TmpLeftSpace->TemporaryStateKyMax + 1;  k < TmpLeftSpace->KyMax; ++k)
	    TmpLeftSpace->TemporaryState[k] = 0;
	  double TmpCoefficient = leftVector[i];
	  Factorial1.SetToOne();
	  Factorial1.Power2Divide(TmpLeftSpace->NbrBosons);
	  for (int k = 0; k <= TmpLeftSpace->TemporaryStateKyMax; ++k)
	    if (TmpLeftSpace->TemporaryState[k] > 1)
	      Factorial1.FactorialDivide(TmpLeftSpace->TemporaryState[k]);
	  
	  if (this->LargeHilbertSpaceDimension >= (1l << 30))
	    {
	      for (long j = 0l; j < TmpRightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  this->FermionToBoson(TmpRightSpace->StateDescription[j], TmpRightSpace->StateKyMax[j] + TmpRightSpace->NbrBosons - 1, 
				       TmpRightSpace->TemporaryState, TmpRightSpace->TemporaryStateKyMax);
		  int k = 0;
		  for (; k <= TmpRightSpace->TemporaryStateKyMax; ++k)
		    {
		      this->TemporaryState[k] = TmpLeftSpace->TemporaryState[k] + TmpRightSpace->TemporaryState[k];
		    }
		  this->TemporaryStateKyMax = TmpRightSpace->TemporaryStateKyMax;
		  if (TmpLeftSpace->TemporaryStateKyMax > TmpRightSpace->TemporaryStateKyMax)
		    {
		      for (; k <= TmpLeftSpace->TemporaryStateKyMax; ++k)
			{
			  this->TemporaryState[k] = TmpLeftSpace->TemporaryState[k];
			}
		      this->TemporaryStateKyMax = TmpLeftSpace->TemporaryStateKyMax;
		    }
		  
		  long TmpPos = this->FindStateLargeIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
		  if (TmpPos < this->LargeHilbertSpaceDimension)
		    {
		      Factorial2 = Factorial1;
		      for (k = 0; k <= TmpRightSpace->TemporaryStateKyMax; ++k)
			if (TmpRightSpace->TemporaryState[k] > 1)
			  Factorial2.FactorialDivide(TmpRightSpace->TemporaryState[k]);
		      for (k = 0; k <= this->TemporaryStateKyMax; ++k)
			if (this->TemporaryState[k] > 1)
			  Factorial2.FactorialMultiply(this->TemporaryState[k]);	
		      
		      symmetrizedVector[TmpPos] += sqrt(Factorial2.GetNumericalValue()) * TmpCoefficient * rightVector[j];
		    }		
		}
	    }
	  else
	    {
	      for (long j = 0l; j < TmpRightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  this->FermionToBoson(TmpRightSpace->StateDescription[j], TmpRightSpace->StateKyMax[j] + TmpRightSpace->NbrBosons - 1, 
				       TmpRightSpace->TemporaryState, TmpRightSpace->TemporaryStateKyMax);
		  int k = 0;
		  for (; k <= TmpRightSpace->TemporaryStateKyMax; ++k)
		    {
		      this->TemporaryState[k] = TmpLeftSpace->TemporaryState[k] + TmpRightSpace->TemporaryState[k];
		    }
		  this->TemporaryStateKyMax = TmpRightSpace->TemporaryStateKyMax;
		  if (TmpLeftSpace->TemporaryStateKyMax > TmpRightSpace->TemporaryStateKyMax)
		    {
		      for (; k <= TmpLeftSpace->TemporaryStateKyMax; ++k)
			{
			  this->TemporaryState[k] = TmpLeftSpace->TemporaryState[k];
			}
		      this->TemporaryStateKyMax = TmpLeftSpace->TemporaryStateKyMax;
		    }
		  
		  int TmpPos = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
		  if (TmpPos < this->HilbertSpaceDimension)
		    {
		      Factorial2 = Factorial1;
		      for (k = 0; k <= TmpRightSpace->TemporaryStateKyMax; ++k)
			if (TmpRightSpace->TemporaryState[k] > 1)
			  Factorial2.FactorialDivide(TmpRightSpace->TemporaryState[k]);
		      for (k = 0; k <= this->TemporaryStateKyMax; ++k)
			if (this->TemporaryState[k] > 1)
			  Factorial2.FactorialMultiply(this->TemporaryState[k]);	
		      
		      symmetrizedVector[TmpPos] += sqrt(Factorial2.GetNumericalValue()) * TmpCoefficient * rightVector[j];
		    }		
		}
	    }
	}  
    }
}

// symmetrized a product of two uncoupled states 
//
// outputVector = reference on the vector which will contain the symmetrozed state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

void BosonOnTorusShort::SymmetrizeU1U1StateCore (ComplexVector& symmetrizedVector, ComplexVector& leftVector, ComplexVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents)
{
  unsigned long LastComponent = firstComponent + nbrComponents;
  
  FactorialCoefficient Factorial1;
  FactorialCoefficient Factorial2;
  if (unnormalizedBasisFlag == true)
    {
      cout << "Unnormalized basis not implemented" << endl;
    }
  else
    {
      BosonOnTorusShort* TmpLeftSpace = (BosonOnTorusShort*) leftSpace;
      BosonOnTorusShort* TmpRightSpace = (BosonOnTorusShort*) rightSpace;
      for (long i = (long) firstComponent; i < (long) LastComponent; ++i)
	{
	  this->FermionToBoson(TmpLeftSpace->StateDescription[i], TmpLeftSpace->StateKyMax[i] + TmpLeftSpace->NbrBosons - 1, 
			       TmpLeftSpace->TemporaryState, TmpLeftSpace->TemporaryStateKyMax);
	  
	  for (int k = TmpLeftSpace->TemporaryStateKyMax + 1;  k < TmpLeftSpace->KyMax; ++k)
	    TmpLeftSpace->TemporaryState[k] = 0;
	  Complex TmpCoefficient = leftVector[i];
	  Factorial1.SetToOne();
	  Factorial1.Power2Divide(TmpLeftSpace->NbrBosons);
	  for (int k = 0; k <= TmpLeftSpace->TemporaryStateKyMax; ++k)
	    if (TmpLeftSpace->TemporaryState[k] > 1)
	      Factorial1.FactorialDivide(TmpLeftSpace->TemporaryState[k]);
	  
	  for (long j = 0l; j < TmpRightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      this->FermionToBoson(TmpRightSpace->StateDescription[j], TmpRightSpace->StateKyMax[j] + TmpRightSpace->NbrBosons - 1, 
				   TmpRightSpace->TemporaryState, TmpRightSpace->TemporaryStateKyMax);
	      int k = 0;
	      for (; k <= TmpRightSpace->TemporaryStateKyMax; ++k)
	      {
		this->TemporaryState[k] = TmpLeftSpace->TemporaryState[k] + TmpRightSpace->TemporaryState[k];
	      }
	      this->TemporaryStateKyMax = TmpRightSpace->TemporaryStateKyMax;
	      if (TmpLeftSpace->TemporaryStateKyMax > TmpRightSpace->TemporaryStateKyMax)
		{
		  for (; k <= TmpLeftSpace->TemporaryStateKyMax; ++k)
		  {
		    this->TemporaryState[k] = TmpLeftSpace->TemporaryState[k];
		  }
		  this->TemporaryStateKyMax = TmpLeftSpace->TemporaryStateKyMax;
		}

	      int TmpPos = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
	      if (TmpPos < this->HilbertSpaceDimension)
		{
		  Factorial2 = Factorial1;
		  for (k = 0; k <= TmpRightSpace->TemporaryStateKyMax; ++k)
		    if (TmpRightSpace->TemporaryState[k] > 1)
		      Factorial2.FactorialDivide(TmpRightSpace->TemporaryState[k]);
		  for (k = 0; k <= this->TemporaryStateKyMax; ++k)
		    if (this->TemporaryState[k] > 1)
		      Factorial2.FactorialMultiply(this->TemporaryState[k]);	
		    
		  symmetrizedVector[TmpPos] += sqrt(Factorial2.GetNumericalValue()) * TmpCoefficient * rightVector[j];
		}
		
	    }
	}  
    }
}

// symmetrize a vector with even number of orbitals 
//
// outputVector = reference on the vector which will contain the symmetrozed state
// leftVector = reference on the vector to be symmetrized
// leftSpace = pointer to the Hilbert space
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

RealVector BosonOnTorusShort::SymmetrizeU1U1SingleState (RealVector& leftVector,  ParticleOnTorus* leftSpace, bool oneInTwoFlag, bool unnormalizedBasisFlag, AbstractArchitecture* architecture)
{
  RealVector SymmetrizedVector (this->LargeHilbertSpaceDimension,true);

//   FQHETorusSymmetrizeU1U1StateOperation Operation (this, leftSpace, rightSpace, &SymmetrizedVector, &leftVector, &rightVector, unnormalizedBasisFlag);
//   Operation.ApplyOperation(architecture);
  unsigned long firstComponent = 0;
  unsigned long nbrComponent = leftSpace->GetHilbertSpaceDimension();
//   timeval TotalStartingTime;
//   gettimeofday (&TotalStartingTime, 0);
//   if (oneInTwoFlag == false)
//     this->SymmetrizeU1U1SingleStateCore (SymmetrizedVector, leftVector, leftSpace, firstComponent, nbrComponent);
//   else
//     this->SymmetrizeU1U1SingleStateOneInTwoCore (SymmetrizedVector, leftVector, leftSpace, firstComponent, nbrComponent);
  
  if ( unnormalizedBasisFlag == false )
  {
    if (SymmetrizedVector.Norm() > std::numeric_limits<double>::epsilon())
      SymmetrizedVector /= SymmetrizedVector.Norm();
  }
  cout << this->TotalKy << " " << SymmetrizedVector.Norm() << endl;
  return SymmetrizedVector;
}


// symmetrize a vector with even number of orbitals 
//
// outputVector = reference on the vector which will contain the symmetrozed state
// leftVector = reference on the vector to be symmetrized
// leftSpace = pointer to the Hilbert space
// groupNeighbouringOrbitals = group neighbouring orbitals instead of grouping orbitals separated by Nphi
// return value = symmetrized state

ComplexVector BosonOnTorusShort::SymmetrizeU1U1SingleState (ComplexVector& leftVector,  ParticleOnTorus* leftSpace, bool groupNeighbouringOrbitals, AbstractArchitecture* architecture)
{
  ComplexVector SymmetrizedVector (this->LargeHilbertSpaceDimension,true);
  Complex TmpCoefficient (0.0, -1.0);
  
  RealVector leftVectorReal(leftVector);
  RealVector leftVectorImaginary(leftSpace->GetLargeHilbertSpaceDimension(), true);
  leftVectorImaginary = TmpCoefficient * leftVector;
  
  RealVector TmpVector(this->LargeHilbertSpaceDimension, true);

//   FQHETorusSymmetrizeU1U1StateOperation Operation (this, leftSpace, rightSpace, &SymmetrizedVector, &leftVector, &rightVector, unnormalizedBasisFlag);
//   Operation.ApplyOperation(architecture);
  unsigned long firstComponent = 0ul;
  unsigned long nbrComponent = leftSpace->GetLargeHilbertSpaceDimension();
//   timeval TotalStartingTime;
//   gettimeofday (&TotalStartingTime, 0);
  if (groupNeighbouringOrbitals == false)
    {
      this->SymmetrizeU1U1SingleStateCore (TmpVector, leftVectorImaginary, leftSpace, firstComponent, nbrComponent);
      SymmetrizedVector += TmpVector;
      SymmetrizedVector *= TmpCoefficient;
      TmpVector.ClearVector(); 
      this->SymmetrizeU1U1SingleStateCore (TmpVector, leftVectorReal, leftSpace, firstComponent, nbrComponent);
      SymmetrizedVector -= TmpVector;
    }
  else
    {
      //      this->SymmetrizeU1U1SingleStateOneInTwoCore (TmpVector, leftVectorImaginary, leftSpace, firstComponent, nbrComponent);
      SymmetrizedVector += TmpVector;
      SymmetrizedVector *= TmpCoefficient;
      TmpVector.ClearVector(); 
      //      this->SymmetrizeU1U1SingleStateOneInTwoCore (TmpVector, leftVectorReal, leftSpace, firstComponent, nbrComponent);
      SymmetrizedVector -= TmpVector;
    }
  
  
  
//   timeval TotalEndingTime;
//   gettimeofday (&TotalEndingTime, 0);
//   double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
//   cout << this->FirstComponent << " " <<  this->NbrComponent << " : " << Dt << "s" << endl;
  if (SymmetrizedVector.Norm() != 0)
    SymmetrizedVector /= SymmetrizedVector.Norm();
  
  return SymmetrizedVector;
}

// symmetrize a vector with even number of orbitals (divides torus aspect ratio by two)
//
// symmetrizedVector = reference on the vector which will contain the symmetrozed state
// leftVector = reference on the vector to be symmetrized
// leftSpace = pointer to the Hilbert space
// first component = index of the first vector component 
// last component = index of the last component

void BosonOnTorusShort::SymmetrizeU1U1SingleStateCore (RealVector& symmetrizedVector, RealVector& leftVector,  ParticleOnTorus* leftSpace, unsigned long firstComponent, unsigned long nbrComponents)
{
  unsigned long LastComponent = firstComponent + nbrComponents;
  
  FactorialCoefficient Factorial1;
  FactorialCoefficient Factorial2;
  BosonOnTorusShort* TmpLeftSpace = (BosonOnTorusShort*) leftSpace;
  for (long i = (long) firstComponent; i < (long) LastComponent; ++i)
    {
      this->FermionToBoson(TmpLeftSpace->StateDescription[i], TmpLeftSpace->StateKyMax[i] + TmpLeftSpace->NbrBosons - 1, 
			   TmpLeftSpace->TemporaryState, TmpLeftSpace->TemporaryStateKyMax);
      
      for (int k = TmpLeftSpace->TemporaryStateKyMax + 1;  k < TmpLeftSpace->KyMax; ++k)
	TmpLeftSpace->TemporaryState[k] = 0;
      double TmpCoefficient = leftVector[i];
      Factorial1.SetToOne();
      Factorial1.Power2Divide(TmpLeftSpace->NbrBosons);
      
      for (int k = 0; k <= TmpLeftSpace->TemporaryStateKyMax; ++k)
	if (TmpLeftSpace->TemporaryState[k] > 1)
	  Factorial1.FactorialDivide(TmpLeftSpace->TemporaryState[k]);
      
      int k = 0;
      for (; k < this->KyMax; ++k)
	{
	  this->TemporaryState[k] = TmpLeftSpace->TemporaryState[k] + TmpLeftSpace->TemporaryState[k + this->KyMax];
	  if (TmpLeftSpace->TemporaryState[k] != 0 || TmpLeftSpace->TemporaryState[k + this->KyMax] != 0)
	    this->TemporaryStateKyMax = k;
	}
      
      int TmpPos = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
      if (TmpPos < this->HilbertSpaceDimension)
	{
	  for (k = 0; k <= this->TemporaryStateKyMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial1.FactorialMultiply(this->TemporaryState[k]);	
	  
	  symmetrizedVector[TmpPos] += sqrt(Factorial1.GetNumericalValue()) * TmpCoefficient;
	}
      
    }
}

// symmetrize a vector by grouping neighbouring orbitals, core part
//
// inputVector = reference on the vector to symmetrize
// nbrOrbitals = number of orbitals to group together
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
// first component = index of the first vector component 
// last component = index of the last component

void BosonOnTorusShort::SymmetrizeSingleStateGroupingNeighbouringOrbitalsCore (ComplexVector& inputVector, ComplexVector* symmetrizedVectors, int nbrOrbitals, unsigned long firstComponent, unsigned long nbrComponents)
{
  long LastComponent = (long) (firstComponent + nbrComponents);
  int TargetSpaceNbrOrbitals = this->KyMax / nbrOrbitals;
  BosonOnTorusShort** TargetSpaces = new BosonOnTorusShort* [TargetSpaceNbrOrbitals];
  unsigned long* TmpState = new unsigned long[TargetSpaceNbrOrbitals];
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      TargetSpaces[i] = 0;
    }  


  FactorialCoefficient Factorial1;
  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      Complex TmpCoefficient = inputVector[i];
      this->FermionToBoson(this->StateDescription[i], this->StateKyMax[i] + this->NbrBosons - 1, 
			   this->TemporaryState, this->TemporaryStateKyMax);
      Factorial1.SetToOne();
      Factorial1.Power2Divide(this->NbrBosons);
      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
	if (this->TemporaryState[k] > 1)
	  Factorial1.FactorialDivide(this->TemporaryState[k]);
      for (int k = this->TemporaryStateKyMax + 1; k < this->KyMax; ++k)
	this->TemporaryState[k] = 0x0ul;
      int TmpTotalKy = 0;
      for (int k = 0; k < TargetSpaceNbrOrbitals; ++k)
	{
	  unsigned long& TmpNbrParticles = TmpState[k];
	  TmpNbrParticles = 0x0ul;
	  for (int l = 0; l < nbrOrbitals; ++l)
	    TmpNbrParticles += this->TemporaryState[k * nbrOrbitals + l];
	  if (TmpNbrParticles > 1)
	    Factorial1.FactorialMultiply(TmpNbrParticles);	
	  TmpTotalKy += k * ((int) TmpNbrParticles);
	}
      TmpTotalKy %= TargetSpaceNbrOrbitals;

      if (TargetSpaces[TmpTotalKy] == 0)
	{
	  TargetSpaces[TmpTotalKy] = new BosonOnTorusShort (this->NbrBosons, TargetSpaceNbrOrbitals, TmpTotalKy);
	  symmetrizedVectors[TmpTotalKy] = ComplexVector(TargetSpaces[TmpTotalKy]->HilbertSpaceDimension, true);
	}	  
      int TmpKyMax = TargetSpaces[TmpTotalKy]->KyMax - 1;
      while (TmpState[TmpKyMax] == 0x0ul)
	--TmpKyMax;
      int TmpPos = TargetSpaces[TmpTotalKy]->FindStateIndex(TargetSpaces[TmpTotalKy]->BosonToFermion(TmpState, TmpKyMax), 
							    TmpKyMax + this->NbrBosons - 1);
      if (TmpPos < TargetSpaces[TmpTotalKy]->HilbertSpaceDimension)
	{
	  symmetrizedVectors[TmpTotalKy][TmpPos] += sqrt(Factorial1.GetNumericalValue()) * TmpCoefficient;
	}
    }
  delete[] TmpState;
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      if (TargetSpaces[i] != 0)
	{
	  delete TargetSpaces[i];
	}
    }
}


// symmetrize a vector by grouping distant and equally separated orbitals, core part
//
// inputVector = reference on the vector to symmetrize
// nbrOrbitals = number of orbitals to group together
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
// first component = index of the first vector component 
// last component = index of the last component

void BosonOnTorusShort::SymmetrizeSingleStateGroupingDistantOrbitalsCore (ComplexVector& inputVector, ComplexVector* symmetrizedVectors, int nbrOrbitals, unsigned long firstComponent, unsigned long nbrComponents, bool twistedTorus)
{
  long LastComponent = (long) (firstComponent + nbrComponents);
  int TargetSpaceNbrOrbitals = this->KyMax / nbrOrbitals;
  BosonOnTorusShort** TargetSpaces = new BosonOnTorusShort* [TargetSpaceNbrOrbitals];
  unsigned long* TmpState = new unsigned long[TargetSpaceNbrOrbitals];
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      TargetSpaces[i] = 0;
    }  

  FactorialCoefficient Factorial1;
  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      Complex TmpCoefficient = inputVector[i];
      this->FermionToBoson(this->StateDescription[i], this->StateKyMax[i] + this->NbrBosons - 1, 
			   this->TemporaryState, this->TemporaryStateKyMax);
      Factorial1.SetToOne();
      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
	if (this->TemporaryState[k] > 1)
	  Factorial1.FactorialDivide(this->TemporaryState[k]);
      for (int k = this->TemporaryStateKyMax + 1; k < this->KyMax; ++k)
	this->TemporaryState[k] = 0x0ul;
      int TmpTotalKy = 0;
      for (int k = 0; k < TargetSpaceNbrOrbitals; ++k)
	{
	  unsigned long& TmpNbrParticles = TmpState[k];
	  TmpNbrParticles = 0x0ul;
	  for (int l = 0; l < nbrOrbitals; ++l)
	  {
	    TmpNbrParticles += this->TemporaryState[k + (TargetSpaceNbrOrbitals * l)];
	    if (twistedTorus == true)
	      TmpCoefficient *= Phase (- (double) (this->TemporaryState[k + (TargetSpaceNbrOrbitals * l)] * (k + (TargetSpaceNbrOrbitals * l)) * (k + (TargetSpaceNbrOrbitals * l)) * nbrOrbitals) * M_PI/ ((double) (this->KyMax)));
	  }
	  if (TmpNbrParticles > 1)
	    Factorial1.FactorialMultiply(TmpNbrParticles);	
	  TmpTotalKy += k * ((int) TmpNbrParticles);
	}
      TmpTotalKy %= TargetSpaceNbrOrbitals;
      if (TargetSpaces[TmpTotalKy] == 0)
	{
	  TargetSpaces[TmpTotalKy] = new BosonOnTorusShort (this->NbrBosons, TargetSpaceNbrOrbitals, TmpTotalKy);
	  symmetrizedVectors[TmpTotalKy] = ComplexVector(TargetSpaces[TmpTotalKy]->HilbertSpaceDimension, true);
	}	  
      int TmpKyMax = TargetSpaces[TmpTotalKy]->KyMax - 1;
      while (TmpState[TmpKyMax] == 0x0ul)
	{
	  --TmpKyMax;
	}
      int TmpPos = TargetSpaces[TmpTotalKy]->FindStateIndex(TargetSpaces[TmpTotalKy]->BosonToFermion(TmpState, TmpKyMax), 
							    TmpKyMax + this->NbrBosons - 1);
      if (TmpPos < TargetSpaces[TmpTotalKy]->HilbertSpaceDimension)
	{
	  symmetrizedVectors[TmpTotalKy][TmpPos] += sqrt(Factorial1.GetNumericalValue()) * TmpCoefficient;
	}
    }
  delete[] TmpState;
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      if (TargetSpaces[i] != 0)
	{
	  delete TargetSpaces[i];
	}
    }
}

// symmetrize a vector by keeping only a subset of equally separated orbitals
//
// inputVector = reference on the vector to symmetrize
// firstOrbitalIndex = index of the first orbital to keep
// symmetrizedVectors = array on the symmetrize states ranging from the smallest Ky to the largest Ky
// periodicity = momentum periodicity (should be a multiple of the number of orbitals)
// firstComponent = first component of the input vector that has to be symmetrized
// nbrComponents = number of components of the input vector that have to be symmetrized
// return value = symmetrized state

void BosonOnTorusShort::SymmetrizeSingleStatePeriodicSubsetOrbitalCore (ComplexVector& inputVector, ComplexVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, 
									unsigned long firstComponent, unsigned long nbrComponents)
{
  bool twistedTorus = true;
  long LastComponent = (long) (firstComponent + nbrComponents);
  int TargetSpaceNbrOrbitals = this->KyMax / periodicity;
  BosonOnTorusShort*** TargetSpaces = new BosonOnTorusShort** [this->NbrBosons + 1];
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      TargetSpaces[i] = new BosonOnTorusShort* [TargetSpaceNbrOrbitals];
      for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
	{
	  TargetSpaces[i][j] = 0;
	}
    }
  unsigned long* TmpState = new unsigned long[TargetSpaceNbrOrbitals];
  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      Complex TmpCoefficient = inputVector[i];
      this->FermionToBoson(this->StateDescription[i], this->StateKyMax[i] + this->NbrBosons - 1, 
			   this->TemporaryState, this->TemporaryStateKyMax);
      for (int k = this->TemporaryStateKyMax + 1; k < this->KyMax; ++k)
	this->TemporaryState[k] = 0x0ul;
      int TmpTotalKy = 0;
      int TmpNbrParticles = 0;
      int Index = 0;
      double phase = 0.0;
      for (int k = firstOrbitalIndex; k < this->KyMax; k += periodicity)
	{
	  unsigned long& TmpNbrParticles2 = this->TemporaryState[k];
	  TmpState[Index] = TmpNbrParticles2;
	  TmpNbrParticles += ((int) TmpNbrParticles2);
	  TmpTotalKy += Index * ((int) TmpNbrParticles2);
// 	  if (twistedTorus == true)
// 	  {
// 	    phase += ((double) (this->TemporaryState[k] * k * k) * 2.0  / ((double) (this->KyMax)));
// 	    TmpCoefficient *= Phase ((double) (this->TemporaryState[k] * k * k) * 2.0 * M_PI / ((double) (this->KyMax)));
// // 	    cout << k << " " << ((double) (this->TemporaryState[k] * k * k) * 2.0  / ((double) (this->KyMax))) <<  " " << this->TemporaryState[k] << endl;
// 	  }
	  ++Index;
	}
      if (TmpNbrParticles > 0)
	{
// 	  if (TmpNbrParticles == this->NbrBosons)
// 	    cout << phase << endl;
	  TmpTotalKy %= TargetSpaceNbrOrbitals;
	  if (TargetSpaces[TmpNbrParticles][TmpTotalKy] == 0)
	    {
	      TargetSpaces[TmpNbrParticles][TmpTotalKy] = new BosonOnTorusShort (TmpNbrParticles, TargetSpaceNbrOrbitals, TmpTotalKy);
	      symmetrizedVectors[TmpNbrParticles][TmpTotalKy] = ComplexVector(TargetSpaces[TmpNbrParticles][TmpTotalKy]->HilbertSpaceDimension, true);
	    }	  
	  int TmpKyMax = TargetSpaces[TmpNbrParticles][TmpTotalKy]->KyMax - 1;
	  while (TmpState[TmpKyMax] == 0x0ul)
	    --TmpKyMax;
	  int TmpPos = TargetSpaces[TmpNbrParticles][TmpTotalKy]->FindStateIndex(TargetSpaces[TmpNbrParticles][TmpTotalKy]->BosonToFermion(TmpState, TmpKyMax), 
										 TmpKyMax + TmpNbrParticles - 1);
	  if (TmpPos < TargetSpaces[TmpNbrParticles][TmpTotalKy]->HilbertSpaceDimension)
	    {
	      symmetrizedVectors[TmpNbrParticles][TmpTotalKy][TmpPos] += TmpCoefficient;
	    }
	}
    }
  delete[] TmpState;
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
	{
	  if (TargetSpaces[i][j] != 0)
	    {
	      delete TargetSpaces[i][j];
	    }
	}
      delete[] TargetSpaces[i];
    }
  delete[] TargetSpaces;
}

// core part of the C4 rotation
//
// inputState = reference on the state that has to be rotated
// inputSpace = Hilbert space associated to the input state
// outputState = reference on the rotated state
// minIndex = minimum index that has to be computed
// nbrIndices = number of indices that have to be computed
// clockwise = the rotation is done clockwise
// return value = reference on the rotated state

ComplexVector& BosonOnTorusShort::CoreC4Rotation (ComplexVector& inputState, ParticleOnTorus* inputSpace, ComplexVector& outputState, int minIndex, int nbrIndices,
						  bool clockwise)
{
  BosonOnTorusShort* TmpInputSpace = (BosonOnTorusShort*) inputSpace;
  unsigned long* TmpInputMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpOutputMonomial = new unsigned long [this->NbrBosons];
  double* LogFactorialCoefficients = new double [this->NbrBosons + 1];
  LogFactorialCoefficients[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    LogFactorialCoefficients[i] = LogFactorialCoefficients[i - 1] + log((double) i);
  int LastIndex = minIndex + nbrIndices;
  Complex Tmp = 0.0;
  Complex Tmp2 = 0.0;
  double TmpLogCoefficient = 0.0;
  double TmpLogCoefficient2 = 0.0;
  double TmpLogCoefficient3 = this->NbrBosons * log((double) this->KyMax);
  double PhaseFactor = 2.0 * M_PI / ((double) this->KyMax);
  if (clockwise == true)
    PhaseFactor *= -1.0;
  ComplexMatrix DeterminantMatrix (this->NbrBosons, this->NbrBosons);
  ComplexMatrix PhaseMatrix (this->KyMax, this->KyMax);
  for (int k = 0; k < this->KyMax; ++k)
    for (int l = 0; l < this->KyMax; ++l)
      PhaseMatrix[k][l] = Phase(PhaseFactor * ((double) (k * l)));
  for (int i = minIndex ; i < LastIndex; ++i)
    {
      this->FermionToBoson(this->StateDescription[i], this->StateKyMax[i] + this->NbrBosons - 1, 
			   this->TemporaryState, this->TemporaryStateKyMax);
      this->ConvertToMonomial(this->StateDescription[i], this->StateKyMax[i] + this->NbrBosons - 1, TmpOutputMonomial);
      Tmp = 0.0;
      TmpLogCoefficient = 0.0;
      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
	TmpLogCoefficient += LogFactorialCoefficients[this->TemporaryState[k]];
      for (int j = 0; j < TmpInputSpace->HilbertSpaceDimension; ++j)
	{
	  this->FermionToBoson(TmpInputSpace->StateDescription[j], TmpInputSpace->StateKyMax[j] + TmpInputSpace->NbrBosons - 1, 
			       TmpInputSpace->TemporaryState, TmpInputSpace->TemporaryStateKyMax);
	  TmpInputSpace->ConvertToMonomial(TmpInputSpace->StateDescription[j], TmpInputSpace->StateKyMax[j] + this->NbrBosons - 1, TmpInputMonomial);
	  for (int k = 0; k < this->NbrBosons; ++k)
	    for (int l = 0; l < this->NbrBosons; ++l)
	      DeterminantMatrix[k][l] = PhaseMatrix[(int) TmpInputMonomial[k]][(int) TmpOutputMonomial[l]];
	  TmpLogCoefficient2 = 0.0;
	  for (int k = 0; k <= TmpInputSpace->TemporaryStateKyMax; ++k)
	    {
	      TmpLogCoefficient2 += LogFactorialCoefficients[TmpInputSpace->TemporaryState[k]];
	    }
	  Tmp += inputState[j] * DeterminantMatrix.Permanent() * exp(0.5 * (-TmpLogCoefficient2 - TmpLogCoefficient - TmpLogCoefficient3));
	}
      outputState[i] = Tmp;      
    }
  delete[] TmpInputMonomial;
  delete[] TmpOutputMonomial;
  delete[] LogFactorialCoefficients;
  return outputState;
}

// symmetrized a product of several uncoupled states 
//
// outputState = reference on the output state
// inputStates = states which will be symmetrized
// inputSpaces = Hilbert spaces attached to each states
// nbrStates = number of states to symmetrize
// firstComponent = first component to symmetrize within the first Hilbert space of inputSpaces
// nbrComponents = number of components to symmetrize within the first Hilbert space of inputSpaces

void BosonOnTorusShort::SymmetrizeU1U1StateCore (ComplexVector& outputState, ComplexVector* inputStates, ParticleOnTorus** inputSpaces, int nbrStates, unsigned long firstComponent, unsigned long nbrComponents)
{
  unsigned long LastComponent = firstComponent + nbrComponents;
  
  FactorialCoefficient* Factorials = new FactorialCoefficient[nbrStates];
  for (long i = (long) firstComponent; i < (long) LastComponent; ++i)
    {
      BosonOnTorusShort* TmpSpace1 = (BosonOnTorusShort*) inputSpaces[0];
      this->FermionToBoson(TmpSpace1->StateDescription[i], TmpSpace1->StateKyMax[i] + TmpSpace1->NbrBosons - 1, 
			   TmpSpace1->TemporaryState, TmpSpace1->TemporaryStateKyMax);
      
      for (int k = TmpSpace1->TemporaryStateKyMax + 1;  k < TmpSpace1->KyMax; ++k)
	TmpSpace1->TemporaryState[k] = 0;
      Complex TmpCoefficient = inputStates[0][i];
      Factorials[0].SetToOne();
      Factorials[0].Power2Divide(TmpSpace1->NbrBosons);
      for (int k = 0; k <= TmpSpace1->TemporaryStateKyMax; ++k)
	if (TmpSpace1->TemporaryState[k] > 1)
	  Factorials[0].FactorialDivide(TmpSpace1->TemporaryState[k]);
      BosonOnTorusShort* TmpSpace2 = (BosonOnTorusShort*) inputSpaces[1];            
      for (long j2 = 0l; j2 < TmpSpace2->LargeHilbertSpaceDimension; ++j2)
	{
	  this->FermionToBoson(TmpSpace2->StateDescription[j2], TmpSpace2->StateKyMax[j2] + TmpSpace2->NbrBosons - 1, 
			       TmpSpace2->TemporaryState, TmpSpace2->TemporaryStateKyMax);
	  if (nbrStates == 2)
	    {	  
	      this->TemporaryStateKyMax = TmpSpace1->TemporaryStateKyMax;
	      if (TmpSpace2->TemporaryStateKyMax > this->TemporaryStateKyMax)
		this->TemporaryStateKyMax = TmpSpace2->TemporaryStateKyMax;
	      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
		this->TemporaryState[k] = 0x0ul;
	      for (int k = 0; k <= TmpSpace1->TemporaryStateKyMax; ++k)
		this->TemporaryState[k] += TmpSpace1->TemporaryState[k];
	      for (int k = 0; k <= TmpSpace2->TemporaryStateKyMax; ++k)
		this->TemporaryState[k] += TmpSpace2->TemporaryState[k];
	      int TmpPos = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
	      if (TmpPos < this->HilbertSpaceDimension)
		{
		  Factorials[1] = Factorials[0];
		  for (int k = 0; k <= TmpSpace2->TemporaryStateKyMax; ++k)
		    if (TmpSpace2->TemporaryState[k] > 1)
		      Factorials[1].FactorialDivide(TmpSpace2->TemporaryState[k]);
		  for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
		    if (this->TemporaryState[k] > 1)
		      Factorials[1].FactorialMultiply(this->TemporaryState[k]);	
		  
		  outputState[TmpPos] += sqrt(Factorials[1].GetNumericalValue()) * TmpCoefficient * inputStates[1][j2];
		}
	    }
	  else
	    {
	      Factorials[1] = Factorials[0];
	      for (int k = 0; k <= TmpSpace2->TemporaryStateKyMax; ++k)
		if (TmpSpace2->TemporaryState[k] > 1)
		  Factorials[1].FactorialDivide(TmpSpace2->TemporaryState[k]);
	      Complex TmpCoefficient2 = TmpCoefficient * inputStates[1][j2];
	      BosonOnTorusShort* TmpSpace3 = (BosonOnTorusShort*) inputSpaces[2];            
	      for (long j3 = 0l; j3 < TmpSpace3->LargeHilbertSpaceDimension; ++j3)
		{
		  this->FermionToBoson(TmpSpace3->StateDescription[j3], TmpSpace3->StateKyMax[j3] + TmpSpace3->NbrBosons - 1, 
				       TmpSpace3->TemporaryState, TmpSpace3->TemporaryStateKyMax);
		  if (nbrStates == 3)
		    {	  
		      this->TemporaryStateKyMax = TmpSpace1->TemporaryStateKyMax;
		      if (TmpSpace2->TemporaryStateKyMax > this->TemporaryStateKyMax)
			this->TemporaryStateKyMax = TmpSpace2->TemporaryStateKyMax;
		      if (TmpSpace3->TemporaryStateKyMax > this->TemporaryStateKyMax)
			this->TemporaryStateKyMax = TmpSpace3->TemporaryStateKyMax;
		      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
			this->TemporaryState[k] = 0x0ul;
		      for (int k = 0; k <= TmpSpace1->TemporaryStateKyMax; ++k)
			this->TemporaryState[k] += TmpSpace1->TemporaryState[k];
		      for (int k = 0; k <= TmpSpace2->TemporaryStateKyMax; ++k)
			this->TemporaryState[k] += TmpSpace2->TemporaryState[k];
		      for (int k = 0; k <= TmpSpace3->TemporaryStateKyMax; ++k)
			this->TemporaryState[k] += TmpSpace3->TemporaryState[k];
		      int TmpPos = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
		      if (TmpPos < this->HilbertSpaceDimension)
			{
			  Factorials[2] = Factorials[1];
			  for (int k = 0; k <= TmpSpace3->TemporaryStateKyMax; ++k)
			    if (TmpSpace3->TemporaryState[k] > 1)
			      Factorials[2].FactorialDivide(TmpSpace3->TemporaryState[k]);
			  for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
			    if (this->TemporaryState[k] > 1)
			      Factorials[2].FactorialMultiply(this->TemporaryState[k]);				  
			  outputState[TmpPos] += sqrt(Factorials[2].GetNumericalValue()) * TmpCoefficient2 * inputStates[2][j3];
			}
		    }
		  else
		    {
		      Factorials[2] = Factorials[1];
		      for (int k = 0; k <= TmpSpace3->TemporaryStateKyMax; ++k)
			if (TmpSpace3->TemporaryState[k] > 1)
			  Factorials[2].FactorialDivide(TmpSpace3->TemporaryState[k]);
		      Complex TmpCoefficient3 = TmpCoefficient2 * inputStates[2][j3];
		      BosonOnTorusShort* TmpSpace4 = (BosonOnTorusShort*) inputSpaces[3];            
		      for (long j4 = 0l; j4 < TmpSpace4->LargeHilbertSpaceDimension; ++j4)
			{
			  this->FermionToBoson(TmpSpace4->StateDescription[j4], TmpSpace4->StateKyMax[j4] + TmpSpace4->NbrBosons - 1, 
					       TmpSpace4->TemporaryState, TmpSpace4->TemporaryStateKyMax);
			  if (nbrStates == 4)
			    {	  
			      this->TemporaryStateKyMax = TmpSpace1->TemporaryStateKyMax;
			      if (TmpSpace2->TemporaryStateKyMax > this->TemporaryStateKyMax)
				this->TemporaryStateKyMax = TmpSpace2->TemporaryStateKyMax;
			      if (TmpSpace3->TemporaryStateKyMax > this->TemporaryStateKyMax)
				this->TemporaryStateKyMax = TmpSpace3->TemporaryStateKyMax;
			      if (TmpSpace4->TemporaryStateKyMax > this->TemporaryStateKyMax)
				this->TemporaryStateKyMax = TmpSpace4->TemporaryStateKyMax;
			      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
				this->TemporaryState[k] = 0x0ul;
			      for (int k = 0; k <= TmpSpace1->TemporaryStateKyMax; ++k)
				this->TemporaryState[k] += TmpSpace1->TemporaryState[k];
			      for (int k = 0; k <= TmpSpace2->TemporaryStateKyMax; ++k)
				this->TemporaryState[k] += TmpSpace2->TemporaryState[k];
			      for (int k = 0; k <= TmpSpace3->TemporaryStateKyMax; ++k)
				this->TemporaryState[k] += TmpSpace3->TemporaryState[k];
			      for (int k = 0; k <= TmpSpace4->TemporaryStateKyMax; ++k)
				this->TemporaryState[k] += TmpSpace4->TemporaryState[k];
			      int TmpPos = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), this->TemporaryStateKyMax + this->NbrBosons - 1);
			      if (TmpPos < this->HilbertSpaceDimension)
				{
				  Factorials[3] = Factorials[2];
				  for (int k = 0; k <= TmpSpace4->TemporaryStateKyMax; ++k)
				    if (TmpSpace4->TemporaryState[k] > 1)
				      Factorials[3].FactorialDivide(TmpSpace4->TemporaryState[k]);
				  for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
				    if (this->TemporaryState[k] > 1)
				      Factorials[3].FactorialMultiply(this->TemporaryState[k]);	
				  outputState[TmpPos] += sqrt(Factorials[3].GetNumericalValue()) * TmpCoefficient3 * inputStates[3][j4];
				}
			    }
			  else
			    {
			      cout << "error, current code cannot symmetrized more than  states" << endl;
			    }
			}
		    }
		}
	    }
	}
    }  
}
  
// create a state from its MPS description
//
// bMatrices = array that gives the B matrices 
// twistMatrix = reference on the twist matrix to insert in the trace
// state = reference to vector that will contain the state description
// mPSSumIndices = diagonal indices that have to be kept in the trace
// nbrMPSSumIndices = number of diagonal indices that have to be kept in the trace
// memory = amount of memory that can be use to precompute matrix multiplications  
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void BosonOnTorusShort::CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, SparseRealMatrix& twistMatrix, RealVector& state, 
						       int* mPSSumIndices, int nbrMPSSumIndices,
						       long memory, long initialIndex, long nbrComponents)
{
  SparseRealMatrix TmpMatrix;
  long MaxIndex = initialIndex + nbrComponents;
  if ((nbrComponents == 0l) || (MaxIndex > this->LargeHilbertSpaceDimension))
    {
      MaxIndex = this->LargeHilbertSpaceDimension;
    }
  long TmpMemory = (((long) bMatrices[0].GetNbrRow()) * ((long)bMatrices[0].GetNbrColumn()));
  if (TmpMemory > (1l << 28))
    TmpMemory = 1l << 28;
  double* TmpMatrixElements = new double [TmpMemory];
  int* TmpColumnIndices = new int [TmpMemory];
  double* TmpElements = new double [bMatrices[0].GetNbrRow()];

  if (memory <= 1l)
    {
      for (long i = initialIndex; i < MaxIndex; ++i)
	{
	  if (((i - initialIndex) % 10000) == 0)
	    cout << "Completed " << (i - initialIndex) << " out of " << (MaxIndex - initialIndex) << endl; 
	  this->FermionToBoson(this->StateDescription[i], this->StateKyMax[i] + this->NbrBosons - 1, 
			       this->TemporaryState, this->TemporaryStateKyMax);
	  TmpMatrix.Copy(twistMatrix);
	  for (int j = 0; j <= this->TemporaryStateKyMax; ++j)
	    {
	      TmpMatrix.Multiply(bMatrices[this->TemporaryState[j]], TmpMatrixElements, TmpColumnIndices, TmpElements);
	    } 
	  for (int j = this->TemporaryStateKyMax + 1; j < this->KyMax; ++j)
	    TmpMatrix.Multiply(bMatrices[0], TmpMatrixElements, TmpColumnIndices, TmpElements);
	  double& TmpComponent = state[i];
	  double Tmp;
	  for (int j = 0; j < nbrMPSSumIndices; ++j)
	    {
	      TmpMatrix.GetMatrixElement(mPSSumIndices[j], mPSSumIndices[j], Tmp);
	      TmpComponent += Tmp;
	    }
	}
    }
  else
    {
//       int PrecalculationBlockSize = (int) memory;
//       unsigned long PrecalculationBlockMask = (0x1ul  << PrecalculationBlockSize) - 0x1ul;
//       int PrecalculationBlockLength  = this->KyMax - (this->KyMax % PrecalculationBlockSize);
//       int RemaingOrbtals = this->KyMax % PrecalculationBlockSize;
//       SparseRealMatrix* TmpBlockbMatrices = new SparseRealMatrix[1 << PrecalculationBlockSize];
//       for (unsigned long TmpState = 0x0ul; TmpState <= PrecalculationBlockMask; ++TmpState)
// 	{
// 	  unsigned long TmpStateDescription = TmpState;
// 	  TmpBlockbMatrices[TmpState].Copy(bMatrices[TmpStateDescription & 0x1ul]);
// 	  TmpStateDescription >>= 1;
// 	  for (int j = 1; j < PrecalculationBlockSize; ++j)
// 	    {
// 	      TmpBlockbMatrices[TmpState].Multiply(bMatrices[TmpStateDescription & 0x1ul], TmpMatrixElements, TmpColumnIndices, TmpElements);
// 	      TmpStateDescription >>= 1;
// 	    }      
// 	}
//       for (long i = initialIndex; i < MaxIndex; ++i)
// 	{
// 	  if (((i - initialIndex) % 10000) == 0)
// 	    cout << "Completed " << (i - initialIndex) << " out of " << (MaxIndex - initialIndex) << endl; 
// 	  unsigned long TmpStateDescription = this->StateDescription[i];
// 	  TmpMatrix.Copy(twistMatrix);
// 	  int j = 0;
// 	  for (; j < PrecalculationBlockLength; j += PrecalculationBlockSize)
// 	    {
// 	      TmpMatrix.Multiply(TmpBlockbMatrices[TmpStateDescription & PrecalculationBlockMask], TmpMatrixElements, TmpColumnIndices, TmpElements);
// 	      TmpStateDescription >>= PrecalculationBlockSize;
// 	    } 
// 	  for (; j < this->KyMax; ++j)
// 	    {
// 	      TmpMatrix.Multiply(bMatrices[TmpStateDescription & 0x1ul], TmpMatrixElements, TmpColumnIndices, TmpElements);
// 	      TmpStateDescription >>= 1;
// 	    } 
// 	  double& TmpComponent = state[i];
// 	  double Tmp;
// 	  for (int j = 0; j < nbrMPSSumIndices; ++j)
// 	    {
// 	      TmpMatrix.GetMatrixElement(mPSSumIndices[j], mPSSumIndices[j], Tmp);
// 	      TmpComponent += Tmp;
// 	    }
// 	}
    }
  delete[] TmpMatrixElements;
  delete[] TmpColumnIndices;
  delete[] TmpElements;
}


// request whether state with given index satisfies a general Pauli exclusion principle
//
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals
// return value = true if teh state satisfies the general Pauli exclusion principle

bool BosonOnTorusShort::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  this->FermionToBoson(this->StateDescription[index], this->StateKyMax[index] + this->NbrBosons - 1, 
		       this->TemporaryState, this->TemporaryStateKyMax);  
  for (int i = this->TemporaryStateKyMax + 1;  i < this->KyMax; ++i)
    this->TemporaryState[i] = 0x0ul;
  unsigned long UnsignedPauliK = (unsigned long) pauliK;
  for (int i = 0; i < this->KyMax; ++i)
    {
      unsigned long TmpOccupation = 0l;
      for (int j = 0; j < pauliR; ++j)
	{
	  TmpOccupation += this->TemporaryState[(i + j) % this->KyMax];
	}
      if (TmpOccupation > UnsignedPauliK)
	return false;
    }
  return true;
}
  
