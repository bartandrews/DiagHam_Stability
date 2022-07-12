////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of particle on a torus                       //
//                                                                            //
//                        last modification : 18/07/2002                      //
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
#include "HilbertSpace/FermionOnTorus.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "HilbertSpace/FermionOnDisk.h"
#include "HilbertSpace/FermionOnSphereFull.h"
#include "Matrix/Matrix.h"
#include "GeneralTools/ArrayTools.h"
#include "MathTools/IntegerAlgebraTools.h"

#include <math.h>
#include <algorithm>
#include <set>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


// basic constructor
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion

FermionOnTorus::FermionOnTorus (int nbrFermions, int maxMomentum)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->KyMax = maxMomentum;
  this->NbrLzValue = this->KyMax + 1;

  this->MomentumModulo = FindGCD(this->NbrFermions, this->KyMax);
  this->StateShift = this->KyMax / this->MomentumModulo;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->KyMax - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->KyMax);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  this->Flag.Initialize();
  this->TotalKy = 0;
  this->TotalKyFlag = false;
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateKyMax = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrFermions, this->KyMax - 1, this->KyMax - 1, 0);
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
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
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// momentumConstraint = index of the momentum orbit

FermionOnTorus::FermionOnTorus (int nbrFermions, int maxMomentum, int momentumConstraint)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->KyMax = maxMomentum;
  this->NbrLzValue = this->KyMax + 1;
  this->TotalKy = momentumConstraint;

  this->MomentumModulo = FindGCD(this->NbrFermions, this->KyMax);
  this->StateShift = this->KyMax / this->MomentumModulo;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->KyMax - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;
  this->TotalKyFlag = true;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->KyMax);
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateKyMax = new int [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->KyMax - 1, this->KyMax - 1, 0, 0);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
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
// fermions = reference on the hilbert space to copy to copy

FermionOnTorus::FermionOnTorus(const FermionOnTorus& fermions)
{
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateKyMax = fermions.StateKyMax;
  this->KyMax = fermions.KyMax;
  this->NbrLzValue = fermions.NbrLzValue;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->TotalKy = fermions.TotalKy;
  this->TotalKyFlag = fermions.TotalKyFlag;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->Flag = fermions.Flag;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
}

// constructor from full datas (with no constraint on the total momentum)
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// hilbertSpaceDimension = Hilbert space dimension
// stateDescription = array describing each state
// stateKyMax = array giving maximum Lz value reached for a fermion in a given state

FermionOnTorus::FermionOnTorus (int nbrFermions, int maxMomentum, int hilbertSpaceDimension, 
				unsigned long* stateDescription, int* stateKyMax)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->KyMax = maxMomentum;
  this->NbrLzValue = this->KyMax + 1;

  this->MomentumModulo = FindGCD(this->NbrFermions, this->KyMax);
  this->StateShift = this->KyMax / this->MomentumModulo;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->KyMax - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->TotalKyFlag = false;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = stateDescription;
  this->StateKyMax = stateKyMax;
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
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
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// momentumConstraint = index of the momentum orbit
// hilbertSpaceDimension = Hilbert space dimension
// stateDescription = array describing each state
// stateKyMax = array giving maximum Lz value reached for a fermion in a given state

FermionOnTorus::FermionOnTorus (int nbrFermions, int maxMomentum, int momentumConstraint, int hilbertSpaceDimension, 
				unsigned long* stateDescription, int* stateKyMax)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->KyMax = maxMomentum;
  this->NbrLzValue = this->KyMax + 1;
  this->TotalKy = momentumConstraint;
  this->TotalKyFlag = true;

  this->MomentumModulo = FindGCD(this->NbrFermions, this->KyMax);
  this->StateShift = this->KyMax / this->MomentumModulo;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->KyMax - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = stateDescription;
  this->StateKyMax = stateKyMax;
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
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


// destructor
//

FermionOnTorus::~FermionOnTorus ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateKyMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
}

// assignment (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnTorus& FermionOnTorus::operator = (const FermionOnTorus& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateKyMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
  if (this->TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateKyMax = fermions.StateKyMax;
  this->KyMax = fermions.KyMax;
  this->NbrLzValue = fermions.NbrLzValue;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->TotalKy = fermions.TotalKy;
  this->TotalKyFlag = fermions.TotalKyFlag;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->Flag = fermions.Flag;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnTorus::Clone()
{
  return new FermionOnTorus(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnTorus::SetTargetSpace(ParticleOnTorus* targetSpace)
{
  printf("Setting target space in FermionOnTorus\n");
  this->TargetSpace = (FermionOnTorus*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int FermionOnTorus::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}


// return number of particles of the target space
//
// return value = Hilbert space dimension

int FermionOnTorus::GetTargetNbrParticles()
{
  return this->TargetSpace->GetNbrParticles();
}


// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnTorus::GetQuantumNumbers ()
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

AbstractQuantumNumber* FermionOnTorus::GetQuantumNumber (int index)
{
  if (this->TotalKyFlag == false)
    {
      return  new PeriodicMomentumQuantumNumber (this->GetMomentumValue(index), this->KyMax);
    }
  else
    return new PeriodicMomentumQuantumNumber (this->TotalKy, this->KyMax);
}

// get momemtum value of a given state
//
// index = state index
// return value = state momentum

int FermionOnTorus::GetMomentumValue(int index)
{
  int StateKyMax = this->StateKyMax[index];
  unsigned long State = this->StateDescription[index];
  int Momentum = 0;
  for (int i = 0; i <= StateKyMax; ++i)
    {
      Momentum += ((State >> i ) & ((unsigned long) 1)) * i;
    }
  return (Momentum % this->KyMax);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnTorus::ExtractSubspace (AbstractQuantumNumber& q, 
						       SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != AbstractQuantumNumber::PeriodicMomentum)
    return 0;
  if (this->TotalKyFlag == true)
    {
      if (this->TotalKy == ((PeriodicMomentumQuantumNumber&) q).GetMomentum())
	return this;
      else 
	return 0;
    }
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
  unsigned long* SubspaceStateDescription = new unsigned long [SubspaceHilbertSpaceDimension];
  int* SubspaceStateKyMax = new int [SubspaceHilbertSpaceDimension];
  int* ConvArray = new int [SubspaceHilbertSpaceDimension];
  for (int i = 0; i < SubspaceHilbertSpaceDimension; ++i)
    {
      ConvArray[i] = TmpConvArray[i];
      SubspaceStateDescription[i] = this->StateDescription[TmpConvArray[i]];
      SubspaceStateKyMax[i] = this->StateKyMax[TmpConvArray[i]];
    }
  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, SubspaceHilbertSpaceDimension, ConvArray);
  return new FermionOnTorus(this->NbrFermions, this->KyMax, Momentum, SubspaceHilbertSpaceDimension, 
			    SubspaceStateDescription, SubspaceStateKyMax);
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2[KyMax])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorus::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateKyMax = this->StateKyMax[index];
  unsigned long State = this->StateDescription[index];
  if ((n1 > StateKyMax) || (n2 > StateKyMax) || ((State & (0x1ul << n1)) == 0) || 
      ((State & (0x1ul << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewKyMax = StateKyMax;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(0x1ul << n2);
  if (NewKyMax == n2)
    {
      if (TmpState != 0x0ul)
	{
	  while ((TmpState >> NewKyMax) == 0)
	    --NewKyMax;
	}
      else
	{
	  NewKyMax = 0;
	}
    }
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(0x1ul << n1);
  if (NewKyMax == n1)
    {
      if (TmpState != 0x0ul)
	{
	  while ((TmpState >> NewKyMax) == 0)
	    --NewKyMax;
	}
      else
	{
	  NewKyMax = 0;
	}
    }
  if ((TmpState & (0x1ul << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewKyMax)
    {
      NewKyMax = m2;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if ((TmpState & (0x1ul << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewKyMax)
    {
      NewKyMax = m1;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewKyMax);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnTorus::ProdA (int index, int* n, int nbrIndices)
{
  this->ProdALzMax = this->StateKyMax[index];
  this->ProdATemporaryState = this->StateDescription[index];
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = n[i];
      if ((this->ProdATemporaryState & (0x1l << Index)) == 0)
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      this->ProdATemporaryState &= ~(0x1l << Index);
    }
  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdALzMax = 0;
      return Coefficient;      
    }
  while (((this->ProdATemporaryState >> this->ProdALzMax) == 0) && (this->ProdALzMax > 0))
    --this->ProdALzMax;

  return Coefficient;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnTorus::AA (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];

  if (((ProdATemporaryState & (((unsigned long) (0x1)) << n1)) == 0) 
      || ((ProdATemporaryState & (((unsigned long) (0x1)) << n2)) == 0) || (n1 == n2))
    return 0.0;

  this->ProdALzMax = this->StateKyMax[index];

  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n1);

  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdALzMax = 0;
      return Coefficient;      
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorus::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = m[i];
      if ((TmpState & (0x1l << Index)) != 0)
	{
	  coefficient = 0.0;
	  return this->TargetSpace->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (0x1l << Index);
    }
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorus::AdAd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if ((TmpState & (((unsigned long) (0x1)) << m2))!= 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m2);
  if ((TmpState & (((unsigned long) (0x1)) << m1))!= 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {      
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m a_m operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for creation operator
// return value =  resulting multiplicative factor 

double FermionOnTorus::AdA (int index, int m)
{
  if (this->StateKyMax[index] < m)
    return 0.0;
  if ((this->StateDescription[index] & (0x1ul << m)) == 0)
    return 0.0;
  else
    return 1.0;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

// attention: check sign returned by this function!
int FermionOnTorus::AdA (int index, int m, int n, double& coefficient)
{
  int StateLzMax = this->StateKyMax[index];
  unsigned long State = this->StateDescription[index];
  if ((n > StateLzMax) || ((State & (((unsigned long) (0x1)) << n)) == 0))
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  int NewLzMax = StateLzMax;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n);
  if ((TmpState != 0x0ul))
    {
      while ((TmpState >> NewLzMax) == 0)
	--NewLzMax;
    }
  else
    NewLzMax = 0;
  if ((TmpState & (((unsigned long) (0x1)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  if (m > NewLzMax)
    {
      NewLzMax = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }

  
  TmpState |= (((unsigned long) (0x1)) << m);
  
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// return matrix representation of the annihilation operator a_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& FermionOnTorus::A (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  int StateKyMax;
  unsigned long State;
  double Coefficient;
  unsigned long GlobalMask = (0x1ul << i);
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)
    {
      StateKyMax = this->StateKyMax[j];
      if (StateKyMax >= i)
	{
	  State = this->StateDescription[j];
	  if ((State & GlobalMask) != 0)
	    {
	      Coefficient = this->SignLookUpTable[(State >> i) & this->SignLookUpTableMask[i]];
	      Coefficient *= this->SignLookUpTable[(State >> (i + 16))  & this->SignLookUpTableMask[i + 16]];
#ifdef  __64_BITS__
	      Coefficient *= this->SignLookUpTable[(State >> (i + 32)) & this->SignLookUpTableMask[i + 32]];
	      Coefficient *= this->SignLookUpTable[(State >> (i + 48)) & this->SignLookUpTableMask[i + 48]];
#endif
	      State &= ~GlobalMask;
	      if (StateKyMax == i)
		while ((State >> StateKyMax) == 0)
		  --StateKyMax;
	      M(this->FindStateIndex(State, StateKyMax), j) = Coefficient;
	    }
	}
    }
  return M;
}

// return matrix representation of the creation operator a^+_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& FermionOnTorus::Ad (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  int StateKyMax;
  unsigned long State;
  double Coefficient;
  unsigned long GlobalMask = (0x1ul << i);
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)
    {
      StateKyMax = this->StateKyMax[j];
      State = this->StateDescription[j];
      if ((State & GlobalMask) == 0)
	{
	  if (i > StateKyMax)
	    {
	      State |= GlobalMask;
	      M(this->FindStateIndex(State, i), j) = 1.0;
	    }
	  else
	    {
	      Coefficient = this->SignLookUpTable[(State >> i) & this->SignLookUpTableMask[i]];
	      Coefficient *= this->SignLookUpTable[(State >> (i + 16))  & this->SignLookUpTableMask[i + 16]];
#ifdef  __64_BITS__
	      Coefficient *= this->SignLookUpTable[(State >> (i + 32)) & this->SignLookUpTableMask[i + 32]];
	      Coefficient *= this->SignLookUpTable[(State >> (i + 48)) & this->SignLookUpTableMask[i + 48]];
#endif
	      State |= GlobalMask;
	      M(this->FindStateIndex(State, StateKyMax), j) = Coefficient;
	    }
	}
    }
  return M;
}


// apply an annihilation operator a_i and return the index in the target space
//
// i = state index
// n = index of annihilation operator
// coefficient = will be multiplied by the prefactor of the bosonic ladder operator
// return value = index in the target space
int FermionOnTorus::A (int index, int n, double &coefficient)
{
  int StateLzMax = this->StateKyMax[index];
  unsigned long State = this->StateDescription[index];
  if ((n > StateLzMax) || ((State & (((unsigned long) (0x1)) << n)) == 0))
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  int NewLzMax = StateLzMax;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n);
  if ((TmpState != 0x0ul))
    {
      while ((TmpState >> NewLzMax) == 0)
	--NewLzMax;
    }
  else
    NewLzMax = 0;  
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply a creation operator a_i and return the index in the target space
//
// i = state index
// m = index of annihilation operator
// coefficient = will be multiplied by the prefactor of the bosonic ladder operator
// return value = index in the target space
int FermionOnTorus::Ad (int index, int m, double &coefficient)
{
  int NewLzMax = this->StateKyMax[index];
  unsigned long TmpState = this->StateDescription[index];
  if ((TmpState & (((unsigned long) (0x1)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  if (m > NewLzMax)
    {
      NewLzMax = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }

  TmpState |= (((unsigned long) (0x1)) << m);  
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}


// find state index
//
// stateDescription = unsigned longeger describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnTorus::FindStateIndex(unsigned long stateDescription, int maxMomentum)
{
//  int Pos = 0;
  int Pos = this->LookUpTable[maxMomentum][stateDescription >> this->LookUpTableShift[maxMomentum]];
  while (this->StateDescription[Pos] != stateDescription)
    ++Pos;
  return Pos;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnTorus::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->KyMax; ++i)
    Str << ((TmpState >> i) & 0x1ul) << " ";
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentKyMax = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int FermionOnTorus::GenerateStates(int nbrFermions, int maxMomentum, int currentKyMax, int pos)
{
  if ((nbrFermions == 0) || (nbrFermions > (currentKyMax + 1)) || (currentKyMax < 0))
    return pos;
  if (nbrFermions == 1)
    {
      for (int i = currentKyMax; i >= 0; --i)
	{
	  this->StateDescription[pos] = 0x1ul << i;
	  this->StateKyMax[pos] = maxMomentum;
	  ++pos;
	}
      return pos;
    }
  int ReducedCurrentKyMax = currentKyMax - 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentKyMax, pos);
  unsigned long Mask = 1 << currentKyMax;
  for (int i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (maxMomentum == currentKyMax)
    return this->GenerateStates(nbrFermions, ReducedCurrentKyMax, ReducedCurrentKyMax, TmpPos);
  else
    return this->GenerateStates(nbrFermions, maxMomentum, ReducedCurrentKyMax, TmpPos);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentKyMax = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

int FermionOnTorus::GenerateStates(int nbrFermions, int maxMomentum, int currentKyMax, int pos, int currentMomentum)
{
  if ((nbrFermions == 0) || (nbrFermions > (currentKyMax + 1)))
    return pos;
  if (nbrFermions == 1)
    {
      int i = this->TotalKy - (currentMomentum % this->KyMax);
      if (i < 0)
	i += this->KyMax;
      for (; i <= currentKyMax; i += this->KyMax)
	{
	  this->StateDescription[pos] = 0x1ul << i;
	  this->StateKyMax[pos] = maxMomentum;
	  ++pos;
	}
      return pos;
    }
  int ReducedCurrentKyMax = currentKyMax - 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentKyMax, pos, currentMomentum + currentKyMax);
  unsigned long Mask = 0x1ul << currentKyMax;
  for (int i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (maxMomentum == currentKyMax)
    return this->GenerateStates(nbrFermions, ReducedCurrentKyMax, ReducedCurrentKyMax, TmpPos, currentMomentum);
  else
    return this->GenerateStates(nbrFermions, maxMomentum, ReducedCurrentKyMax, TmpPos, currentMomentum);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnTorus::GenerateLookUpTable(int memory)
{
  // evaluate look-up table size
  memory /= (4 * this->NbrLzValue);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrLzValue)
    this->MaximumLookUpShift = this->NbrLzValue;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrLzValue];
  this->LookUpTableShift = new int [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize];
  int CurrentKyMax = this->StateKyMax[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentKyMax];
  if (CurrentKyMax < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentKyMax] = 0;
  else
    this->LookUpTableShift[CurrentKyMax] = CurrentKyMax + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentKyMax];
  unsigned long CurrentLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  unsigned long TmpLookUpTableValue;
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentKyMax != this->StateKyMax[i])
	{
 	  CurrentKyMax = this->StateKyMax[i];
	  TmpLookUpTable = this->LookUpTable[CurrentKyMax];
	  if (CurrentKyMax < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentKyMax] = 0;
	  else
	    this->LookUpTableShift[CurrentKyMax] = CurrentKyMax + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentKyMax];
	  CurrentLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      CurrentLookUpTableValue = TmpLookUpTableValue;
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  
  // look-up tables for evaluating sign when applying creation/annihilation operators
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      Count = 0;
      TmpNbr = j;
      while (TmpNbr != 0)
	{
	  if (TmpNbr & 0x1)
	    ++Count;
	  TmpNbr >>= 1;
	}
      if (Count & 1)
	this->SignLookUpTable[j] = -1.0;
      else
	this->SignLookUpTable[j] = 1.0;
    }
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#endif
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// return value = Hilbert space dimension

int FermionOnTorus::EvaluateHilbertSpaceDimension(int nbrFermions, int maxMomentum)
{
  FactorialCoefficient Dimension; 
  Dimension.PartialFactorialMultiply(maxMomentum - nbrFermions + 1, maxMomentum); 
  Dimension.FactorialDivide(nbrFermions);
  return Dimension.GetIntegerValue();
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix FermionOnTorus::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int kySector, RealVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrFermionSector == 0))
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
      if (((kySector % this->KyMax) == this->TotalKy) && (nbrFermionSector == this->NbrFermions))
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

  if (nbrFermionSector == 0)
    {
      if (kySector == 0)
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  double Coefficient = 0.0;
	  unsigned long Mask  = (0x1ul << subsytemSize) - 1ul;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
                Coefficient += groundState[i] * groundState[i];
            }
	  TmpDensityMatrix.SetMatrixElement(0, 0, Coefficient);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  if (nbrFermionSector == this->NbrFermions)
    {
      if ((kySector % this->KyMax) == this->TotalKy)
	{
	  FermionOnDisk TmpDestinationHilbertSpace(nbrFermionSector, kySector, subsytemSize - 1);
	  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
	  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	    {
	      int TmpPos = this->FindStateIndex(TmpDestinationHilbertSpace.StateDescription[i], TmpDestinationHilbertSpace.StateLzMax[i]);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpDensityMatrix.AddToMatrixElement(i, i, groundState[TmpPos] * groundState[TmpPos]);
		  for (int j = i + 1; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
		    {
		      int TmpPos2 = this->FindStateIndex(TmpDestinationHilbertSpace.StateDescription[j], TmpDestinationHilbertSpace.StateLzMax[j]);
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
  FermionOnDisk TmpDestinationHilbertSpace(nbrFermionSector, kySector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrFermions - nbrFermionSector) * subsytemSize);
  int ComplementaryMinTotalKy = ((this->NbrFermions - nbrFermionSector) * (this->NbrFermions - nbrFermionSector - 1)) / 2;
  while (ComplementaryTotalKy < ComplementaryMinTotalKy)
    ComplementaryTotalKy += this->KyMax;
  //  ComplementaryTotalKy = ComplementaryTotalKy % this->KyMax;
  int ComplementaryMaxTotalKy = ((2 * (this->KyMax - subsytemSize - 1) + 1 - (this->NbrFermions - nbrFermionSector)) * (this->NbrFermions - nbrFermionSector)) / 2;    

  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      FermionOnDisk TmpHilbertSpace(this->NbrFermions - nbrFermionSector, ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  int Pos = 0;
	  unsigned long TmpComplementaryState = TmpHilbertSpace.StateDescription[MinIndex] << subsytemSize;
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax + this->NbrFermions - 1;
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
      ComplementaryTotalKy += this->KyMax;
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

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Ky sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnTorus::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int kySector, ComplexVector& groundState)
{
  if (nbrFermionSector == 0)
    {
      if (kySector == 0)
	{
	  HermitianMatrix TmpDensityMatrix(1);
	  Complex Coefficient = 0.0;
	  unsigned long Mask  = (0x1ul << subsytemSize) - 1ul;
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

  if (nbrFermionSector == this->NbrFermions)
    {
      if ((kySector % this->KyMax) == this->TotalKy)
	{
	  FermionOnDisk TmpDestinationHilbertSpace(nbrFermionSector, kySector, subsytemSize - 1);
	  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
	  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	    {
	      int TmpPos = this->FindStateIndex(TmpDestinationHilbertSpace.StateDescription[i], TmpDestinationHilbertSpace.StateLzMax[i]);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpDensityMatrix.AddToMatrixElement(i, i, groundState[TmpPos] * Conj(groundState[TmpPos]));
		  for (int j = i + 1; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
		    {
		      int TmpPos2 = this->FindStateIndex(TmpDestinationHilbertSpace.StateDescription[j], TmpDestinationHilbertSpace.StateLzMax[j]);
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
  FermionOnDisk TmpDestinationHilbertSpace(nbrFermionSector, kySector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrFermions - nbrFermionSector) * subsytemSize);
  int ComplementaryMinTotalKy = ((this->NbrFermions - nbrFermionSector) * (this->NbrFermions - nbrFermionSector - 1)) / 2;
  while (ComplementaryTotalKy < ComplementaryMinTotalKy)
    ComplementaryTotalKy += this->KyMax;
  int ComplementaryMaxTotalKy = ((2 * (this->KyMax - subsytemSize - 1) + 1 - (this->NbrFermions - nbrFermionSector)) * (this->NbrFermions - nbrFermionSector)) / 2;    

  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      FermionOnDisk TmpHilbertSpace(this->NbrFermions - nbrFermionSector, ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  int Pos = 0;
	  unsigned long TmpComplementaryState = TmpHilbertSpace.StateDescription[MinIndex] << subsytemSize;
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax + this->NbrFermions - 1;
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
      ComplementaryTotalKy += this->KyMax;
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

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Ky sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem

RealMatrix FermionOnTorus::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int kySector, RealVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrFermionSector == 0))
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
      if ((kySector == this->TotalKy) && (nbrFermionSector == this->NbrFermions))
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
  if (nbrFermionSector == 0)
    {
      if (kySector == 0)
	{
	  double Coefficient = 0.0;
	  unsigned long Mask  = (0x1ul << subsytemSize) - 1ul;
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
  if (nbrFermionSector == this->NbrFermions)
    {
      if ((kySector % this->KyMax) == this->TotalKy)
	{
	  unsigned long Mask  = ((0x1ul << (this->KyMax - subsytemSize)) - 1ul) << subsytemSize;
	  int TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) ==  0x0ul)
                ++TmpSize;
            }
	  if (TmpSize == 0)
	    {
	      RealMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
	    }
	  RealMatrix TmpEntanglementMatrix(TmpSize, 1, true);
	  TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
		{
 		  TmpEntanglementMatrix.AddToMatrixElement(TmpSize, 0, groundState[i]);
		  ++TmpSize;
		}
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
  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrFermions - nbrFermionSector) * subsytemSize);
  int ComplementaryMinTotalKy = ((this->NbrFermions - nbrFermionSector) * (this->NbrFermions - nbrFermionSector - 1)) / 2;
  while (ComplementaryTotalKy < ComplementaryMinTotalKy)
    ComplementaryTotalKy += this->KyMax;
  int ComplementaryMaxTotalKy = ((2 * (this->KyMax - subsytemSize - 1) + 1 - (this->NbrFermions - nbrFermionSector)) * (this->NbrFermions - nbrFermionSector)) / 2;    
  int NbrComplementarySpaces = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMax;
    }
  if (NbrComplementarySpaces == 0)
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;	  
    }
  FermionOnDisk** ComplementaryHilbertSpaces = new FermionOnDisk*[NbrComplementarySpaces];
  NbrComplementarySpaces = 0;
  ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrFermions - nbrFermionSector) * subsytemSize);
  while (ComplementaryTotalKy < ComplementaryMinTotalKy)
    ComplementaryTotalKy += this->KyMax;
  int TotalComplementaryHilbertSpaceDimension = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ComplementaryHilbertSpaces[NbrComplementarySpaces] = new FermionOnDisk(this->NbrFermions - nbrFermionSector, 
									     ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[NbrComplementarySpaces]->GetHilbertSpaceDimension();
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMax;
    }
 
  FermionOnDisk TmpDestinationHilbertSpace(nbrFermionSector, kySector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,
				   TotalComplementaryHilbertSpaceDimension, true);
  TotalComplementaryHilbertSpaceDimension = 0;
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      for (int MinIndex = 0; MinIndex < ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpComplementaryState = ComplementaryHilbertSpaces[SpaceIndex]->StateDescription[MinIndex] << subsytemSize;
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax;
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
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem

ComplexMatrix FermionOnTorus::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int kySector, ComplexVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrFermionSector == 0))
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
      if ((kySector == this->TotalKy) && (nbrFermionSector == this->NbrFermions))
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
  if (nbrFermionSector == 0)
    {
      if (kySector == 0)
	{
	  unsigned long Mask  = (0x1ul << subsytemSize) - 1ul;
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
            }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  if (nbrFermionSector == this->NbrFermions)
    {
      if ((kySector % this->KyMax) == this->TotalKy)
	{
	  unsigned long Mask  = ((0x1ul << (this->KyMax - subsytemSize)) - 1ul) << subsytemSize;
	  int TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) ==  0x0ul)
                ++TmpSize;
            }
	  if (TmpSize == 0)
	    {
	      ComplexMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
	    }
	  ComplexMatrix TmpEntanglementMatrix(TmpSize, 1, true);
	  TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
		{
 		  TmpEntanglementMatrix.AddToMatrixElement(TmpSize, 0, groundState[i]);
		  ++TmpSize;
		}
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
  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrFermions - nbrFermionSector) * subsytemSize);
  int ComplementaryMinTotalKy = ((this->NbrFermions - nbrFermionSector) * (this->NbrFermions - nbrFermionSector - 1)) / 2;
  while (ComplementaryTotalKy < ComplementaryMinTotalKy)
    ComplementaryTotalKy += this->KyMax;
  int ComplementaryMaxTotalKy = ((2 * (this->KyMax - subsytemSize - 1) + 1 - (this->NbrFermions - nbrFermionSector)) * (this->NbrFermions - nbrFermionSector)) / 2;    
  int NbrComplementarySpaces = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMax;
    }
  if (NbrComplementarySpaces == 0)
    {
      ComplexMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;	  
    }
  FermionOnDisk** ComplementaryHilbertSpaces = new FermionOnDisk*[NbrComplementarySpaces];
  NbrComplementarySpaces = 0;
  ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrFermions - nbrFermionSector) * subsytemSize);
  while (ComplementaryTotalKy < ComplementaryMinTotalKy)
    ComplementaryTotalKy += this->KyMax;
  int TotalComplementaryHilbertSpaceDimension = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ComplementaryHilbertSpaces[NbrComplementarySpaces] = new FermionOnDisk(this->NbrFermions - nbrFermionSector, 
									     ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[NbrComplementarySpaces]->GetHilbertSpaceDimension();
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMax;
    }
 
  FermionOnDisk TmpDestinationHilbertSpace(nbrFermionSector, kySector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  ComplexMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,
				   TotalComplementaryHilbertSpaceDimension, true);
  TotalComplementaryHilbertSpaceDimension = 0;
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      for (int MinIndex = 0; MinIndex < ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpComplementaryState = ComplementaryHilbertSpaces[SpaceIndex]->StateDescription[MinIndex] << subsytemSize;
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax;
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
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem

RealMatrix FermionOnTorus::EvaluatePartialEntanglementMatrixFullKyPartB (int subsytemSize, int nbrFermionSector, int kySector, RealVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrFermionSector == 0))
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
      if ((kySector == this->TotalKy) && (nbrFermionSector == this->NbrFermions))
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
  if (nbrFermionSector == 0)
    {
      if (kySector == 0)
	{
	  unsigned long Mask  = (0x1ul << subsytemSize) - 1ul;
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
	  FermionOnSphereFull TmpFullComplementaryHilbertSpace(this->NbrFermions - nbrFermionSector, this->KyMax - subsytemSize - 1);
	  RealMatrix TmpEntanglementMatrix(1, TmpFullComplementaryHilbertSpace.HilbertSpaceDimension, true);
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
		{
		  unsigned long TmpState = this->StateDescription[i] >> subsytemSize;
		  int TmpKyMax = this->KyMax - subsytemSize - 1;
		  while ((TmpState >> TmpKyMax) == 0x0ul)
		    --TmpKyMax;
 		  TmpEntanglementMatrix.AddToMatrixElement(0, TmpFullComplementaryHilbertSpace.FindStateIndex(TmpState, TmpKyMax), groundState[i]);
		}
            }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  if (nbrFermionSector == this->NbrFermions)
    {
      if ((kySector % this->KyMax) == this->TotalKy)
	{
	  unsigned long Mask  = ((0x1ul << (this->KyMax - subsytemSize)) - 1ul) << subsytemSize;
	  int TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) ==  0x0ul)
                ++TmpSize;
            }
	  if (TmpSize == 0)
	    {
	      RealMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
	    }
	  RealMatrix TmpEntanglementMatrix(TmpSize, 1, true);
	  TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
		{
 		  TmpEntanglementMatrix.AddToMatrixElement(TmpSize, 0, groundState[i]);
		  ++TmpSize;
		}
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
  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrFermions - nbrFermionSector) * subsytemSize);
  int ComplementaryMinTotalKy = ((this->NbrFermions - nbrFermionSector) * (this->NbrFermions - nbrFermionSector - 1)) / 2;
  while (ComplementaryTotalKy < ComplementaryMinTotalKy)
    ComplementaryTotalKy += this->KyMax;
  int ComplementaryMaxTotalKy = ((2 * (this->KyMax - subsytemSize - 1) + 1 - (this->NbrFermions - nbrFermionSector)) * (this->NbrFermions - nbrFermionSector)) / 2;    
  int NbrComplementarySpaces = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMax;
    }
  if (NbrComplementarySpaces == 0)
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;	  
    }
  FermionOnDisk** ComplementaryHilbertSpaces = new FermionOnDisk*[NbrComplementarySpaces];
  NbrComplementarySpaces = 0;
  ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrFermions - nbrFermionSector) * subsytemSize);
  while (ComplementaryTotalKy < ComplementaryMinTotalKy)
    ComplementaryTotalKy += this->KyMax;
  int TotalComplementaryHilbertSpaceDimension = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ComplementaryHilbertSpaces[NbrComplementarySpaces] = new FermionOnDisk(this->NbrFermions - nbrFermionSector, 
									     ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[NbrComplementarySpaces]->GetHilbertSpaceDimension();
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMax;
    }
 
  FermionOnDisk TmpDestinationHilbertSpace(nbrFermionSector, kySector, subsytemSize - 1);
  FermionOnSphereFull TmpFullComplementaryHilbertSpace(this->NbrFermions - nbrFermionSector, this->KyMax - subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,
				   TmpFullComplementaryHilbertSpace.HilbertSpaceDimension, true);
  TotalComplementaryHilbertSpaceDimension = 0;
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      for (int MinIndex = 0; MinIndex < ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpComplementaryState = ComplementaryHilbertSpaces[SpaceIndex]->StateDescription[MinIndex] << subsytemSize;
	  int ComplementaryIndex = TmpFullComplementaryHilbertSpace.FindStateIndex(ComplementaryHilbertSpaces[SpaceIndex]->StateDescription[MinIndex], 
										   ComplementaryHilbertSpaces[SpaceIndex]->StateLzMax[MinIndex]);
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax;
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
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem

ComplexMatrix FermionOnTorus::EvaluatePartialEntanglementMatrixFullKyPartB (int subsytemSize, int nbrFermionSector, int kySector, ComplexVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrFermionSector == 0))
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
      if ((kySector == this->TotalKy) && (nbrFermionSector == this->NbrFermions))
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
  if (nbrFermionSector == 0)
    {
      if (kySector == 0)
	{
	  unsigned long Mask  = (0x1ul << subsytemSize) - 1ul;
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
	  FermionOnSphereFull TmpFullComplementaryHilbertSpace(this->NbrFermions - nbrFermionSector, this->KyMax - subsytemSize - 1);
	  ComplexMatrix TmpEntanglementMatrix(1, TmpFullComplementaryHilbertSpace.HilbertSpaceDimension, true);
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
		{
		  unsigned long TmpState = this->StateDescription[i] >> subsytemSize;
		  int TmpKyMax = this->KyMax - subsytemSize - 1;
		  while ((TmpState >> TmpKyMax) == 0x0ul)
		    --TmpKyMax;
 		  TmpEntanglementMatrix.AddToMatrixElement(0, TmpFullComplementaryHilbertSpace.FindStateIndex(TmpState, TmpKyMax), groundState[i]);
		}
            }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  if (nbrFermionSector == this->NbrFermions)
    {
      if ((kySector % this->KyMax) == this->TotalKy)
	{
	  unsigned long Mask  = ((0x1ul << (this->KyMax - subsytemSize)) - 1ul) << subsytemSize;
	  int TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) ==  0x0ul)
                ++TmpSize;
            }
	  if (TmpSize == 0)
	    {
	      ComplexMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
	    }
	  ComplexMatrix TmpEntanglementMatrix(TmpSize, 1, true);
	  TmpSize = 0;
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              if ((this->StateDescription[i] & Mask) == 0x0ul)
		{
 		  TmpEntanglementMatrix.AddToMatrixElement(TmpSize, 0, groundState[i]);
		  ++TmpSize;
		}
            }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }

  cout << "toto" << endl;
  
  int MinIndex = 0;
  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrFermions - nbrFermionSector) * subsytemSize);
  int ComplementaryMinTotalKy = ((this->NbrFermions - nbrFermionSector) * (this->NbrFermions - nbrFermionSector - 1)) / 2;
  while (ComplementaryTotalKy < ComplementaryMinTotalKy)
    ComplementaryTotalKy += this->KyMax;
  int ComplementaryMaxTotalKy = ((2 * (this->KyMax - subsytemSize - 1) + 1 - (this->NbrFermions - nbrFermionSector)) * (this->NbrFermions - nbrFermionSector)) / 2;    
  int NbrComplementarySpaces = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMax;
    }
  if (NbrComplementarySpaces == 0)
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;	  
    }
  FermionOnDisk** ComplementaryHilbertSpaces = new FermionOnDisk*[NbrComplementarySpaces];
  NbrComplementarySpaces = 0;
  ComplementaryTotalKy = this->TotalKy - kySector - ((this->NbrFermions - nbrFermionSector) * subsytemSize);
  while (ComplementaryTotalKy < ComplementaryMinTotalKy)
    ComplementaryTotalKy += this->KyMax;
  int TotalComplementaryHilbertSpaceDimension = 0;
  while (ComplementaryTotalKy <= ComplementaryMaxTotalKy)
    {
      ComplementaryHilbertSpaces[NbrComplementarySpaces] = new FermionOnDisk(this->NbrFermions - nbrFermionSector, 
									     ComplementaryTotalKy, this->KyMax - subsytemSize - 1);
      TotalComplementaryHilbertSpaceDimension += ComplementaryHilbertSpaces[NbrComplementarySpaces]->GetHilbertSpaceDimension();
      ++NbrComplementarySpaces;
      ComplementaryTotalKy += this->KyMax;
    }
 
  FermionOnDisk TmpDestinationHilbertSpace(nbrFermionSector, kySector, subsytemSize - 1);
  FermionOnSphereFull TmpFullComplementaryHilbertSpace(this->NbrFermions - nbrFermionSector, this->KyMax - subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,
				   TmpFullComplementaryHilbertSpace.HilbertSpaceDimension, true);
  TotalComplementaryHilbertSpaceDimension = 0;
  for (int SpaceIndex = 0; SpaceIndex < NbrComplementarySpaces; ++SpaceIndex)
    {
      for (int MinIndex = 0; MinIndex < ComplementaryHilbertSpaces[SpaceIndex]->HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpComplementaryState = ComplementaryHilbertSpaces[SpaceIndex]->StateDescription[MinIndex] << subsytemSize;
	  int ComplementaryIndex = TmpFullComplementaryHilbertSpace.FindStateIndex(ComplementaryHilbertSpaces[SpaceIndex]->StateDescription[MinIndex], 
										   ComplementaryHilbertSpaces[SpaceIndex]->StateLzMax[MinIndex]);
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[j] | TmpComplementaryState;
	      int TmpKyMax = this->KyMax;
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
      cout << "tata" << endl;
      return TmpEntanglementMatrix;
    }
  else
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Ky sector.
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// kySector = Ky sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix FermionOnTorus::EvaluatePartialDensityMatrixParticlePartition (int nbrFermionSector, int kySector, RealVector& groundState)
{  
  if (nbrFermionSector == 0)
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

  if (nbrFermionSector == this->NbrFermions)
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

  int ComplementaryNbrFermionSector = this->NbrFermions - nbrFermionSector;
  int ComplementaryKySector = this->TotalKy - kySector;
  if (ComplementaryKySector < 0)
    ComplementaryKySector += this->KyMax;
  if (ComplementaryKySector >= this->KyMax)
    ComplementaryKySector -= this->KyMax;


  FermionOnTorus SubsytemSpace (nbrFermionSector, this->KyMax, kySector);
  RealSymmetricMatrix TmpDensityMatrix(SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnTorus ComplementarySpace(ComplementaryNbrFermionSector, this->KyMax, ComplementaryKySector);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.GetHilbertSpaceDimension() << endl;

   
  long TmpNbrNonZeroElements = 0;

  TmpNbrNonZeroElements = this->EvaluatePartialDensityMatrixParticlePartitionCore(0, ComplementarySpace.GetHilbertSpaceDimension(), &ComplementarySpace, &SubsytemSpace, groundState, &TmpDensityMatrix);
  if (TmpNbrNonZeroElements > 0)
    {
      return TmpDensityMatrix;
    }
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnTorus::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorus* complementaryHilbertSpace,  ParticleOnTorus* destinationHilbertSpace,
									RealVector& groundState,  RealSymmetricMatrix* densityMatrix)
{
  FermionOnTorus* TmpHilbertSpace =  (FermionOnTorus*) complementaryHilbertSpace;
  FermionOnTorus* TmpDestinationHilbertSpace =  (FermionOnTorus*) destinationHilbertSpace;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  BinomialCoefficients TmpBinomial (this->NbrFermions);
  double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace->StateDescription[minIndex];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace->StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
 	      int TmpKyMax = this->KyMax;
	      unsigned long TmpState3 = TmpState | TmpState2;
	      while ((TmpState3 >> TmpKyMax) == 0x0ul)
		--TmpKyMax;
	      int TmpPos = this->FindStateIndex(TmpState3, TmpKyMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  double Coefficient = TmpInvBinomial;
		  unsigned long Sign = 0x0ul;
		  int Pos2 = TmpDestinationHilbertSpace->KyMax;
		  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
		    {
		      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
			--Pos2;
		      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
		      TmpState3 ^= TmpState3 >> 32;
#endif	
		      TmpState3 ^= TmpState3 >> 16;
		      TmpState3 ^= TmpState3 >> 8;
		      TmpState3 ^= TmpState3 >> 4;
		      TmpState3 ^= TmpState3 >> 2;
		      TmpState3 ^= TmpState3 >> 1;
		      Sign ^= TmpState3;
		      TmpState2 &= ~(0x1ul << Pos2);
		      --Pos2;
		    }
 		  if ((Sign & 0x1ul) == 0x0ul)		  
 		    Coefficient *= 1.0;
 		  else
 		    Coefficient *= -1.0;
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  TmpStateCoefficient[Pos] = Coefficient;
		  ++Pos;
		}
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
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
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

ComplexVector& FermionOnTorus::CoreC4Rotation (ComplexVector& inputState, ParticleOnTorus* inputSpace, ComplexVector& outputState, int minIndex, int nbrIndices, 
					       bool clockwise)
{
  FermionOnTorus* TmpInputSpace = (FermionOnTorus*) inputSpace;
  unsigned long* TmpInputMonomial = new unsigned long [this->NbrFermions];
  unsigned long* TmpInputMonomial2 = new unsigned long [this->NbrFermions];
  unsigned long* TmpOutputMonomial = new unsigned long [this->NbrFermions];
  int LastIndex = minIndex + nbrIndices;
  Complex Tmp = 0.0;
  Complex Tmp2 = 0.0;
  double TmpCoefficient = pow((double) this->KyMax, -0.5 * ((double) this->NbrFermions));
  double PhaseFactor = 2.0 * M_PI / ((double) this->KyMax);
  if (clockwise == true)
    PhaseFactor *= -1.0;
  ComplexMatrix DeterminantMatrix (this->NbrFermions, this->NbrFermions);
  ComplexMatrix PhaseMatrix (this->KyMax, this->KyMax);
  for (int k = 0; k < this->KyMax; ++k)
    for (int l = 0; l < this->KyMax; ++l)
      PhaseMatrix[k][l] = Phase(PhaseFactor * ((double) (k * l)));
  for (int i = minIndex ; i < LastIndex; ++i)
    {
      this->ConvertToMonomial(this->StateDescription[i], TmpOutputMonomial);
      Tmp = 0.0;
      for (int j = 0; j < TmpInputSpace->HilbertSpaceDimension; ++j)
	{
	  TmpInputSpace->ConvertToMonomial(TmpInputSpace->StateDescription[j], TmpInputMonomial);
	  unsigned long TmpPhase = 0ul;
	  for (int k = 0; k < this->NbrFermions; ++k)
	    for (int l = 0; l < this->NbrFermions; ++l)
	      DeterminantMatrix[k][l] = PhaseMatrix[(int) TmpInputMonomial[k]][(int) TmpOutputMonomial[l]];
#ifdef __LAPACK__
	  Tmp += inputState[j] * TmpCoefficient * DeterminantMatrix.LapackDeterminant();	  
#else
	  Tmp += inputState[j] * TmpCoefficient * DeterminantMatrix.Determinant();
#endif
 	}
      outputState[i] = Tmp;      
    }
  delete[] TmpInputMonomial;
  delete[] TmpInputMonomial2;
  delete[] TmpOutputMonomial;
  return outputState;
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

void FermionOnTorus::SymmetrizeU1U1StateCore (ComplexVector& symmetrizedVector, ComplexVector& leftVector, ComplexVector& rightVector, ParticleOnTorus* leftSpace, ParticleOnTorus* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents)
{
  unsigned long LastComponent = firstComponent + nbrComponents;
  FermionOnTorus* TmpLeftSpace = (FermionOnTorus*) leftSpace;
  FermionOnTorus* TmpRightSpace = (FermionOnTorus*) rightSpace;
  for (long i = (long) firstComponent; i < (long) LastComponent; ++i)
    {
      unsigned long TmpLeftState = TmpLeftSpace->StateDescription[i];
      Complex TmpCoefficient = leftVector[i];
      
      for (long j = 0l; j < TmpRightSpace->LargeHilbertSpaceDimension; ++j)
	{
	  unsigned long TmpRightState = TmpRightSpace->StateDescription[j];
	  if ((TmpLeftState & TmpRightState) == 0x0ul)
	    {
	      int TmpKyMax = this->KyMax;
	      unsigned long TmpState = TmpLeftState | TmpRightState;
	      while ((TmpState >> TmpKyMax) == 0x0ul)
		--TmpKyMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpKyMax);
	      if (TmpPos < this->HilbertSpaceDimension)
		{
		  unsigned long Sign = 0x0ul;
		  int Pos = TmpRightSpace->StateKyMax[j];
		  while ((Pos > 0) && (TmpRightState != 0x0ul))
		    {
		      while (((TmpRightState >> Pos) & 0x1ul) == 0x0ul)
			--Pos;
		      TmpState = TmpLeftState & ((0x1ul << (Pos + 1)) - 1ul);
#ifdef  __64_BITS__
		      TmpState ^= TmpState >> 32;
#endif	
		      TmpState ^= TmpState >> 16;
		      TmpState ^= TmpState >> 8;
		      TmpState ^= TmpState >> 4;
		      TmpState ^= TmpState >> 2;
		      TmpState ^= TmpState >> 1;
		      Sign ^= TmpState;
		      TmpRightState &= ~(0x1ul << Pos);
		      --Pos;
		    }
 		  if ((Sign & 0x1ul) == 0x0ul)		  
		    symmetrizedVector[TmpPos] += TmpCoefficient * rightVector[j];
 		  else
		    symmetrizedVector[TmpPos] -= TmpCoefficient * rightVector[j];
		}
	    }
	}  
    }
}


// symmetrized a product of several uncoupled states 
//
// outputState = reference on the output state
// inputStates = states which will be symmetrized
// inputSpaces = Hilbert spaces attached to each states
// nbrStates = number of states to symmetrize
// firstComponent = first component to symmetrize within the first Hilbert space of inputSpaces
// nbrComponents = number of components to symmetrize within the first Hilbert space of inputSpaces

void FermionOnTorus::SymmetrizeU1U1StateCore (ComplexVector& outputState, ComplexVector* inputStates, ParticleOnTorus** inputSpaces, int nbrStates, unsigned long firstComponent, unsigned long nbrComponents)
{
  unsigned long LastComponent = firstComponent + nbrComponents;
  
  FermionOnTorus* TmpSpace1 = (FermionOnTorus*) inputSpaces[0];
  for (long i = (long) firstComponent; i < (long) LastComponent; ++i)
    {
      Complex TmpCoefficient = inputStates[0][i];
      unsigned long TmpState1 = TmpSpace1->StateDescription[i];
      FermionOnTorus* TmpSpace2 = (FermionOnTorus*) inputSpaces[1];            
      for (long j2 = 0l; j2 < TmpSpace2->LargeHilbertSpaceDimension; ++j2)
	{
	  unsigned long TmpState2 = TmpSpace2->StateDescription[j2];
	  if ((TmpState2 & TmpState1) == 0x0ul)
	    {
	      unsigned long Sign2 = 0x0ul;
	      int Pos = TmpSpace2->StateKyMax[j2];
	      unsigned long TmpState;
	      while ((Pos > 0) && (TmpState2 != 0x0ul))
		{
		  while (((TmpState2 >> Pos) & 0x1ul) == 0x0ul)
		    --Pos;
		  TmpState = TmpState1 & ((0x1ul << (Pos + 1)) - 1ul);
#ifdef  __64_BITS__
		  TmpState ^= TmpState >> 32;
#endif	
		  TmpState ^= TmpState >> 16;
		  TmpState ^= TmpState >> 8;
		  TmpState ^= TmpState >> 4;
		  TmpState ^= TmpState >> 2;
		  TmpState ^= TmpState >> 1;
		  Sign2 ^= TmpState;
		  TmpState2 &= ~(0x1ul << Pos);
		  --Pos;
		}
	      TmpState2 = TmpState1 | TmpSpace2->StateDescription[j2];
	      if (nbrStates == 2)
		{	  
		  int TmpKyMax = this->KyMax;
		  while ((TmpState2 >> TmpKyMax) == 0x0ul)
		    --TmpKyMax;
		  int TmpPos = this->FindStateIndex(TmpState2, TmpKyMax);
		  if (TmpPos < this->HilbertSpaceDimension)
		    {
		      if ((Sign2 & 0x1ul) == 0x0ul)		  
			outputState[TmpPos] += TmpCoefficient * inputStates[1][j2];
		      else
			outputState[TmpPos] -= TmpCoefficient * inputStates[1][j2];
		    }
		}
	      else
		{
		  Complex TmpCoefficient2 = TmpCoefficient * inputStates[1][j2];
		  FermionOnTorus* TmpSpace3 = (FermionOnTorus*) inputSpaces[2];            
		  for (long j3 = 0l; j3 < TmpSpace3->LargeHilbertSpaceDimension; ++j3)
		    {
		      unsigned long TmpState3 = TmpSpace3->StateDescription[j3];
		      if ((TmpState3 & TmpState2) == 0x0ul)
			{
			  unsigned long Sign3 = Sign2;
			  int Pos = TmpSpace3->StateKyMax[j3];
			  unsigned long TmpState;
			  while ((Pos > 0) && (TmpState3 != 0x0ul))
			    {
			      while (((TmpState3 >> Pos) & 0x1ul) == 0x0ul)
				--Pos;
			      TmpState = TmpState2 & ((0x1ul << (Pos + 1)) - 1ul);
#ifdef  __64_BITS__
			      TmpState ^= TmpState >> 32;
#endif	
			      TmpState ^= TmpState >> 16;
			      TmpState ^= TmpState >> 8;
			      TmpState ^= TmpState >> 4;
			      TmpState ^= TmpState >> 2;
			      TmpState ^= TmpState >> 1;
			      Sign3 ^= TmpState;
			      TmpState3 &= ~(0x1ul << Pos);
			      --Pos;
			    }
			  TmpState3 = TmpState2 | TmpSpace3->StateDescription[j3];
			  if (nbrStates == 3)
			    {	  
			      int TmpKyMax = this->KyMax;
			      while ((TmpState3 >> TmpKyMax) == 0x0ul)
				--TmpKyMax;
			      int TmpPos = this->FindStateIndex(TmpState3, TmpKyMax);
			      if (TmpPos < this->HilbertSpaceDimension)
				{
				  if ((Sign3 & 0x1ul) == 0x0ul)		  
				    outputState[TmpPos] += TmpCoefficient2 * inputStates[2][j3];
				  else
				    outputState[TmpPos] -= TmpCoefficient2 * inputStates[2][j3];
				}
			    }
			  else
			    {
			      Complex TmpCoefficient3 = TmpCoefficient2 * inputStates[2][j3];
			      FermionOnTorus* TmpSpace4 = (FermionOnTorus*) inputSpaces[3];            
			      for (long j4 = 0l; j4 < TmpSpace4->LargeHilbertSpaceDimension; ++j4)
				{
				  unsigned long TmpState4 = TmpSpace4->StateDescription[j4];
				  if ((TmpState4 & TmpState3) == 0x0ul)
				    {
				      unsigned long Sign4 = Sign3;
				      int Pos = TmpSpace4->StateKyMax[j4];
				      unsigned long TmpState;
				      while ((Pos > 0) && (TmpState4 != 0x0ul))
					{
					  while (((TmpState4 >> Pos) & 0x1ul) == 0x0ul)
					    --Pos;
					  TmpState = TmpState3 & ((0x1ul << (Pos + 1)) - 1ul);
#ifdef  __64_BITS__
					  TmpState ^= TmpState >> 32;
#endif	
					  TmpState ^= TmpState >> 16;
					  TmpState ^= TmpState >> 8;
					  TmpState ^= TmpState >> 4;
					  TmpState ^= TmpState >> 2;
					  TmpState ^= TmpState >> 1;
					  Sign4 ^= TmpState;
					  TmpState4 &= ~(0x1ul << Pos);
					  --Pos;
					}
				      TmpState4 = TmpState3 | TmpSpace4->StateDescription[j4];
				      if (nbrStates == 4)
					{	  
					  int TmpKyMax = this->KyMax;
					  while ((TmpState4 >> TmpKyMax) == 0x0ul)
					    --TmpKyMax;
					  int TmpPos = this->FindStateIndex(TmpState4, TmpKyMax);
					  if (TmpPos < this->HilbertSpaceDimension)
					    {
					      if ((Sign4 & 0x1ul) == 0x0ul)		  
						outputState[TmpPos] += TmpCoefficient3 * inputStates[3][j4];
					  else
					    outputState[TmpPos] -= TmpCoefficient2 * inputStates[2][j4];
					    }
					}
				      else
					{
					  cout << "error, current code cannot symmetrized more than 4 states" << endl;
					}
				    }
				}
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

void FermionOnTorus::CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, SparseRealMatrix& twistMatrix, RealVector& state, 
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
	  unsigned long TmpStateDescription = this->StateDescription[i];
	  TmpMatrix.Copy(twistMatrix);
	  for (int j = 0; j < this->KyMax; ++j)
	    {
	      TmpMatrix.Multiply(bMatrices[TmpStateDescription & 0x1ul], TmpMatrixElements, TmpColumnIndices, TmpElements);
	      TmpStateDescription >>= 1;
	    } 
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
      int PrecalculationBlockSize = (int) memory;
      unsigned long PrecalculationBlockMask = (0x1ul  << PrecalculationBlockSize) - 0x1ul;
      int PrecalculationBlockLength  = this->KyMax - (this->KyMax % PrecalculationBlockSize);
      int RemaingOrbtals = this->KyMax % PrecalculationBlockSize;
      SparseRealMatrix* TmpBlockbMatrices = new SparseRealMatrix[1 << PrecalculationBlockSize];
      for (unsigned long TmpState = 0x0ul; TmpState <= PrecalculationBlockMask; ++TmpState)
	{
	  unsigned long TmpStateDescription = TmpState;
	  TmpBlockbMatrices[TmpState].Copy(bMatrices[TmpStateDescription & 0x1ul]);
	  TmpStateDescription >>= 1;
	  for (int j = 1; j < PrecalculationBlockSize; ++j)
	    {
	      TmpBlockbMatrices[TmpState].Multiply(bMatrices[TmpStateDescription & 0x1ul], TmpMatrixElements, TmpColumnIndices, TmpElements);
	      TmpStateDescription >>= 1;
	    }      
	}
      for (long i = initialIndex; i < MaxIndex; ++i)
	{
	  if (((i - initialIndex) % 10000) == 0)
	    cout << "Completed " << (i - initialIndex) << " out of " << (MaxIndex - initialIndex) << endl; 
	  unsigned long TmpStateDescription = this->StateDescription[i];
	  TmpMatrix.Copy(twistMatrix);
	  int j = 0;
	  for (; j < PrecalculationBlockLength; j += PrecalculationBlockSize)
	    {
	      TmpMatrix.Multiply(TmpBlockbMatrices[TmpStateDescription & PrecalculationBlockMask], TmpMatrixElements, TmpColumnIndices, TmpElements);
	      TmpStateDescription >>= PrecalculationBlockSize;
	    } 
	  for (; j < this->KyMax; ++j)
	    {
	      TmpMatrix.Multiply(bMatrices[TmpStateDescription & 0x1ul], TmpMatrixElements, TmpColumnIndices, TmpElements);
	      TmpStateDescription >>= 1;
	    } 
	  double& TmpComponent = state[i];
	  double Tmp;
	  for (int j = 0; j < nbrMPSSumIndices; ++j)
	    {
	      TmpMatrix.GetMatrixElement(mPSSumIndices[j], mPSSumIndices[j], Tmp);
	      TmpComponent += Tmp;
	    }
	}
    }
  delete[] TmpMatrixElements;
  delete[] TmpColumnIndices;
  delete[] TmpElements;
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

void FermionOnTorus::CreateStateFromMPSDescription (SparseComplexMatrix* bMatrices, SparseRealMatrix& twistMatrix, ComplexVector& state, 
						    int* mPSSumIndices, int nbrMPSSumIndices,
						    long memory, long initialIndex, long nbrComponents)
{
  SparseComplexMatrix TmpMatrix;
  long MaxIndex = initialIndex + nbrComponents;
  if ((nbrComponents == 0l) || (MaxIndex > this->LargeHilbertSpaceDimension))
    {
      MaxIndex = this->LargeHilbertSpaceDimension;
    }
  long TmpMemory = (((long) bMatrices[0].GetNbrRow()) * ((long)bMatrices[0].GetNbrColumn()));
  if (TmpMemory > (1l << 28))
    TmpMemory = 1l << 28;
  Complex* TmpMatrixElements = new Complex [TmpMemory];
  int* TmpColumnIndices = new int [TmpMemory];
  Complex* TmpElements = new Complex [bMatrices[0].GetNbrRow()];

  if (memory <= 1l)
    {
      for (long i = initialIndex; i < MaxIndex; ++i)
	{
	  if (((i - initialIndex) % 10000) == 0)
	    cout << "Completed " << (i - initialIndex) << " out of " << (MaxIndex - initialIndex) << endl; 
	  unsigned long TmpStateDescription = this->StateDescription[i];
	  TmpMatrix.Copy(twistMatrix);
	  for (int j = 0; j < this->KyMax; ++j)
	    {
	      TmpMatrix.Multiply(bMatrices[TmpStateDescription & 0x1ul], TmpMatrixElements, TmpColumnIndices, TmpElements);
	      TmpStateDescription >>= 1;
	    } 
	  Complex& TmpComponent = state[i];
	  Complex Tmp;
	  for (int j = 0; j < nbrMPSSumIndices; ++j)
	    {
	      TmpMatrix.GetMatrixElement(mPSSumIndices[j], mPSSumIndices[j], Tmp);
	      TmpComponent += Tmp;
	    }
	}
    }
  else
    {
      int PrecalculationBlockSize = (int) memory;
      unsigned long PrecalculationBlockMask = (0x1ul  << PrecalculationBlockSize) - 0x1ul;
      int PrecalculationBlockLength  = this->KyMax - (this->KyMax % PrecalculationBlockSize);
      int RemaingOrbtals = this->KyMax % PrecalculationBlockSize;
      SparseComplexMatrix* TmpBlockbMatrices = new SparseComplexMatrix[1 << PrecalculationBlockSize];
      for (unsigned long TmpState = 0x0ul; TmpState <= PrecalculationBlockMask; ++TmpState)
	{
	  unsigned long TmpStateDescription = TmpState;
	  TmpBlockbMatrices[TmpState].Copy(bMatrices[TmpStateDescription & 0x1ul]);
	  TmpStateDescription >>= 1;
	  for (int j = 1; j < PrecalculationBlockSize; ++j)
	    {
	      TmpBlockbMatrices[TmpState].Multiply(bMatrices[TmpStateDescription & 0x1ul], TmpMatrixElements, TmpColumnIndices, TmpElements);
	      TmpStateDescription >>= 1;
	    }      
	}
      for (long i = initialIndex; i < MaxIndex; ++i)
	{
	  if (((i - initialIndex) % 10000) == 0)
	    cout << "Completed " << (i - initialIndex) << " out of " << (MaxIndex - initialIndex) << endl; 
	  unsigned long TmpStateDescription = this->StateDescription[i];
	  TmpMatrix.Copy(twistMatrix);
	  int j = 0;
	  for (; j < PrecalculationBlockLength; j += PrecalculationBlockSize)
	    {
	      TmpMatrix.Multiply(TmpBlockbMatrices[TmpStateDescription & PrecalculationBlockMask], TmpMatrixElements, TmpColumnIndices, TmpElements);
	      TmpStateDescription >>= PrecalculationBlockSize;
	    } 
	  for (; j < this->KyMax; ++j)
	    {
	      TmpMatrix.Multiply(bMatrices[TmpStateDescription & 0x1ul], TmpMatrixElements, TmpColumnIndices, TmpElements);
	      TmpStateDescription >>= 1;
	    } 
	  Complex& TmpComponent = state[i];
	  Complex Tmp;
	  for (int j = 0; j < nbrMPSSumIndices; ++j)
	    {
	      TmpMatrix.GetMatrixElement(mPSSumIndices[j], mPSSumIndices[j], Tmp);
	      TmpComponent += Tmp;
	    }
	}
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

bool FermionOnTorus::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  unsigned long TmpState = this->StateDescription[index];
  unsigned long RMask = (0x1ul << pauliR) - 0x1ul;
  TmpState |= (TmpState & RMask) << this->KyMax;
  unsigned long UnsignedPauliK = (unsigned long) pauliK;
  for (int i = 0; i < this->KyMax; ++i)
    {
      unsigned long TmpOccupation = 0l;
      unsigned long TmpState2 = TmpState & RMask;
      while (TmpState2 != 0x0ul)
	{
	  TmpOccupation += (TmpState2 & 0x1ul);
	  TmpState2 >>= 1;
	}
      if (TmpOccupation > UnsignedPauliK)
	return false;
      TmpState >>= 1;
    }
  return true;
}



 
// symmetrize a vector by grouping distant and equally separated orbitals, core part
//
// inputVector = reference on the vector to symmetrize
// nbrOrbitals = number of orbitals to group together
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Ky to the largest Ky
// first component = index of the first vector component 
// last component = index of the last component

void FermionOnTorus::SymmetrizeSingleStateGroupingDistantOrbitalsCore (ComplexVector& inputVector, ComplexVector* symmetrizedVectors, int nbrOrbitals, unsigned long firstComponent, unsigned long nbrComponents, bool twistedTorus)
{
  long LastComponent = (long) (firstComponent + nbrComponents);
  int TargetSpaceNbrOrbitals = this->KyMax / nbrOrbitals;
  FermionOnTorus** TargetSpaces = new FermionOnTorus* [TargetSpaceNbrOrbitals];
  unsigned long* TmpState = new unsigned long[TargetSpaceNbrOrbitals];
  
  for (int i = 0; i < TargetSpaceNbrOrbitals; ++i)
    {
      TargetSpaces[i] = 0;
    }  
    
  FermionOnTorus** IntermediaryTargetSpace = new FermionOnTorus* [nbrOrbitals];
  for (int i = 0; i < nbrOrbitals; ++i)
    IntermediaryTargetSpace[i] = 0;
  IntermediaryTargetSpace[0] = new FermionOnTorus(this->NbrFermions, this->KyMax, this->TotalKy);

  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      unsigned long TmpNbrParticles = 0x0ul;
      unsigned long TmpState = this->StateDescription[i];
      unsigned long TmpState1 = TmpState;
      unsigned long TmpState2 = 0x0ul;
      int TmpTotalKy = 0;
      bool OrbitalOverflow = false;
      double sign = 1;
      int TmpIndex = 0;
      int TmpStateIndex = i;
      int TmpIndexTargetSpace = 0;
      int TmpIndexNewTargetSpace = 0;
      for (int k = 0; (k < TargetSpaceNbrOrbitals) && (OrbitalOverflow == false); ++k)
      {
	TmpState = (this->StateDescription[i]) >> k ;
	TmpNbrParticles = 0x0ul;
	for (int l = 0; l < nbrOrbitals; ++l)
	{
	  TmpNbrParticles += TmpState & 0x1ul;
	  if ((TmpState & 0x1ul) == 0x1ul)
	    TmpIndex = l;
	  TmpState >>= TargetSpaceNbrOrbitals;
	}
	if (TmpNbrParticles == 0x1ul)
	{
	  TmpState2 |= 0x1ul << k;
	  TmpTotalKy += k;
	  TmpIndexNewTargetSpace += TmpIndex;
	  TmpIndexNewTargetSpace %= nbrOrbitals;
	  if (IntermediaryTargetSpace[TmpIndexNewTargetSpace] == 0)
	    IntermediaryTargetSpace[TmpIndexNewTargetSpace] = new FermionOnTorus (this->NbrFermions, this->KyMax, ((this->TotalKy - TmpIndexNewTargetSpace*TargetSpaceNbrOrbitals + this->KyMax) % this->KyMax));
	  IntermediaryTargetSpace[TmpIndexTargetSpace]->SetTargetSpace(IntermediaryTargetSpace[TmpIndexNewTargetSpace]);
	  double coefficient;
	  TmpStateIndex = IntermediaryTargetSpace[TmpIndexTargetSpace]->AdA (TmpStateIndex, k, k + TargetSpaceNbrOrbitals*TmpIndex, coefficient);
	  sign *= coefficient;
	  TmpIndexTargetSpace = TmpIndexNewTargetSpace;
	}
	else
	{
	if (TmpNbrParticles > 0x1ul)
	  OrbitalOverflow = true;
	}
      }

      TmpTotalKy %= TargetSpaceNbrOrbitals;
      if (OrbitalOverflow == false)
	{
	  if (TargetSpaces[TmpTotalKy] == 0)
	    {
	      TargetSpaces[TmpTotalKy] = new FermionOnTorus (this->NbrFermions, TargetSpaceNbrOrbitals, TmpTotalKy);
	      symmetrizedVectors[TmpTotalKy] = ComplexVector (TargetSpaces[TmpTotalKy]->HilbertSpaceDimension, true);
	    }	  
	  int TmpKyMax = TargetSpaces[TmpTotalKy]->KyMax - 1;
	  while ((TmpState2 >> TmpKyMax) == 0x0ul)
	    --TmpKyMax;
	  int TmpPos = TargetSpaces[TmpTotalKy]->FindStateIndex(TmpState2, TmpKyMax);
	  if (TmpPos < TargetSpaces[TmpTotalKy]->HilbertSpaceDimension)
	    {
	      symmetrizedVectors[TmpTotalKy][TmpPos] += sign * inputVector[i];
	    }
	}
    }
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

void FermionOnTorus::SymmetrizeSingleStatePeriodicSubsetOrbitalCore (ComplexVector& inputVector, ComplexVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, 
								     unsigned long firstComponent, unsigned long nbrComponents)
{
  bool twistedTorus = true;
  long LastComponent = (long) (firstComponent + nbrComponents);
  int TargetSpaceNbrOrbitals = this->KyMax / periodicity;
  FermionOnTorus** TargetSpaces = new FermionOnTorus* [TargetSpaceNbrOrbitals];
  for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
    {
      TargetSpaces[j] = 0;
    }

  unsigned long* TmpState = new unsigned long[TargetSpaceNbrOrbitals];
  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      Complex TmpCoefficient = inputVector[i];
      unsigned long TmpInputState = this->StateDescription[i];
      int TmpTotalKy = 0;
      int TmpNbrParticles = 0;
      unsigned long TmpOutputState = 0x0ul;
      int CurrentLz = 0;
      for (int k = firstOrbitalIndex; k <= this->KyMax; k += periodicity)
	{
	  if ((TmpInputState & (0x1ul << k)) != 0x0ul)
	    {
	      TmpOutputState |= 0x1ul << CurrentLz;
	      TmpTotalKy += CurrentLz;
	      ++TmpNbrParticles;
	    }
	  ++CurrentLz;
	} 
      TmpTotalKy %= TargetSpaceNbrOrbitals;
      if (TmpNbrParticles == this->NbrFermions)
	{
	  if (TargetSpaces[TmpTotalKy] == 0)
	    {
	      TargetSpaces[TmpTotalKy] = new FermionOnTorus (this->NbrFermions, TargetSpaceNbrOrbitals, TmpTotalKy);
	      symmetrizedVectors[this->NbrFermions][TmpTotalKy] = ComplexVector(TargetSpaces[TmpTotalKy]->HilbertSpaceDimension, true);
	    }	  
	  int TmpKyMax = TargetSpaces[TmpTotalKy]->KyMax - 1;
	  while ((TmpOutputState >> TmpKyMax) == 0x0ul)
	    --TmpKyMax;
	  int TmpPos = TargetSpaces[TmpTotalKy]->FindStateIndex(TmpOutputState, TmpKyMax);
	  if (TmpPos < TargetSpaces[TmpTotalKy]->HilbertSpaceDimension)
	    {
	      symmetrizedVectors[this->NbrFermions][TmpTotalKy][TmpPos] += TmpCoefficient;
	    }
      }
    }
  delete[] TmpState;
  for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
    {
      if (TargetSpaces[j] != 0)
	{
	  delete TargetSpaces[j];
	}
    }
  delete[] TargetSpaces;
}

// symmetrize a vector by keeping only a subset of equally separated orbitals
//
// inputVector = reference on the vector to symmetrize
// firstOrbitalIndex = index of the first orbital to keep
// symmetrizedVectors = array on the symmetrize states ranging from the smallest Ky to the largest Ky
// periodicity = momentum periodicity (should be a multiple of the number of orbitals)
// phase = an optional phase (in pi units) that can be added for each kept and occupied orbital
// firstComponent = first component of the input vector that has to be symmetrized
// nbrComponents = number of components of the input vector that have to be symmetrized
// return value = symmetrized state

void FermionOnTorus::SymmetrizeSingleStatePeriodicSubsetOrbitalCore (ComplexVector& inputVector, ComplexVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, double phase, 
								     unsigned long firstComponent, unsigned long nbrComponents)
{
  bool twistedTorus = true;
  long LastComponent = (long) (firstComponent + nbrComponents);
  int TargetSpaceNbrOrbitals = this->KyMax / periodicity;
  FermionOnTorus** TargetSpaces = new FermionOnTorus* [TargetSpaceNbrOrbitals];
  for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
    {
      TargetSpaces[j] = 0;
    }
  double TmpPhase = -2.0 * M_PI / ((double) this->KyMax);
  if (phase >= (0.5 * M_PI))
    TmpPhase *= -1.0;
  unsigned long* TmpState = new unsigned long[TargetSpaceNbrOrbitals];
  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      Complex TmpCoefficient = inputVector[i];
      unsigned long TmpInputState = this->StateDescription[i];
      int TmpTotalKy = 0;
      int TmpTotalSqrKy = 0;
      int TmpNbrParticles = 0;
      unsigned long TmpOutputState = 0x0ul;
      int CurrentLz = 0;
      for (int k = firstOrbitalIndex; k <= this->KyMax; k += periodicity)
	{
	  if ((TmpInputState & (0x1ul << k)) != 0x0ul)
	    {
	      TmpOutputState |= 0x1ul << CurrentLz;
	      TmpTotalKy += CurrentLz;
	      TmpTotalSqrKy += CurrentLz * CurrentLz;
	      ++TmpNbrParticles;	      
	    }
	  ++CurrentLz;
	} 
      int FullTmpTotalKy = TmpTotalKy;
      TmpTotalKy %= TargetSpaceNbrOrbitals;
      if (TmpNbrParticles == this->NbrFermions)
	{
	  if (TargetSpaces[TmpTotalKy] == 0)
	    {
	      TargetSpaces[TmpTotalKy] = new FermionOnTorus (this->NbrFermions, TargetSpaceNbrOrbitals, TmpTotalKy);
	      symmetrizedVectors[this->NbrFermions][TmpTotalKy] = ComplexVector(TargetSpaces[TmpTotalKy]->HilbertSpaceDimension, true);
	    }	  
	  int TmpKyMax = TargetSpaces[TmpTotalKy]->KyMax - 1;
	  while ((TmpOutputState >> TmpKyMax) == 0x0ul)
	    --TmpKyMax;
	  int TmpPos = TargetSpaces[TmpTotalKy]->FindStateIndex(TmpOutputState, TmpKyMax);
	  if (TmpPos < TargetSpaces[TmpTotalKy]->HilbertSpaceDimension)
	    {
	      symmetrizedVectors[this->NbrFermions][TmpTotalKy][TmpPos] += TmpCoefficient * Phase(TmpPhase * (((double) TmpTotalSqrKy) - ((double) (FullTmpTotalKy * TargetSpaceNbrOrbitals))));
	    }
      }
    }
  delete[] TmpState;
  for (int j = 0; j < TargetSpaceNbrOrbitals; ++j)
    {
      if (TargetSpaces[j] != 0)
	{
	  delete TargetSpaces[j];
	}
    }
  delete[] TargetSpaces;
}

// apply a magnetic translation along x to a given state
//
// index = state index 
// return value = translated state index, negative if a minus phase has to be added

int FermionOnTorus::ApplyXMagneticTranslation(int index)
{
  unsigned long TmpState = this->StateDescription[index];
  TmpState = (TmpState >> this->StateShift) | ((TmpState & this->MomentumMask) << this->ComplementaryStateShift);
  int TmpKyMax = this->KyMax;
  while ((TmpState >> TmpKyMax) == 0x0ul)
    --TmpKyMax;
  return this->FindStateIndex(TmpState, TmpKyMax);
}

// apply a magnetic translation along x to a given state
//
// index = state index 
// sign = additional sign due to the particle statistics
// return value = translated state index

int FermionOnTorus::ApplyXMagneticTranslation(int index, double& sign)
{
  unsigned long TmpState = this->StateDescription[index];
  if (this->NbrFermions & 1)
    {
      sign = 1.0;
      TmpState = (TmpState >> this->StateShift) | ((TmpState & this->MomentumMask) << this->ComplementaryStateShift);
    }
  else
    {
      unsigned long TmpState2 = TmpState & this->MomentumMask;
      TmpState = (TmpState >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
#ifdef  __64_BITS__
      TmpState2 ^= TmpState2 >> 32;
#endif
      TmpState2 ^= TmpState2 >> 16;
      TmpState2 ^= TmpState2 >> 8;
      TmpState2 ^= TmpState2 >> 4;
      TmpState2 ^= TmpState2 >> 2;
      TmpState2 ^= TmpState2 >> 1;
      if (TmpState2 & 0x1ul)
	sign = -1.0;
      else
	sign = -1.0;	
    }
  int TmpKyMax = this->KyMax;
  while ((TmpState >> TmpKyMax) == 0x0ul)
    --TmpKyMax;
  return this->FindStateIndex(TmpState, TmpKyMax);
}

// transform a state expressed on a torus with a given angle to a state expressed on the same trous but a different angle
//
// inputVector = reference on the input vector
// inputAngle = angle (in radian) between the two vectors that span the torus on which the input state is defined
// inputAspectRatio = length ratio of the two vectors that span the torus on which the input state is defined
// outputAngle = angle (in radian) between the two vectors that span the torus on which the output state is defined
// outputAspectRatio = length ratio of the two vectors that span the torus on which the output state is defined
// firstComponent = first component of the input vector that has to be symmetrized
// nbrComponents = number of components of the input vector that have to be symmetrized
// return value = transformed state

ComplexVector FermionOnTorus::ChangeTorusAngle (ComplexVector& inputVector, double inputAngle, double inputAspectRatio, double outputAngle, double outputAspectRatio,
						unsigned long firstComponent, unsigned long nbrComponents)
{
  long LastComponent = (long) (firstComponent + nbrComponents);
  ComplexVector OutputState (this->HilbertSpaceDimension);
  double TmpPhase = - M_PI / ((double) this->KyMax);
  if (inputAngle >= (0.5 * M_PI))
    TmpPhase *= -1.0;
  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      unsigned long TmpInputState = this->StateDescription[i];
      int TmpTotalSqrKy = 0;
      for (int k = 0; k < this->KyMax; ++k)
	{
	  if ((TmpInputState & (0x1ul << k)) != 0x0ul)
	    {
	      TmpTotalSqrKy += k * (k - this->KyMax);
	    }
	} 
      OutputState[i] = inputVector[i] * Phase(TmpPhase * ((double) TmpTotalSqrKy));
    }
  return OutputState;
}

