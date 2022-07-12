////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of fermions on a torus with spin                  //
//                                                                            //
//                        last modification : 10/09/2002                      //
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
#include "HilbertSpace/FermionOnTorusWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;


// constructor with a constraint on total spin momentum and total momentum
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// totalSpinMomentum = twice the total spin momentum to be used as constraint
// momentumConstraint = index of the momentum orbit

FermionOnTorusWithSpin::FermionOnTorusWithSpin (int nbrFermions, int maxMomentum, int totalSpinMomentum, int momentumConstaint)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = maxMomentum;
  this->TotalSpin = totalSpinMomentum;
  this->TotalLz = momentumConstaint;
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrLzValue - 1, 0, this->NbrFermionsUp);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions,  2 * this->NbrLzValue - 1, 2 * this->NbrLzValue - 1, 0, 0, 0);
  this->TargetSpace = this;
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

// constructor with a constraint on total momentum
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// momentumConstraint = index of the momentum orbit

FermionOnTorusWithSpin::FermionOnTorusWithSpin (int nbrFermions, int maxMomentum, int momentumConstaint)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = maxMomentum;
  this->TotalSpin = 0;
  this->TotalLz = momentumConstaint;
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrLzValue - 1, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions,  2 * this->NbrLzValue - 1, 2 * this->NbrLzValue - 1, 0, 0);
  if (((long) this->HilbertSpaceDimension) != this->LargeHilbertSpaceDimension)
    cout << "error, Hilbert space dimension mismatch " << this->HilbertSpaceDimension << " (" << this->LargeHilbertSpaceDimension << ")" << endl;
  this->TargetSpace = this;
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

FermionOnTorusWithSpin::FermionOnTorusWithSpin(const FermionOnTorusWithSpin& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->HighestBit = fermions.HighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnTorusWithSpin::~FermionOnTorusWithSpin ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnTorusWithSpin& FermionOnTorusWithSpin::operator = (const FermionOnTorusWithSpin& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnTorusWithSpin::Clone()
{
  return new FermionOnTorusWithSpin(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnTorusWithSpin::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  List<AbstractQuantumNumber*> TmpList2;
  TmpList2 += new SzQuantumNumber (this->TotalSpin);
  TmpList2 += new PeriodicMomentumQuantumNumber (this->TotalLz, this->NbrLzValue);
  TmpList += new VectorQuantumNumber (TmpList2);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnTorusWithSpin::GetQuantumNumber (int index)
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList +=  new SzQuantumNumber (this->TotalSpin);
  TmpList += new PeriodicMomentumQuantumNumber (this->TotalLz, this->NbrLzValue);
  return new VectorQuantumNumber (TmpList);
}

// get momemtum value of a given state
//
// index = state index
// return value = state momentum

int FermionOnTorusWithSpin::GetMomentumValue(int index)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned int State = this->StateDescription[index];
  int Momentum = 0;
  for (int i = 0; i <= StateHighestBit; ++i)
    {
      Momentum += (((State >> (i << 1)) & 1) + ((State >> ((i << 1) + 1)) & 1)) * i;
    }
  return (Momentum % this->NbrLzValue);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// currentTotalSpinMomentum = total spin momemtum of the fermions that have already been placed
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

int FermionOnTorusWithSpin::GenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos, 
					   int currentTotalSpinMomentum, int currentMomentum)
{
  if ((nbrFermions == 0) || (nbrFermions > (currentMaxMomentum + 1)) || (currentMaxMomentum < 0) || 
      (abs(this->TotalSpin - currentTotalSpinMomentum) > nbrFermions))
    return pos;
  if (nbrFermions == 1)
    {
      if (this->TotalSpin > currentTotalSpinMomentum)
	{
	  int i = this->TotalLz - (currentMomentum % this->NbrLzValue);
	  if (i < 0)
	    i += this->NbrLzValue;
	  i <<= 1;
	  currentMaxMomentum &= ~0x1;
	  for (; i <= currentMaxMomentum; i += (this->NbrLzValue << 1))
	    {
	      this->StateDescription[pos] = 0x1ul << i;
	      this->StateHighestBit[pos] = maxMomentum >> 1;
	      ++pos;
	    }
	}
      else
	{
	  int i = this->TotalLz - (currentMomentum % this->NbrLzValue);
	  if (i < 0)
	    i += this->NbrLzValue;
	  i <<= 1;
	  ++i;
	  if ((currentMaxMomentum & 1) == 0)
	    {
	      --currentMaxMomentum;
	    }
	  for (; i <= currentMaxMomentum; i += (this->NbrLzValue << 1))
	    {
	      this->StateDescription[pos] = 0x1l << i;
	      this->StateHighestBit[pos] = maxMomentum >> 1;
	      ++pos;
	    }
	}
      return pos;
    }
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int ReducedCurrentTotalSpinMomentum = currentTotalSpinMomentum;
  if ((currentMaxMomentum & 1) == 0)
    ReducedCurrentTotalSpinMomentum += 1;
  else
    ReducedCurrentTotalSpinMomentum -= 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentMaxMomentum, pos, ReducedCurrentTotalSpinMomentum, currentMomentum + (currentMaxMomentum >> 1));
  unsigned long Mask = 0x1l << currentMaxMomentum;
  for (int i = pos; i < TmpPos; ++i)
    this->StateDescription[i] |= Mask;
  if (maxMomentum == currentMaxMomentum)
    return this->GenerateStates(nbrFermions, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentTotalSpinMomentum, currentMomentum);
  else
    return this->GenerateStates(nbrFermions, maxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentTotalSpinMomentum, currentMomentum);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

int FermionOnTorusWithSpin::GenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos, 
					   int currentMomentum)
{
  if ((nbrFermions < 0) || (nbrFermions > (currentMaxMomentum + 1)) || ((nbrFermions > 0) && (currentMaxMomentum < 0)))
    return pos;
  if (nbrFermions == 0)
    {
      if ((currentMomentum % this->NbrLzValue) == this->TotalLz)
	{
	  this->StateDescription[pos] = 0x0ul;
	  this->StateHighestBit[pos] = maxMomentum >> 1;
	  ++pos;
	}
      return pos;
    }
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentMaxMomentum, pos, currentMomentum + (currentMaxMomentum >> 1));
  unsigned long Mask = 0x1ul << currentMaxMomentum;
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  if (maxMomentum == currentMaxMomentum)
    return this->GenerateStates(nbrFermions, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentMomentum);
  else
    return this->GenerateStates(nbrFermions, maxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentMomentum);
}

// evaluate Hilbert space dimension for a given total spin momentum
//
// nbrFermions = number of fermions
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// nbrSpinUp = number of particles with spin up
// return value = Hilbert space dimension

long FermionOnTorusWithSpin::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKy, int currentTotalKy, int nbrSpinUp)
{
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrFermions))
    return 0l;

  if (nbrFermions == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->TotalLz)
	return 1l;
      else	
	return 0l;
    }
  if (currentKy < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if (((j + currentTotalKy) % this->NbrLzValue) == this->TotalLz)
	    Count++;
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKy - 1, currentTotalKy + (2 * currentKy), nbrSpinUp - 1);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKy - 1, currentTotalKy + currentKy, nbrSpinUp - 1);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKy - 1, currentTotalKy + currentKy, nbrSpinUp);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKy - 1, currentTotalKy, nbrSpinUp);
  return Count;
}

// evaluate Hilbert space dimension for a given total momentum
//
// nbrFermions = number of fermions
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnTorusWithSpin::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKy, int currentTotalKy)
{
  if (nbrFermions == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->TotalLz)
	return 1l;
      else	
	return 0l;
    }
  if (currentKy < 0)
    return 0l;
  long  Count = 0l;
  if (nbrFermions > 1)
    Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKy - 1, currentTotalKy + (2 * currentKy));
  Count += 2l * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKy - 1, currentTotalKy + currentKy);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKy - 1, currentTotalKy);
  return Count;
}

