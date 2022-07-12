////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on torus with magnetic translations            //
//                                                                            //
//                        last modification : 13/10/2003                      //
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
#include "HilbertSpace/BosonOnTorusWithMagneticTranslations.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "MathTools/IntegerAlgebraTools.h"

#include <math.h>


using std::cout;
using std::endl;


// basic constructor
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// yMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)
// memory = amount of memory that can be used to allocate the look-up table (in bytes)

BosonOnTorusWithMagneticTranslations::BosonOnTorusWithMagneticTranslations (int nbrBosons, int maxMomentum, 
									    int xMomentum, int yMomentum, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = maxMomentum;

  this->NbrMomentum = this->MaxMomentum + 1;
  this->MomentumModulo = FindGCD(this->NbrBosons, this->MaxMomentum);
  this->XMomentum = xMomentum % this->MomentumModulo;
  this->YMomentum = yMomentum % this->MaxMomentum;
 
 this->XMomentumTranslationStep = this->MaxMomentum / this->MomentumModulo;

  this->ReducedNbrState = GetReducedNbrState (this->MaxMomentum);
  this->RemainderNbrState = GetRemainderNbrState (this->MaxMomentum);
  this->TemporaryState.Resize(this->ReducedNbrState);
  this->TemporaryState2.Resize(this->ReducedNbrState);

  this->HilbertSpaceDimension = (2 * this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->MaxMomentum)) / (this->MomentumModulo);
//  cout << this->HilbertSpaceDimension << endl;
  this->HilbertSpaceDimension = this->GenerateStates();
//  cout << this->HilbertSpaceDimension << endl;

  this->HashKeyMask = (unsigned long) 0x1;
  unsigned long TmpMemory = this->MaxMomentum * ((this->HashKeyMask + 1) * ((this->HilbertSpaceDimension * sizeof(int)) + 
									    sizeof(int*) + sizeof(int)) + sizeof(int**) + sizeof(int*));
  while (TmpMemory < memory)
    {
      this->HashKeyMask <<= 1;
      this->HashKeyMask |= (unsigned long) 0x1;
      unsigned long TmpMemory2 = this->MaxMomentum * ((this->HashKeyMask + 1) * ((this->HilbertSpaceDimension * sizeof(int)) + 
										 sizeof(int*) + sizeof(int)) + sizeof(int**) + sizeof(int*));
      if (TmpMemory2 < TmpMemory)
	TmpMemory = memory + 1;
      else
	TmpMemory = TmpMemory2;
    }
  this->HashKeyMask >>= 1;
  this->GenerateLookUpTable();
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

  this->Flag.Initialize();

#ifdef __DEBUG__
/*  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += (this->StateMaxMomentum[i] + 1) * sizeof(int) + sizeof(int*) + sizeof(int);
  cout << "memory requested for Hilbert space = ";
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

BosonOnTorusWithMagneticTranslations::BosonOnTorusWithMagneticTranslations(const BosonOnTorusWithMagneticTranslations& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;

  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->XMomentum = bosons.XMomentum;
  this->YMomentum = bosons.YMomentum;
  this->MomentumModulo = bosons.MomentumModulo;
  this->XMomentumTranslationStep = bosons.XMomentumTranslationStep;

  this->ReducedNbrState = bosons.ReducedNbrState;
  this->RemainderNbrState = bosons.RemainderNbrState;
  this->TemporaryState.Resize(this->ReducedNbrState);
  this->TemporaryState2.Resize(this->ReducedNbrState);

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;

  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  this->LookUpTable = bosons.LookUpTable;
  this->NbrStateInLookUpTable = bosons.NbrStateInLookUpTable;
  this->HashKeyMask = bosons.HashKeyMask;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

  this->Flag = bosons.Flag;
}

// destructor
//

BosonOnTorusWithMagneticTranslations::~BosonOnTorusWithMagneticTranslations ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;

      delete[] this->CompatibilityWithXMomentum;
      delete[] this->CompatibilityWithYMomentum;

      for (int i = 0; i < this->MaxMomentum; ++i)
	{
	  if (this->NbrStateInLookUpTable[i] != 0)
	    {
	      for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
		if (this->NbrStateInLookUpTable[i][j] > 0)
		  delete[] this->LookUpTable[i][j];
	      delete[] this->LookUpTable[i];
	      delete[] this->NbrStateInLookUpTable[i];
	    }
	}
      delete[] this->LookUpTable;
      delete[] this->NbrStateInLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorusWithMagneticTranslations& BosonOnTorusWithMagneticTranslations::operator = (const BosonOnTorusWithMagneticTranslations& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;

      delete[] this->CompatibilityWithXMomentum;
      delete[] this->CompatibilityWithYMomentum;

      for (int i = 0; i < this->MaxMomentum; ++i)
	{
	  if (this->NbrStateInLookUpTable[i] != 0)
	    {
	      for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
		if (this->NbrStateInLookUpTable[i][j] > 0)
		  delete[] this->LookUpTable[i][j];
	      delete[] this->LookUpTable[i];
	      delete[] this->NbrStateInLookUpTable[i];
	    }
	}
      delete[] this->LookUpTable;
      delete[] this->NbrStateInLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }

  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;

  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->XMomentum = bosons.XMomentum;
  this->YMomentum = bosons.YMomentum;
  this->MomentumModulo = bosons.MomentumModulo;
  this->XMomentumTranslationStep = bosons.XMomentumTranslationStep;

  this->ReducedNbrState = bosons.ReducedNbrState;
  this->RemainderNbrState = bosons.RemainderNbrState;
  this->TemporaryState.Resize(this->ReducedNbrState);
  this->TemporaryState2.Resize(this->ReducedNbrState);

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;

  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  this->LookUpTable = bosons.LookUpTable;
  this->NbrStateInLookUpTable = bosons.NbrStateInLookUpTable;
  this->HashKeyMask = bosons.HashKeyMask;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

  this->Flag = bosons.Flag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnTorusWithMagneticTranslations::Clone()
{
  return new BosonOnTorusWithMagneticTranslations(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnTorusWithMagneticTranslations::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->XMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->YMomentum, this->MomentumModulo);
  List<AbstractQuantumNumber*> TmpList2;
  TmpList2 += new VectorQuantumNumber (TmpList);
  return TmpList2;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnTorusWithMagneticTranslations::GetQuantumNumber (int index)
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->XMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->YMomentum, this->MomentumModulo);
  return new VectorQuantumNumber (TmpList);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnTorusWithMagneticTranslations::ExtractSubspace (AbstractQuantumNumber& q, 
									     SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != (AbstractQuantumNumber::Vector | AbstractQuantumNumber::PeriodicMomentum))
    return 0;
  if (((VectorQuantumNumber&) q).GetQuantumNumbers().GetNbrElement() != 2)
    return 0;
  if ((this->XMomentum == ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[0]))->GetMomentum()) &&
      (this->YMomentum == ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[1]))->GetMomentum()))
    return this;
  else 
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
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int BosonOnTorusWithMagneticTranslations::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  BosonOnTorusState& State = this->StateDescription[index];
  if ((n1 > this->StateMaxMomentum[index]) || (n2 > this->StateMaxMomentum[index]) || 
      (State.GetOccupation(n1) == 0) || (State.GetOccupation(n2) == 0) 
      || ((n1 == n2) && (State.GetOccupation(n1) == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }

  this->TemporaryState.Assign(State, this->ReducedNbrState);
  unsigned long Tmp1;
  int Tmp2;
  coefficient = this->PrecalculatedSqrt[this->TemporaryState.DecrementOccupation(n1, Tmp1, Tmp2)];
  coefficient *= this->PrecalculatedSqrt[this->TemporaryState.DecrementOccupation(n2, Tmp1, Tmp2)];  
  coefficient *= this->PrecalculatedSqrt[this->TemporaryState.IncrementOccupation(m1, Tmp1, Tmp2)];
  coefficient *= this->PrecalculatedSqrt[this->TemporaryState.IncrementOccupation(m2, Tmp1, Tmp2)];
//  this->TemporaryState.PrintState(cout, this->ReducedNbrState, this->RemainderNbrState) << endl;
  this->TemporaryState.PutInCanonicalForm(this->TemporaryState2, this->ReducedNbrState, this->RemainderNbrState, 
					  this->MaxMomentum, nbrTranslation, this->XMomentumTranslationStep, Tmp2);
//  this->TemporaryState.PrintState(cout, this->ReducedNbrState, this->RemainderNbrState) << "         " << nbrTranslation << " " << Tmp2 << " " << coefficient << endl;
  if (this->CompatibilityWithXMomentum[Tmp2] == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  Tmp2 = this->FindStateIndex(this->TemporaryState);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[Tmp2]];
  return Tmp2;
}



// find state index
//
// state = state description
// return value = corresponding index

int BosonOnTorusWithMagneticTranslations::FindStateIndex(BosonOnTorusState& state)
{
  int TmpKey = state.GetHashKey(this->ReducedNbrState, this->HashKeyMask);
  int TmpMaxMomentum = state.GetHighestIndex(this->ReducedNbrState);
  int Max = this->NbrStateInLookUpTable[TmpMaxMomentum][TmpKey] - 1;
  int* TmpLookUpTable = this->LookUpTable[TmpMaxMomentum][TmpKey];
  int Min = 0;
  int Pos = 0;
  while ((Max - Min) > 1)
    {
      Pos = (Max + Min) >> 1;
      if (state.Lesser (this->StateDescription[TmpLookUpTable[Pos]], this->ReducedNbrState))
	Max = Pos;
      else
	Min = Pos;
    }
  if (state.Different(this->StateDescription[TmpLookUpTable[Min]], this->ReducedNbrState))
    return TmpLookUpTable[Max];
  else
    return TmpLookUpTable[Min];
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusWithMagneticTranslations::PrintState (ostream& Str, int state)
{
  this->StateDescription[state].PrintState(Str, this->ReducedNbrState, this->RemainderNbrState);
  //  Str  << " position = " << FindStateIndex(this->StateDescription[state]);
/*  int Max = this->StateMaxMomentum[state];
  int i = 0;
  for (; i <= Max; ++i)
    Str << TmpState[i] << " ";
  for (; i < this->MaxMomentum; ++i)
    Str << "0 ";
  Str << " key = " << this->Keys[state] << " max momentum = " << this->MomentumMaxPosition[Max * (this->NbrBosons + 1) + TmpState[Max]]
      << " position = " << FindStateIndex(TmpState, Max) << " momentum = " << this->GetMomentumValue(state);*/
  return Str;
}

// generate all states corresponding to the constraints
// 
// return value = hilbert space dimension

int BosonOnTorusWithMagneticTranslations::GenerateStates()
{
  this->CompatibilityWithXMomentum = new bool [this->NbrMomentum];
  for (int i = 0; i <= this->MomentumModulo; ++i)
    {
      if (((i * this->XMomentum) % this->MomentumModulo) == 0)
	this->CompatibilityWithXMomentum[i * this->XMomentumTranslationStep] = true;
      else
	this->CompatibilityWithXMomentum[i * this->XMomentumTranslationStep] = false;
    }
  int Lim = this->NbrBosons * this->MaxMomentum;
  this->CompatibilityWithYMomentum = new bool [Lim + 1];
  for (int i = 0; i <= Lim; ++i)
    {
      this->CompatibilityWithYMomentum[i]= false;	  
    }
  for (int i = this->YMomentum; i <= Lim; i += this->MomentumModulo)
    {
      this->CompatibilityWithYMomentum[i]= true;	  
    }

  int MaxHilbertSpaceDimension = 0;
  int MaxTmpHilbertSpaceDimension = 0;

  for (int i = 1; i < this->MaxMomentum; ++i)
    for (int j = 1; j <= this->NbrBosons; ++j)
      {
	int Pos = this->RawGenerateStatesDryRun(this->NbrBosons - j, i - 1, 0, j * i);
	if (Pos > MaxTmpHilbertSpaceDimension)
	  MaxTmpHilbertSpaceDimension = Pos;
	MaxHilbertSpaceDimension += Pos;
      }
  cout << MaxTmpHilbertSpaceDimension << " " << MaxHilbertSpaceDimension << endl;

  this->TemporaryStateDescription = new BosonOnTorusState* [MaxTmpHilbertSpaceDimension];
  for (int i = 0; i < MaxTmpHilbertSpaceDimension; ++i)
    this->TemporaryStateDescription[i] = new BosonOnTorusState (this->ReducedNbrState);

  BosonOnTorusState** InternalTemporaryStateDescription = new BosonOnTorusState* [MaxHilbertSpaceDimension];
  int* TmpStateSymmetries = new int [MaxHilbertSpaceDimension];
  BosonOnTorusState* TmpState;
  this->HilbertSpaceDimension = -1;
  int NbrTranslation;
  int StateSymmetry;
  if (this->CompatibilityWithYMomentum[0] == true)
    {
      TmpState = new BosonOnTorusState(this->ReducedNbrState);
      TmpState->SetOccupation (0, this->NbrBosons);
      TmpState->PutInCanonicalForm(this->TemporaryState, this->ReducedNbrState, this->RemainderNbrState, 
				   this->MaxMomentum, NbrTranslation, this->XMomentumTranslationStep);
      StateSymmetry = TmpState->GetStateSymmetry(this->TemporaryState, this->ReducedNbrState, this->RemainderNbrState, 
						 this->MaxMomentum, this->XMomentumTranslationStep);
      if (this->CompatibilityWithXMomentum[StateSymmetry] == true)
	{
	  InternalTemporaryStateDescription[0] = TmpState;
	  TmpStateSymmetries[0] = StateSymmetry;
	  ++this->HilbertSpaceDimension;
	}
      else
	{
	  delete TmpState;
	}
    }
  
  for (int i = 1; i < this->MaxMomentum; ++i)
    {
      for (int j = 1; j <= this->NbrBosons; ++j)
	{
	  int Pos = this->RawGenerateStates(this->NbrBosons - j, i - 1, 0, j * i);
	  for (int k = 0; k < Pos; ++k)
	    {
	      this->TemporaryStateDescription[k]->SetOccupation (i, j);
	      this->TemporaryStateDescription[k]->PutInCanonicalForm(this->TemporaryState, this->ReducedNbrState, 
								     this->RemainderNbrState, this->MaxMomentum, 
								     NbrTranslation, this->XMomentumTranslationStep);
	      StateSymmetry = this->TemporaryStateDescription[k]->GetStateSymmetry(this->TemporaryState, 
										   this->ReducedNbrState, this->RemainderNbrState, 
										   this->MaxMomentum, 
										   this->XMomentumTranslationStep);
	      if ((this->CompatibilityWithXMomentum[StateSymmetry] == true) &&
		  ((this->HilbertSpaceDimension < 0) || 
		   (InternalTemporaryStateDescription[this->HilbertSpaceDimension]->Lesser(*(this->TemporaryStateDescription[k]), 
											   this->ReducedNbrState))))
		{
		  ++this->HilbertSpaceDimension;
		  InternalTemporaryStateDescription[this->HilbertSpaceDimension] = this->TemporaryStateDescription[k];
		  this->TemporaryStateDescription[k] = new BosonOnTorusState (this->ReducedNbrState);
		  TmpStateSymmetries[this->HilbertSpaceDimension] = StateSymmetry;
		}
	      else
		this->TemporaryStateDescription[k]->EmptyState(this->ReducedNbrState);
	    }
	}
    }  

  ++this->HilbertSpaceDimension;
  for (int i = 0; i < MaxTmpHilbertSpaceDimension; ++i)
    delete this->TemporaryStateDescription[i];
  delete[] this->TemporaryStateDescription;
  this->TemporaryStateDescription = 0;
  if (this->HilbertSpaceDimension > 0)
    {
      this->StateDescription = new BosonOnTorusState [this->HilbertSpaceDimension];
      this->NbrStateInOrbit = new int [this->HilbertSpaceDimension];
      this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  this->StateDescription[i].TransfertState(*(InternalTemporaryStateDescription[i]));
	  this->StateMaxMomentum[i] = this->StateDescription[i].GetHighestIndex(this->ReducedNbrState);
	  this->NbrStateInOrbit[i] = TmpStateSymmetries[i];
	}
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    delete InternalTemporaryStateDescription[i];
  delete[] InternalTemporaryStateDescription;
  delete[] TmpStateSymmetries;
  return this->HilbertSpaceDimension;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentMaxMomentum = momentum maximum value for bosons that are still to be placed
// pos = position in StateDescription array where to store states
// currentYMomentum = current value of the momentum in the y direction
// return value = position from which new states have to be stored

int BosonOnTorusWithMagneticTranslations::RawGenerateStates(int nbrBosons, int currentMaxMomentum, 
							    int pos, int currentYMomentum)
{
  if (nbrBosons == 0)
    {
      if (this->CompatibilityWithYMomentum[currentYMomentum] == true)
	{
	  return pos + 1;
	}
      else
	{
 	  return pos;
	}
    }
  if (currentMaxMomentum == 0)
    {
      if (this->CompatibilityWithYMomentum[currentYMomentum] == true)
	{
	  this->TemporaryStateDescription[pos]->SetOccupation (0, nbrBosons);
	  return pos + 1;
	}
      else
	return pos;
    }

  int TmpNbrBosons = 0;
  int DecCurrentMaxMomentum = currentMaxMomentum - 1;
  int TmpPos = pos;
  while (TmpNbrBosons <= nbrBosons)
    {
      TmpPos = this->RawGenerateStates(nbrBosons - TmpNbrBosons, DecCurrentMaxMomentum, pos, currentYMomentum + (TmpNbrBosons * currentMaxMomentum));
      for (int i = pos; i < TmpPos; i++)
	{
	  this->TemporaryStateDescription[i]->SetOccupation (currentMaxMomentum, TmpNbrBosons);
	}      
      ++TmpNbrBosons;
      pos = TmpPos;
    }
  return pos;
}

// generate all states corresponding to the constraints without storage
// 
// nbrBosons = number of bosons
// currentMaxMomentum = momentum maximum value for bosons that are still to be placed
// pos = position in StateDescription array where to store states
// currentYMomentum = current value of the momentum in the y direction
// return value = position from which new states have to be stored

int BosonOnTorusWithMagneticTranslations::RawGenerateStatesDryRun(int nbrBosons, int currentMaxMomentum, 
								  int pos, int currentYMomentum)
{
  if (nbrBosons == 0)
    {
      if (this->CompatibilityWithYMomentum[currentYMomentum] == true)
	{
	  return pos + 1;
	}
      else
	{
 	  return pos;
	}
    }
  if (currentMaxMomentum == 0)
    {
      if (this->CompatibilityWithYMomentum[currentYMomentum] == true)
	{
	  return pos + 1;
	}
      else
	return pos;
    }

  int TmpNbrBosons = 0;
  int DecCurrentMaxMomentum = currentMaxMomentum - 1;
  int TmpPos = pos;
  while (TmpNbrBosons <= nbrBosons)
    {
      TmpPos = this->RawGenerateStatesDryRun(nbrBosons - TmpNbrBosons, DecCurrentMaxMomentum, pos, currentYMomentum + (TmpNbrBosons * currentMaxMomentum));
      ++TmpNbrBosons;
      pos = TmpPos;
    }
  return pos;
}

// generate look-up table associated to current Hilbert space
// 

void BosonOnTorusWithMagneticTranslations::GenerateLookUpTable()
{
  this->RescalingFactors = new double* [this->NbrMomentum];
  for (int i = 1; i <= this->MaxMomentum; ++i)
    {
      this->RescalingFactors[i] = new double [this->NbrMomentum];
      for (int j = 1; j <= this->MaxMomentum; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }

  this->PrecalculatedSqrt = new double [this->IncNbrBosons];
  for (int i = 0; i < this->IncNbrBosons; ++i)
    this->PrecalculatedSqrt[i] = sqrt ((double) i);

  this->LookUpTable = new int** [this->MaxMomentum];
  int CurrentMaxMomentum = this->StateMaxMomentum[0];
  this->LookUpTable[CurrentMaxMomentum] = new int* [this->HashKeyMask + (unsigned long) 1];
  this->NbrStateInLookUpTable = new int* [this->MaxMomentum];
  for (int i = 0; i < this->MaxMomentum; ++i)
    this->NbrStateInLookUpTable[i] = 0;
  this->NbrStateInLookUpTable[CurrentMaxMomentum] = new int [this->HashKeyMask + (unsigned long) 1];
  for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
    this->NbrStateInLookUpTable[CurrentMaxMomentum][j] = 0;
  unsigned long* TmpHashTable = new unsigned long [this->HilbertSpaceDimension];
  unsigned long TmpKey = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentMaxMomentum != this->StateMaxMomentum[i])
	{
	  CurrentMaxMomentum = this->StateMaxMomentum[i];
	  this->LookUpTable[CurrentMaxMomentum] = new int* [this->HashKeyMask + (unsigned long) 1];
	  this->NbrStateInLookUpTable[CurrentMaxMomentum] = new int [this->HashKeyMask + (unsigned long) 1];
	  TmpKey = this->StateDescription[i].GetHashKey(this->ReducedNbrState, this->HashKeyMask);
	  TmpHashTable[i] = TmpKey;
	  for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
	    this->NbrStateInLookUpTable[CurrentMaxMomentum][j] = 0;
	  ++this->NbrStateInLookUpTable[CurrentMaxMomentum][TmpKey];	  
	}
      else
	{
	  TmpKey = this->StateDescription[i].GetHashKey(this->ReducedNbrState, this->HashKeyMask);
	  TmpHashTable[i] = TmpKey;
	  ++this->NbrStateInLookUpTable[CurrentMaxMomentum][TmpKey];
	}
    }
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      if (this->NbrStateInLookUpTable[i] != 0)
	{
	  for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
	    {
	      if (this->NbrStateInLookUpTable[i][j] != 0)
		{
		  this->LookUpTable[i][j] = new int [this->NbrStateInLookUpTable[i][j]];
		  this->NbrStateInLookUpTable[i][j] = 0;
		}
	      else
		{
		  this->LookUpTable[i][j] = 0;
		}
	    }
	}
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    this->LookUpTable[this->StateMaxMomentum[i]][TmpHashTable[i]][this->NbrStateInLookUpTable[this->StateMaxMomentum[i]]
								  [TmpHashTable[i]]++] = i;
  delete[] TmpHashTable;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// return value = Hilbert space dimension

int BosonOnTorusWithMagneticTranslations::EvaluateHilbertSpaceDimension(int nbrBosons, int maxMomentum)
{
  FactorialCoefficient Dimension; 
  Dimension.PartialFactorialMultiply(maxMomentum, maxMomentum + nbrBosons - 1); 
  Dimension.FactorialDivide(nbrBosons);
  return Dimension.GetIntegerValue();
}

