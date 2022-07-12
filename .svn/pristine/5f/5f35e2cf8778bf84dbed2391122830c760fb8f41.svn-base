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
#include "Matrix/Matrix.h"

#include <math.h>


// basic constructor
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion

FermionOnTorus::FermionOnTorus (int nbrFermions, int maxMomentum)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum);
  this->Flag.Initialize();
  this->MomentumConstraint = 0;
  this->MomentumConstraintFlag = false;
  this->StateDescription = new unsigned int [this->HilbertSpaceDimension];
  this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrFermions, this->MaxMomentum - 1, this->MaxMomentum - 1, 0);
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
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->MomentumConstraint = momentumConstraint;
  this->MomentumConstraintFlag = true;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum);
  this->Flag.Initialize();
  this->StateDescription = new unsigned int [this->HilbertSpaceDimension];
  this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->MaxMomentum - 1, this->MaxMomentum - 1, 0, 0);
  cout << this->HilbertSpaceDimension << endl;
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
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;
  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->MomentumConstraint = fermions.MomentumConstraint;
  this->MomentumConstraintFlag = fermions.MomentumConstraintFlag;
  this->Flag = fermions.Flag;
}

// constructor from full datas (with no constraint on the total momentum)
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// hilbertSpaceDimension = Hilbert space dimension
// stateDescription = array describing each state
// stateMaxMomentum = array giving maximum Lz value reached for a fermion in a given state

FermionOnTorus::FermionOnTorus (int nbrFermions, int maxMomentum, int hilbertSpaceDimension, 
				unsigned int* stateDescription, int* stateMaxMomentum)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->MomentumConstraintFlag = false;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = stateDescription;
  this->StateMaxMomentum = stateMaxMomentum;
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
// stateMaxMomentum = array giving maximum Lz value reached for a fermion in a given state

FermionOnTorus::FermionOnTorus (int nbrFermions, int maxMomentum, int momentumConstraint, int hilbertSpaceDimension, 
				unsigned int* stateDescription, int* stateMaxMomentum)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->MomentumConstraint = momentumConstraint;
  this->MomentumConstraintFlag = true;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = stateDescription;
  this->StateMaxMomentum = stateMaxMomentum;
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
      delete[] this->StateMaxMomentum;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnTorus& FermionOnTorus::operator = (const FermionOnTorus& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
    }
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;
  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->MomentumConstraint = fermions.MomentumConstraint;
  this->MomentumConstraintFlag = fermions.MomentumConstraintFlag;
  this->Flag = fermions.Flag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnTorus::Clone()
{
  return new FermionOnTorus(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnTorus::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->MomentumConstraintFlag == false)
    {
      for (int i = 0; i < this->MaxMomentum; ++i)
	TmpList += new PeriodicMomentumQuantumNumber (i, this->MaxMomentum);
    }
  else
    TmpList += new PeriodicMomentumQuantumNumber (this->MomentumConstraint, this->MaxMomentum);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnTorus::GetQuantumNumber (int index)
{
  if (this->MomentumConstraintFlag == false)
    {
      return  new PeriodicMomentumQuantumNumber (this->GetMomentumValue(index), this->MaxMomentum);
    }
  else
    return new PeriodicMomentumQuantumNumber (this->MomentumConstraint, this->MaxMomentum);
}

// get momemtum value of a given state
//
// index = state index
// return value = state momentum

int FermionOnTorus::GetMomentumValue(int index)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  int Momentum = 0;
  for (int i = 0; i <= StateMaxMomentum; ++i)
    {
      Momentum += ((State >> i ) & 1) * i;
    }
  return (Momentum % this->MaxMomentum);
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
  if (this->MomentumConstraintFlag == true)
    if (this->MomentumConstraint == ((PeriodicMomentumQuantumNumber&) q).GetMomentum())
      return this;
    else 
      return 0;
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
  unsigned int* SubspaceStateDescription = new unsigned int [SubspaceHilbertSpaceDimension];
  int* SubspaceStateMaxMomentum = new int [SubspaceHilbertSpaceDimension];
  int* ConvArray = new int [SubspaceHilbertSpaceDimension];
  for (int i = 0; i < SubspaceHilbertSpaceDimension; ++i)
    {
      ConvArray[i] = TmpConvArray[i];
      SubspaceStateDescription[i] = this->StateDescription[TmpConvArray[i]];
      SubspaceStateMaxMomentum[i] = this->StateMaxMomentum[TmpConvArray[i]];
    }
  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, SubspaceHilbertSpaceDimension, ConvArray);
  return new FermionOnTorus(this->NbrFermions, this->MaxMomentum, Momentum, SubspaceHilbertSpaceDimension, 
			    SubspaceStateDescription, SubspaceStateMaxMomentum);
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
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
/*  m1 += this->MaxMomentum;
  m1 >>= 1;
  m2 += this->MaxMomentum;
  m2 >>= 1;
  n1 += this->MaxMomentum;
  n1 >>= 1;
  n2 += this->MaxMomentum;
  n2 >>= 1;*/
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (1 << n1)) == 0) || 
      ((State & (1 << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  int TmpState = State;
  int Mask;
  int SignIndex = NewMaxMomentum - n2;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> n2);
      coefficient = this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> n2) & 0xffff;
      coefficient = this->SignLookUpTable[Mask];
      Mask = (TmpState >> (n2 + 16));
      coefficient *= this->SignLookUpTable[Mask];
    }
  TmpState &= ~(0x1 << n2);
  if (NewMaxMomentum == n2)
    while ((TmpState >> NewMaxMomentum) == 0)
      --NewMaxMomentum;
  SignIndex = NewMaxMomentum - n1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> n1);
      coefficient *= this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> n1) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (n1 + 16));
      coefficient *= this->SignLookUpTable[Mask];
    }
  TmpState &= ~(0x1 << n1);
  if (NewMaxMomentum == n1)
    while ((TmpState >> NewMaxMomentum) == 0)
      --NewMaxMomentum;
  if ((TmpState & (1 << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      TmpState |= (0x1 << m2);
      NewMaxMomentum = m2;
    }
  else
    {
      SignIndex = NewMaxMomentum - m2;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> m2);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> m2) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (m2 + 16));
	  coefficient *= this->SignLookUpTable[Mask];
	}
      TmpState |= (0x1 << m2);
    }
  if ((TmpState & (1 << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      TmpState |= (0x1 << m1);
      NewMaxMomentum = m1;
    }
  else
    {
      SignIndex = NewMaxMomentum - m1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> m1);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> m1) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (m1 + 16));
	  coefficient *= this->SignLookUpTable[Mask];
	}
      TmpState |= (0x1 << m1);
    }
//  cout << hex << TmpState << dec << endl;
  return this->FindStateIndex(TmpState, NewMaxMomentum);
}

// apply a^+_m a_m operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for creation operator
// return value =  resulting multiplicative factor 

double FermionOnTorus::AdA (int index, int m)
{
  if (this->StateMaxMomentum[index] < m)
    return 0.0;
  if ((this->StateDescription[index] & (0x1 << m)) == 0)
    return 0.0;
  else
    return 1.0;
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
  int StateMaxMomentum;
  unsigned int State;
  double Coefficient;
  unsigned int GlobalMask = (0x1 << i);
  int Mask;
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)
    {
      StateMaxMomentum = this->StateMaxMomentum[j];
      if (StateMaxMomentum >= i)
	{
	  State = this->StateDescription[j];
	  if ((State & GlobalMask) != 0)
	    {
	      Mask = (State >> i) & 0xffff;
	      Coefficient = this->SignLookUpTable[Mask];
	      Mask = (State >> (i + 16));
	      Coefficient *= this->SignLookUpTable[Mask];
	      State &= ~GlobalMask;
	      if (StateMaxMomentum == i)
		while ((State >> StateMaxMomentum) == 0)
		  --StateMaxMomentum;
	      M(this->FindStateIndex(State, StateMaxMomentum), j) = Coefficient;
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
  int StateMaxMomentum;
  unsigned int State;
  double Coefficient;
  int Mask;
  unsigned int GlobalMask = (0x1 << i);
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)
    {
      StateMaxMomentum = this->StateMaxMomentum[j];
      State = this->StateDescription[j];
      if ((State & GlobalMask) == 0)
	{
	  if (i > StateMaxMomentum)
	    {
	      State |= GlobalMask;
	      M(this->FindStateIndex(State, i), j) = 1.0;
	    }
	  else
	    {
	      Mask = (State >> i) & 0xffff;
	      Coefficient = this->SignLookUpTable[Mask];
	      Mask = (State >> (i + 16));
	      Coefficient *= this->SignLookUpTable[Mask];
	      State |= GlobalMask;
	      M(this->FindStateIndex(State, StateMaxMomentum), j) = Coefficient;
	    }
	}
    }
  return M;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnTorus::FindStateIndex(unsigned int stateDescription, int maxMomentum)
{
//  int Pos = 0;
  int Pos = this->LookUpTable[maxMomentum][stateDescription >> this->LookUpTableShift[maxMomentum]];
//  cout << maxMomentum << " " << hex << stateDescription << dec << " " << this->LookUpTableShift[maxMomentum] << Pos << endl;
  while (this->StateDescription[Pos] != stateDescription)
    ++Pos;
/*  if (this->StateMaxMomentum[Pos] != maxMomentum)
    cout << "error !!! " << maxMomentum << " " << this->StateMaxMomentum[Pos] << " " << hex << stateDescription << dec << " " << Pos << endl;*/
  return Pos;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnTorus::PrintState (ostream& Str, int state)
{
  int TmpState = this->StateDescription[state];
  for (int i = 0; i < this->MaxMomentum; ++i)
    Str << ((TmpState >> i) & 0x1) << " ";
//  Str << " key = " << this->Keys[state] << " maxMomentum position = " << this->MaxMomentumPosition[Max * (this->NbrFermions + 1) + TmpState[Max]]
  Str << " position = " << FindStateIndex(TmpState, this->StateMaxMomentum[state]) << " momentum = " << this->GetMomentumValue(state);
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int FermionOnTorus::GenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos)
{
  if ((nbrFermions == 0) || (nbrFermions > (currentMaxMomentum + 1)) || (currentMaxMomentum < 0))
    return pos;
  if (nbrFermions == 1)
    {
      for (int i = currentMaxMomentum; i >= 0; --i)
	{
	  this->StateDescription[pos] = 1 << i;
	  this->StateMaxMomentum[pos] = maxMomentum;
	  ++pos;
	}
      return pos;
    }
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentMaxMomentum, pos);
  unsigned int Mask = 1 << currentMaxMomentum;
  for (int i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (maxMomentum == currentMaxMomentum)
    return this->GenerateStates(nbrFermions, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, TmpPos);
  else
    return this->GenerateStates(nbrFermions, maxMomentum, ReducedCurrentMaxMomentum, TmpPos);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

int FermionOnTorus::GenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos, int currentMomentum)
{
  if ((nbrFermions == 0) || (nbrFermions > (currentMaxMomentum + 1)))
    return pos;
//  cout << nbrFermions << " " << maxMomentum << " " << currentMaxMomentum << " " << pos << " " << currentMomentum << endl;
  if (nbrFermions == 1)
    {
      int i = this->MomentumConstraint - (currentMomentum % this->MaxMomentum);
      if (i < 0)
	i += this->MaxMomentum;
      for (; i <= currentMaxMomentum; i += this->MaxMomentum)
	{
	  this->StateDescription[pos] = 1 << i;
	  this->StateMaxMomentum[pos] = maxMomentum;
	  ++pos;
	}
      return pos;
    }
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentMaxMomentum, pos, currentMomentum + currentMaxMomentum);
  unsigned int Mask = 1 << currentMaxMomentum;
  for (int i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (maxMomentum == currentMaxMomentum)
    return this->GenerateStates(nbrFermions, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentMomentum);
  else
    return this->GenerateStates(nbrFermions, maxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentMomentum);
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
  int CurrentMaxMomentum = this->StateMaxMomentum[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
  if (CurrentMaxMomentum < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentMaxMomentum] = 0;
  else
    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
  unsigned int CurrentLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  unsigned int TmpLookUpTableValue;
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentMaxMomentum != this->StateMaxMomentum[i])
	{
 	  CurrentMaxMomentum = this->StateMaxMomentum[i];
	  TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
	  if (CurrentMaxMomentum < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentMaxMomentum] = 0;
	  else
	    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
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

