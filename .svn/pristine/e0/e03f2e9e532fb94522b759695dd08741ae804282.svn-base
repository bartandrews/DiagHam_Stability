////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of bosons on torus with magnetic translations          //
//                         and for system size such that                      //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 30/01/2012                      //
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
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "Architecture/ArchitectureOperation/FQHETorusParticleEntanglementSpectrumOperation.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "MathTools/IntegerAlgebraTools.h"

#include <math.h>
#include <algorithm>
#include <set>


using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
// 

BosonOnTorusWithMagneticTranslationsShort::BosonOnTorusWithMagneticTranslationsShort ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// kxMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// kyMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)

BosonOnTorusWithMagneticTranslationsShort::BosonOnTorusWithMagneticTranslationsShort (int nbrBosons, int maxMomentum, int kxMomentum, int kyMomentum)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = maxMomentum;
  this->FermionicMaxMomentum = this->MaxMomentum + this->NbrBosons - 1;
  this->NbrKyValue = this->MaxMomentum + 1;
 
  this->MomentumModulo = FindGCD(this->NbrBosons, this->MaxMomentum);
  this->KxMomentum = kxMomentum % this->MomentumModulo;
  this->KyMomentum = kyMomentum % this->MaxMomentum;

  this->StateShift = this->MaxMomentum / this->MomentumModulo;
  this->LastMomentumMask = 0x1ul << (this->MaxMomentum + this->NbrBosons - 1);

  this->TemporaryState = new unsigned long [this->NbrKyValue];
  this->ProdATemporaryState = new unsigned long [this->NbrKyValue];

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->MaxMomentum - 1, 0);
  if (this->LargeHilbertSpaceDimension != 0l)
    {
      long TmpLargeHilbertSpaceDimension = this->GenerateStates();
      if (TmpLargeHilbertSpaceDimension != 0l)
	{
	  this->LargeHilbertSpaceDimension = TmpLargeHilbertSpaceDimension;
	  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
	  this->Flag.Initialize();
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
    }
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnTorusWithMagneticTranslationsShort::BosonOnTorusWithMagneticTranslationsShort(const BosonOnTorusWithMagneticTranslationsShort& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrKyValue = bosons.NbrKyValue;
  this->MomentumModulo = bosons.MomentumModulo;
  this->FermionicMaxMomentum = bosons.FermionicMaxMomentum;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMask = bosons.LastMomentumMask;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->Flag = bosons.Flag;
  this->TemporaryState = new unsigned long [this->NbrKyValue];
  this->ProdATemporaryState = new unsigned long [this->NbrKyValue];
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

}

// destructor
//

BosonOnTorusWithMagneticTranslationsShort::~BosonOnTorusWithMagneticTranslationsShort ()
{
  if ((this->LargeHilbertSpaceDimension != 0l) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorusWithMagneticTranslationsShort& BosonOnTorusWithMagneticTranslationsShort::operator = (const BosonOnTorusWithMagneticTranslationsShort& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
    }
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->MaxMomentum = bosons.MaxMomentum;
  this->FermionicMaxMomentum = bosons.FermionicMaxMomentum;
  this->NbrKyValue = bosons.NbrKyValue;
  this->MomentumModulo = bosons.MomentumModulo;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMask = bosons.LastMomentumMask;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->TemporaryState = new unsigned long [this->NbrKyValue];
  this->ProdATemporaryState = new unsigned long [this->NbrKyValue];
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->Flag = bosons.Flag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnTorusWithMagneticTranslationsShort::Clone()
{
  return new BosonOnTorusWithMagneticTranslationsShort(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnTorusWithMagneticTranslationsShort::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->KxMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->KyMomentum, this->MomentumModulo);
  List<AbstractQuantumNumber*> TmpList2;
  TmpList2 += new VectorQuantumNumber (TmpList);
  return TmpList2;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnTorusWithMagneticTranslationsShort::GetQuantumNumber (int index)
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->KxMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->KyMomentum, this->MomentumModulo);
  return new VectorQuantumNumber (TmpList);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnTorusWithMagneticTranslationsShort::ExtractSubspace (AbstractQuantumNumber& q, 
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
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int BosonOnTorusWithMagneticTranslationsShort::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  this->FermionToBoson(this->StateDescription[index], this->FermionicMaxMomentum, this->TemporaryState, this->TemporaryStateKyMax);
  if ((n1 > this->TemporaryStateKyMax) || (n2 > this->TemporaryStateKyMax) || (this->TemporaryState[n1] == 0) || (this->TemporaryState[n2] == 0) || ((n1 == n2) && (this->TemporaryState[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  for (int i = this->TemporaryStateKyMax + 1; i < this->MaxMomentum; ++i)
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
  this->TemporaryStateKyMax =  this->MaxMomentum - 1;
  while (this->TemporaryState[this->TemporaryStateKyMax] == 0x0ul)
    --this->TemporaryStateKyMax;
  unsigned long TmpState = this->FindCanonicalFormAndTestXMomentumConstraint(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), nbrTranslation);
  if (nbrTranslation < 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = this->FermionicMaxMomentum;
  while ((TmpState >> NewMaxMomentum) == 0x0ul)
    --NewMaxMomentum;
  int TmpIndex = this->FindStateIndex(TmpState, NewMaxMomentum);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
  nbrTranslation *= this->StateShift;
  return TmpIndex;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnTorusWithMagneticTranslationsShort::AA (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription[index], this->FermionicMaxMomentum, this->ProdATemporaryState, this->ProdATemporaryStateKyMax);
  if ((n1 > this->ProdATemporaryStateKyMax) || (n2 > this->ProdATemporaryStateKyMax) || (this->ProdATemporaryState[n1] == 0) || 
      (this->ProdATemporaryState[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryState[n1] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryState[n2];
  --this->ProdATemporaryState[n2];
  Coefficient *= this->ProdATemporaryState[n1];
  --this->ProdATemporaryState[n1];
  for (int i = this->ProdATemporaryStateKyMax + 1; i < this->MaxMomentum; ++i)
    this->ProdATemporaryState[i] = 0;
  this->ProdATemporaryStateNbrStateInOrbit = this->NbrStateInOrbit[index];
  return sqrt(Coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnTorusWithMagneticTranslationsShort::ProdA (int index, int* n, int nbrIndices)
{
  this->FermionToBoson(this->StateDescription[index], this->FermionicMaxMomentum, this->ProdATemporaryState, this->ProdATemporaryStateKyMax);
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
  for (int i = this->ProdATemporaryStateKyMax + 1; i < this->MaxMomentum; ++i)
    this->ProdATemporaryState[i] = 0;
  this->ProdATemporaryStateNbrStateInOrbit = this->NbrStateInOrbit[index];
  return sqrt((double) TmpCoefficient);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int BosonOnTorusWithMagneticTranslationsShort::AdAd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  for (int i = 0; i < this->MaxMomentum; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m2];
  coefficient = this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  this->TemporaryStateKyMax = this->MaxMomentum - 1;
  while (this->TemporaryState[this->TemporaryStateKyMax] == 0)
    --this->TemporaryStateKyMax;
  unsigned long TmpState = this->FindCanonicalFormAndTestXMomentumConstraint(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), nbrTranslation);
  if (nbrTranslation < 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = this->FermionicMaxMomentum;
  while ((TmpState >> NewMaxMomentum) == 0x0ul)
    --NewMaxMomentum;
  int TmpIndex = this->FindStateIndex(TmpState, NewMaxMomentum);
  coefficient *= this->RescalingFactors[this->ProdATemporaryStateNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
  nbrTranslation *= this->StateShift;  
  return TmpIndex;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int BosonOnTorusWithMagneticTranslationsShort::ProdAd (int* m, int nbrIndices, double& coefficient, int& nbrTranslation)
{
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      this->TemporaryState[i] = this->ProdATemporaryState[i];
    }
  int TmpCoefficient = 1;
  for (int i = 0; i < nbrIndices; ++i)
    TmpCoefficient *= ++this->TemporaryState[m[i]];
  coefficient = sqrt((double) TmpCoefficient);
  this->TemporaryStateKyMax = this->MaxMomentum - 1;
  while (this->TemporaryState[this->TemporaryStateKyMax] == 0)
    --this->TemporaryStateKyMax;
  unsigned long TmpState = this->FindCanonicalFormAndTestXMomentumConstraint(this->BosonToFermion(this->TemporaryState, this->TemporaryStateKyMax), nbrTranslation);
  if (nbrTranslation < 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = this->FermionicMaxMomentum;
  while ((TmpState >> NewMaxMomentum) == 0x0ul)
    --NewMaxMomentum;
  int TmpIndex = this->FindStateIndex(TmpState, NewMaxMomentum);
  coefficient *= this->RescalingFactors[this->ProdATemporaryStateNbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
  nbrTranslation *= this->StateShift;
  return TmpIndex;
}

// apply a^+_m a_m operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for creation operator
// return value =  resulting multiplicative factor 

double BosonOnTorusWithMagneticTranslationsShort::AdA (int index, int m)
{
  this->FermionToBoson(this->StateDescription[index], this->FermionicMaxMomentum, this->ProdATemporaryState, this->ProdATemporaryStateKyMax);
  return this->ProdATemporaryState[m];
}

// apply a_n operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n = index for annihilation operator
// return value =  multiplicative factor 

double BosonOnTorusWithMagneticTranslationsShort::A (int index, int n)
{
  this->FermionToBoson(this->StateDescription[index], this->FermionicMaxMomentum, this->ProdATemporaryState, this->ProdATemporaryStateKyMax);
  
  if ((n > this->ProdATemporaryStateKyMax) || (this->ProdATemporaryState[n] == 0))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryState[n];
  --this->ProdATemporaryState[n];
  for (int i = this->ProdATemporaryStateKyMax + 1; i < this->MaxMomentum; ++i)
    this->ProdATemporaryState[i] = 0;
  this->ProdATemporaryStateNbrStateInOrbit = this->NbrStateInOrbit[index];
  return Coefficient;
}

// return matrix representation of the annihilation operator a_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& BosonOnTorusWithMagneticTranslationsShort::A (int i, Matrix& M)
{
  return M;
}

// return matrix representation of the creation operator a^+_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& BosonOnTorusWithMagneticTranslationsShort::Ad (int i, Matrix& M)
{
  return M;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusWithMagneticTranslationsShort::PrintStateMonomial (ostream& Str, long state)
{
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  cout << this->StateDescription[state] << endl;
  this->ConvertToMonomial(this->StateDescription[state], this->MaxMomentum + this->NbrBosons - 1, TmpMonomial);
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

int BosonOnTorusWithMagneticTranslationsShort::FindStateIndex(unsigned long stateDescription, int lzmax)
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
 
// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusWithMagneticTranslationsShort::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescription[state], this->MaxMomentum + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
  int i = 0;
  for (; i <= this->TemporaryStateKyMax; ++i)
    Str << this->TemporaryState[i] << " ";
  for (; i < this->MaxMomentum; ++i)
    Str << "0 ";
 return Str;
}


// generate all states with both the kx and ky constraint
// 
// return value = new dimension of the Hilbert space

long BosonOnTorusWithMagneticTranslationsShort::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  long TmpLargeHilbertSpaceDimension = this->RawGenerateStates(this->NbrBosons, this->MaxMomentum - 1, 0l, 0);
  if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
    {
      cout << "error : dimension mismatch while generating the Hilbert space ( is" << TmpLargeHilbertSpaceDimension << ", should be " 
	   << this->LargeHilbertSpaceDimension << ")" << endl;
      return 0l;
    }

  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      int NbrTranslation = 0;
      unsigned long TmpState = this->FindCanonicalFormAndTestXMomentumConstraint(this->StateDescription[i], NbrTranslation);
      if (NbrTranslation == 0)
	{
	  ++TmpLargeHilbertSpaceDimension;
	}
      else
	{
	  this->StateDescription[i] = 0x0ul;
	}
    }

  unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
  this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->StateDescription[i] != 0x0ul)
	{
	  TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	  unsigned long TmpState = this->StateDescription[i];
	  unsigned long TmpReferenceState = TmpState;
	  int TmpOrbitSize = 1;
	  this->ApplySingleTranslation(TmpState);
	  while (TmpState != TmpReferenceState)
	    {
	      this->ApplySingleTranslation(TmpState);
	      ++TmpOrbitSize;
	    }
	  this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = TmpOrbitSize;
	  ++TmpLargeHilbertSpaceDimension;
	}	
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  return TmpLargeHilbertSpaceDimension;
}

// generate all states corresponding to the ky constraint  without taking care of the kx constraint
// 
// nbrBosons = number of bosons
// currentKyMax = momentum maximum value for bosons that are still to be placed
// pos = position in StateDescription array where to store states
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

long BosonOnTorusWithMagneticTranslationsShort::RawGenerateStates(int nbrBosons, int currentKyMax, long pos, int currentMomentum)
{
  if (nbrBosons == 0)
    {
      if ((currentMomentum % this->MaxMomentum) == this->KyMomentum)
	{
	  this->StateDescription[pos] = 0x0ul;
	  return pos + 1l;
	}
      else
	{
	  return pos;
	}
    }

  if (currentKyMax < 0)
    {
      return pos;   
    }

  long TmpPos = pos;
  for (int i = nbrBosons; i > 0; --i)
    {
      TmpPos = this->RawGenerateStates(nbrBosons - i, currentKyMax - 1, pos, currentMomentum + (i * currentKyMax));
      unsigned long Mask = ((0x1ul << i) - 0x1ul) << (currentKyMax + nbrBosons - i);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  return  this->RawGenerateStates(nbrBosons, currentKyMax - 1, pos, currentMomentum);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnTorusWithMagneticTranslationsShort::GenerateLookUpTable(int memory)
{
  // evaluate look-up table size
  int TmpNbrKyValue =  this->MaxMomentum + this->NbrBosons;
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
  this->LookUpTable = new int* [TmpNbrKyValue];
  this->LookUpTableShift = new int [TmpNbrKyValue];
  for (int i = 0; i < TmpNbrKyValue; ++i)
{
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
    this->LookUpTableShift[i] = -1;
}
  int CurrentLzMax = this->MaxMomentum + this->NbrBosons - 1;
  while ((this->StateDescription[0] >> CurrentLzMax) == 0x0ul)
    --CurrentLzMax;
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

  this->RescalingFactors = new double* [this->MaxMomentum + 1];
  for (int i = 1; i <= this->MaxMomentum; ++i)
    {
      this->RescalingFactors[i] = new double [this->MaxMomentum + 1];
      for (int j = 1; j <= this->MaxMomentum; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKyMax = momentum maximum value for bosons that are still to be placed
// currentMomentum = current value of the momentum
// return value = Hilbert space dimension

long BosonOnTorusWithMagneticTranslationsShort::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKyMax, int currentMomentum)
{
  if (nbrBosons == 0)
    {
      if ((currentMomentum % this->MaxMomentum) == this->KyMomentum)
	{
	  return 1l;
	}
      else
	{
	  return 0l;
	}
    }

  if (currentKyMax < 0)
    {
      return 0l;   
    }

  long TmpNbrStates = 0l;
  for (int i = nbrBosons; i >= 0; --i)
    {
      TmpNbrStates += this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKyMax - 1, currentMomentum + (i * currentKyMax));
    }

  return TmpNbrStates;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

// RealSymmetricMatrix BosonOnTorusWithMagneticTranslationsShort::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int kySector, RealVector& groundState)
// {
//   if (subsytemSize <= 0)
//     {
//       if ((kySector == 0) && (nbrBosonSector == 0))
// 	{
// 	  RealSymmetricMatrix TmpDensityMatrix(1);
// 	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
// 	  return TmpDensityMatrix;
// 	}
//       else
// 	{
// 	  RealSymmetricMatrix TmpDensityMatrix;
// 	  return TmpDensityMatrix;	  
// 	}
//     }
//   if (subsytemSize > this->MaxMomentum)
//     {
//       if ((kySector == this->TotalKy) && (nbrBosonSector == this->NbrBosons))
// 	{
// 	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
// 	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
// 	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * groundState[j]);
// 	}
//       else
// 	{
// 	  RealSymmetricMatrix TmpDensityMatrix;
// 	  return TmpDensityMatrix;	  
// 	}
//     }

//   long TmpNbrNonZeroElements = 0;
//   BosonOnTorusWithMagneticTranslationsShort TmpDestinationHilbertSpace(nbrBosonSector, subsytemSize - 1, kySector);
//   cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
//   RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
//   BosonOnTorusWithMagneticTranslationsShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, this->MaxMomentum - subsytemSize, 0);
//   int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];

//   for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
//     {
//       int Pos = 0;
//       //      for (int i = 0; i <= TmpHilbertSpace.StateKyMax[MinIndex]; ++i)
//       //	this->TemporaryState[i] = TmpHilbertSpace.StateDescription[MinIndex][i];
//       for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
// 	{

// // 	  unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
// // 	  int TmpKyMax = this->FermionBasis->MaxMomentum;
// // 	  while (((TmpState >> TmpKyMax) & 0x1ul) == 0x0ul)
// // 	    --TmpKyMax;
// // 	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpKyMax);
// // 	  if (TmpPos != this->HilbertSpaceDimension)
// // 	    {
// // 	      TmpStatePosition[Pos] = TmpPos;
// // 	      TmpStatePosition2[Pos] = j;
// // 	      ++Pos;
// // 	    }
// 	}

//       if (Pos != 0)
// 	{
// 	  ++TmpNbrNonZeroElements;
// 	  for (int j = 0; j < Pos; ++j)
// 	    {
// 	      int Pos2 = TmpStatePosition2[j];
// 	      double TmpValue = groundState[TmpStatePosition[j]];
// 	      for (int k = 0; k < Pos; ++k)
// 		if (TmpStatePosition2[k] >= Pos2)
// 		TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
// 	    }
// 	}
//     }
  
//   delete[] TmpStatePosition2;
//   delete[] TmpStatePosition;
//   if (TmpNbrNonZeroElements > 0)	
//     return TmpDensityMatrix;
//   else
//     {
//       RealSymmetricMatrix TmpDensityMatrixZero;
//       return TmpDensityMatrixZero;
//     }
// }

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given (Kx,Ky) sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// kxSector = Kx sector in which the density matrix has to be evaluated 
// kySector = Ky sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnTorusWithMagneticTranslationsShort::EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kxSector, int kySector, ComplexVector& groundState)
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

//   int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;
//   int ComplementaryKySector = this->TotalKy - kySector;
//   if (ComplementaryKySector < 0)
//     ComplementaryKySector += this->MaxMomentum;
//   if (ComplementaryKySector >= this->MaxMomentum)
//     ComplementaryKySector -= this->MaxMomentum;

//   if (nbrBosonSector == 1)
//     {
//       double TmpValue = 0.0;
//       unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
//       unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
//       cout << "ComplementaryKySector = " << ComplementaryKySector << endl;
//       BosonOnTorusWithMagneticTranslationsShort TmpHilbertSpace(this->NbrBosons - 1, this->MaxMomentum, ComplementaryKySector);
//       FactorialCoefficient Factorial;
// //       for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
// // 	TmpHilbertSpace.PrintState(cout, MinIndex) << " | " << TmpHilbertSpace.StateKyMax[MinIndex] << " | " << hex << TmpHilbertSpace.StateDescription[MinIndex] << dec << endl;
//       for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
// 	{
// 	  TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.StateKyMax[MinIndex] + TmpHilbertSpace.NbrBosons - 1, TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateKyMax);
// 	  TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.TemporaryStateKyMax + TmpHilbertSpace.NbrBosons - 1, TmpMonomial1);
// 	  int TmpIndex2 = 0;
// 	  int TmpIndex4 = 0;
// 	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial1[TmpIndex2] >= kySector))
// 	    {
// 	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
// 	      ++TmpIndex2;
// 	      ++TmpIndex4;		  
// 	    }
// 	  TmpMonomial3[TmpIndex4] = kySector;
// 	  ++TmpIndex4;
// 	  while (TmpIndex2 < ComplementaryNbrBosonSector)
// 	    {
// 	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
// 	      ++TmpIndex2;
// 	      ++TmpIndex4;		  
// 	    }
// 	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
// 	  int TmpPos = this->FindStateIndex(TmpState, TmpMonomial3[0] + this->NbrBosons - 1);
// 	  if (TmpPos != this->HilbertSpaceDimension)
// 	    {
// 	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
// 	      Factorial.SetToOne();
// 	      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateKyMax; ++k)
// 		if (TmpHilbertSpace.TemporaryState[k] > 1)
// 		  Factorial.FactorialDivide(TmpHilbertSpace.TemporaryState[k]);
// 	      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
// 		if (this->TemporaryState[k] > 1)
// 		  Factorial.FactorialMultiply(this->TemporaryState[k]);
// 	      Factorial.BinomialDivide(this->NbrBosons, nbrBosonSector);
// 	      TmpValue += groundState[TmpPos] * groundState[TmpPos] * (Factorial.GetNumericalValue());	
// 	    }
// 	}
      
//       HermitianMatrix TmpDensityMatrix(1);
//       TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
//       return TmpDensityMatrix;
//     }


//   BosonOnTorusWithMagneticTranslationsShort TmpDestinationHilbertSpace(nbrBosonSector, this->MaxMomentum, kySector);
//   cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
//   HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
//   int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   long TmpNbrNonZeroElements = 0;
//   unsigned long* TmpMonomial2 = new unsigned long [nbrBosonSector];
//   unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
//   unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
//   BosonOnTorusWithMagneticTranslationsShort TmpHilbertSpace(ComplementaryNbrBosonSector, this->MaxMomentum, ComplementaryKySector);
//   FactorialCoefficient Factorial;

//   for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
//     {
//       int Pos = 0;
//       TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.StateKyMax[MinIndex] + TmpHilbertSpace.NbrBosons - 1, TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateKyMax);
//       TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex], TmpHilbertSpace.TemporaryStateKyMax + TmpHilbertSpace.NbrBosons - 1 , TmpMonomial1);
//       for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
// 	{
// 	  TmpDestinationHilbertSpace.FermionToBoson(TmpDestinationHilbertSpace.StateDescription[j], TmpDestinationHilbertSpace.StateKyMax[j] + TmpDestinationHilbertSpace.NbrBosons - 1, TmpDestinationHilbertSpace.TemporaryState, TmpDestinationHilbertSpace.TemporaryStateKyMax);
// 	  TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.StateDescription[j], TmpDestinationHilbertSpace.TemporaryStateKyMax + TmpDestinationHilbertSpace.NbrBosons - 1, TmpMonomial2);

// 	  int TmpIndex2 = 0;
// 	  int TmpIndex3 = 0;
// 	  int TmpIndex4 = 0;
// 	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < nbrBosonSector)) 
// 	    {
// 	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
// 		{
// 		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
// 		  ++TmpIndex2;
// 		  ++TmpIndex4;		  
// 		}
// 	      if (TmpIndex2 < ComplementaryNbrBosonSector)
// 		{
// 		  while ((TmpIndex3 < nbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
// 		    {
// 		      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
// 		      ++TmpIndex3;
// 		      ++TmpIndex4;		  
// 		    }
// 		}
// 	    }
// 	  while (TmpIndex2 < ComplementaryNbrBosonSector)
// 	    {
// 	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
// 	      ++TmpIndex2;
// 	      ++TmpIndex4;		  
// 	    }
// 	  while (TmpIndex3 < nbrBosonSector)
// 	    {
// 	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
// 	      ++TmpIndex3;
// 	      ++TmpIndex4;		  
// 	    }

// 	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
// 	  int TmpPos = this->FindStateIndex(TmpState, TmpMonomial3[0] + this->NbrBosons - 1);
// 	  if (TmpPos != this->HilbertSpaceDimension)
// 	    {
//  	      Factorial.SetToOne();
//  	      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateKyMax; ++k)
// 		{
// 		  if (TmpHilbertSpace.TemporaryState[k] > 1)
// 		    Factorial.FactorialDivide(TmpHilbertSpace.TemporaryState[k]);
// 		}
//   	      for (int k = 0; k <= TmpDestinationHilbertSpace.TemporaryStateKyMax; ++k)
// 		{
// 		  if (TmpDestinationHilbertSpace.TemporaryState[k] > 1)
// 		    Factorial.FactorialDivide(TmpDestinationHilbertSpace.TemporaryState[k]);
// 		}
//  	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateKyMax);
//  	      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
// 		if (this->TemporaryState[k] > 1)
//  		  Factorial.FactorialMultiply(this->TemporaryState[k]);
// 	      Factorial.BinomialDivide(this->NbrBosons, nbrBosonSector);
// 	      TmpStatePosition[Pos] = TmpPos;
// 	      TmpStatePosition2[Pos] = j;
// 	      TmpStateCoefficient[Pos] = sqrt(Factorial.GetNumericalValue());
// 	      ++Pos;
// 	    }
// 	}
//       if (Pos != 0)
// 	{
// 	  ++TmpNbrNonZeroElements;
// 	  for (int j = 0; j < Pos; ++j)
// 	    {
// 	      int Pos2 = TmpStatePosition2[j];
// 	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[j];
// 	      for (int k = 0; k < Pos; ++k)
// 		if (TmpStatePosition2[k] >= Pos2)
// 		  {
// 		    TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
// 		  }
// 	    }
// 	}
//     }
//   delete[] TmpStatePosition2;
//   delete[] TmpStatePosition;
//   delete[] TmpStateCoefficient;
//   delete[] TmpMonomial1;
//   delete[] TmpMonomial2;
//   delete[] TmpMonomial3;
//   if (TmpNbrNonZeroElements > 0)	
//     return TmpDensityMatrix;
//   else
  {
    HermitianMatrix TmpDensityMatrixZero;
    return TmpDensityMatrixZero;
  }
}


// convert a state defined in the Ky basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector BosonOnTorusWithMagneticTranslationsShort::ConvertToKxKyBasis(ComplexVector& state, ParticleOnTorus* space)  
{
  BosonOnTorusShort* TmpSpace = (BosonOnTorusShort*) space;
  ComplexVector TmpVector (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      int TmpMaxMomentum = this->FermionicMaxMomentum;
      while ((TmpState >> TmpMaxMomentum) == 0x0ul)
	--TmpMaxMomentum;
      int Pos = TmpSpace->FindStateIndex(TmpState, TmpMaxMomentum);
      if (Pos < TmpSpace->HilbertSpaceDimension)
	{
	  TmpVector[i] =  state[Pos] * sqrt((double) this->NbrStateInOrbit[i]);
	}
    }
  return TmpVector;
}

// convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector BosonOnTorusWithMagneticTranslationsShort::ConvertFromKxKyBasis(ComplexVector& state, ParticleOnTorus* space)
{
  BosonOnTorusShort* TmpSpace = (BosonOnTorusShort*) space;
  ComplexVector TmpVector (TmpSpace->LargeHilbertSpaceDimension, true);
  Complex* FourrierCoefficients = new Complex [this->MomentumModulo];
  for (int i = 0; i < this->MomentumModulo; ++i)
    FourrierCoefficients[i] = Phase (-2.0 * M_PI * ((double) (i * this->KxMomentum)) / ((double) this->MomentumModulo));
  for (long i = 0l; i < TmpSpace->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = TmpSpace->StateDescription[i];
      int NbrTranslation = 0;
      TmpState = this->FindCanonicalFormAndTestXMomentumConstraint(TmpState, NbrTranslation);
      if (NbrTranslation >= 0)
	{
	  int TmpMaxMomentum = this->FermionicMaxMomentum;
	  while ((TmpState >> TmpMaxMomentum) == 0x0ul)
	    --TmpMaxMomentum;
	  int Pos = this->FindStateIndex(TmpState, TmpMaxMomentum);
	  if (Pos < this->HilbertSpaceDimension)
	    {
	      TmpVector[i] =  state[Pos] * FourrierCoefficients[NbrTranslation] / sqrt((double) this->NbrStateInOrbit[Pos]);
	    }
	}
    }
  delete[] FourrierCoefficients;
  return TmpVector;
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

ComplexVector& BosonOnTorusWithMagneticTranslationsShort::CoreC4Rotation (ComplexVector& inputState, ParticleOnTorusWithMagneticTranslations* inputSpace, 
									  ComplexVector& outputState, int minIndex, int nbrIndices, bool clockwise)
{
  BosonOnTorusWithMagneticTranslationsShort* TmpInputSpace = (BosonOnTorusWithMagneticTranslationsShort*) inputSpace;
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
  double TmpLogCoefficient3 = this->NbrBosons * log((double) this->MaxMomentum);
  double PhaseFactor = 2.0 * M_PI / ((double) this->MaxMomentum);
  if (clockwise == true)
    PhaseFactor *= -1.0;
  ComplexMatrix DeterminantMatrix (this->NbrBosons, this->NbrBosons);
  ComplexMatrix PhaseMatrix (this->MaxMomentum, this->MaxMomentum);
  for (int k = 0; k < this->MaxMomentum; ++k)
    for (int l = 0; l < this->MaxMomentum; ++l)
      PhaseMatrix[k][l] = Phase(PhaseFactor * ((double) (k * l)));
  for (int i = minIndex ; i < LastIndex; ++i)
    {
      this->FermionToBoson(this->StateDescription[i], this->FermionicMaxMomentum, 
			   this->TemporaryState, this->TemporaryStateKyMax);
      this->ConvertToMonomial(this->StateDescription[i], this->FermionicMaxMomentum, TmpOutputMonomial);
      Tmp = 0.0;
      TmpLogCoefficient = 0.0;
      for (int k = 0; k <= this->TemporaryStateKyMax; ++k)
	TmpLogCoefficient += LogFactorialCoefficients[this->TemporaryState[k]];
      for (int j = 0; j < TmpInputSpace->HilbertSpaceDimension; ++j)
	{
	  this->FermionToBoson(TmpInputSpace->StateDescription[j], TmpInputSpace->FermionicMaxMomentum, 
			       TmpInputSpace->TemporaryState, TmpInputSpace->TemporaryStateKyMax);
	  TmpInputSpace->ConvertToMonomial(TmpInputSpace->StateDescription[j], TmpInputSpace->FermionicMaxMomentum, TmpInputMonomial);
	  for (int k = 0; k < this->NbrBosons; ++k)
	    for (int l = 0; l < this->NbrBosons; ++l)
	      DeterminantMatrix[k][l] = PhaseMatrix[(int) TmpInputMonomial[k]][(int) TmpOutputMonomial[l]];
	  TmpLogCoefficient2 = 0.0;
	  for (int k = 0; k <= TmpInputSpace->TemporaryStateKyMax; ++k)
	    {
	      TmpLogCoefficient2 += LogFactorialCoefficients[TmpInputSpace->TemporaryState[k]];
	    }
	  Tmp += inputState[j] * DeterminantMatrix.Permanent() * exp(0.5 * (-TmpLogCoefficient2 - TmpLogCoefficient - TmpLogCoefficient3)) * sqrt((double) TmpInputSpace->NbrStateInOrbit[j]);
	}
      outputState[i] = Tmp * sqrt((double) (this->NbrStateInOrbit[i]));      
    }
  delete[] TmpInputMonomial;
  delete[] TmpOutputMonomial;
  delete[] LogFactorialCoefficients;
  return outputState;
}

