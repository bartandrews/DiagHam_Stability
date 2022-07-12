////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of bosons on sphere using                     //
//                       the Lz <-> -Lz symmetry such that                    //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 10/10/2007                      //
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
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
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

BosonOnSphereSymmetricBasisShort::BosonOnSphereSymmetricBasisShort ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// lzMax = maximum Lz value reached by a boson

BosonOnSphereSymmetricBasisShort::BosonOnSphereSymmetricBasisShort (int nbrBosons, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->FermionBasis = 0;
  this->FermionSymmetricBasis = new FermionOnSphereSymmetricBasis(nbrBosons, lzMax + nbrBosons - 1);
  this->HilbertSpaceDimension = this->FermionSymmetricBasis->HilbertSpaceDimension;

  this->StateLzMax = new int [this->HilbertSpaceDimension];
  int TmpLzMax = this->LzMax + this->NbrBosons - 1;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->FermionSymmetricBasis->StateDescription[i];
      while ((TmpState & (0x1ul << TmpLzMax)) == 0x0l)
	--TmpLzMax;
      this->StateLzMax[i] = TmpLzMax;
    }
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->Flag.Initialize();
 
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;

}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSphereSymmetricBasisShort::BosonOnSphereSymmetricBasisShort(const BosonOnSphereSymmetricBasisShort& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = 0;
  this->FermionSymmetricBasis = (FermionOnSphereSymmetricBasis*) bosons.FermionSymmetricBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->StateLzMax =  bosons.StateLzMax;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
}

// destructor
//

BosonOnSphereSymmetricBasisShort::~BosonOnSphereSymmetricBasisShort ()
{
//  if (this->FermionSymmetricBasis != 0)
//    delete this->FermionSymmetricBasis;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    delete[] this->StateLzMax;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereSymmetricBasisShort& BosonOnSphereSymmetricBasisShort::operator = (const BosonOnSphereSymmetricBasisShort& bosons)
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
      delete[] this->StateLzMax;
    }
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateLzMax =  bosons.StateLzMax;
  this->Flag = bosons.Flag;
  this->FermionBasis = 0;
  this->FermionSymmetricBasis = (FermionOnSphereSymmetricBasis*) bosons.FermionSymmetricBasis->Clone();
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereSymmetricBasisShort::Clone()
{
  return new BosonOnSphereSymmetricBasisShort(*this);
}

// convert a given state from Lz-symmetric basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector BosonOnSphereSymmetricBasisShort::ConvertToNbodyBasis(RealVector& state, BosonOnSphereShort& nbodyBasis)
{
  return this->FermionSymmetricBasis->ConvertToNbodyBasis(state, *(nbodyBasis.FermionBasis));
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

int BosonOnSphereSymmetricBasisShort::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  unsigned long State = this->FermionSymmetricBasis->StateDescription[index];
  unsigned long Signature = State & FERMION_SPHERE_SYMMETRIC_BIT;
  State &= FERMION_SPHERE_SYMMETRIC_MASK;
  this->FermionToBoson(State, this->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  if ((n1 > this->TemporaryStateLzMax) || (n2 > this->TemporaryStateLzMax) || (this->TemporaryState[n1] == 0) || (this->TemporaryState[n2] == 0) || ((n1 == n2) && (this->TemporaryState[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  for (int i = this->TemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = 0x0ul;
  coefficient = this->TemporaryState[n2];
  --this->TemporaryState[n2];
  coefficient *= this->TemporaryState[n1];
  --this->TemporaryState[n1];
  ++this->TemporaryState[m2];
  coefficient *= this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  State = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
  State = this->FermionSymmetricBasis->GetSignedCanonicalState(State);
  this->TemporaryStateLzMax += this->NbrBosons;
  --this->TemporaryStateLzMax;
  while (((State & FERMION_SPHERE_SYMMETRIC_MASK) >> this->TemporaryStateLzMax) == 0)
    --this->TemporaryStateLzMax;
  int TmpIndex = this->FermionSymmetricBasis->FindStateIndex(State,  this->TemporaryStateLzMax);
  if ((State & FERMION_SPHERE_SYMMETRIC_BIT) != Signature)
    if (Signature != 0)
      coefficient *= M_SQRT2;
    else
      coefficient *= M_SQRT1_2;
  return TmpIndex;
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereSymmetricBasisShort::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  unsigned long State = this->FermionSymmetricBasis->StateDescription[index];
  unsigned long Signature = State & FERMION_SPHERE_SYMMETRIC_BIT;
  State &= FERMION_SPHERE_SYMMETRIC_MASK;
  this->FermionToBoson(State, this->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  --nbrIndices;
  for (int i = 0; i <= nbrIndices; ++i)
    {
      if ((n[i] > this->TemporaryStateLzMax) || (this->TemporaryState[n[i]] == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  coefficient = 1.0;
  for (int i = nbrIndices; i >= 0; --i)
    {
      if (this->TemporaryState[n[i]] == 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension; 	    
	}
      coefficient *= (double) this->TemporaryState[n[i]];
      --this->TemporaryState[n[i]];
    }
  for (int i = nbrIndices; i >= 0; --i)
    {
      ++this->TemporaryState[m[i]];
      coefficient *= (double) this->TemporaryState[m[i]];
    }
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  State = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
  State = this->FermionSymmetricBasis->GetSignedCanonicalState(State);
  this->TemporaryStateLzMax += this->NbrBosons - 1;
  while (((State & FERMION_SPHERE_SYMMETRIC_MASK) >> this->TemporaryStateLzMax) == 0)
    --this->TemporaryStateLzMax;
  int TmpIndex = this->FermionSymmetricBasis->FindStateIndex(State,  this->TemporaryStateLzMax);
  if ((State & FERMION_SPHERE_SYMMETRIC_BIT) != Signature)
    if (Signature != 0)
      coefficient *= M_SQRT2;
    else
      coefficient *= M_SQRT1_2;
  return TmpIndex;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnSphereSymmetricBasisShort::AA (int index, int n1, int n2)
{
  unsigned long State = this->FermionSymmetricBasis->StateDescription[index];
  this->ProdASignature = State & FERMION_SPHERE_SYMMETRIC_BIT;
  State &= FERMION_SPHERE_SYMMETRIC_MASK;
  this->FermionToBoson(State, this->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  if ((n1 > this->ProdATemporaryStateLzMax) || (n2 > this->ProdATemporaryStateLzMax) || 
      (this->ProdATemporaryState[n1] == 0) || (this->ProdATemporaryState[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryState[n1] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryState[n2];
  --this->ProdATemporaryState[n2];
  Coefficient *= this->ProdATemporaryState[n1];
  --this->ProdATemporaryState[n1];
  for (int i = this->ProdATemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0ul;
  return sqrt(Coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereSymmetricBasisShort::ProdA (int index, int* n, int nbrIndices)
{
  unsigned long State = this->FermionSymmetricBasis->StateDescription[index];
  this->ProdASignature = State & FERMION_SPHERE_SYMMETRIC_BIT;
  State &= FERMION_SPHERE_SYMMETRIC_MASK;
  this->FermionToBoson(State, this->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  int TmpCoefficient = 1;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (n[i] > this->ProdATemporaryStateLzMax)
	return 0.0;
      unsigned long& Tmp = this->ProdATemporaryState[n[i]];
      if (Tmp == 0)
	return 0.0;
      TmpCoefficient *= Tmp;
      --Tmp;
    }
  for (int i = this->ProdATemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt((double) TmpCoefficient);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereSymmetricBasisShort::AdAd (int m1, int m2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m2];
  coefficient = this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  unsigned long State = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
  State = this->FermionSymmetricBasis->GetSignedCanonicalState(State);
  this->TemporaryStateLzMax += this->NbrBosons;
  --this->TemporaryStateLzMax;
  while (((State & FERMION_SPHERE_SYMMETRIC_MASK) >> this->TemporaryStateLzMax) == 0)
    --this->TemporaryStateLzMax;
  int TmpIndex = this->FermionSymmetricBasis->FindStateIndex(State, this->TemporaryStateLzMax);
  if ((State & FERMION_SPHERE_SYMMETRIC_BIT) != this->ProdASignature)
    if (this->ProdASignature != 0)
      coefficient *= M_SQRT2;
    else
      coefficient *= M_SQRT1_2;
  return TmpIndex;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereSymmetricBasisShort::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  int TmpCoefficient = 1;
  for (int i = 0; i < nbrIndices; ++i)
    TmpCoefficient *= ++this->TemporaryState[m[i]];
  coefficient = sqrt((double) TmpCoefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  unsigned long State = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
  State = this->FermionSymmetricBasis->GetSignedCanonicalState(State);
  this->TemporaryStateLzMax += this->NbrBosons - 1;
  while (((State & FERMION_SPHERE_SYMMETRIC_MASK) >> this->TemporaryStateLzMax) == 0)
    --this->TemporaryStateLzMax;
  int TmpIndex = this->FermionSymmetricBasis->FindStateIndex(State,  this->TemporaryStateLzMax);
  if ((State & FERMION_SPHERE_SYMMETRIC_BIT) != this->ProdASignature)
    if (this->ProdASignature != 0)
      coefficient *= M_SQRT2;
    else
      coefficient *= M_SQRT1_2;
  return TmpIndex;
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereSymmetricBasisShort::AdA (int index, int m)
{
  this->FermionToBoson(this->FermionSymmetricBasis->StateDescription[index] & FERMION_SPHERE_SYMMETRIC_MASK, this->StateLzMax[index], 
		       this->TemporaryState, this->TemporaryStateLzMax);
  if (this->TemporaryStateLzMax < m)  
    return 0.0;
  return (double) (this->TemporaryState[m]);  
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereSymmetricBasisShort::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->FermionSymmetricBasis->StateDescription[state] & FERMION_SPHERE_SYMMETRIC_MASK, this->StateLzMax[state], 
		       this->TemporaryState, this->TemporaryStateLzMax);
  int i = 0;
  for (; i <= this->TemporaryStateLzMax; ++i)
    Str << this->TemporaryState[i] << " ";
  for (; i <= this->LzMax; ++i)
    Str << "0 ";
  Str << "   lzmax = " << this->TemporaryStateLzMax;
  return Str;
}


