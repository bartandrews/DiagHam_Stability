////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system size such that            //
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
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereFullShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h" 
#include "MathTools/ClebschGordanCoefficients.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/Permutations.h"

#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereSymmetrizeU1U1StateOperation.h"

#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <map>

//#define VERBOSE_FF


using std::cout;
using std::endl;
using std::dec;
using std::hex;
using std::map;
using std::pair;

// default constructor
//

BosonOnSphereShort::BosonOnSphereShort ()
{
  this->Minors = 0;
  this->KeptCoordinates = 0;
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson

BosonOnSphereShort::BosonOnSphereShort (int nbrBosons, int totalLz, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
//   cout << this->NbrBosons << " " << this->TotalLz << " " << this->LzMax << endl;
  this->NbrLzValue = this->LzMax + 1;
  if (nbrBosons > 0)
    this->FermionBasis = new FermionOnSphere(nbrBosons, totalLz, lzMax + nbrBosons - 1);
  else
    this->FermionBasis = new FermionOnSphere(nbrBosons, totalLz, lzMax);
  this->HilbertSpaceDimension = this->FermionBasis->HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = this->FermionBasis->LargeHilbertSpaceDimension;

  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->Flag.Initialize();
  int TmpLzMax = this->LzMax;
  this->TargetSpace = this;

  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;

}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSphereShort::BosonOnSphereShort(const BosonOnSphereShort& bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
}

// destructor
//

BosonOnSphereShort::~BosonOnSphereShort ()
{
  if (this->FermionBasis != 0)
    delete this->FermionBasis;
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
      if (this->KeptCoordinates != 0)
	delete this->KeptCoordinates;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereShort& BosonOnSphereShort::operator = (const BosonOnSphereShort& bosons)
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
    }
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereShort::Clone()
{
  return new BosonOnSphereShort(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnSphereShort::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->TargetSpace = (BosonOnSphereShort*) targetSpace;
  this->FermionBasis->SetTargetSpace(this->TargetSpace->FermionBasis);
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int BosonOnSphereShort::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereShort::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereShort::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereShort::ExtractSubspace (AbstractQuantumNumber& q, 
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

int BosonOnSphereShort::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  if ((n1 > this->TemporaryStateLzMax) || (n2 > this->TemporaryStateLzMax) || (this->TemporaryState[n1] == 0) || (this->TemporaryState[n2] == 0) || ((n1 == n2) && (this->TemporaryState[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  for (int i = this->TemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = 0ul;
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
  return this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  --nbrIndices;
  for (int i = 0; i <= nbrIndices; ++i)
    {
      if ((n[i] > this->TemporaryStateLzMax) || (this->TemporaryState[n[i]] == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  for (int i = this->TemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = 0ul;
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
  return this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnSphereShort::AA (int index, int n1, int n2)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
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

double BosonOnSphereShort::ProdA (int index, int* n, int nbrIndices)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
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


// apply Prod_i a_mi operator to the state produced using ProdA method (without destroying it)
// use double when calculating normalization factors to avoid overflow
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereShort::ProdAL (int index, int* n, int nbrIndices)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  double TmpCoefficient = 1.0;
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
  return sqrt(TmpCoefficient);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::AdAd (int m1, int m2, double& coefficient)
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
  return this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::ProdAd (int* m, int nbrIndices, double& coefficient)
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
  return this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
// use double when calculating normalization factors to avoid overflow
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::ProdAdL (int* m, int nbrIndices, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  double TmpCoefficient = 1.0;
  for (int i = 0; i < nbrIndices; ++i)
    TmpCoefficient *= ++this->TemporaryState[m[i]];
  coefficient = sqrt(TmpCoefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereShort::AdA (int index, int m)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  if (this->TemporaryStateLzMax < m)  
    return 0.0;
  return (double) (this->TemporaryState[m]);  
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereShort::AdA (long index, int m)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  if (this->TemporaryStateLzMax < m)  
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

int BosonOnSphereShort::AdA (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);

  if ((this->TemporaryStateLzMax < n)  || (this->TemporaryState[n] == 0))
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState[n];
  --this->TemporaryState[n];
  if ((this->TemporaryStateLzMax == n) && (this->TemporaryState[n] == 0))
    {
      while ((this->TemporaryState[this->TemporaryStateLzMax] == 0)&&( this->TemporaryStateLzMax > 0 ))
	--this->TemporaryStateLzMax;
    }
  if (this->TemporaryStateLzMax < m) 
    {
      for (int i = this->TemporaryStateLzMax + 1; i <= m; ++i)
	this->TemporaryState[i] = 0;
      this->TemporaryStateLzMax = m;
    }
  ++this->TemporaryState[m];
  coefficient *= (double) this->TemporaryState[m];
  coefficient = sqrt(coefficient);  
  int TmpIndex = this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
  if (TmpIndex == this->TargetSpace->HilbertSpaceDimension)
    coefficient = 0.0;
  return TmpIndex;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// n = index of the annihilation operator
// return value = double where the multiplicative factor has to be stored

double BosonOnSphereShort::A (int index, int n)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  for (int i = this->ProdATemporaryStateLzMax + 1; i <= this->LzMax ;++i)
    this->ProdATemporaryState[i] = 0;
  
  double coefficient;
  if ((this->ProdATemporaryStateLzMax < n)  || (this->ProdATemporaryState[n] == 0))
    { 
      coefficient = 0.0;
      return coefficient;      
    }
  coefficient = (double) this->ProdATemporaryState[n];
  --this->ProdATemporaryState[n];
  if ((this->ProdATemporaryStateLzMax == n) && (this->ProdATemporaryState[n] == 0))
    {
      while ((this->ProdATemporaryState[this->ProdATemporaryStateLzMax] == 0)&&( this->ProdATemporaryStateLzMax > 0 ))
	--this->ProdATemporaryStateLzMax;
    }
  return coefficient;

}


// apply a^+_m operator to the state produced using the A or Ad method (without destroying it)
//
// m = first index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be store
// return value = index of the destination state 
int BosonOnSphereShort::Ad (int m, double& coefficient)
{
  this->TemporaryStateLzMax = this->ProdATemporaryStateLzMax;
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  if (this->TemporaryStateLzMax < m) 
    {
      for (int i = this->TemporaryStateLzMax + 1; i <= m; ++i)
	this->TemporaryState[i] = 0;
      this->TemporaryStateLzMax = m;
    }
  ++this->TemporaryState[m];
  coefficient *= (double) this->TemporaryState[m];
  coefficient = sqrt(coefficient);  
  int TmpIndex = this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
  if (TmpIndex == this->HilbertSpaceDimension)
    coefficient = 0.0;
  return TmpIndex;
}

// apply a_n  operator to a given state. 
//
// index = index of the state on which the operator has to be applied
// n = index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value =  index of the resulting state 

int BosonOnSphereShort::A (int index, int n, double& coefficient)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  for (int i = this->TemporaryStateLzMax + 1; i <= this->LzMax ;++i)
    this->TemporaryState[i] = 0;
  
  if ((this->TemporaryStateLzMax < n)  || (this->TemporaryState[n] == 0))
    { 
      coefficient = 0.0;
      return coefficient;      
    }
  coefficient = sqrt((double) this->TemporaryState[n]);
  --this->TemporaryState[n];
  if ((this->TemporaryStateLzMax == n) && (this->TemporaryState[n] == 0))
    {
      while ((this->TemporaryState[this->TemporaryStateLzMax] == 0)&&( this->TemporaryStateLzMax > 0 ))
	--this->TemporaryStateLzMax;
    }
  return this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), 
							 this->TemporaryStateLzMax + this->TargetSpace->NbrBosons - 1);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereShort::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state], this->TemporaryState, this->TemporaryStateLzMax);
  int i = 0;
  for (; i <= this->TemporaryStateLzMax; ++i)
    Str << this->TemporaryState[i] << " ";
  for (; i <= this->LzMax; ++i)
    Str << "0 ";
  return Str;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereShort::PrintStateMonomial (ostream& Str, long state)
{
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  this->ConvertToMonomial(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state], TmpMonomial);
  Str << "[";
  if (TmpMonomial[0] != 0)
    Str << TmpMonomial[0];
  for (int i = 1; (i < this->NbrBosons) && (TmpMonomial[i] > 0); ++i)
    Str << "," << TmpMonomial[i];
  Str << "]";
  delete[] TmpMonomial;
  return Str;
}

// print a given State using the monomial notation, with one column per particle (using space as a seperator)
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereShort::PrintColumnFormattedStateMonomial (ostream& Str, long state)
{
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  this->ConvertToMonomial(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state], TmpMonomial);
  if (TmpMonomial[0] != 0)
    Str << TmpMonomial[0];
  for (int i = 1; (i < this->NbrBosons) && (TmpMonomial[i] > 0); ++i)
    Str << " " << TmpMonomial[i];
  delete[] TmpMonomial;
  return Str;
}


// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int BosonOnSphereShort::FindStateIndex(char* stateDescription)
{
   char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != (this->LzMax + 1))
    return -1;
  int TmpNbrParticles = 0;
  int TmpTotalLz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      int Tmp = atoi(TmpDescription[i]);
      this->TemporaryState[i] = Tmp;
      TmpTotalLz += (i * Tmp);
      TmpNbrParticles += Tmp;
      delete[] TmpDescription[i];
    }
  delete[] TmpDescription;
  if ((TmpNbrParticles != this->NbrBosons) || (TmpTotalLz != ((this->TotalLz + this->NbrBosons * this->LzMax) >> 1)))
    return -1;
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
 return this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, NewLzMax), NewLzMax + this->NbrBosons - 1);
}

// find state index from an array of occupied orbitals
//
// stateDescription = array describing the state (stored as k1,k2,k3,...)
// return value = corresponding index, -1 if an error occured
int BosonOnSphereShort::FindStateIndex(int* stateDescription)
{
  for (int i = 0; i <= this->LzMax; ++i)
    this->TemporaryState[i] = 0l;
  for (int i = 0; i < this->NbrBosons; ++i)
    ++this->TemporaryState[stateDescription[i]];
  unsigned long TmpState = this->BosonToFermion(this->TemporaryState, this->LzMax);
  int TmpLzMax = this->FermionBasis->LzMax;
  while ((TmpState >> TmpLzMax) == 0x0ul)
    --TmpLzMax;
  return this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
}

// find state index from unsigned long representation
//
// stateDescription = integer describint the state
// lzMax = the lzmax of the state
// return value = corresponding index, -1 if an error occured

int BosonOnSphereShort::FindStateIndex(unsigned long int stateDescription, int lzMax)
{
  return this->FermionBasis->FindStateIndex(stateDescription, lzMax);
}

// find state index from its occupation number description
//
// stateDescription = array that describes the state in the occupation number basis
// return value = corresponding index, -1 if an error occured

int BosonOnSphereShort::FindStateIndexFromOccupationNumber(unsigned long* stateDescription)
{
  int NewLzMax =  this->LzMax;
  while ((NewLzMax > 0) && (stateDescription[NewLzMax] == 0x0ul))
    --NewLzMax;
  return this->FermionBasis->FindStateIndex(this->BosonToFermion(stateDescription, NewLzMax), NewLzMax + this->NbrBosons - 1);  
}

// evaluate wave function in real space using a given basis and only for a given range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex BosonOnSphereShort::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
						  int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
  ComplexMatrix Perm(this->NbrBosons, this->NbrBosons);
  ComplexMatrix Functions(this->LzMax + 1, this->NbrBosons);
  int NbrComponentPerParticles = position.GetVectorDimension() / this->NbrBosons;
  RealVector TmpCoordinates(NbrComponentPerParticles);
  int* Indices = new int [this->NbrBosons];
  int Pos;
  int Lz;
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      for (int i = 0; i < NbrComponentPerParticles; ++i)
	{
	  TmpCoordinates[i] = position[(j * NbrComponentPerParticles) + i];
	}
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
      this->FermionToBoson(this->FermionBasis->StateDescription[k], this->FermionBasis->StateLzMax[k], this->TemporaryState, this->TemporaryStateLzMax);
      while (Pos < this->NbrBosons)
	{
	  TmpStateDescription = this->TemporaryState[Lz];
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

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  BosonOnSphereShort::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState)
{  
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrBosonSector == 0))
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
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrBosonSector == this->NbrBosons))
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

  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  if ((ShiftedLzComplementarySector < (NbrBosonsComplementarySector * subsytemSize)) || (ShiftedLzComplementarySector > (NbrBosonsComplementarySector * (this->LzMax))))
    {
      RealSymmetricMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }

  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0.0;
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  unsigned long  TmpState2 = 0x0;
	  for (int i = 0; i < nbrBosonSector; ++i)
	    TmpState2 |= 0x1ul << i;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | TmpState2;
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		TmpValue += groundState[TmpPos] * groundState[TmpPos];	
	    }

	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}      
    }
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0;
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		TmpValue += groundState[TmpPos] * groundState[TmpPos];	
	    }	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int MinIndex = 0;
  // int MaxIndex = this->HilbertSpaceDimension - 1;
  if (nbrBosonSector == 1)
    {
      double TmpValue = 0.0;
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | (0x1ul << ShiftedLzSector);
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    TmpValue += groundState[TmpPos] * groundState[TmpPos];	
	}
      RealSymmetricMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  if (NbrBosonsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
      BosonOnSphere TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      double TmpValue;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpValue = groundState[MinIndex + i];
	  for (int j = i; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    TmpDensityMatrix.SetMatrixElement(i, j, TmpValue * groundState[MinIndex + j]);
	}
      return TmpDensityMatrix;
    }


  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
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

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem

RealMatrix BosonOnSphereShort::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState)
{  
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrBosonSector == 0))
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
  
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrBosonSector == this->NbrBosons))
	{
	  RealMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension,1,true);
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
  
  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  if ((ShiftedLzComplementarySector < (NbrBosonsComplementarySector * subsytemSize)) || (ShiftedLzComplementarySector > (NbrBosonsComplementarySector * (this->LzMax))))
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;	  
    }
  long TmpNbrNonZeroElements = 0;
  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  unsigned long  TmpState2 = 0x0;
	  RealMatrix TmpEntanglementMatrix(1,TmpHilbertSpace.HilbertSpaceDimension,true);
	  
	  for (int i = 0; i < nbrBosonSector; ++i)
	    TmpState2 |= 0x1ul << i;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | TmpState2;
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  ++TmpNbrNonZeroElements;
		  TmpEntanglementMatrix.AddToMatrixElement(0,MinIndex,groundState[TmpPos]);
		}
	    }
	  if (TmpNbrNonZeroElements == 0)
	    {
	      RealMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons, 2 * ShiftedLzComplementarySector - ( this->NbrBosons  * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  RealMatrix TmpEntanglementMatrix(1,TmpHilbertSpace.HilbertSpaceDimension,true);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  ++TmpNbrNonZeroElements;
		  TmpEntanglementMatrix.AddToMatrixElement(0,MinIndex,groundState[TmpPos]);
		}
	    }
	  if (TmpNbrNonZeroElements == 0)
	    {
	      RealMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
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
  if (nbrBosonSector == 1)
    {
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
      RealMatrix TmpEntanglementMatrix(1,TmpHilbertSpace.HilbertSpaceDimension,true);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | (0x1ul << ShiftedLzSector);
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      ++TmpNbrNonZeroElements;
	      TmpEntanglementMatrix.AddToMatrixElement(0,MinIndex,groundState[TmpPos]);
	    }
	}
      if (TmpNbrNonZeroElements == 0)
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
      return TmpEntanglementMatrix;
    }
  
  if (NbrBosonsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
      BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,1, true);
      MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpEntanglementMatrix.AddToMatrixElement(i,0,groundState[MinIndex + i]);
	}
      return TmpEntanglementMatrix;
    }
  
  
  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  
  
  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,TmpHilbertSpace.HilbertSpaceDimension, true);
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      ++TmpNbrNonZeroElements;
	      TmpEntanglementMatrix.AddToMatrixElement(j,MinIndex,groundState[TmpPos]);
	    }
	}
    }
  
  if (TmpNbrNonZeroElements > 0)	
    return TmpEntanglementMatrix;
  else
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem

LongRationalMatrix BosonOnSphereShort::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrBosonSector, int lzSector, LongRationalVector& groundState)
{  
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrBosonSector == 0))
	{
	  LongRationalMatrix TmpEntanglementMatrix(1,1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1l);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  LongRationalMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrBosonSector == this->NbrBosons))
	{
	  LongRationalMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension,1,true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[i]);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  LongRationalMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  if ((ShiftedLzComplementarySector < (NbrBosonsComplementarySector * subsytemSize)) || (ShiftedLzComplementarySector > (NbrBosonsComplementarySector * (this->LzMax))))
    {
      LongRationalMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;	  
    }
  long TmpNbrNonZeroElements = 0;
  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  unsigned long  TmpState2 = 0x0;
	  LongRationalMatrix TmpEntanglementMatrix(1,TmpHilbertSpace.HilbertSpaceDimension,true);
	  
	  for (int i = 0; i < nbrBosonSector; ++i)
	    TmpState2 |= 0x1ul << i;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | TmpState2;
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  ++TmpNbrNonZeroElements;
		  TmpEntanglementMatrix.AddToMatrixElement(0,MinIndex,groundState[TmpPos]);
		}
	    }
	  if (TmpNbrNonZeroElements == 0)
	    {
	      LongRationalMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  LongRationalMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons, 2 * ShiftedLzComplementarySector - ( this->NbrBosons  * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  LongRationalMatrix TmpEntanglementMatrix(1,TmpHilbertSpace.HilbertSpaceDimension,true);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  ++TmpNbrNonZeroElements;
		  TmpEntanglementMatrix.AddToMatrixElement(0,MinIndex,groundState[TmpPos]);
		}
	    }
	  if (TmpNbrNonZeroElements == 0)
	    {
	      LongRationalMatrix TmpEntanglementMatrix;
	      return TmpEntanglementMatrix;
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  LongRationalMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }
  
  int MinIndex = 0;
  if (nbrBosonSector == 1)
    {
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
      LongRationalMatrix TmpEntanglementMatrix(1,TmpHilbertSpace.HilbertSpaceDimension,true);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | (0x1ul << ShiftedLzSector);
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      ++TmpNbrNonZeroElements;
	      TmpEntanglementMatrix.AddToMatrixElement(0,MinIndex,groundState[TmpPos]);
	    }
	}
      if (TmpNbrNonZeroElements == 0)
	{
	  LongRationalMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
      return TmpEntanglementMatrix;
    }
  
  if (NbrBosonsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
	  LongRationalMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
      BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      LongRationalMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,1, true);
      MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpEntanglementMatrix.AddToMatrixElement(i,0,groundState[MinIndex + i]);
	}
      return TmpEntanglementMatrix;
    }
  
  
  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  
  
  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
  LongRationalMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension,TmpHilbertSpace.HilbertSpaceDimension, true);
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      ++TmpNbrNonZeroElements;
	      TmpEntanglementMatrix.AddToMatrixElement(j,MinIndex,groundState[TmpPos]);
	    }
	}
    }
  
  if (TmpNbrNonZeroElements > 0)	
    return TmpEntanglementMatrix;
  else
    {
      LongRationalMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
}

// reconstruct a state that contains only a certain subset of Schmidt eigenvalues of the given ground state
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// eigenvalueCut = only keep Schmidt levels that are larger than e^{-eigenvalueCut}
// groundState = reference on the total system ground state
// rebuiltSchmitGroundState = reference on the final state
// diagonalizedDensityMatrix = set of density matrix (Schmidt) eigenvalues
// transformationMatrix = the "truncation" matrix that connects the coefficients  (in the N-body basis) of the ground state and the final (truncated) state
// return value = reconstructed ground state vector

RealVector& BosonOnSphereShort::EvaluatePartialSchmidtDecomposition (int subsytemSize, int nbrBosonSector, int lzSector, double eigenvalueCut, RealVector& groundState, RealVector& rebuiltSchmidtGroundState, RealDiagonalMatrix& diagonalizedDensityMatrix, RealMatrix& transformationMatrix)
{
  if (subsytemSize <= 0)
    {
      cout << "Requested zero orbitals in the subsystem... exiting" << endl;
      return rebuiltSchmidtGroundState;
    }
  if (subsytemSize > this->LzMax)
    {
      cout << "Requested the entire system... exiting" << endl;
      for (long i = 0l; i < groundState.GetLargeVectorDimension(); ++i)
	rebuiltSchmidtGroundState[i] = groundState[i];
      return rebuiltSchmidtGroundState;
    }

  //Evaluate how many eigenvalues to keep 
  int NbrKeptEigenvalues = 0;  
  for (int i = 0; i < diagonalizedDensityMatrix.GetNbrRow(); ++i)
    if (diagonalizedDensityMatrix[i] >=  eigenvalueCut)
      ++NbrKeptEigenvalues;
  cout << "Keeping "<<NbrKeptEigenvalues<<" / "<<diagonalizedDensityMatrix.GetNbrRow()<<" eigenvalues"<<endl;
  if (NbrKeptEigenvalues == 0)
    return rebuiltSchmidtGroundState;
 
 
  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  if ((ShiftedLzComplementarySector < (NbrBosonsComplementarySector * subsytemSize)) || (ShiftedLzComplementarySector > (NbrBosonsComplementarySector * (this->LzMax))))
    {
      return rebuiltSchmidtGroundState;
    }

  //Only one orbital is kept in the subsystem, lzSector has to be 0
  //and reduced density matrix has only one element
  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  unsigned long  TmpState2 = 0x0;
	  for (int i = 0; i < nbrBosonSector; ++i)
	    TmpState2 |= 0x1ul << i;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | TmpState2;
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		  rebuiltSchmidtGroundState[TmpPos] += groundState[TmpPos];
	    }
         return rebuiltSchmidtGroundState;
	}
      else
	{
           return rebuiltSchmidtGroundState;
	}      
    }

  //Subsystem contains 0 bosons, lzSector can therefore only be 0
  //and the reduced density matrix contains only one element
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		  rebuiltSchmidtGroundState[TmpPos] += groundState[TmpPos];	
	    }	  
         return rebuiltSchmidtGroundState;
	}
      else
	{
         return rebuiltSchmidtGroundState;
	}
    }

  //Case when there is only 1 boson in the subsystem A
  //Reduced density matrix only has one element
  if (nbrBosonSector == 1)
    {
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | (0x1ul << ShiftedLzSector);
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
   	      rebuiltSchmidtGroundState[TmpPos] += groundState[TmpPos];	
	}
     return rebuiltSchmidtGroundState;
    }

  //Subsystem A contains all the bosons
  if (NbrBosonsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
         return rebuiltSchmidtGroundState;
	}
      BosonOnSphere TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
      cout << " N_A=N, subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;

      //Build the "truncation" matrix
      RealMatrix TmpMatrix (TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension);
      for (int i = 0; i < TmpMatrix.GetNbrRow(); ++i)
        for (int j = 0; j < TmpMatrix.GetNbrRow(); ++j)
          {
            double Tmp = 0.0;
	    for (int k = 0; k < NbrKeptEigenvalues; ++k)
	      Tmp += transformationMatrix(i, k) * transformationMatrix(j, k);
	    TmpMatrix(i, j) = Tmp;
          }

      int MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  double Tmp = 0.0;
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
              Tmp += groundState[MinIndex + j] * TmpMatrix(i, j);
   	  rebuiltSchmidtGroundState[MinIndex + i] = Tmp;	
	}
      return rebuiltSchmidtGroundState;
    }


  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "N_A<N, subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  RealMatrix TmpMatrix (TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension);
  for (int i = 0; i < TmpMatrix.GetNbrRow(); ++i)
    for (int j = 0; j < TmpMatrix.GetNbrRow(); ++j)
      {
	double Tmp = 0.0;
	for (int k = 0; k < NbrKeptEigenvalues; ++k)
	  Tmp += transformationMatrix(i, k) * transformationMatrix(j, k);
	TmpMatrix(i, j) = Tmp;
      }


  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
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
              double Tmp = 0.0;
	      for (int k = 0; k < Pos; ++k)
                    Tmp += groundState[TmpStatePosition[k]] * TmpMatrix(Pos2, TmpStatePosition2[k]); 
              rebuiltSchmidtGroundState[TmpStatePosition[j]] = Tmp;
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  return rebuiltSchmidtGroundState;  
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem

HermitianMatrix BosonOnSphereShort::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, ComplexVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrBosonSector == 0))
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
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrBosonSector == this->NbrBosons))
	{
	  HermitianMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, Conj(groundState[i]) * groundState[j]);
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  if ((ShiftedLzComplementarySector < (NbrBosonsComplementarySector * subsytemSize)) || (ShiftedLzComplementarySector > (NbrBosonsComplementarySector * (this->LzMax))))
    {
      HermitianMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }

  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  Complex TmpValue = 0.0;
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  unsigned long  TmpState2 = 0x0;
	  for (int i = 0; i < nbrBosonSector; ++i)
	    TmpState2 |= 0x1ul << i;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | TmpState2;
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		TmpValue += Conj(groundState[TmpPos]) * groundState[TmpPos];	
	    }

	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}      
    }
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  Complex TmpValue = 0;
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		TmpValue += Conj(groundState[TmpPos]) * groundState[TmpPos];	
	    }	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int MinIndex = 0;
  // int MaxIndex = this->HilbertSpaceDimension - 1;
  if (nbrBosonSector == 1)
    {
      Complex TmpValue = 0.0;
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | (0x1ul << ShiftedLzSector);
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    TmpValue += Conj(groundState[TmpPos]) * groundState[TmpPos];	
	}
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  if (NbrBosonsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
      BosonOnSphere TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      Complex TmpValue;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpValue = Conj(groundState[MinIndex + i]);
	  for (int j = i; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    TmpDensityMatrix.SetMatrixElement(i, j, TmpValue * groundState[MinIndex + j]);
	}
      return TmpDensityMatrix;
    }


  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
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
	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]);
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
		  }
	    }
	}
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

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particle. The geometrical cut is a stripe.
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax+shitedCut to -Lzmax+shitedCut+subsytemSize-1)
// shiftedCut = first orbital belonging to the subsystem (with angular momentum -Lzmax+shitedCut)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnSphereShort::EvaluateShiftedPartialDensityMatrix (int subsytemSize, int nbrShiftedOrbitals, int nbrBosonSector, int lzSector, RealVector& groundState)
{
  if (nbrShiftedOrbitals == 0)
    return this->EvaluatePartialDensityMatrix(subsytemSize, nbrBosonSector, lzSector, groundState);
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrBosonSector == 0))
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
  if (subsytemSize > this->LzMax)
    {
      RealSymmetricMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }

  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = ((lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1) + (nbrBosonSector * nbrShiftedOrbitals);
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  if ((ShiftedLzComplementarySector < 0) || (ShiftedLzComplementarySector > (NbrBosonsComplementarySector * (this->LzMax))))
    {
      RealSymmetricMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }

  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  for (int NbrBosonsLeft = 0; NbrBosonsLeft <= NbrBosonsComplementarySector; ++NbrBosonsLeft)
    {
      int MaxShiftedLzLeft = NbrBosonsLeft * (nbrShiftedOrbitals - 1);
      for (int ShiftedLzLeft = 0; ShiftedLzLeft <= MaxShiftedLzLeft; ++ShiftedLzLeft)
	{
	  int NbrBosonsRight = NbrBosonsComplementarySector - NbrBosonsLeft;
	  int ShiftedLzRight = ShiftedLzComplementarySector - ShiftedLzLeft;
	  if ((ShiftedLzRight >= (NbrBosonsRight * (nbrShiftedOrbitals + subsytemSize))) && (ShiftedLzRight <= (NbrBosonsRight * this->LzMax)))
	    {
	      BosonOnSphereShort TmpHilbertSpaceLeft(NbrBosonsLeft, 2 * ShiftedLzLeft - (NbrBosonsLeft * (nbrShiftedOrbitals - 1)), 
						     nbrShiftedOrbitals - 1);
	      BosonOnSphereShort TmpHilbertSpaceRight(NbrBosonsRight, 2 * (ShiftedLzRight - (NbrBosonsRight * (nbrShiftedOrbitals + subsytemSize))) - (NbrBosonsRight * (this->LzMax - nbrShiftedOrbitals - subsytemSize)), this->LzMax - nbrShiftedOrbitals - subsytemSize);	      
	      for (int MinIndexLeft = 0; MinIndexLeft < TmpHilbertSpaceLeft.HilbertSpaceDimension; ++MinIndexLeft)    
		{
		  unsigned long TmpComplementaryStateLeft = TmpHilbertSpaceLeft.FermionBasis->StateDescription[MinIndexLeft];
		  for (int MinIndexRight = 0; MinIndexRight < TmpHilbertSpaceRight.HilbertSpaceDimension; ++MinIndexRight)    
		    {
		      int Pos = 0;
		      unsigned long TmpComplementaryState = (TmpHilbertSpaceRight.FermionBasis->StateDescription[MinIndexRight] << (subsytemSize + nbrBosonSector + nbrShiftedOrbitals + NbrBosonsLeft)) | TmpComplementaryStateLeft;
		      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
			{
			  unsigned long TmpState = ((TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] << (nbrShiftedOrbitals + NbrBosonsLeft))
						    | TmpComplementaryState);
			  int TmpLzMax = this->FermionBasis->LzMax;
			  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
			    --TmpLzMax;
			  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
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
		}
	    }
	}
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


// Compute the row and column dimension of the orbital entanglement matrix of 2 cuts counting only those rows/columns that are not completely zero
// Also returns from the set of indices in the reduced density matrix corresponding to rows with atleast one non-zero entry
// Columns contain the hilbert space of B and C which are traced out
// SizeB = number of orbitals in part B, i.e. in the cap around Lzmax/2.
// SizeA = number of orbitals in the bulk of the sphere 
// NbrBosonsA = number of particles that belong to A
// groundState = reference on the total system ground state
// LzA = Lz sector of A in which the density matrix has to be evaluated as measured on a sphere with only A
// return value = pointer with the 1st element being the row dimension
//                2nd element is the column dimension of the oem 
//                3rd element onward gives the positions of the rows in the oem/reduced density matrix which are not filled with 0's
//                (returns 0 if there is a probem/there is no hilbert space)

long* BosonOnSphereShort::Compute2CutEntanglementMatrixDimensions (int SizeB, int SizeA, int NbrBosonsA, int LzA, RealVector& groundState)
{
  //Lz of full ground state on a disk
  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  //Lz of A on a disk
  int ShiftedLzA = (LzA + NbrBosonsA * (SizeA-1) ) >> 1;
  
  
  long DimA;
  int SizeC = this->LzMax +1 -SizeA -SizeB;
  
  BosonOnSphereShort TmpHilbertSpaceA(NbrBosonsA, LzA, SizeA - 1);
  
  if(NbrBosonsA>0)
    {
      DimA = TmpHilbertSpaceA.HilbertSpaceDimension;
    }
  else 
    if(NbrBosonsA==0)
      DimA=1;
    else
      {
	cout<<"The number of bosons in A is negative!";
	return 0; 
      }
  int* TmpIsSqueezedStateA = new int [DimA];
  //Initialize it
  
  long* Dims = new long [2+DimA];
  long PosninDims = 2;
  Dims[0] =0;
  Dims[1] =0;
  for(long i=0; i<DimA; i++)
    {
      TmpIsSqueezedStateA[i] = 0;
      Dims[2+i] =0;
    }
  for(int NbrBosonsB=0; NbrBosonsB <= (this->NbrBosons - NbrBosonsA); NbrBosonsB++)
    {
      int ShiftedLzBMax= (SizeB-1)*NbrBosonsB;
      int NbrBosonsC = this->NbrBosons - NbrBosonsA - NbrBosonsB;
      int ShiftedLzCMin = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - ShiftedLzBMax - (SizeA+SizeB)*NbrBosonsC;
      if(ShiftedLzCMin > (SizeC-1)*NbrBosonsC)
	continue;
      int ShiftedLzCMax = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - (SizeA+SizeB)*NbrBosonsC;
      if(ShiftedLzCMax<0)
	continue;
      int ShiftedLzBTemp = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - (this->LzMax)*NbrBosonsC;
      int ShiftedLzBMin = (ShiftedLzBTemp > 0) ? ShiftedLzBTemp : 0;
      int ShiftedLzBPossMax = (ShiftedLzCMax < ShiftedLzBMax) ? ShiftedLzCMax : ShiftedLzBMax;
      for(int ShiftedLzB=ShiftedLzBMin; ShiftedLzB <= ShiftedLzBPossMax; ShiftedLzB++)
	{
	  //2*Lz of B when only B is on a sphere
	  long DimB, DimC;
	  int SphereLzB = (ShiftedLzB<<1) - ( NbrBosonsB* (SizeB-1) );
	  
	  //Lz of C on a disk  
	  int ShiftedLzC = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - ShiftedLzB - (SizeA+SizeB)*NbrBosonsC;
	  
	  if( (ShiftedLzC<0) || (ShiftedLzC > (SizeC-1)*NbrBosonsC) )
	    {
	      //cout<<"No part C on a disk possible for values of NbrBosonsB="<<NbrBosonsB<<" and ShiftedLzB="<<ShiftedLzB<<endl;
	      continue;
	    }
	  
	  BosonOnSphereShort *TmpHilbertSpaceB=0;
	  if(NbrBosonsB>0)
	    {
	      TmpHilbertSpaceB= new BosonOnSphereShort (NbrBosonsB, SphereLzB, SizeB - 1 );
	      DimB = TmpHilbertSpaceB->GetHilbertSpaceDimension();
	    }
	  else 
	    if(NbrBosonsB==0)
	      DimB=1;
	    else
	      {
		cout<<"The number of bosons in B is negative!";
		return 0; 
	      }
	  
	  
	  //2*Lz of C when only C is on a sphere  
	  int SphereLzC = (ShiftedLzC<<1) - ( NbrBosonsC* (SizeC-1) );
	  BosonOnSphereShort *TmpHilbertSpaceC = 0;  
	  
	  if(NbrBosonsC>0)
	    {
	      TmpHilbertSpaceC = new BosonOnSphereShort (NbrBosonsC, SphereLzC, SizeC - 1 );
	      DimC = TmpHilbertSpaceC->GetHilbertSpaceDimension();
	    }
	  else 
	    if(NbrBosonsC==0)
	      DimC=1;
	    else
	      {
		cout<<"The number of bosons in C is negative!";
		return 0; 
	      }
	  for ( long i=0; i<DimB ; ++i)
	    {
	      for ( long j=0; j<DimC ; ++j)
		{
		  //Shifting the C Hilbert space element left
		  unsigned long TmpStateC, TmpStateB;
		  if(NbrBosonsC>0)
		    TmpStateC= TmpHilbertSpaceC->FermionBasis->StateDescription[j] << (SizeA + SizeB + NbrBosonsA + NbrBosonsB);
		  else
		    TmpStateC = 0;
		  
		  if(NbrBosonsB>0)	
		    TmpStateB = TmpHilbertSpaceB->FermionBasis->StateDescription[i];
		  else
		    TmpStateB = 0;
		  int Pos = 0;
		  for ( long k=0; k<DimA ; ++k)
		    {
		      unsigned long TmpStateA;
		      if(NbrBosonsA>0)	
			TmpStateA = TmpHilbertSpaceA.FermionBasis->StateDescription[k] << (SizeB + NbrBosonsB);
		      else
			TmpStateA = 0;
		      
		      unsigned long TmpFullState = (TmpStateB | TmpStateA) | TmpStateC;
		      
		      int TmpLzMax = this->FermionBasis->LzMax;
		      while (((TmpFullState >> TmpLzMax) & 0x1ul) == 0x0ul)
			--TmpLzMax;
		      int TmpPos = this->FermionBasis->FindStateIndex(TmpFullState, TmpLzMax);
		      
		      if ((TmpPos != this->HilbertSpaceDimension) && ( groundState[TmpPos] != 0) )
			{
			  if(Pos ==0)
			    {
			      Dims[1]++;
			      Pos++;
			    }
			  if(TmpIsSqueezedStateA[k] == 0)
			    {  
			      TmpIsSqueezedStateA[k] = 1;
			      Dims[PosninDims] = k;
			      PosninDims++;
			      Dims[0]++;
			    }
			}
		    }
		  
		}
            }
	  
	  //Finishing the braces of the big loop running over possible NbrParticlesB and ShiftedLzB		
	  if(NbrBosonsC>0)
	    delete TmpHilbertSpaceC;
          if(NbrBosonsB>0)
	    delete TmpHilbertSpaceB;
	  
	}
    }
  
  delete [] TmpIsSqueezedStateA; 
  return Dims;            
}
  			
  			
// evaluate a density matrix with 2 cuts of the whole system described by the RealVector groundState. The reduced density matrix is evaluated for a given Lz sector and number of particles
//
// SizeB = number of orbitals in part B, i.e. in the cap around Lzmax/2.
// SizeA = number of orbitals in the bulk of the sphere 
// NbrBosonsA = number of particles that belong to A
// groundState = reference on the total system ground state
// LzA = Lz sector of A in which the density matrix has to be evaluated as measured on a sphere with only A
// return value = density matrix of the subsytem (return a zero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnSphereShort::Evaluate2CutReducedDensityMatrix (int SizeB, int SizeA, int NbrBosonsA, int LzA, RealVector& groundState)
{
  //Lz of full ground state on a disk
  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  //Lz of A on a disk
  int ShiftedLzA = (LzA + NbrBosonsA * (SizeA-1) ) >> 1;
  //Make a Hilbert Space for A of dimension DimA only if there are some bosons in it
  long DimA;
  int SizeC = this->LzMax +1 -SizeA -SizeB;
  
  BosonOnSphereShort TmpHilbertSpaceA(NbrBosonsA, LzA, SizeA - 1);
  
  if(NbrBosonsA > 0)
    {
      DimA = TmpHilbertSpaceA.HilbertSpaceDimension;
    }
  else 
    if(NbrBosonsA==0)
      DimA=1;
    else
      {
	cout<<"The number of bosons in A is negative!";
	RealSymmetricMatrix TmpDensityMatrixZero;
	return TmpDensityMatrixZero; 
      }
  RealSymmetricMatrix TmpDensityMatrix(DimA, true);
  
  int* TmpFullStatePosition = new int [DimA];
  int* TmpStateAPosition = new int [DimA];
  int TmpNonZeroElements = 0;
  
  for(int NbrBosonsB=0; NbrBosonsB <= (this->NbrBosons - NbrBosonsA); NbrBosonsB++)
    {
      int ShiftedLzBMax= (SizeB-1)*NbrBosonsB;
      int NbrBosonsC = this->NbrBosons - NbrBosonsA - NbrBosonsB;
      int ShiftedLzCMin = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - ShiftedLzBMax - (SizeA+SizeB)*NbrBosonsC;
      if(ShiftedLzCMin > (SizeC-1)*NbrBosonsC)
	continue;
      int ShiftedLzCMax = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - (SizeA+SizeB)*NbrBosonsC;
      if(ShiftedLzCMax<0)
	continue;
      
      for(int ShiftedLzB=0; ShiftedLzB <= ShiftedLzBMax; ShiftedLzB++)
	{
	  //2*Lz of B when only B is on a sphere
	  long DimB, DimC;
	  int SphereLzB = (ShiftedLzB<<1) - ( NbrBosonsB* (SizeB-1) );
	  
	  //Lz of C on a disk  
	  int ShiftedLzC = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - ShiftedLzB - (SizeA+SizeB)*NbrBosonsC;
	  
	  if( (ShiftedLzC<0) || (ShiftedLzC > (SizeC-1)*NbrBosonsC) )
	    {
	      cout<<"No part C on a disk possible for values of NbrBosonsB="<<NbrBosonsB<<" and ShiftedLzB="<<ShiftedLzB<<endl;
	      continue;
	    }
	  BosonOnSphereShort *TmpHilbertSpaceB=0;
	  if(NbrBosonsB>0) 
	    {
	      TmpHilbertSpaceB= new BosonOnSphereShort (NbrBosonsB, SphereLzB, SizeB - 1 );
	      DimB = TmpHilbertSpaceB->GetHilbertSpaceDimension();
	    }
	  else 
	    if(NbrBosonsB==0)
	      DimB=1;
	    else
	      {
		cout<<"The number of bosons in B is negative!";
		RealSymmetricMatrix TmpDensityMatrixZero;
		return TmpDensityMatrixZero; 
	      }
	  
	  //2*Lz of C when only C is on a sphere  
	  int SphereLzC = (ShiftedLzC<<1) - ( NbrBosonsC* (SizeC-1) );
	  BosonOnSphereShort *TmpHilbertSpaceC = 0;  
	  
	  if(NbrBosonsC>0)
	    {
	      TmpHilbertSpaceC = new BosonOnSphereShort (NbrBosonsC, SphereLzC, SizeC - 1 );
	      DimC = TmpHilbertSpaceC->GetHilbertSpaceDimension();
	    }
	  else 
	    if(NbrBosonsC==0)
	      DimC=1;
	    else
	      {
		cout<<"The number of bosons in C is negative!";
		RealSymmetricMatrix TmpDensityMatrixZero;
		return TmpDensityMatrixZero; 
	      }
	  for (long i = 0; i < DimB ; ++i)
	    {
	      for ( long j = 0; j < DimC ; ++j)
		{
		  //Shifting the C Hilbert space element left
		  unsigned long TmpStateC, TmpStateB;
		  if(NbrBosonsC>0)
		    TmpStateC= TmpHilbertSpaceC->FermionBasis->StateDescription[j] << (SizeA + SizeB + NbrBosonsA + NbrBosonsB);
		  else
		    TmpStateC = 0;
		  
		  if(NbrBosonsB>0)	
		    TmpStateB = TmpHilbertSpaceB->FermionBasis->StateDescription[i];
		  else
		    TmpStateB = 0;
		  int Pos=0;
		  
		  for ( long k=0; k<DimA ; ++k)
		    {
		      unsigned long TmpStateA;
		      if(NbrBosonsA>0)	
			TmpStateA = TmpHilbertSpaceA.FermionBasis->StateDescription[k] << (SizeB + NbrBosonsB);
		      else
			TmpStateA = 0;
		      
		      unsigned long TmpFullState = (TmpStateB | TmpStateA) | TmpStateC;
		      
		      int TmpLzMax = this->FermionBasis->LzMax;
		      while (((TmpFullState >> TmpLzMax) & 0x1ul) == 0x0ul)
			--TmpLzMax;
		      int TmpPos = this->FermionBasis->FindStateIndex(TmpFullState, TmpLzMax);
		      
		      if (TmpPos != this->HilbertSpaceDimension)
			{
			  //this->PrintState(cout,TmpPos)<<" "<<groundState[TmpPos]<<" " << k<<endl;
			  TmpFullStatePosition[Pos] = TmpPos;
			  //cout<<TmpPos<<endl;
			  TmpStateAPosition[Pos] = k;
			  ++Pos;
			}
		    }
		  
		  if(Pos!=0)
		    {
		      //Then for the given B and C, we have some contribution to the reduced density matrix
		      ++TmpNonZeroElements;
		      for(long k=0; k<Pos; k++)
			{
			  long RowA=TmpStateAPosition[k];
			  double coeff_k=groundState[TmpFullStatePosition[k]];
			  for(long m=0; m<Pos; m++) 
			    if (TmpStateAPosition[m] >= RowA)
			      TmpDensityMatrix.AddToMatrixElement(RowA,TmpStateAPosition[m],coeff_k*groundState[TmpFullStatePosition[m]]);
			  //cout<<"trying to set something"<<coeff_k;
			}
		    }
		  
		}
	    }
	  //Finishing the braces of the big loop running over possible NbrParticlesB and ShiftedLzB		
  		  if(NbrBosonsC>0)
		    delete TmpHilbertSpaceC;
		  if(NbrBosonsB>0)
		    delete TmpHilbertSpaceB;
		  
	}
    }
  
  delete[] TmpStateAPosition;
  delete[] TmpFullStatePosition;
  if (TmpNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  BosonOnSphereShort::EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState, AbstractArchitecture* architecture)
{  
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
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
      if (lzSector == this->TotalLz)
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

  if (nbrBosonSector == 1)
    {
      double TmpValue = 0.0;
      unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
      unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - 1, this->TotalLz - lzSector, this->LzMax);
      unsigned long ShiftedLzVSector = (lzSector + this->LzMax) >> 1;
      FactorialCoefficient Factorial;
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
	  TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateLzMax);
	  int TmpIndex2 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial1[TmpIndex2] >= ShiftedLzVSector))
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  TmpMonomial3[TmpIndex4] = ShiftedLzVSector;
	  ++TmpIndex4;
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      Factorial.SetToOne();
	      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateLzMax; ++k)
		if (TmpHilbertSpace.TemporaryState[k] > 1)
		  Factorial.FactorialDivide(TmpHilbertSpace.TemporaryState[k]);
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
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

  if (abs(this->TotalLz - lzSector) > (ComplementaryNbrBosonSector * LzMax))
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, this->LzMax);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosonSector, this->TotalLz - lzSector, this->LzMax);

  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &TmpDestinationHilbertSpace, &TmpHilbertSpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
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

long BosonOnSphereShort::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace, RealVector& groundState, RealSymmetricMatrix* densityMatrix)
{
   BosonOnSphereShort* TmpHilbertSpace =  (BosonOnSphereShort*) complementaryHilbertSpace;
   BosonOnSphereShort* TmpDestinationHilbertSpace =  (BosonOnSphereShort*) destinationHilbertSpace;
   int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
   int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
   unsigned long* TmpMonomial2 = new unsigned long [NbrBosonSector];
   unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
   unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
   int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   int MaxIndex = minIndex + nbrIndex;
   long TmpNbrNonZeroElements = 0l;

   double* LogFactorials = new double[this->NbrBosons + 1];
   LogFactorials[0] = 0.0;
   LogFactorials[1] = 0.0;
   for (int i = 2 ; i <= this->NbrBosons; ++i)
     LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
   double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[NbrBosonSector];

   for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->FermionBasis->StateDescription[i], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[i], TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateLzMax);
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateLzMax; ++k)
	TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }
   for (; minIndex < MaxIndex; ++minIndex)    
     {
      int Pos = 0;
      TmpHilbertSpace->ConvertToMonomial(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpMonomial1);
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateLzMax);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace->TemporaryStateLzMax; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  TmpDestinationHilbertSpace->ConvertToMonomial(TmpDestinationHilbertSpace->FermionBasis->StateDescription[j], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[j], TmpMonomial2);
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < NbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonSector)
		{
		  while ((TmpIndex3 < NbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
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
	  while (TmpIndex3 < NbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
	      TmpFactorial *= 0.5; 
	      
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      TmpStateCoefficient[Pos] = exp(TmpFactorial);
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
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		    ++TmpNbrNonZeroElements;
		  }
	    }
	}
     }
   delete[] TmpMonomial2;
   delete[] TmpMonomial1;
   delete[] TmpMonomial3;
   delete[] TmpStatePosition;
   delete[] TmpStatePosition2;
   delete[] TmpStateCoefficient;
   delete[] TmpDestinationLogFactorials;
   return TmpNbrNonZeroElements;
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

long BosonOnSphereShort::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace, ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
   BosonOnSphereShort* TmpHilbertSpace =  (BosonOnSphereShort*) complementaryHilbertSpace;
   BosonOnSphereShort* TmpDestinationHilbertSpace =  (BosonOnSphereShort*) destinationHilbertSpace;
   int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
   int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
   unsigned long* TmpMonomial2 = new unsigned long [NbrBosonSector];
   unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
   unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
   int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   int MaxIndex = minIndex + nbrIndex;
   long TmpNbrNonZeroElements = 0l;

   double* LogFactorials = new double[this->NbrBosons + 1];
   LogFactorials[0] = 0.0;
   LogFactorials[1] = 0.0;
   for (int i = 2 ; i <= this->NbrBosons; ++i)
     LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
   double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[NbrBosonSector];

   for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->FermionBasis->StateDescription[i], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[i], TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateLzMax);
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateLzMax; ++k)
	TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }
   for (; minIndex < MaxIndex; ++minIndex)    
     {
      int Pos = 0;
      TmpHilbertSpace->ConvertToMonomial(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpMonomial1);
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateLzMax);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace->TemporaryStateLzMax; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  TmpDestinationHilbertSpace->ConvertToMonomial(TmpDestinationHilbertSpace->FermionBasis->StateDescription[j], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[j], TmpMonomial2);
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < NbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonSector)
		{
		  while ((TmpIndex3 < NbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
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
	  while (TmpIndex3 < NbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
	      TmpFactorial *= 0.5; 
	      
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      TmpStateCoefficient[Pos] = exp(TmpFactorial);
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]) * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		    ++TmpNbrNonZeroElements;
		  }
	    }
	}
     }
   delete[] TmpMonomial2;
   delete[] TmpMonomial1;
   delete[] TmpMonomial3;
   delete[] TmpStatePosition;
   delete[] TmpStatePosition2;
   delete[] TmpStateCoefficient;
   delete[] TmpDestinationLogFactorials;
   return TmpNbrNonZeroElements;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using real space partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// thetaTop =  inclination angle defining one edge of the cut in degrees
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnSphereShort::EvaluatePartialDensityMatrixRealSpacePartition (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, RealVector& groundState, AbstractArchitecture* architecture)
{
  if ((abs(lzSector) > (nbrBosonSector * this->LzMax)) || (thetaBottom <= thetaTop) || (phiRange <= 0.0))
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

  thetaTop *= M_PI / 180.0;
  thetaBottom *= M_PI / 180.0;
  phiRange /= 360.0;
  
  double* IncompleteBetaThetaTop = 0;
  double* IncompleteBetaThetaBottom = 0;

  this->EvaluatePartialDensityMatrixRealSpacePartitionCoefficient(this->LzMax, thetaTop, thetaBottom, IncompleteBetaThetaTop, IncompleteBetaThetaBottom);

  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  unsigned long* TmpMonomial1 = new unsigned long [this->NbrBosons];
	  double TmpValue = 0.0;
	  for (int MinIndex = 0; MinIndex < this->HilbertSpaceDimension; ++MinIndex)    
	    {
	      this->ConvertToMonomial(this->FermionBasis->StateDescription[MinIndex], this->FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
	      double FormFactor = 0.0;
	      for (int i=0; i < this->NbrBosons; i++)
		FormFactor += log(1.0 - IncompleteBetaThetaBottom[TmpMonomial1[i]] + IncompleteBetaThetaTop[TmpMonomial1[i]] + (1.0 - phiRange) * (IncompleteBetaThetaBottom[TmpMonomial1[i]] - IncompleteBetaThetaTop[TmpMonomial1[i]]) );
	      FormFactor = exp(FormFactor);
	      TmpValue += groundState[MinIndex] * groundState[MinIndex] * FormFactor;	
	    }
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue); 
	  
	  delete[] IncompleteBetaThetaTop;
	  delete[] IncompleteBetaThetaBottom;
	  delete[] TmpMonomial1;
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
      if (lzSector == this->TotalLz)
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension, true);
	  unsigned long* TmpMonomial1 = new unsigned long [this->NbrBosons];
	  double* TmpStateCoefficient = new double [this->HilbertSpaceDimension];
	  for( int i = 0; i < this->HilbertSpaceDimension; i++)
	    {
	      TmpStateCoefficient[i] = 0.5 * this->NbrBosons * log(phiRange);
	      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial1);
	      for( int j=0; j<this->NbrBosons; j++)
		{
		  TmpStateCoefficient [i] += 0.5*log( IncompleteBetaThetaBottom[TmpMonomial1[j]] - IncompleteBetaThetaTop[TmpMonomial1[j]]);
		}
	      TmpStateCoefficient[i] = exp(TmpStateCoefficient[i]);
	    }
	  
	  for(int pos1 = 0; pos1 < this->HilbertSpaceDimension; pos1++)
	    for(int pos2 = pos1; pos2 < this->HilbertSpaceDimension; pos2++)
	      {
		TmpDensityMatrix.SetMatrixElement(pos1, pos2, groundState[pos1]*groundState[pos2]*TmpStateCoefficient[pos1]*TmpStateCoefficient[pos2]);
	      }
	  delete[] TmpMonomial1;
	  delete[] TmpStateCoefficient;
	  delete[] IncompleteBetaThetaTop;
	  delete[] IncompleteBetaThetaBottom;
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;

  if (nbrBosonSector == 1)
    {
      double TmpValue = 0.0;
      unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
      unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - 1, this->TotalLz - lzSector, this->LzMax);
      unsigned long ShiftedLzVSector = (lzSector + this->LzMax) >> 1;
      double TmpStateCoefficient = phiRange * (IncompleteBetaThetaBottom[ShiftedLzVSector] - IncompleteBetaThetaTop[ShiftedLzVSector]);
      FactorialCoefficient Factorial;
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
	  TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateLzMax);
	  double FormFactor = 0.0;
	  for (int i = 0; i < ComplementaryNbrBosonSector; i++)
	    FormFactor += log(1.0 - IncompleteBetaThetaBottom[TmpMonomial1[i]] + IncompleteBetaThetaTop[TmpMonomial1[i]] + (1.0 - phiRange) * (IncompleteBetaThetaBottom[TmpMonomial1[i]] - IncompleteBetaThetaTop[TmpMonomial1[i]]));
	  FormFactor = exp(FormFactor);
	  int TmpIndex2 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial1[TmpIndex2] >= ShiftedLzVSector))
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  TmpMonomial3[TmpIndex4] = ShiftedLzVSector;
	  ++TmpIndex4;
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      Factorial.SetToOne();
	      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateLzMax; ++k)
		if (TmpHilbertSpace.TemporaryState[k] > 1)
		  Factorial.FactorialDivide(TmpHilbertSpace.TemporaryState[k]);
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		if (this->TemporaryState[k] > 1)
		  Factorial.FactorialMultiply(this->TemporaryState[k]);
	      TmpValue += groundState[TmpPos] * groundState[TmpPos] * (Factorial.GetNumericalValue()) * FormFactor * TmpStateCoefficient;	
	    }
	}
      
      RealSymmetricMatrix TmpDensityMatrix(1);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
      delete[] IncompleteBetaThetaTop;
      delete[] IncompleteBetaThetaBottom;
      delete [] TmpMonomial1;
      delete [] TmpMonomial3;
      return TmpDensityMatrix;
    }

  if (abs(this->TotalLz - lzSector) > (ComplementaryNbrBosonSector * LzMax))
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, this->LzMax);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosonSector, this->TotalLz - lzSector, this->LzMax);

   FQHESphereParticleEntanglementSpectrumOperation Operation(this, &TmpDestinationHilbertSpace, &TmpHilbertSpace, groundState, TmpDensityMatrix, IncompleteBetaThetaBottom, IncompleteBetaThetaTop, phiRange);
   Operation.ApplyOperation(architecture);
   if (Operation.GetNbrNonZeroMatrixElements() > 0)	
     return TmpDensityMatrix;
   else
     {
       RealSymmetricMatrix TmpDensityMatrixZero;
       return TmpDensityMatrixZero;
     }
   

//   int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   double* TmpStateCoefficientOccFactorial = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   long TmpNbrNonZeroElements = 0;
//   unsigned long* TmpMonomial2 = 0;
//   unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
//   unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];

//   double* LogFactorials = new double[this->NbrBosons + 1];
//   LogFactorials[0] = 0.0;
//   LogFactorials[1] = 0.0;
//   for (int i = 2 ; i <= this->NbrBosons; ++i)
//     LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 

//   FactorialCoefficient Factorial;
//   unsigned long** TmpDestinationHilbertSpaceOccupationNumbers = new unsigned long* [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   unsigned long** TmpDestinationHilbertSpaceMonomial = new unsigned long* [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   int* TmpDestinationHilbertSpaceLzMax = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
//     {
//       TmpDestinationHilbertSpaceOccupationNumbers[i] = new unsigned long [this->NbrLzValue];
//       TmpDestinationHilbertSpaceMonomial[i] = new unsigned long [nbrBosonSector];
//       TmpDestinationHilbertSpace.FermionToBoson(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpDestinationHilbertSpaceOccupationNumbers[i], TmpDestinationHilbertSpaceLzMax[i]);
//       TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpDestinationHilbertSpaceMonomial[i]);
//       unsigned long* TmpOccupationNumber = TmpDestinationHilbertSpaceOccupationNumbers[i];
//       int TmpLzMax = TmpDestinationHilbertSpaceLzMax[i];
//       double TmpFactor = 0.0;
//       for (int k = 0; k <= TmpLzMax; ++k)
// 	TmpFactor += LogFactorials[TmpOccupationNumber[k]];
//       TmpDestinationLogFactorials[i] =  TmpFactor;
//     }

//   //Compute the coefficients multiplying rhoA in TmpStateCoefficient
//   for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
//     {
//       TmpStateCoefficient [i] = 0.5 * nbrBosonSector * log(phiRange);
      
//       for( int j=0; j<nbrBosonSector; j++)
// 	{
// 	  TmpStateCoefficient [i] += 0.5 * log( IncompleteBetaThetaBottom[TmpDestinationHilbertSpaceMonomial[i][j]] - IncompleteBetaThetaTop[TmpDestinationHilbertSpaceMonomial[i][j]]);
// 	}
//       TmpStateCoefficient[i] = exp(TmpStateCoefficient[i]);
//     }

//   for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
//     {
//       int Pos = 0;
//       TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
//       double FormFactor = 0.0;
//       for (int i=0; i < ComplementaryNbrBosonSector; i++)
// 	FormFactor += log(1.0 - IncompleteBetaThetaBottom[TmpMonomial1[i]] + IncompleteBetaThetaTop[TmpMonomial1[i]] + (1.0 - phiRange) * (IncompleteBetaThetaBottom[TmpMonomial1[i]] - IncompleteBetaThetaTop[TmpMonomial1[i]]) );
//       FormFactor = exp(FormFactor);
//       TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateLzMax);
//       double TmpHilbertSpaceFactorial = 0.0;
//       for (int k = 0; k <= TmpHilbertSpace.TemporaryStateLzMax; ++k)
// 	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace.TemporaryState[k]];
//       for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
// 	{
// 	  TmpMonomial2 = TmpDestinationHilbertSpaceMonomial[j];
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
// 	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
// 	  if (TmpPos != this->HilbertSpaceDimension)
// 	    {	      
// 	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
// 	      double TmpFactorial = 0.0;	      
// 	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
// 		TmpFactorial += LogFactorials[this->TemporaryState[k]];
// 	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j];
// 	      TmpFactorial *= 0.5; 
	      
// 	      TmpStatePosition[Pos] = TmpPos;
// 	      TmpStatePosition2[Pos] = j;
// 	      TmpStateCoefficientOccFactorial[Pos] = exp(TmpFactorial);
// 	      ++Pos;
// 	    }
// 	}
//       if (Pos != 0)
// 	{
// 	  ++TmpNbrNonZeroElements;
// 	  for (int j = 0; j < Pos; ++j)
// 	    {
// 	      int Pos2 = TmpStatePosition2[j];
// 	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[Pos2] * TmpStateCoefficientOccFactorial[j];
// 	      for (int k = 0; k < Pos; ++k)
// 		if (TmpStatePosition2[k] >= Pos2)
// 		  {
// 		    TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k],  FormFactor * TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[TmpStatePosition2[k]] * TmpStateCoefficientOccFactorial[k]);
// 		  }
// 	    }
// 	}
//     }
//   for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
//     {
//       delete[] TmpDestinationHilbertSpaceOccupationNumbers[i];
//       delete[] TmpDestinationHilbertSpaceMonomial[i];
//     }
//   delete[] TmpDestinationHilbertSpaceOccupationNumbers;
//   delete[] TmpDestinationHilbertSpaceLzMax;
//   delete[] TmpDestinationHilbertSpaceMonomial;
//   delete[] TmpStatePosition2;
//   delete[] TmpStatePosition;
//   delete[] TmpStateCoefficient;
//   delete[] TmpMonomial1;
//   delete[] TmpMonomial3;
//   if (TmpNbrNonZeroElements > 0)	
//     return TmpDensityMatrix;
//   else
//     {
//       RealSymmetricMatrix TmpDensityMatrixZero;
//       return TmpDensityMatrixZero;
//     }
}

// core part of the evaluation density matrix real space partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// incompleteBetaThetaTop = pointer to the array where the top part coefficients are stored
// incompleteBetaThetaBotton = pointer on the pointer to the array where the bottom part coefficients are stored
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// return value = number of components that have been added to the density matrix

long BosonOnSphereShort::EvaluatePartialDensityMatrixRealSpacePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
									     RealVector& groundState,  RealSymmetricMatrix* densityMatrix, double* incompleteBetaThetaBottom, double* incompleteBetaThetaTop, double phiRange)
{
  BosonOnSphereShort* TmpHilbertSpace =  (BosonOnSphereShort*) complementaryHilbertSpace;
  BosonOnSphereShort* TmpDestinationHilbertSpace =  (BosonOnSphereShort*) destinationHilbertSpace;
  int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
  int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
  unsigned long* TmpMonomial2 = new unsigned long [NbrBosonSector];
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpStateCoefficientOccFactorial = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;

  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 


  for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->FermionBasis->StateDescription[i], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[i], TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateLzMax);
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateLzMax; ++k)
	TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;

      TmpStateCoefficient[i] = 0.5 * NbrBosonSector * log(phiRange);
      
      TmpDestinationHilbertSpace->ConvertToMonomial(TmpDestinationHilbertSpace->FermionBasis->StateDescription[i], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[i], TmpMonomial2);
     for(int j = 0; j < NbrBosonSector; j++)
	{
	  TmpStateCoefficient[i] += 0.5 * log( incompleteBetaThetaBottom[TmpMonomial2[j]] - incompleteBetaThetaTop[TmpMonomial2[j]]);
	}
      TmpStateCoefficient[i] = exp(TmpStateCoefficient[i]);
    }

  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace->ConvertToMonomial(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpMonomial1);
      double FormFactor = 0.0;
      for (int i=0; i < ComplementaryNbrBosonSector; i++)
	FormFactor += log(1.0 - incompleteBetaThetaBottom[TmpMonomial1[i]] + incompleteBetaThetaTop[TmpMonomial1[i]] + (1.0 - phiRange) * (incompleteBetaThetaBottom[TmpMonomial1[i]] - incompleteBetaThetaTop[TmpMonomial1[i]]) );
      FormFactor = exp(FormFactor);
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateLzMax);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace->TemporaryStateLzMax; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  TmpDestinationHilbertSpace->ConvertToMonomial(TmpDestinationHilbertSpace->FermionBasis->StateDescription[j], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[j], TmpMonomial2);
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < NbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonSector)
		{
		  while ((TmpIndex3 < NbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
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
	  while (TmpIndex3 < NbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {	      
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j];
	      TmpFactorial *= 0.5; 
	      
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      TmpStateCoefficientOccFactorial[Pos] = exp(TmpFactorial);
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[Pos2] * TmpStateCoefficientOccFactorial[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k],  FormFactor * TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[TmpStatePosition2[k]] * TmpStateCoefficientOccFactorial[k]);
		  }
	    }
	}
    }
  delete[] TmpMonomial2;
  delete[] TmpMonomial1;
  delete[] TmpMonomial3;
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  delete[] TmpDestinationLogFactorials;
  delete[] TmpStateCoefficientOccFactorial;
  return TmpNbrNonZeroElements;
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix BosonOnSphereShort::EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState, bool removeBinomialCoefficient)
{	
  if ( abs(lzSector) > (nbrBosonSector * this->LzMax) )
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
  
  if (nbrBosonSector == 0)
    {
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons, this->TotalLz - lzSector, this->LzMax);
      RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  int TmpLzMax = this->LzMax + this->NbrBosons - 1;
	  unsigned long TmpState = this->FermionBasis->StateDescription[i];
	  while ((TmpState >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  
	  TmpEntanglementMatrix.SetMatrixElement(0, TmpHilbertSpace.FermionBasis->FindStateIndex(TmpState, TmpLzMax), groundState[i]);
	}
      return TmpEntanglementMatrix;
    }
  
  
  if (nbrBosonSector == this->NbrBosons)
    {
      if (lzSector == this->TotalLz)
	{
	  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, this->LzMax);
	  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, 1,true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      int TmpLzMax = this->LzMax + this->NbrBosons - 1;
	      unsigned long TmpState = this->FermionBasis->StateDescription[i];
	      while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      
	      TmpEntanglementMatrix.SetMatrixElement(TmpDestinationHilbertSpace.FermionBasis->FindStateIndex(TmpState, TmpLzMax), 0, groundState[i]);
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }

  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;
  
  if ( abs(this->TotalLz - lzSector) > (ComplementaryNbrBosonSector * this->LzMax))
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
  
  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, this->LzMax);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  long TmpNbrNonZeroElements = 0;
  unsigned long* TmpMonomial2 = 0;
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];

  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
  double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[nbrBosonSector];
  if (removeBinomialCoefficient == true)
    TmpLogBinomial = 0.0;

  BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosonSector, this->TotalLz - lzSector, this->LzMax);
  FactorialCoefficient Factorial;
  unsigned long** TmpDestinationHilbertSpaceOccupationNumbers = new unsigned long* [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  unsigned long** TmpDestinationHilbertSpaceMonomial = new unsigned long* [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpDestinationHilbertSpaceLzMax = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpHilbertSpace.HilbertSpaceDimension, true);
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpaceOccupationNumbers[i] = new unsigned long [this->NbrLzValue];
      TmpDestinationHilbertSpaceMonomial[i] = new unsigned long [nbrBosonSector];
      TmpDestinationHilbertSpace.FermionToBoson(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpDestinationHilbertSpaceOccupationNumbers[i], TmpDestinationHilbertSpaceLzMax[i]);
      TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpDestinationHilbertSpaceMonomial[i]);
      unsigned long* TmpOccupationNumber = TmpDestinationHilbertSpaceOccupationNumbers[i];
      int TmpLzMax = TmpDestinationHilbertSpaceLzMax[i];
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpLzMax; ++k)
	TmpFactor += LogFactorials[TmpOccupationNumber[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
      TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateLzMax);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateLzMax; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace.TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  TmpMonomial2 = TmpDestinationHilbertSpaceMonomial[j];
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
	  while ( TmpIndex2 < ComplementaryNbrBosonSector)
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
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
	      TmpFactorial *= 0.5; 	      
	      ++TmpNbrNonZeroElements;
	      double Tmp = exp(TmpFactorial) * groundState[TmpPos];
	      TmpEntanglementMatrix.SetMatrixElement(j, MinIndex, Tmp);
	    }
	}
    }
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      delete[] TmpDestinationHilbertSpaceOccupationNumbers[i];
      delete[] TmpDestinationHilbertSpaceMonomial[i];
    }
  delete[] TmpDestinationHilbertSpaceOccupationNumbers;
  delete[] TmpDestinationHilbertSpaceLzMax;
  delete[] TmpDestinationHilbertSpaceMonomial;
  delete[] TmpMonomial1;
  delete[] TmpMonomial3;
  if (TmpNbrNonZeroElements > 0)
    return TmpEntanglementMatrix;
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// The entanglement matrix is only evaluated in a given Lz sector and both A and B are resticted to a given number of orbitals
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// nbrOrbitalA = number of orbitals that have to be kept for the A part (starting from the leftmost orbital)
// nbrOrbitalA = number of orbitals that have to be kept for the B part (starting from the rightmost orbital)
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix BosonOnSphereShort::EvaluatePartialEntanglementMatrixParticlePartition(int nbrBosonSector, int lzSector, int nbrOrbitalA, int nbrOrbitalB, 
										  RealVector& groundState, bool removeBinomialCoefficient)
{
  int ComplementaryNbrBosonsSector = this->NbrBosons - nbrBosonSector;
  int TotalLzDisk = ConvertLzFromSphereToDisk(this->TotalLz, this->NbrBosons, this->LzMax);
  int LzADisk = ConvertLzFromSphereToDisk(lzSector, nbrBosonSector, nbrOrbitalA - 1);
  if ((LzADisk < 0) || (LzADisk > ((nbrOrbitalA - 1) * nbrBosonSector)))
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
  int LzBDisk = (TotalLzDisk - LzADisk) - ComplementaryNbrBosonsSector * (this->LzMax + 1 - nbrOrbitalB);
  if ((LzBDisk < 0) || (LzBDisk > ((nbrOrbitalB - 1) * ComplementaryNbrBosonsSector)))
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
  int LzBSphere = ConvertLzFromDiskToSphere(LzBDisk, ComplementaryNbrBosonsSector, nbrOrbitalB - 1);

  if (nbrBosonSector == 0)
    {
      BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosonsSector, LzBSphere, nbrOrbitalB - 1);
      RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
      int Shift = (this->LzMax + 1 - nbrOrbitalB);
      for (int i = 0; i < TmpHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  int TmpLzMax = this->LzMax + this->NbrBosons - 1;
	  unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[i] << Shift;
	  while ((TmpState >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  int TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpIndex < this->HilbertSpaceDimension)	    
	    TmpEntanglementMatrix.SetMatrixElement(0, i, groundState[TmpIndex]);
	}
      return TmpEntanglementMatrix;
    }
  
  
  if (nbrBosonSector == this->NbrBosons)
    {
      if (lzSector == this->TotalLz)
	{
	  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, nbrOrbitalA - 1);
	  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, 1,true);
	  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	    {
	      int TmpLzMax = nbrOrbitalA + nbrBosonSector - 1;
	      unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[i];
	      while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpIndex < this->HilbertSpaceDimension)
		TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[TmpIndex]);
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }

  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, nbrOrbitalA - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  unsigned long* TmpMonomial2 = 0;
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonsSector];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
  BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosonsSector, LzBSphere, nbrOrbitalB - 1);
  long TmpNbrNonZeroElements = 0;

  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
  double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonsSector] - LogFactorials[nbrBosonSector];
  if (removeBinomialCoefficient == true)
    TmpLogBinomial = 0.0;

  FactorialCoefficient Factorial;
  unsigned long** TmpDestinationHilbertSpaceOccupationNumbers = new unsigned long* [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  unsigned long** TmpDestinationHilbertSpaceMonomial = new unsigned long* [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpDestinationHilbertSpaceLzMax = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpHilbertSpace.HilbertSpaceDimension, true);
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpaceOccupationNumbers[i] = new unsigned long [this->NbrLzValue];
      TmpDestinationHilbertSpaceMonomial[i] = new unsigned long [nbrBosonSector];
      TmpDestinationHilbertSpace.FermionToBoson(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpDestinationHilbertSpaceOccupationNumbers[i], TmpDestinationHilbertSpaceLzMax[i]);
      TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpDestinationHilbertSpaceMonomial[i]);
      unsigned long* TmpOccupationNumber = TmpDestinationHilbertSpaceOccupationNumbers[i];
      int TmpLzMax = TmpDestinationHilbertSpaceLzMax[i];
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpLzMax; ++k)
	TmpFactor += LogFactorials[TmpOccupationNumber[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
       for (int k = 0; k < ComplementaryNbrBosonsSector; ++k)
	 TmpMonomial1[k] += this->LzMax + 1 - nbrOrbitalA;
      TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateLzMax);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateLzMax; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace.TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  TmpMonomial2 = TmpDestinationHilbertSpaceMonomial[j];
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonsSector) && (TmpIndex3 < nbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonsSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonsSector)
		{
		  while ((TmpIndex3 < nbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
		    {
		      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while ( TmpIndex2 < ComplementaryNbrBosonsSector)
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
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
	      TmpFactorial *= 0.5; 	      
	      ++TmpNbrNonZeroElements;
	      double Tmp = exp(TmpFactorial) * groundState[TmpPos];
	      TmpEntanglementMatrix.SetMatrixElement(j, MinIndex, Tmp);
	    }
	}
    }
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      delete[] TmpDestinationHilbertSpaceOccupationNumbers[i];
      delete[] TmpDestinationHilbertSpaceMonomial[i];
    }
  delete[] TmpDestinationHilbertSpaceOccupationNumbers;
  delete[] TmpDestinationHilbertSpaceLzMax;
  delete[] TmpDestinationHilbertSpaceMonomial;
  delete[] TmpMonomial1;
  delete[] TmpMonomial3;
  if (TmpNbrNonZeroElements > 0)
    return TmpEntanglementMatrix;
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}
  
// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition. The entanglement matrix is only evaluated in a given Lz sector.
// and computed from precalculated particle entanglement matrix
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// thetaTop =  inclination angle defining one edge of the cut in degrees
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& BosonOnSphereShort::EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, RealMatrix& entanglementMatrix)
{
  if ((thetaBottom <= thetaTop) || (phiRange <= 0.0))
    {
      for (int i = 0; i < entanglementMatrix.GetNbrRow(); ++i)
	for (int j = 0; j < entanglementMatrix.GetNbrColumn(); ++j)
	  entanglementMatrix(i, j) = 0.0;
      return entanglementMatrix;
    }
  
  thetaTop *= M_PI / 180.0;
  thetaBottom *= M_PI / 180.0;
  phiRange /= 360.0;
  
  double* IncompleteBetaThetaTop = 0;
  double* IncompleteBetaThetaBottom = 0;

  this->EvaluatePartialDensityMatrixRealSpacePartitionCoefficient(this->LzMax, thetaTop, thetaBottom, IncompleteBetaThetaTop, IncompleteBetaThetaBottom);
  
  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;
  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, this->LzMax);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];

  BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosonSector, this->TotalLz - lzSector, this->LzMax);
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpMonomial3);
      double Tmp = 0.0;
      Tmp = 0.5 * nbrBosonSector * log(phiRange);      
      for( int j = 0; j < nbrBosonSector; j++)
	{
	  Tmp += log( IncompleteBetaThetaBottom[TmpMonomial3[j]] - IncompleteBetaThetaTop[TmpMonomial3[j]]);
	}
      Tmp = exp(0.5 * Tmp);
      for (int j = 0; j < TmpHilbertSpace.HilbertSpaceDimension; ++j)          
	entanglementMatrix(i, j) *= Tmp;      
    }
  
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
      double FormFactor = 0.0;
      for (int i=0; i < ComplementaryNbrBosonSector; i++)
	FormFactor += log(1.0 - IncompleteBetaThetaBottom[TmpMonomial1[i]] + IncompleteBetaThetaTop[TmpMonomial1[i]] + (1.0 - phiRange) * (IncompleteBetaThetaBottom[TmpMonomial1[i]] - IncompleteBetaThetaTop[TmpMonomial1[i]]) );
      FormFactor = exp(0.5 * FormFactor);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	entanglementMatrix(j, MinIndex) *= FormFactor; 
    }
     
  delete[] TmpMonomial1;
  delete[] TmpMonomial3;

  return entanglementMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnSphereShort::EvaluatePartialDensityMatrixGenericRealSpacePartition (int nbrFermionSector, int lzSector, int nbrOrbitalA, double* weightOrbitalA, 
											       int nbrOrbitalB, double* weightOrbitalB, RealVector& groundState, 
											       AbstractArchitecture* architecture)
{
  cout << "error, BosonOnSphereShort::EvaluatePartialDensityMatrixGenericRealSpacePartition is not implemented" << endl;
  return RealSymmetricMatrix();
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& BosonOnSphereShort::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrBosonSector, int lzSector, 
														   int nbrOrbitalA, double* weightOrbitalA, 
														   int nbrOrbitalB, double* weightOrbitalB, RealMatrix& entanglementMatrix)
{
  int ComplementaryNbrBosonsSector = this->NbrBosons - nbrBosonSector;
  int TotalLzDisk = ConvertLzFromSphereToDisk(this->TotalLz, this->NbrBosons, this->LzMax);
  int LzADisk = ConvertLzFromSphereToDisk(lzSector, nbrBosonSector, nbrOrbitalA - 1);
  if ((LzADisk < 0) || (LzADisk > ((nbrOrbitalA - 1) * nbrBosonSector)))
    {
      return entanglementMatrix;	  
    }
  int LzBDisk = (TotalLzDisk - LzADisk) - ComplementaryNbrBosonsSector * (this->LzMax + 1 - nbrOrbitalB);
  if ((LzBDisk < 0) || (LzBDisk > ((nbrOrbitalB - 1) * ComplementaryNbrBosonsSector)))
    {
      return entanglementMatrix;	  
    }
  int LzBSphere = ConvertLzFromDiskToSphere(LzBDisk, ComplementaryNbrBosonsSector, nbrOrbitalB - 1);
  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, nbrOrbitalA - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonsSector];
  unsigned long* TmpMonomial3 = new unsigned long [nbrBosonSector];
  
  BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosonsSector, LzBSphere, nbrOrbitalB - 1);
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpMonomial3);
      double Tmp = 1.0;
      for (int j = 0; j < nbrBosonSector; j++)
	{
	  Tmp *= weightOrbitalA[TmpMonomial3[j]];
	}
      for (int j = 0; j < TmpHilbertSpace.HilbertSpaceDimension; ++j)          
	entanglementMatrix(i, j) *= Tmp;      
    }
  
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
      double FormFactor = 1.0;
      for (int i = 0; i < ComplementaryNbrBosonsSector; i++)
	FormFactor *= weightOrbitalB[TmpMonomial1[i]];
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	entanglementMatrix(j, MinIndex) *= FormFactor; 
    }
  
  delete[] TmpMonomial1;
  delete[] TmpMonomial3;
  
  return entanglementMatrix;
}
  
// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& BosonOnSphereShort::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0;
  if (reference >= 0l)
    {
      Factor /= state[reference];
    }
  else
    {
      reference = 0l;
    }
  this->ConvertToMonomial(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], TmpMonomialReference);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], 
		       this->TemporaryState, this->TemporaryStateLzMax);
  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialDivide(this->TemporaryState[k]);
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	{
	  while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrBosons)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrBosons)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	  ++Index2;
	}
      if (symmetryFactor == true)
	{
	  Factorial = ReferenceFactorial;
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialMultiply(this->TemporaryState[k]);
	  Coefficient *= sqrt(Factorial.GetNumericalValue());
	}
      state[i] *= Coefficient;
    }
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& BosonOnSphereShort::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0;
  if (reference >= 0l)
    Factor = 1.0;
  else
    reference = 0l;
  this->ConvertToMonomial(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], TmpMonomialReference);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
    }
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], 
		       this->TemporaryState, this->TemporaryStateLzMax);
  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialMultiply(this->TemporaryState[k]);

#ifdef VERBOSE_FF
  cout <<"reference monomial: ";
  this->PrintStateMonomial(cout,reference);
#endif
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
#ifdef VERBOSE_FF
      cout <<"target monomial: ";
      this->PrintStateMonomial(cout,i);
      cout << " "<<endl;
#endif
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      if (symmetryFactor == true)
	{
	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	    {
	      while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
#ifdef VERBOSE_FF
		  cout <<"1: "<<InvSqrtCoefficients[TmpMonomialReference[Index1]]<<" ("<<Index1<<") *";
#endif
		  ++Index1;
		  
		}
	      while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
#ifdef VERBOSE_FF
		  cout << " -- 1=2@"<<Index1<<" -- ";
#endif
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
#ifdef VERBOSE_FF
		  cout <<"2: "<<SqrtCoefficients[TmpMonomial[Index2]]<<" ("<<Index2<<")*";
#endif
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrBosons)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
#ifdef VERBOSE_FF
	      cout <<"1: "<<InvSqrtCoefficients[TmpMonomialReference[Index1]]<<" ("<<Index1<<")*";
#endif

	      ++Index1;
	    }
	  while (Index2 < this->NbrBosons)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
#ifdef VERBOSE_FF
	      cout <<"2: "<<SqrtCoefficients[TmpMonomial[Index2]]<<" ("<<Index2<<")*";
#endif
	      ++Index2;
	    }
	  Factorial = ReferenceFactorial;
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialDivide(this->TemporaryState[k]);
	  Coefficient *= sqrt(Factorial.GetNumericalValue());
#ifdef VERBOSE_FF
	  cout <<" fac: "<<sqrt(Factorial.GetNumericalValue())<<endl;
#endif
	}
      else
	{
	  Factorial = ReferenceFactorial;
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialDivide(this->TemporaryState[k]);
	  Coefficient *= sqrt(Factorial.GetNumericalValue());
	}
      state[i] *= Coefficient;
#ifdef VERBOSE_FF
      cout << "total for state: "<<state[i]<<endl;
#endif
    }
  state /= state.Norm();
  return state;
}

// convert a state such that its components are now expressed in the normalized basis, without applying the global normalization to the final state
//
// state = reference to the state to convert
// return value = converted state

RealVector& BosonOnSphereShort::ConvertFromUnnormalizedMonomialNoGlobalNormalization(RealVector& state)
{
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double* SqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt(4.0 * M_PI / (Binomials.GetNumericalCoefficient(this->LzMax, k) * (1.0 + (double) this->LzMax)));
    }
  FactorialCoefficient Factorial;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);

      double Coefficient = 1.0;
      for (int k = 0; k < this->NbrBosons; ++k)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomial[k]];
	}
      Factorial.SetToOne();
      this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	if (this->TemporaryState[k] > 1)
	  Factorial.FactorialDivide(this->TemporaryState[k]);
      Coefficient *= sqrt(Factorial.GetNumericalValue());
      state[i] *= Coefficient;
    }
  delete[] TmpMonomial;
  delete[] SqrtCoefficients;
  return state;
}

// convert a state such that its components, given in the conformal limit,  are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// return value = converted state

RealVector& BosonOnSphereShort::ConvertFromConformalLimit(RealVector& state, long reference)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0;
  if (reference >= 0l)
    Factor = 1.0;
  else
    reference = 0l;
  this->ConvertToMonomial(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], TmpMonomialReference);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;

      while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	{
	  while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrBosons)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrBosons)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	  ++Index2;
	}
      state[i] *= Coefficient;
    }
  state /= state.Norm();
  return state;
}

// normalize Jack with respect to cylinder basis
//
// state = reference to the Jack state to normalize
// aspect = aspect ratio of cylinder
// return value = normalized state

RealVector& BosonOnSphereShort::NormalizeJackToCylinder(RealVector& state, double aspect)
{
  long double Pi_L = 3.14159265358979323846264338328L;
  long double Length = sqrtl((long double)2.0 * Pi_L * (long double)(this->LzMax + 1) * (long double) aspect);
  cout<<"L= "<< Length << " r= "<<aspect<<endl;
  long double kappa = ((long double) 2.0) * Pi_L / Length;
  long double Norm = (long double)0.0;
  long double* LogFactorials = new long double [this->NbrBosons + 1];
  LogFactorials[0] = (long double) 0.0;
  LogFactorials[1] = (long double) 0.0;
  for (int i = 2; i <= this->NbrBosons; ++i)
    {
      LogFactorials[i] = LogFactorials[i - 1] + logl((long double) i);
    }
  for (int i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
      long double Sum2MSquare = (long double) 0.0;
      long double TmpFactorial = (long double) 0.0;
      int j= 0;
      for (int j = 0; j <= this->TemporaryStateLzMax; ++j)
	{
	  if (this->TemporaryState[j] > 0)
	    {
	      Sum2MSquare += (((j - 0.5 * LzMax) * (j - 0.5 * LzMax)) *  this->TemporaryState[j]);
	      TmpFactorial -= LogFactorials[this->TemporaryState[j]]; 
	    }
	}
      state[i] *= expl((((long double)0.5) * kappa * kappa *Sum2MSquare) + (((long double)0.5) * TmpFactorial));      
      Norm += state[i] * state[i];
   }
  delete[] LogFactorials;
  cout << "Norm= "<< Norm << endl;
  state /= sqrtl(Norm); 
  return state;
}

// fuse two states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// leftVector = reference on the vector whose Hilbert space will be fuse to the left
// rightVector = reference on the vector whose Hilbert space will be fuse to the right
// padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
// leftSpace = point to the Hilbert space that will be fuse to the left
// rightSpace = point to the Hilbert space that will be fuse to the right
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

RealVector& BosonOnSphereShort::FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
					    ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
					    bool symmetrizedFlag, double coefficient)
{
  BosonOnSphereShort* LeftSpace = (BosonOnSphereShort*) leftSpace;
  BosonOnSphereShort* RightSpace = (BosonOnSphereShort*) rightSpace;
  if (padding > -2)
    {
      int StateShift = RightSpace->FermionBasis->LzMax + padding + 2;
      for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState1 = LeftSpace->FermionBasis->StateDescription[i] << StateShift;
	  double Coefficient = coefficient * leftVector[i];
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while ((TmpState1 >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  if (symmetrizedFlag == false)
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
		  TmpState2 |= TmpState1;
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
		  outputVector[TmpIndex] = Coefficient2;
		}
	    }
	  else
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
		  TmpState2 |= TmpState1;
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
		  outputVector[TmpIndex] = Coefficient2;
		  unsigned long TmpState3 = this->FermionBasis->GetSymmetricState(TmpState2);
		  if (TmpState3 != TmpState2)
		    {
		      int TmpLzMax2 = this->FermionBasis->LzMax;
		      while ((TmpState3 >> TmpLzMax2) == 0x0ul)
			--TmpLzMax2;
		      TmpIndex = this->FermionBasis->FindStateIndex(TmpState3, TmpLzMax2);
		      outputVector[TmpIndex] = Coefficient2;      
		    }
		}
	    }
	}
    }
  else
    {
      int StateShift = RightSpace->LzMax + padding + 1;
      for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
	{
	  LeftSpace->FermionToBoson(LeftSpace->FermionBasis->StateDescription[i], LeftSpace->FermionBasis->StateLzMax[i], LeftSpace->TemporaryState, LeftSpace->TemporaryStateLzMax);	  
	  this->TemporaryStateLzMax = LeftSpace->TemporaryStateLzMax + StateShift;
	  if (this->TemporaryStateLzMax > this->LzMax)
	    this->TemporaryStateLzMax = this->LzMax;
	  double Coefficient = coefficient * leftVector[i];
	  if (symmetrizedFlag == false)
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  LeftSpace->FermionToBoson(RightSpace->FermionBasis->StateDescription[j], RightSpace->FermionBasis->StateLzMax[j], RightSpace->TemporaryState, RightSpace->TemporaryStateLzMax);
		  for (int k = 0; k <= RightSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k] = RightSpace->TemporaryState[k];
		  for (int k = RightSpace->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
		    this->TemporaryState[k] = 0;
		  for (int k = 0; k <= LeftSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k + StateShift] += LeftSpace->TemporaryState[k];
		  unsigned long TmpState2 = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  int TmpLzMax = this->FermionBasis->LzMax;
		  while ((TmpState2 >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
		  if (TmpIndex == this->HilbertSpaceDimension)
		    {
		      cout << "error while merging components " << i << " and " << j << endl;
		    }
		  outputVector[TmpIndex] = Coefficient2;
		}
	    }
	  else
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  LeftSpace->FermionToBoson(RightSpace->FermionBasis->StateDescription[j], RightSpace->FermionBasis->StateLzMax[j], RightSpace->TemporaryState, RightSpace->TemporaryStateLzMax);
		  for (int k = 0; k <= RightSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k] = RightSpace->TemporaryState[k];
		  for (int k = RightSpace->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
		    this->TemporaryState[k] = 0;
		  for (int k = 0; k <= LeftSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k + StateShift] += LeftSpace->TemporaryState[k];
		  unsigned long TmpState2 = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  int TmpLzMax = this->FermionBasis->LzMax;
		  while ((TmpState2 >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
		  if (TmpIndex == this->HilbertSpaceDimension)
		    {
		      cout << "error while merging components " << i << " and " << j << endl;
		    }
		  outputVector[TmpIndex] = Coefficient2;
		  unsigned long TmpState3 = this->FermionBasis->GetSymmetricState(TmpState2);
		  if (TmpState3 != TmpState2)
		    {
		      int TmpLzMax2 = this->FermionBasis->LzMax;
		      while ((TmpState3 >> TmpLzMax2) == 0x0ul)
			--TmpLzMax2;
		      TmpIndex = this->FermionBasis->FindStateIndex(TmpState3, TmpLzMax2);
		      outputVector[TmpIndex] = Coefficient2;      
		    }
		}
	    }
	}
    }

  return outputVector;
}

// fuse two states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// leftVector = reference on the vector whose Hilbert space will be fuse to the left
// rightVector = reference on the vector whose Hilbert space will be fuse to the right
// padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
// leftSpace = point to the Hilbert space that will be fuse to the left
// rightSpace = point to the Hilbert space that will be fuse to the right
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

LongRationalVector& BosonOnSphereShort::FuseStates (LongRationalVector& outputVector, LongRationalVector& leftVector, LongRationalVector& rightVector, int padding, 
						    ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
						    bool symmetrizedFlag, LongRational& coefficient)
{
  BosonOnSphereShort* LeftSpace = (BosonOnSphereShort*) leftSpace;
  BosonOnSphereShort* RightSpace = (BosonOnSphereShort*) rightSpace;
  if (padding > -2)
    {
      int StateShift = RightSpace->FermionBasis->LzMax + padding + 2;
      for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState1 = LeftSpace->FermionBasis->StateDescription[i] << StateShift;
	  LongRational Coefficient = coefficient * leftVector[i];
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while ((TmpState1 >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  if (symmetrizedFlag == false)
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
		  TmpState2 |= TmpState1;
		  LongRational Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
		  if (TmpIndex == this->HilbertSpaceDimension)
		    {
		      cout << "error while merging components " << i << " and " << j << endl;
		    }
		  outputVector[TmpIndex] = Coefficient2;
		}
	    }
	  else
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
		  TmpState2 |= TmpState1;
		  LongRational Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
		  if (TmpIndex == this->HilbertSpaceDimension)
		    {
		      cout << "error while merging components " << i << " and " << j << endl;
		    }
		  outputVector[TmpIndex] = Coefficient2;
		  unsigned long TmpState3 = this->FermionBasis->GetSymmetricState(TmpState2);
		  if (TmpState3 != TmpState2)
		    {
		      int TmpLzMax2 = this->FermionBasis->LzMax;
		      while ((TmpState3 >> TmpLzMax2) == 0x0ul)
		    --TmpLzMax2;
		      TmpIndex = this->FermionBasis->FindStateIndex(TmpState3, TmpLzMax2);
		      outputVector[TmpIndex] = Coefficient2;      
		    }
		}
	    }
	}
    }
  else
    {
      int StateShift = RightSpace->LzMax + padding + 1;
      cout << "StateShift=" << StateShift << endl;
      for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
	{
	  LeftSpace->FermionToBoson(LeftSpace->FermionBasis->StateDescription[i], LeftSpace->FermionBasis->StateLzMax[i], LeftSpace->TemporaryState, LeftSpace->TemporaryStateLzMax);	  
	  this->TemporaryStateLzMax = LeftSpace->TemporaryStateLzMax + StateShift;
	  if (this->TemporaryStateLzMax > this->LzMax)
	    this->TemporaryStateLzMax = this->LzMax;
	  LongRational Coefficient = coefficient * leftVector[i];
	  if (symmetrizedFlag == false)
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  LeftSpace->FermionToBoson(RightSpace->FermionBasis->StateDescription[j], RightSpace->FermionBasis->StateLzMax[j], RightSpace->TemporaryState, RightSpace->TemporaryStateLzMax);
		  for (int k = 0; k <= RightSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k] = RightSpace->TemporaryState[k];
		  for (int k = RightSpace->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
		    this->TemporaryState[k] = 0;
		  for (int k = 0; k <= LeftSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k + StateShift] += LeftSpace->TemporaryState[k];
		  unsigned long TmpState2 = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
		  LongRational Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  int TmpLzMax = this->FermionBasis->LzMax;
		  while ((TmpState2 >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
		  if (TmpIndex == this->HilbertSpaceDimension)
		    {
		      cout << "error while merging components " << i << " ( ";
		      for (int k = 0; k <= LeftSpace->TemporaryStateLzMax; ++k)
			cout << LeftSpace->TemporaryState[k] << " ";
		      for (int k =  LeftSpace->TemporaryStateLzMax + 1; k <= LeftSpace->LzMax; ++k)
			cout << "0 ";
		      cout << ") and " << j << " ( ";
		      for (int k = 0; k <= RightSpace->TemporaryStateLzMax; ++k)
			cout << RightSpace->TemporaryState[k] << " ";
		      for (int k = RightSpace->TemporaryStateLzMax + 1; k <= RightSpace->LzMax; ++k)
			cout << "0 ";
		      cout << "), merge = ( ";
		      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
			cout << this->TemporaryState[k] << " ";
		      for (int k = this->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
			cout << "0 ";
		      cout << ")" << endl;
		    }
		  outputVector[TmpIndex] = Coefficient2;
		}
	    }
	  else
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  LeftSpace->FermionToBoson(RightSpace->FermionBasis->StateDescription[j], RightSpace->FermionBasis->StateLzMax[j], RightSpace->TemporaryState, RightSpace->TemporaryStateLzMax);
		  for (int k = 0; k <= RightSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k] = RightSpace->TemporaryState[k];
		  for (int k = RightSpace->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
		    this->TemporaryState[k] = 0;
		  for (int k = 0; k <= LeftSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k + StateShift] += LeftSpace->TemporaryState[k];
		  unsigned long TmpState2 = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
		  LongRational Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  int TmpLzMax = this->FermionBasis->LzMax;
		  while ((TmpState2 >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
		  if (TmpIndex == this->HilbertSpaceDimension)
		    {
		      cout << "error while merging components " << i << " ( ";
		      for (int k = 0; k <= RightSpace->TemporaryStateLzMax; ++k)
			cout << RightSpace->TemporaryState[k] << " ";
		      cout << ") and " << j << " ( ";
		      for (int k = 0; k <= LeftSpace->TemporaryStateLzMax; ++k)
			cout << LeftSpace->TemporaryState[k] << " ";
		      cout << "), merge = ( ";
		      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
			cout << this->TemporaryState[k] << " ";
		      cout << ")" << endl;
		    }
		  outputVector[TmpIndex] = Coefficient2;
		  unsigned long TmpState3 = this->FermionBasis->GetSymmetricState(TmpState2);
		  if (TmpState3 != TmpState2)
		    {
		      int TmpLzMax2 = this->FermionBasis->LzMax;
		      while ((TmpState3 >> TmpLzMax2) == 0x0ul)
			--TmpLzMax2;
		      TmpIndex = this->FermionBasis->FindStateIndex(TmpState3, TmpLzMax2);
		      outputVector[TmpIndex] = Coefficient2;      
		    }
		}
	    }
	}
    }
  return outputVector;
}

// fuse multiple states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// nbrInputVectors = number of input vectors
// inputVectors = input vectors whose Hilbert space will be fuse from  left to right
// paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
// inputSpaces = point to the Hilbert space that will be fuse to the left
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

RealVector& BosonOnSphereShort::FuseMultipleStates (RealVector& outputVector, int nbrInputVectors, RealVector* inputVectors, int* paddings, 
						    ParticleOnSphere** inputSpaces, bool symmetrizedFlag, double coefficient)
{
  BosonOnSphereShort** InputSpaces = (BosonOnSphereShort**) inputSpaces;
  for (long i = 0; i <  InputSpaces[0]->LargeHilbertSpaceDimension; ++i)
    {
      this->CoreFuseMultipleStates(outputVector, nbrInputVectors, inputVectors, paddings, InputSpaces, 1, 
				   ((BosonOnSphereShort*) InputSpaces[0])->FermionBasis->StateDescription[i], 
				   ((BosonOnSphereShort*) InputSpaces[0])->FermionBasis->LzMax + paddings[0] + 2,
				   coefficient  * inputVectors[0][i], symmetrizedFlag);	    
      
    }
  return outputVector;
}

// fuse multiple states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// nbrInputVectors = number of input vectors
// inputVectors = input vectors whose Hilbert space will be fuse from  left to right
// paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
// inputSpaces = point to the Hilbert space that will be fuse to the left
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

LongRationalVector& BosonOnSphereShort::FuseMultipleStates (LongRationalVector& outputVector, int nbrInputVectors, LongRationalVector* inputVectors, int* paddings, 
							    ParticleOnSphere** inputSpaces, bool symmetrizedFlag, LongRational& coefficient)
{
  BosonOnSphereShort** InputSpaces = (BosonOnSphereShort**) inputSpaces;
  for (long i = 0; i <  InputSpaces[0]->LargeHilbertSpaceDimension; ++i)
    {
      this->CoreFuseMultipleStates(outputVector, nbrInputVectors, inputVectors, paddings, InputSpaces, 1, 
				   ((BosonOnSphereShort*) InputSpaces[0])->FermionBasis->StateDescription[i], 
				   ((BosonOnSphereShort*) InputSpaces[0])->FermionBasis->LzMax + paddings[0] + 2,
				   coefficient  * inputVectors[0][i], symmetrizedFlag);	    
      
    }
  return outputVector;
}

// core part of multiple state fuse 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// nbrInputVectors = number of input vectors
// inputVectors = input vectors whose Hilbert space will be fuse from  left to right
// paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
// inputSpaces = point to the Hilbert space that will be fuse to the left
// currentPosition = index of the current space to fuse
// currentState = current fermionic state obtained by fusing previous states
// currentCoefficient = current multiplicative coefficient
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry

void BosonOnSphereShort::CoreFuseMultipleStates (RealVector& outputVector, int nbrInputVectors, RealVector* inputVectors, int* paddings, 
						 BosonOnSphereShort** inputSpaces, int currentPosition, unsigned long currentState, int currentPadding, 
						 double currentCoefficient, bool symmetrizedFlag)
{
  if (currentPosition < (nbrInputVectors - 1))
    {
      for (long i = 0; i <  inputSpaces[currentPosition]->LargeHilbertSpaceDimension; ++i)
	{
	  this->CoreFuseMultipleStates(outputVector, nbrInputVectors, inputVectors, paddings, inputSpaces, currentPosition + 1, 
				       currentState | (inputSpaces[currentPosition]->FermionBasis->StateDescription[i] << currentPadding), 
				       currentPadding + inputSpaces[currentPosition]->FermionBasis->LzMax + paddings[currentPosition] + 2,
				       currentCoefficient * inputVectors[currentPosition][i], symmetrizedFlag);	    
	  
	}
    }
  else
    {
      FermionOnSphere* RightSpace = inputSpaces[currentPosition]->FermionBasis;
      RealVector& RightVector = inputVectors[currentPosition];
      if (symmetrizedFlag == false)
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = RightSpace->StateDescription[j] << currentPadding;
	      TmpState |= currentState;
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      outputVector[TmpIndex] = currentCoefficient * RightVector[j];
	    }
	}
      else
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = RightSpace->StateDescription[j] << currentPadding;
	      TmpState |= currentState;
	      double Coefficient = currentCoefficient;
	      Coefficient *= RightVector[j];	  
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      outputVector[TmpIndex] = Coefficient;
	      unsigned long TmpState2 = this->FermionBasis->GetSymmetricState(TmpState);
	      if (TmpState != TmpState2)
		{
		  TmpLzMax = this->FermionBasis->LzMax;
		  while ((TmpState2 >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
		  outputVector[TmpIndex] = Coefficient;      
		}
	    }
	}
    }
}

// core part of multiple state fuse 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// nbrInputVectors = number of input vectors
// inputVectors = input vectors whose Hilbert space will be fuse from  left to right
// paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
// inputSpaces = point to the Hilbert space that will be fuse to the left
// currentPosition = index of the current space to fuse
// currentState = current fermionic state obtained by fusing previous states
// currentCoefficient = current multiplicative coefficient
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry

void BosonOnSphereShort::CoreFuseMultipleStates (LongRationalVector& outputVector, int nbrInputVectors, LongRationalVector* inputVectors, int* paddings, 
						 BosonOnSphereShort** inputSpaces, int currentPosition, unsigned long currentState, int currentPadding, 
						 LongRational currentCoefficient, bool symmetrizedFlag)
{
  if (currentPosition < (nbrInputVectors - 1))
    {
      for (long i = 0; i <  inputSpaces[currentPosition]->LargeHilbertSpaceDimension; ++i)
	{
	  this->CoreFuseMultipleStates(outputVector, nbrInputVectors, inputVectors, paddings, inputSpaces, currentPosition + 1, 
				       currentState | (inputSpaces[currentPosition]->FermionBasis->StateDescription[i] << currentPadding), 
				       currentPadding + inputSpaces[currentPosition]->FermionBasis->LzMax + paddings[currentPosition] + 2,
				       currentCoefficient * inputVectors[currentPosition][i], symmetrizedFlag);	    
	  
	}
    }
  else
    {
      FermionOnSphere* RightSpace = inputSpaces[currentPosition]->FermionBasis;
      LongRationalVector& RightVector = inputVectors[currentPosition];
      if (symmetrizedFlag == false)
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = RightSpace->StateDescription[j] << currentPadding;
	      TmpState |= currentState;
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      outputVector[TmpIndex] = currentCoefficient * RightVector[j];
	    }
	}
      else
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState = RightSpace->StateDescription[j] << currentPadding;
	      TmpState |= currentState;
	      LongRational Coefficient = currentCoefficient;
	      Coefficient *= RightVector[j];	  
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      outputVector[TmpIndex] = Coefficient;
	      unsigned long TmpState2 = this->FermionBasis->GetSymmetricState(TmpState);
	      if (TmpState != TmpState2)
		{
		  TmpLzMax = this->FermionBasis->LzMax;
		  while ((TmpState2 >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
		  outputVector[TmpIndex] = Coefficient;      
		}
	    }
	}
    }
}

// use product rule to produce part of the components of a system from a smaller one
//
// outputVector = reference on the vector which will contain the product rule state  (without zeroing components which do not occur in the fusion)
// inputVector = reference on the vector associated to the smaller system
// inputSpace = pointer to the Hilbert space of the smaller system
// commonPattern = array describing the shared leftmost pattern between the n-body states in both the smaller and larger system sizes
// commonPatterSize = number of elements in the commonPattern array
// addedPattern = array describing the pattern that has to be inserted to go from the smaller system to the larger one
// addedPatterSize = number of elements in the addedPattern array
// coefficient = multiplicqtive fqctor to go fron the component of the smaller system to the larger one
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// return value = reference on the product rule state

RealVector& BosonOnSphereShort::ProductRules (RealVector& outputVector, RealVector& inputVector, ParticleOnSphere* inputSpace, 
					      int* commonPattern, int commonPatterSize, int* addedPattern, int addedPatterSize,
					      double coefficient, bool symmetrizedFlag)
{
  BosonOnSphereShort* InputSpace = (BosonOnSphereShort*) inputSpace;
  int NbrParticlesCommonPattern = 0;
  unsigned long InputPattern = 0x0ul;
  int TmpIndex = 0;
  for (int i = commonPatterSize - 1; i >= 0; --i)
    {
      for (int j = 0; j < commonPattern[i]; ++j)
	{
	  InputPattern |= 0x1ul << TmpIndex;
	  ++TmpIndex;
	}
      ++TmpIndex;
      NbrParticlesCommonPattern += commonPattern[i];
    }
  int NbrParticlesAddedPattern = 0;
  unsigned long OutputPattern = 0x0ul;
  TmpIndex = 0;
  for (int i = addedPatterSize - 1; i >= 0; --i)
    {
      for (int j = 0; j < addedPattern[i]; ++j)
	{
	  OutputPattern |= 0x1ul << TmpIndex;
	  ++TmpIndex;
	}
      ++TmpIndex;
      NbrParticlesAddedPattern += addedPattern[i];
    }
  unsigned long InputMask = ((0x1ul << (commonPatterSize + NbrParticlesCommonPattern)) - 1ul) << (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 1);
  unsigned long InputMask2 = ~InputMask;  
  OutputPattern |= InputPattern << (addedPatterSize + NbrParticlesAddedPattern);
  OutputPattern <<= (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 2);
  InputPattern <<= (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 2); 
  cout << hex << InputMask << " " << OutputPattern << " " << InputPattern << dec << endl;
  int OutputLzMax = this->FermionBasis->LzMax;
  while (((OutputPattern >> OutputLzMax) & 0x1ul) == 0x0ul)
    --OutputLzMax;
  long Count = 0l;
  if (symmetrizedFlag == false)
    {
      for (int i = 0; i <  InputSpace->HilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState1 = InputSpace->FermionBasis->StateDescription[i];
	  //	  cout << i << " " << hex << TmpState1 << dec << endl;
	  if ((TmpState1 & InputMask) == InputPattern)
	    {
	      TmpState1 &= InputMask2;
	      TmpState1 |= OutputPattern;
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState1, OutputLzMax);
	      double& TmpCoef = outputVector[TmpIndex];
	      if (TmpCoef == 0.0)
		++Count;
	      TmpCoef = coefficient * inputVector[i];	  
	    }
	}
    }
  else
    {
      for (int i = 0; i <  InputSpace->HilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState1 = InputSpace->FermionBasis->StateDescription[i];
	  if ((TmpState1 & InputMask) == InputPattern)
	    {
	      TmpState1 &= InputMask2;
	      TmpState1 |= OutputPattern;
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState1, OutputLzMax);
	      double& TmpCoef = outputVector[TmpIndex];
	      if (TmpCoef == 0.0)
		++Count;
	      double TmpCoef3 = coefficient * inputVector[i];
	      TmpCoef = TmpCoef3;	  
	      unsigned long TmpState2 = this->FermionBasis->GetSymmetricState(TmpState1);
	      if (TmpState2 != TmpState1)
		{
		  int TmpLzMax2 = this->FermionBasis->LzMax;
		  while ((TmpState2 >> TmpLzMax2) == 0x0ul)
		    --TmpLzMax2;
		  TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax2);
		  double& TmpCoef2 = outputVector[TmpIndex];
                  if (TmpCoef2 == 0.0)
                    ++Count;
                  TmpCoef2 = TmpCoef3;
		}
	    }

	}
    }
  cout << "nbr of newly added components : " << Count << endl;
  return outputVector;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

double BosonOnSphereShort::JackSqrNormalization (RealVector& outputVector, long minIndex, long nbrComponents)
{
  double SqrNorm = 0.0;
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = ((unsigned long) this->LzMax) >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  for (long i = minIndex; i < MaxIndex; ++i)
    {
      Factorial.SetToOne();
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = 0; k < this->NbrBosons; ++k)
	{
	  if (HalfLzMax < TmpMonomial[k])
	    Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomial[k]);
	  else
	    if (HalfLzMax > TmpMonomial[k])
	      Factorial.PartialFactorialDivide(TmpMonomial[k] + 1, HalfLzMax);
	}	      
      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	if (this->TemporaryState[k] > 1)
	  Factorial.FactorialDivide(this->TemporaryState[k]);
      SqrNorm +=(outputVector[i] * outputVector[i]) * Factorial.GetNumericalValue();
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  cout << endl;
  delete[] TmpMonomial;
  return SqrNorm;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

LongRational BosonOnSphereShort::JackSqrNormalization (LongRationalVector& outputVector, long minIndex, long nbrComponents)
{
  LongRational SqrNorm = 0l;
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = ((unsigned long) this->LzMax) >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  for (long i = minIndex; i < MaxIndex; ++i)
    {
      Factorial.SetToOne();
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = 0; k < this->NbrBosons; ++k)
	{
	  if (HalfLzMax < TmpMonomial[k])
	    Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomial[k]);
	  else
	    if (HalfLzMax > TmpMonomial[k])
	      Factorial.PartialFactorialDivide(TmpMonomial[k] + 1, HalfLzMax);
	}	      
      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	if (this->TemporaryState[k] > 1)
	  Factorial.FactorialDivide(this->TemporaryState[k]);
      SqrNorm += (outputVector[i] * outputVector[i]) * Factorial.GetLongRationalValue();
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  cout << endl;
  delete[] TmpMonomial;
  return SqrNorm;
}

// compute part of the Jack polynomial scalar product in a given range of indices
//
// state1 = reference on the first unnormalized Jack polynomial
// state2 = reference on the second unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

double BosonOnSphereShort::JackScalarProduct (RealVector& state1, RealVector& state2, long minIndex, long nbrComponents)
{
  double SqrNorm = 0.0;
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = ((unsigned long) this->LzMax) >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  for (long i = minIndex; i < MaxIndex; ++i)
    {
      Factorial.SetToOne();
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = 0; k < this->NbrBosons; ++k)
	{
	  if (HalfLzMax < TmpMonomial[k])
	    Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomial[k]);
	  else
	    if (HalfLzMax > TmpMonomial[k])
	      Factorial.PartialFactorialDivide(TmpMonomial[k] + 1, HalfLzMax);
	}	      
      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	if (this->TemporaryState[k] > 1)
	  Factorial.FactorialDivide(this->TemporaryState[k]);
      SqrNorm +=(state1[i] * state2[i]) * Factorial.GetNumericalValue();
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  cout << endl;
  delete[] TmpMonomial;
  return SqrNorm;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state1 = reference on the first unnormalized Jack polynomial
// state2 = reference on the second unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

LongRational BosonOnSphereShort::JackScalarProduct (LongRationalVector& state1, LongRationalVector& state2, long minIndex, long nbrComponents)
{
  LongRational SqrNorm = 0l;
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = ((unsigned long) this->LzMax) >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  for (long i = minIndex; i < MaxIndex; ++i)
    {
      Factorial.SetToOne();
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = 0; k < this->NbrBosons; ++k)
	{
	  if (HalfLzMax < TmpMonomial[k])
	    Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomial[k]);
	  else
	    if (HalfLzMax > TmpMonomial[k])
	      Factorial.PartialFactorialDivide(TmpMonomial[k] + 1, HalfLzMax);
	}	      
      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	if (this->TemporaryState[k] > 1)
	  Factorial.FactorialDivide(this->TemporaryState[k]);
      SqrNorm += (state1[i] * state2[i]) * Factorial.GetLongRationalValue();
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  cout << endl;
  delete[] TmpMonomial;
  return SqrNorm;
}

// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the  component
int BosonOnSphereShort::GetLzValue(int j)
{
  return this->TotalLz;
}

// compute all Kostka coefficients for a given Schur polynomial 
//
// index = index of the partition that describe the Schur polynomial 
// kostkaVector = vector where kostka numbers will be stored
// return value = number of kostka coefficients

long BosonOnSphereShort::KostkaForOneSchur(long index, RealVector& kostkaVector)
{	
  unsigned long * InitialState = new unsigned long[this->NbrLzValue];
  int finalStateLzMax;
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], InitialState, finalStateLzMax);
  int* ReferenceState = new int[this->NbrLzValue];
  for(int i = 0; i <= finalStateLzMax; i++)
    {
      ReferenceState[i] = (int) InitialState[i];
    }
  for(int i = finalStateLzMax + 1; i < this->NbrLzValue; i++)
    {
      ReferenceState[i] = 0;
    }
  BosonOnSphereHaldaneBasisShort* SqueezedSpace = new BosonOnSphereHaldaneBasisShort(this->NbrBosons, this->TotalLz, this->LzMax, 
										     ReferenceState);
  RealVector JackVector(SqueezedSpace->GetHilbertSpaceDimension(), true);
  SqueezedSpace->GenerateJackPolynomial(JackVector, 1.0);
  int Dimension = SqueezedSpace->GetHilbertSpaceDimension();
  
  for (int i = 0; i < Dimension; i++)
    {
      int TmpLzMax = this->LzMax + this->NbrBosons - 1;
      unsigned long TmpState = SqueezedSpace->FermionBasis->StateDescription[i];
      while ((TmpState >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      kostkaVector[this->FermionBasis->FindStateIndex(TmpState, TmpLzMax)] = JackVector[i];
    }
  delete SqueezedSpace;
  delete[] ReferenceState;
  return Dimension;
}

// divide a set of fermionic vectors by a Jastrow factor to get their bosonic counterpart
//
// sourceVector = array of fermionic statesc be divided by a jastrow factor
// outputVector = array of where bosonic states will be stored
// firstComponent = index of the first component to transform 
// nbrComponent = number of components to compute
// nbrStates  = number of states to handle

void BosonOnSphereShort::DivideByJastrow(RealVector* sourceVector, RealVector* OutputVector, long firstComponent,
					 long nbrComponent, int nbrStates)
{
  long IndiceMax = firstComponent + nbrComponent;
  RealVector KostkaVector (this->GetHilbertSpaceDimension(), true);
  for(long i = firstComponent; i < IndiceMax; i++)
    {
      this->KostkaForOneSchur(i,KostkaVector);
      for(long j = i; j < this->GetHilbertSpaceDimension(); j++)
	{
	  for(int k = 0; k < nbrStates; k++)	
	    OutputVector[k][j] += sourceVector[k][i] * KostkaVector[j];
	}
      KostkaVector.ClearVector();
    }
}

void BosonOnSphereShort::FuseParticlesInState(RealVector& firstState, RealVector& outputVector, BosonOnSphereShort* finalSpace, 
					      long minIndex, long nbrComponents)
{
  long MaxIndex= minIndex + nbrComponents;  
  if (nbrComponents == 0l)
    MaxIndex = finalSpace->GetLargeHilbertSpaceDimension();
  long* WeigthVector = new long[finalSpace->GetLargeHilbertSpaceDimension()];
  unsigned long* IndicesVector = new unsigned long[finalSpace->GetLargeHilbertSpaceDimension()];
  for(long Index = minIndex; Index < MaxIndex; Index++)
    {
      int Size = this->FuseParticlesInMonomial(Index, finalSpace, WeigthVector, IndicesVector);
      for(int i = 0; i < Size; i++)
	{
	  int TmpLzMax = finalSpace->LzMax + finalSpace->NbrBosons - 1;
	  while ((IndicesVector[i] >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  outputVector[finalSpace->FermionBasis->FindStateIndex(IndicesVector[i], TmpLzMax)] += firstState[Index] * WeigthVector[i];
	}
    }
}

// fuse particles two by two in a given monomial
//
// index = monomial index
// finalSpace = space where the fused state lies
// weigthVector = weigths of each fused component
// indicesVector = indices of each fused component
// return value = number of generated components when fusing

int BosonOnSphereShort::FuseParticlesInMonomial(long index, BosonOnSphereShort* finalSpace, long* weigthVector, unsigned long* indicesVector)
{
  FactorialCoefficient Coefficient;
  Coefficient.SetToOne();
  Coefficient.Power2Multiply(finalSpace->NbrBosons);
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  unsigned long* TmpMonomial = new unsigned long[this->NbrBosons];
  unsigned long* FinalMonomial = new unsigned long[finalSpace->NbrBosons];
  this->GetMonomial(index,TmpMonomial);
  int Size = 0;
  Size = this->GeneratePairs(TmpMonomial, weigthVector, indicesVector, FinalMonomial, this->NbrBosons, 0, finalSpace);
  for(int l = 0; l <= TemporaryStateLzMax; l++)
    {
      Coefficient.FactorialDivide(this->TemporaryState[l]);		
    }
  FactorialCoefficient Coefficient1;
  for(int i = 0; i < Size; i++)
    {
      Coefficient1 = Coefficient;
      int TmpLzMax = finalSpace->LzMax + finalSpace->NbrBosons - 1;
      while ((indicesVector[i] >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      finalSpace->FermionToBoson(indicesVector[i], TmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
      for(int p = 0; p <= finalSpace->TemporaryStateLzMax; p++)
	{
	  Coefficient1.FactorialMultiply(finalSpace->TemporaryState[p]);
	}
      if(weigthVector[i] < 0)
	{
	  Coefficient1 *= (-weigthVector[i]);
	  weigthVector[i] = (-Coefficient1.GetIntegerValue());
	}
      else
	{
	  Coefficient1 *= weigthVector[i];
	  weigthVector[i] = Coefficient1.GetIntegerValue();
	}
    }
  return Size;
}
	
int BosonOnSphereShort::GeneratePairs(unsigned long* Monomial, long* weigthVector, unsigned long* indicesVector, unsigned long* FinalMonomial, int remainder, int count, BosonOnSphereShort* finalSpace)
{
  if(remainder == 0)
    {
      unsigned long* TMonomial = new unsigned long[finalSpace->NbrBosons];
      for(int pll = 0; pll < finalSpace->NbrBosons; pll++)
	TMonomial[pll] = FinalMonomial[pll];
      SortArrayDownOrdering(TMonomial,finalSpace->NbrBosons);
      unsigned long Result = finalSpace->ConvertFromMonomial(TMonomial);
      for(int k = 0; k < count; k++)
	{
	  if(Result == indicesVector[k])
	    {
	      weigthVector[k]++;
	      return 0;
	    }
	}
      indicesVector[count] = Result;
      weigthVector[count] = 1l;
      return 1;
    }
  for(int i = 1; i < remainder; i++)
    {
      FinalMonomial[finalSpace->NbrBosons - (remainder / 2)] = Monomial[0] + Monomial[i];
      unsigned long* TmpMonomial = new unsigned long[remainder - 2];
      for(int j = 0; j < (i - 1); j++)
	{
	  TmpMonomial[j] = Monomial[j+1];
	}
      for(int j = i; j < (remainder - 1); j++)
	{
	  TmpMonomial[j - 1] = Monomial[j + 1];
	}
      if(remainder == 2)
	count += this->GeneratePairs(TmpMonomial, weigthVector, indicesVector, FinalMonomial, remainder - 2, count, finalSpace);
      else
	count = this->GeneratePairs(TmpMonomial, weigthVector, indicesVector, FinalMonomial, remainder - 2, count, finalSpace);
    }
  return count;
}


// symmetrized a product of two uncoupled states 
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

RealVector BosonOnSphereShort::SymmetrizeU1U1State (RealVector& leftVector, RealVector& rightVector, BosonOnSphereShort* leftSpace, BosonOnSphereShort* rightSpace, bool unnormalizedBasisFlag, AbstractArchitecture* architecture)
{
  RealVector SymmetrizedVector (this->LargeHilbertSpaceDimension,true);

  FQHESphereSymmetrizeU1U1StateOperation Operation (this, leftSpace, rightSpace, &SymmetrizedVector, &leftVector, &rightVector, unnormalizedBasisFlag);
  Operation.ApplyOperation(architecture);

  if ( unnormalizedBasisFlag == false )
    SymmetrizedVector /= SymmetrizedVector.Norm();

  return SymmetrizedVector;
}


// symmetrized a product of two uncoupled states 
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

void BosonOnSphereShort::SymmetrizeU1U1StateCore (RealVector& symmetrizedVector, RealVector& leftVector, RealVector& rightVector, BosonOnSphereShort* leftSpace, BosonOnSphereShort* rightSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents)
{
  long LastComponent = long(firstComponent + nbrComponents);
  
  FactorialCoefficient Factorial1;
  FactorialCoefficient Factorial2;
  if (unnormalizedBasisFlag == true)
    {
      for (long i = (long) firstComponent; i < LastComponent; ++i)
	{
	  this->FermionToBoson(leftSpace->FermionBasis->StateDescription[i], leftSpace->FermionBasis->StateLzMax[i], 
			       leftSpace->TemporaryState, leftSpace->TemporaryStateLzMax);
	  for (int k = leftSpace->TemporaryStateLzMax + 1;  k <= leftSpace->LzMax; ++k)
	    leftSpace->TemporaryState[k] = 0;
	  double TmpCoefficient = leftVector[i];
	  Factorial1.SetToOne();
	  Factorial1.Power2Divide(leftSpace->NbrBosons);
	  for (int k = 0; k <= leftSpace->TemporaryStateLzMax; ++k)
	    if (leftSpace->TemporaryState[k] > 1)
	      Factorial1.FactorialDivide(leftSpace->TemporaryState[k]);
	  
	  for (long j = 0l; j < rightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      this->FermionToBoson(rightSpace->FermionBasis->StateDescription[j], rightSpace->FermionBasis->StateLzMax[j], 
				   rightSpace->TemporaryState, rightSpace->TemporaryStateLzMax);
	      int k = 0;
	      for (; k <= rightSpace->TemporaryStateLzMax; ++k)
		this->TemporaryState[k] = leftSpace->TemporaryState[k] + rightSpace->TemporaryState[k];
	      this->TemporaryStateLzMax = rightSpace->TemporaryStateLzMax;
	      if (leftSpace->TemporaryStateLzMax > rightSpace->TemporaryStateLzMax)
		{
		  for (; k <= leftSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k] = leftSpace->TemporaryState[k];
		  this->TemporaryStateLzMax = leftSpace->TemporaryStateLzMax;
		}
	      int TmpPos = this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
	      if (TmpPos < this->HilbertSpaceDimension)
		{
		  Factorial2 = Factorial1;
		  for (k = 0; k <= rightSpace->TemporaryStateLzMax; ++k)
		    if (rightSpace->TemporaryState[k] > 1)
		      Factorial2.FactorialDivide(rightSpace->TemporaryState[k]);
		  for (k = 0; k <= this->TemporaryStateLzMax; ++k)
		    if (this->TemporaryState[k] > 1)
		      Factorial2.FactorialMultiply(this->TemporaryState[k]);	      
		  symmetrizedVector[TmpPos] += Factorial2.GetNumericalValue() * TmpCoefficient * rightVector[j];
		}
	    }
	}
    }
  else
    {
      for (long i = (long) firstComponent; i < LastComponent; ++i)
	{
	  this->FermionToBoson(leftSpace->FermionBasis->StateDescription[i], leftSpace->FermionBasis->StateLzMax[i], 
			       leftSpace->TemporaryState, leftSpace->TemporaryStateLzMax);
	  for (int k = leftSpace->TemporaryStateLzMax + 1;  k <= leftSpace->LzMax; ++k)
	    leftSpace->TemporaryState[k] = 0;
	  double TmpCoefficient = leftVector[i];
	  Factorial1.SetToOne();
	  Factorial1.Power2Divide(leftSpace->NbrBosons);
	  for (int k = 0; k <= leftSpace->TemporaryStateLzMax; ++k)
	    if (leftSpace->TemporaryState[k] > 1)
	      Factorial1.FactorialDivide(leftSpace->TemporaryState[k]);
	  
	  for (long j = 0l; j < rightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      this->FermionToBoson(rightSpace->FermionBasis->StateDescription[j], rightSpace->FermionBasis->StateLzMax[j], 
				   rightSpace->TemporaryState, rightSpace->TemporaryStateLzMax);
	      int k = 0;
	      for (; k <= rightSpace->TemporaryStateLzMax; ++k)
		this->TemporaryState[k] = leftSpace->TemporaryState[k] + rightSpace->TemporaryState[k];
	      this->TemporaryStateLzMax = rightSpace->TemporaryStateLzMax;
	      if (leftSpace->TemporaryStateLzMax > rightSpace->TemporaryStateLzMax)
		{
		  for (; k <= leftSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k] = leftSpace->TemporaryState[k];
		  this->TemporaryStateLzMax = leftSpace->TemporaryStateLzMax;
		}
	      int TmpPos = this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
	      if (TmpPos < this->HilbertSpaceDimension)
		{
		  Factorial2 = Factorial1;
		  for (k = 0; k <= rightSpace->TemporaryStateLzMax; ++k)
		    if (rightSpace->TemporaryState[k] > 1)
		      Factorial2.FactorialDivide(rightSpace->TemporaryState[k]);
		  for (k = 0; k <= this->TemporaryStateLzMax; ++k)
		    if (this->TemporaryState[k] > 1)
		      Factorial2.FactorialMultiply(this->TemporaryState[k]);	      
		  symmetrizedVector[TmpPos] += sqrt(Factorial2.GetNumericalValue()) * TmpCoefficient * rightVector[j];
		}
	    }
	}     
    }
}

// symmetrized a product of two uncoupled states, using rational input vectors
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// return value = symmetrized state

LongRationalVector BosonOnSphereShort::SymmetrizeU1U1State (LongRationalVector& leftVector, LongRationalVector& rightVector, 
							    BosonOnSphereShort* leftSpace, BosonOnSphereShort* rightSpace, AbstractArchitecture* architecture)
{
  LongRationalVector SymmetrizedVector (this->LargeHilbertSpaceDimension,true);

  FQHESphereSymmetrizeU1U1StateOperation Operation (this, leftSpace, rightSpace, &SymmetrizedVector, &leftVector, &rightVector);
  Operation.ApplyOperation(architecture);

  return SymmetrizedVector;
}


// symmetrized a product of two uncoupled states, using rational input vectors
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// return value = symmetrized state

void BosonOnSphereShort::SymmetrizeU1U1StateCore (LongRationalVector& symmetrizedVector, LongRationalVector& leftVector, LongRationalVector& rightVector, 
						  BosonOnSphereShort* leftSpace, BosonOnSphereShort* rightSpace,
						  unsigned long firstComponent, unsigned long nbrComponents)
{
  long LastComponent = long(firstComponent + nbrComponents);
  
  FactorialCoefficient Factorial1;
  FactorialCoefficient Factorial2;
  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      this->FermionToBoson(leftSpace->FermionBasis->StateDescription[i], leftSpace->FermionBasis->StateLzMax[i], 
			   leftSpace->TemporaryState, leftSpace->TemporaryStateLzMax);
      for (int k = leftSpace->TemporaryStateLzMax + 1;  k <= leftSpace->LzMax; ++k)
	leftSpace->TemporaryState[k] = 0;
      LongRational TmpCoefficient = leftVector[i];
      Factorial1.SetToOne();
      Factorial1.Power2Divide(leftSpace->NbrBosons);
      for (int k = 0; k <= leftSpace->TemporaryStateLzMax; ++k)
 	if (leftSpace->TemporaryState[k] > 1)
 	  Factorial1.FactorialDivide(leftSpace->TemporaryState[k]);
      
      for (long j = 0l; j < rightSpace->LargeHilbertSpaceDimension; ++j)
	{
	  this->FermionToBoson(rightSpace->FermionBasis->StateDescription[j], rightSpace->FermionBasis->StateLzMax[j], 
			       rightSpace->TemporaryState, rightSpace->TemporaryStateLzMax);
	  int k = 0;
	  for (; k <= rightSpace->TemporaryStateLzMax; ++k)
	    this->TemporaryState[k] = leftSpace->TemporaryState[k] + rightSpace->TemporaryState[k];
	  this->TemporaryStateLzMax = rightSpace->TemporaryStateLzMax;
	  if (leftSpace->TemporaryStateLzMax > rightSpace->TemporaryStateLzMax)
	    {
	      for (; k <= leftSpace->TemporaryStateLzMax; ++k)
		this->TemporaryState[k] = leftSpace->TemporaryState[k];
	      this->TemporaryStateLzMax = leftSpace->TemporaryStateLzMax;
	    }
	  int TmpPos = this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
	  if (TmpPos < this->HilbertSpaceDimension)
	    {
 	      Factorial2 = Factorial1;
 	      for (k = 0; k <= rightSpace->TemporaryStateLzMax; ++k)
 		if (rightSpace->TemporaryState[k] > 1)
 		  Factorial2.FactorialDivide(rightSpace->TemporaryState[k]);
 	      for (k = 0; k <= this->TemporaryStateLzMax; ++k)
 		if (this->TemporaryState[k] > 1)
 		  Factorial2.FactorialMultiply(this->TemporaryState[k]);	      
	      symmetrizedVector[TmpPos] += Factorial2.GetLongRationalValue() * TmpCoefficient * rightVector[j];
	    }
	}
    }
}

// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals

bool BosonOnSphereShort::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  int Max = this->TemporaryStateLzMax + 2 - pauliR;
  
  for (int m = 0; m < Max; ++m)
    {
      int Count = TemporaryState[m];
      for (int Offset = 1; Offset < pauliR; ++Offset)
	Count += TemporaryState[m + Offset];
      if (Count > pauliK)
	return false;
    }
  return true;
}

// transform a vector belonging to this vector space in the lz->-lz
//
// finalSpace = the space obtained after the lz->-lz operation
// initialVector = vector on which the operation will be apply
// return value = vector resulting of the operation

RealVector BosonOnSphereShort::GetLzSymmetricVector(ParticleOnSphere* finalSpace, RealVector& initialVector)
{
  BosonOnSphereShort* TmpFinalState = (BosonOnSphereShort*) finalSpace;
  RealVector TmpVector(this->LargeHilbertSpaceDimension, true);
  for(long i = 0 ; i < this->LargeHilbertSpaceDimension ; i++)
    {
      unsigned long Tmp = this->FermionBasis->GetSymmetricState(this->FermionBasis->StateDescription[i]);
      int TmpLzMax = this->FermionBasis->LzMax;
      while (((Tmp >> TmpLzMax) & 0x1ul) == 0x0ul)
	--TmpLzMax;
      TmpVector[TmpFinalState->FermionBasis->FindStateIndex(Tmp, TmpLzMax)] = initialVector[i];
    }
  return TmpVector;
}

// transform a vector belonging to this vector space in the lz->-lz
//
// finalSpace = the space obtained after the lz->-lz operation
// initialVector = vector on which the operation will be apply
// return value = vector resulting of the operation

LongRationalVector BosonOnSphereShort::GetLzSymmetricVector(ParticleOnSphere* finalSpace, LongRationalVector& initialVector)
{
  BosonOnSphereShort* TmpFinalState = (BosonOnSphereShort*) finalSpace;
  LongRationalVector TmpVector(this->LargeHilbertSpaceDimension, true);
  for(long i = 0 ; i < this->LargeHilbertSpaceDimension ; i++)
    {
      unsigned long Tmp = this->FermionBasis->GetSymmetricState(this->FermionBasis->StateDescription[i]);
      int TmpLzMax = this->FermionBasis->LzMax;
      while (((Tmp >> TmpLzMax) & 0x1ul) == 0x0ul)
	--TmpLzMax;
      TmpVector[TmpFinalState->FermionBasis->FindStateIndex(Tmp, TmpLzMax)] = initialVector[i];
    }
  return TmpVector;
}

// Compute the product of two states that belong to the same Hilbert Space
//
// firstState = reference on one of the states whose product will be computed
// secondState = reference on the other state whose product will be computed
// OutputVector = reference on the vector where the result will be stored
// FinalSpace = pointer on the Hilbert Space whose the final state belong

void BosonOnSphereShort::BosonicStateTimeBosonicState(RealVector& firstState, RealVector& secondState, RealVector& outputVector, int minIndex, int nbrComponents, BosonOnSphereShort * finalSpace)
{
  unsigned long * FinalStates = new unsigned long[finalSpace->GetHilbertSpaceDimension()];
  long * Weigth = new long [finalSpace->GetHilbertSpaceDimension()];
  unsigned long * FirstMonomials = new unsigned long[this->NbrBosons];
  unsigned long * SecondMonomials = new unsigned long[this->NbrBosons];
  int MaxIndex = minIndex + nbrComponents;
  
  for (long i = minIndex; i < MaxIndex; i++)
    {
      if(firstState[i] != 0)
	{
	  this->GetMonomial(i,FirstMonomials);
	  this->GetMonomial(i, SecondMonomials);
	  int * EqualPowerIndex = new int[this->NbrBosons];
	  int NbrEqualPower = 0;
	  for (int Index = 0; Index < this->NbrBosons-1; Index++)
	    {
	      if(FirstMonomials[Index] == FirstMonomials[Index+1])
		{
		  EqualPowerIndex[NbrEqualPower] = Index;
		  NbrEqualPower++;
		}
	    }
	  
	  unsigned long NbrStates = this->ProductOfTwoMonomials(FirstMonomials,EqualPowerIndex,NbrEqualPower,SecondMonomials,FinalStates,Weigth,finalSpace);
	  
	  for (long Index = 0l; Index < (long) NbrStates; Index++)
	    {
	      int TmpLzMax = 2 * this->LzMax + this->NbrBosons - 1;
	      while ((FinalStates[Index] >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      
	      outputVector[finalSpace->FermionBasis->FindStateIndex(FinalStates[Index],TmpLzMax)] += firstState[i]*secondState[i]*Weigth[Index];
	    }
	  
	  for (long j = i + 1l; j < this->HilbertSpaceDimension; j++)
	    {
	      if(secondState[j] != 0)
		{
		  this->GetMonomial(j,SecondMonomials);	
		  NbrStates = this->ProductOfTwoMonomials(FirstMonomials,EqualPowerIndex,NbrEqualPower,SecondMonomials,FinalStates,Weigth,finalSpace);
		  for (long Index = 0; Index < (long) NbrStates; Index++)
		    {
		      int TmpLzMax = 2 * this->LzMax + this->NbrBosons - 1;
		      while ((FinalStates[Index] >> TmpLzMax) == 0x0ul)
			--TmpLzMax;
		      outputVector[finalSpace->FermionBasis->FindStateIndex(FinalStates[Index],TmpLzMax)] += (firstState[i]*secondState[j]+firstState[j]*secondState[i])*Weigth[Index];
		    }
		}
	    }
	} 
    }
}


// Compute the product of two states that belong to different Hilbert Spaces
//
// firstState = reference on one of the states whose product will be computed
// secondState = reference on the other state whose product will be computed
// OutputVector = reference on the vector where the result will be stored
// SecondSpace = pointer on the Hilbert Space whose the second state belong
// minIndex = first computed component
// nbrComponents = Nomber of computed components
// FinalSpace = pointer on the Hilbert Space whose the final state belong

void BosonOnSphereShort::BosonicStateTimeBosonicState(RealVector& firstState, RealVector& secondState, RealVector& outputVector,BosonOnSphereShort * secondSpace,int minIndex,int nbrComponents,BosonOnSphereShort * finalSpace)
{
  unsigned long * FinalStates = new unsigned long[finalSpace->GetHilbertSpaceDimension()];
  long * Weigth = new long [finalSpace->GetHilbertSpaceDimension()];
  unsigned long * FirstMonomials = new unsigned long[this->NbrBosons];
  unsigned long * SecondMonomials = new unsigned long[this->NbrBosons];
  int MaxIndex = minIndex + nbrComponents;
  for (long i = minIndex; i < MaxIndex; i++)
    {
      this->GetMonomial(i,FirstMonomials);
      int * EqualPowerIndex = new int[this->NbrBosons];
      int NbrEqualPower = 0;
      for (int Index = 0; Index < this->NbrBosons-1; Index++)
	{
	  if(FirstMonomials[Index] == FirstMonomials[Index+1])
	    {
	      EqualPowerIndex[NbrEqualPower] = Index;
	      NbrEqualPower++;
	    }
	}
      for (long j = 0; j < secondSpace->HilbertSpaceDimension; j++)
	{
	  secondSpace->GetMonomial(j,SecondMonomials);
	  unsigned long NbrStates = this->ProductOfTwoMonomials(FirstMonomials,EqualPowerIndex,NbrEqualPower,SecondMonomials,FinalStates,Weigth,finalSpace);
	  for (unsigned long Index = 0x0ul; Index < NbrStates; Index++)
	    {
	      int TmpLzMax = finalSpace->LzMax + finalSpace->NbrBosons - 1;
	      while ((FinalStates[Index] >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      outputVector[finalSpace->FermionBasis->FindStateIndex(FinalStates[Index],TmpLzMax)]+=firstState[i]*secondState[j]*Weigth[Index];
	    }
	}
    }
}


// Compute the product of two Monomials
//
// firstState = array where the monomial representation of the first state is stored
// secondState = array where the monomial representation of the second state is stored
// finalStates = reference on the array where the fermionic representation of the states product will be stored
// weigth = reference on the array where the coefficient of the states product will be stored
// FinalSpace = pointer on the Hilbert Space whose the final monomials belong

unsigned long BosonOnSphereShort::ProductOfTwoMonomials (unsigned long* firstState,int * equalPowerIndex,const int nbrEqualPower,unsigned long* secondState, unsigned long * & finalStates, long * & weigth,BosonOnSphereShort * finalSpace)
{
  unsigned long NbrStates = 0ul;
  long Coef = 1l;
  unsigned long State [this->NbrBosons];
  
  if(CheckLexiOrder(equalPowerIndex,secondState,nbrEqualPower))
    {
      for (int Index = 0; Index < this->NbrBosons; Index++)
	State[Index] = firstState[Index] + secondState[Index];
	  
      finalStates[0] = finalSpace->ConvertFromMonomial(State);
      weigth[0] = Coef;
      NbrStates++;
    }
  while (std::prev_permutation(secondState,secondState + this->NbrBosons))
    {
      if(CheckLexiOrder(equalPowerIndex,secondState,nbrEqualPower))
	{
	  for (int Index = 0; Index < this->NbrBosons; Index++)
	    State[Index] = firstState[Index] + secondState[Index];
	  Coef = this->ComputeCoefficient(State,firstState);
	  
	  SortArrayDownOrdering(State,this->NbrBosons);
	  NbrStates += SearchInArrayAndSetWeight(finalSpace->ConvertFromMonomial(State),finalStates,weigth,NbrStates,Coef);
	}
    }
  return NbrStates;
}

// Compute the coefficient of a monomial in the decomposition of a product of two monomials
//
// state = array where the monomial representation of the state whose coefficient is to be computed is stored
// firstState = array where the monomial representation of one of the states whose product is computed, is stored
// return value = coefficient of the monomial stored in state


long BosonOnSphereShort::ComputeCoefficient(unsigned long * state,const unsigned long * firstState)
{
 
  int IntegerArray[this->NbrBosons];
  for (int Index = 0; Index < this->NbrBosons; Index++)
    IntegerArray [Index] = Index;

  SortArrayDownOrdering(state,IntegerArray,this->NbrBosons);

  int Index = 1;
  int Indice[this->NbrBosons];
  unsigned long Temp = state[0];
  long Coef = 1;
  Indice[0] = IntegerArray[0];
  long p = 1;
  int SizeIndice = 1;
  while (Index < this->NbrBosons-1)
    {
      if (Temp == state[Index])
	{
	  Indice[SizeIndice] = IntegerArray[Index];
	  SizeIndice++;
	}
      else
	{
	  unsigned long Tableau[SizeIndice];
	  Tableau[0] = firstState [Indice[0]];
	  for (int Index1 = 1; Index1 < SizeIndice; Index1++)
	    {
	      int k = Index1 - 1;
	      bool Compteur = true;
	      while ( (k > -1) && (Compteur) )
		{
		  if(firstState[Indice[Index1]] <= Tableau[k])
		    {
		      Tableau[k+1] = firstState[Indice[Index1]];
		      Compteur = false;
		    }
		  else
		    {
		      Tableau[k+1] = Tableau[k];
		    }
		  k--;
		}
	      if (Compteur == true)
		{
		  Tableau[0] = firstState[Indice[Index1]];
		}
	    }
	  p = 1;
	  while(std::prev_permutation(Tableau,(Tableau + SizeIndice)))
	    p++;
	  Coef *= p;
	  Temp = state[Index];
	  SizeIndice = 1;
	  Indice[0] = IntegerArray[Index];
	}
      Index++;
    }

  if(Temp == state[NbrBosons-1])
    {
      Indice[SizeIndice] = IntegerArray[Index];
      SizeIndice++;
    }

  unsigned long Tableau[SizeIndice];
  Tableau[0] = firstState[Indice[0]];
  for (int Index1 = 1; Index1 < SizeIndice; Index1++)
    {
      int k = Index1 - 1;
      bool Compteur = true;
      while ( (k > -1) && (Compteur) )
	{
	  if(firstState[Indice[Index1]] <= Tableau[k])
	    {
	      Tableau[k+1] = firstState[Indice[Index1]];
	      Compteur = false;
	    }
	  else
	    {
	      Tableau[k+1] = Tableau[k];
	    }
	  k--;
	}
      if (Compteur == true)
	{
	  Tableau[0] = firstState[Indice[Index1]];
	}
    }
  p = 1;
  while(std::prev_permutation(Tableau,Tableau + SizeIndice))
    p++;
  Coef *= p;
  return Coef;
}

// Compute the product of a bosonic state and a fermionic state
//
// firstState = reference on the bosonic state
// secondState = reference on the fermionic state
// outputVector = reference on the vector where the result will be stored
// fermionSpace = pointer on the fermionic Hilbert Space whose the second state belong
// minIndex = first computed component
// nbrComponents = Nomber of computed components
// finalSpace = pointer on the Hilbert Space whose the final state belong

void BosonOnSphereShort::BosonicStateTimeFermionicState(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, FermionOnSphere * fermionSpace, int minIndex, int nbrComponents, FermionOnSphere * finalSpace)
{
  unsigned long * FinalStates = new unsigned long[finalSpace->GetHilbertSpaceDimension()];
  long * Weigth = new long [finalSpace->GetHilbertSpaceDimension()];
  unsigned long * Monomial = new unsigned long[this->NbrBosons];
  unsigned long * Slater = new unsigned long[this->NbrBosons];
  int MaxIndex = minIndex + nbrComponents;
  for (long i = minIndex; i < MaxIndex; i++)
    {
      if(bosonState[i] != 0)
	{
	  this->GetMonomial(i,Monomial);
	  for (long j = 0; j < fermionSpace->HilbertSpaceDimension; j++)
	    {
	      if(fermionState[j] != 0)
		{
		  fermionSpace->GetMonomial(j,Slater);
		  unsigned long NbrStates = this->MonomialsTimesSlater(Slater,Monomial,FinalStates,Weigth,finalSpace);
		  for (unsigned long Index = 0; Index < NbrStates; Index++)
		    {
		      int TmpLzMax = finalSpace->LzMax;
		      while ((FinalStates[Index] >> TmpLzMax) == 0x0ul)
			--TmpLzMax;
		      outputVector[finalSpace->FindStateIndex(FinalStates[Index],TmpLzMax)] += bosonState[i]*fermionState[j]*Weigth[Index];
		    }
		}
	    }
	}
    }
}

// Compute the product of a Monomial and a Slater determinant
//
// firstState = array where the monomial representation of the first state is stored
// secondState = array where the monomial representation of the second state is stored
// finalStates = reference on the array where the fermionic representation of the states product will be stored
// weigth = reference on the array where the coefficient of the states product will be stored
// finalSpace = pointer to the destination Hilbert space
// return value = number of different states generated by the product

unsigned long BosonOnSphereShort::MonomialsTimesSlater (unsigned long* slater, unsigned long* monomial, unsigned long * & finalStates, long * & weigth, FermionOnSphere * finalSpace)
{
  unsigned long NbrStates = 0;
  long Coef = 1;
  unsigned long State[this->NbrBosons];
  bool Bool = true;
  unsigned long Mask;
  unsigned long Sign = 0ul;
  unsigned long TmpState = 0;
  for (int Index = 0; Index < this->NbrBosons; Index++)
    State[Index] = slater[Index] + monomial[Index];
  
  finalStates[0] = finalSpace->ConvertFromMonomial(State);
  weigth[0] = Coef;
  NbrStates++;
  while (std::prev_permutation(monomial,monomial+this->NbrBosons))
    {
      for (int Index = 0; Index < this->NbrBosons; Index++)
	State[Index] = slater[Index] + monomial[Index];
      
      Bool = true;
      TmpState = 0ul;
      Sign = 0ul;
      for(int i = 0; (i < this->NbrBosons) && (Bool == true); i++)
	{
	  Mask = (1ul << State[i]);
	  if ( (TmpState & Mask) != 0ul)
	    Bool = false;
	  unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef _64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
#endif
	  TmpState2 ^= TmpState2 >> 16;
	  TmpState2 ^= TmpState2 >> 8;
	  TmpState2 ^= TmpState2 >> 4;
	  TmpState2 ^= TmpState2 >> 2;
	  TmpState2 ^= TmpState2 >> 1;
	  Sign ^= TmpState2;
	  TmpState |= Mask;
	}
      if(Bool == true)
	{
	  if ((Sign & 0x1ul) == 0ul)
	    {
	      Coef = 1;
	    }
	  else
	    {
	      Coef = -1;
	    }
	  NbrStates += SearchInArrayAndSetWeight(TmpState,finalStates,weigth,NbrStates,Coef);
	}
    }
  return NbrStates;
}

// Compute the product of two fermionic states that belong to different Hilbert Spaces
//
// firstState = reference on one of the states whose product will be computed
// secondState = reference on the other state whose product will be computed
// OutputVector = reference on the vector where the result will be stored
// fermionSpace1 = pointer on the Hilbert Space whose the first state belong to
// fermionSpace2 = pointer on the Hilbert Space whose the second state belong to 
// minIndex = first computed component
// nbrComponents = Nomber of computed components


void BosonOnSphereShort::FermionicStateTimeFermionicState(RealVector& fermionState1, RealVector& fermionState2, RealVector& outputVector, FermionOnSphere * fermionSpace1, FermionOnSphere * fermionSpace2, int minIndex, int nbrComponents)
{
  unsigned long * FinalStates = new unsigned long [this->GetHilbertSpaceDimension()];
  long * Weigth = new long [this->GetHilbertSpaceDimension()];
  unsigned long * Slater1 = new unsigned long[this->NbrBosons];
  unsigned long * Slater2 = new unsigned long[this->NbrBosons];
  int MaxIndex = minIndex + nbrComponents;
  FactorialCoefficient Coefficient;	
  for (long i = minIndex; i < MaxIndex; i++)
    {
      if(fermionState1[i] != 0)
	{
	  fermionSpace1->GetMonomial(i,Slater1);
	  for (long j = 0; j < fermionSpace2->HilbertSpaceDimension; j++)
	    {
	      if(fermionState2[j] != 0)
		{
		  fermionSpace2->GetMonomial(j,Slater2);
		  unsigned long NbrStates = this->SlaterTimesSlater(Slater1,Slater2,FinalStates,Weigth);
		  for (unsigned long Index = 0; Index < NbrStates; Index++)
		    {
		      int TmpLzMax = this->LzMax+this->NbrBosons - 1;
		      while ((FinalStates[Index] >> TmpLzMax) == 0x0ul)
			--TmpLzMax;
		      this->FermionToBoson(FinalStates[Index],TmpLzMax,this->TemporaryState,this->TemporaryStateLzMax);
		      Coefficient.SetToOne();
		      for(int p = 0; p < this->TemporaryStateLzMax + 1; p++)
			{
			  Coefficient.FactorialMultiply(this->TemporaryState[p]);
			}
		      outputVector[this->FermionBasis->FindStateIndex(FinalStates[Index],TmpLzMax)] += fermionState1[i]*fermionState2[j]*Coefficient.GetIntegerValue()*Weigth[Index];
		    }
		}
	    }
	}
    }
}

// Compute the product of two Slater determinants
//
// slater1 = array where the monomial representation of the first slater is stored
// slater2 = array where the monomial representation of the second slater is stored
// finalStates = reference on the array where the fermionic representation of the states product will be stored
// weigth = reference on the array where the coefficient of the states product will be stored

unsigned long BosonOnSphereShort::SlaterTimesSlater (unsigned long* slater1,unsigned long* slater2, unsigned long * & finalStates, long * & weigth)
{
  unsigned long NbrStates = 0;
  long Coef = 1;
  unsigned long State [this->NbrBosons];
  unsigned long TmpState = 0ul;
  unsigned long Sign = 0ul;
  unsigned long Mask = 0ul;
  
  for (int Index = 0; Index < this->NbrBosons; Index++)
    {
      State[Index] = slater1[Index] + slater2[Index];
    }
  
  finalStates[0] = this->ConvertFromMonomial(State);
  weigth[0] = Coef;
  NbrStates++;
  
  while (std::prev_permutation(slater2,slater2+this->NbrBosons))
    {
      for (int Index = 0; Index < this->NbrBosons; Index++)
	{
	  State[Index] = slater1[Index] + slater2[Index];
	}
      
      TmpState = 0ul;
      Sign = 0ul;
      for (int i = 0; i < this->NbrBosons ; i++)
	{
	  Mask = (1ul << slater2[i]);
	  unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
#endif	
	  TmpState2 ^= TmpState2 >> 16;
	  TmpState2 ^= TmpState2 >> 8;
	  TmpState2 ^= TmpState2 >> 4;
	  TmpState2 ^= TmpState2 >> 2;
	  TmpState2 ^= TmpState2 >> 1;
	  Sign ^= TmpState2;
	  TmpState |= Mask;
	}
      SortArrayDownOrdering(State,this->NbrBosons);
      
      if ((Sign & 0x1ul) == 0ul)
	{
	  Coef = 1;
	}
      else
	{
	  Coef = -1;
	}
      NbrStates += SearchInArrayAndSetWeight(this->ConvertFromMonomial(State),finalStates,weigth,NbrStates,Coef);
    }
  return NbrStates;
}

// Compute the product of two fermionic states that belong to the same Hilbert Spaces
//
// fermionState1 = reference on one of the states whose product will be computed
// fermionState2 = reference on the other state whose product will be computed
// OutputVector = reference on the vector where the result will be stored
// fermionSpace = pointer on the fermionic Hilbert Space
// minIndex = first computed component
// nbrComponents = Nomber of computed components

void BosonOnSphereShort::FermionicStateTimeFermionicState(RealVector& fermionState1, RealVector& fermionState2, RealVector& outputVector, FermionOnSphere * fermionSpace, int minIndex, int nbrComponents)
{
  unsigned long * FinalStates= new unsigned long[this->GetHilbertSpaceDimension()];
  long * Weigth = new long [this->GetHilbertSpaceDimension()];
  unsigned long * Slater1 = new unsigned long[this->NbrBosons];
  unsigned long * Slater2 = new unsigned long[this->NbrBosons];
  int MaxIndex = minIndex + nbrComponents;
  FactorialCoefficient Coefficient;	
  for (long i = minIndex; i < MaxIndex; i++)
    {
      if(fermionState1[i] != 0)
	{
	  fermionSpace->GetMonomial(i,Slater1);
	  for (long j = i; j < fermionSpace->HilbertSpaceDimension; j++)
	    {
	      if(fermionState2[j] != 0)
		{
		  fermionSpace->GetMonomial(j,Slater2);
		  unsigned long Limit = this->SlaterTimesSlater(Slater1,Slater2,FinalStates,Weigth);
		  for (unsigned long Index = 0; Index < Limit; Index++)
		    {
		      int TmpLzMax = this->LzMax + this->NbrBosons - 1;
		      while ((FinalStates[Index] >> TmpLzMax) == 0x0ul)
			--TmpLzMax;
		      this->FermionToBoson(FinalStates[Index],TmpLzMax,this->TemporaryState,this->TemporaryStateLzMax);
		      Coefficient.SetToOne();
		      for(int p = 0; p < this->TemporaryStateLzMax + 1; p++)
			{
			  Coefficient.FactorialMultiply(this->TemporaryState[p]);
			}
		      if(i == j)
			outputVector[this->FermionBasis->FindStateIndex(FinalStates[Index],TmpLzMax)] += fermionState1[i]*fermionState2[j]*Coefficient.GetIntegerValue()*Weigth[Index];
		      else
			outputVector[this->FermionBasis->FindStateIndex(FinalStates[Index],TmpLzMax)] += (fermionState1[i]*fermionState2[j]+fermionState1[j]*fermionState2[i])*Coefficient.GetIntegerValue()*Weigth[Index];
		    }
		}
	    }
	}
    }
}

// compute the projection of the product of a bosonic state and the halperin 110 state
//
// bosonState = real vector where the bosonic state is stored
// outputVector = real vector where the result has to be stored
// fermionSpace = pointer to the fermionic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void BosonOnSphereShort::BosonicStateTimePolarizedSlaters(RealVector& bosonState, RealVector& outputVector, FermionOnSphere * fermionSpace , FermionOnSphereWithSpin* finalSpace, int firstComponent,int nbrComponent)
{
  map<unsigned long , double> SortingMap;
  map<unsigned long , double>::iterator It;
  
  BinomialCoefficients Binomial(this->NbrBosons);
  int NbrParticlesPerColor = this->NbrBosons >> 1;
  unsigned long NbrPermutations = Binomial(this->NbrBosons, NbrParticlesPerColor);
  unsigned long* Permutations1 = new unsigned long[NbrPermutations];
  unsigned long* Permutations2 = new unsigned long[NbrPermutations];
  
  EvaluatePermutationsOfSubGroups(NbrPermutations,this->NbrBosons, NbrParticlesPerColor, Permutations1, Permutations2);
  
  unsigned long* Monomial = new unsigned long[this->NbrBosons];
  unsigned long* Slater = new unsigned long[fermionSpace->NbrFermions];
  
  int NbrMax = firstComponent + nbrComponent;
  //  int NbrVariable = 0;
  
  fermionSpace->ConvertToMonomial(fermionSpace->StateDescription[0], Slater);
  
  for (int j = firstComponent; j < NbrMax; j++)
    {
      if(bosonState[j] != 0)
	{
	  this->GetMonomial(j, Monomial);
	  finalSpace->MonomialsTimesPolarizedSlater(Slater, Monomial, SortingMap,NbrPermutations,Permutations1, Permutations2, bosonState[j]);
	}
    }
  
  for ( It = SortingMap.begin() ; It != SortingMap.end(); It++)
    {
      int TmpLzMax = 2 * finalSpace->LzMax + 1;
      while ((( (*It).first >> TmpLzMax) & 0x1ul) == 0x0ul)
	--TmpLzMax;
      outputVector[finalSpace->FindStateIndex((*It).first, TmpLzMax)] += (*It).second;
    }
  
  delete [] Monomial;
  delete [] Slater;
}


void BosonOnSphereShort::EvaluatePartialDensityMatrixMultipartiteParticlePartition( ParticleOnSphere * spaceA, ParticleOnSphere * spaceB, ParticleOnSphere * spaceC,  RealVector groundState, RealSymmetricMatrix* densityMatrix, AbstractArchitecture* architecture)
{
  BosonOnSphereFullShort * SpaceA = ((BosonOnSphereFullShort *) spaceA); 
  BosonOnSphereFullShort * SpaceB = ((BosonOnSphereFullShort *) spaceB);
  BosonOnSphereShort * SpaceC = ((BosonOnSphereShort *) spaceC);
  unsigned long* TmpMonomialA = new unsigned long [SpaceA->NbrBosons];
  unsigned long* TmpMonomialB = new unsigned long [SpaceB->NbrBosons];
  unsigned long* TmpMonomialC = new unsigned long [SpaceC->NbrBosons]; 
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  
  long TensorStateDimension = SpaceA->HilbertSpaceDimension *  SpaceB->HilbertSpaceDimension;
  
  int* TmpStatePositionFinal = new int [TensorStateDimension];
  int* TmpStatePositionA = new int [TensorStateDimension];
  int* TmpStatePositionB = new int [TensorStateDimension];
  double* TmpStateCoefficient = new double [TensorStateDimension];
  long TmpNbrNonZeroElements = 0l;
  long MinIndex = 0l;
  long MaxIndex =  SpaceC->HilbertSpaceDimension;
  double* LogFactorials = new double[this->NbrBosons + 1];
   LogFactorials[0] = 0.0;
   LogFactorials[1] = 0.0;
   for (int i = 2 ; i <= this->NbrBosons; ++i)
     LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
   
   double* TmpLogFactorialsA = new double [SpaceA->HilbertSpaceDimension];
   double* TmpLogFactorialsB = new double [SpaceB->HilbertSpaceDimension];
   double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[SpaceA->NbrBosons] - LogFactorials[SpaceB->NbrBosons] - LogFactorials[SpaceC->NbrBosons];

  for (int i = 0; i < SpaceA->HilbertSpaceDimension; ++i)
    {
      SpaceA->FermionToBoson(SpaceA->FermionBasis->StateDescription[i], SpaceA->FermionBasis->StateLzMax[i], SpaceA->TemporaryState, SpaceA->TemporaryStateLzMax);
      double TmpFactor = 0.0;
      for (int k = 0; k <= SpaceA->TemporaryStateLzMax; ++k)
	TmpFactor += LogFactorials[SpaceA->TemporaryState[k]];
      TmpLogFactorialsA[i] =  TmpFactor;
    }
    
    for (int i = 0; i < SpaceB->HilbertSpaceDimension; ++i)
    {
      SpaceB->FermionToBoson(SpaceB->FermionBasis->StateDescription[i], SpaceB->FermionBasis->StateLzMax[i], SpaceB->TemporaryState, SpaceB->TemporaryStateLzMax);
      double TmpFactor = 0.0;
      for (int k = 0; k <= SpaceB->TemporaryStateLzMax; ++k)
	TmpFactor += LogFactorials[SpaceB->TemporaryState[k]];
      TmpLogFactorialsB[i] =  TmpFactor;
    }
    
    
    
   for (; MinIndex < MaxIndex; ++MinIndex)    
     {
      int Pos = 0;
      SpaceC->ConvertToMonomial(SpaceC->FermionBasis->StateDescription[MinIndex], SpaceC->FermionBasis->StateLzMax[MinIndex], TmpMonomialC);
      SpaceC->FermionToBoson(SpaceC->FermionBasis->StateDescription[MinIndex], SpaceC->FermionBasis->StateLzMax[MinIndex], SpaceC->TemporaryState, SpaceC->TemporaryStateLzMax);
      
      double TmpSpaceCFactorial = 0.0;
      for (int k = 0; k <= SpaceC->TemporaryStateLzMax; ++k)
	TmpSpaceCFactorial += LogFactorials[SpaceC->TemporaryState[k]];
      
      //cout <<"A space dimension = " <<SpaceA->HilbertSpaceDimension<<endl;
      for (int IndexSpaceA = 0; IndexSpaceA < SpaceA->HilbertSpaceDimension; IndexSpaceA++)
	{
	  SpaceA->ConvertToMonomial(SpaceA->FermionBasis->StateDescription[IndexSpaceA], SpaceA->FermionBasis->StateLzMax[IndexSpaceA], TmpMonomialA);
	//cout <<"B space dimension = " <<SpaceB->HilbertSpaceDimension<<endl;
	  for (int IndexSpaceB = 0; IndexSpaceB < SpaceB->HilbertSpaceDimension; IndexSpaceB++)
	{
	  SpaceB->ConvertToMonomial(SpaceB->FermionBasis->StateDescription[IndexSpaceB], SpaceB->FermionBasis->StateLzMax[IndexSpaceB], TmpMonomialB);
	  //cout <<" SpaceA->GetLzValue(IndexSpaceA) = " <<SpaceA->GetLzValue(IndexSpaceA)<<" SpaceB->GetLzValue(IndexSpaceB)"<< SpaceB->GetLzValue(IndexSpaceB)<<endl;
	  //cout <<"this->TotalLz = "<<this->TotalLz<<endl;
	  //cout <<" SpaceC->TotalLz "<< SpaceC->TotalLz<<endl;
	  if ( SpaceA->GetLzValue(IndexSpaceA) + SpaceB->GetLzValue(IndexSpaceB) == this->TotalLz - SpaceC->TotalLz)
	  {
	  int TmpIndexA = 0;
	  int TmpIndexB = 0;
	  int TmpIndexC = 0;
	  int TmpIndexTotal = 0;
	  
	  /*cout <<"A = ";
	  for (int i = 0; i< SpaceA->NbrBosons; i++)
	    {
	     cout << TmpMonomialA[i]<<" ";
	    }
	    cout << "Lz = "<<SpaceA->GetLzValue(IndexSpaceA)<<endl;
	  cout <<"B = ";
	    for (int i = 0; i< SpaceB->NbrBosons; i++)
	    {
	     cout << TmpMonomialB[i]<<" ";
	    }
	    cout << "Lz = "<<SpaceB->GetLzValue(IndexSpaceB)<<endl;
	    cout <<"C = ";
	    for (int i = 0; i< SpaceC->NbrBosons; i++)
	    {
	     cout << TmpMonomialC[i]<<" ";
	    }
	    cout <<endl;*/
	    
	  while ((TmpIndexA < SpaceA->NbrBosons) && (TmpIndexB < SpaceB->NbrBosons) && (TmpIndexC < SpaceC->NbrBosons)) 
	    {
	      while ( (TmpIndexA < SpaceA->NbrBosons) && (TmpMonomialB[TmpIndexB] <= TmpMonomialA[TmpIndexA]) && (TmpMonomialC[TmpIndexC] <= TmpMonomialA[TmpIndexA]))
		{
		  TmpMonomial[TmpIndexTotal] = TmpMonomialA[TmpIndexA];
		  ++TmpIndexA;
		  ++TmpIndexTotal;		  
		}
		
	      if (TmpIndexA < SpaceA->NbrBosons)
		{
		  while ((TmpIndexB < SpaceB->NbrBosons) && (TmpMonomialA[TmpIndexA] <= TmpMonomialB[TmpIndexB])&& (TmpMonomialC[TmpIndexC] <= TmpMonomialB[TmpIndexB]))
		    {
		      TmpMonomial[TmpIndexTotal] = TmpMonomialB[TmpIndexB];
		      ++TmpIndexB;
		      ++TmpIndexTotal;		  
		    }
		}
		
		if ((TmpIndexA < SpaceA->NbrBosons)&&(TmpIndexB < SpaceB->NbrBosons))
		{
		  while ((TmpIndexC < SpaceC->NbrBosons) && (TmpMonomialA[TmpIndexA] <= TmpMonomialC[TmpIndexC])&& (TmpMonomialB[TmpIndexB] <= TmpMonomialC[TmpIndexC]))
		    {
		      TmpMonomial[TmpIndexTotal] = TmpMonomialC[TmpIndexC];
		      ++TmpIndexC;
		      ++TmpIndexTotal;		  
		    }
		}
	    }
	    
	  while ((TmpIndexA < SpaceA->NbrBosons) && (TmpIndexB < SpaceB->NbrBosons))
	    {
	      while ( (TmpIndexA < SpaceA->NbrBosons) && (TmpMonomialB[TmpIndexB] <= TmpMonomialA[TmpIndexA]))
		{
		  TmpMonomial[TmpIndexTotal] = TmpMonomialA[TmpIndexA];
		  ++TmpIndexA;
		  ++TmpIndexTotal;		  
		}
		if (TmpIndexA < SpaceA->NbrBosons)
		{
		  while ((TmpIndexB < SpaceB->NbrBosons) && (TmpMonomialA[TmpIndexA] <= TmpMonomialB[TmpIndexB]))
		    {
		      TmpMonomial[TmpIndexTotal] = TmpMonomialB[TmpIndexB];
		      ++TmpIndexB;
		      ++TmpIndexTotal;		  
		    }
		}
	    }
	    
	    while ((TmpIndexA < SpaceA->NbrBosons) && (TmpIndexC < SpaceC->NbrBosons))
	    {
	      while ( (TmpIndexA < SpaceA->NbrBosons) && (TmpMonomialC[TmpIndexC] <= TmpMonomialA[TmpIndexA]))
		{
		  TmpMonomial[TmpIndexTotal] = TmpMonomialA[TmpIndexA];
		  ++TmpIndexA;
		  ++TmpIndexTotal;		  
		}
		if (TmpIndexA < SpaceA->NbrBosons)
		{
		  while ((TmpIndexC < SpaceC->NbrBosons) && (TmpMonomialA[TmpIndexA] <= TmpMonomialC[TmpIndexC]))
		    {
		      TmpMonomial[TmpIndexTotal] = TmpMonomialC[TmpIndexC];
		      ++TmpIndexC;
		      ++TmpIndexTotal;		  
		    }
		}
	    }
	    
	    while ((TmpIndexB < SpaceB->NbrBosons) && (TmpIndexC < SpaceC->NbrBosons))
	    {
	      while ( (TmpIndexB < SpaceB->NbrBosons) && (TmpMonomialC[TmpIndexC] <= TmpMonomialB[TmpIndexB]))
		{
		  TmpMonomial[TmpIndexTotal] = TmpMonomialB[TmpIndexB];
		  ++TmpIndexB;
		  ++TmpIndexTotal;		  
		}
		if (TmpIndexB < SpaceB->NbrBosons)
		{
		  while ((TmpIndexC < SpaceC->NbrBosons) && (TmpMonomialB[TmpIndexB] <= TmpMonomialC[TmpIndexC]))
		    {
		      TmpMonomial[TmpIndexTotal] = TmpMonomialC[TmpIndexC];
		      ++TmpIndexC;
		      ++TmpIndexTotal;		  
		    }
		}
	    }
	    
	  while (TmpIndexA < SpaceA->NbrBosons)
	    {
	      TmpMonomial[TmpIndexTotal] = TmpMonomialA[TmpIndexA];
	      ++TmpIndexA;
	      ++TmpIndexTotal;		  
	    }
	    
	    while (TmpIndexB < SpaceB->NbrBosons)
	    {
	      TmpMonomial[TmpIndexTotal] = TmpMonomialB[TmpIndexB];
	      ++TmpIndexB;
	      ++TmpIndexTotal;		  
	    }
	    while (TmpIndexC < SpaceC->NbrBosons)
	    {
	      TmpMonomial[TmpIndexTotal] = TmpMonomialC[TmpIndexC];
	      ++TmpIndexC;
	      ++TmpIndexTotal;		  
	    }
	    /*
	    for (int i = 0; i< this->NbrBosons; i++)
	    {
	     cout << TmpMonomial[i]<<" ";
	    }
	    cout <<endl;*/
	    
	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial[0] + this->NbrBosons - 1);
	  if(TmpPos == this->HilbertSpaceDimension)
	  {
	    cout <<"Danger"<<endl;
	    exit(-1);
	  }
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
	      TmpFactorial -= TmpSpaceCFactorial + TmpLogFactorialsA[IndexSpaceA] + TmpLogFactorialsB[IndexSpaceB]+ TmpLogBinomial;
	      TmpFactorial *= 0.5; 
	      TmpStatePositionFinal[Pos] = TmpPos;
	      TmpStatePositionA[Pos] = IndexSpaceA;
	      TmpStatePositionB[Pos] = IndexSpaceB;
	      TmpStateCoefficient[Pos] = exp(TmpFactorial);
	     // cout << TmpPos<<" "<<IndexSpaceA<<" "<<IndexSpaceB<<" "<<TmpStateCoefficient[Pos]<<endl;
	      ++Pos;
	    }
  
  
	  }
	}
	}
	if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int PosA = TmpStatePositionA[j];
	      int PosB = TmpStatePositionB[j];
	      int TensorProductPosL = PosB * SpaceA->HilbertSpaceDimension + PosA;
	      double TmpValue = groundState[TmpStatePositionFinal[j]] * TmpStateCoefficient[j];
	      
	      for (int k = 0; k < Pos; ++k)
	      {
		int TensorProductPosC = TmpStatePositionB[k] * SpaceA->HilbertSpaceDimension + TmpStatePositionA[k];
		if (TensorProductPosC >= TensorProductPosL)
		  {
		    densityMatrix->AddToMatrixElement(TensorProductPosL, TensorProductPosC, TmpValue * groundState[TmpStatePositionFinal[k]] * TmpStateCoefficient[k]);
		    ++TmpNbrNonZeroElements;
		  }
	      }
	    }
	}
     }
  delete [] TmpLogFactorialsA;
  delete [] TmpLogFactorialsB;
  delete [] TmpMonomialA;
  delete [] TmpMonomialB;
  delete [] TmpMonomialC;
  delete [] TmpMonomial;
  delete [] TmpStatePositionFinal; 
  delete [] TmpStatePositionA;
  delete [] TmpStatePositionB;
  delete [] TmpStateCoefficient;
}

// symmetrize a vector with even number of orbitals 
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector to be symmetrized
// leftSpace = pointer to the Hilbert space
// return value = symmetrized state
RealVector BosonOnSphereShort::SymmetrizeU1U1SingleState (RealVector& leftVector, BosonOnSphereShort* leftSpace, bool unnormalizedBasisFlag, AbstractArchitecture* architecture)
{
  RealVector SymmetrizedVector (this->LargeHilbertSpaceDimension,true);
  unsigned long firstComponent = 0;
  unsigned long nbrComponent = leftSpace->GetHilbertSpaceDimension();
  this->SymmetrizeU1U1SingleStateOneInTwoCore (SymmetrizedVector ,leftVector , leftSpace, unnormalizedBasisFlag, firstComponent, nbrComponent);
  /*
  cout << this->TotalLz << " " << SymmetrizedVector.Norm() << endl;
  */
  
  return SymmetrizedVector;
}


// symmetrize a vector with even number of orbitals and rational coefficients
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector to be symmetrized
// leftSpace = pointer to the Hilbert space
// return value = symmetrized state
LongRationalVector BosonOnSphereShort::SymmetrizeU1U1SingleState (LongRationalVector& leftVector, BosonOnSphereShort* leftSpace, bool unnormalizedBasisFlag, AbstractArchitecture* architecture)
{
  LongRationalVector SymmetrizedVector (this->LargeHilbertSpaceDimension,true);
  unsigned long firstComponent = 0;
  unsigned long nbrComponent = leftSpace->GetHilbertSpaceDimension();
  this->SymmetrizeU1U1SingleStateOneInTwoCore (SymmetrizedVector ,leftVector , leftSpace, unnormalizedBasisFlag, firstComponent, nbrComponent);
  
  
  return SymmetrizedVector;
}

// symmetrize a vector with even number of orbitals
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector associated to the first color
// leftSpace = pointer to the Hilbert space of the first color
// return value = symmetrized state
void BosonOnSphereShort::SymmetrizeU1U1SingleStateOneInTwoCore (RealVector& symmetrizedVector, RealVector& leftVector, BosonOnSphereShort* leftSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents)
{
  long LastComponent = long(firstComponent + nbrComponents);
  FactorialCoefficient Factorial1;
  FactorialCoefficient Factorial2;
  if (unnormalizedBasisFlag == false)
  {
   for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      this->FermionToBoson(leftSpace->FermionBasis->StateDescription[i], leftSpace->FermionBasis->StateLzMax[i], 
			       leftSpace->TemporaryState, leftSpace->TemporaryStateLzMax);
	  
      for (int k = leftSpace->TemporaryStateLzMax + 1;  k < leftSpace->NbrLzValue; ++k)
	leftSpace->TemporaryState[k] = 0;
      double TmpCoefficient = leftVector[i];
      Factorial1.SetToOne();
//       Factorial1.Power2Divide(leftSpace->NbrBosons);
      
      for (int k = 0; k <= leftSpace->TemporaryStateLzMax; ++k)
	if (leftSpace->TemporaryState[k] > 1)
	  Factorial1.FactorialDivide(leftSpace->TemporaryState[k]);
//       cout << Factorial1.GetNumericalValue() << endl;
      int k = 0;
      int kyTot = 0;
      for (; k < this->NbrLzValue; ++k)
      {
	this->TemporaryState[k] = leftSpace->TemporaryState[2*k] + leftSpace->TemporaryState[2*k + 1];
	for (int j = 0; j < leftSpace->TemporaryState[2*k]; ++j)
	{
	  Factorial1.FactorialDivide(this->LzMax);
	  Factorial1.FactorialMultiply(leftSpace->LzMax);
	  Factorial1.FactorialDivide(2*k);
	  Factorial1.FactorialMultiply(k);
	  Factorial1.FactorialDivide(leftSpace->LzMax - 2*k);
	  Factorial1.FactorialMultiply(this->LzMax - k);
	}
	for (int j = 0; j < leftSpace->TemporaryState[2*k + 1]; ++j)
	{
	  Factorial1.FactorialDivide(this->LzMax);
	  Factorial1.FactorialMultiply(leftSpace->LzMax);
	  Factorial1.FactorialDivide(2*k + 1);
	  Factorial1.FactorialMultiply(k);
	  Factorial1.FactorialDivide(leftSpace->LzMax - 2*k - 1);
	  Factorial1.FactorialMultiply(this->LzMax - k);
	}
	
	kyTot += this->TemporaryState[k]*k;
	if (leftSpace->TemporaryState[2*k] != 0 || leftSpace->TemporaryState[2*k + 1] != 0)
	  this->TemporaryStateLzMax = k;
      }
	  
      int TmpPos;
      if (kyTot != this->ShiftedTotalLz)
	TmpPos = this->HilbertSpaceDimension;
      else
      {
	TmpPos = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
      }
      if (TmpPos < this->HilbertSpaceDimension)
      {
	for (k = 0; k <= this->TemporaryStateLzMax; ++k)
	  if (this->TemporaryState[k] > 1)
	    Factorial1.FactorialMultiply(this->TemporaryState[k]);	
// 	cout << Factorial1.GetNumericalValue() << endl;
	symmetrizedVector[TmpPos] += sqrt(Factorial1.GetNumericalValue()) * TmpCoefficient;
// 	this->PrintState(cout, TmpPos);
// 	cout << " " << sqrt(Factorial1.GetNumericalValue()) * TmpCoefficient << endl;
      }
    }
  }
  else
  {
    for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      this->FermionToBoson(leftSpace->FermionBasis->StateDescription[i], leftSpace->FermionBasis->StateLzMax[i], 
			       leftSpace->TemporaryState, leftSpace->TemporaryStateLzMax);
	  
      for (int k = leftSpace->TemporaryStateLzMax + 1;  k < leftSpace->NbrLzValue; ++k)
	leftSpace->TemporaryState[k] = 0;
      double TmpCoefficient = leftVector[i];
      Factorial1.SetToOne();
  
      for (int k = 0; k <= leftSpace->TemporaryStateLzMax; ++k)
	if (leftSpace->TemporaryState[k] > 1)
	  Factorial1.FactorialDivide(leftSpace->TemporaryState[k]);
	  
      int k = 0;
      int kyTot = 0;
      for (; k < this->NbrLzValue; ++k)
      {
	this->TemporaryState[k] = leftSpace->TemporaryState[2*k] + leftSpace->TemporaryState[2*k + 1];
	kyTot += this->TemporaryState[k]*k;
	if (leftSpace->TemporaryState[2*k] != 0 || leftSpace->TemporaryState[2*k + 1] != 0)
	  this->TemporaryStateLzMax = k;
      }
	  
      int TmpPos;
      if (kyTot != this->ShiftedTotalLz)
	TmpPos = this->HilbertSpaceDimension;
      else
      {
	TmpPos = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
      }
      if (TmpPos < this->HilbertSpaceDimension)
      {
	for (k = 0; k <= this->TemporaryStateLzMax; ++k)
	  if (this->TemporaryState[k] > 1)
	    Factorial1.FactorialMultiply(this->TemporaryState[k]);	
		    
	symmetrizedVector[TmpPos] += Factorial1.GetNumericalValue() * TmpCoefficient;
      }
    }
  }
}

// symmetrize a vector with even number of orbitals and rational coefficients
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector associated to the first color
// leftSpace = pointer to the Hilbert space of the first color
// return value = symmetrized state
void BosonOnSphereShort::SymmetrizeU1U1SingleStateOneInTwoCore (LongRationalVector& symmetrizedVector, LongRationalVector& leftVector, BosonOnSphereShort* leftSpace, bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents)
{
  long LastComponent = long(firstComponent + nbrComponents);
  
  FactorialCoefficient Factorial1;
  FactorialCoefficient Factorial2;
  if (unnormalizedBasisFlag == false)
    {
      cout << " Error: normalized basis not implemented for rational numbers" << endl;
    }
  else
    {
      for (long i = (long) firstComponent; i < (long) LastComponent; ++i)
	{
	  this->FermionToBoson(leftSpace->FermionBasis->StateDescription[i], leftSpace->FermionBasis->StateLzMax[i], 
			       leftSpace->TemporaryState, leftSpace->TemporaryStateLzMax);
	  
	  for (int k = leftSpace->TemporaryStateLzMax + 1;  k < leftSpace->NbrLzValue; ++k)
	    leftSpace->TemporaryState[k] = 0;
	  LongRational TmpCoefficient = leftVector[i];
// 	  cout << TmpCoefficient << endl;
	  Factorial1.SetToOne();
	  
	  
	  for (int k = 0; k <= leftSpace->TemporaryStateLzMax; ++k)
	    if (leftSpace->TemporaryState[k] > 1)
	      Factorial1.FactorialDivide(leftSpace->TemporaryState[k]);
	  
	  int k = 0;
	  int kyTot = 0;
	  for (; k < this->NbrLzValue; ++k)
	    {
	      this->TemporaryState[k] = leftSpace->TemporaryState[2*k] + leftSpace->TemporaryState[2*k + 1];
	      kyTot += this->TemporaryState[k]*k;
	      if (leftSpace->TemporaryState[2*k] != 0 || leftSpace->TemporaryState[2*k + 1] != 0)
		this->TemporaryStateLzMax = k;
	    }
	  
	  int TmpPos;
	  if (kyTot != this->ShiftedTotalLz)
	    TmpPos = this->HilbertSpaceDimension;
	  else
	    {
	      TmpPos = this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
	    }
	  if (TmpPos < this->HilbertSpaceDimension)
	    {
	      for (k = 0; k <= this->TemporaryStateLzMax; ++k)
		if (this->TemporaryState[k] > 1)
		  Factorial1.FactorialMultiply(this->TemporaryState[k]);
	      
	      symmetrizedVector[TmpPos] += Factorial1.GetLongRationalValue() * TmpCoefficient;
	    }
	}
    }
}

// symmetrize a vector by grouping several orbitals into a single one
//
// inputVector = reference on the vector to symmetrize
// nbrOrbitals = number of orbitals to group together
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Lz to the largest Lz
// lzSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
// return value = number of states that have been generated through the symmetrization procedure

int BosonOnSphereShort::SymmetrizeSingleStateOneIntoManyOrbital (LongRationalVector& inputVector, int nbrOrbitals, LongRationalVector*& symmetrizedVectors, int*& lzSectors)
{
  int TargetSpaceNbrOrbitals = (this->LzMax + 1) / nbrOrbitals;
  int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * this->NbrBosons;
  LongRationalVector* TmpVectors = new LongRationalVector[MaxTotalLz + 1];
  this->SymmetrizeSingleStateOneIntoManyOrbitalCore(inputVector, TmpVectors, nbrOrbitals, 0ul, this->LargeHilbertSpaceDimension);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i <= MaxTotalLz; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  ++NbrGeneratedSectors;
	}     
    }  
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new LongRationalVector[NbrGeneratedSectors];
  lzSectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i <= MaxTotalLz; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i];
	  lzSectors[NbrGeneratedSectors] = 2 * i - ((TargetSpaceNbrOrbitals - 1) * this->NbrBosons);
	  ++NbrGeneratedSectors;	  
	}     
    }  
  return NbrGeneratedSectors;
}
  
// symmetrize a vector by grouping several orbitals into a single one
//
// inputVector = reference on the vector to symmetrize
// symmetrizedVectors = array on the symmetrize states ranging from the smallest Lz to the largest Lz
// nbrOrbitals = number of orbitals to group together
// firstComponent = first component of the input vector that has to be symmetrized
// nbrComponents = number of components of the input vector that have to be symmetrized
// return value = symmetrized state

void BosonOnSphereShort::SymmetrizeSingleStateOneIntoManyOrbitalCore (LongRationalVector& inputVector, LongRationalVector* symmetrizedVectors, int nbrOrbitals, unsigned long firstComponent, unsigned long nbrComponents)
{
  long LastComponent = (long) (firstComponent + nbrComponents);
  int TargetSpaceNbrOrbitals = (this->LzMax + 1) / nbrOrbitals;
  int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * this->NbrBosons;
  BosonOnSphereShort** TargetSpaces = new BosonOnSphereShort* [MaxTotalLz + 1];
  unsigned long* TmpState = new unsigned long[TargetSpaceNbrOrbitals];
  for (int i = 0; i <= MaxTotalLz; ++i)
    {
      TargetSpaces[i] = 0;
    }
  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      LongRational TmpCoefficient = inputVector[i];
      this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	if (this->TemporaryState[k] > 1)
	  TmpCoefficient.FactorialDivide(this->TemporaryState[k]);
      for (int k = this->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
	this->TemporaryState[k] = 0x0ul;
      int TmpTotalLz = 0;
      for (int k = 0; k < TargetSpaceNbrOrbitals; ++k)
	{
	  unsigned long& TmpNbrParticles = TmpState[k];
	  TmpNbrParticles = 0x0ul;
	  for (int l = 0; l < nbrOrbitals; ++l)
	    TmpNbrParticles += this->TemporaryState[k * nbrOrbitals + l];
	  if (TmpNbrParticles > 1)
	    TmpCoefficient.FactorialMultiply(TmpNbrParticles);	    
	  TmpTotalLz += k * ((int) TmpNbrParticles);
	}
       
      if (TargetSpaces[TmpTotalLz] == 0)
	{
	  TargetSpaces[TmpTotalLz] = new BosonOnSphereShort (this->NbrBosons, 2 * TmpTotalLz - ((TargetSpaceNbrOrbitals - 1) * this->NbrBosons), 
							     TargetSpaceNbrOrbitals - 1);
	  symmetrizedVectors[TmpTotalLz] = LongRationalVector(TargetSpaces[TmpTotalLz]->HilbertSpaceDimension, true);
	}	  
      int TmpLzMax = TargetSpaces[TmpTotalLz]->LzMax;
      while (TmpState[TmpLzMax] == 0x0ul)
	--TmpLzMax;
      int TmpPos = TargetSpaces[TmpTotalLz]->FindStateIndex(TargetSpaces[TmpTotalLz]->BosonToFermion(TmpState, TmpLzMax), 
							    TmpLzMax + this->NbrBosons - 1);
      if (TmpPos < TargetSpaces[TmpTotalLz]->HilbertSpaceDimension)
	{
// 	  this->PrintState(cout, i) << " -> ";
// 	  TargetSpaces[TmpTotalLz]->PrintState(cout, TmpPos) << " (lz=" << TmpTotalLz << ") " <<  inputVector[i] << " " << TmpCoefficient << endl;
	  symmetrizedVectors[TmpTotalLz][TmpPos] += TmpCoefficient;
	}
    }
  delete[] TmpState;
  for (int i = 0; i <= MaxTotalLz; ++i)
    {
      if (TargetSpaces[i] != 0)
	{
	  delete TargetSpaces[i];
	}
    }
}
  
// symmetrize a vector by grouping several orbitals that are related by a periodicity condition on their momenta
//
// inputVector = reference on the vector to symmetrize
// periodicity = momentum periodicity (should be a multiple of the number of orbitals)
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Lz to the largest Lz
// lzSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
// return value = number of states that have been generated through the symmetrization procedure

int BosonOnSphereShort::SymmetrizeSingleStatePeriodicOrbitals (LongRationalVector& inputVector, int periodicity, LongRationalVector*& symmetrizedVectors, int*& lzSectors)
{
  int TargetSpaceNbrOrbitals = periodicity;
  int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * this->NbrBosons;
  LongRationalVector* TmpVectors = new LongRationalVector[MaxTotalLz + 1];
  this->SymmetrizeSingleStatePeriodicOrbitalCore(inputVector, TmpVectors, periodicity, 0ul, this->LargeHilbertSpaceDimension);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i <= MaxTotalLz; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  ++NbrGeneratedSectors;
	}     
    }  
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new LongRationalVector[NbrGeneratedSectors];
  lzSectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i <= MaxTotalLz; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i];
	  lzSectors[NbrGeneratedSectors] = 2 * i - ((TargetSpaceNbrOrbitals - 1) * this->NbrBosons);
	  ++NbrGeneratedSectors;	  
	}     
    }  
  return NbrGeneratedSectors;
}
  
// symmetrize a vector by grouping several orbitals that are related by a periodicity condition on their momenta
//
// inputVector = reference on the vector to symmetrize
// symmetrizedVectors = array on the symmetrize states ranging from the smallest Lz to the largest Lz
// periodicity = momentum periodicity (should be a multiple of the number of orbitals)
// firstComponent = first component of the input vector that has to be symmetrized
// nbrComponents = number of components of the input vector that have to be symmetrized
// return value = symmetrized state

void BosonOnSphereShort::SymmetrizeSingleStatePeriodicOrbitalCore (LongRationalVector& inputVector, LongRationalVector* symmetrizedVectors, int periodicity, unsigned long firstComponent, unsigned long nbrComponents)
{
  long LastComponent = (long) (firstComponent + nbrComponents);
  int TargetSpaceNbrOrbitals = periodicity;
  int NbrCopies = (this->LzMax + 1) / periodicity;
  int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * this->NbrBosons;
  BosonOnSphereShort** TargetSpaces = new BosonOnSphereShort* [MaxTotalLz + 1];
  unsigned long* TmpState = new unsigned long[TargetSpaceNbrOrbitals];
  for (int i = 0; i <= MaxTotalLz; ++i)
    {
      TargetSpaces[i] = 0;
    }
  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      LongRational TmpCoefficient = inputVector[i];
      this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	if (this->TemporaryState[k] > 1)
	  TmpCoefficient.FactorialDivide(this->TemporaryState[k]);
      for (int k = this->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
	this->TemporaryState[k] = 0x0ul;
      int TmpTotalLz = 0;
      for (int k = 0; k < TargetSpaceNbrOrbitals; ++k)
	{
	  unsigned long& TmpNbrParticles = TmpState[k];
	  TmpNbrParticles = 0x0ul;
	  for (int l = 0; l < NbrCopies; ++l)
	    TmpNbrParticles += this->TemporaryState[k + l * periodicity];
	  if (TmpNbrParticles > 1)
	    TmpCoefficient.FactorialMultiply(TmpNbrParticles);	    
	  TmpTotalLz += k * ((int) TmpNbrParticles);
	}
       
      if (TargetSpaces[TmpTotalLz] == 0)
	{
	  TargetSpaces[TmpTotalLz] = new BosonOnSphereShort (this->NbrBosons, 2 * TmpTotalLz - ((TargetSpaceNbrOrbitals - 1) * this->NbrBosons), 
							     TargetSpaceNbrOrbitals - 1);
	  symmetrizedVectors[TmpTotalLz] = LongRationalVector(TargetSpaces[TmpTotalLz]->HilbertSpaceDimension, true);
	}	  
      int TmpLzMax = TargetSpaces[TmpTotalLz]->LzMax;
      while (TmpState[TmpLzMax] == 0x0ul)
	--TmpLzMax;
      int TmpPos = TargetSpaces[TmpTotalLz]->FindStateIndex(TargetSpaces[TmpTotalLz]->BosonToFermion(TmpState, TmpLzMax), 
							    TmpLzMax + this->NbrBosons - 1);
      if (TmpPos < TargetSpaces[TmpTotalLz]->HilbertSpaceDimension)
	{
//  	  this->PrintState(cout, i) << " -> ";
//  	  TargetSpaces[TmpTotalLz]->PrintState(cout, TmpPos) << " (lz=" << TmpTotalLz << ") " <<  inputVector[i] << " " << TmpCoefficient << endl;
	  symmetrizedVectors[TmpTotalLz][TmpPos] += TmpCoefficient;
	}
    }
  delete[] TmpState;
  for (int i = 0; i <= MaxTotalLz; ++i)
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
// periodicity = momentum periodicity 
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest number of particles to the largest 
//                      number of partciles and the smallest Lz to the largest Lz
// nbrParticlesSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
// lzSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
// return value = number of states that have been generated through the symmetrization procedure

int BosonOnSphereShort::SymmetrizeSingleStatePeriodicSubsetOrbitals (LongRationalVector& inputVector, int firstOrbitalIndex, int periodicity, 
								     LongRationalVector*& symmetrizedVectors, int*& nbrParticlesSectors, int*& lzSectors)
{
  int TargetSpaceNbrOrbitals = (this->LzMax + 1) / periodicity;
  if ((((this->LzMax + 1) % periodicity) != 0) && firstOrbitalIndex < ((this->LzMax + 1) % periodicity))
    TargetSpaceNbrOrbitals += 1; 
  LongRationalVector** TmpVectors = new LongRationalVector*[this->NbrBosons + 1];
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * i;
      TmpVectors[i] = new LongRationalVector[MaxTotalLz + 1];
    }
  this->SymmetrizeSingleStatePeriodicSubsetOrbitalCore(inputVector, TmpVectors, firstOrbitalIndex, periodicity, 0ul, this->LargeHilbertSpaceDimension);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * i;
      for (int j = 0; j <= MaxTotalLz; ++j)
	{
	  if (TmpVectors[i][j].GetVectorDimension() != 0)	
	    {
	      ++NbrGeneratedSectors;
	    }     
	}
    }  
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new LongRationalVector[NbrGeneratedSectors];
  lzSectors = new int[NbrGeneratedSectors];
  nbrParticlesSectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * i;
      for (int j = 0; j <= MaxTotalLz; ++j)
	{
	  if (TmpVectors[i][j].GetVectorDimension() != 0)	
	    {
	      symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i][j];
	      lzSectors[NbrGeneratedSectors] = 2 * j - ((TargetSpaceNbrOrbitals - 1) * i);
	      nbrParticlesSectors[NbrGeneratedSectors] = i;
	      ++NbrGeneratedSectors;	  
	    }     
	}
    }  
  return NbrGeneratedSectors;
}

// symmetrize a vector by grouping several orbitals that are related by a periodicity condition on their momenta
//
// inputVector = reference on the vector to symmetrize
// firstOrbitalIndex = index of the first orbital to keep
// symmetrizedVectors = array on the symmetrize states ranging from the smallest Lz to the largest Lz
// periodicity = momentum periodicity (should be a multiple of the number of orbitals)
// firstComponent = first component of the input vector that has to be symmetrized
// nbrComponents = number of components of the input vector that have to be symmetrized
// return value = symmetrized state

void BosonOnSphereShort::SymmetrizeSingleStatePeriodicSubsetOrbitalCore (LongRationalVector& inputVector, LongRationalVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, 
									 unsigned long firstComponent, unsigned long nbrComponents)
{
  long LastComponent = (long) (firstComponent + nbrComponents);
  int TargetSpaceNbrOrbitals = (this->LzMax + 1) / periodicity;
  if ((((this->LzMax + 1) % periodicity) != 0) && firstOrbitalIndex < ((this->LzMax + 1) % periodicity))
    TargetSpaceNbrOrbitals += 1; 
  BosonOnSphereShort*** TargetSpaces = new BosonOnSphereShort** [this->NbrBosons + 1];
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * i;
      TargetSpaces[i] = new BosonOnSphereShort* [MaxTotalLz + 1];
      for (int j = 0; j <= MaxTotalLz; ++j)
	{
	  TargetSpaces[i][j] = 0;
	}
    }
  unsigned long* TmpState = new unsigned long[TargetSpaceNbrOrbitals];
  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      LongRational TmpCoefficient = inputVector[i];
      this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = this->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
	this->TemporaryState[k] = 0x0ul;
      int TmpTotalLz = 0;
      int TmpNbrParticles = 0;
      int Index = 0;
      for (int k = firstOrbitalIndex; k <= this->LzMax; k += periodicity)
	{
	  unsigned long& TmpNbrParticles2 = this->TemporaryState[k];
	  TmpState[Index] = TmpNbrParticles2;
	  TmpNbrParticles += ((int) TmpNbrParticles2);
	  TmpTotalLz += Index * ((int) TmpNbrParticles2);
	  ++Index;
	}
       
      if (TmpNbrParticles > 0)
	{
	  if (TargetSpaces[TmpNbrParticles][TmpTotalLz] == 0)
	    {
	      TargetSpaces[TmpNbrParticles][TmpTotalLz] = new BosonOnSphereShort (TmpNbrParticles, 2 * TmpTotalLz - ((TargetSpaceNbrOrbitals - 1) * TmpNbrParticles), 
										  TargetSpaceNbrOrbitals - 1);
	      symmetrizedVectors[TmpNbrParticles][TmpTotalLz] = LongRationalVector(TargetSpaces[TmpNbrParticles][TmpTotalLz]->HilbertSpaceDimension, true);
	    }	  
	  int TmpLzMax = TargetSpaces[TmpNbrParticles][TmpTotalLz]->LzMax;
	  while (TmpState[TmpLzMax] == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = TargetSpaces[TmpNbrParticles][TmpTotalLz]->FindStateIndex(TargetSpaces[TmpNbrParticles][TmpTotalLz]->BosonToFermion(TmpState, TmpLzMax), 
										 TmpLzMax + TmpNbrParticles - 1);
	  if (TmpPos < TargetSpaces[TmpNbrParticles][TmpTotalLz]->HilbertSpaceDimension)
	    {
	      symmetrizedVectors[TmpNbrParticles][TmpTotalLz][TmpPos] += TmpCoefficient;
	    }
	}
    }
  delete[] TmpState;
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * i;
      for (int j = 0; j <= MaxTotalLz; ++j)
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
  
// symmetrize a vector by grouping several orbitals into a single one
//
// inputVector = reference on the vector to symmetrize
// nbrOrbitals = number of orbitals to group together
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest Lz to the largest Lz
// lzSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
// return value = number of states that have been generated through the symmetrization procedure

int BosonOnSphereShort::SymmetrizeSingleStateOneIntoManyOrbital (RealVector& inputVector, int nbrOrbitals, RealVector*& symmetrizedVectors, bool unnormalizedBasisFlag, int*& lzSectors)
{
  int TargetSpaceNbrOrbitals = (this->LzMax + 1) / nbrOrbitals;
  int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * this->NbrBosons;
    
  RealVector* TmpVectors = new RealVector[MaxTotalLz + 1];
  this->SymmetrizeSingleStateOneIntoManyOrbitalCore(inputVector, TmpVectors, unnormalizedBasisFlag, nbrOrbitals, 0ul, this->LargeHilbertSpaceDimension);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i <= MaxTotalLz; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  ++NbrGeneratedSectors;
	}     
    } 
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new RealVector[NbrGeneratedSectors];
  lzSectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i <= MaxTotalLz; ++i)
    {
      if (TmpVectors[i].GetVectorDimension() != 0)	
	{
	  symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i];
	  lzSectors[NbrGeneratedSectors] = 2 * i - ((TargetSpaceNbrOrbitals - 1) * this->NbrBosons);
	  ++NbrGeneratedSectors;	  
	}     
    }  
  return NbrGeneratedSectors;
}

// symmetrize a vector by grouping several orbitals into a single one
//
// inputVector = reference on the vector to symmetrize
// symmetrizedVectors = array on the symmetrize states ranging from the smallest Lz to the largest Lz
// nbrOrbitals = number of orbitals to group together
// firstComponent = first component of the input vector that has to be symmetrized
// nbrComponents = number of components of the input vector that have to be symmetrized
// return value = symmetrized state

void BosonOnSphereShort::SymmetrizeSingleStateOneIntoManyOrbitalCore (RealVector& inputVector, RealVector* symmetrizedVectors, bool unnormalizedBasisFlag, int nbrOrbitals, unsigned long firstComponent, unsigned long nbrComponents)
{
  long LastComponent = (long) (firstComponent + nbrComponents);
  int TargetSpaceNbrOrbitals = (this->LzMax + 1) / nbrOrbitals;
  int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * this->NbrBosons;
  BosonOnSphereShort** TargetSpaces = new BosonOnSphereShort* [MaxTotalLz + 1];
  unsigned long* TmpState = new unsigned long[TargetSpaceNbrOrbitals];
  FactorialCoefficient Factorial1;
  for (int i = 0; i <= MaxTotalLz; ++i)
    {
      TargetSpaces[i] = 0;
    }
  if (unnormalizedBasisFlag == true)
  {
    for (long i = (long) firstComponent; i < LastComponent; ++i)
      {
	double TmpCoefficient = inputVector[i];
	Factorial1.SetToOne();
	this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
	for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	  if (this->TemporaryState[k] > 1)
	    Factorial1.FactorialDivide(this->TemporaryState[k]);
	for (int k = this->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
	  this->TemporaryState[k] = 0x0ul;
	int TmpTotalLz = 0;
	for (int k = 0; k < TargetSpaceNbrOrbitals; ++k)
	  {
	    unsigned long& TmpNbrParticles = TmpState[k];
	    TmpNbrParticles = 0x0ul;
	    for (int l = 0; l < nbrOrbitals; ++l)
	      TmpNbrParticles += this->TemporaryState[k * nbrOrbitals + l];
	    if (TmpNbrParticles > 1)
	      Factorial1.FactorialMultiply(TmpNbrParticles);	    
	    TmpTotalLz += k * ((int) TmpNbrParticles);
	  }
       
	if (TargetSpaces[TmpTotalLz] == 0)
	  {
	    TargetSpaces[TmpTotalLz] = new BosonOnSphereShort (this->NbrBosons, 2 * TmpTotalLz - ((TargetSpaceNbrOrbitals - 1) * this->NbrBosons), 
							     TargetSpaceNbrOrbitals - 1);
	    symmetrizedVectors[TmpTotalLz] = RealVector(TargetSpaces[TmpTotalLz]->HilbertSpaceDimension, true);
	  }	  
	int TmpLzMax = TargetSpaces[TmpTotalLz]->LzMax;
	while (TmpState[TmpLzMax] == 0x0ul)
	  --TmpLzMax;
	int TmpPos = TargetSpaces[TmpTotalLz]->FindStateIndex(TargetSpaces[TmpTotalLz]->BosonToFermion(TmpState, TmpLzMax), 
							    TmpLzMax + this->NbrBosons - 1);
	if (TmpPos < TargetSpaces[TmpTotalLz]->HilbertSpaceDimension)
	  {
// 	  this->PrintState(cout, i) << " -> ";
// 	  TargetSpaces[TmpTotalLz]->PrintState(cout, TmpPos) << " (lz=" << TmpTotalLz << ") " <<  inputVector[i] << " " << TmpCoefficient << endl;
	    symmetrizedVectors[TmpTotalLz][TmpPos] += Factorial1.GetNumericalValue()*TmpCoefficient;
	  }
      }
    }
    else
    {
     for (long i = (long) firstComponent; i < LastComponent; ++i)
      {
	double TmpCoefficient = inputVector[i];
	Factorial1.SetToOne();
	this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
	for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	  if (this->TemporaryState[k] > 1)
	    Factorial1.FactorialDivide(this->TemporaryState[k]);
// 	cout << Factorial1.GetNumericalValue() << endl;
	for (int k = this->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
	  this->TemporaryState[k] = 0x0ul;
	int TmpTotalLz = 0;
	for (int k = 0; k < TargetSpaceNbrOrbitals; ++k)
	  {
	    unsigned long& TmpNbrParticles = TmpState[k];
	    TmpNbrParticles = 0x0ul;
	    for (int l = 0; l < nbrOrbitals; ++l)
	    {
	      TmpNbrParticles += this->TemporaryState[k * nbrOrbitals + l];
	      for (int j = 0; j < this->TemporaryState[k * nbrOrbitals + l]; ++j)
	      {
		Factorial1.FactorialDivide(TargetSpaceNbrOrbitals - 1);
		Factorial1.FactorialMultiply(this->LzMax);
		Factorial1.FactorialDivide(nbrOrbitals*k + l);
		Factorial1.FactorialMultiply(k);
		Factorial1.FactorialDivide(this->LzMax - (nbrOrbitals*k + l));
		Factorial1.FactorialMultiply(TargetSpaceNbrOrbitals - 1 - k);
	      }
	    }
	    
	    if (TmpNbrParticles > 1)
	      Factorial1.FactorialMultiply(TmpNbrParticles);	    
	    TmpTotalLz += k * ((int) TmpNbrParticles);
	  }
       
	if (TargetSpaces[TmpTotalLz] == 0)
	  {
	    TargetSpaces[TmpTotalLz] = new BosonOnSphereShort (this->NbrBosons, 2 * TmpTotalLz - ((TargetSpaceNbrOrbitals - 1) * this->NbrBosons), 
							     TargetSpaceNbrOrbitals - 1);
	    symmetrizedVectors[TmpTotalLz] = RealVector(TargetSpaces[TmpTotalLz]->HilbertSpaceDimension, true);
	  }	  
	int TmpLzMax = TargetSpaces[TmpTotalLz]->LzMax;
	while (TmpState[TmpLzMax] == 0x0ul)
	  --TmpLzMax;
	int TmpPos = TargetSpaces[TmpTotalLz]->FindStateIndex(TargetSpaces[TmpTotalLz]->BosonToFermion(TmpState, TmpLzMax), 
							    TmpLzMax + this->NbrBosons - 1);
	
	if (TmpPos < TargetSpaces[TmpTotalLz]->HilbertSpaceDimension)
	  {
// 	  this->PrintState(cout, i) << " -> ";
// 	  TargetSpaces[TmpTotalLz]->PrintState(cout, TmpPos) << " (lz=" << TmpTotalLz << ") " <<  inputVector[i] << " " << TmpCoefficient << endl;
	    symmetrizedVectors[TmpTotalLz][TmpPos] += sqrt(Factorial1.GetNumericalValue())*TmpCoefficient;
	    
	  }
      } 
    }
  delete[] TmpState;
  for (int i = 0; i <= MaxTotalLz; ++i)
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
// periodicity = momentum periodicity 
// symmetrizedVectors = reference on the array on the symmetrized states ranging from the smallest number of particles to the largest 
//                      number of partciles and the smallest Lz to the largest Lz
// nbrParticlesSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
// lzSectors = reference on the array on twice the Lz sectors that have been generated through the symmetrization procedure
// return value = number of states that have been generated through the symmetrization procedure

int BosonOnSphereShort::SymmetrizeSingleStatePeriodicSubsetOrbitals (RealVector& inputVector, int firstOrbitalIndex, int periodicity, bool unnormalizedBasisFlag,
								     RealVector*& symmetrizedVectors, int*& nbrParticlesSectors, int*& lzSectors)
{
  int TargetSpaceNbrOrbitals = (this->LzMax + 1) / periodicity;
  if ((((this->LzMax + 1) % periodicity) != 0) && firstOrbitalIndex < ((this->LzMax + 1) % periodicity))
    TargetSpaceNbrOrbitals += 1; 
   RealVector** TmpVectors = new RealVector*[this->NbrBosons + 1];
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * i;
      TmpVectors[i] = new RealVector[MaxTotalLz + 1];
    }
  this->SymmetrizeSingleStatePeriodicSubsetOrbitalCore(inputVector, TmpVectors, firstOrbitalIndex, periodicity, unnormalizedBasisFlag, 0ul, this->LargeHilbertSpaceDimension);
  int NbrGeneratedSectors = 0;
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * i;
      for (int j = 0; j <= MaxTotalLz; ++j)
	{
	  if (TmpVectors[i][j].GetVectorDimension() != 0)	
	    {
	      ++NbrGeneratedSectors;
	    }     
	}
    }  
  if (NbrGeneratedSectors == 0)
    return 0;
  symmetrizedVectors = new RealVector[NbrGeneratedSectors];
  lzSectors = new int[NbrGeneratedSectors];
  nbrParticlesSectors = new int[NbrGeneratedSectors];
  NbrGeneratedSectors = 0;
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * i;
      for (int j = 0; j <= MaxTotalLz; ++j)
	{
	  if (TmpVectors[i][j].GetVectorDimension() != 0)	
	    {
	      symmetrizedVectors[NbrGeneratedSectors] = TmpVectors[i][j];
	      lzSectors[NbrGeneratedSectors] = 2 * j - ((TargetSpaceNbrOrbitals - 1) * i);
	      nbrParticlesSectors[NbrGeneratedSectors] = i;
	      ++NbrGeneratedSectors;	  
	    }     
	}
    }  
  return NbrGeneratedSectors;
}

// symmetrize a vector by grouping several orbitals that are related by a periodicity condition on their momenta
//
// inputVector = reference on the vector to symmetrize
// firstOrbitalIndex = index of the first orbital to keep
// symmetrizedVectors = array on the symmetrize states ranging from the smallest Lz to the largest Lz
// periodicity = momentum periodicity (should be a multiple of the number of orbitals)
// firstComponent = first component of the input vector that has to be symmetrized
// nbrComponents = number of components of the input vector that have to be symmetrized
// return value = symmetrized state

void BosonOnSphereShort::SymmetrizeSingleStatePeriodicSubsetOrbitalCore (RealVector& inputVector, RealVector** symmetrizedVectors, int firstOrbitalIndex, int periodicity, bool unnormalizedBasisFlag,
									 unsigned long firstComponent, unsigned long nbrComponents)
{
  long LastComponent = (long) (firstComponent + nbrComponents);
  int TargetSpaceNbrOrbitals = (this->LzMax + 1) / periodicity;
  if ((((this->LzMax + 1) % periodicity) != 0) && firstOrbitalIndex < ((this->LzMax + 1) % periodicity))
    TargetSpaceNbrOrbitals += 1; 
  BosonOnSphereShort*** TargetSpaces = new BosonOnSphereShort** [this->NbrBosons + 1];
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * i;
      TargetSpaces[i] = new BosonOnSphereShort* [MaxTotalLz + 1];
      for (int j = 0; j <= MaxTotalLz; ++j)
	{
	  TargetSpaces[i][j] = 0;
	}
    }
  unsigned long* TmpState = new unsigned long[TargetSpaceNbrOrbitals];
  if (unnormalizedBasisFlag == true)
    {
      for (long i = (long) firstComponent; i < LastComponent; ++i)
	{
	  double TmpCoefficient = inputVector[i];
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = this->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
	    this->TemporaryState[k] = 0x0ul;
	  int TmpTotalLz = 0;
	  int TmpNbrParticles = 0;
	  int Index = 0;
	  for (int k = firstOrbitalIndex; k <= this->LzMax; k += periodicity)
	    {
	      unsigned long& TmpNbrParticles2 = this->TemporaryState[k];
	      TmpState[Index] = TmpNbrParticles2;
	      TmpNbrParticles += ((int) TmpNbrParticles2);
	      TmpTotalLz += Index * ((int) TmpNbrParticles2);
	      ++Index;
	    }
	  
	  if (TmpNbrParticles > 0)
	    {
	      if (TargetSpaces[TmpNbrParticles][TmpTotalLz] == 0)
		{
		  TargetSpaces[TmpNbrParticles][TmpTotalLz] = new BosonOnSphereShort (TmpNbrParticles, 2 * TmpTotalLz - ((TargetSpaceNbrOrbitals - 1) * TmpNbrParticles), 
										      TargetSpaceNbrOrbitals - 1);
		  symmetrizedVectors[TmpNbrParticles][TmpTotalLz] = RealVector(TargetSpaces[TmpNbrParticles][TmpTotalLz]->HilbertSpaceDimension, true);
		}	  
	      int TmpLzMax = TargetSpaces[TmpNbrParticles][TmpTotalLz]->LzMax;
	      while (TmpState[TmpLzMax] == 0x0ul)
		--TmpLzMax;
	      int TmpPos = TargetSpaces[TmpNbrParticles][TmpTotalLz]->FindStateIndex(TargetSpaces[TmpNbrParticles][TmpTotalLz]->BosonToFermion(TmpState, TmpLzMax), 
										     TmpLzMax + TmpNbrParticles - 1);
	      if (TmpPos < TargetSpaces[TmpNbrParticles][TmpTotalLz]->HilbertSpaceDimension)
		{
		  symmetrizedVectors[TmpNbrParticles][TmpTotalLz][TmpPos] += TmpCoefficient;
		}
	    }
	}
    }
  else
    {
      for (long i = (long) firstComponent; i < LastComponent; ++i)
	{
	  FactorialCoefficient Factorial1;
	  Factorial1.SetToOne();
	  double TmpCoefficient = inputVector[i];
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  
	  for (int k = this->TemporaryStateLzMax + 1; k <= this->LzMax; ++k)
	    this->TemporaryState[k] = 0x0ul;
	  
	  int TmpTotalLz = 0;
	  int TmpNbrParticles = 0;
	  int Index = 0;
	  for (int k = firstOrbitalIndex; k <= this->LzMax; k += periodicity)
	    {
	      unsigned long& TmpNbrParticles2 = this->TemporaryState[k];
	      TmpState[Index] = TmpNbrParticles2;
	      TmpNbrParticles += ((int) TmpNbrParticles2);
	      TmpTotalLz += Index * ((int) TmpNbrParticles2);
	      ++Index;
	    }
	  
	  if (TmpNbrParticles > 0)
	    {
	      if (TargetSpaces[TmpNbrParticles][TmpTotalLz] == 0)
		{
		  TargetSpaces[TmpNbrParticles][TmpTotalLz] = new BosonOnSphereShort (TmpNbrParticles, 2 * TmpTotalLz - ((TargetSpaceNbrOrbitals - 1) * TmpNbrParticles), 
										      TargetSpaceNbrOrbitals - 1);
		  symmetrizedVectors[TmpNbrParticles][TmpTotalLz] = RealVector(TargetSpaces[TmpNbrParticles][TmpTotalLz]->HilbertSpaceDimension, true);
		}	  
	      int TmpLzMax = TargetSpaces[TmpNbrParticles][TmpTotalLz]->LzMax;
	      while (TmpState[TmpLzMax] == 0x0ul)
		--TmpLzMax;
	      int TmpPos = TargetSpaces[TmpNbrParticles][TmpTotalLz]->FindStateIndex(TargetSpaces[TmpNbrParticles][TmpTotalLz]->BosonToFermion(TmpState, TmpLzMax), 
										     TmpLzMax + TmpNbrParticles - 1);
	      for (int k = 0; k < TargetSpaceNbrOrbitals; ++k)
		{
		  for (int j = 0; j < this->TemporaryState[k*periodicity + firstOrbitalIndex]; ++j)
		    {
		      Factorial1.FactorialDivide(TargetSpaceNbrOrbitals - 1);
		      Factorial1.FactorialMultiply(k);
		      Factorial1.FactorialMultiply(TargetSpaceNbrOrbitals - 1 - k);
		    }
		  
		  for (int l = 0; l < periodicity; ++l)
		    {
		      if ((k * periodicity + l) <= this->LzMax)
			{
			  for (int j = 0; j < this->TemporaryState[k * periodicity + l]; ++j)
			    {		  
			      Factorial1.FactorialDivide(periodicity*k + l);
			      Factorial1.FactorialMultiply(this->LzMax);		  
			      Factorial1.FactorialDivide(this->LzMax - (periodicity*k + l));
			      if (l != firstOrbitalIndex)
				Factorial1.FactorialDivide(this->TemporaryState[k * periodicity + l]);
			    }
			}
		    }
		}
	      if (TmpPos < TargetSpaces[TmpNbrParticles][TmpTotalLz]->HilbertSpaceDimension)
		{
		  symmetrizedVectors[TmpNbrParticles][TmpTotalLz][TmpPos] += TmpCoefficient*sqrt(Factorial1.GetNumericalValue());
		}
	    }
	}
    }
  delete[] TmpState;
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      int MaxTotalLz = (TargetSpaceNbrOrbitals - 1) * i;
      for (int j = 0; j <= MaxTotalLz; ++j)
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

// Compute the product of two bsonic states, automatically dealing with reverse flux attachement
//
// bosonicState1 = reference on the first bosonic state
// bosonicState2 = reference on the second fermionic state
// outputVector = reference on the vector where the result will be stored
// bosonicSpace1 = pointer on the Hilbert Space associated to the first bosonic state
// bosonicSpace2 = pointer on the Hilbert Space associated to the second bosonic state
// minIndex = first component to compute (refering to the bosonic state)
// nbrComponents = number of components to compute (refering to the first bosonic state)
// unnormalizedFlag = true if the state should be written in the unnormalized basis
// architecture = pointer to the architecture

void BosonOnSphereShort::BosonicStateTimeBosonicState(RealVector& bosonicState1, RealVector& bosonicState2, RealVector& outputVector, 
						      BosonOnSphereShort* bosonicSpace1, BosonOnSphereShort* bosonicSpace2,
						      int minIndex, int nbrComponents, bool unnormalizedFlag, AbstractArchitecture* architecture)
{
  if (nbrComponents == this->GetHilbertSpaceDimension())
    outputVector.ClearVector();
  RealVector FinalState (this->GetHilbertSpaceDimension());
  FactorialCoefficient Coefficient;	
  int MaxIndex = minIndex + nbrComponents;
  double** ThreeOrbitalOverlaps = new double* [this->LzMax + 1];
  unsigned long* TmpSymmetricMonomial1 = new unsigned long[this->NbrBosons];
  unsigned long* TmpSymmetricMonomial2 = new unsigned long[this->NbrBosons];
  double TmpFactorial [this->NbrBosons + 1];
  TmpFactorial[0] = 1;
  TmpFactorial[1] = 1;
  for (int i = 2; i <= this->NbrBosons; ++i)
    TmpFactorial[i] = TmpFactorial[i - 1] / sqrt((double) i);
  
  if (this->LzMax >= bosonicSpace2->LzMax)
    {
      BinomialCoefficients Binomials(this->LzMax);
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  ThreeOrbitalOverlaps[i] = new double [bosonicSpace2->LzMax + 1];
	  double TmpFactor1 = log(((double) ((bosonicSpace2->LzMax + 1) * (bosonicSpace1->LzMax + 1))) / ((double) (this->LzMax + 1)) / (4.0 * M_PI)) - log(Binomials.GetNumericalCoefficient(this->LzMax, i));
	  for (int j = 0; (j <= bosonicSpace2->LzMax) && (j <= i); ++j)
	    {
	      if (unnormalizedFlag == false)
		{
		  ThreeOrbitalOverlaps[i][j] = 0.5 * (TmpFactor1 + log(Binomials.GetNumericalCoefficient(bosonicSpace2->LzMax, j)) 
						      + log(Binomials.GetNumericalCoefficient(bosonicSpace1->LzMax, i - j)));
		}
	      else
		{
		  ThreeOrbitalOverlaps[i][j] = 0.0;
		}
	    }
	}

      double TmpLogFactorial [this->NbrBosons + 1];
      TmpLogFactorial[0] = 0;
      TmpLogFactorial[1] = 0;
      for (int i = 2; i <= this->NbrBosons; ++i)
	TmpLogFactorial[i] = TmpLogFactorial[i - 1] + log((double) i);
      for (int j = minIndex; j < MaxIndex; ++j)
	{
	  if (bosonicState1[j] != 0.0)
	    {
	      int TmpLzMax = bosonicSpace1->FermionBasis->LzMax;
	      while ((bosonicSpace1->FermionBasis->StateDescription[j] >> TmpLzMax) == 0x0ul)
		{
		  --TmpLzMax;
		}
	      bosonicSpace1->ConvertToMonomial(bosonicSpace1->FermionBasis->StateDescription[j], TmpLzMax, TmpSymmetricMonomial1);
	      bosonicSpace1->FermionToBoson(bosonicSpace1->FermionBasis->StateDescription[j], TmpLzMax, bosonicSpace1->TemporaryState, 
					   bosonicSpace1->TemporaryStateLzMax);
	      double TmpFactor = 0.0;
	      for (int p = 0; p <= bosonicSpace1->TemporaryStateLzMax; ++p)
		{
		  TmpFactor -= TmpLogFactorial[bosonicSpace1->TemporaryState[p]];
		}
	      for (int i = 0; i < bosonicSpace2->HilbertSpaceDimension; ++i)
		{
		  if (bosonicState2[i] != 0.0)
		    {
		      int TmpLzMax2 = bosonicSpace2->FermionBasis->LzMax;
		      while ((bosonicSpace2->FermionBasis->StateDescription[i] >> TmpLzMax2) == 0x0ul)
			{
			  --TmpLzMax2;
			}
		      bosonicSpace2->ConvertToMonomial(bosonicSpace2->FermionBasis->StateDescription[i], TmpLzMax2, TmpSymmetricMonomial2);
		      if (unnormalizedFlag == false)
			{
			  bosonicSpace2->FermionToBoson(bosonicSpace2->FermionBasis->StateDescription[i], TmpLzMax2, bosonicSpace2->TemporaryState, 
							bosonicSpace2->TemporaryStateLzMax);
			  this->SymmetricMonomialTimesSymmetricMonomial(TmpSymmetricMonomial1, TmpSymmetricMonomial2, FinalState, ThreeOrbitalOverlaps);
			  for (int Index = 0; Index < FinalState.GetVectorDimension(); ++Index)
			    {
			      if (FinalState[Index] != 0.0)
				{
				  this->FermionToBoson(this->FermionBasis->StateDescription[Index], this->TemporaryState); 
				  Coefficient.SetToOne();
				  for(int p = 0; p <= this->LzMax; ++p)
				    {
				      if (this->TemporaryState[p] > 1)
					Coefficient.FactorialMultiply(this->TemporaryState[p]);
				    }		  
				  for(int p = 0; p <= bosonicSpace1->TemporaryStateLzMax; ++p)
				    {
				      if (bosonicSpace1->TemporaryState[p] > 1)
					Coefficient.FactorialDivide(bosonicSpace1->TemporaryState[p]);
				    }		  
				  for(int p = 0; p <= bosonicSpace2->TemporaryStateLzMax; ++p)
				    {
				      if (bosonicSpace2->TemporaryState[p] > 1)
					Coefficient.FactorialDivide(bosonicSpace2->TemporaryState[p]);
				    }		
				  outputVector[Index] += sqrt(Coefficient.GetNumericalValue()) * bosonicState1[j] * bosonicState2[i] * FinalState[Index];
				}
			    }
			}
		      else
			{
			  this->UnnormalizedSymmetricMonomialTimesSymmetricMonomial(TmpSymmetricMonomial1, TmpSymmetricMonomial2, FinalState);
			  for (int Index = 0; Index < FinalState.GetVectorDimension(); ++Index)
			    {
			      if (FinalState[Index] != 0.0)
				{
				  this->FermionToBoson(this->FermionBasis->StateDescription[Index], this->TemporaryState); 
				  double TmpFactor2 = TmpFactor;
				  for (int p = 0; p <= this->LzMax; ++p)
				    {
				      TmpFactor2 += TmpLogFactorial[this->TemporaryState[p]];
				    }
				  outputVector[Index] += exp(TmpFactor2) * bosonicState1[j] * bosonicState2[i] * FinalState[Index];
				}
			    }
			}
		    }
		}
	    }
	}
    }
  else
    {
      BinomialCoefficients Binomials(bosonicSpace2->LzMax + 1);
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  ThreeOrbitalOverlaps[i] = new double [bosonicSpace2->LzMax + 1];
	}
      int TmpMaxAngularMomentumLambdaLevel = bosonicSpace1->LzMax;
      double TmpPrefactor = sqrt (((double) (bosonicSpace1->LzMax + 1)) * ((double) (TmpMaxAngularMomentumLambdaLevel + 1)) / (4.0 * M_PI * (this->LzMax + 1)));
      ClebschGordanCoefficients Clebsch(TmpMaxAngularMomentumLambdaLevel, bosonicSpace2->LzMax);
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  for (int j = 0; j <= bosonicSpace2->LzMax; ++j)
	    {
	      if (((j - i) >= 0) && ((j - i) <= bosonicSpace1->LzMax))
		{
		  if (unnormalizedFlag == false)
		    {
		      ThreeOrbitalOverlaps[i][j] = (TmpPrefactor * Clebsch.GetCoefficient(- 2 * (j - i) + TmpMaxAngularMomentumLambdaLevel, (2 * j) - bosonicSpace2->LzMax, this->LzMax) 
						    * Clebsch.GetCoefficient(-bosonicSpace1->LzMax, bosonicSpace2->LzMax, this->LzMax));
		    }
		  else
		    {
		      ThreeOrbitalOverlaps[j - i][j] = 1.0;
		    }
		}
	    }
	}
      for (int j = minIndex; j < MaxIndex; ++j)
	{
	  if (bosonicState1[j] != 0.0)
	    {
	      int TmpLzMax = bosonicSpace1->FermionBasis->LzMax;
	      while ((bosonicSpace1->FermionBasis->StateDescription[j] >> TmpLzMax) == 0x0ul)
		{
		  --TmpLzMax;
		}
	      bosonicSpace1->ConvertToMonomial(bosonicSpace1->FermionBasis->StateDescription[j], TmpLzMax, TmpSymmetricMonomial1);
	      bosonicSpace1->FermionToBoson(bosonicSpace1->FermionBasis->StateDescription[j], TmpLzMax, bosonicSpace1->TemporaryState, 
					    bosonicSpace1->TemporaryStateLzMax);
	      double TmpFactor = 1.0;
	      for (int p = 0; p <= bosonicSpace1->TemporaryStateLzMax; ++p)
		{
		  TmpFactor *= TmpFactorial[bosonicSpace1->TemporaryState[p]];
		}
	      for (int i = 0; i < bosonicSpace2->HilbertSpaceDimension; ++i)
		{
		  if (bosonicState2[i] != 0.0)
		    {
		      int TmpLzMax2 = bosonicSpace2->FermionBasis->LzMax;
		      while ((bosonicSpace2->FermionBasis->StateDescription[i] >> TmpLzMax2) == 0x0ul)
			{
			  --TmpLzMax2;
			}
		      bosonicSpace2->ConvertToMonomial(bosonicSpace2->FermionBasis->StateDescription[i], TmpLzMax2, TmpSymmetricMonomial2);
		      bosonicSpace2->FermionToBoson(bosonicSpace2->FermionBasis->StateDescription[i], TmpLzMax2, bosonicSpace2->TemporaryState, 
						    bosonicSpace2->TemporaryStateLzMax);
  		      this->ReverseSymmetricMonomialTimesSymmetricMonomial(TmpSymmetricMonomial1, TmpSymmetricMonomial2, FinalState, ThreeOrbitalOverlaps);
//		      double TmpFactor2 = TmpFactor * bosonicState1[j] * bosonicState2[i];
		      double TmpFactor2 = TmpFactor * bosonicState1[j] * bosonicState2[i];
		      for (int Index = 0; Index < FinalState.GetVectorDimension(); ++Index)
			{
			  if (FinalState[Index] != 0.0)
			    {
			      this->FermionToBoson(this->FermionBasis->StateDescription[Index], this->TemporaryState); 
			      Coefficient.SetToOne();
			      for(int p = 0; p <= this->LzMax; ++p)
				{
				  if (this->TemporaryState[p] > 1)
				    Coefficient.FactorialMultiply(this->TemporaryState[p]);
				}		  
			      for(int p = 0; p <= bosonicSpace1->LzMax; ++p)
				{
				  if (bosonicSpace1->TemporaryState[p] > 1)
				    Coefficient.FactorialDivide(bosonicSpace1->TemporaryState[p]);
				}		  
			      for(int p = 0; p <= bosonicSpace2->LzMax; ++p)
				{
				  if (bosonicSpace2->TemporaryState[p] > 1)
				    Coefficient.FactorialDivide(bosonicSpace2->TemporaryState[p]);
				}		  
			      if (unnormalizedFlag == false)
				{
				  outputVector[Index] += sqrt(Coefficient.GetNumericalValue()) * TmpFactor2 * FinalState[Index];
				}
			      else
				{
				  outputVector[Index] +=  Coefficient.GetNumericalValue() * bosonicState1[j] * bosonicState2[i] * FinalState[Index];
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

// Compute the product of two symmetric monomials
//
// symmetricMonomial1 = first symmetric monomial
// symmetricMonomial2 = second symmetric monomial
// finalState = reference on the vector the produced state will be stored
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored

void BosonOnSphereShort::SymmetricMonomialTimesSymmetricMonomial (unsigned long* symmetricMonomial1, unsigned long* symmetricMonomial2, 
								  RealVector& finalState, double** threeOrbitalOverlaps)
{
  unsigned long TmpNbrStates = 0;
  unsigned long TmpState [this->NbrBosons];
  unsigned long TmpFinalState;
  int TmpLzMax;
  int TmpHeapArray [this->NbrBosons];
  int TmpDim = this->NbrBosons;
  for (int i = 0; i < TmpDim; ++i)
    {
      TmpHeapArray[i] = 0;
    }
  finalState.ClearVector();

  double TmpFactor = 0.0;
  int Tmp = 0;
  for (int i = 0; i < this->NbrBosons; ++i)
    {
      TmpState[i] = symmetricMonomial2[i] + symmetricMonomial1[i];
      TmpFactor += threeOrbitalOverlaps[TmpState[i]][symmetricMonomial2[i]];
    }
  TmpFinalState = this->ConvertFromUnsortedMonomial(TmpState, TmpLzMax);
  int TmpPos = this->FindStateIndex(TmpFinalState, TmpLzMax);
  if (TmpPos != this->HilbertSpaceDimension)
    {
      finalState[TmpPos] += exp(TmpFactor);
    }

  while (std::prev_permutation(symmetricMonomial2, symmetricMonomial2 + this->NbrBosons))
    {
      TmpFactor = 0.0;
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  TmpState[i] = symmetricMonomial2[i] + symmetricMonomial1[i];
	  TmpFactor += threeOrbitalOverlaps[TmpState[i]][symmetricMonomial2[i]];
	}
      TmpFinalState = this->ConvertFromUnsortedMonomial(TmpState, TmpLzMax);
      TmpPos = this->FindStateIndex(TmpFinalState, TmpLzMax);
      if (TmpPos != this->HilbertSpaceDimension)
	{
	  finalState[TmpPos] += exp(TmpFactor);
	}
    }
}

// Compute the product of two symmetric monomials, assuming an unnormalized basis
//
// symmetricMonomial1 = first symmetric monomial
// symmetricMonomial2 = second symmetric monomial
// finalState = reference on the vector the produced state will be stored

void BosonOnSphereShort::UnnormalizedSymmetricMonomialTimesSymmetricMonomial (unsigned long* symmetricMonomial1, unsigned long* symmetricMonomial2, 
									      RealVector& finalState)
{
  int TmpState [this->LzMax + 1];
  int IntSymmetricMonomial1 [this->NbrBosons];
  int IntSymmetricMonomial2 [this->NbrBosons];
  unsigned long TmpFinalState;
  int TmpLzMax;
  finalState.ClearVector();

  for (int i = 0; i <= this->LzMax; ++i)
    TmpState[i] = 0;
  for (int i = 0; i < this->NbrBosons; ++i)
    {
      IntSymmetricMonomial1[i] = (int) symmetricMonomial1[i];
      IntSymmetricMonomial2[i] = (int) symmetricMonomial2[i];
      TmpState[IntSymmetricMonomial1[i] + IntSymmetricMonomial2[i]]++;
    }
  TmpLzMax = this->LzMax;
  while (TmpState[TmpLzMax] == 0)
    {
      --TmpLzMax;
    }
  int Shift = 0;
  TmpFinalState = 0x0ul;
  for (int i = 0; i <= TmpLzMax; ++i)
     {
      TmpFinalState |= ((1ul << TmpState[i]) - 1ul) << Shift;
      Shift += TmpState[i];
      ++Shift;
    }
  TmpLzMax += this->NbrBosons;
  --TmpLzMax;
  int TmpPos = this->FindStateIndex(TmpFinalState, TmpLzMax);
  if (TmpPos != this->HilbertSpaceDimension)
    {
      finalState[TmpPos]++;
    }

  while (std::prev_permutation(IntSymmetricMonomial2, IntSymmetricMonomial2 + this->NbrBosons))
    {
      for (int i = 0; i <= this->LzMax; ++i)
	TmpState[i] = 0;
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  TmpState[IntSymmetricMonomial2[i] + IntSymmetricMonomial1[i]]++;
	}
      TmpLzMax = this->LzMax;
      while (TmpState[TmpLzMax] == 0)
	{
	  --TmpLzMax;
	}
      Shift = 0;
      TmpFinalState = 0x0ul;
      for (int i = 0; i <= TmpLzMax; ++i)
	{
	  TmpFinalState |= ((1ul << TmpState[i]) - 1ul) << Shift;
	  Shift += TmpState[i];
	  ++Shift;
	}
      TmpLzMax += this->NbrBosons;
      --TmpLzMax;
      TmpPos = this->FindStateIndex(TmpFinalState, TmpLzMax);
      if (TmpPos != this->HilbertSpaceDimension)
	{
	  finalState[TmpPos]++;
	}
    }
}
// {
//   unsigned long TmpState [this->NbrBosons];
//   unsigned long TmpFinalState;
//   int TmpLzMax;
//   finalState.ClearVector();

//   for (int i = 0; i < this->NbrBosons; ++i)
//     {
//       TmpState[i] = symmetricMonomial2[i] + symmetricMonomial1[i];
//     }
//   TmpFinalState = this->ConvertFromUnsortedMonomial(TmpState, TmpLzMax);
//   int TmpPos = this->FindStateIndex(TmpFinalState, TmpLzMax);
//   if (TmpPos != this->HilbertSpaceDimension)
//     {
//       finalState[TmpPos]++;
//     }

//   while (std::prev_permutation(symmetricMonomial2, symmetricMonomial2 + this->NbrBosons))
//     {
//       for (int i = 0; i < this->NbrBosons; ++i)
// 	{
// 	  TmpState[i] = symmetricMonomial2[i] + symmetricMonomial1[i];
// 	}
//       TmpFinalState = this->ConvertFromUnsortedMonomial(TmpState, TmpLzMax);
//       TmpPos = this->FindStateIndex(TmpFinalState, TmpLzMax);
//       if (TmpPos != this->HilbertSpaceDimension)
// 	{
// 	  finalState[TmpPos]++;
// 	}
//     }
// }

// Compute the product of two symmetric monomials, assuming an unnormalized basis
//
// symmetricMonomial1 = first symmetric monomial
// symmetricMonomial2 = second symmetric monomial
// finalStateConfigurations = fermionic configurations that are generated durig the multiplication
// finalStateWeights = weights associated to each generated fermionic configurations
// return value = number of generated fermionic configurations

int BosonOnSphereShort::UnnormalizedSymmetricMonomialTimesSymmetricMonomial (unsigned long* symmetricMonomial1, unsigned long* symmetricMonomial2, 
									     unsigned long* finalStateConfigurations, double* finalStateWeights)
{
  int TmpArraySize = 0;
  unsigned long TmpState [this->NbrBosons];
  unsigned long TmpFinalState;
  int TmpLzMax;
  int Tmp = 0;

  for (int i = 0; i < this->NbrBosons; ++i)
    {
      TmpState[i] = symmetricMonomial2[i] + symmetricMonomial1[i];
    }
  TmpFinalState = this->ConvertFromUnsortedMonomial(TmpState);
  TmpArraySize += SearchInArrayAndSetWeight(TmpFinalState, finalStateConfigurations, finalStateWeights, TmpArraySize, 1.0);

  while (std::prev_permutation(symmetricMonomial2, symmetricMonomial2 + this->NbrBosons))
    {
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  TmpState[i] = symmetricMonomial2[i] + symmetricMonomial1[i];
	}
      TmpFinalState = this->ConvertFromUnsortedMonomial(TmpState);
      TmpArraySize += SearchInArrayAndSetWeight(TmpFinalState, finalStateConfigurations, finalStateWeights, TmpArraySize, 1.0);
    }
  return TmpArraySize;
}

// Compute the product of two symmetric monomials, assuming a reverse flux attachment for the first symmetric monomial
//
// symmetricMonomial1 = first symmetric monomial
// symmetricMonomial2 = second symmetric monomial
// finalState = reference on the vector the produced state will be stored
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored

void BosonOnSphereShort::ReverseSymmetricMonomialTimesSymmetricMonomial (unsigned long* symmetricMonomial1, unsigned long* symmetricMonomial2,
									 RealVector& finalState, double** threeOrbitalOverlaps)
{
  cout << "BosonOnSphereShort::ReverseSymmetricMonomialTimesSymmetricMonomial is not implemented" << endl;
}

