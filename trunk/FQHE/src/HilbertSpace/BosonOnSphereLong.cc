////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system size such that            //
//         LzMax + NbrBosons - 1 < 127 or 63 (64 bits or 32bits systems)      //
//                                                                            //
//                        last modification : 14/03/2010                      //
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
#include "HilbertSpace/BosonOnSphereLong.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h" 
#include "GeneralTools/StringTools.h"
#include "GeneralTools/ArrayTools.h"

#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


// default constructor
//

BosonOnSphereLong::BosonOnSphereLong ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson

BosonOnSphereLong::BosonOnSphereLong (int nbrBosons, int totalLz, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->FermionBasis = new FermionOnSphereLong(nbrBosons, totalLz, lzMax + nbrBosons - 1);
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

BosonOnSphereLong::BosonOnSphereLong(const BosonOnSphereLong& bosons)
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
  this->FermionBasis = (FermionOnSphereLong*) bosons.FermionBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
}

// destructor
//

BosonOnSphereLong::~BosonOnSphereLong ()
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
      delete this->KeptCoordinates;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereLong& BosonOnSphereLong::operator = (const BosonOnSphereLong& bosons)
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
  this->FermionBasis = (FermionOnSphereLong*) bosons.FermionBasis->Clone();
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereLong::Clone()
{
  return new BosonOnSphereLong(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnSphereLong::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->TargetSpace = (BosonOnSphereLong*) targetSpace;
  this->FermionBasis->SetTargetSpace(this->TargetSpace->FermionBasis);
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int BosonOnSphereLong::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereLong::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereLong::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereLong::ExtractSubspace (AbstractQuantumNumber& q, 
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

int BosonOnSphereLong::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

int BosonOnSphereLong::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
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

double BosonOnSphereLong::AA (int index, int n1, int n2)
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

double BosonOnSphereLong::ProdA (int index, int* n, int nbrIndices)
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

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereLong::AdAd (int m1, int m2, double& coefficient)
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

int BosonOnSphereLong::ProdAd (int* m, int nbrIndices, double& coefficient)
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

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereLong::AdA (int index, int m)
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

double BosonOnSphereLong::AdA (long index, int m)
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

int BosonOnSphereLong::AdA (int index, int m, int n, double& coefficient)
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
      while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
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
  return this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereLong::PrintState (ostream& Str, int state)
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

ostream& BosonOnSphereLong::PrintStateMonomial (ostream& Str, long state)
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


// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int BosonOnSphereLong::FindStateIndex(char* stateDescription)
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

// evaluate wave function in real space using a given basis and only for a given range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex BosonOnSphereLong::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
						  int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
  ComplexMatrix Perm(this->NbrBosons, this->NbrBosons);
  ComplexMatrix Functions(this->LzMax + 1, this->NbrBosons);
  RealVector TmpCoordinates(2);
  int* Indices = new int [this->NbrBosons];
  int Pos;
  int Lz;
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      TmpCoordinates[0] = position[j << 1];
      TmpCoordinates[1] = position[1 + (j << 1)];
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

RealSymmetricMatrix  BosonOnSphereLong::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState)
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
 	  BosonOnSphereLong TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  ULONGLONG  TmpState2 = ((ULONGLONG) 0x0);
	  for (int i = 0; i < nbrBosonSector; ++i)
	    TmpState2 |= 0x1ul << i;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      ULONGLONG TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | TmpState2;
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
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
 	  BosonOnSphereLong TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      ULONGLONG TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
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

  int TmpComplementarySubsystemLzMax = this->LzMax - subsytemSize;
  int MinIndex = 0;
  int MaxIndex = this->HilbertSpaceDimension - 1;
  if (nbrBosonSector == 1)
    {
      double TmpValue = 0.0;
      BosonOnSphereLong TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  ULONGLONG TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | (0x1ul << ShiftedLzSector);
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
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
      BosonOnSphereLong TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
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


  int TmpNbrBosons;
  int TmpTotalLz;
  int TmpIndex;
  BosonOnSphereLong TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  BosonOnSphereLong TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      ULONGLONG TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  ULONGLONG TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
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

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& BosonOnSphereLong::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0 / state[reference];
  state[reference] = 1.0;
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
    if (i != reference)
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

RealVector& BosonOnSphereLong::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = state[reference];
  state[reference] = 1.0;
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
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    if (i != reference)
      {
	this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
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
	    Factorial = ReferenceFactorial;
	    this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
				 this->TemporaryState, this->TemporaryStateLzMax);
	    for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	      if (this->TemporaryState[k] > 1)
		Factorial.FactorialDivide(this->TemporaryState[k]);
	    Coefficient *= sqrt(Factorial.GetNumericalValue());
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
      }
  state /= state.Norm();
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

RealVector& BosonOnSphereLong::FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
					   ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
					   bool symmetrizedFlag, double coefficient)
{
  BosonOnSphereLong* LeftSpace = (BosonOnSphereLong*) leftSpace;
  BosonOnSphereLong* RightSpace = (BosonOnSphereLong*) rightSpace;
   int StateShift = RightSpace->FermionBasis->LzMax + padding + 2;
  for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
    {
      ULONGLONG TmpState1 = LeftSpace->FermionBasis->StateDescription[i] << StateShift;
      double Coefficient = coefficient * leftVector[i];
      int TmpLzMax = this->FermionBasis->LzMax;
      while ((TmpState1 >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      if (symmetrizedFlag == false)
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      ULONGLONG TmpState2 = RightSpace->FermionBasis->StateDescription[j];
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
	      ULONGLONG TmpState2 = RightSpace->FermionBasis->StateDescription[j];
	      TmpState2 |= TmpState1;
	      double Coefficient2 = Coefficient;
	      Coefficient2 *= rightVector[j];	  
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
	      outputVector[TmpIndex] = Coefficient2;
	      ULONGLONG TmpState3 = this->FermionBasis->GetSymmetricState(TmpState2);
	      if (TmpState3 != TmpState2)
		{
		  int TmpLzMax2 = this->FermionBasis->LzMax;
		  while ((TmpState3 >> TmpLzMax2) == ((ULONGLONG) 0x0ul))
		    --TmpLzMax2;
		  TmpIndex = this->FermionBasis->FindStateIndex(TmpState3, TmpLzMax2);
		  outputVector[TmpIndex] = Coefficient2;      
		}
	    }
	}
    }
  return outputVector;
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

RealVector& BosonOnSphereLong::ProductRules (RealVector& outputVector, RealVector& inputVector, ParticleOnSphere* inputSpace, 
					     int* commonPattern, int commonPatterSize, int* addedPattern, int addedPatterSize,
					     double coefficient, bool symmetrizedFlag)
{
  BosonOnSphereLong* InputSpace = (BosonOnSphereLong*) inputSpace;
  int NbrParticlesCommonPattern = 0;
  ULONGLONG InputPattern = ((ULONGLONG) 0x0ul);
  int TmpIndex = 0;
  for (int i = commonPatterSize - 1; i >= 0; --i)
    {
      for (int j = 0; j < commonPattern[i]; ++j)
	{
	  InputPattern |= ((ULONGLONG) 0x1u) << TmpIndex;
	  ++TmpIndex;
	}
      ++TmpIndex;
      NbrParticlesCommonPattern += commonPattern[i];
    }
  int NbrParticlesAddedPattern = 0;
  ULONGLONG OutputPattern = ((ULONGLONG) 0x0ul);
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
  ULONGLONG InputMask = ((((ULONGLONG) 0x1ul) << (commonPatterSize + NbrParticlesCommonPattern)) - ((ULONGLONG) 1ul)) << (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 1);
  ULONGLONG InputMask2 = ~InputMask;  
  OutputPattern |= InputPattern << (addedPatterSize + NbrParticlesAddedPattern);
  OutputPattern <<= (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 2);
  InputPattern <<= (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 2); 
//  cout << hex << InputMask << " " << OutputPattern << " " << InputPattern << dec << endl;
  int OutputLzMax = this->FermionBasis->LzMax;
  while (((OutputPattern >> OutputLzMax) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
    --OutputLzMax;
  long Count = 0l;
  if (symmetrizedFlag == false)
    {
      for (int i = 0; i <  InputSpace->HilbertSpaceDimension; ++i)
	{
	  ULONGLONG TmpState1 = InputSpace->FermionBasis->StateDescription[i];
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
	  ULONGLONG TmpState1 = InputSpace->FermionBasis->StateDescription[i];
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
		  while ((TmpState2 >> TmpLzMax2) == ((ULONGLONG) 0x0ul))
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

double BosonOnSphereLong::JackSqrNormalization (RealVector& outputVector, long minIndex, long nbrComponents)
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

// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the  component
int BosonOnSphereLong::GetLzValue(int j)
{
  return this->TotalLz;
}

