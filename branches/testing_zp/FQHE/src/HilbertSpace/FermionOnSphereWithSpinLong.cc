////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//       class of fermions on sphere with spin that allow LzMax up to         //
//                  63 (for systems with 128 bit integer support)             //
//               or 31 (on 32 bit systems without 128 bit integer support)    //
//                                                                            //
//                        last modification : 26/09/2008                      //
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
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include <math.h>
#include <bitset>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;

#define WANT_LAPACK

#ifdef __LAPACK__
#ifdef WANT_LAPACK
#define  __USE_LAPACK_HERE__
#endif
#endif


// default constructor
//

FermionOnSphereWithSpinLong::FermionOnSphereWithSpinLong()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twce the total spin value
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinLong::FermionOnSphereWithSpinLong (int nbrFermions, int totalLz, int lzMax, int totalSpin, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->HilbertSpaceDimension = (int) this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
										 (this->TotalSpin + this->NbrFermions) >> 1);
  this->Flag.Initialize();
  this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  

  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						     (this->TotalSpin + this->NbrFermions) >> 1, 0l);
  this->GenerateLookUpTable(memory);
  
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
  UsedMemory = this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
  cout << "memory requested for lookup table = ";
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

FermionOnSphereWithSpinLong::FermionOnSphereWithSpinLong(const FermionOnSphereWithSpinLong& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
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
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
}

// destructor
//

FermionOnSphereWithSpinLong::~FermionOnSphereWithSpinLong ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < (2 * this->NbrLzValue); ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinLong& FermionOnSphereWithSpinLong::operator = (const FermionOnSphereWithSpinLong& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
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
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinLong::Clone()
{
  return new FermionOnSphereWithSpinLong(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnSphereWithSpinLong::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnSphereWithSpinLong::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphereWithSpinLong::ExtractSubspace (AbstractQuantumNumber& q, 
							SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_u_m1 a^+_u_m2 a_u_n1 a_u_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLong::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_d_m1 a^+_d_m2 a_d_n1 a_d_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLong::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

int FermionOnSphereWithSpinLong::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}


// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereWithSpinLong::AduAu (int index, int m)
{
  if ((this->StateDescription[index] & (0x2l << (m << 1))) != 0)
    return 1.0;
  else
    return 0.0;
}

// apply a^+_d_m a_d_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_d_m a_d_m

double FermionOnSphereWithSpinLong::AddAd (int index, int m)
{
  if ((this->StateDescription[index] & (((ULONGLONG) 0x1ul) << (m << 1))) != 0)
    return 1.0;
  else
    return 0.0;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLong::AuAu (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  ++n1;
  n2 <<= 1;
  ++n2;
  if (((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n1)) == 0) || ((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLong::AdAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  n2 <<= 1;
  if (((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n1)) == 0) || ((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLong::AuAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  ++n1;
  n2 <<= 1;
  if (((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n1)) == 0) || ((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n2)) == 0))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLong::AduAdu (int m1, int m2, double& coefficient)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  ++m2;
  if (((TmpState & (((ULONGLONG) 0x1ul) << m1)) != 0) || ((TmpState & (((ULONGLONG) 0x1ul) << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80))  & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80))  & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLong::AddAdd (int m1, int m2, double& coefficient)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  m2 <<= 1;
  if (((TmpState & (((ULONGLONG) 0x1ul) << m1)) != 0) || ((TmpState & (((ULONGLONG) 0x1ul) << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80))  & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80))  & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLong::AduAdd (int m1, int m2, double& coefficient)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  if (((TmpState & (((ULONGLONG) 0x1ul) << m1)) != 0) || ((TmpState & (((ULONGLONG) 0x1ul) << m2)) != 0))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80))  & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80))  & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLong::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  this->ProdALzMax = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (n[i] << 1) + spinIndices[i];
      if ((this->ProdATemporaryState & (((ULONGLONG) ((ULONGLONG) 0x1ul)) << Index)) == 0)
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#ifdef __128_BIT_LONGLONG__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 64)) & this->SignLookUpTableMask[Index + 64]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 80))  & this->SignLookUpTableMask[Index + 80]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 96)) & this->SignLookUpTableMask[Index + 96]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 112)) & this->SignLookUpTableMask[Index + 112]];
#endif
      this->ProdATemporaryState &= ~(((ULONGLONG) ((ULONGLONG) 0x1ul)) << Index);
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;

  return Coefficient;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each annihilation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLong::ProdA (int index, int* n, int spinIndices, int nbrIndices)
{
  this->ProdALzMax = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (n[i] << 1) + ((spinIndices >> i) & 0x1);
      if ((this->ProdATemporaryState & (((ULONGLONG) ((ULONGLONG) 0x1ul)) << Index)) == 0)
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#ifdef __128_BIT_LONGLONG__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 64)) & this->SignLookUpTableMask[Index + 64]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 80))  & this->SignLookUpTableMask[Index + 80]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 96)) & this->SignLookUpTableMask[Index + 96]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 112)) & this->SignLookUpTableMask[Index + 112]];
#endif
      this->ProdATemporaryState &= ~(((ULONGLONG) ((ULONGLONG) 0x1ul)) << Index);
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;

  return Coefficient;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLong::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  ULONGLONG TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (m[i] << 1) + spinIndices[i];
      if ((TmpState & (((ULONGLONG) ((ULONGLONG) 0x1ul)) << Index)) != 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#ifdef __128_BIT_LONGLONG__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 64)) & this->SignLookUpTableMask[Index + 64]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 80))  & this->SignLookUpTableMask[Index + 80]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 96)) & this->SignLookUpTableMask[Index + 96]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 112)) & this->SignLookUpTableMask[Index + 112]];
#endif
	}
      TmpState |= (((ULONGLONG) ((ULONGLONG) 0x1ul)) << Index);
    }
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each creation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLong::ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  ULONGLONG TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (m[i] << 1) + ((spinIndices >> i) & 0x1);
      if ((TmpState & (((ULONGLONG) ((ULONGLONG) 0x1ul)) << Index)) != 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#ifdef __128_BIT_LONGLONG__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 64)) & this->SignLookUpTableMask[Index + 64]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 80))  & this->SignLookUpTableMask[Index + 80]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 96)) & this->SignLookUpTableMask[Index + 96]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 112)) & this->SignLookUpTableMask[Index + 112]];
#endif
	}
      TmpState |= (((ULONGLONG) ((ULONGLONG) 0x1ul)) << Index);
    }
  return this->FindStateIndex(TmpState, NewLzMax);
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereWithSpinLong::FindStateIndex(ULONGLONG stateDescription, int lzmax)
{
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  ULONGLONG CurrentState = this->StateDescription[PosMid];
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

ostream& FermionOnSphereWithSpinLong::PrintState (ostream& Str, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  ULONGLONG Tmp;
  for (int i = this->NbrLzValue-1; i >=0 ; --i)
    {
      Tmp = ((TmpState >> (i << 1)) & ((ULONGLONG) 0x3ul));
      if (Tmp == ((ULONGLONG) 0x1ul))
	Str << "d ";
      else if (Tmp == ((ULONGLONG) 0x2l))
	Str << "u ";
      else if (Tmp == ((ULONGLONG) 0x3l))
	Str << "X ";
      else Str << "0 ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinLong::GenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSpin, long pos)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalSpin > nbrFermions) )
    return pos;
  if ((nbrFermions == 0) && (totalLz == 0) && (totalSpin == 0))
      {
	this->StateDescription[pos] = 0x0ul;
	return (pos + 1l);
      }
    
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz))
    return pos;
    
  if (nbrFermions == 1) 
    if (lzMax >= totalLz)
      {
	this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << ((totalLz << 1) + totalSpin);
	return (pos + 1l);
      }
    else
      return pos;

  if ((lzMax == 0)  && (totalLz != 0))
    return pos;


  long TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 1,  pos);
  ULONGLONG Mask = ((ULONGLONG) 0x3ul) << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1,  pos);
  Mask = ((ULONGLONG) 0x2ul) << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin,  pos);
  Mask = ((ULONGLONG) 0x1ul) << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  return this->GenerateStates(nbrFermions, lzMax - 1, totalLz, totalSpin, pos);
};


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereWithSpinLong::GenerateLookUpTable(unsigned long memory)
{
  // get every highest bit poisition
  ULONGLONG TmpPosition = this->StateDescription[0];
#ifdef __128_BIT_LONGLONG__
  int CurrentHighestBit = 127;
#else
  int CurrentHighestBit = 63;
#endif
  while ((TmpPosition & (((ULONGLONG) 0x1ul) << CurrentHighestBit)) == 0x0ul)
    --CurrentHighestBit;  

  this->StateHighestBit[0] = CurrentHighestBit;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      TmpPosition = this->StateDescription[i];
      while ((TmpPosition & (((ULONGLONG) 0x1ul) << CurrentHighestBit)) == 0x0ul)
	--CurrentHighestBit;  
      this->StateHighestBit[i] = CurrentHighestBit;
   }

  // evaluate look-up table size
  memory /= (sizeof(int*) * 2*this->NbrLzValue);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > 2*this->NbrLzValue)
    this->MaximumLookUpShift = 2*this->NbrLzValue;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [2*this->NbrLzValue];
  this->LookUpTableShift = new int [2*this->NbrLzValue];
  for (int i = 0; i < 2*this->NbrLzValue; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentLargestBit = this->StateHighestBit[0];
  cout << this->NbrLzValue << " " << CurrentLargestBit << endl;
  int* TmpLookUpTable = this->LookUpTable[CurrentLargestBit];
  if (CurrentLargestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentLargestBit] = 0;
  else
    this->LookUpTableShift[CurrentLargestBit] = CurrentLargestBit + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentLargestBit];
  ULONGLONG CurrentLookUpTableValue = this->LookUpTableMemorySize;
  ULONGLONG TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentLargestBit != this->StateHighestBit[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
 	  CurrentLargestBit = this->StateHighestBit[i];
	  TmpLookUpTable = this->LookUpTable[CurrentLargestBit];
	  if (CurrentLargestBit < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentLargestBit] = 0;
	  else
	    this->LookUpTableShift[CurrentLargestBit] = CurrentLargestBit + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentLargestBit];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
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
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
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
#ifdef __128_BIT_LONGLONG__
  this->SignLookUpTableMask = new ULONGLONG [256];
  for (int i = 0; i < 112; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0xffff;
  for (int i = 112; i < 128; ++i)
    this->SignLookUpTableMask[i] = ((ULONGLONG) 0xffff) >> (i - 112);
  for (int i = 128; i < 256; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0;  
#else
  this->SignLookUpTableMask = new ULONGLONG [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((ULONGLONG) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0;
#endif
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// return value = Hilbert space dimension

long FermionOnSphereWithSpinLong::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalSpin > nbrFermions))
    return 0l;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz))
    return 0l;
    
  if (nbrFermions == 1) 
    if (lzMax >= totalLz)
      return 1l;
    else
      return 0l;

  if ((lzMax == 0)  && (totalLz != 0))
    return 0l;

  long Tmp = 0l;  
  if (nbrFermions > 2)    
    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1);
  else
    if ((totalLz == (2 * lzMax)) && (totalSpin == 1))
      ++Tmp;
  return  (Tmp + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, totalSpin));

}


// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphereWithSpinLong::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
							    int firstComponent, int nbrComponent)
{
  Complex Value;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereWithSpinLong::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
  
// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

RealVector FermionOnSphereWithSpinLong::ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace)
{
  RealVector FinalState(this->HilbertSpaceDimension, true);
  for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
    {
      unsigned long TmpUpState = upStateSpace.StateDescription[j];
      int TmpPos = upStateSpace.LzMax;
      while (TmpPos > 0)
	{
	  ULONGLONG Tmp = TmpUpState & (((ULONGLONG) 0x1ul) << TmpPos);
	  TmpUpState |= Tmp << TmpPos;
	  TmpUpState ^= Tmp;
	  --TmpPos;
	}
      TmpUpState <<= 1;
      double TmpComponent = upState[j];
      int Max = 63;
      while ((TmpUpState & (((ULONGLONG) 0x1ul) << Max)) == 0x0ul)
	--Max;
      int Min = 0;
      while ((TmpUpState & (((ULONGLONG) 0x1ul) << Min)) == 0x0ul)
	++Min;
      ULONGLONG TmpUpStateMask = (((ULONGLONG) 0x1ul) << Max) - 1;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	if ((this->StateDescription[i] & TmpUpState) == TmpUpState)
	  {	    
	    ULONGLONG TmpUpState3 = this->StateDescription[i] & TmpUpStateMask;
	    ULONGLONG TmpUpState2 = TmpUpState3;
#ifdef __128_BIT_LONGLONG__
	    TmpUpState3 &= 0x5555555555555555ul;
	    TmpUpState2 &= 0xaaaaaaaaaaaaaaaaul;
#else
	    TmpUpState3 &= 0x55555555ul;
	    TmpUpState2 &= 0xaaaaaaaaul;
#endif	    
	    ULONGLONG Sign = 0x0;
	    int Pos = this->LzMax << 1;
	    while ((Pos > 0) && ((TmpUpState3 & (((ULONGLONG) 0x1ul) << Pos)) == 0x0ul))
	      Pos -= 2;
	    while (Pos > 0)
	      {
		ULONGLONG TmpUpState4 = TmpUpState2 & ((((ULONGLONG) 0x1ul) << Pos) - 1ul);
#ifdef  __64_BITS__
		TmpUpState4 ^= TmpUpState4 >> 32;
#endif	
		TmpUpState4 ^= TmpUpState4 >> 16;
		TmpUpState4 ^= TmpUpState4 >> 8;
		TmpUpState4 ^= TmpUpState4 >> 4;
		TmpUpState4 ^= TmpUpState4 >> 2;
		TmpUpState4 ^= TmpUpState4 >> 1;
		Sign ^= TmpUpState4;
		Pos -= 2;
		while ((Pos > 0) && ((TmpUpState3 & (((ULONGLONG) 0x1ul) << Pos)) == 0x0ul))
		  Pos -= 2;
	      }
	    if ((Sign & ((ULONGLONG) 0x1ul)) == 0x0ul)
	      FinalState[i] = TmpComponent;
	    else
	      FinalState[i] = -TmpComponent;
	  }
    }

  for (int j = 0; j < downStateSpace.HilbertSpaceDimension; ++j)
    {
      ULONGLONG TmpDownState = downStateSpace.StateDescription[j];
      int TmpPos = downStateSpace.LzMax;
      while (TmpPos > 0)
	{
	  ULONGLONG Tmp = TmpDownState & (((ULONGLONG) 0x1ul) << TmpPos);
	  TmpDownState |= Tmp << TmpPos;
	  TmpDownState ^= Tmp;
	  --TmpPos;
	}
      double TmpComponent = downState[j];
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	if ((this->StateDescription[i] & TmpDownState) == TmpDownState)
	  {
	    FinalState[i] *= TmpComponent;
	  }
    }

  return FinalState;
}

// create a U(1) state from an SU(2) state
//
// state = vector describing the SU(2) state
// u1Space = reference on the Hilbert space associated to the U(1) state
// return value = resulting U(1) state

RealVector FermionOnSphereWithSpinLong::ForgeU1FromSU2(RealVector& state, FermionOnSphere& u1Space)
{
  RealVector FinalState(u1Space.GetHilbertSpaceDimension(), true);
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      ULONGLONG TmpState = this->StateDescription[j];
      int TmpPos = this->LzMax << 1;
      while (TmpPos >= 0)
	{
	  if ((((TmpState >> TmpPos) & (TmpState >> (TmpPos + 1))) & ((ULONGLONG) 0x1ul)) != 0x0ul)
	    TmpPos = 1;
	  TmpPos -= 2;
	}
      if (TmpPos != -1)
	{ 
	  TmpPos = 0;
	  ULONGLONG TmpState2 = 0x0ul; 
	  while (TmpPos <= this->LzMax)
	    {
	      TmpState2 |= ((TmpState & ((ULONGLONG) 0x1ul)) | ((TmpState & ((ULONGLONG) 0x2ul)) >> 1)) << TmpPos;
	      TmpState >>= 2;
	      ++TmpPos;
	    }
	  while ((TmpState2 >> TmpPos) == 0x0ul)
	    --TmpPos;
	  FinalState[u1Space.FindStateIndex((unsigned long) TmpState2, TmpPos)] += state[j];
	}
    }
  FinalState /= FinalState.Norm();
  return FinalState;  
}


//Evaluate the Density Matrix of the spin up fermions in a sector with a fixed lzUp 
//
//lzUp twice total momentum of up fermions.
//groundstate = reference on the total system groundstate
//return value = density matrix of the subsystem of spins up fermions.


RealSymmetricMatrix FermionOnSphereWithSpinLong::EvaluatePartialDensityMatrixSpinSeparation (int lzUp,RealVector & groundstate)
{
  if ((NbrFermionsUp==0)&&(lzUp==0))
    {
      RealSymmetricMatrix TmpDensityMatrix(1);
      TmpDensityMatrix.SetMatrixElement(0,0,1.0);
      return TmpDensityMatrix;
    }
  if ((NbrFermionsUp==0)&&(lzUp!=0))
    {
      RealSymmetricMatrix TmpDensityMatrix(1);
      TmpDensityMatrix.SetMatrixElement(0,0,0.0);
      return TmpDensityMatrix;
    }
  int Complementarylz=TotalLz-lzUp;
  int nbrFermionDown=NbrFermions-NbrFermionsUp;
  FermionOnSphereWithSpinLong TmpHilbertSpaceDown (nbrFermionDown,Complementarylz,LzMax,-nbrFermionDown);
  FermionOnSphereWithSpinLong TmpHilbertSpaceUp (NbrFermionsUp,lzUp,LzMax,NbrFermionsUp);
  if ((TmpHilbertSpaceUp.HilbertSpaceDimension==0)||(TmpHilbertSpaceDown.HilbertSpaceDimension==0))
    {
      RealSymmetricMatrix TmpDensityMatrix(1);
      TmpDensityMatrix.SetMatrixElement(0,0,0.0);
      return TmpDensityMatrix;
    }
  RealSymmetricMatrix TmpDensityMatrix(TmpHilbertSpaceUp.HilbertSpaceDimension);
  
  
  unsigned int perm[TmpHilbertSpaceUp.HilbertSpaceDimension][TmpHilbertSpaceDown.HilbertSpaceDimension];
  
  for(int i =0;i< TmpHilbertSpaceUp.HilbertSpaceDimension;i++)
    {
      for (int k=0;k < TmpHilbertSpaceDown.HilbertSpaceDimension;k++)
	{
	  int nbrpermutations=0;
	  int tmpnbr=NbrFermionsUp;
	  for (int p=2*LzMax+1;p>=0;p--)
	    {
	      if (((TmpHilbertSpaceUp.StateDescription[i]>>p)&1)==1)
		{
					tmpnbr--;
		}
	      nbrpermutations+=(((TmpHilbertSpaceDown.StateDescription[k])>>p)&1)*tmpnbr;
	    }
	  perm[i][k]=nbrpermutations;
	}
    }
  
  for( int i =0;i< TmpHilbertSpaceUp.HilbertSpaceDimension;i++)
    {
      for (int j=0;j< i+1; j++)
	{ 
	  double Coefficient=0.0;
	  for (int k=0;k < TmpHilbertSpaceDown.HilbertSpaceDimension;k++)
	    {
	      ULONGLONG TmpStateDescriptionRow=((TmpHilbertSpaceUp.StateDescription[i])|(TmpHilbertSpaceDown.StateDescription[k]));
	      ULONGLONG TmpStateDescriptionColomn=((TmpHilbertSpaceUp.StateDescription[j])|(TmpHilbertSpaceDown.StateDescription[k]));
	      
	      
	      int sign;
	      if(((perm[i][k]+perm[j][k])%2==1))
		{
		  sign=-1;
		}
	      else
		{
		  sign=1;
		}
	      
	      int NewLzMax=2*LzMax+1;
	      while ((TmpStateDescriptionRow >> NewLzMax) == ((ULONGLONG) 0x0ul))
		--NewLzMax;
	      int l=this->FindStateIndex(TmpStateDescriptionRow,NewLzMax);
	      NewLzMax=2*LzMax+1;
	      while ((TmpStateDescriptionColomn >> NewLzMax) == ((ULONGLONG) 0x0ul))
		--NewLzMax;
	      
	      int r=this->FindStateIndex(TmpStateDescriptionColomn,NewLzMax);
	      
	      Coefficient+=groundstate[l]*groundstate[r]*sign;
	    }
	  
	  TmpDensityMatrix.SetMatrixElement(i,j,Coefficient);
	}
    }
  return TmpDensityMatrix;
}

