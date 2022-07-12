////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//    class of fermions on sphere with SU(4) spin that allow LzMax up to      //
//                  31 (for systems with 128 bit integer support)             //
//               or 15 (on 32 bit systems without 128 bit integer support)    //
//                                                                            //
//                        last modification : 23/09/2011                      //
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
#include "HilbertSpace/FermionOnSphereWithSU4SpinLong.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <math.h>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
//

FermionOnSphereWithSU4SpinLong::FermionOnSphereWithSU4SpinLong()
{
  this->HilbertSpaceDimension = 0;
  this->NbrFermions = 0;
  this->IncNbrFermions = 0;
  this->TotalLz = 0;
  this->LzMax = 0;
  this->NbrLzValue = 0;
  this->TotalSpin = 0;
  this->TotalIsospin = 0;
  this->TotalEntanglement = 0;
  this->StateDescription = 0;
  this->StateHighestBit = 0;
  this->NbrFermionsUpPlus = 0;
  this->NbrFermionsDownPlus = 0;
  this->NbrFermionsUpMinus = 0;
  this->NbrFermionsDownMinus = 0;
  this->MaximumLookUpShift = 0;
  this->LookUpTableMemorySize = 0;
  this->LookUpTableShift = 0;
  this->LookUpTable = 0;  
  this->SignLookUpTable = 0;
  this->SignLookUpTableMask = 0;
  this->MaximumSignLookUp = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twice the total spin value
// totalIsospin = twice the total isospin value
// memory = amount of memory granted for precalculations

FermionOnSphereWithSU4SpinLong::FermionOnSphereWithSU4SpinLong (int nbrFermions, int totalLz, int lzMax, int totalSpin, int totalIsospin, 
							unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->TotalIsospin = totalIsospin;
  this->TotalEntanglement = this->IncNbrFermions;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  this->HilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
									   (this->TotalSpin + this->NbrFermions) >> 1, (this->TotalIsospin + this->NbrFermions) >> 1);
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
  this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						       (this->TotalSpin + this->NbrFermions) >> 1, (this->TotalIsospin + this->NbrFermions) >> 1, 0l);
  if (TmpHilbertSpaceDimension != this->HilbertSpaceDimension)
    {
      cout << "Mismatch in State-count and State Generation in FermionOnSphereWithSU4SpinLong!" << endl;
      exit(1);
    }
   this->GenerateLookUpTable(memory);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
//    for (int i = 0 ; i < this->HilbertSpaceDimension; ++i)
//      {
//        cout << i << " = ";
//        this->PrintState(cout, i)  << endl;
//      }
  
// #ifdef __DEBUG__
//   int UsedMemory = 0;
//   UsedMemory += this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
//   cout << "memory requested for Hilbert space = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;
//   UsedMemory = this->NbrLzValue * sizeof(int);
//   UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
//   cout << "memory requested for lookup table = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;

// #endif
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twice the total spin value
// totalIsospin = twice the total isospin value
// totalEntanglement = twice the total entanglement value
// memory = amount of memory granted for precalculations

FermionOnSphereWithSU4SpinLong::FermionOnSphereWithSU4SpinLong (int nbrFermions, int totalLz, int lzMax, int totalSpin, int totalIsospin, int totalEntanglement,
							unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->TotalIsospin = totalIsospin;
  this->TotalEntanglement = totalEntanglement;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  this->HilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
									   (this->TotalSpin + this->NbrFermions) >> 1, (this->TotalIsospin + this->NbrFermions) >> 1,
									   (this->TotalEntanglement + this->NbrFermions) >> 1);
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
  this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						       (this->TotalSpin + this->NbrFermions) >> 1, (this->TotalIsospin + this->NbrFermions) >> 1,
						       (this->TotalEntanglement + this->NbrFermions) >> 1, 0l);
//   if (TmpHilbertSpaceDimension != this->HilbertSpaceDimension)
//     {
//       cout << TmpHilbertSpaceDimension << " " << this->HilbertSpaceDimension << endl;
//       cout << "Mismatch in State-count and State Generation in FermionOnSphereWithSU4SpinLong!" << endl;
//       exit(1);
//     }
  this->HilbertSpaceDimension = TmpHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
   this->GenerateLookUpTable(memory);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  
// #ifdef __DEBUG__
//   int UsedMemory = 0;
//   UsedMemory += this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
//   cout << "memory requested for Hilbert space = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;
//   UsedMemory = this->NbrLzValue * sizeof(int);
//   UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
//   cout << "memory requested for lookup table = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;

// #endif
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereWithSU4SpinLong::FermionOnSphereWithSU4SpinLong(const FermionOnSphereWithSU4SpinLong& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->TotalIsospin = fermions.TotalIsospin;
  this->TotalEntanglement = fermions.TotalEntanglement;
  this->HighestBit = fermions.HighestBit;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnSphereWithSU4SpinLong::~FermionOnSphereWithSU4SpinLong ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      int MaxHighestBit = this->StateHighestBit[0];
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i <= MaxHighestBit; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSU4SpinLong& FermionOnSphereWithSU4SpinLong::operator = (const FermionOnSphereWithSU4SpinLong& fermions)
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
  this->TotalIsospin = fermions.TotalIsospin;
  this->TotalEntanglement = fermions.TotalEntanglement;
  this->HighestBit = fermions.HighestBit;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSU4SpinLong::Clone()
{
  return new FermionOnSphereWithSU4SpinLong(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnSphereWithSU4SpinLong::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnSphereWithSU4SpinLong::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphereWithSU4SpinLong::ExtractSubspace (AbstractQuantumNumber& q, 
							SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m_dp a_m_dp operator to a given state (only spin down isospin plus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_dm a_m_dm

double  FermionOnSphereWithSU4SpinLong::AddpAdp (int index, int m)
{
  if ((this->StateDescription[index] & (((ULONGLONG) 0x2l) << (m << 2))) != ((ULONGLONG) 0x0ul))
    return 1.0;
  else
    return 0.0;
}

// apply a^+_m_up a_m_up operator to a given state  (only spin up isospin plus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_um a_m_um

double FermionOnSphereWithSU4SpinLong::AdupAup (int index, int m)
{
  if ((this->StateDescription[index] & (((ULONGLONG) 0x8l) << (m << 2))) != ((ULONGLONG) 0x0ul))
    return 1.0;
  else
    return 0.0;
}
 
// apply a^+_m_dm a_m_dm operator to a given state (only spin down isospin minus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_dm a_m_dm

double FermionOnSphereWithSU4SpinLong::AddmAdm (int index, int m)
{
  if ((this->StateDescription[index] & (((ULONGLONG) 0x1l) << (m << 2))) != ((ULONGLONG) 0x0ul))
    return 1.0;
  else
    return 0.0;
}

// apply a^+_m_um a_m_um operator to a given state  (only spin up isospin minus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_um a_m_um

double FermionOnSphereWithSU4SpinLong::AdumAum (int index, int m)
{
  if ((this->StateDescription[index] & (((ULONGLONG) 0x4l) << (m << 2))) != ((ULONGLONG) 0x0ul))
    return 1.0;
  else
    return 0.0;
}

// apply a^+_m_up a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdupAup (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 3;
  n = (n << 2) + 3;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_up a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdupAum (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 3;
  n = (n << 2) + 2;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_up a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdupAdp (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 3;
  n = (n << 2) + 1;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_up a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
  
int FermionOnSphereWithSU4SpinLong::AdupAdm (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 3;
  n = (n << 2);
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_um a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdumAup (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 2;
  n = (n << 2) + 3;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_um a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdumAum (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 2;
  n = (n << 2) + 2;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_um a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdumAdp (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 2;
  n = (n << 2) + 1;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_um a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
  
int FermionOnSphereWithSU4SpinLong::AdumAdm (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 2;
  n = (n << 2);
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_dp a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AddpAup (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 1;
  n = (n << 2) + 3;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_dp a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AddpAum (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 1;
  n = (n << 2) + 2;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_dp a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AddpAdp (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 1;
  n = (n << 2) + 1;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_dp a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
  
int FermionOnSphereWithSU4SpinLong::AddpAdm (int index, int m, int n, double& coefficient)
{
  m = (m << 2) + 1;
  n = (n << 2);
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_dm a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AddmAup (int index, int m, int n, double& coefficient)
{
  m = (m << 2);
  n = (n << 2) + 3;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_dm a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AddmAum (int index, int m, int n, double& coefficient)
{
  m = (m << 2);
  n = (n << 2) + 2;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_dm a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AddmAdp (int index, int m, int n, double& coefficient)
{
  m = (m << 2);
  n = (n << 2) + 1;
  return this->GenericAdA(index, m, n, coefficient);
}

// apply a^+_m_dm a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
  
int FermionOnSphereWithSU4SpinLong::AddmAdm (int index, int m, int n, double& coefficient)
{
  m = (m << 2);
  n = (n << 2);
  return this->GenericAdA(index, m, n, coefficient);
}


// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereWithSU4SpinLong::FindStateIndex(ULONGLONG stateDescription, int lzmax)
{
  ULONGLONG CurrentState = stateDescription >> this->LookUpTableShift[lzmax];
  int PosMin = this->LookUpTable[lzmax][CurrentState];
  int PosMax = this->LookUpTable[lzmax][CurrentState+ 1];
  int PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->StateDescription[PosMid];
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

ostream& FermionOnSphereWithSU4SpinLong::PrintState (ostream& Str, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  ULONGLONG Tmp;
  Str << " | ";
  for (int i = this->NbrLzValue-1; i >=0 ; --i)
    {
      Tmp = ((TmpState >> (i << 2)) & ((ULONGLONG) 0xful));
      if (Tmp & ((ULONGLONG) 0x8ul))
	Str << "u+ ";
      else
	Str << "0 ";
      if (Tmp & ((ULONGLONG) 0x4ul))
	Str << "u- ";
      else
	Str << "0 ";
      if (Tmp & ((ULONGLONG) 0x2ul))
	Str << "d+ ";
      else
	Str << "0 ";
      if (Tmp & ((ULONGLONG) 0x1ul))
	Str << "d- ";
      else
	Str << "0 ";
      Str << "| ";
    }
   return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSU4SpinLong::GenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin, long pos)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalIsospin < 0) ||  
      (totalSpin > nbrFermions) || (totalIsospin > nbrFermions))
    return pos;
  if ((nbrFermions == 0) && (totalLz == 0) && (totalSpin == 0) && (totalIsospin == 0))
      {
	this->StateDescription[pos] = ((ULONGLONG) 0x0ul);
	return (pos + 1l);
      }
    
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < totalIsospin) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalIsospin)) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz) || ((((2 * lzMax + nbrFermions + 1 - totalIsospin) * nbrFermions) >> 1) < totalLz))
    return pos;
    
  if (nbrFermions == 1) 
    {
      if (lzMax >= totalLz)
	{
	  this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << ((totalLz << 2) + (totalSpin << 1) + totalIsospin);
	  return (pos + 1l);
	}
      else
	return pos;
    }

  if ((lzMax == 0)  && (totalLz != 0))
    return pos;


  long TmpPos;
  ULONGLONG Mask;
  if (nbrFermions >= 4)
    {
      TmpPos = this->GenerateStates(nbrFermions - 4, lzMax - 1, totalLz - (lzMax << 2), totalSpin - 2, totalIsospin - 2, pos);
      Mask = ((ULONGLONG) 0xful) << ( lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if (nbrFermions >= 3)
    {
      TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 2,  pos);
      Mask = ((ULONGLONG) 0xeul) << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
      TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 1,  pos);
      Mask = ((ULONGLONG) 0xdul) << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 2, totalIsospin - 1,  pos);
  Mask = ((ULONGLONG) 0xcul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  if (nbrFermions >= 3)
    {
      TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 2,  pos);
      Mask = ((ULONGLONG) 0xbul) << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 1, totalIsospin - 2,  pos);
  Mask = ((ULONGLONG) 0xaul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 1, totalIsospin - 1,  pos);
  Mask = ((ULONGLONG) 0x9ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin - 1,  pos);
  Mask = ((ULONGLONG) 0x8ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  if (nbrFermions >= 3)
    {
      TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 1,  pos);
      Mask = ((ULONGLONG) 0x7ul) << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 1, totalIsospin - 1,  pos);
  Mask = ((ULONGLONG) 0x6ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 1, totalIsospin,  pos);
  Mask = ((ULONGLONG) 0x5ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin,  pos);
  Mask = ((ULONGLONG) 0x4ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin, totalIsospin - 1,  pos);
  Mask = ((ULONGLONG) 0x3ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin - 1,  pos);
  Mask = ((ULONGLONG) 0x2ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin,  pos);
  Mask = ((ULONGLONG) 0x1ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  return this->GenerateStates(nbrFermions, lzMax - 1, totalLz, totalSpin, totalIsospin, pos);
};


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// totalEntanglement = number of particles with entanglement plus
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSU4SpinLong::GenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin, int totalEntanglement, long pos)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalIsospin < 0) ||  (totalEntanglement < 0) ||  
      (totalSpin > nbrFermions) || (totalIsospin > nbrFermions) || (totalEntanglement > nbrFermions))
    return pos;
  if ((nbrFermions == 0) && (totalLz == 0) && (totalSpin == 0) && (totalIsospin == 0) && (totalEntanglement == 0))
      {
	this->StateDescription[pos] = ((ULONGLONG) 0x0ul);
	return (pos + 1l);
      }
    
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < totalIsospin) || ((2 * (lzMax + 1)) < totalEntanglement) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalIsospin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions -totalEntanglement )) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz) 
      || ((((2 * lzMax + nbrFermions + 1 - totalIsospin) * nbrFermions) >> 1) < totalLz)
      || ((((2 * lzMax + nbrFermions + 1 - totalEntanglement) * nbrFermions) >> 1) < totalLz))
    return pos;
    
  if (nbrFermions == 1) 
    {
      if ((lzMax >= totalLz) && (totalEntanglement != (totalSpin ^ totalIsospin)))
	{
	  this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << ((totalLz << 2) + (totalSpin << 1) + totalIsospin);
	  return (pos + 1l);
	}
      else
	return pos;
    }

  if ((lzMax == 0)  && (totalLz != 0))
    return pos;


  long TmpPos;
  ULONGLONG Mask;
  if (nbrFermions >= 4)
    {
      TmpPos = this->GenerateStates(nbrFermions - 4, lzMax - 1, totalLz - (lzMax << 2), totalSpin - 2, totalIsospin - 2, totalEntanglement - 2, pos);
      Mask = ((ULONGLONG) 0xful) << ( lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if (nbrFermions >= 3)
    {
      TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 2, totalEntanglement - 1, pos);
      Mask = ((ULONGLONG) 0xeul) << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
      TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 1, totalEntanglement - 2, pos);
      Mask = ((ULONGLONG) 0xdul) << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 2, totalIsospin - 1, totalEntanglement - 1, pos);
  Mask = ((ULONGLONG) 0xcul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  if (nbrFermions >= 3)
    {
      TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 2, totalEntanglement - 2, pos);
      Mask = ((ULONGLONG) 0xbul) << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 1, totalIsospin - 2, totalEntanglement - 1, pos);
  Mask = ((ULONGLONG) 0xaul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 1, totalIsospin - 1, totalEntanglement - 2, pos);
  Mask = ((ULONGLONG) 0x9ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin - 1, totalEntanglement - 1, pos);
  Mask = ((ULONGLONG) 0x8ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  if (nbrFermions >= 3)
    {
      TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 1, totalEntanglement - 1, pos);
      Mask = ((ULONGLONG) 0x7ul) << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 1, totalIsospin - 1, totalEntanglement, pos);
  Mask = ((ULONGLONG) 0x6ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 1, totalIsospin, totalEntanglement - 1, pos);
  Mask = ((ULONGLONG) 0x5ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin, totalEntanglement, pos);
  Mask = ((ULONGLONG) 0x4ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin, totalIsospin - 1, totalEntanglement - 1, pos);
  Mask = ((ULONGLONG) 0x3ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin - 1, totalEntanglement, pos);
  Mask = ((ULONGLONG) 0x2ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin, totalEntanglement - 1, pos);
  Mask = ((ULONGLONG) 0x1ul) << (lzMax << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  return this->GenerateStates(nbrFermions, lzMax - 1, totalLz, totalSpin, totalIsospin, totalEntanglement, pos);
};


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereWithSU4SpinLong::GenerateLookUpTable(unsigned long memory)
{
  // get every highest bit poisition
  ULONGLONG TmpPosition = this->StateDescription[0];
#ifdef __128_BIT_LONGLONG__
  int CurrentHighestBit = 127;
#else
  int CurrentHighestBit = 63;
#endif
  while ((TmpPosition & (((ULONGLONG) 0x1ul) << CurrentHighestBit)) == ((ULONGLONG) 0x0ul))
    --CurrentHighestBit;  
  int MaxHighestBit = CurrentHighestBit;
  this->StateHighestBit[0] = CurrentHighestBit;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      TmpPosition = this->StateDescription[i];
      while ((TmpPosition & (((ULONGLONG) 0x1ul) << CurrentHighestBit)) == ((ULONGLONG) 0x0ul))
	--CurrentHighestBit;  
      this->StateHighestBit[i] = CurrentHighestBit;
   }

  // evaluate look-up table size
  memory /= (sizeof(int*) * MaxHighestBit);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > MaxHighestBit)
    this->MaximumLookUpShift = MaxHighestBit;
  this->LookUpTableMemorySize = ((ULONGLONG) 0x1ul) << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [MaxHighestBit + 1];
  this->LookUpTableShift = new int [MaxHighestBit + 1];
  for (int i = 0; i <= MaxHighestBit; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];

  CurrentHighestBit = this->StateHighestBit[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
  if (CurrentHighestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentHighestBit] = 0;
  else
    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentHighestBit];
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
      if (CurrentHighestBit != this->StateHighestBit[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
 	  CurrentHighestBit = this->StateHighestBit[i];
	  TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
	  if (CurrentHighestBit < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentHighestBit] = 0;
	  else
	    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentHighestBit];
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
  while (CurrentLookUpTableValue > ((ULONGLONG) 0x0ul))
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
    this->SignLookUpTableMask[i] = (ULONGLONG) 0xfffful;
  for (int i = 112; i < 128; ++i)
    this->SignLookUpTableMask[i] = ((ULONGLONG) 0xfffful) >> (i - 112);
  for (int i = 128; i < 256; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0x0ul;  
#else
  this->SignLookUpTableMask = new ULONGLONG [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0xfffful;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((ULONGLONG) 0xfffful) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0x0ul;
#endif
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// return value = Hilbert space dimension

long FermionOnSphereWithSU4SpinLong::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalIsospin < 0) ||  
      (totalSpin > nbrFermions) || (totalIsospin > nbrFermions))
    return 0l;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < totalIsospin) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalIsospin)) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz) || ((((2 * lzMax + nbrFermions + 1 - totalIsospin) * nbrFermions) >> 1) < totalLz))
    return 0l;
    
  if (nbrFermions == 1) 
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }

  if ((lzMax == 0)  && (totalLz != 0))
    return 0l;

  ULONGLONG Tmp = 0l;
  if (nbrFermions >= 3)    
    {
      Tmp += (this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 2, totalIsospin - 1)
	      + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 2)
	      + (2l * this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 1))
	      + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin)
	      + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin, totalIsospin - 1));

      if (nbrFermions > 3)
	{
	  Tmp += (this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 2)
		  + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 1)
		  + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 2)
		  + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 1));
	  if (nbrFermions == 4)
	    {
	      if ((totalLz == (4 * lzMax)) && (totalSpin == 2) && (totalIsospin == 2))
		++Tmp;      
	    }
	  else
	    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 4, lzMax - 1, totalLz - (4 * lzMax), totalSpin - 2, totalIsospin - 2);
	}
      else
	if ((totalLz == (3 * lzMax)) && (((totalSpin == 2) || (totalSpin == 1)) && ((totalIsospin == 2) || (totalIsospin == 1))))
	  ++Tmp;
    }
  else
    if (totalLz == (2 * lzMax))
      {
 	switch (totalSpin)
 	  {
 	  case 2:
	    if (totalIsospin == 1)
	      ++Tmp;
 	    break;
 	  case 1:
	    switch (totalIsospin)
	      {
	      case 2:
		++Tmp;
		break;
	      case 1:
		Tmp += 2l;
		break;
	      case 0:
		++Tmp;
		break;
	      }
	    break;
 	  case 0:
	    if (totalIsospin == 1) 
	      ++Tmp;
	    break; 
 	  }
      }

  return  (Tmp + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin - 1)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin - 1)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, totalSpin, totalIsospin));

}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// totalIsospin = number of particles with isospin plus
// totalEntanglement = number of particles with entanglement plus
// return value = Hilbert space dimension

long FermionOnSphereWithSU4SpinLong::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin, int totalEntanglement)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalIsospin < 0) || (totalEntanglement < 0) ||
      (totalSpin > nbrFermions) || (totalIsospin > nbrFermions) || (totalEntanglement > nbrFermions))
    return 0l;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpin) || ((2 * (lzMax + 1)) < totalIsospin) || ((2 * (lzMax + 1)) < totalEntanglement) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalSpin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalIsospin)) 
      || ((2 * (lzMax + 1)) < (nbrFermions - totalEntanglement)) 
      || ((((2 * lzMax + nbrFermions + 1 - totalSpin) * nbrFermions) >> 1) < totalLz) 
      || ((((2 * lzMax + nbrFermions + 1 - totalIsospin) * nbrFermions) >> 1) < totalLz)
      || ((((2 * lzMax + nbrFermions + 1 - totalEntanglement) * nbrFermions) >> 1) < totalLz))
    return 0l;
    
  if ((nbrFermions == 0) && (totalLz == 0) && (totalSpin == 0) && (totalIsospin == 0) && (totalEntanglement == 0))
    return 1l;
  if (nbrFermions == 1) 
    {
      if ((lzMax >= totalLz) && (totalEntanglement != (totalSpin ^ totalIsospin)))
	return 1l;
      else
	return 0l;
    }

  if ((lzMax == 0)  && (totalLz != 0))
    return 0l;

  ULONGLONG Tmp = 0l;
  if (nbrFermions >= 3)    
    {
      Tmp += (this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 2, totalIsospin - 1, totalEntanglement - 1)
	      + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 2, totalEntanglement - 1)
	      + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 1, totalEntanglement - 2)
	      + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin - 1, totalEntanglement)
	      + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, totalIsospin, totalEntanglement - 1)
	      + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin, totalIsospin - 1, totalEntanglement - 1));
      
      if (nbrFermions > 3)
	{
 	  Tmp += (this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 2, totalEntanglement - 1)
 		  + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 1, totalEntanglement - 1)
 		  + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, totalIsospin - 2, totalEntanglement - 2)
 		  + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, totalIsospin - 1, totalEntanglement - 2));
 	  if (nbrFermions == 4)
 	    {
 	      if ((totalLz == (4 * lzMax)) && (totalSpin == 2) && (totalIsospin == 2) && (totalEntanglement == 2))
 		++Tmp;      
 	    }
 	  else
 	    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 4, lzMax - 1, totalLz - (4 * lzMax), totalSpin - 2, totalIsospin - 2, totalEntanglement - 2);
 	}
      else
	if ((totalLz == (3 * lzMax)) && (((totalSpin * totalIsospin * totalEntanglement) == 4) || ((totalSpin * totalIsospin * totalEntanglement) == 1)))
 	  ++Tmp;
    }
  else
    if (totalLz == (2 * lzMax))
      {
 	switch (totalSpin)
 	  {
 	  case 2:
	    if ((totalIsospin == 1) && (totalEntanglement == 1))
	      ++Tmp;
 	    break;
 	  case 1:
	    switch (totalIsospin)
	      {
	      case 2:
		if (totalEntanglement == 1)
		  ++Tmp;
		break;
	      case 1:
		if (totalEntanglement != 1)
		  ++Tmp;
		break;
	      case 0:
		if (totalEntanglement == 1)
		  ++Tmp;
		break;
	      }
	    break;
 	  case 0:
	    if ((totalIsospin == 1)  && (totalEntanglement == 1))
	      ++Tmp;
	    break; 
 	  }
      }

  return  (Tmp + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin, totalEntanglement)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin - 1, totalEntanglement)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, totalIsospin - 1, totalEntanglement - 1)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, totalIsospin, totalEntanglement - 1)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, totalSpin, totalIsospin, totalEntanglement));

}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = twce the total spin value
// return value = Hilbert space dimension

int FermionOnSphereWithSU4SpinLong::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin, int totalIsospin)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + (nbrFermions * lzMax)) >> 1, 
						    (totalSpin + nbrFermions) >> 1, (totalIsospin + nbrFermions) >> 1);
}


// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphereWithSU4SpinLong::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
							    int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
  ComplexMatrix Slatter(this->NbrFermions, this->NbrFermions);
  ComplexMatrix Functions(this->LzMax + 1, this->NbrFermions);
  RealVector TmpCoordinates(2);
  int* Indices = new int [this->NbrFermions];
  int Pos;
  int Lz;
  for (int j = 0; j < this->NbrFermions; ++j)
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
  double Factor = 1.0;
  for (int i = 2; i <= this->NbrFermions; ++i)
    Factor *= (double) i;
  Factor = 1.0 / sqrt(Factor);
  ULONGLONG TmpStateDescription;
  int LastComponent = firstComponent + nbrComponent;
  for (int k = firstComponent; k < LastComponent; ++k)
    {
      Pos = 0;
      Lz = 0;
      TmpStateDescription = this->StateDescription[k];
      while (Pos < this->NbrFermions)
	{
	  if ((TmpStateDescription & ((ULONGLONG) 0x3l)) != ((ULONGLONG) 0x0l))
	    {
	      Indices[Pos] = Lz;
	      ++Pos;
	    }
	  ++Lz;
	  TmpStateDescription >>= 2;
	}
      for (int i = 0; i < this->NbrFermions; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < this->NbrFermions; ++j)
	    {
	      Slatter[i].Re(j) = TmpColum2.Re(Indices[j]);
	      Slatter[i].Im(j) = TmpColum2.Im(Indices[j]);
	    }
	}
      Complex SlatterDet = Slatter.Determinant();
      Value += SlatterDet * (state[k] * Factor);
    }
  delete[] Indices;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereWithSU4SpinLong::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}

// apply a_n1_up a_n2_up operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereWithSU4SpinLong::AupAup (int index, int n1, int n2)
{
  n1 <<= 2;
  n1 += 3;
  n2 <<= 2;
  n2 += 3;
  return this->GenericAA(index, n1, n2);
}

// apply a_n1_up a_n2_um operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereWithSU4SpinLong::AupAum (int index, int n1, int n2)
{
  n1 <<= 2;
  n1 += 3;
  n2 <<= 2;
  n2 += 2;
  return this->GenericAA(index, n1, n2);
}

// apply a_n1_up a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereWithSU4SpinLong::AupAdp (int index, int n1, int n2)
{
  n1 <<= 2;
  n1 += 3;
  n2 <<= 2;
  ++n2;
  return this->GenericAA(index, n1, n2);
}

// apply a_n1_up a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereWithSU4SpinLong::AupAdm (int index, int n1, int n2)
{
  n1 <<= 2;
  n1 += 3;
  n2 <<= 2;
  return this->GenericAA(index, n1, n2);
}

// apply a_n1_um a_n2_um operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereWithSU4SpinLong::AumAum (int index, int n1, int n2)
{
  n1 <<= 2;
  n1 += 2;
  n2 <<= 2;
  n2 += 2;
  return this->GenericAA(index, n1, n2);
}

// apply a_n1_um a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereWithSU4SpinLong::AumAdp (int index, int n1, int n2)
{
  n1 <<= 2;
  n1 += 2;
  n2 <<= 2;
  ++n2;
  return this->GenericAA(index, n1, n2);
}

// apply a_n1_um a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereWithSU4SpinLong::AumAdm (int index, int n1, int n2)
{
  n1 <<= 2;
  n1 += 2;
  n2 <<= 2;
  return this->GenericAA(index, n1, n2);
}

// apply a_n1_dp a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereWithSU4SpinLong::AdpAdp (int index, int n1, int n2)
{
  n1 <<= 2;
  ++n1;
  n2 <<= 2;
  ++n2;
  return this->GenericAA(index, n1, n2);
}

// apply a_n1_dp a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereWithSU4SpinLong::AdpAdm (int index, int n1, int n2)
{
  n1 <<= 2;
  ++n1;
  n2 <<= 2;
  return this->GenericAA(index, n1, n2);
}

// apply a_n1_dm a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereWithSU4SpinLong::AdmAdm (int index, int n1, int n2)
{
  n1 <<= 2;
  n2 <<= 2;
  return this->GenericAA(index, n1, n2);
}

// apply a^+_m1_up a^+_m2_up operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdupAdup (int m1, int m2, double& coefficient)
{
  m1 <<= 2;
  m1 += 3;
  m2 <<= 2;
  m2 += 3;
  return this->GenericAdAd(m1, m2, coefficient);
}

// apply a^+_m1_up a^+_m2_um operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdupAdum (int m1, int m2, double& coefficient)
{
  m1 <<= 2;
  m1 += 3;
  m2 <<= 2;
  m2 += 2;
  return this->GenericAdAd(m1, m2, coefficient);
}

// apply a^+_m1_up a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdupAddp (int m1, int m2, double& coefficient)
{
  m1 <<= 2;
  m1 += 3;
  m2 <<= 2;
  ++m2;
  return this->GenericAdAd(m1, m2, coefficient);
}

// apply a^+_m1_up a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdupAddm (int m1, int m2, double& coefficient)
{
  m1 <<= 2;
  m1 += 3;
  m2 <<= 2;
  return this->GenericAdAd(m1, m2, coefficient);
}

// apply a^+_m1_um a^+_m2_um operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdumAdum (int m1, int m2, double& coefficient)
{
  m1 <<= 2;
  m1 += 2;
  m2 <<= 2;
  m2 += 2;
  return this->GenericAdAd(m1, m2, coefficient);
}

// apply a^+_m1_um a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdumAddp (int m1, int m2, double& coefficient)
{
  m1 <<= 2;
  m1 += 2;
  m2 <<= 2;
  ++m2;
  return this->GenericAdAd(m1, m2, coefficient);
}

// apply a^+_m1_um a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AdumAddm (int m1, int m2, double& coefficient)
{
  m1 <<= 2;
  m1 += 2;
  m2 <<= 2;
  return this->GenericAdAd(m1, m2, coefficient);
}

// apply a^+_m1_dp a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AddpAddp (int m1, int m2, double& coefficient)
{
  m1 <<= 2;
  ++m1;
  m2 <<= 2;
  ++m2;
  return this->GenericAdAd(m1, m2, coefficient);
}

// apply a^+_m1_dp a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AddpAddm (int m1, int m2, double& coefficient)
{
  m1 <<= 2;
  ++m1;
  m2 <<= 2;
  return this->GenericAdAd(m1, m2, coefficient);
}

// apply a^+_m1_dm a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSU4SpinLong::AddmAddm (int m1, int m2, double& coefficient)
{
  m1 <<= 2;
  m2 <<= 2;
  return this->GenericAdAd(m1, m2, coefficient);
}


// create a U(1) state from an SU(4) state
//
// state = vector describing the SU(4) state
// u1Space = reference on the Hilbert space associated to the U(1) state
// return value = resulting U(1) state

RealVector FermionOnSphereWithSU4SpinLong::ForgeU1FromSU4(RealVector& state, FermionOnSphere& u1Space)
{
  RealVector FinalState(u1Space.GetHilbertSpaceDimension(), true);
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      ULONGLONG TmpState = this->StateDescription[j];
      ULONGLONG TmpState2 = TmpState; 
      int TmpPos = this->LzMax << 2;
      while (TmpPos >=0)
	{
	  ULONGLONG  TmpNbrParticles = TmpState2 & ((ULONGLONG) 0x1ul);
	  TmpState2 >>= 1;
	  TmpNbrParticles += TmpState2 & ((ULONGLONG) 0x1ul);
	  TmpState2 >>= 1;
	  TmpNbrParticles += TmpState2 & ((ULONGLONG) 0x1ul);
	  TmpState2 >>= 1;
	  TmpNbrParticles += TmpState2 & ((ULONGLONG) 0x1ul);
	  TmpState2 >>= 1;
	  if (TmpNbrParticles > ((ULONGLONG) 0x1ul))
	    TmpPos = 1;
	  TmpPos -= 4;
	}
      if (TmpPos != -3)
	{ 
	  TmpPos = 0;
	  TmpState2 = 0x0ul; 
	  while (TmpPos <= this->LzMax)
	    {
	      TmpState2 |= ((TmpState & ((ULONGLONG) 0x1ul)) | ((TmpState & 0x2ul) >> 1) | ((TmpState & 0x4ul) >> 2) | ((TmpState & 0x8ul) >> 3)) << TmpPos;
	      TmpState >>= 4;
	      ++TmpPos;
	    }
	  while ((TmpState2 >> TmpPos) == 0x0ul)
	    --TmpPos;
	  FinalState[u1Space.FindStateIndex(TmpState2, TmpPos)] += state[j];
	}
    }
  FinalState /= FinalState.Norm();
  return FinalState;  
}

// create a SU(2) state from an SU(4) state (fusing same spin values,i.e symmetrizing over the isospin)
//
// state = vector describing the SU(4) state
// su2Space = reference on the Hilbert space associated to the SU(2) state
// return value = resulting SU(2) state

RealVector FermionOnSphereWithSU4SpinLong::ForgeSU2FromSU4(RealVector& state, FermionOnSphereWithSpin& su2Space)
{
  RealVector FinalState(su2Space.GetHilbertSpaceDimension(), true);
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      ULONGLONG TmpState = this->StateDescription[j];
      ULONGLONG TmpState2 = TmpState; 
      int TmpPos = this->LzMax << 2;
      while (TmpPos >=0)
	{
	  if ((((TmpState2 >> 3) & (TmpState2 >> 2) & ((ULONGLONG) 0x1ul)) != 0x0ul) ||
	      (((TmpState2 >> 1) & TmpState2 & ((ULONGLONG) 0x1ul)) != 0x0ul))
	    TmpPos = 1;
	  TmpState2 >>= 4;
	  TmpPos -= 4;
	}
      if (TmpPos != -3)
	{ 
	  TmpPos = 0;
	  TmpState2 = 0x0ul; 
	  int TmpLzMax = this->LzMax << 1;
	  while (TmpPos <= TmpLzMax)
	    {
	      TmpState2 |= ((TmpState | (TmpState >> 1)) & ((ULONGLONG) 0x1ul)) << TmpPos;
	      ++TmpPos;
	      TmpState2 |= (((TmpState >> 2) | (TmpState >> 3)) & ((ULONGLONG) 0x1ul)) << TmpPos;
	      TmpState >>= 4;
	      ++TmpPos;
	    }
	  while ((TmpState2 >> TmpPos) == 0x0ul)
	    --TmpPos;
	  FinalState[su2Space.FindStateIndex(TmpState2, TmpPos)] += state[j];
	}
    }
  FinalState /= FinalState.Norm();
  return FinalState;  
}
