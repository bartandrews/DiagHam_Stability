////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin without            //
//                            sign precalculation table                       //
//                                                                            //
//                        last modification : 12/12/2005                      //
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
#include "HilbertSpace/FermionOnSphereWithSpin.h"
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

FermionOnSphereWithSpin::FermionOnSphereWithSpin()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twce the total spin value
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpin::FermionOnSphereWithSpin (int nbrFermions, int totalLz, int lzMax, int totalSpin, unsigned long memory)
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
//   this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz, this->TotalSpin);
//   long TmpBidule = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
// 							      (this->TotalSpin + this->NbrFermions) >> 1);
  this->LargeHilbertSpaceDimension = (int) this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
										 (this->TotalSpin + this->NbrFermions) >> 1);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
//   if (this->GenerateStates(this->NbrFermions, this->LzMax, this->TotalLz, this->TotalSpin) != this->HilbertSpaceDimension)
//     {
//       cout << "Mismatch in State-count and State Generation in FermionOnSphereWithSpin!" << endl;
//       exit(1);
//     }

  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						     (this->TotalSpin + this->NbrFermions) >> 1, 0l);
  this->GenerateLookUpTable(memory);

//   for (int i=0; i<HilbertSpaceDimension; ++i)
//     PrintState(cout, i)<<endl;
  
#ifdef __DEBUG__
  long UsedMemory = 0;
  UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
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

FermionOnSphereWithSpin::FermionOnSphereWithSpin(const FermionOnSphereWithSpin& fermions)
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

FermionOnSphereWithSpin::~FermionOnSphereWithSpin ()
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
      delete[] this->SignLookUpTableMask;
      delete[] this->SignLookUpTable;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpin& FermionOnSphereWithSpin::operator = (const FermionOnSphereWithSpin& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
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
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpin::Clone()
{
  return new FermionOnSphereWithSpin(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnSphereWithSpin::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnSphereWithSpin::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphereWithSpin::ExtractSubspace (AbstractQuantumNumber& q, 
							SubspaceSpaceConverter& converter)
{
  return 0;
}


// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// m = Lz value of particle to be added
// s = spin index of particle to be added (0=down, 1=up)
// coefficient = reference on the double where the multiplicative factor has to be stored
unsigned long FermionOnSphereWithSpin::Ad (unsigned long state, int m, int s, double& coefficient)
{
  m = (m<<1) + (s&1);
  if ((state & (0x1ul << m)) != 0)
    {
      coefficient=0.0;
      return 0x0l;
    }
  int NewHighestBit = getHighestBit(state)-1;
  coefficient = 1.0;
  if (m > NewHighestBit)
    NewHighestBit = m;
  else
    {
      coefficient *= this->SignLookUpTable[(state >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(state >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(state >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(state >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  state |= (0x1ul << m);

  return state;
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

//#include <bitset>
int FermionOnSphereWithSpin::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  unsigned long signs = 0x0l;
  //bitset<32> tmpB = State;
  n1 = (n1<<1) + 1;
  n2 = (n2<<1) + 1;  
  //cout << "Examining uuuu: " << tmpB << " LzMax: " << StateHighestBit <<" for n's = (" << n1 << ", "<< n2 << ") coeff: " << coefficient << endl;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & (0x1l << n1)) == 0) 
      || ((State & (0x1l << n2)) == 0) || (n1 == n2) || (m1 == m2)) 
    {
      coefficient = 0.0;
      //cout << "First exit" << endl;
      return this->HilbertSpaceDimension;
    }
  // evaluate bit positions corresponding to (m1,up), (m2,up)
  m1 = (m1<<1) + 1;
  m2 = (m2<<1) + 1;
  int NewLargestBit = StateHighestBit;
  //cout << " m's: (" << m1 << ", " << m2 << ")" << endl;
  signs = State & ((0x1l<<n2) -1); // & mask with all bits set at positions right of n2
  State &= ~(0x1l << n2);
  
  /*tmpB = State;
  cout << "Unset n2:  " << tmpB << endl;
  tmpB = signs;
  cout << "signs:     " << tmpB << endl;*/
  signs ^= State & ((0x1l<<n1) -1);
  State &= ~(0x1l << n1);
  /*tmpB = State;
  cout << "Unset n1:  " << tmpB  << endl;
  tmpB = signs;
  cout << "signs:     " << tmpB << endl;*/
  
  // test if possible to create particles at m1, m2:
  if (((State & (0x1l << m2))!= 0)|| ((State & (0x1l << m1)) != 0))
    {
      //cout << "Third exit" << endl;
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  // recalculate NewLargestBit taking into account the above operations:
  
  if ((NewLargestBit == n2)|| (NewLargestBit == n1))
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  signs ^= State & ((0x1l<<m2)-1);
  State |= (0x1l << m2);
  /*tmpB = State;
  cout << "Set m2:    " << tmpB << endl;
  tmpB = signs;
  cout << "signs:     " << tmpB << endl;*/

  if (m1 > NewLargestBit)
    {
      NewLargestBit = m1;
    }

  // in ParticleOnSphereWithSpin Hamiltonian, always m2>m1! -> we can leave these lines out!
  if (m2 > NewLargestBit) 
    { 
      NewLargestBit = m2; 
    } 
  
  
  signs ^= State & ((0x1l<<m1)-1);
  State |= (0x1l << m1);
  /*tmpB = State;
  cout << "Set m1:    " << tmpB << endl;
  tmpB = signs;
  cout << "signs:     " << tmpB << endl;*/
  
  coefficient = ComputeSign (signs);
  //  cout << "Non-zero result found: "  << coefficient << " New Largest Bit: " << NewLargestBit <<endl;
  return this->FindStateIndex(State, NewLargestBit);
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

int FermionOnSphereWithSpin::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  unsigned long signs = 0x0l;
  n1 <<= 1;
  n2 <<= 1;  
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & (0x1l << n1)) == 0) 
      || ((State & (0x1l << n2)) == 0) || (n1 == n2) || (m1 == m2)) // the last two are superflous, though adding some security
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  // evaluate bit positions corresponding to (m1,down), (m2,down)
  m1 <<= 1;
  m2 <<= 1;
  int NewLargestBit = StateHighestBit;
  
  signs = State & ((0x1l<<n2) -1); // & mask with all bits set at positions right of n2
  State &= ~(0x1l << n2);
  
  signs ^= State & ((0x1l<<n1) -1);
  State &= ~(0x1l << n1);

  // test if possible to create particles at m1, m2:
  if (((State & (0x1l << m2))!= 0)|| ((State & (0x1l << m1)) != 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  // recalculate NewLargestBit taking into account the above operations:
  
  if ((NewLargestBit == n2)|| (NewLargestBit == n1))
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  signs ^= State & ((0x1l<<m2)-1);
  State |= (0x1l << m2);
  
  if (m1 > NewLargestBit)
    {
      NewLargestBit = m1;
    }

  // if called from ParticleOnSphereWithSpin... Hamiltonian, always m1>m2!
   if (m2 > NewLargestBit)
    {
      NewLargestBit = m2;
    }


  signs ^= State & ((0x1l<<m1)-1);
  State |= (0x1l << m1);

  coefficient = ComputeSign (signs);
  
  return this->FindStateIndex(State, NewLargestBit);

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

int FermionOnSphereWithSpin::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  unsigned long signs = 0x0l;
  n1 <<= 1;
  n2 = (n2<<1) + 1;  
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & (0x1l << n1)) == 0) 
      || ((State & (0x1l << n2)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  // evaluate bit positions corresponding to (m1,up), (m2,up)
  m1 <<= 1;
  m2 = (m2<<1) + 1;
  int NewLargestBit = StateHighestBit;
  
  signs = State & ((0x1l<<n2) -1); // & mask with all bits set at positions right of n2
  State &= ~(0x1l << n2);
  
  signs ^= State & ((0x1l<<n1) -1);
  State &= ~(0x1l << n1);

  // test if possible to create particles at m1, m2:
  if (((State & (0x1l << m2))!= 0)|| ((State & (0x1l << m1)) != 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  
  // recalculate NewLargestBit taking into account the above operations:
  
  if ((NewLargestBit == n2)|| (NewLargestBit == n1))
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if (m2 > NewLargestBit)
    {
      NewLargestBit = m2;
    }

  signs ^= State & ((0x1l<<m2)-1);
  State |= (0x1l << m2);
  
  if (m1 > NewLargestBit)
    {
      NewLargestBit = m1;
    }

  signs ^= State & ((0x1l<<m1)-1);
  State |= (0x1l << m1);

  coefficient = ComputeSign (signs);
  
  return this->FindStateIndex(State, NewLargestBit);

}


// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereWithSpin::AduAu (int index, int m)
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

double FermionOnSphereWithSpin::AddAd (int index, int m)
{
  if ((this->StateDescription[index] & (0x1l << (m << 1))) != 0)
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

double FermionOnSphereWithSpin::AuAu (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  ++n1;
  n2 <<= 1;
  ++n2;
  if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
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

double FermionOnSphereWithSpin::AdAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  n2 <<= 1;
  if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
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

double FermionOnSphereWithSpin::AuAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  ++n1;
  n2 <<= 1;
  if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ( ((this->ProdATemporaryState >> this->ProdALzMax) == 0) && (this->ProdALzMax>0))
    --this->ProdALzMax;
  return Coefficient;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::AduAdu (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  ++m2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
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
  TmpState |= (0x1ul << m2);
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
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::AddAdd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  m2 <<= 1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
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
  TmpState |= (0x1ul << m2);
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
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::AduAdd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0))
    return this->HilbertSpaceDimension;
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
  TmpState |= (0x1ul << m2);
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
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereWithSpin::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  this->ProdALzMax = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (n[i] << 1) + spinIndices[i];
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

double FermionOnSphereWithSpin::ProdA (int index, int* n, int spinIndices, int nbrIndices)
{
  this->ProdALzMax = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (n[i] << 1) + ((spinIndices >> i) & 0x1);
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

int FermionOnSphereWithSpin::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (m[i] << 1) + spinIndices[i];
      if ((TmpState & (0x1l << Index)) != 0)
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
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (0x1l << Index);
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

int FermionOnSphereWithSpin::ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (m[i] << 1) + ((spinIndices >> i) & 0x1);
      if ((TmpState & (0x1l << Index)) != 0)
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
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (0x1l << Index);
    }
  return this->FindStateIndex(TmpState, NewLzMax);
}

// get the variance of the state
// index = index of state to consider
int FermionOnSphereWithSpin::StateVariance (int index)
{
  int MyStateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  int var=0;
  for (int b=0; b<=MyStateHighestBit; ++b)
    if  ( (State & (0x1ul << b)) != 0x0ul)
      var += (((b>>1)<<1)-LzMax)*(((b>>1)<<1)-LzMax);
  return var;
}



// carefully test whether state is in Hilbert-space and find corresponding state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
// return value = corresponding index, or dimension of space, if not found
int FermionOnSphereWithSpin::CarefulFindStateIndex(unsigned long stateDescription, int highestBit)
{
   if (bitcount(stateDescription)!=this->NbrFermions)
    {
      return this->HilbertSpaceDimension;
    }
  if (highestBit<0)
    {
      highestBit = getHighestBit(stateDescription)-1;
    }
  if (highestBit >= 2*(LzMax+1))
    {
      return this->HilbertSpaceDimension;
    }
  int Index = this->FindStateIndex(stateDescription, highestBit);  
  if (this->StateDescription[Index] == stateDescription)
    return Index;
  else
    {
      cout << "Unexpected situation in CarefulFindStateIndex!"<<endl;
      for (int i=0; i<HilbertSpaceDimension; ++i)
	if (this->StateDescription[i] == stateDescription)
	  cout << "Element now found at i="<<i<<", "<<this->StateDescription[i]
	       <<"="<<stateDescription<<"!"<<endl;      
      return this->HilbertSpaceDimension;
    }

}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereWithSpin::FindStateIndex(unsigned long stateDescription, int lzmax)
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

ostream& FermionOnSphereWithSpin::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  for (int i = this->NbrLzValue-1; i >=0 ; --i)
    {
      Tmp = ((TmpState >> (i << 1)) & ((unsigned long) 0x3));
      if (Tmp == 0x1l)
	Str << "d ";
      else if (Tmp == 0x2l)
	Str << "u ";
      else if (Tmp == 0x3l)
	Str << "X ";
      else Str << "0 ";
    }
//  Str << " value = "<<TmpState;
//  Str << " position = " << this->FindStateIndex(TmpState, this->StateHighestBit[state]);
//   if (state !=  this->FindStateIndex(TmpState, this->StateHighestBit[state]))
//         Str << " error! ";
  return Str;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereWithSpin::PrintStateMonomial (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  Str << "[";
  int i = this->LzMax;
  while (((TmpState >> (i << 1)) & 0x3ul) == 0x0ul)
    --i;
  switch ((TmpState >> (i << 1)) & 0x3ul)
    {
    case 0x1ul:
      Str << i << "d";
      break;
    case 0x2ul:
      Str << i << "u";
      break;
    case 0x3ul:
      Str << i << "u," << i << "d";
      break;
    }
  --i;
  for (; i >=0; --i)
    switch ((TmpState >> (i << 1)) & 0x3ul)
      {
      case 0x1ul:
	Str << "," << i << "d";
	break;
      case 0x2ul:
	Str << "," << i << "u";
	break;
      case 0x3ul:
	Str << "," << i << "u," << i << "d";
	break;
      }
  Str << "]";
  return Str;
}

// print a given State using the monomial notation, separating spin up from spin down
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereWithSpin::PrintStateMonomialSeparatedSpin (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  Str << "[";
  int i = this->LzMax;
  while (((TmpState >> (i << 1)) & 0x1ul) == 0x0ul)
    --i;
  if (((TmpState >> (i << 1)) & 0x1ul) == 0x1ul)
    Str << i << "d";
  --i;
  for (; i >=0; --i)
    if (((TmpState >> (i << 1)) & 0x1ul) == 0x1ul)
      Str << "," << i << "d";
  i = this->LzMax;
  Str << ",";
  while (((TmpState >> (i << 1)) & 0x2ul) == 0x0ul)
    --i;
  if (((TmpState >> (i << 1)) & 0x2ul) == 0x2ul)
    Str << i << "u";
  --i;
  for (; i >=0; --i)
    if (((TmpState >> (i << 1)) & 0x2ul) == 0x2ul)
      Str << "," << i << "u";
  Str << "]";
  return Str;
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// currentLzMax = momentum maximum value for fermions that are still to be placed
// totalLz = momentum total value
// totalSz = spin total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int FermionOnSphereWithSpin::OldGenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSz)
{
  //  codage des etats sur deux bits, -lzMax up down on the lsb's
  
  /*----------------DECLARES---------------*/
  int Is_Lz, Is_Spin;
  unsigned long i, coeff, testLzMax;
  int k, position; 
  int CheckLz;
  int counter;
  int DimOrbit = lzMax+1;
  int currentLargestBit=2*lzMax+1;
  /*-------------INITS---------------------*/
  
  CheckLz = ((totalLz+nbrFermions*lzMax)/2);  //  CheckLz =totalLz+N*S

  i = biggestOne(nbrFermions,2*DimOrbit);


  testLzMax=0x1ul << currentLargestBit;
  counter = 0;        // on exit: dim of subspace
  
  while (i)
    {
      Is_Lz=0;
      Is_Spin=0;
      
      for(k=0;k<DimOrbit;k++)  // k indice va de 0 a 2S
	{  
          position = 2*k;                         // meaning 2*k
          coeff = ((i&(3ul<<position))>>position);
          
          switch(coeff)
	    {
            case 3:
	      {      // neutral to spin!
		Is_Lz += position;
	      }
	      break;
	      
            case 2:
	      { Is_Spin +=1;
	         Is_Lz +=k;}
	      break;
	      
            case 1:
	      { Is_Spin -=1;
	         Is_Lz +=k;}
	      break;
	      
            case 0:
	      // neutral to Spin and Lz
	      break;
	      
            default:
	      printf("severe error in fermion states");
	      break;
	      
	    }
          
          
	}
      
      
      if((Is_Lz == CheckLz) && (Is_Spin == totalSz) ) // project onto fixed spin and Lz
	{
	  this->StateDescription[counter]=i;
	  this->StateHighestBit[counter]=currentLargestBit;
	  counter++;
	}	
      
      i=lastone(i);
      // test if lzMax lowered in next word:
      if (!(i&testLzMax)) 
	{
	  --currentLargestBit;
	  testLzMax=1ul << currentLargestBit;
	}
    }
  return counter;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpin::GenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSpin, long pos)
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
	this->StateDescription[pos] = 0x1ul << ((totalLz << 1) + totalSpin);
	return (pos + 1l);
      }
    else
      return pos;

  if ((lzMax == 0)  && (totalLz != 0))
    return pos;


  long TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), totalSpin - 1,  pos);
  unsigned long Mask = 0x3ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1,  pos);
  Mask = 0x2ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin,  pos);
  Mask = 0x1ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  return this->GenerateStates(nbrFermions, lzMax - 1, totalLz, totalSpin, pos);
};


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereWithSpin::GenerateLookUpTable(unsigned long memory)
{
  // get every highest bit poisition
  unsigned long TmpPosition = this->StateDescription[0];
#ifdef __64_BITS__
  int CurrentHighestBit = 63;
#else
  int CurrentHighestBit = 31;
#endif
  while ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
    --CurrentHighestBit;  

  this->StateHighestBit[0] = CurrentHighestBit;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      TmpPosition = this->StateDescription[i];
      while ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
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
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#endif
}

// compute sign
//
// signs = 
// return value = sign value (+1.0 or -1.0)

double FermionOnSphereWithSpin::ComputeSign(unsigned long signs)
{
  unsigned result=0;
  while(signs) {
    result++;
    signs &= signs-1;
  }
  if (result & 1u) return -1.0;
  else return 1.0;
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = twce the total spin value
// return value = Hilbert space dimension

int FermionOnSphereWithSpin::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin)
{
  //  codage des etats sur deux bits, -lzMax up down on the lsb's
  
  /*----------------DECLARES---------------*/
  
  int Is_Lz, Is_Spin;
  unsigned long i, coeff;
  int k, position; 
  int CheckLz;
  int counter;
  int DimOrbit = lzMax+1;
  /*-------------INITS---------------------*/
  
  CheckLz = ((totalLz+nbrFermions*lzMax)/2);  //  CheckLz =totalLz+N*S

  i = biggestOne(nbrFermions,2*DimOrbit);
  
  counter = 0;        // on exit: dim of subspace
  
  while (i)
    {
      Is_Lz=0;
      Is_Spin=0;
      
      for(k=0;k<DimOrbit;k++)  // k indice va de 0 a 2S
	{  
          position = 2*k;                         // meaning 2*k
          coeff = ((i&(3ul<<position))>>position);
          
          switch(coeff)
	    {
            case 3:
	      {      // neutral to spin!
		Is_Lz += position;
	      }
	      break;
	      
            case 2:
	      { Is_Spin +=1;
	         Is_Lz +=k;}
	      break;
	      
            case 1:
	      { Is_Spin -=1;
	         Is_Lz +=k;}
	      break;
	      
            case 0:
	      // neutral to Spin and Lz
	      break;
	      
            default:
	      printf("severe error in fermion states");
	      break;
	      
	    }
          
          
	}
      
      
      if((Is_Lz == CheckLz) && (Is_Spin == totalSpin) ) // project onto fixed spin and Lz
	counter++;
	
      
      i=lastone(i);
    }
  return counter;
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// return value = Hilbert space dimension

long FermionOnSphereWithSpin::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin)
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

  unsigned long Tmp = 0l;  
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

Complex FermionOnSphereWithSpin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
							    int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
#ifdef __USE_LAPACK_HERE__
  ComplexLapackDeterminant SlaterUp(this->NbrFermionsUp);
  ComplexLapackDeterminant SlaterDown(this->NbrFermionsDown);
#else
  ComplexMatrix SlaterUp(this->NbrFermionsUp, this->NbrFermionsUp);
  ComplexMatrix SlaterDown(this->NbrFermionsDown, this->NbrFermionsDown);
#endif
  ComplexMatrix Functions(this->LzMax + 1, this->NbrFermions);
  RealVector TmpCoordinates(2);
  int* IndicesUp = new int [this->NbrFermionsUp];
  int* IndicesDown = new int [this->NbrFermionsDown];
  int PosUp, PosDown;
  int Lz;
  
  // calculate Basis functions for the given set of coordinates:
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
  // calculate prefactor 1/sqrt(N!);
  double Factor = 1.0;
  for (int i = 2; i <= this->NbrFermions; ++i)
    Factor *= (double) i;
  Factor = 1.0 / sqrt(Factor);
  // get down to the business: adding up terms \sum c_\alpha <r|\alpha>
  unsigned long TmpStateDescription;
  int LastComponent = firstComponent + nbrComponent;
  for (int k = firstComponent; k < LastComponent; ++k)
    {
      PosUp = 0;
      PosDown = 0;
      Lz = 0;
      TmpStateDescription = this->StateDescription[k];
      while ((PosUp < this->NbrFermionsUp)||(PosDown < this->NbrFermionsDown))
	{
	  if ((TmpStateDescription & 0x1l) != 0x0l)
	    {
	      IndicesDown[PosDown] = Lz;
	      ++PosDown;
	    }
	  if ((TmpStateDescription & 0x2l) != 0x0l)
	    {
	      IndicesUp[PosUp] = Lz;
	      ++PosUp;
	    }
	  ++Lz;
	  TmpStateDescription >>= 2;
	}
      for (int i = 0; i < this->NbrFermionsUp; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < this->NbrFermionsUp; ++j)
	    {
#ifdef __USE_LAPACK_HERE__
	      SlaterUp.SetMatrixElement(i,j,TmpColum2.Re(IndicesUp[j]), TmpColum2.Im(IndicesUp[j]));
#else
	      SlaterUp[i].Re(j) = TmpColum2.Re(IndicesUp[j]);
	      SlaterUp[i].Im(j) = TmpColum2.Im(IndicesUp[j]);
#endif
	    }
	}
      for (int i = 0; i < this->NbrFermionsDown; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i+this->NbrFermionsUp];	  
	  for (int j = 0; j < this->NbrFermionsDown; ++j)
	    {
#ifdef __USE_LAPACK_HERE__	      
	      SlaterDown.SetMatrixElement(i,j,TmpColum2.Re(IndicesDown[j]), TmpColum2.Im(IndicesDown[j]));
#else
	      SlaterDown[i].Re(j) = TmpColum2.Re(IndicesDown[j]);
	      SlaterDown[i].Im(j) = TmpColum2.Im(IndicesDown[j]);
#endif
	    }
	}
      Complex SlaterDetUp = SlaterUp.Determinant();
      Complex SlaterDetDown = SlaterDown.Determinant();
      Value += SlaterDetUp * SlaterDetDown * (state[k] * Factor) * this->GetStateSign(k, IndicesDown);
    }
  delete[] IndicesUp;
  delete[] IndicesDown;
  return Value;
}

// compute the sign for permuting all electrons with spin down to the right of those with spin up
// index = index of the state
// return value = sign value (+1.0 or -1.0)
// strategy: match up spin down particles in pairs, then move them out "for free"
// make use of me knowing the location of spin up particles already
double FermionOnSphereWithSpin::GetStateSign(int index, int* IndicesDown)
{ 
  unsigned long State = this->StateDescription[index];
  unsigned long mask, signs;
  int pos1, pos2;
  double result = 1.0;  
  for (int pair=0; pair<NbrFermionsDown/2; ++pair)
    {
      // create mask that highlights fermions between IndicesDown[2*pair] and IndicesDown[2*pair+1]
      pos2 = 2*IndicesDown[2*pair+1];
      pos1 = 2*IndicesDown[2*pair];
      mask = ((0x1l<<pos2) -1); // mask with all bits set at positions right of pos2
      mask ^=  ((0x1l<<(pos1+1)) -1); // ^ mask with all bits set at positions right of pos1+1
      signs = mask & State;
      result *= ComputeSign(signs);
    }
  if (NbrFermionsDown&1) // odd number of fermions...
    {
      pos1 = 2*IndicesDown[NbrFermionsDown-1];
      mask = ((0x1l<<pos1) -1);
      signs = mask & State;
      result *= ComputeSign(signs);
    }
  return result;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereWithSpin::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
  
// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

RealVector FermionOnSphereWithSpin::ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace)
{
  RealVector FinalState(this->HilbertSpaceDimension, true);
  for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
    {
      unsigned long TmpUpState = upStateSpace.StateDescription[j];
// #ifdef  __64_BITS__
//       TmpUpState |= TmpUpState << 32;
//       TmpUpState &= 0xffff0000fffful;
//       TmpUpState |= TmpUpState << 16;
//       TmpUpState &= 0xff00ff00ff00fful;
// #endif
//       TmpUpState |= TmpUpState << 16;
//       TmpUpState &= 0xffff0000fffful;
      
      int TmpPos = upStateSpace.LzMax;
      while (TmpPos > 0)
	{
	  unsigned long Tmp = TmpUpState & (0x1ul << TmpPos);
	  TmpUpState |= Tmp << TmpPos;
	  TmpUpState ^= Tmp;
	  --TmpPos;
	}
      TmpUpState <<= 1;
      double TmpComponent = upState[j];
      int Max = 63;
      while ((TmpUpState & (0x1ul << Max)) == 0x0ul)
	--Max;
      int Min = 0;
      while ((TmpUpState & (0x1ul << Min)) == 0x0ul)
	++Min;
      unsigned long TmpUpStateMask = (0x1ul << Max) - 1;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	if ((this->StateDescription[i] & TmpUpState) == TmpUpState)
	  {	    
	    unsigned long TmpUpState3 = this->StateDescription[i] & TmpUpStateMask;
	    unsigned long TmpUpState2 = TmpUpState3;
#ifdef  __64_BITS__
	    TmpUpState3 &= 0x5555555555555555ul;
	    TmpUpState2 &= 0xaaaaaaaaaaaaaaaaul;
#else
	    TmpUpState3 &= 0x55555555ul;
	    TmpUpState2 &= 0xaaaaaaaaul;
#endif	    
	    unsigned long Sign = 0x0;
	    int Pos = this->LzMax << 1;
	    while ((Pos > 0) && ((TmpUpState3 & (0x1ul << Pos)) == 0x0ul))
	      Pos -= 2;
	    while (Pos > 0)
	      {
		unsigned long TmpUpState4 = TmpUpState2 & ((0x1ul << Pos) - 1ul);
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
		while ((Pos > 0) && ((TmpUpState3 & (0x1ul << Pos)) == 0x0ul))
		  Pos -= 2;
	      }
	    if ((Sign & 0x1ul) == 0x0ul)
	      FinalState[i] = TmpComponent;
	    else
	      FinalState[i] = -TmpComponent;
	  }
    }

  for (int j = 0; j < downStateSpace.HilbertSpaceDimension; ++j)
    {
      unsigned long TmpDownState = downStateSpace.StateDescription[j];
      int TmpPos = downStateSpace.LzMax;
      while (TmpPos > 0)
	{
	  unsigned long Tmp = TmpDownState & (0x1ul << TmpPos);
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

RealVector FermionOnSphereWithSpin::ForgeU1FromSU2(RealVector& state, FermionOnSphere& u1Space)
{
  RealVector FinalState(u1Space.GetHilbertSpaceDimension(), true);
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      unsigned long TmpState = this->StateDescription[j];
      int TmpPos = this->LzMax << 1;
      while (TmpPos >= 0)
	{
	  if ((((TmpState >> TmpPos) & (TmpState >> (TmpPos + 1))) & 0x1ul) != 0x0ul)
	    TmpPos = 1;
	  TmpPos -= 2;
	}
      if (TmpPos != -1)
	{ 
	  TmpPos = 0;
	  unsigned long TmpState2 = 0x0ul; 
	  while (TmpPos <= this->LzMax)
	    {
	      TmpState2 |= ((TmpState & 0x1ul) | ((TmpState & 0x2ul) >> 1)) << TmpPos;
	      TmpState >>= 2;
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

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// return value = converted state

RealVector& FermionOnSphereWithSpin::ConvertToUnnormalizedMonomial(RealVector& state, long reference)
{
  int* TmpMonomialReference = new int [this->NbrFermions];
  int* TmpMonomial = new int [this->NbrFermions];
  double Factor = 1.0 / state[reference];
  state[reference] = 1.0;
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  unsigned long TmpState = this->StateDescription[reference];
  int Index = 0;
  for (int j = this->LzMax; j >= 0; --j)
    {
      switch ((TmpState >> (2 * j)) & 3ul)
	{
	case 0x1ul:
	  TmpMonomialReference[Index++] = j;
	  break;
	case 0x2ul:
	  TmpMonomialReference[Index++] = j;
	  break;
	case 0x3ul:
	  {
	    TmpMonomialReference[Index++] = j;
	    TmpMonomialReference[Index++] = j;
	  }
	  break;
	}
    }
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      Index = 0;
      TmpState = this->StateDescription[i];
      for (int j = this->LzMax; j >= 0; --j)
	{
	  switch ((TmpState >> (2 * j)) & 3ul)
	    {
	    case 0x1ul:
	      TmpMonomial[Index++] = j;
	      break;
	    case 0x2ul:
	      TmpMonomial[Index++] = j;
	      break;
	    case 0x3ul:
	      {
		TmpMonomial[Index++] = j;
		TmpMonomial[Index++] = j;
	      }
	      break;
	    }
	}
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	{
	  while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrFermions)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrFermions)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomialReference[Index2]];
	  ++Index2;
	}
      state[i] *= Coefficient;
    }
  delete[] TmpMonomialReference;
  delete[] TmpMonomial;
  state[reference] = 1.0;
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// return value = converted state

RealVector& FermionOnSphereWithSpin::ConvertFromUnnormalizedMonomial(RealVector& state, long reference)
{
  int* TmpMonomialReference = new int [this->NbrFermions];
  int* TmpMonomial = new int [this->NbrFermions];
  double Factor = 1.0 / state[reference];
  state[reference] = 1.0;
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      SqrtCoefficients[k] =  1.0 / InvSqrtCoefficients[k];
    }
  unsigned long TmpState = this->StateDescription[reference];
  int Index = 0;
  for (int j = this->LzMax; j >= 0; --j)
    {
      switch ((TmpState >> (2 * j)) & 3ul)
	{
	case 0x1ul:
	  TmpMonomialReference[Index++] = j;
	  break;
	case 0x2ul:
	  TmpMonomialReference[Index++] = j;
	  break;
	case 0x3ul:
	  {
	    TmpMonomialReference[Index++] = j;
	    TmpMonomialReference[Index++] = j;
	  }
	  break;
	}
    }
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      Index = 0;
      TmpState = this->StateDescription[i];
      for (int j = this->LzMax; j >= 0; --j)
	{
	  switch ((TmpState >> (2 * j)) & 3ul)
	    {
	    case 0x1ul:
	      TmpMonomial[Index++] = j;
	      break;
	    case 0x2ul:
	      TmpMonomial[Index++] = j;
	      break;
	    case 0x3ul:
	      {
		TmpMonomial[Index++] = j;
		TmpMonomial[Index++] = j;
	      }
	      break;
	    }
	}
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	{
	  while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrFermions)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrFermions)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomialReference[Index2]];
	  ++Index2;
	}
      state[i] *= Coefficient;
    }
  delete[] TmpMonomialReference;
  delete[] TmpMonomial;
  state /= state.Norm();
  return state;
}

// Evaluate the Density Matrix of the spin up fermions in a sector with a fixed lzUp 
//
// lzUp = twice total momentum of up fermions.
// groundstate = reference on the total system groundstate
// return value = density matrix of the subsystem of spins up fermions.

RealSymmetricMatrix FermionOnSphereWithSpin::EvaluatePartialDensityMatrixSpinSeparation (int lzUp, RealVector & groundstate)
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
  FermionOnSphereWithSpin TmpHilbertSpaceDown (nbrFermionDown,Complementarylz,LzMax,-nbrFermionDown);
  FermionOnSphereWithSpin TmpHilbertSpaceUp (NbrFermionsUp,lzUp,LzMax,NbrFermionsUp);
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
	      unsigned long TmpStateDescriptionRow=((TmpHilbertSpaceUp.StateDescription[i])|(TmpHilbertSpaceDown.StateDescription[k]));
	      unsigned long TmpStateDescriptionColomn=((TmpHilbertSpaceUp.StateDescription[j])|(TmpHilbertSpaceDown.StateDescription[k]));
	      
	      
	      int sign;
	      if(((perm[i][k]+perm[j][k])%2==1))
		{
		  sign=-1;
		}
	      else
		{
		  sign=1;
		}
	      int l=CarefulFindStateIndex(TmpStateDescriptionRow,-1);
	      
	      int r=CarefulFindStateIndex(TmpStateDescriptionColomn,-1);
	      
	      Coefficient+=groundstate[l]*groundstate[r]*sign;
	      
	    }
	  TmpDensityMatrix.SetMatrixElement(i,j,Coefficient);
	}
    }
  return TmpDensityMatrix;
}





