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
#include "GeneralTools/ArrayTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/StringTools.h"

#include <math.h>
#include <cstdlib>
#include <fstream>
#include <map>
#include <bitset>
#include <algorithm>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::map;
using std::pair;
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
  if(this->NbrFermions > 0)
    {
      this->LargeHilbertSpaceDimension = (int) this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
											  (this->TotalSpin + this->NbrFermions) >> 1);
      if (this->LargeHilbertSpaceDimension >= (1l << 30))
	this->HilbertSpaceDimension = 0;
      else
	this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
    }
  else
    this->HilbertSpaceDimension = 1l;

  this->Flag.Initialize();
  this->TargetSpace = this;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];  
      
      if (this->NbrFermions > 0)
	{
	  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
							     (this->TotalSpin + this->NbrFermions) >> 1, 0l);
	  this->GenerateLookUpTable(memory);
	}
      else
	{
	  this->StateDescription[0] = 0x0ul; 
	}
    }
  else
    {
      this->StateDescription = 0;
      this->StateHighestBit = 0;
      this->LookUpTableMemorySize = 0;
    }
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
  if (this->NbrFermions > 0)
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

FermionOnSphereWithSpin::~FermionOnSphereWithSpin ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      if (this->NbrFermions > 0)
	{
	  delete[] this->LookUpTableShift;
	  for (int i = 0; i < (2 * this->NbrLzValue); ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;
	  delete[] this->SignLookUpTableMask;
	  delete[] this->SignLookUpTable;
	}
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
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
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

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnSphereWithSpin::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
  this->TargetSpace = (FermionOnSphereWithSpin*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int FermionOnSphereWithSpin::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphereWithSpin::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);
  WriteLittleEndian(File, this->NbrFermions);
  WriteLittleEndian(File, this->LzMax);
  WriteLittleEndian(File, this->TotalLz);
  WriteLittleEndian(File, this->TotalSpin);
  if (this->HilbertSpaceDimension != 0)
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateDescription[i]);
    }
  else
    {
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateDescription[i]);
    }
  File.close();
  return true;
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
  if ((state & (0x1ul << m)) != 0x0ul)
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
      coefficient *= this->SignLookUpTable[(state >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
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
      return this->TargetSpace->HilbertSpaceDimension;
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
  if (((State & (0x1l << m2))!= 0x0ul)|| ((State & (0x1l << m1)) != 0x0ul))
    {
      //cout << "Third exit" << endl;
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
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
  return this->TargetSpace->FindStateIndex(State, NewLargestBit);
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
      return this->TargetSpace->HilbertSpaceDimension;
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
  if (((State & (0x1l << m2))!= 0x0ul)|| ((State & (0x1l << m1)) != 0x0ul))
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
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
  
  return this->TargetSpace->FindStateIndex(State, NewLargestBit);

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
      return this->TargetSpace->HilbertSpaceDimension;
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
  if (((State & (0x1l << m2))!= 0x0ul)|| ((State & (0x1l << m1)) != 0x0ul))
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
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
  
  return this->TargetSpace->FindStateIndex(State, NewLargestBit);
}


// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereWithSpin::AduAu (int index, int m)
{
  if ((this->StateDescription[index] & (0x2l << (m << 1))) != 0x0ul)
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
  if ((this->StateDescription[index] & (0x1l << (m << 1))) != 0x0ul)
    return 1.0;
  else
    return 0.0;
}



// apply a^+_m_u a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::AduAu (int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m = (m<<1) + 1;
  n = (n<<1) + 1;
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(0x1ul << n);
  if (NewLargestBit == n)
    while (((State >> NewLargestBit) == 0) && (NewLargestBit > 0))
      --NewLargestBit;

  if ((State & (0x1ul << m))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (0x1ul << m);
  return this->TargetSpace->FindStateIndex(State, NewLargestBit);
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::AddAd (int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m <<= 1;
  n <<= 1;
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(0x1ul << n);
  if (NewLargestBit == n)
    while (((State >> NewLargestBit) == 0) && (NewLargestBit > 0))
      --NewLargestBit;

  if ((State & (0x1ul << m))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (0x1ul << m);
  //cout << State << " " << NewLargestBit << endl;
  return this->TargetSpace->FindStateIndex(State, NewLargestBit);
}


// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereWithSpin::AduAd (int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m = (m << 1) + 1;
  n <<= 1;
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(0x1ul << n);
  if (NewLargestBit == n)
    while (((State >> NewLargestBit) == 0) && (NewLargestBit > 0))
      --NewLargestBit;

  if ((State & (0x1ul << m))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (0x1ul << m);
  return this->TargetSpace->FindStateIndex(State, NewLargestBit);
}



// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereWithSpin::AddAu (int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m <<= 1;
  n = (n << 1) + 1;  
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0))
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(0x1ul << n);
  if (NewLargestBit == n)
    while (((State >> NewLargestBit) == 0) && (NewLargestBit > 0))
      --NewLargestBit;

  if ((State & (0x1ul << m))!= 0x0ul)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (0x1ul << m);
  return this->TargetSpace->FindStateIndex(State, NewLargestBit);
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
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16)) & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16)) & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  if (this->ProdATemporaryState != 0x0ul)
    {
      while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
	--this->ProdALzMax;
    }
  else
    this->ProdALzMax = 0;
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
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16)) & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16)) & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  if (this->ProdATemporaryState != 0x0ul)
    {
      while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
	--this->ProdALzMax;
    }
  else
    this->ProdALzMax = 0;
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
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16)) & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16)) & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  if (this->ProdATemporaryState != 0x0ul)
    {
      while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
	--this->ProdALzMax;
    }
  else
    this->ProdALzMax = 0;
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
  if (((TmpState & (0x1ul << m1)) != 0x0ul) || ((TmpState & (0x1ul << m2)) != 0x0ul) || (m1 == m2))
    return this->TargetSpace->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16)) & this->SignLookUpTableMask[m2 + 16]];
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
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16)) & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
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
  if (((TmpState & (0x1ul << m1)) != 0x0ul) || ((TmpState & (0x1ul << m2)) != 0x0ul) || (m1 == m2))
    return this->TargetSpace->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16)) & this->SignLookUpTableMask[m2 + 16]];
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
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16)) & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
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
  if (((TmpState & (0x1ul << m1)) != 0x0ul) || ((TmpState & (0x1ul << m2)) != 0x0ul))
    return this->TargetSpace->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16)) & this->SignLookUpTableMask[m2 + 16]];
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
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16)) & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_u a^+_m2_u operator to a state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::AduAdu (int index, int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  ++m2;
  if (((TmpState & (0x1ul << m1)) != 0x0ul) || ((TmpState & (0x1ul << m2)) != 0x0ul) || (m1 == m2))
    return this->TargetSpace->HilbertSpaceDimension;
  int NewLzMax = this->StateHighestBit[index];
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16)) & this->SignLookUpTableMask[m2 + 16]];
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
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16)) & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}
   
// apply a^+_m1_d a^+_m2_d operator to a state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::AddAdd (int index, int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  m1 <<= 1;
  m2 <<= 1;
  if (((TmpState & (0x1ul << m1)) != 0x0ul) || ((TmpState & (0x1ul << m2)) != 0x0ul) || (m1 == m2))
    return this->TargetSpace->HilbertSpaceDimension;
  int NewLzMax = this->StateHighestBit[index];
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16)) & this->SignLookUpTableMask[m2 + 16]];
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
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16)) & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}
  
// apply a^+_m1_u a^+_m2_d operator to a state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::AduAdd (int index, int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  if (((TmpState & (0x1ul << m1)) != 0x0ul) || ((TmpState & (0x1ul << m2)) != 0x0l))
    return this->TargetSpace->HilbertSpaceDimension;
  int NewLzMax = this->StateHighestBit[index];
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16)) & this->SignLookUpTableMask[m2 + 16]];
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
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16)) & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}


// apply a_n1_u a_n2_d operator to a state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::AuAd (int index, int n1, int n2, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  n1 <<= 1;
  ++n1;
  n2 <<= 1;
  if (((TmpState & (0x1ul << n1)) == 0x0ul) || ((TmpState & (0x1ul << n2)) == 0x0ul))
    return this->TargetSpace->HilbertSpaceDimension;
  this->ProdALzMax = this->StateHighestBit[index];
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16)) & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(0x1ul << n2);
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16)) & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(0x1ul << n1);
  int NewLzMax = this->StateHighestBit[index];
  if (TmpState != 0x0ul)
    {
      while ((TmpState >> NewLzMax) == 0)
	--NewLzMax;
    }
  else
    NewLzMax = 0;
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}
  
// apply a_n_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator (spin up)
// return value =  multiplicative factor 

double FermionOnSphereWithSpin::Au (int index, int n)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n <<= 1;
  ++n;
  
  if ((this->ProdATemporaryState & (0x1ul << n)) == 0)
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n) & this->SignLookUpTableMask[n]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n);
  if (this->ProdATemporaryState != 0x0ul)
    {
      while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
	--this->ProdALzMax;
    }
  else
    this->ProdALzMax = 0;
  return Coefficient;
}

// apply a_n_d  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereWithSpin::Ad (int index, int n)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n <<= 1;
  
  if ((this->ProdATemporaryState & (0x1ul << n)) == 0)
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n) & this->SignLookUpTableMask[n]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n);
  if (this->ProdATemporaryState != 0x0ul)
    {
      while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
	--this->ProdALzMax;
    }
  else
    this->ProdALzMax = 0;
  return Coefficient;
}


// apply a^+_m_u operator to the state produced using Au method (without destroying it)
//
// m = first index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::Adu (int m, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m <<= 1;
  ++m;
  
  if ((TmpState & (0x1ul << m)) != 0x0ul)
    return this->TargetSpace->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  
  if (m > NewLzMax)
    NewLzMax = m;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  TmpState |= (0x1ul << m);
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m_d operator to the state produced using Au or Ad method (without destroying it)
//
// m = first index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::Add (int m, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m <<= 1;
  
  if ((TmpState & (0x1ul << m)) != 0x0ul)
    return this->TargetSpace->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  
  if (m > NewLzMax)
    NewLzMax = m;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  TmpState |= (0x1ul << m);
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}


// apply a^+_m_u  operator to a given state. 
//
// index = index of the state on which the operator has to be applied
// m = index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value =  index of the resulting state 

int FermionOnSphereWithSpin::Adu (int index, int m, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  m <<= 1; 
  ++m;
  if ((TmpState & (0x1ul << m)) != 0x0ul)
    return this->TargetSpace->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->StateHighestBit[index];  
  if (m > NewLzMax)
    NewLzMax = m;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  TmpState |= (0x1ul << m);
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m_d  operator to a given state. 
//
// index = index of the state on which the operator has to be applied
// m = index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value =  index of the resulting state 

int FermionOnSphereWithSpin::Add (int index, int m, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  m <<= 1;  
  if ((TmpState & (0x1ul << m)) != 0x0ul)
    return this->TargetSpace->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->StateHighestBit[index];  
  if (m > NewLzMax)
    NewLzMax = m;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  TmpState |= (0x1ul << m);
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply a_n_d operator to a state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// n = second index for annihilation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::Ad (int index, int n, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  n <<= 1;
  if ((TmpState & (0x1ul << n)) == 0x0ul)
    return this->TargetSpace->HilbertSpaceDimension;
  this->ProdALzMax = this->StateHighestBit[index];
  coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  TmpState &= ~(0x1ul << n);
  int NewLzMax = this->StateHighestBit[index];
  if (TmpState != 0x0ul)
    {
      while ((TmpState >> NewLzMax) == 0)
	--NewLzMax;
    }
  else
    NewLzMax = 0;
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}


// apply a_n_u operator to a state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// n = index for annihilation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::Au (int index, int n, double& coefficient)
{
  unsigned long TmpState = this->StateDescription[index];
  n <<= 1;
  n++;
  if ((TmpState & (0x1ul << n)) == 0x0ul)
    return this->TargetSpace->HilbertSpaceDimension;
  this->ProdALzMax = this->StateHighestBit[index];
  coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  TmpState &= ~(0x1ul << n);
  int NewLzMax = this->StateHighestBit[index];
  if (TmpState != 0x0ul)
    {
      while ((TmpState >> NewLzMax) == 0)
	--NewLzMax;
    }
  else
    NewLzMax = 0;

  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
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
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16)) & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      this->ProdATemporaryState &= ~(0x1l << Index);
    }
  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdALzMax = 0;
      return Coefficient;      
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
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16)) & this->SignLookUpTableMask[Index+ 16]];
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
      if ((TmpState & (0x1l << Index)) != 0x0ul)
	{
	  coefficient = 0.0;
	  return this->TargetSpace->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16)) & this->SignLookUpTableMask[Index + 16]];
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (0x1l << Index);
    }
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
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
      if ((TmpState & (0x1l << Index)) != 0x0ul)
	{
	  coefficient = 0.0;
	  return this->TargetSpace->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16)) & this->SignLookUpTableMask[Index + 16]];
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (0x1l << Index);
    }
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// get the variance of the state
//
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

// flip all spins of a given state
// 
// index = index of the state on which the operator has to be applied
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpin::SzToMinusSz (int index, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->StateDescription[index];
  unsigned long TmpState2 = 0x0ul;
  for (int j = 0; j <= this->LzMax; ++j)
    {
      switch (TmpState & 0x3ul)
	{
	case 0x3ul:	  
	  {
	    coefficient *= -1.0;
	    TmpState2 |= 0x3ul << (j << 1);
	  }
	  break;
	case 0x2ul:
	  {
	    TmpState2 |= 0x1ul << (j << 1);
	  }
	  break;
	case 0x1ul:
	  {
	    TmpState2 |= 0x2ul << (j << 1);
	  }
	  break;
	}
      TmpState >>= 2;
    }
  int TmpLzMax = 2 * this->LzMax + 1;
  while (((TmpState2 >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  return this->TargetSpace->CarefulFindStateIndex(TmpState2, TmpLzMax);
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
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    {
      return this->HilbertSpaceDimension;
    }
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
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}  

// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int FermionOnSphereWithSpin::FindStateIndex(char* stateDescription)
{
  char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != (this->LzMax + 1))
    return -1;
  unsigned long TmpState = 0x0ul;
  int TmpNbrParticles = 0;
  int TmpTotalLz = 0;
  int TmpTotalSz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      if (TmpDescription[i][0] == 'u')
	{
	  TmpState |= 0x2ul << (2 * i);
	  TmpTotalLz += i;
	  ++TmpTotalSz;
	  ++TmpNbrParticles;	  
	}
      else
	{
	  if (TmpDescription[i][0] == 'd')
	    {
	      TmpState |= 0x1ul << (2 * i);
	      TmpTotalLz += i;
	      --TmpTotalSz;
	      ++TmpNbrParticles;	  
	    }
	  else
	    {
	      if (TmpDescription[i][0] == 'X')
		{
		  TmpState |= 0x3ul << (2 * i);
		  TmpTotalLz += 2 * i;
		  TmpNbrParticles += 2;	  
		}
	      else
		{
		  if (TmpDescription[i][0] != '0')
		    {
		      return -1;
		    }
		}
	    }
	}
      delete[] TmpDescription[i];
    }
  delete[] TmpDescription;
  if ((TmpNbrParticles != this->NbrFermions) 
      || (TmpTotalLz != ((this->TotalLz + this->NbrFermions * this->LzMax) >> 1))
      || (TmpTotalSz != this->TotalSpin))
    return -1;
  int TmpLzMax = 2 * this->LzMax + 1;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
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
  for (int i = 0; i < this->NbrLzValue; ++i)
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
  return Str;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereWithSpin::PrintStateMonomial (ostream& Str, long state)
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

ostream& FermionOnSphereWithSpin::PrintStateMonomialSeparatedSpin (ostream& Str, long state)
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
    {
      if (lzMax >= totalLz)
	{
	  this->StateDescription[pos] = 0x1ul << ((totalLz << 1) + totalSpin);
	  return (pos + 1l);
	}
      else
	return pos;
    }
  if ((lzMax == 0) && (totalLz != 0))
    {
      return pos;
    }


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
  while (((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul) && (CurrentHighestBit > 0))
    --CurrentHighestBit;  

  if (this->StateHighestBit != 0)
    {
      this->StateHighestBit[0] = CurrentHighestBit;
      for (int i = 1; i < this->HilbertSpaceDimension; ++i)
	{
	  TmpPosition = this->StateDescription[i];
	  while (((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul) && (CurrentHighestBit > 0))
	    --CurrentHighestBit;  
	  this->StateHighestBit[i] = CurrentHighestBit;
	}
      CurrentHighestBit = this->StateHighestBit[0];
    }

  // evaluate look-up table size
  memory /= (sizeof(int*) * 2*this->NbrLzValue);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > (2 * this->NbrLzValue))
    this->MaximumLookUpShift = (2 * this->NbrLzValue);
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [2*this->NbrLzValue];
  this->LookUpTableShift = new int [2*this->NbrLzValue];
  for (int i = 0; i < 2*this->NbrLzValue; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentLargestBit = CurrentHighestBit;
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
      TmpPosition = this->StateDescription[i];
      while (((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul) && (CurrentHighestBit > 0))
	--CurrentHighestBit;  
      if (CurrentLargestBit != CurrentHighestBit)
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  CurrentLargestBit--;
	  while (CurrentLargestBit > CurrentHighestBit)
	    {
	      if (CurrentLargestBit < this->MaximumLookUpShift)
		this->LookUpTableShift[CurrentLargestBit] = 0;
	      else
		this->LookUpTableShift[CurrentLargestBit] = CurrentLargestBit + 1 - this->MaximumLookUpShift;
	      TmpLookUpTable = this->LookUpTable[CurrentLargestBit];
	      CurrentLookUpTableValue = this->LookUpTableMemorySize;
	      while (CurrentLookUpTableValue > 0x0ul)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      CurrentLargestBit--;
	    }
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
  this->GenerateSignLookUpTable();
}

// generate look-up table for sign calculation
// 
void FermionOnSphereWithSpin::GenerateSignLookUpTable()
{
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

long FermionOnSphereWithSpin::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin)
{
  //  codage des etats sur deux bits, -lzMax up down on the lsb's
  
  /*----------------DECLARES---------------*/
  
  int Is_Lz, Is_Spin;
  unsigned long i, coeff;
  int k, position; 
  int CheckLz;
  long counter;
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
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }

  if ((lzMax == 0) && (totalLz != 0))
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

RealVector FermionOnSphereWithSpin::ForgeSU2FromU1(RealVector& upState, ParticleOnSphere* upStateSpace, RealVector& downState, ParticleOnSphere* downStateSpace)
{
  return this->ForgeSU2FromU1(upState, *((FermionOnSphere*) upStateSpace), downState, *((FermionOnSphere*) downStateSpace));
}

// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

ComplexVector FermionOnSphereWithSpin::ForgeSU2FromU1(ComplexVector& upState, ParticleOnSphere* upStateSpace, ComplexVector& downState, ParticleOnSphere* downStateSpace)
{
  return this->ForgeSU2FromU1(upState, *((FermionOnSphere*) upStateSpace), downState, *((FermionOnSphere*) downStateSpace));
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
      unsigned long TmpUpStateMask = (0x1ul << (Max + 1)) - 0x1ul;
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
	{
	  if ((this->StateDescription[i] & TmpDownState) == TmpDownState)
	    {
	      FinalState[i] *= TmpComponent;
	    }
	}
    }

  return FinalState;
}

// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

ComplexVector FermionOnSphereWithSpin::ForgeSU2FromU1(ComplexVector& upState, FermionOnSphere& upStateSpace, ComplexVector& downState, FermionOnSphere& downStateSpace)
{
  ComplexVector FinalState(this->HilbertSpaceDimension, true);
  for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
    {
      unsigned long TmpUpState = upStateSpace.StateDescription[j];      
      int TmpPos = upStateSpace.LzMax;
      while (TmpPos > 0)
	{
	  unsigned long Tmp = TmpUpState & (0x1ul << TmpPos);
	  TmpUpState |= Tmp << TmpPos;
	  TmpUpState ^= Tmp;
	  --TmpPos;
	}
      TmpUpState <<= 1;
      Complex TmpComponent = upState[j];
      int Max = 63;
      while ((TmpUpState & (0x1ul << Max)) == 0x0ul)
	--Max;
      int Min = 0;
      while ((TmpUpState & (0x1ul << Min)) == 0x0ul)
	++Min;
      unsigned long TmpUpStateMask = (0x1ul << (Max + 1)) - 0x1ul;
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
      Complex TmpComponent = downState[j];
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  if ((this->StateDescription[i] & TmpDownState) == TmpDownState)
	    {
	      FinalState[i] *= TmpComponent;
	    }
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

// convert a given state from a generic basis to the current Sz subspace basis
//
// state = reference on the vector to convert
// basis = reference on the basis associated to state
// return value = converted vector

RealVector FermionOnSphereWithSpin::ConvertFromNbodyBasis(RealVector& state, FermionOnSphereWithSpin& basis)
{
  RealVector TmpVector (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    TmpVector[i] = state[basis.FindStateIndex(this->StateDescription[i], this->StateHighestBit[i])];
  return TmpVector;
}

// convert a given state from a generic basis to the current Sz subspace basis
//
// state = reference on the vector to convert
// basis = reference on the basis associated to state
// return value = converted vector

ComplexVector FermionOnSphereWithSpin::ConvertFromNbodyBasis(ComplexVector& state, FermionOnSphereWithSpin& basis)
{
  ComplexVector TmpVector (basis.GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    TmpVector[basis.FindStateIndex(this->StateDescription[i], this->StateHighestBit[i])] = state[i];
  
  return TmpVector;
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& FermionOnSphereWithSpin::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
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
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& FermionOnSphereWithSpin::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  int* TmpMonomialReference = new int [this->NbrFermions];
  int* TmpMonomial = new int [this->NbrFermions];
  double Factor = 1.0;
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
  
  
  unsigned int **perm;

  perm = new unsigned int*[TmpHilbertSpaceUp.HilbertSpaceDimension];
  for (int i=0; i<TmpHilbertSpaceUp.HilbertSpaceDimension; ++i)
    perm[i] = new unsigned int[TmpHilbertSpaceDown.HilbertSpaceDimension];
  
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

  for (int i=0; i<TmpHilbertSpaceUp.HilbertSpaceDimension; ++i)
    delete [] perm[i];
  delete [] perm;
  
  return TmpDensityMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix FermionOnSphereWithSpin::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, int szSector, RealVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrFermionSector == 0) && (szSector == 0))
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
      if ((lzSector == this->TotalLz) && (nbrFermionSector == this->NbrFermions) && (szSector == this->TotalSpin))
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

  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int ShiftedTotalLz = (this->TotalLz + this->NbrFermions * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrFermionSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector - (NbrFermionsComplementarySector * subsytemSize);
  int SzComplementarySector = this->TotalSpin - szSector;

  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  int NbrFermionsSectorUp = (nbrFermionSector + szSector) >> 1;
	  int NbrFermionsSectorDown = (nbrFermionSector - szSector) >> 1;
	  double TmpValue = 0.0;
 	  FermionOnSphereWithSpin TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize, SzComplementarySector);
	  unsigned long  TmpState2 = 0x0;
	  if (NbrFermionsSectorUp != 0)
	    TmpState2 |= 0x2ul;
	  if (NbrFermionsSectorDown != 0)
	    TmpState2 |= 0x1ul;
	  
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1)) | TmpState2;
	      int TmpLzMax = 2 * this->LzMax + 1;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
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
  if (nbrFermionSector == 0)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0;
 	  FermionOnSphereWithSpin TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize, SzComplementarySector);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1);
	      int TmpLzMax = 2 * this->LzMax + 1;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
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

  if (nbrFermionSector == 1)
    {
      double TmpValue = 0.0;
      FermionOnSphereWithSpin TmpHilbertSpace(NbrFermionsComplementarySector, (2 * ShiftedLzComplementarySector) - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize, SzComplementarySector);
      unsigned long TmpMask = 0x2l;
      if (szSector < 0)
	TmpMask = 0x1l;
      TmpMask <<= ShiftedLzSector << 1;
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1)) | TmpMask;
	  int TmpLzMax = 2 * this->LzMax + 1;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    TmpValue += groundState[TmpPos] * groundState[TmpPos];	
	}
      RealSymmetricMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  if (NbrFermionsComplementarySector == 0)
    {
      if ((ShiftedLzComplementarySector != 0) || (SzComplementarySector != 0))
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
      FermionOnSphereWithSpin TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1, szSector);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      int MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      double TmpValue;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpValue = groundState[MinIndex + i];
	  for (int j = i; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    TmpDensityMatrix.SetMatrixElement(i, j, TmpValue * groundState[MinIndex + j]);
	}
      return TmpDensityMatrix;
    }

  FermionOnSphereWithSpin TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1, szSector);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  if (TmpDestinationHilbertSpace.HilbertSpaceDimension == 0)
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  FermionOnSphereWithSpin TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize, SzComplementarySector);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = 2 * this->LzMax + 1;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
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

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealMatrix FermionOnSphereWithSpin::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int lzSector, int szSector, RealVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrFermionSector == 0) && (szSector == 0))
	{
	  RealMatrix TmpEntanglementMatrix(1, 1);
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
      if ((lzSector == this->TotalLz) && (nbrFermionSector == this->NbrFermions) && (szSector == this->TotalSpin))
	{
	  RealMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	      TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[i]);
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }

  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int ShiftedTotalLz = (this->TotalLz + this->NbrFermions * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrFermionSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector - (NbrFermionsComplementarySector * subsytemSize);
  int SzComplementarySector = this->TotalSpin - szSector;

  long TmpNbrNonZeroElements = 0;

  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  int NbrFermionsSectorUp = (nbrFermionSector + szSector) >> 1;
	  int NbrFermionsSectorDown = (nbrFermionSector - szSector) >> 1;
	  double TmpValue = 0.0;
 	  FermionOnSphereWithSpin TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize, SzComplementarySector);
	  unsigned long  TmpState2 = 0x0;
	  if (NbrFermionsSectorUp != 0)
	    TmpState2 |= 0x2ul;
	  if (NbrFermionsSectorDown != 0)
	    TmpState2 |= 0x1ul;
	  
          RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1)) | TmpState2;
	      int TmpLzMax = 2 * this->LzMax + 1;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
                {
                   TmpNbrNonZeroElements++;
		   TmpEntanglementMatrix.AddToMatrixElement(0, MinIndex, groundState[TmpPos]);	
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
  if (nbrFermionSector == 0)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0;
 	  FermionOnSphereWithSpin TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize, SzComplementarySector);
          RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1);
	      int TmpLzMax = 2 * this->LzMax + 1;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
                {
                   TmpNbrNonZeroElements++;
		   TmpEntanglementMatrix.AddToMatrixElement(0, MinIndex, groundState[TmpPos]);	
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

  if (nbrFermionSector == 1)
    {
      double TmpValue = 0.0;
      FermionOnSphereWithSpin TmpHilbertSpace(NbrFermionsComplementarySector, (2 * ShiftedLzComplementarySector) - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize, SzComplementarySector);
      RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
      unsigned long TmpMask = 0x2l;
      if (szSector < 0)
	TmpMask = 0x1l;
      TmpMask <<= ShiftedLzSector << 1;
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1)) | TmpMask;
	  int TmpLzMax = 2 * this->LzMax + 1;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
             {
                TmpNbrNonZeroElements++;
	        TmpEntanglementMatrix.AddToMatrixElement(0, MinIndex, groundState[TmpPos]);	
              }	
	}
       if (TmpNbrNonZeroElements == 0)
        {
          RealMatrix TmpEntanglementMatrix;
          return TmpEntanglementMatrix;
        }
       return TmpEntanglementMatrix;    
     }


  if (NbrFermionsComplementarySector == 0)
    {
      if ((ShiftedLzComplementarySector != 0) || (SzComplementarySector != 0))
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
      FermionOnSphereWithSpin TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1, szSector);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, 1, true);
      int MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	    TmpEntanglementMatrix.AddToMatrixElement(i, 0, groundState[MinIndex + i]);
	}
      return TmpEntanglementMatrix;
    }

  FermionOnSphereWithSpin TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1, szSector);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
 
  FermionOnSphereWithSpin TmpHilbertSpace(NbrFermionsComplementarySector, 2 * ShiftedLzComplementarySector - (NbrFermionsComplementarySector * (this->LzMax - subsytemSize)), this->LzMax - subsytemSize, SzComplementarySector);
 
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpHilbertSpace.HilbertSpaceDimension, true);
  
  TmpNbrNonZeroElements = 0;

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.StateDescription[MinIndex] << (subsytemSize << 1);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = 2 * this->LzMax + 1;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
              TmpNbrNonZeroElements++;
              TmpEntanglementMatrix.AddToMatrixElement(j, MinIndex, groundState[TmpPos]);
	    }
	}

     }

  if (TmpNbrNonZeroElements == 0)
   {
     RealMatrix TmpEntanglementMatrix;
     return TmpEntanglementMatrix;
   }
  return TmpEntanglementMatrix;    
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  FermionOnSphereWithSpin::EvaluatePartialDensityMatrixParticlePartition (int nbrFermionSector, int lzSector, int szSector, RealVector& groundState)
{  
  if (nbrFermionSector == 0)
    {
      if ((lzSector == 0) && (szSector == 0))
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

  if (nbrFermionSector == this->NbrFermions)
    {
      if ((lzSector == this->TotalLz) && (szSector == this->TotalSpin))
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

  int ComplementaryNbrFermionSector = this->NbrFermions - nbrFermionSector;
  BinomialCoefficients TmpBinomial (this->NbrFermions);
  double TmpInvBinomial = 1.0 / ( TmpBinomial(this->NbrFermions, nbrFermionSector) );
  /*
  if (nbrFermionSector == 1)
    {
      if (abs(szSector) == 1)
	{
	  double TmpValue = 0.0;
	  FermionOnSphereWithSpin TmpHilbertSpace(this->NbrFermions - 1, this->TotalLz - lzSector, this->LzMax, this->TotalSpin - szSector);
	  unsigned long ShiftedLzVSector = (lzSector + this->LzMax) >> 1;
	  unsigned long TmpMask = 0x1ul << ((ShiftedLzVSector << 1) + ((szSector + 1) >> 1));
	  unsigned long TmpMask2 = TmpMask - 1ul;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex];
	      if ((TmpState & TmpMask) == 0x0ul)
		{
		  TmpState |= TmpMask;
		  int TmpLzMax = (this->LzMax << 1) + 1;
		  while ((TmpState >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
		  if (TmpPos != this->HilbertSpaceDimension)
		    {
		      TmpValue += groundState[TmpPos] * groundState[TmpPos] * TmpInvBinomial;	
		    }
		}
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
  */
  
  FermionOnSphereWithSpin TmpDestinationHilbertSpace(nbrFermionSector, lzSector, this->LzMax, szSector);
  cout <<"nbrFermionSector = " << nbrFermionSector<< " lzSector = " << lzSector<< " LzMax =  "<< LzMax<< " szSector = "<< szSector<<endl;
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;
  FermionOnSphereWithSpin TmpHilbertSpace(ComplementaryNbrFermionSector, this->TotalLz - lzSector, this->LzMax, this->TotalSpin - szSector);
  cout <<"ComplementaryNbrFermionSector = " << ComplementaryNbrFermionSector << " lzSector = " << ( this->TotalLz - lzSector  )<< " LzMax =  "<< LzMax<< " szSector = "<< ( this->TotalSpin - szSector )<<endl;
  cout << "complementary Hilbert space dimension = " << TmpHilbertSpace.HilbertSpaceDimension << endl;
  TmpInvBinomial = sqrt(TmpInvBinomial);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex];
      cout << "TmpState = "<< TmpState <<endl;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace.StateDescription[j];
	  cout << "TmpState2 = "<< TmpState2 <<endl;
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
	      int TmpLzMax = (this->LzMax << 1) +1;
	      unsigned long TmpState3 = TmpState | TmpState2;
	      cout << "TmpState3 = "<< TmpState3 <<endl;
	      cout <<"TmpLzMax = "<<TmpLzMax<<endl;
	      while ((TmpState3 >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      cout <<"TmpLzMax = "<<TmpLzMax<<endl;
	      int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
 		  double Coefficient = TmpInvBinomial;
		  unsigned long Sign = 0x0ul;
		  int Pos2 = (TmpDestinationHilbertSpace.LzMax << 1) + 1;
		  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
		    {
		      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
			--Pos2;
		      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
		      TmpState3 ^= TmpState3 >> 32;
#endif	
		      TmpState3 ^= TmpState3 >> 16;
		      TmpState3 ^= TmpState3 >> 8;
		      TmpState3 ^= TmpState3 >> 4;
		      TmpState3 ^= TmpState3 >> 2;
		      TmpState3 ^= TmpState3 >> 1;
		      Sign ^= TmpState3;
		      TmpState2 &= ~(0x1ul << Pos2);
		      --Pos2;
		    }
 		  if ((Sign & 0x1ul) == 0x0ul)		  
 		    Coefficient *= 1.0;//TmpInvBinomial;
 		  else
 		    Coefficient *= -1.0;//TmpInvBinomial;
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  TmpStateCoefficient[Pos] = Coefficient;
		  ++Pos;
		}
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
		    TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  delete[] TmpStateCoefficient;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector and a gien Sz sector.
//
// nbrFermionSector = number of particles that belong to the subsytem
// lzSector = Lz sector in which the density matrix has to be evaluated
// szSector = Sz sector in which the density matrix has to be evaluated
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix FermionOnSphereWithSpin::EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int lzSector, int szSector, RealVector& groundState, bool removeBinomialCoefficient)
{
  if (nbrFermionSector == 0)
    {
      if ((lzSector == 0) && (szSector == 0))
        {
          FermionOnSphereWithSpin TmpHilbertSpace(this->NbrFermions, this->TotalLz - lzSector, this->LzMax, this->TotalSpin);
          RealMatrix TmpEntanglementMatrix(1, TmpHilbertSpace.HilbertSpaceDimension, true);
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              int TmpLzMax = (this->LzMax <<1) + 1;
              unsigned long TmpState = this->StateDescription[i];
              while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
             
              TmpEntanglementMatrix.SetMatrixElement(0, TmpHilbertSpace.FindStateIndex(TmpState, TmpLzMax), groundState[i]);
            }
          return TmpEntanglementMatrix;
        }
      else
        {
          RealMatrix TmpEntanglementMatrix;
          return TmpEntanglementMatrix;
        }
    }
 
 
  if (nbrFermionSector == this->NbrFermions)
    {
      if ((lzSector == this->TotalLz) && (szSector == this->TotalSpin))
        {
          FermionOnSphereWithSpin TmpDestinationHilbertSpace(nbrFermionSector, lzSector, this->LzMax, this->TotalSpin);
          RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, 1,true);
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              int TmpLzMax = (this->LzMax << 1) + 1;
              unsigned long TmpState = this->StateDescription[i];
              while ((TmpState >> TmpLzMax) == 0x0ul)
                --TmpLzMax;
	      
              TmpEntanglementMatrix.SetMatrixElement(TmpDestinationHilbertSpace.FindStateIndex(TmpState, TmpLzMax), 0, groundState[i]);
            }
          return TmpEntanglementMatrix;
        }
      else
        {
          RealMatrix TmpEntanglementMatrix;
          return TmpEntanglementMatrix;  
        }
    }
 
  int ComplementaryNbrFermionSector = this->NbrFermions - nbrFermionSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
 
  FermionOnSphereWithSpin TmpDestinationHilbertSpace(nbrFermionSector, lzSector, this->LzMax, szSector);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  long TmpNbrNonZeroElements = 0;

  double TmpInvBinomial = 1.0;
  if(removeBinomialCoefficient == false )
    {
      BinomialCoefficients TmpBinomial (this->NbrFermions);
      TmpInvBinomial = sqrt(1.0 / (TmpBinomial( this->NbrFermions , nbrFermionSector))); 
    }
 
 
  FermionOnSphereWithSpin TmpHilbertSpace(ComplementaryNbrFermionSector, this->TotalLz - lzSector, this->LzMax, ComplementarySzSector);
 
  FactorialCoefficient Factorial;
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpHilbertSpace.HilbertSpaceDimension, true);
 
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex];
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
        {
          unsigned long TmpState2 = TmpDestinationHilbertSpace.StateDescription[j];
          if ((TmpState & TmpState2) == 0x0ul)
            {
              int TmpLzMax = (this->LzMax << 1) + 1;
              unsigned long TmpState3 = TmpState | TmpState2;
              while ((TmpState3 >> TmpLzMax) == 0x0ul)
                --TmpLzMax;
              int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
              if (TmpPos != this->HilbertSpaceDimension)
                {
                  double Coefficient = TmpInvBinomial;
                  unsigned long Sign = 0x0ul;
                  int Pos2 = (TmpDestinationHilbertSpace.LzMax << 1) + 1;
                  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
                    {
                      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
                        --Pos2;
                      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
                      TmpState3 ^= TmpState3 >> 32;
#endif 
                      TmpState3 ^= TmpState3 >> 16;
                      TmpState3 ^= TmpState3 >> 8;
                      TmpState3 ^= TmpState3 >> 4;
                      TmpState3 ^= TmpState3 >> 2;
                      TmpState3 ^= TmpState3 >> 1;
                      Sign ^= TmpState3;
                      TmpState2 &= ~(0x1ul << Pos2);
                      --Pos2;
                    }
                  if ((Sign & 0x1ul) == 0x0ul)           
                    Coefficient *= 1.0;
                  else
                    Coefficient *= -1.0;
                 
                  TmpNbrNonZeroElements++;
                  TmpEntanglementMatrix.SetMatrixElement(j, MinIndex, Coefficient*groundState[TmpPos]);
                }
            }
        }
    }
 
  if (TmpNbrNonZeroElements > 0)
    return TmpEntanglementMatrix;
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnSphereWithSpin::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
										 RealVector& groundState, RealSymmetricMatrix* densityMatrix)
{
  FermionOnSphereWithSpin* TmpHilbertSpace = (FermionOnSphereWithSpin*) complementaryHilbertSpace;
  FermionOnSphereWithSpin* TmpDestinationHilbertSpace = (FermionOnSphereWithSpin*) destinationHilbertSpace;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;
  BinomialCoefficients TmpBinomial (this->NbrFermions);
  double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace->HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace->StateDescription[MinIndex];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace->StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
	      int TmpLzMax = (this->LzMax << 1) +1;
	      unsigned long TmpState3 = TmpState | TmpState2;
	      while ((TmpState3 >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
 		  double Coefficient = TmpInvBinomial;
		  unsigned long Sign = 0x0ul;
		  int Pos2 = (TmpDestinationHilbertSpace->LzMax << 1) + 1;
		  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
		    {
		      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
			--Pos2;
		      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
		      TmpState3 ^= TmpState3 >> 32;
#endif	
		      TmpState3 ^= TmpState3 >> 16;
		      TmpState3 ^= TmpState3 >> 8;
		      TmpState3 ^= TmpState3 >> 4;
		      TmpState3 ^= TmpState3 >> 2;
		      TmpState3 ^= TmpState3 >> 1;
		      Sign ^= TmpState3;
		      TmpState2 &= ~(0x1ul << Pos2);
		      --Pos2;
		    }
 		  if ((Sign & 0x1ul) == 0x0ul)		  
 		    Coefficient *= 1.0;
 		  else
 		    Coefficient *= -1.0;
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  TmpStateCoefficient[Pos] = Coefficient;
		  ++Pos;
		}
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
		  }
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnSphereWithSpin::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
										 ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
  FermionOnSphereWithSpin* TmpHilbertSpace = (FermionOnSphereWithSpin*) complementaryHilbertSpace;
  FermionOnSphereWithSpin* TmpDestinationHilbertSpace = (FermionOnSphereWithSpin*) destinationHilbertSpace;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  Complex* TmpStateCoefficient = new Complex [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;
  BinomialCoefficients TmpBinomial (this->NbrFermions);
  double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace->HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace->StateDescription[MinIndex];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace->StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
	      int TmpLzMax = (this->LzMax << 1) +1;
	      unsigned long TmpState3 = TmpState | TmpState2;
	      while ((TmpState3 >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
 		  double Coefficient = TmpInvBinomial;
		  unsigned long Sign = 0x0ul;
		  int Pos2 = (TmpDestinationHilbertSpace->LzMax << 1) + 1;
		  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
		    {
		      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
			--Pos2;
		      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
		      TmpState3 ^= TmpState3 >> 32;
#endif	
		      TmpState3 ^= TmpState3 >> 16;
		      TmpState3 ^= TmpState3 >> 8;
		      TmpState3 ^= TmpState3 >> 4;
		      TmpState3 ^= TmpState3 >> 2;
		      TmpState3 ^= TmpState3 >> 1;
		      Sign ^= TmpState3;
		      TmpState2 &= ~(0x1ul << Pos2);
		      --Pos2;
		    }
 		  if ((Sign & 0x1ul) == 0x0ul)		  
 		    Coefficient *= 1.0;
 		  else
 		    Coefficient *= -1.0;
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  TmpStateCoefficient[Pos] = Coefficient;
		  ++Pos;
		}
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
		  }
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
}


// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition. The entanglement matrix is only evaluated in a given Lz sector.
// and computed from precalculated particle entanglement matrix
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// thetaTop =  inclination angle defining one edge of the cut in degrees
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& FermionOnSphereWithSpin::EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, int szSector, double thetaTop, double thetaBottom, double phiRange, RealMatrix& entanglementMatrix)
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
  

  int ComplementaryNbrFermionsSector = this->NbrFermions - nbrFermionSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  FermionOnSphereWithSpin TmpDestinationHilbertSpace(nbrFermionSector, lzSector, this->LzMax, szSector);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  int* TmpMonomialUp = new int [ComplementaryNbrFermionsSector];
  int* TmpMonomialDown = new int [ComplementaryNbrFermionsSector];
  int* TmpDestinationMonomialUp = new int [this->NbrFermions];
  int* TmpDestinationMonomialDown = new int [this->NbrFermions];
  
  FermionOnSphereWithSpin TmpHilbertSpace(ComplementaryNbrFermionsSector, this->TotalLz - lzSector, this->LzMax, ComplementarySzSector);
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.StateDescription[i],TmpDestinationMonomialUp,TmpDestinationMonomialDown);
      double Tmp = 0.0;
      Tmp = 0.5 * nbrFermionSector * log(phiRange);
      for( int j = 0; j < TmpDestinationHilbertSpace.NbrFermionsUp; j++)
	{
	  Tmp += log( IncompleteBetaThetaBottom[TmpDestinationMonomialUp[j]] - IncompleteBetaThetaTop[TmpDestinationMonomialUp[j]]);
	}
      for( int j = 0; j < TmpDestinationHilbertSpace.NbrFermionsDown; j++)
	{
	  Tmp += log( IncompleteBetaThetaBottom[TmpDestinationMonomialDown[j]] - IncompleteBetaThetaTop[TmpDestinationMonomialDown[j]]);
	}
      Tmp = exp(0.5 * Tmp);
      for (int j = 0; j < TmpHilbertSpace.HilbertSpaceDimension; ++j)          
	entanglementMatrix(i, j) *= Tmp;      
    }
  
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex], TmpMonomialUp, TmpMonomialDown);
      double FormFactor = 0.0;
      for (int i = 0; i < TmpHilbertSpace.NbrFermionsUp; i++)
	FormFactor += log(1.0 - IncompleteBetaThetaBottom[TmpMonomialUp[i]] + IncompleteBetaThetaTop[TmpMonomialUp[i]] + (1.0 - phiRange) * (IncompleteBetaThetaBottom[TmpMonomialUp[i]] - IncompleteBetaThetaTop[TmpMonomialUp[i]]) );
      
      for (int i = 0; i <  TmpHilbertSpace.NbrFermionsDown; i++)
	FormFactor += log(1.0 - IncompleteBetaThetaBottom[TmpMonomialDown[i]] + IncompleteBetaThetaTop[TmpMonomialDown[i]] + (1.0 - phiRange) * (IncompleteBetaThetaBottom[TmpMonomialDown[i]] - IncompleteBetaThetaTop[TmpMonomialDown[i]]) );

      FormFactor = exp(0.5 * FormFactor);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	entanglementMatrix(j, MinIndex) *= FormFactor; 
    }
  
  delete[] TmpMonomialUp;
  delete[] TmpMonomialDown;  
  delete[] TmpDestinationMonomialDown;
  delete[] TmpDestinationMonomialUp;

  return entanglementMatrix;
}


// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization

double FermionOnSphereWithSpin::JackSqrNormalization (RealVector& outputVector, long minIndex, long nbrComponents)
{
  double SqrNorm = 0.0;
  int* TmpMonomialUp = new int [this->NbrFermionsUp];
  int* TmpMonomialDown = new int [this->NbrFermionsDown];
  FactorialCoefficient Factorial;
  int HalfLzMax = this->LzMax >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  for (long i = minIndex; i < MaxIndex; ++i)
    MaxIndex = this->LargeHilbertSpaceDimension;
  for (long i = minIndex; i < MaxIndex; ++i)
    if (outputVector[i] != 0.0)
      {
	Factorial.SetToOne();
	this->ConvertToMonomial(this->StateDescription[i], TmpMonomialUp, TmpMonomialDown);
	for (int k = 0; k < this->NbrFermionsUp; ++k)
	  {
	    if (HalfLzMax < TmpMonomialUp[k])
	      Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomialUp[k]);
	    else
	      if (HalfLzMax > TmpMonomialUp[k])
		Factorial.PartialFactorialDivide(TmpMonomialUp[k] + 1, HalfLzMax);
	  }
	for (int k = 0; k < this->NbrFermionsDown; ++k)
	  {
	    if (HalfLzMax < TmpMonomialDown[k])
	      Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomialDown[k]);
	    else
	      if (HalfLzMax > TmpMonomialDown[k])
		Factorial.PartialFactorialDivide(TmpMonomialDown[k] + 1, HalfLzMax);
	  }
	SqrNorm +=(outputVector[i] * outputVector[i]) * Factorial.GetNumericalValue();
	if ((i & 0x3fffl) == 0l)
	  {
	    cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	    cout.flush();
	  }
      }
  cout << endl;
  delete[] TmpMonomialUp;
  delete[] TmpMonomialDown;
  return SqrNorm;
}

// particle hole conjugate the spin down electrons, only (valid for N_phi=2N-1)
// source: input state vector
// target: output state vector
void FermionOnSphereWithSpin::ParticleHoleConjugateDownSpins(RealVector &source, RealVector &target)
{
  if (this->NbrLzValue!=2*NbrFermionsDown)
    {
      cout << "ParticleHoleConjugateDownSpins only defined for half filled down spins"<<endl;
      exit(1);
    }
  if (source.GetVectorDimension()!=this->HilbertSpaceDimension)
    {
      cout << "Input vector does not match Hilbert space dimension"<<endl;
      exit(1);
    }
  target.ResizeAndClean(this->HilbertSpaceDimension);
  unsigned long DownMask=0x0ul;
  unsigned long TmpState;
  int CurrentHighestBit, Index;
  for (int i=0; i<NbrLzValue; ++i)
    DownMask |= 0x1ul<<(i<<1);
  for (int i=0; i<this->HilbertSpaceDimension; ++i)
    {
      TmpState = this->StateDescription[i];
      TmpState ^= DownMask;
      CurrentHighestBit = (LzMax<<1) + 1;
      while ((TmpState & (0x1ul << CurrentHighestBit)) == 0x0ul)
	--CurrentHighestBit;  
      Index = this->FindStateIndex(TmpState,CurrentHighestBit);
      target[Index]=source[i];
    }
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

RealVector& FermionOnSphereWithSpin::FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
						 ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
						 bool symmetrizedFlag, double coefficient)
{
  FermionOnSphereWithSpin* LeftSpace = (FermionOnSphereWithSpin*) leftSpace;
  FermionOnSphereWithSpin* RightSpace = (FermionOnSphereWithSpin*) rightSpace;
  int StateShift = (RightSpace->LzMax + padding + 1) << 1;
  for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long TmpState1 = LeftSpace->StateDescription[i] << StateShift;
      double Coefficient = coefficient * leftVector[i];
      int TmpLzMax = 2 * this->LzMax + 1;
      while ((TmpState1 >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      if (symmetrizedFlag == false)
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState2 = RightSpace->StateDescription[j];
	      TmpState2 |= TmpState1;
	      double Coefficient2 = Coefficient;
	      Coefficient2 *= rightVector[j];	  
	      int TmpIndex = this->FindStateIndex(TmpState2, TmpLzMax);
	      outputVector[TmpIndex] = Coefficient2;
	    }
	}
      else
	{
// 	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
// 	    {
// 	      unsigned long TmpState2 = RightSpace->StateDescription[j];
// 	      TmpState2 |= TmpState1;
// 	      double Coefficient2 = Coefficient;
// 	      Coefficient2 *= rightVector[j];	  
// 	      int TmpIndex = this->FindStateIndex(TmpState2, TmpLzMax);
// 	      outputVector[TmpIndex] = Coefficient2;
// 	      unsigned long TmpState3 = this->GetSymmetricState(TmpState2);
// 	      if (TmpState3 != TmpState2)
// 		{
// 		  int TmpLzMax2 = this->LzMax;
// 		  while ((TmpState3 >> TmpLzMax2) == 0x0ul)
// 		    --TmpLzMax2;
// 		  TmpIndex = this->FindStateIndex(TmpState3, TmpLzMax2);
// 		  outputVector[TmpIndex] = Coefficient2;      
// 		}
// 	    }
	}
    }
  return outputVector;
}

// compute the product of a monomial and the halperin 110 state
//
// slater = array where the monomial representation of the slater determinant for half the number of particles is stored
// monomial = array where the monomial representation is stored
// sortingMap = map in which the generated states and their coefficient will be stored
// nbrPermutations = number of different permutations
// permutations1 = array where are stored the permutations of the spin up
// permutations2 = array where are stored the permutations of the spin down
// initialCoef = inital coefficient in front of the monomial

void FermionOnSphereWithSpin::MonomialsTimesPolarizedSlater(unsigned long * slater, unsigned long * monomial ,map<unsigned long , double> & sortingMap, unsigned long nbrPermutations , unsigned long * permutations1, unsigned long * permutations2,double initialCoef)
{
  unsigned long* State = new unsigned long[this->NbrFermions];
  pair <map <unsigned long, double>::iterator, bool> InsertionResult;
  
  int HalfNbrParticles = this->NbrFermions>>1;
  unsigned long * HalfMonomialsUp = new unsigned long[HalfNbrParticles];
  unsigned long * HalfMonomialsDown = new unsigned long[HalfNbrParticles];
  
  unsigned long TmpState = 0ul;
  unsigned long Mask = 0ul;
  unsigned long Sign = 0ul;
 
  double MonomialFact = initialCoef / (double) MultiplicitiesFactorial(monomial,this->NbrFermions);
  double TmpCoef;
  
  double CoefInitial;
  
  for (unsigned long IndexPermutations = 0; IndexPermutations < nbrPermutations ;IndexPermutations++)
    {
      unsigned long TmpPermUp = permutations1[IndexPermutations];
      unsigned long TmpPermDown = permutations2[IndexPermutations];
      
      for (int i = 0; i < HalfNbrParticles ; i++)
	{
	  HalfMonomialsUp[i] = monomial[(TmpPermUp >> (i * 5)) & 0x1ful];
	  HalfMonomialsDown[i] = monomial[(TmpPermDown >> (i * 5)) & 0x1ful];
	}
      
      CoefInitial =  ((double)MultiplicitiesFactorial(HalfMonomialsUp,HalfNbrParticles) * MultiplicitiesFactorial(HalfMonomialsDown,HalfNbrParticles)) * MonomialFact;
      
      do
	{
	  
	  for (int i = 0; i < HalfNbrParticles ; i++)
	    State[i] = HalfMonomialsUp[i] + slater[i];
	  
	  do
	    {
	      TmpCoef = CoefInitial;
	      for (int i = 0; i < HalfNbrParticles ; i++)
		State[i+HalfNbrParticles] = HalfMonomialsDown[i] + slater[i];
	      
	      TmpState = 0ul;
	      Sign = 0ul;
	      bool Bool = true;
	      
	      for (int i = 0; (i < HalfNbrParticles )&& (Bool); i++)
		{
		  Mask = (1ul << ((State[i]<<1) +1));
		  if((TmpState & Mask) != 0x0ul)
		    Bool = false;
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
	      
	      for (int i = HalfNbrParticles; (i < this->NbrFermions)&& (Bool); i++)
		{
		  Mask = (1ul << ((State[i]<<1)));
		  if((TmpState & Mask) != 0x0ul)
		    Bool = false;
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
	      
	      if (Bool)
		{
		  if ((Sign & 0x1ul) != 0ul)
		    TmpCoef *= -1l;
		  
		  InsertionResult = sortingMap.insert (pair <unsigned long, double> (TmpState, TmpCoef));
		  
		  if (InsertionResult.second == false)
		    {
		      InsertionResult.first->second += TmpCoef;
		    }	
		}
	    }
	  while (std::prev_permutation(HalfMonomialsDown, HalfMonomialsDown + HalfNbrParticles));
	}
      while (std::prev_permutation(HalfMonomialsUp, HalfMonomialsUp + HalfNbrParticles));
    }
  delete [] State;
}

// compute the projection of the product of a monomial in the two lowest LL and the halperin 110 state
//
// slater = array where the monomial representation of the slater determinant for half the number of particles is stored
// monomial = array where the monomial representation is stored
// sortingMap = map in which the generated states and their coefficient will be stored
// nbrPermutations = number of different permutations
// permutations1 = array where are stored the permutations of the spin up
// permutations2 = array where are stored the permutations of the spin down
// initialCoef = inital coefficient in front of the monomial
/*
void FermionOnSphereWithSpin::MonomialsTimesPolarizedSlaterProjection(unsigned long * slater, unsigned long * monomial, map<unsigned long , double> & sortingMap, unsigned long nbrPermutations , unsigned long * permutations1, unsigned long * permutations2, double initialCoef)
{
  unsigned long* State = new unsigned long[this->NbrFermions];
  pair <map <unsigned long, double>::iterator, bool> InsertionResult;
  
  int HalfNbrParticles = this->NbrFermions>>1;
  unsigned long * HalfMonomialsUp = new unsigned long[HalfNbrParticles];
  unsigned long * HalfMonomialsDown = new unsigned long[HalfNbrParticles];
  double CoefUp = 1.0;
  double CoefDown = 1.0;
  unsigned long TmpState = 0ul;
  unsigned long Mask = 0ul;
  unsigned long Sign = 0ul;
	
  long TmpLzMaxUp = this->LzMax - HalfNbrParticles + 3;
  long TmpFinalLzMaxUp = 2l + this->LzMax;
  double InverseFactor = 1.0 / (((double) TmpLzMaxUp) * ((double) TmpFinalLzMaxUp));
  double CoefInitial;
  double MonomialFact = initialCoef / (double) MultiplicitiesFactorial(monomial,this->NbrFermions);       
	
  for (unsigned long IndexPermutations = 0; IndexPermutations < nbrPermutations ; IndexPermutations++)
    {
      unsigned long TmpPermUp = permutations1[IndexPermutations];
      unsigned long TmpPermDown = permutations2[IndexPermutations];
		
      for (int i = 0; i < HalfNbrParticles ; i++)
	{
	  HalfMonomialsUp[i] = monomial[(TmpPermUp >> (i * 5)) & 0x1ful];
	  HalfMonomialsDown[i] = monomial[(TmpPermDown >> (i * 5)) & 0x1ful];
	}
      
      CoefInitial =  ((double)MultiplicitiesFactorial(HalfMonomialsUp,HalfNbrParticles) * MultiplicitiesFactorial(HalfMonomialsDown,HalfNbrParticles)) * MonomialFact;
      
      
      do
	{	
	  CoefUp = CoefInitial;
	  	  
	  for(int k = 0 ; k < HalfNbrParticles ; k++)
	    {
	      State[k] = (HalfMonomialsUp[k]>>1) + slater[k];
	      if ((HalfMonomialsUp[k] & 0x1ul) != 0ul) //not zero so have to project
		{
		  long Numerator = -((HalfMonomialsUp[k]>>1) * TmpFinalLzMaxUp) + (State[k] * TmpLzMaxUp);
		  if (Numerator == 0l)
		    { 
		      CoefUp = 0.0;		      
		      break;
		    }
		  else
		    CoefUp *= ((double) Numerator) * InverseFactor;
		}
	      State[k]--;
	    }
	  
	  if (CoefUp != 0.0)
	    {
	      do
		{
		  CoefDown = 1.0;		    
		  for(int k = 0 ; k < HalfNbrParticles ; k++)
		    {
		      State[k+HalfNbrParticles] = (HalfMonomialsDown[k]>>1) + slater[k];
		      if ((HalfMonomialsDown[k] & 0x1ul) != 0ul)
			{
			  long Numerator = -((HalfMonomialsDown[k]>>1) * TmpFinalLzMaxUp) + (State[HalfNbrParticles + k] * TmpLzMaxUp);
			  if (Numerator == 0l)
			    {
			      CoefDown = 0.0;
			      break;
			    }
			  else
			    CoefDown *= ((double) Numerator) * InverseFactor;
			}
		      State[HalfNbrParticles + k]--;
		    }
		  
		  if (CoefDown != 0.0)
		    {
		      
		      TmpState = 0ul;
		      Sign = 0ul;
		      bool Bool = true;
		      
		      for (int i = 0; i < HalfNbrParticles ; i++)
			{
			  Mask = (1ul << ((State[i]<<1) +1));
			  if((TmpState & Mask) != 0x0ul)
			    {
			      Bool = false;
			      break;
			    }
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
		      
		      if ( Bool ) 
			{
			  for (int i = HalfNbrParticles; i < this->NbrFermions ; i++)
			    {
			      Mask = (1ul << ((State[i]<<1)));
			      if((TmpState & Mask) != 0x0ul)
				{
				  Bool = false;
				  break;
				}
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
			}
		      		      
		      if (Bool)
			{
			  if ((Sign & 0x1ul) != 0ul)
			    CoefDown *= -1l;
			  
			  InsertionResult = sortingMap.insert (pair <unsigned long, double> (TmpState , CoefDown*CoefUp));
			  
			  if (InsertionResult.second == false)
			    {
			      InsertionResult.first->second += CoefDown*CoefUp;
			    }
			}
		    }
		}
	      while (std::prev_permutation(HalfMonomialsDown, HalfMonomialsDown + HalfNbrParticles));
	    }
	}
      while (std::prev_permutation(HalfMonomialsUp, HalfMonomialsUp + HalfNbrParticles));
    }
  delete [] State;
}*/


// compute the projection of the product of a monomial in the two lowest LL and the halperin 110 state
//
// slater = array where the monomial representation of the slater determinant for half the number of particles is stored
// monomial = array where the monomial representation is stored
// sortingMap = map in which the generated states and their coefficient will be stored
// nbrPermutations = number of different permutations
// permutations1 = array where are stored the permutations of the spin up
// permutations2 = array where are stored the permutations of the spin down
// initialCoef = inital coefficient in front of the monomial

void FermionOnSphereWithSpin::MonomialsTimesPolarizedSlaterProjection(unsigned long * slater, unsigned long * monomial, map<unsigned long , double> & sortingMap, unsigned long nbrPermutations , unsigned long * permutations1, unsigned long * permutations2, double initialCoef)
{
  unsigned long* State = new unsigned long[this->NbrFermions];
  pair <map <unsigned long, double>::iterator, bool> InsertionResult;
  
  int HalfNbrParticles = this->NbrFermions>>1;
  unsigned long * HalfMonomialsUp = new unsigned long[HalfNbrParticles];
  unsigned long * HalfMonomialsDown = new unsigned long[HalfNbrParticles];
  double CoefUp = 1.0;
  double CoefDown = 1.0;
  unsigned long TmpState = 0ul;
  unsigned long Mask = 0ul;
  unsigned long Sign = 0ul;
	
  long TmpLzMaxUp = this->LzMax - HalfNbrParticles + 3;
  long TmpFinalLzMaxUp = 2l + this->LzMax;
  double InverseFactor = 1.0 / (((double) TmpLzMaxUp) * ((double) TmpFinalLzMaxUp));
  double CoefInitial;
  double MonomialFact = initialCoef / (double) MultiplicitiesFactorial(monomial,this->NbrFermions);
    
  FactorialCoefficient *FactCoef = new FactorialCoefficient();
  FactCoef->FactorialMultiply(HalfNbrParticles);
  unsigned long** SlaterPermutations = new unsigned long*[FactCoef->GetIntegerValue()];  
  double *SlaterSigns = new double[FactCoef->GetIntegerValue()];
  delete FactCoef;
  unsigned long* TmpSlaterPermutation = new unsigned long[HalfNbrParticles];
  unsigned long* TmpSlaterPermutation2 = new unsigned long[HalfNbrParticles];
  unsigned long* TmpSlaterPermutation3 = new unsigned long[HalfNbrParticles];
  memcpy(TmpSlaterPermutation, slater, sizeof(unsigned long) * HalfNbrParticles);
  int NumPermutations;
  
  int NbrSlaterPermutations = 0;
  do 
    {
      SlaterPermutations[NbrSlaterPermutations] = new unsigned long[this->NbrFermions];      
      memcpy(SlaterPermutations[NbrSlaterPermutations], TmpSlaterPermutation, sizeof(unsigned long) * HalfNbrParticles);
      
      memcpy(TmpSlaterPermutation2, slater, sizeof(unsigned long) * HalfNbrParticles);
      memcpy(TmpSlaterPermutation3, SlaterPermutations[NbrSlaterPermutations], sizeof(unsigned long) * HalfNbrParticles);
      NumPermutations = 0;
      SortArrayDownOrdering(TmpSlaterPermutation3,  TmpSlaterPermutation2, HalfNbrParticles, NumPermutations);
      if ( (NumPermutations & 0x1) == 1 )
	{
	  SlaterSigns[NbrSlaterPermutations] = -1.0;
	}
      else
	{
	  SlaterSigns[NbrSlaterPermutations] = 1.0;
	}
      NbrSlaterPermutations++;
    }
  while (std::prev_permutation(TmpSlaterPermutation, TmpSlaterPermutation + HalfNbrParticles));
  delete [] TmpSlaterPermutation;
  delete [] TmpSlaterPermutation2;
  delete [] TmpSlaterPermutation3;
	
  for (unsigned long IndexPermutations = 0; IndexPermutations < nbrPermutations ; IndexPermutations++)
    {
      unsigned long TmpPermUp = permutations1[IndexPermutations];
      unsigned long TmpPermDown = permutations2[IndexPermutations];
		
      for (int i = 0; i < HalfNbrParticles ; i++)
	{
	  HalfMonomialsUp[i] = monomial[(TmpPermUp >> (i * 5)) & 0x1ful];
	  HalfMonomialsDown[i] = monomial[(TmpPermDown >> (i * 5)) & 0x1ful];
	}
      
      //CoefInitial =  ((double)MultiplicitiesFactorial(HalfMonomialsUp,HalfNbrParticles) * MultiplicitiesFactorial(HalfMonomialsDown,HalfNbrParticles)) * MonomialFact;
      
      CoefInitial =  MonomialFact;
      for ( int SlaterPermIndexUp = 0 ; SlaterPermIndexUp < NbrSlaterPermutations; SlaterPermIndexUp++ )
	{
	  
	  CoefUp = CoefInitial * SlaterSigns[SlaterPermIndexUp];
	  
	  for(int k = 0 ; k < HalfNbrParticles ; k++)
	    {
	      State[k] = (HalfMonomialsUp[k]>>1) + SlaterPermutations[SlaterPermIndexUp][k];
	      if ((HalfMonomialsUp[k] & 0x1ul) != 0ul) //not zero so have to project
		{
		  long Numerator = -((HalfMonomialsUp[k]>>1) * TmpFinalLzMaxUp) + (State[k] * TmpLzMaxUp);
		  if (Numerator == 0l)
		    { 
		      CoefUp = 0.0;		      
		      break;
		    }
		  else
		    CoefUp *= ((double) Numerator) * InverseFactor;
		}
	      State[k]--;
	    }
	  
	  if (CoefUp != 0.0)
	    {
	      for ( int SlaterPermIndex = 0 ; SlaterPermIndex < NbrSlaterPermutations; SlaterPermIndex++ )
		{
		  CoefDown = SlaterSigns[SlaterPermIndex];		    
		  for(int k = 0 ; k < HalfNbrParticles ; k++)
		    {
		      State[k+HalfNbrParticles] = (HalfMonomialsDown[k]>>1) + SlaterPermutations[SlaterPermIndex][k];
		      if ((HalfMonomialsDown[k] & 0x1ul) != 0ul)
			{
			  long Numerator = -((HalfMonomialsDown[k]>>1) * TmpFinalLzMaxUp) + (State[HalfNbrParticles + k] * TmpLzMaxUp);
			  if (Numerator == 0l)
			    {
			      CoefDown = 0.0;
			      break;
			    }
			  else
			    CoefDown *= ((double) Numerator) * InverseFactor;
			}
		      State[HalfNbrParticles + k]--;
		    }
		  
		  if (CoefDown != 0.0)
		    {
		      
		      TmpState = 0ul;
		      Sign = 0ul;
		      bool Bool = true;
		      
		      for (int i = 0; i < HalfNbrParticles ; i++)
			{
			  Mask = (1ul << ((State[i]<<1) +1));
			  if((TmpState & Mask) != 0x0ul)
			    {
			      Bool = false;
			      break;
			    }
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

                      if ( Bool ) 
			{
			  for (int i = HalfNbrParticles; i < this->NbrFermions ; i++)
			    {
			      Mask = (1ul << ((State[i]<<1)));
			      if((TmpState & Mask) != 0x0ul)
				{
				  Bool = false;
				  break;
				}
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
			}
		      		      
		      if (Bool)
			{
			  if ((Sign & 0x1ul) != 0ul)
			    CoefDown *= -1l;
			  
			  InsertionResult = sortingMap.insert (pair <unsigned long, double> (TmpState , CoefDown*CoefUp));
			  
			  if (InsertionResult.second == false)
			    {
			      InsertionResult.first->second += CoefDown*CoefUp;
			    }
			}
		    }
		}
	      //while (std::prev_permutation(HalfMonomialsDown, HalfMonomialsDown + HalfNbrParticles));
	    }
	}
      //while (std::prev_permutation(HalfMonomialsUp, HalfMonomialsUp + HalfNbrParticles));
    }
  for ( int SlaterPermIndex = 0 ; SlaterPermIndex < NbrSlaterPermutations; SlaterPermIndex++ )
    {
      delete [] SlaterPermutations[SlaterPermIndex];
    }
  delete [] SlaterPermutations;
  delete [] SlaterSigns;
  delete [] State;
}


// compute the projection of the product of a monomial in the two lowest LL and the halperin 110 state
//
// slaterPermutations = array of arrays where the monomial representation of the slater determinant for half the number of particles for each permutation are stored
// slaterSigns = array of the signs for each permutation of 
// nbrSlaterPermutations = number of permutations of the slater monomial in the array
// monomial = array where the monomial representation is stored
// sortingMap = map in which the generated states and their coefficient will be stored
// nbrPermutations = number of different permutations
// permutations1 = array where are stored the permutations of the spin up
// permutations2 = array where are stored the permutations of the spin down
// initialCoef = inital coefficient in front of the monomial

void FermionOnSphereWithSpin::MonomialsTimesPolarizedSlaterProjection(unsigned long ** slaterPermutations, double *slaterSigns, int nbrSlaterPermutations, unsigned long * monomial, map<unsigned long , double> & sortingMap, unsigned long nbrPermutations , unsigned long * permutations1, unsigned long * permutations2, double initialCoef)
{
  unsigned long* State = new unsigned long[this->NbrFermions];
  pair <map <unsigned long, double>::iterator, bool> InsertionResult;  
  int HalfNbrParticles = this->NbrFermions>>1;
  unsigned long * HalfMonomialsUp = new unsigned long[HalfNbrParticles];
  unsigned long * HalfMonomialsDown = new unsigned long[HalfNbrParticles];
  double CoefUp = 1.0;
  double CoefDown = 1.0;
  unsigned long TmpState = 0ul;
  unsigned long Mask = 0ul;
  unsigned long Mask2 = 0ul;
  unsigned long Sign = 0ul;
	
  long TmpLzMaxUp = this->LzMax - HalfNbrParticles + 3;
  long TmpFinalLzMaxUp = 2l + this->LzMax;
  double InverseFactor = 1.0 / (((double) TmpLzMaxUp) * ((double) TmpFinalLzMaxUp));
  double CoefInitial;
  double MonomialFact = initialCoef / (double) MultiplicitiesFactorial(monomial,this->NbrFermions);
    
  	
  for (unsigned long IndexPermutations = 0; IndexPermutations < nbrPermutations ; IndexPermutations++)
    {
      unsigned long TmpPermUp = permutations1[IndexPermutations];
      unsigned long TmpPermDown = permutations2[IndexPermutations];
		
      for (int i = 0; i < HalfNbrParticles ; i++)
	{
	  HalfMonomialsUp[i] = monomial[(TmpPermUp >> (i * 5)) & 0x1ful];
	  HalfMonomialsDown[i] = monomial[(TmpPermDown >> (i * 5)) & 0x1ful];
	}            
      
      CoefInitial =  MonomialFact;
      for ( int SlaterPermIndexUp = 0 ; SlaterPermIndexUp < nbrSlaterPermutations; SlaterPermIndexUp++ )
	{
	  
	  CoefUp = CoefInitial * slaterSigns[SlaterPermIndexUp];
	  
	  for(int k = 0 ; k < HalfNbrParticles ; k++)
	    {
	      State[k] = (HalfMonomialsUp[k]>>1) + slaterPermutations[SlaterPermIndexUp][k];
	      if ((HalfMonomialsUp[k] & 0x1ul) != 0ul) //not zero so have to project
		{
		  long Numerator = -((HalfMonomialsUp[k]>>1) * TmpFinalLzMaxUp) + (State[k] * TmpLzMaxUp);
		  if (Numerator == 0l)
		    { 
		      CoefUp = 0.0;		      
		      break;
		    }
		  else
		    CoefUp *= ((double) Numerator) * InverseFactor;
		}
	      State[k]--;
	    }
	  
	  if (CoefUp != 0.0)
	    {
	      for ( int SlaterPermIndex = 0 ; SlaterPermIndex < nbrSlaterPermutations; SlaterPermIndex++ )
		{
		  CoefDown = slaterSigns[SlaterPermIndex];		    
		  for(int k = 0 ; k < HalfNbrParticles ; k++)
		    {
		      State[k+HalfNbrParticles] = (HalfMonomialsDown[k]>>1) + slaterPermutations[SlaterPermIndex][k];
		      if ((HalfMonomialsDown[k] & 0x1ul) != 0ul)
			{
			  long Numerator = -((HalfMonomialsDown[k]>>1) * TmpFinalLzMaxUp) + (State[HalfNbrParticles + k] * TmpLzMaxUp);
			  if (Numerator == 0l)
			    {
			      CoefDown = 0.0;
			      break;
			    }
			  else
			    CoefDown *= ((double) Numerator) * InverseFactor;
			}
		      State[HalfNbrParticles + k]--;
		    }
		  
		  if (CoefDown != 0.0)
		    {
		      
		      TmpState = 0ul;
		      Sign = 0ul;
		      bool Bool = true;
		      
		     for (int i = 0; i < HalfNbrParticles ; i++)
			{
			  Mask = (1ul << ((State[i]<<1) +1));
			  if((TmpState & Mask) != 0x0ul)
			    {
			      Bool = false;
			      break;
			    }
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

                      if ( Bool ) 
			{
			  for (int i = HalfNbrParticles; i < this->NbrFermions ; i++)
			    {
			      Mask = (1ul << ((State[i]<<1)));
			      if((TmpState & Mask) != 0x0ul)
				{
				  Bool = false;
				  break;
				}
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
			}
                     
		      if (Bool)
			{
			  if ((Sign & 0x1ul) != 0ul)
			    CoefDown *= -1l;
			  
			  InsertionResult = sortingMap.insert (pair <unsigned long, double> (TmpState , CoefDown*CoefUp));
			  
			  if (InsertionResult.second == false)
			    {
			      InsertionResult.first->second += CoefDown*CoefUp;
			    }
			}
		    }
		}	      
	    }
	}     
    }  
  delete [] State;
}


// compute the product of a monomial and the halperin 110 state
//
// slater = array where the monomial representation of the slater determinant for half the number of particles is stored
// monomial = array where the monomial representation is stored
// sortingMap = map in which the generated states and their coefficient will be stored
// nbrPermutations = number of different permutations
// permutations1 = array where are stored the permutations of the spin up
// permutations2 = array where are stored the permutations of the spin down
// initialCoef = inital coefficient in front of the monomial

void FermionOnSphereWithSpin::MonomialsTimesPolarizedSlater(unsigned long * slaterUp, unsigned long * slaterDown, unsigned long * monomial ,map<unsigned long , double> & sortingMap, unsigned long nbrPermutations , unsigned long * permutations1, unsigned long * permutations2,double initialCoef)
{
  unsigned long* State = new unsigned long[this->NbrFermions];
  pair <map <unsigned long, double>::iterator, bool> InsertionResult;
  
  unsigned long * HalfMonomialsUp = new unsigned long[this->NbrFermionsUp];
  unsigned long * HalfMonomialsDown = new unsigned long[this->NbrFermionsDown];
  
  unsigned long TmpState = 0ul;
  unsigned long Mask = 0ul;
  unsigned long Sign = 0ul;
  long TmpLzMaxUp = this->LzMax - this->NbrFermionsUp + 3;
  long TmpFinalLzMaxUp = this->LzMax + 2;
  double InverseFactor = 1.0 / (((double) TmpLzMaxUp) * ((double) TmpFinalLzMaxUp));	
  
  double MonomialFact = initialCoef / (double) MultiplicitiesFactorial(monomial,this->NbrFermions);
  double TmpCoef;
  double CoefUp = 1.0;
  double CoefDown = 1.0;
	double CoefInitial ;
  
  for (unsigned long IndexPermutations = 0; IndexPermutations < nbrPermutations ;IndexPermutations++)
    {
      unsigned long TmpPermUp = permutations1[IndexPermutations];
      unsigned long TmpPermDown = permutations2[IndexPermutations];
      
      
      for (int i = 0; i < this->NbrFermionsUp ; i++)
	{
	  HalfMonomialsUp[i] = monomial[(TmpPermUp >> (i * 5)) & 0x1ful];
	}
      for (int i = 0; i < this->NbrFermionsDown ; i++)
	{
	  HalfMonomialsDown[i] = monomial[(TmpPermDown >> (i * 5)) & 0x1ful];
	}
      
      CoefInitial =  ((double)MultiplicitiesFactorial(HalfMonomialsUp,this->NbrFermionsUp) * MultiplicitiesFactorial(HalfMonomialsDown,this->NbrFermionsDown)) * MonomialFact;
      
      do
	{
	  CoefUp = CoefInitial;
	  
	  for(int k = 0 ; k < this->NbrFermionsUp ; k++)
	    {
	      State[k] = (HalfMonomialsUp[k]>>1) + slaterUp[k];
	      if ((HalfMonomialsUp[k] & 0x1ul) != 0ul) //not zero so have to project
		{
		  long Numerator = -((HalfMonomialsUp[k]>>1) * TmpFinalLzMaxUp) + (State[k] * TmpLzMaxUp);
		  if (Numerator == 0l)
		    { 
		      CoefUp = 0.0;		      
		      break;
		    }
		  else
		    CoefUp *= ((double) Numerator) * InverseFactor;
		}
	      State[k]--;
	    }
	    
	  if(CoefUp != 0.0)
	    {
	      do
		{
		  TmpCoef = CoefInitial;
		  CoefDown = 1.0;		    
		  
		  for(int k = 0 ; k < this->NbrFermionsDown ; k++)
		    {
		      State[k+this->NbrFermionsUp] = (HalfMonomialsDown[k]>>1) + slaterDown[k];
		      if ((HalfMonomialsDown[k] & 0x1ul) != 0ul)
			{
			  long Numerator = -((HalfMonomialsDown[k]>>1) * TmpFinalLzMaxUp) + (State[this->NbrFermionsUp + k] * TmpLzMaxUp);
			  if (Numerator == 0l)
			    {
			      CoefDown = 0.0;
			      break;
			    }
			  else
			    CoefDown *= ((double) Numerator) * InverseFactor;
			}
		      State[this->NbrFermionsUp + k]--;
		    }
		  
		  
		  if(CoefDown != 0.0)
		    {
		      TmpState = 0ul;
		      Sign = 0ul;
		      bool Bool = true;
		      
		      for (int i = 0; (i < this->NbrFermionsUp )&& (Bool); i++)
			{
			  Mask = (1ul << ((State[i]<<1) +1));
			  if((TmpState & Mask) != 0x0ul)
			    Bool = false;
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
	      
	      for (int i = this->NbrFermionsUp; (i < this->NbrFermions)&& (Bool); i++)
		{
		  Mask = (1ul << ((State[i]<<1)));
		  if((TmpState & Mask) != 0x0ul)
		    Bool = false;
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
	      
	      if (Bool)
		{
	  
		  if ((Sign & 0x1ul) != 0ul)
		    CoefDown *= -1l;
		  
		  InsertionResult = sortingMap.insert (pair <unsigned long, double> (TmpState, CoefDown*CoefUp));
		  
		  if (InsertionResult.second == false)
		    {
		      InsertionResult.first->second += CoefDown*CoefUp;
		    }	
		}
				}
	    }
			while (std::prev_permutation(HalfMonomialsDown, HalfMonomialsDown + this->NbrFermionsDown));
		}
	}
      while (std::prev_permutation(HalfMonomialsUp, HalfMonomialsUp + this->NbrFermionsUp));
    }
  delete [] State;
}

// apply a Gutzwiller projection (in the orbital space) to a given state
//
// state = reference on the state to project
// space = pointer to the Hilbert space where state is defined
// return value = Gutzwiller projected state

ComplexVector FermionOnSphereWithSpin::GutzwillerProjection(ComplexVector& state, ParticleOnSphere* space)
{
  FermionOnSphereWithSpin* TmpSpace = (FermionOnSphereWithSpin*) space;
  ComplexVector ProjectedState (this->LargeHilbertSpaceDimension, true);
  for (long i = 0l; i < TmpSpace->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = TmpSpace->StateDescription[i];
#ifdef  __64_BITS__
      if ((((TmpState & 0xaaaaaaaaaaaaaaaaul) >> 1) & (TmpState & 0x5555555555555555ul)) == 0x0ul)
#else
      if ((((TmpState & 0xaaaaaaaaul) >> 1) & (TmpState & 0x55555555ul)) == 0x0ul)
#endif	    
	{
	  int TmpIndex = this->FindStateIndex(TmpState, TmpSpace->StateHighestBit[i]);
	  if (TmpIndex < this->HilbertSpaceDimension)
	    {
	      ProjectedState[TmpIndex] = state[i];
	    }
	}
    }
  return ProjectedState;
}

// convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void FermionOnSphereWithSpin::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis,
						    long firstComponent, long nbrComponents)
{
  int* TmpMomentumIndices = new int [this->NbrFermions];
  int* TmpSU2Indices = new int [this->NbrFermions];
  int* TmpSU2Indices2 = new int [this->NbrFermions];
  targetState.ClearVector();
  long LastComponent = firstComponent + nbrComponents;
  if (nbrComponents == 0)
    LastComponent = this->LargeHilbertSpaceDimension;
  for (long i = firstComponent; i < LastComponent; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      unsigned long Tmp;
      int TmpLzMax = this->NbrLzValue << 1;
      int TmpIndex = 0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  Tmp = (TmpState >> (j << 1)) & 0x3ul;;
	  if ((Tmp & 0x2ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 1;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x1ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 0;
	      ++TmpIndex;
	    }	  
	}
      this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpMomentumIndices, TmpSU2Indices, TmpSU2Indices2, oneBodyBasis);
    }
  delete[] TmpMomentumIndices;
  delete[] TmpSU2Indices;
  delete[] TmpSU2Indices2;
}

// compute the transformation matrix from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// return value = transformation matrix

ComplexMatrix FermionOnSphereWithSpin::TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis)
{
  int* TmpMomentumIndices = new int [this->NbrFermions];
  int* TmpSU2Indices = new int [this->NbrFermions];
  int* TmpSU2Indices2 = new int [this->NbrFermions];
  ComplexMatrix TmpMatrix(this->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      unsigned long Tmp;
      int TmpLzMax = this->NbrLzValue << 1;
      int TmpIndex = 0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  Tmp = (TmpState >> (j << 1)) & 0x3ul;;
	  if ((Tmp & 0x2ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 1;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x1ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 0;
	      ++TmpIndex;
	    }	  
	}
      this->TransformOneBodyBasisRecursive(TmpMatrix[i], 1.0, 0, TmpMomentumIndices, TmpSU2Indices, TmpSU2Indices2, oneBodyBasis);
    }
  delete[] TmpMomentumIndices;
  delete[] TmpSU2Indices;
  delete[] TmpSU2Indices2;
  return TmpMatrix;
}

// recursive part of the convertion from a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// targetState = vector where the transformed state has to be stored
// coefficient = current coefficient to assign
// position = current particle consider in the n-body state
// momentumIndices = array that gives the momentum partition of the initial n-body state
// initialSU2Indices = array that gives the spin dressing the initial n-body state
// currentSU2Indices = array that gives the spin dressing the current transformed n-body state
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector

void FermionOnSphereWithSpin::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
							     int position, int* momentumIndices, int* initialSU2Indices,
							     int* currentSU2Indices, ComplexMatrix* oneBodyBasis) 
{
  if (position == this->NbrFermions)
    {
      unsigned long TmpState = 0x0ul;
      unsigned long TmpState2;
      unsigned long Mask = 0x0ul;
      unsigned long MaskSign = 0x0ul;
      for (int i = 0; i < this->NbrFermions; ++i)
	{
	  Mask = 0x1ul << ((momentumIndices[i] << 1) + currentSU2Indices[i]);
	  if ((TmpState & Mask) != 0x0ul)
	    return;
	  TmpState2 = TmpState & (Mask - 0x1ul);
#ifdef __64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
#endif
	  TmpState2 ^= (TmpState2 >> 16);
	  TmpState2 ^= (TmpState2 >> 8);
	  TmpState2 ^= (TmpState2 >> 4);
	  TmpState2 ^= (TmpState2 >> 2);
	  MaskSign ^= (TmpState2 ^ (TmpState2 >> 1)) & 0x1ul;
	  TmpState |= Mask;
	}
      int TmpLzMax = this->NbrLzValue << 1;
      while ((TmpState >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      int Index = this->FindStateIndex(TmpState, TmpLzMax);
      if (Index < this->HilbertSpaceDimension)
	{
	  if (MaskSign == 0ul)
	    {
	      targetState[Index] += coefficient;
	    }
	  else
	    {
	      targetState[Index] -= coefficient;
	    }
	}
      return;      
    }
  else
    {
      currentSU2Indices[position] = 0;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSU2Indices[position]][1]), position + 1, momentumIndices, initialSU2Indices, currentSU2Indices, oneBodyBasis);
      currentSU2Indices[position] = 1;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSU2Indices[position]][0]), position + 1, momentumIndices, initialSU2Indices, currentSU2Indices, oneBodyBasis);
    }
}
