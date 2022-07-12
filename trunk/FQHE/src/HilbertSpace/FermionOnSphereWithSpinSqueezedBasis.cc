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
//                        class author: Gunnar Möller                         //
//                        last modification : 20/12/2007                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinSqueezedBasis.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/Endian.h"

#include <math.h>
#include <bitset>
#include <cstdlib>

using std::cout;
using std::endl;
using std::ios;
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

FermionOnSphereWithSpinSqueezedBasis::FermionOnSphereWithSpinSqueezedBasis()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twce the total spin value
// referenceState = array that describes the reference state to start from (each entry has two bits corresponding to up 2l / down 1l)
// memory = amount of memory granted for precalculations
FermionOnSphereWithSpinSqueezedBasis::FermionOnSphereWithSpinSqueezedBasis (int nbrFermions, int &totalLz, int lzMax, int totalSpin,  int* referenceState, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalSpin = totalSpin;
  if ((NbrFermions&1) !=  (TotalSpin&1))
    {
      cout << "NbrFermions and TotalSpin need to have the same parity"<<endl;
      exit(1);
    }
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;  
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  // generate binary representation of reference state:
  this->ReferenceState = 0x0ul;
  int ReferenceStateHighestBit = 0;
  this->TotalLz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      this->ReferenceState |= ((unsigned long) (referenceState[i] & 3)) << 2*i;
      switch(referenceState[i]&3)
	{
	case 1:
	  ReferenceStateHighestBit = 2*i;
	  this->TotalLz += i;
	case 2:
	  ReferenceStateHighestBit = 2*i+1;
	  this->TotalLz += i;
	case 3:
	  ReferenceStateHighestBit = 2*i+1;
	  this->TotalLz += 2*i;
	}
    }
  this->TotalLz = ((this->TotalLz << 1) - (this->LzMax * this->NbrFermions)) >> 1;
  totalLz = this->TotalLz;
  // calculate full HilbertSpace dimension
  this->HilbertSpaceDimension = (int) this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
										 (this->TotalSpin + this->NbrFermions) >> 1);
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];
  // assign bit-valued Flags of whether to keep states:
#ifdef  __64_BITS__
  int ReducedHilbertSpaceDimension = (this->HilbertSpaceDimension >> 6) + 1;
#else
  int ReducedHilbertSpaceDimension = (this->HilbertSpaceDimension >> 5) + 1;
#endif
  this->KeepStateFlag = new unsigned long [ReducedHilbertSpaceDimension];
  for (int i = 0; i < ReducedHilbertSpaceDimension; ++i)
    this->KeepStateFlag[i] = 0x0l;
  // generate full Hilbert-space
  this->HilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, this->LzMax,
							(this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
							(this->TotalSpin + this->NbrFermions) >> 1, 0l);
  cout << "Full Hilbert-space dimension: " <<HilbertSpaceDimension<<endl;
  this->GenerateLookUpTable(memory);
  
  int MaxSweeps = (this->NbrFermionsUp * (this->NbrFermionsUp - 1))/2 // max number of squeezes within up-spins
    + (this->NbrFermionsDown * (this->NbrFermionsDown - 1))/2 // max number of squeezes within up-spins
    + (this->NbrFermionsUp * this->NbrFermionsDown)*2; // max number of generalized squeezes between different spins
  
  this->TmpGeneratedStates =  new unsigned long [MaxSweeps * 1000];
  this->TmpGeneratedStatesHighestBit = new int [MaxSweeps * 1000];
  long Memory = 0l;

  // mark reference state to be kept:
  int TmpIndex = this->FindStateIndex(this->ReferenceState, ReferenceStateHighestBit);
#ifdef  __64_BITS__
  this->KeepStateFlag[TmpIndex >> 6] = 0x1l << (TmpIndex & 0x3f);
#else
  this->KeepStateFlag[TmpIndex >> 5] = 0x1l << (TmpIndex & 0x1f);
#endif
  // recursively, determine descendents:
  this->GenerateDescendingStates(ReferenceStateHighestBit, this->ReferenceState, Memory);

  for (int i=0; i<HilbertSpaceDimension; ++i)
    {
      this->PrintState(cout, i);
      cout << " " << StateHighestBit[i];
      cout << endl;
    }

  int NewHilbertSpaceDimension = 0;
  unsigned long TmpKeepStateFlag;
  int TmpNbrOne[] = {  
  0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
  4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
  for (int i = 0; i < ReducedHilbertSpaceDimension; ++i)
    {
      TmpKeepStateFlag = this->KeepStateFlag[i];
      NewHilbertSpaceDimension += TmpNbrOne[TmpKeepStateFlag & 0xffl];
      NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 8) & 0xffl];
      NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 16) & 0xffl];
      NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 24) & 0xffl];
#ifdef  __64_BITS__
      NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 32) & 0xffl];
      NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 40) & 0xffl];
      NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 48) & 0xffl];
      NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 56) & 0xffl];      
#endif
    }

  delete[] this->SignLookUpTable;
  delete[] this->SignLookUpTableMask;
  delete[] this->LookUpTableShift;
  for (int i = 0; i < this->NbrLzValue; ++i)
    delete[] this->LookUpTable[i];
  delete[] this->LookUpTable;
  unsigned long* TmpStateDescription = new unsigned long [NewHilbertSpaceDimension];
  int* TmpStateHighestBit = new int [NewHilbertSpaceDimension];
  NewHilbertSpaceDimension = 0;
  int TotalIndex = 0;
#ifdef  __64_BITS__
  if ((this->HilbertSpaceDimension & 0x3f) != 0)
#else
  if ((this->HilbertSpaceDimension & 0x1f) != 0)
#endif
    --ReducedHilbertSpaceDimension;
  for (int i = 0; i < ReducedHilbertSpaceDimension; ++i)
    {
      TmpKeepStateFlag = this->KeepStateFlag[i];
#ifdef  __64_BITS__
      for (int j = 0; j < 64; ++j)
#else
      for (int j = 0; j < 32; ++j)
#endif
	{
	  if ((TmpKeepStateFlag >> j) & 0x1l)
	    {
	      TmpStateDescription[NewHilbertSpaceDimension] =  this->StateDescription[TotalIndex];
	      TmpStateHighestBit[NewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
	      ++NewHilbertSpaceDimension;
	    }
	  else
	    {	      
	      cout << "discarding ";
	      this->PrintState(cout, TotalIndex);
	      cout << endl;
	    }
	  ++TotalIndex;
	}
    }
#ifdef  __64_BITS__
  this->HilbertSpaceDimension &= 0x3f;
 #else
  this->HilbertSpaceDimension &= 0x1f;
 #endif
  if (this->HilbertSpaceDimension != 0)
    {
      TmpKeepStateFlag = this->KeepStateFlag[ReducedHilbertSpaceDimension];
      for (int j = 0; j < this->HilbertSpaceDimension; ++j)
	{
	  if ((TmpKeepStateFlag >> j) & 0x1l)
	    {
	      TmpStateDescription[NewHilbertSpaceDimension] =  this->StateDescription[TotalIndex];
	      TmpStateHighestBit[NewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
	      ++NewHilbertSpaceDimension;
	    }
	  else
	    {	      
	      cout << "discarding ";
	      this->PrintState(cout, TotalIndex);
	      cout << endl;
	    }
	  ++TotalIndex;
	}
    }
  
  delete[] this->StateDescription;
  delete[] this->StateHighestBit;
  delete[] this->KeepStateFlag;
  this->StateDescription = TmpStateDescription;
  this->StateHighestBit = TmpStateHighestBit;
  this->HilbertSpaceDimension = NewHilbertSpaceDimension;
  cout << "Reduced Hilbert-space dimension: " <<HilbertSpaceDimension<<endl;
  delete[] this->TmpGeneratedStates;
  delete[] this->TmpGeneratedStatesHighestBit;

  for (int i=0; i<HilbertSpaceDimension; ++i)
    {
      this->PrintState(cout, i);
      cout << " " << StateHighestBit[i];
      cout << endl;
    }

  this->GenerateLookUpTable(memory);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  

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

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinSqueezedBasis::FermionOnSphereWithSpinSqueezedBasis (char* fileName, unsigned long memory)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      this->HilbertSpaceDimension = 0;
      return;
    }
  ReadLittleEndian(File, this->HilbertSpaceDimension);
  ReadLittleEndian(File, this->NbrFermions);
  ReadLittleEndian(File, this->LzMax);
  ReadLittleEndian(File, this->TotalLz);
  ReadLittleEndian(File, this->TotalSpin);
  ReadLittleEndian(File, this->ReferenceState);
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    ReadLittleEndian(File, this->StateDescription[i]);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    ReadLittleEndian(File, this->StateHighestBit[i]);

  File.close();

  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;  
  this->Flag.Initialize();

  this->GenerateLookUpTable(memory);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  unsigned long UsedMemory = 0l;
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << this->MaximumSignLookUp) * sizeof(double);
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

FermionOnSphereWithSpinSqueezedBasis::FermionOnSphereWithSpinSqueezedBasis(const FermionOnSphereWithSpinSqueezedBasis& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->ReferenceState = fermions.ReferenceState;
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnSphereWithSpinSqueezedBasis::~FermionOnSphereWithSpinSqueezedBasis ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinSqueezedBasis& FermionOnSphereWithSpinSqueezedBasis::operator = (const FermionOnSphereWithSpinSqueezedBasis& fermions)
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
  this->ReferenceState = fermions.ReferenceState;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
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

AbstractHilbertSpace* FermionOnSphereWithSpinSqueezedBasis::Clone()
{
  return new FermionOnSphereWithSpinSqueezedBasis(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphereWithSpinSqueezedBasis::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->NbrFermions);
  WriteLittleEndian(File, this->LzMax);
  WriteLittleEndian(File, this->TotalLz);
  WriteLittleEndian(File, this->TotalSpin);
  WriteLittleEndian(File, this->ReferenceState);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    WriteLittleEndian(File, this->StateDescription[i]);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    WriteLittleEndian(File, this->StateHighestBit[i]);
  
  File.close();
  return true;
}


// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnSphereWithSpinSqueezedBasis::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnSphereWithSpinSqueezedBasis::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphereWithSpinSqueezedBasis::ExtractSubspace (AbstractQuantumNumber& q, 
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

//#include <bitset>
int FermionOnSphereWithSpinSqueezedBasis::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

int FermionOnSphereWithSpinSqueezedBasis::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

int FermionOnSphereWithSpinSqueezedBasis::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

double FermionOnSphereWithSpinSqueezedBasis::AduAu (int index, int m)
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

double FermionOnSphereWithSpinSqueezedBasis::AddAd (int index, int m)
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

double FermionOnSphereWithSpinSqueezedBasis::AuAu (int index, int n1, int n2)
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

double FermionOnSphereWithSpinSqueezedBasis::AdAd (int index, int n1, int n2)
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

double FermionOnSphereWithSpinSqueezedBasis::AuAd (int index, int n1, int n2)
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

int FermionOnSphereWithSpinSqueezedBasis::AduAdu (int m1, int m2, double& coefficient)
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

int FermionOnSphereWithSpinSqueezedBasis::AddAdd (int m1, int m2, double& coefficient)
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

int FermionOnSphereWithSpinSqueezedBasis::AduAdd (int m1, int m2, double& coefficient)
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

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereWithSpinSqueezedBasis::FindStateIndex(unsigned long stateDescription, int lzmax)
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

ostream& FermionOnSphereWithSpinSqueezedBasis::PrintState (ostream& Str, int state)
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
//   Str << " position = " << this->FindStateIndex(TmpState, this->StateHighestBit[state]);
//   if (state !=  this->FindStateIndex(TmpState, this->StateHighestBit[state]))
//         Str << " error! ";
  return Str;
}


// generate all descendents of the given referenceState by squeezing operations
// 
// highestBit = highest non-zero bit in reference state
// referenceState = state whose direct descendents should be calculated
// memory = index of memory slot to be used
void FermionOnSphereWithSpinSqueezedBasis::GenerateDescendingStates(int highestBit, unsigned long referenceState, long& memory)
{
  int MaxSweepsUpUp = (this->NbrFermionsUp * (this->NbrFermionsUp - 1)) >> 1;
  int MaxSweepsDownDown = (this->NbrFermionsDown * (this->NbrFermionsDown - 1)) >> 1;
  int MaxSweepsUpDown = this->NbrFermionsUp * this->NbrFermionsDown;
  int MaxSweeps = MaxSweepsUpUp + MaxSweepsDownDown + (MaxSweepsUpDown << 1);
  // assign temporary memory avoiding overwriting in recursion
  unsigned long* TmpGeneratedStates2 = this->TmpGeneratedStates + (MaxSweeps * memory);
  int* TmpHighestBit = this->TmpGeneratedStatesHighestBit  + (MaxSweeps  * memory);
  memory += 1;

  // up-up squeezing operations
  int TmpCurrentLzMax = 2;
  int TmpCurrentLzMax2;
  int TmpMaxLz, TmpMaxLzUp;
  // calculate highest Lz of spin up fermion less 1:
  if (highestBit&1)
    TmpMaxLz = (highestBit>>1) - 1;
  else
    {
      TmpMaxLz = highestBit-1;
      while ((TmpMaxLz>0) && ((referenceState & (1l<<TmpMaxLz)) == 0)) 
	TmpMaxLz-=2;
      TmpMaxLz = (highestBit>>1) - 1;
    }
  TmpMaxLzUp = TmpMaxLz;
  bitset<32> test;
  test = referenceState;
  // cout << "Reference-state =    " << test << " TmpMaxLz=" <<TmpMaxLz<<", TmpCurrentLzMax="<<TmpCurrentLzMax<<" highestBit="<<highestBit << endl;
  // cout << "up-up descendents"<<endl;
  int NbrEntries = 0;
  unsigned long TmpReferenceState;  
  while (TmpCurrentLzMax < TmpMaxLz)
    {
      while ((TmpCurrentLzMax < TmpMaxLz) && (((referenceState >> (TmpCurrentLzMax<<1)) & 0xal) != 0x8l))
	++TmpCurrentLzMax;      
      if (TmpCurrentLzMax < TmpMaxLz)
	{
	  TmpReferenceState = (referenceState & ~(0xal << (TmpCurrentLzMax<<1))) | (0x2l << (TmpCurrentLzMax<<1));
	  TmpCurrentLzMax2 = TmpCurrentLzMax - 2;
	  while (TmpCurrentLzMax2 >= 0)
	    {
	      while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> (TmpCurrentLzMax2<<1)) & 0xal) != 0x2l))
		--TmpCurrentLzMax2;
	      if (TmpCurrentLzMax2 >= 0)
		{
		  TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0xal << (TmpCurrentLzMax2<<1))) | (0x8l << (TmpCurrentLzMax2<<1));
		  
		  test = TmpGeneratedStates2[NbrEntries];
		  // cout << "Descending state " <<NbrEntries << " = " << test << endl;		  
		  TmpHighestBit[NbrEntries] = highestBit;
		  ++NbrEntries;
		  --TmpCurrentLzMax2;
		}	      
	    }
	  ++TmpCurrentLzMax;
	}
    }
  if (((referenceState >> (TmpCurrentLzMax<<1)) & 0xal) == 0x8l)
    {
      TmpReferenceState = (referenceState & ~(0xal << (TmpCurrentLzMax<<1))) | (0x2l << (TmpCurrentLzMax<<1));
      TmpCurrentLzMax2 = TmpCurrentLzMax - 2;
      int NewHighestBit;      
      if ( (TmpReferenceState & (1l << (highestBit-1))) != 0l )
	NewHighestBit = highestBit - 1;
      else NewHighestBit = highestBit - 2;      
      while (TmpCurrentLzMax2 >= 0)
	{
	  while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> (TmpCurrentLzMax2<<1)) & 0xal) != 0x2l))
	    --TmpCurrentLzMax2;
	  if (TmpCurrentLzMax2 >= 0)
	    {
	      TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0xal << (TmpCurrentLzMax2<<1))) | (0x8l << (TmpCurrentLzMax2<<1));	      
	      TmpHighestBit[NbrEntries] = NewHighestBit;
	      test = TmpGeneratedStates2[NbrEntries];
	      // cout << "Descending state " <<NbrEntries << " = " << test << " (highest bit " << TmpHighestBit[NbrEntries] << ")" << endl;		  
	      ++NbrEntries;
	      --TmpCurrentLzMax2;
	    }
	}      
    }

  // down-down squeezing operations
  TmpCurrentLzMax = 2;
  // calculate highest Lz of spin down fermion less 1:
  if (highestBit&1)
    {
      TmpMaxLz = highestBit-1;
      while ((TmpMaxLz>0) && ((referenceState & (1l<<TmpMaxLz)) == 0)) TmpMaxLz-=2;
      TmpMaxLz = (highestBit>>1) - 1;
    }
  else
    TmpMaxLz = (highestBit>>1) - 1;
  int TmpMaxLzDown = TmpMaxLz;
  // cout << "down down descendents" << endl;
  while (TmpCurrentLzMax < TmpMaxLz)
    {
      while ((TmpCurrentLzMax < TmpMaxLz) && (((referenceState >> (TmpCurrentLzMax<<1)) & 0x5l) != 0x4l))
	++TmpCurrentLzMax;
      if (TmpCurrentLzMax < TmpMaxLz)
	{
	  TmpReferenceState = (referenceState & ~(0x5l << (TmpCurrentLzMax<<1))) | (0x1l << (TmpCurrentLzMax<<1));
	  TmpCurrentLzMax2 = TmpCurrentLzMax - 2;
	  while (TmpCurrentLzMax2 >= 0)
	    {
	      while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> (TmpCurrentLzMax2<<1)) & 0x5l) != 0x1l))
		--TmpCurrentLzMax2;
	      if (TmpCurrentLzMax2 >= 0)
		{
		  TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0x5l << (TmpCurrentLzMax2<<1))) | (0x4l << (TmpCurrentLzMax2<<1));
		  
		  test = TmpGeneratedStates2[NbrEntries];
		  // cout << "Descending state " <<NbrEntries << " = " << test << endl;		  
		  TmpHighestBit[NbrEntries] = highestBit;
		  ++NbrEntries;
		  --TmpCurrentLzMax2;
		}	      
	    }
	  ++TmpCurrentLzMax;
	}
    }
  if (((referenceState >> (TmpCurrentLzMax<<1)) & 0x5l) == 0x4l)
    {
      TmpReferenceState = (referenceState & ~(0x5l << (TmpCurrentLzMax<<1))) | (0x1l << (TmpCurrentLzMax<<1));
      TmpCurrentLzMax2 = TmpCurrentLzMax - 2;
      int NewHighestBit;
      // cout << "final dd highestBit=" << highestBit << endl;
      if (highestBit&1)
	NewHighestBit = highestBit;
      else
	{
	  NewHighestBit = highestBit-1;
	  if ((TmpReferenceState & (1l << NewHighestBit)) == 0l ) --NewHighestBit;
	}
      while (TmpCurrentLzMax2 >= 0)
	{
	  while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> (TmpCurrentLzMax2<<1)) & 0x5l) != 0x1l))
	    --TmpCurrentLzMax2;
	  if (TmpCurrentLzMax2 >= 0)
	    {
	      TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0x5l << (TmpCurrentLzMax2<<1))) | (0x4l << (TmpCurrentLzMax2<<1));	      
	      TmpHighestBit[NbrEntries] = NewHighestBit;
	      test = TmpGeneratedStates2[NbrEntries];
	      // cout << "Descending state " <<NbrEntries << " = " << test << " (highest bit " << TmpHighestBit[NbrEntries] << ")" << endl;		  
	      ++NbrEntries;
	      --TmpCurrentLzMax2;
	    }
	}      
    }


  // up-down squeezing operations
  TmpCurrentLzMax = 2;
  TmpMaxLz = TmpMaxLzUp;
  int TmpMaxLz2 = TmpMaxLzDown+1;
  if (TmpMaxLz2>=this->LzMax)
    {
      TmpMaxLz2 = this->LzMax-1;
      while ((referenceState & (1l << (TmpMaxLz2<<1))) == 0)
	--TmpMaxLz2;
    }
  test = referenceState;
  // cout << "Reference-state =    " << test << " TmpMaxLz=" <<TmpMaxLz<<", TmpCurrentLzMax="<<TmpCurrentLzMax<<" highestBit="<<highestBit <<", TmpMaxLz2="<<TmpMaxLz2<< endl;
  // cout << "up-down descendents"<<endl;
  while (TmpCurrentLzMax < TmpMaxLz)
    {
      while ((TmpCurrentLzMax < TmpMaxLz) && (((referenceState >> (TmpCurrentLzMax<<1)) & 0xal) != 0x8l))
	++TmpCurrentLzMax;
      if (TmpCurrentLzMax < TmpMaxLz)
	{
	  TmpReferenceState = (referenceState & ~(0xal << (TmpCurrentLzMax<<1))) | (0x2l << (TmpCurrentLzMax<<1));
	  TmpCurrentLzMax2 = 0;
	  while (TmpCurrentLzMax2 <= TmpMaxLz2)
	    {
	      while ((TmpCurrentLzMax2 <= TmpMaxLz2) && (((referenceState >> (TmpCurrentLzMax2<<1)) & 0x5l) != 0x1l))
		++TmpCurrentLzMax2;
	      if (TmpCurrentLzMax2 <= TmpMaxLz2)
		{
		  TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0x5l << (TmpCurrentLzMax2<<1))) | (0x4l << (TmpCurrentLzMax2<<1));
		  
		  test = TmpGeneratedStates2[NbrEntries];
		  // cout << "Descending state " <<NbrEntries << " = " << test;
		  if (TmpCurrentLzMax2<(highestBit>>1))
		    TmpHighestBit[NbrEntries] = highestBit;
		  else
		    TmpHighestBit[NbrEntries] = ((TmpCurrentLzMax2+1)<<1);
		  // cout << ", highestBit = " << TmpHighestBit[NbrEntries] << " at TmpCurrentLzMax2="<<TmpCurrentLzMax2<< endl;
		  ++NbrEntries;
		  ++TmpCurrentLzMax2;
		}
	    }
	  ++TmpCurrentLzMax;
	}
    }
  if (((referenceState >> (TmpCurrentLzMax<<1)) & 0xal) == 0x8l)
    {
      TmpReferenceState = (referenceState & ~(0xal << (TmpCurrentLzMax<<1))) | (0x2l << (TmpCurrentLzMax<<1));
      TmpCurrentLzMax2 = 0;
      int NewHighestBit;      
      if ( (TmpReferenceState & (1l << (highestBit-1))) != 0l )
	NewHighestBit = highestBit - 1;
      else NewHighestBit = highestBit - 2;
      while (TmpCurrentLzMax2 <= TmpMaxLz2)
	{
	  while ((TmpCurrentLzMax2 <= TmpMaxLz2) && (((referenceState >> (TmpCurrentLzMax2<<1)) & 0xal) != 0x2l))
	    ++TmpCurrentLzMax2;
	  if (TmpCurrentLzMax2 <= TmpMaxLz2)
	    {
	      TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0x5l << (TmpCurrentLzMax2<<1))) | (0x4l << (TmpCurrentLzMax2<<1));
	      if ((TmpCurrentLzMax2<<1)<highestBit-2)
		TmpHighestBit[NbrEntries] = NewHighestBit;
	      else TmpHighestBit[NbrEntries] = ((TmpCurrentLzMax2+1)<<1);
	      test = TmpGeneratedStates2[NbrEntries];
	      // cout << "Descending state " <<NbrEntries << " = " << test << " (highest bit " << TmpHighestBit[NbrEntries] << " at TmpCurrentLzMax="<<TmpCurrentLzMax<< " at TmpCurrentLzMax2="<<TmpCurrentLzMax2<<")" << endl;		  
	      ++NbrEntries;
	      ++TmpCurrentLzMax2;
	    }
	}      
    }

  // down-up squeezing operations
  TmpCurrentLzMax = 2;
  TmpMaxLz = TmpMaxLzDown;
  TmpMaxLz2 = TmpMaxLzUp+1;
  if (TmpMaxLz2>=this->LzMax)
    {
      TmpMaxLz2 = this->LzMax-1;
      while ((referenceState & (2l << (TmpMaxLz2<<1))) == 0)
	--TmpMaxLz2;
    }
  test = referenceState;
  // cout << "Reference-state =    " << test << " TmpMaxLz=" <<TmpMaxLz<<", TmpCurrentLzMax="<<TmpCurrentLzMax<<" highestBit="<<highestBit <<", TmpMaxLz2="<<TmpMaxLz2<< endl;
  // cout << "down-up descendents"<<endl;
  while (TmpCurrentLzMax < TmpMaxLz)
    {
      while ((TmpCurrentLzMax < TmpMaxLz) && (((referenceState >> (TmpCurrentLzMax<<1)) & 0x5l) != 0x4l))
	++TmpCurrentLzMax;
      if (TmpCurrentLzMax < TmpMaxLz)
	{
	  TmpReferenceState = (referenceState & ~(0x5l << (TmpCurrentLzMax<<1))) | (0x1l << (TmpCurrentLzMax<<1));
	  TmpCurrentLzMax2 = 0;
	  while (TmpCurrentLzMax2 <= TmpMaxLz2)
	    {
	      while ((TmpCurrentLzMax2 <= TmpMaxLz2) && (((referenceState >> (TmpCurrentLzMax2<<1)) & 0xal) != 0x2l))
		++TmpCurrentLzMax2;
	      if (TmpCurrentLzMax2 <= TmpMaxLz2)
		{
		  TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0xal << (TmpCurrentLzMax2<<1))) | (0x8l << (TmpCurrentLzMax2<<1));
		  
		  test = TmpGeneratedStates2[NbrEntries];
		  // cout << "Descending state " <<NbrEntries << " = " << test;
		  if ((TmpCurrentLzMax2<<1)+3 < highestBit)
		    TmpHighestBit[NbrEntries] = highestBit;
		  else
		    TmpHighestBit[NbrEntries] = ((TmpCurrentLzMax2<<1) + 3);
		  // cout << ", highestBit = " << TmpHighestBit[NbrEntries] << " at TmpCurrentLzMax2="<<TmpCurrentLzMax2<< endl;
		  ++NbrEntries;
		  ++TmpCurrentLzMax2;
		}
	    }
	  ++TmpCurrentLzMax;
	}
    }
  if (((referenceState >> (TmpCurrentLzMax<<1)) & 0x5l) == 0x4l)
    {
      TmpReferenceState = (referenceState & ~(0x5l << (TmpCurrentLzMax<<1))) | (0x1l << (TmpCurrentLzMax<<1));
      TmpCurrentLzMax2 = 0;
      int NewHighestBit;
      if (highestBit&1)
	NewHighestBit = highestBit;
      else
	{
	  NewHighestBit = highestBit-1;
	  if ((TmpReferenceState & (1l << NewHighestBit)) == 0l ) --NewHighestBit;
	}      
      while (TmpCurrentLzMax2 <= TmpMaxLz2)
	{
	  while ((TmpCurrentLzMax2 <= TmpMaxLz2) && (((referenceState >> (TmpCurrentLzMax2<<1)) & 0x5l) != 0x1l))
	    ++TmpCurrentLzMax2;
	  if (TmpCurrentLzMax2 <= TmpMaxLz2)
	    {
	      TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0xal << (TmpCurrentLzMax2<<1))) | (0x8l << (TmpCurrentLzMax2<<1));
	      if ((TmpCurrentLzMax2<<1)+3 < NewHighestBit)
		TmpHighestBit[NbrEntries] = NewHighestBit;
	      else
		TmpHighestBit[NbrEntries] = ((TmpCurrentLzMax2<<1) + 3);
	      test = TmpGeneratedStates2[NbrEntries];
	      // cout << "Descending state " <<NbrEntries << " = " << test << " (highest bit " << TmpHighestBit[NbrEntries] << " at TmpCurrentLzMax="<<TmpCurrentLzMax<< " at TmpCurrentLzMax2="<<TmpCurrentLzMax2<<")" << endl;		  
	      ++NbrEntries;
	      ++TmpCurrentLzMax2;
	    }
	}      
    }
  
  int TmpIndex;
  int NbrNewEntries = 0;
  for (int i = 0; i < NbrEntries; ++i)
    {
      TmpIndex = this->FindStateIndex(TmpGeneratedStates2[i], TmpHighestBit[i]);
#ifdef __64_BITS__
      if ((this->KeepStateFlag[TmpIndex >> 6] >> (TmpIndex & 0x3f)) & 0x1l)
	{
	  TmpGeneratedStates2[i] = 0x0l;
	}
      else
	{
	  this->KeepStateFlag[TmpIndex >> 6] |= 0x1l << (TmpIndex & 0x3f);
	  ++NbrNewEntries;
	}
#else
      if ((this->KeepStateFlag[TmpIndex >> 5] >> (TmpIndex & 0x1f)) & 0x1l)
	{
	  TmpGeneratedStates2[i] = 0x0l;
	}
      else
	{
	  this->KeepStateFlag[TmpIndex >> 5] |= 0x1l << (TmpIndex & 0x1f);
	  ++NbrNewEntries;
	}
#endif
    }

  if (NbrNewEntries > 0)
    for (int i = 0; i < NbrEntries; ++i)
      if (TmpGeneratedStates2[i] != 0x0l)
	this->GenerateDescendingStates(TmpHighestBit[i], TmpGeneratedStates2[i], memory);
  
  memory -= 1;  
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinSqueezedBasis::RawGenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSpin, long pos)
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

void FermionOnSphereWithSpinSqueezedBasis::GenerateLookUpTable(unsigned long memory)
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

double FermionOnSphereWithSpinSqueezedBasis::ComputeSign(unsigned long signs)
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

long FermionOnSphereWithSpinSqueezedBasis::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin)
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

long FermionOnSphereWithSpinSqueezedBasis::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin)
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

Complex FermionOnSphereWithSpinSqueezedBasis::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
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
double FermionOnSphereWithSpinSqueezedBasis::GetStateSign(int index, int* IndicesDown)
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

void FermionOnSphereWithSpinSqueezedBasis::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
  
// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

RealVector FermionOnSphereWithSpinSqueezedBasis::ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace)
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
