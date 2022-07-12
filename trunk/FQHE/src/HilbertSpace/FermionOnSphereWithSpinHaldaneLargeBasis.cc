////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin                    //
//    including the Haldane squeezing technique for large system size         //
//                                                                            //
//                        last modification : 05/09/2009                      //
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
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneLargeBasis.h"
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
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;
using std::ifstream;
using std::ios;

#define WANT_LAPACK

#ifdef __LAPACK__
#ifdef WANT_LAPACK
#define  __USE_LAPACK_HERE__
#endif
#endif


// default constructor
//

FermionOnSphereWithSpinHaldaneLargeBasis::FermionOnSphereWithSpinHaldaneLargeBasis()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twce the total spin value
// rootPartitions = array of root partitions describing the squeezed basis
// nbrRootPartitions = number of root partitions
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinHaldaneLargeBasis::FermionOnSphereWithSpinHaldaneLargeBasis (int nbrFermions, int& totalLz, int lzMax, int& totalSpin, int** rootPartitions, int nbrRootPartitions, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;
  this->LzMax = lzMax;
  this->HighestBit = 2 * this->LzMax + 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;

  this->NbrRootPartitions = nbrRootPartitions;
  this->RootPartitions = new unsigned long [this->NbrRootPartitions];
  for (int j = 0; j < this->NbrRootPartitions; ++j)
    {
      this->RootPartitions[j] = 0x0ul;
      int TmpTotalLz = 0;
      int TmpTotalSpin = 0;
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->RootPartitions[j] |= ((unsigned long) (rootPartitions[j][i] & 3)) << (2 * i);
	  if (rootPartitions[j][i] != 0)
	    {
	      switch (rootPartitions[j][i])
		{
		case 1:
		  {
		    TmpTotalLz += i;
		    --TmpTotalSpin;
		  }
		  break;
		case 2:
		  {
		    TmpTotalLz += i;
		    ++TmpTotalSpin;
		  }
		  break;
		case 3:
		  {
		    TmpTotalLz += i;
		    TmpTotalLz += i;
		  }
		  break;
		}
	    }
	}
      if (j == 0)
	{
	  this->TotalLz = TmpTotalLz;
	  this->TotalSpin = TmpTotalSpin;
	}
      else
	{
	  if (this->TotalLz != TmpTotalLz)
	    cout << "warning : root partition " << j << "does not have the same TotalLz as root partition 0" << endl;
	  if (this->TotalSpin != TmpTotalSpin)
	    cout << "warning : root partition " << j << "does not have the same TotalSpin as root partition 0" << endl;
	}
    }
  this->TotalLz = ((this->TotalLz << 1) - (this->LzMax * this->NbrFermions));
  totalLz = this->TotalLz;
  totalSpin = this->TotalSpin;

  cout << "computing Hilbert space dimension" << endl;
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
										(this->TotalSpin + this->NbrFermions) >> 1);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateHighestBit = 0;  
  cout << "generating Hilbert space (" << this->LargeHilbertSpaceDimension << ")" << endl;
  this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
							  (this->TotalSpin + this->NbrFermions) >> 1, 0l);
  cout << "generating look-up table " << endl;
  this->GenerateLookUpTable(memory);
#ifdef  __64_BITS__
  long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 6) + 1;
#else
  long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 5) + 1;
#endif
  this->KeepStateFlag = new unsigned long [ReducedHilbertSpaceDimension + 1];
  for (long i = 0; i < ReducedHilbertSpaceDimension; ++i)
    this->KeepStateFlag[i] = 0x0l;
  int MaxSweeps = (this->NbrFermions * (this->NbrFermions - 1)) >> 1;  
  this->TmpGeneratedStates =  new unsigned long [MaxSweeps * 1000];
  this->TmpGeneratedStatesLzMax = new int [MaxSweeps * 1000];

  for (int j = 0; j < this->NbrRootPartitions; ++j)
    {
      long Memory = 0l;
      int TmpLzMax = this->HighestBit;
      while (((this->RootPartitions[j] >> TmpLzMax) & 0x1ul) == 0x0ul)
	--TmpLzMax;
      long TmpIndex = this->FindStateIndex(this->RootPartitions[j]);
#ifdef  __64_BITS__
      this->KeepStateFlag[TmpIndex >> 6] |= 0x1l << (TmpIndex & 0x3f);
#else
      this->KeepStateFlag[TmpIndex >> 5] |= 0x1l << (TmpIndex & 0x1f);
#endif
      cout << "generating squeezed states from root partition " << j << endl;
      this->GenerateSqueezedStates(TmpLzMax, this->RootPartitions[j], 1, Memory);  
    }

  long NewHilbertSpaceDimension = 0;
  unsigned long TmpKeepStateFlag;
  long TmpNbrOne[] = {  
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
  for (long i = 0; i < ReducedHilbertSpaceDimension; ++i)
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
  for (int i = 0; i < (2 * this->NbrLzValue); ++i)
    delete[] this->LargeLookUpTable[i];
  delete[] this->LargeLookUpTable;
  unsigned long* TmpStateDescription = new unsigned long [NewHilbertSpaceDimension];
  NewHilbertSpaceDimension = 0l;
  long TotalIndex = 0l;
#ifdef  __64_BITS__
  if ((this->LargeHilbertSpaceDimension & 0x3fl) != 0)
#else
  if ((this->LargeHilbertSpaceDimension & 0x1fl) != 0)
#endif
    --ReducedHilbertSpaceDimension;
  for (long i = 0; i < ReducedHilbertSpaceDimension; ++i)
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
	      ++NewHilbertSpaceDimension;
	    }
	  ++TotalIndex;
	}
    }
#ifdef  __64_BITS__
  this->LargeHilbertSpaceDimension &= 0x3fl;
 #else
  this->LargeHilbertSpaceDimension &= 0x1fl;
 #endif
  if (this->LargeHilbertSpaceDimension != 0l)
    {
      TmpKeepStateFlag = this->KeepStateFlag[ReducedHilbertSpaceDimension];
      for (long j = 0; j < this->LargeHilbertSpaceDimension; ++j)
	{
	  if ((TmpKeepStateFlag >> j) & 0x1l)
	    {
	      TmpStateDescription[NewHilbertSpaceDimension] =  this->StateDescription[TotalIndex];
	      ++NewHilbertSpaceDimension;
	    }
	  ++TotalIndex;
	}
    }
  
  delete[] this->StateDescription;
  delete[] this->KeepStateFlag;
  this->StateDescription = TmpStateDescription;
  this->LargeHilbertSpaceDimension = NewHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

  delete[] this->TmpGeneratedStates;
  delete[] this->TmpGeneratedStatesLzMax;

  this->GenerateLookUpTable(memory);

  
#ifdef __DEBUG__
  long UsedMemory = 0;
  UsedMemory += this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024l)
    if (UsedMemory >= 1048576l)
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

FermionOnSphereWithSpinHaldaneLargeBasis::FermionOnSphereWithSpinHaldaneLargeBasis (char* fileName, unsigned long memory)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      this->HilbertSpaceDimension = 0;
      return;
    }
  File.seekg (0l, ios::end);
  unsigned long FileSize = File.tellg ();
  File.close();

  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      this->HilbertSpaceDimension = 0;
      return;
    }
 

  ReadLittleEndian(File, this->HilbertSpaceDimension);
  ReadLittleEndian(File, this->LargeHilbertSpaceDimension);
  ReadLittleEndian(File, this->NbrFermions);
  ReadLittleEndian(File, this->LzMax);
  ReadLittleEndian(File, this->TotalLz);
  ReadLittleEndian(File, this->TotalSpin);
  ReadLittleEndian(File, this->NbrRootPartitions);
  this->RootPartitions = new unsigned long [this->NbrRootPartitions];
  for (int j = 0; j < this->NbrRootPartitions; ++j)
    ReadLittleEndian(File, this->RootPartitions[j]);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    ReadLittleEndian(File, this->StateDescription[i]);
  File.close();

  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;
  this->HighestBit = 2 * this->LzMax + 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  this->StateHighestBit = 0;  

  this->GenerateLookUpTable(memory);

#ifdef __DEBUG__
  long UsedMemory = 0;
  UsedMemory += this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024l)
    if (UsedMemory >= 1048576l)
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

FermionOnSphereWithSpinHaldaneLargeBasis::FermionOnSphereWithSpinHaldaneLargeBasis(const FermionOnSphereWithSpinHaldaneLargeBasis& fermions)
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
  this->LargeLookUpTable = fermions.LargeLookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrRootPartitions = fermions.NbrRootPartitions;
  this->RootPartitions = new unsigned long [this->NbrRootPartitions];
  for (int i = 0; i < this->NbrRootPartitions; ++i)
    this->RootPartitions[i] = fermions.RootPartitions[i];
}

// destructor
//

FermionOnSphereWithSpinHaldaneLargeBasis::~FermionOnSphereWithSpinHaldaneLargeBasis ()
{
  delete[] this->RootPartitions;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinHaldaneLargeBasis& FermionOnSphereWithSpinHaldaneLargeBasis::operator = (const FermionOnSphereWithSpinHaldaneLargeBasis& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
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
  this->NbrRootPartitions = fermions.NbrRootPartitions;
  delete[] this->RootPartitions;
  this->RootPartitions = new unsigned long [this->NbrRootPartitions];
  for (int i = 0; i < this->NbrRootPartitions; ++i)
    this->RootPartitions[i] = fermions.RootPartitions[i];
  this->LargeLookUpTable = fermions.LargeLookUpTable;  
  return *this;
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphereWithSpinHaldaneLargeBasis::WriteHilbertSpace (char* fileName)
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
  WriteLittleEndian(File, this->NbrRootPartitions);
  for (int j = 0; j < this->NbrRootPartitions; ++j)
    WriteLittleEndian(File, this->RootPartitions[j]);
  File.close();
  return true;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinHaldaneLargeBasis::Clone()
{
  return new FermionOnSphereWithSpinHaldaneLargeBasis(*this);
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereWithSpinHaldaneLargeBasis::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  return -1;
}


// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

long FermionOnSphereWithSpinHaldaneLargeBasis::FindStateIndex(unsigned long stateDescription)
{
  int TmpLzMax = this->HighestBit;
  while ((stateDescription >> TmpLzMax) == 0x0ul)
    --TmpLzMax;
  long PosMax = stateDescription >> this->LookUpTableShift[TmpLzMax];
  long PosMin = this->LargeLookUpTable[TmpLzMax][PosMax];
  PosMax = this->LargeLookUpTable[TmpLzMax][PosMax + 1];
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
      return this->LargeHilbertSpaceDimension;
    else
      return PosMin;
}


// generate all squeezed states from a root partition
// 
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinHaldaneLargeBasis::GenerateSqueezedStates(int lzMax, unsigned long referenceState, long pos, long& memory)
{
  int MaxSweeps = (this->NbrFermions * (this->NbrFermions - 1)) >> 1;  
  unsigned long* TmpGeneratedStates2 = this->TmpGeneratedStates + (MaxSweeps * memory);
  int* TmpLzMax = this->TmpGeneratedStatesLzMax  + (MaxSweeps  * memory);  
  memory += 1;
  int TmpCurrentLzMax = 1;
  int TmpCurrentLzMax2;
  int TmpMax = lzMax - 2;
  int NbrEntries = 0;
  unsigned long TmpReferenceState;  
  while (TmpCurrentLzMax < TmpMax)
    {
      while ((TmpCurrentLzMax < TmpMax) && (((referenceState >> TmpCurrentLzMax) & 0x5l) != 0x4l))
	++TmpCurrentLzMax;
      if (TmpCurrentLzMax < TmpMax)
	{
	  TmpReferenceState = (referenceState & ~(0x5l << TmpCurrentLzMax)) | (0x1l << TmpCurrentLzMax);
	  TmpCurrentLzMax2 = TmpCurrentLzMax - 1;
	  while (TmpCurrentLzMax2 >= 0)
	    {
	      while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> TmpCurrentLzMax2) & 0x5l) != 0x1l))
		--TmpCurrentLzMax2;
	      if (TmpCurrentLzMax2 >= 0)
		{
		  TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0x5l << TmpCurrentLzMax2)) | (0x4l << TmpCurrentLzMax2);
		  TmpLzMax[NbrEntries] = lzMax;
		  ++NbrEntries;
		  --TmpCurrentLzMax2;
		}	      
	    }
	  ++TmpCurrentLzMax;
	}
    }
  if (((referenceState >> TmpCurrentLzMax) & 0x5l) == 0x4l)
    {
      TmpReferenceState = (referenceState & ~(0x5l << TmpCurrentLzMax)) | (0x1l << TmpCurrentLzMax);
      TmpCurrentLzMax2 = TmpCurrentLzMax - 3;
      while (TmpCurrentLzMax2 >= 0)
	{
	  while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> TmpCurrentLzMax2) & 0x5l) != 0x1l))
	    --TmpCurrentLzMax2;
	  if (TmpCurrentLzMax2 >= 0)
	    {
	      TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0x5l << TmpCurrentLzMax2)) | (0x4l << TmpCurrentLzMax2);
	      if ((referenceState & (0x1ul << (lzMax - 1))) == 0x0ul)
		TmpLzMax[NbrEntries] = lzMax - 2;
	      else
		TmpLzMax[NbrEntries] = lzMax - 1;
	      ++NbrEntries;
	      --TmpCurrentLzMax2;
	    }
	}      
    }

  long TmpIndex;
  int NbrNewEntries = 0;
  for (int i = 0; i < NbrEntries; ++i)
    {
      unsigned long& TmpState = TmpGeneratedStates2[i];
      TmpIndex = this->FindStateIndex(TmpState);
#ifdef __64_BITS__
      if ((this->KeepStateFlag[TmpIndex >> 6] >> (TmpIndex & 0x3fl)) & 0x1l)
	{
	  TmpState = 0x0l;
	}
      else
	{
	  this->KeepStateFlag[TmpIndex >> 6] |= 0x1l << (TmpIndex & 0x3fl);
	  ++NbrNewEntries;
	}
#else
      if ((this->KeepStateFlag[TmpIndex >> 5] >> (TmpIndex & 0x1fl)) & 0x1l)
	{
	  TmpState = 0x0l;
	}
      else
	{
	  this->KeepStateFlag[TmpIndex >> 5] |= 0x1l << (TmpIndex & 0x1fl);
	  ++NbrNewEntries;
	}      
#endif
    }

  if (NbrNewEntries > 0)
    for (int i = 0; i < NbrEntries; ++i)
      if (TmpGeneratedStates2[i] != 0x0l)
	pos = this->GenerateSqueezedStates(TmpLzMax[i], TmpGeneratedStates2[i], pos, memory);

  memory -= 1;
  return pos;
}


// convert a gien state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

// RealVector FermionOnSphereWithSpinHaldaneLargeBasis::ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis)
// {
//   RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
//   for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//     TmpVector[nbodyBasis.FindStateIndex(this->StateDescription[i])] = state[i];
//   return TmpVector;
// }

// convert a given state from the usual n-body basis to the Haldane basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

// RealVector FermionOnSphereWithSpinHaldaneLargeBasis::ConvertFromNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis)
// {
//   RealVector TmpVector (this->HilbertSpaceDimension, true);
//   for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//     TmpVector[i] = state[nbodyBasis.FindStateIndex(this->StateDescription[i])];
//   TmpVector /= TmpVector.Norm();
//   return TmpVector;
// }

// create the Jack polynomial decomposition corresponding to the root partition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& FermionOnSphereWithSpinHaldaneLargeBasis::GenerateJackPolynomial(RealVector& jack, double alpha)
{
  for (int i = 0; i < this->NbrRootPartitions; ++i)
    {
      int TmpLzMax = 2 * this->LzMax + 1;
      while (((this->RootPartitions[i] >> TmpLzMax) & 0x1ul) == 0x0ul)
	--TmpLzMax;
      long TmpIndex = this->FindStateIndex(this->RootPartitions[i]);
      if (TmpIndex < this->LargeHilbertSpaceDimension)
	jack[TmpIndex] = 1.0;
    }
  double InvAlpha= -2.0;
  int ReducedNbrFermionsUp = this->NbrFermionsUp - 1;
  int ReducedNbrFermionsDown = this->NbrFermionsDown - 1;


  unsigned long MaxRoot = this->StateDescription[0];
  int* TmpMonomialUp = new int [this->NbrFermionsUp];
  int* TmpMonomialDown = new int [this->NbrFermionsDown];
  int* TmpMonomialUp2 = new int [this->NbrFermionsUp];
  int* TmpMonomialDown2 = new int [this->NbrFermionsDown];
  this->ConvertToMonomial(MaxRoot, TmpMonomialUp, TmpMonomialDown);
  double RhoRoot = 0.0;
  for (int j = 0; j < this->NbrFermionsUp; ++j)
    {
      double Tmp = TmpMonomialUp[j];
      RhoRoot += Tmp * (Tmp - InvAlpha * ((double) (j + 1)));
      for (int k = 0; k < this->NbrFermionsDown; ++k)
	if (Tmp > TmpMonomialDown[k])
	  RhoRoot -= (double) (Tmp - TmpMonomialDown[k]);
    }
  for (int j = 0; j < this->NbrFermionsDown; ++j)
    {
      double Tmp = TmpMonomialDown[j];
      RhoRoot += Tmp * (Tmp - InvAlpha * ((double) (j + 1)));
      for (int k = 0; k < this->NbrFermionsUp; ++k)
	if (Tmp > TmpMonomialUp[k])
	  RhoRoot -= (double) (Tmp - TmpMonomialUp[k]);
    }

  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    if (jack[i] == 0.0)
      {      
	unsigned long CurrentPartition = this->StateDescription[i];
	this->ConvertToMonomial(CurrentPartition, TmpMonomialUp, TmpMonomialDown);
	double Rho = 0.0;
	for (int j = 0; j < this->NbrFermionsUp; ++j)
	  {
	    double Tmp = TmpMonomialUp[j];
	    Rho += Tmp * (Tmp - InvAlpha * ((double) (j + 1)));
	    for (int k = 0; k < this->NbrFermionsDown; ++k)
	      if (Tmp > TmpMonomialDown[k])
		Rho -= (double) (Tmp - TmpMonomialDown[k]);
	  }
	for (int j = 0; j < this->NbrFermionsDown; ++j)
	  {
	    double Tmp = TmpMonomialDown[j];
	    Rho += Tmp * (Tmp - InvAlpha * ((double) (j + 1)));
	    for (int k = 0; k < this->NbrFermionsUp; ++k)
	      if (Tmp > TmpMonomialUp[k])
		Rho -= (double) (Tmp - TmpMonomialUp[k]);
	  }
	
	double Factor = 1.0;
#ifdef __64_BITS__
	if (((CurrentPartition >> 1) & 0x5555555555555555ul) == ((CurrentPartition & 0x5555555555555555ul)))
#else
	  if (((CurrentPartition >> 1) & 0x55555555ul) == ((CurrentPartition & 0x55555555ul)))
#endif
	    Factor = 0.5;
	
	double Coefficient = 0.0;
	
	for (int j1 = 0; j1 < ReducedNbrFermionsUp; ++j1)
	  for (int j2 = j1 + 1; j2 < this->NbrFermionsUp; ++j2)
	    {
	      double Diff = (double) (TmpMonomialUp[j1] - TmpMonomialUp[j2]);
	      unsigned int Max = TmpMonomialUp[j2];
	      unsigned long TmpState = CurrentPartition;
	      int Tmpj1 = j1;
	      int Tmpj2 = j2;
	      for (int l = 0; l < this->NbrFermionsUp; ++l)
		TmpMonomialUp2[l] = TmpMonomialUp[l];	    
	      double Sign = 1.0;
	      for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		{
		  ++TmpMonomialUp2[Tmpj1];
		  --TmpMonomialUp2[Tmpj2];
		  Diff += 2.0;
		  while ((Tmpj1 > 0) && (TmpMonomialUp2[Tmpj1] >= TmpMonomialUp2[Tmpj1 - 1]))
		    {
		      unsigned long Tmp = TmpMonomialUp2[Tmpj1 - 1];
		      TmpMonomialUp2[Tmpj1 - 1] = TmpMonomialUp2[Tmpj1];
		      TmpMonomialUp2[Tmpj1] = Tmp;
		      --Tmpj1;
		      Sign *= -1.0; 
		    }
		  while ((Tmpj2 < ReducedNbrFermionsUp) && (TmpMonomialUp2[Tmpj2] <= TmpMonomialUp2[Tmpj2 + 1]))
		    {
		      unsigned long Tmp = TmpMonomialUp2[Tmpj2 + 1];
		      TmpMonomialUp2[Tmpj2 + 1] = TmpMonomialUp2[Tmpj2];
		      TmpMonomialUp2[Tmpj2] = Tmp;
		      ++Tmpj2;
		      Sign *= -1.0; 
		    }
		  if ((TmpMonomialUp2[Tmpj1] != TmpMonomialUp2[Tmpj1 + 1]) && (TmpMonomialUp2[Tmpj2] != TmpMonomialUp2[Tmpj2 - 1]))
		    {
		      TmpState = this->ConvertFromMonomial(TmpMonomialUp2, TmpMonomialDown);
		      if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
			{
			  int TmpLzMax = TmpMonomialUp2[0];
			  if (TmpLzMax < TmpMonomialDown[0])
			    TmpLzMax = TmpMonomialDown[0] << 1;
			  else
			    {
			      TmpLzMax <<= 1;
			      ++TmpLzMax;
			  }
			  long TmpIndex = this->FindStateIndex(TmpState);
			  if (TmpIndex < this->LargeHilbertSpaceDimension)
			    {
			      Coefficient += Sign * Diff * jack[TmpIndex];
			    }
			}
		    }
		}
	    }
	
	for (int j1 = 0; j1 < ReducedNbrFermionsDown; ++j1)
	  for (int j2 = j1 + 1; j2 < this->NbrFermionsDown; ++j2)
	    {
	      double Diff = (double) (TmpMonomialDown[j1] - TmpMonomialDown[j2]);
	      unsigned int Max = TmpMonomialDown[j2];
	      unsigned long TmpState = CurrentPartition;
	      int Tmpj1 = j1;
	      int Tmpj2 = j2;
	      for (int l = 0; l < this->NbrFermionsDown; ++l)
		TmpMonomialDown2[l] = TmpMonomialDown[l];	    
	      double Sign = 1.0;
	      for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		{
		  ++TmpMonomialDown2[Tmpj1];
		  --TmpMonomialDown2[Tmpj2];
		  Diff += 2.0;
		  while ((Tmpj1 > 0) && (TmpMonomialDown2[Tmpj1] >= TmpMonomialDown2[Tmpj1 - 1]))
		    {
		      unsigned long Tmp = TmpMonomialDown2[Tmpj1 - 1];
		      TmpMonomialDown2[Tmpj1 - 1] = TmpMonomialDown2[Tmpj1];
		      TmpMonomialDown2[Tmpj1] = Tmp;
		      --Tmpj1;
		      Sign *= -1.0; 
		    }
		  while ((Tmpj2 < ReducedNbrFermionsDown) && (TmpMonomialDown2[Tmpj2] <= TmpMonomialDown2[Tmpj2 + 1]))
		    {
		      unsigned long Tmp = TmpMonomialDown2[Tmpj2 + 1];
		      TmpMonomialDown2[Tmpj2 + 1] = TmpMonomialDown2[Tmpj2];
		      TmpMonomialDown2[Tmpj2] = Tmp;
		      ++Tmpj2;
		      Sign *= -1.0; 
		    }
		  if ((TmpMonomialDown2[Tmpj1] != TmpMonomialDown2[Tmpj1 + 1]) && (TmpMonomialDown2[Tmpj2] != TmpMonomialDown2[Tmpj2 - 1]))
		    {
		      TmpState = this->ConvertFromMonomial(TmpMonomialUp, TmpMonomialDown2);
		      if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
			{
			  int TmpLzMax = TmpMonomialUp[0];
			  if (TmpLzMax < TmpMonomialDown2[0])
			    TmpLzMax = TmpMonomialDown2[0] << 1;
			  else
			    {
			      TmpLzMax <<= 1;
			      ++TmpLzMax;
			    }
			  long TmpIndex = this->FindStateIndex(TmpState);
			  if (TmpIndex < this->LargeHilbertSpaceDimension)
			    {
			      Coefficient += Sign * Diff * jack[TmpIndex];
			    }
			}
		    }
		}
	    }
	
	for (int j1 = 0; j1 < this->NbrFermionsUp; ++j1)
	  {
	    int j2 = 0;
	    while ((j2 < this->NbrFermionsDown) && (TmpMonomialDown[j2] > TmpMonomialUp[j1]))
	      ++j2;
	    for (; j2 < this->NbrFermionsDown; ++j2)
	      {
		double Diff = (double) (TmpMonomialUp[j1] - TmpMonomialDown[j2]);
		unsigned int Max = TmpMonomialDown[j2];
		unsigned long TmpState = CurrentPartition;
		int Tmpj1 = j1;
		int Tmpj2 = j2;
		for (int l = 0; l < this->NbrFermionsUp; ++l)
		  TmpMonomialUp2[l] = TmpMonomialUp[l];	    
		for (int l = 0; l < this->NbrFermionsDown; ++l)
		  TmpMonomialDown2[l] = TmpMonomialDown[l];	    
		double Sign = 1.0;
		for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		  {
		    ++TmpMonomialUp2[Tmpj1];
		    --TmpMonomialDown2[Tmpj2];
		    Diff += 2.0;
		    while ((Tmpj1 > 0) && (TmpMonomialUp2[Tmpj1] >= TmpMonomialUp2[Tmpj1 - 1]))
		      {
			unsigned long Tmp = TmpMonomialUp2[Tmpj1 - 1];
			TmpMonomialUp2[Tmpj1 - 1] = TmpMonomialUp2[Tmpj1];
			TmpMonomialUp2[Tmpj1] = Tmp;
			--Tmpj1;
			Sign *= -1.0; 
		      }
		    while ((Tmpj2 < ReducedNbrFermionsDown) && (TmpMonomialDown2[Tmpj2] <= TmpMonomialDown2[Tmpj2 + 1]))
		      {
			unsigned long Tmp = TmpMonomialDown2[Tmpj2 + 1];
			TmpMonomialDown2[Tmpj2 + 1] = TmpMonomialDown2[Tmpj2];
			TmpMonomialDown2[Tmpj2] = Tmp;
			++Tmpj2;
			Sign *= -1.0; 
		      }
		    if (((Tmpj1 == ReducedNbrFermionsUp) || (TmpMonomialUp2[Tmpj1] != TmpMonomialUp2[Tmpj1 + 1])) &&
			((Tmpj2 == 0) || (TmpMonomialDown2[Tmpj2] != TmpMonomialDown2[Tmpj2 - 1])))
		      {
			TmpState = this->ConvertFromMonomial(TmpMonomialUp2, TmpMonomialDown2);
			if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
			  {
			    int TmpLzMax = TmpMonomialUp2[0];
			    if (TmpLzMax < TmpMonomialDown2[0])
			      TmpLzMax = TmpMonomialDown2[0] << 1;
			    else
			      {
				TmpLzMax <<= 1;
				++TmpLzMax;
			      }
			    long TmpIndex = this->FindStateIndex(TmpState);
			    if (TmpIndex < this->LargeHilbertSpaceDimension)
			      {
				Coefficient += Sign * Diff * jack[TmpIndex];
			      }
			  }
		      }
		  }
	      }
	  }
	
	for (int j1 = 0; j1 < this->NbrFermionsDown; ++j1)
	  {
	    int j2 = 0;
	    while ((j2 < this->NbrFermionsUp) && (TmpMonomialUp[j2] >= TmpMonomialDown[j1]))
	      ++j2;
	    for (; j2 < this->NbrFermionsUp; ++j2)
	      {
		double Diff = (double) (TmpMonomialDown[j1] - TmpMonomialUp[j2]);
		unsigned int Max = TmpMonomialUp[j2];
		unsigned long TmpState = CurrentPartition;
		int Tmpj1 = j1;
		int Tmpj2 = j2;
		for (int l = 0; l < this->NbrFermionsDown; ++l)
		  TmpMonomialDown2[l] = TmpMonomialDown[l];	    
		for (int l = 0; l < this->NbrFermionsUp; ++l)
		  TmpMonomialUp2[l] = TmpMonomialUp[l];	    
		double Sign = 1.0;
		for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		  {
		    ++TmpMonomialDown2[Tmpj1];
		    --TmpMonomialUp2[Tmpj2];
		    Diff += 2.0;
		    while ((Tmpj1 > 0) && (TmpMonomialDown2[Tmpj1] >= TmpMonomialDown2[Tmpj1 - 1]))
		      {
			unsigned long Tmp = TmpMonomialDown2[Tmpj1 - 1];
			TmpMonomialDown2[Tmpj1 - 1] = TmpMonomialDown2[Tmpj1];
			TmpMonomialDown2[Tmpj1] = Tmp;
			--Tmpj1;
			Sign *= -1.0; 
		      }
		    while ((Tmpj2 < ReducedNbrFermionsUp) && (TmpMonomialUp2[Tmpj2] <= TmpMonomialUp2[Tmpj2 + 1]))
		      {
			unsigned long Tmp = TmpMonomialUp2[Tmpj2 + 1];
			TmpMonomialUp2[Tmpj2 + 1] = TmpMonomialUp2[Tmpj2];
			TmpMonomialUp2[Tmpj2] = Tmp;
			++Tmpj2;
			Sign *= -1.0; 
		      }
		    if (((Tmpj1 == ReducedNbrFermionsDown) || (TmpMonomialDown2[Tmpj1] != TmpMonomialDown2[Tmpj1 + 1]))
			&& ((Tmpj2 == 0) || (TmpMonomialUp2[Tmpj2] != TmpMonomialUp2[Tmpj2 - 1])))
		      {
			TmpState = this->ConvertFromMonomial(TmpMonomialUp2, TmpMonomialDown2);
			if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
			  {
			    int TmpLzMax = TmpMonomialUp2[0];
			    if (TmpLzMax < TmpMonomialDown2[0])
			      TmpLzMax = TmpMonomialDown2[0] << 1;
			    else
			      {
				TmpLzMax <<= 1;
				++TmpLzMax;
			      }
			    long TmpIndex = this->FindStateIndex(TmpState);
			    if (TmpIndex < this->LargeHilbertSpaceDimension)
			      {
				Coefficient += Sign * Diff * jack[TmpIndex];
			      }
			  }
		      }
		  }
	      }
	  }
	
	jack[i] = Coefficient * InvAlpha / (RhoRoot - Rho);
	if ((i & 0xffl) == 0l)
	  {
	    cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	  }
      }
  

  double GlobalSign = this->GetSpinSeparationSign(this->StateDescription[0]);  
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      jack[i] *= GlobalSign * this->GetSpinSeparationSign(this->StateDescription[i]);
    }
  return jack;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereWithSpinHaldaneLargeBasis::GenerateLookUpTable(unsigned long memory)
{
  // get every highest bit poisition
  unsigned long TmpState = this->StateDescription[0];
#ifdef __64_BITS__
  int CurrentHighestBit = 63;
#else
  int CurrentHighestBit = 31;
#endif
  while ((TmpState & (0x1ul << CurrentHighestBit)) == 0x0ul)
    --CurrentHighestBit;  

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
  this->LargeLookUpTable = new long* [2*this->NbrLzValue];
  this->LookUpTableShift = new int [2*this->NbrLzValue];
  for (int i = 0; i < 2*this->NbrLzValue; ++i)
    this->LargeLookUpTable[i] = new long [this->LookUpTableMemorySize + 1];
  int CurrentLargestBit = CurrentHighestBit;
  long* TmpLookUpTable = this->LargeLookUpTable[CurrentLargestBit];
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
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {     
      TmpState = this->StateDescription[i];
      while ((TmpState & (0x1ul << CurrentHighestBit)) == 0x0ul)
	--CurrentHighestBit;  
      if (CurrentLargestBit != CurrentHighestBit)
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
 	  CurrentLargestBit = CurrentHighestBit;
	  TmpLookUpTable = this->LargeLookUpTable[CurrentLargestBit];
	  if (CurrentLargestBit < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentLargestBit] = 0;
	  else
	    this->LookUpTableShift[CurrentLargestBit] = CurrentLargestBit + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentLargestBit];
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
      TmpLookUpTable[CurrentLookUpTableValue] = this->LargeHilbertSpaceDimension - 1l;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->LargeHilbertSpaceDimension - 1l;

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


