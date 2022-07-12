////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//   that allow LzMax up to 127 (for systems with 128 bit integer support)    //
//               or 63 (on 32 bit systems without 128 bit integer support)    //
//                                                                            //
//                        last modification : 22/09/2007                      //
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
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/Endian.h"
#include "Polynomial/RationalPolynomial.h"
#include "Polynomial/LongRationalPolynomial.h"
#include "Architecture/ArchitectureOperation/FQHESphereJackGeneratorSumRationalPolynomialOperation.h"

#include <math.h>
#include <fstream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
//

FermionOnSphereHaldaneBasisLong::FermionOnSphereHaldaneBasisLong()
{
  this->HilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = reference on twice the momentum total value, totalLz will be recomputed from referenceState and stored in totalLz
// lzMax = twice the maximum Lz value reached by a fermion
// memory = amount of memory granted for precalculations
// referenceState = array that describes the reference state to start from

FermionOnSphereHaldaneBasisLong::FermionOnSphereHaldaneBasisLong (int nbrFermions, int& totalLz, int lzMax, int* referenceState,
								  unsigned long memory)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->ReferenceState = ((ULONGLONG) 0x0ul);
  int ReferenceStateLzMax = 0;
  this->TotalLz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      this->ReferenceState |= ((ULONGLONG) (referenceState[i] & 1)) << i;
      if (referenceState[i] == 1)
	{
	  ReferenceStateLzMax = i;
	  this->TotalLz += i;
	}
    }
  this->TotalLz = ((this->TotalLz << 1) - (this->LzMax * this->NbrFermions));
  totalLz = this->TotalLz;
#ifdef __128_BIT_LONGLONG__
  this->InvertShift = 64 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 32 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();

  this->StateDescription = new ULONGLONG [this->LargeHilbertSpaceDimension];
  this->StateLzMax = new int [this->LargeHilbertSpaceDimension];
#ifdef  __64_BITS__
  int ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 6) + 1;
#else
  int ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 5) + 1;
#endif
  this->KeepStateFlag = new unsigned long [ReducedHilbertSpaceDimension];
  for (int i = 0; i < ReducedHilbertSpaceDimension; ++i)
    this->KeepStateFlag[i] = 0x0l;
  this->RawGenerateStates(this->NbrFermions, this->LzMax, this->LzMax, (this->TotalLz + this->NbrFermions * this->LzMax) >> 1, 0);
  this->GenerateLookUpTable(memory);
  int MaxSweeps = (this->NbrFermions * (this->NbrFermions - 1)) >> 1;  
  this->TmpGeneratedStates =  new ULONGLONG [MaxSweeps * 1000];
  this->TmpGeneratedStatesLzMax = new int [MaxSweeps * 1000];
  long Memory = 0l;

  int TmpIndex = this->FindStateIndex(this->ReferenceState, ReferenceStateLzMax);
#ifdef  __64_BITS__
  this->KeepStateFlag[TmpIndex >> 6] = 0x1l << (TmpIndex & 0x3f);
#else
  this->KeepStateFlag[TmpIndex >> 5] = 0x1l << (TmpIndex & 0x1f);
#endif
  this->GenerateStates(ReferenceStateLzMax, this->ReferenceState, 1, Memory);  

  long NewHilbertSpaceDimension = 0l;
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
  ULONGLONG* TmpStateDescription = new ULONGLONG [NewHilbertSpaceDimension];
  int* TmpStateLzMax = new int [NewHilbertSpaceDimension];
  NewHilbertSpaceDimension = 0l;
  int TotalIndex = 0;
#ifdef  __64_BITS__
  if ((this->LargeHilbertSpaceDimension & 0x3fl) != 0)
#else
  if ((this->LargeHilbertSpaceDimension & 0x1fl) != 0)
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
	      TmpStateLzMax[NewHilbertSpaceDimension] = this->StateLzMax[TotalIndex];
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
      for (long j = 0l; j < this->LargeHilbertSpaceDimension; ++j)
	{
	  if ((TmpKeepStateFlag >> j) & 0x1l)
	    {
	      TmpStateDescription[NewHilbertSpaceDimension] =  this->StateDescription[TotalIndex];
	      TmpStateLzMax[NewHilbertSpaceDimension] = this->StateLzMax[TotalIndex];
	      ++NewHilbertSpaceDimension;
	    }
	  ++TotalIndex;
	}
    }
  
  delete[] this->StateDescription;
  delete[] this->StateLzMax;
  delete[] this->KeepStateFlag;
  this->StateDescription = TmpStateDescription;
  this->StateLzMax = TmpStateLzMax;
  this->LargeHilbertSpaceDimension = NewHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

  delete[] this->TmpGeneratedStates;
  delete[] this->TmpGeneratedStatesLzMax;

  this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
  unsigned long UsedMemory = 0l;
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(ULONGLONG) + sizeof(int));
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

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

FermionOnSphereHaldaneBasisLong::FermionOnSphereHaldaneBasisLong (char* fileName, unsigned long memory)
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
  ReadLittleEndian(File, this->LargeHilbertSpaceDimension);
  ReadLittleEndian(File, this->NbrFermions);
  ReadLittleEndian(File, this->LzMax);
  ReadLittleEndian(File, this->TotalLz);
  ReadLittleEndian(File, this->ReferenceState);
  this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    ReadLittleEndian(File, this->StateDescription[i]);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    ReadLittleEndian(File, this->StateLzMax[i]);

  File.close();

  this->TargetSpace = this;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->Flag.Initialize();
#ifdef __128_BIT_LONGLONG__
  this->InvertShift = 64 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 32 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;

  this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
  unsigned long UsedMemory = 0l;
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(ULONGLONG) + sizeof(int));
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

FermionOnSphereHaldaneBasisLong::FermionOnSphereHaldaneBasisLong(const FermionOnSphereHaldaneBasisLong& fermions)
{
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->ReferenceState = fermions.ReferenceState;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
}

// destructor
//

FermionOnSphereHaldaneBasisLong::~FermionOnSphereHaldaneBasisLong ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereHaldaneBasisLong& FermionOnSphereHaldaneBasisLong::operator = (const FermionOnSphereHaldaneBasisLong& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->ReferenceState = fermions.ReferenceState;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereHaldaneBasisLong::Clone()
{
  return new FermionOnSphereHaldaneBasisLong(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphereHaldaneBasisLong::WriteHilbertSpace (char* fileName)
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
  WriteLittleEndian(File, this->ReferenceState);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    WriteLittleEndian(File, this->StateDescription[i]);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    WriteLittleEndian(File, this->StateLzMax[i]);
  File.close();
  return true;
}

// convert a gien state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereHaldaneBasisLong::ConvertToNbodyBasis(RealVector& state, FermionOnSphereLong& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    TmpVector[nbodyBasis.FindStateIndex(this->StateDescription[i], this->StateLzMax[i])] = state[i];
  return TmpVector;
}


// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereHaldaneBasisLong::FindStateIndex(ULONGLONG stateDescription, int lzmax)
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
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}


// generate all states corresponding to the constraints
// 
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int FermionOnSphereHaldaneBasisLong::GenerateStates(int lzMax, ULONGLONG referenceState, int pos, long& memory)
{
  int MaxSweeps = (this->NbrFermions * (this->NbrFermions - 1)) >> 1;  
  ULONGLONG* TmpGeneratedStates2 = this->TmpGeneratedStates + (MaxSweeps * memory);
  int* TmpLzMax = this->TmpGeneratedStatesLzMax  + (MaxSweeps  * memory);  
  memory += 1;
  int TmpCurrentLzMax = 2;
  int TmpCurrentLzMax2;
  int TmpMax = lzMax - 1;
  int NbrEntries = 0;
  ULONGLONG TmpReferenceState;
  
  while (TmpCurrentLzMax < TmpMax)
    {
      while ((TmpCurrentLzMax < TmpMax) && (((referenceState >> TmpCurrentLzMax) & ((ULONGLONG) 0x3ul)) != ((ULONGLONG) 0x2ul)))
	++TmpCurrentLzMax;
      if (TmpCurrentLzMax < TmpMax)
	{
	  TmpReferenceState = (referenceState & ~(((ULONGLONG) 0x3ul) << TmpCurrentLzMax)) | (((ULONGLONG) 0x1ul) << TmpCurrentLzMax);
	  TmpCurrentLzMax2 = TmpCurrentLzMax - 2;
	  while (TmpCurrentLzMax2 >= 0)
	    {
	      while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> TmpCurrentLzMax2) & ((ULONGLONG) 0x3ul)) != ((ULONGLONG) 0x1ul)))
		--TmpCurrentLzMax2;
	      if (TmpCurrentLzMax2 >= 0)
		{
		  TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(((ULONGLONG) 0x3ul) << TmpCurrentLzMax2)) | (((ULONGLONG) 0x2ul) << TmpCurrentLzMax2);
		  TmpLzMax[NbrEntries] = lzMax;
		  ++NbrEntries;
		  --TmpCurrentLzMax2;
		}	      
	    }
	  ++TmpCurrentLzMax;
	}
    }
  if (((referenceState >> TmpCurrentLzMax) & ((ULONGLONG) 0x3ul)) == ((ULONGLONG) 0x2ul))
    {
      TmpReferenceState = (referenceState & ~(((ULONGLONG) 0x3ul) << TmpCurrentLzMax)) | (((ULONGLONG) 0x1ul) << TmpCurrentLzMax);
      TmpCurrentLzMax2 = TmpCurrentLzMax - 2;
      while (TmpCurrentLzMax2 >= 0)
	{
	  while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> TmpCurrentLzMax2) & ((ULONGLONG) 0x3ul)) != ((ULONGLONG) 0x1ul)))
	    --TmpCurrentLzMax2;
	  if (TmpCurrentLzMax2 >= 0)
	    {
	      TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(((ULONGLONG) 0x3ul) << TmpCurrentLzMax2)) | (((ULONGLONG) 0x2ul) << TmpCurrentLzMax2);
	      TmpLzMax[NbrEntries] = lzMax - 1;
	      ++NbrEntries;
	      --TmpCurrentLzMax2;
	    }
	}      
    }

  int TmpIndex;
  int NbrNewEntries = 0;
  for (int i = 0; i < NbrEntries; ++i)
    {
      TmpIndex = this->FindStateIndex(TmpGeneratedStates2[i], TmpLzMax[i]);
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
	pos = this->GenerateStates(TmpLzMax[i], TmpGeneratedStates2[i], pos, memory);

  memory -= 1;
  return pos;
}


// generate all states (i.e. all possible skew symmetric polynomials with fixed Lz)
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// currentLzMax = momentum maximum value for fermions that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int FermionOnSphereHaldaneBasisLong::RawGenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int pos)
{
  if ((nbrFermions == 0) || (totalLz < 0) || (currentLzMax < (nbrFermions - 1)))
    return pos;
  int LzTotalMax = ((2 * currentLzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (totalLz > LzTotalMax)
    return pos;
  if ((nbrFermions == 1) && (currentLzMax >= totalLz))
    {
      this->StateDescription[pos] = ((ULONGLONG) 0x1) << totalLz;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }
  if (LzTotalMax == totalLz)
    {
      ULONGLONG Mask = 0;
      for (int i = currentLzMax - nbrFermions + 1; i <= currentLzMax; ++i)
	Mask |= (((ULONGLONG) 1) << i);
      this->StateDescription[pos] = Mask;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }

  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = this->RawGenerateStates(nbrFermions - 1, lzMax, ReducedCurrentLzMax, totalLz - currentLzMax, pos);
  ULONGLONG Mask = ((ULONGLONG) 1) << currentLzMax;
  for (int i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (lzMax == currentLzMax)
    return this->RawGenerateStates(nbrFermions, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, TmpPos);
  else
    return this->RawGenerateStates(nbrFermions, lzMax, ReducedCurrentLzMax, totalLz, TmpPos);
}

// create the Jack polynomial decomposition corresponding to the root partition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& FermionOnSphereHaldaneBasisLong::GenerateJackPolynomial(RealVector& jack, double alpha)
{
  jack[0] = 1.0;
  double InvAlpha =  2.0 * (1.0 - alpha) / alpha;

  unsigned long* TmpMonomial = new unsigned long [this->NbrFermions];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrFermions];

  double RhoRoot = 0.0;
  ULONGLONG MaxRoot = this->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, TmpMonomial);
  for (int j = 0; j < this->NbrFermions; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - InvAlpha * ((double) j));
  int ReducedNbrFermions = this->NbrFermions - 1;

  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      double Rho = 0.0;
      ULONGLONG CurrentPartition = this->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, TmpMonomial);
      for (int j = 0; j < this->NbrFermions; ++j)
	Rho += TmpMonomial[j] * (TmpMonomial[j] - InvAlpha * ((double) j));
      double Coefficient = 0.0;
      for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
	  {
	    double Diff = (double) (TmpMonomial[j1] - TmpMonomial[j2]);
	    unsigned int Max = TmpMonomial[j2];
	    ULONGLONG TmpState = (ULONGLONG) 0x0ul;
	    int Tmpj1 = j1;
	    int Tmpj2 = j2;
	    for (int l = 0; l < this->NbrFermions; ++l)
	      TmpMonomial2[l] = TmpMonomial[l];	    
	    double Sign = 1.0;
	    for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
	      {
		++TmpMonomial2[Tmpj1];
		--TmpMonomial2[Tmpj2];
		while ((Tmpj1 > 0) && (TmpMonomial2[Tmpj1] >= TmpMonomial2[Tmpj1 - 1]))
		  {
		    unsigned long Tmp = TmpMonomial2[Tmpj1 - 1];
		    TmpMonomial2[Tmpj1 - 1] = TmpMonomial2[Tmpj1];
		    TmpMonomial2[Tmpj1] = Tmp;
		    --Tmpj1;
		    Sign *= -1.0; 
		  }
                while ((Tmpj2 < ReducedNbrFermions) && (TmpMonomial2[Tmpj2] <= TmpMonomial2[Tmpj2 + 1]))
                  {
                    unsigned long Tmp = TmpMonomial2[Tmpj2 + 1];
                    TmpMonomial2[Tmpj2 + 1] = TmpMonomial2[Tmpj2];
                    TmpMonomial2[Tmpj2] = Tmp;
                    ++Tmpj2;
 		    Sign *= -1.0; 
                 }
		if ((TmpMonomial2[Tmpj1] != TmpMonomial2[Tmpj1 + 1]) && (TmpMonomial2[Tmpj2] != TmpMonomial2[Tmpj2 - 1]))
		  {
		    TmpState = this->ConvertFromMonomial(TmpMonomial2);
		    if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		      {
			long TmpIndex = this->FindStateIndex(TmpState, TmpMonomial2[0]);
			if (TmpIndex < this->HilbertSpaceDimension)
			  Coefficient += Sign * Diff * jack[TmpIndex];
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
  delete[] TmpMonomial;
  cout << endl;

  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& FermionOnSphereHaldaneBasisLong::GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha)
{
  jack[0] = 1.0;
  double InvAlpha =  2.0 * (1.0 - alpha) / alpha;

  unsigned long* TmpMonomial = new  unsigned long[this->NbrFermions];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrFermions];

  double RhoRoot = 0.0;
  ULONGLONG MaxRoot = this->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, TmpMonomial);
  for (int j = 0; j < this->NbrFermions; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - InvAlpha * ((double) j));
  int ReducedNbrFermions = this->NbrFermions - 1;  
  double SymSign = 1.0;
  if ((((this->NbrFermions * ReducedNbrFermions) >> 1) & 1) != 0)
    SymSign = -1.0;
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      double Rho = 0.0;
      ULONGLONG CurrentPartition = this->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, TmpMonomial);
      for (int j = 0; j < this->NbrFermions; ++j)
	Rho += TmpMonomial[j] * (TmpMonomial[j] - InvAlpha * ((double) j));
      double Coefficient = 0.0;
      for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
	  {
	    double Diff = (double) (TmpMonomial[j1] - TmpMonomial[j2]);
	    unsigned int Max = TmpMonomial[j2];
	    ULONGLONG TmpState = (ULONGLONG) 0x0ul;
	    int Tmpj1 = j1;
	    int Tmpj2 = j2;
	    for (int l = 0; l < this->NbrFermions; ++l)
	      TmpMonomial2[l] = TmpMonomial[l];	    
	    double Sign = 1.0;
	    for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
	      {
		++TmpMonomial2[Tmpj1];
		--TmpMonomial2[Tmpj2];
		while ((Tmpj1 > 0) && (TmpMonomial2[Tmpj1] >= TmpMonomial2[Tmpj1 - 1]))
		  {
		    unsigned long Tmp = TmpMonomial2[Tmpj1 - 1];
		    TmpMonomial2[Tmpj1 - 1] = TmpMonomial2[Tmpj1];
		    TmpMonomial2[Tmpj1] = Tmp;
		    --Tmpj1;
		    Sign *= -1.0; 
		  }
                while ((Tmpj2 < ReducedNbrFermions) && (TmpMonomial2[Tmpj2] <= TmpMonomial2[Tmpj2 + 1]))
                  {
                    unsigned long Tmp = TmpMonomial2[Tmpj2 + 1];
                    TmpMonomial2[Tmpj2 + 1] = TmpMonomial2[Tmpj2];
                    TmpMonomial2[Tmpj2] = Tmp;
                    ++Tmpj2;
 		    Sign *= -1.0; 
                 }
		if ((TmpMonomial2[Tmpj1] != TmpMonomial2[Tmpj1 + 1]) && (TmpMonomial2[Tmpj2] != TmpMonomial2[Tmpj2 - 1]))
		  {
		    TmpState = this->ConvertFromMonomial(TmpMonomial2);
		    if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		      {
			long TmpIndex = this->FindStateIndex(TmpState, TmpMonomial2[0]);
			if (TmpIndex < this->HilbertSpaceDimension)
			  Coefficient += Sign * Diff * jack[TmpIndex];
		      }
		  }
	      }
	  }

      unsigned long TmpSymState = this->GetSymmetricState(CurrentPartition);
      Coefficient *= InvAlpha;
      Coefficient /= (RhoRoot - Rho);
      if (TmpSymState < CurrentPartition)
	{
	  long TmpIndex = this->FindStateIndex(TmpSymState, this->LzMax - TmpMonomial[ReducedNbrFermions]);
	  jack[TmpIndex] = SymSign *Coefficient;
	}
      jack[i] = Coefficient;
      if ((i & 0xffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  delete[] TmpMonomial;
  cout << endl;

  return jack;
}

// check partitions that may lead to singular coefficient in a given Jack polynomial decomposition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// error = error when comparing two rho values
// return value = vector with non-zero component being rho factor of possible singular coefficients

RealVector& FermionOnSphereHaldaneBasisLong::CheckPossibleSingularCoefficientsInJackPolynomial(RealVector& jack, double alpha, double error)
{
  double InvAlpha =  2.0 * (1.0 - alpha) / alpha;

  unsigned long* TmpMonomial = new unsigned long [this->NbrFermions];

  double RhoRoot = 0.0;
  ULONGLONG MaxRoot = this->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, TmpMonomial);
  for (int j = 0; j < this->NbrFermions; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - InvAlpha * ((double) j));
  jack[0] = RhoRoot;
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      double Rho = 0.0;
      ULONGLONG CurrentPartition = this->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, TmpMonomial);
      for (int j = 0; j < this->NbrFermions; ++j)
	Rho += TmpMonomial[j] * (TmpMonomial[j] - InvAlpha * ((double) j));
      if ((fabs(RhoRoot - Rho) < error) || (fabs(RhoRoot - Rho) < (error * fabs(RhoRoot))))
	jack[i] = Rho;
      else
	jack[i] = 0.0;
    }
  delete[] TmpMonomial;
  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alphaNumerator = numerator of the Jack polynomial alpha coefficient
// alphaDenominator = numerator of the Jack polynomial alpha coefficient
// architecture = architecture to use for precalculation
// symbolicDepth = use symbolic calculation to solve singular values if non zero, if greater than zero it will use symbolic calculation up to a given depth (below that depth it will rely on numerical calculation),
//                 -1 if the symbolic calculation has to be done up to the point the singular value problem has been solved
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// fileName = optional file name to store temporary calculations
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

LongRationalVector& FermionOnSphereHaldaneBasisLong::GenerateJackPolynomial(LongRationalVector& jack, long alphaNumerator, long alphaDenominator, AbstractArchitecture* architecture, int symbolicDepth, long minIndex, long maxIndex, char* fileName)
{
  jack[0] = 1l;
  LongRational InvAlpha (2l * (alphaDenominator - alphaNumerator), alphaNumerator);

  int ReducedNbrFermions = this->NbrFermions - 1;
  unsigned long* TmpMonomial = new unsigned long [this->NbrFermions];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrFermions];
  long* ConnectedIndices = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients  = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];
  long* ConnectedIndices2 = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients2  = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];

  long CounterMask = 0xffffl;
  if (symbolicDepth > 1)
    {
      CounterMask = 0xffl;
    }

  LongRational RhoRoot = 0l;
  LongRational Rho = 0l;
  ULONGLONG MaxRoot = this->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, TmpMonomial);
  long TmpRhoRootInvAlphaCoef = 0l;
  long TmpRhoRootConstCoef  = 0l;
  for (int j = 0; j < this->NbrFermions; ++j)
    {
      TmpRhoRootInvAlphaCoef -= TmpMonomial[j] * ((long) j);
      TmpRhoRootConstCoef += TmpMonomial[j] * (TmpMonomial[j] - 1l);
      RhoRoot += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
    }
  LongRational RhoRootInvAlphaCoef(TmpRhoRootInvAlphaCoef);
  LongRational RhoRootConstCoef(TmpRhoRootConstCoef);

  LongRationalPolynomial* TmpNumerators = 0;
  LongRationalPolynomial* TmpDenominators = 0;
  int* EvaluatedCoeffcients = 0;
  if (symbolicDepth != 0)
    {
      TmpNumerators = new LongRationalPolynomial[this->LargeHilbertSpaceDimension];
      TmpDenominators = new LongRationalPolynomial[this->LargeHilbertSpaceDimension];		  
      EvaluatedCoeffcients = new int[this->LargeHilbertSpaceDimension];
      
      this->GenerateSingleJackPolynomialCoefficient(0, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2,
						    RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture);
    }
  LongRational Coefficient = 0l;
  LongRational Coefficient2 = 0l;
  
  if (minIndex < 1l)
    minIndex = 1l;
  if ((maxIndex <= minIndex) || (maxIndex > this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension;

  for (long i = minIndex; i < maxIndex; ++i)
    {
      if (symbolicDepth > 1)
	{
	  this->GenerateSingleJackPolynomialCoefficient(i, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2,
							RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture);
	  if (symbolicDepth == 3)
	    {
	      cout << TmpNumerators[i] << endl;
	      cout << TmpDenominators[i] << endl;
	    }
	}
      if (jack[i] == 0l)
	{
	  Rho = 0l;
	  ULONGLONG CurrentPartition = this->StateDescription[i];
	  this->ConvertToMonomial(CurrentPartition, TmpMonomial);
	  for (int j = 0; j < this->NbrFermions; ++j)
	    Rho += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
	  if (Rho == RhoRoot)
	    {
 	      if (symbolicDepth == 0)
 		{
 		  cout << "warning : singular value detected at position " << i << ", skipping the rest of the calculation" << endl;
 		  return jack;
 		}
	      if (symbolicDepth != 3)
		{
		  cout << "singular value detected at position " << i << ", using symbolic calculation" << endl;
		}
	      if (fileName != 0)
		jack.WriteVector(fileName);
	      for (long j = 0l; j < this->LargeHilbertSpaceDimension; ++j)
		{
		  if (TmpNumerators[j].Defined())
		    EvaluatedCoeffcients[j] = 1;
		  else
		    EvaluatedCoeffcients[j] = 0;
		}
	      long NbrComponentToEvaluate = this->GenerateSingleJackPolynomialCoefficientCountOnly(i, EvaluatedCoeffcients, TmpMonomial, TmpMonomial2, MaxRoot);
	      cout << "number of components to evalute using symbolic calculation : " << NbrComponentToEvaluate << endl;
	      this->GenerateSingleJackPolynomialCoefficient(i, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2,
							    RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture);
	      TmpNumerators[i].PolynomialEvaluate(InvAlpha, Coefficient);
	      Coefficient /= TmpDenominators[i].PolynomialEvaluate(InvAlpha);
	      if (symbolicDepth != 3)
		{
		  cout << "----------------------------------------------------------------" << endl
		       << "result = " << endl;
		  cout << TmpNumerators[i] << endl;
		  cout << TmpDenominators[i] << endl;
		  cout << Coefficient << endl;
		}		  
	      jack[i] = Coefficient;
	    }
	  else
	    {
	      Coefficient = 0l;
	      int Pos = 0;
	      for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
		for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
		  {
		    long Diff = (long) (TmpMonomial[j1] - TmpMonomial[j2]);
		    unsigned int Max = TmpMonomial[j2];
		    ULONGLONG TmpState = 0x0ul;
		    int Tmpj1 = j1;
		    int Tmpj2 = j2;
		    for (int l = 0; l < this->NbrFermions; ++l)
		      TmpMonomial2[l] = TmpMonomial[l];	    
		    long Sign = 1l;
		    for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		      {
			++TmpMonomial2[Tmpj1];
			--TmpMonomial2[Tmpj2];
			while ((Tmpj1 > 0) && (TmpMonomial2[Tmpj1] > TmpMonomial2[Tmpj1 - 1]))
			  {
			    unsigned long Tmp = TmpMonomial2[Tmpj1 - 1];
			    TmpMonomial2[Tmpj1 - 1] = TmpMonomial2[Tmpj1];
			    TmpMonomial2[Tmpj1] = Tmp;
			    --Tmpj1;
			    Sign *= -1l; 
			  }
			while ((Tmpj2 < ReducedNbrFermions) && (TmpMonomial2[Tmpj2] < TmpMonomial2[Tmpj2 + 1]))
			  {
			    unsigned long Tmp = TmpMonomial2[Tmpj2 + 1];
			    TmpMonomial2[Tmpj2 + 1] = TmpMonomial2[Tmpj2];
			    TmpMonomial2[Tmpj2] = Tmp;
			    ++Tmpj2;
			    Sign *= -1l; 
			  }
			if ((TmpMonomial2[Tmpj1] != TmpMonomial2[Tmpj1 + 1]) && (TmpMonomial2[Tmpj2] != TmpMonomial2[Tmpj2 - 1]))
			  {
			    TmpState = this->ConvertFromMonomial(TmpMonomial2);
			    if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
			      {
				long TmpIndex = this->FindStateIndex(TmpState, TmpMonomial2[0]);
				if (TmpIndex < this->HilbertSpaceDimension)
				  {
				    ConnectedIndices[Pos] = TmpIndex;
				    ConnectedCoefficients[Pos] = Sign * Diff;
				    ++Pos;
				  }
 			      }
			  }
		      }
		  }
	      int NbrConnected = 1l;
	      if (Pos > 1)
		{
		  SortArrayDownOrdering<long>(ConnectedIndices, ConnectedCoefficients, Pos);
		  int TmpIndex = 1;
		  while (TmpIndex < Pos)
		    {
		      while ((TmpIndex < Pos) && (ConnectedIndices[TmpIndex] == ConnectedIndices[TmpIndex - 1]))
			++TmpIndex;
		      if (TmpIndex < Pos)
			++NbrConnected;
		      ++TmpIndex;
		    }
		  ConnectedIndices2[0] = ConnectedIndices[0];
		  ConnectedCoefficients2[0] = ConnectedCoefficients[0];
		  TmpIndex = 1;
		  NbrConnected = 1;
		  while (TmpIndex < Pos)
		    {
		      while ((TmpIndex < Pos) && (ConnectedIndices[TmpIndex] == ConnectedIndices[TmpIndex - 1]))
			{
			  ConnectedCoefficients2[NbrConnected - 1] += ConnectedCoefficients[TmpIndex];
			  ++TmpIndex;
			}
		      if (TmpIndex < Pos)
			{
			  ConnectedIndices2[NbrConnected] = ConnectedIndices[TmpIndex];
			  ConnectedCoefficients2[NbrConnected] = ConnectedCoefficients[TmpIndex];	   
			  ++NbrConnected;
			}
		      ++TmpIndex;
		    }
		}
	      else
		{
		  ConnectedIndices2[0] = ConnectedIndices[0];
		  ConnectedCoefficients2[0] = ConnectedCoefficients[0];
		}
 	      Coefficient = ConnectedCoefficients2[0];	  
	      Coefficient *= jack[ConnectedIndices2[0]];
	      for (int j = 1; j < NbrConnected; ++j)
		{
		  Coefficient2 = ConnectedCoefficients2[j];
		  Coefficient2 *= jack[ConnectedIndices2[j]];
		  Coefficient += Coefficient2;
		}
	      Coefficient *= InvAlpha;
	      Rho -= RhoRoot;
	      Rho.Neg();
	      Coefficient /= Rho;
	      jack[i] = Coefficient;
	    }
	}
      if ((i & CounterMask) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  delete[] TmpMonomial;
  delete[] TmpMonomial2;
  delete[] ConnectedIndices;
  delete[] ConnectedCoefficients;
  delete[] ConnectedIndices2;
  delete[] ConnectedCoefficients2;
  cout << endl;
  if (symbolicDepth != 0)
    {
      delete[] TmpNumerators;
      delete[] TmpDenominators;		  
      delete[] EvaluatedCoeffcients;
    }
  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur and the resulting state is invariant under the Lz<->-Lz symmetry
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alphaNumerator = numerator of the Jack polynomial alpha coefficient
// alphaDenominator = numerator of the Jack polynomial alpha coefficient
// symbolicDepth = use symbolic calculation to solve singular values if non zero, if greater than zero it will use symbolic calculation up to a given depth (below that depth it will rely on numerical calculation),
//                 -1 if the symbolic calculation has to be done up to the point the singular value problem has been solved
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// fileName = optional file name to store temporary calculations
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

LongRationalVector& FermionOnSphereHaldaneBasisLong::GenerateSymmetrizedJackPolynomial(LongRationalVector& jack, long alphaNumerator, long alphaDenominator, AbstractArchitecture* architecture, int symbolicDepth, long minIndex, long maxIndex, char* fileName)
{
  jack[0] = 1l;
  LongRational InvAlpha (2l * (alphaDenominator - alphaNumerator), alphaNumerator);

  int ReducedNbrFermions = this->NbrFermions - 1;
  unsigned long* TmpMonomial = new unsigned long [this->NbrFermions];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrFermions];
  long* ConnectedIndices = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients  = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];
  long* ConnectedIndices2 = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients2  = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];

  long CounterMask = 0xffffl;
  if (symbolicDepth > 1)
    {
      CounterMask = 0xffl;
    }

  LongRational RhoRoot = 0l;
  LongRational Rho = 0l;
  ULONGLONG MaxRoot = this->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, TmpMonomial);
  long TmpRhoRootInvAlphaCoef = 0l;
  long TmpRhoRootConstCoef  = 0l;
  for (int j = 0; j < this->NbrFermions; ++j)
    {
      TmpRhoRootInvAlphaCoef -= TmpMonomial[j] * ((long) j);
      TmpRhoRootConstCoef += TmpMonomial[j] * (TmpMonomial[j] - 1l);
      RhoRoot += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
    }
  LongRational RhoRootInvAlphaCoef(TmpRhoRootInvAlphaCoef);
  LongRational RhoRootConstCoef(TmpRhoRootConstCoef);
  
  LongRationalPolynomial* TmpNumerators = 0;
  LongRationalPolynomial* TmpDenominators = 0;
  int* EvaluatedCoeffcients = 0 ;
  if (symbolicDepth != 0)
    {
      TmpNumerators = new LongRationalPolynomial[this->LargeHilbertSpaceDimension];
      TmpDenominators = new LongRationalPolynomial[this->LargeHilbertSpaceDimension];		  
      EvaluatedCoeffcients = new int[this->LargeHilbertSpaceDimension];

      this->GenerateSingleJackPolynomialCoefficient(0, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2, 
						    RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture);
    }
  LongRational Coefficient = 0l;
  LongRational Coefficient2 = 0l;
  long SymSign = 1l;
  if ((((this->NbrFermions * ReducedNbrFermions) >> 1) & 1) != 0)
    SymSign = -1l;

  if (minIndex < 1l)
    minIndex = 1l;
  if ((maxIndex <= minIndex) || (maxIndex > this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension;
  
  for (long i = minIndex; i < maxIndex; ++i)
    {
      if (jack[i] == 0l)
	{
	  Rho = 0l;
	  ULONGLONG CurrentPartition = this->StateDescription[i];
	  this->ConvertToMonomial(CurrentPartition, TmpMonomial);
	  for (int j = 0; j < this->NbrFermions; ++j)
	    Rho += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
	  if (symbolicDepth > 1)
	    {
	      this->GenerateSingleJackPolynomialCoefficient(i, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2, 
							    RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture, 0, true);
	      long TmpIndex = this->FindStateIndex(this->GetSymmetricState(CurrentPartition), this->LzMax - TmpMonomial[ReducedNbrFermions]);
	      if (i < TmpIndex)
		{
		  TmpNumerators[TmpIndex] = TmpNumerators[i];
		  TmpDenominators[TmpIndex] = TmpDenominators[i];
		}
	      if (symbolicDepth == 3)
		{
		  cout << TmpNumerators[i] << endl;
		  cout << TmpDenominators[i] << endl;
		}
	    }
	  if (Rho == RhoRoot)
	    {
 	      if (symbolicDepth == 0)
 		{
 		  cout << "warning : singular value detected at position " << i << ", skipping the rest of the calculation" << endl;
 		  return jack;
 		}
	      if (symbolicDepth != 3)
		{
		  cout << "singular value detected at position " << i << ", using symbolic calculation" << endl;
		}
	      if (fileName != 0)
		jack.WriteVector(fileName);
	      for (long j = 0l; j < this->LargeHilbertSpaceDimension; ++j)
		{
		  if (TmpNumerators[j].Defined())
		    EvaluatedCoeffcients[j] = 1;
		  else
		    EvaluatedCoeffcients[j] = 0;
		}
	      long NbrComponentToEvaluate = this->GenerateSingleJackPolynomialCoefficientCountOnly(i, EvaluatedCoeffcients, TmpMonomial, TmpMonomial2, MaxRoot);
	      cout << "number of components to evalute using symbolic calculation : " << NbrComponentToEvaluate << endl;
	      this->GenerateSingleJackPolynomialCoefficient(i, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2, 
							    RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture, 0, true);
	      LongRational Tmp = TmpNumerators[i].PolynomialEvaluate(InvAlpha);
	      Tmp /= TmpDenominators[i].PolynomialEvaluate(InvAlpha);
	      if (symbolicDepth != 3)
		{
		  cout << "----------------------------------------------------------------" << endl
		       << "result = " << endl;
		  cout << TmpNumerators[i] << endl;
		  cout << TmpDenominators[i] << endl;
		  cout << Tmp << endl;
		}		  
	      long TmpIndex = this->FindStateIndex(this->GetSymmetricState(CurrentPartition), this->LzMax - TmpMonomial[ReducedNbrFermions]);
	      jack[i] = Tmp;
	      if (i < TmpIndex)
		{
		  jack[TmpIndex] = Tmp;
		  TmpNumerators[TmpIndex] = TmpNumerators[i];
		  TmpDenominators[TmpIndex] = TmpDenominators[i];
		}
	    }
	  else
	    {
	      Coefficient = 0l;
	      int Pos = 0;
	      for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
		for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
		  {
		    long Diff = (long) (TmpMonomial[j1] - TmpMonomial[j2]);
		    unsigned int Max = TmpMonomial[j2];
		    ULONGLONG TmpState = 0x0ul;
		    int Tmpj1 = j1;
		    int Tmpj2 = j2;
		    for (int l = 0; l < this->NbrFermions; ++l)
		      TmpMonomial2[l] = TmpMonomial[l];	 
		    long Sign = 1l;   
		    for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		      {
			++TmpMonomial2[Tmpj1];
			--TmpMonomial2[Tmpj2];
			while ((Tmpj1 > 0) && (TmpMonomial2[Tmpj1] > TmpMonomial2[Tmpj1 - 1]))
			  {
			    unsigned long Tmp = TmpMonomial2[Tmpj1 - 1];
			    TmpMonomial2[Tmpj1 - 1] = TmpMonomial2[Tmpj1];
			    TmpMonomial2[Tmpj1] = Tmp;
			    Sign *= -1.0; 
			    --Tmpj1;
			  }
			while ((Tmpj2 < ReducedNbrFermions) && (TmpMonomial2[Tmpj2] < TmpMonomial2[Tmpj2 + 1]))
			  {
			    unsigned long Tmp = TmpMonomial2[Tmpj2 + 1];
			    TmpMonomial2[Tmpj2 + 1] = TmpMonomial2[Tmpj2];
			    TmpMonomial2[Tmpj2] = Tmp;
			    Sign *= -1.0; 
			    ++Tmpj2;
			  }
			if ((TmpMonomial2[Tmpj1] != TmpMonomial2[Tmpj1 + 1]) && (TmpMonomial2[Tmpj2] != TmpMonomial2[Tmpj2 - 1]))
			  {
			    TmpState = this->ConvertFromMonomial(TmpMonomial2);
			    if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
			      {
				long TmpIndex = this->FindStateIndex(TmpState, TmpMonomial2[0]);
				if (TmpIndex < this->HilbertSpaceDimension)
				  {
				    long TmpIndex2 = this->FindStateIndex(this->GetSymmetricState(TmpState), this->LzMax - TmpMonomial2[ReducedNbrFermions]);
				    if (TmpIndex < TmpIndex2)
				      {
					ConnectedIndices[Pos] = TmpIndex;
					ConnectedCoefficients[Pos] = Sign * Diff;
				      }
				    else
				      {
					ConnectedIndices[Pos] = TmpIndex2;
					ConnectedCoefficients[Pos] = SymSign * Sign * Diff;
				      }
				    ++Pos;
				  }
			      }
			  }
		      }
		  }
	      int NbrConnected = 1l;
	      if (Pos > 1)
		{
		  SortArrayDownOrdering<long>(ConnectedIndices, ConnectedCoefficients, Pos);
		  int TmpIndex = 1;
		  while (TmpIndex < Pos)
		    {
		      while ((TmpIndex < Pos) && (ConnectedIndices[TmpIndex] == ConnectedIndices[TmpIndex - 1]))
			++TmpIndex;
		      if (TmpIndex < Pos)
			++NbrConnected;
		      ++TmpIndex;
		    }
		  ConnectedIndices2[0] = ConnectedIndices[0];
		  ConnectedCoefficients2[0] = ConnectedCoefficients[0];
		  TmpIndex = 1;
		  NbrConnected = 1;
		  while (TmpIndex < Pos)
		    {
		      while ((TmpIndex < Pos) && (ConnectedIndices[TmpIndex] == ConnectedIndices[TmpIndex - 1]))
			{
			  ConnectedCoefficients2[NbrConnected - 1] += ConnectedCoefficients[TmpIndex];
			  ++TmpIndex;
			}
		      if (TmpIndex < Pos)
			{
			  ConnectedIndices2[NbrConnected] = ConnectedIndices[TmpIndex];
			  ConnectedCoefficients2[NbrConnected] = ConnectedCoefficients[TmpIndex];	   
			  ++NbrConnected;
			}
		      ++TmpIndex;
		    }
		}
	      else
		{
		  ConnectedIndices2[0] = ConnectedIndices[0];
		  ConnectedCoefficients2[0] = ConnectedCoefficients[0];
		}
 	      Coefficient = ConnectedCoefficients2[0];	  
	      Coefficient *= jack[ConnectedIndices2[0]];
	      for (int j = 1; j < NbrConnected; ++j)
		{
		  Coefficient2 = ConnectedCoefficients2[j];
		  Coefficient2 *= jack[ConnectedIndices2[j]];
		  Coefficient += Coefficient2;
		}
	      long TmpIndex = this->FindStateIndex(this->GetSymmetricState(CurrentPartition), this->LzMax - TmpMonomial[ReducedNbrFermions]);
	      Coefficient *= InvAlpha;
	      Rho -= RhoRoot;
	      Rho.Neg();
	      Coefficient /= Rho;
	      if (i < TmpIndex)
		jack[TmpIndex] = Coefficient;
	      jack[i] = Coefficient;
	    }
	}
      if ((i & CounterMask) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  delete[] TmpMonomial;
  delete[] TmpMonomial2;
  delete[] ConnectedIndices;
  delete[] ConnectedCoefficients;
  delete[] ConnectedIndices2;
  delete[] ConnectedCoefficients2;
  if (symbolicDepth != 0)
    {
      delete[] TmpNumerators;
      delete[] TmpDenominators;		  
      delete[] EvaluatedCoeffcients;
    }
  cout << endl;
  return jack;
}

// compute a single coefficient of the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur and using (partial symbolic calculation)
//
// index = index of the component to compute
// numerators = array of polynomials attached to each coefficient numerator
// denominators = array of polynomials attached to each coefficient denominator
// tmpMonomial = temporary array for monomial description
// tmpMonomial2 = temporary array for monomial description
// rhoRootInvAlphaCoef = coefficient in front of inv alpha in the rho for the root partition
// rhoRootConstCoef = constant coefficient in the rho for the root partition
// maxRoot = fermionic expression for the root partition
// architecture = architecture to use for precalculation
// currentNbrComponents = current number of components computed for the symbolic of the index-th compoment
// symmetryFlag = true iof the Lz <-> -Lz symmetry has to be used
// return value = total number of components computed for the symbolic of the index-th compoment

long FermionOnSphereHaldaneBasisLong::GenerateSingleJackPolynomialCoefficient(long index, LongRationalPolynomial* numerators, LongRationalPolynomial* denominators, long* connectedIndices, long* connectedCoefficients, unsigned long* tmpMonomial, unsigned long* tmpMonomial2, LongRational& rhoRootInvAlphaCoef, LongRational& rhoRootConstCoef, ULONGLONG maxRoot, AbstractArchitecture* architecture, long currentNbrComponents, bool symmetryFlag)
{
  if (numerators[index].Defined())
    return 0l;
  if (index == 0l)
    {
      if (!numerators[0l].Defined())
	{
	  numerators[0l] = LongRationalPolynomial(0);
	  denominators[0l] = LongRationalPolynomial(0);
	  numerators[0l][0] = 1l;
	  denominators[0l][0] = 1l;
	  return 0l;
	}
    }


  int ReducedNbrFermions = this->NbrFermions - 1;  

  ULONGLONG CurrentPartition = this->StateDescription[index];
  this->ConvertToMonomial(CurrentPartition, tmpMonomial);
  long TmpRhoInvAlphaCoef = 0l;
  long TmpRhoConstCoef  = 0l;
  for (int j = 0; j < this->NbrFermions; ++j)
    {
      TmpRhoInvAlphaCoef -= tmpMonomial[j] * ((long) j);
      TmpRhoConstCoef += tmpMonomial[j] * (tmpMonomial[j] - 1l);
    }
  LongRational RhoInvAlphaCoef(TmpRhoInvAlphaCoef);
  LongRational RhoConstCoef (TmpRhoConstCoef);
  double SymSign = 1l;
  if ((((this->NbrFermions * ReducedNbrFermions) >> 1) & 1) != 0)
    SymSign = -1l;
  
  int Pos = 0;
  for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
    for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
      {
	long Diff = (long) (tmpMonomial[j1] - tmpMonomial[j2]);
	unsigned int Max = tmpMonomial[j2];
	ULONGLONG TmpState = 0x0ul;
	int Tmpj1 = j1;
	int Tmpj2 = j2;
	for (int l = 0; l < this->NbrFermions; ++l)
	  tmpMonomial2[l] = tmpMonomial[l];	    
	long Sign = 1l;
	for (unsigned int k = 1; (k <= Max) && (TmpState < maxRoot); ++k)
	  {
	    ++tmpMonomial2[Tmpj1];
	    --tmpMonomial2[Tmpj2];
	    while ((Tmpj1 > 0) && (tmpMonomial2[Tmpj1] > tmpMonomial2[Tmpj1 - 1]))
	      {
		unsigned long Tmp = tmpMonomial2[Tmpj1 - 1];
		tmpMonomial2[Tmpj1 - 1] = tmpMonomial2[Tmpj1];
		tmpMonomial2[Tmpj1] = Tmp;
		--Tmpj1;
		Sign *= -1l; 
	      }
	    while ((Tmpj2 < ReducedNbrFermions) && (tmpMonomial2[Tmpj2] < tmpMonomial2[Tmpj2 + 1]))
	      {
		unsigned long Tmp = tmpMonomial2[Tmpj2 + 1];
		tmpMonomial2[Tmpj2 + 1] = tmpMonomial2[Tmpj2];
		tmpMonomial2[Tmpj2] = Tmp;
		++Tmpj2;
		Sign *= -1l; 
	      }
	    if ((tmpMonomial2[Tmpj1] != tmpMonomial2[Tmpj1 + 1]) && (tmpMonomial2[Tmpj2] != tmpMonomial2[Tmpj2 - 1]))
	      {
		TmpState = this->ConvertFromMonomial(tmpMonomial2);
		if ((TmpState <= maxRoot) && (TmpState > CurrentPartition))
		  {
		    long TmpIndex = this->FindStateIndex(TmpState, tmpMonomial2[0]);
		    if (TmpIndex < this->HilbertSpaceDimension)
		      {
			if (symmetryFlag == false)
			  {
			    connectedIndices[Pos] = TmpIndex;
			    connectedCoefficients[Pos] = Sign * Diff;
			    ++Pos;
			  }
			else
			  {
			    long TmpIndex2 = this->FindStateIndex(this->GetSymmetricState(TmpState), this->LzMax - tmpMonomial2[ReducedNbrFermions]);
			    if (TmpIndex < TmpIndex2)
			      {
				connectedIndices[Pos] = TmpIndex;
				connectedCoefficients[Pos] = Sign * Diff;
			      }
			    else
			      {
				connectedIndices[Pos] = TmpIndex2;
				connectedCoefficients[Pos] = SymSign * Sign * Diff;
			      }
			    ++Pos;
			  }
		      }
		  }
	      }
	  }
      }
  long* ConnectedIndices2 = 0;
  long* ConnectedCoefficients2  = 0;
  int NbrConnected = 1l;
  if (Pos > 1)
    {
      SortArrayDownOrdering<long>(connectedIndices, connectedCoefficients, Pos);
      int TmpIndex = 1;
      while (TmpIndex < Pos)
	{
	  while ((TmpIndex < Pos) && (connectedIndices[TmpIndex] == connectedIndices[TmpIndex - 1]))
	    ++TmpIndex;
	  if (TmpIndex < Pos)
	    ++NbrConnected;
	  ++TmpIndex;
	}
      ConnectedIndices2 = new long[NbrConnected];
      ConnectedCoefficients2  = new long[NbrConnected];
      ConnectedIndices2[0] = connectedIndices[0];
      ConnectedCoefficients2[0] = connectedCoefficients[0];
      TmpIndex = 1;
      NbrConnected = 1;
      while (TmpIndex < Pos)
	{
	  while ((TmpIndex < Pos) && (connectedIndices[TmpIndex] == connectedIndices[TmpIndex - 1]))
	    {
	      ConnectedCoefficients2[NbrConnected - 1] += connectedCoefficients[TmpIndex];
	      ++TmpIndex;
	    }
	  if (TmpIndex < Pos)
	    {
	      ConnectedIndices2[NbrConnected] = connectedIndices[TmpIndex];
	      ConnectedCoefficients2[NbrConnected] = connectedCoefficients[TmpIndex];	   
	      ++NbrConnected;
	    }
	  ++TmpIndex;
	}
    }
  else
    {
      ConnectedIndices2 = new long[NbrConnected];
      ConnectedCoefficients2  = new long[NbrConnected];
      ConnectedIndices2[0] = connectedIndices[0];
      ConnectedCoefficients2[0] = connectedCoefficients[0];
    }

  int TmpNbrComponents = 0;
  for (int i = 0; i < NbrConnected; ++i)
    {
      if (!numerators[ConnectedIndices2[i]].Defined())
        {
	  TmpNbrComponents += this->GenerateSingleJackPolynomialCoefficient(ConnectedIndices2[i], numerators, denominators, connectedIndices, connectedCoefficients, 
									    tmpMonomial, tmpMonomial2, rhoRootInvAlphaCoef, rhoRootConstCoef, maxRoot, architecture, currentNbrComponents + TmpNbrComponents, symmetryFlag);
	}
    }

  cout << (currentNbrComponents + TmpNbrComponents) << " symbolic components computed                  \r";
  cout.flush();

  LongRational Tmp4;
  if (NbrConnected > 1)
    {
      FQHESphereJackGeneratorSumRationalPolynomialOperation Operation(index, numerators, denominators, ConnectedIndices2, ConnectedCoefficients2, NbrConnected);
      Operation.ApplyOperation(architecture);
    }
  else
    {
      numerators[index]= numerators[ConnectedIndices2[0]];
      denominators[index] = denominators[ConnectedIndices2[0]];
      numerators[index] *= ConnectedCoefficients2[0];
    }
  numerators[index].ShiftPowers(1);

  LongRational Tmp2 = rhoRootInvAlphaCoef;
  Tmp2 -= RhoInvAlphaCoef;
  if (Tmp2 != 0l)
    {
      numerators[index] /= Tmp2;
      LongRational Tmp3 = rhoRootConstCoef;
      Tmp3 -= RhoConstCoef;
      Tmp3 /= -Tmp2;
      numerators[index].PolynomialEvaluate(Tmp3, Tmp4);
      if (Tmp4 == 0l)
	{
	  numerators[index].LocalMonomialDivision(Tmp3);
	}
      else
	{
	  denominators[index].LocalMonomialMultiplication(Tmp3);
	}
    }
  else
    {
      Tmp2 = rhoRootConstCoef;
      Tmp2 -= RhoConstCoef;
      if (Tmp2 != 0)
	{
	  numerators[index] /= Tmp2;
	}
    }

  delete[] ConnectedIndices2;
  delete[] ConnectedCoefficients2;

  return (TmpNbrComponents + 1);
}


// compute of many coefficients have to be computed to get a single coefficient of the Jack polynomial decomposition 
//
// index = index of the component to compute
// evaluatedCoeffcients = that indicates which coefficients have already been computed
// tmpMonomial = temporary array for monomial description
// tmpMonomial2 = temporary array for monomial description
// maxRoot = fermionic expression for the root partition
// return value = true if a fully symbolic calculation has been performed

long FermionOnSphereHaldaneBasisLong::GenerateSingleJackPolynomialCoefficientCountOnly(long index, int* evaluatedCoeffcients, unsigned long* tmpMonomial, unsigned long* tmpMonomial2, ULONGLONG maxRoot)
{
  if (index == 0l)
    return 1l;
  int ReducedNbrFermions = this->NbrFermions - 1;  

  ULONGLONG CurrentPartition = this->StateDescription[index];
  this->ConvertToMonomial(CurrentPartition, tmpMonomial);

  long* ConnectedIndices2 = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];
  int NbrConnected  = 0;
  for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
    for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
      {
	unsigned int Max = tmpMonomial[j2];
	ULONGLONG TmpState = 0x0ul;
	int Tmpj1 = j1;
	int Tmpj2 = j2;
	for (int l = 0; l < this->NbrFermions; ++l)
	  tmpMonomial2[l] = tmpMonomial[l];	    
	for (unsigned int k = 1; (k <= Max) && (TmpState < maxRoot); ++k)
	  {
	    ++tmpMonomial2[Tmpj1];
	    --tmpMonomial2[Tmpj2];
	    while ((Tmpj1 > 0) && (tmpMonomial2[Tmpj1] > tmpMonomial2[Tmpj1 - 1]))
	      {
		unsigned long Tmp = tmpMonomial2[Tmpj1 - 1];
		tmpMonomial2[Tmpj1 - 1] = tmpMonomial2[Tmpj1];
		tmpMonomial2[Tmpj1] = Tmp;
		--Tmpj1;
	      }
	    while ((Tmpj2 < ReducedNbrFermions) && (tmpMonomial2[Tmpj2] < tmpMonomial2[Tmpj2 + 1]))
	      {
		unsigned long Tmp = tmpMonomial2[Tmpj2 + 1];
		tmpMonomial2[Tmpj2 + 1] = tmpMonomial2[Tmpj2];
		tmpMonomial2[Tmpj2] = Tmp;
		++Tmpj2;
	      }
	    if ((tmpMonomial2[Tmpj1] != tmpMonomial2[Tmpj1 + 1]) && (tmpMonomial2[Tmpj2] != tmpMonomial2[Tmpj2 - 1]))
	      {
		TmpState = this->ConvertFromMonomial(tmpMonomial2);
		if ((TmpState <= maxRoot) && (TmpState > CurrentPartition))
		  {
		    long TmpIndex = this->FindStateIndex(TmpState, tmpMonomial2[0]);
		    if (TmpIndex < this->HilbertSpaceDimension)
		      {
			ConnectedIndices2[NbrConnected] = TmpIndex;
			++NbrConnected;
		      }
		  }
	      }
	  }
      }

  long TmpNbrConnected = 1;
  for (int i = 0; i < NbrConnected; ++i)
    {
      if (evaluatedCoeffcients[ConnectedIndices2[i]] == 0)
        {
 	  TmpNbrConnected += this->GenerateSingleJackPolynomialCoefficientCountOnly(ConnectedIndices2[i], evaluatedCoeffcients, 
 										    tmpMonomial, tmpMonomial2, maxRoot);
	  evaluatedCoeffcients[ConnectedIndices2[i]] = 1;
	}
    }

  delete[] ConnectedIndices2;

  return TmpNbrConnected;
}


