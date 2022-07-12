////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//                           with Lz <-> -Lz symmetry                         //
//   that allow LzMax up to 127 (for systems with 128 bit integer support)    //
//               or 63 (on 32 bit systems without 128 bit integer support)    //
//                                                                            //
//                        last modification : 23/09/2007                      //
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
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/Endian.h"

#include <math.h>
#include <fstream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// basic constructor
// 
// nbrFermions = number of fermions
// lzMax = twice the maximum Lz value reached by a fermion
// memory = amount of memory granted for precalculations
// referenceState = array that describes the reference state to start from

FermionOnSphereHaldaneSymmetricBasisLong::FermionOnSphereHaldaneSymmetricBasisLong (int nbrFermions, int lzMax, int* referenceState,
										    unsigned long memory)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
#ifdef __128_BIT_LONGLONG__
  this->InvertShift = 64 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 32 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;
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
  this->TotalLz = ((this->TotalLz << 1) - (this->LzMax * this->NbrFermions)) >> 1;
  this->TotalLz = 0;

  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  this->Flag.Initialize();

  this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
#ifdef  __64_BITS__
  int ReducedHilbertSpaceDimension = (this->HilbertSpaceDimension >> 6) + 1;
#else
  int ReducedHilbertSpaceDimension = (this->HilbertSpaceDimension >> 5) + 1;
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
  NewHilbertSpaceDimension = 0;
#ifdef  __64_BITS__
  if ((this->HilbertSpaceDimension & 0x3f) == 0)
#else
  if ((this->HilbertSpaceDimension & 0x1f) == 0)
#endif
    --ReducedHilbertSpaceDimension;
  int TotalIndex = 0;
  for (int i = 0; i < ReducedHilbertSpaceDimension; ++i)
    {
      TmpKeepStateFlag = this->KeepStateFlag[i];
#ifdef  __64_BITS__
      for (int j = 0; j < 64; ++j)
#else
      for (int j = 0; j < 32; ++j)
#endif
	{
	  if ((TmpKeepStateFlag >> j) & 0x1ul)
	    {
	      if (this->StateDescription[TotalIndex] == this->GetCanonicalState(this->StateDescription[TotalIndex]))
		++NewHilbertSpaceDimension;
	      else
		this->KeepStateFlag[i] &= ~(0x1ul << j);
	    }
	  ++TotalIndex;
	}
    }
  ULONGLONG* TmpStateDescription = new ULONGLONG [NewHilbertSpaceDimension];
  int* TmpStateLzMax = new int [NewHilbertSpaceDimension];
  NewHilbertSpaceDimension = 0;
  TotalIndex = 0;
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
  
  delete[] this->StateDescription;
  delete[] this->StateLzMax;
  delete[] this->KeepStateFlag;
  this->StateDescription = TmpStateDescription;
  this->StateLzMax = TmpStateLzMax;
  this->HilbertSpaceDimension = NewHilbertSpaceDimension;

  delete[] this->TmpGeneratedStates;
  delete[] this->TmpGeneratedStatesLzMax;

  this->GenerateLookUpTable(memory);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    this->GetStateSymmetry(this->StateDescription[i]);
  delete[] this->StateLzMax;
  this->StateLzMax = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
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

FermionOnSphereHaldaneSymmetricBasisLong::FermionOnSphereHaldaneSymmetricBasisLong (char* fileName, unsigned long memory)
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
  ReadLittleEndian(File, this->ReferenceState);
  this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    ReadLittleEndian(File, this->StateDescription[i]);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    ReadLittleEndian(File, this->StateLzMax[i]);

  File.close();

#ifdef __128_BIT_LONGLONG__
  this->InvertShift = 64 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 32 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;
  this->TargetSpace = this;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->Flag.Initialize();

  this->StateLzMax = new int [this->HilbertSpaceDimension];
  int TmpLzMax = this->LzMax;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->StateDescription[i] &= FERMION_SPHERE_LONG_SYMMETRIC_MASK;
      while ((this->StateDescription[i] >> TmpLzMax) == ((ULONGLONG) 0))
	--TmpLzMax;
      this->StateLzMax[i] = TmpLzMax;
    }
   this->GenerateLookUpTable(memory);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    this->GetStateSymmetry(this->StateDescription[i]);
  delete[] this->StateLzMax;
  this->StateLzMax = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  unsigned long UsedMemory = 0l;
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * sizeof(ULONGLONG);
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

FermionOnSphereHaldaneSymmetricBasisLong::FermionOnSphereHaldaneSymmetricBasisLong(const FermionOnSphereHaldaneSymmetricBasisLong& fermions)
{
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->ReferenceState = fermions.ReferenceState;
  this->Flag = fermions.Flag;
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

FermionOnSphereHaldaneSymmetricBasisLong::~FermionOnSphereHaldaneSymmetricBasisLong ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereHaldaneSymmetricBasisLong& FermionOnSphereHaldaneSymmetricBasisLong::operator = (const FermionOnSphereHaldaneSymmetricBasisLong& fermions)
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereHaldaneSymmetricBasisLong::Clone()
{
  return new FermionOnSphereHaldaneSymmetricBasisLong(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphereHaldaneSymmetricBasisLong::WriteHilbertSpace (char* fileName)
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
  WriteLittleEndian(File, this->ReferenceState);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    WriteLittleEndian(File, this->StateDescription[i]);
  File.close();
  return true;
}

// convert a gien state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereHaldaneSymmetricBasisLong::ConvertToNbodyBasis(RealVector& state, FermionOnSphereLong& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  ULONGLONG TmpState;
  ULONGLONG Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      TmpState = this->GetSignedCanonicalState(nbodyBasis.StateDescription[i]);
      NewLzMax = this->LzMax;
      Signature = TmpState & FERMION_SPHERE_LONG_SYMMETRIC_BIT;
      TmpState &= FERMION_SPHERE_LONG_SYMMETRIC_MASK;
      while ((TmpState >> NewLzMax) == ((ULONGLONG) 0x0ul))
	--NewLzMax;
      if (Signature != ((ULONGLONG) 0x0ul))	
	TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2;
      else
	TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)];
    }
  return TmpVector;
}

// convert a gien state from Haldane basis to the usual symmetric n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the  symmetric nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereHaldaneSymmetricBasisLong::ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereSymmetricBasisLong& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      int NewLzMax = this->LzMax;
      ULONGLONG TmpState = this->StateDescription[i] & FERMION_SPHERE_LONG_SYMMETRIC_MASK;
      while ((TmpState >> NewLzMax) == ((ULONGLONG) 0x0ul))
	--NewLzMax;
      TmpVector[nbodyBasis.FindStateIndex(this->StateDescription[i], NewLzMax)] = state[i];
    }
  return TmpVector;
}



// convert a gien state from Lz-symmetric Haldane basis to the usual Haldane n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereHaldaneSymmetricBasisLong::ConvertToHaldaneNbodyBasis(RealVector& state, FermionOnSphereHaldaneBasisLong& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  ULONGLONG TmpState;
  ULONGLONG Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      TmpState = this->GetSignedCanonicalState(nbodyBasis.StateDescription[i]);
      NewLzMax = this->LzMax;
      Signature = TmpState & FERMION_SPHERE_LONG_SYMMETRIC_BIT;
      TmpState &= FERMION_SPHERE_LONG_SYMMETRIC_MASK;
      while ((TmpState >> NewLzMax) == ((ULONGLONG) 0x0ul))
	--NewLzMax;
      if (Signature != 0x0ul)	
	TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2;
      else
	TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)];
    }
  return TmpVector;
}

// convert a given state from the usual n-body Haldane basis to the Lz-symmetric Haldane basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereHaldaneSymmetricBasisLong::ConvertToSymmetricHaldaneNbodyBasis(RealVector& state, FermionOnSphereLong& nbodyBasis)
{
  RealVector TmpVector (this->GetHilbertSpaceDimension(), true);
  ULONGLONG TmpState;
  ULONGLONG Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      TmpState = this->GetSignedCanonicalState(nbodyBasis.StateDescription[i]);
      NewLzMax = this->LzMax;
      Signature = TmpState & FERMION_SPHERE_LONG_SYMMETRIC_BIT;
      TmpState &= FERMION_SPHERE_LONG_SYMMETRIC_MASK;
      while ((TmpState >> NewLzMax) == ((ULONGLONG) 0x0ul))
	--NewLzMax;
      if (Signature != ((ULONGLONG) 0x0ul))	
	TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += state[i] * M_SQRT1_2;
      else
	TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = state[i];
    }
  return TmpVector;  
}

// generate all states corresponding to the constraints
// 
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int FermionOnSphereHaldaneSymmetricBasisLong::GenerateStates(int lzMax, ULONGLONG referenceState, int pos, long& memory)
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

int FermionOnSphereHaldaneSymmetricBasisLong::RawGenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int pos)
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

