////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//        that has to handle (temporary) Hilbert space larger than 2^31       //
//                                                                            //
//                        last modification : 14/05/2009                      //
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
#include "HilbertSpace/FermionOnSphereHaldaneLargeBasis.h"
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


// default constructor
//

FermionOnSphereHaldaneLargeBasis::FermionOnSphereHaldaneLargeBasis()
{
  this->HilbertSpaceDimension = 0;
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = reference on twice the momentum total value, totalLz will be recomputed from referenceState and stored in totalLz
// lzMax = twice the maximum Lz value reached by a fermion
// memory = amount of memory granted for precalculations
// referenceState = array that describes the reference state to start from

FermionOnSphereHaldaneLargeBasis::FermionOnSphereHaldaneLargeBasis (int nbrFermions, int& totalLz, int lzMax, int* referenceState,
								    unsigned long memory)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->ReferenceState = 0x0ul;
  int ReferenceStateLzMax = 0;  
  this->StateLzMax = 0;
  this->LookUpTableShift = 0;
  this->LookUpTable = 0;
  this->TotalLz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      this->ReferenceState |= ((unsigned long) (referenceState[i] & 1)) << i;
      if (referenceState[i] == 1)
	{
	  ReferenceStateLzMax = i;
	  this->TotalLz += i;
	}
    }
  this->TotalLz = ((this->TotalLz << 1) - (this->LzMax * this->NbrFermions));
  totalLz = this->TotalLz;

#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->LzMax + 1 ) >> 1);
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

  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
#ifdef  __64_BITS__
  long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 6) + 1;
#else
  long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 5) + 1;
#endif
  this->KeepStateFlag = new unsigned long [ReducedHilbertSpaceDimension];
  for (int i = 0; i < ReducedHilbertSpaceDimension; ++i)
    this->KeepStateFlag[i] = 0x0l;
  this->RawGenerateStates(this->NbrFermions, this->LzMax, this->LzMax, (this->TotalLz + this->NbrFermions * this->LzMax) >> 1, 0);
  this->GenerateLookUpTable(memory);
  cout << "B-Tree look-up table depth = " << this->BTreeLookUpTableDepth << endl;
//   for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//     {
//       if (i != this->FindStateIndex(this->StateDescription[i]))
// 	  cout << "Error : " << i << " " << this->FindStateIndex(this->StateDescription[i]) << endl;
//     }
//    for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//      {
//        cout << i << " : " << StateDescription[i] << endl;
//      }
  
//   int truc = 0;
//   for (int i = 0; i <= this->BTreeLookUpTableDepth; ++i)
//     {
//       cout << "level " << i << endl;
//       for (int j = 0; j < (1 << i); ++j)
// 	{
// 	  cout << "  " << j << " -> " << this->BTreeLookUpTable[truc] << endl;
// 	  ++truc;
// 	}
//     }
//   cout << "hil : " << endl;
//   for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//     {
//       cout << this->StateDescription[i] << endl;
//     }

  int MaxSweeps = (this->NbrFermions * (this->NbrFermions - 1)) >> 1;  
  this->TmpGeneratedStates =  new unsigned long [MaxSweeps * 1000];
  this->TmpGeneratedStatesLzMax = new int [MaxSweeps * 1000];
  long Memory = 0l;

  long TmpIndex = this->FindStateIndex(this->ReferenceState);
#ifdef  __64_BITS__
  this->KeepStateFlag[TmpIndex >> 6] = 0x1l << (TmpIndex & 0x3f);
#else
  this->KeepStateFlag[TmpIndex >> 5] = 0x1l << (TmpIndex & 0x1f);
#endif
  this->GenerateStates(ReferenceStateLzMax, this->ReferenceState, 1, Memory);  

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

  unsigned long* TmpStateDescription = new unsigned long [NewHilbertSpaceDimension];
  NewHilbertSpaceDimension = 0l;
  long TotalIndex = 0;
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
  if (this->LargeHilbertSpaceDimension >= (1l << 31))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

  delete[] this->TmpGeneratedStates;
  delete[] this->TmpGeneratedStatesLzMax;

  this->GenerateLookUpTable(memory);

#ifdef __DEBUG__
  unsigned long UsedMemory = 0l;
  UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * sizeof(unsigned long);
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += (1l << (this->BTreeLookUpTableDepth + 1)) * sizeof(unsigned long);
  UsedMemory +=  (1l << this->MaximumSignLookUp) * sizeof(double);
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

FermionOnSphereHaldaneLargeBasis::FermionOnSphereHaldaneLargeBasis (char* fileName, unsigned long memory)
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
  FileSize -= sizeof (int);
  ReadLittleEndian(File, this->LargeHilbertSpaceDimension);
  FileSize -= sizeof (unsigned long);
  cout << this->LargeHilbertSpaceDimension << endl;
  ReadLittleEndian(File, this->NbrFermions);
  FileSize -= sizeof (int);
  ReadLittleEndian(File, this->LzMax);
  FileSize -= sizeof (int);
  ReadLittleEndian(File, this->TotalLz);
  FileSize -= sizeof (int);
  ReadLittleEndian(File, this->ReferenceState);
  FileSize -= sizeof (unsigned long);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    ReadLittleEndian(File, this->StateDescription[i]);
  File.close();

#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->LzMax + 1 ) >> 1);
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

  this->GenerateLookUpTable(memory);
  this->LookUpTableShift = 0;
  this->LookUpTable = 0;
#ifdef __DEBUG__
  unsigned long UsedMemory = 0l;
  UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(unsigned long));
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += (1l << (this->BTreeLookUpTableDepth + 1)) * sizeof(unsigned long);
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

FermionOnSphereHaldaneLargeBasis::FermionOnSphereHaldaneLargeBasis(const FermionOnSphereHaldaneLargeBasis& fermions)
{
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = 0;
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
}

// destructor
//

FermionOnSphereHaldaneLargeBasis::~FermionOnSphereHaldaneLargeBasis ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereHaldaneLargeBasis& FermionOnSphereHaldaneLargeBasis::operator = (const FermionOnSphereHaldaneLargeBasis& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
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
  this->InitializeWaveFunctionEvaluation();
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereHaldaneLargeBasis::Clone()
{
  return new FermionOnSphereHaldaneLargeBasis(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphereHaldaneLargeBasis::WriteHilbertSpace (char* fileName)
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
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    WriteLittleEndian(File, this->StateDescription[i]);
  File.close();
  return true;
}

// convert a gien state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereHaldaneLargeBasis::ConvertToNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  int TmpLzMax = this->LzMax;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      while ((TmpState>> TmpLzMax) != 0ul)
	--TmpLzMax;
      TmpVector[nbodyBasis.FindStateIndex(TmpState, TmpLzMax)] = state[i];
    }
  return TmpVector;
}

// convert a given state from the usual n-body basis to the Haldane basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereHaldaneLargeBasis::ConvertFromNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis)
{
  RealVector TmpVector (this->HilbertSpaceDimension, true);
  int TmpLzMax = this->LzMax;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      while ((TmpState>> TmpLzMax) != 0ul)
	--TmpLzMax;
      TmpVector[i] = state[nbodyBasis.FindStateIndex(TmpState, TmpLzMax)];
    }
  TmpVector /= TmpVector.Norm();
  return TmpVector;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// return value = corresponding index

long FermionOnSphereHaldaneLargeBasis::FindStateIndex(unsigned long stateDescription)
{
  long PosMax = 0l;
  long PosMin = 0l;
  for (int i = 0; i <= this->BTreeLookUpTableDepth; ++i)
    {
      if (stateDescription <= this->BTreeLookUpTable[PosMin + PosMax])
	{
	  PosMax <<= 1;
	  ++PosMax;
	}
      else
	PosMax <<= 1;
      PosMin += 1l << i;      
    }
  PosMin = PosMax;
  ++PosMin;
  PosMax *= this->BTreeLookUpTableNbrIndices;
  if (PosMin == (1l << (this->BTreeLookUpTableDepth + 1)))
    PosMin = this->LargeHilbertSpaceDimension - 1l;
  else
    PosMin *= this->BTreeLookUpTableNbrIndices;

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

// generate all states (i.e. all possible skew symmetric polynomials with fixed Lz)
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// currentLzMax = momentum maximum value for fermions that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereHaldaneLargeBasis::RawGenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, long pos)
{
  if ((nbrFermions == 0) || (totalLz < 0) || (currentLzMax < (nbrFermions - 1)))
    return pos;
  int LzTotalMax = ((2 * currentLzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (totalLz > LzTotalMax)
    return pos;
  if ((nbrFermions == 1) && (currentLzMax >= totalLz))
    {
      this->StateDescription[pos] = ((unsigned long) 0x1) << totalLz;
      return pos + 1l;
    }
  if (LzTotalMax == totalLz)
    {
      unsigned long Mask = 0;
      for (int i = currentLzMax - nbrFermions + 1; i <= currentLzMax; ++i)
	Mask |= (((unsigned long) 1) << i);
      this->StateDescription[pos] = Mask;
      return pos + 1l;
    }

  int ReducedCurrentLzMax = currentLzMax - 1;
  long TmpPos = this->RawGenerateStates(nbrFermions - 1, lzMax, ReducedCurrentLzMax, totalLz - currentLzMax, pos);
  unsigned long Mask = ((unsigned long) 1) << currentLzMax;
  for (long i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (lzMax == currentLzMax)
    return this->RawGenerateStates(nbrFermions, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, TmpPos);
  else
    return this->RawGenerateStates(nbrFermions, lzMax, ReducedCurrentLzMax, totalLz, TmpPos);
}

// generate all states corresponding to the constraints
// 
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereHaldaneLargeBasis::GenerateStates(int lzMax, unsigned long referenceState, long pos, long& memory)
{
  int MaxSweeps = (this->NbrFermions * (this->NbrFermions - 1)) >> 1;  
  unsigned long* TmpGeneratedStates2 = this->TmpGeneratedStates + (MaxSweeps * memory);
  int* TmpLzMax = this->TmpGeneratedStatesLzMax  + (MaxSweeps  * memory);  
  memory += 1;
  int TmpCurrentLzMax = 2;
  int TmpCurrentLzMax2;
  int TmpMax = lzMax - 1;
  int NbrEntries = 0;
  unsigned long TmpReferenceState;
  
  while (TmpCurrentLzMax < TmpMax)
    {
      while ((TmpCurrentLzMax < TmpMax) && (((referenceState >> TmpCurrentLzMax) & 0x3l) != 0x2l))
	++TmpCurrentLzMax;
      if (TmpCurrentLzMax < TmpMax)
	{
	  TmpReferenceState = (referenceState & ~(0x3l << TmpCurrentLzMax)) | (0x1l << TmpCurrentLzMax);
	  TmpCurrentLzMax2 = TmpCurrentLzMax - 2;
	  while (TmpCurrentLzMax2 >= 0)
	    {
	      while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> TmpCurrentLzMax2) & 0x3l) != 0x1l))
		--TmpCurrentLzMax2;
	      if (TmpCurrentLzMax2 >= 0)
		{
		  TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0x3l << TmpCurrentLzMax2)) | (0x2l << TmpCurrentLzMax2);
		  TmpLzMax[NbrEntries] = lzMax;
		  ++NbrEntries;
		  --TmpCurrentLzMax2;
		}	      
	    }
	  ++TmpCurrentLzMax;
	}
    }
  if (((referenceState >> TmpCurrentLzMax) & 0x3l) == 0x2l)
    {
      TmpReferenceState = (referenceState & ~(0x3l << TmpCurrentLzMax)) | (0x1l << TmpCurrentLzMax);
      TmpCurrentLzMax2 = TmpCurrentLzMax - 2;
      while (TmpCurrentLzMax2 >= 0)
	{
	  while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> TmpCurrentLzMax2) & 0x3l) != 0x1l))
	    --TmpCurrentLzMax2;
	  if (TmpCurrentLzMax2 >= 0)
	    {
	      TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0x3l << TmpCurrentLzMax2)) | (0x2l << TmpCurrentLzMax2);
	      TmpLzMax[NbrEntries] = lzMax - 1;
	      ++NbrEntries;
	      --TmpCurrentLzMax2;
	    }
	}      
    }

  unsigned long TmpIndex;
  int NbrNewEntries = 0;
  for (int i = 0; i < NbrEntries; ++i)
    {
      unsigned long& TmpState = TmpGeneratedStates2[i];
      TmpIndex = this->FindStateIndex(TmpState);
#ifdef __64_BITS__
      if ((this->KeepStateFlag[TmpIndex >> 6] >> (TmpIndex & 0x3f)) & 0x1l)
	{
	  TmpGeneratedStates2[i] = 0x0l;
	}
      else
	{
	  this->KeepStateFlag[TmpIndex >> 6] |= 0x1l << (TmpIndex & 0x3f);
	  ++NbrNewEntries;
	  if (this->TotalLz == 0)
	    {
	      unsigned long TmpSymmetricState = this->GetSymmetricState (TmpState);
	      int TmpSymLzMax = this->LzMax;
	      while (((TmpSymmetricState >> TmpSymLzMax) & 0x1ul) == 0x0ul)
		--TmpSymLzMax;
	      TmpIndex = this->FindStateIndex(TmpSymmetricState);
	      this->KeepStateFlag[TmpIndex >> 6] |= 0x1l << (TmpIndex & 0x3f);	      
	    }
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
	  if (this->TotalLz == 0)
	    {
	      unsigned long TmpSymmetricState = this->GetSymmetricState (TmpState);
	      int TmpSymLzMax = this->LzMax;
	      while (((TmpSymmetricState >> TmpSymLzMax) & 0x1ul) == 0x0ul)
		--TmpSymLzMax;
	      TmpIndex = this->FindStateIndex(TmpSymmetricState);
	      this->KeepStateFlag[TmpIndex >> 5] |= 0x1l << (TmpIndex & 0x1f);	      
	    }
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


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereHaldaneLargeBasis::GenerateLookUpTable(unsigned long memory)
{

  // build B-tree look-up table
  memory >>= 3;
  this->BTreeLookUpTableDepth = 0;
  while (((this->LargeHilbertSpaceDimension >> (this->BTreeLookUpTableDepth + 1)) != 0l) && ((1ul << this->BTreeLookUpTableDepth) < memory))
    ++this->BTreeLookUpTableDepth;
  --this->BTreeLookUpTableDepth;
  this->BTreeLookUpTable = new unsigned long [(1l << (this->BTreeLookUpTableDepth + 1)) - 1l];
  long TmpIndex = 1l;
  this->BTreeLookUpTable[0] = this->StateDescription[(this->LargeHilbertSpaceDimension >> (this->BTreeLookUpTableDepth + 1)) << this->BTreeLookUpTableDepth];
  for (int i = 1; i <= this->BTreeLookUpTableDepth; ++i)
    {
      long TmpJump = (this->LargeHilbertSpaceDimension >> (this->BTreeLookUpTableDepth + 1)) << (this->BTreeLookUpTableDepth - i + 1);
      long TmpIndex2 = TmpJump >> 1;
      for (int j = 0; j < (1 << i); ++j)
	{
	  this->BTreeLookUpTable[TmpIndex] = this->StateDescription[TmpIndex2];
	  ++TmpIndex;
	  TmpIndex2 += TmpJump;
	}
    }
  this->BTreeLookUpTableNbrIndices = this->LargeHilbertSpaceDimension >> (this->BTreeLookUpTableDepth + 1);

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
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#endif
}

// create the Jack polynomial decomposition corresponding to the root partition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& FermionOnSphereHaldaneLargeBasis::GenerateJackPolynomial(RealVector& jack, double alpha)
{
  jack[0] = 1.0;
  double InvAlpha =  2.0 * (1.0 - alpha) / alpha;

  unsigned long* TmpMonomial = new unsigned long [this->NbrFermions];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrFermions];

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, TmpMonomial);
  for (int j = 0; j < this->NbrFermions; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - InvAlpha * ((double) j));
  int ReducedNbrFermions = this->NbrFermions - 1;

  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      double Rho = 0.0;
      unsigned long CurrentPartition = this->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, TmpMonomial);
      for (int j = 0; j < this->NbrFermions; ++j)
	Rho += TmpMonomial[j] * (TmpMonomial[j] - InvAlpha * ((double) j));
      double Coefficient = 0.0;
      for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
	  {
	    double Diff = (double) (TmpMonomial[j1] - TmpMonomial[j2]);
	    unsigned int Max = TmpMonomial[j2];
	    unsigned long TmpState = 0x0ul;
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
			long TmpIndex = this->FindStateIndex(TmpState);
			if (TmpIndex < this->LargeHilbertSpaceDimension)
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

RealVector& FermionOnSphereHaldaneLargeBasis::GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha)
{
  jack[0l] = 1.0;
  double InvAlpha =  2.0 * (1.0 - alpha) / alpha;

  unsigned long* TmpMonomial = new unsigned long [this->NbrFermions];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrFermions];

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, TmpMonomial);
  for (int j = 0; j < this->NbrFermions; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - InvAlpha * ((double) j));
  int ReducedNbrFermions = this->NbrFermions - 1;  
  double SymSign = 1.0;
  if ((((this->NbrFermions * ReducedNbrFermions) >> 1) & 1) != 0)
    SymSign = -1.0;
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (jack[i] == 0.0)
{
      double Rho = 0.0;
      unsigned long CurrentPartition = this->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, TmpMonomial);
      for (int j = 0; j < this->NbrFermions; ++j)
	Rho += TmpMonomial[j] * (TmpMonomial[j] - InvAlpha * ((double) j));
      double Coefficient = 0.0;
      for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
	  {
	    double Diff = (double) (TmpMonomial[j1] - TmpMonomial[j2]);
	    unsigned int Max = TmpMonomial[j2];
	    unsigned long TmpState = 0x0ul;
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
			long TmpIndex = this->FindStateIndex(TmpState);
			if (TmpIndex < this->LargeHilbertSpaceDimension)
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
	  long TmpIndex = this->FindStateIndex(TmpSymState);
	  jack[TmpIndex] = SymSign *Coefficient;
	}
      jack[i] = Coefficient;
}
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

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix FermionOnSphereHaldaneLargeBasis::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState)
{  
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrFermionSector == 0))
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
      if ((lzSector == this->TotalLz) && (nbrFermionSector == this->NbrFermions))
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
	  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	    for (long j = i; j < this->LargeHilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * groundState[j]);
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;  
	}
    }

  unsigned long TmpMask = (((0x1ul << (this->LzMax + 2)) - 1) >> subsytemSize) << subsytemSize;
  unsigned long TmpSubsystemMask = (0x1ul << subsytemSize) - 1;
  long TmpIndex;
  int ShiftedTotalLz = (this->TotalLz + this->NbrFermions * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrFermionSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int TmpStateMaxLz = ShiftedLzComplementarySector - (((NbrFermionsComplementarySector - 2 + (subsytemSize << 1)) * (NbrFermionsComplementarySector - 1)) >> 1);
  long MinIndex = 0;
  long MaxIndex = this->LargeHilbertSpaceDimension - 1l;
//   if ((NbrFermionsComplementarySector > 0) && ((NbrFermionsComplementarySector + subsytemSize - 2) > this->StateLzMax[MaxIndex]))
//     MaxIndex = this->LookUpTable[NbrFermionsComplementarySector + subsytemSize - 2][0];
//   if ((TmpStateMaxLz < this->StateLzMax[0]) && ((TmpStateMaxLz + 1) >  this->StateLzMax[MaxIndex]) && (TmpStateMaxLz >= subsytemSize))
//     MinIndex = this->LookUpTable[TmpStateMaxLz + 1][0];
  
  unsigned long TmpComplementarySubsystem;
  int TmpNbrFermions;
  int TmpTotalLz;
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
  int TmpSumOccupation[] = {
    0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6,
    4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 10, 10,
    5, 5, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 10, 10, 11, 11,
    9, 9, 10, 10, 11, 11, 12, 12, 12, 12, 13, 13, 14, 14, 15, 15,
    6, 6, 7, 7, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 12, 12,
    10, 10, 11, 11, 12, 12, 13, 13, 13, 13, 14, 14, 15, 15, 16, 16,
    11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17,
    15, 15, 16, 16, 17, 17, 18, 18, 18, 18, 19, 19, 20, 20, 21, 21,
    7, 7, 8, 8, 9, 9, 10, 10, 10, 10, 11, 11, 12, 12, 13, 13,
    11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17,
    12, 12, 13, 13, 14, 14, 15, 15, 15, 15, 16, 16, 17, 17, 18, 18,
    16, 16, 17, 17, 18, 18, 19, 19, 19, 19, 20, 20, 21, 21, 22, 22,
    13, 13, 14, 14, 15, 15, 16, 16, 16, 16, 17, 17, 18, 18, 19, 19,
    17, 17, 18, 18, 19, 19, 20, 20, 20, 20, 21, 21, 22, 22, 23, 23,
    18, 18, 19, 19, 20, 20, 21, 21, 21, 21, 22, 22, 23, 23, 24, 24,
    22, 22, 23, 23, 24, 24, 25, 25, 25, 25, 26, 26, 27, 27, 28, 28};
  int TmpPartialNbrOne;
  if (nbrFermionSector <= 1)
    {
      unsigned long Key = 0x0ul;
      if (nbrFermionSector == 1)
	Key = 0x1ul << ShiftedLzSector;
      double TmpValue = 0.0;
      while (MinIndex <= MaxIndex)
	{
	  if ((this->StateDescription[MinIndex] & TmpSubsystemMask) == Key)
	    TmpValue += groundState[MinIndex] * groundState[MinIndex];	    
	  ++MinIndex;
	}
      RealSymmetricMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  else
    {
      FermionOnSphere TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
      RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      long TmpNbrNonZeroElements = 0;
      while (MinIndex <= MaxIndex)
	{
	  TmpIndex = MinIndex;
	  TmpComplementarySubsystem = this->StateDescription[TmpIndex] & TmpMask;
	  ++TmpIndex;
	  while ((TmpIndex <= MaxIndex) && ((this->StateDescription[TmpIndex] & TmpMask) == TmpComplementarySubsystem))
	    ++TmpIndex;
	  TmpPartialNbrOne = TmpNbrOne[TmpComplementarySubsystem & 0xffl];
	  TmpNbrFermions = TmpPartialNbrOne;
	  TmpTotalLz = TmpSumOccupation[TmpComplementarySubsystem & 0xffl];
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 8) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 8) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 3;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 16) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 16) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 4;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 24) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 24) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 24;
#ifdef  __64_BITS__
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 32) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 32) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 5;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 40) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 40) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 40;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 48) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 48) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 48;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 56) & 0xffl];      
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 56) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 56;
#endif
	  if ((TmpNbrFermions == NbrFermionsComplementarySector) && (ShiftedLzComplementarySector == TmpTotalLz))
	    {
	      int Pos = 0;
	      for (long i = MinIndex; i < TmpIndex; ++i)
		{
		  unsigned long TmpState = this->StateDescription[i] & TmpSubsystemMask;
		  int TmpLzMax = subsytemSize - 1;
		  while ((TmpState & (0x1ul << TmpLzMax)) == 0x0ul)
		    --TmpLzMax;
		  TmpStatePosition[Pos] = TmpDestinationHilbertSpace.FindStateIndex(TmpState, TmpLzMax);
		  ++Pos;
		}
	      int Pos2;
	      Pos = 0;
	      int Pos3;
	      double TmpValue;
	      for (long i = MinIndex; i < TmpIndex; ++i)
		{
		  Pos2 = 0;
		  Pos3 = TmpStatePosition[Pos];
		  TmpValue = groundState[i];
		  for (int j = MinIndex; j < TmpIndex; ++j)
		    {
		      if (Pos3 <=  TmpStatePosition[Pos2])
			{
			  TmpDensityMatrix.AddToMatrixElement(Pos3, TmpStatePosition[Pos2], TmpValue * groundState[j]);
			  ++TmpNbrNonZeroElements;
			}
		      ++Pos2;
		    }
		  ++Pos;
		}
	    }
	  MinIndex = TmpIndex;
	}
      delete[] TmpStatePosition;
      if (TmpNbrNonZeroElements > 0)	
	return TmpDensityMatrix;
      else
	{
	  RealSymmetricMatrix TmpDensityMatrixZero;
	  return TmpDensityMatrixZero;
	}
    }
}
