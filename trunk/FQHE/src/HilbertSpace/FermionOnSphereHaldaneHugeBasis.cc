////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//                                                                            //
//                        last modification : 06/07/2006                      //
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
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/ListIterator.h"
#include "GeneralTools/Endian.h"
#include "MathTools/FactorialCoefficient.h" 
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHESphereJackGeneratorOperation.h"


#include <math.h>
#include <cstring>
#include <fstream>
#include <sys/time.h>


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
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// maxFileSize = maximum file size (in MBytes)
// memory = amount of memory granted for precalculations
// referenceState = array that describes the reference state to start from
// symmetricFlag = indicate if a symmetric basis has to be used (only available if totalLz = 0)
// fullDimension = provide the full (i.e. without squeezing) Hilbert space dimension (0 if it has to be computed)

FermionOnSphereHaldaneHugeBasis::FermionOnSphereHaldaneHugeBasis (int nbrFermions, int totalLz, int lzMax, unsigned long maxFileSize, int* referenceState, unsigned long memory,
								  bool symmetricFlag, long fullDimension)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->SizeLimit = maxFileSize << 17;
  this->HilbertSpaceDimension = 0;

  this->NbrPrefixSector = 0;
  this->NbrRootSuffix = 0;
  this->RootSuffixShift = 0;
  this->RootPrefixMask = 0x0ul;
  this->PrefixSectors = 0;
  this->RootSuffix = 0;
  this->RootSuffixSectorPositions = 0;
  this->RootSuffixOffset = 0;
  this->RootSuffixSectorSize = 0;
  this->SuffixLookUpTable = 0;
  this->SuffixLookUpTableSize = 0;
  this->SuffixLookUpTableShift = 0;
  this->TemporaryMonomial = new unsigned long [this->NbrFermions];
  this->TemporaryMonomial2 = new unsigned long [this->NbrFermions];

  this->ReferenceState = 0x0l;
  this->SymmetricFlag = symmetricFlag;
  this->MaximumSignLookUp = 16;
  int ReferenceStateLzMax = 0;
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
  this->Flag.Initialize();

#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;

  if (fullDimension == 0l)
    this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  else 
    this->LargeHilbertSpaceDimension = fullDimension;
  cout << "total dimension = " << this->LargeHilbertSpaceDimension << endl;
  ListIterator<unsigned long> IterFileSizes(this->FileSizes);
  unsigned long* TmpFileSize;
  int MinCommonLz = this->LzMax;
  while ((TmpFileSize = IterFileSizes()))
    {
      if (MinCommonLz > ((int) ((*TmpFileSize) & 0xffl)))
	MinCommonLz = ((*TmpFileSize) & 0xffl);
    }  
  int MaxPartialNbrFermions = this->LzMax - MinCommonLz;  
  this->StateHighestLzShift = MinCommonLz + 1;
  this->LargeHilbertSpaceDimension = 0l;
  int ShiftedTotalLz = (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1;
  this->PartialHilbertSpaceDimension = 0l;
  this->NbrFiles = 0;
  unsigned long LargestFile = 0;
  unsigned long DiskSpace = 0;
  for (int i = MaxPartialNbrFermions; i >= 0; --i)
    {
      int TmpMin = ((i * (i - 1)) >> 1);
      int TmpMax = ((MaxPartialNbrFermions - 1) * i) - ((i * (i - 1)) >> 1);
      int ShiftedTotalLz2 = ShiftedTotalLz - (i * this->StateHighestLzShift);
      int ShiftedTotalLz3  = ShiftedTotalLz2 - (((this->NbrFermions - i) * ((2 * MinCommonLz) - (this->NbrFermions - i) + 1)) >> 1);
      int ShiftedTotalLz4 = ShiftedTotalLz2 - (((this->NbrFermions - i) * (this->NbrFermions - i - 1)) >> 1);
      if (TmpMin < ShiftedTotalLz3)
	TmpMin = ShiftedTotalLz3;
      if (TmpMax > ShiftedTotalLz4)	
	TmpMax = ShiftedTotalLz4;
      if (TmpMax >= TmpMin)
	this->NbrFiles += (TmpMax - TmpMin) + 1;
    }
  this->StateDescriptionFileSizes = new unsigned long[this->NbrFiles];
  this->NbrFiles = 0;
  for (int i = MaxPartialNbrFermions; i >= 0; --i)
    {
      int TmpMin = ((i * (i - 1)) >> 1);
      int TmpMax = ((MaxPartialNbrFermions - 1) * i) - ((i * (i - 1)) >> 1);
      int ShiftedTotalLz2 = ShiftedTotalLz - (i * this->StateHighestLzShift);
      int ShiftedTotalLz3  = ShiftedTotalLz2 - (((this->NbrFermions - i) * ((2 * MinCommonLz) - (this->NbrFermions - i) + 1)) >> 1);
      int ShiftedTotalLz4 = ShiftedTotalLz2 - (((this->NbrFermions - i) * (this->NbrFermions - i - 1)) >> 1);
      if (TmpMin < ShiftedTotalLz3)
	TmpMin = ShiftedTotalLz3;
      if (TmpMax > ShiftedTotalLz4)	
	TmpMax = ShiftedTotalLz4;
      for (; TmpMax >= TmpMin; --TmpMax)
	{
          unsigned long TmpPartialDimension1 = this->ShiftedEvaluateHilbertSpaceDimension2(i, MaxPartialNbrFermions - 1, TmpMax);	  
	  unsigned long TmpPartialDimension2 = this->ShiftedEvaluateHilbertSpaceDimension2(this->NbrFermions - i, MinCommonLz, ShiftedTotalLz2 - TmpMax);
	  if (TmpPartialDimension1 == 0l)
	    TmpPartialDimension1 = 1l;
	  if (LargestFile < TmpPartialDimension2)
	    LargestFile = TmpPartialDimension2;
	  this->LargeHilbertSpaceDimension += TmpPartialDimension1 * TmpPartialDimension2;
	  this->PartialHilbertSpaceDimension += TmpPartialDimension1;
	  DiskSpace += TmpPartialDimension2;
	  this->StateDescriptionFileSizes[this->NbrFiles] = TmpPartialDimension2;
	  ++this->NbrFiles;
	}
    }
  cout << "total dimension = " << this->LargeHilbertSpaceDimension << "  (" << this->PartialHilbertSpaceDimension << ")" << endl;
  cout << "total requested disk space = " << (DiskSpace >> 17) <<  " Mbytes splitted in " << this->NbrFiles << " files " << endl;
  cout << "largest file = " << (LargestFile >> 7) << " KBytes (" << LargestFile<< " elements)" << endl;
  if (memory > (DiskSpace << 3))
    {
      cout << "disk storage is useless, keep everything in memory" << endl;
      this->NoDiskFlag = true;
    }
  else
    this->NoDiskFlag = false;
  this->StateDescription = new unsigned long [this->PartialHilbertSpaceDimension];
  this->StateDescriptionFileIndex = new int [this->PartialHilbertSpaceDimension];
  this->StateDescriptionIndexShift = new long [this->PartialHilbertSpaceDimension];
  this->PartialHilbertSpaceDimension = 0l;
  this->NbrFiles = 0;
  this->StateDescriptionIndexShift[0] = 0l;
  this->StateHighestLzToIndex = new int [1 << ((this->LzMax - MinCommonLz) + 2)];
  for (int i = MaxPartialNbrFermions; i >= 0; --i)
    {
      long TmpPartialHilbertSpaceDimension;
      int TmpMin = ((i * (i - 1)) >> 1);
      int TmpMax = ((MaxPartialNbrFermions - 1) * i) - ((i * (i - 1)) >> 1);
      int ShiftedTotalLz2 = ShiftedTotalLz - (i * this->StateHighestLzShift);
      int ShiftedTotalLz3  = ShiftedTotalLz2 - (((this->NbrFermions - i) * ((2 * MinCommonLz) - (this->NbrFermions - i) + 1)) >> 1);
      int ShiftedTotalLz4 = ShiftedTotalLz2 - (((this->NbrFermions - i) * (this->NbrFermions - i - 1)) >> 1);
      if (TmpMin < ShiftedTotalLz3)
	TmpMin = ShiftedTotalLz3;
      if (TmpMax > ShiftedTotalLz4)	
	TmpMax = ShiftedTotalLz4;
      for (; TmpMax >= TmpMin; --TmpMax)	  
	{
	  TmpPartialHilbertSpaceDimension = this->RawGenerateStates(i, MaxPartialNbrFermions - 1, MaxPartialNbrFermions - 1, TmpMax, this->PartialHilbertSpaceDimension);
	  for (; this->PartialHilbertSpaceDimension < TmpPartialHilbertSpaceDimension; ++this->PartialHilbertSpaceDimension)
	    this->StateDescriptionFileIndex[this->PartialHilbertSpaceDimension] = this->NbrFiles;
	  ++this->NbrFiles;	  
	}      
    }

  SortArrayDownOrdering(this->StateDescription, this->StateDescriptionFileIndex, this->PartialHilbertSpaceDimension);
  this->StateDescriptionIndexShift[0] = 0l;
  this->StateHighestLzToIndex[this->StateDescription[0]] = 0;
  for (long i = 1; i < this->PartialHilbertSpaceDimension; ++i)
    {
      this->StateDescriptionIndexShift[i] = this->StateDescriptionIndexShift[i - 1] + this->StateDescriptionFileSizes[this->StateDescriptionFileIndex[i - 1]];
      this->StateHighestLzToIndex[this->StateDescription[i]] = i;
    }
  this->HighestLzStateMask = ~((0x1l << this->StateHighestLzShift) - 1);
  this->NbrBuffers = memory / (LargestFile << 3);
  if (this->NbrBuffers < 2)
    this->NbrBuffers = 2;
  this->StateDescriptionBuffers = new unsigned long* [this->NbrBuffers];
  this->StateDescriptionLzSectorBuffers = new unsigned long* [this->NbrBuffers];
  for (int i = 0; i < this->NbrBuffers; ++i)
    {
      this->StateDescriptionBuffers[i] = new unsigned long[LargestFile];
      this->StateDescriptionLzSectorBuffers[i] = new unsigned long [this->StateHighestLzShift + 1];
    }
  this->StateDescriptionFileNames = new char*[this->NbrFiles];
  int CurrentFileIndex = 0;
  unsigned long* TmpStateDescription = this->StateDescription;
  for (int i = MaxPartialNbrFermions; i >= 0; --i)
    {
      int TmpMin = ((i * (i - 1)) >> 1);
      int TmpMax = ((MaxPartialNbrFermions - 1) * i) - ((i * (i - 1)) >> 1);
      int ShiftedTotalLz2 = ShiftedTotalLz - (i * this->StateHighestLzShift);
      int ShiftedTotalLz3  = ShiftedTotalLz2 - (((this->NbrFermions - i) * ((2 * MinCommonLz) - (this->NbrFermions - i) + 1)) >> 1);
      int ShiftedTotalLz4 = ShiftedTotalLz2 - (((this->NbrFermions - i) * (this->NbrFermions - i - 1)) >> 1);
      if (TmpMin < ShiftedTotalLz3)
	TmpMin = ShiftedTotalLz3;
      if (TmpMax > ShiftedTotalLz4)	
	TmpMax = ShiftedTotalLz4;
      for (; TmpMax >= TmpMin; --TmpMax)	  
	{
	  this->StateDescription = this->StateDescriptionBuffers[0];
	  this->RawGenerateStates(this->NbrFermions - i, MinCommonLz, MinCommonLz, ShiftedTotalLz2 - TmpMax, 0);
	  this->FindLzMaxSectors(this->StateDescription, this->StateDescriptionFileSizes[CurrentFileIndex], this->StateDescriptionLzSectorBuffers[0], this->StateHighestLzShift - 1);
	  this->StateDescriptionFileNames[CurrentFileIndex] = new char[256];
	  sprintf (this->StateDescriptionFileNames[CurrentFileIndex], "fermions_sphere_n_%d_2s_%d.%d.hb", this->NbrFermions, this->LzMax, CurrentFileIndex);	  
	  ofstream File; 
	  File.open(this->StateDescriptionFileNames[CurrentFileIndex], ios::binary | ios::out);
	  for (int i = 0; i <= this->StateHighestLzShift; ++i)
	    WriteLittleEndian(File, this->StateDescriptionLzSectorBuffers[0][i]);
	  for (unsigned long i = 0; i < this->StateDescriptionFileSizes[CurrentFileIndex]; ++i)
	    WriteLittleEndian(File, this->StateDescription[i]);
	  File.close();
	  ++CurrentFileIndex;
	}      
    }
  this->StateDescription = TmpStateDescription;
  this->BufferIndices = new int [this->NbrBuffers];
  this->BufferAges = new int [this->NbrBuffers];
  for (int i = 1; i < this->NbrBuffers; ++i)
    {
      this->BufferIndices[i] = -1;
      this->BufferAges[i] = this->NbrBuffers + 1;
    }
  this->BufferIndices[0] = this->NbrFiles - 1;
  this->BufferAges[0] = 1;
  this->FileToBuffer = new int [this->NbrFiles];
  for (int i = 0; i < this->NbrFiles; ++i)
    this->FileToBuffer[i] = -1;
  this->FileToBuffer[this->NbrFiles - 1] = 0;

#ifdef  __64_BITS__
   unsigned long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 6) + 1l;
#else
   unsigned long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 5) + 1l;
#endif
   this->KeepStateFlag = new unsigned long [ReducedHilbertSpaceDimension];
   for (unsigned long i = 0; i < ReducedHilbertSpaceDimension; ++i)
     this->KeepStateFlag[i] = 0x0l;

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
   this->GenerateStates(ReferenceStateLzMax, this->ReferenceState, 1l, Memory);

   if ((this->SymmetricFlag == false) || (this->TotalLz != 0))
     {
       unsigned long NewHilbertSpaceDimension = 0;
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
       for (unsigned long i = 0; i < ReducedHilbertSpaceDimension; ++i)
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
       cout << "Haldane space dimension = " << NewHilbertSpaceDimension << endl;
       this->FullLargeHilbertSpaceDimension =  this->LargeHilbertSpaceDimension;
       this->LargeHilbertSpaceDimension = NewHilbertSpaceDimension;
       if (this->LargeHilbertSpaceDimension >= (1l << 30))
	 this->HilbertSpaceDimension = 0;
       else
	 this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
     }
   else
     {
       cout << "symmetrizing..." << endl;
       unsigned long NewHilbertSpaceDimension = 0;
       unsigned long NewTotalHilbertSpaceDimension = 0;
       long TotalIndex = 0;
       long LocalIndex = 0;
       unsigned long TmpStateHighestLz = this->StateDescription[0] << this->StateHighestLzShift;
       int TmpFileIndex = this->StateDescriptionFileIndex[this->StateHighestLzToIndex[this->StateDescription[0]]];
       long TmpLimit = this->StateDescriptionFileSizes[TmpFileIndex];
       int TmpBufferIndex = this->FileToBuffer[TmpFileIndex];
       if (TmpBufferIndex < 0)
	 TmpBufferIndex = this->LoadLowestLzBuffer(TmpFileIndex);
       unsigned long* TmpStateDescriptionBuffers = this->StateDescriptionBuffers[TmpBufferIndex];
       for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	 {
#ifdef  __64_BITS__
	   unsigned long& TmpKeepStateFlag = this->KeepStateFlag[i >> 6];	   
#else
	   unsigned long& TmpKeepStateFlag = this->KeepStateFlag[i >> 5];	   
#endif	   
	   if (this->IsCanonicalState(TmpStateHighestLz | TmpStateDescriptionBuffers[LocalIndex]))
	     {
	       ++NewTotalHilbertSpaceDimension;
#ifdef  __64_BITS__
	       if (((TmpKeepStateFlag >> (i & 0x3ful)) & 0x1ul) != 0ul)
		 ++NewHilbertSpaceDimension;
	       else
		 TmpKeepStateFlag &= ~(0x1ul << (i & 0x3ful));
#else
	       if (((TmpKeepStateFlag >> (i & 0x1ful)) & 0x1ul) != 0ul)
		 ++NewHilbertSpaceDimension;
	       else
		 TmpKeepStateFlag &= ~(0x1ul << (i & 0x1ful));
#endif
	     }
	   ++LocalIndex;
	   if (LocalIndex == TmpLimit)
	     {
	       ++TotalIndex;
	       if (TotalIndex < this->PartialHilbertSpaceDimension)
		 {
		   TmpStateHighestLz = this->StateDescription[TotalIndex] << this->StateHighestLzShift;
		   TmpFileIndex = this->StateDescriptionFileIndex[this->StateHighestLzToIndex[this->StateDescription[TotalIndex]]];
		   TmpLimit = this->StateDescriptionFileSizes[TmpFileIndex];
		   TmpBufferIndex = this->FileToBuffer[TmpFileIndex];
		   if (TmpBufferIndex < 0)
		     TmpBufferIndex = this->LoadLowestLzBuffer(TmpFileIndex);
		   TmpStateDescriptionBuffers = this->StateDescriptionBuffers[TmpBufferIndex];
		   LocalIndex = 0ul;
		 }
	     }
	 }
       cout << "total symmetric space dimension = " << NewTotalHilbertSpaceDimension << endl;
       cout << "Haldane symmetric space dimension = " << NewHilbertSpaceDimension << endl;
       this->FullLargeHilbertSpaceDimension =  this->LargeHilbertSpaceDimension;
       this->LargeHilbertSpaceDimension = NewHilbertSpaceDimension;
       if (this->LargeHilbertSpaceDimension >= (1l << 30))
	 this->HilbertSpaceDimension = 0;
       else
	 this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
     }

  cout << "final Hilbert space dimension " << this->LargeHilbertSpaceDimension << endl;
}

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memoryHilbert = amount of memory granted to store the Hilbert space (in Mbytes)
// memory = amount of memory allowed for precalculations

FermionOnSphereHaldaneHugeBasis::FermionOnSphereHaldaneHugeBasis(char* fileName, long memoryHilbert, long memory)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      this->HilbertSpaceDimension = 0;
      return;
    }
  this->NbrBuffers = 0;
  this->NbrFiles = 0;
  this->StateDescriptionBuffers = 0;
  this->StateDescriptionFileNames = 0;
  this->StateDescriptionFileSizes = 0;
  this->KeepStateFlag = 0;
  ReadLittleEndian(File, this->HilbertSpaceDimension);
  ReadLittleEndian(File, this->LargeHilbertSpaceDimension);
  ReadLittleEndian(File, this->NbrFermions);
  ReadLittleEndian(File, this->LzMax);
  ReadLittleEndian(File, this->TotalLz);
  ReadLittleEndian(File, this->ReferenceState);
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
  if ((this->LargeHilbertSpaceDimension << 3) < (memoryHilbert << 20))
    {
      this->NbrPrefixSector = 0;
      this->NbrRootSuffix = 0;
      this->RootSuffixShift = 0;
      this->RootPrefixMask = 0x0ul;
      this->PrefixSectors = 0;
      this->RootSuffix = 0;
      this->RootSuffixSectorPositions = 0;
      this->RootSuffixOffset = 0;
      this->RootSuffixSectorSize = 0;
      this->SuffixLookUpTable = 0;
      this->SuffixLookUpTableSize = 0;
      this->SuffixLookUpTableShift = 0;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	ReadLittleEndian(File, this->StateDescription[i]);
      File.close();
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      unsigned long UsedMemory = 0l;
      UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
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
  else 
    {
      File.close();
      this->HilbertSpaceFileName = new char[strlen(fileName) + 1];
      strcpy (this->HilbertSpaceFileName, fileName);
      this->FileHeaderSize = (4 * sizeof(int)) + (2 * sizeof(long));
      this->StateDescription = 0;
      this->SparseHilbertSpaceDimension = memoryHilbert << 17;
      unsigned long CurrentPartition;
      unsigned long Mask = ((1l << ((this->LzMax / 2) + 2)) - 1l) << (this->LzMax / 2);
      this->RootSuffixShift = LzMax / 2;
      this->RootPrefixMask = (1l << this->RootSuffixShift) - 1l;
      CurrentPartition &= Mask;
      long Count = 1l;
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
      this->NbrRootSuffix = 1l;
      long** SectorSize = new long* [this->NbrFermions + 1];
      unsigned int*** Sectors = new unsigned int**[this->NbrFermions + 1];
      int TmpMaxTotalLz = ((this->LzMax / 2) + 1) * this->NbrFermions - ((this->NbrFermions * (this->NbrFermions - 1)) / 2);
      for (int i = 0; i <= this->NbrFermions; ++i)
	{
	  SectorSize[i] = new long [TmpMaxTotalLz + 1];
	  Sectors[i] = new unsigned int*[TmpMaxTotalLz +1];
	  for (int j = 0; j <= TmpMaxTotalLz; ++j)
	    {
	      SectorSize[i][j] = 0l;
	      Sectors[i][j] = 0;
	    }
	}

      ifstream FileHilbert;      
      FileHilbert.open(this->HilbertSpaceFileName, ios::binary | ios::in);
      FileHilbert.seekg (this->FileHeaderSize, ios::beg);
      ReadLittleEndian(FileHilbert, CurrentPartition);	  
      CurrentPartition >>= this->RootSuffixShift;
      int TmpNbrFermions = 0;
      int TmpTotalLz = 0;
      int TmpPartialNbrOne = 0;
      TmpPartialNbrOne = TmpNbrOne[CurrentPartition & 0xffl];
      TmpNbrFermions = TmpPartialNbrOne;
      TmpTotalLz = TmpSumOccupation[CurrentPartition & 0xffl];
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 8) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 8) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne << 3;
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 16) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 16) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne << 4;
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 24) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 24) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne * 24;
#ifdef  __64_BITS__
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 32) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 32) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne << 5;
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 40) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 40) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne * 40;
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 48) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 48) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne * 48;
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 56) & 0xffl];      
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 56) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne * 56;
#endif
      for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  unsigned long MaxRoot;
	  ReadLittleEndian(FileHilbert, MaxRoot);	  
	  if (CurrentPartition == (MaxRoot >> this->RootSuffixShift))
	    {
	      Count++;
	    }
	  else
	    {
	      ++this->NbrRootSuffix;
	      if (SectorSize[TmpNbrFermions][TmpTotalLz] == 0l)
		{
		  SectorSize[TmpNbrFermions][TmpTotalLz] = Count;
		}
	      CurrentPartition = (MaxRoot >> this->RootSuffixShift);
	      TmpPartialNbrOne = TmpNbrOne[CurrentPartition & 0xffl];
	      TmpNbrFermions = TmpPartialNbrOne;
	      TmpTotalLz = TmpSumOccupation[CurrentPartition & 0xffl];
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 8) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 8) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne << 3;
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 16) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 16) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne << 4;
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 24) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 24) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne * 24;
#ifdef  __64_BITS__
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 32) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 32) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne << 5;
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 40) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 40) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne * 40;
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 48) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 48) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne * 48;
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 56) & 0xffl];      
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 56) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne * 56;
#endif
	      Count = 1l;
	    }	  
	}
      if (SectorSize[TmpNbrFermions][TmpTotalLz] == 0l)
	{
	  SectorSize[TmpNbrFermions][TmpTotalLz] = Count;
	}
      FileHilbert.close();
      cout << "TotalCount = "  << this->NbrRootSuffix << endl;      

      long SumSector = 0l;
      this->NbrPrefixSector = 0l;
      for (int i = 0; i <= this->NbrFermions; ++i)
	{
	  for (int j = 0; j <= TmpMaxTotalLz; ++j)
	    {
	      SumSector += SectorSize[i][j];
	      if (SectorSize[i][j] != 0l)
		{
		  Sectors[i][j] = new unsigned int [SectorSize[i][j]];
		  ++this->NbrPrefixSector;
		}
	      else
		Sectors[i][j] = 0;
	      SectorSize[i][j] = 0l;
	    }
	}
      this->PrefixSectors = new unsigned int* [this->NbrPrefixSector];
      this->RootSuffix = new unsigned int[this->NbrRootSuffix];
      this->RootSuffixSectorPositions = new unsigned int*[this->NbrRootSuffix];
      this->RootSuffixOffset = new long [this->NbrRootSuffix + 1];
      this->RootSuffixSectorSize = new long [this->NbrRootSuffix];
      FileHilbert.open(this->HilbertSpaceFileName, ios::binary | ios::in);
      FileHilbert.seekg (this->FileHeaderSize, ios::beg);
      ReadLittleEndian(FileHilbert, CurrentPartition);	  
      unsigned long MaxRoot = CurrentPartition;
      CurrentPartition >>= this->RootSuffixShift;
      this->RootSuffix[0] = (unsigned int) CurrentPartition;
      this->RootSuffixOffset[0] = 0l;
      TmpPartialNbrOne = TmpNbrOne[CurrentPartition & 0xffl];
      TmpNbrFermions = TmpPartialNbrOne;
      TmpTotalLz = TmpSumOccupation[CurrentPartition & 0xffl];
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 8) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 8) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne << 3;
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 16) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 16) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne << 4;
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 24) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 24) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne * 24;
#ifdef  __64_BITS__
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 32) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 32) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne << 5;
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 40) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 40) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne * 40;
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 48) & 0xffl];
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 48) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne * 48;
      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 56) & 0xffl];      
      TmpNbrFermions += TmpPartialNbrOne;
      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 56) & 0xffl];
      TmpTotalLz += TmpPartialNbrOne * 56;
#endif
      unsigned int* CurrentSector = Sectors[TmpNbrFermions][TmpTotalLz];
      CurrentSector[0] = MaxRoot & this->RootPrefixMask;
      this->RootSuffixSectorPositions[0] = CurrentSector;
      this->NbrRootSuffix = 0;
      Count = 1;
      for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  ReadLittleEndian(FileHilbert, MaxRoot);	  
	  if (CurrentPartition == (MaxRoot >> this->RootSuffixShift))
	    {
	      if (CurrentSector != 0)
		{
		  CurrentSector[Count] = (unsigned int) (MaxRoot & this->RootPrefixMask);
		  ++Count;
		}
	    }
	  else
	    {
	      if (CurrentSector == 0)
		{
		  this->RootSuffixSectorSize[this->NbrRootSuffix] = SectorSize[TmpNbrFermions][TmpTotalLz];
		}
	      else
		{
		  this->RootSuffixSectorSize[this->NbrRootSuffix] = Count;
		}
	      ++this->NbrRootSuffix;
	      this->RootSuffix[this->NbrRootSuffix] = (unsigned int) (MaxRoot >> this->RootSuffixShift);
	      if (CurrentSector != 0)
		{
		  SectorSize[TmpNbrFermions][TmpTotalLz] = Count;
		}
	      CurrentPartition = (MaxRoot >> this->RootSuffixShift);
	      TmpPartialNbrOne = TmpNbrOne[CurrentPartition & 0xffl];
	      TmpNbrFermions = TmpPartialNbrOne;
	      TmpTotalLz = TmpSumOccupation[CurrentPartition & 0xffl];
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 8) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 8) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne << 3;
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 16) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 16) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne << 4;
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 24) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 24) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne * 24;
#ifdef  __64_BITS__
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 32) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 32) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne << 5;
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 40) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 40) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne * 40;
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 48) & 0xffl];
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 48) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne * 48;
	      TmpPartialNbrOne = TmpNbrOne[(CurrentPartition >> 56) & 0xffl];      
	      TmpNbrFermions += TmpPartialNbrOne;
	      TmpTotalLz += TmpSumOccupation[(CurrentPartition >> 56) & 0xffl];
	      TmpTotalLz += TmpPartialNbrOne * 56;
#endif
	      Count = 0;
	      if (SectorSize[TmpNbrFermions][TmpTotalLz] == 0)
		{
		  CurrentSector = Sectors[TmpNbrFermions][TmpTotalLz];
		  CurrentSector[Count] = (unsigned int) (MaxRoot & this->RootPrefixMask);
		  ++Count;
		}
	      else
		{
		  CurrentSector = 0;
		}
	      this->RootSuffixSectorPositions[this->NbrRootSuffix] = Sectors[TmpNbrFermions][TmpTotalLz];
	      this->RootSuffixOffset[this->NbrRootSuffix] = i;
	    }	  
	}
      if (CurrentSector != 0)
	{
	  SectorSize[TmpNbrFermions][TmpTotalLz] = Count;
	}
      this->RootSuffixSectorSize[this->NbrRootSuffix] = SectorSize[TmpNbrFermions][TmpTotalLz];
      ++this->NbrRootSuffix;
      this->RootSuffixOffset[this->NbrRootSuffix] = this->RootSuffixOffset[this->NbrRootSuffix - 1l];
      FileHilbert.close();
      this->NbrPrefixSector = 0l;
      for (int i = 0; i <= this->NbrFermions; ++i)
	{
	  for (int j = 0; j <= TmpMaxTotalLz; ++j)
	    {
	      if (SectorSize[i][j] != 0l)
		{
		  this->PrefixSectors[this->NbrPrefixSector] = Sectors[i][j];
		  this->NbrPrefixSector++;
		}
	    }
	  delete[] Sectors[i];
	  delete[] SectorSize[i];
	}
      delete[] Sectors;
      delete[] SectorSize;	  
      cout << "Factorization ratio = " << (((SumSector + this->NbrRootSuffix) * 100.0) / ((double) this->LargeHilbertSpaceDimension)) << "%" << endl;
      this->GenerateLookUpTableFactorized();
    }
  this->TemporaryMonomial = new unsigned long [this->NbrFermions];
  this->TemporaryMonomial2 = new unsigned long [this->NbrFermions];

}



// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereHaldaneHugeBasis::FermionOnSphereHaldaneHugeBasis(const FermionOnSphereHaldaneHugeBasis& fermions)
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
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->Flag = fermions.Flag;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->KeepStateFlag = fermions.KeepStateFlag;
  this->NbrBuffers = fermions.NbrBuffers;
  this->NbrFiles = fermions.NbrFiles;
  this->StateDescriptionBuffers = fermions.StateDescriptionBuffers;
  this->StateDescriptionFileNames = fermions.StateDescriptionFileNames;
  this->StateDescriptionFileSizes = fermions.StateDescriptionFileSizes;
  this->NbrPrefixSector  = fermions.NbrPrefixSector;
  this->NbrRootSuffix  = fermions.NbrRootSuffix;
  this->RootSuffixShift  = fermions.RootSuffixShift;
  this->RootPrefixMask  = fermions.RootPrefixMask;
  this->PrefixSectors  = fermions.PrefixSectors;
  this->RootSuffix  = fermions.RootSuffix;
  this->RootSuffixSectorPositions  = fermions.RootSuffixSectorPositions;
  this->RootSuffixOffset  = fermions.RootSuffixOffset;
  this->RootSuffixSectorSize  = fermions.RootSuffixSectorSize;
  this->SuffixLookUpTable  = fermions.SuffixLookUpTable;
  this->SuffixLookUpTableSize  = fermions.SuffixLookUpTableSize;
  this->SuffixLookUpTableShift  = fermions.SuffixLookUpTableShift;
  this->TemporaryMonomial = new unsigned long [this->NbrFermions];
  this->TemporaryMonomial2 = new unsigned long [this->NbrFermions];
}

// destructor
//

FermionOnSphereHaldaneHugeBasis::~FermionOnSphereHaldaneHugeBasis ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      if (NbrBuffers > 0)
	{
	  for (int i = 0; i < this->NbrBuffers; ++i)
	    delete[] this->StateDescriptionBuffers[i];
	  delete[] this->StateDescriptionBuffers;
	}
      if  (this->NbrFiles > 0)
	{
	  for (int i = 0; i < this->NbrFiles; ++i)
	    delete[] this->StateDescriptionFileNames[i];
	  delete[] this->StateDescriptionFileNames;
	  delete[] this->StateDescriptionFileSizes;
	}
      if (this->KeepStateFlag != 0)
	delete[] this->KeepStateFlag;
    }
  delete[] this->TemporaryMonomial;
  delete[] this->TemporaryMonomial2;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereHaldaneHugeBasis& FermionOnSphereHaldaneHugeBasis::operator = (const FermionOnSphereHaldaneHugeBasis& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
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
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->Flag = fermions.Flag;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->KeepStateFlag = fermions.KeepStateFlag;
  this->NbrPrefixSector  = fermions.NbrPrefixSector;
  this->NbrRootSuffix  = fermions.NbrRootSuffix;
  this->RootSuffixShift  = fermions.RootSuffixShift;
  this->RootPrefixMask  = fermions.RootPrefixMask;
  this->PrefixSectors  = fermions.PrefixSectors;
  this->RootSuffix  = fermions.RootSuffix;
  this->RootSuffixSectorPositions  = fermions.RootSuffixSectorPositions;
  this->RootSuffixOffset  = fermions.RootSuffixOffset;
  this->RootSuffixSectorSize  = fermions.RootSuffixSectorSize;
  this->SuffixLookUpTable  = fermions.SuffixLookUpTable;
  this->SuffixLookUpTableSize  = fermions.SuffixLookUpTableSize;
  this->SuffixLookUpTableShift  = fermions.SuffixLookUpTableShift;
  this->TemporaryMonomial = new unsigned long [this->NbrFermions];
  this->TemporaryMonomial2 = new unsigned long [this->NbrFermions];
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereHaldaneHugeBasis::Clone()
{
  return new FermionOnSphereHaldaneHugeBasis(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphereHaldaneHugeBasis::WriteHilbertSpace (char* fileName)
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

  long TotalIndex = 0;
  long LocalIndex = 0;
  unsigned long TmpStateHighestLz = this->StateDescription[0] << this->StateHighestLzShift;
  int TmpFileIndex = this->StateDescriptionFileIndex[this->StateHighestLzToIndex[this->StateDescription[0]]];
  long TmpLimit = this->StateDescriptionFileSizes[TmpFileIndex];
  int TmpBufferIndex = this->FileToBuffer[TmpFileIndex];
  if (TmpBufferIndex < 0)
    TmpBufferIndex = this->LoadLowestLzBuffer(TmpFileIndex);
  unsigned long* TmpStateDescriptionBuffers = this->StateDescriptionBuffers[TmpBufferIndex];
  long Count = 0;;
  for (long i = 0; i < this->FullLargeHilbertSpaceDimension; ++i)
    {
#ifdef  __64_BITS__
      unsigned long& TmpKeepStateFlag = this->KeepStateFlag[i >> 6];	   
#else
      unsigned long& TmpKeepStateFlag = this->KeepStateFlag[i >> 5];	   
#endif	   
#ifdef  __64_BITS__
      if (((TmpKeepStateFlag >> (i & 0x3ful)) & 0x1ul) != 0ul)
	{
	  unsigned long TmpState = TmpStateHighestLz | TmpStateDescriptionBuffers[LocalIndex];
	  WriteLittleEndian(File, TmpState);
	  Count++;
	}
#else
      if (((TmpKeepStateFlag >> (i & 0x1ful)) & 0x1ul) != 0ul)
	{
	  unsigned long TmpState = TmpStateHighestLz | TmpStateDescriptionBuffers[LocalIndex];
	  WriteLittleEndian(File, TmpState);
	  Count++;
	}
#endif
      ++LocalIndex;
      if (LocalIndex == TmpLimit)
	{
	  ++TotalIndex;
	  if (TotalIndex < this->PartialHilbertSpaceDimension)
	    {
	      TmpStateHighestLz = this->StateDescription[TotalIndex] << this->StateHighestLzShift;
	      TmpFileIndex = this->StateDescriptionFileIndex[this->StateHighestLzToIndex[this->StateDescription[TotalIndex]]];
	      TmpLimit = this->StateDescriptionFileSizes[TmpFileIndex];
	      TmpBufferIndex = this->FileToBuffer[TmpFileIndex];
	      if (TmpBufferIndex < 0)
		TmpBufferIndex = this->LoadLowestLzBuffer(TmpFileIndex);
	      TmpStateDescriptionBuffers = this->StateDescriptionBuffers[TmpBufferIndex];
	      LocalIndex = 0ul;
	    }
	}
    }  
  File.close();
  cout << "Nbr Saved states = " << Count << endl;
  return true;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnSphereHaldaneHugeBasis::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnSphereHaldaneHugeBasis::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphereHaldaneHugeBasis::ExtractSubspace (AbstractQuantumNumber& q, 
							SubspaceSpaceConverter& converter)
{
  return 0;
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnSphereHaldaneHugeBasis::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->TargetSpace = (FermionOnSphereHaldaneHugeBasis*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int FermionOnSphereHaldaneHugeBasis::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// convert a given state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereHaldaneHugeBasis::ConvertToNbodyBasis(RealVector& state, FermionOnSphereHaldaneHugeBasis& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetLargeHilbertSpaceDimension(), true);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      TmpVector[nbodyBasis.FindStateIndexFactorized(this->GetStateFactorized(i))] = state[i];
    }
  return TmpVector;
}

// convert a given state from the usual n-body basis to the Haldane basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereHaldaneHugeBasis::ConvertFromNbodyBasis(RealVector& state, FermionOnSphereHaldaneHugeBasis& nbodyBasis)
{
  RealVector TmpVector (this->LargeHilbertSpaceDimension, true);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    TmpVector[i] = state[nbodyBasis.FindStateIndexFactorized(this->GetStateFactorized(i))];
  TmpVector /= TmpVector.Norm();
  return TmpVector;
}

// convert a given state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

LongRationalVector FermionOnSphereHaldaneHugeBasis::ConvertToNbodyBasis(LongRationalVector& state, FermionOnSphereHaldaneHugeBasis& nbodyBasis)
{
  LongRationalVector TmpVector (nbodyBasis.GetLargeHilbertSpaceDimension(), true);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      TmpVector[nbodyBasis.FindStateIndexFactorized(this->GetStateFactorized(i))] = state[i];
    }
  return TmpVector;
}

// convert a given state from the usual n-body basis to the Haldane basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

LongRationalVector FermionOnSphereHaldaneHugeBasis::ConvertFromNbodyBasis(LongRationalVector& state, FermionOnSphereHaldaneHugeBasis& nbodyBasis)
{
  LongRationalVector TmpVector (this->LargeHilbertSpaceDimension, true);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    TmpVector[i] = state[nbodyBasis.FindStateIndexFactorized(this->GetStateFactorized(i))];
  TmpVector /= TmpVector.Norm();
  return TmpVector;
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

int FermionOnSphereHaldaneHugeBasis::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return 0;
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereHaldaneHugeBasis::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  return 0;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereHaldaneHugeBasis::ProdA (int index, int* n, int nbrIndices)
{
  return 0;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereHaldaneHugeBasis::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  return 0;
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereHaldaneHugeBasis::AdA (int index, int m)
{
  if ((this->GetStateFactorized((long) index) & (0x1ul << m)) != 0ul)
    return 1.0;
  else
    return 0.0;
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereHaldaneHugeBasis::AdA (long index, int m)
{
  if ((this->GetStateFactorized(index) & (0x1ul << m)) != 0ul)
    return 1.0;
  else
    return 0.0;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
  
int FermionOnSphereHaldaneHugeBasis::AdA (int index, int m, int n, double& coefficient)
{
  return 0;
}

// find state index assuming the Hilbert space is stored using the sparse technique
//
// stateDescription = unsigned integer describing the state
// return value = corresponding index

long FermionOnSphereHaldaneHugeBasis::FindStateIndex(unsigned long stateDescription)
{
  unsigned long TmpHighestLz = stateDescription >> this->StateHighestLzShift;
  int TmpBufferIndex =  this->StateHighestLzToIndex[TmpHighestLz];
  TmpHighestLz = this->StateDescriptionIndexShift[TmpBufferIndex];
  int TmpFileIndex = this->StateDescriptionFileIndex[TmpBufferIndex];
  TmpBufferIndex = this->FileToBuffer[TmpFileIndex];
  if (TmpBufferIndex < 0)
    TmpBufferIndex = this->LoadLowestLzBuffer(TmpFileIndex);

  stateDescription &= ~this->HighestLzStateMask;
  int TmpLzMax = this->StateHighestLzShift;
  while (((stateDescription >> TmpLzMax) & 0x1l) == 0l)
    --TmpLzMax;
  long PosMax = this->StateDescriptionLzSectorBuffers[TmpBufferIndex][TmpLzMax + 1];
  long PosMin = this->StateDescriptionLzSectorBuffers[TmpBufferIndex][TmpLzMax] - 1l;
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long* TmpStateDescriptionBuffers = this->StateDescriptionBuffers[TmpBufferIndex];
  unsigned long CurrentState = TmpStateDescriptionBuffers[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	PosMax = PosMid;
      else
	PosMin = PosMid;
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = TmpStateDescriptionBuffers[PosMid];
    }
  if (CurrentState == stateDescription)
    return (TmpHighestLz + PosMid);
  else
    return (TmpHighestLz + PosMin);
}

// load a lowest Lz part file into memory
//
// fileIndex = index of the file to read
// return file = index of the buffer that has been loaded

int FermionOnSphereHaldaneHugeBasis::LoadLowestLzBuffer(int fileIndex)
{
  int Index = 0;  
  while (this->BufferAges[Index] < this->NbrBuffers)
    ++Index;
  if (this->BufferIndices[Index] >= 0)
    this->FileToBuffer[this->BufferIndices[Index]] = -1;
  unsigned long* TmpBuffer = this->StateDescriptionBuffers[Index];
  unsigned long* TmpBuffer2 = this->StateDescriptionLzSectorBuffers[Index];
  unsigned long TmpSize = this->StateDescriptionFileSizes[fileIndex];
  ifstream File;
  File.open(this->StateDescriptionFileNames[fileIndex], ios::binary | ios::in);
  for (int i = 0; i <= this->StateHighestLzShift; ++i)
    ReadLittleEndian(File, TmpBuffer2[i]);
  for (unsigned long i = 0; i < TmpSize; ++i)
    ReadLittleEndian(File, TmpBuffer[i]);
  File.close();
  this->BufferAges[Index] = 0;
  this->BufferIndices[Index] = fileIndex;
  for (int i = 0; i < this->NbrBuffers; ++i)
    ++this->BufferAges[i];
  this->FileToBuffer[this->BufferIndices[Index]] = Index;
  return Index;
}

// find state index assuming the whole Hilbert space is stored in memory
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

long FermionOnSphereHaldaneHugeBasis::FindStateIndexMemory(unsigned long stateDescription, int lzmax)
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
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->LargeHilbertSpaceDimension;
    else
      return PosMin;
}

// find state index when hilbert space storage is based on factorized algorithm
//
// stateDescription = unsigned integer describing the state
// return value = corresponding index

long FermionOnSphereHaldaneHugeBasis::FindStateIndexFactorized(unsigned long stateDescription)
{
  unsigned int TmpSuffix = (unsigned int) (stateDescription >> this->RootSuffixShift);
  unsigned int CurrentState = TmpSuffix >> this->SuffixLookUpTableShift;
  if (CurrentState >= this->SuffixLookUpTableSize)
    return this->LargeHilbertSpaceDimension;
  long PosMax = this->SuffixLookUpTable[CurrentState + 1];//0l;
  long PosMin = this->SuffixLookUpTable[CurrentState];//this->NbrRootSuffix - 1l;
  long PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->RootSuffix[PosMid];
  while ((PosMax != PosMid) && (CurrentState != TmpSuffix))
    {
      if (CurrentState > TmpSuffix)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->RootSuffix[PosMid];
    }
  long TmpSuffixPos = PosMin;
  if (CurrentState == TmpSuffix)
    TmpSuffixPos = PosMid;
  else
    if ((this->RootSuffix[PosMin] != TmpSuffix) && (this->RootSuffix[PosMax] != TmpSuffix))
      return this->LargeHilbertSpaceDimension;
  TmpSuffix = (unsigned int) (stateDescription & this->RootPrefixMask);
  PosMax = 0l;
  PosMin = this->RootSuffixSectorSize[TmpSuffixPos] - 1l;
  PosMid = (PosMin + PosMax) >> 1;
  unsigned int* TmpPrefixSector = this->RootSuffixSectorPositions[TmpSuffixPos];
  CurrentState = TmpPrefixSector[PosMid];
  while ((PosMax != PosMid) && (CurrentState != TmpSuffix))
    {
      if (CurrentState > TmpSuffix)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = TmpPrefixSector[PosMid];
    }
  if (CurrentState == TmpSuffix)
    return PosMid + this->RootSuffixOffset[TmpSuffixPos];
  else
    if ((TmpPrefixSector[PosMin] != TmpSuffix) && (TmpPrefixSector[PosMax] != TmpSuffix))
      return this->LargeHilbertSpaceDimension;
    else
      return PosMin + this->RootSuffixOffset[TmpSuffixPos];      
}

// find multiple state indices when hilbert space storage is based on factorized algorithm
//
// stateDescriptions = array of unsigned integer describing the states (state has to be sorted from the largest to the smallest one)
// nbrStates = number of states to process
// indices = array where state indices will be stored

void FermionOnSphereHaldaneHugeBasis::FindMultipleStateIndexFactorized(unsigned long* stateDescriptions, int nbrStates, long* indices)
{
  unsigned int TmpSuffix = (unsigned int) (stateDescriptions[0] >> this->RootSuffixShift);
  unsigned int CurrentState = TmpSuffix >> this->SuffixLookUpTableShift;
  long PosMax = this->SuffixLookUpTable[CurrentState + 1];
  long PosMin = this->SuffixLookUpTable[CurrentState];
  long PosMid = (PosMin + PosMax) >> 1;
  long PosPrefixMax = 0l;
  long PosPrefixMin = 0l;
  long PosPrefixMid = 0l;
  unsigned int* TmpPrefixSector = 0;
  unsigned int TmpPrefix = 0;
  CurrentState = this->RootSuffix[PosMid];
  while ((PosMax != PosMid) && (CurrentState != TmpSuffix))
    {
      if (CurrentState > TmpSuffix)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->RootSuffix[PosMid];
    }
  long TmpSuffixPos = PosMin;
  if ((this->RootSuffix[PosMin] != TmpSuffix) && (this->RootSuffix[PosMax] != TmpSuffix))
    {
      TmpSuffixPos = PosMid;
      indices[0] = this->LargeHilbertSpaceDimension;
    }
  else
    {
      if (CurrentState == TmpSuffix)
	TmpSuffixPos = PosMid;
      TmpPrefix = (unsigned int) (stateDescriptions[0] & this->RootPrefixMask);
      PosPrefixMax = 0l;
      PosPrefixMin = this->RootSuffixSectorSize[TmpSuffixPos] - 1l;
      PosPrefixMid = (PosPrefixMin + PosPrefixMax) >> 1;
      TmpPrefixSector = this->RootSuffixSectorPositions[TmpSuffixPos];
      CurrentState = TmpPrefixSector[PosPrefixMid];
      while ((PosPrefixMax != PosPrefixMid) && (CurrentState != TmpPrefix))
	{
	  if (CurrentState > TmpPrefix)
	    {
	      PosPrefixMax = PosPrefixMid;
	    }
	  else
	    {
	      PosPrefixMin = PosPrefixMid;
	    } 
	  PosPrefixMid = (PosPrefixMin + PosPrefixMax) >> 1;
	  CurrentState = TmpPrefixSector[PosPrefixMid];
	}
      if (CurrentState == TmpPrefix)
	{
	  indices[0] = PosPrefixMid + this->RootSuffixOffset[TmpSuffixPos];
	}
      else
	if ((TmpPrefixSector[PosPrefixMin] != TmpPrefix) && (TmpPrefixSector[PosPrefixMax] != TmpPrefix))
	  indices[0] = this->LargeHilbertSpaceDimension;
	else
	  {
	    PosPrefixMax = PosPrefixMin;
	    indices[0] = PosPrefixMin + this->RootSuffixOffset[TmpSuffixPos];      
	  }
    }
  for (int i = 1; i < nbrStates; ++i)
    {
      unsigned int TmpSuffix2 = (unsigned int) (stateDescriptions[i] >> this->RootSuffixShift);
      if (TmpSuffix == TmpSuffix2)
	{
	  if (TmpPrefixSector != 0)
	    {
	      PosPrefixMin = this->RootSuffixSectorSize[TmpSuffixPos] - 1l;
	      PosPrefixMid = (PosPrefixMin + PosPrefixMax) >> 1;
	      TmpPrefix = (unsigned int) (stateDescriptions[i] & this->RootPrefixMask);
	      CurrentState = TmpPrefixSector[PosPrefixMid];
	      while ((PosPrefixMax != PosPrefixMid) && (CurrentState != TmpPrefix))
		{
		  if (CurrentState > TmpPrefix)
		    {
		      PosPrefixMax = PosPrefixMid;
		    }
		  else
		    {
		      PosPrefixMin = PosPrefixMid;
		    } 
		  PosPrefixMid = (PosPrefixMin + PosPrefixMax) >> 1;
		  CurrentState = TmpPrefixSector[PosPrefixMid];
		}
	      if (CurrentState == TmpPrefix)
		{
		  indices[i] = PosPrefixMid + this->RootSuffixOffset[TmpSuffixPos];
		}
	      else
		{
		  PosPrefixMax = PosPrefixMid;
		  if ((TmpPrefixSector[PosPrefixMin] != TmpPrefix) && (TmpPrefixSector[PosPrefixMax] != TmpPrefix))
		    {
		      indices[i] = this->LargeHilbertSpaceDimension;
		    }
		  else
		    {
		      indices[i] = PosPrefixMin + this->RootSuffixOffset[TmpSuffixPos];      
		    }	  
		}
	    }
	  else
	    indices[i] = this->LargeHilbertSpaceDimension;
	}
      else
	{
	  TmpSuffix = TmpSuffix2;
	  CurrentState = TmpSuffix >> this->SuffixLookUpTableShift;
	  PosMin = this->SuffixLookUpTable[CurrentState + 1];
	  if (PosMin > PosMax)
	    PosMax = PosMin;
	  PosMin = this->SuffixLookUpTable[CurrentState];
	  PosMid = (PosMin + PosMax) >> 1;
	  CurrentState = this->RootSuffix[PosMid];
	  while ((PosMax != PosMid) && (CurrentState != TmpSuffix))
	    {
	      if (CurrentState > TmpSuffix)
		{
		  PosMax = PosMid;
		}
	      else
		{
		  PosMin = PosMid;
		} 
	      PosMid = (PosMin + PosMax) >> 1;
	      CurrentState = this->RootSuffix[PosMid];
	    }
	  TmpSuffixPos = PosMin;
	  if ((this->RootSuffix[PosMin] != TmpSuffix) && (this->RootSuffix[PosMax] != TmpSuffix))
	    {
	      TmpSuffixPos = PosMid;
	      TmpPrefixSector = 0;
	      indices[i] = this->LargeHilbertSpaceDimension;
	    }
	  else
	    {
	      if (CurrentState == TmpSuffix)
		TmpSuffixPos = PosMid;
	      TmpPrefix = (unsigned int) (stateDescriptions[i] & this->RootPrefixMask);
	      PosPrefixMax = 0l;
	      PosPrefixMin = this->RootSuffixSectorSize[TmpSuffixPos] - 1l;
	      PosPrefixMid = (PosPrefixMin + PosPrefixMax) >> 1;
	      TmpPrefixSector = this->RootSuffixSectorPositions[TmpSuffixPos];
	      CurrentState = TmpPrefixSector[PosPrefixMid];
	      while ((PosPrefixMax != PosPrefixMid) && (CurrentState != TmpPrefix))
		{
		  if (CurrentState > TmpPrefix)
		    {
		      PosPrefixMax = PosPrefixMid;
		    }
		  else
		    {
		      PosPrefixMin = PosPrefixMid;
		    } 
		  PosPrefixMid = (PosPrefixMin + PosPrefixMax) >> 1;
		  CurrentState = TmpPrefixSector[PosPrefixMid];
		}
	      if (CurrentState == TmpPrefix)
		{
		  indices[i] = PosPrefixMid + this->RootSuffixOffset[TmpSuffixPos];
		}
	      else
		if ((TmpPrefixSector[PosPrefixMin] != TmpPrefix) && (TmpPrefixSector[PosPrefixMax] != TmpPrefix))
		  indices[i] = this->LargeHilbertSpaceDimension;
		else
		  {
		    PosPrefixMax = PosPrefixMin;
		    indices[i] = PosPrefixMin + this->RootSuffixOffset[TmpSuffixPos];      
		  }
	    }
	}
    }
}

// get a state description from its index when hilbert space storage is based on factorized algorithm
//
// index = state index
// return value = unsigned integer describing the state

unsigned long FermionOnSphereHaldaneHugeBasis::GetStateFactorized(long index)
{
  long PosMax = 0l;
  long PosMin = this->NbrRootSuffix;
  long PosMid = (PosMin + PosMax) >> 1;
  long CurrentIndex = this->RootSuffixOffset[PosMid];
  while ((PosMin - PosMax) > 1l)
    {
      if (CurrentIndex <= index)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentIndex = this->RootSuffixOffset[PosMid];
    }
//  cout << PosMin << " " << PosMid << " " << PosMax << endl;
  return ((((unsigned long) this->RootSuffix[PosMax]) << this->RootSuffixShift) | ((unsigned long) this->RootSuffixSectorPositions[PosMax][index - this->RootSuffixOffset[PosMax]]));
}

// find state index when hilbert space storage is based on sparse algorithm
//
// stateDescription = unsigned integer describing the state
// return value = corresponding index

long FermionOnSphereHaldaneHugeBasis::FindStateIndexSparse(unsigned long stateDescription)
{
  long PosMax = this->SparseHilbertSpaceDimension;
  long PosMin = 0l;
  long PosMid = (PosMin + PosMax) >> 1;
  cout << PosMin << " " << PosMax << " " << hex << stateDescription << dec << endl;
  unsigned long CurrentState;
  while ((PosMax - PosMin) > 1)
    {
      CurrentState = this->SparseHilbertSpaceDescription[PosMid];
      if (CurrentState > stateDescription)
	PosMin = PosMid;
      else
	PosMax = PosMid;
      PosMid = (PosMin + PosMax) >> 1;
    }
  unsigned long* TmpBuffer;
  cout << PosMin << " " << PosMax << " " << hex << stateDescription << dec << endl;
  this->LoadSparseBuffer(PosMin, TmpBuffer, PosMax);
  PosMin *= this->SparseHilbertSpaceChunckSize;
  long PosMin2 = 0l;
  --PosMax;
  PosMid = (PosMin2 + PosMax) >> 1;
  CurrentState = TmpBuffer[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	PosMin2 = PosMid;
      else
	PosMax = PosMid;
      PosMid = (PosMin2 + PosMax) >> 1;
      CurrentState = TmpBuffer[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid + PosMin;
  else
    return PosMin2 + PosMin;
}

// load a part of the Hilbert space associated to one element of the sparse Hilbert basis
//
// sparseIndex = index associated to the element of the sparse Hilbert space 
// buffer = reference on the pointer to the part of the total space 
// bufferSize = reference on the size of the buffer
 
void FermionOnSphereHaldaneHugeBasis::LoadSparseBuffer(long sparseIndex, unsigned long*& buffer, long& bufferSize)
{
  int TmpIndex = this->SparseHilbertSpaceBufferIndices[sparseIndex];
  if (sparseIndex == (SparseHilbertSpaceDimension - 1))
    bufferSize = this->SparseHilbertSpaceRemainderChunckSize;
  else
    bufferSize = this->SparseHilbertSpaceChunckSize;
  if (TmpIndex >= 0)
    {
      buffer = this->SparseBuffers[TmpIndex];
      for (int i = 0; i < this->NbrBuffers; ++i)
	--this->BufferAges[i];
      ++this->BufferAges[TmpIndex];
      return;
    }
  int TmpMinAge = this->BufferAges[0];
  int TmpIndex2 = 0l;
  for (int i = 1; i < this->NbrBuffers; ++i)
    if (this->BufferAges[i] < TmpMinAge)
      {
	TmpIndex2 = i;
	TmpMinAge = this->BufferAges[i];
      }
  for (int i = 0; i < this->NbrBuffers; ++i)
    --this->BufferAges[i];
  this->BufferAges[TmpIndex2] = 0;
  ifstream File;
  File.open(this->HilbertSpaceFileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << this->HilbertSpaceFileName << endl;
      return;
    }
  File.seekg (this->FileHeaderSize, ios::beg);
  File.seekg (this->SparseHilbertSpaceChunckSize * sparseIndex, ios::beg);
  unsigned long* TmpBuffer = this->SparseBuffers[TmpIndex2];;
  for (long i = 0; i < bufferSize; ++i)
    ReadLittleEndian(File, TmpBuffer[i]);
  File.close();
  
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereHaldaneHugeBasis::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->NbrLzValue; ++i)
    Str << ((TmpState >> i) & ((unsigned long) 0x1)) << " ";
  return Str;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereHaldaneHugeBasis::PrintStateMonomial (ostream& Str, long state)
{
  unsigned long TmpState = this->StateDescription[state];
  Str << "[";
  int i = this->LzMax;
  while (((TmpState >> i) & 0x1ul) == 0x0ul)
    --i;
  Str << i;
  --i;
  for (; i >=0; --i)
    if (((TmpState >> i) & 0x1ul) != 0x0ul)
      Str << "," << i;
  Str << "]";
  return Str;
}

// print a given State using the monomial notation, with one column per particle (using space as a seperator)
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereHaldaneHugeBasis::PrintColumnFormattedStateMonomial (ostream& Str, long state)
{
  unsigned long TmpState = this->GetStateFactorized(state);
  int i = this->LzMax;
  while (((TmpState >> i) & 0x1ul) == 0x0ul)
    --i;
  Str << i;
  --i;
  for (; i >=0; --i)
    if (((TmpState >> i) & 0x1ul) != 0x0ul)
      Str << " " << i;
  return Str;
}

// generate all states corresponding to the constraints
// 
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereHaldaneHugeBasis::GenerateStates(int lzMax, unsigned long referenceState, long pos, long& memory)
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


// generate all states (i.e. all possible skew symmetric polynomials with fixed Lz)
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// currentLzMax = momentum maximum value for fermions that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereHaldaneHugeBasis::RawGenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, long pos)
{
  if ((nbrFermions == 0) && (totalLz == 0))
    {
      this->StateDescription[pos] = 0l;
      return pos + 1l;
    }
  if ((nbrFermions == 0) || (totalLz < 0) || (currentLzMax < (nbrFermions - 1)))
    return pos;
  int LzTotalMax = ((2 * currentLzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (totalLz > LzTotalMax)
    return pos;
  if ((nbrFermions == 1) && (currentLzMax >= totalLz))
    {
      this->StateDescription[pos] = (0x1l << totalLz);
      return pos + 1l;
    }
  if (LzTotalMax == totalLz)
    {
      unsigned long Mask = 0l;
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

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereHaldaneHugeBasis::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->NbrLzValue);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrLzValue)
    this->MaximumLookUpShift = this->NbrLzValue;
  this->LookUpTableMemorySize = 1l << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new long* [this->NbrLzValue];
  this->LookUpTableShift = new int [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->LookUpTable[i] = new long [this->LookUpTableMemorySize + 1];
  int CurrentLzMax = this->LzMax;
  unsigned long CurrentState = this->StateDescription[0];
  while (((CurrentState >> CurrentLzMax) & 0x1ul) == 0x0ul)
    --CurrentLzMax;
  int TmpLzMax = CurrentLzMax;
  long* TmpLookUpTable = this->LookUpTable[CurrentLzMax];
  if (CurrentLzMax < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentLzMax] = 0;
  else
    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentLzMax];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = CurrentState >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      CurrentState = this->StateDescription[i];
      while (((CurrentState >> TmpLzMax) & 0x1ul) == 0x0ul)
	--TmpLzMax;
      if (CurrentLzMax != TmpLzMax)
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  /*	  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
	    cout << TmpLookUpTable[j] << " ";
	    cout << endl << "-------------------------------------------" << endl;*/
 	  CurrentLzMax = TmpLzMax;
	  TmpLookUpTable = this->LookUpTable[CurrentLzMax];
	  if (CurrentLzMax < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentLzMax] = 0;
	  else
	    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentLzMax];
	  TmpLookUpTableValue = CurrentState >> CurrentShift;
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
	  TmpLookUpTableValue = CurrentState >> CurrentShift;
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
  this->GenerateSignLookUpTable();
}

// generate look-up table associated to current Hilbert space in factorized mode
// 

void FermionOnSphereHaldaneHugeBasis::GenerateLookUpTableFactorized()
{
  int TmpLookUpSize = 20;
  this->SuffixLookUpTableShift = this->LzMax + 1 - this->RootSuffixShift - TmpLookUpSize;
  if (this->SuffixLookUpTableShift < 0)
    {
      this->SuffixLookUpTableShift = 0;
      this->SuffixLookUpTableSize = 1 << (TmpLookUpSize + this->SuffixLookUpTableShift);
    }
  else
    {
      this->SuffixLookUpTableSize = 1 << TmpLookUpSize;
    }
  cout << "look up size " <<  this->SuffixLookUpTableSize << endl;
  this->SuffixLookUpTable = new long [this->SuffixLookUpTableSize + 1];
  unsigned int TmpSuffix = this->RootSuffix[0] >> this->SuffixLookUpTableShift;  
  unsigned int CurrentSuffix = this->SuffixLookUpTableSize;
  while (CurrentSuffix > TmpSuffix)
    {
      this->SuffixLookUpTable[CurrentSuffix] = 0l;
      --CurrentSuffix;
    }
  this->SuffixLookUpTable[CurrentSuffix] = 0l;
  for (long i = 0l; i < this->NbrRootSuffix; ++i)
    {
      TmpSuffix = this->RootSuffix[i] >> this->SuffixLookUpTableShift;
      if (TmpSuffix != CurrentSuffix)
	{
	  while (CurrentSuffix > TmpSuffix)
	    {
	      this->SuffixLookUpTable[CurrentSuffix] = i;
	      --CurrentSuffix;
	    }
	  this->SuffixLookUpTable[CurrentSuffix] = i;
	}
    }
  while (CurrentSuffix > 0)
    {
      this->SuffixLookUpTable[CurrentSuffix] = this->NbrRootSuffix - 1l;
      --CurrentSuffix;
    }
  this->SuffixLookUpTable[0] = this->NbrRootSuffix - 1l;
  this->GenerateSignLookUpTable();
}

// generate look-up table for sign calculation
// 

void FermionOnSphereHaldaneHugeBasis::GenerateSignLookUpTable()
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


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereHaldaneHugeBasis::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + nbrFermions * lzMax) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long FermionOnSphereHaldaneHugeBasis::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)))
    return 0l;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return 0l;
  if ((nbrFermions == 1) && (lzMax >= totalLz))
    return 1l;
  if (LzTotalMax == totalLz)
    return 1l;
  long Tmp1 = this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax);
  long Tmp2 = this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz);  
  if ((Tmp1 < this->SizeLimit) && (Tmp2 < this->SizeLimit) && ((Tmp1 + Tmp2) > this->SizeLimit))
    {
      this->FileSizes += ((((unsigned long) (nbrFermions - 1)) << 8) | ((unsigned long) (lzMax - 1)) | (((unsigned long) (totalLz - lzMax)) << 16));
      this->FileSizes += ((((unsigned long) nbrFermions) << 8) | ((unsigned long) (lzMax - 1)) | (((unsigned long) totalLz) << 16));
    }
  return (Tmp1 + Tmp2);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long FermionOnSphereHaldaneHugeBasis::ShiftedEvaluateHilbertSpaceDimension2(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)))
    return 0l;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return 0l;
  if ((nbrFermions == 1) && (lzMax >= totalLz))
    return 1l;
  if (LzTotalMax == totalLz)
    return 1l;
  return (this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax) +
	  this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz));  
}

// find position of each sector with a given maximum Lz
// 
// stateDescription = array that contains state description
// dimension = dimension of the stateDescription array
// lzSectors = array that where the position of each sector with a given maximum Lz will be stored
// lzMax = maximum momentum value for a fermion

void FermionOnSphereHaldaneHugeBasis::FindLzMaxSectors(unsigned long* stateDescription, long dimension, unsigned long* lzSectors, int lzMax)
{
  for (int i = 0; i <= lzMax; ++i)
    lzSectors[i] = 0xffffffffl;
  unsigned long TmpStateDescription;
  int CurrentLzMax = lzMax;
  int TmpLzMax;
  long Position = 0l;
  lzSectors[lzMax + 1] = 0l;
  while (Position < dimension)
    {
      TmpLzMax = CurrentLzMax;
      TmpStateDescription = stateDescription[Position];
      while (((TmpStateDescription >> TmpLzMax) & 0x1l) == 0x0l)
	--TmpLzMax;
      if (CurrentLzMax != TmpLzMax)
	{
	  lzSectors[CurrentLzMax] = Position;
	  lzSectors[TmpLzMax + 1] = Position;
	  CurrentLzMax = TmpLzMax;
	}
      ++Position;
    }
  lzSectors[CurrentLzMax] = dimension;
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphereHaldaneHugeBasis::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
							       int firstComponent, int nbrComponent)
{
  Complex Value;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereHaldaneHugeBasis::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
  
// create the Jack polynomial decomposition corresponding to the root partition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// partialSave = save partial results in a given vector file
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& FermionOnSphereHaldaneHugeBasis::GenerateJackPolynomial(RealVector& jack, double alpha, long minIndex, long maxIndex, char* partialSave)
{
  jack[0] = 1.0;
  double InvAlpha =  2.0 * (1.0 - alpha) / alpha;

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->TemporaryMonomial);
  for (int j = 0; j < this->NbrFermions; ++j)
    RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - InvAlpha * ((double) j));
  int ReducedNbrFermions = this->NbrFermions - 1;

  if (minIndex <= 0)
    minIndex = 1;
  if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension - 1l;
  for (long i = minIndex; i <= maxIndex; ++i)
    {
      double Rho = 0.0;
      unsigned long CurrentPartition = this->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, this->TemporaryMonomial);
      for (int j = 0; j < this->NbrFermions; ++j)
	Rho += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - InvAlpha * ((double) j));
      double Coefficient = 0.0;
      for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
	  {
	    double Diff = (double) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
	    unsigned int Max = this->TemporaryMonomial[j2];
	    unsigned long TmpState = 0x0ul;
	    int Tmpj1 = j1;
	    int Tmpj2 = j2;
	    for (int l = 0; l < this->NbrFermions; ++l)
	      this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
	    double Sign = 1.0;
	    for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
	      {
		++this->TemporaryMonomial2[Tmpj1];
		--this->TemporaryMonomial2[Tmpj2];
		while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] >= this->TemporaryMonomial2[Tmpj1 - 1]))
		  {
		    unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
		    this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
		    this->TemporaryMonomial2[Tmpj1] = Tmp;
		    --Tmpj1;
		    Sign *= -1.0; 
		  }
                while ((Tmpj2 < ReducedNbrFermions) && (this->TemporaryMonomial2[Tmpj2] <= this->TemporaryMonomial2[Tmpj2 + 1]))
                  {
                    unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
                    this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
                    this->TemporaryMonomial2[Tmpj2] = Tmp;
                    ++Tmpj2;
 		    Sign *= -1.0; 
                 }
		if ((this->TemporaryMonomial2[Tmpj1] != this->TemporaryMonomial2[Tmpj1 + 1]) && (this->TemporaryMonomial2[Tmpj2] != this->TemporaryMonomial2[Tmpj2 - 1]))
		  {
		    TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
		    if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		      {
			long TmpIndex = this->FindStateIndexMemory(TmpState, this->TemporaryMonomial2[0]);
			if (TmpIndex < this->LargeHilbertSpaceDimension)
			  Coefficient += Sign * Diff * jack[TmpIndex];
		      }
		  }
	      }
	  }
      jack[i] = Coefficient * InvAlpha / (RhoRoot - Rho);
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	  if ((partialSave != 0) && ((i & 0xffffffl) == 0l))
	    jack.WriteVector(partialSave);
	}
    }
  cout << endl;

  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// partialSave = save partial results in a given vector file
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& FermionOnSphereHaldaneHugeBasis::GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha, long minIndex, long maxIndex, char* partialSave)
{
  jack[0] = 1.0;
  double InvAlpha =  2.0 * (1.0 - alpha) / alpha;

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->TemporaryMonomial);
  for (int j = 0; j < this->NbrFermions; ++j)
    RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - InvAlpha * ((double) j));
  int ReducedNbrFermions = this->NbrFermions - 1;  
  double SymSign = 1.0;
  if ((((this->NbrFermions * ReducedNbrFermions) >> 1) & 1) != 0)
    SymSign = -1.0;
  if (minIndex <= 0)
    minIndex = 1;
  if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension - 1l;
  for (long i = minIndex; i <= maxIndex; ++i)
    {
      double Rho = 0.0;
      unsigned long CurrentPartition = this->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, this->TemporaryMonomial);
      for (int j = 0; j < this->NbrFermions; ++j)
	Rho += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - InvAlpha * ((double) j));
      double Coefficient = 0.0;
      for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
	  {
	    double Diff = (double) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
	    unsigned int Max = this->TemporaryMonomial[j2];
	    unsigned long TmpState = 0x0ul;
	    int Tmpj1 = j1;
	    int Tmpj2 = j2;
	    for (int l = 0; l < this->NbrFermions; ++l)
	      this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
	    double Sign = 1.0;
	    for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
	      {
		++this->TemporaryMonomial2[Tmpj1];
		--this->TemporaryMonomial2[Tmpj2];
		while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] >= this->TemporaryMonomial2[Tmpj1 - 1]))
		  {
		    unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
		    this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
		    this->TemporaryMonomial2[Tmpj1] = Tmp;
		    --Tmpj1;
		    Sign *= -1.0; 
		  }
                while ((Tmpj2 < ReducedNbrFermions) && (this->TemporaryMonomial2[Tmpj2] <= this->TemporaryMonomial2[Tmpj2 + 1]))
                  {
                    unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
                    this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
                    this->TemporaryMonomial2[Tmpj2] = Tmp;
                    ++Tmpj2;
 		    Sign *= -1.0; 
                 }
		if ((this->TemporaryMonomial2[Tmpj1] != this->TemporaryMonomial2[Tmpj1 + 1]) && (this->TemporaryMonomial2[Tmpj2] != this->TemporaryMonomial2[Tmpj2 - 1]))
		  {
		    TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
		    if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		      {
			long TmpIndex = this->FindStateIndexMemory(TmpState, this->TemporaryMonomial2[0]);
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
	  long TmpIndex = this->FindStateIndexMemory(TmpSymState, this->LzMax - this->TemporaryMonomial[ReducedNbrFermions]);
	  jack[TmpIndex] = SymSign *Coefficient;
	}
      jack[i] = Coefficient;
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	  if ((partialSave != 0) && ((i & 0xffffffl) == 0l))
	    jack.WriteVector(partialSave);
	}
    }
  cout << endl;

  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alphaNumerator = numerator of the Jack polynomial alpha coefficient
// alphaDenominator = numerator of the Jack polynomial alpha coefficient
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// partialSave = save partial results in a given vector file
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

LongRationalVector& FermionOnSphereHaldaneHugeBasis::GenerateJackPolynomial(LongRationalVector& jack, long alphaNumerator, long alphaDenominator, long minIndex, long maxIndex, char* partialSave)
{
  jack[0l] = 1l;
  LongRational InvAlpha (2l * (alphaDenominator - alphaNumerator), alphaNumerator);

  int ReducedNbrFermions = this->NbrFermions - 1;
  long* ConnectedIndices = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients  = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];
  long* ConnectedIndices2 = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients2  = new long [((this->NbrFermions * ReducedNbrFermions) >> 1) * (this->LzMax + 1)];

  LongRational RhoRoot = 0l;
  LongRational Rho = 0l;
  unsigned long MaxRoot = this->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->TemporaryMonomial);
  for (int j = 0; j < this->NbrFermions; ++j)
    RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - InvAlpha * ((long) j));

  LongRational Coefficient = 0l;
  LongRational Coefficient2 = 0l;

  if (minIndex <= 0)
    minIndex = 1;
  if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension - 1l;
  for (long i = minIndex; i <= maxIndex; ++i)
    {
      Rho = 0l;
      unsigned long CurrentPartition = this->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, this->TemporaryMonomial);
      for (int j = 0; j < this->NbrFermions; ++j)
	Rho += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - InvAlpha * ((long) j));
      if (Rho == RhoRoot)
	{
	  cout << "warning : singular value detected at position " << i << ", skipping the rest of the calculation" << endl;
	  return jack;
	}
      else
	{
	  Coefficient = 0l;
	  int Pos = 0;
	  for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
	    for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
	      {
		long Diff = (long) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
		unsigned int Max = this->TemporaryMonomial[j2];
		unsigned long TmpState = 0x0ul;
		int Tmpj1 = j1;
		int Tmpj2 = j2;
		for (int l = 0; l < this->NbrFermions; ++l)
		  this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
		long Sign = 1l;
		for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		  {
		    ++this->TemporaryMonomial2[Tmpj1];
		    --this->TemporaryMonomial2[Tmpj2];
		    while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] >= this->TemporaryMonomial2[Tmpj1 - 1]))
		      {
			unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
			this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
			this->TemporaryMonomial2[Tmpj1] = Tmp;
			--Tmpj1;
			Sign *= -1l; 
		      }
		    while ((Tmpj2 < ReducedNbrFermions) && (this->TemporaryMonomial2[Tmpj2] <= this->TemporaryMonomial2[Tmpj2 + 1]))
		      {
			unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
			this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
			this->TemporaryMonomial2[Tmpj2] = Tmp;
			++Tmpj2;
			Sign *= -1l; 
		      }
		    if ((this->TemporaryMonomial2[Tmpj1] != this->TemporaryMonomial2[Tmpj1 + 1]) && (this->TemporaryMonomial2[Tmpj2] != this->TemporaryMonomial2[Tmpj2 - 1]))
		      {
			TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
			if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
			  {
			    long TmpIndex = this->FindStateIndexMemory(TmpState, this->TemporaryMonomial2[0]);
			    if (TmpIndex < this->LargeHilbertSpaceDimension)
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
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	  if ((partialSave != 0) && ((i & 0xffffffl) == 0l))
	    jack.WriteVector(partialSave);
	}
    }
  cout << endl;

  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition and using sparse storage
//
// alpha = value of the Jack polynomial alpha coefficient
// architecture = architecture to use for precalculation
// partialSave = save partial results in a given vector file
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// memory = amount of memory (in bytes) allowed for temporary vector storage (0 if the whole vector has to be stored in memory)
// memoryBlock = amount of memory (in bytes) allowed for precomputing state indices
// resumeFlag = true if the calculation has to be resumed from a previous one (assuming partialSave contains already computed components)

void FermionOnSphereHaldaneHugeBasis::GenerateJackPolynomialSparse(double alpha, AbstractArchitecture* architecture, char* partialSave, long minIndex, long maxIndex, long memory, long memoryBlock, bool resumeFlag)
{
  if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension - 1l;
  double TmpComponent = 1.0;
  long FileShift = 4l;
  if (this->HilbertSpaceDimension <= 0)
    FileShift = 12l;
  if ((minIndex <= 0) && (resumeFlag == false))
    {
      ofstream File;
      File.open(partialSave, ios::binary | ios::out);
      WriteLittleEndian(File, this->HilbertSpaceDimension);  
      if (this->HilbertSpaceDimension <= 0)
	{
	  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);  
	}
      WriteLittleEndian(File, TmpComponent);  
      File.close();
    }
  TmpComponent = 0.0;
   
  memory >>= 3;
  if ((memory > this->LargeHilbertSpaceDimension) || (memory <= 0))
    {
      cout << "vector does not require temporary disk storage" << endl;
      memory = this->LargeHilbertSpaceDimension;
    }
  double* TmpVectorBuffer = new double [memory];
  long BufferGlobalIndex = 0l;

  double InvAlpha =  2.0 * (1.0 - alpha) / alpha;

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->GetStateFactorized(0l);;
  this->ConvertToMonomial(MaxRoot, this->TemporaryMonomial);
  for (int j = 0; j < this->NbrFermions; ++j)
    RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - InvAlpha * ((double) j));
  int ReducedNbrFermions = this->NbrFermions - 1;  

  if (minIndex <= 0)
    minIndex = 1;
  TmpVectorBuffer[0l] = 1.0;
  int MaxArraySize = ((this->NbrFermions * (this->NbrFermions - 1)) / 2) * (this->LzMax + 1);
  long NbrBlocks =  memoryBlock / (((2l* sizeof (long)) + (sizeof(double))) * MaxArraySize);
  if (NbrBlocks == 0)
    {
      NbrBlocks = 100000;
    }
  if (NbrBlocks > this->LargeHilbertSpaceDimension)
    {
      NbrBlocks = (this->LargeHilbertSpaceDimension >> 1) + 1;
    }
  cout << "number of precalculation blocks = " << NbrBlocks << endl;
  long DisplayStep = (this->LargeHilbertSpaceDimension / (1000 * NbrBlocks)) * NbrBlocks;
  long** TmpIndexArray = new long* [NbrBlocks];
  double** TmpComponentArray = new double* [NbrBlocks];
  unsigned long** TmpStateArray = new unsigned long* [NbrBlocks]; 
  double* TmpRhoArray = new double [NbrBlocks];
  int* TmpNbrComputedComponentArray = new int [NbrBlocks];
  for (int j = 0; j < NbrBlocks; ++j)
    {
      TmpIndexArray[j] = new long [MaxArraySize];
      TmpComponentArray[j] = new double [MaxArraySize];
      TmpStateArray[j] = new unsigned long [MaxArraySize];
    }
  FQHESphereJackGeneratorOperation Operation(this, InvAlpha, MaxRoot, TmpIndexArray, TmpStateArray, TmpComponentArray, TmpRhoArray, TmpNbrComputedComponentArray, true, false);

  if (resumeFlag == true)
    {
      ifstream File;
      File.open(partialSave, ios::binary | ios::in);
      File.seekg (0, ios::end);
      long TmpResumePos = File.tellg();
      File.close();
      TmpResumePos -= FileShift;
      TmpResumePos /= sizeof(double); 	      
      long TmpResumeMinPos = TmpResumePos - NbrBlocks;
      long LimNbrBlocks = NbrBlocks;
      if (TmpResumeMinPos < 0l)
	{
	  TmpResumeMinPos = 0l;
	  LimNbrBlocks = TmpResumePos - TmpResumeMinPos + 1;
	}
      long TmpMaxIndex = TmpResumeMinPos + NbrBlocks - 1l;
      if (TmpMaxIndex > maxIndex)
	{
	  LimNbrBlocks = NbrBlocks - (TmpMaxIndex - maxIndex);
	  TmpMaxIndex = maxIndex;
	}
      if (LimNbrBlocks > 0)
	{
	  cout << "consistency check, " << TmpResumePos << " components have already been computed, checking the last " << LimNbrBlocks << " ones" << endl;      
	  Operation.SetIndicesRange(TmpResumeMinPos, LimNbrBlocks);
	  Operation.ApplyOperation(architecture);
	  ifstream OutputFile;
	  OutputFile.open(partialSave, ios::binary | ios::in);
	  double RefCoefficient = 0.0;

	  for (long k = 0l; k < LimNbrBlocks; ++k)
	    {
	      OutputFile.seekg (FileShift + (TmpResumeMinPos * sizeof(double)), ios::beg);
	      ReadLittleEndian(OutputFile, RefCoefficient);
	      double Coefficient = 0.0;
	      if (TmpNbrComputedComponentArray[k] >= 0)
		{
		  for (int j = 0; j < TmpNbrComputedComponentArray[k]; ++j)
		    {
		      long TmpIndex = TmpIndexArray[k][j];
		      if (TmpIndex < this->LargeHilbertSpaceDimension)
			{		  
			  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
			  ReadLittleEndian (OutputFile, TmpComponent);
			  Coefficient += TmpComponentArray[k][j] * TmpComponent;
			}	      	    
		    }		  
		  Coefficient *= InvAlpha;
		  Coefficient /= (RhoRoot - TmpRhoArray[k]);
		}
	      else
		{
		  long TmpIndex = TmpIndexArray[k][0];
		  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
		  ReadLittleEndian (OutputFile, Coefficient);
		}
	      if (Coefficient != RefCoefficient)
		{
		  cout << "error, invalid Jack : component " << TmpResumeMinPos << " is " << RefCoefficient << ", should be " << Coefficient << endl;
		  OutputFile.close();
		  return;
		}
	      ++TmpResumeMinPos;
	    }
	  TmpResumeMinPos = TmpResumePos - memory;
	  if (TmpResumeMinPos < 0l)
	    TmpResumeMinPos = 0l;
	  BufferGlobalIndex = TmpResumeMinPos;
	  OutputFile.seekg ((TmpResumeMinPos * sizeof(double)) + FileShift, ios::beg);
	  double TmpComponent;
	  for (; TmpResumeMinPos < TmpResumePos; ++TmpResumeMinPos)
	    {	      
	      ReadLittleEndian (OutputFile, TmpComponent);
	      TmpVectorBuffer[TmpResumeMinPos % memory] = TmpComponent;	      
	    }
	  OutputFile.close();
	}
      cout << "consistency check done, resuming calculation now" << endl;
      minIndex = TmpResumePos;
    }

  fstream OutputFile;
  OutputFile.open(partialSave, ios::in | ios::binary | ios::out);

  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);


  for (long i = minIndex; i <= maxIndex;)
    {
      long TmpMaxIndex = i + NbrBlocks - 1l;
      long LimNbrBlocks = NbrBlocks;
      if (TmpMaxIndex > maxIndex)
	{
	  LimNbrBlocks = NbrBlocks - (TmpMaxIndex - maxIndex);
	  TmpMaxIndex = maxIndex;
	}
      Operation.SetIndicesRange(i, LimNbrBlocks);
      Operation.ApplyOperation(architecture);
 
      for (long k = 0l; k < LimNbrBlocks; ++k)
	{
	  if (TmpNbrComputedComponentArray[k] >= 0)
	    {
	      double Coefficient = 0.0;
	      for (int j = 0; j < TmpNbrComputedComponentArray[k]; ++j)
		{
		  long TmpIndex = TmpIndexArray[k][j];
		  if (TmpIndex < this->LargeHilbertSpaceDimension)
		    {		  
		      if (TmpIndex < BufferGlobalIndex)
			{
			  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
			  ReadLittleEndian (OutputFile, TmpComponent);
			}
		      else
			{
			  TmpComponent = TmpVectorBuffer[TmpIndex % memory];
			}
		      Coefficient += TmpComponentArray[k][j] * TmpComponent;
		    }	      	    
		}
 
	      Coefficient *= InvAlpha;
	      Coefficient /= (RhoRoot - TmpRhoArray[k]);
	      OutputFile.seekg (0, ios::end);
	      WriteLittleEndian(OutputFile, Coefficient);
	      if (i >= memory)
		++BufferGlobalIndex;
	      TmpVectorBuffer[i % memory] = Coefficient;
	      ++i;
	    }
	  else
	    {
	      long TmpIndex = TmpIndexArray[k][0];
	      if (TmpIndex < BufferGlobalIndex)
		{
		  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
		  ReadLittleEndian (OutputFile, TmpComponent);
		}
	      else
		{
		  TmpComponent = TmpVectorBuffer[TmpIndex % memory];
		}
	      OutputFile.seekg (0, ios::end);
	      WriteLittleEndian(OutputFile, TmpComponent); 	  
	      if (i >= memory)
		++BufferGlobalIndex;
	      TmpVectorBuffer[i % memory] = TmpComponent;
	      ++i;
	    }
	}
      if ((i & DisplayStep) == 0l)
      	{
     	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
      	  cout.flush();
      	}
    }
  OutputFile.close();
  delete[] TmpStateArray;
  delete[] TmpIndexArray;
  delete[] TmpComponentArray;
  cout << endl;
}


// create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry and using sparse storage
//
// alpha = value of the Jack polynomial alpha coefficient
// architecture = architecture to use for precalculation
// partialSave = save partial results in a given vector file
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// memory = amount of memory (in bytes) allowed for temporary vector storage (0 if the whole vector has to be stored in memory)
// memoryBlock = amount of memory (in bytes) allowed for precomputing state indices
// resumeFlag = true if the calculation has to be resumed from a previous one (assuming partialSave contains already computed components)

void FermionOnSphereHaldaneHugeBasis::GenerateSymmetrizedJackPolynomialSparse(double alpha, AbstractArchitecture* architecture, char* partialSave, long minIndex, long maxIndex, long memory, long memoryBlock, bool resumeFlag)
{
  if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension - 1l;
  double TmpComponent = 1.0;
  long FileShift = 4l;
  if (this->HilbertSpaceDimension <= 0)
    FileShift = 12l;
  if ((minIndex <= 0) && (resumeFlag == false))
    {
      ofstream File;
      File.open(partialSave, ios::binary | ios::out);
      WriteLittleEndian(File, this->HilbertSpaceDimension);  
      if (this->HilbertSpaceDimension <= 0)
	{
	  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);  
	}
      WriteLittleEndian(File, TmpComponent);  
      File.close();
    }
  TmpComponent = 0.0;
   
  memory >>= 3;
  if ((memory > this->LargeHilbertSpaceDimension) || (memory <= 0))
    {
      cout << "vector does not require temporary disk storage" << endl;
      memory = this->LargeHilbertSpaceDimension;
    }
  double* TmpVectorBuffer = new double [memory];
  long BufferGlobalIndex = 0l;

  double InvAlpha =  2.0 * (1.0 - alpha) / alpha;

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->GetStateFactorized(0l);;
  this->ConvertToMonomial(MaxRoot, this->TemporaryMonomial);
  for (int j = 0; j < this->NbrFermions; ++j)
    RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - InvAlpha * ((double) j));
  int ReducedNbrFermions = this->NbrFermions - 1;  
  double SymSign = 1.0;
  if ((((this->NbrFermions * ReducedNbrFermions) >> 1) & 1) != 0)
    SymSign = -1.0;

  if (minIndex <= 0)
    minIndex = 1;
  TmpVectorBuffer[0l] = 1.0;
  int MaxArraySize = ((this->NbrFermions * (this->NbrFermions - 1)) / 2) * (this->LzMax + 1);
  long NbrBlocks =  memoryBlock / (((2l* sizeof (long)) + (sizeof(double))) * MaxArraySize);
  if (NbrBlocks == 0)
    {
      NbrBlocks = 100000;
    }
  if (NbrBlocks > this->LargeHilbertSpaceDimension)
    {
      NbrBlocks = (this->LargeHilbertSpaceDimension >> 1) + 1;
    }
  cout << "number of precalculation blocks = " << NbrBlocks << endl;
  long DisplayStep = (this->LargeHilbertSpaceDimension / (1000 * NbrBlocks)) * NbrBlocks;
  long** TmpIndexArray = new long* [NbrBlocks];
  double** TmpComponentArray = new double* [NbrBlocks];
  unsigned long** TmpStateArray = new unsigned long* [NbrBlocks]; 
  double* TmpRhoArray = new double [NbrBlocks];
  int* TmpNbrComputedComponentArray = new int [NbrBlocks];
  for (int j = 0; j < NbrBlocks; ++j)
    {
      TmpIndexArray[j] = new long [MaxArraySize];
      TmpComponentArray[j] = new double [MaxArraySize];
      TmpStateArray[j] = new unsigned long [MaxArraySize];
    }
  FQHESphereJackGeneratorOperation Operation(this, InvAlpha, MaxRoot, TmpIndexArray, TmpStateArray, TmpComponentArray, TmpRhoArray, TmpNbrComputedComponentArray, true, true);

  if (resumeFlag == true)
    {
      ifstream File;
      File.open(partialSave, ios::binary | ios::in);
      File.seekg (0, ios::end);
      long TmpResumePos = File.tellg();
      File.close();
      TmpResumePos -= FileShift;
      TmpResumePos /= sizeof(double); 	      
      long TmpResumeMinPos = TmpResumePos - NbrBlocks;
      long LimNbrBlocks = NbrBlocks;
      if (TmpResumeMinPos < 0l)
	{
	  TmpResumeMinPos = 0l;
	  LimNbrBlocks = TmpResumePos - TmpResumeMinPos + 1;
	}
      long TmpMaxIndex = TmpResumeMinPos + NbrBlocks - 1l;
      if (TmpMaxIndex > maxIndex)
	{
	  LimNbrBlocks = NbrBlocks - (TmpMaxIndex - maxIndex);
	  TmpMaxIndex = maxIndex;
	}
      if (LimNbrBlocks > 0)
	{
	  cout << "consistency check, " << TmpResumePos << " components have already been computed, checking the last " << LimNbrBlocks << " ones" << endl;      
	  Operation.SetIndicesRange(TmpResumeMinPos, LimNbrBlocks);
	  Operation.ApplyOperation(architecture);
	  ifstream OutputFile;
	  OutputFile.open(partialSave, ios::binary | ios::in);
	  double RefCoefficient = 0.0;

	  for (long k = 0l; k < LimNbrBlocks; ++k)
	    {
	      OutputFile.seekg (FileShift + (TmpResumeMinPos * sizeof(double)), ios::beg);
	      ReadLittleEndian(OutputFile, RefCoefficient);
	      double Coefficient = 0.0;
	      if (TmpNbrComputedComponentArray[k] >= 0)
		{
		  for (int j = 0; j < TmpNbrComputedComponentArray[k]; ++j)
		    {
		      long TmpIndex = TmpIndexArray[k][j];
		      if (TmpIndex < this->LargeHilbertSpaceDimension)
			{		  
			  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
			  ReadLittleEndian (OutputFile, TmpComponent);
			  Coefficient += TmpComponentArray[k][j] * TmpComponent;
			}	      	    
		    }		  
		  Coefficient *= InvAlpha;
		  Coefficient /= (RhoRoot - TmpRhoArray[k]);
		}
	      else
		{
		  long TmpIndex = TmpIndexArray[k][0];
		  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
		  ReadLittleEndian (OutputFile, Coefficient);
		  Coefficient *= SymSign;
		}
	      if (Coefficient != RefCoefficient)
		{
		  cout << "error, invalid Jack : component " << TmpResumeMinPos << " is " << RefCoefficient << ", should be " << Coefficient << endl;
		  OutputFile.close();
		  return;
		}
	      ++TmpResumeMinPos;
	    }
	  TmpResumeMinPos = TmpResumePos - memory;
	  if (TmpResumeMinPos < 0l)
	    TmpResumeMinPos = 0l;
	  BufferGlobalIndex = TmpResumeMinPos;
	  OutputFile.seekg ((TmpResumeMinPos * sizeof(double)) + FileShift, ios::beg);
	  double TmpComponent;
	  for (; TmpResumeMinPos < TmpResumePos; ++TmpResumeMinPos)
	    {	      
	      ReadLittleEndian (OutputFile, TmpComponent);
	      TmpVectorBuffer[TmpResumeMinPos % memory] = TmpComponent;	      
	    }
	  OutputFile.close();
	}
      cout << "consistency check done, resuming calculation now" << endl;
      minIndex = TmpResumePos;
    }

  fstream OutputFile;
  OutputFile.open(partialSave, ios::in | ios::binary | ios::out);

  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);


  for (long i = minIndex; i <= maxIndex;)
    {
      long TmpMaxIndex = i + NbrBlocks - 1l;
      long LimNbrBlocks = NbrBlocks;
      if (TmpMaxIndex > maxIndex)
	{
	  LimNbrBlocks = NbrBlocks - (TmpMaxIndex - maxIndex);
	  TmpMaxIndex = maxIndex;
	}
      Operation.SetIndicesRange(i, LimNbrBlocks);
      Operation.ApplyOperation(architecture);
 
      for (long k = 0l; k < LimNbrBlocks; ++k)
	{
	  if (TmpNbrComputedComponentArray[k] >= 0)
	    {
	      double Coefficient = 0.0;
	      for (int j = 0; j < TmpNbrComputedComponentArray[k]; ++j)
		{
		  long TmpIndex = TmpIndexArray[k][j];
		  if (TmpIndex < this->LargeHilbertSpaceDimension)
		    {		  
		      if (TmpIndex < BufferGlobalIndex)
			{
			  long TmpIndex2 = this->FindStateIndexFactorized(this->GetSymmetricState(TmpStateArray[k][j]));
			  if ((TmpIndex2 >= BufferGlobalIndex) && (TmpIndex2 < i))
			    {
			      TmpComponent = TmpVectorBuffer[TmpIndex2 % memory];
			      TmpComponent *= SymSign;
			    }
			  else
			    {
			      OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
			      ReadLittleEndian (OutputFile, TmpComponent);
			    }
			}
		      else
			{
			  TmpComponent = TmpVectorBuffer[TmpIndex % memory];
			}
		      Coefficient += TmpComponentArray[k][j] * TmpComponent;
		    }	      	    
		}
 
	      Coefficient *= InvAlpha;
	      Coefficient /= (RhoRoot - TmpRhoArray[k]);
	      OutputFile.seekg (0, ios::end);
	      WriteLittleEndian(OutputFile, Coefficient);
	      if (i >= memory)
		++BufferGlobalIndex;
	      TmpVectorBuffer[i % memory] = Coefficient;
	      ++i;
	    }
	  else
	    {
	      long TmpIndex = TmpIndexArray[k][0];
	      if (TmpIndex < BufferGlobalIndex)
		{
		  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
		  ReadLittleEndian (OutputFile, TmpComponent);
		}
	      else
		{
		  TmpComponent = TmpVectorBuffer[TmpIndex % memory];
		}
	      TmpComponent *= SymSign;
	      OutputFile.seekg (0, ios::end);
	      WriteLittleEndian(OutputFile, TmpComponent); 	  
	      if (i >= memory)
		++BufferGlobalIndex;
	      TmpVectorBuffer[i % memory] = TmpComponent;
	      ++i;
	    }
	}
      if ((i & DisplayStep) == 0l)
      	{
     	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
      	  cout.flush();
      	}
    }
  OutputFile.close();
  delete[] TmpStateArray;
  delete[] TmpIndexArray;
  delete[] TmpComponentArray;
  cout << endl;
}

// core part of the Jack generator using the Lz<->-Lz symmetry and the factorized algorithm
//
// invAlpha = inverse of the Jack polynomial alpha coefficient
// maxRoot = root partition (in fermionic binary representation)
// partialSave = save partial results in a given vector file
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// indexArray = array where state indices are stored
// stateArray = array use to store computed state description
// componentArray = array where computed component numerical factors are stored
// nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
// rhoArray = rho factor associated to each state

void FermionOnSphereHaldaneHugeBasis::GenerateSymmetrizedJackPolynomialFactorizedCore(double invAlpha, unsigned long maxRoot, long minIndex, long maxIndex, unsigned long** stateArray, double** componentArray, long** indexArray, int* nbrComputedComponents, double* rhoArray)
{
  int ReducedNbrFermions = this->NbrFermions - 1;
  for (long i = minIndex; i <= maxIndex; ++i)
    {
      unsigned long* TmpStateArray = stateArray[i - minIndex];
      double* TmpComponentArray = componentArray[i - minIndex];
      long* TmpIndexArray = indexArray[i - minIndex];
      double Rho = 0.0;
      unsigned long CurrentPartition = 0x0ul;
      CurrentPartition = this->GetStateFactorized(i);
      unsigned long TmpSymState = this->GetSymmetricState(CurrentPartition);
      if (TmpSymState <= CurrentPartition)
	{
	  this->ConvertToMonomial(CurrentPartition, this->TemporaryMonomial);
	  for (int j = 0; j < this->NbrFermions; ++j)
	    Rho += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - invAlpha * ((double) j));
	  rhoArray[i - minIndex] = Rho;
	  int NbrComputedComponents = 0;
	  for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
	    for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
	      {
		double Diff = (double) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
		unsigned int Max = this->TemporaryMonomial[j2];
		unsigned long TmpState = 0x0ul;
		int Tmpj1 = j1;
		int Tmpj2 = j2;
		for (int l = 0; l < this->NbrFermions; ++l)
		  this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
		double Sign = 1.0;
		for (unsigned int k = 1; (k <= Max) && (TmpState < maxRoot); ++k)
		  {
		    ++this->TemporaryMonomial2[Tmpj1];
		    --this->TemporaryMonomial2[Tmpj2];
		    while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] >= this->TemporaryMonomial2[Tmpj1 - 1]))
		      {
			unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
			this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
			this->TemporaryMonomial2[Tmpj1] = Tmp;
			--Tmpj1;
			Sign *= -1.0; 
		      }
		    while ((Tmpj2 < ReducedNbrFermions) && (this->TemporaryMonomial2[Tmpj2] <= this->TemporaryMonomial2[Tmpj2 + 1]))
		      {
			unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
			this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
			this->TemporaryMonomial2[Tmpj2] = Tmp;
			++Tmpj2;
			Sign *= -1.0; 
		      }
		    if ((this->TemporaryMonomial2[Tmpj1] != this->TemporaryMonomial2[Tmpj1 + 1]) && (this->TemporaryMonomial2[Tmpj2] != this->TemporaryMonomial2[Tmpj2 - 1]))
		      {
			TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
			if ((TmpState <= maxRoot) && (TmpState > CurrentPartition))
			  {
			    TmpComponentArray[NbrComputedComponents] = Sign * Diff;
			    TmpStateArray[NbrComputedComponents] = TmpState;
			    ++NbrComputedComponents;
			  }
		      }
		  }
	      }
	  nbrComputedComponents[i - minIndex] = NbrComputedComponents;
	  for (int j = 0; j < NbrComputedComponents; ++j)
	    {
	      TmpIndexArray[j] = this->FindStateIndexFactorized(TmpStateArray[j]);
	    }
	}
      else
	{
	  nbrComputedComponents[i - minIndex] = -1;
	   TmpIndexArray[0] = this->FindStateIndexFactorized(TmpSymState);
	}
    }
}

// core part of the Jack generator using the factorized algorithm
//
// invAlpha = inverse of the Jack polynomial alpha coefficient
// maxRoot = root partition (in fermionic binary representation)
// partialSave = save partial results in a given vector file
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// indexArray = array where state indices are stored
// stateArray = array use to store computed state description
// componentArray = array where computed component numerical factors are stored
// nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
// rhoArray = rho factor associated to each state

void FermionOnSphereHaldaneHugeBasis::GenerateJackPolynomialFactorizedCore(double invAlpha, unsigned long maxRoot, long minIndex, long maxIndex, unsigned long** stateArray, double** componentArray, long** indexArray, int* nbrComputedComponents, double* rhoArray)
{
  int ReducedNbrFermions = this->NbrFermions - 1;
  for (long i = minIndex; i <= maxIndex; ++i)
    {
      unsigned long* TmpStateArray = stateArray[i - minIndex];
      double* TmpComponentArray = componentArray[i - minIndex];
      long* TmpIndexArray = indexArray[i - minIndex];
      double Rho = 0.0;
      unsigned long CurrentPartition = 0x0ul;
      CurrentPartition = this->GetStateFactorized(i);
      int NbrComputedComponents = 0;
      this->ConvertToMonomial(CurrentPartition, this->TemporaryMonomial);
      for (int j = 0; j < this->NbrFermions; ++j)
	Rho += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - invAlpha * ((double) j));
      rhoArray[i - minIndex] = Rho;
      for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
	  {
	    double Diff = (double) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
	    unsigned int Max = this->TemporaryMonomial[j2];
	    unsigned long TmpState = 0x0ul;
	    int Tmpj1 = j1;
	    int Tmpj2 = j2;
	    for (int l = 0; l < this->NbrFermions; ++l)
	      this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
	    double Sign = 1.0;
	    for (unsigned int k = 1; (k <= Max) && (TmpState < maxRoot); ++k)
	      {
		++this->TemporaryMonomial2[Tmpj1];
		--this->TemporaryMonomial2[Tmpj2];
		while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] >= this->TemporaryMonomial2[Tmpj1 - 1]))
		  {
		    unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
		    this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
		    this->TemporaryMonomial2[Tmpj1] = Tmp;
		    --Tmpj1;
		    Sign *= -1.0; 
		  }
                while ((Tmpj2 < ReducedNbrFermions) && (this->TemporaryMonomial2[Tmpj2] <= this->TemporaryMonomial2[Tmpj2 + 1]))
                  {
                    unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
                    this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
                    this->TemporaryMonomial2[Tmpj2] = Tmp;
                    ++Tmpj2;
 		    Sign *= -1.0; 
                 }
		if ((this->TemporaryMonomial2[Tmpj1] != this->TemporaryMonomial2[Tmpj1 + 1]) && (this->TemporaryMonomial2[Tmpj2] != this->TemporaryMonomial2[Tmpj2 - 1]))
		  {
		    TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
		    if ((TmpState <= maxRoot) && (TmpState > CurrentPartition))
		      {
			long TmpIndex = this->FindStateIndexFactorized(TmpState);
			if (TmpIndex < this->LargeHilbertSpaceDimension)
			  {
			    TmpComponentArray[NbrComputedComponents] = Sign * Diff;
			    TmpStateArray[NbrComputedComponents] = TmpState;
			    ++NbrComputedComponents;
			  }
		      }
		  }
	      }
	  }
      nbrComputedComponents[i - minIndex] = NbrComputedComponents;
      for (int j = 0; j < NbrComputedComponents; ++j)
	{
	  TmpIndexArray[j] = this->FindStateIndexFactorized(TmpStateArray[j]);
	}
    }
}

// create the Jack polynomial decomposition corresponding to the root partition, using an optimized version of the code
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// partialSave = save partial results in a given vector file
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& FermionOnSphereHaldaneHugeBasis::OptimizedGenerateJackPolynomial(RealVector& jack, double alpha, long minIndex, long maxIndex, char* partialSave)
{
  jack[0] = 1.0;
  double InvAlpha =  2.0 * (1.0 - alpha) / alpha;

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->TemporaryMonomial);
  for (int j = 0; j < this->NbrFermions; ++j)
    RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - InvAlpha * ((double) j));

  long ConnectedPartitionsBufferSize = 400;
  int* ConnectedPartitionsBufferIndices = new int [ConnectedPartitionsBufferSize];
  int* NbrConnectedPartitions = new int [ConnectedPartitionsBufferSize];
  long** ConnectedPartitionIndices = new long* [ConnectedPartitionsBufferSize];
  double** Factors = new double* [ConnectedPartitionsBufferSize];
  int MaxConnected = (((this->NbrFermions - 1) * this->NbrFermions) / 2) * this->LzMax;
  for (long i = 0; i < ConnectedPartitionsBufferSize; ++i)
    {
      ConnectedPartitionIndices[i] = new long[MaxConnected];
      Factors[i] = new double[MaxConnected];
    }

  if (minIndex <= 0)
    minIndex = 1;
  if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension - 1l;
  for (long i = minIndex; i <= maxIndex; ++i)
    {
      double Rho = 0.0;
      long Limit = i;
      long NbrIndices = 0;
      while ((NbrIndices < ConnectedPartitionsBufferSize) && (Limit < this->LargeHilbertSpaceDimension) && (jack[Limit]))
	{
	  if (jack[Limit] == 0.0)
	    {
	      ConnectedPartitionsBufferIndices[NbrIndices] = Limit;
	      ++NbrIndices;
	    }
	  ++Limit;
	}
      if (NbrIndices > 0)
	{
	  int CurrentIndex = 0;	  
	  for (long k = i; k < Limit; ++k)
	    {
	      double& Coefficient = jack[k];
	      if (Coefficient == 0.0)
		{
		  unsigned long CurrentPartition = this->StateDescription[k];
		  this->ConvertToMonomial(CurrentPartition, this->TemporaryMonomial);
		  for (int j = 0; j < this->NbrFermions; ++j)
		    Rho += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - InvAlpha * ((double) j));
		  double Coefficient2 = 0.0;
		  int TmpNbrConnectedPartitions = NbrConnectedPartitions[CurrentIndex];
		  long* TmpConnectedPartitionIndices = ConnectedPartitionIndices[CurrentIndex];
		  double* TmpFactors = Factors[CurrentIndex];
		  for (int j = 0; j < TmpNbrConnectedPartitions; ++j)
		    Coefficient2 += TmpFactors[j] * jack[TmpConnectedPartitionIndices[j]];
		  Coefficient = Coefficient2;
		  ++CurrentIndex;
		}
	    }
	}
      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
      cout.flush();
      i = Limit;
    }
  cout << endl;

  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry, using an optimized version of the code
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// partialSave = save partial results in a given vector file
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& FermionOnSphereHaldaneHugeBasis::OptimizedGenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha, long minIndex, long maxIndex, char* partialSave)
{
  return jack;
}

// find squeezed partitions that are connected throught the Jack calculation algorithm
//
// nbrPartitions = number of partitions whose connection has to be computed
// partitionIndices = indices of the partitions whose connection has to be computed
// nbrConnectedPartitions = array where the number of partitions connected to a given one will be stored
// connectedPartitionIndices = array where the index of the partitions connected to a given one will be stored
// factors = numerical factor that relates two connected partitions
// rootPartition = Jack root partition

void FermionOnSphereHaldaneHugeBasis::GetConnectedSqueezedPartitions(long nbrPartitions, long* partitionIndices, int* nbrConnectedPartitions, 
								     long** connectedPartitionIndices, double** factors, unsigned long rootPartition)
{
  int ReducedNbrFermions = this->NbrFermions - 1;  

  for (long i = 0; i < nbrPartitions; ++i)
    {
      long* TmpConnectedPartitionIndices = connectedPartitionIndices[i];
      double* TmpFactors = factors[i];
      int CurrentIndex = 0;
      unsigned long CurrentPartition = this->StateDescription[partitionIndices[i]];
      this->ConvertToMonomial(CurrentPartition, this->TemporaryMonomial);
      for (int j1 = 0; j1 < ReducedNbrFermions; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrFermions; ++j2)
	  {
	    double Diff = (double) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
	    unsigned int Max = this->TemporaryMonomial[j2];
	    unsigned long TmpState = 0x0ul;
	    int Tmpj1 = j1;
	    int Tmpj2 = j2;
	    for (int l = 0; l < this->NbrFermions; ++l)
	      this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
	    double Sign = 1.0;
	    for (unsigned int k = 1; (k <= Max) && (TmpState < rootPartition); ++k)
	      {
		++this->TemporaryMonomial2[Tmpj1];
		--this->TemporaryMonomial2[Tmpj2];
		while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] >= this->TemporaryMonomial2[Tmpj1 - 1]))
		  {
		    unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
		    this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
		    this->TemporaryMonomial2[Tmpj1] = Tmp;
		    --Tmpj1;
		    Sign *= -1.0; 
		  }
                while ((Tmpj2 < ReducedNbrFermions) && (this->TemporaryMonomial2[Tmpj2] <= this->TemporaryMonomial2[Tmpj2 + 1]))
                  {
                    unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
                    this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
                    this->TemporaryMonomial2[Tmpj2] = Tmp;
                    ++Tmpj2;
 		    Sign *= -1.0; 
                 }
		if ((this->TemporaryMonomial2[Tmpj1] != this->TemporaryMonomial2[Tmpj1 + 1]) && (this->TemporaryMonomial2[Tmpj2] != this->TemporaryMonomial2[Tmpj2 - 1]))
		  {
		    TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
		    if ((TmpState <= rootPartition) && (TmpState > CurrentPartition))
		      {
			long TmpIndex = this->FindStateIndexMemory(TmpState, this->TemporaryMonomial2[0]);
			if (TmpIndex < this->LargeHilbertSpaceDimension)
			  {
			    TmpConnectedPartitionIndices[CurrentIndex] = TmpIndex;
			    TmpFactors[CurrentIndex] = Sign * Diff;
			    ++CurrentIndex;
			  }
		      }
		  }
	      }
	  }
      nbrConnectedPartitions[i] = CurrentIndex;
    }
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

double FermionOnSphereHaldaneHugeBasis::JackSqrNormalization (RealVector& outputVector, long minIndex, long nbrComponents)
{
  double SqrNorm = 0.0;
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = this->LzMax >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  if (this->NbrRootSuffix == 0)
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  Factorial.SetToOne();
	  this->ConvertToMonomial(this->StateDescription[i], this->TemporaryMonomial);
	  for (int k = 0; k < this->NbrFermions; ++k)
	    {
	      if (HalfLzMax < this->TemporaryMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, this->TemporaryMonomial[k]);
	      else
		if (HalfLzMax > this->TemporaryMonomial[k])
		  Factorial.PartialFactorialDivide(this->TemporaryMonomial[k] + 1, HalfLzMax);
	    }	      
	  SqrNorm +=(outputVector[i] * outputVector[i]) * Factorial.GetNumericalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  else
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  Factorial.SetToOne();
	  this->ConvertToMonomial(this->GetStateFactorized(i), this->TemporaryMonomial);
	  for (int k = 0; k < this->NbrFermions; ++k)
	    {
	      if (HalfLzMax < this->TemporaryMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, this->TemporaryMonomial[k]);
	      else
		if (HalfLzMax > this->TemporaryMonomial[k])
		  Factorial.PartialFactorialDivide(this->TemporaryMonomial[k] + 1, HalfLzMax);
	    }	      
	  SqrNorm +=(outputVector[i] * outputVector[i]) * Factorial.GetNumericalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  cout << endl;
  return SqrNorm;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

LongRational FermionOnSphereHaldaneHugeBasis::JackSqrNormalization (LongRationalVector& outputVector, long minIndex, long nbrComponents)
{
  LongRational SqrNorm = 0l;
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = this->LzMax >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  if (this->NbrRootSuffix == 0)
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  Factorial.SetToOne();
	  this->ConvertToMonomial(this->StateDescription[i], this->TemporaryMonomial);
	  for (int k = 0; k < this->NbrFermions; ++k)
	    {
	      if (HalfLzMax < this->TemporaryMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, this->TemporaryMonomial[k]);
	      else
		if (HalfLzMax > this->TemporaryMonomial[k])
		  Factorial.PartialFactorialDivide(this->TemporaryMonomial[k] + 1, HalfLzMax);
	    }	      
	  SqrNorm +=(outputVector[i] * outputVector[i]) * Factorial.GetLongRationalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  else
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  Factorial.SetToOne();
	  this->ConvertToMonomial(this->GetStateFactorized(i), this->TemporaryMonomial);
	  for (int k = 0; k < this->NbrFermions; ++k)
	    {
	      if (HalfLzMax < this->TemporaryMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, this->TemporaryMonomial[k]);
	      else
		if (HalfLzMax > this->TemporaryMonomial[k])
		  Factorial.PartialFactorialDivide(this->TemporaryMonomial[k] + 1, HalfLzMax);
	    }	      
	  SqrNorm +=(outputVector[i] * outputVector[i]) * Factorial.GetLongRationalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  cout << endl;
  return SqrNorm;
}

// compute part of the Jack polynomial scalar product in a given range of indices
//
// state1 = reference on the first unnormalized Jack polynomial
// state2 = reference on the second unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

double FermionOnSphereHaldaneHugeBasis::JackScalarProduct (RealVector& state1, RealVector& state2, long minIndex, long nbrComponents)
{
  double SqrNorm = 0.0;
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = this->LzMax >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  if (this->NbrRootSuffix == 0)
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  Factorial.SetToOne();
	  this->ConvertToMonomial(this->StateDescription[i], this->TemporaryMonomial);
	  for (int k = 0; k < this->NbrFermions; ++k)
	    {
	      if (HalfLzMax < this->TemporaryMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, this->TemporaryMonomial[k]);
	      else
		if (HalfLzMax > this->TemporaryMonomial[k])
		  Factorial.PartialFactorialDivide(this->TemporaryMonomial[k] + 1, HalfLzMax);
	    }	      
	  SqrNorm +=(state1[i] * state2[i]) * Factorial.GetNumericalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  else
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  Factorial.SetToOne();
	  this->ConvertToMonomial(this->GetStateFactorized(i), this->TemporaryMonomial);
	  for (int k = 0; k < this->NbrFermions; ++k)
	    {
	      if (HalfLzMax < this->TemporaryMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, this->TemporaryMonomial[k]);
	      else
		if (HalfLzMax > this->TemporaryMonomial[k])
		  Factorial.PartialFactorialDivide(this->TemporaryMonomial[k] + 1, HalfLzMax);
	    }	      
	  SqrNorm +=(state1[i] * state2[i]) * Factorial.GetNumericalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  cout << endl;
  return SqrNorm;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state1 = reference on the first unnormalized Jack polynomial
// state2 = reference on the second unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

LongRational FermionOnSphereHaldaneHugeBasis::JackScalarProduct (LongRationalVector& state1, LongRationalVector& state2, long minIndex, long nbrComponents)
{
  LongRational SqrNorm = 0l;
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = this->LzMax >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  if (this->NbrRootSuffix == 0)
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  Factorial.SetToOne();
	  this->ConvertToMonomial(this->StateDescription[i], this->TemporaryMonomial);
	  for (int k = 0; k < this->NbrFermions; ++k)
	    {
	      if (HalfLzMax < this->TemporaryMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, this->TemporaryMonomial[k]);
	      else
		if (HalfLzMax > this->TemporaryMonomial[k])
		  Factorial.PartialFactorialDivide(this->TemporaryMonomial[k] + 1, HalfLzMax);
	    }	      
	  SqrNorm +=(state1[i] * state2[i]) * Factorial.GetLongRationalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  else
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  Factorial.SetToOne();
	  this->ConvertToMonomial(this->GetStateFactorized(i), this->TemporaryMonomial);
	  for (int k = 0; k < this->NbrFermions; ++k)
	    {
	      if (HalfLzMax < this->TemporaryMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, this->TemporaryMonomial[k]);
	      else
		if (HalfLzMax > this->TemporaryMonomial[k])
		  Factorial.PartialFactorialDivide(this->TemporaryMonomial[k] + 1, HalfLzMax);
	    }	      
	  SqrNorm +=(state1[i] * state2[i]) * Factorial.GetLongRationalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  cout << endl;
  return SqrNorm;
}

// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the  component
int FermionOnSphereHaldaneHugeBasis::GetLzValue(int j)
{
  return this->TotalLz;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  FermionOnSphereHaldaneHugeBasis::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState)
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
	  RealSymmetricMatrix TmpDensityMatrix(this->LargeHilbertSpaceDimension);
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

  int ShiftedTotalLz = (this->TotalLz + this->NbrFermions * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrFermionSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  if ((ShiftedLzComplementarySector < (NbrFermionsComplementarySector * subsytemSize)) || (ShiftedLzComplementarySector > (NbrFermionsComplementarySector * (this->LzMax))))
    {
      RealSymmetricMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }

  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0.0;
 	  FermionOnSphere TmpHilbertSpace(this->NbrFermions - nbrFermionSector, 2 * ShiftedLzComplementarySector - ((this->NbrFermions - nbrFermionSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  unsigned long  TmpState2 = 0x0;
	  for (int i = 0; i < nbrFermionSector; ++i)
	    TmpState2 |= 0x1ul << i;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << subsytemSize) | TmpState2;
	      long TmpPos = this->FindStateIndexFactorized(TmpState);
	      if (TmpPos != this->LargeHilbertSpaceDimension)
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
 	  FermionOnSphere TmpHilbertSpace(this->NbrFermions - nbrFermionSector, 2 * ShiftedLzComplementarySector - ((this->NbrFermions - nbrFermionSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << subsytemSize;
	      long TmpPos = this->FindStateIndexFactorized(TmpState);
	      if (TmpPos != this->LargeHilbertSpaceDimension)
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
  long MinIndex = 0l;
  long MaxIndex = this->LargeHilbertSpaceDimension - 1l;
  if (nbrFermionSector == 1)
    {
      double TmpValue = 0.0;
      FermionOnSphere TmpHilbertSpace(this->NbrFermions - nbrFermionSector, 2 * ShiftedLzComplementarySector - ((this->NbrFermions - nbrFermionSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << subsytemSize) | (0x1ul << ShiftedLzSector);
	  long TmpPos = this->FindStateIndexFactorized(TmpState);
	  if (TmpPos != this->LargeHilbertSpaceDimension)
	    TmpValue += groundState[TmpPos] * groundState[TmpPos];	
	}
      RealSymmetricMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  if (NbrFermionsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
      FermionOnSphere TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      MinIndex = this->LargeHilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      double TmpValue;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpValue = groundState[MinIndex + i];
	  for (int j = i; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    TmpDensityMatrix.SetMatrixElement(i, j, TmpValue * groundState[MinIndex + j]);
	}
      return TmpDensityMatrix;
    }


  int TmpNbrFermions;
  int TmpTotalLz;
  int TmpIndex;
  FermionOnSphere TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  long* TmpStatePosition = new long [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  FermionOnSphere TmpHilbertSpace(this->NbrFermions - nbrFermionSector, 2 * ShiftedLzComplementarySector - ((this->NbrFermions - nbrFermionSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.StateDescription[MinIndex] << subsytemSize;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.StateDescription[j] | TmpComplementaryState;
	  long TmpPos = this->FindStateIndexFactorized(TmpState);
	  if (TmpPos != this->LargeHilbertSpaceDimension)
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

RealVector& FermionOnSphereHaldaneHugeBasis::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  double* SqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  if (reference >= 0l)
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0 / state[reference];
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	  InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
	}
      unsigned long TmpState = this->StateDescription[reference];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & 1ul) != 0ul)
	  TmpMonomialReference[Index++] = j;
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      TmpMonomial[Index++] = j;
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
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      state[reference] = 1.0;
      delete[] InvSqrtCoefficients;
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
   }
  else
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0;
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	  InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
	}
      unsigned long TmpState = this->StateDescription[0l];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & 1ul) != 0ul)
	  TmpMonomialReference[Index++] = j;
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      TmpMonomial[Index++] = j;
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
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      delete[] InvSqrtCoefficients;
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
   }
  delete[] SqrtCoefficients;
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& FermionOnSphereHaldaneHugeBasis::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  BinomialCoefficients Binomials(this->LzMax);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  if ( this->StateDescription != 0)
    {
      if (reference >= 0l)
	{
	  int* TmpMonomialReference = new int [this->NbrFermions];
	  int* TmpMonomial = new int [this->NbrFermions];
	  double Factor = 1.0;
	  double* SqrtCoefficients = new double [this->LzMax + 1];
	  double* InvSqrtCoefficients = new double [this->LzMax + 1];
	  for (int k = 0; k <= this->LzMax; ++k)
	    {
	      InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	      SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
	    }
	  unsigned long TmpState = this->StateDescription[reference];
	  int Index = 0;
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      TmpMonomialReference[Index++] = j;
	  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      Index = 0;
	      TmpState = this->StateDescription[i];
	      for (int j = this->LzMax; j >= 0; --j)
		if (((TmpState >> j) & 1ul) != 0ul)
		  TmpMonomial[Index++] = j;
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
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}
	      state[i] *= Coefficient;
	    }
	  delete[] TmpMonomialReference;
	  delete[] TmpMonomial;
	  delete[] InvSqrtCoefficients;
	  state /= state.Norm();
	}
      else
	{
	  int* TmpMonomialReference = new int [this->NbrFermions];
	  int* TmpMonomial = new int [this->NbrFermions];
	  double Factor = 1.0;
	  double* SqrtCoefficients = new double [this->LzMax + 1];
	  double* InvSqrtCoefficients = new double [this->LzMax + 1];
	  for (int k = 0; k <= this->LzMax; ++k)
	    {
	      InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	      SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
	    }
	  unsigned long TmpState = this->StateDescription[0l];
	  int Index = 0;
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      TmpMonomialReference[Index++] = j;
	  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      Index = 0;
	      TmpState = this->StateDescription[i];
	      for (int j = this->LzMax; j >= 0; --j)
		if (((TmpState >> j) & 1ul) != 0ul)
		  TmpMonomial[Index++] = j;
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
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}
	      state[i] *= Coefficient;
	    }
	  delete[] TmpMonomialReference;
	  delete[] TmpMonomial;
	  delete[] InvSqrtCoefficients;
	}
    }
  else
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0;
      double* SqrtCoefficients = new double [this->LzMax + 1];
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	  SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
	}
      unsigned long TmpState = this->GetStateFactorized(reference);
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & 1ul) != 0ul)
	  TmpMonomialReference[Index++] = j;
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->GetStateFactorized(i);
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      TmpMonomial[Index++] = j;
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
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
      state /= state.Norm();
    }
  delete[] SqrtCoefficients;
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

RealVector& FermionOnSphereHaldaneHugeBasis::FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
							 ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
							 bool symmetrizedFlag, double coefficient)
{
  FermionOnSphere* LeftSpace = (FermionOnSphere*) leftSpace;
  FermionOnSphere* RightSpace = (FermionOnSphere*) rightSpace;
  int StateShift = RightSpace->LzMax + padding + 1;
  long Count = 0l;
  for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long TmpState1 = LeftSpace->StateDescription[i] << StateShift;
      int TmpLzMax = this->LzMax;
      while ((TmpState1 >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      double Coefficient = coefficient * leftVector[i];
      if (symmetrizedFlag == false)
	{
	  if (this->CheckDiskStorage() == true)
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->StateDescription[j];
		  TmpState2 |= TmpState1;
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  long TmpIndex = this->FindStateIndexFactorized(TmpState2);
		  double& TmpCoef = outputVector[TmpIndex];
		  if (TmpCoef == 0.0)
		    ++Count;
		  TmpCoef = Coefficient2;
		}
	    }
	  else
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->StateDescription[j];
		  TmpState2 |= TmpState1;
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  long TmpIndex = this->FindStateIndexMemory(TmpState2, TmpLzMax);
		  double& TmpCoef = outputVector[TmpIndex];
		  if (TmpCoef == 0.0)
		    ++Count;
		  TmpCoef = Coefficient2;
		}
	    }
	}
      else
	{
	  if (this->CheckDiskStorage() == true)
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->StateDescription[j];
		  TmpState2 |= TmpState1;
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  long TmpIndex = this->FindStateIndexFactorized(TmpState2);
		  double& TmpCoef = outputVector[TmpIndex];
		  if (TmpCoef == 0.0)
		    ++Count;
		  TmpCoef = Coefficient2;
		  unsigned long TmpState3 = this->GetSymmetricState(TmpState2);
		  if (TmpState3 != TmpState2)
		    {
		      TmpIndex = this->FindStateIndexFactorized(TmpState3);
		      double& TmpCoef2 = outputVector[TmpIndex];
		      if (TmpCoef2 == 0.0)
			++Count;
		      TmpCoef2 = Coefficient2;
		    }
		}
	    }
	  else
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->StateDescription[j];
		  TmpState2 |= TmpState1;
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  long TmpIndex = this->FindStateIndexMemory(TmpState2, TmpLzMax);
		  double& TmpCoef = outputVector[TmpIndex];
		  if (TmpCoef == 0.0)
		    ++Count;
		  TmpCoef = Coefficient2;
		  unsigned long TmpState3 = this->GetSymmetricState(TmpState2);
		  if (TmpState3 != TmpState2)
		    {
		      int TmpLzMax2 = this->LzMax;
		      while ((TmpState3 >> TmpLzMax2) == 0x0ul)
			--TmpLzMax2;
		      TmpIndex = this->FindStateIndexMemory(TmpState3, TmpLzMax2);
		      double& TmpCoef2 = outputVector[TmpIndex];
		      if (TmpCoef2 == 0.0)
			++Count;
		      TmpCoef2 = Coefficient2;
		    }
		}
	    }
	}
    }
  cout << "nbr of newly added components : " << Count << endl;
  return outputVector;
}

// normalize Jack with respect to cylinder basis
//
// state = reference to the Jack state to normalize
// aspect = aspect ratio of cylinder
// return value = normalized state

RealVector& FermionOnSphereHaldaneHugeBasis::NormalizeJackToCylinder(RealVector& state, double aspect)
{
  unsigned long TmpState;
  long double Pi_L = 3.14159265358979323846264338328L;
  long double Length = sqrtl((long double)2.0 * Pi_L * (long double)(this->LzMax + 1) * (long double)aspect);
  cout<<"L= "<<Length<<" r= "<<aspect<<endl;
  long double kappa = (long double)2.0 * Pi_L/Length;
  long double Norm = (long double)0.0;

  if (this->StateDescription != 0)
    {
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  TmpState = this->StateDescription[i];
	  long double Sum2MSquare = (long double)0.0;
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      Sum2MSquare += (j - 0.5*LzMax) * (j - 0.5*LzMax);
	  
	  state[i] *= expl((long double)0.5 * kappa * kappa * Sum2MSquare); 
	  
	  Norm += state[i] * state[i];
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  else
    {
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  TmpState = this->GetStateFactorized(i);
	  long double Sum2MSquare = (long double)0.0;
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      Sum2MSquare += (j - 0.5*LzMax) * (j - 0.5*LzMax);
	  
	  state[i] *= expl((long double)0.5 * kappa * kappa * Sum2MSquare); 
	  
	  Norm += state[i] * state[i];
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  cout<<"Norm= "<<Norm<<endl;
  state /= sqrtl(Norm);
 
  return state;
}
