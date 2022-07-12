////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin                    //
//                  including the Haldane squeezing technique                 //
//                  and the  Lz<->-Lz and Sz<->-Sz symmetries                 //
//                                                                            //
//                        last modification : 26/06/2009                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneLzSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/ArrayTools.h"
#include <math.h>
#include <bitset>
#include <algorithm>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;
using std::pair;
using std::map;


#define WANT_LAPACK

#ifdef __LAPACK__
#ifdef WANT_LAPACK
#define  __USE_LAPACK_HERE__
#endif
#endif


// default constructor
//

FermionOnSphereWithSpinHaldaneLzSzSymmetry::FermionOnSphereWithSpinHaldaneLzSzSymmetry()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// lzMax = twice the maximum Lz value reached by a fermion
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
// rootPartitions = array of root partitions describing the squeezed basis
// nbrRootPartitions = number of root partitions
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinHaldaneLzSzSymmetry::FermionOnSphereWithSpinHaldaneLzSzSymmetry (int nbrFermions, int lzMax, bool minusSzParity, bool minusLzParity, 
											int** rootPartitions, int nbrRootPartitions, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
#ifdef __64_BITS__
  if ((this->LzMax & 1) == 0)
    {
      this->InvertShift = 32 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 2;
    }
  else
    {
      this->InvertShift = 32 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
#else
  if ((this->LzMax & 1) != 0)
    {
      this->InvertShift = 16 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
  else
    {
      this->InvertShift = 16 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 1;
    }
#endif

  this->SzParitySign = 1.0;
  if (minusSzParity == true)
    this->SzParitySign = -1.0;
  this->LzParitySign = 1.0;
  if (minusLzParity == true)
    this->LzParitySign = -1.0;
  if (minusLzParity == minusSzParity)
    this->LzSzSameParityFlag = true;
  else
    this->LzSzSameParityFlag = false;

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

  this->LargeHilbertSpaceDimension = (int) this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
										 (this->TotalSpin + this->NbrFermions) >> 1);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						     (this->TotalSpin + this->NbrFermions) >> 1, 0l);
  this->GenerateLookUpTable(memory);
#ifdef  __64_BITS__
  long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 6) + 1;
#else
  long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 5) + 1;
#endif
  this->KeepStateFlag = new unsigned long [ReducedHilbertSpaceDimension];
  for (int i = 0; i < ReducedHilbertSpaceDimension; ++i)
    this->KeepStateFlag[i] = 0x0l;
  int MaxSweeps = (this->NbrFermions * (this->NbrFermions - 1)) >> 1;  
  this->TmpGeneratedStates =  new unsigned long [MaxSweeps * 1000];
  this->TmpGeneratedStatesLzMax = new int [MaxSweeps * 1000];

  for (int j = 0; j < this->NbrRootPartitions; ++j)
    {
      long Memory = 0l;
      int TmpLzMax = 2 * this->LzMax + 1;
      while (((this->RootPartitions[j] >> TmpLzMax) & 0x1ul) == 0x0ul)
	--TmpLzMax;
      int TmpIndex = this->FindStateIndex(this->RootPartitions[j], TmpLzMax);
#ifdef  __64_BITS__
      this->KeepStateFlag[TmpIndex >> 6] |= 0x1l << (TmpIndex & 0x3f);
#else
      this->KeepStateFlag[TmpIndex >> 5] |= 0x1l << (TmpIndex & 0x1f);
#endif
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
  for (int i = 0; i < (2 * this->NbrLzValue); ++i)
    delete[] this->LookUpTable[i];
  delete[] this->LookUpTable;
  unsigned long* TmpStateDescription = new unsigned long [NewHilbertSpaceDimension];
  int* TmpStateHighestBit = new int [NewHilbertSpaceDimension];
  NewHilbertSpaceDimension = 0l;
  int TotalIndex = 0;
#ifdef  __64_BITS__
  if ((this->LargeHilbertSpaceDimension & 0x3fl) != 0)
#else
  if ((this->LargeHilbertSpaceDimension & 0x1fl) != 0)
#endif
    --ReducedHilbertSpaceDimension;
  int TmpHilbertSpaceDimension = 0;
  if (minusSzParity == minusLzParity)
    {
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
		    unsigned long TmpState = this->StateDescription[TotalIndex];
		    if (this->GetCanonicalState(TmpState) == TmpState)
		      {
			this->GetStateSymmetry(TmpState);
			if ((TmpState & FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT) == FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT)
			  {
			    TmpStateDescription[NewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
			    TmpStateHighestBit[NewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
			    ++NewHilbertSpaceDimension;
			  }
			else
			  {
			    unsigned long TmpStateParity = this->StateDescription[TotalIndex];
			    this->GetStateSingletParity(TmpStateParity);
			    if ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusLzParity == false))
				|| (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusLzParity == true)))
			      {
				TmpStateDescription[NewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
				TmpStateHighestBit[NewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
				++NewHilbertSpaceDimension;
			      }
			  }
		      }
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
		  unsigned long TmpState = this->StateDescription[TotalIndex];
		  if (this->GetCanonicalState(TmpState) == TmpState)
		    {
		      this->GetStateSymmetry(TmpState);
		      if ((TmpState & FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT) == FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT)
			{
			  TmpStateDescription[NewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
			  TmpStateHighestBit[NewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
			  ++NewHilbertSpaceDimension;
			}
		      else
			{
			  unsigned long TmpStateParity = this->StateDescription[TotalIndex];
			  this->GetStateSingletParity(TmpStateParity);
			  if ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusLzParity == false))
			      || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusLzParity == true)))
			    {
			      TmpStateDescription[NewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
			      TmpStateHighestBit[NewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
			      ++NewHilbertSpaceDimension;
			    }
			}
		    }
		}
	      ++TotalIndex;
	    }
	}
    }
  else
    {
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
		    unsigned long TmpState = this->StateDescription[TotalIndex];
		    if (this->GetCanonicalState(TmpState) == TmpState)
		      {
			this->GetStateSymmetry(TmpState);
			if (((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) != 0x0ul) && ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) !=  FERMION_SPHERE_SU2_LZSZ_SYMMETRIC_BIT))		      
			  {
			    if ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) == FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT)
			      {
				TmpStateDescription[NewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
				TmpStateHighestBit[NewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
				++NewHilbertSpaceDimension;
			      }
			    else
			      {
				unsigned long TmpStateParity = this->StateDescription[TotalIndex];
				this->GetStateSingletParity(TmpStateParity);
				if ((((TmpState & FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT) == 0x0ul) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusLzParity == false))
												     || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusLzParity == true))))
				    || (((TmpState & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusSzParity == false))
													|| (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusSzParity == true)))))
				  {
				    TmpStateDescription[NewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
				    TmpStateHighestBit[NewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
				    ++NewHilbertSpaceDimension;
				  }
			      }
			  }
		      }
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
		  unsigned long TmpState = this->StateDescription[TotalIndex];
		  if (this->GetCanonicalState(TmpState) == TmpState)
		    {
		      this->GetStateSymmetry(TmpState);
		      if (((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) != 0x0ul) && ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) !=  FERMION_SPHERE_SU2_LZSZ_SYMMETRIC_BIT))		      
			{
			  if ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) == FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT)
			    {
			      TmpStateDescription[NewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
			      TmpStateHighestBit[NewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
			      ++NewHilbertSpaceDimension;
			    }
			  else
			    {
			      unsigned long TmpStateParity = this->StateDescription[TotalIndex];
			      this->GetStateSingletParity(TmpStateParity);
			      if ((((TmpState & FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT) == 0x0ul) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusLzParity == false))
												   || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusLzParity == true))))
				  || (((TmpState & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusSzParity == false))
												      || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusSzParity == true)))))
				{
				  TmpStateDescription[NewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
				  TmpStateHighestBit[NewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
				  ++NewHilbertSpaceDimension;
				}
			    }
			}
		    }
		  
		}
	      ++TotalIndex;
	    }
	}
    }
  delete[] this->StateDescription;
  delete[] this->StateHighestBit;
  delete[] this->KeepStateFlag;
  this->StateDescription = TmpStateDescription;
  this->StateHighestBit = TmpStateHighestBit;
  this->LargeHilbertSpaceDimension = NewHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

  delete[] this->TmpGeneratedStates;
  delete[] this->TmpGeneratedStatesLzMax;

  this->GenerateLookUpTable(memory);
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	this->GetStateSymmetry(this->StateDescription[i]);

  
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

// textureless constructor
// 
// nbrFermions = number of fermions
// lzMax = twice the maximum Lz value reached by a fermion
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
// texturelessRootPartition = root partition describing the squeezed basis, spin texture has to be added on top of it   
// nbrRootPartitions = number of root partitions
// texturelessFlag = flag to indicate textureless squeezing.
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinHaldaneLzSzSymmetry::FermionOnSphereWithSpinHaldaneLzSzSymmetry (int nbrFermions, int lzMax, bool minusSzParity, bool minusLzParity, 
											int** texturelessRootPartition, int nbrRootPartitions, bool texturelessFlag, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  
#ifdef __64_BITS__
  if ((this->LzMax & 1) == 0)
    {
      this->InvertShift = 32 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 2;
    }
  else
    {
      this->InvertShift = 32 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
#else
  if ((this->LzMax & 1) != 0)
    {
      this->InvertShift = 16 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
  else
    {
      this->InvertShift = 16 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 1;
    }
#endif

  this->SzParitySign = 1.0;
  if (minusSzParity == true)
    this->SzParitySign = -1.0;
  this->LzParitySign = 1.0;
  if (minusLzParity == true)
    this->LzParitySign = -1.0;
  if (minusLzParity == minusSzParity)
    this->LzSzSameParityFlag = true;
  else
    this->LzSzSameParityFlag = false;    
  
  this->NbrRootPartitions = nbrRootPartitions;
  this->RootPartitions = new unsigned long [this->NbrRootPartitions];  
  
  for ( int j = 0 ; j < nbrRootPartitions; j++ )
    {
      this->RootPartitions[j] = 0x0ul;
      int TmpTotalLz = 0;
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  if (texturelessRootPartition[j][i] != 0)
	    {
	      switch (texturelessRootPartition[j][i])
		{
		case 1:
		  {
		    this->RootPartitions[j] |= 0x1ul << (2 * i);
		    TmpTotalLz += i;
		  }
		  break;
		case 2:
		  {
		    this->RootPartitions[j] |= 0x2ul << (2 * i);
		    TmpTotalLz += i;
		    TmpTotalLz += i;
		  }
		  break;
		}
	    }
	}
      if ( j == 0 ) 
	{
	  this->TotalLz = TmpTotalLz;
	}
      else
	{
	  if ( this->TotalLz != TmpTotalLz)
	    cout << "warning : root partition " << j << "does not have the same TotalLz as root partition 0" << endl;
	}
    }    

  this->TotalLz = ((this->TotalLz << 1) - (this->LzMax * this->NbrFermions));

  BosonOnSphereShort *Space = new BosonOnSphereShort(this->NbrFermions, 0, this->LzMax);
 
  this->LargeHilbertSpaceDimension = Space->GetTargetHilbertSpaceDimension();
										 
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    {
      this->HilbertSpaceDimension = 0;
    }
  else
    {
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
    }
  this->Flag.Initialize();      
 
#ifdef  __64_BITS__
  long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 6) + 1;
#else
  long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 5) + 1;
#endif
  this->KeepStateFlag = new unsigned long [ReducedHilbertSpaceDimension];
  for (int i = 0; i < ReducedHilbertSpaceDimension; ++i)
    this->KeepStateFlag[i] = 0x0l;
  
  int MaxSweeps = (this->NbrFermions * (this->NbrFermions - 1)) >> 1;  
  this->TmpGeneratedStates =  new unsigned long [MaxSweeps * 1000];
  this->TmpGeneratedStatesLzMax = new int [MaxSweeps * 1000];
  
  for ( int j = 0 ; j < this->NbrRootPartitions ; j++ )
    {
      long Memory = 0l;      
      unsigned long TmpState = 0x0ul;
      int NbrBosonsPlaced = 0;   
      for ( int CurrentLz = 0 ; CurrentLz <= this->LzMax  ; CurrentLz++)
	{
	  if ( ((this->RootPartitions[j] >> (CurrentLz << 1)) & 0x03l ) == 0x01l )
	    {
	      TmpState += 0x01l << (CurrentLz + NbrBosonsPlaced);
	      NbrBosonsPlaced++;
	    }
	  else if ( ((this->RootPartitions[j] >> (CurrentLz << 1)) & 0x03l ) == 0x02l )
	    {
	      TmpState += 0x03l << ( CurrentLz + NbrBosonsPlaced);
	      NbrBosonsPlaced+= 2;
	    }
	}            
      int TmpLzMax = this->LzMax + this->NbrFermions;
      while ( ((TmpState >> TmpLzMax) & 0x01ul)  == 0x0ul ) 
	{
	  TmpLzMax--;
	}  
      int TmpIndex = Space->FindStateIndex(TmpState, TmpLzMax);
      TmpLzMax = 2 * this->LzMax + 1;    
#ifdef  __64_BITS__
      this->KeepStateFlag[TmpIndex >> 6] |= 0x1l << (TmpIndex & 0x3f);
#else
      this->KeepStateFlag[TmpIndex >> 5] |= 0x1l << (TmpIndex & 0x1f);
#endif  
      this->GenerateSqueezedTexturelessStates(TmpLzMax, this->RootPartitions[j], 1, Space, Memory);  
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
  
  unsigned long* TmpStateDescription = new unsigned long [NewHilbertSpaceDimension];
  int* TmpStateHighestBit = new int [NewHilbertSpaceDimension];
  NewHilbertSpaceDimension = 0l;
  int TotalIndex = 0;

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
	      TmpStateDescription[NewHilbertSpaceDimension] =  Space->FermionBasis->StateDescription[TotalIndex];
	      TmpStateHighestBit[NewHilbertSpaceDimension] = Space->FermionBasis->StateLzMax[TotalIndex];
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
	      TmpStateDescription[NewHilbertSpaceDimension] =  Space->FermionBasis->StateDescription[TotalIndex];
	      TmpStateHighestBit[NewHilbertSpaceDimension] = Space->FermionBasis->StateLzMax[TotalIndex];
	      ++NewHilbertSpaceDimension;
	    }
	  ++TotalIndex;
	}
    }
    
  delete Space;  
  delete []this->KeepStateFlag;  
  delete []this->TmpGeneratedStates;
  delete []this->TmpGeneratedStatesLzMax;

  this->LargeHilbertSpaceDimension = (int) this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
										 (this->TotalSpin + this->NbrFermions) >> 1);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1,  (this->TotalSpin + this->NbrFermions) >> 1, 0l);
  this->GenerateLookUpTable(memory);
 
#ifdef  __64_BITS__
  ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 6) + 1;
#else
  ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 5) + 1;
#endif
  this->KeepStateFlag = new unsigned long [ReducedHilbertSpaceDimension];
  for (int i = 0; i < ReducedHilbertSpaceDimension; ++i)
    this->KeepStateFlag[i] = 0x0l;
  
  this->EvaluatePermutations();
  
  int NbrSingles;
  int *Positions = new int[this->NbrFermions];
  unsigned long TmpDressedState, TmpDressedState2;
  int BosonsPassed;
  int NewNewHilbertSpaceDimension = 0;
  //now dress with spins. First set 2s to be Xs (3 or 11) and use the permutations to set the ups (2 or 10)and downs (1 or 01).
  for ( int i = 0 ; i < NewHilbertSpaceDimension ; i++ )
    {
      TmpDressedState = 0x0ul;
      BosonsPassed = 0;
      NbrSingles = 0;
      for ( int Lz = 0 ; Lz <= this->LzMax  ; Lz ++ ) 
	{
	  if  ( ((TmpStateDescription[i] >> (Lz + BosonsPassed)) & 0x03ul) ==  0x03ul ) 
	    {
	      TmpDressedState += (0x03ul << (Lz << 1)); 
	      BosonsPassed += 2;
	    }
	  else if ( ((TmpStateDescription[i] >> (Lz + BosonsPassed)) & 0x03ul) ==  0x01ul ) 
	    {
	      Positions[NbrSingles++] = Lz;
	      BosonsPassed += 1;
	    }
	}   
	
      
      for ( int j = 0 ; j < this->NbrPermutations[NbrSingles >> 1]; j++ )
	{
	  TmpDressedState2 = TmpDressedState;
	  for ( int k = 0 ; k < NbrSingles; k++ )
	    {
	      if ( ((this->Permutations[NbrSingles >> 1][j]  >>  k)  & 0x01ul)  == 0x01ul ) 
		{
		  TmpDressedState2 += (0x02ul << (Positions[k] << 1));
		}
	      else
		{
		  TmpDressedState2 += (0x01ul << (Positions[k]<<1));
		}
	    }
	  int TmpLzMax = 2*this->LzMax + 1;
	  while ( ((TmpDressedState2 >> TmpLzMax ) & 0x1ul) == 0x0ul)
	      --TmpLzMax;
	  int TmpIndex = this->FindStateIndex(TmpDressedState2, TmpLzMax);
#ifdef __64_BITS__
	  if (((this->KeepStateFlag[TmpIndex >> 6] >> (TmpIndex & 0x3f)) & 0x1l) == 0x0l ) 
	    {	      
	      this->KeepStateFlag[TmpIndex >> 6] |= 0x1l << (TmpIndex & 0x3f);
	      NewNewHilbertSpaceDimension++;
	    }
#else
	  if (((this->KeepStateFlag[TmpIndex >> 5] >> (TmpIndex & 0x1f)) & 0x1l) == 0x0l)	    
	    {
	      this->KeepStateFlag[TmpIndex >> 5] |= 0x1l << (TmpIndex & 0x1f);	      
	      NewNewHilbertSpaceDimension++;
	    }      
#endif	  	  	  
	}
    }
  
  
  delete[] this->NbrPermutations;
  for (int i = 0 ; i <= this->NbrFermionsUp; ++i)
    delete[] this->Permutations[i];
  delete[] this->Permutations;
  delete []TmpStateDescription;
  delete []TmpStateHighestBit;
  delete []Positions;    
  
  TmpStateDescription = new unsigned long[NewNewHilbertSpaceDimension];
  TmpStateHighestBit = new int[NewNewHilbertSpaceDimension];
   
  NewNewHilbertSpaceDimension = 0;
  TotalIndex = 0;
  
  #ifdef  __64_BITS__
  if ((this->LargeHilbertSpaceDimension & 0x3fl) != 0)
#else
  if ((this->LargeHilbertSpaceDimension & 0x1fl) != 0)
#endif
    --ReducedHilbertSpaceDimension;
  int TmpHilbertSpaceDimension = 0;
  if (minusSzParity == minusLzParity)
    {
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
		    unsigned long TmpState = this->StateDescription[TotalIndex];
		    if (this->GetCanonicalState(TmpState) == TmpState)
		      {
			this->GetStateSymmetry(TmpState);
			if ((TmpState & FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT) == FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT)
			  {
			    TmpStateDescription[NewNewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
			    TmpStateHighestBit[NewNewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
			    ++NewNewHilbertSpaceDimension;
			  }
			else
			  {
			    unsigned long TmpStateParity = this->StateDescription[TotalIndex];
			    this->GetStateSingletParity(TmpStateParity);
			    if ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusLzParity == false))
				|| (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusLzParity == true)))
			      {
				TmpStateDescription[NewNewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
				TmpStateHighestBit[NewNewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
				++NewNewHilbertSpaceDimension;
			      }
			  }
		      }
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
		  unsigned long TmpState = this->StateDescription[TotalIndex];
		  if (this->GetCanonicalState(TmpState) == TmpState)
		    {
		      this->GetStateSymmetry(TmpState);
		      if ((TmpState & FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT) == FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT)
			{
			  TmpStateDescription[NewNewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
			  TmpStateHighestBit[NewNewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
			  ++NewNewHilbertSpaceDimension;
			}
		      else
			{
			  unsigned long TmpStateParity = this->StateDescription[TotalIndex];
			  this->GetStateSingletParity(TmpStateParity);
			  if ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusLzParity == false))
			      || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusLzParity == true)))
			    {
			      TmpStateDescription[NewNewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
			      TmpStateHighestBit[NewNewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
			      ++NewNewHilbertSpaceDimension;
			    }
			}
		    }
		}
	      ++TotalIndex;
	    }
	}
    }
  else
    {
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
		    unsigned long TmpState = this->StateDescription[TotalIndex];
		    if (this->GetCanonicalState(TmpState) == TmpState)
		      {
			this->GetStateSymmetry(TmpState);
			if (((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) != 0x0ul) && ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) !=  FERMION_SPHERE_SU2_LZSZ_SYMMETRIC_BIT))		      
			  {
			    if ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) == FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT)
			      {
				TmpStateDescription[NewNewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
				TmpStateHighestBit[NewNewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
				++NewNewHilbertSpaceDimension;
			      }
			    else
			      {
				unsigned long TmpStateParity = this->StateDescription[TotalIndex];
				this->GetStateSingletParity(TmpStateParity);
				if ((((TmpState & FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT) == 0x0ul) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusLzParity == false))
												     || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusLzParity == true))))
				    || (((TmpState & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusSzParity == false))
													|| (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusSzParity == true)))))
				  {
				    TmpStateDescription[NewNewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
				    TmpStateHighestBit[NewNewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
				    ++NewNewHilbertSpaceDimension;
				  }
			      }
			  }
		      }
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
		  unsigned long TmpState = this->StateDescription[TotalIndex];
		  if (this->GetCanonicalState(TmpState) == TmpState)
		    {
		      this->GetStateSymmetry(TmpState);
		      if (((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) != 0x0ul) && ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) !=  FERMION_SPHERE_SU2_LZSZ_SYMMETRIC_BIT))		      
			{
			  if ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) == FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT)
			    {
			      TmpStateDescription[NewNewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
			      TmpStateHighestBit[NewNewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
			      ++NewNewHilbertSpaceDimension;
			    }
			  else
			    {
			      unsigned long TmpStateParity = this->StateDescription[TotalIndex];
			      this->GetStateSingletParity(TmpStateParity);
			      if ((((TmpState & FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT) == 0x0ul) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusLzParity == false))
												   || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusLzParity == true))))
				  || (((TmpState & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusSzParity == false))
												      || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusSzParity == true)))))
				{
				  TmpStateDescription[NewNewHilbertSpaceDimension] =  TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK;
				  TmpStateHighestBit[NewNewHilbertSpaceDimension] = this->StateHighestBit[TotalIndex];
				  ++NewNewHilbertSpaceDimension;
				}
			    }
			}
		    }
		  
		}
	      ++TotalIndex;
	    }
	}
    }
    
  delete[] this->SignLookUpTable;
  delete[] this->SignLookUpTableMask;
  delete[] this->LookUpTableShift;
  for (int i = 0; i < (2 * this->NbrLzValue); ++i)
    delete[] this->LookUpTable[i];
  delete[] this->LookUpTable;

  
  delete[] this->StateDescription;
  delete[] this->StateHighestBit;
  delete[] this->KeepStateFlag;
  this->StateDescription = TmpStateDescription;
  this->StateHighestBit = TmpStateHighestBit;
  this->LargeHilbertSpaceDimension = NewNewHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

  this->GenerateLookUpTable(memory);
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	this->GetStateSymmetry(this->StateDescription[i]);
  
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

FermionOnSphereWithSpinHaldaneLzSzSymmetry::FermionOnSphereWithSpinHaldaneLzSzSymmetry(const FermionOnSphereWithSpinHaldaneLzSzSymmetry& fermions)
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
  this->NbrRootPartitions = fermions.NbrRootPartitions;
  this->RootPartitions = new unsigned long [this->NbrRootPartitions];
  for (int i = 0; i < this->NbrRootPartitions; ++i)
    this->RootPartitions[i] = fermions.RootPartitions[i];
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->LzParitySign = fermions.LzParitySign;
  this->SzParitySign = fermions.SzParitySign;
  this->LzSzSameParityFlag = fermions.LzSzSameParityFlag;
}

// destructor
//

FermionOnSphereWithSpinHaldaneLzSzSymmetry::~FermionOnSphereWithSpinHaldaneLzSzSymmetry ()
{
  delete[] this->RootPartitions;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinHaldaneLzSzSymmetry& FermionOnSphereWithSpinHaldaneLzSzSymmetry::operator = (const FermionOnSphereWithSpinHaldaneLzSzSymmetry& fermions)
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
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->LzParitySign = fermions.LzParitySign;
  this->SzParitySign = fermions.SzParitySign;
  this->LzSzSameParityFlag = fermions.LzSzSameParityFlag;
  this->NbrRootPartitions = fermions.NbrRootPartitions;
  delete[] this->RootPartitions;
  this->RootPartitions = new unsigned long [this->NbrRootPartitions];
  for (int i = 0; i < this->NbrRootPartitions; ++i)
    this->RootPartitions[i] = fermions.RootPartitions[i];
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinHaldaneLzSzSymmetry::Clone()
{
  return new FermionOnSphereWithSpinHaldaneLzSzSymmetry(*this);
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereWithSpinHaldaneLzSzSymmetry::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  stateDescription &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = (this->StateDescription[PosMid] & FERMION_SPHERE_SU2_SYMMETRIC_MASK);
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
      CurrentState = (this->StateDescription[PosMid] & FERMION_SPHERE_SU2_SYMMETRIC_MASK);
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}


// generate all squeezed states from a root partition
// 
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinHaldaneLzSzSymmetry::GenerateSqueezedStates(int lzMax, unsigned long referenceState, long pos, long& memory)
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

  int TmpIndex;
  int NbrNewEntries = 0;
  for (int i = 0; i < NbrEntries; ++i)
    {
      unsigned long& TmpState = TmpGeneratedStates2[i];
      TmpIndex = this->FindStateIndex(TmpState, TmpLzMax[i]);
      //      if (TmpIndex >= this->HilbertSpaceDimension)
	//	cout << "error " << TmpIndex << " " << hex  << TmpState << dec << " " << TmpLzMax[i] << endl;
#ifdef __64_BITS__
      if ((this->KeepStateFlag[TmpIndex >> 6] >> (TmpIndex & 0x3f)) & 0x1l)
	{
	  TmpState = 0x0l;
	}
      else
	{
	  this->KeepStateFlag[TmpIndex >> 6] |= 0x1l << (TmpIndex & 0x3f);
	  ++NbrNewEntries;
	}
#else
      if ((this->KeepStateFlag[TmpIndex >> 5] >> (TmpIndex & 0x1f)) & 0x1l)
	{
	  TmpState = 0x0l;
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
	pos = this->GenerateSqueezedStates(TmpLzMax[i], TmpGeneratedStates2[i], pos, memory);

  memory -= 1;
  return pos;
}

// convert a gien state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereWithSpinHaldaneLzSzSymmetry::ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  int TmpIndex;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      Signature = nbodyBasis.StateDescription[i];
      TmpState = this->GetSignedCanonicalState(Signature);
      Signature = TmpState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
      unsigned long TmpState2 = TmpState;
      TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
      NewLzMax = 1 + (this->LzMax << 1);
      while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;	 
      switch (Signature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT)
	{
	case 0x0ul:
	  {
	    if (this->LzSzSameParityFlag == true)
	      {
		Signature = TmpState;
		this->GetStateSingletParity(Signature);
		if (((1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul))) * this->SzParitySign) > 0.0)
		  {
		    TmpIndex = this->FindStateIndex(TmpState, NewLzMax);
		    if ( TmpIndex < this->HilbertSpaceDimension ) 
		      {		    
			TmpVector[i] = state[TmpIndex];
		      }
		  }
	      }
	  }
	  break;
	case FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT :
	  {
	    if (this->LzSzSameParityFlag == true)
	      {
		Signature = TmpState;
		this->GetStateSingletParity(Signature);
		double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
		TmpIndex = this->FindStateIndex(TmpState, NewLzMax);
		if ( TmpIndex < this->HilbertSpaceDimension ) 
		  {		    		  
		    TmpVector[i] = ((1.0 + ((double) (((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) 
						   | (TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT)) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
				* state[TmpIndex] * M_SQRT1_2);
		  }
	      }		
	  }
	  break;
	case FERMION_SPHERE_SU2_SZ_SYMMETRIC_TEST :
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
	    if ((TmpSign * this->LzParitySign) > 0.0)
	      {
		TmpIndex = this->FindStateIndex(TmpState, NewLzMax);
		if ( TmpIndex < this->HilbertSpaceDimension ) 
		  {		    		  
		    TmpVector[i] = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
				* state[TmpIndex] * M_SQRT1_2);
		  }
	      }
	  }
	  break;
	case FERMION_SPHERE_SU2_LZ_SYMMETRIC_TEST :
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
	    if ((TmpSign * this->SzParitySign) > 0.0)
	      {
		TmpIndex = this->FindStateIndex(TmpState, NewLzMax);
		if ( TmpIndex < this->HilbertSpaceDimension ) 
		  {		    		  
		    TmpVector[i] = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->LzParitySign) - 1.0))
				* state[TmpIndex] * M_SQRT1_2);
		  }
	      }
	  }
	  break;
	default:
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
	    TmpIndex = this->FindStateIndex(TmpState, NewLzMax);
	    if ( TmpIndex < this->HilbertSpaceDimension ) 
	      {		    		  
		TmpVector[i] = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->LzParitySign) - 1.0))
			    * (1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
			    * state[TmpIndex] * 0.5);
	      }
	  }
	  break;
	}
    }
  return TmpVector;
}

// convert a given state from the usual n-body basis to the symmetric basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereWithSpinHaldaneLzSzSymmetry::ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis)
{
  RealVector TmpVector (this->HilbertSpaceDimension, true);
//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//     TmpVector[i] = state[nbodyBasis.FindStateIndex(this->StateDescription[i], this->StateHighestBit[i])];
//   TmpVector /= TmpVector.Norm();
  return TmpVector;
}

// generate all squeezed states from a textureless root partition
// 
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinHaldaneLzSzSymmetry::GenerateSqueezedTexturelessStates(int lzMax, unsigned long referenceState, long pos, BosonOnSphereShort *bosonSpace, long& memory)
{
  int MaxSweeps = (this->NbrFermions * (this->NbrFermions - 1)) >> 1;  
  unsigned long* TmpGeneratedStates2 = this->TmpGeneratedStates + (MaxSweeps * memory);
  memory += 1;
  int TmpCurrentLeftLz = (lzMax - 1) >> 1;
  int TmpCurrentRightLz;
  int NbrEntries = 0;
  unsigned long TmpReferenceState;
  unsigned long TmpReferenceState2;  
  //int TwiceLzMax = this->LzMax << 1;
          
  while (TmpCurrentLeftLz >= 2 )
    {
      while ( (TmpCurrentLeftLz >= 2)  &&  ((((referenceState >> (TmpCurrentLeftLz << 1) ) & 0x3l) == 0x0l) || (((referenceState >> ((TmpCurrentLeftLz-1) << 1) ) & 0x2l) != 0x0l)) )
	TmpCurrentLeftLz -= 1;      
      
      if (TmpCurrentLeftLz >= 2) 
	{
	  if  (((referenceState >> (TmpCurrentLeftLz << 1) ) & 0x3l) == 0x1l )
	    {
	      TmpReferenceState = referenceState - (0x1ul << (TmpCurrentLeftLz << 1));	      
	    }
	  else if (((referenceState >> (TmpCurrentLeftLz << 1) ) & 0x3l) == 0x2l ) 
	    {
	      TmpReferenceState = referenceState - (0x2ul << (TmpCurrentLeftLz << 1)) + (0x1ul << (TmpCurrentLeftLz << 1));	      	      
	    }
	  if (((referenceState >> ((TmpCurrentLeftLz-1) << 1) ) & 0x1l) == 0x1l  )
	    {
	      TmpReferenceState = TmpReferenceState - (0x1ul << ((TmpCurrentLeftLz-1) << 1)) + (0x2ul << ((TmpCurrentLeftLz-1) << 1)) ;	      	      
	    }
	  else
	    {
	      TmpReferenceState = TmpReferenceState + (0x1ul << ((TmpCurrentLeftLz-1) << 1)) ;	      	      
	    }	  
	  TmpCurrentRightLz = TmpCurrentLeftLz - 2;
	  while ( TmpCurrentRightLz >= 0 )
	    {
	      while ( (TmpCurrentRightLz >= 0) &&  ((((TmpReferenceState >> (TmpCurrentRightLz << 1) ) & 0x3l) == 0x0l) || (((TmpReferenceState >> ((TmpCurrentRightLz+1) << 1) ) & 0x2l) != 0x0l)) )
		TmpCurrentRightLz -= 1;
	      
	      if (TmpCurrentRightLz >= 0) 
		{
		  if  (((TmpReferenceState >> (TmpCurrentRightLz << 1) ) & 0x3l) == 0x1l )
		    {
		      TmpReferenceState2 = TmpReferenceState - (0x1ul << (TmpCurrentRightLz << 1));	      
		    }
		  else if (((TmpReferenceState >> (TmpCurrentRightLz << 1) ) & 0x3l) == 0x2l ) 
		    {
		      TmpReferenceState2 = TmpReferenceState - (0x2ul << (TmpCurrentRightLz << 1)) + (0x1ul << (TmpCurrentRightLz << 1));	      	      
		    }
		  if (((TmpReferenceState >> ((TmpCurrentRightLz-1) << 1) ) & 0x1l) == 0x1l  )
		    {
		      TmpReferenceState2 = TmpReferenceState2 - (0x1ul << ((TmpCurrentRightLz+1) << 1)) + (0x2ul << ((TmpCurrentRightLz+1) << 1)) ;	      	      
		    }
		  else
		    {
		      TmpReferenceState2 = TmpReferenceState2 + (0x1ul << ((TmpCurrentRightLz+1) << 1)) ;	      	      
		    }
		  TmpGeneratedStates2[NbrEntries] = TmpReferenceState2;
		  ++NbrEntries;
		}
		TmpCurrentRightLz -= 1;
	    }
	  TmpCurrentLeftLz -= 1;
	}
    }	    	    	    	    
            
  int NbrNewEntries =0;
  unsigned long *NewEntries = new unsigned long[NbrEntries];
  for (int i = 0; i < NbrEntries; ++i)
    {
      unsigned long TmpState = 0x0ul;
      
      int NbrBosonsPlaced = 0;
      for ( int CurrentLz = 0 ; CurrentLz <= ((lzMax - 1) >> 1) ; CurrentLz++)
	{
	  if ( ((TmpGeneratedStates2[i] >> (CurrentLz << 1)) & 0x03l ) == 0x01l )
	    {
	      TmpState += 0x01l << (CurrentLz + NbrBosonsPlaced);
	      NbrBosonsPlaced++;
	    }
	  else if ( ((TmpGeneratedStates2[i] >> (CurrentLz << 1)) & 0x03l ) == 0x02l )
	    {
	      TmpState += 0x03l << ( CurrentLz + NbrBosonsPlaced);
	      NbrBosonsPlaced+= 2;
	    }
	}            
      int TmpLzMax = this->LzMax + this->NbrFermions;
      while ( ((TmpState >> TmpLzMax) & 0x01ul)  == 0x0ul ) 
	{
	  TmpLzMax--;
	}
      int TmpIndex = bosonSpace->FindStateIndex(TmpState, TmpLzMax);
#ifdef __64_BITS__
      if (((this->KeepStateFlag[TmpIndex >> 6] >> (TmpIndex & 0x3f)) & 0x1l) == 0x0ul)
	{
	  this->KeepStateFlag[TmpIndex >> 6] |= 0x1l << (TmpIndex & 0x3f);
	  NewEntries[NbrNewEntries] = TmpGeneratedStates2[i];	  
	  ++NbrNewEntries;
	}
#else
      if (((this->KeepStateFlag[TmpIndex >> 5] >> (TmpIndex & 0x1f)) & 0x1l) == 0x0ul)
	{
	  this->KeepStateFlag[TmpIndex >> 5] |= 0x1l << (TmpIndex & 0x1f);
	  NewEntries[NbrNewEntries] = TmpGeneratedStates2[i];
	  ++NbrNewEntries;
	}          
#endif    
    }

  
  for (int i = 0; i < NbrNewEntries; ++i)
    {
      if ( NewEntries[i] != 0x0l)
	{
	  int TmpLzMax = 0;
	  pos = this->GenerateSqueezedTexturelessStates(lzMax, NewEntries[i], pos, bosonSpace, memory);
	}
	memory -= 1;
    }
  delete [] NewEntries;
  return pos;
}

// evaluate all permutations requested to sapply spin texture
//
void FermionOnSphereWithSpinHaldaneLzSzSymmetry::EvaluatePermutations()
{
  BinomialCoefficients Binomial(this->NbrFermions);
  this->NbrPermutations = new int [this->NbrFermionsUp + 1];
  this->Permutations = new unsigned long*[this->NbrFermionsUp + 1];
  Binomial(this->NbrFermions, this->NbrFermionsUp);
  int TmpNbrFermions = 0;
  for (int i = 0; i <= this->NbrFermionsUp; ++i)
    {
      this->NbrPermutations[i] = Binomial(TmpNbrFermions, i); 
      this->Permutations[i] = new unsigned long[this->NbrPermutations[i]];
      unsigned long MinValue = (0x1ul << i) - 0x1ul;
      unsigned long MaxValue = MinValue << (TmpNbrFermions - i);
      unsigned long* TmpPermutations = this->Permutations[i];
      int TmpNbrPermutations = 0;
      for (; MinValue <= MaxValue; ++MinValue)
	{
	  int Count = 0;
	  int Pos = 0;
	  while ((Pos < TmpNbrFermions) && (Count <= i))
	    {
	      if (((MinValue >> Pos) & 0x1ul) != 0x0ul)
		++Count;
	      ++Pos;
	    }
	  if (Count == i)
	    {
	      TmpPermutations[TmpNbrPermutations] = MinValue;
	      ++TmpNbrPermutations;
	    }
	}
      TmpNbrFermions += 2;
    }
  return;
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

void FermionOnSphereWithSpinHaldaneLzSzSymmetry::MonomialsTimesPolarizedSlaterProjection(unsigned long * slater, unsigned long * monomial, map<unsigned long , double> & sortingMap, unsigned long nbrPermutations , unsigned long * permutations1, unsigned long * permutations2, double initialCoef)
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
	  
	  
	  for(int k = 0 ; (k < HalfNbrParticles) && (CoefUp != 0.0); k++)
	    {
	      State[k] = (HalfMonomialsUp[k]>>1) + slater[k];
	      if ((HalfMonomialsUp[k] & 0x1ul) != 0ul)
		{
		  long Numerator = -((HalfMonomialsUp[k]>>1) * TmpFinalLzMaxUp) + (State[k] * TmpLzMaxUp);
		  if (Numerator == 0l)
		    CoefUp = 0.0;
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
		    
		  for(int k = 0 ; (k < HalfNbrParticles) && (CoefDown != 0.0); k++)
		    {
		      State[k+HalfNbrParticles] = (HalfMonomialsDown[k]>>1) + slater[k];
		      if ((HalfMonomialsDown[k] & 0x1ul) != 0ul)
			{
			  long Numerator = -((HalfMonomialsDown[k]>>1) * TmpFinalLzMaxUp) + (State[HalfNbrParticles + k] * TmpLzMaxUp);
			  if (Numerator == 0l)
			    CoefDown = 0.0;
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
		      
		      for (int i = 0; (i < HalfNbrParticles )&& (Bool); i++)
			{
			  Mask = (1ul << ((State[i]<<1) +1));
			  if((TmpState & Mask) != 0)
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
			  if((TmpState & Mask) != 0)
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
			  
			  double TmpCoef = CoefDown*CoefUp;
			  unsigned long NewTmpState = TmpState;
			  //int TmpIndex = this->SymmetrizeAdAdResult(NewTmpState, TmpCoef);
			  int TmpIndex =  this->ConvertToSymmetricNbodyBasis(TmpCoef, NewTmpState);
			  if ( TmpIndex < this->HilbertSpaceDimension )
			    {
			      //cout << TmpState << " -> " << NewTmpState << " (" << TmpIndex << "), Coeff: " << CoefDown*CoefUp << " -> " << TmpCoef << endl;			      
			      //InsertionResult = sortingMap.insert (pair <unsigned long, double> (TmpState , CoefDown*CoefUp));
			      InsertionResult = sortingMap.insert (pair <unsigned long, double> (NewTmpState , TmpCoef));
			      
			      if (InsertionResult.second == false)
				{
				  InsertionResult.first->second += TmpCoef;
				}
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
}

// convert a given coefficient of state from the usual n-body basis to the symmetric basis
//
// coefficient = coefficient of configuration in usual n-body basis
// stateDescription = configuration in usual n-body basis
// return value = index in this basis.

int FermionOnSphereWithSpinHaldaneLzSzSymmetry::ConvertToSymmetricNbodyBasis(double& coefficient, unsigned long& stateDescription)
{
  //RealVector TmpVector (this->GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  
  Signature = stateDescription;
  TmpState = this->GetSignedCanonicalState(Signature);  
  Signature = TmpState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  unsigned long TmpState2 = TmpState;
  TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  NewLzMax = 1 + (this->LzMax << 1);
  while ((TmpState >> NewLzMax) == 0x0ul)
   --NewLzMax;	 
  stateDescription = TmpState;
  //switch (Signature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT)
  switch (Signature & FERMION_SPHERE_SU2_SYMMETRIC_BIT)
    {
    case 0x0ul:
      {
	if (this->LzSzSameParityFlag == true)
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    if (((1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul))) * this->SzParitySign) > 0.0)
	      {		
		return this->FindStateIndex(TmpState, NewLzMax);		
	      }
	  }
      }
      break;
    case FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT :
      {
	if (this->LzSzSameParityFlag == true)
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));	    
// 	    TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) (((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) 
// 											| (TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT)) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
// 								    * state[i] * M_SQRT1_2);
	    coefficient = ((1.0 + ((double) (((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) 
											| (TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT)) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
								    * coefficient * M_SQRT1_2);
	    return this->FindStateIndex(TmpState, NewLzMax);
	  }		
      }
      break;
    case FERMION_SPHERE_SU2_SZ_SYMMETRIC_TEST :
      {
	Signature = TmpState;
	this->GetStateSingletParity(Signature);
	double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
	if ((TmpSign * this->LzParitySign) > 0.0)
	  {
// 	    TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
// 								    * state[i] * M_SQRT1_2);
	    coefficient = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
								    * coefficient * M_SQRT1_2);
	    return this->FindStateIndex(TmpState, NewLzMax);
	  }
      }
      break;
    case FERMION_SPHERE_SU2_LZ_SYMMETRIC_TEST :
      {
	Signature = TmpState;
	this->GetStateSingletParity(Signature);
	double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
	if ((TmpSign * this->SzParitySign) > 0.0)
	  {
	    //TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->LzParitySign) - 1.0))
	//							     * state[i] * M_SQRT1_2);
	    coefficient = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->LzParitySign) - 1.0))
								    * coefficient * M_SQRT1_2);
	    return this->FindStateIndex(TmpState, NewLzMax);
	  }
      }
      break;
    default:
      {
	Signature = TmpState;
	this->GetStateSingletParity(Signature);
	double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
// 	TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->LzParitySign) - 1.0))
// 								* (1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
// 								* state[i] * 0.5);
	coefficient = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->LzParitySign) - 1.0))
								* (1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
								* coefficient * 0.5);
	return this->FindStateIndex(TmpState, NewLzMax);
      }
      break;
    }    
  return this->HilbertSpaceDimension;  
}


