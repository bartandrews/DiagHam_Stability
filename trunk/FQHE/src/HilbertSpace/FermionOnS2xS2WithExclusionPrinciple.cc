////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                  class of fermions on the 4D manifold S2 x S2              //
//               with an exclusion principle along within each sphere         //
//                                                                            //
//                        last modification : 08/12/2016                      //
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
#include "HilbertSpace/FermionOnS2xS2WithExclusionPrinciple.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <math.h>
#include <cstdlib>
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

FermionOnS2xS2WithExclusionPrinciple::FermionOnS2xS2WithExclusionPrinciple ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrFluxQuanta1 = number of flux quanta for the first sphere
// nbrFluxQuanta2 = number of flux quanta for the second sphere
// totalLz1 = total angular momentum for the first sphere
// totalLz2 = total angular momentum for the second sphere
// memory = amount of memory granted for precalculations
  
FermionOnS2xS2WithExclusionPrinciple::FermionOnS2xS2WithExclusionPrinciple (int nbrFermions, int nbrFluxQuanta1, int nbrFluxQuanta2, int totalLz1, int totalLz2, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrFluxQuanta1 = nbrFluxQuanta1;
  this->NbrFluxQuanta2 = nbrFluxQuanta2;
  this->TotalLz1 = totalLz1;
  this->TotalLz2 = totalLz2;
  this->NbrSiteX =  this->NbrFluxQuanta1 + 1;
  this->NbrSiteY =  this->NbrFluxQuanta2 + 1;
  this->KxMomentum = (this->TotalLz1 + (this->NbrFermions * this->NbrFluxQuanta1)) >> 1;
  this->KyMomentum = (this->TotalLz2 + (this->NbrFermions * this->NbrFluxQuanta2)) >> 1;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateLzMax = new int [this->HilbertSpaceDimension];  
      this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l);
      if (this->NbrFermions == 0)
	{
	  this->StateLzMax = new int [this->LargeHilbertSpaceDimension];
	  this->StateDescription[0] = 0x0ul; 
	  this->StateLzMax[0] = 0;
	}
      else
	{
	  unsigned long TmpMask1 = (0x1ul << this->NbrSiteY) | 0x1ul;
	  unsigned long TmpMask2 = (0x1ul << (2 * this->NbrSiteY)) | 0x1ul;
	  long TmpLargeHilbertSpaceDimension = 0l;
	  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      int TmpIndex = 0;
	      while ((TmpIndex < this->NbrLzValue) && (((this->StateDescription[i] >> TmpIndex) & TmpMask1) != TmpMask1) && (((this->StateDescription[i] >> TmpIndex) & TmpMask2) != TmpMask2))
		{
		  ++TmpIndex;		  
		}
	      if (TmpIndex < this->NbrLzValue)
		{
		  this->StateDescription[i] = 0x0ul;
		}
	      else
		{
		  unsigned long TmpMask3 = (0x3ul << this->NbrSiteY) | 0x3ul;
		  unsigned long TmpPattern0110 = (0x1ul << this->NbrSiteY) | 0x2ul;
		  unsigned long TmpPattern1001 = (0x2ul << this->NbrSiteY) | 0x1ul;
		  for (int k = 1;  (k < this->NbrSiteX) && (this->StateDescription[i] != 0x0ul); ++k)
		    {
		      for (int j = 1;  j < (this->NbrSiteY) && (this->StateDescription[i] != 0x0ul); ++j)
			{
			  if (((this->StateDescription[i] & TmpMask3) == TmpPattern1001) || ((this->StateDescription[i] & TmpMask3) == TmpPattern0110))
			    {
			      this->StateDescription[i] = 0x0ul;
			    }
			  TmpMask3 <<= 1;
			  TmpPattern0110 <<= 1;
			  TmpPattern1001 <<= 1;
			}
		      TmpMask3 <<= 1;
		      TmpPattern0110 <<= 1;
		      TmpPattern1001 <<= 1;
		    }
		  if (this->StateDescription[i] != 0x0ul)
		    {
		      ++TmpLargeHilbertSpaceDimension;
		    }
		}
	    }	  
	  if (TmpLargeHilbertSpaceDimension > 0l)
	    {
	      unsigned long* TmpStateDescriptions = new unsigned long[TmpLargeHilbertSpaceDimension];
	      TmpLargeHilbertSpaceDimension = 0l;
	      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
		{
		  if (this->StateDescription[i] != 0x0ul)
		    {
		      TmpStateDescriptions[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
		      ++TmpLargeHilbertSpaceDimension;
		    }
		}
	      delete[] this->StateDescription;
	      this->StateDescription = TmpStateDescriptions;
	      this->LargeHilbertSpaceDimension = TmpLargeHilbertSpaceDimension;	      
	      if (this->LargeHilbertSpaceDimension >= (1l << 30))
		this->HilbertSpaceDimension = 0;
	      else
		this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
	      this->StateLzMax = new int [this->LargeHilbertSpaceDimension];
	      int TmpLzMax = this->LzMax;
	      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
		{
		  while ((this->StateDescription[i] >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  this->StateLzMax[i] = TmpLzMax;
		}
	      this->TemporaryCanonicalArray = new unsigned long [this->LargeHilbertSpaceDimension];
	      this->TemporaryCanonicalArray3x3Block = new unsigned long [this->LargeHilbertSpaceDimension];
	    }
	  else
	    {
	      this->LargeHilbertSpaceDimension = 0l;
	      this->HilbertSpaceDimension = 0;
	    }
	}
      this->GenerateLookUpTable(memory);      
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
}
 
// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnS2xS2WithExclusionPrinciple::FermionOnS2xS2WithExclusionPrinciple(const FermionOnS2xS2WithExclusionPrinciple& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrFluxQuanta1 = fermions.NbrFluxQuanta1;
  this->NbrFluxQuanta2 = fermions.NbrFluxQuanta2;
  this->TotalLz1 = fermions.TotalLz1;
  this->TotalLz2 = fermions.TotalLz2;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
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

FermionOnS2xS2WithExclusionPrinciple::~FermionOnS2xS2WithExclusionPrinciple ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnS2xS2WithExclusionPrinciple& FermionOnS2xS2WithExclusionPrinciple::operator = (const FermionOnS2xS2WithExclusionPrinciple& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
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
  this->NbrFluxQuanta1 = fermions.NbrFluxQuanta1;
  this->NbrFluxQuanta2 = fermions.NbrFluxQuanta2;
  this->TotalLz1 = fermions.TotalLz1;
  this->TotalLz2 = fermions.TotalLz2;
  this->LzMax = fermions.LzMax;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnS2xS2WithExclusionPrinciple::Clone()
{
  return new FermionOnS2xS2WithExclusionPrinciple(*this);
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnS2xS2WithExclusionPrinciple::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions == 0)
    {
      if ((currentTotalKx == this->KxMomentum) && (currentTotalKy == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 3, currentTotalKx + currentKx, currentTotalKy + currentKy);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnS2xS2WithExclusionPrinciple::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions == 0)
    {
      if ((currentTotalKx == this->KxMomentum) && (currentTotalKy == this->KyMomentum))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  long TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 3, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
  unsigned long Mask = 0x1ul << ((currentKx * this->NbrSiteY) + currentKy);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos);
}


// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals

bool FermionOnS2xS2WithExclusionPrinciple::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  this->TemporaryCanonicalArray[0] = this->StateDescription[index];
  int Tmp = FindCanonical(this->StateDescription[index], 1);
  if (Tmp <= 0)
    {
      return false;
    }
  if (this->TemporaryCanonicalArray[Tmp - 1] == this->StateDescription[index])
    {
      return true;
    }
  else
    {
      return false;
    }
//   if (this->StateDescription[index] != this->FindCanonical(this->StateDescription[index], 0, 0))
//     {
//       return false;
//     }
  return true;
}

// test if a configuration satisfies the core exclusion principle 
//
// state = configuration to test
// return value = true if the configuration satisfies the core exclusion principle

bool FermionOnS2xS2WithExclusionPrinciple::CheckCoreExclusion(unsigned long state)
{
  unsigned long TmpMask = (0x3ul << this->NbrSiteY) | 0x3ul;
  unsigned long TmpPattern0110 = (0x1ul << this->NbrSiteY) | 0x2ul;
  unsigned long TmpPattern1001 = (0x2ul << this->NbrSiteY) | 0x1ul;
  int TmpNbrPatterns = 0;
  for (int i = 1;  i < this->NbrSiteX; ++i)
    {
      for (int j = 1;  j < this->NbrSiteY; ++j)
	{
	  if (((state & TmpMask) == TmpPattern0110) || ((state & TmpMask) == TmpPattern1001))
	    {
	      return false;
	    }
	  TmpMask <<= 1;
	  TmpPattern0110 <<= 1;
	  TmpPattern1001 <<= 1;
	}
      TmpMask <<= 1;
      TmpPattern0110 <<= 1;
      TmpPattern1001 <<= 1;
    }
  TmpMask = (0x1ul << (2 * this->NbrSiteY)) | (0x1ul << this->NbrSiteY) | 0x1ul;
  unsigned long TmpPattern101 = (0x1ul << (2 * this->NbrSiteY)) | 0x1ul;
  unsigned long TmpPattern110 = (0x1ul << this->NbrSiteY) | 0x1ul;
  for (int i = 2;  i < this->NbrSiteX; ++i)
    {
      for (int j = 0;  j < this->NbrSiteY; ++j)
	{
	  if (((state & TmpMask) == TmpPattern101) || ((state & TmpMask) == TmpPattern110))
	    {
	      return false;
	    }
	  TmpMask <<= 1;
	  TmpPattern101 <<= 1;
	  TmpPattern110 <<= 1;
	}
    }
  TmpMask = 0x7ul;
  TmpPattern101 = 0x5ul;
  TmpPattern110 = 0x3ul;
  for (int i = 0;  i < this->NbrSiteX; ++i)
    {
      for (int j = 2;  j < this->NbrSiteY; ++j)
	{
	  if (((state & TmpMask) == TmpPattern101) || ((state & TmpMask) == TmpPattern110))
	    {
	      return false;
	    }
	  TmpMask <<= 1;
	  TmpPattern101 <<= 1;
	  TmpPattern110 <<= 1;
	}
      TmpMask <<= 2;
      TmpPattern101 <<= 2;
      TmpPattern110 <<= 2;
    }
  return true;
}

// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals

unsigned long FermionOnS2xS2WithExclusionPrinciple::FindCanonical(unsigned long state, int xPosition, int yPosition)
{
  if (yPosition >= (this->NbrSiteY - 1))
    {
      yPosition = 0;
      xPosition++;
    }
  if (this->CheckCoreExclusion(state) == false)
    {
      return 0x0ul;
    }
  if (xPosition >= (this->NbrSiteX - 1))
    {
      return state;
    }
  int TmpShift = (xPosition *  this->NbrSiteY) + yPosition;
  unsigned long TmpMaskVertical = (0x3ul << this->NbrSiteY) | (0x3ul << (2 * this->NbrSiteY)) | 0x3ul;
  unsigned long TmpPattern01_00_10 = (0x1ul << (2 * this->NbrSiteY)) | 0x2ul;
  unsigned long TmpPattern10_00_01 = (0x2ul << (2 * this->NbrSiteY)) | 0x1ul;
  TmpMaskVertical <<= TmpShift;
  TmpPattern01_00_10 <<= TmpShift;
  TmpPattern10_00_01 <<= TmpShift;
  if ((state & TmpMaskVertical) == TmpPattern01_00_10)
    {
      unsigned long TmpState1 = this->FindCanonical(state, xPosition, yPosition + 1);
      unsigned long TmpState2 = this->FindCanonical((state & ~TmpMaskVertical) | TmpPattern10_00_01, xPosition, yPosition + 1);
      if ((TmpState1 == 0x0ul) || (TmpState2 == 0x0ul))
	return 0x0ul;
      if (TmpState1 < TmpState2)
	return TmpState2;
      else
	return TmpState1;
    }
  if ((state & TmpMaskVertical) == TmpPattern10_00_01)
    {
      unsigned long TmpState1 = this->FindCanonical(state, xPosition, yPosition + 1);
      unsigned long TmpState2 = this->FindCanonical((state & ~TmpMaskVertical) | TmpPattern01_00_10, xPosition, yPosition + 1);
      if ((TmpState1 == 0x0ul) || (TmpState2 == 0x0ul))
	return 0x0ul;
      if (TmpState1 < TmpState2)
	return TmpState2;
      else
	return TmpState1;      
    }
  if (yPosition < (this->NbrSiteY - 2))
    {
      unsigned long TmpMaskHorizontal = (0x7ul << this->NbrSiteY) | 0x7ul;
      unsigned long TmpPattern100_001 = (0x4ul << this->NbrSiteY) | 0x1ul;
      unsigned long TmpPattern001_100 = (0x1ul << this->NbrSiteY) | 0x4ul;
      TmpMaskHorizontal <<= TmpShift;
      TmpPattern100_001 <<= TmpShift;
      TmpPattern001_100 <<= TmpShift;
      if ((state & TmpMaskHorizontal) == TmpPattern100_001)
	{
	  unsigned long TmpState1 = this->FindCanonical(state, xPosition, yPosition + 1);
	  unsigned long TmpState2 = this->FindCanonical((state & ~TmpMaskHorizontal) | TmpPattern001_100, xPosition, yPosition + 1);
	  if ((TmpState1 == 0x0ul) || (TmpState2 == 0x0ul))
	    return 0x0ul;
	  if (TmpState1 < TmpState2)
	    return TmpState2;
	  else
	    return TmpState1;
	}
      if ((state & TmpMaskHorizontal) == TmpPattern001_100)
	{
	  unsigned long TmpState1 = this->FindCanonical(state, xPosition, yPosition + 1);
	  unsigned long TmpState2 = this->FindCanonical((state & ~TmpMaskHorizontal) | TmpPattern100_001, xPosition, yPosition + 1);
	  if ((TmpState1 == 0x0ul) || (TmpState2 == 0x0ul))
	    return 0x0ul;
	  if (TmpState1 < TmpState2)
	    return TmpState2;
	  else
	    return TmpState1;      
	}      
    }
  return this->FindCanonical(state, xPosition, yPosition + 1);  
}

// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals

int FermionOnS2xS2WithExclusionPrinciple::FindCanonical(unsigned long state, int currentPosition)
{
  if (this->CheckCoreExclusion(state) == false)
    {
      return -1;
    }  
  this->TemporaryCanonicalArray3x3Block[0] = state;
  if (this->FindCanonical3x3Block(state, 1) < 0)
    {
      return -1;
    }
//   if (this->TemporaryCanonicalArray3x3Block[0] != state)
//     {
//       return -1;
//     }

  unsigned long TmpMaskVertical = (0x3ul << this->NbrSiteY) | (0x3ul << (2 * this->NbrSiteY)) | 0x3ul;
  unsigned long TmpPattern01_00_10 = (0x1ul << (2 * this->NbrSiteY)) | 0x2ul;
  unsigned long TmpPattern10_00_01 = (0x2ul << (2 * this->NbrSiteY)) | 0x1ul;
  for (int i = 2;  i < this->NbrSiteX; ++i)
    {
      for (int j = 1;  j < this->NbrSiteY; ++j)
	{
	  if ((state & TmpMaskVertical) == TmpPattern01_00_10)
	    {
	      unsigned long TmpState = (state & ~TmpMaskVertical) | TmpPattern10_00_01;
	      int Tmp = SearchInSortedArrayAndInsert<unsigned long>(TmpState, this->TemporaryCanonicalArray, currentPosition);
	      if (Tmp == 1)
		{
		  currentPosition = this->FindCanonical(TmpState, currentPosition + 1);
		  if (currentPosition < 0)
		    return -1;		  
		}
	    }
	  else
	    {
	      if ((state & TmpMaskVertical) == TmpPattern10_00_01)
		{
		  unsigned long TmpState = (state & ~TmpMaskVertical) | TmpPattern01_00_10;
		  int Tmp = SearchInSortedArrayAndInsert<unsigned long>(TmpState, this->TemporaryCanonicalArray, currentPosition);
		  if (Tmp == 1)
		    {
		      currentPosition = this->FindCanonical(TmpState, currentPosition + 1);
		      if (currentPosition < 0)
			return -1;		      
		    }
		}
	    }
	  TmpMaskVertical <<= 1;
	  TmpPattern10_00_01 <<= 1;
	  TmpPattern01_00_10 <<= 1;
	}
      TmpMaskVertical <<= 1;
      TmpPattern10_00_01 <<= 1;
      TmpPattern01_00_10 <<= 1;
    }
  unsigned long TmpMaskHorizontal = (0x7ul << this->NbrSiteY) | 0x7ul;
  unsigned long TmpPattern100_001 = (0x4ul << this->NbrSiteY) | 0x1ul;
  unsigned long TmpPattern001_100 = (0x1ul << this->NbrSiteY) | 0x4ul;
  for (int i = 1;  i < this->NbrSiteX; ++i)
    {
      for (int j = 2;  j < this->NbrSiteY; ++j)
	{
	  if ((state & TmpMaskHorizontal) == TmpPattern100_001)
	    {
	      unsigned long TmpState = (state & ~TmpMaskHorizontal) | TmpPattern001_100;
	      int Tmp = SearchInSortedArrayAndInsert<unsigned long>(TmpState, this->TemporaryCanonicalArray, currentPosition);
	      if (Tmp == 1)
		{
		  currentPosition = this->FindCanonical(TmpState, currentPosition + 1);
		  if (currentPosition < 0)
		    return -1;		      
		}
	    }
	  if ((state & TmpMaskHorizontal) == TmpPattern001_100)
	    {
	      unsigned long TmpState = (state & ~TmpMaskHorizontal) | TmpPattern100_001;
	      int Tmp = SearchInSortedArrayAndInsert<unsigned long>(TmpState, this->TemporaryCanonicalArray, currentPosition);
	      if (Tmp == 1)
		{
		  currentPosition = this->FindCanonical(TmpState, currentPosition + 1);
		  if (currentPosition < 0)
		    return -1;		      
		}
	    }      
	  TmpMaskHorizontal <<= 1;
	  TmpPattern100_001 <<= 1;
	  TmpPattern001_100 <<= 1;
	}
      TmpMaskHorizontal <<= 2;
      TmpPattern100_001 <<= 2;
      TmpPattern001_100 <<= 2;
    }
//   unsigned long TmpMaskBlock = (0x7ul << (2 * this->NbrSiteY)) | (0x7ul << this->NbrSiteY) | 0x7ul;
//   unsigned long TmpPattern100_000_001 = (0x4ul << (2 *this->NbrSiteY)) | 0x1ul;
//   unsigned long TmpPattern001_000_100 = (0x1ul << (2 *this->NbrSiteY)) | 0x4ul;
//   for (int i = 2;  i < this->NbrSiteX; ++i)
//     {
//       for (int j = 2;  j < this->NbrSiteY; ++j)
// 	{
// 	  if ((state & TmpMaskBlock) == TmpPattern100_000_001)
// 	    {
// 	      unsigned long TmpState = (state & ~TmpMaskBlock) | TmpPattern001_000_100;
// 	      int Tmp = SearchInSortedArrayAndInsert<unsigned long>(TmpState, this->TemporaryCanonicalArray, currentPosition);
// 	      if (Tmp == 1)
// 		{
// 		  currentPosition = this->FindCanonical(TmpState, currentPosition + 1);
// 		  if (currentPosition < 0)
// 		    return -1;		      
// 		}
// 	    }
// 	  if ((state & TmpMaskBlock) == TmpPattern001_000_100)
// 	    {
// 	      unsigned long TmpState = (state & ~TmpMaskBlock) | TmpPattern100_000_001;
// 	      int Tmp = SearchInSortedArrayAndInsert<unsigned long>(TmpState, this->TemporaryCanonicalArray, currentPosition);
// 	      if (Tmp == 1)
// 		{
// 		  currentPosition = this->FindCanonical(TmpState, currentPosition + 1);
// 		  if (currentPosition < 0)
// 		    return -1;		      
// 		}
// 	    }      
// 	  TmpMaskBlock <<= 1;
// 	  TmpPattern100_000_001 <<= 1;
// 	  TmpPattern001_000_100 <<= 1;
// 	}
//       TmpMaskBlock <<= 2;
//       TmpPattern100_000_001 <<= 2;
//       TmpPattern001_000_100 <<= 2;
//     }
  return currentPosition;
}

// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals

int FermionOnS2xS2WithExclusionPrinciple::FindCanonical3x3Block(unsigned long state, int currentPosition)
{
  unsigned long TmpMaskBlock = (0x7ul << (2 * this->NbrSiteY)) | (0x7ul << this->NbrSiteY) | 0x7ul;
  unsigned long TmpPattern100_000_001 = (0x4ul << (2 *this->NbrSiteY)) | 0x1ul;
  unsigned long TmpPattern001_000_100 = (0x1ul << (2 *this->NbrSiteY)) | 0x4ul;

  unsigned long TmpPattern101_000_001 = (0x5ul << (2 *this->NbrSiteY)) | 0x1ul;
  unsigned long TmpPattern100_000_101 = (0x4ul << (2 *this->NbrSiteY)) | 0x5ul;
  unsigned long TmpPattern001_000_101 = (0x1ul << (2 *this->NbrSiteY)) | 0x5ul;
  unsigned long TmpPattern101_000_100 = (0x5ul << (2 *this->NbrSiteY)) | 0x4ul;
  for (int i = 2;  i < this->NbrSiteX; ++i)
    {
      for (int j = 2;  j < this->NbrSiteY; ++j)
	{
	  if (((state & TmpMaskBlock) == TmpPattern101_000_001) || ((state & TmpMaskBlock) == TmpPattern101_000_100) ||
	      ((state & TmpMaskBlock) == TmpPattern100_000_101) || ((state & TmpMaskBlock) == TmpPattern001_000_101))
	    {
	      return -1;
	    }
	  if ((state & TmpMaskBlock) == TmpPattern100_000_001)
	    {
	      unsigned long TmpState = (state & ~TmpMaskBlock) | TmpPattern001_000_100;
	      int Tmp = SearchInSortedArrayAndInsert<unsigned long>(TmpState, this->TemporaryCanonicalArray3x3Block, currentPosition);
	      if (Tmp == 1)
		{
		  currentPosition = this->FindCanonical3x3Block(TmpState, currentPosition + 1);
		  if (currentPosition < 0)
		    return -1;		      
		}
	    }
	  if ((state & TmpMaskBlock) == TmpPattern001_000_100)
	    {
	      unsigned long TmpState = (state & ~TmpMaskBlock) | TmpPattern100_000_001;
	      int Tmp = SearchInSortedArrayAndInsert<unsigned long>(TmpState, this->TemporaryCanonicalArray3x3Block, currentPosition);
	      if (Tmp == 1)
		{
		  currentPosition = this->FindCanonical3x3Block(TmpState, currentPosition + 1);
		  if (currentPosition < 0)
		    return -1;		      
		}
	    }      
	  TmpMaskBlock <<= 1;
	  TmpPattern100_000_001 <<= 1;
	  TmpPattern001_000_100 <<= 1; 
	  TmpPattern101_000_001 <<= 1;
	  TmpPattern100_000_101 <<= 1;
	  TmpPattern001_000_101 <<= 1;
	  TmpPattern101_000_100 <<= 1;
	}
      TmpMaskBlock <<= 2;
      TmpPattern100_000_001 <<= 2;
      TmpPattern001_000_100 <<= 2;
      TmpPattern101_000_001 <<= 2;
      TmpPattern100_000_101 <<= 2;
      TmpPattern001_000_101 <<= 2;
      TmpPattern101_000_100 <<= 2;
    }
  return currentPosition;
}
