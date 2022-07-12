////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere including two                //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 19/05/2009                      //
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
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"
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

FermionOnSphereTwoLandauLevels::FermionOnSphereTwoLandauLevels()
{
}

// basic constructor with no contraint on the number of particles per spin component
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMaxUp = twice the maximum Lz value reached by a fermion with a spin up
// lzMaxDown = twice the maximum Lz value reached by a fermion with a spin down
// memory = amount of memory granted for precalculations

FermionOnSphereTwoLandauLevels::FermionOnSphereTwoLandauLevels (int nbrFermions, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->LzMaxUp = lzMaxUp;
  this->LzMaxDown = lzMaxDown;
  if (this->LzMaxUp >= this->LzMaxDown)
    {
      this->LzMax = this->LzMaxUp;
      this->LzShiftUp = 0;
      this->LzShiftDown = (this->LzMaxUp - this->LzMaxDown) >> 1;
    }
  else
    {
      this->LzMax = this->LzMaxDown;
      this->LzShiftDown = 0;
      this->LzShiftUp = (this->LzMaxDown - this->LzMaxUp) >> 1;
    }
  this->LzTotalShift = this->LzMaxDown + this->LzMaxUp;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateFullHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  int TmpDimension = this->GenerateFullStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 0);
  if (TmpDimension != this->HilbertSpaceDimension)
    {
      cout << "Mismatch in State-count and State Generation in FermionOnSphereTwoLandauLevels! " << this->HilbertSpaceDimension << " " << TmpDimension  << endl;
  for (int i = 0; i < TmpDimension; ++i)
    this->PrintState(cout, i) << endl;
       exit(1);
    }

  //  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermionsUp, this->NbrFermionsDown, this->LzMaxUp, this->LzMaxDown, );
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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereTwoLandauLevels::FermionOnSphereTwoLandauLevels(const FermionOnSphereTwoLandauLevels& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->LzMaxUp = fermions.LzMaxUp;
  this->LzMaxDown = fermions.LzMaxDown;
  this->LzShiftUp = fermions.LzShiftUp;
  this->LzShiftDown = fermions.LzShiftDown;
  this->LzTotalShift = fermions.LzTotalShift;
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

FermionOnSphereTwoLandauLevels::~FermionOnSphereTwoLandauLevels ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereTwoLandauLevels& FermionOnSphereTwoLandauLevels::operator = (const FermionOnSphereTwoLandauLevels& fermions)
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
  this->LzMaxUp = fermions.LzMaxUp;
  this->LzMaxDown = fermions.LzMaxDown;
  this->LzShiftUp = fermions.LzShiftUp;
  this->LzShiftDown = fermions.LzShiftDown;
  this->LzTotalShift = fermions.LzTotalShift;
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

AbstractHilbertSpace* FermionOnSphereTwoLandauLevels::Clone()
{
  return new FermionOnSphereTwoLandauLevels(*this);
}

// generate all states corresponding to the constraints
// 
// nbrFermionsUp = number of fermions with spin up
// nbrFermionsDown = number of fermions with spin down
// lzMaxUp = momentum maximum value for a fermion
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereTwoLandauLevels::GenerateStates(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz, long pos)
{
  cout << "warning : untested code" << endl;
  if ((nbrFermionsUp < 0) || (nbrFermionsDown < 0) || (totalLz < 0))
    return pos;
  if ((nbrFermionsUp == 0) && (totalLz == 0) && (nbrFermionsDown == 0))
      {
	this->StateDescription[pos] = 0x0ul;
	return (pos + 1l);
      }
    
  if (lzMax < 0)
    return pos;
  
  if ((nbrFermionsUp == 1) && (nbrFermionsDown == 0)) 
    if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
      {
	this->StateDescription[pos] = 0x2ul << (totalLz << 1);
	return (pos + 1l);
      }
    else
      return pos;
  if ((nbrFermionsDown == 1) && (nbrFermionsUp == 0)) 
    if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
      {
	this->StateDescription[pos] = 0x1ul << (totalLz << 1);
	return (pos + 1l);
      }
    else
      return pos;

  if ((lzMax == 0) && (totalLz != 0))
    return pos;


  if (((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp)) &&
      ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown)))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp - 1, nbrFermionsDown - 1, lzMax - 1, totalLz - (2 * lzMax), pos);
      unsigned long Mask = 0x3ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp - 1, nbrFermionsDown, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x2ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp, nbrFermionsDown  - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x1ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  return this->GenerateStates(nbrFermionsUp, nbrFermionsDown, lzMax - 1, totalLz, pos);
};

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereTwoLandauLevels::GenerateFullStates(int nbrFermions, int lzMax, int totalLz, long pos)
{
  if (nbrFermions < 0)
    return pos;
  if ((nbrFermions == 0) && (totalLz == 0))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
    
  if (lzMax < 0)
    return pos;
  
  if (nbrFermions == 1) 
     {
       if (lzMax >= totalLz)
	 {
	   if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
	     {
	       this->StateDescription[pos] = 0x2ul << (totalLz << 1);
	       ++pos;
	     }
	   if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
	     {
	       this->StateDescription[pos] = 0x1ul << (totalLz << 1);
	       ++pos;
	     }
	 }
       return pos;
     }

  if ((lzMax == 0) && (totalLz != 0))
    return pos;

  if (((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp)) &&
      ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown)))
    {
      long TmpPos = this->GenerateFullStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), pos);
      unsigned long Mask = 0x3ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      long TmpPos = this->GenerateFullStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x2ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
    {
      long TmpPos = this->GenerateFullStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x1ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }

  return this->GenerateFullStates(nbrFermions, lzMax - 1, totalLz, pos);
}

// evaluate Hilbert space dimension
//
// nbrFermionsUp = number of fermions with spin up
// nbrFermionsDown = number of fermions with spin down
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereTwoLandauLevels::ShiftedEvaluateHilbertSpaceDimension(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz)
{
  cout << "warning : untested code" << endl;
  if ((nbrFermionsUp < 0) || (nbrFermionsDown < 0) || (totalLz < 0))
    return 0l;
  if ((nbrFermionsUp == 0) && (nbrFermionsDown == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if ((nbrFermionsUp == 1) && (nbrFermionsDown == 0)) 
    if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
      return 1l;
    else
      return 0l;
  if ((nbrFermionsUp == 0) && (nbrFermionsDown == 1)) 
    if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
      return 1l;
    else
      return 0l;

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;

  long Tmp = 0l;
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp - 1, nbrFermionsDown - 1, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp - 1, nbrFermionsDown, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp, nbrFermionsDown - 1, lzMax - 1, totalLz - lzMax);
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp, nbrFermionsDown, lzMax- 1, totalLz);
  return Tmp;
}

// evaluate Hilbert space dimension without constraint on the number of particles per level
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereTwoLandauLevels::ShiftedEvaluateFullHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions < 0) || (totalLz < 0))
    return 0l;
  if ((nbrFermions == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if (nbrFermions == 1) 
    {
      long Tmp = 0l;
      if (lzMax >= totalLz)
	{
	  if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
	    ++Tmp;
	  if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
	    ++Tmp;
	}
      return Tmp;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;

  long Tmp = 0l;
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
	Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions -1, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))    
    Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax);
  Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz);
  return Tmp;
}


// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

RealVector FermionOnSphereTwoLandauLevels::ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace)
{
  RealVector FinalState(this->HilbertSpaceDimension, true);
  for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
    {
      unsigned long TmpUpState = upStateSpace.StateDescription[j] << this->LzShiftUp;
      int TmpPos = upStateSpace.LzMax + this->LzShiftUp;
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
      unsigned long TmpDownState = downStateSpace.StateDescription[j] << this->LzShiftDown;
      int TmpPos = downStateSpace.LzMax + this->LzShiftDown;
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

