////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//        class of fermions on sphere with spin and two Landau levels         //
//                                                                            //
//                        last modification : 09/12/2016                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <math.h>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
//

FermionOnSphereWithSpinTwoLandauLevels::FermionOnSphereWithSpinTwoLandauLevels()
{
  this->LzMaxLLL = 0;
  this->LzMax2LL = 0;
  this->LzShift2LL = 0;
  this->LzShift2LL = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion in the lowest Landau level
// totalSpin = twice the total spin value
// totalIsospin = twice the total isospin value
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinTwoLandauLevels::FermionOnSphereWithSpinTwoLandauLevels (int nbrFermions, int totalLz, int lzMax, int totalSpin, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->TotalIsospin = 0;
  this->TotalEntanglement = 0;
  this->LzMax = lzMax + 2;
  this->NbrLzValue = this->LzMax + 1;
  this->LzMaxLLL = this->LzMax - 2;
  this->LzMax2LL = this->LzMax;
  this->LzShift2LL = 0;
  this->LzShiftLLL = 1;
  this->NbrFermionsUp = (this->TotalSpin + this->NbrFermions) >> 1;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrFermionsUpPlus = this->NbrFermionsUp;
  this->NbrFermionsDownPlus = this->NbrFermionsDown;
  this->NbrFermionsUpMinus = this->NbrFermionsUp;
  this->NbrFermionsDownMinus = this->NbrFermionsDown;

  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
										(this->TotalSpin + this->NbrFermions) >> 1);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;  
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						       (this->TotalSpin + this->NbrFermions) >> 1, 0l);
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << "Mismatch in State-count (" << this->LargeHilbertSpaceDimension << ") and State Generation (" 
	   << TmpHilbertSpaceDimension << ") in FermionOnSphereWithSpinTwoLandauLevels!" << endl;
      exit(1);
    }
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      long UsedMemory = 0;
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
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024l)
	if (UsedMemory >= 1048576l)
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

FermionOnSphereWithSpinTwoLandauLevels::FermionOnSphereWithSpinTwoLandauLevels(const FermionOnSphereWithSpinTwoLandauLevels& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->TotalIsospin = fermions.TotalIsospin;
  this->TotalEntanglement = fermions.TotalEntanglement;
  this->LzMaxLLL = fermions.LzMaxLLL;
  this->LzMax2LL = fermions.LzMax2LL;
  this->LzShift2LL = fermions.LzShift2LL;
  this->LzShiftLLL = fermions.LzShiftLLL;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermionsUpPlus = fermions.NbrFermionsUpPlus;
  this->NbrFermionsDownPlus = fermions.NbrFermionsDownPlus;
  this->NbrFermionsUpMinus = fermions.NbrFermionsUpMinus;
  this->NbrFermionsDownMinus = fermions.NbrFermionsDownMinus;
  this->HighestBit = fermions.HighestBit;
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

FermionOnSphereWithSpinTwoLandauLevels::~FermionOnSphereWithSpinTwoLandauLevels ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinTwoLandauLevels& FermionOnSphereWithSpinTwoLandauLevels::operator = (const FermionOnSphereWithSpinTwoLandauLevels& fermions)
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
  this->TotalIsospin = fermions.TotalIsospin;
  this->TotalEntanglement = fermions.TotalEntanglement;
  this->LzMaxLLL = fermions.LzMaxLLL;
  this->LzMax2LL = fermions.LzMax2LL;
  this->LzShift2LL = fermions.LzShift2LL;
  this->LzShiftLLL = fermions.LzShiftLLL;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermionsUpPlus = fermions.NbrFermionsUpPlus;
  this->NbrFermionsDownPlus = fermions.NbrFermionsDownPlus;
  this->NbrFermionsUpMinus = fermions.NbrFermionsUpMinus;
  this->NbrFermionsDownMinus = fermions.NbrFermionsDownMinus;
  this->HighestBit = fermions.HighestBit;
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

AbstractHilbertSpace* FermionOnSphereWithSpinTwoLandauLevels::Clone()
{
  return new FermionOnSphereWithSpinTwoLandauLevels(*this);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinTwoLandauLevels::GenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSpin, long pos)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalSpin > nbrFermions))
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
	  if (((this->LzMax2LL + this->LzShift2LL) >= totalLz) && (totalLz >= this->LzShift2LL))
	    {
	      this->StateDescription[pos] = 0x1ul << ((4 * totalLz) + 1 + (2 * totalSpin));
	      ++pos;
	    }
	  if (((this->LzMaxLLL + this->LzShiftLLL) >= totalLz) && (totalLz >= this->LzShiftLLL))
	    {
	      this->StateDescription[pos] = 0x1ul << ((4 * totalLz) + (2 * totalSpin));
	      ++pos;
	    }
	}
      return pos;
    }

  long TmpPos;
  unsigned long TmpMask;

  if ((lzMax <= (this->LzMax2LL + this->LzShift2LL)) && (lzMax >= this->LzShift2LL))
    {
      if ((lzMax <= (this->LzMaxLLL + this->LzShiftLLL)) && (lzMax >= this->LzShiftLLL))
	{
	  TmpMask = 0xful << (lzMax << 2);
	  TmpPos = this->GenerateStates(nbrFermions - 4, lzMax - 1, totalLz - (4 * lzMax), totalSpin - 2, pos);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= TmpMask;
	  TmpMask = 0xeul << (lzMax << 2);
	  TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, pos);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= TmpMask;
	  TmpMask = 0xdul << (lzMax << 2);
	  TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2, pos);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= TmpMask;
	  TmpMask = 0xcul << (lzMax << 2);
	  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 2, pos);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= TmpMask;
 	  TmpMask = 0xbul << (lzMax << 2);
	  TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, pos);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= TmpMask;
	}
      TmpMask = 0xaul << (lzMax << 2);
      TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= TmpMask;
      if ((lzMax <= (this->LzMaxLLL + this->LzShiftLLL)) && (lzMax >= this->LzShiftLLL))
	{
	  TmpMask = 0x9ul << (lzMax << 2);
	  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, pos);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= TmpMask;
	}
      TmpMask = 0x8ul << (lzMax << 2);
      TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= TmpMask;
      if ((lzMax <= (this->LzMaxLLL + this->LzShiftLLL)) && (lzMax >= this->LzShiftLLL))
	{
	  TmpMask = 0x7ul << (lzMax << 2);
	  TmpPos = this->GenerateStates(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1, pos);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= TmpMask;
	  TmpMask = 0x6ul << (lzMax << 2);
	  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, pos);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= TmpMask;
	}
    }
  if ((lzMax <= (this->LzMaxLLL + this->LzShiftLLL)) && (lzMax >= this->LzShiftLLL))
    {
      TmpMask = 0x5ul << (lzMax << 2);
      TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= TmpMask;
      TmpMask = 0x4ul << (lzMax << 2);
      TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= TmpMask;
   }
  if ((lzMax <= (this->LzMax2LL + this->LzShift2LL)) && (lzMax >= this->LzShift2LL))
    {
      if ((lzMax <= (this->LzMaxLLL + this->LzShiftLLL)) && (lzMax >= this->LzShiftLLL))
	{
	  TmpMask = 0x3ul << (lzMax << 2);
	  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin, pos);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= TmpMask;
	}
      TmpMask = 0x2ul << (lzMax << 2);
      TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= TmpMask;
    }
  if ((lzMax <= (this->LzMaxLLL + this->LzShiftLLL)) && (lzMax >= this->LzShiftLLL))
    {
      TmpMask = 0x1ul << (lzMax << 2);
      TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin, pos);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= TmpMask;
    }

  return this->GenerateStates(nbrFermions, lzMax - 1, totalLz, totalSpin, pos);
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// return value = Hilbert space dimension

long FermionOnSphereWithSpinTwoLandauLevels::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalSpin > nbrFermions))
    return 0l;

  if ((nbrFermions == 0) && (totalLz == 0))
    {
      return 1l;
    }

  if (lzMax < 0)
    return 0l;
  
  if (nbrFermions == 1) 
    {
      long Tmp = 0l;
      if (lzMax >= totalLz)
	{
	  if (((this->LzMax2LL + this->LzShift2LL) >= totalLz) && (totalLz >= this->LzShift2LL))
	    ++Tmp;
	  if (((this->LzMaxLLL + this->LzShiftLLL) >= totalLz) && (totalLz >= this->LzShiftLLL))
	    ++Tmp;
	}
      return Tmp;
    }


  long Tmp = 0l;
  if ((lzMax <= (this->LzMax2LL + this->LzShift2LL)) && (lzMax >= this->LzShift2LL))
    {
      if ((lzMax <= (this->LzMaxLLL + this->LzShiftLLL)) && (lzMax >= this->LzShiftLLL))
	{
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 4, lzMax - 1, totalLz - (4 * lzMax), totalSpin - 2);
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2);
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1);
	}
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1);
      if ((lzMax <= (this->LzMaxLLL + this->LzShiftLLL)) && (lzMax >= this->LzShiftLLL))
	{
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 2);
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 2);
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1);
	}
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1);
      if ((lzMax <= (this->LzMaxLLL + this->LzShiftLLL)) && (lzMax >= this->LzShiftLLL))
	{
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, lzMax - 1, totalLz - (3 * lzMax), totalSpin - 1);
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1);
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin);
	}
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin);
    }
  if ((lzMax <= (this->LzMaxLLL + this->LzShiftLLL)) && (lzMax >= this->LzShiftLLL))
    {
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), totalSpin - 1);
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin - 1);
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, totalSpin);
    }
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, totalSpin);
  return Tmp;
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = twce the total spin value
// return value = Hilbert space dimension

int FermionOnSphereWithSpinTwoLandauLevels::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + (nbrFermions * lzMax)) >> 1, 
						    (totalSpin + nbrFermions) >> 1);
}


