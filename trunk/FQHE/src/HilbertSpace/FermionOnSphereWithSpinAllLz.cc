////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin and                //
//                  without any contraint on the momentum along z             //
//                                                                            //
//                        last modification : 19/09/2016                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinAllLz.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/StringTools.h"

#include <cmath>
#include <bitset>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;


// default constructor
//

FermionOnSphereWithSpinAllLz::FermionOnSphereWithSpinAllLz()
{
}

// basic constructor
// 
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twice the total spin value
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinAllLz::FermionOnSphereWithSpinAllLz (int nbrFermions, int lzMax, int totalSpin, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin)/2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;

  this->Flag.Initialize();
  this->TargetSpace = this;
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];  

  this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, 0, 0l);

  this->GenerateLookUpTable(memory);


#ifdef __DEBUG__
  long UsedMemory = 0;
  UsedMemory += this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
  UsedMemory = this->NbrLzValue * sizeof(int);
  if (this->NbrFermions > 0)
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

FermionOnSphereWithSpinAllLz::FermionOnSphereWithSpinAllLz(const FermionOnSphereWithSpinAllLz& fermions)
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
  this->HighestBit = fermions.HighestBit;
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

FermionOnSphereWithSpinAllLz::~FermionOnSphereWithSpinAllLz ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinAllLz& FermionOnSphereWithSpinAllLz::operator = (const FermionOnSphereWithSpinAllLz& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
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
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinAllLz::Clone()
{
  return new FermionOnSphereWithSpinAllLz(*this);
}


  
// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalSpin = twice the total spin
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinAllLz::GenerateStates(int nbrFermions, int lzMax, int totalSpin, long pos)
{
  if (nbrFermions < 0)
    {
      return 0l;
    }
  if (lzMax < 0)
    {
      if ((totalSpin == this->TotalSpin) && (nbrFermions == 0))
	{
	  this->StateDescription[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else
	{
	  return pos;
	}
    }
  if (nbrFermions == 0)
    {
      if (totalSpin == this->TotalSpin)
	{
	  this->StateDescription[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else
	{
	  return pos;
	}
    }
  long TmpPos;
  unsigned long Mask;
  TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalSpin,  pos);
  Mask = 0x3ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalSpin + 1,  pos);
  Mask = 0x2ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalSpin - 1,  pos);
  Mask = 0x1ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, lzMax - 1, totalSpin,  pos);
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// posMax = highest position for next particle to be placed
// totalSpin = twice the total spin
// return value = Hilbert space dimension

long FermionOnSphereWithSpinAllLz::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalSpin)
{
  if (nbrFermions < 0)
    {
      return 0l;
    }
  if (lzMax < 0)
    {
      if ((totalSpin == this->TotalSpin) && (nbrFermions == 0))
	return 1l;
      else
	return 0l;
    }
  if (nbrFermions == 0)
    {
      if (totalSpin == this->TotalSpin)
	return 1l;
      else
	return 0l;
      
    }
  return  (this->EvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalSpin)
           + this->EvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalSpin + 1)
	   + this->EvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalSpin - 1)
	   + this->EvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalSpin));
}


