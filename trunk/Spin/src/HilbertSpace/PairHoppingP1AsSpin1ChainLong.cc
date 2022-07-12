////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of pair hopping p=1 (a.k.a. PXP model)             //
//                     Hilbert space written as spin 1 chain                  //
//                             for more than 32 spins                         //
//                                                                            //
//                        last modification : 18/03/2019                      //
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


#include "HilbertSpace/PairHoppingP1AsSpin1ChainLong.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;



// default constructor
//

PairHoppingP1AsSpin1ChainLong::PairHoppingP1AsSpin1ChainLong () 
{
  this->PeriodicBoundaryConditions = false;
  this->UseForEntanglementMatrix = false;
}

// constructor for complete Hilbert space
//
// chainLength = number of spin / group of 2p+1 orbitals
// periodicBoundaryConditions = true if the system uses periodic boundary conditions
// memorySize = memory size in bytes allowed for look-up table
// useForEntanglementMatrix = true if the hilbert space has to be generated for the entanglement matrix calculations

PairHoppingP1AsSpin1ChainLong::PairHoppingP1AsSpin1ChainLong (int chainLength, bool periodicBoundaryConditions, int memorySize, bool useForEntanglementMatrix) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->UseForEntanglementMatrix = useForEntanglementMatrix;
  
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(0, 0);

  this->StateDescription = new ULONGLONG [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, 0, 0);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory +=  this->LargeHilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->ChainLength * sizeof(int);
      UsedMemory += this->ChainLength * this->LookUpTableMemorySize * sizeof(int);
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
// chain = reference on chain to copy

PairHoppingP1AsSpin1ChainLong::PairHoppingP1AsSpin1ChainLong (const PairHoppingP1AsSpin1ChainLong& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->LookUpTable = chain.LookUpTable;
      this->MaximumLookUpShift = chain.MaximumLookUpShift;
      this->LookUpTableMemorySize = chain.LookUpTableMemorySize;
      this->LookUpTableShift = chain.LookUpTableShift;

      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->PeriodicBoundaryConditions = chain.PeriodicBoundaryConditions;
      this->UseForEntanglementMatrix = chain.UseForEntanglementMatrix;
    }
  else
    {
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0l;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->PeriodicBoundaryConditions = false;
      this->UseForEntanglementMatrix = false;

      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;
    }
}

// destructor
//

PairHoppingP1AsSpin1ChainLong::~PairHoppingP1AsSpin1ChainLong () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

PairHoppingP1AsSpin1ChainLong& PairHoppingP1AsSpin1ChainLong::operator = (const PairHoppingP1AsSpin1ChainLong& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTableShift;
      int TmpMaxBitPosition = 2 * this->ChainLength;
      for (int i = 0; i < TmpMaxBitPosition; ++i)
	this->LookUpTable[i];
      delete[] this->LookUpTable;
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->StateDescription = chain.StateDescription;
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->Sz = chain.Sz;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->PeriodicBoundaryConditions = chain.PeriodicBoundaryConditions;
      this->UseForEntanglementMatrix = chain.UseForEntanglementMatrix;

      this->LookUpTable = chain.LookUpTable;
      this->MaximumLookUpShift = chain.MaximumLookUpShift;
      this->LookUpTableMemorySize = chain.LookUpTableMemorySize;
      this->LookUpTableShift = chain.LookUpTableShift;
    }
  else
    {
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0l;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->PeriodicBoundaryConditions = false;
      this->UseForEntanglementMatrix = false;
 
      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* PairHoppingP1AsSpin1ChainLong::Clone()
{
  return new PairHoppingP1AsSpin1ChainLong (*this);
}

// evaluate Hilbert space dimension
//
// previousSpin = value of the previous spin (0 for -1, 1 for 0 and 2 for +1)
// sitePosition = site on chain where spin has to be changed
// return value = Hilbert space dimension

long PairHoppingP1AsSpin1ChainLong::EvaluateHilbertSpaceDimension(int previousSpin, int sitePosition)
{
  if (sitePosition == this->ChainLength)
    {
      return 1l;	  
    }
  long TmpDimension = 0l;
  if (sitePosition != 0)
    {
      if (previousSpin == 2)
	{
	  TmpDimension += this->EvaluateHilbertSpaceDimension(0, sitePosition + 1);
	}
      else
	{
	  TmpDimension += this->EvaluateHilbertSpaceDimension(2, sitePosition + 1);
	  TmpDimension += this->EvaluateHilbertSpaceDimension(1, sitePosition + 1);
	}
    }
  else
    {
      TmpDimension += this->EvaluateHilbertSpaceDimension(2, 1);
      TmpDimension += this->EvaluateHilbertSpaceDimension(1, 1);
      TmpDimension += this->EvaluateHilbertSpaceDimension(0, 1);
    }
  return TmpDimension;
}

// generate all states with neither constraint from boundary conditions nor discrete symmtry constraint
//
// statePosition = position for the new states
// previousSpin = value of the previous spin (0 for -1, 1 for 0 and 2 for +1)
// sitePosition = site on chain where spin has to be changed
// return value = number of generated states

long PairHoppingP1AsSpin1ChainLong::RawGenerateStates(long statePosition, int previousSpin, int sitePosition) 
{
  if (sitePosition == this->ChainLength)
    {
      this->StateDescription[statePosition] = ((ULONGLONG) 0x0ul);
      return (statePosition + 1l);
    }
  if (sitePosition != 0)
    {
      if (previousSpin == 2)
	{
	  statePosition = this->RawGenerateStates(statePosition, 0, sitePosition + 1);
	}
      else
	{
	  ULONGLONG TmpMask = ((ULONGLONG) 0x3ul) << (sitePosition * 2);
	  long TmpPosition = this->RawGenerateStates(statePosition, 2, sitePosition + 1);
	  for (; statePosition < TmpPosition; ++statePosition)
	    this->StateDescription[statePosition] |= TmpMask;
	  TmpMask = ((ULONGLONG) 0x2ul) << (sitePosition * 2);
	  TmpPosition = this->RawGenerateStates(statePosition, 1, sitePosition + 1);
	  for (; statePosition < TmpPosition; ++statePosition)
	    this->StateDescription[statePosition] |= TmpMask;
	}
    }
  else
    {
      ULONGLONG TmpMask = ((ULONGLONG) 0x3ul) << (sitePosition * 2);
      long TmpPosition = this->RawGenerateStates(statePosition, 2, 1);
      for (; statePosition < TmpPosition; ++statePosition)
	this->StateDescription[statePosition] |= TmpMask;
      TmpMask = ((ULONGLONG) 0x2ul) << (sitePosition * 2);
      TmpPosition = this->RawGenerateStates(statePosition, 1, 1);
      for (; statePosition < TmpPosition; ++statePosition)
	this->StateDescription[statePosition] |= TmpMask;
      statePosition = this->RawGenerateStates(statePosition, 0, 1);
    }
  return statePosition;
}


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long PairHoppingP1AsSpin1ChainLong::GenerateStates()
{
  long TmpHilbertSpaceDimension = 0l;
  ULONGLONG TmpDiscardMask = ~((ULONGLONG) 0x0ul);
  if (this->PeriodicBoundaryConditions == true)
    {
      int TmpLastUnitCellShift = (this->ChainLength - 2) * 2;
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  ULONGLONG Tmp = (this->StateDescription[i] & ((ULONGLONG) 0x3ul)) | ((this->StateDescription[i] >> TmpLastUnitCellShift) & 0xcul);
	  if ((Tmp == 0xful) || (Tmp == 0xeul) || (Tmp == 0x8ul) || (Tmp == ((ULONGLONG) 0x0ul)))
	    {
	      this->StateDescription[i] = TmpDiscardMask;
	    }
	  else
	    {
	      ++TmpHilbertSpaceDimension;
	    }
	}
    }
  else
    {
       if (this->UseForEntanglementMatrix == false)
	{
	  int TmpLastUnitCellShift = (this->ChainLength - 1) * 2;
	  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      if (((this->StateDescription[i] & ((ULONGLONG) 0x3ul)) == ((ULONGLONG) 0x0ul)) ||
		  (((this->StateDescription[i] >> TmpLastUnitCellShift) & ((ULONGLONG) 0x3ul)) == ((ULONGLONG) 0x3ul)))
		{
		  this->StateDescription[i] = TmpDiscardMask;
		}
	      else
		{
		  ++TmpHilbertSpaceDimension;
		}
	    }
	}
       else
	 {
	   TmpHilbertSpaceDimension = this->LargeHilbertSpaceDimension;
	 }
    }
  if (TmpHilbertSpaceDimension > 0l)
    {
      ULONGLONG* TmpStateDescription = new ULONGLONG [TmpHilbertSpaceDimension];
      TmpHilbertSpaceDimension = 0l;
      for  (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if ( this->StateDescription[i] != TmpDiscardMask)
	    {
	      TmpStateDescription[TmpHilbertSpaceDimension] = this->StateDescription[i];
	      ++TmpHilbertSpaceDimension;
	    }
	}
      SortArrayDownOrdering<ULONGLONG>(TmpStateDescription, TmpHilbertSpaceDimension);
      delete[] this->StateDescription;
      this->StateDescription = TmpStateDescription;
    }
  return TmpHilbertSpaceDimension;
}

// apply the swap operator within the unit cell
//
// unitCellCoordinate = coordinate of the unit cell
// siteIndex = index of the leftmost site within the unit cell
// state = index of the state on which the operator has to be applied
// return value = index of resulting state 

int PairHoppingP1AsSpin1ChainLong::SwapOperator (int unitCellCoordinate, int siteIndex, int state)
{
  return this->HilbertSpaceDimension;
}

// apply the operator coupling neighboring unit cells
//
// leftUnitCellCoordinate = coordinate of the left unit cell
// rightUnitCellCoordinate = coordinate of the right unit cell
// state = index of the state on which the operator has to be applied
// return value = index of resulting state 

int PairHoppingP1AsSpin1ChainLong::PlusMinusOperator (int leftUnitCellCoordinate, int rightUnitCellCoordinate, int state)
{
  ULONGLONG TmpStateDescription = this->StateDescription[state];
  leftUnitCellCoordinate <<= 1;
  rightUnitCellCoordinate <<= 1;
  ULONGLONG Tmp = ((TmpStateDescription >> leftUnitCellCoordinate) & ((ULONGLONG) 0x3ul)) | (((TmpStateDescription >> rightUnitCellCoordinate) & ((ULONGLONG) 0x3ul)) << 2);
  if (Tmp == ((ULONGLONG) 0xaul))
    {
      TmpStateDescription &= ~(((ULONGLONG) 0x3ul) << leftUnitCellCoordinate);
      TmpStateDescription |= ((ULONGLONG) 0x3ul) << leftUnitCellCoordinate;
      TmpStateDescription &= ~(((ULONGLONG) 0x3ul) << rightUnitCellCoordinate);
      TmpStateDescription |= ((ULONGLONG) 0x0ul) << rightUnitCellCoordinate;
      return this->FindStateIndex(TmpStateDescription);
   }
  if (Tmp == ((ULONGLONG) 0x3ul))
    {
      TmpStateDescription &= ~(((ULONGLONG) 0x3ul) << leftUnitCellCoordinate);
      TmpStateDescription |= ((ULONGLONG) 0x2ul) << leftUnitCellCoordinate;
      TmpStateDescription &= ~(((ULONGLONG) 0x3ul) << rightUnitCellCoordinate);
      TmpStateDescription |= ((ULONGLONG) 0x2ul) << rightUnitCellCoordinate;            
      return this->FindStateIndex(TmpStateDescription);
    }
  return this->HilbertSpaceDimension;
}

