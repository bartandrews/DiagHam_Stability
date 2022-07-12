////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of pair hopping with generic p                  //
//                     Hilbert space written as spin 1 chain                  //
//                                                                            //
//                        last modification : 03/11/2019                      //
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


#include "HilbertSpace/PairHoppingGenericPAsSpin1Chain.h"
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

PairHoppingGenericPAsSpin1Chain::PairHoppingGenericPAsSpin1Chain () 
{
  this->PeriodicBoundaryConditions = false;
  this->NbrMinusLeft = -1;
  this->NbrPlusRight = -1;
  this->UnitCellSize = 3;
  this->FirstUnitCellMask = (0x1ul << (this->UnitCellSize << 1)) - 0x1ul;
  this->PossibleUnitCellsPerNbrMinus = 0;
  this->NbrPossibleUnitCellsPerNbrMinus = 0;
  this->PossibleUnitCellsNbrPlusPerNbrMinus = 0;
  this->PossibleUnitCells = 0;
  this->PossibleUnitCellsNbrPlus = 0;
}

// constructor for complete Hilbert space
//
// chainLength = number of spin / group of 2p+1 orbitals
// pValue = p value
// periodicBoundaryConditions = true if the system uses periodic boundary conditions
// memorySize = memory size in bytes allowed for look-up table
// useForEntanglementMatrix = true if the hilbert space has to be generated for the entanglement matrix calculations

PairHoppingGenericPAsSpin1Chain::PairHoppingGenericPAsSpin1Chain (int chainLength, int pValue, bool periodicBoundaryConditions, int memorySize, bool useForEntanglementMatrix) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->PeriodicBoundaryConditions = periodicBoundaryConditions;
  this->UseForEntanglementMatrix = useForEntanglementMatrix;
  this->NbrMinusLeft = -1;
  this->NbrPlusRight = -1;
  this->UnitCellSize = pValue;
  this->FirstUnitCellMask = (0x1ul << (this->UnitCellSize << 1)) - 0x1ul;

  this->GenerateAllPossibleUnitCells();
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(0, 0);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, 0, 0);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory +=  this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
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

// constructor for the Hilbert space, fixing the number of pluses in the rightmost unit cell and the number of minuses in the leftmost unit cell
//
// chainLength = number of spin / group of 2p+1 orbitals
// nbrPlusRight = number of pluses in the rightmost unit cell
// nbrMinusLeft = number of minuses in the lefttmost unit cell

PairHoppingGenericPAsSpin1Chain::PairHoppingGenericPAsSpin1Chain (int chainLength, int nbrPlusRight, int nbrMinusLeft, int memorySize)
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->PeriodicBoundaryConditions = false;
  this->UseForEntanglementMatrix = true;
  this->NbrMinusLeft = nbrMinusLeft;
  this->NbrPlusRight = nbrPlusRight;
  this->UnitCellSize = 3;
  this->FirstUnitCellMask = (0x1ul << (this->UnitCellSize << 1)) - 0x1ul;
  this->GenerateAllPossibleUnitCells();
  
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(0, 0);
  
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, 0, 0);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory +=  this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
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

PairHoppingGenericPAsSpin1Chain::PairHoppingGenericPAsSpin1Chain (const PairHoppingGenericPAsSpin1Chain& chain) 
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

      this->UnitCellSize =  chain.UnitCellSize;
      this->FirstUnitCellMask = chain.FirstUnitCellMask;
      this->GenerateAllPossibleUnitCells();
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

      this->UnitCellSize = 0;
      this->FirstUnitCellMask = 0x0ul;
      this->PossibleUnitCellsPerNbrMinus = 0;
      this->NbrPossibleUnitCellsPerNbrMinus = 0;      
      this->PossibleUnitCellsNbrPlusPerNbrMinus = 0;
      this->PossibleUnitCells = 0;
      this->PossibleUnitCellsNbrPlus = 0;
    }
}

// destructor
//

PairHoppingGenericPAsSpin1Chain::~PairHoppingGenericPAsSpin1Chain () 
{
  if (this->NbrPossibleUnitCellsPerNbrMinus != 0)
    {
      for (int i = 0; i <= this->UnitCellSize; ++i)
	{
	  if (this->NbrPossibleUnitCellsPerNbrMinus[i] != 0)
	    {
	      delete[] this->PossibleUnitCellsPerNbrMinus[i];
	      delete[] this->PossibleUnitCellsNbrPlusPerNbrMinus[i];
	    }
	}
      delete[] this->PossibleUnitCellsPerNbrMinus;
      delete[] this->NbrPossibleUnitCellsPerNbrMinus;
      delete[] this->PossibleUnitCellsNbrPlusPerNbrMinus;
      delete[] this->PossibleUnitCells;
      delete[] this->PossibleUnitCellsNbrPlus;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

PairHoppingGenericPAsSpin1Chain& PairHoppingGenericPAsSpin1Chain::operator = (const PairHoppingGenericPAsSpin1Chain& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTableShift;
      int TmpMaxBitPosition = 2 * this->ChainLength;
      for (int i = 0; i < TmpMaxBitPosition; ++i)
	this->LookUpTable[i];
      delete[] this->LookUpTable;
      if (this->NbrPossibleUnitCellsPerNbrMinus != 0)
	{
	  for (int i = 0; i <= this->UnitCellSize; ++i)
	    {
	      if (this->NbrPossibleUnitCellsPerNbrMinus[i] != 0)
		{
		  delete[] this->PossibleUnitCellsPerNbrMinus[i];
		  delete[] this->PossibleUnitCellsNbrPlusPerNbrMinus[i];
		}
	    }
	  delete[] this->PossibleUnitCellsPerNbrMinus;
	  delete[] this->NbrPossibleUnitCellsPerNbrMinus;
	  delete[] this->PossibleUnitCellsNbrPlusPerNbrMinus;
	  delete[] this->PossibleUnitCells;
	  delete[] this->PossibleUnitCellsNbrPlus;
	}
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

      this->UnitCellSize =  chain.UnitCellSize;
      this->FirstUnitCellMask = chain.FirstUnitCellMask;
      this->GenerateAllPossibleUnitCells();
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

      this->UnitCellSize = 0;
      this->FirstUnitCellMask = 0x0ul;
      this->PossibleUnitCellsPerNbrMinus = 0;
      this->NbrPossibleUnitCellsPerNbrMinus = 0;      
      this->PossibleUnitCellsNbrPlusPerNbrMinus = 0;
      this->PossibleUnitCells = 0;
      this->PossibleUnitCellsNbrPlus = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* PairHoppingGenericPAsSpin1Chain::Clone()
{
  return new PairHoppingGenericPAsSpin1Chain (*this);
}

// evaluate Hilbert space dimension
//
// previousNbrPlus = number of +1 spins in the previous unit cell
// sitePosition = site on chain where spin has to be changed
// return value = Hilbert space dimension

long PairHoppingGenericPAsSpin1Chain::EvaluateHilbertSpaceDimension(int previousNbrPlus, int sitePosition)
{
  if (sitePosition == this->ChainLength)
    {
      return 1l;	  
    }
  long TmpDimension = 0l;
  if (sitePosition != 0)
    {
      for (int i = 0; i < this->NbrPossibleUnitCellsPerNbrMinus[previousNbrPlus]; ++i)
	{
	  TmpDimension += this->EvaluateHilbertSpaceDimension(this->PossibleUnitCellsNbrPlusPerNbrMinus[previousNbrPlus][i], sitePosition + this->UnitCellSize);	  
	}
    }
  else
    {
      for (int i = 0; i < this->NbrPossibleUnitCells; ++i)
	{
	  TmpDimension += this->EvaluateHilbertSpaceDimension(this->PossibleUnitCellsNbrPlus[i], this->UnitCellSize);	  
	}
    }
  return TmpDimension;
}

// generate all states with neither constraint from boundary conditions nor discrete symmtry constraint
//
// statePosition = position for the new states
// previousNbrPlus = number of +1 spins in the previous unit cell
// sitePosition = site on chain where spin has to be changed
// return value = number of generated states

long PairHoppingGenericPAsSpin1Chain::RawGenerateStates(long statePosition, int previousNbrPlus, int sitePosition) 
{
  if (sitePosition == this->ChainLength)
    {
      this->StateDescription[statePosition] = 0x0ul;
      return (statePosition + 1l);
    }
  if (sitePosition != 0)
    {
      for (int i = 0; i < this->NbrPossibleUnitCellsPerNbrMinus[previousNbrPlus]; ++i)
	{
	  unsigned long TmpMask = this->PossibleUnitCellsPerNbrMinus[previousNbrPlus][i] << (sitePosition * 2);
	  long TmpPosition = this->RawGenerateStates(statePosition, this->PossibleUnitCellsNbrPlusPerNbrMinus[previousNbrPlus][i], sitePosition + this->UnitCellSize);
	  for (; statePosition < TmpPosition; ++statePosition)
	    this->StateDescription[statePosition] |= TmpMask;	  
	}
    }
  else
    {
      for (int i = 0; i < this->NbrPossibleUnitCells; ++i)
	{
	  unsigned long TmpMask = this->PossibleUnitCells[i];
	  long TmpPosition = this->RawGenerateStates(statePosition, this->PossibleUnitCellsNbrPlus[i], this->UnitCellSize);
	  for (; statePosition < TmpPosition; ++statePosition)
	    this->StateDescription[statePosition] |= TmpMask;	  
	}
    }
  return statePosition;
}


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long PairHoppingGenericPAsSpin1Chain::GenerateStates()
{
  cout << "p=" << this->UnitCellSize << " has " << this->NbrPossibleUnitCells << " possible unit cells" << endl;
  long TmpHilbertSpaceDimension = 0l;
  unsigned long TmpDiscardMask = ~0x0ul;
  int TmpNbrBitsPerUnitCell = 1 << (this->UnitCellSize << 1);
  int* TmpNbrMinusTable = new int[TmpNbrBitsPerUnitCell];
  int* TmpNbrPlusTable = new int[TmpNbrBitsPerUnitCell];
  int TmpLastUnitCellShift = (this->ChainLength - this->UnitCellSize) * 2;
  for (int i = 0; i < TmpNbrBitsPerUnitCell; ++i)
    {
      TmpNbrMinusTable[i] = 0;
      TmpNbrPlusTable[i] = 0;
      for (int j = 0; j < this->UnitCellSize; ++j)
	{
	  int Tmp = (i >> (2 * j)) & 0x3;
	  if (Tmp == 0)
	    {
	      TmpNbrMinusTable[i]++;
	    }
	  else
	    {
	      if (Tmp == 3)
		{
		  TmpNbrPlusTable[i]++;
		}
	    }
	}      
    }
  if (this->PeriodicBoundaryConditions == true)
    {
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if (TmpNbrPlusTable[(this->StateDescription[i] >> TmpLastUnitCellShift) & this->FirstUnitCellMask] != TmpNbrMinusTable[this->StateDescription[i]  & this->FirstUnitCellMask])
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
	  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      if ((TmpNbrPlusTable[(this->StateDescription[i] >> TmpLastUnitCellShift) & this->FirstUnitCellMask] != 0) ||
		  (TmpNbrMinusTable[this->StateDescription[i]  & this->FirstUnitCellMask] != 0))
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
	  if ((this->NbrMinusLeft != - 1) && (this->NbrPlusRight != -1))
	    {
	      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
		{
		  if ((TmpNbrPlusTable[(this->StateDescription[i] >> TmpLastUnitCellShift) & this->FirstUnitCellMask] != this->NbrPlusRight) ||
		      (TmpNbrMinusTable[this->StateDescription[i]  & this->FirstUnitCellMask] != this->NbrMinusLeft))
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
    }
  if (TmpHilbertSpaceDimension > 0l)
    {
      unsigned long* TmpStateDescription = new unsigned long [TmpHilbertSpaceDimension];
      TmpHilbertSpaceDimension = 0l;
      for  (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if ( this->StateDescription[i] != TmpDiscardMask)
	    {
	      TmpStateDescription[TmpHilbertSpaceDimension] = this->StateDescription[i];
	      ++TmpHilbertSpaceDimension;
	    }
	}
      SortArrayDownOrdering<unsigned long>(TmpStateDescription, TmpHilbertSpaceDimension);
      delete[] this->StateDescription;
      this->StateDescription = TmpStateDescription;
    }
  delete[] TmpNbrMinusTable;
  delete[] TmpNbrPlusTable;

  return TmpHilbertSpaceDimension;
}

// apply the swap operator within the unit cell
//
// unitCellCoordinate = coordinate of the unit cell
// siteIndex = index of the leftmost site within the unit cell
// state = index of the state on which the operator has to be applied
// return value = index of resulting state 

int PairHoppingGenericPAsSpin1Chain::SwapOperator (int unitCellCoordinate, int siteIndex, int state)
{
  unsigned long TmpStateDescription = this->StateDescription[state];
  unitCellCoordinate *= this->UnitCellSize;
  unitCellCoordinate += siteIndex;
  unitCellCoordinate *= 2;
  unsigned long Tmp = (TmpStateDescription >> unitCellCoordinate) & this->FirstUnitCellMask;
  switch (Tmp)
    {
    case 0xbul:
      {
	TmpStateDescription &= ~(0xful << unitCellCoordinate);
	TmpStateDescription |= 0xeul << unitCellCoordinate;
	return this->FindStateIndex(TmpStateDescription);
      }
      break;
    case 0xeul:
      {
	TmpStateDescription &= ~(0xful << unitCellCoordinate);
	TmpStateDescription |= 0xbul << unitCellCoordinate;
	return this->FindStateIndex(TmpStateDescription);
      }
      break;
    case 0x2ul:
      {
	TmpStateDescription &= ~(0xful << unitCellCoordinate);
	TmpStateDescription |= 0x8ul << unitCellCoordinate;
	return this->FindStateIndex(TmpStateDescription);
      }
      break;
    case 0x8ul:
      {
	TmpStateDescription &= ~(0xful << unitCellCoordinate);
	TmpStateDescription |= 0x2ul << unitCellCoordinate;
	return this->FindStateIndex(TmpStateDescription);
      }
      break;
    }  
  return this->HilbertSpaceDimension;
}

// apply the swap operator within the unit cell with a constraint of the unit cell parity
//
// unitCellCoordinate = coordinate of the unit cell
// siteIndex = index of the leftmost site within the unit cell
// state = index of the state on which the operator has to be applied
// return value = index of resulting state 

int PairHoppingGenericPAsSpin1Chain::SwapOperatorPlus (int unitCellCoordinate, int siteIndex, int state)
{
  unsigned long TmpStateDescription = this->StateDescription[state];
  if ((unitCellCoordinate & 1) == 0)
    {
      unitCellCoordinate *= this->UnitCellSize;
      unitCellCoordinate += siteIndex;
      unitCellCoordinate *= 2;
      unsigned long Tmp = (TmpStateDescription >> unitCellCoordinate) & this->FirstUnitCellMask;
      switch (Tmp)
	{
	case 0xbul:
	  {
	    TmpStateDescription &= ~(0xful << unitCellCoordinate);
	    TmpStateDescription |= 0xeul << unitCellCoordinate;
	    return this->FindStateIndex(TmpStateDescription);
	  }
	  break;
	case 0x8ul:
	  {
	    TmpStateDescription &= ~(0xful << unitCellCoordinate);
	    TmpStateDescription |= 0x2ul << unitCellCoordinate;
	    return this->FindStateIndex(TmpStateDescription);
	  }
	  break;
	}
    }
  else
    {
      unitCellCoordinate *= this->UnitCellSize;
      unitCellCoordinate += siteIndex;
      unitCellCoordinate *= 2;
      unsigned long Tmp = (TmpStateDescription >> unitCellCoordinate) & this->FirstUnitCellMask;
      switch (Tmp)
	{
	case 0xeul:
	  {
	    TmpStateDescription &= ~(0xful << unitCellCoordinate);
	    TmpStateDescription |= 0xbul << unitCellCoordinate;
	    return this->FindStateIndex(TmpStateDescription);
	  }
	  break;
	case 0x2ul:
	  {
	    TmpStateDescription &= ~(0xful << unitCellCoordinate);
	    TmpStateDescription |= 0x8ul << unitCellCoordinate;
	    return this->FindStateIndex(TmpStateDescription);
	  }
	  break;
	}
    }
  return this->HilbertSpaceDimension;
}

// apply the operator coupling neighboring unit cells
//
// leftUnitCellCoordinate = coordinate of the left unit cell
// rightUnitCellCoordinate = coordinate of the right unit cell
// state = index of the state on which the operator has to be applied
// return value = index of resulting state 

int PairHoppingGenericPAsSpin1Chain::PlusMinusOperator (int leftUnitCellCoordinate, int rightUnitCellCoordinate, int state)
{
  unsigned long TmpStateDescription = this->StateDescription[state];
  leftUnitCellCoordinate *= 2 * this->UnitCellSize;
  leftUnitCellCoordinate += 2 * (this->UnitCellSize - 1);
  rightUnitCellCoordinate *= 2 * this->UnitCellSize;
  unsigned long Tmp = ((TmpStateDescription >> leftUnitCellCoordinate) & 0x3ul) | (((TmpStateDescription >> rightUnitCellCoordinate) & 0x3ul) << 2);
  if (Tmp == 0xaul)
    {
      TmpStateDescription &= ~(0x3ul << leftUnitCellCoordinate);
      TmpStateDescription |= 0x3ul << leftUnitCellCoordinate;
      TmpStateDescription &= ~(0x3ul << rightUnitCellCoordinate);
      TmpStateDescription |= 0x0ul << rightUnitCellCoordinate;
      return this->FindStateIndex(TmpStateDescription);
   }
  if (Tmp == 0x3ul)
    {
      TmpStateDescription &= ~(0x3ul << leftUnitCellCoordinate);
      TmpStateDescription |= 0x2ul << leftUnitCellCoordinate;
      TmpStateDescription &= ~(0x3ul << rightUnitCellCoordinate);
      TmpStateDescription |= 0x2ul << rightUnitCellCoordinate;            
      return this->FindStateIndex(TmpStateDescription);
    }
  return this->HilbertSpaceDimension;
}

// apply the operator coupling neighboring unit cells with a constraint of the unit cell parity
//
// leftUnitCellCoordinate = coordinate of the left unit cell
// rightUnitCellCoordinate = coordinate of the right unit cell
// state = index of the state on which the operator has to be applied
// return value = index of resulting state 

int PairHoppingGenericPAsSpin1Chain::PlusMinusOperatorPlus (int leftUnitCellCoordinate, int rightUnitCellCoordinate, int state)
{
  unsigned long TmpStateDescription = this->StateDescription[state];
  if ((rightUnitCellCoordinate & 1) == 0)
    {
      leftUnitCellCoordinate *= 2 * this->UnitCellSize;
      leftUnitCellCoordinate += 2 * (this->UnitCellSize - 1);
      rightUnitCellCoordinate *= 2 * this->UnitCellSize;
      unsigned long Tmp = ((TmpStateDescription >> leftUnitCellCoordinate) & 0x3ul) | (((TmpStateDescription >> rightUnitCellCoordinate) & 0x3ul) << 2);
      if (Tmp == 0xaul)
	{
	  TmpStateDescription &= ~(0x3ul << leftUnitCellCoordinate);
	  TmpStateDescription |= 0x3ul << leftUnitCellCoordinate;
	  TmpStateDescription &= ~(0x3ul << rightUnitCellCoordinate);
	  TmpStateDescription |= 0x0ul << rightUnitCellCoordinate;
	  return this->FindStateIndex(TmpStateDescription);
	}
    }
  else
    {
      leftUnitCellCoordinate *= this->UnitCellSize << 1;
      leftUnitCellCoordinate += (this->UnitCellSize - 1) << 1;
      rightUnitCellCoordinate *= this->UnitCellSize << 1;
      unsigned long Tmp = ((TmpStateDescription >> leftUnitCellCoordinate) & 0x3ul) | (((TmpStateDescription >> rightUnitCellCoordinate) & 0x3ul) << 2);
      if (Tmp == 0x3ul)
	{
	  TmpStateDescription &= ~(0x3ul << leftUnitCellCoordinate);
	  TmpStateDescription |= 0x2ul << leftUnitCellCoordinate;
	  TmpStateDescription &= ~(0x3ul << rightUnitCellCoordinate);
	  TmpStateDescription |= 0x2ul << rightUnitCellCoordinate;            
	  return this->FindStateIndex(TmpStateDescription);
	}
    }
  return this->HilbertSpaceDimension;
}


// generate all the posible unit cells 
//

void PairHoppingGenericPAsSpin1Chain::GenerateAllPossibleUnitCells()
{
  unsigned long TmpMaxNbrPossibleUnitCells = 0x3ul;
  for (int i = 1; i < this->UnitCellSize; ++i)
    {
      TmpMaxNbrPossibleUnitCells *= 0x3ul;
    }
  this->NbrPossibleUnitCells = 0; 
  this->PossibleUnitCells = new unsigned long[TmpMaxNbrPossibleUnitCells];
  this->PossibleUnitCellsNbrPlus = new int[TmpMaxNbrPossibleUnitCells];
  int* TmpPossibleUnitCellsNbrMinus = new int[TmpMaxNbrPossibleUnitCells];

  for (long i = ((0x1l << (this->UnitCellSize << 1)) - 0x1l); i >= 0x0l; --i)
    {
      bool PlusFlag = false;
      bool MinusFlag = false;
      bool DiscardFlag = false;
      int Tmp = 0;
      while ((PlusFlag == false) && (DiscardFlag == false) && (Tmp < this->UnitCellSize))
	{
	  long Tmp2 = ((i >> (Tmp << 1)) & 0x3l);
	  if (Tmp2 == 0x3l)
	    {
	      PlusFlag = true;
	    }
	  else
	    {
	      if (Tmp2 == 0x1l)
		{
		  DiscardFlag = true;
		}
	    }
	  ++Tmp;
	}
      if (PlusFlag == true)
	{
	  while ((MinusFlag == false) && (DiscardFlag == false) && (Tmp < this->UnitCellSize))
	    {
	      long Tmp2 = ((i >> (Tmp << 1)) & 0x3l);
	      if (Tmp2 == 0x0l)
		{
		  MinusFlag = true;
		}
	      else
		{
		  if (Tmp2 == 0x1l)
		    {
		      DiscardFlag = true;
		    }
		}
	      ++Tmp;
	    }
	}
      if ((MinusFlag == false) && (DiscardFlag == false))
	{
	  this->PossibleUnitCells[this->NbrPossibleUnitCells] = (unsigned long) i;
	  this->PossibleUnitCellsNbrPlus[this->NbrPossibleUnitCells] = 0;
	  TmpPossibleUnitCellsNbrMinus[this->NbrPossibleUnitCells] = 0;
	  for (Tmp = 0; Tmp < this->UnitCellSize; ++Tmp)
	    {
	      long Tmp2 = ((i >> (Tmp << 1)) & 0x3l);
	      if (Tmp2 == 0x3l)
		{
		  this->PossibleUnitCellsNbrPlus[this->NbrPossibleUnitCells]++;		  
		}
	      else
		{
		  if (Tmp2 == 0x0l)
		    {
		      TmpPossibleUnitCellsNbrMinus[this->NbrPossibleUnitCells]++;
		    }
		}
	    }
	  this->NbrPossibleUnitCells++;
	}
    }

  this->NbrPossibleUnitCellsPerNbrMinus = new int[this->UnitCellSize + 1];
  for (int i = 0; i <= this->UnitCellSize; ++i)
    {
      this->NbrPossibleUnitCellsPerNbrMinus[i] = 0;      
    }
  for (int i = 0; i < this->NbrPossibleUnitCells; ++i)
    {
      this->NbrPossibleUnitCellsPerNbrMinus[TmpPossibleUnitCellsNbrMinus[i]]++;
    }
  this->PossibleUnitCellsPerNbrMinus = new unsigned long*[this->UnitCellSize + 1];
  this->PossibleUnitCellsNbrPlusPerNbrMinus = new int*[this->UnitCellSize + 1];
  for (int i = 0; i <= this->UnitCellSize; ++i)
    {
      if (this->NbrPossibleUnitCellsPerNbrMinus[i] != 0)
	{
	  this->PossibleUnitCellsPerNbrMinus[i] = new unsigned long [this->NbrPossibleUnitCellsPerNbrMinus[i]];
	  this->PossibleUnitCellsNbrPlusPerNbrMinus[i] = new int [this->NbrPossibleUnitCellsPerNbrMinus[i]];
	}
      else
	{
	  this->PossibleUnitCellsPerNbrMinus[i] = 0;
	}
      this->NbrPossibleUnitCellsPerNbrMinus[i] = 0;     
    }
  for (int i = 0; i < this->NbrPossibleUnitCells; ++i)
    {
      this->PossibleUnitCellsPerNbrMinus[TmpPossibleUnitCellsNbrMinus[i]][this->NbrPossibleUnitCellsPerNbrMinus[TmpPossibleUnitCellsNbrMinus[i]]] = this->PossibleUnitCells[i];
      this->PossibleUnitCellsNbrPlusPerNbrMinus[TmpPossibleUnitCellsNbrMinus[i]][this->NbrPossibleUnitCellsPerNbrMinus[TmpPossibleUnitCellsNbrMinus[i]]] = this->PossibleUnitCellsNbrPlus[i];
      this->NbrPossibleUnitCellsPerNbrMinus[TmpPossibleUnitCellsNbrMinus[i]]++;
    }
  delete[] TmpPossibleUnitCellsNbrMinus;
}
