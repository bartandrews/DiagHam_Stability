////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of pair hopping p=2                        //
//            Hilbert space written as spin 1 chain with translations         //
//                                                                            //
//                        last modification : 13/03/2019                      //
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


#include "HilbertSpace/PairHoppingP2AsSpin1ChainWithTranslations.h"
#include "HilbertSpace/PairHoppingP2AsSpin1Chain.h"
#include "HilbertSpace/Spin1Chain.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <math.h>

using std::cout;
using std::endl;
using std::dec;
using std::hex;


// default constructor
//

PairHoppingP2AsSpin1ChainWithTranslations::PairHoppingP2AsSpin1ChainWithTranslations () 
{
}

// constructor for Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// memory = amount of memory granted for precalculations

PairHoppingP2AsSpin1ChainWithTranslations::PairHoppingP2AsSpin1ChainWithTranslations (int chainLength, int momentum, unsigned long memory) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->Momentum = momentum;

  this->MaxXMomentum = this->ChainLength >> 1;
  this->StateXShift = 2 * (this->ChainLength / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this-> ChainLength) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(0, 0, 0);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, 0, 0, 0);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memory);
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

PairHoppingP2AsSpin1ChainWithTranslations::PairHoppingP2AsSpin1ChainWithTranslations (const PairHoppingP2AsSpin1ChainWithTranslations& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->MaxXMomentum = chain.MaxXMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;

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
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;

      this->MaxXMomentum = 0;
      this->StateXShift = 0;
      this->ComplementaryStateXShift = 0;
      this->XMomentumMask = 0x0ul;

      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;
    }
}

// destructor
//

PairHoppingP2AsSpin1ChainWithTranslations::~PairHoppingP2AsSpin1ChainWithTranslations () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

PairHoppingP2AsSpin1ChainWithTranslations& PairHoppingP2AsSpin1ChainWithTranslations::operator = (const PairHoppingP2AsSpin1ChainWithTranslations& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  delete[] this->StateDescription;
	  delete[] this->LookUpTable;
	  for (int i = 1; i <= this->ChainLength; ++i)
	    {
	      delete[] this->RescalingFactors[i];
	    } 
	  delete[] this->RescalingFactors;
	  delete[] this->NbrStateInOrbit;
	}
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->MaxXMomentum = chain.MaxXMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;

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
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;

      this->MaxXMomentum = 0;
      this->StateXShift = 0;
      this->ComplementaryStateXShift = 0;
      this->XMomentumMask = 0x0ul;

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

AbstractHilbertSpace* PairHoppingP2AsSpin1ChainWithTranslations::Clone()
{
  return new PairHoppingP2AsSpin1ChainWithTranslations(*this);
}

// apply the swap operator within the unit cell
//
// unitCellCoordinate = coordinate of the unit cell
// siteIndex = index of the leftmost site within the unit cell
// state = index of the state on which the operator has to be applied
// return value = index of resulting state 

int PairHoppingP2AsSpin1ChainWithTranslations::SwapOperator (int unitCellCoordinate, int siteIndex, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpStateDescription = this->StateDescription[state];
  unitCellCoordinate <<= 2;
  unsigned long Tmp = (TmpStateDescription >> unitCellCoordinate) & 0xful;
  switch (Tmp)
    {
    case 0xbul:
      {
	TmpStateDescription &= ~(0xful << unitCellCoordinate);
	TmpStateDescription |= 0xeul << unitCellCoordinate;
	coefficient = 1.0;
	return this->SymmetrizeResult(TmpStateDescription, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    case 0xeul:
      {
	TmpStateDescription &= ~(0xful << unitCellCoordinate);
	TmpStateDescription |= 0xbul << unitCellCoordinate;
	coefficient = 1.0;
	return this->SymmetrizeResult(TmpStateDescription, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    case 0x2ul:
      {
	TmpStateDescription &= ~(0xful << unitCellCoordinate);
	TmpStateDescription |= 0x8ul << unitCellCoordinate;
	coefficient = 1.0;
	return this->SymmetrizeResult(TmpStateDescription, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    case 0x8ul:
      {
	TmpStateDescription &= ~(0xful << unitCellCoordinate);
	TmpStateDescription |= 0x2ul << unitCellCoordinate;
	coefficient = 1.0;
	return this->SymmetrizeResult(TmpStateDescription, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    }  
  return this->HilbertSpaceDimension;
}

// apply the operator coupling neighboring unit cells
//
// leftUnitCellCoordinate = coordinate of the left unit cell
// rightUnitCellCoordinate = coordinate of the right unit cell
// state = index of the state on which the operator has to be applied
// return value = index of resulting state 

int PairHoppingP2AsSpin1ChainWithTranslations::PlusMinusOperator (int leftUnitCellCoordinate, int rightUnitCellCoordinate, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpStateDescription = this->StateDescription[state];
  leftUnitCellCoordinate <<= 2;
  leftUnitCellCoordinate += 2;
  rightUnitCellCoordinate <<= 2;
  unsigned long Tmp = ((TmpStateDescription >> leftUnitCellCoordinate) & 0x3ul) | (((TmpStateDescription >> rightUnitCellCoordinate) & 0x3ul) << 2);
  if (Tmp == 0xaul)
    {
      TmpStateDescription &= ~(0x3ul << leftUnitCellCoordinate);
      TmpStateDescription |= 0x3ul << leftUnitCellCoordinate;
      TmpStateDescription &= ~(0x3ul << rightUnitCellCoordinate);
      TmpStateDescription |= 0x0ul << rightUnitCellCoordinate;
      coefficient = 1.0;
      return this->SymmetrizeResult(TmpStateDescription, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
   }
  if (Tmp == 0x3ul)
    {
      TmpStateDescription &= ~(0x3ul << leftUnitCellCoordinate);
      TmpStateDescription |= 0x2ul << leftUnitCellCoordinate;
      TmpStateDescription &= ~(0x3ul << rightUnitCellCoordinate);
      TmpStateDescription |= 0x2ul << rightUnitCellCoordinate;            
      coefficient = 1.0;
      return this->SymmetrizeResult(TmpStateDescription, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
    }
  return this->HilbertSpaceDimension;
}

// evaluate Hilbert space dimension
//
// previousNbrPlus = number of +1 spins in the previous unit cell
// initialNbrMinus = number of -1 spins in the first unit cell
// sitePosition = site on chain where spin has to be changed
// return value = Hilbert space dimension

long PairHoppingP2AsSpin1ChainWithTranslations::EvaluateHilbertSpaceDimension(int previousNbrPlus, int initialNbrMinus, int sitePosition)
{
  if (sitePosition == this->ChainLength)
    {
      if (previousNbrPlus == initialNbrMinus)
	{
	  return 1l;
	}
      else
	{
	  return 0l;
	}
    }
  long TmpDimension = 0l;
  if (sitePosition != 0)
    {
      if (previousNbrPlus == 2)
	{
	  TmpDimension += this->EvaluateHilbertSpaceDimension(0, initialNbrMinus, sitePosition + 2);
	}
      else
	{
	  if (previousNbrPlus == 1)
	    {
	      TmpDimension += this->EvaluateHilbertSpaceDimension(0, initialNbrMinus, sitePosition + 2);
	      TmpDimension += this->EvaluateHilbertSpaceDimension(0, initialNbrMinus, sitePosition + 2);
	      TmpDimension += this->EvaluateHilbertSpaceDimension(1, initialNbrMinus, sitePosition + 2);
	    }
	  else
	    {
	      TmpDimension += this->EvaluateHilbertSpaceDimension(0, initialNbrMinus, sitePosition + 2);
	      TmpDimension += this->EvaluateHilbertSpaceDimension(1, initialNbrMinus, sitePosition + 2);
	      TmpDimension += this->EvaluateHilbertSpaceDimension(1, initialNbrMinus, sitePosition + 2);
	      TmpDimension += this->EvaluateHilbertSpaceDimension(2, initialNbrMinus, sitePosition + 2);
	    }
	}
    }
  else
    {
      TmpDimension += this->EvaluateHilbertSpaceDimension(0, 2, 2);
      TmpDimension += this->EvaluateHilbertSpaceDimension(0, 1, 2);
      TmpDimension += this->EvaluateHilbertSpaceDimension(0, 1, 2);
      TmpDimension += this->EvaluateHilbertSpaceDimension(0, 0, 2);
      TmpDimension += this->EvaluateHilbertSpaceDimension(1, 1, 2);
      TmpDimension += this->EvaluateHilbertSpaceDimension(1, 0, 2);
      TmpDimension += this->EvaluateHilbertSpaceDimension(1, 0, 2);
      TmpDimension += this->EvaluateHilbertSpaceDimension(2, 0, 2);
    }
  return TmpDimension;
}

// generate all states with neither constraint from boundary conditions nor discrete symmtry constraint
//
// statePosition = position for the new states
// previousNbrPlus = number of +1 spins in the previous unit cell
// initialNbrMinus = number of -1 spins in the first unit cell
// sitePosition = site on chain where spin has to be changed
// return value = number of generated states

long PairHoppingP2AsSpin1ChainWithTranslations::RawGenerateStates(long statePosition, int previousNbrPlus, int initialNbrMinus, int sitePosition) 
{
  if (sitePosition == this->ChainLength)
    {
      if (previousNbrPlus == initialNbrMinus)
	{
	  this->StateDescription[statePosition] = 0x0ul;
	  return (statePosition + 1l);
	}
      else
	{
	  return statePosition;
	}
    }
  if (sitePosition != 0)
    {
      if (previousNbrPlus == 2)
	{
	  statePosition = this->RawGenerateStates(statePosition, 0, initialNbrMinus, sitePosition + 2);
	}
      else
	{
	  if (previousNbrPlus == 1)
	    {
	      unsigned long TmpMask = 0x2ul << (sitePosition * 2);
	      long TmpPosition = this->RawGenerateStates(statePosition, 0, initialNbrMinus, sitePosition + 2);
	      for (; statePosition < TmpPosition; ++statePosition)
		this->StateDescription[statePosition] |= TmpMask;
	      TmpMask = 0x8ul << (sitePosition * 2);
	      TmpPosition = this->RawGenerateStates(statePosition, 0, initialNbrMinus, sitePosition + 2);
	      for (; statePosition < TmpPosition; ++statePosition)
		this->StateDescription[statePosition] |= TmpMask;
	      TmpMask = 0xcul << (sitePosition * 2);
	      TmpPosition = this->RawGenerateStates(statePosition, 1, initialNbrMinus, sitePosition + 2);
	      for (; statePosition < TmpPosition; ++statePosition)
		this->StateDescription[statePosition] |= TmpMask;
	    }
	  else
	    {
	      unsigned long TmpMask = 0xaul << (sitePosition * 2);
	      long TmpPosition = this->RawGenerateStates(statePosition, 0, initialNbrMinus, sitePosition + 2);
	      for (; statePosition < TmpPosition; ++statePosition)
		this->StateDescription[statePosition] |= TmpMask;
	      TmpMask = 0xbul << (sitePosition * 2);
	      TmpPosition = this->RawGenerateStates(statePosition, 1, initialNbrMinus, sitePosition + 2);
	      for (; statePosition < TmpPosition; ++statePosition)
		this->StateDescription[statePosition] |= TmpMask;
	      TmpMask = 0xeul << (sitePosition * 2);
	      TmpPosition = this->RawGenerateStates(statePosition, 1, initialNbrMinus, sitePosition + 2);
	      for (; statePosition < TmpPosition; ++statePosition)
		this->StateDescription[statePosition] |= TmpMask;
	      TmpMask = 0xful << (sitePosition * 2);
	      TmpPosition = this->RawGenerateStates(statePosition, 2, initialNbrMinus, sitePosition + 2);
	      for (; statePosition < TmpPosition; ++statePosition)
		this->StateDescription[statePosition] |= TmpMask;
	    }
	}
    }
  else
    {
      statePosition = this->RawGenerateStates(statePosition, 0, 2, 2);
      unsigned long TmpMask = 0x2ul << (sitePosition * 2);
      long TmpPosition = this->RawGenerateStates(statePosition, 0, 1, 2);
      for (; statePosition < TmpPosition; ++statePosition)
	this->StateDescription[statePosition] |= TmpMask;
      TmpMask = 0x8ul << (sitePosition * 2);
      TmpPosition = this->RawGenerateStates(statePosition, 0, 1, 2);
      for (; statePosition < TmpPosition; ++statePosition)
	this->StateDescription[statePosition] |= TmpMask;
      TmpMask = 0xaul << (sitePosition * 2);
      TmpPosition = this->RawGenerateStates(statePosition, 0, 0, 2);
      for (; statePosition < TmpPosition; ++statePosition)
	this->StateDescription[statePosition] |= TmpMask;
      TmpMask = 0xbul << (sitePosition * 2);
      TmpPosition = this->RawGenerateStates(statePosition, 1, 0, 2);
      for (; statePosition < TmpPosition; ++statePosition)
	this->StateDescription[statePosition] |= TmpMask;
      TmpMask = 0xcul << (sitePosition * 2);
      TmpPosition = this->RawGenerateStates(statePosition, 1, 1, 2);
      for (; statePosition < TmpPosition; ++statePosition)
	this->StateDescription[statePosition] |= TmpMask;
      TmpMask = 0xeul << (sitePosition * 2);
      TmpPosition = this->RawGenerateStates(statePosition, 1, 0, 2);
      for (; statePosition < TmpPosition; ++statePosition)
	this->StateDescription[statePosition] |= TmpMask;
      TmpMask = 0xful << (sitePosition * 2);
      TmpPosition = this->RawGenerateStates(statePosition, 2, 0, 2);
      for (; statePosition < TmpPosition; ++statePosition)
	this->StateDescription[statePosition] |= TmpMask;
    }
  return statePosition;
}

// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long PairHoppingP2AsSpin1ChainWithTranslations::GenerateStates()
{
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  double InversionSzSymmetrySign;
  unsigned long TmpDiscardMask = ~0x0ul;

  int TmpNbrMinusTable[] = {2, -1, 1, -1, -1, -1, -1, -1, 1, -1, 0, 0, 1, -1, 0, 0};
  int TmpNbrPlusTable[] =  {0, -1, 0, -1, -1, -1, -1, -1, 0, -1, 0, 1, 1, -1, 1, 2};
  int TmpLastUnitCellShift = (this->ChainLength - 2) * 2;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (TmpNbrPlusTable[(this->StateDescription[i] >> TmpLastUnitCellShift) & 0xf] != TmpNbrMinusTable[this->StateDescription[i]  & 0xf])
	{
	  this->StateDescription[i] = TmpDiscardMask;
	}
      else
	{
	  if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, InversionSzSymmetrySign) == this->StateDescription[i]))
	    {
	      if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
		{
		  ++TmpLargeHilbertSpaceDimension;
		}
	      else
		{
		  this->StateDescription[i] = TmpDiscardMask;
		}
	    }
	  else
	    {
	      this->StateDescription[i] = TmpDiscardMask;
	    }
	}
    }
  
  if (TmpLargeHilbertSpaceDimension > 0l)
    {
      unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];
      this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
      TmpLargeHilbertSpaceDimension = 0l;
      for  (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if ( this->StateDescription[i] != TmpDiscardMask)
	    {
	      TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	      ++TmpLargeHilbertSpaceDimension;
	    }
	}
      SortArrayDownOrdering<unsigned long>(TmpStateDescription, TmpLargeHilbertSpaceDimension);     
      delete[] this->StateDescription;
      this->StateDescription = TmpStateDescription;
      for  (long i = 0; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  this->NbrStateInOrbit[i] = this->FindOrbitSize(this->StateDescription[i]);
	}
    }
  return TmpLargeHilbertSpaceDimension;
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated  (fixing the boundray unit cells to nbr left minuses= szSector / (pValue + 1), nbr right pluses = szSector % (pValue + 1))
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix PairHoppingP2AsSpin1ChainWithTranslations::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, 1);
	  Complex Tmp(1.0, 0.0);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  int TmpNbrMinusLeft = szSector / 3;
  int TmpNbrPlusRight = szSector % 3;
  PairHoppingP2AsSpin1Chain TmpDestinationHilbertSpace(nbrSites, TmpNbrPlusRight, TmpNbrMinusLeft, 1000000);
  PairHoppingP2AsSpin1Chain TmpHilbertSpace(this->ChainLength - nbrSites, TmpNbrMinusLeft, TmpNbrPlusRight, 1000000);

  RealMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = 2 * nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  int TmpNbrTranslation;
  int TmpNbrTranslationToIdentity;
  double* TmpPhases = new double [2 * this->MaxXMomentum];
  double Coef = 2.0 * M_PI * ((double) this->Momentum) / ((double) this->MaxXMomentum);
  for (int i = 0; i < (2 * this->MaxXMomentum); ++i)
    {
      TmpPhases[i] = cos(Coef * ((double) i));
    }

  unsigned long Mask1 = (0x1ul << Shift) - 0x1ul;
  unsigned long Mask2;
#ifdef  __64_BITS__
  if (this->ChainLength < 32)
#else
  if (this->ChainLength < 16)   
#endif
    {
      Mask2 = (0x1ul << (2 * this->ChainLength)) - 0x1ul;      
    }
  else
    {
      Mask2 = ~0x0ul;
    }
  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = (TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask1)) & Mask2;
	  double Coefficient = 1.0;
	  int TmpPos = this->SymmetrizeResult(TmpState2, 1, Coefficient, TmpNbrTranslation);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos] * TmpPhases[TmpNbrTranslation]  * Coefficient);
	    }
	}
    }
  delete[] TmpPhases;
  return TmpEntanglementMatrix;
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated  (fixing the boundray unit cells to nbr left minuses= szSector / (pValue + 1), nbr right pluses = szSector % (pValue + 1))
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix PairHoppingP2AsSpin1ChainWithTranslations::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, 1);
	  Complex Tmp(1.0, 0.0);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  int TmpNbrMinusLeft = szSector / 3;
  int TmpNbrPlusRight = szSector % 3;
  PairHoppingP2AsSpin1Chain TmpDestinationHilbertSpace(nbrSites, TmpNbrPlusRight, TmpNbrMinusLeft, 1000000);
  PairHoppingP2AsSpin1Chain TmpHilbertSpace(this->ChainLength - nbrSites, TmpNbrMinusLeft, TmpNbrPlusRight, 1000000);

  ComplexMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = 2 * nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  int TmpNbrTranslation;
  int TmpNbrTranslationToIdentity;
  Complex* TmpPhases = new Complex [2 * this->MaxXMomentum];
  double Coef = 2.0 * M_PI * ((double) this->Momentum) / ((double) this->MaxXMomentum);
  for (int i = 0; i < (2 * this->MaxXMomentum); ++i)
    {
      TmpPhases[i] = Phase(Coef * ((double) i));
    }

  unsigned long Mask1 = (0x1ul << Shift) - 0x1ul;
  unsigned long Mask2;
#ifdef  __64_BITS__
  if (this->ChainLength < 32)
#else
  if (this->ChainLength < 16)   
#endif
    {
      Mask2 = (0x1ul << (2 * this->ChainLength)) - 0x1ul;      
    }
  else
    {
      Mask2 = ~0x0ul;
    }
  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = (TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask1)) & Mask2;
	  double Coefficient = 1.0;
	  int TmpPos = this->SymmetrizeResult(TmpState2, 1, Coefficient, TmpNbrTranslation);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos] * TmpPhases[TmpNbrTranslation]  * Coefficient);
	    }
	}
    }
  delete[] TmpPhases;
  return TmpEntanglementMatrix;
}
