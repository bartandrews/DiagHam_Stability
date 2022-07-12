////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of pair hopping with generic p                  //
//            Hilbert space written as spin 1 chain with translations         //
//                   and the combination of inversion and spin flip           //
//                                                                            //
//                        last modification : 16/10/2019                      //
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


#include "HilbertSpace/PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry.h"
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

PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry::PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry () 
{
  this->InversionSector = 0.0;
  this->InversionShift = 0;
  this->InversionUnshift = 0;
  this->SzSymmetryMask = 0x0ul;
}

// constructor for Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// inversionSector = inversion symmetry sector (can be either +1 or -1)
// pValue = p value
// memory = amount of memory granted for precalculations

PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry::PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry (int chainLength, int momentum, int inversionSector, int pValue, unsigned long memory) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->Momentum = momentum;
  this->UnitCellSize = pValue;
  this->FirstUnitCellMask = (0x1ul << (this->UnitCellSize << 1)) - 0x1ul;
  this->GenerateAllPossibleUnitCells();

  if (inversionSector == 1)
    {
      this->InversionSector = 1.0;
    }
  else
    {
      this->InversionSector = -1.0;
    }
#ifdef __64_BITS__
  this->InversionShift = 32 - ((this->ChainLength >> 1) << 1);
#else
  this->InversionShift = 16 - ((this->ChainLength >> 1) << 1);
#endif
  if ((this->ChainLength & 1) == 0)
    this->InversionUnshift = this->InversionShift;
  else
    this->InversionUnshift = this->InversionShift - 2;
#ifdef __64_BITS__
  if (this-> ChainLength < 32)
#else
  if (this-> ChainLength < 16)
#endif
    {
      this->SzSymmetryMask = (0x1ul << (2 * this-> ChainLength)) - 0x1ul;
    }
  else
    {
      this->SzSymmetryMask  = ~0x0ul;
    }

  this->MaxXMomentum = this->ChainLength / this->UnitCellSize;
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

PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry::PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry (const PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry& chain) 
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

      this->InversionSector = chain.InversionSector;
      this->InversionShift = chain.InversionShift;
      this->InversionUnshift = chain.InversionUnshift;
      this->SzSymmetryMask = chain.SzSymmetryMask;
 
      this->MaxXMomentum = chain.MaxXMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;

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
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;

      this->InversionSector = 0.0;
      this->InversionShift = 0;
      this->InversionUnshift = 0;
      this->SzSymmetryMask = 0x0ul;

      this->MaxXMomentum = 0;
      this->StateXShift = 0;
      this->ComplementaryStateXShift = 0;
      this->XMomentumMask = 0x0ul;

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
      this->PossibleUnitCellsNbrMinus = 0;
    }
}

// destructor
//

PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry::~PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry& PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry::operator = (const PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry& chain)
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
	      delete[] this->PossibleUnitCellsNbrMinus;
	    }
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

      this->InversionSector = chain.InversionSector;
      this->InversionShift = chain.InversionShift;
      this->InversionUnshift = chain.InversionUnshift;
      this->SzSymmetryMask = chain.SzSymmetryMask;
 
      this->MaxXMomentum = chain.MaxXMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;

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
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;

      this->InversionSector = 0.0;
      this->InversionShift = 0;
      this->InversionUnshift = 0;
      this->SzSymmetryMask = 0x0ul;

      this->MaxXMomentum = 0;
      this->StateXShift = 0;
      this->ComplementaryStateXShift = 0;
      this->XMomentumMask = 0x0ul;

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
      this->PossibleUnitCellsNbrMinus = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry::Clone()
{
  return new PairHoppingGenericPAsSpin1ChainWithTranslationsAndInversionSzSymmetry (*this);
}

