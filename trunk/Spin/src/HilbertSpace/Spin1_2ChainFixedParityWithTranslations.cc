////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with translastion invariance         //  
//                             and fixed parity                               //
//                                                                            //
//                        last modification : 06/06/2014                      //
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


#include "HilbertSpace/Spin1_2ChainFixedParityWithTranslations.h"
#include "MathTools/Complex.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "GeneralTools/ArrayTools.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;
using std::hex;
using std::dec;


#define M_SQRT3 1.73205080756888
#define M_SQRT6 2.44948974278318


// default constructor
//

Spin1_2ChainFixedParityWithTranslations::Spin1_2ChainFixedParityWithTranslations () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->ChainLength = 0;
  this->Momentum = 0;
  this->StateMask = 0x0ul;
  this->StateShift = 0;
  this->ComplementaryStateShift = 0;
  this->Sz = 0;
  this->SzParity = 0;
  this->FixedSpinProjectionFlag = false;
  this->CompatibilityWithMomentum = 0;
  this->RescalingFactors = 0;
  this->NbrStateInOrbit = 0;
}

// constructor for Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// translationStep = indicates the step for an elementary translation
// parity = parity of the total (Sz + 1/2) (can be 0 or 1)
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

Spin1_2ChainFixedParityWithTranslations::Spin1_2ChainFixedParityWithTranslations (int chainLength, int momentum, int translationStep, int parity, 
										  int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->SzParity = parity;
  this->Sz = 0;
  this->FixedSpinProjectionFlag = false;
  this->ComplementaryStateShift = this->ChainLength - translationStep;
  this->StateMask = (0x1ul << translationStep) - 1ul;
  this->StateShift = translationStep;
  this->Momentum = momentum;

  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->CreatePrecalculationTable();

  long TmpHilbertSpaceDimension = 1l << (this->ChainLength - 1);
  this->StateDescription = new unsigned long [TmpHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = 0l;
  unsigned long TmpState;
  unsigned long TmpState2;
  int NbrTranslation;
  int CurrentNbrStateInOrbit;
  unsigned long DicardFlag = ~0x0ul;
  unsigned long TmpMax = 1ul << this->ChainLength;
  unsigned long TmpSzParity = this->SzParity ^ (this->ChainLength & 1);
  for (unsigned long i = 0ul; i <= TmpMax; ++i)
    {
      unsigned long TmpState = i;
#ifdef __64_BITS__
      TmpState ^= TmpState >> 32;
#endif
      TmpState ^= TmpState >> 16;
      TmpState ^= TmpState >> 8;
      TmpState ^= TmpState >> 4;
      TmpState ^= TmpState >> 2;
      TmpState ^= TmpState >> 1;
      if ((TmpState & 0x1ul) == TmpSzParity)
	{
	  TmpState = i;
	  TmpState2 = this->FindCanonicalForm(TmpState, NbrTranslation);
	  if (TmpState2 == TmpState)
	    {
	      CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpState2);
	      if (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true)
		{
		  this->StateDescription[this->LargeHilbertSpaceDimension] = TmpState2;
		  ++this->LargeHilbertSpaceDimension;
		}
	    }
	}
    }

  unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int [this->LargeHilbertSpaceDimension];
  for (long i = this->LargeHilbertSpaceDimension - 1l; i >= 0l; --i)
    {
      TmpStateDescription[i] = this->StateDescription[i];
      CurrentNbrStateInOrbit = this->FindNumberTranslation(this->StateDescription[i]);
      this->NbrStateInOrbit[i] = CurrentNbrStateInOrbit;
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;  
  if (this->HilbertSpaceDimension > 0)
    this->CreateLookUpTable();
}

// constructor from pre-constructed datas
//
// hilbertSpaceDimension = Hilbert space dimension
// chainDescription = array describing states
// chainLength = number of spin 1
// momemtum = total momentum of each state
// parity = parity of the total (Sz + 1/2) (can be 0 or 1)
// lookUpTableShift = shift to apply to a state to obtain an index to the look-up table 
// complementaryStateShift = shift to apply to move the spin from one end to the other one

Spin1_2ChainFixedParityWithTranslations::Spin1_2ChainFixedParityWithTranslations (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
										  int momentum, int parity, int lookUpTableShift, 
										  int complementaryStateShift)
{
  this->Flag.Initialize();
  this->ComplementaryStateShift = complementaryStateShift;
  this->LookUpTableShift = lookUpTableShift;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->StateDescription = chainDescription;
  this->Momentum = momentum;
  this->Sz = 0;
  this->SzParity = parity;
  this->FixedSpinProjectionFlag = false;
  this->ChainLength = chainLength;
  this->CreatePrecalculationTable();
  this->CreateLookUpTable();
}
  
// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainFixedParityWithTranslations::Spin1_2ChainFixedParityWithTranslations (const Spin1_2ChainFixedParityWithTranslations& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->SzParity = chain.SzParity;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->StateMask = chain.StateMask;
      this->StateShift = chain.StateShift;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->Sz = 0;
      this->SzParity = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->StateMask = 0x0ul;
      this->StateShift = 0;
      this->ComplementaryStateShift = 0;
    }
}

// destructor
//

Spin1_2ChainFixedParityWithTranslations::~Spin1_2ChainFixedParityWithTranslations () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainFixedParityWithTranslations& Spin1_2ChainFixedParityWithTranslations::operator = (const Spin1_2ChainFixedParityWithTranslations& chain)
{
  if ((this->ChainLength != 0) && (this->HilbertSpaceDimension > 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTable;
      delete[] this->CompatibilityWithMomentum;
      int TmpPeriodicity = this->ChainLength / this->StateShift;
      for (int i = 1; i <= TmpPeriodicity; ++i)
	{
	  delete[] this->RescalingFactors[i];
	} 
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->SzParity = chain.SzParity;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->StateMask = chain.StateMask;
      this->StateShift = chain.StateShift;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
   }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->Sz = 0;
      this->SzParity = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->StateMask = 0x0ul;
      this->StateShift = 0;
      this->ComplementaryStateShift = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainFixedParityWithTranslations::Clone()
{
  return new Spin1_2ChainFixedParityWithTranslations (*this);
}

