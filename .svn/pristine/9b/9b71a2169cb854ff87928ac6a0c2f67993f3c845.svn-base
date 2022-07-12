////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with translation invariance          //
//                           and the inversion symmetry                       //
//                                                                            //
//                        last modification : 27/07/2016                      //
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


#include "HilbertSpace/Spin1_2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1_2Chain.h"
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


// default constructor
//

Spin1_2ChainWithTranslationsAndInversionSymmetry::Spin1_2ChainWithTranslationsAndInversionSymmetry () 
{
  this->InversionSector = 0.0;
  this->InversionShift = 0;
  this->InversionUnshift = 0;
}

// constructor for Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// translationStep = indicates the step for an elementary translation
// inversionSector = inversion symmetry sector (can be either +1 or -1)
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

Spin1_2ChainWithTranslationsAndInversionSymmetry::Spin1_2ChainWithTranslationsAndInversionSymmetry (int chainLength, int momentum, int translationStep, int inversionSector, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->ComplementaryStateShift = this->ChainLength - translationStep;
  this->StateMask = (0x1ul << translationStep) - 1ul;
  this->StateShift = translationStep;
  this->Momentum = momentum;
  this->MaxXMomentum = chainLength;
  if (inversionSector == 1)
    {
      this->InversionSector = 1.0;
    }
  else
    {
      this->InversionSector = -1.0;
    }
#ifdef __64_BITS__
  this->InversionShift = 32 - ((((this->ChainLength + 1) >> 1)) * this->StateShift);
#else
  this->InversionShift = 16 - ((((this->ChainLength + 1) >> 1) - 1) * this->StateShift);
#endif
  if ((this->ChainLength & 1) == 0)
    this->InversionUnshift = this->InversionShift -  this->StateShift;
  else
    this->InversionUnshift = this->InversionShift;
  this->SzSymmetrySector = 0.0;
  this->SzSymmetryMask = (0x1ul << this->ChainLength) - 0x1ul;
  this->ComputeInversionTable();

  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->CreatePrecalculationTable();

  cout << "warning : untested code" << endl;
  long TmpHilbertSpaceDimension = (1l <<  this->ChainLength );
  this->LargeHilbertSpaceDimension = 0l;
  unsigned long TmpState;
  unsigned long TmpState2;
  int NbrTranslation;
  int CurrentNbrStateInOrbit;
  double TmpSign;
  for (long i = TmpHilbertSpaceDimension - 1l; i >= 0l; --i)
    {
      TmpState = (unsigned long) i;
      if (this->FindCanonicalForm(TmpState, NbrTranslation, TmpSign) == TmpState)
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      ++this->LargeHilbertSpaceDimension;
	    }
	}
    }
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int [this->LargeHilbertSpaceDimension];
  
  this->LargeHilbertSpaceDimension = 0l;
  for (long i = TmpHilbertSpaceDimension - 1l; i >= 0l; --i)
    {
      TmpState = (unsigned long) i;
      if (this->FindCanonicalForm(TmpState, NbrTranslation, TmpSign) == TmpState)
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      this->StateDescription[this->LargeHilbertSpaceDimension] = TmpState;
	      this->NbrStateInOrbit[this->LargeHilbertSpaceDimension] = this->FindOrbitSize(TmpState);
	      ++this->LargeHilbertSpaceDimension;
	    }
	}
    }
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;  

  if (this->HilbertSpaceDimension > 0)
    this->CreateLookUpTable();
}

// constructor for Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1/2
// momemtum = total momentum of each state
// translationStep = indicates the step for an elementary translation
// inversionSector = inversion symmetry sector (can be either +1 or -1)
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

Spin1_2ChainWithTranslationsAndInversionSymmetry::Spin1_2ChainWithTranslationsAndInversionSymmetry (int chainLength, int momentum, int translationStep, int inversionSector, int sz, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz;
  this->FixedSpinProjectionFlag = true;
  this->ComplementaryStateShift = this->ChainLength - translationStep;
  this->StateMask = (0x1ul << translationStep) - 1ul;
  this->StateShift = translationStep;
  this->Momentum = momentum;
  this->MaxXMomentum = chainLength;
  if (inversionSector == 1)
    {
      this->InversionSector = 1.0;
    }
  else
    {
      this->InversionSector = -1.0;
    }
#ifdef __64_BITS__
  this->InversionShift = 32 - ((((this->ChainLength + 1) >> 1)) * this->StateShift);
#else
  this->InversionShift = 16 - ((((this->ChainLength + 1) >> 1) - 1) * this->StateShift);
#endif
  if ((this->ChainLength & 1) == 0)
    this->InversionUnshift = this->InversionShift -  this->StateShift;
  else
    this->InversionUnshift = this->InversionShift;
  this->SzSymmetrySector = 0.0;
  this->SzSymmetryMask = (0x1ul << this->ChainLength) - 0x1ul;
  this->ComputeInversionTable();

  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->CreatePrecalculationTable();

  this->StateDescription = new unsigned long [this->EvaluateHilbertSpaceDimension(this->ChainLength, this->Sz)];
  long TmpHilbertSpaceDimension = this->GenerateStates(0l, this->ChainLength - 1, (this->ChainLength + this->Sz) >> 1);
  this->LargeHilbertSpaceDimension = 0l;
  unsigned long TmpState;
  unsigned long TmpState2;
  int NbrTranslation;
  int CurrentNbrStateInOrbit;
  unsigned long DicardFlag = ~0x0ul;
  double TmpSign;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      TmpState = this->StateDescription[i];
      TmpState2 = this->FindCanonicalForm(TmpState, NbrTranslation, TmpSign);
      if (TmpState2 == TmpState)
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      ++this->LargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescription[i] = DicardFlag;
	    }
	}
      else
	{
	  this->StateDescription[i] = DicardFlag;
	}
    }
  
  unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int [this->LargeHilbertSpaceDimension];
  
  this->LargeHilbertSpaceDimension = 0l;
  for (long i = TmpHilbertSpaceDimension - 1l; i >= 0l; --i)
    {
      if (this->StateDescription[i] != DicardFlag)
	{
	  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->StateDescription[i];
	  this->NbrStateInOrbit[this->LargeHilbertSpaceDimension] = this->FindOrbitSize(this->StateDescription[i]);
	  ++this->LargeHilbertSpaceDimension;
	}
    }
  
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;  
  this->LookUpTable =0;
  if (this->HilbertSpaceDimension > 0)
    this->CreateLookUpTable();
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainWithTranslationsAndInversionSymmetry::Spin1_2ChainWithTranslationsAndInversionSymmetry (const Spin1_2ChainWithTranslationsAndInversionSymmetry& chain) 
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
      this->Momentum = chain.Momentum;
      this->MaxXMomentum = chain.MaxXMomentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->StateMask = chain.StateMask;
      this->StateShift = chain.StateShift;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->SzSymmetrySector = chain.SzSymmetrySector;
      this->SzSymmetryMask = chain.SzSymmetryMask;
      this->InversionSector = chain.InversionSector;
      this->InversionShift = chain.InversionShift;
      this->InversionUnshift = chain.InversionUnshift;
      this->ComputeInversionTable();
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
      this->MaxXMomentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->StateMask = 0x0ul;
      this->StateShift = 0;
      this->ComplementaryStateShift = 0;
      this->SzSymmetrySector = 0.0;
      this->SzSymmetryMask = 0x0ul;
      this->InversionSector = 0.0;
      this->InversionShift = 0;
      this->InversionUnshift = 0;
    }
}

// destructor
//

Spin1_2ChainWithTranslationsAndInversionSymmetry::~Spin1_2ChainWithTranslationsAndInversionSymmetry () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainWithTranslationsAndInversionSymmetry& Spin1_2ChainWithTranslationsAndInversionSymmetry::operator = (const Spin1_2ChainWithTranslationsAndInversionSymmetry& chain)
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
      this->Momentum = chain.Momentum;
      this->MaxXMomentum = chain.MaxXMomentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->StateMask = chain.StateMask;
      this->StateShift = chain.StateShift;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->SzSymmetrySector = chain.SzSymmetrySector;
      this->SzSymmetryMask = chain.SzSymmetryMask;
      this->InversionSector = chain.InversionSector;
      this->InversionShift = chain.InversionShift;
      this->InversionUnshift = chain.InversionUnshift;
      this->ComputeInversionTable();
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
      this->MaxXMomentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->StateMask = 0x0ul;
      this->StateShift = 0;
      this->ComplementaryStateShift = 0;
      this->SzSymmetrySector = 0.0;
      this->SzSymmetryMask = 0x0ul;
      this->InversionSector = 0.0;
      this->InversionShift = 0;
      this->InversionUnshift = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainWithTranslationsAndInversionSymmetry::Clone()
{
  return new Spin1_2ChainWithTranslationsAndInversionSymmetry (*this);
}

// compute the inversion table
//

void Spin1_2ChainWithTranslationsAndInversionSymmetry::ComputeInversionTable()
{
  switch (this->StateShift)
    {
    case 1:
      {
	for (unsigned long i = 0x0ul; i <= 0xfful; ++i)
	  {
	    this->InversionTable[i]  = (i & 0x01ul) << 7;
	    this->InversionTable[i] |= (i & 0x02ul) << 5;
	    this->InversionTable[i] |= (i & 0x04ul) << 3;
	    this->InversionTable[i] |= (i & 0x08ul) << 1;
	    this->InversionTable[i] |= (i & 0x10ul) >> 1;
	    this->InversionTable[i] |= (i & 0x20ul) >> 3;
	    this->InversionTable[i] |= (i & 0x40ul) >> 5;
	    this->InversionTable[i] |= (i & 0x80ul) >> 7;
	  }
      }
      break;
    case 2:
      {
	for (unsigned long i = 0x0ul; i <= 0xfful; ++i)
	  {
	    this->InversionTable[i]  = (i & 0x03ul) << 6;
	    this->InversionTable[i] |= (i & 0x0cul) << 2;
	    this->InversionTable[i] |= (i & 0x30ul) >> 2;
	    this->InversionTable[i] |= (i & 0xc0ul) >> 6;
	  }
      }
      break;
    case 4:
      {
	for (unsigned long i = 0x0ul; i <= 0xfful; ++i)
	  {
	    this->InversionTable[i]  = (i & 0x0ful) << 4;
	    this->InversionTable[i]  = (i & 0xf0ul) >> 4;
	  }
      }
      break;
    default:
      {
	cout << "warning, Spin1_2ChainWithTranslationsAndInversionSymmetry does not support a block size of " << this->StateShift << endl;
	for (unsigned long i = 0x0ul; i <= 0xfful; ++i)
	  {
	    this->InversionTable[i]  = 0x0ul;
	  }	
      }
    }
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix Spin1_2ChainWithTranslationsAndInversionSymmetry::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
	{
	  RealMatrix TmpEntanglementMatrix(1, 1);
          double Tmp = 1.0;
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
      
    }
  if (nbrSites == this->ChainLength)
    {
      if (szSector == this->Sz)
	{
	  RealMatrix TmpEntanglementMatrix(1, 1);
	  double Tmp = 1.0;
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}      
    }
  Spin1_2Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1_2Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);

  RealMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  int TmpNbrTranslation;
  int TmpNbrTranslationToIdentity;
  double* TmpPhases = new double [2 * this->ChainLength];
  TmpPhases[0] = 1.0;
  if (this->Momentum == 0)
    {
      for (int i = 1; i < (2 * this->ChainLength); ++i)
	{
	  TmpPhases[i] = 1.0;
	}
    }
  else
    {
      for (int i = 1; i < (2 * this->ChainLength); ++i)
	{
	  TmpPhases[i] = -TmpPhases[i - 1];
	}
    }
  unsigned long Mask1 = (0x1ul << Shift) - 0x1ul;
  unsigned long Mask2 = (0x1ul << this->ChainLength) - 0x1ul;
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
	      TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos] * TmpPhases[TmpNbrTranslation] * Coefficient);
	    }
	}
    }
  delete[] TmpPhases;
  return TmpEntanglementMatrix;
}



