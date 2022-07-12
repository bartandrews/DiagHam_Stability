////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of Potts-3 chain with the translation symmetry            //
//                            and inversion symmetry                          //
//                                                                            //
//                        last modification : 10/12/2019                      //
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


#include "HilbertSpace/Potts3ChainWithTranslationsAndInversion.h"
#include "HilbertSpace/Potts3ChainWithTranslations.h"
#include "HilbertSpace/Potts3Chain.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;



// default constructor
//

Potts3ChainWithTranslationsAndInversion::Potts3ChainWithTranslationsAndInversion () 
{
  this->InversionSector = 1.0;
  this->InversionShift = 0;
  this->InversionUnshift = 0;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = momentum sector
// inversionSector = inversion symmetry sector (can be either +1 or -1)
// memorySize = memory size in bytes allowed for look-up table

Potts3ChainWithTranslationsAndInversion::Potts3ChainWithTranslationsAndInversion (int chainLength, int momentum, int inversionSector, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->ComplementaryStateShift = (this->ChainLength - 1) << 1;
  this->FixedQuantumNumberFlag = false;
  this->Momentum = momentum;
  this->MaxXMomentum = this->ChainLength;
  this->LargeHilbertSpaceDimension = 3l;
  for (int i = 1; i < chainLength; i++)
    this->LargeHilbertSpaceDimension *= 3l;
  if (inversionSector == 1)
    {
      this->InversionSector = 1.0;
    }
  else
    {
      this->InversionSector = -1.0;
    }
#ifdef __64_BITS__
  this->InversionShift = 32 - ((((this->ChainLength + 1) >> 1)) << 1);
#else
  this->InversionShift = 16 - ((((this->ChainLength + 1) >> 1) - 1) << 1);
#endif
  if ((this->ChainLength & 1) == 0)
    this->InversionUnshift = this->InversionShift - 2;
  else
    this->InversionUnshift = this->InversionShift;
  

  this->LookUpPosition = 0;
  this->LookUpTableSize = 4;
  memorySize >>= 2;
  this->LookUpTableMask = 0xfffffffc;
  while ((this->LookUpPosition <= this->ChainLength) && (memorySize >=  4))
    {
      this->LookUpTableMask <<= 2;
      this->LookUpTableSize <<= 2;
      memorySize >>= 2;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];

  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = this->RawGenerateStates (this->ChainLength - 1, 0);
  this->GenerateStates();
}

// constructor for complete Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// sz = twice the value of total Sz component
// momemtum = momentum sector
// inversionSector = inversion symmetry sector (can be either +1 or -1)
// memorySize = memory size in bytes allowed for look-up table

Potts3ChainWithTranslationsAndInversion::Potts3ChainWithTranslationsAndInversion (int chainLength, int sz, int momentum, int inversionSector, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->ComplementaryStateShift = (this->ChainLength - 1) << 1;
  this->Sz = sz % 3;
  this->FixedQuantumNumberFlag = true;
  this->Momentum = momentum;
  this->MaxXMomentum = this->ChainLength;
  if (inversionSector == 1)
    {
      this->InversionSector = 1.0;
    }
  else
    {
      this->InversionSector = -1.0;
    }
#ifdef __64_BITS__
  this->InversionShift = 32 - ((((this->ChainLength + 1) >> 1)) << 1);
#else
  this->InversionShift = 16 - ((((this->ChainLength + 1) >> 1) - 1) << 1);
#endif
  if ((this->ChainLength & 1) == 0)
    this->InversionUnshift = this->InversionShift - 2;
  else
    this->InversionUnshift = this->InversionShift;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->ChainLength - 1, this->Sz);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->LookUpPosition = 0;
  this->LookUpTableSize = 4;
  memorySize >>= 1;
  this->LookUpTableMask = 0xfffffffc;
  while ((this->LookUpPosition < this->ChainLength) && (memorySize >=  4))
    {
      this->LookUpTableMask <<= 2;
      this->LookUpTableSize <<= 2;
      memorySize >>= 2;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];

  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = this->RawGenerateStates (this->ChainLength - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->GenerateStates ();
}


// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Potts3ChainWithTranslationsAndInversion::Potts3ChainWithTranslationsAndInversion (const Potts3ChainWithTranslationsAndInversion& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->ChainDescription = chain.ChainDescription;
      this->Sz = chain.Sz;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->Momentum = chain.Momentum;
      this->MaxXMomentum = chain.MaxXMomentum;
      this->InversionSector = chain.InversionSector;
      this->InversionShift = chain.InversionShift;
      this->InversionUnshift = chain.InversionUnshift;
      this->RescalingFactors = chain.RescalingFactors;
   }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->NbrStateInOrbit = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
      this->LargeHilbertSpaceDimension = 0l;
      this->ComplementaryStateShift = 0;
      this->Momentum = 0;
      this->MaxXMomentum = 0;
      this->InversionSector = 0.0;
      this->InversionShift = 0;
      this->InversionUnshift = 0;
    }
}

// destructor
//

Potts3ChainWithTranslationsAndInversion::~Potts3ChainWithTranslationsAndInversion () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Potts3ChainWithTranslationsAndInversion& Potts3ChainWithTranslationsAndInversion::operator = (const Potts3ChainWithTranslationsAndInversion& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->ChainDescription;
      delete[] this->LookUpTable;
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainDescription = chain.ChainDescription;
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->LookUpPosition = chain.LookUpPosition;
      this->Sz = chain.Sz;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->Momentum = chain.Momentum;
      this->MaxXMomentum = chain.MaxXMomentum;
      this->InversionSector = chain.InversionSector;
      this->InversionShift = chain.InversionShift;
      this->InversionUnshift = chain.InversionUnshift;
      this->RescalingFactors = chain.RescalingFactors;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->NbrStateInOrbit = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
      this->LargeHilbertSpaceDimension = 0l;
      this->ComplementaryStateShift = 0;
      this->Momentum = 0;
      this->MaxXMomentum = 0;
      this->InversionSector = 0.0;
      this->InversionShift = 0;
      this->InversionUnshift = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Potts3ChainWithTranslationsAndInversion::Clone()
{
  return new Potts3ChainWithTranslationsAndInversion (*this);
}

// generate all states corresponding to a given momnetum
//

void Potts3ChainWithTranslationsAndInversion::GenerateStates()
{
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
#ifdef  __64_BITS__
  unsigned long Discard = 0xfffffffffffffffful;
#else
  unsigned long Discard = 0xfffffffful;
#endif
  double TmpSign;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->ChainDescription[i], NbrTranslationX, TmpSign) == this->ChainDescription[i]))
	{
	  if (this->TestMomentumConstraint(this->ChainDescription[i]) == true)
	    {
	      ++TmpLargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->ChainDescription[i] = Discard;
	    }
	}
      else
	{
	  this->ChainDescription[i] = Discard;
	}
    }
  if (TmpLargeHilbertSpaceDimension > 0)
    {
      unsigned long* TmpChainDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
      this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
      TmpLargeHilbertSpaceDimension = 0l;
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if (this->ChainDescription[i] != Discard)
	    {
	      TmpChainDescription[TmpLargeHilbertSpaceDimension] = this->ChainDescription[i];
	      this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = this->FindOrbitSize(this->ChainDescription[i]);
	      ++TmpLargeHilbertSpaceDimension;
	    }
	}
      delete[] this->ChainDescription;
      this->ChainDescription = TmpChainDescription;
    }
  else
    {
      delete[] this->ChainDescription;
    }

  this->LargeHilbertSpaceDimension = TmpLargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->RescalingFactors = new double* [2 * this->ChainLength + 1];
  for (int i = 1; i <= (2 * this->ChainLength); ++i)
    {
      this->RescalingFactors[i] = new double [2 * this->ChainLength + 1];
      for (int j = 1; j <= (2 * this->ChainLength); ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
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

RealMatrix Potts3ChainWithTranslationsAndInversion::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
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
  Potts3Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  int QB = this->Sz - szSector;
  if (QB < 0)
    {
      QB += 3;
    }
  Potts3Chain TmpHilbertSpace(this->ChainLength - nbrSites, QB, 1000000);

  RealMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = 2 * nbrSites;
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
  unsigned long Mask2 = (0x1ul << (2 * this->ChainLength)) - 0x1ul;
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

