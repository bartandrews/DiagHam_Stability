////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of spin 2 chain with translations                  //
//                          and the inversion symmetry                        //
//                                                                            //
//                        last modification : 12/11/2016                      //
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


#include "HilbertSpace/Spin2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin2Chain.h"
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

Spin2ChainWithTranslationsAndInversionSymmetry::Spin2ChainWithTranslationsAndInversionSymmetry () 
{
  this->InversionSector = 1.0;
  this->InversionShift = 0;
  this->InversionUnshift = 0;
}

// constructor for Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// inversionSector = inversion symmetry sector (can be either +1 or -1)
// memory = amount of memory granted for precalculations

Spin2ChainWithTranslationsAndInversionSymmetry::Spin2ChainWithTranslationsAndInversionSymmetry (int chainLength, int momentum, int inversionSector, unsigned long memory) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->Momentum = momentum;
  this->SzSymmetrySector = 0.0;
  if (inversionSector == 1)
    {
      this->InversionSector = 1.0;
    }
  else
    {
      this->InversionSector = -1.0;
    }
  this->InversionShift = 0;
  this->InversionUnshift = 60 - (this->ChainLength * 3);

  this->MaxXMomentum = this->ChainLength;
  this->StateXShift = 3 * (this->ChainLength / this->MaxXMomentum);
  this->ComplementaryStateXShift = (3 * this-> ChainLength) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;
  this->SzSymmetryMask = (0x1ul << (3 * this-> ChainLength)) - 0x1ul;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->ChainLength);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->BuildProjectors();
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

// constructor for Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component (should be equal to zero)
// inversionSector = inversion symmetry sector (can be either +1 or -1)
// memory = amount of memory granted for precalculations

Spin2ChainWithTranslationsAndInversionSymmetry::Spin2ChainWithTranslationsAndInversionSymmetry (int chainLength, int momentum, int inversionSector, int sz, unsigned long memory) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz;
  this->FixedSpinProjectionFlag = true;
  this->Momentum = momentum;
  this->SzSymmetrySector = 0.0;
  if (inversionSector == 1)
    {
      this->InversionSector = 1.0;
    }
  else
    {
      this->InversionSector = -1.0;
    }
  this->InversionShift = 0;
  this->InversionUnshift = 60 - (this->ChainLength * 3);

  this->SzSymmetryMask = (0x1ul << (3 * this-> ChainLength)) - 0x1ul;

  this->MaxXMomentum = this->ChainLength;
  this->StateXShift = 3 * (this->ChainLength / this->MaxXMomentum);
  this->ComplementaryStateXShift = (3 * this-> ChainLength) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->ChainLength);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1, 0);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->BuildProjectors();
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

Spin2ChainWithTranslationsAndInversionSymmetry::Spin2ChainWithTranslationsAndInversionSymmetry (const Spin2ChainWithTranslationsAndInversionSymmetry& chain) 
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

      this->SzSymmetryMask = chain.SzSymmetryMask;
      this->SzSymmetrySector = chain.SzSymmetrySector;
      this->InversionSector = chain.InversionSector;
      this->InversionShift = chain.InversionShift;
      this->InversionUnshift = chain.InversionUnshift;

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
 
      this->SzSymmetryMask = 0x0ul;
      this->SzSymmetrySector = 0.0;
      this->InversionSector = 0.0;
      this->InversionShift = 0;
      this->InversionUnshift = 0;

      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

Spin2ChainWithTranslationsAndInversionSymmetry::~Spin2ChainWithTranslationsAndInversionSymmetry () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin2ChainWithTranslationsAndInversionSymmetry& Spin2ChainWithTranslationsAndInversionSymmetry::operator = (const Spin2ChainWithTranslationsAndInversionSymmetry& chain)
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

      this->SzSymmetryMask = chain.SzSymmetryMask;
      this->SzSymmetrySector = chain.SzSymmetrySector;
      this->InversionSector = chain.InversionSector;
      this->InversionShift = chain.InversionShift;
      this->InversionUnshift = chain.InversionUnshift;

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

      this->SzSymmetryMask = 0x0ul;
      this->SzSymmetrySector = 0.0;
      this->InversionSector = 0.0;
      this->InversionShift = 0;
      this->InversionUnshift = 0;
     
      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin2ChainWithTranslationsAndInversionSymmetry::Clone()
{
  return new Spin2ChainWithTranslationsAndInversionSymmetry (*this);
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix Spin2ChainWithTranslationsAndInversionSymmetry::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
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
  Spin2Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin2Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);

  RealMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = 3 * nbrSites;
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
  unsigned long Mask2 = (0x1ul << (3 * this->ChainLength)) - 0x1ul;
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

