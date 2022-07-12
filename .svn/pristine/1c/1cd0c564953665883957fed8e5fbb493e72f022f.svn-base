////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of fermions on lattice with spin  and Gutzwiller          //
//   projection in real space with translation invariance in two directions   //
//                            and Sz<->-Sz symmetry                           //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 03/11/2014                      //
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
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "GeneralTools/StringTools.h"

#include <math.h>
#include <cstdlib>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// xMomentum = momentum sector in the x direction
// maxXMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation (int nbrFermions, int nbrSite, bool minusSzParity, 
																		    int xMomentum, int maxXMomentum,
																		    int yMomentum, int maxYMomentum,
																		    unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalSpin = 0;
  this->SzParitySign = 1.0;
  this->SzParity = 0x0ul;
  if (minusSzParity == true)
    {
      this->SzParitySign = -1.0;
      this->SzParity = 0x1ul;
    }
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum = maxXMomentum;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = this->MaxXMomentum;

  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = 2 * (this->NbrSite / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this->MaxMomentum) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->MaxYMomentum =  maxYMomentum;
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = (2 * this->NbrSite) / this->StateXShift;
  this->StateYShift = 2 * (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
  this->YMomentumFullMask = 0x0ul;
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
    }
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 
  this->NbrFermionsParity = (~((unsigned long) this->NbrFermions)) & 0x1ul;

  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions);
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 31))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->GenerateSignLookUpTable();
      this->LargeHilbertSpaceDimension  = this->GenerateStates();
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  this->GenerateLookUpTable(memory);	  
//	  this->CheckHilbertSpace();
#ifdef __DEBUG__
	  long UsedMemory = 0;
	  UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
	  cout << "memory requested for Hilbert space = ";
	  if (UsedMemory >= 1024)
	    if (UsedMemory >= 1048576)
	      cout << (UsedMemory >> 20) << "Mo" << endl;
	    else
	      cout << (UsedMemory >> 10) << "ko" <<  endl;
	  else
	    cout << UsedMemory << endl;
	  UsedMemory = this->NbrMomentum * sizeof(int);
	  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
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
}


// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// xMomentum = momentum sector in the x direction
// maxXMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation (int nbrFermions, int totalSpin, int nbrSite, bool minusSzParity, 
																		    int xMomentum, int maxXMomentum,
																		    int yMomentum, int maxYMomentum,
																		    unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalSpin = 0;
  this->SzParitySign = 1.0;
  this->SzParity = 0x0ul;
  if (minusSzParity == true)
    {
      this->SzParitySign = -1.0;
      this->SzParity = 0x1ul;
    }
  this->NbrFermionsUp = (totalSpin + this->NbrFermions) >> 1;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum = maxXMomentum;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = this->MaxXMomentum;

  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = 2 * (this->NbrSite / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this->MaxMomentum) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->MaxYMomentum =  maxYMomentum;
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = (2 * this->NbrSite) / this->StateXShift;
  this->StateYShift = 2 * (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
  this->YMomentumFullMask = 0x0ul;
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
    }
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 
  this->NbrFermionsParity = (~((unsigned long) this->NbrFermions)) & 0x1ul;

  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFermionsUp);
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 31))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->GenerateSignLookUpTable();
      this->LargeHilbertSpaceDimension  = this->GenerateStates();
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  this->GenerateLookUpTable(memory);	  
//	  this->CheckHilbertSpace();
#ifdef __DEBUG__
	  long UsedMemory = 0;
	  UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
	  cout << "memory requested for Hilbert space = ";
	  if (UsedMemory >= 1024)
	    if (UsedMemory >= 1048576)
	      cout << (UsedMemory >> 20) << "Mo" << endl;
	    else
	      cout << (UsedMemory >> 10) << "ko" <<  endl;
	  else
	    cout << UsedMemory << endl;
	  UsedMemory = this->NbrMomentum * sizeof(int);
	  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
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
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation(const FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation& fermions)
{
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->SzParitySign = fermions.SzParitySign;
  this->SzParity = fermions.SzParity;
  this->NbrSite = fermions.NbrSite;

  this->MaxXMomentum = fermions.MaxXMomentum;
  this->XMomentum = fermions.XMomentum;
  this->StateXShift = fermions.StateXShift;
  this->ComplementaryStateXShift = fermions.ComplementaryStateXShift;
  this->XMomentumMask = fermions.XMomentumMask;
  this->MaxYMomentum = fermions.MaxYMomentum;
  this->YMomentum = fermions.YMomentum;
  this->NbrYMomentumBlocks = fermions.NbrYMomentumBlocks;
  this->StateYShift = fermions.StateYShift;
  this->YMomentumBlockSize = fermions.YMomentumBlockSize;
  this->ComplementaryStateYShift = fermions.ComplementaryStateYShift;
  this->YMomentumMask = fermions.YMomentumMask;
  this->YMomentumBlockMask = fermions.YMomentumBlockMask;  
  this->YMomentumFullMask = fermions.YMomentumFullMask;
  this->ComplementaryYMomentumFullMask = fermions.ComplementaryYMomentumFullMask; 

  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::~FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation& FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::operator = (const FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->SzParitySign = fermions.SzParitySign;
  this->SzParity = fermions.SzParity;
  this->NbrSite = fermions.NbrSite;

  this->MaxXMomentum = fermions.MaxXMomentum;
  this->XMomentum = fermions.XMomentum;
  this->StateXShift = fermions.StateXShift;
  this->ComplementaryStateXShift = fermions.ComplementaryStateXShift;
  this->XMomentumMask = fermions.XMomentumMask;
  this->MaxYMomentum = fermions.MaxYMomentum;
  this->YMomentum = fermions.YMomentum;
  this->NbrYMomentumBlocks = fermions.NbrYMomentumBlocks;
  this->StateYShift = fermions.StateYShift;
  this->YMomentumBlockSize = fermions.YMomentumBlockSize;
  this->ComplementaryStateYShift = fermions.ComplementaryStateYShift;
  this->YMomentumMask = fermions.YMomentumMask;
  this->YMomentumBlockMask = fermions.YMomentumBlockMask;  
  this->YMomentumFullMask = fermions.YMomentumFullMask;
  this->ComplementaryYMomentumFullMask = fermions.ComplementaryYMomentumFullMask; 

  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::Clone()
{
  return new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation(*this);
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::RawGenerateStates(int nbrFermions, int currentSite, long pos)
{
  return this->RawGenerateStatesWithHoleCounting(nbrFermions, currentSite, currentSite - nbrFermions + 1, pos);
}

// generate all states corresponding to the constraints with Sz conserved
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrSpinUp = number of fermions with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::RawGenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos)
{
  return this->RawGenerateStatesWithHoleCounting(nbrFermions, currentSite, currentSite - nbrFermions + 1, nbrSpinUp, pos);
}

// generate all states corresponding to the constraints, knowing the number of holes
// 
// nbrFermions = number of fermions
// currentSite = current site index in real state
// nbrHoles = number of unoccupied sites
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::RawGenerateStatesWithHoleCounting(int nbrFermions, int currentSite, int nbrHoles, long pos)
{
  if (nbrFermions == 0)
    {
      if (nbrHoles == (currentSite + 1))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else
	{
	  return pos;
	}
    }
  if (currentSite < 0)
    return pos;
  long TmpPos = this->RawGenerateStatesWithHoleCounting(nbrFermions - 1, currentSite - 1, nbrHoles, pos);
  unsigned long Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStatesWithHoleCounting(nbrFermions - 1, currentSite - 1, nbrHoles, pos);
   Mask = 0x1ul << ((currentSite) << 1);
   for (; pos < TmpPos; ++pos)
     this->StateDescription[pos] |= Mask;
   if (nbrHoles == 0)
     return pos;
   else
     return this->RawGenerateStatesWithHoleCounting(nbrFermions, currentSite - 1, nbrHoles - 1, pos);
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrHoles = number of unoccupied sites
// nbrSpinUp = number of fermions with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::RawGenerateStatesWithHoleCounting(int nbrFermions, int currentSite, int nbrHoles, int nbrSpinUp, long pos)
{
 if ((nbrFermions == 0) && (nbrSpinUp == 0))
    {
      if (nbrHoles == (currentSite + 1))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else
	{
	  return pos;
	}
    }
   if ((currentSite < 0) || (nbrFermions < 0) || (nbrSpinUp > nbrFermions) || (nbrSpinUp < 0))
    return pos;
   
  long TmpPos = this->RawGenerateStatesWithHoleCounting(nbrFermions - 1, currentSite - 1, nbrHoles, nbrSpinUp - 1, pos);
  unsigned long Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStatesWithHoleCounting(nbrFermions - 1, currentSite - 1, nbrHoles, nbrSpinUp, pos);
   Mask = 0x1ul << ((currentSite) << 1);
   for (; pos < TmpPos; ++pos)
     this->StateDescription[pos] |= Mask;
   if (nbrHoles == 0)
     return pos;
   else
     return this->RawGenerateStatesWithHoleCounting(nbrFermions, currentSite - 1, nbrHoles - 1, nbrSpinUp, pos); 
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::EvaluateHilbertSpaceDimension(int nbrFermions)
{
  BinomialCoefficients binomials(this->NbrSite);
  int NbrHoles = this->NbrSite - this->NbrFermions;
  long dimension = binomials(this->NbrSite, NbrHoles);
  for (int i = 0; i < this->NbrFermions; ++i)
    dimension *= 2l;
  return dimension;
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// nbrSpinUp = number of fermions with spin up
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation::EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp)
{
  BinomialCoefficients binomials(this->NbrSite);
  int NbrHoles = this->NbrSite - this->NbrFermions;
  long dimension = binomials(this->NbrSite, NbrHoles);
  BinomialCoefficients binomials1(this->NbrFermions);
  dimension *= binomials1(this->NbrFermions, this->NbrFermionsUp);
  return dimension;
}


