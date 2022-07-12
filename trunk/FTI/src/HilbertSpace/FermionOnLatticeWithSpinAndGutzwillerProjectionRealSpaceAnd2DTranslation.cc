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
//                                                                            //
//                       class author: Nicolas Regnault                       //
//                                                                            //
//                        last modification : 20/08/2014                      //
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
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
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
#include "Architecture/ArchitectureOperation/FQHETorusParticleEntanglementSpectrumOperation.h"
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

FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation ()
{
  this->NbrFermions = 0;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalSpin = 0;
  this->MaxMomentum = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = 0;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = 0;
  this->XMomentum = 0; 
  this->YMomentum = 0;
  this->StateShift = 2 * (this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((unsigned long) 1);
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateHighestBit = 0;  
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// xMomentum = momentum sector in the x direction
// maxXMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (int nbrFermions, int nbrSite, 
																		    int xMomentum, int maxXMomentum,
																		    int yMomentum, int maxYMomentum,
																		    unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalSpin = 0;
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
  this->NbrSitePerUnitCell = this->NbrSite /  (this->MaxYMomentum * this->MaxXMomentum);
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
  this->TargetSpace = this;
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
// xMomentum = momentum sector in the x direction
// maxXMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation (int nbrFermions, int totalSpin, int nbrSite, 
																		    int xMomentum, int maxXMomentum,
																		    int yMomentum, int maxYMomentum,
																		    unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalSpin = 0;
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
  this->NbrSitePerUnitCell = this->NbrSite /  (this->MaxYMomentum * this->MaxXMomentum);
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
  this->TargetSpace = this;
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

FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation(const FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation& fermions)
{
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrSite = fermions.NbrSite;

  this->NbrSitePerUnitCell = fermions.NbrSitePerUnitCell;
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
  if (fermions.TargetSpace != &fermions)
    {
      this->TargetSpace = fermions.TargetSpace;
    }
  else
    {
      this->TargetSpace = this;
    }

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

FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::~FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation& FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::operator = (const FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation& fermions)
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
  this->NbrSite = fermions.NbrSite;

  this->NbrSitePerUnitCell = fermions.NbrSitePerUnitCell;
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
  if (fermions.TargetSpace != &fermions)
    {
      this->TargetSpace = fermions.TargetSpace;
    }
  else
    {
      this->TargetSpace = this;
    }

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

AbstractHilbertSpace* FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::Clone()
{
  return new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation(*this);
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::RawGenerateStates(int nbrFermions, int currentSite, long pos)
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

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::RawGenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos)
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

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::RawGenerateStatesWithHoleCounting(int nbrFermions, int currentSite, int nbrHoles, long pos)
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

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::RawGenerateStatesWithHoleCounting(int nbrFermions, int currentSite, int nbrHoles, int nbrSpinUp, long pos)
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

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::EvaluateHilbertSpaceDimension(int nbrFermions)
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

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp)
{
  BinomialCoefficients binomials(this->NbrSite);
  int NbrHoles = this->NbrSite - this->NbrFermions;
  long dimension = binomials(this->NbrSite, NbrHoles);
  BinomialCoefficients binomials1(this->NbrFermions);
  dimension *= binomials1(this->NbrFermions, this->NbrFermionsUp);
  return dimension;
}


// convert a state defined in the real space basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::ConvertToKxKyBasis(ComplexVector& state, ParticleOnSphere* space)  
{
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace* TmpSpace = (FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) space;
  ComplexVector TmpVector (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      int TmpMaxMomentum = this->StateHighestBit[i];
      int Pos = TmpSpace->FindStateIndex(TmpState, TmpMaxMomentum);
      if (Pos < TmpSpace->HilbertSpaceDimension)
	{
	  TmpVector[i] =  state[Pos] * sqrt((double) this->NbrStateInOrbit[i]);
	}
    }
  return TmpVector;
}


// convert a state defined in the (Kx,Ky) basis into a state in the real space basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space)
{
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace* TmpSpace = (FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) space;
  ComplexVector TmpVector (TmpSpace->HilbertSpaceDimension, true);
  Complex** FourrierCoefficients = new Complex* [this->MomentumModulo];
  for (int i = 0; i < this->MaxXMomentum; ++i)
  {
    FourrierCoefficients[i] = new Complex [this->MaxYMomentum];
    for (int j = 0; j < this->MaxYMomentum; ++j)
    {
      FourrierCoefficients[i][j] = Phase (-2.0 * M_PI * ((double) (i * this->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * this->YMomentum) / ((double) this->MaxYMomentum)));
    }
  }
  for (int i = 0; i < TmpSpace->HilbertSpaceDimension; ++i)
    {
      int nbrTranslationX;
      int nbrTranslationY;
      unsigned long TmpState = TmpSpace->StateDescription[i];
      TmpState = this->FindCanonicalForm(TmpState, nbrTranslationX, nbrTranslationY);
      nbrTranslationX = (this->MaxXMomentum - nbrTranslationX) % this->MaxXMomentum;
      nbrTranslationY = (this->MaxYMomentum - nbrTranslationY) % this->MaxYMomentum;
      int TmpMaxMomentum = 2 * this->NbrSite + 1;
      while ((TmpState >> TmpMaxMomentum) == 0x0ul)
	--TmpMaxMomentum;
      
      int Pos = this->FindStateIndex(TmpState, TmpMaxMomentum);
      if (Pos < this->HilbertSpaceDimension)
	{
	  TmpVector[i] =  (state[Pos] * (1.0 - (2.0 * ((double) ((this->ReorderingSign[Pos] >> ((nbrTranslationY * this->MaxXMomentum) + nbrTranslationX)) & 0x1ul)))) * 
			       FourrierCoefficients[nbrTranslationX][nbrTranslationY] / sqrt((double) this->NbrStateInOrbit[Pos]));
	}
      }
  delete[] FourrierCoefficients;
  return TmpVector;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = subsystem momentum along the x direction
// kySector = subsystem momentum along the x direction
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, 
														  ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      if ((kxSector == this->XMomentum) && (kySector == this->YMomentum))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementaryKxMomentum = (this->XMomentum - kxSector);
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->MaxXMomentum;
  int ComplementaryKyMomentum = (this->YMomentum - kySector);
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->MaxYMomentum;
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation SubsytemSpace (nbrParticleSector, this->NbrSite, kxSector, this->MaxXMomentum, kySector, this->MaxYMomentum);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;
  FQHETorusParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, (ParticleOnTorusWithSpinAndMagneticTranslations*) &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  cout << "nbr matrix elements non zero = " << Operation.GetNbrNonZeroMatrixElements() << endl;
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}
  
  
// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in given momentum and Sz sectors.
//
// nbrParticleSector = number of particles that belong to the subsytem
// kxSector = subsystem momentum along the x direction
// kySector = subsystem momentum along the x direction
// szSector  = twice the total Sz value of the subsytem
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int szSector, int kxSector, int kySector, 
																	 ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0) && (szSector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      if ((kxSector == this->XMomentum) && (kySector == this->YMomentum) && (this->TotalSpin == szSector))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  if (abs(ComplementarySzSector) > ComplementaryNbrParticles)
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

  int ComplementaryKxMomentum = (this->XMomentum - kxSector);
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->MaxXMomentum;
  int ComplementaryKyMomentum = (this->YMomentum - kySector);
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->MaxYMomentum;
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation SubsytemSpace (nbrParticleSector, szSector, this->NbrSite, kxSector, this->MaxXMomentum, kySector, this->MaxYMomentum);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace ComplementarySpace (ComplementaryNbrParticles, ComplementarySzSector, this->NbrSite);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;
  FQHETorusParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, (ParticleOnTorusWithSpinAndMagneticTranslations*) &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  cout << "nbr matrix elements non zero = " << Operation.GetNbrNonZeroMatrixElements() << endl;
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}
  
// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, 
													   ParticleOnTorusWithSpinAndMagneticTranslations* complementaryHilbertSpace,  
													   ParticleOnTorusWithSpinAndMagneticTranslations* destinationHilbertSpace,
													   ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation* TmpDestinationHilbertSpace =  (FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation*) destinationHilbertSpace;
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace* TmpHilbertSpace = (FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) complementaryHilbertSpace;
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace* TmpDestinationFullHilbertSpace = 0;
  if (TmpDestinationHilbertSpace->SzFlag == false)
    {
      TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace(TmpDestinationHilbertSpace->NbrFermions,
												    TmpDestinationHilbertSpace->NbrSite);
    }
  else
    {
      TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace(TmpDestinationHilbertSpace->NbrFermions,
												    TmpDestinationHilbertSpace->TotalSpin,
												    TmpDestinationHilbertSpace->NbrSite);
    }
  int* TmpStatePosition = new int [TmpHilbertSpace->GetLargeHilbertSpaceDimension()];
  int* TmpStatePosition2 = new int [TmpHilbertSpace->GetLargeHilbertSpaceDimension()];
  Complex* TmpStateCoefficient = new Complex [TmpHilbertSpace->GetLargeHilbertSpaceDimension()];
  Complex** FourrierCoefficientsDestination = new Complex* [this->MomentumModulo];
  Complex** FourrierCoefficients = new Complex* [this->MomentumModulo];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      FourrierCoefficientsDestination[i] = new Complex [this->MaxYMomentum];
      FourrierCoefficients[i] = new Complex [this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{
	  FourrierCoefficientsDestination[i][j] = Phase (2.0 * M_PI * ((double) (i * TmpDestinationHilbertSpace->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * TmpDestinationHilbertSpace->YMomentum) / ((double) this->MaxYMomentum)));
	  FourrierCoefficients[i][j] = Phase (2.0 * M_PI * ((double) (i * this->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * this->YMomentum) / ((double) this->MaxYMomentum)));
	}
    }
  
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  BinomialCoefficients TmpBinomial (this->NbrFermions);
  double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace->StateDescription[minIndex];
      for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationFullHilbertSpace->StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
	      int TmpDestinationNbrTranslationX;
	      int TmpDestinationNbrTranslationY;
	      unsigned long TmpCanonicalState2 = TmpDestinationHilbertSpace->FindCanonicalForm(TmpState2, TmpDestinationNbrTranslationX, TmpDestinationNbrTranslationY);
	      int TmpDestinationLzMax = 2 * TmpDestinationHilbertSpace->NbrSite - 1;
	      while ((TmpCanonicalState2 >> TmpDestinationLzMax) == 0x0ul)
		--TmpDestinationLzMax;
	      int RealDestinationIndex = TmpDestinationHilbertSpace->FindStateIndex(TmpCanonicalState2, TmpDestinationLzMax);
	      if (RealDestinationIndex < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
		{
		  int TmpLzMax = 2 * this->NbrSite - 1;
		  int TmpNbrTranslationX;
		  int TmpNbrTranslationY;
		  unsigned long TmpState3 = this->FindCanonicalForm((TmpState | TmpState2), TmpNbrTranslationX, TmpNbrTranslationY);
		  while ((TmpState3 >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
		  if (TmpPos < this->HilbertSpaceDimension)
		    {      
		      int NbrTranslationX = (this->MaxXMomentum - TmpNbrTranslationX) % this->MaxXMomentum;
		      int NbrTranslationY = (this->MaxYMomentum - TmpNbrTranslationY) % this->MaxYMomentum;
		      int DestinationNbrTranslationX = (this->MaxXMomentum - TmpDestinationNbrTranslationX) % this->MaxXMomentum;
		      int DestinationNbrTranslationY = (this->MaxYMomentum - TmpDestinationNbrTranslationY) % this->MaxYMomentum;
		      Complex Coefficient = TmpInvBinomial * FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY] * FourrierCoefficientsDestination[TmpDestinationNbrTranslationX][TmpDestinationNbrTranslationY] / sqrt ((double) (this->NbrStateInOrbit[TmpPos] * TmpDestinationHilbertSpace->NbrStateInOrbit[RealDestinationIndex]));
		      unsigned long Sign =  ((this->ReorderingSign[TmpPos] >> ((NbrTranslationY * this->MaxXMomentum) + NbrTranslationX))
					     ^ (TmpDestinationHilbertSpace->ReorderingSign[RealDestinationIndex] >> (DestinationNbrTranslationY * TmpDestinationHilbertSpace->MaxXMomentum + DestinationNbrTranslationX))) & 0x1ul;
		      int Pos2 = 2 * TmpDestinationHilbertSpace->NbrSite - 1;
		      unsigned long TmpState22 = TmpState2;
		      while ((Pos2 > 0) && (TmpState22 != 0x0ul))
			{
			  while (((TmpState22 >> Pos2) & 0x1ul) == 0x0ul)
			    --Pos2;
			  TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
			  TmpState3 ^= TmpState3 >> 32;
#endif	
			  TmpState3 ^= TmpState3 >> 16;
			  TmpState3 ^= TmpState3 >> 8;
			  TmpState3 ^= TmpState3 >> 4;
			  TmpState3 ^= TmpState3 >> 2;
			  TmpState3 ^= TmpState3 >> 1;
			  Sign ^= TmpState3;
			  TmpState22 &= ~(0x1ul << Pos2);
			  --Pos2;
			}
		      if ((Sign & 0x1ul) != 0x0ul)		  
			Coefficient *= -1.0;
		      TmpStatePosition[Pos] = TmpPos;
		      TmpStatePosition2[Pos] = RealDestinationIndex;
		      TmpStateCoefficient[Pos] = Coefficient;
		      ++Pos;
		    }
		}
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]] * TmpStateCoefficient[j]);
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      delete[] FourrierCoefficientsDestination[i];
      delete[] FourrierCoefficients[i];
    }
  delete[] FourrierCoefficientsDestination;
  delete[] FourrierCoefficients;
  delete TmpDestinationFullHilbertSpace;
  return TmpNbrNonZeroElements;
}

// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

ComplexVector FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::ForgeSU2FromU1(ComplexVector& upState, ParticleOnSphere* upStateSpace,
												       ComplexVector& downState, ParticleOnSphere* downStateSpace)
{
  return this->ForgeSU2FromU1(upState, *((FermionOnLatticeRealSpaceAnd2DTranslation*) upStateSpace),
			      downState, *((FermionOnLatticeRealSpaceAnd2DTranslation*) downStateSpace));
}

// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

ComplexVector FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation::ForgeSU2FromU1(ComplexVector& upState, FermionOnLatticeRealSpaceAnd2DTranslation& upStateSpace,
												       ComplexVector& downState, FermionOnLatticeRealSpaceAnd2DTranslation& downStateSpace)
{
  ComplexVector FinalState(this->HilbertSpaceDimension, true);
  Complex** FourrierCoefficients = new Complex* [this->MomentumModulo];
  Complex** FourrierCoefficientsDownState = new Complex* [this->MomentumModulo];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      FourrierCoefficientsDownState[i] = new Complex [this->MaxYMomentum];
      FourrierCoefficients[i] = new Complex [this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{
	  FourrierCoefficientsDownState[i][j] = Phase (2.0 * M_PI * ((double) (i * downStateSpace.XMomentum) / ((double) this->MaxXMomentum) + (double) (j * downStateSpace.YMomentum) / ((double) this->MaxYMomentum)));
	  FourrierCoefficients[i][j] = Phase (2.0 * M_PI * ((double) (i * this->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * this->YMomentum) / ((double) this->MaxYMomentum)));
	}
    }

   if ((upStateSpace.NbrFermions + downStateSpace.NbrFermions) != this->NbrSite)
     {
       cout << "cases where the number of sites is not equal to the total number of fermions is not supported" << endl;
     }
   else
     {
      unsigned long TmpSpinlessMask = (0x1ul << this->NbrSite) - 0x1ul;
      for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpUpState = upStateSpace.StateDescription[j];
	  unsigned long TmpDownState = (~TmpUpState) & TmpSpinlessMask;
	  int TmpDownStateNbrTranslationX;
	  int TmpDownStateNbrTranslationY;
	  unsigned long TmpCanonicalDownState = downStateSpace.FindCanonicalForm(TmpDownState, TmpDownStateNbrTranslationX, TmpDownStateNbrTranslationY);
	  int TmpPos = downStateSpace.NbrSite - 1;
	  while ((TmpPos > 0) && ((TmpCanonicalDownState & (0x1ul << TmpPos)) == 0x0ul))
	    {
	      --TmpPos;
	    }
	  int TmpSpinDownIndex = downStateSpace.FindStateIndex(TmpCanonicalDownState, TmpPos);
	  if (TmpSpinDownIndex != downStateSpace.HilbertSpaceDimension)
	    {
	      TmpPos = upStateSpace.NbrSite - 1;
	      while (TmpPos > 0)
		{
		  unsigned long Tmp = TmpUpState & (0x1ul << TmpPos);
		  TmpUpState |= Tmp << TmpPos;
		  TmpUpState ^= Tmp;
		  --TmpPos;
		}
	      TmpUpState <<= 1;
	      TmpPos = downStateSpace.NbrSite - 1;
	      while (TmpPos > 0)
		{
		  unsigned long Tmp = TmpDownState & (0x1ul << TmpPos);
		  TmpDownState |= Tmp << TmpPos;
		  TmpDownState ^= Tmp;
		  --TmpPos;
		}
	      TmpUpState |= TmpDownState;
	      int TmpNbrTranslationX;
	      int TmpNbrTranslationY;
	      unsigned long TmpSpinfulCanonical = this->FindCanonicalForm(TmpUpState, TmpNbrTranslationX, TmpNbrTranslationY);
#ifdef  __64_BITS__
	      int Max = 63;
#else
	      int Max = 31;
#endif
	      while ((TmpSpinfulCanonical & (0x1ul << Max)) == 0x0ul)
		{
		  --Max;
		}
	      int TmpSpinfulIndex = this->FindStateIndex(TmpSpinfulCanonical, Max);
	      if (TmpSpinfulIndex != this->HilbertSpaceDimension)
		{
		  unsigned long TmpUpState3 = TmpUpState;
		  unsigned long TmpUpState2 = TmpUpState3;
#ifdef  __64_BITS__
		  TmpUpState3 &= 0x5555555555555555ul;
		  TmpUpState2 &= 0xaaaaaaaaaaaaaaaaul;
#else
		  TmpUpState3 &= 0x55555555ul;
		  TmpUpState2 &= 0xaaaaaaaaul;
#endif	    
		  unsigned long Sign = 0x0;
		  int Pos = 2 * this->NbrSite - 2;
		  while ((Pos > 0) && ((TmpUpState3 & (0x1ul << Pos)) == 0x0ul))
		    Pos -= 2;
		  while (Pos > 0)
		    {
		      unsigned long TmpUpState4 = TmpUpState2 & ((0x1ul << Pos) - 1ul);
#ifdef  __64_BITS__
		      TmpUpState4 ^= TmpUpState4 >> 32;
#endif	
		      TmpUpState4 ^= TmpUpState4 >> 16;
		      TmpUpState4 ^= TmpUpState4 >> 8;
		      TmpUpState4 ^= TmpUpState4 >> 4;
		      TmpUpState4 ^= TmpUpState4 >> 2;
		      TmpUpState4 ^= TmpUpState4 >> 1;
		      Sign ^= TmpUpState4;
		      Pos -= 2;
		      while ((Pos > 0) && ((TmpUpState3 & (0x1ul << Pos)) == 0x0ul))
			Pos -= 2;
		    }
		  int NbrTranslationX = (this->MaxXMomentum - TmpNbrTranslationX) % this->MaxXMomentum;
		  int NbrTranslationY = (this->MaxYMomentum - TmpNbrTranslationY) % this->MaxYMomentum;
		  int DownStateNbrTranslationX = (this->MaxXMomentum - TmpDownStateNbrTranslationX) % this->MaxXMomentum;
		  int DownStateNbrTranslationY = (this->MaxYMomentum - TmpDownStateNbrTranslationY) % this->MaxYMomentum;
		  Sign ^= (downStateSpace.ReorderingSign[TmpSpinDownIndex] >> (TmpDownStateNbrTranslationY * downStateSpace.MaxXMomentum + TmpDownStateNbrTranslationX)) & 0x1ul;	      
		  Sign ^= (this->ReorderingSign[TmpSpinfulIndex] >> (NbrTranslationY * this->MaxXMomentum + NbrTranslationX)) & 0x1ul;
		  Complex Coefficient = ((FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY] *
					 FourrierCoefficientsDownState[DownStateNbrTranslationX][DownStateNbrTranslationY]) *
					 sqrt (((double) this->NbrStateInOrbit[TmpSpinfulIndex]) / ((double) (upStateSpace.NbrStateInOrbit[j] * downStateSpace.NbrStateInOrbit[TmpSpinDownIndex]))));
		  if ((Sign & 0x1ul) == 0x0ul)
		    FinalState[TmpSpinfulIndex] = upState[j] * downState[TmpSpinDownIndex] * Coefficient;
		  else
		    FinalState[TmpSpinfulIndex] = -upState[j] * downState[TmpSpinDownIndex] * Coefficient;		  
		}
	    }
	}
     }
   
   for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      delete[] FourrierCoefficients[i];
      delete[] FourrierCoefficientsDownState[i];
    } 
  delete[] FourrierCoefficients;
  delete[] FourrierCoefficientsDownState;
  return FinalState;
}
