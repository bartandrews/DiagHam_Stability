////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Nicolas Regnault                      //
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
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"
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

FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinRealSpaceAnd2DTranslation ()
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
  this->MomentumModulo = 1;
  this->XMomentum = 0; 
  this->YMomentum = 0;
  this->NbrSitePerUnitCell = 0;
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
  this->RescalingFactors = 0;
  this->TargetSpace = this;
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

FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (int nbrFermions, int nbrSite, int xMomentum, int maxXMomentum,
												      int yMomentum, int  maxYMomentum, unsigned long memory)
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
  this->MaxYMomentum = maxYMomentum;
  this->NbrSitePerUnitCell = this->NbrSite /  (this->MaxYMomentum * this->MaxXMomentum);
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = this->MaxXMomentum;

  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = 2 * (this->NbrSite / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this->MaxMomentum) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;
  this->TargetSpace = this;

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
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
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
	}
#endif
    }
}



// basic constructor when Sz is preserved
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// totalSpin = twice the total spin value
// xMomentum = momentum sector in the x direction
// maxXMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (int nbrFermions, int totalSpin, int nbrSite, int xMomentum, int maxXMomentum,
												      int yMomentum, int  maxYMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (totalSpin + this->NbrFermions) >> 1;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
  this->NbrSitePerUnitCell = this->NbrSite /  (this->MaxYMomentum * this->MaxXMomentum);
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
  this->TargetSpace = this;

  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFermionsUp);
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
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
	}
#endif
    }
}



// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::FermionOnLatticeWithSpinRealSpaceAnd2DTranslation(const FermionOnLatticeWithSpinRealSpaceAnd2DTranslation& fermions)
{
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrSite = fermions.NbrSite;
  this->NbrFermionsParity = fermions.NbrFermionsParity;

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

FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::~FermionOnLatticeWithSpinRealSpaceAnd2DTranslation ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeWithSpinRealSpaceAnd2DTranslation& FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::operator = (const FermionOnLatticeWithSpinRealSpaceAnd2DTranslation& fermions)
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
  this->NbrFermionsParity = fermions.NbrFermionsParity;

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
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;

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

AbstractHilbertSpace* FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::Clone()
{
  return new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      Tmp = (TmpState >> (i << 1)) & 03ul;
      switch (Tmp)
	{
	case 0x0ul:
	  Str << "0 ";
	  break;
	case 0x1ul:
	  Str << "d ";
	  break;
	case 0x2ul:
	  Str << "u ";
	  break;
	case 0x3ul:
	  Str << "X ";
	  break;
	}
    }
//    Str << " " << this->FindStateIndex(this->StateDescription[state], this->StateHighestBit[state]);
//    if (this->FindStateIndex(this->StateDescription[state], this->StateHighestBit[state]) != state)
//      {
//        Str << "  error";
//      }
//    Str << " orb=" << this->NbrStateInOrbit[state];
  return Str;
}


// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
  this->TargetSpace = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) targetSpace;
}

// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  if (this->SzFlag == false)
    this->RawGenerateStates(this->NbrFermions, this->NbrSite - 1, 0l);
  else
    this->RawGenerateStates(this->NbrFermions, this->NbrSite - 1, this->NbrFermionsUp, 0l);
  return this->CoreGenerateStates();
}


// generate all states corresponding to the constraints (core part of the method)
//
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::CoreGenerateStates()
{
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY) == this->StateDescription[i]))
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      ++TmpLargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescription[i] = 0x0ul;
	    }
	}
      else
	{
	  this->StateDescription[i] = 0x0ul;
	}
    }
  if (TmpLargeHilbertSpaceDimension == 0l)
    return 0l;
  unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
  this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
  this->ReorderingSign = new unsigned long [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->StateDescription[i] != 0x0ul)
	{
	  TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	  this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = this->FindOrbitSize(this->StateDescription[i]);
	  unsigned long& TmpSign = this->ReorderingSign[TmpLargeHilbertSpaceDimension];
	  TmpSign = 0x0ul;	  
	  int Index = 1;
	  unsigned long TmpState =  this->StateDescription[i];
	  for (int m = 0; m < this->MaxYMomentum; ++m)
	    {
	      unsigned long TmpState2 = TmpState;
	      for (int n = 1; n < this->MaxXMomentum; ++n)
		{
		  TmpSign |= (this->GetSignAndApplySingleXTranslation(TmpState2) << Index) ^ ((TmpSign & (0x1ul << (Index - 1))) << 1);
		  ++Index;
		}
	      TmpSign |= ((this->GetSignAndApplySingleYTranslation(TmpState) << Index) 
			  ^ ((TmpSign & (0x1ul << (Index - this->MaxXMomentum))) << this->MaxXMomentum));
	      ++Index;
	    }
	  ++TmpLargeHilbertSpaceDimension;
	}
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;

  this->StateHighestBit = new int [TmpLargeHilbertSpaceDimension];  
  int CurrentMaxMomentum = (2 * this->MaxMomentum) + 1;
  while (((this->StateDescription[0] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
    --CurrentMaxMomentum;
  this->StateHighestBit[0] = CurrentMaxMomentum;
  for (long i = 1l; i < TmpLargeHilbertSpaceDimension; ++i)
    {
      while (((this->StateDescription[i] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
	--CurrentMaxMomentum;
      this->StateHighestBit[i] = CurrentMaxMomentum;
    }
  return TmpLargeHilbertSpaceDimension;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::RawGenerateStates(int nbrFermions, int currentSite, long pos)
{
  if (nbrFermions == 0)
    {
      this->StateDescription[pos] = 0x0ul;	  
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0))
    return pos;
  if (nbrFermions == 1)
    {
      for (int j = currentSite; j >= 0; --j)
	{
	  this->StateDescription[pos] = 0x2ul << (j << 1);
	  ++pos;
	  this->StateDescription[pos] = 0x1ul << (j << 1);
	  ++pos;
	}
      return pos;
    }
  long TmpPos = this->RawGenerateStates(nbrFermions - 2, currentSite - 1, pos);
  unsigned long Mask = 0x3ul << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, pos);
  Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
   TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, pos);
   Mask = 0x1ul << (currentSite << 1);
   for (; pos < TmpPos; ++pos)
     this->StateDescription[pos] |= Mask;
  return this->RawGenerateStates(nbrFermions, currentSite - 1, pos);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrSpinUp = number of fermions with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::RawGenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos)
{
  if ((nbrFermions == 0) && (nbrSpinUp == 0))
    {
      this->StateDescription[pos] = 0x0ul;	  
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0) || (nbrSpinUp > nbrFermions) || (nbrSpinUp < 0))
    return pos;
  if (nbrFermions == 1)
    {
      if (nbrSpinUp == 1)
	{
	  for (int j = currentSite; j >= 0; --j)
	    {
	      this->StateDescription[pos] = 0x2ul << (j << 1);
	      ++pos;
	    }
	}
      else
	{
	  for (int j = currentSite; j >= 0; --j)
	    {
	      this->StateDescription[pos] = 0x1ul << (j << 1);
	      ++pos;
	    }
	}
      return pos;
    }
  long TmpPos = this->RawGenerateStates(nbrFermions - 2, currentSite - 1, nbrSpinUp - 1, pos);
  unsigned long Mask = 0x3ul << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, nbrSpinUp - 1, pos);
  Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, nbrSpinUp, pos);
  Mask = 0x1ul << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->RawGenerateStates(nbrFermions, currentSite - 1, nbrSpinUp, pos);   
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::EvaluateHilbertSpaceDimension(int nbrFermions)
{
  BinomialCoefficients binomials(2*this->NbrSite);
  long dimension = binomials(2*this->NbrSite, this->NbrFermions);
  return dimension;
}

// evaluate Hilbert space dimension with a fixed number of fermions with spin up
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrSpinUp = number of fermions with spin up
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::EvaluateHilbertSpaceDimension(int nbrFermions,int nbrSpinUp)
{
  BinomialCoefficients binomials(this->NbrSite);
  long Dimension = binomials(this->NbrSite, this->NbrFermionsUp);
  Dimension *= binomials(this->NbrSite, this->NbrFermionsDown);
  return Dimension;
}


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::GenerateLookUpTable(unsigned long memory)
{
  // get every highest bit poisition
  unsigned long TmpPosition = this->StateDescription[0];
#ifdef __64_BITS__
  int CurrentHighestBit = 63;
#else
  int CurrentHighestBit = 31;
#endif
  while ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
    --CurrentHighestBit;  

  // evaluate look-up table size
  memory /= (sizeof(int*) * this->NbrFermionStates);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrFermionStates)
    this->MaximumLookUpShift = this->NbrFermionStates;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrFermionStates];
  this->LookUpTableShift = new int [this->NbrFermionStates];
  for (int i = 0; i < this->NbrFermionStates; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentLargestBit = CurrentHighestBit;
  int* TmpLookUpTable = this->LookUpTable[CurrentLargestBit];
  if (CurrentLargestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentLargestBit] = 0;
  else
    this->LookUpTableShift[CurrentLargestBit] = CurrentLargestBit + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentLargestBit];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {     
      TmpPosition = this->StateDescription[i];
      while ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
	--CurrentHighestBit;  
      if (CurrentLargestBit != CurrentHighestBit)
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  --CurrentLargestBit;
	  while (CurrentLargestBit > CurrentHighestBit)
	    {
	      this->LookUpTableShift[CurrentLargestBit] = -1;
	      --CurrentLargestBit;
	    }
 	  CurrentLargestBit = CurrentHighestBit;
	  TmpLookUpTable = this->LookUpTable[CurrentLargestBit];
	  if (CurrentLargestBit < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentLargestBit] = 0;
	  else
	    this->LookUpTableShift[CurrentLargestBit] = CurrentLargestBit + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentLargestBit];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;

  this->RescalingFactors = new double* [this->NbrMomentum];
  for (int i = 1; i <= this->MaxMomentum; ++i)
    {
      this->RescalingFactors[i] = new double [this->NbrMomentum];
      for (int j = 1; j <= this->MaxMomentum; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}

// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::FindStateIndex(char* stateDescription)
{
  char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != this->MaxMomentum)
    return -1;
  unsigned long TmpState = 0x0ul;
  int TmpNbrParticles = 0;
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      if (TmpDescription[i][0] == 'u')
	{
	  TmpState |= 0x2ul << (2 * i);
	  ++TmpNbrParticles;	  
	}
      else
	{
	  if (TmpDescription[i][0] == 'd')
	    {
	      TmpState |= 0x1ul << (2 * i);
	      ++TmpNbrParticles;	  
	    }
	  else
	    {
	      if (TmpDescription[i][0] == 'X')
		{
		  TmpState |= 0x3ul << (2 * i);
		  TmpNbrParticles += 2;	  
		}
	      else
		{
		  if (TmpDescription[i][0] != '0')
		    {
		      return -1;
		    }
		}
	    }
	}
      delete[] TmpDescription[i];
    }
  delete[] TmpDescription;
  if (TmpNbrParticles != this->NbrFermions)
    return -1;
  int TmpLzMax = 2 * this->MaxMomentum + 1;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
}

// find state index
//
// stateDescription = unsigned longeger describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::FindStateIndex(unsigned long stateDescription, int maxMomentum)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    return this->HilbertSpaceDimension;
  if (this->LookUpTableShift[maxMomentum] < 0)
    return this->HilbertSpaceDimension;
  long PosMax = stateDescription >> this->LookUpTableShift[maxMomentum];
  long PosMin = this->LookUpTable[maxMomentum][PosMax];
  PosMax = this->LookUpTable[maxMomentum][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }

  if (CurrentState == stateDescription)
    return PosMid;
  else
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}

// apply a^+_m_sigma a_n_sigma operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator including the orbital and the spin index
// n = index of the annihilation operator including the orbital and the spin index
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::AdsigmaAsigma (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long State = this->StateDescription[index];
  if ((State & (0x1ul << n)) == 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(0x1ul << n);
  if ((State & (0x1ul << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
  coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
  coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
  State |= (0x1ul << m);
  this->ProdATemporaryNbrStateInOrbit =  this->NbrStateInOrbit[index];
  return this->SymmetrizeAdAdResult(State, coefficient, nbrTranslationX, nbrTranslationY);
}
  
// apply a_n1_sigma a_n2_sigma operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdsigmaAdsigma call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 

double FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::AsigmaAsigma (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || 
      ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;  
  double coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  ProdATemporaryState &= ~(0x1ul << n2);
  coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  ProdATemporaryState &= ~(0x1ul << n1);
  this->ProdATemporaryNbrStateInOrbit =  this->NbrStateInOrbit[index];
  return coefficient;
}


// apply a^+_m1_sigma a^+_m2_sigma operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations to applied in the y direction to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::AdsigmaAdsigma (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
  TmpState |= (0x1ul << m2);
  
  coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
  TmpState |= (0x1ul << m1);
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}

// apply a^+_m1_sigma a^+_m2_sigma operator to the state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::AdsigmaAdsigma (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long TmpState = this->StateDescription[index];
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
  TmpState |= (0x1ul << m2);
  
  coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
  TmpState |= (0x1ul << m1);
  this->ProdATemporaryNbrStateInOrbit =  this->NbrStateInOrbit[index];
  return this->SymmetrizeAdAdResultTarget(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}

// apply a^+_m operator to the state produced using Au method (without destroying it)
//
// m = first index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::Adu (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m <<= 1;
  ++m;
  
  if ((TmpState & (0x1ul << m)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  
  coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
  TmpState |= (0x1ul << m);
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}

// apply a^+_m operator to the state produced using Au method (without destroying it)
//
// m = first index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::Add (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m <<= 1;
    
  if ((TmpState & (0x1ul << m)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  
  coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
  TmpState |= (0x1ul << m);
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}

// convert a state defined in the real space basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::ConvertToKxKyBasis(ComplexVector& state, ParticleOnSphere* space)  
{
  FermionOnLatticeWithSpinRealSpace* TmpSpace = (FermionOnLatticeWithSpinRealSpace*) space;
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

ComplexVector FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space)
{
  FermionOnLatticeWithSpinRealSpace* TmpSpace = (FermionOnLatticeWithSpinRealSpace*) space;
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

// convert a given state from a given  n-body basis basis to another one
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis where state is defined
// return value = converted vector

ComplexVector FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::ConvertToNbodyBasis(ComplexVector& state, FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* nbodyBasis)
{
  ComplexVector TmpVector (this->LargeHilbertSpaceDimension, true);
  if (nbodyBasis->LargeHilbertSpaceDimension < this->LargeHilbertSpaceDimension)
    {
      for (long i= 0l; i < nbodyBasis->LargeHilbertSpaceDimension; ++i)
	{
	  int TmpLzMax = 2 * nbodyBasis->NbrSite - 1;
	  unsigned long TmpState = nbodyBasis->StateDescription[i];
	  while ((TmpState >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  int TmpIndex = this->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpIndex < this->HilbertSpaceDimension)
	    {
	      TmpVector[TmpIndex] = state[i];
	    }
	}
    }
  else
    {
      for (long i= 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  int TmpLzMax = 2 * this->NbrSite - 1;
	  unsigned long TmpState = this->StateDescription[i];
	  while ((TmpState >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  int TmpIndex = nbodyBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpIndex < nbodyBasis->HilbertSpaceDimension)
	    {
	      TmpVector[i] = state[TmpIndex];
	    }
	}
    }
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

HermitianMatrix FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, 
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
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation SubsytemSpace (nbrParticleSector, this->NbrSite, kxSector, this->MaxXMomentum, kySector, this->MaxYMomentum);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;
  architecture->SetDimension(ComplementarySpace.GetHilbertSpaceDimension());
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

HermitianMatrix FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int szSector, int kxSector, int kySector, 
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
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation SubsytemSpace (nbrParticleSector, szSector, this->NbrSite, kxSector, this->MaxXMomentum, kySector, this->MaxYMomentum);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation ComplementarySpace (ComplementaryNbrParticles, ComplementarySzSector, this->NbrSite, 
  									ComplementaryKxMomentum, this->MaxXMomentum, ComplementaryKyMomentum, this->MaxYMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;
  architecture->SetDimension(ComplementarySpace.GetHilbertSpaceDimension());
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

// evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = subsystem momentum along the x direction
// kySector = subsystem momentum along the x direction
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, 
														  int nbrGroundStates, ComplexVector* groundStates, double* weights, 
														  AbstractArchitecture* architecture)
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
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation SubsytemSpace (nbrParticleSector, this->NbrSite, kxSector, this->MaxXMomentum, kySector, this->MaxYMomentum);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;
  architecture->SetDimension(ComplementarySpace.GetHilbertSpaceDimension());
  FQHETorusParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, (ParticleOnTorusWithSpinAndMagneticTranslations*) &ComplementarySpace, 
							   nbrGroundStates,  groundStates, weights, TmpDensityMatrix);
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
  
// evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. The density matrix is only evaluated in given momentum and Sz sectors.
//
// nbrParticleSector = number of particles that belong to the subsytem
// kxSector = subsystem momentum along the x direction
// kySector = subsystem momentum along the x direction
// szSector  = twice the total Sz value of the subsytem
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// architecture = pointer to the architecture to use parallelized algorithm
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int szSector, int kxSector, int kySector, 
														  int nbrGroundStates, ComplexVector* groundStates, double* weights, 
														  AbstractArchitecture* architecture)
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
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation SubsytemSpace (nbrParticleSector, szSector, this->NbrSite, kxSector, this->MaxXMomentum, kySector, this->MaxYMomentum);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation ComplementarySpace (ComplementaryNbrParticles, ComplementarySzSector, this->NbrSite, 
  									ComplementaryKxMomentum, this->MaxXMomentum, ComplementaryKyMomentum, this->MaxYMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;
  architecture->SetDimension(ComplementarySpace.GetHilbertSpaceDimension());
  FQHETorusParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, (ParticleOnTorusWithSpinAndMagneticTranslations*) &ComplementarySpace, 
							   nbrGroundStates,  groundStates, weights, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  cout << "nbr matrix elements non zero = " << Operation.GetNbrNonZeroMatrixElements() << endl;
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    {
      if (Operation.GetMatrix().GetNbrRow() != TmpDensityMatrix.GetNbrRow())
	{
	  TmpDensityMatrix = Operation.GetMatrix();
	}
      return TmpDensityMatrix;
    }
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

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, 
													   ParticleOnTorusWithSpinAndMagneticTranslations* complementaryHilbertSpace,  
													   ParticleOnTorusWithSpinAndMagneticTranslations* destinationHilbertSpace,
													   ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* TmpDestinationHilbertSpace =  (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) destinationHilbertSpace;
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* TmpHilbertSpace = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) complementaryHilbertSpace;
  FermionOnLatticeWithSpinRealSpace* TmpDestinationFullHilbertSpace = 0;
  if (TmpDestinationHilbertSpace->SzFlag == false)
    {
      TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpace(TmpDestinationHilbertSpace->NbrFermions,
									     TmpDestinationHilbertSpace->NbrSite);
    }
  else
    {
      TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpace(TmpDestinationHilbertSpace->NbrFermions,
									     TmpDestinationHilbertSpace->TotalSpin,
									     TmpDestinationHilbertSpace->NbrSite);
    }
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
  ComplexMatrix TmpEntanglementMatrix (nbrIndex, TmpDestinationHilbertSpace->GetHilbertSpaceDimension(), true);
  
  int TmpMinIndex = minIndex;
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
 		      TmpEntanglementMatrix[RealDestinationIndex][minIndex - TmpMinIndex] += groundState[TmpPos] * Coefficient;
		    }
		}
	    }
	}
    }
  for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      for (int j = i; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  ++TmpNbrNonZeroElements;
	  densityMatrix->UnsafeAddToMatrixElement(i, j, TmpEntanglementMatrix[j] * TmpEntanglementMatrix[i]);
	}
    }
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

// core part of the evaluation density matrix particle partition calculation involving a sum of projectors
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, 
													   ParticleOnTorusWithSpinAndMagneticTranslations* complementaryHilbertSpace,  
													   ParticleOnTorusWithSpinAndMagneticTranslations* destinationHilbertSpace,
													   int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix* densityMatrix)
{
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* TmpDestinationHilbertSpace =  (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) destinationHilbertSpace;
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* TmpHilbertSpace = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) complementaryHilbertSpace;
  FermionOnLatticeWithSpinRealSpace* TmpDestinationFullHilbertSpace = 0;
  if (TmpDestinationHilbertSpace->SzFlag == false)
    {
      TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpace(TmpDestinationHilbertSpace->NbrFermions,
									     TmpDestinationHilbertSpace->NbrSite);
    }
  else
    {
      TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpace(TmpDestinationHilbertSpace->NbrFermions,
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
  Complex* TmpValues = new Complex[nbrGroundStates];

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
	      for (int l = 0; l < nbrGroundStates; ++l)
		TmpValues[l] = weights[l] * Conj(groundStates[l][TmpStatePosition[j]] * TmpStateCoefficient[j]);
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    Complex Tmp = 0.0;
		    for (int l = 0; l < nbrGroundStates; ++l)
		      Tmp += TmpValues[l] * groundStates[l][TmpStatePosition[k]] * TmpStateCoefficient[k];
		    densityMatrix->UnsafeAddToMatrixElement(Pos2, TmpStatePosition2[k], Tmp);		
		  }
	    }
	}
    }
  delete[] TmpValues;
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

// generate an eta pairing state
// 
// return value = vector corresponding to the  eta pairing state

ComplexVector FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::GenerateEtaPairingState()
{
  ComplexVector TmpVector (this->LargeHilbertSpaceDimension, true);
  double* TmpFactors = new double[this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
    {
      int TmpX;
      int TmpY;
      int TmpOrbitalIndex;
      this->GetLinearizedIndex(i, TmpX, TmpY, TmpOrbitalIndex);
      TmpFactors[i] = (double) (1 - (((TmpX ^ TmpY) & 1) << 1));
    }
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      double TmpFactor = 1.0;
      unsigned long TmpState = this->StateDescription[i];
      for (int j = 0; (j < this->NbrSite) && (TmpFactor != 0.0); ++j)
	{	  
	  if ((TmpState & 0x3ul) == 0x3ul)
	    {
	      TmpFactor *= TmpFactors[j];
	    }
	  else
	    {
	      if ((TmpState & 0x3ul) != 0x0ul)
		TmpFactor = 0.0;
	    }
	  TmpState >>= 2;
	}
      if (TmpFactor != 0.0)
	TmpVector[i] = TmpFactor * sqrt((double) this->NbrStateInOrbit[i]);
    }
  delete[] TmpFactors;
  TmpVector /= TmpVector.Norm();
  return TmpVector;
}

// generate a state near the eta pairing state (with out broken pair)
// 
// xDistance = distance along x between the two particles of the broken pair
// yDistance = distance along y between the two particles of the broken pair
// return value = vector corresponding to the  eta pairing state

ComplexVector FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::GenerateEtaPairingNearbyState(int xDistance, int yDistance)
{
//   if ((2 * xDistance) > this->MaxXMomentum)
//     xDistance = this->MaxXMomentum - xDistance;
//   if ((2 * yDistance) > this->MaxYMomentum)
//     yDistance = this->MaxYMomentum - yDistance;
  if ((xDistance == 0) && (yDistance == 0))
    {
      return this->GenerateEtaPairingState();
    }
  ComplexVector TmpVector (this->LargeHilbertSpaceDimension, true);
  double* TmpFactors = new double[this->NbrSite];
  int TmpX1;
  int TmpY1;
  int TmpOrbitalIndex1; 
  for (int i = 0; i < this->NbrSite; ++i)
    {
      this->GetLinearizedIndex(i, TmpX1, TmpY1, TmpOrbitalIndex1);
      TmpFactors[i] = (double) (1 - (((TmpX1 ^ TmpY1) & 1) << 1));
    }
  int TmpX2;
  int TmpY2;
  int TmpOrbitalIndex2;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      double TmpFactor = 1.0;
      unsigned long TmpState = this->StateDescription[i];
      int TwiceNbrBrokenPairs = 0;
      int BrokenPairCoordinates[3];
      for (int j = 0; (j < this->NbrSite) && (TwiceNbrBrokenPairs <= 2); ++j)
	{	  
	  if ((TmpState & 0x3ul) == 0x3ul)
	    {
	      TmpFactor *= TmpFactors[j];
	    }
	  else
	    {
	      if ((TmpState & 0x3ul) != 0x0ul)
		{
		  BrokenPairCoordinates[TwiceNbrBrokenPairs] = j;
		  ++TwiceNbrBrokenPairs;
		}
	    }
	  TmpState >>= 2;
	}
      if (TwiceNbrBrokenPairs == 2)
	{
	  this->GetLinearizedIndex(BrokenPairCoordinates[0], TmpX1, TmpY1, TmpOrbitalIndex1);
	  this->GetLinearizedIndex(BrokenPairCoordinates[1], TmpX2, TmpY2, TmpOrbitalIndex2);
	  if (((this->StateDescription[i] >> (BrokenPairCoordinates[0] << 1)) & 0x3ul) == 0x1ul)
	    {
	      TmpX1 += xDistance;
	      TmpY1 += yDistance;
	      TmpX1 %= this->MaxXMomentum;
	      TmpY1 %= this->MaxYMomentum;
	      if ((TmpX1 == TmpX2) && (TmpY1 == TmpY2))
		{
		  TmpVector[i] = TmpFactor * TmpFactors[BrokenPairCoordinates[0]] * sqrt((double) this->NbrStateInOrbit[i]);
		}
	    }
	  else
	    {
	      TmpX2 += xDistance;
	      TmpY2 += yDistance;
	      TmpX2 %= this->MaxXMomentum;
	      TmpY2 %= this->MaxYMomentum;
	      if ((TmpX1 == TmpX2) && (TmpY1 == TmpY2))
		{
		  TmpVector[i] = -TmpFactor * TmpFactors[BrokenPairCoordinates[1]] * sqrt((double) this->NbrStateInOrbit[i]);
		}
	    }
// 	  TmpX1 = abs(TmpX1 - TmpX2);
// 	  TmpY1 = abs(TmpY1 - TmpY2);
// 	  if ((2 * TmpX1) > this->MaxXMomentum)
// 	    TmpX1 = this->MaxXMomentum - TmpX1;
// 	  if ((2 * TmpY1) > this->MaxYMomentum)
// 	    TmpY1 = this->MaxYMomentum - TmpY1;
// 	  if ((TmpX1 == xDistance) && (TmpY1 == yDistance))
// 	    {
// 	      if (((this->StateDescription[i] >> (BrokenPairCoordinates[0] << 1)) & 0x3ul) == 0x1ul)
// 		{
// 		  TmpVector[i] = TmpFactor * TmpFactors[BrokenPairCoordinates[0]] * sqrt((double) this->NbrStateInOrbit[i]);
// 		}
// 	      else
// 		{
// 		  TmpVector[i] = TmpFactor * TmpFactors[BrokenPairCoordinates[1]] * sqrt((double) this->NbrStateInOrbit[i]);
// 		}
// 	    }	  
	}
    }
  delete[] TmpFactors;
  TmpVector /= TmpVector.Norm();
  return TmpVector;
}

// generate an eta pairing state, using a generic state instead of the vacuum
// 
// initialSpace = space where the generic state is defined
// initialState = generic state
// return value = vector corresponding to the eta pairing state built on top on the generic state

ComplexVector FermionOnLatticeWithSpinRealSpaceAnd2DTranslation::GenerateEtaPairingState(FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* initialSpace, ComplexVector& initialState)
{
  int NbrPairs = (this->NbrFermions - initialSpace->NbrFermions) >> 1;
  FermionOnLatticeRealSpace PairSpace (NbrPairs, this->NbrSite);

  ComplexVector TmpVector (this->LargeHilbertSpaceDimension, true);
  double* TmpFactors = new double[this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
    {
      int TmpX;
      int TmpY;
      int TmpOrbitalIndex;
      this->GetLinearizedIndex(i, TmpX, TmpY, TmpOrbitalIndex);
      TmpFactors[i] = (double) (1 - (((TmpX ^ TmpY) & 1) << 1));
    }
  Complex** FourrierCoefficients = new Complex* [this->MomentumModulo];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      FourrierCoefficients[i] = new Complex [this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{
	  FourrierCoefficients[i][j] = Phase (-2.0 * M_PI * ((double) (i * this->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * this->YMomentum) / ((double) this->MaxYMomentum)));
	}
    }
  unsigned long* ExpendedStateDescription = new unsigned long[PairSpace.LargeHilbertSpaceDimension];
  double* ExpendedFactors = new double[PairSpace.LargeHilbertSpaceDimension];
  for (long j = 0; j < PairSpace.LargeHilbertSpaceDimension; ++j)
    {
      unsigned long TmpState = 0x0ul;
      double TmpFactor = 1.0;
      unsigned long TmpState2 = PairSpace.StateDescription[j];
      for (int i = 0; i < this->NbrSite; ++i)
	{
	  if (((TmpState2 >> i) & 0x1ul) != 0x0ul)
	    {
	      TmpFactor *= TmpFactors[i];
	      TmpState |= 0x3ul << (i << 1);
	    }
	}
      ExpendedStateDescription[j] = TmpState;
      ExpendedFactors[j] = TmpFactor;
    }

  for (long i = 0l; i < initialSpace->LargeHilbertSpaceDimension; ++i)
    {
      double TmpFactor = 1.0;
      unsigned long TmpState = initialSpace->StateDescription[i];
      for (long j = 0; j < PairSpace.LargeHilbertSpaceDimension; ++j)
	{
	  if ((ExpendedStateDescription[j] & TmpState) == 0x0ul)
	    {
	      int TmpNbrTranslationX;
	      int TmpNbrTranslationY;
	      unsigned long TmpCanonicalState = this->FindCanonicalForm(ExpendedStateDescription[j] | TmpState, TmpNbrTranslationX, TmpNbrTranslationY);
	      int TmpLzMax = 2 * this->NbrSite - 1;
	      while ((TmpCanonicalState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpIndex = this->FindStateIndex(TmpCanonicalState, TmpLzMax);
	      if (TmpIndex < this->LargeHilbertSpaceDimension)
		{
  		  TmpVector[TmpIndex] = ((ExpendedFactors[j] * sqrt(((double) this->NbrStateInOrbit[TmpIndex]) / ((double) initialSpace->NbrStateInOrbit[i]))) 
  					 * FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY] * initialState[i]);
//  		  TmpVector[TmpIndex] = ((ExpendedFactors[j] * sqrt(((double) initialSpace->NbrStateInOrbit[i]) / ((double) this->NbrStateInOrbit[TmpIndex])))
//  					 * FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY] * initialState[i]);
		  TmpNbrTranslationX = (this->MaxXMomentum - TmpNbrTranslationX) % this->MaxXMomentum;
		  TmpNbrTranslationY = (this->MaxYMomentum - TmpNbrTranslationY) % this->MaxYMomentum;
		  TmpVector[TmpIndex] *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> ((TmpNbrTranslationY * this->MaxXMomentum) + TmpNbrTranslationX)) & 0x1ul))); 
		}
	    }
	}
    }

  delete[] ExpendedStateDescription;
  delete[] ExpendedFactors;
  delete[] FourrierCoefficients;
  delete[] TmpFactors;
  TmpVector /= TmpVector.Norm();
  return TmpVector;  
}

