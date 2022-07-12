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
//                        class author: Cecile Repellin                       //
//                                                                            //
//                        last modification : 02/04/2018                      //
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
#include "HilbertSpace/FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation.h"
#include "HilbertSpace/FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusion.h"
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

FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation::FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation ()
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

FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation::FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation (int nbrFermions, int totalSpin, int nbrSite, int xMomentum, int maxXMomentum, int yMomentum, int  maxYMomentum, unsigned long memory)
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

  
  this->InitializeHexagonArrays();
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFermionsUp, this->NbrSite - 1, 0x0ul);
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

FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation::FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation(const FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation& fermions)
{
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrSite = fermions.NbrSite;
  this->NbrFermionsParity = fermions.NbrFermionsParity;
  this->ListIndicesPerPlaquette = new int*[(this->NbrSite >> 1)];
  for (int i = 0; i < (this->NbrSite >> 1); ++i)
  {
    this->ListIndicesPerPlaquette[i] = new int[6];
    for (int j = 0; j < 6; ++j)
      this->ListIndicesPerPlaquette[i][j] = fermions.ListIndicesPerPlaquette[i][j];   
  }
  
  this->NbrLargerIndicesInPlaquette = new int*[this->NbrSite];
  this->LargerIndicesInPlaquette = new int**[this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
  {
    this->NbrLargerIndicesInPlaquette[i] = new int[3];
    this->LargerIndicesInPlaquette[i] = new int*[3];
    for (int j = 0; j < 3; ++j)
    {
      this->NbrLargerIndicesInPlaquette[i][j] = fermions.NbrLargerIndicesInPlaquette[i][j];
      if (this->NbrLargerIndicesInPlaquette[i][j] != 0)
	this->LargerIndicesInPlaquette[i][j] = new int[this->NbrLargerIndicesInPlaquette[i][j]];
      else
	this->LargerIndicesInPlaquette[i][j] = 0;
      for (int k = 0; k < this->NbrLargerIndicesInPlaquette[i][j]; ++k)
	this->LargerIndicesInPlaquette[i][j][k] = fermions.LargerIndicesInPlaquette[i][j][k];
    }
  }

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

FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation::~FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation ()
{
  for (int i = 0; i < (this->NbrSite >> 1); ++i)
    delete[] this->ListIndicesPerPlaquette[i];
  delete[] this->ListIndicesPerPlaquette;
  
  for (int i = 0; i < this->NbrSite; ++i)
  {
    delete[] this->NbrLargerIndicesInPlaquette[i];
    for (int j = 0; j < 3; ++j)
      delete[] this->LargerIndicesInPlaquette[i][j];
    delete[] this->LargerIndicesInPlaquette[i];
  }
  delete[] this->NbrLargerIndicesInPlaquette;
  delete[] this->LargerIndicesInPlaquette;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation& FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation::operator = (const FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation& fermions)
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
  this->ListIndicesPerPlaquette = new int*[(this->NbrSite >> 1)];
  for (int i = 0; i < (this->NbrSite >> 1); ++i)
  {
    this->ListIndicesPerPlaquette[i] = new int[6];
    for (int j = 0; j < 6; ++j)
      this->ListIndicesPerPlaquette[i][j] = fermions.ListIndicesPerPlaquette[i][j];   
  }
  
  this->NbrLargerIndicesInPlaquette = new int*[this->NbrSite];
  this->LargerIndicesInPlaquette = new int**[this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
  {
    this->NbrLargerIndicesInPlaquette[i] = new int[3];
    this->LargerIndicesInPlaquette[i] = new int*[3];
    for (int j = 0; j < 3; ++j)
    {
      this->NbrLargerIndicesInPlaquette[i][j] = fermions.NbrLargerIndicesInPlaquette[i][j];
      if (this->NbrLargerIndicesInPlaquette[i][j] != 0)
	this->LargerIndicesInPlaquette[i][j] = new int[this->NbrLargerIndicesInPlaquette[i][j]];
      else
	this->LargerIndicesInPlaquette[i][j] = 0;
      for (int k = 0; k < this->NbrLargerIndicesInPlaquette[i][j]; ++k)
	this->LargerIndicesInPlaquette[i][j][k] = fermions.LargerIndicesInPlaquette[i][j][k];
    }
  }

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

AbstractHilbertSpace* FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation::Clone()
{
  return new FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation(*this);
}


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(this->NbrFermions, this->NbrFermionsUp, this->NbrSite - 1, 0x0ul, 0l);
  return this->CoreGenerateStates();
}




// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrSpinUp = number of fermions with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation::RawGenerateStates(int nbrFermions, int nbrSpinUp, int currentSite, unsigned long currentConfiguration, long pos)
{
   if ((nbrFermions == 0) && (nbrSpinUp == 0))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0) || (nbrSpinUp > nbrFermions) || (nbrSpinUp < 0))
    {
      return pos;
    }

  long TmpPos;
  unsigned long Mask;
  int TmpCharge = this->FindMaximumChargeSurroundingPlaquettes(currentSite, currentConfiguration);
  if (TmpCharge > 6)
    return pos;
  if (TmpCharge < 5)
  {
    Mask =  0x3ul << (currentSite << 1);
    TmpPos = this->RawGenerateStates(nbrFermions - 2, nbrSpinUp - 1, currentSite - 1, currentConfiguration | Mask, pos);
    for (; pos < TmpPos; ++pos)
      this->StateDescription[pos] |= Mask;
  }
  if (TmpCharge < 6)
  {
    Mask =  0x2ul << (currentSite << 1);
    TmpPos = this->RawGenerateStates(nbrFermions - 1, nbrSpinUp - 1, currentSite - 1, currentConfiguration | Mask, pos);
    for (; pos < TmpPos; ++pos)
      this->StateDescription[pos] |= Mask;
    
    Mask =  0x1ul << (currentSite << 1);
    TmpPos = this->RawGenerateStates(nbrFermions - 1, nbrSpinUp, currentSite - 1, currentConfiguration | Mask, pos);
    for (; pos < TmpPos; ++pos)
      this->StateDescription[pos] |= Mask;
  }  
  return this->RawGenerateStates(nbrFermions, nbrSpinUp, currentSite - 1, currentConfiguration, pos);
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// nbrSpinUp = number of fermions with spin up
// currentSite = current site linearized index
// currentConfiguration = configuraton at the current x coordinate
// return value = Hilbert space dimension

long FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation::EvaluateHilbertSpaceDimension(int nbrFermions, int nbrSpinUp, int currentSite, unsigned long currentConfiguration)
{
  if ((nbrFermions == 0) && (nbrSpinUp == 0))
    {
      return 1l;
    }
    
  if ((currentSite < 0) || (nbrFermions < 0) || (nbrSpinUp > nbrFermions) || (nbrSpinUp < 0))
    {
      return 0l;
    }

  long TmpDimension = 0l;
  int TmpCharge = this->FindMaximumChargeSurroundingPlaquettes(currentSite, currentConfiguration);
  if (TmpCharge > 6)
    return TmpDimension;
  if (TmpCharge < 5)
      TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, nbrSpinUp - 1, currentSite - 1, currentConfiguration | (0x3ul << (currentSite << 1)));
  if (TmpCharge < 6)
  {
    TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, nbrSpinUp - 1, currentSite - 1, currentConfiguration | (0x2ul << (currentSite << 1)));
    TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, nbrSpinUp, currentSite - 1, currentConfiguration | (0x1ul << (currentSite << 1)));
  }  
  TmpDimension += this->EvaluateHilbertSpaceDimension(nbrFermions, nbrSpinUp, currentSite - 1, currentConfiguration);
  return TmpDimension;
}


// initialize all of the arrays that will be used to implement hexagon exclusion conditions
//

void FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation::InitializeHexagonArrays()
{
  this->ListIndicesPerPlaquette = new int*[this->NbrSite >> 1];
  for (int i = 0; i < this->MaxXMomentum; ++i)
  {
    for (int j = 0; j < this->MaxYMomentum; ++j)
    {
      int TmpIndex = this->MaxYMomentum * i + j;
      this->ListIndicesPerPlaquette[TmpIndex] = new int[6];      
      this->ListIndicesPerPlaquette[TmpIndex][0] = this->GetLinearizedIndexSafe(i, j, 0);
      this->ListIndicesPerPlaquette[TmpIndex][1] = this->GetLinearizedIndexSafe(i, j, 1);
      this->ListIndicesPerPlaquette[TmpIndex][2] = this->GetLinearizedIndexSafe(i + 1, j, 0);
      this->ListIndicesPerPlaquette[TmpIndex][3] = this->GetLinearizedIndexSafe(i + 1, j - 1, 1);
      this->ListIndicesPerPlaquette[TmpIndex][4] = this->GetLinearizedIndexSafe(i + 1, j - 1 , 0);
      this->ListIndicesPerPlaquette[TmpIndex][5] = this->GetLinearizedIndexSafe(i, j - 1, 1);
      
    }
  }
  
  this->NbrLargerIndicesInPlaquette = new int*[this->NbrSite];
  this->LargerIndicesInPlaquette = new int**[this->NbrSite];
  for (int i = 0; i < this->NbrSite; ++i)
  {
    this->NbrLargerIndicesInPlaquette[i] = new int[3];
    this->LargerIndicesInPlaquette[i] = new int*[3];
    for (int j = 0; j < 3; ++j)
      this->NbrLargerIndicesInPlaquette[i][j] = 0;
    
    int j = 0;
    bool tmpFlag;
    int tmp = 0;
    while (tmp < 3)
    {
      int tmpNbrIndex = 0;
      tmpFlag = false;
      for (int k = 0; k < 6; ++k)
      {
	if (this->ListIndicesPerPlaquette[j][k] > i)
	  tmpNbrIndex += 1;
	if (this->ListIndicesPerPlaquette[j][k] == i)
	  tmpFlag = true;
      }
      ++j;
      
      if (tmpFlag)
      {
	this->NbrLargerIndicesInPlaquette[i][tmp] = tmpNbrIndex;
	this->LargerIndicesInPlaquette[i][tmp] = new int[this->NbrLargerIndicesInPlaquette[i][tmp]];
	++tmp;
      }
    }
    
    j = 0;
    tmp = 0;
    while (tmp < 3)
    {
      int tmpNbrIndex = 0;
      tmpFlag = false;
      for (int k = 0; k < 6; ++k)
      {
	if ((this->ListIndicesPerPlaquette[j][k] > i) && (tmpNbrIndex < this->NbrLargerIndicesInPlaquette[i][tmp]))
	{
	  this->LargerIndicesInPlaquette[i][tmp][tmpNbrIndex] = this->ListIndicesPerPlaquette[j][k];
	  tmpNbrIndex += 1;
	}
	if (this->ListIndicesPerPlaquette[j][k] == i)
	  tmpFlag = true;
      }
      ++j;
      
      if (tmpFlag)
	++tmp;
    }
  }
}

// compute the charge on all three hexagons surrounding one site and find the largest one
//
// siteIndex = site index
// configuraton =hilbert space configuraton
// return value = maximum charge
  
int FermionOnHoneycombLatticeWithSpinRealSpacePlaquetteExclusionAnd2DTranslation::FindMaximumChargeSurroundingPlaquettes (int siteIndex, unsigned long configuration)
{
  int TmpIndex;
  int TmpCharge = 0;
  int TmpChargePlaquette;
  for (int i = 0; i < 3; ++i)
  {
    TmpChargePlaquette = 0;
    if (this->NbrLargerIndicesInPlaquette[siteIndex][i] > 2)
    {
      for (int j = 0; j < this->NbrLargerIndicesInPlaquette[siteIndex][i]; ++j)
      {
	TmpIndex = this->LargerIndicesInPlaquette[siteIndex][i][j];
	TmpChargePlaquette += ((configuration >> (TmpIndex << 1)) & 0x1ul);
	TmpChargePlaquette += ((configuration >> ((TmpIndex << 1) + 1)) & 0x1ul);
      }
      
      if (TmpChargePlaquette > 4)
	TmpChargePlaquette += ((configuration >> (siteIndex << 1)) & 0x1ul);
	TmpChargePlaquette += ((configuration >> ((siteIndex << 1) + 1)) & 0x1ul);
      
      if (TmpChargePlaquette > TmpCharge)
	TmpCharge = TmpChargePlaquette;
    }
  }
  return TmpCharge;
}