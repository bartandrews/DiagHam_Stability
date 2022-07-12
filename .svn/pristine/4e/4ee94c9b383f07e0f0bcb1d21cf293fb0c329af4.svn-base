////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//       in real space with translation invariance in two directions and      //
//                  Sz<->-Sz symmetry supporting up to 64 sites               //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 15/07/2016                      //
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
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceLong.h"
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

FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong ()
{
  this->SzParitySign = 1.0;
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

FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong (int nbrFermions, int nbrSite, bool minusSzParity, 
															  int xMomentum, int maxXMomentum,
															  int yMomentum, int  maxYMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->SzParitySign = 1.0;
  this->SzParity = ((ULONGLONG) 0x0ul);
  if (minusSzParity == true)
    {
      this->SzParitySign = -1.0;
      this->SzParity = ((ULONGLONG) 0x1ul);
    }
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = this->MaxXMomentum;

  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (((ULONGLONG) 0x1ul) << this->StateShift) - ((ULONGLONG) 0x1ul);

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = 2 * (this->NbrSite / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this->MaxMomentum) - this->StateXShift;
  this->XMomentumMask = (((ULONGLONG) 0x1ul) << this->StateXShift) - ((ULONGLONG) 0x1ul);

  this->MaxYMomentum =  maxYMomentum;
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = (2 * this->NbrSite) / this->StateXShift;
  this->StateYShift = 2 * (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (((ULONGLONG) 0x1ul) << this->StateYShift) - ((ULONGLONG) 0x1ul);
  this->YMomentumBlockMask = (((ULONGLONG) 0x1ul) << this->YMomentumBlockSize) - ((ULONGLONG) 0x1ul);  
  this->YMomentumFullMask = ((ULONGLONG) 0x0ul);
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
    }
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 

  this->NbrFermionsParity = (~((ULONGLONG) this->NbrFermions)) & ((ULONGLONG) 0x1ul);


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
	  UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
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
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// xMomentum = momentum sector in the x direction
// maxXMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong (int nbrFermions, int totalSpin, int nbrSite, bool minusSzParity, 
															  int xMomentum, int maxXMomentum,
															  int yMomentum, int  maxYMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalSpin = totalSpin;
  this->SzParitySign = 1.0;
  this->SzParity = ((ULONGLONG) 0x0ul);
  if (minusSzParity == true)
    {
      this->SzParitySign = -1.0;
      this->SzParity = ((ULONGLONG) 0x1ul);
    }
  this->NbrFermionsUp = (totalSpin + this->NbrFermions) >> 1;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = this->MaxXMomentum;

  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (((ULONGLONG) 0x1ul) << this->StateShift) - ((ULONGLONG) 0x1ul);

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = 2 * (this->NbrSite / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this->MaxMomentum) - this->StateXShift;
  this->XMomentumMask = (((ULONGLONG) 0x1ul) << this->StateXShift) - ((ULONGLONG) 0x1ul);

  this->MaxYMomentum =  maxYMomentum;
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = (2 * this->NbrSite) / this->StateXShift;
  this->StateYShift = 2 * (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (((ULONGLONG) 0x1ul) << this->StateYShift) - ((ULONGLONG) 0x1ul);
  this->YMomentumBlockMask = (((ULONGLONG) 0x1ul) << this->YMomentumBlockSize) - ((ULONGLONG) 0x1ul);  
  this->YMomentumFullMask = ((ULONGLONG) 0x0ul);
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
    }
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 

  this->NbrFermionsParity = (~((ULONGLONG) this->NbrFermions)) & ((ULONGLONG) 0x1ul);


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
	  UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
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

FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong(const FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong& fermions)
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
  this->NbrFermionsParity = fermions.NbrFermionsParity;

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

FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::~FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong& FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::operator = (const FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong& fermions)
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
  this->NbrFermionsParity = fermions.NbrFermionsParity;

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

AbstractHilbertSpace* FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::Clone()
{
  return new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong(*this);
}

// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::GenerateStates()
{
  this->StateDescription = new ULONGLONG [this->LargeHilbertSpaceDimension];
  if (this->SzFlag == false)
    this->RawGenerateStates(this->NbrFermions, this->NbrSite - 1, 0l);
  else
    this->RawGenerateStates(this->NbrFermions, this->NbrSite - 1, this->NbrFermionsUp, 0l);
  return this->CoreGenerateStates();
}


// generate all states corresponding to the constraints (core part of the method)
//
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::CoreGenerateStates()
{
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
  int NbrSzSymmetry;
  double AdditionalSign;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY, NbrSzSymmetry) == this->StateDescription[i]))
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      ++TmpLargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescription[i] = ((ULONGLONG) 0x0ul);
	    }
	}
      else
	{
	  this->StateDescription[i] = ((ULONGLONG) 0x0ul);
	}
    }
  if (TmpLargeHilbertSpaceDimension == 0l)
    return 0l;
  ULONGLONG* TmpStateDescription = new ULONGLONG [TmpLargeHilbertSpaceDimension];  
  this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
  this->ReorderingSign = new unsigned long [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->StateDescription[i] != ((ULONGLONG) 0x0ul))
	{
	  TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	  this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = this->FindOrbitSize(this->StateDescription[i]);
	  unsigned long& TmpSign = this->ReorderingSign[TmpLargeHilbertSpaceDimension];
	  TmpSign = ((ULONGLONG) 0x0ul);
	  int Index = 0;
	  for (int m = 0; m < this->MaxYMomentum; ++m)
	    {
	      int TmpNbrTranslationY = (this->MaxYMomentum - m) % this->MaxYMomentum;
	      for (int n = 0; n < this->MaxXMomentum; ++n)
		{
		  ULONGLONG TmpState =  this->StateDescription[i];
		  int TmpNbrTranslationX = (this->MaxXMomentum - n) % this->MaxXMomentum;
		  for (int n2 = 0; n2 < TmpNbrTranslationX; ++n2)
		    {
		      this->ApplySingleXTranslation(TmpState); 
		    }
		  for (int m2 = 0; m2 < TmpNbrTranslationY; ++m2)
		    {
		      this->ApplySingleYTranslation(TmpState); 
		    }		  
		  TmpSign |=  this->FindReorderingSign(TmpState, n, m, 0) << Index;
		  ++Index;
		  this->ApplySzSymmetry(TmpState);
		  TmpSign |=  this->FindReorderingSign(TmpState, n, m, 1) << Index;
		  ++Index;		  
		}
	    }
	  ++TmpLargeHilbertSpaceDimension;
	}
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;

  this->StateHighestBit = new int [TmpLargeHilbertSpaceDimension];  
  int CurrentMaxMomentum = (2 * this->MaxMomentum) + 1;
  while (((this->StateDescription[0] >> CurrentMaxMomentum) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
    --CurrentMaxMomentum;
  this->StateHighestBit[0] = CurrentMaxMomentum;
  for (long i = 1l; i < TmpLargeHilbertSpaceDimension; ++i)
    {
      while (((this->StateDescription[i] >> CurrentMaxMomentum) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
	--CurrentMaxMomentum;
      this->StateHighestBit[i] = CurrentMaxMomentum;
    }
  return TmpLargeHilbertSpaceDimension;
}


// convert a given state from the n-body basis with a fized Sz parity to the full n-body basis
//
// state = reference on the vector to convert
// targetNbodyBasis = reference on the nbody-basis where the final state will be expressed
// return value = converted vector

ComplexVector FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::ConvertToNbodyBasis(ComplexVector& state, FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong* targetNbodyBasis)
{
  ComplexVector TmpVector (targetNbodyBasis->LargeHilbertSpaceDimension, true);
  int TmpNbrTranslationX;
  int TmpNbrTranslationY;
  double TmpCoefficient;
  Complex** FourrierCoefficients = new Complex* [this->MomentumModulo];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      FourrierCoefficients[i] = new Complex [this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{
	  FourrierCoefficients[i][j] = Phase (2.0 * M_PI * ((double) (i * this->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * this->YMomentum) / ((double) this->MaxYMomentum)));
	}
    }
  for (long i= 0l; i < targetNbodyBasis->LargeHilbertSpaceDimension; ++i)
    {
      ULONGLONG TmpState = targetNbodyBasis->StateDescription[i];
      TmpCoefficient = 1.0;
      this->ProdATemporaryNbrStateInOrbit = targetNbodyBasis->NbrStateInOrbit[i];
      int TmpIndex = this->SymmetrizeAdAdResult(TmpState, TmpCoefficient, TmpNbrTranslationX, TmpNbrTranslationY);
      if (TmpIndex < this->HilbertSpaceDimension)
	{
	  TmpVector[i] = state[TmpIndex] * FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY] * TmpCoefficient;
	}
    }
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      delete[] FourrierCoefficients[i];
    }
  delete[] FourrierCoefficients;
  return TmpVector;
}

// convert a given state from the full n-body basis to the current n-body basis with a fized Sz parity 
//
// state = reference on the vector to convert
// inputNbodyBasis = reference on the nbody-basis where the inital state is expressed
// return value = converted vector

ComplexVector FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::ConvertFromNbodyBasis(ComplexVector& state, FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong* inputNbodyBasis)
{
  ComplexVector TmpVector (this->LargeHilbertSpaceDimension, true);
  int TmpNbrTranslationX;
  int TmpNbrTranslationY;
  double TmpCoefficient;
  Complex** FourrierCoefficients = new Complex* [this->MomentumModulo];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      FourrierCoefficients[i] = new Complex [this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{
	  FourrierCoefficients[i][j] = Phase (2.0 * M_PI * ((double) (i * this->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * this->YMomentum) / ((double) this->MaxYMomentum)));
	}
    }
  for (long i= 0l; i < inputNbodyBasis->LargeHilbertSpaceDimension; ++i)
    {
      ULONGLONG TmpState = inputNbodyBasis->StateDescription[i];
      TmpCoefficient = 1.0;
      this->ProdATemporaryNbrStateInOrbit = inputNbodyBasis->NbrStateInOrbit[i];
      int TmpIndex = this->SymmetrizeAdAdResult(TmpState, TmpCoefficient, TmpNbrTranslationX, TmpNbrTranslationY);
      if (TmpIndex < this->HilbertSpaceDimension)
	{
	  TmpVector[TmpIndex] += state[i] * FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY] * TmpCoefficient;
	}
    }
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      delete[] FourrierCoefficients[i];
    }
  delete[] FourrierCoefficients;
  return TmpVector;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in given momentum and Sz sectors.
//
// nbrParticleSector = number of particles that belong to the subsytem
// szParitySector = Sz parity sector (can be either -1 or +1)
// kxSector = subsystem momentum along the x direction
// kySector = subsystem momentum along the x direction
// szSector  = twice the total Sz value of the subsytem
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int szSector, int szParitySector, 
															    int kxSector, int kySector, 
															    ComplexVector& groundState, AbstractArchitecture* architecture)
{
  int TmpCurrentSzParitySector = 1;
  if (this->SzParitySign < 0.0)
    {
      TmpCurrentSzParitySector = -1;
    }
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0) && (szSector == 0) && (szParitySector == 1))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      if ((kxSector == this->XMomentum) && (kySector == this->YMomentum) && (this->TotalSpin == szSector) & (szParitySector == TmpCurrentSzParitySector))
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
  int ComplementarySzParitySector = TmpCurrentSzParitySector * szParitySector;
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong SubsytemSpace (nbrParticleSector, szSector, this->NbrSite, (szParitySector == -1), 
									     kxSector, this->MaxXMomentum, kySector, this->MaxYMomentum);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong ComplementarySpace (ComplementaryNbrParticles, ComplementarySzSector, this->NbrSite, 
										  (ComplementarySzParitySector == -1), ComplementaryKxMomentum, this->MaxXMomentum, 
										  ComplementaryKyMomentum , this->MaxYMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;
  architecture->SetDimension(ComplementarySpace.GetHilbertSpaceDimension());
  FQHETorusParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, (ParticleOnTorusWithSpinAndMagneticTranslations*) &ComplementarySpace, groundState, TmpDensityMatrix);
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
  
// evaluate a density matrix of a subsystem of the whole system described by a given sum of projectors, using particle partition. The density matrix is only evaluated in given momentum and Sz sectors.
//
// nbrParticleSector = number of particles that belong to the subsytem
// szSector  = twice the total Sz value of the subsytem
// szParitySector = Sz parity sector (can be either -1 or +1)
// kxSector = subsystem momentum along the x direction
// kySector = subsystem momentum along the x direction
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// architecture = pointer to the architecture to use parallelized algorithm
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int szSector, int szParitySector, 
															    int kxSector, int kySector, 
															    int nbrGroundStates, ComplexVector* groundStates, 
															    double* weights, 
															    AbstractArchitecture* architecture)
{
  int TmpCurrentSzParitySector = 1;
  if (this->SzParitySign < 0.0)
    {
      TmpCurrentSzParitySector = -1;
    }
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0) && (szSector == 0) && (szParitySector == 1))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      if ((kxSector == this->XMomentum) && (kySector == this->YMomentum) && (this->TotalSpin == szSector) & (szParitySector == TmpCurrentSzParitySector))
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
  int ComplementarySzParitySector = TmpCurrentSzParitySector * szParitySector;
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong SubsytemSpace (nbrParticleSector, szSector, this->NbrSite, (szParitySector == -1), 
									     kxSector, this->MaxXMomentum, kySector, this->MaxYMomentum);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong ComplementarySpace (ComplementaryNbrParticles, ComplementarySzSector, this->NbrSite, 
										  (ComplementarySzParitySector == -1), ComplementaryKxMomentum, this->MaxXMomentum, 
										  ComplementaryKyMomentum , this->MaxYMomentum);
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
  
// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, 
														     ParticleOnTorusWithSpinAndMagneticTranslations* complementaryHilbertSpace,  
														     ParticleOnTorusWithSpinAndMagneticTranslations* destinationHilbertSpace,
														     ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
  Complex** FourrierCoefficientsDestination = new Complex* [this->MomentumModulo];
  Complex** FourrierCoefficients = new Complex* [this->MomentumModulo];
  long TmpNbrNonZeroElements = 0l;
  if ((destinationHilbertSpace->GetTotalSpin() != 0) || (this->TotalSpin != 0))
    {
      FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong* TmpDestinationHilbertSpace =  (FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong*) destinationHilbertSpace;
      FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong* TmpHilbertSpace = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong*) complementaryHilbertSpace;
      FermionOnLatticeWithSpinRealSpaceLong* TmpDestinationFullHilbertSpace = 0;
      if (TmpDestinationHilbertSpace->SzFlag == false)
	{
	  TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpaceLong(TmpDestinationHilbertSpace->NbrFermions,
										     TmpDestinationHilbertSpace->NbrSite);
	}
      else
	{
	  TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpaceLong(TmpDestinationHilbertSpace->NbrFermions,
										     TmpDestinationHilbertSpace->TotalSpin,
										     TmpDestinationHilbertSpace->NbrSite);
	}
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
      BinomialCoefficients TmpBinomial (this->NbrFermions);
      double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
      this->ProdATemporaryNbrStateInOrbit = 1;
      TmpDestinationHilbertSpace->ProdATemporaryNbrStateInOrbit = 1;

      Complex* TmpDestinationFactors = new Complex[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
      int* TmpDestinationRealIndices = new int[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
      for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  ULONGLONG TmpCanonicalState2 = TmpDestinationFullHilbertSpace->StateDescription[j];
	  int TmpDestinationNbrTranslationX;
	  int TmpDestinationNbrTranslationY;
	  double TmpDestinationCoefficient = 1.0;
	  int RealDestinationIndex = TmpDestinationHilbertSpace->SymmetrizeAdAdResult(TmpCanonicalState2, TmpDestinationCoefficient, TmpDestinationNbrTranslationX, TmpDestinationNbrTranslationY);
	  TmpDestinationRealIndices[j] = RealDestinationIndex;
	  if (RealDestinationIndex < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
	    {
	      TmpDestinationFactors[j] = (TmpDestinationCoefficient * TmpInvBinomial) * FourrierCoefficientsDestination[TmpDestinationNbrTranslationX][TmpDestinationNbrTranslationY];
	    }
	}

      for (; minIndex < MaxIndex; ++minIndex)    
	{
	  int Pos = 0;
	  ULONGLONG TmpState = TmpHilbertSpace->StateDescription[minIndex];
	  double TmpRescalingFactor = sqrt((double)  TmpHilbertSpace->NbrStateInOrbit[minIndex]);
	  for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
	    {
	      ULONGLONG TmpState2 = TmpDestinationFullHilbertSpace->StateDescription[j];
	      if ((TmpState & TmpState2) == ((ULONGLONG) 0x0ul))
		{
		  if (TmpDestinationRealIndices[j] < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
		    {
		      int TmpNbrTranslationX;
		      int TmpNbrTranslationY;
		      double TmpCoefficient = TmpRescalingFactor;
		      ULONGLONG TmpState3 = TmpState | TmpState2;
		      int TmpPos = this->SymmetrizeAdAdResult(TmpState3, TmpCoefficient, TmpNbrTranslationX, TmpNbrTranslationY);
		      
		      if (TmpPos < this->HilbertSpaceDimension)
			{      
 			  Complex Coefficient = TmpCoefficient * TmpDestinationFactors[j] * FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY];
			  ULONGLONG Sign = ((ULONGLONG) 0x0ul);
			  int Pos2 = 2 * TmpDestinationHilbertSpace->NbrSite - 1;
			  ULONGLONG TmpState22 = TmpState2;
			  while ((Pos2 > 0) && (TmpState22 != ((ULONGLONG) 0x0ul)))
			    {
			      while (((TmpState22 >> Pos2) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
				--Pos2;
			      TmpState3 = TmpState & (( ((ULONGLONG) 0x1ul) << (Pos2 + 1)) - ((ULONGLONG) 0x1ul));
#ifdef __128_BIT_LONGLONG__
			      TmpState3 ^= TmpState3 >> 64;
#endif	
			      TmpState3 ^= TmpState3 >> 32;
			      TmpState3 ^= TmpState3 >> 16;
			      TmpState3 ^= TmpState3 >> 8;
			      TmpState3 ^= TmpState3 >> 4;
			      TmpState3 ^= TmpState3 >> 2;
			      TmpState3 ^= TmpState3 >> 1;
			      Sign ^= TmpState3;
			      TmpState22 &= ~(((ULONGLONG) 0x1ul) << Pos2);
			      --Pos2;
			    }
			  if ((Sign & ((ULONGLONG) 0x1ul)) != ((ULONGLONG) 0x0ul))		  
			    Coefficient *= -1.0;
			  TmpEntanglementMatrix[TmpDestinationRealIndices[j]][minIndex - TmpMinIndex] += groundState[TmpPos] * Coefficient;
			}
		    }
		}
	    }
	}
      for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
	{
	  for (int j = i; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	    {
	      densityMatrix->UnsafeAddToMatrixElement(i, j, TmpEntanglementMatrix[j] * TmpEntanglementMatrix[i]);
	      ++TmpNbrNonZeroElements;
	    }
	}
      delete TmpDestinationFullHilbertSpace;
      delete[] TmpDestinationFactors;
      delete[] TmpDestinationRealIndices;
    }
  else
    {
      FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong* TmpDestinationHilbertSpace =  (FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong*) destinationHilbertSpace;
      FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong* TmpHilbertSpace = (FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong*) complementaryHilbertSpace;
      FermionOnLatticeWithSpinRealSpaceLong* TmpDestinationFullHilbertSpace = 0;
      if (TmpDestinationHilbertSpace->SzFlag == false)
	{
	  TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpaceLong(TmpDestinationHilbertSpace->NbrFermions,
										     TmpDestinationHilbertSpace->NbrSite);
	}
      else
	{
	  TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpaceLong(TmpDestinationHilbertSpace->NbrFermions,
										     TmpDestinationHilbertSpace->TotalSpin,
										     TmpDestinationHilbertSpace->NbrSite);
	}
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
      BinomialCoefficients TmpBinomial (this->NbrFermions);
      double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
      this->ProdATemporaryNbrStateInOrbit = 1;
      TmpDestinationHilbertSpace->ProdATemporaryNbrStateInOrbit = 1;

      Complex* TmpDestinationFactors = new Complex[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
      int* TmpDestinationRealIndices = new int[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
      for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  ULONGLONG TmpCanonicalState2 = TmpDestinationFullHilbertSpace->StateDescription[j];
	  int TmpDestinationNbrTranslationX;
	  int TmpDestinationNbrTranslationY;
	  double TmpDestinationCoefficient = 1.0;
	  int RealDestinationIndex = TmpDestinationHilbertSpace->SymmetrizeAdAdResult(TmpCanonicalState2, TmpDestinationCoefficient, TmpDestinationNbrTranslationX, TmpDestinationNbrTranslationY);
	  TmpDestinationRealIndices[j] = RealDestinationIndex;
	  if (RealDestinationIndex < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
	    {
	      TmpDestinationFactors[j] = (TmpDestinationCoefficient * TmpInvBinomial) * FourrierCoefficientsDestination[TmpDestinationNbrTranslationX][TmpDestinationNbrTranslationY];
	    }
	}

      for (; minIndex < MaxIndex; ++minIndex)    
	{
	  int Pos = 0;
	  ULONGLONG TmpState = TmpHilbertSpace->StateDescription[minIndex];
	  double TmpRescalingFactor = sqrt((double)  TmpHilbertSpace->NbrStateInOrbit[minIndex]);
	  for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
	    {
	      ULONGLONG TmpState2 = TmpDestinationFullHilbertSpace->StateDescription[j];
	      if ((TmpState & TmpState2) == ((ULONGLONG) 0x0ul))
		{
		  if (TmpDestinationRealIndices[j] < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
		    {
		      int TmpNbrTranslationX;
		      int TmpNbrTranslationY;
		      double TmpCoefficient = TmpRescalingFactor;
		      ULONGLONG TmpState3 = TmpState | TmpState2;
		      int TmpPos = this->SymmetrizeAdAdResult(TmpState3, TmpCoefficient, TmpNbrTranslationX, TmpNbrTranslationY);
		      
		      if (TmpPos < this->HilbertSpaceDimension)
			{      
 			  Complex Coefficient = TmpCoefficient * TmpDestinationFactors[j] * FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY];
			  ULONGLONG Sign = ((ULONGLONG) 0x0ul);
			  int Pos2 = 2 * TmpDestinationHilbertSpace->NbrSite - 1;
			  ULONGLONG TmpState22 = TmpState2;
			  while ((Pos2 > 0) && (TmpState22 != ((ULONGLONG) 0x0ul)))
			    {
			      while (((TmpState22 >> Pos2) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
				--Pos2;
			      TmpState3 = TmpState & ((((ULONGLONG) 0x1ul) << (Pos2 + 1)) - ((ULONGLONG) 0x1ul));
#ifdef __128_BIT_LONGLONG__
			      TmpState3 ^= TmpState3 >> 64;
#endif	
			      TmpState3 ^= TmpState3 >> 32;
			      TmpState3 ^= TmpState3 >> 16;
			      TmpState3 ^= TmpState3 >> 8;
			      TmpState3 ^= TmpState3 >> 4;
			      TmpState3 ^= TmpState3 >> 2;
			      TmpState3 ^= TmpState3 >> 1;
			      Sign ^= TmpState3;
			      TmpState22 &= ~(((ULONGLONG) 0x1ul) << Pos2);
			      --Pos2;
			    }
			  if ((Sign & ((ULONGLONG) 0x1ul)) != ((ULONGLONG) 0x0ul))		  
			    Coefficient *= -1.0;
			  TmpEntanglementMatrix[TmpDestinationRealIndices[j]][minIndex - TmpMinIndex] += groundState[TmpPos] * Coefficient;
			}
		    }
		}
	    }
	}
      for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
	{
	  for (int j = i; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	    {
	      densityMatrix->UnsafeAddToMatrixElement(i, j, TmpEntanglementMatrix[j] * TmpEntanglementMatrix[i]);
	      ++TmpNbrNonZeroElements;
	    }
	}
      delete TmpDestinationFullHilbertSpace;
      delete[] TmpDestinationFactors;
      delete[] TmpDestinationRealIndices;
    }
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      delete[] FourrierCoefficientsDestination[i];
      delete[] FourrierCoefficients[i];
    }
  delete[] FourrierCoefficientsDestination;
  delete[] FourrierCoefficients;
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

long FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorusWithSpinAndMagneticTranslations* complementaryHilbertSpace,
														     ParticleOnTorusWithSpinAndMagneticTranslations* destinationHilbertSpace,
														     int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix* densityMatrix)
{
  Complex** FourrierCoefficientsDestination = new Complex* [this->MomentumModulo];
  Complex** FourrierCoefficients = new Complex* [this->MomentumModulo];
  long TmpNbrNonZeroElements = 0l;
  if ((destinationHilbertSpace->GetTotalSpin() != 0) || (this->TotalSpin != 0))
    {
      FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong* TmpDestinationHilbertSpace =  (FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong*) destinationHilbertSpace;
      FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong* TmpHilbertSpace = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong*) complementaryHilbertSpace;
      FermionOnLatticeWithSpinRealSpaceLong* TmpDestinationFullHilbertSpace = 0;
      if (TmpDestinationHilbertSpace->SzFlag == false)
	{
	  TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpaceLong(TmpDestinationHilbertSpace->NbrFermions,
										     TmpDestinationHilbertSpace->NbrSite);
	}
      else
	{
	  TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpaceLong(TmpDestinationHilbertSpace->NbrFermions,
										     TmpDestinationHilbertSpace->TotalSpin,
										     TmpDestinationHilbertSpace->NbrSite);
	}
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
      
      ComplexMatrix* TmpEntanglementMatrices = new ComplexMatrix[nbrGroundStates];
      for (int k = 0; k < nbrGroundStates; ++k)
	TmpEntanglementMatrices[k] = ComplexMatrix(nbrIndex, TmpDestinationHilbertSpace->GetHilbertSpaceDimension(), true);
      int TmpMinIndex = minIndex;
      int MaxIndex = minIndex + nbrIndex;
      BinomialCoefficients TmpBinomial (this->NbrFermions);
      double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
      this->ProdATemporaryNbrStateInOrbit = 1;
      TmpDestinationHilbertSpace->ProdATemporaryNbrStateInOrbit = 1;

      Complex* TmpDestinationFactors = new Complex[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
      int* TmpDestinationRealIndices = new int[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
      for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  ULONGLONG TmpCanonicalState2 = TmpDestinationFullHilbertSpace->StateDescription[j];
	  int TmpDestinationNbrTranslationX;
	  int TmpDestinationNbrTranslationY;
	  double TmpDestinationCoefficient = 1.0;
	  int RealDestinationIndex = TmpDestinationHilbertSpace->SymmetrizeAdAdResult(TmpCanonicalState2, TmpDestinationCoefficient, TmpDestinationNbrTranslationX, TmpDestinationNbrTranslationY);
	  TmpDestinationRealIndices[j] = RealDestinationIndex;
	  if (RealDestinationIndex < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
	    {
	      TmpDestinationFactors[j] = (TmpDestinationCoefficient * TmpInvBinomial) * FourrierCoefficientsDestination[TmpDestinationNbrTranslationX][TmpDestinationNbrTranslationY];
	    }
	}

      for (; minIndex < MaxIndex; ++minIndex)    
	{
	  int Pos = 0;
	  ULONGLONG TmpState = TmpHilbertSpace->StateDescription[minIndex];
	  double TmpRescalingFactor = sqrt((double)  TmpHilbertSpace->NbrStateInOrbit[minIndex]);
	  for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
	    {
	      ULONGLONG TmpState2 = TmpDestinationFullHilbertSpace->StateDescription[j];
	      if ((TmpState & TmpState2) == ((ULONGLONG) 0x0ul))
		{
		  if (TmpDestinationRealIndices[j] < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
		    {
		      int TmpNbrTranslationX;
		      int TmpNbrTranslationY;
		      double TmpCoefficient = TmpRescalingFactor;
		      ULONGLONG TmpState3 = TmpState | TmpState2;
		      int TmpPos = this->SymmetrizeAdAdResult(TmpState3, TmpCoefficient, TmpNbrTranslationX, TmpNbrTranslationY);
		      
		      if (TmpPos < this->HilbertSpaceDimension)
			{      
 			  Complex Coefficient = TmpCoefficient * TmpDestinationFactors[j] * FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY];
			  ULONGLONG Sign = ((ULONGLONG) 0x0ul);
			  int Pos2 = 2 * TmpDestinationHilbertSpace->NbrSite - 1;
			  ULONGLONG TmpState22 = TmpState2;
			  while ((Pos2 > 0) && (TmpState22 != ((ULONGLONG) 0x0ul)))
			    {
			      while (((TmpState22 >> Pos2) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
				--Pos2;
			      TmpState3 = TmpState & ((((ULONGLONG) 0x1ul) << (Pos2 + 1)) - ((ULONGLONG) 0x1ul));
#ifdef __128_BIT_LONGLONG__
			      TmpState3 ^= TmpState3 >> 64;
#endif	
			      TmpState3 ^= TmpState3 >> 32;
			      TmpState3 ^= TmpState3 >> 16;
			      TmpState3 ^= TmpState3 >> 8;
			      TmpState3 ^= TmpState3 >> 4;
			      TmpState3 ^= TmpState3 >> 2;
			      TmpState3 ^= TmpState3 >> 1;
			      Sign ^= TmpState3;
			      TmpState22 &= ~(((ULONGLONG) 0x1ul) << Pos2);
			      --Pos2;
			    }
			  if ((Sign & ((ULONGLONG) 0x1ul)) != ((ULONGLONG) 0x0ul))		  
			    Coefficient *= -1.0;
			  for (int k = 0; k < nbrGroundStates; ++k)
			    TmpEntanglementMatrices[k][TmpDestinationRealIndices[j]][minIndex - TmpMinIndex] += groundStates[k][TmpPos] * Coefficient;
			}
		    }
		}
	    }
	}
      for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
	{
	  for (int j = i; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	    {
	      Complex Tmp = 0.0;
	      for (int k = 0; k < nbrGroundStates; ++k)
		Tmp += weights[k] * TmpEntanglementMatrices[k][j] * TmpEntanglementMatrices[k][i];
	      densityMatrix->UnsafeAddToMatrixElement(i, j, Tmp);
	      ++TmpNbrNonZeroElements;
	    }
	}
      delete TmpDestinationFullHilbertSpace;
      delete[] TmpDestinationFactors;
      delete[] TmpDestinationRealIndices;
      delete[] TmpEntanglementMatrices;
    }
  else
    {
      FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong* TmpDestinationHilbertSpace =  (FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong*) destinationHilbertSpace;
      FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong* TmpHilbertSpace = (FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong*) complementaryHilbertSpace;
      FermionOnLatticeWithSpinRealSpaceLong* TmpDestinationFullHilbertSpace = 0;
      if (TmpDestinationHilbertSpace->SzFlag == false)
	{
	  TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpaceLong(TmpDestinationHilbertSpace->NbrFermions,
										     TmpDestinationHilbertSpace->NbrSite);
	}
      else
	{
	  TmpDestinationFullHilbertSpace = new FermionOnLatticeWithSpinRealSpaceLong(TmpDestinationHilbertSpace->NbrFermions,
										     TmpDestinationHilbertSpace->TotalSpin,
										     TmpDestinationHilbertSpace->NbrSite);
	}
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
      
      ComplexMatrix* TmpEntanglementMatrices = new ComplexMatrix[nbrGroundStates];
      for (int k = 0; k < nbrGroundStates; ++k)
	TmpEntanglementMatrices[k] = ComplexMatrix(nbrIndex, TmpDestinationHilbertSpace->GetHilbertSpaceDimension(), true);
      int TmpMinIndex = minIndex;
      int MaxIndex = minIndex + nbrIndex;
      BinomialCoefficients TmpBinomial (this->NbrFermions);
      double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
      this->ProdATemporaryNbrStateInOrbit = 1;
      TmpDestinationHilbertSpace->ProdATemporaryNbrStateInOrbit = 1;

      Complex* TmpDestinationFactors = new Complex[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
      int* TmpDestinationRealIndices = new int[TmpDestinationFullHilbertSpace->GetHilbertSpaceDimension()];
      for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  ULONGLONG TmpCanonicalState2 = TmpDestinationFullHilbertSpace->StateDescription[j];
	  int TmpDestinationNbrTranslationX;
	  int TmpDestinationNbrTranslationY;
	  double TmpDestinationCoefficient = 1.0;
	  int RealDestinationIndex = TmpDestinationHilbertSpace->SymmetrizeAdAdResult(TmpCanonicalState2, TmpDestinationCoefficient, TmpDestinationNbrTranslationX, TmpDestinationNbrTranslationY);
	  TmpDestinationRealIndices[j] = RealDestinationIndex;
	  if (RealDestinationIndex < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
	    {
	      TmpDestinationFactors[j] = (TmpDestinationCoefficient * TmpInvBinomial) * FourrierCoefficientsDestination[TmpDestinationNbrTranslationX][TmpDestinationNbrTranslationY];
	    }
	}

      for (; minIndex < MaxIndex; ++minIndex)    
	{
	  int Pos = 0;
	  ULONGLONG TmpState = TmpHilbertSpace->StateDescription[minIndex];
	  double TmpRescalingFactor = sqrt((double)  TmpHilbertSpace->NbrStateInOrbit[minIndex]);
	  for (int j = 0; j < TmpDestinationFullHilbertSpace->HilbertSpaceDimension; ++j)
	    {
	      ULONGLONG TmpState2 = TmpDestinationFullHilbertSpace->StateDescription[j];
	      if ((TmpState & TmpState2) == ((ULONGLONG) 0x0ul))
		{
		  if (TmpDestinationRealIndices[j] < TmpDestinationHilbertSpace->GetHilbertSpaceDimension())
		    {
		      int TmpNbrTranslationX;
		      int TmpNbrTranslationY;
		      double TmpCoefficient = TmpRescalingFactor;
		      ULONGLONG TmpState3 = TmpState | TmpState2;
		      int TmpPos = this->SymmetrizeAdAdResult(TmpState3, TmpCoefficient, TmpNbrTranslationX, TmpNbrTranslationY);
		      
		      if (TmpPos < this->HilbertSpaceDimension)
			{      
 			  Complex Coefficient = TmpCoefficient * TmpDestinationFactors[j] * FourrierCoefficients[TmpNbrTranslationX][TmpNbrTranslationY];
			  ULONGLONG Sign = ((ULONGLONG) 0x0ul);
			  int Pos2 = 2 * TmpDestinationHilbertSpace->NbrSite - 1;
			  ULONGLONG TmpState22 = TmpState2;
			  while ((Pos2 > 0) && (TmpState22 != ((ULONGLONG) 0x0ul)))
			    {
			      while (((TmpState22 >> Pos2) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
				--Pos2;
			      TmpState3 = TmpState & ((((ULONGLONG) 0x1ul) << (Pos2 + 1)) -  ((ULONGLONG) 0x1ul));
#ifdef __128_BIT_LONGLONG__
			      TmpState3 ^= TmpState3 >> 64;
#endif	
			      TmpState3 ^= TmpState3 >> 32;
			      TmpState3 ^= TmpState3 >> 16;
			      TmpState3 ^= TmpState3 >> 8;
			      TmpState3 ^= TmpState3 >> 4;
			      TmpState3 ^= TmpState3 >> 2;
			      TmpState3 ^= TmpState3 >> 1;
			      Sign ^= TmpState3;
			      TmpState22 &= ~(((ULONGLONG) 0x1ul) << Pos2);
			      --Pos2;
			    }
			  if ((Sign & ((ULONGLONG) 0x1ul)) != ((ULONGLONG) 0x0ul))		  
			    Coefficient *= -1.0;
			  for (int k = 0; k < nbrGroundStates; ++k)
			    TmpEntanglementMatrices[k][TmpDestinationRealIndices[j]][minIndex - TmpMinIndex] += groundStates[k][TmpPos] * Coefficient;
			}
		    }
		}
	    }
	}
      for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
	{
	  for (int j = i; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	    {
	      Complex Tmp = 0.0;
	      for (int k = 0; k < nbrGroundStates; ++k)
		Tmp += weights[k] * TmpEntanglementMatrices[k][j] * TmpEntanglementMatrices[k][i];
	      densityMatrix->UnsafeAddToMatrixElement(i, j, Tmp);
	      ++TmpNbrNonZeroElements;
	    }
	}
      delete TmpDestinationFullHilbertSpace;
      delete[] TmpDestinationFactors;
      delete[] TmpDestinationRealIndices;
      delete[] TmpEntanglementMatrices;
    }
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      delete[] FourrierCoefficientsDestination[i];
      delete[] FourrierCoefficients[i];
    }
  delete[] FourrierCoefficientsDestination;
  delete[] FourrierCoefficients;
  return TmpNbrNonZeroElements;
}

