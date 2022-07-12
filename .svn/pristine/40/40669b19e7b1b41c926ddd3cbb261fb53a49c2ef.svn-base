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
//                   and an exclusion rule for neighboring sites              //
//                                                                            //
//                                                                            //
//                        last modification : 11/02/2018                      //
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
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceWithExclusion.h"
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

FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion::FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion ()
{
  this->NbrExcludedSiteMasks = 0;
  this->ExcludedSiteMasks = 0;
  this->ExcludedSiteCenters = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// xMomentum = momentum sector in the x direction
// xTranslation = translation that has to be applied on the site index to connect two sites with a translation in the x direction
// yMomentum = momentum sector in the y direction
// yPeriodicity = periodicity in the y direction with respect to site numbering 
// excludedSites = array of the linearized coordinates of the excluded sites (for the unit cell at (0,0), first index being the orbital index)
// nbrExcludedSites = number of excluded sites per orbital within the unit cell at (0,0)
// memory = amount of memory granted for precalculations

FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion::FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion (int nbrFermions, int nbrSite, int xMomentum,  int maxXMomentum,
														int yMomentum,  int maxYMomentum, int** excludedSites, 
														int* nbrExcludedSites, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum =  maxXMomentum;
  this->MomentumModulo =  this->MaxXMomentum;

  this->MomentumIncrement = (this->NbrFermions * this->StateShift / 2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = this->NbrSite / this->MaxXMomentum;
  this->ComplementaryStateXShift = this->MaxMomentum - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->MaxYMomentum =  maxYMomentum;
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = this->NbrSite / this->StateXShift;
  this->StateYShift = (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
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

  this->NbrExcludedSiteMasks = this->StateYShift * this->MaxYMomentum * this->MaxXMomentum;
  this->ExcludedSiteMasks = new unsigned long [this->NbrExcludedSiteMasks];
  this->ExcludedSiteCenters = new unsigned long [this->NbrExcludedSiteMasks];
  this->NbrExcludedSiteMasks = 0;
  for (int k = 0; k < this->StateYShift; ++k)
    {
      this->ExcludedSiteCenters[this->NbrExcludedSiteMasks] = 0x1ul << k;
      this->ExcludedSiteMasks[this->NbrExcludedSiteMasks] = 0x0ul;
      for (int l = 0; l < nbrExcludedSites[k]; ++l)
	{
	  this->ExcludedSiteMasks[this->NbrExcludedSiteMasks] |= 0x1ul << excludedSites[k][l];
	}
      ++this->NbrExcludedSiteMasks;
    }
  this->NbrExcludedSiteMasks = 0;
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{
	  int TmpIndex = 0;
	  for (int k = 0; k < this->StateYShift; ++k)
	    {
	      unsigned long TmpExcludedSiteCenter = this->ExcludedSiteCenters[k];
	      unsigned long TmpExcludedSiteMask = this->ExcludedSiteMasks[k];
	      for (int i2 = 0; i2 < i; ++i2)
		{
		  this->ApplySingleXTranslation(TmpExcludedSiteCenter);
		  this->ApplySingleXTranslation(TmpExcludedSiteMask);
		}
	      for (int j2 = 0; j2 < j; ++j2)
		{
		  this->ApplySingleYTranslation(TmpExcludedSiteCenter);
		  this->ApplySingleYTranslation(TmpExcludedSiteMask);
		}
	      this->ExcludedSiteCenters[this->NbrExcludedSiteMasks] = TmpExcludedSiteCenter;
	      this->ExcludedSiteMasks[this->NbrExcludedSiteMasks] = TmpExcludedSiteMask;
	      ++this->NbrExcludedSiteMasks;
	    }
	}      
    }

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
	  this->StateMaxMomentum = new int [this->LargeHilbertSpaceDimension];  
	  int CurrentMaxMomentum = this->MaxMomentum;
	  while (((this->StateDescription[0] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
	    --CurrentMaxMomentum;
	  this->StateMaxMomentum[0] = CurrentMaxMomentum;
	  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      while (((this->StateDescription[i] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
		--CurrentMaxMomentum;
	      this->StateMaxMomentum[i] = CurrentMaxMomentum;
	    }
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
#endif
	}
    }
}

// constructor from prebuilt exclusion rules
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// xMomentum = momentum sector in the x direction
// xTranslation = translation that has to be applied on the site index to connect two sites with a translation in the x direction
// yMomentum = momentum sector in the y direction
// yPeriodicity = periodicity in the y direction with respect to site numbering 
// excludedSiteMasks = masks used to detected excluded sites around a given position
// excludedSiteCenters = masks used to indicate the site around which a given exclusion rule is defined
// nbrExcludedSiteMasks = number of masks in ExcludedSiteMasks
// memory = amount of memory granted for precalculations

FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion::FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion (int nbrFermions, int nbrSite, int xMomentum,  int maxXMomentum,
														int yMomentum,  int maxYMomentum, unsigned long* excludedSiteMasks,
														unsigned long* excludedSiteCenters, int nbrExcludedSiteMasks,
														unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum =  maxXMomentum;
  this->MomentumModulo =  this->MaxXMomentum;

  this->MomentumIncrement = (this->NbrFermions * this->StateShift / 2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = this->NbrSite / this->MaxXMomentum;
  this->ComplementaryStateXShift = this->MaxMomentum - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->MaxYMomentum =  maxYMomentum;
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = this->NbrSite / this->StateXShift;
  this->StateYShift = (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
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

  this->NbrExcludedSiteMasks = nbrExcludedSiteMasks;
  this->ExcludedSiteMasks = new unsigned long [this->NbrExcludedSiteMasks];
  this->ExcludedSiteCenters = new unsigned long [this->NbrExcludedSiteMasks];
  for (int k = 0; k < this->NbrExcludedSiteMasks; ++k)
    {
      this->ExcludedSiteCenters[k] = excludedSiteCenters[k];
      this->ExcludedSiteMasks[k] = excludedSiteMasks[k];
    }

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
	  this->StateMaxMomentum = new int [this->LargeHilbertSpaceDimension];  
	  int CurrentMaxMomentum = this->MaxMomentum;
	  while (((this->StateDescription[0] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
	    --CurrentMaxMomentum;
	  this->StateMaxMomentum[0] = CurrentMaxMomentum;
	  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      while (((this->StateDescription[i] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
		--CurrentMaxMomentum;
	      this->StateMaxMomentum[i] = CurrentMaxMomentum;
	    }
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
#endif
	}
    }
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion::FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion(const FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion& fermions)
{
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
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

  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->NbrExcludedSiteMasks = fermions.NbrExcludedSiteMasks;
  this->ExcludedSiteMasks = new unsigned long [this->NbrExcludedSiteMasks];
  this->ExcludedSiteCenters = new unsigned long [this->NbrExcludedSiteMasks];
  for (int i = 0; i < this->NbrExcludedSiteMasks; ++i)
    {
      this->ExcludedSiteMasks[i] = fermions.ExcludedSiteMasks[i];
      this->ExcludedSiteCenters[i] = fermions.ExcludedSiteCenters[i];
    }

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;

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

}

// destructor
//

FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion::~FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion ()
{
  delete[] this->ExcludedSiteMasks;
  delete[] this->ExcludedSiteCenters;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion& FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion::operator = (const FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;

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

      delete[] this->ExcludedSiteMasks;
      delete[] this->ExcludedSiteCenters;
    }
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
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

  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->NbrExcludedSiteMasks = fermions.NbrExcludedSiteMasks;
  this->ExcludedSiteMasks = new unsigned long [this->NbrExcludedSiteMasks];
  this->ExcludedSiteCenters = new unsigned long [this->NbrExcludedSiteMasks];
  for (int i = 0; i < this->NbrExcludedSiteMasks; ++i)
    {
      this->ExcludedSiteMasks[i] = fermions.ExcludedSiteMasks[i];
      this->ExcludedSiteCenters[i] = fermions.ExcludedSiteCenters[i];
    }

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = this->LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;

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

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion::Clone()
{
  return new FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion(*this);
}


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(this->NbrFermions, this->NbrSite - 1, 0l);
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      int TmpPosition = 0;
      unsigned long TmpState = this->StateDescription[i];
      while ((TmpPosition < this->NbrExcludedSiteMasks) && (((TmpState & this->ExcludedSiteCenters[TmpPosition]) == 0x0ul)
							    || ((TmpState & this->ExcludedSiteMasks[TmpPosition]) == 0x0ul)))
	{
	  ++TmpPosition;
	}
      if ((TmpPosition == this->NbrExcludedSiteMasks) &&
	  (this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY) == this->StateDescription[i]))
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
  return TmpLargeHilbertSpaceDimension;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = kx sector in which the density matrix has to be evaluated 
// kySector = kx sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture)
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
  FermionOnLatticeRealSpaceAnd2DTranslationWithExclusion SubsytemSpace (nbrParticleSector, this->NbrSite, kxSector, this->MaxXMomentum, kySector, this->MaxYMomentum,
									this->ExcludedSiteMasks, this->ExcludedSiteCenters, this->NbrExcludedSiteMasks);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeRealSpaceWithExclusion ComplementarySpace (ComplementaryNbrParticles, this->NbrSite,
							     this->ExcludedSiteMasks, this->ExcludedSiteCenters, this->NbrExcludedSiteMasks);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;
  FQHETorusParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, (ParticleOnTorusWithMagneticTranslations*) &ComplementarySpace, groundState, TmpDensityMatrix);
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


