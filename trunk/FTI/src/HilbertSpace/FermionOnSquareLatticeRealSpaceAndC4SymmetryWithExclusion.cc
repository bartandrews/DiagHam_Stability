////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice in real space           //
//                with C4 symmetry and generic neighbor exclusion             //
//                                                                            //
//                        last modification : 25/02/2018                      //
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
#include "HilbertSpace/FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion.h"
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
#include <bitset>


using std::cout;
using std::endl;
using std::bitset;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion::FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion ()
{
  this->NbrExcludedSiteMasks = 0;
  this->ExcludedSiteMasks = 0;
  this->ExcludedSiteCenters = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSitesX = number of sites along the x (or y direction)
// symmetrySector = C4 symmetry sector (either 0, 1, 2 or 3)
  // excludedSites = array of the linearized coordinates of the excluded sites around each site
  // nbrExcludedSites = number of excluded sites around each site
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion::FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion (int nbrFermions, int nbrSitesX, 
														      int symmetrySector, int** excludedSites,
														      int* nbrExcludedSites, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->NbrSite = nbrSitesX * nbrSitesX;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum =  4;
  this->MomentumModulo =  this->MaxXMomentum;
  this->NbrSitePerUnitCell = this->NbrSite / this->MaxXMomentum;
  this->StateShift = this->MaxMomentum / this->MomentumModulo;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift / 2) % this->MomentumModulo;
  this->ComplementaryStateShift =  2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->XMomentum = symmetrySector % this->MaxXMomentum;
  this->StateXShift = this->NbrSite / this->MaxXMomentum;

  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->NbrFermionsParity = (~((unsigned long) this->NbrFermions)) & 0x1ul;
  if ((this->NbrSite & 1) == 0)
    {
      this->ComplementaryStateXShift = this->MaxMomentum - this->StateXShift;
      this->CenterSiteMask = 0x0ul;
    }
  else
    {
      this->ComplementaryStateXShift = (this->MaxMomentum - 1) - this->StateXShift;
      this->CenterSiteMask = 0x1ul << (this->NbrSite - 1);
    }
  this->ComplementaryCenterSiteMask = ~this->CenterSiteMask;
  this->CenterSitePosition = this->MaxMomentum - 1;

  this->NbrExcludedSiteMasks = 0;
  for (int i = 0; i < this->NbrSite; ++i)
    {
      if (nbrExcludedSites[i] > 0)
	{
	  this->NbrExcludedSiteMasks++;
	}
    }
  this->ExcludedSiteMasks = new unsigned long [this->NbrExcludedSiteMasks];
  this->ExcludedSiteCenters = new unsigned long [this->NbrExcludedSiteMasks];
  this->NbrExcludedSiteMasks = 0;
  for (int i = 0; i < this->NbrSite; ++i)
    {
      if (nbrExcludedSites[i] > 0)
	{
	  this->ExcludedSiteCenters[this->NbrExcludedSiteMasks] = 0x1ul << i;
	  this->ExcludedSiteMasks[this->NbrExcludedSiteMasks] = 0x0ul;
	  for (int l = 0; l < nbrExcludedSites[i]; ++l)
	    {
	      this->ExcludedSiteMasks[this->NbrExcludedSiteMasks] |= 0x1ul << excludedSites[i][l];
	    }
	  ++this->NbrExcludedSiteMasks;
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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion::FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion(const FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion& fermions)
{
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->NbrSite = fermions.NbrSite;
  this->NbrSitePerUnitCell = fermions.NbrSitePerUnitCell;
  this->MaxXMomentum = fermions.MaxXMomentum;
  this->XMomentum = fermions.XMomentum;
  this->StateXShift = fermions.StateXShift;
  this->ComplementaryStateXShift = fermions.ComplementaryStateXShift;
  this->XMomentumMask = fermions.XMomentumMask;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;
  this->CenterSiteMask = fermions.CenterSiteMask;
  this->ComplementaryCenterSiteMask = fermions.ComplementaryCenterSiteMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;

  this->NbrExcludedSiteMasks = fermions.NbrExcludedSiteMasks;
  this->ExcludedSiteMasks = new unsigned long [this->NbrExcludedSiteMasks];
  this->ExcludedSiteCenters = new unsigned long [this->NbrExcludedSiteMasks];
  for (int i = 0; i < this->NbrExcludedSiteMasks; ++i)
    {
      this->ExcludedSiteMasks[i] = fermions.ExcludedSiteMasks[i];
      this->ExcludedSiteCenters[i] = fermions.ExcludedSiteCenters[i];
    }

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

FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion::~FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion ()
{
  delete[] this->ExcludedSiteMasks;
  delete[] this->ExcludedSiteCenters;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion& FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion::operator = (const FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion& fermions)
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
    }
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->NbrSite = fermions.NbrSite;
  this->NbrSitePerUnitCell = fermions.NbrSitePerUnitCell;
 
  this->MaxXMomentum = fermions.MaxXMomentum;
  this->XMomentum = fermions.XMomentum;
  this->StateXShift = fermions.StateXShift;
  this->ComplementaryStateXShift = fermions.ComplementaryStateXShift;
  this->XMomentumMask = fermions.XMomentumMask;

  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;
  this->CenterSiteMask = fermions.CenterSiteMask;
  this->ComplementaryCenterSiteMask = fermions.ComplementaryCenterSiteMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = this->LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;

  this->NbrExcludedSiteMasks = fermions.NbrExcludedSiteMasks;
  this->ExcludedSiteMasks = new unsigned long [this->NbrExcludedSiteMasks];
  this->ExcludedSiteCenters = new unsigned long [this->NbrExcludedSiteMasks];
  for (int i = 0; i < this->NbrExcludedSiteMasks; ++i)
    {
      this->ExcludedSiteMasks[i] = fermions.ExcludedSiteMasks[i];
      this->ExcludedSiteCenters[i] = fermions.ExcludedSiteCenters[i];
    }

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

AbstractHilbertSpace* FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion::Clone()
{
  return new FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion(*this);
}


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long FermionOnSquareLatticeRealSpaceAndC4SymmetryWithExclusion::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(this->NbrFermions, this->NbrSite - 1, 0l);
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;

  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX) == this->StateDescription[i]))
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      int TmpPosition = 0;
	      unsigned long TmpState = this->StateDescription[i];
	      while ((TmpPosition < this->NbrExcludedSiteMasks) && (((TmpState & this->ExcludedSiteCenters[TmpPosition]) == 0x0ul)
								    || ((TmpState & this->ExcludedSiteMasks[TmpPosition]) == 0x0ul)))
		{
		  ++TmpPosition;
		}
	      if (TmpPosition == this->NbrExcludedSiteMasks)
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
	  unsigned long TmpState2 = TmpState;
	  for (int n = 1; n < this->MaxXMomentum; ++n)
	    {
	      TmpSign |= (this->GetSignAndApplySingleXTranslation(TmpState2) << Index) ^ ((TmpSign & (0x1ul << (Index - 1))) << 1);
	      ++Index;
	    }
	  ++TmpLargeHilbertSpaceDimension;
	}
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  return TmpLargeHilbertSpaceDimension;
}

