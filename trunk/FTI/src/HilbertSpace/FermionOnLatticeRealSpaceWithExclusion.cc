////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of fermions on lattice in real space              //
//                   and an exclusion rule for neighboring sites              //
//                                                                            //
//                        last modification : 17/02/2018                      //
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
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "GeneralTools/StringTools.h"

#include <math.h>
#include <cstdlib>
#include <fstream>
#include <bitset>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

FermionOnLatticeRealSpaceWithExclusion::FermionOnLatticeRealSpaceWithExclusion ()
{
  this->NbrExcludedSiteMasks = 0;
  this->ExcludedSiteMasks = 0;
  this->ExcludedSiteCenters = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// excludedSites = array of the linearized coordinates of the excluded sites around each site
// nbrExcludedSites = number of excluded sites around each site
// memory = amount of memory granted for precalculations

FermionOnLatticeRealSpaceWithExclusion::FermionOnLatticeRealSpaceWithExclusion (int nbrFermions, int nbrSite, int** excludedSites,
										int* nbrExcludedSites, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;

  this->NbrExcludedSiteMasks = this->NbrSite;
  this->ExcludedSiteMasks = new unsigned long [this->NbrExcludedSiteMasks];
  this->ExcludedSiteCenters = new unsigned long [this->NbrExcludedSiteMasks];
  for (int k = 0; k < this->NbrExcludedSiteMasks; ++k)
    {
      this->ExcludedSiteCenters[k] = 0x1ul << k;
      this->ExcludedSiteMasks[k] = 0x0ul;
      for (int l = 0; l < nbrExcludedSites[k]; ++l)
	{
	  this->ExcludedSiteMasks[k] |= 0x1ul << excludedSites[k][l];
	}
    }

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions);
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      this->Flag.Initialize();
      this->TargetSpace = this;
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
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
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

// constructor from prebuilt exclusion rules
// 
// nbrFermions = number of fermions
// nbrSite = number of sites
// excludedSiteMasks = masks used to detected excluded sites around a given position
// excludedSiteCenters = masks used to indicate the site around which a given exclusion rule is defined
// nbrExcludedSiteMasks = number of masks in ExcludedSiteMasks

FermionOnLatticeRealSpaceWithExclusion::FermionOnLatticeRealSpaceWithExclusion (int nbrFermions, int nbrSite, unsigned long* excludedSiteMasks,
										unsigned long* excludedSiteCenters, int nbrExcludedSiteMasks, 
										unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;

  this->NbrExcludedSiteMasks = nbrExcludedSiteMasks;
  this->ExcludedSiteMasks = new unsigned long [this->NbrExcludedSiteMasks];
  this->ExcludedSiteCenters = new unsigned long [this->NbrExcludedSiteMasks];
  for (int k = 0; k < this->NbrExcludedSiteMasks; ++k)
    {
      this->ExcludedSiteCenters[k] = excludedSiteCenters[k];
      this->ExcludedSiteMasks[k] = excludedSiteMasks[k];
    }

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions);
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      this->Flag.Initialize();
      this->TargetSpace = this;
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
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
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
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeRealSpaceWithExclusion::FermionOnLatticeRealSpaceWithExclusion(const FermionOnLatticeRealSpaceWithExclusion& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSite = fermions.NbrSite;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;

  this->NbrExcludedSiteMasks = fermions.NbrExcludedSiteMasks;
  this->ExcludedSiteMasks = new unsigned long [this->NbrExcludedSiteMasks];
  this->ExcludedSiteCenters = new unsigned long [this->NbrExcludedSiteMasks];
  for (int i = 0; i < this->NbrExcludedSiteMasks; ++i)
    {
      this->ExcludedSiteMasks[i] = fermions.ExcludedSiteMasks[i];
      this->ExcludedSiteCenters[i] = fermions.ExcludedSiteCenters[i];
    }

  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnLatticeRealSpaceWithExclusion::~FermionOnLatticeRealSpaceWithExclusion ()
{
  delete[] this->ExcludedSiteMasks;
  delete[] this->ExcludedSiteCenters;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeRealSpaceWithExclusion& FermionOnLatticeRealSpaceWithExclusion::operator = (const FermionOnLatticeRealSpaceWithExclusion& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->ExcludedSiteMasks;
      delete[] this->ExcludedSiteCenters;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrSite = fermions.NbrSite;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  

  this->NbrExcludedSiteMasks = fermions.NbrExcludedSiteMasks;
  this->ExcludedSiteMasks = new unsigned long [this->NbrExcludedSiteMasks];
  this->ExcludedSiteCenters = new unsigned long [this->NbrExcludedSiteMasks];
  for (int i = 0; i < this->NbrExcludedSiteMasks; ++i)
    {
      this->ExcludedSiteMasks[i] = fermions.ExcludedSiteMasks[i];
      this->ExcludedSiteCenters[i] = fermions.ExcludedSiteCenters[i];
    }

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnLatticeRealSpaceWithExclusion::Clone()
{
  return new FermionOnLatticeRealSpaceWithExclusion(*this);
}

// generate all states corresponding to the constraints
// 
// return value = Hilbert space dimension

long FermionOnLatticeRealSpaceWithExclusion::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->FermionOnLatticeRealSpace::GenerateStates(this->NbrFermions, this->NbrSite - 1, 0l);
  long TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
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
  unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->StateDescription[i] != 0x0ul)
	{
	  TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	  ++TmpLargeHilbertSpaceDimension;
	}
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  this->LargeHilbertSpaceDimension = TmpLargeHilbertSpaceDimension;

  this->StateLzMax = new int [this->LargeHilbertSpaceDimension];  
  int CurrentLzMax = this->NbrLzValue;
  while ((((this->StateDescription[0] >> CurrentLzMax) & 0x1ul) == 0x0ul)&&(CurrentLzMax>=0))
    --CurrentLzMax;
  this->StateLzMax[0] = CurrentLzMax;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while (((this->StateDescription[i] >> CurrentLzMax) & 0x1ul) == 0x0ul)
	--CurrentLzMax;
      this->StateLzMax[i] = CurrentLzMax;
    }

  return this->LargeHilbertSpaceDimension;
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = kx sector in which the density matrix has to be evaluated 
// kySector = kx sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeRealSpaceWithExclusion::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  FermionOnLatticeRealSpaceWithExclusion SubsytemSpace (nbrParticleSector, this->NbrSite, 
							this->ExcludedSiteMasks, this->ExcludedSiteCenters, this->NbrExcludedSiteMasks);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeRealSpaceWithExclusion ComplementarySpace (ComplementaryNbrParticles, this->NbrSite, 
							     this->ExcludedSiteMasks, this->ExcludedSiteCenters, this->NbrExcludedSiteMasks);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;

  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

