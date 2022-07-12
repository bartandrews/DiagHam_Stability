////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//                  in real space supporting up to 64 sites                   //
//                                                                            //
//                       class author: Cecile Repellin                        //
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
#include "HilbertSpace/FermionOnSphereWithSpin.h"
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

FermionOnLatticeWithSpinRealSpaceLong::FermionOnLatticeWithSpinRealSpaceLong ()
{
  this->NbrFermions = 0;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = 0;
  this->LzMax = this->NbrSite;
  this->NbrLzValue = this->LzMax + 1;
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
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinRealSpaceLong::FermionOnLatticeWithSpinRealSpaceLong (int nbrFermions, int nbrSite, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSite - 1, 0l);
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space, " << TmpLargeHilbertSpaceDimension << " generated states, should be " << this->LargeHilbertSpaceDimension << endl;
	}
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

// basic constructor when Sz is preserved
// 
// nbrFermions = number of fermions
// totalSpin = twice the total spin value
// nbrSite = number of sites in the x direction
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinRealSpaceLong::FermionOnLatticeWithSpinRealSpaceLong (int nbrFermions, int totalSpin, int nbrSite, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalLz = 0;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (totalSpin + this->NbrFermions) >> 1;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFermionsUp);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSite - 1, this->NbrFermionsUp, 0l);      
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space, " << TmpLargeHilbertSpaceDimension << " generated states, should be " << this->LargeHilbertSpaceDimension << endl;
	}
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

FermionOnLatticeWithSpinRealSpaceLong::FermionOnLatticeWithSpinRealSpaceLong(const FermionOnLatticeWithSpinRealSpaceLong& fermions)
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
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnLatticeWithSpinRealSpaceLong::~FermionOnLatticeWithSpinRealSpaceLong ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeWithSpinRealSpaceLong& FermionOnLatticeWithSpinRealSpaceLong::operator = (const FermionOnLatticeWithSpinRealSpaceLong& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
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
  this->SzFlag = fermions.SzFlag;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnLatticeWithSpinRealSpaceLong::Clone()
{
  return new FermionOnLatticeWithSpinRealSpaceLong(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnLatticeWithSpinRealSpaceLong::PrintState (ostream& Str, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  ULONGLONG Tmp;
  for (int i = 0; i < this->LzMax; ++i)
    {
      Tmp = (TmpState >> (i << 1)) & ((ULONGLONG) 0x3ul);
      switch (Tmp)
	{
	case ((ULONGLONG) 0x0ul):
	  Str << "0 ";
	  break;
	case ((ULONGLONG) 0x1ul):
	  Str << "d ";
	  break;
	case ((ULONGLONG) 0x2ul):
	  Str << "u ";
	  break;
	case ((ULONGLONG) 0x3ul):
	  Str << "X ";
	  break;
	}
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinRealSpaceLong::GenerateStates(int nbrFermions, int currentSite, long pos)
{
  if (nbrFermions == 0)
    {
      this->StateDescription[pos] = ((ULONGLONG) 0x0ul);	  
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0))
    return pos;
  if (nbrFermions == 1)
    {
      for (int j = currentSite; j >= 0; --j)
	{
	  this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << (j << 1);
	  ++pos;
	  this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << (j << 1);
	  ++pos;
	}
      return pos;
    }
  long TmpPos = this->GenerateStates(nbrFermions - 2, currentSite - 1, pos);
  ULONGLONG Mask = ((ULONGLONG) 0x3ul) << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, pos);
  Mask = ((ULONGLONG) 0x2ul) << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
   TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, pos);
   Mask = ((ULONGLONG) 0x1ul) << (currentSite << 1);
   for (; pos < TmpPos; ++pos)
     this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentSite - 1, pos);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrSpinUp = number of fermions with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinRealSpaceLong::GenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos)
{
  if ((nbrFermions == 0) && (nbrSpinUp == 0))
    {
      this->StateDescription[pos] = ((ULONGLONG) 0x0ul);	  
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
	      this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << (j << 1);
	      ++pos;
	    }
	}
      else
	{
	  for (int j = currentSite; j >= 0; --j)
	    {
	      this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << (j << 1);
	      ++pos;
	    }
	}
      return pos;
    }
  long TmpPos = this->GenerateStates(nbrFermions - 2, currentSite - 1, nbrSpinUp - 1, pos);
  ULONGLONG Mask = ((ULONGLONG) 0x3ul) << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, nbrSpinUp - 1, pos);
  Mask = ((ULONGLONG) 0x2ul) << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, nbrSpinUp, pos);
  Mask = ((ULONGLONG) 0x1ul) << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentSite - 1, nbrSpinUp, pos);   
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// return value = Hilbert space dimension
long FermionOnLatticeWithSpinRealSpaceLong::EvaluateHilbertSpaceDimension(int nbrFermions)
{
  BinomialCoefficients binomials(2 * this->NbrSite);
  long Dimension = binomials(2 * this->NbrSite, this->NbrFermions);
  return Dimension;
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

long FermionOnLatticeWithSpinRealSpaceLong::EvaluateHilbertSpaceDimension(int nbrFermions,int nbrSpinUp)
{
  BinomialCoefficients binomials(this->NbrSite);
  long Dimension = binomials(this->NbrSite, this->NbrFermionsUp);
  Dimension *= binomials(this->NbrSite, this->NbrFermionsDown);
  return Dimension;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = kx sector in which the density matrix has to be evaluated 
// kySector = kx sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeWithSpinRealSpaceLong::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture)
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
  FermionOnLatticeWithSpinRealSpaceLong SubsytemSpace (nbrParticleSector, this->NbrSite);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinRealSpaceLong ComplementarySpace (ComplementaryNbrParticles, this->NbrSite);
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

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in given momentum and Sz sectors.
//
// nbrParticleSector = number of particles that belong to the subsytem
// szSector  = twice the total Sz value of the subsytem
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeWithSpinRealSpaceLong::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if (szSector == 0)
        {
          HermitianMatrix TmpDensityMatrix(1, true);
          TmpDensityMatrix(0, 0) = 1.0;
          return TmpDensityMatrix;
        }
      else
        {
          HermitianMatrix TmpDensityMatrixZero;
          return TmpDensityMatrixZero;
        }
    }
  if (nbrParticleSector == this->NbrFermions)
    {
      if (this->TotalSpin == szSector)
        {
          HermitianMatrix TmpDensityMatrix(1, true);
          TmpDensityMatrix(0, 0) = 1.0;
          return TmpDensityMatrix;
        }
      else
        {
          HermitianMatrix TmpDensityMatrixZero;
          return TmpDensityMatrixZero;
        }
    }
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  int ComplementarySzSector =  this->TotalSpin - szSector;
  if (abs(ComplementarySzSector) > ComplementaryNbrParticles)
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
   }
   
  FermionOnLatticeWithSpinRealSpaceLong SubsytemSpace (nbrParticleSector, szSector, this->NbrSite, 10000000);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinRealSpaceLong ComplementarySpace (ComplementaryNbrParticles, ComplementarySzSector, this->NbrSite, 10000000);
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

// evaluate the orbital cut entanglement matrix. The entanglement matrix is only evaluated for fixed number of particles
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// keptOrbitals = array of orbitals that have to be kept, should be sorted from the smallest index to the largest index 
// nbrKeptOrbitals = array of orbitals that have to be kept
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem

ComplexMatrix FermionOnLatticeWithSpinRealSpaceLong::EvaluatePartialEntanglementMatrix (int nbrParticleSector, int nbrKeptOrbitals, int* keptOrbitals, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  if ((nbrParticleSector > (2 * nbrKeptOrbitals)) || 
      (ComplementaryNbrParticles > (2 * (this->NbrSite - nbrKeptOrbitals))))
    {
      ComplexMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
  if (nbrKeptOrbitals == 0)
    {
      if (nbrParticleSector == 0)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, this->HilbertSpaceDimension, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      TmpEntanglementMatrix[i][0] = groundState[i];
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  if (nbrKeptOrbitals == this->NbrSite)
    {
      if (nbrParticleSector == this->NbrFermions)
	{
	  ComplexMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      TmpEntanglementMatrix[0][i] = groundState[i];
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  this->KeptOrbitals = new int [nbrKeptOrbitals];
  for (int i = 0 ; i < nbrKeptOrbitals; ++i) 
    this->KeptOrbitals[i] = keptOrbitals[i];
  FermionOnLatticeWithSpinRealSpaceLong SubsytemSpace (nbrParticleSector, nbrKeptOrbitals);
  FermionOnLatticeWithSpinRealSpaceLong ComplementarySpace (ComplementaryNbrParticles, this->NbrSite - nbrKeptOrbitals);
  ComplexMatrix TmpEntanglementMatrix (SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySpace.HilbertSpaceDimension, true);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;

  long TmpEntanglementMatrixZero = this->EvaluatePartialEntanglementMatrixCore(0, ComplementarySpace.HilbertSpaceDimension, &ComplementarySpace, &SubsytemSpace, groundState, &TmpEntanglementMatrix);
//   FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpEntanglementMatrix);
//   Operation.ApplyOperation(architecture);
//   if (Operation.GetNbrNonZeroMatrixElements() > 0)	
  if (TmpEntanglementMatrixZero > 0)
     return TmpEntanglementMatrix;
   else
     {
       ComplexMatrix TmpEntanglementMatrixZero;
       return TmpEntanglementMatrixZero;
     }    
}


// evaluate the orbital cut entanglement matrix. The entanglement matrix is only evaluated for fixed number of particles and Sz
//
// nbrParticleSector = number of particles that belong to the subsytem
// szSector  = twice the total Sz value of the subsytem
// groundState = reference on the total system ground state
// keptOrbitals = array of orbitals that have to be kept, should be sorted from the smallest index to the largest index
// nbrKeptOrbitals = array of orbitals that have to be kept
// architecture = pointer to the architecture to use parallelized algorithm
// return value = entanglement matrix of the subsytem
 
ComplexMatrix FermionOnLatticeWithSpinRealSpaceLong::EvaluatePartialEntanglementMatrix (int nbrParticleSector, int szSector, int nbrKeptOrbitals, int* keptOrbitals, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  if ((nbrParticleSector > (2 * nbrKeptOrbitals)) ||
      (ComplementaryNbrParticles > (2 * (this->NbrSite - nbrKeptOrbitals))))
    {
      ComplexMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
  if (nbrKeptOrbitals == 0)
    {
      if ((nbrParticleSector == 0) && (szSector == 0))
        {
          ComplexMatrix TmpEntanglementMatrix(1, this->HilbertSpaceDimension, true);
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              TmpEntanglementMatrix[i][0] = groundState[i];
            }
          return TmpEntanglementMatrix;
        }
      else
        {
          ComplexMatrix TmpEntanglementMatrix;
          return TmpEntanglementMatrix;
        }
    }
  if (nbrKeptOrbitals == this->NbrSite)
    {
      if ((nbrParticleSector == this->NbrFermions) && (this->TotalSpin == szSector))
        {
          ComplexMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
          for (int i = 0; i < this->HilbertSpaceDimension; ++i)
            {
              TmpEntanglementMatrix[0][i] = groundState[i];
            }
          return TmpEntanglementMatrix;
        }
      else
        {
          ComplexMatrix TmpEntanglementMatrix;
          return TmpEntanglementMatrix;
        }
    }
  int ComplementarySzSector =  this->TotalSpin - szSector;
  if (abs(ComplementarySzSector) > ComplementaryNbrParticles)
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
   }
  this->KeptOrbitals = new int [nbrKeptOrbitals];
  for (int i = 0 ; i < nbrKeptOrbitals; ++i)
    this->KeptOrbitals[i] = keptOrbitals[i];
  FermionOnLatticeWithSpinRealSpaceLong SubsytemSpace (nbrParticleSector, szSector, nbrKeptOrbitals);
  FermionOnLatticeWithSpinRealSpaceLong ComplementarySpace (ComplementaryNbrParticles, ComplementarySzSector, this->NbrSite - nbrKeptOrbitals);
  ComplexMatrix TmpEntanglementMatrix (SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySpace.HilbertSpaceDimension, true);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;

  long TmpEntanglementMatrixZero = this->EvaluatePartialEntanglementMatrixCore(0, ComplementarySpace.HilbertSpaceDimension, &ComplementarySpace, &SubsytemSpace, groundState, &TmpEntanglementMatrix);
//   FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpEntanglementMatrix);
//   Operation.ApplyOperation(architecture);
//   if (Operation.GetNbrNonZeroMatrixElements() > 0)  
  if (TmpEntanglementMatrixZero > 0)
     return TmpEntanglementMatrix;
   else
     {
       ComplexMatrix TmpEntanglementMatrixZero;
       return TmpEntanglementMatrixZero;
     }    
}
  
// core part of the evaluation orbital cut entanglement matrix calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnLatticeWithSpinRealSpaceLong::EvaluatePartialEntanglementMatrixCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
									       ComplexVector& groundState, ComplexMatrix* entanglementMatrix)
{
  FermionOnLatticeWithSpinRealSpaceLong* TmpHilbertSpace = (FermionOnLatticeWithSpinRealSpaceLong*) complementaryHilbertSpace;
  FermionOnLatticeWithSpinRealSpaceLong* TmpDestinationHilbertSpace = (FermionOnLatticeWithSpinRealSpaceLong*) destinationHilbertSpace;
  long TmpNbrNonZeroElements = 0;
  int* TraceOutOrbitals = new int [this->LzMax - TmpDestinationHilbertSpace->LzMax];
  int MaxIndex = minIndex + nbrIndex;
  int TmpIndex = 0;
  for (int i = 0; i < this->LzMax; ++i)
    {
      if (SearchInArray<int>(i, this->KeptOrbitals, TmpDestinationHilbertSpace->LzMax) < 0)
	{
	  TraceOutOrbitals[TmpIndex++] = i;
	}
    }
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      ULONGLONG TmpStateCompact = TmpHilbertSpace->StateDescription[minIndex];
      ULONGLONG TmpState = ((ULONGLONG) 0x0ul);
      for (int i = 0 ; i < TmpHilbertSpace->LzMax; ++i)
	TmpState |= ((TmpStateCompact >> (2 * i)) & ((ULONGLONG) 0x3ul)) << (2 * TraceOutOrbitals[i]);
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  ULONGLONG TmpStateCompact2 = TmpDestinationHilbertSpace->StateDescription[j];
	  ULONGLONG TmpState2 = ((ULONGLONG) 0x0ul);
	  for (int i = 0 ; i < TmpDestinationHilbertSpace->LzMax; ++i)
	    TmpState2 |= ((TmpStateCompact2 >> (2 * i)) & ((ULONGLONG) 0x3ul)) << (2 * this->KeptOrbitals[i]);
	  ULONGLONG TmpState3 = TmpState | TmpState2;
	  int TmpLzMax = (this->LzMax << 1) + 1; 
	  while ((TmpState3 >> TmpLzMax) == ((ULONGLONG) 0x0ul))
	    --TmpLzMax;
	  int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
	  if ((TmpPos != this->HilbertSpaceDimension) && ((groundState[TmpPos].Re != 0.0) || (groundState[TmpPos].Im != 0.0)))
	    {
	      double Coefficient = 1.0;
	      ULONGLONG Sign = ((ULONGLONG) 0x0ul);
	      int Pos2 = (this->LzMax << 1) + 1;
	      while ((Pos2 > 0) && (TmpState2 != ((ULONGLONG) 0x0ul)))
		{
		  while (((TmpState2 >> Pos2) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
		    --Pos2;
		  TmpState3 = TmpState & ((((ULONGLONG) 0x1ul) << (Pos2 + 1)) - 1ul);
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
		  TmpState2 &= ~(((ULONGLONG) 0x1ul) << Pos2);
		  --Pos2;
		}
	      if ((Sign & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))		  
		Coefficient *= 1.0;
	      else
		Coefficient *= -1.0;
	      entanglementMatrix->AddToMatrixElement(j, minIndex, Coefficient * groundState[TmpPos]);
	      ++TmpNbrNonZeroElements;
	    }
	}
    }
  delete[] TraceOutOrbitals;
  return TmpNbrNonZeroElements;
}


// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnLatticeWithSpinRealSpaceLong::FindStateIndex(ULONGLONG stateDescription, int lzmax)
{
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  ULONGLONG CurrentState = this->StateDescription[PosMid];
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
    return PosMin;
}  

// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int FermionOnLatticeWithSpinRealSpaceLong::FindStateIndex(char* stateDescription)
{
  char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != this->LzMax)
    return -1;
  ULONGLONG TmpState = ((ULONGLONG) 0x0ul);
  int TmpNbrParticles = 0;
  for (int i = 0; i < this->LzMax; ++i)
    {
      if (TmpDescription[i][0] == 'u')
	{
	  TmpState |= ((ULONGLONG) 0x2ul) << (2 * i);
	  ++TmpNbrParticles;	  
	}
      else
	{
	  if (TmpDescription[i][0] == 'd')
	    {
	      TmpState |= ((ULONGLONG) 0x1ul) << (2 * i);
	      ++TmpNbrParticles;	  
	    }
	  else
	    {
	      if (TmpDescription[i][0] == 'X')
		{
		  TmpState |= ((ULONGLONG) 0x3ul) << (2 * i);
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
  int TmpLzMax = 2 * this->LzMax + 1;
  while (((TmpState >> TmpLzMax) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
}
