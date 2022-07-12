////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of fermions on lattice with spin  and Gutzwiller          //
//             projection in real space with the Sz<->-Sz symmetry            //
//                                                                            //
//                        last modification : 18/08/2014                      //
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


// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace (int nbrFermions, int nbrSite, bool minusSzParity, unsigned long memory)
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
  this->NbrLzValue = this->LzMax;
  this->MaximumSignLookUp = 16;
  this->SzParitySign = 1.0;
  if (minusSzParity == true)
    this->SzParitySign = -1.0;
  if (this->NbrFermions > 0)
    {
      this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions);
      cout << "Intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension >= (1l << 30))
	this->HilbertSpaceDimension = 0;
      else
	this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  this->Flag.Initialize();
	  this->TargetSpace = this;
	  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
	  this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSite - 1, this->NbrSite - this->NbrFermions, 0l);
	  this->LargeHilbertSpaceDimension = this->GenerateStatesWithSzSymmetry(minusSzParity);
	  if (this->LargeHilbertSpaceDimension >= (1l << 30))
	    this->HilbertSpaceDimension = 0;
	  else
	    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
	  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
	  if (this->LargeHilbertSpaceDimension > 0l)
	    {
	      this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];  
	      this->GenerateLookUpTable(memory);
	      for (long i = 0l; i < this->HilbertSpaceDimension; ++i)
		this->GetStateSymmetry(this->StateDescription[i]);	      
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
    }
  else
    {
      this->LargeHilbertSpaceDimension = 1l;
      this->HilbertSpaceDimension = 1;
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      this->StateDescription[0] = 0x0ul;
      this->StateHighestBit[0] = 0;
      //      this->GenerateLookUpTable(memory);
    }    
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalSpin = twice the Sz projection
// nbrSite = number of sites
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace (int nbrFermions, int totalSpin, int nbrSite, bool minusSzParity, unsigned long memory)
{
  cout << "untested code, please fix it!" << endl;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (totalSpin + this->NbrFermions) >> 1;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite;
  this->NbrLzValue = this->LzMax;
  this->MaximumSignLookUp = 16;
  this->SzParitySign = 1.0;
  if (minusSzParity == true)
    this->SzParitySign = -1.0;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFermionsUp);
  cout << "Intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      if (this->LargeHilbertSpaceDimension >= (1l << 30))
	this->HilbertSpaceDimension = 0;
      else
	this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  this->Flag.Initialize();
	  this->TargetSpace = this;
	  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
	  // this should be fixed
	  this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSite - 1, this->NbrFermionsUp, 0l);
	  this->LargeHilbertSpaceDimension = this->GenerateStatesWithSzSymmetry(minusSzParity);
	  if (this->LargeHilbertSpaceDimension >= (1l << 30))
	    this->HilbertSpaceDimension = 0;
	  else
	    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
	  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
	  if (this->LargeHilbertSpaceDimension > 0l)
	    {
	      this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];  
	      this->GenerateLookUpTable(memory);
	      for (long i = 0l; i < this->HilbertSpaceDimension; ++i)
		this->GetStateSymmetry(this->StateDescription[i]);	      
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
    }
  else
    {
      this->LargeHilbertSpaceDimension = 1l;
      this->HilbertSpaceDimension = 1;
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      this->StateDescription[0] = 0x0ul;
      this->StateHighestBit[0] = 0;
    }    
}


// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace(const FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace& fermions)
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
  this->HighestBit = fermions.HighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->SzParitySign = fermions.SzParitySign;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::~FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace& FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::operator = (const FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace& fermions)
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
  this->SzParitySign = fermions.SzParitySign;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::Clone()
{
  return new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace(*this);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrHoles = number of unoccupied sites
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::GenerateStates(int nbrFermions, int currentSite, int nbrHoles, long pos)
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
  long TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, nbrHoles, pos);
  unsigned long Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, nbrHoles, pos);
   Mask = 0x1ul << ((currentSite) << 1);
   for (; pos < TmpPos; ++pos)
     this->StateDescription[pos] |= Mask;
   if (nbrHoles == 0)
     return pos;
   else
     return this->GenerateStates(nbrFermions, currentSite - 1, nbrHoles - 1, pos);
};

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrHoles = number of unoccupied sites
// nbrSpinUp = number of fermions with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::GenerateStates(int nbrFermions, int currentSite, int nbrHoles, int nbrSpinUp, long pos)
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
   
  long TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, nbrHoles, nbrSpinUp - 1, pos);
  unsigned long Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, nbrHoles, nbrSpinUp, pos);
   Mask = 0x1ul << ((currentSite) << 1);
   for (; pos < TmpPos; ++pos)
     this->StateDescription[pos] |= Mask;
   if (nbrHoles == 0)
     return pos;
   else
     return this->GenerateStates(nbrFermions, currentSite - 1, nbrHoles - 1, nbrSpinUp, pos); 
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::EvaluateHilbertSpaceDimension(int nbrFermions)
{
  BinomialCoefficients binomials(this->NbrSite);
  int NbrHoles = this->NbrSite - this->NbrFermions;
  long dimension = binomials(this->NbrSite, NbrHoles);
  for (int i = 0; i < this->NbrFermions; ++i)
    dimension *= 2l;
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

long FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::EvaluateHilbertSpaceDimension(int nbrFermions,int nbrSpinUp)
{
  BinomialCoefficients binomials(this->NbrSite);
  int NbrHoles = this->NbrSite - this->NbrFermions;
  long dimension = binomials(this->NbrSite, NbrHoles);
  BinomialCoefficients binomials1(this->NbrFermions);
  dimension *= binomials1(this->NbrFermions, this->NbrFermionsUp);
  return dimension;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  stateDescription &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  if ((stateDescription > (this->StateDescription[0] & FERMION_SPHERE_SU2_SYMMETRIC_MASK)) || 
      (stateDescription < (this->StateDescription[this->HilbertSpaceDimension - 1] & FERMION_SPHERE_SU2_SYMMETRIC_MASK)))
    return this->HilbertSpaceDimension;

  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = (this->StateDescription[PosMid] & FERMION_SPHERE_SU2_SYMMETRIC_MASK);
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
      CurrentState = (this->StateDescription[PosMid] & FERMION_SPHERE_SU2_SYMMETRIC_MASK);
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    if (((this->StateDescription[PosMin] & FERMION_SPHERE_SU2_SYMMETRIC_MASK) != stateDescription)
	&& ((this->StateDescription[PosMax] & FERMION_SPHERE_SU2_SYMMETRIC_MASK) != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}

