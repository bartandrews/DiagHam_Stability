////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//                                  in real space                             //
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
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpace.h"
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

FermionOnLatticeWithSpinSzSymmetryRealSpace::FermionOnLatticeWithSpinSzSymmetryRealSpace ()
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
  this->SzParitySign = 0.0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinSzSymmetryRealSpace::FermionOnLatticeWithSpinSzSymmetryRealSpace (int nbrFermions, int nbrSite, bool minusSzParity, unsigned long memory)
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
  this->SzParitySign = 1.0;
  if (minusSzParity == true)
    this->SzParitySign = -1.0;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions);
  cout << "Intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSite - 1, 0l);
      cout << "No symmetry" << endl;
      for (int i = 0; i < this->LargeHilbertSpaceDimension; ++i)
      {
	this->PrintState(cout, i);
	cout <<  endl;
      }
      
      
      this->LargeHilbertSpaceDimension = this->GenerateStatesWithSzSymmetry(minusSzParity);
      cout << "Sz -> -Sz symmetry" << endl;
      for (int i = 0; i < this->LargeHilbertSpaceDimension; ++i)
      {
	this->PrintState(cout, i);
	cout <<  endl;
      }
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


// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinSzSymmetryRealSpace::FermionOnLatticeWithSpinSzSymmetryRealSpace (int nbrFermions, int totalSpin, int nbrSite, bool minusSzParity, unsigned long memory)
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
  this->SzParitySign = 1.0;
  if (minusSzParity == true)
    this->SzParitySign = -1.0;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrFermionsUp);
  cout << "Intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
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


// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeWithSpinSzSymmetryRealSpace::FermionOnLatticeWithSpinSzSymmetryRealSpace(const FermionOnLatticeWithSpinSzSymmetryRealSpace& fermions)
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

FermionOnLatticeWithSpinSzSymmetryRealSpace::~FermionOnLatticeWithSpinSzSymmetryRealSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeWithSpinSzSymmetryRealSpace& FermionOnLatticeWithSpinSzSymmetryRealSpace::operator = (const FermionOnLatticeWithSpinSzSymmetryRealSpace& fermions)
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

AbstractHilbertSpace* FermionOnLatticeWithSpinSzSymmetryRealSpace::Clone()
{
  return new FermionOnLatticeWithSpinSzSymmetryRealSpace(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnLatticeWithSpinSzSymmetryRealSpace::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  for (int i = 0; i < this->LzMax; ++i)
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
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinSzSymmetryRealSpace::GenerateStates(int nbrFermions, int currentSite, long pos)
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
  long TmpPos = this->GenerateStates(nbrFermions - 2, currentSite - 1, pos);
  unsigned long Mask = 0x3ul << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, pos);
  Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
   TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, pos);
   Mask = 0x1ul << (currentSite << 1);
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

long FermionOnLatticeWithSpinSzSymmetryRealSpace::GenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos)
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
  long TmpPos = this->GenerateStates(nbrFermions - 2, currentSite - 1, nbrSpinUp - 1, pos);
  unsigned long Mask = 0x3ul << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, nbrSpinUp - 1, pos);
  Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, nbrSpinUp, pos);
  Mask = 0x1ul << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentSite - 1, nbrSpinUp, pos);   
}

// generate all states corresponding to the constraints
// 
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// return value = number of generated states

long FermionOnLatticeWithSpinSzSymmetryRealSpace::GenerateStatesWithSzSymmetry(bool minusSzParity)
{
  long TmpHilbertSpaceDimension = 0l;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->GetCanonicalState(this->StateDescription[i]) != this->StateDescription[i])
	this->StateDescription[i] = 0x0ul;
      else
	{
	  this->GetStateSymmetry(this->StateDescription[i]);
	  if ((this->StateDescription[i] & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul)
	    {
	      unsigned long TmpStateParity = this->StateDescription[i];
	      this->GetStateSingletParity(TmpStateParity);
	      if ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0) && (minusSzParity == false))
		  || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0) && (minusSzParity == true)))
		++TmpHilbertSpaceDimension;
	      else
		this->StateDescription[i] = 0x0ul;
	    }
	  else
	    {
	      ++TmpHilbertSpaceDimension;
	      this->StateDescription[i] &= ~FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT;
	    }
	}
    }
  unsigned long* TmpStateDescription = new unsigned long [TmpHilbertSpaceDimension];
  TmpHilbertSpaceDimension = 0;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    if (this->StateDescription[i] != 0x0ul)
      {
	TmpStateDescription[TmpHilbertSpaceDimension] = this->StateDescription[i];
	++TmpHilbertSpaceDimension;
      }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  return TmpHilbertSpaceDimension;
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinSzSymmetryRealSpace::EvaluateHilbertSpaceDimension(int nbrFermions)
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

long FermionOnLatticeWithSpinSzSymmetryRealSpace::EvaluateHilbertSpaceDimension(int nbrFermions,int nbrSpinUp)
{
  BinomialCoefficients binomials(this->NbrSite);
  long Dimension = binomials(this->NbrSite, this->NbrFermionsUp);
  Dimension *= binomials(this->NbrSite, this->NbrFermionsDown);
  return Dimension;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnLatticeWithSpinSzSymmetryRealSpace::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  stateDescription &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
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
    return PosMin;
}  

// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int FermionOnLatticeWithSpinSzSymmetryRealSpace::FindStateIndex(char* stateDescription)
{
  char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != this->LzMax)
    return -1;
  unsigned long TmpState = 0x0ul;
  int TmpNbrParticles = 0;
  for (int i = 0; i < this->LzMax; ++i)
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
  int TmpLzMax = 2 * this->LzMax + 1;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
}
