////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of fermions on lattice with spin  and Gutzwiller          //
//                          projection in real space                          //
//                                                                            //
//                       class author: Nicolas Regnault                       //
//                                                                            //
//                        last modification : 17/06/2014                      //
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
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnSphere.h"
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
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace (int nbrFermions, int nbrSite, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  if (this->NbrFermions > 0)
    {
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
	  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
	  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
	  long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSite - 1, this->NbrSite - this->NbrFermions, 0l);
	  if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	    {
	      cout << "error while generating the Hilbert space, " << TmpLargeHilbertSpaceDimension << " generated states, should be " << this->LargeHilbertSpaceDimension << endl;
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

// basic constructor when Sz is preserved
// 
// nbrFermions = number of fermions
//nbrSite = number of sites
// totalSpin = twice the total spin value
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace(int nbrFermions, int totalSpin, int nbrSite, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalLz = 0;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (totalSpin + this->NbrFermions) >> 1;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  if (this->NbrFermions > 0)
    {
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
	  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
	  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
	  long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSite - 1, this->NbrSite - this->NbrFermions, this->NbrFermionsUp, 0l);
	  if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	    {
	      cout << "error while generating the Hilbert space, " << TmpLargeHilbertSpaceDimension << " generated states, should be " << this->LargeHilbertSpaceDimension << endl;
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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace(const FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace& fermions)
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
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::~FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace& FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::operator = (const FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace& fermions)
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

AbstractHilbertSpace* FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::Clone()
{
  return new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace(*this);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrHoles = number of unoccupied sites
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::GenerateStates(int nbrFermions, int currentSite, int nbrHoles, long pos)
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

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::GenerateStates(int nbrFermions, int currentSite, int nbrHoles, int nbrSpinUp, long pos)
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

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::EvaluateHilbertSpaceDimension(int nbrFermions)
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

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::EvaluateHilbertSpaceDimension(int nbrFermions,int nbrSpinUp)
{ 
  BinomialCoefficients binomials(this->NbrSite);
  int NbrHoles = this->NbrSite - this->NbrFermions;
  long dimension = binomials(this->NbrSite, NbrHoles);
  BinomialCoefficients binomials1(this->NbrFermions);
  dimension *= binomials1(this->NbrFermions, this->NbrFermionsUp);
  return dimension;
}

// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

RealVector FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::ForgeSU2FromU1(RealVector& upState, ParticleOnSphere* upStateSpace, RealVector& downState, ParticleOnSphere* downStateSpace)
{
  return this->ForgeSU2FromU1(upState, *((FermionOnSphere*) upStateSpace), downState, *((FermionOnSphere*) downStateSpace));
}

// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

ComplexVector FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::ForgeSU2FromU1(ComplexVector& upState, ParticleOnSphere* upStateSpace, ComplexVector& downState, ParticleOnSphere* downStateSpace)
{
  return this->ForgeSU2FromU1(upState, *((FermionOnSphere*) upStateSpace), downState, *((FermionOnSphere*) downStateSpace));
}

// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

RealVector FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace)
{
  RealVector FinalState(this->HilbertSpaceDimension, true);
  if ((upStateSpace.NbrFermions + downStateSpace.NbrFermions) != this->NbrSite)
    {
      for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpUpState = upStateSpace.StateDescription[j];      
	  int TmpPos = upStateSpace.LzMax;
	  while (TmpPos > 0)
	    {
	      unsigned long Tmp = TmpUpState & (0x1ul << TmpPos);
	      TmpUpState |= Tmp << TmpPos;
	      TmpUpState ^= Tmp;
	      --TmpPos;
	    }
	  TmpUpState <<= 1;
	  double TmpComponent = upState[j];
	  int Max = 63;
	  while ((TmpUpState & (0x1ul << Max)) == 0x0ul)
	    --Max;
	  int Min = 0;
	  while ((TmpUpState & (0x1ul << Min)) == 0x0ul)
	    ++Min;
	  unsigned long TmpUpStateMask = (0x1ul << (Max + 1)) - 0x1ul;
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if ((this->StateDescription[i] & TmpUpState) == TmpUpState)
	      {	    
		unsigned long TmpUpState3 = this->StateDescription[i] & TmpUpStateMask;
		unsigned long TmpUpState2 = TmpUpState3;
#ifdef  __64_BITS__
		TmpUpState3 &= 0x5555555555555555ul;
		TmpUpState2 &= 0xaaaaaaaaaaaaaaaaul;
#else
		TmpUpState3 &= 0x55555555ul;
		TmpUpState2 &= 0xaaaaaaaaul;
#endif	    
		unsigned long Sign = 0x0;
		int Pos = this->LzMax << 1;
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
		if ((Sign & 0x1ul) == 0x0ul)
		  FinalState[i] = TmpComponent;
		else
		  FinalState[i] = -TmpComponent;
	      }
	}
      
      for (int j = 0; j < downStateSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpDownState = downStateSpace.StateDescription[j];
	  int TmpPos = downStateSpace.LzMax;
	  while (TmpPos > 0)
	    {
	      unsigned long Tmp = TmpDownState & (0x1ul << TmpPos);
	      TmpDownState |= Tmp << TmpPos;
	      TmpDownState ^= Tmp;
	      --TmpPos;
	    }
	  double TmpComponent = downState[j];
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      if ((this->StateDescription[i] & TmpDownState) == TmpDownState)
		{
		  FinalState[i] *= TmpComponent;
		}
	    }
	}
    }
  else
    {
      unsigned long TmpSpinlessMask = (0x1ul << this->NbrSite) - 0x1ul;
      for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpUpState = upStateSpace.StateDescription[j];
	  unsigned long TmpDownState = (~TmpUpState) & TmpSpinlessMask;
	  int TmpPos = downStateSpace.LzMax;
	  while ((TmpPos > 0) && ((TmpDownState & (0x1ul << TmpPos)) == 0x0ul))
	    {
	      --TmpPos;
	    }
	  int TmpSpinDownIndex = downStateSpace.FindStateIndex(TmpDownState, TmpPos);
	  if (TmpSpinDownIndex != downStateSpace.HilbertSpaceDimension)
	    {
	      TmpPos = upStateSpace.LzMax;
	      while (TmpPos > 0)
		{
		  unsigned long Tmp = TmpUpState & (0x1ul << TmpPos);
		  TmpUpState |= Tmp << TmpPos;
		  TmpUpState ^= Tmp;
		  --TmpPos;
		}
	      TmpUpState <<= 1;
	      TmpPos = downStateSpace.LzMax;
	      while (TmpPos > 0)
		{
		  unsigned long Tmp = TmpDownState & (0x1ul << TmpPos);
		  TmpDownState |= Tmp << TmpPos;
		  TmpDownState ^= Tmp;
		  --TmpPos;
		}
#ifdef  __64_BITS__
	      int Max = 63;
#else
	      int Max = 31;
#endif
	      TmpUpState |= TmpDownState;
	      while ((TmpUpState & (0x1ul << Max)) == 0x0ul)
		{
		  --Max;
		}
	      int TmpSpinfulIndex = this->FindStateIndex(TmpUpState, Max);
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
		  int Pos = this->LzMax << 1;
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
		  if ((Sign & 0x1ul) == 0x0ul)
		    FinalState[TmpSpinfulIndex] = upState[j] * downState[TmpSpinDownIndex];
		  else
		    FinalState[TmpSpinfulIndex] = -upState[j] * downState[TmpSpinDownIndex];		  
		}
	    }
	}
    }

  return FinalState;
}

// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

ComplexVector FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::ForgeSU2FromU1(ComplexVector& upState, FermionOnSphere& upStateSpace, ComplexVector& downState, FermionOnSphere& downStateSpace)
{
  ComplexVector FinalState(this->HilbertSpaceDimension, true);
  if ((upStateSpace.NbrFermions + downStateSpace.NbrFermions) != this->NbrSite)
    {
      for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpUpState = upStateSpace.StateDescription[j];
	  
	  int TmpPos = upStateSpace.LzMax;
	  while (TmpPos > 0)
	    {
	      unsigned long Tmp = TmpUpState & (0x1ul << TmpPos);
	      TmpUpState |= Tmp << TmpPos;
	      TmpUpState ^= Tmp;
	      --TmpPos;
	    }
	  TmpUpState <<= 1;
	  Complex TmpComponent = upState[j];
	  int Max = 63;
	  while ((TmpUpState & (0x1ul << Max)) == 0x0ul)
	    --Max;
	  int Min = 0;
	  while ((TmpUpState & (0x1ul << Min)) == 0x0ul)
	    ++Min;
	  unsigned long TmpUpStateMask = (0x1ul << (2 * this->LzMax + 1)) - 0x1ul;
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if ((this->StateDescription[i] & TmpUpState) == TmpUpState)
	      {	    
		unsigned long TmpUpState3 = this->StateDescription[i] & TmpUpStateMask;
		unsigned long TmpUpState2 = TmpUpState3;
#ifdef  __64_BITS__
		TmpUpState3 &= 0x5555555555555555ul;
		TmpUpState2 &= 0xaaaaaaaaaaaaaaaaul;
#else
		TmpUpState3 &= 0x55555555ul;
		TmpUpState2 &= 0xaaaaaaaaul;
#endif	    
		unsigned long Sign = 0x0;
		int Pos = this->LzMax << 1;
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
		if ((Sign & 0x1ul) == 0x0ul)
		  FinalState[i] = TmpComponent;
		else
		  FinalState[i] = -TmpComponent;
	      }
	}
      
      for (int j = 0; j < downStateSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpDownState = downStateSpace.StateDescription[j];
	  int TmpPos = downStateSpace.LzMax;
	  while (TmpPos > 0)
	    {
	      unsigned long Tmp = TmpDownState & (0x1ul << TmpPos);
	      TmpDownState |= Tmp << TmpPos;
	      TmpDownState ^= Tmp;
	      --TmpPos;
	    }
	  Complex TmpComponent = downState[j];
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      if ((this->StateDescription[i] & TmpDownState) == TmpDownState)
		{
		  FinalState[i] *= TmpComponent;
		}
	    }
	}
    }
  else
    {
      unsigned long TmpSpinlessMask = (0x1ul << this->NbrSite) - 0x1ul;
      for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpUpState = upStateSpace.StateDescription[j];
	  unsigned long TmpDownState = (~TmpUpState) & TmpSpinlessMask;
	  int TmpPos = downStateSpace.LzMax;
	  while ((TmpPos > 0) && ((TmpDownState & (0x1ul << TmpPos)) == 0x0ul))
	    {
	      --TmpPos;
	    }
	  int TmpSpinDownIndex = downStateSpace.FindStateIndex(TmpDownState, TmpPos);
	  if (TmpSpinDownIndex != downStateSpace.HilbertSpaceDimension)
	    {
	      TmpPos = upStateSpace.LzMax;
	      while (TmpPos > 0)
		{
		  unsigned long Tmp = TmpUpState & (0x1ul << TmpPos);
		  TmpUpState |= Tmp << TmpPos;
		  TmpUpState ^= Tmp;
		  --TmpPos;
		}
	      TmpUpState <<= 1;
	      TmpPos = downStateSpace.LzMax;
	      while (TmpPos > 0)
		{
		  unsigned long Tmp = TmpDownState & (0x1ul << TmpPos);
		  TmpDownState |= Tmp << TmpPos;
		  TmpDownState ^= Tmp;
		  --TmpPos;
		}
#ifdef  __64_BITS__
	      int Max = 63;
#else
	      int Max = 31;
#endif
	      TmpUpState |= TmpDownState;
	      while ((TmpUpState & (0x1ul << Max)) == 0x0ul)
		{
		  --Max;
		}
	      int TmpSpinfulIndex = this->FindStateIndex(TmpUpState, Max);
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
		  int Pos = this->LzMax << 1;
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
		  if ((Sign & 0x1ul) == 0x0ul)
		    FinalState[TmpSpinfulIndex] = upState[j] * downState[TmpSpinDownIndex];
		  else
		    FinalState[TmpSpinfulIndex] = -upState[j] * downState[TmpSpinDownIndex];		  
		}
	    }
	}
   }

  return FinalState;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = kx sector in which the density matrix has to be evaluated 
// kySector = kx sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture)
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
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace SubsytemSpace (nbrParticleSector, this->NbrSite);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite);
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

HermitianMatrix FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
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
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace SubsytemSpace (nbrParticleSector, szSector, this->NbrSite);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  FermionOnLatticeWithSpinRealSpace ComplementarySpace (ComplementaryNbrParticles, ComplementarySzSector, this->NbrSite);
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
// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
														  ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace* TmpHilbertSpace = (FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) complementaryHilbertSpace;
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace* TmpDestinationHilbertSpace = (FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace*) destinationHilbertSpace;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  Complex* TmpStateCoefficient = new Complex [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;
  BinomialCoefficients TmpBinomial (this->NbrFermions);
  double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrFermions, TmpDestinationHilbertSpace->NbrFermions));
  int MaxIndex = minIndex + nbrIndex;
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace->StateDescription[minIndex];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace->StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
	      unsigned long TmpState3 = TmpState | TmpState2;
	      int TmpLzMax = (this->LzMax << 1);
	      int TmpIndex = 0;
	      while ((TmpIndex < TmpLzMax) && (((TmpState3 >> TmpIndex) & 0x3ul) != 0x3ul))
		{
		  TmpIndex += 2;
		}

	      if (TmpIndex >= TmpLzMax)
		{
		  TmpLzMax = (this->LzMax << 1) + 1; 
		  while ((TmpState3 >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
		  if (TmpPos != this->HilbertSpaceDimension)
		    {
		      double Coefficient = TmpInvBinomial;
		      unsigned long Sign = 0x0ul;
		      int Pos2 = (TmpDestinationHilbertSpace->LzMax << 1) + 1;
		      while ((Pos2 > 0) && (TmpState2 != 0x0ul))
			{
			  while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
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
			  TmpState2 &= ~(0x1ul << Pos2);
			  --Pos2;
			}
		      if ((Sign & 0x1ul) == 0x0ul)		  
			Coefficient *= 1.0;
		      else
			Coefficient *= -1.0;
		      TmpStatePosition[Pos] = TmpPos;
		      TmpStatePosition2[Pos] = j;
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
	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]) * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
}

// evaluate the orbital cut entanglement matrix. The entanglement matrix is only evaluated for fixed number of particles
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// keptOrbitals = array of orbitals that have to be kept, should be sorted from the smallest index to the largest index 
// nbrKeptOrbitals = array of orbitals that have to be kept
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem

ComplexMatrix FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::EvaluatePartialEntanglementMatrix (int nbrParticleSector, int nbrKeptOrbitals, int* keptOrbitals, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  int ComplementaryNbrParticles = this->NbrFermions - nbrParticleSector;
  if ((nbrParticleSector > nbrKeptOrbitals) || (ComplementaryNbrParticles > (this->NbrSite - nbrKeptOrbitals)))
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
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace SubsytemSpace (nbrParticleSector, nbrKeptOrbitals);
  FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite - nbrKeptOrbitals);
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

  
// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    return this->HilbertSpaceDimension;

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
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}

