////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere including four               //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 14/10/2010                      //
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
#include "HilbertSpace/FermionOnSphereFourLandauLevels.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/ArrayTools.h"


#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <map>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::pair;
using std::map;


// default constructor
// 

FermionOnSphereFourLandauLevels::FermionOnSphereFourLandauLevels ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// memory = amount of memory granted for precalculations

FermionOnSphereFourLandauLevels::FermionOnSphereFourLandauLevels (int nbrFermions, int totalLz, int lzMax, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  this->TotalIsospin = 0;
  this->TotalEntanglement = 0;
  this->LzMax = lzMax + 6;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  this->HilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax - 6, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1);
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax - 6, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 0l);
  if (TmpHilbertSpaceDimension != this->HilbertSpaceDimension)
    {
      cout << "Mismatch in State-count (" << this->HilbertSpaceDimension << ") and State Generation (" << TmpHilbertSpaceDimension << ") in FermionOnSphereFourLandauLevels!" << endl;
      exit(1);
    }
  this->HilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
  this->GenerateLookUpTable(memory);
//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)	
//     this->PrintState(cout, i) << endl;
// #ifdef __DEBUG__
//   int UsedMemory = 0;
//   UsedMemory += this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
//   cout << "memory requested for Hilbert space = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;
//   UsedMemory = this->NbrLzValue * sizeof(int);
//   UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
//   cout << "memory requested for lookup table = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;

// #endif
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereFourLandauLevels::FermionOnSphereFourLandauLevels(const FermionOnSphereFourLandauLevels& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->TotalIsospin = fermions.TotalIsospin;
  this->TotalEntanglement = fermions.TotalEntanglement;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnSphereFourLandauLevels::~FermionOnSphereFourLandauLevels ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < 2*this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereFourLandauLevels& FermionOnSphereFourLandauLevels::operator = (const FermionOnSphereFourLandauLevels& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->TotalSpin = fermions.TotalSpin;
  this->TotalIsospin = fermions.TotalIsospin;
  this->TotalEntanglement = fermions.TotalEntanglement;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereFourLandauLevels::Clone()
{
  return new FermionOnSphereFourLandauLevels(*this);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// nbrFluxQuanta = number of flux quanta
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereFourLandauLevels::GenerateStates(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz, long pos)
{
  if ((nbrFermions < 0) || (totalLz < 0))
    return pos;
  if ((nbrFermions == 0) && (totalLz == 0))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
  if (lzMax < 0) 
    return pos;
    
  if (nbrFermions == 1) 
    {
      if (lzMax >= totalLz)
	{
	  if (((nbrFluxQuanta + 6) >= totalLz) && (totalLz >= 0))
	    {
	      this->StateDescription[pos] = 0x8ul << (totalLz << 2);
	      ++pos;
	    }
	  if (((nbrFluxQuanta + 5) >= totalLz) && (totalLz >= 1))
	    {
	      this->StateDescription[pos] = 0x4ul << (totalLz << 2);
	      ++pos;
	    }
	  if (((nbrFluxQuanta + 4) >= totalLz) && (totalLz >= 2))
	    {
	      this->StateDescription[pos] = 0x2ul << (totalLz << 2);
	      ++pos;
	    }
	  if (((nbrFluxQuanta + 3) >= totalLz) && (totalLz >= 3))
	    {
	      this->StateDescription[pos] = 0x1ul << (totalLz << 2);
	      ++pos;
	    }
	}
      return pos;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return pos;

  long TmpPos = pos;
  unsigned long Mask = 0x0ul;
  if ((lzMax <= (nbrFluxQuanta + 6)) && (lzMax >= 0))
    {
      if ((lzMax <= (nbrFluxQuanta + 5)) && (lzMax >= 1))
	{
	  if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 2))
	    {
	      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))
		{
		  TmpPos = this->GenerateStates(nbrFermions - 4, nbrFluxQuanta, lzMax - 1, totalLz - (4 * lzMax), pos);
		  Mask = 0xful << (lzMax << 2);
		  for (; pos < TmpPos; ++pos)
		    this->StateDescription[pos] |= Mask;
		}
	      TmpPos = this->GenerateStates(nbrFermions - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax), pos);
	      Mask = 0xeul << (lzMax << 2);
	      for (; pos < TmpPos; ++pos)
		this->StateDescription[pos] |= Mask;
	    }
	  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))
	    {
	      TmpPos = this->GenerateStates(nbrFermions - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax), pos);
	      Mask = 0xdul << (lzMax << 2);
	      for (; pos < TmpPos; ++pos)
		this->StateDescription[pos] |= Mask;
	    }
	  TmpPos = this->GenerateStates(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
	  Mask = 0xcul << (lzMax << 2);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= Mask;
	}
      if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 2))
	{
	  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))
	    {
	      TmpPos = this->GenerateStates(nbrFermions - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax), pos);
	      Mask = 0xbul << (lzMax << 2);
	      for (; pos < TmpPos; ++pos)
		this->StateDescription[pos] |= Mask;
	    }
	  TmpPos = this->GenerateStates(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
	  Mask = 0xaul << (lzMax << 2);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= Mask;
	}
      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))    
	{
	  TmpPos = this->GenerateStates(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
	  Mask = 0x9ul << (lzMax << 2);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= Mask;
	}  
      TmpPos = this->GenerateStates(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax, pos);
      Mask = 0x8ul << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;  
    }
  if ((lzMax <= (nbrFluxQuanta + 5)) && (lzMax >= 1))
    {
      if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 2))
	{
	  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))
	    {
	      TmpPos = this->GenerateStates(nbrFermions - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax), pos);
	      Mask = 0x7ul << (lzMax << 2);
	      for (; pos < TmpPos; ++pos)
		this->StateDescription[pos] |= Mask;
	    }
	  TmpPos = this->GenerateStates(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
	  Mask = 0x6ul << (lzMax << 2);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= Mask;
	}
      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))
	{
	  TmpPos = this->GenerateStates(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
	  Mask = 0x5ul << (lzMax << 2);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= Mask;
	}
      TmpPos = this->GenerateStates(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax, pos);
      Mask = 0x4ul << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 2))
    {
      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))
	{
	  TmpPos = this->GenerateStates(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
	  Mask = 0x3ul << (lzMax << 2);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= Mask;
	}
      TmpPos = this->GenerateStates(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax, pos);
      Mask = 0x2ul << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))    
    {
      TmpPos = this->GenerateStates(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax, pos);
      Mask = 0x1ul << (lzMax << 2);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  return this->GenerateStates(nbrFermions, nbrFluxQuanta, lzMax - 1, totalLz, pos);
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// nbrFluxQuanta = number of flux quanta
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereFourLandauLevels::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz)
{
  if ((nbrFermions < 0) || (totalLz < 0))
    return 0l;
  if ((nbrFermions == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if (nbrFermions == 1) 
    {
      long Tmp = 0l;
      if (lzMax >= totalLz)
	{
	  if (((nbrFluxQuanta + 6) >= totalLz) && (totalLz >= 0))
	    ++Tmp;
	  if (((nbrFluxQuanta + 5) >= totalLz) && (totalLz >= 1))
	    ++Tmp;
	  if (((nbrFluxQuanta + 4) >= totalLz) && (totalLz >= 2))
	    ++Tmp;
	  if (((nbrFluxQuanta + 3) >= totalLz) && (totalLz >= 3))
	    ++Tmp;
	}
      return Tmp;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;

  long Tmp = 0l;
  if ((lzMax <= (nbrFluxQuanta + 6)) && (lzMax >= 0))
    {
      if ((lzMax <= (nbrFluxQuanta + 5)) && (lzMax >= 1))
	{
	  if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 2))
	    {
	      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))    
		{
		  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 4, nbrFluxQuanta, lzMax - 1, totalLz - (4 * lzMax));
		}	  
	      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax));
	    }
	  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))    
	    {
	      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax));
	    }
	   Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
	}
      if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 2))
	{
	  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))    
	    {
	      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax));
	    }	  
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
	}
      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (nbrFluxQuanta + 5)) && (lzMax >= 1))
    {
      if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 2))
	{
	  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))    
	    {
	      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax));
	    }	  
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
	}
      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))    
	{
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
	}
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 2))    
    {
      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))    
	{
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
	}
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 3))    
    {
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
    }
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, nbrFluxQuanta, lzMax - 1, totalLz);
  return Tmp;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereFourLandauLevels::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << " | ";
  for (int i = this->NbrLzValue-1; i >=0 ; --i)
    {
      Tmp = ((TmpState >> (i << 2)) & ((unsigned long) 0xful));
      if (Tmp & 0x8ul)
	Str << "4 ";
      else
	Str << "0 ";
      if (Tmp & 0x4ul)
	Str << "3 ";
      else
	Str << "0 ";
      if (Tmp & 0x2ul)
	Str << "2 ";
      else
	Str << "0 ";
      if (Tmp & 0x1ul)
	Str << "1 ";
      else
	Str << "0 ";
      Str << "| ";
    }
  return Str;
}

// compute the projection of the product of a bosonic state and a fermionic state
//
// bosonState = real vector where the bosonic state is stored
// fermionState = real vector where the fermionic state is stored
// outputVector = real vector where the result has to be stored
// bosonSpace = pointer to the bosonic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed
// reverseFluxFlag = true if it a reverse flux attachment

void FermionOnSphereFourLandauLevels::BosonicStateTimeFermionicState(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, BosonOnSphereShort* bosonSpace, FermionOnSphere* finalSpace, int firstComponent,int nbrComponent, bool reverseFluxFlag)
{
  map<unsigned long , double> SortingMap;
  map<unsigned long , double>::iterator It;
  
  
  unsigned long* Monomial = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent+nbrComponent;
  int NbrVariable = 0;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  
  BinomialCoefficients * TmpBinomials = 0;
  
  if (reverseFluxFlag == true)
    TmpBinomials = new BinomialCoefficients(finalSpace->LzMax + this->LzMax -4);
  
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j] != 0)
	{
	  if ( reverseFluxFlag == false)
	    this->ConvertToMonomialVariable(this->StateDescription[j], Slater,NbrVariable,Variable);
	  else
	    this->ConvertToMonomialLandau(this->StateDescription[j], Slater,Variable);
	  
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if(bosonState[i] != 0)
		{
		  bosonSpace->GetMonomial(i,Monomial);
		  if ( reverseFluxFlag == false)
		    this->MonomialsTimesSlaterProjection(Slater,Monomial,Variable,NbrVariable,SortingMap,finalSpace);
		  else
		    this->MonomialsTimesSlaterProjectionReverse(Slater,Monomial,Variable,SortingMap,*TmpBinomials,finalSpace);
		  
		  for ( It = SortingMap.begin(); It != SortingMap.end(); It++)
		    {
		      int TmpLzMax = finalSpace->LzMax;
		      while ( ( ( (*It).first  >> TmpLzMax) & 0x1ul ) == 0x0ul)
			--TmpLzMax;
		      outputVector[finalSpace->FindStateIndex( (*It).first, TmpLzMax)] += bosonState[i] * fermionState[j] * (*It).second;
		    }
		    SortingMap.clear();
		}
	    }
	}
    }
}

// compute the product and the projection of a Slater determinant and a monomial 
// 
// slater = array where the slater is stored in its monomial representation
// monomial = array where the monomial is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace
// return value = number of different obtained states

void FermionOnSphereFourLandauLevels::MonomialsTimesSlaterProjection(unsigned long* slater, unsigned long* monomial, unsigned long* variable, int nbrVariable, map <unsigned long, double> & sortingMap, FermionOnSphere* finalSpace)
{
  pair <map <unsigned long, double>::iterator, bool> InsertionResult;
  
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0;
  bool Bool = true;
  double Coef = 1.0;
  long PowerIn;
  long PowerOut;
  long Numerator;
  long AlphaIn = (this->LzMax-0x2ul) * (this->LzMax-0x3ul);
  long AlphaOut = (finalSpace->LzMax+0x4ul) * (finalSpace->LzMax+0x3ul);
  long Denominator = (this->LzMax*(this->LzMax-2l)*(this->LzMax-1l)*(finalSpace->LzMax+3l)*(finalSpace->LzMax+4l)*(finalSpace->LzMax+2l));
  for (int i = 0; i < this->NbrFermions ; i++)
    State[i]=slater[i]+monomial[i];
  for(int k = 0 ; (k < nbrVariable) && (Coef != 0.0); k++)
    {
      PowerIn = (long) slater[variable[k]>>2];
      PowerOut = (long) State[variable[k]>>2];
      switch(variable[k] & 0x3ul)
	{
	case 0ul:
	  {
	    Numerator = (PowerIn-0x2ul)*(0x2ul+finalSpace->LzMax)-(PowerOut-0x2ul)*(this->LzMax-0x4ul);
	    if(Numerator == 0x0l)
	      Coef = 0.0;
	    else
	      Coef *= ((double)Numerator / (double)((this->LzMax-0x4ul)*(0x2ul+finalSpace->LzMax)));
	    break;
	  }
	case 1ul:
	  {
	    Numerator = ((AlphaOut*(PowerIn-0x1ul)*(PowerIn-0x2ul)-AlphaIn*(PowerOut-0x1ul)*(PowerOut-0x2ul))*(0x2ul+finalSpace->LzMax)-((this->LzMax-0x3ul)*(0x2ul*(PowerIn-0x1ul)-(this->LzMax-0x2ul))*AlphaOut-(finalSpace->LzMax+0x3ul)*(0x2ul*(PowerOut-0x1ul)-(finalSpace->LzMax+0x4ul))*AlphaIn)*(PowerOut-0x2ul));
	    if(Numerator == 0x0l)
	      Coef= 0.0;
	    else
	      Coef*=((double)Numerator/(double)(AlphaOut*(0x2ul+finalSpace->LzMax)));
	    break;
	  }
	case 2ul : 
	  {
	    Numerator = (-(2l + finalSpace->LzMax)*(3l + finalSpace->LzMax)*(4l + finalSpace->LzMax)*PowerIn*PowerIn*PowerIn + 3l *(3l + finalSpace->LzMax)* (4l + finalSpace->LzMax)* PowerIn*PowerIn* (6l + finalSpace->LzMax + PowerOut *(-2l + this->LzMax) - 2l* this->LzMax) + (-2l + PowerOut) *(-1l + PowerOut) *PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) *this->LzMax - (4l +finalSpace->LzMax)* PowerIn *(3l *PowerOut *(6l + finalSpace->LzMax - 3l* this->LzMax)* (-2l + this->LzMax) + 3l *PowerOut*PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) + 2l *(30l + 11l* finalSpace->LzMax + finalSpace->LzMax*finalSpace->LzMax - 3l *(6l + finalSpace->LzMax) *this->LzMax + 3l *this->LzMax*this->LzMax)));
	    
	    if(Numerator == 0x0l)
	      Coef = 0.0;
	    else
	      Coef *= ((double)Numerator/((double)(Denominator)));
	    break;
	  }
	  
	}
    }
  unsigned long Mask;
  unsigned long Sign = 0ul;
  
  if (Coef != 0.0)
    {
      for(int i = 0 ; (i < this->NbrFermions)&&(Bool); i++)
	{
	  Mask = (1ul << (State[i]-0x3ul));
	  if ( (TmpState&Mask) != 0ul)
	    Bool = false;
	  unsigned long TmpState2 = TmpState&(Mask-0x1ul);
#ifdef _64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
#endif
	  TmpState2 ^= TmpState2 >> 16;
	  TmpState2 ^= TmpState2 >> 8;
	  TmpState2 ^= TmpState2 >> 4;
	  TmpState2 ^= TmpState2 >> 2;
	  TmpState2 ^= TmpState2 >> 1;
	  Sign ^= TmpState2;
	  TmpState|=Mask;
	}
      if(Bool == true)
	{
	  if((Sign & 0x1ul) != 0ul)
	    Coef *= -1.0;

	  InsertionResult = sortingMap.insert (pair<unsigned long,double> (TmpState, Coef));
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += Coef;
	    }
	}
    }
  while (std::prev_permutation(monomial,monomial+this->NbrFermions))
    {
      Coef = 1.0;
      for(int i = 0 ; i<this->NbrFermions;i++)
	State[i] = slater[i] + monomial[i];
      
      for(int k = 0 ; (k < nbrVariable) && (Coef!= 0.0); k++)
	{
	  PowerIn = (long) slater[variable[k]>>2];
	  PowerOut = (long) State[variable[k]>>2];
	  switch(variable[k] & 0x3ul)
	    {
	    case 0ul:
	      {
		Numerator = (PowerIn-0x2ul)*(0x2ul+finalSpace->LzMax)-(PowerOut-0x2ul)*(this->LzMax-0x4ul);
		if(Numerator == 0x0l)
		  Coef = 0.0;
		else
		  Coef *= ((double)Numerator / (double)((this->LzMax-0x4ul)*(0x2ul+finalSpace->LzMax)));
		break;
	      }
	    case 1ul:
	      {
		Numerator= ((AlphaOut*(PowerIn-0x1ul)*(PowerIn-0x2ul)-AlphaIn*(PowerOut-0x1ul)*(PowerOut-0x2ul))*(0x2ul+finalSpace->LzMax)-((this->LzMax-0x3ul)*(0x2ul*(PowerIn-0x1ul)-(this->LzMax-0x2ul))*AlphaOut-(finalSpace->LzMax+0x3ul)*(0x2ul*(PowerOut-0x1ul)-(finalSpace->LzMax+0x4ul))*AlphaIn)*(PowerOut-0x2ul));
		if(Numerator == 0x0l)
		  Coef = 0.0;
		else
		  Coef *= ((double)Numerator/(double)(AlphaOut*(0x2ul+finalSpace->LzMax)));
		break;
	      }
	    case 2ul : 
	      {
		Numerator = (-(2l + finalSpace->LzMax)*(3l + finalSpace->LzMax)*(4l + finalSpace->LzMax)*PowerIn*PowerIn*PowerIn + 3l *(3l + finalSpace->LzMax)* (4l + finalSpace->LzMax)* PowerIn*PowerIn* (6l + finalSpace->LzMax + PowerOut *(-2l + this->LzMax) - 2l* this->LzMax) + (-2l + PowerOut) *(-1l + PowerOut) *PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) *this->LzMax - (4l +finalSpace->LzMax)* PowerIn *(3l *PowerOut *(6l + finalSpace->LzMax - 3l* this->LzMax)* (-2l + this->LzMax) + 3l *PowerOut*PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) + 2l *(30l + 11l* finalSpace->LzMax + finalSpace->LzMax*finalSpace->LzMax - 3l *(6l + finalSpace->LzMax) *this->LzMax + 3l *this->LzMax*this->LzMax)));
		if(Numerator == 0x0l)
		  Coef= 0.0;
		else
		  Coef *=((double)Numerator/((double)(Denominator)));
		break;
	      }
	      
	    }
	}
      
      if (Coef != 0.0)
	{
	  Bool = true;
	  TmpState = 0;
	  Sign = 0ul;
	  for (int i=0; (i < this->NbrFermions)&&(Bool);i++)
	    {
	      Mask = (1ul << (State[i]-0x3ul));
	      if((TmpState&Mask) != 0)
		Bool = false;
	      unsigned long TmpState2 = TmpState&(Mask-0x1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  if(Bool == true )
	    {
	      
	      if((Sign & 0x1ul) == 1ul)
		Coef *= -1.0;
	      
	      InsertionResult = sortingMap.insert (pair<unsigned long,double> (TmpState, Coef));
	      if (InsertionResult.second == false)
		{
		  InsertionResult.first->second += Coef;
		}
	    }
	}
    }
  delete [] State;
}

// compute the projection of the product of a fermionic state in the LLL and a fermionic state in four LL
//
// lllFermionState = real vector where the lowest Landau level fermionic state is stored
// fermionState = real vector where the two Landau level fermionic state is stored
// outputVector = real vector where the result has to be stored
// lllFermionSpace = pointer to the lowest Landau level Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereFourLandauLevels::LLLFermionicStateTimeFermionicState(RealVector& lllFermionState, RealVector& fermionState, RealVector& outputVector, FermionOnSphere* lllFermionSpace, BosonOnSphereShort* finalSpace, int firstComponent,int nbrComponent)
{
  map<unsigned long , double> SortingMap;
  map<unsigned long , double>::iterator It;
  
  unsigned long* LLLSlater = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent+nbrComponent;
  int NbrVariable = 0;
  FactorialCoefficient Coefficient;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j] != 0)
	{
	  NbrVariable = 0;
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if(lllFermionState[i] != 0)
		{
		  lllFermionSpace->GetMonomial(i, LLLSlater);
		  this->SlaterTimesSlaterProjection(Slater,LLLSlater,Variable,NbrVariable, SortingMap,finalSpace);
		  for ( It = SortingMap.begin(); It != SortingMap.end(); It++)
		    {
		      int FTmpLzMax = finalSpace->LzMax + this->NbrFermions - 1;
		      
		      while ( ( (  (*It).first  >> FTmpLzMax) & 0x1ul) == 0x0ul)
			--FTmpLzMax;
		      finalSpace->FermionToBoson( (*It).first ,FTmpLzMax,finalSpace->TemporaryState,finalSpace->TemporaryStateLzMax);
		      Coefficient.SetToOne();
		      for(int p = 0; p <= finalSpace->TemporaryStateLzMax; p++)
			{
			  Coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			}
		      outputVector[finalSpace->FermionBasis->FindStateIndex( (*It).first, FTmpLzMax)] += lllFermionState[i] * fermionState[j] *  (*It).second * Coefficient.GetIntegerValue();
		    }
		    SortingMap.clear();
		}
	    }
	}
    }
}


// compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in four Landau levels
//
// slater = array where the slater determinant in the two landau levels is stored in its monomial representation
// lllslater = array where the slater determinant in the LLL is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace

void FermionOnSphereFourLandauLevels::SlaterTimesSlaterProjection(unsigned long* slater,unsigned long* lllslater,unsigned long * variable,int nbrVariable, map <unsigned long, double> & sortingMap, BosonOnSphereShort* finalSpace)
{
  pair <map <unsigned long, double>::iterator, bool> InsertionResult;
  
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0;
  
  double Coef = 1.0;
  long PowerIn;
  long PowerOut;
  long Numerator;
  
  long AlphaIn = (this->LzMax-0x2ul) * (this->LzMax-0x3ul);
  long AlphaOut = (finalSpace->LzMax+0x4ul) * (finalSpace->LzMax+0x3ul);
  long Denominator = (this->LzMax*(this->LzMax-2l)*(this->LzMax-1l)*(finalSpace->LzMax+3l)*(finalSpace->LzMax+4l)*(finalSpace->LzMax+2l));
  
  for (int i = 0; i < this->NbrFermions ; i++)
    State[i] = slater[i] + lllslater[i];
  
  for(int k = 0 ; (k < nbrVariable) && (Coef != 0.0); k++)
    {
      PowerIn = (long) slater[variable[k]>>2];
      PowerOut = (long) State[variable[k]>>2];
      switch(variable[k] & 0x3ul)
	{
	case 0ul:
	  {
	    Numerator = (PowerIn-0x2ul)*(0x2ul+finalSpace->LzMax)-(PowerOut-0x2ul)*(this->LzMax-0x4ul);
	    if(Numerator == 0x0l)
	      Coef = 0.0;
	    else
	      Coef *= ((double)Numerator / (double)((this->LzMax-0x4ul)*(0x2ul+finalSpace->LzMax)));
	    break;
	  }
	case 1ul:
	  {
	    Numerator = ((AlphaOut*(PowerIn-0x1ul)*(PowerIn-0x2ul)-AlphaIn*(PowerOut-0x1ul)*(PowerOut-0x2ul))*(0x2ul+finalSpace->LzMax)-((this->LzMax-0x3ul)*(0x2ul*(PowerIn-0x1ul)-(this->LzMax-0x2ul))*AlphaOut-(finalSpace->LzMax+0x3ul)*(0x2ul*(PowerOut-0x1ul)-(finalSpace->LzMax+0x4ul))*AlphaIn)*(PowerOut-0x2ul));
	    
	    if(Numerator == 0x0l)
	      Coef = 0.0;
	    else
	      Coef *= ((double)Numerator/(double)(AlphaOut*(0x2ul+finalSpace->LzMax)));
	    break;
	  }
	case 2ul : 
	  {
	    Numerator = (-(2l + finalSpace->LzMax)*(3l + finalSpace->LzMax)*(4l + finalSpace->LzMax)*PowerIn*PowerIn*PowerIn + 3l *(3l + finalSpace->LzMax)* (4l + finalSpace->LzMax)* PowerIn*PowerIn* (6l + finalSpace->LzMax + PowerOut *(-2l + this->LzMax) - 2l* this->LzMax) + (-2l + PowerOut) *(-1l + PowerOut) *PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) *this->LzMax - (4l +finalSpace->LzMax)* PowerIn *(3l *PowerOut *(6l + finalSpace->LzMax - 3l* this->LzMax)* (-2l + this->LzMax) + 3l *PowerOut*PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) + 2l *(30l + 11l* finalSpace->LzMax + finalSpace->LzMax*finalSpace->LzMax - 3l *(6l + finalSpace->LzMax) *this->LzMax + 3l *this->LzMax*this->LzMax)));
	    
	    if(Numerator == 0x0l)
	      Coef = 0.0;
	    else
	      Coef *=((double)Numerator/((double)(Denominator)));
	    break;
	  }
	  
	}
    }
  
  unsigned long Mask=0ul;
  unsigned long Sign = 0ul;
  if(Coef != 0.0)
    {
      for (int i=0; (i < this->NbrFermions);i++)
	State[i] -= 3;
	
      InsertionResult = sortingMap.insert (pair<unsigned long,double> (finalSpace->ConvertFromMonomial(State), Coef));
      if (InsertionResult.second == false)
	{
	  InsertionResult.first->second += Coef;
	}
    }
  while (std::prev_permutation(lllslater, lllslater + this->NbrFermions))
    {
      Coef = 1.0;
      for(int i = 0 ; i < this->NbrFermions; i++)
	State[i] = slater[i] + lllslater[i];
      for(int k = 0 ; (k < nbrVariable) && (Coef != 0.0); k++)
	{
	  PowerIn = (long) slater[variable[k]>>2];
	  PowerOut = (long) State[variable[k]>>2];
	  switch(variable[k] & 0x3ul)
	    {
	    case 0ul:
	      {
		Numerator = (PowerIn-0x2ul)*(0x2ul+finalSpace->LzMax)-(PowerOut-0x2ul)*(this->LzMax-0x4ul);
		if(Numerator == 0x0l)
		  Coef = 0.0;
		else
		  Coef *= ((double)Numerator / (double)((this->LzMax-0x4ul)*(0x2ul+finalSpace->LzMax)));
		break;
	      }
	    case 1ul:
	      {
		Numerator = ((AlphaOut*(PowerIn-0x1ul)*(PowerIn-0x2ul)-AlphaIn*(PowerOut-0x1ul)*(PowerOut-0x2ul))*(0x2ul+finalSpace->LzMax)-((this->LzMax-0x3ul)*(0x2ul*(PowerIn-0x1ul)-(this->LzMax-0x2ul))*AlphaOut-(finalSpace->LzMax+0x3ul)*(0x2ul*(PowerOut-0x1ul)-(finalSpace->LzMax+0x4ul))*AlphaIn)*(PowerOut-0x2ul));
		
		if(Numerator == 0x0l)
		  Coef = 0.0;
		else
		  Coef  *=((double)Numerator/(double)(AlphaOut*(0x2ul+finalSpace->LzMax)));
		break;
	      }
	    case 2ul : 
	      {
		
		Numerator = (-(2l + finalSpace->LzMax)*(3l + finalSpace->LzMax)*(4l + finalSpace->LzMax)*PowerIn*PowerIn*PowerIn + 3l *(3l + finalSpace->LzMax)* (4l + finalSpace->LzMax)* PowerIn*PowerIn* (6l + finalSpace->LzMax + PowerOut *(-2l + this->LzMax) - 2l* this->LzMax) + (-2l + PowerOut) *(-1l + PowerOut) *PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) *this->LzMax - (4l +finalSpace->LzMax)* PowerIn *(3l *PowerOut *(6l + finalSpace->LzMax - 3l* this->LzMax)* (-2l + this->LzMax) + 3l *PowerOut*PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) + 2l *(30l + 11l* finalSpace->LzMax + finalSpace->LzMax*finalSpace->LzMax - 3l *(6l + finalSpace->LzMax) *this->LzMax + 3l *this->LzMax*this->LzMax)));
		if(Numerator == 0x0l)
		  Coef= 0.0;
		else
		  Coef *= ((double)Numerator/((double)(Denominator)));
		break;
	      }
	      
	    }
	}
      if( Coef != 0.0 )
	{
	  TmpState = 0ul;
	  Sign = 0ul;
	  for (int i = 0; i < this->NbrFermions ; i++)
	    {
	      State[i] -= 3;
	      Mask = (1ul << lllslater[i]);
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  SortArrayDownOrdering(State,this->NbrFermions);
	  
	  if((Sign & 0x1ul) != 0ul)
	    Coef *= -1.0;
	  
	  InsertionResult = sortingMap.insert (pair<unsigned long,double> (finalSpace->ConvertFromMonomial(State), Coef));
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += Coef;
	    }
	}
    }
  delete [] State;
}


// compute the number of particles in each Landau level
//
// state = ID of the state to handle
// LLOccupationConfiguration = array where the decomposition will be store
void  FermionOnSphereFourLandauLevels::LandauLevelOccupationNumber(int state, int * lLOccupationConfiguration)
{
  unsigned long State = this->StateDescription[state];
  for(int i = this->LzMax; i >= 0; i--)
    {
      switch ((State >> (i << 2)) & 0xful)
	{
	case 0x1ul: 
	  {
	    lLOccupationConfiguration[0]++;
	  }
	  break;
	case 0x2ul:
	  {
	    lLOccupationConfiguration[1]++;
	  }
	  break;
	case 0x3ul:
	  {
	    lLOccupationConfiguration[0]++;
	    lLOccupationConfiguration[1]++;
	  }
	  break;
	case 0x4ul: 
	  {
	    lLOccupationConfiguration[2]++;
	  }
	  break;
	case 0x5ul:
	  {
	    lLOccupationConfiguration[0]++;
	    lLOccupationConfiguration[2]++;
	  }
	  break;
	case 0x6ul:
	  {
	    lLOccupationConfiguration[1]++;
	    lLOccupationConfiguration[2]++;
	  }
	  break;
	case 0x7ul:
	  {
	    lLOccupationConfiguration[0]++;
	    lLOccupationConfiguration[1]++;
	    lLOccupationConfiguration[2]++;
	  }
	  break;
	case 0x8ul:
	  {
	    lLOccupationConfiguration[3]++;
	  }
	  break;
	case 0x9ul:
	  {
	    lLOccupationConfiguration[0]++;
	    lLOccupationConfiguration[3]++;
	  }
	  break;
	case 0xaul:
	  {
	    lLOccupationConfiguration[1]++;
	    lLOccupationConfiguration[3]++;
	  }
	  break;
	case 0xbul:
	  {
	    lLOccupationConfiguration[0]++;
	    lLOccupationConfiguration[1]++;
	    lLOccupationConfiguration[3]++;
	  }
	  break;
	case 0xcul:
	  {
	    lLOccupationConfiguration[2]++;
	    lLOccupationConfiguration[3]++;
	  }
	  break;
	case 0xdul:
	  {
	    lLOccupationConfiguration[0]++;
	    lLOccupationConfiguration[2]++;
	    lLOccupationConfiguration[3]++;
	  }
	  break;
	case 0xeul:
	  {
	    lLOccupationConfiguration[1]++;
	    lLOccupationConfiguration[2]++;
	    lLOccupationConfiguration[3]++;
	  }
	  break;
	case 0xful:
	  {
	    lLOccupationConfiguration[0]++;
	    lLOccupationConfiguration[1]++;
	    lLOccupationConfiguration[2]++;
	    lLOccupationConfiguration[3]++;
	  }
	  break;
	default : 
	  {
	    break;
	  }
	}
    }
}

// compute the projection of the product of a bosonic state and a fermionic state
//
// bosonState = real vector where the bosonic state is stored
// fermionState = real vector where the fermionic state is stored
// outputVector = real vector where the result has to be stored
// bosonSpace = pointer to the bosonic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereFourLandauLevels::BosonicStateTimeFermionicState(LongRationalVector & bosonState, LongRationalVector& fermionState, LongRationalVector& outputVector, BosonOnSphereShort* bosonSpace, FermionOnSphere* finalSpace, int firstComponent,int nbrComponent)
{
  map<unsigned long , LongRational> SortingMap;
  map<unsigned long , LongRational>::iterator It;
  
  unsigned long* Monomial = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent+nbrComponent;
  int NbrVariable = 0;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j].IsZero() == false )
	{
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater,NbrVariable,Variable);
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if(bosonState[i].IsZero() == false)
		{
		  bosonSpace->GetMonomial(i,Monomial);
		  this->MonomialsTimesSlaterProjection(Slater,Monomial,Variable,NbrVariable,SortingMap,finalSpace);
		  for ( It = SortingMap.begin(); It != SortingMap.end(); It++)
		    {
		      int TmpLzMax = finalSpace->LzMax;
		      while ( ( ( (*It).first >> TmpLzMax) & 0x1ul) == 0x0ul)
			--TmpLzMax;
		      outputVector[finalSpace->FindStateIndex( (*It).first , TmpLzMax)] += bosonState[i] * fermionState[j] *  (*It).second;
		    }
		    SortingMap.clear();
		}
	    }
	}
    }
}

// compute the product and the projection of a Slater determinant and a monomial 
// 
// slater = array where the slater is stored in its monomial representation
// monomial = array where the monomial is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace

void FermionOnSphereFourLandauLevels::MonomialsTimesSlaterProjection(unsigned long* slater, unsigned long* monomial, unsigned long* variable, int nbrVariable, map <unsigned long, LongRational> & sortingMap, FermionOnSphere* finalSpace)
{

  pair <map <unsigned long, LongRational>::iterator, bool> InsertionResult;

  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0;

  bool Bool = true;
  LongRational Coef = 1l;
  
  long PowerIn;
  long PowerOut;
  long Numerator;
  long AlphaIn = (this->LzMax-0x2ul) * (this->LzMax-0x3ul);
  long AlphaOut = (finalSpace->LzMax+0x4ul) * (finalSpace->LzMax+0x3ul);
  long Denominator = (this->LzMax*(this->LzMax-2l)*(this->LzMax-1l)*(finalSpace->LzMax+3l)*(finalSpace->LzMax+4l)*(finalSpace->LzMax+2l));
  
  for (int i = 0; i < this->NbrFermions ; i++)
    State[i]=slater[i]+monomial[i];
  
  for(int k = 0 ; (k < nbrVariable) && (Coef.IsZero() == false); k++)
    {
      PowerIn = (long) slater[variable[k]>>2];
      PowerOut = (long) State[variable[k]>>2];
      switch(variable[k] & 0x3ul)
	{
	case 0ul:
	  {
	    Numerator = (PowerIn-0x2ul)*(0x2ul+finalSpace->LzMax)-(PowerOut-0x2ul)*(this->LzMax-0x4ul);
	    if(Numerator == 0x0l)
	      Coef = 0l;
	    else
	      {
		Coef *= Numerator;
		Coef /= (this->LzMax-0x4ul)*(0x2ul+finalSpace->LzMax);
	      }
	    break;
	  }
	case 1ul:
	  {
	    Numerator = ((AlphaOut*(PowerIn-0x1ul)*(PowerIn-0x2ul)-AlphaIn*(PowerOut-0x1ul)*(PowerOut-0x2ul))*(0x2ul+finalSpace->LzMax)-((this->LzMax-0x3ul)*(0x2ul*(PowerIn-0x1ul)-(this->LzMax-0x2ul))*AlphaOut-(finalSpace->LzMax+0x3ul)*(0x2ul*(PowerOut-0x1ul)-(finalSpace->LzMax+0x4ul))*AlphaIn)*(PowerOut-0x2ul));
	    if(Numerator == 0x0l)
	      Coef= 0l;
	    else
	      {
		Coef *= Numerator;
		Coef /= AlphaOut*(0x2ul+finalSpace->LzMax);
	      }
	    break;
	  }
	case 2ul : 
	  {
	    Numerator = (-(2l + finalSpace->LzMax)*(3l + finalSpace->LzMax)*(4l + finalSpace->LzMax)*PowerIn*PowerIn*PowerIn + 3l *(3l + finalSpace->LzMax)* (4l + finalSpace->LzMax)* PowerIn*PowerIn* (6l + finalSpace->LzMax + PowerOut *(-2l + this->LzMax) - 2l* this->LzMax) + (-2l + PowerOut) *(-1l + PowerOut) *PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) *this->LzMax - (4l +finalSpace->LzMax)* PowerIn *(3l *PowerOut *(6l + finalSpace->LzMax - 3l* this->LzMax)* (-2l + this->LzMax) + 3l *PowerOut*PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) + 2l *(30l + 11l* finalSpace->LzMax + finalSpace->LzMax*finalSpace->LzMax - 3l *(6l + finalSpace->LzMax) *this->LzMax + 3l *this->LzMax*this->LzMax)));
	    
	    if(Numerator == 0x0l)
	      Coef = 0l;
	    else
	      {
		Coef *= Numerator;
		Coef /= Denominator;
	      }
	    break;
	  }
	  
	}
    }
  unsigned long Mask;
  unsigned long Sign = 0ul;
  
  if (Coef.IsZero() == false)
    {
      for(int i = 0 ; (i < this->NbrFermions)&&(Bool); i++)
	{
	  Mask = (1ul << (State[i]-0x3ul));
	  if ( (TmpState&Mask) != 0ul)
	    Bool = false;
	  unsigned long TmpState2 = TmpState&(Mask-0x1ul);
#ifdef _64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
#endif
	  TmpState2 ^= TmpState2 >> 16;
	  TmpState2 ^= TmpState2 >> 8;
	  TmpState2 ^= TmpState2 >> 4;
	  TmpState2 ^= TmpState2 >> 2;
	  TmpState2 ^= TmpState2 >> 1;
	  Sign ^= TmpState2;
	  TmpState|=Mask;
	}
      if(Bool == true)
	{
	  if((Sign & 0x1ul) != 0ul)
	    Coef *= -1l;
      
	  InsertionResult = sortingMap.insert (pair<unsigned long, LongRational> ( TmpState , Coef));
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += Coef;
	    }
	}
    }
  while (std::prev_permutation(monomial,monomial+this->NbrFermions))
    {
      Coef = 1l;
      
      for(int i = 0 ; i<this->NbrFermions;i++)
	State[i] = slater[i] + monomial[i];
      
      for(int k = 0 ; (k < nbrVariable) && (Coef.IsZero() == false); k++)
	{
	  PowerIn = (long) slater[variable[k]>>2];
	  PowerOut = (long) State[variable[k]>>2];
	  switch(variable[k] & 0x3ul)
	    {
	    case 0ul:
	      {
		Numerator = (PowerIn-0x2ul)*(0x2ul+finalSpace->LzMax)-(PowerOut-0x2ul)*(this->LzMax-0x4ul);
		if(Numerator == 0x0l)
		  Coef = 0l;
		else
		  {
		    Coef *= Numerator;
		    Coef /= (this->LzMax-0x4ul)*(0x2ul+finalSpace->LzMax);
		  }
		break;
	      }
	    case 1ul:
	      {
		Numerator = ((AlphaOut*(PowerIn-0x1ul)*(PowerIn-0x2ul)-AlphaIn*(PowerOut-0x1ul)*(PowerOut-0x2ul))*(0x2ul+finalSpace->LzMax)-((this->LzMax-0x3ul)*(0x2ul*(PowerIn-0x1ul)-(this->LzMax-0x2ul))*AlphaOut-(finalSpace->LzMax+0x3ul)*(0x2ul*(PowerOut-0x1ul)-(finalSpace->LzMax+0x4ul))*AlphaIn)*(PowerOut-0x2ul));
		if(Numerator == 0x0l)
		  Coef = 0l;
		else
		  {
		    Coef *= Numerator;
		    Coef /= AlphaOut*(0x2ul+finalSpace->LzMax);
		  }
		break;
	      }
	    case 2ul : 
	      {
		Numerator = (-(2l + finalSpace->LzMax)*(3l + finalSpace->LzMax)*(4l + finalSpace->LzMax)*PowerIn*PowerIn*PowerIn + 3l *(3l + finalSpace->LzMax)* (4l + finalSpace->LzMax)* PowerIn*PowerIn* (6l + finalSpace->LzMax + PowerOut *(-2l + this->LzMax) - 2l* this->LzMax) + (-2l + PowerOut) *(-1l + PowerOut) *PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) *this->LzMax - (4l +finalSpace->LzMax)* PowerIn *(3l *PowerOut *(6l + finalSpace->LzMax - 3l* this->LzMax)* (-2l + this->LzMax) + 3l *PowerOut*PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) + 2l *(30l + 11l* finalSpace->LzMax + finalSpace->LzMax*finalSpace->LzMax - 3l *(6l + finalSpace->LzMax) *this->LzMax + 3l *this->LzMax*this->LzMax)));
		
		if(Numerator == 0x0l)
		  Coef= 0l;
		else
		  {
		    Coef *= Numerator;
		    Coef /= Denominator;
		  }
		break;
	      }
	    }
	}
      
      if (Coef.IsZero() == false)
	{
	  Bool = true;
	  TmpState = 0;
	  Sign = 0ul;
	  for (int i= 0; (i < this->NbrFermions) && ( Bool == true ); i++)
	    {
		  Mask = (1ul << (State[i]-0x3ul));
		  if((TmpState&Mask) != 0)
		    Bool = false;
		  unsigned long TmpState2 = TmpState&(Mask-0x1ul);
#ifdef  __64_BITS__
		  TmpState2 ^= TmpState2 >> 32;
#endif
		  TmpState2 ^= TmpState2 >> 16;
		  TmpState2 ^= TmpState2 >> 8;
		  TmpState2 ^= TmpState2 >> 4;
		  TmpState2 ^= TmpState2 >> 2;
		  TmpState2 ^= TmpState2 >> 1;
		  Sign ^= TmpState2;
		  TmpState |= Mask;
		}
	      
	      if(Bool == true)
		{
		  if((Sign & 0x1ul) == 1ul)
		    Coef *= -1l;
	  
		  InsertionResult = sortingMap.insert (pair<unsigned long, LongRational> ( TmpState , Coef));
		  if (InsertionResult.second == false)
		    {
		      InsertionResult.first->second += Coef;
		    }
		}
	}
    }
  delete [] State;
}

  
// compute the projection of the product of a fermionic state in the LLL and a fermionic state in four LL
//
// lllFermionState = real vector where the lowest Landau level fermionic state is stored
// fermionState = real vector where the two Landau level fermionic state is stored
// outputVector = real vector where the result has to be stored
// lllFermionSpace = pointer to the lowest Landau level Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereFourLandauLevels::LLLFermionicStateTimeFermionicState(LongRationalVector& lllFermionState, LongRationalVector& fermionState, LongRationalVector & outputVector, FermionOnSphere* lllFermionSpace, BosonOnSphereShort* finalSpace, int firstComponent, int nbrComponent)
{ 
  map<unsigned long , LongRational> SortingMap;
  map<unsigned long , LongRational>::iterator It;
  
  unsigned long* LLLSlater = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  
  int NbrMax = firstComponent + nbrComponent;
  int NbrVariable = 0;
  
  FactorialCoefficient Coefficient;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j].IsZero() == false)
	{
	  NbrVariable = 0;
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if(lllFermionState[i].IsZero() == false)
		{
		  lllFermionSpace->GetMonomial(i, LLLSlater);
		  this->SlaterTimesSlaterProjection(Slater, LLLSlater, Variable, NbrVariable, SortingMap, finalSpace);
		 
		  for ( It = SortingMap.begin(); It != SortingMap.end(); It++)		      
		    {
		      int FTmpLzMax = finalSpace->LzMax+this->NbrFermions-1;
		      
		      while (( ( (*It).first >> FTmpLzMax) & 0x1ul) == 0x0ul)
			--FTmpLzMax;
		      finalSpace->FermionToBoson( (*It).first ,FTmpLzMax,finalSpace->TemporaryState,finalSpace->TemporaryStateLzMax);
		      Coefficient.SetToOne();
		      for(int p = 0; p <= finalSpace->TemporaryStateLzMax; p++)
			{
			  Coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			}
		      
		      outputVector[finalSpace->FermionBasis->FindStateIndex( (*It).first,FTmpLzMax)] += lllFermionState[i] * fermionState[j] * (*It).second * Coefficient.GetIntegerValue();
		    }
		    SortingMap.clear();
		}
	    }
	}
    }
}


// compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in four Landau levels
//
// slater = array where the slater determinant in the two landau levels is stored in its monomial representation
// lllslater = array where the slater determinant in the LLL is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// sortingMap = map in which the generated states and their coefficient will be stored
// finalSpace = pointer to the final HilbertSpace

void FermionOnSphereFourLandauLevels::SlaterTimesSlaterProjection(unsigned long* slater, unsigned long* lllslater, unsigned long * variable, int nbrVariable, map <unsigned long, LongRational > & sortingMap, BosonOnSphereShort* finalSpace)
{
  pair <map <unsigned long, LongRational>::iterator, bool> InsertionResult;
    
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0;
  LongRational Coef = 1l;
  
  long PowerIn;
  long PowerOut;
  long Numerator;
  long AlphaIn = (this->LzMax-0x2ul) * (this->LzMax-0x3ul);
  long AlphaOut = (finalSpace->LzMax+0x4ul) * (finalSpace->LzMax+0x3ul);
  long Denominator = (this->LzMax*(this->LzMax-2l)*(this->LzMax-1l)*(finalSpace->LzMax+3l)*(finalSpace->LzMax+4l)*(finalSpace->LzMax+2l));
  
  for (int i = 0; i < this->NbrFermions ; i++)
    State[i] = slater[i] + lllslater[i];
  
  for(int k = 0 ; (k < nbrVariable) && (Coef.IsZero() == false); k++)
    {
      PowerIn = (long) slater[variable[k]>>2];
      PowerOut = (long) State[variable[k]>>2];
      switch(variable[k] & 0x3ul)
	{
	case 0ul:
	  {
	    Numerator = (PowerIn-0x2ul)*(0x2ul+finalSpace->LzMax)-(PowerOut-0x2ul)*(this->LzMax-0x4ul);
	    if(Numerator == 0x0l)
	      Coef = 0l;
	    else
	      {
		Coef *= Numerator;
		Coef /= (this->LzMax-0x4ul)*(0x2ul+finalSpace->LzMax);
	      }
	    break;
	  }
	case 1ul:
	  {
	    Numerator = ((AlphaOut*(PowerIn-0x1ul)*(PowerIn-0x2ul)-AlphaIn*(PowerOut-0x1ul)*(PowerOut-0x2ul))*(0x2ul+finalSpace->LzMax)-((this->LzMax-0x3ul)*(0x2ul*(PowerIn-0x1ul)-(this->LzMax-0x2ul))*AlphaOut-(finalSpace->LzMax+0x3ul)*(0x2ul*(PowerOut-0x1ul)-(finalSpace->LzMax+0x4ul))*AlphaIn)*(PowerOut-0x2ul));
	    if(Numerator == 0x0l)
	      Coef = 0l;
	    else
	      {
		Coef *= Numerator;
		Coef /= AlphaOut*(0x2ul+finalSpace->LzMax);
	      }
	    break;
	  }
	case 2ul : 
	  {
	    Numerator = (-(2l + finalSpace->LzMax)*(3l + finalSpace->LzMax)*(4l + finalSpace->LzMax)*PowerIn*PowerIn*PowerIn + 3l *(3l + finalSpace->LzMax)* (4l + finalSpace->LzMax)* PowerIn*PowerIn* (6l + finalSpace->LzMax + PowerOut *(-2l + this->LzMax) - 2l* this->LzMax) + (-2l + PowerOut) *(-1l + PowerOut) *PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) *this->LzMax - (4l +finalSpace->LzMax)* PowerIn *(3l *PowerOut *(6l + finalSpace->LzMax - 3l* this->LzMax)* (-2l + this->LzMax) + 3l *PowerOut*PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) + 2l *(30l + 11l* finalSpace->LzMax + finalSpace->LzMax*finalSpace->LzMax - 3l *(6l + finalSpace->LzMax) *this->LzMax + 3l *this->LzMax*this->LzMax)));
	    
	    if(Numerator == 0x0l)
	      Coef= 0l;
	    else
	      {
		Coef *= Numerator;
		Coef /= Denominator;
	      }
	    break;
	  }
	}
    }
  
  unsigned long Mask=0ul;
  unsigned long Sign = 0ul;
  if ( Coef.IsZero() == false)
    {
      for (int i = 0; i < this->NbrFermions; i++)
	State[i] -= 3;
      
      InsertionResult = sortingMap.insert ( pair <unsigned long, LongRational> (finalSpace->ConvertFromMonomial(State), Coef));
      
      if (InsertionResult.second == false)
	{
	  InsertionResult.first->second += Coef;
	}
      
    }
  while (std::prev_permutation(lllslater, lllslater + this->NbrFermions))
    {
      Coef = 1l;
      for(int i = 0 ; i < this->NbrFermions; i++)
	State[i] = slater[i] + lllslater[i];
      for(int k = 0 ; (k < nbrVariable) && (Coef.IsZero() == false); k++)
	{
	  PowerIn = (long) slater[variable[k]>>2];
	  PowerOut = (long) State[variable[k]>>2];
	  switch(variable[k] & 0x3ul)
	    {
	    case 0ul:
	      {
		Numerator = (PowerIn-0x2ul)*(0x2ul+finalSpace->LzMax)-(PowerOut-0x2ul)*(this->LzMax-0x4ul);
		if(Numerator == 0x0l)
		  Coef = 0l;
		else
		  {
		    Coef *= Numerator;
		    Coef /= (this->LzMax-0x4ul)*(0x2ul+finalSpace->LzMax);
		  }
		break;
	      }
	    case 1ul:
	      {
		Numerator = ((AlphaOut*(PowerIn-0x1ul)*(PowerIn-0x2ul)-AlphaIn*(PowerOut-0x1ul)*(PowerOut-0x2ul))*(0x2ul+finalSpace->LzMax)-((this->LzMax-0x3ul)*(0x2ul*(PowerIn-0x1ul)-(this->LzMax-0x2ul))*AlphaOut-(finalSpace->LzMax+0x3ul)*(0x2ul*(PowerOut-0x1ul)-(finalSpace->LzMax+0x4ul))*AlphaIn)*(PowerOut-0x2ul));
		if(Numerator == 0x0l)
		  Coef = 0l;
		else
		  {
		    Coef *= Numerator;
		    Coef /= AlphaOut*(0x2ul+finalSpace->LzMax);
		  }
		break;
	      }
	    case 2ul : 
	      {
		Numerator = (-(2l + finalSpace->LzMax)*(3l + finalSpace->LzMax)*(4l + finalSpace->LzMax)*PowerIn*PowerIn*PowerIn + 3l *(3l + finalSpace->LzMax)* (4l + finalSpace->LzMax)* PowerIn*PowerIn* (6l + finalSpace->LzMax + PowerOut *(-2l + this->LzMax) - 2l* this->LzMax) + (-2l + PowerOut) *(-1l + PowerOut) *PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) *this->LzMax - (4l +finalSpace->LzMax)* PowerIn *(3l *PowerOut *(6l + finalSpace->LzMax - 3l* this->LzMax)* (-2l + this->LzMax) + 3l *PowerOut*PowerOut *(-2l + this->LzMax) *(-1l + this->LzMax) + 2l *(30l + 11l* finalSpace->LzMax + finalSpace->LzMax*finalSpace->LzMax - 3l *(6l + finalSpace->LzMax) *this->LzMax + 3l *this->LzMax*this->LzMax)));
		
		if(Numerator == 0x0l)
		  Coef= 0l;
		else
		  {
		    Coef *= Numerator;
		    Coef /= Denominator;
		  }
		break;
	      }
	    }
	}
      if ( Coef.IsZero() == false)
	{
	  TmpState = 0ul;
	  Sign = 0ul;
	  for (int i = 0; i < this->NbrFermions ; i++)
	    {
	      State[i] -= 3;
	      Mask = (1ul << lllslater[i]);
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  SortArrayDownOrdering(State,this->NbrFermions);
	  
	  if((Sign & 0x1ul) != 0ul)
	    Coef *= -1l;
	  
	  InsertionResult = sortingMap.insert ( pair <unsigned long, LongRational> (finalSpace->ConvertFromMonomial(State), Coef));
	  
	  if (InsertionResult.second == false)
	    {
	      InsertionResult.first->second += Coef;
	    }
	}
    }
  delete [] State;
}


// compute the product and the projection of a Slater determinant and a monomial with reverse flux attachment
// 
// slater = array where the slater is stored in its monomial representation
// monomial = array where the monomial is stored in its monomial representation
// landau =  array where the landau level of fermions is stored
// sortingMap = map in which the generated states and their coefficient will be stored
// binomialsCoefficient = binomials coefficient needed in the computation
// finalSpace = pointer to the final HilbertSpace

void FermionOnSphereFourLandauLevels::MonomialsTimesSlaterProjectionReverse(unsigned long* slater, unsigned long* monomial, unsigned long* landau, map <unsigned long, double> & sortingMap,BinomialCoefficients& binomialsCoefficient, FermionOnSphere* finalSpace)
{
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0;
  
  pair <map <unsigned long, double>::iterator, bool> InsertionResult;

  bool Bool = true;
  double Coef = 1.0;
  
  double TmpCoef;
  double TmpCoef2;
  int LzMaxDown = this->LzMax - 6;
  long OtherSpaceNphi =  finalSpace->LzMax + LzMaxDown;
	
  do
    {
      Coef = 1.0;
      for(int k = 0 ; (k < this->NbrFermions) && (Coef != 0.0); k++)
	{
	  TmpCoef = 0.0;
	  
	  if((monomial[k] + ((long) slater[k]) >= LzMaxDown + 3l)&&( monomial[k] + ((long) slater[k]) <= (LzMaxDown + 3l +finalSpace->LzMax)))
	    {
	      State[k] = monomial[k] + ((long) slater[k]) - LzMaxDown - 3l;
	      
	      for (int s = 0; s <= landau[k]; s++)
		{
		  if (( slater[k]  > LzMaxDown + 3l + s )||( slater[k] + landau[k] < s + 3 ))
		    {
		    }
		  else
		    {
		      TmpCoef2 = binomialsCoefficient.GetNumericalCoefficient(landau[k] , s);
		      TmpCoef2  *= binomialsCoefficient.GetNumericalCoefficient(LzMaxDown + landau[k], slater[k] + landau[k] - s - 3);
		      TmpCoef2  /= binomialsCoefficient.GetNumericalCoefficient(OtherSpaceNphi + landau[k], monomial[k] + s);
		      TmpCoef2/=(OtherSpaceNphi + 1 + landau[k]);
		      if ( (s & 0x1ul) != 0ul)
			TmpCoef -= TmpCoef2 ;
		      else
			TmpCoef += TmpCoef2;				
		    }
		}
	      TmpCoef *= binomialsCoefficient.GetNumericalCoefficient(finalSpace->LzMax,State[k]);
	      TmpCoef *= finalSpace->LzMax + 1;
	      Coef *= TmpCoef;
	    }
	  else
	    {
	      Coef = 0.0;
	    }
	}
      if(Coef != 0.0)
	{			
	  TmpState = 0ul;
	  unsigned long Mask;
	  unsigned long Sign = 0ul;
	  Bool = true;
	  for(int i = 0; ( i < this->NbrFermions ) && ( Bool == true ); i++)
	    {
	      Mask= 1ul << State[i];
	      if ( (TmpState & Mask) != 0ul)
		Bool = false;
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef _64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  if (Bool == true)
	    {
	      if ((Sign & 0x1ul) != 0ul)
		Coef *= -1.0;
	      
	      InsertionResult = sortingMap.insert (pair <unsigned long, double> (TmpState, Coef));
	      
	      if (InsertionResult.second == false)
		{
		  InsertionResult.first->second += Coef;
		}
	    }
	}
    }
  while (std::prev_permutation(monomial,monomial+this->NbrFermions));
  delete [] State;
}

