////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with SU(8) spin              //
//                                                                            //
//                        last modification : 11/10/2006                      //
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
#include "HilbertSpace/FermionOnSphereWithSU8Spin.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <math.h>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
//

FermionOnSphereWithSU8Spin::FermionOnSphereWithSU8Spin()
{
  this->HilbertSpaceDimension = 0;
  this->NbrFermions = 0;
  this->IncNbrFermions = 0;
  this->TotalLz = 0;
  this->LzMax = 0;
  this->NbrLzValue = 0;
  this->NbrFermions1 = 0;
  this->NbrFermions2 = 0;
  this->NbrFermions3 = 0;
  this->NbrFermions4 = 0;
  this->NbrFermions5 = 0;
  this->NbrFermions6 = 0;
  this->NbrFermions7 = 0;
  this->NbrFermions8 = 0;
  this->StateDescription = 0;
  this->StateHighestBit = 0;
  this->MaximumLookUpShift = 0;
  this->LookUpTableMemorySize = 0;
  this->LookUpTableShift = 0;
  this->LookUpTable = 0;  
  this->SignLookUpTable = 0;
  this->SignLookUpTableMask = 0;
  this->MaximumSignLookUp = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// nbrParticleSigma = array that provides the number of particles with a given internal degree of freedom
// memory = amount of memory granted for precalculations

FermionOnSphereWithSU8Spin::FermionOnSphereWithSU8Spin (int nbrFermions, int totalLz, int lzMax, int* nbrParticleSigma,
							unsigned long memory)
{
  cout << "FermionOnSphereWithSU8Spin is some heavily untested code, use at your risk and debug it" << endl;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->NbrFermions1 = nbrParticleSigma[0];
  this->NbrFermions2 = nbrParticleSigma[1];
  this->NbrFermions3 = nbrParticleSigma[2];
  this->NbrFermions4 = nbrParticleSigma[3];
  this->NbrFermions5 = nbrParticleSigma[4];
  this->NbrFermions6 = nbrParticleSigma[5];
  this->NbrFermions7 = nbrParticleSigma[6];
  this->NbrFermions8 = nbrParticleSigma[7];
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
										this->NbrFermions1, this->NbrFermions2, this->NbrFermions3, this->NbrFermions4,
										this->NbrFermions5, this->NbrFermions6, this->NbrFermions7, this->NbrFermions8);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;  
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						       this->NbrFermions1, this->NbrFermions2, this->NbrFermions3, this->NbrFermions4,
									   this->NbrFermions5, this->NbrFermions6, this->NbrFermions7, this->NbrFermions8, 0l);
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
       cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
       cout << "Mismatch in State-count and State Generation in FermionOnSphereWithSU8Spin!" << endl;
       exit(1);
     }
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
  this->GenerateLookUpTable(memory);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  
#ifdef __DEBUG__
  long UsedMemory = 0l;
  UsedMemory += this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024l)
    if (UsedMemory >= 1048576l)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
  UsedMemory = this->NbrLzValue * sizeof(int);
  UsedMemory +=((long) (this->NbrLzValue * this->LookUpTableMemorySize)) * sizeof(int);
  cout << "memory requested for lookup table = ";
  if (UsedMemory >= 1024l)
    if (UsedMemory >= 1048576l)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereWithSU8Spin::FermionOnSphereWithSU8Spin(const FermionOnSphereWithSU8Spin& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->NbrFermions1 = fermions.NbrFermions1;
  this->NbrFermions2 = fermions.NbrFermions2;
  this->NbrFermions3 = fermions.NbrFermions3;
  this->NbrFermions4 = fermions.NbrFermions4;
  this->NbrFermions5 = fermions.NbrFermions5;
  this->NbrFermions6 = fermions.NbrFermions6;
  this->NbrFermions7 = fermions.NbrFermions7;
  this->NbrFermions8 = fermions.NbrFermions8;
  this->HighestBit = fermions.HighestBit;
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

FermionOnSphereWithSU8Spin::~FermionOnSphereWithSU8Spin ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      int MaxHighestBit = this->StateHighestBit[0];
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i <= MaxHighestBit; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSU8Spin& FermionOnSphereWithSU8Spin::operator = (const FermionOnSphereWithSU8Spin& fermions)
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
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->NbrFermions1 = fermions.NbrFermions1;
  this->NbrFermions2 = fermions.NbrFermions2;
  this->NbrFermions3 = fermions.NbrFermions3;
  this->NbrFermions4 = fermions.NbrFermions4;
  this->NbrFermions5 = fermions.NbrFermions5;
  this->NbrFermions6 = fermions.NbrFermions6;
  this->NbrFermions7 = fermions.NbrFermions7;
  this->NbrFermions8 = fermions.NbrFermions8;
  this->HighestBit = fermions.HighestBit;
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

AbstractHilbertSpace* FermionOnSphereWithSU8Spin::Clone()
{
  return new FermionOnSphereWithSU8Spin(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnSphereWithSU8Spin::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnSphereWithSU8Spin::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphereWithSU8Spin::ExtractSubspace (AbstractQuantumNumber& q, 
							SubspaceSpaceConverter& converter)
{
  return 0;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereWithSU8Spin::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  unsigned long CurrentState = stateDescription >> this->LookUpTableShift[lzmax];
//  cout << hex << stateDescription << dec << " " << lzmax << " " << endl;
//  cout << this->LookUpTableShift[lzmax] << endl;
  int PosMin = this->LookUpTable[lzmax][CurrentState];
  int PosMax = this->LookUpTable[lzmax][CurrentState+ 1];
  int PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->StateDescription[PosMid];
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



// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereWithSU8Spin::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << " | ";
  for (int i = this->NbrLzValue-1; i >=0 ; --i)
    {
      Tmp = ((TmpState >> (i << 3)) & ((unsigned long) 0xfful));
      for (int j = 7; j >= 0; --j)
	{
	  if (Tmp & (0x1ul << 7))
	    {
	      Str << j << " ";
	    }
	  else
	    {
	      Str << "- ";
	    }
	}
      Str << "| ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// nbrParticles1 = number of particles with sigma=1
// nbrParticles2 = number of particles with sigma=2
// nbrParticles3 = number of particles with sigma=3
// nbrParticles4 = number of particles with sigma=4
// nbrParticles5 = number of particles with sigma=5
// nbrParticles6 = number of particles with sigma=6
// nbrParticles7 = number of particles with sigma=7
// nbrParticles8 = number of particles with sigma=8
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSU8Spin::GenerateStates(int nbrFermions, int lzMax, int totalLz, int nbrParticles1, int nbrParticles2,
						int nbrParticles3, int nbrParticles4, int nbrParticles5, int nbrParticles6,
						int nbrParticles7, int nbrParticles8, long pos)
{
  if ((nbrFermions < 0) || (totalLz < 0)
      || (nbrParticles1 < 0) || (nbrParticles2 < 0) || (nbrParticles3 < 0) || (nbrParticles4 < 0)
      || (nbrParticles5 < 0) || (nbrParticles6 < 0) || (nbrParticles7 < 0) || (nbrParticles8 < 0)
      || (nbrParticles1 > nbrFermions) || (nbrParticles2 > nbrFermions) || (nbrParticles3 > nbrFermions) || (nbrParticles4 > nbrFermions)
      || (nbrParticles5 > nbrFermions) || (nbrParticles6 > nbrFermions) || (nbrParticles7 > nbrFermions) || (nbrParticles8 > nbrFermions))
    {
      return pos;
    }
    
  if ((nbrFermions == 0) && (totalLz == 0) && (nbrParticles1 == 0) && (nbrParticles2 == 0) && (nbrParticles3 == 0) && (nbrParticles4 == 0)
      && (nbrParticles5 == 0) && (nbrParticles6 == 0) && (nbrParticles7 == 0) && (nbrParticles8 == 0))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
    
  if ((lzMax < 0) || (lzMax < nbrParticles1) || (lzMax < nbrParticles2) || (lzMax < nbrParticles3) || (lzMax < nbrParticles4) 
      || (lzMax < nbrParticles5) || (lzMax < nbrParticles6) || (lzMax < nbrParticles7) || (lzMax < nbrParticles8))
    {
      return pos;
    }
    
  if (nbrFermions == 1) 
    {
      if (lzMax >= totalLz)
	{
	  this->StateDescription[pos] = ((unsigned long) (nbrParticles1 | (nbrParticles2 << 1) | (nbrParticles3 << 2) | (nbrParticles4 << 3)
							  | (nbrParticles5 << 4) | (nbrParticles6 << 5) | (nbrParticles7 << 6)
							  | (nbrParticles8 << 7)) << (totalLz << 3));
	  return (pos + 1l);
	}
      else
	return pos;
    }

  if ((lzMax == 0)  && (totalLz != 0))
    return pos;


  long TmpPos;
  unsigned long Mask;
  int TmpNbrParticles;
  
  for (int i = 255; i >= 0; --i)
    {
      TmpNbrParticles = ((i & 1) + ((i & 2) >> 1) + ((i & 4) >> 2) + ((i & 8) >> 3)
			 + ((i & 16) >> 4) + ((i & 32) >> 5) + ((i & 64) >> 6) + ((i & 128) >> 7));
      TmpPos = this->GenerateStates(nbrFermions - TmpNbrParticles, lzMax - 1, totalLz - (TmpNbrParticles * lzMax),
				    nbrParticles1 - (i & 1), nbrParticles2 - ((i & 2) >> 1),
				    nbrParticles3 - ((i & 4) >> 2), nbrParticles4 - ((i & 8) >> 3),
				    nbrParticles5 - ((i & 16) >> 4), nbrParticles6 - ((i & 32) >> 5),
				    nbrParticles7 - ((i & 64) >> 6), nbrParticles8 - ((i & 128) >> 7), pos);
      Mask = ((unsigned long) i) << (lzMax << 3);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
   }
  return pos;
};


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereWithSU8Spin::GenerateLookUpTable(unsigned long memory)
{
  // get every highest bit poisition
  unsigned long TmpPosition = this->StateDescription[0];
#ifdef __64_BITS__
  int CurrentHighestBit = 63;
#else
  int CurrentHighestBit = 31;
#endif
  while ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
    --CurrentHighestBit;  
  int MaxHighestBit = CurrentHighestBit;
  this->StateHighestBit[0] = CurrentHighestBit;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      TmpPosition = this->StateDescription[i];
      while ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
	--CurrentHighestBit;  
      this->StateHighestBit[i] = CurrentHighestBit;
   }

  // evaluate look-up table size
  memory /= (sizeof(int*) * MaxHighestBit);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > MaxHighestBit)
    this->MaximumLookUpShift = MaxHighestBit;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [MaxHighestBit + 1];
  this->LookUpTableShift = new int [MaxHighestBit + 1];
  for (int i = 0; i <= MaxHighestBit; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];

  CurrentHighestBit = this->StateHighestBit[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
  if (CurrentHighestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentHighestBit] = 0;
  else
    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentHighestBit];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentHighestBit != this->StateHighestBit[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
 	  CurrentHighestBit = this->StateHighestBit[i];
	  TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
	  if (CurrentHighestBit < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentHighestBit] = 0;
	  else
	    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentHighestBit];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;

  // look-up tables for evaluating sign when applying creation/annihilation operators
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      Count = 0;
      TmpNbr = j;
      while (TmpNbr != 0)
	{
	  if (TmpNbr & 0x1)
	    ++Count;
	  TmpNbr >>= 1;
	}
      if (Count & 1)
	this->SignLookUpTable[j] = -1.0;
      else
	this->SignLookUpTable[j] = 1.0;
    }
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#endif
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// nbrParticles1 = number of particles with sigma=1
// nbrParticles2 = number of particles with sigma=2
// nbrParticles3 = number of particles with sigma=3
// nbrParticles4 = number of particles with sigma=4
// nbrParticles5 = number of particles with sigma=5
// nbrParticles6 = number of particles with sigma=6
// nbrParticles7 = number of particles with sigma=7
// nbrParticles8 = number of particles with sigma=8
// return value = Hilbert space dimension

long FermionOnSphereWithSU8Spin::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int nbrParticles1, int nbrParticles2,
								      int nbrParticles3, int nbrParticles4, int nbrParticles5, int nbrParticles6,
								      int nbrParticles7, int nbrParticles8)
{
  if ((nbrFermions < 0) || (totalLz < 0)
      || (nbrParticles1 < 0) || (nbrParticles2 < 0) || (nbrParticles3 < 0) || (nbrParticles4 < 0)
      || (nbrParticles5 < 0) || (nbrParticles6 < 0) || (nbrParticles7 < 0) || (nbrParticles8 < 0)
      || (nbrParticles1 > nbrFermions) || (nbrParticles2 > nbrFermions) || (nbrParticles3 > nbrFermions) || (nbrParticles4 > nbrFermions)
      || (nbrParticles5 > nbrFermions) || (nbrParticles6 > nbrFermions) || (nbrParticles7 > nbrFermions) || (nbrParticles8 > nbrFermions))
    return 0l;
  if ((lzMax < 0) || (lzMax < nbrParticles1) || (lzMax < nbrParticles2) || (lzMax < nbrParticles3) || (lzMax < nbrParticles4) 
      || (lzMax < nbrParticles5) || (lzMax < nbrParticles6) || (lzMax < nbrParticles7) || (lzMax < nbrParticles8))
    return 0l;
    
  if ((nbrFermions == 0) && (totalLz == 0) && (nbrParticles1 == 0) && (nbrParticles2 == 0) && (nbrParticles3 == 0) && (nbrParticles4 == 0)
      && (nbrParticles5 == 0) && (nbrParticles6 == 0) && (nbrParticles7 == 0) && (nbrParticles8 == 0))
    return 1l;
  if (nbrFermions == 1) 
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }

  if ((lzMax == 0)  && (totalLz != 0))
    return 0l;

  unsigned long Tmp = 0l;
  for (int i = 255; i >= 0; --i)
    {
      int TmpNbrParticles = ((i & 1) + ((i & 2) >> 1) + ((i & 4) >> 2) + ((i & 8) >> 3)
			     + ((i & 16) >> 4) + ((i & 32) >> 5) + ((i & 64) >> 6) + ((i & 128) >> 7));
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - TmpNbrParticles, lzMax - 1, totalLz - (TmpNbrParticles * lzMax),
							nbrParticles1 - (i & 1), nbrParticles2 - ((i & 2) >> 1),
							nbrParticles3 - ((i & 4) >> 2), nbrParticles4 - ((i & 8) >> 3),
							nbrParticles3 - ((i & 16) >> 4), nbrParticles4 - ((i & 32) >> 5),
							nbrParticles3 - ((i & 64) >> 6), nbrParticles4 - ((i & 128) >> 7));
    }
  return Tmp;
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphereWithSU8Spin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
							    int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
  ComplexMatrix Slatter(this->NbrFermions, this->NbrFermions);
  ComplexMatrix Functions(this->LzMax + 1, this->NbrFermions);
  RealVector TmpCoordinates(2);
  int* Indices = new int [this->NbrFermions];
  int Pos;
  int Lz;
  for (int j = 0; j < this->NbrFermions; ++j)
    {
      TmpCoordinates[0] = position[j << 1];
      TmpCoordinates[1] = position[1 + (j << 1)];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	  Functions[j].Re(i) = Tmp.Re;
	  Functions[j].Im(i) = Tmp.Im;
	}
    }
  double Factor = 1.0;
  for (int i = 2; i <= this->NbrFermions; ++i)
    Factor *= (double) i;
  Factor = 1.0 / sqrt(Factor);
  unsigned long TmpStateDescription;
  int LastComponent = firstComponent + nbrComponent;
  for (int k = firstComponent; k < LastComponent; ++k)
    {
      Pos = 0;
      Lz = 0;
      TmpStateDescription = this->StateDescription[k];
      while (Pos < this->NbrFermions)
	{
	  if ((TmpStateDescription & 0x3l) != 0x0l)
	    {
	      Indices[Pos] = Lz;
	      ++Pos;
	    }
	  ++Lz;
	  TmpStateDescription >>= 2;
	}
      for (int i = 0; i < this->NbrFermions; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < this->NbrFermions; ++j)
	    {
	      Slatter[i].Re(j) = TmpColum2.Re(Indices[j]);
	      Slatter[i].Im(j) = TmpColum2.Im(Indices[j]);
	    }
	}
      Complex SlatterDet = Slatter.Determinant();
      Value += SlatterDet * (state[k] * Factor);
    }
  delete[] Indices;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereWithSU8Spin::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}

// create a U(1) state from an SU(8) state
//
// state = vector describing the SU(8) state
// u1Space = reference on the Hilbert space associated to the U(1) state
// return value = resulting U(1) state

RealVector FermionOnSphereWithSU8Spin::ForgeU1FromSU8(RealVector& state, FermionOnSphere& u1Space)
{
  RealVector FinalState(u1Space.GetHilbertSpaceDimension(), true);
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      unsigned long TmpState = this->StateDescription[j];
      unsigned long TmpState2 = TmpState; 
      int TmpPos = this->LzMax << 3;
      while (TmpPos >=0)
	{
	  unsigned long  TmpNbrParticles = TmpState2 & 0x1ul;
	  TmpState2 >>= 1;
	  TmpNbrParticles += TmpState2 & 0x1ul;
	  TmpState2 >>= 1;
	  TmpNbrParticles += TmpState2 & 0x1ul;
	  TmpState2 >>= 1;
	  TmpNbrParticles += TmpState2 & 0x1ul;
	  TmpState2 >>= 1;
	  TmpState2 >>= 1;
	  TmpNbrParticles += TmpState2 & 0x1ul;
	  TmpState2 >>= 1;
	  TmpNbrParticles += TmpState2 & 0x1ul;
	  TmpState2 >>= 1;
	  TmpNbrParticles += TmpState2 & 0x1ul;
	  TmpState2 >>= 1;
	  if (TmpNbrParticles > 0x1ul)
	    TmpPos = 1;
	  TmpPos -= 8;
	}
      if (TmpPos != -7)
	{ 
	  TmpPos = 0;
	  TmpState2 = 0x0ul; 
	  while (TmpPos <= this->LzMax)
	    {
	      TmpState2 |= ((TmpState & 0x1ul) | ((TmpState & 0x2ul) >> 1) | ((TmpState & 0x4ul) >> 2) | ((TmpState & 0x8ul) >> 3)
			    | ((TmpState & 0x10ul) >> 4) | ((TmpState & 0x20ul) >> 5) | ((TmpState & 0x40ul) >> 6) | ((TmpState & 0x80ul) >> 7)) << TmpPos;
	      TmpState >>= 8;
	      ++TmpPos;
	    }
	  while ((TmpState2 >> TmpPos) == 0x0ul)
	    --TmpPos;
	  FinalState[u1Space.FindStateIndex(TmpState2, TmpPos)] += state[j];
	}
    }
  FinalState /= FinalState.Norm();
  return FinalState;  
}

// create a SU(2) state from an SU(8) state (fusing same spin values,i.e symmetrizing over the isospin)
//
// state = vector describing the SU(8) state
// su2Space = reference on the Hilbert space associated to the SU(2) state
// return value = resulting SU(2) state

RealVector FermionOnSphereWithSU8Spin::ForgeSU2FromSU8(RealVector& state, FermionOnSphereWithSpin& su2Space)
{
  RealVector FinalState(su2Space.GetHilbertSpaceDimension(), true);
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      unsigned long TmpState = this->StateDescription[j];
      unsigned long TmpState2 = TmpState; 
      int TmpPos = this->LzMax << 3;
      while (TmpPos >=0)
	{
	  if ((((TmpState2 >> 3) & (TmpState2 >> 2) & 0x1ul) != 0x0ul) ||
	      (((TmpState2 >> 1) & TmpState2 & 0x1ul) != 0x0ul))
	    TmpPos = 1;
	  TmpState2 >>= 8;
	  TmpPos -= 8;
	}
      if (TmpPos != -7)
	{ 
	  TmpPos = 0;
	  TmpState2 = 0x0ul; 
	  int TmpLzMax = this->LzMax << 1;
	  while (TmpPos <= TmpLzMax)
	    {
	      TmpState2 |= ((TmpState | (TmpState >> 1)) & 0x1ul) << TmpPos;
	      ++TmpPos;
	      TmpState2 |= (((TmpState >> 2) | (TmpState >> 3)) & 0x1ul) << TmpPos;
	      TmpState >>= 4;
	      ++TmpPos;
	    }
	  while ((TmpState2 >> TmpPos) == 0x0ul)
	    --TmpPos;
	  FinalState[su2Space.FindStateIndex(TmpState2, TmpPos)] += state[j];
	}
    }
  FinalState /= FinalState.Norm();
  return FinalState;  
}

// convert a state from one SU(8) basis to another, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void FermionOnSphereWithSU8Spin::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, long firstComponent, long nbrComponents)
{
  int* TmpMomentumIndices = new int [this->NbrFermions];
  int* TmpSU8Indices = new int [this->NbrFermions];
  int* TmpSU8Indices2 = new int [this->NbrFermions];
  targetState.ClearVector();
  long LastComponent = firstComponent + nbrComponents;
  if (nbrComponents == 0)
    LastComponent = this->LargeHilbertSpaceDimension;
  for (long i = firstComponent; i < LastComponent; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      unsigned long Tmp;
      int TmpLzMax = this->NbrLzValue << 3;
      int TmpIndex = 0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  Tmp = (TmpState >> (j << 3)) & 0xfful;;
	  if ((Tmp & 0x80ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 7;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x40ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 6;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x20ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 5;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x10ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 4;
	      ++TmpIndex;
	    }	  
	  if ((Tmp & 0x8ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 3;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x4ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 2;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x2ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 1;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x1ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 0;
	      ++TmpIndex;
	    }	  
	}
      this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpMomentumIndices, TmpSU8Indices, TmpSU8Indices2, oneBodyBasis);
    }
  delete[] TmpMomentumIndices;
  delete[] TmpSU8Indices;
  delete[] TmpSU8Indices2;
}

// compute the transformation matrix from one SU(8) basis to another, transforming the one body basis in each momentum sector
//
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// return value = transformation matrix

ComplexMatrix FermionOnSphereWithSU8Spin::TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis)
{
  int* TmpMomentumIndices = new int [this->NbrFermions];
  int* TmpSU8Indices = new int [this->NbrFermions];
  int* TmpSU8Indices2 = new int [this->NbrFermions];
  ComplexMatrix TmpMatrix(this->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      unsigned long Tmp;
      int TmpLzMax = this->NbrLzValue << 3;
      int TmpIndex = 0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  Tmp = (TmpState >> (j << 3)) & 0xfful;;
	  if ((Tmp & 0x80ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 7;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x40ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 6;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x20ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 5;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x10ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 4;
	      ++TmpIndex;
	    }	  
	  if ((Tmp & 0x8ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 3;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x4ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 2;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x2ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 1;
	      ++TmpIndex;
	    }
	  if ((Tmp & 0x1ul) != 0x0ul)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU8Indices[TmpIndex] = 0;
	      ++TmpIndex;
	    }	  
	}
      this->TransformOneBodyBasisRecursive(TmpMatrix[i], 1.0, 0, TmpMomentumIndices, TmpSU8Indices, TmpSU8Indices2, oneBodyBasis);
    }
  delete[] TmpMomentumIndices;
  delete[] TmpSU8Indices;
  delete[] TmpSU8Indices2;
  return TmpMatrix;
}

// recursive part of the convertion from a state from one SU(8) basis to another, transforming the one body basis in each momentum sector
//
// targetState = vector where the transformed state has to be stored
// coefficient = current coefficient to assign
// position = current particle consider in the n-body state
// momentumIndices = array that gives the momentum partition of the initial n-body state
// initialSU8Indices = array that gives the spin dressing the initial n-body state
// currentSU8Indices = array that gives the spin dressing the current transformed n-body state
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector

void FermionOnSphereWithSU8Spin::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
								int position, int* momentumIndices, int* initialSU8Indices, int* currentSU8Indices, ComplexMatrix* oneBodyBasis) 
{
//   cout << position << " : " << endl;
//   for (int i = 0; i < position; ++i)
//     cout << currentSU8Indices[i] << " ";
//   cout << endl;
  if (position == this->NbrFermions)
    {
      unsigned long TmpState = 0x0ul;
      unsigned long TmpState2;
      unsigned long Mask = 0x0ul;
      unsigned long MaskSign = 0x0ul;
      for (int i = 0; i < this->NbrFermions; ++i)
	{
	  Mask = 0x1ul << ((momentumIndices[i] << 3) + currentSU8Indices[i]);
	  if ((TmpState & Mask) != 0x0ul)
	    return;
	  TmpState2 = TmpState & (Mask - 0x1ul);
#ifdef __64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
#endif
	  TmpState2 ^= (TmpState2 >> 16);
	  TmpState2 ^= (TmpState2 >> 8);
	  TmpState2 ^= (TmpState2 >> 4);
	  TmpState2 ^= (TmpState2 >> 2);
	  MaskSign ^= (TmpState2 ^ (TmpState2 >> 1)) & 0x1ul;
	  TmpState |= Mask;
	}
      int TmpLzMax = this->NbrLzValue << 3;
      while ((TmpState >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      int Index = this->FindStateIndex(TmpState, TmpLzMax);
      if (Index < this->HilbertSpaceDimension)
	{
	  if (MaskSign == 0ul)
	    {
	      targetState[Index] += coefficient;
	    }
	  else
	    {
	      targetState[Index] -= coefficient;
	    }
	}
      return;      
    }
  else
    {
      for (int i = 0; i < 8; ++i)
	{
	  currentSU8Indices[position] = i;
	  this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][7 - initialSU8Indices[position]][7 - i]), position + 1, momentumIndices, initialSU8Indices, currentSU8Indices, oneBodyBasis);
	}
    }
}

// compute the projection matrix from the SU(8) Hilbert space to an SU(2) Hilbert space
// 
// targetSpace = pointer to the SU(2) Hilbert space
// spinUp = index of the component that has to be consider as a spin up
// spinDown = index of the component that has to be consider as a spin down
// return value = projection matrix

ComplexMatrix FermionOnSphereWithSU8Spin::TransformationMatrixSU8ToSU2(ParticleOnSphereWithSpin* targetSpace, int spinUp, int spinDown)
{
  FermionOnSphereWithSpin* TmpTargetSpace = (FermionOnSphereWithSpin*) targetSpace;
  ComplexMatrix TmpMatrix (TmpTargetSpace->HilbertSpaceDimension, this->HilbertSpaceDimension, true) ;
  unsigned long MaskUp = 0x1ul << spinUp;
  unsigned long MaskDown = 0x1ul << spinDown;
  unsigned long Mask = (~(MaskUp | MaskDown)) & 0xfful;
  double ReorderingSign = 1.0;
  if (spinUp > spinDown)
    ReorderingSign = -1.0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState1 = this->StateDescription[i];
      unsigned long TmpState2 = 0x0ul;
      bool Flag = true;
      double Sign = 1.0;
      for (int j = 0; (j <= this->LzMax) && (Flag == true); ++j)
	{
	  if ((TmpState1 & (Mask << (j << 3))) == 0x0ul)
	    {
	      if ((TmpState1 & (MaskUp << (j << 3))) != 0x0ul)
		{
		  if ((TmpState1 & (MaskDown << (j << 3))) != 0x0ul)
		    {
		      TmpState2 |= 0x3ul << (j << 1);
		      Sign *= ReorderingSign;
		    }
		  else
		    {
		      TmpState2 |= 0x1ul << (j << 1);
		    }
		}
	      else
		{
		  if ((TmpState1 & (MaskDown << (j << 3))) != 0x0ul)
		    {
		      TmpState2 |= 0x2ul << (j << 1);
		    }
		}
	    }
	  else
	    {
	      Flag = false;
	    }
	}
      if (Flag == true)
	{
	  int TmpLzMax = (this->LzMax << 1) + 1;
	  while ((TmpState2 >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  int Index = TmpTargetSpace->FindStateIndex(TmpState2, TmpLzMax);
	  if (Index < TmpTargetSpace->HilbertSpaceDimension)
	    {
	      TmpMatrix[i][Index] = Sign;
	    }
	}
    }
  return TmpMatrix;
}

