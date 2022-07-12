////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of fermions on disk with no restriction on the         //
//               number of reachable states or the number of fermions         //
//                                                                            //
//                        last modification : 03/03/2004                      //
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
#include "HilbertSpace/FermionOnDiskUnlimited.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"

#include <math.h>


using std::cout;
using std::endl;
using std::hex;
using std::dec;


// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = momentum total value
// lzMax = maximum angular momentum that a single particle can reach (negative if it has to be deduced from nbrFermions and totalLz)
// memory = amount of memory granted for precalculations

FermionOnDiskUnlimited::FermionOnDiskUnlimited (int nbrFermions, int totalLz, int lzMax, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  if (lzMax < nbrFermions)
    this->LzMax = this->TotalLz - (((this->NbrFermions - 1) * (this->NbrFermions - 2)) / 2);
  else
    this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  this->Flag.Initialize();
  this->StateDescription = new FermionOnSphereLongState [this->HilbertSpaceDimension];
  this->TemporaryStateReducedNbrState = FermionOnSphereLongStateGetReducedNbrState(this->NbrLzValue);
  this->TemporaryState.Resize(this->TemporaryStateReducedNbrState);
  this->ProdATemporaryStateReducedNbrState = this->TemporaryStateReducedNbrState;
  this->ProdATemporaryState.Resize(this->ProdATemporaryStateReducedNbrState);
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->ReducedNbrState = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrFermions, this->LzMax, this->LzMax, this->TotalLz, 0);
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  long UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += ((this->ReducedNbrState[i] + 1) * sizeof(unsigned long));
  UsedMemory += this->HilbertSpaceDimension * (3 * sizeof(int));
  UsedMemory += this->NbrLzValue * sizeof(int);
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      if (this->NbrStateInLookUpTable[i] != 0)
	{
	  for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
	    if (this->NbrStateInLookUpTable[i][j] > 0)      
	      UsedMemory += this->NbrStateInLookUpTable[i][j] * sizeof(int);
	}
    }
  UsedMemory += this->NbrLzValue * sizeof(int*);
  UsedMemory +=  (1 << this->MaximumSignLookUp) * sizeof(double);
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
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

FermionOnDiskUnlimited::FermionOnDiskUnlimited(const FermionOnDiskUnlimited& fermions)
{
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;

  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->ReducedNbrState = fermions.ReducedNbrState;

  this->LookUpTable = fermions.LookUpTable;
  this->NbrStateInLookUpTable = fermions.NbrStateInLookUpTable;
  this->HashKeyMask = fermions.HashKeyMask;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;

  this->TemporaryStateReducedNbrState = FermionOnSphereLongStateGetReducedNbrState(this->NbrLzValue);
  this->TemporaryState.Resize(this->TemporaryStateReducedNbrState);
  this->ProdATemporaryStateReducedNbrState = this->TemporaryStateReducedNbrState;
  this->ProdATemporaryState.Resize(this->ProdATemporaryStateReducedNbrState);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnDiskUnlimited::~FermionOnDiskUnlimited ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnDiskUnlimited& FermionOnDiskUnlimited::operator = (const FermionOnDiskUnlimited& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->ReducedNbrState;
      for (int i = 0; i < this->NbrLzValue; ++i)
	{
	  if (this->NbrStateInLookUpTable[i] != 0)
	    {
	      for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
		if (this->NbrStateInLookUpTable[i][j] > 0)
		  delete[] this->LookUpTable[i][j];
	      delete[] this->LookUpTable[i];
	      delete[] this->NbrStateInLookUpTable[i];
	    }
	}
      delete[] this->LookUpTable;
      delete[] this->NbrStateInLookUpTable;
    }
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
 
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->ReducedNbrState = fermions.ReducedNbrState;

  this->LookUpTable = fermions.LookUpTable;
  this->NbrStateInLookUpTable = fermions.NbrStateInLookUpTable;
  this->HashKeyMask = fermions.HashKeyMask;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;

  this->TemporaryStateReducedNbrState = FermionOnSphereLongStateGetReducedNbrState(this->NbrLzValue);
  this->TemporaryState.Resize(this->TemporaryStateReducedNbrState);
  this->ProdATemporaryStateReducedNbrState = this->TemporaryStateReducedNbrState;
  this->ProdATemporaryState.Resize(this->ProdATemporaryStateReducedNbrState);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}


