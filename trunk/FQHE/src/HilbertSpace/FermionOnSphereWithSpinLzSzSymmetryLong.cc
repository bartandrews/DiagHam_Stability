////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of fermions on sphere with spin with             //
//                           Lz<->-Lz and Sz<->-Sz symmetries                 //
//    that allow LzMax up to 60 (for systems with 128 bit integer support)    //
//               or 27 (on 32 bit systems without 128 bit integer support)    //
//                                                                            //
//                                                                            //
//                        last modification : 26/09/2008                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/Endian.h"

#include <math.h>
#include <bitset>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;
using std::ofstream;
using std::ifstream;
using std::ios;

#define WANT_LAPACK

#ifdef __LAPACK__
#ifdef WANT_LAPACK
#define  __USE_LAPACK_HERE__
#endif
#endif


// default constructor 
//

FermionOnSphereWithSpinLzSzSymmetryLong::FermionOnSphereWithSpinLzSzSymmetryLong ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// lzMax = twice the maximum Lz value reached by a fermion
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinLzSzSymmetryLong::FermionOnSphereWithSpinLzSzSymmetryLong (int nbrFermions, int lzMax, bool minusSzParity, bool minusLzParity, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin) / 2;
  this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin) / 2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
#ifdef __128_BIT_LONGLONG__
  if ((this->LzMax & 1) == 0)
    {
      this->InvertShift = 64 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 2;
    }
  else
    {
      this->InvertShift = 64 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
#else
  if ((this->LzMax & 1) != 0)
    {
      this->InvertShift = 32 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
  else
    {
      this->InvertShift = 32 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 1;
    }
#endif
  this->HilbertSpaceDimension = (int) this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
										 (this->TotalSpin + this->NbrFermions) >> 1);
  this->SzParitySign = 1.0;
  if (minusSzParity == true)
    this->SzParitySign = -1.0;
  this->LzParitySign = 1.0;
  if (minusLzParity == true)
    this->LzParitySign = -1.0;
  if (minusLzParity == minusSzParity)
    this->LzSzSameParityFlag = true;
  else
    this->LzSzSameParityFlag = false;
  this->Flag.Initialize();
  this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						     (this->TotalSpin + this->NbrFermions) >> 1, 0l);
  int TmpHilbertSpaceDimension = 0;
  if (minusSzParity == minusLzParity)
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	if (this->GetCanonicalState(this->StateDescription[i]) != this->StateDescription[i])
	  this->StateDescription[i] = ((ULONGLONG) 0x0ul);
	else
	  {
	    ULONGLONG TmpState = this->StateDescription[i];
	    this->GetStateSymmetry(TmpState);
	    if ((TmpState & FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT_LONG) == FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT_LONG)
	      ++TmpHilbertSpaceDimension;
	    else
	      {
		ULONGLONG TmpStateParity = this->StateDescription[i];
		this->GetStateSingletParity(TmpStateParity);
		if ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) == ((ULONGLONG) 0x0ul)) && (minusLzParity == false))
		    || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) != ((ULONGLONG) 0x0ul)) && (minusLzParity == true)))
		  ++TmpHilbertSpaceDimension;
		else
		  this->StateDescription[i] = ((ULONGLONG) 0x0ul);
	      }
	  }
    }
  else
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	if (this->GetCanonicalState(this->StateDescription[i]) != this->StateDescription[i])
	  this->StateDescription[i] = ((ULONGLONG) 0x0ul);
	else
	  {
	    ULONGLONG TmpState = this->StateDescription[i];
	    this->GetStateSymmetry(TmpState);
	    if (((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) != ((ULONGLONG) 0x0ul)) && ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) !=  FERMION_SPHERE_SU2_LZSZ_SYMMETRIC_BIT_LONG))		      
	      {
		if ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG) == FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG)
		  {
		    ++TmpHilbertSpaceDimension;
		  }
		else
		  {
		    ULONGLONG TmpStateParity = TmpState;
		    this->GetStateSingletParity(TmpStateParity);
		    if ((((TmpState & FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT_LONG) == ((ULONGLONG) 0x0ul)) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) == ((ULONGLONG) 0x0ul)) && (minusLzParity == false))
											 || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) != ((ULONGLONG) 0x0ul)) && (minusLzParity == true))))
			|| (((TmpState & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT_LONG) == ((ULONGLONG) 0x0ul)) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) == ((ULONGLONG) 0x0ul)) && (minusSzParity == false))
											    || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) != ((ULONGLONG) 0x0ul)) && (minusSzParity == true)))))
		      ++TmpHilbertSpaceDimension;
		    else
		      this->StateDescription[i] = ((ULONGLONG) 0x0ul);		    
		  }
	      }
	    else
	      this->StateDescription[i] = ((ULONGLONG) 0x0ul);
	  }
    }
  cout << "dim = " << TmpHilbertSpaceDimension << endl;
  ULONGLONG* TmpStateDescription = new ULONGLONG [TmpHilbertSpaceDimension];
  TmpHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    if (this->StateDescription[i] != ((ULONGLONG) 0x0ul))
      {
	TmpStateDescription[TmpHilbertSpaceDimension] = this->StateDescription[i];
	++TmpHilbertSpaceDimension;
      }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  this->HilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

  if (this->HilbertSpaceDimension > 0)
    {
      this->StateHighestBit =  new int [TmpHilbertSpaceDimension];
      this->GenerateLookUpTable(memory);
      delete[] this->StateHighestBit;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	this->GetStateSymmetry(this->StateDescription[i]);
      this->StateHighestBit = 0;
    }
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += this->HilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereWithSpinLzSzSymmetryLong::FermionOnSphereWithSpinLzSzSymmetryLong(const FermionOnSphereWithSpinLzSzSymmetryLong& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
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
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->LzParitySign = fermions.LzParitySign;
  this->SzParitySign = fermions.SzParitySign;
  this->LzSzSameParityFlag = fermions.LzSzSameParityFlag;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinLzSzSymmetryLong::FermionOnSphereWithSpinLzSzSymmetryLong (char* fileName, unsigned long memory)
{
  this->ReadHilbertSpace(fileName);
  this->IncNbrFermions = this->NbrFermions + 1;
  this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin) / 2;
  this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin) / 2;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
#ifdef __128_BIT_LONGLONG__
  if ((this->LzMax & 1) == 0)
    {
      this->InvertShift = 64 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 2;
    }
  else
    {
      this->InvertShift = 64 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
#else
  if ((this->LzMax & 1) != 0)
    {
      this->InvertShift = 32 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
  else
    {
      this->InvertShift = 32 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 1;
    }
#endif
  if (this->HilbertSpaceDimension > 0)
    {
      this->GenerateLookUpTable(memory);
      delete[] this->StateHighestBit;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	this->GetStateSymmetry(this->StateDescription[i]);
      this->StateHighestBit = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += this->HilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
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

// destructor
//

FermionOnSphereWithSpinLzSzSymmetryLong::~FermionOnSphereWithSpinLzSzSymmetryLong ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinLzSzSymmetryLong& FermionOnSphereWithSpinLzSzSymmetryLong::operator = (const FermionOnSphereWithSpinLzSzSymmetryLong& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->LzParitySign = fermions.LzParitySign;
  this->SzParitySign = fermions.SzParitySign;
  this->LzSzSameParityFlag = fermions.LzSzSameParityFlag;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinLzSzSymmetryLong::Clone()
{
  return new FermionOnSphereWithSpinLzSzSymmetryLong(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphereWithSpinLzSzSymmetryLong::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->NbrFermions);
  WriteLittleEndian(File, this->LzMax);
  WriteLittleEndian(File, this->TotalLz);
  WriteLittleEndian(File, this->TotalSpin);
  WriteLittleEndian(File, this->SzParitySign);
  WriteLittleEndian(File, this->LzParitySign);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    WriteLittleEndian(File, this->StateDescription[i]);
  File.close();
  return true;
}

// read Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description is stored
// return value = true if no error occured

bool FermionOnSphereWithSpinLzSzSymmetryLong::ReadHilbertSpace (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      this->HilbertSpaceDimension = 0;
      return false;
    }
  ReadLittleEndian(File, this->HilbertSpaceDimension);
  ReadLittleEndian(File, this->NbrFermions);
  ReadLittleEndian(File, this->LzMax);
  ReadLittleEndian(File, this->TotalLz);
  ReadLittleEndian(File, this->TotalSpin);
  ReadLittleEndian(File, this->SzParitySign);
  ReadLittleEndian(File, this->LzParitySign);
  if (this->SzParitySign == this->LzParitySign)
    this->LzSzSameParityFlag = true;
  else
    this->LzSzSameParityFlag = false;
  this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    ReadLittleEndian(File, this->StateDescription[i]);
  File.close();
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  int NewLzMax = 1 + (2 * this->LzMax);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->StateDescription[i] &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
      ULONGLONG TmpState = this->StateDescription[i];
      while (((TmpState >> NewLzMax) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
	--NewLzMax;
      this->StateHighestBit[i] = NewLzMax;
    }
  return true;
}
  
// convert a given state from symmetric basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector  

RealVector FermionOnSphereWithSpinLzSzSymmetryLong::ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpinLong& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  ULONGLONG TmpState;
  ULONGLONG Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      Signature = nbodyBasis.StateDescription[i];
      TmpState = this->GetSignedCanonicalState(Signature);
      Signature = TmpState & FERMION_SPHERE_SU2_SYMMETRIC_BIT_LONG;
      ULONGLONG TmpState2 = TmpState;
      TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
      NewLzMax = 1 + (this->LzMax << 1);
      while ((TmpState >> NewLzMax) == ((ULONGLONG) 0x0ul))
	--NewLzMax;	 
      switch (Signature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG)
	{
	case ((ULONGLONG) 0x0ul):
	  {
	    if (this->LzSzSameParityFlag == true)
	      {
		Signature = TmpState;
		this->GetStateSingletParity(Signature);
		if (((1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul)))) * this->SzParitySign) > 0.0)
		  TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)];
	      }
	  }
	  break;
	case FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT_LONG :
	  {
	    if (this->LzSzSameParityFlag == true)
	      {
		Signature = TmpState;
		this->GetStateSingletParity(Signature);
		double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul))));
		TmpVector[i] = ((1.0 + ((double) (((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT_LONG) 
						   | (TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT_LONG)) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->SzParitySign) - 1.0))
				* state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2);
	      }		
	  }
	  break;
	case FERMION_SPHERE_SU2_SZ_SYMMETRIC_TEST_LONG :
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul))));
	    if ((TmpSign * this->LzParitySign) > 0.0)
	      {
		TmpVector[i] = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->SzParitySign) - 1.0))
				* state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2);
	      }
	  }
	  break;
	case FERMION_SPHERE_SU2_LZ_SYMMETRIC_TEST_LONG :
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul))));
	    if ((TmpSign * this->SzParitySign) > 0.0)
	      {
		TmpVector[i] = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->LzParitySign) - 1.0))
				* state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2);
	      }
	  }
	  break;
	default:
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul))));
	    TmpVector[i] = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->LzParitySign) - 1.0))
			    * (1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->SzParitySign) - 1.0))
			    * state[this->FindStateIndex(TmpState, NewLzMax)] * 0.5);
	  }
	  break;
	}
    }
  return TmpVector;  
}

// convert a given state from the usual n-body basis to the symmetric basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereWithSpinLzSzSymmetryLong::ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSpinLong& nbodyBasis)
{
  RealVector TmpVector (this->GetHilbertSpaceDimension(), true);
  ULONGLONG TmpState;
  ULONGLONG Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      Signature = nbodyBasis.StateDescription[i];
      TmpState = this->GetSignedCanonicalState(Signature);
      Signature = TmpState & FERMION_SPHERE_SU2_SYMMETRIC_BIT_LONG;
      ULONGLONG TmpState2 = TmpState;
      TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
      NewLzMax = 1 + (this->LzMax << 1);
      while ((TmpState >> NewLzMax) == ((ULONGLONG) 0x0ul))
	--NewLzMax;	 
      switch (Signature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT_LONG)
	{
	case ((ULONGLONG) 0x0ul):
	  {
	    if (this->LzSzSameParityFlag == true)
	      {
		Signature = TmpState;
		this->GetStateSingletParity(Signature);
		if (((1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul)))) * this->SzParitySign) > 0.0)
		  TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = state[i];
	      }
	  }
	  break;
	case FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT_LONG :
	  {
	    if (this->LzSzSameParityFlag == true)
	      {
		Signature = TmpState;
		this->GetStateSingletParity(Signature);
		double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul))));
		TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) (((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT_LONG) 
											   | (TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT_LONG)) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->SzParitySign) - 1.0))
									* state[i] * M_SQRT1_2);
	      }		
	  }
	  break;
	case FERMION_SPHERE_SU2_SZ_SYMMETRIC_TEST_LONG :
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul))));
	    if ((TmpSign * this->LzParitySign) > 0.0)
	      {
		TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->SzParitySign) - 1.0))
									* state[i] * M_SQRT1_2);
	      }
	  }
	  break;
	case FERMION_SPHERE_SU2_LZ_SYMMETRIC_TEST_LONG :
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul))));
	    if ((TmpSign * this->SzParitySign) > 0.0)
	      {
		TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->LzParitySign) - 1.0))
									* state[i] * M_SQRT1_2);
	      }
	  }
	  break;
	default:
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul))));
	    TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->LzParitySign) - 1.0))
								    * (1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT_LONG) & ((ULONGLONG) 0x1ul))) * ((TmpSign * this->SzParitySign) - 1.0))
								    * state[i] * 0.5);
	  }
	  break;
	}
    }
  return TmpVector;  
}

// apply a^+_u_m1 a^+_u_m2 a_u_n1 a_u_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetryLong::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  cout << "WARNING: using deprecated method AduAduAuAu" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_d_m1 a^+_d_m2 a_d_n1 a_d_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetryLong::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  cout << "WARNING: using deprecated method AddAddAdAd" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_d_m1 a^+_u_m2 a_d_n1 a_u_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetryLong::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  cout << "WARNING: using deprecated method AddAduAdAu" << endl;
  return this->HilbertSpaceDimension;
}


// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereWithSpinLzSzSymmetryLong::AduAu (int index, int m)
{
  if ((this->StateDescription[index] & (((ULONGLONG) 0x2ul) << (m << 1))) != ((ULONGLONG) 0x0ul))
    return 1.0;
  else
    return 0.0;
}

// apply a^+_d_m a_d_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_d_m a_d_m

double FermionOnSphereWithSpinLzSzSymmetryLong::AddAd (int index, int m)
{
  if ((this->StateDescription[index] & (((ULONGLONG) 0x1ul) << (m << 1))) != ((ULONGLONG) 0x0ul))
    return 1.0;
  else
    return 0.0;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLzSzSymmetryLong::AuAu (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  ++n1;
  n2 <<= 1;
  ++n2;
  ULONGLONG TmpMask = (((ULONGLONG) 0x1ul) << n1) ^ (((ULONGLONG) 0x1ul) << n2);
  if (((this->ProdATemporaryState & TmpMask) ^ TmpMask) || (n1 == n2))
    return 0.0;
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT_LONG;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;

  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  return Coefficient;
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLzSzSymmetryLong::AdAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  n2 <<= 1;

  ULONGLONG TmpMask = (((ULONGLONG) 0x1ul) << n1) ^ (((ULONGLONG) 0x1ul) << n2);
  if (((this->ProdATemporaryState & TmpMask) ^ TmpMask) || (n1 == n2))
    return 0.0;
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT_LONG;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  return Coefficient;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLzSzSymmetryLong::AuAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  ++n1;
  n2 <<= 1;

  ULONGLONG TmpMask = (((ULONGLONG) 0x1ul) << n1) ^ (((ULONGLONG) 0x1ul) << n2);
  if ((this->ProdATemporaryState & TmpMask) ^ TmpMask)
    return 0.0;
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT_LONG;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  return Coefficient;
}


// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetryLong::AduAdu (int m1, int m2, double& coefficient)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  ++m2;
  if (((TmpState & (((ULONGLONG) 0x1ul) << m1)) != ((ULONGLONG) 0x0ul)) || ((TmpState & (((ULONGLONG) 0x1ul) << m2)) != ((ULONGLONG) 0x0ul)) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80))  & this->SignLookUpTableMask[m2 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80))  & this->SignLookUpTableMask[m1 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  
  return this->SymmetrizeAdAdResult(TmpState, coefficient);
}
  
// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetryLong::AddAdd (int m1, int m2, double& coefficient)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  m2 <<= 1;
  if (((TmpState & (((ULONGLONG) 0x1ul) << m1)) != ((ULONGLONG) 0x0ul)) || ((TmpState & (((ULONGLONG) 0x1ul) << m2)) != ((ULONGLONG) 0x0ul)) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80))  & this->SignLookUpTableMask[m2 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80))  & this->SignLookUpTableMask[m1 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
  TmpState |= (((ULONGLONG) 0x1ul) << m1);

  return this->SymmetrizeAdAdResult(TmpState, coefficient);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetryLong::AduAdd (int m1, int m2, double& coefficient)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  if (((TmpState & (((ULONGLONG) 0x1ul) << m1)) != ((ULONGLONG) 0x0ul)) || ((TmpState & (((ULONGLONG) 0x1ul) << m2)) != ((ULONGLONG) 0x0ul)) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80))  & this->SignLookUpTableMask[m2 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80))  & this->SignLookUpTableMask[m1 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
  TmpState |= (((ULONGLONG) 0x1ul) << m1);

  return this->SymmetrizeAdAdResult(TmpState, coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLzSzSymmetryLong::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  this->ProdALzMax = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT_LONG;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (n[i] << 1) + spinIndices[i];
      if ((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << Index)) == ((ULONGLONG) 0x0ul))
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#ifdef __128_BIT_LONGLONG__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 64)) & this->SignLookUpTableMask[Index + 64]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 80))  & this->SignLookUpTableMask[Index + 80]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 96)) & this->SignLookUpTableMask[Index + 96]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 112)) & this->SignLookUpTableMask[Index + 112]];
#endif
      this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << Index);
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == ((ULONGLONG) 0x0ul))
    --this->ProdALzMax;

  return Coefficient;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each annihilation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLzSzSymmetryLong::ProdA (int index, int* n, int spinIndices, int nbrIndices)
{
  this->ProdALzMax = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT_LONG;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (n[i] << 1) + ((spinIndices >> i) & 0x1);
      if ((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << Index)) == ((ULONGLONG) 0x0ul))
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#ifdef __128_BIT_LONGLONG__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 64)) & this->SignLookUpTableMask[Index + 64]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 80))  & this->SignLookUpTableMask[Index + 80]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 96)) & this->SignLookUpTableMask[Index + 96]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 112)) & this->SignLookUpTableMask[Index + 112]];
#endif
      this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << Index);
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == ((ULONGLONG) 0x0ul))
    --this->ProdALzMax;

  return Coefficient;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetryLong::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  ULONGLONG TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (m[i] << 1) + spinIndices[i];
      if ((TmpState & (((ULONGLONG) 0x1ul) << Index)) != ((ULONGLONG) 0x0ul))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#ifdef __128_BIT_LONGLONG__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 64)) & this->SignLookUpTableMask[Index + 64]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 80))  & this->SignLookUpTableMask[Index + 80]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 96)) & this->SignLookUpTableMask[Index + 96]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 112)) & this->SignLookUpTableMask[Index + 112]];
#endif
	}
      TmpState |= (((ULONGLONG) 0x1ul) << Index);
    }
  return this->SymmetrizeAdAdResult(TmpState, coefficient);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each creation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetryLong::ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  ULONGLONG TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (m[i] << 1) + ((spinIndices >> i) & 0x1);
      if ((TmpState & (((ULONGLONG) 0x1ul) << Index)) != ((ULONGLONG) 0x0ul))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#ifdef __128_BIT_LONGLONG__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 64)) & this->SignLookUpTableMask[Index + 64]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 80))  & this->SignLookUpTableMask[Index + 80]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 96)) & this->SignLookUpTableMask[Index + 96]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 112)) & this->SignLookUpTableMask[Index + 112]];
#endif
	}
      TmpState |= (((ULONGLONG) 0x1ul) << Index);
    }
  return this->SymmetrizeAdAdResult(TmpState, coefficient);
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereWithSpinLzSzSymmetryLong::FindStateIndex(ULONGLONG stateDescription, int lzmax)
{
  stateDescription &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  ULONGLONG CurrentState = (this->StateDescription[PosMid] & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG);
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
      CurrentState = (this->StateDescription[PosMid] & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG);
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    return PosMin;
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphereWithSpinLzSzSymmetryLong::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
								   int firstComponent, int nbrComponent)
{
  Complex Value;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereWithSpinLzSzSymmetryLong::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
  
