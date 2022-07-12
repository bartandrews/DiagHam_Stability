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
//                                                                            //
//                        last modification : 12/08/2007                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
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

FermionOnSphereWithSpinLzSzSymmetry::FermionOnSphereWithSpinLzSzSymmetry ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// lzMax = twice the maximum Lz value reached by a fermion
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// minusLzParity = select the  Lz <-> -Lz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinLzSzSymmetry::FermionOnSphereWithSpinLzSzSymmetry (int nbrFermions, int lzMax, bool minusSzParity, bool minusLzParity, unsigned long memory)
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
#ifdef __64_BITS__
  if ((this->LzMax & 1) == 0)
    {
      this->InvertShift = 32 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 2;
    }
  else
    {
      this->InvertShift = 32 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
#else
  if ((this->LzMax & 1) != 0)
    {
      this->InvertShift = 16 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
  else
    {
      this->InvertShift = 16 - this->LzMax;
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
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						     (this->TotalSpin + this->NbrFermions) >> 1, 0l);
  int TmpHilbertSpaceDimension = 0;
  if (minusSzParity == minusLzParity)
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	if (this->GetCanonicalState(this->StateDescription[i]) != this->StateDescription[i])
	  this->StateDescription[i] = 0x0ul;
	else
	  {
//	    cout << "1 : " << hex << this->StateDescription[i] << dec << endl;
	    unsigned long TmpState = this->StateDescription[i];
	    this->GetStateSymmetry(TmpState);
	    if ((TmpState & FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT) == FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT)
	      ++TmpHilbertSpaceDimension;
	    else
	      {
		unsigned long TmpStateParity = this->StateDescription[i];
		this->GetStateSingletParity(TmpStateParity);
//		cout << "2 : " << hex << TmpStateParity << dec << endl;
		if ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusLzParity == false))
		    || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusLzParity == true)))
		  ++TmpHilbertSpaceDimension;
		else
		  this->StateDescription[i] = 0x0ul;
	      }
	  }
    }
  else
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	if (this->GetCanonicalState(this->StateDescription[i]) != this->StateDescription[i])
	  this->StateDescription[i] = 0x0ul;
	else
	  {
	    unsigned long TmpState = this->StateDescription[i];
	    this->GetStateSymmetry(TmpState);
//	    if ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) > FERMION_SPHERE_SU2_LZSZ_SYMMETRIC_BIT)
	    if (((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) != 0x0ul) && ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) !=  FERMION_SPHERE_SU2_LZSZ_SYMMETRIC_BIT))		      
	      {
		if ((TmpState & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT) == FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT)
		  {
		    ++TmpHilbertSpaceDimension;
		  }
		else
		  {
		    unsigned long TmpStateParity = TmpState;
		    this->GetStateSingletParity(TmpStateParity);
		    if ((((TmpState & FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT) == 0x0ul) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusLzParity == false))
											 || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusLzParity == true))))
			|| (((TmpState & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul) && ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0x0ul) && (minusSzParity == false))
											    || (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0x0ul) && (minusSzParity == true)))))
		      ++TmpHilbertSpaceDimension;
		    else
		      this->StateDescription[i] = 0x0ul;		    
		  }
	      }
	    else
	      this->StateDescription[i] = 0x0ul;
	  }
    }
  cout << "dim = " << TmpHilbertSpaceDimension << endl;
  unsigned long* TmpStateDescription = new unsigned long [TmpHilbertSpaceDimension];
  TmpHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    if (this->StateDescription[i] != 0x0ul)
      {
	TmpStateDescription[TmpHilbertSpaceDimension] = this->StateDescription[i];
	++TmpHilbertSpaceDimension;
      }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  this->HilbertSpaceDimension = TmpHilbertSpaceDimension;

  if (this->HilbertSpaceDimension > 0)
    {
      this->StateHighestBit =  new int [TmpHilbertSpaceDimension];
      this->GenerateLookUpTable(memory);
      delete[] this->StateHighestBit;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	this->GetStateSymmetry(this->StateDescription[i]);
      this->StateHighestBit = 0;
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	this->PrintState(cout, i) << endl;
    }
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
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

FermionOnSphereWithSpinLzSzSymmetry::FermionOnSphereWithSpinLzSzSymmetry(const FermionOnSphereWithSpinLzSzSymmetry& fermions)
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
}

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinLzSzSymmetry::FermionOnSphereWithSpinLzSzSymmetry (char* fileName, unsigned long memory)
{
  this->ReadHilbertSpace(fileName);
  this->IncNbrFermions = this->NbrFermions + 1;
  this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin) / 2;
  this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin) / 2;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
#ifdef __64_BITS__
  if ((this->LzMax & 1) == 0)
    {
      this->InvertShift = 32 - this->LzMax;
      this->InvertUnshift = this->InvertShift - 2;
    }
  else
    {
      this->InvertShift = 32 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
#else
  if ((this->LzMax & 1) != 0)
    {
      this->InvertShift = 16 - (this->LzMax + 1);
      this->InvertUnshift = this->InvertShift;
    }
  else
    {
      this->InvertShift = 16 - this->LzMax;
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
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
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

FermionOnSphereWithSpinLzSzSymmetry::~FermionOnSphereWithSpinLzSzSymmetry ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinLzSzSymmetry& FermionOnSphereWithSpinLzSzSymmetry::operator = (const FermionOnSphereWithSpinLzSzSymmetry& fermions)
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
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinLzSzSymmetry::Clone()
{
  return new FermionOnSphereWithSpinLzSzSymmetry(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphereWithSpinLzSzSymmetry::WriteHilbertSpace (char* fileName)
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

bool FermionOnSphereWithSpinLzSzSymmetry::ReadHilbertSpace (char* fileName)
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
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    ReadLittleEndian(File, this->StateDescription[i]);
  File.close();
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  int NewLzMax = 1 + (2 * this->LzMax);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->StateDescription[i] &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
      unsigned long TmpState = this->StateDescription[i];
      while (((TmpState >> NewLzMax) & 0x1ul) == 0x0ul)
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

RealVector FermionOnSphereWithSpinLzSzSymmetry::ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      Signature = nbodyBasis.StateDescription[i];
      TmpState = this->GetSignedCanonicalState(Signature);
      Signature = TmpState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
      unsigned long TmpState2 = TmpState;
      TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
      NewLzMax = 1 + (this->LzMax << 1);
      while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;	 
      switch (Signature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT)
	{
	case 0x0ul:
	  {
	    if (this->LzSzSameParityFlag == true)
	      {
		Signature = TmpState;
		this->GetStateSingletParity(Signature);
		if (((1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul))) * this->SzParitySign) > 0.0)
		  TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)];
	      }
	  }
	  break;
	case FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT :
	  {
	    if (this->LzSzSameParityFlag == true)
	      {
		Signature = TmpState;
		this->GetStateSingletParity(Signature);
		double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
		TmpVector[i] = ((1.0 + ((double) (((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) 
						   | (TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT)) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
				* state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2);
	      }		
	  }
	  break;
	case FERMION_SPHERE_SU2_SZ_SYMMETRIC_TEST :
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
	    if ((TmpSign * this->LzParitySign) > 0.0)
	      {
		TmpVector[i] = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
				* state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2);
	      }
	  }
	  break;
	case FERMION_SPHERE_SU2_LZ_SYMMETRIC_TEST :
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
	    if ((TmpSign * this->SzParitySign) > 0.0)
	      {
		TmpVector[i] = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->LzParitySign) - 1.0))
				* state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2);
	      }
	  }
	  break;
	default:
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
	    TmpVector[i] = ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->LzParitySign) - 1.0))
			    * (1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
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

RealVector FermionOnSphereWithSpinLzSzSymmetry::ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis)
{
  RealVector TmpVector (this->GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      Signature = nbodyBasis.StateDescription[i];
      TmpState = this->GetSignedCanonicalState(Signature);
      Signature = TmpState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
      unsigned long TmpState2 = TmpState;
      TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
      NewLzMax = 1 + (this->LzMax << 1);
      while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;	 
      switch (Signature & FERMION_SPHERE_SU2_FULLY_SYMMETRIC_BIT)
	{
	case 0x0ul:
	  {
	    if (this->LzSzSameParityFlag == true)
	      {
		Signature = TmpState;
		this->GetStateSingletParity(Signature);
		if (((1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul))) * this->SzParitySign) > 0.0)
		  TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = state[i];
	      }
	  }
	  break;
	case FERMION_SPHERE_SU2_LZ_SZ_SYMMETRIC_BIT :
	  {
	    if (this->LzSzSameParityFlag == true)
	      {
		Signature = TmpState;
		this->GetStateSingletParity(Signature);
		double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
		TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) (((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) 
											   | (TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT)) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
									* state[i] * M_SQRT1_2);
	      }		
	  }
	  break;
	case FERMION_SPHERE_SU2_SZ_SYMMETRIC_TEST :
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
	    if ((TmpSign * this->LzParitySign) > 0.0)
	      {
		TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
									* state[i] * M_SQRT1_2);
	      }
	  }
	  break;
	case FERMION_SPHERE_SU2_LZ_SYMMETRIC_TEST :
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
	    if ((TmpSign * this->SzParitySign) > 0.0)
	      {
		TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->LzParitySign) - 1.0))
									* state[i] * M_SQRT1_2);
	      }
	  }
	  break;
	default:
	  {
	    Signature = TmpState;
	    this->GetStateSingletParity(Signature);
	    double TmpSign = (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)));
	    TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += ((1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_LZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->LzParitySign) - 1.0))
								    * (1.0 + ((double) ((TmpState2 >> FERMION_SPHERE_SU2_SZ_TRANSFORMATION_SHIFT) & 0x1ul)) * ((TmpSign * this->SzParitySign) - 1.0))
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

int FermionOnSphereWithSpinLzSzSymmetry::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  cout << "WARNING: using deprecated method AduAduAuAu" << endl;
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  unsigned long signs = 0x0l;
  //bitset<32> tmpB = State;
  n1 = (n1<<1) + 1;
  n2 = (n2<<1) + 1;  
  //cout << "Examining uuuu: " << tmpB << " LzMax: " << StateHighestBit <<" for n's = (" << n1 << ", "<< n2 << ") coeff: " << coefficient << endl;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & (0x1l << n1)) == 0) 
      || ((State & (0x1l << n2)) == 0) || (n1 == n2) || (m1 == m2)) 
    {
      coefficient = 0.0;
      //cout << "First exit" << endl;
      return this->HilbertSpaceDimension;
    }
  // evaluate bit positions corresponding to (m1,up), (m2,up)
  m1 = (m1<<1) + 1;
  m2 = (m2<<1) + 1;
  int NewLargestBit = StateHighestBit;
  //cout << " m's: (" << m1 << ", " << m2 << ")" << endl;
  signs = State & ((0x1l<<n2) -1); // & mask with all bits set at positions right of n2
  State &= ~(0x1l << n2);
  
  /*tmpB = State;
  cout << "Unset n2:  " << tmpB << endl;
  tmpB = signs;
  cout << "signs:     " << tmpB << endl;*/
  signs ^= State & ((0x1l<<n1) -1);
  State &= ~(0x1l << n1);
  /*tmpB = State;
  cout << "Unset n1:  " << tmpB  << endl;
  tmpB = signs;
  cout << "signs:     " << tmpB << endl;*/
  
  // test if possible to create particles at m1, m2:
  if (((State & (0x1l << m2))!= 0)|| ((State & (0x1l << m1)) != 0))
    {
      //cout << "Third exit" << endl;
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  // recalculate NewLargestBit taking into account the above operations:
  
  if ((NewLargestBit == n2)|| (NewLargestBit == n1))
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  signs ^= State & ((0x1l<<m2)-1);
  State |= (0x1l << m2);
  /*tmpB = State;
  cout << "Set m2:    " << tmpB << endl;
  tmpB = signs;
  cout << "signs:     " << tmpB << endl;*/

  if (m1 > NewLargestBit)
    {
      NewLargestBit = m1;
    }

  // in ParticleOnSphereWithSpin Hamiltonian, always m2>m1! -> we can leave these lines out!
  if (m2 > NewLargestBit) 
    { 
      NewLargestBit = m2; 
    } 
  
  
  signs ^= State & ((0x1l<<m1)-1);
  State |= (0x1l << m1);
  /*tmpB = State;
  cout << "Set m1:    " << tmpB << endl;
  tmpB = signs;
  cout << "signs:     " << tmpB << endl;*/
  
  coefficient = ComputeSign (signs);
  //  cout << "Non-zero result found: "  << coefficient << " New Largest Bit: " << NewLargestBit <<endl;
  return this->FindStateIndex(State, NewLargestBit);
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

int FermionOnSphereWithSpinLzSzSymmetry::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  cout << "WARNING: using deprecated method AddAddAdAd" << endl;
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  unsigned long signs = 0x0l;
  n1 <<= 1;
  n2 <<= 1;  
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & (0x1l << n1)) == 0) 
      || ((State & (0x1l << n2)) == 0) || (n1 == n2) || (m1 == m2)) // the last two are superflous, though adding some security
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  // evaluate bit positions corresponding to (m1,down), (m2,down)
  m1 <<= 1;
  m2 <<= 1;
  int NewLargestBit = StateHighestBit;
  
  signs = State & ((0x1l<<n2) -1); // & mask with all bits set at positions right of n2
  State &= ~(0x1l << n2);
  
  signs ^= State & ((0x1l<<n1) -1);
  State &= ~(0x1l << n1);

  // test if possible to create particles at m1, m2:
  if (((State & (0x1l << m2))!= 0)|| ((State & (0x1l << m1)) != 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  // recalculate NewLargestBit taking into account the above operations:
  
  if ((NewLargestBit == n2)|| (NewLargestBit == n1))
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  signs ^= State & ((0x1l<<m2)-1);
  State |= (0x1l << m2);
  
  if (m1 > NewLargestBit)
    {
      NewLargestBit = m1;
    }

  // if called from ParticleOnSphereWithSpin... Hamiltonian, always m1>m2!
   if (m2 > NewLargestBit)
    {
      NewLargestBit = m2;
    }


  signs ^= State & ((0x1l<<m1)-1);
  State |= (0x1l << m1);

  coefficient = ComputeSign (signs);
  
  return this->FindStateIndex(State, NewLargestBit);

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

int FermionOnSphereWithSpinLzSzSymmetry::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  cout << "WARNING: using deprecated method AddAduAdAu" << endl;
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  unsigned long signs = 0x0l;
  n1 <<= 1;
  n2 = (n2<<1) + 1;  
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & (0x1l << n1)) == 0) 
      || ((State & (0x1l << n2)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  // evaluate bit positions corresponding to (m1,up), (m2,up)
  m1 <<= 1;
  m2 = (m2<<1) + 1;
  int NewLargestBit = StateHighestBit;
  
  signs = State & ((0x1l<<n2) -1); // & mask with all bits set at positions right of n2
  State &= ~(0x1l << n2);
  
  signs ^= State & ((0x1l<<n1) -1);
  State &= ~(0x1l << n1);

  // test if possible to create particles at m1, m2:
  if (((State & (0x1l << m2))!= 0)|| ((State & (0x1l << m1)) != 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  
  // recalculate NewLargestBit taking into account the above operations:
  
  if ((NewLargestBit == n2)|| (NewLargestBit == n1))
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if (m2 > NewLargestBit)
    {
      NewLargestBit = m2;
    }

  signs ^= State & ((0x1l<<m2)-1);
  State |= (0x1l << m2);
  
  if (m1 > NewLargestBit)
    {
      NewLargestBit = m1;
    }

  signs ^= State & ((0x1l<<m1)-1);
  State |= (0x1l << m1);

  coefficient = ComputeSign (signs);
  
  return this->FindStateIndex(State, NewLargestBit);

}


// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereWithSpinLzSzSymmetry::AduAu (int index, int m)
{
  if ((this->StateDescription[index] & (0x2l << (m << 1))) != 0)
    return 1.0;
  else
    return 0.0;
}

// apply a^+_d_m a_d_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_d_m a_d_m

double FermionOnSphereWithSpinLzSzSymmetry::AddAd (int index, int m)
{
  if ((this->StateDescription[index] & (0x1l << (m << 1))) != 0)
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

double FermionOnSphereWithSpinLzSzSymmetry::AuAu (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  ++n1;
  n2 <<= 1;
  ++n2;
  unsigned long TmpMask = (0x1ul << n1) ^ (0x1ul << n2);
  if (((this->ProdATemporaryState & TmpMask) ^ TmpMask) || (n1 == n2))
    return 0.0;
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;

  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  return Coefficient;
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLzSzSymmetry::AdAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  n2 <<= 1;

  unsigned long TmpMask = (0x1ul << n1) ^ (0x1ul << n2);
  if (((this->ProdATemporaryState & TmpMask) ^ TmpMask) || (n1 == n2))
    return 0.0;
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  return Coefficient;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLzSzSymmetry::AuAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 <<= 1;
  ++n1;
  n2 <<= 1;

  unsigned long TmpMask = (0x1ul << n1) ^ (0x1ul << n2);
  if ((this->ProdATemporaryState & TmpMask) ^ TmpMask)
    return 0.0;
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  return Coefficient;
}


// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetry::AduAdu (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  ++m2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
  TmpState |= (0x1ul << m2);
  coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
  TmpState |= (0x1ul << m1);
  
  return this->SymmetrizeAdAdResult(TmpState, coefficient);
}
  
// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetry::AddAdd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  m2 <<= 1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
  TmpState |= (0x1ul << m2);
  coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
  TmpState |= (0x1ul << m1);

  return this->SymmetrizeAdAdResult(TmpState, coefficient);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetry::AduAdd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
  TmpState |= (0x1ul << m2);
  coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
  TmpState |= (0x1ul << m1);

  return this->SymmetrizeAdAdResult(TmpState, coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLzSzSymmetry::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (n[i] << 1) + spinIndices[i];
      if ((this->ProdATemporaryState & (0x1l << Index)) == 0)
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      this->ProdATemporaryState &= ~(0x1l << Index);
    }
  return Coefficient;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each annihilation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereWithSpinLzSzSymmetry::ProdA (int index, int* n, int spinIndices, int nbrIndices)
{
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (n[i] << 1) + ((spinIndices >> i) & 0x1);
      if ((this->ProdATemporaryState & (0x1l << Index)) == 0)
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      this->ProdATemporaryState &= ~(0x1l << Index);
    }
  return Coefficient;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinLzSzSymmetry::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (m[i] << 1) + spinIndices[i];
      if ((TmpState & (0x1l << Index)) != 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      TmpState |= (0x1l << Index);
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

int FermionOnSphereWithSpinLzSzSymmetry::ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = (m[i] << 1) + ((spinIndices >> i) & 0x1);
      if ((TmpState & (0x1l << Index)) != 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      TmpState |= (0x1l << Index);
    }
  return this->SymmetrizeAdAdResult(TmpState, coefficient);
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereWithSpinLzSzSymmetry::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  stateDescription &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
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

Complex FermionOnSphereWithSpinLzSzSymmetry::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
								   int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
#ifdef __USE_LAPACK_HERE__
  ComplexLapackDeterminant SlaterUp(this->NbrFermionsUp);
  ComplexLapackDeterminant SlaterDown(this->NbrFermionsDown);
#else
  ComplexMatrix SlaterUp(this->NbrFermionsUp, this->NbrFermionsUp);
  ComplexMatrix SlaterDown(this->NbrFermionsDown, this->NbrFermionsDown);
#endif
  ComplexMatrix Functions(this->LzMax + 1, this->NbrFermions);
  RealVector TmpCoordinates(2);
  int* IndicesUp = new int [this->NbrFermionsUp];
  int* IndicesDown = new int [this->NbrFermionsDown];
  int PosUp, PosDown;
  int Lz;
  
  // calculate Basis functions for the given set of coordinates:
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
  // calculate prefactor 1/sqrt(N!);
  double Factor = 1.0;
  for (int i = 2; i <= this->NbrFermions; ++i)
    Factor *= (double) i;
  Factor = 1.0 / sqrt(Factor);
  // get down to the business: adding up terms \sum c_\alpha <r|\alpha>
  unsigned long TmpStateDescription;
  int LastComponent = firstComponent + nbrComponent;
  for (int k = firstComponent; k < LastComponent; ++k)
    {
      PosUp = 0;
      PosDown = 0;
      Lz = 0;
      TmpStateDescription = this->StateDescription[k];
      while ((PosUp < this->NbrFermionsUp)||(PosDown < this->NbrFermionsDown))
	{
	  if ((TmpStateDescription & 0x1l) != 0x0l)
	    {
	      IndicesDown[PosDown] = Lz;
	      ++PosDown;
	    }
	  if ((TmpStateDescription & 0x2l) != 0x0l)
	    {
	      IndicesUp[PosUp] = Lz;
	      ++PosUp;
	    }
	  ++Lz;
	  TmpStateDescription >>= 2;
	}
      for (int i = 0; i < this->NbrFermionsUp; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < this->NbrFermionsUp; ++j)
	    {
#ifdef __USE_LAPACK_HERE__
	      SlaterUp.SetMatrixElement(i,j,TmpColum2.Re(IndicesUp[j]), TmpColum2.Im(IndicesUp[j]));
#else
	      SlaterUp[i].Re(j) = TmpColum2.Re(IndicesUp[j]);
	      SlaterUp[i].Im(j) = TmpColum2.Im(IndicesUp[j]);
#endif
	    }
	}
      for (int i = 0; i < this->NbrFermionsDown; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i+this->NbrFermionsUp];	  
	  for (int j = 0; j < this->NbrFermionsDown; ++j)
	    {
#ifdef __USE_LAPACK_HERE__	      
	      SlaterDown.SetMatrixElement(i,j,TmpColum2.Re(IndicesDown[j]), TmpColum2.Im(IndicesDown[j]));
#else
	      SlaterDown[i].Re(j) = TmpColum2.Re(IndicesDown[j]);
	      SlaterDown[i].Im(j) = TmpColum2.Im(IndicesDown[j]);
#endif
	    }
	}
      Complex SlaterDetUp = SlaterUp.Determinant();
      Complex SlaterDetDown = SlaterDown.Determinant();
      Value += SlaterDetUp * SlaterDetDown * (state[k] * Factor) * this->GetStateSign(k, IndicesDown);
    }
  delete[] IndicesUp;
  delete[] IndicesDown;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereWithSpinLzSzSymmetry::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
  
