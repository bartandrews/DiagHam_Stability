////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//              class of fermions on sphere with SU(3) spin including         //
//                      the Z3 and Tz discrete symmetries                     //
//                                                                            //
//                        last modification : 20/02/2008                      //
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
#include "HilbertSpace/FermionOnSphereWithSU3SpinTzZ3Symmetry.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/Endian.h"

#include <math.h>
#include <fstream>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

FermionOnSphereWithSU3SpinTzZ3Symmetry::FermionOnSphereWithSU3SpinTzZ3Symmetry()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// minusTzParity = select the  Tz <-> -Tz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnSphereWithSU3SpinTzZ3Symmetry::FermionOnSphereWithSU3SpinTzZ3Symmetry (int nbrFermions, int totalLz, int lzMax, bool minusTzParity, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalTz = 0;
  this->TotalY = 0;
  this->TzParitySign = 1.0;
  this->LzParitySign = 1.0;
  this->YParitySign = 1.0;
  if (minusTzParity == true)
    this->TzParitySign = -1.0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  int N1 = (2 * this->NbrFermions) + this->TotalY + (3 * this->TotalTz);
  int N2 = (2 * this->NbrFermions) + this->TotalY - (3 * this->TotalTz);
  int N3 = this->NbrFermions - this->TotalY;
  if ((N1 < 0) || (N2 < 0) || (N3 < 0) || ((N1 % 6) != 0) || ((N2 % 6) != 0) || ((N3 % 3) != 0))
    this->HilbertSpaceDimension = 0;
  else
    {
      N1 /= 6;
      N2 /= 6;
      N3 /= 3;
      this->HilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, N1, N2, N3);
    }
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						       N1, N2, N3, 0l);
  if (TmpHilbertSpaceDimension != this->HilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->HilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in FermionOnSphereWithSU3SpinTzZ3Symmetry!" << endl;
      exit(1);
    }

  TmpHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (this->GetCanonicalState(this->StateDescription[i]) != this->StateDescription[i])
	this->StateDescription[i] = 0x0ul;
      else
	{
	  if ((this->GetStateSymmetry(this->StateDescription[i]) & FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT) != 0x0ul)
	    {
	      if ((this->GetStateDoubletTripletParity(this->StateDescription[i]) * this->TzParitySign) > 0.0)
		++TmpHilbertSpaceDimension;
	      else
		this->StateDescription[i] = 0x0ul;
	    }
	  else
	    ++TmpHilbertSpaceDimension;
	}
    }
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  if (this->HilbertSpaceDimension > 0)
    {
      this->StateHighestBit =  new int [TmpHilbertSpaceDimension];
      this->GenerateLookUpTable(memory);
      delete[] this->StateHighestBit;
      this->StateHighestBit = 0;
      cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
      //      for (int i = 0; i < this->HilbertSpaceDimension; ++i)	
      //	{	  
      //	  this->PrintState(cout, i) << hex << this->StateDescription[i] << dec << " sign = " << this->GetStateSymmetry(this->StateDescription[i]) << endl;
      //	}
// " r = " << this->GetStateRotationSign(this->StateDescription[i], 0x0ul) 
// 				   << " l = " << this->GetStateRotationSign(this->StateDescription[i], FERMION_SPHERE_SU3_Z3LEFTROTATION_BIT) << endl;
#ifdef __DEBUG__
      long UsedMemory = 0;
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
    }
#endif
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereWithSU3SpinTzZ3Symmetry::FermionOnSphereWithSU3SpinTzZ3Symmetry(const FermionOnSphereWithSU3SpinTzZ3Symmetry& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalTz = fermions.TotalTz;
  this->TotalY = fermions.TotalY;
  this->LzParitySign = fermions.LzParitySign;
  this->TzParitySign = fermions.TzParitySign;
  this->YParitySign = fermions.YParitySign;
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

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

FermionOnSphereWithSU3SpinTzZ3Symmetry::FermionOnSphereWithSU3SpinTzZ3Symmetry (char* fileName, unsigned long memory)
{
  this->ReadHilbertSpace(fileName);
  this->IncNbrFermions = this->NbrFermions + 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  if (this->HilbertSpaceDimension > 0)
    {
      this->GenerateLookUpTable(memory);
      delete[] this->StateHighestBit;
      this->StateHighestBit = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  long UsedMemory = 0;
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

FermionOnSphereWithSU3SpinTzZ3Symmetry::~FermionOnSphereWithSU3SpinTzZ3Symmetry ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSU3SpinTzZ3Symmetry& FermionOnSphereWithSU3SpinTzZ3Symmetry::operator = (const FermionOnSphereWithSU3SpinTzZ3Symmetry& fermions)
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
  this->TotalTz = fermions.TotalTz;
  this->TotalY = fermions.TotalY;
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

AbstractHilbertSpace* FermionOnSphereWithSU3SpinTzZ3Symmetry::Clone()
{
  return new FermionOnSphereWithSU3SpinTzZ3Symmetry(*this);
}

// convert a given state from symmetric basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector  

RealVector FermionOnSphereWithSU3SpinTzZ3Symmetry::ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSU3Spin& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      TmpState = nbodyBasis.StateDescription[i];
      unsigned long TmpState2 = TmpState;
      Signature = this->GetSignedCanonicalState(TmpState);
      NewLzMax = 2 + (this->LzMax * 3);
      while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;	 
      switch (Signature & FERMION_SPHERE_SU3_TZZ3_SYMMETRIC_BIT)
	{
	case  FERMION_SPHERE_SU3_TZZ3_SYMMETRIC_BIT :
	  {
	    if ((this->GetStateDoubletTripletParity(TmpState2) * this->TzParitySign) > 0.0)
	      {
		TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)];
	      }
	  }
	  break;
	case FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT :
	  {
	    if ((this->GetStateDoubletTripletParity(TmpState2) * this->TzParitySign) > 0.0)
	      {
		if (TmpState2 != TmpState)
		  TmpVector[i] = this->GetStateRotationSign(TmpState2, Signature) * state[this->FindStateIndex(TmpState, NewLzMax)] * MSQRT1_3;
		else
		  TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)] * MSQRT1_3;
	      }
	  }
	  break;
	default:
	  {
	    double TmpCoefficient = 1.0;
	    if ((Signature & FERMION_SPHERE_SU3_TZ_FLIP_BIT) != 0x0ul)
	      {
		TmpCoefficient *= this->GetStateSingletTzParity(TmpState2);
		TmpCoefficient *= this->TzParitySign;
	      }
	    if ((Signature & FERMION_SPHERE_SU3_Z3ROTATION_BIT) != 0x0ul)
	      {
		TmpCoefficient *= this->GetStateRotationSign(TmpState2, Signature);
	      }
	    TmpVector[i] = TmpCoefficient * state[this->FindStateIndex(TmpState, NewLzMax)] * MSQRT1_6;
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

RealVector FermionOnSphereWithSU3SpinTzZ3Symmetry::ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSU3Spin& nbodyBasis)
{
  RealVector TmpVector (this->GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      TmpState = nbodyBasis.StateDescription[i];
      unsigned long TmpState2 = TmpState;
      Signature = this->GetSignedCanonicalState(TmpState);
      NewLzMax = 2 + (this->LzMax * 3);
      while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;	 
      switch (Signature & FERMION_SPHERE_SU3_TZZ3_SYMMETRIC_BIT)
	{
	case  FERMION_SPHERE_SU3_TZZ3_SYMMETRIC_BIT :
	  {
	    if ((this->GetStateDoubletTripletParity(TmpState2) * this->TzParitySign) > 0.0)
	      {
		TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = state[i];
	      }
	  }
	  break;
	case FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT :
	  {
	    if ((this->GetStateDoubletTripletParity(TmpState2) * this->TzParitySign) > 0.0)
	      {
		if (TmpState2 != TmpState)
		  TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = this->GetStateRotationSign(TmpState2, Signature) * state[i] * MSQRT1_3;
		else
		  TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = state[i] * MSQRT1_3;
	      }
	  }
	  break;
	default:
	  {
	    double TmpCoefficient = 1.0;
	    if ((Signature & FERMION_SPHERE_SU3_TZ_FLIP_BIT) != 0x0ul)
	      {
		TmpCoefficient *= this->GetStateSingletTzParity(TmpState2);
		TmpCoefficient *= this->TzParitySign;
	      }
	    if ((Signature & FERMION_SPHERE_SU3_Z3ROTATION_BIT) != 0x0ul)
	      {
		TmpCoefficient *= this->GetStateRotationSign(TmpState2, Signature);
	      }
	    TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = TmpCoefficient * state[i] * MSQRT1_6;
	  }
	  break;
	}
    }
  return TmpVector;  
}
