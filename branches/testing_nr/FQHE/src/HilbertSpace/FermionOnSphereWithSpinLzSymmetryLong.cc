////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of fermions on sphere with spin with             //
//                                Lz<->-Lz symmetry                           //
//    that allow LzMax up to 60 (for systems with 128 bit integer support)    //
//               or 27 (on 32 bit systems without 128 bit integer support)    //
//                                                                            //
//                        last modification : 28/09/2008                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include <math.h>
#include <bitset>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;

#define WANT_LAPACK

#ifdef __LAPACK__
#ifdef WANT_LAPACK
#define  __USE_LAPACK_HERE__
#endif
#endif


// default constructor 
//

FermionOnSphereWithSpinLzSymmetryLong::FermionOnSphereWithSpinLzSymmetryLong ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// lzMax = twice the maximum Lz value reached by a fermion
// totalSz = twce the total spin value
// minusParity = select the  Lz <-> -Lz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinLzSymmetryLong::FermionOnSphereWithSpinLzSymmetryLong (int nbrFermions, int lzMax, int totalSz, bool minusParity, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->TotalSpin = totalSz;
  this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin) / 2;
  this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin) / 2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
#ifdef __64_BITS__
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
  this->LzParitySign = 1.0;
  if (minusParity == true)
    this->LzParitySign = -1.0;
  this->Flag.Initialize();
  this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						     (this->TotalSpin + this->NbrFermions) >> 1, 0l);
  int TmpHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    if (this->GetCanonicalState(this->StateDescription[i]) != this->StateDescription[i])
      this->StateDescription[i] = ((ULONGLONG) 0x0ul);
    else
      {
	this->GetStateSymmetry(this->StateDescription[i]);
	if ((this->StateDescription[i] & FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT_LONG) == ((ULONGLONG) 0x0ul))
	  {
	    ULONGLONG TmpStateParity = this->StateDescription[i];
	    this->GetStateSingletParity(TmpStateParity);
	    if ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) == 0) && (minusParity == false))
		|| (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) != 0) && (minusParity == true)))	
	      ++TmpHilbertSpaceDimension;
	    else
	      this->StateDescription[i] = ((ULONGLONG) 0x0ul);
	  }
	else
	  {
	     ++TmpHilbertSpaceDimension;
	     this->StateDescription[i] &= ~FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT_LONG;
	  }	
      }
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
  if (this->HilbertSpaceDimension > 0)
    {
      this->StateHighestBit =  new int [TmpHilbertSpaceDimension];
      this->GenerateLookUpTable(memory);
      delete[] this->StateHighestBit;
      this->StateHighestBit = 0;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	this->GetStateSymmetry(this->StateDescription[i]);
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
    }
#endif
}

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinLzSymmetryLong::FermionOnSphereWithSpinLzSymmetryLong (char* fileName, unsigned long memory)
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

FermionOnSphereWithSpinLzSymmetryLong::FermionOnSphereWithSpinLzSymmetryLong(const FermionOnSphereWithSpinLzSymmetryLong& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
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
  this->LzParitySign = fermions.LzParitySign;
}

// destructor
//

FermionOnSphereWithSpinLzSymmetryLong::~FermionOnSphereWithSpinLzSymmetryLong ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinLzSymmetryLong& FermionOnSphereWithSpinLzSymmetryLong::operator = (const FermionOnSphereWithSpinLzSymmetryLong& fermions)
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
  this->LzParitySign = fermions.LzParitySign;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinLzSymmetryLong::Clone()
{
  return new FermionOnSphereWithSpinLzSymmetryLong(*this);
}

// convert a given state from symmetric basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector  

RealVector FermionOnSphereWithSpinLzSymmetryLong::ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpinLong& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  ULONGLONG TmpState;
  ULONGLONG Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      Signature = nbodyBasis.StateDescription[i];
      TmpState = this->GetSignedCanonicalState(Signature);
      if ((TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG) == Signature)
	{
	  Signature = TmpState & FERMION_SPHERE_SU2_SYMMETRIC_BIT_LONG;
	  TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
	  NewLzMax = 1 + (this->LzMax << 1);
	  while ((TmpState >> NewLzMax) == ((ULONGLONG) 0x0ul))
	    --NewLzMax;
	  if (Signature != ((ULONGLONG) 0x0ul))	
	    TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2;
	  else
	    {
	      Signature = TmpState;
	      this->GetStateSingletParity(Signature);
	      if ((((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) == ((ULONGLONG) 0x0ul)) && (this->LzParitySign > 0.0))
		  || (((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) != ((ULONGLONG) 0x0ul)) && (this->LzParitySign < 0.0)))
		TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)];
	    }
	}
      else
	{
	  TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
	  NewLzMax = 1 + (this->LzMax << 1);
	  while ((TmpState >> NewLzMax) == ((ULONGLONG) 0x0ul))
	    --NewLzMax;
	  Signature = TmpState;
	  this->GetStateSingletParity(Signature);
	  TmpVector[i] =  (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul)))) * this->LzParitySign * state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2;
	}
    }
  return TmpVector;  
}

// convert a given state from the usual n-body basis to the symmetric basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereWithSpinLzSymmetryLong::ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSpinLong& nbodyBasis)
{
  RealVector TmpVector (this->GetHilbertSpaceDimension(), true);
  ULONGLONG TmpState;
  ULONGLONG Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      Signature = nbodyBasis.StateDescription[i];
      TmpState = this->GetSignedCanonicalState(Signature);
      if ((TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG) == Signature)
	{
	  Signature = TmpState & FERMION_SPHERE_SU2_SYMMETRIC_BIT_LONG;
	  TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
	  NewLzMax = 1 + (this->LzMax << 1);
	  while ((TmpState >> NewLzMax) == ((ULONGLONG) 0x0ul))
	    --NewLzMax;
	  if (Signature != ((ULONGLONG) 0x0ul))	
	    TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += state[i] * M_SQRT1_2;
	  else
	    {
	      Signature = TmpState;
	      this->GetStateSingletParity(Signature);
	      if ((((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) == ((ULONGLONG) 0x0ul)) && (this->LzParitySign > 0.0))
		  || (((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT_LONG) != ((ULONGLONG) 0x0ul)) && (this->LzParitySign < 0.0)))
		TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = state[i];
	    }
	}
      else
	{
	  TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK_LONG;
	  NewLzMax = 1 + (this->LzMax << 1);
	  while ((TmpState >> NewLzMax) == ((ULONGLONG) 0x0ul))
	    --NewLzMax;
	  Signature = TmpState;
	  this->GetStateSingletParity(Signature);
	  TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT_LONG) & ((ULONGLONG) 0x1ul)))) * this->LzParitySign * state[i] * M_SQRT1_2;
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

int FermionOnSphereWithSpinLzSymmetryLong::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

int FermionOnSphereWithSpinLzSymmetryLong::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

int FermionOnSphereWithSpinLzSymmetryLong::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  cout << "WARNING: using deprecated method AddAduAdAu" << endl;
  return this->HilbertSpaceDimension;
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphereWithSpinLzSymmetryLong::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
								   int firstComponent, int nbrComponent)
{
  Complex Value;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereWithSpinLzSymmetryLong::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
  
