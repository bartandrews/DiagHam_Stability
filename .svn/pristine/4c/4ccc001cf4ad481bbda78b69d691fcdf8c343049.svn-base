////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of fermions on sphere with spin with             //
//                                Sz<->-Sz symmetry                           //
//                                                                            //
//                        last modification : 13/08/2007                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
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

FermionOnSphereWithSpinSzSymmetry::FermionOnSphereWithSpinSzSymmetry ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// minusParity = select the  Sz <-> -Sz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinSzSymmetry::FermionOnSphereWithSpinSzSymmetry (int nbrFermions, int totalLz, int lzMax, bool minusParity, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin) / 2;
  this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin) / 2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->HilbertSpaceDimension = (int) this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
										 (this->TotalSpin + this->NbrFermions) >> 1);
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
						     (this->TotalSpin + this->NbrFermions) >> 1, 0l);
  this->SzParitySign = 1.0;
  if (minusParity == true)
    this->SzParitySign = -1.0;
  int TmpHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    if (this->GetCanonicalState(this->StateDescription[i]) != this->StateDescription[i])
      this->StateDescription[i] = 0x0ul;
    else
      {
	this->GetStateSymmetry(this->StateDescription[i]);
	if ((this->StateDescription[i] & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul)
	  {
	    unsigned long TmpStateParity = this->StateDescription[i];
	    this->GetStateSingletParity(TmpStateParity);
	    if ((((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0) && (minusParity == false))
		|| (((TmpStateParity & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0) && (minusParity == true)))
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  if (this->HilbertSpaceDimension>0)
    {
      this->StateHighestBit =  new int [TmpHilbertSpaceDimension];
      this->GenerateLookUpTable(memory);
      delete[] this->StateHighestBit;
      this->StateHighestBit = 0;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	this->GetStateSymmetry(this->StateDescription[i]);
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
  else cout << "Hilbert space empty" << endl;
}

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinSzSymmetry::FermionOnSphereWithSpinSzSymmetry (char* fileName, unsigned long memory)
{
  this->ReadHilbertSpace(fileName);
  this->IncNbrFermions = this->NbrFermions + 1;
  this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin) / 2;
  this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin) / 2;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereWithSpinSzSymmetry::FermionOnSphereWithSpinSzSymmetry(const FermionOnSphereWithSpinSzSymmetry& fermions)
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
  this->SzParitySign = fermions.SzParitySign;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnSphereWithSpinSzSymmetry::~FermionOnSphereWithSpinSzSymmetry ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinSzSymmetry& FermionOnSphereWithSpinSzSymmetry::operator = (const FermionOnSphereWithSpinSzSymmetry& fermions)
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
  this->SzParitySign = fermions.SzParitySign;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinSzSymmetry::Clone()
{
  return new FermionOnSphereWithSpinSzSymmetry(*this);
}

// convert a given state from symmetric basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector  

RealVector FermionOnSphereWithSpinSzSymmetry::ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      Signature = nbodyBasis.StateDescription[i];
      TmpState = this->GetSignedCanonicalState(Signature);
      if ((TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK) == Signature)
        {
          Signature = TmpState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
          TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
	  NewLzMax = 1 + (this->LzMax << 1);
          while ((TmpState >> NewLzMax) == 0x0ul)
            --NewLzMax;
          if (Signature != 0x0ul)	
            TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2;
          else
	    {
	      Signature = TmpState;
	      this->GetStateSingletParity(Signature);
	      if ((((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0) && (this->SzParitySign > 0.0))
		  || (((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0) && (this->SzParitySign < 0.0)))
		TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)];
	    }
        }
      else
        {
          TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
	  NewLzMax = 1 + (this->LzMax << 1);
          while ((TmpState >> NewLzMax) == 0x0ul)
            --NewLzMax;
	  Signature = TmpState;
	  this->GetStateSingletParity(Signature);
          TmpVector[i] =  (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul))) * this->SzParitySign * state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2;
        }
    }
  return TmpVector;  
}

// convert a given state from the usual n-body basis to the symmetric basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereWithSpinSzSymmetry::ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSpin& nbodyBasis)
{
  RealVector TmpVector (this->GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      Signature = nbodyBasis.StateDescription[i];
      TmpState = this->GetSignedCanonicalState(Signature);
      if ((TmpState & FERMION_SPHERE_SU2_SYMMETRIC_MASK) == Signature)
	{	  
	  Signature = TmpState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
	  TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
	  NewLzMax = 1 + (this->LzMax << 1);
	  while ((TmpState >> NewLzMax) == 0x0ul)
	    --NewLzMax;
	  if ((Signature & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) == 0x0ul)	
	    {
	      Signature = TmpState;
	      this->GetStateSingletParity(Signature);
	      if ((((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0) && (this->SzParitySign > 0.0))
		  || (((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0) && (this->SzParitySign < 0.0)))
		TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = state[i];
	    }
	  else
	    TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += state[i] * M_SQRT1_2;
	}
      else
	{
	  TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
	  NewLzMax = 1 + (this->LzMax << 1);
	  while ((TmpState >> NewLzMax) == 0x0ul)
	    --NewLzMax;
	  Signature = TmpState;
	  this->GetStateSingletParity(Signature);
	  TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul)))* this->SzParitySign * state[i] * M_SQRT1_2;
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

int FermionOnSphereWithSpinSzSymmetry::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

int FermionOnSphereWithSpinSzSymmetry::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

int FermionOnSphereWithSpinSzSymmetry::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphereWithSpinSzSymmetry::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
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

void FermionOnSphereWithSpinSzSymmetry::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
  
// apply a_n_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator (spin up)
// return value =  multiplicative factor 

double FermionOnSphereWithSpinSzSymmetry::Au (int index, int n)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n <<= 1;
  ++n;
  
  if ((this->ProdATemporaryState & (0x1ul << n)) == 0)
    return 0.0;
  
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  
//   this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n) & this->SignLookUpTableMask[n]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n);
//   if (this->ProdATemporaryState != 0x0ul)
//     {
//       while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
// 	--this->ProdALzMax;
//     }
//   else
//     this->ProdALzMax = 0;
  return Coefficient;
}


// apply a_n_d  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereWithSpinSzSymmetry::Ad (int index, int n)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n <<= 1;
  
  if ((this->ProdATemporaryState & (0x1ul << n)) == 0)
    return 0.0;
  
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  
//   this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n) & this->SignLookUpTableMask[n]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n);
//   if (this->ProdATemporaryState != 0x0ul)
//     {
//       while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
// 	--this->ProdALzMax;
//     }
//   else
//     this->ProdALzMax = 0;
  return Coefficient;
}

// apply a^+_m_u operator to the state produced using Au method (without destroying it)
//
// m = first index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinSzSymmetry::Adu (int m, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m <<= 1;
  ++m;
  
  if ((TmpState & (0x1ul << m)) != 0)
    return this->HilbertSpaceDimension;
//   int NewLzMax = this->ProdALzMax;
  
//   if (m > NewLzMax)
//     NewLzMax = m;
//   else
//     {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
//     }
  TmpState |= (0x1ul << m);
  
  
//   cout <<  "ad" << m << " -> " ;
//   for (int i = 0; i < 2*this->NbrLzValue; ++i)
//     cout << ((TmpState >> i) & 0x1ul) << " " ;
//   cout << (TmpState & FERMION_SPHERE_SU2_SZ_SYMMETRIC_BIT) << endl;
  
  return this->SymmetrizeAdAdResult(TmpState, coefficient);
}

// apply a^+_m_d operator to the state produced using Au or Ad method (without destroying it)
//
// m = first index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereWithSpinSzSymmetry::Add (int m, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m <<= 1;
  
  if ((TmpState & (0x1ul << m)) != 0)
    return this->TargetSpace->HilbertSpaceDimension;
//   int NewLzMax = this->ProdALzMax;
  
//   if (m > NewLzMax)
//     NewLzMax = m;
//   else
//     {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
//     }
  TmpState |= (0x1ul << m);
  return this->SymmetrizeAdAdResult(TmpState, coefficient);
}