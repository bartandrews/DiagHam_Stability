/////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of fermions on sphere with spin with             //
//                                Lz<->-Lz symmetry                           //
//                                                                            //
//                        last modification : 15/08/2007                      //
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
#include "HilbertSpace/FermionOnSphereWithSpinAllSzLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
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

FermionOnSphereWithSpinAllSzLzSymmetry::FermionOnSphereWithSpinAllSzLzSymmetry ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// lzMax = twice the maximum Lz value reached by a fermion
// totalSz = twce the total spin value
// minusParity = select the  Lz <-> -Lz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinAllSzLzSymmetry::FermionOnSphereWithSpinAllSzLzSymmetry (int nbrFermions, int lzMax, bool minusParity, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  //this->TotalSpin = totalSz;
  //this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin) / 2;
  //this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin) / 2;
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

  // temporary space to accelerate state generation
  this->MaxTotalLz = new int*[2*NbrLzValue];
  for (int i=0; i<2*NbrLzValue; ++i)
    {
      MaxTotalLz[i] = new int[NbrFermions+1];
      for (int nbrFermions=0; nbrFermions<=NbrFermions; ++nbrFermions)
	{
	  MaxTotalLz[i][nbrFermions] = 0;
	  for (int n=0; (n<nbrFermions); ++n)	  
	    MaxTotalLz[i][nbrFermions] += (i-n)>>1;
 	}
    }

  this->HilbertSpaceDimension = (int) this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, (this->LzMax<<1)+1, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1);


  this->LzParitySign = 1.0;
  if (minusParity == true)
    this->LzParitySign = -1.0;

  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, (this->LzMax<<1)+1, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 0x0l);
  delete[] this->StateHighestBit;
  cout<<"Done with raw gen"<<endl;

  int TmpHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    if (this->GetCanonicalState(this->StateDescription[i]) != this->StateDescription[i])
      this->StateDescription[i] = 0x0ul;
    else
      {
	this->GetStateSymmetry(this->StateDescription[i]);
	if ((this->StateDescription[i] & FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT) == 0x0ul)
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
	     this->StateDescription[i] &= ~FERMION_SPHERE_SU2_LZ_SYMMETRIC_BIT;
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
  if (this->HilbertSpaceDimension > 0)
    {
      this->StateHighestBit =  new int [TmpHilbertSpaceDimension];
      this->GenerateLookUpTable(memory);
      delete[] this->StateHighestBit;
      this->StateHighestBit = 0;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	this->GetStateSymmetry(this->StateDescription[i]);

/*
//......................................................................................................
cout<<"DIm="<<this->HilbertSpaceDimension<<endl;
for (int i=0;i<this->HilbertSpaceDimension;i++)
	{
		cout<<"Vect "<<i<<" "; this->PrintState(cout,i); cout<<endl;
	}

cout<<"...................................................."<<endl;
for (int i=0;i<this->HilbertSpaceDimension;i++)
	{
		cout<<i<<"_____"; this->PrintState(cout,i); cout<<endl;
		cout<<"Acting with "<<endl;
		for (int j=0;j<=this->LzMax;j++)
			{
				double Coefficient;
				int ind=this->AduAd (i, j, j, Coefficient);
				if (ind>=this->HilbertSpaceDimension) cout<<"Warning"<<endl;
				cout<<j<<" outindex="<<ind<<" "; this->PrintState(cout,ind); cout<<endl;
			}	
	}




//......................................................................................................
*/

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
    }

#endif
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

  // clean up temporary space
  for (int i=0; i<2*NbrLzValue; ++i)
    delete [] MaxTotalLz[i];
  delete [] MaxTotalLz;


}

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinAllSzLzSymmetry::FermionOnSphereWithSpinAllSzLzSymmetry (char* fileName, unsigned long memory)
{
  this->ReadHilbertSpace(fileName);
  this->IncNbrFermions = this->NbrFermions + 1;
  //this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin) / 2;
  //this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin) / 2;
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
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

FermionOnSphereWithSpinAllSzLzSymmetry::FermionOnSphereWithSpinAllSzLzSymmetry(const FermionOnSphereWithSpinAllSzLzSymmetry& fermions)
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
  //this->TotalSpin = fermions.TotalSpin;
  //this->NbrFermionsUp = fermions.NbrFermionsUp;
 // this->NbrFermionsDown = fermions.NbrFermionsDown;
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnSphereWithSpinAllSzLzSymmetry::~FermionOnSphereWithSpinAllSzLzSymmetry ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinAllSzLzSymmetry& FermionOnSphereWithSpinAllSzLzSymmetry::operator = (const FermionOnSphereWithSpinAllSzLzSymmetry& fermions)
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
  //this->TotalSpin = fermions.TotalSpin;
  //this->NbrFermionsUp = fermions.NbrFermionsUp;
  //this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->LzParitySign = fermions.LzParitySign;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinAllSzLzSymmetry::Clone()
{
  return new FermionOnSphereWithSpinAllSzLzSymmetry(*this);
}


// convert a given state from symmetric basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector  

RealVector FermionOnSphereWithSpinAllSzLzSymmetry::ConvertToNbodyBasis(RealVector& state, FermionOnSphereWithSpinAllSz& nbodyBasis)
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
	      if ((((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0) && (this->LzParitySign > 0.0))
		  || (((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0) && (this->LzParitySign < 0.0)))
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
	  TmpVector[i] =  (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul))) * this->LzParitySign * state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2;
	}
    }
  return TmpVector;  
}

// convert a given state from the usual n-body basis to the symmetric basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereWithSpinAllSzLzSymmetry::ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphereWithSpinAllSz& nbodyBasis)
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
	  if (Signature != 0x0ul)	
	    TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += state[i] * M_SQRT1_2;
	  else
	    {
	      Signature = TmpState;
	      this->GetStateSingletParity(Signature);
	      if ((((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) == 0) && (this->LzParitySign > 0.0))
		  || (((Signature & FERMION_SPHERE_SU2_SINGLETPARITY_BIT) != 0) && (this->LzParitySign < 0.0)))
		TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = state[i];
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
	  TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += (1.0 - 2.0 * ((double) ((Signature >> FERMION_SPHERE_SU2_SINGLETPARITY_SHIFT) & 0x1ul))) * this->LzParitySign * state[i] * M_SQRT1_2;
	}
    }
  return TmpVector;  
}



// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereWithSpinAllSzLzSymmetry::AduAd (int index, int m1, int n2, double& Coefficient)
{

  this->ProdATemporaryState = this->StateDescription[index];
  n2 <<= 1;
  unsigned long TmpMask = (0x1ul << n2);
  if ((this->ProdATemporaryState & TmpMask) ^ TmpMask)
    return this->HilbertSpaceDimension;
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);

  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 1;
  ++m1;

  if ((TmpState & (0x1ul << m1)) != 0) 
    return this->HilbertSpaceDimension;

  Coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
  TmpState |= (0x1ul << m1);
  
  return this->SymmetrizeAdAdResult(TmpState, Coefficient);

}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereWithSpinAllSzLzSymmetry::AddAu (int index, int m1, int n2, double& Coefficient)
{

  this->ProdATemporaryState = this->StateDescription[index];
  n2 <<= 1;
  ++n2;

  unsigned long TmpMask = (0x1ul << n2);
  if ((this->ProdATemporaryState & TmpMask) ^ TmpMask)
    return this->HilbertSpaceDimension;
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SU2_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
  Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);

  unsigned long TmpState = this->ProdATemporaryState;
  m1 <<= 1;

  if ((TmpState & (0x1ul << m1)) != 0) 
    return this->HilbertSpaceDimension;

  Coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  Coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
  TmpState |= (0x1ul << m1);
  
  return this->SymmetrizeAdAdResult(TmpState, Coefficient);
}
  
// Project the state from the tunneling space (all Sz's)
// to the space with the fixed projection of Sz (given by SzValue)
//
// state = state that needs to be projected
// su2Space = the subspace onto which the projection is carried out
// SzValue = the desired value of Sz

RealVector FermionOnSphereWithSpinAllSzLzSymmetry::ForgeSU2FromTunneling(RealVector& state, FermionOnSphereWithSpinLzSymmetry& su2Space, int SzValue)
{
  RealVector FinalState(su2Space.GetHilbertSpaceDimension(), true);
  int counter=0;
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      unsigned long TmpState = this->StateDescription[j];
      TmpState= this->GetSignedCanonicalState(TmpState);
	 TmpState &= FERMION_SPHERE_SU2_SYMMETRIC_MASK;
      int TmpPos = this->LzMax << 1;
      int TmpSzValue=0;
      while (TmpPos >= 0)
	{
	  TmpSzValue+=((TmpState>>(TmpPos+1)) & 0x1ul) - ( (TmpState>>TmpPos) & 0x1ul );
	  TmpPos -= 2;
	}
      if (TmpSzValue == SzValue)
	{ 
       ++counter;
       int NewLzMax = 1 + (this->LzMax << 1);
       while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;
       FinalState[su2Space.FindStateIndex(TmpState, NewLzMax)] += state[j];
       cout<<"su2 "<<su2Space.FindStateIndex(TmpState, NewLzMax)<<" "<<TmpState<<endl;
	}
    }
  cout<<"Nbr of stored components = "<<counter<<endl;
  //FinalState /= FinalState.Norm();
  return FinalState;  
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
long FermionOnSphereWithSpinAllSzLzSymmetry::GenerateStates(int nbrFermions, int posMax, int totalLz, long pos)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (posMax < (nbrFermions - 1)))
    return pos;
  
//   if ((nbrFermions == 0) && (totalLz == 0)) // this could occur only if we assigned two fermions at once
//     {
//       cout << "0 fermions"<<endl;
//       this->StateDescription[pos] = 0x0ul;
//       this->StateHighestBit[pos++] = 0;
//       return pos;
//     }

  int LzTotalMax = this->MaxTotalLz[posMax][nbrFermions];  
//   int LzTotalMax = 0;
//   for (int n=0; n<nbrFermions; ++n) LzTotalMax += (posMax-n)/2;
  
  if (LzTotalMax < totalLz)
    {
      return pos;
    }

  if (nbrFermions == 1)
    {
      if ((posMax>>1) >= totalLz)
	{
	  if ( ((posMax>>1) > totalLz) || (((posMax>>1) == totalLz) && (posMax&1)))
	    {
	      this->StateHighestBit[pos] = (totalLz << 1) + 1;
	      this->StateDescription[pos++] = 0x1ul << ((totalLz << 1) + 1);
	    }
	  this->StateHighestBit[pos] = (totalLz << 1);
	  this->StateDescription[pos++] = 0x1ul << (totalLz << 1);
	}
      return pos;
    }

  if (((posMax>>1) == 0) && (totalLz != 0))
    return pos;
  
  long TmpPos;
  unsigned long Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, posMax - 1, totalLz - (posMax>>1),  pos);
  Mask = 0x1ul << posMax;
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, posMax - 1, totalLz, pos);
};


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// posMax = highest position for next particle to be placed
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereWithSpinAllSzLzSymmetry::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int posMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (posMax < (nbrFermions - 1)))
    return 0;
  
  int LzTotalMax = this->MaxTotalLz[posMax][nbrFermions];
//   int LzTotalMax = 0;
//   for (int n=0; n<nbrFermions; ++n) LzTotalMax += (posMax-n)>>1;
//   cout << "LzTotalMax="<<LzTotalMax<<", LzTotalMax2="<<LzTotalMax2<<endl;

  if (LzTotalMax < totalLz)
    return 0;
  if ((nbrFermions == 1) && ((posMax>>1) >= totalLz))
    {
      if ((posMax>>1) >totalLz) return 2;
      else return 1+(posMax&1);
    }
  return  (this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, posMax - 1, totalLz - (posMax>>1))
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, posMax - 1, totalLz));
}

