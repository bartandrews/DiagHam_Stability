////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of fermions on sphere using                     //
//                            the Lz <-> -Lz symmetry                         //
//                                                                            //
//                        last modification : 10/03/2007                      //
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
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/Endian.h"

#include <math.h>
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
// lzMax = twice the maximum Lz value reached by a fermion
// memory = amount of memory granted for precalculations

FermionOnSphereSymmetricBasis::FermionOnSphereSymmetricBasis (int nbrFermions, int lzMax, unsigned long memory)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrFermions, this->LzMax, this->LzMax, (this->NbrFermions * this->LzMax) >> 1, 0);
  int TmpHilbertSpaceDimension = 0;
   for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    if (this->GetCanonicalState(this->StateDescription[i]) != this->StateDescription[i])
      this->StateDescription[i] = 0x0ul;
    else
      ++TmpHilbertSpaceDimension;
  unsigned long* TmpStateDescription = new unsigned long [TmpHilbertSpaceDimension];
  int* TmpStateLzMax = new int [TmpHilbertSpaceDimension];
  TmpHilbertSpaceDimension = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    if (this->StateDescription[i] != 0x0ul)
      {
	TmpStateDescription[TmpHilbertSpaceDimension] = this->StateDescription[i];
	TmpStateLzMax[TmpHilbertSpaceDimension] = this->StateLzMax[i];
	++TmpHilbertSpaceDimension;
      }
  delete[] this->StateDescription;
  delete[] this->StateLzMax;
  this->StateDescription = TmpStateDescription;
  this->StateLzMax = TmpStateLzMax;
  this->HilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);
  delete[] this->StateLzMax;
  this->StateLzMax = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    this->GetStateSymmetry(this->StateDescription[i]);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * sizeof(unsigned long);
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
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

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

FermionOnSphereSymmetricBasis::FermionOnSphereSymmetricBasis (char* fileName, unsigned long memory)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      this->HilbertSpaceDimension = 0;
      return;
    }
  ReadLittleEndian(File, this->HilbertSpaceDimension);
  ReadLittleEndian(File, this->NbrFermions);
  ReadLittleEndian(File, this->LzMax);
  ReadLittleEndian(File, this->TotalLz);
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    ReadLittleEndian(File, this->StateDescription[i]);
  File.close();

  this->TargetSpace = this;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->Flag.Initialize();
#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;

  this->StateLzMax = new int [this->HilbertSpaceDimension];
  int TmpLzMax = this->LzMax;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->StateDescription[i] &= FERMION_SPHERE_SYMMETRIC_MASK;
      while ((this->StateDescription[i] >> TmpLzMax) == 0)
	--TmpLzMax;
      this->StateLzMax[i] = TmpLzMax;
    }
  this->GenerateLookUpTable(memory);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    this->GetStateSymmetry(this->StateDescription[i]);
  delete[] this->StateLzMax;
  this->StateLzMax = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

#ifdef __DEBUG__
  unsigned long UsedMemory = 0l;
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * sizeof(unsigned long);
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
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

FermionOnSphereSymmetricBasis::FermionOnSphereSymmetricBasis(const FermionOnSphereSymmetricBasis& fermions)
{
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->Flag = fermions.Flag;
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

FermionOnSphereSymmetricBasis::~FermionOnSphereSymmetricBasis ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
//       delete[] this->StateDescription;
//       delete[] this->StateLzMax;
//       delete[] this->SignLookUpTable;
//       delete[] this->SignLookUpTableMask;
//       delete[] this->LookUpTableShift;
//       for (int i = 0; i < this->NbrLzValue; ++i)
// 	delete[] this->LookUpTable[i];
//       delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereSymmetricBasis& FermionOnSphereSymmetricBasis::operator = (const FermionOnSphereSymmetricBasis& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateLzMax != 0)
	delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->InvertShift = fermions.InvertShift;
  this->InvertUnshift = fermions.InvertUnshift;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereSymmetricBasis::Clone()
{
  return new FermionOnSphereSymmetricBasis(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphereSymmetricBasis::WriteHilbertSpace (char* fileName)
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
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    WriteLittleEndian(File, this->StateDescription[i]);
  File.close();
  return true;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphereSymmetricBasis::ExtractSubspace (AbstractQuantumNumber& q, 
								      SubspaceSpaceConverter& converter)
{
  return 0;
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnSphereSymmetricBasis::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->TargetSpace = (FermionOnSphereSymmetricBasis*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int FermionOnSphereSymmetricBasis::GetTargetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}

// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id

int FermionOnSphereSymmetricBasis::GetHilbertSpaceAdditionalSymmetry()
{
  return ParticleOnSphere::LzMinusLzSymmetry;
}

// convert a given state from Lz-symmetric basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereSymmetricBasis::ConvertToNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      TmpState = this->GetSignedCanonicalState(nbodyBasis.StateDescription[i]);
      NewLzMax = this->LzMax;
      Signature = TmpState & FERMION_SPHERE_SYMMETRIC_BIT;
      TmpState &= FERMION_SPHERE_SYMMETRIC_MASK;
      while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;
      if (Signature != 0x0ul)	
	TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)] * M_SQRT1_2;
      else
	TmpVector[i] = state[this->FindStateIndex(TmpState, NewLzMax)];
    }
  return TmpVector;  
}

// convert a given state from the usual n-body basis to the Lz-symmetric basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector FermionOnSphereSymmetricBasis::ConvertToSymmetricNbodyBasis(RealVector& state, FermionOnSphere& nbodyBasis)
{
  RealVector TmpVector (this->GetHilbertSpaceDimension(), true);
  unsigned long TmpState;
  unsigned long Signature;  
  int NewLzMax;
  for (int i = 0; i < nbodyBasis.GetHilbertSpaceDimension(); ++i)
    {
      TmpState = this->GetSignedCanonicalState(nbodyBasis.StateDescription[i]);
      NewLzMax = this->LzMax;
      Signature = TmpState & FERMION_SPHERE_SYMMETRIC_BIT;
      TmpState &= FERMION_SPHERE_SYMMETRIC_MASK;
      while ((TmpState >> NewLzMax) == 0x0ul)
	--NewLzMax;
      if (Signature != 0x0ul)	
	TmpVector[this->FindStateIndex(TmpState, NewLzMax)] += state[i] * M_SQRT1_2;
      else
	TmpVector[this->FindStateIndex(TmpState, NewLzMax)] = state[i];
    }
  return TmpVector;  
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereSymmetricBasis::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  unsigned long State = this->StateDescription[index];
  unsigned long Signature = State & FERMION_SPHERE_SYMMETRIC_BIT;
  State &= FERMION_SPHERE_SYMMETRIC_MASK;
  if (((State & (((unsigned long) (0x1)) << n1)) == 0) 
      || ((State & (((unsigned long) (0x1)) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n2);
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n1);
  int NewLzMax = this->LzMax;
  while ((TmpState >> NewLzMax) == 0)
    --NewLzMax;
  if ((TmpState & (((unsigned long) (0x1)) << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewLzMax)
    {
      NewLzMax = m2;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m2);
  if ((TmpState & (((unsigned long) (0x1)) << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewLzMax)
    {
      NewLzMax = m1;
    }
  else
    {      
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m1);
  TmpState = this->GetSignedCanonicalState(TmpState);
  NewLzMax = this->LzMax;
  while (((TmpState & FERMION_SPHERE_SYMMETRIC_MASK) >> NewLzMax) == 0)
    --NewLzMax;
  int TmpIndex = this->FindStateIndex(TmpState, NewLzMax);
  if ((TmpState & FERMION_SPHERE_SYMMETRIC_BIT) != Signature)
    {
      if (Signature != 0)
	coefficient *= M_SQRT2;
      else
	coefficient *= M_SQRT1_2;
    }
  return TmpIndex;
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereSymmetricBasis::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  unsigned long State = this->StateDescription[index];
  unsigned long Signature = State & FERMION_SPHERE_SYMMETRIC_BIT;
  State &= FERMION_SPHERE_SYMMETRIC_MASK;
  --nbrIndices;
  for (int i = 0; i < nbrIndices; ++i)
    {
      if ((State & (((unsigned long) (0x1)) << n[i])) == 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      for (int j = i + 1; j <= nbrIndices; ++j)
	if ((n[i] == n[j]) || (m[i] == m[j]))
	  {
	    coefficient = 0.0;
	    return this->HilbertSpaceDimension; 	    
	  }
    }
  if ((State & (((unsigned long) (0x1)) << n[nbrIndices])) == 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }

  unsigned long TmpState = State;

  int Index;
  coefficient = 1.0;
  for (int i = nbrIndices; i >= 0; --i)
    {
      Index = n[i];
      coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      TmpState &= ~(((unsigned long) (0x1)) << Index);
    }
  int NewLzMax = this->LzMax;
  while ((TmpState >> NewLzMax) == 0)
    --NewLzMax;  
  for (int i = nbrIndices; i >= 0; --i)
    {
      Index = m[i];
      if ((TmpState & (((unsigned long) (0x1)) << Index))!= 0)
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
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (((unsigned long) (0x1)) << Index);
    }
  TmpState = this->GetSignedCanonicalState(TmpState);
  NewLzMax = this->LzMax;
  while (((TmpState & FERMION_SPHERE_SYMMETRIC_MASK) >> NewLzMax) == 0)
    --NewLzMax;
  int TmpIndex = this->FindStateIndex(TmpState, NewLzMax);
  if ((TmpState & FERMION_SPHERE_SYMMETRIC_BIT) != Signature)
    {
      if (Signature != 0)
	coefficient *= M_SQRT2;
      else
	coefficient *= M_SQRT1_2;
    }
  return TmpIndex;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereSymmetricBasis::ProdA (int index, int* n, int nbrIndices)
{
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SYMMETRIC_MASK;
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = n[i];
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
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;

  return Coefficient;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereSymmetricBasis::AA (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  unsigned long TmpMask = (0x1l << n1) | (0x1l << n2);
  if ((this->ProdATemporaryState & TmpMask) ^ TmpMask)
    return 0.0;
  this->ProdASignature = this->ProdATemporaryState & FERMION_SPHERE_SYMMETRIC_BIT;
  this->ProdATemporaryState &= FERMION_SPHERE_SYMMETRIC_MASK;
  this->ProdATemporaryState &= ~TmpMask;
  TmpMask = this->ProdATemporaryState;
  TmpMask &= ((0x1l << n1) - 1ul);
  TmpMask &= ~((0x1l << n2) - 1ul);
#ifdef  __64_BITS__
  TmpMask ^= TmpMask >> 32;
#endif
  TmpMask ^= TmpMask >> 16;
  TmpMask ^= TmpMask >> 8;
  TmpMask ^= TmpMask >> 4;
  TmpMask ^= TmpMask >> 2;
  TmpMask ^= TmpMask << 1;
  TmpMask &= 0x2ul;
  return (1.0 - ((double) TmpMask));
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereSymmetricBasis::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = m[i];
      if ((TmpState & (0x1l << Index)) != 0)
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
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (0x1l << Index);
    }

  TmpState = this->GetSignedCanonicalState(TmpState);
  NewLzMax = this->LzMax;
  while (((TmpState & FERMION_SPHERE_SYMMETRIC_MASK) >> NewLzMax) == 0)
    --NewLzMax;
  int TmpIndex = this->FindStateIndex(TmpState, NewLzMax);
  if ((TmpState & FERMION_SPHERE_SYMMETRIC_BIT) != this->ProdASignature)
    {
      if (this->ProdASignature != 0)
	coefficient *= M_SQRT2;
      else
	coefficient *= M_SQRT1_2;
    }
  return TmpIndex;
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereSymmetricBasis::AdAd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  unsigned long TmpMask = (0x1l << m1) | (0x1l << m2);
  if (TmpState & TmpMask)
    return  this->HilbertSpaceDimension;
  TmpState |= TmpMask;
  TmpMask = TmpState;
  TmpMask &= ((0x1l << m1) - 1ul);
  TmpMask &= ~((0x1l << m2) - 1ul);
#ifdef  __64_BITS__
  TmpMask ^= TmpMask >> 32;
#endif
  TmpMask ^= TmpMask >> 16;
  TmpMask ^= TmpMask >> 8;
  TmpMask ^= TmpMask >> 4;
  TmpMask ^= TmpMask >> 2;
  TmpMask ^= TmpMask << 1;
  TmpMask &= 0x2ul;
  coefficient = (1.0 - ((double) TmpMask));
  TmpState = this->GetSignedCanonicalState(TmpState);
  int NewLzMax = this->LzMax;
  TmpMask = TmpState & FERMION_SPHERE_SYMMETRIC_MASK;
  while ((TmpMask >> NewLzMax) == 0)
    --NewLzMax;
  NewLzMax = this->FindStateIndex(TmpState, NewLzMax);
  if ((TmpState & FERMION_SPHERE_SYMMETRIC_BIT) != this->ProdASignature)
    {
      if (this->ProdASignature != 0)
	coefficient *= M_SQRT2;
      else
	coefficient *= M_SQRT1_2;
    }
   return NewLzMax;
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereSymmetricBasis::AdA (int index, int m)
{
  if ((this->StateDescription[index] & (((unsigned long) (0x1)) << m)) != 0)
    return 1.0;
  else
    return 0.0;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
  
int FermionOnSphereSymmetricBasis::AdA (int index, int m, int n, double& coefficient)
{
  unsigned long State = this->StateDescription[index];
  unsigned long Signature = State & FERMION_SPHERE_SYMMETRIC_BIT;
  State &= FERMION_SPHERE_SYMMETRIC_MASK;
  if ((State & (((unsigned long) (0x1)) << n)) == 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n);
  int NewLzMax = this->LzMax;
  while ((TmpState >> NewLzMax) == 0)
    --NewLzMax;
  if ((TmpState & (((unsigned long) (0x1)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLzMax)
    {
      NewLzMax = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m);
  TmpState = this->GetSignedCanonicalState(TmpState);
  NewLzMax = this->LzMax;
  while (((TmpState & FERMION_SPHERE_SYMMETRIC_MASK) >> NewLzMax) == 0)
    --NewLzMax;
  int TmpIndex = this->FindStateIndex(TmpState, NewLzMax);
  if ((TmpState & FERMION_SPHERE_SYMMETRIC_BIT) != Signature)
    {
      if (Signature != 0)
	coefficient *= M_SQRT2;
      else
	coefficient *= M_SQRT1_2;
    }
  return TmpIndex;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereSymmetricBasis::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->NbrLzValue; ++i)
    Str << ((TmpState >> i) & ((unsigned long) 0x1)) << " ";
  return Str;
}


// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphereSymmetricBasis::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
							     int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
  ComplexMatrix Slater(this->NbrFermions, this->NbrFermions);
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
	  if ((TmpStateDescription & ((unsigned long) 1)) == ((unsigned long) 1))
	    {
	      Indices[Pos] = Lz;
	      ++Pos;
	    }
	  ++Lz;
	  TmpStateDescription >>= 1;
	}
      for (int i = 0; i < this->NbrFermions; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < this->NbrFermions; ++j)
	    {
	      Slater[i].Re(j) = TmpColum2.Re(Indices[j]);
	      Slater[i].Im(j) = TmpColum2.Im(Indices[j]);
	    }
	}
      Complex SlaterDet = Slater.Determinant();
      Value += SlaterDet * (state[k] * Factor);
    }
  delete[] Indices;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereSymmetricBasis::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
  
// compute particule-hole symmetric state from a given state
//
// state = vector corresponding to the state to symmetrize
// holeBasis = n-body basis on which the symmetrized state has to be expressed

RealVector  FermionOnSphereSymmetricBasis::ParticleHoleSymmetrize (RealVector& state, FermionOnSphere& holeBasis)
{
  RealVector TmpVector(holeBasis.HilbertSpaceDimension, true);
  unsigned long TmpMask = (0x1ul << (this->LzMax + 1)) - 1;
  FermionOnSphereSymmetricBasis& TmpHoleBasis = (FermionOnSphereSymmetricBasis&) holeBasis;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = TmpHoleBasis.GetSignedCanonicalState((~this->StateDescription[i]) & TmpMask);
      int TmpLzMax = this->LzMax;
      while ((TmpState & (0x1ul << TmpLzMax)) == 0x0l)
	--TmpLzMax;
      TmpHoleBasis.GetSignedCanonicalState(TmpState);
      TmpVector[TmpHoleBasis.FindStateIndex(TmpState, TmpLzMax)] = state[i];
    }
  return TmpVector;
}

