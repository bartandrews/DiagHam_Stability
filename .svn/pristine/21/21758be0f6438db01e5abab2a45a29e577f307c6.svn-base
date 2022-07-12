////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of fermions on sphere                       //
//                                                                            //
//                        last modification : 24/06/2002                      //
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
#include "HilbertSpace/FermionOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/StringTools.h"
#include "MathTools/FactorialCoefficient.h" 

#include <math.h>
#include <stdlib.h>
#include <fstream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constuctor
//

FermionOnSphere::FermionOnSphere()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// memory = amount of memory granted for precalculations

FermionOnSphere::FermionOnSphere (int nbrFermions, int totalLz, int lzMax, unsigned long memory)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrFermions, this->LzMax, this->LzMax, (this->TotalLz + this->NbrFermions * this->LzMax) >> 1, 0);
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphere::FermionOnSphere(const FermionOnSphere& fermions)
{
  this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
}

// destructor
//

FermionOnSphere::~FermionOnSphere ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateLzMax != 0)
	delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      if (this->LookUpTableShift != 0)
	{
	  delete[] this->LookUpTableShift;
	  for (int i = 0; i < this->NbrLzValue; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;
	}
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphere& FermionOnSphere::operator = (const FermionOnSphere& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
  if (this->TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->InitializeWaveFunctionEvaluation();
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphere::Clone()
{
  return new FermionOnSphere(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphere::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);
  WriteLittleEndian(File, this->NbrFermions);
  WriteLittleEndian(File, this->LzMax);
  WriteLittleEndian(File, this->TotalLz);
  if (this->HilbertSpaceDimension != 0)
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateDescription[i]);
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateLzMax[i]);
    }
  else
    {
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateDescription[i]);
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->StateLzMax[i]);
    }
  File.close();
  return true;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnSphere::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnSphere::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphere::ExtractSubspace (AbstractQuantumNumber& q, 
							SubspaceSpaceConverter& converter)
{
  return 0;
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnSphere::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->TargetSpace = (FermionOnSphere*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int FermionOnSphere::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
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

int FermionOnSphere::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  unsigned long State = this->StateDescription[index];
  if ((n1 > StateLzMax) || (n2 > StateLzMax) || ((State & (((unsigned long) (0x1)) << n1)) == 0) 
      || ((State & (((unsigned long) (0x1)) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = StateLzMax;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n2);
  if (NewLzMax == n2)
    while ((TmpState >> NewLzMax) == 0)
      --NewLzMax;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n1);
  if (NewLzMax == n1)
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
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphere::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  unsigned long State = this->StateDescription[index];
  --nbrIndices;
  for (int i = 0; i < nbrIndices; ++i)
    {
      if ((n[i] > StateLzMax) || ((State & (((unsigned long) (0x1)) << n[i])) == 0))
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
  if (n[nbrIndices] > StateLzMax)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }

  int NewLzMax = StateLzMax;
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
      if (NewLzMax == Index)
	while ((TmpState >> NewLzMax) == 0)
	  --NewLzMax;
    }
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
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphere::ProdA (int index, int* n, int nbrIndices)
{
  this->ProdALzMax = this->StateLzMax[index];
  this->ProdATemporaryState = this->StateDescription[index];
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

double FermionOnSphere::AA (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];

  if (((ProdATemporaryState & (((unsigned long) (0x1)) << n1)) == 0) 
      || ((ProdATemporaryState & (((unsigned long) (0x1)) << n2)) == 0) || (n1 == n2))
    return 0.0;

  this->ProdALzMax = this->StateLzMax[index];

  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n1);

  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphere::ProdAd (int* m, int nbrIndices, double& coefficient)
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
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphere::AdAd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if ((TmpState & (((unsigned long) (0x1)) << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
  if (m2 > NewLzMax)
    NewLzMax = m2;
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
    NewLzMax = m1;
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
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphere::AdA (int index, int m)
{
  if ((this->StateDescription[index] & (0x1ul << m)) != 0ul)
    return 1.0;
  else
    return 0.0;
}


// get the variance of the state
// index = index of state to consider
int FermionOnSphere::StateVariance (int index)
{
  int MyStateLzMax = this->StateLzMax[index];
  unsigned long State = this->StateDescription[index];
  int var=0;
  for (int m=0; m<=MyStateLzMax; ++m)
    if  ( (State & (0x1ul << m)) != 0x0ul)
      var += ((m<<1)-LzMax)*((m<<1)-LzMax);
  return var;
}


// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

// attention: check sign returned by this function!
int FermionOnSphere::AdA (int index, int m, int n, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  unsigned long State = this->StateDescription[index];
  if ((n > StateLzMax) || ((State & (((unsigned long) (0x1)) << n)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = StateLzMax;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n);
  if (NewLzMax == n)
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
  return this->TargetSpace->FindStateIndex(TmpState, NewLzMax);
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphere::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
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

// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int FermionOnSphere::FindStateIndex(char* stateDescription)
{
  char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != (this->LzMax + 1))
    return -1;
  unsigned long TmpState = 0x0ul;
  int TmpNbrParticles = 0;
  int TmpTotalLz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      int Tmp = atoi(TmpDescription[i]);
      TmpState |= ((unsigned long)  Tmp) << i;
      TmpTotalLz += (i * Tmp);
      TmpNbrParticles += Tmp;
      delete[] TmpDescription[i];
    }
  delete[] TmpDescription;
  if ((TmpNbrParticles != this->NbrFermions) || (TmpTotalLz != ((this->TotalLz + this->NbrFermions * this->LzMax) >> 1)))
    return -1;
  int TmpLzMax = this->LzMax;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  return this->FindStateIndex(TmpState, TmpLzMax);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphere::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->NbrLzValue; ++i)
    Str << ((TmpState >> i) & ((unsigned long) 0x1)) << " ";
//  Str << " key = " << this->Keys[state] << " lzmax position = " << this->LzMaxPosition[Max * (this->NbrFermions + 1) + TmpState[Max]]
//  Str << " position = " << this->FindStateIndex(TmpState, this->StateLzMax[state]);
//  if (state !=  this->FindStateIndex(TmpState, this->StateLzMax[state]))
//    Str << " error! ";
  return Str;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphere::PrintStateMonomial (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  Str << "[";
  int i = this->LzMax;
  while (((TmpState >> i) & 0x1ul) == 0x0ul)
    --i;
  Str << i;
  --i;
  for (; i >=0; --i)
    if (((TmpState >> i) & 0x1ul) != 0x0ul)
      Str << "," << i;
  Str << "]";
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// currentLzMax = momentum maximum value for fermions that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int FermionOnSphere::GenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int pos)
{
  if ((nbrFermions == 0) || (totalLz < 0) || (currentLzMax < (nbrFermions - 1)))
    return pos;
  int LzTotalMax = ((2 * currentLzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (totalLz > LzTotalMax)
    return pos;
  if ((nbrFermions == 1) && (currentLzMax >= totalLz))
    {
      this->StateDescription[pos] = ((unsigned long) 0x1) << totalLz;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }
  if (LzTotalMax == totalLz)
    {
      unsigned long Mask = 0;
      for (int i = currentLzMax - nbrFermions + 1; i <= currentLzMax; ++i)
	Mask |= (((unsigned long) 1) << i);
      this->StateDescription[pos] = Mask;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }

  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, lzMax, ReducedCurrentLzMax, totalLz - currentLzMax, pos);
  unsigned long Mask = ((unsigned long) 1) << currentLzMax;
  for (int i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (lzMax == currentLzMax)
    return this->GenerateStates(nbrFermions, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, TmpPos);
  else
    return this->GenerateStates(nbrFermions, lzMax, ReducedCurrentLzMax, totalLz, TmpPos);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphere::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->NbrLzValue);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrLzValue)
    this->MaximumLookUpShift = this->NbrLzValue;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrLzValue];
  this->LookUpTableShift = new int [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentLzMax = this->StateLzMax[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentLzMax];
  if (CurrentLzMax < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentLzMax] = 0;
  else
    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentLzMax];
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
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  /*	  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
	    cout << TmpLookUpTable[j] << " ";
	    cout << endl << "-------------------------------------------" << endl;*/
 	  CurrentLzMax = this->StateLzMax[i];
	  TmpLookUpTable = this->LookUpTable[CurrentLzMax];
	  if (CurrentLzMax < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentLzMax] = 0;
	  else
	    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentLzMax];
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
//	      CurrentLookUpTableValue = TmpLookUpTableValue;
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
  /*  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
    cout << TmpLookUpTable[j] << " ";
    cout << endl << "-------------------------------------------" << endl;*/
  
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
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#endif
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphere::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + nbrFermions * lzMax) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long FermionOnSphere::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)))
    return 0l;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return 0l;
  if ((nbrFermions == 1) && (lzMax >= totalLz))
    return 1l;
  if (LzTotalMax == totalLz)
    return 1l;
  return  (this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz));
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphere::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
					       int firstComponent, int nbrComponent)
{
  // fields used by EvaluateWaveFunction
#ifdef __LAPACK__
  ComplexLapackDeterminant Slater(this->NbrFermions);
#else
  ComplexMatrix Slater(this->NbrFermions, this->NbrFermions);
#endif
  ComplexMatrix Functions(this->LzMax + 1, this->NbrFermions);
  
  // temporary array used to stored indices when evaluating wave function
  int* Indices = new int [this->NbrFermions];
  
  Complex Value;
  Complex Tmp;
  RealVector TmpCoordinates(2);
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
#ifdef __LAPACK__
	      Slater.SetMatrixElement(i,j,TmpColum2.Re(Indices[j]), TmpColum2.Im(Indices[j]));
#else
	      Slater[i].Re(j) = TmpColum2.Re(Indices[j]);
	      Slater[i].Im(j) = TmpColum2.Im(Indices[j]);
#endif
	    }
	}
      //cout << Slater << endl;

      // can calculate with lapack for a regular ComplexMatrix by Complex SlaterDet = Slater.LapackDeterminant();

      Complex SlaterDet = Slater.Determinant();
      Value += SlaterDet * (state[k] * Factor);
    }
  delete [] Indices;
  return Value;
}

// evaluate wave functions in real space using a given basis and only for agiven range of components
//
// states = array of vector corresponding to the state in the Fock basis
// nbrStates = number of states in the states array
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// waveFuntions = array where the  wave function values at the given location will be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate

void FermionOnSphere::EvaluateWaveFunctions (RealVector* states, int nbrStates, RealVector& position, AbstractFunctionBasis& basis,
					     Complex* waveFuntions, int firstComponent, int nbrComponent)
{
  // fields used by EvaluateWaveFunction
#ifdef __LAPACK__
  ComplexLapackDeterminant Slater(this->NbrFermions);
#else
  ComplexMatrix Slater(this->NbrFermions, this->NbrFermions);
#endif
  ComplexMatrix Functions(this->LzMax + 1, this->NbrFermions);
  
  // temporary array used to stored indices when evaluating wave function
  int* Indices = new int [this->NbrFermions];
  
  for (int i = 0; i < nbrStates; ++i)
    waveFuntions[i] = 0.0;
  Complex Tmp;
  RealVector TmpCoordinates(2);
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
#ifdef __LAPACK__
	      Slater.SetMatrixElement(i,j,TmpColum2.Re(Indices[j]), TmpColum2.Im(Indices[j]));
#else
	      Slater[i].Re(j) = TmpColum2.Re(Indices[j]);
	      Slater[i].Im(j) = TmpColum2.Im(Indices[j]);
#endif
	    }
	}
      //cout << Slater << endl;

      // can calculate with lapack for a regular ComplexMatrix by Complex SlaterDet = Slater.LapackDeterminant();

      Complex SlaterDet = Slater.Determinant();
      SlaterDet *= Factor;
      for (int i = 0; i < nbrStates; ++i)
	waveFuntions[i] += SlaterDet * states[i][k];
    }
  delete [] Indices;
}
                               
  
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphere::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}

  
// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  FermionOnSphere::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState)
{  
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrFermionSector == 0))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrFermionSector == this->NbrFermions))
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * groundState[j]);
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;  
	}
    }

  unsigned long TmpMask = (((0x1ul << (this->LzMax + 2)) - 1) >> subsytemSize) << subsytemSize;
  unsigned long TmpSubsystemMask = (0x1ul << subsytemSize) - 1;
  int TmpIndex;
  int ShiftedTotalLz = (this->TotalLz + this->NbrFermions * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrFermionSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int TmpStateMaxLz = ShiftedLzComplementarySector - (((NbrFermionsComplementarySector - 2 + (subsytemSize << 1)) * (NbrFermionsComplementarySector - 1)) >> 1);
  int MinIndex = 0;
  int MaxIndex = this->HilbertSpaceDimension - 1;
  if ((NbrFermionsComplementarySector > 0) && ((NbrFermionsComplementarySector + subsytemSize - 2) > this->StateLzMax[MaxIndex]))
    MaxIndex = this->LookUpTable[NbrFermionsComplementarySector + subsytemSize - 2][0];
  if ((TmpStateMaxLz < this->StateLzMax[0]) && ((TmpStateMaxLz + 1) >  this->StateLzMax[MaxIndex]) && (TmpStateMaxLz >= subsytemSize))
    MinIndex = this->LookUpTable[TmpStateMaxLz + 1][0];
  
  unsigned long TmpComplementarySubsystem;
  int TmpNbrFermions;
  int TmpTotalLz;
  int TmpNbrOne[] = {  
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
  int TmpSumOccupation[] = {
    0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6,
    4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 10, 10,
    5, 5, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 10, 10, 11, 11,
    9, 9, 10, 10, 11, 11, 12, 12, 12, 12, 13, 13, 14, 14, 15, 15,
    6, 6, 7, 7, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 12, 12,
    10, 10, 11, 11, 12, 12, 13, 13, 13, 13, 14, 14, 15, 15, 16, 16,
    11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17,
    15, 15, 16, 16, 17, 17, 18, 18, 18, 18, 19, 19, 20, 20, 21, 21,
    7, 7, 8, 8, 9, 9, 10, 10, 10, 10, 11, 11, 12, 12, 13, 13,
    11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17,
    12, 12, 13, 13, 14, 14, 15, 15, 15, 15, 16, 16, 17, 17, 18, 18,
    16, 16, 17, 17, 18, 18, 19, 19, 19, 19, 20, 20, 21, 21, 22, 22,
    13, 13, 14, 14, 15, 15, 16, 16, 16, 16, 17, 17, 18, 18, 19, 19,
    17, 17, 18, 18, 19, 19, 20, 20, 20, 20, 21, 21, 22, 22, 23, 23,
    18, 18, 19, 19, 20, 20, 21, 21, 21, 21, 22, 22, 23, 23, 24, 24,
    22, 22, 23, 23, 24, 24, 25, 25, 25, 25, 26, 26, 27, 27, 28, 28};
  int TmpPartialNbrOne;
  if (nbrFermionSector <= 1)
    {
      unsigned long Key = 0x0ul;
      if (nbrFermionSector == 1)
	Key = 0x1ul << ShiftedLzSector;
      double TmpValue = 0.0;
      while (MinIndex <= MaxIndex)
	{
	  if ((this->StateDescription[MinIndex] & TmpSubsystemMask) == Key)
	    TmpValue += groundState[MinIndex] * groundState[MinIndex];	    
	  ++MinIndex;
	}
      RealSymmetricMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  else
    {
      FermionOnSphere TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
      RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      long TmpNbrNonZeroElements = 0;
      while (MinIndex <= MaxIndex)
	{
	  TmpIndex = MinIndex;
	  TmpComplementarySubsystem = this->StateDescription[TmpIndex] & TmpMask;
	  ++TmpIndex;
	  while ((TmpIndex <= MaxIndex) && ((this->StateDescription[TmpIndex] & TmpMask) == TmpComplementarySubsystem))
	    ++TmpIndex;
	  TmpPartialNbrOne = TmpNbrOne[TmpComplementarySubsystem & 0xffl];
	  TmpNbrFermions = TmpPartialNbrOne;
	  TmpTotalLz = TmpSumOccupation[TmpComplementarySubsystem & 0xffl];
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 8) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 8) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 3;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 16) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 16) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 4;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 24) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 24) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 24;
#ifdef  __64_BITS__
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 32) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 32) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 5;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 40) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 40) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 40;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 48) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 48) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 48;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 56) & 0xffl];      
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 56) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 56;
#endif
	  if ((TmpNbrFermions == NbrFermionsComplementarySector) && (ShiftedLzComplementarySector == TmpTotalLz))
	    {
	      int Pos = 0;
	      for (int i = MinIndex; i < TmpIndex; ++i)
		{
		  unsigned long TmpState = this->StateDescription[i] & TmpSubsystemMask;
		  int TmpLzMax = subsytemSize - 1;
		  while ((TmpState & (0x1ul << TmpLzMax)) == 0x0ul)
		    --TmpLzMax;
		  TmpStatePosition[Pos] = TmpDestinationHilbertSpace.FindStateIndex(TmpState, TmpLzMax);
		  ++Pos;
		}
	      int Pos2;
	      Pos = 0;
	      int Pos3;
	      double TmpValue;
	      for (int i = MinIndex; i < TmpIndex; ++i)
		{
		  Pos2 = 0;
		  Pos3 = TmpStatePosition[Pos];
		  TmpValue = groundState[i];
		  for (int j = MinIndex; j < TmpIndex; ++j)
		    {
		      if (Pos3 <=  TmpStatePosition[Pos2])
			{
			  TmpDensityMatrix.AddToMatrixElement(Pos3, TmpStatePosition[Pos2], TmpValue * groundState[j]);
			  ++TmpNbrNonZeroElements;
			}
		      ++Pos2;
		    }
		  ++Pos;
		}
	    }
	  MinIndex = TmpIndex;
	}
      delete[] TmpStatePosition;
      if (TmpNbrNonZeroElements > 0)	
	return TmpDensityMatrix;
      else
	{
	  RealSymmetricMatrix TmpDensityMatrixZero;
	  return TmpDensityMatrixZero;
	}
    }
}

// evaluate a density matrix of a shited subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// subsystemShift = shift (in number of states) of the subsytem with repect to the leftmost state (i.e -LzMax)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem (return a zero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  FermionOnSphere::EvaluateShiftedPartialDensityMatrix (int subsytemSize, int subsystemShift, int nbrFermionSector, int lzSector, RealVector& groundState)
{  
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrFermionSector == 0))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrFermionSector == this->NbrFermions))
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * groundState[j]);
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  unsigned long TmpMask = ((0x1ul << (this->LzMax + 2)) - 1);
  unsigned long TmpSubsystemMask = ((0x1ul << subsytemSize) - 1) << subsystemShift;
  TmpMask &= ~TmpSubsystemMask;
  int TmpIndex;
  int ShiftedTotalLz = (this->TotalLz + this->NbrFermions * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrFermionSector * (subsytemSize + (subsystemShift << 1) - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  //  int TmpStateMaxLz = ShiftedLzComplementarySector - (((NbrFermionsComplementarySector - 2 + (subsytemSize << 1)) * (NbrFermionsComplementarySector - 1)) >> 1);
  int MinIndex = 0;
  int MaxIndex = this->HilbertSpaceDimension - 1;
//   if ((NbrFermionsComplementarySector > 0) && ((NbrFermionsComplementarySector + subsytemSize - 2) > this->StateLzMax[MaxIndex]))
//     MaxIndex = this->LookUpTable[NbrFermionsComplementarySector + subsytemSize - 2][0];
//   if ((TmpStateMaxLz < this->StateLzMax[0]) && ((TmpStateMaxLz + 1) >  this->StateLzMax[MaxIndex]) && (TmpStateMaxLz >= subsytemSize))
//     MinIndex = this->LookUpTable[TmpStateMaxLz + 1][0];
  
  unsigned long TmpComplementarySubsystem;
  int TmpNbrFermions;
  int TmpTotalLz;
  int TmpNbrOne[] = {  
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
  int TmpSumOccupation[] = {
    0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6,
    4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 10, 10,
    5, 5, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 10, 10, 11, 11,
    9, 9, 10, 10, 11, 11, 12, 12, 12, 12, 13, 13, 14, 14, 15, 15,
    6, 6, 7, 7, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 12, 12,
    10, 10, 11, 11, 12, 12, 13, 13, 13, 13, 14, 14, 15, 15, 16, 16,
    11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17,
    15, 15, 16, 16, 17, 17, 18, 18, 18, 18, 19, 19, 20, 20, 21, 21,
    7, 7, 8, 8, 9, 9, 10, 10, 10, 10, 11, 11, 12, 12, 13, 13,
    11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17,
    12, 12, 13, 13, 14, 14, 15, 15, 15, 15, 16, 16, 17, 17, 18, 18,
    16, 16, 17, 17, 18, 18, 19, 19, 19, 19, 20, 20, 21, 21, 22, 22,
    13, 13, 14, 14, 15, 15, 16, 16, 16, 16, 17, 17, 18, 18, 19, 19,
    17, 17, 18, 18, 19, 19, 20, 20, 20, 20, 21, 21, 22, 22, 23, 23,
    18, 18, 19, 19, 20, 20, 21, 21, 21, 21, 22, 22, 23, 23, 24, 24,
    22, 22, 23, 23, 24, 24, 25, 25, 25, 25, 26, 26, 27, 27, 28, 28};
  int TmpPartialNbrOne;
  if (nbrFermionSector <= 1)
    {
      unsigned long Key = 0x0ul;
      if (nbrFermionSector == 1)
	Key = 0x1ul << ShiftedLzSector;
      double TmpValue = 0.0;
      while (MinIndex <= MaxIndex)
	{
	  if ((this->StateDescription[MinIndex] & TmpSubsystemMask) == Key)
	    TmpValue += groundState[MinIndex] * groundState[MinIndex];	    
	  ++MinIndex;
	}
      RealSymmetricMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  else
    {
      FermionOnSphere TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
      RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      long TmpNbrNonZeroElements = 0;
      while (MinIndex <= MaxIndex)
	{
	  TmpIndex = MinIndex;
	  TmpComplementarySubsystem = this->StateDescription[TmpIndex] & TmpMask;
	  ++TmpIndex;
	  while ((TmpIndex <= MaxIndex) && ((this->StateDescription[TmpIndex] & TmpMask) == TmpComplementarySubsystem))
	    ++TmpIndex;
	  TmpPartialNbrOne = TmpNbrOne[TmpComplementarySubsystem & 0xffl];
	  TmpNbrFermions = TmpPartialNbrOne;
	  TmpTotalLz = TmpSumOccupation[TmpComplementarySubsystem & 0xffl];
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 8) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 8) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 3;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 16) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 16) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 4;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 24) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 24) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 24;
#ifdef  __64_BITS__
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 32) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 32) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 5;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 40) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 40) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 40;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 48) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 48) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 48;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 56) & 0xffl];      
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 56) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 56;
#endif
	  if ((TmpNbrFermions == NbrFermionsComplementarySector) && (ShiftedLzComplementarySector == TmpTotalLz))
	    {
	      int Pos = 0;
	      for (int i = MinIndex; i < TmpIndex; ++i)
		{
		  unsigned long TmpState = (this->StateDescription[i] & TmpSubsystemMask) >> subsystemShift;
		  int TmpLzMax = subsytemSize - 1;
		  while ((TmpState & (0x1ul << TmpLzMax)) == 0x0ul)
		    --TmpLzMax;
		  TmpStatePosition[Pos] = TmpDestinationHilbertSpace.FindStateIndex(TmpState, TmpLzMax);
		  ++Pos;
		}
	      int Pos2;
	      Pos = 0;
	      int Pos3;
	      double TmpValue;
	      for (int i = MinIndex; i < TmpIndex; ++i)
		{
		  Pos2 = 0;
		  Pos3 = TmpStatePosition[Pos];
		  TmpValue = groundState[i];
		  for (int j = MinIndex; j < TmpIndex; ++j)
		    {
		      if (Pos3 <=  TmpStatePosition[Pos2])
			{
			  TmpDensityMatrix.AddToMatrixElement(Pos3, TmpStatePosition[Pos2], TmpValue * groundState[j]);
			  ++TmpNbrNonZeroElements;
			}
		      ++Pos2;
		    }
		  ++Pos;
		}
	    }
	  MinIndex = TmpIndex;
	}
      delete[] TmpStatePosition;
      if (TmpNbrNonZeroElements > 0)	
	return TmpDensityMatrix;
      else
	{
	  RealSymmetricMatrix TmpDensityMatrixZero;
	  return TmpDensityMatrixZero;
	}
    }
}

// compute part of the Schmidt decomposition, allowing cut in the reduced denisty matrix eigenvalue space
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// eigenvalueCut = discard all contribution from the reduced density matrix whose eigenvalues is lower than eigenvalueCut
// rebuiltSchmidtGroundState = reference on the state to whose current sector contribution to the Schmidt decomposition will be added 
// diagonalizedDensityMatrix = reference on the diagonalized reduced density matrix associated to the current sector (with down ordered diagonal elements)
// transformationMatrix =  reference on the transformation matric that diagonalizes the reduced density matrix
// return value = reference on rebuiltSchmidtGroundState

RealVector& FermionOnSphere::EvaluatePartialSchmidtDecomposition(int subsytemSize, int nbrFermionSector, int lzSector, double eigenvalueCut,
								 RealVector& groundState, RealVector& rebuiltSchmidtGroundState,
								 RealDiagonalMatrix& diagonalizedDensityMatrix, RealMatrix& transformationMatrix)
{
  if (subsytemSize <= 0)
    {
      return rebuiltSchmidtGroundState;
    }
  if (subsytemSize > this->LzMax)
    {
      for (long i = 0l; i < groundState.GetLargeVectorDimension(); ++i)
	rebuiltSchmidtGroundState[i] = groundState[i];
      return rebuiltSchmidtGroundState;
    }
  int NbrKeptEigenvalues = 0;  
  for (int i = 0; i < diagonalizedDensityMatrix.GetNbrRow(); ++i)
    if (diagonalizedDensityMatrix[i] >=  eigenvalueCut)
      ++NbrKeptEigenvalues;
  if (NbrKeptEigenvalues == 0)
    return rebuiltSchmidtGroundState;

  unsigned long TmpMask = (((0x1ul << (this->LzMax + 2)) - 1) >> subsytemSize) << subsytemSize;
  unsigned long TmpSubsystemMask = (0x1ul << subsytemSize) - 1;
  int TmpIndex;
  int ShiftedTotalLz = (this->TotalLz + this->NbrFermions * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrFermionSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrFermionsComplementarySector = this->NbrFermions - nbrFermionSector;
  int TmpStateMaxLz = ShiftedLzComplementarySector - (((NbrFermionsComplementarySector - 2 + (subsytemSize << 1)) * (NbrFermionsComplementarySector - 1)) >> 1);
  int MinIndex = 0;
  int MaxIndex = this->HilbertSpaceDimension - 1;
  if ((NbrFermionsComplementarySector > 0) && ((NbrFermionsComplementarySector + subsytemSize - 2) > this->StateLzMax[MaxIndex]))
    MaxIndex = this->LookUpTable[NbrFermionsComplementarySector + subsytemSize - 2][0];
  if ((TmpStateMaxLz < this->StateLzMax[0]) && ((TmpStateMaxLz + 1) >  this->StateLzMax[MaxIndex]) && (TmpStateMaxLz >= subsytemSize))
    MinIndex = this->LookUpTable[TmpStateMaxLz + 1][0];
  
  unsigned long TmpComplementarySubsystem;
  int TmpNbrFermions;
  int TmpTotalLz;
  int TmpNbrOne[] = {  
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
  int TmpSumOccupation[] = {
    0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6,
    4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 9, 9, 10, 10,
    5, 5, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 10, 10, 11, 11,
    9, 9, 10, 10, 11, 11, 12, 12, 12, 12, 13, 13, 14, 14, 15, 15,
    6, 6, 7, 7, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 12, 12,
    10, 10, 11, 11, 12, 12, 13, 13, 13, 13, 14, 14, 15, 15, 16, 16,
    11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17,
    15, 15, 16, 16, 17, 17, 18, 18, 18, 18, 19, 19, 20, 20, 21, 21,
    7, 7, 8, 8, 9, 9, 10, 10, 10, 10, 11, 11, 12, 12, 13, 13,
    11, 11, 12, 12, 13, 13, 14, 14, 14, 14, 15, 15, 16, 16, 17, 17,
    12, 12, 13, 13, 14, 14, 15, 15, 15, 15, 16, 16, 17, 17, 18, 18,
    16, 16, 17, 17, 18, 18, 19, 19, 19, 19, 20, 20, 21, 21, 22, 22,
    13, 13, 14, 14, 15, 15, 16, 16, 16, 16, 17, 17, 18, 18, 19, 19,
    17, 17, 18, 18, 19, 19, 20, 20, 20, 20, 21, 21, 22, 22, 23, 23,
    18, 18, 19, 19, 20, 20, 21, 21, 21, 21, 22, 22, 23, 23, 24, 24,
    22, 22, 23, 23, 24, 24, 25, 25, 25, 25, 26, 26, 27, 27, 28, 28};
  int TmpPartialNbrOne;
  
  if (nbrFermionSector <= 1)
    {
      unsigned long Key = 0x0ul;
      if (nbrFermionSector == 1)
	Key = 0x1ul << ShiftedLzSector;
      while (MinIndex <= MaxIndex)
	{
	  if ((this->StateDescription[MinIndex] & TmpSubsystemMask) == Key)
	    rebuiltSchmidtGroundState[MinIndex] += groundState[MinIndex];
	  ++MinIndex;
	}
      return rebuiltSchmidtGroundState;
    }
  else
    {
      FermionOnSphere TmpDestinationHilbertSpace(nbrFermionSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
      RealMatrix TmpMatrix (TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension);
      for (int i = 0; i < TmpMatrix.GetNbrRow(); ++i)
	for (int j = 0; j < TmpMatrix.GetNbrRow(); ++j)
	  {
	    double Tmp = 0.0;
	    for (int k = 0; k < NbrKeptEigenvalues; ++k)
	      Tmp += transformationMatrix(i, k) * transformationMatrix(j, k);
	    TmpMatrix(i, j) = Tmp;
	  }
	
      while (MinIndex <= MaxIndex)
	{
	  TmpIndex = MinIndex;
	  TmpComplementarySubsystem = this->StateDescription[TmpIndex] & TmpMask;
	  ++TmpIndex;
	  while ((TmpIndex <= MaxIndex) && ((this->StateDescription[TmpIndex] & TmpMask) == TmpComplementarySubsystem))
	    ++TmpIndex;
	  TmpPartialNbrOne = TmpNbrOne[TmpComplementarySubsystem & 0xffl];
	  TmpNbrFermions = TmpPartialNbrOne;
	  TmpTotalLz = TmpSumOccupation[TmpComplementarySubsystem & 0xffl];
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 8) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 8) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 3;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 16) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 16) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 4;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 24) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 24) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 24;
#ifdef  __64_BITS__
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 32) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 32) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne << 5;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 40) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 40) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 40;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 48) & 0xffl];
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 48) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 48;
	  TmpPartialNbrOne = TmpNbrOne[(TmpComplementarySubsystem >> 56) & 0xffl];      
	  TmpNbrFermions += TmpPartialNbrOne;
	  TmpTotalLz += TmpSumOccupation[(TmpComplementarySubsystem >> 56) & 0xffl];
	  TmpTotalLz += TmpPartialNbrOne * 56;
#endif
	  if ((TmpNbrFermions == NbrFermionsComplementarySector) && (ShiftedLzComplementarySector == TmpTotalLz))
	    {
	      int Pos = 0;
	      for (int i = MinIndex; i < TmpIndex; ++i)
		{
		  unsigned long TmpState = this->StateDescription[i] & TmpSubsystemMask;
		  int TmpLzMax = subsytemSize - 1;
		  while ((TmpState & (0x1ul << TmpLzMax)) == 0x0ul)
		    --TmpLzMax;
		  TmpStatePosition[Pos] = TmpDestinationHilbertSpace.FindStateIndex(TmpState, TmpLzMax);
		  ++Pos;
		}
	      int Pos2;
	      Pos = 0;
	      double TmpValue;
	      for (int i = MinIndex; i < TmpIndex; ++i)
		{
		  Pos2 = 0;
		  double Tmp = 0.0;
		  for (int j = MinIndex; j < TmpIndex; ++j)
		    {
		      Tmp += groundState[j] * TmpMatrix(TmpStatePosition[Pos], TmpStatePosition[Pos2]);
		      ++Pos2;
		    }
		  rebuiltSchmidtGroundState[i] = Tmp;
		  ++Pos;
		}
	    }
	  MinIndex = TmpIndex;
	}
//       for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
// 	{
// 	  unsigned long TmpPattern = TmpDestinationHilbertSpace.StateDescription[i];
// 	  double TmpNorm = 0.0;
// 	  int TmpIndex = 0;
// 	  for (int j = MinIndex; j <= MaxIndex; ++j)
// 	    {
// 	      if ((this->StateDescription[j] & TmpSubsystemMask) == TmpPattern)
// 		{
// 		  TmpNorm += groundState[j] * groundState[j];
// 		  TmpStatePosition[TmpIndex] = j;
// 		  ++TmpIndex;
// 		}
// 	    }
// 	  if ((TmpIndex > 0) && (TmpNorm > MACHINE_PRECISION))
// 	    {
// 	      TmpNorm = 1.0 / sqrt(TmpNorm);
// 	      double Sign = 0.0;
// 	      for (int j = 0; j <  NbrKeptEigenvalues; ++j)
// 		{
// 		  double TmpEigenvalue =  diagonalizedDensityMatrix[j];
// 		  double Tmp = (sqrt(TmpEigenvalue) * TmpNorm) * transformationMatrix(j, i);
// 		  double Tmp2 = groundState[TmpStatePosition[0]];
// 		  if (TmpEigenvalue >  eigenvalueCut)
// 		    Sign += Tmp * Tmp2;
// 		}
// 	      if ((Sign * groundState[TmpStatePosition[0]]) >= 0.0)
// 		Sign = 1.0;
// 	      else
// 		Sign = -1.0;
// 	      for (int j = 0; j <  NbrKeptEigenvalues; ++j)
// 		{
// 		  double TmpEigenvalue =  diagonalizedDensityMatrix[j];
// 		  if (TmpEigenvalue >  eigenvalueCut)
// 		    {
// 		      TmpEigenvalue = Sign * (sqrt(TmpEigenvalue) * TmpNorm) * transformationMatrix(j, i);
// 		      for (int k = 0; k < TmpIndex; ++k)
// 			rebuiltSchmidtGroundState[TmpStatePosition[k]] += TmpEigenvalue * groundState[TmpStatePosition[k]];
// 		    }
// 		}
// 	    }
// 	}

      delete[] TmpStatePosition;
    }
  return rebuiltSchmidtGroundState;
}

// compute particule-hole symmetric state from a given state
//
// state = vector corresponding to the state to symmetrize
// holeBasis = n-body basis on which the symmetrized state has to be expressed

RealVector FermionOnSphere::ParticleHoleSymmetrize (RealVector& state, FermionOnSphere& holeBasis)
{
  RealVector TmpVector(holeBasis.HilbertSpaceDimension, true);
  unsigned long TmpMask = (0x1ul << (this->LzMax + 1)) - 1;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = (~this->StateDescription[i]) & TmpMask;
      int TmpLzMax = this->LzMax;
      while ((TmpState & (0x1ul << TmpLzMax)) == 0x0l)
	--TmpLzMax;
      TmpVector[holeBasis.FindStateIndex(TmpState, TmpLzMax)] = state[i];
    }
  return TmpVector;
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// return value = converted state

RealVector& FermionOnSphere::ConvertToUnnormalizedMonomial(RealVector& state, long reference)
{
  double* SqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  if (reference >= 0l)
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0 / state[reference];
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	  InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
	}
      unsigned long TmpState = this->StateDescription[reference];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & 1ul) != 0ul)
	  TmpMonomialReference[Index++] = j;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      TmpMonomial[Index++] = j;
	  int Index1 = 0;
	  int Index2 = 0;
	  double Coefficient = Factor;
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	    {
	      while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomialReference[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      state[reference] = 1.0;
      delete[] InvSqrtCoefficients;
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
   }
  else
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0;
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	  InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
	}
      unsigned long TmpState = this->StateDescription[0l];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & 1ul) != 0ul)
	  TmpMonomialReference[Index++] = j;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      TmpMonomial[Index++] = j;
	  int Index1 = 0;
	  int Index2 = 0;
	  double Coefficient = Factor;
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	    {
	      while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomialReference[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      delete[] InvSqrtCoefficients;
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
//       for (int k = 0; k <= this->LzMax; ++k)
// 	SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k) * ((double) (this->LzMax + 1)));
//       unsigned long TmpState;
//       int Index = 0;
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	{
// 	  Index = 0;
// 	  TmpState = this->StateDescription[i];
// 	  double Coefficient = 1.0;
// 	  for (int j = this->LzMax; j >= 0; --j)
// 	    if (((TmpState >> j) & 1ul) != 0ul)
// 	      Coefficient /= SqrtCoefficients[j];
// 	  state[i] *= Coefficient;
// 	}
   }
  delete[] SqrtCoefficients;
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// return value = converted state

RealVector& FermionOnSphere::ConvertFromUnnormalizedMonomial(RealVector& state, long reference)
{
  BinomialCoefficients Binomials(this->LzMax);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  if (reference >= 0l)
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0 / state[reference];
      state[reference] = 1.0;
      double* SqrtCoefficients = new double [this->LzMax + 1];
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	  SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
	}
      unsigned long TmpState = this->StateDescription[reference];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & 1ul) != 0ul)
	  TmpMonomialReference[Index++] = j;
      for (int i = 1; i < this->HilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      TmpMonomial[Index++] = j;
	  int Index1 = 0;
	  int Index2 = 0;
	  double Coefficient = Factor;
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	    {
	      while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomialReference[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
      delete[] InvSqrtCoefficients;
      state /= state.Norm();
    }
  else
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0;
      double* SqrtCoefficients = new double [this->LzMax + 1];
      double* InvSqrtCoefficients = new double [this->LzMax + 1];
      for (int k = 0; k <= this->LzMax; ++k)
	{
	  InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
	  SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
	}
      unsigned long TmpState = this->StateDescription[0l];
      int Index = 0;
      for (int j = this->LzMax; j >= 0; --j)
	if (((TmpState >> j) & 1ul) != 0ul)
	  TmpMonomialReference[Index++] = j;
      for (int i = 1; i < this->HilbertSpaceDimension; ++i)
	{
	  Index = 0;
	  TmpState = this->StateDescription[i];
	  for (int j = this->LzMax; j >= 0; --j)
	    if (((TmpState >> j) & 1ul) != 0ul)
	      TmpMonomial[Index++] = j;
	  int Index1 = 0;
	  int Index2 = 0;
	  double Coefficient = Factor;
	  while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions))
	    {
	      while ((Index1 < this->NbrFermions) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomialReference[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
      delete[] InvSqrtCoefficients;
//       for (int k = 0; k <= this->LzMax; ++k)
// 	SqrtCoefficients[k] = 1.0 / sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k)* ((double) (this->LzMax + 1)));
//       unsigned long TmpState;
//       int Index = 0;
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	{
// 	  Index = 0;
// 	  TmpState = this->StateDescription[i];
// 	  double Coefficient = 1.0;
// 	  for (int j = this->LzMax; j >= 0; --j)
// 	    if (((TmpState >> j) & 1ul) != 0ul)
// 	      Coefficient /= SqrtCoefficients[j];
// 	  state[i] *= Coefficient;
// 	}
    }
  delete[] SqrtCoefficients;
  return state;
}

// fuse two states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// leftVector = reference on the vector whose Hilbert space will be fuse to the left
// rightVector = reference on the vector whose Hilbert space will be fuse to the right
// padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
// leftSpace = point to the Hilbert space that will be fuse to the left
// rightSpace = point to the Hilbert space that will be fuse to the right
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// return value = reference on the fused state

RealVector& FermionOnSphere::FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
				 ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace, bool symmetrizedFlag)
{
  FermionOnSphere* LeftSpace = (FermionOnSphere*) leftSpace;
  FermionOnSphere* RightSpace = (FermionOnSphere*) rightSpace;
  int StateShift = RightSpace->LzMax + padding + 1;
  for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long TmpState1 = LeftSpace->StateDescription[i] << StateShift;
      double Coefficient = leftVector[i];
      int TmpLzMax = this->LzMax;
      while ((TmpState1 >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      if (symmetrizedFlag == false)
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState2 = RightSpace->StateDescription[j];
	      TmpState2 |= TmpState1;
	      double Coefficient2 = Coefficient;
	      Coefficient2 *= rightVector[j];	  
	      int TmpIndex = this->FindStateIndex(TmpState2, TmpLzMax);
	      outputVector[TmpIndex] = Coefficient2;
	    }
	}
      else
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState2 = RightSpace->StateDescription[j];
	      TmpState2 |= TmpState1;
	      double Coefficient2 = Coefficient;
	      Coefficient2 *= rightVector[j];	  
	      int TmpIndex = this->FindStateIndex(TmpState2, TmpLzMax);
	      outputVector[TmpIndex] = Coefficient2;
	      unsigned long TmpState3 = this->GetSymmetricState(TmpState2);
	      if (TmpState3 != TmpState2)
		{
		  int TmpLzMax2 = this->LzMax;
		  while ((TmpState3 >> TmpLzMax2) == 0x0ul)
		    --TmpLzMax2;
		  TmpIndex = this->FindStateIndex(TmpState3, TmpLzMax2);
		  outputVector[TmpIndex] = Coefficient2;      
		}
	    }
	}
    }
  return outputVector;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

double FermionOnSphere::JackSqrNormalization (RealVector& outputVector, long minIndex, long nbrComponents)
{
  double SqrNorm = 0.0;
  int* TmpMonomial = new int [this->NbrFermions];
  FactorialCoefficient Factorial;
  int HalfLzMax = this->LzMax >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  for (long i = minIndex; i < MaxIndex; ++i)
    {
      Factorial.SetToOne();
      this->ConvertToMonomial(this->StateDescription[i], TmpMonomial);
      for (int k = 0; k < this->NbrFermions; ++k)
	{
	  if (HalfLzMax < TmpMonomial[k])
	    Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomial[k]);
	  else
	    if (HalfLzMax > TmpMonomial[k])
	      Factorial.PartialFactorialDivide(TmpMonomial[k] + 1, HalfLzMax);
	}	      
      SqrNorm +=(outputVector[i] * outputVector[i]) * Factorial.GetNumericalValue();
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  cout << endl;
  delete[] TmpMonomial;
  return SqrNorm;
}

// remove part of each Fock state, discarding component if the Fock state does not a given pattern
//
// inputVector = state to truncate
// reducedSpace = Hilbert space where the truncated state will lie
// pattern = array describing the pattern 
// patternSize = pattern size
// patternShift = indicate where the pattern has to be applied
// return value = trucated state

RealVector FermionOnSphere::TruncateStateWithPatternConstraint(RealVector& inputVector, ParticleOnSphere* reducedSpace, int* pattern, int patternSize, int patternShift)
{
  FermionOnSphere* TmpSpace = (FermionOnSphere*) reducedSpace; 
  unsigned long PatternMask =  (0x1ul << patternSize) - 0x1ul;
  PatternMask <<= patternShift;
  unsigned long Pattern = 0x0ul;
  for (int i = 0; i < patternSize; ++i)
    {
      if (pattern[i] == 1)
	Pattern |= 0x1ul << i;
    }
  Pattern <<= patternShift;
  unsigned long Mask = (0x1ul << patternShift) - 0x1ul;
  RealVector TmpVector (TmpSpace->LargeHilbertSpaceDimension, true);
  unsigned long Tmp;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {      
      Tmp = this->StateDescription[i]; 
      if ((Tmp & PatternMask) == Pattern)
	{	  	  
	  Tmp = ((Tmp >> patternSize) & (~Mask)) | (Tmp & Mask);
	  int TmpLzMax = TmpSpace->LzMax;
	  while (((Tmp >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpIndex = TmpSpace->FindStateIndex(Tmp, TmpLzMax);
	  if (TmpIndex < TmpSpace->HilbertSpaceDimension)
	    TmpVector[TmpIndex] = inputVector[i];
	}
    }
  return TmpVector;
}
 
// remove part of each Fock state, discarding component if the Fock state does not a given pattern
//
// inputVector = state to truncate
// reducedSpace = Hilbert space where the truncated state will lie
// pattern = array describing the pattern 
// patternSize = pattern size
// patternShift = indicate where the pattern has to be applied
// return value = trucated state

// RealVector FermionOnSphere::TruncateStateWithPatternConstraint(RealVector& inputVector, ParticleOnSphere* reducedSpace, int* pattern, int patternSize, int patternShift)
// {
//   FermionOnSphere* TmpSpace = (FermionOnSphere*) reducedSpace; 
//   unsigned long PatternMask =  (0x1ul << patternSize) - 0x1ul;
//   PatternMask <<= patternShift;
//   unsigned long Pattern = 0x0ul;
//   for (int i = 0; i < patternSize; ++i)
//     {
//       if (pattern[i] == 1)
// 	Pattern |= 0x1ul << i;
//     }
//   Pattern <<= patternShift;
//   unsigned long Mask = (0x1ul << patternShift) - 0x1ul;
//   RealVector TmpVector (TmpSpace->LargeHilbertSpaceDimension, true);
//   unsigned long Tmp;
//   for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//     {      
//       Tmp = this->StateDescription[i]; 
//       if ((Tmp & PatternMask) == Pattern)
// 	{	  	  
// 	  Tmp = ((Tmp >> patternSize) & (~Mask)) | (Tmp & Mask);
// 	  int TmpLzMax = TmpSpace->LzMax;
// 	  while (((Tmp >> TmpLzMax) & 0x1ul) == 0x0ul)
// 	    --TmpLzMax;
// 	  int TmpIndex = TmpSpace->FindStateIndex(Tmp, TmpLzMax);
// 	  if (TmpIndex < TmpSpace->HilbertSpaceDimension)
// 	    TmpVector[TmpIndex] = inputVector[i];
// 	}
//     }
//   return TmpVector;
// }
 
