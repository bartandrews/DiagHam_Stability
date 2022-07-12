////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of fermions on disk with restriction on the            //
//              number of reachable states or the number of fermions          //
//                                                                            //
//                        last modification : 01/02/2004                      //
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
#include "HilbertSpace/FermionOnDisk.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"

#include <math.h>

#ifdef __GSL__                                                                                                                                
#include <gsl/gsl_sf_gamma.h>                                                                                                                 
#endif                                                                                                                                        


using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constuctor
//

FermionOnDisk::FermionOnDisk()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = momentum total value
// lzMax = maximum angular momentum that a single particle can reach (negative if it has to be deduced from nbrFermions and totalLz)
// memory = amount of memory granted for precalculations

FermionOnDisk::FermionOnDisk (int nbrFermions, int totalLz, int lzMax, unsigned long memory)
{
  this->TargetSpace = this;
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
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrFermions, this->LzMax, this->LzMax, this->TotalLz, 0);
  this->ShiftedTotalLz = this->TotalLz;
  this->TotalLz = (this->TotalLz << 1) - (this->NbrFermions * this->LzMax);
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
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

FermionOnDisk::FermionOnDisk(const FermionOnDisk& fermions)
{
  this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->ShiftedTotalLz = fermions.ShiftedTotalLz;
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnDisk::~FermionOnDisk ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnDisk& FermionOnDisk::operator = (const FermionOnDisk& fermions)
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
  this->ShiftedTotalLz = fermions.ShiftedTotalLz;
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
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& FermionOnDisk::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  SqrtCoefficients[0] = 1.0;
  InvSqrtCoefficients[0] = 1.0;
  for (int k = 1; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt((double) k) * SqrtCoefficients[k - 1];
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  if (reference >= 0l)
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0 / state[reference];
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
   }
  delete[] SqrtCoefficients;
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& FermionOnDisk::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  InvSqrtCoefficients[0] = 1.0;
  SqrtCoefficients[0] = 1.0;
  for (int k = 1; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt((double) k) * SqrtCoefficients[k - 1];
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  this->ConvertFromUnnormalizedMonomialCore(state, SqrtCoefficients, InvSqrtCoefficients, reference, symmetryFactor);
  delete[] InvSqrtCoefficients;
  delete[] SqrtCoefficients;
  return state;
}

// convert a state such that its components are now expressed in the normalized basis, shifting all orbitals
//
// state = reference to the state to convert
// shift = shift to apply to each orbitals
// reference = set which component has been normalized to 1
// return value = converted state

RealVector& FermionOnDisk::ShiftedConvertFromUnnormalizedMonomial(RealVector& state, int shift, long reference)
{
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  InvSqrtCoefficients[0] = 1.0;
  SqrtCoefficients[0] = 1.0;
  for (int k = 1; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt((double) (shift + k)) * SqrtCoefficients[k - 1];
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  this->ConvertFromUnnormalizedMonomialCore(state, SqrtCoefficients, InvSqrtCoefficients, reference, false);
  delete[] InvSqrtCoefficients;
  delete[] SqrtCoefficients;
  return state;
}

// core part of the convertion of a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// sqrtCoefficients = array that contains the normalization coefficients
// invSqrtCoefficients = array that contains the inverts of the normalization coefficients
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& FermionOnDisk::ConvertFromUnnormalizedMonomialCore(RealVector& state, double* sqrtCoefficients, double* invSqrtCoefficients, long reference, bool symmetryFactor)
{
  if (reference >= 0l)
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0;
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
		  Coefficient *= invSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= sqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= invSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= sqrtCoefficients[TmpMonomialReference[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
      state /= state.Norm();
    }
  else
    {
      int* TmpMonomialReference = new int [this->NbrFermions];
      int* TmpMonomial = new int [this->NbrFermions];
      double Factor = 1.0;
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
		  Coefficient *= invSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrFermions) && (Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrFermions) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= sqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrFermions)
	    {
	      Coefficient *= invSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrFermions)
	    {
	      Coefficient *= sqrtCoefficients[TmpMonomialReference[Index2]];
	      ++Index2;
	    }
	  state[i] *= Coefficient;
	}
      delete[] TmpMonomialReference;
      delete[] TmpMonomial;
    }
  return state;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using real space partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// radius = radius of the A disk
// groundState = reference on the total system ground state
// shift = shift to apply to each orbitals
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix FermionOnDisk::EvaluatePartialDensityMatrixRealSpacePartition (int nbrFermionSector, int lzSector,  double radius , RealVector& groundState,int shift, AbstractArchitecture* architecture)
{
#ifdef __GSL__	
  double ReducedRadius = radius * radius / 2.0;
  double* IncompleteNormalizedGamma = new double[this->LzMax + 1];
  
  for (int i = 0; i <= this->LzMax ; i++)
    {
      IncompleteNormalizedGamma[i] = gsl_sf_gamma_inc_Q( ((double) i +shift) + 1.0 , ReducedRadius);
    }
  if (nbrFermionSector == 0)
    {
      if (lzSector == 0)
	{
	  unsigned long* TmpMonomial1 = new unsigned long [this->NbrFermions];
	  double TmpValue = 0.0;
	  for (int MinIndex = 0; MinIndex < this->HilbertSpaceDimension; ++MinIndex)    
	    {
	      this->ConvertToMonomial(this->StateDescription[MinIndex], TmpMonomial1);
	      double FormFactor = 0.0;
	      for (int i=0; i < this->NbrFermions; i++)
		FormFactor += log(1.0 - IncompleteNormalizedGamma[TmpMonomial1[i]]);
	      FormFactor = exp(FormFactor);
	      TmpValue += groundState[MinIndex] * groundState[MinIndex] * FormFactor;	
	    }
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue); 
	  
	  delete[] IncompleteNormalizedGamma;
	  delete[] TmpMonomial1;
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  
  if (nbrFermionSector == this->NbrFermions)
    {
      if (lzSector == this->ShiftedTotalLz)
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension, true);
	  unsigned long* TmpMonomial1 = new unsigned long [this->NbrFermions];
	  double* TmpStateCoefficient = new double [this->HilbertSpaceDimension];
	  for( int i = 0; i < this->HilbertSpaceDimension; i++)
	    {
	      TmpStateCoefficient[i] = 0.0;
	      this->ConvertToMonomial(this->StateDescription[i], TmpMonomial1);
	      for( int j=0; j<this->NbrFermions; j++)
		{
		  TmpStateCoefficient [i] += 0.5*log( IncompleteNormalizedGamma[TmpMonomial1[j]]);
		}
	      TmpStateCoefficient[i] = exp(TmpStateCoefficient[i]);
	    }
	  
	  for(int pos1 = 0; pos1 < this->HilbertSpaceDimension; pos1++)
	    for(int pos2 = pos1; pos2 < this->HilbertSpaceDimension; pos2++)
	      {
		TmpDensityMatrix.SetMatrixElement(pos1, pos2, groundState[pos1]*groundState[pos2]*TmpStateCoefficient[pos1]*TmpStateCoefficient[pos2]);
	      }
	  delete[] TmpMonomial1;
	  delete[] TmpStateCoefficient;
	  delete[] IncompleteNormalizedGamma;
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  
  int ComplementaryNbrFermionSector = this->NbrFermions - nbrFermionSector;
  
  if (nbrFermionSector == 1)
    {
      double TmpValue = 0.0;
			FermionOnDisk TmpHilbertSpace(this->NbrFermions - 1, this->ShiftedTotalLz - lzSector, this->LzMax);
      unsigned long TmpMask = 0x1ul << lzSector;
      double TmpStateCoefficient = IncompleteNormalizedGamma[lzSector];
      unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrFermionSector];
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex], TmpMonomial1);
	  double FormFactor = 0.0;
	  for (int i = 0; i < ComplementaryNbrFermionSector; i++)
	    FormFactor += log(1.0 - IncompleteNormalizedGamma[TmpMonomial1[i]]);
	  FormFactor = exp(FormFactor);
	  unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex];
	  if ((TmpState & TmpMask) == 0x0ul)
	    {
	      TmpState |= TmpMask;
	      int TmpLzMax = this->LzMax;
	      while ((TmpState >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  TmpValue += groundState[TmpPos] * groundState[TmpPos] *FormFactor*TmpStateCoefficient;	
		}
	    }
	}
      RealSymmetricMatrix TmpDensityMatrix(1);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
      delete[] TmpMonomial1;
      delete[] IncompleteNormalizedGamma;
      return TmpDensityMatrix;
    }
  
  
  FermionOnDisk TmpDestinationHilbertSpace(nbrFermionSector, lzSector, this->LzMax);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrFermionSector];
  unsigned long* TmpMonomial2 = new unsigned long [nbrFermionSector];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
	
  double* TmpStateCoefficient_Sign = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;
	
	if ( this->ShiftedTotalLz < lzSector)
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
    FermionOnDisk TmpHilbertSpace(ComplementaryNbrFermionSector, this->ShiftedTotalLz - lzSector, this->LzMax);
	
  //Compute the coefficients multiplying rhoA in TmpStateCoefficient
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      TmpStateCoefficient [i] = 0.0;
      TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.StateDescription[i], TmpMonomial2);
      for( int j=0; j<nbrFermionSector; j++)
	{
	  TmpStateCoefficient [i] += 0.5 * log( IncompleteNormalizedGamma[TmpMonomial2[j]] );
	}
      TmpStateCoefficient[i] = exp(TmpStateCoefficient[i]);
    }
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex];
		
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex], TmpMonomial1);
      double FormFactor = 0.0;
      for (int i=0; i < ComplementaryNbrFermionSector; i++)
	FormFactor += log(1.0 - IncompleteNormalizedGamma[TmpMonomial1[i]] );
      FormFactor = exp(FormFactor);
		
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpDestinationHilbertSpace.StateDescription[j];
	  if ((TmpState & TmpState2) == 0x0ul)
	    {
	      int TmpLzMax = this->LzMax;
	      unsigned long TmpState3 = TmpState | TmpState2;
	      while ((TmpState3 >> TmpLzMax) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
				
	      if (TmpPos != this->HilbertSpaceDimension)
		{
		  double Coefficient = 1.0;
		  unsigned long Sign = 0x0ul;
		  int Pos2 = TmpDestinationHilbertSpace.LzMax;
		  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
		    {
		      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
			--Pos2;
		      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
#ifdef  __64_BITS__
		      TmpState3 ^= TmpState3 >> 32;
#endif	
		      TmpState3 ^= TmpState3 >> 16;
		      TmpState3 ^= TmpState3 >> 8;
		      TmpState3 ^= TmpState3 >> 4;
		      TmpState3 ^= TmpState3 >> 2;
		      TmpState3 ^= TmpState3 >> 1;
		      Sign ^= TmpState3;
		      TmpState2 &= ~(0x1ul << Pos2);
		      --Pos2;
		    }
		  if ((Sign & 0x1ul) == 0x0ul)		  
		    Coefficient *= 1.0;
		  else
		    Coefficient *= -1.0;
		  TmpStatePosition[Pos] = TmpPos;
		  TmpStatePosition2[Pos] = j;
		  TmpStateCoefficient_Sign[Pos] = Coefficient;
		  ++Pos;
		}
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient_Sign[j] *TmpStateCoefficient[Pos2];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], FormFactor* TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient_Sign[k] *TmpStateCoefficient[TmpStatePosition2[k]]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  delete[] TmpStateCoefficient;
  delete[] TmpStateCoefficient_Sign;
  delete[] TmpMonomial1;
  delete[] TmpMonomial2;
  delete[] IncompleteNormalizedGamma;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
#else
  RealSymmetricMatrix TmpDensityMatrixZero;
  return TmpDensityMatrixZero;
#endif
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition. The entanglement matrix is only evaluated in a given Lz sector.
// and computed from precalculated particle entanglement matrix
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// radius = radius of the A disk
// shift = shift to apply to each orbitals
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& FermionOnDisk::EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, double radius, int shift, RealMatrix& entanglementMatrix)
{
#ifdef __GSL__
  double ReducedRadius = radius * radius / 2.0;
  double* IncompleteNormalizedGamma = new double[this->LzMax + 1];
  
  int ComplementaryNbrFermionSector = this->NbrFermions - nbrFermionSector;
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrFermionSector];
  unsigned long* TmpMonomial3 = new unsigned long [nbrFermionSector];
  
  
  for (int i = 0; i <= this->LzMax ; i++)
    {
      IncompleteNormalizedGamma[i] = gsl_sf_gamma_inc_Q( ((double) i + shift) + 1.0 , ReducedRadius);
    }
  
  if (nbrFermionSector == 0)
    {
      FermionOnDisk TmpHilbertSpace(this->NbrFermions, this->ShiftedTotalLz - lzSector , this->LzMax);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex] , TmpMonomial1);
	  double FormFactor = 0.0;
	  for (int i = 0; i < ComplementaryNbrFermionSector; i++)
	    FormFactor += log(1.0 - IncompleteNormalizedGamma[TmpMonomial1[i]]);
	  FormFactor = exp(0.5 * FormFactor);
	  entanglementMatrix(0, MinIndex) *= FormFactor; 
	}
      return entanglementMatrix;
    }
  
  if (ComplementaryNbrFermionSector == 0)
    {
      FermionOnDisk TmpDestinationHilbertSpace(nbrFermionSector, lzSector, this->LzMax);
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.StateDescription[i],  TmpMonomial3);
	  double Tmp = 0.0;     
	  for( int j = 0; j < nbrFermionSector; j++)
	    {
	      Tmp += log( IncompleteNormalizedGamma[TmpMonomial3[j]] );
	    }
	  Tmp = exp(0.5 * Tmp);        
	  entanglementMatrix(i, 0) *= Tmp;      
	}
      return entanglementMatrix;
    }
  
  FermionOnDisk TmpDestinationHilbertSpace(nbrFermionSector, lzSector, this->LzMax);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  int ComplementaryLz = (this->TotalLz + this->LzMax * this->NbrFermions)>>1;
  FermionOnDisk TmpHilbertSpace( ComplementaryNbrFermionSector, this->ShiftedTotalLz - lzSector, this->LzMax);
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.StateDescription[i],  TmpMonomial3);
      double Tmp = 0.0;     
      for( int j = 0; j < nbrFermionSector; j++)
	{
	  Tmp += log( IncompleteNormalizedGamma[TmpMonomial3[j]] );
	}
      Tmp = exp(0.5 * Tmp);
      for (int j = 0; j < TmpHilbertSpace.HilbertSpaceDimension; ++j)          
	entanglementMatrix(i, j) *= Tmp;      
    }
  
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.StateDescription[MinIndex] , TmpMonomial1);
      double FormFactor = 0.0;
      for (int i = 0; i < ComplementaryNbrFermionSector; i++)
	FormFactor += log(1.0 - IncompleteNormalizedGamma[TmpMonomial1[i]]);
      FormFactor = exp(0.5 * FormFactor);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	entanglementMatrix(j, MinIndex) *= FormFactor; 
    }
  
  delete[] TmpMonomial1;
  delete[] TmpMonomial3;
  
  return entanglementMatrix;
  
#else
  return entanglementMatrix;
#endif
}

