////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of bosons on sphere  using the Haldane basis           //
//                            for system size such that                       //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 08/07/2008                      //
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
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/ArrayTools.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h" 
#include "GeneralTools/Endian.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHESphereJackGeneratorOperation.h"

#include <math.h>
#include <fstream>
#include <sys/time.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::fstream;
using std::ios;


// default constructor
//

BosonOnSphereHaldaneHugeBasisShort::BosonOnSphereHaldaneHugeBasisShort ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson
// maxFileSize = maximum file size (in MBytes)
// referenceState = array that describes the reference state to start from
// memory = amount of memory granted for precalculations
// referenceState = array that describes the reference state to start from
// symmetricFlag = indicate if a symmetric basis has to be used (only available if totalLz = 0)
// fullDimension = provide the full (i.e. without squeezing) Hilbert space dimension (0 if it has to be computed)

BosonOnSphereHaldaneHugeBasisShort::BosonOnSphereHaldaneHugeBasisShort (int nbrBosons, int totalLz, int lzMax, unsigned long maxFileSize, int* referenceState, unsigned long memory, bool symmetricFlag, long fullDimension)
{
  this->TargetSpace = this;
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  int ShiftedLzMax = this->LzMax + nbrBosons - 1;
  int* TmpReferenceState = new int[ShiftedLzMax + 1];
  int TmpIndex = 0;
  for (int i = 0; i <= ShiftedLzMax; ++i)
    TmpReferenceState[i] = 0;   
  for (int i = 0; i <= this->LzMax; ++i)
    {
      for (int j = 0; j < referenceState[i]; ++j)
	{
	  TmpReferenceState[TmpIndex] = 1;
	  ++TmpIndex;
	}
      ++TmpIndex;
    }
  this->FermionBasis = 0;
  this->FermionHugeBasis = new FermionOnSphereHaldaneHugeBasis(nbrBosons, totalLz, this->LzMax + nbrBosons - 1, maxFileSize, TmpReferenceState, memory);
  delete[] TmpReferenceState;
  this->HilbertSpaceDimension = this->FermionHugeBasis->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->FermionHugeBasis->GetLargeHilbertSpaceDimension();

  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->TemporaryMonomial = new unsigned long [this->NbrBosons];
  this->TemporaryMonomial2 = new unsigned long [this->NbrBosons];

  this->Flag.Initialize();
  int TmpLzMax = this->LzMax;


  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;

}

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memoryHilbert = amount of memory granted to store the Hilbert space (in Mbytes)

BosonOnSphereHaldaneHugeBasisShort::BosonOnSphereHaldaneHugeBasisShort (char* fileName, unsigned long memoryHilbert)
{
  timeval TotalStartingTime;
  gettimeofday (&(TotalStartingTime), 0);

  this->FermionBasis = 0;
  this->FermionHugeBasis = new FermionOnSphereHaldaneHugeBasis(fileName, memoryHilbert);
  this->HilbertSpaceDimension = this->FermionHugeBasis->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->FermionHugeBasis->GetLargeHilbertSpaceDimension();
  
  timeval TotalEndingTime;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
  cout << "Hilbert space factorization done in " << Dt << "s" <<endl;
 
//  gettimeofday (&(TotalStartingTime), 0);

//   ifstream FileHilbert;
//   FileHilbert.open(this->FermionHugeBasis->HilbertSpaceFileName, ios::binary | ios::in);
//   FileHilbert.seekg (this->FermionHugeBasis->FileHeaderSize, ios::beg);
//   unsigned long CurrentPartition;
//   for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//     {
//       ReadLittleEndian(FileHilbert, CurrentPartition);
//       if (i != this->FermionHugeBasis->FindStateIndexFactorized(CurrentPartition))
// //      if (CurrentPartition != this->FermionHugeBasis->GetStateFactorized(i))
// 	cout << i << " " << this->FermionHugeBasis->FindStateIndexFactorized(CurrentPartition) << " " << hex << CurrentPartition << " " << this->FermionHugeBasis->GetStateFactorized(i) << dec << endl;
//     }
//   FileHilbert.close();

//   long Test = 5;
//   long TmpShift = 470;
//   long TmpStep = 5;
//   unsigned long* TmpStateArray = new unsigned long [Test];
//   long* TmpIndexArray = new long[Test];
//   for (long i = 0; i < Test; ++i)
//     TmpStateArray[i] = this->FermionHugeBasis->GetStateFactorized((i * TmpStep) + TmpShift);
//   this->FermionHugeBasis->FindMultipleStateIndexFactorized(TmpStateArray, Test, TmpIndexArray);
//   for (long i = 0; i < Test; ++i)
//     if (((i * TmpStep) + TmpShift) != TmpIndexArray[i])
//       cout << ((i * TmpStep) + TmpShift) << " : " << TmpIndexArray[i] << " " << endl;
 
//   gettimeofday (&(TotalEndingTime), 0);
//   Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
//     ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
//   cout << "Hilbert space consistency check done in " << Dt << "s" <<endl;
//  exit(0) ;

  this->NbrBosons = this->FermionHugeBasis->NbrFermions;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = this->FermionHugeBasis->TotalLz;
  this->LzMax = this->FermionHugeBasis->LzMax - this->NbrBosons + 1;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;

  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->TemporaryMonomial = new unsigned long [this->NbrBosons];
  this->TemporaryMonomial2 = new unsigned long [this->NbrBosons];
  this->Flag.Initialize();
  int TmpLzMax = this->LzMax;

  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;
}


// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSphereHaldaneHugeBasisShort::BosonOnSphereHaldaneHugeBasisShort(const BosonOnSphereHaldaneHugeBasisShort& bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = 0;
  this->FermionHugeBasis = (FermionOnSphereHaldaneHugeBasis*) bosons.FermionHugeBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
  this->TemporaryMonomial = new unsigned long [this->NbrBosons];
  this->TemporaryMonomial2 = new unsigned long [this->NbrBosons];
}

// destructor
//

BosonOnSphereHaldaneHugeBasisShort::~BosonOnSphereHaldaneHugeBasisShort ()
{
  if (this->FermionHugeBasis != 0)
    delete this->FermionHugeBasis;
  delete[] this->TemporaryMonomial;
  delete[] this->TemporaryMonomial2;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereHaldaneHugeBasisShort& BosonOnSphereHaldaneHugeBasisShort::operator = (const BosonOnSphereHaldaneHugeBasisShort& bosons)
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = 0;
  this->FermionHugeBasis = (FermionOnSphereHaldaneHugeBasis*) bosons.FermionHugeBasis->Clone();
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->TemporaryMonomial = new unsigned long [this->NbrBosons];
  this->TemporaryMonomial2 = new unsigned long [this->NbrBosons];

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereHaldaneHugeBasisShort::Clone()
{
  return new BosonOnSphereHaldaneHugeBasisShort(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool BosonOnSphereHaldaneHugeBasisShort::WriteHilbertSpace (char* fileName)
{
  return this->FermionHugeBasis->WriteHilbertSpace(fileName);
}

// convert a given state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector BosonOnSphereHaldaneHugeBasisShort::ConvertToNbodyBasis(RealVector& state, BosonOnSphereHaldaneHugeBasisShort& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetLargeHilbertSpaceDimension(), true);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      TmpVector[nbodyBasis.FermionHugeBasis->FindStateIndexFactorized(this->FermionHugeBasis->GetStateFactorized(i))] = state[i];
    }
  return TmpVector;
}

// convert a given state from the usual n-body basis to the Haldane basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector BosonOnSphereHaldaneHugeBasisShort::ConvertFromNbodyBasis(RealVector& state, BosonOnSphereHaldaneHugeBasisShort& nbodyBasis)
{
  RealVector TmpVector (this->LargeHilbertSpaceDimension, true);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    TmpVector[i] = state[nbodyBasis.FermionHugeBasis->FindStateIndexFactorized(this->FermionHugeBasis->GetStateFactorized(i))];
  TmpVector /= TmpVector.Norm();
  return TmpVector;
}

// convert a given state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

LongRationalVector BosonOnSphereHaldaneHugeBasisShort::ConvertToNbodyBasis(LongRationalVector& state, BosonOnSphereHaldaneHugeBasisShort& nbodyBasis)
{
  LongRationalVector TmpVector (nbodyBasis.GetLargeHilbertSpaceDimension(), true);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      TmpVector[nbodyBasis.FermionHugeBasis->FindStateIndexFactorized(this->FermionHugeBasis->GetStateFactorized(i))] = state[i];
    }
  return TmpVector;
}

// convert a given state from the usual n-body basis to the Haldane basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

LongRationalVector BosonOnSphereHaldaneHugeBasisShort::ConvertFromNbodyBasis(LongRationalVector& state, BosonOnSphereHaldaneHugeBasisShort& nbodyBasis)
{
  LongRationalVector TmpVector (this->LargeHilbertSpaceDimension, true);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    TmpVector[i] = state[nbodyBasis.FermionHugeBasis->FindStateIndexFactorized(this->FermionHugeBasis->GetStateFactorized(i))];
  TmpVector /= TmpVector.Norm();
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

int BosonOnSphereHaldaneHugeBasisShort::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  unsigned long TmpState = this->FermionHugeBasis->GetStateFactorized((long) index);
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  this->FermionToBoson(TmpState, TmpLzMax, 
		       this->TemporaryState, this->TemporaryStateLzMax);
  if ((n1 > this->TemporaryStateLzMax) || (n2 > this->TemporaryStateLzMax) || (this->TemporaryState[n1] == 0) || (this->TemporaryState[n2] == 0) || ((n1 == n2) && (this->TemporaryState[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  for (int i = this->TemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = 0ul;
  coefficient = this->TemporaryState[n2];
  --this->TemporaryState[n2];
  coefficient *= this->TemporaryState[n1];
  --this->TemporaryState[n1];
  ++this->TemporaryState[m2];
  coefficient *= this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->FermionHugeBasis->FindStateIndexFactorized(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax));
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

long BosonOnSphereHaldaneHugeBasisShort::AdAdAA (long index, int m1, int m2, int n1, int n2, double& coefficient)
{
  unsigned long TmpState = this->FermionHugeBasis->GetStateFactorized(index);
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  this->FermionToBoson(TmpState, TmpLzMax, 
		       this->TemporaryState, this->TemporaryStateLzMax);
  if ((n1 > this->TemporaryStateLzMax) || (n2 > this->TemporaryStateLzMax) || (this->TemporaryState[n1] == 0) || (this->TemporaryState[n2] == 0) || ((n1 == n2) && (this->TemporaryState[n1] == 1)))
    {
      coefficient = 0.0;
      return this->LargeHilbertSpaceDimension;
    }
  for (int i = this->TemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = 0ul;
  coefficient = this->TemporaryState[n2];
  --this->TemporaryState[n2];
  coefficient *= this->TemporaryState[n1];
  --this->TemporaryState[n1];
  ++this->TemporaryState[m2];
  coefficient *= this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->FermionHugeBasis->FindStateIndexFactorized(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax));
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereHaldaneHugeBasisShort::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  return 0;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereHaldaneHugeBasisShort::ProdA (int index, int* n, int nbrIndices)
{
  return 0;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereHaldaneHugeBasisShort::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  return 0;
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereHaldaneHugeBasisShort::AdA (int index, int m)
{
  unsigned long TmpState = this->FermionHugeBasis->GetStateFactorized((long) index);
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  this->FermionToBoson(TmpState, TmpLzMax, 
		       this->TemporaryState, this->TemporaryStateLzMax);
  if (this->TemporaryStateLzMax < m)  
    return 0.0;
  return (double) (this->TemporaryState[m]);  
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereHaldaneHugeBasisShort::AdA (long index, int m)
{
  unsigned long TmpState = this->FermionHugeBasis->GetStateFactorized(index);
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  this->FermionToBoson(TmpState, TmpLzMax, 
		       this->TemporaryState, this->TemporaryStateLzMax);
  if (this->TemporaryStateLzMax < m)  
    return 0.0;
  return (double) (this->TemporaryState[m]);  
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
  
int BosonOnSphereHaldaneHugeBasisShort::AdA (int index, int m, int n, double& coefficient)
{
  unsigned long TmpState = this->FermionHugeBasis->GetStateFactorized((long) index);
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  this->FermionToBoson(TmpState, TmpLzMax, 
		       this->TemporaryState, this->TemporaryStateLzMax);
  if ((this->TemporaryStateLzMax < n)  || (this->TemporaryState[n] == 0))
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState[n];
  --this->TemporaryState[n];
  if ((this->TemporaryStateLzMax == n) && (this->TemporaryState[n] == 0))
    {
      while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
	--this->TemporaryStateLzMax;
    }
  if (this->TemporaryStateLzMax < m) 
    {
      for (int i = this->TemporaryStateLzMax + 1; i <= m; ++i)
	this->TemporaryState[i] = 0;
      this->TemporaryStateLzMax = m;
    }
  ++this->TemporaryState[m];
  coefficient *= (double) this->TemporaryState[m];
  coefficient = sqrt(coefficient);  
  return this->FermionHugeBasis->FindStateIndexFactorized(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax));
  return 0;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereHaldaneHugeBasisShort::PrintStateMonomial (ostream& Str, long state)
{
  unsigned long TmpState = this->FermionHugeBasis->StateDescription[state];
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  this->ConvertToMonomial(TmpState, TmpLzMax, this->TemporaryMonomial);
  Str << "[";
  if (this->TemporaryMonomial[0] != 0)
    Str << this->TemporaryMonomial[0];
  for (int i = 1; (i < this->NbrBosons) && (this->TemporaryMonomial[i] > 0); ++i)
    Str << "," << this->TemporaryMonomial[i];
  Str << "]";
  return Str;
}

// print a given State using the monomial notation, with one column per particle (using space as a seperator)
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereHaldaneHugeBasisShort::PrintColumnFormattedStateMonomial (ostream& Str, long state)
{
  unsigned long TmpState = this->FermionHugeBasis->GetStateFactorized(state);
  int i = this->LzMax;
  while (((TmpState >> i) & 0x1ul) == 0x0ul)
    --i;
  Str << i;
  --i;
  for (; i >=0; --i)
    if (((TmpState >> i) & 0x1ul) != 0x0ul)
      Str << " " << i;
  return Str;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  BosonOnSphereHaldaneHugeBasisShort::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState)
{  
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrBosonSector == 0))
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
      if ((lzSector == this->TotalLz) && (nbrBosonSector == this->NbrBosons))
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->LargeHilbertSpaceDimension);
	  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	    for (long j = i; j < this->LargeHilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * groundState[j]);
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  if ((ShiftedLzComplementarySector < (NbrBosonsComplementarySector * subsytemSize)) || (ShiftedLzComplementarySector > (NbrBosonsComplementarySector * (this->LzMax))))
    {
      RealSymmetricMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }

  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0.0;
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  unsigned long  TmpState2 = 0x0;
	  for (int i = 0; i < nbrBosonSector; ++i)
	    TmpState2 |= 0x1ul << i;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | TmpState2;
	      long TmpPos = this->FermionHugeBasis->FindStateIndexFactorized(TmpState);
	      if (TmpPos != this->LargeHilbertSpaceDimension)
		TmpValue += groundState[TmpPos] * groundState[TmpPos];	
	    }

	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}      
    }
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0;
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	      long TmpPos = this->FermionHugeBasis->FindStateIndexFactorized(TmpState);
	      if (TmpPos != this->LargeHilbertSpaceDimension)
		TmpValue += groundState[TmpPos] * groundState[TmpPos];	
	    }	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int TmpComplementarySubsystemLzMax = this->LzMax - subsytemSize;
  long MinIndex = 0l;
  long MaxIndex = this->LargeHilbertSpaceDimension - 1l;
  if (nbrBosonSector == 1)
    {
      double TmpValue = 0.0;
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | (0x1ul << ShiftedLzSector);
	  long TmpPos = this->FermionHugeBasis->FindStateIndexFactorized(TmpState);
	  if (TmpPos != this->LargeHilbertSpaceDimension)
	    TmpValue += groundState[TmpPos] * groundState[TmpPos];	
	}
      RealSymmetricMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  if (NbrBosonsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
      BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      MinIndex = this->LargeHilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      double TmpValue;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpValue = groundState[MinIndex + i];
	  for (int j = i; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    TmpDensityMatrix.SetMatrixElement(i, j, TmpValue * groundState[MinIndex + j]);
	}
      return TmpDensityMatrix;
    }


  int TmpNbrBosons;
  int TmpTotalLz;
  int TmpIndex;
  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  long* TmpStatePosition = new long [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	  long TmpPos = this->FermionHugeBasis->FindStateIndexFactorized(TmpState);
	  if (TmpPos != this->LargeHilbertSpaceDimension)
	    {
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// create the Jack polynomial decomposition corresponding to the root partition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// partialSave = save partial results in a given vector file
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& BosonOnSphereHaldaneHugeBasisShort::GenerateJackPolynomial(RealVector& jack, double alpha, long minIndex, long maxIndex, char* partialSave)
{
  jack[0l] = 1.0;
  double InvAlpha =  2.0 / alpha;

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->FermionHugeBasis->StateDescription[0l];
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((MaxRoot >> TmpLzMax) & 0x1ul) == 0ul)
    --TmpLzMax;
  this->ConvertToMonomial(MaxRoot, TmpLzMax, this->TemporaryMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - 1.0 - InvAlpha * ((double) j));
  int ReducedNbrBosons = this->NbrBosons - 1;

  if (minIndex <= 0)
    minIndex = 1;
  if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension - 1l;
  for (long i = minIndex; i <= maxIndex; ++i)
    {
      double Rho = 0.0;
      unsigned long CurrentPartition = this->FermionHugeBasis->StateDescription[i];
      while (((CurrentPartition >> TmpLzMax) & 0x1ul) == 0ul)
	--TmpLzMax;
      this->ConvertToMonomial(CurrentPartition, TmpLzMax, this->TemporaryMonomial);
      for (int j = 0; j < this->NbrBosons; ++j)
	Rho += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - 1.0 - InvAlpha * ((double) j));
      double Coefficient = 0.0;
      for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	  {
	    double Diff = (double) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
	    unsigned int Max = this->TemporaryMonomial[j2];
	    unsigned long TmpState = 0x0ul;
	    int Tmpj1 = j1;
	    int Tmpj2 = j2;
	    for (int l = 0; l < this->NbrBosons; ++l)
	      this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
	    for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
	      {
		++this->TemporaryMonomial2[Tmpj1];
		--this->TemporaryMonomial2[Tmpj2];
		Diff += 2.0;
		while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] > this->TemporaryMonomial2[Tmpj1 - 1]))
		  {
		    unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
		    this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
		    this->TemporaryMonomial2[Tmpj1] = Tmp;
		    --Tmpj1;
		  }
                while ((Tmpj2 < ReducedNbrBosons) && (this->TemporaryMonomial2[Tmpj2] < this->TemporaryMonomial2[Tmpj2 + 1]))
                  {
                    unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
                    this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
                    this->TemporaryMonomial2[Tmpj2] = Tmp;
                    ++Tmpj2;
                  }
		TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
		if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		  {
		    long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState, this->TemporaryMonomial2[0] + ReducedNbrBosons);
		    if (TmpIndex < this->LargeHilbertSpaceDimension)
		      Coefficient += Diff * jack[TmpIndex];
		  }
	      }
	  }
      jack[i] = Coefficient * InvAlpha / (RhoRoot - Rho);
      if ((i & 0xffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	  if ((partialSave != 0) && ((i & 0xfffffffl) == 0l))
	    jack.WriteVector(partialSave);
	}
    }
  cout << endl;

  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alphaNumerator = numerator of the Jack polynomial alpha coefficient
// alphaDenominator = numerator of the Jack polynomial alpha coefficient
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// partialSave = save partial results in a given vector file
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

LongRationalVector& BosonOnSphereHaldaneHugeBasisShort::GenerateJackPolynomial(LongRationalVector& jack, long alphaNumerator, long alphaDenominator, long minIndex, long maxIndex, char* partialSave)
{
  jack[0l] = 1l;
  LongRational InvAlpha (2l * alphaDenominator, alphaNumerator);

  int ReducedNbrBosons = this->NbrBosons - 1;
  long* ConnectedIndices = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients  = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedIndices2 = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients2  = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];

  LongRational RhoRoot = 0l;
  LongRational Rho = 0l;
  unsigned long MaxRoot = this->FermionHugeBasis->StateDescription[0l];
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((MaxRoot >> TmpLzMax) & 0x1ul) == 0ul)
    --TmpLzMax;
  this->ConvertToMonomial(MaxRoot, TmpLzMax, this->TemporaryMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - 1l - InvAlpha * ((long) j));

  LongRational Coefficient = 0l;
  LongRational Coefficient2 = 0l;

  if (minIndex <= 0)
    minIndex = 1;
  if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension - 1l;
  for (long i = minIndex; i <= maxIndex; ++i)
    {
      Rho = 0l;
      unsigned long CurrentPartition = this->FermionHugeBasis->StateDescription[i];
      TmpLzMax = this->FermionHugeBasis->LzMax;
      while (((MaxRoot >> TmpLzMax) & 0x1ul) == 0ul)
	--TmpLzMax;
      this->ConvertToMonomial(CurrentPartition, TmpLzMax, this->TemporaryMonomial);
      for (int j = 0; j < this->NbrBosons; ++j)
	Rho += this->TemporaryMonomial[j] * ((this->TemporaryMonomial[j] - 1l) - InvAlpha * ((long) j));
      if (Rho == RhoRoot)
	{
	  cout << "warning : singular value detected at position " << i << ", skipping the rest of the calculation" << endl;
	  return jack;
	}
      else
	{
	  Coefficient = 0l;
	  int Pos = 0;
	  for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	    for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	      {
		long Diff = (long) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
		unsigned int Max = this->TemporaryMonomial[j2];
		unsigned long TmpState = 0x0ul;
		int Tmpj1 = j1;
		int Tmpj2 = j2;
		for (int l = 0; l < this->NbrBosons; ++l)
		  this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
		for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		  {
		    ++this->TemporaryMonomial2[Tmpj1];
		    --this->TemporaryMonomial2[Tmpj2];
		    Diff += 2l;
		    while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] > this->TemporaryMonomial2[Tmpj1 - 1]))
		      {
			unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
			this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
			this->TemporaryMonomial2[Tmpj1] = Tmp;
			--Tmpj1;
		      }
		    while ((Tmpj2 < ReducedNbrBosons) && (this->TemporaryMonomial2[Tmpj2] < this->TemporaryMonomial2[Tmpj2 + 1]))
		      {
			unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
			this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
			this->TemporaryMonomial2[Tmpj2] = Tmp;
			++Tmpj2;
		      }
		    TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
		    if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		      {
			long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState, this->TemporaryMonomial2[0] + ReducedNbrBosons);
			if (TmpIndex < this->LargeHilbertSpaceDimension)
			  {
			    ConnectedIndices[Pos] = TmpIndex;
			    ConnectedCoefficients[Pos] = Diff;
			    ++Pos;
			  }
		      }
		  }
	      }
	      int NbrConnected = 1l;
	      if (Pos > 1)
		{
		  SortArrayDownOrdering<long>(ConnectedIndices, ConnectedCoefficients, Pos);
		  int TmpIndex = 1;
		  while (TmpIndex < Pos)
		    {
		      while ((TmpIndex < Pos) && (ConnectedIndices[TmpIndex] == ConnectedIndices[TmpIndex - 1]))
			++TmpIndex;
		      if (TmpIndex < Pos)
			++NbrConnected;
		      ++TmpIndex;
		    }
		  ConnectedIndices2[0] = ConnectedIndices[0];
		  ConnectedCoefficients2[0] = ConnectedCoefficients[0];
		  TmpIndex = 1;
		  NbrConnected = 1;
		  while (TmpIndex < Pos)
		    {
		      while ((TmpIndex < Pos) && (ConnectedIndices[TmpIndex] == ConnectedIndices[TmpIndex - 1]))
			{
			  ConnectedCoefficients2[NbrConnected - 1] += ConnectedCoefficients[TmpIndex];
			  ++TmpIndex;
			}
		      if (TmpIndex < Pos)
			{
			  ConnectedIndices2[NbrConnected] = ConnectedIndices[TmpIndex];
			  ConnectedCoefficients2[NbrConnected] = ConnectedCoefficients[TmpIndex];	   
			  ++NbrConnected;
			}
		      ++TmpIndex;
		    }
		}
	      else
		{
		  ConnectedIndices2[0] = ConnectedIndices[0];
		  ConnectedCoefficients2[0] = ConnectedCoefficients[0];
		}
 	      Coefficient = ConnectedCoefficients2[0];	  
	      Coefficient *= jack[ConnectedIndices2[0]];
	      for (int j = 1; j < NbrConnected; ++j)
		{
		  Coefficient2 = ConnectedCoefficients2[j];
		  Coefficient2 *= jack[ConnectedIndices2[j]];
		  Coefficient += Coefficient2;
		}
	      Coefficient *= InvAlpha;
	      Rho -= RhoRoot;
	      Rho.Neg();
	      Coefficient /= Rho;
	      jack[i] = Coefficient;
	    }
      if ((i & 0xfffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	  if ((partialSave != 0) && ((i & 0xfffffffl) == 0l))
	    jack.WriteVector(partialSave);
	}
    }
  cout << endl;
  delete[] ConnectedIndices;
  delete[] ConnectedCoefficients;
  delete[] ConnectedIndices2;
  delete[] ConnectedCoefficients2;
  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// partialSave = save partial results in a given vector file
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& BosonOnSphereHaldaneHugeBasisShort::GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha, long minIndex, long maxIndex, char* partialSave)
{
  jack[0l] = 1.0;
  double InvAlpha =  2.0 / alpha;

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->FermionHugeBasis->StateDescription[0l];
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((MaxRoot >> TmpLzMax) & 0x1ul) == 0ul)
    --TmpLzMax;
  this->ConvertToMonomial(MaxRoot, TmpLzMax, this->TemporaryMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - 1.0 - InvAlpha * ((double) j));
  int ReducedNbrBosons = this->NbrBosons - 1;
  if (minIndex <= 0)
    minIndex = 1;
  if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension - 1l;
  for (long i = minIndex; i <= maxIndex; ++i)
    if (jack[i] == 0.0)
      {
	double Rho = 0.0;
	unsigned long CurrentPartition = this->FermionHugeBasis->StateDescription[i];
	while (((CurrentPartition >> TmpLzMax) & 0x1ul) == 0ul)
	  --TmpLzMax;
	this->ConvertToMonomial(CurrentPartition, TmpLzMax, this->TemporaryMonomial);
// 	if (i != this->FermionHugeBasis->FindStateIndexMemory(CurrentPartition, this->TemporaryMonomial[0] + ReducedNbrBosons))
// 	  cout << i << " " << this->FermionHugeBasis->FindStateIndexMemory(CurrentPartition, this->TemporaryMonomial[0] + ReducedNbrBosons) << endl;
	for (int j = 0; j < this->NbrBosons; ++j)
	  Rho += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - 1.0 - InvAlpha * ((double) j));
	double Coefficient = 0.0;
	for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	  for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	    {
	      double Diff = (double) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
	      unsigned int Max = this->TemporaryMonomial[j2];
	      unsigned long TmpState = 0x0ul;
	      int Tmpj1 = j1;
	      int Tmpj2 = j2;
	      for (int l = 0; l < this->NbrBosons; ++l)
		this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
	      for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		{
		  ++this->TemporaryMonomial2[Tmpj1];
		  --this->TemporaryMonomial2[Tmpj2];
		  Diff += 2.0;
		  while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] > this->TemporaryMonomial2[Tmpj1 - 1]))
		    {
		      unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
		      this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
		      this->TemporaryMonomial2[Tmpj1] = Tmp;
		      --Tmpj1;
		    }
		  while ((Tmpj2 < ReducedNbrBosons) && (this->TemporaryMonomial2[Tmpj2] < this->TemporaryMonomial2[Tmpj2 + 1]))
		    {
		      unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
		      this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
		      this->TemporaryMonomial2[Tmpj2] = Tmp;
		      ++Tmpj2;
		    }
		  TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
		  if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		    {
		      long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState, this->TemporaryMonomial2[0] + ReducedNbrBosons);
		      if (TmpIndex < this->LargeHilbertSpaceDimension)
			Coefficient += Diff * jack[TmpIndex];
		    }
		}
	    }
	
	long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(this->FermionHugeBasis->GetSymmetricState(CurrentPartition), (this->LzMax - this->TemporaryMonomial2[ReducedNbrBosons]) + ReducedNbrBosons);
	Coefficient *= InvAlpha;
	Coefficient /= (RhoRoot - Rho);
	if (i < TmpIndex)
	  jack[TmpIndex] = Coefficient;
	jack[i] = Coefficient;
	if ((i & 0xffl) == 0l)
	  {
	    cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	    cout.flush();
	    if ((partialSave != 0) && ((i & 0xfffffffl) == 0l))
	      jack.WriteVector(partialSave);
	  }
      }
  cout << endl;

  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition and using sparse storage
//
// alpha = value of the Jack polynomial alpha coefficient
// architecture = architecture to use for precalculation
// partialSave = save partial results in a given vector file
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// memory = amount of memory (in bytes) allowed for temporary vector storage (0 if the whole vector has to be stored in memory)
// memoryBlock = amount of memory (in bytes) allowed for precomputing state indices
// resumeFlag = true if the calculation has to be resumed from a previous one (assuming partialSave contains already computed components)

void  BosonOnSphereHaldaneHugeBasisShort::GenerateJackPolynomialSparse(double alpha, AbstractArchitecture* architecture, char* partialSave, long minIndex, long maxIndex, long memory, long memoryBlock, bool resumeFlag)
{
  if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension - 1l;
  double TmpComponent = 1.0;
  long FileShift = 4l;
  if (this->HilbertSpaceDimension <= 0)
    FileShift = 12l;
  if ((minIndex <= 0) && (resumeFlag == false))
    {
      ofstream File;
      File.open(partialSave, ios::binary | ios::out);
      WriteLittleEndian(File, this->HilbertSpaceDimension);  
      if (this->HilbertSpaceDimension <= 0)
	{
	  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);  
	}
      WriteLittleEndian(File, TmpComponent);  
      File.close();
    }
  TmpComponent = 0.0;
  
  memory >>= 3;
  if ((memory > this->LargeHilbertSpaceDimension) || (memory <= 0))
    {
      cout << "vector does not require temporary disk storage" << endl;
      memory = this->LargeHilbertSpaceDimension;
    }
  double* TmpVectorBuffer = new double [memory];
  long BufferGlobalIndex = 0l;


  double InvAlpha =  2.0 / alpha;

  double RhoRoot = 0.0;
  unsigned long MaxRoot = 0x0ul;
  MaxRoot = this->FermionHugeBasis->GetStateFactorized(0l);
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((MaxRoot >> TmpLzMax) & 0x1ul) == 0ul)
    --TmpLzMax;
  this->ConvertToMonomial(MaxRoot, TmpLzMax, this->TemporaryMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - 1.0 - InvAlpha * ((double) j));
  if (minIndex <= 0)
    minIndex = 1;
  TmpVectorBuffer[0l] = 1.0;
  int MaxArraySize = ((this->NbrBosons * (this->NbrBosons - 1)) / 2) * (this->LzMax + 1);
  long NbrBlocks =  memoryBlock / (((2l* sizeof (long)) + (sizeof(double))) * MaxArraySize);
  if (NbrBlocks == 0)
    {
      NbrBlocks = 100000;
    }
  if (NbrBlocks > this->LargeHilbertSpaceDimension)
    {
      NbrBlocks = (this->LargeHilbertSpaceDimension >> 1) + 1;
    }
  cout << "number of precalculation blocks = " << NbrBlocks << endl;
  long DisplayStep = (this->LargeHilbertSpaceDimension / (1000 * NbrBlocks)) * NbrBlocks;
  long** TmpIndexArray = new long* [NbrBlocks];
  double** TmpComponentArray = new double* [NbrBlocks];
  unsigned long** TmpStateArray = new unsigned long* [NbrBlocks]; 
  double* TmpRhoArray = new double [NbrBlocks];
  int* TmpNbrComputedComponentArray = new int [NbrBlocks];
  for (int j = 0; j < NbrBlocks; ++j)
    {
      TmpIndexArray[j] = new long [MaxArraySize];
      TmpComponentArray[j] = new double [MaxArraySize];
      TmpStateArray[j] = new unsigned long [MaxArraySize];
    }
  FQHESphereJackGeneratorOperation Operation(this, InvAlpha, MaxRoot, TmpIndexArray, TmpStateArray, TmpComponentArray, TmpRhoArray, TmpNbrComputedComponentArray, false, false);

  if (resumeFlag == true)
    {
      ifstream File;
      File.open(partialSave, ios::binary | ios::in);
      File.seekg (0, ios::end);
      long TmpResumePos = File.tellg();
      File.close();
      TmpResumePos -= FileShift;
      TmpResumePos /= sizeof(double); 	      
      long TmpResumeMinPos = TmpResumePos - NbrBlocks;
      long LimNbrBlocks = NbrBlocks;
      if (TmpResumeMinPos < 0l)
	{
	  TmpResumeMinPos = 0l;
	  LimNbrBlocks = TmpResumePos - TmpResumeMinPos + 1;
	}
      long TmpMaxIndex = TmpResumeMinPos + NbrBlocks - 1l;
      if (TmpMaxIndex > maxIndex)
	{
	  LimNbrBlocks = NbrBlocks - (TmpMaxIndex - maxIndex);
	  TmpMaxIndex = maxIndex;
	}
      if (LimNbrBlocks > 0)
	{
	  cout << "consistency check, " << TmpResumePos << " components have already been computed, checking the last " << LimNbrBlocks << " ones" << endl;      
	  Operation.SetIndicesRange(TmpResumeMinPos, LimNbrBlocks);
	  Operation.ApplyOperation(architecture);
	  ifstream OutputFile;
	  OutputFile.open(partialSave, ios::binary | ios::in);
	  double RefCoefficient = 0.0;

	  for (long k = 0l; k < LimNbrBlocks; ++k)
	    {
	      OutputFile.seekg (FileShift + (TmpResumeMinPos * sizeof(double)), ios::beg);
	      ReadLittleEndian(OutputFile, RefCoefficient);
	      double Coefficient = 0.0;
	      if (TmpNbrComputedComponentArray[k] >= 0)
		{
		  for (int j = 0; j < TmpNbrComputedComponentArray[k]; ++j)
		    {
		      long TmpIndex = TmpIndexArray[k][j];
		      if (TmpIndex < this->LargeHilbertSpaceDimension)
			{		  
			  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
			  ReadLittleEndian (OutputFile, TmpComponent);
			  Coefficient += TmpComponentArray[k][j] * TmpComponent;
			}	      	    
		    }		  
		  Coefficient *= InvAlpha;
		  Coefficient /= (RhoRoot - TmpRhoArray[k]);
		}
	      else
		{
		  long TmpIndex = TmpIndexArray[k][0];
		  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
		  ReadLittleEndian (OutputFile, Coefficient);
		}
	      if (Coefficient != RefCoefficient)
		{
		  cout << "error, invalid Jack : component " << TmpResumeMinPos << " is " << RefCoefficient << ", should be " << Coefficient << endl;
		  OutputFile.close();
		  return;
		}
	      ++TmpResumeMinPos;
	    }
	  TmpResumeMinPos = TmpResumePos - memory;
	  if (TmpResumeMinPos < 0l)
	    TmpResumeMinPos = 0l;
	  BufferGlobalIndex = TmpResumeMinPos;
	  OutputFile.seekg ((TmpResumeMinPos * sizeof(double)) + FileShift, ios::beg);
	  double TmpComponent;
	  for (; TmpResumeMinPos < TmpResumePos; ++TmpResumeMinPos)
	    {	      
	      ReadLittleEndian (OutputFile, TmpComponent);
	      TmpVectorBuffer[TmpResumeMinPos % memory] = TmpComponent;	      
	    }
	  OutputFile.close();
	}
      cout << "consistency check done, resuming calculation now" << endl;
      minIndex = TmpResumePos;
    }

  fstream OutputFile;
  OutputFile.open(partialSave, ios::in | ios::binary | ios::out);

  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);

  for (long i = minIndex; i <= maxIndex;)
    {
      long TmpMaxIndex = i + NbrBlocks - 1l;
      long LimNbrBlocks = NbrBlocks;
      if (TmpMaxIndex > maxIndex)
	{
	  LimNbrBlocks = NbrBlocks - (TmpMaxIndex - maxIndex);
	  TmpMaxIndex = maxIndex;
	}
      Operation.SetIndicesRange(i, LimNbrBlocks);
      Operation.ApplyOperation(architecture);
      for (long k = 0l; k < LimNbrBlocks; ++k)
	{
	  if (TmpNbrComputedComponentArray[k] >= 0)
	    {
	      double Coefficient = 0.0;
	      for (int j = 0; j < TmpNbrComputedComponentArray[k]; ++j)
		{
		  long TmpIndex = TmpIndexArray[k][j];
		  if (TmpIndex < this->LargeHilbertSpaceDimension)
		    {		  
		      if (TmpIndex < BufferGlobalIndex)
			{
			  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
			  ReadLittleEndian (OutputFile, TmpComponent);
			}
		      else
			{
			  TmpComponent = TmpVectorBuffer[TmpIndex % memory];
			}
		      Coefficient += TmpComponentArray[k][j] * TmpComponent;
		    }	      	    
		}
 
	      Coefficient *= InvAlpha;
	      Coefficient /= (RhoRoot - TmpRhoArray[k]);
	      OutputFile.seekg (0, ios::end);
	      WriteLittleEndian(OutputFile, Coefficient);
	      if (i >= memory)
		++BufferGlobalIndex;
	      TmpVectorBuffer[i % memory] = Coefficient;
	      ++i;
	    }
	  else
	    {
	      long TmpIndex = TmpIndexArray[k][0];
	      if (TmpIndex < BufferGlobalIndex)
		{
		  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
		  ReadLittleEndian (OutputFile, TmpComponent);
		}
	      else
		{
		  TmpComponent = TmpVectorBuffer[TmpIndex % memory];
		}
	      OutputFile.seekg (0, ios::end);
	      WriteLittleEndian(OutputFile, TmpComponent); 	  
	      if (i >= memory)
		++BufferGlobalIndex;
	      TmpVectorBuffer[i % memory] = TmpComponent;
	      ++i;
	    }
	}
      if ((i & DisplayStep) == 0l)
      	{
     	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
      	  cout.flush();
      	}
    }
  OutputFile.close();
  delete[] TmpStateArray;
  delete[] TmpIndexArray;
  delete[] TmpComponentArray;
  cout << endl;
}

// create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry and using sparse storage
//
// alpha = value of the Jack polynomial alpha coefficient
// architecture = architecture to use for precalculation
// partialSave = save partial results in a given vector file
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// memory = amount of memory (in bytes) allowed for temporary vector storage (0 if the whole vector has to be stored in memory)
// memoryBlock = amount of memory (in bytes) allowed for precomputing state indices
// resumeFlag = true if the calculation has to be resumed from a previous one (assuming partialSave contains already computed components)

void BosonOnSphereHaldaneHugeBasisShort::GenerateSymmetrizedJackPolynomialSparse(double alpha, AbstractArchitecture* architecture, char* partialSave, long minIndex, long maxIndex, long memory, long memoryBlock, bool resumeFlag)
{
  if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension - 1l;
  double TmpComponent = 1.0;
  long FileShift = 4l;
  if (this->HilbertSpaceDimension <= 0)
    FileShift = 12l;
  if ((minIndex <= 0) && (resumeFlag == false))
    {
      ofstream File;
      File.open(partialSave, ios::binary | ios::out);
      WriteLittleEndian(File, this->HilbertSpaceDimension);  
      if (this->HilbertSpaceDimension <= 0)
	{
	  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);  
	}
      WriteLittleEndian(File, TmpComponent);  
      File.close();
    }
  TmpComponent = 0.0;
  
  memory >>= 3;
  if ((memory > this->LargeHilbertSpaceDimension) || (memory <= 0))
    {
      cout << "vector does not require temporary disk storage" << endl;
      memory = this->LargeHilbertSpaceDimension;
    }
  double* TmpVectorBuffer = new double [memory];
  long BufferGlobalIndex = 0l;


  double InvAlpha =  2.0 / alpha;

  double RhoRoot = 0.0;
  unsigned long MaxRoot = 0x0ul;
  MaxRoot = this->FermionHugeBasis->GetStateFactorized(0l);
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((MaxRoot >> TmpLzMax) & 0x1ul) == 0ul)
    --TmpLzMax;
  this->ConvertToMonomial(MaxRoot, TmpLzMax, this->TemporaryMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - 1.0 - InvAlpha * ((double) j));
  if (minIndex <= 0)
    minIndex = 1;
  TmpVectorBuffer[0l] = 1.0;
  int MaxArraySize = ((this->NbrBosons * (this->NbrBosons - 1)) / 2) * (this->LzMax + 1);
  long NbrBlocks =  memoryBlock / (((2l* sizeof (long)) + (sizeof(double))) * MaxArraySize);
  if (NbrBlocks == 0)
    {
      NbrBlocks = 100000;
    }
  if (NbrBlocks > this->LargeHilbertSpaceDimension)
    {
      NbrBlocks = (this->LargeHilbertSpaceDimension >> 1) + 1;
    }
  cout << "number of precalculation blocks = " << NbrBlocks << endl;
  long DisplayStep = (this->LargeHilbertSpaceDimension / (1000 * NbrBlocks)) * NbrBlocks;
  long** TmpIndexArray = new long* [NbrBlocks];
  double** TmpComponentArray = new double* [NbrBlocks];
  unsigned long** TmpStateArray = new unsigned long* [NbrBlocks]; 
  double* TmpRhoArray = new double [NbrBlocks];
  int* TmpNbrComputedComponentArray = new int [NbrBlocks];
  for (int j = 0; j < NbrBlocks; ++j)
    {
      TmpIndexArray[j] = new long [MaxArraySize];
      TmpComponentArray[j] = new double [MaxArraySize];
      TmpStateArray[j] = new unsigned long [MaxArraySize];
    }
  FQHESphereJackGeneratorOperation Operation(this, InvAlpha, MaxRoot, TmpIndexArray, TmpStateArray, TmpComponentArray, TmpRhoArray, TmpNbrComputedComponentArray, false, true);

  if (resumeFlag == true)
    {
      ifstream File;
      File.open(partialSave, ios::binary | ios::in);
      File.seekg (0, ios::end);
      long TmpResumePos = File.tellg();
      File.close();
      TmpResumePos -= FileShift;
      TmpResumePos /= sizeof(double); 	      
      long TmpResumeMinPos = TmpResumePos - NbrBlocks;
      long LimNbrBlocks = NbrBlocks;
      if (TmpResumeMinPos < 0l)
	{
	  TmpResumeMinPos = 0l;
	  LimNbrBlocks = TmpResumePos - TmpResumeMinPos + 1;
	}
      long TmpMaxIndex = TmpResumeMinPos + NbrBlocks - 1l;
      if (TmpMaxIndex > maxIndex)
	{
	  LimNbrBlocks = NbrBlocks - (TmpMaxIndex - maxIndex);
	  TmpMaxIndex = maxIndex;
	}
      if (LimNbrBlocks > 0)
	{
	  cout << "consistency check, " << TmpResumePos << " components have already been computed, checking the last " << LimNbrBlocks << " ones" << endl;      
	  Operation.SetIndicesRange(TmpResumeMinPos, LimNbrBlocks);
	  Operation.ApplyOperation(architecture);
	  ifstream OutputFile;
	  OutputFile.open(partialSave, ios::binary | ios::in);
	  double RefCoefficient = 0.0;

	  for (long k = 0l; k < LimNbrBlocks; ++k)
	    {
	      OutputFile.seekg (FileShift + (TmpResumeMinPos * sizeof(double)), ios::beg);
	      ReadLittleEndian(OutputFile, RefCoefficient);
	      double Coefficient = 0.0;
	      if (TmpNbrComputedComponentArray[k] >= 0)
		{
		  for (int j = 0; j < TmpNbrComputedComponentArray[k]; ++j)
		    {
		      long TmpIndex = TmpIndexArray[k][j];
		      if (TmpIndex < this->LargeHilbertSpaceDimension)
			{		  
			  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
			  ReadLittleEndian (OutputFile, TmpComponent);
			  Coefficient += TmpComponentArray[k][j] * TmpComponent;
			}	      	    
		    }		  
		  Coefficient *= InvAlpha;
		  Coefficient /= (RhoRoot - TmpRhoArray[k]);
		}
	      else
		{
		  long TmpIndex = TmpIndexArray[k][0];
		  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
		  ReadLittleEndian (OutputFile, Coefficient);
		}
	      if (Coefficient != RefCoefficient)
		{
		  cout << "error, invalid Jack : component " << TmpResumeMinPos << " is " << RefCoefficient << ", should be " << Coefficient << endl;
		  OutputFile.close();
		  return;
		}
	      ++TmpResumeMinPos;
	    }
	  TmpResumeMinPos = TmpResumePos - memory;
	  if (TmpResumeMinPos < 0l)
	    TmpResumeMinPos = 0l;
	  BufferGlobalIndex = TmpResumeMinPos;
	  OutputFile.seekg ((TmpResumeMinPos * sizeof(double)) + FileShift, ios::beg);
	  double TmpComponent;
	  for (; TmpResumeMinPos < TmpResumePos; ++TmpResumeMinPos)
	    {	      
	      ReadLittleEndian (OutputFile, TmpComponent);
	      TmpVectorBuffer[TmpResumeMinPos % memory] = TmpComponent;	      
	    }
	  OutputFile.close();
	}
      cout << "consistency check done, resuming calculation now" << endl;
      minIndex = TmpResumePos;
    }

  fstream OutputFile;
  OutputFile.open(partialSave, ios::in | ios::binary | ios::out);

  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);

  for (long i = minIndex; i <= maxIndex;)
    {
      long TmpMaxIndex = i + NbrBlocks - 1l;
      long LimNbrBlocks = NbrBlocks;
      if (TmpMaxIndex > maxIndex)
	{
	  LimNbrBlocks = NbrBlocks - (TmpMaxIndex - maxIndex);
	  TmpMaxIndex = maxIndex;
	}
      Operation.SetIndicesRange(i, LimNbrBlocks);
      Operation.ApplyOperation(architecture);
      for (long k = 0l; k < LimNbrBlocks; ++k)
	{
	  if (TmpNbrComputedComponentArray[k] >= 0)
	    {
	      double Coefficient = 0.0;
	      for (int j = 0; j < TmpNbrComputedComponentArray[k]; ++j)
		{
		  long TmpIndex = TmpIndexArray[k][j];
		  if (TmpIndex < this->LargeHilbertSpaceDimension)
		    {		  
		      if (TmpIndex < BufferGlobalIndex)
			{
			  long TmpIndex2 = this->FermionHugeBasis->FindStateIndexFactorized(this->FermionHugeBasis->GetSymmetricState(TmpStateArray[k][j]));
			  if ((TmpIndex2 >= BufferGlobalIndex) && (TmpIndex2 < i))
			    {
			      TmpComponent = TmpVectorBuffer[TmpIndex2 % memory];
			    }
			  else
			    {
			      OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
			      ReadLittleEndian (OutputFile, TmpComponent);
			    }
			}
		      else
			{
			  TmpComponent = TmpVectorBuffer[TmpIndex % memory];
			}
		      Coefficient += TmpComponentArray[k][j] * TmpComponent;
		    }	      	    
		}
 
	      Coefficient *= InvAlpha;
	      Coefficient /= (RhoRoot - TmpRhoArray[k]);
	      OutputFile.seekg (0, ios::end);
	      WriteLittleEndian(OutputFile, Coefficient);
	      if (i >= memory)
		++BufferGlobalIndex;
	      TmpVectorBuffer[i % memory] = Coefficient;
	      ++i;
	    }
	  else
	    {
	      long TmpIndex = TmpIndexArray[k][0];
	      if (TmpIndex < BufferGlobalIndex)
		{
		  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
		  ReadLittleEndian (OutputFile, TmpComponent);
		}
	      else
		{
		  TmpComponent = TmpVectorBuffer[TmpIndex % memory];
		}
	      OutputFile.seekg (0, ios::end);
	      WriteLittleEndian(OutputFile, TmpComponent); 	  
	      if (i >= memory)
		++BufferGlobalIndex;
	      TmpVectorBuffer[i % memory] = TmpComponent;
	      ++i;
	    }
	}
      if ((i & DisplayStep) == 0l)
      	{
	  //	  gettimeofday (&(TotalEndingTime), 0);
     	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
      	  cout.flush();
	  //	  double Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	  //	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  //	  cout << "done in " << Dt << " s" << endl;
	  //	  gettimeofday (&(TotalStartingTime), 0);
      	}
    }
  OutputFile.close();
  delete[] TmpStateArray;
  delete[] TmpIndexArray;
  delete[] TmpComponentArray;
  cout << endl;
}

// create the Jack polynomial decomposition corresponding to the root partition and using sparse storage
//
// alpha = value of the Jack polynomial alpha coefficient
// architecture = architecture to use for precalculation
// partialSave = save partial results in a given vector file
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// memory = amount of memory (in bytes) allowed for temporary vector storage (0 if the whole vector has to be stored in memory)
// memoryBlock = amount of memory (in bytes) allowed for precomputing state indices
// resumeFlag = true if the calculation has to be resumed from a previous one (assuming partialSave contains already computed components)

void  BosonOnSphereHaldaneHugeBasisShort::GenerateJackPolynomialSparse(long alphaNumerator, long alphaDenominator, AbstractArchitecture* architecture, char* partialSave, long minIndex, long maxIndex, long memory, long memoryBlock, bool resumeFlag)
{
//   if ((maxIndex <= 0) || (maxIndex >= this->LargeHilbertSpaceDimension))
//     maxIndex = this->LargeHilbertSpaceDimension - 1l;
//   double TmpComponent = 1.0;
//   long FileShift = 4l;
//   if (this->HilbertSpaceDimension <= 0)
//     FileShift = 12l;
//   if ((minIndex <= 0) && (resumeFlag == false))
//     {
//       ofstream File;
//       File.open(partialSave, ios::binary | ios::out);
//       WriteLittleEndian(File, this->HilbertSpaceDimension);  
//       if (this->HilbertSpaceDimension <= 0)
// 	{
// 	  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);  
// 	}
//       WriteLittleEndian(File, TmpComponent);  
//       File.close();
//     }
//   TmpComponent = 0.0;
  
//   memory >>= 3;
//   if ((memory > this->LargeHilbertSpaceDimension) || (memory <= 0))
//     {
//       cout << "vector does not require temporary disk storage" << endl;
//       memory = this->LargeHilbertSpaceDimension;
//     }
//   LongRational* TmpVectorBuffer = new LongRational [memory];
//   long BufferGlobalIndex = 0l;


//   LongRational InvAlpha =  2.0 / alpha;

//   LongRational RhoRoot = 0.0;
//   unsigned long MaxRoot = 0x0ul;
//   MaxRoot = this->FermionHugeBasis->GetStateFactorized(0l);
//   int TmpLzMax = this->FermionHugeBasis->LzMax;
//   while (((MaxRoot >> TmpLzMax) & 0x1ul) == 0ul)
//     --TmpLzMax;
//   this->ConvertToMonomial(MaxRoot, TmpLzMax, this->TemporaryMonomial);
//   for (int j = 0; j < this->NbrBosons; ++j)
//     RhoRoot += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - 1l - InvAlpha * ((long) j));
//   if (minIndex <= 0)
//     minIndex = 1;
//   TmpVectorBuffer[0l] = 1l;
//   int MaxArraySize = ((this->NbrBosons * (this->NbrBosons - 1)) / 2) * (this->LzMax + 1);
//   long NbrBlocks =  memoryBlock / (((2l* sizeof (long)) + (sizeof(double))) * MaxArraySize);
//   if (NbrBlocks == 0)
//     {
//       NbrBlocks = 100000;
//     }
//   if (NbrBlocks > this->LargeHilbertSpaceDimension)
//     {
//       NbrBlocks = (this->LargeHilbertSpaceDimension >> 1) + 1;
//     }
//   cout << "number of precalculation blocks = " << NbrBlocks << endl;
//   long DisplayStep = (this->LargeHilbertSpaceDimension / (1000 * NbrBlocks)) * NbrBlocks;
//   long** TmpIndexArray = new long* [NbrBlocks];
//   double** TmpComponentArray = new LongRational* [NbrBlocks];
//   unsigned long** TmpStateArray = new unsigned long* [NbrBlocks]; 
//   double* TmpRhoArray = new LongRational [NbrBlocks];
//   int* TmpNbrComputedComponentArray = new int [NbrBlocks];
//   for (int j = 0; j < NbrBlocks; ++j)
//     {
//       TmpIndexArray[j] = new long [MaxArraySize];
//       TmpComponentArray[j] = new LongRational [MaxArraySize];
//       TmpStateArray[j] = new unsigned long [MaxArraySize];
//     }
//   FQHESphereJackGeneratorOperation Operation(this, InvAlpha, MaxRoot, TmpIndexArray, TmpStateArray, TmpComponentArray, TmpRhoArray, TmpNbrComputedComponentArray, false, false);

//   if (resumeFlag == true)
//     {
// //       ifstream File;
// //       File.open(partialSave, ios::binary | ios::in);
// //       File.seekg (0, ios::end);
// //       long TmpResumePos = File.tellg();
// //       File.close();
// //       TmpResumePos -= FileShift;
// //       TmpResumePos /= sizeof(double); 	      
// //       long TmpResumeMinPos = TmpResumePos - NbrBlocks;
// //       long LimNbrBlocks = NbrBlocks;
// //       if (TmpResumeMinPos < 0l)
// // 	{
// // 	  TmpResumeMinPos = 0l;
// // 	  LimNbrBlocks = TmpResumePos - TmpResumeMinPos + 1;
// // 	}
// //       long TmpMaxIndex = TmpResumeMinPos + NbrBlocks - 1l;
// //       if (TmpMaxIndex > maxIndex)
// // 	{
// // 	  LimNbrBlocks = NbrBlocks - (TmpMaxIndex - maxIndex);
// // 	  TmpMaxIndex = maxIndex;
// // 	}
// //       if (LimNbrBlocks > 0)
// // 	{
// // 	  cout << "consistency check, " << TmpResumePos << " components have already been computed, checking the last " << LimNbrBlocks << " ones" << endl;      
// // 	  Operation.SetIndicesRange(TmpResumeMinPos, LimNbrBlocks);
// // 	  Operation.ApplyOperation(architecture);
// // 	  ifstream OutputFile;
// // 	  OutputFile.open(partialSave, ios::binary | ios::in);
// // 	  double RefCoefficient = 0.0;

// // 	  for (long k = 0l; k < LimNbrBlocks; ++k)
// // 	    {
// // 	      OutputFile.seekg (FileShift + (TmpResumeMinPos * sizeof(double)), ios::beg);
// // 	      ReadLittleEndian(OutputFile, RefCoefficient);
// // 	      double Coefficient = 0.0;
// // 	      if (TmpNbrComputedComponentArray[k] >= 0)
// // 		{
// // 		  for (int j = 0; j < TmpNbrComputedComponentArray[k]; ++j)
// // 		    {
// // 		      long TmpIndex = TmpIndexArray[k][j];
// // 		      if (TmpIndex < this->LargeHilbertSpaceDimension)
// // 			{		  
// // 			  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
// // 			  ReadLittleEndian (OutputFile, TmpComponent);
// // 			  Coefficient += TmpComponentArray[k][j] * TmpComponent;
// // 			}	      	    
// // 		    }		  
// // 		  Coefficient *= InvAlpha;
// // 		  Coefficient /= (RhoRoot - TmpRhoArray[k]);
// // 		}
// // 	      else
// // 		{
// // 		  long TmpIndex = TmpIndexArray[k][0];
// // 		  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
// // 		  ReadLittleEndian (OutputFile, Coefficient);
// // 		}
// // 	      if (Coefficient != RefCoefficient)
// // 		{
// // 		  cout << "error, invalid Jack : component " << TmpResumeMinPos << " is " << RefCoefficient << ", should be " << Coefficient << endl;
// // 		  OutputFile.close();
// // 		  return;
// // 		}
// // 	      ++TmpResumeMinPos;
// // 	    }
// // 	  TmpResumeMinPos = TmpResumePos - memory;
// // 	  if (TmpResumeMinPos < 0l)
// // 	    TmpResumeMinPos = 0l;
// // 	  BufferGlobalIndex = TmpResumeMinPos;
// // 	  OutputFile.seekg ((TmpResumeMinPos * sizeof(double)) + FileShift, ios::beg);
// // 	  double TmpComponent;
// // 	  for (; TmpResumeMinPos < TmpResumePos; ++TmpResumeMinPos)
// // 	    {	      
// // 	      ReadLittleEndian (OutputFile, TmpComponent);
// // 	      TmpVectorBuffer[TmpResumeMinPos % memory] = TmpComponent;	      
// // 	    }
// // 	  OutputFile.close();
// // 	}
// //       cout << "consistency check done, resuming calculation now" << endl;
// //       minIndex = TmpResumePos;
//     }

//   fstream OutputFile;
//   OutputFile.open(partialSave, ios::in | ios::binary | ios::out);

//   timeval TotalStartingTime;
//   timeval TotalEndingTime;
//   gettimeofday (&(TotalStartingTime), 0);

//   for (long i = minIndex; i <= maxIndex;)
//     {
//       long TmpMaxIndex = i + NbrBlocks - 1l;
//       long LimNbrBlocks = NbrBlocks;
//       if (TmpMaxIndex > maxIndex)
// 	{
// 	  LimNbrBlocks = NbrBlocks - (TmpMaxIndex - maxIndex);
// 	  TmpMaxIndex = maxIndex;
// 	}
//       Operation.SetIndicesRange(i, LimNbrBlocks);
//       Operation.ApplyOperation(architecture);
//       for (long k = 0l; k < LimNbrBlocks; ++k)
// 	{
// 	  if (TmpNbrComputedComponentArray[k] >= 0)
// 	    {
// 	      double Coefficient = 0.0;
// 	      for (int j = 0; j < TmpNbrComputedComponentArray[k]; ++j)
// 		{
// 		  long TmpIndex = TmpIndexArray[k][j];
// // 		  if (TmpIndex < this->LargeHilbertSpaceDimension)
// // 		    {		  
// // 		      if (TmpIndex < BufferGlobalIndex)
// // 			{
// // 			  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
// // 			  ReadLittleEndian (OutputFile, TmpComponent);
// // 			}
// // 		      else
// // 			{
// 			  TmpComponent = TmpVectorBuffer[TmpIndex % memory];
// 			  //			}
// 		      Coefficient += TmpComponentArray[k][j] * TmpComponent;
// 		    }	      	    
// 		}
 
// 	      Coefficient *= InvAlpha;
// 	      Coefficient /= (RhoRoot - TmpRhoArray[k]);
// 	      OutputFile.seekg (0, ios::end);
// 	      WriteLittleEndian(OutputFile, Coefficient);
// 	      if (i >= memory)
// 		++BufferGlobalIndex;
// 	      TmpVectorBuffer[i % memory] = Coefficient;
// 	      ++i;
// 	    }
// 	  else
// 	    {
// 	      long TmpIndex = TmpIndexArray[k][0];
// 	      if (TmpIndex < BufferGlobalIndex)
// 		{
// 		  OutputFile.seekg ((TmpIndex * sizeof(double)) + FileShift, ios::beg);
// 		  ReadLittleEndian (OutputFile, TmpComponent);
// 		}
// 	      else
// 		{
// 		  TmpComponent = TmpVectorBuffer[TmpIndex % memory];
// 		}
// 	      OutputFile.seekg (0, ios::end);
// 	      WriteLittleEndian(OutputFile, TmpComponent); 	  
// 	      if (i >= memory)
// 		++BufferGlobalIndex;
// 	      TmpVectorBuffer[i % memory] = TmpComponent;
// 	      ++i;
// 	    }
// 	}
//       if ((i & DisplayStep) == 0l)
//       	{
//      	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
//       	  cout.flush();
//       	}
//     }
//   OutputFile.close();
//   delete[] TmpStateArray;
//   delete[] TmpIndexArray;
//   delete[] TmpComponentArray;
//   cout << endl;
}

// core part of the Jack generator using the Lz<->-Lz symmetry and the factorized algorithm
//
// invAlpha = inverse of the Jack polynomial alpha coefficient
// maxRoot = root partition (in fermionic binary representation)
// partialSave = save partial results in a given vector file
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// indexArray = array where state indices are stored
// stateArray = array use to store computed state description
// componentArray = array where computed component numerical factors are stored
// nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
// rhoArray = rho factor associated to each state

void BosonOnSphereHaldaneHugeBasisShort::GenerateSymmetrizedJackPolynomialFactorizedCore(double invAlpha, unsigned long maxRoot, long minIndex, long maxIndex, unsigned long** stateArray, double** componentArray, long** indexArray, int* nbrComputedComponents, double* rhoArray)
{
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  int ReducedNbrBosons = this->NbrBosons - 1;
  for (long i = minIndex; i <= maxIndex; ++i)
    {
      unsigned long* TmpStateArray = stateArray[i - minIndex];
      double* TmpComponentArray = componentArray[i - minIndex];
      long* TmpIndexArray = indexArray[i - minIndex];
      double Rho = 0.0;
      unsigned long CurrentPartition = 0x0ul;
      CurrentPartition = this->FermionHugeBasis->GetStateFactorized(i);
      unsigned long TmpSymState = this->FermionHugeBasis->GetSymmetricState(CurrentPartition);
      if (TmpSymState <= CurrentPartition)
	{
	  while (((CurrentPartition >> TmpLzMax) & 0x1ul) == 0ul)
	    --TmpLzMax;
	  this->ConvertToMonomial(CurrentPartition, TmpLzMax, this->TemporaryMonomial);
	  for (int j = 0; j < this->NbrBosons; ++j)
	    Rho += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - 1.0 - invAlpha * ((double) j));
	  rhoArray[i - minIndex] = Rho;
	  int NbrComputedComponents = 0;
	  for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	    for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	      {
		double Diff = (double) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
		unsigned int Max = this->TemporaryMonomial[j2];
		unsigned long TmpState = 0x0ul;
		int Tmpj1 = j1;
		int Tmpj2 = j2;
		for (int l = 0; l < this->NbrBosons; ++l)
		  this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
		for (unsigned int k = 1; (k <= Max) && (TmpState < maxRoot); ++k)
		  {
		    ++this->TemporaryMonomial2[Tmpj1];
		    --this->TemporaryMonomial2[Tmpj2];
		    Diff += 2.0;
		    while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] > this->TemporaryMonomial2[Tmpj1 - 1]))
		      {
			unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
			this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
			this->TemporaryMonomial2[Tmpj1] = Tmp;
			--Tmpj1;
		      }
		    while ((Tmpj2 < ReducedNbrBosons) && (this->TemporaryMonomial2[Tmpj2] < this->TemporaryMonomial2[Tmpj2 + 1]))
		      {
			unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
			this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
			this->TemporaryMonomial2[Tmpj2] = Tmp;
			++Tmpj2;
		      }
		    TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
		    if ((TmpState <= maxRoot) && (TmpState > CurrentPartition))
		      {
			TmpComponentArray[NbrComputedComponents] = Diff;
			TmpStateArray[NbrComputedComponents] = TmpState;
			++NbrComputedComponents;
		      }
		  }
	      }
	  nbrComputedComponents[i - minIndex] = NbrComputedComponents;
	  for (int j = 0; j < NbrComputedComponents; ++j)
	    {
	      TmpIndexArray[j] = this->FermionHugeBasis->FindStateIndexFactorized(TmpStateArray[j]);
	    }
	}
      else
	{
	  nbrComputedComponents[i - minIndex] = -1;
	   TmpIndexArray[0] = this->FermionHugeBasis->FindStateIndexFactorized(TmpSymState);
	}
    }
}

// core part of the Jack generator using the factorized algorithm
//
// invAlpha = inverse of the Jack polynomial alpha coefficient
// maxRoot = root partition (in fermionic binary representation)
// partialSave = save partial results in a given vector file
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// indexArray = array where state indices are stored
// stateArray = array use to store computed state description
// componentArray = array where computed component numerical factors are stored
// nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
// rhoArray = rho factor associated to each state

void BosonOnSphereHaldaneHugeBasisShort::GenerateJackPolynomialFactorizedCore(double invAlpha, unsigned long maxRoot, long minIndex, long maxIndex, unsigned long** stateArray, double** componentArray, long** indexArray, int* nbrComputedComponents, double* rhoArray)
{
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  int ReducedNbrBosons = this->NbrBosons - 1;
  for (long i = minIndex; i <= maxIndex; ++i)
    {
      unsigned long* TmpStateArray = stateArray[i - minIndex];
      double* TmpComponentArray = componentArray[i - minIndex];
      long* TmpIndexArray = indexArray[i - minIndex];
      double Rho = 0.0;
      unsigned long CurrentPartition = 0x0ul;
      CurrentPartition = this->FermionHugeBasis->GetStateFactorized(i);
      while (((CurrentPartition >> TmpLzMax) & 0x1ul) == 0ul)
	--TmpLzMax;
      this->ConvertToMonomial(CurrentPartition, TmpLzMax, this->TemporaryMonomial);
      for (int j = 0; j < this->NbrBosons; ++j)
	Rho += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - 1.0 - invAlpha * ((double) j));
      rhoArray[i - minIndex] = Rho;
      int NbrComputedComponents = 0;
      for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	  {
	    double Diff = (double) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
	    unsigned int Max = this->TemporaryMonomial[j2];
	    unsigned long TmpState = 0x0ul;
	    int Tmpj1 = j1;
	    int Tmpj2 = j2;
	    for (int l = 0; l < this->NbrBosons; ++l)
	      this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
	    for (unsigned int k = 1; (k <= Max) && (TmpState < maxRoot); ++k)
	      {
		++this->TemporaryMonomial2[Tmpj1];
		--this->TemporaryMonomial2[Tmpj2];
		Diff += 2.0;
		while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] > this->TemporaryMonomial2[Tmpj1 - 1]))
		  {
		    unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
		    this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
		    this->TemporaryMonomial2[Tmpj1] = Tmp;
		    --Tmpj1;
		  }
		while ((Tmpj2 < ReducedNbrBosons) && (this->TemporaryMonomial2[Tmpj2] < this->TemporaryMonomial2[Tmpj2 + 1]))
		  {
		    unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
		    this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
		    this->TemporaryMonomial2[Tmpj2] = Tmp;
		    ++Tmpj2;
		  }
		TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
		if ((TmpState <= maxRoot) && (TmpState > CurrentPartition))
		  {
		    TmpComponentArray[NbrComputedComponents] = Diff;
		    TmpStateArray[NbrComputedComponents] = TmpState;
		    ++NbrComputedComponents;
		  }
	      }
	  }
      nbrComputedComponents[i - minIndex] = NbrComputedComponents;
      for (int j = 0; j < NbrComputedComponents; ++j)
	{
	  TmpIndexArray[j] = this->FermionHugeBasis->FindStateIndexFactorized(TmpStateArray[j]);
	}
    }
}

// core part of the Jack generator using the factorized algorithm
//
// invAlpha = inverse of the Jack polynomial alpha coefficient
// maxRoot = root partition (in fermionic binary representation)
// partialSave = save partial results in a given vector file
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// indexArray = array where state indices are stored
// stateArray = array use to store computed state description
// componentArray = array where computed component numerical factors are stored
// nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
// rhoArray = rho factor associated to each state

void BosonOnSphereHaldaneHugeBasisShort::GenerateJackPolynomialFactorizedCore(LongRational invAlpha, unsigned long maxRoot, long minIndex, long maxIndex, unsigned long** stateArray, LongRational** componentArray, long** indexArray, int* nbrComputedComponents, LongRational* rhoArray)
{
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  int ReducedNbrBosons = this->NbrBosons - 1;
  for (long i = minIndex; i <= maxIndex; ++i)
    {
      unsigned long* TmpStateArray = stateArray[i - minIndex];
      LongRational* TmpComponentArray = componentArray[i - minIndex];
      long* TmpIndexArray = indexArray[i - minIndex];
      LongRational Rho = 0l;
      unsigned long CurrentPartition = 0x0ul;
      CurrentPartition = this->FermionHugeBasis->GetStateFactorized(i);
      while (((CurrentPartition >> TmpLzMax) & 0x1ul) == 0ul)
	--TmpLzMax;
      this->ConvertToMonomial(CurrentPartition, TmpLzMax, this->TemporaryMonomial);
      for (int j = 0; j < this->NbrBosons; ++j)
	Rho += this->TemporaryMonomial[j] * (this->TemporaryMonomial[j] - 1l - invAlpha * ((long) j));
      rhoArray[i - minIndex] = Rho;
      int NbrComputedComponents = 0;
      for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	  {
	    long Diff = (long) (this->TemporaryMonomial[j1] - this->TemporaryMonomial[j2]);
	    unsigned int Max = this->TemporaryMonomial[j2];
	    unsigned long TmpState = 0x0ul;
	    int Tmpj1 = j1;
	    int Tmpj2 = j2;
	    for (int l = 0; l < this->NbrBosons; ++l)
	      this->TemporaryMonomial2[l] = this->TemporaryMonomial[l];	    
	    for (unsigned int k = 1; (k <= Max) && (TmpState < maxRoot); ++k)
	      {
		++this->TemporaryMonomial2[Tmpj1];
		--this->TemporaryMonomial2[Tmpj2];
		Diff += 2.0;
		while ((Tmpj1 > 0) && (this->TemporaryMonomial2[Tmpj1] > this->TemporaryMonomial2[Tmpj1 - 1]))
		  {
		    unsigned long Tmp = this->TemporaryMonomial2[Tmpj1 - 1];
		    this->TemporaryMonomial2[Tmpj1 - 1] = this->TemporaryMonomial2[Tmpj1];
		    this->TemporaryMonomial2[Tmpj1] = Tmp;
		    --Tmpj1;
		  }
		while ((Tmpj2 < ReducedNbrBosons) && (this->TemporaryMonomial2[Tmpj2] < this->TemporaryMonomial2[Tmpj2 + 1]))
		  {
		    unsigned long Tmp = this->TemporaryMonomial2[Tmpj2 + 1];
		    this->TemporaryMonomial2[Tmpj2 + 1] = this->TemporaryMonomial2[Tmpj2];
		    this->TemporaryMonomial2[Tmpj2] = Tmp;
		    ++Tmpj2;
		  }
		TmpState = this->ConvertFromMonomial(this->TemporaryMonomial2);
		if ((TmpState <= maxRoot) && (TmpState > CurrentPartition))
		  {
		    TmpComponentArray[NbrComputedComponents] = Diff;
		    TmpStateArray[NbrComputedComponents] = TmpState;
		    ++NbrComputedComponents;
		  }
	      }
	  }
      nbrComputedComponents[i - minIndex] = NbrComputedComponents;
      for (int j = 0; j < NbrComputedComponents; ++j)
	{
	  TmpIndexArray[j] = this->FermionHugeBasis->FindStateIndexFactorized(TmpStateArray[j]);
	}
    }
}

	


// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& BosonOnSphereHaldaneHugeBasisShort::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  double Factor = 1.0 / state[reference];
  state[reference] = 1.0;
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  unsigned long TmpState = 0x0ul;
  if (this->FermionHugeBasis->CheckDiskStorage() == false)
    {
      TmpState = this->FermionHugeBasis->StateDescription[reference];
    }
  else
    {
      TmpState = this->FermionHugeBasis->GetStateFactorized(reference);
    }
   while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  this->ConvertToMonomial(TmpState, TmpLzMax, this->TemporaryMonomial2);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(TmpState, TmpLzMax, 
		       this->TemporaryState, this->TemporaryStateLzMax);
  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialDivide(this->TemporaryState[k]);
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    if (i != reference)
      {
	if (this->FermionHugeBasis->CheckDiskStorage() == false)
	  {
	    TmpState = this->FermionHugeBasis->StateDescription[i];
	  }
	else
	  {
	    TmpState = this->FermionHugeBasis->GetStateFactorized(i);
	  }
	while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	  --TmpLzMax;
	this->ConvertToMonomial(TmpState, TmpLzMax, this->TemporaryMonomial);
	int Index1 = 0;
	int Index2 = 0;
	double Coefficient = Factor;
	while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	  {
	    while ((Index1 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] > this->TemporaryMonomial[Index2]))
	      {
		Coefficient *= InvSqrtCoefficients[this->TemporaryMonomial2[Index1]];
		++Index1;
	      }
	    while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] == this->TemporaryMonomial[Index2]))
	      {
		++Index1;
		++Index2;
	      }
	    while ((Index2 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] < this->TemporaryMonomial[Index2]))
	      {
		Coefficient *= SqrtCoefficients[this->TemporaryMonomial[Index2]];
		++Index2;
	      }	  
	  }
	while (Index1 < this->NbrBosons)
	  {
	    Coefficient *= InvSqrtCoefficients[this->TemporaryMonomial2[Index1]];
	    ++Index1;
	  }
	while (Index2 < this->NbrBosons)
	  {
	    Coefficient *= SqrtCoefficients[this->TemporaryMonomial2[Index2]];
	    ++Index2;
	  }
	if (symmetryFactor == true)
	  {
	    Factorial = ReferenceFactorial;
	    this->FermionToBoson(TmpState, TmpLzMax, 
				 this->TemporaryState, this->TemporaryStateLzMax);
	    for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	      if (this->TemporaryState[k] > 1)
		Factorial.FactorialMultiply(this->TemporaryState[k]);
	    Coefficient *= sqrt(Factorial.GetNumericalValue());
	  }
	state[i] *= Coefficient;
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
      }
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& BosonOnSphereHaldaneHugeBasisShort::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  double Factor = 1.0;
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  unsigned long TmpState = 0x0ul;
  if (this->FermionHugeBasis->CheckDiskStorage() == false)
    {
      TmpState = this->FermionHugeBasis->StateDescription[reference];
    }
  else
    {
      TmpState = this->FermionHugeBasis->GetStateFactorized(reference);
    }
   while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
   this->ConvertToMonomial(TmpState, TmpLzMax, this->TemporaryMonomial2);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
    }
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(TmpState, TmpLzMax, 
		       this->TemporaryState, this->TemporaryStateLzMax);
  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialMultiply(this->TemporaryState[k]);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    if (i != reference)
      {
	if (this->FermionHugeBasis->CheckDiskStorage() == false)
	  {
	    TmpState = this->FermionHugeBasis->StateDescription[i];
	  }
	else
	  {
	    TmpState = this->FermionHugeBasis->GetStateFactorized(i);
	  }
	while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	  --TmpLzMax;
	this->ConvertToMonomial(TmpState, TmpLzMax, this->TemporaryMonomial);
	int Index1 = 0;
	int Index2 = 0;
	double Coefficient = Factor;
	if (symmetryFactor == true)
	  {
	    while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	      {
		while ((Index1 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] > this->TemporaryMonomial[Index2]))
		  {
		    Coefficient *= InvSqrtCoefficients[this->TemporaryMonomial2[Index1]];
		    ++Index1;
		  }
		while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] == this->TemporaryMonomial[Index2]))
		  {
		    ++Index1;
		    ++Index2;
		  }
		while ((Index2 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] < this->TemporaryMonomial[Index2]))
		  {
		    Coefficient *= SqrtCoefficients[this->TemporaryMonomial[Index2]];
		    ++Index2;
		  }	  
	      }
	    while (Index1 < this->NbrBosons)
	      {
		Coefficient *= InvSqrtCoefficients[this->TemporaryMonomial2[Index1]];
		++Index1;
	      }
	    while (Index2 < this->NbrBosons)
	      {
		Coefficient *= SqrtCoefficients[this->TemporaryMonomial2[Index2]];
		++Index2;
	      }
	    Factorial = ReferenceFactorial;
	    this->FermionToBoson(TmpState, TmpLzMax, 
				 this->TemporaryState, this->TemporaryStateLzMax);
	    for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	      if (this->TemporaryState[k] > 1)
		Factorial.FactorialDivide(this->TemporaryState[k]);
	    Coefficient *= sqrt(Factorial.GetNumericalValue());
	  }
	else
	  {
	    Factorial = ReferenceFactorial;
	    this->FermionToBoson(TmpState, TmpLzMax, 
				 this->TemporaryState, this->TemporaryStateLzMax);
	    for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	      if (this->TemporaryState[k] > 1)
		Factorial.FactorialDivide(this->TemporaryState[k]);
	    Coefficient *= sqrt(Factorial.GetNumericalValue());
	  }
	state[i] *= Coefficient;
	if ((i & 0x3fffl) == 0l)
	  {
	    cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	    cout.flush();
	  }
      }
  state /= state.Norm();
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
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

RealVector& BosonOnSphereHaldaneHugeBasisShort::FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
							    ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
							    bool symmetrizedFlag, double coefficient)
{
  BosonOnSphereShort* LeftSpace = (BosonOnSphereShort*) leftSpace;
  BosonOnSphereShort* RightSpace = (BosonOnSphereShort*) rightSpace;
  int StateShift = RightSpace->FermionBasis->LzMax + padding + 2;
  long Count = 0l;
  for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long TmpState1 = LeftSpace->FermionBasis->StateDescription[i] << StateShift;
      int TmpLzMax = this->FermionHugeBasis->LzMax;
      while ((TmpState1 >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      double Coefficient = coefficient * leftVector[i];
      if (symmetrizedFlag == false)
	{
	  if (this->FermionHugeBasis->CheckDiskStorage() == true)
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
		  TmpState2 |= TmpState1;
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  long TmpIndex = this->FermionHugeBasis->FindStateIndexFactorized(TmpState2);
		  double& TmpCoef = outputVector[TmpIndex];
		  if (TmpCoef == 0.0)
		    ++Count;
		  TmpCoef = Coefficient2;
		}
	    }
	  else
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
		  TmpState2 |= TmpState1;
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState2, TmpLzMax);
		  double& TmpCoef = outputVector[TmpIndex];
		  if (TmpCoef == 0.0)
		    ++Count;
		  TmpCoef = Coefficient2;
		}
	    }
	}
      else
	{
	  if (this->FermionHugeBasis->CheckDiskStorage() == true)
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
		  TmpState2 |= TmpState1;
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  long TmpIndex = this->FermionHugeBasis->FindStateIndexFactorized(TmpState2);
		  double& TmpCoef = outputVector[TmpIndex];
		  if (TmpCoef == 0.0)
		    ++Count;
		  TmpCoef = Coefficient2;
		  unsigned long TmpState3 = this->FermionHugeBasis->GetSymmetricState(TmpState2);
		  if (TmpState3 != TmpState2)
		    {
		      TmpIndex = this->FermionHugeBasis->FindStateIndexFactorized(TmpState3);
		      double& TmpCoef2 = outputVector[TmpIndex];
		      if (TmpCoef2 == 0.0)
			++Count;
		      TmpCoef2 = Coefficient2;
		    }
		}
	    }
	  else
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
		  TmpState2 |= TmpState1;
		  double Coefficient2 = Coefficient;
		  Coefficient2 *= rightVector[j];	  
		  long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState2, TmpLzMax);
		  double& TmpCoef = outputVector[TmpIndex];
		  if (TmpCoef == 0.0)
		    ++Count;
		  TmpCoef = Coefficient2;
		  unsigned long TmpState3 = this->FermionHugeBasis->GetSymmetricState(TmpState2);
		  if (TmpState3 != TmpState2)
		    {
		      int TmpLzMax2 = this->FermionHugeBasis->LzMax;
		      while ((TmpState3 >> TmpLzMax2) == 0x0ul)
			--TmpLzMax2;
		      TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState3, TmpLzMax2);
		      double& TmpCoef2 = outputVector[TmpIndex];
		      if (TmpCoef2 == 0.0)
			++Count;
		      TmpCoef2 = Coefficient2;
		    }
		}
	    }
	}
    }
  cout << "nbr of newly added components : " << Count << endl;
  return outputVector;
}

// fuse multiple states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// nbrInputVectors = number of input vectors
// inputVectors = input vectors whose Hilbert space will be fuse from  left to right
// paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
// inputSpaces = point to the Hilbert space that will be fuse to the left
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

RealVector& BosonOnSphereHaldaneHugeBasisShort::FuseMultipleStates (RealVector& outputVector, int nbrInputVectors, RealVector* inputVectors, int* paddings, 
								    ParticleOnSphere** inputSpaces, bool symmetrizedFlag, double coefficient)
{
  BosonOnSphereShort** InputSpaces = (BosonOnSphereShort**) inputSpaces;
  for (long i = 0; i <  InputSpaces[0]->LargeHilbertSpaceDimension; ++i)
    {
      this->CoreFuseMultipleStates(outputVector, nbrInputVectors, inputVectors, paddings, InputSpaces, 1, 
				   ((BosonOnSphereShort*) InputSpaces[0])->FermionBasis->StateDescription[i], 
				   ((BosonOnSphereShort*) InputSpaces[0])->FermionBasis->LzMax + paddings[0] + 2,
				   coefficient  * inputVectors[0][i], symmetrizedFlag);	    
      
    }
  return outputVector;
}

// core part of multiple state fuse 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// nbrInputVectors = number of input vectors
// inputVectors = input vectors whose Hilbert space will be fuse from  left to right
// paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
// inputSpaces = point to the Hilbert space that will be fuse to the left
// currentPosition = index of the current space to fuse
// currentState = current fermionic state obtained by fusing previous states
// currentCoefficient = current multiplicative coefficient
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry

void BosonOnSphereHaldaneHugeBasisShort::CoreFuseMultipleStates (RealVector& outputVector, int nbrInputVectors, RealVector* inputVectors, int* paddings, 
								 BosonOnSphereShort** inputSpaces, int currentPosition, unsigned long currentState, int currentPadding, 
								 double currentCoefficient, bool symmetrizedFlag)
{
  if (currentPosition < (nbrInputVectors - 1))
    {
      for (long i = 0; i <  inputSpaces[currentPosition]->LargeHilbertSpaceDimension; ++i)
	{
	  this->CoreFuseMultipleStates(outputVector, nbrInputVectors, inputVectors, paddings, inputSpaces, currentPosition + 1, 
				       currentState | (inputSpaces[currentPosition]->FermionBasis->StateDescription[i] << currentPadding), 
				       currentPadding + inputSpaces[currentPosition]->FermionBasis->LzMax + paddings[currentPosition] + 2,
				       currentCoefficient * inputVectors[currentPosition][i], symmetrizedFlag);	    
	  
	}
    }
  else
    {
      FermionOnSphere* RightSpace = inputSpaces[currentPosition]->FermionBasis;
      RealVector& RightVector = inputVectors[currentPosition];
      if (symmetrizedFlag == false)
	{
	  if (this->FermionHugeBasis->CheckDiskStorage() == true)
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState = RightSpace->StateDescription[j] << currentPadding;
		  TmpState |= currentState;
		  long TmpIndex = this->FermionHugeBasis->FindStateIndexFactorized(TmpState);
		  outputVector[TmpIndex] = currentCoefficient * RightVector[j];
		}
	    }
	  else
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState = RightSpace->StateDescription[j] << currentPadding;
		  TmpState |= currentState;
		  int TmpLzMax = this->FermionBasis->LzMax;
		  while ((TmpState >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState, TmpLzMax);
		  outputVector[TmpIndex] = currentCoefficient * RightVector[j];
		}
	    }
	}
      else
	{
	  if (this->FermionHugeBasis->CheckDiskStorage() == true)
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState = RightSpace->StateDescription[j] << currentPadding;
		  TmpState |= currentState;
		  double Coefficient = currentCoefficient;
		  Coefficient *= RightVector[j];	  
		  long TmpIndex = this->FermionHugeBasis->FindStateIndexFactorized(TmpState);
		  outputVector[TmpIndex] = Coefficient;
		  unsigned long TmpState2 = this->FermionBasis->GetSymmetricState(TmpState);
		  if (TmpState != TmpState2)
		    {
		      TmpIndex = this->FermionHugeBasis->FindStateIndexFactorized(TmpState2);
		      outputVector[TmpIndex] = Coefficient;      
		    }
		}
	    }
	  else
	    {
	      for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
		{
		  unsigned long TmpState = RightSpace->StateDescription[j] << currentPadding;
		  TmpState |= currentState;
		  double Coefficient = currentCoefficient;
		  Coefficient *= RightVector[j];	  
		  int TmpLzMax = this->FermionBasis->LzMax;
		  while ((TmpState >> TmpLzMax) == 0x0ul)
		    --TmpLzMax;
		  long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState, TmpLzMax);
		  outputVector[TmpIndex] = Coefficient;
		  unsigned long TmpState2 = this->FermionBasis->GetSymmetricState(TmpState);
		  if (TmpState != TmpState2)
		    {
		      TmpLzMax = this->FermionBasis->LzMax;
		      while ((TmpState2 >> TmpLzMax) == 0x0ul)
			--TmpLzMax;
		      TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState2, TmpLzMax);
		      outputVector[TmpIndex] = Coefficient;      
		    }
		}
	    }
	}
    }
}

// use product rule to produce part of the components of a system from a smaller one
//
// outputVector = reference on the vector which will contain the product rule state  (without zeroing components which do not occur in the fusion)
// inputVector = reference on the vector associated to the smaller system
// inputSpace = pointer to the Hilbert space of the smaller system
// commonPattern = array describing the shared leftmost pattern between the n-body states in both the smaller and larger system sizes
// commonPatterSize = number of elements in the commonPattern array
// addedPattern = array describing the pattern that has to be inserted to go from the smaller system to the larger one
// addedPatterSize = number of elements in the addedPattern array
// coefficient = multiplicqtive fqctor to go fron the component of the smaller system to the larger one
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// return value = reference on the product rule state

RealVector& BosonOnSphereHaldaneHugeBasisShort::ProductRules (RealVector& outputVector, RealVector& inputVector, ParticleOnSphere* inputSpace, 
							      int* commonPattern, int commonPatterSize, int* addedPattern, int addedPatterSize,
							      double coefficient, bool symmetrizedFlag)
{
  BosonOnSphereShort* InputSpace = (BosonOnSphereShort*) inputSpace;
  int NbrParticlesCommonPattern = 0;
  unsigned long InputPattern = 0x0ul;
  int TmpIndex = 0;
  for (int i = commonPatterSize - 1; i >= 0; --i)
    {
      for (int j = 0; j < commonPattern[i]; ++j)
	{
	  InputPattern |= 0x1ul << TmpIndex;
	  ++TmpIndex;
	}
      ++TmpIndex;
      NbrParticlesCommonPattern += commonPattern[i];
    }
  int NbrParticlesAddedPattern = 0;
  unsigned long OutputPattern = 0x0ul;
  TmpIndex = 0;
  for (int i = addedPatterSize - 1; i >= 0; --i)
    {
      for (int j = 0; j < addedPattern[i]; ++j)
	{
	  OutputPattern |= 0x1ul << TmpIndex;
	  ++TmpIndex;
	}
      ++TmpIndex;
      NbrParticlesAddedPattern += addedPattern[i];
    }
  unsigned long InputMask = ((0x1ul << (commonPatterSize + NbrParticlesCommonPattern)) - 1ul) << (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 1);
  unsigned long InputMask2 = ~InputMask;  
  OutputPattern |= InputPattern << (addedPatterSize + NbrParticlesAddedPattern);
  OutputPattern <<= (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 2);
  InputPattern <<= (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 2); 
  int OutputLzMax = this->FermionHugeBasis->LzMax;
  while (((OutputPattern >> OutputLzMax) & 0x1ul) == 0x0ul)
    --OutputLzMax;
  long Count = 0l;
  if (symmetrizedFlag == false)
    {
      for (long i = 0; i <  InputSpace->LargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState1 = InputSpace->FermionBasis->StateDescription[i];
	  if ((TmpState1 & InputMask) == InputPattern)
	    {
	      TmpState1 &= InputMask2;
	      TmpState1 |= OutputPattern;
	      long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState1, OutputLzMax);
	      double& TmpCoef = outputVector[TmpIndex];
	      if (TmpCoef == 0.0)
		++Count;
	      TmpCoef = coefficient * inputVector[i];	  
	    }
	}
    }
  else
    {
      for (long i = 0; i <  InputSpace->LargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState1 = InputSpace->FermionBasis->StateDescription[i];
	  if ((TmpState1 & InputMask) == InputPattern)
	    {
	      TmpState1 &= InputMask2;
	      TmpState1 |= OutputPattern;
	      long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState1, OutputLzMax);
	      double& TmpCoef = outputVector[TmpIndex];
	      if (TmpCoef == 0.0)
		++Count;
	      double TmpCoef3 = coefficient * inputVector[i];
	      TmpCoef = TmpCoef3;	  
	      unsigned long TmpState2 = this->FermionHugeBasis->GetSymmetricState(TmpState1);
	      if (TmpState2 != TmpState1)
		{
		  int TmpLzMax2 = this->FermionHugeBasis->LzMax;
		  while ((TmpState2 >> TmpLzMax2) == 0x0ul)
		    --TmpLzMax2;
		  TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState2, TmpLzMax2);
		  double& TmpCoef2 = outputVector[TmpIndex];
                  if (TmpCoef2 == 0.0)
                    ++Count;
                  TmpCoef2 = TmpCoef3;
		}
	    }
	}
    }
  cout << "nbr of newly added components : " << Count << endl;
  return outputVector;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

double BosonOnSphereHaldaneHugeBasisShort::JackSqrNormalization (RealVector& outputVector, long minIndex, long nbrComponents)
{
  double SqrNorm = 0.0;
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = ((unsigned long) this->LzMax) >> 1;
  int TmpLzMax= this->FermionHugeBasis->LzMax;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  if (this->FermionHugeBasis->NbrRootSuffix == 0)
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  unsigned long TmpState = this->FermionHugeBasis->StateDescription[i];
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  Factorial.SetToOne();
	  this->ConvertToMonomial(TmpState, TmpLzMax, this->TemporaryMonomial);
	  this->FermionToBoson(TmpState, TmpLzMax, 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k < this->NbrBosons; ++k)
	    {
	      if (HalfLzMax < this->TemporaryMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, this->TemporaryMonomial[k]);
	      else
		if (HalfLzMax > this->TemporaryMonomial[k])
		  Factorial.PartialFactorialDivide(this->TemporaryMonomial[k] + 1, HalfLzMax);
	    }	      
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialDivide(this->TemporaryState[k]);
	  SqrNorm += (outputVector[i] * outputVector[i]) * Factorial.GetNumericalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  else
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  unsigned long TmpState = this->FermionHugeBasis->GetStateFactorized(i);
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  Factorial.SetToOne();
	  this->ConvertToMonomial(TmpState, TmpLzMax, this->TemporaryMonomial);
	  this->FermionToBoson(TmpState, TmpLzMax, 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k < this->NbrBosons; ++k)
	    {
	      if (HalfLzMax < this->TemporaryMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, this->TemporaryMonomial[k]);
	      else
		if (HalfLzMax > this->TemporaryMonomial[k])
		  Factorial.PartialFactorialDivide(this->TemporaryMonomial[k] + 1, HalfLzMax);
	    }	      
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialDivide(this->TemporaryState[k]);
	  SqrNorm += (outputVector[i] * outputVector[i]) * Factorial.GetNumericalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  cout << endl;
  return SqrNorm;
}


// compute part of the Jack polynomial scalar product in a given range of indices
//
// state1 = reference on the first unnormalized Jack polynomial
// state2 = reference on the second unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

double BosonOnSphereHaldaneHugeBasisShort::JackScalarProduct (RealVector& state1, RealVector& state2, long minIndex, long nbrComponents)
{
  double SqrNorm = 0.0;
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = ((unsigned long) this->LzMax) >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  int TmpLzMax= this->FermionHugeBasis->LzMax;
  if (this->FermionHugeBasis->NbrRootSuffix == 0)
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  unsigned long TmpState = this->FermionHugeBasis->StateDescription[i];
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  Factorial.SetToOne();
	  this->ConvertToMonomial(TmpState, TmpLzMax, TmpMonomial);
	  this->FermionToBoson(TmpState, TmpLzMax,
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k < this->NbrBosons; ++k)
	    {
	      if (HalfLzMax < TmpMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomial[k]);
	      else
		if (HalfLzMax > TmpMonomial[k])
		  Factorial.PartialFactorialDivide(TmpMonomial[k] + 1, HalfLzMax);
	    }	      
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialDivide(this->TemporaryState[k]);
	  SqrNorm +=(state1[i] * state2[i]) * Factorial.GetNumericalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  else
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  unsigned long TmpState = this->FermionHugeBasis->GetStateFactorized(i);
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  Factorial.SetToOne();
	  this->ConvertToMonomial(TmpState, TmpLzMax, TmpMonomial);
	  this->FermionToBoson(TmpState, TmpLzMax,
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k < this->NbrBosons; ++k)
	    {
	      if (HalfLzMax < TmpMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomial[k]);
	      else
		if (HalfLzMax > TmpMonomial[k])
		  Factorial.PartialFactorialDivide(TmpMonomial[k] + 1, HalfLzMax);
	    }	      
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialDivide(this->TemporaryState[k]);
	  SqrNorm +=(state1[i] * state2[i]) * Factorial.GetNumericalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  cout << endl;
  delete[] TmpMonomial;
  return SqrNorm;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state1 = reference on the first unnormalized Jack polynomial
// state2 = reference on the second unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

LongRational BosonOnSphereHaldaneHugeBasisShort::JackScalarProduct (LongRationalVector& state1, LongRationalVector& state2, long minIndex, long nbrComponents)
{
  LongRational SqrNorm = 0l;
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = ((unsigned long) this->LzMax) >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  int TmpLzMax= this->FermionHugeBasis->LzMax;
  if (this->FermionHugeBasis->NbrRootSuffix == 0)
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  unsigned long TmpState = this->FermionHugeBasis->StateDescription[i];
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  Factorial.SetToOne();
	  this->ConvertToMonomial(TmpState, TmpLzMax, TmpMonomial);
	  this->FermionToBoson(TmpState, TmpLzMax, 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k < this->NbrBosons; ++k)
	    {
	      if (HalfLzMax < TmpMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomial[k]);
	      else
		if (HalfLzMax > TmpMonomial[k])
		  Factorial.PartialFactorialDivide(TmpMonomial[k] + 1, HalfLzMax);
	    }	      
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialDivide(this->TemporaryState[k]);
	  SqrNorm += (state1[i] * state2[i]) * Factorial.GetLongRationalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  else
    {
      for (long i = minIndex; i < MaxIndex; ++i)
	{
	  unsigned long TmpState = this->FermionHugeBasis->GetStateFactorized(i);
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  Factorial.SetToOne();
	  this->ConvertToMonomial(TmpState, TmpLzMax, TmpMonomial);
	  this->FermionToBoson(TmpState, TmpLzMax, 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k < this->NbrBosons; ++k)
	    {
	      if (HalfLzMax < TmpMonomial[k])
		Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomial[k]);
	      else
		if (HalfLzMax > TmpMonomial[k])
		  Factorial.PartialFactorialDivide(TmpMonomial[k] + 1, HalfLzMax);
	    }	      
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialDivide(this->TemporaryState[k]);
	  SqrNorm += (state1[i] * state2[i]) * Factorial.GetLongRationalValue();
	  if ((i & 0x3fffl) == 0l)
	    {
	      cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	      cout.flush();
	    }
	}
    }
  cout << endl;
  delete[] TmpMonomial;
  return SqrNorm;
}

