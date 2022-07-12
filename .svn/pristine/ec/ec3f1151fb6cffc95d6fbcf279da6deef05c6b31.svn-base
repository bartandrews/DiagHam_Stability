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
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/ArrayTools.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h" 

#include <math.h>


using std::cout;
using std::endl;


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
  this->FermionBasis = 0;
  this->FermionHugeBasis = new FermionOnSphereHaldaneHugeBasis(fileName, memoryHilbert);
  this->HilbertSpaceDimension = this->FermionHugeBasis->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->FermionHugeBasis->GetLargeHilbertSpaceDimension();

  this->NbrBosons = this->FermionHugeBasis->NbrFermions;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = this->FermionHugeBasis->TotalLz;
  this->LzMax = this->FermionHugeBasis->LzMax - this->NbrBosons + 1;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;

  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
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
}

// destructor
//

BosonOnSphereHaldaneHugeBasisShort::~BosonOnSphereHaldaneHugeBasisShort ()
{
  if (this->FermionHugeBasis != 0)
    delete this->FermionHugeBasis;
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

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereHaldaneHugeBasisShort::PrintStateMonomial (ostream& Str, int state)
{
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long TmpState = this->FermionHugeBasis->StateDescription[state];
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  this->ConvertToMonomial(TmpState, TmpLzMax, TmpMonomial);
  Str << "[";
  if (TmpMonomial[0] != 0)
    Str << TmpMonomial[0];
  for (int i = 1; (i < this->NbrBosons) && (TmpMonomial[i] > 0); ++i)
    Str << "," << TmpMonomial[i];
  Str << "]";
  delete[] TmpMonomial;
  return Str;
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

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->FermionHugeBasis->StateDescription[0l];
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((MaxRoot >> TmpLzMax) & 0x1ul) == 0ul)
    --TmpLzMax;
  this->ConvertToMonomial(MaxRoot, TmpLzMax, TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
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
      this->ConvertToMonomial(CurrentPartition, TmpLzMax, TmpMonomial);
      for (int j = 0; j < this->NbrBosons; ++j)
	Rho += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
      double Coefficient = 0.0;
      for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	  {
	    double Diff = (double) (TmpMonomial[j1] - TmpMonomial[j2]);
	    unsigned int Max = TmpMonomial[j2];
	    unsigned long TmpState = 0x0ul;
	    int Tmpj1 = j1;
	    int Tmpj2 = j2;
	    for (int l = 0; l < this->NbrBosons; ++l)
	      TmpMonomial2[l] = TmpMonomial[l];	    
	    for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
	      {
		++TmpMonomial2[Tmpj1];
		--TmpMonomial2[Tmpj2];
		Diff += 2.0;
		while ((Tmpj1 > 0) && (TmpMonomial2[Tmpj1] > TmpMonomial2[Tmpj1 - 1]))
		  {
		    unsigned long Tmp = TmpMonomial2[Tmpj1 - 1];
		    TmpMonomial2[Tmpj1 - 1] = TmpMonomial2[Tmpj1];
		    TmpMonomial2[Tmpj1] = Tmp;
		    --Tmpj1;
		  }
                while ((Tmpj2 < ReducedNbrBosons) && (TmpMonomial2[Tmpj2] < TmpMonomial2[Tmpj2 + 1]))
                  {
                    unsigned long Tmp = TmpMonomial2[Tmpj2 + 1];
                    TmpMonomial2[Tmpj2 + 1] = TmpMonomial2[Tmpj2];
                    TmpMonomial2[Tmpj2] = Tmp;
                    ++Tmpj2;
                  }
		TmpState = this->ConvertFromMonomial(TmpMonomial2);
		if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		  {
		    long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
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
  delete[] TmpMonomial;
  cout << endl;

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

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->FermionHugeBasis->StateDescription[0l];
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  while (((MaxRoot >> TmpLzMax) & 0x1ul) == 0ul)
    --TmpLzMax;
  this->ConvertToMonomial(MaxRoot, TmpLzMax, TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
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
	this->ConvertToMonomial(CurrentPartition, TmpLzMax, TmpMonomial);
	if (i != this->FermionHugeBasis->FindStateIndexMemory(CurrentPartition, TmpMonomial[0] + ReducedNbrBosons))
	  cout << i << " " << this->FermionHugeBasis->FindStateIndexMemory(CurrentPartition, TmpMonomial[0] + ReducedNbrBosons) << endl;
	for (int j = 0; j < this->NbrBosons; ++j)
	  Rho += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
	double Coefficient = 0.0;
	for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	  for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	    {
	      double Diff = (double) (TmpMonomial[j1] - TmpMonomial[j2]);
	      unsigned int Max = TmpMonomial[j2];
	      unsigned long TmpState = 0x0ul;
	      int Tmpj1 = j1;
	      int Tmpj2 = j2;
	      for (int l = 0; l < this->NbrBosons; ++l)
		TmpMonomial2[l] = TmpMonomial[l];	    
	      for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		{
		  ++TmpMonomial2[Tmpj1];
		  --TmpMonomial2[Tmpj2];
		  Diff += 2.0;
		  while ((Tmpj1 > 0) && (TmpMonomial2[Tmpj1] > TmpMonomial2[Tmpj1 - 1]))
		    {
		      unsigned long Tmp = TmpMonomial2[Tmpj1 - 1];
		      TmpMonomial2[Tmpj1 - 1] = TmpMonomial2[Tmpj1];
		      TmpMonomial2[Tmpj1] = Tmp;
		      --Tmpj1;
		    }
		  while ((Tmpj2 < ReducedNbrBosons) && (TmpMonomial2[Tmpj2] < TmpMonomial2[Tmpj2 + 1]))
		    {
		      unsigned long Tmp = TmpMonomial2[Tmpj2 + 1];
		      TmpMonomial2[Tmpj2 + 1] = TmpMonomial2[Tmpj2];
		      TmpMonomial2[Tmpj2] = Tmp;
		      ++Tmpj2;
		    }
		  TmpState = this->ConvertFromMonomial(TmpMonomial2);
		  if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		    {
		      long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
		      if (TmpIndex < this->LargeHilbertSpaceDimension)
			Coefficient += Diff * jack[TmpIndex];
		    }
		}
	    }
	
	long TmpIndex = this->FermionHugeBasis->FindStateIndexMemory(this->FermionHugeBasis->GetSymmetricState(CurrentPartition), (this->LzMax - TmpMonomial2[ReducedNbrBosons]) + ReducedNbrBosons);
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
  delete[] TmpMonomial;
  cout << endl;

  return jack;
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// return value = converted state

RealVector& BosonOnSphereHaldaneHugeBasisShort::ConvertToUnnormalizedMonomial(RealVector& state, long reference)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0 / state[reference];
  state[reference] = 1.0;
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  unsigned long TmpState = this->FermionHugeBasis->StateDescription[reference];
   while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  this->ConvertToMonomial(TmpState, TmpLzMax, TmpMonomialReference);
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
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      TmpState = this->FermionHugeBasis->StateDescription[i];
      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	--TmpLzMax;
      this->ConvertToMonomial(TmpState, TmpLzMax, TmpMonomial);
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	{
	  while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrBosons)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrBosons)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomialReference[Index2]];
	  ++Index2;
	}
      Factorial = ReferenceFactorial;
      this->FermionToBoson(TmpState, TmpLzMax, 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	if (this->TemporaryState[k] > 1)
	  Factorial.FactorialMultiply(this->TemporaryState[k]);
      Coefficient *= sqrt(Factorial.GetNumericalValue());
      state[i] *= Coefficient;
    }
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// return value = converted state

RealVector& BosonOnSphereHaldaneHugeBasisShort::ConvertFromUnnormalizedMonomial(RealVector& state, long reference)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = state[reference];
  state[reference] = 1.0;
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  unsigned long TmpState = this->FermionHugeBasis->StateDescription[reference];
   while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
   this->ConvertToMonomial(TmpState, TmpLzMax, TmpMonomialReference);
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
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      TmpState = this->FermionHugeBasis->StateDescription[i];
      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	--TmpLzMax;
      this->ConvertToMonomial(TmpState, TmpLzMax, TmpMonomial);
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	{
	  while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrBosons)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
 	  ++Index1;
	}
      while (Index2 < this->NbrBosons)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomialReference[Index2]];
	  ++Index2;
	}
      Factorial = ReferenceFactorial;
      this->FermionToBoson(TmpState, TmpLzMax, 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	if (this->TemporaryState[k] > 1)
	  Factorial.FactorialDivide(this->TemporaryState[k]);
      Coefficient *= sqrt(Factorial.GetNumericalValue());
      state[i] *= Coefficient;
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
// return value = reference on the fused state

RealVector& BosonOnSphereHaldaneHugeBasisShort::FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
							    ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
							    bool symmetrizedFlag)
{
  BosonOnSphereShort* LeftSpace = (BosonOnSphereShort*) leftSpace;
  BosonOnSphereShort* RightSpace = (BosonOnSphereShort*) rightSpace;
  int StateShift = RightSpace->FermionBasis->LzMax + padding + 2;
  long Count = 0l;
  for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long TmpState1 = LeftSpace->FermionBasis->StateDescription[i] << StateShift;
      double Coefficient = leftVector[i];
      int TmpLzMax = this->FermionHugeBasis->LzMax;
      while ((TmpState1 >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      if (symmetrizedFlag == false)
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
  cout << "nbr of newly added components : " << Count << endl;
  return outputVector;
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
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = ((unsigned long) this->LzMax) >> 1;
  int TmpLzMax= this->FermionHugeBasis->LzMax;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
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
      SqrNorm += (outputVector[i] * outputVector[i]) * Factorial.GetNumericalValue();
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

