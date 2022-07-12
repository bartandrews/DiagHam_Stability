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
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/ArrayTools.h"
#include "Polynomial/RationalPolynomial.h"
#include "Polynomial/LongRationalPolynomial.h"
#include "Architecture/ArchitectureOperation/FQHESphereJackGeneratorSumRationalPolynomialOperation.h"

#include <math.h>


using std::cout;
using std::endl;


// default constructor
//

BosonOnSphereHaldaneBasisShort::BosonOnSphereHaldaneBasisShort ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson
// referenceState = array that describes the reference state to start from

BosonOnSphereHaldaneBasisShort::BosonOnSphereHaldaneBasisShort (int nbrBosons, int& totalLz, int lzMax, int* referenceState)
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
  this->TotalLz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      for (int j = 0; j < referenceState[i]; ++j)
	{
	  TmpReferenceState[TmpIndex] = 1;
	  ++TmpIndex;
	  this->TotalLz += i;
	}
      ++TmpIndex;
    }
  this->TotalLz = (2 * this->TotalLz) - (this->NbrBosons *  this->LzMax);
  this->FermionBasis = new FermionOnSphereHaldaneBasis(nbrBosons, totalLz, this->LzMax + nbrBosons - 1, TmpReferenceState);
  totalLz = this->TotalLz;
  delete[] TmpReferenceState;
  this->HilbertSpaceDimension = this->FermionBasis->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->FermionBasis->GetLargeHilbertSpaceDimension();

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
// memory = amount of memory granted for precalculations

BosonOnSphereHaldaneBasisShort::BosonOnSphereHaldaneBasisShort (char* fileName)
{
  this->FermionBasis = new FermionOnSphereHaldaneBasis(fileName);
  this->HilbertSpaceDimension = this->FermionBasis->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->FermionBasis->GetLargeHilbertSpaceDimension();

  this->NbrBosons = this->FermionBasis->NbrFermions;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = this->FermionBasis->TotalLz;
  this->LzMax = this->FermionBasis->LzMax - this->NbrBosons + 1;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;

  this->TargetSpace = this;
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

BosonOnSphereHaldaneBasisShort::BosonOnSphereHaldaneBasisShort(const BosonOnSphereHaldaneBasisShort& bosons)
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
  this->FermionBasis = (FermionOnSphereHaldaneBasis*) bosons.FermionBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
}

// destructor
//

BosonOnSphereHaldaneBasisShort::~BosonOnSphereHaldaneBasisShort ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereHaldaneBasisShort& BosonOnSphereHaldaneBasisShort::operator = (const BosonOnSphereHaldaneBasisShort& bosons)
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
  this->FermionBasis = (FermionOnSphereHaldaneBasis*) bosons.FermionBasis->Clone();
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];

  return *this;
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool BosonOnSphereHaldaneBasisShort::WriteHilbertSpace (char* fileName)
{
  return this->FermionBasis->WriteHilbertSpace(fileName);
}

// convert a given state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector BosonOnSphereHaldaneBasisShort::ConvertToNbodyBasis(RealVector& state, BosonOnSphereShort& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      TmpVector[nbodyBasis.FermionBasis->FindStateIndex(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i])] = state[i];
    }
  return TmpVector;
}

// convert a given state from the usual n-body basis to the Haldane basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector BosonOnSphereHaldaneBasisShort::ConvertFromNbodyBasis(RealVector& state, BosonOnSphereShort& nbodyBasis)
{
  RealVector TmpVector (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    TmpVector[i] = state[nbodyBasis.FermionBasis->FindStateIndex(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i])];
  TmpVector /= TmpVector.Norm();
  return TmpVector;
}

// convert a given state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

LongRationalVector BosonOnSphereHaldaneBasisShort::ConvertToNbodyBasis(LongRationalVector& state, BosonOnSphereShort& nbodyBasis)
{
  LongRationalVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      TmpVector[nbodyBasis.FermionBasis->FindStateIndex(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i])] = state[i];
    }
  return TmpVector;
}

// convert a given state from the usual n-body basis to the Haldane basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

LongRationalVector BosonOnSphereHaldaneBasisShort::ConvertFromNbodyBasis(LongRationalVector& state, BosonOnSphereShort& nbodyBasis)
{
  LongRationalVector TmpVector (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    TmpVector[i] = state[nbodyBasis.FermionBasis->FindStateIndex(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i])];
  TmpVector /= TmpVector.Norm();
  return TmpVector;
}

// create the Jack polynomial decomposition corresponding to the root partition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& BosonOnSphereHaldaneBasisShort::GenerateJackPolynomial(RealVector& jack, double alpha)
{
  jack[0] = 1.0;
  double InvAlpha =  2.0 / alpha;

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
  int ReducedNbrBosons = this->NbrBosons - 1;

  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (jack[i] == 0.0)
	{
	  double Rho = 0.0;
	  unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
	  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
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
			long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
			if (TmpIndex < this->HilbertSpaceDimension)
			  Coefficient += Diff * jack[TmpIndex];
		      }
		  }
	      }
	  jack[i] = Coefficient * InvAlpha / (RhoRoot - Rho);
	}
      else
	{
	  cout << "toto" << endl;
	}
      if ((i & 0xffl) == 0l)
	{
	  //	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  //	  cout.flush();
	}
    }
  delete[] TmpMonomial;
  cout << endl;

  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alphaNumerator = numerator of the Jack polynomial alpha coefficient
// alphaDenominator = numerator of the Jack polynomial alpha coefficient
// symbolicDepth = use symbolic calculation to solve singular values if non zero, if greater than zero it will use symbolic calculation up to a given depth (below that depth it will rely on numerical calculation),
//                 -1 if the symbolic calculation has to be done up to the point the singular value problem has been solved
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RationalVector& BosonOnSphereHaldaneBasisShort::GenerateJackPolynomial(RationalVector& jack, long alphaNumerator, long alphaDenominator, int symbolicDepth)
{
  jack[0] = 1l;
  Rational InvAlpha (2l * alphaDenominator, alphaNumerator);

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];

  Rational RhoRoot = 0l;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      RhoRoot += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
    }
  int ReducedNbrBosons = this->NbrBosons - 1;

  RationalPolynomial* TmpNumerators = 0;
  RationalPolynomial* TmpDenominators = 0;
  if (symbolicDepth != 0)
    {
      TmpNumerators = new RationalPolynomial[this->LargeHilbertSpaceDimension];
      TmpDenominators = new RationalPolynomial[this->LargeHilbertSpaceDimension];		  
    }
  long CounterMask = 0xffffl;
  if (symbolicDepth > 1)
    {
      CounterMask = 0xffl;
    }

  if (symbolicDepth != 0)
    {
      this->GenerateSingleJackPolynomialCoefficient(jack, 0, TmpNumerators, TmpDenominators);
    }
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (symbolicDepth > 1)
	{
	  this->GenerateSingleJackPolynomialCoefficient(jack, i, TmpNumerators, TmpDenominators);
	  if (symbolicDepth == 3)
	    {
	      cout << TmpNumerators[i] << endl;
	      cout << TmpDenominators[i] << endl;
	    }
	}
      if (jack[i].Num() == 0l)
	{
	  Rational Rho = 0l;
	  unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
	  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
	  for (int j = 0; j < this->NbrBosons; ++j)
	    Rho += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
	  if (Rho == RhoRoot)
	    {
 	      if (symbolicDepth == 0)
 		{
 		  cout << "warning : singular value detected at position " << i << ", skipping the rest of the calculation" << endl;
 		  return jack;
 		}
	      if (symbolicDepth != 3)
		{
		  cout << "singular value detected at position " << i << ", using symbolic calculation" << endl;
		}
	      this->GenerateSingleJackPolynomialCoefficient(jack, i, TmpNumerators, TmpDenominators);
	      Rational Tmp = TmpNumerators[i].PolynomialEvaluate(InvAlpha);
	      Tmp /= TmpDenominators[i].PolynomialEvaluate(InvAlpha);
	      if (symbolicDepth != 3)
		{
		  cout << "----------------------------------------------------------------" << endl
		       << "result = " << endl;
		  cout << TmpNumerators[i] << endl;
		  cout << TmpDenominators[i] << endl;
		  cout << Tmp << endl;
		}		  
	      jack[i] = Tmp;
	    }
	  else
	    {
	      Rational Coefficient = 0;
	      for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
		for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
		  {
		    long Diff = (long) (TmpMonomial[j1] - TmpMonomial[j2]);
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
			Diff += 2l;
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
			    long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
			    if (TmpIndex < this->HilbertSpaceDimension)
			      Coefficient += Diff * jack[TmpIndex];
			  }
		      }
		  }
	      jack[i] = (Coefficient * InvAlpha) / (RhoRoot - Rho);
	      //	      Rational TmpR = TmpNumerators[i].PolynomialEvaluate(InvAlpha) / TmpDenominators[i].PolynomialEvaluate(InvAlpha);
	      //	      cout << jack[i] << " " << TmpR << endl;
	      // cout << jack[i] << " ";
	      // this->PrintStateMonomial(cout, i) << endl;
	    }
	}
      if ((i & CounterMask) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  delete[] TmpMonomial;
  if (symbolicDepth != 0)
    {
      delete[] TmpNumerators;
      delete[] TmpDenominators;		  
    }
  cout << endl;
  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alphaNumerator = numerator of the Jack polynomial alpha coefficient
// alphaDenominator = numerator of the Jack polynomial alpha coefficient
// architecture = architecture to use for precalculation
// symbolicDepth = use symbolic calculation to solve singular values if non zero, if greater than zero it will use symbolic calculation up to a given depth (below that depth it will rely on numerical calculation),
//                 -1 if the symbolic calculation has to be done up to the point the singular value problem has been solved
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// fileName = optional file name to store temporary calculations
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

LongRationalVector& BosonOnSphereHaldaneBasisShort::GenerateJackPolynomial(LongRationalVector& jack, long alphaNumerator, long alphaDenominator, AbstractArchitecture* architecture, int symbolicDepth, long minIndex, long maxIndex, char* fileName)
{
  jack[0] = 1l;
  LongRational InvAlpha (2l * alphaDenominator, alphaNumerator);

  int ReducedNbrBosons = this->NbrBosons - 1;
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];
  long* ConnectedIndices = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients  = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedIndices2 = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients2  = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];

  long CounterMask = 0xffffl;
  if (symbolicDepth > 1)
    {
      CounterMask = 0xffl;
    }

  LongRational RhoRoot = 0l;
  LongRational Rho = 0l;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  long TmpRhoRootInvAlphaCoef = 0l;
  long TmpRhoRootConstCoef  = 0l;
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      TmpRhoRootInvAlphaCoef -= TmpMonomial[j] * ((long) j);
      TmpRhoRootConstCoef += TmpMonomial[j] * (TmpMonomial[j] - 1l);
      RhoRoot += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
    }
  LongRational RhoRootInvAlphaCoef(TmpRhoRootInvAlphaCoef);
  LongRational RhoRootConstCoef(TmpRhoRootConstCoef);

  LongRationalPolynomial* TmpNumerators = 0;
  LongRationalPolynomial* TmpDenominators = 0;		  
  int* EvaluatedCoeffcients = 0;
  if (symbolicDepth != 0)
    {
      TmpNumerators = new LongRationalPolynomial[this->LargeHilbertSpaceDimension];
      TmpDenominators = new LongRationalPolynomial[this->LargeHilbertSpaceDimension];		  
      EvaluatedCoeffcients = new int[this->LargeHilbertSpaceDimension];
      this->GenerateSingleJackPolynomialCoefficient(0, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2,
						    RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture);
    }
  LongRational Coefficient = 0l;
  LongRational Coefficient2 = 0l;
  
  if (minIndex < 1l)
    minIndex = 1l;
  if ((maxIndex <= minIndex) || (maxIndex > this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension;

  for (long i = minIndex; i < maxIndex; ++i)
    {
      if (symbolicDepth > 1)
	{
	  this->GenerateSingleJackPolynomialCoefficient(i, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2,
							RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture);
	  if (symbolicDepth == 3)
	    {
	      cout << TmpNumerators[i] << endl;
	      cout << TmpDenominators[i] << endl;
	    }
	}
      if (jack[i] == 0l)
	{
	  Rho = 0l;
	  unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
	  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
	  for (int j = 0; j < this->NbrBosons; ++j)
	    Rho += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
	  if (Rho == RhoRoot)
	    {
 	      if (symbolicDepth == 0)
 		{
 		  cout << "warning : singular value detected at position " << i << ", skipping the rest of the calculation" << endl;
 		  return jack;
 		}
	      if (symbolicDepth != 3)
		{
		  cout << "singular value detected at position " << i << ", using symbolic calculation" << endl;
		}
	      if (fileName != 0)
		jack.WriteVector(fileName);
	      for (long j = 0l; j < this->LargeHilbertSpaceDimension; ++j)
		{
		  if (TmpNumerators[j].Defined())
		    EvaluatedCoeffcients[j] = 1;
		  else
		    EvaluatedCoeffcients[j] = 0;
		}
	      long NbrComponentToEvaluate = this->GenerateSingleJackPolynomialCoefficientCountOnly(i, EvaluatedCoeffcients, TmpMonomial, TmpMonomial2, MaxRoot);
	      cout << "number of components to evalute using symbolic calculation : " << NbrComponentToEvaluate << endl;
	      this->GenerateSingleJackPolynomialCoefficient(i, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2,
							    RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture);
	      TmpNumerators[i].PolynomialEvaluate(InvAlpha, Coefficient);
	      Coefficient /= TmpDenominators[i].PolynomialEvaluate(InvAlpha);
	      if (symbolicDepth != 3)
		{
		  cout << "----------------------------------------------------------------" << endl
		       << "result = " << endl;
		  cout << TmpNumerators[i] << endl;
		  cout << TmpDenominators[i] << endl;
		  cout << Coefficient << endl;
		}		  
	      jack[i] = Coefficient;
	    }
	  else
	    {
	      Coefficient = 0l;
	      int Pos = 0;
	      for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
		for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
		  {
		    long Diff = (long) (TmpMonomial[j1] - TmpMonomial[j2]);
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
			Diff += 2l;
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
 			    long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
 			    if (TmpIndex < this->HilbertSpaceDimension)
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
	}
      if ((i & CounterMask) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  delete[] TmpMonomial;
  delete[] TmpMonomial2;
  delete[] ConnectedIndices;
  delete[] ConnectedCoefficients;
  delete[] ConnectedIndices2;
  delete[] ConnectedCoefficients2;
  cout << endl;
  if (symbolicDepth != 0)
    {
      delete[] TmpNumerators;
      delete[] TmpDenominators;		  
      delete[] EvaluatedCoeffcients;
    }
  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur and the resulting state is invariant under the Lz<->-Lz symmetry
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alphaNumerator = numerator of the Jack polynomial alpha coefficient
// alphaDenominator = numerator of the Jack polynomial alpha coefficient
// symbolicDepth = use symbolic calculation to solve singular values if non zero, if greater than zero it will use symbolic calculation up to a given depth (below that depth it will rely on numerical calculation),
//                 -1 if the symbolic calculation has to be done up to the point the singular value problem has been solved
// minIndex = start computing the Jack polynomial from the minIndex-th component
// maxIndex = stop  computing the Jack polynomial up to the maxIndex-th component (0 if it has to be computed up to the end)
// fileName = optional file name to store temporary calculations
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

LongRationalVector& BosonOnSphereHaldaneBasisShort::GenerateSymmetrizedJackPolynomial(LongRationalVector& jack, long alphaNumerator, long alphaDenominator, AbstractArchitecture* architecture, int symbolicDepth, long minIndex, long maxIndex, char* fileName)
{
  jack[0] = 1l;
  LongRational InvAlpha (2l * alphaDenominator, alphaNumerator);

  int ReducedNbrBosons = this->NbrBosons - 1;
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];
  long* ConnectedIndices = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients  = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedIndices2 = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients2  = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];

  long CounterMask = 0xffffl;
  if (symbolicDepth > 1)
    {
      CounterMask = 0xffl;
    }

  LongRational RhoRoot = 0l;
  LongRational Rho = 0l;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  long TmpRhoRootInvAlphaCoef = 0l;
  long TmpRhoRootConstCoef  = 0l;
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      TmpRhoRootInvAlphaCoef -= TmpMonomial[j] * ((long) j);
      TmpRhoRootConstCoef += TmpMonomial[j] * (TmpMonomial[j] - 1l);
      RhoRoot += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
    }
  LongRational RhoRootInvAlphaCoef(TmpRhoRootInvAlphaCoef);
  LongRational RhoRootConstCoef(TmpRhoRootConstCoef);

  LongRationalPolynomial* TmpNumerators = 0;
  LongRationalPolynomial* TmpDenominators = 0;
  int* EvaluatedCoeffcients = 0;
  if (symbolicDepth != 0)
    {
      TmpNumerators = new LongRationalPolynomial[this->LargeHilbertSpaceDimension];
      TmpDenominators = new LongRationalPolynomial[this->LargeHilbertSpaceDimension];		  
      EvaluatedCoeffcients = new int[this->LargeHilbertSpaceDimension];
      this->GenerateSingleJackPolynomialCoefficient(0, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2, 
						    RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture);
    }
  LongRational Coefficient = 0l;
  LongRational Coefficient2 = 0l;

  if (minIndex < 1l)
    minIndex = 1l;
  if ((maxIndex <= minIndex) || (maxIndex > this->LargeHilbertSpaceDimension))
    maxIndex = this->LargeHilbertSpaceDimension;
  
  for (long i = minIndex; i < maxIndex; ++i)
    {
      if (jack[i] == 0l)
	{
	  Rho = 0l;
	  unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
	  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
	  for (int j = 0; j < this->NbrBosons; ++j)
	    Rho += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
	  if (symbolicDepth > 1)
	    {
	      this->GenerateSingleJackPolynomialCoefficient(i, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2, 
							    RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture, 0, true);
	      long TmpIndex = this->FermionBasis->FindStateIndex(((FermionOnSphereHaldaneBasis*) this->FermionBasis)->GetSymmetricState(CurrentPartition), (this->LzMax - TmpMonomial[ReducedNbrBosons]) + ReducedNbrBosons);
	      if (i < TmpIndex)
		{
		  TmpNumerators[TmpIndex] = TmpNumerators[i];
		  TmpDenominators[TmpIndex] = TmpDenominators[i];
		}
	      if (symbolicDepth == 3)
		{
		  cout << TmpNumerators[i] << endl;
		  cout << TmpDenominators[i] << endl;
		}
	    }
	  if (Rho == RhoRoot)
	    {
 	      if (symbolicDepth == 0)
 		{
 		  cout << "warning : singular value detected at position " << i << ", skipping the rest of the calculation" << endl;
 		  return jack;
 		}
	      if (symbolicDepth != 3)
		{
		  cout << "singular value detected at position " << i << ", using symbolic calculation" << endl;
		}
	      if (fileName != 0)
		jack.WriteVector(fileName);
	      for (long j = 0l; j < this->LargeHilbertSpaceDimension; ++j)
		{
		  if (TmpNumerators[j].Defined())
		    EvaluatedCoeffcients[j] = 1;
		  else
		    EvaluatedCoeffcients[j] = 0;
		}
	      long NbrComponentToEvaluate = this->GenerateSingleJackPolynomialCoefficientCountOnly(i, EvaluatedCoeffcients, TmpMonomial, TmpMonomial2, MaxRoot);
	      cout << "number of components to evalute using symbolic calculation : " << NbrComponentToEvaluate << endl;
	      this->GenerateSingleJackPolynomialCoefficient(i, TmpNumerators, TmpDenominators, ConnectedIndices, ConnectedCoefficients, TmpMonomial, TmpMonomial2, 
							    RhoRootInvAlphaCoef, RhoRootConstCoef, MaxRoot, architecture, 0, true);
	      LongRational Tmp = TmpNumerators[i].PolynomialEvaluate(InvAlpha);
	      Tmp /= TmpDenominators[i].PolynomialEvaluate(InvAlpha);
	      if (symbolicDepth != 3)
		{
		  cout << "----------------------------------------------------------------" << endl
		       << "result = " << endl;
		  cout << TmpNumerators[i] << endl;
		  cout << TmpDenominators[i] << endl;
		  cout << Tmp << endl;
		}		  
	      long TmpIndex = this->FermionBasis->FindStateIndex(((FermionOnSphereHaldaneBasis*) this->FermionBasis)->GetSymmetricState(CurrentPartition), (this->LzMax - TmpMonomial[ReducedNbrBosons]) + ReducedNbrBosons);
	      jack[i] = Tmp;
	      if (i < TmpIndex)
		{
		  jack[TmpIndex] = Tmp;
		  TmpNumerators[TmpIndex] = TmpNumerators[i];
		  TmpDenominators[TmpIndex] = TmpDenominators[i];
		}
	    }
	  else
	    {
	      Coefficient = 0l;
	      int Pos = 0;
	      for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
		for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
		  {
		    long Diff = (long) (TmpMonomial[j1] - TmpMonomial[j2]);
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
			Diff += 2l;
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
			    long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
			    if (TmpIndex < this->HilbertSpaceDimension)
			      {
				long TmpIndex2 = this->FermionBasis->FindStateIndex(((FermionOnSphereHaldaneBasis*) this->FermionBasis)->GetSymmetricState(TmpState), (this->LzMax - TmpMonomial2[ReducedNbrBosons]) + ReducedNbrBosons);
				if (TmpIndex < TmpIndex2)
				  ConnectedIndices[Pos] = TmpIndex;
				else
				  ConnectedIndices[Pos] = TmpIndex2;
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
	      long TmpIndex = this->FermionBasis->FindStateIndex(((FermionOnSphereHaldaneBasis*) this->FermionBasis)->GetSymmetricState(CurrentPartition), (this->LzMax - TmpMonomial[ReducedNbrBosons]) + ReducedNbrBosons);
	      Coefficient *= InvAlpha;
	      Rho -= RhoRoot;
	      Rho.Neg();
	      Coefficient /= Rho;
	      if (i < TmpIndex)
		jack[TmpIndex] = Coefficient;
	      jack[i] = Coefficient;
	    }
	}
      if ((i & CounterMask) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  delete[] TmpMonomial;
  delete[] TmpMonomial2;
  delete[] ConnectedIndices;
  delete[] ConnectedCoefficients;
  delete[] ConnectedIndices2;
  delete[] ConnectedCoefficients2;
  if (symbolicDepth != 0)
    {
      delete[] TmpNumerators;
      delete[] TmpDenominators;		  
      delete[] EvaluatedCoeffcients;
    }
  cout << endl;
  return jack;
}

// compute a single coefficient of the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur and using (partial symbolic calculation)
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// index = index of the component to compute
// numerators = array of polynomials attached to each coefficient numerator
// denominators = array of polynomials attached to each coefficient denominator
// return value = true if a fully symbolic calculation has been performed

bool BosonOnSphereHaldaneBasisShort::GenerateSingleJackPolynomialCoefficient(RationalVector& jack, long index, RationalPolynomial* numerators, RationalPolynomial* denominators)
{
  if (numerators[index].Defined())
    return true;
  if (index == 0l)
    {
      if (!numerators[0l].Defined())
	{
	  numerators[0l] = RationalPolynomial(0);
	  denominators[0l] = RationalPolynomial(0);
	  numerators[0l][0] = 1l;
	  denominators[0l][0] = 1l;
	  return true;
	}
    }

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];
  
  Rational RhoRootInvAlphaCoef = 0l;
  Rational RhoRootConstCoef = 0l;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      RhoRootInvAlphaCoef -= TmpMonomial[j] * ((long) j);
      RhoRootConstCoef += TmpMonomial[j] * (TmpMonomial[j] - 1l);
    }
  int ReducedNbrBosons = this->NbrBosons - 1;
  
  Rational RhoInvAlphaCoef = 0l;
  Rational RhoConstCoef = 0l;

  unsigned long CurrentPartition = this->FermionBasis->StateDescription[index];
  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[index], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      RhoInvAlphaCoef -= TmpMonomial[j] * ((long) j);
      RhoConstCoef += TmpMonomial[j] * (TmpMonomial[j] - 1l);
    }
  
  int Pos = 0;
  long* ConnectedIndices = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients  = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
    for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
      {
	long Diff = (long) (TmpMonomial[j1] - TmpMonomial[j2]);
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
	    Diff += 2l;
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
		long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
		if (TmpIndex < this->HilbertSpaceDimension)
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
      long* TmpConnectedIndices = new long [NbrConnected];
      long* TmpConnectedCoefficients  = new long [NbrConnected];
      TmpConnectedIndices[0] = ConnectedIndices[0];
      TmpConnectedCoefficients[0] = ConnectedCoefficients[0];
      TmpIndex = 1;
      NbrConnected = 1;
      while (TmpIndex < Pos)
	{
	  while ((TmpIndex < Pos) && (ConnectedIndices[TmpIndex] == ConnectedIndices[TmpIndex - 1]))
	    {
	      TmpConnectedCoefficients[NbrConnected - 1] += ConnectedCoefficients[TmpIndex];
	      ++TmpIndex;
	    }
	  if (TmpIndex < Pos)
	    {
	      TmpConnectedIndices[NbrConnected] = ConnectedIndices[TmpIndex];
	      TmpConnectedCoefficients[NbrConnected] = ConnectedCoefficients[TmpIndex];	   
	      ++NbrConnected;
	    }
	  ++TmpIndex;
	}
      delete[] ConnectedIndices;
      delete[] ConnectedCoefficients;
      ConnectedIndices = TmpConnectedIndices;
      ConnectedCoefficients = TmpConnectedCoefficients;
    }
  bool SymbolicFlag = false;

  SymbolicFlag = true;
  for (int i = 0; i < NbrConnected; ++i)
    {
      if (!numerators[ConnectedIndices[i]].Defined())
        {
	  this->GenerateSingleJackPolynomialCoefficient(jack, ConnectedIndices[i], numerators, denominators);
	}

    }


  numerators[index] = numerators[ConnectedIndices[0]];
  denominators[index] = denominators[ConnectedIndices[0]];
  numerators[index] *= ConnectedCoefficients[0];
  for (int i = 1; i < NbrConnected; ++i)
    {
      RationalPolynomial TmpNumerator = numerators[ConnectedIndices[i]];
      TmpNumerator *= denominators[index];
      TmpNumerator *= ConnectedCoefficients[i];
      numerators[index] *= denominators[ConnectedIndices[i]];
      numerators[index] += TmpNumerator;
      denominators[index] *= denominators[ConnectedIndices[i]];
      for (int j = 0; j < denominators[index].GetPolynomialDegree(); ++j)
	{
	  Rational Tmp4 = numerators[index].PolynomialEvaluate(denominators[index].PolynomialRationalRoot(j));	
	  if (Tmp4.Num() == 0l)
	    {
	      numerators[index] = numerators[index].MonomialDivision(denominators[index].PolynomialRationalRoot(j));
	      denominators[index] = denominators[index].MonomialDivision(denominators[index].PolynomialRationalRoot(j));
	      --j;
	    }
 	}
    }

  numerators[index].ShiftPowers(1);

  Rational Tmp2 = (RhoRootInvAlphaCoef - RhoInvAlphaCoef);

  if (Tmp2.Num() != 0l)
    {
      numerators[index] /= Tmp2;
      Rational Tmp3 = RhoRootConstCoef - RhoConstCoef;
      Tmp3 /= -Tmp2;
      Rational Tmp4 = numerators[index].PolynomialEvaluate(Tmp3);
      if (Tmp4.Num() == 0l)
	{
	  numerators[index] = numerators[index].MonomialDivision(Tmp3);
	}
      else
	{
	  denominators[index] = denominators[index].MonomialMultiplication(Tmp3);
	}
    }
  else
    {
      Tmp2 = (RhoRootConstCoef - RhoConstCoef);
      if (Tmp2.Num() != 0)
	{
	  numerators[index] /= Tmp2;
	}
    }

  delete[] ConnectedIndices;
  delete[] ConnectedCoefficients;
  delete[] TmpMonomial;
  delete[] TmpMonomial2;
  return SymbolicFlag;
}

// compute a single coefficient of the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur and using (partial symbolic calculation)
//
// index = index of the component to compute
// numerators = array of polynomials attached to each coefficient numerator
// denominators = array of polynomials attached to each coefficient denominator
// tmpMonomial = temporary array for monomial description
// tmpMonomial2 = temporary array for monomial description
// rhoRootInvAlphaCoef = coefficient in front of inv alpha in the rho for the root partition
// rhoRootConstCoef = constant coefficient in the rho for the root partition
// maxRoot = fermionic expression for the root partition
// architecture = architecture to use for precalculation
// currentNbrComponents = current number of components computed for the symbolic of the index-th compoment
// symmetryFlag = true iof the Lz <-> -Lz symmetry has to be used
// return value = total number of components computed for the symbolic of the index-th compoment

long BosonOnSphereHaldaneBasisShort::GenerateSingleJackPolynomialCoefficient(long index, LongRationalPolynomial* numerators, LongRationalPolynomial* denominators, long* connectedIndices, long* connectedCoefficients, unsigned long* tmpMonomial, unsigned long* tmpMonomial2, LongRational& rhoRootInvAlphaCoef, LongRational& rhoRootConstCoef, unsigned long maxRoot, AbstractArchitecture* architecture, long currentNbrComponents, bool symmetryFlag)
{
  if (numerators[index].Defined())
    return 0l;
  if (index == 0l)
    {
      if (!numerators[0l].Defined())
	{
	  numerators[0l] = LongRationalPolynomial(0);
	  denominators[0l] = LongRationalPolynomial(0);
	  numerators[0l][0] = 1l;
	  denominators[0l][0] = 1l;
	  return 0l;
	}
    }


  int ReducedNbrBosons = this->NbrBosons - 1;  

  unsigned long CurrentPartition = this->FermionBasis->StateDescription[index];
  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[index], tmpMonomial);
  long TmpRhoInvAlphaCoef = 0l;
  long TmpRhoConstCoef  = 0l;
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      TmpRhoInvAlphaCoef -= tmpMonomial[j] * ((long) j);
      TmpRhoConstCoef += tmpMonomial[j] * (tmpMonomial[j] - 1l);
    }
  LongRational RhoInvAlphaCoef(TmpRhoInvAlphaCoef);
  LongRational RhoConstCoef (TmpRhoConstCoef);
  
  int Pos = 0;
  for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
    for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
      {
	long Diff = (long) (tmpMonomial[j1] - tmpMonomial[j2]);
	unsigned int Max = tmpMonomial[j2];
	unsigned long TmpState = 0x0ul;
	int Tmpj1 = j1;
	int Tmpj2 = j2;
	for (int l = 0; l < this->NbrBosons; ++l)
	  tmpMonomial2[l] = tmpMonomial[l];	    
	for (unsigned int k = 1; (k <= Max) && (TmpState < maxRoot); ++k)
	  {
	    ++tmpMonomial2[Tmpj1];
	    --tmpMonomial2[Tmpj2];
	    Diff += 2l;
	    while ((Tmpj1 > 0) && (tmpMonomial2[Tmpj1] > tmpMonomial2[Tmpj1 - 1]))
	      {
		unsigned long Tmp = tmpMonomial2[Tmpj1 - 1];
		tmpMonomial2[Tmpj1 - 1] = tmpMonomial2[Tmpj1];
		tmpMonomial2[Tmpj1] = Tmp;
		--Tmpj1;
	      }
	    while ((Tmpj2 < ReducedNbrBosons) && (tmpMonomial2[Tmpj2] < tmpMonomial2[Tmpj2 + 1]))
	      {
		unsigned long Tmp = tmpMonomial2[Tmpj2 + 1];
		tmpMonomial2[Tmpj2 + 1] = tmpMonomial2[Tmpj2];
		tmpMonomial2[Tmpj2] = Tmp;
		++Tmpj2;
	      }
	    TmpState = this->ConvertFromMonomial(tmpMonomial2);
	    if ((TmpState <= maxRoot) && (TmpState > CurrentPartition))
	      {
		long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, tmpMonomial2[0] + ReducedNbrBosons);
		if (TmpIndex < this->HilbertSpaceDimension)
		  {
		    if (symmetryFlag == false)
		      {
			connectedIndices[Pos] = TmpIndex;
			connectedCoefficients[Pos] = Diff;
			++Pos;
		      }
		    else
		      {
			long TmpIndex2 = this->FermionBasis->FindStateIndex(((FermionOnSphereHaldaneBasis*) this->FermionBasis)->GetSymmetricState(TmpState), (this->LzMax - tmpMonomial2[ReducedNbrBosons]) + ReducedNbrBosons);
			if (TmpIndex < TmpIndex2)
			  connectedIndices[Pos] = TmpIndex;
			else
			  connectedIndices[Pos] = TmpIndex2;
			connectedCoefficients[Pos] = Diff;
			++Pos;
		      }
		  }
	      }
	  }
      }
  long* ConnectedIndices2 = 0;
  long* ConnectedCoefficients2  = 0;
  int NbrConnected = 1l;
  if (Pos > 1)
    {
      SortArrayDownOrdering<long>(connectedIndices, connectedCoefficients, Pos);
      int TmpIndex = 1;
      while (TmpIndex < Pos)
	{
	  while ((TmpIndex < Pos) && (connectedIndices[TmpIndex] == connectedIndices[TmpIndex - 1]))
	    ++TmpIndex;
	  if (TmpIndex < Pos)
	    ++NbrConnected;
	  ++TmpIndex;
	}
      ConnectedIndices2 = new long[NbrConnected];
      ConnectedCoefficients2  = new long[NbrConnected];
      ConnectedIndices2[0] = connectedIndices[0];
      ConnectedCoefficients2[0] = connectedCoefficients[0];
      TmpIndex = 1;
      NbrConnected = 1;
      while (TmpIndex < Pos)
	{
	  while ((TmpIndex < Pos) && (connectedIndices[TmpIndex] == connectedIndices[TmpIndex - 1]))
	    {
	      ConnectedCoefficients2[NbrConnected - 1] += connectedCoefficients[TmpIndex];
	      ++TmpIndex;
	    }
	  if (TmpIndex < Pos)
	    {
	      ConnectedIndices2[NbrConnected] = connectedIndices[TmpIndex];
	      ConnectedCoefficients2[NbrConnected] = connectedCoefficients[TmpIndex];	   
	      ++NbrConnected;
	    }
	  ++TmpIndex;
	}
    }
  else
    {
      ConnectedIndices2 = new long[NbrConnected];
      ConnectedCoefficients2  = new long[NbrConnected];
      ConnectedIndices2[0] = connectedIndices[0];
      ConnectedCoefficients2[0] = connectedCoefficients[0];
    }

  int TmpNbrComponents = 0;
  for (int i = 0; i < NbrConnected; ++i)
    {
      if (!numerators[ConnectedIndices2[i]].Defined())
        {
	  TmpNbrComponents += this->GenerateSingleJackPolynomialCoefficient(ConnectedIndices2[i], numerators, denominators, connectedIndices, connectedCoefficients, 
									    tmpMonomial, tmpMonomial2, rhoRootInvAlphaCoef, rhoRootConstCoef, maxRoot, architecture, currentNbrComponents + TmpNbrComponents, symmetryFlag);
	}
    }

  cout << (currentNbrComponents + TmpNbrComponents) << " symbolic components computed                  \r";
  cout.flush();

  LongRational Tmp4;
  if (NbrConnected > 1)
    {
      FQHESphereJackGeneratorSumRationalPolynomialOperation Operation(index, numerators, denominators, ConnectedIndices2, ConnectedCoefficients2, NbrConnected);
      Operation.ApplyOperation(architecture);
    }
  else
    {
      numerators[index]= numerators[ConnectedIndices2[0]];
      denominators[index] = denominators[ConnectedIndices2[0]];
      numerators[index] *= ConnectedCoefficients2[0];
    }
  numerators[index].ShiftPowers(1);

  LongRational Tmp2 = rhoRootInvAlphaCoef;
  Tmp2 -= RhoInvAlphaCoef;
  if (Tmp2 != 0l)
    {
      numerators[index] /= Tmp2;
      LongRational Tmp3 = rhoRootConstCoef;
      Tmp3 -= RhoConstCoef;
      Tmp3 /= -Tmp2;
      numerators[index].PolynomialEvaluate(Tmp3, Tmp4);
      if (Tmp4 == 0l)
	{
	  numerators[index].LocalMonomialDivision(Tmp3);
	}
      else
	{
	  denominators[index].LocalMonomialMultiplication(Tmp3);
	}
    }
  else
    {
      Tmp2 = rhoRootConstCoef;
      Tmp2 -= RhoConstCoef;
      if (Tmp2 != 0)
	{
	  numerators[index] /= Tmp2;
	}
    }

  delete[] ConnectedIndices2;
  delete[] ConnectedCoefficients2;

  return (TmpNbrComponents + 1);
}

// compute of many coefficients have to be computed to get a single coefficient of the Jack polynomial decomposition 
//
// index = index of the component to compute
// evaluatedCoeffcients = that indicates which coefficients have already been computed
// tmpMonomial = temporary array for monomial description
// tmpMonomial2 = temporary array for monomial description
// maxRoot = fermionic expression for the root partition
// return value = true if a fully symbolic calculation has been performed

long BosonOnSphereHaldaneBasisShort::GenerateSingleJackPolynomialCoefficientCountOnly(long index, int* evaluatedCoeffcients, unsigned long* tmpMonomial, unsigned long* tmpMonomial2, unsigned long maxRoot)
{
  if (index == 0l)
    return 1l;
  int ReducedNbrBosons = this->NbrBosons - 1;  

  unsigned long CurrentPartition = this->FermionBasis->StateDescription[index];
  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[index], tmpMonomial);

  long* ConnectedIndices2 = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  int NbrConnected  = 0;
  for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
    for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
      {
	unsigned int Max = tmpMonomial[j2];
	unsigned long TmpState = 0x0ul;
	int Tmpj1 = j1;
	int Tmpj2 = j2;
	for (int l = 0; l < this->NbrBosons; ++l)
	  tmpMonomial2[l] = tmpMonomial[l];	    
	for (unsigned int k = 1; (k <= Max) && (TmpState < maxRoot); ++k)
	  {
	    ++tmpMonomial2[Tmpj1];
	    --tmpMonomial2[Tmpj2];
	    while ((Tmpj1 > 0) && (tmpMonomial2[Tmpj1] > tmpMonomial2[Tmpj1 - 1]))
	      {
		unsigned long Tmp = tmpMonomial2[Tmpj1 - 1];
		tmpMonomial2[Tmpj1 - 1] = tmpMonomial2[Tmpj1];
		tmpMonomial2[Tmpj1] = Tmp;
		--Tmpj1;
	      }
	    while ((Tmpj2 < ReducedNbrBosons) && (tmpMonomial2[Tmpj2] < tmpMonomial2[Tmpj2 + 1]))
	      {
		unsigned long Tmp = tmpMonomial2[Tmpj2 + 1];
		tmpMonomial2[Tmpj2 + 1] = tmpMonomial2[Tmpj2];
		tmpMonomial2[Tmpj2] = Tmp;
		++Tmpj2;
	      }
	    TmpState = this->ConvertFromMonomial(tmpMonomial2);
	    if ((TmpState <= maxRoot) && (TmpState > CurrentPartition))
	      {
		long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, tmpMonomial2[0] + ReducedNbrBosons);
		if (TmpIndex < this->HilbertSpaceDimension)
		  {
		    ConnectedIndices2[NbrConnected] = TmpIndex;
		    ++NbrConnected;
		  }
	      }
	  }
      }

  long TmpNbrConnected = 1;
  for (int i = 0; i < NbrConnected; ++i)
    {
      if (evaluatedCoeffcients[ConnectedIndices2[i]] == 0)
        {
 	  TmpNbrConnected += this->GenerateSingleJackPolynomialCoefficientCountOnly(ConnectedIndices2[i], evaluatedCoeffcients, 
 										    tmpMonomial, tmpMonomial2, maxRoot);
	  evaluatedCoeffcients[ConnectedIndices2[i]] = 1;
	}
    }

  delete[] ConnectedIndices2;

  return TmpNbrConnected;
}



// create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& BosonOnSphereHaldaneBasisShort::GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha)
{
  jack[0] = 1.0;
  double InvAlpha =  2.0 / alpha;

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
  int ReducedNbrBosons = this->NbrBosons - 1;
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (jack[i] == 0.0)
	{
	  double Rho = 0.0;
	  unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
	  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
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
			long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
			if (TmpIndex < this->HilbertSpaceDimension)
			  Coefficient += Diff * jack[TmpIndex];
		      }
		  }
	      }
	  
	  long TmpIndex = this->FermionBasis->FindStateIndex(((FermionOnSphereHaldaneBasis*) this->FermionBasis)->GetSymmetricState(CurrentPartition), (this->LzMax - TmpMonomial[ReducedNbrBosons]) + ReducedNbrBosons);
	  Coefficient *= InvAlpha;
	  Coefficient /= (RhoRoot - Rho);
	  if (i < TmpIndex)
	    jack[TmpIndex] = Coefficient;
	  jack[i] = Coefficient;
	}
      if ((i & 0xffffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
      }
  delete[] TmpMonomial;
  cout << endl;

  return jack;
}

// check partitions that may lead to singular coefficient in a given Jack polynomial decomposition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// error = error when comparing two rho values
// return value = vector with non-zero component being rho factor of possible singular coefficients

RealVector& BosonOnSphereHaldaneBasisShort::CheckPossibleSingularCoefficientsInJackPolynomial(RealVector& jack, double alpha, double error)
{
  double InvAlpha =  2.0 / alpha;

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));

  jack[0] = RhoRoot;
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      double Rho = 0.0;
      unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
      for (int j = 0; j < this->NbrBosons; ++j)
	Rho += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
      if ((fabs(RhoRoot - Rho) < error) || (fabs(RhoRoot - Rho) < (error * fabs(RhoRoot))))
	jack[i] = Rho;
      else
	jack[i] = 0.0;
      }
  delete[] TmpMonomial;

  return jack;
}
  
// check partitions that may lead to singular coefficient in a given Jack polynomial decomposition, assuming only rational numbers occur
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alphaNumerator = numerator of the Jack polynomial alpha coefficient
// alphaDenominator = numerator of the Jack polynomial alpha coefficient
// checkConnectivity = if true, compute how many componets are involved in the calculation of a given singular coefficients
// return value = vector with non-zero component being rho factor of possible singular coefficients

RationalVector& BosonOnSphereHaldaneBasisShort::CheckPossibleSingularCoefficientsInJackPolynomial(RationalVector& jack, long alphaNumerator, long alphaDenominator, bool checkConnectivity)
{
  Rational InvAlpha (2l * alphaDenominator, alphaNumerator);

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];
  int* EvaluatedCoeffcients = 0;
  if (checkConnectivity == true)
    {
      EvaluatedCoeffcients = new int[this->LargeHilbertSpaceDimension];
    }

  Rational RhoRoot = 0;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));

  jack[0] = RhoRoot;
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      Rational Rho = 0l;
      unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
      for (int j = 0; j < this->NbrBosons; ++j)
	Rho += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
      if (RhoRoot == Rho)
	{
	  if (checkConnectivity == false)
	    {
	      jack[i] = Rho;
	    }
	  else
	    {
	      for (long j = 0l; j < this->LargeHilbertSpaceDimension; ++j)
		EvaluatedCoeffcients[j] = 0;
	      jack[i] = this->GenerateSingleJackPolynomialCoefficientCountOnly(i, EvaluatedCoeffcients, TmpMonomial, TmpMonomial2, MaxRoot);	      
	    }
	}
      else
	jack[i] = 0l;
      }
  delete[] TmpMonomial;
  delete[] TmpMonomial2;
  if (EvaluatedCoeffcients != 0)
    {
      delete[] EvaluatedCoeffcients;
    }
  return jack;
}
  
// check partitions that may lead to singular coefficient in a given Jack polynomial decomposition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// error = error when comparing two rho values
// return value = vector with non-zero component being rho factor of possible singular coefficients

void BosonOnSphereHaldaneBasisShort::CheckMaximumConnectedStateInJackPolynomial()
{
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  int ReducedNbrBosons = this->NbrBosons - 1;
  long LargestDistance = 0;
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
      long MinIndex = this->LargeHilbertSpaceDimension;
      for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	  {
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
		    long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
		    if ((TmpIndex < this->HilbertSpaceDimension) && (TmpIndex < MinIndex))
		      {
			MinIndex = TmpIndex;
		      }
		  }
	      }
	  }
      //      this->PrintState(cout, i) << " : " << MinIndex << endl;
      if ((MinIndex < this->LargeHilbertSpaceDimension) && (LargestDistance < (i - MinIndex)))
	{
	  LargestDistance = (i - MinIndex);
	}
    }
  cout << "largest distance = " << LargestDistance << "   (" << ((LargestDistance * 100.0) / ((double) this->LargeHilbertSpaceDimension))<< "%,  dim = " << this->LargeHilbertSpaceDimension << ")" << endl; 
}
  
