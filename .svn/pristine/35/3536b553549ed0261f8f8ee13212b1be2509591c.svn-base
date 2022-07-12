////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of bosons on disk for system size such that             //
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
#include "HilbertSpace/BosonOnDiskShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/FactorialCoefficient.h" 

#include <math.h>
#include <iostream>

#ifdef __GSL__                                                                                                                                
#include <gsl/gsl_sf_gamma.h>                                                                                                                 
#endif

using std::cout;
using std::endl;


// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum angular momentum that a single particle can reach (negative if it has to be deduced from nbrBosons and totalLz)

BosonOnDiskShort::BosonOnDiskShort (int nbrBosons, int totalLz, int lzMax)
{
  this->TargetSpace = this;
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  if (lzMax < 0)
    this->LzMax = totalLz;
  else
    this->LzMax = lzMax;
  this->ShiftedTotalLz = this->TotalLz;
  this->TotalLz = (this->TotalLz << 1) - (this->NbrBosons * this->LzMax);
  this->NbrLzValue = this->LzMax + 1;
  this->FermionBasis = new FermionOnSphere(nbrBosons, this->TotalLz, this->LzMax  + nbrBosons - 1);
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

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnDiskShort::BosonOnDiskShort(const BosonOnDiskShort& bosons)
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
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
}

// copy constructor, preversing only some specific states 
// 
// bosons = reference on the hilbert space to copy to copy
// nbrPreservedStates = number of preserved states
// preservedStates = array of flags that indicates if the corresponding state has to be preserved 
//                   (dimension of the array should the one of the original Hilbert space)

BosonOnDiskShort::BosonOnDiskShort (const BosonOnDiskShort& bosons, int nbrPreservedStates, bool* preservedStates)
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
  this->HilbertSpaceDimension = nbrPreservedStates;
  this->LargeHilbertSpaceDimension = (long) nbrPreservedStates;
  this->Flag.Initialize();
  this->FermionBasis = new FermionOnSphere(*(bosons.FermionBasis), nbrPreservedStates, preservedStates);
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
}

// destructor
//

BosonOnDiskShort::~BosonOnDiskShort ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnDiskShort& BosonOnDiskShort::operator = (const BosonOnDiskShort& bosons)
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
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  return *this;
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& BosonOnDiskShort::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0 / state[reference];
  state[reference] = 1.0;
  this->ConvertToMonomial(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], TmpMonomialReference);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  SqrtCoefficients[0] = 1.0;
  InvSqrtCoefficients[0] = 1.0;
  for (int k = 1; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt((double) k) * SqrtCoefficients[k - 1];
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], 
		       this->TemporaryState, this->TemporaryStateLzMax);
  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialDivide(this->TemporaryState[k]);
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
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
      if (symmetryFactor == true)
	{
	  Factorial = ReferenceFactorial;
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialMultiply(this->TemporaryState[k]);
	  Coefficient *= sqrt(Factorial.GetNumericalValue());
	}
      state[i] *= Coefficient;
    }
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& BosonOnDiskShort::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
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
  this->ConvertFromUnnormalizedMonomialCore(state, SqrtCoefficients, InvSqrtCoefficients, reference, symmetryFactor);
  delete[] SqrtCoefficients;
  delete[] InvSqrtCoefficients;
  return state;
}

// convert a state such that its components are now expressed in the normalized basis, shifting all orbitals
//
// state = reference to the state to convert
// shift = shift to apply to each orbitals
// reference = set which component has been normalized to 1
// return value = converted state

RealVector& BosonOnDiskShort::ShiftedConvertFromUnnormalizedMonomial(RealVector& state, int shift, long reference)
{
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  SqrtCoefficients[0] = 1.0;
  InvSqrtCoefficients[0] = 1.0;
  for (int k = 1; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt((double) (shift + k)) * SqrtCoefficients[k - 1];
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  this->ConvertFromUnnormalizedMonomialCore(state, SqrtCoefficients, InvSqrtCoefficients, reference, true);
  delete[] SqrtCoefficients;
  delete[] InvSqrtCoefficients;
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

RealVector& BosonOnDiskShort::ConvertFromUnnormalizedMonomialCore(RealVector& state, double* sqrtCoefficients, double* invSqrtCoefficients, long reference, bool symmetryFactor)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0;
  this->ConvertToMonomial(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], TmpMonomialReference);
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], 
		       this->TemporaryState, this->TemporaryStateLzMax);
  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialMultiply(this->TemporaryState[k]);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	{
	  while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
	    {
	      Coefficient *= invSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
	    {
	      Coefficient *= sqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrBosons)
	{
	  Coefficient *= invSqrtCoefficients[TmpMonomialReference[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrBosons)
	{
	  Coefficient *= sqrtCoefficients[TmpMonomialReference[Index2]];
	  ++Index2;
	}
      if (symmetryFactor == true)
	{
	  Factorial = ReferenceFactorial;
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialDivide(this->TemporaryState[k]);
	  Coefficient *= sqrt(Factorial.GetNumericalValue());
	}
      state[i] *= Coefficient;
    }
  state /= state.Norm();
  delete[] TmpMonomialReference;
  delete[] TmpMonomial;
  return state;
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition. The entanglement matrix is only evaluated in a given Lz sector.
// and computed from precalculated particle entanglement matrix
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// radius = radius of the A disk
// shift = shift to apply to each orbitals
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& BosonOnDiskShort::EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrBosonSector, int lzSector, double radius, int shift, RealMatrix& entanglementMatrix)
{
#ifdef __GSL__
  double ReducedRadius = radius * radius / 2.0;
  double* IncompleteNormalizedGamma = new double[this->LzMax + 1];
  
  for (int i = 0; i <= this->LzMax ; i++)
    {
      IncompleteNormalizedGamma[i] = gsl_sf_gamma_inc_Q( ((double) i + shift) + 1.0 , ReducedRadius);
    }
  
  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;
  BosonOnDiskShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, this->LzMax);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
	
  int ComplementaryLz = (this->TotalLz + this->LzMax * this->NbrBosons)>>1;
  BosonOnDiskShort TmpHilbertSpace(ComplementaryNbrBosonSector, ComplementaryLz - lzSector, this->LzMax);
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpMonomial3);
      double Tmp = 0.0;     
      for( int j = 0; j < nbrBosonSector; j++)
	{
	  Tmp += log( IncompleteNormalizedGamma[TmpMonomial3[j]] );
	}
      Tmp = exp(0.5 * Tmp);
      for (int j = 0; j < TmpHilbertSpace.HilbertSpaceDimension; ++j)          
	entanglementMatrix(i, j) *= Tmp;      
    }
  
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
      double FormFactor = 0.0;
      for (int i=0; i < ComplementaryNbrBosonSector; i++)
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
