////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system size such that            //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                        without a fixed total Lz value                      //
//                                                                            //
//                        last modification : 19/04/2010                      //
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
#include "HilbertSpace/BosonOnSphereFullShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h" 
#include "GeneralTools/StringTools.h"
#include "GeneralTools/ArrayTools.h"

#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


// default constructor
//

BosonOnSphereFullShort::BosonOnSphereFullShort ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// lzMax = maximum Lz value reached by a boson

BosonOnSphereFullShort::BosonOnSphereFullShort (int nbrBosons, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->FermionBasis = new FermionOnSphereFull(nbrBosons, lzMax + nbrBosons - 1);
  this->HilbertSpaceDimension = this->FermionBasis->HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = this->FermionBasis->LargeHilbertSpaceDimension;

  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->Flag.Initialize();
  int TmpLzMax = this->LzMax;
  this->TargetSpace = this;

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

BosonOnSphereFullShort::BosonOnSphereFullShort(const BosonOnSphereFullShort& bosons)
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

// destructor
//

BosonOnSphereFullShort::~BosonOnSphereFullShort ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereFullShort& BosonOnSphereFullShort::operator = (const BosonOnSphereFullShort& bosons)
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

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereFullShort::Clone()
{
  return new BosonOnSphereFullShort(*this);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereFullShort::AdA (int index, int m)
{
  return 0.0;
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereFullShort::AdA (long index, int m)
{
  return 0.0;
}

// evaluate a density matrix of a subsystem of the whole system described |left><right|. The density matrix is only evaluated for a fixed number of particles
// 
// densityMatrix = reference on the temporary storage for the reduced density matrix
// leftSpace = Hilbert space associated to the left state 
// leftState = reference on the left state
// rightSpace = Hilbert space associated to the right state 
// rightState = reference on the right state
// return value = reference on the reduced density matrix of the subsytem 

RealMatrix& BosonOnSphereFullShort::EvaluatePartialDensityMatrix (RealMatrix& densityMatrix, BosonOnSphereShort* leftSpace, RealVector& leftState, BosonOnSphereShort* rightSpace, RealVector& rightState)
{
  densityMatrix.ClearMatrix();
  if (this->LzMax <= 0)
    {
      if (this->NbrBosons == 0)
	{
	  densityMatrix.SetMatrixElement(0, 0, 1.0);
	  return densityMatrix;
	}
      else
	{
	  return densityMatrix;	  
	}
    }
  FermionOnSphereFull* TmpSpace = (FermionOnSphereFull*) this->FermionBasis;
  if (this->LzMax >= leftSpace->LzMax)
    {
      if (this->NbrBosons == leftSpace->NbrBosons)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (TmpSpace->TotalLzValues[i] == leftSpace->TotalLz)
	      {
		for (int j = 0; j < this->HilbertSpaceDimension; ++j)
		  {
		    if (TmpSpace->TotalLzValues[j] == rightSpace->TotalLz)
		      densityMatrix.SetMatrixElement(i, j, leftState[i] * rightState[j]);
		  }
	      }
	  return densityMatrix;	  
	}
      else
	{
	  return densityMatrix;	  
	}
    }
  
  long TmpNbrNonZeroElements = 0;
  int ShiftedLeftSpaceTotalLz = (leftSpace->TotalLz + (leftSpace->LzMax * leftSpace->NbrBosons)) >> 1;
  int ShiftedRightSpaceTotalLz = (rightSpace->TotalLz + (rightSpace->LzMax * rightSpace->NbrBosons)) >> 1;
  for (int CurrentLeftLzSector = 0; CurrentLeftLzSector < 0; CurrentLeftLzSector += 2)
    {
      int UnshiftedCurrentLeftLzSector = (CurrentLeftLzSector << 1) - (this->NbrBosons * this->LzMax);
      cout << "processing Lz=" <<  UnshiftedCurrentLeftLzSector << " sector" << endl;
      BosonOnSphereShort TmpLeftDestinationHilbertSpace(this->NbrBosons, UnshiftedCurrentLeftLzSector, this->LzMax);
      int ComplementarySubsystemLzMax = leftSpace->LzMax - this->LzMax - 1;
      int ComplementaryNbrBosons = leftSpace->NbrBosons - this->NbrBosons;
      BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosons, ((ShiftedLeftSpaceTotalLz - CurrentLeftLzSector) << 1) - (ComplementarySubsystemLzMax * ComplementaryNbrBosons), ComplementarySubsystemLzMax);
      int ShiftedCurrentRightLzSector = ShiftedRightSpaceTotalLz - CurrentLeftLzSector;
      BosonOnSphereShort TmpRightDestinationHilbertSpace(this->NbrBosons, (ShiftedCurrentRightLzSector << 1) - (this->NbrBosons * this->LzMax), this->LzMax);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  int Pos = 0;
	  unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (this->LzMax + 1 + this->NbrBosons);
	  for (int j = 0; j < TmpLeftDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpLeftState = TmpLeftDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	      int TmpLzMax = TmpLeftDestinationHilbertSpace.FermionBasis->LzMax;
	      while (((TmpLeftState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = TmpLeftDestinationHilbertSpace.FermionBasis->FindStateIndex(TmpLeftState, TmpLzMax);
	      if (TmpPos != TmpLeftDestinationHilbertSpace.HilbertSpaceDimension)
		{
		  double TmpCoefficient = leftState[TmpPos];
		  int TmpIndex = this->FermionBasis->FindStateIndex(TmpLeftDestinationHilbertSpace.FermionBasis->StateDescription[j], TmpLeftDestinationHilbertSpace.FermionBasis->StateLzMax[j]);
		  for (int k = 0; k < TmpRightDestinationHilbertSpace.HilbertSpaceDimension; ++k)
		    {
		      unsigned long TmpRightState = TmpRightDestinationHilbertSpace.FermionBasis->StateDescription[k] | TmpComplementaryState;
		      TmpLzMax = TmpRightDestinationHilbertSpace.FermionBasis->LzMax;
		      while (((TmpRightState >> TmpLzMax) & 0x1ul) == 0x0ul)
			--TmpLzMax;
		      int TmpPos2 = TmpRightDestinationHilbertSpace.FermionBasis->FindStateIndex(TmpRightState, TmpLzMax);
		      if (TmpPos2 != TmpRightDestinationHilbertSpace.HilbertSpaceDimension)
			{
			  int TmpIndex2 = this->FermionBasis->FindStateIndex(TmpRightDestinationHilbertSpace.FermionBasis->StateDescription[k], TmpRightDestinationHilbertSpace.FermionBasis->StateLzMax[k]);
			  densityMatrix.AddToMatrixElement(TmpIndex, TmpIndex2, TmpCoefficient * rightState[TmpPos2]);
			  ++TmpNbrNonZeroElements;
			}
		    }
		}
 	    }
 	}
     }
  return densityMatrix;
}
