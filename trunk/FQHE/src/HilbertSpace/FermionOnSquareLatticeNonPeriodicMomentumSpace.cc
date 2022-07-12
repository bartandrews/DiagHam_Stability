////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                        class of fermions on square lattice                 //
//                  in momentum space with open boundary conditions           //
//                                                                            //
//                        last modification : 28/02/2011                      //
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
#include "HilbertSpace/FermionOnSquareLatticeNonPeriodicMomentumSpace.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"

#include <math.h>
#include <cstdlib>
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
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrAllowedKx = number of kx momenta allowed in the space
// nbrAllowedKy = number of kx momenta allowed in the space
// minKx = minimum value of kx momenta allowed in the space
// minKy = minimum value of ky momenta allowed in the space
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeNonPeriodicMomentumSpace::FermionOnSquareLatticeNonPeriodicMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int nbrAllowedKx, int nbrAllowedKy, int minKx, int minKy, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrAllowedKx = nbrAllowedKx;
  this->NbrAllowedKy = nbrAllowedKy;
  this->MinKx=minKx;
  this->MinKy=minKy;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrAllowedKx - 1, this->NbrAllowedKy - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateLzMax = new int [this->HilbertSpaceDimension];  
      this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrAllowedKx - 1, this->NbrAllowedKy - 1, 0, 0, 0l);
      int TmpLzMax = this->LzMax;
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  while ((this->StateDescription[i] >> TmpLzMax) == 0x0ul && TmpLzMax>0)
	    --TmpLzMax;
	  this->StateLzMax[i] = TmpLzMax;
	}
      this->GenerateLookUpTable(memory);      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
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
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSquareLatticeNonPeriodicMomentumSpace::FermionOnSquareLatticeNonPeriodicMomentumSpace(const FermionOnSquareLatticeNonPeriodicMomentumSpace& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->NbrAllowedKx = fermions.NbrAllowedKx;
  this->NbrAllowedKy = fermions.NbrAllowedKy;
  this->MinKx=fermions.MinKx;
  this->MinKy=fermions.MinKy;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnSquareLatticeNonPeriodicMomentumSpace::~FermionOnSquareLatticeNonPeriodicMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeNonPeriodicMomentumSpace& FermionOnSquareLatticeNonPeriodicMomentumSpace::operator = (const FermionOnSquareLatticeNonPeriodicMomentumSpace& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->NbrAllowedKx = fermions.NbrAllowedKx;
  this->NbrAllowedKy = fermions.NbrAllowedKy;
  this->MinKx=fermions.MinKx;
  this->MinKy=fermions.MinKy;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSquareLatticeNonPeriodicMomentumSpace::Clone()
{
  return new FermionOnSquareLatticeNonPeriodicMomentumSpace(*this);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKxShifted = current momentum along x for a single particle, shifted to fit in [0,NbrAllowedKx)
// currentKyShifted = current momentum along y for a single particle, shifted to fit in [0,NbrAllowedKy)
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeNonPeriodicMomentumSpace::GenerateStates(int nbrFermions, int currentKxShifted, int currentKyShifted, int currentTotalKx, int currentTotalKy, long pos)
{
  if (currentKyShifted < 0)
    {
      currentKyShifted = this->NbrAllowedKy - 1;
      currentKxShifted--;
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKxShifted < 0)
    return pos;
  int currentKx = (currentKxShifted + this->MinKx) % this->NbrSiteX;
  int currentKy = (currentKyShifted + this->MinKy) % this->NbrSiteY;
  if (nbrFermions == 1)
    {
      for (int jShifted = currentKyShifted; jShifted >= 0; --jShifted)
	{
          int j = (jShifted + this->MinKy) % this->NbrSiteY;
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      this->StateDescription[pos] = 0x1ul << ((currentKx * this->NbrSiteY) + j);
	      ++pos;
	    }
	}
      for (int iShifted = currentKxShifted - 1; iShifted >= 0; --iShifted)
	{
          int i = (iShifted + this->MinKx) % this->NbrSiteX;
	  for (int jShifted = this->NbrAllowedKy - 1; jShifted >= 0; --jShifted)
	    {
              int j = (jShifted + this->MinKy) % this->NbrSiteY;
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
 		  this->StateDescription[pos] = 0x1ul << ((i * this->NbrSiteY) + j);
 		  ++pos;
		}
	    }
	}
      return pos;
    }
  long TmpPos = this->GenerateStates(nbrFermions - 1, currentKxShifted, currentKyShifted - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
  unsigned long Mask = 0x1ul << ((currentKx * this->NbrSiteY) + currentKy);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentKxShifted, currentKyShifted - 1, currentTotalKx, currentTotalKy, pos);
};


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKxShifted = current momentum along x for a single particle, shifted to fit in [0,NbrAllowedKx)
// currentKyShifted = current momentum along y for a single particle, shifted to fit in [0,NbrAllowedKy)
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnSquareLatticeNonPeriodicMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKxShifted, int currentKyShifted, int currentTotalKx, int currentTotalKy)
{
  if (currentKyShifted < 0)
    {
      currentKyShifted = this->NbrAllowedKy - 1;
      currentKxShifted--;
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKxShifted < 0)
    return 0l;
  long Count = 0;
  int currentKx = (currentKxShifted + this->MinKx) % this->NbrSiteX;
  int currentKy = (currentKyShifted + this->MinKy) % this->NbrSiteY;
  if (nbrFermions == 1)
    {
      for (int jShifted = currentKyShifted; jShifted >= 0; --jShifted)
	{
          int j = (jShifted + this->MinKy) % this->NbrSiteY;
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    ++Count;
	}
      for (int iShifted = currentKxShifted - 1; iShifted >= 0; --iShifted)
	{
          int i = (iShifted + this->MinKx) % this->NbrSiteX;
	  for (int jShifted = this->NbrAllowedKy - 1; jShifted >= 0; --jShifted)
	    {
              int j = (jShifted + this->MinKy) % this->NbrSiteY;
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		++Count;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKxShifted, currentKyShifted - 1, currentTotalKx + currentKx, currentTotalKy + currentKy);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKxShifted, currentKyShifted - 1, currentTotalKx, currentTotalKy);
  return Count;
}

