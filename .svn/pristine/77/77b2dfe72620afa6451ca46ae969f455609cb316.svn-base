////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on hyper cubic 4D lattice with spin         //
//                                  in momentum space                         //
//                                                                            //
//                        last modification : 12/07/2011                      //
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
#include "HilbertSpace/FermionOnHyperCubicLatticeWithSpinMomentumSpace.h"
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
// nbrSiteZ = number of sites in the z direction
// nbrSiteT = number of sites in the t direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// ktMomentum = momentum along the t direction
// memory = amount of memory granted for precalculations

FermionOnHyperCubicLatticeWithSpinMomentumSpace::FermionOnHyperCubicLatticeWithSpinMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT, int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteT = nbrSiteT;
  this->NbrSiteYZT = this->NbrSiteZ * this->NbrSiteY * this->NbrSiteT;
  this->NbrSiteZT = this->NbrSiteZ * this->NbrSiteT;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->KzMomentum = kzMomentum;
  this->KtMomentum = ktMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ * this->NbrSiteT;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ -1, this->NbrSiteT - 1, 0, 0, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, this->NbrSiteT - 1, 0, 0, 0, 0, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
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

// basic constructor when Sz is preserved
// 
// nbrFermions = number of fermions
// nbrSpinUp = number of particles with spin up
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// nbrSiteT = number of sites in the t direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// ktMomentum = momentum along the t direction
// memory = amount of memory granted for precalculations

FermionOnHyperCubicLatticeWithSpinMomentumSpace::FermionOnHyperCubicLatticeWithSpinMomentumSpace (int nbrFermions, int nbrSpinUp, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int nbrSiteT, int kxMomentum, int kyMomentum, int kzMomentum, int ktMomentum, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = nbrSpinUp;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteT = nbrSiteT;
  this->NbrSiteYZT = this->NbrSiteZ * this->NbrSiteY * this->NbrSiteT;
  this->NbrSiteZT = this->NbrSiteZ * this->NbrSiteT;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->KzMomentum = kzMomentum;
  this->KtMomentum = ktMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, this->NbrSiteT - 1, 0, 0, 0, 0, this->NbrFermionsUp);
  cout << "dim1 " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, this->NbrSiteT - 1, 0, 0, 0, 0, this->NbrFermionsUp, 0l);
      cout << "dim2 " << this->LargeHilbertSpaceDimension << endl;
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

FermionOnHyperCubicLatticeWithSpinMomentumSpace::FermionOnHyperCubicLatticeWithSpinMomentumSpace(const FermionOnHyperCubicLatticeWithSpinMomentumSpace& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->NbrSiteZ = fermions.NbrSiteZ;
  this->NbrSiteT = fermions.NbrSiteT;
  this->NbrSiteYZT = fermions.NbrSiteYZT;
  this->NbrSiteZT = fermions.NbrSiteZT;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->KzMomentum = fermions.KzMomentum;
  this->KtMomentum = fermions.KtMomentum;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->HighestBit = fermions.HighestBit;
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

FermionOnHyperCubicLatticeWithSpinMomentumSpace::~FermionOnHyperCubicLatticeWithSpinMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnHyperCubicLatticeWithSpinMomentumSpace& FermionOnHyperCubicLatticeWithSpinMomentumSpace::operator = (const FermionOnHyperCubicLatticeWithSpinMomentumSpace& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
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
  this->NbrSiteZ = fermions.NbrSiteZ;
  this->NbrSiteYZT = fermions.NbrSiteYZT;
  this->NbrSiteZT = fermions.NbrSiteZT;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->KzMomentum = fermions.KzMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->SzFlag = fermions.SzFlag;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnHyperCubicLatticeWithSpinMomentumSpace::Clone()
{
  return new FermionOnHyperCubicLatticeWithSpinMomentumSpace(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnHyperCubicLatticeWithSpinMomentumSpace::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> (i << 1));
      int TmpKx = i / this->NbrSiteYZT;
      int TmpKy = i % this->NbrSiteYZT;
      int TmpKz = TmpKy % this->NbrSiteZT;
      TmpKy /= this->NbrSiteZT;
      int TmpKt = TmpKz % this->NbrSiteT;
      TmpKz /= this->NbrSiteT;
      if ((Tmp & 0x1l) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << "," << TmpKt << ",+)";
      if ((Tmp & 0x2l) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << "," << TmpKt << ",-)";
    }
  Str << "]";
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along t for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnHyperCubicLatticeWithSpinMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, long pos)
{
  if (currentKt < 0)
    {
      currentKt = this->NbrSiteT - 1;
      currentKz--;
      if (currentKz < 0)
	{
	  currentKz = this->NbrSiteZ - 1;
	  currentKy--;
	  if (currentKy < 0)
	    {
	      currentKy = this->NbrSiteY - 1;
	      currentKx--;
	    }
	}
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum)
	   && ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum)&& ((currentTotalKt % this->NbrSiteT) == this->KtMomentum))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      for (int l = currentKt; l >= 0; --l)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((currentKy + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && 
	      (((currentKz + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum) && (((l + currentTotalKt) % this->NbrSiteT) == this->KtMomentum))
	    {
	      this->StateDescription[pos] = 0x2ul << (((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) * this->NbrSiteT + l) << 1);
	      ++pos;
	      this->StateDescription[pos] = 0x1ul << (((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) * this->NbrSiteT + l) << 1);
	      ++pos;
	    }
	}
      for (int k = currentKz; k >= 0; --k)
	{
	  for (int l = this->NbrSiteT - 1; l >= 0; --l)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((currentKy + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && 
		  (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum) && (((l + currentTotalKt) % this->NbrSiteT) == this->KtMomentum))
		{
		  this->StateDescription[pos] = 0x2ul << (((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + k) * this->NbrSiteT + l) << 1);
		  ++pos;
		  this->StateDescription[pos] = 0x1ul << (((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + k) * this->NbrSiteT + l) << 1);
		  ++pos;
		}
	    }
	}
      for (int j = currentKy - 1; j >= 0; --j)
	{
	  for (int k = this->NbrSiteZ - 1; k >= 0; --k)
	    {
	      for (int l = this->NbrSiteT - 1; l >= 0; --l)
		{
		  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum)
		      && (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum) && (((l + currentTotalKt) % this->NbrSiteT) == this->KtMomentum))
		    {
		      this->StateDescription[pos] = 0x2ul << (((((currentKx * this->NbrSiteY) + j) * this->NbrSiteZ + k) * this->NbrSiteT + l) << 1);
		      ++pos;
		      this->StateDescription[pos] = 0x1ul << (((((currentKx * this->NbrSiteY) + j) * this->NbrSiteZ + k) * this->NbrSiteT + l) << 1);
		      ++pos;
		    }
		}
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      for (int k = this->NbrSiteZ - 1; k >= 0; --k)
		{
		  for (int l = this->NbrSiteT - 1; l >= 0; --l)
		    {
		      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum)
			  && (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum) && (((l + currentTotalKt) % this->NbrSiteT) == this->KtMomentum))
			{
			  this->StateDescription[pos] = 0x2ul << ((((i * this->NbrSiteY) + j) * this->NbrSiteZ + k) << 1);
			  ++pos;
			  this->StateDescription[pos] = 0x1ul << ((((i * this->NbrSiteY) + j) * this->NbrSiteZ + k) << 1);
			  ++pos;
			}
		    }
		}
	    }
	}
      return pos;
    }
  long TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz), currentTotalKt + (2 * currentKt), pos);
  unsigned long Mask = 0x3ul << (((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) * this->NbrSiteT + currentKt) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, currentTotalKt + currentKt, pos);
  Mask = 0x2ul << (((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) * this->NbrSiteT + currentKt) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, currentTotalKt + currentKt, pos);
  Mask = 0x1ul << (((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) * this->NbrSiteT + currentKt) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx, currentTotalKy, currentTotalKz, currentTotalKt, pos);
};

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along t for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// nbrSpinUp = number of fermions with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnHyperCubicLatticeWithSpinMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, int nbrSpinUp, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }

  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrFermions))
    return 0l;

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
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      if (nbrSpinUp == 1)
	{
	  for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 1);
		  ++pos;
		}
	      for (int i = currentKx - 1; i >= 0; --i)
		{
		  for (int j = this->NbrSiteY - 1; j >= 0; --j)
		    {
		      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
			{
			  this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) << 1);
			  ++pos;
			}
		    }
		}
	    }
	}
      else
	{
	  for (int j = currentKy; j >= 0; --j)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 1);
		  ++pos;
		}
	      for (int i = currentKx - 1; i >= 0; --i)
		{
		  for (int j = this->NbrSiteY - 1; j >= 0; --j)
		    {
		      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
			{
			  this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) << 1);
			  ++pos;
			}
		    }
		}
	    }
	}
      return pos;
    }
  long TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz), currentTotalKt + (2 * currentKt), nbrSpinUp - 1, pos);
  unsigned long Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, currentTotalKt + currentKt, nbrSpinUp - 1, pos);
  Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, currentTotalKt + currentKt, nbrSpinUp, pos);
  Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx, currentTotalKy, currentTotalKz, currentTotalKt, nbrSpinUp, pos);  
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along t for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// return value = Hilbert space dimension

long FermionOnHyperCubicLatticeWithSpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt)
{
  if (currentKt < 0)
    {
      currentKt = this->NbrSiteT - 1;
      currentKz--;
      if (currentKz < 0)
	{
	  currentKz = this->NbrSiteZ - 1;
	  currentKy--;
	  if (currentKy < 0)
	    {
	      currentKy = this->NbrSiteY - 1;
	      currentKx--;
	    }
	}
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum)
	  && ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum) && ((currentTotalKt % this->NbrSiteT) == this->KtMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int l = currentKt; l >= 0; --l)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((currentKy + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && 
	      (((currentKz + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum) && (((l + currentTotalKt) % this->NbrSiteT) == this->KtMomentum))
	    Count += 2l;
	}
      for (int k = currentKz; k >= 0; --k)
	{
	  for (int l = this->NbrSiteT - 1; l >= 0; --l)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((currentKy + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && 
		  (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum) && (((l + currentTotalKt) % this->NbrSiteT) == this->KtMomentum))
		Count += 2l;
	    }
	}
      for (int j = currentKy - 1; j >= 0; --j)
	{
	  for (int k = this->NbrSiteZ - 1; k >= 0; --k)
	    {
	      for (int l = this->NbrSiteT - 1; l >= 0; --l)
		{
		  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum)
		      && (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum) && (((l + currentTotalKt) % this->NbrSiteT) == this->KtMomentum))
		    Count += 2l;
		}
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      for (int k = this->NbrSiteZ - 1; k >= 0; --k)
		{
		  for (int l = this->NbrSiteT - 1; l >= 0; --l)
		    {
		      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum)
			  && (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum) && (((l + currentTotalKt) % this->NbrSiteT) == this->KtMomentum))
			Count += 2l;
		    }
		}
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz), currentTotalKt + (2 * currentKt));
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, currentTotalKt + currentKt));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx, currentTotalKy, currentTotalKz, currentTotalKt);
  return Count;
}


// evaluate Hilbert space dimension with a fixed number of fermions with spin up
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentKt = current momentum along t for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentTotalKt = current total momentum along t
// nbrSpinUp = number of fermions with spin up
// return value = Hilbert space dimension

long FermionOnHyperCubicLatticeWithSpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentKt, int currentTotalKx, int currentTotalKy, int currentTotalKz, int currentTotalKt, int nbrSpinUp)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrFermions))
    return 0l;

  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Count++;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count++;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz), currentTotalKt + (2 * currentKt), nbrSpinUp - 1);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, currentTotalKt + currentKt, nbrSpinUp);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, currentTotalKt + currentKt, nbrSpinUp - 1);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy, currentKz, currentKt - 1, currentTotalKx, currentTotalKy, currentTotalKz, currentTotalKt, nbrSpinUp);
  return Count;
}

