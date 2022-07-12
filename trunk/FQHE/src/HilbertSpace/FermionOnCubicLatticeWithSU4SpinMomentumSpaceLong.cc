////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                class of fermions on cubic lattice with SU(4) spin          //
//               in momentum space that allows a number of sites up to        //
//                  31 (for systems with 128 bit integer support)             //
//               or 15 (on 32 bit systems without 128 bit integer support)    //                        
//                                                                            //
//                        last modification : 24/09/2011                      //
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
#include "HilbertSpace/FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong.h"
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
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// memory = amount of memory granted for precalculations

FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong::FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int kxMomentum, int kyMomentum, int kzMomentum, unsigned long memory)
{ 
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->PzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUpPlus = 0;
  this->NbrFermionsDownPlus = 0;
  this->NbrFermionsUpMinus = 0;
  this->NbrFermionsDownMinus = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteYZ = this->NbrSiteZ * this->NbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->KzMomentum = kzMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ -1, 0, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, 0, 0, 0, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      this->GenerateLookUpTable(memory);
      double coef = 0.0;
/*
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
      {
	for (int k1 = 0; k1 < 2; k1++)
	  {
	    for (int k2 = 0; k2 < 2; k2++)
	    {
	      double coef1 = this->AsigmaAsigma(i,k1,k2, 0, 0);
	      if (coef1 != 0)
	      {
		if (this->AdsigmaAdsigma(k1,k2,0,0,coef) != this->HilbertSpaceDimension)
		{
		  cout << i ;
		  this->PrintState(cout, i) << " : " ;
		  cout << this->AdsigmaAdsigma(k1,k2,0,0,coef) ;
		  this->PrintState(cout,this->AdsigmaAdsigma(k1,k2,0,0,coef)) << " " << coef1 << " " << coef << endl;
		}
	    }
	  }
	  }
      }*/
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(ULONGLONG) + sizeof(int));
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

FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong::FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong(const FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong& fermions)
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
  this->NbrSiteYZ = fermions.NbrSiteYZ;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->KzMomentum = fermions.KzMomentum;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->PzFlag = fermions.PzFlag;
  this->NbrFermionsUpPlus = fermions.NbrFermionsUpPlus;
  this->NbrFermionsDownPlus = fermions.NbrFermionsDownPlus;
  this->NbrFermionsUpMinus = fermions.NbrFermionsUpMinus;
  this->NbrFermionsDownMinus = fermions.NbrFermionsDownMinus;
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
}

// destructor
//

FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong::~FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong& FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong::operator = (const FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
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
  this->NbrSiteYZ = fermions.NbrSiteYZ;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->KzMomentum = fermions.KzMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->SzFlag = fermions.SzFlag;
  this->PzFlag = fermions.PzFlag;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUpPlus = fermions.NbrFermionsUpPlus;
  this->NbrFermionsDownPlus = fermions.NbrFermionsDownPlus;
  this->NbrFermionsUpMinus = fermions.NbrFermionsUpMinus;
  this->NbrFermionsDownMinus = fermions.NbrFermionsDownMinus;
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

AbstractHilbertSpace* FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong::Clone()
{
  return new FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong::PrintState (ostream& Str, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  ULONGLONG Tmp;
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> (i << 2));
      int TmpKx = i / this->NbrSiteYZ;
      int TmpKy = i % this->NbrSiteYZ;
      int TmpKz = TmpKy % this->NbrSiteZ;
      TmpKy /= this->NbrSiteZ;
      if ((Tmp & ((ULONGLONG) 0x8ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << ",u,A)";
      if ((Tmp & ((ULONGLONG) 0x4ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << ",u,B)";
      if ((Tmp & ((ULONGLONG) 0x2ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << ",d,A)";
      if ((Tmp & ((ULONGLONG) 0x1ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << ",d,B)";
    }
  Str << "]";
//      Str  << this->StateHighestBit[state] << " " << hex << ((unsigned long) (TmpState >> 64)) << "." << ((unsigned long) TmpState) << dec << endl;

//   Str << " " << state << " " << this->FindStateIndex (TmpState, this->StateHighestBit[state]);
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz, long pos)
{
  //  cout << nbrFermions << " " << currentKx << " " << currentKy << " " << currentKz << " " << currentTotalKx << " " << currentTotalKy << " " << currentTotalKz << endl;
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
  if (nbrFermions < 0)
    return pos;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum)
	   && ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum))
	{
	  this->StateDescription[pos] = ((ULONGLONG) 0x0ul);	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      for (int k = currentKz; k >= 0; --k)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((currentKy + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && 
	      (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum))
	    {
 	      this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + k) << 2);
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + k) << 2);
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + k) << 2);
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + k) << 2);
	      ++pos;
	    }
	}
      for (int j = currentKy - 1; j >= 0; --j)
	{
	  for (int k = this->NbrSiteZ - 1; k >= 0; --k)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum)
		  && (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum))
		{
		  this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << ((((currentKx * this->NbrSiteY) + j) * this->NbrSiteZ + k) << 2);
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << ((((currentKx * this->NbrSiteY) + j) * this->NbrSiteZ + k) << 2);
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << ((((currentKx * this->NbrSiteY) + j) * this->NbrSiteZ + k) << 2);
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << ((((currentKx * this->NbrSiteY) + j) * this->NbrSiteZ + k) << 2);
		  ++pos;
		}
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      for (int k = this->NbrSiteZ - 1; k >= 0; --k)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum)
		      && (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum))
		    {
		      this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << ((((i * this->NbrSiteY) + j) * this->NbrSiteZ + k) << 2);
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << ((((i * this->NbrSiteY) + j) * this->NbrSiteZ + k) << 2);
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << ((((i * this->NbrSiteY) + j) * this->NbrSiteZ + k) << 2);
		      ++pos;
		      this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << ((((i * this->NbrSiteY) + j) * this->NbrSiteZ + k) << 2);
		      ++pos;
		    }
		}
	    }
	}
      return pos;
    }

  long TmpPos = this->GenerateStates(nbrFermions - 4, currentKx, currentKy, currentKz - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), currentTotalKz + (4 * currentKz), pos);
  ULONGLONG Mask = ((ULONGLONG) 0xful) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy, currentKz - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), currentTotalKz + (3 * currentKz), pos);
  Mask = ((ULONGLONG) 0xeul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy, currentKz - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), currentTotalKz + (3 * currentKz), pos);
  Mask = ((ULONGLONG) 0xdul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy, currentKz - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz), pos);
  Mask = ((ULONGLONG) 0xcul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy, currentKz - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), currentTotalKz + (3 * currentKz), pos);
  Mask = ((ULONGLONG) 0xbul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy, currentKz - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz), pos);
  Mask = ((ULONGLONG) 0xaul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy, currentKz - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz), pos);
  Mask = ((ULONGLONG) 0x9ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy, currentKz - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, pos);
  Mask = ((ULONGLONG) 0x8ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy, currentKz - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), currentTotalKz + (3 * currentKz), pos);
  Mask = ((ULONGLONG) 0x7ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy, currentKz - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz), pos);
  Mask = ((ULONGLONG) 0x6ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy, currentKz - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz), pos);
  Mask = ((ULONGLONG) 0x5ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy, currentKz - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, pos);
  Mask = ((ULONGLONG) 0x4ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy, currentKz - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz), pos);
  Mask = ((ULONGLONG) 0x3ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy, currentKz - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, pos);
  Mask = ((ULONGLONG) 0x2ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy, currentKz - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz, pos);
  Mask = ((ULONGLONG) 0x1ul) << ((((currentKx * this->NbrSiteY) + currentKy) * this->NbrSiteZ + currentKz) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;


  return this->GenerateStates(nbrFermions, currentKx, currentKy, currentKz - 1, currentTotalKx, currentTotalKy, currentTotalKz, pos);
};


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// return value = Hilbert space dimension

long FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentKz, int currentTotalKx, int currentTotalKy, int currentTotalKz)
{
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
  if (nbrFermions < 0)
    return 0l;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum)
	  && ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int k = currentKz; k >= 0; --k)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((currentKy + currentTotalKy) % this->NbrSiteY) == this->KyMomentum) && 
	      (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum))
	    Count += 4l;
	}
      for (int j = currentKy - 1; j >= 0; --j)
	{
	  for (int k = this->NbrSiteZ - 1; k >= 0; --k)
	    {
	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum)
		  && (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum))
		Count += 4l;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      for (int k = this->NbrSiteZ - 1; k >= 0; --k)
		{
		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum)
		      && (((k + currentTotalKz) % this->NbrSiteZ) == this->KzMomentum))
		    Count += 4l;
		}
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy, currentKz - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), currentTotalKz + (4 * currentKz));
  Count += (4 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy, currentKz - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), currentTotalKz + (3 * currentKz)));
  Count += (6 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy, currentKz - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), currentTotalKz + (2 * currentKz)));
  Count += (4 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy, currentKz - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, currentTotalKz + currentKz));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy, currentKz - 1, currentTotalKx, currentTotalKy, currentTotalKz);
  return Count;
}


