////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(4) spin         //
//               in momentum space that allows a number of sites up to        //
//                  31 (for systems with 128 bit integer support)             //
//               or 15 (on 32 bit systems without 128 bit integer support)    //
//                                                                            //
//                        last modification : 26/09/2011                      //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"
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
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong::FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
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
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	this->PrintState(cout, i) << " " << hex << this->StateDescription[i] << dec << endl;
      this->GenerateLookUpTable(memory);
      
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


// constructor when conserving spin and isospin
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong::FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, int totalSpin, int totalIsospin, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->PzFlag = true;
  this->TotalLz = 0;
  this->TotalSpin = totalSpin;
  this->TotalIsospin = totalIsospin;
  this->NbrFermionsUpPlus = ((this->NbrFermions+this->TotalSpin)/2 + this->TotalIsospin)/2;
  this->NbrFermionsDownPlus = ((this->NbrFermions-this->TotalSpin)/2 + this->TotalIsospin)/2;
  this->NbrFermionsUpMinus = ((this->NbrFermions+this->TotalSpin)/2 - this->TotalIsospin)/2;
  this->NbrFermionsDownMinus = ((this->NbrFermions-this->TotalSpin)/2 - this->TotalIsospin)/2;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, (this->NbrFermions+this->TotalSpin)/2, (this->NbrFermions+this->TotalIsospin)/2);
  
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l, (this->NbrFermions+this->TotalSpin)/2, (this->NbrFermions+this->TotalIsospin)/2);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	this->PrintState(cout, i) << " " << hex << this->StateDescription[i] << dec << endl;
      this->GenerateLookUpTable(memory);
      
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

FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong::FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong(const FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
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

FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong::~FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong& FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong::operator = (const FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong& fermions)
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
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
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

AbstractHilbertSpace* FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong::Clone()
{
  return new FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong::PrintState (ostream& Str, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  ULONGLONG Tmp;
  cout << "lzMax = " << (this->NbrLzValue) << " : " ;
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> (i << 2));
      int TmpKx = i / this->NbrSiteY;
      int TmpKy = i % this->NbrSiteY;
      if ((Tmp & ((ULONGLONG) 0x8ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",u,A)";
      if ((Tmp & ((ULONGLONG) 0x4ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",u,B)";
      if ((Tmp & ((ULONGLONG) 0x2ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",d,A)";
      if ((Tmp & ((ULONGLONG) 0x1ul)) != ((ULONGLONG) 0x0ul))
	Str << "(" << TmpKx << "," << TmpKy << ",d,B)";
    }
  Str << "]";
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

long FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions < 0)
    return pos;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
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
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << (((currentKx * this->NbrSiteY) + j) << 2);
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << (((currentKx * this->NbrSiteY) + j) << 2);
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << (((currentKx * this->NbrSiteY) + j) << 2);
	      ++pos;
	      this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << (((currentKx * this->NbrSiteY) + j) << 2);
	      ++pos;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << (((i * this->NbrSiteY) + j) << 2);
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << (((i * this->NbrSiteY) + j) << 2);
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << (((i * this->NbrSiteY) + j) << 2);
		  ++pos;
		  this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << (((i * this->NbrSiteY) + j) << 2);
		  ++pos;
		}
	    }
	}
      return pos;
    }

  long TmpPos = this->GenerateStates(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), pos);
  ULONGLONG Mask = ((ULONGLONG) 0xful) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos);
  Mask = ((ULONGLONG) 0xeul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos);
  Mask = ((ULONGLONG) 0xdul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = ((ULONGLONG) 0xcul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos);
  Mask = ((ULONGLONG) 0xbul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = ((ULONGLONG) 0xaul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = ((ULONGLONG) 0x9ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
  Mask = ((ULONGLONG) 0x8ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos);
  Mask = ((ULONGLONG) 0x7ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = ((ULONGLONG) 0x6ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = ((ULONGLONG) 0x5ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
  Mask = ((ULONGLONG) 0x4ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
  Mask = ((ULONGLONG) 0x3ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
  Mask = ((ULONGLONG) 0x2ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
  Mask = ((ULONGLONG) 0x1ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;


  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos);
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

long FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos, int nbrFermionsUp, int nbrFermionsPlus)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0) || (nbrFermionsUp < 0) || (nbrFermionsPlus < 0) ||  
      (nbrFermionsUp > nbrFermions) || (nbrFermionsPlus > nbrFermions))
    return pos;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
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
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
            switch (nbrFermionsUp + 2*nbrFermionsPlus)
            {
                case 3:
                    this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << (((currentKx * this->NbrSiteY) + j) << 2);
                    break;
                case 1:
                    this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << (((currentKx * this->NbrSiteY) + j) << 2);
                    break;
                case 2:
                    this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << (((currentKx * this->NbrSiteY) + j) << 2);
                    break;
                case 0:
                    this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << (((currentKx * this->NbrSiteY) + j) << 2);
                    break;
            }
	      ++pos;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
            {
                switch (nbrFermionsUp + 2*nbrFermionsPlus)
                {
                    case 3:
                        this->StateDescription[pos] = ((ULONGLONG) 0x8ul) << (((i * this->NbrSiteY) + j) << 2);
                        break;
                    case 1:
                        this->StateDescription[pos] = ((ULONGLONG) 0x4ul) << (((i * this->NbrSiteY) + j) << 2);
                        break;
                    case 2:
                        this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << (((i * this->NbrSiteY) + j) << 2);
                        break;
                    case 0:
                        this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << (((i * this->NbrSiteY) + j) << 2);
                        break;
                }
                ++pos;
            }
	    }
	}
      return pos;
    }

  long TmpPos = this->GenerateStates(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), pos, nbrFermionsUp - 2, nbrFermionsPlus - 2);
  ULONGLONG Mask = ((ULONGLONG) 0xful) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 2, nbrFermionsPlus - 2);
  Mask = ((ULONGLONG) 0xeul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 2, nbrFermionsPlus - 1);
  Mask = ((ULONGLONG) 0xdul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 2, nbrFermionsPlus - 1);
  Mask = ((ULONGLONG) 0xcul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 1, nbrFermionsPlus - 2);
  Mask = ((ULONGLONG) 0xbul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1, nbrFermionsPlus - 2);
  Mask = ((ULONGLONG) 0xaul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1, nbrFermionsPlus - 1);
  Mask = ((ULONGLONG) 0x9ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp - 1, nbrFermionsPlus - 1);
  Mask = ((ULONGLONG) 0x8ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), pos, nbrFermionsUp - 1, nbrFermionsPlus - 1);
  Mask = ((ULONGLONG) 0x7ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp - 1, nbrFermionsPlus - 1);
  Mask = ((ULONGLONG) 0x6ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp -  1, nbrFermionsPlus);
  Mask = ((ULONGLONG) 0x5ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp - 1, nbrFermionsPlus);
  Mask = ((ULONGLONG) 0x4ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos, nbrFermionsUp, nbrFermionsPlus - 1);
  Mask = ((ULONGLONG) 0x3ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp, nbrFermionsPlus - 1);
  Mask = ((ULONGLONG) 0x2ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;

  TmpPos = this->GenerateStates(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos, nbrFermionsUp, nbrFermionsPlus);
  Mask = ((ULONGLONG) 0x1ul) << (((currentKx * this->NbrSiteY) + currentKy) << 2);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;


  return this->GenerateStates(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos, nbrFermionsUp, nbrFermionsPlus);
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions < 0)
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
	    Count += 4l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 4l;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy));
  Count += (4 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  Count += (6 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy)));
  Count += (4 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}



// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrFermionsUp, int nbrFermionsPlus)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)|| (nbrFermionsUp < 0) || (nbrFermionsPlus < 0) ||  
      (nbrFermionsUp > nbrFermions) || (nbrFermionsPlus > nbrFermions))
    return 0l;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  return 1l;
	}
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
	    Count += 1l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 1l;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy), nbrFermionsUp - 2, nbrFermionsPlus - 2);
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 2, nbrFermionsPlus - 2);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 2, nbrFermionsPlus - 1);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 1, nbrFermionsPlus - 2);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy), nbrFermionsUp - 1, nbrFermionsPlus - 1);
  
  Count += (2 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 1, nbrFermionsPlus -1));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 2, nbrFermionsPlus - 1);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp -1, nbrFermionsPlus - 2);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp - 1, nbrFermionsPlus);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrFermionsUp, nbrFermionsPlus - 1);
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp, nbrFermionsPlus);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp - 1, nbrFermionsPlus);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp, nbrFermionsPlus - 1);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrFermionsUp - 1, nbrFermionsPlus - 1);
  
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, nbrFermionsUp, nbrFermionsPlus);
  return Count;
}
