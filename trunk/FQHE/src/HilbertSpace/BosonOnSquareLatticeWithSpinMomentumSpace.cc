////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons on square lattice with spin              //
//                                in momentum space                           //
//                                                                            //
//                        last modification : 02/12/2011                      //
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
#include "HilbertSpace/BosonOnSquareLatticeWithSpinMomentumSpace.h"
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
// nbrBosons = number of bosons
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

BosonOnSquareLatticeWithSpinMomentumSpace::BosonOnSquareLatticeWithSpinMomentumSpace (int nbrBosons, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrBosonsUp = 0;
  this->NbrBosonsDown = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->TemporaryState = new unsigned [this->NbrLzValue];
      this->ProdATemporaryState = new unsigned [this->NbrLzValue];
      this->Flag.Initialize();
      this->StateDescription = new unsigned* [this->HilbertSpaceDimension];
      this->StateLzMaxUp = new unsigned [this->HilbertSpaceDimension];
      this->StateLzMaxDown = new unsigned [this->HilbertSpaceDimension];
      int TmpLzMax = this->LzMax;
      if (this->ShiftedTotalLz < TmpLzMax)
	{
	  TmpLzMax = this->ShiftedTotalLz;	  
	}

      long TmpDim = GenerateStates(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l);
      
      if (TmpDim!=this->HilbertSpaceDimension)
	cout << "Count inconsistent: "<<TmpDim<<" vs " << this->HilbertSpaceDimension<<endl;

      this->GenerateLookUpTable(memory);
      this->CoherenceFactors=new double[NbrBosons * NbrBosons + 1];
      for (int i=0; i< (NbrBosons * NbrBosons + 1); ++i)
	this->CoherenceFactors[i] = sqrt((double)i);
      this->KeptCoordinates = new int;
      (*(this->KeptCoordinates)) = -1;
      this->Minors = 0;
      
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
// nbrBosons = number of bosons
// nbrSpinUp = number of particles with spin up
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

BosonOnSquareLatticeWithSpinMomentumSpace::BosonOnSquareLatticeWithSpinMomentumSpace (int nbrBosons, int nbrSpinUp, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->NbrBosonsUp = nbrSpinUp;
  this->NbrBosonsDown = this->NbrBosons - this->NbrBosons;
  this->TotalSpin = this->NbrBosonsUp - this->NbrBosonsDown;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->NbrBosonsUp);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->TemporaryState = new unsigned [this->NbrLzValue];
      this->ProdATemporaryState = new unsigned [this->NbrLzValue];
      this->Flag.Initialize();
      this->StateDescription = new unsigned* [this->HilbertSpaceDimension];
      this->StateLzMaxUp = new unsigned [this->HilbertSpaceDimension];
      this->StateLzMaxDown = new unsigned [this->HilbertSpaceDimension];
      int TmpLzMax = this->LzMax;
      if (this->ShiftedTotalLz < TmpLzMax)
	{
	  TmpLzMax = this->ShiftedTotalLz;	  
	}

      long TmpDim = GenerateStates(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->NbrBosonsUp, 0l);
      
      if (TmpDim!=this->HilbertSpaceDimension)
	cout << "Count inconsistent: "<<TmpDim<<" vs " << this->HilbertSpaceDimension<<endl;

      this->GenerateLookUpTable(memory);
      this->CoherenceFactors=new double[NbrBosons * NbrBosons + 1];
      for (int i=0; i< (NbrBosons * NbrBosons + 1); ++i)
	this->CoherenceFactors[i] = sqrt((double)i);
      this->KeptCoordinates = new int;
      (*(this->KeptCoordinates)) = -1;
      this->Minors = 0;
      
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


  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->SzFlag = true;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrBosonsUp = nbrSpinUp;
  this->NbrBosonsDown = this->NbrBosons - this->NbrBosonsUp;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->NbrBosonsUp);
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
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->NbrBosonsUp, 0l);
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

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSquareLatticeWithSpinMomentumSpace::BosonOnSquareLatticeWithSpinMomentumSpace(const BosonOnSquareLatticeWithSpinMomentumSpace& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->TotalSpin = bosons.TotalSpin;
  this->SzFlag = bosons.SzFlag;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->HighestBit = bosons.HighestBit;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  this->SignLookUpTable = bosons.SignLookUpTable;
  this->SignLookUpTableMask = bosons.SignLookUpTableMask;
  this->MaximumSignLookUp = bosons.MaximumSignLookUp;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

BosonOnSquareLatticeWithSpinMomentumSpace::~BosonOnSquareLatticeWithSpinMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSquareLatticeWithSpinMomentumSpace& BosonOnSquareLatticeWithSpinMomentumSpace::operator = (const BosonOnSquareLatticeWithSpinMomentumSpace& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->NbrLzValue = bosons.NbrLzValue;
  this->SzFlag = bosons.SzFlag;
  this->TotalSpin = bosons.TotalSpin;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSquareLatticeWithSpinMomentumSpace::Clone()
{
  return new BosonOnSquareLatticeWithSpinMomentumSpace(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSquareLatticeWithSpinMomentumSpace::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> (i << 1));
      int TmpKx = i / this->NbrSiteY;
      int TmpKy = i % this->NbrSiteY;
      if ((Tmp & 0x2l) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << ",+)";
      if ((Tmp & 0x1l) != 0ul)
	Str << "(" << TmpKx << "," << TmpKy << ",-)";
    }
  Str << "]";
//   Str << " " << TmpState; 
//   Str << " " << hex << TmpState << dec; 
//   Str << " " << state;
//   int TmpLzMax = (this->LzMax << 1) + 1;
//   while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
//     --TmpLzMax;
//   Str << " " << this->FindStateIndex(TmpState, TmpLzMax);
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSquareLatticeWithSpinMomentumSpace::GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos)
{
//   if (currentKy < 0)
//     {
//       currentKy = this->NbrSiteY - 1;
//       currentKx--;
//     }
//   if (nbrBosons == 0)
//     {
//       if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
// 	{
// 	  this->StateDescription[pos] = 0x0ul;	  
// 	  return (pos + 1l);
// 	}
//       else	
// 	return pos;
//     }
//   if (currentKx < 0)
//     return pos;
//   if (nbrBosons == 1)
//     {
//       for (int j = currentKy; j >= 0; --j)
// 	{
// 	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
// 	    {
// 	      this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 1);
// 	      ++pos;
// 	      this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 1);
// 	      ++pos;
// 	    }
// 	}
//       for (int i = currentKx - 1; i >= 0; --i)
// 	{
// 	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
// 	    {
// 	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
// 		{
// 		  this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) << 1);
// 		  ++pos;
//  		  this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) << 1);
//  		  ++pos;
// 		}
// 	    }
// 	}
//       return pos;
//     }
//   long TmpPos = this->GenerateStates(nbrBosons - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), pos);
//   unsigned long Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
//   for (; pos < TmpPos; ++pos)
//     this->StateDescription[pos] |= Mask;
//   TmpPos = this->GenerateStates(nbrBosons - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
//   Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
//   for (; pos < TmpPos; ++pos)
//     this->StateDescription[pos] |= Mask;
//    TmpPos = this->GenerateStates(nbrBosons - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, pos);
//    Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
//    for (; pos < TmpPos; ++pos)
//      this->StateDescription[pos] |= Mask;
   return this->GenerateStates(nbrBosons, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, pos);
};

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrSpinUp = number of bosons with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSquareLatticeWithSpinMomentumSpace::GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrSpinUp, long pos)
{
//   if (currentKy < 0)
//     {
//       currentKy = this->NbrSiteY - 1;
//       currentKx--;
//     }

//   if ((nbrSpinUp < 0) || (nbrSpinUp > nbrBosons))
//     return 0l;

//   if (nbrBosons == 0)
//     {
//       if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
// 	{
// 	  this->StateDescription[pos] = 0x0ul;	  
// 	  return (pos + 1l);
// 	}
//       else	
// 	return pos;
//     }
//   if (currentKx < 0)
//     return pos;
//   if (nbrBosons == 1)
//     {
//       if (nbrSpinUp == 1)
// 	{
// 	  for (int j = currentKy; j >= 0; --j)
// 	    {
// 	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
// 		{
// 		  this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 1);
// 		  ++pos;
// 		}
// 	    }
// 	  for (int i = currentKx - 1; i >= 0; --i)
// 	    {
// 	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
// 		{
// 		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
// 		    {
// 		      this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) << 1);
// 		      ++pos;
// 			}
// 		}
// 	    }
// 	}
//       else
// 	{
// 	  for (int j = currentKy; j >= 0; --j)
// 	    {
// 	      if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
// 		{
// 		  this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 1);
// 		  ++pos;
// 		}
// 	    }
// 	  for (int i = currentKx - 1; i >= 0; --i)
// 	    {
// 	      for (int j = this->NbrSiteY - 1; j >= 0; --j)
// 		{
// 		  if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
// 		    {
// 		      this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) << 1);
// 		      ++pos;
// 		    }
// 		}
// 	    }
// 	}
//       return pos;
//     }
//   long TmpPos = this->GenerateStates(nbrBosons - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy), nbrSpinUp - 1, pos);
//   unsigned long Mask = 0x3ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
//   for (; pos < TmpPos; ++pos)
//     this->StateDescription[pos] |= Mask;
//   TmpPos = this->GenerateStates(nbrBosons - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrSpinUp - 1, pos);
//   Mask = 0x2ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
//   for (; pos < TmpPos; ++pos)
//     this->StateDescription[pos] |= Mask;
//   TmpPos = this->GenerateStates(nbrBosons - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy, nbrSpinUp, pos);
//   Mask = 0x1ul << (((currentKx * this->NbrSiteY) + currentKy) << 1);
//   for (; pos < TmpPos; ++pos)
//     this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrBosons, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, nbrSpinUp, pos);  
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonOnSquareLatticeWithSpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrBosons == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Count += 2l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 2l;
	    }
	}
      return Count;
    }
  for (int i = nbrBosons; i >= 0; --i)
    Count += ((long) i + 1l) * this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKx, currentKy - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy));
  return Count;
}


// evaluate Hilbert space dimension with a fixed number of bosons with spin up
//
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrSpinUp = number of bosons with spin up
// return value = Hilbert space dimension

long BosonOnSquareLatticeWithSpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrSpinUp)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrBosons))
    return 0l;

  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrBosons == 1)
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
  for (int i = nbrBosons; i >= 0; --i)
    for (int j = i; j >= 0; --j)
      Count += this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKx, currentKy - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy), nbrSpinUp - j);
  return Count;
}


