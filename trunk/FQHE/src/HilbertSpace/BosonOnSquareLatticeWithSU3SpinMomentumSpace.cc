////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//               class of bosons on a square lattice with SU(3) spin          //
//                                in momentum space                           //
//                                                                            //
//                        last modification : 06/12/2011                      //
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
#include "HilbertSpace/BosonOnSquareLatticeWithSU3SpinMomentumSpace.h"
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
#include "GeneralTools/ArrayTools.h"

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

BosonOnSquareLatticeWithSU3SpinMomentumSpace::BosonOnSquareLatticeWithSU3SpinMomentumSpace (int nbrBosons, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->TotalTz = 0;
  this->TotalY = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryState1 = new unsigned long[this->NbrLzValue];
  this->TemporaryState2 = new unsigned long[this->NbrLzValue];
  this->TemporaryState3 = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryState1;
  this->TemporaryStateSigma[1] = this->TemporaryState2;
  this->TemporaryStateSigma[2] = this->TemporaryState3;
  this->ProdATemporaryState1 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState2 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState3 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryState1;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryState2;
  this->ProdATemporaryStateSigma[2] = this->ProdATemporaryState3;

  this->N1LzMax = this->LzMax + this->NbrBosons - 1;
  this->N2LzMax = this->LzMax + this->NbrBosons - 1;
  this->N3LzMax = this->LzMax + this->NbrBosons - 1;
  this->FermionicLzMax = this->N1LzMax;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescription1 = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescription2 = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescription3 = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionSigma[0] = this->StateDescription1;
      this->StateDescriptionSigma[1] = this->StateDescription2;
      this->StateDescriptionSigma[2] = this->StateDescription3;
      this->Flag.Initialize();
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      for (long i = 0; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->StateDescription1[i];
	  unsigned long Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescription1[i] >>= this->NbrBosons - Tmp; 
	  TmpState = this->StateDescription2[i];
	  Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescription2[i] >>= this->NbrBosons - Tmp; 
	  TmpState = this->StateDescription3[i];
	  Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescription3[i] >>= this->NbrBosons - Tmp; 
	}
      SortTripleElementArrayDownOrdering<unsigned long>(this->StateDescription1, this->StateDescription2, this->StateDescription3, TmpLargeHilbertSpaceDimension);
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (4 * sizeof(unsigned long));
      cout << "memory requested for Hilbert space = ";
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

BosonOnSquareLatticeWithSU3SpinMomentumSpace::BosonOnSquareLatticeWithSU3SpinMomentumSpace(const BosonOnSquareLatticeWithSU3SpinMomentumSpace& bosons)
{
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->N1LzMax = bosons.N1LzMax;
  this->N2LzMax = bosons.N2LzMax;
  this->N3LzMax = bosons.N3LzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryState1 = new unsigned long[this->NbrLzValue];
  this->TemporaryState2 = new unsigned long[this->NbrLzValue];
  this->TemporaryState3 = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryState1;
  this->TemporaryStateSigma[1] = this->TemporaryState2;
  this->TemporaryStateSigma[2] = this->TemporaryState3;
  this->ProdATemporaryState1 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState2 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState3 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryState1;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryState2;
  this->ProdATemporaryStateSigma[2] = this->ProdATemporaryState3;
  this->StateDescription1 = bosons.StateDescription1;
  this->StateDescription2 = bosons.StateDescription2;
  this->StateDescription3 = bosons.StateDescription3;
  this->StateDescriptionSigma[0] = this->StateDescription1;
  this->StateDescriptionSigma[1] = this->StateDescription2;
  this->StateDescriptionSigma[2] = this->StateDescription3;
  this->NbrUniqueStateDescription1 = bosons.NbrUniqueStateDescription1;
  this->UniqueStateDescription1 = bosons.UniqueStateDescription1;
  this->UniqueStateDescriptionSubArraySize1 = bosons.UniqueStateDescriptionSubArraySize1;
  this->NbrUniqueStateDescription2 = bosons.NbrUniqueStateDescription2;
  this->UniqueStateDescription2 = bosons.UniqueStateDescription2;
  this->UniqueStateDescriptionSubArraySize2 = bosons.UniqueStateDescriptionSubArraySize2;
  this->FirstIndexUniqueStateDescription2 = bosons.FirstIndexUniqueStateDescription2;
}

// destructor
//

BosonOnSquareLatticeWithSU3SpinMomentumSpace::~BosonOnSquareLatticeWithSU3SpinMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSquareLatticeWithSU3SpinMomentumSpace& BosonOnSquareLatticeWithSU3SpinMomentumSpace::operator = (const BosonOnSquareLatticeWithSU3SpinMomentumSpace& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription1;
      delete[] this->StateDescription2;
      delete[] this->StateDescription3;
      delete[] this->UniqueStateDescription1;
      delete[] this->UniqueStateDescriptionSubArraySize1;
      delete[] this->NbrUniqueStateDescription2;
      for (long i = 0l; i < this->NbrUniqueStateDescription1; ++i)
	{
	  delete[] this->UniqueStateDescription2[i];
	  delete[] this->UniqueStateDescriptionSubArraySize2[i];
	  delete[] this->FirstIndexUniqueStateDescription2[i];
	}
      delete[] this->UniqueStateDescription2;
      delete[] this->UniqueStateDescriptionSubArraySize2;
      delete[] this->FirstIndexUniqueStateDescription2;
    }
  delete[] this->TemporaryState1;
  delete[] this->TemporaryState2;
  delete[] this->TemporaryState3;
  delete[] this->ProdATemporaryState1;
  delete[] this->ProdATemporaryState2;
  delete[] this->ProdATemporaryState3;
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->N1LzMax = bosons.N1LzMax;
  this->N2LzMax = bosons.N2LzMax;
  this->N3LzMax = bosons.N3LzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryState1 = new unsigned long[this->NbrLzValue];
  this->TemporaryState2 = new unsigned long[this->NbrLzValue];
  this->TemporaryState3 = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryState1;
  this->TemporaryStateSigma[1] = this->TemporaryState2;
  this->TemporaryStateSigma[2] = this->TemporaryState3;
  this->ProdATemporaryState1 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState2 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryState3 = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryState1;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryState2;
  this->ProdATemporaryStateSigma[2] = this->ProdATemporaryState3;
  this->StateDescription1 = bosons.StateDescription1;
  this->StateDescription2 = bosons.StateDescription2;
  this->StateDescription3 = bosons.StateDescription3;
  this->StateDescriptionSigma[0] = this->StateDescription1;
  this->StateDescriptionSigma[1] = this->StateDescription2;
  this->StateDescriptionSigma[2] = this->StateDescription3;
  this->NbrUniqueStateDescription1 = bosons.NbrUniqueStateDescription1;
  this->UniqueStateDescription1 = bosons.UniqueStateDescription1;
  this->UniqueStateDescriptionSubArraySize1 = bosons.UniqueStateDescriptionSubArraySize1;
  this->NbrUniqueStateDescription2 = bosons.NbrUniqueStateDescription2;
  this->UniqueStateDescription2 = bosons.UniqueStateDescription2;
  this->UniqueStateDescriptionSubArraySize2 = bosons.UniqueStateDescriptionSubArraySize2;
  this->FirstIndexUniqueStateDescription2 = bosons.FirstIndexUniqueStateDescription2;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSquareLatticeWithSU3SpinMomentumSpace::Clone()
{
  return new BosonOnSquareLatticeWithSU3SpinMomentumSpace(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSquareLatticeWithSU3SpinMomentumSpace::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescription1[state], this->StateDescription2[state], this->StateDescription3[state],
		       TemporaryState1, TemporaryState2, TemporaryState3); 

  unsigned long Tmp;
  Str << "[";
  for (int i = 0; i <= this->LzMax; ++i)
    {
      int TmpKx = i / this->NbrSiteY;
      int TmpKy = i % this->NbrSiteY;
      if (this->TemporaryState1[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryState1[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << ",A)";
	}
      if (this->TemporaryState2[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryState2[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << ",B)";
	}
      if (this->TemporaryState3[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryState3[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << ",C)";
	}
    }
  Str << "]";
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentFermionicPosition1 = current fermionic position within the state description for the type 1 particles
// currentFermionicPosition2 = current fermionic position within the state description for the type 2 particles
// currentFermionicPosition3 = current fermionic position within the state description for the type 3 particles
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSquareLatticeWithSU3SpinMomentumSpace::GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentFermionicPosition1, int currentFermionicPosition2, int currentFermionicPosition3, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
 	  this->StateDescription1[pos] = 0x0ul;
 	  this->StateDescription2[pos] = 0x0ul;
 	  this->StateDescription3[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  for (int i = nbrBosons; i >= 0; --i)
    {
      unsigned long Mask1 = ((0x1ul << i) - 0x1ul) << (currentFermionicPosition1 - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  unsigned long Mask2 = ((0x1ul << j) - 0x1ul) << (currentFermionicPosition2 - j - 1);
	  for (int k = nbrBosons - i - j; k >= 0; --k)
	    {
	      long TmpPos = this->GenerateStates(nbrBosons - i - j - k, currentKx, currentKy - 1, currentTotalKx + ((i + j + k) * currentKx), currentTotalKy + ((i + j + k) * currentKy), currentFermionicPosition1 - i - 1, currentFermionicPosition2 - j - 1, currentFermionicPosition3 - k - 1, pos);
	      unsigned long Mask3 = ((0x1ul << k) - 0x1ul) << (currentFermionicPosition3 - k - 1);
	      for (; pos < TmpPos; ++pos)
		{
 		  this->StateDescription1[pos] |= Mask1;
 		  this->StateDescription2[pos] |= Mask2;
 		  this->StateDescription3[pos] |= Mask3;
		}
	    }
	}
    }
  return pos;
};


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonOnSquareLatticeWithSU3SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
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
	    Count += 3l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 3l;
	    }
	}
      return Count;
    }
  for (int i = nbrBosons; i >= 0; --i)
    Count += ((((long) i + 1l) * ((long) i + 2l)) / 2l) * this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKx, currentKy - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy));
  return Count;
}
