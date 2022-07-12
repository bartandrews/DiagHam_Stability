////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//               class of bosons on a square lattice with SU(4) spin          //
//                                in momentum space                           //
//                                                                            //
//                        last modification : 21/12/2011                      //
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
#include "HilbertSpace/BosonOnSquareLatticeWithSU4SpinMomentumSpace.h"
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

BosonOnSquareLatticeWithSU4SpinMomentumSpace::BosonOnSquareLatticeWithSU4SpinMomentumSpace (int nbrBosons, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->TotalIsospin = 0;
  this->TotalEntanglement = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUpPlus;
  this->TemporaryStateSigma[1] = this->TemporaryStateUpMinus;
  this->TemporaryStateSigma[2] = this->TemporaryStateDownPlus;
  this->TemporaryStateSigma[3] = this->TemporaryStateDownMinus;
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUpPlus;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateUpMinus;
  this->ProdATemporaryStateSigma[2] = this->ProdATemporaryStateDownPlus;
  this->ProdATemporaryStateSigma[3] = this->ProdATemporaryStateDownMinus;

  this->NUpPlusLzMax = this->LzMax + this->NbrBosons - 1;
  this->NUpMinusLzMax = this->LzMax + this->NbrBosons - 1;
  this->NDownPlusLzMax = this->LzMax + this->NbrBosons - 1;
  this->NDownMinusLzMax = this->LzMax + this->NbrBosons - 1;
  this->FermionicLzMax = this->NUpPlusLzMax;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUpPlus = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionUpMinus = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDownPlus = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDownMinus = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionSigma[0] = this->StateDescriptionUpPlus;
      this->StateDescriptionSigma[1] = this->StateDescriptionUpMinus;
      this->StateDescriptionSigma[2] = this->StateDescriptionDownPlus;
      this->StateDescriptionSigma[3] = this->StateDescriptionDownMinus;
      this->Flag.Initialize();
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space count: " << this->LargeHilbertSpaceDimension << " generate: " << TmpLargeHilbertSpaceDimension << endl;
	}
      for (long i = 0; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->StateDescriptionUpPlus[i];
	  unsigned long Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescriptionUpPlus[i] >>= this->NbrBosons - Tmp; 
	  TmpState = this->StateDescriptionUpMinus[i];
	  Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescriptionUpMinus[i] >>= this->NbrBosons - Tmp; 
	  TmpState = this->StateDescriptionDownPlus[i];
	  Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescriptionDownPlus[i] >>= this->NbrBosons - Tmp; 
	  TmpState = this->StateDescriptionDownMinus[i];
	  Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescriptionDownMinus[i] >>= this->NbrBosons - Tmp; 
	}
      SortQuadElementArrayDownOrdering<unsigned long>(this->StateDescriptionUpPlus, this->StateDescriptionUpMinus, 
					   this->StateDescriptionDownPlus, this->StateDescriptionDownMinus, TmpLargeHilbertSpaceDimension);
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

BosonOnSquareLatticeWithSU4SpinMomentumSpace::BosonOnSquareLatticeWithSU4SpinMomentumSpace(const BosonOnSquareLatticeWithSU4SpinMomentumSpace& bosons)
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
  this->TotalSpin = bosons.TotalSpin;
  this->TotalIsospin = bosons.TotalIsospin;
  this->TotalEntanglement = bosons.TotalEntanglement;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpPlusLzMax = bosons.NUpPlusLzMax;
  this->NUpMinusLzMax = bosons.NUpMinusLzMax;
  this->NDownPlusLzMax = bosons.NDownPlusLzMax;
  this->NDownMinusLzMax = bosons.NDownMinusLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUpPlus;
  this->TemporaryStateSigma[1] = this->TemporaryStateUpMinus;
  this->TemporaryStateSigma[2] = this->TemporaryStateDownPlus;
  this->TemporaryStateSigma[3] = this->TemporaryStateDownMinus;
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUpPlus;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateUpMinus;
  this->ProdATemporaryStateSigma[2] = this->ProdATemporaryStateDownPlus;
  this->ProdATemporaryStateSigma[3] = this->ProdATemporaryStateDownMinus;
  this->StateDescriptionUpPlus = bosons.StateDescriptionUpPlus;
  this->StateDescriptionUpMinus = bosons.StateDescriptionUpMinus;
  this->StateDescriptionDownPlus = bosons.StateDescriptionDownPlus;
  this->StateDescriptionDownMinus = bosons.StateDescriptionDownMinus;
  this->StateDescriptionSigma[0] = this->StateDescriptionUpPlus;
  this->StateDescriptionSigma[1] = this->StateDescriptionUpMinus;
  this->StateDescriptionSigma[2] = this->StateDescriptionDownPlus;
  this->StateDescriptionSigma[3] = this->StateDescriptionDownMinus;
  this->NbrUniqueStateDescriptionUpPlus = bosons.NbrUniqueStateDescriptionUpPlus;
  this->UniqueStateDescriptionUpPlus = bosons.UniqueStateDescriptionUpPlus;
  this->UniqueStateDescriptionSubArraySizeUpPlus = bosons.UniqueStateDescriptionSubArraySizeUpPlus;
  this->NbrUniqueStateDescriptionUpMinus = bosons.NbrUniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionUpMinus = bosons.UniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionSubArraySizeUpMinus = bosons.UniqueStateDescriptionSubArraySizeUpMinus;
  this->FirstIndexUniqueStateDescriptionUpMinus = bosons.FirstIndexUniqueStateDescriptionUpMinus;
  this->NbrUniqueStateDescriptionDownPlus = bosons.NbrUniqueStateDescriptionDownPlus;
  this->UniqueStateDescriptionDownPlus = bosons.UniqueStateDescriptionDownPlus;
  this->UniqueStateDescriptionSubArraySizeDownPlus = bosons.UniqueStateDescriptionSubArraySizeDownPlus;
  this->FirstIndexUniqueStateDescriptionDownPlus = bosons.FirstIndexUniqueStateDescriptionDownPlus;
}

// destructor
//

BosonOnSquareLatticeWithSU4SpinMomentumSpace::~BosonOnSquareLatticeWithSU4SpinMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSquareLatticeWithSU4SpinMomentumSpace& BosonOnSquareLatticeWithSU4SpinMomentumSpace::operator = (const BosonOnSquareLatticeWithSU4SpinMomentumSpace& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescriptionUpPlus;
      delete[] this->StateDescriptionUpMinus;
      delete[] this->StateDescriptionDownPlus;
      delete[] this->StateDescriptionDownMinus;
      delete[] this->UniqueStateDescriptionUpPlus;
      delete[] this->UniqueStateDescriptionSubArraySizeUpPlus;
      delete[] this->NbrUniqueStateDescriptionUpMinus;
      for (long i = 0l; i < this->NbrUniqueStateDescriptionUpPlus; ++i)
	{
	  delete[] this->UniqueStateDescriptionUpMinus[i];
	  delete[] this->UniqueStateDescriptionSubArraySizeUpMinus[i];
	  delete[] this->FirstIndexUniqueStateDescriptionUpMinus[i];
	}
      delete[] this->UniqueStateDescriptionUpMinus;
      delete[] this->UniqueStateDescriptionSubArraySizeUpMinus;
      delete[] this->FirstIndexUniqueStateDescriptionUpMinus;
    }
  delete[] this->TemporaryStateUpPlus;
  delete[] this->TemporaryStateUpMinus;
  delete[] this->TemporaryStateDownPlus;
  delete[] this->TemporaryStateDownMinus;
  delete[] this->ProdATemporaryStateUpPlus;
  delete[] this->ProdATemporaryStateUpMinus;
  delete[] this->ProdATemporaryStateDownPlus;
  delete[] this->ProdATemporaryStateDownMinus;
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->TotalIsospin = bosons.TotalIsospin;
  this->TotalEntanglement = bosons.TotalEntanglement;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpPlusLzMax = bosons.NUpPlusLzMax;
  this->NUpMinusLzMax = bosons.NUpMinusLzMax;
  this->NDownPlusLzMax = bosons.NDownPlusLzMax;
  this->NDownMinusLzMax = bosons.NDownMinusLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUpPlus;
  this->TemporaryStateSigma[1] = this->TemporaryStateUpMinus;
  this->TemporaryStateSigma[2] = this->TemporaryStateDownPlus;
  this->TemporaryStateSigma[3] = this->TemporaryStateDownMinus;
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUpPlus;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateUpMinus;
  this->ProdATemporaryStateSigma[2] = this->ProdATemporaryStateDownPlus;
  this->ProdATemporaryStateSigma[3] = this->ProdATemporaryStateDownMinus;
  this->StateDescriptionUpPlus = bosons.StateDescriptionUpPlus;
  this->StateDescriptionUpMinus = bosons.StateDescriptionUpMinus;
  this->StateDescriptionDownPlus = bosons.StateDescriptionDownPlus;
  this->StateDescriptionDownMinus = bosons.StateDescriptionDownMinus;
  this->StateDescriptionSigma[0] = this->StateDescriptionUpPlus;
  this->StateDescriptionSigma[1] = this->StateDescriptionUpMinus;
  this->StateDescriptionSigma[2] = this->StateDescriptionDownPlus;
  this->StateDescriptionSigma[3] = this->StateDescriptionDownMinus;
  this->NbrUniqueStateDescriptionUpPlus = bosons.NbrUniqueStateDescriptionUpPlus;
  this->UniqueStateDescriptionUpPlus = bosons.UniqueStateDescriptionUpPlus;
  this->UniqueStateDescriptionSubArraySizeUpPlus = bosons.UniqueStateDescriptionSubArraySizeUpPlus;
  this->NbrUniqueStateDescriptionUpMinus = bosons.NbrUniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionUpMinus = bosons.UniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionSubArraySizeUpMinus = bosons.UniqueStateDescriptionSubArraySizeUpMinus;
  this->FirstIndexUniqueStateDescriptionUpMinus = bosons.FirstIndexUniqueStateDescriptionUpMinus;
  this->NbrUniqueStateDescriptionDownPlus = bosons.NbrUniqueStateDescriptionDownPlus;
  this->UniqueStateDescriptionDownPlus = bosons.UniqueStateDescriptionDownPlus;
  this->UniqueStateDescriptionSubArraySizeDownPlus = bosons.UniqueStateDescriptionSubArraySizeDownPlus;
  this->FirstIndexUniqueStateDescriptionDownPlus = bosons.FirstIndexUniqueStateDescriptionDownPlus;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSquareLatticeWithSU4SpinMomentumSpace::Clone()
{
  return new BosonOnSquareLatticeWithSU4SpinMomentumSpace(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSquareLatticeWithSU4SpinMomentumSpace::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[state], this->StateDescriptionUpMinus[state], 
		       this->StateDescriptionDownPlus[state], this->StateDescriptionDownMinus[state],
		       TemporaryStateUpPlus, TemporaryStateUpMinus, TemporaryStateDownPlus, TemporaryStateDownMinus); 
  
  unsigned long Tmp;
  Str << "[";
  for (int i = 0; i <= this->LzMax; ++i)
    {
      int TmpKx = i / this->NbrSiteY;
      int TmpKy = i % this->NbrSiteY;
      if (this->TemporaryStateUpPlus[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateUpPlus[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << ",u,+)";
	}
      if (this->TemporaryStateUpMinus[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateUpMinus[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << ",u,-)";
	}
      if (this->TemporaryStateDownPlus[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateDownPlus[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << ",d,+)";
	}
      if (this->TemporaryStateDownMinus[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateDownMinus[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << ",d,-)";
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
// currentFermionicPositionUpMinus = current fermionic position within the state description for the type up-plus particles
// currentFermionicPositionUpPlus = current fermionic position within the state description for the type up-minus particles
// currentFermionicPositionDownPlus = current fermionic position within the state description for the type down-plus particles
// currentFermionicPositionDownMinus = current fermionic position within the state description for the type down-minus particles
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSquareLatticeWithSU4SpinMomentumSpace::GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, 
								  int currentFermionicPositionUpPlus, int currentFermionicPositionUpMinus, 
								  int currentFermionicPositionDownPlus, int currentFermionicPositionDownMinus, long pos)
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
 	  this->StateDescriptionUpPlus[pos] = 0x0ul;
 	  this->StateDescriptionUpMinus[pos] = 0x0ul;
 	  this->StateDescriptionDownPlus[pos] = 0x0ul;
 	  this->StateDescriptionDownMinus[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  for (int i = nbrBosons; i >= 0; --i)
    {
      unsigned long Mask1 = ((0x1ul << i) - 0x1ul) << (currentFermionicPositionUpPlus - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  unsigned long Mask2 = ((0x1ul << j) - 0x1ul) << (currentFermionicPositionUpMinus - j - 1);
	  for (int k = nbrBosons - i - j; k >= 0; --k)
	    {
	      unsigned long Mask3 = ((0x1ul << k) - 0x1ul) << (currentFermionicPositionDownPlus - k - 1);
	      for (int l = nbrBosons - i - j - k; l >= 0; --l)
		{
		  long TmpPos = this->GenerateStates(nbrBosons - i - j - k - l, currentKx, currentKy - 1, currentTotalKx + ((i + j + k + l) * currentKx), currentTotalKy + ((i + j + k + l) * currentKy), currentFermionicPositionUpPlus - i - 1, currentFermionicPositionUpMinus - j - 1, currentFermionicPositionDownPlus - k - 1, currentFermionicPositionDownMinus - l - 1, pos);
		  unsigned long Mask4 = ((0x1ul << l) - 0x1ul) << (currentFermionicPositionDownMinus - l - 1);
		  for (; pos < TmpPos; ++pos)
		    {
		      this->StateDescriptionUpPlus[pos] |= Mask1;
		      this->StateDescriptionUpMinus[pos] |= Mask2;
		      this->StateDescriptionDownPlus[pos] |= Mask3;
		      this->StateDescriptionDownMinus[pos] |= Mask4;
		    }
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

long BosonOnSquareLatticeWithSU4SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
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
  for (int i = nbrBosons; i >= 0; --i)
    Count += ((((long) i + 1l) * ((long) i + 2l) * ((long) i + 3l)) / 6l) * this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKx, currentKy - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy));
  return Count;
}
