////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                 class of bosons on a cubic lattice with SU(4) spin         //
//                                in momentum space                           //
//                                                                            //
//                        last modification : 19/03/2012                      //
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
#include "HilbertSpace/BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong.h"
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
// nbrSiteZ = number of sites in the z direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// kzMomentum = momentum along the z direction
// memory = amount of memory granted for precalculations

BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong::BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong (int nbrBosons, int nbrSiteX, int nbrSiteY, int nbrSiteZ, int kxMomentum, int kyMomentum, int kzMomentum, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->TotalIsospin = 0;
  this->TotalEntanglement = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteYZ = this->NbrSiteZ * this->NbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->KzMomentum = kzMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
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
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, 0, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUpPlus = new ULONGLONG [this->LargeHilbertSpaceDimension];
      this->StateDescriptionUpMinus = new ULONGLONG [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDownPlus = new ULONGLONG [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDownMinus = new ULONGLONG [this->LargeHilbertSpaceDimension];
      this->StateDescriptionSigma[0] = this->StateDescriptionUpPlus;
      this->StateDescriptionSigma[1] = this->StateDescriptionUpMinus;
      this->StateDescriptionSigma[2] = this->StateDescriptionDownPlus;
      this->StateDescriptionSigma[3] = this->StateDescriptionDownMinus;
      this->Flag.Initialize();
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, this->NbrSiteZ - 1, 0, 0, 0, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      for (long i = 0; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  ULONGLONG TmpState = this->StateDescriptionUpPlus[i];
	  ULONGLONG Tmp = (ULONGLONG)0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & (ULONGLONG)0x1ul;
	  this->StateDescriptionUpPlus[i] >>= this->NbrBosons - Tmp; 
	  TmpState = this->StateDescriptionUpMinus[i];
	  Tmp = (ULONGLONG)0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & (ULONGLONG)0x1ul;
	  this->StateDescriptionUpMinus[i] >>= this->NbrBosons - Tmp; 
	  TmpState = this->StateDescriptionDownPlus[i];
	  Tmp = (ULONGLONG)0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & (ULONGLONG)0x1ul;
	  this->StateDescriptionDownPlus[i] >>= this->NbrBosons - Tmp; 
	  TmpState = this->StateDescriptionDownMinus[i];
	  Tmp = (ULONGLONG)0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & (ULONGLONG)0x1ul;
	  this->StateDescriptionDownMinus[i] >>= this->NbrBosons - Tmp; 
	}
      SortQuadElementArrayDownOrdering<ULONGLONG>(this->StateDescriptionUpPlus, this->StateDescriptionUpMinus, 
					   this->StateDescriptionDownPlus, this->StateDescriptionDownMinus, TmpLargeHilbertSpaceDimension);
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (4 * sizeof(ULONGLONG));
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

BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong::BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong(const BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong& bosons)
{
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->NbrSiteZ = bosons.NbrSiteZ;
  this->NbrSiteYZ = bosons.NbrSiteYZ;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->KzMomentum = bosons.KzMomentum;
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

BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong::~BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong& BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong::operator = (const BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong& bosons)
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
  this->NbrSiteZ = bosons.NbrSiteZ;
  this->NbrSiteYZ = bosons.NbrSiteYZ;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->KzMomentum = bosons.KzMomentum;
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

AbstractHilbertSpace* BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong::Clone()
{
  return new BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[state], this->StateDescriptionUpMinus[state], 
		       this->StateDescriptionDownPlus[state], this->StateDescriptionDownMinus[state],
		       this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, this->TemporaryStateDownPlus, this->TemporaryStateDownMinus); 
  
  Str << "[";
  for (int i = 0; i <= this->LzMax; ++i)
    {
      int TmpKx = i / this->NbrSiteYZ;
      int TmpKy = i % this->NbrSiteYZ;
      int TmpKz = TmpKy % this->NbrSiteZ;
      TmpKy /= this->NbrSiteZ;
      if (this->TemporaryStateUpPlus[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateUpPlus[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << ",u,A)";
	}
      if (this->TemporaryStateUpMinus[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateUpMinus[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << ",u,B)";
	}
      if (this->TemporaryStateDownPlus[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateDownPlus[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << ",d,A)";
	}
      if (this->TemporaryStateDownMinus[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateDownMinus[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << "," << TmpKz << ",d,B)";
	}
// 	Str << this->TemporaryStateUpPlus[i] << " " << this->TemporaryStateUpMinus[i] << " " << this->TemporaryStateDownPlus[i] << " " << this->TemporaryStateDownMinus[i] ;
    }
  Str << "]";
  
//   Str << "   " << state << " " <<  this->FindStateIndex(this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, 
// 							this->TemporaryStateDownPlus, this->TemporaryStateDownMinus)
//       << " " << this->FindStateIndex(this->StateDescriptionUpPlus[state], this->StateDescriptionUpMinus[state], 
// 				     this->StateDescriptionDownPlus[state], this->StateDescriptionDownMinus[state]);
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// currentFermionicPositionUpMinus = current fermionic position within the state description for the type up-plus particles
// currentFermionicPositionUpPlus = current fermionic position within the state description for the type up-minus particles
// currentFermionicPositionDownPlus = current fermionic position within the state description for the type down-plus particles
// currentFermionicPositionDownMinus = current fermionic position within the state description for the type down-minus particles
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong::GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentKz, 
								 int currentTotalKx, int currentTotalKy, int currentTotalKz, 
								 int currentFermionicPositionUpPlus, int currentFermionicPositionUpMinus, 
								 int currentFermionicPositionDownPlus, int currentFermionicPositionDownMinus, long pos)
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
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum) && 
	  ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum))
	{
 	  this->StateDescriptionUpPlus[pos] = (ULONGLONG)0x0ul;
 	  this->StateDescriptionUpMinus[pos] = (ULONGLONG)0x0ul;
 	  this->StateDescriptionDownPlus[pos] = (ULONGLONG)0x0ul;
 	  this->StateDescriptionDownMinus[pos] = (ULONGLONG)0x0ul;
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  for (int i = nbrBosons; i >= 0; --i)
    {
      ULONGLONG Mask1 = (((ULONGLONG)0x1ul << i) - (ULONGLONG)0x1ul) << (currentFermionicPositionUpPlus - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  ULONGLONG Mask2 = (((ULONGLONG)0x1ul << j) - (ULONGLONG)0x1ul) << (currentFermionicPositionUpMinus - j - 1);
	  for (int k = nbrBosons - i - j; k >= 0; --k)
	    {
	      ULONGLONG Mask3 = (((ULONGLONG)0x1ul << k) - (ULONGLONG)0x1ul) << (currentFermionicPositionDownPlus - k - 1);
	      for (int l = nbrBosons - i - j - k; l >= 0; --l)
		{
		  long TmpPos = this->GenerateStates(nbrBosons - i - j - k - l, currentKx, currentKy, currentKz - 1, currentTotalKx + ((i + j + k + l) * currentKx), currentTotalKy + ((i + j + k + l) * currentKy), currentTotalKz + ((i + j + k + l) * currentKz), currentFermionicPositionUpPlus - i - 1, currentFermionicPositionUpMinus - j - 1, currentFermionicPositionDownPlus - k - 1, currentFermionicPositionDownMinus - l - 1, pos);
		  ULONGLONG Mask4 = (((ULONGLONG)0x1ul << l) - (ULONGLONG)0x1ul) << (currentFermionicPositionDownMinus - l - 1);
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
// currentKz = current momentum along z for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentTotalKz = current total momentum along z
// return value = Hilbert space dimension

long BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentKz,
										int currentTotalKx, int currentTotalKy, int currentTotalKz)
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
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum) &&
	  ((currentTotalKz % this->NbrSiteZ) == this->KzMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  for (int i = nbrBosons; i >= 0; --i)
    Count += ((((long) i + 1l) * ((long) i + 2l) * ((long) i + 3l)) / 6l) * this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKx, currentKy, currentKz - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy), currentTotalKz + (i * currentKz));
  return Count;
}
