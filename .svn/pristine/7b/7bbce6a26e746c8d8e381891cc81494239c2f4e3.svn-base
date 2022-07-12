////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin without            //
//                            sign precalculation table                       //
//        supporting system partial polarization (i.e. an infinite Zeeman     //
//                            energy on some orbitals)                        //
//                                                                            //
//                        last modification : 03/10/2018                      //
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
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpinPartialPolarization.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/ArrayTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/StringTools.h"

#include <math.h>
#include <cstdlib>
#include <fstream>
#include <map>
#include <bitset>
#include <algorithm>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::map;
using std::pair;
using std::bitset;

#define WANT_LAPACK

#ifdef __LAPACK__
#ifdef WANT_LAPACK
#define  __USE_LAPACK_HERE__
#endif
#endif


// default constructor
//

FermionOnSphereWithSpinPartialPolarization::FermionOnSphereWithSpinPartialPolarization()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twce the total spin value
// nbrSpinPolarizedOrbitals = number of orbitals (from the lowest momenta) that are fully spin polarized (i.e. where only spin up are allowed)
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinPartialPolarization::FermionOnSphereWithSpinPartialPolarization (int nbrFermions, int totalLz, int lzMax, int totalSpin, int nbrSpinPolarizedOrbitals, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->NbrSpinPolarizedOrbitals = nbrSpinPolarizedOrbitals;
  if(this->NbrFermions > 0)
    {
      this->LargeHilbertSpaceDimension = (int) this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 
											  (this->NbrFermions + this->TotalSpin) >> 1, (this->NbrFermions - this->TotalSpin) >> 1);
      if (this->LargeHilbertSpaceDimension >= (1l << 30))
	this->HilbertSpaceDimension = 0;
      else
	this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
    }
  else
    {
      this->HilbertSpaceDimension = 1l;
    }
  
  this->Flag.Initialize();
  this->TargetSpace = this;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateHighestBit = new int [this->LargeHilbertSpaceDimension];  
      
      if (this->NbrFermions > 0)
	{
	  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1,
							     (this->NbrFermions + this->TotalSpin) >> 1, (this->NbrFermions - this->TotalSpin) >> 1, 0l);
	  this->GenerateLookUpTable(memory);
	}
      else
	{
	  this->StateDescription[0] = 0x0ul; 
	}
    }
  else
    {
      this->StateDescription = 0;
      this->StateHighestBit = 0;
      this->LookUpTableMemorySize = 0;
    }

  
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
  if (this->NbrFermions > 0)
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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereWithSpinPartialPolarization::FermionOnSphereWithSpinPartialPolarization(const FermionOnSphereWithSpinPartialPolarization& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
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
  this->NbrSpinPolarizedOrbitals = fermions.NbrSpinPolarizedOrbitals;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnSphereWithSpinPartialPolarization::~FermionOnSphereWithSpinPartialPolarization ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereWithSpinPartialPolarization& FermionOnSphereWithSpinPartialPolarization::operator = (const FermionOnSphereWithSpinPartialPolarization& fermions)
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
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->NbrSpinPolarizedOrbitals = fermions.NbrSpinPolarizedOrbitals;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereWithSpinPartialPolarization::Clone()
{
  return new FermionOnSphereWithSpinPartialPolarization(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnSphereWithSpinPartialPolarization::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
  this->TargetSpace = (FermionOnSphereWithSpinPartialPolarization*) targetSpace;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// nbrNUp = number of particles with quantum number up
// nbrNDown = number of particles with quantum number down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereWithSpinPartialPolarization::GenerateStates(int nbrFermions, int lzMax, int totalLz, int nbrNUp, int nbrNDown, long pos)
{
  if ((nbrFermions < 0) || (totalLz < 0)  || (nbrNUp < 0) || (nbrNDown < 0))
    return pos;
  if ((nbrFermions == 0) && (totalLz == 0) && (nbrNUp == 0) && (nbrNDown == 0))
      {
	this->StateDescription[pos] = 0x0ul;
	return (pos + 1l);
      }
    
  if ((lzMax < 0) || (nbrNUp > (lzMax + 1)) || (nbrNDown > (lzMax + 1)) || ((lzMax < this->NbrSpinPolarizedOrbitals) && (nbrNDown > this->NbrSpinPolarizedOrbitals)))
    return pos;
    
  long TmpPos = 0l;
  if (this->NbrSpinPolarizedOrbitals <= lzMax)
    {
      TmpPos = this->GenerateStates(nbrFermions - 2, lzMax - 1, totalLz - (lzMax << 1), nbrNUp - 1, nbrNDown - 1,  pos);
      unsigned long Mask = 0x3ul << (lzMax << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
      TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrNUp - 1,  nbrNDown, pos);
      Mask = 0x2ul << (lzMax << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
      TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrNUp,  nbrNDown - 1, pos);
      Mask = 0x1ul << (lzMax << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  else
    {
      TmpPos = this->GenerateStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrNUp - 1, nbrNDown ,  pos);
      unsigned long Mask = 0x2ul << (lzMax << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  return this->GenerateStates(nbrFermions, lzMax - 1, totalLz, nbrNUp,  nbrNDown, pos);
};


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// nbrNUp = number of particles with quantum number up
// nbrNDown = number of particles with quantum number down
// return value = Hilbert space dimension

long FermionOnSphereWithSpinPartialPolarization::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int nbrNUp, int nbrNDown)
{
  if ((nbrFermions < 0) || (totalLz < 0) || (nbrNUp < 0) || (nbrNDown < 0))
    return 0l;
  if ((nbrFermions == 0) && (totalLz == 0))
    return 1l;
  if ((lzMax < 0) || (nbrNUp > (lzMax + 1)) || (nbrNDown > (lzMax + 1)) || ((lzMax < this->NbrSpinPolarizedOrbitals) && (nbrNDown > this->NbrSpinPolarizedOrbitals)))
    return 0l;
    
  unsigned long Tmp = 0l;  
  if (this->NbrSpinPolarizedOrbitals <= lzMax)
    {
     if (nbrFermions >= 2)
       Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), nbrNUp - 1, nbrNDown - 1);
     Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrNUp - 1, nbrNDown);
     Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrNUp, nbrNDown - 1);
     Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, nbrNUp, nbrNDown);
    }
  else
    {
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax, nbrNUp - 1, nbrNDown);
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz, nbrNUp, nbrNDown);
    }
  return Tmp;
}



