////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                    class of bosons on sphere including three               //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 15/12/2011                      //
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
#include "HilbertSpace/BosonOnSphereThreeLandauLevels.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/ArrayTools.h"
#include "MathTools/FactorialCoefficient.h"

#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <map>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::map;
using std::pair;


// default constructor
// 

BosonOnSphereThreeLandauLevels::BosonOnSphereThreeLandauLevels ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// memory = amount of memory granted for precalculations

BosonOnSphereThreeLandauLevels::BosonOnSphereThreeLandauLevels (int nbrBosons, int totalLz, int lzMax, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalY = 0;
  this->TotalTz = 0;
  this->LzMax = lzMax + 4;
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->N1LzMax = this->LzMax + this->NbrBosons - 1;
  this->N2LzMax = this->LzMax + this->NbrBosons - 1;
  this->N3LzMax = this->LzMax + this->NbrBosons - 1;
  this->FermionicLzMax = this->N1LzMax;
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

  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax - 4, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1);
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
      //            long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, 0l);
   long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax - 4, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 0l);

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

BosonOnSphereThreeLandauLevels::BosonOnSphereThreeLandauLevels(const BosonOnSphereThreeLandauLevels& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->N1LzMax = bosons.LzMax;
  this->N2LzMax = bosons.LzMax;
  this->N3LzMax = bosons.LzMax;
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

BosonOnSphereThreeLandauLevels::~BosonOnSphereThreeLandauLevels ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereThreeLandauLevels& BosonOnSphereThreeLandauLevels::operator = (const BosonOnSphereThreeLandauLevels& bosons)
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

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->N1LzMax = bosons.LzMax;
  this->N2LzMax = bosons.LzMax;
  this->N3LzMax = bosons.LzMax;
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

AbstractHilbertSpace* BosonOnSphereThreeLandauLevels::Clone()
{
  return new BosonOnSphereThreeLandauLevels(*this);
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// nbrFluxQuanta = number of flux quanta
// lzMax = momentum maximum value for a boson in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereThreeLandauLevels::GenerateStates(int nbrBosons, int nbrFluxQuanta, int lzMax, int totalLz, long pos)
{
  if ((nbrBosons < 0) || (totalLz < 0))
    return pos;
  if ((nbrBosons == 0) && (totalLz == 0))
    {
      //      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
  if (lzMax < 0) 
    return pos;
    
//   if (nbrBosons == 1) 
//     {
//       if (lzMax >= totalLz)
// 	{
// 	  if (((nbrFluxQuanta + 4) >= totalLz) && (totalLz >= 0))
// 	    {
// 	      this->StateDescription[pos] = 0x4ul << (totalLz * 3);
// 	      ++pos;
// 	    }
// 	  if (((nbrFluxQuanta + 3) >= totalLz) && (totalLz >= 1))
// 	    {
// 	      this->StateDescription[pos] = 0x2ul << (totalLz * 3);
// 	      ++pos;
// 	    }
// 	  if (((nbrFluxQuanta + 2) >= totalLz) && (totalLz >= 2))
// 	    {
// 	      this->StateDescription[pos] = 0x1ul << (totalLz * 3);
// 	      ++pos;
// 	    }
// 	}
//       return pos;
//     }

  if ((lzMax == 0) && (totalLz != 0))
    return pos;

  long TmpPos = pos;
//   unsigned long Mask = 0x0ul;
//   if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 0))
//     {
//       if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 1))
// 	{
// 	  if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
// 	    {
// 	      TmpPos = this->GenerateStates(nbrBosons - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax), pos);
// 	      Mask = 0x7ul << (lzMax * 3);
// 	      for (; pos < TmpPos; ++pos)
// 		this->StateDescription[pos] |= Mask;
// 	    }
// 	  TmpPos = this->GenerateStates(nbrBosons - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
// 	  Mask = 0x6ul << (lzMax * 3);
// 	  for (; pos < TmpPos; ++pos)
// 	    this->StateDescription[pos] |= Mask;
// 	}
//       if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
// 	{
// 	  TmpPos = this->GenerateStates(nbrBosons - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
// 	  Mask = 0x5ul << (lzMax * 3);
// 	  for (; pos < TmpPos; ++pos)
// 	    this->StateDescription[pos] |= Mask;
// 	}
//       TmpPos = this->GenerateStates(nbrBosons - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax, pos);
//       Mask = 0x4ul << (lzMax * 3);
//       for (; pos < TmpPos; ++pos)
// 	this->StateDescription[pos] |= Mask;
//     }
//   if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 1))
//     {
//       if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
// 	{
// 	  TmpPos = this->GenerateStates(nbrBosons - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
// 	  Mask = 0x3ul << (lzMax * 3);
// 	  for (; pos < TmpPos; ++pos)
// 	    this->StateDescription[pos] |= Mask;
// 	}
//       TmpPos = this->GenerateStates(nbrBosons - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax, pos);
//       Mask = 0x2ul << (lzMax * 3);
//       for (; pos < TmpPos; ++pos)
// 	this->StateDescription[pos] |= Mask;
//     }
//   if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))    
//     {
//       TmpPos = this->GenerateStates(nbrBosons - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax, pos);
//       Mask = 0x1ul << (lzMax * 3);
//       for (; pos < TmpPos; ++pos)
// 	this->StateDescription[pos] |= Mask;
//     }
  return this->GenerateStates(nbrBosons, nbrFluxQuanta, lzMax - 1, totalLz, pos);
}


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// nbrFluxQuanta = number of flux quanta
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

long BosonOnSphereThreeLandauLevels::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int nbrFluxQuanta, int lzMax, int totalLz)
{
  if ((nbrBosons < 0) || (totalLz < 0))
    return 0l;
  if ((nbrBosons == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if (nbrBosons == 1) 
    {
      long Tmp = 0l;
      if (lzMax >= totalLz)
	{
	  if (((nbrFluxQuanta + 4) >= totalLz) && (totalLz >= 0))
	    ++Tmp;
	  if (((nbrFluxQuanta + 3) >= totalLz) && (totalLz >= 1))
	    ++Tmp;
	  if (((nbrFluxQuanta + 2) >= totalLz) && (totalLz >= 2))
	    ++Tmp;
	}
      return Tmp;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;

  long Tmp = 0l;
  if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 0))
    {
      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 1))
	{
	  if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax));
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
	}
      if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 1))
    {
      if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))    
    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, nbrFluxQuanta, lzMax - 1, totalLz);
  return Tmp;
}

