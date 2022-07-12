////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of bosons on a torus with spin                   //
//                                                                            //
//                        last modification : 03/04/2012                      //
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
#include "HilbertSpace/BosonOnTorusWithSU3Spin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "HilbertSpace/BosonOnTorusShort.h" 
#include "GeneralTools/ArrayTools.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"

#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;


// constructor with a constraint on total momentum
// 
// nbrBosons = number of bosons
// totalTz = twice the total Tz value
// totalY = three time the total Y value
// maxMomentum = momentum maximum value for a boson
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

BosonOnTorusWithSU3Spin::BosonOnTorusWithSU3Spin (int nbrBosons, int maxMomentum, int kyMomentum, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = maxMomentum;
  this->TotalY = 0;
  this->TotalTz = 0;
  this->KyMomentum = kyMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = maxMomentum;
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

  this->N1LzMax = this->LzMax + this->NbrBosons;
  this->N2LzMax = this->LzMax + this->NbrBosons;
  this->N3LzMax = this->LzMax + this->NbrBosons;
  this->FermionicLzMax = this->N1LzMax;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, 0);

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
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, 0, 
								this->LzMax + this->NbrBosons + 1, this->LzMax + this->NbrBosons + 1, 
								this->LzMax + this->NbrBosons + 1, 0l);
      cout  << "Dimension = " << this->LargeHilbertSpaceDimension << endl;
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
      UsedMemory += (long) this->HilbertSpaceDimension * (3 * sizeof(unsigned long));
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

// constructor with a constraint on total spin momentum and total momentum
// 
// nbrBosons = number of bosons
// totalTz = twice the total Tz value
// totalY = three time the total Y value
// maxMomentum = momentum maximum value for a boson
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

BosonOnTorusWithSU3Spin::BosonOnTorusWithSU3Spin (int nbrBosons, int totalTz, int totalY, int maxMomentum, int kyMomentum, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = maxMomentum;
  this->TotalY = totalY;
  this->TotalTz = totalTz;
  this->KyMomentum = kyMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = maxMomentum;
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

  int N1 = (2 * nbrBosons) + totalY + (3 * totalTz);
  int N2 = (2 * nbrBosons) + totalY - (3 * totalTz);
  int N3 = nbrBosons - totalY;
  if ((N1 < 0) || (N2 < 0) || (N3 < 0) || ((N1 % 6) != 0) || ((N2 % 6) != 0) || ((N3 % 3) != 0))
    this->LargeHilbertSpaceDimension = 0l;
  else
    {
      N1 /= 6;
      N2 /= 6;
      N3 /= 3;
      this->N1LzMax = this->LzMax + N1 - 1;
      this->N2LzMax = this->LzMax + N2 - 1;
      this->N3LzMax = this->LzMax + N3 - 1;
      this->FermionicLzMax = this->N1LzMax;
      if (this->N2LzMax > this->FermionicLzMax)
	this->FermionicLzMax = this->N2LzMax;
      if (this->N3LzMax > this->FermionicLzMax)
	this->FermionicLzMax = this->N3LzMax;
      this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, 0, N1, N2, N3);
    }
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
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, this->LzMax, this->LzMax, 0, 
								N1, N2, N3, 0l);
      cout  << "Dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      SortTripleElementArrayDownOrdering<unsigned long>(this->StateDescription1, this->StateDescription2, this->StateDescription3, TmpLargeHilbertSpaceDimension);
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (3 * sizeof(unsigned long));
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

BosonOnTorusWithSU3Spin::BosonOnTorusWithSU3Spin(const BosonOnTorusWithSU3Spin& bosons)
{
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

BosonOnTorusWithSU3Spin::~BosonOnTorusWithSU3Spin ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorusWithSU3Spin& BosonOnTorusWithSU3Spin::operator = (const BosonOnTorusWithSU3Spin& bosons)
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

AbstractHilbertSpace* BosonOnTorusWithSU3Spin::Clone()
{
  return new BosonOnTorusWithSU3Spin(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusWithSU3Spin::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescription1[state], this->StateDescription2[state], this->StateDescription3[state],
		       this->TemporaryState1, this->TemporaryState2, this->TemporaryState3); 

  unsigned long Tmp;
  Str << " | ";
  for (int i = 0; i <= this->LzMax; ++i)
    {
      Str << "(" << this->TemporaryState1[i] << "," << this->TemporaryState2[i] << "," << this->TemporaryState3[i] << ") | ";
    }
  //  Str << " : " << this->FindStateIndex(this->StateDescription1[state], this->StateDescription2[state], this->StateDescription3[state]) << " = " << state;
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKy1 = current momentum along y for a single type 1 particle
// currentKy2 = current momentum along y for a single type 2 particle
// currentKy3 = current momentum along y for a single type 3 particle
// currentTotalKy = current total momentum along y
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnTorusWithSU3Spin::GenerateStates(int nbrBosons, int currentKy1, int currentKy2, int currentKy3, int currentTotalKy, 
					     int nbrN1, int nbrN2, int nbrN3, long pos)
{
  if ((nbrBosons < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0))
    return pos;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	{
	  this->StateDescription1[pos] = 0x0ul;
	  this->StateDescription2[pos] = 0x0ul;
	  this->StateDescription3[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else
	return pos;
    }
  if ((currentKy1 < 0) || (currentKy2 < 0) || (currentKy3 < 0))
    return pos;

  long TmpPos;
  for (int i = nbrN1; i >= 0; --i)    
    {
      unsigned long Mask1 = ((0x1ul << i) - 1ul) << (currentKy1 + nbrN1 - i);
      for (int j = nbrN2; j >= 0; --j)
	{
    	  unsigned long Mask2 = ((0x1ul << j) - 1ul) << (currentKy2 + nbrN2 - j);	  
	  for (int k = nbrN3; k >= 0; --k)
	    {
	      unsigned long Mask3 = ((0x1ul << k) - 1ul) << (currentKy3 + nbrN3 - k);	  
	      TmpPos = this->GenerateStates(nbrBosons - i - j - k, currentKy1 - 1, currentKy2 - 1, currentKy3 - 1, currentTotalKy + (currentKy1 * i) + (currentKy2 * j) + (currentKy3 * k), 
					    nbrN1 - i, nbrN2 - j, nbrN3 - k, pos); 
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

// generate all states corresponding to the constraints
//! 
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle 
// currentTotalKy = current total momentum along y
// currentFermionicPositionKy1 = current fermionic position within the state description for the type 1 particles
// currentFermionicPositionKy2 = current fermionic position within the state description for the type 2 particles
// currentFermionicPositionKy3 = current fermionic position within the state description for the type 3 particles
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnTorusWithSU3Spin::GenerateStates(int nbrBosons, int currentKy, int currentTotalKy, int currentFermionicPositionKy1,
					     int currentFermionicPositionKy2, int currentFermionicPositionKy3, long pos)
{
  if (nbrBosons < 0)
    return pos;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	{
	  this->StateDescription1[pos] = 0x0ul;
	  this->StateDescription2[pos] = 0x0ul;
	  this->StateDescription3[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else
	return pos;
    }
  if (currentKy < 0)
    return pos;

  long TmpPos;
  for (int i = nbrBosons; i >= 0; --i)    
    {
      unsigned long Mask1 = ((0x1ul << i) - 0x1ul) << (currentFermionicPositionKy1 - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
    	  unsigned long Mask2 = ((0x1ul << j) - 0x1ul) << (currentFermionicPositionKy2 - j - 1);	  
	  for (int k = nbrBosons - i - j; k >= 0; --k)
	    {
	      unsigned long Mask3 = ((0x1ul << k) - 0x1ul) << (currentFermionicPositionKy3 - k - 1);	  
	      TmpPos = this->GenerateStates(nbrBosons - i - j - k, currentKy - 1, currentTotalKy + (currentKy * (i + j +k)), 
					    currentFermionicPositionKy1 - i - 1, currentFermionicPositionKy2 - j -1,
					    currentFermionicPositionKy3 - k - 1, pos); 
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
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// return value = Hilbert space dimension

long BosonOnTorusWithSU3Spin::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrN1, int nbrN2, int nbrN3)
{
  if ((nbrBosons < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0))
    return 0l;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	return 1l;
      else	
	return 0l;
    }
  if (currentKy < 0)
    return 0l;
  long Tmp = 0l;
  if (nbrBosons == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if (((j + currentTotalKy) % this->NbrLzValue) == this->KyMomentum)
	    ++Tmp;
	}
      return Tmp;
    }
  for (int i = nbrN1; i >= 0; --i)
    for (int j = nbrN2; j >= 0; --j)
      for (int k = nbrN3; k >= 0; --k)
	Tmp += this->EvaluateHilbertSpaceDimension(nbrBosons - (i + j + k), currentKy - 1, currentTotalKy + (currentKy * (i + j + k)), 
						   nbrN1 - i, nbrN2 - j, nbrN3 - k);
  return  Tmp;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonOnTorusWithSU3Spin::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy)
{
  if (nbrBosons < 0)
    return 0l;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	return 1l;
      else	
	return 0l;
    }
  if (currentKy < 0)
    return 0l;
  long Tmp = 0l;
  if (nbrBosons == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if (((j + currentTotalKy) % this->NbrLzValue) == this->KyMomentum)
	    Tmp += 3;
	}
      return Tmp;
    }
  for (int i = nbrBosons; i >= 0; --i)
    Tmp += ((((long) i + 1l) * ((long) i + 2l)) / 2l) * this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKy - 1, currentTotalKy + (currentKy * i)); 
  return  Tmp;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrN1Sector = number of type 1 particles  that belong to the subsytem 
// nbrN2Sector = number of type 1 particles  that belong to the subsytem 
// nbrN3Sector = number of type 1 particles  that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnTorusWithSU3Spin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
											    int nbrN1Sector, int nbrN2Sector, int nbrN3Sector, RealVector& groundState, AbstractArchitecture* architecture)
{
  int N1 = (2 * this->NbrBosons) + this->TotalY + (3 * this->TotalTz);
  int N2 = (2 * this->NbrBosons) + this->TotalY - (3 * this->TotalTz);
  int N3 = this->NbrBosons -  this->TotalY;
  N1 /= 6;
  N2 /= 6;
  N3 /= 3;
 if (nbrParticleSector == 0)
    {
      if ((lzSector == 0) && (nbrN1Sector == 0) && (nbrN2Sector == 0) && (nbrN3Sector == 0))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrBosons)
    {
       if ((lzSector == this->KyMomentum) && (nbrN1Sector == N1) && (nbrN2Sector == N2) && (nbrN3Sector == N3))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementaryKyMomentum = (this->KyMomentum - lzSector) % this->NbrLzValue;
  int ComplementaryN1 = N1 - nbrN1Sector;
  int ComplementaryN2 = N2 - nbrN2Sector;
  int ComplementaryN3 = N3 - nbrN3Sector;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrLzValue;
  int SubsytemTotalTz = nbrN1Sector - nbrN2Sector;
  int SubsytemTotalY = nbrN1Sector +nbrN2Sector  - (2 * nbrN3Sector);
  int ComplementaryTotalTz = ComplementaryN1 - ComplementaryN2;
  int ComplementaryTotalY = ComplementaryN1 + ComplementaryN2 - (2 * ComplementaryN3);
  cout << "ky = " << this->KyMomentum << " " << lzSector << " " << ComplementaryKyMomentum << endl;
  cout << "n1 = " << N1 << " " << nbrN1Sector << " " << ComplementaryN1 << endl;
  cout << "n2 = " << N2 << " " << nbrN2Sector << " " << ComplementaryN2 << endl;
  cout << "n3 = " << N3 << " " << nbrN3Sector << " " << ComplementaryN3 << endl;
  BosonOnTorusWithSU3Spin SubsytemSpace (nbrParticleSector, SubsytemTotalTz, SubsytemTotalY, this->NbrLzValue, lzSector);
  RealSymmetricMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnTorusWithSU3Spin ComplementarySpace (ComplementaryNbrParticles, ComplementaryTotalTz, ComplementaryTotalY,  this->NbrLzValue, ComplementaryKyMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnTorusWithSU3Spin::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
										 RealVector& groundState, RealSymmetricMatrix* densityMatrix)
{
  BosonOnTorusWithSU3Spin* TmpHilbertSpace =  (BosonOnTorusWithSU3Spin*) complementaryHilbertSpace;
  BosonOnTorusWithSU3Spin* TmpDestinationHilbertSpace =  (BosonOnTorusWithSU3Spin*) destinationHilbertSpace;
  int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
  int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  
  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
  double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[NbrBosonSector];
  for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescription1[i], TmpDestinationHilbertSpace->StateDescription2[i], TmpDestinationHilbertSpace->StateDescription3[i], TmpDestinationHilbertSpace->TemporaryState1, TmpDestinationHilbertSpace->TemporaryState2, TmpDestinationHilbertSpace->TemporaryState3); 

      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->LzMax; ++k)
	{
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState1[k]];
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState2[k]];
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState3[k]];
	}
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }

  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->StateDescription1[minIndex], TmpHilbertSpace->StateDescription2[minIndex], TmpHilbertSpace->StateDescription3[minIndex], 
				      TmpHilbertSpace->TemporaryState1, TmpHilbertSpace->TemporaryState2, TmpHilbertSpace->TemporaryState3);
       double TmpHilbertSpaceFactorial = 0.0;
       for (int k = 0; k <= TmpHilbertSpace->LzMax; ++k)
	 {
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState1[k]];
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState2[k]];
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState3[k]];
	 }
       for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	 {
	   TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescription1[j], TmpDestinationHilbertSpace->StateDescription2[j], TmpDestinationHilbertSpace->StateDescription3[j], 
						      TmpDestinationHilbertSpace->TemporaryState1, TmpDestinationHilbertSpace->TemporaryState2, TmpDestinationHilbertSpace->TemporaryState3);
	   for (int k = 0; k <=  TmpDestinationHilbertSpace->LzMax; ++k)
	     {
	       this->TemporaryState1[k] = TmpDestinationHilbertSpace->TemporaryState1[k];
	       this->TemporaryState2[k] = TmpDestinationHilbertSpace->TemporaryState2[k];
	       this->TemporaryState3[k] = TmpDestinationHilbertSpace->TemporaryState3[k];
	     }
	   for (int k = TmpDestinationHilbertSpace->LzMax + 1; k <=  this->LzMax; ++k)
	     {
	       this->TemporaryState1[k] = 0x0ul;
	       this->TemporaryState2[k] = 0x0ul;
	       this->TemporaryState3[k] = 0x0ul;
	     }	   
	   for (int k = 0; k <=  TmpHilbertSpace->LzMax; ++k)
	     {
	       this->TemporaryState1[k] += TmpHilbertSpace->TemporaryState1[k];
	       this->TemporaryState2[k] += TmpHilbertSpace->TemporaryState2[k];
	       this->TemporaryState3[k] += TmpHilbertSpace->TemporaryState3[k];
	     }

	   int TmpPos = this->FindStateIndex(this->TemporaryState1, this->TemporaryState2, this->TemporaryState3);
	   if (TmpPos != this->HilbertSpaceDimension)
	     {
	       double TmpFactorial = 0.0;	      
	       for (int k = 0; k <= this->LzMax; ++k)
		 {
		   TmpFactorial += LogFactorials[this->TemporaryState1[k]];
		   TmpFactorial += LogFactorials[this->TemporaryState2[k]];
		   TmpFactorial += LogFactorials[this->TemporaryState3[k]];
		 }
	       TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
	       TmpFactorial *= 0.5; 
	       
	       TmpStatePosition[Pos] = TmpPos;
	       TmpStatePosition2[Pos] = j;
	       TmpStateCoefficient[Pos] = exp(TmpFactorial);
	       ++Pos;
	     }
	 }
       if (Pos != 0)
 	{
 	  ++TmpNbrNonZeroElements;
 	  for (int j = 0; j < Pos; ++j)
 	    {
 	      int Pos2 = TmpStatePosition2[j];
 	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[j];
 	      for (int k = 0; k < Pos; ++k)
 		if (TmpStatePosition2[k] >= Pos2)
 		  {
 		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
 		  }
 	    }
 	}
     }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
  
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrN1Sector = number of type 1 particles  that belong to the subsytem 
// nbrN2Sector = number of type 1 particles  that belong to the subsytem 
// nbrN3Sector = number of type 1 particles  that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnTorusWithSU3Spin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
											    ComplexVector& groundState, AbstractArchitecture* architecture)
{
 if (nbrParticleSector == 0)
    {
      if (lzSector == 0)
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrBosons)
    {
       if (lzSector == this->KyMomentum)
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
    
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementaryKyMomentum = (this->KyMomentum - lzSector) % this->NbrLzValue;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrLzValue;
  cout << "ky = " << this->KyMomentum << " " << lzSector << " " << ComplementaryKyMomentum << endl;
  BosonOnTorusWithSU3Spin SubsytemSpace (nbrParticleSector, this->NbrLzValue, lzSector);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnTorusWithSU3Spin ComplementarySpace (ComplementaryNbrParticles, this->NbrLzValue, ComplementaryKyMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnTorusWithSU3Spin::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
										 ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
  BosonOnTorusWithSU3Spin* TmpHilbertSpace =  (BosonOnTorusWithSU3Spin*) complementaryHilbertSpace;
  BosonOnTorusWithSU3Spin* TmpDestinationHilbertSpace =  (BosonOnTorusWithSU3Spin*) destinationHilbertSpace;
  int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
  int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  
  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
  double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[NbrBosonSector];
  for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescription1[i], TmpDestinationHilbertSpace->StateDescription2[i], TmpDestinationHilbertSpace->StateDescription3[i], TmpDestinationHilbertSpace->TemporaryState1, TmpDestinationHilbertSpace->TemporaryState2, TmpDestinationHilbertSpace->TemporaryState3); 

      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->LzMax; ++k)
	{
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState1[k]];
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState2[k]];
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState3[k]];
	}
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }

  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->StateDescription1[minIndex], TmpHilbertSpace->StateDescription2[minIndex], TmpHilbertSpace->StateDescription3[minIndex], 
				      TmpHilbertSpace->TemporaryState1, TmpHilbertSpace->TemporaryState2, TmpHilbertSpace->TemporaryState3);
       double TmpHilbertSpaceFactorial = 0.0;
       for (int k = 0; k <= TmpHilbertSpace->LzMax; ++k)
	 {
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState1[k]];
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState2[k]];
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState3[k]];
	 }
       for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	 {
	   TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescription1[j], TmpDestinationHilbertSpace->StateDescription2[j], TmpDestinationHilbertSpace->StateDescription3[j], 
						      TmpDestinationHilbertSpace->TemporaryState1, TmpDestinationHilbertSpace->TemporaryState2, TmpDestinationHilbertSpace->TemporaryState3);
	   for (int k = 0; k <=  TmpDestinationHilbertSpace->LzMax; ++k)
	     {
	       this->TemporaryState1[k] = TmpDestinationHilbertSpace->TemporaryState1[k];
	       this->TemporaryState2[k] = TmpDestinationHilbertSpace->TemporaryState2[k];
	       this->TemporaryState3[k] = TmpDestinationHilbertSpace->TemporaryState3[k];
	     }
	   for (int k = TmpDestinationHilbertSpace->LzMax + 1; k <=  this->LzMax; ++k)
	     {
	       this->TemporaryState1[k] = 0x0ul;
	       this->TemporaryState2[k] = 0x0ul;
	       this->TemporaryState3[k] = 0x0ul;
	     }	   
	   for (int k = 0; k <=  TmpHilbertSpace->LzMax; ++k)
	     {
	       this->TemporaryState1[k] += TmpHilbertSpace->TemporaryState1[k];
	       this->TemporaryState2[k] += TmpHilbertSpace->TemporaryState2[k];
	       this->TemporaryState3[k] += TmpHilbertSpace->TemporaryState3[k];
	     }

	   int TmpPos = this->FindStateIndex(this->TemporaryState1, this->TemporaryState2, this->TemporaryState3);
	   if (TmpPos != this->HilbertSpaceDimension)
	     {
	       double TmpFactorial = 0.0;	      
	       for (int k = 0; k <= this->LzMax; ++k)
		 {
		   TmpFactorial += LogFactorials[this->TemporaryState1[k]];
		   TmpFactorial += LogFactorials[this->TemporaryState2[k]];
		   TmpFactorial += LogFactorials[this->TemporaryState3[k]];
		 }
	       TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
	       TmpFactorial *= 0.5; 
	       
	       TmpStatePosition[Pos] = TmpPos;
	       TmpStatePosition2[Pos] = j;
	       TmpStateCoefficient[Pos] = exp(TmpFactorial);
	       ++Pos;
	     }
	 }
       if (Pos != 0)
 	{
 	  ++TmpNbrNonZeroElements;
 	  for (int j = 0; j < Pos; ++j)
 	    {
 	      int Pos2 = TmpStatePosition2[j];
 	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]) * TmpStateCoefficient[j];
 	      for (int k = 0; k < Pos; ++k)
 		if (TmpStatePosition2[k] >= Pos2)
 		  {
 		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
 		  }
 	    }
 	}
     }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
  
}
