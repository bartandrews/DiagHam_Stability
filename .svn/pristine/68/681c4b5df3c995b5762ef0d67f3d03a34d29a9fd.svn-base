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
#include "HilbertSpace/BosonOnTorusWithSpin.h"
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


// constructor with a constraint on the total Ky momentum
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

BosonOnTorusWithSpin::BosonOnTorusWithSpin (int nbrBosons, int maxMomentum, int kyMomentum, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->SzFlag = false;
  this->TotalLz = maxMomentum;
  this->TotalSpin = 0;
  this->NbrBosonsUp = 0;
  this->NbrBosonsDown = 0;
  this->KyMomentum = kyMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = maxMomentum;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;

  this->NUpLzMax = this->LzMax + this->NbrBosons;
  this->NDownLzMax = this->LzMax + this->NbrBosons;
  this->FermionicLzMax = this->NUpLzMax;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->TargetSpace = this;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionSigma[0] = this->StateDescriptionUp;
      this->StateDescriptionSigma[1] = this->StateDescriptionDown;
      this->Flag.Initialize();
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, 0, this->LzMax + this->NbrBosons + 1, this->LzMax + this->NbrBosons + 1, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      cout  << "Dimension = " << this->LargeHilbertSpaceDimension << endl;
      for (long i = 0; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->StateDescriptionUp[i];
	  unsigned long Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescriptionUp[i] >>= this->NbrBosons - Tmp; 
	  TmpState = this->StateDescriptionDown[i];
	  Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescriptionDown[i] >>= this->NbrBosons - Tmp; 
	}
      SortDoubleElementArrayDownOrdering<unsigned long>(this->StateDescriptionUp, this->StateDescriptionDown, TmpLargeHilbertSpaceDimension);
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

// constructor with a constraint on total spin momentum and total momentum
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// totalSpin = twice the total spin along z
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

BosonOnTorusWithSpin::BosonOnTorusWithSpin (int nbrBosons, int maxMomentum, int totalSpin, int kyMomentum, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->SzFlag = true;
  this->TotalLz = maxMomentum;
  this->TotalSpin = totalSpin;
  this->NbrBosonsUp = (this->NbrBosons + this->TotalSpin) / 2;
  this->NbrBosonsDown = (this->NbrBosons - this->TotalSpin) / 2;
  this->KyMomentum = kyMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = maxMomentum;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;

  this->NUpLzMax = this->LzMax + this->NbrBosons;
  this->NDownLzMax = this->LzMax + this->NbrBosons;
  this->FermionicLzMax = this->NUpLzMax;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, 0, this->NbrBosonsUp);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->TargetSpace = this;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionSigma[0] = this->StateDescriptionUp;
      this->StateDescriptionSigma[1] = this->StateDescriptionDown;
      this->Flag.Initialize();
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, 0, this->LzMax + this->NbrBosonsUp + 1, this->LzMax + this->NbrBosonsDown + 1, this->NbrBosonsUp, 0l);
      cout  << "Dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      SortDoubleElementArrayDownOrdering<unsigned long>(this->StateDescriptionUp, this->StateDescriptionDown, TmpLargeHilbertSpaceDimension);
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

BosonOnTorusWithSpin::BosonOnTorusWithSpin(const BosonOnTorusWithSpin& bosons)
{
  this->KyMomentum = bosons.KyMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->SzFlag = bosons.SzFlag;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

BosonOnTorusWithSpin::~BosonOnTorusWithSpin ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorusWithSpin& BosonOnTorusWithSpin::operator = (const BosonOnTorusWithSpin& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescriptionUp;
      delete[] this->StateDescriptionDown;
      delete[] this->UniqueStateDescriptionUp;
      delete[] this->UniqueStateDescriptionSubArraySizeUp;
      delete[] this->FirstIndexUniqueStateDescriptionUp;
    }
  delete[] this->TemporaryStateUp;
  delete[] this->TemporaryStateDown;
  delete[] this->ProdATemporaryStateUp;
  delete[] this->ProdATemporaryStateDown;
  this->KyMomentum = bosons.KyMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->SzFlag = bosons.SzFlag;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnTorusWithSpin::Clone()
{
  return new BosonOnTorusWithSpin(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusWithSpin::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUp[state], this->StateDescriptionDown[state],
		       this->TemporaryStateUp, this->TemporaryStateDown); 

  unsigned long Tmp;
  Str << " | ";
  for (int i = 0; i <= this->LzMax; ++i)
    {
      Str << "(" << this->TemporaryStateUp[i] << "," << this->TemporaryStateDown[i] << ") | ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnTorusWithSpin::GenerateStates(int nbrBosons, int currentKy, int currentTotalKy, int currentFermionicPositionUp, int currentFermionicPositionDown, long pos)
{
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	{
 	  this->StateDescriptionUp[pos] = 0x0ul;
 	  this->StateDescriptionDown[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKy < 0)
    return pos;
  for (int i = nbrBosons; i >= 0; --i)
    {
      unsigned long MaskUp = ((0x1ul << i) - 0x1ul) << (currentFermionicPositionUp - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  long TmpPos = this->GenerateStates(nbrBosons - i - j, currentKy - 1, currentTotalKy + ((i + j) * currentKy), currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, pos);
	  unsigned long MaskDown = ((0x1ul << j) - 0x1ul) << (currentFermionicPositionDown - j - 1);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescriptionUp[pos] |= MaskUp;
	      this->StateDescriptionDown[pos] |= MaskDown;
	    }
	}
    }
  return pos;
};

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// nbrSpinUp = number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnTorusWithSpin::GenerateStates(int nbrBosons, int currentKy, int currentTotalKy, 
					  int currentFermionicPositionUp, int currentFermionicPositionDown, int nbrSpinUp, long pos)
{
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrBosons))
    return 0l;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	{
 	  this->StateDescriptionUp[pos] = 0x0ul;
 	  this->StateDescriptionDown[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKy < 0)
    return pos;
  for (int i = nbrSpinUp; i >= 0; --i)
    {
      unsigned long MaskUp = ((0x1ul << i) - 0x1ul) << (currentFermionicPositionUp - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  long TmpPos = this->GenerateStates(nbrBosons - i - j, currentKy - 1, currentTotalKy + ((i + j) * currentKy), currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, nbrSpinUp - i, pos);
	  unsigned long MaskDown = ((0x1ul << j) - 0x1ul) << (currentFermionicPositionDown - j - 1);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescriptionUp[pos] |= MaskUp;
	      this->StateDescriptionDown[pos] |= MaskDown;
	    }
	}
    }
  return pos;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonOnTorusWithSpin::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy)
{
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	return 1l;
      else	
	return 0l;
    }
  if (currentKy < 0)
    return 0l;
  long Count = 0;
  if (nbrBosons == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if (((j + currentTotalKy) % this->NbrLzValue) == this->KyMomentum)
	    Count += 2l;
	}
      return Count;
    }
  for (int i = nbrBosons; i >= 0; --i)
    Count += (((long) i) + 1l) * this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKy - 1, currentTotalKy + (i * currentKy));
  return Count;
}

// evaluate Hilbert space dimension for a given total spin momentum
//
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// nbrSpinUp = number of particles with spin up
// return value = Hilbert space dimension

long BosonOnTorusWithSpin::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrSpinUp)
{
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrBosons))
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
  long Count = 0;
  if (nbrBosons == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if (((j + currentTotalKy) % this->NbrLzValue) == this->KyMomentum)
	    Count++;
	}
      return Count;
    }
  for (int i = nbrBosons; i >= 0; --i)
    for (int j = i; j >= 0; --j)
      Count += this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKy - 1, currentTotalKy + (i * currentKy), nbrSpinUp -j);
  return Count;
}


// project out any configurations that have particles on levels other than lll
//
// inputVector = vector to apply the projection to
// outputVector = projected vector
// finalSpace = reference to output vector space

void BosonOnTorusWithSpin::ProjectionInTheLowestLevel(RealVector &inputVector, RealVector & outputVector, BosonOnTorusShort *finalSpace)
{
  unsigned long Etat;
  int Idx; 
  for(int i = 0 ; i < finalSpace->GetHilbertSpaceDimension() ; i++)
    {
      Etat = finalSpace->StateDescription[i]; 
      Idx = this->FindStateIndex(0,Etat);
      if ( Idx < this->HilbertSpaceDimension ) 
	{
	  outputVector[i] = inputVector[Idx];
	}
    }
  cout <<"Norm after projection" <<outputVector.Norm()<<endl;
}


void BosonOnTorusWithSpin::ProjectionInTheLowestLevel(ComplexVector &inputVector, ComplexVector & outputVector, BosonOnTorusShort *finalSpace)
{
  unsigned long Etat;
  int Idx; 
  for(int i = 0 ; i < finalSpace->GetHilbertSpaceDimension() ; i++)
    {
      Etat = finalSpace->StateDescription[i]; 
      Idx = this->FindStateIndex(Etat,0);
      if ( Idx < this->HilbertSpaceDimension ) 
	{
	  cout <<i<< " "<<Idx<<endl;
	  outputVector[i] = inputVector[Idx];
	}
    }
  cout <<"Norm after projection" <<outputVector.Norm()<<endl;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnTorusWithSpin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
											 RealVector& groundState, AbstractArchitecture* architecture)
{
 if (nbrParticleSector == 0)
    {
      if (lzSector == 0)
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
       if (lzSector == this->KyMomentum)
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
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrLzValue;
  cout << "ky = " << this->KyMomentum << " " << lzSector << " " << ComplementaryKyMomentum << endl;
  BosonOnTorusWithSpin SubsytemSpace (nbrParticleSector, this->NbrLzValue, lzSector);
  RealSymmetricMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnTorusWithSpin ComplementarySpace (ComplementaryNbrParticles, this->NbrLzValue, ComplementaryKyMomentum);
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

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnTorusWithSpin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
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
  BosonOnTorusWithSpin SubsytemSpace (nbrParticleSector, this->NbrLzValue, lzSector);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnTorusWithSpin ComplementarySpace (ComplementaryNbrParticles, this->NbrLzValue, ComplementaryKyMomentum);
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


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrNUpSector = number of spin up  that belong to the subsytem 
// nbrNDownSector = number of spin down  that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnTorusWithSpin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
											 int nbrNUpSector, int nbrNDownSector, RealVector& groundState, AbstractArchitecture* architecture)
{
 if (nbrParticleSector == 0)
    {
      if ((lzSector == 0) && (nbrNUpSector == 0) && (nbrNDownSector == 0))
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
       if ((lzSector == this->KyMomentum) && (nbrNUpSector == this->NbrBosonsUp) && (nbrNDownSector == this->NbrBosonsDown))
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
  int ComplementaryNUp = this->NbrBosonsUp - nbrNUpSector;
  int ComplementaryNDown = this->NbrBosonsDown - nbrNDownSector;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrLzValue;
  int SubsytemTotalSz = nbrNUpSector - nbrNDownSector;
  int ComplementaryTotalSz = ComplementaryNUp - ComplementaryNDown;
  cout << "ky = " << this->KyMomentum << " " << lzSector << " " << ComplementaryKyMomentum << endl;
  cout << "nup = " << this->NbrBosonsUp << " " << nbrNUpSector << " " << ComplementaryNUp << endl;
  cout << "ndown = " << this->NbrBosonsDown << " " << nbrNDownSector << " " << ComplementaryNDown << endl;
  BosonOnTorusWithSpin SubsytemSpace (nbrParticleSector, this->NbrLzValue, SubsytemTotalSz, lzSector);
  RealSymmetricMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnTorusWithSpin ComplementarySpace (ComplementaryNbrParticles,  this->NbrLzValue, ComplementaryTotalSz, ComplementaryKyMomentum);
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

long BosonOnTorusWithSpin::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
									      RealVector& groundState, RealSymmetricMatrix* densityMatrix)
{
  BosonOnTorusWithSpin* TmpHilbertSpace =  (BosonOnTorusWithSpin*) complementaryHilbertSpace;
  BosonOnTorusWithSpin* TmpDestinationHilbertSpace =  (BosonOnTorusWithSpin*) destinationHilbertSpace;
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
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescriptionUp[i], TmpDestinationHilbertSpace->StateDescriptionDown[i],
						 TmpDestinationHilbertSpace->TemporaryStateUp, TmpDestinationHilbertSpace->TemporaryStateDown); 

      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->LzMax; ++k)
	{
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryStateUp[k]];
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryStateDown[k]];
	}
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }

  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->StateDescriptionUp[minIndex], TmpHilbertSpace->StateDescriptionDown[minIndex],
				      TmpHilbertSpace->TemporaryStateUp, TmpHilbertSpace->TemporaryStateDown);
       double TmpHilbertSpaceFactorial = 0.0;
       for (int k = 0; k <= TmpHilbertSpace->LzMax; ++k)
	 {
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryStateUp[k]];
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryStateDown[k]];
	 }
       for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	 {
	   TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescriptionUp[j], TmpDestinationHilbertSpace->StateDescriptionDown[j], 
						      TmpDestinationHilbertSpace->TemporaryStateUp, TmpDestinationHilbertSpace->TemporaryStateDown);
	   for (int k = 0; k <=  TmpDestinationHilbertSpace->LzMax; ++k)
	     {
	       this->TemporaryStateUp[k] = TmpDestinationHilbertSpace->TemporaryStateUp[k];
	       this->TemporaryStateDown[k] = TmpDestinationHilbertSpace->TemporaryStateDown[k];
	     }
	   for (int k = TmpDestinationHilbertSpace->LzMax + 1; k <=  this->LzMax; ++k)
	     {
	       this->TemporaryStateUp[k] = 0x0ul;
	       this->TemporaryStateDown[k] = 0x0ul;
	     }	   
	   for (int k = 0; k <=  TmpHilbertSpace->LzMax; ++k)
	     {
	       this->TemporaryStateUp[k] += TmpHilbertSpace->TemporaryStateUp[k];
	       this->TemporaryStateDown[k] += TmpHilbertSpace->TemporaryStateDown[k];
	     }

	   int TmpPos = this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
	   if (TmpPos != this->HilbertSpaceDimension)
	     {
	       double TmpFactorial = 0.0;	      
	       for (int k = 0; k <= this->LzMax; ++k)
		 {
		   TmpFactorial += LogFactorials[this->TemporaryStateUp[k]];
		   TmpFactorial += LogFactorials[this->TemporaryStateDown[k]];
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

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnTorusWithSpin::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
									      ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
  BosonOnTorusWithSpin* TmpHilbertSpace =  (BosonOnTorusWithSpin*) complementaryHilbertSpace;
  BosonOnTorusWithSpin* TmpDestinationHilbertSpace =  (BosonOnTorusWithSpin*) destinationHilbertSpace;
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
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescriptionUp[i], TmpDestinationHilbertSpace->StateDescriptionDown[i],
						 TmpDestinationHilbertSpace->TemporaryStateUp, TmpDestinationHilbertSpace->TemporaryStateDown); 

      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->LzMax; ++k)
	{
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryStateUp[k]];
	  TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryStateDown[k]];
	}
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }

  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->StateDescriptionUp[minIndex], TmpHilbertSpace->StateDescriptionDown[minIndex],
				      TmpHilbertSpace->TemporaryStateUp, TmpHilbertSpace->TemporaryStateDown);
       double TmpHilbertSpaceFactorial = 0.0;
       for (int k = 0; k <= TmpHilbertSpace->LzMax; ++k)
	 {
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryStateUp[k]];
	   TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryStateDown[k]];
	 }
       for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	 {
	   TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->StateDescriptionUp[j], TmpDestinationHilbertSpace->StateDescriptionDown[j], 
						      TmpDestinationHilbertSpace->TemporaryStateUp, TmpDestinationHilbertSpace->TemporaryStateDown);
	   for (int k = 0; k <=  TmpDestinationHilbertSpace->LzMax; ++k)
	     {
	       this->TemporaryStateUp[k] = TmpDestinationHilbertSpace->TemporaryStateUp[k];
	       this->TemporaryStateDown[k] = TmpDestinationHilbertSpace->TemporaryStateDown[k];
	     }
	   for (int k = TmpDestinationHilbertSpace->LzMax + 1; k <=  this->LzMax; ++k)
	     {
	       this->TemporaryStateUp[k] = 0x0ul;
	       this->TemporaryStateDown[k] = 0x0ul;
	     }	   
	   for (int k = 0; k <=  TmpHilbertSpace->LzMax; ++k)
	     {
	       this->TemporaryStateUp[k] += TmpHilbertSpace->TemporaryStateUp[k];
	       this->TemporaryStateDown[k] += TmpHilbertSpace->TemporaryStateDown[k];
	     }

	   int TmpPos = this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
	   if (TmpPos != this->HilbertSpaceDimension)
	     {
	       double TmpFactorial = 0.0;	      
	       for (int k = 0; k <= this->LzMax; ++k)
		 {
		   TmpFactorial += LogFactorials[this->TemporaryStateUp[k]];
		   TmpFactorial += LogFactorials[this->TemporaryStateDown[k]];
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
