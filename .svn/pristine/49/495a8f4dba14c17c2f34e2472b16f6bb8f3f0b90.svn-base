////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of bosons on sphere                        //
//                                                                            //
//                        last modification : 05/07/2002                      //
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
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/StringTools.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"

#include "GeneralTools/List.h"
#include "GeneralTools/OrderedList.h"
#include "GeneralTools/SmallIntegerArray.h"
#include "GeneralTools/Permutations.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

// testing flag - switches debut output
// #define TESTING

// default constructor
//

BosonOnSphereWithSpin::BosonOnSphereWithSpin ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson
// totalSpin = spin projection
// memory = amount of memory granted for precalculations
//
BosonOnSphereWithSpin::BosonOnSphereWithSpin (int nbrBosons, int totalLz, int lzMax, int totalSpin, unsigned long memory)
{
  cout << "BosonOnSphereWithSpin"<<endl;
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->TotalSpin = totalSpin;
  this->NbrBosonsUp = (this->NbrBosons + this->TotalSpin) >> 1;
  this->NbrBosonsDown = (this->NbrBosons - this->TotalSpin) >> 1;

  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, this->TotalLz, this->TotalSpin);
  cout << "dim = " << this->HilbertSpaceDimension << endl;
  this->TemporaryState = new unsigned [this->NbrLzValue];
  this->TemporaryStateUp = new unsigned long [this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned [this->NbrLzValue];
  this->TemporaryMonomials = new unsigned long [this->NbrBosons];
  this->Flag.Initialize();
  this->TargetSpace = this;
  this->StateDescription = new unsigned* [this->HilbertSpaceDimension];
  this->StateLzMaxUp = new unsigned [this->HilbertSpaceDimension];
  this->StateLzMaxDown = new unsigned [this->HilbertSpaceDimension];
  int TmpLzMax = this->LzMax;
  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  //int TmpDim = this->GenerateStates(this->NbrBosonsUp, this->NbrBosonsDown, TmpLzMax, TmpLzMax, this->ShiftedTotalLz, 0);
  long TmpDim = GenerateStates(this->NbrBosonsUp, this->NbrBosonsDown, this->LzMax, this->TotalLz);
  
  if (TmpDim != this->HilbertSpaceDimension)
    cout << "Count inconsistent: "<<TmpDim<<" vs " << this->HilbertSpaceDimension<<endl;
  if (this->HilbertSpaceDimension > 0)
    {
      this->GenerateLookUpTable(memory);
      this->CoherenceFactors=new double[NbrBosons*NbrBosons+1];
      for (int i=0; i<NbrBosons*NbrBosons+1; ++i)
	this->CoherenceFactors[i]=sqrt((double)i);
      this->KeptCoordinates = new int;
      (*(this->KeptCoordinates)) = -1;
      this->Minors = 0;
      this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

      //   for (int i=0; i<HilbertSpaceDimension; ++i)
      //     {
      //       PrintState(cout,i)<<endl;      
      //     }
      
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(unsigned long);  // StateDescriptionUp/Down
      UsedMemory += this->HilbertSpaceDimension * sizeof(unsigned); // StateInfo
      UsedMemory += (0x1l<<(this->LzMax+this->NbrBosonsUp)) * sizeof(unsigned long);  // LookUpTableUp 
      UsedMemory += (0x1l<<(this->LzMax+this->NbrBosonsDown)) * sizeof(unsigned long); // LookUpTableDown
      cout << "HS using ";
      PrintMemorySize(cout, UsedMemory)<<endl;
#endif
    }
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSphereWithSpin::BosonOnSphereWithSpin(const BosonOnSphereWithSpin& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->TotalSpin = bosons.TotalSpin;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->StateLzMaxUp = bosons.StateLzMaxUp;
  this->StateLzMaxDown = bosons.StateLzMaxDown;
  this->StateInfo = bosons.StateInfo;
  this->LookUpTableUp = bosons.LookUpTableUp;
  this->LookUpTableDown = bosons.LookUpTableDown;
  this->CoherenceFactors = bosons.CoherenceFactors;
  this->Flag = bosons.Flag;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned [this->NbrLzValue];
  this->TemporaryStateUp = new unsigned long [this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned [this->NbrLzValue];
  this->TemporaryMonomials = new unsigned long [this->NbrBosons];
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;

}

// destructor
//

BosonOnSphereWithSpin::~BosonOnSphereWithSpin ()
{
  delete[] this->TemporaryState;
  delete[] this->TemporaryStateUp;
  delete[] this->TemporaryStateDown;
  delete[] this->TemporaryMonomials;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->StateDescription!=NULL)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    delete[] this->StateDescription[i];
	  delete[] this->StateDescription;
	}
      delete[] this->StateDescriptionUp;
      delete[] this->StateDescriptionDown;
      delete[] this->StateInfo;
      delete[] this->LookUpTableUp;
      delete[] this->LookUpTableDown;
      delete[] this->CoherenceFactors;
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereWithSpin& BosonOnSphereWithSpin::operator = (const BosonOnSphereWithSpin& bosons)
{
  delete[] this->TemporaryState;
  delete[] this->TemporaryStateUp;
  delete[] this->TemporaryStateDown;
  delete[] this->TemporaryMonomials;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->StateDescription!=NULL)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    delete[] this->StateDescription[i];
	  delete[] this->StateDescription;
	}
      delete[] this->StateDescriptionUp;
      delete[] this->StateDescriptionDown;
      delete[] this->StateInfo;
      delete[] this->LookUpTableUp;
      delete[] this->LookUpTableDown;
      delete[] this->CoherenceFactors;
      if (this->StateLzMaxUp!=NULL)
	delete[] this->StateLzMaxUp;
      if (this->StateLzMaxDown!=NULL)
	delete[] this->StateLzMaxDown;
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->TotalSpin = bosons.TotalSpin;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->StateLzMaxUp = bosons.StateLzMaxUp;
  this->StateLzMaxDown = bosons.StateLzMaxDown;
  this->StateInfo = bosons.StateInfo;
  this->LookUpTableUp = bosons.LookUpTableUp;
  this->LookUpTableDown = bosons.LookUpTableDown;
  this->CoherenceFactors = bosons.CoherenceFactors;
  this->Flag = bosons.Flag;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned [this->NbrLzValue];
  this->TemporaryStateUp = new unsigned long [this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned [this->NbrLzValue];
  this->TemporaryMonomials = new unsigned long [this->NbrBosons];
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereWithSpin::Clone()
{
  return new BosonOnSphereWithSpin(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereWithSpin::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereWithSpin::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnSphereWithSpin::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
  this->TargetSpace = (BosonOnSphereWithSpin*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int BosonOnSphereWithSpin::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}


//Get the fermionic state description
//stateDescriptionUp = unsigned long where state description of up spins will be stored
//stateDescriptionDown = unsigned long where state description of down spins will be stored
void BosonOnSphereWithSpin::GetBosonicDescription(int state, unsigned long * & stateDescriptionUp, unsigned long * & stateDescriptionDown)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[state], this->StateDescriptionDown[state], this->StateInfo[state],
		       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
  unsigned* TmpState = TemporaryState; // this->StateDescription[state];
  for (int i=0; i <= LzMax; ++i)
    {
      stateDescriptionUp[i] = (TmpState[i]>>16);
      stateDescriptionDown[i] = (TmpState[i]&0xffff);
    }
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereWithSpin::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m1_d a^+_m2_d a_n1_d a_n2_d operator to a given state (with m1+m2=n1+n2, only spin down)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// apply a^+_m1_u a^+_m2_u a_n1_u a_n2_u operator to a given state (with m1+m2=n1+n2, only spin up)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// apply a^+_d_m1 a^+_u_m2 a_d_n1 a_u_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AduAdu (int m1, int m2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  this->TemporaryState[m2] += 0x10000;
  //coefficient = (this->TemporaryState[m2] >> 16);
  int CoherenceIndex = (this->TemporaryState[m2] >> 16);
  this->TemporaryState[m1] += 0x10000;
  //coefficient *= (this->TemporaryState[m1] >> 16);
  CoherenceIndex *= (this->TemporaryState[m1] >> 16);
  coefficient = this->CoherenceFactors[CoherenceIndex];
//   int NewLzSzMax = this->LzMax;
//   while (this->TemporaryState[NewLzSzMax] == 0)
//     --NewLzSzMax;
//   if (this->TemporaryState[NewLzSzMax]&0xffff0000)
//     {
//       NewLzSzMax<<=1;
//       NewLzSzMax+=1;
//     }
//   else NewLzSzMax<<=1;
  return this->TargetSpace->FindStateIndex(this->TemporaryState);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AdAd method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AddAdd (int m1, int m2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m2];
  //coefficient = (this->TemporaryState[m2] & 0xffff);
  int CoherenceIndex = (this->TemporaryState[m2] & 0xffff);
  ++this->TemporaryState[m1];
  //coefficient *= (this->TemporaryState[m1] & 0xffff);
  CoherenceIndex *= (this->TemporaryState[m1] & 0xffff);
  coefficient = this->CoherenceFactors[CoherenceIndex];

  return this->TargetSpace->FindStateIndex(this->TemporaryState);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAd method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::AduAdd (int m1, int m2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m2];
  //coefficient = (this->TemporaryState[m2] & 0xffff);
  int CoherenceIndex = (this->TemporaryState[m2] & 0xffff);
  this->TemporaryState[m1] += 0x10000;
  //coefficient *= (this->TemporaryState[m1] >> 16);
  CoherenceIndex *= (this->TemporaryState[m1] >> 16);
  coefficient = this->CoherenceFactors[CoherenceIndex];
  return this->TargetSpace->FindStateIndex(this->TemporaryState);
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 
double BosonOnSphereWithSpin::AuAu (int index, int n1, int n2)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       ProdATemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
  if ((n1 > CurrentLzMaxUp) || (n2 > CurrentLzMaxUp) || ((ProdATemporaryState[n1] >> 16) == 0)
      || ((ProdATemporaryState[n2] >> 16) == 0) || ((n1 == n2) && ((ProdATemporaryState[n1] >> 16) == 1)))
    return 0.0;
  //double Coefficient = (this->ProdATemporaryState[n2] >> 16);
  int CoherenceIndex = (this->ProdATemporaryState[n2] >> 16);
  this->ProdATemporaryState[n2] -= 0x10000;
  //Coefficient *= (this->ProdATemporaryState[n1] >> 16);
  CoherenceIndex *= (this->ProdATemporaryState[n1] >> 16);
  this->ProdATemporaryState[n1] -= 0x10000;

  return this->CoherenceFactors[CoherenceIndex];
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double BosonOnSphereWithSpin::AdAd (int index, int n1, int n2)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       ProdATemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);

  if ((n1 > CurrentLzMaxDown) || (n2 > CurrentLzMaxDown) || ((ProdATemporaryState[n1] & 0xffff) == 0)
      || ((ProdATemporaryState[n2] & 0xffff) == 0) || ((n1 == n2) && ((ProdATemporaryState[n1] & 0xffff) == 1)))
    return 0.0;
  //double Coefficient = (this->ProdATemporaryState[n2] & 0xffff);
  int CoherenceIndex = (this->ProdATemporaryState[n2] & 0xffff);
  --this->ProdATemporaryState[n2];
  //Coefficient *= (this->ProdATemporaryState[n1] & 0xffff);
  CoherenceIndex *= (this->ProdATemporaryState[n1] & 0xffff);
  --this->ProdATemporaryState[n1];
  return this->CoherenceFactors[CoherenceIndex];
}

// apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double BosonOnSphereWithSpin::AuAd (int index, int n1, int n2)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       ProdATemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
  if ((n1 > CurrentLzMaxUp) || (n2 > CurrentLzMaxDown) || ((ProdATemporaryState[n1] >> 16) == 0)
      || ((ProdATemporaryState[n2] & 0xffff) == 0))
    return 0.0;
  //double Coefficient = (this->ProdATemporaryState[n2] & 0xffff);
  int CoherenceIndex = (this->ProdATemporaryState[n2] & 0xffff);
  --this->ProdATemporaryState[n2];
  //Coefficient *= (this->ProdATemporaryState[n1] >> 16);
  CoherenceIndex *= (this->ProdATemporaryState[n1] >> 16);
  this->ProdATemporaryState[n1] -= 0x10000;
  return this->CoherenceFactors[CoherenceIndex];
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereWithSpin::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       ProdATemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
  
  int TmpCoefficient = 1;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      unsigned& Tmp = this->ProdATemporaryState[n[i]];
      if (spinIndices[i] == 0)
	{
	  if ((n[i]>CurrentLzMaxDown)||((Tmp & 0xffff) == 0))
	    return 0.0;
	  TmpCoefficient *= (Tmp & 0xffff);
	  --Tmp;
	}
      else
	{
	  if ((n[i]>CurrentLzMaxUp)||((Tmp >> 16) == 0))
	    return 0.0;
	  TmpCoefficient *= (Tmp >> 16);
	  Tmp -= 0x10000;
	}
    }
  return sqrt((double) TmpCoefficient);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpin::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  int i = 0;
  for (; i < this->NbrLzValue; ++i)
    {
      this->TemporaryState[i] = this->ProdATemporaryState[i];
    }
  int TmpCoefficient = 1;
  for (i = 0; i < nbrIndices; ++i)
    if (spinIndices[i] == 0)
      {
	++this->TemporaryState[m[i]];
	TmpCoefficient *= (this->TemporaryState[m[i]] &0xffff);
      }
    else
      {
	this->TemporaryState[m[i]] += 0x10000;
	TmpCoefficient *= (this->TemporaryState[m[i]] >> 16);
      }
  coefficient = sqrt((double) TmpCoefficient);
  return this->TargetSpace->FindStateIndex(this->TemporaryState);
}


// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereWithSpin::AduAu (int index, int m)
{
  unsigned Info = this->StateInfo[index];
  int CurrentLzMaxUp = ((Info>>10)&0x3ff) - NbrBosonsUp + (NbrBosonsUp!=0);
  if (CurrentLzMaxUp < m)
    return 0.0;
  else
    {
      int CurrentLzMaxDown;
      this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], Info,
			   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
      return (double) ((this->TemporaryState[m] >> 16));
    }
}

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereWithSpin::AddAd (int index, int m)
{
  unsigned Info = this->StateInfo[index];
  int CurrentLzMaxDown = (Info&0x3ff) - NbrBosonsDown + (NbrBosonsDown!=0);
  if (CurrentLzMaxDown < m)
    {
      return 0.0;
    }
  else
    {
      int CurrentLzMaxUp;
      this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], Info,
			   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
      return (double) ((this->TemporaryState[m] & 0xffff));
    }
}

// apply a^+_m_u a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnSphereWithSpin::AduAu (int index, int m, int n, double& coefficient)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);

  if ((n > CurrentLzMaxUp) || ((TemporaryState[n] >> 16) == 0)) // || ((State[n] & 0xffff) == 0)) // shift for up, mask for down
    {
      return 0.0;
    }
  //double TmpCoefficient = (this->TemporaryState[n] >> 16);
  int CoherenceIndex = (this->TemporaryState[n] >> 16);
  this->TemporaryState[n] -= 0x10000;
  
  this->TemporaryState[m] += 0x10000;
  //TmpCoefficient *= (this->TemporaryState[m] >> 16);
  CoherenceIndex *= (this->TemporaryState[m] >> 16);  
  coefficient *= this->CoherenceFactors[CoherenceIndex];
  return this->TargetSpace->FindStateIndex(this->TemporaryState);
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnSphereWithSpin::AddAd (int index, int m, int n, double& coefficient)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);

  if ((n > CurrentLzMaxDown) || ((TemporaryState[n] & 0xffff) == 0)) // || ((State[n] & 0xffff) == 0)) // shift for up, mask for down
    {
      return 0.0;
    }
  //double TmpCoefficient = (this->TemporaryState[n] & 0xffff);
  int CoherenceIndex = (this->TemporaryState[n] & 0xffff);
  --this->TemporaryState[n];
  ++this->TemporaryState[m];
  //TmpCoefficient *= (this->TemporaryState[m] & 0xffff);
  CoherenceIndex *= (this->TemporaryState[m] & 0xffff);
  coefficient *= this->CoherenceFactors[CoherenceIndex];
  return this->TargetSpace->FindStateIndex(this->TemporaryState);
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnSphereWithSpin::AduAd (int index, int m, int n, double& coefficient)
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnSphereWithSpin::AddAu (int index, int m, int n, double& coefficient)
{
  return this->TargetSpace->HilbertSpaceDimension;
}


// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int BosonOnSphereWithSpin::FindStateIndex(char* stateDescription)
{
  char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != (this->LzMax + 1))
    return -1;
  int TmpNbrParticles = 0;
  int TmpTotalLz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      int Tmp = atoi(TmpDescription[i]);
      this->TemporaryState[i] = Tmp;
      TmpTotalLz += (i * Tmp);
      TmpNbrParticles += Tmp;
      delete[] TmpDescription[i];
    }
  delete[] TmpDescription;
  if ((TmpNbrParticles != this->NbrBosons) || (TmpTotalLz != ((this->TotalLz + this->NbrBosons * this->LzMax) >> 1)))
    return -1;
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState);
}



// sort an array and reflect permutations in auxiliary array
//
// length = length of arrays
// sortArray = array to be sorted
// auxArray = auxiliary array
//
void BosonOnSphereWithSpin::ShellSortAux(unsigned length, unsigned long* sortArray, unsigned *auxArray, int *auxArray2)
{
  unsigned inc = round(length/2.0);
  unsigned long tmpI;
  unsigned tmpA;
  int tmpA2;
  while (inc > 0)
    {
      for (unsigned i = inc; i< length; ++i)
	{
	  tmpI = sortArray[i];
	  tmpA = auxArray[i];
	  tmpA2 = auxArray2[i];
	  unsigned j = i;
	  while ((j>=inc) && ( sortArray[j-inc] < tmpI) )
	    {
	      sortArray[j] = sortArray[j - inc];
	      auxArray[j] = auxArray[j - inc];
	      auxArray2[j] = auxArray2[j - inc];
	      j = j - inc;
	    }
	  sortArray[j] = tmpI;
	  auxArray[j] = tmpA;
	  auxArray2[j] = tmpA2;
	}
      inc = round(inc / 2.2);
    }
}
  

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSpin::PrintState (ostream& Str, int state)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[state], this->StateDescriptionDown[state], this->StateInfo[state],
		       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
  unsigned* TmpState = TemporaryState; // this->StateDescription[state];
  int i = 0;
  Str << state << " = " << "|";
  Str << (TmpState[i]>>16) << "u " << (TmpState[i]&0xffff)<<"d";
  for (i=1; i <= LzMax; ++i)
    Str << " | "<< (TmpState[i]>>16) << "u " << (TmpState[i]&0xffff)<<"d";
  //  for (; i <= this->LzMax; ++i)
  //    Str << " | 0u 0d";
  //  Str << "  lzmaxU = " <<((this->StateInfo[state]>>20)&0x03ffu)<< "  lzMaxD = "<< ((this->StateInfo[state]>>10)&0x03ffu)
  //      << "  nbrUp = " << (this->StateInfo[state]&0x03ffu) <<"  position = ";
  //  Str << FindStateIndex(TmpState, this->StateInfo[state]&0x03ffu);
  return Str;
}


// print a given State
//
// Str = reference on current output stream 
// myState = explicit form of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSpin::PrintState (ostream& Str, unsigned *myState)
{
  int i = 0;
  Str << (myState[i]>>16) << "u " << (myState[i]&0xffff)<<"d";
  for (i=1; i <= LzMax; ++i)
    Str << " | "<< (myState[i]>>16) << "u " << (myState[i]&0xffff)<<"d";
  //Str << " key = " << this->Keys[state] << " lzmax  = " << this->StateLzMax[state]<< " position = " << FindStateIndex(myState, Max);
//   Str << " key = " << this->Keys[state] << " lzmax position = " << this->LzMaxPosition[Max * (this->NbrBosons + 1) + TmpState[Max]]
//       << " position = " << FindStateIndex(TmpState, Max);
  return Str;
}



// generate all states corresponding to the constraints
// 
// nbrBosonsUp = number of bosons with spin up
// nbrBosonsDown = number of bosons with spin down
// lzMax = momentum maximum value for a boson in the state
// totalLz = momentum total value
// return value = position from which new states have to be stored

long BosonOnSphereWithSpin::GenerateStates(int nbrBosonsUp, int nbrBosonsDown, int lzMax, int totalLz)
{
  long Pos=0;
  
  for (int TotalLzUp=totalLz-nbrBosonsDown*lzMax; TotalLzUp<=totalLz+nbrBosonsDown*lzMax; TotalLzUp+=2)
    {
      int ShiftedTotalLzUp = (TotalLzUp + nbrBosonsUp * this->LzMax) >> 1;
      int ShiftedTotalLzDown = (totalLz-TotalLzUp + nbrBosonsDown * this->LzMax) >> 1;
      if ((ShiftedTotalLzUp>=0)&&(ShiftedTotalLzDown>=0))
	Pos=this->GenerateStatesWithConstraint(nbrBosonsUp, nbrBosonsDown, lzMax, lzMax, lzMax, ShiftedTotalLzUp, ShiftedTotalLzDown, Pos, /*level*/ 0);
    }
  return Pos;
}

// generate all states corresponding to the constraints
// 
// nbrBosonsUp = number of bosons with spin up
// nbrBosonsDown = number of bosons with spin down
// lzMax = momentum maximum value for a boson in the state
// lzMaxUp = momentum maximum value for a spin-up boson in the state
// currentLzMax = momentum maximum value for bosons that are still to be placed
// totalLz = momentum total value
// totalLzUp = momentum total value for spin-up bosons
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSpin::GenerateStatesWithConstraint(int nbrBosonsUp, int nbrBosonsDown, int lzMax, int lzMaxUp, int currentLzMax, int totalLzUp, int totalLzDown, long pos, int level)
{
#ifdef TESTING
  for (int i=0; i<level; ++i) cout << " ";
  cout << "GenerateStates(Up:"<< nbrBosonsUp<<", Down:"<< nbrBosonsDown<<", lzMax:"<<lzMax<<", lzMaxUp="<<lzMaxUp<<", currLz:"<<currentLzMax<< ", totalLzUp:" <<totalLzUp<<", totalLzDown:" <<totalLzDown<<", pos:"<< pos<<")"<<endl;
#endif

  // place up spin bosons, first, given the additional constraint on the angular momentum
  if (nbrBosonsUp>0)
    {
      if (((nbrBosonsUp * currentLzMax) < totalLzUp) || (pos == this->HilbertSpaceDimension))
	{
	  return pos;
	}
      if ((nbrBosonsUp * currentLzMax) == totalLzUp)
	{
	  long TmpPos;
	  if (nbrBosonsDown>0)
	    TmpPos = this->GenerateStatesWithConstraint(0, nbrBosonsDown, lzMax, 0, lzMax, 0, totalLzDown, pos, level+1);
	  else
	    {
	      TmpPos = pos+1;
	      this->StateDescription[pos] = new unsigned [this->LzMax + 1];
	      this->StateLzMaxDown[pos] = 0;
	      unsigned* TmpState = this->StateDescription[pos];
	      for (int i = 0; i <= this->LzMax; ++i)
		TmpState[i] = 0;
	    }
	  for (long i = pos; i < TmpPos; i++)
	    {
	      this->StateDescription[i][currentLzMax] |= nbrBosonsUp<<16;
	      this->StateLzMaxUp[i] = lzMaxUp;
	    }
	  return TmpPos;
	}
      if ((currentLzMax == 0) || (totalLzUp == 0))
	{
	  long TmpPos;
	  if (nbrBosonsDown>0)
	    TmpPos = this->GenerateStatesWithConstraint(0, nbrBosonsDown, lzMax, 0, lzMax, 0, totalLzDown, pos, level+1);
	  else
	    {
	      TmpPos = pos+1;
	      this->StateDescription[pos] = new unsigned [this->LzMax + 1];
	      this->StateLzMaxDown[pos] = 0;
	      unsigned* TmpState = this->StateDescription[pos];
	      for (int i = 0; i <= this->LzMax; ++i)
		TmpState[i] = 0;
	    }
	  for (long i = pos; i < TmpPos; i++)
	    {
	      this->StateDescription[i][0] |= nbrBosonsUp<<16;
	      this->StateLzMaxUp[i] = lzMaxUp;
	    }
	  return TmpPos;
	}

      int TmpTotalLzUp = totalLzUp / currentLzMax;
      int TmpNbrBosonsUp = nbrBosonsUp - TmpTotalLzUp;
      TmpTotalLzUp = totalLzUp - TmpTotalLzUp * currentLzMax;
      int ReducedCurrentLzMax = currentLzMax - 1;
      long TmpPos = pos;
      while (TmpNbrBosonsUp < nbrBosonsUp)
	{
	  TmpPos = this->GenerateStatesWithConstraint(TmpNbrBosonsUp, nbrBosonsDown, lzMax, lzMaxUp, ReducedCurrentLzMax, TmpTotalLzUp, totalLzDown, pos, level+1);
	  for (long i = pos; i < TmpPos; i++)
	    this->StateDescription[i][currentLzMax] |= (nbrBosonsUp - TmpNbrBosonsUp)<<16;
	  ++TmpNbrBosonsUp;
	  pos = TmpPos;
	  TmpTotalLzUp += currentLzMax;
	}
      if (lzMaxUp == currentLzMax)
	return this->GenerateStatesWithConstraint(nbrBosonsUp, nbrBosonsDown, lzMax, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLzUp, totalLzDown, pos, level+1);
      else
	return this->GenerateStatesWithConstraint(nbrBosonsUp, nbrBosonsDown, lzMax, lzMaxUp, ReducedCurrentLzMax, totalLzUp, totalLzDown, pos, level+1);
    }
  else // place down spins
    {
      if ((nbrBosonsDown == 0) || ((nbrBosonsDown * currentLzMax) < totalLzDown) || (pos == this->HilbertSpaceDimension))
	{
	  return pos;
	}
      if ((nbrBosonsDown * currentLzMax) == totalLzDown)
	{
	  this->StateDescription[pos] = new unsigned [this->LzMax + 1];
	  unsigned* TmpState = this->StateDescription[pos];
	  for (int i = 0; i <= this->LzMax; ++i)
	    TmpState[i] = 0;
	  TmpState[currentLzMax] = nbrBosonsDown;
	  this->StateLzMaxDown[pos] = lzMax;
	  this->StateLzMaxUp[pos] = 0;
	  return pos + 1;
	}
      if ((currentLzMax == 0) || (totalLzDown == 0))
	{
	  this->StateDescription[pos] = new unsigned [this->LzMax + 1];
	  unsigned* TmpState = this->StateDescription[pos];
	  for (int i = 1; i <= this->LzMax; ++i)
	    TmpState[i] = 0;
	  TmpState[0] = nbrBosonsDown;
	  this->StateLzMaxDown[pos] = lzMax;
	  this->StateLzMaxUp[pos] = 0;
	  return pos + 1;
	}
      int TmpTotalLzDown = totalLzDown / currentLzMax;
      int TmpNbrBosonsDown = nbrBosonsDown - TmpTotalLzDown;
      TmpTotalLzDown = totalLzDown - TmpTotalLzDown * currentLzMax;
      int ReducedCurrentLzMax = currentLzMax - 1;
      long TmpPos = pos;
      while (TmpNbrBosonsDown < nbrBosonsDown)
	{
	  TmpPos = this->GenerateStatesWithConstraint(0, TmpNbrBosonsDown, lzMax, 0, ReducedCurrentLzMax, 0, TmpTotalLzDown, pos, level+1);
	  for (long i = pos; i < TmpPos; i++)
	    this->StateDescription[i][currentLzMax] = nbrBosonsDown - TmpNbrBosonsDown;
	  ++TmpNbrBosonsDown;
	  pos = TmpPos;
	  TmpTotalLzDown += currentLzMax;
	}
      if (lzMax == currentLzMax)
	return this->GenerateStatesWithConstraint(0, nbrBosonsDown, ReducedCurrentLzMax, 0, ReducedCurrentLzMax, 0, totalLzDown, pos, level+1);
      else
	return this->GenerateStatesWithConstraint(0, nbrBosonsDown, lzMax, 0, ReducedCurrentLzMax, 0, totalLzDown, pos, level+1);
    }
}


/* conventional Generate States, no tensored index
// generate all states corresponding to the constraints
// 
// nbrBosonsUp = number of bosons with spin up
// nbrBosonsDown = number of bosons with spin down
// lzMax = momentum maximum value for a boson in the state
// currentLzMax = momentum maximum value for bosons that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int BosonOnSphereWithSpin::GenerateStates(int nbrBosonsUp, int nbrBosonsDown, int lzMax, int currentLzMax, int totalLz, int pos, int level)
{
#ifdef TESTING
  for (int i=0; i<level; ++i) cout << " ";
  cout << "GenerateStates(Up:"<< nbrBosonsUp<<", Down:"<< nbrBosonsDown<<", lzMax:"<<lzMax<<", currLz:"<<currentLzMax<< ", totalLz:" <<totalLz<<", pos:"<< pos<<")"<<endl;
#endif
  if ((nbrBosonsUp < 0) || (nbrBosonsDown < 0) || (totalLz < 0) || (currentLzMax < 0) || (((nbrBosonsUp + nbrBosonsDown) * currentLzMax) < totalLz))
    return pos;
  if ((nbrBosonsUp + nbrBosonsDown) == 0)
    {
      if (totalLz == 0)
	{	  
	  this->StateDescription[pos] = new unsigned [lzMax + 1];
	  unsigned* TmpState = this->StateDescription[pos];
	  for (int i = 0; i <= lzMax; ++i)
	    TmpState[i] = 0;
	  this->StateLzMaxUp[pos] = 0;
	  this->StateLzMaxDown[pos] = 0;
	  return pos + 1;
	}
      else return pos;
    }
  if (((nbrBosonsUp + nbrBosonsDown) * currentLzMax) == totalLz)
    {
      this->StateDescription[pos] = new unsigned [lzMax + 1];
      unsigned* TmpState = this->StateDescription[pos];
      for (int i = 0; i <= lzMax; ++i)
	TmpState[i] = 0;
      TmpState[currentLzMax] = (nbrBosonsUp << 16) | nbrBosonsDown;
      if (nbrBosonsUp>0)
	this->StateLzMaxUp[pos] = currentLzMax;
      else
	this->StateLzMaxUp[pos] = 0;
      if (nbrBosonsDown>0)
	this->StateLzMaxDown[pos] = currentLzMax;
      else
	this->StateLzMaxDown[pos] = 0;
      return pos + 1;
    }
  if ((nbrBosonsDown + nbrBosonsUp) == 1)
    if (lzMax >= totalLz)
      {
	this->StateDescription[pos] = new unsigned [lzMax + 1];
	unsigned* TmpState = this->StateDescription[pos];
	for (int i = 0; i <= lzMax; ++i)
	  TmpState[i] = 0;
	if (nbrBosonsUp == 1)
	  {
	    TmpState[totalLz] = 1 << 16; // place spin up
	    this->StateLzMaxUp[pos] = totalLz;
	    this->StateLzMaxDown[pos] = 0;
	  }
	else
	  {
	    TmpState[totalLz] = 1; // place spin down
	    this->StateLzMaxUp[pos] = 0;
	    this->StateLzMaxDown[pos] = totalLz;
	  }
	return pos + 1;	
      }
    else return pos;
  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = pos;
  for (int UpToPlace=nbrBosonsUp; UpToPlace>=0; --UpToPlace)
    for (int DownToPlace=nbrBosonsDown; DownToPlace>=0; --DownToPlace)    
      {
	TmpPos = this->GenerateStates(nbrBosonsUp-UpToPlace, nbrBosonsDown-DownToPlace, lzMax, ReducedCurrentLzMax, totalLz-(UpToPlace+DownToPlace)*currentLzMax, pos, level+1);
	for (int i = pos; i < TmpPos; i++)
	  {
	    this->StateDescription[i][currentLzMax] = (UpToPlace << 16) | DownToPlace;
	    if (UpToPlace>0)
	      this->StateLzMaxUp[i] = currentLzMax;
	    if (DownToPlace>0)
	      this->StateLzMaxDown[i] = currentLzMax;
	  }	
	pos = TmpPos;
      }
  return pos;
}

*/

/*
int BosonOnSphereWithSpin::GenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int pos)
{
  // cout << "GenerateStates(n:"<< nbrBosons<<", currLz:"<<currentLzMax<< ", totalLz:" <<totalLz<<", pos:"<< pos<<")"<<endl;
  int CurrentMomentum = currentLzMax%this->NbrLzValue;
  int CurrentSpinUp = currentLzMax/this->NbrLzValue;
  if ((nbrBosons < 0) || (totalLz < 0) || (currentLzMax < 0) || ((CurrentSpinUp==0)&&((nbrBosons * CurrentMomentum) < totalLz)))
    return pos;
  if (nbrBosons == 0)
    {
      if (totalLz == 0)
	{
	  this->StateDescription[pos] = new unsigned [lzMax + 1];
	  unsigned* TmpState = this->StateDescription[pos];
	  for (int i = 0; i <= lzMax; ++i)
	    TmpState[i] = 0;
	  this->StateNbrUp[pos] = 0;
	  this->StateLzMaxUp[pos] = 0;
	  this->StateLzMaxDown[pos] = 0;
	  return pos + 1;
	}
      else return pos;
    }
  if ((CurrentSpinUp==0)&&(nbrBosons == 1))
    {
      if (currentLzMax >= totalLz)
	{
	  this->StateDescription[pos] = new unsigned [lzMax + 1];
	  unsigned* TmpState = this->StateDescription[pos];
	  for (int i = 0; i <= lzMax; ++i)
	    TmpState[i] = 0;	    	
	  TmpState[totalLz] = 1; // place spin down
	  this->StateNbrUp[pos] = 0;
	  this->StateLzMaxUp[pos] = 0;
	  this->StateLzMaxDown[pos] = totalLz;
	  return pos + 1;	
	}
      else return pos;
    }
  if ((CurrentSpinUp==0)&&((nbrBosons * CurrentMomentum) == totalLz))
    {
      this->StateDescription[pos] = new unsigned [lzMax + 1];
      unsigned* TmpState = this->StateDescription[pos];
      for (int i = 0; i <= lzMax; ++i)
	TmpState[i] = 0;	    	
      TmpState[CurrentMomentum] = nbrBosons; // place spin down
      this->StateNbrUp[pos] = 0;
      this->StateLzMaxUp[pos] = 0;
      this->StateLzMaxDown[pos] = CurrentMomentum;
      return pos + 1;	
    }
  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = pos;
  int TmpNbrUp;
  for (int ToPlace=nbrBosons; ToPlace>=0; --ToPlace)
    {
      TmpPos = this->GenerateStates(nbrBosons-ToPlace, lzMax, ReducedCurrentLzMax, totalLz-ToPlace*CurrentMomentum, pos);
      TmpNbrUp = ToPlace*CurrentSpinUp;
      for (int i = pos; i < TmpPos; i++)
	{
	  this->StateDescription[i][CurrentMomentum] |= ToPlace << (16*CurrentSpinUp);
	  this->StateNbrUp[i] += TmpNbrUp;
	  if (ToPlace>0)
	    {
	      if (CurrentSpinUp>0)
		this->StateLzMaxUp[i] = CurrentMomentum;
	      else
		this->StateLzMaxDown[i] = CurrentMomentum;
	    }
	}
      pos = TmpPos;
    }
  return pos;
}*/



// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereWithSpin::GenerateLookUpTable(int memory)
{
  // re-code storage of LzMaxUp etc in StateInfo & convert all states into packed storage
  this->StateInfo = new unsigned[this->HilbertSpaceDimension];
  this->StateDescriptionUp = new unsigned long[this->HilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long[this->HilbertSpaceDimension];
  int LzMaxUp, LzMaxDown;
  for (int i=0; i<this->HilbertSpaceDimension; ++i)
    {
      LzMaxUp = this->StateLzMaxUp[i] + NbrBosonsUp;
      if (NbrBosonsUp!=0) --LzMaxUp;
      LzMaxDown = this->StateLzMaxDown[i] + NbrBosonsDown;
      if (NbrBosonsUp!=NbrBosons) --LzMaxDown;

      this->StateInfo[i]=(LzMaxUp<<10) | LzMaxDown;
      this->BosonToFermion(this->StateDescriptionUp[i], this->StateDescriptionDown[i], LzMaxUp, LzMaxDown,
			   this->StateDescription[i]);
#ifdef TESTING
      cout << i << " ";
      this->PrintState(cout, StateDescription[i]) <<" "<<this->StateDescriptionUp[i]<<" "<< this->StateDescriptionDown[i] << endl;
#endif
      this->GetMonomial(i,TemporaryMonomials);
      int TotalLzUp=0;
      for(int k=0; k<this->NbrBosonsUp; ++k) 
	TotalLzUp += TemporaryMonomials[k];
      this->StateInfo[i]|=(TotalLzUp<<20);
    }
  
  // size of lookup tables
  unsigned long LookUpTableSizeUp = 0x1ul<<(this->LzMax+this->NbrBosonsUp);
  unsigned long LookUpTableSizeDown = 0x1ul<<(this->LzMax+this->NbrBosonsDown);
  // look-up table for up/down spins
  this->LookUpTableUp = new unsigned long[LookUpTableSizeUp];
  this->LookUpTableDown = new unsigned long[LookUpTableSizeDown];

  unsigned long LastUpSpins=this->StateDescriptionUp[0];
  this->LookUpTableUp[LastUpSpins]=0;
  this->LookUpTableDown[this->StateDescriptionDown[0]]=0;
  unsigned SectorCount=1;
  for (int i=1; i<this->HilbertSpaceDimension; ++i)
    {
      if (LastUpSpins!=this->StateDescriptionUp[i])
	{
	  LastUpSpins=this->StateDescriptionUp[i];
	  this->LookUpTableUp[LastUpSpins]=i;
	  this->LookUpTableDown[this->StateDescriptionDown[i]]=0;
	  SectorCount=1;
	}
      else
	{
	  this->LookUpTableDown[this->StateDescriptionDown[i]]=SectorCount;
	  ++SectorCount;
	}
    }
  
  for (int i=0; i<this->HilbertSpaceDimension; ++i)
    delete [] this->StateDescription[i];
  delete [] this->StateDescription;
  this->StateDescription=NULL;
  delete [] this->StateLzMaxUp;
  this->StateLzMaxUp=NULL;
  delete [] this->StateLzMaxDown;
  this->StateLzMaxDown=NULL;

#ifdef TESTING
  unsigned Info;
  int Index;
  for (int i=0; i<this->HilbertSpaceDimension; ++i)
    {
      Info = this->StateInfo[i];
      LzMaxUp = ((Info>>10)&0x3ff); // +CurrentStateNbrUp-(CurrentStateNbrUp!=0);
      LzMaxDown = (Info&0x3ff); // + NbrBosons - CurrentStateNbrUp - (CurrentStateNbrUp!=NbrBosons);
      Index = this->FindStateIndex(this->StateDescriptionUp[i], this->StateDescriptionDown[i]);
      if (Index != i)
	cout << "Error in state "<<i<<" base: " << this->LookUpTableUp[this->StateDescriptionUp[i]]
	     << " Offset: "<<this->LookUpTableDown[this->StateDescriptionDown[i]] << endl;
      else
	cout << "State "<<i<<" found"<<endl;
    }
#endif

}


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// totalSpin = twice the total spin value
// return value = Hilbert space dimension

long BosonOnSphereWithSpin::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalSpin)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, 
						    (totalLz + lzMax * nbrBosons) >> 1, (totalSpin + nbrBosons) >> 1);
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// totalSpin = number of particles with spin up
// return value = Hilbert space dimension      

long BosonOnSphereWithSpin::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalSpin)
{
  if ((nbrBosons < 0) || (totalLz < 0)  || (totalSpin < 0) || (totalSpin > nbrBosons) || ((lzMax * nbrBosons) < totalLz))
    return 0l;
    
  if (nbrBosons == 1) 
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }
  if (totalLz == 0)
    return 1l;

  unsigned long Tmp = 0l;  
  for (int i = totalSpin; i >= 0; --i)
    for (int j = (nbrBosons - totalSpin); j >= 0; --j)
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - (i + j), lzMax - 1, totalLz - (lzMax * (i + j)), totalSpin - i);
  return Tmp;
}




// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex BosonOnSphereWithSpin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
{
  return this->EvaluateWaveFunction(state, position, basis, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// return value = wave function evaluated at the given location

Complex BosonOnSphereWithSpin::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
							      int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, this->HilbertSpaceDimension);
}



// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex BosonOnSphereWithSpin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
					     int firstComponent, int nbrComponent)
{
  int IncMaxNbrBosons = (NbrBosonsUp>NbrBosonsDown? NbrBosonsUp : NbrBosonsDown) + 1;
  Complex Value(0.0, 0.0), TmpValue(0.0, 0.0);
  Complex Tmp;
  ComplexMatrix PermUp(NbrBosonsUp, NbrBosonsUp);
  ComplexMatrix FunctionsUp(this->LzMax + 1, NbrBosonsUp);
  ComplexMatrix PermDown(NbrBosonsDown, NbrBosonsDown);  
  ComplexMatrix FunctionsDown(this->LzMax + 1, NbrBosonsDown);
  RealVector TmpCoordinates(2);
  int* Indices = new int [this->NbrBosons];
  int Pos;
  int Lz;
  for (int j = 0; j < NbrBosonsUp; ++j)
    {
      TmpCoordinates[0] = position[j << 1];
      TmpCoordinates[1] = position[1 + (j << 1)];
      for (int i = 0; i <= this->LzMax; ++i)
	{	  
	  basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	  FunctionsUp[j].Re(i) = Tmp.Re;
	  FunctionsUp[j].Im(i) = Tmp.Im;
	}
    }
  for (int j = NbrBosonsUp; j < this->NbrBosons; ++j)
    {
      TmpCoordinates[0] = position[j << 1];
      TmpCoordinates[1] = position[1 + (j << 1)];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	  FunctionsDown[j-NbrBosonsUp].Re(i) = Tmp.Re;
	  FunctionsDown[j-NbrBosonsUp].Im(i) = Tmp.Im;
	}
    }
  double* Factors = new double [IncMaxNbrBosons];  
  Factors[0] = 1.0;
  Factors[1] = 1.0;
  for (int i = 2; i < IncMaxNbrBosons; ++i)
    Factors[i] = Factors[i - 1] / sqrt((double) i);
  double TmpFactor, TmpComponent;
  int TmpStateDescription;
  int LastComponent = firstComponent + nbrComponent;  
  int* ChangeBitSignUp;
  int* ChangeBitUp;
  PermUp.EvaluateFastPermanentPrecalculationArray(ChangeBitUp, ChangeBitSignUp);
  int* ChangeBitSignDown;
  int* ChangeBitDown;
  int CurrentLzMaxUp, CurrentLzMaxDown;
  PermDown.EvaluateFastPermanentPrecalculationArray(ChangeBitDown, ChangeBitSignDown);  
  for (int k = firstComponent; k < LastComponent; ++k)
    {
      TmpComponent = state[k];
      if (TmpComponent!=0.0)
	{
	  TmpFactor = Factors[NbrBosonsUp];
	  Pos = 0;
	  Lz = 0;
	  this->FermionToBoson(this->StateDescriptionUp[k], this->StateDescriptionDown[k], this->StateInfo[k],
			       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
	  while (Pos < NbrBosonsUp)
	    {	      
	      TmpStateDescription = (TemporaryState[Lz])>>16;
	      if (TmpStateDescription != 0)
		{
		  TmpFactor *= Factors[TmpStateDescription];
		  for (int j = 0; j < TmpStateDescription; ++j)
		    {
		      Indices[Pos] = Lz;
		      ++Pos;
		    }
		}
	      ++Lz;
	    }
	  for (int i = 0; i < NbrBosonsUp; ++i)
	    {
	      ComplexVector& TmpColum2 = FunctionsUp[i];	  
	      for (int j = 0; j < NbrBosonsUp; ++j)
		{
		  PermUp[i].Re(j) = TmpColum2.Re(Indices[j]);
		  PermUp[i].Im(j) = TmpColum2.Im(Indices[j]);
		}
	    }
	  TmpValue = PermUp.FastPermanent(ChangeBitUp, ChangeBitSignUp) * TmpFactor;
	  Pos = 0;
	  Lz = 0;
	  TmpFactor = Factors[NbrBosonsDown];
	  while (Pos < NbrBosonsDown)
	    {
	      TmpStateDescription = (TemporaryState[Lz])&0xffff;
	      if (TmpStateDescription != 0)
		{
		  TmpFactor *= Factors[TmpStateDescription];
		  for (int j = 0; j < TmpStateDescription; ++j)
		    {
		      Indices[Pos] = Lz;
		      ++Pos;
		    }
		}
	      ++Lz;
	    }
	  for (int i = 0; i < NbrBosonsDown; ++i)
	    {
	      ComplexVector& TmpColum2 = FunctionsDown[i];	  
	      for (int j = 0; j < NbrBosonsDown; ++j)
		{
		  PermDown[i].Re(j) = TmpColum2.Re(Indices[j]);
		  PermDown[i].Im(j) = TmpColum2.Im(Indices[j]);
		}
	    }
	  TmpValue *= PermDown.FastPermanent(ChangeBitDown, ChangeBitSignDown) * TmpFactor;
	  Value+=TmpComponent*TmpValue;
	}
    }
  delete[] ChangeBitSignUp;
  delete[] ChangeBitUp;
  delete[] ChangeBitSignDown;
  delete[] ChangeBitDown;
  delete[] Factors;
  delete[] Indices;
  return Value;
}

// evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex BosonOnSphereWithSpin::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
							      int nextCoordinates, int firstComponent, int nbrComponent)
{
  cout << "Attention, TimeCoherence not fully implemented, yet!"<<endl;
  double* Factors = new double [this->NbrBosons + 1];
  Factors[0] = 1.0;
  Factors[1] = 1.0;
  for (int i = 2; i <= this->NbrBosons; ++i)
    Factors[i] = Factors[i - 1] / sqrt((double) i);
  double TmpFactor;
  Complex Value;
  Complex Tmp;
  Complex TmpPerm;
  int* Indices = new int [this->NbrBosons];
  int Pos;
  int Lz;
  int TmpStateDescription;
  int CurrentLzMaxUp, CurrentLzMaxDown;
  int LastComponent = firstComponent + nbrComponent;
  if ((*(this->KeptCoordinates)) == -1)
    {
      ComplexMatrix Perm(this->NbrBosons, this->NbrBosons);
      ComplexMatrix Functions(this->LzMax + 1, this->NbrBosons);
      RealVector TmpCoordinates(2);
      for (int j = 0; j < this->NbrBosons; ++j)
	{
	  TmpCoordinates[0] = position[j << 1];
	  TmpCoordinates[1] = position[1 + (j << 1)];
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	      Functions[j].Re(i) = Tmp.Re;
	      Functions[j].Im(i) = Tmp.Im;
	    }
	}
      int* ChangeBitSign;
      int* ChangeBit;
      Perm.EvaluateFastPermanentPrecalculationArray(ChangeBit, ChangeBitSign, true);
      for (int k = firstComponent; k < LastComponent; ++k)
	{
	  Pos = 0;
	  Lz = 0;
	  TmpFactor = state[k] * Factors[this->NbrBosons];
	  this->FermionToBoson(this->StateDescriptionUp[k], this->StateDescriptionDown[k], this->StateInfo[k],
			       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
	  while (Pos < this->NbrBosons)
	    {
	      TmpStateDescription = this->TemporaryState[Lz];
	      if (TmpStateDescription != 0)
		{
		  TmpFactor *= Factors[TmpStateDescription];
		  for (int j = 0; j < TmpStateDescription; ++j)
		    {
		      Indices[Pos] = Lz;
		      ++Pos;
		    }
		}
	      ++Lz;
	    }
	  for (int i = 0; i < this->NbrBosons; ++i)
	    {
	      ComplexVector& TmpColum2 = Functions[i];	  
	      for (int j = 0; j < this->NbrBosons; ++j)
		{
		  Perm[i].Re(j) = TmpColum2.Re(Indices[j]);
		  Perm[i].Im(j) = TmpColum2.Im(Indices[j]);
		}
	    }
	  if (this->Minors[k] == 0)
	    {
	      this->Minors[k] = new Complex [this->NbrBosons];
	    }
	  Perm.FastPermanentMinorDevelopment(ChangeBit, ChangeBitSign, nextCoordinates, this->Minors[k]);
	  TmpPerm = 0.0;
	  for (int i = 0; i < this->NbrBosons; ++i)
	    TmpPerm += this->Minors[k][i] * Complex (Perm[nextCoordinates].Re(i), 
						     Perm[nextCoordinates].Im(i));
	  Value += TmpPerm * TmpFactor;
	}
      delete[] ChangeBitSign;
      delete[] ChangeBit;
      (*(this->KeptCoordinates)) = nextCoordinates;
    }
  else
    {
      Complex* Functions = new Complex[this->LzMax + 1];
      RealVector TmpCoordinates(2);
      TmpCoordinates[0] = position[(*(this->KeptCoordinates)) << 1];
      TmpCoordinates[1] = position[1 + ((*(this->KeptCoordinates)) << 1)];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  basis.GetFunctionValue(TmpCoordinates, Functions[i], i);
	}
      for (int k = firstComponent; k < LastComponent; ++k)
	{
	  Pos = 0;
	  Lz = 0;
	  TmpFactor = Factors[this->NbrBosons] * state[k];
	  this->FermionToBoson(this->StateDescriptionUp[k], this->StateDescriptionDown[k], this->StateInfo[k],
			       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
	  while (Pos < this->NbrBosons)
	    {
	      TmpStateDescription = this->TemporaryState[Lz];
	      if (TmpStateDescription != 0)
		{
		  TmpFactor *= Factors[TmpStateDescription];
		  for (int j = 0; j < TmpStateDescription; ++j)
		    {
		      Indices[Pos] = Lz;
		      ++Pos;
		    }
		}
	      ++Lz;
	    }
	  Complex* TmpMinors = this->Minors[k];
	  TmpPerm = 0.0;
	  for (int i = 0; i < this->NbrBosons; ++i)
	    TmpPerm += TmpMinors[i] * Functions[Indices[i]];
	  Value += TmpPerm * TmpFactor;
	}
      delete[] Functions;
      (*(this->KeptCoordinates)) = -1;
    }
  delete[] Factors;
  delete[] Indices;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void BosonOnSphereWithSpin::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
  if ((timeCoherence == true) && (this->Minors == 0))
    {
      this->Minors = new Complex* [this->HilbertSpaceDimension];
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	this->Minors[i] = 0;
    }
}



// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  BosonOnSphereWithSpin::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState)
{
  RealSymmetricMatrix TmpDensityMatrix;
  return TmpDensityMatrix;	  

  /*
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrBosonSector == 0))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrBosonSector == this->NbrBosons))
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * groundState[j]);
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  
  int CurrentLzMaxUp, CurrentLzMaxDown;
  unsigned TemporaryStateNbrUp;
  
  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0;
	  for (int MinIndex = 0; MinIndex < this->HilbertSpaceDimension; ++MinIndex)
	    {
	      this->FermionToBoson(this->StateDescriptionUp[MinIndex], this->StateDescriptionDown[MinIndex], this->StateInfo[MinIndex],
				   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
	      if (this->TemporaryState[0] == (unsigned)nbrBosonSector)
		TmpValue += groundState[MinIndex] * groundState[MinIndex];
	    }
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}      
    }
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0;
	  int MinIndex = 0;
	  while ((MinIndex < this->HilbertSpaceDimension) && (this->GetStateLzMax(MinIndex) >= subsytemSize))
	    {
	      int TmpPos = 0;
	      this->FermionToBoson(this->StateDescriptionUp[MinIndex], this->StateDescriptionDown[MinIndex], this->StateInfo[MinIndex],
				   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
	      while ((TmpPos < subsytemSize) && (TemporaryState[TmpPos] == 0))
		++TmpPos;
	      if (TmpPos == subsytemSize)
		{
		  TmpValue += groundState[MinIndex] * groundState[MinIndex];
		}
	      ++MinIndex;
	    }
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  if (ShiftedLzComplementarySector < 0)
    {
      RealSymmetricMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  int MinIndex = 0;
  int MaxIndex = this->HilbertSpaceDimension - 1;
  if (nbrBosonSector == 1)
    {
      double TmpValue = 0.0;
      int TmpLzMax = this->GetStateLzMax(MinIndex);
      while ((MinIndex <= MaxIndex) && (subsytemSize <= TmpLzMax))
	{
	  this->FermionToBoson(this->StateDescriptionUp[MinIndex], this->StateDescriptionDown[MinIndex], this->StateInfo[MinIndex],
			       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
	  if (TemporaryState[ShiftedLzSector] == 1)
	    {	      
	      int TmpPos = 0;
	      int TmpNbrBosons = 0;
	      while (TmpPos < subsytemSize)
		TmpNbrBosons += TemporaryState[TmpPos++];
	      if (TmpNbrBosons == 1)
		TmpValue += groundState[MinIndex] * groundState[MinIndex];	    
	    }
	  ++MinIndex;
	  if (MinIndex <= MaxIndex)
	    TmpLzMax = this->GetStateLzMax(MinIndex);
	}
      RealSymmetricMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  if (NbrBosonsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
      BosonOnSphereWithSpin TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      double TmpValue;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpValue = groundState[MinIndex + i];
	  for (int j = i; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    TmpDensityMatrix.SetMatrixElement(i, j, TmpValue * groundState[MinIndex + j]);
	}
      return TmpDensityMatrix;
    }


  int TmpNbrBosons;
  int TmpTotalLz;
  int TmpIndex;
  BosonOnSphereWithSpin TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  long TmpNbrNonZeroElements = 0;
  int TmpComplementarySubsystemLzMax = this->GetStateLzMax(MinIndex);
  unsigned* TmpComplementarySubsystem = new unsigned[NbrLzValue];
  while ((MinIndex <= MaxIndex) && (TmpComplementarySubsystemLzMax >= subsytemSize))
    {
      this->FermionToBoson(this->StateDescriptionUp[MinIndex], this->StateDescriptionDown[MinIndex], this->StateInfo[MinIndex],
			   TmpComplementarySubsystem, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
      TmpIndex = MinIndex + 1;
      int TmpPos = TmpComplementarySubsystemLzMax;	  
      int TmpLzMax = TmpPos;
      while ((TmpIndex <= MaxIndex) && (TmpPos == TmpLzMax))
	{
	  if (TmpLzMax == this->GetStateLzMax(TmpIndex))
	    {
	      TmpPos = subsytemSize;
	      this->FermionToBoson(this->StateDescriptionUp[TmpIndex], this->StateDescriptionDown[TmpIndex], this->StateInfo[TmpIndex],
				   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);

	      while ((TmpPos <= TmpLzMax) && (this->TemporaryState[TmpPos] == TmpComplementarySubsystem[TmpPos]))
		++TmpPos;
	      if (TmpPos > TmpLzMax)
		{
		  ++TmpIndex;
		  --TmpPos;
		}
	      else
		{
		  TmpPos = -1;
		}
	    }
	  else
	    {
	      TmpPos = -1;
	    }
	}
      TmpNbrBosons = 0;
      TmpTotalLz = 0;
      TmpPos = subsytemSize;	  
      while (TmpPos <= TmpComplementarySubsystemLzMax)
	{
	  TmpNbrBosons += TmpComplementarySubsystem[TmpPos];
	  TmpTotalLz += TmpComplementarySubsystem[TmpPos] * TmpPos;
	  ++TmpPos;
	}
      if ((TmpNbrBosons == NbrBosonsComplementarySector) && (ShiftedLzComplementarySector == TmpTotalLz))
	{
	  int Pos = 0;
	  for (int i = MinIndex; i < TmpIndex; ++i)
	    {
	      this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i], this->StateInfo[i],
				   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);
	      int TmpLzMax = subsytemSize - 1;
	      while (TemporaryState[TmpLzMax] == 0) 
		--TmpLzMax;
	      TmpStatePosition[Pos] = TmpDestinationHilbertSpace.FindStateIndex(TemporaryState);
	      ++Pos;
	    }
	  int Pos2;
	  Pos = 0;
	  int Pos3;
	  double TmpValue;
	  for (int i = MinIndex; i < TmpIndex; ++i)
	    {
	      Pos2 = 0;
	      Pos3 = TmpStatePosition[Pos];
	      TmpValue = groundState[i];
	      for (int j = MinIndex; j < TmpIndex; ++j)
		{
		  if (Pos3 <=  TmpStatePosition[Pos2])
		    {
		      TmpDensityMatrix.AddToMatrixElement(Pos3, TmpStatePosition[Pos2], TmpValue * groundState[j]);
		      ++TmpNbrNonZeroElements;
		    }
		  ++Pos2;
		}
	      ++Pos;
	    }
	}
      MinIndex = TmpIndex;
      if (MinIndex <= MaxIndex)
	TmpComplementarySubsystemLzMax = this->GetStateLzMax(MinIndex);
    }
  delete[] TmpStatePosition;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
  */
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrParticleSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnSphereWithSpin::EvaluatePartialDensityMatrix (int subsytemSize, int nbrParticleSector, int lzSector,  int szSector, RealVector& groundState)
{
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySystemSize = this->LzMax + 1 - subsytemSize;
  int ComplementaryTotalSpin = this->TotalSpin - szSector;
  int ComplementaryTotalLz = ((this->TotalLz + (this->NbrBosons * this->LzMax)) 
			      - (lzSector + (nbrParticleSector * (subsytemSize - 1))));
  ComplementaryTotalLz -= (subsytemSize * ComplementaryNbrParticles) << 1;
  ComplementaryTotalLz -= ComplementaryNbrParticles * (ComplementarySystemSize - 1);
  cout << "Lz = " << this->TotalLz << " " << lzSector << " " << ComplementaryTotalLz << endl;
  cout << "Sz = " << this->TotalSpin << " " << szSector << " " << ComplementaryTotalSpin << endl;
  if ((nbrParticleSector == 0) && (lzSector == 0))
    {
      double Tmp = 0.0;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  unsigned long Count = 0ul;
	  this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i],
			       this->TemporaryStateUp, this->TemporaryStateDown);
	  for (int j = 0; j < subsytemSize; ++j)
	    {
	      Count += this->TemporaryStateUp[j];
	      Count += this->TemporaryStateDown[j];
	    }
	  if (Count == 0ul)
	    {
	      Tmp += groundState[i] * groundState[i];
	    }
	}
      if (Tmp == 0.0)
	{
	  RealSymmetricMatrix TmpDensityMatrixZero;
	  return TmpDensityMatrixZero;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix.AddToMatrixElement(0, 0, Tmp);
	  return TmpDensityMatrix;
	}
    }

  BosonOnSphereWithSpin SubsytemSpace (nbrParticleSector, lzSector, subsytemSize - 1, szSector);
  if (SubsytemSpace.GetHilbertSpaceDimension() == 0)
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
  BosonOnSphereWithSpin ComplementarySpace (ComplementaryNbrParticles, ComplementaryTotalLz, ComplementarySystemSize - 1, ComplementaryTotalSpin);
  if (ComplementarySpace.GetHilbertSpaceDimension() == 0)
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
  RealSymmetricMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


//   FQHESphereEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
//   Operation.ApplyOperation(architecture);
//   if (Operation.GetNbrNonZeroMatrixElements() > 0)	
//     return TmpDensityMatrix;
//   else
//     {
//       RealSymmetricMatrix TmpDensityMatrixZero;
//       return TmpDensityMatrixZero;
//     }

  long TmpNbrNonZeroMatrixElements = this->EvaluatePartialDensityMatrixCore(0, ComplementarySpace.GetHilbertSpaceDimension(), &ComplementarySpace, &SubsytemSpace, groundState, &TmpDensityMatrix);
  if (TmpNbrNonZeroMatrixElements > 0)
    {	
      return TmpDensityMatrix;
    }
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
// nbrNUpSector = number of spin up  that belong to the subsytem 
// nbrNDownSector = number of spin down  that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnSphereWithSpin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
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
       if ((lzSector == this->TotalLz) && (nbrNUpSector == this->NbrBosonsUp) && (nbrNDownSector == this->NbrBosonsDown))
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
  int ComplementaryTotalLz = (this->TotalLz - lzSector);
  int ComplementaryNUp = this->NbrBosonsUp - nbrNUpSector;
  int ComplementaryNDown = this->NbrBosonsDown - nbrNDownSector;
  int SubsytemTotalSz = nbrNUpSector - nbrNDownSector;
  int ComplementaryTotalSz = ComplementaryNUp - ComplementaryNDown;
  cout << "Lz = " << this->TotalLz << " " << lzSector << " " << ComplementaryTotalLz << endl;
  cout << "nup = " << this->NbrBosonsUp << " " << nbrNUpSector << " " << ComplementaryNUp << endl;
  cout << "ndown = " << this->NbrBosonsDown << " " << nbrNDownSector << " " << ComplementaryNDown << endl;
  BosonOnSphereWithSpin SubsytemSpace (nbrParticleSector, lzSector, this->LzMax, SubsytemTotalSz);
  if (SubsytemSpace.GetHilbertSpaceDimension() == 0)
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
  BosonOnSphereWithSpin ComplementarySpace (ComplementaryNbrParticles, ComplementaryTotalLz, this->LzMax, ComplementaryTotalSz);
  if (ComplementarySpace.GetHilbertSpaceDimension() == 0)
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
  RealSymmetricMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
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

// core part of the evaluation density matrix calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnSphereWithSpin::EvaluatePartialDensityMatrixCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
							      RealVector& groundState, RealSymmetricMatrix* densityMatrix)
{
  BosonOnSphereWithSpin* TmpHilbertSpace =  (BosonOnSphereWithSpin*) complementaryHilbertSpace;
  BosonOnSphereWithSpin* TmpDestinationHilbertSpace =  (BosonOnSphereWithSpin*) destinationHilbertSpace;
  int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
  int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  int OrbitalShift = TmpDestinationHilbertSpace->LzMax + 1;
  
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->StateDescriptionUp[minIndex], TmpHilbertSpace->StateDescriptionDown[minIndex],
				      TmpHilbertSpace->TemporaryStateUp, TmpHilbertSpace->TemporaryStateDown);
      for (int k = 0; k <=  TmpHilbertSpace->LzMax; ++k)
	{
	  this->TemporaryStateUp[k + OrbitalShift] = TmpHilbertSpace->TemporaryStateUp[k];
	  this->TemporaryStateDown[k + OrbitalShift] = TmpHilbertSpace->TemporaryStateDown[k];
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
	   int TmpPos = this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
	   if (TmpPos != this->HilbertSpaceDimension)
	     {
	       TmpStatePosition[Pos] = TmpPos;
	       TmpStatePosition2[Pos] = j;
	       ++Pos;
	     }
	 }
       if (Pos != 0)
 	{
 	  ++TmpNbrNonZeroElements;
 	  for (int j = 0; j < Pos; ++j)
 	    {
 	      int Pos2 = TmpStatePosition2[j];
 	      double TmpValue = groundState[TmpStatePosition[j]];
 	      for (int k = 0; k < Pos; ++k)
 		if (TmpStatePosition2[k] >= Pos2)
 		  {
 		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
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

long BosonOnSphereWithSpin::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
									      RealVector& groundState, RealSymmetricMatrix* densityMatrix)
{
  BosonOnSphereWithSpin* TmpHilbertSpace =  (BosonOnSphereWithSpin*) complementaryHilbertSpace;
  BosonOnSphereWithSpin* TmpDestinationHilbertSpace =  (BosonOnSphereWithSpin*) destinationHilbertSpace;
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
      for (int k = 0; k <=  TmpHilbertSpace->LzMax; ++k)
	{
	  this->TemporaryStateUp[k] += TmpHilbertSpace->TemporaryStateUp[k];
	  this->TemporaryStateDown[k] += TmpHilbertSpace->TemporaryStateDown[k];
	}
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

// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

RealVector BosonOnSphereWithSpin::ForgeSU2FromU1(RealVector& upState, BosonOnSphere& upStateSpace, RealVector& downState, BosonOnSphere& downStateSpace)
{
  RealVector FinalState(this->HilbertSpaceDimension, true);
  return FinalState;
}


// Project the state from the su2 space
// to the U(1) space (u1Space)
//
// state = state that needs to be projected
// u1Space = the subspace onto which the projection is carried out
RealVector BosonOnSphereWithSpin::ForgeU1FromSU2(RealVector& state, BosonOnSphere& u1Space)
{
  int Dim2=u1Space.GetHilbertSpaceDimension();
  RealVector FinalState(u1Space.GetHilbertSpaceDimension(), true);
  int Rejected=0, Index, LzMax, CurrentLzMaxUp, CurrentLzMaxDown;
  int *TmpState = new int[NbrLzValue];
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      this->FermionToBoson(this->StateDescriptionUp[j], this->StateDescriptionDown[j], this->StateInfo[j],
			   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown);

      for (int i=0; i<NbrLzValue; ++i)
	TmpState[i] = (this->TemporaryState[i]&0xffff)+(this->TemporaryState[i]>>16);
      LzMax = std::max(CurrentLzMaxUp,CurrentLzMaxDown);
      if ((Index=u1Space.FindStateIndex(TmpState, LzMax))<Dim2)
	{
	  FinalState[Index] += state[j];
	}
      else
	++Rejected;
    } //End loop over HilbertSpace

  if (Rejected>0)
    cout<<"Attention: ForgeU1 rejected "<<Rejected<<" components"<<endl; 
  FinalState /= FinalState.Norm();
  delete[] TmpState;
  return FinalState;
}

// Calculate mean value <Sx> in a given state
//
// state = given state
double BosonOnSphereWithSpin::MeanSxValue(RealVector& state)
{
  double Coefficient;
  int Dim = this->HilbertSpaceDimension;
  
  RealVector FinalState(Dim, true);
  
  for (int i = 0; i < Dim; ++i)    
    for (int j=0; j<=this->LzMax; ++j)
      {
	int Index=this->AduAd(i,j,j,Coefficient);
	if ( (Index<Dim) && (Coefficient != 0.0))
	  FinalState[Index]+=state[i]*Coefficient*0.5;
	Index=this->AddAu(i,j,j,Coefficient);
	if ( (Index<Dim) && (Coefficient != 0.0))
	  FinalState[Index]+=state[i]*Coefficient*0.5;
      }
  return (FinalState*state);
}


// Calculate mean value <Sz> in a given state
//
// state = given state
double BosonOnSphereWithSpin::MeanSzValue(RealVector& state)
{
  double Result = 0.0;
  int Dim = this->HilbertSpaceDimension;

  for (int i = 0; i < Dim; ++i)    
    Result+=state[i]*state[i]*0.5*(2*this->StateInfo[i]-this->NbrBosons);

  return Result;
}



// Compute the product of two states that belong to different Hilbert Spaces
//
// firstState = reference on one of the states whose product will be computed
// polarizedState = reference on the other state whose product will be computed
// OutputVector = reference on the vector where the result will be stored
// PolarizedSpace = pointer on the Hilbert Space whose the polarized state belong
// minIndex = first computed component
// nbrComponents = Number of computed components
// FinalSpace = pointer on the Hilbert Space whose the final state belong

//My attempt

void BosonOnSphereWithSpin::BosonicStateWithSpinTimesBosonicState( RealVector& spinfulState, RealVector& polarizedState, RealVector& outputVector, BosonOnSphereShort * polarizedSpace, int minIndex, int nbrComponents, BosonOnSphereWithSpin * finalSpace)
{
  //cout << "BosonOnSphereWithSpin::BosonicStateWithSpinTimesBosonicState("<<minIndex<<", "<< nbrComponents<<")"<<endl;
  unsigned long * FirstMonomials = new unsigned long[this->NbrBosons];
  unsigned long * FirstPolarizedMonomials = new unsigned long[this->NbrBosons];
  unsigned long * FirstPolarizedMonomialsOccupationBasis = new unsigned long[polarizedSpace->LzMax+1];  
  unsigned * FirstPolarizedMonomialsOccupationBasisInt = new unsigned[polarizedSpace->LzMax+1];

  OrderedList< SmallIntegerArray > SpinUpOccupationBasisPartitions;
  List< SmallIntegerArray > SpinDownOccupationBasisPartitions;

  unsigned long * SpinUpPartitionExponents = new unsigned long[this->NbrBosonsUp];
  unsigned long * SpinDownPartitionExponents = new unsigned long[this->NbrBosonsDown];

  unsigned long * SpinfulUpOccupationNos = new unsigned long [this->LzMax + 1];
  unsigned long * SpinfulDownOccupationNos = new unsigned long [this->LzMax + 1];

  unsigned long * upConversionFromSmallIntegerArray = new unsigned long [finalSpace->LzMax+1];
  unsigned long * downConversionFromSmallIntegerArray = new unsigned long [finalSpace->LzMax+1];
  
  unsigned long * FinalStateUpOccupationNos = new unsigned long [finalSpace->LzMax + 1];
  unsigned long * FinalStateUpExponents = new unsigned long [this->NbrBosons];
  
  unsigned long * FinalStateDownOccupationNos = new unsigned long [finalSpace->LzMax + 1];
  unsigned long * FinalStateDownExponents = new unsigned long [this->NbrBosons];

  int * EqualPowerIndexUp = new int[this->NbrBosonsUp];
  int * EqualPowerIndexDown = new int[this->NbrBosonsDown];

  int MaxIndex = minIndex + nbrComponents;
  for (long i = minIndex; i < MaxIndex; i++)
   {
    if( spinfulState[i] != 0 ) {
      
      this->GetMonomial(i, FirstMonomials);
      
      unsigned long * FirstUpMonomials = &FirstMonomials[0];
      unsigned long * FirstDownMonomials = &FirstMonomials[this->NbrBosonsUp];

      for(int lz=0; lz<= this->LzMax; lz++)
	{
	  SpinfulUpOccupationNos[lz] = 0;
	  SpinfulDownOccupationNos[lz] = 0;
	}
    
      this->ConvertPolarizedMonomialToOccupationBasis( FirstUpMonomials, SpinfulUpOccupationNos, this->NbrBosonsUp );
      this->ConvertPolarizedMonomialToOccupationBasis( FirstDownMonomials, SpinfulDownOccupationNos, this->NbrBosonsDown );
      /*
      cout << "Next spinful state to lift\n";
      for(int lz=0; lz<= this->LzMax; lz++)
	{
	  cout << SpinfulUpOccupationNos[lz] << " u " << SpinfulDownOccupationNos[lz] << " d | "; 
	}
      cout << "\n";
      */
      
      int NbrEqualPowerUp = 0;
      int NbrEqualPowerDown = 0;
      
      for (int Index = 0; Index < this->NbrBosonsUp-1; Index++)
	{
	  if(FirstUpMonomials[Index] == FirstUpMonomials[Index+1])
	    {
	      EqualPowerIndexUp[NbrEqualPowerUp] = Index;
	      NbrEqualPowerUp++;
	    }
	}
      
      for (int Index = 0 ; Index < this->NbrBosonsDown-1; Index++) 
	{
	  if(FirstDownMonomials[Index] == FirstDownMonomials[Index+1]) 
	    {
	      EqualPowerIndexDown[NbrEqualPowerDown] = Index;
	      NbrEqualPowerDown++;
	    }
	}
      
      
      for(int j = 0; j < polarizedSpace->GetHilbertSpaceDimension(); j++)
	{
	  if( polarizedState[j] != 0 ) {       
	    polarizedSpace->GetMonomial( j, FirstPolarizedMonomials );
	    /*
	    cout << "Lifting initial states \n polarized ";
	    for(int particle = 0; particle < this->NbrBosons; particle++)
	      {
		cout << FirstPolarizedMonomials[particle] << " ";
	      }
	    cout << "\tspin up ";
	    for(int particle = 0; particle < this->NbrBosonsUp; particle++)
	      {
		cout << FirstUpMonomials[particle] << " ";
	      }
	    cout << "\tspin down ";
	    for(int particle = 0; particle < this->NbrBosonsDown; particle++)
	      {
		cout << FirstDownMonomials[particle] << " ";
	      }
	    cout << "\n";
	    */
	    //Create the set of all length NbrBosonsUp subsets of FirstPolarizedMonomials

	    SpinUpOccupationBasisPartitions.DeleteList();
	    SpinDownOccupationBasisPartitions.DeleteList();

	    for( int lzvalue = 0; lzvalue <= polarizedSpace->LzMax; lzvalue++)
	      {
		FirstPolarizedMonomialsOccupationBasis[ lzvalue ] = 0x0ul;
	      }
	    
	    this->ConvertPolarizedMonomialToOccupationBasis( FirstPolarizedMonomials, FirstPolarizedMonomialsOccupationBasis, polarizedSpace->NbrBosons );
	    for( int lzvalue = 0; lzvalue <= polarizedSpace->LzMax; lzvalue++ )
	      {
		FirstPolarizedMonomialsOccupationBasisInt[ lzvalue ] = (int) FirstPolarizedMonomialsOccupationBasis[ lzvalue ];
	      }

	    //generate all possible ways of partitioning the exponents into up and down

	    int NbrSpinPartitions = AllGivenSizeSubsets( FirstPolarizedMonomialsOccupationBasisInt, this->NbrBosons, polarizedSpace->LzMax + 1, this->NbrBosonsUp, SpinUpOccupationBasisPartitions, SpinDownOccupationBasisPartitions );

	    for(int p=0; p<NbrSpinPartitions; p++) 
	      {
		cout << "\t";
		for(int lz=0;  lz <= polarizedSpace->LzMax; lz++) 
		  {
		    cout << SpinUpOccupationBasisPartitions[p].GetElement(lz) << " ";
		  } 
	       	cout << "\t";
		for(int lz=0; lz <= polarizedSpace->LzMax; lz++)
		  {
		    cout << SpinDownOccupationBasisPartitions[p].GetElement(lz) << " ";
		  }
		cout << "\n";
	      }

	    for (int spinPartition = 0; spinPartition < NbrSpinPartitions; spinPartition++) 
	      {
		//convert each occupation basis to the monomial basis so they can be multiplied
		//calculate LzMax and TotalLz for up and down spaces and convert to unsigned longs
		int SpinUpMaxLzOfPartition = polarizedSpace->LzMax;
		int SpinDownMaxLzOfPartition = polarizedSpace->LzMax;

		int SpinUpTmpLzOfPartition = polarizedSpace->LzMax;
		int SpinDownTmpLzOfPartition = polarizedSpace->LzMax;
		while(SpinUpTmpLzOfPartition > 0 && SpinUpOccupationBasisPartitions[spinPartition].GetElement( SpinUpTmpLzOfPartition ) == 0)
		  SpinUpTmpLzOfPartition--;
		while(SpinDownTmpLzOfPartition > 0 && SpinDownOccupationBasisPartitions[spinPartition].GetElement( SpinDownTmpLzOfPartition ) == 0)
		  SpinDownTmpLzOfPartition--;
  
		SpinUpMaxLzOfPartition = SpinUpTmpLzOfPartition;
		SpinDownMaxLzOfPartition = SpinDownTmpLzOfPartition;

		int SpinUpTotalLzOfPartition = 0;
		int SpinDownTotalLzOfPartition = 0;

		//unsigned long * upConversionFromSmallIntegerArray = new unsigned long [SpinUpMaxLzOfPartition+1];
		//unsigned long * downConversionFromSmallIntegerArray = new unsigned long [SpinDownMaxLzOfPartition+1];

		//cout << "Spin up partition lz occupation values\n\t";
		for(int lzvalue = 0; lzvalue <= SpinUpMaxLzOfPartition && this->NbrBosonsUp!=0; lzvalue++) 
		  {
		    upConversionFromSmallIntegerArray[ lzvalue ] = (unsigned long) SpinUpOccupationBasisPartitions[spinPartition].GetElement( lzvalue );
		    SpinUpTotalLzOfPartition += SpinUpOccupationBasisPartitions[spinPartition].GetElement( lzvalue ) * (2*lzvalue - polarizedSpace->LzMax);
		    //cout << upConversionFromSmallIntegerArray[ lzvalue ] << " ";
		  }
		//cout << "Spin down partition ";
		for(int lzvalue = 0; lzvalue <= SpinDownMaxLzOfPartition && this->NbrBosonsDown!=0; lzvalue++)
		  {
		    downConversionFromSmallIntegerArray[ lzvalue ] = (unsigned long) SpinDownOccupationBasisPartitions[spinPartition].GetElement( lzvalue );
		    SpinDownTotalLzOfPartition += SpinDownOccupationBasisPartitions[spinPartition].GetElement( lzvalue ) * (2*lzvalue - polarizedSpace->LzMax);
		    //cout << downConversionFromSmallIntegerArray[ lzvalue ] << " ";
		  }
		cout << "\n";

		if(this->NbrBosonsUp==0)
		  {
		    BosonOnSphereShort * InitialSpinDownSpace = new BosonOnSphereShort(this->NbrBosonsDown, SpinDownTotalLzOfPartition, SpinDownMaxLzOfPartition);
		    cout << "InitialSpinDownSpace initialised with NbrBosonsDown " << this->NbrBosonsDown << " TotalLz " << SpinDownTotalLzOfPartition << " LzMax " << SpinDownMaxLzOfPartition << "\n";
		    //cout << "this->GetStateLzMaxDown(" << i << ") = " << this->GetStateLzMaxDown(i) << "\n";
		    BosonOnSphereShort * FinalSpinDownSpace = new BosonOnSphereShort(this->NbrBosonsDown, this->GetTotalLzDown(i) + SpinDownTotalLzOfPartition, this->GetStateLzMaxDown(i) + SpinDownMaxLzOfPartition );
		    cout << "FinalSpinDownSpace initialised with NbrBosonsDown " << this->NbrBosonsDown << " TotalLz " << this->GetTotalLzDown(i) + SpinDownTotalLzOfPartition << " LzMax " << this->GetStateLzMaxDown(i) + SpinDownMaxLzOfPartition << "\n";

		    unsigned long * FinalDownStates = new unsigned long [ FinalSpinDownSpace->GetHilbertSpaceDimension() ];
		    long * DownWeight = new long [ FinalSpinDownSpace->GetHilbertSpaceDimension() ];

		    InitialSpinDownSpace->ConvertToMonomial( downConversionFromSmallIntegerArray, SpinDownMaxLzOfPartition, SpinDownPartitionExponents );

		    unsigned long NbrDownStates = InitialSpinDownSpace->ProductOfTwoMonomials( FirstDownMonomials, EqualPowerIndexDown, NbrEqualPowerDown, SpinDownPartitionExponents, FinalDownStates, DownWeight, FinalSpinDownSpace);

		    //	    cout << "NbrDownStates = " << NbrDownStates << "\n";

		    for(unsigned downIndex = 0; downIndex < NbrDownStates; downIndex++)
		      {
			//unsigned long * FinalStateDownOccupationNos = new unsigned long [finalSpace->LzMax + 1];
			for( int lz = 0; lz <= finalSpace->LzMax; lz++ )
			  FinalStateDownOccupationNos[lz] = 0x0ul;
			//	cout <<  "FinalSpinDownSpace->LzMax " << FinalSpinDownSpace->LzMax << "\n";
			//unsigned long * FinalStateDownExponents = new unsigned long [this->NbrBosonsDown];
			for(int p=0; p<this->NbrBosonsDown; p++)
			  FinalStateDownExponents[p] = 0x0ul;
			this->ConvertToMonomial( 0x0ul, FinalDownStates[downIndex], FinalSpinDownSpace->LzMax, FinalStateDownExponents );
			/*
			cout << "Exponents ";
			for(int p = 0; p < this->NbrBosonsDown; p++)
			  cout << FinalStateDownExponents[p] << " ";
			cout << "\n";
			*/
			this->ConvertPolarizedMonomialToOccupationBasis( FinalStateDownExponents, FinalStateDownOccupationNos, this->NbrBosonsDown );

			double geometricCorrectionFactorDown = this->GeometricCorrectionFactor( FirstMonomials, FirstPolarizedMonomials, FinalStateDownExponents, polarizedSpace->LzMax, this->LzMax );
			double occupationCorrectionFactorDown = this->OccupationCorrectionFactor(FirstPolarizedMonomialsOccupationBasis, SpinfulDownOccupationNos, FinalStateDownOccupationNos, SpinDownMaxLzOfPartition, this->GetStateLzMaxDown(i) );

			unsigned long EmptyUpState = 0x0ul;
			outputVector[ finalSpace->FindStateIndex( EmptyUpState, FinalDownStates[downIndex] ) ] += spinfulState[i]*polarizedState[j]*DownWeight[downIndex]*geometricCorrectionFactorDown*occupationCorrectionFactorDown;
		      }
		    delete InitialSpinDownSpace;
		    delete FinalSpinDownSpace;
		    delete [] FinalDownStates;
		    delete [] DownWeight;
		  }
		else if (this->NbrBosonsDown==0)
		  {
		    BosonOnSphereShort * InitialSpinUpSpace = new BosonOnSphereShort(this->NbrBosonsUp, SpinUpTotalLzOfPartition, SpinUpMaxLzOfPartition);
		    BosonOnSphereShort * FinalSpinUpSpace = new BosonOnSphereShort(this->NbrBosonsUp, this->GetTotalLzUp(i) + SpinUpTotalLzOfPartition, this->GetStateLzMaxUp(i) + SpinUpMaxLzOfPartition );
		  
		    unsigned long * FinalUpStates = new unsigned long [ FinalSpinUpSpace->GetHilbertSpaceDimension() ];
		    long * UpWeight = new long [FinalSpinUpSpace->GetHilbertSpaceDimension() ];

		    InitialSpinUpSpace->ConvertToMonomial( upConversionFromSmallIntegerArray, SpinUpMaxLzOfPartition, SpinUpPartitionExponents );

		    unsigned long NbrUpStates = InitialSpinUpSpace->ProductOfTwoMonomials( FirstUpMonomials, EqualPowerIndexUp, NbrEqualPowerUp, SpinUpPartitionExponents, FinalUpStates, UpWeight, FinalSpinUpSpace);
		    
		    //	    cout << "NbrUpStates = " << NbrUpStates << "\n";

		    for( unsigned upIndex = 0; upIndex < NbrUpStates; upIndex++ )
		      {
			
			//unsigned long * FinalStateUpOccupationNos = new unsigned long [finalSpace->LzMax + 1];
			for( int lz = 0; lz <= finalSpace->LzMax; lz++ )
			  FinalStateUpOccupationNos[lz] = 0x0ul;
			//unsigned long * FinalStateUpExponents = new unsigned long [this->NbrBosonsUp];
			for(int p=0; p < this->NbrBosonsUp; p++)
			  FinalStateUpExponents[p] = 0x0ul;
			this->ConvertToMonomial( FinalUpStates[upIndex], 0x0ul, FinalSpinUpSpace->LzMax, FinalStateUpExponents );
			/*
			cout << "Exponents ";
			for(int p = 0; p < this->NbrBosonsUp; p++)
			  cout << FinalStateUpExponents[p] << " ";
			cout << "\n";
			*/
			this->ConvertPolarizedMonomialToOccupationBasis( FinalStateUpExponents, FinalStateUpOccupationNos, this->NbrBosonsUp );

			double geometricCorrectionFactorUp = this->GeometricCorrectionFactor( FirstPolarizedMonomials, FirstMonomials, FinalStateUpExponents, polarizedSpace->LzMax, this->LzMax );
			double occupationCorrectionFactorUp = this->OccupationCorrectionFactor( FirstPolarizedMonomialsOccupationBasis, SpinfulUpOccupationNos, FinalStateUpOccupationNos, SpinUpMaxLzOfPartition, this->GetStateLzMaxUp(i) );
			
			//cout << "Weight " << UpWeight[upIndex] << "\n";

			unsigned long EmptyDownState = 0x0ul;
			outputVector[ finalSpace->FindStateIndex( FinalUpStates[upIndex], EmptyDownState ) ] += spinfulState[i]*polarizedState[j]*UpWeight[upIndex]*occupationCorrectionFactorUp*geometricCorrectionFactorUp;
			//cout << "UpWeight[" << upIndex << "] = " << UpWeight[upIndex] << "\n";
		      }
		    delete InitialSpinUpSpace;
		    delete FinalSpinUpSpace;
		    delete [] FinalUpStates;
		    delete [] UpWeight;
		  }
		else
		  {
		    BosonOnSphereShort * InitialSpinUpSpace = new BosonOnSphereShort(this->NbrBosonsUp, SpinUpTotalLzOfPartition, polarizedSpace->LzMax );
		    BosonOnSphereShort * InitialSpinDownSpace = new BosonOnSphereShort(this->NbrBosonsDown, SpinDownTotalLzOfPartition, polarizedSpace->LzMax );

		    cout << "InitialSpinUpSpace initialised with NbrBosonsUp="<<this->NbrBosonsUp<<", TotalLz " << SpinUpTotalLzOfPartition << " LzMax " << SpinUpMaxLzOfPartition << "\n";
		    cout << "InitialSpinDownSpace initialised with NbrBosonsDown="<<this->NbrBosonsDown<<", TotalLz " << SpinDownTotalLzOfPartition << " LzMax " << SpinDownMaxLzOfPartition << "\n";

		   
		    BosonOnSphereShort * FinalSpinUpSpace = new BosonOnSphereShort(this->NbrBosonsUp, this->GetTotalLzUp(i) + SpinUpTotalLzOfPartition, this->LzMax + polarizedSpace->LzMax );
		    BosonOnSphereShort * FinalSpinDownSpace = new BosonOnSphereShort(this->NbrBosonsDown, this->GetTotalLzDown(i) + SpinDownTotalLzOfPartition, this->LzMax + polarizedSpace->LzMax );
		    
		    int lzmaxup = this->GetStateLzMaxUp(i);
		    int lzmaxdown = this->GetStateLzMaxDown(i);
		    cout << "Creating FinalSpinUpSpace NbrBosonsUp " << this->NbrBosonsUp << " TotalLz " << this->GetTotalLzUp(i) << " + " << SpinUpTotalLzOfPartition << " = " << this->GetTotalLzUp(i) + SpinUpTotalLzOfPartition << " lzmax " << lzmaxup << " + " << SpinUpMaxLzOfPartition << " = " << SpinUpMaxLzOfPartition + lzmaxup << "\n";
		    cout << "Arguments used: "<<this->NbrBosonsUp<<", "<<this->GetTotalLzUp(i) + SpinUpTotalLzOfPartition<<", "<<this->LzMax + polarizedSpace->LzMax<<endl;
		    cout << "Creating FinalSpinDownSpace NbrBosonsDown " << this->NbrBosonsDown << " TotalLz " << this->GetTotalLzDown(i) << " + " << SpinDownTotalLzOfPartition << " = " << this->GetTotalLzDown(i) + SpinDownTotalLzOfPartition << " lzmax " << lzmaxdown << " + " << SpinDownMaxLzOfPartition << " = "<< SpinDownMaxLzOfPartition + lzmaxdown << "\n";
		    cout << "Arguments used: "<<this->NbrBosonsDown<<", "<<this->GetTotalLzDown(i) + SpinDownTotalLzOfPartition<<", "<<this->LzMax + polarizedSpace->LzMax<<endl;

		    
		    
		    unsigned long * FinalUpStates = new unsigned long [FinalSpinUpSpace->GetHilbertSpaceDimension() ];
		    unsigned long * FinalDownStates = new unsigned long [FinalSpinDownSpace->GetHilbertSpaceDimension() ];
		    long * UpWeight = new long [FinalSpinUpSpace->GetHilbertSpaceDimension() ];
		    long * DownWeight = new long [FinalSpinDownSpace->GetHilbertSpaceDimension() ];
    
		    InitialSpinUpSpace->ConvertToMonomial( upConversionFromSmallIntegerArray, SpinUpMaxLzOfPartition, SpinUpPartitionExponents );
		    InitialSpinDownSpace->ConvertToMonomial( downConversionFromSmallIntegerArray, SpinDownMaxLzOfPartition, SpinDownPartitionExponents );
		    
		    //multiply monomials within each spin group
		    unsigned long NbrUpStates = InitialSpinUpSpace->ProductOfTwoMonomials( FirstUpMonomials, EqualPowerIndexUp, NbrEqualPowerUp, SpinUpPartitionExponents, FinalUpStates, UpWeight, FinalSpinUpSpace); 
		    unsigned long NbrDownStates = InitialSpinDownSpace->ProductOfTwoMonomials( FirstDownMonomials, EqualPowerIndexDown, NbrEqualPowerDown, SpinDownPartitionExponents, FinalDownStates, DownWeight, FinalSpinDownSpace);
		    
		    //find all the combinations of spin products and add the result to the output vector
		    for( unsigned long upIndex = 0x0ul; upIndex < NbrUpStates; upIndex++) 
		      {
			//unsigned long * FinalStateUpOccupationNos = new unsigned long [finalSpace->LzMax + 1];
			for( int lz = 0; lz <= finalSpace->LzMax; lz++ )
			  FinalStateUpOccupationNos[lz] = 0x0ul;
			//	cout <<  "FinalSpinUpSpace->LzMax " << FinalSpinUpSpace->LzMax << "\n";
			//unsigned long * FinalStateUpExponents = new unsigned long [this->NbrBosonsUp];
			for(int p=0; p < this->NbrBosonsUp; p++)
			  FinalStateUpExponents[p] = 0x0ul;
			this->ConvertToMonomial( FinalUpStates[upIndex], 0x0ul, FinalSpinUpSpace->LzMax, FinalStateUpExponents );
		
			this->ConvertPolarizedMonomialToOccupationBasis( FinalStateUpExponents, FinalStateUpOccupationNos, this->NbrBosonsUp );
					
			for( unsigned long downIndex = 0x0ul; downIndex < NbrDownStates; downIndex++ ) 
			  {
			    // unsigned long * FinalStateDownOccupationNos = new unsigned long [finalSpace->LzMax + 1];
			    for( int lz = 0; lz <= finalSpace->LzMax; lz++ )
			      FinalStateDownOccupationNos[lz] = 0x0ul;
			    // unsigned long * FinalStateDownExponents = new unsigned long [this->NbrBosonsDown];
			    for(int p=0; p<this->NbrBosonsDown; p++)
			      FinalStateDownExponents[p] = 0x0ul;
			    this->ConvertToMonomial( 0x0ul, FinalDownStates[downIndex], FinalSpinDownSpace->LzMax, FinalStateDownExponents );
			    this->ConvertPolarizedMonomialToOccupationBasis( FinalStateDownExponents, FinalStateDownOccupationNos, this->NbrBosonsDown );
			    
			    
			    double geometricCorrectionFactor = this->GeometricCorrectionFactor( FirstPolarizedMonomials, FirstMonomials, FinalStateUpExponents, FinalStateDownExponents, polarizedSpace->LzMax, this->LzMax );
			    double occupationCorrectionFactor = this->OccupationCorrectionFactor( FirstPolarizedMonomialsOccupationBasis, SpinfulUpOccupationNos, SpinfulDownOccupationNos, FinalStateUpOccupationNos, FinalStateDownOccupationNos, polarizedSpace->LzMax, this->LzMax );
			    
			    outputVector[ finalSpace->FindStateIndex( FinalUpStates[upIndex], FinalDownStates[downIndex] ) ] += spinfulState[i]*polarizedState[j]*UpWeight[upIndex]*DownWeight[downIndex]*geometricCorrectionFactor*occupationCorrectionFactor;
			  }
		      }
		    
		    delete InitialSpinUpSpace;
		    delete InitialSpinDownSpace;
		    delete FinalSpinUpSpace;
		    delete FinalSpinDownSpace;
		    delete [] FinalUpStates;
		    delete [] FinalDownStates;
		    delete [] UpWeight;
		    delete [] DownWeight;
		  }
	      }
	  }
	}
    }
  }
  
  delete [] FinalStateUpOccupationNos;
  delete [] FinalStateUpExponents;

  delete [] FinalStateDownOccupationNos;
  delete [] FinalStateDownExponents;

  delete [] upConversionFromSmallIntegerArray;
  delete [] downConversionFromSmallIntegerArray;

  delete [] EqualPowerIndexUp;
  delete [] EqualPowerIndexDown;
  
  delete [] SpinfulUpOccupationNos;
  delete [] SpinfulDownOccupationNos;

  delete [] FirstMonomials;
  delete [] FirstPolarizedMonomials;
  delete [] FirstPolarizedMonomialsOccupationBasis;
  delete [] FirstPolarizedMonomialsOccupationBasisInt;

  delete [] SpinUpPartitionExponents;
  delete [] SpinDownPartitionExponents;
}



void BosonOnSphereWithSpin::SymmetriseOverGroupsAndAddToVector(unsigned long * & PlusStateUp, unsigned long * & MinusStateUp, unsigned long * & PlusStateDown, unsigned long * & MinusStateDown, double coefficient, RealVector & OutputVector )
{
  unsigned long * FinalUpState = new unsigned long [this->LzMax+1];
  unsigned long * FinalDownState = new unsigned long [this->LzMax+1];
  FactorialCoefficient occupationCorrection( 1l );
  for(int lz=0; lz <= this->LzMax; lz++)
    {
      FinalUpState[lz] = PlusStateUp[lz] + MinusStateUp[lz];
      occupationCorrection.FactorialMultiply( FinalUpState[lz] );
      occupationCorrection.FactorialDivide( PlusStateUp[lz] );
      occupationCorrection.FactorialDivide( MinusStateUp[lz] );

      FinalDownState[lz] = PlusStateDown[lz] + MinusStateDown[lz];
      occupationCorrection.FactorialMultiply( FinalDownState[lz] );
      occupationCorrection.FactorialDivide( PlusStateDown[lz] );
      occupationCorrection.FactorialDivide( MinusStateDown[lz] );
    }
  double occCorrFactor = sqrt( occupationCorrection.GetNumericalValue() );

  long index = this->FindStateIndex( FinalUpState, FinalDownState );
  OutputVector[ index ] += occCorrFactor*coefficient;
  delete [] FinalUpState;
  delete [] FinalDownState;

}

//Compute the geometric correction factor for a given product state when multiplying two monomials and working with second quantised forms on the sphere for a given spin value for the case of a fully polarized state
//
//firstState = reference on array where monomial representation of first state stored
//secondState = reference on array where monomial representation of second state stored
//productState = reference on array where monomial representation of a given final state in the product is stored
//lzMaxOne = twice maximum lz value for a boson in first state
//lzMaxTwo = twice maximum lz value for a boson in second state
double BosonOnSphereWithSpin::GeometricCorrectionFactor(unsigned long * firstMonomial, unsigned long * secondMonomial, unsigned long * productMonomial, int lzMaxOne, int lzMaxTwo)
{
  //cout << "N_phi polarized " << lzMaxOne << " N_phi spinful " << lzMaxTwo << " N_phi final " << lzMaxOne + lzMaxTwo << "\n";
  FactorialCoefficient factor(1l);
  for(int particle=0; particle<this->NbrBosons; particle++)
    {
      // cout << productMonomial[particle] << "!/" << firstMonomial[particle] << "!" << secondMonomial[particle] << "!\n";
      // cout << lzMaxOne+lzMaxTwo-productMonomial[particle] << "!/" << lzMaxOne-firstMonomial[particle] << "!" << lzMaxTwo-secondMonomial[particle] << "!\n";
      factor.FactorialMultiply( productMonomial[particle] );
      factor.FactorialDivide( firstMonomial[particle] );
      factor.FactorialDivide( secondMonomial[particle] );
      factor.FactorialMultiply( lzMaxOne + lzMaxTwo - productMonomial[particle] );
      factor.FactorialDivide( lzMaxOne - firstMonomial[particle] );
      factor.FactorialDivide( lzMaxTwo - secondMonomial[particle] );
    }
  //cout << "GeometricCorrectionFactor squared " << factor << "\n";
  double SquaredResult = factor.GetNumericalValue();
  //cout << "GeometricCorrectionFactor " << sqrt(SquaredResult) << "\n";
  return sqrt(SquaredResult);
}

//Compute the geometric correction factor for a given product state when multiplying two monomials and working with second quantised forms on the sphere for a given spin value for the case of a non-fully polarized spin state
//
//firstState = reference on array where monomial representation of first state stored
//secondState = reference on array where monomial representation of second state stored
//productState = reference on array where monomial representation of a given final state in the product is stored
//lzMaxOne = twice maximum lz value for a boson in first state
//lzMaxTwo = twice maximum lz value for a boson in second state
double BosonOnSphereWithSpin::GeometricCorrectionFactor(unsigned long * firstMonomial, unsigned long * secondMonomial, unsigned long * productMonomialUp, unsigned long * productMonomialDown, int lzMaxOne, int lzMaxTwo)
{
  //cout << "N_phi polarized " << lzMaxOne << " N_phi spinful " << lzMaxTwo << " N_phi final " << lzMaxOne + lzMaxTwo << "\n";

  FactorialCoefficient factor(1l);
  for(int particle=0; particle < this->NbrBosonsUp; particle++)
    {
      factor.FactorialMultiply( productMonomialUp[particle] );
      factor.FactorialDivide( firstMonomial[particle] );
      factor.FactorialDivide( secondMonomial[particle] );
      factor.FactorialMultiply( lzMaxOne + lzMaxTwo - productMonomialUp[particle] );
      factor.FactorialDivide( lzMaxOne - firstMonomial[particle] );
      factor.FactorialDivide( lzMaxTwo - secondMonomial[particle] );
    }
  for(int particle=0; particle < this->NbrBosonsDown; particle++)
    {
      factor.FactorialMultiply( productMonomialDown[particle] );
      factor.FactorialDivide( firstMonomial[this->NbrBosonsUp + particle] );
      factor.FactorialDivide( secondMonomial[this->NbrBosonsUp + particle] );
      factor.FactorialMultiply( lzMaxOne + lzMaxTwo - productMonomialDown[particle] );
      factor.FactorialDivide( lzMaxOne - firstMonomial[this->NbrBosonsUp + particle] );
      factor.FactorialDivide( lzMaxTwo - secondMonomial[this->NbrBosonsDown + particle] );
    }
  //cout << "GeometricCorrectionFactor squared " << factor << "\n";

  double SquaredResult = factor.GetNumericalValue();
  //cout << "GeometricCorrectionFactor " << sqrt(SquaredResult) << "\n";
  return sqrt(SquaredResult);
}

//Compute the occupation correction factor for a given product state when multiplying two monomials and working with second quantised forms on the sphere for the case of fully polarized states
//
//firstState = reference on array where bosonic representation of first state stored
//secondState = reference on array where bosonic representation of second state stored
//productState = reference on array where bosonic representation of a given final state in the product is stored
//lzMaxOne = twice maximum lz value for a boson in first state
//lzMaxTwo = twice maximum lz value for a boson in second state
double BosonOnSphereWithSpin::OccupationCorrectionFactor(unsigned long * firstState, unsigned long * secondState, unsigned long * productState, int lzMaxOne, int lzMaxTwo)
{
  FactorialCoefficient factor(1l);
  int lz = 0;
  if( lzMaxOne < lzMaxTwo ) 
    {
      while( lz <= lzMaxOne )
	{
	  factor.FactorialMultiply( firstState[lz] );
	  factor.FactorialMultiply( secondState[lz] );
	  factor.FactorialDivide( productState[lz] );
	  lz++;
	}
      while( lz <= lzMaxTwo )
	{
	  factor.FactorialMultiply( secondState[lz] );
	  factor.FactorialDivide( productState[lz] );
	  lz++;
	}
      while(lz <= lzMaxOne + lzMaxTwo)
	{
	  factor.FactorialDivide( productState[lz] );
	  lz++;
	}
    }
  else
    {
      while( lz <= lzMaxTwo )
	{
	  factor.FactorialMultiply( firstState[lz] );
	  factor.FactorialMultiply( secondState[lz] );
	  factor.FactorialDivide( productState[lz] );
	  lz++;
	}
      while( lz <= lzMaxOne )
	{
	  factor.FactorialMultiply( firstState[lz] );
	  factor.FactorialDivide( productState[lz] );
	  lz++;
	}
      while( lz <= lzMaxOne + lzMaxTwo )
	{
	  factor.FactorialDivide( productState[lz] );
	  lz++;
	}
    }
  factor.FactorialDivide( this->NbrBosons );
  double squaredResult = factor.GetNumericalValue();
  //cout << "OccupationCorrectionFactor " << sqrt(squaredResult) << "\n"; 
  return sqrt( squaredResult );
}


//Compute the occupation correction factor for a given product state when multiplying two monomials and working with second quantised forms on the sphere for the case of spinful states
//
//firstState = reference on array where bosonic representation of first state stored
//secondState = reference on array where bosonic representation of second state stored
//productState = reference on array where bosonic representation of a given final state in the product is stored
//lzMaxOne = twice maximum lz value for a boson in first state
//lzMaxTwo = twice maximum lz value for a boson in second state

double BosonOnSphereWithSpin::OccupationCorrectionFactor(unsigned long * firstPolarizedState, unsigned long * secondUpState, 
							 unsigned long * secondDownState, unsigned long * productUpState, unsigned long * productDownState, 
							 int lzMaxOne, int lzMaxTwo)
{
  FactorialCoefficient factor(1l);
  int lz = 0;
  if( lzMaxOne < lzMaxTwo ) 
    {
      while( lz <= lzMaxOne )
	{
	  factor.FactorialMultiply( firstPolarizedState[lz] );
	  factor.FactorialMultiply( secondUpState[lz] );
	  factor.FactorialMultiply( secondDownState[lz] );
	  factor.FactorialDivide( productUpState[lz] );
	  factor.FactorialDivide( productDownState[lz] );

	  lz++;
	}
      while( lz <= lzMaxTwo )
	{
	  factor.FactorialMultiply( secondUpState[lz] );
	  factor.FactorialMultiply( secondDownState[lz] );
	  factor.FactorialDivide( productUpState[lz] );
	  factor.FactorialDivide( productDownState[lz] );
	  lz++;
	}
      while(lz <= lzMaxOne + lzMaxTwo)
	{
	  factor.FactorialDivide( productUpState[lz] );
	  factor.FactorialDivide( productDownState[lz] );
	  lz++;
	}
    }
  else
    {
      while( lz <= lzMaxTwo )
	{
	  factor.FactorialMultiply( firstPolarizedState[lz] );
	  factor.FactorialMultiply( secondUpState[lz] );
	  factor.FactorialMultiply( secondDownState[lz] );
	  factor.FactorialDivide( productUpState[lz] );
	  factor.FactorialDivide( productDownState[lz] );
	  lz++;
	}
      while( lz <= lzMaxOne )
	{
	  factor.FactorialMultiply( firstPolarizedState[lz] );
	  factor.FactorialDivide( productUpState[lz] );
	  factor.FactorialDivide( productDownState[lz] );
	  lz++;
	}
      while( lz <= lzMaxOne + lzMaxTwo )
	{
	  factor.FactorialDivide( productUpState[lz] );
	  factor.FactorialDivide( productDownState[lz] );
	  lz++;
	}
    }
  factor.FactorialDivide(this->NbrBosons);
  double squaredResult = factor.GetNumericalValue();
  //cout << "OccupationCorrectionFactor " << sqrt(squaredResult) << "\n"; 
  return sqrt( squaredResult );
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& BosonOnSphereWithSpin::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  int TemporaryStateLzMaxUp;
  int TemporaryStateLzMaxDown;
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0 / state[reference];
  state[reference] = 1.0;
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  cout << hex << this->StateDescriptionUp[reference] << " " << this->StateDescriptionDown[reference] << dec << " " << ((this->StateInfo[reference]  >> 10) & 0x3ffu) << endl;
  this->ConvertToMonomial(this->StateDescriptionUp[reference], this->StateDescriptionDown[reference], (this->StateInfo[reference] >> 10) & 0x3ffu, TmpMonomialReference);
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(this->StateDescriptionUp[reference], this->StateDescriptionDown[reference], this->LzMax, 
		       this->TemporaryState, TemporaryStateLzMaxUp, TemporaryStateLzMaxDown);
  for (int k = 0; k <= TemporaryStateLzMaxUp; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialMultiply(this->TemporaryState[k] >> 16);
  for (int k = 0; k <= TemporaryStateLzMaxDown; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialMultiply(this->TemporaryState[k] & 0xffff);
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->StateDescriptionUp[i], this->StateDescriptionDown[i], this->LzMax, TmpMonomial);
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	{
	  while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrBosons)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrBosons)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	  ++Index2;
	}
      if (symmetryFactor == true)
	{
	  Factorial = ReferenceFactorial;
	  this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i], this->StateInfo[i], 
			       this->TemporaryState, TemporaryStateLzMaxUp, TemporaryStateLzMaxDown);
	  for (int k = 0; k <= TemporaryStateLzMaxUp; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialMultiply(this->TemporaryState[k] >> 16);
	  for (int k = 0; k <= TemporaryStateLzMaxDown; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialMultiply(this->TemporaryState[k] & 0xffff);
	  Coefficient *= sqrt(Factorial.GetNumericalValue());
	}
      state[i] *= Coefficient;
    }
  delete[] TmpMonomialReference;
  delete[] TmpMonomial;
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& BosonOnSphereWithSpin::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  int* TmpMonomialReference = new int [this->NbrBosons];
  int* TmpMonomial = new int [this->NbrBosons];
  double Factor = 1.0;
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      SqrtCoefficients[k] =  1.0 / InvSqrtCoefficients[k];
    }
//   unsigned long TmpState = this->StateDescription[reference];
//   int Index = 0;
//   for (int j = this->LzMax; j >= 0; --j)
//     {
//       switch ((TmpState >> (2 * j)) & 3ul)
// 	{
// 	case 0x1ul:
// 	  TmpMonomialReference[Index++] = j;
// 	  break;
// 	case 0x2ul:
// 	  TmpMonomialReference[Index++] = j;
// 	  break;
// 	case 0x3ul:
// 	  {
// 	    TmpMonomialReference[Index++] = j;
// 	    TmpMonomialReference[Index++] = j;
// 	  }
// 	  break;
// 	}
//     }
//   for (int i = 1; i < this->HilbertSpaceDimension; ++i)
//     {
//       Index = 0;
//       TmpState = this->StateDescription[i];
//       for (int j = this->LzMax; j >= 0; --j)
// 	{
// 	  switch ((TmpState >> (2 * j)) & 3ul)
// 	    {
// 	    case 0x1ul:
// 	      TmpMonomial[Index++] = j;
// 	      break;
// 	    case 0x2ul:
// 	      TmpMonomial[Index++] = j;
// 	      break;
// 	    case 0x3ul:
// 	      {
// 		TmpMonomial[Index++] = j;
// 		TmpMonomial[Index++] = j;
// 	      }
// 	      break;
// 	    }
// 	}
//       int Index1 = 0;
//       int Index2 = 0;
//       double Coefficient = Factor;
//       while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
// 	{
// 	  while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
// 	    {
// 	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
// 	      ++Index1;
// 	    }
// 	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
// 	    {
// 	      ++Index1;
// 	      ++Index2;
// 	    }
// 	  while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
// 	    {
// 	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
// 	      ++Index2;
// 	    }	  
// 	}
//       while (Index1 < this->NbrBosons)
// 	{
// 	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
// 	  ++Index1;
// 	}
//       while (Index2 < this->NbrBosons)
// 	{
// 	  Coefficient *= SqrtCoefficients[TmpMonomialReference[Index2]];
// 	  ++Index2;
// 	}
//       state[i] *= Coefficient;
//     }
//   delete[] TmpMonomialReference;
//   delete[] TmpMonomial;
  state /= state.Norm();
  return state;
}

// normalize Jack with respect to cylinder basis
//
// state = reference to the Jack state to normalize
// aspect = aspect ratio of cylinder
// return value = normalized state

RealVector& BosonOnSphereWithSpin::NormalizeJackToCylinder(RealVector& state, double aspect)
{
//   long double Pi_L = 3.14159265358979323846264338328L;
//   long double Length = sqrtl((long double)2.0 * Pi_L * (long double)(this->LzMax + 1) * (long double) aspect);
//   cout<<"L= "<< Length << " r= "<<aspect<<endl;
//   long double kappa = ((long double) 2.0) * Pi_L / Length;
//   long double Norm = (long double)0.0;
//   long double* LogFactorials = new long double [this->NbrBosons + 1];
//   LogFactorials[0] = (long double) 0.0;
//   LogFactorials[1] = (long double) 0.0;
//   for (int i = 2; i <= this->NbrBosons; ++i)
//     {
//       LogFactorials[i] = LogFactorials[i - 1] + logl((long double) i);
//     }
//   for (int i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//     {
//       this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
// 			   this->TemporaryState, this->TemporaryStateLzMax);
//       long double Sum2MSquare = (long double) 0.0;
//       long double TmpFactorial = (long double) 0.0;
//       int j= 0;
//       for (int j = 0; j <= this->TemporaryStateLzMax; ++j)
// 	{
// 	  if (this->TemporaryState[j] > 0)
// 	    {
// 	      Sum2MSquare += (((j - 0.5 * LzMax) * (j - 0.5 * LzMax)) *  this->TemporaryState[j]);
// 	      TmpFactorial -= LogFactorials[this->TemporaryState[j]]; 
// 	    }
// 	}
//       state[i] *= expl((((long double)0.5) * kappa * kappa *Sum2MSquare) + (((long double)0.5) * TmpFactorial));      
//       Norm += state[i] * state[i];
//    }
//   delete[] LogFactorials;
//   cout << "Norm= " << Norm << endl;
//   state /= sqrtl(Norm);
 
  return state;
}

