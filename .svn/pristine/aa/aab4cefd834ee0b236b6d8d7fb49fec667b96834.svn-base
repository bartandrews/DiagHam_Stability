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
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/StringTools.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

// testing flag - switches debut output
//#define TESTING

// default constructor
//

BosonOnSphereWithSpinAllSz::BosonOnSphereWithSpinAllSz ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson
// memory = amount of memory granted for precalculations
//
BosonOnSphereWithSpinAllSz::BosonOnSphereWithSpinAllSz (int nbrBosons, int totalLz, int lzMax, unsigned long memory)
{
  cout << "BosonOnSphereWithSpinAllSz"<<endl;
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->SzParity = -1;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, this->TotalLz);
  cout << "dim = " << this->HilbertSpaceDimension << endl;
  this->TemporaryState = new unsigned [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned [this->NbrLzValue];
  this->Flag.Initialize();
  this->TargetSpace = this;
  this->StateDescription = new unsigned* [this->HilbertSpaceDimension];
  this->StateLzMaxUp = new unsigned [this->HilbertSpaceDimension];
  this->StateLzMaxDown = new unsigned [this->HilbertSpaceDimension];
  this->StateNbrUp = new unsigned [this->HilbertSpaceDimension];
  this->SubspaceRestriction = true;
  this->SubspaceSz = NbrBosons & 1;
  int TmpLzMax = this->LzMax;
  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  int StartLzMax = (LzMax<<1)+1;
  int TmpDim = this->GenerateStates(this->NbrBosons, LzMax, StartLzMax, this->ShiftedTotalLz, 0);
  if (TmpDim!=this->HilbertSpaceDimension)
    cout << "Count inconsistent: "<<TmpDim<<" vs " << this->HilbertSpaceDimension<<endl;
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
  UsedMemory += 2*this->HilbertSpaceDimension * sizeof(unsigned long);  // StateDescriptionUp/Down
  UsedMemory += this->HilbertSpaceDimension * sizeof(unsigned); // StateInfo
  // 
  UsedMemory += this->NbrTensoredElements * (sizeof(unsigned long) + 2*sizeof(unsigned) + sizeof(int));  // TensoredElements, UpSpinLookUpTable  , DownSpinLookUpTable, TensoredLzMax
  UsedMemory += (this->TotalLz + this->NbrBosons) * (LookUpTableMemorySize+1) * sizeof(int); // LookUpTableShift, LookUpTable
  cout << "HS using ";
  PrintMemorySize(cout, UsedMemory)<<endl;
#endif
}


// basic constructor with pair parity
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson
// szParity = parity of the total Sz : (-1)^Sz == 1 (even parity 0 generalized for all N by: N_\up mod 2 == \lgauss N/2 \rgauss mod 2)
// memory = amount of memory granted for precalculations
//
BosonOnSphereWithSpinAllSz::BosonOnSphereWithSpinAllSz (int nbrBosons, int totalLz, int lzMax, int szParity, unsigned long memory)
{
  cout << "BosonOnSphereWithSpinAllSz with pair parity"<<endl;
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedParity = ((nbrBosons>>1) + szParity) & 0x1;
  this->NbrLzValue = this->LzMax + 1;
  this->SzParity = szParity;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, this->TotalLz, this->SzParity);
  cout << "dim = " << this->HilbertSpaceDimension << endl;
  this->TemporaryState = new unsigned [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned [this->NbrLzValue];
  this->Flag.Initialize();
  this->TargetSpace = this;
  this->StateDescription = new unsigned* [this->HilbertSpaceDimension];
  this->StateLzMaxUp = new unsigned [this->HilbertSpaceDimension];
  this->StateLzMaxDown = new unsigned [this->HilbertSpaceDimension];
  this->StateNbrUp = new unsigned [this->HilbertSpaceDimension];
  this->SubspaceRestriction = true;
  this->SubspaceSz = NbrBosons & 1;
  int TmpLzMax = this->LzMax;
  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  int StartLzMax = (LzMax<<1)+1;
  int TmpDim = this->GenerateStates(this->NbrBosons, LzMax, StartLzMax, this->ShiftedTotalLz, 0, ShiftedParity, 0);
  if (TmpDim!=this->HilbertSpaceDimension)
    cout << "Count inconsistent: "<<TmpDim<<" vs " << this->HilbertSpaceDimension<<endl;
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
  UsedMemory += 2*this->HilbertSpaceDimension * sizeof(unsigned long);  // StateDescriptionUp/Down
  UsedMemory += this->HilbertSpaceDimension * sizeof(unsigned); // StateInfo
  // 
  UsedMemory += this->NbrTensoredElements * (sizeof(unsigned long) + 2*sizeof(unsigned) + sizeof(int));  // TensoredElements, UpSpinLookUpTable  , DownSpinLookUpTable, TensoredLzMax
  UsedMemory += (this->TotalLz + this->NbrBosons) * (LookUpTableMemorySize+1) * sizeof(int); // LookUpTableShift, LookUpTable
  cout << "HS using ";
  PrintMemorySize(cout, UsedMemory)<<endl;
#endif
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSphereWithSpinAllSz::BosonOnSphereWithSpinAllSz(const BosonOnSphereWithSpinAllSz& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->SzParity = bosons.SzParity;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->StateLzMaxUp = bosons.StateLzMaxUp;
  this->StateLzMaxDown = bosons.StateLzMaxDown;
  this->StateNbrUp = bosons.StateNbrUp;
  this->StateInfo = bosons.StateInfo;
  this->UpSpinLookUpTable = bosons.UpSpinLookUpTable;
  this->DownSpinLookUpTable = bosons.DownSpinLookUpTable;
  this->TensoredElements = bosons.TensoredElements;
  this->TensoredLzMax = bosons.TensoredLzMax;
  this->NbrTensoredElements = bosons.NbrTensoredElements;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->CoherenceFactors = bosons.CoherenceFactors;
  this->Flag = bosons.Flag;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned [this->NbrLzValue];
  this->SubspaceRestriction=bosons.SubspaceRestriction;
  this->SubspaceSz=bosons.SubspaceSz;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

BosonOnSphereWithSpinAllSz::~BosonOnSphereWithSpinAllSz ()
{
  delete[] this->TemporaryState;
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
      delete[] this->TensoredElements;
      delete[] this->TensoredLzMax;
      delete[] this->UpSpinLookUpTable;
      delete[] this->DownSpinLookUpTable;
      if (this->StateLzMaxUp!=NULL)
	delete[] this->StateLzMaxUp;
      if (this->StateLzMaxDown!=NULL)
	delete[] this->StateLzMaxDown;
      if (this->StateNbrUp!=NULL)
	delete[] this->StateNbrUp;
      if (this->LookUpTableShift != 0)
	{
	  delete[] this->LookUpTableShift;
	  for (int i = 0; i < NbrBosons + LzMax; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;
	}
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

BosonOnSphereWithSpinAllSz& BosonOnSphereWithSpinAllSz::operator = (const BosonOnSphereWithSpinAllSz& bosons)
{
  delete[] this->TemporaryState;
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
      delete[] this->TensoredElements;
      delete[] this->TensoredLzMax;
      delete[] this->UpSpinLookUpTable;
      delete[] this->DownSpinLookUpTable;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
      delete[] this->CoherenceFactors;
      if (this->StateLzMaxUp!=NULL)
	delete[] this->StateLzMaxUp;
      if (this->StateLzMaxDown!=NULL)
	delete[] this->StateLzMaxDown;
      if (this->StateNbrUp!=NULL)
	delete[] this->StateNbrUp;
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
  this->SzParity = bosons.SzParity;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->StateLzMaxUp = bosons.StateLzMaxUp;
  this->StateLzMaxDown = bosons.StateLzMaxDown;
  this->StateNbrUp = bosons.StateNbrUp;
  this->StateInfo = bosons.StateInfo;
  this->UpSpinLookUpTable = bosons.UpSpinLookUpTable;
  this->DownSpinLookUpTable = bosons.DownSpinLookUpTable;
  this->TensoredElements = bosons.TensoredElements;
  this->TensoredLzMax = bosons.TensoredLzMax;
  this->NbrTensoredElements = bosons.NbrTensoredElements;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;
  this->CoherenceFactors = bosons.CoherenceFactors;
  this->Flag = bosons.Flag;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned [this->NbrLzValue];
  this->SubspaceRestriction=bosons.SubspaceRestriction;
  this->SubspaceSz=bosons.SubspaceSz;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereWithSpinAllSz::Clone()
{
  return new BosonOnSphereWithSpinAllSz(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereWithSpinAllSz::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereWithSpinAllSz::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereWithSpinAllSz::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnSphereWithSpinAllSz::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
  this->TargetSpace = (BosonOnSphereWithSpinAllSz*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int BosonOnSphereWithSpinAllSz::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
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

int BosonOnSphereWithSpinAllSz::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

int BosonOnSphereWithSpinAllSz::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
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

int BosonOnSphereWithSpinAllSz::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpinAllSz::AduAdu (int m1, int m2, double& coefficient)
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
  return this->TargetSpace->FindStateIndex(this->TemporaryState, ProdATemporaryStateNbrUp + 2);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AdAd method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpinAllSz::AddAdd (int m1, int m2, double& coefficient)
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

  return this->TargetSpace->FindStateIndex(this->TemporaryState, ProdATemporaryStateNbrUp);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAd method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSpinAllSz::AduAdd (int m1, int m2, double& coefficient)
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
  return this->TargetSpace->FindStateIndex(this->TemporaryState, ProdATemporaryStateNbrUp + 1);
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 

double BosonOnSphereWithSpinAllSz::AuAu (int index, int n1, int n2)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       ProdATemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, ProdATemporaryStateNbrUp);
  if ((n1 > CurrentLzMaxUp) || (n2 > CurrentLzMaxUp) || ((ProdATemporaryState[n1] >> 16) == 0)
      || ((ProdATemporaryState[n2] >> 16) == 0) || ((n1 == n2) && ((ProdATemporaryState[n1] >> 16) == 1)))
    return 0.0;
  //double Coefficient = (this->ProdATemporaryState[n2] >> 16);
  int CoherenceIndex = (this->ProdATemporaryState[n2] >> 16);
  this->ProdATemporaryState[n2] -= 0x10000;
  //Coefficient *= (this->ProdATemporaryState[n1] >> 16);
  CoherenceIndex *= (this->ProdATemporaryState[n1] >> 16);
  this->ProdATemporaryState[n1] -= 0x10000;
  ProdATemporaryStateNbrUp-=2;
  return this->CoherenceFactors[CoherenceIndex];
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double BosonOnSphereWithSpinAllSz::AdAd (int index, int n1, int n2)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       ProdATemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, ProdATemporaryStateNbrUp);

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

double BosonOnSphereWithSpinAllSz::AuAd (int index, int n1, int n2)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       ProdATemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, ProdATemporaryStateNbrUp);
  if ((n1 > CurrentLzMaxUp) || (n2 > CurrentLzMaxDown) || ((ProdATemporaryState[n1] >> 16) == 0)
      || ((ProdATemporaryState[n2] & 0xffff) == 0))
    return 0.0;
  //double Coefficient = (this->ProdATemporaryState[n2] & 0xffff);
  int CoherenceIndex = (this->ProdATemporaryState[n2] & 0xffff);
  --this->ProdATemporaryState[n2];
  //Coefficient *= (this->ProdATemporaryState[n1] >> 16);
  CoherenceIndex *= (this->ProdATemporaryState[n1] >> 16);
  this->ProdATemporaryState[n1] -= 0x10000;
  --this->ProdATemporaryStateNbrUp;
  return this->CoherenceFactors[CoherenceIndex];
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereWithSpinAllSz::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       ProdATemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, ProdATemporaryStateNbrUp);
  
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
	  --this->ProdATemporaryStateNbrUp;
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

int BosonOnSphereWithSpinAllSz::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  int i = 0;
  for (; i < this->NbrLzValue; ++i)
    {
      this->TemporaryState[i] = this->ProdATemporaryState[i];
    }
  int TmpCoefficient = 1;
  int TemporaryStateNbrUp = ProdATemporaryStateNbrUp;
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
	++TemporaryStateNbrUp;
      }
  coefficient = sqrt((double) TmpCoefficient);
  return this->TargetSpace->FindStateIndex(this->TemporaryState, TemporaryStateNbrUp);
}


// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereWithSpinAllSz::AduAu (int index, int m)
{
  unsigned Info = this->StateInfo[index];
  unsigned TemporaryStateNbrUp = Info&0x03ff;
  int CurrentLzMaxUp = (Info>>20) - TemporaryStateNbrUp + (TemporaryStateNbrUp!=0);
  if (CurrentLzMaxUp < m)
    return 0.0;
  else
    {
      int CurrentLzMaxDown;
      this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], Info,
			   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
      return (double) ((this->TemporaryState[m] >> 16));
    }
}

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereWithSpinAllSz::AddAd (int index, int m)
{
  unsigned Info = this->StateInfo[index];
  unsigned TemporaryStateNbrUp = Info & 0x03ff;	
  int CurrentLzMaxDown = ((Info>>10)&0x03ff) - NbrBosons + TemporaryStateNbrUp + (TemporaryStateNbrUp!=(unsigned)NbrBosons);
  if (CurrentLzMaxDown < m)
    {
      return 0.0;
    }
  else
    {
      int CurrentLzMaxUp;
      this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], Info,
			   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
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
int BosonOnSphereWithSpinAllSz::AduAu (int index, int m, int n, double& coefficient)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  unsigned TemporaryStateNbrUp;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);

  if ((n > CurrentLzMaxUp) || ((TemporaryState[n] >> 16) == 0)) // || ((State[n] & 0xffff) == 0)) // shift for up, mask for down
    {
      coefficient=0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  //double TmpCoefficient = (this->TemporaryState[n] >> 16);
  int CoherenceIndex = (this->TemporaryState[n] >> 16);
  this->TemporaryState[n] -= 0x10000;
  
  this->TemporaryState[m] += 0x10000;
  //TmpCoefficient *= (this->TemporaryState[m] >> 16);
  CoherenceIndex *= (this->TemporaryState[m] >> 16);  
  coefficient = this->CoherenceFactors[CoherenceIndex];
  return this->TargetSpace->FindStateIndex(this->TemporaryState, TemporaryStateNbrUp);
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnSphereWithSpinAllSz::AddAd (int index, int m, int n, double& coefficient)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  unsigned TemporaryStateNbrUp;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);

  if ((n > CurrentLzMaxDown) || ((TemporaryState[n] & 0xffff) == 0)) // || ((State[n] & 0xffff) == 0)) // shift for up, mask for down
    {
      coefficient=0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  //double TmpCoefficient = (this->TemporaryState[n] & 0xffff);
  int CoherenceIndex = (this->TemporaryState[n] & 0xffff);
  --this->TemporaryState[n];
  ++this->TemporaryState[m];
  //TmpCoefficient *= (this->TemporaryState[m] & 0xffff);
  CoherenceIndex *= (this->TemporaryState[m] & 0xffff);
  coefficient = this->CoherenceFactors[CoherenceIndex];
  return this->TargetSpace->FindStateIndex(this->TemporaryState, TemporaryStateNbrUp);
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnSphereWithSpinAllSz::AduAd (int index, int m, int n, double& coefficient)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  unsigned TemporaryStateNbrUp;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);

  if ((n > CurrentLzMaxDown) || ((TemporaryState[n] & 0xffff) == 0)) // || ((State[n] & 0xffff) == 0)) // shift for up, mask for down
    {
      coefficient=0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  //double TmpCoefficient = (this->TemporaryState[n] & 0xffff);
  int CoherenceIndex = (this->TemporaryState[n] & 0xffff);
  --this->TemporaryState[n];
  
  this->TemporaryState[m] += 0x10000;
  //TmpCoefficient *= (this->TemporaryState[m] >> 16);
  CoherenceIndex *= (this->TemporaryState[m] >> 16);
  coefficient = this->CoherenceFactors[CoherenceIndex];
  return this->TargetSpace->FindStateIndex(this->TemporaryState, TemporaryStateNbrUp + 1);
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnSphereWithSpinAllSz::AddAu (int index, int m, int n, double& coefficient)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  unsigned TemporaryStateNbrUp;
  this->FermionToBoson(this->StateDescriptionUp[index], this->StateDescriptionDown[index], this->StateInfo[index],
		       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
  if ((n > CurrentLzMaxUp) || ((TemporaryState[n] >> 16) == 0)) // || ((State[n] & 0xffff) == 0)) // shift for up, mask for down
    {
      coefficient=0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  //double TmpCoefficient = (this->TemporaryState[n] >> 16);
  int CoherenceIndex = (this->TemporaryState[n] >> 16);
  this->TemporaryState[n] -= 0x10000;
  
  ++this->TemporaryState[m];
  //TmpCoefficient *= (this->TemporaryState[m]) & 0xffff;
  CoherenceIndex *= (this->TemporaryState[m]) & 0xffff;
  // coefficient *= this->CoherenceFactors[(this->TemporaryState[n] >> 16)*(this->TemporaryState[m]& 0xffff)];
  coefficient = this->CoherenceFactors[CoherenceIndex];
  return this->TargetSpace->FindStateIndex(this->TemporaryState, TemporaryStateNbrUp - 1);
}


// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnSphereWithSpinAllSz::FindStateIndex(unsigned* stateDescription, unsigned nbrUp)
{
  unsigned long FinalStateUp;
  unsigned long FinalStateDown;
  int LzMaxUp, LzMaxDown;
  this->BosonToFermion(FinalStateUp, FinalStateDown, LzMaxUp, LzMaxDown, stateDescription, nbrUp);
  //  cout << "up: "<< hex << FinalStateUp << dec << " deduced lzmax="<<LzMaxUp<<endl;
  //  cout << "down: "<< hex << FinalStateDown << dec << " deduced lzmax="<<LzMaxDown<<endl;
  return FindStateIndex(FinalStateUp, FinalStateDown, LzMaxUp, LzMaxDown);
}



// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int BosonOnSphereWithSpinAllSz::FindStateIndex(char* stateDescription)
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
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}



// sort an array and reflect permutations in auxiliary array
//
// length = length of arrays
// sortArray = array to be sorted
// auxArray = auxiliary array
//
void BosonOnSphereWithSpinAllSz::ShellSortAux(unsigned length, unsigned long* sortArray, unsigned *auxArray, int *auxArray2)
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
  
// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the Lz component
int BosonOnSphereWithSpinAllSz::GetLzValue(int j)
{
  return this->TotalLz;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSpinAllSz::PrintState (ostream& Str, int state)
{
  int CurrentLzMaxUp, CurrentLzMaxDown;
  unsigned TemporaryStateNbrUp;
  this->FermionToBoson(this->StateDescriptionUp[state], this->StateDescriptionDown[state], this->StateInfo[state],
		       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
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

ostream& BosonOnSphereWithSpinAllSz::PrintState (ostream& Str, unsigned *myState)
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
// nbrBosons = number of bosons remaining to be placed
// lzMax = momentum maximum value for a boson in the state 
// currentLzMax = momentum maximum value for bosons that are still to be placed
//                (counting from 0 to lzMax for down and lzMax+1 to 2LzMax+1 for up)
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

/*
  if ((nbrBosons < 0) || (totalLz < 0) || (currentLzMax < 0))// || ((CurrentSpinUp==0)&&((nbrBosons * CurrentMomentum) < totalLz)))
    {
      return 0;
    }
  if (nbrBosons == 0)
    {
      if (totalLz == 0)
	{
	  return 1;
	}
      else
	{
	  return 0;
	}
    }
  if ((CurrentSpinUp==0)&&(nbrBosons == 1))
    if (currentLzMax >= totalLz)
      {
	return 1;
      }
  if ((CurrentSpinUp==0)&&((nbrBosons * CurrentMomentum) == totalLz))
    {
      return 1;
    }
*/

int BosonOnSphereWithSpinAllSz::GenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int pos)
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
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// currentLzMax = momentum maximum value for bosons that are still to be placed
// totalLz = momentum total value
// totalNbrUp = number of up-particles placed so far
// szParity = target parity of totalNbrUp
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
int BosonOnSphereWithSpinAllSz::GenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int totalNbrUp, int szParity, int pos)
{
  // cout << "GenerateStates(n:"<< nbrBosons<<", currLz:"<<currentLzMax<< ", totalLz:" <<totalLz<<", pos:"<< pos<<")"<<endl;
  int CurrentMomentum = currentLzMax%this->NbrLzValue;
  int CurrentSpinUp = currentLzMax/this->NbrLzValue;
  if ((nbrBosons < 0) || (totalLz < 0) || (currentLzMax < 0) || ((CurrentSpinUp==0)&&((nbrBosons * CurrentMomentum) < totalLz)))
    return pos;
  if (nbrBosons == 0)
    {
      if ((totalLz == 0)&&((totalNbrUp & 0x1) == szParity))

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
      if ((currentLzMax >= totalLz)&&((totalNbrUp & 0x1) == szParity))
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
  if (((CurrentSpinUp==0)&&((nbrBosons * CurrentMomentum) == totalLz))&&((totalNbrUp & 0x1) == szParity))
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
      TmpPos = this->GenerateStates(nbrBosons-ToPlace, lzMax, ReducedCurrentLzMax, totalLz-ToPlace*CurrentMomentum, totalNbrUp+ToPlace*CurrentSpinUp, szParity, pos);
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
}


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereWithSpinAllSz::GenerateLookUpTable(int memory)
{
  // re-code storage of LzMaxUp etc in StateInfo & convert all states into packed storage
  this->StateInfo = new unsigned[this->HilbertSpaceDimension];
  this->StateDescriptionUp = new unsigned long[this->HilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long[this->HilbertSpaceDimension];
  unsigned NbrStatesUp=0;
  int LzMaxUp, LzMaxDown, CurrentStateNbrUp;
  unsigned long LastUp=~0x0ul;
  unsigned *TmpBlockStart = new unsigned[this->HilbertSpaceDimension];
  for (int i=0; i<this->HilbertSpaceDimension; ++i)
    {
      CurrentStateNbrUp = this->StateNbrUp[i];
      LzMaxUp = this->StateLzMaxUp[i] + CurrentStateNbrUp;
      if (CurrentStateNbrUp!=0) --LzMaxUp;
      LzMaxDown = this->StateLzMaxDown[i] + NbrBosons - CurrentStateNbrUp;
      if (CurrentStateNbrUp!=NbrBosons) --LzMaxDown;

      this->StateInfo[i]=(LzMaxUp<<20) | (LzMaxDown<<10) | this->StateNbrUp[i];
      this->BosonToFermion(this->StateDescriptionUp[i], this->StateDescriptionDown[i], LzMaxUp, LzMaxDown,
			   this->StateDescription[i], (this->StateNbrUp[i]));
#ifdef TESTING
      cout << i << " ";
      this->PrintState(cout, StateDescription[i])<<" "<<this->StateDescriptionUp[i]<<" "<< this->StateDescriptionDown[i] << endl;
#endif
      if (this->StateDescriptionUp[i]!=LastUp)
	{
	  LastUp=this->StateDescriptionUp[i];
	  TmpBlockStart[NbrStatesUp]=i;
	  ++NbrStatesUp;
	}
    }
  for (int i=0; i<this->HilbertSpaceDimension; ++i)
    delete [] this->StateDescription[i];
  delete [] this->StateDescription;
  this->StateDescription=NULL;
  delete [] this->StateLzMaxDown;  
  this->StateLzMaxDown=NULL;

  this->UpSpinLookUpTable = new unsigned[NbrStatesUp];  
  this->DownSpinLookUpTable = new unsigned[NbrStatesUp];
  this->NbrTensoredElements = NbrStatesUp;
  this->TensoredElements = new unsigned long[this->NbrTensoredElements];
  this->TensoredLzMax = new int[this->NbrTensoredElements];
  for (unsigned i=0; i<NbrStatesUp; ++i)
    {
      this->UpSpinLookUpTable[i]= TmpBlockStart[i];
      this->DownSpinLookUpTable[i] = 0;
      this->TensoredElements[i] = this->StateDescriptionUp[TmpBlockStart[i]];
      this->TensoredLzMax[i] = this->StateLzMaxUp[TmpBlockStart[i]]+this->StateNbrUp[TmpBlockStart[i]];
      if (this->TensoredLzMax[i]>0) --this->TensoredLzMax[i];
    }
  delete [] this->StateLzMaxUp;
  delete [] this->StateNbrUp;
  this->StateLzMaxUp=NULL;
  this->StateNbrUp=NULL;
  delete [] TmpBlockStart;

  this->ShellSortAux(NbrStatesUp, this->TensoredElements, this->UpSpinLookUpTable, this->TensoredLzMax);

#ifdef TESTING
  for (unsigned i=0; i<NbrStatesUp; ++i)
    cout << "T["<<i<<"]="<<this->TensoredElements[i]<<" "<<std::hex<<this->TensoredElements[i]<<std::dec<<" "<<TensoredLzMax[i]<<endl;
#endif
  
  //exit(1);
  // generate fast lookup table for TensoredIndex
  // evaluate look-up table size
  int PolarizedMaxBits = NbrBosons + LzMax;
  memory /= (sizeof(int*) * PolarizedMaxBits);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > PolarizedMaxBits)
    this->MaximumLookUpShift = PolarizedMaxBits;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [PolarizedMaxBits];
  this->LookUpTableShift = new int [PolarizedMaxBits];
  for (int i = 0; i < PolarizedMaxBits; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentLzMax = this->TensoredLzMax[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentLzMax];
  if (CurrentLzMax < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentLzMax] = 0;
  else
    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentLzMax];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->TensoredElements[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (unsigned i = 0; i < this->NbrTensoredElements; ++i)
    {
      if (CurrentLzMax != this->TensoredLzMax[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  /*
	  cout << "CurrentLzMax = "<<CurrentLzMax<<endl;
	  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
	    cout << TmpLookUpTable[j] << " ";
	  cout << endl << "-------------------------------------------" << endl;
	  */
 	  CurrentLzMax = this->TensoredLzMax[i];
	  TmpLookUpTable = this->LookUpTable[CurrentLzMax];
	  if (CurrentLzMax < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentLzMax] = 0;
	  else
	    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentLzMax];
	  TmpLookUpTableValue = this->TensoredElements[i] >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->TensoredElements[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
//	      CurrentLookUpTableValue = TmpLookUpTableValue;
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->NbrTensoredElements - 1;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->NbrTensoredElements - 1;
  /*
  cout << "CurrentLzMax = 0"<<endl;
  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
    cout << TmpLookUpTable[j] << " ";
  cout << endl << "-------------------------------------------" << endl;
  */

  // deduce final lookup structure
  unsigned Index;
  int DownElementLzMax;
  // int LzOffset;
  LastUp = ~0x0ul;
  unsigned SectorCount=0;
  for (int i=0; i<this->HilbertSpaceDimension; ++i)
    {
      DownElementLzMax = (this->StateInfo[i] >> 10) & 0x03ffu;
      //      cout << "Call search "<< std::hex << this->StateDescriptionDown[i] <<std::dec<< " with lzmax=" << DownElementLzMax << endl;
      Index = this->FindTensoredIndex(this->StateDescriptionDown[i], DownElementLzMax);
#ifdef TESTING
      if (TensoredElements[Index]==this->StateDescriptionDown[i])
	{
	  cout << "Element "<<TensoredElements[Index]<<" found"<<endl;
	}
      else
	cout << "Error: Element "<<this->StateDescriptionDown[i]<<" NOT found"<<endl;
#endif
      if (LastUp!=this->StateDescriptionUp[i])
	{
	  LastUp=this->StateDescriptionUp[i];
	  // this->DownSpinLookUpTable[Index]=0;
	  SectorCount=1;
	}
      else
	{
#ifdef TESTING
	  if ((this->DownSpinLookUpTable[Index]!=0) && (this->DownSpinLookUpTable[Index]!=SectorCount))
	    {
	      cout << hex << this->StateDescriptionDown[i] << " vs " << TensoredElements[Index]<<std::dec<<endl;
	      cout << "Error: space does not seem to factorize well at state "<<i<<" - new index "<<SectorCount<<" old: "<<
		this->DownSpinLookUpTable[Index]<<endl;
	    }
#endif
	  this->DownSpinLookUpTable[Index]=SectorCount;
	  ++SectorCount;
	}
      //++this->DownSpinLookUpTable[Index];
    }
#ifdef TESTING
  unsigned Info;
  for (unsigned i=0; i<(unsigned)this->HilbertSpaceDimension; ++i)
    {
      Info = this->StateInfo[i];
      CurrentStateNbrUp = Info&0x03ff;
      LzMaxUp = (Info>>20); // +CurrentStateNbrUp-(CurrentStateNbrUp!=0);
      LzMaxDown = ((Info>>10)&0x03ff); // + NbrBosons - CurrentStateNbrUp - (CurrentStateNbrUp!=NbrBosons);
      Index = this->FindStateIndex(this->StateDescriptionUp[i], this->StateDescriptionDown[i], LzMaxUp, LzMaxDown);
      if (Index != i)
	cout << "Error in state "<<i<<" base: " << this->UpSpinLookUpTable[this->FindTensoredIndex(StateDescriptionUp[i], LzMaxUp)]
	     << " Offset: "<<this->DownSpinLookUpTable[this->FindTensoredIndex(StateDescriptionDown[i], LzMaxDown)] << endl;
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
// return value = Hilbert space dimension

int BosonOnSphereWithSpinAllSz::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, 2*lzMax+1, (totalLz + lzMax * nbrBosons) >> 1,0);
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// szParity = parity of the total Sz : (-1)^Sz == 1?
// return value = Hilbert space dimension

int BosonOnSphereWithSpinAllSz::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int szParity)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, 2*lzMax+1, (totalLz + lzMax * nbrBosons) >> 1, 0, ((nbrBosons>>1) + szParity)&0x1, 0);
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons remaining to be placed
// lzMax = momentum maximum value for a boson in the state 
// currentLzMax = momentum maximum value for bosons that are still to be placed
//                (counting from 0 to lzMax for down and lzMax+1 to 2LzMax+1 for up)
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int BosonOnSphereWithSpinAllSz::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int currentLzMax, int totalLz, int level)
{
  int CurrentMomentum = currentLzMax%(this->NbrLzValue);
  int CurrentSpinUp = currentLzMax/(this->NbrLzValue);
  //for (int i=0; i<level; ++i) cout << "  ";
  //cout << "SEV(n="<<nbrBosons<<", lz="<<currentLzMax<<", "<<totalLz<<") - Lz="<<CurrentMomentum<<endl;
  if ((nbrBosons < 0) || (totalLz < 0) || (currentLzMax < 0) || ((CurrentSpinUp==0)&&((nbrBosons * CurrentMomentum) < totalLz)))
    {
//       for (int i=0; i<level; ++i) cout << "  ";
//       cout << "add 0"<<endl;
      return 0;
    }
  if (nbrBosons == 0)
    {
      if (totalLz == 0)
	{
// 	  for (int i=0; i<level; ++i) cout << "  ";
// 	  cout << "add 1"<<endl;
	  return 1;
	}
      else
	{
// 	  for (int i=0; i<level; ++i) cout << "  ";	  
// 	  cout << "add 0"<<endl;
	  return 0;
	}
    }
  if ((CurrentSpinUp==0)&&(nbrBosons == 1))
    {
      if (currentLzMax >= totalLz)
	{
// 	  for (int i=0; i<level; ++i) cout << "  ";
// 	  cout << "add 1"<<endl;
	  return 1;
	}
      else return 0;
    }
  if ((CurrentSpinUp==0)&&((nbrBosons * CurrentMomentum) == totalLz))
    {
//       for (int i=0; i<level; ++i) cout << "  ";
//       cout << "add 1"<<endl;
      return 1;
    }
  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpDim = 0;
  for (int ToPlace=nbrBosons; ToPlace>=0; --ToPlace)
    TmpDim += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons-ToPlace, ReducedCurrentLzMax, totalLz-ToPlace*CurrentMomentum, level+1);
//   for (int i=0; i<level; ++i) cout << "  ";
//   cout << "add: "<<TmpDim<<endl;
  return TmpDim;
}


// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons remaining to be placed
// lzMax = momentum maximum value for a boson in the state 
// currentLzMax = momentum maximum value for bosons that are still to be placed
//                (counting from 0 to lzMax for down and lzMax+1 to 2LzMax+1 for up)
// totalLz = momentum total value
// totalNbrUp = number of up-particles placed so far
// szParity = target parity of the total number of up particles
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int BosonOnSphereWithSpinAllSz::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int currentLzMax, int totalLz, int totalNbrUp, int szParity, int level)
{
  int CurrentMomentum = currentLzMax%(this->NbrLzValue);
  int CurrentSpinUp = currentLzMax/(this->NbrLzValue);
  //for (int i=0; i<level; ++i) cout << "  ";
  //cout << "SEV(n="<<nbrBosons<<", lz="<<currentLzMax<<", "<<totalLz<<") - Lz="<<CurrentMomentum<<endl;
  if ((nbrBosons < 0) || (totalLz < 0) || (currentLzMax < 0) || ((CurrentSpinUp==0)&&((nbrBosons * CurrentMomentum) < totalLz)))
    {
//       for (int i=0; i<level; ++i) cout << "  ";
//       cout << "add 0"<<endl;
      return 0;
    }
  if (nbrBosons == 0)
    {
      if (totalLz == 0)
	{
// 	  for (int i=0; i<level; ++i) cout << "  ";
// 	  cout << "add 1"<<endl;
	  if ((totalNbrUp & 0x1) == szParity)
	    return 1;
	  else
	    return 0;
	}
      else
	{
// 	  for (int i=0; i<level; ++i) cout << "  ";	  
// 	  cout << "add 0"<<endl;
	  return 0;
	}
    }
  if ((CurrentSpinUp==0)&&(nbrBosons == 1))
    {
      if (currentLzMax >= totalLz)
	{
// 	  for (int i=0; i<level; ++i) cout << "  ";
// 	  cout << "add 1"<<endl;
	  if ((totalNbrUp & 0x1) == szParity)
	    return 1;
	  else
	    return 0;
	}
      else return 0;
    }
  if ((CurrentSpinUp==0)&&((nbrBosons * CurrentMomentum) == totalLz))
    {
//       for (int i=0; i<level; ++i) cout << "  ";
//       cout << "add 1"<<endl;
      if ((totalNbrUp & 0x1) == szParity)
	return 1;
    }
  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpDim = 0;
  for (int ToPlace=nbrBosons; ToPlace>=0; --ToPlace)
    TmpDim += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons-ToPlace, ReducedCurrentLzMax, totalLz-ToPlace*CurrentMomentum, totalNbrUp+ToPlace*CurrentSpinUp, szParity, level+1);
//   for (int i=0; i<level; ++i) cout << "  ";
//   cout << "add: "<<TmpDim<<endl;
  return TmpDim;
}


// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex BosonOnSphereWithSpinAllSz::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
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

Complex BosonOnSphereWithSpinAllSz::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
							      int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, this->HilbertSpaceDimension);
}


// set a restriction to evaluate the wavefunction on a given subspace of fixed sz
// twoSz = subspace to restrict to
// restriction = flag whether restriction shall be set (true) or deleted (false)
bool BosonOnSphereWithSpinAllSz::WaveFunctionSubSpace(int twoSz, bool restriction)
{
  if (restriction == false)
    {
      this->SubspaceRestriction = false;
      this->SubspaceSz = 0;
      return true;
    }
  
  if (((twoSz & 1) != ( this->NbrBosons &1)) || (twoSz > this->NbrBosons) || (twoSz < -this->NbrBosons))
    return false;
  this->SubspaceRestriction = true;
  this->SubspaceSz = twoSz;
  return true;
}

// get weight of wavefunction in current subspace
// state = vector to be considered
//
// return = weight
double BosonOnSphereWithSpinAllSz::GetWeightInSubSpace(RealVector& state)
{
  int NbrBosonsUp = (this->NbrBosons + this->SubspaceSz) >> 1;
  double Weight = 0.0;
  for (int k = 0; k < this->HilbertSpaceDimension; ++k)
    if ((StateInfo[k]&0x03ffu) == (unsigned)NbrBosonsUp)
      Weight+=state[k]*state[k];
  return Weight;
}


// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex BosonOnSphereWithSpinAllSz::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
					     int firstComponent, int nbrComponent)
{
  if (SubspaceRestriction==false)
    {
      cout << "Wavefunction evaluation in BosonOnSphereWithSpinAllSz defined only for subspaces of fixed sz."<<endl;
      cout << "Indicate a subspace first, using BosonOnSphereWithSpinAllSz::WaveFunctionSubSpace(int twoSz)"<<endl;
    }
  int NbrBosonsUp = (this->NbrBosons + this->SubspaceSz) >> 1;
  int NbrBosonsDown = (this->NbrBosons - this->SubspaceSz) >> 1;
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
  unsigned TemporaryStateNbrUp;
  PermDown.EvaluateFastPermanentPrecalculationArray(ChangeBitDown, ChangeBitSignDown);  
  for (int k = firstComponent; k < LastComponent; ++k)
    if ((StateInfo[k]&0x03ffu) == (unsigned)NbrBosonsUp)
      {
      TmpComponent = state[k];
      if (TmpComponent!=0.0)
	{
	  TmpFactor = Factors[NbrBosonsUp];
	  Pos = 0;
	  Lz = 0;
	  this->FermionToBoson(this->StateDescriptionUp[k], this->StateDescriptionDown[k], this->StateInfo[k],
			       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
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

Complex BosonOnSphereWithSpinAllSz::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
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
  unsigned TemporaryStateNbrUp;
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
			       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
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
			       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
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

void BosonOnSphereWithSpinAllSz::InitializeWaveFunctionEvaluation (bool timeCoherence)
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

RealSymmetricMatrix  BosonOnSphereWithSpinAllSz::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState)
{  
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
			       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
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
      BosonOnSphereWithSpinAllSz TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
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
  BosonOnSphereWithSpinAllSz TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  
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
				   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);

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
				   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
	      int TmpLzMax = subsytemSize - 1;
	      while (TemporaryState[TmpLzMax] == 0) 
		--TmpLzMax;
	      TmpStatePosition[Pos] = TmpDestinationHilbertSpace.FindStateIndex(TemporaryState, TmpLzMax);
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
}




// Project the state from the tunneling space (all Sz's)
// to the space with the fixed projection of Sz (given by SzValue)
//
// state = state that needs to be projected
// su2Space = the subspace onto which the projection is carried out
// SzValue = the desired value of Sz

RealVector BosonOnSphereWithSpinAllSz::ForgeSU2FromTunneling(RealVector& state, BosonOnSphereWithSpin& su2Space, int SzValue)
{
  RealVector FinalState(su2Space.GetHilbertSpaceDimension(), true);
  if ((this->NbrBosons &1)!=(SzValue&1))
    {
      cout << "NbrBosons and SzValue need to have the same parity"<<endl;
      return FinalState;
    }
  unsigned TargetNbrUp = (SzValue+this->NbrBosons)>>1;
  int counter=0;
  double Weight=0.0;
  int TotalLzMax, CurrentLzMaxUp, CurrentLzMaxDown;
  //int LzSzMax;
  unsigned TemporaryStateNbrUp;
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      if ((StateInfo[j]&0x03ffu) == TargetNbrUp)
	{
	  this->FermionToBoson(this->StateDescriptionUp[j], this->StateDescriptionDown[j], this->StateInfo[j],
			       TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);
	  // changed interface for BosonOnSphereWithSpin - no LzSzMax value required for FindStateIndex!
	  //LzSzMax = std::max(2*CurrentLzMaxUp+1,2*CurrentLzMaxDown);
	  TotalLzMax = std::max(CurrentLzMaxUp,CurrentLzMaxDown);
	  for (int i=TotalLzMax+1; i<NbrLzValue; ++i)
	    TemporaryState[i]=0;
	  ++counter;
	  FinalState[su2Space.FindStateIndex(TemporaryState/*, LzSzMax */)] += state[j];
	  Weight+=state[j]*state[j];
	}
    }
  cout<<"Nbr of stored components = "<<counter<<endl;
  cout<<"Weight = "<<Weight<<endl;
  if (fabs(Weight)>1e-12)
    FinalState /= FinalState.Norm();
  return FinalState;  

}

// Project the state from the tunneling space (all Sz's)
// to the U(1) space (u1Space)
//
// state = state that needs to be projected
// u1Space = the subspace onto which the projection is carried out
RealVector BosonOnSphereWithSpinAllSz::ForgeU1FromTunneling(RealVector& state, BosonOnSphere& u1Space)
{
  int Dim2=u1Space.GetHilbertSpaceDimension();
  RealVector FinalState(u1Space.GetHilbertSpaceDimension(), true);
  int Rejected=0, Index, LzMax, CurrentLzMaxUp, CurrentLzMaxDown;
  unsigned TemporaryStateNbrUp;
  int *TmpState = new int[NbrLzValue];
  for (int j = 0; j < this->HilbertSpaceDimension; ++j)    
    {
      this->FermionToBoson(this->StateDescriptionUp[j], this->StateDescriptionDown[j], this->StateInfo[j],
			   TemporaryState, CurrentLzMaxUp, CurrentLzMaxDown, TemporaryStateNbrUp);

      for (int i=0; i<NbrLzValue; ++i)
	TmpState[i] = (this->TemporaryState[i]&0x03ff)+(this->TemporaryState[i]>>16);
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
double BosonOnSphereWithSpinAllSz::MeanSxValue(RealVector& state)
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
double BosonOnSphereWithSpinAllSz::MeanSzValue(RealVector& state)
{
  double Result = 0.0;
  int Dim = this->HilbertSpaceDimension;

  for (int i = 0; i < Dim; ++i)    
    Result+=state[i]*state[i]*0.5*(2*this->StateInfo[i]-this->NbrBosons);

  return Result;
}
