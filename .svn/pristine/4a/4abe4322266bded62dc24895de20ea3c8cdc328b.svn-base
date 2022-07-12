////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//      altenate version of the class for bosons on sphere with SU(2) spin    //
//                                                                            //
//                        last modification : 26/01/2012                      //
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
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinAllLz.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinTwoLandauLevels.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/Permutations.h"
#include "GeneralTools/Endian.h"
#include "Architecture/ArchitectureOperation/FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <map>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
// 

BosonOnSphereWithSU2Spin::BosonOnSphereWithSU2Spin ()
{
}

// basic constructor without any constraint on Sz
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson

BosonOnSphereWithSU2Spin::BosonOnSphereWithSU2Spin (int nbrBosons, int totalLz, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;

  this->NbrBosonsUp = 0;
  this->NbrBosonsDown = 0;
  this->NUpLzMax = this->LzMax + this->NbrBosons - 1;
  this->NDownLzMax = this->LzMax + this->NbrBosons - 1;
  this->FermionicLzMax = this->NUpLzMax;
  if (this->NDownLzMax > this->FermionicLzMax)
    this->FermionicLzMax = this->NDownLzMax;
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1);
  
  this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 
						       this->LzMax + this->NbrBosons + 1, this->LzMax + this->NbrBosons + 1, 0l);

  for (long i = 0; i < TmpHilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescriptionUp[i];
      int TmpNbrUps = this->ComputeNbrParticles(this->StateDescriptionUp[i]);
      this->StateDescriptionUp[i] >>= this->NbrBosons - TmpNbrUps; 
      this->StateDescriptionDown[i] >>= TmpNbrUps; 
    }
  SortDoubleElementArrayDownOrdering<unsigned long>(this->StateDescriptionUp, this->StateDescriptionDown, TmpHilbertSpaceDimension);
  
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2Spin!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->TargetSpace = this;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  



  this->GenerateLookUpTable(10000000);

#ifdef __DEBUG__
   long UsedMemory = 0l;
   UsedMemory += this->LargeHilbertSpaceDimension * (2l * sizeof(unsigned long));
   cout << "memory requested for Hilbert space = ";
   if (UsedMemory >= 1024l)
    if (UsedMemory >= 1048576l)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// totalSpin = twice the total spin value
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU2Spin::BosonOnSphereWithSU2Spin (int nbrBosons, int totalLz, int lzMax, int totalSpin, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;

  this->NbrBosonsUp = this->NbrBosons + this->TotalSpin;
  this->NbrBosonsDown = this->NbrBosons - this->TotalSpin;
  if ((this->NbrBosonsUp < 0) || ((this->NbrBosonsUp & 0x1) != 0) ||
      (this->NbrBosonsDown < 0) || ((this->NbrBosonsDown & 0x1) != 0))
    this->LargeHilbertSpaceDimension = 0l;
  else
    {
      this->NbrBosonsUp >>= 1;
      this->NbrBosonsDown >>= 1;
      this->NUpLzMax = this->LzMax + this->NbrBosonsUp - 1;
      this->NDownLzMax = this->LzMax + this->NbrBosonsDown - 1;
      this->FermionicLzMax = this->NUpLzMax;
      if (this->NDownLzMax > this->FermionicLzMax)
	this->FermionicLzMax = this->NDownLzMax;
      this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, this->NbrBosonsUp, this->NbrBosonsDown);
    }
  this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 
						       this->NbrBosonsUp, this->NbrBosonsDown, 0l);
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2Spin!" << endl;
      exit(1);
    }
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

  this->TargetSpace = this;

  this->GenerateLookUpTable(memory);
  //   for (int i = 0; i < this->HilbertSpaceDimension; ++i)	
  //      {
  //        cout << i << " : ";
  //        this->PrintState(cout, i);
  //        cout << this->FindStateIndex(this->StateDescriptionUp[i], this->StateDescriptionDown[i]);
  //        cout << endl;
//        unsigned long Tmp1;
//        unsigned long Tmp2;
//        this->FermionToBoson(this->StateDescriptionUp[i], this->NUpLzMax, this->TemporaryStateUp);
//        this->FermionToBoson(this->StateDescriptionDown[i], this->NDownLzMax, this->TemporaryStateDown);
//        this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, Tmp1, Tmp2);
//        cout << this->StateDescriptionUp[i] << " " << Tmp1 << " & " << this->StateDescriptionDown[i] << " " << Tmp2 << endl;
//        for (int j = 0; j < this->NbrLzValue; ++j)
// 	 {
// 	   cout << " (" << this->TemporaryStateUp[j] << "," << this->TemporaryStateDown[j] << ")";
// 	 }
//        cout << endl;
//     }
  
#ifdef __DEBUG__
  long UsedMemory = 0;
  UsedMemory += this->LargeHilbertSpaceDimension * (2l * sizeof(unsigned long));
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024l)
    if (UsedMemory >= 1048576l)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
}

// constructor from a binary file that describes the Hilbert space
// 
// fileName = name of the binary file
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU2Spin::BosonOnSphereWithSU2Spin (char* fileName, unsigned long memory)
{
  this->ReadHilbertSpace(fileName);

  this->NbrLzValue = this->LzMax + 1;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->Flag.Initialize();
 
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;

  this->NbrBosonsUp = this->NbrBosons + this->TotalSpin;
  this->NbrBosonsDown = this->NbrBosons - this->TotalSpin;
  this->NbrBosonsUp >>= 1;
  this->NbrBosonsDown >>= 1;
  this->NUpLzMax = this->LzMax + this->NbrBosonsUp - 1;
  this->NDownLzMax = this->LzMax + this->NbrBosonsDown - 1;
  this->FermionicLzMax = this->NUpLzMax;
  if (this->NDownLzMax > this->FermionicLzMax)
    this->FermionicLzMax = this->NDownLzMax;
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  this->TargetSpace = this;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += this->LargeHilbertSpaceDimension * (4l * sizeof(unsigned long));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024l)
	if (UsedMemory >= 1048576l)
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

BosonOnSphereWithSU2Spin::BosonOnSphereWithSU2Spin(const BosonOnSphereWithSU2Spin& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->TotalSpin = bosons.TotalSpin;
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
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

BosonOnSphereWithSU2Spin::~BosonOnSphereWithSU2Spin ()
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
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereWithSU2Spin& BosonOnSphereWithSU2Spin::operator = (const BosonOnSphereWithSU2Spin& bosons)
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
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
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
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
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

AbstractHilbertSpace* BosonOnSphereWithSU2Spin::Clone()
{
  return new BosonOnSphereWithSU2Spin(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereWithSU2Spin::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereWithSU2Spin::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereWithSU2Spin::ExtractSubspace (AbstractQuantumNumber& q, 
								 SubspaceSpaceConverter& converter)
{
  return 0;
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnSphereWithSU2Spin::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
  this->TargetSpace = (BosonOnSphereWithSU2Spin*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int BosonOnSphereWithSU2Spin::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// apply a^+_m_u a_m_u operator to a given state  (only spin up isospin plus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_um a_m_um

double  BosonOnSphereWithSU2Spin::AduAu (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  return (double) (this->TemporaryStateUp[m]);  
}

// apply a^+_m_d a_m_d operator to a given state (only spin down isospin plus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_dm a_m_dm

double BosonOnSphereWithSU2Spin::AddAd (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  return (double) (this->TemporaryStateDown[m]);  
}

// find state index
//
// stateDescriptionUp = unsigned integer describing the fermionic state for type up particles
// stateDescriptionDown = unsigned integer describing the fermionic state for type down particles
// return value = corresponding index

int BosonOnSphereWithSU2Spin::FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
{
  int PosMin = 0;
  int PosMax = this->NbrUniqueStateDescriptionUp - 1;
  if ((stateDescriptionUp > this->UniqueStateDescriptionUp[PosMin]) || (stateDescriptionUp < this->UniqueStateDescriptionUp[PosMax]))
    {
      return this->HilbertSpaceDimension;
    }
  int PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->UniqueStateDescriptionUp[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescriptionUp))
    {
       if (CurrentState > stateDescriptionUp)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->UniqueStateDescriptionUp[PosMid];
    }
  if (CurrentState != stateDescriptionUp)
    {
      PosMid = PosMax;
      if (this->UniqueStateDescriptionUp[PosMid] != stateDescriptionUp)
	return this->HilbertSpaceDimension;
     }

  PosMin = this->FirstIndexUniqueStateDescriptionUp[PosMid];
  PosMax = PosMin + this->UniqueStateDescriptionSubArraySizeUp[PosMid] - 1;
  PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->StateDescriptionDown[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescriptionDown))
    {
       if (CurrentState > stateDescriptionDown)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->StateDescriptionDown[PosMid];
    }
  if (CurrentState == stateDescriptionDown)
    return PosMid;
  if ( this->StateDescriptionDown[PosMax] != stateDescriptionDown)
    return this->HilbertSpaceDimension;
  else
    return PosMax;
}



// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSU2Spin::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUp[state], this->StateDescriptionDown[state],
		       this->TemporaryStateUp, this->TemporaryStateDown); 

  unsigned long Tmp;
  Str << " | ";
  for (int i = this->LzMax; i >=0 ; --i)
    {
      Str << "(" << this->TemporaryStateUp[i] << "," << this->TemporaryStateDown[i] << ") | ";
    }
  return Str;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSU2Spin::PrintStateMonomial (ostream& Str, long state)
{
  unsigned long* TmpMonomialUp = new unsigned long [this->NbrBosonsUp];
  unsigned long* TmpMonomialDown = new unsigned long [this->NbrBosonsDown];
  this->ConvertToMonomial(this->StateDescriptionUp[state], this->StateDescriptionDown[state], TmpMonomialUp, TmpMonomialDown);
  Str << "[";
  for (int i = 0; i < (this->NbrBosonsUp - 1); ++i)
    Str << TmpMonomialUp[i] << ",";
  Str << TmpMonomialUp[this->NbrBosonsUp - 1];
  Str << "][";
  for (int i = 0; i < (this->NbrBosonsDown - 1); ++i)
    Str << TmpMonomialDown[i] << ",";
  Str << TmpMonomialDown[this->NbrBosonsDown - 1];
  Str << "]";
  delete[] TmpMonomialUp;
  delete[] TmpMonomialDown;
  return Str;
}
 
// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMaxUp = momentum maximum value for a boson in the state with up
// lzMaxDown = momentum maximum value for a boson in the state with down
// totalLz = momentum total value
// nbrNUp = number of particles with quantum up
// nbrNDown = number of particles with quantum down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSU2Spin::GenerateStates(int nbrBosons, int lzMaxUp, int lzMaxDown, int totalLz, 
					      int nbrNUp, int nbrNDown, long pos)
{
  if ((nbrBosons < 0) || (totalLz < 0) || (nbrNUp < 0) || (nbrNDown < 0))
    return pos;
  if ((nbrBosons == 0) && (totalLz == 0))
    {
      this->StateDescriptionUp[pos] = 0x0ul;
      this->StateDescriptionDown[pos] = 0x0ul;
      return (pos + 1l);
    }
  if ((lzMaxUp < 0) || (lzMaxDown < 0))
    return pos;

  long TmpPos;
  unsigned long MaskUp;
  unsigned long MaskDown;

  if (nbrNUp == 0)
    {
      for (int l = nbrNDown; l > 0; --l)
	{
	  TmpPos = this->GenerateStates(nbrBosons - l, 0, lzMaxDown - 1, totalLz - (lzMaxDown * l), 
					0, nbrNDown - l, pos); 
	  MaskDown = ((0x1ul << l) - 1ul) << (lzMaxDown + nbrNDown - l);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescriptionDown[pos] |= MaskDown;
	    }
	}
      pos = this->GenerateStates(nbrBosons, 0, lzMaxDown - 1, totalLz, 0, nbrNDown, pos);
      return pos;
    }
  
  TmpPos = this->GenerateStates(nbrBosons - nbrNUp, 0, lzMaxDown, totalLz - (lzMaxUp * nbrNUp), 0, nbrNDown, pos); 
  MaskUp = ((0x1ul << nbrNUp) - 1ul) << lzMaxUp;
  for (; pos < TmpPos; ++pos)
    this->StateDescriptionUp[pos] |= MaskUp;
  for (int i = nbrNUp - 1; i > 0; --i)
    {
      TmpPos = this->GenerateStates(nbrBosons - i, lzMaxUp - 1, lzMaxDown, totalLz - (lzMaxUp * i), nbrNUp - i, nbrNDown, pos); 
      MaskUp = ((0x1ul << i) - 1ul) << (lzMaxUp + nbrNUp - i);
      for (; pos < TmpPos; ++pos)
	{
	  this->StateDescriptionUp[pos] |= MaskUp;
	}
    }
  pos = this->GenerateStates(nbrBosons, lzMaxUp - 1, lzMaxDown, totalLz, nbrNUp, nbrNDown, pos);
  return pos;
};


// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSU2Spin::GenerateStates(int nbrBosons, int lzMax, int totalLz, 
					      int currentFermionicPositionUp, int currentFermionicPositionDown, long pos)
{
  if (nbrBosons == 0)
    {
      if (totalLz == 0)
	{
	  this->StateDescriptionUp[pos] = 0x0ul;
	  this->StateDescriptionDown[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else
	{
	  return pos;
	}
    }
  if (lzMax < 0)
    return pos;

  for (int i = nbrBosons; i >= 0; --i)
    {
      unsigned long MaskUp = ((0x1ul << i) - 0x1ul)  << (currentFermionicPositionUp - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  long TmpPos = this->GenerateStates(nbrBosons - i - j, lzMax - 1, totalLz - (lzMax * (i + j)), 
					     currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, pos); 
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


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereWithSU2Spin::GenerateLookUpTable(unsigned long memory)
{  
  long TmpUniquePartition = 1l;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescriptionUp[i - 1] == this->StateDescriptionUp[i]))
	{
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	++TmpUniquePartition;
    }

  this->NbrUniqueStateDescriptionUp = TmpUniquePartition;
  this->UniqueStateDescriptionUp = new unsigned long [this->NbrUniqueStateDescriptionUp];
  this->UniqueStateDescriptionSubArraySizeUp = new int [this->NbrUniqueStateDescriptionUp];
  this->FirstIndexUniqueStateDescriptionUp = new int [this->NbrUniqueStateDescriptionUp];
  TmpUniquePartition = 0l;
  this->UniqueStateDescriptionUp[0l] = this->StateDescriptionUp[0l];
  this->UniqueStateDescriptionSubArraySizeUp[0] = 1;
  this->FirstIndexUniqueStateDescriptionUp[0] = 0;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescriptionUp[i - 1] == this->StateDescriptionUp[i]))
	{
	  ++this->UniqueStateDescriptionSubArraySizeUp[TmpUniquePartition];
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	{
	  ++TmpUniquePartition;
	  this->UniqueStateDescriptionUp[TmpUniquePartition] = this->StateDescriptionUp[i];
	  this->UniqueStateDescriptionSubArraySizeUp[TmpUniquePartition] = 1; 
	  this->FirstIndexUniqueStateDescriptionUp[TmpUniquePartition] = i;
	}
    }
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrNUp = number of particles with quantum number up
// nbrNDown = number of particles with quantum number down
// return value = Hilbert space dimension

long BosonOnSphereWithSU2Spin::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, 
								    int nbrNUp, int nbrNDown)
{
  if ((nbrBosons < 0) || (totalLz < 0) || (nbrNUp < 0) || (nbrNDown < 0))
    return 0l;
  if ((nbrBosons == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0)
    return 0l;
  if (nbrBosons == 1)
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }
  long Tmp = 0l;
  for (int i = nbrNUp; i >= 0; --i)
    for (int j = nbrNDown; j >= 0; --j)
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - (i + j), lzMax - 1, totalLz - (lzMax * (i + j)), 
							nbrNUp - i, nbrNDown - j);
  return  Tmp;
}

// evaluate Hilbert space dimension without the Sz constraint
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

long BosonOnSphereWithSU2Spin::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons < 0) || (totalLz < 0))
    return 0l;
  if ((nbrBosons == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0)
    return 0l;
  if (nbrBosons == 1)
    {
      if (lzMax >= totalLz)
	return 2l;
      else
	return 0l;
    }
  long Tmp = 0l;
  for (int i = nbrBosons; i >= 0; --i)
    Tmp += (i + 1) * this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - i, lzMax - 1, totalLz - (lzMax * i));
  return Tmp;
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2Spin::AduAu (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  if (this->TemporaryStateUp[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUp[n];
  --this->TemporaryStateUp[n];
  ++this->TemporaryStateUp[m];
  coefficient *= (double) this->TemporaryStateUp[m];
  coefficient = sqrt(coefficient);  
  return this->TargetSpace->FindStateIndex(this->BosonToFermion(this->TemporaryStateUp), this->StateDescriptionDown[index]);  
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2Spin::AduAd (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  
  if (this->TemporaryStateDown[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  /*cout<<"Index = "<<index<<" m = "<<m << " n " <<n<<endl;
  cout<<" Up = ";
  for (int i = 0; i <=  this->LzMax; i++ )
    {
      cout <<this->TemporaryStateUp[i]<<" ";
    }
  cout <<endl;
  cout <<"Down = ";
  for (int i = 0; i <=  this->LzMax; i++ )
    {
      cout <<this->TemporaryStateDown[i]<<" ";
    }
  cout <<endl;*/
  coefficient = (double) this->TemporaryStateDown[n];
  --this->TemporaryStateDown[n];
  ++this->TemporaryStateUp[m];
  coefficient *= (double) this->TemporaryStateUp[m];
  coefficient = sqrt(coefficient); 
  return this->TargetSpace->FindStateIndex(this->BosonToFermion(this->TemporaryStateUp), this->BosonToFermion(this->TemporaryStateDown));
  /*unsigned long Tmp = this->FindStateIndex(this->BosonToFermion(this->TemporaryStateUp), this->BosonToFermion(this->TemporaryStateDown));
  cout <<"Tmp = "<<Tmp<<endl;
  cout<<" Up = ";
  for (int i = 0; i <=  this->LzMax; i++ )
    {
      cout <<this->TemporaryStateUp[i]<<" ";
    }
  cout <<endl;
  cout <<"Down = ";
  for (int i = 0; i <=  this->LzMax; i++ )
    {
      cout <<this->TemporaryStateDown[i]<<" ";
    } 
  cout <<endl;
  cout <<"end"<<endl;
  
  return Tmp;  */
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2Spin::AddAu (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  if (this->TemporaryStateUp[n] == 0)
    { 
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUp[n];
  --this->TemporaryStateUp[n];
  ++this->TemporaryStateDown[m];
  coefficient *= (double) this->TemporaryStateDown[m];
  coefficient = sqrt(coefficient);  
  return this->TargetSpace->FindStateIndex(this->TargetSpace->BosonToFermion(this->TemporaryStateUp), this->TargetSpace->BosonToFermion(this->TemporaryStateDown));  
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2Spin::AddAd (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  if (this->TemporaryStateDown[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDown[n];
  --this->TemporaryStateDown[n];
  ++this->TemporaryStateDown[m];
  coefficient *= (double) this->TemporaryStateDown[m];
  coefficient = sqrt(coefficient);  
  return this->TargetSpace->FindStateIndex(this->StateDescriptionUp[index], this->BosonToFermion(this->TemporaryStateDown));  
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU2Spin::AuAu (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  if ((this->ProdATemporaryStateUp[n1] == 0) || (this->ProdATemporaryStateUp[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateUp[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  double Coefficient = this->ProdATemporaryStateUp[n2];
  --this->ProdATemporaryStateUp[n2];
  Coefficient *= this->ProdATemporaryStateUp[n1];
  --this->ProdATemporaryStateUp[n1];
  return sqrt(Coefficient);
}

// apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU2Spin::AuAd (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateUp[n1] == 0) || (this->ProdATemporaryStateDown[n2] == 0))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryStateDown[n2];
  --this->ProdATemporaryStateDown[n2];
  Coefficient *= this->ProdATemporaryStateUp[n1];
  --this->ProdATemporaryStateUp[n1];
  return sqrt(Coefficient);
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU2Spin::AdAd (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateDown[n1] == 0) || (this->ProdATemporaryStateDown[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateDown[n1] == 1)))    
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  double Coefficient = this->ProdATemporaryStateDown[n2];
  --this->ProdATemporaryStateDown[n2];
  Coefficient *= this->ProdATemporaryStateDown[n1];
  --this->ProdATemporaryStateDown[n1];
  return sqrt(Coefficient);
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU2Spin::AduAdu (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUp, this->TemporaryStateUp, coefficient);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU2Spin::AduAdd (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUp, this->TemporaryStateDown, coefficient);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU2Spin::AddAdd (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDown, this->TemporaryStateDown, coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereWithSU2Spin::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  double Coefficient = 1.0;
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (spinIndices[i] == 0)
	{
	  unsigned long& Tmp = this->ProdATemporaryStateDown[n[i]];
	  if (Tmp == 0x0ul)
	    {
	      return 0.0;
	    }
	  Coefficient *= Tmp;
	  --Tmp;
	}
      else
	{
	  unsigned long& Tmp = this->ProdATemporaryStateUp[n[i]];
	  if (Tmp == 0x0ul)
	    {
	      return 0.0;
	    }
	  Coefficient *= Tmp;
	  --Tmp;
	}
    }
  return sqrt(Coefficient);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2Spin::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      this->TemporaryStateUp[i] = this->ProdATemporaryStateUp[i];
      this->TemporaryStateDown[i] = this->ProdATemporaryStateDown[i];
    }
  coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (spinIndices[i] == 0)
	{
	  unsigned long& Tmp = this->TemporaryStateDown[m[i]];
	  ++Tmp;
	  coefficient *= Tmp;
	}
      else
	{
	  unsigned long& Tmp = this->TemporaryStateUp[m[i]];
	  ++Tmp;
	  coefficient *= Tmp;
	}
    }
  coefficient = sqrt(coefficient);
  return this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool BosonOnSphereWithSU2Spin::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->NbrBosons);
  WriteLittleEndian(File, this->LzMax);
  WriteLittleEndian(File, this->TotalLz);
  WriteLittleEndian(File, this->TotalSpin);
  WriteBlockLittleEndian(File, this->StateDescriptionUp, this->LargeHilbertSpaceDimension);
  WriteBlockLittleEndian(File, this->StateDescriptionDown, this->LargeHilbertSpaceDimension);
  File.close();
  return true;
}

// read Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description is stored
// return value = true if no error occured

bool BosonOnSphereWithSU2Spin::ReadHilbertSpace (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      this->LargeHilbertSpaceDimension = 0l;
      this->HilbertSpaceDimension;
      return false;
    }
  ReadLittleEndian(File, this->LargeHilbertSpaceDimension);
  ReadLittleEndian(File, this->HilbertSpaceDimension);
  ReadLittleEndian(File, this->NbrBosons);
  ReadLittleEndian(File, this->LzMax);
  ReadLittleEndian(File, this->TotalLz);
  ReadLittleEndian(File, this->TotalSpin);
  this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
  ReadBlockLittleEndian(File, this->StateDescriptionUp, this->LargeHilbertSpaceDimension);
  ReadBlockLittleEndian(File, this->StateDescriptionDown, this->LargeHilbertSpaceDimension);
  File.close();
  return true;
}
  
// convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void BosonOnSphereWithSU2Spin::TransformOneBodyBasis(RealVector& initialState, RealVector& targetState, RealMatrix* oneBodyBasis, long firstComponent, long nbrComponents)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU2Indices = new int [this->NbrBosons];
  int* TmpSU2Indices2 = new int [this->NbrBosons];
  double* OccupationCoefficientArray = new double [this->NbrBosons + 1];
  OccupationCoefficientArray[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    OccupationCoefficientArray[i] = OccupationCoefficientArray[i - 1] + 0.5 * log((double) i);
  targetState.ClearVector();
  long LastComponent = firstComponent + nbrComponents;
  if (nbrComponents == 0)
    LastComponent = this->LargeHilbertSpaceDimension;
  for (long i = firstComponent; i < LastComponent; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i],
			   this->TemporaryStateUp, this->TemporaryStateDown); 
      double OccupationCoefficient = 0.0;
      int TmpIndex = 0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  for (int l = 0; l < this->TemporaryStateDown[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 1;
	      ++TmpIndex;	      
	    }
	  for (int l = 0; l < this->TemporaryStateUp[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 0;
	      ++TmpIndex;	      
	    }
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUp[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDown[j]];
	}
      this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpMomentumIndices, TmpSU2Indices, TmpSU2Indices2, oneBodyBasis, OccupationCoefficient, OccupationCoefficientArray);
    }
  delete[] OccupationCoefficientArray;
  delete[] TmpMomentumIndices;
  delete[] TmpSU2Indices;
  delete[] TmpSU2Indices2;
}

// convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void BosonOnSphereWithSU2Spin::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, long firstComponent, long nbrComponents)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU2Indices = new int [this->NbrBosons];
  int* TmpSU2Indices2 = new int [this->NbrBosons];
  double* OccupationCoefficientArray = new double [this->NbrBosons + 1];
  OccupationCoefficientArray[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    OccupationCoefficientArray[i] = OccupationCoefficientArray[i - 1] + 0.5 * log((double) i);
  targetState.ClearVector();
  long LastComponent = firstComponent + nbrComponents;
  if (nbrComponents == 0)
    LastComponent = this->LargeHilbertSpaceDimension;
  for (long i = firstComponent; i < LastComponent; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i],
			   this->TemporaryStateUp, this->TemporaryStateDown); 
      double OccupationCoefficient = 0.0;
      int TmpIndex = 0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  for (int l = 0; l < this->TemporaryStateDown[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 1;
	      ++TmpIndex;	      
	    }
	  for (int l = 0; l < this->TemporaryStateUp[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 0;
	      ++TmpIndex;	      
	    }
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUp[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDown[j]];
	}
      this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpMomentumIndices, TmpSU2Indices, TmpSU2Indices2, oneBodyBasis, OccupationCoefficient, OccupationCoefficientArray);
    }
  delete[] OccupationCoefficientArray;
  delete[] TmpMomentumIndices;
  delete[] TmpSU2Indices;
  delete[] TmpSU2Indices2;
}

// compute the transformation matrix from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// return value = transformation matrix

ComplexMatrix BosonOnSphereWithSU2Spin::TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU2Indices = new int [this->NbrBosons];
  int* TmpSU2Indices2 = new int [this->NbrBosons];
  ComplexMatrix TmpMatrix(this->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  double* OccupationCoefficientArray = new double [this->NbrBosons + 1];
  OccupationCoefficientArray[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    OccupationCoefficientArray[i] = OccupationCoefficientArray[i - 1] + 0.5 * log((double) i);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i], 
			   this->TemporaryStateUp, this->TemporaryStateDown); 
      int TmpIndex = 0;
      double OccupationCoefficient = 0.0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  for (int l = 0; l < this->TemporaryStateDown[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 1;
	      ++TmpIndex;	      
	    }
	  for (int l = 0; l < this->TemporaryStateUp[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 0;
	      ++TmpIndex;	      
	    }
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUp[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDown[j]];
	}
      this->TransformOneBodyBasisRecursive(TmpMatrix[i], 1.0, 0, TmpMomentumIndices, TmpSU2Indices, TmpSU2Indices2, oneBodyBasis,
					   OccupationCoefficient, OccupationCoefficientArray);
    }
  delete[] TmpMomentumIndices;
  delete[] TmpSU2Indices;
  delete[] TmpSU2Indices2;
  delete[] OccupationCoefficientArray;
  return TmpMatrix;
}

// recursive part of the convertion from a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// targetState = vector where the transformed state has to be stored
// coefficient = current coefficient to assign
// position = current particle consider in the n-body state
// momentumIndices = array that gives the momentum partition of the initial n-body state
// initialSU2Indices = array that gives the spin dressing the initial n-body state
// currentSU2Indices = array that gives the spin dressing the current transformed n-body state
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
// occupationCoefficientArray = array that provides 1/2 ln (N!)

void BosonOnSphereWithSU2Spin::TransformOneBodyBasisRecursive(RealVector& targetState, double coefficient,
							      int position, int* momentumIndices, int* initialSU2Indices, int* currentSU2Indices, RealMatrix* oneBodyBasis,
							      double occupationCoefficient, double* occupationCoefficientArray) 
{
  if (position == this->NbrBosons)
    {
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->TemporaryStateUp[i] = 0ul;
	  this->TemporaryStateDown[i] = 0ul;
	}
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  switch (currentSU2Indices[i])
	    {
	    case 0:
	      this->TemporaryStateUp[momentumIndices[i]]++;
	      break;
	    case 1:
	      this->TemporaryStateDown[momentumIndices[i]]++;
	      break;
	    }
	}
      int Index = this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
      if (Index < this->HilbertSpaceDimension)
	{
	  
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateUp[i]];
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateDown[i]];
	    }
	  targetState[Index] += coefficient * exp (occupationCoefficient);
	}
      return;      
    }
  else
    {
      currentSU2Indices[position] = 0;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSU2Indices[position]][1]), position + 1, momentumIndices, initialSU2Indices, currentSU2Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
      currentSU2Indices[position] = 1;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSU2Indices[position]][0]), position + 1, momentumIndices, initialSU2Indices, currentSU2Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
    }
}

// recursive part of the convertion from a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// targetState = vector where the transformed state has to be stored
// coefficient = current coefficient to assign
// position = current particle consider in the n-body state
// momentumIndices = array that gives the momentum partition of the initial n-body state
// initialSU2Indices = array that gives the spin dressing the initial n-body state
// currentSU2Indices = array that gives the spin dressing the current transformed n-body state
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
// occupationCoefficientArray = array that provides 1/2 ln (N!)

void BosonOnSphereWithSU2Spin::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
							      int position, int* momentumIndices, int* initialSU2Indices, int* currentSU2Indices, ComplexMatrix* oneBodyBasis,
							      double occupationCoefficient, double* occupationCoefficientArray) 
{
  if (position == this->NbrBosons)
    {
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->TemporaryStateUp[i] = 0ul;
	  this->TemporaryStateDown[i] = 0ul;
	}
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  switch (currentSU2Indices[i])
	    {
	    case 0:
	      this->TemporaryStateUp[momentumIndices[i]]++;
	      break;
	    case 1:
	      this->TemporaryStateDown[momentumIndices[i]]++;
	      break;
	    }
	}
      int Index = this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
      if (Index < this->HilbertSpaceDimension)
	{
	  
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateUp[i]];
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateDown[i]];
	    }
	  targetState[Index] += coefficient * exp (occupationCoefficient);
	}
      return;      
    }
  else
    {
      currentSU2Indices[position] = 0;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSU2Indices[position]][1]), position + 1, momentumIndices, initialSU2Indices, currentSU2Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
      currentSU2Indices[position] = 1;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSU2Indices[position]][0]), position + 1, momentumIndices, initialSU2Indices, currentSU2Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
    }
}

// compute the projection matrix from the SU(2) Hilbert space to an U(1) Hilbert space
// 
// targetSpace = pointer to the U(1) Hilbert space
// type = type of particles that has to be kept (0 for type up, 1 for type down)
// return value = projection matrix

ComplexMatrix BosonOnSphereWithSU2Spin::TransformationMatrixSU2ToU1(BosonOnSphereShort* targetSpace, int type)
{
  ComplexMatrix TmpMatrix (targetSpace->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  unsigned long* TmpStateDescription;
  unsigned long* TmpStateDescriptionOther;
  if (type == 0)
    {
      TmpStateDescription = this->StateDescriptionUp;
      TmpStateDescriptionOther = this->StateDescriptionDown;
    }
  else
    {
      TmpStateDescription = this->StateDescriptionDown;
      TmpStateDescriptionOther = this->StateDescriptionUp;
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (TmpStateDescriptionOther[i] == 0x0ul)
	{
	  unsigned long TmpState = TmpStateDescription[i];
	  int TmpLzMax = this->FermionicLzMax;
	  while ((TmpState >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  int Index = targetSpace->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (Index < targetSpace->HilbertSpaceDimension)
	    {
	      TmpMatrix[i][Index] = 1.0;
	    }
	}
    }
  return TmpMatrix;
}

// convert state of a SU(2) Hilbert space with fixed Sz to a SU(2) space with all sz sectors
//
// state = state that needs to be projected
// su2space = SU(2) space with fixed sz of the input state
// return value = input state expression in the SU(2) basis

RealVector BosonOnSphereWithSU2Spin::SU2ToSU2AllSz(RealVector& state, ParticleOnSphereWithSpin* su2space)
{
  RealVector TmpVector (this->LargeHilbertSpaceDimension, true);
  BosonOnSphereWithSU2Spin* InputSpace = (BosonOnSphereWithSU2Spin*) su2space;
  unsigned long TmpStateUp;
  unsigned long TmpStateDown;
  for (long i = 0; i < InputSpace->LargeHilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(InputSpace->StateDescriptionUp[i], InputSpace->NUpLzMax, this->ProdATemporaryStateUp);
      this->FermionToBoson(InputSpace->StateDescriptionDown[i], InputSpace->NDownLzMax, this->ProdATemporaryStateDown);
      this->TargetSpace->BosonToFermion(this->ProdATemporaryStateUp, this->ProdATemporaryStateDown, TmpStateUp, TmpStateDown);
      int TmpIndex = this->TargetSpace->FindStateIndex(TmpStateUp, TmpStateDown);
      if (TmpIndex < this->HilbertSpaceDimension)
	{
	  TmpVector[TmpIndex] = state[i];
	}
    }
  return TmpVector;
}

// convert state of a SU(2) Hilbert space with fixed Sz to a SU(2) space with all sz sectors
//
// state = state that needs to be projected
// su2space = SU(2) space with fixed sz of the input state
// return value = input state expression in the SU(2) basis

ComplexVector BosonOnSphereWithSU2Spin::SU2ToSU2AllSz(ComplexVector& state, ParticleOnSphereWithSpin* su2space)
{
  ComplexVector TmpVector (this->LargeHilbertSpaceDimension, true);
  BosonOnSphereWithSU2Spin* InputSpace = (BosonOnSphereWithSU2Spin*) su2space;
  unsigned long TmpStateUp;
  unsigned long TmpStateDown;
  for (long i = 0; i < InputSpace->LargeHilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(InputSpace->StateDescriptionUp[i], InputSpace->NUpLzMax, this->ProdATemporaryStateUp);
      this->FermionToBoson(InputSpace->StateDescriptionDown[i], InputSpace->NDownLzMax, this->ProdATemporaryStateDown);
      this->TargetSpace->BosonToFermion(this->ProdATemporaryStateUp, this->ProdATemporaryStateDown, TmpStateUp, TmpStateDown);
      int TmpIndex = this->TargetSpace->FindStateIndex(TmpStateUp, TmpStateDown);
      if (TmpIndex < this->HilbertSpaceDimension)
	{
	  TmpVector[TmpIndex] = state[i];
	}
    }
  return TmpVector;
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& BosonOnSphereWithSU2Spin::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  unsigned long* TmpMonomialReferenceUp = new unsigned long [this->NbrBosonsUp];
  unsigned long* TmpMonomialReferenceDown = new unsigned long [this->NbrBosonsDown];
  unsigned long* TmpMonomialUp = new unsigned long [this->NbrBosonsUp];
  unsigned long* TmpMonomialDown = new unsigned long [this->NbrBosonsDown];
  double Factor = 1.0 / state[reference];
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  this->ConvertToMonomial(this->StateDescriptionUp[reference], this->StateDescriptionDown[reference], TmpMonomialReferenceUp, TmpMonomialReferenceDown);
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(this->StateDescriptionUp[reference], this->StateDescriptionDown[reference], this->TemporaryStateUp, this->TemporaryStateDown);
  for (int k = 0; k <= this->LzMax; ++k)
    if (this->TemporaryStateUp[k] > 1)
      ReferenceFactorial.FactorialMultiply(this->TemporaryStateUp[k]);
  for (int k = 0; k <= this->LzMax; ++k)
    if (this->TemporaryStateDown[k] > 1)
      ReferenceFactorial.FactorialMultiply(this->TemporaryStateDown[k]);
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->StateDescriptionUp[i], this->StateDescriptionDown[i], TmpMonomialUp, TmpMonomialDown);
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      while ((Index1 < this->NbrBosonsUp) && (Index2 < this->NbrBosonsUp))
	{
	  while ((Index1 < this->NbrBosonsUp) && (TmpMonomialReferenceUp[Index1] > TmpMonomialUp[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReferenceUp[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrBosonsUp) && (Index2 < this->NbrBosonsUp) && (TmpMonomialReferenceUp[Index1] == TmpMonomialUp[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrBosonsUp) && (TmpMonomialReferenceUp[Index1] < TmpMonomialUp[Index2]))
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomialUp[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrBosonsUp)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReferenceUp[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrBosonsUp)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomialUp[Index2]];
	  ++Index2;
	}
      Index1 = 0;
      Index2 = 0;
      while ((Index1 < this->NbrBosonsDown) && (Index2 < this->NbrBosonsDown))
	{
	  while ((Index1 < this->NbrBosonsDown) && (TmpMonomialReferenceDown[Index1] > TmpMonomialDown[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReferenceDown[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrBosonsDown) && (Index2 < this->NbrBosonsDown) && (TmpMonomialReferenceDown[Index1] == TmpMonomialDown[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrBosonsDown) && (TmpMonomialReferenceDown[Index1] < TmpMonomialDown[Index2]))
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomialDown[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrBosonsDown)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReferenceDown[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrBosonsDown)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomialDown[Index2]];
	  ++Index2;
	}
      if (symmetryFactor == true)
	{
	  Factorial = ReferenceFactorial;
	  this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i], this->TemporaryStateUp, this->TemporaryStateDown);
	  for (int k = 0; k <= LzMax; ++k)
	    if (this->TemporaryStateUp[k] > 1)
	      Factorial.FactorialMultiply(this->TemporaryStateUp[k]);
	  for (int k = 0; k <= LzMax; ++k)
	    if (this->TemporaryStateDown[k] > 1)
	      Factorial.FactorialMultiply(this->TemporaryStateDown[k]);
	  Coefficient *= sqrt(Factorial.GetNumericalValue());
	}
      state[i] *= Coefficient;
    }
  delete[] TmpMonomialReferenceUp;
  delete[] TmpMonomialReferenceDown;
  delete[] TmpMonomialUp;
  delete[] TmpMonomialDown;
  return state;
}

// normalize from the cylinder geometry to the Jack normalization
//
// state = reference to the state to unnormalize
// aspect = cylinder aspect ratio
// reference = set which component as to be normalized to 1
// return value = unnormalized state

RealVector& BosonOnSphereWithSU2Spin::NormalizeCylinderToJack(RealVector& state, double aspect, long reference)
{
  long double Pi_L = 3.14159265358979323846264338328L;
  double Perimeter = sqrt(2.0 * M_PI * (double)(this->LzMax + 1) / aspect);
  double KappaSquare = 2.0 * M_PI / Perimeter;
  KappaSquare *= KappaSquare;
  double TmpShift = 0.5 * ((double) this->LzMax);
  long ReferenceIndex = reference;
  double* LogFactorials = new double [this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2; i <= this->NbrBosons; ++i)
    {
      LogFactorials[i] = LogFactorials[i - 1] + log((double) i);
    }


  double ReferenceExpFactor = 0.0;
  this->FermionToBoson(this->StateDescriptionUp[ReferenceIndex], this->StateDescriptionDown[ReferenceIndex], this->TemporaryStateUp, this->TemporaryStateDown);
  double SumMSquare = 0.0;
  double TmpFactorial = 0.0;
  for (int j = 0; j <= this->LzMax; ++j)
    {
      if (this->TemporaryStateUp[j] > 0)
	{
	  SumMSquare -= (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateUp[j]);
	  TmpFactorial += LogFactorials[this->TemporaryStateUp[j]]; 
	}
      if (this->TemporaryStateDown[j] > 0)
	{
	  SumMSquare -= (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateDown[j]);
	  TmpFactorial += LogFactorials[this->TemporaryStateDown[j]]; 
	}
    }  
  ReferenceExpFactor = (KappaSquare * SumMSquare) + TmpFactorial;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i], this->TemporaryStateUp, this->TemporaryStateDown);
      SumMSquare = 0.0;
      TmpFactorial = 0.0;
      for (int j = 0; j <= this->LzMax; ++j)
	{
	  if (this->TemporaryStateUp[j] > 0)
	    {
	      SumMSquare -= (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateUp[j]);
	      TmpFactorial += LogFactorials[this->TemporaryStateUp[j]]; 
	    }
	  if (this->TemporaryStateDown[j] > 0)
	    {
	      SumMSquare -= (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateDown[j]);
	      TmpFactorial += LogFactorials[this->TemporaryStateDown[j]]; 
	    }
	}
      state[i] *= exp(0.5 * ((KappaSquare * SumMSquare) + TmpFactorial - ReferenceExpFactor));      
    }
  state /= state[reference];

  delete[] LogFactorials;
  return state;
}

// normalize Jack with respect to cylinder basis
//
// state = reference to the Jack state to normalize
// aspect = aspect ratio of cylinder
// return value = normalized state

RealVector& BosonOnSphereWithSU2Spin::NormalizeJackToCylinder(RealVector& state, double aspect)
{
  long double Pi_L = 3.14159265358979323846264338328L;
  double Perimeter = sqrt(2.0 * M_PI * (double)(this->LzMax + 1) / aspect);
  double KappaSquare = 2.0 * M_PI / Perimeter;
  KappaSquare *= KappaSquare;
  double TmpShift = 0.5 * ((double) this->LzMax);
  long ReferenceIndex = 0l;
  double Norm = 0.0;//state[0l] * state[0l];
  double* LogFactorials = new double [this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2; i <= this->NbrBosons; ++i)
    {
      LogFactorials[i] = LogFactorials[i - 1] + log((double) i);
    }


  double ReferenceExpFactor = 0.0;
  this->FermionToBoson(this->StateDescriptionUp[ReferenceIndex], this->StateDescriptionDown[ReferenceIndex], this->TemporaryStateUp, this->TemporaryStateDown);
  double SumMSquare = 0.0;
  double TmpFactorial = 0.0;
  for (int j = 0; j <= this->LzMax; ++j)
    {
      if (this->TemporaryStateUp[j] > 0)
	{
	  SumMSquare += (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateUp[j]);
	  TmpFactorial -= LogFactorials[this->TemporaryStateUp[j]]; 
	}
      if (this->TemporaryStateDown[j] > 0)
	{
	  SumMSquare += (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateDown[j]);
	  TmpFactorial -= LogFactorials[this->TemporaryStateDown[j]]; 
	}
    }  
  ReferenceExpFactor = (KappaSquare * SumMSquare) + TmpFactorial;
  for (int i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i], this->TemporaryStateUp, this->TemporaryStateDown);
      SumMSquare = 0.0;
      TmpFactorial = 0.0;
      for (int j = 0; j <= this->LzMax; ++j)
	{
	  if (this->TemporaryStateUp[j] > 0)
	    {
	      SumMSquare += (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateUp[j]);
	      TmpFactorial -= LogFactorials[this->TemporaryStateUp[j]]; 
	    }
	  if (this->TemporaryStateDown[j] > 0)
	    {
	      SumMSquare += (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateDown[j]);
	      TmpFactorial -= LogFactorials[this->TemporaryStateDown[j]]; 
	    }
	}
      state[i] *= exp(0.5 * ((KappaSquare * SumMSquare) + TmpFactorial - ReferenceExpFactor));      
      Norm += state[i] * state[i];
    }
  cout << "Norm= "<< Norm << endl;
  state /= sqrt(Norm);
  delete[] LogFactorials;
  return state;
}

// normalize a state defined on the sphere geometry with respect to cylinder basis
//
// state = reference to the state to normalize
// aspect = aspect ratio of cylinder
// return value = normalized state

RealVector& BosonOnSphereWithSU2Spin::NormalizeSphereToCylinder(RealVector& state, double aspect)
{
  double Perimeter = sqrt(2.0 * M_PI * (double)(this->LzMax + 1) / aspect);
  double KappaSquare = 2.0 * M_PI / Perimeter;
  KappaSquare *= KappaSquare;
  double TmpShift = 0.5 * ((double) this->LzMax);
  long ReferenceIndex = 0l;
  double Norm = 0.0;
  double* LogSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      LogSqrtCoefficients[k] = log(Binomials.GetNumericalCoefficient(this->LzMax, k));
    }
  double* LogFactorials = new double [this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2; i <= this->NbrBosons; ++i)
    {
      LogFactorials[i] = LogFactorials[i - 1] + log((double) i);
    }

  double ReferenceExpFactor = 0.0;
  this->FermionToBoson(this->StateDescriptionUp[ReferenceIndex], this->StateDescriptionDown[ReferenceIndex], this->TemporaryStateUp, this->TemporaryStateDown);
  double SumMSquare = 0.0;
  double TmpSphereFactor = 0.0;
  for (int j = 0; j <= this->LzMax; ++j)
    {
      if (this->TemporaryStateUp[j] > 0)
	{
	  SumMSquare += (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateUp[j]);
	  TmpSphereFactor += this->TemporaryStateUp[j] * LogSqrtCoefficients[j];
	}
      if (this->TemporaryStateDown[j] > 0)
	{
	  SumMSquare += (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateDown[j]);
	  TmpSphereFactor += this->TemporaryStateDown[j] * LogSqrtCoefficients[j];
	}
    }  
  ReferenceExpFactor = (KappaSquare * SumMSquare) + TmpSphereFactor;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i], this->TemporaryStateUp, this->TemporaryStateDown);
      SumMSquare = 0.0;
      TmpSphereFactor = 0.0;
      for (int j = 0; j <= this->LzMax; ++j)
	{
	  if (this->TemporaryStateUp[j] > 0)
	    {
	      SumMSquare += (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateUp[j]);
	      TmpSphereFactor += this->TemporaryStateUp[j] * LogSqrtCoefficients[j];
	    }
	  if (this->TemporaryStateDown[j] > 0)
	    {
	      SumMSquare += (((((double) j) - TmpShift) * (((double) j) - TmpShift)) *  this->TemporaryStateDown[j]);
	      TmpSphereFactor += this->TemporaryStateDown[j] * LogSqrtCoefficients[j];
	    }
	}
      state[i] *= exp(0.5 * ((KappaSquare * SumMSquare) + TmpSphereFactor - ReferenceExpFactor));      
      Norm += state[i] * state[i];
    }
  cout << "Norm= "<< Norm << endl;
  state /= sqrt(Norm);
  return state;
}

// Compute the product of a spinful fermionic state with a Van der Monde determinant 
//
// fermionicState = reference on the spinful fermionic state
// outputVector = reference on the vector where the result will be stored
// fermionicSpace = pointer on the Hilbert Space associated to the spinful fermionic state
// minIndex = first component to compute
// nbrComponents = number of components to compute
// unnormalizedFlag = true if the state should be written in the unnormalized basis
// cylinderFlag = true if the state should be written on the cylinder geometry
// cylinderPerimeter = cylinder perimeter
// architecture = pointer to the architecture

void BosonOnSphereWithSU2Spin::SlaterTimeSpinfulFermionicState(RealVector& fermionicState, RealVector& outputVector, FermionOnSphereWithSpin* fermionicSpace, 
							       int minIndex, int nbrComponents, bool unnormalizedFlag, bool cylinderFlag, double cylinderPerimeter,
							       AbstractArchitecture* architecture)
{
  outputVector.ClearVector();
  RealVector FinalState (this->GetHilbertSpaceDimension());
  unsigned long* SlaterUp = new unsigned long[this->NbrBosons];
  unsigned long* SlaterDown = new unsigned long[this->NbrBosons];
  int MaxIndex = minIndex + nbrComponents;
  FactorialCoefficient Coefficient;	
  double** ThreeOrbitalOverlaps = new double* [this->LzMax + 1];

  if ((this->LzMax + 1) >= this->NbrBosons)
    {
      if (cylinderFlag == false)
	{
	  BinomialCoefficients Binomials(this->LzMax);
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      ThreeOrbitalOverlaps[i] = new double [this->NbrBosons];
	      double TmpFactor1 = log(((double) ((fermionicSpace->LzMax + 1) * this->NbrBosons)) / ((double) (this->LzMax + 1)) / (4.0 * M_PI)) - log(Binomials.GetNumericalCoefficient(this->LzMax, i));
	      for (int j = 0; (j < this->NbrBosons) && (j <= i); ++j)
		{
		  if (unnormalizedFlag == false)
		    {
		      ThreeOrbitalOverlaps[i][j] = 0.5 * (TmpFactor1 + log(Binomials.GetNumericalCoefficient(this->NbrBosons - 1, j)) 
							  + log(Binomials.GetNumericalCoefficient(fermionicSpace->LzMax, i - j)));
		    }
		  else
		    {
		      ThreeOrbitalOverlaps[i][j] = 0.0;
		    }
		}
	    }
	}
      else
	{
	  double TmpShift1 = 0.5 * ((double) this->LzMax);
	  double TmpShift2 = 0.5 * ((double) (this->NbrBosons - 1));
	  double TmpShift3 = 0.5 * ((double) fermionicSpace->LzMax);
	  double SquareMagneticLengthRescaling2 = ((double) this->LzMax) / ((double) (this->NbrBosons - 1));
	  double SquareMagneticLengthRescaling3 = ((double) this->LzMax) / ((double) fermionicSpace->LzMax);
	  double TmpFactor = 1.0 / (1.0 + (1.0 / SquareMagneticLengthRescaling2) + (1.0 / SquareMagneticLengthRescaling3));
	  
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      double Kappa2Factor =  2.0 * M_PI / cylinderPerimeter;
	      Kappa2Factor *= Kappa2Factor;
	      ThreeOrbitalOverlaps[i] = new double [this->NbrBosons];
	      for (int j = 0; (j < this->NbrBosons) && (j <= i); ++j)
		{
		  double TmpK1 = ((double) i) - TmpShift1;
		  double TmpK2 = ((double) j) - TmpShift2;
		  double TmpK3 = ((double) (i - j)) - TmpShift3;
		  ThreeOrbitalOverlaps[i][j] = log(this->ComputeIntegralPhi0Phi0Phi0OnCylinder(TmpK1, 1.0, TmpK2, SquareMagneticLengthRescaling2, TmpK3, SquareMagneticLengthRescaling3, cylinderPerimeter));
// 		  ThreeOrbitalOverlaps[i][j] = -0.5 * Kappa2Factor * ((TmpK1 * TmpK1) + SquareMagneticLengthRescaling2 * (TmpK2 * TmpK2)
// 								      + SquareMagneticLengthRescaling3 * (TmpK3 * TmpK3)
// 								      - TmpFactor * (TmpK1 + TmpK2 + TmpK3) * (TmpK1 + TmpK2 + TmpK3));
		}
	    }
	}
      for (long i = minIndex; i < MaxIndex; i++)
	{
	  if (fermionicState[i] != 0.0)
	    {
	      fermionicSpace->ConvertToMonomial(fermionicSpace->StateDescription[i], SlaterUp, SlaterDown);
	      this->VanDerMondeTimesSlater(SlaterUp, SlaterDown, FinalState, ThreeOrbitalOverlaps);
	      for (int Index = 0; Index < FinalState.GetVectorDimension(); ++Index)
		{
		  if (FinalState[Index] != 0.0)
		    {
		      this->FermionToBoson(this->StateDescriptionUp[Index], this->StateDescriptionDown[Index],
					   this->TemporaryStateUp, this->TemporaryStateDown); 
		      Coefficient.SetToOne();
		      for(int p = 0; p <= this->LzMax; ++p)
			{
			  if (this->TemporaryStateUp[p] > 1)
			    Coefficient.FactorialMultiply(this->TemporaryStateUp[p]);
			  if (this->TemporaryStateDown[p] > 1)
			    Coefficient.FactorialMultiply(this->TemporaryStateDown[p]);
			}		  
 		      if (unnormalizedFlag == false)
 			{
 			  outputVector[Index] += sqrt(Coefficient.GetNumericalValue()) * fermionicState[i] * FinalState[Index];
 			}
 		      else
			{
			  outputVector[Index] += Coefficient.GetNumericalValue() * fermionicState[i] * FinalState[Index];
			}
		    }
		}
	    }
	}
    }
  else
    {
      if (cylinderFlag == false)
	{
	  BinomialCoefficients Binomials(this->NbrBosons);
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      ThreeOrbitalOverlaps[i] = new double [this->NbrBosons];
	    }
	  for (int i = 0; i <= fermionicSpace->LzMax; ++i)
	    {
	      for (int j = 0; j < this->NbrBosons; ++j)
		{
		  if (((j - i) >= 0) && ((j - i) <= this->LzMax))
		    {
		      if (unnormalizedFlag == false)
			{
			  ThreeOrbitalOverlaps[j - i][j] = -0.5 * (log(((double) ((fermionicSpace->LzMax + 1) * this->NbrBosons)) / ((double) (this->LzMax + 1)) / (4.0 * M_PI)) 
								   - log(Binomials.GetNumericalCoefficient(this->LzMax, j - i))
								   + log(Binomials.GetNumericalCoefficient(this->NbrBosons - 1, j)) 
								   + log(Binomials.GetNumericalCoefficient(fermionicSpace->LzMax, i)));
			}
		      else
			{
			  ThreeOrbitalOverlaps[j - i][j] = 0.0;
			}
		    }
		}
	    }
	}
      else
	{
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      ThreeOrbitalOverlaps[i] = new double [this->NbrBosons];
	    }
	  double Kappa2Factor =  2.0 * M_PI / cylinderPerimeter;
	  Kappa2Factor *= Kappa2Factor;
	  double TmpShift1 = 0.5 * ((double) this->LzMax);
	  double TmpShift2 = 0.5 * ((double) (this->NbrBosons - 1));
	  double TmpShift3 = 0.5 * ((double) fermionicSpace->LzMax);
// 	  double SquareMagneticLengthRescaling2 = ((double) this->LzMax) / ((double) (this->NbrBosons - 1));
// 	  double SquareMagneticLengthRescaling3 = ((double) this->LzMax) / ((double) fermionicSpace->LzMax);
// 	  double TmpFactor = 1.0 / (1.0 + (1.0 / SquareMagneticLengthRescaling2) + (1.0 / SquareMagneticLengthRescaling3));
	  
// 	  for (int i = 0; i <= fermionicSpace->LzMax; ++i)
// 	    {
// 	      for (int j = 0; j < this->NbrBosons; ++j)
// 		{
// 		  if (((j - i) >= 0) && ((j - i) <= this->LzMax))
// 		    {
// 		      double TmpK1 = ((double) j - i) - TmpShift1;
// 		      double TmpK2 = ((double) j) - TmpShift2;
// 		      double TmpK3 = (((double) i) - TmpShift3);
// 		      ThreeOrbitalOverlaps[j - i][j] = -0.5 * Kappa2Factor * ((TmpK1 * TmpK1) + SquareMagneticLengthRescaling2 * (TmpK2 * TmpK2)
// 									      + SquareMagneticLengthRescaling3 * (TmpK3 * TmpK3)
// 									      - TmpFactor * (TmpK1 + TmpK2 + TmpK3) * (TmpK1 + TmpK2 + TmpK3));
// 		    }
// 		}
// 	    }
	  
	  Kappa2Factor /= ((double) (fermionicSpace->LzMax + this->NbrBosons - 1)) / ((double) (this->LzMax + 1));
	  double SquareMagneticLengthRescaling2 = ((double) this->LzMax) / ((double) (this->NbrBosons - 1));
	  double SquareMagneticLengthRescaling3 = ((double) this->LzMax) / ((double) fermionicSpace->LzMax);
	  double TmpFactor = 1.0 / (1.0 + (1.0 / SquareMagneticLengthRescaling2) - (1.0 / SquareMagneticLengthRescaling3));
	  cout << TmpFactor << " " << fermionicSpace->LzMax << " " << (this->NbrBosons - 1) << endl;
	  for (int i = 0; i <= fermionicSpace->LzMax; ++i)
	    {
	      for (int j = 0; j < this->NbrBosons; ++j)
		{
		  if (((j - i) >= 0) && ((j - i) <= this->LzMax))
		    {
		      double TmpK1 = ((double) j - i) - TmpShift1;
		      double TmpK2 = ((double) j) - TmpShift2;
		      double TmpK3 = (((double) i) - TmpShift3);
		      ThreeOrbitalOverlaps[j - i][j] = log(this->ComputeIntegralPhi0Phi0Phi0OnCylinder(TmpK2, SquareMagneticLengthRescaling2, TmpK1, 1.0, -TmpK3, SquareMagneticLengthRescaling3, cylinderPerimeter));
// 		      ThreeOrbitalOverlaps[j - i][j] = -0.5 * Kappa2Factor * ((TmpK1 * TmpK1) + SquareMagneticLengthRescaling2 * (TmpK2 * TmpK2)
// 									      + SquareMagneticLengthRescaling3 * (TmpK3 * TmpK3)
// 									      - TmpFactor * (TmpK1 + TmpK2 + TmpK3) * (TmpK1 + TmpK2 + TmpK3));
		    }
		}
	    }
	}

      for (long i = minIndex; i < MaxIndex; i++)
	{
	  if (fermionicState[i] != 0.0)
	    {
	      fermionicSpace->ConvertToMonomial(fermionicSpace->StateDescription[i], SlaterUp, SlaterDown);

	      FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation TmpOperation(this, true, SlaterUp, SlaterDown, ThreeOrbitalOverlaps);
	      TmpOperation.ApplyOperation(architecture);
	      FinalState = TmpOperation.GetDestinationVector();
	      
//	      this->ReverseVanDerMondeTimesSlater(SlaterUp, SlaterDown, FinalState, ThreeOrbitalOverlaps);

// 	      for (int j = 0; j < this->NbrBosons; ++j)
// 		this->ReverseVanDerMondeTimesSlater(SlaterUp, SlaterDown, FinalState, ThreeOrbitalOverlaps, j);

	      for (int Index = 0; Index < FinalState.GetVectorDimension(); ++Index)
		{
		  if (FinalState[Index] != 0.0)
		    {
		      this->FermionToBoson(this->StateDescriptionUp[Index], this->StateDescriptionDown[Index],
					   this->TemporaryStateUp, this->TemporaryStateDown); 
		      Coefficient.SetToOne();
		      for(int p = 0; p <= this->LzMax; ++p)
			{
			  if (this->TemporaryStateUp[p] > 1)
			    Coefficient.FactorialMultiply(this->TemporaryStateUp[p]);
			  if (this->TemporaryStateDown[p] > 1)
			    Coefficient.FactorialMultiply(this->TemporaryStateDown[p]);
			}		  
		      if (unnormalizedFlag == false)
			{
			  outputVector[Index] += sqrt(Coefficient.GetNumericalValue()) * fermionicState[i] * FinalState[Index];
			}
		      else
			{
			  outputVector[Index] += Coefficient.GetNumericalValue() * fermionicState[i] * FinalState[Index];
			}
		    }
		}
	    }
	}
    }
  for (int i = 0; i <= this->LzMax; ++i)
    {
      delete[] ThreeOrbitalOverlaps[i];
    }  
  delete[] ThreeOrbitalOverlaps;
}

// Compute the product of a spinful fermionic state with a Van der Monde determinant 
//
// fermionicState = reference on the spinful fermionic state
// outputVector = reference on the vector where the result will be stored
// fermionicSpace = pointer on the Hilbert Space associated to the spinful fermionic state
// minIndex = first component to compute
// nbrComponents = number of components to compute
// unnormalizedFlag = true if the state should be written in the unnormalized basis
// cylinderFlag = true if the state should be written on the cylinder geometry
// cylinderPerimeter = cylinder perimeter
// architecture = pointer to the architecture

void BosonOnSphereWithSU2Spin::SlaterTimeSpinfulFermionicState(RealVector& fermionicState, RealVector& outputVector, FermionOnSphereWithSpinTwoLandauLevels* fermionicSpace, 
							       int minIndex, int nbrComponents, bool unnormalizedFlag, bool cylinderFlag, double cylinderPerimeter, 
							       AbstractArchitecture* architecture)
{
  outputVector.ClearVector();
  RealVector FinalState (this->GetHilbertSpaceDimension());
  unsigned long* SlaterLLLUp = new unsigned long[this->NbrBosons];
  unsigned long* Slater2LLUp = new unsigned long[this->NbrBosons];
  unsigned long* SlaterLLLDown = new unsigned long[this->NbrBosons];
  unsigned long* Slater2LLDown = new unsigned long[this->NbrBosons];
  int MaxIndex = minIndex + nbrComponents;
  FactorialCoefficient Coefficient;	
  double*** ThreeOrbitalOverlaps = new double** [2];
  int MaxLambdaLevelIndex = 1;
  for (int LambdaLevel = 0; LambdaLevel <= MaxLambdaLevelIndex; ++LambdaLevel)
    {      
      ThreeOrbitalOverlaps[LambdaLevel] = new double* [this->LzMax + 1];
    }

  if ((this->LzMax + 1) >= this->NbrBosons)
    {
      for (int LambdaLevel = 0; LambdaLevel <= MaxLambdaLevelIndex; ++LambdaLevel)
	{
	  int TmpMaxAngularMomentumLambdaLevel = fermionicSpace->LzMax + 2 * (LambdaLevel - MaxLambdaLevelIndex);
	  double TmpPrefactor = (1.0 - (2.0 * ((double) (LambdaLevel & 1)))) * sqrt (((double) this->NbrBosons) * ((double) (TmpMaxAngularMomentumLambdaLevel + 1)) / (4.0 * M_PI * (this->LzMax + 1)));
	  
	  ClebschGordanCoefficients Clebsch(TmpMaxAngularMomentumLambdaLevel, this->NbrBosons - 1);
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      ThreeOrbitalOverlaps[LambdaLevel][i] = new double [this->NbrBosons];
	      for (int j = 0; j < this->NbrBosons; ++j)
		{
		  if (((i - j + LambdaLevel) >= 0) && ((i - j + LambdaLevel) <= TmpMaxAngularMomentumLambdaLevel))
		    {
		      if (unnormalizedFlag == false)
			{
//			  cout << i << " " << j << " :" << (2 * (i - j + LambdaLevel) - TmpMaxAngularMomentumLambdaLevel) << " " << ((2 * j) - (this->NbrBosons - 1)) <<  endl;		      
			  ThreeOrbitalOverlaps[LambdaLevel][i][j] = (TmpPrefactor * Clebsch.GetCoefficient(2 * (i - j + LambdaLevel) - TmpMaxAngularMomentumLambdaLevel, (2 * j) - (this->NbrBosons - 1), this->LzMax) 
								     * Clebsch.GetCoefficient(fermionicSpace->LzMax - 2 * MaxLambdaLevelIndex, this->NbrBosons - 1, this->LzMax));
//			  cout << LambdaLevel << " " << i << " " << j << " : " << ThreeOrbitalOverlaps[LambdaLevel][i][j] << endl;

			}
		      else
			{
			  ThreeOrbitalOverlaps[LambdaLevel][i][j] = 1.0;
			}
		    }
		  else
		    {
		      ThreeOrbitalOverlaps[LambdaLevel][i][j] = 0.0;
		    }		  
		}
	    }
	}
      for (long i = minIndex; i < MaxIndex; i++)
	{
	  if (fermionicState[i] != 0.0)
	    {
	      int TmpNbrBosonsLLLUp;
	      int TmpNbrBosons2LLUp;
	      int TmpNbrBosonsLLLDown;
	      int TmpNbrBosons2LLDown;
	      fermionicSpace->ConvertToMonomial(fermionicSpace->StateDescription[i], SlaterLLLUp, Slater2LLUp, SlaterLLLDown, Slater2LLDown, TmpNbrBosonsLLLUp, TmpNbrBosonsLLLDown);	      
	      this->VanDerMondeTimesSlater(SlaterLLLUp, Slater2LLUp, SlaterLLLDown, Slater2LLDown, TmpNbrBosonsLLLUp, TmpNbrBosonsLLLDown, 
					   FinalState, ThreeOrbitalOverlaps);
	      for (int Index = 0; Index < FinalState.GetVectorDimension(); ++Index)
		{
		  if (FinalState[Index] != 0.0)
		    {
		      this->FermionToBoson(this->StateDescriptionUp[Index], this->StateDescriptionDown[Index],
					   this->TemporaryStateUp, this->TemporaryStateDown); 
		      Coefficient.SetToOne();
		      for(int p = 0; p <= this->LzMax; ++p)
			{
			  if (this->TemporaryStateUp[p] > 1)
			    Coefficient.FactorialMultiply(this->TemporaryStateUp[p]);
			  if (this->TemporaryStateDown[p] > 1)
			    Coefficient.FactorialMultiply(this->TemporaryStateDown[p]);
			}		  
		      if (unnormalizedFlag == false)
			{
			  outputVector[Index] += sqrt(Coefficient.GetNumericalValue()) * fermionicState[i] * FinalState[Index];
			}
		      else
			{
			  outputVector[Index] += Coefficient.GetNumericalValue() * fermionicState[i] * FinalState[Index];
			}
		    }
		}
	    }
	}
    }
  else
    {
      for (int LambdaLevel = 0; LambdaLevel <= MaxLambdaLevelIndex; ++LambdaLevel)
	{
	  int TmpMaxAngularMomentumLambdaLevel = fermionicSpace->LzMax + 2 * (LambdaLevel - MaxLambdaLevelIndex);
	  double TmpPrefactor = (1.0 - (2.0 * ((double) (LambdaLevel & 1)))) * sqrt (((double) this->NbrBosons) * ((double) (TmpMaxAngularMomentumLambdaLevel + 1)) / (4.0 * M_PI * (this->LzMax + 1)));
	  
	  ClebschGordanCoefficients Clebsch(TmpMaxAngularMomentumLambdaLevel, this->NbrBosons - 1);
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      ThreeOrbitalOverlaps[LambdaLevel][i] = new double [this->NbrBosons];
	      for (int j = 0; j < this->NbrBosons; ++j)
		{
		  if (((j - i + LambdaLevel) >= 0) && ((j - i + LambdaLevel) <= TmpMaxAngularMomentumLambdaLevel))
		    {
		      if (unnormalizedFlag == false)
			{
//			  cout << i << " " << j << " :" << (-2 * (j - i + LambdaLevel) + TmpMaxAngularMomentumLambdaLevel) << " " << ((2 * j) - (this->NbrBosons - 1)) <<  " " << this->LzMax << endl;		      
			  ThreeOrbitalOverlaps[LambdaLevel][i][j] = (TmpPrefactor * Clebsch.GetCoefficient(- 2 * (j - i + LambdaLevel) + TmpMaxAngularMomentumLambdaLevel, (2 * j) - (this->NbrBosons - 1), this->LzMax) 
								     * Clebsch.GetCoefficient(-(fermionicSpace->LzMax - 2 * MaxLambdaLevelIndex), this->NbrBosons - 1, this->LzMax));
// 			  cout << LambdaLevel << " " << i << " " << j << " : " << ThreeOrbitalOverlaps[LambdaLevel][i][j] << " " << TmpPrefactor 
// 			       << " " << Clebsch.GetCoefficient(-(fermionicSpace->LzMax - 2 * MaxLambdaLevelIndex), this->NbrBosons - 1, this->LzMax) << endl;

			}
		      else
			{
			  ThreeOrbitalOverlaps[LambdaLevel][i][j] = 1.0;
			}
		    }
		  else
		    {
		      ThreeOrbitalOverlaps[LambdaLevel][i][j] = 0.0;
		    }		  
		}
	    }
	}
       for (long i = minIndex; i < MaxIndex; i++)
 	{
 	  if (fermionicState[i] != 0.0)
 	    {
	      int TmpNbrBosonsLLLUp;
	      int TmpNbrBosons2LLUp;
	      int TmpNbrBosonsLLLDown;
	      int TmpNbrBosons2LLDown;
	      fermionicSpace->ConvertToMonomial(fermionicSpace->StateDescription[i], SlaterLLLUp, Slater2LLUp, SlaterLLLDown, Slater2LLDown, TmpNbrBosonsLLLUp, TmpNbrBosonsLLLDown);	      

 	      FQHESphereWithSU2SpinVanDerMondeTimesSlaterOperation TmpOperation(this, true,  SlaterLLLUp, Slater2LLUp, SlaterLLLDown, Slater2LLDown, TmpNbrBosonsLLLUp, TmpNbrBosonsLLLDown, ThreeOrbitalOverlaps);
 	      TmpOperation.ApplyOperation(architecture);
 	      FinalState = TmpOperation.GetDestinationVector();
	      
// 	      this->ReverseVanDerMondeTimesSlater(SlaterLLLUp, Slater2LLUp, SlaterLLLDown, Slater2LLDown, TmpNbrBosonsLLLUp, TmpNbrBosonsLLLDown, 
// 						  FinalState, ThreeOrbitalOverlaps);

//  	      for (int j = 0; j < this->NbrBosons; ++j)
//  		this->ReverseVanDerMondeTimesSlater(SlaterLLLUp, Slater2LLUp, SlaterLLLDown, Slater2LLDown, TmpNbrBosonsLLLUp, TmpNbrBosonsLLLDown, FinalState, ThreeOrbitalOverlaps, j);

 	      for (int Index = 0; Index < FinalState.GetVectorDimension(); ++Index)
 		{
 		  if (FinalState[Index] != 0.0)
 		    {
 		      this->FermionToBoson(this->StateDescriptionUp[Index], this->StateDescriptionDown[Index],
 					   this->TemporaryStateUp, this->TemporaryStateDown); 
 		      Coefficient.SetToOne();
 		      for(int p = 0; p <= this->LzMax; ++p)
 			{
 			  if (this->TemporaryStateUp[p] > 1)
 			    Coefficient.FactorialMultiply(this->TemporaryStateUp[p]);
 			  if (this->TemporaryStateDown[p] > 1)
 			    Coefficient.FactorialMultiply(this->TemporaryStateDown[p]);
 			}		  
 		      if (unnormalizedFlag == false)
 			{
 			  outputVector[Index] += sqrt(Coefficient.GetNumericalValue()) * fermionicState[i] * FinalState[Index];
 			}
		      else
			{
			  outputVector[Index] += Coefficient.GetNumericalValue() * fermionicState[i] * FinalState[Index];
			}
		    }
		}
	    }
	}
    }
  for (int LambdaLevel = 0; LambdaLevel <= MaxLambdaLevelIndex; ++LambdaLevel)
    {      
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  delete[] ThreeOrbitalOverlaps[LambdaLevel][i];
	}  
      delete[] ThreeOrbitalOverlaps[LambdaLevel];
    }
  delete[] ThreeOrbitalOverlaps;
}

// compute the integral \int phi_{k1,0}^* phi_{k2,0} phi_{k3,0} on the cylinder geometry where phi_{k,n} in wavefunction with momentum k and int the n-th LL
// 
// k1 = momentum of the first orbital for the phi_{k1,n} wavefunction
// l1 = square value of magnetic length for the phi_{k1,n} wavefunction
// k2 = momentum of the first orbital for the phi_{k2,n} wavefunction
// l2 = square value of magnetic length for the phi_{k2,n} wavefunction
// k3 = momentum of the first orbital for the phi_{k3,n} wavefunction
// l3 = square value of magnetic length for the phi_{k3,n} wavefunction
// perimeter = cylinder perimeter
// return value = corresponding integral

double BosonOnSphereWithSU2Spin::ComputeIntegralPhi0Phi0Phi0OnCylinder(double k1, double l1, double k2, double l2, double k3, double l3, double perimeter)
{
  return  exp (-2.0 * M_PI * M_PI / (perimeter * perimeter) * ((k1 * k1 * l1) + (k2 * k2 * l2) + (k3 * k3 * l3)
							       - 0.5 * ((k1 + k2 + k3) * (k1 + k2 + k3) * l1))) * sqrt (l1 / perimeter);
}

// Compute the product of a spinful Slater determinant with a Van der Monde determinant
//
// slaterUp = monomial representation of the Slater spin up part
// slaterDown = monomial representation of the Slater spin up part
// finalState = reference on the vector the produced state will be stored
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored

void BosonOnSphereWithSU2Spin::VanDerMondeTimesSlater (unsigned long* slaterUp, unsigned long* slaterDown, RealVector& finalState, 
						       double** threeOrbitalOverlaps)
{
  unsigned long TmpNbrStates = 0;
  long CoefUp = 1;
  long CoefDown = 1;
  unsigned long StateUp [this->NbrBosonsUp];
  unsigned long StateDown [this->NbrBosonsDown];
  unsigned long TmpFinalStateUp;
  unsigned long TmpFinalStateDown;
  unsigned long TmpState = 0ul;
  double Sign = 1.0;
  unsigned long Mask = 0ul;
  unsigned long VanDerMonde [this->NbrBosons];
  unsigned long TmpHeapArray [this->NbrBosons];
  unsigned long TmpDim = ((unsigned long) this->NbrBosons);
  for (unsigned long i = 0ul; i < TmpDim; ++i)
    {
      VanDerMonde[i] = i;
      TmpHeapArray[i] = 0ul;
    }
  finalState.ClearVector();


  double TmpFactorUp = 0.0;
  bool ChangedUp = false;
  for (int i = 0; i < this->NbrBosonsUp; ++i)
    {
      StateUp[i] = slaterUp[i] + VanDerMonde[i];
      TmpFactorUp += threeOrbitalOverlaps[StateUp[i]][VanDerMonde[i]];
    }
  double TmpFactorDown = 0.0;
  bool ChangedDown = false;
  for (int i = 0; i < this->NbrBosonsDown; ++i)
    {
      StateDown[i] = slaterDown[i] + VanDerMonde[this->NbrBosonsUp + i];
      TmpFactorDown += threeOrbitalOverlaps[StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
    }
  
  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
  if (TmpPos != this->HilbertSpaceDimension)
    {
      finalState[TmpPos] += Sign * exp(TmpFactorUp + TmpFactorDown);
    }
  ++TmpNbrStates;
  unsigned long Tmp = 0;
  while (Tmp < TmpDim)
    {
      if (TmpHeapArray[Tmp] < Tmp)
	{
	  if ((Tmp & 0x1ul) == 0x0ul)
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[0];
	      VanDerMonde[0] = Tmp2;
	      if (Tmp >= this->NbrBosonsUp)
		ChangedDown = true;
	      ChangedUp = true;
	      
	    }
	  else
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[TmpHeapArray[Tmp]];
	      VanDerMonde[TmpHeapArray[Tmp]] = Tmp2;
	      if (TmpHeapArray[Tmp] >= this->NbrBosonsUp)
		{
		  ChangedDown = true;
		}
	      else
		{
		  ChangedUp = true;		  
		  if (Tmp >= this->NbrBosonsUp)
		    {
		      ChangedDown = true;
		    }
		}
	    }
	  if (ChangedUp == true)
	    {
	      TmpFactorUp = 0.0;
	      for (int i = 0; i < this->NbrBosonsUp; ++i)
		{
		  StateUp[i] = slaterUp[i] + VanDerMonde[i];
		  TmpFactorUp += threeOrbitalOverlaps[StateUp[i]][VanDerMonde[i]];
		}
	      this->ConvertFromMonomial(StateUp, this->TemporaryStateUp, this->NbrBosonsUp);
	      TmpFinalStateUp = this->BosonToFermion(this->TemporaryStateUp);
	      ChangedUp = false;
	    }
	  if (ChangedDown == true)
	    {
	      TmpFactorDown = 0.0;
	      for (int i = 0; i < this->NbrBosonsDown; ++i)
		{
		  StateDown[i] = slaterDown[i] + VanDerMonde[this->NbrBosonsUp + i];
		  TmpFactorDown += threeOrbitalOverlaps[StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
		}
	      this->ConvertFromMonomial(StateDown, this->TemporaryStateDown, this->NbrBosonsDown);
	      TmpFinalStateDown = this->BosonToFermion(this->TemporaryStateDown);
	      ChangedDown = false;	      
	    }
	  Sign *= -1.0;
	  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      finalState[TmpPos] += Sign* exp(TmpFactorUp + TmpFactorDown);
	    }
	  ++TmpHeapArray[Tmp];
	  Tmp = 0ul;
	}
      else
	{
	  TmpHeapArray[Tmp]= 0ul;
	  ++Tmp;
	}
    }
}

// Compute the product of a spinful Slater determinant with a Van der Monde determinant, assuming a reverse flux attachment
//
// slaterUp = monomial representation of the Slater spin up part
// slaterDown = monomial representation of the Slater spin up part
// finalState = reference on the vector the produced state will be stored
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored

void BosonOnSphereWithSU2Spin::ReverseVanDerMondeTimesSlater (unsigned long* slaterUp, unsigned long* slaterDown, RealVector& finalState, 
							      double** threeOrbitalOverlaps)
{
  unsigned long TmpNbrStates = 0;
  long CoefUp = 1;
  long CoefDown = 1;
  unsigned long StateUp [this->NbrBosonsUp];
  unsigned long StateDown [this->NbrBosonsDown];
  unsigned long TmpFinalStateUp;
  unsigned long TmpFinalStateDown;
  unsigned long TmpState = 0ul;
  double Sign = 1.0;
  unsigned long Mask = 0ul;
  unsigned long VanDerMonde [this->NbrBosons];
  int TmpHeapArray [this->NbrBosons];
  int TmpDim = this->NbrBosons;
  for (int i = 0; i < TmpDim; ++i)
    {
      VanDerMonde[i] = (unsigned long) i;
      TmpHeapArray[i] = 0;
    }
  finalState.ClearVector();


  double TmpFactor = 0.0;
  int Tmp = 0;
  bool DiscardFlag = false;
  for (int i = 0; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
    {
      StateUp[i] = VanDerMonde[i] - slaterUp[i];
      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
	{
	  TmpFactor += threeOrbitalOverlaps[StateUp[i]][VanDerMonde[i]];
	}
      else
	{
	  DiscardFlag = true;
	}
    }
  if (DiscardFlag == false)
    {
      for (int i = 0; (i < this->NbrBosonsDown) && (DiscardFlag == false); ++i)
	{
	  StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slaterDown[i];
	  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
	    {
	      TmpFactor += threeOrbitalOverlaps[StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
	    }
	  else
	    {
	      DiscardFlag = true;
	    }
	}
      if (DiscardFlag == false)
	{
	  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
	  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
	  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      finalState[TmpPos] += Sign * exp(TmpFactor);
	    }
	}
    }

  while (Tmp < TmpDim)
    {
      if (TmpHeapArray[Tmp] < Tmp)
	{
	  if ((Tmp & 0x1ul) == 0x0ul)
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[0];
	      VanDerMonde[0] = Tmp2;
	      
	    }
	  else
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[TmpHeapArray[Tmp]];
	      VanDerMonde[TmpHeapArray[Tmp]] = Tmp2;
	    }
	  Sign *= -1.0;
	  DiscardFlag = false;
	  TmpFactor = 0.0;
	  for (int i = 0; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
	    {
	      StateUp[i] = VanDerMonde[i] - slaterUp[i];
	      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
		{
		  TmpFactor += threeOrbitalOverlaps[StateUp[i]][VanDerMonde[i]];
		}
	      else
		{
		  DiscardFlag = true;
		}
	    }
	  if (DiscardFlag == false)
	    {
	      for (int i = 0; (i < this->NbrBosonsDown) && (DiscardFlag == false); ++i)
		{
		  StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slaterDown[i];
		  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
		    {
		      TmpFactor += threeOrbitalOverlaps[StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
		    }
		  else
		    {
		      DiscardFlag = true;
		    }
		}
	      if (DiscardFlag == false)
		{
		  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
		  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
		  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
		  if (TmpPos != this->HilbertSpaceDimension)
		    {
		      finalState[TmpPos] += Sign* exp(TmpFactor);
		    }
		}
	    }
	  ++TmpHeapArray[Tmp];
	  Tmp = 0;
	}
      else
	{
	  TmpHeapArray[Tmp]= 0;
	  ++Tmp;
	}
    }
}

// Compute the product of a spinful Slater determinant with a Van der Monde determinant, assuming a reverse flux attachment and peforming only (N-1)! permutations
//
// slaterUp = monomial representation of the Slater spin up part
// slaterDown = monomial representation of the Slater spin up part
// finalState = reference on the vector the produced state will be stored
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored
// position = perform a swap between the last element the one at given position 

void BosonOnSphereWithSU2Spin::ReverseVanDerMondeTimesSlater (unsigned long* slaterUp, unsigned long* slaterDown, RealVector& finalState, 
							      double** threeOrbitalOverlaps, int position)
{
  unsigned long TmpNbrStates = 0;
  long CoefUp = 1;
  long CoefDown = 1;
  unsigned long StateUp [this->NbrBosonsUp];
  unsigned long StateDown [this->NbrBosonsDown];
  unsigned long TmpFinalStateUp;
  unsigned long TmpFinalStateDown;
  unsigned long TmpState = 0ul;
  double Sign = 1.0;
  unsigned long Mask = 0ul;
  unsigned long VanDerMonde [this->NbrBosons];
  int TmpHeapArray [this->NbrBosons];
  int TmpDim = this->NbrBosons - 1;
  for (int i = 0; i < TmpDim; ++i)
    {
      VanDerMonde[i] = (unsigned long) i;
      TmpHeapArray[i] = 0;
    }
  VanDerMonde[TmpDim] = position;
  VanDerMonde[position] = TmpDim;
  if (position != TmpDim)
    {
      Sign *= -1.0;
    }
  double TmpFactor = 0.0;
  int Tmp = 0;
  bool DiscardFlag = false;
  for (int i = 0; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
    {
      StateUp[i] = VanDerMonde[i] - slaterUp[i];
      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
	{
	  TmpFactor += threeOrbitalOverlaps[StateUp[i]][VanDerMonde[i]];
	}
      else
	{
	  DiscardFlag = true;
	}
    }
  if (DiscardFlag == false)
    {
      for (int i = 0; (i < this->NbrBosonsDown) && (DiscardFlag == false); ++i)
	{
	  StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slaterDown[i];
	  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
	    {
	      TmpFactor += threeOrbitalOverlaps[StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
	    }
	  else
	    {
	      DiscardFlag = true;
	    }
	}
      if (DiscardFlag == false)
	{
	  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
	  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
	  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      finalState[TmpPos] += Sign * exp(TmpFactor);
	    }
	}
    }

  while (Tmp < TmpDim)
    {
      if (TmpHeapArray[Tmp] < Tmp)
	{
	  if ((Tmp & 0x1) == 0x0)
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[0];
	      VanDerMonde[0] = Tmp2;
	      
	    }
	  else
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[TmpHeapArray[Tmp]];
	      VanDerMonde[TmpHeapArray[Tmp]] = Tmp2;
	    }
	  Sign *= -1.0;
	  DiscardFlag = false;
	  TmpFactor = 0.0;
	  for (int i = 0; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
	    {
	      StateUp[i] = VanDerMonde[i] - slaterUp[i];
	      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
		{
		  TmpFactor += threeOrbitalOverlaps[StateUp[i]][VanDerMonde[i]];
		}
	      else
		{
		  DiscardFlag = true;
		}
	    }
	  if (DiscardFlag == false)
	    {
	      for (int i = 0; (i < this->NbrBosonsDown) && (DiscardFlag == false); ++i)
		{
		  StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slaterDown[i];
		  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
		    {
		      TmpFactor += threeOrbitalOverlaps[StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
		    }
		  else
		    {
		      DiscardFlag = true;
		    }
		}
	      if (DiscardFlag == false)
		{
		  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
		  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
		  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
		  if (TmpPos != this->HilbertSpaceDimension)
		    {
		      finalState[TmpPos] += Sign* exp(TmpFactor);
		    }
		}
	    }
	  ++TmpHeapArray[Tmp];
	  Tmp = 0;
	}
      else
	{
	  TmpHeapArray[Tmp]= 0;
	  ++Tmp;
	}
    }
}

// Compute the product of a spinful Slater determinant in two Landau levels with a Van der Monde determinant
//
// slaterLLLUp = monomial representation of the lowest Landau part of the Slater spin up part
// slater2LLUp = monomial representation of the second Landau part of the Slater spin up part
// slaterLLLDown = monomial representation of the lowest Landau part  of the Slater spin down part
// slater2LLDown = monomial representation of the second Landau part of the Slater spin down part
// nbrBosonsLLLUp - number of spin up bosons in the lowest Landau level
// nbrBosonsLLLDown - number of spin down bosons in the lowest Landau level
// finalState = reference on the vector the produced state will be stored
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored

void BosonOnSphereWithSU2Spin::VanDerMondeTimesSlater (unsigned long* slaterLLLUp, unsigned long* slater2LLUp, 
						       unsigned long* slaterLLLDown, unsigned long* slater2LLDown, 
						       int nbrBosonsLLLUp, int nbrBosonsLLLDown, 
						       RealVector& finalState, double*** threeOrbitalOverlaps)
{
  unsigned long TmpNbrStates = 0;
  long CoefUp = 1;
  long CoefDown = 1;
  unsigned long StateUp [this->NbrBosonsUp];
  unsigned long StateDown [this->NbrBosonsDown];
  unsigned long TmpFinalStateUp;
  unsigned long TmpFinalStateDown;
  unsigned long TmpState = 0ul;
  double Sign = 1.0;
  unsigned long Mask = 0ul;
  unsigned long VanDerMonde [this->NbrBosons];
  int TmpHeapArray [this->NbrBosons];
  int TmpDim = this->NbrBosons;
  int Tmp = 0;
  for (int i = 0; i < TmpDim; ++i)
    {
      VanDerMonde[i] = (unsigned long) i;
      TmpHeapArray[i] = 0;
    }
  finalState.ClearVector();

  double TmpFactor = 1.0;
  bool DiscardFlag = false;
  for (int i = 0; (i < nbrBosonsLLLUp) && (DiscardFlag == false); ++i)
    {
      StateUp[i] = slaterLLLUp[i] + VanDerMonde[i] - 1;
      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
	{
	  TmpFactor *= threeOrbitalOverlaps[0][StateUp[i]][VanDerMonde[i]];
	}
      else
	{
	  DiscardFlag = true;
	}
    }
  if (DiscardFlag == false)
    {
      for (int i = nbrBosonsLLLUp; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
	{
	  StateUp[i] = slater2LLUp[i - nbrBosonsLLLUp] + VanDerMonde[i] - 1;
	  if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
	    {
	      TmpFactor *= threeOrbitalOverlaps[1][StateUp[i]][VanDerMonde[i]];
	    }
	  else
	    {
	      DiscardFlag = true;
	    }
	}
      if (DiscardFlag == false)
	{
	  for (int i = 0; i < (nbrBosonsLLLDown) && (DiscardFlag == false); ++i)
	    {
	      StateDown[i] = slaterLLLDown[i] + VanDerMonde[this->NbrBosonsUp + i] - 1;
	      if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
		{
		  TmpFactor *= threeOrbitalOverlaps[0][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
		}
	      else
		{
		  DiscardFlag = true;
		}
	    }
	  if (DiscardFlag == false)
	    {
	      for (int i = nbrBosonsLLLDown; i < (this->NbrBosonsDown) && (DiscardFlag == false); ++i)
		{
		  StateDown[i] = slater2LLDown[i - nbrBosonsLLLDown] + VanDerMonde[this->NbrBosonsUp + i] - 1;
		  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
		    {
		      TmpFactor *= threeOrbitalOverlaps[1][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
		    }
		  else
		    {
		      DiscardFlag = true;
		    }
		}
	      if (DiscardFlag == false)
		{
		  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
		  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
		  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
		  if (TmpPos != this->HilbertSpaceDimension)
		    {
		      finalState[TmpPos] += Sign * TmpFactor;
		    }
		}
	    }
	}
    }
  
  while (Tmp < TmpDim)
    {
      if (TmpHeapArray[Tmp] < Tmp)
	{
	  if ((Tmp & 0x1) == 0x0)
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[0];
	      VanDerMonde[0] = Tmp2;
	      
	    }
	  else
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[TmpHeapArray[Tmp]];
	      VanDerMonde[TmpHeapArray[Tmp]] = Tmp2;
	    }
	  Sign *= -1.0;
	  DiscardFlag = false;
	  TmpFactor = 1.0;
	  for (int i = 0; (i < nbrBosonsLLLUp) && (DiscardFlag == false); ++i)
	    {
	      StateUp[i] = slaterLLLUp[i] + VanDerMonde[i] - 1;
	      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
		{
		  TmpFactor *= threeOrbitalOverlaps[0][StateUp[i]][VanDerMonde[i]];
		}
	      else
		{
		  DiscardFlag = true;
		}
	    }
	  if (DiscardFlag == false)
	    {
	      for (int i = nbrBosonsLLLUp; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
		{
		  StateUp[i] = slater2LLUp[i - nbrBosonsLLLUp] + VanDerMonde[i] - 1;
		  if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
		    {
		      TmpFactor *= threeOrbitalOverlaps[1][StateUp[i]][VanDerMonde[i]];
		    }
		  else
		    {
		      DiscardFlag = true;
		    }
		}
	      if (DiscardFlag == false)
		{
		  for (int i = 0; (i < nbrBosonsLLLDown) && (DiscardFlag == false); ++i)
		    {
		      StateDown[i] = slaterLLLDown[i] + VanDerMonde[this->NbrBosonsUp + i] - 1;
		      if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
			{
			  TmpFactor *= threeOrbitalOverlaps[0][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
			}
		      else
			{
			  DiscardFlag = true;
			}
		    }
		  if (DiscardFlag == false)
		    {
		      for (int i = nbrBosonsLLLDown; (i < this->NbrBosonsDown) && (DiscardFlag == false); ++i)
			{
			  StateDown[i] = slater2LLDown[i - nbrBosonsLLLDown] + VanDerMonde[this->NbrBosonsUp + i] - 1;
			  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
			    {
			      TmpFactor *= threeOrbitalOverlaps[1][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
			    }
			  else
			    {
			      DiscardFlag = true;
			    }
			}
		      if (DiscardFlag == false)
			{
			  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
			  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
			  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
			  if (TmpPos != this->HilbertSpaceDimension)
			    {
			      finalState[TmpPos] += Sign* TmpFactor;
			    }
			}
		    }
		}
	    }
	  ++TmpHeapArray[Tmp];
	  Tmp = 0;
	}
      else
	{
	  TmpHeapArray[Tmp]= 0;
	  ++Tmp;
	}
    }
}

// Compute the product of a spinful Slater determinant in two Landau levels with a Van der Monde determinant, assuming a reverse flux attachment
//
// slaterLLLUp = monomial representation of the lowest Landau part of the Slater spin up part
// slater2LLUp = monomial representation of the second Landau part of the Slater spin up part
// slaterLLLDown = monomial representation of the lowest Landau part  of the Slater spin down part
// slater2LLDown = monomial representation of the second Landau part of the Slater spin down part
// nbrBosonsLLLUp - number of spin up bosons in the lowest Landau level
// nbrBosonsLLLDown - number of spin down bosons in the lowest Landau level
// finalState = reference on the vector the produced state will be stored
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored

void BosonOnSphereWithSU2Spin::ReverseVanDerMondeTimesSlater (unsigned long* slaterLLLUp, unsigned long* slater2LLUp, unsigned long* slaterLLLDown, unsigned long* slater2LLDown, 
							      int nbrBosonsLLLUp, int nbrBosonsLLLDown, RealVector& finalState, double*** threeOrbitalOverlaps)
{
  unsigned long TmpNbrStates = 0;
  long CoefUp = 1;
  long CoefDown = 1;
  unsigned long StateUp [this->NbrBosonsUp];
  unsigned long StateDown [this->NbrBosonsDown];
  unsigned long TmpFinalStateUp;
  unsigned long TmpFinalStateDown;
  unsigned long TmpState = 0ul;
  double Sign = 1.0;
  unsigned long Mask = 0ul;
  unsigned long VanDerMonde [this->NbrBosons];
  int TmpHeapArray [this->NbrBosons];
  int TmpDim = this->NbrBosons;
  int Tmp = 0;
  for (int i = 0; i < TmpDim; ++i)
    {
      VanDerMonde[i] = (unsigned long) i;
      TmpHeapArray[i] = 0;
    }
  finalState.ClearVector();

  double TmpFactor = 1.0;
  bool DiscardFlag = false;
  for (int i = 0; (i < nbrBosonsLLLUp) && (DiscardFlag == false); ++i)
    {
      StateUp[i] = VanDerMonde[i] -  slaterLLLUp[i] + 1;
      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
	{
	  TmpFactor *= threeOrbitalOverlaps[0][StateUp[i]][VanDerMonde[i]];
	}
      else
	{
	  DiscardFlag = true;
	}
    }
  if (DiscardFlag == false)
    {
      for (int i = nbrBosonsLLLUp; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
	{
	  StateUp[i] = VanDerMonde[i] - slater2LLUp[i - nbrBosonsLLLUp] + 1;
	  if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
	    {
	      TmpFactor *= threeOrbitalOverlaps[1][StateUp[i]][VanDerMonde[i]];
	    }
	  else
	    {
	      DiscardFlag = true;
	    }
	}
      if (DiscardFlag == false)
	{
	  for (int i = 0; i < (nbrBosonsLLLDown) && (DiscardFlag == false); ++i)
	    {
	      StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slaterLLLDown[i] + 1;
	      if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
		{
		  TmpFactor *= threeOrbitalOverlaps[0][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
		}
	      else
		{
		  DiscardFlag = true;
		}
	    }
	  if (DiscardFlag == false)
	    {
	      for (int i = nbrBosonsLLLDown; i < (this->NbrBosonsDown) && (DiscardFlag == false); ++i)
		{
		  StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slater2LLDown[i - nbrBosonsLLLDown] + 1;
		  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
		    {
		      TmpFactor *= threeOrbitalOverlaps[1][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
		    }
		  else
		    {
		      DiscardFlag = true;
		    }
		}
	      if (DiscardFlag == false)
		{
		  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
		  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
		  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
		  if (TmpPos != this->HilbertSpaceDimension)
		    {
		      finalState[TmpPos] += Sign * TmpFactor;
		    }
		}
	    }
	}
    }
  
  while (Tmp < TmpDim)
    {
      if (TmpHeapArray[Tmp] < Tmp)
	{
	  if ((Tmp & 0x1) == 0x0)
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[0];
	      VanDerMonde[0] = Tmp2;
	      
	    }
	  else
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[TmpHeapArray[Tmp]];
	      VanDerMonde[TmpHeapArray[Tmp]] = Tmp2;
	    }
	  Sign *= -1.0;
	  DiscardFlag = false;
	  TmpFactor = 1.0;
	  for (int i = 0; (i < nbrBosonsLLLUp) && (DiscardFlag == false); ++i)
	    {
	      StateUp[i] = VanDerMonde[i] - slaterLLLUp[i] + 1;
	      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
		{
		  TmpFactor *= threeOrbitalOverlaps[0][StateUp[i]][VanDerMonde[i]];
		}
	      else
		{
		  DiscardFlag = true;
		}
	    }
	  if (DiscardFlag == false)
	    {
	      for (int i = nbrBosonsLLLUp; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
		{
		  StateUp[i] = VanDerMonde[i] - slater2LLUp[i - nbrBosonsLLLUp] + 1;
		  if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
		    {
		      TmpFactor *= threeOrbitalOverlaps[1][StateUp[i]][VanDerMonde[i]];
		    }
		  else
		    {
		      DiscardFlag = true;
		    }
		}
	      if (DiscardFlag == false)
		{
		  for (int i = 0; (i < nbrBosonsLLLDown) && (DiscardFlag == false); ++i)
		    {
		      StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slaterLLLDown[i] + 1;
		      if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
			{
			  TmpFactor *= threeOrbitalOverlaps[0][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
			}
		      else
			{
			  DiscardFlag = true;
			}
		    }
		  if (DiscardFlag == false)
		    {
		      for (int i = nbrBosonsLLLDown; (i < this->NbrBosonsDown) && (DiscardFlag == false); ++i)
			{
			  StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slater2LLDown[i - nbrBosonsLLLDown] + 1;
			  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
			    {
			      TmpFactor *= threeOrbitalOverlaps[1][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
			    }
			  else
			    {
			      DiscardFlag = true;
			    }
			}
		      if (DiscardFlag == false)
			{
			  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
			  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
			  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
			  if (TmpPos != this->HilbertSpaceDimension)
			    {
			      finalState[TmpPos] += Sign* TmpFactor;
			    }
			}
		    }
		}
	    }
	  ++TmpHeapArray[Tmp];
	  Tmp = 0;
	}
      else
	{
	  TmpHeapArray[Tmp]= 0;
	  ++Tmp;
	}
    }
}
  
// Compute the product of a spinful Slater determinant in two Landau levels with a Van der Monde determinant, assuming a reverse flux attachment
//
// slaterLLLUp = monomial representation of the lowest Landau part of the Slater spin up part
// slater2LLUp = monomial representation of the second Landau part of the Slater spin up part
// slaterLLLDown = monomial representation of the lowest Landau part  of the Slater spin down part
// slater2LLDown = monomial representation of the second Landau part of the Slater spin down part
// nbrBosonsLLLUp - number of spin up bosons in the lowest Landau level
// nbrBosonsLLLDown - number of spin down bosons in the lowest Landau level
// finalState = reference on the vector the produced state will be stored
// threeOrbitalOverlaps = array where the integrals of the three orbital product are stored
// position = perform a swap between the last element the one at given position 

void BosonOnSphereWithSU2Spin::ReverseVanDerMondeTimesSlater (unsigned long* slaterLLLUp, unsigned long* slater2LLUp, unsigned long* slaterLLLDown, unsigned long* slater2LLDown, 
							      int nbrBosonsLLLUp, int nbrBosonsLLLDown, RealVector& finalState, double*** threeOrbitalOverlaps, int position)
{
  unsigned long TmpNbrStates = 0;
  long CoefUp = 1;
  long CoefDown = 1;
  unsigned long StateUp [this->NbrBosonsUp];
  unsigned long StateDown [this->NbrBosonsDown];
  unsigned long TmpFinalStateUp;
  unsigned long TmpFinalStateDown;
  unsigned long TmpState = 0ul;
  double Sign = 1.0;
  unsigned long Mask = 0ul;
  unsigned long VanDerMonde [this->NbrBosons];
  int TmpHeapArray [this->NbrBosons];
  int TmpDim = this->NbrBosons - 1;
  int Tmp = 0;
  for (int i = 0; i < TmpDim; ++i)
    {
      VanDerMonde[i] = (unsigned long) i;
      TmpHeapArray[i] = 0;
    }
  VanDerMonde[TmpDim] = position;
  VanDerMonde[position] = TmpDim;
  if (position != TmpDim)
    {
      Sign *= -1.0;
    }

  double TmpFactor = 1.0;
  bool DiscardFlag = false;
  for (int i = 0; (i < nbrBosonsLLLUp) && (DiscardFlag == false); ++i)
    {
      StateUp[i] = VanDerMonde[i] -  slaterLLLUp[i] + 1;
      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
	{
	  TmpFactor *= threeOrbitalOverlaps[0][StateUp[i]][VanDerMonde[i]];
	}
      else
	{
	  DiscardFlag = true;
	}
    }
  if (DiscardFlag == false)
    {
      for (int i = nbrBosonsLLLUp; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
	{
	  StateUp[i] = VanDerMonde[i] - slater2LLUp[i - nbrBosonsLLLUp] + 1;
	  if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
	    {
	      TmpFactor *= threeOrbitalOverlaps[1][StateUp[i]][VanDerMonde[i]];
	    }
	  else
	    {
	      DiscardFlag = true;
	    }
	}
      if (DiscardFlag == false)
	{
	  for (int i = 0; i < (nbrBosonsLLLDown) && (DiscardFlag == false); ++i)
	    {
	      StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slaterLLLDown[i] + 1;
	      if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
		{
		  TmpFactor *= threeOrbitalOverlaps[0][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
		}
	      else
		{
		  DiscardFlag = true;
		}
	    }
	  if (DiscardFlag == false)
	    {
	      for (int i = nbrBosonsLLLDown; i < (this->NbrBosonsDown) && (DiscardFlag == false); ++i)
		{
		  StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slater2LLDown[i - nbrBosonsLLLDown] + 1;
		  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
		    {
		      TmpFactor *= threeOrbitalOverlaps[1][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
		    }
		  else
		    {
		      DiscardFlag = true;
		    }
		}
	      if (DiscardFlag == false)
		{
		  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
		  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
		  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
		  if (TmpPos != this->HilbertSpaceDimension)
		    {
		      finalState[TmpPos] += Sign * TmpFactor;
		    }
		}
	    }
	}
    }
  
  while (Tmp < TmpDim)
    {
      if (TmpHeapArray[Tmp] < Tmp)
	{
	  if ((Tmp & 0x1) == 0x0)
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[0];
	      VanDerMonde[0] = Tmp2;
	      
	    }
	  else
	    {
	      unsigned long Tmp2 = VanDerMonde[Tmp];
	      VanDerMonde[Tmp] = VanDerMonde[TmpHeapArray[Tmp]];
	      VanDerMonde[TmpHeapArray[Tmp]] = Tmp2;
	    }
	  Sign *= -1.0;
	  DiscardFlag = false;
	  TmpFactor = 1.0;
	  for (int i = 0; (i < nbrBosonsLLLUp) && (DiscardFlag == false); ++i)
	    {
	      StateUp[i] = VanDerMonde[i] - slaterLLLUp[i] + 1;
	      if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
		{
		  TmpFactor *= threeOrbitalOverlaps[0][StateUp[i]][VanDerMonde[i]];
		}
	      else
		{
		  DiscardFlag = true;
		}
	    }
	  if (DiscardFlag == false)
	    {
	      for (int i = nbrBosonsLLLUp; (i < this->NbrBosonsUp) && (DiscardFlag == false); ++i)
		{
		  StateUp[i] = VanDerMonde[i] - slater2LLUp[i - nbrBosonsLLLUp] + 1;
		  if ((StateUp[i] >= 0) && (StateUp[i] <= this->LzMax))
		    {
		      TmpFactor *= threeOrbitalOverlaps[1][StateUp[i]][VanDerMonde[i]];
		    }
		  else
		    {
		      DiscardFlag = true;
		    }
		}
	      if (DiscardFlag == false)
		{
		  for (int i = 0; (i < nbrBosonsLLLDown) && (DiscardFlag == false); ++i)
		    {
		      StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slaterLLLDown[i] + 1;
		      if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
			{
			  TmpFactor *= threeOrbitalOverlaps[0][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
			}
		      else
			{
			  DiscardFlag = true;
			}
		    }
		  if (DiscardFlag == false)
		    {
		      for (int i = nbrBosonsLLLDown; (i < this->NbrBosonsDown) && (DiscardFlag == false); ++i)
			{
			  StateDown[i] = VanDerMonde[this->NbrBosonsUp + i] - slater2LLDown[i - nbrBosonsLLLDown] + 1;
			  if ((StateDown[i] >= 0) && (StateDown[i] <= this->LzMax))
			    {
			      TmpFactor *= threeOrbitalOverlaps[1][StateDown[i]][VanDerMonde[this->NbrBosonsUp + i]];
			    }
			  else
			    {
			      DiscardFlag = true;
			    }
			}
		      if (DiscardFlag == false)
			{
			  this->ConvertFromMonomial(StateUp, StateDown, this->TemporaryStateUp, this->TemporaryStateDown);
			  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
			  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
			  if (TmpPos != this->HilbertSpaceDimension)
			    {
			      finalState[TmpPos] += Sign* TmpFactor;
			    }
			}
		    }
		}
	    }
	  ++TmpHeapArray[Tmp];
	  Tmp = 0;
	}
      else
	{
	  TmpHeapArray[Tmp]= 0;
	  ++Tmp;
	}
    }
}
  
// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnSphereWithSU2Spin::EvaluatePartialDensityMatrix (int subsytemSize, int nbrParticleSector, int lzSector, int szSector, RealVector& groundState)
{
  RealMatrix TmpDensityMatrix;
  return TmpDensityMatrix;
 }

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealMatrix BosonOnSphereWithSU2Spin::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrParticleSector, int lzSector, int szSector, RealVector& groundState)
{
  int nbrOrbitalA = subsytemSize;
  int nbrOrbitalB = this->LzMax + 1 - nbrOrbitalA;  
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;


  int TotalLzDisk = ConvertLzFromSphereToDisk(this->TotalLz, this->NbrBosons, this->LzMax);  
  int LzSectorDisk = ConvertLzFromSphereToDisk(lzSector, nbrParticleSector, nbrOrbitalA - 1);
  int ComplementaryLzSectorDisk = (TotalLzDisk - LzSectorDisk) - (ComplementaryNbrParticles * nbrOrbitalA);
  int ComplementaryLzSector = ConvertLzFromDiskToSphere(ComplementaryLzSectorDisk, ComplementaryNbrParticles, nbrOrbitalB - 1);

  if ((abs(ComplementarySzSector) > ComplementaryNbrParticles) || (abs(ComplementaryLzSector) > ((nbrOrbitalB - 1) * ComplementaryNbrParticles)))
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;  
    }
  BosonOnSphereWithSU2Spin SubsytemSpace(nbrParticleSector, lzSector, nbrOrbitalA - 1, szSector);
  BosonOnSphereWithSU2Spin ComplementarySubsytemSpace(ComplementaryNbrParticles, ComplementaryLzSector, nbrOrbitalB - 1, ComplementarySzSector);  
  RealMatrix TmpEntanglementMatrix(SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySubsytemSpace.GetHilbertSpaceDimension(), true);
  long TmpNbrNonZeroElements = 0l;
  unsigned long** TmpSubsytemSpaceOccupationNumbersUp = new unsigned long* [SubsytemSpace.HilbertSpaceDimension];
  unsigned long** TmpSubsytemSpaceOccupationNumbersDown = new unsigned long* [SubsytemSpace.HilbertSpaceDimension];
  unsigned long TmpStateUp;
  unsigned long TmpStateDown;
  for (int i = 0; i < SubsytemSpace.HilbertSpaceDimension; ++i)
    {
      TmpSubsytemSpaceOccupationNumbersUp[i] = new unsigned long [SubsytemSpace.NbrLzValue];
      TmpSubsytemSpaceOccupationNumbersDown[i] = new unsigned long [SubsytemSpace.NbrLzValue];
      SubsytemSpace.FermionToBoson(SubsytemSpace.StateDescriptionUp[i], SubsytemSpace.StateDescriptionDown[i], 
				   TmpSubsytemSpaceOccupationNumbersUp[i], TmpSubsytemSpaceOccupationNumbersDown[i]);
    }

  for (int MinIndex = 0; MinIndex < ComplementarySubsytemSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      ComplementarySubsytemSpace.FermionToBoson(ComplementarySubsytemSpace.StateDescriptionUp[MinIndex], ComplementarySubsytemSpace.StateDescriptionDown[MinIndex],  
						ComplementarySubsytemSpace.TemporaryStateUp, ComplementarySubsytemSpace.TemporaryStateDown);
      for (int j = 0; j < SubsytemSpace.HilbertSpaceDimension; ++j)
	{
	  for (int i = 0; i <= SubsytemSpace.LzMax; ++i)
	    {
	      this->TemporaryStateUp[i] = TmpSubsytemSpaceOccupationNumbersUp[j][i];
	      this->TemporaryStateDown[i] = TmpSubsytemSpaceOccupationNumbersDown[j][i];
	    }
	  for (int i = 0; i <= ComplementarySubsytemSpace.LzMax; ++i)
	    {
	      this->TemporaryStateUp[i + nbrOrbitalA] = ComplementarySubsytemSpace.TemporaryStateUp[i];
	      this->TemporaryStateDown[i + nbrOrbitalA] = ComplementarySubsytemSpace.TemporaryStateDown[i];
	    }
	  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpStateUp, TmpStateDown);
	  int TmpPos = this->FindStateIndex(TmpStateUp, TmpStateDown);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpEntanglementMatrix.SetMatrixElement(j, MinIndex, groundState[TmpPos]);
	      ++TmpNbrNonZeroElements;
	    }
	}
    }

  for (int i = 0; i < SubsytemSpace.HilbertSpaceDimension; ++i)
    {
      delete[] TmpSubsytemSpaceOccupationNumbersUp[i];
      delete[] TmpSubsytemSpaceOccupationNumbersDown[i];
    }
  delete[] TmpSubsytemSpaceOccupationNumbersUp;
  delete[] TmpSubsytemSpaceOccupationNumbersDown;
  if (TmpNbrNonZeroElements > 0l)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}


// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix BosonOnSphereWithSU2Spin::EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, RealVector& groundState,
											 bool removeBinomialCoefficient, AbstractArchitecture* architecture)
{
  int nbrOrbitalA = this->LzMax + 1;
  int nbrOrbitalB = this->LzMax + 1;  
  if (nbrParticleSector == 0)
    {
      if ((lzSector == 0) && (szSector == 0))
	{
	  RealMatrix TmpEntanglementMatrix(1, this->HilbertSpaceDimension, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    TmpEntanglementMatrix.SetMatrixElement(0, i, groundState[i]);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  if (nbrParticleSector == this->NbrBosons)
    {
      if ((lzSector == this->TotalLz) && (szSector == this->TotalSpin))
	{
	  RealMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[i]);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    } 
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  int ComplementaryLzSector = this->TotalLz - lzSector;
  if ((abs(ComplementarySzSector) > ComplementaryNbrParticles) || (abs(ComplementaryLzSector) > (this->LzMax * ComplementaryNbrParticles)))
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;  
    }

  BosonOnSphereWithSU2Spin SubsytemSpace(nbrParticleSector, lzSector, this->LzMax, szSector);
  BosonOnSphereWithSU2Spin ComplementarySubsytemSpace(ComplementaryNbrParticles, ComplementaryLzSector, this->LzMax, ComplementarySzSector);  
  RealMatrix TmpEntanglementMatrix(SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySubsytemSpace.GetHilbertSpaceDimension(), true);

  long TmpNbrNonZeroElements = this->EvaluatePartialEntanglementMatrixParticlePartitionCore(0, ComplementarySubsytemSpace.GetHilbertSpaceDimension(),
											    &ComplementarySubsytemSpace, &SubsytemSpace, 
											    groundState, &TmpEntanglementMatrix, removeBinomialCoefficient);

  if (TmpNbrNonZeroElements > 0l)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}
   
// core part of the entanglement matrix evaluation for the particle partition
// 
// minIndex = first index to consider in the complementary Hilbert space
// nbrIndex = number of indices to consider in the complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// entanglementMatrix = pointer to entanglement matrix
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = number of components that have been added to the entanglement matrix

long BosonOnSphereWithSU2Spin::EvaluatePartialEntanglementMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace, ParticleOnSphere* destinationHilbertSpace, 
										       RealVector& groundState, RealMatrix* entanglementMatrix, bool removeBinomialCoefficient)
{

  BosonOnSphereWithSU2Spin* SubsytemSpace = (BosonOnSphereWithSU2Spin*) ((BosonOnSphereWithSU2Spin*) destinationHilbertSpace)->Clone();
  BosonOnSphereWithSU2Spin* ComplementarySubsytemSpace = (BosonOnSphereWithSU2Spin*) ((BosonOnSphereWithSU2Spin*) complementaryHilbertSpace)->Clone();

  int NbrOrbitalA = SubsytemSpace->LzMax + 1;
  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
  double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementarySubsytemSpace->NbrBosons] - LogFactorials[SubsytemSpace->NbrBosons];
  if (removeBinomialCoefficient == true)
    TmpLogBinomial = 0.0;

  long TmpNbrNonZeroElements = 0l;
  unsigned long** TmpSubsytemSpaceOccupationNumbersUp = new unsigned long* [SubsytemSpace->HilbertSpaceDimension];
  unsigned long** TmpSubsytemSpaceOccupationNumbersDown = new unsigned long* [SubsytemSpace->HilbertSpaceDimension];
  unsigned long** TmpSubsytemSpaceMonomialUp = new unsigned long* [SubsytemSpace->HilbertSpaceDimension];
  unsigned long** TmpSubsytemSpaceMonomialDown = new unsigned long* [SubsytemSpace->HilbertSpaceDimension];
  double* TmpSubsytemLogFactorials = new double [SubsytemSpace->HilbertSpaceDimension];
  unsigned long* TmpMonomialUp1 = new unsigned long [ComplementarySubsytemSpace->NbrBosons];
  unsigned long* TmpMonomialDown1 = new unsigned long [ComplementarySubsytemSpace->NbrBosons];
  unsigned long* TmpMonomialUp3 = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomialDown3 = new unsigned long [this->NbrBosons];

  for (int i = 0; i < SubsytemSpace->HilbertSpaceDimension; ++i)
    {
      TmpSubsytemSpaceOccupationNumbersUp[i] = new unsigned long [this->NbrLzValue];
      TmpSubsytemSpaceOccupationNumbersDown[i] = new unsigned long [this->NbrLzValue];
      TmpSubsytemSpaceMonomialUp[i] = new unsigned long [SubsytemSpace->NbrBosons];
      TmpSubsytemSpaceMonomialDown[i] = new unsigned long [SubsytemSpace->NbrBosons];
      SubsytemSpace->FermionToBoson(SubsytemSpace->StateDescriptionUp[i], SubsytemSpace->StateDescriptionDown[i], 
				      TmpSubsytemSpaceOccupationNumbersUp[i], TmpSubsytemSpaceOccupationNumbersDown[i]);
      SubsytemSpace->ConvertToMonomial(SubsytemSpace->StateDescriptionUp[i], SubsytemSpace->StateDescriptionDown[i], 
					 TmpSubsytemSpaceMonomialUp[i], TmpSubsytemSpaceMonomialDown[i]);
      unsigned long* TmpOccupationNumberUp = TmpSubsytemSpaceOccupationNumbersUp[i];
      unsigned long* TmpOccupationNumberDown = TmpSubsytemSpaceOccupationNumbersDown[i];
      double TmpFactor = 0.0;
      for (int k = 0; k <= SubsytemSpace->LzMax; ++k)
	{
	  TmpFactor += LogFactorials[TmpOccupationNumberUp[k]];
	  TmpFactor += LogFactorials[TmpOccupationNumberDown[k]];
	}
      TmpSubsytemLogFactorials[i] = TmpFactor;      
    }
  int MaxIndex = minIndex + nbrIndex;
  for (int MinIndex = minIndex; MinIndex < MaxIndex; ++MinIndex)    
    {
      ComplementarySubsytemSpace->ConvertToMonomial(ComplementarySubsytemSpace->StateDescriptionUp[MinIndex], ComplementarySubsytemSpace->StateDescriptionDown[MinIndex], 
						   TmpMonomialUp1, TmpMonomialDown1);
       for (int k = 0; k < ComplementarySubsytemSpace->NbrBosons; ++k)
	 {
	   TmpMonomialUp1[k] += this->LzMax + 1 - NbrOrbitalA;
	   TmpMonomialDown1[k] += this->LzMax + 1 - NbrOrbitalA;
	 }   
      ComplementarySubsytemSpace->FermionToBoson(ComplementarySubsytemSpace->StateDescriptionUp[MinIndex], ComplementarySubsytemSpace->StateDescriptionDown[MinIndex],  
						 ComplementarySubsytemSpace->TemporaryStateUp, ComplementarySubsytemSpace->TemporaryStateDown);
      double ComplementarySubsytemSpaceFactorial = 0.0;
      for (int k = 0; k <= ComplementarySubsytemSpace->LzMax; ++k)
	{
	  ComplementarySubsytemSpaceFactorial += LogFactorials[ComplementarySubsytemSpace->TemporaryStateUp[k]];
	  ComplementarySubsytemSpaceFactorial += LogFactorials[ComplementarySubsytemSpace->TemporaryStateDown[k]];
	}
      for (int j = 0; j < SubsytemSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long* TmpMonomialUp2 = TmpSubsytemSpaceMonomialUp[j];
	  unsigned long* TmpMonomialDown2 = TmpSubsytemSpaceMonomialDown[j];
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementarySubsytemSpace->NbrBosonsUp) && (TmpIndex3 < SubsytemSpace->NbrBosonsUp)) 
	    {
	      while ((TmpIndex2 < ComplementarySubsytemSpace->NbrBosonsUp) && (TmpMonomialUp2[TmpIndex3] <= TmpMonomialUp1[TmpIndex2]))
		{
		  TmpMonomialUp3[TmpIndex4] = TmpMonomialUp1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementarySubsytemSpace->NbrBosonsUp)
		{
		  while ((TmpIndex3 < SubsytemSpace->NbrBosonsUp) && (TmpMonomialUp1[TmpIndex2] <= TmpMonomialUp2[TmpIndex3]))
		    {
		      TmpMonomialUp3[TmpIndex4] = TmpMonomialUp2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementarySubsytemSpace->NbrBosonsUp)
	    {
	      TmpMonomialUp3[TmpIndex4] = TmpMonomialUp1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < SubsytemSpace->NbrBosonsUp)
	    {
	      TmpMonomialUp3[TmpIndex4] = TmpMonomialUp2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }
	  TmpIndex2 = 0;
	  TmpIndex3 = 0;
	  TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementarySubsytemSpace->NbrBosonsDown) && (TmpIndex3 < SubsytemSpace->NbrBosonsDown)) 
	    {
	      while ((TmpIndex2 < ComplementarySubsytemSpace->NbrBosonsDown) && (TmpMonomialDown2[TmpIndex3] <= TmpMonomialDown1[TmpIndex2]))
		{
		  TmpMonomialDown3[TmpIndex4] = TmpMonomialDown1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementarySubsytemSpace->NbrBosonsDown)
		{
		  while ((TmpIndex3 < SubsytemSpace->NbrBosonsDown) && (TmpMonomialDown1[TmpIndex2] <= TmpMonomialDown2[TmpIndex3]))
		    {
		      TmpMonomialDown3[TmpIndex4] = TmpMonomialDown2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementarySubsytemSpace->NbrBosonsDown)
	    {
	      TmpMonomialDown3[TmpIndex4] = TmpMonomialDown1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < SubsytemSpace->NbrBosonsDown)
	    {
	      TmpMonomialDown3[TmpIndex4] = TmpMonomialDown2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpStateUp;
	  unsigned long TmpStateDown;
	  this->ConvertFromMonomial(TmpMonomialUp3, TmpMonomialDown3, TmpStateUp, TmpStateDown);
	  int TmpPos = this->FindStateIndex(TmpStateUp, TmpStateDown);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpStateUp, TmpStateDown, this->TemporaryStateUp, this->TemporaryStateDown);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->LzMax; ++k)
		{
		  TmpFactorial += LogFactorials[this->TemporaryStateUp[k]];
		  TmpFactorial += LogFactorials[this->TemporaryStateDown[k]];
		}
	      TmpFactorial -= ComplementarySubsytemSpaceFactorial + TmpSubsytemLogFactorials[j] + TmpLogBinomial;
	      TmpFactorial *= 0.5; 	      
	      ++TmpNbrNonZeroElements;
	      double Tmp = exp(TmpFactorial) * groundState[TmpPos];
	      entanglementMatrix->SetMatrixElement(j, MinIndex, Tmp);
	    }
	}
    }

  for (int i = 0; i < SubsytemSpace->HilbertSpaceDimension; ++i)
    {
      delete[] TmpSubsytemSpaceOccupationNumbersUp[i];
      delete[] TmpSubsytemSpaceOccupationNumbersDown[i];
      delete[] TmpSubsytemSpaceMonomialUp[i];
      delete[] TmpSubsytemSpaceMonomialDown[i];
    }
  delete[] TmpSubsytemSpaceOccupationNumbersUp;
  delete[] TmpSubsytemSpaceOccupationNumbersDown;
  delete[] TmpSubsytemSpaceMonomialUp;
  delete[] TmpSubsytemSpaceMonomialDown;
  delete[] LogFactorials;
  delete SubsytemSpace;
  delete ComplementarySubsytemSpace;
  return TmpNbrNonZeroElements;
}
   

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
// weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
// weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& BosonOnSphereWithSU2Spin::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int lzSector, int szSector,
															 int nbrOrbitalA, double* weightOrbitalAUp, double* weightOrbitalADown, 
															 int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown, RealMatrix& entanglementMatrix)
{
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  int TotalLzDisk = ConvertLzFromSphereToDisk(this->TotalLz, this->NbrBosons, this->LzMax);
  int LzADisk = ConvertLzFromSphereToDisk(lzSector, nbrParticleSector, nbrOrbitalA - 1);
  if ((LzADisk < 0) || (LzADisk > ((nbrOrbitalA - 1) * nbrParticleSector)))
    {
      return entanglementMatrix;	  
    }
  int LzBDisk = (TotalLzDisk - LzADisk) - ComplementaryNbrParticles * (this->LzMax + 1 - nbrOrbitalB);
  if ((LzBDisk < 0) || (LzBDisk > ((nbrOrbitalB - 1) * ComplementaryNbrParticles)))
    {
      return entanglementMatrix;	  
    }
  int ComplementaryLzSector = ConvertLzFromDiskToSphere(LzBDisk, ComplementaryNbrParticles, nbrOrbitalB - 1);
  BosonOnSphereWithSU2Spin SubsytemSpace(nbrParticleSector, lzSector, nbrOrbitalA - 1, szSector);
  BosonOnSphereWithSU2Spin ComplementarySubsytemSpace(ComplementaryNbrParticles, ComplementaryLzSector, nbrOrbitalB - 1, ComplementarySzSector);  

  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;
  unsigned long* TmpMonomialUp1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpMonomialDown1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpMonomialUp3 = new unsigned long [nbrParticleSector];
  unsigned long* TmpMonomialDown3 = new unsigned long [nbrParticleSector];
  
  for (int i = 0; i < SubsytemSpace.HilbertSpaceDimension; ++i)
    {
      SubsytemSpace.ConvertToMonomial(SubsytemSpace.StateDescriptionUp[i], SubsytemSpace.StateDescriptionDown[i], TmpMonomialUp3, TmpMonomialDown3);
      double Tmp = 1.0;
      for (int j = 0; j < SubsytemSpace.NbrBosonsUp; j++)
	{
	  Tmp *= weightOrbitalAUp[TmpMonomialUp3[j]];
	}
      for (int j = 0; j < SubsytemSpace.NbrBosonsDown; j++)
	{
	  Tmp *= weightOrbitalADown[TmpMonomialDown3[j]];
	}
      for (int j = 0; j < ComplementarySubsytemSpace.HilbertSpaceDimension; ++j)          
	entanglementMatrix(i, j) *= Tmp;      
    }
  
  for (int MinIndex = 0; MinIndex < ComplementarySubsytemSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      ComplementarySubsytemSpace.ConvertToMonomial(ComplementarySubsytemSpace.StateDescriptionUp[MinIndex], ComplementarySubsytemSpace.StateDescriptionDown[MinIndex], TmpMonomialUp1, TmpMonomialDown1);
      double FormFactor = 1.0;
      for (int i = 0; i < ComplementarySubsytemSpace.NbrBosonsUp; i++)
	FormFactor *= weightOrbitalBUp[TmpMonomialUp1[i]];
      for (int i = 0; i < ComplementarySubsytemSpace.NbrBosonsDown; i++)
	FormFactor *= weightOrbitalBDown[TmpMonomialDown1[i]];
      for (int j = 0; j < SubsytemSpace.HilbertSpaceDimension; ++j)
	entanglementMatrix(j, MinIndex) *= FormFactor; 
    }
  
  delete[] TmpMonomialUp1;
  delete[] TmpMonomialDown1;
  delete[] TmpMonomialUp3;
  delete[] TmpMonomialDown3;
  
  return entanglementMatrix;
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition breaking the momentum conservation. 
// The entanglement matrix is computed from precalculated particle entanglement matrices in each momentum sector
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// nbrConnectedOrbitalAUp = number of orbitals connected to a given one by the A part real space cut (for the spin up)
// nbrConnectedOrbitalADown = number of orbitals connected to a given one by the A part real space cut (for the spin down)
// connectedOrbitalAUp = orbitals taht connected to a given one by the A part real space cut (for the spin up)
// connectedOrbitalADown = orbitals taht connected to a given one by the A part real space cut (for the spin down)
// weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
// weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// nbrConnectedOrbitalBUp = number of orbitals connected to a given one by the B part real space cut (for the spin up)
// nbrConnectedOrbitalBDown = number of orbitals connected to a given one by the B part real space cut (for the spin down)
// connectedOrbitalBUp = orbitals taht connected to a given one by the B part real space cut (for the spin up)
// connectedOrbitalBDown = orbitals taht connected to a given one by the B part real space cut (for the spin down)
// weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
// weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
// nbrEntanglementMatrices = number of available entanglement matrices with a fixed moementum
// entanglementMatrixLzSectors = momentum sector of each entanglement matrix
// entanglementMatrices = array containing the entanglement matrices with a fixed moementum
// return value = real space entanglement matrix

RealMatrix BosonOnSphereWithSU2Spin::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int szSector,
															int nbrOrbitalA, int* nbrConnectedOrbitalAUp, int* nbrConnectedOrbitalADown,
															int** connectedOrbitalAUp, int** connectedOrbitalADown, 
															double** weightOrbitalAUp, double** weightOrbitalADown, 
															int nbrOrbitalB, int* nbrConnectedOrbitalBUp, int* nbrConnectedOrbitalBDown, 
															int** connectedOrbitalBUp, int** connectedOrbitalBDown, 
															double** weightOrbitalBUp, double** weightOrbitalBDown, 
															int nbrEntanglementMatrices, int* entanglementMatrixLzSectors,
															RealMatrix* entanglementMatrices)
{
  int TotalNbrRow = 0;
  int TotalNbrColumn = 0;

  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;


  BosonOnSphereWithSU2SpinAllLz* TotalSubsystemSpace = new BosonOnSphereWithSU2SpinAllLz(nbrParticleSector, nbrOrbitalA - 1, szSector);
  BosonOnSphereWithSU2SpinAllLz* TotalComplementarySubsystemSpace = new BosonOnSphereWithSU2SpinAllLz(ComplementaryNbrParticles, nbrOrbitalB - 1, ComplementarySzSector);

  cout << "size of the full entanglement matrix " << TotalSubsystemSpace->GetHilbertSpaceDimension() << " x " << TotalComplementarySubsystemSpace->GetHilbertSpaceDimension()
       << " (requiring  " << (((double) TotalSubsystemSpace->GetHilbertSpaceDimension()) * ((double) TotalComplementarySubsystemSpace->GetHilbertSpaceDimension())/ 131072.0) << " Mb )" << endl;

  BosonOnSphereWithSU2Spin** SubsystemSpaces = new BosonOnSphereWithSU2Spin* [nbrEntanglementMatrices];
  BosonOnSphereWithSU2Spin** ComplementarySubsystemSpaces = new BosonOnSphereWithSU2Spin* [nbrEntanglementMatrices];
  int TotalLzDisk = ConvertLzFromSphereToDisk(this->TotalLz, this->NbrBosons, this->LzMax);  
  
  for (int i = 0; i < nbrEntanglementMatrices; ++i)
    {
      int LzADisk = ConvertLzFromSphereToDisk(entanglementMatrixLzSectors[i], nbrParticleSector, nbrOrbitalA - 1);
      int LzBDisk = (TotalLzDisk - LzADisk) - ComplementaryNbrParticles * (this->LzMax + 1 - nbrOrbitalB);
      int ComplementaryLzSector = ConvertLzFromDiskToSphere(LzBDisk, ComplementaryNbrParticles, nbrOrbitalB - 1);
      SubsystemSpaces[i] = new BosonOnSphereWithSU2Spin(nbrParticleSector, entanglementMatrixLzSectors[i], nbrOrbitalA - 1, szSector);
      ComplementarySubsystemSpaces[i] = new BosonOnSphereWithSU2Spin(ComplementaryNbrParticles, ComplementaryLzSector, nbrOrbitalB - 1, ComplementarySzSector);
    }

  RealMatrix TmpEntanglementMatrix (TotalSubsystemSpace->GetHilbertSpaceDimension(), TotalComplementarySubsystemSpace->GetHilbertSpaceDimension(), true);

  unsigned long* TmpComplementaryMonomialUp1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpComplementaryMonomialDown1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpComplementaryMonomialUp2 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpComplementaryMonomialDown2 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpMonomialUp1 = new unsigned long [nbrParticleSector];
  unsigned long* TmpMonomialDown1 = new unsigned long [nbrParticleSector];
  unsigned long* TmpMonomialUp2 = new unsigned long [nbrParticleSector];
  unsigned long* TmpMonomialDown2 = new unsigned long [nbrParticleSector];

  int SubsystemNbrUp = SubsystemSpaces[0]->NbrBosonsUp;
  int SubsystemNbrDown = SubsystemSpaces[0]->NbrBosonsDown;
  int ComplementarySubsystemNbrUp = ComplementarySubsystemSpaces[0]->NbrBosonsUp;
  int ComplementarySubsystemNbrDown = ComplementarySubsystemSpaces[0]->NbrBosonsDown;
  RealMatrix TmpSubsystemPermanentMatrixUp (SubsystemNbrUp, SubsystemNbrUp);
  RealMatrix TmpSubsystemPermanentMatrixDown (SubsystemNbrDown, SubsystemNbrDown);
  RealMatrix TmpComplementarySubsystemPermanentMatrixUp (ComplementarySubsystemNbrUp, ComplementarySubsystemNbrUp);
  RealMatrix TmpComplementarySubsystemPermanentMatrixDown (ComplementarySubsystemNbrDown, ComplementarySubsystemNbrDown);

  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 

  double* TotalSubsystemOccupationFactors = new double[TotalSubsystemSpace->GetHilbertSpaceDimension()];
  unsigned long* TmpSubsystemSpaceOccupationNumbersUp = new unsigned long [TotalSubsystemSpace->NbrLzValue];
  unsigned long* TmpSubsystemSpaceOccupationNumbersDown = new unsigned long [TotalSubsystemSpace->NbrLzValue];
  for (int i = 0; i < TotalSubsystemSpace->HilbertSpaceDimension; ++i)
    {
      TotalSubsystemSpace->FermionToBoson(TotalSubsystemSpace->StateDescriptionUp[i], TotalSubsystemSpace->StateDescriptionDown[i], 
					  TmpSubsystemSpaceOccupationNumbersUp, TmpSubsystemSpaceOccupationNumbersDown);
      double TmpFactor = 0.0;
      for (int k = 0; k <= TotalSubsystemSpace->LzMax; ++k)
	{
	  TmpFactor += LogFactorials[TmpSubsystemSpaceOccupationNumbersUp[k]];
	  TmpFactor += LogFactorials[TmpSubsystemSpaceOccupationNumbersDown[k]];
	}
      TotalSubsystemOccupationFactors[i] = TmpFactor;      
    }
  double* TotalComplementarySubsystemOccupationFactors = new double[TotalComplementarySubsystemSpace->GetHilbertSpaceDimension()];
  unsigned long* TmpComplementarySubsystemSpaceOccupationNumbersUp = new unsigned long [TotalComplementarySubsystemSpace->NbrLzValue];
  unsigned long* TmpComplementarySubsystemSpaceOccupationNumbersDown = new unsigned long [TotalComplementarySubsystemSpace->NbrLzValue];
  for (int i = 0; i <  TotalComplementarySubsystemSpace->GetHilbertSpaceDimension(); ++i)
    {
      TotalComplementarySubsystemSpace->FermionToBoson(TotalComplementarySubsystemSpace->StateDescriptionUp[i], TotalComplementarySubsystemSpace->StateDescriptionDown[i], 
						      TmpComplementarySubsystemSpaceOccupationNumbersUp, TmpComplementarySubsystemSpaceOccupationNumbersDown);
      double TmpFactor = 0.0;
      for (int k = 0; k <= TotalComplementarySubsystemSpace->LzMax; ++k)
	{
	  TmpFactor += LogFactorials[TmpComplementarySubsystemSpaceOccupationNumbersUp[k]];
	  TmpFactor += LogFactorials[TmpComplementarySubsystemSpaceOccupationNumbersDown[k]];
	}
      TotalComplementarySubsystemOccupationFactors[i] = TmpFactor;      
    }
  double** TotalAllSubsystemOccupationFactors = new double*[nbrEntanglementMatrices];
  double** TotalAllComplementarySubsystemOccupationFactors = new double*[nbrEntanglementMatrices];
  for (int i = 0; i < nbrEntanglementMatrices; ++i)
    {
      TotalAllComplementarySubsystemOccupationFactors[i] = new double[ComplementarySubsystemSpaces[i]->GetHilbertSpaceDimension()];
      for (int TmpComplementarySubsystemIndex = 0; TmpComplementarySubsystemIndex <  ComplementarySubsystemSpaces[i]->GetHilbertSpaceDimension(); ++TmpComplementarySubsystemIndex)
	{
	  ComplementarySubsystemSpaces[i]->FermionToBoson(ComplementarySubsystemSpaces[i]->StateDescriptionUp[TmpComplementarySubsystemIndex], 
							 ComplementarySubsystemSpaces[i]->StateDescriptionDown[TmpComplementarySubsystemIndex], 
							 TmpComplementarySubsystemSpaceOccupationNumbersUp, TmpComplementarySubsystemSpaceOccupationNumbersDown);
	  double TmpFactor = 0.0;
	  for (int k = 0; k <= ComplementarySubsystemSpaces[i]->LzMax; ++k)
	    {
	      TmpFactor += LogFactorials[TmpComplementarySubsystemSpaceOccupationNumbersUp[k]];
	      TmpFactor += LogFactorials[TmpComplementarySubsystemSpaceOccupationNumbersDown[k]];
	    }
	  TotalAllComplementarySubsystemOccupationFactors[i][TmpComplementarySubsystemIndex] = TmpFactor;      
	}
      TotalAllSubsystemOccupationFactors[i] = new double[SubsystemSpaces[i]->GetHilbertSpaceDimension()];
      for (int TmpSubsystemIndex = 0; TmpSubsystemIndex <  SubsystemSpaces[i]->GetHilbertSpaceDimension(); ++TmpSubsystemIndex)
	{
	  SubsystemSpaces[i]->FermionToBoson(SubsystemSpaces[i]->StateDescriptionUp[TmpSubsystemIndex], 
					    SubsystemSpaces[i]->StateDescriptionDown[TmpSubsystemIndex], 
					    TmpSubsystemSpaceOccupationNumbersUp, TmpSubsystemSpaceOccupationNumbersDown);
	  double TmpFactor = 0.0;
	  for (int k = 0; k <= SubsystemSpaces[i]->LzMax; ++k)
	    {
	      TmpFactor += LogFactorials[TmpSubsystemSpaceOccupationNumbersUp[k]];
	      TmpFactor += LogFactorials[TmpSubsystemSpaceOccupationNumbersDown[k]];
	    }
	  TotalAllSubsystemOccupationFactors[i][TmpSubsystemIndex] = TmpFactor;      
	}
    }



  for (int i = 0; i < nbrEntanglementMatrices; ++i)
    {
      cout << "size of the B tranformation matrix " << ComplementarySubsystemSpaces[i]->GetHilbertSpaceDimension() << " x " << TotalComplementarySubsystemSpace->GetHilbertSpaceDimension()
	   << " (requiring  " << (((double) ComplementarySubsystemSpaces[i]->GetHilbertSpaceDimension()) * ((double) TotalComplementarySubsystemSpace->GetHilbertSpaceDimension())/ 131072.0) << " Mb )" << endl;
      RealMatrix TmpComplementaryTransformationMatrix (ComplementarySubsystemSpaces[i]->GetHilbertSpaceDimension(), TotalComplementarySubsystemSpace->GetHilbertSpaceDimension(), true);
      for (int TmpComplementarySubsystemIndex = 0; TmpComplementarySubsystemIndex <  ComplementarySubsystemSpaces[i]->GetHilbertSpaceDimension(); ++TmpComplementarySubsystemIndex)
	{
	  ComplementarySubsystemSpaces[i]->ConvertToMonomial(ComplementarySubsystemSpaces[i]->StateDescriptionUp[TmpComplementarySubsystemIndex],
							    ComplementarySubsystemSpaces[i]->StateDescriptionDown[TmpComplementarySubsystemIndex],
							    TmpComplementaryMonomialUp1, TmpComplementaryMonomialDown1);
	  for (int TmpTotalComplementarySubsystemIndex = 0; TmpTotalComplementarySubsystemIndex <  TotalComplementarySubsystemSpace->GetHilbertSpaceDimension(); ++TmpTotalComplementarySubsystemIndex)
	    {
	      TotalComplementarySubsystemSpace->ConvertToMonomial(TotalComplementarySubsystemSpace->StateDescriptionUp[TmpTotalComplementarySubsystemIndex],
								 TotalComplementarySubsystemSpace->StateDescriptionDown[TmpTotalComplementarySubsystemIndex],
								 TmpComplementaryMonomialUp2, TmpComplementaryMonomialDown2);
	      double Tmp = 1.0;
	      if (ComplementarySubsystemNbrUp > 0)
		{
		  TmpComplementarySubsystemPermanentMatrixUp.ClearMatrix();
		  for (int j = 0; j < ComplementarySubsystemNbrUp; ++j)
		    {
		      for (int k = 0; k < ComplementarySubsystemNbrUp; ++k)
			{
			  int TmpIndex = SearchInArray<int>(TmpComplementaryMonomialUp2[k], connectedOrbitalBUp[TmpComplementaryMonomialUp1[j]], 
							    nbrConnectedOrbitalBUp[TmpComplementaryMonomialUp1[j]]);
			  if (TmpIndex >= 0)
			    {
			      TmpComplementarySubsystemPermanentMatrixUp.SetMatrixElement(j, k, weightOrbitalBUp[TmpComplementaryMonomialUp1[j]][TmpIndex]);
			    }
			}
		    }
		  Tmp *= TmpComplementarySubsystemPermanentMatrixUp.Permanent();
		}
	      if (ComplementarySubsystemNbrDown > 0)
		{
		  TmpComplementarySubsystemPermanentMatrixDown.ClearMatrix();
		  for (int j = 0; j < ComplementarySubsystemNbrDown; ++j)
		    {
		      for (int k = 0; k < ComplementarySubsystemNbrDown; ++k)
			{
			  int TmpIndex = SearchInArray<int>(TmpComplementaryMonomialDown2[k], connectedOrbitalBDown[TmpComplementaryMonomialDown1[j]], 
							    nbrConnectedOrbitalBDown[TmpComplementaryMonomialDown1[j]]);
			  if (TmpIndex >= 0)
			    {
			      TmpComplementarySubsystemPermanentMatrixDown.SetMatrixElement(j, k, weightOrbitalBDown[TmpComplementaryMonomialDown1[j]][TmpIndex]);
			    }
			}
		    }
		   Tmp *= TmpComplementarySubsystemPermanentMatrixDown.Permanent();
		}
	      Tmp *= exp(-0.5 * (TotalComplementarySubsystemOccupationFactors[TmpTotalComplementarySubsystemIndex] 
				 +  + TotalAllComplementarySubsystemOccupationFactors[i][TmpComplementarySubsystemIndex]));
	      TmpComplementaryTransformationMatrix.SetMatrixElement(TmpComplementarySubsystemIndex, TmpTotalComplementarySubsystemIndex, Tmp);
	    }
	}
      cout << "size of the A tranformation matrix " << SubsystemSpaces[i]->GetHilbertSpaceDimension() << " x " << TotalSubsystemSpace->GetHilbertSpaceDimension()
	   << " (requiring  " << (((double) SubsystemSpaces[i]->GetHilbertSpaceDimension()) * ((double) TotalSubsystemSpace->GetHilbertSpaceDimension())/ 131072.0) << " Mb )" << endl;
      RealMatrix TmpTransformationMatrix (TotalSubsystemSpace->GetHilbertSpaceDimension(), SubsystemSpaces[i]->GetHilbertSpaceDimension(), true);
      for (int TmpSubsystemIndex = 0; TmpSubsystemIndex <  SubsystemSpaces[i]->GetHilbertSpaceDimension(); ++TmpSubsystemIndex)
	{
	  SubsystemSpaces[i]->ConvertToMonomial(SubsystemSpaces[i]->StateDescriptionUp[TmpSubsystemIndex],
					      SubsystemSpaces[i]->StateDescriptionDown[TmpSubsystemIndex],
					      TmpMonomialUp1, TmpMonomialDown1);
	  for (int TmpTotalSubsystemIndex = 0; TmpTotalSubsystemIndex <  TotalSubsystemSpace->GetHilbertSpaceDimension(); ++TmpTotalSubsystemIndex)
	    {
	      TotalSubsystemSpace->ConvertToMonomial(TotalSubsystemSpace->StateDescriptionUp[TmpTotalSubsystemIndex],
						    TotalSubsystemSpace->StateDescriptionDown[TmpTotalSubsystemIndex],
						    TmpMonomialUp2, TmpMonomialDown2);
	      double Tmp = 1.0;
	      if (SubsystemNbrUp > 0)
		{
		  TmpSubsystemPermanentMatrixUp.ClearMatrix();
		  for (int j = 0; j < SubsystemNbrUp; ++j)
		    {
		      for (int k = 0; k < SubsystemNbrUp; ++k)
			{
			  int TmpIndex = SearchInArray<int>(TmpMonomialUp2[k], connectedOrbitalAUp[TmpMonomialUp1[j]], 
							    nbrConnectedOrbitalAUp[TmpMonomialUp1[j]]);
			  if (TmpIndex >= 0)
			    {
			      TmpSubsystemPermanentMatrixUp.SetMatrixElement(j, k, weightOrbitalAUp[TmpMonomialUp1[j]][TmpIndex]);
			    }
			}
		    }
		  Tmp *= TmpSubsystemPermanentMatrixUp.Permanent();
		}
	      if (SubsystemNbrDown > 0)
		{
		  TmpSubsystemPermanentMatrixDown.ClearMatrix();
		  for (int j = 0; j < SubsystemNbrDown; ++j)
		    {
		      for (int k = 0; k < SubsystemNbrDown; ++k)
			{
			  int TmpIndex = SearchInArray<int>(TmpMonomialDown2[k], connectedOrbitalADown[TmpMonomialDown1[j]], 
							    nbrConnectedOrbitalADown[TmpMonomialDown1[j]]);
			  if (TmpIndex >= 0)
			    {
			      TmpSubsystemPermanentMatrixDown.SetMatrixElement(j, k, weightOrbitalADown[TmpMonomialDown1[j]][TmpIndex]);
			    }
			}
		    }
		  Tmp *= TmpSubsystemPermanentMatrixDown.Permanent();
		}
	      Tmp *= exp(-0.5 * (TotalSubsystemOccupationFactors[TmpTotalSubsystemIndex] + TotalAllSubsystemOccupationFactors[i][TmpSubsystemIndex]));
	      TmpTransformationMatrix.SetMatrixElement(TmpTotalSubsystemIndex, TmpSubsystemIndex, Tmp);
	    }
	}
//       cout << "processing " << i << endl;
//       cout << entanglementMatrices[i] << endl;
//       cout << TmpComplementaryTransformationMatrix << endl;
//       cout << TmpTransformationMatrix << endl;
      RealMatrix TmpMatrix = entanglementMatrices[i] * TmpComplementaryTransformationMatrix;
      TmpEntanglementMatrix.AddMultiply (TmpTransformationMatrix, TmpMatrix);
    }
  
  for (int i = 0; i < nbrEntanglementMatrices; ++i)
    {
      delete SubsystemSpaces[i];
      delete ComplementarySubsystemSpaces[i];
      delete[] TotalAllSubsystemOccupationFactors[i];
      delete[] TotalAllComplementarySubsystemOccupationFactors[i];
    }
  delete[] SubsystemSpaces;
  delete[] ComplementarySubsystemSpaces;
  delete[] TotalSubsystemOccupationFactors;
  delete[] TotalComplementarySubsystemOccupationFactors;
  delete[] TotalAllSubsystemOccupationFactors;
  delete[] TotalAllComplementarySubsystemOccupationFactors;
  delete[] LogFactorials;

  return TmpEntanglementMatrix;
}



// convert a given state from a generic basis from the current Sz subspace basis
//
// state = reference on the vector to convert
// space = reference on the basis associated to state
// return value = converted vector

RealVector BosonOnSphereWithSU2Spin::ConvertToNbodyBasis(RealVector& state, ParticleOnSphereWithSpin* space)
{
  RealVector TmpVector (state, true);
  return TmpVector;
}
  
// convert a given state from a generic basis to the current Sz subspace basis
//
// state = reference on the vector to convert
// space = reference on the basis associated to state
// return value = converted vector

RealVector BosonOnSphereWithSU2Spin::ConvertFromNbodyBasis(RealVector& state, ParticleOnSphereWithSpin* space)
{
  RealVector TmpVector (state, true);
  return TmpVector;
}
  
// symmetrized a product of two decoupled states, core part
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// firstComponent = index of the first component
// nbrComponents = number of components to symmetrize
// return value = symmetrized state

void BosonOnSphereWithSU2Spin::SymmetrizeSU2SU2StateCore (RealVector& symmetrizedVector, RealVector& leftVector, RealVector& rightVector, 
							  ParticleOnSphereWithSpin* leftSpace, ParticleOnSphereWithSpin* rightSpace, 
							  bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents)
{
  BosonOnSphereWithSU2Spin* TmpRightSpace = (BosonOnSphereWithSU2Spin*) rightSpace;
  BosonOnSphereWithSU2Spin* TmpLeftSpace = (BosonOnSphereWithSU2Spin*) leftSpace;
  long LastComponent = long(firstComponent + nbrComponents);
  double* LogFactorials = new double [this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2; i <= this->NbrBosons; ++i)
    {
      LogFactorials[i] = LogFactorials[i - 1] + log((double) i);
    }
  unsigned long TmpFinalStateUp;
  unsigned long TmpFinalStateDown;

  for (long i = (long) firstComponent; i < LastComponent; ++i)
    {
      TmpLeftSpace->FermionToBoson(TmpLeftSpace->StateDescriptionUp[i], TmpLeftSpace->StateDescriptionDown[i], 
				   TmpLeftSpace->TemporaryStateUp, TmpLeftSpace->TemporaryStateDown);
      double TmpLogFactor1 = 0.0;
      for (int k = 0; k <= TmpLeftSpace->LzMax; ++k)
	{
	  TmpLogFactor1 -= LogFactorials[TmpLeftSpace->TemporaryStateUp[k]];
	  TmpLogFactor1 -= LogFactorials[TmpLeftSpace->TemporaryStateDown[k]];
	}
      double TmpCoefficient = leftVector[i];
      for (long j = 0l; j < TmpRightSpace->LargeHilbertSpaceDimension; ++j)
	{
	  TmpRightSpace->FermionToBoson(TmpRightSpace->StateDescriptionUp[j], TmpRightSpace->StateDescriptionDown[j], 
					TmpRightSpace->TemporaryStateUp, TmpRightSpace->TemporaryStateDown);
	  for (int k = 0; k <= TmpRightSpace->LzMax; ++k)
	    {
	      this->TemporaryStateUp[k] = TmpRightSpace->TemporaryStateUp[k] + TmpLeftSpace->TemporaryStateUp[k];
	      this->TemporaryStateDown[k] = TmpRightSpace->TemporaryStateDown[k] + TmpLeftSpace->TemporaryStateDown[k];
	    }
	  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpFinalStateUp, TmpFinalStateDown);
	  int TmpPos = this->FindStateIndex(TmpFinalStateUp, TmpFinalStateDown);
	  if (TmpPos < this->HilbertSpaceDimension)
	    {
	      double TmpLogFactor2 = 0.0;
	      for (int k = 0; k <= TmpRightSpace->LzMax; ++k)
		{
		  TmpLogFactor2 -= LogFactorials[TmpRightSpace->TemporaryStateUp[k]];
		  TmpLogFactor2 -= LogFactorials[TmpRightSpace->TemporaryStateDown[k]];
		}
	      for (int k = 0; k <= this->LzMax; ++k)
		{
		  TmpLogFactor2 += LogFactorials[this->TemporaryStateUp[k]];
		  TmpLogFactor2 += LogFactorials[this->TemporaryStateDown[k]];
		}
	      symmetrizedVector[TmpPos] += exp(0.5 * (TmpLogFactor1 + TmpLogFactor2)) * TmpCoefficient * rightVector[j];
	    }
	}
    }
  delete[] LogFactorials;
}

