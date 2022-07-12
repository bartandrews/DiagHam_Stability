////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                    class of bosons on sphere with SU(4) spin               //
//                                                                            //
//                        last modification : 19/12/2011                      //
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
#include "HilbertSpace/BosonOnSphereWithSU4SpinLong.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereLong.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <math.h>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
// 

BosonOnSphereWithSU4SpinLong::BosonOnSphereWithSU4SpinLong ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// totalSpin = twice the total spin value
// totalIsospin = twice the total isospin value
// totalEntanglement = twice the total entanglement value
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU4SpinLong::BosonOnSphereWithSU4SpinLong (int nbrBosons, int totalLz, int lzMax, int totalSpin, int totalIsospin, int totalEntanglement,
						    unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->TotalIsospin = totalIsospin;
  this->TotalEntanglement = totalEntanglement;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUpPlus;
  this->TemporaryStateSigma[1] = this->TemporaryStateUpMinus;
  this->TemporaryStateSigma[2] = this->TemporaryStateDownPlus;
  this->TemporaryStateSigma[3] = this->TemporaryStateDownMinus;
  this->ProdATemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUpPlus;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateUpMinus;
  this->ProdATemporaryStateSigma[2] = this->ProdATemporaryStateDownPlus;
  this->ProdATemporaryStateSigma[3] = this->ProdATemporaryStateDownMinus;

  int NUpPlus = this->NbrBosons + this->TotalSpin + this->TotalIsospin + this->TotalEntanglement;
  int NUpMinus = this->NbrBosons + this->TotalSpin - this->TotalIsospin - this->TotalEntanglement;
  int NDownPlus = this->NbrBosons - this->TotalSpin + this->TotalIsospin - this->TotalEntanglement;
  int NDownMinus = this->NbrBosons - this->TotalSpin - this->TotalIsospin + this->TotalEntanglement;
  this->NUpPlusLzMax = this->LzMax + NUpPlus - 1;
  this->NUpMinusLzMax = this->LzMax + NUpMinus - 1;
  this->NDownPlusLzMax = this->LzMax + NDownPlus - 1;
  this->NDownMinusLzMax = this->LzMax + NDownMinus - 1;
  this->FermionicLzMax = this->NUpPlusLzMax;
  if (this->NUpMinusLzMax > this->FermionicLzMax)
    this->FermionicLzMax = this->NUpMinusLzMax;
  if (this->NDownPlusLzMax > this->FermionicLzMax)
    this->FermionicLzMax = this->NDownPlusLzMax;
  if (this->NDownMinusLzMax > this->FermionicLzMax)
    this->FermionicLzMax = this->NDownMinusLzMax;
  if ((NUpPlus < 0) || ((NUpPlus & 0x3) != 0) || (NUpMinus < 0) || ((NUpMinus & 0x3) != 0) ||
      (NDownPlus < 0) || ((NDownPlus & 0x3) != 0) || (NDownMinus < 0) || ((NDownMinus & 0x3) != 0))
    this->LargeHilbertSpaceDimension = 0l;
  else
    {
      NUpPlus >>= 2;
      NUpMinus >>= 2;
      NDownPlus >>= 2;
      NDownMinus >>= 2;
      this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, NUpPlus, NUpMinus, NDownPlus, NDownMinus);
    }
  this->StateDescriptionUpPlus = new ULONGLONG [this->LargeHilbertSpaceDimension];
  this->StateDescriptionUpMinus = new ULONGLONG [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDownPlus = new ULONGLONG [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDownMinus = new ULONGLONG [this->LargeHilbertSpaceDimension];
  this->StateDescriptionSigma[0] = this->StateDescriptionUpPlus;
  this->StateDescriptionSigma[1] = this->StateDescriptionUpMinus;
  this->StateDescriptionSigma[2] = this->StateDescriptionDownPlus;
  this->StateDescriptionSigma[3] = this->StateDescriptionDownMinus;
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, this->LzMax, this->LzMax, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 
						       NUpPlus, NUpMinus, NDownPlus, NDownMinus, 0l);
//   for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
//     cout << i << " : " << hex << this->StateDescriptionUpPlus[i] << " " << this->StateDescriptionUpMinus[i] << " " << this->StateDescriptionDownPlus[i] << dec << endl;
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU4Spin!" << endl;
      exit(1);
    }
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  



  this->GenerateLookUpTable(memory);
//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)	
//     {
//       cout << i << " : ";
//       this->PrintState(cout, i);
//       cout << this->FindStateIndex(this->StateDescriptionUpPlus[i], this->StateDescriptionUpMinus[i], this->StateDescriptionDownPlus[i]);
//       cout << endl;
//     }
#ifdef __DEBUG__
   int UsedMemory = 0;
   UsedMemory += this->HilbertSpaceDimension * (4 * sizeof(ULONGLONG));
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

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSphereWithSU4SpinLong::BosonOnSphereWithSU4SpinLong(const BosonOnSphereWithSU4SpinLong& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpPlusLzMax = bosons.LzMax;
  this->NUpMinusLzMax = bosons.LzMax;
  this->NDownPlusLzMax = bosons.LzMax;
  this->NDownMinusLzMax = bosons.LzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->TotalSpin = bosons.TotalSpin;
  this->TotalIsospin = bosons.TotalIsospin;
  this->TotalEntanglement = bosons.TotalEntanglement;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUpPlus;
  this->TemporaryStateSigma[1] = this->TemporaryStateUpMinus;
  this->TemporaryStateSigma[2] = this->TemporaryStateDownPlus;
  this->TemporaryStateSigma[3] = this->TemporaryStateDownMinus;
  this->ProdATemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
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

BosonOnSphereWithSU4SpinLong::~BosonOnSphereWithSU4SpinLong ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescriptionUpPlus;
      delete[] this->StateDescriptionUpMinus;
      delete[] this->StateDescriptionDownPlus;
      delete[] this->StateDescriptionDownMinus;
      delete[] this->UniqueStateDescriptionUpPlus;
      delete[] this->UniqueStateDescriptionSubArraySizeUpPlus;
      for (long i = 0l; i < this->NbrUniqueStateDescriptionUpPlus; ++i)
	{
	  for (long j = 0l; j < this->NbrUniqueStateDescriptionUpMinus[i]; ++j)
	    {
	      delete[] this->UniqueStateDescriptionDownPlus[i][j];
	      delete[] this->UniqueStateDescriptionSubArraySizeDownPlus[i][j];
	      delete[] this->FirstIndexUniqueStateDescriptionDownPlus[i][j];
	    }
	  delete[] this->UniqueStateDescriptionUpMinus[i];
	  delete[] this->UniqueStateDescriptionSubArraySizeUpMinus[i];
	  delete[] this->FirstIndexUniqueStateDescriptionUpMinus[i];
	  delete[] this->NbrUniqueStateDescriptionDownPlus[i];
	  delete[] this->UniqueStateDescriptionDownPlus[i];
	  delete[] this->UniqueStateDescriptionSubArraySizeDownPlus[i];
	  delete[] this->FirstIndexUniqueStateDescriptionDownPlus[i];
	}
      delete[] this->NbrUniqueStateDescriptionUpMinus;
      delete[] this->UniqueStateDescriptionUpMinus;
      delete[] this->UniqueStateDescriptionSubArraySizeUpMinus;
      delete[] this->FirstIndexUniqueStateDescriptionUpMinus;
      delete[] this->NbrUniqueStateDescriptionDownPlus;
      delete[] this->UniqueStateDescriptionDownPlus;
      delete[] this->UniqueStateDescriptionSubArraySizeDownPlus;
      delete[] this->FirstIndexUniqueStateDescriptionDownPlus;

    }
  delete[] this->TemporaryStateUpPlus;
  delete[] this->TemporaryStateUpMinus;
  delete[] this->TemporaryStateDownPlus;
  delete[] this->TemporaryStateDownMinus;
  delete[] this->ProdATemporaryStateUpPlus;
  delete[] this->ProdATemporaryStateUpMinus;
  delete[] this->ProdATemporaryStateDownPlus;
  delete[] this->ProdATemporaryStateDownMinus;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereWithSU4SpinLong& BosonOnSphereWithSU4SpinLong::operator = (const BosonOnSphereWithSU4SpinLong& bosons)
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
  this->NUpPlusLzMax = bosons.LzMax;
  this->NUpMinusLzMax = bosons.LzMax;
  this->NDownPlusLzMax = bosons.LzMax;
  this->NDownMinusLzMax = bosons.LzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUpPlus;
  this->TemporaryStateSigma[1] = this->TemporaryStateUpMinus;
  this->TemporaryStateSigma[2] = this->TemporaryStateDownPlus;
  this->TemporaryStateSigma[3] = this->TemporaryStateDownMinus;
  this->ProdATemporaryStateUpPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->NbrLzValue];
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

AbstractHilbertSpace* BosonOnSphereWithSU4SpinLong::Clone()
{
  return new BosonOnSphereWithSU4SpinLong(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereWithSU4SpinLong::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereWithSU4SpinLong::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereWithSU4SpinLong::ExtractSubspace (AbstractQuantumNumber& q, 
								 SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m_up a_m_up operator to a given state  (only spin up isospin plus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_um a_m_um

double  BosonOnSphereWithSU4SpinLong::AdupAup (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->TemporaryStateUpPlus);
  return (double) (this->TemporaryStateUpPlus[m]);  
}

// apply a^+_m_um a_m_um operator to a given state  (only spin up isospin minus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_um a_m_um

double BosonOnSphereWithSU4SpinLong::AdumAum (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->TemporaryStateUpMinus);
  return (double) (this->TemporaryStateUpMinus[m]);  
}
 
// apply a^+_m_dp a_m_dp operator to a given state (only spin down isospin plus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_dm a_m_dm

double BosonOnSphereWithSU4SpinLong::AddpAdp (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->TemporaryStateDownPlus);
  return (double) (this->TemporaryStateDownPlus[m]);  
}

// apply a^+_m_dm a_m_dm operator to a given state (only spin down isospin minus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_dm a_m_dm

double BosonOnSphereWithSU4SpinLong::AddmAdm (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->TemporaryStateDownMinus);
  return (double) (this->TemporaryStateDownMinus[m]);  
}

// find state index
//
// stateDescriptionUpPlus = unsigned integer describing the fermionic state for type up-plus particles
// stateDescriptionUpMinus = unsigned integer describing the fermionic state for type up-minus particles
// stateDescriptionDownPlus = unsigned integer describing the fermionic state for type down-plus particles
// stateDescriptionDownMinus = unsigned integer describing the fermionic state for type down-minus particles
// return value = corresponding index

int BosonOnSphereWithSU4SpinLong::FindStateIndex(ULONGLONG stateDescriptionUpPlus, ULONGLONG stateDescriptionUpMinus, 
					     ULONGLONG stateDescriptionDownPlus, ULONGLONG stateDescriptionDownMinus)
{
  int PosMin = 0;
  int PosMax = this->NbrUniqueStateDescriptionUpPlus - 1;
  int PosMidUpPlus = (PosMin + PosMax) >> 1;
  ULONGLONG CurrentState = this->UniqueStateDescriptionUpPlus[PosMidUpPlus];
  while ((PosMax > PosMin) && (CurrentState != stateDescriptionUpPlus))
    {
       if (CurrentState > stateDescriptionUpPlus)
	 {
	   PosMin = PosMidUpPlus + 1;
	 }
       else
 	{
 	  PosMax = PosMidUpPlus - 1;
	} 
       PosMidUpPlus = (PosMin + PosMax) >> 1;
       CurrentState = this->UniqueStateDescriptionUpPlus[PosMidUpPlus];
    }
  if (CurrentState != stateDescriptionUpPlus)
    PosMidUpPlus = PosMax;

  ULONGLONG* TmpStateDescriptionArray = this->UniqueStateDescriptionUpMinus[PosMidUpPlus];
  PosMin = 0;
  PosMax = this->NbrUniqueStateDescriptionUpMinus[PosMidUpPlus] - 1;
  int PosMidUpMinus = (PosMin + PosMax) >> 1;
  CurrentState = TmpStateDescriptionArray[PosMidUpMinus];
  while ((PosMax > PosMin) && (CurrentState != stateDescriptionUpMinus))
    {
       if (CurrentState > stateDescriptionUpMinus)
	 {
	   PosMin = PosMidUpMinus + 1;
	 }
       else
 	{
 	  PosMax = PosMidUpMinus - 1;
	} 
       PosMidUpMinus = (PosMin + PosMax) >> 1;
       CurrentState = TmpStateDescriptionArray[PosMidUpMinus];
    }
  if (CurrentState != stateDescriptionUpMinus)
    PosMidUpMinus = PosMax;

  TmpStateDescriptionArray = this->UniqueStateDescriptionDownPlus[PosMidUpPlus][PosMidUpMinus];
  PosMin = 0;
  PosMax = this->NbrUniqueStateDescriptionDownPlus[PosMidUpPlus][PosMidUpMinus] - 1;
  int PosMid = (PosMin + PosMax) >> 1;
  CurrentState = TmpStateDescriptionArray[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescriptionDownPlus))
    {
       if (CurrentState > stateDescriptionDownPlus)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = TmpStateDescriptionArray[PosMid];
    }
  if (CurrentState != stateDescriptionDownPlus)
    PosMid = PosMax;

  PosMin = this->FirstIndexUniqueStateDescriptionDownPlus[PosMidUpPlus][PosMidUpMinus][PosMid];
  PosMax = PosMin + this->UniqueStateDescriptionSubArraySizeDownPlus[PosMidUpPlus][PosMidUpMinus][PosMid] - 1;
  PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->StateDescriptionDownMinus[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescriptionDownMinus))
    {
       if (CurrentState > stateDescriptionDownMinus)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->StateDescriptionDownMinus[PosMid];
    }

  if (CurrentState != stateDescriptionDownMinus)
    return PosMax;
  else
    return PosMid;
}



// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSU4SpinLong::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[state], this->StateDescriptionUpMinus[state], 
		       this->StateDescriptionDownPlus[state], this->StateDescriptionDownMinus[state],
		       this->TemporaryStateUpPlus, this->TemporaryStateUpMinus,
		       this->TemporaryStateDownPlus, this->TemporaryStateDownMinus); 

  Str << " | ";
  for (int i = this->LzMax; i >=0 ; --i)
    {
      Str << "(" << this->TemporaryStateUpPlus[i] << "," << this->TemporaryStateUpMinus[i] << "," 
	  << this->TemporaryStateDownPlus[i] << "," << this->TemporaryStateDownMinus[i] << ") | ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax1 = momentum maximum value for a boson in the state with up-plus
// lzMax2 = momentum maximum value for a boson in the state with up-minus
// lzMax3 = momentum maximum value for a boson in the state with down-plus
// lzMax4 = momentum maximum value for a boson in the state with down-minus
// totalLz = momentum total value
// nbrNUpPlus = number of particles with quantum up-plus
// nbrNUpMinus = number of particles with quantum up-minus
// nbrNDownPlus = number of particles with quantum down-plus
// nbrNDownMinus = number of particles with quantum down-minus
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSU4SpinLong::GenerateStates(int nbrBosons, int lzMax1, int lzMax2, int lzMax3, int lzMax4, int totalLz, 
					      int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus, long pos)
{
  if ((nbrBosons < 0) || (totalLz < 0) || (nbrNUpPlus < 0) || (nbrNUpMinus < 0) || (nbrNDownPlus < 0) || (nbrNDownMinus < 0))
    return pos;
  if ((nbrBosons == 0) && (totalLz == 0))
    {
      this->StateDescriptionUpPlus[pos] = (ULONGLONG)0x0ul;
      this->StateDescriptionUpMinus[pos] = (ULONGLONG)0x0ul;
      this->StateDescriptionDownPlus[pos] = (ULONGLONG)0x0ul;
      this->StateDescriptionDownMinus[pos] = (ULONGLONG)0x0ul;
      return (pos + 1l);
    }
  if ((lzMax1 < 0) || (lzMax2 < 0) || (lzMax3 < 0) || (lzMax4 < 0))
    return pos;

  if (nbrBosons == 1) 
    {
      if ((nbrNDownMinus == 1) && (lzMax4 >= totalLz))
	{
	  this->StateDescriptionUpPlus[pos] = (ULONGLONG)0x0ul;
	  this->StateDescriptionUpMinus[pos] = (ULONGLONG)0x0ul;
	  this->StateDescriptionDownPlus[pos] = (ULONGLONG)0x0ul;
	  this->StateDescriptionDownMinus[pos] = (ULONGLONG)0x1ul << totalLz;
	  return (pos + 1l);
	}
      return pos;
    }

  long TmpPos;
  ULONGLONG Mask1;
  ULONGLONG Mask2;
  ULONGLONG Mask3;
  ULONGLONG Mask4;

  if (nbrNUpPlus == 0)
    {
      if (nbrNUpMinus == 0)
	{
	  if (nbrNDownPlus == 0)
	    {
	      for (int l = nbrNDownMinus; l > 0; --l)
		{
		  TmpPos = this->GenerateStates(nbrBosons - l, 0, 0, 0, lzMax4 - 1, totalLz - (lzMax4 * l), 
						0, 0, 0, nbrNDownMinus - l, pos); 
		  Mask4 = (((ULONGLONG)0x1ul << l) - (ULONGLONG)1ul) << (lzMax4 + nbrNDownMinus - l);
		  for (; pos < TmpPos; ++pos)
		    {
		      this->StateDescriptionDownMinus[pos] |= Mask4;
		    }
		}
	      pos = this->GenerateStates(nbrBosons, 0, 0, 0, lzMax4 - 1, totalLz, 0, 0, 0, nbrNDownMinus, pos);
	      return pos;
	    }
	  for (int k = nbrNDownPlus; k > 0; --k)
	    {
	      TmpPos = this->GenerateStates(nbrBosons - k, 0, 0, lzMax3 - 1, lzMax4, totalLz - (lzMax3 * k), 
					    0, 0, nbrNDownPlus - k, nbrNDownMinus, pos); 
	      Mask3 = (((ULONGLONG)0x1ul << k) - (ULONGLONG)1ul) << (lzMax3 + nbrNDownPlus - k);
	      for (; pos < TmpPos; ++pos)
		{
		  this->StateDescriptionDownPlus[pos] |= Mask3;
		}
	    }
	  pos = this->GenerateStates(nbrBosons, 0, 0, lzMax3 - 1, lzMax4, totalLz, 0, 0, nbrNDownPlus, nbrNDownMinus, pos);
	  return pos;
	}
      TmpPos = this->GenerateStates(nbrBosons - nbrNUpMinus, 0, 0, lzMax3, lzMax4, totalLz - (lzMax2 * nbrNUpMinus), 
				    0, 0, nbrNDownPlus, nbrNDownMinus, pos); 
      Mask2 = (((ULONGLONG)0x1ul << nbrNUpMinus) - (ULONGLONG)1ul) << lzMax2;
      for (; pos < TmpPos; ++pos)
	this->StateDescriptionUpMinus[pos] |= Mask2;
      for (int j = nbrNUpMinus - 1; j > 0; --j)
	{
	  TmpPos = this->GenerateStates(nbrBosons - j, 0, lzMax2 - 1, lzMax3, lzMax4, totalLz - (lzMax2 * j), 
					0, nbrNUpMinus - j, nbrNDownPlus, nbrNDownMinus, pos); 
	  Mask2 = (((ULONGLONG)0x1ul << j) - (ULONGLONG)1ul) << (lzMax2 + nbrNUpMinus - j);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescriptionUpMinus[pos] |= Mask2;
	}
      pos = this->GenerateStates(nbrBosons, 0, lzMax2 - 1, lzMax3, lzMax4, totalLz, 
				 0, nbrNUpMinus, nbrNDownPlus, nbrNDownMinus, pos);
      return pos;
    }
  
  TmpPos = this->GenerateStates(nbrBosons - nbrNUpPlus, 0, lzMax2, lzMax3, lzMax4, totalLz - (lzMax1 * nbrNUpPlus), 
				0, nbrNUpMinus, nbrNDownPlus, nbrNDownMinus, pos); 
  Mask1 = (((ULONGLONG)0x1ul << nbrNUpPlus) - (ULONGLONG)1ul) << lzMax1;
  for (; pos < TmpPos; ++pos)
    this->StateDescriptionUpPlus[pos] |= Mask1;
  for (int i = nbrNUpPlus - 1; i > 0; --i)
    {
      TmpPos = this->GenerateStates(nbrBosons - i, lzMax1 - 1, lzMax2, lzMax3, lzMax4, totalLz - (lzMax1 * i), 
				    nbrNUpPlus - i, nbrNUpMinus, nbrNDownPlus, nbrNDownMinus, pos); 
      Mask1 = (((ULONGLONG)0x1ul << i) - (ULONGLONG)1ul) << (lzMax1 + nbrNUpPlus - i);
      for (; pos < TmpPos; ++pos)
	{
	  this->StateDescriptionUpPlus[pos] |= Mask1;
	}
    }
  pos = this->GenerateStates(nbrBosons, lzMax1 - 1, lzMax2, lzMax3, lzMax4, totalLz, 
			     nbrNUpPlus, nbrNUpMinus, nbrNDownPlus, nbrNDownMinus, pos);
  return pos;
};


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereWithSU4SpinLong::GenerateLookUpTable(unsigned long memory)
{  
  long TmpUniquePartition = 1l;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescriptionUpPlus[i - 1] == this->StateDescriptionUpPlus[i]))
	{
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	++TmpUniquePartition;
    }

  this->NbrUniqueStateDescriptionUpPlus = TmpUniquePartition;
  this->UniqueStateDescriptionUpPlus = new ULONGLONG [this->NbrUniqueStateDescriptionUpPlus];
  this->UniqueStateDescriptionSubArraySizeUpPlus = new int [this->NbrUniqueStateDescriptionUpPlus];
  TmpUniquePartition = 0l;
  this->UniqueStateDescriptionUpPlus[0l] = this->StateDescriptionUpPlus[0l];
  this->UniqueStateDescriptionSubArraySizeUpPlus[0l] = 1;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescriptionUpPlus[i - 1] == this->StateDescriptionUpPlus[i]))
	{
	  ++this->UniqueStateDescriptionSubArraySizeUpPlus[TmpUniquePartition];
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	{
	  ++TmpUniquePartition;
	  this->UniqueStateDescriptionUpPlus[TmpUniquePartition] = this->StateDescriptionUpPlus[i];
	  this->UniqueStateDescriptionSubArraySizeUpPlus[TmpUniquePartition] = 1; 
	}
    }

  this->NbrUniqueStateDescriptionUpMinus = new int [this->NbrUniqueStateDescriptionUpPlus];
  TmpUniquePartition = 0;
  long TmpIndex = 0l;
  while (TmpIndex < this->LargeHilbertSpaceDimension)
    {
      long Lim = TmpIndex + this->UniqueStateDescriptionSubArraySizeUpPlus[TmpUniquePartition];
      this->NbrUniqueStateDescriptionUpMinus[TmpUniquePartition] = 1;
      ++TmpIndex;
      while (TmpIndex < Lim)
	{
	  while ((TmpIndex < Lim) && (this->StateDescriptionUpMinus[TmpIndex - 1] == this->StateDescriptionUpMinus[TmpIndex]))
	    ++TmpIndex;
	  if (TmpIndex < Lim)
	    {
	      ++this->NbrUniqueStateDescriptionUpMinus[TmpUniquePartition];
	      ++TmpIndex;
	    }
	}
      ++TmpUniquePartition;
    }
  this->UniqueStateDescriptionUpMinus = new ULONGLONG* [this->NbrUniqueStateDescriptionUpPlus];
  this->UniqueStateDescriptionSubArraySizeUpMinus = new int* [this->NbrUniqueStateDescriptionUpPlus];
  this->FirstIndexUniqueStateDescriptionUpMinus = new int* [this->NbrUniqueStateDescriptionUpPlus];
  this->NbrUniqueStateDescriptionDownPlus = new int* [this->NbrUniqueStateDescriptionUpPlus];
  this->UniqueStateDescriptionDownPlus = new ULONGLONG** [this->NbrUniqueStateDescriptionUpPlus];
  this->UniqueStateDescriptionSubArraySizeDownPlus = new int** [this->NbrUniqueStateDescriptionUpPlus];
  this->FirstIndexUniqueStateDescriptionDownPlus = new int** [this->NbrUniqueStateDescriptionUpPlus];
  for (long i = 0l; i < this->NbrUniqueStateDescriptionUpPlus; ++i)
    {
      this->UniqueStateDescriptionUpMinus[i] = new ULONGLONG [this->NbrUniqueStateDescriptionUpMinus[i]];
      this->UniqueStateDescriptionSubArraySizeUpMinus[i] = new int [this->NbrUniqueStateDescriptionUpMinus[i]];
      this->FirstIndexUniqueStateDescriptionUpMinus[i] = new int [this->NbrUniqueStateDescriptionUpMinus[i]];
      this->NbrUniqueStateDescriptionDownPlus[i] = new int [this->NbrUniqueStateDescriptionUpMinus[i]]; 
      this->UniqueStateDescriptionDownPlus[i] = new ULONGLONG* [this->NbrUniqueStateDescriptionUpMinus[i]];
      this->UniqueStateDescriptionSubArraySizeDownPlus[i] = new int* [this->NbrUniqueStateDescriptionUpMinus[i]];
      this->FirstIndexUniqueStateDescriptionDownPlus[i] = new int* [this->NbrUniqueStateDescriptionUpMinus[i]];
   }

  TmpUniquePartition = 0;
  TmpIndex = 0l;
  while (TmpIndex < this->LargeHilbertSpaceDimension)
    {
      long Lim = TmpIndex + this->UniqueStateDescriptionSubArraySizeUpPlus[TmpUniquePartition];
      int TmpUniquePartition2 = 0;
      this->UniqueStateDescriptionUpMinus[TmpUniquePartition][TmpUniquePartition2] = this->StateDescriptionUpMinus[TmpIndex];
      this->UniqueStateDescriptionSubArraySizeUpMinus[TmpUniquePartition][TmpUniquePartition2] = 1;
      this->FirstIndexUniqueStateDescriptionUpMinus[TmpUniquePartition][TmpUniquePartition2] = TmpIndex;
      ++TmpIndex;
      while (TmpIndex < Lim)
	{
	  while ((TmpIndex < Lim) && (this->StateDescriptionUpMinus[TmpIndex - 1] == this->StateDescriptionUpMinus[TmpIndex]))
	    {
	      ++this->UniqueStateDescriptionSubArraySizeUpMinus[TmpUniquePartition][TmpUniquePartition2];	      
	      ++TmpIndex;
	    }
	  if (TmpIndex < Lim)
	    {
	      ++TmpUniquePartition2;
	      this->UniqueStateDescriptionUpMinus[TmpUniquePartition][TmpUniquePartition2] = this->StateDescriptionUpMinus[TmpIndex];
	      this->UniqueStateDescriptionSubArraySizeUpMinus[TmpUniquePartition][TmpUniquePartition2] = 1;
	      this->FirstIndexUniqueStateDescriptionUpMinus[TmpUniquePartition][TmpUniquePartition2] = TmpIndex;
	      ++TmpIndex;
	    }
	}
      ++TmpUniquePartition;
    }

  for (long i = 0l; i < this->NbrUniqueStateDescriptionUpPlus; ++i)
    {
      for (long j = 0; j < this->NbrUniqueStateDescriptionUpMinus[i]; ++j)
	{
	  TmpIndex = this->FirstIndexUniqueStateDescriptionUpMinus[i][j];
	  long Lim = TmpIndex + this->UniqueStateDescriptionSubArraySizeUpMinus[i][j];
	  long TmpUnique = 1l;
	  ++TmpIndex;
	  while (TmpIndex < Lim)
	    {
	      while ((TmpIndex < Lim) && (this->StateDescriptionDownPlus[TmpIndex - 1] == this->StateDescriptionDownPlus[TmpIndex]))
		++TmpIndex;
	      if (TmpIndex < Lim)
		{
		  ++TmpIndex;
		  TmpUnique++;
		}
	    }	  
	  this->NbrUniqueStateDescriptionDownPlus[i][j] = TmpUnique;
	  this->UniqueStateDescriptionDownPlus[i][j] = new ULONGLONG [TmpUnique];
	  this->UniqueStateDescriptionSubArraySizeDownPlus[i][j] = new int [TmpUnique];
	  this->FirstIndexUniqueStateDescriptionDownPlus[i][j] = new int [TmpUnique];
	}
    }

  for (long i = 0l; i < this->NbrUniqueStateDescriptionUpPlus; ++i)
    {
      for (long j = 0; j < this->NbrUniqueStateDescriptionUpMinus[i]; ++j)
	{
	  TmpIndex = this->FirstIndexUniqueStateDescriptionUpMinus[i][j];
	  long Lim = TmpIndex + this->UniqueStateDescriptionSubArraySizeUpMinus[i][j];
	  long Tmp = 0l;
	  this->UniqueStateDescriptionDownPlus[i][j][Tmp] = this->StateDescriptionDownPlus[TmpIndex];	
	  this->UniqueStateDescriptionSubArraySizeDownPlus[i][j][Tmp] = 1;
	  this->FirstIndexUniqueStateDescriptionDownPlus[i][j][Tmp] = TmpIndex;	  
	  long TmpUnique = 1l;
	  ++TmpIndex;
	  while (TmpIndex < Lim)
	    {
	      while ((TmpIndex < Lim) && (this->StateDescriptionDownPlus[TmpIndex - 1] == this->StateDescriptionDownPlus[TmpIndex]))
		{
		  ++TmpIndex;
		  ++this->UniqueStateDescriptionSubArraySizeDownPlus[i][j][Tmp];
		}
	      if (TmpIndex < Lim)
		{
		  ++Tmp;
		  this->UniqueStateDescriptionDownPlus[i][j][Tmp] = this->StateDescriptionDownPlus[TmpIndex];	
		  this->UniqueStateDescriptionSubArraySizeDownPlus[i][j][Tmp] = 1;
		  this->FirstIndexUniqueStateDescriptionDownPlus[i][j][Tmp] = TmpIndex;	  
		  ++TmpIndex;
		}
	    }	  
	}
    }
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrNUpPlus = number of particles with quantum number up-plus
// nbrNUpMinus = number of particles with quantum number up-minus
// nbrNDownPlus = number of particles with quantum number down-plus
// nbrNDownMinus = number of particles with quantum number down-minus
// return value = Hilbert space dimension

long BosonOnSphereWithSU4SpinLong::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, 
								    int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus)
{
  if ((nbrBosons < 0) || (totalLz < 0) || (nbrNUpPlus < 0) || (nbrNUpMinus < 0) || (nbrNDownPlus < 0) || (nbrNDownMinus < 0))
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
  for (int i = nbrNUpPlus; i >= 0; --i)
    for (int j = nbrNUpMinus; j >= 0; --j)
      for (int k = nbrNDownPlus; k >= 0; --k)
	for (int l = nbrNDownMinus; l >= 0; --l)
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - (i + j + k + l), lzMax - 1, totalLz - (lzMax * (i + j + k + l)), 
							    nbrNUpPlus - i, nbrNUpMinus - j, nbrNDownPlus - k, nbrNDownMinus - l);
  return  Tmp;
}

// apply a^+_m_up a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AdupAup (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->TemporaryStateUpPlus);
  if (this->TemporaryStateUpPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpPlus[n];
  --this->TemporaryStateUpPlus[n];
  ++this->TemporaryStateUpPlus[m];
  coefficient *= (double) this->TemporaryStateUpPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStateUpPlus), this->StateDescriptionUpMinus[index], 
			      this->StateDescriptionDownPlus[index], this->StateDescriptionDownMinus[index]);  
}

// apply a^+_m_up a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AdupAum (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->TemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->TemporaryStateUpMinus);
  if (this->TemporaryStateUpMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpMinus[n];
  --this->TemporaryStateUpMinus[n];
  ++this->TemporaryStateUpPlus[m];
  coefficient *= (double) this->TemporaryStateUpPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStateUpPlus), this->BosonToFermion(this->TemporaryStateUpMinus), 
			      this->StateDescriptionDownPlus[index], this->StateDescriptionDownMinus[index]);  
}

// apply a^+_m_up a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AdupAdp (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->TemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->TemporaryStateDownPlus);
  if (this->TemporaryStateDownPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownPlus[n];
  --this->TemporaryStateDownPlus[n];
  ++this->TemporaryStateUpPlus[m];
  coefficient *= (double) this->TemporaryStateUpPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStateUpPlus), this->StateDescriptionUpMinus[index],
			      this->BosonToFermion(this->TemporaryStateDownPlus), this->StateDescriptionDownMinus[index]);  
}

// apply a^+_m_up a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AdupAdm (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->TemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->TemporaryStateDownMinus);
  if (this->TemporaryStateDownMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownMinus[n];
  --this->TemporaryStateDownMinus[n];
  ++this->TemporaryStateUpPlus[m];
  coefficient *= (double) this->TemporaryStateUpPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStateUpPlus), this->StateDescriptionUpMinus[index],
			      this->StateDescriptionDownPlus[index], this->BosonToFermion(this->TemporaryStateDownMinus));  
}

// apply a^+_m_um a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AdumAup (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->TemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->TemporaryStateUpMinus);
  if (this->TemporaryStateUpPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpPlus[n];
  --this->TemporaryStateUpPlus[n];
  ++this->TemporaryStateUpMinus[m];
  coefficient *= (double) this->TemporaryStateUpMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStateUpPlus), this->BosonToFermion(this->TemporaryStateUpMinus), 
			      this->StateDescriptionDownPlus[index], this->StateDescriptionDownMinus[index]);  
}

// apply a^+_m_um a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AdumAum (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->TemporaryStateUpMinus);
  if (this->TemporaryStateUpMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpMinus[n];
  --this->TemporaryStateUpMinus[n];
  ++this->TemporaryStateUpMinus[m];
  coefficient *= (double) this->TemporaryStateUpMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionUpPlus[index], this->BosonToFermion(this->TemporaryStateUpMinus), 
			      this->StateDescriptionDownPlus[index], this->StateDescriptionDownMinus[index]);  
}

// apply a^+_m_um a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AdumAdp (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->TemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->TemporaryStateUpMinus);
  if (this->TemporaryStateDownPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownPlus[n];
  --this->TemporaryStateDownPlus[n];
  ++this->TemporaryStateUpMinus[m];
  coefficient *= (double) this->TemporaryStateUpMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionUpPlus[index], this->BosonToFermion(this->TemporaryStateUpMinus), 
			      this->BosonToFermion(this->TemporaryStateDownPlus), this->StateDescriptionDownMinus[index]);  
}

// apply a^+_m_um a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AdumAdm (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->TemporaryStateDownMinus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->TemporaryStateUpMinus);
  if (this->TemporaryStateDownMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownMinus[n];
  --this->TemporaryStateDownMinus[n];
  ++this->TemporaryStateUpMinus[m];
  coefficient *= (double) this->TemporaryStateUpMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionUpPlus[index], this->BosonToFermion(this->TemporaryStateUpMinus), 
			      this->StateDescriptionDownPlus[index], this->BosonToFermion(this->TemporaryStateDownMinus));
}

// apply a^+_m_dp a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AddpAup (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->TemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->TemporaryStateDownPlus);
  if (this->TemporaryStateUpPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpPlus[n];
  --this->TemporaryStateUpPlus[n];
  ++this->TemporaryStateDownPlus[m];
  coefficient *= (double) this->TemporaryStateDownPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStateUpPlus), this->StateDescriptionUpMinus[index],
			      this->BosonToFermion(this->TemporaryStateDownPlus), this->StateDescriptionDownMinus[index]);  
}

// apply a^+_m_dp a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AddpAum (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->TemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->TemporaryStateUpMinus);
  if (this->TemporaryStateUpMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpMinus[n];
  --this->TemporaryStateUpMinus[n];
  ++this->TemporaryStateDownPlus[m];
  coefficient *= (double) this->TemporaryStateDownPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionUpPlus[index], this->BosonToFermion(this->TemporaryStateUpMinus), 
			      this->BosonToFermion(this->TemporaryStateDownPlus), this->StateDescriptionDownMinus[index]);  
}

// apply a^+_m_dp a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AddpAdp (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->TemporaryStateDownPlus);
  if (this->TemporaryStateDownPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownPlus[n];
  --this->TemporaryStateDownPlus[n];
  ++this->TemporaryStateDownPlus[m];
  coefficient *= (double) this->TemporaryStateDownPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionUpPlus[index], this->StateDescriptionUpMinus[index], 
			      this->BosonToFermion(this->TemporaryStateDownPlus), this->StateDescriptionDownMinus[index]);  
}

// apply a^+_m_dp a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AddpAdm (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->TemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->TemporaryStateDownMinus);
  if (this->TemporaryStateDownMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownMinus[n];
  --this->TemporaryStateDownMinus[n];
  ++this->TemporaryStateDownPlus[m];
  coefficient *= (double) this->TemporaryStateDownPlus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionUpPlus[index], this->StateDescriptionUpMinus[index], 
			      this->BosonToFermion(this->TemporaryStateDownPlus), this->BosonToFermion(this->TemporaryStateDownMinus));  
}

// apply a^+_m_dm a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AddmAup (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->TemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->TemporaryStateDownMinus);
  if (this->TemporaryStateUpPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpPlus[n];
  --this->TemporaryStateUpPlus[n];
  ++this->TemporaryStateDownMinus[m];
  coefficient *= (double) this->TemporaryStateDownMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryStateUpPlus), this->StateDescriptionUpMinus[index],
			      this->StateDescriptionDownPlus[index], this->BosonToFermion(this->TemporaryStateDownMinus));  
}

// apply a^+_m_dm a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AddmAum (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->TemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->TemporaryStateDownMinus);
  if (this->TemporaryStateUpMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUpMinus[n];
  --this->TemporaryStateUpMinus[n];
  ++this->TemporaryStateDownMinus[m];
  coefficient *= (double) this->TemporaryStateDownMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionUpPlus[index], this->BosonToFermion(this->TemporaryStateUpMinus), 
			      this->StateDescriptionDownPlus[index], this->BosonToFermion(this->TemporaryStateDownMinus));  
}

// apply a^+_m_dm a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AddmAdp (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->TemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->TemporaryStateDownMinus);
  if (this->TemporaryStateDownPlus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownPlus[n];
  --this->TemporaryStateDownPlus[n];
  ++this->TemporaryStateDownMinus[m];
  coefficient *= (double) this->TemporaryStateDownMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionUpPlus[index], this->StateDescriptionUpMinus[index], 
			      this->BosonToFermion(this->TemporaryStateDownPlus), this->BosonToFermion(this->TemporaryStateDownMinus));  
}

// apply a^+_m_dm a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU4SpinLong::AddmAdm (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->TemporaryStateDownMinus);
  if (this->TemporaryStateDownMinus[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDownMinus[n];
  --this->TemporaryStateDownMinus[n];
  ++this->TemporaryStateDownMinus[m];
  coefficient *= (double) this->TemporaryStateDownMinus[m];
  coefficient = sqrt(coefficient);  
  return this->FindStateIndex(this->StateDescriptionUpPlus[index], this->StateDescriptionUpMinus[index], 
			      this->StateDescriptionDownPlus[index], this->BosonToFermion(this->TemporaryStateDownMinus));  
}

// apply a_n1_up a_n2_up operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinLong::AupAup (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->ProdATemporaryStateUpPlus);
  if ((this->ProdATemporaryStateUpPlus[n1] == 0) || (this->ProdATemporaryStateUpPlus[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateUpPlus[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->ProdATemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->ProdATemporaryStateDownMinus);
  double Coefficient = this->ProdATemporaryStateUpPlus[n2];
  --this->ProdATemporaryStateUpPlus[n2];
  Coefficient *= this->ProdATemporaryStateUpPlus[n1];
  --this->ProdATemporaryStateUpPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_up a_n2_um operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinLong::AupAum (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->ProdATemporaryStateUpMinus);
  if ((this->ProdATemporaryStateUpPlus[n1] == 0) || (this->ProdATemporaryStateUpMinus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->ProdATemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->ProdATemporaryStateDownMinus);
  double Coefficient = this->ProdATemporaryStateUpMinus[n2];
  --this->ProdATemporaryStateUpMinus[n2];
  Coefficient *= this->ProdATemporaryStateUpPlus[n1];
  --this->ProdATemporaryStateUpPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_up a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinLong::AupAdp (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->ProdATemporaryStateDownPlus);
  if ((this->ProdATemporaryStateUpPlus[n1] == 0) || (this->ProdATemporaryStateDownPlus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->ProdATemporaryStateDownMinus);
  double Coefficient = this->ProdATemporaryStateDownPlus[n2];
  --this->ProdATemporaryStateDownPlus[n2];
  Coefficient *= this->ProdATemporaryStateUpPlus[n1];
  --this->ProdATemporaryStateUpPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_up a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinLong::AupAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->ProdATemporaryStateDownMinus);
  if ((this->ProdATemporaryStateUpPlus[n1] == 0) || (this->ProdATemporaryStateDownMinus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->ProdATemporaryStateDownPlus);
  double Coefficient = this->ProdATemporaryStateDownMinus[n2];
  --this->ProdATemporaryStateDownMinus[n2];
  Coefficient *= this->ProdATemporaryStateUpPlus[n1];
  --this->ProdATemporaryStateUpPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_um a_n2_um operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinLong::AumAum (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->ProdATemporaryStateUpMinus);
  if ((this->ProdATemporaryStateUpMinus[n1] == 0) || (this->ProdATemporaryStateUpMinus[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateUpMinus[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->ProdATemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->ProdATemporaryStateDownMinus);
  double Coefficient = this->ProdATemporaryStateUpMinus[n2];
  --this->ProdATemporaryStateUpMinus[n2];
  Coefficient *= this->ProdATemporaryStateUpMinus[n1];
  --this->ProdATemporaryStateUpMinus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_um a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinLong::AumAdp (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->ProdATemporaryStateDownPlus);
  if ((this->ProdATemporaryStateUpMinus[n1] == 0) || (this->ProdATemporaryStateDownPlus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->ProdATemporaryStateDownMinus);
  double Coefficient = this->ProdATemporaryStateDownPlus[n2];
  --this->ProdATemporaryStateDownPlus[n2];
  Coefficient *= this->ProdATemporaryStateUpMinus[n1];
  --this->ProdATemporaryStateUpMinus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_um a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinLong::AumAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->ProdATemporaryStateDownMinus);
  if ((this->ProdATemporaryStateUpMinus[n1] == 0) || (this->ProdATemporaryStateDownMinus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->ProdATemporaryStateDownPlus);
  double Coefficient = this->ProdATemporaryStateDownMinus[n2];
  --this->ProdATemporaryStateDownMinus[n2];
  Coefficient *= this->ProdATemporaryStateUpMinus[n1];
  --this->ProdATemporaryStateUpMinus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_dp a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinLong::AdpAdp (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->ProdATemporaryStateDownPlus);
  if ((this->ProdATemporaryStateDownPlus[n1] == 0) || (this->ProdATemporaryStateDownPlus[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateDownPlus[n1] == 1)))    
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->ProdATemporaryStateDownMinus);
  double Coefficient = this->ProdATemporaryStateDownPlus[n2];
  --this->ProdATemporaryStateDownPlus[n2];
  Coefficient *= this->ProdATemporaryStateDownPlus[n1];
  --this->ProdATemporaryStateDownPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_dp a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinLong::AdpAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->ProdATemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->ProdATemporaryStateDownMinus);
  if ((this->ProdATemporaryStateDownPlus[n1] == 0) || (this->ProdATemporaryStateDownMinus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->ProdATemporaryStateUpMinus);
  double Coefficient = this->ProdATemporaryStateDownMinus[n2];
  --this->ProdATemporaryStateDownMinus[n2];
  Coefficient *= this->ProdATemporaryStateDownPlus[n1];
  --this->ProdATemporaryStateDownPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_dm a_n2_dm operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU4SpinLong::AdmAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusLzMax, this->ProdATemporaryStateDownMinus);
  if ((this->ProdATemporaryStateDownMinus[n1] == 0) || (this->ProdATemporaryStateDownMinus[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateDownMinus[n1] == 1)))    
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusLzMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusLzMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusLzMax, this->ProdATemporaryStateDownPlus);
  double Coefficient = this->ProdATemporaryStateDownMinus[n2];
  --this->ProdATemporaryStateDownMinus[n2];
  Coefficient *= this->ProdATemporaryStateDownMinus[n1];
  --this->ProdATemporaryStateDownMinus[n1];
  return sqrt(Coefficient);
}

// apply a^+_m1_up a^+_m2_up operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinLong::AdupAdup (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateUpPlus, coefficient);
}

// apply a^+_m1_up a^+_m2_um operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinLong::AdupAdum (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, coefficient);
}

// apply a^+_m1_up a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinLong::AdupAddp (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateDownPlus, coefficient);
}

// apply a^+_m1_up a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinLong::AdupAddm (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateDownMinus, coefficient);
}

// apply a^+_m1_um a^+_m2_um operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinLong::AdumAdum (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpMinus, this->TemporaryStateUpMinus, coefficient);
}

// apply a^+_m1_um a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinLong::AdumAddp (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpMinus, this->TemporaryStateDownPlus, coefficient);
}

// apply a^+_m1_um a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinLong::AdumAddm (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpMinus, this->TemporaryStateDownMinus, coefficient);
}

// apply a^+_m1_dp a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinLong::AddpAddp (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDownPlus, this->TemporaryStateDownPlus, coefficient);
}

// apply a^+_m1_dp a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinLong::AddpAddm (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDownPlus, this->TemporaryStateDownMinus, coefficient);
}

// apply a^+_m1_dm a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU4SpinLong::AddmAddm (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDownMinus, this->TemporaryStateDownMinus, coefficient);
}

// convert a state from one SU(4) basis to another, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void BosonOnSphereWithSU4SpinLong::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, 
							 long firstComponent, long nbrComponents)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU4Indices = new int [this->NbrBosons];
  int* TmpSU4Indices2 = new int [this->NbrBosons];
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
      this->FermionToBoson(this->StateDescriptionUpPlus[i], this->StateDescriptionUpMinus[i], 
			   this->StateDescriptionDownPlus[i], this->StateDescriptionDownMinus[i],
			   this->TemporaryStateUpPlus, this->TemporaryStateUpMinus,
			   this->TemporaryStateDownPlus, this->TemporaryStateDownMinus); 
      double OccupationCoefficient = 0.0;
      int TmpIndex = 0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  for (unsigned l = 0; l < this->TemporaryStateDownMinus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 3;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateDownPlus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 2;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateUpMinus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 1;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateUpPlus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 0;
	      ++TmpIndex;	      
	    }
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUpPlus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUpMinus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDownPlus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDownMinus[j]];
	}
      this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpMomentumIndices, TmpSU4Indices, TmpSU4Indices2, oneBodyBasis, OccupationCoefficient, OccupationCoefficientArray);
    }
  delete[] OccupationCoefficientArray;
  delete[] TmpMomentumIndices;
  delete[] TmpSU4Indices;
  delete[] TmpSU4Indices2;
}

// compute the transformation matrix from one SU(4) basis to another, transforming the one body basis in each momentum sector
//
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// return value = transformation matrix

ComplexMatrix BosonOnSphereWithSU4SpinLong::TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU4Indices = new int [this->NbrBosons];
  int* TmpSU4Indices2 = new int [this->NbrBosons];
  ComplexMatrix TmpMatrix(this->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  double* OccupationCoefficientArray = new double [this->NbrBosons + 1];
  OccupationCoefficientArray[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    OccupationCoefficientArray[i] = OccupationCoefficientArray[i - 1] + 0.5 * log((double) i);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUpPlus[i], this->StateDescriptionUpMinus[i], 
			   this->StateDescriptionDownPlus[i], this->StateDescriptionDownMinus[i],
			   this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, 
			   this->TemporaryStateDownPlus, this->TemporaryStateDownMinus); 
      int TmpIndex = 0;
      double OccupationCoefficient = 0.0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  for (unsigned l = 0; l < this->TemporaryStateDownMinus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 3;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateDownPlus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 2;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateUpMinus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 1;
	      ++TmpIndex;	      
	    }
	  for (unsigned l = 0; l < this->TemporaryStateUpPlus[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU4Indices[TmpIndex] = 0;
	      ++TmpIndex;	      
	    }
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUpPlus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUpMinus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDownPlus[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDownMinus[j]];
	}
      this->TransformOneBodyBasisRecursive(TmpMatrix[i], 1.0, 0, TmpMomentumIndices, TmpSU4Indices, TmpSU4Indices2, oneBodyBasis,
					   OccupationCoefficient, OccupationCoefficientArray);
    }
  delete[] TmpMomentumIndices;
  delete[] TmpSU4Indices;
  delete[] TmpSU4Indices2;
  delete[] OccupationCoefficientArray;
  return TmpMatrix;
}

// recursive part of the convertion from a state from one SU(4) basis to another, transforming the one body basis in each momentum sector
//
// targetState = vector where the transformed state has to be stored
// coefficient = current coefficient to assign
// position = current particle consider in the n-body state
// momentumIndices = array that gives the momentum partition of the initial n-body state
// initialSU4Indices = array that gives the spin dressing the initial n-body state
// currentSU4Indices = array that gives the spin dressing the current transformed n-body state
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
// occupationCoefficientArray = array that provides 1/2 ln (N!)

void BosonOnSphereWithSU4SpinLong::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
							      int position, int* momentumIndices, int* initialSU4Indices, int* currentSU4Indices, ComplexMatrix* oneBodyBasis,
							      double occupationCoefficient, double* occupationCoefficientArray) 
{
  if (position == this->NbrBosons)
    {
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->TemporaryStateUpPlus[i] = 0ul;
	  this->TemporaryStateUpMinus[i] = 0ul; 
	  this->TemporaryStateDownPlus[i] = 0ul;
	  this->TemporaryStateDownMinus[i] = 0ul;
	}
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  switch (currentSU4Indices[i])
	    {
	    case 0:
	      this->TemporaryStateUpPlus[momentumIndices[i]]++;
	      break;
	    case 1:
	      this->TemporaryStateUpMinus[momentumIndices[i]]++;
	      break;
	    case 2:
	      this->TemporaryStateDownPlus[momentumIndices[i]]++;
	      break;
	    case 3:
	      this->TemporaryStateDownMinus[momentumIndices[i]]++;
	      break;
	    }
	}
      int Index = this->FindStateIndex(this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, 
				       this->TemporaryStateDownPlus, this->TemporaryStateDownMinus);
      if (Index < this->HilbertSpaceDimension)
	{
	  
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateUpPlus[i]];
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateUpMinus[i]];
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateDownPlus[i]];
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateDownMinus[i]];
	    }
	  targetState[Index] += coefficient * exp (occupationCoefficient);
	}
      return;      
    }
  else
    {
      currentSU4Indices[position] = 0;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][3 - initialSU4Indices[position]][3]), position + 1, momentumIndices, initialSU4Indices, currentSU4Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
      currentSU4Indices[position] = 1;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][3 - initialSU4Indices[position]][2]), position + 1, momentumIndices, initialSU4Indices, currentSU4Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
      currentSU4Indices[position] = 2;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][3 - initialSU4Indices[position]][1]), position + 1, momentumIndices, initialSU4Indices, currentSU4Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
      currentSU4Indices[position] = 3;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][3 - initialSU4Indices[position]][0]), position + 1, momentumIndices, initialSU4Indices, currentSU4Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
    }
}

// compute the projection matrix from the SU(4) Hilbert space to an U(1) Hilbert space
// 
// targetSpace = pointer to the U(1) Hilbert space
// type = type of particles that has to be kept (0 for type up-plus, 1 for type up-minus, 2 for type down-plus, 3 for type down-minus)
// return value = projection matrix

ComplexMatrix BosonOnSphereWithSU4SpinLong::TransformationMatrixSU4ToU1(BosonOnSphereLong* targetSpace, int type)
{
  ComplexMatrix TmpMatrix (targetSpace->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  ULONGLONG* TmpStateDescription=0;
  ULONGLONG* TmpStateDescriptionOther1=0;
  ULONGLONG* TmpStateDescriptionOther2=0;
  ULONGLONG* TmpStateDescriptionOther3=0;
  switch (type)
    {
    case 0:
      {
	TmpStateDescription = this->StateDescriptionUpPlus;
	TmpStateDescriptionOther1 = this->StateDescriptionUpMinus;
	TmpStateDescriptionOther2 = this->StateDescriptionDownPlus;
	TmpStateDescriptionOther3 = this->StateDescriptionDownMinus;
      }
      break;
    case 1:
      {
	TmpStateDescription = this->StateDescriptionUpMinus;
	TmpStateDescriptionOther1 = this->StateDescriptionUpPlus;
	TmpStateDescriptionOther2 = this->StateDescriptionDownPlus;
	TmpStateDescriptionOther3 = this->StateDescriptionDownMinus;
      }
      break;
    case 2:
      {
	TmpStateDescription = this->StateDescriptionDownPlus;
	TmpStateDescriptionOther1 = this->StateDescriptionUpMinus;
	TmpStateDescriptionOther2 = this->StateDescriptionUpPlus;
	TmpStateDescriptionOther3 = this->StateDescriptionDownMinus;
      }
      break;
    case 3:
      {
	TmpStateDescription = this->StateDescriptionDownMinus;
	TmpStateDescriptionOther1 = this->StateDescriptionUpMinus;
	TmpStateDescriptionOther2 = this->StateDescriptionUpPlus;
	TmpStateDescriptionOther3 = this->StateDescriptionDownPlus;
      }
      break;
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if ((TmpStateDescriptionOther1[i] == (ULONGLONG)0x0ul) && (TmpStateDescriptionOther2[i] == (ULONGLONG)0x0ul) && 
	  (TmpStateDescriptionOther3[i] == (ULONGLONG)0x0ul))
	{
	  ULONGLONG TmpState = TmpStateDescription[i];
	  int TmpLzMax = this->FermionicLzMax;
	  while ((TmpState >> TmpLzMax) == (ULONGLONG)0x0ul)
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

// // // compute the projection matrix from the SU(4) Hilbert space to an SU(2) Hilbert space
// // // 
// // // targetSpace = pointer to the SU(2) Hilbert space
// // // spinUp = index of the component that has to be consider as a spin up
// // // spinDown = index of the component that has to be consider as a spin down
// // // return value = projection matrix
// // 
// // ComplexMatrix BosonOnSphereWithSU4SpinLong::TransformationMatrixSU4ToSU2(ParticleOnSphereWithSpin* targetSpace, int spinUp, int spinDown)
// // {
// //   BosonOnSphereWithSU2Spin* TmpTargetSpace = (BosonOnSphereWithSU2Spin*) targetSpace;
// //   ComplexMatrix TmpMatrix (TmpTargetSpace->HilbertSpaceDimension, this->HilbertSpaceDimension, true) ;
// // 
// //   ULONGLONG* SpinUpComponent = 0;
// //   ULONGLONG* SpinDownComponent = 0;
// //   ULONGLONG* UnwantedComponent1 = 0;
// //   ULONGLONG* UnwantedComponent2 = 0;
// //   
// //   if (spinUp == 0)
// //     {
// //       SpinUpComponent = this->StateDescriptionUpPlus;
// //       switch (spinDown)
// // 	{
// // 	case 1:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionUpMinus;
// // 	    UnwantedComponent1 = this->StateDescriptionDownPlus;
// // 	    UnwantedComponent2 = this->StateDescriptionDownMinus;
// // 	  }
// // 	  break;
// // 	case 2:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionDownPlus;
// // 	    UnwantedComponent1 = this->StateDescriptionUpMinus;
// // 	    UnwantedComponent2 = this->StateDescriptionDownMinus;
// // 	  }
// // 	  break;
// // 	case 3:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionDownMinus;
// // 	    UnwantedComponent1 = this->StateDescriptionUpMinus;
// // 	    UnwantedComponent2 = this->StateDescriptionDownPlus;
// // 	  }
// // 	  break;
// // 	}
// //     }
// //   if (spinUp == 1)
// //     {
// //       SpinUpComponent = this->StateDescriptionUpMinus;
// //       switch (spinDown)
// // 	{
// // 	case 0:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionUpPlus;
// // 	    UnwantedComponent1 = this->StateDescriptionDownPlus;
// // 	    UnwantedComponent2 = this->StateDescriptionDownMinus;
// // 	  }
// // 	  break;
// // 	case 2:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionDownPlus;
// // 	    UnwantedComponent1 = this->StateDescriptionUpPlus;
// // 	    UnwantedComponent2 = this->StateDescriptionDownMinus;
// // 	  }
// // 	  break;
// // 	case 3:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionDownMinus;
// // 	    UnwantedComponent1 = this->StateDescriptionUpPlus;
// // 	    UnwantedComponent2 = this->StateDescriptionDownPlus;
// // 	  }
// // 	  break;
// // 	}
// //     }
// // 
// //   if (spinUp == 2)
// //     {
// //       SpinUpComponent = this->StateDescriptionDownPlus;
// //       switch (spinDown)
// // 	{
// // 	case 0:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionUpPlus;
// // 	    UnwantedComponent1 = this->StateDescriptionUpMinus;
// // 	    UnwantedComponent2 = this->StateDescriptionDownMinus;
// // 	  }
// // 	  break;
// // 	case 1:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionUpMinus;
// // 	    UnwantedComponent1 = this->StateDescriptionUpPlus;
// // 	    UnwantedComponent2 = this->StateDescriptionDownMinus;
// // 	  }
// // 	  break;
// // 	case 3:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionDownMinus;
// // 	    UnwantedComponent1 = this->StateDescriptionUpPlus;
// // 	    UnwantedComponent2 = this->StateDescriptionUpMinus;
// // 	  }
// // 	  break;
// // 	}
// //     }
// // 
// //   if (spinUp == 3)
// //     {
// //       SpinUpComponent = this->StateDescriptionDownMinus;
// //       switch (spinDown)
// // 	{
// // 	case 0:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionUpPlus;
// // 	    UnwantedComponent1 = this->StateDescriptionUpMinus;
// // 	    UnwantedComponent2 = this->StateDescriptionDownPlus;
// // 	  }
// // 	  break;
// // 	case 1:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionUpMinus;
// // 	    UnwantedComponent1 = this->StateDescriptionUpPlus;
// // 	    UnwantedComponent2 = this->StateDescriptionDownPlus;
// // 	  }
// // 	  break;
// // 	case 2:
// // 	  {
// // 	    SpinDownComponent = this->StateDescriptionDownPlus;
// // 	    UnwantedComponent1 = this->StateDescriptionUpPlus;
// // 	    UnwantedComponent2 = this->StateDescriptionUpMinus;
// // 	  }
// // 	  break;
// // 	}
// //     }
// // 
// //   for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// //     {
// //       if ((UnwantedComponent1[i] == (ULONGLONG)0x0ul) && (UnwantedComponent2[i] == (ULONGLONG)0x0ul))
// // 	{
// // 	  int Index = TmpTargetSpace->FindStateIndex(SpinUpComponent[i], SpinDownComponent[i]);
// // 	  if (Index < TmpTargetSpace->HilbertSpaceDimension)
// // 	    {
// // 	      TmpMatrix[i][Index] = 1.0;
// // 	    }
// // 	}
// //     }
// //   return TmpMatrix;
// // }
