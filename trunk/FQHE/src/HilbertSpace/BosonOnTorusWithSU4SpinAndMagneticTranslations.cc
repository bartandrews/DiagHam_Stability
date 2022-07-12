////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//               class of boson with SU4 spin on a torus taking               //
//                    into account magnetic translations                      //
//                                                                            //
//                        last modification : 21/06/2012                      //
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
#include "HilbertSpace/BosonOnTorusWithSU4SpinAndMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "Matrix/Matrix.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ArrayTools.h"

#include <cmath>
#include <cstdlib>


using std::cout;
using std::endl;
using std::dec;
using std::hex;
using std::abs;


// basic constructor
// 
// nbrBosons= number of bosons
// totalSpin = twice the total spin value
// totalIsospin = twice the total isospin value
// totalEntanglement = twice the total entanglement value
// maxMomentum = momentum maximum value for a boson
// kxMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// kyMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)

BosonOnTorusWithSU4SpinAndMagneticTranslations::BosonOnTorusWithSU4SpinAndMagneticTranslations (int nbrBosons, int totalSpin, int totalIsospin, int totalEntanglement, int maxMomentum, 
												int kxMomentum, int kyMomentum)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalSpin = totalSpin;
  this->TotalIsospin = totalIsospin;
  this->TotalEntanglement = totalEntanglement;
  this->NbrBosonsUpPlus = this->NbrBosons + this->TotalSpin + this->TotalIsospin + this->TotalEntanglement;
  this->NbrBosonsUpMinus = this->NbrBosons + this->TotalSpin - this->TotalIsospin - this->TotalEntanglement;
  this->NbrBosonsDownPlus = this->NbrBosons - this->TotalSpin + this->TotalIsospin - this->TotalEntanglement;
  this->NbrBosonsDownMinus = this->NbrBosons - this->TotalSpin - this->TotalIsospin + this->TotalEntanglement;
  this->NbrBosonsUpPlus >>= 2;
  this->NbrBosonsUpMinus >>= 2;
  this->NbrBosonsDownPlus  >>= 2;
  this->NbrBosonsDownMinus  >>= 2;
  
//  cout << this->NbrBosonsUpPlus << " " << this->NbrBosonsUpMinus << " " << this->NbrBosonsDownPlus << " : " << this->TotalTz << " " << this->TotalY << endl;
  this->MaxMomentum = maxMomentum;  
  this->KyMax = this->MaxMomentum - 1;
  this->MomentumModulo = FindGCD(this->NbrBosons, this->MaxMomentum);
  this->KxMomentum = kxMomentum % this->MomentumModulo;
  this->KyMomentum = kyMomentum % this->MaxMomentum;

  this->TemporaryStateUpPlus = new unsigned long[this->MaxMomentum];
  this->TemporaryStateUpMinus = new unsigned long[this->MaxMomentum];
  this->TemporaryStateDownPlus = new unsigned long[this->MaxMomentum];
  this->TemporaryStateDownMinus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateUpPlus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->MaxMomentum];
  this->NUpPlusKyMax = this->KyMax + this->NbrBosons;
  this->NUpMinusKyMax = this->KyMax + this->NbrBosons;
  this->NDownPlusKyMax = this->KyMax + this->NbrBosons;
  this->NDownMinusKyMax = this->KyMax + this->NbrBosons;
  this->FermionicKyMax = this->NUpPlusKyMax;

  this->StateShift = (this->MaxMomentum / this->MomentumModulo);
  this->LastMomentumMaskUpPlus = 0x1ul << (this->MaxMomentum + this->NbrBosonsUpPlus - 1);
  this->LastMomentumMaskUpMinus = 0x1ul << (this->MaxMomentum + this->NbrBosonsUpMinus - 1);
  this->LastMomentumMaskDownPlus = 0x1ul << (this->MaxMomentum + this->NbrBosonsDownPlus - 1);
  this->LastMomentumMaskDownMinus = 0x1ul << (this->MaxMomentum + this->NbrBosonsDownMinus - 1);

  this->MomentumIncrement = (this->NbrBosons * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((unsigned long) 1);
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= ((unsigned long) 1);
    }

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->KyMax, 0, this->NbrBosonsUpPlus, this->NbrBosonsUpMinus, 
									 this->NbrBosonsDownPlus, this->NbrBosonsDownMinus);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUpPlus = new unsigned long[this->LargeHilbertSpaceDimension];
      this->StateDescriptionUpMinus = new unsigned long[this->LargeHilbertSpaceDimension];
      this->StateDescriptionDownPlus = new unsigned long[this->LargeHilbertSpaceDimension];
      this->StateDescriptionDownMinus = new unsigned long[this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->RawGenerateStates(this->NbrBosons, this->KyMax, this->KyMax, this->KyMax, this->KyMax, 0, 
								   this->NbrBosonsUpPlus, this->NbrBosonsUpMinus, this->NbrBosonsDownPlus, this->NbrBosonsDownMinus, 0l);     
      this->LargeHilbertSpaceDimension = this->GenerateStates();
      cout  << "Dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  if (this->LargeHilbertSpaceDimension >= (1l << 30))
	    this->HilbertSpaceDimension = 0;
	  else
	    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
	  this->GenerateLookUpTable(1000000);
	  this->Flag.Initialize();
#ifdef __DEBUG__
	  long UsedMemory = 0l;
	  UsedMemory += this->LargeHilbertSpaceDimension * (3 * sizeof(unsigned long) + sizeof(int));
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
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnTorusWithSU4SpinAndMagneticTranslations::BosonOnTorusWithSU4SpinAndMagneticTranslations(const BosonOnTorusWithSU4SpinAndMagneticTranslations& bosons)
{
  this->MomentumModulo = bosons.MomentumModulo;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;

  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalSpin = bosons.TotalSpin;
  this->TotalIsospin = bosons.TotalIsospin;
  this->TotalEntanglement = bosons.TotalEntanglement;
  this->NbrBosonsUpPlus = bosons.NbrBosonsUpPlus;
  this->NbrBosonsUpMinus = bosons.NbrBosonsUpMinus;
  this->NbrBosonsDownPlus = bosons.NbrBosonsDownPlus;
  this->MaxMomentum = bosons.MaxMomentum;  
  this->KyMax = bosons.KyMax;
  this->NUpPlusKyMax = bosons.NUpPlusKyMax;
  this->NUpMinusKyMax = bosons.NUpMinusKyMax;
  this->NDownPlusKyMax = bosons.NDownPlusKyMax;
  this->NDownMinusKyMax = bosons.NDownMinusKyMax;
  this->FermionicKyMax = bosons.FermionicKyMax;

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;

  this->TemporaryStateUpPlus = new unsigned long[this->MaxMomentum];
  this->TemporaryStateUpMinus = new unsigned long[this->MaxMomentum];
  this->TemporaryStateDownPlus = new unsigned long[this->MaxMomentum];
  this->TemporaryStateDownMinus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateUpPlus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->MaxMomentum];

  this->StateDescriptionUpPlus = bosons.StateDescriptionUpPlus;
  this->StateDescriptionUpMinus = bosons.StateDescriptionUpMinus;
  this->StateDescriptionDownPlus = bosons.StateDescriptionDownPlus;
  this->StateDescriptionDownMinus = bosons.StateDescriptionDownMinus;

  this->UniqueStateDescriptionUpPlus = bosons.UniqueStateDescriptionUpPlus;
  this->UniqueStateDescriptionSubArraySizeUpPlus = bosons.UniqueStateDescriptionSubArraySizeUpPlus;
  this->NbrUniqueStateDescriptionUpMinus = bosons.NbrUniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionUpMinus = bosons.UniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionSubArraySizeUpMinus = bosons.UniqueStateDescriptionSubArraySizeUpMinus;
  this->FirstIndexUniqueStateDescriptionUpMinus = bosons.FirstIndexUniqueStateDescriptionUpMinus;

  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMaskUpPlus = bosons.LastMomentumMaskUpPlus;
  this->LastMomentumMaskUpMinus = bosons.LastMomentumMaskUpMinus;
  this->LastMomentumMaskDownPlus = bosons.LastMomentumMaskDownPlus;
  this->LastMomentumMaskDownMinus = bosons.LastMomentumMaskDownMinus;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;
  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  this->Flag = bosons.Flag;
}

// destructor
//

BosonOnTorusWithSU4SpinAndMagneticTranslations::~BosonOnTorusWithSU4SpinAndMagneticTranslations ()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
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

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
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

BosonOnTorusWithSU4SpinAndMagneticTranslations& BosonOnTorusWithSU4SpinAndMagneticTranslations::operator = (const BosonOnTorusWithSU4SpinAndMagneticTranslations& bosons)
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

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  delete[] this->TemporaryStateUpPlus;
  delete[] this->TemporaryStateUpMinus;
  delete[] this->TemporaryStateDownPlus;
  delete[] this->TemporaryStateDownMinus;
  delete[] this->ProdATemporaryStateUpPlus;
  delete[] this->ProdATemporaryStateUpMinus;
  delete[] this->ProdATemporaryStateDownPlus;
  delete[] this->ProdATemporaryStateDownMinus;

  this->MomentumModulo = bosons.MomentumModulo;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;

  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalSpin = bosons.TotalSpin;
  this->TotalIsospin = bosons.TotalIsospin;
  this->TotalEntanglement = bosons.TotalEntanglement;
  this->NbrBosonsUpPlus = bosons.NbrBosonsUpPlus;
  this->NbrBosonsUpMinus = bosons.NbrBosonsUpMinus;
  this->NbrBosonsDownPlus = bosons.NbrBosonsDownPlus;
  this->MaxMomentum = bosons.MaxMomentum;  
  this->KyMax = bosons.KyMax;
  this->NUpPlusKyMax = bosons.NUpPlusKyMax;
  this->NUpMinusKyMax = bosons.NUpMinusKyMax;
  this->NDownPlusKyMax = bosons.NDownPlusKyMax;
  this->NDownMinusKyMax = bosons.NDownMinusKyMax;
  this->FermionicKyMax = bosons.FermionicKyMax;

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;

  this->TemporaryStateUpPlus = new unsigned long[this->MaxMomentum];
  this->TemporaryStateUpMinus = new unsigned long[this->MaxMomentum];
  this->TemporaryStateDownPlus = new unsigned long[this->MaxMomentum];
  this->TemporaryStateDownMinus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateUpPlus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateUpMinus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateDownPlus = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateDownMinus = new unsigned long[this->MaxMomentum];

  this->StateDescriptionUpPlus = bosons.StateDescriptionUpPlus;
  this->StateDescriptionUpMinus = bosons.StateDescriptionUpMinus;
  this->StateDescriptionDownPlus = bosons.StateDescriptionDownPlus;
  this->StateDescriptionDownMinus = bosons.StateDescriptionDownMinus;

  this->UniqueStateDescriptionUpPlus = bosons.UniqueStateDescriptionUpPlus;
  this->UniqueStateDescriptionSubArraySizeUpPlus = bosons.UniqueStateDescriptionSubArraySizeUpPlus;
  this->NbrUniqueStateDescriptionUpMinus = bosons.NbrUniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionUpMinus = bosons.UniqueStateDescriptionUpMinus;
  this->UniqueStateDescriptionSubArraySizeUpMinus = bosons.UniqueStateDescriptionSubArraySizeUpMinus;
  this->FirstIndexUniqueStateDescriptionUpMinus = bosons.FirstIndexUniqueStateDescriptionUpMinus;

  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMaskUpPlus = bosons.LastMomentumMaskUpPlus;
  this->LastMomentumMaskUpMinus = bosons.LastMomentumMaskUpMinus;
  this->LastMomentumMaskDownPlus = bosons.LastMomentumMaskDownPlus;
  this->LastMomentumMaskDownMinus = bosons.LastMomentumMaskDownMinus;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;
  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  this->Flag = bosons.Flag;

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnTorusWithSU4SpinAndMagneticTranslations::Clone()
{
  return new BosonOnTorusWithSU4SpinAndMagneticTranslations(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnTorusWithSU4SpinAndMagneticTranslations::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->KxMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->KyMomentum, this->MomentumModulo);
  List<AbstractQuantumNumber*> TmpList2;
  TmpList2 += new VectorQuantumNumber (TmpList);
  return TmpList2;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnTorusWithSU4SpinAndMagneticTranslations::GetQuantumNumber (int index)
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->KxMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->KyMomentum, this->MomentumModulo);
  return new VectorQuantumNumber (TmpList);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnTorusWithSU4SpinAndMagneticTranslations::ExtractSubspace (AbstractQuantumNumber& q, 
										    SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != (AbstractQuantumNumber::Vector | AbstractQuantumNumber::PeriodicMomentum))
    return 0;
  if (((VectorQuantumNumber&) q).GetQuantumNumbers().GetNbrElement() != 2)
    return 0;
  if ((this->KxMomentum == ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[0]))->GetMomentum()) &&
      (this->KyMomentum == ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[1]))->GetMomentum()))
    return this;
  else 
    return 0;
}


// apply a^+_m_up a_m_up operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_up a_m_up

double BosonOnTorusWithSU4SpinAndMagneticTranslations::AdupAup (int index, int m)
{
  cout << "warning : AdupAup not implemented" << endl;
  return 0.0;
}

// apply a^+_m_um a_m_um operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_um a_m_um

double BosonOnTorusWithSU4SpinAndMagneticTranslations::AdumAum (int index, int m)
{
  cout << "warning : AdumAum not implemented" << endl;
  return 0.0;
}

// apply a^+_m_dp a_m_dp operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_dp a_m_dp

double BosonOnTorusWithSU4SpinAndMagneticTranslations::AddpAdp (int index, int m)
{
  cout << "warning : AddpAdp not implemented" << endl;
  return 0.0;
}

// apply a^+_m_dm a_m_dm operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_dm a_m_dm

double BosonOnTorusWithSU4SpinAndMagneticTranslations::AddmAdm (int index, int m)
{
  cout << "warning : AddmAdm not implemented" << endl;
  return 0.0;
}

// apply a^+_m_up a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdupAup (int index, int m, int n, double& coefficient)
{
  cout << "warning : AdupAup not implemented" << endl;
  return 0.0;
}

// apply a^+_m_up a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdupAum (int index, int m, int n, double& coefficient)
{
  cout << "warning : AdupAum not implemented" << endl;
  return 0.0;
}

// apply a^+_m_up a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdupAdp (int index, int m, int n, double& coefficient)
{
  cout << "warning : AdupAdp not implemented" << endl;
  return 0.0;
}

// apply a^+_m_up a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdupAdm (int index, int m, int n, double& coefficient)
{
  cout << "warning : AdupAdm not implemented" << endl;
  return 0.0;
}


// apply a^+_m_um a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdumAup (int index, int m, int n, double& coefficient)
{
  cout << "warning : AdumAup not implemented" << endl;
  return 0.0;
}

// apply a^+_m_um a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdumAum (int index, int m, int n, double& coefficient)
{
  cout << "warning : AdumAum not implemented" << endl;
  return 0.0;
}

// apply a^+_m_um a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdumAdp (int index, int m, int n, double& coefficient)
{
  cout << "warning : AdumAdp not implemented" << endl;
  return 0.0;
}

// apply a^+_m_um a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdumAdm (int index, int m, int n, double& coefficient)
{
  cout << "warning : AdumAdm not implemented" << endl;
  return 0.0;
}

// apply a^+_m_dp a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AddpAup (int index, int m, int n, double& coefficient)
{
  cout << "warning : AddpAup not implemented" << endl;
  return 0.0;
}

// apply a^+_m_dp a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AddpAum (int index, int m, int n, double& coefficient)
{
  cout << "warning : AddpAum not implemented" << endl;
  return 0.0;
}

// apply a^+_m_dp a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AddpAdp (int index, int m, int n, double& coefficient)
{
  cout << "warning : AddpAdp not implemented" << endl;
  return 0.0;
}

// apply a^+_m_dp a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AddpAdm (int index, int m, int n, double& coefficient)
{
  cout << "warning : AddpAdm not implemented" << endl;
  return 0.0;
}

// apply a^+_m_dm a_n_up operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AddmAup (int index, int m, int n, double& coefficient)
{
  cout << "warning : AddmAup not implemented" << endl;
  return 0.0;
}

// apply a^+_m_dm a_n_um operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AddmAum (int index, int m, int n, double& coefficient)
{
  cout << "warning : AddmAum not implemented" << endl;
  return 0.0;
}

// apply a^+_m_dm a_n_dp operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AddmAdp (int index, int m, int n, double& coefficient)
{
  cout << "warning : AddmAdp not implemented" << endl;
  return 0.0;
}

// apply a^+_m_dm a_n_dm operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AddmAdm (int index, int m, int n, double& coefficient)
{
  cout << "warning : AddmAdm not implemented" << endl;
  return 0.0;
}

// apply a_n1_up a_n2_up operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnTorusWithSU4SpinAndMagneticTranslations::AupAup (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusKyMax, this->ProdATemporaryStateUpPlus);
  if ((this->ProdATemporaryStateUpPlus[n1] == 0) || (this->ProdATemporaryStateUpPlus[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateUpPlus[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusKyMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusKyMax, this->ProdATemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusKyMax, this->ProdATemporaryStateDownMinus);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
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
  
double BosonOnTorusWithSU4SpinAndMagneticTranslations::AupAum (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusKyMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusKyMax, this->ProdATemporaryStateUpMinus);
  if ((this->ProdATemporaryStateUpPlus[n1] == 0) || (this->ProdATemporaryStateUpMinus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusKyMax, this->ProdATemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusKyMax, this->ProdATemporaryStateDownMinus);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
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
  
double BosonOnTorusWithSU4SpinAndMagneticTranslations::AupAdp (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusKyMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusKyMax, this->ProdATemporaryStateDownPlus);
  if ((this->ProdATemporaryStateUpPlus[n1] == 0) || (this->ProdATemporaryStateDownPlus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusKyMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusKyMax, this->ProdATemporaryStateDownMinus);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
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
  
double BosonOnTorusWithSU4SpinAndMagneticTranslations::AupAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusKyMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusKyMax, this->ProdATemporaryStateDownMinus);
  if ((this->ProdATemporaryStateUpPlus[n1] == 0) || (this->ProdATemporaryStateDownMinus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusKyMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusKyMax, this->ProdATemporaryStateDownPlus);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryStateDownMinus[n2];
  --this->ProdATemporaryStateDownPlus[n2];
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

double BosonOnTorusWithSU4SpinAndMagneticTranslations::AumAum (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusKyMax, this->ProdATemporaryStateUpMinus);
  if ((this->ProdATemporaryStateUpMinus[n1] == 0) || (this->ProdATemporaryStateUpMinus[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateUpPlus[n2] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusKyMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusKyMax, this->ProdATemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusKyMax, this->ProdATemporaryStateDownMinus);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
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
  
double BosonOnTorusWithSU4SpinAndMagneticTranslations::AumAdp (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusKyMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusKyMax, this->ProdATemporaryStateDownPlus);
  if ((this->ProdATemporaryStateUpMinus[n1] == 0) || (this->ProdATemporaryStateDownPlus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusKyMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusKyMax, this->ProdATemporaryStateDownMinus);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
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
  
double BosonOnTorusWithSU4SpinAndMagneticTranslations::AumAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusKyMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusKyMax, this->ProdATemporaryStateDownMinus);
  if ((this->ProdATemporaryStateUpMinus[n1] == 0) || (this->ProdATemporaryStateDownMinus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusKyMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusKyMax, this->ProdATemporaryStateDownPlus);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
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

double BosonOnTorusWithSU4SpinAndMagneticTranslations::AdpAdp (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusKyMax, this->ProdATemporaryStateDownPlus);
  if ((this->ProdATemporaryStateDownPlus[n1] == 0) || (this->ProdATemporaryStateDownPlus[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateDownPlus[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusKyMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusKyMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusKyMax, this->ProdATemporaryStateDownMinus);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryStateDownPlus[n2];
  --this->ProdATemporaryStateDownPlus[n2];
  Coefficient *= this->ProdATemporaryStateDownPlus[n1];
  --this->ProdATemporaryStateDownPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_dp a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnTorusWithSU4SpinAndMagneticTranslations::AdpAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusKyMax, this->ProdATemporaryStateDownPlus);
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusKyMax, this->ProdATemporaryStateDownMinus);
  if ((this->ProdATemporaryStateDownPlus[n1] == 0) || (this->ProdATemporaryStateDownMinus[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusKyMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusKyMax, this->ProdATemporaryStateUpMinus);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryStateDownMinus[n2];
  --this->ProdATemporaryStateDownMinus[n2];
  Coefficient *= this->ProdATemporaryStateDownPlus[n1];
  --this->ProdATemporaryStateDownPlus[n1];
  return sqrt(Coefficient);
}

// apply a_n1_dp a_n2_dp operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnTorusWithSU4SpinAndMagneticTranslations::AdmAdm (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionDownMinus[index], this->NDownMinusKyMax, this->ProdATemporaryStateDownMinus);
  if ((this->ProdATemporaryStateDownMinus[n1] == 0) || (this->ProdATemporaryStateDownMinus[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateDownMinus[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUpPlus[index], this->NUpPlusKyMax, this->ProdATemporaryStateUpPlus);
  this->FermionToBoson(this->StateDescriptionUpMinus[index], this->NUpMinusKyMax, this->ProdATemporaryStateUpMinus);
  this->FermionToBoson(this->StateDescriptionDownPlus[index], this->NDownPlusKyMax, this->ProdATemporaryStateDownPlus);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
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
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdupAdup (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateUpPlus, coefficient, nbrTranslation);
}

// apply a^+_m1_up a^+_m2_um operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 

int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdupAdum (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, coefficient, nbrTranslation);
}

// apply a^+_m1_up a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
  
int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdupAddp (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateDownPlus, coefficient, nbrTranslation);
}

// apply a^+_m1_up a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
  
int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdupAddm (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpPlus, this->TemporaryStateDownMinus, coefficient, nbrTranslation);
}

// apply a^+_m1_um a^+_m2_um operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
  
int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdumAdum (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpMinus, this->TemporaryStateUpMinus, coefficient, nbrTranslation);
}

// apply a^+_m1_um a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
 
int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdumAddp (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpMinus, this->TemporaryStateDownPlus, coefficient, nbrTranslation);
}

// apply a^+_m1_um a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
 
int BosonOnTorusWithSU4SpinAndMagneticTranslations::AdumAddm (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUpMinus, this->TemporaryStateDownMinus, coefficient, nbrTranslation);
}

// apply a^+_m1_dp a^+_m2_dp operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
  
int BosonOnTorusWithSU4SpinAndMagneticTranslations::AddpAddp (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDownPlus, this->TemporaryStateDownPlus, coefficient, nbrTranslation);
}

// apply a^+_m1_dp a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
  
int BosonOnTorusWithSU4SpinAndMagneticTranslations::AddpAddm (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDownPlus, this->TemporaryStateDownMinus, coefficient, nbrTranslation);
}

// apply a^+_m1_dm a^+_m2_dm operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
  
int BosonOnTorusWithSU4SpinAndMagneticTranslations::AddmAddm (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDownMinus, this->TemporaryStateDownMinus, coefficient, nbrTranslation);
}

// find state index
//
// stateDescriptionUpPlus = unsigned integer describing the fermionic state for type up-plus particles
// stateDescriptionUpMinus = unsigned integer describing the fermionic state for up-minus particles
// stateDescriptionDownPlus = unsigned integer describing the fermionic state for down-plus particles
// stateDescriptionDownMinus = unsigned integer describing the fermionic state for down-minus particles
// return value = corresponding index

int BosonOnTorusWithSU4SpinAndMagneticTranslations::FindStateIndex(unsigned long stateDescriptionUpPlus, unsigned long stateDescriptionUpMinus,
								   unsigned long stateDescriptionDownPlus, unsigned long stateDescriptionDownMinus)
{
  int PosMin = 0;
  int PosMax = this->NbrUniqueStateDescriptionUpPlus - 1;
  int PosMidUpPlus = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->UniqueStateDescriptionUpPlus[PosMidUpPlus];
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

  unsigned long* TmpStateDescriptionArray = this->UniqueStateDescriptionUpMinus[PosMidUpPlus];
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

ostream& BosonOnTorusWithSU4SpinAndMagneticTranslations::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUpPlus[state], this->StateDescriptionUpMinus[state], this->StateDescriptionDownPlus[state], this->StateDescriptionDownMinus[state],
		       this->TemporaryStateUpPlus, this->TemporaryStateUpMinus, this->TemporaryStateDownPlus, this->TemporaryStateDownMinus); 

  unsigned long Tmp;
  Str << " | ";
  for (int i = 0; i <= this->KyMax; ++i)
    {
      Str << "(" << this->TemporaryStateUpPlus[i] << "," << this->TemporaryStateUpMinus[i] << "," << this->TemporaryStateDownPlus[i] << "," << this->TemporaryStateDownMinus[i] << ") | ";
    }
  return Str;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnTorusWithSU4SpinAndMagneticTranslations::GenerateLookUpTable(unsigned long memory)
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
  this->UniqueStateDescriptionUpPlus = new unsigned long [this->NbrUniqueStateDescriptionUpPlus];
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
  this->UniqueStateDescriptionUpMinus = new unsigned long* [this->NbrUniqueStateDescriptionUpPlus];
  this->UniqueStateDescriptionSubArraySizeUpMinus = new int* [this->NbrUniqueStateDescriptionUpPlus];
  this->FirstIndexUniqueStateDescriptionUpMinus = new int* [this->NbrUniqueStateDescriptionUpPlus];
  this->NbrUniqueStateDescriptionDownPlus = new int* [this->NbrUniqueStateDescriptionUpPlus];
  this->UniqueStateDescriptionDownPlus = new unsigned long** [this->NbrUniqueStateDescriptionUpPlus];
  this->UniqueStateDescriptionSubArraySizeDownPlus = new int** [this->NbrUniqueStateDescriptionUpPlus];
  this->FirstIndexUniqueStateDescriptionDownPlus = new int** [this->NbrUniqueStateDescriptionUpPlus];
  for (long i = 0l; i < this->NbrUniqueStateDescriptionUpPlus; ++i)
    {
      this->UniqueStateDescriptionUpMinus[i] = new unsigned long [this->NbrUniqueStateDescriptionUpMinus[i]];
      this->UniqueStateDescriptionSubArraySizeUpMinus[i] = new int [this->NbrUniqueStateDescriptionUpMinus[i]];
      this->FirstIndexUniqueStateDescriptionUpMinus[i] = new int [this->NbrUniqueStateDescriptionUpMinus[i]];
      this->NbrUniqueStateDescriptionDownPlus[i] = new int [this->NbrUniqueStateDescriptionUpMinus[i]]; 
      this->UniqueStateDescriptionDownPlus[i] = new unsigned long* [this->NbrUniqueStateDescriptionUpMinus[i]];
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
	  this->UniqueStateDescriptionDownPlus[i][j] = new unsigned long [TmpUnique];
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

  this->RescalingFactors = new double* [2 * this->MaxMomentum];
  for (int i = 1; i <= this->MaxMomentum; ++i)
    {
      this->RescalingFactors[i] = new double [2 * this->MaxMomentum];
      for (int j = 1; j <= this->MaxMomentum; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}

// generate all states corresponding to the constraints
// 
// return value = hilbert space dimension

long BosonOnTorusWithSU4SpinAndMagneticTranslations::GenerateStates()
{
  unsigned long* TmpStateDescriptionUpPlus = new unsigned long[this->LargeHilbertSpaceDimension];
  unsigned long* TmpStateDescriptionUpMinus = new unsigned long[this->LargeHilbertSpaceDimension];
  unsigned long* TmpStateDescriptionDownPlus = new unsigned long[this->LargeHilbertSpaceDimension];
  unsigned long* TmpStateDescriptionDownMinus = new unsigned long[this->LargeHilbertSpaceDimension];
  long TmpDimension = 0l;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      int NbrTranslation = 0;
      if ((this->FindCanonicalFormAndTestXMomentumConstraint(this->StateDescriptionUpPlus[i], this->StateDescriptionUpMinus[i], 
							     this->StateDescriptionDownPlus[i], this->StateDescriptionDownMinus[i], NbrTranslation) == true) && (NbrTranslation == 0))
	{
	  TmpStateDescriptionUpPlus[TmpDimension] = this->StateDescriptionUpPlus[i];
	  TmpStateDescriptionUpMinus[TmpDimension] = this->StateDescriptionUpMinus[i];
	  TmpStateDescriptionDownPlus[TmpDimension] = this->StateDescriptionDownPlus[i];
	  TmpStateDescriptionDownMinus[TmpDimension] = this->StateDescriptionDownMinus[i];
	  ++TmpDimension;
	}
    }
  delete[] this->StateDescriptionUpPlus;
  delete[] this->StateDescriptionUpMinus;
  delete[] this->StateDescriptionDownPlus;
  delete[] this->StateDescriptionDownMinus;
  this->LargeHilbertSpaceDimension = TmpDimension;
  this->StateDescriptionUpPlus = new unsigned long[this->LargeHilbertSpaceDimension];
  this->StateDescriptionUpMinus = new unsigned long[this->LargeHilbertSpaceDimension];
  this->StateDescriptionDownPlus = new unsigned long[this->LargeHilbertSpaceDimension];
  this->StateDescriptionDownMinus = new unsigned long[this->LargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int[this->LargeHilbertSpaceDimension];
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->StateDescriptionUpPlus[i] = TmpStateDescriptionUpPlus[i];
      this->StateDescriptionUpMinus[i] = TmpStateDescriptionUpMinus[i];
      this->StateDescriptionDownPlus[i] = TmpStateDescriptionDownPlus[i];
      this->StateDescriptionDownMinus[i] = TmpStateDescriptionDownMinus[i];
      unsigned long TmpStateUpPlus = this->StateDescriptionUpPlus[i];
      unsigned long TmpStateUpMinus = this->StateDescriptionUpMinus[i];
      unsigned long TmpStateDownPlus = this->StateDescriptionDownPlus[i];
      unsigned long TmpStateDownMinus = this->StateDescriptionDownMinus[i];
      unsigned long TmpReferenceStateUpPlus = TmpStateUpPlus;
      unsigned long TmpReferenceStateUpMinus = TmpStateUpMinus;
      unsigned long TmpReferenceStateDownPlus = TmpStateDownPlus;
      unsigned long TmpReferenceStateDownMinus = TmpStateDownMinus;
      int TmpOrbitSize = 1;
      this->ApplySingleTranslation(TmpStateUpPlus, TmpStateUpMinus, TmpStateDownPlus, TmpStateDownMinus);
      while ((TmpStateUpPlus != TmpReferenceStateUpPlus) || (TmpStateUpMinus != TmpReferenceStateUpMinus) || 
	     (TmpStateDownPlus != TmpReferenceStateDownPlus) || (TmpStateDownMinus != TmpReferenceStateDownMinus))
	{
	  this->ApplySingleTranslation(TmpStateUpPlus, TmpStateUpMinus, TmpStateDownPlus, TmpStateDownMinus);
	  ++TmpOrbitSize;
	}
      this->NbrStateInOrbit[i] = TmpOrbitSize;
    } 
  delete[] TmpStateDescriptionUpPlus;
  delete[] TmpStateDescriptionUpMinus;
  delete[] TmpStateDescriptionDownPlus;
  delete[] TmpStateDescriptionDownMinus;
  return TmpDimension;
}


// generate all states corresponding to the constraints without the mangetic translations
// 
// nbrBosons = number of bosons
// currentKyUpPlus = current momentum along y for a single type up-plus particle
// currentKyUpMinus = current momentum along y for a single type up-minus particle
// currentKyDownPlus = current momentum along y for a single type down-plus particle
// currentKyDownMinus = current momentum along y for a single type down-minus particle
// currentTotalKy = current total momentum along y
// nbrNUpPlus = number of particles with quantum number up-plus
// nbrNUpMinus = number of particles with quantum number up-minus
// nbrNDownPlus = number of particles with quantum number down-plus
// nbrNDownMinus = number of particles with quantum number down-minus
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnTorusWithSU4SpinAndMagneticTranslations::RawGenerateStates(int nbrBosons, int currentKyUpPlus, int currentKyUpMinus, int currentKyDownPlus, int currentKyDownMinus, 
								       int currentTotalKy, int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus, long pos)
{
//  cout << nbrBosons << " | " << currentKy1 << " " << currentKy2 << " " << currentKy3 << " | " << currentTotalKy << "  | " << nbrNUpPlus << " " << nbrNUpMinus << " " << nbrNDownPlus << " | " << pos << endl;
  if ((nbrBosons < 0) || (nbrNUpPlus < 0) || (nbrNUpMinus < 0) || (nbrNDownPlus < 0) || (nbrNDownMinus < 0))
    return pos;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->MaxMomentum) == this->KyMomentum)
	{
//	  cout << "OK" << endl;
	  this->StateDescriptionUpPlus[pos] = 0x0ul;
	  this->StateDescriptionUpMinus[pos] = 0x0ul;
	  this->StateDescriptionDownPlus[pos] = 0x0ul;
	  this->StateDescriptionDownMinus[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else
	return pos;
    }
  if ((currentKyUpPlus < 0) || (currentKyUpMinus < 0) || (currentKyDownPlus < 0) || (currentKyDownMinus < 0))
    return pos;

  long TmpPos;
  for (int i = nbrNUpPlus; i >= 0; --i)    
    {
      unsigned long MaskUpPlus = ((0x1ul << i) - 1ul) << (currentKyUpPlus + nbrNUpPlus - i);
      for (int j = nbrNUpMinus; j >= 0; --j)
	{
    	  unsigned long MaskUpMinus = ((0x1ul << j) - 1ul) << (currentKyUpMinus + nbrNUpMinus - j);	  
	  for (int k = nbrNDownPlus; k >= 0; --k)
	    {
	      unsigned long MaskDownPlus = ((0x1ul << k) - 1ul) << (currentKyDownPlus + nbrNDownPlus - k);	  
	      for (int l = nbrNDownMinus; l >= 0; --l)
		{
		  unsigned long MaskDownMinus = ((0x1ul << l) - 1ul) << (currentKyDownMinus + nbrNDownMinus - l);	  
		  TmpPos = this->RawGenerateStates(nbrBosons - i - j - k - l, currentKyUpPlus - 1, currentKyUpMinus - 1, currentKyDownPlus - 1, currentKyDownMinus - 1, 
						   currentTotalKy + (currentKyUpPlus * i) + (currentKyUpMinus * j) + (currentKyDownPlus * k) + (currentKyDownMinus * l), 
						   nbrNUpPlus - i, nbrNUpMinus - j, nbrNDownPlus - k, nbrNDownMinus - l, pos); 
		  for (; pos < TmpPos; ++pos)
		    {
		      this->StateDescriptionUpPlus[pos] |= MaskUpPlus;
		      this->StateDescriptionUpMinus[pos] |= MaskUpMinus;
		      this->StateDescriptionDownPlus[pos] |= MaskDownPlus;
		      this->StateDescriptionDownMinus[pos] |= MaskDownMinus;
		    }
		}
	    }
	}
    }
  return pos;
};


// evaluate Hilbert space dimension for a given total spin momentum
//
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// nbrNUpPlus = number of particles with quantum number up-plus
// nbrNUpMinus = number of particles with quantum number up-minus
// nbrNDownPlus = number of particles with quantum number down-plus
// nbrNDownMinus = number of particles with quantum number down-minus
// return value = Hilbert space dimension

long BosonOnTorusWithSU4SpinAndMagneticTranslations::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrNUpPlus, int nbrNUpMinus, int nbrNDownPlus, int nbrNDownMinus)
{
  if ((nbrBosons < 0) || (nbrNUpPlus < 0) || (nbrNUpMinus < 0) || (nbrNDownPlus < 0) || (nbrNDownMinus < 0))
    return 0l;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->MaxMomentum) == this->KyMomentum)
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
	  if (((j + currentTotalKy) % this->MaxMomentum) == this->KyMomentum)
	    ++Tmp;
	}
      return Tmp;
    }
  for (int i = nbrNUpPlus; i >= 0; --i)
    for (int j = nbrNUpMinus; j >= 0; --j)
      for (int k = nbrNDownPlus; k >= 0; --k)
	for (int l = nbrNDownMinus; l >= 0; --l)
	  Tmp += this->EvaluateHilbertSpaceDimension(nbrBosons - (i + j + k + l), currentKy - 1, currentTotalKy + (currentKy * (i + j + k + l)), 
						     nbrNUpPlus - i, nbrNUpMinus - j, nbrNDownPlus - k, nbrNDownMinus - l);
  return  Tmp;
}



