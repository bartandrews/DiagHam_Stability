////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//               class of boson with SU3 spin on a torus taking               //
//                    into account magnetic translations                      //
//                                                                            //
//                        last modification : 18/06/2012                      //
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
#include "HilbertSpace/BosonOnTorusWithSU3SpinAndMagneticTranslations.h"
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
// totalTz = twice the total Tz value
// totalY = three time the total Y value
// maxMomentum = momentum maximum value for a boson
// kxMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// kyMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)

BosonOnTorusWithSU3SpinAndMagneticTranslations::BosonOnTorusWithSU3SpinAndMagneticTranslations (int nbrBosons, int totalTz, int totalY, int maxMomentum, 
												int kxMomentum, int kyMomentum)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalY = totalY;
  this->TotalTz = totalTz;
  this->NbrBosons1 = (2 * nbrBosons) + totalY + (3 * totalTz);
  this->NbrBosons2 = (2 * nbrBosons) + totalY - (3 * totalTz);
  this->NbrBosons3 = nbrBosons - totalY;
  this->NbrBosons1 /= 6;
  this->NbrBosons2 /= 6;
  this->NbrBosons3 /= 3;
  
//  cout << this->NbrBosons1 << " " << this->NbrBosons2 << " " << this->NbrBosons3 << " : " << this->TotalTz << " " << this->TotalY << endl;
  this->MaxMomentum = maxMomentum;  
  this->KyMax = this->MaxMomentum - 1;
  this->MomentumModulo = FindGCD(this->NbrBosons, this->MaxMomentum);
  this->KxMomentum = kxMomentum % this->MomentumModulo;
  this->KyMomentum = kyMomentum % this->MaxMomentum;

  this->TemporaryState1 = new unsigned long[this->MaxMomentum];
  this->TemporaryState2 = new unsigned long[this->MaxMomentum];
  this->TemporaryState3 = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryState1 = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryState2 = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryState3 = new unsigned long[this->MaxMomentum];
  this->N1KyMax = this->KyMax + this->NbrBosons;
  this->N2KyMax = this->KyMax + this->NbrBosons;
  this->N3KyMax = this->KyMax + this->NbrBosons;
  this->FermionicKyMax = this->N1KyMax;

  this->StateShift = (this->MaxMomentum / this->MomentumModulo);
  this->LastMomentumMask1 = 0x1ul << (this->MaxMomentum + this->NbrBosons1 - 1);
  this->LastMomentumMask2 = 0x1ul << (this->MaxMomentum + this->NbrBosons2 - 1);
  this->LastMomentumMask3 = 0x1ul << (this->MaxMomentum + this->NbrBosons3 - 1);

  this->MomentumIncrement = (this->NbrBosons * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((unsigned long) 1);
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= ((unsigned long) 1);
    }

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->KyMax, 0, this->NbrBosons1, this->NbrBosons2, this->NbrBosons3);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescription1 = new unsigned long[this->LargeHilbertSpaceDimension];
      this->StateDescription2 = new unsigned long[this->LargeHilbertSpaceDimension];
      this->StateDescription3 = new unsigned long[this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->RawGenerateStates(this->NbrBosons, this->KyMax, this->KyMax, this->KyMax, 0, this->NbrBosons1, this->NbrBosons2, this->NbrBosons3, 0l);     
      SortTripleElementArrayDownOrdering<unsigned long>(this->StateDescription1, this->StateDescription2, this->StateDescription3, this->LargeHilbertSpaceDimension);
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

BosonOnTorusWithSU3SpinAndMagneticTranslations::BosonOnTorusWithSU3SpinAndMagneticTranslations(const BosonOnTorusWithSU3SpinAndMagneticTranslations& bosons)
{
  this->MomentumModulo = bosons.MomentumModulo;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;

  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->NbrBosons1 = bosons.NbrBosons1;
  this->NbrBosons2 = bosons.NbrBosons2;
  this->NbrBosons3 = bosons.NbrBosons3;
  this->MaxMomentum = bosons.MaxMomentum;  
  this->KyMax = bosons.KyMax;
  this->N1KyMax = bosons.N1KyMax;
  this->N2KyMax = bosons.N2KyMax;
  this->N3KyMax = bosons.N3KyMax;
  this->FermionicKyMax = bosons.FermionicKyMax;

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;

  this->TemporaryState1 = new unsigned long[this->MaxMomentum];
  this->TemporaryState2 = new unsigned long[this->MaxMomentum];
  this->TemporaryState3 = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryState1 = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryState2 = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryState3 = new unsigned long[this->MaxMomentum];

  this->StateDescription1 = bosons.StateDescription1;
  this->StateDescription2 = bosons.StateDescription2;
  this->StateDescription3 = bosons.StateDescription3;

  this->NbrUniqueStateDescription1 = bosons.NbrUniqueStateDescription1;
  this->UniqueStateDescription1 = bosons.UniqueStateDescription1;
  this->UniqueStateDescriptionSubArraySize1 = bosons.UniqueStateDescriptionSubArraySize1;
  this->NbrUniqueStateDescription2 = bosons.NbrUniqueStateDescription2;
  this->UniqueStateDescription2 = bosons.UniqueStateDescription2;
  this->UniqueStateDescriptionSubArraySize2 = bosons.UniqueStateDescriptionSubArraySize2;
  this->FirstIndexUniqueStateDescription2 = bosons.FirstIndexUniqueStateDescription2;

  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMask1 = bosons.LastMomentumMask1;
  this->LastMomentumMask2 = bosons.LastMomentumMask2;
  this->LastMomentumMask3 = bosons.LastMomentumMask3;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;
  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  this->Flag = bosons.Flag;
}

// destructor
//

BosonOnTorusWithSU3SpinAndMagneticTranslations::~BosonOnTorusWithSU3SpinAndMagneticTranslations ()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
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
      delete[] this->ProdATemporaryState1;
      delete[] this->ProdATemporaryState2;
      delete[] this->ProdATemporaryState3;
      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorusWithSU3SpinAndMagneticTranslations& BosonOnTorusWithSU3SpinAndMagneticTranslations::operator = (const BosonOnTorusWithSU3SpinAndMagneticTranslations& bosons)
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
      delete[] this->ProdATemporaryState1;
      delete[] this->ProdATemporaryState2;
      delete[] this->ProdATemporaryState3;
      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }

  this->MomentumModulo = bosons.MomentumModulo;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;

  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->NbrBosons1 = bosons.NbrBosons1;
  this->NbrBosons2 = bosons.NbrBosons2;
  this->NbrBosons3 = bosons.NbrBosons3;
  this->MaxMomentum = bosons.MaxMomentum;  
  this->KyMax = bosons.KyMax;
  this->N1KyMax = bosons.N1KyMax;
  this->N2KyMax = bosons.N2KyMax;
  this->N3KyMax = bosons.N3KyMax;
  this->FermionicKyMax = bosons.FermionicKyMax;

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;

  this->TemporaryState1 = new unsigned long[this->MaxMomentum];
  this->TemporaryState2 = new unsigned long[this->MaxMomentum];
  this->TemporaryState3 = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryState1 = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryState2 = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryState3 = new unsigned long[this->MaxMomentum];

  this->StateDescription1 = bosons.StateDescription1;
  this->StateDescription2 = bosons.StateDescription2;
  this->StateDescription3 = bosons.StateDescription3;

  this->NbrUniqueStateDescription1 = bosons.NbrUniqueStateDescription1;
  this->UniqueStateDescription1 = bosons.UniqueStateDescription1;
  this->UniqueStateDescriptionSubArraySize1 = bosons.UniqueStateDescriptionSubArraySize1;
  this->NbrUniqueStateDescription2 = bosons.NbrUniqueStateDescription2;
  this->UniqueStateDescription2 = bosons.UniqueStateDescription2;
  this->UniqueStateDescriptionSubArraySize2 = bosons.UniqueStateDescriptionSubArraySize2;
  this->FirstIndexUniqueStateDescription2 = bosons.FirstIndexUniqueStateDescription2;

  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMask1 = bosons.LastMomentumMask1;
  this->LastMomentumMask2 = bosons.LastMomentumMask2;
  this->LastMomentumMask3 = bosons.LastMomentumMask3;
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

AbstractHilbertSpace* BosonOnTorusWithSU3SpinAndMagneticTranslations::Clone()
{
  return new BosonOnTorusWithSU3SpinAndMagneticTranslations(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnTorusWithSU3SpinAndMagneticTranslations::GetQuantumNumbers ()
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

AbstractQuantumNumber* BosonOnTorusWithSU3SpinAndMagneticTranslations::GetQuantumNumber (int index)
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

AbstractHilbertSpace* BosonOnTorusWithSU3SpinAndMagneticTranslations::ExtractSubspace (AbstractQuantumNumber& q, 
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


// apply a^+_m_1 a_m_1 operator to a given state (only state 1 Tz=+1/2, Y=+1/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_1 a_m_1

double BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad1A1 (int index, int m)
{
  cout << "warning : Ad1A1 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_2 a_m_2 operator to a given state (only state 2 Tz=-1/2, Y=+1/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_2 a_m_2

double BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad2A2 (int index, int m)
{
  cout << "warning : Ad2A2 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_3 a_m_3 operator to a given state (only state 3 Tz=0, Y=-2/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_3 a_m_3

double BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad3A3 (int index, int m)
{
  cout << "warning : Ad3A3 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_1 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad1A1 (int index, int m, int n, double& coefficient)
{
  cout << "warning : Ad1A1 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_1 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad1A2 (int index, int m, int n, double& coefficient)
{
  cout << "warning : Ad1A2 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_1 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad1A3 (int index, int m, int n, double& coefficient)
{
  cout << "warning : Ad1A2 not implemented" << endl;
  return 0.0;
}


// apply a^+_m_2 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad2A1 (int index, int m, int n, double& coefficient)
{
  cout << "warning : Ad2A1 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_2 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad2A2 (int index, int m, int n, double& coefficient)
{
  cout << "warning : Ad2A2 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_2 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad2A3 (int index, int m, int n, double& coefficient)
{
  cout << "warning : Ad2A3 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_3 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad3A1 (int index, int m, int n, double& coefficient)
{
  cout << "warning : Ad3A1 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_3 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad3A2 (int index, int m, int n, double& coefficient)
{
  cout << "warning : Ad3A2 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_3 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad3A3 (int index, int m, int n, double& coefficient)
{
  cout << "warning : Ad3A3 not implemented" << endl;
  return 0.0;
}

// apply a_n1_1 a_n2_1 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnTorusWithSU3SpinAndMagneticTranslations::A1A1 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1KyMax, this->ProdATemporaryState1);
  if ((this->ProdATemporaryState1[n1] == 0) || (this->ProdATemporaryState1[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryState1[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription2[index], this->N2KyMax, this->ProdATemporaryState2);
  this->FermionToBoson(this->StateDescription3[index], this->N3KyMax, this->ProdATemporaryState3);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryState1[n2];
  --this->ProdATemporaryState1[n2];
  Coefficient *= this->ProdATemporaryState1[n1];
  --this->ProdATemporaryState1[n1];
  return sqrt(Coefficient);
}

// apply a_n1_1 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 
  
double BosonOnTorusWithSU3SpinAndMagneticTranslations::A1A2 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1KyMax, this->ProdATemporaryState1);
  this->FermionToBoson(this->StateDescription2[index], this->N2KyMax, this->ProdATemporaryState2);
  if ((this->ProdATemporaryState1[n1] == 0) || (this->ProdATemporaryState2[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription3[index], this->N3KyMax, this->ProdATemporaryState3);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryState2[n2];
  --this->ProdATemporaryState2[n2];
  Coefficient *= this->ProdATemporaryState1[n1];
  --this->ProdATemporaryState1[n1];
  return sqrt(Coefficient);
}

// apply a_n1_1 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 
  
double BosonOnTorusWithSU3SpinAndMagneticTranslations::A1A3 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription1[index], this->N1KyMax, this->ProdATemporaryState1);
  this->FermionToBoson(this->StateDescription3[index], this->N3KyMax, this->ProdATemporaryState3);
  if ((this->ProdATemporaryState1[n1] == 0) || (this->ProdATemporaryState3[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription2[index], this->N2KyMax, this->ProdATemporaryState2);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryState3[n2];
  --this->ProdATemporaryState3[n2];
  Coefficient *= this->ProdATemporaryState1[n1];
  --this->ProdATemporaryState1[n1];
  return sqrt(Coefficient);
}

// apply a_n1_2 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnTorusWithSU3SpinAndMagneticTranslations::A2A2 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription2[index], this->N2KyMax, this->ProdATemporaryState2);
  if ((this->ProdATemporaryState2[n1] == 0) || (this->ProdATemporaryState2[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryState1[n2] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription1[index], this->N1KyMax, this->ProdATemporaryState1);
  this->FermionToBoson(this->StateDescription3[index], this->N3KyMax, this->ProdATemporaryState3);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryState2[n2];
  --this->ProdATemporaryState2[n2];
  Coefficient *= this->ProdATemporaryState2[n1];
  --this->ProdATemporaryState2[n1];
  return sqrt(Coefficient);
}

// apply a_n1_2 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 
  
double BosonOnTorusWithSU3SpinAndMagneticTranslations::A2A3 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription2[index], this->N2KyMax, this->ProdATemporaryState2);
  this->FermionToBoson(this->StateDescription3[index], this->N3KyMax, this->ProdATemporaryState3);
  if ((this->ProdATemporaryState2[n1] == 0) || (this->ProdATemporaryState3[n2] == 0))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription1[index], this->N1KyMax, this->ProdATemporaryState1);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryState3[n2];
  --this->ProdATemporaryState3[n2];
  Coefficient *= this->ProdATemporaryState2[n1];
  --this->ProdATemporaryState2[n1];
  return sqrt(Coefficient);
}

// apply a_n1_3 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnTorusWithSU3SpinAndMagneticTranslations::A3A3 (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescription3[index], this->N3KyMax, this->ProdATemporaryState3);
  if ((this->ProdATemporaryState3[n1] == 0) || (this->ProdATemporaryState3[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryState3[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescription1[index], this->N1KyMax, this->ProdATemporaryState1);
  this->FermionToBoson(this->StateDescription2[index], this->N2KyMax, this->ProdATemporaryState2);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryState3[n2];
  --this->ProdATemporaryState3[n2];
  Coefficient *= this->ProdATemporaryState3[n1];
  --this->ProdATemporaryState3[n1];
  return sqrt(Coefficient);
}

// apply a^+_m1_1 a^+_m2_1 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 

int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad1Ad1 (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryState1, this->TemporaryState1, coefficient, nbrTranslation);
}

// apply a^+_m1_1 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 

int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad1Ad2 (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryState1, this->TemporaryState2, coefficient, nbrTranslation);
}

// apply a^+_m1_1 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
  
int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad1Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryState1, this->TemporaryState3, coefficient, nbrTranslation);
}

// apply a^+_m1_2 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
  
int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad2Ad2 (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryState2, this->TemporaryState2, coefficient, nbrTranslation);
}

// apply a^+_m1_2 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
 
int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad2Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryState2, this->TemporaryState3, coefficient, nbrTranslation);
}

// apply a^+_m1_3 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = number of translations needed to find the canonical state
// return value = index of the destination state 
  
int BosonOnTorusWithSU3SpinAndMagneticTranslations::Ad3Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryState3, this->TemporaryState3, coefficient, nbrTranslation);
}

// find state index
//
// stateDescription1 = unsigned integer describing the fermionic state for type 1 particles
// stateDescription2 = unsigned integer describing the fermionic state for type 2 particles
// stateDescription3 = unsigned integer describing the fermionic state for type 3 particles
// return value = corresponding index

int BosonOnTorusWithSU3SpinAndMagneticTranslations::FindStateIndex(unsigned long stateDescription1, unsigned long stateDescription2, unsigned long stateDescription3)
{
  int PosMin = 0;
  int PosMax = this->NbrUniqueStateDescription1 - 1;
  int PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->UniqueStateDescription1[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescription1))
    {
       if (CurrentState > stateDescription1)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->UniqueStateDescription1[PosMid];
    }
  if (CurrentState != stateDescription1)
    PosMid = PosMax;
  unsigned long* TmpStateDescriptionArray = this->UniqueStateDescription2[PosMid];
  int* TmpFirstIndexUniqueStateDescription2 = this->FirstIndexUniqueStateDescription2[PosMid];
  int* TmpUniqueStateDescriptionSubArraySize2 = this->UniqueStateDescriptionSubArraySize2[PosMid];
  PosMin = 0;
  PosMax = this->NbrUniqueStateDescription2[PosMid] - 1;
  PosMid = (PosMin + PosMax) >> 1;
  CurrentState = TmpStateDescriptionArray[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescription2))
    {
       if (CurrentState > stateDescription2)
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
  if (CurrentState != stateDescription2)
    PosMid = PosMax;
  PosMin = TmpFirstIndexUniqueStateDescription2[PosMid];
  PosMax = PosMin + TmpUniqueStateDescriptionSubArraySize2[PosMid] - 1;
  PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->StateDescription3[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescription3))
    {
       if (CurrentState > stateDescription3)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->StateDescription3[PosMid];
    }
  if (CurrentState != stateDescription3)
    return PosMax;
  else
    return PosMid;
}
// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusWithSU3SpinAndMagneticTranslations::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescription1[state], this->StateDescription2[state], this->StateDescription3[state],
		       this->TemporaryState1, this->TemporaryState2, this->TemporaryState3); 

  unsigned long Tmp;
  Str << " | ";
  for (int i = 0; i <= this->KyMax; ++i)
    {
      Str << "(" << this->TemporaryState1[i] << "," << this->TemporaryState2[i] << "," << this->TemporaryState3[i] << ") | ";
    }
  int TmpIndex = this->FindStateIndex(this->StateDescription1[state], this->StateDescription2[state], this->StateDescription3[state]);
  Str << state << " " << TmpIndex;
  if (state != TmpIndex)
    Str << " error";
  Str << " nbr states in orbit " << this->NbrStateInOrbit[state] << endl;
  return Str;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnTorusWithSU3SpinAndMagneticTranslations::GenerateLookUpTable(unsigned long memory)
{  
  long TmpUniquePartition = 1l;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescription1[i - 1] == this->StateDescription1[i]))
	{
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	++TmpUniquePartition;
    }

  this->NbrUniqueStateDescription1 = TmpUniquePartition;
  this->UniqueStateDescription1 = new unsigned long [this->NbrUniqueStateDescription1];
  this->UniqueStateDescriptionSubArraySize1 = new int [this->NbrUniqueStateDescription1];
  TmpUniquePartition = 0l;
  this->UniqueStateDescription1[0l] = this->StateDescription1[0l];
  this->UniqueStateDescriptionSubArraySize1[0l] = 1;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescription1[i - 1] == this->StateDescription1[i]))
	{
	  ++this->UniqueStateDescriptionSubArraySize1[TmpUniquePartition];
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	{
	  ++TmpUniquePartition;
	  this->UniqueStateDescription1[TmpUniquePartition] = this->StateDescription1[i];
	  this->UniqueStateDescriptionSubArraySize1[TmpUniquePartition] = 1; 
	}
    }

  this->NbrUniqueStateDescription2 = new int [this->NbrUniqueStateDescription1];
  TmpUniquePartition = 0;
  long TmpIndex = 0l;
  while (TmpIndex < this->LargeHilbertSpaceDimension)
    {
      long Lim = TmpIndex + this->UniqueStateDescriptionSubArraySize1[TmpUniquePartition];
      this->NbrUniqueStateDescription2[TmpUniquePartition] = 1;
      ++TmpIndex;
      while (TmpIndex < Lim)
	{
	  while ((TmpIndex < Lim) && (this->StateDescription2[TmpIndex - 1] == this->StateDescription2[TmpIndex]))
	    ++TmpIndex;
	  if (TmpIndex < Lim)
	    {
	      ++this->NbrUniqueStateDescription2[TmpUniquePartition];
	      ++TmpIndex;
	    }
	}
      ++TmpUniquePartition;
    }
  this->UniqueStateDescription2 = new unsigned long* [this->NbrUniqueStateDescription1];
  this->UniqueStateDescriptionSubArraySize2 = new int* [this->NbrUniqueStateDescription1];
  this->FirstIndexUniqueStateDescription2 = new int* [this->NbrUniqueStateDescription1];
  for (long i = 0l; i < this->NbrUniqueStateDescription1; ++i)
    {
      this->UniqueStateDescription2[i] = new unsigned long [this->NbrUniqueStateDescription2[i]];
      this->UniqueStateDescriptionSubArraySize2[i] = new int [this->NbrUniqueStateDescription2[i]];
      this->FirstIndexUniqueStateDescription2[i] = new int [this->NbrUniqueStateDescription2[i]];
    }

  TmpUniquePartition = 0;
  TmpIndex = 0l;
  while (TmpIndex < this->LargeHilbertSpaceDimension)
    {
      long Lim = TmpIndex + this->UniqueStateDescriptionSubArraySize1[TmpUniquePartition];
      int TmpUniquePartition2 = 0;
      this->UniqueStateDescription2[TmpUniquePartition][TmpUniquePartition2] = this->StateDescription2[TmpIndex];
      this->UniqueStateDescriptionSubArraySize2[TmpUniquePartition][TmpUniquePartition2] = 1;
      this->FirstIndexUniqueStateDescription2[TmpUniquePartition][TmpUniquePartition2] = TmpIndex;
      ++TmpIndex;
      while (TmpIndex < Lim)
	{
	  while ((TmpIndex < Lim) && (this->StateDescription2[TmpIndex - 1] == this->StateDescription2[TmpIndex]))
	    {
	      ++this->UniqueStateDescriptionSubArraySize2[TmpUniquePartition][TmpUniquePartition2];	      
	      ++TmpIndex;
	    }
	  if (TmpIndex < Lim)
	    {
	      ++TmpUniquePartition2;
	      this->UniqueStateDescription2[TmpUniquePartition][TmpUniquePartition2] = this->StateDescription2[TmpIndex];
	      this->UniqueStateDescriptionSubArraySize2[TmpUniquePartition][TmpUniquePartition2] = 1;
	      this->FirstIndexUniqueStateDescription2[TmpUniquePartition][TmpUniquePartition2] = TmpIndex;
	      ++TmpIndex;
	    }
	}
      ++TmpUniquePartition;
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

long BosonOnTorusWithSU3SpinAndMagneticTranslations::GenerateStates()
{
  unsigned long* TmpStateDescription1 = new unsigned long[this->LargeHilbertSpaceDimension];
  unsigned long* TmpStateDescription2 = new unsigned long[this->LargeHilbertSpaceDimension];
  unsigned long* TmpStateDescription3 = new unsigned long[this->LargeHilbertSpaceDimension];
  long TmpDimension = 0l;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      int NbrTranslation = 0;
      if ((this->FindCanonicalFormAndTestXMomentumConstraint(this->StateDescription1[i], this->StateDescription2[i], this->StateDescription3[i], NbrTranslation) == true) && (NbrTranslation == 0))
	{
	  TmpStateDescription1[TmpDimension] = this->StateDescription1[i];
	  TmpStateDescription2[TmpDimension] = this->StateDescription2[i];
	  TmpStateDescription3[TmpDimension] = this->StateDescription3[i];
	  ++TmpDimension;
	}
    }
  delete[] this->StateDescription1;
  delete[] this->StateDescription2;
  delete[] this->StateDescription3;
  this->LargeHilbertSpaceDimension = TmpDimension;
  this->StateDescription1 = new unsigned long[this->LargeHilbertSpaceDimension];
  this->StateDescription2 = new unsigned long[this->LargeHilbertSpaceDimension];
  this->StateDescription3 = new unsigned long[this->LargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int[this->LargeHilbertSpaceDimension];
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->StateDescription1[i] = TmpStateDescription1[i];
      this->StateDescription2[i] = TmpStateDescription2[i];
      this->StateDescription3[i] = TmpStateDescription3[i];
      unsigned long TmpState1 = this->StateDescription1[i];
      unsigned long TmpState2 = this->StateDescription2[i];
      unsigned long TmpState3 = this->StateDescription3[i];
      unsigned long TmpReferenceState1 = TmpState1;
      unsigned long TmpReferenceState2 = TmpState2;
      unsigned long TmpReferenceState3 = TmpState3;
      int TmpOrbitSize = 1;
      this->ApplySingleTranslation(TmpState1, TmpState2, TmpState3);
      while ((TmpState1 != TmpReferenceState1) || (TmpState2 != TmpReferenceState2) || (TmpState3 != TmpReferenceState3))
	{
	  this->ApplySingleTranslation(TmpState1, TmpState2, TmpState3);
	  ++TmpOrbitSize;
	}
      this->NbrStateInOrbit[i] = TmpOrbitSize;
    } 
  delete[] TmpStateDescription1;
  delete[] TmpStateDescription2;
  delete[] TmpStateDescription3;
  return TmpDimension;
}


// generate all states corresponding to the constraints without the mangetic translations
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

long BosonOnTorusWithSU3SpinAndMagneticTranslations::RawGenerateStates(int nbrBosons, int currentKy1, int currentKy2, int currentKy3, int currentTotalKy, 
								       int nbrN1, int nbrN2, int nbrN3, long pos)
{
  if ((nbrBosons < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0))
    return pos;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->MaxMomentum) == this->KyMomentum)
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
	      TmpPos = this->RawGenerateStates(nbrBosons - i - j - k, currentKy1 - 1, currentKy2 - 1, currentKy3 - 1, currentTotalKy + (currentKy1 * i) + (currentKy2 * j) + (currentKy3 * k), 
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





//   long TmpPos;
//   unsigned long Mask1;
//   unsigned long Mask2;
//   unsigned long Mask3;

//   if (nbrN1 == 0)
//     {
//       if (nbrN2 == 0)
// 	{
// 	  for (int k = nbrN3; k > 0; --k)
// 	    {
// 	      TmpPos = this->RawGenerateStates(nbrBosons - k, 0, 0, currentKy3 - 1, currentTotalKy + (currentKy3 * k), 
// 					       0, 0, nbrN3 - k, pos); 
// 	      Mask3 = ((0x1ul << k) - 1ul) << (currentKy3 + nbrN3 - k);
// 	      for (; pos < TmpPos; ++pos)
// 		{
// 		  this->StateDescription3[pos] |= Mask3;
// 		}
// 	    }
// 	  pos = this->RawGenerateStates(nbrBosons, 0, 0, currentKy3 - 1, currentTotalKy, 0, 0, nbrN3, pos);
// 	  return pos;
// 	}
//       TmpPos = this->RawGenerateStates(nbrBosons - nbrN2, 0, 0, currentKy3, currentTotalKy + (currentKy2 * nbrN2), 
// 				    0, 0, nbrN3, pos); 
//       Mask2 = ((0x1ul << nbrN2) - 1ul) << currentKy2;
//       for (; pos < TmpPos; ++pos)
// 	this->StateDescription2[pos] |= Mask2;
//       for (int j = nbrN2 - 1; j > 0; --j)
// 	{
// 	  TmpPos = this->RawGenerateStates(nbrBosons - j, 0, currentKy2 - 1, currentKy3, currentTotalKy + (currentKy2 * j), 
// 					   0, nbrN2 - j, nbrN3, pos); 
// 	  Mask2 = ((0x1ul << j) - 1ul) << (currentKy2 + nbrN2 - j);
// 	  for (; pos < TmpPos; ++pos)
// 	    this->StateDescription2[pos] |= Mask2;
// 	}
//       pos = this->RawGenerateStates(nbrBosons, 0, currentKy2 - 1, currentKy3, currentTotalKy, 0, nbrN2, nbrN3, pos);
//       return pos;
//     }
  
//   TmpPos = this->RawGenerateStates(nbrBosons - nbrN1, 0, currentKy2, currentKy3, currentTotalKy + (currentKy1 * nbrN1), 
// 				   0, nbrN2, nbrN3, pos); 
//   Mask1 = ((0x1ul << nbrN1) - 1ul) << currentKy1;
//   for (; pos < TmpPos; ++pos)
//     this->StateDescription1[pos] |= Mask1;
//   for (int i = nbrN1 - 1; i > 0; --i)
//     {
//       TmpPos = this->RawGenerateStates(nbrBosons - i, currentKy1 - 1, currentKy2, currentKy3, currentTotalKy + (currentKy1 * i), 
// 				       nbrN1 - i, nbrN2, nbrN3, pos); 
//       Mask1 = ((0x1ul << i) - 1ul) << (currentKy1 + nbrN1 - i);
//       for (; pos < TmpPos; ++pos)
// 	{
// 	  this->StateDescription1[pos] |= Mask1;
// 	}
//     }
//   pos = this->RawGenerateStates(nbrBosons, currentKy1 - 1, currentKy2, currentKy3, currentTotalKy, nbrN1, nbrN2, nbrN3, pos);
//   return pos;
};


// evaluate Hilbert space dimension for a given total spin momentum
//
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// return value = Hilbert space dimension

long BosonOnTorusWithSU3SpinAndMagneticTranslations::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrN1, int nbrN2, int nbrN3)
{
  if ((nbrBosons < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0))
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
  for (int i = nbrN1; i >= 0; --i)
    for (int j = nbrN2; j >= 0; --j)
      for (int k = nbrN3; k >= 0; --k)
	Tmp += this->EvaluateHilbertSpaceDimension(nbrBosons - (i + j + k), currentKy - 1, currentTotalKy + (currentKy * (i + j + k)), 
						   nbrN1 - i, nbrN2 - j, nbrN3 - k);
  return  Tmp;
}



