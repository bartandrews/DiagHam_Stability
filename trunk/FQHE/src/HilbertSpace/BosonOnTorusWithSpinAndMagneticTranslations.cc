////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                  class of boson with spin on a torus taking                //
//                    into account magnetic translations                      //
//                                                                            //
//                        last modification : 26/04/2012                      //
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
#include "HilbertSpace/BosonOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "HilbertSpace/BosonOnTorusWithSpin.h"
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


// default constructor
// 

BosonOnTorusWithSpinAndMagneticTranslations::BosonOnTorusWithSpinAndMagneticTranslations ()
{
}

// basic constructor
// 
// nbrBosons= number of bosons
// totalSpin = twice the total spin value
// maxMomentum = momentum maximum value for a boson
// kxMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// kyMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)

BosonOnTorusWithSpinAndMagneticTranslations::BosonOnTorusWithSpinAndMagneticTranslations (int nbrBosons, int totalSpin, int maxMomentum, 
											  int kxMomentum, int kyMomentum)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->SzFlag = true;
  this->TotalSpin = totalSpin;
  this->NbrBosonsUp = (this->NbrBosons + this->TotalSpin) / 2;
  this->NbrBosonsDown = (this->NbrBosons - this->TotalSpin) / 2;

  this->MaxMomentum = maxMomentum;  
  this->KyMax = this->MaxMomentum - 1;
  this->MomentumModulo = FindGCD(this->NbrBosons, this->MaxMomentum);
  this->KxMomentum = kxMomentum % this->MomentumModulo;
  this->KyMomentum = kyMomentum % this->MaxMomentum;

  this->TemporaryStateUp = new unsigned long[this->MaxMomentum];
  this->TemporaryStateDown = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateUp = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateDown = new unsigned long[this->MaxMomentum];
  this->NUpKyMax = this->KyMax + this->NbrBosons;
  this->NDownKyMax = this->KyMax + this->NbrBosons;
  this->FermionicKyMax = this->NUpKyMax;

  this->StateShift = (this->MaxMomentum / this->MomentumModulo);
  this->LastMomentumMaskUp = 0x1ul << (this->MaxMomentum + this->NbrBosonsUp - 1);
  this->LastMomentumMaskDown = 0x1ul << (this->MaxMomentum + this->NbrBosonsDown - 1);

  this->MomentumIncrement = (this->NbrBosons * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((unsigned long) 1);
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= ((unsigned long) 1);
    }

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->KyMax, 0, this->NbrBosonsUp);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->TargetSpace = this;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      cout  << "Dimension without Kx = " << this->LargeHilbertSpaceDimension << endl;
      this->StateDescriptionUp = new unsigned long[this->LargeHilbertSpaceDimension];
      this->StateDescriptionDown = new unsigned long[this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->RawGenerateStates(this->NbrBosons, this->KyMax, 0, this->KyMax + this->NbrBosonsUp + 1, this->KyMax + this->NbrBosonsDown + 1, this->NbrBosonsUp, 0l);     
      cout  << "Sorting temporary Hilbert space" << endl;
      SortDoubleElementArrayDownOrdering<unsigned long>(this->StateDescriptionUp, this->StateDescriptionDown, TmpLargeHilbertSpaceDimension);
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
	  UsedMemory += this->LargeHilbertSpaceDimension * (2l * sizeof(unsigned long) + sizeof(int));
	  UsedMemory += this->NbrUniqueStateDescriptionUp * (2l * sizeof(int) + sizeof(unsigned long));
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

BosonOnTorusWithSpinAndMagneticTranslations::BosonOnTorusWithSpinAndMagneticTranslations(const BosonOnTorusWithSpinAndMagneticTranslations& bosons)
{
  this->MomentumModulo = bosons.MomentumModulo;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;

  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalSpin = bosons.TotalSpin;
  this->SzFlag = bosons.SzFlag;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->MaxMomentum = bosons.MaxMomentum;  
  this->KyMax = bosons.KyMax;
  this->NUpKyMax = bosons.NUpKyMax;
  this->NDownKyMax = bosons.NDownKyMax;
  this->FermionicKyMax = bosons.FermionicKyMax;

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;

  this->TemporaryStateUp = new unsigned long[this->MaxMomentum];
  this->TemporaryStateDown = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateUp = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateDown = new unsigned long[this->MaxMomentum];

  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;

  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMaskUp = bosons.LastMomentumMaskUp;
  this->LastMomentumMaskDown = bosons.LastMomentumMaskDown;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;
  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->Flag = bosons.Flag;
}

// destructor
//

BosonOnTorusWithSpinAndMagneticTranslations::~BosonOnTorusWithSpinAndMagneticTranslations ()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescriptionUp;
      delete[] this->StateDescriptionDown;
      delete[] this->UniqueStateDescriptionUp;
      delete[] this->UniqueStateDescriptionSubArraySizeUp;
      delete[] this->FirstIndexUniqueStateDescriptionUp;
      delete[] this->TemporaryStateUp;
      delete[] this->TemporaryStateDown;
      delete[] this->ProdATemporaryStateUp;
      delete[] this->ProdATemporaryStateDown;
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

BosonOnTorusWithSpinAndMagneticTranslations& BosonOnTorusWithSpinAndMagneticTranslations::operator = (const BosonOnTorusWithSpinAndMagneticTranslations& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescriptionUp;
      delete[] this->StateDescriptionDown;
      delete[] this->UniqueStateDescriptionUp;
      delete[] this->UniqueStateDescriptionSubArraySizeUp;
      delete[] this->FirstIndexUniqueStateDescriptionUp;
      delete[] this->TemporaryStateUp;
      delete[] this->TemporaryStateDown;
      delete[] this->ProdATemporaryStateUp;
      delete[] this->ProdATemporaryStateDown;
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
  this->TotalSpin = bosons.TotalSpin;
  this->SzFlag = bosons.SzFlag;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->MaxMomentum = bosons.MaxMomentum;  
  this->KyMax = bosons.KyMax;
  this->NUpKyMax = bosons.NUpKyMax;
  this->NDownKyMax = bosons.NDownKyMax;
  this->FermionicKyMax = bosons.FermionicKyMax;

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;

  this->TemporaryStateUp = new unsigned long[this->MaxMomentum];
  this->TemporaryStateDown = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateUp = new unsigned long[this->MaxMomentum];
  this->ProdATemporaryStateDown = new unsigned long[this->MaxMomentum];

  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;

  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->LastMomentumMaskUp = bosons.LastMomentumMaskUp;
  this->LastMomentumMaskDown = bosons.LastMomentumMaskDown;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;
  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->Flag = bosons.Flag;

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnTorusWithSpinAndMagneticTranslations::Clone()
{
  return new BosonOnTorusWithSpinAndMagneticTranslations(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnTorusWithSpinAndMagneticTranslations::GetQuantumNumbers ()
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

AbstractQuantumNumber* BosonOnTorusWithSpinAndMagneticTranslations::GetQuantumNumber (int index)
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

AbstractHilbertSpace* BosonOnTorusWithSpinAndMagneticTranslations::ExtractSubspace (AbstractQuantumNumber& q, 
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


// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnTorusWithSpinAndMagneticTranslations::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
  this->TargetSpace = (BosonOnTorusWithSpinAndMagneticTranslations*) targetSpace;
}

// apply a^+_(d,m1) a^+_(d,m2) a_(d,n1) a_(d,n2) operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 
int BosonOnTorusWithSpinAndMagneticTranslations::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  cout << "warning : AddAdAdAd not implemented" << endl;
  return -1;
}

// apply a^+_(u,m1) a^+_(u,m2) a_(u,n1) a_(u,n2) operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 
int BosonOnTorusWithSpinAndMagneticTranslations::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  cout << "warning : AduAuAuAu not implemented" << endl;
  return -1;
}

// apply a^+_(u,m1) a^+_(d,m2) a_(d,n1) a_(u,n2) operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to be applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int BosonOnTorusWithSpinAndMagneticTranslations::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  cout << "warning : AddAuAdAu not implemented" << endl;
  return -1;
}


// apply a^+_m_d a_m_d operator to a given state (only spin down)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnTorusWithSpinAndMagneticTranslations::AddAd (int index, int m)
{
  cout << "warning : AddAd not implemented" << endl;
  return 0.0;
}
  

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnTorusWithSpinAndMagneticTranslations::AduAu (int index, int m)
{
  cout << "warning : AduAu not implemented" << endl;
  return 0.0;
}

// apply a^+_m_d a_m_u operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int BosonOnTorusWithSpinAndMagneticTranslations::AddAu (int index, int m, double& coefficient, int& nbrTranslation)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpKyMax, this->ProdATemporaryStateUp);
  if (this->ProdATemporaryStateUp[m] == 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownKyMax, this->ProdATemporaryStateDown);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  coefficient = this->ProdATemporaryStateUp[m];
  --this->ProdATemporaryStateUp[m];
  ++this->ProdATemporaryStateDown[m];
  coefficient *= this->ProdATemporaryStateDown[m];
  unsigned long TmpStateUp;
  unsigned long TmpStateDown;
  this->TargetSpace->BosonToFermion(this->ProdATemporaryStateUp, this->ProdATemporaryStateDown, TmpStateUp, TmpStateDown);
  if (this->TargetSpace->FindCanonicalFormAndTestXMomentumConstraint(TmpStateUp, TmpStateDown, nbrTranslation) == false)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  coefficient = sqrt(coefficient);
  int TmpIndex = this->TargetSpace->FindStateIndex(TmpStateUp, TmpStateDown);
  coefficient *= this->ProdANbrStateInOrbit[this->TargetSpace->NbrStateInOrbit[TmpIndex]];
  nbrTranslation *= this->TargetSpace->StateShift;
  return TmpIndex;
}

// apply a^+_m_u a_m_d operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int BosonOnTorusWithSpinAndMagneticTranslations::AduAd (int index, int m, double& coefficient, int& nbrTranslation)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownKyMax, this->ProdATemporaryStateDown);
  if (this->ProdATemporaryStateDown[m] == 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpKyMax, this->ProdATemporaryStateUp);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  coefficient = this->ProdATemporaryStateDown[m];
  --this->ProdATemporaryStateDown[m];
  ++this->ProdATemporaryStateUp[m];
  coefficient *= this->ProdATemporaryStateUp[m];
  unsigned long TmpStateUp;
  unsigned long TmpStateDown;
  this->TargetSpace->BosonToFermion(this->ProdATemporaryStateUp, this->ProdATemporaryStateDown, TmpStateUp, TmpStateDown);
  if (this->TargetSpace->FindCanonicalFormAndTestXMomentumConstraint(TmpStateUp, TmpStateDown, nbrTranslation) == false)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  coefficient = sqrt(coefficient);
  int TmpIndex = this->TargetSpace->FindStateIndex(TmpStateUp, TmpStateDown);
  coefficient *= this->ProdANbrStateInOrbit[this->TargetSpace->NbrStateInOrbit[TmpIndex]];
  nbrTranslation *= this->TargetSpace->StateShift;
  return TmpIndex;
}

// apply a^+_m_d a_n_u operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int BosonOnTorusWithSpinAndMagneticTranslations::AddAu (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpKyMax, this->ProdATemporaryStateUp);
  if (this->ProdATemporaryStateUp[n] == 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownKyMax, this->ProdATemporaryStateDown);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  coefficient = this->ProdATemporaryStateUp[n];
  --this->ProdATemporaryStateUp[n];
  ++this->ProdATemporaryStateDown[m];
  coefficient *= this->ProdATemporaryStateDown[m];
  unsigned long TmpStateUp;
  unsigned long TmpStateDown;
  this->TargetSpace->BosonToFermion(this->ProdATemporaryStateUp, this->ProdATemporaryStateDown, TmpStateUp, TmpStateDown);
  if (this->TargetSpace->FindCanonicalFormAndTestXMomentumConstraint(TmpStateUp, TmpStateDown, nbrTranslation) == false)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  coefficient = sqrt(coefficient);
  int TmpIndex = this->TargetSpace->FindStateIndex(TmpStateUp, TmpStateDown);
  coefficient *= this->ProdANbrStateInOrbit[this->TargetSpace->NbrStateInOrbit[TmpIndex]];
  nbrTranslation *= this->TargetSpace->StateShift;
  return TmpIndex;
}

// apply a^+_m_u a_n_d operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int BosonOnTorusWithSpinAndMagneticTranslations::AduAd (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownKyMax, this->ProdATemporaryStateDown);
  if (this->ProdATemporaryStateDown[n] == 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpKyMax, this->ProdATemporaryStateUp);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  coefficient = this->ProdATemporaryStateDown[n];
  --this->ProdATemporaryStateDown[n];
  ++this->ProdATemporaryStateUp[m];
  coefficient *= this->ProdATemporaryStateUp[m];
  unsigned long TmpStateUp;
  unsigned long TmpStateDown;
  this->TargetSpace->BosonToFermion(this->ProdATemporaryStateUp, this->ProdATemporaryStateDown, TmpStateUp, TmpStateDown);
  if (this->TargetSpace->FindCanonicalFormAndTestXMomentumConstraint(TmpStateUp, TmpStateDown, nbrTranslation) == false)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  coefficient = sqrt(coefficient);
  int TmpIndex = this->TargetSpace->FindStateIndex(TmpStateUp, TmpStateDown);
  coefficient *= this->ProdANbrStateInOrbit[this->TargetSpace->NbrStateInOrbit[TmpIndex]];
  nbrTranslation *= this->TargetSpace->StateShift;
  return TmpIndex;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 
double BosonOnTorusWithSpinAndMagneticTranslations::AuAu (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpKyMax, this->ProdATemporaryStateUp);
  if ((this->ProdATemporaryStateUp[n1] == 0) || (this->ProdATemporaryStateUp[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateUp[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownKyMax, this->ProdATemporaryStateDown);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryStateUp[n2];
  --this->ProdATemporaryStateUp[n2];
  Coefficient *= this->ProdATemporaryStateUp[n1];
  --this->ProdATemporaryStateUp[n1];
  return sqrt(Coefficient);
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 
double BosonOnTorusWithSpinAndMagneticTranslations::AdAd (int index, int n1, int n2)
{  
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownKyMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateDown[n1] == 0) || (this->ProdATemporaryStateDown[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateDown[n1] == 1)))    
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpKyMax, this->ProdATemporaryStateUp);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryStateDown[n2];
  --this->ProdATemporaryStateDown[n2];
  Coefficient *= this->ProdATemporaryStateDown[n1];
  --this->ProdATemporaryStateDown[n1];
  return sqrt(Coefficient);
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double BosonOnTorusWithSpinAndMagneticTranslations::AuAd (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpKyMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownKyMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateUp[n1] == 0) || (this->ProdATemporaryStateDown[n2] == 0))
    {
      return 0.0;
    }
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
  double Coefficient = this->ProdATemporaryStateDown[n2];
  --this->ProdATemporaryStateDown[n2];
  Coefficient *= this->ProdATemporaryStateUp[n1];
  --this->ProdATemporaryStateUp[n1];
  return sqrt(Coefficient);
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnTorusWithSpinAndMagneticTranslations::AduAdu (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUp, this->TemporaryStateUp, coefficient, nbrTranslation);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnTorusWithSpinAndMagneticTranslations::AddAdd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDown, this->TemporaryStateDown, coefficient, nbrTranslation);
}

  

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnTorusWithSpinAndMagneticTranslations::AduAdd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUp, this->TemporaryStateDown, coefficient, nbrTranslation);
}


// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnTorusWithSpinAndMagneticTranslations::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  double Coefficient = 1.0;
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpKyMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownKyMax, this->ProdATemporaryStateDown);
  this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[index]];
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
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int BosonOnTorusWithSpinAndMagneticTranslations::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient, int& nbrTranslation)
{
  for (int i = 0; i < this->MaxMomentum; ++i)
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
  unsigned long TmpStateUp;
  unsigned long TmpStateDown;
  this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpStateUp, TmpStateDown);
  if (this->FindCanonicalFormAndTestXMomentumConstraint(TmpStateUp, TmpStateDown, nbrTranslation) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpStateUp, TmpStateDown);
  coefficient = sqrt(coefficient);
  coefficient *= this->ProdANbrStateInOrbit[this->NbrStateInOrbit[TmpIndex]];
  nbrTranslation *= this->StateShift;
  return TmpIndex;
}

// find state index
//
// stateDescriptionUp = unsigned integer describing the fermionic state for type up particles
// stateDescriptionDown = unsigned integer describing the fermionic state for type down particles
// return value = corresponding index

int BosonOnTorusWithSpinAndMagneticTranslations::FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
{
  int PosMin = 0;
  int PosMax = this->NbrUniqueStateDescriptionUp - 1;
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
    PosMid = PosMax;

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
  if (CurrentState != stateDescriptionDown)
    return PosMax;
  else
    return PosMid;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusWithSpinAndMagneticTranslations::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUp[state], this->StateDescriptionDown[state],
		       this->TemporaryStateUp, this->TemporaryStateDown); 

  unsigned long Tmp;
  Str << " | ";
  for (int i = 0; i <= this->KyMax; ++i)
    {
      Str << "(" << this->TemporaryStateUp[i] << "," << this->TemporaryStateDown[i] << ") | ";
    }
  return Str;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnTorusWithSpinAndMagneticTranslations::GenerateLookUpTable(unsigned long memory)
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

long BosonOnTorusWithSpinAndMagneticTranslations::GenerateStates()
{
  unsigned long* TmpStateDescriptionUp = new unsigned long[this->LargeHilbertSpaceDimension];
  unsigned long* TmpStateDescriptionDown = new unsigned long[this->LargeHilbertSpaceDimension];
  long TmpDimension = 0l;
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      int NbrTranslation = 0;
      if ((this->FindCanonicalFormAndTestXMomentumConstraint(this->StateDescriptionUp[i], this->StateDescriptionDown[i], NbrTranslation) == true) && (NbrTranslation == 0))
	{
	  TmpStateDescriptionUp[TmpDimension] = this->StateDescriptionUp[i];
	  TmpStateDescriptionDown[TmpDimension] = this->StateDescriptionDown[i];
	  ++TmpDimension;
	}
    }
  delete[] this->StateDescriptionUp;
  delete[] this->StateDescriptionDown;
  this->LargeHilbertSpaceDimension = TmpDimension;
  this->StateDescriptionUp = new unsigned long[this->LargeHilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long[this->LargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int[this->LargeHilbertSpaceDimension];
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->StateDescriptionUp[i] = TmpStateDescriptionUp[i];
      this->StateDescriptionDown[i] = TmpStateDescriptionDown[i];
      unsigned long TmpStateUp = this->StateDescriptionUp[i];
      unsigned long TmpStateDown = this->StateDescriptionDown[i];
      unsigned long TmpReferenceStateUp = TmpStateUp;
      unsigned long TmpReferenceStateDown = TmpStateDown;
      int TmpOrbitSize = 1;
      this->ApplySingleTranslation(TmpStateUp, TmpStateDown);
      while ((TmpStateUp != TmpReferenceStateUp) || (TmpStateDown != TmpReferenceStateDown))
	{
	  this->ApplySingleTranslation(TmpStateUp, TmpStateDown);
	  ++TmpOrbitSize;
	}
      this->NbrStateInOrbit[i] = TmpOrbitSize;
    } 
  delete[] TmpStateDescriptionUp;
  delete[] TmpStateDescriptionDown;
  return TmpDimension;
}


// generate all states corresponding to the constraints without the mangetic translations
// 
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnTorusWithSpinAndMagneticTranslations::RawGenerateStates(int nbrBosons, int currentKy, int currentTotalKy, 
								    int currentFermionicPositionUp, int currentFermionicPositionDown, long pos)
{
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->MaxMomentum) == this->KyMomentum)
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
	  long TmpPos = this->RawGenerateStates(nbrBosons - i - j, currentKy - 1, currentTotalKy + ((i + j) * currentKy), currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, pos);
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

// generate all states corresponding to the constraints without the mangetic translations
// 
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// nbrSpinUp = number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnTorusWithSpinAndMagneticTranslations::RawGenerateStates(int nbrBosons, int currentKy, int currentTotalKy, 
								    int currentFermionicPositionUp, int currentFermionicPositionDown, int nbrSpinUp, long pos)
{
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrBosons))
    return 0l;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->MaxMomentum) == this->KyMomentum)
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
	  long TmpPos = this->RawGenerateStates(nbrBosons - i - j, currentKy - 1, currentTotalKy + ((i + j) * currentKy), currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, nbrSpinUp - i, pos);
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

long BosonOnTorusWithSpinAndMagneticTranslations::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy)
{
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->MaxMomentum) == this->KyMomentum)
	return 1l;
      else	
	return 0l;
    }
  if (currentKy < 0)
    return 0l;
  long Count = 0l;
  if (nbrBosons == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if (((j + currentTotalKy) % this->MaxMomentum) == this->KyMomentum)
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

long BosonOnTorusWithSpinAndMagneticTranslations::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrSpinUp)
{
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrBosons))
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
  long Count = 0l;
  if (nbrBosons == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if (((j + currentTotalKy) % this->MaxMomentum) == this->KyMomentum)
	    Count++;
	}
      return Count;
    }
  for (int i = nbrBosons; i >= 0; --i)
    for (int j = i; j >= 0; --j)
      Count += this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKy - 1, currentTotalKy + (i * currentKy), nbrSpinUp -j);
  return Count;
}

// convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector BosonOnTorusWithSpinAndMagneticTranslations::ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space)
{
  BosonOnTorusWithSpin* TmpSpace = (BosonOnTorusWithSpin*) space;
  ComplexVector TmpVector (TmpSpace->HilbertSpaceDimension, true);
  Complex* FourrierCoefficients = new Complex [this->MomentumModulo];
  for (int i = 0; i < this->MomentumModulo; ++i)
    FourrierCoefficients[i] = Phase (-2.0 * M_PI * ((double) (i * this->KxMomentum)) / ((double) this->MomentumModulo));
  for (int i = 0; i < TmpSpace->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpStateUp = TmpSpace->StateDescriptionUp[i];
      unsigned long TmpStateDown = TmpSpace->StateDescriptionDown[i];
      int NbrTranslation = 0;
      if (this->FindCanonicalFormAndTestXMomentumConstraint(TmpStateUp, TmpStateDown, NbrTranslation) == true)
	{
	  int Pos = this->FindStateIndex(TmpStateUp, TmpStateDown);
	  if (Pos < this->HilbertSpaceDimension)
	    {
	      TmpVector[i] =  state[Pos] * FourrierCoefficients[NbrTranslation] / sqrt((double) this->NbrStateInOrbit[Pos]);
	    }
	}
    }
  delete[] FourrierCoefficients;
  return TmpVector;
}


