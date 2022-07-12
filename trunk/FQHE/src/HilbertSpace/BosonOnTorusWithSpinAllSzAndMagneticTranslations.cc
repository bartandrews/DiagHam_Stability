////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                  class of boson with spin on a torus taking                //
//            into account magnetic translations and all Sz sectors           //
//                                                                            //
//                        last modification : 21/10/2015                      //
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
#include "HilbertSpace/BosonOnTorusWithSpinAllSzAndMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "HilbertSpace/BosonOnTorusWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "Matrix/Matrix.h"
#include "Matrix/ComplexMatrix.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ArrayTools.h"

#include <cmath>
#include <cstdlib>


using std::cout;
using std::endl;
using std::dec;
using std::hex;
using std::abs;

// constructor when Sz is not conserved
// 
// nbrBosons= number of bosons
// maxMomentum = momentum maximum value for a boson
// kxMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// kyMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)

BosonOnTorusWithSpinAllSzAndMagneticTranslations::BosonOnTorusWithSpinAllSzAndMagneticTranslations (int nbrBosons, int maxMomentum, int kxMomentum, int kyMomentum)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->SzFlag = false;
  this->TotalSpin = 0;
  this->NbrBosonsUp = (this->NbrBosons + this->TotalSpin) / 2;
  this->NbrBosonsDown = (this->NbrBosons - this->TotalSpin) / 2;

  cout << "sizeof " << sizeof(short int) << endl;

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

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->KyMax, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->TargetSpace = this;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUp = new unsigned long[this->LargeHilbertSpaceDimension];
      this->StateDescriptionDown = new unsigned long[this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->RawGenerateStates(this->NbrBosons, this->KyMax, 0, this->KyMax + this->NbrBosons + 1, this->KyMax + this->NbrBosons + 1, 0l);     
      cout  << "Dimension without kx = " << TmpLargeHilbertSpaceDimension << endl;
      for (long i = 0; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->StateDescriptionUp[i];
	  int TmpNbrUps = this->ComputeNbrParticles(this->StateDescriptionUp[i]);
	  this->StateDescriptionUp[i] >>= this->NbrBosons - TmpNbrUps; 
	  this->StateDescriptionDown[i] >>= TmpNbrUps; 
	}
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

BosonOnTorusWithSpinAllSzAndMagneticTranslations::BosonOnTorusWithSpinAllSzAndMagneticTranslations(const BosonOnTorusWithSpinAllSzAndMagneticTranslations& bosons)
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

BosonOnTorusWithSpinAllSzAndMagneticTranslations::~BosonOnTorusWithSpinAllSzAndMagneticTranslations ()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorusWithSpinAllSzAndMagneticTranslations& BosonOnTorusWithSpinAllSzAndMagneticTranslations::operator = (const BosonOnTorusWithSpinAllSzAndMagneticTranslations& bosons)
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

AbstractHilbertSpace* BosonOnTorusWithSpinAllSzAndMagneticTranslations::Clone()
{
  return new BosonOnTorusWithSpinAllSzAndMagneticTranslations(*this);
}


// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnTorusWithSpinAllSzAndMagneticTranslations::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
  this->TargetSpace = (BosonOnTorusWithSpinAllSzAndMagneticTranslations*) targetSpace;
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
int BosonOnTorusWithSpinAllSzAndMagneticTranslations::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
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
int BosonOnTorusWithSpinAllSzAndMagneticTranslations::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
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

int BosonOnTorusWithSpinAllSzAndMagneticTranslations::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  cout << "warning : AddAuAdAu not implemented" << endl;
  return -1;
}


// apply a^+_m_d a_m_d operator to a given state (only spin down)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnTorusWithSpinAllSzAndMagneticTranslations::AddAd (int index, int m)
{
  cout << "warning : AddAd not implemented" << endl;
  return 0.0;
}
  

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnTorusWithSpinAllSzAndMagneticTranslations::AduAu (int index, int m)
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

int BosonOnTorusWithSpinAllSzAndMagneticTranslations::AddAu (int index, int m, double& coefficient, int& nbrTranslation)
{
  int TmpNbrUps = this->ComputeNbrParticles(this->StateDescriptionUp[index]);
  this->FermionToBoson(this->StateDescriptionUp[index], this->MaxMomentum + TmpNbrUps - 1, this->ProdATemporaryStateUp);
  if (this->ProdATemporaryStateUp[m] == 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1, this->ProdATemporaryStateDown);
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

int BosonOnTorusWithSpinAllSzAndMagneticTranslations::AduAd (int index, int m, double& coefficient, int& nbrTranslation)
{
  int TmpNbrUps = this->ComputeNbrParticles(this->StateDescriptionUp[index]);
  this->FermionToBoson(this->StateDescriptionDown[index], this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1, this->ProdATemporaryStateDown);
  if (this->ProdATemporaryStateDown[m] == 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->MaxMomentum + TmpNbrUps - 1, this->ProdATemporaryStateUp);
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

int BosonOnTorusWithSpinAllSzAndMagneticTranslations::AddAu (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  int TmpNbrUps = this->ComputeNbrParticles(this->StateDescriptionUp[index]);
  this->FermionToBoson(this->StateDescriptionUp[index], this->MaxMomentum + TmpNbrUps - 1, this->ProdATemporaryStateUp);
  if (this->ProdATemporaryStateUp[n] == 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1, this->ProdATemporaryStateDown);
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

int BosonOnTorusWithSpinAllSzAndMagneticTranslations::AduAd (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  int TmpNbrUps = this->ComputeNbrParticles(this->StateDescriptionUp[index]);
  this->FermionToBoson(this->StateDescriptionDown[index], this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1, this->ProdATemporaryStateDown);
  if (this->ProdATemporaryStateDown[n] == 0)
    {
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->MaxMomentum + TmpNbrUps - 1, this->ProdATemporaryStateUp);
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
double BosonOnTorusWithSpinAllSzAndMagneticTranslations::AuAu (int index, int n1, int n2)
{
  int TmpNbrUps = this->ComputeNbrParticles(this->StateDescriptionUp[index]);
  this->FermionToBoson(this->StateDescriptionUp[index], this->MaxMomentum + TmpNbrUps - 1, this->ProdATemporaryStateUp);
  if ((this->ProdATemporaryStateUp[n1] == 0) || (this->ProdATemporaryStateUp[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateUp[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1, this->ProdATemporaryStateDown);
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
double BosonOnTorusWithSpinAllSzAndMagneticTranslations::AdAd (int index, int n1, int n2)
{  
  int TmpNbrUps = this->ComputeNbrParticles(this->StateDescriptionUp[index]);
  this->FermionToBoson(this->StateDescriptionDown[index], this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateDown[n1] == 0) || (this->ProdATemporaryStateDown[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateDown[n1] == 1)))    
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->MaxMomentum + TmpNbrUps - 1, this->ProdATemporaryStateUp);
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

double BosonOnTorusWithSpinAllSzAndMagneticTranslations::AuAd (int index, int n1, int n2)
{
  int TmpNbrUps = this->ComputeNbrParticles(this->StateDescriptionUp[index]);
  this->FermionToBoson(this->StateDescriptionUp[index], this->MaxMomentum + TmpNbrUps - 1, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1, this->ProdATemporaryStateDown);
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

// generate all states corresponding to the constraints
// 
// return value = hilbert space dimension

long BosonOnTorusWithSpinAllSzAndMagneticTranslations::GenerateStates()
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
      int TmpNbrUps = this->ComputeNbrParticles(this->StateDescriptionUp[i]);
      unsigned long TmpLastMomentumMaskUp = 0x1ul << (this->MaxMomentum + TmpNbrUps - 1);
      unsigned long TmpLastMomentumMaskDown = 0x1ul << (this->MaxMomentum + this->NbrBosons - TmpNbrUps - 1);
      this->ApplySingleTranslation(TmpStateUp, TmpStateDown, TmpLastMomentumMaskUp, TmpLastMomentumMaskDown);
      while ((TmpStateUp != TmpReferenceStateUp) || (TmpStateDown != TmpReferenceStateDown))
	{
	  this->ApplySingleTranslation(TmpStateUp, TmpStateDown, TmpLastMomentumMaskUp, TmpLastMomentumMaskDown);
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

long BosonOnTorusWithSpinAllSzAndMagneticTranslations::RawGenerateStates(int nbrBosons, int currentKy, int currentTotalKy, 
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

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonOnTorusWithSpinAllSzAndMagneticTranslations::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy)
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
  long Count = 0;
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

// convert state of a SU(2) Hilbert space with fixed Sz to a SU(2) space with all sz sectors
//
// state = state that needs to be projected
// su2space = SU(2) space with fixed sz of the input state
// return value = input state expression in the SU(2) basis

ComplexVector BosonOnTorusWithSpinAllSzAndMagneticTranslations::SU2ToSU2AllSz(ComplexVector& state, ParticleOnSphereWithSpin* su2space)
{
  ComplexVector TmpVector (this->LargeHilbertSpaceDimension, true);
  BosonOnTorusWithSpinAndMagneticTranslations* InputSpace = (BosonOnTorusWithSpinAndMagneticTranslations*) su2space;
  unsigned long TmpStateUp;
  unsigned long TmpStateDown;
  for (long i = 0; i < InputSpace->LargeHilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(InputSpace->StateDescriptionUp[i], InputSpace->NUpKyMax, this->ProdATemporaryStateUp);
      this->FermionToBoson(InputSpace->StateDescriptionDown[i], InputSpace->NDownKyMax, this->ProdATemporaryStateDown);
      this->TargetSpace->BosonToFermion(this->ProdATemporaryStateUp, this->ProdATemporaryStateDown, TmpStateUp, TmpStateDown);
      int TmpIndex = this->TargetSpace->FindStateIndex(TmpStateUp, TmpStateDown);
      if (TmpIndex < this->HilbertSpaceDimension)
	{
	  TmpVector[TmpIndex] = state[i];
	}
    }
  return TmpVector;
}

// convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void BosonOnTorusWithSpinAllSzAndMagneticTranslations::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, 
									     long firstComponent, long nbrComponents)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU2Indices = new int [this->NbrBosons];
  int* TmpSU2Indices2 = new int [this->NbrBosons];
  double* OccupationCoefficientArray = new double [this->NbrBosons + 1];
  OccupationCoefficientArray[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    OccupationCoefficientArray[i] = OccupationCoefficientArray[i - 1] + 0.5 * log((double) i);
  targetState.ClearVector();
  Complex* FourrierCoefficients = new Complex [this->MomentumModulo];
  for (int i = 0; i < this->MomentumModulo; ++i)
    FourrierCoefficients[i] = Phase (2.0 * M_PI * ((double) (i * this->KxMomentum)) / ((double) this->MomentumModulo));
  long LastComponent = firstComponent + nbrComponents;
  if (nbrComponents == 0)
    LastComponent = this->LargeHilbertSpaceDimension;
  for (long i = firstComponent; i < LastComponent; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i],
			   this->TemporaryStateUp, this->TemporaryStateDown); 
      this->ProdANbrStateInOrbit = this->RescalingFactors[this->NbrStateInOrbit[i]];
      double OccupationCoefficient = 0.0;
      int TmpIndex = 0;
      for (int j = this->KyMax; j >= 0; --j)
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
      this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpMomentumIndices, TmpSU2Indices, TmpSU2Indices2, oneBodyBasis, OccupationCoefficient, OccupationCoefficientArray, FourrierCoefficients);
    }
  delete[] OccupationCoefficientArray;
  delete[] TmpMomentumIndices;
  delete[] TmpSU2Indices;
  delete[] TmpSU2Indices2;
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
// fourrierCoefficients = array of Fourrier coefficients for the kx momentum

void BosonOnTorusWithSpinAllSzAndMagneticTranslations::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
										      int position, int* momentumIndices, int* initialSU2Indices, 
										      int* currentSU2Indices, ComplexMatrix* oneBodyBasis,
										      double occupationCoefficient, double* occupationCoefficientArray,
										      Complex* fourrierCoefficients) 
{
  if (position == this->NbrBosons)
    {
      for (int i = 0; i <= this->KyMax; ++i)
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
      int TmpNbrTranslation = 0;
      for (int i = 0; i <= this->KyMax; ++i)
	{
	  occupationCoefficient += occupationCoefficientArray[this->TemporaryStateUp[i]];
	  occupationCoefficient += occupationCoefficientArray[this->TemporaryStateDown[i]];
	}
      unsigned long TmpStateUp;
      unsigned long TmpStateDown;
      this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, TmpStateUp, TmpStateDown);
      if (this->TargetSpace->FindCanonicalFormAndTestXMomentumConstraint(TmpStateUp, TmpStateDown, TmpNbrTranslation) == true)
	{
	  int Index = this->TargetSpace->FindStateIndex(TmpStateUp, TmpStateDown);
	  if (Index < this->HilbertSpaceDimension)
	    {
	      targetState[Index] += (coefficient * this->ProdANbrStateInOrbit[this->TargetSpace->NbrStateInOrbit[Index]] * exp (occupationCoefficient)) * fourrierCoefficients[TmpNbrTranslation];
	    }
	}
      return;      
    }
  else
    {
      currentSU2Indices[position] = 0;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSU2Indices[position]][1]), position + 1, momentumIndices, initialSU2Indices, currentSU2Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray, fourrierCoefficients);
      currentSU2Indices[position] = 1;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSU2Indices[position]][0]), position + 1, momentumIndices, initialSU2Indices, currentSU2Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray, fourrierCoefficients);
    }
}

