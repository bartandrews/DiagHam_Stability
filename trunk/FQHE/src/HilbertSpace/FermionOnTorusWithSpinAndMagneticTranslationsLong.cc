////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of fermion with spin on a torus taking               //
//     into account magnetic translations involving more than 32 orbitals     //
//                                                                            //
//                        last modification : 16/06/2016                      //
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
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslationsLong.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorusWithSpin.h"
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

FermionOnTorusWithSpinAndMagneticTranslationsLong::FermionOnTorusWithSpinAndMagneticTranslationsLong ()
{
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrFermions = 0;
  this->TotalSpin = 0;
  this->IncNbrFermions = 0;

  this->MaxMomentum = 0;
  this->NbrMomentum = 0;
  this->NbrFermionStates = 0;
  this->MomentumModulo = 0;
  this->XMomentum = 0;
  this->YMomentum = 0;

  this->MomentumIncrement = 0;
  this->StateShift = 0;
  this->ComplementaryStateShift = 0;
  this->MomentumMask = ((ULONGLONG) 0x0ul);

  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateHighestBit = 0;

  this->MaximumLookUpShift = 0;
  this->LookUpTableMemorySize = 0;
  this->LookUpTableShift = 0;
  this->LookUpTable = 0;

  this->SignLookUpTable = 0;
  this->SignLookUpTableMask = 0;
  this->MaximumSignLookUp = 0;
  this->NbrParticleLookUpTable = 0;

  this->RescalingFactors = 0;
  this->NbrStateInOrbit = 0;

  this->ReorderingSign = 0;

  this->LargeHilbertSpaceDimension = 0l;
}  


// basic constructor
// 
// nbrFermions= number of fermions
// totalSpin = twice the total spin value
// maxMomentum = momentum maximum value for a fermion
// xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
// yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)

FermionOnTorusWithSpinAndMagneticTranslationsLong::FermionOnTorusWithSpinAndMagneticTranslationsLong (int nbrFermions, int totalSpin, int maxMomentum, 
											      int xMomentum, int yMomentum)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;  
  this->MaxMomentum = maxMomentum;  

  this->NbrMomentum = this->MaxMomentum + 1;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = FindGCD(this->NbrFermions, this->MaxMomentum);
  cout << "MomentumModulo=" << MomentumModulo<<endl;
  this->XMomentum = xMomentum % this->MomentumModulo;
  this->YMomentum = yMomentum % this->MaxMomentum;

  this->StateShift = 2 * (this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrFermions * this->StateShift / 2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((ULONGLONG) 0x1ul);
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= ((ULONGLONG) 0x1ul);
    }

  this->MaximumSignLookUp = 16;
  this->GenerateSignLookUpTable();
  this->HilbertSpaceDimension = ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum - 1, 0, this->NbrFermionsUp);
  cout << "intermediate Hilbert space dimension = " << HilbertSpaceDimension << endl;
  this->HilbertSpaceDimension = this->GenerateStates();
  cout << "Hilbert space dimension = "<< HilbertSpaceDimension << endl;
  this->Flag.Initialize();

  if (this->HilbertSpaceDimension !=0)
    this->GenerateLookUpTable(1000000);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrMomentum * sizeof(int);
  //  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
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

// basic constructor without constraint on Sz
// 
// nbrFermions= number of fermions
// maxMomentum = momentum maximum value for a fermion
// xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
// yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)  

FermionOnTorusWithSpinAndMagneticTranslationsLong::FermionOnTorusWithSpinAndMagneticTranslationsLong (int nbrFermions, int maxMomentum, int xMomentum, int yMomentum)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalSpin = 0;
  this->NbrFermionsUp = this->NbrFermions;
  this->NbrFermionsDown = 0;  
  this->MaxMomentum = maxMomentum;  

  this->NbrMomentum = this->MaxMomentum + 1;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = FindGCD(this->NbrFermions, this->MaxMomentum);
  cout << "MomentumModulo=" << MomentumModulo<<endl;
  this->XMomentum = xMomentum % this->MomentumModulo;
  this->YMomentum = yMomentum % this->MaxMomentum;

  this->StateShift = 2*(this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2*this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((ULONGLONG) 0x1ul);
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= ((ULONGLONG) 0x1ul);
    }

  this->MaximumSignLookUp = 16;
  this->GenerateSignLookUpTable();
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum - 1, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  this->HilbertSpaceDimension = this->GenerateStates(true);
  cout << "Hilbert space dimension = "<< HilbertSpaceDimension << endl;
  this->Flag.Initialize();
  if (this->HilbertSpaceDimension != 0)
    this->GenerateLookUpTable(1000000);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrMomentum * sizeof(int);
  //  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
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
// fermions = reference on the hilbert space to copy to copy

FermionOnTorusWithSpinAndMagneticTranslationsLong::FermionOnTorusWithSpinAndMagneticTranslationsLong(const FermionOnTorusWithSpinAndMagneticTranslationsLong& fermions)
{
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->IncNbrFermions = fermions.IncNbrFermions;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnTorusWithSpinAndMagneticTranslationsLong::~FermionOnTorusWithSpinAndMagneticTranslationsLong ()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
      delete[] this->NbrStateInOrbit;      
      delete[] this->ReorderingSign;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;
      delete[] this->SignLookUpTableMask;      
  
      if (this->HilbertSpaceDimension != 0)
	{
	  delete[] this->LookUpTableShift;
	  for (int i = 0; i < this->NbrFermionStates; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;	  
	  
	  if (this->RescalingFactors != 0)
	    {
	      for (int i = 1; i <= this->MaxMomentum ; ++i)
		delete[] this->RescalingFactors[i];
	      delete[] this->RescalingFactors;	  
	    }
	}
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnTorusWithSpinAndMagneticTranslationsLong& FermionOnTorusWithSpinAndMagneticTranslationsLong::operator = (const FermionOnTorusWithSpinAndMagneticTranslationsLong& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnTorusWithSpinAndMagneticTranslationsLong::Clone()
{
  return new FermionOnTorusWithSpinAndMagneticTranslationsLong(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnTorusWithSpinAndMagneticTranslationsLong::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->XMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->YMomentum, this->MomentumModulo);
  List<AbstractQuantumNumber*> TmpList2;
  TmpList2 += new VectorQuantumNumber (TmpList);
  return TmpList2;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnTorusWithSpinAndMagneticTranslationsLong::GetQuantumNumber (int index)
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->XMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->YMomentum, this->MomentumModulo);
  return new VectorQuantumNumber (TmpList);
}

// get momemtum value in the y direction of a given state
//
// index = state index
// return value = state momentum in the y direction

int FermionOnTorusWithSpinAndMagneticTranslationsLong::GetYMomentumValue(int index)
{
  int StateHighestBit = this->StateHighestBit[index];
  ULONGLONG State = this->StateDescription[index];
  int Momentum = 0;
  for (int i = 0; i <= StateHighestBit; ++i)
    {
      Momentum += ((int) ((State >> (2*i) ) & ((ULONGLONG) 0x1ul) + (State >> (2*i+1) ) & ((ULONGLONG) 0x1ul))) * i;
    }
  return (Momentum % this->MomentumModulo);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnTorusWithSpinAndMagneticTranslationsLong::ExtractSubspace (AbstractQuantumNumber& q, 
									       SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != (AbstractQuantumNumber::Vector | AbstractQuantumNumber::PeriodicMomentum))
    return 0;
  if (((VectorQuantumNumber&) q).GetQuantumNumbers().GetNbrElement() != 2)
    return 0;
  if ((this->XMomentum == ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[0]))->GetMomentum()) &&
      (this->YMomentum == ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[1]))->GetMomentum()))
    return this;
  else 
    return 0;
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
int FermionOnTorusWithSpinAndMagneticTranslationsLong::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  ULONGLONG State = this->StateDescription[index];
  n1<<=1;
  n2<<=1;  
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((((ULONGLONG) 0x1ul)) << n1)) == 0) || 
      ((State & ((((ULONGLONG) 0x1ul)) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  ULONGLONG TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 64))  & this->SignLookUpTableMask[n2 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 80))  & this->SignLookUpTableMask[n2 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 96))  & this->SignLookUpTableMask[n2 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 112))  & this->SignLookUpTableMask[n2 + 112]];
#endif
  TmpState &= ~(((ULONGLONG) 0x1ul) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 64))  & this->SignLookUpTableMask[n1 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 80))  & this->SignLookUpTableMask[n1 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 96))  & this->SignLookUpTableMask[n1 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 112))  & this->SignLookUpTableMask[n1 + 112]];
#endif
  TmpState &= ~(((ULONGLONG) 0x1ul) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == ((ULONGLONG) 0x0ul))
      --NewHighestBit;
  m1<<=1;
  m2<<=1;
  if ((TmpState & (((ULONGLONG) 0x1ul) << m2))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewHighestBit)
    {
      NewHighestBit = m2;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80)) & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  if ((TmpState & (((ULONGLONG) 0x1ul) << m1))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    {
      NewHighestBit = m1;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80)) & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  TmpState = this->FindCanonicalForm(TmpState, NewHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, NewHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
  nbrTranslation *= this->StateShift/2;
  return TmpIndex;
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

int FermionOnTorusWithSpinAndMagneticTranslationsLong::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  ULONGLONG State = this->StateDescription[index];
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((((ULONGLONG) 0x1ul)) << n1)) == ((ULONGLONG) 0x0ul)) || 
      ((State & ((((ULONGLONG) 0x1ul)) << n2)) == ((ULONGLONG) 0x0ul)) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  ULONGLONG TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  TmpState &= ~(((ULONGLONG) 0x1ul) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  TmpState &= ~(((ULONGLONG) 0x1ul) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;  
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  if ((TmpState & (((ULONGLONG) 0x1ul) << m2))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewHighestBit)
    {
      NewHighestBit = m2;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80)) & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  if ((TmpState & (((ULONGLONG) 0x1ul) << m1))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    {
      NewHighestBit = m1;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80)) & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  TmpState = this->FindCanonicalForm(TmpState, NewHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, NewHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
  nbrTranslation *= this->StateShift/2;
  return TmpIndex;
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
int FermionOnTorusWithSpinAndMagneticTranslationsLong::AduAduAuAuV (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  ULONGLONG State = this->StateDescription[index];
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((((ULONGLONG) 0x1ul)) << n1)) == 0) || 
      ((State & ((((ULONGLONG) 0x1ul)) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  ULONGLONG TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  TmpState &= ~(((ULONGLONG) 0x1ul) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  TmpState &= ~(((ULONGLONG) 0x1ul) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == ((ULONGLONG) 0x0ul))
      --NewHighestBit;
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  if ((TmpState & (((ULONGLONG) 0x1ul) << m2))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewHighestBit)
    {
      NewHighestBit = m2;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80)) & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  if ((TmpState & (((ULONGLONG) 0x1ul) << m1))!= ((ULONGLONG) 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    {
      NewHighestBit = m1;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80)) & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  TmpState = this->FindCanonicalForm(TmpState, NewHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, NewHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
  nbrTranslation *= this->StateShift/2;
  return TmpIndex;
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
int FermionOnTorusWithSpinAndMagneticTranslationsLong::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  ULONGLONG State = this->StateDescription[index];
  n1<<=1;
  n2<<=1;
  n2+=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((((ULONGLONG) 0x1ul)) << n1)) == 0) || ((State & ((((ULONGLONG) 0x1ul)) << n2)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  ULONGLONG TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  TmpState &= ~(((ULONGLONG) 0x1ul) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  TmpState &= ~(((ULONGLONG) 0x1ul) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == ((ULONGLONG) 0x0ul))
      --NewHighestBit;
  m1<<=1;
  m2<<=1;
  m2+=1;
  if ((TmpState & (((ULONGLONG) 0x1ul) << m2))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewHighestBit)
    {
      NewHighestBit = m2;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80)) & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  if ((TmpState & (((ULONGLONG) 0x1ul) << m1))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    {
      NewHighestBit = m1;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80)) & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  TmpState = this->FindCanonicalForm(TmpState, NewHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, NewHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
  nbrTranslation *= this->StateShift/2;
  return TmpIndex;
}


// apply a^+_m_d a_m_d operator to a given state (only spin down)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
double FermionOnTorusWithSpinAndMagneticTranslationsLong::AddAd (int index, int m)
{
  if (this->StateHighestBit[index] < 2*m)
    return 0.0;
  if ((this->StateDescription[index] & (((ULONGLONG) 0x1ul) << (m<<1))) == ((ULONGLONG) 0x0ul))
    return 0.0;
  else
    return 1.0;
}
  

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
double FermionOnTorusWithSpinAndMagneticTranslationsLong::AduAu (int index, int m)
{
  if (this->StateHighestBit[index] < 2*m+1)
    return 0.0;
  if ((this->StateDescription[index] & (((ULONGLONG) 0x1ul) << (2*m+1))) == ((ULONGLONG) 0x0ul))
    return 0.0;
  else
    return 1.0;
}

// apply a^+_m_u a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithSpinAndMagneticTranslationsLong::AduAu (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  ULONGLONG State = this->StateDescription[index];
  m = (m << 1) + 1;
  n = (n << 1) + 1;
  if ((n > StateHighestBit) || ((State & (((ULONGLONG) 0x1ul) << n)) == ((ULONGLONG) 0x0ul)) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(State >> (n + 64)) & this->SignLookUpTableMask[n + 64]];
  coefficient *= this->SignLookUpTable[(State >> (n + 80)) & this->SignLookUpTableMask[n + 80]];
  coefficient *= this->SignLookUpTable[(State >> (n + 96)) & this->SignLookUpTableMask[n + 96]];
  coefficient *= this->SignLookUpTable[(State >> (n + 112)) & this->SignLookUpTableMask[n + 112]];
#endif
  State &= ~(((ULONGLONG) 0x1ul) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == ((ULONGLONG) 0x0ul))
      --NewLargestBit;

  if ((State & (((ULONGLONG) 0x1ul) << m))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(State >> (m + 64)) & this->SignLookUpTableMask[m + 64]];
      coefficient *= this->SignLookUpTable[(State >> (m + 80)) & this->SignLookUpTableMask[m + 80]];
      coefficient *= this->SignLookUpTable[(State >> (m + 96)) & this->SignLookUpTableMask[m + 96]];
      coefficient *= this->SignLookUpTable[(State >> (m + 112)) & this->SignLookUpTableMask[m + 112]];
#endif
    }
  State |= (((ULONGLONG) 0x1ul) << m);
  State = this->FindCanonicalForm(State, NewLargestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(State, NewLargestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(State, NewLargestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift / 2;
    }
  return TmpIndex;
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithSpinAndMagneticTranslationsLong::AddAd (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  ULONGLONG State = this->StateDescription[index];
  m <<= 1;
  n <<= 1;
  if ((n > StateHighestBit) || ((State & (((ULONGLONG) 0x1ul) << n)) == ((ULONGLONG) 0x0ul)) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(State >> (n + 64)) & this->SignLookUpTableMask[n + 64]];
  coefficient *= this->SignLookUpTable[(State >> (n + 80)) & this->SignLookUpTableMask[n + 80]];
  coefficient *= this->SignLookUpTable[(State >> (n + 96)) & this->SignLookUpTableMask[n + 96]];
  coefficient *= this->SignLookUpTable[(State >> (n + 112)) & this->SignLookUpTableMask[n + 112]];
#endif
  State &= ~(((ULONGLONG) 0x1ul) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == ((ULONGLONG) 0x0ul))
      --NewLargestBit;

  if ((State & (((ULONGLONG) 0x1ul) << m))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(State >> (m + 64)) & this->SignLookUpTableMask[m + 64]];
      coefficient *= this->SignLookUpTable[(State >> (m + 80)) & this->SignLookUpTableMask[m + 80]];
      coefficient *= this->SignLookUpTable[(State >> (m + 96)) & this->SignLookUpTableMask[m + 96]];
      coefficient *= this->SignLookUpTable[(State >> (m + 112)) & this->SignLookUpTableMask[m + 112]];
#endif
    }
  State |= (((ULONGLONG) 0x1ul) << m);
  State = this->FindCanonicalForm(State, NewLargestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(State, NewLargestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(State, NewLargestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift / 2;
    }
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

int FermionOnTorusWithSpinAndMagneticTranslationsLong::AduAd (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  ULONGLONG State = this->StateDescription[index];
  m = (m << 1) + 1;
  n <<= 1;
  if ((n > StateHighestBit) || ((State & (((ULONGLONG) 0x1ul) << n)) == ((ULONGLONG) 0x0ul)) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(State >> (n + 64)) & this->SignLookUpTableMask[n + 64]];
  coefficient *= this->SignLookUpTable[(State >> (n + 80)) & this->SignLookUpTableMask[n + 80]];
  coefficient *= this->SignLookUpTable[(State >> (n + 96)) & this->SignLookUpTableMask[n + 96]];
  coefficient *= this->SignLookUpTable[(State >> (n + 112)) & this->SignLookUpTableMask[n + 112]];
#endif
  State &= ~(((ULONGLONG) 0x1ul) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == ((ULONGLONG) 0x0ul))
      --NewLargestBit;

  if ((State & (((ULONGLONG) 0x1ul) << m))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(State >> (m + 64)) & this->SignLookUpTableMask[m + 64]];
      coefficient *= this->SignLookUpTable[(State >> (m + 80)) & this->SignLookUpTableMask[m + 80]];
      coefficient *= this->SignLookUpTable[(State >> (m + 96)) & this->SignLookUpTableMask[m + 96]];
      coefficient *= this->SignLookUpTable[(State >> (m + 112)) & this->SignLookUpTableMask[m + 112]];
#endif
    }
  State |= (((ULONGLONG) 0x1ul) << m);
  State = this->FindCanonicalForm(State, NewLargestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(State, NewLargestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(State, NewLargestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift / 2;
    }
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

int FermionOnTorusWithSpinAndMagneticTranslationsLong::AddAu (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  ULONGLONG State = this->StateDescription[index];
  m <<= 1;
  n = (n << 1) + 1;  
  if ((n > StateHighestBit) || ((State & (((ULONGLONG) 0x1ul) << n)) == ((ULONGLONG) 0x0ul)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(State >> (n + 64)) & this->SignLookUpTableMask[n + 64]];
  coefficient *= this->SignLookUpTable[(State >> (n + 80)) & this->SignLookUpTableMask[n + 80]];
  coefficient *= this->SignLookUpTable[(State >> (n + 96)) & this->SignLookUpTableMask[n + 96]];
  coefficient *= this->SignLookUpTable[(State >> (n + 112)) & this->SignLookUpTableMask[n + 112]];
#endif
  State &= ~(((ULONGLONG) 0x1ul) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == ((ULONGLONG) 0x0ul))
      --NewLargestBit;

  if ((State & (((ULONGLONG) 0x1ul) << m))!= ((ULONGLONG) 0x0ul))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(State >> (m + 64)) & this->SignLookUpTableMask[m + 64]];
      coefficient *= this->SignLookUpTable[(State >> (m + 80)) & this->SignLookUpTableMask[m + 80]];
      coefficient *= this->SignLookUpTable[(State >> (m + 96)) & this->SignLookUpTableMask[m + 96]];
      coefficient *= this->SignLookUpTable[(State >> (m + 112)) & this->SignLookUpTableMask[m + 112]];
#endif
    }
  State |= (((ULONGLONG) 0x1ul) << m);
  State = this->FindCanonicalForm(State, NewLargestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(State, NewLargestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(State, NewLargestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & ((ULONGLONG) 0x1ul))));
      nbrTranslation *= this->StateShift / 2;
    }
  return TmpIndex;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinAndMagneticTranslationsLong::AuAu (int index, int n1, int n2)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1; 
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((((ULONGLONG) 0x1ul)) << n1)) == ((ULONGLONG) 0x0ul)) || 
      ((ProdATemporaryState & ((((ULONGLONG) 0x1ul)) << n2)) == ((ULONGLONG) 0x0ul)) || (n1 == n2))
    return 0.0;  
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == ((ULONGLONG) 0x0ul))
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == ((ULONGLONG) 0x0ul))
      --ProdAHighestBit;
  this->ProdATemporaryNbrStateInOrbit =  this->NbrStateInOrbit[index];
  return coefficient;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinAndMagneticTranslationsLong::AuAuV (int index, int n1, int n2)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((((ULONGLONG) 0x1ul)) << n1)) == ((ULONGLONG) 0x0ul)) || 
      ((ProdATemporaryState & ((((ULONGLONG) 0x1ul)) << n2)) == ((ULONGLONG) 0x0ul)) || (n1 == n2))
    return 0.0;  
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == ((ULONGLONG) 0x0ul))
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == ((ULONGLONG) 0x0ul))
      --ProdAHighestBit;
  this->ProdATemporaryNbrStateInOrbit =  this->NbrStateInOrbit[index];
  return coefficient;
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinAndMagneticTranslationsLong::AdAd (int index, int n1, int n2)
{  
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n2<<=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((((ULONGLONG) 0x1ul)) << n1)) == ((ULONGLONG) 0x0ul)) || 
      ((ProdATemporaryState & ((((ULONGLONG) 0x1ul)) << n2)) == ((ULONGLONG) 0x0ul)) || (n1 == n2))
    return 0.0;
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == ((ULONGLONG) 0x0ul))
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == ((ULONGLONG) 0x0ul))
      --ProdAHighestBit;
  this->ProdATemporaryNbrStateInOrbit =  this->NbrStateInOrbit[index];
  return coefficient;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinAndMagneticTranslationsLong::AuAd (int index, int n1, int n2)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n1+=1;
  n2<<=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((((ULONGLONG) 0x1ul)) << n1)) == ((ULONGLONG) 0x0ul)) || 
      ((ProdATemporaryState & ((((ULONGLONG) 0x1ul)) << n2)) == ((ULONGLONG) 0x0ul)) || (n1 == n2))
    return 0.0;
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == ((ULONGLONG) 0x0ul))
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == ((ULONGLONG) 0x0ul))
      --ProdAHighestBit;
  this->ProdATemporaryNbrStateInOrbit =  this->NbrStateInOrbit[index];
  return coefficient;
}

// apply a_n1_d a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinAndMagneticTranslationsLong::AdAu (int index, int n1, int n2)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n2<<=1;
  n2+=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((((ULONGLONG) 0x1ul)) << n1)) == ((ULONGLONG) 0x0ul)) || 
      ((ProdATemporaryState & ((((ULONGLONG) 0x1ul)) << n2)) == ((ULONGLONG) 0x0ul)) || (n1 == n2))
    return 0.0;
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 64)) & this->SignLookUpTableMask[n2 + 64]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 80)) & this->SignLookUpTableMask[n2 + 80]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 96)) & this->SignLookUpTableMask[n2 + 96]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 112)) & this->SignLookUpTableMask[n2 + 112]];
#endif
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == ((ULONGLONG) 0x0ul))
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#ifdef __128_BIT_LONGLONG__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 64)) & this->SignLookUpTableMask[n1 + 64]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 80)) & this->SignLookUpTableMask[n1 + 80]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 96)) & this->SignLookUpTableMask[n1 + 96]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 112)) & this->SignLookUpTableMask[n1 + 112]];
#endif
  ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == ((ULONGLONG) 0x0ul))
      --ProdAHighestBit;
  this->ProdATemporaryNbrStateInOrbit =  this->NbrStateInOrbit[index];
  return coefficient;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnTorusWithSpinAndMagneticTranslationsLong::AduAdu (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  int TmpAHighestBit = ProdAHighestBit;
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  if (((TmpState & (((ULONGLONG) 0x1ul) << m1)) != ((ULONGLONG) 0x0ul)) || ((TmpState & (((ULONGLONG) 0x1ul) << m2)) != ((ULONGLONG) 0x0ul)) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  if (m2 > TmpAHighestBit)
    TmpAHighestBit = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80)) & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  
  if (m1 > TmpAHighestBit)
    TmpAHighestBit = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80)) & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  TmpState = this->FindCanonicalForm(TmpState, TmpAHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, TmpAHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, TmpAHighestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[ProdAIndex]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift/2;
    }
  return TmpIndex;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnTorusWithSpinAndMagneticTranslationsLong::AduAduV (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  int TmpAHighestBit = ProdAHighestBit;
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  if (((TmpState & (((ULONGLONG) 0x1ul) << m1)) != ((ULONGLONG) 0x0ul)) || ((TmpState & (((ULONGLONG) 0x1ul) << m2)) != ((ULONGLONG) 0x0ul)) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  if (m2 > TmpAHighestBit)
    TmpAHighestBit = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80)) & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  if (m1 > TmpAHighestBit)
    TmpAHighestBit = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80)) & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  TmpState = this->FindCanonicalForm(TmpState, TmpAHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, TmpAHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, TmpAHighestBit);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[ProdAIndex]][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
  nbrTranslation *= this->StateShift/2;
  return TmpIndex;
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnTorusWithSpinAndMagneticTranslationsLong::AddAdd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  int TmpAHighestBit = ProdAHighestBit;
  m1<<=1;
  m2<<=1;
  if (((TmpState & (((ULONGLONG) 0x1ul) << m1)) != ((ULONGLONG) 0x0ul)) || ((TmpState & (((ULONGLONG) 0x1ul) << m2)) != ((ULONGLONG) 0x0ul)) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  if (m2 > TmpAHighestBit)
    TmpAHighestBit = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80)) & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  
  if (m1 > TmpAHighestBit)
    TmpAHighestBit = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80)) & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  TmpState = this->FindCanonicalForm(TmpState, TmpAHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, TmpAHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, TmpAHighestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[ProdAIndex]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift/2;
    }
  return TmpIndex;
}

  

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnTorusWithSpinAndMagneticTranslationsLong::AduAdd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  ULONGLONG TmpState = this->ProdATemporaryState;
  int TmpAHighestBit = ProdAHighestBit;
  m1<<=1;
  m1+=1;
  m2<<=1;
  if (((TmpState & (((ULONGLONG) 0x1ul) << m1)) != ((ULONGLONG) 0x0ul)) || ((TmpState & (((ULONGLONG) 0x1ul) << m2)) != ((ULONGLONG) 0x0ul)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  if (m2 > TmpAHighestBit)
    TmpAHighestBit = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 64)) & this->SignLookUpTableMask[m2 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 80)) & this->SignLookUpTableMask[m2 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 96)) & this->SignLookUpTableMask[m2 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 112)) & this->SignLookUpTableMask[m2 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m2);
  
  if (m1 > TmpAHighestBit)
    TmpAHighestBit = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#ifdef __128_BIT_LONGLONG__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 64)) & this->SignLookUpTableMask[m1 + 64]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 80)) & this->SignLookUpTableMask[m1 + 80]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 96)) & this->SignLookUpTableMask[m1 + 96]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 112)) & this->SignLookUpTableMask[m1 + 112]];
#endif
    }
  TmpState |= (((ULONGLONG) 0x1ul) << m1);
  TmpState = this->FindCanonicalForm(TmpState, TmpAHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, TmpAHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, TmpAHighestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[ProdAIndex]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift/2;
    }
  return TmpIndex;
}


// apply a_n_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator (spin up)
// return value =  multiplicative factor 

double FermionOnTorusWithSpinAndMagneticTranslationsLong::Au (int index, int n)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
//   this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  n <<= 1;
  ++n;
  
  if ((n >  this->ProdAHighestBit) || ((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n)) == ((ULONGLONG) 0x0ul)))
    {
      return 0.0;
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n) & this->SignLookUpTableMask[n]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 64)) & this->SignLookUpTableMask[n + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 80)) & this->SignLookUpTableMask[n + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 96)) & this->SignLookUpTableMask[n + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 112)) & this->SignLookUpTableMask[n + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n);
  if (this->ProdATemporaryState == ((ULONGLONG) 0x0ul))
    {
      this->ProdAHighestBit = 0;
    }
  else
    {
      if (this->ProdAHighestBit == n)
	while ((this->ProdATemporaryState >> this->ProdAHighestBit) == ((ULONGLONG) 0x0ul))
	  --this->ProdAHighestBit;
    }
  return Coefficient;
}

// apply a_n_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next Ad call
//
// index = index of the state on which the operator has to be applied
// n = first index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnTorusWithSpinAndMagneticTranslationsLong::Ad (int index, int n)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  n <<= 1;
  
  if ((n >  this->ProdAHighestBit) || ((this->ProdATemporaryState & (((ULONGLONG) 0x1ul) << n)) == ((ULONGLONG) 0x0ul)))
    {
      return 0.0;
    }
  this->ProdATemporaryNbrStateInOrbit = this->NbrStateInOrbit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n) & this->SignLookUpTableMask[n]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#ifdef __128_BIT_LONGLONG__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 64)) & this->SignLookUpTableMask[n + 64]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 80)) & this->SignLookUpTableMask[n + 80]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 96)) & this->SignLookUpTableMask[n + 96]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n + 112)) & this->SignLookUpTableMask[n + 112]];
#endif
  this->ProdATemporaryState &= ~(((ULONGLONG) 0x1ul) << n);
  if (this->ProdATemporaryState == ((ULONGLONG) 0x0ul))
    {
      this->ProdAHighestBit = 0;
    }
  else
    {
      if (this->ProdAHighestBit == n)
	while ((this->ProdATemporaryState >> this->ProdAHighestBit) == ((ULONGLONG) 0x0ul))
	  --this->ProdAHighestBit;
    }
  return Coefficient;
}


// convert a state defined in the Ky basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

// ComplexVector FermionOnTorusWithSpinAndMagneticTranslationsLong::ConvertToKxKyBasis(ComplexVector& state, ParticleOnSphere* space)  
// {
//   FermionOnTorusWithSpinLong* TmpSpace = (FermionOnTorusWithSpinLong*) space;
//   ComplexVector TmpVector (this->HilbertSpaceDimension, true);
//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//     {
//       ULONGLONG TmpState = this->StateDescription[i];
//       int TmpMaxMomentum = this->StateHighestBit[i];
//       int Pos = TmpSpace->FindStateIndex(TmpState, TmpMaxMomentum);
//       if (Pos < TmpSpace->HilbertSpaceDimension)
// 	{
// 	  TmpVector[i] =  state[Pos] * sqrt((double) this->NbrStateInOrbit[i]);
// 	}
//     }
//   return TmpVector;
// }

// convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

// ComplexVector FermionOnTorusWithSpinAndMagneticTranslationsLong::ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space)
// {
//   FermionOnTorusWithSpinLong* TmpSpace = (FermionOnTorusWithSpinLong*) space;
//   ComplexVector TmpVector (TmpSpace->HilbertSpaceDimension, true);
//   Complex* FourrierCoefficients = new Complex [this->MomentumModulo];
//   for (int i = 0; i < this->MomentumModulo; ++i)
//     FourrierCoefficients[i] = Phase (-2.0 * M_PI * ((double) (i * this->XMomentum)) / ((double) this->MomentumModulo));
//   for (int i = 0; i < TmpSpace->HilbertSpaceDimension; ++i)
//     {
//       ULONGLONG TmpState = TmpSpace->StateDescription[i];
//       int NbrTranslation = 0;
//       int TmpMaxMomentum = TmpSpace->StateHighestBit[i];
//       TmpState = this->FindCanonicalFormAndTestXMomentumConstraint(TmpState, TmpMaxMomentum, NbrTranslation);
//       if (NbrTranslation >= 0)
// 	{
// 	  int Pos = this->FindStateIndex(TmpState, TmpMaxMomentum);
// 	  if (Pos < this->HilbertSpaceDimension)
// 	    {
// 	      TmpVector[i] =  (state[Pos] * (1.0 - (2.0 * ((double) ((this->ReorderingSign[Pos] >> NbrTranslation) & 0x1ul)))) * 
// 			       FourrierCoefficients[NbrTranslation] / sqrt((double) this->NbrStateInOrbit[Pos]));
// 	    }
// 	}
//     }
//   delete[] FourrierCoefficients;
//   return TmpVector;
// }


// find canonical form of a state description
//
// stateDescription = unsigned integer describing the state
// stateHighestBit = reference on the highest non-zero bit in the canonical representation
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description

ULONGLONG FermionOnTorusWithSpinAndMagneticTranslationsLong::FindCanonicalForm(ULONGLONG stateDescription, int& stateHighestBit, int& nbrTranslation)
{
  nbrTranslation = 0;
  ULONGLONG CanonicalState = stateDescription;
  ULONGLONG stateDescriptionReference = stateDescription;
  int index = 1;
  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);  
  while (stateDescriptionReference != stateDescription)
    {
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = index;
	}
      stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      ++index;
    }
  stateHighestBit = 2*this->MaxMomentum-1;
  stateDescription = ((ULONGLONG) 0x1ul) << stateHighestBit;
  while ((CanonicalState & stateDescription) == 0)
    {
      --stateHighestBit;
      stateDescription >>= 1;
    }
  if (nbrTranslation != 0)
    nbrTranslation = index - nbrTranslation;
  return CanonicalState;
}

// find how many translations on the x direction are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

int  FermionOnTorusWithSpinAndMagneticTranslationsLong::FindNumberXTranslation(ULONGLONG stateDescription)
{
  ULONGLONG TmpState = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
  int index = 1;  
  while (TmpState != stateDescription)
    {
      TmpState = (TmpState >> this->StateShift) | ((TmpState & this->MomentumMask) << this->ComplementaryStateShift);
      ++index;
    }
  return index;
}

// test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescription = unsigned integer describing the state
// maxMomentum = maximum momentum value that can be reached by a fermion in the stateDescription state
// return value = true if the state satisfy the x momentum constraint

bool FermionOnTorusWithSpinAndMagneticTranslationsLong::TestXMomentumConstraint(ULONGLONG stateDescription, int maxMomentum)
{
  if (this->NbrFermions & 1)
    {
      if (this->XMomentum == 0)
	return true;
      ULONGLONG TmpState = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      int index = 1;  
      while (TmpState != stateDescription)
	{
	  TmpState = (TmpState >> this->StateShift) | ((TmpState & this->MomentumMask) << this->ComplementaryStateShift);
	  ++index;
	}
      if (((this->XMomentum * index) % this->MomentumModulo) == 0)
	return true;
      else
	return false;
    }
  else
    {
      ULONGLONG TmpState2 = stateDescription & this->MomentumMask;
      ULONGLONG TmpState = (stateDescription >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
      int TmpSignature = 0;
      int index = 1;  
#ifdef __128_BIT_LONGLONG__
      TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & ((ULONGLONG) 0xfffful)] 
		       + this->NbrParticleLookUpTable[(TmpState2 >> 16) & ((ULONGLONG) 0xfffful)]
		       + this->NbrParticleLookUpTable[(TmpState2 >> 32) & ((ULONGLONG) 0xfffful)] 
		       + this->NbrParticleLookUpTable[(TmpState2 >> 48) & ((ULONGLONG) 0xfffful)]
		       + this->NbrParticleLookUpTable[(TmpState2 >> 64) & ((ULONGLONG) 0xfffful)]
		       + this->NbrParticleLookUpTable[(TmpState2 >> 80) & ((ULONGLONG) 0xfffful)]
		       + this->NbrParticleLookUpTable[(TmpState2 >> 96) & ((ULONGLONG) 0xfffful)] 
		       + this->NbrParticleLookUpTable[(TmpState2 >> 112) & ((ULONGLONG) 0xfffful)]) & 1;
#else
      TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & ((ULONGLONG) 0xfffful)] 
		       + this->NbrParticleLookUpTable[(TmpState2 >> 16) & ((ULONGLONG) 0xfffful)]
		       + this->NbrParticleLookUpTable[(TmpState2 >> 32) & ((ULONGLONG) 0xfffful)] 
		       + this->NbrParticleLookUpTable[(TmpState2 >> 48) & ((ULONGLONG) 0xfffful)]) & 1;
#endif

      while (TmpState != stateDescription)
	{
	  TmpState2 = TmpState & this->MomentumMask;
	  TmpState = (TmpState >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
#ifdef __128_BIT_LONGLONG__
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & ((ULONGLONG) 0xfffful)] 
			   + this->NbrParticleLookUpTable[(TmpState2 >> 16) & ((ULONGLONG) 0xfffful)]
			   + this->NbrParticleLookUpTable[(TmpState2 >> 32) & ((ULONGLONG) 0xfffful)] 
			   + this->NbrParticleLookUpTable[(TmpState2 >> 48) & ((ULONGLONG) 0xfffful)]) & 1;
#else
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & ((ULONGLONG) 0xfffful)] 
			   + this->NbrParticleLookUpTable[(TmpState2 >> 16) & ((ULONGLONG) 0xfffful)]
			   + this->NbrParticleLookUpTable[(TmpState2 >> 32) & ((ULONGLONG) 0xfffful)] 
			   + this->NbrParticleLookUpTable[(TmpState2 >> 48) & ((ULONGLONG) 0xfffful)]
			   + this->NbrParticleLookUpTable[(TmpState2 >> 64) & ((ULONGLONG) 0xfffful)]
			   + this->NbrParticleLookUpTable[(TmpState2 >> 80) & ((ULONGLONG) 0xfffful)]
			   + this->NbrParticleLookUpTable[(TmpState2 >> 96) & ((ULONGLONG) 0xfffful)] 
			   + this->NbrParticleLookUpTable[(TmpState2 >> 112) & ((ULONGLONG) 0xfffful)]) & 1;
#endif
	  ++index;
	}
      if ((((this->XMomentum * index) - ((this->MomentumModulo * TmpSignature) >> 1)) % this->MomentumModulo) == 0)
	return true;
      else
	return false;
    }
}

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescription = unsigned integer describing the state
// stateHighestBit = reference on the highest non-zero bit in the canonical representation
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description and -1 in nbrTranslation if the state does not fit the x momentum constraint

ULONGLONG FermionOnTorusWithSpinAndMagneticTranslationsLong::FindCanonicalFormAndTestXMomentumConstraint(ULONGLONG stateDescription, int& stateHighestBit, int& nbrTranslation)
{
  nbrTranslation = 0;
  ULONGLONG CanonicalState = stateDescription;
  ULONGLONG stateDescriptionReference = stateDescription;
  int index = 1;  
  if (this->NbrFermions & 1)
    {
      stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      while (stateDescriptionReference != stateDescription)
	{
	  if (stateDescription < CanonicalState)
	    {
	      CanonicalState = stateDescription;
	      nbrTranslation = index;
	    }
	  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
	  ++index;
	}
      if (((this->XMomentum * index) % this->MomentumModulo) != 0)
	{
	  nbrTranslation = -1;
	  return CanonicalState;
	}
    }
  else
    {
      ULONGLONG TmpState = stateDescription & this->MomentumMask;
      stateDescription = (stateDescription >> this->StateShift) | (TmpState << this->ComplementaryStateShift);
      int TmpSignature = 0;
#ifdef __128_BIT_LONGLONG__
      TmpSignature += (this->NbrParticleLookUpTable[TmpState & ((ULONGLONG) 0xfffful)] 
		       + this->NbrParticleLookUpTable[(TmpState >> 16) & ((ULONGLONG) 0xfffful)]
		       + this->NbrParticleLookUpTable[(TmpState >> 32) & ((ULONGLONG) 0xfffful)] 
		       + this->NbrParticleLookUpTable[(TmpState >> 48) & ((ULONGLONG) 0xfffful)]) & 1;
#else
      TmpSignature += (this->NbrParticleLookUpTable[TmpState & ((ULONGLONG) 0xfffful)] 
		       + this->NbrParticleLookUpTable[(TmpState >> 16) & ((ULONGLONG) 0xfffful)]
		       + this->NbrParticleLookUpTable[(TmpState >> 32) & ((ULONGLONG) 0xfffful)] 
		       + this->NbrParticleLookUpTable[(TmpState >> 48) & ((ULONGLONG) 0xfffful)]
		       + this->NbrParticleLookUpTable[(TmpState >> 64) & ((ULONGLONG) 0xfffful)]
		       + this->NbrParticleLookUpTable[(TmpState >> 80) & ((ULONGLONG) 0xfffful)]
		       + this->NbrParticleLookUpTable[(TmpState >> 96) & ((ULONGLONG) 0xfffful)] 
		       + this->NbrParticleLookUpTable[(TmpState >> 112) & ((ULONGLONG) 0xfffful)]) & 1;
#endif

      while (stateDescription != stateDescriptionReference)
	{
	  if (stateDescription < CanonicalState)
	    {
	      CanonicalState = stateDescription;
	      nbrTranslation = index;
	    }
	  TmpState = stateDescription & this->MomentumMask;
	  stateDescription = (stateDescription >> this->StateShift) | (TmpState << this->ComplementaryStateShift);
#ifdef __128_BIT_LONGLONG__
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState & ((ULONGLONG) 0xfffful)] 
			   + this->NbrParticleLookUpTable[(TmpState >> 16) & ((ULONGLONG) 0xfffful)]
			   + this->NbrParticleLookUpTable[(TmpState >> 32) & ((ULONGLONG) 0xfffful)] 
			   + this->NbrParticleLookUpTable[(TmpState >> 48) & ((ULONGLONG) 0xfffful)]) & 1;
#else
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState & ((ULONGLONG) 0xfffful)] 
			   + this->NbrParticleLookUpTable[(TmpState >> 16) & ((ULONGLONG) 0xfffful)]
			   + this->NbrParticleLookUpTable[(TmpState >> 32) & ((ULONGLONG) 0xfffful)] 
			   + this->NbrParticleLookUpTable[(TmpState >> 48) & ((ULONGLONG) 0xfffful)]
			   + this->NbrParticleLookUpTable[(TmpState >> 64) & ((ULONGLONG) 0xfffful)]
			   + this->NbrParticleLookUpTable[(TmpState >> 80) & ((ULONGLONG) 0xfffful)]
			   + this->NbrParticleLookUpTable[(TmpState >> 96) & ((ULONGLONG) 0xfffful)] 
			   + this->NbrParticleLookUpTable[(TmpState >> 112) & ((ULONGLONG) 0xfffful)]) & 1;
#endif
	  ++index;
	}
      if ((((this->XMomentum * index) - ((this->MomentumModulo * TmpSignature) >> 1)) % this->MomentumModulo) != 0)
	{
	  nbrTranslation = -1;
	  return CanonicalState;
	}
    }
  if (nbrTranslation != 0)
    {
      stateHighestBit = 2*this->MaxMomentum-1;
      stateDescription = ((ULONGLONG) 0x1ul) << stateHighestBit;
      while ((CanonicalState & stateDescription) == ((ULONGLONG) 0x0ul))      
	{
	  --stateHighestBit;
	  stateDescription >>= 1;
	}
      nbrTranslation = index - nbrTranslation;
    }
  return CanonicalState;
}

// find state index
//
// stateDescription = unsigned longeger describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnTorusWithSpinAndMagneticTranslationsLong::FindStateIndex(ULONGLONG stateDescription, int maxMomentum)
{
  long PosMax = stateDescription >> this->LookUpTableShift[maxMomentum];
  long PosMin = this->LookUpTable[maxMomentum][PosMax];
  PosMax = this->LookUpTable[maxMomentum][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  ULONGLONG CurrentState = this->StateDescription[PosMid];
  while ((PosMin != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    return PosMax;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnTorusWithSpinAndMagneticTranslationsLong::PrintState (ostream& Str, int state)
{
  ULONGLONG TmpState = this->StateDescription[state];
  unsigned long Tmp;
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      Tmp = (unsigned long) ((TmpState >> (i << 1)) & ((ULONGLONG) 0x3ul));
      switch (Tmp)
	{
	case 0x1ul:
	  Str << "d ";
	  break;
	case 0x2ul:
	  Str << "u ";
	  break;
	case 0x3ul:
	  Str << "X ";
	  break;
	default:
	  Str << "0 ";
	  break;
	}
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// fullSzFlag = if true, does not apply the Sz contraint
// fullKyFlag = if true, does not apply the Ky contraint
// return value = hilbert space dimension

int FermionOnTorusWithSpinAndMagneticTranslationsLong::GenerateStates(bool fullSzFlag, bool fullKyFlag)
{
  int maxHighestBit = 2*(this->MaxMomentum + 1);
  this->StateDescription = new ULONGLONG [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  if (fullKyFlag == false)
    {
      if (fullSzFlag == false)
	this->HilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, this->MaxMomentum - 1, 0, this->NbrFermionsUp, 0l);
      else
	this->HilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, 2 * this->MaxMomentum - 1, 2 * this->MaxMomentum - 1, 0l, 0);  
    }
  else
    {
      if (fullSzFlag == false)
      {
	this->HilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, this->MaxMomentum - 1, this->NbrFermionsUp, 0l);
      }
      else
      {
	this->HilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, this->MaxMomentum - 1, 0l);
      }
      int TmpMaxHighestBit = 2*(this->MaxMomentum + 1);
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  while (((this->StateDescription[i] >> TmpMaxHighestBit) & ((ULONGLONG) 0x1ul)) == ((ULONGLONG) 0x0ul))
	    --TmpMaxHighestBit;
	  this->StateHighestBit[i] = TmpMaxHighestBit;
	}
    }
  int* TmpNbrStateDescription = new int [maxHighestBit];  
  for (int i = 0; i < maxHighestBit; ++i)
    {
      TmpNbrStateDescription[i] = 0;
    }
  int NbrTranslation;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->StateDescription[i] = this->FindCanonicalForm(this->StateDescription[i], this->StateHighestBit[i], NbrTranslation);
      ++TmpNbrStateDescription[this->StateHighestBit[i]];
    }  
  ULONGLONG** TmpStateDescription = new ULONGLONG* [maxHighestBit];  
  bool** TmpStateDescriptionFlag = new bool* [maxHighestBit];  
  for (int i = 0; i < maxHighestBit; ++i)
    {
      if (TmpNbrStateDescription[i] != 0)
	{
	  TmpStateDescription[i] = new ULONGLONG [TmpNbrStateDescription[i]];
	  TmpStateDescriptionFlag[i] = new bool [TmpNbrStateDescription[i]];
	}
      else
	{
	  TmpStateDescription[i] = 0;
	  TmpStateDescriptionFlag[i] = 0;
	}
      TmpNbrStateDescription[i] = 0;
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    TmpStateDescription[this->StateHighestBit[i]][TmpNbrStateDescription[this->StateHighestBit[i]]++] = this->StateDescription[i];      
  int TmpHilbertSpaceDimension = 0;
  for (int i = 0; i < maxHighestBit; ++i) 
    {
      int CurrentNbrState = TmpNbrStateDescription[i];
      if (CurrentNbrState > 0)
	{
	  ULONGLONG* TmpStateArray = TmpStateDescription[i];
	  bool* TmpStateArrayFlag = TmpStateDescriptionFlag[i];
	  SortArrayUpOrdering(TmpStateArray, CurrentNbrState);
	  int TmpNbr = CurrentNbrState;
	  if (this->TestXMomentumConstraint(TmpStateArray[0], i) == true)
	    {
	      TmpStateArrayFlag[0] = true;
	    }
	  else
	    {
	      TmpStateArrayFlag[0] = false;
	      --TmpNbr;
	    }
	  for (int j = 1; j < CurrentNbrState; ++j)
	    {
	      while ((j < CurrentNbrState) && (TmpStateArray[j - 1] == TmpStateArray[j]))
		{
		  TmpStateArrayFlag[j] = false;
		  --TmpNbr;
		  ++j;
		}
	      if (j < CurrentNbrState)
		{
		  if (this->TestXMomentumConstraint(TmpStateArray[j], i) == true)
		    {
		      TmpStateArrayFlag[j] = true;
		    }
		  else
		    {
		      TmpStateArrayFlag[j] = false;
		      --TmpNbr;	      
		    }
		}
	    }
	  TmpHilbertSpaceDimension += TmpNbr;
	}
    }
  delete[] this->StateHighestBit;
  delete[] this->StateDescription;
  if (TmpHilbertSpaceDimension > 0)
    {
      this->StateDescription = new ULONGLONG [TmpHilbertSpaceDimension];
      this->StateHighestBit = new int [TmpHilbertSpaceDimension];
      this->NbrStateInOrbit = new int [TmpHilbertSpaceDimension];
      this->ReorderingSign = new unsigned long [TmpHilbertSpaceDimension];
      int Pos = 0;
      ULONGLONG TmpState;
      ULONGLONG TmpState2;
      int TmpSignature;
      unsigned long TmpReorderingSignUp;
      unsigned long TmpReorderingSignDown;
      unsigned long TmpReorderingSign;
      int TmpNbrParticleUp;
      int TmpNbrParticleDown;
      int TmpNbrParticle;
      for (int i = 0; i <maxHighestBit; ++i) 
	{
	  int CurrentNbrState = TmpNbrStateDescription[i];
	  if (CurrentNbrState > 0)
	    {
	      ULONGLONG* TmpStateArray = TmpStateDescription[i];
	      bool* TmpStateArrayFlag = TmpStateDescriptionFlag[i];
	      for (int j = 0; j < CurrentNbrState; ++j)
		{
		  // old code:
		  if (TmpStateArrayFlag[j] == true)
		    {
		      this->StateDescription[Pos] = TmpStateArray[j];
		      this->StateHighestBit[Pos] = i;
		      this->NbrStateInOrbit[Pos] = this->FindNumberXTranslation(this->StateDescription[Pos]);		  
		      TmpState = this->StateDescription[Pos];
		      TmpSignature = 0;
		      TmpReorderingSign = 0x0ul;
		      if ((this->NbrFermions & 1) == 0)
			{
			  for (int k = 1; k <= this->NbrStateInOrbit[Pos]; ++k)
			    {
			      TmpState2 = TmpState & this->MomentumMask;
			      TmpState =  (TmpState >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
			      TmpNbrParticle = 0;
			      for (int l = 0; l < this->StateShift; ++l)
				{
				  if (TmpState2 & ((ULONGLONG) 0x1ul))
				    ++TmpNbrParticle;
				  TmpState2 >>= 1;
				}
			      if (TmpNbrParticle & 1)
				{			 
				  TmpReorderingSign |= (((TmpReorderingSign << 1) & (0x1ul << k))) ^ (0x1ul << k);
				  ++TmpSignature;
				}
			      else
				{
				  TmpReorderingSign |= ((TmpReorderingSign << 1) & (0x1ul << k));
				}			  
			    }
			}
		      // do not assign this result, just keep for comparison
		      // this->ReorderingSign[Pos] = TmpReorderingSign;
		      // ++Pos;
		    }
		  // end old code
		  // new code:
		  if (TmpStateArrayFlag[j] == true)
		    {
		      this->StateDescription[Pos] = TmpStateArray[j];
		      this->StateHighestBit[Pos] = i;
		      this->NbrStateInOrbit[Pos] = this->FindNumberXTranslation(this->StateDescription[Pos]);		  
		      TmpState = this->StateDescription[Pos];
		      TmpSignature = 0;
		      TmpReorderingSignUp = 0x0ul;
		      TmpReorderingSignDown = 0x0ul;
		      if (((this->NbrFermionsUp & 1) == 0)||((this->NbrFermionsDown & 1) == 0))
			{
			  for (int k = 1; k <= this->NbrStateInOrbit[Pos]; ++k)
			    {
			      TmpState2 = TmpState & this->MomentumMask;
			      TmpState =  (TmpState >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
			      TmpNbrParticleUp = 0;
			      TmpNbrParticleDown = 0;
			      for (int l = 0; l < this->StateShift; l+=2)
				{
				  if ((TmpState2 & ((ULONGLONG) 0x1ul)) != ((ULONGLONG) 0x0ul))
				    ++TmpNbrParticleDown;
				  if ((TmpState2 & ((ULONGLONG) 0x2ul)) != ((ULONGLONG) 0x0ul))
				    ++TmpNbrParticleUp;
				  TmpState2 >>= 2;
				}
			      if ((TmpNbrParticleUp & 1)&&((this->NbrFermionsUp & 1) == 0))
				{			 
				  TmpReorderingSignUp |= (((TmpReorderingSignUp << 1) & (0x1ul << k))) ^ (0x1ul << k);
				  ++TmpSignature;
				}
			      else
				{
				  TmpReorderingSignUp |= ((TmpReorderingSignUp << 1) & (0x1ul << k));
				}
			      if ((TmpNbrParticleDown & 1)&&((this->NbrFermionsDown & 1) == 0))
				{			 
				  TmpReorderingSignDown |= (((TmpReorderingSignDown << 1) & (0x1ul << k))) ^ (0x1ul << k);
				  ++TmpSignature;
				}
			      else
				{
				  TmpReorderingSignDown |= ((TmpReorderingSignDown << 1) & (0x1ul << k));
				}			  
			      
			    }
			}
		      this->ReorderingSign[Pos] = TmpReorderingSign;// TmpReorderingSignUp^TmpReorderingSignDown;
		      ++Pos;
		    }
		  
		}
	      delete[] TmpStateArray;
	      delete[] TmpStateArrayFlag;
	    }
	}
      delete[] TmpStateDescription;
      delete[] TmpStateDescriptionFlag;
      delete[] TmpNbrStateDescription;
    }

  return TmpHilbertSpaceDimension;
}

// evaluate Hilbert space dimension without the translation symmmetry along x, Sz or Ky quantu numbers
//
// nbrFermions = number of fermions
// currentMomentum = current one-body momentum 
// pos = position in StateDescription array where to store states
// return value = Hilbert space dimension

long FermionOnTorusWithSpinAndMagneticTranslationsLong::RawGenerateStates(int nbrFermions, int currentMomentum, long pos)
{
  if (nbrFermions == 0)
    {
      this->StateDescription[pos] == ((ULONGLONG) 0x0ul);	  
      return (pos + 1l);
    }
  if ((currentMomentum < 0) || (nbrFermions < 0))
    return pos;
  if (nbrFermions == 1)
    {
      for (int j = currentMomentum; j >= 0; --j)
	{
	  this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << (j << 1);
	  ++pos;
	  this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << (j << 1);
	  ++pos;
	}
      return pos;
    }
  long TmpPos = this->RawGenerateStates(nbrFermions - 2, currentMomentum - 1, pos);
ULONGLONG Mask = ((ULONGLONG) 0x3ul) << (currentMomentum << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStates(nbrFermions - 1, currentMomentum - 1, pos);
  Mask = ((ULONGLONG) 0x2ul) << ((currentMomentum) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
   TmpPos = this->RawGenerateStates(nbrFermions - 1, currentMomentum - 1, pos);
   Mask = ((ULONGLONG) 0x1ul) << (currentMomentum << 1);
   for (; pos < TmpPos; ++pos)
     this->StateDescription[pos] |= Mask;
  return this->RawGenerateStates(nbrFermions, currentMomentum - 1, pos);
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalMomentum = momentum total value
// totalSpinUp = number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = Hilbert space dimension

long FermionOnTorusWithSpinAndMagneticTranslationsLong::RawGenerateStates(int nbrFermions, int lzMax, int totalMomentum, int totalSpinUp, long pos)
{  
  if ((nbrFermions < 0) || (totalSpinUp < 0) || (totalSpinUp > nbrFermions))
    return pos;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpinUp) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpinUp)) )
    return pos;
  if ((nbrFermions == 0) && (lzMax == 0) && (totalSpinUp == 0) && ((totalMomentum % MaxMomentum)== YMomentum))
    {
      this->StateDescription[pos] = ((ULONGLONG) 0x0ul);
      return (pos + 1l);
    }
  
  if (nbrFermions == 1)
    {
      for (int k = 0; k <= lzMax; ++k)
	if (((k+totalMomentum) % MaxMomentum) == YMomentum)
	  {	    
	    this->StateDescription[pos++] = ((ULONGLONG) 0x1ul) << ((k << 1) + totalSpinUp);
	  }
      return (pos);
    }
  
  if ((lzMax == 0)  && ( (totalMomentum % MaxMomentum) != YMomentum))
    return pos;

  // put last two particles:
  if ((lzMax==0) && (nbrFermions == 2) && (totalSpinUp == 1))
    {
      this->StateDescription[pos] = ((ULONGLONG) 0x3ul);
      return (pos + 1l);
    }

  // enter recursion, here:
  // put two particles
  long TmpPos = this->RawGenerateStates(nbrFermions - 2, lzMax - 1, totalMomentum + (2 * lzMax), totalSpinUp - 1,  pos);
  ULONGLONG Mask = ((ULONGLONG) 0x3ul) << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  // put one particle with spin up
  TmpPos = this->RawGenerateStates(nbrFermions - 1, lzMax - 1, totalMomentum + lzMax, totalSpinUp - 1,  pos);
  Mask = ((ULONGLONG) 0x2ul) << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  // put one particle with spin down
  TmpPos = this->RawGenerateStates(nbrFermions - 1, lzMax - 1, totalMomentum + lzMax, totalSpinUp,  pos);
  Mask = ((ULONGLONG) 0x1ul) << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  // do not put any particles
  return this->RawGenerateStates(nbrFermions, lzMax - 1, totalMomentum, totalSpinUp, pos);  
  
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

long FermionOnTorusWithSpinAndMagneticTranslationsLong::RawGenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, long pos, 
								     int currentMomentum)
{
  if ((nbrFermions < 0) || (nbrFermions > (currentMaxMomentum + 1)) || ((nbrFermions > 0) && (currentMaxMomentum < 0)))
    return pos;
  if (nbrFermions == 0)
    {
      if ((currentMomentum % this->MaxMomentum) == this->YMomentum)
	{
	  this->StateDescription[pos] = ((ULONGLONG) 0x0ul);
	  this->StateHighestBit[pos] = maxMomentum >> 1;
	  ++pos;
	}
      return pos;
    }
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  long TmpPos = this->RawGenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentMaxMomentum, pos, currentMomentum + (currentMaxMomentum >> 1));
  ULONGLONG Mask = ((ULONGLONG) 0x1ul) << currentMaxMomentum;
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  if (maxMomentum == currentMaxMomentum)
    return this->RawGenerateStates(nbrFermions, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentMomentum);
  else
    return this->RawGenerateStates(nbrFermions, maxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentMomentum);
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrSpinUp = number of fermions with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnTorusWithSpinAndMagneticTranslationsLong::RawGenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos)
{
  if ((nbrFermions == 0) && (nbrSpinUp == 0))
    {
      this->StateDescription[pos] = ((ULONGLONG) 0x0ul);	  
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0) || (nbrSpinUp > nbrFermions) || (nbrSpinUp < 0))
    return pos;
  if (nbrFermions == 1)
    {
      if (nbrSpinUp == 1)
	{
	  for (int j = currentSite; j >= 0; --j)
	    {
	      this->StateDescription[pos] = ((ULONGLONG) 0x2ul) << (j << 1);
	      ++pos;
	    }
	}
      else
	{
	  for (int j = currentSite; j >= 0; --j)
	    {
	      this->StateDescription[pos] = ((ULONGLONG) 0x1ul) << (j << 1);
	      ++pos;
	    }
	}
      return pos;
    }
  long TmpPos = this->RawGenerateStates(nbrFermions - 2, currentSite - 1, nbrSpinUp - 1, pos);
  ULONGLONG Mask = ((ULONGLONG) 0x3ul) << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, nbrSpinUp - 1, pos);
  Mask = ((ULONGLONG) 0x2ul) << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, nbrSpinUp, pos);
  Mask = ((ULONGLONG) 0x1ul) << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->RawGenerateStates(nbrFermions, currentSite - 1, nbrSpinUp, pos);   
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnTorusWithSpinAndMagneticTranslationsLong::GenerateLookUpTable(int memory)
{
  // evaluate look-up table size
  memory /= (4 * this->NbrFermionStates);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrFermionStates)
    this->MaximumLookUpShift = this->NbrFermionStates;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrFermionStates];
  this->LookUpTableShift = new int [this->NbrFermionStates];
  for (int i = 0; i < this->NbrFermionStates; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentMaxMomentum = this->StateHighestBit[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
  if (CurrentMaxMomentum < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentMaxMomentum] = 0;
  else
    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
  long CurrentLookUpTableValue = 0;
  long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue < TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      ++CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentMaxMomentum != this->StateHighestBit[i])
	{
	  ++CurrentLookUpTableValue;
	  while (CurrentLookUpTableValue <= this->LookUpTableMemorySize)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      ++CurrentLookUpTableValue;
	    }
 	  CurrentMaxMomentum = this->StateHighestBit[i];
	  TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
	  if (CurrentMaxMomentum < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentMaxMomentum] = 0;
	  else
	    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  CurrentLookUpTableValue = 0;
	  while (CurrentLookUpTableValue < TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      ++CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      ++CurrentLookUpTableValue;
	      while (CurrentLookUpTableValue < TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  ++CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  ++CurrentLookUpTableValue;
  while (CurrentLookUpTableValue <= this->LookUpTableMemorySize)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      ++CurrentLookUpTableValue;
    }
  
  this->RescalingFactors = new double* [this->NbrFermionStates];
  for (int i = 1; i <= this->MaxMomentum; ++i)
    {
      this->RescalingFactors[i] = new double [this->NbrFermionStates];
      for (int j = 1; j <= this->MaxMomentum; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}

// generate look-up table associated to sign calculations
// 

void FermionOnTorusWithSpinAndMagneticTranslationsLong::GenerateSignLookUpTable()
{
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  this->NbrParticleLookUpTable = new int [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      Count = 0;
      TmpNbr = j;
      while (TmpNbr != 0)
	{
	  if (TmpNbr & 0x1)
	    ++Count;
	  TmpNbr >>= 1;
	}
      if (Count & 1)
	{
	  this->SignLookUpTable[j] = -1.0;
	  this->NbrParticleLookUpTable[j] = 1;
	}
      else
	{
	  this->SignLookUpTable[j] = 1.0;
	  this->NbrParticleLookUpTable[j] = 0;
	}
    }
#ifdef __128_BIT_LONGLONG__
  this->SignLookUpTableMask = new ULONGLONG [256];
  for (int i = 0; i < 112; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0xffff;
  for (int i = 112; i < 128; ++i)
    this->SignLookUpTableMask[i] = ((ULONGLONG) 0xffff) >> (i - 112);
  for (int i = 128; i < 256; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0;  
#else
  this->SignLookUpTableMask = new ULONGLONG [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((ULONGLONG) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (ULONGLONG) 0;
#endif
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalMomentum = momentum total value
// totalSpinUp = number of particles with spin up
// return value = Hilbert space dimension

long FermionOnTorusWithSpinAndMagneticTranslationsLong::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalMomentum, int totalSpinUp)
{  
  if ((nbrFermions < 0) || (totalSpinUp < 0) || (totalSpinUp > nbrFermions))
    return 0l;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpinUp) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpinUp)) )
    return 0l;

  if (nbrFermions == 1)
    {
      long Tmp = 0;
      for (int k = 0; k <= lzMax; ++k)
	if (((k+totalMomentum) % MaxMomentum) == YMomentum)
	  ++Tmp;
      return Tmp;
    }
  

  if ((lzMax == 0)  && ( (totalMomentum % MaxMomentum) != YMomentum))
    return 0l;

  unsigned long Tmp = 0l;  
  if (nbrFermions > 2)
    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalMomentum + (2 * lzMax), totalSpinUp - 1);
  else
    if ((totalSpinUp == 1) && (((totalMomentum + 2*lzMax) % MaxMomentum) == YMomentum) )
      ++Tmp;
  return  (Tmp + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalMomentum + lzMax, totalSpinUp - 1)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalMomentum + lzMax, totalSpinUp)
	   + this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalMomentum, totalSpinUp));
}

// evaluate Hilbert space dimension for a given total momentum
//
// nbrFermions = number of fermions
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnTorusWithSpinAndMagneticTranslationsLong::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKy, int currentTotalKy)
{
  if (nbrFermions == 0)
    {
      if ((currentTotalKy % this->MaxMomentum) == this->YMomentum)
	return 1l;
      else	
	return 0l;
    }
  if (currentKy < 0)
    return 0l;
  long  Count = 0l;
  if (nbrFermions > 1)
    Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKy - 1, currentTotalKy + (2 * currentKy));
  Count += 2l * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKy - 1, currentTotalKy + currentKy);
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKy - 1, currentTotalKy);
  return Count;
}



// convert a given state from a generic basis to the current Sz subspace basis
//
// state = reference on the vector to convert
// basis = reference on the basis associated to state
// return value = converted vector

ComplexVector FermionOnTorusWithSpinAndMagneticTranslationsLong::ConvertFromNbodyBasis(ComplexVector& state, FermionOnTorusWithSpinAndMagneticTranslationsLong& basis)
{
  ComplexVector TmpVector (basis.GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
  {
    TmpVector[basis.FindStateIndex(this->StateDescription[i], this->StateHighestBit[i])] = state[i];
  }
  return TmpVector;
}

// apply a Gutzwiller projection (in the orbital space) to a given state
//
// state = reference on the state to project
// space = pointer to the Hilbert space where state is defined
// return value = Gutzwiller projected state

ComplexVector FermionOnTorusWithSpinAndMagneticTranslationsLong::GutzwillerProjection(ComplexVector& state, ParticleOnSphere* space)
{
  FermionOnTorusWithSpinAndMagneticTranslationsLong* TmpSpace = (FermionOnTorusWithSpinAndMagneticTranslationsLong*) space;
  ComplexVector ProjectedState (this->LargeHilbertSpaceDimension, true);
#ifdef __128_BIT_LONGLONG__
  ULONGLONG TmpSpinUpMask = (((ULONGLONG) 0xaaaaaaaaaaaaaaaaul) << 64) | ((ULONGLONG) 0xaaaaaaaaaaaaaaaaul);
  ULONGLONG TmpSpinDownMask = (((ULONGLONG) 0x5555555555555555ul) << 64) | ((ULONGLONG) 0x5555555555555555ul);  
#else
  ULONGLONG TmpSpinUpMask   = (((ULONGLONG) 0xaaaaaaaaul) << 64) | ((ULONGLONG) 0xaaaaaaaaul);
  ULONGLONG TmpSpinDownMask = (((ULONGLONG) 0x55555555ul) << 64) | ((ULONGLONG) 0x55555555ul);  
#endif	    
  for (long i = 0l; i < TmpSpace->LargeHilbertSpaceDimension; ++i)
    {
      ULONGLONG TmpState = TmpSpace->StateDescription[i];
      if ((((TmpState & TmpSpinUpMask) >> 1) & (TmpState & TmpSpinDownMask)) == ((ULONGLONG)0x0ul))
	{
	  int TmpIndex = this->FindStateIndex(TmpState, TmpSpace->StateHighestBit[i]);
	  if (TmpIndex < this->HilbertSpaceDimension)
	    {
	      ProjectedState[TmpIndex] = state[i];
	    }
	}
    }
  return ProjectedState;
}

