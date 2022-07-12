////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                   class author: Gunnar Möller                              //
//                                                                            //
//   class of fermion on a torus taking into account magnetic translations    //
//                                                                            //
//                        last modification : 26/11/2007                      //
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
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnTorusWithSpinNew.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
//#include "QuantumNumber/VectorQuantumNumber.h"
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
// nbrFermions= number of fermions
// totalSpin = twice the total spin value
// maxMomentum = momentum maximum value for a fermion
// yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)

FermionOnTorusWithSpinNew::FermionOnTorusWithSpinNew (int nbrFermions, int totalSpin, int maxMomentum, 
						      int yMomentum)
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
  this->YMomentum = yMomentum % this->MaxMomentum;

  this->StateShift = 2*(this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2*this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((unsigned long) 1);
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= ((unsigned long) 1);
    }

  this->MaximumSignLookUp = 16;
  this->GenerateSignLookUpTable();
  // testing the count of states:
  if (false)
    {
      int tmpYMomentum=YMomentum;
      long sum=0;
      for ( int i=0; i<this->MaxMomentum; ++i)
	{
	  this->YMomentum=i;
	  this->HilbertSpaceDimension = ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum - 1, 0, this->NbrFermionsUp);
	  sum += this->HilbertSpaceDimension;
	  cout << "YMomentum="<<YMomentum << "\tD="<<this->HilbertSpaceDimension<<endl;
	}
      cout << "sum D= "<<sum<<endl;
      cout << "A priori = " << this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->TotalSpin, this->MaxMomentum)<<endl;
      YMomentum=tmpYMomentum;

      if (abs(TotalSpin)==NbrFermions)
	{
	  cout << "comparing:" << endl;
	  FermionOnTorus *test=new FermionOnTorus	    (NbrFermions, MaxMomentum, yMomentum);
	  test->GetHilbertSpaceDimension();
	  cout << "done comparing:" << endl;
	}
    }
  // end test section
  cout << "YMomentum="<<YMomentum <<endl;
  this->HilbertSpaceDimension = ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum - 1, 0, this->NbrFermionsUp);
  cout << "Full dimension: "<<HilbertSpaceDimension<<endl;
  this->HilbertSpaceDimension = this->GenerateStates();
  this->Flag.Initialize();

  if (false)
    {
      for (int i=0; i<this->HilbertSpaceDimension; ++i)
	PrintState(cout,i)<<endl;
    }

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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnTorusWithSpinNew::FermionOnTorusWithSpinNew(const FermionOnTorusWithSpinNew& fermions)
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

  this->Flag = fermions.Flag;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnTorusWithSpinNew::~FermionOnTorusWithSpinNew ()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;
      delete[] this->SignLookUpTableMask;      
  
      if (this->HilbertSpaceDimension != 0)
	{
	  delete[] this->LookUpTableShift;
	  for (int i = 0; i < this->NbrFermionStates; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;	  
	  
	}
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnTorusWithSpinNew& FermionOnTorusWithSpinNew::operator = (const FermionOnTorusWithSpinNew& fermions)
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

  this->Flag = fermions.Flag;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnTorusWithSpinNew::Clone()
{
  return new FermionOnTorusWithSpinNew(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnTorusWithSpinNew::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->YMomentum, this->MomentumModulo);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnTorusWithSpinNew::GetQuantumNumber (int index)
{
  return new PeriodicMomentumQuantumNumber (this->YMomentum, this->MomentumModulo);
}

// get momemtum value in the y direction of a given state
//
// index = state index
// return value = state momentum in the y direction

int FermionOnTorusWithSpinNew::GetYMomentumValue(int index)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  int Momentum = 0;
  for (int i = 0; i <= StateHighestBit; ++i)
    {
      Momentum += ((State >> (2*i) ) & ((unsigned long) 1)+(State >> (2*i+1) ) & ((unsigned long) 1)) * i;
    }
  return (Momentum % this->MomentumModulo);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnTorusWithSpinNew::ExtractSubspace (AbstractQuantumNumber& q, 
									       SubspaceSpaceConverter& converter)
{
  return 0;
}


// apply sum_m au^+_m au_m operator to a given state
//
// index = index of the state on which the operator has to be applied
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpinNew::SumAduAu (int index, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned int State = this->StateDescription[index];
  coefficient = 0.0;
  for (int m = 1; m <= StateHighestBit; m+=2)
    if ((State & (1 << m)) != 0)
      {
	coefficient += 1.0;
      }
  return index;
}


// apply sum_m ad^+_m ad_m operator to a given state
//
// index = index of the state on which the operator has to be applied
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpinNew::SumAddAd (int index, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned int State = this->StateDescription[index];
  coefficient = 0.0;
  for (int m = 0; m < StateHighestBit; m+=2)
    if ((State & (1 << m)) != 0)
      {
	coefficient += 1.0;
      }
  return index;
}

// apply au^+_m au_m operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for density operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpinNew::AduAu (int index, int m, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned int State = this->StateDescription[index];
  if ((m > StateHighestBit) || ((State & (2 << (m << 1))) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  return index;
}

// apply ad^+_m ad_m operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for density operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpinNew::AddAd (int index, int m, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned int State = this->StateDescription[index];
  if ((m > StateHighestBit) || ((State & (1 << (m << 1))) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  return index;
}


// apply a^+_(d,m1) a^+_(d,m2) a_(d,n1) a_(d,n2) operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n2<<=1;  
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || 
      ((State & ((0x1ul) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  m1<<=1;
  m2<<=1;
  if ((TmpState & (((unsigned long) 1) << m2))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  if ((TmpState & (((unsigned long) 0x1ul) << m1))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
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
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || 
      ((State & ((0x1ul) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;  
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  if ((TmpState & (((unsigned long) 0x1) << m2))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  if ((TmpState & (((unsigned long) 0x1) << m1))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
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
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AduAduAuAuV (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || 
      ((State & ((0x1ul) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n2);
  cout << "Sign after n2="<<n2<<" in AduAduAuAu:" << coefficient << endl;
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n1);
  cout << "Sign after n1="<<n1<<" in AduAduAuAu:" << coefficient << endl;
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  cout << "TmpState in AduAduAuAu=" << TmpState << endl;
  if ((TmpState & (((unsigned long) 1) << m2))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  cout << "Sign after m2="<<m2<<" in AduAduAuAu:" << coefficient << endl;
  if ((TmpState & (((unsigned long) 1) << m1))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  cout << "Sign after m1="<<m1<<" in AduAduAuAu:" << coefficient << endl;
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
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
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AduAddAuAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n1+=1;
  n2<<=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || ((State & ((0x1ul) << n2)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  m1<<=1;
  m1+=1;
  m2<<=1;
  if ((TmpState & (((unsigned long) 1) << m2))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  if ((TmpState & (((unsigned long) 0x1) << m1))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
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
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AddAddAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || 
      ((State & ((0x1ul) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;  
  m1<<=1;
  m2<<=1;
  if ((TmpState & (((unsigned long) 0x1) << m2))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  if ((TmpState & (((unsigned long) 0x1) << m1))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
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
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AduAduAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n2<<=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || 
      ((State & ((0x1ul) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;  
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  if ((TmpState & (((unsigned long) 0x1) << m2))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  if ((TmpState & (((unsigned long) 0x1) << m1))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
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
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AduAduAuAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n1+=1;
  n2<<=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || 
      ((State & ((0x1ul) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;  
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  if ((TmpState & (((unsigned long) 0x1) << m2))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  if ((TmpState & (((unsigned long) 0x1) << m1))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
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
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AddAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || 
      ((State & ((0x1ul) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;  
  m1<<=1;
  m2<<=1;
  m2+=1;
  if ((TmpState & (((unsigned long) 0x1) << m2))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  if ((TmpState & (((unsigned long) 0x1) << m1))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
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
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AddAddAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n2<<=1;
  n2+=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || 
      ((State & ((0x1ul) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;  
  m1<<=1;
  m2<<=1;
  if ((TmpState & (((unsigned long) 0x1) << m2))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  if ((TmpState & (((unsigned long) 0x1) << m1))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
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
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AduAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n2<<=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || 
      ((State & ((0x1ul) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;  
  m1<<=1;
  m1+=1;
  m2<<=1;
  if ((TmpState & (((unsigned long) 0x1) << m2))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  if ((TmpState & (((unsigned long) 0x1) << m1))!= 0)
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
  return TmpIndex;
}


// apply a^+_m_d a_m_d operator to a given state (only spin down)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
double FermionOnTorusWithSpinNew::AddAd (int index, int m)
{
  if (this->StateHighestBit[index] < 2*m)
    return 0.0;
  if ((this->StateDescription[index] & (((unsigned long) 0x1) << (m<<1))) == 0)
    return 0.0;
  else
    return 1.0;
}
  

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
double FermionOnTorusWithSpinNew::AduAu (int index, int m)
{
  if (this->StateHighestBit[index] < 2*m+1)
    return 0.0;
  if ((this->StateDescription[index] & (((unsigned long) 0x1) << (2*m+1))) == 0)
    return 0.0;
  else
    return 1.0;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinNew::AuAu (int index, int n1, int n2)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((0x1ul) << n1)) == 0) || 
      ((ProdATemporaryState & ((0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;  
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  ProdATemporaryState &= ~(((unsigned long) 0x1) << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  ProdATemporaryState &= ~(((unsigned long) 0x1) << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  return coefficient;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinNew::AuAuV (int index, int n1, int n2)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((0x1ul) << n1)) == 0) || 
      ((ProdATemporaryState & ((0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;  
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  ProdATemporaryState &= ~(((unsigned long) 0x1) << n2);
  cout << "Sign after n2="<<n2<<" in AuAu:" << coefficient << endl;
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  ProdATemporaryState &= ~(((unsigned long) 0x1) << n1);
  cout << "Sign after n1="<<n1<<" in AuAu:" << coefficient << endl;
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  return coefficient;
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinNew::AdAd (int index, int n1, int n2)
{  
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n2<<=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((0x1ul) << n1)) == 0) || 
      ((ProdATemporaryState & ((0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  ProdATemporaryState &= ~(((unsigned long) 0x1) << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  ProdATemporaryState &= ~(((unsigned long) 0x1) << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  return coefficient;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinNew::AuAd (int index, int n1, int n2)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n1+=1;
  n2<<=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((0x1ul) << n1)) == 0) || 
      ((ProdATemporaryState & ((0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  ProdATemporaryState &= ~(((unsigned long) 0x1) << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  ProdATemporaryState &= ~(((unsigned long) 0x1) << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  return coefficient;
}

// apply a_n1_d a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinNew::AdAu (int index, int n1, int n2)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n2<<=1;
  n2+=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((0x1ul) << n1)) == 0) || 
      ((ProdATemporaryState & ((0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  ProdATemporaryState &= ~(((unsigned long) 0x1) << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  ProdATemporaryState &= ~(((unsigned long) 0x1) << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  return coefficient;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AduAdu (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  int TmpAHighestBit = ProdAHighestBit;
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  
  if (m1 > TmpAHighestBit)
    TmpAHighestBit = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, TmpAHighestBit);
  return TmpIndex;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AduAduV (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  int TmpAHighestBit = ProdAHighestBit;
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  cout << "TmpState in AduAduV=" << TmpState << endl;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  cout << "Sign after m2="<<m2<<" in AduAdu:" << coefficient << endl;
  if (m1 > TmpAHighestBit)
    TmpAHighestBit = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  cout << "Sign after m1="<<m1<<" in AduAdu:" << coefficient << endl;
  int TmpIndex = this->FindStateIndex(TmpState, TmpAHighestBit);
  return TmpIndex;
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AddAdd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  int TmpAHighestBit = ProdAHighestBit;
  m1<<=1;
  m2<<=1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  
  if (m1 > TmpAHighestBit)
    TmpAHighestBit = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, TmpAHighestBit);
  return TmpIndex;
}

  

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnTorusWithSpinNew::AduAdd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  int TmpAHighestBit = ProdAHighestBit;
  m1<<=1;
  m1+=1;
  m2<<=1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0))
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
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  
  if (m1 > TmpAHighestBit)
    TmpAHighestBit = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  int TmpIndex = this->FindStateIndex(TmpState, TmpAHighestBit);
  return TmpIndex;
}


// find state index
//
// stateDescription = unsigned longeger describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnTorusWithSpinNew::FindStateIndex(unsigned long stateDescription, int maxMomentum)
{
  long PosMax = stateDescription >> this->LookUpTableShift[maxMomentum];
  long PosMin = this->LookUpTable[maxMomentum][PosMax];
  PosMax = this->LookUpTable[maxMomentum][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
//  cout << PosMin << " " << PosMid << " " << PosMax << " " << hex << CurrentState << " " << stateDescription << dec << endl;
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

ostream& FermionOnTorusWithSpinNew::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->MaxMomentum; ++i)
    Str << ((TmpState >> (2 * i)) & 0x1ul) << ((TmpState >> ((2 * i) + 1)) & 0x1ul) << " ";
//   Str << "  (" << hex << this->ReorderingSign[state] << dec << ")";
   Str << " " << this->FindStateIndex(this->StateDescription[state], this->StateHighestBit[state]);
//   if (this->FindStateIndex(this->StateDescription[state], this->StateHighestBit[state]) != state)
//     {
//       Str << "  error";
//     }
  return Str;
}

// generate all states corresponding to the constraints
// 
// return value = hilbert space dimension

int FermionOnTorusWithSpinNew::GenerateStates()
{
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  this->HilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, this->MaxMomentum - 1, 0, this->NbrFermionsUp, 0);
  return this->HilbertSpaceDimension;
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalMomentum = momentum total value
// totalSpinUp = number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = Hilbert space dimension

long FermionOnTorusWithSpinNew::RawGenerateStates(int nbrFermions, int lzMax, int totalMomentum, int totalSpinUp, long pos)
{  
  if ((nbrFermions < 0) || (totalSpinUp < 0) || (totalSpinUp > nbrFermions))
    return pos;
  if ((lzMax < 0) || ((2 * (lzMax + 1)) < totalSpinUp) || ((2 * (lzMax + 1)) < (nbrFermions - totalSpinUp)) )
    return pos;
  if ((nbrFermions == 0) && (lzMax == 0) && (totalSpinUp == 0) && ((totalMomentum % MaxMomentum)== YMomentum))
    {
      this->StateDescription[pos] = 0x0ul;
      this->StateHighestBit[pos] = 0;
      return (pos + 1l);
    }
  
  if (nbrFermions == 1)
    {
      for (int k = 0; k <= lzMax; ++k)
	if (((k+totalMomentum) % MaxMomentum) == YMomentum)
	  {
	    this->StateHighestBit[pos] = ((k << 1) + totalSpinUp);
	    this->StateDescription[pos++] = 0x1ul << ((k << 1) + totalSpinUp);
	  }
      return (pos);
    }
  
  if ((lzMax == 0)  && ( (totalMomentum % MaxMomentum) != YMomentum))
    return pos;

  // put last two particles:
  if ((lzMax==0) && (nbrFermions == 2) && (totalSpinUp == 1))
    {
      this->StateHighestBit[pos] = 1;
      this->StateDescription[pos] = 0x3ul;
      return (pos + 1l);
    }

  // enter recursion, here:
  // put two particles
  long TmpPos = this->RawGenerateStates(nbrFermions - 2, lzMax - 1, totalMomentum + (2 * lzMax), totalSpinUp - 1,  pos);
  unsigned long Mask = 0x3ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    {
      this->StateHighestBit[pos] = (lzMax << 1)+1;
      this->StateDescription[pos] |= Mask;
    }
  // put one particle with spin up
  TmpPos = this->RawGenerateStates(nbrFermions - 1, lzMax - 1, totalMomentum + lzMax, totalSpinUp - 1,  pos);
  Mask = 0x2ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    {
      this->StateHighestBit[pos] = (lzMax << 1)+1;
      this->StateDescription[pos] |= Mask;
    }
  // put one particle with spin down
  TmpPos = this->RawGenerateStates(nbrFermions - 1, lzMax - 1, totalMomentum + lzMax, totalSpinUp,  pos);
  Mask = 0x1ul << (lzMax << 1);
  for (; pos < TmpPos; ++pos)
    {
      this->StateHighestBit[pos] = (lzMax << 1);
      this->StateDescription[pos] |= Mask;
    }
  // do not put any particles
  return this->RawGenerateStates(nbrFermions, lzMax - 1, totalMomentum, totalSpinUp, pos);  
  
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnTorusWithSpinNew::GenerateLookUpTable(int memory)
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
}

// generate look-up table associated to sign calculations
// 

void FermionOnTorusWithSpinNew::GenerateSignLookUpTable()
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
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#endif
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// totalSpin = twice z-projection of total spin
// maxMomentum = momentum maximum value for a fermion
// return value = Hilbert space dimension
int FermionOnTorusWithSpinNew::EvaluateHilbertSpaceDimension(int nbrFermions, int totalSpin, int maxMomentum)
{
  FactorialCoefficient Dimension;
  int nbrFermionsUp = (nbrFermions+totalSpin)/2;
  int nbrFermionsDown = (nbrFermions-totalSpin)/2;  
  Dimension.PartialFactorialMultiply(maxMomentum - nbrFermionsUp + 1, maxMomentum); 
  Dimension.FactorialDivide(nbrFermionsUp);
  Dimension.PartialFactorialMultiply(maxMomentum - nbrFermionsDown + 1, maxMomentum); 
  Dimension.FactorialDivide(nbrFermionsDown);
  return Dimension.GetIntegerValue();
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalMomentum = momentum total value
// totalSpinUp = number of particles with spin up
// return value = Hilbert space dimension

long FermionOnTorusWithSpinNew::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalMomentum, int totalSpinUp)
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
