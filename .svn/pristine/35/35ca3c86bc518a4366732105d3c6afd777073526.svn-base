////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of fermions on a torus with spin                  //
//                                                                            //
//                        last modification : 10/09/2002                      //
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
#include "HilbertSpace/FermionOnTorusWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include <math.h>


// basic constructor
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion

FermionOnTorusWithSpin::FermionOnTorusWithSpin (int nbrFermions, int maxMomentum)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->TotalSpinMomentumFlag = false;
  this->MomentumConstraint = 0;
  this->MomentumConstraintFlag = false;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum);
  cout << this->HilbertSpaceDimension << endl;
  this->Flag.Initialize();
  this->StateDescription = new unsigned int [this->HilbertSpaceDimension];
  this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions,  2 * this->MaxMomentum - 1, 2 * this->MaxMomentum - 1, 0);
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
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

// constructor with a constraint on total spin momentum and total momentum
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// totalSpinMomentum = twice the total spin momentum to be used as constraint
// momentumConstraint = index of the momentum orbit

FermionOnTorusWithSpin::FermionOnTorusWithSpin (int nbrFermions, int maxMomentum, int totalSpinMomentum, int momentumConstaint)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->TotalSpinMomentum = totalSpinMomentum;
  this->TotalSpinMomentumFlag = true;
  this->MomentumConstraint = momentumConstaint;
  this->MomentumConstraintFlag = true;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum, totalSpinMomentum);
  cout << this->HilbertSpaceDimension << endl;
  this->Flag.Initialize();
  this->StateDescription = new unsigned int [this->HilbertSpaceDimension];
  this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermions,  2 * this->MaxMomentum - 1, 2 * this->MaxMomentum - 1, 0, 0, 0);
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
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

// constructor with a constraint on total spin momentum
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// totalSpinMomentum = twice the total spin momentum to be used as constraint

FermionOnTorusWithSpin::FermionOnTorusWithSpin (int nbrFermions, int maxMomentum, int totalSpinMomentum)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->TotalSpinMomentum = totalSpinMomentum;
  this->TotalSpinMomentumFlag = true;
  this->MomentumConstraint = 0;
  this->MomentumConstraintFlag = false;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum, totalSpinMomentum);
  cout << this->HilbertSpaceDimension << endl;
  this->Flag.Initialize();
  this->StateDescription = new unsigned int [this->HilbertSpaceDimension];
  this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrFermions,  2 * this->MaxMomentum - 1, 2 * this->MaxMomentum - 1, 0, 0);
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
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

// constructor from full datas
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// totalSpinMomentum = twice the total spin momentum to be used as constraint
// hilbertSpaceDimension = Hilbert space dimension
// stateDescription = array describing each state
// stateMaxMomentum = array giving maximum Lz value reached for a fermion in a given state
// momentumConstraintFlag = flag for momementum constraint
// momentumConstraint = index of the momentum orbit

FermionOnTorusWithSpin::FermionOnTorusWithSpin (int nbrFermions, int maxMomentum, int totalSpinMomentum, int hilbertSpaceDimension, 
						unsigned int* stateDescription, int* StateMaxMomentum, bool momentumConstraintFlag, int momentumConstraint)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->TotalSpinMomentum = totalSpinMomentum;
  this->TotalSpinMomentumFlag = true;
  this->MomentumConstraint = momentumConstraint;
  this->MomentumConstraintFlag = momentumConstraintFlag;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = stateDescription;
  this->StateMaxMomentum = StateMaxMomentum;
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
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

FermionOnTorusWithSpin::FermionOnTorusWithSpin(const FermionOnTorusWithSpin& fermions)
{
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;
  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpinMomentum = fermions.TotalSpinMomentum;
  this->TotalSpinMomentumFlag = fermions.TotalSpinMomentumFlag;
  this->Flag = fermions.Flag;
}

// destructor
//

FermionOnTorusWithSpin::~FermionOnTorusWithSpin ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnTorusWithSpin& FermionOnTorusWithSpin::operator = (const FermionOnTorusWithSpin& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
    }
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateMaxMomentum = fermions.StateMaxMomentum;
  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnTorusWithSpin::Clone()
{
  return new FermionOnTorusWithSpin(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnTorusWithSpin::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->TotalSpinMomentumFlag == true)
    {
      if (this->MomentumConstraintFlag == true)
	{
	  List<AbstractQuantumNumber*> TmpList2;
	  TmpList2 += new SzQuantumNumber (this->TotalSpinMomentum);
	  TmpList2 += new PeriodicMomentumQuantumNumber (this->MomentumConstraint, this->MaxMomentum);
	  TmpList += new VectorQuantumNumber (TmpList2);
	}
      else
	{
	  for (int i = 0; i < this->MaxMomentum; ++i)
	    {
	      List<AbstractQuantumNumber*> TmpList2;
	      TmpList2 += new SzQuantumNumber (this->TotalSpinMomentum);
	      TmpList2 += new PeriodicMomentumQuantumNumber (i, this->MaxMomentum);
	      TmpList += new VectorQuantumNumber (TmpList2);
	    }
	}
    }
  else
    {
      int Max;
      if (this->NbrFermions > this->MaxMomentum)
	Max = this->MaxMomentum;
      else
	Max = this->NbrFermions;
      for (int i = -Max; i <= Max; i += 2)
	for (int j = 0; j < this->MaxMomentum; ++j)
	  {
	    List<AbstractQuantumNumber*> TmpList2;
	    TmpList2 += new SzQuantumNumber (i);
	    TmpList2 += new PeriodicMomentumQuantumNumber (j, this->MaxMomentum);
	    TmpList += new VectorQuantumNumber (TmpList2);
	  }	
    }
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnTorusWithSpin::GetQuantumNumber (int index)
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->TotalSpinMomentumFlag == true)
    {
      TmpList +=  new SzQuantumNumber (this->TotalSpinMomentum);
    }
  else
    {
      double Coef;
      this->SumAudAu(index, Coef);
      TmpList += new SzQuantumNumber (2 * this->NbrFermions - ((int) Coef));
    }    
  if (this->MomentumConstraintFlag == true)
    {
      TmpList += new PeriodicMomentumQuantumNumber (this->MomentumConstraint, this->MaxMomentum);
    }
  else
    {
      TmpList += new PeriodicMomentumQuantumNumber (this->GetMomentumValue(index), this->MaxMomentum);
    }    
  return new VectorQuantumNumber (TmpList);
}

// get momemtum value of a given state
//
// index = state index
// return value = state momentum

int FermionOnTorusWithSpin::GetMomentumValue(int index)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  int Momentum = 0;
  for (int i = 0; i <= StateMaxMomentum; ++i)
    {
      Momentum += (((State >> (i << 1)) & 1) + ((State >> ((i << 1) + 1)) & 1)) * i;
    }
  return (Momentum % this->MaxMomentum);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnTorusWithSpin::ExtractSubspace (AbstractQuantumNumber& q, 
							       SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() == AbstractQuantumNumber::Sz)
    {
      if (this->TotalSpinMomentumFlag == true)
	if (this->TotalSpinMomentum == ((SzQuantumNumber&) q).GetSz())
	  return this;
	else 
	  return 0;
      double Sz = ((double) (((SzQuantumNumber&) q).GetSz() + this->NbrFermions)) * 0.5;
      double Coef;
      int SubspaceHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum, ((SzQuantumNumber&) q).GetSz());
      unsigned int* SubspaceStateDescription = new unsigned int [SubspaceHilbertSpaceDimension];
      int* SubspaceStateMaxMomentum = new int [SubspaceHilbertSpaceDimension];
      int* ConvArray = new int [SubspaceHilbertSpaceDimension];
      int Count = 0;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  this->SumAudAu(i, Coef);
	  if (Coef == Sz)
	    {
	      SubspaceStateDescription[Count] = this->StateDescription[i];
	      SubspaceStateMaxMomentum[Count] = this->StateMaxMomentum[i];
	      ConvArray[Count] = i;
	      ++Count;
	    }
	}
      converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, SubspaceHilbertSpaceDimension, ConvArray);
      return new FermionOnTorusWithSpin(this->NbrFermions, this->MaxMomentum, ((SzQuantumNumber&) q).GetSz(), SubspaceHilbertSpaceDimension, 
					SubspaceStateDescription, SubspaceStateMaxMomentum);
    }
  if ((q.GetQuantumNumberType() == AbstractQuantumNumber::PeriodicMomentum) && (this->TotalSpinMomentumFlag == true))
    {
      if (this->MomentumConstraintFlag == true)
	if (this->MomentumConstraint == ((PeriodicMomentumQuantumNumber&) q).GetMomentum())
	  return this;
	else 
	  return 0;
      int Momentum = ((PeriodicMomentumQuantumNumber&) q).GetMomentum();
      int* TmpConvArray = new int [this->HilbertSpaceDimension];
      int SubspaceHilbertSpaceDimension = 0;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  if (this->GetMomentumValue(i) == Momentum)
	    {
	      TmpConvArray[SubspaceHilbertSpaceDimension] = i;
	      ++SubspaceHilbertSpaceDimension;
	    }
	}
      unsigned int* SubspaceStateDescription = new unsigned int [SubspaceHilbertSpaceDimension];
      int* SubspaceStateMaxMomentum = new int [SubspaceHilbertSpaceDimension];
      int* ConvArray = new int [SubspaceHilbertSpaceDimension];
      for (int i = 0; i < SubspaceHilbertSpaceDimension; ++i)
	{
	  ConvArray[i] = TmpConvArray[i];
	  SubspaceStateDescription[i] = this->StateDescription[TmpConvArray[i]];
	  SubspaceStateMaxMomentum[i] = this->StateMaxMomentum[TmpConvArray[i]];
	}
      delete[] TmpConvArray;
      converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, SubspaceHilbertSpaceDimension, ConvArray);
      return new FermionOnTorusWithSpin(this->NbrFermions, this->MaxMomentum, this->TotalSpinMomentum, SubspaceHilbertSpaceDimension, 
					SubspaceStateDescription, SubspaceStateMaxMomentum, true, Momentum);
    }
  if ((q.GetQuantumNumberType() & AbstractQuantumNumber::Vector) && 
      (((VectorQuantumNumber&) q)[1]->GetQuantumNumberType() == AbstractQuantumNumber::PeriodicMomentum) && 
      (((VectorQuantumNumber&) q)[0]->GetQuantumNumberType() == AbstractQuantumNumber::Sz))
    {
      PeriodicMomentumQuantumNumber* MomentumNumber = (PeriodicMomentumQuantumNumber*) ((VectorQuantumNumber&) q)[1];
      SzQuantumNumber* SzNumber = (SzQuantumNumber*) ((VectorQuantumNumber&) q)[0];
      if ((this->MomentumConstraintFlag == false) && (this->TotalSpinMomentumFlag == false))
	{
	  int Momentum = MomentumNumber->GetMomentum();
	  double Sz = ((double) (SzNumber->GetSz() + this->NbrFermions)) * 0.5;
	  double Coef;
	  int* TmpConvArray = new int [this->HilbertSpaceDimension];
	  int SubspaceHilbertSpaceDimension = 0;
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      this->SumAudAu(i, Coef);
	      if ((Coef == Sz) && (this->GetMomentumValue(i) == Momentum))
		{
		  TmpConvArray[SubspaceHilbertSpaceDimension] = i;
		  ++SubspaceHilbertSpaceDimension;
		}
	    }
	  unsigned int* SubspaceStateDescription = new unsigned int [SubspaceHilbertSpaceDimension];
	  int* SubspaceStateMaxMomentum = new int [SubspaceHilbertSpaceDimension];
	  int* ConvArray = new int [SubspaceHilbertSpaceDimension];
	  for (int i = 0; i < SubspaceHilbertSpaceDimension; ++i)
	    {
	      ConvArray[i] = TmpConvArray[i];
	      SubspaceStateDescription[i] = this->StateDescription[TmpConvArray[i]];
	      SubspaceStateMaxMomentum[i] = this->StateMaxMomentum[TmpConvArray[i]];
	    }
	  delete[] TmpConvArray;
	  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, SubspaceHilbertSpaceDimension, ConvArray);
	  return new FermionOnTorusWithSpin(this->NbrFermions, this->MaxMomentum, SzNumber->GetSz(), SubspaceHilbertSpaceDimension, 
					    SubspaceStateDescription, SubspaceStateMaxMomentum, true, Momentum);
	}
      if ((this->MomentumConstraintFlag == true) && (this->TotalSpinMomentumFlag == true))
	if ((MomentumNumber->GetMomentum() == this->MomentumConstraint) && (this->TotalSpinMomentum == SzNumber->GetSz()))
	  return this;
	else
	  return 0;
      if ((this->MomentumConstraintFlag == false) && (this->TotalSpinMomentumFlag == true) && (this->TotalSpinMomentum == SzNumber->GetSz()))
	{
	  int Momentum = MomentumNumber->GetMomentum();
	  int* TmpConvArray = new int [this->HilbertSpaceDimension];
	  int SubspaceHilbertSpaceDimension = 0;
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      if (this->GetMomentumValue(i) == Momentum)
		{
		  TmpConvArray[SubspaceHilbertSpaceDimension] = i;
		  ++SubspaceHilbertSpaceDimension;
		}
	    }
	  unsigned int* SubspaceStateDescription = new unsigned int [SubspaceHilbertSpaceDimension];
	  int* SubspaceStateMaxMomentum = new int [SubspaceHilbertSpaceDimension];
	  int* ConvArray = new int [SubspaceHilbertSpaceDimension];
	  for (int i = 0; i < SubspaceHilbertSpaceDimension; ++i)
	    {
	      ConvArray[i] = TmpConvArray[i];
	      SubspaceStateDescription[i] = this->StateDescription[TmpConvArray[i]];
	      SubspaceStateMaxMomentum[i] = this->StateMaxMomentum[TmpConvArray[i]];
	    }
	  delete[] TmpConvArray;
	  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, SubspaceHilbertSpaceDimension, ConvArray);
	  return new FermionOnTorusWithSpin(this->NbrFermions, this->MaxMomentum, this->TotalSpinMomentum, SubspaceHilbertSpaceDimension, 
					    SubspaceStateDescription, SubspaceStateMaxMomentum, true, Momentum);
	}
      return 0;
    }
  return 0;
}

// apply sum_m au^+_m au_m operator to a given state
//
// index = index of the state on which the operator has to be applied
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpin::SumAudAu (int index, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  coefficient = 0.0;
  for (int m = 0; m <= StateMaxMomentum; ++m)
    if ((State & (1 << (m << 1))) != 0)
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

int FermionOnTorusWithSpin::SumAddAd (int index, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  coefficient = 0.0;
  for (int m = 0; m < StateMaxMomentum; ++m)
    if ((State & (2 << (m << 1))) != 0)
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

int FermionOnTorusWithSpin::AudAu (int index, int m, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((m > StateMaxMomentum) || ((State & (1 << (m << 1))) == 0))
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

int FermionOnTorusWithSpin::AddAd (int index, int m, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((m > StateMaxMomentum) || ((State & (2 << (m << 1))) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  return index;
}

// apply au^+_m1 au^+_m2 au_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpin::AudAudAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (1 << (n1 << 1))) == 0) || 
      ((State & (1 << (n2 << 1))) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  int TmpState = State;
  int Mask;
  int SignIndex = (NewMaxMomentum - n2) << 1;
  int Twice = n2 << 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient = this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient = this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n2)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = n1 << 1;
  SignIndex = (NewMaxMomentum - n1) << 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient *= this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n1)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = m2 << 1;
  if ((TmpState & (1 << Twice)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m2;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m2) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  Twice = m1 << 1;
  if ((TmpState & (1 << Twice))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m1;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m1) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  return this->FindStateIndex(TmpState, NewMaxMomentum);
}

// apply ad^+_m1 ad^+_m2 ad_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpin::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (2 << (n1 << 1))) == 0) || 
      ((State & (2 << (n2 << 1))) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  int TmpState = State;
  int Mask;
  int SignIndex = (NewMaxMomentum - n2) << 1;
  int Twice = (n2 << 1) + 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient = this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient = this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n2)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (n1 << 1) + 1;
  SignIndex = (NewMaxMomentum - n1) << 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient *= this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n1)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (m2 << 1) + 1;
  if ((TmpState & (1 << Twice)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m2;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m2) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  Twice = (m1 << 1) + 1;
  if ((TmpState & (1 << Twice))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m1;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m1) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  return this->FindStateIndex(TmpState, NewMaxMomentum);
}

// apply ad^+_m1 ad^+_m2 au_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpin::AddAddAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (1 << (n1 << 1))) == 0) || 
      ((State & (1 << (n2 << 1))) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  int TmpState = State;
  int Mask;
  int SignIndex = (NewMaxMomentum - n2) << 1;
  int Twice = (n2 << 1);
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient = this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient = this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n2)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (n1 << 1);
  SignIndex = (NewMaxMomentum - n1) << 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient *= this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n1)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (m2 << 1) + 1;
  if ((TmpState & (1 << Twice)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m2;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m2) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  Twice = (m1 << 1) + 1;
  if ((TmpState & (1 << Twice))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m1;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m1) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  return this->FindStateIndex(TmpState, NewMaxMomentum);
}

// apply au^+_m1 au^+_m2 ad_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpin::AudAudAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (2 << (n1 << 1))) == 0) || 
      ((State & (2 << (n2 << 1))) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  int TmpState = State;
  int Mask;
  int SignIndex = (NewMaxMomentum - n2) << 1;
  int Twice = (n2 << 1) + 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient = this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient = this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n2)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (n1 << 1) + 1;
  SignIndex = (NewMaxMomentum - n1) << 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient *= this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n1)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (m2 << 1);
  if ((TmpState & (1 << Twice)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m2;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m2) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  Twice = (m1 << 1);
  if ((TmpState & (1 << Twice))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m1;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m1) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  return this->FindStateIndex(TmpState, NewMaxMomentum);
}

// apply au^+_m1 au^+_m2 au_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpin::AudAudAuAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (1 << (n1 << 1))) == 0) || 
      ((State & (2 << (n2 << 1))) == 0) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  int TmpState = State;
  int Mask;
  int SignIndex = (NewMaxMomentum - n2) << 1;
  int Twice = (n2 << 1) + 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient = this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient = this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n2)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (n1 << 1);
  SignIndex = (NewMaxMomentum - n1) << 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient *= this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n1)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (m2 << 1);
  if ((TmpState & (1 << Twice)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m2;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m2) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  Twice = (m1 << 1);
  if ((TmpState & (1 << Twice))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m1;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m1) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  return this->FindStateIndex(TmpState, NewMaxMomentum);
}

// apply ad^+_m1 au^+_m2 au_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpin::AddAudAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (1 << (n1 << 1))) == 0) || 
      ((State & (1 << (n2 << 1))) == 0) || (n1 == n2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  int TmpState = State;
  int Mask;
  int SignIndex = (NewMaxMomentum - n2) << 1;
  int Twice = (n2 << 1);
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient = this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient = this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n2)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (n1 << 1);
  SignIndex = (NewMaxMomentum - n1) << 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient *= this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n1)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (m2 << 1);
  if ((TmpState & (1 << Twice)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m2;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m2) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  Twice = (m1 << 1) + 1;
  if ((TmpState & (1 << Twice))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m1;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m1) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  return this->FindStateIndex(TmpState, NewMaxMomentum);
}

// apply ad^+_m1 ad^+_m2 ad_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpin::AddAddAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (2 << (n1 << 1))) == 0) || 
      ((State & (1 << (n2 << 1))) == 0) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  int TmpState = State;
  int Mask;
  int SignIndex = (NewMaxMomentum - n2) << 1;
  int Twice = (n2 << 1);
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient = this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient = this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n2)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (n1 << 1) + 1;
  SignIndex = (NewMaxMomentum - n1) << 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient *= this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n1)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (m2 << 1) + 1;
  if ((TmpState & (1 << Twice)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m2;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m2) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  Twice = (m1 << 1) + 1;
  if ((TmpState & (1 << Twice))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m1;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m1) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  return this->FindStateIndex(TmpState, NewMaxMomentum);
}

// apply au^+_m1 ad^+_m2 ad_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpin::AudAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (2 << (n1 << 1))) == 0) || 
      ((State & (2 << (n2 << 1))) == 0) || (n1 == n2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  int TmpState = State;
  int Mask;
  int SignIndex = (NewMaxMomentum - n2) << 1;
  int Twice = (n2 << 1) + 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient = this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient = this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n2)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (n1 << 1) + 1;
  SignIndex = (NewMaxMomentum - n1) << 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient *= this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n1)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (m2 << 1) + 1;
  if ((TmpState & (1 << Twice)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m2;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m2) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  Twice = (m1 << 1);
  if ((TmpState & (1 << Twice))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m1;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m1) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  return this->FindStateIndex(TmpState, NewMaxMomentum);
}

// apply au^+_m1 ad^+_m2 au_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpin::AudAddAuAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned int State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (1 << (n1 << 1))) == 0) || 
      ((State & (2 << (n2 << 1))) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  int TmpState = State;
  int Mask;
  int SignIndex = (NewMaxMomentum - n2) << 1;
  int Twice = (n2 << 1) + 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient = this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient = this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n2)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (n1 << 1);
  SignIndex = (NewMaxMomentum - n1) << 1;
  if (SignIndex < 16)
    {
      Mask = (TmpState >> Twice);
      coefficient *= this->SignLookUpTable[Mask];
    }
  else
    {
      Mask = (TmpState >> Twice) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
      Mask = (TmpState >> (Twice + 16)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 32)) & 0xffff;
      coefficient *= this->SignLookUpTable[Mask];
      Mask = (TmpState >> (Twice + 48));
      coefficient *= this->SignLookUpTable[Mask];
#else
      Mask = (TmpState >> (Twice + 16));
      coefficient *= this->SignLookUpTable[Mask];
#endif
    }
  TmpState &= ~(0x1 << Twice);
  if (NewMaxMomentum == n1)
    while ((TmpState >> (NewMaxMomentum << 1)) == 0)
      --NewMaxMomentum;
  Twice = (m2 << 1) + 1;
  if ((TmpState & (1 << Twice)) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m2;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m2) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  Twice = (m1 << 1);
  if ((TmpState & (1 << Twice))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      TmpState |= (0x1 << Twice);
      NewMaxMomentum = m1;
    }
  else
    {
      SignIndex = (NewMaxMomentum - m1) << 1;
      if (SignIndex < 16)
	{
	  Mask = (TmpState >> Twice);
	  coefficient *= this->SignLookUpTable[Mask];
	}
      else
	{
	  Mask = (TmpState >> Twice) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
#ifdef  __64_BITS__
	  Mask = (TmpState >> (Twice + 16)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 32)) & 0xffff;
	  coefficient *= this->SignLookUpTable[Mask];
	  Mask = (TmpState >> (Twice + 48));
	  coefficient *= this->SignLookUpTable[Mask];
#else
	  Mask = (TmpState >> (Twice + 16));
	  coefficient *= this->SignLookUpTable[Mask];
#endif
	}
      TmpState |= (0x1 << Twice);
    }
  return this->FindStateIndex(TmpState, NewMaxMomentum);
}

// find state index
//
// stateDescription = unsigned integer describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnTorusWithSpin::FindStateIndex(unsigned int stateDescription, int maxMomentum)
{
//  int Pos = 0;
  int Pos = this->LookUpTable[maxMomentum][stateDescription >> this->LookUpTableShift[maxMomentum]];
  while ((Pos < this->HilbertSpaceDimension) && (this->StateDescription[Pos] != stateDescription))
    ++Pos;
  return Pos;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnTorusWithSpin::PrintState (ostream& Str, int state)
{
  int TmpState = this->StateDescription[state];
  for (int i = 0; i < this->MaxMomentum; ++i)
    {
      switch ((TmpState >> (i * 2)) & 0x3)
	{
	case 0x01:
	  Str << " + ";
	  break;
	case 0x02:
	  Str << " - ";
	  break;
	case 0x03:
	  Str << "+- ";
	  break;
	default:
	  Str << " 0 "; 
	}
    }
//  Str << " key = " << this->Keys[state] << " maxMomentum position = " << this->MaxMomentumPosition[Max * (this->NbrFermions + 1) + TmpState[Max]]
  Str << " position = " << FindStateIndex(TmpState, this->StateMaxMomentum[state]) << " P = " << this->GetMomentumValue(state);
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int FermionOnTorusWithSpin::GenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos)
{
  if ((nbrFermions == 0) || (nbrFermions > (currentMaxMomentum + 1)) || (currentMaxMomentum < 0))
    return pos;
  if (nbrFermions == 1)
    {
      for (int i = currentMaxMomentum; i >= 0; --i)
	{
	  this->StateDescription[pos] = 1 << i;
	  this->StateMaxMomentum[pos] = maxMomentum >> 1;
	  ++pos;
	}
      return pos;
    }
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentMaxMomentum, pos);
  unsigned int Mask = 1 << currentMaxMomentum;
  for (int i = pos; i < TmpPos; ++i)
    this->StateDescription[i] |= Mask;
  if (maxMomentum == currentMaxMomentum)
    return this->GenerateStates(nbrFermions, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, TmpPos);
  else
    return this->GenerateStates(nbrFermions, maxMomentum, ReducedCurrentMaxMomentum, TmpPos);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// currentTotalSpinMomentum = total spin momemtum of the fermions that have already been placed
// return value = position from which new states have to be stored

int FermionOnTorusWithSpin::GenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos, int currentTotalSpinMomentum)
{
  if ((nbrFermions == 0) || (nbrFermions > (currentMaxMomentum + 1)) || (currentMaxMomentum < 0) || 
      (fabs(this->TotalSpinMomentum - currentTotalSpinMomentum) > nbrFermions))
    return pos;
  if (nbrFermions == 1)
    {
      if (this->TotalSpinMomentum > currentTotalSpinMomentum)
	{
	  currentMaxMomentum &= ~0x1;
	  for (; currentMaxMomentum >= 0; currentMaxMomentum -= 2)
	    {
	      this->StateDescription[pos] = 1 << currentMaxMomentum;
	      this->StateMaxMomentum[pos] = maxMomentum >> 1;
	      ++pos;
	    }
	}
      else
	{
	  if ((currentMaxMomentum & 1) == 0)
	    --currentMaxMomentum;
	  for (; currentMaxMomentum >= 0; currentMaxMomentum -= 2)
	    {
	      this->StateDescription[pos] = 1 << currentMaxMomentum;
	      this->StateMaxMomentum[pos] = maxMomentum >> 1;
	      ++pos;
	    }
	}
      return pos;
    }
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int ReducedCurrentTotalSpinMomentum = currentTotalSpinMomentum;
  if ((currentMaxMomentum & 1) == 0)
    ReducedCurrentTotalSpinMomentum += 1;
  else
    ReducedCurrentTotalSpinMomentum -= 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentMaxMomentum, pos, ReducedCurrentTotalSpinMomentum);
  unsigned int Mask = 1 << currentMaxMomentum;
  for (int i = pos; i < TmpPos; ++i)
    this->StateDescription[i] |= Mask;
  if (maxMomentum == currentMaxMomentum)
    return this->GenerateStates(nbrFermions, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentTotalSpinMomentum);
  else
    return this->GenerateStates(nbrFermions, maxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentTotalSpinMomentum);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// currentTotalSpinMomentum = total spin momemtum of the fermions that have already been placed
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

int FermionOnTorusWithSpin::GenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos, int currentTotalSpinMomentum, int currentMomentum)
{
  if ((nbrFermions == 0) || (nbrFermions > (currentMaxMomentum + 1)) || (currentMaxMomentum < 0) || 
      (fabs(this->TotalSpinMomentum - currentTotalSpinMomentum) > nbrFermions))
    return pos;
  if (nbrFermions == 1)
    {
      if (this->TotalSpinMomentum > currentTotalSpinMomentum)
	{
	  int i = this->MomentumConstraint - (currentMomentum % this->MaxMomentum);
	  if (i < 0)
	    i += this->MaxMomentum;
	  i <<= 1;
	  currentMaxMomentum &= ~0x1;
	  for (; i <= currentMaxMomentum; i += (this->MaxMomentum << 1))
	    {
	      this->StateDescription[pos] = 1 << i;
	      this->StateMaxMomentum[pos] = maxMomentum >> 1;
	      ++pos;
	    }
	}
      else
	{
	  int i = this->MomentumConstraint - (currentMomentum % this->MaxMomentum);
	  if (i < 0)
	    i += this->MaxMomentum;
	  i <<= 1;
	  ++i;
	  if ((currentMaxMomentum & 1) == 0)
	    {
	      --currentMaxMomentum;
	    }
	  for (; i <= currentMaxMomentum; i += (this->MaxMomentum << 1))
	    {
	      this->StateDescription[pos] = 1 << i;
	      this->StateMaxMomentum[pos] = maxMomentum >> 1;
	      ++pos;
	    }
	}
      return pos;
    }
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int ReducedCurrentTotalSpinMomentum = currentTotalSpinMomentum;
  if ((currentMaxMomentum & 1) == 0)
    ReducedCurrentTotalSpinMomentum += 1;
  else
    ReducedCurrentTotalSpinMomentum -= 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentMaxMomentum, pos, ReducedCurrentTotalSpinMomentum, currentMomentum + (currentMaxMomentum >> 1));
  unsigned int Mask = 1 << currentMaxMomentum;
  for (int i = pos; i < TmpPos; ++i)
    this->StateDescription[i] |= Mask;
  if (maxMomentum == currentMaxMomentum)
    return this->GenerateStates(nbrFermions, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentTotalSpinMomentum, currentMomentum);
  else
    return this->GenerateStates(nbrFermions, maxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentTotalSpinMomentum, currentMomentum);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnTorusWithSpin::GenerateLookUpTable(int memory)
{
  // evaluate look-up table size
  memory /= (4 * this->NbrLzValue);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 2;
      this->MaximumLookUpShift += 2;
    }
  if (this->MaximumLookUpShift > (2 * this->NbrLzValue))
    {
      this->MaximumLookUpShift = 2 * this->NbrLzValue;
      this->LookUpTableMemorySize = 1 << (this->MaximumLookUpShift + 1);
    }
  else
    {
      this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;
      --this->MaximumLookUpShift;
    }

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrLzValue];
  this->LookUpTableShift = new int [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize];
  int CurrentMaxMomentum = this->StateMaxMomentum[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
  if ((2 * CurrentMaxMomentum) < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentMaxMomentum] = 0;
  else
    this->LookUpTableShift[CurrentMaxMomentum] = 2 * (CurrentMaxMomentum + 1) - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
  unsigned int CurrentLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  unsigned int TmpLookUpTableValue;
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentMaxMomentum != this->StateMaxMomentum[i])
	{
 	  CurrentMaxMomentum = this->StateMaxMomentum[i];
	  TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
	  if ((2 * CurrentMaxMomentum) < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentMaxMomentum] = 0;
	  else
	    this->LookUpTableShift[CurrentMaxMomentum] = 2 * (CurrentMaxMomentum + 1) - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
	  CurrentLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      CurrentLookUpTableValue = TmpLookUpTableValue;
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  
  // look-up tables for evaluating sign when applying creation/annihilation operators
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
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
	this->SignLookUpTable[j] = -1.0;
      else
	this->SignLookUpTable[j] = 1.0;
    }
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// return value = Hilbert space dimension

int FermionOnTorusWithSpin::EvaluateHilbertSpaceDimension(int nbrFermions, int maxMomentum)
{
  FactorialCoefficient Dimension; 
  Dimension.PartialFactorialMultiply(2 * maxMomentum - nbrFermions + 1, 2 * maxMomentum); 
  Dimension.FactorialDivide(nbrFermions);
  return (Dimension.GetIntegerValue());
}

// evaluate Hilbert space dimension for a given total spin momentum
//
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// spinMomemtum = twice the total spin momentum
// return value = Hilbert space dimension

int FermionOnTorusWithSpin::EvaluateHilbertSpaceDimension(int nbrFermions, int maxMomentum, int spinMomemtum)
{
  FactorialCoefficient Dimension; 
  int NbrUpFermion = (spinMomemtum + nbrFermions) >> 1;
  Dimension.PartialFactorialMultiply(maxMomentum - NbrUpFermion + 1, maxMomentum); 
  Dimension.FactorialDivide(NbrUpFermion);
  NbrUpFermion = nbrFermions - NbrUpFermion;
  Dimension.PartialFactorialMultiply(maxMomentum - NbrUpFermion + 1, maxMomentum); 
  Dimension.FactorialDivide(NbrUpFermion);
  return (Dimension.GetIntegerValue());
}

