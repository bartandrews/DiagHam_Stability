////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of doubled spin 0 +1/2 chain with translations           //
//                                                                            //
//                        last modification : 21/01/2016                      //
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


#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <math.h>

using std::cout;
using std::endl;


// default constructor
//

DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry::DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry () 
{
  this->ZEigenvalueBra = 0;
  this->ZEigenvalueKet = 0;
}

 

// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry::DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry (int chainLength, int diffSz,  int zEigenvalueBra, int zEigenvalueKet, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->ZEigenvalueBra=zEigenvalueBra;
  this->ZEigenvalueKet=zEigenvalueKet;
  this->ComplementaryStateShift = 2*(this->ChainLength - 1);
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->PowerD = new int [this->ChainLength];
  this->PowerD[0]=1;
  for(int i = 1 ; i <this->ChainLength;i++)
    this->PowerD[i]=this->PowerD[i-1]*9;
  
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, this->DiffSz);
  this->ShiftNegativeDiffSz =   this->LargeHilbertSpaceDimension;
  if (this->DiffSz !=0 )
    this->LargeHilbertSpaceDimension += this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, -this->DiffSz);
  
  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, this->DiffSz, 0l);
  if (this->DiffSz != 0)
    TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, -this->DiffSz, TmpHilbertSpaceDimension);
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in DoubledSpin0_1_2_ChainWithTranslations!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  
  this->LargeHilbertSpaceDimension = 0l;

  unsigned long DicardFlag = ~0x0ul;
  
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      if ((this->ComputeZValueBra(this->ChainDescription[i]) == this->ZEigenvalueBra) &&(this->ComputeZValueKet(this->ChainDescription[i]) == this->ZEigenvalueKet) )
	{
	  ++this->LargeHilbertSpaceDimension;
	}
      else
	{
	  this->ChainDescription[i] = DicardFlag;
	}
    }
  
  
  unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];

  
  this->LargeHilbertSpaceDimension = 0l;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      if ( this->ChainDescription[i] != DicardFlag)
	{
	  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->ChainDescription[i];
	  ++this->LargeHilbertSpaceDimension;
	}
    }
  
  delete[]  this->ChainDescription;

  this->ChainDescription = TmpStateDescription;
  
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;  
  this->LookUpTable =0;
  SortArrayDownOrdering(this->ChainDescription ,this->HilbertSpaceDimension);
  
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
    }
  
  this->RescalingFactors = 0;

}
 


// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry::DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry (int chainLength, int momentum, int diffSz, int zEigenvalueBra, int zEigenvalueKet, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->Momentum = momentum;
  this->ZEigenvalueBra=zEigenvalueBra;
  this->ZEigenvalueKet=zEigenvalueKet;
  this->ComplementaryStateShift = 2*(this->ChainLength - 1);
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->PowerD = new int [this->ChainLength];
  this->PowerD[0]=1;
  for(int i = 1 ; i <this->ChainLength;i++)
    this->PowerD[i]=this->PowerD[i-1]*9;
  
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, this->DiffSz);
  this->ShiftNegativeDiffSz =   this->LargeHilbertSpaceDimension;
  if (this->DiffSz !=0 )
    this->LargeHilbertSpaceDimension += this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, -this->DiffSz);

  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
    
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, this->DiffSz, 0l);
  if (this->DiffSz != 0)
    TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, -this->DiffSz, TmpHilbertSpaceDimension);
  
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2Spin!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = 0l;
  unsigned long TmpCanonicalState;

  int NbrTranslation;
  int CurrentNbrStateInOrbit;
  unsigned long DicardFlag = ~0x0ul;
  unsigned long TmpState;
  
  this->CreatePrecalculationTable();

  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      TmpState = this->ChainDescription[i];
      this->FindCanonicalForm(TmpState,TmpCanonicalState, NbrTranslation);

      if (TmpState  == TmpCanonicalState)
	{
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpCanonicalState);
	  if ((this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true)&&(this->ComputeZValueBra(TmpCanonicalState) == this->ZEigenvalueBra) &&(this->ComputeZValueKet(TmpCanonicalState) == this->ZEigenvalueKet) )
	    {
	      ++this->LargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->ChainDescription[i] = DicardFlag;
	    }
	}
      else
	{
	  this->ChainDescription[i] = DicardFlag;
	}
    }
  
  unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  
  this->NbrStateInOrbit = new int [this->LargeHilbertSpaceDimension];

  this->LargeHilbertSpaceDimension = 0l;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      if ( this->ChainDescription[i] != DicardFlag)
	{
	  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->ChainDescription[i];
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(this->ChainDescription[i]);
	  this->NbrStateInOrbit[this->LargeHilbertSpaceDimension] = CurrentNbrStateInOrbit;
	  ++this->LargeHilbertSpaceDimension;
	}
    }
  
  delete[]  this->ChainDescription;
  this->ChainDescription = TmpStateDescription;
  
  
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;  
  this->LookUpTable =0;

  if (this->HilbertSpaceDimension > 0)
    this->GenerateLookUpTable(memorySize);
  this->EvaluateExponentialFactors();
  
}


// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry::DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry (const DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry & chain)
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescription = chain.ChainDescription;
      this->DiffSz = chain.DiffSz;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->UniqueStateDescriptionBra = chain.UniqueStateDescriptionBra;
      this->UniqueStateDescriptionSubArraySizeBra = chain.UniqueStateDescriptionSubArraySizeBra;
      this->NbrUniqueStateDescriptionBra = chain.NbrUniqueStateDescriptionBra;
      this->FirstIndexUniqueStateDescriptionBra = chain.FirstIndexUniqueStateDescriptionBra;
      this->ShiftNegativeDiffSz = chain.ShiftNegativeDiffSz;
      this->BraShiftNegativeSz =  chain.BraShiftNegativeSz;
      this->PowerD = chain.PowerD;
      this->TranslationPhase=chain.TranslationPhase;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainDescriptionBra = 0;
      this->ChainDescriptionKet = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->DiffSz = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->ShiftNegativeDiffSz =0;
      this->BraShiftNegativeSz =0;
      this->PowerD = 0;
      this->TranslationPhase=0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry::~DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry & DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry::operator = (const DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry & chain)
{
  DoubledSpin0_1_2_ChainWithTranslationsStaggered::operator =(chain);
  this->ZEigenvalueBra=chain.ZEigenvalueBra;
  this->ZEigenvalueKet=chain.ZEigenvalueKet;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry::Clone()
{
  return new DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry (*this);
}


