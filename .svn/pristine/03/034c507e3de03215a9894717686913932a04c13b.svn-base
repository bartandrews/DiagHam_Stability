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
//                        last modification : 25/02/2016                      //
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


#include "HilbertSpace/Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;


// default constructor
//

Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers:: Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers ()
{
  this->SubLatticeDifference = 0;
}


// constructor for Hilbert space corresponding no constaint on Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers:: Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers (int chainLength, int subLatticeDifference, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->ComplementaryStateShift = 2*(this->ChainLength - 1);
  this->SubLatticeDifference =  subLatticeDifference;
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  memorySize/=8;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;
  
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1);
  this->ShiftNegativeDiffSz = this->LargeHilbertSpaceDimension;
  
  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, 0l);
  
  SortArrayDownOrdering(this->ChainDescription ,TmpHilbertSpaceDimension);
  
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in Spin0_1_2_ChainWithTranslations!" << endl;
      exit(1);
    } 
  
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = 0l;
  unsigned long DicardFlag = ~0x0ul;
  int SublatticeNumber;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      this->ComputeDifferenceSubLatticeNumberZero(this->ChainDescription[i],  SublatticeNumber);
      if (  SublatticeNumber == this->SubLatticeDifference  )
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
  
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable();
    }
  this->RescalingFactors = 0;
  
  for(int i=0; i < this->LargeHilbertSpaceDimension; i++)
    {
      if ( i !=  this->FindStateIndex(this->ChainDescription[i]) )
	{
	  cout <<"Problem in Find state index "<< i<< " "<< this->FindStateIndex(this->ChainDescription[i])<<endl;
	}
    }
  cout <<"Hilbert Space dimension = "<< this->GetHilbertSpaceDimension()<<endl;  
}


// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers:: Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers (int chainLength, int diffSz, int subLatticeDifference, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->ComplementaryStateShift = 2*(this->ChainLength - 1);
  this->SubLatticeDifference =  subLatticeDifference;
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  memorySize/=8;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;
  
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->DiffSz);
  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->DiffSz, 0l);
  
  SortArrayDownOrdering(this->ChainDescription ,TmpHilbertSpaceDimension);
  
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in Spin0_1_2_ChainWithTranslations!" << endl;
      exit(1);
    } 
  
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  
  this->LargeHilbertSpaceDimension = 0l;
  unsigned long DicardFlag = ~0x0ul;
  int SublatticeNumber;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      this->ComputeDifferenceSubLatticeNumberZero(this->ChainDescription[i],  SublatticeNumber);
      if (  SublatticeNumber == this->SubLatticeDifference  )
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
  
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable();
    }
  
  this->RescalingFactors = 0;
  
  for(int i=0; i < this->LargeHilbertSpaceDimension; i++)
    {
      if ( i !=  this->FindStateIndex(this->ChainDescription[i]) )
	{
	  cout <<"Problem in Find state index "<< i<< " "<< this->FindStateIndex(this->ChainDescription[i])<<endl;
	}
    }
}
 

// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers:: Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers (int chainLength, int momentum, int translationStep, int diffSz,  int subLatticeDifference, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->Momentum = momentum;
  this->SubLatticeDifference =  subLatticeDifference;
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;
  
  this->StateMask = (0x1ul << (translationStep<<1)) - 1ul;
  this->StateShift = translationStep<<1;
  this->MaxXMomentum = this->ChainLength/ translationStep;
  this->ComplementaryStateShift = 2*(this->ChainLength - translationStep);

  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->DiffSz);

  
  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->DiffSz, 0l);

  SortArrayDownOrdering(this->ChainDescription ,TmpHilbertSpaceDimension);
  
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in  Spin0_1_2_ChainWithTranslations!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = 0l;
  unsigned long TmpCanonicalState;
  int NbrTranslation;
  int CurrentNbrStateInOrbit;
  unsigned long DicardFlag = ~0x0ul;
  unsigned long TmpState;
  this->CreatePrecalculationTable();
  int TmpSublatticeNumber;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      TmpState = this->ChainDescription[i];
      
      this->ComputeDifferenceSubLatticeNumberZero(this->ChainDescription[i],  TmpSublatticeNumber);
      if (  TmpSublatticeNumber == this->SubLatticeDifference  )
	{
	  this->FindCanonicalForm(TmpState,TmpCanonicalState,NbrTranslation);
	  
	  
	  if (TmpState  == TmpCanonicalState)
	    {
	      CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpCanonicalState);
	      if (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true)
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
  this->ShiftNegativeDiffSz++;
  delete[]  this->ChainDescription;
  this->ChainDescription = TmpStateDescription;

  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;  
  this->LookUpTable =0;

  if (this->HilbertSpaceDimension > 0)
    this->GenerateLookUpTable();
}


// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers::Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers (const Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers & chain) :   Spin0_1_2_ChainWithTranslations(chain)
{
  this->SubLatticeDifference =  chain.SubLatticeDifference;
}

// destructor
//

Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers::~Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers & Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers::operator = (const Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers & chain)
{
  Spin0_1_2_ChainWithTranslations::operator =(chain);
  this->SubLatticeDifference =  chain.SubLatticeDifference;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers::Clone()
{
  return new Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers (*this);
}

