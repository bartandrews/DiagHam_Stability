////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            /
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

#include "HilbertSpace/Spin0_1_2_ChainWithTranslations.h"
#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslations.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/ArrayTools.h"
#include <iostream>
#include <math.h>

using std::cout;
using std::endl;


// default constructor
//

DoubledSpin0_1_2_ChainWithTranslations::DoubledSpin0_1_2_ChainWithTranslations () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->ChainLength = 0;
  this->Momentum = 0;
  this->ComplementaryStateShift = 0;
  this->DiffSz = 0;
  this->FixedSpinProjectionFlag = false;
  this->CompatibilityWithMomentum = 0;
  this->RescalingFactors = 0;
  this->NbrStateInOrbit = 0;
  this->ShiftNegativeDiffSz = 0;
  this->BraShiftNegativeSz = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  this->PowerD = 0;
  this->TranslationPhase=0;
}

// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin0_1_2_ChainWithTranslations::DoubledSpin0_1_2_ChainWithTranslations (int chainLength, int diffSz, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->ComplementaryStateShift = 2*(this->ChainLength - 1);
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;
  this->TranslationPhase=0;

  this->PowerD = new int [this->ChainLength];
  this->PowerD[0]=1;
  for(int i = 1 ; i <this->ChainLength;i++)
    this->PowerD[i]=this->PowerD[i-1]*9;
  
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, this->DiffSz);
  this->ShiftNegativeDiffSz = this->LargeHilbertSpaceDimension;

/*  if (this->DiffSz !=0 )
    this->LargeHilbertSpaceDimension += this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, -this->DiffSz);*/

  this->ChainDescriptionBra = 0;
  this->ChainDescriptionKet = 0;
  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, this->DiffSz, 0l);

/*  if (this->DiffSz != 0)
    TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, -this->DiffSz, TmpHilbertSpaceDimension);*/

  SortArrayDownOrdering(this->ChainDescription ,TmpHilbertSpaceDimension);

  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in DoubledSpin0_1_2_ChainWithTranslations!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
    }
  this->RescalingFactors = 0;
  int MemoryCost = sizeof(unsigned long) * this->LargeHilbertSpaceDimension + this->ChainLength *  sizeof(int);
  cout <<" Memory cost of the Hilbert Space " <<MemoryCost <<endl;
}
 


// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin0_1_2_ChainWithTranslations::DoubledSpin0_1_2_ChainWithTranslations (int chainLength, int momentum,  int translationStep, int diffSz, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;


  this->MaxXMomentum = this->ChainLength/ translationStep;
  this->ComplementaryStateShift = 2*(this->ChainLength - translationStep);
  this->Momentum = momentum %  this->MaxXMomentum;
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

/*  if (this->DiffSz !=0 )
    this->LargeHilbertSpaceDimension += this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, -this->DiffSz);*/


  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, this->DiffSz, 0l);

/*  if (this->DiffSz != 0)
    TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, -this->DiffSz, TmpHilbertSpaceDimension);*/
  
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
  SortArrayDownOrdering(this->ChainDescription ,TmpHilbertSpaceDimension);
  
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      TmpState = this->ChainDescription[i];
      this->FindCanonicalForm(TmpState,TmpCanonicalState, NbrTranslation);

      if ((TmpState  == TmpCanonicalState))
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

  unsigned long MemoryCost = (sizeof(unsigned long) + sizeof(int)) * this->LargeHilbertSpaceDimension + this->ChainLength *  sizeof(int) + this->MaxXMomentum*sizeof(Complex);
  cout <<" Memory cost of the Hilbert Space " <<MemoryCost <<endl;
}


// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

DoubledSpin0_1_2_ChainWithTranslations::DoubledSpin0_1_2_ChainWithTranslations (const DoubledSpin0_1_2_ChainWithTranslations & chain)
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
      this->ChainDescriptionBra = chain.ChainDescriptionBra;
      this->ChainDescriptionKet = chain.ChainDescriptionKet;
      this->DiffSz = chain.DiffSz;
      this->Momentum = chain.Momentum;
      this->MaxXMomentum = chain.MaxXMomentum;
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
      this->ChainDescriptionBra = 0;
      this->ChainDescriptionKet = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->DiffSz = 0;
      this->ChainDescription = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->ShiftNegativeDiffSz =0;
      this->BraShiftNegativeSz =0;
      this->PowerD = 0;
      this->TranslationPhase=0;
      this->MaxXMomentum = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

DoubledSpin0_1_2_ChainWithTranslations::~DoubledSpin0_1_2_ChainWithTranslations () 
{
  delete [] this->PowerD;
  delete [] this->ChainDescription;
  delete [] TranslationPhase;
  this->LargeHilbertSpaceDimension = 0;
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

DoubledSpin0_1_2_ChainWithTranslations & DoubledSpin0_1_2_ChainWithTranslations::operator = (const DoubledSpin0_1_2_ChainWithTranslations & chain)
{
  AbstractDoubledSpinChainWithTranslations::operator =(chain);
  this->ShiftNegativeDiffSz = chain.ShiftNegativeDiffSz;
  this->BraShiftNegativeSz = chain.BraShiftNegativeSz;
  this->PowerD = chain.PowerD;
  this->ChainDescription = chain.ChainDescription;
  this->TranslationPhase = chain.TranslationPhase;
  this->MaxXMomentum = chain.MaxXMomentum;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* DoubledSpin0_1_2_ChainWithTranslations::Clone()
{
  return new DoubledSpin0_1_2_ChainWithTranslations (*this);
}

// return value of twice spin projection on (Oz) for a given state
//
// stateDescription = state to which the spin projection has to be evaluated
// return value = twice spin projection on (Oz)

inline int DoubledSpin0_1_2_ChainWithTranslations::GetTotalSz (unsigned long stateDescriptionBra, unsigned long stateDescriptionKet)
{
  int TmpSz = 0;
  int Sign;
  for (int i = 0; i < this->ChainLength; i++)
    {
      switch (stateDescriptionBra & 0x3ul)
	{
	case 0x1:
	  TmpSz += 1;
	  break;
	case 0x0:
	  TmpSz -= 1;
	  break;
	}
      stateDescriptionBra >>= 2;
      switch (stateDescriptionKet & 0x3ul)
	{
	case 0x1:
	  TmpSz -= 1;
	  break;
	case 0x0:
	  TmpSz += 1;
	  break;
	}
      stateDescriptionKet >>= 2;
    }
  return TmpSz;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& DoubledSpin0_1_2_ChainWithTranslations::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned int tmpBra,tmpKet;
  int BraNumber, KetNumber;
  unsigned long StateDescription = this->ChainDescription[state];  

  this->ComputeDiffereenceSubLatticeNumberZero(StateDescription, BraNumber, KetNumber);

  Str << this->FindStateIndex(StateDescription) << " : " << StateDescription<< " : ";
  for (int j = this->ChainLength; j >0; j--)
    {
      this->GetBraAndKetIndicesFromCommonIndex(tmpBra,tmpKet, StateDescription%9);
      StateDescription/=9;
      Str << "(";
      if (tmpBra == 0)
	Str << "d ";
      else
	if (tmpBra == 1)
	  Str << "u ";
	else
	  Str << "0 ";
      Str << ",";
      if (tmpKet == 0)
	Str << "d ";
      else
	if (tmpKet == 1)
	  Str << "u ";
	else
	  Str << "0 ";
      Str << ") ";
    }
  Str<<" SublatticeZeroNumberDifference = (" << BraNumber<<"," << KetNumber<< ")";
  return Str;
}


// generate all states corresponding to the constraints
// 
// lengthBra = length of the chain to be decided for bra spins
// lengthBra = length of the chain to be decided for ket spins
// diffSz = difference of spin projection between bra and ket chain
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
long DoubledSpin0_1_2_ChainWithTranslations::GenerateStates(int lengthBra, int lengthKet, int diffSz, long pos)
{
  if ((lengthKet == 0) && (lengthBra == 0))
    {
      if (diffSz == 0) 
	{
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(2,2);
	  pos++;
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(1,1);
	  pos++;
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(0,0);
	  pos++;
	  return pos;
	}
      if (diffSz == 1) 
	{
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(1,2);
	  pos++;
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(2,0);
	  pos++;
	  return pos;
	}
      if (diffSz == -1) 
	{
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(2,1);
	  pos++;
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(0,2);
	  pos++;
	  return pos;
	}

      if (diffSz == 2) 
	{
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(1,0);
	  pos++;
	  return pos;
	}
      if (diffSz == -2) 
	{
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(0,1);
	  pos++;
	  return pos;
	}
      return pos;
    }

  long TmpPos;
  unsigned long MaskBra;
  unsigned long MaskKet;
  
  if(lengthKet > 0)
    { 

      TmpPos = this->GenerateStates(lengthBra,lengthKet-1, diffSz+1, pos); 
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescription[pos] += (GetCommonIndexFromBraAndKetIndices(0,1)*this->PowerD[lengthKet]);
	}
      TmpPos = this->GenerateStates(lengthBra,lengthKet-1, diffSz, pos); 
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescription[pos] += GetCommonIndexFromBraAndKetIndices(0,2)*this->PowerD[lengthKet];
	}
      TmpPos = this->GenerateStates(lengthBra,lengthKet-1, diffSz-1, pos); 
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescription[pos] += GetCommonIndexFromBraAndKetIndices(0,0)*this->PowerD[lengthKet];
	}
      return pos;
     }
  
  if (lengthKet == 0)
    {

      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz-1, pos); 
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescription[pos] += GetCommonIndexFromBraAndKetIndices(1,0)*this->PowerD[lengthBra];
	}
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz, pos); 
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescription[pos] += GetCommonIndexFromBraAndKetIndices(2,0)*this->PowerD[lengthBra];
	} 
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz+1, pos); 
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescription[pos] += GetCommonIndexFromBraAndKetIndices(0,0)*this->PowerD[lengthBra];
	}
      return pos;
    }
  return pos;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrNUp = number of particles with quantum number up
// nbrNDown = number of particles with quantum number down
// return value = Hilbert space dimension

long DoubledSpin0_1_2_ChainWithTranslations::ShiftedEvaluateHilbertSpaceDimension(int lengthBra, int lengthKet, int diffSz)
{
  if ((lengthBra < 0) || (lengthKet < 0))
    return 0;
  
  if ((lengthBra == 0) && (lengthKet == 0))
    {
      if (diffSz == 0) 
	{
	  return 3;
	}
      if (diffSz == 1) 
	{
	  return 2;
	}
      if (diffSz == -1) 
	{
	  return 2;
	}

      if (diffSz == 2) 
	{
	  return 1;
	}
      if (diffSz == -2) 
	{
	  return 1;
	}
      return 0;
    }  
  long Tmp=0;
  
  if (lengthBra == 0)
    {
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra,lengthKet-1, diffSz-1);   
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra,lengthKet-1, diffSz); 
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra,lengthKet-1, diffSz+1); 
      return Tmp;
    }

  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra-1,lengthKet, diffSz+1); 
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra-1,lengthKet, diffSz); 
  Tmp +=  this->ShiftedEvaluateHilbertSpaceDimension(lengthBra-1,lengthKet, diffSz-1);
  return Tmp;
}



// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void DoubledSpin0_1_2_ChainWithTranslations::GenerateLookUpTable(unsigned long memory)
{  
  this->LookUpTable = new long [82];

  int CurrentBeginning = this->ChainDescription[0] / this->PowerD[this->ChainLength-2];
  for (int i =CurrentBeginning+1; i <82;i++)
    {
      this->LookUpTable[i]=0;
    }
  for(int i =0; i < this->HilbertSpaceDimension; i++)
    {
      if( this->ChainDescription[i]/this->PowerD[this->ChainLength-2] >= CurrentBeginning)
	{

	}
      else
	{
	  for(int p = this->ChainDescription[i] / this->PowerD[this->ChainLength-2] + 1 ;p <= CurrentBeginning;p++)
	    this->LookUpTable[p] = i;
	  CurrentBeginning = this->ChainDescription[i] / this->PowerD[this->ChainLength-2];
	}
      
    }
  for(int p = 0 ;p <= CurrentBeginning;p++)
  this->LookUpTable[p] = this->HilbertSpaceDimension;
}
 

int DoubledSpin0_1_2_ChainWithTranslations::FindStateIndex(unsigned long stateDescription)
{
  int PosMax = this->LookUpTable[stateDescription/this->PowerD[this->ChainLength-2]];
  int PosMin = this->LookUpTable[stateDescription/this->PowerD[this->ChainLength-2]+1];
  int PosMid = (PosMin + PosMax) >>1;
  while( (PosMax - PosMin) >1 ) 
    {
      if (this->ChainDescription[PosMid]  >  stateDescription)
	{
	  PosMin = PosMid;
	}
      else
	{
	  if (this->ChainDescription[PosMid]  <  stateDescription)
	    {
	      PosMax = PosMid;
	    }
	  else
	    return PosMid;
	}
      
      PosMid = (PosMin + PosMax) >>1;
    }
  
  if (this->ChainDescription[PosMid] == stateDescription)
    return PosMid;
  else
    {
      if ((PosMax == this->HilbertSpaceDimension) ||(this->ChainDescription[PosMax] != stateDescription))
	return this->HilbertSpaceDimension;
      else
	return PosMax;
    }
}
 
// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix DoubledSpin0_1_2_ChainWithTranslations::EvaluatePartialDensityMatrix (int szSector, RealVector& groundState)
{
  Spin0_1_2_ChainWithTranslations TmpDestinationHilbertSpace(this->ChainLength, szSector,10000,10000);
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  unsigned long TmpBra,TmpKet,TmpState;

  for(int i=0;i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      for(int j=0;j < TmpDestinationHilbertSpace.HilbertSpaceDimension;j++)
	{
	  TmpState=0;
	  TmpBra = TmpDestinationHilbertSpace.ChainDescription[i];
	  TmpKet = TmpDestinationHilbertSpace.ChainDescription[j];
	  for (int p = 0;p <this->ChainLength;p++)
	    {
	      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpBra &0x3ul, TmpKet &0x3ul);
	      TmpBra>>=2;
	      TmpKet>>=2;
	    }

	  TmpDensityMatrix.SetMatrixElement(i,j,groundState[this->FindStateIndex (TmpState)]);
	}
    }
  return TmpDensityMatrix;
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix DoubledSpin0_1_2_ChainWithTranslations::EvaluatePartialDensityMatrix (int szSector, ComplexVector& groundState)
{
  Spin0_1_2_ChainWithTranslations TmpDestinationHilbertSpace(this->ChainLength, szSector,10000,10000);
  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  unsigned long TmpBra,TmpKet,TmpState;


  ComplexMatrix HRep(TmpDestinationHilbertSpace.HilbertSpaceDimension,TmpDestinationHilbertSpace.HilbertSpaceDimension,true);
  for(int i=0;i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      for(int j=0;j < TmpDestinationHilbertSpace.HilbertSpaceDimension;j++)
	{
	  TmpState=0;
	  TmpBra = TmpDestinationHilbertSpace.ChainDescription[i];
	  TmpKet = TmpDestinationHilbertSpace.ChainDescription[j];
	  for (int p = 0;p <this->ChainLength;p++)
	    {
	      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpBra &0x3ul, TmpKet &0x3ul);
	      TmpBra>>=2;
	      TmpKet>>=2;
	    }
	  int Index = this->FindStateIndex (TmpState);
	  if (Index < this->HilbertSpaceDimension ) 
	    HRep.SetMatrixElement(i,j,groundState[Index]);
	}
    }

  ComplexMatrix SquareRho( TmpDestinationHilbertSpace.HilbertSpaceDimension,  TmpDestinationHilbertSpace.HilbertSpaceDimension,true);
  SquareRho = HRep*HRep;
  
  Complex Trace = 0.0;
  Complex Tmp = 0.0;
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      SquareRho.GetMatrixElement(i,i,Tmp);
      Trace+=Tmp;
    }
  cout <<"Trace "<< Trace<<endl;
//  HRep/= sqrt(Norm(Trace));
  

  Complex Tmp1;
  Complex Tmp2;
  cout << "check hermiticity" << endl;

  double Error = 5e-8;
  
  for (int i = 0; i <  TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    for (int j = i; j <  TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
      {
	HRep.GetMatrixElement(i, j, Tmp1);
	HRep.GetMatrixElement(j, i, Tmp2);
	if (Norm(Tmp1 - Conj(Tmp2)) > Error )
	  {
	    cout << "error at " << i << " " << j << " : " << Tmp1 << " " << Tmp2 << " " << Norm(Tmp1 - Conj(Tmp2)) << " (should be lower than " << (Error ) << ")" << endl;
	  }
	HRep.GetMatrixElement(i,j,Tmp);
	TmpDensityMatrix.SetMatrixElement(i,j,Tmp);
      }  
  return TmpDensityMatrix;
}



void DoubledSpin0_1_2_ChainWithTranslations::ConvertToGeneralSpace(ComplexVector vSource,ComplexVector & vDestination)
{
  for(int i =0; i <this->HilbertSpaceDimension; i++)
    {
      vDestination[(long) this->ChainDescription[i]] = vSource[i];
    }
}

void DoubledSpin0_1_2_ChainWithTranslations::AddConvertFromGeneralSpace(ComplexVector vSource,ComplexVector & vDestination)
{
  for(int i =0; i <this->HilbertSpaceDimension; i++)
    {
      vDestination[i] = vSource[(long) this->ChainDescription[i]];
    }
}


void DoubledSpin0_1_2_ChainWithTranslations::ConvertToGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination)
{
  for(int i =0; i <this->HilbertSpaceDimension; i++)
    {
      vDestination[(long) this->ChainDescription[i]] = vSource[i]* sqrt ( ((double) this->NbrStateInOrbit[i]));
    }
}

void DoubledSpin0_1_2_ChainWithTranslations::AddConvertFromGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination)
{
  unsigned long TmpState;
  for(int i =0; i <this->HilbertSpaceDimension; i++)
    {
      TmpState = (unsigned long) this->ChainDescription[i];
      for(int p =0 ;p <	  this->NbrStateInOrbit[i];p++)
	{
	  vDestination[i] += this->TranslationPhase[p] * vSource[(int ) TmpState] / sqrt ( ((double) this->NbrStateInOrbit[i]));
 	  this->ApplySingleXTranslation(TmpState);
	}
    }
}

// evaluate all exponential factors
//   

void  DoubledSpin0_1_2_ChainWithTranslations::EvaluateExponentialFactors()
{
  this->TranslationPhase = new Complex[this->MaxXMomentum];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      this->TranslationPhase[i] = Phase(2.0 * M_PI * ((this->Momentum * ((double) i) / ((double) this->MaxXMomentum))));
    }
}

/*
// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix DoubledSpin0_1_2_ChainWithTranslations::EvaluatePartialDensityMatrix (int szSector, int momentumSector, ComplexVector& groundState)
{
  Spin0_1_2_ChainWithTranslations TmpDestinationHilbertSpace(this->ChainLength, momentumSector ,szSector,10000,10000);

  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  unsigned long TmpState,TmpBra,TmpKet,TmpCanonical;
  int NbrTranslation;
  for(int i=0;i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      for(int j=i;j < TmpDestinationHilbertSpace.HilbertSpaceDimension;j++)
	{
	  TmpState=0;
	  TmpBra = TmpDestinationHilbertSpace.ChainDescription[i];
	  TmpKet = TmpDestinationHilbertSpace.ChainDescription[j];
	  for (int l=0; l < TmpDestinationHilbertSpace.NbrStateInOrbit[i]; l++)
	    {
	      TmpDestinationHilbertSpace.ApplySingleXTranslation(TmpBra);
	      for (int t=0; t < TmpDestinationHilbertSpace.NbrStateInOrbit[j]; t++)
		{
		  TmpDestinationHilbertSpace.ApplySingleXTranslation(TmpKet);
		  TmpState=0;
		  for (int p = 0 ; p <this->ChainLength ; p++)
		    {
		      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpBra &0x3ul, TmpKet &0x3ul);
		      TmpBra>>=2;
		      TmpKet>>=2;
		    }
		  
		  this->FindCanonicalForm(TmpState,TmpCanonical,NbrTranslation);
		  
		  int Index = this->FindStateIndex (TmpState);
		  if (Index < this->HilbertSpaceDimension ) 
		    {
		      double TmpFactor=sqrt( (double) (TmpDestinationHilbertSpace.NbrStateInOrbit[i] * TmpDestinationHilbertSpace.NbrStateInOrbit[j])  /( double) this->NbrStateInOrbit[Index]); 
		      TmpDensityMatrix.SetMatrixElement(i,j,groundState[Index]*this->TranslationPhase[t]*this->TranslationPhase[l]*Conj(this->TranslationPhase[NbrTranslation])*TmpFactor );
		    }
		  
		}
	    }
	}
    }
  return TmpDensityMatrix;
}
*/


HermitianMatrix DoubledSpin0_1_2_ChainWithTranslations::EvaluatePartialDensityMatrix (int szSector, int momentumSector, ComplexVector& groundState)
{
  Spin0_1_2_ChainWithTranslations TmpDestinationHilbertSpace(this->ChainLength, momentumSector ,this->ChainLength / this->MaxXMomentum, szSector,10000,10000);

  int ComplementaryKSector = (this->Momentum - momentumSector) %  this->MaxXMomentum;
  if (ComplementaryKSector < 0)
    ComplementaryKSector +=  this->MaxXMomentum;
  
  //  int ComplementaryKSector = momentumSector;
  
  Spin0_1_2_ChainWithTranslations TmpHilbertSpace(this->ChainLength,ComplementaryKSector,this->ChainLength / this->MaxXMomentum,szSector,10000,10000);
  
  int MaxDimension = TmpDestinationHilbertSpace.HilbertSpaceDimension;
  
  if ( MaxDimension < TmpHilbertSpace.HilbertSpaceDimension )
    {
      MaxDimension = TmpHilbertSpace.HilbertSpaceDimension;
    }
  
  HermitianMatrix TmpDensityMatrix(MaxDimension, true);
  unsigned long TmpState,TmpBra,TmpKet,TmpCanonicalState,ReferenceBra,ReferenceKet; 
  int NbrTranslation;

  ComplexMatrix HRep( MaxDimension, MaxDimension,true);
  
  for(int i = 0 ; i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      for(int j = 0 ; j < TmpHilbertSpace.HilbertSpaceDimension;j++)
	{
	  ReferenceBra = TmpDestinationHilbertSpace.ChainDescription[i];
	  ReferenceKet = TmpHilbertSpace.ChainDescription[j];
	  
	  for(int t = 0; t < TmpDestinationHilbertSpace.NbrStateInOrbit[i]; t++)
	    {
	      for(int k = 0; k < TmpHilbertSpace.NbrStateInOrbit[j]; k++)
		{
		  TmpState = 0;
		  TmpBra = ReferenceBra;
		  TmpKet = ReferenceKet;
		  for (int p = 0 ; p < this->ChainLength ; p++)
		    {
		      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpBra &0x3ul, TmpKet &0x3ul);
		      TmpBra>>=2;
		      TmpKet>>=2;
		    }
		  
		  this->FindCanonicalForm(TmpState,TmpCanonicalState, NbrTranslation);
		  
		  int Index = this->FindStateIndex (TmpCanonicalState);
		  if (Index < this->HilbertSpaceDimension ) 
		    {
		      double TmpFactor = sqrt( (double) (TmpDestinationHilbertSpace.NbrStateInOrbit[i] * TmpHilbertSpace.NbrStateInOrbit[j] * ( double) this->NbrStateInOrbit[Index])); 
		      //	      double TmpFactor=1.0;
		      double Argument =  2.0 * M_PI * (NbrTranslation * this->Momentum  - t * momentumSector  - k * ComplementaryKSector) /  this->MaxXMomentum ;
		      HRep.AddToMatrixElement(i,j,groundState[Index]/TmpFactor*Phase(Argument));
		    }
		  TmpHilbertSpace.ApplySingleXTranslation(ReferenceKet);
		}
	      TmpDestinationHilbertSpace.ApplySingleXTranslation(ReferenceBra);
	    }
	}
    }

  ComplexMatrix SquareRho( MaxDimension, MaxDimension,true);
  SquareRho = HRep*HRep;

  Complex Trace = 0.0;
  Complex Tmp = 0.0;
  for (int i = 0; i < MaxDimension;i++)
    {
      SquareRho.GetMatrixElement(i,i,Tmp);
      Trace+=Tmp;
    }
  cout <<"Trace "<< Trace<<endl;
  HRep/= Phase(0.5*Arg(Trace));
  Complex Tmp1;
  Complex Tmp2;
  cout << "check hermiticity" << endl;
  
  double Error = 5e-8;
  
  for (int i = 0; i <   MaxDimension; ++i)
    for (int j = i; j <  MaxDimension; ++j)
      {
	HRep.GetMatrixElement(i, j, Tmp1);
	HRep.GetMatrixElement(j, i, Tmp2);
	if (Norm(Tmp1 - Conj(Tmp2)) > Error )
	  {
	    cout << "error at " << i << " " << j << " : " << Tmp1 << " " << Tmp2 << " " << Norm(Tmp1 - Conj(Tmp2)) << " (should be lower than " << (Error ) << ")" << endl;
	  }
	HRep.GetMatrixElement(i,j,Tmp);
	TmpDensityMatrix.SetMatrixElement(i,j,Tmp);
      }
  
  return TmpDensityMatrix;
}


void  DoubledSpin0_1_2_ChainWithTranslations::ApplyInversionSymmetry(ComplexVector & sourceVector,  ComplexVector & destinationVector)
{
//  cout <<"using void  DoubledSpin0_1_2_ChainWithTranslations::ApplyInversionSymmetry(ComplexVector & sourceVector,  ComplexVector & destinationVector)"<<endl;
  RealDiagonalMatrix ZMatrix(9,true);
  double  Tmp[3];
  Tmp[0] = -1.0;Tmp[1] = -1.0;Tmp[2] = 1.0;
  for(int i =0;i <9; i++)
    {
      ZMatrix.SetMatrixElement(i,i,Tmp[i%3] *Tmp[i/3] ); 
    }
//  cout <<ZMatrix<<endl;
  double TmpFactor;
  unsigned long OldState;
  int TmpIndice;
  unsigned long  TmpStateDescription;


  int tmp = this->FindNextInversionSymetricIndice(0);

  if((tmp ==this->HilbertSpaceDimension) )
    {
      tmp=0;
      while(Norm(sourceVector[tmp]) < 1e-8 )
	tmp++;
    }
  
  while(Norm(sourceVector[tmp]) < 1e-8 )
    {
      tmp=this->FindNextInversionSymetricIndice(tmp);
    }
  
  cout <<tmp<<" "<< Norm(sourceVector[tmp])<<endl;
//  sourceVector*=Conj(sourceVector[tmp])/Norm(sourceVector[tmp]);
  sourceVector/=Phase(Arg(sourceVector[tmp]));

//  this->NormalizeDensityMatrix(sourceVector);
  for (int i =0; i <this->HilbertSpaceDimension; i++)
    {
      cout <<" i = " <<i<< " : "<<this->ChainDescription[i]<<endl;;
      TmpFactor = 1.0;
      OldState = 0ul;
      TmpStateDescription =  this->ChainDescription[i]; 
      for (int p = 0 ; p <this->ChainLength ; p++)
	{
	  TmpIndice = (TmpStateDescription/this->PowerD[p]) % this->PowerD[1];
	  cout <<TmpIndice<<" ";
	  if  (p%2 == 0 ) 
	    TmpFactor*=ZMatrix(TmpIndice,TmpIndice);
	  OldState+= TmpIndice*this->PowerD[this->ChainLength-p-1];
	}

      cout <<endl<<" "<<OldState<< " "<<this->FindStateIndex(OldState)<<" "<<sourceVector[this->FindStateIndex(OldState)]<<endl;
      destinationVector[i] = Conj(TmpFactor*sourceVector[this->FindStateIndex(OldState)]) ;
    }
}


int DoubledSpin0_1_2_ChainWithTranslations::FindNextInversionSymetricIndice(int previousOne)
{
  unsigned long OldState;
  int TmpIndice;
  unsigned long  TmpStateDescription;
  for (int i =previousOne+1; i <this->HilbertSpaceDimension; i++)
    {
      OldState = 0ul;
      TmpStateDescription =  this->ChainDescription[i]; 
      for (int p = 0 ; p <this->ChainLength ; p++)
	{
	  TmpIndice = (TmpStateDescription/this->PowerD[p]) % this->PowerD[1];
	  OldState+= TmpIndice*this->PowerD[this->ChainLength-p-1];
	} 
      if (OldState ==  this->ChainDescription[i]) 
	return i;
    } 
  return this->HilbertSpaceDimension;
}


void  DoubledSpin0_1_2_ChainWithTranslations::ApplyInversionSymmetry(ComplexVector & sourceVector,  ComplexVector & destinationVector, bool translationFlag)
{
  if (translationFlag == false ) 
    {
      this->ApplyInversionSymmetry(sourceVector, destinationVector) ;
    }
  else
    {
      //  cout <<"using void  DoubledSpin0_1_2_ChainWithTranslations::ApplyInversionSymmetry(ComplexVector & sourceVector,  ComplexVector & destinationVector)"<<endl;
      RealDiagonalMatrix ZMatrix(9,true);
      double  Tmp[3];
      Tmp[0] = -1.0;Tmp[1] = -1.0;Tmp[2] = 1.0;
      for(int i =0;i <9; i++)
	{
	  ZMatrix.SetMatrixElement(i,i,Tmp[i%3] *Tmp[i/3] ); 
	}
      //  cout <<ZMatrix<<endl;
      double TmpFactor;
      unsigned long OldState;
      int TmpIndice;
      unsigned long  TmpStateDescription;
      unsigned long  TmpCanonicalState;
      int NbrTranslation;

      int tmp = this->FindNextInversionSymetricIndice(0);
      
      if((tmp ==this->HilbertSpaceDimension) )
	{
	  tmp=0;
	  while(Norm(sourceVector[tmp]) < 1e-8 )
	    tmp++;
	}
      
      while(Norm(sourceVector[tmp]) < 1e-8 )
	{
	  tmp=this->FindNextInversionSymetricIndice(tmp);
	}
      

      sourceVector*=Conj(sourceVector[tmp])/Norm(sourceVector[tmp]);
      for (int i =0; i <this->HilbertSpaceDimension; i++)
	{
	  cout <<" i = " <<i<< " : "<<this->ChainDescription[i]<<endl;;
	  TmpFactor = 1.0;
	  OldState = 0ul;
	  TmpStateDescription =  this->ChainDescription[i]; 
	  for (int p = 0 ; p <this->ChainLength ; p++)
	    {
	      TmpIndice = (TmpStateDescription/this->PowerD[p]) % this->PowerD[1];
	      cout <<TmpIndice<<" ";
	      if  (p%2 == 0 ) 
		TmpFactor*=ZMatrix(TmpIndice,TmpIndice);
	      OldState+= TmpIndice*this->PowerD[this->ChainLength-p-1];
	    }

	  this->FindCanonicalForm(OldState,TmpCanonicalState, NbrTranslation);	  
//	  cout <<endl<<" "<<OldState<< " "<<" "<<this->FindStateIndex(TmpCanonicalState)<<" "<< NbrTranslation <<" "<<sourceVector[this->FindStateIndex(TmpCanonicalState)]<<endl;


	  destinationVector[i] = Conj(TmpFactor*sourceVector[this->FindStateIndex(TmpCanonicalState)]*this->TranslationPhase[NbrTranslation]) ;
	}
    }
}


void  DoubledSpin0_1_2_ChainWithTranslations::NormalizeDensityMatrix(ComplexVector & sourceVector)
{
  unsigned long SourceState,TmpState;
  unsigned int TmpBra,TmpKet;
  Complex Trace = 0.0;
  for(int i = 0 ; i < this->HilbertSpaceDimension ; i++)
    {
      if ( Norm(sourceVector[i]) > 1e-8 ) 
	{
	  TmpState=0;
	  SourceState = this->ChainDescription[i];
	  for (int p = 0; p < this->ChainLength; p++)
	    {
	      this->GetBraAndKetIndicesFromCommonIndex(TmpBra,TmpKet, SourceState%9);
	      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpKet, TmpBra );
	      SourceState/=9;
	    } 
	  int Index = this->FindStateIndex (TmpState);
//	  cout <<endl<<"Index = " << Index<<endl;
	  if (Index < this->HilbertSpaceDimension ) 
	    {
	      sourceVector/=Phase(0.5*(Arg(sourceVector[i]) + Arg(sourceVector[Index])));
	      return;
	    }
	}
    }
}

void DoubledSpin0_1_2_ChainWithTranslations::ComputeDiffereenceSubLatticeNumberZero( unsigned long stateDescription , int & braNumber, int & ketNumber )
{
  unsigned long SourceState =  stateDescription;
  braNumber = 0;
  ketNumber = 0;
  unsigned  int TmpBra, TmpKet;
  for (int p = 0;p <this->ChainLength;p++)
    {
      this->GetBraAndKetIndicesFromCommonIndex(TmpBra,TmpKet, SourceState%9);
      if( p%2 == 0 )
	{
	  if (TmpBra == 2)
	    braNumber++;
	  if (TmpKet == 2) 
	    ketNumber++;

	}
      else
	{
	  if (TmpBra == 2)
	    braNumber--;
	  if (TmpKet == 2) 
	    ketNumber--;
	}
      SourceState/=9;
    }
}
