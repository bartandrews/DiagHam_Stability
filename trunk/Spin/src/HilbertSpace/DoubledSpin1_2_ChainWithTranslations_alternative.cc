////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of doubled spin 1/2 chain with translations              //
//                                                                            //
//                        last modification : 01/05/2016                      //
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

#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/DoubledSpin1_2_ChainWithTranslations_alternative.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "GeneralTools/ArrayTools.h"
#include <iostream>
#include <math.h>

using std::cout;
using std::endl;


// default constructor
//

DoubledSpin1_2_ChainWithTranslations_alternative::DoubledSpin1_2_ChainWithTranslations_alternative () 
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
// chainLength = number of spins
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin1_2_ChainWithTranslations_alternative::DoubledSpin1_2_ChainWithTranslations_alternative (int chainLength, int momentum, int diffSz, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->Momentum = momentum;
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
    this->PowerD[i] = this->PowerD[i-1] * 4;
  
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, this->DiffSz);
  this->ShiftNegativeDiffSz = this->LargeHilbertSpaceDimension;
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
}


// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

DoubledSpin1_2_ChainWithTranslations_alternative::DoubledSpin1_2_ChainWithTranslations_alternative (const DoubledSpin1_2_ChainWithTranslations_alternative & chain)
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
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

DoubledSpin1_2_ChainWithTranslations_alternative::~DoubledSpin1_2_ChainWithTranslations_alternative () 
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

DoubledSpin1_2_ChainWithTranslations_alternative & DoubledSpin1_2_ChainWithTranslations_alternative::operator = (const DoubledSpin1_2_ChainWithTranslations_alternative & chain)
{
  AbstractDoubledSpinChainWithTranslations::operator =(chain);
  this->ShiftNegativeDiffSz = chain.ShiftNegativeDiffSz;
  this->BraShiftNegativeSz = chain.BraShiftNegativeSz;
  this->PowerD = chain.PowerD;
  this->ChainDescription = chain.ChainDescription;
  this->TranslationPhase = chain.TranslationPhase;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* DoubledSpin1_2_ChainWithTranslations_alternative::Clone()
{
  return new DoubledSpin1_2_ChainWithTranslations_alternative (*this);
}

// return value of twice spin projection on (Oz) for a given state
//
// stateDescription = state to which the spin projection has to be evaluated
// return value = twice spin projection on (Oz)

inline int DoubledSpin1_2_ChainWithTranslations_alternative::GetTotalSz (unsigned long stateDescriptionBra, unsigned long stateDescriptionKet)
{
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      switch (stateDescriptionBra & 0x1ul)
	{
	case 0x1:
	  TmpSz += 1;
	  break;
	case 0x0:
	  TmpSz -= 1;
	  break;
	}
      stateDescriptionBra >>= 1;
      switch (stateDescriptionKet & 0x1ul)
	{
	case 0x1ul:
	  TmpSz -= 1;
	  break;
	case 0x0:
	  TmpSz += 1;
	  break;
	}
      stateDescriptionKet >>= 1;
    }
  return TmpSz;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& DoubledSpin1_2_ChainWithTranslations_alternative::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned int tmpBra,tmpKet;
  unsigned long StateDescription = this->ChainDescription[state];  

  Str << this->FindStateIndex(StateDescription) << " : " << StateDescription<< " : ";
  for (int j = this->ChainLength; j > 0; j--)
    {
      tmpBra = ((StateDescription >> (2 * j - 1) ) & 0x1ul);
      tmpKet =  ((StateDescription >> (2 * j - 2)) & 0x1ul);
      
      Str << "(";
      if (tmpBra == 0)
	Str << "d ";
      else
	Str << "u ";
      Str << ",";
      if (tmpKet == 0)
	Str << "d ";
      else
	Str << "u ";
      Str << ") ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// lengthBra = length of the chain to be decided for bra spins
// lengthBra = length of the chain to be decided for ket spins
// diffSz = difference of spin projection between bra and ket chain
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
long DoubledSpin1_2_ChainWithTranslations_alternative::GenerateStates(int lengthBra, int lengthKet, int diffSz, long pos)
{
  if ((lengthKet == 0) && (lengthBra == 0))
    {
      if (diffSz == 0) 
	{
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(1,1);
	  pos++;
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(0,0);
	  pos++;
	  return pos;
	}
      if (diffSz == 1) 
	{
	  this->ChainDescription[pos] = GetCommonIndexFromBraAndKetIndices(1,0);
	  pos++;
	  return pos;
	}
      if (diffSz == -1) 
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
      TmpPos = this->GenerateStates(lengthBra,lengthKet-1, diffSz, pos); 
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescription[pos] += GetCommonIndexFromBraAndKetIndices(0,1)*this->PowerD[lengthKet];
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
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz, pos); 
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescription[pos] += GetCommonIndexFromBraAndKetIndices(1,0)*this->PowerD[lengthBra];
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

long DoubledSpin1_2_ChainWithTranslations_alternative::ShiftedEvaluateHilbertSpaceDimension(int lengthBra, int lengthKet, int diffSz)
{
  if ((lengthBra < 0) || (lengthKet < 0))
    return 0;
  
  if ((lengthBra == 0) && (lengthKet == 0))
    {
      if (diffSz == 0) 
	{
	  return 2;
	}
      if (diffSz == 1) 
	{
	  return 1;
	}
      if (diffSz == -1) 
	{
	  return 1;
	}
      return 0;
    }  
  long Tmp=0;
  
  if (lengthBra == 0)
    {
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra,lengthKet-1, diffSz); 
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra,lengthKet-1, diffSz-1); 
      return Tmp;
    }
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra-1,lengthKet, diffSz); 
  Tmp +=  this->ShiftedEvaluateHilbertSpaceDimension(lengthBra-1,lengthKet, diffSz+1);
  return Tmp;
}



// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void DoubledSpin1_2_ChainWithTranslations_alternative::GenerateLookUpTable(unsigned long memory)
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
 

int DoubledSpin1_2_ChainWithTranslations_alternative::FindStateIndex(unsigned long stateDescription)
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

RealSymmetricMatrix DoubledSpin1_2_ChainWithTranslations_alternative::EvaluatePartialDensityMatrix (int szSector, RealVector& groundState)
{
  Spin1_2Chain TmpDestinationHilbertSpace(this->ChainLength, szSector, 10000);
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  unsigned long TmpBra,TmpKet,TmpState;

  for (int i=0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; i++)
    {
      for (int j=0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; j++)
	{
	  TmpState=0;
	  TmpBra = TmpDestinationHilbertSpace.StateDescription[i];
	  TmpKet = TmpDestinationHilbertSpace.StateDescription[j];
	  for (int p = 0;p <this->ChainLength;p++)
	    {
	      TmpState += this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpBra &0x1ul, TmpKet &0x1ul);
	      TmpBra >>= 1;
	      TmpKet >>= 1;
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

HermitianMatrix DoubledSpin1_2_ChainWithTranslations_alternative::EvaluatePartialDensityMatrix (int szSector, ComplexVector& groundState)
{
  Spin1_2Chain TmpDestinationHilbertSpace(this->ChainLength, szSector, 10000);
  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  ComplexMatrix HRep(TmpDestinationHilbertSpace.HilbertSpaceDimension,TmpDestinationHilbertSpace.HilbertSpaceDimension,true);

  unsigned long TmpBra,TmpKet,TmpState;

  for(int i=0;i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      for(int j=0;j < TmpDestinationHilbertSpace.HilbertSpaceDimension;j++)
	{
	  TmpState=0;
	  TmpBra = TmpDestinationHilbertSpace.StateDescription[i];
	  TmpKet = TmpDestinationHilbertSpace.StateDescription[j];
	  for (int p = 0;p <this->ChainLength;p++)
	    {
	      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpBra &0x1ul, TmpKet &0x1ul);
	      TmpBra >>= 1;
	      TmpKet >>= 1;
	    }
//	  cout<<i<< " "<<j<<" "<<TmpState <<endl;
	  int Index = this->FindStateIndex (TmpState);
	  if (Index < this->HilbertSpaceDimension ) 
	    HRep.SetMatrixElement(i,j,groundState[Index]);
	}
    }
  
  Complex Tmp1;
  Complex Tmp2;
  cout << "check hermiticity" << endl;
  double AverageNorm = 0.1;
  double Error = MACHINE_PRECISION;
  
  
  for (int i = 0; i <  TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    for (int j = i; j <  TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
      {
	HRep.GetMatrixElement(i, j, Tmp1);
	HRep.GetMatrixElement(j, i, Tmp2);
	if (Norm(Tmp1 - Conj(Tmp2)) > Error )
	  {
	    cout << "error at " << i << " " << j << " : " << Tmp1 << " " << Tmp2 << " " << Norm(Tmp1 - Conj(Tmp2)) << " (should be lower than " << (Error * AverageNorm) << ")" << endl;
	  }
      }  
  
  return TmpDensityMatrix;
}


void DoubledSpin1_2_ChainWithTranslations_alternative::ConvertToGeneralSpace(ComplexVector vSource,ComplexVector & vDestination)
{
  for(int i = 0; i <this->HilbertSpaceDimension; i++)
    {
      vDestination[(long) this->ChainDescription[i]] = vSource[i];
    }
}

void DoubledSpin1_2_ChainWithTranslations_alternative::AddConvertFromGeneralSpace(ComplexVector vSource,ComplexVector & vDestination)
{
  for(int i = 0; i <this->HilbertSpaceDimension; i++)
    {
      vDestination[i] = vSource[(long) this->ChainDescription[i]];
    }
}


void DoubledSpin1_2_ChainWithTranslations_alternative::ConvertToGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination)
{
  for(int i = 0; i <this->HilbertSpaceDimension; i++)
    {
      vDestination[(long) this->ChainDescription[i]] = vSource[i]* sqrt ( ((double) this->NbrStateInOrbit[i]));
    }
}

void DoubledSpin1_2_ChainWithTranslations_alternative::AddConvertFromGeneralSpaceWithMomentum(ComplexVector vSource,ComplexVector & vDestination)
{
  int TmpState;
  for(int i =0; i <this->HilbertSpaceDimension; i++)
    {
      TmpState = (long) this->ChainDescription[i];
      for(int p =0 ;p <	  this->NbrStateInOrbit[i];p++)
	{
	  vDestination[i] += this->TranslationPhase[p] * vSource[TmpState] / sqrt ( ((double) this->NbrStateInOrbit[i]));
	  TmpState = ((long) TmpState/4) + ((long) TmpState%4)*this->PowerD[this->ChainLength-1];
	}
    }
}

// evaluate all exponential factors
//   

void  DoubledSpin1_2_ChainWithTranslations_alternative::EvaluateExponentialFactors()
{
  this->TranslationPhase = new Complex[this->ChainLength];
  for (int i = 0; i < this->ChainLength; ++i)
    { 
      this->TranslationPhase[i] = Phase(2.0 * M_PI * ((this->Momentum * ((double) i) / ((double) this->ChainLength))));
    }
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix DoubledSpin1_2_ChainWithTranslations_alternative::EvaluatePartialDensityMatrix (int szSector, int momentumSector, ComplexVector& groundState)
{
  Spin1_2ChainWithTranslations TmpDestinationHilbertSpace(this->ChainLength, momentumSector ,szSector,10000,10000);

  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  unsigned long TmpState,TmpBra,TmpKet,TmpCanonical;
  int NbrTranslation;
  for(int i=0;i < TmpDestinationHilbertSpace.HilbertSpaceDimension;i++)
    {
      for(int j=i;j < TmpDestinationHilbertSpace.HilbertSpaceDimension;j++)
	{
	  TmpState=0;
	  TmpBra = TmpDestinationHilbertSpace.StateDescription[i];
	  TmpKet = TmpDestinationHilbertSpace.StateDescription[j];
	  for (int l=0; l < TmpDestinationHilbertSpace.NbrStateInOrbit[i]; l++)
	    {
	      TmpDestinationHilbertSpace.ApplySingleXTranslation(TmpBra);
	      for (int t=0; t < TmpDestinationHilbertSpace.NbrStateInOrbit[j]; t++)
		{
		  TmpDestinationHilbertSpace.ApplySingleXTranslation(TmpKet);
		  TmpState=0;
		  for (int p = 0 ; p <this->ChainLength ; p++)
		    {
		      TmpState+= this->PowerD[p] * this->GetCommonIndexFromBraAndKetIndices(TmpBra &0x1ul, TmpKet &0x1ul);
		      TmpBra >>= 1;
		      TmpKet >>= 1;
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
