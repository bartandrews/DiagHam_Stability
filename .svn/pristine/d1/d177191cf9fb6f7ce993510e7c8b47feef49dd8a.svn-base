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

#include <bitset>
#include "HilbertSpace/Spin0_1_2_ChainWithTranslations.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;


// default constructor
//

Spin0_1_2_ChainWithTranslations::Spin0_1_2_ChainWithTranslations () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->ChainDescription = 0;
  this->ChainLength = 0;
  this->Momentum = 0;
  this->ComplementaryStateShift = 0;
  this->DiffSz = 0;
  this->FixedSpinProjectionFlag = false;
  this->CompatibilityWithMomentum = 0;
  this->RescalingFactors = 0;
  this->NbrStateInOrbit = 0;
  this->ShiftNegativeDiffSz = 0;
  this->ShiftLookUpTableNegativeSz=0;
  this->ShiftNegativeSz = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}




// constructor for Hilbert space corresponding no constaint on Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

Spin0_1_2_ChainWithTranslations::Spin0_1_2_ChainWithTranslations (int chainLength, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->ComplementaryStateShift = 2*(this->ChainLength - 1);
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
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

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

Spin0_1_2_ChainWithTranslations::Spin0_1_2_ChainWithTranslations (int chainLength, int diffSz, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->ComplementaryStateShift = 2*(this->ChainLength - 1);
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
  this->ShiftNegativeDiffSz = this->LargeHilbertSpaceDimension;

/*  if (this->DiffSz !=0 )
    this->LargeHilbertSpaceDimension += this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, -this->DiffSz);*/
  
  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->DiffSz, 0l);
/*  if (this->DiffSz != 0)
    TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1,-this->DiffSz, TmpHilbertSpaceDimension);*/

  SortArrayDownOrdering(this->ChainDescription ,TmpHilbertSpaceDimension);
  
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in Spin0_1_2_ChainWithTranslations!" << endl;
      exit(1);
    } 
  
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

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

Spin0_1_2_ChainWithTranslations::Spin0_1_2_ChainWithTranslations (int chainLength, int momentum, int translationStep, int diffSz, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->Momentum = momentum;

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
  this->ShiftNegativeDiffSz = this->LargeHilbertSpaceDimension;

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

  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      TmpState = this->ChainDescription[i];
      this->FindCanonicalForm(TmpState,TmpCanonicalState,NbrTranslation);
     
      if (TmpState  == TmpCanonicalState)
	{
	  //CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpCanonicalState)
	  if (this->TestMomentumConstraint(TmpCanonicalState) == true){
	      //	  if (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true){
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
	  if(this->GetTotalSz(this->ChainDescription[i]) >=0 )
	    this->ShiftNegativeDiffSz=this->LargeHilbertSpaceDimension;
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

Spin0_1_2_ChainWithTranslations::Spin0_1_2_ChainWithTranslations (const Spin0_1_2_ChainWithTranslations & chain)
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
      this->ShiftNegativeDiffSz = chain.ShiftNegativeDiffSz;
      this->ShiftNegativeSz =  chain.ShiftNegativeSz;
      this->ShiftLookUpTableNegativeSz = chain.ShiftLookUpTableNegativeSz;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->DiffSz = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->ShiftNegativeDiffSz =0;
      this->ShiftNegativeSz =0;
      this->ShiftLookUpTableNegativeSz = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

Spin0_1_2_ChainWithTranslations::~Spin0_1_2_ChainWithTranslations () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin0_1_2_ChainWithTranslations & Spin0_1_2_ChainWithTranslations::operator = (const Spin0_1_2_ChainWithTranslations & chain)
{
  this->Flag = chain.Flag;
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
  this->ShiftNegativeDiffSz = chain.ShiftNegativeDiffSz;
  this->ShiftNegativeSz = chain.ShiftNegativeSz;
  this->ShiftLookUpTableNegativeSz = chain.ShiftLookUpTableNegativeSz;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin0_1_2_ChainWithTranslations::Clone()
{
  return new Spin0_1_2_ChainWithTranslations (*this);
}

// return value of twice spin projection on (Oz) for a given state
//
// stateDescription = state to which the spin projection has to be evaluated
// return value = twice spin projection on (Oz)

inline int Spin0_1_2_ChainWithTranslations::GetTotalSz (unsigned long stateDescription)
{
  int TmpSz = 0;

  for (int i = 0; i < this->ChainLength; i++)
    {
      switch (stateDescription & 0x3ul)
	{
	case 0x1:
	  TmpSz += 1;
	  break;
	case 0x0:
	  TmpSz -= 1;
	  break;
	}
      stateDescription >>= 2;
    }
  return TmpSz;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin0_1_2_ChainWithTranslations::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long tmpState;
  unsigned long StateDescription = this->ChainDescription[state];  
  Str << this->FindStateIndex(StateDescription) << " : "; 
  for (int j = this->ChainLength; j >0; j--)
    {
      tmpState = ((StateDescription >> (( j -1) << 1)) & 0x3ul);
      
      if (tmpState == 0x0ul)
	Str << "d";
      else
	if (tmpState == 0x1ul)
	  Str << "u";
	else
	  Str << "0";
      Str << " ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// length = length of the chain to be decided for bra spins
// diffSz = difference of spin projection between bra and ket chain
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
long Spin0_1_2_ChainWithTranslations::GenerateStates(int length, int diffSz, long pos)
{
  if (length == 0)
    {
      if (diffSz == 0) 
	{
	  this->ChainDescription[pos] = 0x1ul<<1;
	  pos ++;
	  return pos;
	}
      if (diffSz == 1) 
	{
	  this->ChainDescription[pos] = 0x1ul;
	  pos ++;
	  return pos;
	}
      if (diffSz == -1) 
	{
	  this->ChainDescription[pos] = 0x0ul;
	  pos ++;
	  return pos;
	}
      return pos;
    }
  
  long TmpPos;
  unsigned long Mask;
  
  if(length > 0)
    { 
      TmpPos = this->GenerateStates(length-1, diffSz-1, pos); 
      Mask = (((0x1ul)) << ((length<<1)));
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescription[pos] |= Mask;
	}
      TmpPos = this->GenerateStates(length-1,diffSz, pos); 
      Mask = (((0x1ul << 1)) << ((length<<1)));
      for (; pos < TmpPos; ++pos)
	this->ChainDescription[pos] |= Mask;
      TmpPos = this->GenerateStates(length-1, diffSz+1, pos); 
      Mask = (((0x0ul)) << ((length<<1)));
      for (; pos < TmpPos; ++pos)
	this->ChainDescription[pos] |= Mask;
      return pos;
    }
}


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrNUp = number of particles with quantum number up
// nbrNDown = number of particles with quantum number down
// return value = Hilbert space dimension

long Spin0_1_2_ChainWithTranslations::ShiftedEvaluateHilbertSpaceDimension(int length, int diffSz)
{
  if (length == 0)
    {
      if ((diffSz == 0)||(diffSz == 1)||(diffSz == -1))
	{
	  return 1;
	}
      return 0;
    } 
  
  long Tmp=0;
  
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(length-1, diffSz+1); 
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(length-1, diffSz); 
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(length-1, diffSz-1);
  return Tmp;
}


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrNUp = number of particles with quantum number up
// nbrNDown = number of particles with quantum number down
// return value = Hilbert space dimension

long Spin0_1_2_ChainWithTranslations::ShiftedEvaluateHilbertSpaceDimension(int length)
{
  if (length == 0)
    {
      return 3;
    } 
  
  return 3*this->ShiftedEvaluateHilbertSpaceDimension(length-1); 
}


// generate all states corresponding to the constraints
// 
// length = length of the chain to be decided for bra spins
// diffSz = difference of spin projection between bra and ket chain
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
long Spin0_1_2_ChainWithTranslations::GenerateStates(int length,  long pos)
{
  if (length == 0)
    {
      this->ChainDescription[pos] = 0x1ul<<1;
      pos ++;
      this->ChainDescription[pos] = 0x1ul;
      pos ++;
      this->ChainDescription[pos] = 0x0ul;
      pos ++;
      return pos;
    }
  
  long TmpPos;
  unsigned long Mask;
  
  if(length > 0)
    { 
      TmpPos = this->GenerateStates(length-1, pos); 
      Mask = (((0x1ul)) << ((length<<1)));
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescription[pos] |= Mask;
	}
      TmpPos = this->GenerateStates(length-1, pos); 
      Mask = (((0x1ul << 1)) << ((length<<1)));
      for (; pos < TmpPos; ++pos)
	this->ChainDescription[pos] |= Mask;
      TmpPos = this->GenerateStates(length-1, pos); 
      Mask = (((0x0ul)) << ((length<<1)));
      for (; pos < TmpPos; ++pos)
	this->ChainDescription[pos] |= Mask;
      return pos;
    }
}


double Spin0_1_2_ChainWithTranslations::TotalSzSz (int index)
{
  cout <<"Calling undefined function double Spin0_1_2_ChainWithTranslations::TotalSzSz (int index) "<<endl;
  return 0.0;
}

double Spin0_1_2_ChainWithTranslations::SziSzj (int i, int j, int state)
{
  int TmpState = this->ChainDescription[state];
  int TmpI =  (( TmpState >> (i << 1) ) & 0x03ul );
  int TmpJ =  (( TmpState >> (j << 1) ) & 0x03ul );

  if ( ( TmpI == 2 ) ||  ( TmpJ == 2 ) )
    return 0.0;
  else
    {
      if ( TmpI == TmpJ )
	return 0.25;      
      else
	return -0.25;      
    }
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

int Spin0_1_2_ChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient)
{
  int TmpState = this->ChainDescription[state];
  int TmpI =  (( TmpState >> (i << 1) ) & 0x03ul );
  int TmpJ =  (( TmpState >> (j << 1) ) & 0x03ul );
  
  if ( ( TmpI == 1 ) && ( TmpJ == 0 ) )
    {
//      std::bitset<8> x (TmpState);
      coefficient = 1.0;
      unsigned long Mask = (0x1ul << (2*i)) | (0x1ul << (2*j));
      TmpState ^= Mask;
//      std::bitset<8> t (Mask);
//      std::bitset<8> y (TmpState);
//      cout <<x<<" " <<y<<" "<<t<<endl;
      return this->FindStateIndex(TmpState);
    }
  else
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension; 
    }
  
}




int Spin0_1_2_ChainWithTranslations::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int Spin0_1_2_ChainWithTranslations::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int Spin0_1_2_ChainWithTranslations::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int Spin0_1_2_ChainWithTranslations::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation) "<<endl;
  return 0;
}

int Spin0_1_2_ChainWithTranslations::SpiSpi (int i, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int Spin0_1_2_ChainWithTranslations::SpiSpi (int i, int state, double& coefficient, int& nbrTranslation) "<<endl;
  return 0;
}

int Spin0_1_2_ChainWithTranslations::SmiSmi (int i, int state, double& coefficient, int& nbrTranslation) 
{
  cout <<"Calling undefined function int Spin0_1_2_ChainWithTranslations::SmiSmi (int i, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int Spin0_1_2_ChainWithTranslations::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int Spin0_1_2_ChainWithTranslations::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int Spin0_1_2_ChainWithTranslations::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int Spin0_1_2_ChainWithTranslations::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int Spin0_1_2_ChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int Spin0_1_2_ChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int Spin0_1_2_ChainWithTranslations::Spi (int i, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int Spin0_1_2_ChainWithTranslations::Spi (int i, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int Spin0_1_2_ChainWithTranslations::Smi (int i, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int  Spin0_1_2_ChainWithTranslations::Smi (int i, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

void Spin0_1_2_ChainWithTranslations::CreatePrecalculationTable()
{
  int TmpPeriodicity = this->MaxXMomentum;
  this->CompatibilityWithMomentum = new bool [TmpPeriodicity + 1];
  for (int i = 0; i <= TmpPeriodicity; ++i)
    if (((i * this->Momentum) % TmpPeriodicity) == 0)
      this->CompatibilityWithMomentum[i] = true;
    else
      this->CompatibilityWithMomentum[i] = false;
 
  this->RescalingFactors = new double* [TmpPeriodicity + 1];
  for (int i = 1; i <= TmpPeriodicity; ++i)
    {
      this->RescalingFactors[i] = new double [TmpPeriodicity + 1];
      for (int j = 1; j <= TmpPeriodicity; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void Spin0_1_2_ChainWithTranslations::GenerateLookUpTable()
{  
  long LowPos;
  long MidPos;
  long HighPos;
  unsigned long Max2 = (this->ChainDescription[0]) >> this->LookUpTableShift;
  unsigned long Max3 = (this->ChainDescription[this->HilbertSpaceDimension-1]) >> this->LookUpTableShift;
  
  this->LookUpTable = new long [Max2 - Max3 +2];
  this->ShiftLookUpTable = Max3;
  
  for (long i = Max3; i <=Max2 ; i++)
    {
      LowPos = 0;
      HighPos = this->HilbertSpaceDimension - 1;
      while ((HighPos - LowPos) > 1)
	{
	  MidPos = (HighPos + LowPos) >> 1;
	  if (this->ChainDescription[MidPos] <= (i << this->LookUpTableShift))
	    HighPos = MidPos;
	  else
	    LowPos = MidPos;
	}      
      this->LookUpTable[i-this->ShiftLookUpTable] = HighPos;
    }
  
  this->LookUpTable[Max2 - Max3 +1] = 0;
}
 

int Spin0_1_2_ChainWithTranslations::FindStateIndex(unsigned long state)
{
  unsigned long MidPos = (state >> this->LookUpTableShift)-this->ShiftLookUpTable;  
  
  unsigned long LowPos = this->LookUpTable[MidPos+1];
  unsigned long HighPos = this->LookUpTable[MidPos];
  while ( ( HighPos - LowPos ) > 1)
    {
      MidPos = (HighPos + LowPos) >> 1;
      if (this->ChainDescription[MidPos] <= state)
	HighPos = MidPos;
      else
	LowPos = MidPos;
    }

  if (this->ChainDescription[LowPos] == state ) 
    return LowPos;
  if (this->ChainDescription[HighPos] == state ) 
    return HighPos;   
  return this->HilbertSpaceDimension;
}



//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// state = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint
bool Spin0_1_2_ChainWithTranslations::TestMomentumConstraint(unsigned long state) 
{
  return this->CompatibilityWithMomentum[this->FindNumberTranslation(state)];
}

// get the value of the spin (i.e. S) at a given site
// 
// site = site index
// return value = twice the spin

int Spin0_1_2_ChainWithTranslations::GetLocalSpin (int site, int state)
{
  int TmpI = (( this->ChainDescription[state] >> (site << 1) ) & 0x03ul );
  if (TmpI == 2)
    return 0.0;
  else
    return 1.0;

}
