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


#include "HilbertSpace/VirtualSpacePEPSWithTranslations.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;


// default constructor
//

VirtualSpacePEPSWithTranslations::VirtualSpacePEPSWithTranslations () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->ChainDescription = 0;
  this->ChainLength = 0;
  this->Momentum = 0;
  this->CompatibilityWithMomentum = 0;
  this->RescalingFactors = 0;
  this->NbrStateInOrbit = 0;
  this->BondDimension=0;
  this->PowerD = 0;
  this->TranslationPhase = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}




// constructor for Hilbert space corresponding no constaint on Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

VirtualSpacePEPSWithTranslations::VirtualSpacePEPSWithTranslations (int chainLength, int bondDimension, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->BondDimension= bondDimension;
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  memorySize/=8;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->PowerD = new int [this->ChainLength];
  this->PowerD[0]=1;
  this->PowerD[1]=  this->BondDimension;
  for(int i = 2 ; i <this->ChainLength;i++)
    this->PowerD[i]=this->PowerD[i-1]*this->PowerD[1];
  
  this->LargeHilbertSpaceDimension =  this->PowerD[this->ChainLength-1]*this->PowerD[1];
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

VirtualSpacePEPSWithTranslations::VirtualSpacePEPSWithTranslations (int chainLength, int bondDimension, int momentum, int translationStep,  int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Momentum = momentum;
  this->BondDimension = bondDimension;

  this->PowerD = new int [this->ChainLength];
  this->PowerD[0]=1;
  this->PowerD[1]=  this->BondDimension;
  for(int i = 2 ; i <this->ChainLength;i++)
    this->PowerD[i]=this->PowerD[i-1]*this->PowerD[1] ;
  
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;
  
  this->MaxXMomentum = this->ChainLength/ translationStep;

  this->LargeHilbertSpaceDimension =  this->PowerD[this->ChainLength-1]*this->PowerD[1];

  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, 0l);
  
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
    this->GenerateLookUpTable();
}


// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

VirtualSpacePEPSWithTranslations::VirtualSpacePEPSWithTranslations (const VirtualSpacePEPSWithTranslations & chain)
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescription = chain.ChainDescription;
      this->Momentum = chain.Momentum;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->BondDimension = chain.BondDimension;
      this->PowerD = chain.PowerD;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->BondDimension = 0;
      this->Momentum = 0;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->PowerD = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

VirtualSpacePEPSWithTranslations::~VirtualSpacePEPSWithTranslations () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

VirtualSpacePEPSWithTranslations & VirtualSpacePEPSWithTranslations::operator = (const VirtualSpacePEPSWithTranslations & chain)
{
  this->Flag = chain.Flag;
  this->ChainLength = chain.ChainLength;
  this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
  this->LookUpTable = chain.LookUpTable;
  this->LookUpTableShift = chain.LookUpTableShift;
  this->ChainDescription = chain.ChainDescription;
  this->Momentum = chain.Momentum;
  this->BondDimension = chain.BondDimension;
  this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
  this->RescalingFactors = chain.RescalingFactors;
  this->NbrStateInOrbit = chain.NbrStateInOrbit;
  this->PowerD = chain.PowerD;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* VirtualSpacePEPSWithTranslations::Clone()
{
  return new VirtualSpacePEPSWithTranslations (*this);
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& VirtualSpacePEPSWithTranslations::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long tmpState;
  unsigned long StateDescription = this->ChainDescription[state];  
  Str << this->FindStateIndex(StateDescription) << " : "; 
  for (int j = this->ChainLength; j >0; j--)
    {
      tmpState =  StateDescription%this->PowerD[1];
      StateDescription/=this->PowerD[1];
      Str << tmpState<< " ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// length = length of the chain to be decided for bra spins
// diffSz = difference of spin projection between bra and ket chain
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
long VirtualSpacePEPSWithTranslations::GenerateStates(int length,  long pos)
{
  if (length == 0 )
    {
      for (int i =0; i < this->PowerD[1]; i++)
	{
	  this->ChainDescription[pos] = i;
	  pos++;
	}
      return pos;
    }
  
  if(length > 0)
    { 
      long TmpPos;
      for (int i = 0; i < this->PowerD[1]; i++)
	{
	  TmpPos = this->GenerateStates(length-1, pos); 
	  for (; pos < TmpPos; ++pos)
	    {
	      this->ChainDescription[pos] += i*this->PowerD[length];
	    }
	}
      return pos;
    }
  return pos;
}


double VirtualSpacePEPSWithTranslations::TotalSzSz (int index)
{
  cout <<"Calling undefined function double VirtualSpacePEPSWithTranslations::TotalSzSz (int index) "<<endl;
  return 0.0;
}

double VirtualSpacePEPSWithTranslations::SziSzj (int i, int j, int state)
{
  cout <<"Calling undefined function double VirtualSpacePEPSWithTranslations::SziSzj (int i, int j, int state)"<<endl;
  return 0.0;
}

int VirtualSpacePEPSWithTranslations::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int VirtualSpacePEPSWithTranslations::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int VirtualSpacePEPSWithTranslations::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int VirtualSpacePEPSWithTranslations::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation) "<<endl;
  return 0;
}

int VirtualSpacePEPSWithTranslations::SpiSpi (int i, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int VirtualSpacePEPSWithTranslations::SpiSpi (int i, int state, double& coefficient, int& nbrTranslation) "<<endl;
  return 0;
}

int VirtualSpacePEPSWithTranslations::SmiSmi (int i, int state, double& coefficient, int& nbrTranslation) 
{
  cout <<"Calling undefined function int VirtualSpacePEPSWithTranslations::SmiSmi (int i, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int VirtualSpacePEPSWithTranslations::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int VirtualSpacePEPSWithTranslations::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int VirtualSpacePEPSWithTranslations::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int VirtualSpacePEPSWithTranslations::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int VirtualSpacePEPSWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int VirtualSpacePEPSWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int VirtualSpacePEPSWithTranslations::Spi (int i, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int VirtualSpacePEPSWithTranslations::Spi (int i, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int VirtualSpacePEPSWithTranslations::Smi (int i, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int VirtualSpacePEPSWithTranslations::Smi (int i, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

void VirtualSpacePEPSWithTranslations::CreatePrecalculationTable()
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

void VirtualSpacePEPSWithTranslations::GenerateLookUpTable()
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
 

int VirtualSpacePEPSWithTranslations::FindStateIndex(unsigned long state)
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
