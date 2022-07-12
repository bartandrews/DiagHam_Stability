////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of Abstract spin chain with translations              //
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


#include "HilbertSpace/AbstractDoubledSpinChain.h"

#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <math.h>

using std::cout;
using std::endl;


#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif


// destructor
//

AbstractDoubledSpinChain::~AbstractDoubledSpinChain () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  delete[] this->ChainDescriptionBra;
	  delete[] this->ChainDescriptionKet;
	  delete[] this->UniqueStateDescriptionBra;
	  delete[] this->UniqueStateDescriptionSubArraySizeBra;
	  delete[] this->FirstIndexUniqueStateDescriptionBra;
	}
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

AbstractDoubledSpinChain & AbstractDoubledSpinChain::operator = (const AbstractDoubledSpinChain & chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->LargeHilbertSpaceDimension > 0l)
	{
 	  delete[] this->ChainDescriptionBra;
	  delete[] this->ChainDescriptionKet;
	  delete[] this->UniqueStateDescriptionBra;
	  delete[] this->UniqueStateDescriptionSubArraySizeBra;
	  delete[] this->FirstIndexUniqueStateDescriptionBra;
	}
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescriptionBra = chain.ChainDescriptionBra;
      this->ChainDescriptionKet = chain.ChainDescriptionKet;
      this->DiffSz = chain.DiffSz;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
   }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescriptionBra = 0;
      this->ChainDescriptionKet = 0;
      this->ChainLength = 0;
      this->DiffSz = 0;
      this->FixedSpinProjectionFlag = false;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// return Hilbert space dimension
//
// return value = Hilbert space dimension

int AbstractDoubledSpinChain::GetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}

// return value of twice spin projection of the Bra - the one of the ketfor a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int AbstractDoubledSpinChain::TotalSz (int index)
{
  if (this->FixedSpinProjectionFlag == true)
    return this->DiffSz;
  return this->GetTotalSz (this->ChainDescriptionBra[index],this->ChainDescriptionKet[index]);
}


// find state index
//
// state = state description
// return value = corresponding index

int AbstractDoubledSpinChain::FindStateIndex(unsigned long stateBra,unsigned long stateKet) 
{
  int PosMin = 0;
  int PosMax = this->NbrUniqueStateDescriptionBra - 1;
  int PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->UniqueStateDescriptionBra[PosMid];
  while (( (PosMax  - PosMin ) > 1 ) && (CurrentState != stateBra))
    {
       if (CurrentState > stateBra)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->UniqueStateDescriptionBra[PosMid];
    }
  
  if (CurrentState != stateBra)
    PosMid = PosMax;

  if (this->UniqueStateDescriptionBra[PosMid] != stateBra)
    return this->HilbertSpaceDimension;

  PosMin = this->FirstIndexUniqueStateDescriptionBra[PosMid];
  PosMax = PosMin + this->UniqueStateDescriptionSubArraySizeBra[PosMid] - 1;
  PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->ChainDescriptionKet[PosMid];
  while ( ( (PosMax  - PosMin ) > 1) && (CurrentState != stateKet))
    {
      if (CurrentState > stateKet)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->ChainDescriptionKet[PosMid];
    }
  if (this->ChainDescriptionKet[PosMax] == stateKet)
    return PosMax;
  if (this->ChainDescriptionKet[PosMid] == stateKet)
    return PosMid;
  return this->HilbertSpaceDimension;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void AbstractDoubledSpinChain::GenerateLookUpTable(unsigned long memory)
{  
  long TmpUniquePartition = 1l;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->ChainDescriptionBra[i - 1] == this->ChainDescriptionBra[i]))
	{
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	++TmpUniquePartition;
    }
  
  this->NbrUniqueStateDescriptionBra = TmpUniquePartition;
  this->UniqueStateDescriptionBra = new unsigned long [this->NbrUniqueStateDescriptionBra];
  this->UniqueStateDescriptionSubArraySizeBra = new int [this->NbrUniqueStateDescriptionBra];
  this->FirstIndexUniqueStateDescriptionBra = new int [this->NbrUniqueStateDescriptionBra];
  TmpUniquePartition = 0l;
  this->UniqueStateDescriptionBra[0l] = this->ChainDescriptionBra[0l];
  this->UniqueStateDescriptionSubArraySizeBra[0] = 1;
  this->FirstIndexUniqueStateDescriptionBra[0] = 0;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->ChainDescriptionBra[i - 1] == this->ChainDescriptionBra[i]))
	{
	  ++this->UniqueStateDescriptionSubArraySizeBra[TmpUniquePartition];
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	{
	  ++TmpUniquePartition;
	  this->UniqueStateDescriptionBra[TmpUniquePartition] = this->ChainDescriptionBra[i];
	  this->UniqueStateDescriptionSubArraySizeBra[TmpUniquePartition] = 1; 
	  this->FirstIndexUniqueStateDescriptionBra[TmpUniquePartition] = i;
	}
    }
}

/*
double AbstractDoubledSpinChainWithTranslations::TotalSzSz (int index)
{
  cout <<"Calling undefined function double AbstractDoubledSpinChainWithTranslations::TotalSzSz (int index) "<<endl;
  return 0.0;
}

double AbstractDoubledSpinChainWithTranslations::SziSzj (int i, int j, int state)
{
  cout <<"Calling undefined function double AbstractDoubledSpinChainWithTranslations::SziSzj (int i, int j, int state)"<<endl;
  return 0.0;
}

int AbstractDoubledSpinChainWithTranslations::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int AbstractDoubledSpinChainWithTranslations::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int AbstractDoubledSpinChainWithTranslations::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int AbstractDoubledSpinChainWithTranslations::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation) "<<endl;
  return 0;
}

int AbstractDoubledSpinChainWithTranslations::SpiSpi (int i, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int AbstractDoubledSpinChainWithTranslations::SpiSpi (int i, int state, double& coefficient, int& nbrTranslation) "<<endl;
  return 0;
}

int AbstractDoubledSpinChainWithTranslations::SmiSmi (int i, int state, double& coefficient, int& nbrTranslation) 
{
  cout <<"Calling undefined function int AbstractDoubledSpinChainWithTranslations::SmiSmi (int i, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int AbstractDoubledSpinChainWithTranslations::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int AbstractDoubledSpinChainWithTranslations::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}


int AbstractDoubledSpinChainWithTranslations::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int AbstractDoubledSpinChainWithTranslations::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int AbstractDoubledSpinChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int AbstractDoubledSpinChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int AbstractDoubledSpinChainWithTranslations::Spi (int i, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int AbstractDoubledSpinChainWithTranslations::Spi (int i, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}

int AbstractDoubledSpinChainWithTranslations::Smi (int i, int state, double& coefficient, int& nbrTranslation)
{
  cout <<"Calling undefined function int  AbstractDoubledSpinChainWithTranslations::Smi (int i, int state, double& coefficient, int& nbrTranslation)"<<endl;
  return 0;
}
*/

int AbstractDoubledSpinChain::FindStateIndexFromLinearizedIndex(unsigned long linearizedState)
{
  cout <<"Calling undefined function int  AbstractDoubledSpinChainWithTranslations::FindStateIndexFindStateIndexFromLinearizedIndex(unsigned long linearizedState)"<<endl;
  return 0;
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractDoubledSpinChain::Smi (int i, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractDoubledSpinChain::Smi" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i operator on a given state
//
// i = position of Sz operator
// state = index of the state to be applied on Sz_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

int AbstractDoubledSpinChain::Szi (int i, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractDoubledSpinChain::Szi" << endl;
  return this->HilbertSpaceDimension;
}
  
// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractDoubledSpinChain::SpiSpj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractDoubledSpinChain::SpiSpj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractDoubledSpinChain::SmiSmj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractDoubledSpinChain::SmiSmj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractDoubledSpinChain::SpiSzj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractSpinChainWithTranslations::SpiSzj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractDoubledSpinChain::SmiSzj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractDoubledSpinChain::SmiSzj" << endl;
  return this->HilbertSpaceDimension;
}
  
// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

int AbstractDoubledSpinChain::SmiSpj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractDoubledSpinChain::SmiSpj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractDoubledSpinChain::Spi (int i, int state, double& coefficient)
{
  cout << "warning, using dummy method AbstractDoubledSpinChain::Spi" << endl;
  return this->HilbertSpaceDimension;
}

double AbstractDoubledSpinChain::SziSzj (int i, int j, int state)
{
  cout <<"Calling undefined function double AbstractDoubledSpinChainWithTranslations::SziSzj (int i, int j, int state)"<<endl;
  return 0.0;
}


void AbstractDoubledSpinChain::GetChainDescriptionInCondensedForm(unsigned long * OldHilbertSpace)
{
  cout <<"Calling undefined function void AbstractDoubledSpinChain::GetChainDescriptionInCondensedForm(unsigned long * OldHilbertSpace)"<<endl;
}
