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


#include "HilbertSpace/DoubledSpin1_2_Chain.h"
#include "Matrix/HermitianMatrix.h"
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


// default constructor
//

DoubledSpin1_2_Chain::DoubledSpin1_2_Chain () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->ChainDescriptionBra = 0;
  this->ChainDescriptionKet = 0;
  this->ChainLength = 0;
  this->ComplementaryStateShift = 0;
  this->DiffSz = 0;
  this->FixedSpinProjectionFlag = false;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}



// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin1_2_Chain::DoubledSpin1_2_Chain (int chainLength, int diffSz, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->ComplementaryStateShift = (this->ChainLength - 1) << 1;
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;
  
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, this->DiffSz);
  
  this->ChainDescriptionBra = new unsigned long [this->LargeHilbertSpaceDimension];
  this->ChainDescriptionKet = new unsigned long [this->LargeHilbertSpaceDimension];
  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, this->DiffSz, 0l);
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2Spin!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
    }
}
  
// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

DoubledSpin1_2_Chain::DoubledSpin1_2_Chain (const DoubledSpin1_2_Chain & chain)
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescriptionBra = chain.ChainDescriptionBra;
      this->ChainDescriptionKet = chain.ChainDescriptionKet;
      this->DiffSz = chain.DiffSz;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->UniqueStateDescriptionBra = chain.UniqueStateDescriptionBra;
      this->UniqueStateDescriptionSubArraySizeBra = chain.UniqueStateDescriptionSubArraySizeBra;
      this->NbrUniqueStateDescriptionBra = chain.NbrUniqueStateDescriptionBra;
      this->FirstIndexUniqueStateDescriptionBra = chain.FirstIndexUniqueStateDescriptionBra;

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
      this->DiffSz = 0;
      this->FixedSpinProjectionFlag = false;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

DoubledSpin1_2_Chain::~DoubledSpin1_2_Chain () 
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

DoubledSpin1_2_Chain & DoubledSpin1_2_Chain::operator = (const DoubledSpin1_2_Chain & chain)
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
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
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
      this->ComplementaryStateShift = 0;
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

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* DoubledSpin1_2_Chain::Clone()
{
  return new DoubledSpin1_2_Chain (*this);
}

// return value of twice spin projection of the Bra - the one of the ket for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int DoubledSpin1_2_Chain::TotalSz (int index)
{
  if (this->FixedSpinProjectionFlag == true)
    return this->DiffSz;
  return this->GetTotalSz (this->ChainDescriptionBra[index],this->ChainDescriptionKet[index]);
}

// return value of twice spin projection on (Oz) for a given state
//
// stateDescription = state to which the spin projection has to be evaluated
// return value = twice spin projection on (Oz)

inline int DoubledSpin1_2_Chain::GetTotalSz (unsigned long stateDescriptionBra,unsigned long stateDescriptionKet)
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

// find state index
//
// state = state description
// return value = corresponding index

int DoubledSpin1_2_Chain::FindStateIndex(unsigned long stateBra,unsigned long stateKet) 
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

void DoubledSpin1_2_Chain::GenerateLookUpTable(unsigned long memory)
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


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& DoubledSpin1_2_Chain::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long tmpBra,tmpKet;
  unsigned long StateDescriptionBra = this->ChainDescriptionBra[state];  
  unsigned long StateDescriptionKet = this->ChainDescriptionKet[state];  
  Str << this->FindStateIndex(StateDescriptionBra,StateDescriptionKet) << " : ";
  for (int j = this->ChainLength - 1; j >=0; j--)
    {
      tmpBra = ((StateDescriptionBra >> j ) & 0x1ul);
      tmpKet =  ((StateDescriptionKet >> j) & 0x1ul);

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

long DoubledSpin1_2_Chain::GenerateStates(int lengthBra, int lengthKet, int diffSz, long pos)
{
  if (lengthKet == 0)
    {
      if ( diffSz == -1 )
	{
	  this->ChainDescriptionKet[pos] = 0x1ul;
	  pos ++;
	  return pos;
	}
      if ( diffSz == 1)
	{
	  this->ChainDescriptionKet[pos] = 0x0ul;
	  pos ++;
	  return pos;
	}
      return pos;
    }

 if(lengthBra == 0)
    { 
      long TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz-1, pos); 
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescriptionBra[pos] = 0x1ul;
	}
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz+1, pos); 
      for (; pos < TmpPos; ++pos)
	this->ChainDescriptionBra[pos] = 0x0ul;
      return pos;
    }

  long TmpPos;
  unsigned long MaskBra;
  unsigned long MaskKet;
  
  if(lengthBra > 0)
    { 
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz-1, pos); 
      MaskBra = 0x1ul  << lengthBra;
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescriptionBra[pos] |= MaskBra;
	}      
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz+1, pos); 
      return TmpPos;
    }
  
  if (lengthKet > 0)
    {
      TmpPos = this->GenerateStates(lengthBra,lengthKet-1, diffSz+1, pos); 
      MaskKet = 0x1ul << lengthKet;
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescriptionKet[pos] |= MaskKet;
	}
      TmpPos = this->GenerateStates(lengthBra,lengthKet-1, diffSz-1, pos); 
      return  TmpPos;
    }
}


double DoubledSpin1_2_Chain::TotalSzSz (int index)
{
  cout <<"Calling undefined function double DoubledSpin1_2_chain::TotalSzSz (int index)"<<endl;
  return 0.0;
}

double DoubledSpin1_2_Chain::SziSzj (int i, int j, int state)
{
  cout <<"Calling undefined function double DoubledSpin1_2_chain::SziSzj (int i, int j, int state)"<<endl;
  return 0.0;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrNUp = number of particles with quantum number up
// nbrNDown = number of particles with quantum number down
// return value = Hilbert space dimension

long DoubledSpin1_2_Chain::ShiftedEvaluateHilbertSpaceDimension(int lengthBra, int lengthKet, int diffSz)
{
  if ((lengthBra < 0) || (lengthKet < 0))
    return 0;
  
  if ((lengthBra == 0) && (lengthKet == 0))
    {
      if (diffSz == 0) 
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
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra,lengthKet-1, diffSz+1); 
      return Tmp;
    }

  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra-1,lengthKet, diffSz+1); 
  Tmp +=  this->ShiftedEvaluateHilbertSpaceDimension(lengthBra-1,lengthKet, diffSz-1);
  return Tmp;
}




// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int DoubledSpin1_2_Chain::Smi (int i, int state, double& coefficient)
{
  cout << "warning, using dummy method DoubledSpin1_2_chain::Smi" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i operator on a given state
//
// i = position of Sz operator
// state = index of the state to be applied on Sz_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

int DoubledSpin1_2_Chain::Szi (int i, int state, double& coefficient)
{
  cout << "warning, using dummy method DoubledSpin1_2_chain::Szi" << endl;
  return this->HilbertSpaceDimension;
}
  
// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int DoubledSpin1_2_Chain::SpiSpj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method DoubledSpin1_2_chain::SpiSpj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int DoubledSpin1_2_Chain::SmiSmj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method DoubledSpin1_2_chain::SmiSmj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int DoubledSpin1_2_Chain::SpiSzj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method DoubledSpin1_2_chain::SpiSzj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int DoubledSpin1_2_Chain::SmiSzj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method DoubledSpin1_2_chain::SmiSzj" << endl;
  return this->HilbertSpaceDimension;
}
  
// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

int DoubledSpin1_2_Chain::SmiSpj (int i, int j, int state, double& coefficient)
{
  cout << "warning, using dummy method DoubledSpin1_2_chain::SmiSpj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int DoubledSpin1_2_Chain::Spi (int i, int state, double& coefficient)
{
  cout << "warning, using dummy method DoubledSpin1_2_chain::Spi" << endl;
  return this->HilbertSpaceDimension;
}
