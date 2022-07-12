////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of spin 2 chain with translations                  //
//                                                                            //
//                        last modification : 09/11/2016                      //
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


#include "HilbertSpace/Spin2ChainWithTranslations.h"
#include "HilbertSpace/Spin2Chain.h"
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
using std::dec;
using std::hex;


#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif
#ifndef M_SQRT6
#define M_SQRT6	2.4494897427831781
#endif


// default constructor
//

Spin2ChainWithTranslations::Spin2ChainWithTranslations () 
{
  this->Spin3ProjectorNbrCoefficients = 0;
  this->Spin4ProjectorNbrCoefficients = 0;
  this->Spin3ProjectorCoefficients = 0;
  this->Spin4ProjectorCoefficients = 0;
  this->Spin3ProjectorStates = 0;
  this->Spin4ProjectorStates = 0;
}

// constructor for Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// memory = amount of memory granted for precalculations

Spin2ChainWithTranslations::Spin2ChainWithTranslations (int chainLength, int momentum, unsigned long memory) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->Momentum = momentum;

  this->MaxXMomentum = this->ChainLength;
  this->StateXShift = 3 * (this->ChainLength / this->MaxXMomentum);
  this->ComplementaryStateXShift = (3 * this-> ChainLength) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->Spin3ProjectorNbrCoefficients = 0;
  this->Spin4ProjectorNbrCoefficients = 0;
  this->Spin3ProjectorCoefficients = 0;
  this->Spin4ProjectorCoefficients = 0;
  this->Spin3ProjectorStates = 0;
  this->Spin4ProjectorStates = 0;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->ChainLength);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->BuildProjectors();
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory +=  this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->ChainLength * sizeof(int);
      UsedMemory += this->ChainLength * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}

// constructor for Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memory = amount of memory granted for precalculations

Spin2ChainWithTranslations::Spin2ChainWithTranslations (int chainLength, int momentum, int sz, unsigned long memory) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz;
  this->FixedSpinProjectionFlag = true;
  this->Momentum = momentum;

  this->MaxXMomentum = this->ChainLength;
  this->StateXShift = 3 * (this->ChainLength / this->MaxXMomentum);
  this->ComplementaryStateXShift = (3 * this-> ChainLength) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->Spin3ProjectorNbrCoefficients = 0;
  this->Spin4ProjectorNbrCoefficients = 0;
  this->Spin3ProjectorCoefficients = 0;
  this->Spin4ProjectorCoefficients = 0;
  this->Spin3ProjectorStates = 0;
  this->Spin4ProjectorStates = 0;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->ChainLength);
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1, 0);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->BuildProjectors();
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory +=  this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
      UsedMemory = this->ChainLength * sizeof(int);
      UsedMemory += this->ChainLength * this->LookUpTableMemorySize * sizeof(int);
      cout << "memory requested for lookup table = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin2ChainWithTranslations::Spin2ChainWithTranslations (const Spin2ChainWithTranslations& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->MaxXMomentum = chain.MaxXMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;

      this->LookUpTable = chain.LookUpTable;
      this->MaximumLookUpShift = chain.MaximumLookUpShift;
      this->LookUpTableMemorySize = chain.LookUpTableMemorySize;
      this->LookUpTableShift = chain.LookUpTableShift;

      this->Spin3ProjectorNbrCoefficients = chain.Spin3ProjectorNbrCoefficients;
      this->Spin4ProjectorNbrCoefficients = chain.Spin4ProjectorNbrCoefficients;
      this->Spin3ProjectorCoefficients = chain.Spin3ProjectorCoefficients;
      this->Spin4ProjectorCoefficients = chain.Spin4ProjectorCoefficients;
      this->Spin3ProjectorStates = chain.Spin3ProjectorStates;
      this->Spin4ProjectorStates = chain.Spin4ProjectorStates;
    }
  else
    {
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0l;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;

      this->MaxXMomentum = 0;
      this->StateXShift = 0;
      this->ComplementaryStateXShift = 0;
      this->XMomentumMask = 0x0ul;

      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;

      this->Spin3ProjectorNbrCoefficients = 0;
      this->Spin4ProjectorNbrCoefficients = 0;
      this->Spin3ProjectorCoefficients = 0;
      this->Spin4ProjectorCoefficients = 0;
      this->Spin3ProjectorStates = 0;
      this->Spin4ProjectorStates = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

Spin2ChainWithTranslations::~Spin2ChainWithTranslations () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Spin3ProjectorNbrCoefficients != 0)
	{
	  for (int i = 0; i < 64; ++i)
	    {
	      if (this->Spin3ProjectorNbrCoefficients[i] > 0)
		{
		  delete[] this->Spin3ProjectorCoefficients[i];
		  delete[] this->Spin3ProjectorStates[i];
		}
	    }
	  delete[] this->Spin3ProjectorNbrCoefficients;
	  delete[] this->Spin3ProjectorCoefficients;
	  delete[] this->Spin3ProjectorStates;
	}
      if (this->Spin4ProjectorNbrCoefficients != 0)
	{
	  for (int i = 0; i < 64; ++i)
	    {
	      if (this->Spin4ProjectorNbrCoefficients[i] > 0)
		{
		  delete[] this->Spin4ProjectorCoefficients[i];
		  delete[] this->Spin4ProjectorStates[i];
		}
	    }
	  delete[] this->Spin4ProjectorNbrCoefficients;
	  delete[] this->Spin4ProjectorCoefficients;
	  delete[] this->Spin4ProjectorStates;
	}
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin2ChainWithTranslations& Spin2ChainWithTranslations::operator = (const Spin2ChainWithTranslations& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  delete[] this->StateDescription;
	  delete[] this->LookUpTable;
	  for (int i = 1; i <= this->ChainLength; ++i)
	    {
	      delete[] this->RescalingFactors[i];
	    } 
	  delete[] this->RescalingFactors;
	  delete[] this->NbrStateInOrbit;
	}
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->MaxXMomentum = chain.MaxXMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;

      this->LookUpTable = chain.LookUpTable;
      this->MaximumLookUpShift = chain.MaximumLookUpShift;
      this->LookUpTableMemorySize = chain.LookUpTableMemorySize;
      this->LookUpTableShift = chain.LookUpTableShift;
   }
  else
    {
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0l;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;

      this->MaxXMomentum = 0;
      this->StateXShift = 0;
      this->ComplementaryStateXShift = 0;
      this->XMomentumMask = 0x0ul;

      this->LookUpTable = 0;
      this->MaximumLookUpShift = 0;
      this->LookUpTableMemorySize = 0;
      this->LookUpTableShift = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin2ChainWithTranslations::Clone()
{
  return new Spin2ChainWithTranslations (*this);
}

// get the normalization factor in front of each basis state (i.e. 1/sqrt(orbit size))
//
// return value = pointer to normalization factors

double* Spin2ChainWithTranslations::GetBasisNormalization()
{
  double* TmpNorm = new double[this->HilbertSpaceDimension];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      TmpNorm[i] = 1.0 / sqrt((double) this->NbrStateInOrbit[i]);
    }
  return TmpNorm;
}
  

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin2ChainWithTranslations::TotalSz (int index)
{
  if (this->FixedSpinProjectionFlag == true)
    return this->Sz;
  unsigned long State = this->StateDescription[index];
  int TmpSz = 0;
  unsigned long TmpState;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpState = State & 0x7ul;
      switch (TmpState)
	{
	case 0x4ul:
	  TmpSz += 4;
	  break;
	case 0x3ul:
	  TmpSz += 2;
	  break;
	case 0x1ul:
	  TmpSz -= 2;
	  break;
	case 0x0ul:
	  TmpSz -= 4;
	  break;
	}
      State >>= 3;
    }
  return TmpSz;
}

// return value of the value of the sum of the square of spin projection on (Oz) 
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

double Spin2ChainWithTranslations::TotalSzSz (int index)
{  
  unsigned long State = this->StateDescription[index];
  double TmpSzSz = 0.0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      if ((State & 0x3ul) != 0x2ul)
	TmpSzSz += 1.0;
      State >>= 3;
    }
  return TmpSzSz;
}

// return value of twice spin projection on (Oz) for a given state
//
// stateDescription = state to which the spin projection has to be evaluated
// return value = twice spin projection on (Oz)

inline int Spin2ChainWithTranslations::GetTotalSz (unsigned long stateDescription)
{
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      switch (stateDescription & 0x7ul)
	{
	case 0x4ul:
	  TmpSz += 4;
	  break;
	case 0x3ul:
	  TmpSz += 2;
	  break;
	case 0x1ul:
	  TmpSz -= 2;
	  break;
	case 0x0ul:
	  TmpSz -= 4;
	  break;
	}
      stateDescription >>= 3;
    }
  return TmpSz;
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double Spin2ChainWithTranslations::SziSzj (int i, int j, int state)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = (tmpState >> (j * 3)) & 0x7ul;
  tmpState >>= i * 3;
  tmpState &= 0x7ul;
  switch (tmpState | (tmpState2 << 4))
    {
    case 0x44ul:
      return 4.0;
    case 0x00ul:
      return 4.0;
    case 0x40ul:
      return -4.0;
    case 0x04ul:
      return -4.0;      
    case 0x33ul:
      return 1.0;
    case 0x11ul:
      return 1.0;
    case 0x31ul:
      return -1.0;
    case 0x13ul:
      return -1.0;      
    case 0x43ul:
      return 2.0;
    case 0x34ul:
      return 2.0;
    case 0x10ul:
      return 2.0;
    case 0x01ul:
      return 2.0;
    case 0x41ul:
      return -2.0;
    case 0x14ul:
      return -2.0;      
    case 0x30ul:
      return -2.0;
    case 0x03ul:
      return -2.0;      
    default: 
      return 0.0;
    }
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin2ChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i *= 3;
  tmpState >>= i;
  tmpState &= 0x7ul;
  switch (tmpState)
    {
    case 0x4ul:
      {
	coefficient = 2.0;
	State &= ~(0x4ul << i);
	State |= (0x3ul << i);
      }
      break;
    case 0x3ul:
      {
	coefficient = M_SQRT6;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT6;
	State &= ~(0x2ul << i);
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  j *= 3;
  tmpState2 >>= j;
  tmpState2 &= 0x7ul;
  switch (tmpState2)
    {
    case 0x4ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x3ul:
      {
	coefficient *= 2.0;
	State &= ~(0x3ul << j);
	State |= (0x4ul << j);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT6;
	State |= (0x1ul << j);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT6;
	State &= ~(0x1ul << j);
	State |= (0x2ul << j);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    case 0x0ul:
      {
	coefficient *= 2.0;
	State |= (0x1ul << j);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i operator on a given state (only valid if there is no constraint on total Sz)
//
// i = operator position
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin2ChainWithTranslations::Spi (int i, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  i *= 3;
  tmpState >>= i;
  tmpState &= 0x7ul;
  switch (tmpState)
    {
    case 0x4ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x3ul:
      {
	coefficient = 2.0;
	State &= ~(0x3ul << i);
	State |= (0x4ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT6;
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient = M_SQRT6;
	State &= ~(0x1ul << i);
	State |= (0x2ul << i);
      }
      break;
    case 0x0ul:
      {
	coefficient = 2.0;
	State |= (0x1ul << i);
      }
      break;
    }	  
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}

// return index of resulting state from application of S-_i operator on a given state (only valid if there is no constraint on total Sz)
//
// i = operator position
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin2ChainWithTranslations::Smi (int i, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  i *= 3;
  tmpState >>= i;
  tmpState &= 0x7ul;
  switch (tmpState)
    {
    case 0x4ul:
      {
	coefficient = 2.0;
	State &= ~(0x4ul << i);
	State |= (0x3ul << i);
      }
      break;
    case 0x3ul:
      {
	coefficient = M_SQRT6;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT6;
	State &= ~(0x2ul << i);
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}
    
// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin2ChainWithTranslations::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long State = TmpState;
  i *= 3;
  TmpState >>= i;
  TmpState &= 0x7ul;
  switch (TmpState)
    {
    case 0x4ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x3ul:
      {
	coefficient = 2.0;
	State &= ~(0x3ul << i);
	State |= (0x4ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT6;
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient = M_SQRT6;
	State &= ~(0x1ul << i);
	State |= (0x2ul << i);
      }
      break;
    case 0x0ul:
      {
	coefficient = 2.0;
	State |= (0x1ul << i);
      }
      break;
    }	  
  TmpState = State;
  j *= 3;
  TmpState >>= j;
  TmpState &= 0x7ul;
  switch (TmpState)
    {
    case 0x4ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x3ul:
      {
	coefficient *= 2.0;
	State &= ~(0x3ul << j);
	State |= (0x4ul << j);
	return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT6;
	State |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT6;
	State &= ~(0x1ul << j);
	State |= (0x2ul << j);
      }
      break;
    case 0x0ul:
      {
	coefficient *= 2.0;
	State |= (0x1ul << j);
      }
      break;
    }	  
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin2ChainWithTranslations::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long State = TmpState;
  i *= 3;
  TmpState >>= i;
  TmpState &= 0x7ul;
  switch (TmpState)
    {
    case 0x4ul:
      {
	coefficient = 2.0;
	State &= ~(0x4ul << i);
	State |= (0x3ul << i);
      }
      break;
    case 0x3ul:
      {
	coefficient = M_SQRT6;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT6;
	State &= ~(0x2ul << i);
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  TmpState = State;
  j *= 3;
  TmpState >>= j;
  TmpState &= 0x7ul;
  switch (TmpState)
    {
    case 0x4ul:
      {
	coefficient *= 2.0;
	State &= ~(0x4ul << j);
	State |= (0x3ul << j);
      }
      break;
    case 0x3ul:
      {
	coefficient *= M_SQRT6;
	State &= ~(0x1ul << j);
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT6;
	State &= ~(0x2ul << j);
	State |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	State &= ~(0x1ul << j);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}

// return index of resulting state from application of S+_i S+_i operator on a given state
//
// i = position of first S+ operator
// state = index of the state to be applied on S+_i S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin2ChainWithTranslations::SpiSpi (int i, int state, double& coefficient, int& nbrTranslation)
{
  return this->SpiSpj (i, i, state, coefficient, nbrTranslation);
}

// return index of resulting state from application of S-_i S-_i operator on a given state
//
// i = position of the S- operator
// state = index of the state to be applied on S-_i S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin2ChainWithTranslations::SmiSmi (int i, int state, double& coefficient, int& nbrTranslation)
{
  return this->SmiSmj (i, i, state, coefficient, nbrTranslation);
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin2ChainWithTranslations::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long State = TmpState;
  TmpState >>= (j * 3);
  TmpState &= 0x7ul;
  switch (TmpState)
    {
    case 0x4ul:
      coefficient *= 2.0;
      break;
    case 0x2ul:
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x1ul:
      coefficient *= -1.0;
      break;
    case 0x0ul:
      coefficient *= -2.0;
      break;
    }
  TmpState = State;
  i *= 3;
  TmpState >>= i;
  TmpState &= 0x7ul;
  switch (TmpState)
    {
    case 0x4ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x3ul:
      {
	coefficient *= 2.0;
	State &= ~(0x3ul << i);
	State |= (0x4ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT6;
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT6;
	State &= ~(0x1ul << i);
	State |= (0x2ul << i);
      }
      break;
    case 0x0ul:
      {
	coefficient *= 2.0;
	State |= (0x1ul << i);
      }
      break;
    }	  
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin2ChainWithTranslations::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long State = TmpState;
  TmpState >>= (j * 3);
  TmpState &= 0x7ul;
  switch (TmpState)
    {
    case 0x4ul:
      coefficient *= 2.0;
      break;
    case 0x2ul:
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x1ul:
      coefficient *= -1.0;
      break;
    case 0x0ul:
      coefficient *= -2.0;
      break;
    }
  TmpState = State;
  i *= 3;
  TmpState >>= i;
  TmpState &= 0x7ul;
  switch (TmpState)
    {
    case 0x4ul:
      {
	coefficient *= 2.0;
	State &= ~(0x4ul << j);
	State |= (0x3ul << j);
      }
      break;
    case 0x3ul:
      {
	coefficient *= M_SQRT6;
	State &= ~(0x1ul << j);
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT6;
	State &= ~(0x2ul << j);
	State |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	State &= ~(0x1ul << j);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}

// return index of resulting state from application of S-_i1 S+_j1 S-_i2 S+_j2 operator on a given state
//
// i1 = position of leftmost S- operator
// j1 = position of leftmost S+ operator
// i2 = position of rightmost S- operator
// j2 = position of rightmost S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state (orbit index)

int Spin2ChainWithTranslations::SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient, int& nbrTranslation)
{
  cout << "Spin2ChainWithTranslations::SmiSpjSmiSpj is not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i1 Sz_j1 S-_i2 S+_j2 operator on a given state
//
// i1 = position of first Sz operator
// j1 = position of second Sz operator
// i2 = position of S- operator
// j2 = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state (orbit index)

int Spin2ChainWithTranslations::SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient, int& nbrTranslation)
{
  cout << "Spin2ChainWithTranslations::SziSzjSmiSpj is not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin2ChainWithTranslations::PrintState (ostream& Str, int state)
{
  unsigned long tmp;
  unsigned long TmpStateDescription = this->StateDescription[state];  
  for (int j = 0; j < this->ChainLength; j++)
    {
      tmp = TmpStateDescription & 0x7ul;
      switch (TmpStateDescription & 0x7ul)
	{
	case 0x0ul:
	  Str << "-2 ";
	  break;
	case 0x1ul:
	  Str << "-1 ";
	  break;
	case 0x2ul:
	  Str << "0 ";
	  break;
	case 0x3ul:
	  Str << "1 ";
	  break;
	case 0x4ul:
	  Str << "2 ";
	  break;
	}
      TmpStateDescription >>= 3;
    }
  Str << " (" << this->NbrStateInOrbit[state] << ")";
  return Str;
}


// evaluate Hilbert space dimension with no constraint on the total Sz
//
// nbrSites = number of sites
// return value = Hilbert space dimension

long Spin2ChainWithTranslations::EvaluateHilbertSpaceDimension(int nbrSites)
{
  long Tmp = 5l;
  for (int i = 1; i < this->ChainLength; ++i)
    {
      Tmp *= 5l;
    }
  return Tmp;
}


// evaluate Hilbert space dimension
//
// sz = twice the Sz value
// nbrSites = number of sites
// return value = Hilbert space dimension

long Spin2ChainWithTranslations::EvaluateHilbertSpaceDimension(int sz, int nbrSites)
{
  if (nbrSites == 0)
    {
      if (sz == 0)
	{
	  return 1l;	  
	}
      else
	{
	  return 0l;	  
	}
    }
  long TmpDimension = this->EvaluateHilbertSpaceDimension(sz - 4, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(sz - 2, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(sz, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(sz + 2, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(sz + 4, nbrSites - 1);
  return TmpDimension;
}

// generate all states with no constraint on total Sz and no discrete symmtry constraint
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// return value = number of generated states

long Spin2ChainWithTranslations::RawGenerateStates(long statePosition, int sitePosition) 
{
  if (sitePosition < 0)
    {
      this->StateDescription[statePosition] = 0x0ul;
      return (statePosition + 1l);
    }

  unsigned long TmpMask = 0x4ul << (sitePosition * 3);
  long TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = 0x3ul << (sitePosition * 3);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = 0x2ul << (sitePosition * 3);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = 0x1ul << (sitePosition * 3);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  return this->RawGenerateStates(statePosition, sitePosition - 1);  
}

// generate all states corresponding to a given total Sz and no discrete symmtry constraint
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentSz = total Sz value of current state
// return value = number of generated states

long Spin2ChainWithTranslations::RawGenerateStates(long statePosition, int sitePosition, int currentSz) 
{
  if (sitePosition < 0)
    {
      if (currentSz == this->Sz)
	{
	  this->StateDescription[statePosition] = 0x0ul;
	  return (statePosition + 1l);
	}
      else
	{
	  return statePosition;
	}
    }
  unsigned long TmpMask = 0x4ul << (sitePosition * 3);
  long TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1, currentSz + 4);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = 0x3ul << (sitePosition * 3);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1, currentSz + 2);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = 0x2ul << (sitePosition * 3);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1, currentSz);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  TmpMask = 0x1ul << (sitePosition * 3);
  TmpPosition = this->RawGenerateStates(statePosition, sitePosition - 1, currentSz - 2);  
  for (; statePosition < TmpPosition; ++statePosition)
    this->StateDescription[statePosition] |= TmpMask;
  return this->RawGenerateStates(statePosition, sitePosition - 1, currentSz - 4);  
}


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void Spin2ChainWithTranslations::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->ChainLength);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  int TmpMaxBitPosition = 3 * this->ChainLength;
  if (this->MaximumLookUpShift > TmpMaxBitPosition)
    this->MaximumLookUpShift = TmpMaxBitPosition;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [TmpMaxBitPosition];
  this->LookUpTableShift = new int [TmpMaxBitPosition];
  for (int i = 0; i < TmpMaxBitPosition; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentMaxMomentum = TmpMaxBitPosition;
  while (((this->StateDescription[0] >> CurrentMaxMomentum) == 0x0ul) && (CurrentMaxMomentum > 0))
    --CurrentMaxMomentum;
  int* TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
  if (CurrentMaxMomentum < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentMaxMomentum] = 0;
  else
    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      int TmpMaxMomentum = CurrentMaxMomentum;
      while (((this->StateDescription[i] >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
	--TmpMaxMomentum;
      if (CurrentMaxMomentum != TmpMaxMomentum)
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  --CurrentMaxMomentum;
	  while (CurrentMaxMomentum > TmpMaxMomentum)
	    {
	      this->LookUpTableShift[CurrentMaxMomentum] = -1;
	      --CurrentMaxMomentum;
	    }
 	  CurrentMaxMomentum = TmpMaxMomentum;
	  TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
	  if (CurrentMaxMomentum < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentMaxMomentum] = 0;
	  else
	    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;
  this->ComputeRescalingFactors();
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix Spin2ChainWithTranslations::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, 1);
          Complex Tmp(1.0, 0.0);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
      
    }
  if (nbrSites == this->ChainLength)
    {
      if (szSector == this->Sz)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, 1);
          Complex Tmp(1.0, 0.0);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}      
    }
  Spin2Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin2Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);

  ComplexMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = 3 * nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  int TmpNbrTranslation;
  int TmpNbrTranslationToIdentity;
  Complex* TmpPhases = new Complex [2 * this->ChainLength];
  double Coef = 2.0 * M_PI * ((double) this->Momentum) / ((double) this->ChainLength);
  for (int i = 0; i < (2 * this->ChainLength); ++i)
    {
      TmpPhases[i] = Phase(Coef * ((double) i));
    }

  unsigned long Mask1 = (0x1ul << Shift) - 0x1ul;
  unsigned long Mask2 = (0x1ul << (3 * this->ChainLength)) - 0x1ul;
  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = (TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask1)) & Mask2;
	  double Coefficient = 1.0;
	  int TmpPos = this->SymmetrizeResult(TmpState2, 1, Coefficient, TmpNbrTranslation);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos] * TmpPhases[TmpNbrTranslation] / sqrt((double) this->NbrStateInOrbit[TmpPos]));
	    }
	}
    }
  delete[] TmpPhases;
  return TmpEntanglementMatrix;
}

// build the projectors of two spins on a given total spin sector
//

void Spin2ChainWithTranslations::BuildProjectors()
{
  this->Spin3ProjectorNbrCoefficients = new int[64];
  this->Spin4ProjectorNbrCoefficients = new int[64];
  RealMatrix TmpS4 (25, 9, true);
  RealMatrix TmpS3 (25, 7, true);
  for (int i = 0; i < 64; ++i)
    {
      this->Spin3ProjectorNbrCoefficients[i] = 0;
      this->Spin4ProjectorNbrCoefficients[i] = 0;      
    }

  TmpS4.SetMatrixElement((0 * 5) + 0, 0, 1.0);

  TmpS4.SetMatrixElement((0 * 5) + 1, 1, sqrt(1.0 / 2.0));
  TmpS4.SetMatrixElement((1 * 5) + 0, 1, sqrt(1.0 / 2.0));

  TmpS4.SetMatrixElement((0 * 5) + 2, 2, sqrt(3.0 / 14.0));
  TmpS4.SetMatrixElement((1 * 5) + 1, 2, sqrt(4.0 / 7.0));
  TmpS4.SetMatrixElement((2 * 5) + 0, 2, sqrt(3.0 / 14.0));

  TmpS4.SetMatrixElement((0 * 5) + 3, 3, sqrt(1.0 / 14.0));
  TmpS4.SetMatrixElement((1 * 5) + 2, 3, sqrt(3.0 / 7.0));
  TmpS4.SetMatrixElement((2 * 5) + 1, 3, sqrt(3.0 / 7.0));
  TmpS4.SetMatrixElement((3 * 5) + 0, 3, sqrt(1.0 / 14.0));

  TmpS4.SetMatrixElement((0 * 5) + 4, 3, sqrt(1.0 / 70.0));
  TmpS4.SetMatrixElement((1 * 5) + 3, 3, sqrt(8.0 / 35.0));
  TmpS4.SetMatrixElement((2 * 5) + 2, 3, sqrt(18.0 / 35.0));
  TmpS4.SetMatrixElement((3 * 5) + 1, 3, sqrt(8.0 / 35.0));
  TmpS4.SetMatrixElement((4 * 5) + 0, 3, sqrt(1.0 / 70.0));

  TmpS4.SetMatrixElement((4 * 5) + 1, 5, sqrt(1.0 / 14.0));
  TmpS4.SetMatrixElement((3 * 5) + 2, 5, sqrt(3.0 / 7.0));
  TmpS4.SetMatrixElement((2 * 5) + 3, 5, sqrt(3.0 / 7.0));
  TmpS4.SetMatrixElement((1 * 5) + 4, 5, sqrt(1.0 / 14.0));

  TmpS4.SetMatrixElement((4 * 5) + 2, 6, sqrt(3.0 / 14.0));
  TmpS4.SetMatrixElement((3 * 5) + 3, 6, sqrt(4.0 / 7.0));
  TmpS4.SetMatrixElement((2 * 5) + 4, 6, sqrt(3.0 / 14.0));

  TmpS4.SetMatrixElement((4 * 5) + 3, 7, sqrt(1.0 / 2.0));
  TmpS4.SetMatrixElement((3 * 5) + 4, 7, sqrt(1.0 / 2.0));

  TmpS4.SetMatrixElement((4 * 5) + 4, 8, 1.0);


  TmpS3.SetMatrixElement((0 * 5) + 1, 0, sqrt(1.0 / 2.0));
  TmpS3.SetMatrixElement((1 * 5) + 0, 0, -sqrt(1.0 / 2.0));

  TmpS3.SetMatrixElement((0 * 5) + 2, 1, sqrt(1.0 / 2.0));
  TmpS3.SetMatrixElement((2 * 5) + 0, 1, -sqrt(1.0 / 2.0));

  TmpS3.SetMatrixElement((0 * 5) + 3, 2, sqrt(3.0 / 10.0));
  TmpS3.SetMatrixElement((1 * 5) + 2, 2, sqrt(1.0 / 5.0));
  TmpS3.SetMatrixElement((2 * 5) + 1, 2, -sqrt(1.0 / 5.0));
  TmpS3.SetMatrixElement((3 * 5) + 0, 2, -sqrt(3.0 / 10.0));

  TmpS3.SetMatrixElement((0 * 5) + 4, 3, sqrt(1.0 / 10.0));
  TmpS3.SetMatrixElement((1 * 5) + 3, 3, sqrt(2.0 / 5.0));
  TmpS3.SetMatrixElement((3 * 5) + 1, 3, -sqrt(2.0 / 5.0));
  TmpS3.SetMatrixElement((4 * 5) + 0, 3, -sqrt(1.0 / 10.0));

  TmpS3.SetMatrixElement((4 * 5) + 1, 4, -sqrt(3.0 / 10.0));
  TmpS3.SetMatrixElement((3 * 5) + 2, 4, -sqrt(1.0 / 5.0));
  TmpS3.SetMatrixElement((2 * 5) + 3, 4, sqrt(1.0 / 5.0));
  TmpS3.SetMatrixElement((1 * 5) + 4, 4, sqrt(3.0 / 10.0));

  TmpS3.SetMatrixElement((4 * 5) + 2, 5, -sqrt(1.0 / 2.0));
  TmpS3.SetMatrixElement((2 * 5) + 4, 5, sqrt(1.0 / 2.0));

  TmpS3.SetMatrixElement((4 * 5) + 3, 6, -sqrt(1.0 / 2.0));
  TmpS3.SetMatrixElement((3 * 5) + 4, 6, sqrt(1.0 / 2.0));

  RealMatrix TmpS4Transpose = TmpS4.DuplicateAndTranspose();
  RealMatrix TmpS3Transpose = TmpS3.DuplicateAndTranspose();
  RealMatrix TmpS4Projector = TmpS4 * TmpS4Transpose;
  RealMatrix TmpS3Projector = TmpS3 * TmpS3Transpose;
  
//  cout <<  TmpS3Projector << endl;
  for (int i = 0; i < 5; ++i)
    {      
      for (int j = 0; j < 5; ++j)
	{
	  for (int k = 0; k < 25; ++k)
	    {
	      if (TmpS3Projector[k][(i * 5) + j] != 0.0)
		{
		  this->Spin3ProjectorNbrCoefficients[(i << 3) | j]++;
		}
	      if (TmpS4Projector[k][(i * 5) + j] != 0.0)
		{
		  this->Spin4ProjectorNbrCoefficients[(i << 3) | j]++;
		}
	    }
	}
    }
  this->Spin3ProjectorCoefficients = new double*[64];
  this->Spin3ProjectorStates = new unsigned long*[64];
  this->Spin4ProjectorCoefficients = new double*[64];
  this->Spin4ProjectorStates = new unsigned long*[64];
  for (int i = 0; i < 64; ++i)
    {
      if (this->Spin3ProjectorNbrCoefficients[i] > 0)
	{
	  this->Spin3ProjectorCoefficients[i] = new double[this->Spin3ProjectorNbrCoefficients[i]];
	  this->Spin3ProjectorStates[i] = new unsigned long[this->Spin3ProjectorNbrCoefficients[i]];
	  this->Spin3ProjectorNbrCoefficients[i] = 0;
	}
      if (this->Spin4ProjectorNbrCoefficients[i] > 0)
	{
	  this->Spin4ProjectorCoefficients[i] = new double[this->Spin4ProjectorNbrCoefficients[i]];
	  this->Spin4ProjectorStates[i] = new unsigned long[this->Spin4ProjectorNbrCoefficients[i]];
	  this->Spin4ProjectorNbrCoefficients[i] = 0;
	}
    }
  for (int i = 0; i < 5; ++i)
    {      
      for (int j = 0; j < 5; ++j)
	{
	  for (int k = 0; k < 5; ++k)
	    {
	      for (int l = 0; l < 5; ++l)
		{
		  if (TmpS3Projector[(k * 5) + l][(i * 5) + j] != 0.0)
		    {
		      this->Spin3ProjectorCoefficients[(i << 3) | j][this->Spin3ProjectorNbrCoefficients[(i << 3) | j]] = TmpS3Projector[(k * 5) + l][(i * 5) + j];
		      this->Spin3ProjectorStates[(i << 3) | j][this->Spin3ProjectorNbrCoefficients[(i << 3) | j]] = (unsigned long) ((l << 3) | k);
		      this->Spin3ProjectorNbrCoefficients[(i << 3) | j]++;
		    }
		  if (TmpS4Projector[(k * 5) + l][(i * 5) + j] != 0.0)
		    {
		      this->Spin4ProjectorCoefficients[(i << 3) | j][this->Spin4ProjectorNbrCoefficients[(i << 3) | j]] = TmpS4Projector[(k * 5) + l][(i * 5) + j];
		      this->Spin4ProjectorStates[(i << 3) | j][this->Spin4ProjectorNbrCoefficients[(i << 3) | j]] = (unsigned long) ((l << 3) | k);
		      this->Spin4ProjectorNbrCoefficients[(i << 3) | j]++;
		    }
		}
	    }
	}
    }

}

// compute all the states connected to a single one by a two site spin projector
//
// i = index of the first site
// j = index of the first site
// state = index of state on whcih the projector should be applied
// indices = pointer to the array where the connected state indeices will be stored
// coefficients = pointer to the array where coefficients related to each connected states will be stored
// nbrTranslations  = pointer to the array where number of translations related to each connected states will be stored
// spinProjectorNbrCoefficients = number of coefficients per component for the spin projector
// spinProjectorCoefficients = coefficients per component for the spin projector
// spinProjectorStates = configuration per component for the spin projector
// return value = number of connected states

int Spin2ChainWithTranslations::SpinProjector (int  i, int j, int state, int* indices, double* coefficients, int* nbrTranslations,
					       int* spinProjectorNbrCoefficients,  double** spinProjectorCoefficients, unsigned long** spinProjectorStates)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long TmpIndex = (((TmpState >> (3 * i)) & 0x7) << 3) | ((TmpState >> (3 * j)) & 0x7);
  TmpState &= ~(0x7ul << (3 * i));
  TmpState &= ~(0x7ul << (3 * j));
  unsigned long TmpState2;
  int TmpNbrIndex = 0;
  for (int k = 0; k < spinProjectorNbrCoefficients[TmpIndex]; ++k)
    {
      unsigned long TmpIndex2 = spinProjectorStates[TmpIndex][k]; 
      TmpState2 = TmpState;
      TmpState2 |= (TmpIndex2 & 0x7ul) << (3 * i);
      TmpState2 |= (TmpIndex2 >> 3) << (3 * j);
      double TmpCoefficient = spinProjectorCoefficients[TmpIndex][k];
      int TmpNbrTranslation = 0;
      int TmpIndex3 = this->SymmetrizeResult(TmpState2, this->NbrStateInOrbit[state], TmpCoefficient, TmpNbrTranslation);
      if (TmpIndex3 < this->HilbertSpaceDimension)
	{
	  indices[TmpNbrIndex] = TmpIndex3;
	  coefficients[TmpNbrIndex] = TmpCoefficient;
	  nbrTranslations[TmpNbrIndex] = TmpNbrTranslation;
	  ++TmpNbrIndex;
	}
    }
  return TmpNbrIndex;
}
