////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of spin 3/2 chain with Sz contraint                   //
//                               and 2d translations                          //
//                                                                            //
//                        last modification : 11/07/2018                      //
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


#include "HilbertSpace/Spin3_2ChainAnd2DTranslation.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"

#include <math.h>
#include <iostream>


using std::cout;
using std::dec;
using std::hex;
using std::endl;


#ifndef M_SQRT3
#define M_SQRT3 1.73205080756887729
#endif


// default constructor
//

Spin3_2ChainAnd2DTranslation::Spin3_2ChainAnd2DTranslation ()
{
  this->Flag.Initialize();
  this->ChainLength = 0;
  this->HilbertSpaceDimension = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->Sz = 0;
  this->StateDescription = 0;
  this->LookUpTable = 0;
  this->LookUpTableMask = 0;
  this->LookUpPosition = 0;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// xMomentum = momentum along the x direction
// maxXMomentum = number of sites in the x direction
// yMomentum = momentum along the y direction
// maxYMomentum = number of sites in the y direction
// memory = amount of memory granted for precalculations

Spin3_2ChainAnd2DTranslation::Spin3_2ChainAnd2DTranslation (int nbrSite, int sz, int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, unsigned long memory) 
{
  this->Flag.Initialize();
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
  this->NbrSite = nbrSite;
  this->ChainLength = this->NbrSite;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;

  this->Sz = sz;
  this->FixedQuantumNumberFlag = true;
  
  this->StateXShift = 2 * (this->NbrSite / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this->NbrSite) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->NbrYMomentumBlocks = (2 * this->NbrSite) / this->StateXShift;
  this->StateYShift = 2 * (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
  this->YMomentumFullMask = 0x0ul;
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
    }
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->NbrSite);
  this->LargeHilbertSpaceDimension = this->GenerateStates();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
 
  if (this->LargeHilbertSpaceDimension > 0)
    {
      this->GenerateLookUpTable(memory);
      this->ComputeRescalingFactors();

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
      UsedMemory = this->NbrSite * sizeof(int);
      UsedMemory += this->NbrSite * this->LookUpTableMemorySize * sizeof(int);
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

Spin3_2ChainAnd2DTranslation::Spin3_2ChainAnd2DTranslation (const Spin3_2ChainAnd2DTranslation& chain)
{
  this->Flag = chain.Flag;
  if (chain.NbrSite != 0)
    {
      this->NbrSite = chain.NbrSite;
      this->ChainLength = chain.NbrSite;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->MaxXMomentum = chain.MaxXMomentum;
      this->XMomentum = chain.XMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;
      this->MaxYMomentum = chain.MaxYMomentum;
      this->YMomentum = chain.YMomentum;
      this->NbrYMomentumBlocks = chain.NbrYMomentumBlocks;
      this->StateYShift = chain.StateYShift;
      this->YMomentumBlockSize = chain.YMomentumBlockSize;
      this->ComplementaryStateYShift = chain.ComplementaryStateYShift;
      this->YMomentumMask = chain.YMomentumMask;
      this->YMomentumBlockMask = chain.YMomentumBlockMask;  
      this->YMomentumFullMask = chain.YMomentumFullMask;
      this->ComplementaryYMomentumFullMask = chain.ComplementaryYMomentumFullMask; 

      this->Sz = chain.Sz;

      this->StateDescription = chain.StateDescription;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->LookUpTableShift = chain.LookUpTableShift;

      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
    }
  else
    {
      this->NbrSite = 0;
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableShift = 0;
    }
}

// destructor
//

Spin3_2ChainAnd2DTranslation::~Spin3_2ChainAnd2DTranslation () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin3_2ChainAnd2DTranslation& Spin3_2ChainAnd2DTranslation::operator = (const Spin3_2ChainAnd2DTranslation& chain)
{
  if ((this->Flag.Used() == true) && (this->NbrSite != 0) && (this->Flag.Shared() == false))
    {
      delete[] this->StateDescription;
    }
  this->Flag = chain.Flag;
  if (chain.NbrSite != 0)
    {
      this->NbrSite = chain.NbrSite;
      this->ChainLength = chain.NbrSite;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;

      this->MaxXMomentum = chain.MaxXMomentum;
      this->XMomentum = chain.XMomentum;
      this->StateXShift = chain.StateXShift;
      this->ComplementaryStateXShift = chain.ComplementaryStateXShift;
      this->XMomentumMask = chain.XMomentumMask;
      this->MaxYMomentum = chain.MaxYMomentum;
      this->YMomentum = chain.YMomentum;
      this->NbrYMomentumBlocks = chain.NbrYMomentumBlocks;
      this->StateYShift = chain.StateYShift;
      this->YMomentumBlockSize = chain.YMomentumBlockSize;
      this->ComplementaryStateYShift = chain.ComplementaryStateYShift;
      this->YMomentumMask = chain.YMomentumMask;
      this->YMomentumBlockMask = chain.YMomentumBlockMask;  
      this->YMomentumFullMask = chain.YMomentumFullMask;
      this->ComplementaryYMomentumFullMask = chain.ComplementaryYMomentumFullMask; 

      this->Sz = chain.Sz;

      this->StateDescription = chain.StateDescription;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;      
    }
  else
    {
      this->NbrSite = 0;
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Sz = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin3_2ChainAnd2DTranslation::Clone()
{
  return new Spin3_2ChainAnd2DTranslation (*this);
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int Spin3_2ChainAnd2DTranslation::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i *= 2;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      {
	coefficient = M_SQRT3;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient = 2.0;
	State &= ~(0x2ul << i);
	State |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient = M_SQRT3;
	State &= ~(0x1ul << i);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  j *= 2;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  switch (tmpState2)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	State |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	State &= ~(0x1ul << j);
	State |= (0x2ul << j);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	State |= (0x1ul << j);
      }
      break;
    }	  
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
}

// return index of resulting state from application of Sz_i Sz_j S-_k S+_l operator on a given state
//
// i = position of the first Sz operator
// j = position of the second Sz operator
// k = position of S- operator
// l = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int Spin3_2ChainAnd2DTranslation::SziSzjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{  
  unsigned long TmpState = this->StateDescription[state];
  l *= 2;
  switch ((TmpState >> l) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << l);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	TmpState &= ~(0x1ul << l);
	TmpState |= (0x2ul << l);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << l);
      }
      break;
    }	  
  k *= 2;
  switch ((TmpState >> k) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << k);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << k);
	TmpState |= (0x1ul << k);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << k);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	   
  j *= 2;
  switch ((TmpState >> j) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }  
  i *= 2;
  switch ((TmpState >> i) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }  
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
}

// return index of resulting state from application of S-_i S+_j S-_k S+_l operator on a given state
//
// i = position of the first S- operator
// j = position of the first S+ operator
// k = position of the second S- operator
// l = position of the second S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int Spin3_2ChainAnd2DTranslation::SmiSpjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{  
  unsigned long TmpState = this->StateDescription[state];
  l *= 2;
  switch ((TmpState >> l) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << l);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	TmpState &= ~(0x1ul << l);
	TmpState |= (0x2ul << l);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << l);
      }
      break;
    }	  
  k *= 2;
  switch ((TmpState >> k) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << k);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << k);
	TmpState |= (0x1ul << k);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << k);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	   
  j *= 2;
  switch ((TmpState >> j) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x1ul << j);
	TmpState |= (0x2ul << j);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << j);
      }
      break;
    }	  
  i *= 2;
  switch ((TmpState >> i) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << i);
	TmpState |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
}


// return index of resulting state from application of S-z_i Sz_j Sz_k Sz_l S-_m S+_n operator on a given state
//
// i = position of the first Sz operator
// j = position of the second Sz operator
// k = position of the third Sz operator
// l = position of the fourth Sz operator
// m = position of the first S- operator
// n = position of the first S+ operator
// state = index of the state to be applied on the Sz_i Sz_j Sz_k Sz_l S-_m S+_n operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int Spin3_2ChainAnd2DTranslation::SziSzjSzkSzlSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
						      int& nbrTranslationX, int& nbrTranslationY)
{  
  unsigned long TmpState = this->StateDescription[state];
  n *= 2;
  switch ((TmpState >> n) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << n);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	TmpState &= ~(0x1ul << n);
	TmpState |= (0x2ul << n);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << n);
      }
      break;
    }	  
  m *= 2;
  switch ((TmpState >> m) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << m);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << m);
	TmpState |= (0x1ul << m);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << m);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	   
  l *= 2;
  switch ((TmpState >> l) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }  
  k *= 2;
  switch ((TmpState >> k) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }  
  j *= 2;
  switch ((TmpState >> j) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }  
  i *= 2;
  switch ((TmpState >> i) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }  
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
}
  
// return index of resulting state from application of Sz_i Sz_j S-_k S+_l S-_m S+_n operator on a given state
//
// i = position of the first Sz operator
// j = position of the second Sz operator
// k = position of the first S- operator
// l = position of the first S+ operator
// m = position of the second S- operator
// n = position of the second S+ operator
// state = index of the state to be applied on the Sz_i Sz_j S-_k S+_l S-_m S+_n operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int Spin3_2ChainAnd2DTranslation::SziSzjSmkSplSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
						      int& nbrTranslationX, int& nbrTranslationY)
{  
  unsigned long TmpState = this->StateDescription[state];
  n *= 2;
  switch ((TmpState >> n) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << n);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	TmpState &= ~(0x1ul << n);
	TmpState |= (0x2ul << n);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << n);
      }
      break;
    }	  
  m *= 2;
  switch ((TmpState >> m) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << m);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << m);
	TmpState |= (0x1ul << m);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << m);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	   
  l *= 2;
  switch ((TmpState >> l) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << l);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x1ul << l);
	TmpState |= (0x2ul << l);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << l);
      }
      break;
    }	  
  k *= 2;
  switch ((TmpState >> k) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << k);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << k);
	TmpState |= (0x1ul << k);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << k);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	   
  j *= 2;
  switch ((TmpState >> j) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }  
  i *= 2;
  switch ((TmpState >> i) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }  
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
}

// return index of resulting state from application of S-_i S+_j Sz_k Sz_l S-_m S+_n operator on a given state
//
// i = position of the first S- operator
// j = position of the first S+ operator
// k = position of the first Sz operator
// l = position of the second Sz operator
// m = position of the second S- operator
// n = position of the second S+ operator
// state = index of the state to be applied on the S-_i S+_j Sz_k Sz_l S-_m S+_n operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int Spin3_2ChainAnd2DTranslation::SmiSpjSzkSzlSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
						      int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long TmpState = this->StateDescription[state];
  n *= 2;
  switch ((TmpState >> n) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << n);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	TmpState &= ~(0x1ul << n);
	TmpState |= (0x2ul << n);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << n);
      }
      break;
    }	  
  m *= 2;
  switch ((TmpState >> m) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << m);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << m);
	TmpState |= (0x1ul << m);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << m);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	   
  l *= 2;
  switch ((TmpState >> l) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }  
  k *= 2;
  switch ((TmpState >> k) & 0x3ul)
    {
    case 0x3ul:
      coefficient *= 1.5;
      break;
    case 0x2ul:
      coefficient *= 0.5;
      break;
    case 0x1ul:
      coefficient *= -0.5;
      break;
    case 0x0ul:
      coefficient *= -1.5;
      break;
    }  
  j *= 2;
  switch ((TmpState >> j) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x1ul << j);
	TmpState |= (0x2ul << j);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << j);
      }
      break;
    }	  
  i *= 2;
  switch ((TmpState >> i) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << i);
	TmpState |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
}
  
// return index of resulting state from application of S-_i S+_j S-_k S+_l S-_m S+_n operator on a given state
//
// i = position of the first S- operator
// j = position of the first S+ operator
// k = position of the second S- operator
// l = position of the second S+ operator
// m = position of the third S- operator
// n = position of the third S+ operator
// state = index of the state to be applied on S-_i S+_j S-_k S+_l S-_m S+_n operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int Spin3_2ChainAnd2DTranslation::SmiSpjSmkSplSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
						      int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long TmpState = this->StateDescription[state];
  n *= 2;
  switch ((TmpState >> n) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << n);
      }
      break;
    case 0x1ul:
      {
	coefficient = 2.0;
	TmpState &= ~(0x1ul << n);
	TmpState |= (0x2ul << n);
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT3;
	TmpState |= (0x1ul << n);
      }
      break;
    }	  
  m *= 2;
  switch ((TmpState >> m) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << m);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << m);
	TmpState |= (0x1ul << m);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << m);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	   
  l *= 2;
  switch ((TmpState >> l) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << l);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x1ul << l);
	TmpState |= (0x2ul << l);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << l);
      }
      break;
    }	  
  k *= 2;
  switch ((TmpState >> k) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << k);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << k);
	TmpState |= (0x1ul << k);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << k);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	   
  j *= 2;
  switch ((TmpState >> j) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x1ul << j);
	TmpState |= (0x2ul << j);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << j);
      }
      break;
    }	  
  i *= 2;
  switch ((TmpState >> i) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << i);
	TmpState |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << i);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
}

// return index of resulting state from application of S+_i S-_j Sz_k operator on a given state
//
// i = position of S+ operator
// j = position of S- operator
// k = position of Sz operator
// state = index of the state to be applied on S+_i S-_j Sz_k operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int Spin3_2ChainAnd2DTranslation::SpiSmjSzk (int i, int j, int k, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long TmpState = this->StateDescription[state];
  k *= 2;
  switch ((TmpState >> k) & 0x3ul)
    {
    case 0x3ul:
      coefficient = 1.5;
      break;
    case 0x2ul:
      coefficient = 0.5;
      break;
    case 0x1ul:
      coefficient = -0.5;
      break;
    case 0x0ul:
      coefficient = -1.5;
      break;
    }  
  j *= 2;
  switch ((TmpState >> j) & 0x3ul)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << j);
      }
      break;
    case 0x2ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x2ul << j);
	TmpState |= (0x1ul << j);
      }
      break;
    case 0x1ul:
      {
	coefficient *= M_SQRT3;
	TmpState &= ~(0x1ul << j);
      }
      break;
    case 0x0ul:
      {
	return this->HilbertSpaceDimension;
      }
    }	  
  i *= 2;
  switch ((TmpState >> i) & 0x3ul)
    {
    case 0x3ul:
      {
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << i);
      }
      break;
    case 0x1ul:
      {
	coefficient *= 2.0;
	TmpState &= ~(0x1ul << i);
	TmpState |= (0x2ul << i);
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT3;
	TmpState |= (0x1ul << i);
      }
      break;
    }	  
  return this->SymmetrizeResult(TmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
}


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long Spin3_2ChainAnd2DTranslation::GenerateStates()
{
  
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->NbrSite - 1, 0);
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
#ifdef  __64_BITS__
  unsigned long Discard = 0xfffffffffffffffful;
#else
  unsigned long Discard = 0xfffffffful;
#endif
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY) == this->StateDescription[i]))
	{
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      ++TmpLargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescription[i] = Discard;
	    }
	}
      else
	{
	  this->StateDescription[i] = Discard;
	}
    }
  if (TmpLargeHilbertSpaceDimension > 0)
    {
      unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
      this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
      TmpLargeHilbertSpaceDimension = 0l;
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if (this->StateDescription[i] != Discard)
	    {
	      TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	      this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = this->FindOrbitSize(this->StateDescription[i]);
	      ++TmpLargeHilbertSpaceDimension;
	    }
	}
      delete[] this->StateDescription;
      this->StateDescription = TmpStateDescription;
    }
  else
    {
      delete[] this->StateDescription;
    }
  return TmpLargeHilbertSpaceDimension;
}

// compute the rescaling factors
//

void Spin3_2ChainAnd2DTranslation::ComputeRescalingFactors()
{
  this->RescalingFactors = new double* [this->NbrSite + 1];
  for (int i = 1; i <= this->NbrSite; ++i)
    {
      this->RescalingFactors[i] = new double [this->NbrSite + 1];
      for (int j = 1; j <= this->NbrSite; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}

// convert a state defined in the real space basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector Spin3_2ChainAnd2DTranslation::ConvertToKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  Spin3_2Chain* TmpSpace = (Spin3_2Chain*) space;
  ComplexVector TmpVector (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      int Pos = TmpSpace->FindStateIndex(TmpState);
      if (Pos < TmpSpace->HilbertSpaceDimension)
	{
	  TmpVector[i] =  state[Pos] * sqrt((double) this->NbrStateInOrbit[i]);
	}
    }
  return TmpVector;
}

// convert a state defined in the (Kx,Ky) basis into a state in the real space basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector Spin3_2ChainAnd2DTranslation::ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  Spin3_2Chain* TmpSpace = (Spin3_2Chain*) space;
  ComplexVector TmpVector (TmpSpace->HilbertSpaceDimension, true);
  Complex** FourrierCoefficients = new Complex* [this->MaxXMomentum];
  
  for (int i = 0; i < this->MaxXMomentum; ++i)
    {
      FourrierCoefficients[i] = new Complex [this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{
	  FourrierCoefficients[i][j] = Phase (-2.0 * M_PI * ((double) (i * this->XMomentum) / ((double) this->MaxXMomentum) + (double) (j * this->YMomentum) / ((double) this->MaxYMomentum)));
	}
    }
  
  for (int i = 0; i < TmpSpace->HilbertSpaceDimension; ++i)
    {
      int NbrTranslationX;
      int NbrTranslationY;
      unsigned long TmpState = TmpSpace->StateDescription[i];
      TmpState = this->FindCanonicalForm(TmpState, NbrTranslationX, NbrTranslationY);
      NbrTranslationX = (this->MaxXMomentum - NbrTranslationX) % this->MaxXMomentum;
      NbrTranslationY = (this->MaxYMomentum - NbrTranslationY) % this->MaxYMomentum;
      int TmpMaxMomentum = this->NbrSite;
      while (((TmpState >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
	--TmpMaxMomentum;
      
      int Pos = this->FindStateIndex(TmpState, TmpMaxMomentum);
      if (Pos < this->HilbertSpaceDimension)
	{
	  TmpVector[i] =  (state[Pos] * FourrierCoefficients[NbrTranslationX][NbrTranslationY] / sqrt((double) this->NbrStateInOrbit[Pos]));
	}
    }
  delete[] FourrierCoefficients;
  return TmpVector;
}
  

