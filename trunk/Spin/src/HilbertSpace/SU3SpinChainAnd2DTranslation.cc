////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                        Class author Cecile Repellin                        //
//                                                                            //
//                                                                            //
//                class of SU(3) spin chain with Sz contraint                 //
//                               and 2d translations                          //
//                                                                            //
//                        last modification : 06/02/2018                      //
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


#include "HilbertSpace/SU3SpinChainAnd2DTranslation.h"
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



// default constructor
//

SU3SpinChainAnd2DTranslation::SU3SpinChainAnd2DTranslation ()
{
  this->Flag.Initialize();
  this->ChainLength = 0;
  this->HilbertSpaceDimension = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->Tz = 0;
  this->Y = 0;
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

SU3SpinChainAnd2DTranslation::SU3SpinChainAnd2DTranslation (int nbrSite, int tz, int y, int xMomentum, int maxXMomentum, int yMomentum, int maxYMomentum, unsigned long memory) 
{
  this->Flag.Initialize();
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
  this->ChainLength = nbrSite;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;

  this->Tz = tz;
  this->Y = y;
  this->FixedSpinProjectionFlag = true;
    
  this->StateXShift = 2 * (this->ChainLength / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this-> ChainLength) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->NbrYMomentumBlocks = 2 * this->ChainLength / this->StateXShift;
  this->StateYShift = 2 * (this->ChainLength / (this->MaxYMomentum * this->MaxXMomentum));
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

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(0, 0, this->ChainLength);
  cout << "Hilbert space dimension before taking into account translation symmetries = " << (this->LargeHilbertSpaceDimension) << endl;
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
    
//      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//       {
// 	unsigned long state = this->StateDescription[i];
// 	int TmpMaxMomentum = 2*this->ChainLength;
// 	while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
// 	  --TmpMaxMomentum;
// 	cout << "Space : " ;
// 	this->PrintState(cout, i) << " " << i << " " << this->FindStateIndex(state, TmpMaxMomentum) << endl;
//       }
}

// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

SU3SpinChainAnd2DTranslation::SU3SpinChainAnd2DTranslation (const SU3SpinChainAnd2DTranslation& chain)
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
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

      this->Tz = chain.Tz;
      this->Y = chain.Y;

      this->StateDescription = chain.StateDescription;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;      
    }
  else
    {
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->LargeHilbertSpaceDimension = 0;
      this->Tz = 0;
      this->Y = 0;
      this->StateDescription = 0;
      this->LookUpTable = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->LookUpTableSize = 0;
    }
}

// destructor
//

SU3SpinChainAnd2DTranslation::~SU3SpinChainAnd2DTranslation () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

SU3SpinChainAnd2DTranslation& SU3SpinChainAnd2DTranslation::operator = (const SU3SpinChainAnd2DTranslation& chain)
{
  if ((this->Flag.Used() == true) && (this->ChainLength != 0) && (this->Flag.Shared() == false))
    {
      delete[] this->StateDescription;
    }
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
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

      this->Tz = chain.Tz;
      this->Y = chain.Y;

      this->StateDescription = chain.StateDescription;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;

      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;      
    }
  else
    {
      this->ChainLength = 0;
      this->HilbertSpaceDimension = 0;
      this->Tz = 0;
      this->Y = 0;
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

AbstractHilbertSpace* SU3SpinChainAnd2DTranslation::Clone()
{
  return new SU3SpinChainAnd2DTranslation (*this);
}

// return index of resulting state from application of Pij operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied onPij operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int SU3SpinChainAnd2DTranslation::Pij (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long TmpMask1 = ((tmpState >> (i << 1)) & 0x3ul) ;
  unsigned long TmpMask2 = ((tmpState >> (j << 1)) & 0x3ul);
  coefficient = 1.0;
  nbrTranslationX = 0;
  nbrTranslationY = 0;
  
  if (TmpMask1 == TmpMask2)
    return state;
  
  unsigned long TmpMask3 = 0x3ul << (i << 1);
  unsigned long TmpMask4 = 0x3ul << (j << 1);
  
  
  tmpState &= (~(TmpMask3));
  tmpState &= (~(TmpMask4));
  
  tmpState |= (TmpMask1 << (j << 1));
  tmpState |= (TmpMask2 << (i << 1));
  
  return this->SymmetrizeResult(tmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
 
}


// return index of resulting state from application of four-site exchange operator on a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = index of resulting state

int SU3SpinChainAnd2DTranslation::Pijkl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long TmpMask1 = ((tmpState >> (i << 1)) & 0x3ul) ;
  unsigned long TmpMask2 = ((tmpState >> (j << 1)) & 0x3ul);
  unsigned long TmpMask3 = ((tmpState >> (k << 1)) & 0x3ul) ;
  unsigned long TmpMask4 = ((tmpState >> (l << 1)) & 0x3ul);
  coefficient = 1.0;
  nbrTranslationX = 0;
  nbrTranslationY = 0;
  
  if ((TmpMask1 == TmpMask2) && (TmpMask2 == TmpMask3) && (TmpMask3 == TmpMask4))
    return state;
  
  
  unsigned long TmpMask = (0x3ul << (i << 1)) | (0x3ul << (j << 1)) | (0x3ul << (k << 1)) | (0x3ul << (l << 1)); 
  
  tmpState &= (~TmpMask);
  
  tmpState |= (TmpMask1 << (j << 1));
  tmpState |= (TmpMask2 << (k << 1));
  tmpState |= (TmpMask3 << (l << 1));
  tmpState |= (TmpMask4 << (i << 1));
  
  double TmpCoefficient;
  return this->SymmetrizeResult(tmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslationX, nbrTranslationY);
}

// find state index
//
// stateDescription = unsigned longeger describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int SU3SpinChainAnd2DTranslation::FindStateIndex(unsigned long stateDescription, int maxMomentum)
{
  if ((stateDescription > this->StateDescription[0]) || (stateDescription < this->StateDescription[this->HilbertSpaceDimension - 1]))
    return this->HilbertSpaceDimension;
  if (this->LookUpTableShift[maxMomentum] < 0)
    return this->HilbertSpaceDimension;
  long PosMax = stateDescription >> this->LookUpTableShift[maxMomentum];
  long PosMin = this->LookUpTable[maxMomentum][PosMax];
  PosMax = this->LookUpTable[maxMomentum][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }

  if (CurrentState == stateDescription)
    return PosMid;
  else
    if ((this->StateDescription[PosMin] != stateDescription) && (this->StateDescription[PosMax] != stateDescription))
      return this->HilbertSpaceDimension;
    else
      return PosMin;
}

// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long SU3SpinChainAnd2DTranslation::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(0l, this->ChainLength - 1, 0, 0);
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

void SU3SpinChainAnd2DTranslation::ComputeRescalingFactors()
{
  this->RescalingFactors = new double* [this->ChainLength + 1];
  for (int i = 1; i <= this->ChainLength; ++i)
    {
      this->RescalingFactors[i] = new double [this->ChainLength + 1];
      for (int j = 1; j <= this->ChainLength; ++j)
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

ComplexVector SU3SpinChainAnd2DTranslation::ConvertToKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  SU3SpinChain* TmpSpace = (SU3SpinChain*) space;
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

ComplexVector SU3SpinChainAnd2DTranslation::ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  SU3SpinChain* TmpSpace = (SU3SpinChain*) space;
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
      int TmpMaxMomentum = 2*this->ChainLength;
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
  

